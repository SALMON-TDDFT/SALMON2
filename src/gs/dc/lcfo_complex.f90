!
!  Copyright 2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!

#include "config.h"
module lcfo_complex
  implicit none

  private
  public :: check_lcfo_complex_options, dc_lcfo_complex

  real(8), parameter :: lcfo_tol = 1d-10

  type :: s_lcfo_complex_halo
    integer :: id_src = -1
    integer :: id_dst = -1
    integer :: ifrag_src = 0
    integer :: dvec(3) = 0
    integer :: length(3) = 0
    integer :: dsp_send(3) = 0
    integer :: dsp_recv(3) = 0
    complex(8), allocatable :: buf_send(:,:,:,:,:)
    complex(8), allocatable :: buf_recv(:,:,:,:,:)
    complex(8), allocatable :: mat_H_local(:,:,:)
  end type s_lcfo_complex_halo

contains

  subroutine check_lcfo_complex_options(system,info,dc)
    use ieee_arithmetic, only: ieee_is_finite
    use salmon_global, only: energy_cut,lambda_cut,lcfo_eigensolver, &
         num_fragment,yn_dc_lcfo_diag,yn_spinorbit,theory
    use structures, only: s_dcdft,s_dft_system,s_parallel_info
    implicit none
    type(s_dft_system), intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    type(s_dcdft), intent(in) :: dc
    integer :: ik,n

#if defined(USE_OPENACC) || defined(USE_CUDA)
    stop "DC-LCFO complex: GPU/OpenACC/CUDA is unsupported."
#endif
    if (yn_spinorbit == 'y') stop "DC-LCFO complex: spin-orbit and noncollinear calculations are unsupported."
    if (trim(theory) /= 'dft') stop "DC-LCFO complex: theory must be dft."
    if (trim(lcfo_eigensolver) /= 'lapack') &
      stop "DC-LCFO complex: lcfo_eigensolver must be 'lapack'."
    if (yn_dc_lcfo_diag /= 'y') &
      stop "DC-LCFO complex: yn_dc_lcfo_diag='y' is required."
    if (system%if_real_orbital) stop "DC-LCFO complex: complex orbitals are required."
    if (system%nspin < 1 .or. system%nspin > 2) stop "DC-LCFO complex: unsupported spin count."
    if (.not.ieee_is_finite(system%hvol) .or. system%hvol <= 0d0) &
      stop "DC-LCFO complex: invalid grid volume."
    if (system%nk < 1 .or. dc%nstate_tot < 1 .or. dc%nstate_frag < 1) &
      stop "DC-LCFO complex: invalid k-point or state count."
    if (info%npk < 1 .or. info%npk > system%nk .or. info%nporbital < 1 .or. &
        info%nporbital > system%no) stop "DC-LCFO complex: invalid MPI distribution."
    if (info%numk < 1) stop "DC-LCFO complex: every k-point group must own a k point."
    if (dc%system_tot%nk /= system%nk) stop "DC-LCFO complex: total and fragment nk differ."
    if (.not.allocated(system%vec_k) .or. .not.allocated(system%wtk)) &
      stop "DC-LCFO complex: k-point data is not allocated."
    if (.not.allocated(dc%system_tot%vec_k) .or. .not.allocated(dc%system_tot%wtk)) &
      stop "DC-LCFO complex: total k-point data is not allocated."
    if (maxval(abs(system%vec_k-dc%system_tot%vec_k)) > lcfo_tol .or. &
        maxval(abs(system%wtk-dc%system_tot%wtk)) > lcfo_tol) &
      stop "DC-LCFO complex: fragment and total k-point lists differ."
    if (.not.ieee_is_finite(energy_cut) .or. .not.ieee_is_finite(lambda_cut) .or. &
        lambda_cut <= 0d0) stop "DC-LCFO complex: invalid LCFO cutoff."
    if (any(dc%nxyz_domain < 1) .or. any(dc%nxyz_buffer < 0) .or. &
        any(dc%nxyz_buffer > dc%nxyz_domain)) stop "DC-LCFO complex: invalid domain/buffer."
    if ((dc%id_frag == 0) .neqv. (info%id_rko == 0)) &
      stop "DC-LCFO complex: fragment root and rko root do not agree."
    if (dc%n_frag /= product(num_fragment)) stop "DC-LCFO complex: fragment count mismatch."
    if (dc%nstate_frag /= system%no) stop "DC-LCFO complex: nstate_frag and system%no differ."
    do n=1,3
      if (num_fragment(n) > 1) then
        if (any(abs(system%vec_k(n,1:system%nk)) > lcfo_tol)) &
          stop "DC-LCFO complex: k component along a fragment-split direction must be zero."
      else if (dc%nxyz_buffer(n) /= 0) then
        stop "DC-LCFO complex: buffer is allowed only along split directions."
      end if
    end do
    do ik=1,system%nk
      if (.not.ieee_is_finite(system%wtk(ik))) stop "DC-LCFO complex: non-finite k weight."
    end do
  end subroutine check_lcfo_complex_options

  subroutine dc_lcfo_complex(lg,mg,system,info,stencil,ppg,energy,v_local,spsi,shpsi,sttpsi,srg,dc)
    use communication, only: comm_bcast,comm_irecv,comm_isend,comm_summation,comm_wait_all
    use eigen_subdiag_sub, only: eigen_zheev
    use hamiltonian, only: hpsi
    use ieee_arithmetic, only: ieee_is_finite
    use filesystem, only: get_filehandle
    use math_constants, only: pi
    use salmon_global, only: energy_cut,lambda_cut,sysname, &
         yn_dc_lcfo_diag
    use structures
    implicit none
    type(s_rgrid),         intent(in) :: lg,mg
    type(s_dft_system),    intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    type(s_stencil),       intent(in) :: stencil
    type(s_pp_grid),       intent(in) :: ppg
    type(s_dft_energy),    intent(in) :: energy
    type(s_scalar),        intent(in) :: v_local(system%nspin)
    type(s_orbital),       intent(in) :: spsi
    type(s_orbital)                  :: shpsi,sttpsi
    type(s_sendrecv_grid)             :: srg
    type(s_dcdft)                     :: dc

    type(s_lcfo_complex_halo) :: halo(26)
    integer, allocatable :: id_array(:),n_basis(:,:,:),n_mat(:,:),index_basis(:,:,:,:)
    integer, allocatable :: req_send(:),req_recv(:)
    real(8), allocatable :: esp_tot(:,:,:)
    complex(8), allocatable :: f_basis(:,:,:,:,:),work_basis(:,:,:,:,:)
    complex(8), allocatable :: hf(:,:,:,:,:),work_hf(:,:,:,:,:)
    complex(8), allocatable :: mat_h_local(:,:,:)
    complex(8), allocatable :: hsend(:,:),hmat(:,:),vmat(:,:)
    real(8), allocatable :: eval(:)
    integer :: n_halo,ik,ispin,io,ix,iy,iz,n,nactive,istat,global_status
    integer :: nspin,nk,m,pass,j,ifrag,i,jj,src_frag,tag,ios
    integer, allocatable :: nb(:)
    real(8) :: hvol,herm_err,basis_err,ortho_err,resid_err,norm_h
    logical :: owns_ik

    call check_lcfo_complex_options(system,info,dc)
    if (.not.allocated(spsi%zwf) .or. .not.allocated(sttpsi%zwf) .or. &
        .not.allocated(shpsi%zwf)) stop "DC-LCFO complex: complex orbital scratch is not allocated."
    if (.not.allocated(energy%esp)) stop "DC-LCFO complex: SCF eigenvalues are not allocated."
    if (.not.ieee_is_finite(system%mu)) stop "DC-LCFO complex: chemical potential is not finite."
    if (any(shape(energy%esp) /= [system%no,system%nk,system%nspin])) &
      stop "DC-LCFO complex: SCF eigenvalue shape is inconsistent."
    if (any(.not.ieee_is_finite(energy%esp))) &
      stop "DC-LCFO complex: SCF eigenvalues are not finite."

    if (dc%id_tot == 0) write(*,*) "start DC-LCFO complex"
    hvol = system%hvol
    nspin = system%nspin
    nk = system%nk
    m = dc%nstate_frag

    allocate(id_array(dc%n_frag))
    call init_fragment_roots(id_array)
    call init_halo(n_halo,halo,id_array)
    allocate(n_basis(dc%n_frag,nspin,nk),n_mat(nspin,nk))
    allocate(index_basis(m,dc%n_frag,nspin,nk),esp_tot(dc%nstate_tot,nspin,nk))
    n_basis = 0
    n_mat = 0
    index_basis = 0
    esp_tot = 0d0

    do ik=1,nk
      owns_ik = (info%ik_s <= ik .and. ik <= info%ik_e)
      allocate(f_basis(dc%nxyz_domain(1),dc%nxyz_domain(2),dc%nxyz_domain(3),nspin,m))
      allocate(nb(nspin))
      call build_basis(ik,f_basis,nb,basis_err)
      call collect_basis_dimensions(ik,nb)
      if (any(n_mat(:,ik) < dc%nstate_tot)) &
        stop "DC-LCFO complex: nstate_tot exceeds a k-dependent Hamiltonian dimension."

      sttpsi%zwf = (0d0,0d0)
      shpsi%zwf = (0d0,0d0)
      if (owns_ik) then
        do io=info%io_s,info%io_e
          if (io <= m) then
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              if (ix <= dc%nxyz_domain(1) .and. iy <= dc%nxyz_domain(2) .and. &
                  iz <= dc%nxyz_domain(3)) &
                sttpsi%zwf(ix,iy,iz,1:nspin,io,ik,1) = f_basis(ix,iy,iz,1:nspin,io)
            end do
            end do
            end do
          end if
        end do
      end if

      allocate(hf(lg%num(1),lg%num(2),lg%num(3),nspin,m))
      allocate(work_hf(lg%num(1),lg%num(2),lg%num(3),nspin,m))
      call hpsi(sttpsi,shpsi,info,mg,v_local,system,stencil,srg,ppg)
      work_hf = (0d0,0d0)
      if (owns_ik) then
        do io=info%io_s,info%io_e
          if (io <= m) then
            do ispin=1,nspin
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              work_hf(ix,iy,iz,ispin,io) = shpsi%zwf(ix,iy,iz,ispin,io,ik,1)
            end do
            end do
            end do
            end do
          end if
        end do
      end if
      call comm_summation(work_hf,hf,size(hf),info%icomm_rko)
      deallocate(work_hf)

      allocate(mat_h_local(m,m,nspin))
      call assemble_halo_hamiltonian(ik,f_basis,hf,mat_h_local,n_halo,halo,nactive, &
           req_send,req_recv)
      deallocate(hf,f_basis)

      do ispin=1,nspin
        n = n_mat(ispin,ik)
        allocate(hsend(n,n),hmat(n,n),vmat(n,n),eval(n))
        hsend = (0d0,0d0)
        if (dc%id_frag == 0) then
          ifrag = dc%i_frag
          do io=1,n_basis(ifrag,ispin,ik)
            i = index_basis(io,ifrag,ispin,ik)
            do jj=1,n_basis(ifrag,ispin,ik)
              j = index_basis(jj,ifrag,ispin,ik)
              hsend(i,j) = mat_h_local(io,jj,ispin)
            end do
          end do
          do j=1,n_halo
            src_frag = halo(j)%ifrag_src
            if (.not.allocated(halo(j)%mat_h_local)) cycle
            do io=1,n_basis(ifrag,ispin,ik)
              i = index_basis(io,ifrag,ispin,ik)
              do jj=1,n_basis(src_frag,ispin,ik)
                n = index_basis(jj,src_frag,ispin,ik)
                hsend(n,i) = hsend(n,i) + 0.5d0*halo(j)%mat_h_local(jj,io,ispin)
                hsend(i,n) = hsend(i,n) + &
                     0.5d0*conjg(halo(j)%mat_h_local(jj,io,ispin))
              end do
            end do
          end do
        end if
        call comm_summation(hsend,hmat,size(hmat),dc%icomm_tot)
        herm_err = maxval(abs(hmat-conjg(transpose(hmat)))) / &
             max(1d0,maxval(abs(hmat)))
        istat = 0
        if (.not.ieee_is_finite(herm_err) .or. herm_err > lcfo_tol) istat = 1
        call check_collective_status(istat,"Hermitian Hamiltonian",ik,ispin)
        hmat = 0.5d0*(hmat+conjg(transpose(hmat)))
        do i=1,n
          hmat(i,i) = dcmplx(real(hmat(i,i),8),0d0)
        end do
        vmat = (0d0,0d0)
        istat = 0
        if (dc%id_tot == 0) then
          call eigen_zheev(hmat,eval,vmat,lapack_info=istat)
          if (istat /= 0) istat = 1
          if (istat == 0) call check_eigensystem(hmat,vmat,eval,dc%nstate_tot, &
               ortho_err,resid_err)
          if (istat == 0 .and. (ortho_err > lcfo_tol .or. resid_err > lcfo_tol)) istat = 1
          if (dc%id_tot == 0) write(*,'(a,2i6,4(1x,es12.4))') &
               "complex LCFO k/spin, hermitian/basis/eigen/residual:",ik,ispin, &
               herm_err,basis_err,ortho_err,resid_err
        end if
        call check_collective_status(istat,"ZHEEV",ik,ispin)
        call comm_bcast(eval,dc%icomm_tot,0)
        call comm_bcast(vmat,dc%icomm_tot,0)
        esp_tot(:,ispin,ik) = eval(1:dc%nstate_tot)
        deallocate(hsend,hmat,vmat,eval)
      end do
      deallocate(mat_h_local)
      call deallocate_halo_buffers(n_halo,halo)
      deallocate(nb)
    end do

    call write_complex_eigenvalues(esp_tot,n_basis,n_mat)
    deallocate(esp_tot,index_basis,n_mat,n_basis,id_array)
    if (dc%id_tot == 0) write(*,*) "end DC-LCFO complex"

  contains

    subroutine init_fragment_roots(ids)
      use communication, only: comm_summation
      implicit none
      integer, intent(out) :: ids(:)
      integer :: local(size(ids))
      local = 0
      if (dc%id_frag == 0) local(dc%i_frag) = dc%id_tot + 1
      call comm_summation(local,ids,size(ids),dc%icomm_tot)
      ids = ids - 1
      if (any(ids < 0)) stop "DC-LCFO complex: fragment root map is incomplete."
    end subroutine init_fragment_roots

    subroutine init_halo(nh,halos,ids)
      use salmon_global, only: num_fragment
      implicit none
      integer, intent(out) :: nh
      type(s_lcfo_complex_halo), intent(out) :: halos(:)
      integer, intent(in) :: ids(:)
      integer :: lx,ly,lz,h,ifg,n,d(3),ir1(3),ir2(3),delta(3)
      integer :: nhdir(3)
      nhdir = 0
      do n=1,3
        if (num_fragment(n) > 1) nhdir(n) = 1
      end do
      h = 0
      do lx=-nhdir(1),nhdir(1)
      do ly=-nhdir(2),nhdir(2)
      do lz=-nhdir(3),nhdir(3)
        if (lx == 0 .and. ly == 0 .and. lz == 0) cycle
        h = h + 1
        halos(h)%dvec = [lx,ly,lz]
        halos(h)%id_dst = -1
        halos(h)%id_src = -1
        do ifg=1,dc%n_frag
          ir1 = dc%ixyz_frag(:,ifg)
          ir2 = dc%ixyz_frag(:,dc%i_frag) + halos(h)%dvec*dc%nxyz_domain
          delta = mod(ir1-ir2,dc%lg_tot%num)
          if (all(delta == 0) .and. halos(h)%id_dst < 0) halos(h)%id_dst = ids(ifg)
          ir2 = dc%ixyz_frag(:,dc%i_frag) - halos(h)%dvec*dc%nxyz_domain
          delta = mod(ir1-ir2,dc%lg_tot%num)
          if (all(delta == 0) .and. halos(h)%id_src < 0) then
            halos(h)%id_src = ids(ifg)
            halos(h)%ifrag_src = ifg
          end if
        end do
        if (halos(h)%id_dst < 0 .or. halos(h)%id_src < 0) &
          stop "DC-LCFO complex: halo neighbor was not found."
        do n=1,3
          select case(halos(h)%dvec(n))
          case(0)
            halos(h)%length(n) = dc%nxyz_domain(n)
            halos(h)%dsp_send(n) = 0
            halos(h)%dsp_recv(n) = 0
          case(1)
            halos(h)%length(n) = dc%nxyz_buffer(n)
            halos(h)%dsp_send(n) = dc%nxyz_domain(n)-dc%nxyz_buffer(n)
            halos(h)%dsp_recv(n) = dc%nxyz_domain(n)+dc%nxyz_buffer(n)
          case(-1)
            halos(h)%length(n) = dc%nxyz_buffer(n)
            halos(h)%dsp_send(n) = 0
            halos(h)%dsp_recv(n) = dc%nxyz_domain(n)
          end select
        end do
      end do
      end do
      end do
      nh = h
    end subroutine init_halo

    subroutine build_basis(ik0,basis,nb0,basis_err0)
      use communication, only: comm_bcast,comm_summation
      use eigen_subdiag_sub, only: eigen_zheev
      implicit none
      integer, intent(in) :: ik0
      complex(8), intent(out) :: basis(:,:,:,:,:)
      integer, intent(out) :: nb0(:)
      real(8), intent(out) :: basis_err0
      complex(8), allocatable :: local(:,:,:,:,:),work(:,:,:,:,:),smat(:,:,:),umat(:,:,:)
      real(8), allocatable :: lambda(:,:)
      integer :: io0,jo0,isp0,ix0,iy0,iz0,ieig,ib,pass0,status0
      complex(8) :: coeff
      real(8) :: norm2

      allocate(local(size(basis,1),size(basis,2),size(basis,3),nspin,m))
      allocate(work(size(basis,1),size(basis,2),size(basis,3),nspin,m))
      allocate(smat(m,m,nspin),umat(m,m,nspin),lambda(m,nspin))
      local = (0d0,0d0)
      do io0=info%io_s,info%io_e
        if (io0 <= m .and. info%ik_s <= ik0 .and. ik0 <= info%ik_e) then
          do isp0=1,nspin
          do iz0=mg%is(3),mg%ie(3)
          do iy0=mg%is(2),mg%ie(2)
          do ix0=mg%is(1),mg%ie(1)
            if (ix0 <= dc%nxyz_domain(1) .and. iy0 <= dc%nxyz_domain(2) .and. &
                iz0 <= dc%nxyz_domain(3)) then
              if (energy%esp(io0,ik0,isp0)-system%mu < energy_cut) &
                local(ix0,iy0,iz0,isp0,io0) = spsi%zwf(ix0,iy0,iz0,isp0,io0,ik0,1)
            end if
          end do
          end do
          end do
          end do
        end if
      end do
      call comm_summation(local,basis,size(basis),info%icomm_rko)
      deallocate(local)
      smat = (0d0,0d0)
      umat = (0d0,0d0)
      lambda = 0d0
      status0 = 0
      if (dc%id_frag == 0) then
        status0 = 0
        do isp0=1,nspin
        do io0=1,m
        do jo0=1,m
          smat(io0,jo0,isp0) = hvol*sum(conjg(basis(:,:,:,isp0,io0))*basis(:,:,:,isp0,jo0))
        end do
        end do
        end do
        status0 = 0
        do isp0=1,nspin
          call eigen_zheev(smat(:,:,isp0),lambda(:,isp0),umat(:,:,isp0),lapack_info=status0)
          if (status0 /= 0) exit
          if (any(.not.ieee_is_finite(lambda(:,isp0))) .or. &
              minval(lambda(:,isp0)) < -lcfo_tol*max(1d0,maxval(abs(lambda(:,isp0))))) then
            status0 = 1
            exit
          end if
        end do
      end if
      call check_collective_status(status0,"fragment overlap ZHEEV",ik0,0)
      work = basis
      basis = (0d0,0d0)
      nb0 = 0
      if (dc%id_frag == 0) then
        do isp0=1,nspin
          ib = 0
          do ieig=m,1,-1
            if (lambda(ieig,isp0) > lambda_cut) then
              ib = ib + 1
              do jo0=1,m
                basis(:,:,:,isp0,ib) = basis(:,:,:,isp0,ib) + &
                     work(:,:,:,isp0,jo0)*umat(jo0,ieig,isp0)/sqrt(lambda(ieig,isp0))
              end do
            end if
          end do
          nb0(isp0) = ib
        end do
        do pass0=1,2
          do isp0=1,nspin
            ib = count(lambda(:,isp0) > lambda_cut)
            do io0=1,ib
              do jo0=1,io0-1
                coeff = hvol*sum(conjg(basis(:,:,:,isp0,jo0))*basis(:,:,:,isp0,io0))
                basis(:,:,:,isp0,io0) = basis(:,:,:,isp0,io0) - &
                     basis(:,:,:,isp0,jo0)*coeff
              end do
              norm2 = hvol*sum(abs(basis(:,:,:,isp0,io0))**2)
              if (.not.ieee_is_finite(norm2) .or. norm2 <= 100d0*epsilon(1d0)) then
                status0 = 1
                exit
              end if
              basis(:,:,:,isp0,io0) = basis(:,:,:,isp0,io0)/sqrt(norm2)
            end do
          end do
        end do
      end if
      call check_collective_status(status0,"Gram-Schmidt",ik0,0)
      basis_err0 = 0d0
      if (dc%id_frag == 0) then
        do isp0=1,nspin
          do io0=1,nb0(isp0)
          do jo0=1,nb0(isp0)
            coeff = hvol*sum(conjg(basis(:,:,:,isp0,io0))*basis(:,:,:,isp0,jo0))
            if (io0 == jo0) coeff = coeff-(1d0,0d0)
            basis_err0 = max(basis_err0,abs(coeff))
          end do
          end do
        end do
      end if
      call comm_bcast(basis_err0,info%icomm_rko,0)
      status0 = 0
      if (basis_err0 > lcfo_tol .or. .not.ieee_is_finite(basis_err0)) status0 = 1
      call check_collective_status(status0,"fragment basis orthogonality",ik0,0)
      call comm_bcast(nb0,info%icomm_rko,0)
      call comm_bcast(basis,info%icomm_rko,0)
      deallocate(lambda,umat,smat,work)
    end subroutine build_basis

    subroutine collect_basis_dimensions(ik0,nb0)
      use communication, only: comm_summation
      implicit none
      integer, intent(in) :: ik0,nb0(:)
      integer :: local_n(dc%n_frag,nspin),isp0,ifg0,io0,idx0
      local_n = 0
      if (dc%id_frag == 0) local_n(dc%i_frag,:) = nb0
      call comm_summation(local_n,n_basis(:,:,ik0),dc%n_frag*nspin,dc%icomm_tot)
      do isp0=1,nspin
        idx0 = 0
        do ifg0=1,dc%n_frag
          do io0=1,n_basis(ifg0,isp0,ik0)
            idx0 = idx0 + 1
            index_basis(io0,ifg0,isp0,ik0) = idx0
          end do
        end do
        n_mat(isp0,ik0) = idx0
      end do
    end subroutine collect_basis_dimensions

    subroutine assemble_halo_hamiltonian(ik0,basis,hf0,diag_h,nh,halos,nact,rs,rr)
      use communication, only: comm_irecv,comm_isend,comm_wait_all
      implicit none
      integer, intent(in) :: ik0,nh
      complex(8), intent(in) :: basis(:,:,:,:,:),hf0(:,:,:,:,:)
      complex(8), intent(out) :: diag_h(:,:,:)
      type(s_lcfo_complex_halo), intent(inout) :: halos(:)
      integer, intent(out) :: nact
      integer, allocatable, intent(out) :: rs(:),rr(:)
      integer :: h,ia,ib,ic,isp0,io0,jo0,jj0,k0,tag_send,tag_recv
      integer :: l(3),d(3),nreq

      diag_h = (0d0,0d0)
      if (dc%id_frag == 0) then
        do isp0=1,nspin
        do io0=1,n_basis(dc%i_frag,isp0,ik0)
        do jo0=1,n_basis(dc%i_frag,isp0,ik0)
          diag_h(io0,jo0,isp0) = hvol*sum(conjg(basis(:,:,:,isp0,io0))* &
               hf0(1:dc%nxyz_domain(1),1:dc%nxyz_domain(2),1:dc%nxyz_domain(3),isp0,jo0))
        end do
        end do
        end do
      end if

      nact = 0
      do h=1,nh
        if (all(halos(h)%length > 0)) nact = nact + 1
      end do
      allocate(rs(max(1,nact)),rr(max(1,nact)))
      nreq = 0
      if (dc%id_frag == 0) then
        do h=1,nh
          if (.not.all(halos(h)%length > 0)) cycle
          nreq = nreq + 1
          l = halos(h)%length
          d = halos(h)%dsp_send
          allocate(halos(h)%buf_send(l(1),l(2),l(3),nspin,m), &
               halos(h)%buf_recv(l(1),l(2),l(3),nspin,m))
          do ic=1,l(3)
          do ib=1,l(2)
          do ia=1,l(1)
            halos(h)%buf_send(ia,ib,ic,:,:) = basis(d(1)+ia,d(2)+ib,d(3)+ic,:,:)
          end do
          end do
          end do
          tag_send = halo_tag(dc%i_frag,halos(h)%dvec)
          tag_recv = halo_tag(halos(h)%ifrag_src,halos(h)%dvec)
          rr(nreq) = comm_irecv(halos(h)%buf_recv,halos(h)%id_src,tag_recv,dc%icomm_tot)
          rs(nreq) = comm_isend(halos(h)%buf_send,halos(h)%id_dst,tag_send,dc%icomm_tot)
        end do
        if (nreq > 0) then
          call comm_wait_all(rr(1:nreq))
          call comm_wait_all(rs(1:nreq))
        end if
        nreq = 0
        do h=1,nh
          if (.not.all(halos(h)%length > 0)) cycle
          nreq = nreq + 1
          l = halos(h)%length
          d = halos(h)%dsp_recv
          allocate(halos(h)%mat_h_local(m,m,nspin))
          halos(h)%mat_h_local = (0d0,0d0)
          do isp0=1,nspin
          do jo0=1,m
          do io0=1,m
            do ic=1,l(3)
            do ib=1,l(2)
            do ia=1,l(1)
              halos(h)%mat_h_local(jo0,io0,isp0) = halos(h)%mat_h_local(jo0,io0,isp0) + &
                   hvol*conjg(halos(h)%buf_recv(ia,ib,ic,isp0,jo0))* &
                   hf0(d(1)+ia,d(2)+ib,d(3)+ic,isp0,io0)
            end do
            end do
            end do
          end do
          end do
          end do
          deallocate(halos(h)%buf_send,halos(h)%buf_recv)
        end do
      end if
      deallocate(rs,rr)
    end subroutine assemble_halo_hamiltonian

    integer function halo_tag(ifg,dvec) result(tag0)
      integer, intent(in) :: ifg,dvec(3)
      tag0 = (ifg-1)*27 + (dvec(1)+1)*9 + (dvec(2)+1)*3 + dvec(3) + 2
    end function halo_tag

    subroutine deallocate_halo_buffers(nh,halos)
      implicit none
      integer, intent(in) :: nh
      type(s_lcfo_complex_halo), intent(inout) :: halos(:)
      integer :: h
      do h=1,nh
        if (allocated(halos(h)%buf_send)) deallocate(halos(h)%buf_send)
        if (allocated(halos(h)%buf_recv)) deallocate(halos(h)%buf_recv)
        if (allocated(halos(h)%mat_h_local)) deallocate(halos(h)%mat_h_local)
      end do
    end subroutine deallocate_halo_buffers

    subroutine check_eigensystem(h0,v0,e0,nt,oe,re)
      implicit none
      complex(8), intent(in) :: h0(:,:),v0(:,:)
      real(8), intent(in) :: e0(:)
      integer, intent(in) :: nt
      real(8), intent(out) :: oe,re
      complex(8), allocatable :: gram(:,:),res(:)
      integer :: row,col,k
      complex(8) :: dot
      real(8) :: denom
      allocate(gram(size(v0,2),size(v0,2)))
      gram = (0d0,0d0)
      do row=1,size(v0,2)
      do col=1,size(v0,2)
        dot = (0d0,0d0)
        do k=1,size(v0,1)
          dot = dot + conjg(v0(k,row))*v0(k,col)
        end do
        gram(row,col) = dot
      end do
      end do
      do col=1,size(gram,1)
        gram(col,col) = gram(col,col)-(1d0,0d0)
      end do
      oe = maxval(abs(gram))
      re = 0d0
      denom = max(1d0,sqrt(sum(abs(h0)**2)))
      allocate(res(size(h0,1)))
      do col=1,nt
        do row=1,size(h0,1)
          dot = (0d0,0d0)
          do k=1,size(h0,2)
            dot = dot + h0(row,k)*v0(k,col)
          end do
          res(row) = dot-e0(col)*v0(row,col)
        end do
        re = max(re,sqrt(sum(abs(res)**2))/max(denom,abs(e0(col)),1d0))
      end do
      deallocate(gram,res)
      if (.not.ieee_is_finite(oe) .or. .not.ieee_is_finite(re) .or. &
          oe > lcfo_tol .or. re > lcfo_tol) then
        ! The caller converts the failed diagnostic into a collective status.
        re = max(re,2d0*lcfo_tol)
      end if
    end subroutine check_eigensystem

    subroutine check_collective_status(local_status,stage,ik0,isp0)
      use communication, only: comm_summation
      implicit none
      integer, intent(in) :: local_status,ik0,isp0
      character(*), intent(in) :: stage
      integer :: total_status
      call comm_summation(local_status,total_status,dc%icomm_tot)
      if (total_status /= 0) then
        if (dc%id_tot == 0) write(*,'(a,2i6,1x,a)') &
             "DC-LCFO complex failure at k/spin:",ik0,isp0,trim(stage)
        stop "DC-LCFO complex: collective failure."
      end if
    end subroutine check_collective_status

    subroutine write_complex_eigenvalues(values,nbasis,nmat0)
      use communication, only: comm_summation
      use filesystem, only: get_filehandle
      implicit none
      real(8), intent(in) :: values(:,:,:)
      integer, intent(in) :: nbasis(:,:,:),nmat0(:,:)
      integer :: iu,ik0,isp0,ib0,ifg0,local_status,total_status
      character(256) :: filename
      real(8) :: qred(3)
      local_status = 0
      if (dc%id_tot == 0) then
        iu = get_filehandle()
        filename = trim(dc%base_directory)//trim(sysname)//'_lcfo_complex_eigen.data'
        open(iu,file=filename,status='replace',form='formatted',iostat=local_status)
        if (local_status == 0) then
          write(iu,'(a)') '# DC-LCFO complex scalar eigenvalues; format_version=1'
          write(iu,'(a,1x,i0,1x,a,1x,i0,1x,a,1x,i0)') '# nk=',nk,'nspin=',nspin, &
               'nstate=',dc%nstate_tot
          write(iu,'(a)') '# columns: ik ispin iband eigenvalue_Ha'
          do ik0=1,nk
            qred = [dot_product(dc%system_tot%primitive_a(:,1),dc%system_tot%vec_k(:,ik0)), &
                 dot_product(dc%system_tot%primitive_a(:,2),dc%system_tot%vec_k(:,ik0)), &
                 dot_product(dc%system_tot%primitive_a(:,3),dc%system_tot%vec_k(:,ik0))]/(2d0*pi)
            write(iu,'(a,i0,a,3(1x,es26.16e3),a,3(1x,es26.16e3),a,es26.16e3)') &
                 '# k=',ik0,' reduced=',qred,' cartesian=',dc%system_tot%vec_k(:,ik0), &
                 ' weight=',dc%system_tot%wtk(ik0)
            write(iu,'(a,i0,a,*(1x,i0))') '# k=',ik0,' nmat=',nmat0(:,ik0)
            do ifg0=1,dc%n_frag
              write(iu,'(a,i0,a,i0,a,*(1x,i0))') '# k=',ik0,' fragment=',ifg0, &
                   ' n_basis=',nbasis(ifg0,:,ik0)
            end do
            do isp0=1,nspin
              do ib0=1,dc%nstate_tot
                write(iu,'(3(i0,1x),es26.16e3)') ik0,isp0,ib0,values(ib0,isp0,ik0)
              end do
            end do
          end do
          write(iu,'(a)') '# end DC-LCFO complex eigenvalues'
          close(iu)
        end if
      end if
      call comm_summation(local_status,total_status,dc%icomm_tot)
      if (total_status /= 0) stop "DC-LCFO complex: failed to write eigenvalue output."
    end subroutine write_complex_eigenvalues

  end subroutine dc_lcfo_complex

end module lcfo_complex
