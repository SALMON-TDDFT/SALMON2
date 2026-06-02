!
!  Copyright 2019-2024 SALMON developers
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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90-------100-------110-------120-------130

#include "config.h"
module lcfo_soi
  use dc_fragment_geometry, only: get_fragment_domain
  implicit none

  private
  public :: dc_lcfo_soi

  character(32),parameter :: binfile_wf_soi = "wavefunctions_soi.bin", &
  &                          binfile_rg_soi = "rgrid_index_soi.bin", &
  &                          binfile_bf_soi = "basis_functions_soi.bin", &
  &                          binfile_bfb_soi = "basis_functions_buffer_soi.bin", &
  &                          binfile_hl_soi = "hamiltonian_local_soi.bin"
  integer, parameter :: basis_buffer_magic_soi = -22022213
  integer, parameter :: basis_buffer_version_soi = 1

contains

  subroutine dc_lcfo_soi(lg,mg,system,info,stencil,ppg,energy,v_local,spsi,shpsi,sttpsi,srg,dc)
    use communication, only: comm_summation
    use salmon_global, only: yn_dc_lcfo_diag, yn_spinorbit
    use structures
    implicit none
    type(s_rgrid),        intent(in) :: lg,mg
    type(s_dft_system),   intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_stencil),      intent(in) :: stencil
    type(s_pp_grid),      intent(in) :: ppg
    type(s_dft_energy),   intent(in) :: energy
    type(s_scalar),       intent(in) :: v_local(system%nspin)
    type(s_orbital),      intent(in) :: spsi
    type(s_orbital)                  :: shpsi,sttpsi
    type(s_sendrecv_grid)            :: srg
    type(s_dcdft)                    :: dc
    type halo_info
      integer :: id_src,id_dst,ifrag_src,dvec(3),length(3),dsp_send(3),dsp_recv(3)
      complex(8),allocatable :: buf_send(:,:,:,:,:),buf_recv(:,:,:,:,:),mat_H_local(:,:,:)
    end type halo_info
    type(halo_info) :: halo(26)
    integer :: nspin,n_halo
    integer :: id_array(dc%n_frag)
    integer :: n_basis(dc%n_frag,system%nspin), n_mat(system%nspin)
    integer :: index_basis(dc%nstate_frag,dc%n_frag,system%nspin)
    real(8) :: hvol
    complex(8),allocatable :: f_basis(:,:,:,:,:),hf(:,:,:,:,:),wrk_array(:,:,:,:,:) &
    & ,mat_H_local(:,:,:),coef_wf(:,:,:),basis_transform(:,:,:)
    real(8),allocatable :: esp_tot(:,:)
    logical :: use_spinor_basis
    integer :: i,j,n,ix,iy,iz,io,jo,ispin,ifrag,jfrag,i_halo

    if(.not. allocated(spsi%zwf)) stop "DC-LCFO(SOI): spsi%zwf is required."

    if(dc%id_tot==0) write(*,*) "start DC-LCFO(SOI)"
    hvol = system%hvol
    nspin = system%nspin
    use_spinor_basis = (yn_spinorbit=='y' .and. nspin>1)

    call init_lcfo
    call calc_basis
    call hpsi_basis
    if(dc%id_tot==0) write(*,*) "basis functions operation: done"

    call calc_hamiltonian_matrix
    if(dc%id_tot==0) write(*,*) "Hamiltonian matrix: done"

    if(yn_dc_lcfo_diag=='y') then
      allocate(esp_tot(maxval(n_mat),nspin))
      if(dc%id_frag==0) allocate(coef_wf(dc%nstate_frag,dc%nstate_tot,nspin))
      call diag_hermitian
      if(dc%id_tot==0) write(*,*) "diagonalization: done"
    end if

    call output

    if(allocated(coef_wf)) deallocate(coef_wf)
    if(allocated(f_basis)) deallocate(f_basis)
    if(allocated(esp_tot)) deallocate(esp_tot)
    if(allocated(mat_H_local)) deallocate(mat_H_local)
    do i=1,n_halo
      if(allocated(halo(i)%mat_H_local)) deallocate(halo(i)%mat_H_local)
    end do
    if(dc%id_tot==0) write(*,*) "end DC-LCFO(SOI)"

  contains

    subroutine init_lcfo
      use salmon_global, only: num_fragment
      implicit none
      integer :: lx,ly,lz
      integer :: nxyz_domain(3)
      integer,dimension(3) :: nh,ir1,ir2,d
      integer :: id_tmp(dc%n_frag)

      id_tmp = 0
      if(dc%id_frag==0) id_tmp(dc%i_frag) = dc%id_tot + 1
      call comm_summation(id_tmp,id_array,dc%n_frag,dc%icomm_tot)
      id_array = id_array - 1
      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)

      nh = 0
      do n=1,3
        if(dc%nxyz_buffer(n) > nxyz_domain(n)) stop "DC-LCFO(SOI): buffer > domain"
        if(num_fragment(n) > 1) nh(n) = 1
      end do

      i = 0
      do lx=-nh(1),nh(1)
      do ly=-nh(2),nh(2)
      do lz=-nh(3),nh(3)
        if(lx==0 .and. ly==0 .and. lz==0) cycle
        i = i + 1
        halo(i)%dvec(1:3) = [lx, ly, lz]
        halo(i)%id_dst = -1
        halo(i)%id_src = -1
        do ifrag=1,dc%n_frag
          ir1(1:3) = dc%ixyz_frag(1:3,ifrag)
          ir2(1:3) = dc%ixyz_frag(1:3,dc%i_frag) + halo(i)%dvec(1:3)*nxyz_domain(1:3)
          d(1:3) = mod(ir1(1:3) - ir2(1:3), dc%lg_tot%num(1:3))
          if(d(1)==0 .and. d(2)==0 .and. d(3)==0 .and. halo(i)%id_dst < 0) then
            halo(i)%id_dst = id_array(ifrag)
          end if

          ir2(1:3) = dc%ixyz_frag(1:3,dc%i_frag) - halo(i)%dvec(1:3)*nxyz_domain(1:3)
          d(1:3) = mod(ir1(1:3) - ir2(1:3), dc%lg_tot%num(1:3))
          if(d(1)==0 .and. d(2)==0 .and. d(3)==0 .and. halo(i)%id_src < 0) then
            halo(i)%id_src = id_array(ifrag)
            halo(i)%ifrag_src = ifrag
          end if
        end do
        if(halo(i)%id_dst < 0 .or. halo(i)%id_src < 0) stop "DC-LCFO(SOI): dst, src"
      end do
      end do
      end do
      n_halo = i

      do i=1,n_halo
        do n=1,3
          select case (halo(i)%dvec(n))
          case(0)
            halo(i)%length(n) = nxyz_domain(n)
            halo(i)%dsp_send(n) = 0
            halo(i)%dsp_recv(n) = 0
          case(1)
            halo(i)%length(n) = dc%nxyz_buffer(n)
            halo(i)%dsp_send(n) = nxyz_domain(n) - dc%nxyz_buffer(n)
            halo(i)%dsp_recv(n) = nxyz_domain(n) + dc%nxyz_buffer(n)
          case(-1)
            halo(i)%length(n) = dc%nxyz_buffer(n)
            halo(i)%dsp_send(n) = 0
            halo(i)%dsp_recv(n) = nxyz_domain(n)
          end select
        end do
      end do
    end subroutine init_lcfo

    subroutine calc_basis
      use eigen_subdiag_sub, only: eigen_zheev
      use salmon_global, only: energy_cut,lambda_cut
      implicit none
      integer :: nb(nspin),itmp(dc%n_frag,nspin)
      integer :: nxyz_domain(3)
      logical :: include_orb
      complex(8),dimension(dc%nstate_frag,dc%nstate_frag,system%nspin) :: mat_S,mat_U
      complex(8),dimension(dc%nstate_frag,dc%nstate_frag) :: mat_S_spinor,mat_U_spinor
      real(8),dimension(dc%nstate_frag,system%nspin) :: lambda
      real(8),dimension(dc%nstate_frag) :: lambda_spinor
      complex(8) :: ztmp
      real(8) :: norm_basis
      
      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)

      allocate(f_basis  (nxyz_domain(1),nxyz_domain(2),nxyz_domain(3),nspin,dc%nstate_frag))
      allocate(wrk_array(nxyz_domain(1),nxyz_domain(2),nxyz_domain(3),nspin,dc%nstate_frag))
      if(.not. allocated(basis_transform)) allocate(basis_transform(dc%nstate_frag,dc%nstate_frag,nspin))
      basis_transform = (0d0,0d0)

      wrk_array = (0d0,0d0)
      if(use_spinor_basis) then
        do io=info%io_s,info%io_e
          include_orb = any(energy%esp(io,1,1:nspin) - system%mu < energy_cut)
          if(.not. include_orb) cycle
          do ispin=1,nspin
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            if(ix <= nxyz_domain(1) .and. iy <= nxyz_domain(2) .and. iz <= nxyz_domain(3)) then
              wrk_array(ix,iy,iz,ispin,io) = spsi%zwf(ix,iy,iz,ispin,io,1,1)
            end if
          end do
          end do
          end do
          end do
        end do
      else
        do io=info%io_s,info%io_e
        do ispin=1,nspin
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          if(ix <= nxyz_domain(1) .and. iy <= nxyz_domain(2) .and. iz <= nxyz_domain(3) &
          & .and. energy%esp(io,1,ispin) - system%mu < energy_cut) then
            wrk_array(ix,iy,iz,ispin,io) = spsi%zwf(ix,iy,iz,ispin,io,1,1)
          end if
        end do
        end do
        end do
        end do
        end do
      end if
      call comm_summation(wrk_array,f_basis,product(nxyz_domain)*nspin*dc%nstate_frag,info%icomm_rko)

      if(use_spinor_basis) then
        do io=1,dc%nstate_frag
        do jo=1,dc%nstate_frag
          ztmp = (0d0,0d0)
          do ispin=1,nspin
            ztmp = ztmp + sum(conjg(f_basis(:,:,:,ispin,io))*f_basis(:,:,:,ispin,jo)) * hvol
          end do
          mat_S_spinor(io,jo) = ztmp
        end do
        end do
        call eigen_zheev(mat_S_spinor,lambda_spinor,mat_U_spinor)
      else
        do ispin=1,nspin
        do io=1,dc%nstate_frag
        do jo=1,dc%nstate_frag
          ztmp = sum(conjg(f_basis(:,:,:,ispin,io))*f_basis(:,:,:,ispin,jo)) * hvol
          mat_S(io,jo,ispin) = ztmp
        end do
        end do
        end do
        do ispin=1,nspin
          call eigen_zheev(mat_S(:,:,ispin),lambda(:,ispin),mat_U(:,:,ispin))
        end do
      end if

      wrk_array = f_basis
      f_basis = (0d0,0d0)
      if(use_spinor_basis) then
        i = 0
        do io=dc%nstate_frag,1,-1
          if(lambda_spinor(io) > lambda_cut) then
            i = i + 1
            do ispin=1,nspin
            do jo=1,dc%nstate_frag
              f_basis(:,:,:,ispin,i) = f_basis(:,:,:,ispin,i) &
              & + wrk_array(:,:,:,ispin,jo) * mat_U_spinor(jo,io) / sqrt(lambda_spinor(io))
              basis_transform(jo,i,ispin) = mat_U_spinor(jo,io) / sqrt(lambda_spinor(io))
            end do
            end do
          end if
        end do
        nb(1:nspin) = i
      else
        do ispin=1,nspin
          i = 0
          do io=dc%nstate_frag,1,-1
            if(lambda(io,ispin) > lambda_cut) then
              i = i + 1
              do jo=1,dc%nstate_frag
                f_basis(:,:,:,ispin,i) = f_basis(:,:,:,ispin,i) &
                & + wrk_array(:,:,:,ispin,jo) * mat_U(jo,io,ispin) / sqrt(lambda(io,ispin))
                basis_transform(jo,i,ispin) = mat_U(jo,io,ispin) / sqrt(lambda(io,ispin))
              end do
            end if
          end do
          nb(ispin) = i
        end do
      end if

      wrk_array = f_basis
      if(use_spinor_basis) then
        do io=1,nb(1)
          do jo=1,io-1
            ztmp = (0d0,0d0)
            do ispin=1,nspin
              ztmp = ztmp + sum(conjg(f_basis(:,:,:,ispin,jo))*wrk_array(:,:,:,ispin,io)) * hvol
            end do
            ztmp = ztmp / sum(conjg(f_basis(:,:,:,1:nspin,jo))*f_basis(:,:,:,1:nspin,jo)) / hvol
            do ispin=1,nspin
              wrk_array(:,:,:,ispin,io) = wrk_array(:,:,:,ispin,io) &
              & - f_basis(:,:,:,ispin,jo) * ztmp
              basis_transform(:,io,ispin) = basis_transform(:,io,ispin) &
              & - basis_transform(:,jo,ispin) * ztmp
            end do
          end do
          ztmp = (0d0,0d0)
          do ispin=1,nspin
            ztmp = ztmp + sum(conjg(wrk_array(:,:,:,ispin,io))*wrk_array(:,:,:,ispin,io)) * hvol
          end do
          do ispin=1,nspin
            wrk_array(:,:,:,ispin,io) = wrk_array(:,:,:,ispin,io) / sqrt(real(ztmp))
            basis_transform(:,io,ispin) = basis_transform(:,io,ispin) / sqrt(real(ztmp))
          end do
        end do
      else
        do ispin=1,nspin
          do io=1,nb(ispin)
            do jo=1,io-1
              ztmp = sum(conjg(f_basis(:,:,:,ispin,jo))*wrk_array(:,:,:,ispin,io)) * hvol
              wrk_array(:,:,:,ispin,io) = wrk_array(:,:,:,ispin,io) &
              & - f_basis(:,:,:,ispin,jo) * ztmp / (sum(conjg(f_basis(:,:,:,ispin,jo))*f_basis(:,:,:,ispin,jo)) * hvol)
              basis_transform(:,io,ispin) = basis_transform(:,io,ispin) &
              & - basis_transform(:,jo,ispin) * ztmp / &
              &   (sum(conjg(f_basis(:,:,:,ispin,jo))*f_basis(:,:,:,ispin,jo)) * hvol)
            end do
            norm_basis = sqrt(real(sum(conjg(wrk_array(:,:,:,ispin,io))*wrk_array(:,:,:,ispin,io))) * hvol)
            wrk_array(:,:,:,ispin,io) = wrk_array(:,:,:,ispin,io) / norm_basis
            basis_transform(:,io,ispin) = basis_transform(:,io,ispin) / norm_basis
          end do
        end do
      end if
      f_basis = wrk_array

      sttpsi%zwf = (0d0,0d0)
      do io=info%io_s,info%io_e
      do ispin=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        if(ix <= dc%nxyz_domain(1) .and. iy <= dc%nxyz_domain(2) .and. iz <= dc%nxyz_domain(3)) then
          sttpsi%zwf(ix,iy,iz,ispin,io,1,1) = f_basis(ix,iy,iz,ispin,io)
        end if
      end do
      end do
      end do
      end do
      end do

      itmp = 0
      if(dc%id_frag==0) itmp(dc%i_frag,1:nspin) = nb(1:nspin)
      call comm_summation(itmp,n_basis,dc%n_frag*nspin,dc%icomm_tot)
      index_basis = 0
      do ispin=1,nspin
        i = 0
        do ifrag=1,dc%n_frag
          do io=1,n_basis(ifrag,ispin)
            i = i + 1
            index_basis(io,ifrag,ispin) = i
          end do
        end do
        n_mat(ispin) = i
      end do

      deallocate(wrk_array)
    end subroutine calc_basis

    subroutine hpsi_basis
      use hamiltonian, only: hpsi
      implicit none

      allocate(hf       (lg%num(1),lg%num(2),lg%num(3),nspin,dc%nstate_frag))
      allocate(wrk_array(lg%num(1),lg%num(2),lg%num(3),nspin,dc%nstate_frag))

      call hpsi(sttpsi,shpsi,info,mg,v_local,system,stencil,srg,ppg)

      wrk_array = (0d0,0d0)
      do io=info%io_s,info%io_e
      do ispin=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        wrk_array(ix,iy,iz,ispin,io) = shpsi%zwf(ix,iy,iz,ispin,io,1,1)
      end do
      end do
      end do
      end do
      end do
      call comm_summation(wrk_array,hf,product(lg%num)*nspin*dc%nstate_frag,info%icomm_rko)

      deallocate(wrk_array)
    end subroutine hpsi_basis

    subroutine calc_hamiltonian_matrix
      use communication, only: comm_isend, comm_irecv, comm_wait_all
      implicit none
      integer :: d(3),l(3)
      integer :: itag_send,itag_recv
      integer,dimension(n_halo) :: ireq_send,ireq_recv
      complex(8) :: ztmp

      allocate(mat_H_local(dc%nstate_frag,dc%nstate_frag,nspin))
      mat_H_local = (0d0,0d0)
      l = dc%nxyz_domain
      if(use_spinor_basis) then
        do io=1,n_basis(dc%i_frag,1)
        do jo=1,n_basis(dc%i_frag,1)
          ztmp = (0d0,0d0)
          do ispin=1,nspin
            ztmp = ztmp + sum(conjg(f_basis(1:l(1),1:l(2),1:l(3),ispin,io)) &
            & * hf(1:l(1),1:l(2),1:l(3),ispin,jo)) * hvol
          end do
          mat_H_local(io,jo,1:nspin) = ztmp
        end do
        end do
      else
        do ispin=1,nspin
        do io=1,n_basis(dc%i_frag,ispin)
        do jo=1,n_basis(dc%i_frag,ispin)
          mat_H_local(io,jo,ispin) = sum(conjg(f_basis(1:l(1),1:l(2),1:l(3),ispin,io)) &
          & * hf(1:l(1),1:l(2),1:l(3),ispin,jo)) * hvol
        end do
        end do
        end do
      end if

      if(dc%id_frag==0) then
        do i_halo=1,n_halo
          l = halo(i_halo)%length
          d = halo(i_halo)%dsp_send
          allocate(halo(i_halo)%buf_send(l(1),l(2),l(3),nspin,dc%nstate_frag))
          allocate(halo(i_halo)%buf_recv(l(1),l(2),l(3),nspin,dc%nstate_frag))
          do iz=1,l(3)
          do iy=1,l(2)
          do ix=1,l(1)
            halo(i_halo)%buf_send(ix,iy,iz,:,:) = f_basis(d(1)+ix,d(2)+iy,d(3)+iz,:,:)
          end do
          end do
          end do

          itag_send = dc%i_frag
          ireq_send(i_halo) = comm_isend(halo(i_halo)%buf_send,halo(i_halo)%id_dst,itag_send,dc%icomm_tot)
          itag_recv = halo(i_halo)%ifrag_src
          ireq_recv(i_halo) = comm_irecv(halo(i_halo)%buf_recv,halo(i_halo)%id_src,itag_recv,dc%icomm_tot)
        end do
        call comm_wait_all(ireq_recv)
        call comm_wait_all(ireq_send)

        do i_halo=1,n_halo
          l = halo(i_halo)%length
          d = halo(i_halo)%dsp_recv
          allocate(halo(i_halo)%mat_H_local(dc%nstate_frag,dc%nstate_frag,nspin))
          halo(i_halo)%mat_H_local = (0d0,0d0)
          if(use_spinor_basis) then
            do io=1,dc%nstate_frag
            do jo=1,dc%nstate_frag
              ztmp = (0d0,0d0)
              do ispin=1,nspin
                do iz=1,l(3)
                do iy=1,l(2)
                do ix=1,l(1)
                  ztmp = ztmp + conjg(halo(i_halo)%buf_recv(ix,iy,iz,ispin,io)) &
                  & * hf(d(1)+ix,d(2)+iy,d(3)+iz,ispin,jo) * hvol
                end do
                end do
                end do
              end do
              halo(i_halo)%mat_H_local(io,jo,1:nspin) = ztmp
            end do
            end do
          else
            do ispin=1,nspin
            do io=1,dc%nstate_frag
            do jo=1,dc%nstate_frag
              do iz=1,l(3)
              do iy=1,l(2)
              do ix=1,l(1)
                halo(i_halo)%mat_H_local(io,jo,ispin) = halo(i_halo)%mat_H_local(io,jo,ispin) &
                & + conjg(halo(i_halo)%buf_recv(ix,iy,iz,ispin,io)) * hf(d(1)+ix,d(2)+iy,d(3)+iz,ispin,jo) * hvol
              end do
              end do
              end do
            end do
            end do
            end do
          end if
          deallocate(halo(i_halo)%buf_send,halo(i_halo)%buf_recv)
          if(dc%id_tot==0) write(*,*) "Halo communication #",i_halo,": done"
        end do
      end if
      deallocate(hf)
    end subroutine calc_hamiltonian_matrix

    subroutine diag_hermitian
      use eigen_subdiag_sub, only: eigen_zheev
#ifdef USE_SCALAPACK
      use eigen_subdiag_sub, only: eigen_pzheevd
      use scalapack_module, only: create_gridmap, init_blacs
#endif
      use salmon_global, only: yn_scalapack, yn_eigenexa
      implicit none
      complex(8) :: wrk
      complex(8),allocatable :: mat_H(:,:),mat_V(:,:)
      logical :: use_scalapack_diag

#ifdef USE_SCALAPACK
      use_scalapack_diag = (yn_scalapack=='y')
      if (use_scalapack_diag .and. dc%isize_tot>1 .and. dc%info_tot%isize_o==1) then
        use_scalapack_diag = .false.
        if(dc%id_tot==0) write(*,*) '[WARN] LCFO-SOI: disable ScaLAPACK when isize_o==1 in multi-rank DC run (fallback to ZHEEV).'
      end if
      if(use_scalapack_diag) then
        if(.not. allocated(dc%info_tot%gridmap)) call create_gridmap(dc%info_tot)
        if(.not. dc%info_tot%flag_blacs_gridinit .or. maxval(n_mat) > dc%info_tot%desca(3)) then
          if(dc%id_tot==0 .and. dc%info_tot%flag_blacs_gridinit .and. maxval(n_mat) > dc%info_tot%desca(3)) then
            write(*,'(1x,a,i0,a,i0,a)') '[LCFO-SOI] reinit BLACS descriptor: old_n=', dc%info_tot%desca(3), &
              ', required_n=', maxval(n_mat), '.'
          end if
          call init_blacs(dc%info_tot,maxval(n_mat))
        end if
      end if
#else
      use_scalapack_diag = .false.
#endif

      do ispin=1,nspin
        if(dc%id_tot==0) write(*,*) "Hermitian diag (SOI), #dim=",n_mat(ispin)
        n = n_mat(ispin)
        if(n<=0) then
          if(dc%id_tot==0) write(*,*) "[WARN] skip diagonalization in LCFO(SOI): n_mat<=0 for spin=",ispin
          cycle
        end if
        allocate(mat_H(n,n),mat_V(n,n))
        mat_V = (0d0,0d0)
        if(dc%id_frag==0) then
          ifrag = dc%i_frag
          do io=1,n_basis(ifrag,ispin) ; i = index_basis(io,ifrag,ispin)
          do jo=1,n_basis(ifrag,ispin) ; j = index_basis(jo,ifrag,ispin)
            mat_V(i,j) = mat_H_local(io,jo,ispin)
          end do
          end do

          do i_halo=1,n_halo ; jfrag = halo(i_halo)%ifrag_src
            do jo=1,n_basis(jfrag,ispin) ; j = index_basis(jo,jfrag,ispin)
            do io=1,n_basis(ifrag,ispin) ; i = index_basis(io,ifrag,ispin)
              wrk = 0.5d0 * halo(i_halo)%mat_H_local(jo,io,ispin)
              mat_V(j,i) = mat_V(j,i) + wrk
              mat_V(i,j) = mat_V(i,j) + conjg(wrk)
            end do
            end do
          end do
        end if
        call comm_summation(mat_V,mat_H,n*n,dc%icomm_tot)
        if(use_scalapack_diag) then
#ifdef USE_SCALAPACK
          if(dc%id_tot==0) write(*,*) "[LCFO-SOI] diagonalizer = ScaLAPACK PZHEEVD"
          if(.not. dc%info_tot%flag_blacs_gridinit) then
            if(dc%id_tot==0) write(*,*) "[WARN] BLACS grid is not initialized; fallback to ZHEEV."
            call eigen_zheev(mat_H,esp_tot(1:n,ispin),mat_V)
          else
            call eigen_pzheevd(dc%info_tot,mat_H,esp_tot(1:n,ispin),mat_V)
          end if
#else
          if(dc%id_tot==0) write(*,*) "[WARN] yn_scalapack='y' but ScaLAPACK is not enabled; fallback to ZHEEV."
          call eigen_zheev(mat_H,esp_tot(1:n,ispin),mat_V)
#endif
        else
          if(dc%id_tot==0) write(*,*) "[LCFO-SOI] diagonalizer = LAPACK ZHEEV"
          if(dc%id_tot==0 .and. yn_eigenexa=='y') then
            write(*,*) "[WARN] Complex EigenExa path is unavailable in LCFO(SOI); using ZHEEV."
          end if
          call eigen_zheev(mat_H,esp_tot(1:n,ispin),mat_V)
        end if

        if(dc%id_frag==0) then
          ifrag = dc%i_frag
          coef_wf(1:dc%nstate_frag,1:dc%nstate_tot,ispin) = (0d0,0d0)
          if (n < dc%nstate_tot .and. dc%id_tot==0) then
            write(*,'(1x,a,i0,a,i0,a)') '[WARN] n_mat < nstate_tot in diag_hermitian(SOI): n_mat=', n, &
              ' nstate_tot=', dc%nstate_tot, ' (truncating coefficients)'
          end if
          do i=1,min(n,dc%nstate_tot)
          do jo=1,n_basis(ifrag,ispin) ; j = index_basis(jo,ifrag,ispin)
            coef_wf(jo,i,ispin) = mat_V(j,i)
          end do
          end do
        end if
        deallocate(mat_H,mat_V)
      end do
    end subroutine diag_hermitian

    subroutine output
      use salmon_global, only: base_directory, sysname, unit_energy
      use filesystem, only: get_filehandle
      use inputoutput, only: uenergy_from_au
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      implicit none
      integer :: iunit,i_halo
      integer :: nxyz_domain(3), nxyz_buffer_seed(3), nxyz_box(3)
      integer :: ibx, iby, ibz, sx, sy, sz, raw_io, ibasis
      integer :: lb_zwf(3), ub_zwf(3), io_lb, io_ub
      character(256) :: filename
      real(8) :: esp_out
      complex(8) :: coef_val
      complex(8), allocatable :: phi_box_local(:,:,:), phi_box_sum(:,:,:)

      if(dc%id_tot==0 .and. yn_dc_lcfo_diag=='y') then
        iunit = get_filehandle()
        filename = trim(dc%base_directory)//trim(sysname)//"_eigen_soi.data"
        open(iunit,file=filename)
        write(iunit,'("#esp: single-particle energies (eigen energies) calculated by DC-LCFO(SOI) method")')
        write(iunit,'("#io: orbital index")')
        select case(unit_energy)
        case('au','a.u.')
          write(iunit,'("# 1:io, 2:esp[a.u.]")')
        case('ev','eV')
          write(iunit,'("# 1:io, 2:esp[eV]")')
        end select
        do ispin=1,nspin
          write(iunit,'("# spin=",1x,i5)') ispin
          if(n_mat(ispin) < dc%nstate_tot) then
            write(iunit,'("# warning: only first",1x,i0,1x,"states are available in LCFO(SOI)")') n_mat(ispin)
          end if
          do i=1,dc%nstate_tot
            if(i <= n_mat(ispin)) then
              esp_out = esp_tot(i,ispin)
            else
              ! States beyond n_mat are unavailable; use NaN as an explicit sentinel
              ! so that readers can detect missing data rather than silently using zero.
              esp_out = ieee_value(esp_out, ieee_quiet_nan)
            end if
            write(iunit,'(1x,i5,e26.16e3)') i,esp_out*uenergy_from_au
          end do
        end do
        close(iunit)
      end if

      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)

      if(dc%id_frag==0) then
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_rg_soi
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) lg%num(1:3), dc%lg_tot%num(1:3)
        do n=1,3
          write(iunit) dc%jxyz_tot(1:lg%num(n),n)
        end do
        close(iunit)

        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_bf_soi
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) nxyz_domain(1:3),nspin,dc%nstate_frag
        write(iunit) n_basis(dc%i_frag,1:nspin)
        write(iunit) f_basis(1:nxyz_domain(1),1:nxyz_domain(2),1:nxyz_domain(3),1:nspin,1:dc%nstate_frag)
        close(iunit)
      end if

      ! DG surface flux needs the same SOI basis on the DC buffer box, not
      ! only the core element.  Apply the core orthonormalization transform
      ! to the buffered spinor wavefunctions and write an unwrapped box.
      nxyz_buffer_seed(1:3) = dc%nxyz_buffer(1:3)
      nxyz_box(1:3) = nxyz_domain(1:3) + 2 * nxyz_buffer_seed(1:3)
      lb_zwf(1) = lbound(spsi%zwf, 1)
      lb_zwf(2) = lbound(spsi%zwf, 2)
      lb_zwf(3) = lbound(spsi%zwf, 3)
      ub_zwf(1) = ubound(spsi%zwf, 1)
      ub_zwf(2) = ubound(spsi%zwf, 2)
      ub_zwf(3) = ubound(spsi%zwf, 3)
      io_lb = lbound(spsi%zwf, 5)
      io_ub = ubound(spsi%zwf, 5)
      allocate(phi_box_local(nxyz_box(1), nxyz_box(2), nxyz_box(3)))
      allocate(phi_box_sum(nxyz_box(1), nxyz_box(2), nxyz_box(3)))
      if(dc%id_frag==0) then
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_bfb_soi
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) basis_buffer_magic_soi, basis_buffer_version_soi
        write(iunit) nxyz_domain(1:3), nxyz_buffer_seed(1:3), nspin, dc%nstate_frag
        write(iunit) n_basis(dc%i_frag,1:nspin)
      end if
      do ispin=1,nspin
        do ibasis=1,dc%nstate_frag
          phi_box_local = (0d0,0d0)
          if(ibasis <= n_basis(dc%i_frag,ispin)) then
            do raw_io=max(1,io_lb),min(dc%nstate_frag,io_ub)
              coef_val = basis_transform(raw_io,ibasis,ispin)
              if(abs(coef_val) <= 0d0) cycle
              do ibz=1,nxyz_box(3)
                sz = dc_buffer_box_to_local_index_soi(ibz,nxyz_domain(3),nxyz_buffer_seed(3))
                if(sz < lb_zwf(3) .or. sz > ub_zwf(3)) cycle
                do iby=1,nxyz_box(2)
                  sy = dc_buffer_box_to_local_index_soi(iby,nxyz_domain(2),nxyz_buffer_seed(2))
                  if(sy < lb_zwf(2) .or. sy > ub_zwf(2)) cycle
                  do ibx=1,nxyz_box(1)
                    sx = dc_buffer_box_to_local_index_soi(ibx,nxyz_domain(1),nxyz_buffer_seed(1))
                    if(sx < lb_zwf(1) .or. sx > ub_zwf(1)) cycle
                    phi_box_local(ibx,iby,ibz) = phi_box_local(ibx,iby,ibz) &
                    & + coef_val * spsi%zwf(sx,sy,sz,ispin,raw_io,1,1)
                  end do
                end do
              end do
            end do
          end if
          call comm_summation(phi_box_local,phi_box_sum,product(nxyz_box),dc%icomm_frag)
          if(dc%id_frag==0) write(iunit) phi_box_sum(1:nxyz_box(1),1:nxyz_box(2),1:nxyz_box(3))
        end do
      end do
      if(dc%id_frag==0) close(iunit)
      deallocate(phi_box_local,phi_box_sum)

      if(dc%id_frag==0) then
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_hl_soi
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) mat_H_local(1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        write(iunit) n_halo
        do i_halo=1,n_halo
          write(iunit) halo(i_halo)%mat_H_local(1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        end do
        close(iunit)

        if(yn_dc_lcfo_diag=='y') then
          iunit = get_filehandle()
          filename = trim(base_directory)//binfile_wf_soi
          open(iunit,file=filename,form='unformatted',access='stream')
          write(iunit) dc%n_frag, nspin, dc%nstate_frag, dc%nstate_tot
          write(iunit) n_mat(1:nspin)
          write(iunit) n_basis(1:dc%n_frag,1:nspin)
          write(iunit) index_basis(1:dc%nstate_frag,1:dc%n_frag,1:nspin)
          write(iunit) coef_wf(1:dc%nstate_frag,1:dc%nstate_tot,1:nspin)
          close(iunit)
        end if
      end if
    end subroutine output

    integer function dc_buffer_box_to_local_index_soi(ibox, ndom, nbuf) result(iloc)
      implicit none
      integer, intent(in) :: ibox, ndom, nbuf

      if(nbuf <= 0) then
        iloc = ibox
      else if(ibox <= nbuf) then
        iloc = ndom + nbuf + ibox
      else
        iloc = ibox - nbuf
      end if
    end function dc_buffer_box_to_local_index_soi

  end subroutine dc_lcfo_soi

end module lcfo_soi
