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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

! DC-LCFO method [Phys. Rev. B 95, 045106 (2017).]

#include "config.h"
module lcfo_flux
  use dc_fragment_geometry, only: get_fragment_domain
  implicit none

  private
  public :: dc_lcfo_flux

  character(32),parameter :: binfile_wf = "wavefunctions.bin", &
  &                          binfile_rg = "rgrid_index.bin", &
  &                          binfile_bf = "basis_functions.bin", &
  &                          binfile_bfb = "basis_functions_buffer.bin", &
  &                          binfile_hl = "hamiltonian_local.bin", &
  &                          binfile_vl = "velocity_local.bin"
  integer, parameter :: basis_buffer_magic = -22022212
  integer, parameter :: basis_buffer_version = 2

contains

  subroutine dc_lcfo_flux(lg,mg,system,info,stencil,ppg,energy,v_local,spsi,shpsi,sttpsi,srg,dc)
    use communication, only: comm_summation
    use salmon_global, only: yn_dc_lcfo_diag
    use structures
    implicit none
    type(s_rgrid),        intent(in) :: lg,mg
    type(s_dft_system),   intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_stencil),      intent(in) :: stencil
    type(s_pp_grid),      intent(in) :: ppg
    type(s_dft_energy),   intent(in) :: energy
    type(s_scalar),       intent(in) :: V_local(system%nspin)
    type(s_orbital),      intent(in) :: spsi
    type(s_orbital)                  :: shpsi,sttpsi
    type(s_sendrecv_grid)            :: srg
    type(s_dcdft)                    :: dc
    !
    type halo_info
      integer :: id_src,id_dst,ifrag_src,dvec(3),length(3),dsp_send(3),dsp_recv(3),axis
      real(8),allocatable :: mat_H_local(:,:,:)
      real(8),allocatable :: mat_V_local(:,:,:,:)
      real(8),allocatable :: trace_send(:,:,:),trace_recv(:,:,:)
    end type halo_info
    !
    type(halo_info) :: halo(26) ! 26 = 3^3-1
    integer :: nspin,n_halo
    integer :: stencil_radius(3)
    integer :: id_array(dc%n_frag)
    integer :: n_basis(dc%n_frag,system%nspin), n_mat(system%nspin)
    integer :: index_basis(dc%nstate_frag,dc%n_frag,system%nspin)
    real(8) :: hvol
    real(8),allocatable :: f_basis(:,:,:,:,:),hf(:,:,:,:,:),wrk_array(:,:,:,:,:) &
    & ,esp_tot(:,:),mat_H_local(:,:,:),mat_V_local(:,:,:,:),coef_wf(:,:,:),basis_transform(:,:,:)
    !
    integer :: i,j,n,ix,iy,iz,io,jo,ispin,ifrag,jfrag,i_halo

    if(dc%id_tot==0) write(*,*) "start DC-LCFO-Flux"
    hvol = system%hvol
    nspin = system%nspin
    call init_lcfo
    call calc_basis
    call hpsi_basis
    if(dc%id_tot==0) write(*,*) "basis functions operation: done"

    call calc_hamiltonian_matrix
    if(dc%id_tot==0) write(*,*) "Hamiltonian matrix: done"

    if(yn_dc_lcfo_diag=='y') then
      allocate(esp_tot(maxval(n_mat),nspin))
      if(dc%id_frag==0) allocate(coef_wf(dc%nstate_frag,dc%nstate_tot,nspin))
      if(allocated(coef_wf)) coef_wf = 0d0
#ifdef USE_EIGENEXA
      call diag_eigenexa
#else
      stop "DC-LCFO-Flux requires EigenExa; dense LAPACK fallback is disabled."
#endif
      if(dc%id_tot==0) write(*,*) "diagonalization: done"
    end if

    call output

    if(allocated(coef_wf)) deallocate(coef_wf)
    if(allocated(f_basis)) deallocate(f_basis)
    if(allocated(basis_transform)) deallocate(basis_transform)
    if(allocated(esp_tot)) deallocate(esp_tot)
    if(allocated(mat_H_local)) deallocate(mat_H_local)
    if(allocated(mat_V_local)) deallocate(mat_V_local)
    do i=1,n_halo
      if(allocated(halo(i)%mat_H_local)) deallocate(halo(i)%mat_H_local)
      if(allocated(halo(i)%mat_V_local)) deallocate(halo(i)%mat_V_local)
    end do
    if(dc%id_tot==0) write(*,*) "end DC-LCFO-Flux"

  contains

    subroutine init_lcfo
      use salmon_global, only: num_fragment
      implicit none
      integer :: lx,ly,lz,nonzero_dirs
      integer :: nxyz_domain(3)
      integer,dimension(3) :: nh,ir1,ir2,d
      integer :: id_tmp(dc%n_frag)

      id_tmp = 0
      if(.not. stencil%if_orthogonal) stop "DC-LCFO-Flux: nonorthogonal lattice is not supported"
      if(dc%optimized_fragment_geometry) stop "DC-LCFO-Flux: optimized fragment geometry is not supported"
      if(mod(dc%isize_tot,dc%n_frag) /= 0) stop "DC-LCFO-Flux: MPI size must be divisible by number of fragments"
      if(dc%id_frag==0) id_tmp(dc%i_frag) = dc%id_tot + 1
      call comm_summation(id_tmp,id_array,dc%n_frag,dc%icomm_tot)
      id_array = id_array - 1
      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)
      do n=1,3
        stencil_radius(n) = active_laplacian_radius(n)
      end do

      nh = 0
      do n=1,3 ! x,y,z
        if(dc%nxyz_buffer(n) > nxyz_domain(n)) stop "DC-LCFO: buffer > domain"
        if(num_fragment(n) > 1 .and. dc%nxyz_buffer(n) < stencil_radius(n)) &
        & stop "DC-LCFO-Flux: buffer is smaller than the active stencil radius"
        if(num_fragment(n) > 1) nh(n) = 1
      end do

      i = 0
      do lx=-nh(1),nh(1)
      do ly=-nh(2),nh(2)
      do lz=-nh(3),nh(3)
        if(lx==0 .and. ly==0 .and. lz==0) cycle
        nonzero_dirs = count([lx,ly,lz] /= 0)
        if(nonzero_dirs /= 1) cycle
        i = i + 1
        halo(i)%dvec(1:3) = [lx, ly, lz]
        halo(i)%axis = 0
        do n=1,3
          if(halo(i)%dvec(n) /= 0) halo(i)%axis = n
        end do
        halo(i)%id_dst = -1
        halo(i)%id_src = -1
        do ifrag=1,dc%n_frag
        ! dc%ixyz_frag: r-grid index of the fragment origin
          ir1(1:3) = dc%ixyz_frag(1:3,ifrag) ! position of fragment ifrag
        ! dst neighbor (+)
          ir2(1:3) = dc%ixyz_frag(1:3,dc%i_frag) + halo(i)%dvec(1:3)*nxyz_domain(1:3) ! neighbor fragment
          d(1:3) = mod( ir1(1:3) - ir2(1:3) , dc%lg_tot%num(1:3) )
          if(d(1)==0 .and. d(2)==0 .and. d(3)==0 .and. halo(i)%id_dst < 0) then
            halo(i)%id_dst = id_array(ifrag) ! process ID of the communication destination
          end if
        ! src neighbor (-)
          ir2(1:3) = dc%ixyz_frag(1:3,dc%i_frag) - halo(i)%dvec(1:3)*nxyz_domain(1:3) ! neighbor fragment
          d(1:3) = mod( ir1(1:3) - ir2(1:3) , dc%lg_tot%num(1:3) )
          if(d(1)==0 .and. d(2)==0 .and. d(3)==0 .and. halo(i)%id_src < 0) then
            halo(i)%id_src = id_array(ifrag) ! process ID of the communication source
            halo(i)%ifrag_src = ifrag
          end if
        end do ! ifrag
        if(halo(i)%id_dst < 0 .or. halo(i)%id_src < 0) stop "DC-LCFO: dst, src"
      end do
      end do
      end do
      n_halo = i ! # of the halo regions (neighbor fragments)

      do i=1,n_halo
        do n=1,3 ! x,y,z
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
      use eigen_subdiag_sub, only: eigen_dsyev
      use salmon_global, only: energy_cut,lambda_cut
      implicit none
      integer, parameter :: n_spectrum_thr = 7
      real(8), parameter :: spectrum_thr(n_spectrum_thr) = &
        [1.0d-3, 1.0d-6, 1.0d-8, 1.0d-10, 1.0d-12, 1.0d-14, 1.0d-16]
      integer :: nb(nspin),itmp(dc%n_frag,nspin)
      integer :: nxyz_domain(3),ix1,ix2,iy1,iy2,iz1,iz2
      integer :: ithr, count_by_frag(n_spectrum_thr,dc%n_frag,nspin)
      integer :: count_all(n_spectrum_thr,dc%n_frag,nspin)
      real(8),dimension(dc%nstate_frag,dc%nstate_frag,system%nspin) :: mat_S,mat_U
      real(8),dimension(dc%nstate_frag,system%nspin) :: lambda
      real(8) :: alpha_gs, norm_basis
      logical :: active_state

      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)
      ix1 = max(mg%is(1),1); ix2 = min(mg%ie(1),nxyz_domain(1))
      iy1 = max(mg%is(2),1); iy2 = min(mg%ie(2),nxyz_domain(2))
      iz1 = max(mg%is(3),1); iz2 = min(mg%ie(3),nxyz_domain(3))
      if(dc%id_tot==0 .and. energy_cut <= 0d0) then
        write(*,'(1x,a)') &
          "[DC-LCFO-FLUX-BASIS] warning: energy_cut<=0 keeps only local states below mu; RT f-sum may be incomplete"
      end if

      allocate(f_basis  (nxyz_domain(1),nxyz_domain(2),nxyz_domain(3),nspin,dc%nstate_frag))
      allocate(wrk_array(nxyz_domain(1),nxyz_domain(2),nxyz_domain(3),nspin,dc%nstate_frag))
      if(.not. allocated(basis_transform)) allocate(basis_transform(dc%nstate_frag,dc%nstate_frag,nspin))
      basis_transform = 0d0

    ! f_basis <-- | \bar{\phi} > (projected fragment orbitals)
      wrk_array = 0d0
!$omp parallel do collapse(2) private(io,ispin,iz,iy,ix,active_state) schedule(static)
      do io=info%io_s,info%io_e
      do ispin=1,nspin
        active_state = energy%esp(io,1,ispin) - system%mu < energy_cut
        if(active_state) then
      do iz=iz1,iz2
      do iy=iy1,iy2
      do ix=ix1,ix2
          wrk_array(ix,iy,iz,ispin,io) = spsi%rwf(ix,iy,iz,ispin,io,1,1) ! | \phi > @ core domain
      end do
      end do
      end do
        end if
      end do
      end do
!$omp end parallel do
      call comm_summation(wrk_array,f_basis,product(nxyz_domain)*nspin*dc%nstate_frag,info%icomm_rko)

    ! mat_S <-- S_{ij} = < \bar{\phi}_i | \bar{\phi}_j > (overlap matrix)
!$omp parallel do collapse(3) private(ispin,io,jo) schedule(static)
      do ispin=1,nspin
      do io=1,dc%nstate_frag
      do jo=1,dc%nstate_frag
        mat_S(io,jo,ispin) = sum(f_basis(:,:,:,ispin,io)*f_basis(:,:,:,ispin,jo)) * hvol
      end do
      end do
      end do
!$omp end parallel do

    ! diagonalize mat_S
      do ispin=1,nspin
        call eigen_dsyev(mat_S(:,:,ispin),lambda(:,ispin),mat_U(:,:,ispin))
      end do
      count_by_frag = 0
      if(dc%id_frag==0) then
        do ispin=1,nspin
          do ithr=1,n_spectrum_thr
            count_by_frag(ithr,dc%i_frag,ispin) = count(lambda(:,ispin) > spectrum_thr(ithr))
          end do
        end do
      end if
      call comm_summation(count_by_frag,count_all,n_spectrum_thr*dc%n_frag*nspin,dc%icomm_tot)
      if(dc%id_tot==0) then
        do ispin=1,nspin
          write(*,'(1x,a,i0,7(a,1pe9.1,a,i0,a,i0))') &
            "[DC-LCFO-FLUX-BASIS-SPECTRUM] ispin=", ispin, &
            " thr=", spectrum_thr(1), " min/max=", minval(count_all(1,:,ispin)), "/", maxval(count_all(1,:,ispin)), &
            " thr=", spectrum_thr(2), " min/max=", minval(count_all(2,:,ispin)), "/", maxval(count_all(2,:,ispin)), &
            " thr=", spectrum_thr(3), " min/max=", minval(count_all(3,:,ispin)), "/", maxval(count_all(3,:,ispin)), &
            " thr=", spectrum_thr(4), " min/max=", minval(count_all(4,:,ispin)), "/", maxval(count_all(4,:,ispin)), &
            " thr=", spectrum_thr(5), " min/max=", minval(count_all(5,:,ispin)), "/", maxval(count_all(5,:,ispin)), &
            " thr=", spectrum_thr(6), " min/max=", minval(count_all(6,:,ispin)), "/", maxval(count_all(6,:,ispin)), &
            " thr=", spectrum_thr(7), " min/max=", minval(count_all(7,:,ispin)), "/", maxval(count_all(7,:,ispin))
        end do
      end if

    ! f_basis <-- | lambda > (basis functions)
      wrk_array = f_basis
      f_basis = 0d0
      do ispin=1,nspin
        i = 0 ! count # of basis functions
        do io=dc%nstate_frag,1,-1
          if( lambda(io,ispin) > lambda_cut ) then ! cutoff for the eigenvalues of the overlap matrix
            i = i + 1 ! count # of basis functions
            do jo=1,dc%nstate_frag
              f_basis(:,:,:,ispin,i) = f_basis(:,:,:,ispin,i) &
              & + wrk_array(:,:,:,ispin,jo) * mat_U(jo,io,ispin) / sqrt(lambda(io,ispin))
              basis_transform(jo,i,ispin) = mat_U(jo,io,ispin) / sqrt(lambda(io,ispin))
            end do
          end if
        end do ! io
        nb(ispin) = i ! # of basis functions
      end do ! ispin

    ! Gram–Schmidt orthonormalization
      wrk_array = f_basis
      do ispin=1,nspin
        do io=1,nb(ispin)
          do jo=1,io-1
            alpha_gs = (sum(f_basis(:,:,:,ispin,jo)*wrk_array(:,:,:,ispin,io)) * hvol) &
            & / (sum(f_basis(:,:,:,ispin,jo)*f_basis(:,:,:,ispin,jo)) * hvol)
            wrk_array(:,:,:,ispin,io) = wrk_array(:,:,:,ispin,io) &
            & - f_basis(:,:,:,ispin,jo) * alpha_gs
            basis_transform(:,io,ispin) = basis_transform(:,io,ispin) &
            & - basis_transform(:,jo,ispin) * alpha_gs
          end do
          norm_basis = sqrt( sum(wrk_array(:,:,:,ispin,io)*wrk_array(:,:,:,ispin,io)) * hvol )
          if(norm_basis <= 1.0d-12) stop "DC-LCFO-Flux: null basis after core-S cleanup"
          basis_transform(:,io,ispin) = basis_transform(:,io,ispin) &
          & / norm_basis
          wrk_array(:,:,:,ispin,io) = wrk_array(:,:,:,ispin,io) &
          & / norm_basis
        end do
      end do ! ispin
      f_basis = wrk_array

    ! sttpsi <-- f_basis == | lambda > (basis functions)
      sttpsi%rwf = 0d0
!$omp parallel do collapse(2) private(io,ispin,iz,iy,ix) schedule(static)
      do io=info%io_s,info%io_e
      do ispin=1,nspin
      do iz=iz1,iz2
      do iy=iy1,iy2
      do ix=ix1,ix2
        sttpsi%rwf(ix,iy,iz,ispin,io,1,1) = f_basis(ix,iy,iz,ispin,io)
      end do
      end do
      end do
      end do
      end do
!$omp end parallel do

    ! n_basis: # of basis functions
      itmp = 0
      if(dc%id_frag==0) itmp(dc%i_frag,1:nspin) = nb(1:nspin)
      call comm_summation(itmp,n_basis,dc%n_frag*nspin,dc%icomm_tot)
      index_basis = 0
      do ispin=1,nspin
        i = 0
        do ifrag=1,dc%n_frag
          do io=1,n_basis(ifrag,ispin)
            i = i + 1
            index_basis(io,ifrag,ispin) = i ! index_basis: index for the total matrix
          end do
        end do
        n_mat(ispin) = i ! n_mat: dimension of the total matrix
      end do

      deallocate(wrk_array)

    end subroutine calc_basis

    subroutine hpsi_basis
      use hamiltonian, only: hpsi
      implicit none

      allocate(hf       (lg%num(1),lg%num(2),lg%num(3),nspin,dc%nstate_frag))
      allocate(wrk_array(lg%num(1),lg%num(2),lg%num(3),nspin,dc%nstate_frag))

    ! shpsi <-- H | lambda > (Hamiltonian operation)
      call hpsi(sttpsi,shpsi,info,mg,v_local,system,stencil,srg,ppg)

    ! hf <-- shpsi == H | lambda >
      wrk_array = 0d0
!$omp parallel do collapse(2) private(io,ispin,iz,iy,ix) schedule(static)
      do io=info%io_s,info%io_e
      do ispin=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        wrk_array(ix,iy,iz,ispin,io) = shpsi%rwf(ix,iy,iz,ispin,io,1,1)
      end do
      end do
      end do
      end do
      end do
!$omp end parallel do
      call comm_summation(wrk_array,hf,product(lg%num)*nspin*dc%nstate_frag,info%icomm_rko)

      deallocate(wrk_array)

    end subroutine hpsi_basis

    subroutine exchange_surface_trace_halo()
      use communication, only: comm_isend, comm_irecv, comm_wait_all
      implicit none
      integer :: axis, side, send_side, face_pt, npts
      integer :: ix, iy, iz, io, ispin
      integer :: itag_send, itag_recv, itag_dir
      integer, dimension(n_halo) :: ireq_send, ireq_recv

      if(dc%id_frag /= 0) return

      do i_halo=1,n_halo
        axis = halo(i_halo)%axis
        send_side = halo(i_halo)%dvec(axis)
        ! The receiver uses dnu_r with the receiver-side outward normal.
        ! Convert the local face derivative before sending so the surface
        ! formula is identical on both endpoints of the face.
        side = -send_side
        npts = face_point_count(axis)
        if(allocated(halo(i_halo)%trace_send)) deallocate(halo(i_halo)%trace_send)
        if(allocated(halo(i_halo)%trace_recv)) deallocate(halo(i_halo)%trace_recv)
        allocate(halo(i_halo)%trace_send(npts,nspin,2*dc%nstate_frag))
        allocate(halo(i_halo)%trace_recv(npts,nspin,2*dc%nstate_frag))
        halo(i_halo)%trace_send = 0d0
        halo(i_halo)%trace_recv = 0d0

        face_pt = 0
        do iz=1,dc%nxyz_domain(3)
        do iy=1,dc%nxyz_domain(2)
        do ix=1,dc%nxyz_domain(1)
          if(face_axis_index([ix,iy,iz], axis) /= face_coord(axis, send_side)) cycle
          face_pt = face_pt + 1
          do ispin=1,nspin
            do io=1,n_basis(dc%i_frag,ispin)
              halo(i_halo)%trace_send(face_pt,ispin,io) = &
              & local_basis_value(ix,iy,iz,ispin,io)
              halo(i_halo)%trace_send(face_pt,ispin,dc%nstate_frag+io) = &
              & local_basis_dn(ix,iy,iz,ispin,io,axis,side)
            end do
          end do
        end do
        end do
        end do

        itag_dir = halo_tag_offset(halo(i_halo)%dvec)
        itag_send = dc%i_frag + itag_dir*dc%n_frag
        ireq_send(i_halo) = comm_isend(halo(i_halo)%trace_send,halo(i_halo)%id_dst,itag_send,dc%icomm_tot)
        itag_recv = halo(i_halo)%ifrag_src + itag_dir*dc%n_frag
        ireq_recv(i_halo) = comm_irecv(halo(i_halo)%trace_recv,halo(i_halo)%id_src,itag_recv,dc%icomm_tot)
      end do
      call comm_wait_all(ireq_recv)
      call comm_wait_all(ireq_send)

    end subroutine exchange_surface_trace_halo

    subroutine release_surface_trace_halo()
      implicit none

      if(dc%id_frag /= 0) return
      do i_halo=1,n_halo
        if(allocated(halo(i_halo)%trace_send)) deallocate(halo(i_halo)%trace_send)
        if(allocated(halo(i_halo)%trace_recv)) deallocate(halo(i_halo)%trace_recv)
      end do

    end subroutine release_surface_trace_halo

    integer function halo_tag_offset(dvec) result(offset)
      implicit none
      integer,intent(in) :: dvec(3)

      offset = (dvec(1) + 1)*9 + (dvec(2) + 1)*3 + (dvec(3) + 1)
    end function halo_tag_offset

    subroutine calc_hamiltonian_matrix
      implicit none
      integer :: axis, side, face_pt, ix_face, iy_face, iz_face
      integer :: l(3), idir
      real(8), parameter :: surface_penalty_factor = 10.0d0
      real(8) :: area_weight, alpha, u_l, v_l, dnu_l, dnv_l, u_r, dnu_r
      real(8) :: term_sum, term_face, term_vsum, term_vface, pavg
      real(8), allocatable :: trace_local(:,:,:)
      real(8), allocatable :: grad_work(:,:,:)

    ! diagonal block < lambda_{ifrag,io} | H | lambda_{ifrag,jo} >
      allocate(mat_H_local(dc%nstate_frag,dc%nstate_frag,nspin))
      allocate(mat_V_local(3,dc%nstate_frag,dc%nstate_frag,nspin))
      mat_H_local = 0d0
      mat_V_local = 0d0
      l = dc%nxyz_domain
      do ispin=1,nspin
!$omp parallel do private(io,jo,idir) schedule(static)
      do io=1,n_basis(dc%i_frag,ispin)
      do jo=1,n_basis(dc%i_frag,ispin)
        mat_H_local(io,jo,ispin) = &
        & + sum(f_basis(1:l(1),1:l(2),1:l(3),ispin,io)*hf(1:l(1),1:l(2),1:l(3),ispin,jo)) * hvol
      end do
      end do
!$omp end parallel do
      end do
      allocate(grad_work(l(1),l(2),l(3)))
      do ispin=1,nspin
      do idir=1,3
      do jo=1,n_basis(dc%i_frag,ispin)
!$omp parallel do collapse(3) private(ix,iy,iz) schedule(static)
        do iz=1,l(3)
        do iy=1,l(2)
        do ix=1,l(1)
          grad_work(ix,iy,iz) = local_basis_grad(ix,iy,iz,ispin,jo,idir)
        end do
        end do
        end do
!$omp end parallel do
!$omp parallel do private(io) schedule(static)
        do io=1,n_basis(dc%i_frag,ispin)
          mat_V_local(idir,io,jo,ispin) = &
          & sum(f_basis(1:l(1),1:l(2),1:l(3),ispin,io)*grad_work(1:l(1),1:l(2),1:l(3))) * hvol
        end do
!$omp end parallel do
      end do
      end do
      end do
      deallocate(grad_work)

    ! DG surface jump/average/penalty terms on face-neighbor blocks.
      if(dc%id_frag==0) then
        call exchange_surface_trace_halo()
        do i_halo=1,n_halo
          axis = halo(i_halo)%axis
          side = -halo(i_halo)%dvec(axis)
          jfrag = halo(i_halo)%ifrag_src
          allocate(halo(i_halo)%mat_H_local(dc%nstate_frag,dc%nstate_frag,nspin))
          allocate(halo(i_halo)%mat_V_local(3,dc%nstate_frag,dc%nstate_frag,nspin))
          halo(i_halo)%mat_H_local = 0d0
          halo(i_halo)%mat_V_local = 0d0
          area_weight = system%hvol / system%hgs(axis)
          alpha = surface_penalty_factor / system%hgs(axis)
          allocate(trace_local(face_point_count(axis),nspin,2*dc%nstate_frag))
          trace_local = 0d0
          face_pt = 0
          do iz_face=1,dc%nxyz_domain(3)
          do iy_face=1,dc%nxyz_domain(2)
          do ix_face=1,dc%nxyz_domain(1)
            if(face_axis_index([ix_face,iy_face,iz_face], axis) /= face_coord(axis, side)) cycle
            face_pt = face_pt + 1
            do ispin=1,nspin
              do io=1,n_basis(dc%i_frag,ispin)
                trace_local(face_pt,ispin,io) = local_basis_value(ix_face,iy_face,iz_face,ispin,io)
                trace_local(face_pt,ispin,dc%nstate_frag+io) = &
                & local_basis_dn(ix_face,iy_face,iz_face,ispin,io,axis,side)
              end do
            end do
          end do
          end do
          end do
          do ispin=1,nspin
!$omp parallel do private(io,jo,face_pt,ix_face,iy_face,iz_face,u_l,v_l,dnu_l,dnv_l, &
!$omp& term_sum,term_face,term_vsum,term_vface) schedule(static)
          do io=1,n_basis(dc%i_frag,ispin)
          do jo=1,n_basis(dc%i_frag,ispin)
            term_sum = 0d0
            term_vsum = 0d0
            face_pt = 0
            do iz_face=1,dc%nxyz_domain(3)
            do iy_face=1,dc%nxyz_domain(2)
            do ix_face=1,dc%nxyz_domain(1)
              if(face_axis_index([ix_face,iy_face,iz_face], axis) /= face_coord(axis, side)) cycle
              face_pt = face_pt + 1
              v_l = trace_local(face_pt,ispin,io)
              dnv_l = trace_local(face_pt,ispin,dc%nstate_frag+io)
              u_l = trace_local(face_pt,ispin,jo)
              dnu_l = trace_local(face_pt,ispin,dc%nstate_frag+jo)
              term_face = (-0.25d0 * v_l * dnu_l - 0.25d0 * dnv_l * u_l + alpha * v_l * u_l) * area_weight
              term_vface = -0.5d0 * real(side,8) * v_l * u_l * area_weight
              term_sum = term_sum + term_face
              term_vsum = term_vsum + term_vface
            end do
            end do
            end do
            mat_H_local(io,jo,ispin) = mat_H_local(io,jo,ispin) + term_sum
            mat_V_local(axis,io,jo,ispin) = mat_V_local(axis,io,jo,ispin) + term_vsum
          end do
          end do
!$omp end parallel do
          end do
          ! halo%mat_H_local(io_remote, jo_local) stores the transpose of
          ! the local-to-remote surface block, i.e. H(remote,local).  The
          ! term below is evaluated as H(local,remote) with receiver-normal
          ! traces and is placed transposed by construction for EigenExa.
          do ispin=1,nspin
!$omp parallel do private(io,jo,face_pt,ix_face,iy_face,iz_face,v_l,dnv_l,u_r,dnu_r, &
!$omp& term_sum,term_face,term_vsum,term_vface) schedule(static)
          do io=1,n_basis(jfrag,ispin)
          do jo=1,n_basis(dc%i_frag,ispin)
            term_sum = 0d0
            term_vsum = 0d0
            face_pt = 0
            do iz_face=1,dc%nxyz_domain(3)
            do iy_face=1,dc%nxyz_domain(2)
            do ix_face=1,dc%nxyz_domain(1)
              if(face_axis_index([ix_face,iy_face,iz_face], axis) /= face_coord(axis, side)) cycle
              face_pt = face_pt + 1
              v_l = trace_local(face_pt,ispin,jo)
              dnv_l = trace_local(face_pt,ispin,dc%nstate_frag+jo)
              u_r = halo(i_halo)%trace_recv(face_pt,ispin,io)
              dnu_r = halo(i_halo)%trace_recv(face_pt,ispin,dc%nstate_frag+io)
              term_face = (-0.25d0 * v_l * dnu_r + 0.25d0 * dnv_l * u_r - alpha * v_l * u_r) * area_weight
              term_vface = 0.5d0 * real(side,8) * v_l * u_r * area_weight
              term_sum = term_sum + term_face
              term_vsum = term_vsum + term_vface
            end do
            end do
            end do
            halo(i_halo)%mat_H_local(io,jo,ispin) = halo(i_halo)%mat_H_local(io,jo,ispin) + term_sum
            halo(i_halo)%mat_V_local(axis,io,jo,ispin) = halo(i_halo)%mat_V_local(axis,io,jo,ispin) + term_vsum
          end do
          end do
!$omp end parallel do
          end do
          deallocate(trace_local)
          if(dc%id_tot==0) write(*,*) "Halo communication #",i_halo,": done"
        end do
        call release_surface_trace_halo()
      end if ! dc%id_frag==0
      do ispin=1,nspin
        do idir=1,3
          do io=1,n_basis(dc%i_frag,ispin)
            mat_V_local(idir,io,io,ispin) = 0d0
            do jo=io+1,n_basis(dc%i_frag,ispin)
              pavg = 0.5d0 * (mat_V_local(idir,io,jo,ispin) - mat_V_local(idir,jo,io,ispin))
              mat_V_local(idir,io,jo,ispin) = pavg
              mat_V_local(idir,jo,io,ispin) = -pavg
            end do
          end do
        end do
      end do
      deallocate(hf)

    end subroutine calc_hamiltonian_matrix

    integer function face_point_count(axis) result(npts)
      implicit none
      integer,intent(in) :: axis

      select case(axis)
      case(1)
        npts = dc%nxyz_domain(2) * dc%nxyz_domain(3)
      case(2)
        npts = dc%nxyz_domain(1) * dc%nxyz_domain(3)
      case(3)
        npts = dc%nxyz_domain(1) * dc%nxyz_domain(2)
      case default
        npts = 0
      end select
    end function face_point_count

    integer function face_coord(axis, side) result(idx)
      implicit none
      integer,intent(in) :: axis, side

      if(side > 0) then
        idx = dc%nxyz_domain(axis)
      else
        idx = 1
      end if
    end function face_coord

    integer function face_axis_index(idx3, axis) result(idx)
      implicit none
      integer,intent(in) :: idx3(3), axis

      idx = idx3(axis)
    end function face_axis_index

    real(8) function local_basis_value(ix,iy,iz,ispin,ibasis) result(val)
      implicit none
      integer,intent(in) :: ix,iy,iz,ispin,ibasis
      integer :: raw_io, io_lb, io_ub, sx, sy, sz
      real(8) :: coef_val

      val = 0d0
      if(ibasis < 1 .or. ibasis > dc%nstate_frag) return
      if(ispin < 1 .or. ispin > nspin) return
      if(ix >= 1 .and. ix <= dc%nxyz_domain(1) .and. &
         iy >= 1 .and. iy <= dc%nxyz_domain(2) .and. &
         iz >= 1 .and. iz <= dc%nxyz_domain(3)) then
        val = f_basis(ix,iy,iz,ispin,ibasis)
        return
      end if
      sx = buffered_local_index(ix, 1)
      sy = buffered_local_index(iy, 2)
      sz = buffered_local_index(iz, 3)
      if(sx < lbound(spsi%rwf,1) .or. sx > ubound(spsi%rwf,1)) return
      if(sy < lbound(spsi%rwf,2) .or. sy > ubound(spsi%rwf,2)) return
      if(sz < lbound(spsi%rwf,3) .or. sz > ubound(spsi%rwf,3)) return
      io_lb = lbound(spsi%rwf,5)
      io_ub = ubound(spsi%rwf,5)
      do raw_io=max(1,io_lb),min(dc%nstate_frag,io_ub)
        coef_val = basis_transform(raw_io,ibasis,ispin)
        if(abs(coef_val) <= 0d0) cycle
        val = val + coef_val * spsi%rwf(sx,sy,sz,ispin,raw_io,1,1)
      end do
    end function local_basis_value

    integer function buffered_local_index(idx, axis) result(mapped)
      implicit none
      integer,intent(in) :: idx, axis

      if(idx < 1) then
        mapped = dc%nxyz_domain(axis) + dc%nxyz_buffer(axis) + (1 - idx)
      else
        mapped = idx
      end if
    end function buffered_local_index

    real(8) function local_basis_grad(ix,iy,iz,ispin,ibasis,axis) result(grad_axis)
      implicit none
      integer,intent(in) :: ix,iy,iz,ispin,ibasis,axis
      integer :: dist

      grad_axis = 0d0
      do dist=1,size(stencil%coef_nab,1)
        select case(axis)
        case(1)
          grad_axis = grad_axis + stencil%coef_nab(dist,axis) * &
          & (local_basis_value(ix+dist,iy,iz,ispin,ibasis) - local_basis_value(ix-dist,iy,iz,ispin,ibasis))
        case(2)
          grad_axis = grad_axis + stencil%coef_nab(dist,axis) * &
          & (local_basis_value(ix,iy+dist,iz,ispin,ibasis) - local_basis_value(ix,iy-dist,iz,ispin,ibasis))
        case(3)
          grad_axis = grad_axis + stencil%coef_nab(dist,axis) * &
          & (local_basis_value(ix,iy,iz+dist,ispin,ibasis) - local_basis_value(ix,iy,iz-dist,ispin,ibasis))
        end select
      end do
    end function local_basis_grad

    real(8) function local_basis_dn(ix,iy,iz,ispin,ibasis,axis,side) result(dn)
      implicit none
      integer,intent(in) :: ix,iy,iz,ispin,ibasis,axis,side

      dn = real(side,8) * local_basis_grad(ix,iy,iz,ispin,ibasis,axis)
    end function local_basis_dn

    real(8) function volume_velocity_integral(io,jo,ispin,axis,l) result(val)
      implicit none
      integer,intent(in) :: io,jo,ispin,axis,l(3)
      integer :: ix,iy,iz

      val = 0d0
      do iz=1,l(3)
        do iy=1,l(2)
          do ix=1,l(1)
            val = val + f_basis(ix,iy,iz,ispin,io) * local_basis_grad(ix,iy,iz,ispin,jo,axis)
          end do
        end do
      end do
      val = val * hvol
    end function volume_velocity_integral

    integer function active_laplacian_radius(axis) result(radius)
      implicit none
      integer,intent(in) :: axis
      integer :: dist

      radius = 0
      do dist=1,size(stencil%coef_lap,1)
        if(abs(stencil%coef_lap(dist,axis)) > 0d0) radius = dist
      end do
    end function active_laplacian_radius

#ifdef USE_EIGENEXA
    subroutine diag_eigenexa
      use communication, only: comm_bcast, comm_summation
      use eigen_subdiag_sub, only: eigen_dsyev
      use eigen_libs_mod
      implicit none
      integer, parameter :: coef_gather_target_elems = 2000000
      integer, parameter :: velocity_diag_max_dim = 1024
      integer :: n,nx,ny,ix_s,ix_e,iy_s,iy_e,ix_loc,iy_loc,ifrag_x,ifrag_y,io_x,io_y
      integer :: nnod,x_nnod,y_nnod,inod,x_inod,y_inod
      integer :: jfrag_halo(n_halo)
      integer, allocatable :: io_array(:),ifrag_array(:)
      integer :: c0, c1, ncol, nstate_chunk
      integer :: target_frag, nbasis_diag, level, max_level, state_col, n_entry
      integer :: pass, best, tmp_i
      integer, parameter :: nsample_max = 3
      integer :: nsample, isample, sample_state(nsample_max)
      integer :: nstate_use, nocc_diag, occ, virt, idir_diag
      real(8) :: eps, hval, rel_col, rel_row, rel_col_max, rel_row_max
      real(8) :: occ_weight, gap, amp, strength_pair
      real(8) :: strength_total, strength_max, sum_gap_weighted, sum_inv_gap, occ_sum
      real(8) :: mean_gap, fsum_ratio
      real(8), allocatable :: h_div(:,:), h_ref_div(:,:), v_div(:,:), h(:,:,:)
      real(8), allocatable :: h_block(:,:,:), h_local_diag(:,:), evec_local(:,:), eval_local(:)
      real(8), allocatable :: eval_list(:)
      real(8), allocatable :: coef_state_norm(:), coef_state_norm_alt(:)
      real(8), allocatable :: v_tmp1(:,:), v_tmp2(:,:)
      integer, allocatable :: frag_list(:), level_list(:)
      real(8), allocatable :: v_col_local(:,:), v_col(:,:), hv_col_local(:,:), hv_col(:,:)
      real(8), allocatable :: v_row_local(:,:), v_row(:,:), hv_row_local(:,:), hv_row(:,:)
      real(8), allocatable :: eigvec_local(:,:), eigvec(:,:), vop(:,:), vwork(:)
      logical, allocatable :: repair_state(:)
      logical :: use_transpose_export, block_diag_h

      allocate(h(dc%nstate_frag,dc%nstate_frag,0:n_halo))
      block_diag_h = use_block_diag_hamiltonian_mode()
      if(dc%id_tot==0 .and. block_diag_h) then
        write(*,'(1x,a)') "[DC-LCFO-FLUX] block-diag H diagonalization: enabled"
      end if
      do ispin=1,nspin
        if(dc%id_tot==0) write(*,*) "eigenexa diag, #dim=",n_mat(ispin)
        n = n_mat(ispin)

        allocate(io_array(n),ifrag_array(n))
        do ifrag=1,dc%n_frag
          do io=1,n_basis(ifrag,ispin) ; i = index_basis(io,ifrag,ispin)
            io_array(i) = io
            ifrag_array(i) = ifrag
          end do
        end do

        call eigen_init(dc%icomm_tot)
        call eigen_get_matdims( n, nx, ny )
        call eigen_get_procs( nnod, x_nnod, y_nnod )
        call eigen_get_id   ( inod, x_inod, y_inod )
        allocate( h_div(nx,ny), h_ref_div(nx,ny), v_div(nx,ny) )
        ix_s = eigen_loop_start( 1, x_nnod, x_inod )
        ix_e = eigen_loop_end  ( n, x_nnod, x_inod )
        iy_s = eigen_loop_start( 1, y_nnod, y_inod )
        iy_e = eigen_loop_end  ( n, y_nnod, y_inod )
        nstate_chunk = max(1, min(dc%nstate_tot, &
          max(1, coef_gather_target_elems / max(1, dc%nstate_frag))))
        allocate(v_tmp1(dc%nstate_frag,nstate_chunk))
        allocate(v_tmp2(dc%nstate_frag,nstate_chunk))
        allocate(coef_state_norm(dc%nstate_tot), coef_state_norm_alt(dc%nstate_tot))
        allocate(repair_state(dc%nstate_tot))

        h_div = 0d0
        do ifrag=1,dc%n_frag
          if(ifrag==dc%i_frag .and. dc%id_frag==0) then
            h(:,:,0) = mat_H_local(:,:,ispin)
            do i_halo=1,n_halo
              jfrag_halo(i_halo) = halo(i_halo)%ifrag_src ! src fragment (recv)
              h(:,:,i_halo) = halo(i_halo)%mat_H_local(:,:,ispin)
            end do
          end if
          call comm_bcast( h, dc%icomm_tot, id_array(ifrag) )
          call comm_bcast( jfrag_halo, dc%icomm_tot, id_array(ifrag) )
          do iy_loc=iy_s,iy_e
            iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
            if(iy > n) cycle
            ifrag_y = ifrag_array(iy)
            io_y = io_array(iy)
            do ix_loc=ix_s,ix_e
              ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
              if(ix > n) cycle
              ifrag_x = ifrag_array(ix)
              io_x = io_array(ix)
              if(block_diag_h) then
                if(ifrag_x == ifrag_y) then
                  if(ifrag_x == ifrag) then
                    h_div(ix_loc,iy_loc) = h_div(ix_loc,iy_loc) + 0.5d0 * (h(io_x,io_y,0) + h(io_y,io_x,0))
                    do i_halo=1,n_halo
                      h_div(ix_loc,iy_loc) = h_div(ix_loc,iy_loc) &
                      & + 0.25d0 * (h(io_y,io_x,i_halo) + h(io_x,io_y,i_halo))
                    end do
                  end if
                  do i_halo=1,n_halo
                    if(ifrag_x == jfrag_halo(i_halo)) then
                      h_div(ix_loc,iy_loc) = h_div(ix_loc,iy_loc) &
                      & + 0.25d0 * (h(io_x,io_y,i_halo) + h(io_y,io_x,i_halo))
                    end if
                  end do
                end if
              else
                if(ifrag_x == ifrag .and. ifrag_y == ifrag) then
                  h_div(ix_loc,iy_loc) = h(io_x,io_y,0)
                end if
                do i_halo=1,n_halo
                  if( ifrag_x == jfrag_halo(i_halo) .and. ifrag_y == ifrag ) then
                    h_div(ix_loc,iy_loc) = h_div(ix_loc,iy_loc) &
                    & + 0.5d0 * h(io_x,io_y,i_halo) ! 0.5d0: avoid double counting face-neighbor pairs
                  else if( ifrag_x == ifrag .and. ifrag_y == jfrag_halo(i_halo) ) then
                    h_div(ix_loc,iy_loc) = h_div(ix_loc,iy_loc) &
                    & + 0.5d0 * h(io_y,io_x,i_halo) ! 0.5d0: avoid double counting face-neighbor pairs
                  end if
                end do ! i_halo
              end if
            end do ! ix_loc
          end do ! iy_loc
        end do ! ifrag
        if(dc%id_tot==0) write(*,*) "h_div: done"

        h_ref_div = h_div
        call eigen_sx(n, n, h_div, nx, esp_tot(1:n,ispin), v_div, nx)
        if(dc%id_tot==0) write(*,*) "eigen_sx: done"

        nsample = 1
        sample_state(1) = 1
        i = min(n, max(1, dc%nstate_tot / 2))
        if(i /= sample_state(1)) then
          nsample = nsample + 1
          sample_state(nsample) = i
        end if
        i = min(n, dc%nstate_tot)
        if(all(sample_state(1:nsample) /= i)) then
          nsample = nsample + 1
          sample_state(nsample) = i
        end if

        allocate(v_col_local(n,nsample), v_col(n,nsample), hv_col_local(n,nsample), hv_col(n,nsample))
        allocate(v_row_local(n,nsample), v_row(n,nsample), hv_row_local(n,nsample), hv_row(n,nsample))
        v_col_local = 0d0
        v_row_local = 0d0
        do iy_loc=iy_s,iy_e
          iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
          if(iy > n) cycle
          do ix_loc=ix_s,ix_e
            ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
            if(ix > n) cycle
            do isample=1,nsample
              if(iy == sample_state(isample)) v_col_local(ix,isample) = v_div(ix_loc,iy_loc)
              if(ix == sample_state(isample)) v_row_local(iy,isample) = v_div(ix_loc,iy_loc)
            end do
          end do
        end do
        call comm_summation(v_col_local, v_col, n*nsample, dc%icomm_tot)
        call comm_summation(v_row_local, v_row, n*nsample, dc%icomm_tot)

        hv_col_local = 0d0
        hv_row_local = 0d0
        do iy_loc=iy_s,iy_e
          iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
          if(iy > n) cycle
          do ix_loc=ix_s,ix_e
            ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
            if(ix > n) cycle
            hval = h_ref_div(ix_loc,iy_loc)
            if(hval == 0d0) cycle
            do isample=1,nsample
              hv_col_local(ix,isample) = hv_col_local(ix,isample) + hval * v_col(iy,isample)
              hv_row_local(ix,isample) = hv_row_local(ix,isample) + hval * v_row(iy,isample)
            end do
          end do
        end do
        call comm_summation(hv_col_local, hv_col, n*nsample, dc%icomm_tot)
        call comm_summation(hv_row_local, hv_row, n*nsample, dc%icomm_tot)

        rel_col_max = 0d0
        rel_row_max = 0d0
        do isample=1,nsample
          eps = esp_tot(sample_state(isample),ispin)
          rel_col = sqrt(sum((hv_col(:,isample) - eps * v_col(:,isample))**2) &
            / max(1.0d-300, sum(hv_col(:,isample)**2)))
          rel_row = sqrt(sum((hv_row(:,isample) - eps * v_row(:,isample))**2) &
            / max(1.0d-300, sum(hv_row(:,isample)**2)))
          rel_col_max = max(rel_col_max, rel_col)
          rel_row_max = max(rel_row_max, rel_row)
        end do
        use_transpose_export = (rel_row_max < rel_col_max * 1.0d-3)
        if(dc%id_tot==0) then
          write(*,'(1x,a,i0,2(a,1pe12.5),a,l1)') &
            "[DC-LCFO-FLUX] EigenExa vector check: ispin=", ispin, &
            " rel_col=", rel_col_max, " rel_row=", rel_row_max, &
            " transpose_export=", use_transpose_export
        end if
        deallocate(v_col_local, v_col, hv_col_local, hv_col)
        deallocate(v_row_local, v_row, hv_row_local, hv_row)

        nstate_use = min(dc%nstate_tot, n)
        if(n <= velocity_diag_max_dim .and. nstate_use <= velocity_diag_max_dim) then
          idir_diag = 3
          allocate(eigvec_local(n,nstate_use), eigvec(n,nstate_use))
          allocate(vop(n,n), vwork(n))
          eigvec_local = 0d0
          do iy_loc=iy_s,iy_e
            iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
            if(iy < 1 .or. iy > nstate_use) cycle
            do ix_loc=ix_s,ix_e
              ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
              if(ix > n) cycle
              eigvec_local(ix,iy) = v_div(ix_loc,iy_loc)
            end do
          end do
          call comm_summation(eigvec_local, eigvec, n*nstate_use, dc%icomm_tot)

          vop = 0d0
          do ifrag=1,dc%n_frag
            if(ifrag==dc%i_frag .and. dc%id_frag==0) then
              h(:,:,0) = mat_V_local(idir_diag,:,:,ispin)
              do i_halo=1,n_halo
                jfrag_halo(i_halo) = halo(i_halo)%ifrag_src
                h(:,:,i_halo) = halo(i_halo)%mat_V_local(idir_diag,:,:,ispin)
              end do
            end if
            call comm_bcast( h, dc%icomm_tot, id_array(ifrag) )
            call comm_bcast( jfrag_halo, dc%icomm_tot, id_array(ifrag) )
            do iy=1,n
              ifrag_y = ifrag_array(iy)
              io_y = io_array(iy)
              do ix=1,n
                ifrag_x = ifrag_array(ix)
                io_x = io_array(ix)
                if(ifrag_x == ifrag .and. ifrag_y == ifrag) then
                  if(block_diag_h) then
                    vop(ix,iy) = vop(ix,iy) + h(io_x,io_y,0)
                  else
                    vop(ix,iy) = h(io_x,io_y,0)
                  end if
                end if
                if(block_diag_h) then
                  if(ifrag_x == ifrag_y) then
                    if(ifrag_x == ifrag) then
                      do i_halo=1,n_halo
                        vop(ix,iy) = vop(ix,iy) + 0.5d0 * h(io_y,io_x,i_halo)
                      end do
                    end if
                    do i_halo=1,n_halo
                      if(ifrag_x == jfrag_halo(i_halo)) then
                        vop(ix,iy) = vop(ix,iy) + 0.5d0 * h(io_x,io_y,i_halo)
                      end if
                    end do
                  end if
                else
                  do i_halo=1,n_halo
                    if(ifrag_x == jfrag_halo(i_halo) .and. ifrag_y == ifrag) then
                      vop(ix,iy) = vop(ix,iy) + 0.5d0 * h(io_x,io_y,i_halo)
                    else if(ifrag_x == ifrag .and. ifrag_y == jfrag_halo(i_halo)) then
                      vop(ix,iy) = vop(ix,iy) + 0.5d0 * h(io_y,io_x,i_halo)
                    end if
                  end do
                end if
              end do
            end do
          end do

          nocc_diag = 0
          occ_sum = 0d0
          do occ=1,nstate_use
            occ_weight = max(0d0, system%rocc(occ,1,ispin))
            if(occ_weight > 1d-12) nocc_diag = occ
            occ_sum = occ_sum + occ_weight
          end do
          strength_total = 0d0
          strength_max = 0d0
          sum_gap_weighted = 0d0
          sum_inv_gap = 0d0
          do virt=nocc_diag+1,nstate_use
            vwork(:) = matmul(vop, eigvec(:,virt))
            do occ=1,nocc_diag
              occ_weight = max(0d0, system%rocc(occ,1,ispin))
              if(occ_weight <= 1d-12) cycle
              gap = esp_tot(virt,ispin) - esp_tot(occ,ispin)
              amp = sum(eigvec(:,occ) * vwork(:))
              strength_pair = occ_weight * amp * amp
              strength_total = strength_total + strength_pair
              strength_max = max(strength_max, strength_pair)
              sum_gap_weighted = sum_gap_weighted + strength_pair * gap
              if(gap > 1d-12) sum_inv_gap = sum_inv_gap + strength_pair / gap
            end do
          end do
          if(strength_total > 0d0) then
            mean_gap = sum_gap_weighted / strength_total
          else
            mean_gap = 0d0
          end if
          if(occ_sum > 0d0) then
            fsum_ratio = 2d0 * sum_inv_gap / occ_sum
          else
            fsum_ratio = 0d0
          end if
          if(dc%id_tot==0) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,4(a,1pe13.5))') &
              "[DC-LCFO-FLUX-TRANSITION] idir=", idir_diag, " ispin=", ispin, &
              " nocc=", nocc_diag, " nvirt=", max(0,nstate_use-nocc_diag), &
              " total=", strength_total, " max_pair=", strength_max, &
              " mean_gap_eV=", mean_gap * 27.211386245988d0, &
              " fsum_ratio=", fsum_ratio
          end if
          deallocate(eigvec_local, eigvec, vop, vwork)
        else if(dc%id_tot==0) then
          write(*,'(1x,a,i0,a,i0,a,i0)') &
            "[DC-LCFO-FLUX-TRANSITION] skipped: n=", n, &
            " nstate=", nstate_use, " max_auto=", velocity_diag_max_dim
        end if

        coef_state_norm = 0d0
        do ifrag=1,dc%n_frag
          do c0=1,dc%nstate_tot,nstate_chunk
            c1 = min(dc%nstate_tot, c0+nstate_chunk-1)
            ncol = c1-c0+1
            v_tmp1(:,1:ncol) = 0d0
            if(use_transpose_export) then
              do iy_loc=iy_s,iy_e
                iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
                if(iy > n) cycle
                ifrag_y = ifrag_array(iy)
                if(ifrag_y /= ifrag) cycle
                io_y = io_array(iy)
                do ix_loc=ix_s,ix_e
                  ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
                  if(ix < c0 .or. ix > c1 .or. ix > n) cycle
                  v_tmp1(io_y,ix-c0+1) = v_div(ix_loc,iy_loc)
                end do
              end do
            else
              do iy_loc=iy_s,iy_e
                iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
                if(iy < c0 .or. iy > c1 .or. iy > n) cycle
                do ix_loc=ix_s,ix_e
                  ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
                  if(ix > n) cycle
                  ifrag_x = ifrag_array(ix)
                  io_x = io_array(ix)
                  if(ifrag_x == ifrag) then
                    v_tmp1(io_x,iy-c0+1) = v_div(ix_loc,iy_loc)
                  end if
                end do
              end do
            end if
            call comm_summation(v_tmp1(:,1:ncol),v_tmp2(:,1:ncol), &
              dc%nstate_frag*ncol,dc%icomm_tot)
            do i=1,ncol
              coef_state_norm(c0+i-1) = coef_state_norm(c0+i-1) + sum(v_tmp2(:,i)**2)
            end do
            if(ifrag==dc%i_frag .and. dc%id_frag==0) then
              coef_wf(:,c0:c1,ispin) = v_tmp2(:,1:ncol)
            end if
          end do
        end do ! ifrag

        repair_state(:) = (coef_state_norm(:) <= 1.0d-12)
        if(any(repair_state)) then
          if(dc%id_tot==0) then
            write(*,'(1x,a,i0,a,1pe13.5)') &
              "[DC-LCFO-FLUX] repairing near-zero exported eigenvector columns: count=", &
              count(repair_state), " min_norm2=", minval(coef_state_norm)
          end if
          coef_state_norm_alt = 0d0
          do ifrag=1,dc%n_frag
            do c0=1,dc%nstate_tot,nstate_chunk
              c1 = min(dc%nstate_tot, c0+nstate_chunk-1)
              ncol = c1-c0+1
              if(.not. any(repair_state(c0:c1))) cycle
              v_tmp1(:,1:ncol) = 0d0
              if(use_transpose_export) then
                do iy_loc=iy_s,iy_e
                  iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
                  if(iy < c0 .or. iy > c1 .or. iy > n) cycle
                  do ix_loc=ix_s,ix_e
                    ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
                    if(ix > n) cycle
                    ifrag_x = ifrag_array(ix)
                    io_x = io_array(ix)
                    if(ifrag_x == ifrag) then
                      v_tmp1(io_x,iy-c0+1) = v_div(ix_loc,iy_loc)
                    end if
                  end do
                end do
              else
                do iy_loc=iy_s,iy_e
                  iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
                  if(iy > n) cycle
                  ifrag_y = ifrag_array(iy)
                  if(ifrag_y /= ifrag) cycle
                  io_y = io_array(iy)
                  do ix_loc=ix_s,ix_e
                    ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
                    if(ix < c0 .or. ix > c1 .or. ix > n) cycle
                    v_tmp1(io_y,ix-c0+1) = v_div(ix_loc,iy_loc)
                  end do
                end do
              end if
              call comm_summation(v_tmp1(:,1:ncol),v_tmp2(:,1:ncol), &
                dc%nstate_frag*ncol,dc%icomm_tot)
              do i=1,ncol
                if(.not. repair_state(c0+i-1)) cycle
                coef_state_norm_alt(c0+i-1) = coef_state_norm_alt(c0+i-1) + sum(v_tmp2(:,i)**2)
                if(ifrag==dc%i_frag .and. dc%id_frag==0) then
                  coef_wf(:,c0+i-1,ispin) = v_tmp2(:,i)
                end if
              end do
            end do
          end do
          if(dc%id_tot==0) then
            write(*,'(1x,a,1pe13.5)') &
              "[DC-LCFO-FLUX] repaired eigenvector column min_norm2=", &
              minval(coef_state_norm_alt, mask=repair_state)
          end if
        end if

        if(block_diag_h) then
          if(dc%id_tot==0) then
            write(*,'(1x,a)') &
              "[DC-LCFO-FLUX] exporting fragment-local block eigenvectors for block-diag H"
          end if
          allocate(h_block(dc%nstate_frag,dc%nstate_frag,dc%n_frag))
          allocate(eval_list(n), frag_list(n), level_list(n))
          h_block = 0d0
          eval_list = huge(1.0d0)
          frag_list = 0
          level_list = 0
          do ifrag=1,dc%n_frag
            if(ifrag==dc%i_frag .and. dc%id_frag==0) then
              h(:,:,0) = mat_H_local(:,:,ispin)
              do i_halo=1,n_halo
                jfrag_halo(i_halo) = halo(i_halo)%ifrag_src
                h(:,:,i_halo) = halo(i_halo)%mat_H_local(:,:,ispin)
              end do
            end if
            call comm_bcast( h, dc%icomm_tot, id_array(ifrag) )
            call comm_bcast( jfrag_halo, dc%icomm_tot, id_array(ifrag) )
            do target_frag=1,dc%n_frag
              if(target_frag == ifrag) then
                do io_x=1,n_basis(target_frag,ispin)
                do io_y=1,n_basis(target_frag,ispin)
                  h_block(io_x,io_y,target_frag) = h_block(io_x,io_y,target_frag) &
                  & + 0.5d0 * (h(io_x,io_y,0) + h(io_y,io_x,0))
                  do i_halo=1,n_halo
                    h_block(io_x,io_y,target_frag) = h_block(io_x,io_y,target_frag) &
                    & + 0.25d0 * (h(io_y,io_x,i_halo) + h(io_x,io_y,i_halo))
                  end do
                end do
                end do
              end if
              do i_halo=1,n_halo
                if(target_frag == jfrag_halo(i_halo)) then
                  do io_x=1,n_basis(target_frag,ispin)
                  do io_y=1,n_basis(target_frag,ispin)
                    h_block(io_x,io_y,target_frag) = h_block(io_x,io_y,target_frag) &
                    & + 0.25d0 * (h(io_x,io_y,i_halo) + h(io_y,io_x,i_halo))
                  end do
                  end do
                end if
              end do
            end do
          end do
          if(dc%id_frag==0) coef_wf(:,:,ispin) = 0d0
          n_entry = 0
          do ifrag=1,dc%n_frag
            nbasis_diag = n_basis(ifrag,ispin)
            if(nbasis_diag <= 0) cycle
            allocate(h_local_diag(nbasis_diag,nbasis_diag))
            allocate(evec_local(nbasis_diag,nbasis_diag))
            allocate(eval_local(nbasis_diag))
            h_local_diag(:,:) = h_block(1:nbasis_diag,1:nbasis_diag,ifrag)
            call eigen_dsyev(h_local_diag, eval_local, evec_local)
            max_level = min(nbasis_diag, dc%nstate_tot)
            do level=1,max_level
              if(n_entry >= n) exit
              n_entry = n_entry + 1
              eval_list(n_entry) = eval_local(level)
              frag_list(n_entry) = ifrag
              level_list(n_entry) = level
              h_block(1:nbasis_diag,level,ifrag) = evec_local(1:nbasis_diag,level)
            end do
            deallocate(h_local_diag,evec_local,eval_local)
          end do
          do pass=1,max(0,n_entry-1)
            best = pass
            do level=pass+1,n_entry
              if(eval_list(level) < eval_list(best)) best = level
            end do
            if(best /= pass) then
              eps = eval_list(pass)
              eval_list(pass) = eval_list(best)
              eval_list(best) = eps
              tmp_i = frag_list(pass)
              frag_list(pass) = frag_list(best)
              frag_list(best) = tmp_i
              tmp_i = level_list(pass)
              level_list(pass) = level_list(best)
              level_list(best) = tmp_i
            end if
          end do
          do state_col=1,min(n_entry,dc%nstate_tot)
            ifrag = frag_list(state_col)
            level = level_list(state_col)
            if(ifrag < 1 .or. ifrag > dc%n_frag) cycle
            nbasis_diag = n_basis(ifrag,ispin)
            if(level < 1 .or. level > nbasis_diag) cycle
            esp_tot(state_col,ispin) = eval_list(state_col)
            if(ifrag==dc%i_frag .and. dc%id_frag==0) then
              coef_wf(1:nbasis_diag,state_col,ispin) = h_block(1:nbasis_diag,level,ifrag)
            end if
          end do
          deallocate(h_block,eval_list,frag_list,level_list)
        end if

        deallocate(h_div,h_ref_div,v_div,v_tmp1,v_tmp2,io_array,ifrag_array)
        deallocate(coef_state_norm,coef_state_norm_alt,repair_state)
        call eigen_free()
      end do ! ispin

      deallocate(h)
    end subroutine diag_eigenexa
#endif

    logical function use_block_diag_hamiltonian_mode() result(enabled)
      implicit none
      character(16) :: env_value
      integer :: env_status

      enabled = .false.
      env_value = ''
      call get_environment_variable('SALMON_DG_BLOCK_DIAG_H', env_value, status=env_status)
      if(env_status == 0) then
        select case(trim(adjustl(env_value)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enabled = .true.
        case default
          enabled = .false.
        end select
      end if
    end function use_block_diag_hamiltonian_mode

    subroutine output
      use salmon_global, only: base_directory, sysname, unit_energy
      use filesystem, only: get_filehandle
      use inputoutput, only: uenergy_from_au
      implicit none
      integer :: iunit,i_halo
      integer :: nxyz_domain(3)
      integer :: nxyz_box(3), nxyz_buffer_seed(3)
      integer :: lb_rwf(3), ub_rwf(3), io_lb, io_ub
      integer :: ibasis, raw_io, ibx, iby, ibz, sx, sy, sz
      character(256) :: filename
      real(8) :: coef_val
      real(8), allocatable :: phi_box_local(:,:,:), phi_box_sum(:,:,:)

    ! total system data
      if(dc%id_tot==0 .and. yn_dc_lcfo_diag=='y') then
      ! eigen.data
        iunit = get_filehandle()
        filename = trim(dc%base_directory)//trim(sysname)//"_eigen.data" ! @ ./data_dcdft/total/
        open(iunit,file=filename)
        write(iunit,'("#esp: single-particle energies (eigen energies) calculated by DC-LCFO method")')
        write(iunit,'("#io: orbital index")')
        select case(unit_energy)
        case('au','a.u.')
          write(iunit,'("# 1:io, 2:esp[a.u.]")')
        case('ev','eV')
          write(iunit,'("# 1:io, 2:esp[eV]")')
        end select
        do ispin=1,nspin
          write(iunit,'("# spin=",1x,i5)') ispin
          do i=1,dc%nstate_tot
            write(iunit,'(1x,i5,e26.16e3)') i,esp_tot(i,ispin)*uenergy_from_au
          end do
        end do
        close(iunit)
      end if

    ! fragment data
      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)
      if(dc%id_frag==0) then
      ! r-grid index
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_rg ! base_directory==./data_dcdft/fragments/dc%i_frag/
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) lg%num(1:3), dc%lg_tot%num(1:3)
        do n=1,3 ! x,y,z
          write(iunit) dc%jxyz_tot(1:lg%num(n),n)
        end do
        close(iunit)
      ! basis functions | lambda >
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_bf ! base_directory==./data_dcdft/fragments/dc%i_frag/
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) nxyz_domain(1:3),nspin,dc%nstate_frag
        write(iunit) n_basis(dc%i_frag,1:nspin) ! # of basis functions
        write(iunit) f_basis(1:nxyz_domain(1),1:nxyz_domain(2),1:nxyz_domain(3) &
        & ,1:nspin,1:dc%nstate_frag) ! basis functions | lambda >
        close(iunit)
      end if

      ! buffered basis functions for DG surface Flux/Flow terms
        nxyz_buffer_seed(1:3) = dc%nxyz_buffer(1:3)
        nxyz_box(1:3) = nxyz_domain(1:3) + 2 * nxyz_buffer_seed(1:3)
        lb_rwf(1) = lbound(spsi%rwf, 1)
        lb_rwf(2) = lbound(spsi%rwf, 2)
        lb_rwf(3) = lbound(spsi%rwf, 3)
        ub_rwf(1) = ubound(spsi%rwf, 1)
        ub_rwf(2) = ubound(spsi%rwf, 2)
        ub_rwf(3) = ubound(spsi%rwf, 3)
        io_lb = lbound(spsi%rwf, 5)
        io_ub = ubound(spsi%rwf, 5)
        allocate(phi_box_local(nxyz_box(1), nxyz_box(2), nxyz_box(3)))
        allocate(phi_box_sum(nxyz_box(1), nxyz_box(2), nxyz_box(3)))
        if(dc%id_frag==0) then
          iunit = get_filehandle()
          filename = trim(base_directory)//binfile_bfb
          open(iunit,file=filename,form='unformatted',access='stream')
          write(iunit) basis_buffer_magic, basis_buffer_version
          write(iunit) nxyz_domain(1:3), nxyz_buffer_seed(1:3), nspin, dc%nstate_frag
          write(iunit) n_basis(dc%i_frag,1:nspin)
        end if
        do ispin=1,nspin
          do ibasis=1,dc%nstate_frag
            phi_box_local(:,:,:) = 0d0
            if(ibasis <= n_basis(dc%i_frag,ispin)) then
              do raw_io=max(1,io_lb),min(dc%nstate_frag,io_ub)
                coef_val = basis_transform(raw_io,ibasis,ispin)
                if(abs(coef_val) <= 0d0) cycle
                do ibz=1,nxyz_box(3)
                  sz = dc_buffer_box_to_local_index(ibz,nxyz_domain(3),nxyz_buffer_seed(3))
                  if(sz < lb_rwf(3) .or. sz > ub_rwf(3)) cycle
                  do iby=1,nxyz_box(2)
                    sy = dc_buffer_box_to_local_index(iby,nxyz_domain(2),nxyz_buffer_seed(2))
                    if(sy < lb_rwf(2) .or. sy > ub_rwf(2)) cycle
                    do ibx=1,nxyz_box(1)
                      sx = dc_buffer_box_to_local_index(ibx,nxyz_domain(1),nxyz_buffer_seed(1))
                      if(sx < lb_rwf(1) .or. sx > ub_rwf(1)) cycle
                      phi_box_local(ibx,iby,ibz) = phi_box_local(ibx,iby,ibz) &
                      & + coef_val * spsi%rwf(sx,sy,sz,ispin,raw_io,1,1)
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
      ! local hamiltonian matrix
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_hl
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) mat_H_local(1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        write(iunit) n_halo
        do i_halo=1,n_halo
          write(iunit) halo(i_halo)%mat_H_local(1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        end do
        close(iunit)
      ! local velocity matrix dH/dA at A=0, built from the same Flux-LCFO basis
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_vl
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) mat_V_local(1:3,1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        write(iunit) n_halo
        do i_halo=1,n_halo
          write(iunit) halo(i_halo)%mat_V_local(1:3,1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        end do
        close(iunit)
        if(yn_dc_lcfo_diag=='y') then
        ! coefficients of the wavefunctions
          iunit = get_filehandle()
          filename = trim(base_directory)//binfile_wf
          open(iunit,file=filename,form='unformatted',access='stream')
          write(iunit) dc%n_frag, nspin, dc%nstate_frag, dc%nstate_tot
          write(iunit) n_mat(1:nspin)
          write(iunit) n_basis(1:dc%n_frag,1:nspin)
          write(iunit) index_basis(1:dc%nstate_frag,1:dc%n_frag,1:nspin)
          write(iunit) coef_wf(1:dc%nstate_frag,1:dc%nstate_tot,1:nspin)
          write(iunit) esp_tot(1:dc%nstate_tot,1:nspin)
          close(iunit)
        end if
      end if

    end subroutine output

  end subroutine dc_lcfo_flux

  integer function dc_buffer_box_to_local_index(ibox, ndom, nbuf) result(iloc)
    implicit none
    integer, intent(in) :: ibox, ndom, nbuf

    if (nbuf <= 0) then
      iloc = ibox
    else if (ibox <= nbuf) then
      iloc = ndom + nbuf + ibox
    else
      iloc = ibox - nbuf
    end if
  end function dc_buffer_box_to_local_index

end module lcfo_flux
