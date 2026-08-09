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
module lcfo
  implicit none
  
  private
  public :: dc_lcfo, init_conventional_from_dcdft
  
  character(32),parameter :: binfile_wf = "wavefunctions.bin", &
  &                          binfile_rg = "rgrid_index.bin", &
  &                          binfile_bf = "basis_functions.bin", &
  &                          binfile_hl = "hamiltonian_local.bin"
  
contains

  subroutine dc_lcfo(lg,mg,system,info,stencil,ppg,energy,v_local,spsi,shpsi,sttpsi,srg,dc)
    use communication, only: comm_summation
    use salmon_global, only: yn_dc_lcfo_diag, lcfo_eigensolver
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
      integer :: id_src,id_dst,ifrag_src,dvec(3),length(3),dsp_send(3),dsp_recv(3)
      real(8),allocatable :: buf_send(:,:,:,:,:),buf_recv(:,:,:,:,:),mat_H_local(:,:,:)
    end type halo_info
    !
    type(halo_info) :: halo(26) ! 26 = 3^3-1
    integer :: nspin,n_halo
    integer :: id_array(dc%n_frag)
    integer :: n_basis(dc%n_frag,system%nspin), n_mat(system%nspin)
    integer :: index_basis(dc%nstate_frag,dc%n_frag,system%nspin)
    real(8) :: hvol
    real(8),allocatable :: f_basis(:,:,:,:,:),hf(:,:,:,:,:),wrk_array(:,:,:,:,:) &
    & ,esp_tot(:,:),mat_H_local(:,:,:),coef_wf(:,:,:)
    !
    integer :: i,j,n,ix,iy,iz,io,jo,ispin,ifrag,jfrag,i_halo
    
    if(dc%id_tot==0) write(*,*) "start DC-LCFO"
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
      if(dc%id_frag==0) then
        allocate(coef_wf(dc%nstate_frag,dc%nstate_tot,nspin))
        coef_wf = 0d0
      end if
      if(any(n_mat < dc%nstate_tot)) then
        stop "DC-LCFO: nstate_tot exceeds the Hamiltonian matrix dimension."
      end if
      select case(trim(lcfo_eigensolver))
      case('lapack')
        call diag_lapack
      case('eigenexa')
#ifdef USE_EIGENEXA
        call diag_eigenexa
#else
        stop "DC-LCFO: EigenExa support is not enabled."
#endif
      case('slepc')
#ifdef USE_SLEPC
        call diag_slepc
#else
        stop "DC-LCFO: SLEPc support is not enabled."
#endif
      case('chefsi')
#ifdef USE_SCALAPACK
        call diag_chefsi_driver
#else
        stop "DC-LCFO: ScaLAPACK support is not enabled."
#endif
      case default
        stop "DC-LCFO: invalid lcfo_eigensolver."
      end select
      if(dc%id_tot==0) write(*,*) "diagonalization: done"
!      call test_write_psi
    end if
  
    call output

    if(allocated(coef_wf)) deallocate(coef_wf)
    if(allocated(f_basis)) deallocate(f_basis)
    if(allocated(esp_tot)) deallocate(esp_tot)
    if(allocated(mat_H_local)) deallocate(mat_H_local)
    do i=1,n_halo
      if(allocated(halo(i)%mat_H_local)) deallocate(halo(i)%mat_H_local)
    end do
    if(dc%id_tot==0) write(*,*) "end DC-LCFO"
    
  contains
  
    subroutine init_lcfo
      use salmon_global, only: num_fragment
      implicit none
      integer :: lx,ly,lz
      integer,dimension(3) :: nh,ir1,ir2,d
      integer :: id_tmp(dc%n_frag)
      
      id_tmp = 0
      if(dc%id_frag==0) id_tmp(dc%i_frag) = dc%id_tot + 1
      call comm_summation(id_tmp,id_array,dc%n_frag,dc%icomm_tot)
      id_array = id_array - 1
      
      nh = 0
      do n=1,3 ! x,y,z
        if(dc%nxyz_buffer(n) > dc%nxyz_domain(n)) stop "DC-LCFO: buffer > domain"
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
        ! dc%ixyz_frag: r-grid index of the fragment origin
          ir1(1:3) = dc%ixyz_frag(1:3,ifrag) ! position of fragment ifrag
        ! dst neighbor (+)
          ir2(1:3) = dc%ixyz_frag(1:3,dc%i_frag) + halo(i)%dvec(1:3)*dc%nxyz_domain(1:3) ! neighbor fragment
          d(1:3) = mod( ir1(1:3) - ir2(1:3) , dc%lg_tot%num(1:3) )
          if(d(1)==0 .and. d(2)==0 .and. d(3)==0 .and. halo(i)%id_dst < 0) then
            halo(i)%id_dst = id_array(ifrag) ! process ID of the communication destination
          end if
        ! src neighbor (-)
          ir2(1:3) = dc%ixyz_frag(1:3,dc%i_frag) - halo(i)%dvec(1:3)*dc%nxyz_domain(1:3) ! neighbor fragment
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
            halo(i)%length(n) = dc%nxyz_domain(n)
            halo(i)%dsp_send(n) = 0
            halo(i)%dsp_recv(n) = 0
          case(1)
            halo(i)%length(n) = dc%nxyz_buffer(n)
            halo(i)%dsp_send(n) = dc%nxyz_domain(n) - dc%nxyz_buffer(n)
            halo(i)%dsp_recv(n) = dc%nxyz_domain(n) + dc%nxyz_buffer(n)
          case(-1)
            halo(i)%length(n) = dc%nxyz_buffer(n)
            halo(i)%dsp_send(n) = 0
            halo(i)%dsp_recv(n) = dc%nxyz_domain(n)
          end select
        end do
      end do
    
    end subroutine init_lcfo
    
    subroutine calc_basis
      use eigen_subdiag_sub, only: eigen_dsyev
      use salmon_global, only: energy_cut,lambda_cut
      implicit none
      integer :: nb(nspin),itmp(dc%n_frag,nspin)
      real(8),dimension(dc%nstate_frag,dc%nstate_frag,system%nspin) :: mat_S,mat_U
      real(8),dimension(dc%nstate_frag,system%nspin) :: lambda
      
      allocate(f_basis  (dc%nxyz_domain(1),dc%nxyz_domain(2),dc%nxyz_domain(3),nspin,dc%nstate_frag))
      allocate(wrk_array(dc%nxyz_domain(1),dc%nxyz_domain(2),dc%nxyz_domain(3),nspin,dc%nstate_frag))
      
    ! f_basis <-- | \bar{\phi} > (projected fragment orbitals)
      wrk_array = 0d0
      do io=info%io_s,info%io_e
      do ispin=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        if( ix <= dc%nxyz_domain(1) .and. iy <= dc%nxyz_domain(2) .and. iz <= dc%nxyz_domain(3) &
        & .and. energy%esp(io,1,ispin) - system%mu < energy_cut ) then ! energy cutoff
          wrk_array(ix,iy,iz,ispin,io) = spsi%rwf(ix,iy,iz,ispin,io,1,1) ! | \phi > @ core domain
        end if
      end do
      end do
      end do
      end do
      end do
      call comm_summation(wrk_array,f_basis,product(dc%nxyz_domain)*nspin*dc%nstate_frag,info%icomm_rko)
      
    ! mat_S <-- S_{ij} = < \bar{\phi}_i | \bar{\phi}_j > (overlap matrix)
      do ispin=1,nspin
      do io=1,dc%nstate_frag
      do jo=1,dc%nstate_frag
        mat_S(io,jo,ispin) = sum(f_basis(:,:,:,ispin,io)*f_basis(:,:,:,ispin,jo)) * hvol
      end do
      end do
      end do
      
    ! diagonalize mat_S
      do ispin=1,nspin
        call eigen_dsyev(mat_S(:,:,ispin),lambda(:,ispin),mat_U(:,:,ispin))
      end do
      
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
            wrk_array(:,:,:,ispin,io) = wrk_array(:,:,:,ispin,io) &
            & - f_basis(:,:,:,ispin,jo) * sum(f_basis(:,:,:,ispin,jo)*wrk_array(:,:,:,ispin,io)) &
            & / sum(f_basis(:,:,:,ispin,jo)*f_basis(:,:,:,ispin,jo))
          end do
          wrk_array(:,:,:,ispin,io) = wrk_array(:,:,:,ispin,io) &
          & / sqrt( sum(wrk_array(:,:,:,ispin,io)*wrk_array(:,:,:,ispin,io)) * hvol )
        end do
      end do ! ispin
      f_basis = wrk_array
      
    ! sttpsi <-- f_basis == | lambda > (basis functions)
      sttpsi%rwf = 0d0
      do io=info%io_s,info%io_e
      do ispin=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        if( ix <= dc%nxyz_domain(1) .and. iy <= dc%nxyz_domain(2) .and. iz <= dc%nxyz_domain(3) ) then
          sttpsi%rwf(ix,iy,iz,ispin,io,1,1) = f_basis(ix,iy,iz,ispin,io)
        end if
      end do
      end do
      end do
      end do
      end do
      
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
      call comm_summation(wrk_array,hf,product(lg%num)*nspin*dc%nstate_frag,info%icomm_rko)
      
      deallocate(wrk_array)
      
    end subroutine hpsi_basis
    
    subroutine calc_hamiltonian_matrix
      use communication, only: comm_isend, comm_irecv, comm_wait_all
      implicit none
      integer :: d(3),l(3)
      integer :: itag_send,itag_recv
      integer,dimension(n_halo) :: ireq_send,ireq_recv
      
    ! diagonal block < lambda_{ifrag,io} | H | lambda_{ifrag,jo} >
      allocate(mat_H_local(dc%nstate_frag,dc%nstate_frag,nspin))
      mat_H_local = 0d0
      l = dc%nxyz_domain
      do ispin=1,nspin
      do io=1,n_basis(dc%i_frag,ispin)
      do jo=1,n_basis(dc%i_frag,ispin)
        mat_H_local(io,jo,ispin) = &
        & + sum(f_basis(1:l(1),1:l(2),1:l(3),ispin,io)*hf(1:l(1),1:l(2),1:l(3),ispin,jo)) * hvol
      end do
      end do
      end do
            
    ! off-diagonal block < lambda_{jfrag,jo} | H | lambda_{ifrag,io} >
      if(dc%id_frag==0) then
      ! halo communication for f_basis == | lambda >
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

        ! MPI_ISEND
          itag_send = dc%i_frag
          ireq_send(i_halo) = comm_isend(halo(i_halo)%buf_send,halo(i_halo)%id_dst,itag_send,dc%icomm_tot)
        ! MPI_IRECV
          itag_recv = halo(i_halo)%ifrag_src
          ireq_recv(i_halo) = comm_irecv(halo(i_halo)%buf_recv,halo(i_halo)%id_src,itag_recv,dc%icomm_tot)
        end do ! i_halo=1,n_halo
        call comm_wait_all(ireq_recv)
        call comm_wait_all(ireq_send)
      ! halo(i_halo)%mat_H_local: local Hamiltonian matrix < lambda | H | lambda > (off-diagonal block)
        do i_halo=1,n_halo
          l = halo(i_halo)%length
          d = halo(i_halo)%dsp_recv
          allocate(halo(i_halo)%mat_H_local(dc%nstate_frag,dc%nstate_frag,nspin))
          halo(i_halo)%mat_H_local = 0d0
          do ispin=1,nspin
          do io=1,dc%nstate_frag
          do jo=1,dc%nstate_frag
            do iz=1,l(3)
            do iy=1,l(2)
            do ix=1,l(1)
              halo(i_halo)%mat_H_local(io,jo,ispin) = halo(i_halo)%mat_H_local(io,jo,ispin) &
              & + halo(i_halo)%buf_recv(ix,iy,iz,ispin,io) * hf(d(1)+ix,d(2)+iy,d(3)+iz,ispin,jo) * hvol
            end do
            end do
            end do
          end do
          end do
          end do
          deallocate(halo(i_halo)%buf_send,halo(i_halo)%buf_recv)
          if(dc%id_tot==0) write(*,*) "Halo communication #",i_halo,": done"
        end do
      end if ! dc%id_frag==0
      deallocate(hf)
            
    end subroutine calc_hamiltonian_matrix
    
    subroutine diag_lapack
      use eigen_subdiag_sub, only: eigen_dsyev
      implicit none
      real(8) :: wrk
      real(8),allocatable :: mat_H(:,:),mat_V(:,:)
      
      do ispin=1,nspin
        if(dc%id_tot==0) write(*,*) "lapack diag, #dim=",n_mat(ispin)
        n = n_mat(ispin)
      ! mat_H <-- total Hamiltonian matrix < lambda | H | lambda >
        allocate(mat_H(n,n),mat_V(n,n))
        mat_V = 0d0
        if(dc%id_frag==0) then
          ifrag = dc%i_frag
        ! diagonal block < lambda_{ifrag,io} | H | lambda_{ifrag,jo} >
          do io=1,n_basis(ifrag,ispin) ; i = index_basis(io,ifrag,ispin)
          do jo=1,n_basis(ifrag,ispin) ; j = index_basis(jo,ifrag,ispin)
            mat_V(i,j) = mat_H_local(io,jo,ispin)
          end do
          end do
        ! off-diagonal block < lambda_{jfrag,jo} | H | lambda_{ifrag,io} >
          do i_halo=1,n_halo ; jfrag = halo(i_halo)%ifrag_src ! src fragment (recv)
            do jo=1,n_basis(jfrag,ispin) ; j = index_basis(jo,jfrag,ispin)
            do io=1,n_basis(ifrag,ispin) ; i = index_basis(io,ifrag,ispin)
            ! mat_H_local(jo,io) == < lambda_{jfrag,jo} | H | lambda_{ifrag,io} >
              wrk = 0.5d0* halo(i_halo)%mat_H_local(jo,io,ispin) ! 0.5d0* : for double counting
              mat_V(j,i) = mat_V(j,i) + wrk
              mat_V(i,j) = mat_V(i,j) + wrk
            end do
            end do
          end do
        end if ! dc%id_frag==0
        call comm_summation(mat_V,mat_H,n*n,dc%icomm_tot)
        call eigen_dsyev(mat_H,esp_tot(1:n,ispin),mat_V)
        if(dc%id_frag==0) then
          ifrag = dc%i_frag
          do i=1,dc%nstate_tot
          do jo=1,n_basis(ifrag,ispin) ; j = index_basis(jo,ifrag,ispin)
            coef_wf(jo,i,ispin) = mat_V(j,i) ! coefficients of the wavefunctions
          end do
          end do
        end if
        deallocate(mat_H,mat_V)
      end do ! ispin
      
    end subroutine diag_lapack

#ifdef USE_SCALAPACK
    subroutine diag_chefsi_driver
      use lcfo_diag_chefsi, only: diag_chefsi
      use salmon_global, only: lcfo_diag_chefsi_filter_degree
      implicit none
      integer :: h,frag
      integer, allocatable :: halo_src(:),halo_dst(:)
      integer, allocatable :: halo_root_src(:),halo_dvec(:,:)
      real(8), allocatable :: h_halo(:,:,:,:)

      allocate(halo_src(n_halo),halo_dst(n_halo))
      allocate(halo_root_src(n_halo),halo_dvec(3,n_halo))
      allocate(h_halo(dc%nstate_frag,dc%nstate_frag,nspin,n_halo))
      h_halo = 0d0
      do h=1,n_halo
        halo_src(h) = halo(h)%ifrag_src
        halo_root_src(h) = halo(h)%id_src
        halo_dvec(:,h) = halo(h)%dvec
        halo_dst(h) = 0
        do frag=1,dc%n_frag
          if(id_array(frag)==halo(h)%id_dst) halo_dst(h) = frag
        end do
        if(halo_dst(h)==0) stop "DC-LCFO CheFSI: destination fragment not found."
        if(dc%id_frag==0) h_halo(:,:,:,h) = halo(h)%mat_H_local
      end do

      call diag_chefsi(dc,nspin,lcfo_diag_chefsi_filter_degree,n_basis,n_mat, &
      & n_halo,halo_src, &
      & halo_dst,halo_root_src,halo_dvec,mat_H_local,h_halo, &
      & esp_tot,coef_wf)

      deallocate(h_halo,halo_dvec,halo_root_src,halo_dst,halo_src)
    end subroutine diag_chefsi_driver
#endif

#ifdef USE_EIGENEXA
    subroutine diag_eigenexa
      use communication, only: comm_bcast
      use eigen_libs_mod
      implicit none
      integer :: n,nx,ny,ix_s,ix_e,iy_s,iy_e,ix_loc,iy_loc,ifrag_x,ifrag_y,io_x,io_y
      integer :: nnod,x_nnod,y_nnod,inod,x_inod,y_inod
      integer :: jfrag_halo(n_halo)
      integer, allocatable :: io_array(:),ifrag_array(:)
      real(8), allocatable :: h_div(:,:), v_div(:,:), h(:,:,:), v_tmp1(:,:), v_tmp2(:,:)
      
      allocate(h(dc%nstate_frag,dc%nstate_frag,0:n_halo))
      allocate(v_tmp1(dc%nstate_frag,dc%nstate_tot))
      allocate(v_tmp2(dc%nstate_frag,dc%nstate_tot))
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
        allocate( h_div(nx,ny), v_div(nx,ny) )
        ix_s = eigen_loop_start( 1, x_nnod, x_inod )
        ix_e = eigen_loop_end  ( n, x_nnod, x_inod )
        iy_s = eigen_loop_start( 1, y_nnod, y_inod )
        iy_e = eigen_loop_end  ( n, y_nnod, y_inod )
        
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
            ifrag_y = ifrag_array(iy)
            io_y = io_array(iy)
            do ix_loc=ix_s,ix_e
              ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
              ifrag_x = ifrag_array(ix)
              io_x = io_array(ix)
              if(ifrag_x == ifrag .and. ifrag_y == ifrag) then
                h_div(ix_loc,iy_loc) = h(io_x,io_y,0)
              end if
              do i_halo=1,n_halo
                if( ifrag_x == jfrag_halo(i_halo) .and. ifrag_y == ifrag ) then
                  h_div(ix_loc,iy_loc) = h_div(ix_loc,iy_loc) &
                  & + 0.5d0* h(io_x,io_y,i_halo) ! 0.5d0* : for double counting
                else if( ifrag_x == ifrag .and. ifrag_y == jfrag_halo(i_halo) ) then
                  h_div(ix_loc,iy_loc) = h_div(ix_loc,iy_loc) &
                  & + 0.5d0* h(io_y,io_x,i_halo) ! 0.5d0* : for double counting
                end if
              end do ! i_halo
            end do ! ix_loc
          end do ! iy_loc
        end do ! ifrag
        if(dc%id_tot==0) write(*,*) "h_div: done"
        
        call eigen_sx(n, n, h_div, nx, esp_tot(1:n,ispin), v_div, nx)
        if(dc%id_tot==0) write(*,*) "eigen_sx: done"
        
        do ifrag=1,dc%n_frag
          v_tmp1 = 0d0
          do iy_loc=iy_s,iy_e
            iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
            do ix_loc=ix_s,ix_e
              ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
              ifrag_x = ifrag_array(ix)
              io_x = io_array(ix)
              if(iy <= dc%nstate_tot .and. ifrag_x == ifrag) then
                v_tmp1(io_x,iy) = v_div(ix_loc,iy_loc)
              end if
            end do
          end do
          call comm_summation(v_tmp1,v_tmp2,dc%nstate_frag*dc%nstate_tot,dc%icomm_tot)
          if(ifrag==dc%i_frag .and. dc%id_frag==0) then
            coef_wf(:,:,ispin) = v_tmp2
          end if
        end do ! ifrag
        
        deallocate(h_div,v_div,io_array,ifrag_array)
        call eigen_free()
      end do ! ispin
      
      deallocate(h,v_tmp1,v_tmp2)
    end subroutine diag_eigenexa
#endif

#ifdef USE_SLEPC
    subroutine diag_slepc
      use communication, only: comm_create_group, comm_free_group
      use slepceps
      implicit none
      type(tMat) :: mat_h
      type(tVec) :: eigvec
      type(tEPS) :: eps
      integer(PETSC_INT_KIND) :: n_petsc,nev,n_local,row_start,row_end
      integer(PETSC_INT_KIND) :: nconv,niter,ieig,k,gidx,expected_start
      integer(PETSC_INT_KIND) :: row(1),ncols,col_start,col_end,diag_count
      integer :: ierr
      real(PETSC_REAL_KIND) :: eig_r,eig_i
      real(PETSC_REAL_KIND) :: rel_error,max_error
      integer(PETSC_INT_KIND),allocatable :: d_nnz(:),o_nnz(:),cols(:),eig_indices(:)
      real(PETSC_REAL_KIND),allocatable :: vals(:),eig_values(:)
      integer :: icomm_slepc,key_slepc,n_local_default,local_offset
      integer :: n_before,n_connected,io_local
      integer,allocatable :: adjacency_local(:,:),adjacency(:,:)
      real(8),allocatable :: coef_local(:,:),coef_fragment(:,:)

      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      call check_slepc(ierr,'SlepcInitialize')

      key_slepc = (dc%i_frag-1)*dc%isize_tot + dc%id_frag
      icomm_slepc = comm_create_group(dc%icomm_tot,0,key_slepc)

      allocate(adjacency_local(dc%n_frag,dc%n_frag))
      allocate(adjacency(dc%n_frag,dc%n_frag))
      adjacency_local = 0
      if(dc%id_frag==0) then
        adjacency_local(dc%i_frag,dc%i_frag) = 1
        do i_halo=1,n_halo
          adjacency_local(dc%i_frag,halo(i_halo)%ifrag_src) = 1
        end do
      end if
      call comm_summation(adjacency_local,adjacency,dc%n_frag*dc%n_frag,dc%icomm_tot)
      adjacency = max(adjacency,transpose(adjacency))
      deallocate(adjacency_local)

      allocate(coef_local(dc%nstate_frag,dc%nstate_tot))
      allocate(coef_fragment(dc%nstate_frag,dc%nstate_tot))
      allocate(cols(max(1,dc%nstate_frag)),vals(max(1,dc%nstate_frag)))

      do ispin=1,nspin
        if(dc%id_tot==0) write(*,*) "SLEPc Krylov-Schur diag, #dim=",n_mat(ispin)
        n = n_mat(ispin)
        n_petsc = int(n,PETSC_INT_KIND)
        nev = int(dc%nstate_tot,PETSC_INT_KIND)

        n_local_default = n_basis(dc%i_frag,ispin)/dc%isize_frag
        if(dc%id_frag < mod(n_basis(dc%i_frag,ispin),dc%isize_frag)) then
          n_local_default = n_local_default + 1
        end if
        local_offset = dc%id_frag*(n_basis(dc%i_frag,ispin)/dc%isize_frag) &
        & + min(dc%id_frag,mod(n_basis(dc%i_frag,ispin),dc%isize_frag))
        n_before = sum(n_basis(1:dc%i_frag-1,ispin))
        expected_start = int(n_before+local_offset,PETSC_INT_KIND)
        n_local = int(n_local_default,PETSC_INT_KIND)

        call MatCreate(icomm_slepc,mat_h,ierr)
        call check_slepc(ierr,'MatCreate')
        call MatSetSizes(mat_h,n_local,n_local,n_petsc,n_petsc,ierr)
        call check_slepc(ierr,'MatSetSizes')
        call MatSetType(mat_h,MATMPIAIJ,ierr)
        call check_slepc(ierr,'MatSetType')
        call MatGetOwnershipRange(mat_h,row_start,row_end,ierr)
        call check_slepc(ierr,'MatGetOwnershipRange')
        if(row_start /= expected_start .or. row_end-row_start /= n_local) then
          stop "DC-LCFO: unexpected PETSc row ownership."
        end if

        allocate(d_nnz(n_local_default),o_nnz(n_local_default))
        diag_count = 0_PETSC_INT_KIND
        n_connected = 0
        do jfrag=1,dc%n_frag
          if(adjacency(dc%i_frag,jfrag)==0) cycle
          col_start = int(sum(n_basis(1:jfrag-1,ispin)),PETSC_INT_KIND)
          col_end = col_start + int(n_basis(jfrag,ispin),PETSC_INT_KIND)
          diag_count = diag_count + max(0_PETSC_INT_KIND, &
          & min(col_end,row_end)-max(col_start,row_start))
          n_connected = n_connected + n_basis(jfrag,ispin)
        end do
        d_nnz = diag_count
        o_nnz = int(n_connected,PETSC_INT_KIND)-diag_count
        call MatMPIAIJSetPreallocation(mat_h,0_PETSC_INT_KIND,d_nnz, &
        & 0_PETSC_INT_KIND,o_nnz,ierr)
        call check_slepc(ierr,'MatMPIAIJSetPreallocation')
        deallocate(d_nnz,o_nnz)
        call MatSetOption(mat_h,MAT_SYMMETRIC,PETSC_TRUE,ierr)
        call check_slepc(ierr,'MatSetOption(MAT_SYMMETRIC)')
        call MatSetOption(mat_h,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr)
        call check_slepc(ierr,'MatSetOption(MAT_SYMMETRY_ETERNAL)')

        if(dc%id_frag==0) then
          ifrag = dc%i_frag
          ncols = int(n_basis(ifrag,ispin),PETSC_INT_KIND)
          do io=1,n_basis(ifrag,ispin)
            row(1) = int(index_basis(io,ifrag,ispin)-1,PETSC_INT_KIND)
            do jo=1,n_basis(ifrag,ispin)
              cols(jo) = int(index_basis(jo,ifrag,ispin)-1,PETSC_INT_KIND)
              if(io <= jo) then
                vals(jo) = mat_H_local(io,jo,ispin)
              else
                vals(jo) = mat_H_local(jo,io,ispin)
              end if
            end do
            call MatSetValues(mat_h,1_PETSC_INT_KIND,row,ncols,cols,vals,ADD_VALUES,ierr)
            call check_slepc(ierr,'MatSetValues(diagonal block)')
          end do

          do i_halo=1,n_halo
            jfrag = halo(i_halo)%ifrag_src
            ncols = int(n_basis(jfrag,ispin),PETSC_INT_KIND)
            do io=1,n_basis(ifrag,ispin)
              row(1) = int(index_basis(io,ifrag,ispin)-1,PETSC_INT_KIND)
              do jo=1,n_basis(jfrag,ispin)
                cols(jo) = int(index_basis(jo,jfrag,ispin)-1,PETSC_INT_KIND)
                vals(jo) = 0.5d0*halo(i_halo)%mat_H_local(jo,io,ispin)
              end do
              call MatSetValues(mat_h,1_PETSC_INT_KIND,row,ncols,cols,vals,ADD_VALUES,ierr)
              call check_slepc(ierr,'MatSetValues(off-diagonal block)')
            end do

            ncols = int(n_basis(ifrag,ispin),PETSC_INT_KIND)
            do jo=1,n_basis(jfrag,ispin)
              row(1) = int(index_basis(jo,jfrag,ispin)-1,PETSC_INT_KIND)
              do io=1,n_basis(ifrag,ispin)
                cols(io) = int(index_basis(io,ifrag,ispin)-1,PETSC_INT_KIND)
                vals(io) = 0.5d0*halo(i_halo)%mat_H_local(jo,io,ispin)
              end do
              call MatSetValues(mat_h,1_PETSC_INT_KIND,row,ncols,cols,vals,ADD_VALUES,ierr)
              call check_slepc(ierr,'MatSetValues(transposed off-diagonal block)')
            end do
          end do
        end if

        call MatAssemblyBegin(mat_h,MAT_FINAL_ASSEMBLY,ierr)
        call check_slepc(ierr,'MatAssemblyBegin')
        call MatAssemblyEnd(mat_h,MAT_FINAL_ASSEMBLY,ierr)
        call check_slepc(ierr,'MatAssemblyEnd')

        call EPSCreate(icomm_slepc,eps,ierr)
        call check_slepc(ierr,'EPSCreate')
        call EPSSetOperators(eps,mat_h,PETSC_NULL_MAT,ierr)
        call check_slepc(ierr,'EPSSetOperators')
        call EPSSetProblemType(eps,EPS_HEP,ierr)
        call check_slepc(ierr,'EPSSetProblemType')
        call EPSSetType(eps,EPSKRYLOVSCHUR,ierr)
        call check_slepc(ierr,'EPSSetType')
        call EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL,ierr)
        call check_slepc(ierr,'EPSSetWhichEigenpairs')
        call EPSSetDimensions(eps,nev,PETSC_DETERMINE,PETSC_DETERMINE,ierr)
        call check_slepc(ierr,'EPSSetDimensions')
        call EPSSetOptionsPrefix(eps,'dc_lcfo_',ierr)
        call check_slepc(ierr,'EPSSetOptionsPrefix')
        call EPSSetFromOptions(eps,ierr)
        call check_slepc(ierr,'EPSSetFromOptions')
        call EPSSolve(eps,ierr)
        call check_slepc(ierr,'EPSSolve')
        call EPSGetConverged(eps,nconv,ierr)
        call check_slepc(ierr,'EPSGetConverged')
        if(nconv < nev) then
          if(dc%id_tot==0) write(*,*) "SLEPc converged eigenpairs:",nconv," requested:",nev
          stop "DC-LCFO: SLEPc did not converge all requested eigenpairs."
        end if

        call MatCreateVecs(mat_h,eigvec,PETSC_NULL_VEC,ierr)
        call check_slepc(ierr,'MatCreateVecs')
        allocate(eig_indices(n_local_default),eig_values(n_local_default))
        do k=0_PETSC_INT_KIND,n_local-1_PETSC_INT_KIND
          eig_indices(k+1) = row_start+k
        end do
        coef_local = 0d0
        max_error = 0d0
        do ieig=0_PETSC_INT_KIND,nev-1_PETSC_INT_KIND
          call EPSGetEigenpair(eps,ieig,eig_r,eig_i,eigvec,PETSC_NULL_VEC,ierr)
          call check_slepc(ierr,'EPSGetEigenpair')
          esp_tot(ieig+1,ispin) = real(eig_r,8)
          call VecGetValues(eigvec,n_local,eig_indices,eig_values,ierr)
          call check_slepc(ierr,'VecGetValues')
          do k=0_PETSC_INT_KIND,n_local-1_PETSC_INT_KIND
            gidx = row_start+k
            io_local = int(gidx-int(n_before,PETSC_INT_KIND))+1
            coef_local(io_local,ieig+1) = real(eig_values(k+1),8)
          end do
          call EPSComputeError(eps,ieig,EPS_ERROR_RELATIVE,rel_error,ierr)
          call check_slepc(ierr,'EPSComputeError')
          max_error = max(max_error,rel_error)
        end do
        call comm_summation(coef_local,coef_fragment, &
        & dc%nstate_frag*dc%nstate_tot,dc%icomm_frag)
        if(dc%id_frag==0) coef_wf(:,:,ispin) = coef_fragment

        call EPSGetIterationNumber(eps,niter,ierr)
        call check_slepc(ierr,'EPSGetIterationNumber')
        if(dc%id_tot==0) then
          write(*,*) "SLEPc iterations:",niter," max relative residual:",max_error
        end if
        deallocate(eig_indices,eig_values)
        call VecDestroy(eigvec,ierr)
        call check_slepc(ierr,'VecDestroy')
        call EPSDestroy(eps,ierr)
        call check_slepc(ierr,'EPSDestroy')
        call MatDestroy(mat_h,ierr)
        call check_slepc(ierr,'MatDestroy')
      end do

      deallocate(adjacency,coef_local,coef_fragment,cols,vals)
      call comm_free_group(icomm_slepc)
      call SlepcFinalize(ierr)
      call check_slepc(ierr,'SlepcFinalize')
    end subroutine diag_slepc

    subroutine check_slepc(ierr,where)
      use slepceps
      implicit none
      integer,intent(in) :: ierr
      character(*),intent(in) :: where
      if(ierr /= 0) then
        if(dc%id_tot==0) write(*,*) "SLEPc/PETSc error in ",trim(where),": ",ierr
        stop "DC-LCFO: SLEPc/PETSc call failed."
      end if
    end subroutine check_slepc
#endif
    
    subroutine output
      use salmon_global, only: base_directory, sysname, unit_energy
      use filesystem, only: get_filehandle
      use inputoutput, only: uenergy_from_au
      implicit none
      integer :: iunit,i_halo
      character(256) :: filename
      
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
        write(iunit) dc%nxyz_domain(1:3),nspin,dc%nstate_frag
        write(iunit) n_basis(dc%i_frag,1:nspin) ! # of basis functions
        write(iunit) f_basis(1:dc%nxyz_domain(1),1:dc%nxyz_domain(2),1:dc%nxyz_domain(3) &
        & ,1:nspin,1:dc%nstate_frag) ! basis functions | lambda >
        close(iunit)
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
          close(iunit)
        end if
      end if
      
    end subroutine output
    
!+++++++++++++++++
    subroutine test_write_psi
      use salmon_global, only: natom, kion, rion, base_directory
      use write_file3d
      implicit none
      character(60) :: suffix='./psi_test'
      character(30) :: phys_quantity='psi'
      real(8),dimension(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3)) :: wrk1,wrk2
      integer :: ix_tot,iy_tot,iz_tot
      integer,parameter :: io_out=1025
      character(256) :: dir_tmp
      
      if(n_mat(1) < io_out) return

      wrk1 = 0d0
      if(dc%id_frag==0) then
        ispin = 1
        io = io_out
        jfrag = dc%i_frag
        do jo=1,n_basis(jfrag,ispin) ; j = index_basis(jo,jfrag,ispin)
        do iz=1,dc%nxyz_domain(3); iz_tot = dc%jxyz_tot(iz,3)
        do iy=1,dc%nxyz_domain(2); iy_tot = dc%jxyz_tot(iy,2)
        do ix=1,dc%nxyz_domain(1); ix_tot = dc%jxyz_tot(ix,1)
          wrk1(ix_tot,iy_tot,iz_tot) = wrk1(ix_tot,iy_tot,iz_tot) &
          & + f_basis(ix,iy,iz,ispin,jo) * coef_wf(jo,io,ispin)
        end do
        end do
        end do
        end do
      end if
      call comm_summation(wrk1,wrk2,product(dc%lg_tot%num(1:3)),dc%icomm_tot)
        
    ! override (fragment --> total)
      natom = dc%system_tot%nion
      deallocate(kion,rion)
      allocate(kion(natom),rion(3,natom))
      kion = dc%system_tot%kion
      rion = dc%system_tot%rion
      dir_tmp = base_directory
      base_directory = dc%base_directory

      call write_cube(dc%lg_tot,103,suffix,phys_quantity,wrk2,dc%system_tot)
    
    ! override (total --> fragment)
      natom = system%nion
      deallocate(kion,rion)
      allocate(kion(natom),rion(3,natom))
      kion = system%kion
      rion = system%rion
      base_directory = dir_tmp

    end subroutine test_write_psi
!++++++++++++++++
  
  end subroutine dc_lcfo

!===================================================================================================================================

! yn_conventional_from_dcdft==y : conventional calculation but wavefunctions are reconstructed from DC-LCFO data
  subroutine init_conventional_from_dcdft(lg,mg,system,info,spsi)
    use communication, only: comm_is_root, comm_summation, comm_bcast
    use filesystem, only: get_filehandle
    use salmon_global,only: num_fragment
    use structures
    implicit none
    type(s_rgrid),        intent(in) :: lg,mg
    type(s_dft_system),   intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_orbital)                  :: spsi
    !
    character(32),parameter :: bdir_frag='./data_dcdft/fragments/'
    character(256) :: filename
    integer :: iunit, n_frag, nspin, nstate_frag, nstate_tot
    integer :: i,j,jfrag,ispin,io,jo,ix,iy,iz,ix_tot,iy_tot,iz_tot,n
    integer,dimension(3) :: lgnum_frag,lgnum_tmp,nxyz_domain
    !
    integer,allocatable :: n_mat(:),n_basis(:,:),index_basis(:,:,:),jxyz_tot(:,:)
    real(8),allocatable :: f_basis(:,:,:,:,:),coef_wf(:,:,:),wrk1(:,:,:),wrk2(:,:,:)
    
    nspin = system%nspin
    n_frag = product(num_fragment)
    nstate_tot = 0 ! initial
    
    if(comm_is_root(info%id_rko)) then
      write(*,*) "start init_conventional_from_dcdft"
      write(*,*) "yn_conventional_from_dcdft==y : conventional calculation but wavefunctions are reconstructed from DC-LCFO data"
      write(*,*) "read from ./data_dcdft directory"
    end if
    
  ! read fragment data ./data_dcdft/fragments/000001, 000002, ...  
    if(info%isize_rko < n_frag) stop "yn_conventional_from_dcdft=y: MPI size is too small."
    n = info%isize_rko / n_frag
    jfrag = -1
    do j=1,n_frag
      if( j*n == (info%id_rko+1) ) then
        jfrag = j
        exit
      end if
    end do
    if(jfrag > 0) then ! myrank (info%id_rko) <--> jfrag
    ! coefficients of the wavefunctions
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), jfrag, '/', binfile_wf
      open(iunit,file=filename,form='unformatted',access='stream')
      read(iunit) n_frag, nspin, nstate_frag, nstate_tot
      if( n_frag /= product(num_fragment) .or. nspin /= system%nspin ) stop "data_dcdft: input mismatch"
      allocate(n_mat(nspin))
      allocate(n_basis(n_frag,nspin))
      allocate(index_basis(nstate_frag,n_frag,nspin))
      allocate(coef_wf(nstate_frag,nstate_tot,nspin))
      read(iunit) n_mat(1:nspin)
      read(iunit) n_basis(1:n_frag,1:nspin)
      read(iunit) index_basis(1:nstate_frag,1:n_frag,1:nspin)
      read(iunit) coef_wf(1:nstate_frag,1:nstate_tot,1:nspin)
      close(iunit)    
    ! r-grid index
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), jfrag, '/', binfile_rg
      open(iunit,file=filename,form='unformatted',access='stream')
      read(iunit) lgnum_frag(1:3), lgnum_tmp(1:3)
      if( any( lgnum_tmp /= lg%num ) ) stop "data_dcdft: input mismatch (lg)"
      allocate(jxyz_tot(maxval(lgnum_frag),3))
      do n=1,3 ! x,y,z
        read(iunit) jxyz_tot(1:lgnum_frag(n),n)
      end do
      close(iunit)
    ! basis functions | lambda >
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), jfrag, '/', binfile_bf
      open(iunit,file=filename,form='unformatted',access='stream')
      read(iunit) nxyz_domain(1:3),i,j ! i,j: dummy
      read(iunit) lgnum_tmp(1:nspin) ! dummy
      if(i /= nspin .or. j /= nstate_frag .or. any( lgnum_tmp(1:nspin) /= n_basis(jfrag,1:nspin) ) ) then
        stop "data_dcdft: input mismatch (basis_functions.bin)"
      end if
      allocate(f_basis(1:nxyz_domain(1),1:nxyz_domain(2),1:nxyz_domain(3),1:nspin,1:nstate_frag))
      read(iunit) f_basis(1:nxyz_domain(1),1:nxyz_domain(2),1:nxyz_domain(3),1:nspin,1:nstate_frag)
      close(iunit)
    end if
    
  ! r-grid wavefunctions
    allocate(wrk1(lg%num(1),lg%num(2),lg%num(3)))
    allocate(wrk2(lg%num(1),lg%num(2),lg%num(3)))
    do ispin=1,nspin
    do io=1,system%no
      wrk1 = 0d0
      if(jfrag > 0 .and. io <= nstate_tot) then ! myrank (info%id_rko) <--> jfrag
        do jo=1,n_basis(jfrag,ispin) ; j = index_basis(jo,jfrag,ispin)
        do iz=1,nxyz_domain(3); iz_tot = jxyz_tot(iz,3)
        do iy=1,nxyz_domain(2); iy_tot = jxyz_tot(iy,2)
        do ix=1,nxyz_domain(1); ix_tot = jxyz_tot(ix,1)
          wrk1(ix_tot,iy_tot,iz_tot) = wrk1(ix_tot,iy_tot,iz_tot) &
          & + f_basis(ix,iy,iz,ispin,jo) * coef_wf(jo,io,ispin)
        end do
        end do
        end do
        end do
      end if
      call comm_summation(wrk1,wrk2,product(lg%num(1:3)),info%icomm_rko)
      if(info%io_s <= io .and. io <= info%io_e) then
        if(allocated(spsi%rwf)) then
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            spsi%rwf(ix,iy,iz,ispin,io,1,1) = wrk2(ix,iy,iz)
          end do
          end do
          end do
        else if(allocated(spsi%zwf)) then
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            spsi%zwf(ix,iy,iz,ispin,io,1,1) = dcmplx(wrk2(ix,iy,iz))
          end do
          end do
          end do
        end if
      end if
    end do
    end do

    if(comm_is_root(info%id_rko)) then
      write(*,*) "end init_conventional_from_dcdft"
    end if
    
    if(jfrag > 0) deallocate(n_mat,n_basis,index_basis,jxyz_tot,coef_wf,f_basis)
    deallocate(wrk1,wrk2)
  end subroutine init_conventional_from_dcdft
  
end module lcfo
