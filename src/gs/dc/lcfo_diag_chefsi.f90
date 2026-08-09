!
!  Copyright 2019-2026 SALMON developers
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

module lcfo_diag_chefsi
  implicit none
  private

  integer, parameter :: dense_block_size = 64
  integer, parameter :: max_filter_cycle = 200
  integer, parameter :: lanczos_max_steps = 30
  real(8), parameter :: residual_tolerance = 1d-7
  real(8), parameter :: lanczos_residual_factor = 10d0
  real(8), parameter :: lanczos_relative_margin = 0.1d0
  real(8), parameter :: filter_growth_margin = 1d6

  type :: s_blacs_grid
    integer :: ictxt = -1
    integer :: nprow = 0
    integer :: npcol = 0
    integer :: myrow = -1
    integer :: mycol = -1
  end type s_blacs_grid

  type :: s_matrix_layout
    integer :: desc(9) = 0
    integer :: nrow_local = 0
    integer :: ncol_local = 0
    integer :: lda = 1
  end type s_matrix_layout

  type :: s_hamiltonian_workspace
    integer, allocatable :: request_send(:)
    integer, allocatable :: request_recv(:)
    real(8), allocatable :: x_halo(:,:,:)
  end type s_hamiltonian_workspace

  type :: s_chefsi_workspace
    type(s_hamiltonian_workspace) :: hamiltonian
    real(8), allocatable :: q0(:,:),q1(:,:),q2(:,:),hx(:,:)
    real(8), allocatable :: xd(:,:),hxd(:,:),xrot(:,:)
    real(8), allocatable :: projected(:,:),ritz_vector(:,:)
    real(8), allocatable :: residual_local(:),residual(:)
    real(8), allocatable :: eigen_work(:)
    integer, allocatable :: eigen_iwork(:)
  end type s_chefsi_workspace

  public :: diag_chefsi

contains

  subroutine diag_chefsi(dc,nspin,filter_degree,n_basis,n_mat,n_halo,halo_src, &
  & halo_dst,halo_root_src,halo_dvec,h_diag,h_halo,esp_tot,coef_wf)
    use communication, only: comm_bcast, comm_create_group, &
    & comm_free_group, comm_get_max, comm_isend, comm_irecv, &
    & comm_summation, comm_wait_all
    use eigen_subdiag_sub, only: eigen_dsyev
    use structures, only: s_dcdft
    use timer, only: LOG_CHEFSI_SETUP,LOG_CHEFSI_TOTAL, &
    & timer_begin,timer_end
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: nspin,filter_degree,n_basis(:,:),n_mat(:),n_halo
    integer, intent(in) :: halo_src(:),halo_dst(:),halo_root_src(:)
    integer, intent(in) :: halo_dvec(:,:)
    real(8), intent(in) :: h_diag(:,:,:),h_halo(:,:,:,:)
    real(8), intent(inout) :: esp_tot(:,:)
    real(8), allocatable, intent(inout) :: coef_wf(:,:,:)
    type(s_blacs_grid) :: grid_natural,grid_dense
    integer :: icomm_row,ispin,nsub,nbuffer
    integer :: max_basis,npadded
    real(8), allocatable :: h_row(:,:,:,:),h_diag_sym(:,:,:)

    call timer_begin(LOG_CHEFSI_TOTAL)
    call timer_begin(LOG_CHEFSI_SETUP)
    max_basis = size(h_diag,1)
    npadded = max_basis*dc%n_frag
    nbuffer = ceiling(0.05d0*real(dc%nstate_tot,8))
    nbuffer = min(nbuffer,minval(n_mat)-dc%nstate_tot)
    nbuffer = max(0,nbuffer)
    nsub = dc%nstate_tot+nbuffer

    if(dc%isize_tot /= dc%n_frag*dc%isize_frag) then
      stop "DC-LCFO CheFSI: inconsistent fragment process grid."
    end if
    if(nsub < 1 .or. any(n_mat < nsub)) then
      stop "DC-LCFO CheFSI: invalid subspace dimension."
    end if

    call initialize_blacs_grids(dc,grid_natural,grid_dense)
    icomm_row = comm_create_group(dc%icomm_tot,dc%id_frag,dc%i_frag)

    allocate(h_row(max_basis,max_basis,nspin,n_halo))
    allocate(h_diag_sym(max_basis,max_basis,nspin))
    call prepare_hamiltonian_blocks
    call timer_end(LOG_CHEFSI_SETUP)

    if(dc%id_tot==0) then
      write(*,*) "CheFSI subspace:",nsub," target:",dc%nstate_tot
      write(*,*) "CheFSI polynomial degree:",filter_degree
    end if

    do ispin=1,nspin
      call solve_one_spin(ispin)
    end do

    deallocate(h_diag_sym,h_row)
    call comm_free_group(icomm_row)
    call finalize_blacs_grids(grid_natural,grid_dense)
    call timer_end(LOG_CHEFSI_TOTAL)

  contains

    subroutine prepare_hamiltonian_blocks
      implicit none
      integer :: h,i,j,s,send_tag,recv_tag
      integer, allocatable :: request_send(:),request_recv(:)
      real(8), allocatable :: h_recv(:,:,:,:)

      h_row = 0d0
      h_diag_sym = 0d0
      do s=1,nspin
        do i=1,max_basis
          do j=i,max_basis
            h_diag_sym(i,j,s) = h_diag(i,j,s)
            h_diag_sym(j,i,s) = h_diag(i,j,s)
          end do
        end do
      end do

      allocate(h_recv(max_basis,max_basis,nspin,n_halo))
      h_recv = 0d0
      if(dc%id_frag==0) then
        allocate(request_send(n_halo),request_recv(n_halo))
        do h=1,n_halo
          send_tag = direction_tag(-halo_dvec(:,h))
          recv_tag = direction_tag( halo_dvec(:,h))
          request_send(h) = comm_isend(h_halo(:,:,:,h), &
          & halo_root_src(h),send_tag,dc%icomm_tot)
          request_recv(h) = comm_irecv(h_recv(:,:,:,h), &
          & halo_root_src(h),recv_tag,dc%icomm_tot)
        end do
        call comm_wait_all(request_recv)
        call comm_wait_all(request_send)
        deallocate(request_recv,request_send)

        do h=1,n_halo
          do s=1,nspin
            do j=1,n_basis(halo_src(h),s)
              do i=1,n_basis(dc%i_frag,s)
                h_row(i,j,s,h) = 0.5d0*(h_halo(j,i,s,h) &
                & +h_recv(i,j,s,h))
              end do
            end do
          end do
        end do
      end if
      call comm_bcast(h_row,dc%icomm_frag,0)
      deallocate(h_recv)
    end subroutine prepare_hamiltonian_blocks

    subroutine solve_one_spin(s)
      implicit none
      integer, intent(in) :: s
      type(s_matrix_layout) :: layout_natural,layout_dense,layout_small
      type(s_chefsi_workspace) :: workspace
      integer :: cycle,lanczos_steps,ncol,nstate
      real(8) :: gershgorin_upper,lambda_upper,lambda_cut,max_residual
      real(8) :: lanczos_residual,lanczos_ritz
      real(8), allocatable :: x(:,:),eigenvalue(:)

      nstate = dc%nstate_tot
      call initialize_layouts(npadded,nsub,max_basis,grid_natural, &
      & grid_dense,layout_natural,layout_dense,layout_small)
      ncol = layout_natural%ncol_local
      allocate(x(layout_natural%lda,max(1,ncol)))
      allocate(eigenvalue(nsub))
      call initialize_workspace(layout_natural,layout_dense, &
      & layout_small,nstate,eigenvalue,workspace)
      x = 0d0

      call initialize_subspace(s,layout_natural,x)
      call rayleigh_ritz(s,layout_natural,layout_dense,layout_small, &
      & x,eigenvalue,workspace)
      call calculate_residual(s,layout_natural,x,eigenvalue, &
      & nstate,max_residual,workspace)
      lambda_cut = eigenvalue(nsub)
      call estimate_upper_bound(s,layout_natural,workspace, &
      & eigenvalue(1),lambda_cut,gershgorin_upper,lambda_upper, &
      & lanczos_ritz,lanczos_residual,lanczos_steps)

      if(dc%id_tot==0) then
        write(*,*) "CheFSI diag, #dim=",n_mat(s)
        write(*,*) "CheFSI Gershgorin upper bound:",gershgorin_upper
        write(*,*) "CheFSI Lanczos largest Ritz value:",lanczos_ritz
        write(*,*) "CheFSI Lanczos residual estimate:",lanczos_residual
        write(*,*) "CheFSI Lanczos steps:",lanczos_steps
        write(*,*) "CheFSI selected upper spectral bound:",lambda_upper
        write(*,*) "CheFSI initial maximum residual:",max_residual
      end if

      do cycle=1,max_filter_cycle
        if(max_residual <= residual_tolerance) exit
        if(nsub==n_mat(s)) exit
        lambda_cut = eigenvalue(nsub)
        call chebyshev_filter(s,layout_natural,x,eigenvalue(1), &
        & lambda_cut,lambda_upper,gershgorin_upper,workspace)
        call orthonormalize(layout_natural,layout_dense, &
        & layout_small,x,workspace)
        call rayleigh_ritz(s,layout_natural,layout_dense, &
        & layout_small,x,eigenvalue,workspace)
        call calculate_residual(s,layout_natural,x,eigenvalue, &
        & nstate,max_residual,workspace)
        if(dc%id_tot==0) then
          write(*,*) "CheFSI cycle:",cycle," maximum residual:", &
          & max_residual
        end if
      end do

      if(max_residual > residual_tolerance) then
        if(dc%id_tot==0) write(*,*) &
        & "CheFSI failed to converge; maximum residual:",max_residual
        stop "DC-LCFO CheFSI: eigensolver did not converge."
      end if

      esp_tot(1:nstate,s) = eigenvalue(1:nstate)
      call export_coefficients(s,layout_natural,x,nstate)
      call finalize_workspace(workspace)
      deallocate(eigenvalue,x)
    end subroutine solve_one_spin

    subroutine initialize_workspace(layout_n,layout_d,layout_g,nstate, &
    & eigenvalue,workspace)
      use timer, only: LOG_CHEFSI_SETUP,timer_begin,timer_end
      implicit none
      type(s_matrix_layout), intent(in) :: layout_n,layout_d,layout_g
      integer, intent(in) :: nstate
      real(8), intent(inout) :: eigenvalue(:)
      type(s_chefsi_workspace), intent(out) :: workspace
      integer :: info,liwork,lwork,lwork_min,ncol
      real(8) :: work_query(1)

      call timer_begin(LOG_CHEFSI_SETUP)
      ncol = max(1,layout_n%ncol_local)
      allocate(workspace%hamiltonian%request_send(max(1,n_halo)))
      allocate(workspace%hamiltonian%request_recv(max(1,n_halo)))
      allocate(workspace%hamiltonian%x_halo(max_basis,ncol,max(1,n_halo)))
      allocate(workspace%q0(layout_n%lda,ncol))
      allocate(workspace%q1(layout_n%lda,ncol))
      allocate(workspace%q2(layout_n%lda,ncol))
      allocate(workspace%hx(layout_n%lda,ncol))
      allocate(workspace%xd(layout_d%lda,max(1,layout_d%ncol_local)))
      allocate(workspace%hxd(layout_d%lda,max(1,layout_d%ncol_local)))
      allocate(workspace%xrot(layout_d%lda,max(1,layout_d%ncol_local)))
      allocate(workspace%projected(layout_g%lda, &
      & max(1,layout_g%ncol_local)))
      allocate(workspace%ritz_vector(layout_g%lda, &
      & max(1,layout_g%ncol_local)))
      allocate(workspace%residual_local(nstate),workspace%residual(nstate))

      liwork = max(1,2+7*nsub+8*grid_dense%npcol)
      allocate(workspace%eigen_iwork(liwork))
      call pdsyevd('V','L',nsub,workspace%projected,1,1,layout_g%desc, &
      & eigenvalue,workspace%ritz_vector,1,1,layout_g%desc, &
      & work_query,-1,workspace%eigen_iwork,liwork,info)
      if(info/=0) then
        stop "DC-LCFO CheFSI: PDSYEVD workspace query failed."
      end if
      lwork_min = max(1+6*nsub+2*layout_g%nrow_local* &
      & layout_g%ncol_local,3*nsub+max(dense_block_size* &
      & (layout_g%nrow_local+1),3*dense_block_size))
      lwork = max(10*lwork_min,10*nint(work_query(1)))
      allocate(workspace%eigen_work(lwork))
      call timer_end(LOG_CHEFSI_SETUP)
    end subroutine initialize_workspace

    subroutine finalize_workspace(workspace)
      use timer, only: LOG_CHEFSI_SETUP,timer_begin,timer_end
      implicit none
      type(s_chefsi_workspace), intent(inout) :: workspace

      call timer_begin(LOG_CHEFSI_SETUP)
      deallocate(workspace%eigen_iwork,workspace%eigen_work)
      deallocate(workspace%residual,workspace%residual_local)
      deallocate(workspace%ritz_vector,workspace%projected)
      deallocate(workspace%xrot,workspace%hxd,workspace%xd)
      deallocate(workspace%hx,workspace%q2,workspace%q1,workspace%q0)
      deallocate(workspace%hamiltonian%x_halo)
      deallocate(workspace%hamiltonian%request_recv)
      deallocate(workspace%hamiltonian%request_send)
      call timer_end(LOG_CHEFSI_SETUP)
    end subroutine finalize_workspace

    subroutine initialize_subspace(s,layout,x)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      integer, intent(in) :: s
      type(s_matrix_layout), intent(in) :: layout
      real(8), intent(out) :: x(:,:)
      integer :: i,j,k,nb,offset,owner,local_column,local_state
      integer, allocatable :: order(:)
      real(8), allocatable :: block(:,:),value(:),vector(:,:)
      real(8), allocatable :: value_local(:),value_global(:)

      nb = n_basis(dc%i_frag,s)
      allocate(block(nb,nb),value(nb),vector(nb,nb))
      do j=1,nb
        do i=1,nb
          block(i,j) = h_diag_sym(i,j,s)
        end do
      end do
      call eigen_dsyev(block,value,vector)

      allocate(value_local(n_mat(s)),value_global(n_mat(s)))
      value_local = 0d0
      offset = sum(n_basis(1:dc%i_frag-1,s))
      if(dc%id_frag==0) value_local(offset+1:offset+nb) = value
      call comm_summation(value_local,value_global,n_mat(s),dc%icomm_tot)
      if(any(.not.ieee_is_finite(value_global))) then
        stop "DC-LCFO CheFSI: non-finite initial eigenvalue."
      end if
      allocate(order(n_mat(s)))
      do i=1,n_mat(s)
        order(i) = i
      end do
      call sort_eigenvalue_index(value_global,order,1,n_mat(s))

      x = 0d0
      do k=1,nsub
        if(order(k)<=offset .or. order(k)>offset+nb) cycle
        owner = (k-1)/natural_column_block(nsub,grid_natural%npcol)
        if(owner/=grid_natural%mycol) cycle
        local_column = k-owner*natural_column_block( &
        & nsub,grid_natural%npcol)
        local_state = order(k)-offset
        x(1:nb,local_column) = vector(1:nb,local_state)
      end do

      deallocate(order,value_global,value_local,vector,value,block)
    end subroutine initialize_subspace

    subroutine apply_hamiltonian(s,layout,x,y,workspace)
      use timer, only: LOG_CHEFSI_H_APPLY,LOG_CHEFSI_H_COMM_POST, &
      & LOG_CHEFSI_H_DIAG,LOG_CHEFSI_H_HALO, &
      & LOG_CHEFSI_H_RECV_WAIT,LOG_CHEFSI_H_SEND_WAIT, &
      & timer_begin,timer_end
      implicit none
      integer, intent(in) :: s
      type(s_matrix_layout), intent(in) :: layout
      real(8), intent(in) :: x(:,:)
      real(8), intent(out) :: y(:,:)
      type(s_hamiltonian_workspace), intent(inout) :: workspace
      integer :: h,nb,ncol,source_rank,destination_rank

      call timer_begin(LOG_CHEFSI_H_APPLY)
      nb = n_basis(dc%i_frag,s)
      ncol = layout%ncol_local
      y = 0d0
      if(ncol==0) then
        call timer_end(LOG_CHEFSI_H_APPLY)
        return
      end if

      call timer_begin(LOG_CHEFSI_H_COMM_POST)
      do h=1,n_halo
        source_rank = halo_src(h)-1
        destination_rank = halo_dst(h)-1
        workspace%request_send(h) = comm_isend( &
        & x(:,1:ncol),destination_rank, &
        & dc%i_frag,icomm_row)
        workspace%request_recv(h) = comm_irecv( &
        & workspace%x_halo(:,1:ncol,h),source_rank,halo_src(h), &
        & icomm_row)
      end do
      call timer_end(LOG_CHEFSI_H_COMM_POST)

      call timer_begin(LOG_CHEFSI_H_DIAG)
      call dgemm('N','N',nb,ncol,nb,1d0,h_diag_sym(:,:,s), &
      & max_basis,x,layout%lda,0d0,y,layout%lda)
      call timer_end(LOG_CHEFSI_H_DIAG)
      if(n_halo>0) then
        call timer_begin(LOG_CHEFSI_H_RECV_WAIT)
        call comm_wait_all(workspace%request_recv(1:n_halo))
        call timer_end(LOG_CHEFSI_H_RECV_WAIT)
      end if
      call timer_begin(LOG_CHEFSI_H_HALO)
      do h=1,n_halo
        call dgemm('N','N',nb,ncol,n_basis(halo_src(h),s),1d0, &
        & h_row(:,:,s,h),max_basis,workspace%x_halo(:,:,h),max_basis, &
        & 1d0,y,layout%lda)
      end do
      call timer_end(LOG_CHEFSI_H_HALO)
      if(n_halo>0) then
        call timer_begin(LOG_CHEFSI_H_SEND_WAIT)
        call comm_wait_all(workspace%request_send(1:n_halo))
        call timer_end(LOG_CHEFSI_H_SEND_WAIT)
      end if
      call timer_end(LOG_CHEFSI_H_APPLY)
    end subroutine apply_hamiltonian

    subroutine chebyshev_filter(s,layout,x,lower,cutoff,upper, &
    & fallback_upper,workspace)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use timer, only: LOG_CHEFSI_FILTER,timer_begin,timer_end
      implicit none
      integer, intent(in) :: s
      type(s_matrix_layout), intent(in) :: layout
      real(8), intent(inout) :: x(:,:)
      real(8), intent(in) :: lower,cutoff,fallback_upper
      real(8), intent(inout) :: upper
      type(s_chefsi_workspace), intent(inout) :: workspace
      integer :: attempt,degree,ncol,unstable
      real(8) :: center,log_limit,max_amplitude,radius,scaled_lower

      ncol = layout%ncol_local
      if(ncol==0) return
      call timer_begin(LOG_CHEFSI_FILTER)
      do attempt=1,2
        center = 0.5d0*(upper+cutoff)
        radius = 0.5d0*(upper-cutoff)
        if(radius<=epsilon(1d0)*max(1d0,abs(center))) then
          stop "DC-LCFO CheFSI: invalid filter interval."
        end if

        workspace%q0(:,1:ncol) = x(:,1:ncol)
        call apply_hamiltonian(s,layout,workspace%q0,workspace%hx, &
        & workspace%hamiltonian)
        workspace%q1(:,1:ncol) = (workspace%hx(:,1:ncol) &
        & -center*workspace%q0(:,1:ncol))/radius
        do degree=2,filter_degree
          call apply_hamiltonian(s,layout,workspace%q1,workspace%hx, &
          & workspace%hamiltonian)
          workspace%q2(:,1:ncol) = 2d0*(workspace%hx(:,1:ncol) &
          & -center*workspace%q1(:,1:ncol))/radius &
          & -workspace%q0(:,1:ncol)
          workspace%q0(:,1:ncol) = workspace%q1(:,1:ncol)
          workspace%q1(:,1:ncol) = workspace%q2(:,1:ncol)
        end do

        unstable = 0
        if(any(.not.ieee_is_finite(workspace%q1(:,1:ncol)))) then
          unstable = 1
        else
          max_amplitude = maxval(abs(workspace%q1(:,1:ncol)))
          scaled_lower = abs((lower-center)/radius)
          log_limit = log(filter_growth_margin)
          if(scaled_lower>1d0) then
            log_limit = log_limit+real(filter_degree,8)* &
            & log(scaled_lower+sqrt(scaled_lower**2-1d0))
          end if
          log_limit = min(log(huge(1d0))-log(10d0),log_limit)
          if(max_amplitude>0d0) then
            if(log(max_amplitude)>log_limit) unstable = 1
          end if
        end if
        call comm_get_max(unstable,dc%icomm_tot)
        if(unstable==0) exit
        if(upper>=fallback_upper*(1d0-10d0*epsilon(1d0))) then
          stop "DC-LCFO CheFSI: unstable Chebyshev recurrence."
        end if
        upper = fallback_upper
        if(dc%id_tot==0) write(*,*) &
        & "CheFSI filter fallback to Gershgorin upper bound:",upper
      end do
      x(:,1:ncol) = workspace%q1(:,1:ncol)
      call timer_end(LOG_CHEFSI_FILTER)
    end subroutine chebyshev_filter

    subroutine orthonormalize(layout_n,layout_d,layout_g,x,workspace)
      use timer, only: LOG_CHEFSI_ORTHO,timer_begin,timer_end
      implicit none
      type(s_matrix_layout), intent(in) :: layout_n,layout_d,layout_g
      real(8), intent(inout) :: x(:,:)
      type(s_chefsi_workspace), intent(inout) :: workspace
      integer :: info,pass

      call timer_begin(LOG_CHEFSI_ORTHO)
      workspace%xd = 0d0
      call redistribute_to_dense(layout_n,x,layout_d,workspace%xd)
      do pass=1,2
        workspace%projected = 0d0
        call pdsyrk('L','T',nsub,npadded,1d0,workspace%xd,1,1, &
        & layout_d%desc,0d0,workspace%projected,1,1,layout_g%desc)
        call pdpotrf('L',nsub,workspace%projected,1,1, &
        & layout_g%desc,info)
        if(info/=0) then
          call distributed_qr(layout_d,workspace%xd)
          exit
        end if
        call pdtrsm('R','L','T','N',npadded,nsub,1d0, &
        & workspace%projected,1,1,layout_g%desc,workspace%xd, &
        & 1,1,layout_d%desc)
      end do
      call redistribute_to_natural(layout_d,workspace%xd,layout_n,x)
      call timer_end(LOG_CHEFSI_ORTHO)
    end subroutine orthonormalize

    subroutine distributed_qr(layout,xd)
      implicit none
      type(s_matrix_layout), intent(in) :: layout
      real(8), intent(inout) :: xd(:,:)
      integer :: info,lwork
      real(8) :: work_query(1)
      real(8), allocatable :: tau(:),work(:)

      allocate(tau(max(1,layout%ncol_local+dense_block_size)))
      call pdgeqrf(npadded,nsub,xd,1,1,layout%desc,tau, &
      & work_query,-1,info)
      if(info/=0) stop "DC-LCFO CheFSI: PDGEQRF workspace query failed."
      lwork = max(1,nint(work_query(1)))
      deallocate(tau)
      allocate(tau(max(1,layout%ncol_local+dense_block_size)))
      allocate(work(lwork))
      call pdgeqrf(npadded,nsub,xd,1,1,layout%desc,tau, &
      & work,lwork,info)
      if(info/=0) stop "DC-LCFO CheFSI: PDGEQRF failed."
      deallocate(work)

      call pdorgqr(npadded,nsub,nsub,xd,1,1,layout%desc,tau, &
      & work_query,-1,info)
      if(info/=0) stop "DC-LCFO CheFSI: PDORGQR workspace query failed."
      lwork = max(1,nint(work_query(1)))
      allocate(work(lwork))
      call pdorgqr(npadded,nsub,nsub,xd,1,1,layout%desc,tau, &
      & work,lwork,info)
      if(info/=0) stop "DC-LCFO CheFSI: PDORGQR failed."
      deallocate(work,tau)
    end subroutine distributed_qr

    subroutine rayleigh_ritz(s,layout_n,layout_d,layout_g,x,eigenvalue, &
    & workspace)
      use timer, only: LOG_CHEFSI_PROJECT,LOG_CHEFSI_RAYLEIGH_RITZ, &
      & LOG_CHEFSI_ROTATE,timer_begin,timer_end
      implicit none
      integer, intent(in) :: s
      type(s_matrix_layout), intent(in) :: layout_n,layout_d,layout_g
      real(8), intent(inout) :: x(:,:)
      real(8), intent(out) :: eigenvalue(:)
      type(s_chefsi_workspace), intent(inout) :: workspace

      call timer_begin(LOG_CHEFSI_RAYLEIGH_RITZ)
      call apply_hamiltonian(s,layout_n,x,workspace%hx, &
      & workspace%hamiltonian)
      workspace%xd = 0d0
      workspace%hxd = 0d0
      call redistribute_to_dense(layout_n,x,layout_d,workspace%xd)
      call redistribute_to_dense(layout_n,workspace%hx,layout_d, &
      & workspace%hxd)
      workspace%projected = 0d0
      call timer_begin(LOG_CHEFSI_PROJECT)
      call pdgemm('T','N',nsub,nsub,npadded,1d0,workspace%xd,1,1, &
      & layout_d%desc,workspace%hxd,1,1,layout_d%desc,0d0, &
      & workspace%projected,1,1,layout_g%desc)
      call timer_end(LOG_CHEFSI_PROJECT)
      call diagonalize_projected(layout_g,workspace%projected, &
      & eigenvalue,workspace%ritz_vector,workspace%eigen_work, &
      & workspace%eigen_iwork)
      workspace%xrot = 0d0
      call timer_begin(LOG_CHEFSI_ROTATE)
      call pdgemm('N','N',npadded,nsub,nsub,1d0,workspace%xd,1,1, &
      & layout_d%desc,workspace%ritz_vector,1,1,layout_g%desc,0d0, &
      & workspace%xrot,1,1,layout_d%desc)
      call timer_end(LOG_CHEFSI_ROTATE)
      call redistribute_to_natural(layout_d,workspace%xrot,layout_n,x)
      call timer_end(LOG_CHEFSI_RAYLEIGH_RITZ)
    end subroutine rayleigh_ritz

    subroutine diagonalize_projected(layout,projected,eigenvalue, &
    & vector,work,iwork)
      use timer, only: LOG_CHEFSI_PROJECT_EIGEN,timer_begin,timer_end
      implicit none
      type(s_matrix_layout), intent(in) :: layout
      real(8), intent(inout) :: projected(:,:)
      real(8), intent(out) :: eigenvalue(:),vector(:,:)
      real(8), intent(inout) :: work(:)
      integer, intent(inout) :: iwork(:)
      integer :: info

      call timer_begin(LOG_CHEFSI_PROJECT_EIGEN)
      call pdsyevd('V','L',nsub,projected,1,1,layout%desc, &
      & eigenvalue,vector,1,1,layout%desc,work,size(work), &
      & iwork,size(iwork),info)
      if(info/=0) then
        if(dc%id_tot==0) write(*,*) "PDSYEVD error:",info
        stop "DC-LCFO CheFSI: projected eigensolver failed."
      end if
      call timer_end(LOG_CHEFSI_PROJECT_EIGEN)
    end subroutine diagonalize_projected

    subroutine calculate_residual(s,layout,x,eigenvalue,nstate, &
    & max_residual,workspace)
      use timer, only: LOG_CHEFSI_RESIDUAL,timer_begin,timer_end
      implicit none
      integer, intent(in) :: s,nstate
      type(s_matrix_layout), intent(in) :: layout
      real(8), intent(in) :: x(:,:),eigenvalue(:)
      real(8), intent(out) :: max_residual
      type(s_chefsi_workspace), intent(inout) :: workspace
      integer :: global_column,local_column,nb

      call timer_begin(LOG_CHEFSI_RESIDUAL)
      call apply_hamiltonian(s,layout,x,workspace%hx, &
      & workspace%hamiltonian)
      workspace%residual_local = 0d0
      nb = n_basis(dc%i_frag,s)
      do local_column=1,layout%ncol_local
        global_column = natural_global_column(local_column, &
        & grid_natural%mycol,nsub,grid_natural%npcol)
        if(global_column>nstate) cycle
        workspace%residual_local(global_column) = sum(( &
        & workspace%hx(1:nb,local_column) &
        & -eigenvalue(global_column)*x(1:nb,local_column))**2)
      end do
      call comm_summation(workspace%residual_local,workspace%residual, &
      & nstate,dc%icomm_tot)
      do global_column=1,nstate
        workspace%residual(global_column) = sqrt(max(0d0, &
        & workspace%residual(global_column))) &
        & /max(1d0,abs(eigenvalue(global_column)))
      end do
      max_residual = maxval(workspace%residual)
      call timer_end(LOG_CHEFSI_RESIDUAL)
    end subroutine calculate_residual

    subroutine estimate_upper_bound(s,layout,workspace,lower,cutoff, &
    & gershgorin_upper,upper,lanczos_ritz,lanczos_residual,nsteps)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use timer, only: LOG_CHEFSI_LANCZOS,timer_begin,timer_end
      implicit none
      integer, intent(in) :: s
      type(s_matrix_layout), intent(in) :: layout
      type(s_chefsi_workspace), intent(inout) :: workspace
      real(8), intent(in) :: lower,cutoff
      real(8), intent(out) :: gershgorin_upper,upper
      real(8), intent(out) :: lanczos_ritz,lanczos_residual
      integer, intent(out) :: nsteps
      real(8) :: candidate

      call timer_begin(LOG_CHEFSI_LANCZOS)
      gershgorin_upper = gershgorin_upper_bound(s)
      if(nsub==n_mat(s)) then
        lanczos_ritz = gershgorin_upper
        lanczos_residual = 0d0
        nsteps = 0
        upper = gershgorin_upper
        call timer_end(LOG_CHEFSI_LANCZOS)
        return
      end if

      call lanczos_upper_bound(s,layout,workspace,lanczos_ritz, &
      & lanczos_residual,nsteps)
      candidate = lanczos_ritz+max( &
      & lanczos_residual_factor*lanczos_residual, &
      & lanczos_relative_margin*max(1d0,lanczos_ritz-lower))
      if(.not.ieee_is_finite(candidate) .or. &
      & candidate<=cutoff+epsilon(1d0)*max(1d0,abs(cutoff))) then
        upper = gershgorin_upper
      else
        upper = min(gershgorin_upper,candidate)
      end if
      call timer_end(LOG_CHEFSI_LANCZOS)
    end subroutine estimate_upper_bound

    function gershgorin_upper_bound(s) result(upper)
      implicit none
      integer, intent(in) :: s
      integer :: h,i,nb
      real(8) :: local_upper,upper
      real(8) :: send_value(1),recv_value(1)

      nb = n_basis(dc%i_frag,s)
      local_upper = -huge(1d0)
      do i=1,nb
        send_value(1) = h_diag_sym(i,i,s) &
        & +sum(abs(h_diag_sym(i,1:nb,s))) &
        & -abs(h_diag_sym(i,i,s))
        do h=1,n_halo
          send_value(1) = send_value(1) &
          & +sum(abs(h_row(i,1:n_basis(halo_src(h),s),s,h)))
        end do
        local_upper = max(local_upper,send_value(1))
      end do
      send_value(1) = local_upper
      call comm_get_max(send_value,recv_value,1,dc%icomm_tot)
      upper = recv_value(1)+max(1d-10,1d-8*abs(recv_value(1)))
    end function gershgorin_upper_bound

    subroutine lanczos_upper_bound(s,layout,workspace,ritz,residual, &
    & nsteps)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      integer, intent(in) :: s
      type(s_matrix_layout), intent(in) :: layout
      type(s_chefsi_workspace), intent(inout) :: workspace
      real(8), intent(out) :: ritz,residual
      integer, intent(out) :: nsteps
      type(s_matrix_layout) :: vector_layout
      integer :: global_row,i,j,k,nb,offset,pass
      integer(8) :: hash_value
      real(8) :: beta,global_value,local_value
      real(8) :: send_value(2),recv_value(2)
      real(8), allocatable :: alpha(:),beta_value(:),basis(:,:)
      real(8), allocatable :: eigenvalue(:)
      real(8), allocatable :: overlap_global(:),overlap_local(:)
      real(8), allocatable :: tridiagonal(:,:),vector(:,:)

      nb = n_basis(dc%i_frag,s)
      offset = sum(n_basis(1:dc%i_frag-1,s))
      vector_layout = layout
      vector_layout%ncol_local = 1
      allocate(alpha(lanczos_max_steps))
      allocate(beta_value(lanczos_max_steps))
      allocate(basis(layout%lda,lanczos_max_steps))
      allocate(overlap_local(lanczos_max_steps))
      allocate(overlap_global(lanczos_max_steps))
      alpha = 0d0
      beta_value = 0d0
      basis = 0d0
      do i=1,nb
        global_row = offset+i
        hash_value = modulo(int(global_row,8)+ &
        & int(grid_natural%mycol,8)*int(n_mat(s),8)+12345_8,1048573_8)
        hash_value = modulo(104729_8*hash_value**2+12345_8,1048573_8)
        basis(i,1) = real(hash_value,8)/1048573d0-0.5d0
      end do
      local_value = sum(basis(1:nb,1)**2)
      call comm_summation(local_value,global_value,icomm_row)
      if(global_value<=tiny(1d0)) then
        basis = 0d0
        if(offset==0) basis(1,1) = 1d0
        global_value = 1d0
      end if
      basis(1:nb,1) = basis(1:nb,1)/sqrt(global_value)

      beta = 0d0
      nsteps = lanczos_max_steps
      do k=1,lanczos_max_steps
        call apply_hamiltonian(s,vector_layout,basis(:,k:k), &
        & workspace%hx(:,1:1),workspace%hamiltonian)
        if(k>1) then
          workspace%hx(1:nb,1) = workspace%hx(1:nb,1) &
          & -beta*basis(1:nb,k-1)
        end if
        local_value = dot_product(basis(1:nb,k), &
        & workspace%hx(1:nb,1))
        call comm_summation(local_value,alpha(k),icomm_row)
        workspace%hx(1:nb,1) = workspace%hx(1:nb,1) &
        & -alpha(k)*basis(1:nb,k)

        do pass=1,2
          overlap_local = 0d0
          do j=1,k
            overlap_local(j) = dot_product(basis(1:nb,j), &
            & workspace%hx(1:nb,1))
          end do
          call comm_summation(overlap_local(1:k), &
          & overlap_global(1:k),k,icomm_row)
          do j=1,k
            workspace%hx(1:nb,1) = workspace%hx(1:nb,1) &
            & -overlap_global(j)*basis(1:nb,j)
          end do
        end do
        local_value = sum(workspace%hx(1:nb,1)**2)
        call comm_summation(local_value,global_value,icomm_row)
        beta = sqrt(max(0d0,global_value))
        beta_value(k) = beta
        if(k==lanczos_max_steps .or. beta<=100d0*epsilon(1d0)) then
          nsteps = k
          exit
        end if
        basis(1:nb,k+1) = workspace%hx(1:nb,1)/beta
      end do

      allocate(tridiagonal(nsteps,nsteps),eigenvalue(nsteps))
      allocate(vector(nsteps,nsteps))
      tridiagonal = 0d0
      do k=1,nsteps
        tridiagonal(k,k) = alpha(k)
      end do
      do k=1,nsteps-1
        tridiagonal(k,k+1) = beta_value(k)
        tridiagonal(k+1,k) = beta_value(k)
      end do
      call eigen_dsyev(tridiagonal,eigenvalue,vector)
      ritz = eigenvalue(nsteps)
      residual = beta*abs(vector(nsteps,nsteps))
      if(.not.ieee_is_finite(ritz) .or. &
      & .not.ieee_is_finite(residual)) then
        ritz = -huge(1d0)
        residual = huge(1d0)
      end if
      send_value = [ritz,residual]
      call comm_get_max(send_value,recv_value,2,dc%icomm_tot)
      ritz = recv_value(1)
      residual = recv_value(2)
      deallocate(vector,eigenvalue,tridiagonal)
      deallocate(overlap_global,overlap_local,basis,beta_value,alpha)
    end subroutine lanczos_upper_bound

    subroutine export_coefficients(s,layout,x,nstate)
      use timer, only: LOG_CHEFSI_EXPORT,timer_begin,timer_end
      implicit none
      integer, intent(in) :: s,nstate
      type(s_matrix_layout), intent(in) :: layout
      real(8), intent(in) :: x(:,:)
      integer :: global_column,local_column,nb
      real(8), allocatable :: local_coef(:,:),fragment_coef(:,:)

      call timer_begin(LOG_CHEFSI_EXPORT)
      nb = n_basis(dc%i_frag,s)
      allocate(local_coef(max_basis,nstate),fragment_coef(max_basis,nstate))
      local_coef = 0d0
      do local_column=1,layout%ncol_local
        global_column = natural_global_column(local_column, &
        & grid_natural%mycol,nsub,grid_natural%npcol)
        if(global_column>nstate) cycle
        local_coef(1:nb,global_column) = x(1:nb,local_column)
      end do
      call comm_summation(local_coef,fragment_coef,max_basis*nstate, &
      & dc%icomm_frag)
      if(dc%id_frag==0) coef_wf(:,:,s) = fragment_coef
      deallocate(fragment_coef,local_coef)
      call timer_end(LOG_CHEFSI_EXPORT)
    end subroutine export_coefficients

    subroutine redistribute_to_dense(layout_n,xn,layout_d,xd)
      use timer, only: LOG_CHEFSI_REDISTRIBUTE,timer_begin,timer_end
      implicit none
      type(s_matrix_layout), intent(in) :: layout_n,layout_d
      real(8), intent(in) :: xn(:,:)
      real(8), intent(inout) :: xd(:,:)
      call timer_begin(LOG_CHEFSI_REDISTRIBUTE)
      call pdgemr2d(npadded,nsub,xn,1,1,layout_n%desc,xd,1,1, &
      & layout_d%desc,grid_natural%ictxt)
      call timer_end(LOG_CHEFSI_REDISTRIBUTE)
    end subroutine redistribute_to_dense

    subroutine redistribute_to_natural(layout_d,xd,layout_n,xn)
      use timer, only: LOG_CHEFSI_REDISTRIBUTE,timer_begin,timer_end
      implicit none
      type(s_matrix_layout), intent(in) :: layout_d,layout_n
      real(8), intent(in) :: xd(:,:)
      real(8), intent(inout) :: xn(:,:)
      call timer_begin(LOG_CHEFSI_REDISTRIBUTE)
      call pdgemr2d(npadded,nsub,xd,1,1,layout_d%desc,xn,1,1, &
      & layout_n%desc,grid_natural%ictxt)
      call timer_end(LOG_CHEFSI_REDISTRIBUTE)
    end subroutine redistribute_to_natural

  end subroutine diag_chefsi

  subroutine initialize_blacs_grids(dc,natural,dense)
    use communication, only: comm_summation
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    type(s_blacs_grid), intent(out) :: natural,dense
    integer :: iam,nprocs,i,j
    integer, allocatable :: map_local(:,:),map_global(:,:)
    integer, allocatable :: rank_local(:),rank_global(:),dense_map(:,:)

    call blacs_pinfo(iam,nprocs)
    if(nprocs<dc%isize_tot) then
      stop "DC-LCFO CheFSI: BLACS process count is too small."
    end if

    natural%nprow = dc%n_frag
    natural%npcol = dc%isize_frag
    allocate(map_local(natural%nprow,natural%npcol))
    allocate(map_global(natural%nprow,natural%npcol))
    map_local = 0
    map_local(dc%i_frag,dc%id_frag+1) = iam+1
    call comm_summation(map_local,map_global,size(map_local),dc%icomm_tot)
    map_global = map_global-1
    call blacs_get(0,0,natural%ictxt)
    call blacs_gridmap(natural%ictxt,map_global,natural%nprow, &
    & natural%nprow,natural%npcol)
    call blacs_gridinfo(natural%ictxt,natural%nprow,natural%npcol, &
    & natural%myrow,natural%mycol)
    if(natural%myrow/=dc%i_frag-1 .or. &
    & natural%mycol/=dc%id_frag) then
      stop "DC-LCFO CheFSI: unexpected natural BLACS coordinates."
    end if
    deallocate(map_global,map_local)

    dense%nprow = int(sqrt(real(dc%isize_tot,8)))
    do while(dense%nprow>1)
      if(mod(dc%isize_tot,dense%nprow)==0) exit
      dense%nprow = dense%nprow-1
    end do
    dense%npcol = dc%isize_tot/dense%nprow
    allocate(rank_local(dc%isize_tot),rank_global(dc%isize_tot))
    rank_local = 0
    rank_local(dc%id_tot+1) = iam+1
    call comm_summation(rank_local,rank_global,dc%isize_tot,dc%icomm_tot)
    rank_global = rank_global-1
    allocate(dense_map(dense%nprow,dense%npcol))
    do j=1,dense%npcol
      do i=1,dense%nprow
        dense_map(i,j) = rank_global((j-1)*dense%nprow+i)
      end do
    end do
    call blacs_get(0,0,dense%ictxt)
    call blacs_gridmap(dense%ictxt,dense_map,dense%nprow, &
    & dense%nprow,dense%npcol)
    call blacs_gridinfo(dense%ictxt,dense%nprow,dense%npcol, &
    & dense%myrow,dense%mycol)
    deallocate(dense_map,rank_global,rank_local)
  end subroutine initialize_blacs_grids

  subroutine finalize_blacs_grids(natural,dense)
    implicit none
    type(s_blacs_grid), intent(inout) :: natural,dense
    call blacs_gridexit(dense%ictxt)
    call blacs_gridexit(natural%ictxt)
    dense%ictxt = -1
    natural%ictxt = -1
  end subroutine finalize_blacs_grids

  subroutine initialize_layouts(npadded,nsub,max_basis,natural,dense, &
  & layout_natural,layout_dense,layout_small)
    implicit none
    integer, intent(in) :: npadded,nsub,max_basis
    type(s_blacs_grid), intent(in) :: natural,dense
    type(s_matrix_layout), intent(out) :: layout_natural
    type(s_matrix_layout), intent(out) :: layout_dense,layout_small
    integer :: info,nb_natural
    integer :: numroc

    nb_natural = natural_column_block(nsub,natural%npcol)
    layout_natural%nrow_local = numroc(npadded,max_basis, &
    & natural%myrow,0,natural%nprow)
    layout_natural%ncol_local = numroc(nsub,nb_natural, &
    & natural%mycol,0,natural%npcol)
    layout_natural%lda = max(1,layout_natural%nrow_local)
    call descinit(layout_natural%desc,npadded,nsub,max_basis, &
    & nb_natural,0,0,natural%ictxt,layout_natural%lda,info)
    if(info/=0) stop "DC-LCFO CheFSI: natural DESCINIT failed."

    layout_dense%nrow_local = numroc(npadded,dense_block_size, &
    & dense%myrow,0,dense%nprow)
    layout_dense%ncol_local = numroc(nsub,dense_block_size, &
    & dense%mycol,0,dense%npcol)
    layout_dense%lda = max(1,layout_dense%nrow_local)
    call descinit(layout_dense%desc,npadded,nsub,dense_block_size, &
    & dense_block_size,0,0,dense%ictxt,layout_dense%lda,info)
    if(info/=0) stop "DC-LCFO CheFSI: dense DESCINIT failed."

    layout_small%nrow_local = numroc(nsub,dense_block_size, &
    & dense%myrow,0,dense%nprow)
    layout_small%ncol_local = numroc(nsub,dense_block_size, &
    & dense%mycol,0,dense%npcol)
    layout_small%lda = max(1,layout_small%nrow_local)
    call descinit(layout_small%desc,nsub,nsub,dense_block_size, &
    & dense_block_size,0,0,dense%ictxt,layout_small%lda,info)
    if(info/=0) stop "DC-LCFO CheFSI: projected DESCINIT failed."
  end subroutine initialize_layouts

  integer function natural_column_block(nsub,npcol)
    implicit none
    integer, intent(in) :: nsub,npcol
    natural_column_block = max(1,(nsub+npcol-1)/npcol)
  end function natural_column_block

  integer function natural_global_column(local_column,mycol,nsub,npcol)
    implicit none
    integer, intent(in) :: local_column,mycol,nsub,npcol
    natural_global_column = mycol*natural_column_block(nsub,npcol) &
    & +local_column
  end function natural_global_column

  integer function direction_tag(dvec)
    implicit none
    integer, intent(in) :: dvec(3)
    direction_tag = (dvec(1)+1)*9+(dvec(2)+1)*3+dvec(3)+2
  end function direction_tag

  recursive subroutine sort_eigenvalue_index(value,index,left,right)
    implicit none
    real(8), intent(inout) :: value(:)
    integer, intent(inout) :: index(:)
    integer, intent(in) :: left,right
    integer :: i,j,index_tmp
    real(8) :: pivot,value_tmp

    if(left>=right) return
    pivot = value((left+right)/2)
    i = left
    j = right
    do
      do
        if(i>right) exit
        if(value(i)>=pivot) exit
        i = i+1
      end do
      do
        if(j<left) exit
        if(value(j)<=pivot) exit
        j = j-1
      end do
      if(i>j) exit
      value_tmp = value(i)
      value(i) = value(j)
      value(j) = value_tmp
      index_tmp = index(i)
      index(i) = index(j)
      index(j) = index_tmp
      i = i+1
      j = j-1
    end do
    if(left<j) call sort_eigenvalue_index(value,index,left,j)
    if(i<right) call sort_eigenvalue_index(value,index,i,right)
  end subroutine sort_eigenvalue_index

end module lcfo_diag_chefsi
