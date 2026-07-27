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
  integer, parameter :: filter_degree = 20
  integer, parameter :: max_filter_cycle = 8
  real(8), parameter :: residual_tolerance = 1d-8

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

  public :: diag_chefsi

contains

  subroutine diag_chefsi(dc,nspin,n_basis,n_mat,n_halo,halo_src, &
  & halo_dst,halo_root_src,halo_dvec,h_diag,h_halo,esp_tot,coef_wf)
    use communication, only: comm_bcast, comm_create_group, &
    & comm_free_group, comm_get_max, comm_isend, comm_irecv, &
    & comm_summation, comm_wait_all
    use eigen_subdiag_sub, only: eigen_dsyev
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: nspin,n_basis(:,:),n_mat(:),n_halo
    integer, intent(in) :: halo_src(:),halo_dst(:),halo_root_src(:)
    integer, intent(in) :: halo_dvec(:,:)
    real(8), intent(in) :: h_diag(:,:,:),h_halo(:,:,:,:)
    real(8), intent(inout) :: esp_tot(:,:)
    real(8), allocatable, intent(inout) :: coef_wf(:,:,:)
    type(s_blacs_grid) :: grid_natural,grid_dense
    integer :: icomm_row,ispin,nsub,nbuffer
    integer :: max_basis,npadded
    real(8), allocatable :: h_row(:,:,:,:),h_diag_sym(:,:,:)

    max_basis = size(h_diag,1)
    npadded = max_basis*dc%n_frag
    nbuffer = min(256,max(8,dc%nstate_tot/100))
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
      integer :: cycle,ncol,nstate
      real(8) :: lambda_upper,lambda_cut,max_residual
      real(8), allocatable :: x(:,:),eigenvalue(:)

      nstate = dc%nstate_tot
      call initialize_layouts(npadded,nsub,max_basis,grid_natural, &
      & grid_dense,layout_natural,layout_dense,layout_small)
      ncol = layout_natural%ncol_local
      allocate(x(layout_natural%lda,max(1,ncol)))
      allocate(eigenvalue(nsub))
      x = 0d0

      call initialize_subspace(s,layout_natural,x)
      call rayleigh_ritz(s,layout_natural,layout_dense,layout_small, &
      & x,eigenvalue)
      call calculate_residual(s,layout_natural,x,eigenvalue, &
      & nstate,max_residual)
      lambda_upper = estimate_upper_bound(s)

      if(dc%id_tot==0) then
        write(*,*) "CheFSI diag, #dim=",n_mat(s)
        write(*,*) "CheFSI upper spectral bound:",lambda_upper
        write(*,*) "CheFSI initial maximum residual:",max_residual
      end if

      do cycle=1,max_filter_cycle
        if(max_residual <= residual_tolerance) exit
        if(nsub==n_mat(s)) exit
        lambda_cut = eigenvalue(nsub)
        call chebyshev_filter(s,layout_natural,x,lambda_cut, &
        & lambda_upper)
        call orthonormalize(layout_natural,layout_dense, &
        & layout_small,x)
        call rayleigh_ritz(s,layout_natural,layout_dense, &
        & layout_small,x,eigenvalue)
        call calculate_residual(s,layout_natural,x,eigenvalue, &
        & nstate,max_residual)
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
      deallocate(eigenvalue,x)
    end subroutine solve_one_spin

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

    subroutine apply_hamiltonian(s,layout,x,y)
      implicit none
      integer, intent(in) :: s
      type(s_matrix_layout), intent(in) :: layout
      real(8), intent(in) :: x(:,:)
      real(8), intent(out) :: y(:,:)
      integer :: h,nb,ncol,source_rank,destination_rank
      integer, allocatable :: request_send(:),request_recv(:)
      real(8), allocatable :: x_halo(:,:,:)

      nb = n_basis(dc%i_frag,s)
      ncol = layout%ncol_local
      y = 0d0
      if(ncol==0) return

      allocate(x_halo(max_basis,ncol,n_halo))
      allocate(request_send(n_halo),request_recv(n_halo))
      do h=1,n_halo
        source_rank = halo_src(h)-1
        destination_rank = halo_dst(h)-1
        request_send(h) = comm_isend(x(:,1:ncol),destination_rank, &
        & dc%i_frag,icomm_row)
        request_recv(h) = comm_irecv(x_halo(:,:,h),source_rank, &
        & halo_src(h),icomm_row)
      end do

      call dgemm('N','N',nb,ncol,nb,1d0,h_diag_sym(:,:,s), &
      & max_basis,x,layout%lda,0d0,y,layout%lda)
      call comm_wait_all(request_recv)
      do h=1,n_halo
        call dgemm('N','N',nb,ncol,n_basis(halo_src(h),s),1d0, &
        & h_row(:,:,s,h),max_basis,x_halo(:,:,h),max_basis, &
        & 1d0,y,layout%lda)
      end do
      call comm_wait_all(request_send)
      deallocate(request_recv,request_send,x_halo)
    end subroutine apply_hamiltonian

    subroutine chebyshev_filter(s,layout,x,cutoff,upper)
      implicit none
      integer, intent(in) :: s
      type(s_matrix_layout), intent(in) :: layout
      real(8), intent(inout) :: x(:,:)
      real(8), intent(in) :: cutoff,upper
      integer :: degree,ncol
      real(8) :: center,radius
      real(8), allocatable :: q0(:,:),q1(:,:),q2(:,:),work(:,:)

      ncol = layout%ncol_local
      if(ncol==0) return
      center = 0.5d0*(upper+cutoff)
      radius = 0.5d0*(upper-cutoff)
      if(radius<=epsilon(1d0)*max(1d0,abs(center))) then
        stop "DC-LCFO CheFSI: invalid filter interval."
      end if

      allocate(q0(layout%lda,ncol),q1(layout%lda,ncol))
      allocate(q2(layout%lda,ncol),work(layout%lda,ncol))
      q0 = x(:,1:ncol)
      call apply_hamiltonian(s,layout,q0,work)
      q1 = (work-center*q0)/radius
      do degree=2,filter_degree
        call apply_hamiltonian(s,layout,q1,work)
        q2 = 2d0*(work-center*q1)/radius-q0
        q0 = q1
        q1 = q2
      end do
      x(:,1:ncol) = q1
      deallocate(work,q2,q1,q0)
    end subroutine chebyshev_filter

    subroutine orthonormalize(layout_n,layout_d,layout_g,x)
      implicit none
      type(s_matrix_layout), intent(in) :: layout_n,layout_d,layout_g
      real(8), intent(inout) :: x(:,:)
      integer :: info,pass
      real(8), allocatable :: xd(:,:),gram(:,:)

      allocate(xd(layout_d%lda,max(1,layout_d%ncol_local)))
      allocate(gram(layout_g%lda,max(1,layout_g%ncol_local)))
      xd = 0d0
      call redistribute_to_dense(layout_n,x,layout_d,xd)
      do pass=1,2
        gram = 0d0
        call pdsyrk('L','T',nsub,npadded,1d0,xd,1,1, &
        & layout_d%desc,0d0,gram,1,1,layout_g%desc)
        call pdpotrf('L',nsub,gram,1,1,layout_g%desc,info)
        if(info/=0) then
          call distributed_qr(layout_d,xd)
          exit
        end if
        call pdtrsm('R','L','T','N',npadded,nsub,1d0,gram, &
        & 1,1,layout_g%desc,xd,1,1,layout_d%desc)
      end do
      call redistribute_to_natural(layout_d,xd,layout_n,x)
      deallocate(gram,xd)
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

    subroutine rayleigh_ritz(s,layout_n,layout_d,layout_g,x,eigenvalue)
      implicit none
      integer, intent(in) :: s
      type(s_matrix_layout), intent(in) :: layout_n,layout_d,layout_g
      real(8), intent(inout) :: x(:,:)
      real(8), intent(out) :: eigenvalue(:)
      real(8), allocatable :: hx(:,:),xd(:,:),hxd(:,:),xrot(:,:)
      real(8), allocatable :: projected(:,:),ritz_vector(:,:)

      allocate(hx(layout_n%lda,max(1,layout_n%ncol_local)))
      allocate(xd(layout_d%lda,max(1,layout_d%ncol_local)))
      allocate(hxd(layout_d%lda,max(1,layout_d%ncol_local)))
      allocate(xrot(layout_d%lda,max(1,layout_d%ncol_local)))
      allocate(projected(layout_g%lda,max(1,layout_g%ncol_local)))
      allocate(ritz_vector(layout_g%lda,max(1,layout_g%ncol_local)))

      call apply_hamiltonian(s,layout_n,x,hx)
      xd = 0d0
      hxd = 0d0
      call redistribute_to_dense(layout_n,x,layout_d,xd)
      call redistribute_to_dense(layout_n,hx,layout_d,hxd)
      projected = 0d0
      call pdgemm('T','N',nsub,nsub,npadded,1d0,xd,1,1, &
      & layout_d%desc,hxd,1,1,layout_d%desc,0d0,projected, &
      & 1,1,layout_g%desc)
      call diagonalize_projected(layout_g,projected,eigenvalue, &
      & ritz_vector)
      xrot = 0d0
      call pdgemm('N','N',npadded,nsub,nsub,1d0,xd,1,1, &
      & layout_d%desc,ritz_vector,1,1,layout_g%desc,0d0, &
      & xrot,1,1,layout_d%desc)
      call redistribute_to_natural(layout_d,xrot,layout_n,x)

      deallocate(ritz_vector,projected,xrot,hxd,xd,hx)
    end subroutine rayleigh_ritz

    subroutine diagonalize_projected(layout,projected,eigenvalue,vector)
      implicit none
      type(s_matrix_layout), intent(in) :: layout
      real(8), intent(inout) :: projected(:,:)
      real(8), intent(out) :: eigenvalue(:),vector(:,:)
      integer :: info,lwork,liwork,lwork_min
      integer, allocatable :: iwork(:)
      real(8) :: work_query(1)
      real(8), allocatable :: work(:)

      liwork = max(1,2+7*nsub+8*grid_dense%npcol)
      allocate(iwork(liwork))
      call pdsyevd('V','L',nsub,projected,1,1,layout%desc, &
      & eigenvalue,vector,1,1,layout%desc,work_query,-1, &
      & iwork,liwork,info)
      if(info/=0) stop "DC-LCFO CheFSI: PDSYEVD workspace query failed."
      lwork_min = max(1+6*nsub+2*layout%nrow_local* &
      & layout%ncol_local,3*nsub+max(dense_block_size* &
      & (layout%nrow_local+1),3*dense_block_size))
      ! PDSYEVD partitions WORK internally before calling PDLASRT.  Some
      ! ScaLAPACK versions under-report that nested workspace in a query.
      lwork = max(10*lwork_min,10*nint(work_query(1)))
      allocate(work(lwork))
      call pdsyevd('V','L',nsub,projected,1,1,layout%desc, &
      & eigenvalue,vector,1,1,layout%desc,work,lwork, &
      & iwork,liwork,info)
      if(info/=0) then
        if(dc%id_tot==0) write(*,*) "PDSYEVD error:",info
        stop "DC-LCFO CheFSI: projected eigensolver failed."
      end if
      deallocate(work,iwork)
    end subroutine diagonalize_projected

    subroutine calculate_residual(s,layout,x,eigenvalue,nstate,max_residual)
      implicit none
      integer, intent(in) :: s,nstate
      type(s_matrix_layout), intent(in) :: layout
      real(8), intent(in) :: x(:,:),eigenvalue(:)
      real(8), intent(out) :: max_residual
      integer :: global_column,local_column,nb
      real(8), allocatable :: hx(:,:),residual_local(:),residual(:)

      allocate(hx(layout%lda,max(1,layout%ncol_local)))
      allocate(residual_local(nstate),residual(nstate))
      call apply_hamiltonian(s,layout,x,hx)
      residual_local = 0d0
      nb = n_basis(dc%i_frag,s)
      do local_column=1,layout%ncol_local
        global_column = natural_global_column(local_column, &
        & grid_natural%mycol,nsub,grid_natural%npcol)
        if(global_column>nstate) cycle
        residual_local(global_column) = sum((hx(1:nb,local_column) &
        & -eigenvalue(global_column)*x(1:nb,local_column))**2)
      end do
      call comm_summation(residual_local,residual,nstate,dc%icomm_tot)
      do global_column=1,nstate
        residual(global_column) = sqrt(max(0d0,residual(global_column))) &
        & /max(1d0,abs(eigenvalue(global_column)))
      end do
      max_residual = maxval(residual)
      deallocate(residual,residual_local,hx)
    end subroutine calculate_residual

    function estimate_upper_bound(s) result(upper)
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
    end function estimate_upper_bound

    subroutine export_coefficients(s,layout,x,nstate)
      implicit none
      integer, intent(in) :: s,nstate
      type(s_matrix_layout), intent(in) :: layout
      real(8), intent(in) :: x(:,:)
      integer :: global_column,local_column,nb
      real(8), allocatable :: local_coef(:,:),fragment_coef(:,:)

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
    end subroutine export_coefficients

    subroutine redistribute_to_dense(layout_n,xn,layout_d,xd)
      implicit none
      type(s_matrix_layout), intent(in) :: layout_n,layout_d
      real(8), intent(in) :: xn(:,:)
      real(8), intent(inout) :: xd(:,:)
      call pdgemr2d(npadded,nsub,xn,1,1,layout_n%desc,xd,1,1, &
      & layout_d%desc,grid_natural%ictxt)
    end subroutine redistribute_to_dense

    subroutine redistribute_to_natural(layout_d,xd,layout_n,xn)
      implicit none
      type(s_matrix_layout), intent(in) :: layout_d,layout_n
      real(8), intent(in) :: xd(:,:)
      real(8), intent(inout) :: xn(:,:)
      call pdgemr2d(npadded,nsub,xd,1,1,layout_d%desc,xn,1,1, &
      & layout_n%desc,grid_natural%ictxt)
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
