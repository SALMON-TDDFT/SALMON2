! Column-streamed linear algebra on a Cartesian spatial/orbital process grid.
! Mesh arrays contain local orbital columns; ScaLAPACK ACE metrics use tiles.
module exx_orbitals
  use exx_sparse_orbitals, only: s_sparse_orbitals,sparse_dot,sparse_valid
  use communication, only: comm_get_groupinfo,comm_get_max,comm_summation,comm_bcast
  use exx_ace, only: s_exx_ace,exx_ace_clear
  use exx_distributed_metric, only: distributed_metric_available,distributed_metric_build, &
    distributed_metric_rotate,distributed_metric_apply
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: orbital_layout,orbital_check,orbital_overlap,orbital_rotate
  public :: orbital_ace_build,orbital_ace_apply,orbital_hermitian_action
  ! Rank-specific loops keep IEEE inquiries scalar on Fujitsu compilers.
  private :: salmon_all_finite,finite_real_1d,finite_real_2d,finite_real_3d
  interface salmon_all_finite
    module procedure finite_real_1d,finite_real_2d,finite_real_3d
  end interface
contains
  subroutine orbital_check(bad,comm_r,comm_o)
    implicit none
    integer,intent(inout) :: bad
    integer,intent(in) :: comm_r,comm_o
    call comm_get_max(bad,comm_r)
    call comm_get_max(bad,comm_o)
  end subroutine

  subroutine orbital_layout(nlocal,comm_r,comm_o,counts,first,status)
    implicit none
    integer,intent(in) :: nlocal,comm_r,comm_o
    integer,allocatable,intent(out) :: counts(:)
    integer,intent(out) :: first,status
    integer :: rank,np,owner,largest
    call comm_get_groupinfo(comm_o,rank,np)
    allocate(counts(0:np-1))
    do owner=0,np-1
      counts(owner)=nlocal
      call comm_bcast(counts(owner),comm_o,owner)
    enddo
    status=0
    do owner=0,np-1
      largest=counts(owner)
      call comm_get_max(largest,comm_r)
      if(largest/=counts(owner).or.counts(owner)<0)status=1
    enddo
    call orbital_check(status,comm_r,comm_o)
    first=1+sum(counts(:rank-1))
  end subroutine

  subroutine orbital_overlap(left,right,dv,comm_r,comm_o,counts,first,matrix,phase,sparse_left)
    ! left/right have the same contiguous column ownership. Optional phase is
    ! applied to right grid rows. Global small matrix is identical on all peers.
    implicit none
    complex(8),intent(in) :: left(:,:),right(:,:)
    real(8),intent(in) :: dv
    integer,intent(in) :: comm_r,comm_o,counts(0:),first
    complex(8),intent(out) :: matrix(:,:)
    complex(8),intent(in),optional :: phase(:)
    type(s_sparse_orbitals),intent(in),optional :: sparse_left
    complex(8),allocatable :: column(:),total(:,:)
    integer :: owner,rank,np,j,k,i
    call comm_get_groupinfo(comm_o,rank,np)
    allocate(column(size(right,1)),total(size(matrix,1),size(matrix,2)))
    matrix=0d0;k=0
    do owner=0,np-1
      do j=1,counts(owner)
        k=k+1
        if(rank==owner)column=right(:,j)
        call comm_bcast(column,comm_o,owner)
        if(present(phase))column=column*phase
        if(present(sparse_left))then
          do i=1,size(left,2)
            matrix(first+i-1,k)=sparse_dot(sparse_left,i,column)*dv
          enddo
        else
          matrix(first:first+size(left,2)-1,k)=matmul(conjg(transpose(left)),column)*dv
        endif
      enddo
    enddo
    call comm_summation(matrix,total,size(matrix),comm_r)
    call comm_summation(total,matrix,size(matrix),comm_o)
  end subroutine

  subroutine orbital_rotate(input,matrix,comm_o,counts,first,output,weights)
    implicit none
    complex(8),intent(in) :: input(:,:),matrix(:,:)
    integer,intent(in) :: comm_o,counts(0:),first
    complex(8),intent(out) :: output(:,:)
    real(8),intent(in),optional :: weights(:)
    complex(8),allocatable :: column(:)
    integer :: rank,np,owner,i,j,k,g
    call comm_get_groupinfo(comm_o,rank,np)
    allocate(column(size(input,1)))
    output=0d0;k=0
    do owner=0,np-1
      do i=1,counts(owner)
        k=k+1
        if(rank==owner)then
          column=input(:,i)
          if(present(weights))column=column*weights(i)
        endif
        call comm_bcast(column,comm_o,owner)
        ! Communication stays outside OMP; each worker owns distinct output elements.
!$omp parallel do collapse(2) default(none) schedule(static) &
!$omp shared(output,column,matrix,k,first) private(j,g)
        do j=1,size(output,2)
          do g=1,size(output,1)
            output(g,j)=output(g,j)+column(g)*matrix(k,first+j-1)
          enddo
        enddo
!$omp end parallel do
      enddo
    enddo
  end subroutine

  subroutine orbital_hermitian_action(u,w,dv,comm_r,comm_o,budget,correction_norm,status)
    ! Complete the projected metric without assuming perfectly orthonormal U.
    ! A=U^dagger W, G=U^dagger U, D=U G^-1 (A^dagger-A)/2.
    ! Measure ||D|| on distributed rows/columns before accepting W+D.
    implicit none
    complex(8),intent(in) :: u(:,:,:)
    complex(8),intent(inout) :: w(:,:,:)
    real(8),intent(in) :: dv,budget
    integer,intent(in) :: comm_r,comm_o
    real(8),intent(out) :: correction_norm
    integer,intent(out) :: status
    integer,allocatable :: counts(:)
    complex(8),allocatable :: gram(:,:),metric(:,:),delta(:,:)
    real(8) :: local_norm,spatial_norm,total_norm,scale,scale_max(1),scale_global(1)
    integer :: first,n,bad
    bad=0;status=1;correction_norm=0d0
    if(any(shape(u)/=shape(w)).or.size(u,3)/=1.or.dv<=0d0.or..not.ieee_is_finite(dv))bad=1
    if(.not.ieee_is_finite(budget).or.budget<0d0)bad=1
    if(.not.salmon_all_finite(real(u)).or..not.salmon_all_finite(aimag(u)))bad=1
    if(.not.salmon_all_finite(real(w)).or..not.salmon_all_finite(aimag(w)))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    call orbital_layout(size(u,2),comm_r,comm_o,counts,first,bad)
    if(bad/=0.or.sum(counts)<1)return
    n=sum(counts);allocate(gram(n,n),metric(n,n),delta(size(u,1),size(u,2)))
    call orbital_overlap(u(:,:,1),u(:,:,1),dv,comm_r,comm_o,counts,first,gram)
    call orbital_overlap(u(:,:,1),w(:,:,1),dv,comm_r,comm_o,counts,first,metric)
    gram=.5d0*(gram+conjg(transpose(gram)))
    metric=.5d0*(conjg(transpose(metric))-metric)
    if(.not.salmon_all_finite(real(gram)).or..not.salmon_all_finite(aimag(gram)))bad=1
    if(.not.salmon_all_finite(real(metric)).or..not.salmon_all_finite(aimag(metric)))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    call zposv('U',n,n,gram,n,metric,n,bad)
    if(bad/=0)bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    call orbital_rotate(u(:,:,1),metric,comm_o,counts,first,delta)
    scale=0d0
    if(size(delta)>0)scale=maxval(abs(delta))
    call comm_get_max([scale],scale_max,1,comm_r)
    call comm_get_max(scale_max,scale_global,1,comm_o)
    scale=scale_global(1)
    local_norm=0d0
    if(scale>0d0)local_norm=sum((abs(delta)/scale)**2)*dv
    call comm_summation(local_norm,spatial_norm,comm_r)
    call comm_summation(spatial_norm,total_norm,comm_o)
    correction_norm=scale*sqrt(total_norm)*(1d0+128d0*epsilon(1d0))
    if(.not.ieee_is_finite(correction_norm).or.correction_norm>budget)bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    w(:,:,1)=w(:,:,1)+delta
    status=0
  end subroutine

  subroutine orbital_ace_build(ace,u,w,dv,comm_r,comm_o,status,packed,comm_matrix,sparse_u)
    use exx_blas_threads, only: exx_blas_thread_control
!$  use omp_lib, only: omp_get_max_threads,omp_in_parallel
    implicit none
    type(s_exx_ace),intent(inout) :: ace
    complex(8),intent(in) :: u(:,:,:),w(:,:,:)
    real(8),intent(in) :: dv
    integer,intent(in) :: comm_r,comm_o
    integer,intent(out) :: status
    integer,intent(in),optional :: comm_matrix
    logical,intent(in),optional :: packed
    type(s_sparse_orbitals),intent(in),optional :: sparse_u
    integer :: workers,thread_state(3)
    workers=1;thread_state=0
!$  workers=omp_get_max_threads()
!$  if(omp_in_parallel())workers=1
    if(workers>1)call exx_blas_thread_control(workers,thread_state)
    ! The core has early exits and an internal packer; keep cleanup in this wrapper.
    call orbital_ace_build_core(ace,u,w,dv,comm_r,comm_o,status,packed,comm_matrix,sparse_u)
    call exx_blas_thread_control(0,thread_state)
  end subroutine orbital_ace_build

  subroutine orbital_ace_build_core(ace,u,w,dv,comm_r,comm_o,status,packed,comm_matrix,sparse_u)
    implicit none
    type(s_exx_ace),intent(inout) :: ace
    complex(8),intent(in) :: u(:,:,:),w(:,:,:)
    real(8),intent(in) :: dv
    integer,intent(in) :: comm_r,comm_o
    integer,intent(out) :: status
    integer,intent(in),optional :: comm_matrix
    logical,intent(in),optional :: packed
    type(s_sparse_orbitals),intent(in),optional :: sparse_u
    logical :: store_packed
    integer,allocatable :: counts(:)
    complex(8),allocatable :: metric(:,:),work(:)
    real(8),allocatable :: e(:),rwork(:)
    integer :: first,n,j,bad,nonzero
    real(8) :: scale
    status=1;bad=0
    call exx_ace_clear(ace)
    store_packed=.false.
    if(present(packed))store_packed=packed
    if(any(shape(u)/=shape(w)).or.size(u,3)/=1.or.dv<=0d0.or..not.ieee_is_finite(dv))bad=1
    ! With sparse_u, u supplies only the dimensions (callers may pass w twice).
    if(present(sparse_u))then
      if(.not.sparse_valid(sparse_u,size(w,1),size(w,2)))bad=1
    else
      if(.not.salmon_all_finite(real(u)).or..not.salmon_all_finite(aimag(u)))bad=1
    endif
    if(.not.salmon_all_finite(real(w)).or..not.salmon_all_finite(aimag(w)))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    call orbital_layout(size(u,2),comm_r,comm_o,counts,first,bad)
    if(bad/=0.or.sum(counts)<1)return
    n=sum(counts)
    nonzero=0
    if(any(w/=(0d0,0d0)))nonzero=1
    call orbital_check(nonzero,comm_r,comm_o)
    ace%dv=dv;ace%condition=0d0
    if(present(comm_matrix))then
      if(distributed_metric_available(comm_matrix))then
        call distributed_metric_build(ace,u(:,:,1),w(:,:,1),dv,comm_o,comm_matrix,counts,first,nonzero,bad,sparse_u)
        if(bad/=0)then
          call exx_ace_clear(ace)
          return
        endif
        if(store_packed)then
          call pack_action()
        else
          allocate(ace%factors(size(u,1),size(u,2),1))
          call distributed_metric_rotate(ace,w(:,:,1),comm_o,counts,first,ace%factors(:,:,1))
          deallocate(ace%metric_factor,ace%metric_rows,ace%metric_cols)
          ace%metric_distributed=.false.
        endif
        status=0
        return
      endif
    endif
    allocate(metric(n,n),work(2*n),e(n),rwork(max(1,3*n-2)))
    if(nonzero==0)then
      if(store_packed)then
        allocate(ace%metric_factor(n,n));ace%metric_factor=0d0
        call pack_action()
      else
        allocate(ace%factors(size(u,1),size(u,2),1));ace%factors=0d0
      endif
      status=0;return
    endif
    call orbital_overlap(u(:,:,1),w(:,:,1),-dv,comm_r,comm_o,counts,first,metric,sparse_left=sparse_u)
    if(.not.salmon_all_finite(real(metric)).or..not.salmon_all_finite(aimag(metric)))bad=1
    scale=sqrt(sum(abs(metric)**2))
    if(scale==0d0.or.sqrt(sum(abs(metric-transpose(conjg(metric)))**2))>1d-10*scale)bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    metric=.5d0*(metric+transpose(conjg(metric)))
    call zheev('V','U',n,metric,n,e,work,size(work),rwork,bad)
    if(bad/=0)bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    if(e(n)<=0d0.or.e(1)<=1d-12*e(n))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    ace%condition=e(n)/e(1)
    do j=1,n
      metric(:,j)=metric(:,j)/sqrt(e(j))
    enddo
    if(store_packed)then
      ace%metric_factor=metric
      call pack_action()
    else
      allocate(ace%factors(size(u,1),size(u,2),1))
      call orbital_rotate(w(:,:,1),metric,comm_o,counts,first,ace%factors(:,:,1))
    endif
    status=0
  contains
    subroutine pack_action()
      implicit none
      integer :: column,g,k
      allocate(ace%offset(size(w,2)+1),ace%row(count(w/=(0d0,0d0))),ace%values(count(w/=(0d0,0d0))))
      k=1
      do column=1,size(w,2)
        ace%offset(column)=k
        do g=1,size(w,1)
          if(w(g,column,1)==(0d0,0d0))cycle
          ace%row(k)=g;ace%values(k)=w(g,column,1);k=k+1
        enddo
      enddo
      ace%offset(size(w,2)+1)=k;ace%grid_rows=size(w,1);ace%packed=.true.
    end subroutine
  end subroutine

  subroutine orbital_ace_apply(ace,target,action,comm_r,comm_o,status)
    implicit none
    type(s_exx_ace),intent(in) :: ace
    complex(8),intent(in) :: target(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    integer,intent(in) :: comm_r,comm_o
    integer,intent(out) :: status
    integer,allocatable :: counts(:),target_counts(:)
    complex(8),allocatable :: column(:),overlap(:),total(:)
    integer :: first,rank,np,owner,i,j,g,bad,packed_mode
    complex(8) :: dot
    status=1;action=0d0;bad=0
    packed_mode=merge(1,0,ace%packed)
    call comm_get_max(packed_mode,comm_r);call comm_get_max(packed_mode,comm_o)
    if((packed_mode==1).neqv.ace%packed)bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    if(ace%packed)then
      call packed_ace_apply(ace,target,action,comm_r,comm_o,status)
      return
    endif
    if(.not.allocated(ace%factors))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    if(size(ace%factors,1)/=size(target,1).or.size(ace%factors,3)/=1.or.size(target,3)/=1)bad=1
    if(any(shape(target)/=shape(action)).or.ace%dv<=0d0.or..not.ieee_is_finite(ace%dv))bad=1
    if(.not.salmon_all_finite(real(target)).or..not.salmon_all_finite(aimag(target)))bad=1
    if(.not.salmon_all_finite(real(ace%factors)).or..not.salmon_all_finite(aimag(ace%factors)))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    call orbital_layout(size(target,2),comm_r,comm_o,target_counts,first,bad)
    if(bad/=0)return
    call orbital_layout(size(ace%factors,2),comm_r,comm_o,counts,first,bad)
    if(bad/=0)return
    allocate(column(size(target,1)),overlap(size(target,2)),total(size(target,2)))
    call comm_get_groupinfo(comm_o,rank,np)
    do owner=0,np-1
      do i=1,counts(owner)
        if(rank==owner)column=ace%factors(:,i,1)
        call comm_bcast(column,comm_o,owner)
        ! Keep grid summation order within each independent target orbital.
!$omp parallel do default(none) private(j,g,dot) shared(target,column,overlap,ace)
        do j=1,size(target,2)
          dot=0d0
          do g=1,size(target,1)
            dot=dot+conjg(column(g))*target(g,j,1)
          enddo
          overlap(j)=dot*ace%dv
        enddo
!$omp end parallel do
        call comm_summation(overlap,total,size(overlap),comm_r)
!$omp parallel do default(none) private(j) shared(target,action,column,total)
        do j=1,size(target,2)
          action(:,j,1)=action(:,j,1)-column*total(j)
        enddo
!$omp end parallel do
      enddo
    enddo
    status=0
  end subroutine
  subroutine packed_ace_apply(ace,target,action,comm_r,comm_o,status)
    implicit none
    type(s_exx_ace),intent(in) :: ace
    complex(8),intent(in) :: target(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    integer,intent(in) :: comm_r,comm_o
    integer,intent(out) :: status
    integer,allocatable :: counts(:),target_counts(:)
    complex(8),allocatable :: column(:),partial_action(:),total_action(:),overlap(:),overlap_o(:),overlap_r(:),coeff(:),rotated(:)
    integer :: bad,ng,n,first,target_first,rank,np,owner,i,j,k,start,finish,distributed_mode
    status=1;action=0d0;bad=0;ng=size(target,1)
    distributed_mode=merge(1,0,ace%metric_distributed)
    call orbital_check(distributed_mode,comm_r,comm_o)
    if((distributed_mode==1).neqv.ace%metric_distributed)bad=1
    if(.not.allocated(ace%offset).or..not.allocated(ace%row).or. &
       .not.allocated(ace%values).or..not.allocated(ace%metric_factor))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    if(ng/=ace%grid_rows.or.size(target,3)/=1.or.any(shape(target)/=shape(action)))bad=1
    if(ace%dv<=0d0.or..not.ieee_is_finite(ace%dv))bad=1
    if(size(ace%offset)<1.or.size(ace%row)/=size(ace%values))bad=1
    if(.not.salmon_all_finite(real(target)).or..not.salmon_all_finite(aimag(target)))bad=1
    if(.not.salmon_all_finite(real(ace%values)).or..not.salmon_all_finite(aimag(ace%values)))bad=1
    if(.not.salmon_all_finite(real(ace%metric_factor)).or..not.salmon_all_finite(aimag(ace%metric_factor)))bad=1
    if(any(ace%row<1).or.any(ace%row>ng))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    if(ace%offset(1)/=1.or.ace%offset(size(ace%offset))/=size(ace%values)+1)bad=1
    if(any(ace%offset(2:)<ace%offset(:size(ace%offset)-1)))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    call orbital_layout(size(ace%offset)-1,comm_r,comm_o,counts,first,bad)
    if(bad/=0)return
    n=sum(counts)
    if(n<1)bad=1
    if(ace%metric_distributed)then
      if(.not.allocated(ace%metric_rows).or..not.allocated(ace%metric_cols))bad=1
      call orbital_check(bad,comm_r,comm_o)
      if(bad/=0)return
      if(ace%metric_order/=n)bad=1
      if(any(shape(ace%metric_factor)/=[size(ace%metric_rows),size(ace%metric_cols)]))bad=1
      if(any(ace%metric_rows<0).or.any(ace%metric_rows>n))bad=1
      if(any(ace%metric_cols<0).or.any(ace%metric_cols>n))bad=1
    else
      if(any(shape(ace%metric_factor)/=[n,n]))bad=1
    endif
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    call orbital_layout(size(target,2),comm_r,comm_o,target_counts,target_first,bad)
    if(bad/=0)return
    call comm_get_groupinfo(comm_o,rank,np)
    allocate(column(ng),partial_action(ng),total_action(ng),overlap(n),overlap_o(n),overlap_r(n),coeff(n),rotated(n))
    ! Stream each target; only a column and band vectors are communicated.
    do owner=0,np-1
      do i=1,target_counts(owner)
        if(rank==owner)column=target(:,i,1)
        call comm_bcast(column,comm_o,owner)
        overlap=0d0
!$omp parallel do default(none) private(j,start,finish) shared(ace,column,overlap,first)
        do j=1,size(ace%offset)-1
          start=ace%offset(j);finish=ace%offset(j+1)-1
          overlap(first+j-1)=sum(conjg(ace%values(start:finish))*column(ace%row(start:finish)))*ace%dv
        enddo
!$omp end parallel do
        call comm_summation(overlap,overlap_o,n,comm_o)
        call comm_summation(overlap_o,overlap_r,n,comm_r)
        ! Do not form A A^H explicitly: near the accepted conditioning limit,
        ! that inverse loses cancellation accuracy relative to dense ACE factors.
        if(ace%metric_distributed)then
          call distributed_metric_apply(ace,overlap_r,rotated,coeff,overlap)
        else
          call zgemv('C',n,n,(1d0,0d0),ace%metric_factor,n,overlap_r,1,(0d0,0d0),rotated,1)
          call zgemv('N',n,n,(1d0,0d0),ace%metric_factor,n,rotated,1,(0d0,0d0),coeff,1)
        endif
        partial_action=0d0
        do j=1,size(ace%offset)-1
          do k=ace%offset(j),ace%offset(j+1)-1
            partial_action(ace%row(k))=partial_action(ace%row(k))-ace%values(k)*coeff(first+j-1)
          enddo
        enddo
        call comm_summation(partial_action,total_action,ng,comm_o)
        if(rank==owner)action(:,i,1)=total_action
      enddo
    enddo
    status=0
  end subroutine


  pure logical function finite_real_1d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:)
    real(8) :: value
    integer :: i
    finite=.false.
    do i=1,size(values,1)
      value=values(i)
      if(.not.ieee_is_finite(value))return
    enddo
    finite=.true.
  end function

  pure logical function finite_real_2d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:,:)
    real(8) :: value
    integer :: i,j
    finite=.false.
    do j=1,size(values,2)
      do i=1,size(values,1)
        value=values(i,j)
        if(.not.ieee_is_finite(value))return
      enddo
    enddo
    finite=.true.
  end function

  pure logical function finite_real_3d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:,:,:)
    real(8) :: value
    integer :: i,j,k
    finite=.false.
    do k=1,size(values,3)
      do j=1,size(values,2)
        do i=1,size(values,1)
          value=values(i,j,k)
          if(.not.ieee_is_finite(value))return
        enddo
      enddo
    enddo
    finite=.true.
  end function
end module
