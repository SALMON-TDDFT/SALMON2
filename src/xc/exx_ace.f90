! Adaptively compressed exchange: fixed occupied source state, arbitrary targets.
module exx_ace
  use iso_c_binding, only: c_double,c_double_complex
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: s_exx_ace,exx_ace_build,exx_ace_apply,exx_ace_average,exx_ace_clear,exx_ace_ready
  ! Spatial peers must share orbital/k dimensions, dv and collective call order.
  ! The callback sums a small matrix in place; grid rows are never gathered.
  abstract interface
    subroutine blas_thread_control(requested,state)
      implicit none
      integer,intent(in) :: requested
      integer,intent(inout) :: state(3)
    end subroutine
    subroutine grid_sum(matrix)
      import c_double_complex
      implicit none
      complex(c_double_complex),intent(inout) :: matrix(:,:)
    end subroutine
  end interface
  type s_exx_ace
    complex(c_double_complex),allocatable :: factors(:,:,:)
    real(c_double) :: dv=0d0,condition=0d0
    ! Exact-nonzero W and A=V/sqrt(e) with A A^H=(-U^H W dv)^-1.
    logical :: packed=.false.
    ! No BLACS context is retained: ordinary assignment safely copies this state.
    logical :: metric_distributed=.false.
    integer :: metric_comm=0,metric_order=0
    integer,allocatable :: metric_rows(:),metric_cols(:)
    integer :: grid_rows=0
    integer,allocatable :: offset(:),row(:)
    complex(c_double_complex),allocatable :: values(:),metric_factor(:,:)
  end type
contains
  ! Scalar temporaries avoid the array IEEE-expression crash seen with frtpx.
  logical function finite_orbitals(values) result(finite)
    implicit none
    complex(c_double_complex),intent(in) :: values(:,:,:)
    real(c_double) :: component
    integer :: i,j,k
    finite=.false.
    do k=1,size(values,3)
    do j=1,size(values,2)
      do i=1,size(values,1)
        component=real(values(i,j,k),c_double)
        if(.not.ieee_is_finite(component))return
        component=aimag(values(i,j,k))
        if(.not.ieee_is_finite(component))return
      enddo
    enddo
    enddo
    finite=.true.
  end function

  ! Scalar temporaries avoid the array IEEE-expression crash seen with frtpx.
  logical function finite_matrix(values) result(finite)
    implicit none
    complex(c_double_complex),intent(in) :: values(:,:)
    real(c_double) :: component
    integer :: i,j
    finite=.false.
    do j=1,size(values,2)
      do i=1,size(values,1)
        component=real(values(i,j),c_double)
        if(.not.ieee_is_finite(component))return
        component=aimag(values(i,j))
        if(.not.ieee_is_finite(component))return
      enddo
    enddo
    finite=.true.
  end function

  subroutine exx_ace_clear(ace)
    implicit none
    type(s_exx_ace),intent(inout) :: ace
    type(s_exx_ace) :: empty
    ace=empty
  end subroutine

  logical function exx_ace_ready(ace)
    implicit none
    type(s_exx_ace),intent(in) :: ace
    exx_ace_ready=allocated(ace%factors)
    if(ace%packed)exx_ace_ready=allocated(ace%values).and.allocated(ace%metric_factor).and. &
      allocated(ace%offset).and.allocated(ace%row)
  end function

  subroutine exx_ace_average(left,right,average,ierr)
    implicit none
    type(s_exx_ace),intent(in) :: left,right
    type(s_exx_ace),intent(out) :: average
    integer,intent(out) :: ierr
    integer :: ng,nl,nr,nk
    ierr=1
    if(.not.allocated(left%factors).or..not.allocated(right%factors))return
    ng=size(left%factors,1);nl=size(left%factors,2);nk=size(left%factors,3);nr=size(right%factors,2)
    if(size(right%factors,1)/=ng.or.size(right%factors,3)/=nk.or.left%dv/=right%dv)return
    allocate(average%factors(ng,nl+nr,nk))
    average%factors(:,:nl,:)=left%factors/sqrt(2d0)
    average%factors(:,nl+1:,:)=right%factors/sqrt(2d0)
    average%dv=left%dv
    ! This representation applies (K_left+K_right)/2, not an average of orbitals.
    average%condition=max(left%condition,right%condition)
    ierr=0
  end subroutine

  subroutine exx_ace_build(ace,u,w,dv,ierr,sum_grid,thread_control)
!$  use omp_lib, only: omp_get_max_threads,omp_in_parallel
    implicit none
    type(s_exx_ace),intent(inout) :: ace
    complex(c_double_complex),intent(in) :: u(:,:,:),w(:,:,:)
    real(c_double),intent(in) :: dv
    integer,intent(out) :: ierr
    procedure(grid_sum),optional :: sum_grid
    procedure(blas_thread_control),optional :: thread_control
    complex(c_double_complex) :: check(1,1)
    complex(c_double_complex),allocatable :: metric(:,:),work(:)
    real(c_double),allocatable :: e(:),rwork(:)
    integer :: n,ng,nk,ik,status,workers,bad,thread_state(3),worker_state(3)
    real(c_double) :: condition,one_condition
    logical :: parallel_k
    complex(c_double_complex),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    external :: zgemm,zheev
    ierr=1
    call exx_ace_clear(ace)
    check=0d0
    if(any(shape(u)/=shape(w)).or.dv<=0.or..not.ieee_is_finite(dv))check=1d0
    if(.not.finite_orbitals(u))check=1d0
    if(.not.finite_orbitals(w))check=1d0
    ng=size(u,1);n=size(u,2);nk=size(u,3)
    if(min(n,nk)<1)check=1d0
    if(ng<1.and..not.present(sum_grid))check=1d0
    if(present(sum_grid))call sum_grid(check)
    if(real(check(1,1))/=0d0)return
    workers=1;thread_state=0
!$  workers=omp_get_max_threads()
!$  if(omp_in_parallel())workers=1
    parallel_k=nk>1.and..not.present(sum_grid).and.workers>1
    if(present(thread_control).and.workers>1)then
      if(parallel_k)then
        call thread_control(1,thread_state)
      else
        call thread_control(workers,thread_state)
      endif
    endif
    parallel_k=parallel_k.and.thread_state(1)>0
    allocate(ace%factors(ng,n,nk))
    ace%dv=dv;ace%condition=0d0
    bad=0;condition=0d0
    if(parallel_k)then
      workers=min(workers,nk)
!$omp parallel default(none) num_threads(workers) &
!$omp shared(n,nk,ng,u,w,dv,ace) private(metric,work,e,rwork,ik,status,one_condition,worker_state) &
!$omp reduction(max:bad,condition)
      worker_state=0
      call thread_control(1,worker_state)
      allocate(metric(n,n),work(max(1,2*n)),e(n),rwork(max(1,3*n-2)))
!$omp do schedule(static)
      do ik=1,nk
        call build_one(ik,metric,work,e,rwork,one_condition,status)
        condition=max(condition,one_condition)
        bad=max(bad,status)
      enddo
!$omp end do
      deallocate(metric,work,e,rwork)
      call thread_control(0,worker_state)
!$omp end parallel
    else
      allocate(metric(n,n),work(max(1,2*n)),e(n),rwork(max(1,3*n-2)))
      do ik=1,nk
        call build_one(ik,metric,work,e,rwork,one_condition,status,sum_grid)
        condition=max(condition,one_condition)
        bad=max(bad,status)
        if(bad/=0)exit
      enddo
    endif
    if(present(thread_control))call thread_control(0,thread_state)
    if(bad/=0)then
      deallocate(ace%factors)
      return
    endif
    ace%condition=condition;ierr=0
  contains
    subroutine build_one(k,metric,work,e,rwork,kcondition,error,reduce_grid)
      implicit none
      integer,intent(in) :: k
      complex(c_double_complex),intent(inout) :: metric(:,:),work(:)
      real(c_double),intent(inout) :: e(:),rwork(:)
      real(c_double),intent(out) :: kcondition
      integer,intent(out) :: error
      procedure(grid_sum),optional :: reduce_grid
      complex(c_double_complex) :: nonzero(1,1)
      real(c_double) :: norm
      integer :: column,lapack_status
      error=1;nonzero=0d0;kcondition=0d0
      if(any(w(:,:,k)/=zero))nonzero=1d0
      if(present(reduce_grid))call reduce_grid(nonzero)
      if(real(nonzero(1,1))==0d0)then
        ace%factors(:,:,k)=zero;error=0
        return
      endif
      metric=zero
      if(ng>0)call zgemm('C','N',n,n,ng,-one*dv,u(1,1,k),ng,w(1,1,k),ng,zero,metric(1,1),n)
      if(present(reduce_grid))call reduce_grid(metric)
      if(.not.finite_matrix(metric))return
      norm=sqrt(sum(abs(metric)**2))
      if(norm==0d0.or.sqrt(sum(abs(metric-transpose(conjg(metric)))**2))>1d-10*norm)return
      metric=.5d0*(metric+transpose(conjg(metric)))
      call zheev('V','U',n,metric,n,e,work,size(work),rwork,lapack_status)
      if(lapack_status/=0)return
      if(e(n)<=0.or.e(1)<=1d-12*e(n))return
      kcondition=e(n)/e(1)
      do column=1,n
        metric(:,column)=metric(:,column)/sqrt(e(column))
      enddo
      if(ng>0)call zgemm('N','N',ng,n,n,one,w(1,1,k),ng,metric(1,1),n,zero,ace%factors(1,1,k),ng)
      error=0
    end subroutine
  end subroutine

  subroutine exx_ace_apply(ace,target,action,ierr,sum_grid,thread_control)
!$  use omp_lib, only: omp_get_max_threads,omp_in_parallel
    implicit none
    type(s_exx_ace),intent(in) :: ace
    complex(c_double_complex),intent(in) :: target(:,:,:)
    complex(c_double_complex),intent(out) :: action(:,:,:)
    integer,intent(out) :: ierr
    procedure(grid_sum),optional :: sum_grid
    procedure(blas_thread_control),optional :: thread_control
    complex(c_double_complex) :: check(1,1)
    complex(c_double_complex),allocatable :: overlap(:,:)
    complex(c_double_complex),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    integer :: ng,no,nt,nk,ik,workers,thread_state(3),worker_state(3)
    logical :: parallel_k
    external :: zgemm
    ierr=1;action=zero
    check=0d0
    if(.not.allocated(ace%factors))check=1d0
    if(present(sum_grid))call sum_grid(check)
    if(real(check(1,1))/=0d0)return
    ng=size(ace%factors,1);no=size(ace%factors,2);nk=size(ace%factors,3);nt=size(target,2)
    check=0d0
    if(size(target,1)/=ng.or.size(target,3)/=nk.or.nt<1.or.any(shape(target)/=shape(action)))check=1d0
    if(.not.finite_orbitals(target))check=1d0
    if(.not.finite_orbitals(ace%factors))check=1d0
    if(present(sum_grid))call sum_grid(check)
    if(real(check(1,1))/=0d0)return
    workers=1;thread_state=0
!$  workers=omp_get_max_threads()
!$  if(omp_in_parallel())workers=1
    parallel_k=nk>1.and..not.present(sum_grid).and.workers>1
    if(present(thread_control).and.workers>1)then
      if(parallel_k)then
        call thread_control(1,thread_state)
      else
        call thread_control(workers,thread_state)
      endif
    endif
    parallel_k=parallel_k.and.thread_state(1)>0
    if(parallel_k)then
      workers=min(workers,nk)
!$omp parallel default(none) num_threads(workers) &
!$omp shared(no,nt,nk,ng,ace,target,action) private(overlap,ik,worker_state)
      worker_state=0
      call thread_control(1,worker_state)
      allocate(overlap(no,nt))
!$omp do schedule(static)
      do ik=1,nk
        call apply_one(ik,overlap)
      enddo
!$omp end do
      deallocate(overlap)
      call thread_control(0,worker_state)
!$omp end parallel
    else
      allocate(overlap(no,nt))
      do ik=1,nk
        call apply_one(ik,overlap,sum_grid)
      enddo
      deallocate(overlap)
    endif
    if(present(thread_control))call thread_control(0,thread_state)
    ierr=0
  contains
    subroutine apply_one(k,overlap,reduce_grid)
      implicit none
      integer,intent(in) :: k
      complex(c_double_complex),intent(inout) :: overlap(:,:)
      procedure(grid_sum),optional :: reduce_grid
      overlap=zero
      if(ng>0)call zgemm('C','N',no,nt,ng,one*ace%dv,ace%factors(1,1,k),ng,target(1,1,k),ng,zero,overlap(1,1),no)
      if(present(reduce_grid))call reduce_grid(overlap)
      if(ng>0)call zgemm('N','N',ng,nt,no,-one,ace%factors(1,1,k),ng,overlap(1,1),no,zero,action(1,1,k),ng)
    end subroutine
  end subroutine
end module
