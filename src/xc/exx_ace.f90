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

  subroutine exx_ace_build(ace,u,w,dv,ierr,sum_grid)
    implicit none
    type(s_exx_ace),intent(inout) :: ace
    complex(c_double_complex),intent(in) :: u(:,:,:),w(:,:,:)
    real(c_double),intent(in) :: dv
    integer,intent(out) :: ierr
    procedure(grid_sum),optional :: sum_grid
    complex(c_double_complex) :: check(1,1)
    complex(c_double_complex),allocatable :: metric(:,:),work(:)
    real(c_double),allocatable :: e(:),rwork(:)
    real(c_double) :: scale
    integer :: n,ng,nk,ik,j,status
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
    allocate(metric(n,n),work(max(1,2*n)),e(n),rwork(max(1,3*n-2)),ace%factors(ng,n,nk))
    ace%dv=dv;ace%condition=0
    do ik=1,nk
      ! A fragment without electrons has exactly zero exchange, a valid operator.
      check=0d0
      if(any(w(:,:,ik)/=zero))check=1d0
      if(present(sum_grid))call sum_grid(check)
      if(real(check(1,1))==0d0)then
        ace%factors(:,:,ik)=zero
        cycle
      endif
      metric=zero
      if(ng>0)call zgemm('C','N',n,n,ng,-one*dv,u(1,1,ik),ng,w(1,1,ik),ng,zero,metric(1,1),n)
      if(present(sum_grid))call sum_grid(metric)
      if(.not.finite_matrix(metric))goto 900
      scale=sqrt(sum(abs(metric)**2))
      if(scale==0d0.or.sqrt(sum(abs(metric-transpose(conjg(metric)))**2))>1d-10*scale)goto 900
      metric=.5d0*(metric+transpose(conjg(metric)))
      call zheev('V','U',n,metric,n,e,work,size(work),rwork,status)
      if(status/=0)goto 900
      if(e(n)<=0.or.e(1)<=1d-12*e(n))goto 900
      ace%condition=max(ace%condition,e(n)/e(1))
      do j=1,n;metric(:,j)=metric(:,j)/sqrt(e(j));enddo
      if(ng>0)call zgemm('N','N',ng,n,n,one,w(1,1,ik),ng,metric(1,1),n,zero,ace%factors(1,1,ik),ng)
    enddo
    ierr=0;return
900 deallocate(ace%factors)
  end subroutine

  subroutine exx_ace_apply(ace,target,action,ierr,sum_grid)
    implicit none
    type(s_exx_ace),intent(in) :: ace
    complex(c_double_complex),intent(in) :: target(:,:,:)
    complex(c_double_complex),intent(out) :: action(:,:,:)
    integer,intent(out) :: ierr
    procedure(grid_sum),optional :: sum_grid
    complex(c_double_complex) :: check(1,1)
    complex(c_double_complex),allocatable :: overlap(:,:)
    complex(c_double_complex),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    integer :: ng,no,nt,nk,ik
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
    allocate(overlap(no,nt))
    do ik=1,nk
      overlap=zero
      if(ng>0)call zgemm('C','N',no,nt,ng,one*ace%dv,ace%factors(1,1,ik),ng,target(1,1,ik),ng,zero,overlap(1,1),no)
      if(present(sum_grid))call sum_grid(overlap)
      if(ng>0)call zgemm('N','N',ng,nt,no,-one,ace%factors(1,1,ik),ng,overlap(1,1),no,zero,action(1,1,ik),ng)
    enddo
    deallocate(overlap)
    ierr=0
  end subroutine
end module
