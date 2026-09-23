! Adaptively compressed exchange: fixed occupied source state, arbitrary targets.
module hse_ace
  use iso_c_binding, only: c_double,c_double_complex
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: hse_ace_state,hse_ace_build,hse_ace_apply
  type hse_ace_state
    complex(c_double_complex),allocatable :: factors(:,:,:)
    real(c_double) :: dv=0d0,condition=0d0
  end type
contains
  subroutine hse_ace_build(ace,u,w,dv,ierr)
    type(hse_ace_state),intent(inout) :: ace
    complex(c_double_complex),intent(in) :: u(:,:,:),w(:,:,:)
    real(c_double),intent(in) :: dv
    integer,intent(out) :: ierr
    complex(c_double_complex),allocatable :: metric(:,:),work(:)
    real(c_double),allocatable :: e(:),rwork(:)
    real(c_double) :: scale
    integer :: n,ng,nk,ik,j,status
    complex(c_double_complex),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    external :: zgemm,zheev
    ierr=1
    if(allocated(ace%factors))deallocate(ace%factors)
    if(any(shape(u)/=shape(w)).or.dv<=0.or..not.ieee_is_finite(dv))return
    if(.not.all(ieee_is_finite(real(u))).or..not.all(ieee_is_finite(aimag(u))))return
    if(.not.all(ieee_is_finite(real(w))).or..not.all(ieee_is_finite(aimag(w))))return
    ng=size(u,1);n=size(u,2);nk=size(u,3)
    if(min(ng,n,nk)<1)return
    allocate(metric(n,n),work(max(1,2*n)),e(n),rwork(max(1,3*n-2)),ace%factors(ng,n,nk))
    ace%dv=dv;ace%condition=0
    do ik=1,nk
      call zgemm('C','N',n,n,ng,-one*dv,u(1,1,ik),ng,w(1,1,ik),ng,zero,metric(1,1),n)
      scale=sqrt(sum(abs(metric)**2))
      if(scale==0d0.or.sqrt(sum(abs(metric-transpose(conjg(metric)))**2))>1d-10*scale)goto 900
      metric=.5d0*(metric+transpose(conjg(metric)))
      call zheev('V','U',n,metric,n,e,work,size(work),rwork,status)
      if(status/=0)goto 900
      if(e(n)<=0.or.e(1)<=1d-12*e(n))goto 900
      ace%condition=max(ace%condition,e(n)/e(1))
      do j=1,n;metric(:,j)=metric(:,j)/sqrt(e(j));enddo
      call zgemm('N','N',ng,n,n,one,w(1,1,ik),ng,metric(1,1),n,zero,ace%factors(1,1,ik),ng)
    enddo
    ierr=0;return
900 deallocate(ace%factors)
  end subroutine

  subroutine hse_ace_apply(ace,target,action,ierr)
    type(hse_ace_state),intent(in) :: ace
    complex(c_double_complex),intent(in) :: target(:,:,:)
    complex(c_double_complex),intent(out) :: action(:,:,:)
    integer,intent(out) :: ierr
    complex(c_double_complex),allocatable :: overlap(:,:)
    complex(c_double_complex),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    integer :: ng,no,nt,nk,ik
    external :: zgemm
    ierr=1;action=zero
    if(.not.allocated(ace%factors))return
    ng=size(ace%factors,1);no=size(ace%factors,2);nk=size(ace%factors,3);nt=size(target,2)
    if(size(target,1)/=ng.or.size(target,3)/=nk.or.nt<1.or.any(shape(target)/=shape(action)))return
    allocate(overlap(no,nt))
    do ik=1,nk
      call zgemm('C','N',no,nt,ng,one*ace%dv,ace%factors(1,1,ik),ng,target(1,1,ik),ng,zero,overlap(1,1),no)
      call zgemm('N','N',ng,nt,no,-one,ace%factors(1,1,ik),ng,overlap(1,1),no,zero,action(1,1,ik),ng)
    enddo
    ierr=0
  end subroutine
end module
