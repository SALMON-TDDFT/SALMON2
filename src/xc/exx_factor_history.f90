! Experimental local-row ACE factor prediction: three accepted endpoints only.
module exx_factor_history
 use iso_fortran_env,only:real64,int64
 use,intrinsic::ieee_arithmetic
 implicit none
 private
 public::s_factor_history,history_accept,history_predict,factor_align
 type s_factor_history
  complex(real64),allocatable::x(:,:,:)
  integer::interval=8
  real(real64)::g(2,2)=0,b(2)=0,coeff(2)=[1d0,0d0]
  integer::count=0,last=-1,teachers=0
  logical::ready=.false.
 end type
 abstract interface
  subroutine sum_callback(a)
   import real64
   complex(real64),intent(inout)::a(:,:)
  end subroutine
 end interface
contains

 logical function scalar_finite_real(a) result(ok)
  real(8),intent(in)::a(:)
  integer::i
  ok=.false.
  do i=1,size(a)
   if(.not.ieee_is_finite(a(i)))return
  enddo
  ok=.true.
 end function

 logical function scalar_finite_complex(a) result(ok)
  complex(8),intent(in)::a(:,:)
  integer::i,j
  ok=.false.
  do j=1,size(a,2)
   do i=1,size(a,1)
    if(.not.ieee_is_finite(real(a(i,j),8)))return
    if(.not.ieee_is_finite(aimag(a(i,j))))return
   enddo
  enddo
  ok=.true.
 end function

 subroutine factor_align(x,ref,sumgrid,status)
  implicit none
  complex(real64),intent(inout)::x(:,:)
  complex(real64),intent(in)::ref(:,:)
  procedure(sum_callback)::sumgrid
  integer,intent(out)::status
  complex(real64),allocatable::m(:,:),u(:,:),vh(:,:),r(:,:),work(:),row(:)
  real(real64),allocatable::s(:),rw(:)
  complex(real64)::v,bad(1,1)
  integer::n,ng,i,j,k,info
  status=1;ng=size(x,1);n=size(x,2)
  if(any(shape(x)/=shape(ref)).or.n<1)return
  bad=0
  if(.not.scalar_finite_complex(x))bad=1
  call sumgrid(bad)
  if(abs(bad(1,1))>0)return
  allocate(m(n,n),u(n,n),vh(n,n),r(n,n),work(max(1,8*n*n)),s(n),rw(5*n))
!$omp parallel do collapse(2) private(i,v)
  do k=1,n
  do j=1,n
   v=0
   do i=1,ng
    v=v+conjg(x(i,j))*ref(i,k)
   enddo
   m(j,k)=v
  enddo
  enddo
!$omp end parallel do
  call sumgrid(m)
  call zgesvd('A','A',n,n,m,n,s,u,n,vh,n,work,size(work),rw,info)
  if(info/=0)return
  if(minval(s)<=maxval(s)*1d-12)return
  r=matmul(u,vh)
  ! Thread-private row workspace is only n elements; no extra local-grid frame.
!$omp parallel private(row,i,j)
  allocate(row(n))
!$omp do
  do i=1,ng
   row=x(i,:)
   do j=1,n
    x(i,j)=sum(row*r(:,j))
   enddo
  enddo
!$omp end do
  deallocate(row)
!$omp end parallel
  status=0
 end subroutine
 subroutine history_accept(h,x,dv,step,sumgrid,status)
  implicit none
  type(s_factor_history),intent(inout)::h
  complex(real64),intent(inout)::x(:,:)
  real(real64),intent(in)::dv
  integer,intent(in)::step
  procedure(sum_callback)::sumgrid
  integer,intent(out)::status
  complex(real64)::d(2),y,acc(6,1)
  real(real64)::g(2,2),b(2),gg(2,2),rhs(2,1),ridge
  real(real64)::af(2,2),scales(2),solution(2,1),rcond,ferr(1),berr(1),solve_work(6)
  integer::solve_iwork(2)
  character::equed
  external::dposvx
  integer::i,j,k,l,info
  status=1
  if(step<0.or.dv<=0.or.h%interval<1)return
  if(.not.scalar_finite_complex(x))return
  if(h%count>0)then
   if(step-h%last/=h%interval)return
   if(any(shape(x)/=shape(h%x(:,:,1))))return
   call factor_align(x,h%x(:,:,1),sumgrid,status)
   if(status/=0)return
  else
   allocate(h%x(size(x,1),size(x,2),3));h%x=0
  endif
  status=1
  if(h%count==3)then
   g=0;b=0
!$omp parallel do collapse(2) private(d,y,k,l) reduction(+:g,b)
   do j=1,size(x,2)
   do i=1,size(x,1)
    d=0
    do k=1,2
     d(k)=h%x(i,j,k)-h%x(i,j,k+1)
    enddo
    y=x(i,j)-h%x(i,j,1)
    do k=1,2
     b(k)=b(k)+real(conjg(d(k))*y,real64)
     do l=1,2
      g(k,l)=g(k,l)+real(conjg(d(k))*d(l),real64)
     enddo
    enddo
   enddo
   enddo
!$omp end parallel do
   acc(:,1)=cmplx([reshape(g,[4]),b]*dv,0d0,real64);call sumgrid(acc)
   h%g=h%g+reshape(real(acc(1:4,1),real64),[2,2]);h%b=h%b+real(acc(5:6,1),real64)
   h%teachers=h%teachers+1;ridge=0
   do k=1,2
    ridge=ridge+h%g(k,k)
   enddo
   ridge=1d-8*ridge;gg=h%g;rhs(:,1)=h%b
   do k=1,2
    gg(k,k)=gg(k,k)+ridge
   enddo
   h%ready=.false.;info=1
   ! Solve the regularized Gram system directly; avoid determinant underflow.
   call dposvx('E','U',2,1,gg,2,af,2,equed,scales,rhs,2,solution,2, &
     rcond,ferr,berr,solve_work,solve_iwork,info)
   if(info==0)then
    rhs=solution
    if(scalar_finite_real(rhs(:,1)).and.sum(abs(rhs(:,1)))<32d0)then
     h%coeff=rhs(:,1);h%ready=h%teachers>=8
    endif
   endif
  endif
!$omp parallel do collapse(2) private(k)
  do j=1,size(x,2)
  do i=1,size(x,1)
   do k=3,2,-1
    h%x(i,j,k)=h%x(i,j,k-1)
   enddo
   h%x(i,j,1)=x(i,j)
  enddo
  enddo
!$omp end parallel do
  h%count=min(h%count+1,3);h%last=step;status=0
 end subroutine
 subroutine history_predict(h,horizon,x,status)
  implicit none
  type(s_factor_history),intent(in)::h
  integer,intent(in)::horizon
  complex(real64),intent(out)::x(:,:)
  integer,intent(out)::status
  real(real64)::f,a(2)
  integer::i,j,k
  status=1
  if(.not.h%ready.or.horizon<1.or.horizon>h%interval)return
  if(any(shape(x)/=shape(h%x(:,:,1))))return
  f=real(horizon,real64)/real(h%interval,real64)
  a=f*f*h%coeff;a(1)=f+f*f*(h%coeff(1)-1d0)
!$omp parallel do collapse(2) private(k)
  do j=1,size(x,2)
  do i=1,size(x,1)
   x(i,j)=h%x(i,j,1)
   do k=1,2
    x(i,j)=x(i,j)+a(k)*(h%x(i,j,k)-h%x(i,j,k+1))
   enddo
  enddo
  enddo
!$omp end parallel do
  if(.not.scalar_finite_complex(x))return
  status=0
 end subroutine
end module
