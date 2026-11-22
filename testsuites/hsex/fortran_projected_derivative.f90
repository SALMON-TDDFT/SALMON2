program test_projected_derivative
  use hsex_projected_derivative, only: projected_derivative,rank1_second
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  integer,parameter :: ng=6,na=2
  integer :: box(ng),i,status
  real(8) :: r(ng),hvol,a,delta
  complex(8) :: phase(ng),x(ng),s0(na,na),s1(na,na),s2(na,na),sigma(na,na),d1(na,na),d2(na,na)
  complex(8) :: velocity(ng),mass(ng),hp(ng),hm(ng),h0(ng)
  integer :: indices(6)
  real(8) :: coords(6)
  complex(8) :: beta(6),ppx(4),ppmass(4),ppplus(4),ppminus(4),ppzero(4)
  box=[1,1,2,2,2,1];r=[-.4d0,.2d0,.5d0,.9d0,1.2d0,-.7d0];hvol=.37d0;a=.17d0;delta=1d-4
  s0=reshape([cmplx(1d0,0d0,8),cmplx(.2d0,-.1d0,8),cmplx(.2d0,.1d0,8),cmplx(.7d0,0d0,8)],[na,na])
  s1=reshape([cmplx(.3d0,0d0,8),cmplx(-.1d0,-.2d0,8),cmplx(-.1d0,.2d0,8),cmplx(.1d0,0d0,8)],[na,na])
  s2=reshape([cmplx(.2d0,0d0,8),cmplx(.05d0,-.03d0,8),cmplx(.05d0,.03d0,8),cmplx(-.2d0,0d0,8)],[na,na])
  do i=1,ng
    x(i)=cmplx(sin(dble(i)),cos(dble(i)),8)
  enddo
  phase=sqrt(hvol/3d0)*exp(cmplx(0d0,1d0,8)*a*r)
  sigma=s0+cos(1.3d0*a)*s1+sin(1.3d0*a)*s2
  d1=1.3d0*(-sin(1.3d0*a)*s1+cos(1.3d0*a)*s2)
  d2=-1.3d0**2*(cos(1.3d0*a)*s1+sin(1.3d0*a)*s2)
  call action(a+delta,hp);call action(a-delta,hm);call action(a,h0)
  call projected_derivative(box,phase,r,sigma,d1,hvol,x,velocity,status,d2,mass)
  if(status/=0)error stop 'projected derivative rejected finite input'
  if(maxval(abs(velocity-(hp-hm)/(2d0*delta)))>1d-8)error stop 'HSEX velocity differs from finite-A H action'
  if(maxval(abs(mass-(hp-2d0*h0+hm)/delta**2))>1d-6)error stop 'HSEX mass differs from finite-A H action'
  indices=[1,2,2,3,4,4];coords=[-.4d0,.2d0,1.3d0,.9d0,-.2d0,1.1d0]
  do i=1,6
    beta(i)=cmplx(sin(dble(i)),cos(dble(2*i)),8)
  enddo
  ppx=x(:4)
  call pp_action(a+delta,ppplus);call pp_action(a-delta,ppminus);call pp_action(a,ppzero)
  call rank1_second(indices,coords,beta*exp(cmplx(0d0,-a,8)*coords),.7d0,ppx,ppmass,status)
  if(status/=0)error stop 'nonlocal projector mass rejected'
  if(maxval(abs(ppmass-(ppplus-2d0*ppzero+ppminus)/delta**2))>2d-6) &
    error stop 'nonlocal mass image/scatter/phase error'
  beta(1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call rank1_second(indices,coords,beta,.7d0,ppx,ppmass,status)
  if(status==0)error stop 'NaN nonlocal projector accepted'
  x(1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call projected_derivative(box,phase,r,sigma,d1,hvol,x,velocity,status,d2,mass)
  if(status==0)error stop 'NaN projected derivative accepted'
  x(1)=cmplx(ieee_value(0d0,ieee_positive_inf),0d0,8)
  call projected_derivative(box,phase,r,sigma,d1,hvol,x,velocity,status,d2,mass)
  if(status==0)error stop 'infinite projected derivative accepted'
  print *, 'PROJECTED_DERIVATIVE_PASS'
contains
  subroutine pp_action(value,y)
    implicit none
    real(8),intent(in) :: value
    complex(8),intent(out) :: y(4)
    complex(8) :: f(4)
    integer :: g
    f=0d0
    do g=1,6
      f(indices(g))=f(indices(g))+beta(g)*exp(cmplx(0d0,-value*coords(g),8))
    enddo
    y=.7d0*f*dot_product(f,ppx)
  end subroutine
  subroutine action(value,y)
    implicit none
    real(8),intent(in) :: value
    complex(8),intent(out) :: y(ng)
    complex(8) :: f(ng),s(na,na),px(na)
    integer :: g
    f=sqrt(hvol/3d0)*exp(cmplx(0d0,1d0,8)*value*r)
    s=s0+cos(1.3d0*value)*s1+sin(1.3d0*value)*s2;px=0d0
    do g=1,ng
      px(box(g))=px(box(g))+f(g)*x(g)
    enddo
    px=matmul(s,px)
    y=conjg(f)*px(box)/hvol
  end subroutine
end program
