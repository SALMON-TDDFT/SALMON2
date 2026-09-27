! Periodic rVV10 nonlocal correlation in Hartree atomic units.
! Sabatini, Gorni and de Gironcoli, PRB 87, 041108(R) (2013).
! Spline-channel convolution; analytic 3D Fourier transform of the rational
! kernel avoids a radial table, finite radial cutoff, and periodic image sums.
! The q saturation is the standard 12-term rVV10 interpolation regularization.
module rvv10
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: rvv10_evaluate,rvv10_kernel_fourier
  include 'fftw3.f03'
contains
  real(8) function rvv10_kernel_fourier(q1,q2,g) result(value)
    real(8),intent(in) :: q1,q2,g
    real(8) :: a,b,c,pi,coef(3),r(3)
    pi=acos(-1d0)
    a=1/sqrt(q1);b=1/sqrt(q2);c=sqrt(2/(q1+q2))
    if(q1==q2)then
      value=-.75d0*pi*pi/(4*q1**1.5d0)*(1+g*a)*exp(-g*a)
    else if(g==0d0)then
      value=-3*pi*pi/(q1*q2*(q1+q2)*(a+b)*(a+c)*(b+c))
    else
      r=[a,b,c]
      coef(1)=1/((b*b-a*a)*(c*c-a*a))
      coef(2)=1/((a*a-b*b)*(c*c-b*b))
      coef(3)=1/((a*a-c*c)*(b*b-c*c))
      value=-3*pi*pi/(q1*q2*(q1+q2)*g)*sum(coef*exp(-g*r))
    endif
  end function

  ! Return energy per volume, partial dE/d(rho) and dE/d(sigma).
  ! Caller adds -div(2*vsigma*grad rho) using the adjoint of its density gradient.
  subroutine rvv10_evaluate(n,h,rho,sigma,b,c,nq,energy,vrho,vsigma,status)
    integer,intent(in) :: n(3),nq
    real(8),intent(in) :: h(3),rho(:),sigma(:),b,c
    real(8),intent(out) :: energy(:),vrho(:),vsigma(:)
    integer,intent(out) :: status
    integer :: ng,i,j,a,d,x,y,z,ig,p(3)
    real(8) :: pi,beta,q,qn,qs,amp,fac,v1,v2,kappa,w,wg,exponent,ds,t,power,g,phi
    real(8),allocatable :: mesh(:),second(:,:),basis(:,:),deriv(:,:),qn_all(:),qs_all(:),amplitude(:)
    complex(c_double_complex),allocatable :: theta(:,:),u(:,:),work(:)
    type(c_ptr) :: forward,backward
    status=1;ng=product(n)
    if(any(n<1).or.nq<8.or.nq>128)return
    if(any(h<=0d0).or..not.all(ieee_is_finite(h)))return
    if(.not.ieee_is_finite(b).or..not.ieee_is_finite(c).or.b<=0d0.or.c<0d0)return
    if(size(rho)/=ng.or.size(sigma)/=ng.or.size(energy)/=ng.or.size(vrho)/=ng.or.size(vsigma)/=ng)return
    if(any(rho<0d0).or.any(sigma<0d0))return
    if(.not.all(ieee_is_finite(rho)).or..not.all(ieee_is_finite(sigma)))return
    pi=acos(-1d0);beta=(3/b**2)**.75d0/32
    allocate(mesh(nq),second(nq,nq),basis(ng,nq),deriv(ng,nq),qn_all(ng),qs_all(ng),amplitude(ng))
    allocate(theta(ng,nq),u(ng,nq),work(ng))
    do a=1,nq
      mesh(a)=exp(log(1d-4)+real(a-1,8)/(nq-1)*log(.5d0/1d-4))
    enddo
    call spline_second(mesh,second)
    do i=1,ng
      q=.5d0;qn=0d0;qs=0d0;amp=0d0
      if(rho(i)>1d-18)then
        kappa=1.5d0*b*pi*(rho(i)/(9*pi))**(1d0/6)
        wg=c*(sigma(i)/rho(i)**2)**2
        w=sqrt(4*pi*rho(i)/3+wg)
        q=w/kappa
        qn=( (4*pi/3-4*wg/rho(i))/(2*w)-w/(6*rho(i)) )/kappa
        qs=c*sigma(i)/(rho(i)**4*w*kappa)
        t=q/.5d0
        if(t>=2d0)then
          q=.5d0;qn=0d0;qs=0d0
        else
          exponent=0d0;ds=0d0;power=1d0
          do j=1,12
            ds=ds+power;power=power*t;exponent=exponent+power/j
          enddo
          fac=exp(-exponent)
          q=.5d0*(1-fac);qn=qn*ds*fac;qs=qs*ds*fac
        endif
        if(q<mesh(1))then
          q=mesh(1);qn=0d0;qs=0d0
        endif
        amp=rho(i)/kappa**1.5d0
      endif
      call spline_basis(mesh,second,q,basis(i,:),deriv(i,:))
      amplitude(i)=amp;qn_all(i)=qn;qs_all(i)=qs
      theta(i,:)=amp*basis(i,:)
    enddo
    forward=fftw_plan_dft_3d(n(3),n(2),n(1),work,work,FFTW_FORWARD,FFTW_ESTIMATE)
    backward=fftw_plan_dft_3d(n(3),n(2),n(1),work,work,FFTW_BACKWARD,FFTW_ESTIMATE)
    if(.not.c_associated(forward).or..not.c_associated(backward))then
      if(c_associated(forward))call fftw_destroy_plan(forward)
      if(c_associated(backward))call fftw_destroy_plan(backward)
      return
    endif
    do a=1,nq
      work=theta(:,a);call fftw_execute_dft(forward,work,work);theta(:,a)=work
    enddo
    u=0d0
    ! O(nq**2 * ng) work, O(nq * ng) storage; no ng**2 pair matrix.
!$omp parallel do collapse(3) private(x,y,z,p,ig,g,a,d,phi)
    do z=0,n(3)-1;do y=0,n(2)-1;do x=0,n(1)-1
      p=[x,y,z];where(p>=(n+1)/2)p=p-n
      ig=1+x+n(1)*(y+n(2)*z);g=sqrt(sum((2*pi*p/(n*h))**2))
      do a=1,nq
        do d=1,a
          phi=rvv10_kernel_fourier(mesh(a),mesh(d),g)
          u(ig,a)=u(ig,a)+phi*theta(ig,d)
          if(a/=d)u(ig,d)=u(ig,d)+phi*theta(ig,a)
        enddo
      enddo
    enddo;enddo;enddo
!$omp end parallel do
    do a=1,nq
      work=u(:,a);call fftw_execute_dft(backward,work,work);u(:,a)=work/ng
    enddo
    call fftw_destroy_plan(forward);call fftw_destroy_plan(backward)
    energy=beta*rho;vrho=beta;vsigma=0d0
    do i=1,ng
      if(rho(i)<=1d-18)cycle
      v1=sum(basis(i,:)*real(u(i,:),8));v2=sum(deriv(i,:)*real(u(i,:),8))
      amp=amplitude(i)
      energy(i)=energy(i)+.5d0*amp*v1
      vrho(i)=vrho(i)+.75d0*amp/rho(i)*v1+amp*qn_all(i)*v2
      vsigma(i)=amp*qs_all(i)*v2
    enddo
    if(.not.all(ieee_is_finite(energy)).or..not.all(ieee_is_finite(vrho)).or. &
       .not.all(ieee_is_finite(vsigma)))return
    status=0
  end subroutine

  subroutine spline_second(x,second)
    real(8),intent(in) :: x(:)
    real(8),intent(out) :: second(:,:)
    real(8) :: lower(size(x)),diag(size(x)),upper(size(x)),rhs(size(x)),y(size(x)),factor
    integer :: n,j,i
    n=size(x);second=0d0
    do j=1,n
      y=0d0;y(j)=1d0;diag=1d0;lower=0d0;upper=0d0;rhs=0d0
      do i=2,n-1
        lower(i)=x(i)-x(i-1);upper(i)=x(i+1)-x(i)
        diag(i)=2*(lower(i)+upper(i))
        rhs(i)=6*((y(i+1)-y(i))/upper(i)-(y(i)-y(i-1))/lower(i))
      enddo
      do i=2,n
        factor=lower(i)/diag(i-1)
        diag(i)=diag(i)-factor*upper(i-1);rhs(i)=rhs(i)-factor*rhs(i-1)
      enddo
      second(n,j)=rhs(n)/diag(n)
      do i=n-1,1,-1
        second(i,j)=(rhs(i)-upper(i)*second(i+1,j))/diag(i)
      enddo
    enddo
  end subroutine

  subroutine spline_basis(x,second,q,basis,deriv)
    real(8),intent(in) :: x(:),second(:,:),q
    real(8),intent(out) :: basis(:),deriv(:)
    real(8) :: h,a,b
    integer :: lo,hi,j
    lo=1;hi=size(x)
    do while(hi-lo>1)
      j=(lo+hi)/2
      if(q<x(j))then
        hi=j
      else
        lo=j
      endif
    enddo
    h=x(hi)-x(lo);a=(x(hi)-q)/h;b=1-a
    basis=((a**3-a)*second(lo,:)+(b**3-b)*second(hi,:))*h*h/6
    deriv=(-(3*a*a-1)*second(lo,:)+(3*b*b-1)*second(hi,:))*h/6
    basis(lo)=basis(lo)+a;basis(hi)=basis(hi)+b
    deriv(lo)=deriv(lo)-1/h;deriv(hi)=deriv(hi)+1/h
  end subroutine
end module
