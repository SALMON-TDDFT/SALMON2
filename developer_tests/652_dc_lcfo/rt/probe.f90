program probe
  use, intrinsic :: ieee_arithmetic
  use lcfo_rt_core
  implicit none
  complex(8) :: h(3,3),u(3,3),c(3,2),out(3,2),back(3,2),exact(3,2),rho(3,3),r0(3,3)
  complex(8) :: a(3,2),b(3,2),ref(3,2),bad(3,3),ident(3,3)
  real(8) :: eig(3),f(2),dt,pi,err1,err2,res,density(3),potential(3),energy_grid,energy_coeff,last_time
  complex(8) :: projected(3,3),rotated(3,2)
  integer :: i,j,k,status,it
  pi=acos(-1d0); eig=[-0.7d0,0.2d0,1.1d0]; f=[2d0,0.4d0]
  do i=1,3
    do j=1,3
      u(i,j)=exp(cmplx(0d0,2*pi*(i-1)*(j-1)/3,8))/sqrt(3d0)
    end do
  end do
  h=(0d0,0d0);ident=(0d0,0d0)
  do i=1,3
    h=h+eig(i)*spread(u(:,i),2,3)*spread(conjg(u(:,i)),1,3)
    ident(i,i)=1d0
  end do
  c=ident(:,:2);dt=0.08d0
  call lcfo_cayley_step(h,c,dt,out,status)
  call require(status==0,'Cayley accepted')
  exact=(0d0,0d0)
  do i=1,3
    exact=exact+exp(cmplx(0d0,-dt*eig(i),8))* &
      matmul(spread(u(:,i),2,3)*spread(conjg(u(:,i)),1,3),c)
  end do
  err1=maxval(abs(out-exact))
  call lcfo_cayley_step(h,c,dt/2,a,status)
  call lcfo_cayley_step(h,a,dt/2,b,status)
  err2=maxval(abs(b-exact))
  call require(err1/err2>3.9d0.and.err1/err2<4.1d0,'second order vs analytic phases')
  call require(maxval(abs(matmul(conjg(transpose(out)),out)-ident(:2,:2)))<2d-14,'orthogonality')
  call lcfo_cayley_step(h,out,-dt,back,status)
  call require(maxval(abs(back-c))<2d-14,'time reversal')
  call lcfo_density(out,f,rho,status)
  call require(status==0,'fractional density accepted')
  call require(abs(sum([(real(rho(i,i),8),i=1,3)])-sum(f))<2d-14,'particle number')
  call require(maxval(abs(rho-conjg(transpose(rho))))<2d-14,'density Hermitian')
  call require(maxval(abs(aimag(rho)))>1d-4,'complex coherence retained')
  c=u(:,:2)
  call lcfo_density(c,f,r0,status)
  call lcfo_cayley_step(h,c,1d0,out,status)
  call lcfo_density(out,f,rho,status)
  call require(maxval(abs(rho-r0))<2d-14,'stationary fractional density')
  bad=h;bad(1,2)=bad(1,2)+0.01d0
  call lcfo_cayley_step(bad,c,dt,out,status)
  call require(status/=0.and.maxval(abs(out-c))==0d0,'reject nonhermitian without modifying state')
  bad=h;bad(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call lcfo_cayley_step(bad,c,dt,out,status)
  call require(status/=0.and.maxval(abs(out-c))==0d0,'reject nonfinite')
  call lcfo_density(c,[-0.1d0,1d0],rho,status)
  call require(status/=0,'reject negative occupations')
  c=ident(:,:2)
  call lcfo_midpoint_step(c,f,0d0,dt,nonlinear,1d-13,100,out,res,it,status)
  call require(status==0.and.it>1.and.res<1d-13,'nonlinear midpoint convergence')
  call lcfo_midpoint_step(out,f,dt,-dt,nonlinear,1d-13,100,back,res,it,status)
  call require(status==0.and.maxval(abs(back-c))<1d-12,'nonlinear time reversal')
  call lcfo_midpoint_step(c,f,0d0,dt,nonlinear,1d-15,1,back,res,it,status)
  call require(status/=0.and.maxval(abs(back-c))==0d0,'failed convergence preserves state')
  call lcfo_midpoint_step(c,f,0d0,dt,broken,1d-13,100,back,res,it,status)
  call require(status/=0.and.maxval(abs(back-c))==0d0,'callback failure preserves state')
  ref=c
  do k=1,32
    call lcfo_midpoint_step(ref,f,(k-1)*dt/32,dt/32,nonlinear,1d-13,100,a,res,it,status)
    call require(status==0,'reference converges');ref=a
  end do
  call lcfo_midpoint_step(c,f,0d0,dt/2,nonlinear,1d-13,100,a,res,it,status)
  call lcfo_midpoint_step(a,f,dt/2,dt/2,nonlinear,1d-13,100,b,res,it,status)
  err1=maxval(abs(out-ref));err2=maxval(abs(b-ref))
  call require(err1/err2>3.8d0.and.err1/err2<4.2d0,'nonlinear second order')
  call lcfo_midpoint_step(c,f,1.2d0,dt,timed,1d-13,100,out,res,it,status)
  call require(status==0.and.abs(last_time-(1.2d0+dt/2))<1d-15,'callback evaluated at midpoint time')
  call lcfo_cayley_step(h*(1d0+1.2d0+dt/2),c,dt,exact,status)
  call require(maxval(abs(out-exact))<2d-14,'explicit time dependence')
  call lcfo_midpoint_step(c,f,0d0,dt,bad_h,1d-13,100,out,res,it,status)
  call require(status==2.and.maxval(abs(out-c))==0d0,'nonhermitian callback rejected')
  call lcfo_midpoint_step(c,f,0d0,dt,nan_h,1d-13,100,out,res,it,status)
  call require(status/=0.and.maxval(abs(out-c))==0d0,'nonfinite callback rejected')
  c=ident(:,:2);potential=[-0.4d0,0.7d0,1.3d0]
  call lcfo_grid_density(u/sqrt(0.5d0),c,f,density,status)
  call require(status==0.and.abs(sum(density)*0.5d0-sum(f))<2d-14,'grid density charge')
  call lcfo_project_potential(u/sqrt(0.5d0),potential,0.5d0,projected,status)
  call require(status==0,'potential projection accepted')
  energy_grid=sum(density*potential)*0.5d0
  energy_coeff=sum(f*real(sum(conjg(c)*matmul(projected,c),dim=1),8))
  call require(abs(energy_grid-energy_coeff)<2d-14,'grid/projected potential expectation parity')
  call require(maxval(abs(projected-conjg(transpose(projected))))<2d-14,'projected potential hermiticity')
  call lcfo_project_potential(u,potential,-1d0,projected,status)
  call require(status/=0,'reject invalid quadrature volume')
  print *, 'PASS lcfo_rt_core: analytic, conservation, reversal, nonlinear and rejection checks'
contains
  subroutine require(ok,name)
    logical,intent(in)::ok
    character(*),intent(in)::name
    if(.not.ok)then
      print *, 'FAIL ',name
      stop 1
    end if
  end subroutine
  subroutine nonlinear(p,t,hh,istat)
    complex(8),intent(in)::p(:,:)
    real(8),intent(in)::t
    complex(8),intent(out)::hh(:,:)
    integer,intent(out)::istat
    integer::n
    hh=h+0.2d0*p
    do n=1,3
      hh(n,n)=hh(n,n)+0.3d0*real(p(n,n),8)
    end do
    istat=0
  end subroutine
  subroutine timed(p,t,hh,istat)
    complex(8),intent(in)::p(:,:)
    real(8),intent(in)::t
    complex(8),intent(out)::hh(:,:)
    integer,intent(out)::istat
    last_time=t;hh=h*(1d0+t);istat=0
  end subroutine
  subroutine bad_h(p,t,hh,istat)
    complex(8),intent(in)::p(:,:)
    real(8),intent(in)::t
    complex(8),intent(out)::hh(:,:)
    integer,intent(out)::istat
    hh=h;hh(1,2)=hh(1,2)+.1d0;istat=0
  end subroutine
  subroutine nan_h(p,t,hh,istat)
    complex(8),intent(in)::p(:,:)
    real(8),intent(in)::t
    complex(8),intent(out)::hh(:,:)
    integer,intent(out)::istat
    hh=h;hh(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8);istat=0
  end subroutine
  subroutine broken(p,t,hh,istat)
    complex(8),intent(in)::p(:,:)
    real(8),intent(in)::t
    complex(8),intent(out)::hh(:,:)
    integer,intent(out)::istat
    hh=0;istat=1
  end subroutine
end program
