! Passive isotropic LD fitting. Local workspaces only; no solver/MPI lifetime.
module lorentz_drude_fit
  use iso_c_binding
  use,intrinsic :: ieee_arithmetic,only:ieee_is_finite
  use salmonpy_convergence_kernel,only:scaled_square
  implicit none
contains
  logical function residuals(n,p,x,target,weights,par,r,cost) result(ok)
    integer,intent(in)::n,p
    real(c_double),intent(in)::x(n),weights(n),par(1+3*p)
    complex(c_double_complex),intent(in)::target(n)
    real(c_double),intent(out)::r(2*n),cost
    complex(c_double_complex)::z,den
    integer::i,q
    ok=.false.;cost=0d0
    do q=1,1+3*p
      if(.not.ieee_is_finite(par(q)))return
    enddo
    do i=1,n
      z=cmplx(par(1),0d0,kind=c_double)
      do q=1,p
        if(par(3*q-1)==0d0)cycle
        if(abs(x(i))==par(3*q).and.(par(3*q+1)==0d0.or.x(i)==0d0))return
        den=cmplx(par(3*q)**2-x(i)**2,-par(3*q+1)*x(i),kind=c_double)
        if(.not.ieee_is_finite(real(den,kind=c_double)).or..not.ieee_is_finite(aimag(den)))return
        if(abs(den)==0d0)return
        z=z+par(3*q-1)/den
      enddo
      r(i)=real(z-target(i),kind=c_double)*weights(i)
      r(n+i)=aimag(z-target(i))*weights(i)
      if(.not.ieee_is_finite(r(i)).or..not.ieee_is_finite(r(n+i)))return
      cost=cost+r(i)**2+r(n+i)**2
      if(.not.ieee_is_finite(cost))return
    enddo
    ok=.true.
  end function

  logical function solve_spd(a,b,m) result(ok)
    integer,intent(in)::m
    real(c_double),intent(inout)::a(m,m),b(m)
    real(c_double)::v
    integer::i,j,k
    ok=.false.
    do j=1,m
      v=a(j,j)-sum(a(j,1:j-1)**2)
      if(.not.ieee_is_finite(v).or.v<=0d0)return
      a(j,j)=sqrt(v)
      do i=j+1,m
        a(i,j)=(a(i,j)-sum(a(i,1:j-1)*a(j,1:j-1)))/a(j,j)
        if(.not.ieee_is_finite(a(i,j)))return
      enddo
    enddo
    do i=1,m
      b(i)=(b(i)-sum(a(i,1:i-1)*b(1:i-1)))/a(i,i)
      if(.not.ieee_is_finite(b(i)))return
    enddo
    do i=m,1,-1
      b(i)=(b(i)-sum(a(i+1:m,i)*b(i+1:m)))/a(i,i)
      if(.not.ieee_is_finite(b(i)))return
    enddo
    ok=.true.
  end function

  integer(c_int) function sp_ld_fit(n,p,fit_eps,drude,maxiter,starts,has_initial,fixed_eps, &
                                    omega,target,initial,parameters,stats,details) bind(C) result(status)
    integer(c_int),value::n,p,fit_eps,drude,maxiter,starts,has_initial
    real(c_double),value::fixed_eps
    real(c_double),intent(in)::omega(*),initial(*)
    complex(c_double_complex),intent(in)::target(*)
    real(c_double),intent(out)::parameters(*),stats(4)
    integer(c_int),intent(out)::details(3)
    real(c_double),allocatable::x(:),weights(:),par(:),trial(:),best(:),r(:),rt(:),jac(:,:),a(:,:),b(:)
    logical,allocatable::free(:)
    real(c_double)::scale,cost,newcost,bestcost,lambda,h,oldvalue,es,eq,rs,rq,fast,magnitude
    integer::m,i,j,q,seed,iter,ios,bestiter,bestconv,conv
    logical::ok,columns_ok
    status=1;stats=0d0;details=0
    if(n<2.or.n>huge(n)/2.or.p<1.or.p>(huge(p)-1)/3.or.maxiter<1.or.starts<1)return
    if(.not.ieee_is_finite(fixed_eps).or.fixed_eps<=0d0)return
    m=1+3*p;scale=0d0
    ! Reject address-size overflow before workspace products or allocations.
    if(int(m,c_size_t)>huge(0_c_size_t)/8_c_size_t/int(m,c_size_t))return
    if(int(n,c_size_t)>huge(0_c_size_t)/16_c_size_t/int(m,c_size_t))return
    do i=1,n
      if(.not.ieee_is_finite(omega(i)).or.omega(i)<0d0)return
      if(drude/=0.and.omega(i)==0d0)return
      if(.not.ieee_is_finite(real(target(i),kind=c_double)).or..not.ieee_is_finite(aimag(target(i))))return
      if(.not.ieee_is_finite(abs(target(i))))return
      scale=max(scale,omega(i))
    enddo
    if(scale<=0d0)return
    status=3
    allocate(x(n),weights(n),par(m),trial(m),best(m),r(2*n),rt(2*n),jac(2*n,m),a(m,m),b(m),free(m),stat=ios)
    if(ios/=0)return
    x=omega(1:n)/scale
    do i=1,n
      weights(i)=1d0/max(1d0,abs(target(i)))
    enddo
    free=.true.;free(1)=fit_eps/=0
    if(drude/=0)free(3)=.false.
    bestcost=huge(1d0);bestiter=0;bestconv=0
    do seed=1,starts
      par(1)=fixed_eps
      do q=1,p
        par(3*q-1)=.15d0/dble(p)
        par(3*q)=.2d0+.7d0*(dble(q)-.5d0)/dble(p)+.025d0*dble(seed-1)
        par(3*q+1)=.03d0+.08d0*dble(seed-1)
      enddo
      if(drude/=0)then
        par(3)=0d0
        par(2)=max(.005d0,abs(target(1))*x(1)**2)
      endif
      if(has_initial/=0.and.seed==1)then
        par(1)=initial(1)
        do q=1,p
          par(3*q-1)=initial(3*q)*(initial(2)/scale)**2
          par(3*q)=initial(3*q+1)/scale
          par(3*q+1)=initial(3*q+2)/scale
        enddo
      endif
      if(fit_eps==0)par(1)=fixed_eps
      if(drude/=0)par(3)=0d0
      ok=residuals(n,p,x,target(1:n),weights,par,r,cost)
      if(.not.ok)cycle
      lambda=.01d0;conv=0
      do iter=1,maxiter
        jac=0d0;columns_ok=.true.
        do j=1,m
          if(.not.free(j))cycle
          oldvalue=par(j);h=1d-6*(1d0+abs(oldvalue));par(j)=oldvalue+h
          ok=residuals(n,p,x,target(1:n),weights,par,rt,newcost)
          par(j)=oldvalue
          if(.not.ok)then
            columns_ok=.false.;exit
          endif
          jac(:,j)=(rt-r)/h
        enddo
        if(.not.columns_ok)exit
        a=matmul(transpose(jac),jac);b=-matmul(transpose(jac),r)
        if(maxval(abs(b))<1d-10)then
          conv=1;exit
        endif
        do j=1,m
          a(j,j)=a(j,j)+lambda*(a(j,j)+1d-12)
          if(.not.free(j))a(j,j)=1d0
        enddo
        ok=solve_spd(a,b,m)
        if(ok)then
          trial=max(0d0,par+b)
          if(free(1))trial(1)=max(1d-9,trial(1))
          if(.not.free(1))trial(1)=fixed_eps
          if(drude/=0)trial(3)=0d0
          ok=residuals(n,p,x,target(1:n),weights,trial,rt,newcost)
        endif
        if(ok)then
          if(newcost<cost)then
            par=trial;r=rt;cost=newcost;lambda=max(1d-12,lambda*.3d0)
            if(cost<1d-18*dble(n))then
              conv=1;exit
            endif
            cycle
          endif
        endif
        lambda=lambda*10d0
        if(lambda>1d16)exit
      enddo
      if(cost<bestcost)then
        bestcost=cost;best=par;bestiter=min(iter,maxiter);bestconv=conv
      endif
    enddo
    status=2
    if(bestcost==huge(1d0))return
    ok=residuals(n,p,x,target(1:n),weights,best,r,cost)
    if(.not.ok)return
    es=0d0;eq=1d0;rs=0d0;rq=1d0
    do i=1,n
      call scaled_square(r(i)/weights(i),es,eq)
      call scaled_square(r(n+i)/weights(i),es,eq)
      call scaled_square(real(target(i),kind=c_double),rs,rq)
      call scaled_square(aimag(target(i)),rs,rq)
    enddo
    stats(1)=huge(1d0)
    if(rs>0d0)stats(1)=(es/rs)*sqrt(eq/rq)
    if(es==0d0)stats(1)=0d0
    if(.not.ieee_is_finite(stats(1)))stats(1)=huge(1d0)
    stats(2)=sqrt(cost/dble(n));stats(3)=cost
    parameters(1)=best(1);magnitude=0d0;fast=0d0
    do q=1,p
      parameters(3*q-1)=best(3*q-1)*scale**2
      parameters(3*q)=best(3*q)*scale
      parameters(3*q+1)=best(3*q+1)*scale
      magnitude=magnitude+parameters(3*q-1)
      if(parameters(3*q-1)>0d0)fast=max(fast,parameters(3*q),parameters(3*q+1))
    enddo
    do i=1,m
      if(.not.ieee_is_finite(parameters(i)).or.parameters(i)<0d0)return
    enddo
    fast=max(fast,sqrt(magnitude))
    if(fast>0d0)stats(4)=.1d0/fast
    do i=2,4
      if(.not.ieee_is_finite(stats(i)))return
    enddo
    details=[bestiter,bestconv,starts];status=0
  end function

  integer(c_int) function sp_ld_nk(n,nreal,k,eps) bind(C) result(status)
    integer(c_int),value::n
    real(c_double),intent(in)::nreal(*),k(*)
    complex(c_double_complex),intent(out)::eps(*)
    integer::i
    status=1
    if(n<1)return
    do i=1,n
      if(.not.ieee_is_finite(nreal(i)).or..not.ieee_is_finite(k(i)))return
      if(nreal(i)<0d0.or.k(i)<0d0)return
      eps(i)=cmplx(nreal(i),k(i),kind=c_double)**2
      if(.not.ieee_is_finite(real(eps(i),kind=c_double)).or..not.ieee_is_finite(aimag(eps(i))))return
    enddo
    status=0
  end function
  integer(c_int) function sp_ld_band_fraction(n,m,omega,values,low,high,stats) bind(C) result(status)
    integer(c_int),value::n,m
    real(c_double),intent(in)::omega(*)
    complex(c_double_complex),intent(in)::values(m,*)
    real(c_double),value::low,high
    real(c_double),intent(out)::stats(3)
    real(c_double)::peak,total,outside,weight
    integer::i,j
    status=1;stats=0d0;peak=0d0
    if(n<1.or.m<1)return
    if(.not.ieee_is_finite(low).or..not.ieee_is_finite(high))return
    if(low<0d0.or.high<=low)return
    do i=1,n
      if(.not.ieee_is_finite(omega(i)).or.omega(i)<0d0)return
      do j=1,m
        if(.not.ieee_is_finite(real(values(j,i),kind=c_double)))return
        if(.not.ieee_is_finite(aimag(values(j,i))))return
        if(.not.ieee_is_finite(abs(values(j,i))))return
        peak=max(peak,abs(values(j,i)))
      enddo
    enddo
    total=0d0;outside=0d0
    if(peak>0d0)then
      do i=1,n
        do j=1,m
          weight=(abs(values(j,i))/peak)**2
          total=total+weight
          if(omega(i)<low.or.omega(i)>high)outside=outside+weight
        enddo
      enddo
      stats(1)=outside/total
    endif
    stats(2)=total;stats(3)=peak;status=0
  end function
  integer(c_int) function sp_ld_wavelength(n,wavelength,omega,energy) bind(C) result(status)
    integer(c_int),value::n
    real(c_double),intent(in)::wavelength(*)
    real(c_double),intent(out)::omega(*),energy(*)
    ! Same CODATA-2018 values as SALMON math/phys_constants.f90.
    real(c_double),parameter::c_au=137.035999084d0,a0_um=5.29177210903d-5
    real(c_double),parameter::hartree_ev=4.3597447222071d-18/1.602176634d-19
    integer::i
    status=1
    if(n<1)return
    do i=1,n
      if(.not.ieee_is_finite(wavelength(i)).or.wavelength(i)<=0d0)return
      omega(i)=2d0*acos(-1d0)*c_au*a0_um/wavelength(i)
      energy(i)=omega(i)*hartree_ev
      if(.not.ieee_is_finite(omega(i)).or..not.ieee_is_finite(energy(i)))return
    enddo
    status=0
  end function

  integer(c_int) function sp_ld_fit_seed(n,p,fit_eps,drude,seed,fixed_eps,omega,target,par) bind(C) result(status)
    integer(c_int),value::n,p,fit_eps,drude,seed
    real(c_double),value::fixed_eps
    real(c_double),intent(in)::omega(n)
    complex(c_double_complex),intent(in)::target(n)
    real(c_double),intent(out)::par(1+3*p)
    real(c_double)::scale
    integer::i,q
    status=1
    if(n<2.or.p<1.or.seed<1)return
    scale=0d0
    do i=1,n
      if(.not.ieee_is_finite(omega(i)).or.omega(i)<0d0)return
      scale=max(scale,omega(i))
    enddo
    if(scale<=0d0)return
    par(1)=fixed_eps
    do q=1,p
      par(3*q-1)=.15d0/dble(p)
      par(3*q)=.2d0+.7d0*(dble(q)-.5d0)/dble(p)+.025d0*dble(seed-1)
      par(3*q+1)=.03d0+.08d0*dble(seed-1)
    enddo
    if(drude/=0)then
      par(3)=0d0
      par(2)=max(.005d0,abs(target(1))*(omega(1)/scale)**2)
    endif
    status=0
  end function

  integer(c_int) function sp_ld_fit_evaluate(n,p,omega,target,par,r,parameters,stats) bind(C) result(status)
    integer(c_int),value::n,p
    real(c_double),intent(in)::omega(n),par(1+3*p)
    complex(c_double_complex),intent(in)::target(n)
    real(c_double),intent(out)::r(2*n),parameters(1+3*p),stats(4)
    real(c_double)::x(n),weights(n),scale,cost,es,eq,rs,rq,fast,magnitude
    integer::i,q
    status=1;stats=0d0;scale=0d0
    if(n<2.or.p<1)return
    do i=1,n
      if(.not.ieee_is_finite(omega(i)).or.omega(i)<0d0)return
      if(.not.ieee_is_finite(real(target(i),kind=c_double)))return
      if(.not.ieee_is_finite(aimag(target(i))))return
      scale=max(scale,omega(i))
      weights(i)=1d0/max(1d0,abs(target(i)))
    enddo
    if(scale<=0d0)return
    do i=1,n
      x(i)=omega(i)/scale
    enddo
    if(.not.residuals(n,p,x,target,weights,par,r,cost))return
    es=0d0;eq=0d0;rs=0d0;rq=0d0
    do i=1,n
      call scaled_square(r(i)/weights(i),es,eq)
      call scaled_square(r(n+i)/weights(i),es,eq)
      call scaled_square(real(target(i),kind=c_double),rs,rq)
      call scaled_square(aimag(target(i)),rs,rq)
    enddo
    stats(1)=huge(1d0)
    if(rs>0d0)stats(1)=(es/rs)*sqrt(eq/rq)
    if(es==0d0)stats(1)=0d0
    stats(2)=sqrt(cost/dble(n));stats(3)=cost
    parameters(1)=par(1);fast=0d0;magnitude=0d0
    do q=1,p
      parameters(3*q-1)=par(3*q-1)*scale**2
      parameters(3*q)=par(3*q)*scale
      parameters(3*q+1)=par(3*q+1)*scale
      magnitude=magnitude+parameters(3*q-1)
      if(parameters(3*q-1)>0d0)fast=max(fast,parameters(3*q),parameters(3*q+1))
    enddo
    fast=max(fast,sqrt(magnitude))
    if(fast>0d0)stats(4)=.1d0/fast
    do i=1,4
      if(.not.ieee_is_finite(stats(i)))return
    enddo
    status=0
  end function

  integer(c_int) function sp_ld_tensor(n,scalar,tensor) bind(C) result(status)
    integer(c_int),value::n
    complex(c_double_complex),intent(in)::scalar(n)
    complex(c_double_complex),intent(out)::tensor(3,3,n)
    integer::i,j,k
    status=1
    if(n<1)return
    do k=1,n
      if(.not.ieee_is_finite(real(scalar(k),kind=c_double)))return
      if(.not.ieee_is_finite(aimag(scalar(k))))return
      do j=1,3
        do i=1,3
          tensor(i,j,k)=cmplx(0d0,0d0,kind=c_double)
          if(i==j)tensor(i,j,k)=scalar(k)
        enddo
      enddo
    enddo
    status=0
  end function
end module
