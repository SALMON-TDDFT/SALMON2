! Reciprocal passive tensor LD observables; all work is caller/handle owned.
module lorentz_drude_tensor
  use iso_c_binding
  use,intrinsic :: ieee_arithmetic,only:ieee_is_finite
  implicit none
contains
  logical function positive_matrix(a,strict) result(ok)
    real(c_double),intent(in)::a(3,3)
    logical,intent(in)::strict
    real(c_double)::b(3,3),scale,tol,minor,det
    integer::i,j
    ok=.false.;scale=0d0;tol=64d0*epsilon(1d0)
    do j=1,3
      do i=1,3
        if(.not.ieee_is_finite(a(i,j)))return
        scale=max(scale,abs(a(i,j)))
      enddo
    enddo
    if(scale==0d0)then
      ok=.not.strict;return
    endif
    b=a/scale
    do j=1,3
      do i=1,3
        if(abs(b(i,j)-b(j,i))>tol)return
      enddo
    enddo
    do i=1,3
      if(b(i,i)<-tol)return
    enddo
    do j=2,3
      do i=1,j-1
        minor=b(i,i)*b(j,j)-b(i,j)*b(j,i)
        if(minor < -tol)return
      enddo
    enddo
    det=b(1,1)*(b(2,2)*b(3,3)-b(2,3)*b(3,2)) &
       -b(1,2)*(b(2,1)*b(3,3)-b(2,3)*b(3,1)) &
       +b(1,3)*(b(2,1)*b(3,2)-b(2,2)*b(3,1))
    if(det < -tol)return
    if(strict)then
      if(b(1,1)<=0d0.or.b(1,1)*b(2,2)-b(1,2)*b(2,1)<=0d0.or.det<=0d0)return
    endif
    ok=.true.
  end function

  integer(c_int) function sp_ld_tensor_validate(n,wp,f,w,g,eps) bind(C) result(status)
    integer(c_int),value::n
    real(c_double),value::wp
    real(c_double),intent(in)::f(3,3,n),w(n),g(n),eps(3,3)
    integer::q
    status=1
    if(n<1.or.n>huge(n)/9)return
    if(.not.ieee_is_finite(wp).or.wp<0d0)return
    if(.not.positive_matrix(eps,.true.))return
    do q=1,n
      if(.not.ieee_is_finite(w(q)).or..not.ieee_is_finite(g(q)))return
      if(w(q)<0d0.or.g(q)<0d0)return
      if(.not.positive_matrix(f(:,:,q),.false.))return
    enddo
    status=0
  end function

  integer(c_int) function sp_ld_tensor_permittivity(n,m,wp,eps,f,w,g,omega,result) bind(C) result(status)
    integer(c_int),value::n,m
    real(c_double),value::wp
    real(c_double),intent(in)::eps(3,3),f(3,3,n),w(n),g(n),omega(m)
    complex(c_double_complex),intent(out)::result(3,3,m)
    complex(c_double_complex)::den,value
    integer::i,j,k,q
    status=sp_ld_tensor_validate(n,wp,f,w,g,eps)
    if(status/=0)return
    status=1
    if(m<0)return
    do k=1,m
      if(.not.ieee_is_finite(omega(k)))return
      result(:,:,k)=cmplx(eps,0d0,kind=c_double)
      do q=1,n
        if(wp==0d0.or.maxval(abs(f(:,:,q)))==0d0)cycle
        if(abs(omega(k))==w(q).and.(g(q)==0d0.or.omega(k)==0d0))return
        den=cmplx(w(q)*w(q)-omega(k)*omega(k),-g(q)*omega(k),kind=c_double)
        if(.not.ieee_is_finite(real(den,kind=c_double)).or..not.ieee_is_finite(aimag(den)))return
        if(abs(den)==0d0)return
        do j=1,3
          do i=1,3
            value=result(i,j,k)+f(i,j,q)*wp**2/den
            if(.not.ieee_is_finite(real(value,kind=c_double)).or..not.ieee_is_finite(aimag(value)))return
            result(i,j,k)=value
          enddo
        enddo
      enddo
    enddo
    status=0
  end function

  integer(c_int) function sp_ld_tensor_background(m,eps,e,p,h,d,u) bind(C) result(status)
    integer(c_int),value::m
    real(c_double),intent(in)::eps(3,3),e(3,m),p(3,m),h(3,m)
    real(c_double),intent(out)::d(3,m),u(m)
    real(c_double)::v,energy,pi
    integer::i,j,k
    status=1;pi=acos(-1d0)
    if(m<0.or..not.positive_matrix(eps,.true.))return
    do k=1,m
      energy=0d0
      do i=1,3
        if(.not.ieee_is_finite(e(i,k)).or..not.ieee_is_finite(p(i,k)))return
        if(.not.ieee_is_finite(h(i,k)))return
        v=0d0
        do j=1,3
          v=v+eps(i,j)*e(j,k)
        enddo
        d(i,k)=v+4d0*pi*p(i,k)
        energy=energy+e(i,k)*v+h(i,k)*h(i,k)
        if(.not.ieee_is_finite(d(i,k)).or..not.ieee_is_finite(energy))return
      enddo
      u(k)=energy/(8d0*pi)
    enddo
    status=0
  end function

  subroutine normalized_stokes(z,s,ok)
    complex(c_double_complex),intent(in)::z(2)
    real(c_double),intent(out)::s(3)
    logical,intent(out)::ok
    complex(c_double_complex)::a,b,cross
    real(c_double)::scale,norm
    scale=max(abs(z(1)),abs(z(2)));ok=.false.
    if(.not.ieee_is_finite(scale).or.scale==0d0)return
    a=z(1)/scale;b=z(2)/scale;norm=abs(a)**2+abs(b)**2
    cross=a*conjg(b)
    s(1)=(abs(a)**2-abs(b)**2)/norm
    s(2)=2d0*real(cross,kind=c_double)/norm;s(3)=2d0*aimag(cross)/norm
    ok=.true.
  end subroutine

  integer(c_int) function sp_waveplate_diagnostic(n,time,e,omega,nx,ny,length,phasor,stokes,expected,stats) &
      bind(C) result(status)
    integer(c_int),value::n
    real(c_double),intent(in)::time(n),e(3,n)
    real(c_double),value::omega,nx,ny,length
    complex(c_double_complex),intent(out)::phasor(2)
    real(c_double),intent(out)::stokes(3),expected(3),stats(4)
    complex(c_double_complex)::phase,previous,transfer(2),den
    real(c_double)::dt,index,argument
    integer::k,j
    logical::ok
    status=1;phasor=cmplx(0d0,0d0,kind=c_double);stats=0d0
    if(n<2)return
    if(.not.ieee_is_finite(omega).or..not.ieee_is_finite(length))return
    if(.not.ieee_is_finite(nx).or..not.ieee_is_finite(ny))return
    if(omega<=0d0.or.length<=0d0.or.nx<=0d0.or.ny<=0d0)return
    do k=1,n
      if(.not.ieee_is_finite(time(k)))return
      do j=1,3
        if(.not.ieee_is_finite(e(j,k)))return
      enddo
      if(k>1)then
        dt=time(k)-time(k-1)
        if(.not.ieee_is_finite(dt).or.dt<=0d0)return
        argument=omega*time(k)
        if(.not.ieee_is_finite(argument))return
        phase=exp(cmplx(0d0,argument,kind=c_double))
        previous=exp(cmplx(0d0,omega*time(k-1),kind=c_double))
        do j=1,2
          phasor(j)=phasor(j)+.5d0*dt*(e(j,k)*phase+e(j,k-1)*previous)
        enddo
      endif
    enddo
    call normalized_stokes(phasor,stokes,ok)
    if(.not.ok)return
    do j=1,2
      index=nx
      if(j==2)index=ny
      argument=omega*index*length/137.035999084d0
      if(.not.ieee_is_finite(argument))return
      phase=exp(cmplx(0d0,argument,kind=c_double))
      den=(1d0+index)**2-(index-1d0)**2*phase**2
      if(abs(den)==0d0)return
      transfer(j)=4d0*index*phase/den
    enddo
    call normalized_stokes(transfer,expected,ok)
    if(.not.ok)return
    stats(1)=omega*(ny-nx)*length/137.035999084d0
    phase=phasor(2)*conjg(phasor(1))
    stats(2)=atan2(aimag(phase),real(phase,kind=c_double))
    phase=transfer(2)*conjg(transfer(1))
    stats(3)=atan2(aimag(phase),real(phase,kind=c_double))
    stats(4)=sqrt(sum((stokes-expected)**2))
    do j=1,4
      if(.not.ieee_is_finite(stats(j)))return
    enddo
    status=0
  end function
  ! Each column is one independent experiment, not a polarization assumption.
  integer(c_int) function sp_optical_tensor(n,omega,j,e,rtol,sigma,eps,mask) bind(C) result(status)
    integer(c_int),value::n
    real(c_double),intent(in)::omega(n)
    complex(c_double_complex),intent(in)::j(3,3,n),e(3,3,n)
    real(c_double),value::rtol
    complex(c_double_complex),intent(out)::sigma(3,3,n),eps(3,3,n)
    integer(c_int),intent(out)::mask(n)
    complex(c_double_complex)::a(3,6),row(6),pivot
    real(c_double)::scale,peak,pivsize
    integer::k,i,l,r,best
    logical::valid
    status=1
    if(n<1.or..not.ieee_is_finite(rtol))return
    if(rtol<0d0.or.rtol>=1d0)return
    peak=0d0
    do k=1,n
      if(.not.ieee_is_finite(omega(k)))return
      do l=1,3
        do i=1,3
          if(.not.ieee_is_finite(real(j(i,l,k),c_double)))return
          if(.not.ieee_is_finite(aimag(j(i,l,k))))return
          if(.not.ieee_is_finite(real(e(i,l,k),c_double)))return
          if(.not.ieee_is_finite(aimag(e(i,l,k))))return
          peak=max(peak,abs(e(i,l,k)))
        enddo
      enddo
    enddo
    ! Invalid bins are masked. Python supplies NaNs after this native operation.
    sigma=cmplx(0d0,0d0,c_double);eps=sigma;mask=0
    do k=1,n
      if(omega(k)<=0d0)cycle
      scale=maxval(abs(e(:,:,k)))
      if(scale<=rtol*peak.or.scale==0d0)cycle
      a(:,1:3)=e(:,:,k)/scale;a(:,4:6)=cmplx(0d0,0d0,c_double)
      do i=1,3
        a(i,i+3)=cmplx(1d0,0d0,c_double)
      enddo
      valid=.true.
      do l=1,3
        best=l;pivsize=abs(a(l,l))
        do r=l+1,3
          if(abs(a(r,l))>pivsize)then
            best=r;pivsize=abs(a(r,l))
          endif
        enddo
        if(pivsize<=max(rtol,64d0*epsilon(1d0)))then
          valid=.false.;exit
        endif
        row=a(l,:);a(l,:)=a(best,:);a(best,:)=row
        pivot=a(l,l);a(l,:)=a(l,:)/pivot
        do r=1,3
          if(r==l)cycle
          pivot=a(r,l);a(r,:)=a(r,:)-pivot*a(l,:)
        enddo
      enddo
      if(.not.valid)cycle
      sigma(:,:,k)=matmul(j(:,:,k)/scale,a(:,4:6))
      eps(:,:,k)=cmplx(0d0,4d0*acos(-1d0)/omega(k),c_double)*sigma(:,:,k)
      do i=1,3
        eps(i,i,k)=eps(i,i,k)+1d0
      enddo
      do l=1,3
        do i=1,3
          if(.not.ieee_is_finite(real(eps(i,l,k),c_double)))valid=.false.
          if(.not.ieee_is_finite(aimag(eps(i,l,k))))valid=.false.
          if(.not.ieee_is_finite(real(sigma(i,l,k),c_double)))valid=.false.
          if(.not.ieee_is_finite(aimag(sigma(i,l,k))))valid=.false.
        enddo
      enddo
      if(valid)mask(k)=1
    enddo
    status=0
  end function
  integer(c_int) function sp_ld_tensor_fit(n,p,omega,target,eps,par,residual,f,w,g) bind(C) result(status)
    integer(c_int),value::n,p
    real(c_double),intent(in)::omega(n),eps(3,3),par(8*p)
    complex(c_double_complex),intent(in)::target(3,3,n)
    real(c_double),intent(out)::residual(18*n),f(3,3,p),w(p),g(p)
    real(c_double)::l(3,3),weight
    complex(c_double_complex)::z(3,3),den,diff
    integer::q,k,i,j,b,t
    status=1
    if(n<1.or.p<1)return
    do i=1,8*p
      if(.not.ieee_is_finite(par(i)))return
    enddo
    do q=1,p
      b=8*(q-1);l=0d0
      l(1,1)=par(b+1);l(2,1)=par(b+2);l(2,2)=par(b+3)
      l(3,1)=par(b+4);l(3,2)=par(b+5);l(3,3)=par(b+6)
      f(:,:,q)=matmul(l,transpose(l));w(q)=par(b+7);g(q)=par(b+8)
      if(w(q)<0d0.or.g(q)<0d0)return
    enddo
    do k=1,n
      if(.not.ieee_is_finite(omega(k)).or.omega(k)<=0d0)return
      z=cmplx(eps,0d0,c_double);weight=0d0
      do j=1,3
        do i=1,3
          if(.not.ieee_is_finite(real(target(i,j,k),c_double)))return
          if(.not.ieee_is_finite(aimag(target(i,j,k))))return
          weight=weight+abs(target(i,j,k))**2
        enddo
      enddo
      weight=1d0/max(1d0,sqrt(weight))
      do q=1,p
        den=cmplx(w(q)**2-omega(k)**2,-g(q)*omega(k),c_double)
        if(abs(den)==0d0)return
        z=z+f(:,:,q)/den
      enddo
      do j=1,3
        do i=1,3
          t=18*(k-1)+2*(3*(j-1)+i)-1
          diff=(z(i,j)-target(i,j,k))*weight
          residual(t)=real(diff,c_double);residual(t+1)=aimag(diff)
          if(.not.ieee_is_finite(residual(t)))return
          if(.not.ieee_is_finite(residual(t+1)))return
        enddo
      enddo
    enddo
    status=0
  end function
  integer(c_int) function sp_waveplate_rotated(n,time,e,omega,nx,ny,length,angle,phasor,stokes,expected,stats) &
      bind(C) result(status)
    integer(c_int),value::n
    real(c_double),intent(in)::time(n),e(3,n)
    real(c_double),value::omega,nx,ny,length,angle
    complex(c_double_complex),intent(out)::phasor(2)
    real(c_double),intent(out)::stokes(3),expected(3),stats(4)
    real(c_double),allocatable::principal(:,:)
    real(c_double)::cs,sn,reference(3)
    complex(c_double_complex)::z(2)
    logical::ok
    integer::k
    status=1
    if(n<2.or..not.ieee_is_finite(angle))return
    cs=cos(angle);sn=sin(angle);allocate(principal(3,n))
    do k=1,n
      principal(1,k)=cs*e(1,k)+sn*e(2,k)
      principal(2,k)=-sn*e(1,k)+cs*e(2,k)
      principal(3,k)=e(3,k)
    enddo
    status=sp_waveplate_diagnostic(n,time,principal,omega,nx,ny,length,z,stokes,reference,stats)
    if(status/=0)return
    phasor(1)=cs*z(1)-sn*z(2);phasor(2)=sn*z(1)+cs*z(2)
    call normalized_stokes(phasor,stokes,ok)
    if(.not.ok)then
      status=1;return
    endif
    expected(1)=cos(2d0*angle)*reference(1)-sin(2d0*angle)*reference(2)
    expected(2)=sin(2d0*angle)*reference(1)+cos(2d0*angle)*reference(2)
    expected(3)=reference(3)
    stats(4)=sqrt(sum((stokes-expected)**2))
  end function
end module
