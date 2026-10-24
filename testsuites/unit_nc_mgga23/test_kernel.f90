program test_kernel
  use nc_mgga23_kernel, only: nc_mgga23_eval
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
  implicit none
  real(8) :: scratch(36),a(36),d(36),ep,em,e,ex,ec,h,err,n(4),g(4,3),l(4),t(4),j(4,3),pi,r,expected,rs,bt,tw,zc,gm2,ml
  real(8) :: dn(4),dg(4,3),dl(4),dt(4),dj(4,3),e0,shift(3)
  complex(8):: psi(2,2),gp(2,2,3),lp(2,2),zspin(2,2)
  integer :: k,ierr,itest
  a=0d0
  a(1:4)=[1.2d0,.2d0,-.1d0,.3d0]
  a(5:16)=[.1d0,.03d0,-.02d0,.01d0,-.06d0,.04d0,.02d0,-.03d0,.03d0,-.01d0,.02d0,.04d0]
  a(17:20)=[-.15d0,.03d0,.02d0,-.01d0]
  a(21:24)=[2d0,.1d0,-.02d0,.15d0]
  a(25:36)=.01d0
  call evaluate(a,e,d,ierr)
  if(ierr/=0) stop 1
  do k=1,36
    h=1d-6
    a(k)=a(k)+h
    call evaluate(a,ep,scratch,ierr)
    a(k)=a(k)-2*h
    call evaluate(a,em,scratch,ierr)
    a(k)=a(k)+h
    err=abs((ep-em)/(2*h)-d(k))
    if(err>2d-7) then
      print *, 'derivative failed', k,err,d(k),(ep-em)/(2*h)
      stop 2
    endif
  enddo
  e0=e
  ! A proper global spin rotation: x -> y, y -> z, z -> x.
  call rotate(a(1:4)); call rotate(a(17:20)); call rotate(a(21:24))
  do k=1,3
    call rotate(a(4+4*(k-1)+1:4+4*k))
    call rotate(a(24+4*(k-1)+1:24+4*k))
  enddo
  call evaluate(a,e,d,ierr)
  if(abs(e-e0)>1d-13.or.ierr/=0) stop 3
  do k=1,36
    h=a(k)
    a(k)=ieee_value(0d0,ieee_quiet_nan)
    call evaluate(a,e,d,ierr)
    if(ierr==0) stop 4
    a(k)=ieee_value(0d0,ieee_positive_inf)
    call evaluate(a,e,d,ierr)
    if(ierr==0) stop 5
    a(k)=h
  enddo
  a=0d0
  call evaluate(a,e,d,ierr)
  if(ierr/=0.or.e/=0d0.or.maxval(abs(d))/=0d0) stop 6
  ! Exact hydrogenic exchange-hole model, fully polarized real 1s orbital.
  pi=acos(-1d0)
  do k=1,3
    r=.5d0*k; h=exp(-2*r)/pi
    n=[h,0d0,0d0,h]; g=0d0; g(:,1)=-2d0*n
    l=(4d0-4d0/r)*n; t=n; j=0d0
    call nc_mgga23_eval(n,g,l,t,j,e,dn,dg,dl,dt,dj,ierr,ex,ec)
    expected=-h/(2d0*r)*(1d0-exp(-2d0*r)*(1d0+r))
    if(ierr/=0.or.abs(ex-expected)>1d-13.or.abs(ec)>1d-14) stop 7
  enddo
  ! Curvatures on both sides and exactly at Q=0; all partials remain smooth.
  do k=-1,1
    a=0d0; a(1)=1d0; a(21)=2d0; a(17)=8d0/3d0+real(k,8)
    call evaluate(a,e,d,ierr)
    if(ierr/=0) stop 8
    do itest=1,36
      h=1d-6
      a(itest)=a(itest)+h; call evaluate(a,ep,scratch,ierr)
      a(itest)=a(itest)-2*h; call evaluate(a,em,scratch,ierr)
      a(itest)=a(itest)+h
      if(abs((ep-em)/(2*h)-d(itest))>3d-7) stop 9
    enddo
  enddo
  ! Independent explicit evaluation of corrected correlation, mixed polarization.
  n=[1.2d0,.2d0,-.1d0,.3d0]; t=[2d0,.1d0,.05d0,.15d0]
  g=.02d0; l=[-.1d0,.02d0,-.03d0,.04d0]; j=.01d0
  call nc_mgga23_eval(n,g,l,t,j,e,dn,dg,dl,dt,dj,ierr,ex,ec)
  gm2=sum(g(2:4,:)**2); ml=sum(n(2:4)*l(2:4)); tw=sum(g(1,:)**2)/(4*n(1))
  bt=(sum(n*t)/2-sum(j*j)/2+(n(1)*l(1)-ml)/4+(sum(g(1,:)**2)-gm2)/8)/n(1)
  zc=1-sum(n(2:4)**2)/n(1)**2; rs=(4*pi*n(1)/3)**(-1d0/3d0)
  expected=-.04918d0*n(1)*zc/(1+.349d0*(4*pi/3)**(1d0/3d0)*rs) &
    *(1+.132d0*(4*pi/3)**(8d0/3d0)*rs**5/2 &
    *(2*bt+gm2/(4*n(1))-3*tw+(ml-n(1)*l(1))/(4*n(1))) &
    *exp(-.2533d0*(4*pi/3)**(1d0/3d0)*rs))
  if(ierr/=0.or.abs(ec-expected)>1d-14) stop 10
  ! Unaffected original Eqs. (14),(15), independently contracted matrices.
  bt=(sum(n*t)/2-sum(j*j)/2-sum(n*l)/8)/n(1) ! tau_HF, not bar_tau
  expected=-.04918d0*n(1)*zc/(1+.349d0*(4*pi/3)**(1d0/3d0)*rs) &
    *(1+.132d0*(4*pi/3)**(8d0/3d0)*rs**5*(bt-tw+l(1)/4) &
    *exp(-.2533d0*(4*pi/3)**(1d0/3d0)*rs))
  if(abs(ec-expected)>1d-14) stop 14
  ! Common local U(1) phase gradient: J'=J+A*N, T'=T+2 A.J+A^2 N.
  e0=e; shift=[.17d0,-.23d0,.1d0]
  do k=1,3
    t=t+2*shift(k)*j(:,k)+shift(k)**2*n
    j(:,k)=j(:,k)+shift(k)*n
  enddo
  call nc_mgga23_eval(n,g,l,t,j,e,dn,dg,dl,dt,dj,ierr)
  if(ierr/=0.or.abs(e-e0)>1d-13) stop 11
  ! Local SU(2) rotation U(r)=exp(i q.r sigma_z), U(0)=I, including Laplacian terms.
  psi=reshape([cmplx(.7d0,.1d0,8),cmplx(.2d0,-.3d0,8), &
    cmplx(.1d0,.4d0,8),cmplx(.6d0,.2d0,8)],[2,2])
  do k=1,3
    gp(:,:,k)=cmplx(.1d0*k,.03d0,8)*psi+cmplx(.02d0,-.01d0*k,8)
  enddo
  lp=cmplx(-.2d0,.1d0,8)*psi
  call orbital_fields(psi,gp,lp,n,g,l,t,j)
  call nc_mgga23_eval(n,g,l,t,j,e,dn,dg,dl,dt,dj,ierr)
  if(ierr/=0) stop 12
  e0=e
  zspin=psi; zspin(2,:)=-zspin(2,:)
  do k=1,3
    lp(1,:)=lp(1,:)+cmplx(0d0,2*shift(k),8)*gp(1,:,k)
    lp(2,:)=lp(2,:)-cmplx(0d0,2*shift(k),8)*gp(2,:,k)
    lp=lp-shift(k)**2*psi
    gp(:,:,k)=gp(:,:,k)+cmplx(0d0,shift(k),8)*zspin
  enddo
  call orbital_fields(psi,gp,lp,n,g,l,t,j)
  call nc_mgga23_eval(n,g,l,t,j,e,dn,dg,dl,dt,dj,ierr)
  if(ierr/=0.or.abs(e-e0)>1d-12) then
    print *, 'SU(2) failed',e,e0,ierr
    stop 13
  endif
  print *, 'PASS: partials, rotation, NaN/Inf, vacuum, hydrogen, Q=0, corrected CS, U(1), SU(2)'
contains
  subroutine orbital_fields(p,gp,lp,n,g,l,t,j)
    implicit none
    complex(8),intent(in)::p(2,2),gp(2,2,3),lp(2,2)
    real(8),intent(out)::n(4),g(4,3),l(4),t(4),j(4,3)
    complex(8)::nn(2,2),gg(2,2),ll(2,2),tt(2,2),jj(2,2)
    integer::k
    nn=matmul(p,conjg(transpose(p)))
    ll=matmul(lp,conjg(transpose(p)))+matmul(p,conjg(transpose(lp)))
    tt=0d0
    do k=1,3
      gg=matmul(gp(:,:,k),conjg(transpose(p)))+matmul(p,conjg(transpose(gp(:,:,k))))
      jj=(matmul(gp(:,:,k),conjg(transpose(p)))-matmul(p,conjg(transpose(gp(:,:,k)))))/cmplx(0d0,2d0,8)
      call pauli(gg,g(:,k)); call pauli(jj,j(:,k))
      tt=tt+matmul(gp(:,:,k),conjg(transpose(gp(:,:,k))))
    enddo
    ll=ll+2*tt
    call pauli(nn,n); call pauli(ll,l); call pauli(tt,t)
  end subroutine
  subroutine pauli(v,a)
    implicit none
    complex(8),intent(in)::v(2,2)
    real(8),intent(out)::a(4)
    a=[real(v(1,1)+v(2,2),8),2*real(v(1,2),8),-2*aimag(v(1,2)),real(v(1,1)-v(2,2),8)]
  end subroutine
  subroutine rotate(v)
    implicit none
    real(8),intent(inout)::v(4)
    v(2:4)=[v(4),v(2),v(3)]
  end subroutine
  subroutine evaluate(v,e,d,ierr)
    implicit none
    real(8),intent(in)::v(36)
    real(8),intent(out)::e,d(36)
    integer,intent(out)::ierr
    real(8)::dn(4),dg(4,3),dl(4),dt(4),dj(4,3)
    call nc_mgga23_eval(v(1:4),reshape(v(5:16),[4,3]),v(17:20),v(21:24), &
      reshape(v(25:36),[4,3]),e,dn,dg,dl,dt,dj,ierr)
    d(1:4)=dn; d(5:16)=reshape(dg,[12]); d(17:20)=dl; d(21:24)=dt; d(25:36)=reshape(dj,[12])
  end subroutine
end program
