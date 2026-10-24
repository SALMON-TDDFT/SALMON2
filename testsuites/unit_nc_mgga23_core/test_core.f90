program test_core
  use nc_mgga23_core,only: core23_fields,core23_energy_deda
  use nc_mgga23_kernel,only: nc_mgga23_eval
  use,intrinsic:: ieee_arithmetic,only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  real(8)::n(4),g(4,3),l(4),t(4),j(4,3),nc(4),tc(4),jc(4,3),ac(3),jref(4,3),deda(3)
  real(8)::dn(4),dg(4,3),dl(4),dt(4),dj(4,3),e,e0,ep,em,h,rho,tau,shift(3),r(17),base(17)
  real(8)::nv(4),gv(4,3),lv(4),tv(4),jv(4,3),gc(4,3),lc(4),ncheck(4),gcheck(4,3),lcheck(4), &
    tcheck(4),jcheck(4,3)
  complex(8)::p(2,2),dp(2,2,3),lp(2,2),cp(2,2),cdp(2,2,3),clp(2,2),zspin(2,2),ci
  integer::k,ierr,i
  ci=cmplx(0d0,1d0,8)
  cp=0d0;cp(1,1)=sqrt(.2d0);cp(2,2)=cp(1,1)
  cdp=0d0
  do k=1,3
    cdp(1,1,k)=.03d0*k;cdp(2,2,k)=.03d0*k
  enddo
  ! Spin-current-carrying Kramers-paired reference, with scalar N and T.
  cdp(1,2,1)=ci*.04d0;cdp(2,1,1)=ci*.04d0
  cdp(1,2,2)=.02d0;cdp(2,1,2)=-.02d0
  cdp(1,1,3)=cdp(1,1,3)+ci*.03d0;cdp(2,2,3)=cdp(2,2,3)-ci*.03d0
  clp=-.08d0*cp
  call orbital_fields(cp,cdp,clp,nc,gc,lc,tc,jref)
  rho=nc(1);tau=tc(1)
  if(maxval(abs(nc(2:4)))>1d-15.or.maxval(abs(tc(2:4)))>1d-15) stop 1
  ac=[.11d0,-.07d0,.03d0]
  call core23_fields(rho,tau,jref,ac,nc,tc,jc,ierr)
  if(ierr/=0) stop 2
  do k=1,3
    clp=clp-2d0*ci*ac(k)*cdp(:,:,k)-ac(k)**2*cp
    cdp(:,:,k)=cdp(:,:,k)-ci*ac(k)*cp
  enddo
  call orbital_fields(cp,cdp,clp,ncheck,gcheck,lcheck,tcheck,jcheck)
  if(maxval(abs(nc-ncheck))+maxval(abs(tc-tcheck))+maxval(abs(jc-jcheck))>1d-14) stop 3
  if(maxval(abs(gc-gcheck))+maxval(abs(lc-lcheck))>1d-14) stop 4
  p=reshape([cmplx(.7d0,.1d0,8),cmplx(.2d0,-.3d0,8), &
    cmplx(.1d0,.4d0,8),cmplx(.6d0,.2d0,8)],[2,2])
  do k=1,3
    dp(:,:,k)=cmplx(.1d0*k,.03d0,8)*p+cmplx(.02d0,-.01d0*k,8)
  enddo
  lp=cmplx(-.2d0,.1d0,8)*p
  call orbital_fields(p,dp,lp,nv,gv,lv,tv,jv)
  n=nv+nc;g=gv+gc;l=lv+lc;t=tv+tc;j=jv+jc
  call nc_mgga23_eval(n,g,l,t,j,e,dn,dg,dl,dt,dj,ierr)
  if(ierr/=0) stop 5
  call core23_energy_deda(nc,jref,ac,dt,dj,deda)
  e0=e
  h=1d-6
  do k=1,3
    ac(k)=ac(k)+h
    call core23_fields(rho,tau,jref,ac,nc,tc,jc,ierr)
    call nc_mgga23_eval(n,g,l,tv+tc,jv+jc,ep,dn,dg,dl,dt,dj,ierr)
    ac(k)=ac(k)-2*h
    call core23_fields(rho,tau,jref,ac,nc,tc,jc,ierr)
    call nc_mgga23_eval(n,g,l,tv+tc,jv+jc,em,dn,dg,dl,dt,dj,ierr)
    ac(k)=ac(k)+h
    if(abs((ep-em)/(2*h)-deda(k))>2d-8) stop 6
  enddo
  ! Co-transform BOTH core and valence spinors with a local SU(2) rotation.
  shift=[.17d0,-.23d0,.1d0]
  zspin=p;zspin(2,:)=-zspin(2,:)
  do k=1,3
    lp(1,:)=lp(1,:)+2*ci*shift(k)*dp(1,:,k)
    lp(2,:)=lp(2,:)-2*ci*shift(k)*dp(2,:,k)
    lp=lp-shift(k)**2*p
    dp(:,:,k)=dp(:,:,k)+ci*shift(k)*zspin
  enddo
  zspin=cp;zspin(2,:)=-zspin(2,:)
  do k=1,3
    clp(1,:)=clp(1,:)+2*ci*shift(k)*cdp(1,:,k)
    clp(2,:)=clp(2,:)-2*ci*shift(k)*cdp(2,:,k)
    clp=clp-shift(k)**2*cp
    cdp(:,:,k)=cdp(:,:,k)+ci*shift(k)*zspin
  enddo
  call orbital_fields(p,dp,lp,nv,gv,lv,tv,jv)
  call orbital_fields(cp,cdp,clp,nc,gc,lc,tc,jc)
  call nc_mgga23_eval(nv+nc,gv+gc,lv+lc,tv+tc,jv+jc,e,dn,dg,dl,dt,dj,ierr)
  if(ierr/=0.or.abs(e-e0)>1d-12) stop 7
  ! Common U(1) phase gradient; inhomogeneous amplitudes retained.
  do k=1,3
    lp=lp+2*ci*shift(k)*dp(:,:,k)-shift(k)**2*p
    dp(:,:,k)=dp(:,:,k)+ci*shift(k)*p
    clp=clp+2*ci*shift(k)*cdp(:,:,k)-shift(k)**2*cp
    cdp(:,:,k)=cdp(:,:,k)+ci*shift(k)*cp
  enddo
  call orbital_fields(p,dp,lp,nv,gv,lv,tv,jv)
  call orbital_fields(cp,cdp,clp,nc,gc,lc,tc,jc)
  call nc_mgga23_eval(nv+nc,gv+gc,lv+lc,tv+tc,jv+jc,e,dn,dg,dl,dt,dj,ierr)
  if(ierr/=0.or.abs(e-e0)>1d-12) stop 8
  jref=0d0
  call core23_fields(0d0,0d0,jref,ac,nc,tc,jc,ierr)
  if(ierr/=0.or.maxval(abs(nc))+maxval(abs(tc))+maxval(abs(jc))>0d0) stop 9
  base=0d0;base(1:2)=[.4d0,.3d0];base(15:17)=ac
  do i=1,17
    r=base;r(i)=ieee_value(0d0,ieee_quiet_nan)
    call core23_fields(r(1),r(2),reshape(r(3:14),[4,3]),r(15:17),nc,tc,jc,ierr)
    if(ierr==0) stop 10
    r=base;r(i)=ieee_value(0d0,ieee_positive_inf)
    call core23_fields(r(1),r(2),reshape(r(3:14),[4,3]),r(15:17),nc,tc,jc,ierr)
    if(ierr==0) stop 11
  enddo
  print *, 'PASS core phase response, dE/dAc finite differences, core+valence U1/SU2 covariance, vacuum, IEEE'
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
end program
