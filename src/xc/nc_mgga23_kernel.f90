! Noncollinear BR exchange and corrected CS correlation (Hartree atomic units).
! N,T,J=(x0 I+x.sigma)/2; T=sum_i grad(psi_i) grad(psi_i)^dagger (NO 1/2).
! PRB 107,165111 (2023), Eqs. (5)-(8), App. B; PRB 110,199902 (2024), Eq. (6).
! Forward automatic differentiation includes the implicit BR root derivative.
module nc_mgga23_kernel
  implicit none
  private
  public :: nc_mgga23_eval
  integer,parameter :: nv=36
  type :: dual
    real(8) :: v=0d0
    real(8) :: d(nv)=0d0
  end type
  interface operator(+)
    module procedure adddd,addrd,adddr
  end interface
  interface operator(-)
    module procedure subdd,subrd,subdr,negate
  end interface
  interface operator(*)
    module procedure muldd,mulrd,muldr
  end interface
  interface operator(/)
    module procedure divdd,divrd,divdr
  end interface
contains
  pure function negate(a) result(c)
    implicit none
    type(dual),intent(in)::a
    type(dual)::c
    c%v=-a%v; c%d=-a%d
  end function
  subroutine nc_mgga23_eval(n4,grad4,lap4,tau4,j4,e,dn,dgrad,dlap,dtau,dj,ierr,ex,ec)
    use,intrinsic :: ieee_arithmetic,only: ieee_is_finite
    implicit none
    real(8),intent(in):: n4(4),grad4(4,3),lap4(4),tau4(4),j4(4,3)
    real(8),intent(out):: e,dn(4),dgrad(4,3),dlap(4),dtau(4),dj(4,3)
    integer,intent(out):: ierr
    real(8),intent(out),optional:: ex,ec
    real(8):: raw(nv),pi,ca,cb,cc,cd,lo,hi,x,mid,r,fp,ev
    type(dual):: a(nv),n(4),g(4,3),l(4),t(4),j(4,3),m2,gn2,gm2,j2,nt,ml
    type(dual):: tw,bt,zx,zc,h,q,ratio,xx,rs,fx,exd,ecd,tot,bracket
    integer:: i,k,it
    e=0d0; dn=0d0; dgrad=0d0; dlap=0d0; dtau=0d0; dj=0d0; ierr=0
    if(present(ex)) ex=0d0
    if(present(ec)) ec=0d0
    raw(1:4)=n4; raw(5:16)=reshape(grad4,[12]); raw(17:20)=lap4
    raw(21:24)=tau4; raw(25:36)=reshape(j4,[12])
    do i=1,nv
      if(.not.ieee_is_finite(raw(i))) then
        ierr=1
        return
      endif
    enddo
    if(n4(1)<0d0.or.sqrt(sum(n4(2:4)**2))>n4(1)*(1d0+1d-10)+1d-25) then
      ierr=2
      return
    endif
    ! Explicit vacuum convention; no positive-density clipping or alpha clipping.
    if(n4(1)<=1d-25) return
    do i=1,nv
      a(i)%v=raw(i); a(i)%d(i)=1d0
    enddo
    n=a(1:4); g=reshape(a(5:16),[4,3]); l=a(17:20)
    t=a(21:24); j=reshape(a(25:36),[4,3])
    m2=constant(0d0); gn2=constant(0d0); gm2=constant(0d0)
    j2=constant(0d0); nt=n(1)*t(1); ml=constant(0d0)
    do i=2,4
      m2=m2+n(i)*n(i); nt=nt+n(i)*t(i); ml=ml+n(i)*l(i)
    enddo
    do k=1,3
      gn2=gn2+g(1,k)*g(1,k)
      do i=1,4
        j2=j2+j(i,k)*j(i,k)
        if(i>1) gm2=gm2+g(i,k)*g(i,k)
      enddo
    enddo
    ! Eq. (6) in Pauli trace coordinates; commutator trace is identically zero.
    bt=(0.5d0*(nt-j2)+0.25d0*(n(1)*l(1)-ml)+0.125d0*(gn2-gm2))/n(1)
    tw=gn2/(4d0*n(1))
    zx=0.5d0*(1d0+m2/(n(1)*n(1)))
    zc=1d0-m2/(n(1)*n(1))
    h=n(1)*zx
    q=(l(1)-1.6d0*(bt-tw))/6d0
    pi=acos(-1d0)
    ratio=6d0*q/(4d0*pi**(2d0/3d0)*power(h,5d0/3d0))
    r=ratio%v
    if(.not.ieee_is_finite(r)) then
      ierr=3
      return
    endif
    ! Monotone equation (x-2) exp(2x/3)/x=ratio, continuous through Q=0, x=2.
    lo=1d-150; hi=4d0
    do while(br_value(hi)<r)
      hi=hi*2d0
      if(hi>512d0) then
        ierr=3
        return
      endif
    enddo
    do it=1,600
      mid=0.5d0*(lo+hi)
      if(br_value(mid)>r) then
        hi=mid
      else
        lo=mid
      endif
      if(hi-lo<2d-15*max(mid,1d-150)) exit
    enddo
    x=0.5d0*(lo+hi)
    ev=exp(2d0*x/3d0)
    fp=ev*(2d0/3d0-4d0/(3d0*x)+2d0/(x*x))
    xx%v=x; xx%d=ratio%d/fp
    ! Stable [1-exp(-x)(1+x/2)]/x, avoiding cancellation as x -> 0.
    if(x<1d-3) then
      fx=0.5d0-xx*xx/12d0+power(xx,3d0)/24d0-power(xx,4d0)/80d0 &
         +power(xx,5d0)/360d0-power(xx,6d0)/2016d0
    else
      fx=(1d0-exponential(-1d0*xx)*(1d0+0.5d0*xx))/xx
    endif
    exd=-1d0*n(1)*power(pi*h,1d0/3d0)*exponential(xx/3d0)*fx
    ca=.04918d0; cb=.132d0*(4d0*pi/3d0)**(8d0/3d0)
    cc=.2533d0*(4d0*pi/3d0)**(1d0/3d0); cd=.349d0*(4d0*pi/3d0)**(1d0/3d0)
    rs=power((4d0*pi/3d0)*n(1),-1d0/3d0)
    ! Erratum Eq. (6): both the 2*bar_tau coefficient and Laplacian term matter.
    bracket=2d0*bt+gm2/(4d0*n(1))-3d0*tw+(ml-n(1)*l(1))/(4d0*n(1))
    ecd=-ca*n(1)*zc*(1d0+0.5d0*cb*power(rs,5d0)*bracket*exponential(-cc*rs))/(1d0+cd*rs)
    tot=exd+ecd
    if(.not.ieee_is_finite(tot%v)) then
      ierr=4
      return
    endif
    do i=1,nv
      if(.not.ieee_is_finite(tot%d(i))) then
        ierr=4
        return
      endif
    enddo
    e=tot%v; dn=tot%d(1:4); dgrad=reshape(tot%d(5:16),[4,3]); dlap=tot%d(17:20)
    dtau=tot%d(21:24); dj=reshape(tot%d(25:36),[4,3])
    if(present(ex)) ex=exd%v
    if(present(ec)) ec=ecd%v
  end subroutine
  pure real(8) function br_value(x) result(y)
    implicit none
    real(8),intent(in)::x
    y=(1d0-2d0/x)*exp(2d0*x/3d0)
  end function
  pure function constant(v) result(a)
    implicit none
    real(8),intent(in)::v
    type(dual)::a
    a%v=v; a%d=0d0
  end function
  pure function power(a,p) result(c)
    implicit none
    type(dual),intent(in)::a
    real(8),intent(in)::p
    type(dual)::c
    c%v=a%v**p; c%d=p*a%v**(p-1d0)*a%d
  end function
  pure function exponential(a) result(c)
    implicit none
    type(dual),intent(in)::a
    type(dual)::c
    c%v=exp(a%v); c%d=c%v*a%d
  end function
  pure function adddd(a,b) result(c)
    implicit none
    type(dual),intent(in)::a,b
    type(dual)::c
    c%v=a%v+b%v; c%d=a%d+b%d
  end function
  pure function addrd(a,b) result(c)
    implicit none
    real(8),intent(in)::a
    type(dual),intent(in)::b
    type(dual)::c
    c=adddd(constant(a),b)
  end function
  pure function adddr(a,b) result(c)
    implicit none
    type(dual),intent(in)::a
    real(8),intent(in)::b
    type(dual)::c
    c=adddd(a,constant(b))
  end function
  pure function subdd(a,b) result(c)
    implicit none
    type(dual),intent(in)::a,b
    type(dual)::c
    c%v=a%v-b%v; c%d=a%d-b%d
  end function
  pure function subrd(a,b) result(c)
    implicit none
    real(8),intent(in)::a
    type(dual),intent(in)::b
    type(dual)::c
    c=subdd(constant(a),b)
  end function
  pure function subdr(a,b) result(c)
    implicit none
    type(dual),intent(in)::a
    real(8),intent(in)::b
    type(dual)::c
    c=subdd(a,constant(b))
  end function
  pure function muldd(a,b) result(c)
    implicit none
    type(dual),intent(in)::a,b
    type(dual)::c
    c%v=a%v*b%v; c%d=a%d*b%v+a%v*b%d
  end function
  pure function mulrd(a,b) result(c)
    implicit none
    real(8),intent(in)::a
    type(dual),intent(in)::b
    type(dual)::c
    c=muldd(constant(a),b)
  end function
  pure function muldr(a,b) result(c)
    implicit none
    type(dual),intent(in)::a
    real(8),intent(in)::b
    type(dual)::c
    c=muldd(a,constant(b))
  end function
  pure function divdd(a,b) result(c)
    implicit none
    type(dual),intent(in)::a,b
    type(dual)::c
    c%v=a%v/b%v; c%d=(a%d-c%v*b%d)/b%v
  end function
  pure function divrd(a,b) result(c)
    implicit none
    real(8),intent(in)::a
    type(dual),intent(in)::b
    type(dual)::c
    c=divdd(constant(a),b)
  end function
  pure function divdr(a,b) result(c)
    implicit none
    type(dual),intent(in)::a
    real(8),intent(in)::b
    type(dual)::c
    c=divdd(a,constant(b))
  end function
end module nc_mgga23_kernel
