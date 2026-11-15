! Fixed-temperature, core-weighted DC occupations and their fixed-charge response.
! Temperature denotes kBT (Hartree); weights include core norms and k weights.
module dc_thermal
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer,parameter,public :: dc_thermal_invalid=1,dc_thermal_capacity=2, &
    dc_thermal_unconverged=3,dc_thermal_singular=4
  public :: solve_dc_thermal,response_dc_thermal,fermi_entropy
  ! Rank-specific loops keep IEEE inquiries scalar on Fujitsu compilers.
  private :: salmon_all_finite,finite_real_1d
  interface salmon_all_finite
    module procedure finite_real_1d
  end interface
contains
  pure real(8) function fermi_entropy(f) result(s)
    real(8),intent(in) :: f
    s=0d0
    if(f>0d0.and.f<1d0)s=-f*log(f)-(1d0-f)*log(1d0-f)
  end function

  pure logical function valid_inputs(e,w,t,g) result(valid)
    real(8),intent(in) :: e(:),w(:),t,g
    valid=.false.
    if(size(e)==0.or.size(w)/=size(e))return
    if(.not.salmon_all_finite(e).or..not.salmon_all_finite(w))return
    if(.not.ieee_is_finite(t).or..not.ieee_is_finite(g))return
    if(t<=0d0.or.g<=0d0.or.any(w<0d0))return
    valid=sum(w)>0d0
  end function

  pure subroutine evaluate(e,t,mu,f,b,s)
    real(8),intent(in) :: e(:),t,mu
    real(8),intent(out) :: f(:),b(:),s(:)
    real(8) :: x,z,l
    integer :: i
    do i=1,size(e)
      x=(e(i)-mu)/t
      z=exp(-abs(x))
      if(x>=0d0)then
        f(i)=z/(1d0+z)
      else
        f(i)=1d0/(1d0+z)
      endif
      ! Keep minority tails even when the majority occupation rounds to one.
      b(i)=z/(1d0+z)**2
      if(z<1d-8)then
        l=z*(1d0-z*.5d0)
      else
        l=log(1d0+z)
      endif
      s(i)=l+abs(x)*z/(1d0+z)
    enddo
  end subroutine

  pure subroutine solve_dc_thermal(e,w,t,g,n,mu,f,ts,status)
    real(8),intent(in) :: e(:),w(:),t,g,n
    real(8),intent(out) :: mu,f(:),ts
    integer,intent(out) :: status
    real(8) :: lo,hi,charge,capacity,tol,endpoint_tol,b(size(e)),s(size(e))
    integer :: iter
    status=dc_thermal_invalid;mu=0d0;f=0d0;ts=0d0
    if(.not.valid_inputs(e,w,t,g).or.size(f)/=size(e))return
    if(.not.ieee_is_finite(n))return
    capacity=g*sum(w)
    ! Endpoint allowance is only for normalization roundoff in occupied-only runs.
    endpoint_tol=1d-10+64d0*epsilon(n)*capacity
    status=dc_thermal_capacity
    if(n<0d0.or.n>capacity+endpoint_tol)return
    lo=minval(e)-40d0*t;hi=maxval(e)+40d0*t
    if(n==0d0.or.n>=capacity)then
      mu=lo
      if(n>0d0)mu=hi
      call evaluate(e,t,mu,f,b,s)
      ts=g*t*sum(w*s);status=0
      return
    endif
    tol=32d0*epsilon(n)*max(1d0,capacity)
    status=dc_thermal_unconverged
    do iter=1,256
      mu=.5d0*lo+.5d0*hi
      call evaluate(e,t,mu,f,b,s)
      charge=g*sum(w*f)
      if(abs(charge-n)<=tol)then
        status=0
        exit
      endif
      if(mu==lo.or.mu==hi)exit
      if(charge<n)then
        lo=mu
      else
        hi=mu
      endif
    enddo
    ts=g*t*sum(w*s)
  end subroutine

  pure subroutine response_dc_thermal(e,w,t,g,mu,de,dw,dmu,df,dts,status)
    real(8),intent(in) :: e(:),w(:),t,g,mu,de(:),dw(:)
    real(8),intent(out) :: dmu,df(:),dts
    integer,intent(out) :: status
    real(8) :: f(size(e)),b(size(e)),s(size(e)),susceptibility
    status=dc_thermal_invalid;dmu=0d0;df=0d0;dts=0d0
    if(.not.valid_inputs(e,w,t,g))return
    if(size(de)/=size(e).or.size(dw)/=size(e).or.size(df)/=size(e))return
    if(.not.ieee_is_finite(mu))return
    if(.not.salmon_all_finite(de).or..not.salmon_all_finite(dw))return
    call evaluate(e,t,mu,f,b,s)
    susceptibility=sum(w*b)
    ! Do not divide by a numerically unresolved charge susceptibility.
    status=dc_thermal_singular
    if(susceptibility<=128d0*epsilon(t)*sum(w))return
    dmu=(sum(w*b*de)-t*sum(f*dw))/susceptibility
    df=b*((dmu-de)/t)
    dts=g*sum(t*s*dw+(e-mu)*w*df)
    status=0
  end subroutine

  pure logical function finite_real_1d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:)
    real(8) :: value
    integer :: i
    finite=.false.
    do i=1,size(values,1)
      value=values(i)
      if(.not.ieee_is_finite(value))return
    enddo
    finite=.true.
  end function

end module
