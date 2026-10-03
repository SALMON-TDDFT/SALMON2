! Serial fixed-Q inference controller. No automatic activation or calibration.
module exx_surrogate_controller
  use iso_fortran_env,only:real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use exx_surrogate_model,only:surrogate_predict
  use exx_surrogate_guard,only:surrogate_psd,surrogate_action_check
  implicit none
  private
  type,public::s_surrogate_controller
    integer::epoch=-1,p=0,period=0,count=0,last_exact=-1
    logical::bundle_ready=.false.
    real(real64)::dt=0,dv=0
    complex(real64),allocatable::q(:,:),history(:,:,:)
    real(real64),allocatable::coefficients(:,:)
    integer,allocatable::steps(:)
    logical,allocatable::certified(:)
  end type
  public::surrogate_initialize,surrogate_record,surrogate_set_class,surrogate_evaluate,surrogate_certify_bundle
contains
  logical function valid_matrix(b)result(ok)
    complex(real64),intent(in)::b(:,:)
    integer::i,j
    ok=.false.
    do j=1,size(b,2)
      do i=1,size(b,1)
        if(.not.ieee_is_finite(real(b(i,j),real64)).or..not.ieee_is_finite(aimag(b(i,j))))return
      enddo
    enddo
    ok=.true.
  end function

  subroutine surrogate_initialize(state,q,epoch,p,period,dt,dv,status)
    type(s_surrogate_controller),intent(inout)::state
    complex(real64),intent(in)::q(:,:)
    integer,intent(in)::epoch,p,period
    real(real64),intent(in)::dt,dv
    integer,intent(out)::status
    type(s_surrogate_controller)::fresh
    complex(real64),allocatable::gram(:,:)
    integer::k,i
    status=1;k=size(q,2)
    if(epoch<0.or.p<2.or.p>4.or.period<1.or.k<1)return
    if(.not.all(ieee_is_finite([dt,dv])).or.dt<=0.or.dv<=0.or..not.valid_matrix(q))return
    gram=matmul(conjg(transpose(q)),q)*dv
    do i=1,k
      gram(i,i)=gram(i,i)-1
    enddo
    if(maxval(abs(gram))>1d-10)return
    fresh%q=q;fresh%epoch=epoch;fresh%p=p;fresh%period=period;fresh%dt=dt;fresh%dv=dv
    allocate(fresh%history(k,k,p),fresh%steps(p),fresh%coefficients(p-1,period),fresh%certified(period))
    fresh%history=0;fresh%steps=-1;fresh%coefficients=0;fresh%certified=.false.
    state=fresh;status=0
  end subroutine

  subroutine surrogate_record(state,step,b,strict,accepted,status)
    type(s_surrogate_controller),intent(inout)::state
    integer,intent(in)::step
    complex(real64),intent(in)::b(:,:)
    logical,intent(in)::strict,accepted
    integer,intent(out)::status
    integer::k,p
    status=1
    if(.not.allocated(state%q).or..not.strict.or..not.accepted)return
    if(step<0.or.step<=state%last_exact)return
    k=size(state%q,2);p=state%p
    if(any(shape(b)/=[k,k]).or..not.valid_matrix(b))return
    if(maxval(abs(b-conjg(transpose(b))))>1d-12*max(maxval(abs(b)),1d0))return
    ! Nonperiodic exact endpoints are excluded, with recertification required.
    state%last_exact=step
    if(mod(step,state%period)/=0)then
      state%certified=.false.;state%bundle_ready=.false.;status=0;return
    endif
    if(state%count>0)then
      if(step-state%steps(state%count)/=state%period)then
        state%count=0;state%certified=.false.;state%bundle_ready=.false.
      endif
    endif
    if(state%count==p)then
      state%history(:,:,:p-1)=state%history(:,:,2:p)
      state%steps(:p-1)=state%steps(2:p);state%count=p-1
    endif
    state%count=state%count+1;state%history(:,:,state%count)=b;state%steps(state%count)=step
    status=0
  end subroutine

  subroutine surrogate_set_class(state,horizon,coeff,validated,status)
    type(s_surrogate_controller),intent(inout)::state
    integer,intent(in)::horizon
    real(real64),intent(in)::coeff(:)
    logical,intent(in)::validated
    integer,intent(out)::status
    status=1
    if(.not.allocated(state%q).or..not.validated)return
    if(horizon<1.or.horizon>state%period.or.size(coeff)/=state%p-1)return
    if(state%certified(horizon))return
    if(.not.all(ieee_is_finite(coeff)))return
    ! Caller is responsible for offline/shadow certification, never RT self-training.
    state%coefficients(:,horizon)=coeff;state%certified(horizon)=.true.;status=0
  end subroutine

  subroutine surrogate_certify_bundle(state,coeff,validated,status)
    type(s_surrogate_controller),intent(inout)::state
    real(real64),intent(in)::coeff(:,:)
    logical,intent(in)::validated
    integer,intent(out)::status
    status=1
    if(.not.allocated(state%q).or..not.validated.or.state%bundle_ready)return
    if(any(shape(coeff)/=[state%p-1,state%period]))return
    if(.not.all(ieee_is_finite(coeff)))return
    if(state%count/=state%p)return
    if(mod(state%steps(state%p),state%period)/=0)return
    if(state%last_exact/=state%steps(state%p))return
    ! Whole bundle installed only after independent caller-side shadow/test checks.
    state%coefficients=coeff;state%certified=.true.;state%bundle_ready=.true.;status=0
  end subroutine

  subroutine surrogate_evaluate(state,step,psi,exact,correction_max,atol,rtol,floor, &
    predicted,correction,absolute,relative,accepted,status)
    type(s_surrogate_controller),intent(in)::state
    integer,intent(in)::step
    complex(real64),intent(in)::psi(:,:),exact(:,:)
    real(real64),intent(in)::correction_max,atol,rtol,floor
    complex(real64),intent(out)::predicted(:,:)
    real(real64),intent(out)::correction,absolute,relative
    logical,intent(out)::accepted
    integer,intent(out)::status
    integer::h
    status=1;accepted=.false.;predicted=0
    correction=huge(1d0);absolute=huge(1d0);relative=huge(1d0)
    if(.not.allocated(state%q).or.state%count/=state%p)return
    h=step-state%steps(state%p)
    if(h<1.or.h>state%period)return
    if(.not.state%bundle_ready.or..not.state%certified(h))return
    call surrogate_predict(state%history,real(state%steps,real64)*state%dt,step*state%dt, &
      state%coefficients(:,h),predicted,status)
    if(status/=0)return
    call surrogate_psd(predicted,correction_max,floor,correction,status)
    if(status/=0)return
    call surrogate_action_check(state%q,predicted,psi,exact,state%dv,atol,rtol,floor, &
      absolute,relative,accepted,status)
  end subroutine
end module
