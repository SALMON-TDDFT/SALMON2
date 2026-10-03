! Strict-shadow action coverage only; not a physical/spectral certification.
module exx_surrogate_shadow
  use iso_fortran_env,only:real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  integer,parameter,public::shadow_beginning=1,shadow_predictor=2,shadow_endpoint=3
  type,public::s_surrogate_shadow
    integer::period=0,min_count=0,last_step=-1,last_stage=0,last_horizon=-1
    real(real64)::atol=0,rtol=0
    integer,allocatable::counts(:,:),failures(:,:)
    real(real64),allocatable::max_absolute(:,:),max_relative(:,:)
  end type
  public::shadow_initialize,shadow_observe,shadow_complete
contains
  subroutine shadow_initialize(state,period,min_count,atol,rtol,status)
    type(s_surrogate_shadow),intent(inout)::state
    integer,intent(in)::period,min_count
    real(real64),intent(in)::atol,rtol
    integer,intent(out)::status
    type(s_surrogate_shadow)::fresh
    status=1
    if(period<1.or.min_count<1)return
    if(.not.all(ieee_is_finite([atol,rtol])).or.atol<0.or.rtol<0)return
    fresh%period=period;fresh%min_count=min_count;fresh%atol=atol;fresh%rtol=rtol
    allocate(fresh%counts(0:period,3),fresh%failures(0:period,3), &
      fresh%max_absolute(0:period,3),fresh%max_relative(0:period,3))
    fresh%counts=0;fresh%failures=0;fresh%max_absolute=0;fresh%max_relative=0
    state=fresh;status=0
  end subroutine

  subroutine shadow_observe(state,step,stage,horizon,absolute,relative,strict_trajectory,status)
    type(s_surrogate_shadow),intent(inout)::state
    integer,intent(in)::step,stage,horizon
    real(real64),intent(in)::absolute,relative
    logical,intent(in)::strict_trajectory
    integer,intent(out)::status
    integer::h
    status=1
    if(.not.allocated(state%counts).or..not.strict_trajectory)return
    if(step<1.or.stage<1.or.stage>3)return
    if(horizon<0.or.horizon>state%period)return
    if(stage==1.and.horizon==state%period)return
    if(stage>1.and.horizon==0)return
    if(.not.all(ieee_is_finite([absolute,relative])).or.absolute<0.or.relative<0)return
    if(stage==shadow_beginning)then
      if(state%last_stage/=0.and.state%last_stage/=shadow_endpoint)return
      if(step<=state%last_step)return
    else
      if(state%last_stage/=stage-1.or.step/=state%last_step)return
      if(stage==shadow_predictor.and.horizon/=state%last_horizon+1)return
      if(stage==shadow_endpoint.and.horizon/=state%last_horizon)return
    endif
    h=horizon
    state%counts(h,stage)=state%counts(h,stage)+1
    state%max_absolute(h,stage)=max(state%max_absolute(h,stage),absolute)
    state%max_relative(h,stage)=max(state%max_relative(h,stage),relative)
    if(absolute>state%atol.or.relative>state%rtol)state%failures(h,stage)=state%failures(h,stage)+1
    state%last_step=step;state%last_stage=stage;state%last_horizon=horizon;status=0
  end subroutine

  logical function shadow_complete(state)result(complete)
    type(s_surrogate_shadow),intent(in)::state
    complete=.false.
    if(.not.allocated(state%counts))return
    if(state%last_stage/=shadow_endpoint)return
    complete=all(state%counts(0:state%period-1,1)>=state%min_count).and. &
      all(state%counts(1:state%period,shadow_predictor)>=state%min_count).and. &
      all(state%counts(1:state%period,shadow_endpoint)>=state%min_count).and.all(state%failures==0)
  end function
end module
