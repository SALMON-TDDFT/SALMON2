! Stage bookkeeping only. No RT operator, density or wavefunction mutation.
module exx_surrogate_stage
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer,parameter,public :: surrogate_start=0,surrogate_predictor=1
  integer,parameter,public :: surrogate_corrector_start=2,surrogate_endpoint=3
  type,public :: s_surrogate_stage
    integer :: accepted_step=0,pending_step=-1,stage=surrogate_start,retries=0
    logical :: register_strict=.false.
  end type
  public :: surrogate_begin,surrogate_visit,surrogate_rollback,surrogate_accept
contains
  subroutine surrogate_begin(state,step,ierr)
    type(s_surrogate_stage),intent(inout) :: state
    integer,intent(in) :: step
    integer,intent(out) :: ierr
    ierr=1
    if(state%pending_step/=-1.or.step/=state%accepted_step+1)return
    state%pending_step=step;state%stage=surrogate_start;state%retries=0
    state%register_strict=.false.;ierr=0
  end subroutine

  subroutine surrogate_visit(state,stage,dt,time,ierr)
    type(s_surrogate_stage),intent(inout) :: state
    integer,intent(in) :: stage
    real(real64),intent(in) :: dt
    real(real64),intent(out) :: time
    integer,intent(out) :: ierr
    ierr=1;time=0
    if(state%pending_step==-1.or..not.ieee_is_finite(dt).or.dt<=0)return
    if(stage/=state%stage+1.or.stage>surrogate_endpoint)return
    state%stage=stage
    if(stage==surrogate_corrector_start)then
      time=state%accepted_step*dt
    else
      time=state%pending_step*dt
    endif
    ierr=0
  end subroutine

  subroutine surrogate_rollback(state,ierr)
    type(s_surrogate_stage),intent(inout) :: state
    integer,intent(out) :: ierr
    ierr=1
    if(state%pending_step==-1.or.state%retries/=0)return
    ! Production adapter must restore all physical transaction state separately.
    state%retries=1;state%stage=surrogate_start;state%register_strict=.false.;ierr=0
  end subroutine

  subroutine surrogate_accept(state,strict,ierr)
    type(s_surrogate_stage),intent(inout) :: state
    logical,intent(in) :: strict
    integer,intent(out) :: ierr
    ierr=1;state%register_strict=.false.
    if(state%pending_step==-1.or.state%stage/=surrogate_endpoint)return
    if(state%retries/=0.and..not.strict)return
    state%accepted_step=state%pending_step;state%pending_step=-1
    state%stage=surrogate_start;state%retries=0
    state%register_strict=strict;ierr=0
  end subroutine
end module
