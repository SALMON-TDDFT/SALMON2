#include "config.h"
module dg_hybrid_continuation_controller
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi, only: MPI_Allreduce, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_INTEGER8, MPI_MAX, MPI_MIN, MPI_SUCCESS
#endif
  implicit none
  private
  integer,parameter::residual_channel_count=4

  type,public::s_dg_hybrid_controller_controls
    real(real64)::initial_step=0.125d0,minimum_step=0.015625d0,maximum_step=0.5d0
    real(real64)::growth_factor=1.5d0,shrink_factor=0.5d0,residual_growth_limit=4d0
    real(real64)::density_damping=0.5d0,minimum_projector_overlap=0.9d0
    integer::maximum_rollbacks=8,iteration_limit=100
    real(real64)::intermediate_tolerance(residual_channel_count)=[1d-5,1d-5,1d-5,1d-6]
    real(real64)::final_tolerance(residual_channel_count)=[1d-8,1d-8,1d-8,1d-10]
  end type s_dg_hybrid_controller_controls

  type,public::s_dg_hybrid_trial_state
    real(real64),allocatable::density(:),potential(:),occupations(:),eigenvalues(:),mixing_history(:)
    complex(real64),allocatable::projector(:,:),trace(:,:)
    integer::density_epoch=-1,operator_epoch=-1,projector_epoch=-1,trace_epoch=-1,derived_epoch=-1
    integer(int64)::operator_structure_fingerprint=0_int64,operator_value_fingerprint=0_int64
    logical::trace_cache_valid=.false.
  end type s_dg_hybrid_trial_state

  type,public::s_dg_hybrid_stage_report
    real(real64)::residuals(residual_channel_count)=huge(1d0)
    real(real64)::tolerances(residual_channel_count)=0d0
    real(real64)::projector_overlap=0d0
    integer::iteration=0
    logical::electron_ok=.false.,occupation_ok=.false.,hermitian_ok=.false.,symmetry_ok=.false.
    logical::real_space_ok=.false.,finite_ok=.false.,gap_shrinking=.false.
  end type s_dg_hybrid_stage_report

  type,public::s_dg_hybrid_stage_schedule
    integer::ordinary_limit=0,ordinary_solve_count=0,refresh_solve_count=0
    logical::refresh_pending=.false.,refresh_active=.false.
  end type s_dg_hybrid_stage_schedule

  type,public::s_dg_hybrid_controller
    logical::valid=.false.,trial_active=.false.,trace_valid=.false.,rollback_since_accept=.false.
    real(real64)::accepted_lambda=0d0,trial_lambda=0d0,step=0d0
    integer::rollback_count=0,inner_sample_count=0,growth_streak=0
    type(s_dg_hybrid_controller_controls)::controls
    type(s_dg_hybrid_trial_state)::accepted_state
    real(real64)::previous_inner(residual_channel_count)=0d0
    real(real64),allocatable::face_lambda(:)
  end type s_dg_hybrid_controller

  public::default_dg_hybrid_controller_controls,dg_hybrid_stage_tolerances,&
    validate_dg_hybrid_controller_contract,&
    initialize_dg_hybrid_controller,propose_dg_hybrid_trial,observe_dg_hybrid_inner_residuals,&
    decide_dg_hybrid_stage,reject_dg_hybrid_trial,initialize_dg_hybrid_stage_schedule,&
    begin_dg_hybrid_stage_solve,schedule_dg_hybrid_candidate_checks,complete_dg_hybrid_stage_solve,&
    dg_hybrid_continuation_state_count
contains
  pure subroutine dg_hybrid_continuation_state_count(occupations,basis_count,solve_count,meaningful_gap,ok)
    real(real64),intent(in)::occupations(:)
    integer,intent(in)::basis_count
    integer,intent(out)::solve_count
    logical,intent(out)::meaningful_gap,ok
    real(real64),parameter::occupation_floor=64d0*epsilon(1d0)
    ok=size(occupations)>0.and.basis_count>=size(occupations).and.&
      all(ieee_is_finite(occupations)).and.all(occupations>=0d0).and.all(occupations<=2d0)
    solve_count=0;meaningful_gap=.false.
    if(.not.ok)return
    meaningful_gap=occupations(size(occupations))>occupation_floor.and.basis_count>size(occupations)
    solve_count=size(occupations)+merge(1,0,meaningful_gap)
  end subroutine dg_hybrid_continuation_state_count

  pure subroutine initialize_dg_hybrid_stage_schedule(iteration_limit,schedule)
    integer,intent(in)::iteration_limit
    type(s_dg_hybrid_stage_schedule),intent(out)::schedule
    schedule=s_dg_hybrid_stage_schedule()
    schedule%ordinary_limit=max(0,iteration_limit)
  end subroutine initialize_dg_hybrid_stage_schedule

  pure subroutine begin_dg_hybrid_stage_solve(schedule,run_solve,iteration)
    type(s_dg_hybrid_stage_schedule),intent(inout)::schedule
    logical,intent(out)::run_solve
    integer,intent(out)::iteration
    run_solve=.false.;iteration=schedule%ordinary_solve_count+schedule%refresh_solve_count
    if(schedule%refresh_pending)then
      schedule%refresh_pending=.false.;schedule%refresh_active=.true.
      schedule%refresh_solve_count=schedule%refresh_solve_count+1
      iteration=schedule%ordinary_limit+schedule%refresh_solve_count;run_solve=.true.
    else if(.not.schedule%refresh_active.and.schedule%ordinary_solve_count<schedule%ordinary_limit)then
      schedule%ordinary_solve_count=schedule%ordinary_solve_count+1
      iteration=schedule%ordinary_solve_count;run_solve=.true.
    endif
  end subroutine begin_dg_hybrid_stage_solve

  pure subroutine schedule_dg_hybrid_candidate_checks(schedule,cheap_candidate,run_expensive)
    type(s_dg_hybrid_stage_schedule),intent(in)::schedule
    logical,intent(in)::cheap_candidate
    logical,intent(out)::run_expensive
    run_expensive=cheap_candidate
  end subroutine schedule_dg_hybrid_candidate_checks

  pure subroutine complete_dg_hybrid_stage_solve(schedule,stage_converged,trial_lambda,&
      refresh_scheduled,final_refresh_performed)
    type(s_dg_hybrid_stage_schedule),intent(inout)::schedule
    logical,intent(in)::stage_converged
    real(real64),intent(in)::trial_lambda
    logical,intent(out)::refresh_scheduled,final_refresh_performed
    refresh_scheduled=.false.;final_refresh_performed=.false.
    if(schedule%refresh_active)then
      final_refresh_performed=stage_converged;schedule%refresh_active=.false.
    else if(stage_converged.and.trial_lambda==1d0.and.schedule%refresh_solve_count==0)then
      schedule%refresh_pending=.true.;refresh_scheduled=.true.
    endif
  end subroutine complete_dg_hybrid_stage_solve

  subroutine validate_dg_hybrid_controller_contract(icomm,controls,ok,message)
    integer,intent(in)::icomm;type(s_dg_hybrid_controller_controls),intent(in)::controls
    logical,intent(out)::ok;character(*),intent(out)::message
#ifdef USE_MPI
    integer::local_bad,global_bad,ierr
    integer(int64)::local_hash,minimum_hash,maximum_hash
    local_bad=merge(0,1,valid_controls(controls))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ok=.false.;message='invalid continuation controller controls';return;endif
    local_hash=controller_controls_fingerprint(controls)
    call MPI_Allreduce(local_hash,minimum_hash,1,MPI_INTEGER8,MPI_MIN,icomm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_hash,maximum_hash,1,MPI_INTEGER8,MPI_MAX,icomm,ierr)
    ok=ierr==MPI_SUCCESS.and.minimum_hash==maximum_hash
    if(ok)then;message='';else;message='rank-disagreeing continuation controller controls';endif
#else
    ok=.false.;message='continuation controller validation requires MPI'
#endif
  end subroutine validate_dg_hybrid_controller_contract

  subroutine default_dg_hybrid_controller_controls(controls)
    type(s_dg_hybrid_controller_controls),intent(out)::controls
    controls=s_dg_hybrid_controller_controls()
  end subroutine default_dg_hybrid_controller_controls

  subroutine dg_hybrid_stage_tolerances(controls,lambda,tolerances)
    type(s_dg_hybrid_controller_controls),intent(in)::controls
    real(real64),intent(in)::lambda
    real(real64),intent(out)::tolerances(residual_channel_count)
    real(real64)::bounded_lambda
    bounded_lambda=min(1d0,max(0d0,lambda))
    tolerances=max(controls%final_tolerance,(1d0-bounded_lambda)*controls%intermediate_tolerance+&
      bounded_lambda*controls%final_tolerance)
  end subroutine dg_hybrid_stage_tolerances

  subroutine initialize_dg_hybrid_controller(icomm,controls,accepted_lambda,accepted_state,face_count,controller,ok,message)
    integer,intent(in)::icomm,face_count
    type(s_dg_hybrid_controller_controls),intent(in)::controls
    real(real64),intent(in)::accepted_lambda
    type(s_dg_hybrid_trial_state),intent(in)::accepted_state
    type(s_dg_hybrid_controller),intent(out)::controller
    logical,intent(out)::ok;character(*),intent(out)::message
#ifdef USE_MPI
    integer::local_bad,global_bad,ierr
    integer(int64)::local_hash,minimum_hash,maximum_hash
    real(real64)::minimum_lambda,maximum_lambda
    local_bad=merge(0,1,valid_controls(controls).and.accepted_lambda>=0d0.and.accepted_lambda<=1d0.and.&
      ieee_is_finite(accepted_lambda).and.face_count>=0.and.valid_state(accepted_state).and.&
      accepted_state%trace_cache_valid)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ok=.false.;message='invalid continuation controller initialization';return;endif
    call MPI_Allreduce(accepted_lambda,minimum_lambda,1,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(accepted_lambda,maximum_lambda,1,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_lambda/=maximum_lambda)then
      ok=.false.;message='rank-disagreeing accepted lambda';return
    endif
    local_hash=controller_controls_fingerprint(controls)
    call MPI_Allreduce(local_hash,minimum_hash,1,MPI_INTEGER8,MPI_MIN,icomm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_hash,maximum_hash,1,MPI_INTEGER8,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      ok=.false.;message='rank-disagreeing continuation controller controls';return
    endif
    controller%controls=controls;controller%accepted_lambda=accepted_lambda;controller%trial_lambda=accepted_lambda
    controller%step=controls%initial_step;controller%accepted_state=accepted_state
    allocate(controller%face_lambda(face_count));controller%face_lambda=accepted_lambda
    controller%valid=.true.;controller%trial_active=.false.
    controller%trace_valid=controller%accepted_state%trace_cache_valid;ok=.true.;message=''
#else
    ok=.false.;message='continuation controller requires MPI'
#endif
  end subroutine initialize_dg_hybrid_controller

  subroutine propose_dg_hybrid_trial(icomm,controller,state,ok,message)
    integer,intent(in)::icomm
    type(s_dg_hybrid_controller),intent(inout)::controller
    type(s_dg_hybrid_trial_state),intent(inout)::state
    logical,intent(out)::ok;character(*),intent(out)::message
#ifdef USE_MPI
    integer::local_bad,global_bad,ierr
    local_bad=merge(0,1,controller%valid.and..not.controller%trial_active.and.controller%accepted_lambda<1d0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ok=.false.;message='cannot propose continuation trial';return;endif
    state=controller%accepted_state
    state%trace_cache_valid=.false.
    controller%trial_lambda=min(1d0,controller%accepted_lambda+controller%step)
    controller%face_lambda=controller%trial_lambda;controller%trial_active=.true.;controller%trace_valid=.false.
    controller%inner_sample_count=0;controller%growth_streak=0;controller%previous_inner=0d0
    ok=.true.;message=''
#else
    ok=.false.;message='continuation trial proposal requires MPI'
#endif
  end subroutine propose_dg_hybrid_trial

  subroutine observe_dg_hybrid_inner_residuals(icomm,controller,residuals,reject_requested,ok,message)
    integer,intent(in)::icomm
    type(s_dg_hybrid_controller),intent(inout)::controller
    real(real64),intent(in)::residuals(residual_channel_count)
    logical,intent(out)::reject_requested,ok;character(*),intent(out)::message
#ifdef USE_MPI
    integer::local_bad,global_bad,ierr
    real(real64)::growth,numerical_floor,global_residuals(residual_channel_count)
    local_bad=merge(0,1,controller%valid.and.controller%trial_active.and.all(ieee_is_finite(residuals)).and.&
      all(residuals>=0d0))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      reject_requested=.false.;ok=.false.;message='invalid inner continuation residual sample';return
    endif
    call MPI_Allreduce(residuals,global_residuals,residual_channel_count,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS)then;reject_requested=.false.;ok=.false.;message='inner residual canonicalization failed';return;endif
    numerical_floor=sqrt(tiny(1d0));reject_requested=.false.
    if(controller%inner_sample_count>0)then
      growth=maxval(global_residuals/max(controller%previous_inner,numerical_floor))
      if(growth>controller%controls%residual_growth_limit)then
        controller%growth_streak=controller%growth_streak+1
      else;controller%growth_streak=0
      endif
      reject_requested=controller%growth_streak>=2
    endif
    controller%previous_inner=global_residuals;controller%inner_sample_count=controller%inner_sample_count+1
    ok=.true.;message=''
#else
    reject_requested=.false.;ok=.false.;message='inner residual observation requires MPI'
#endif
  end subroutine observe_dg_hybrid_inner_residuals

  subroutine decide_dg_hybrid_stage(icomm,controller,state,report,accept,ok,message)
    integer,intent(in)::icomm
    type(s_dg_hybrid_controller),intent(inout)::controller
    type(s_dg_hybrid_trial_state),intent(in)::state
    type(s_dg_hybrid_stage_report),intent(in)::report
    logical,intent(out)::accept,ok;character(*),intent(out)::message
#ifdef USE_MPI
    logical::local_accept
    integer::local_integer,global_integer,ierr,global_iteration,gap_integer,global_gap_integer
    real(real64)::stage_tolerances(residual_channel_count)
    call dg_hybrid_stage_tolerances(controller%controls,controller%trial_lambda,stage_tolerances)
    local_accept=controller%valid.and.controller%trial_active.and.valid_state(state).and.state%trace_cache_valid.and.&
      report%iteration>=1.and.&
      report%iteration<=controller%controls%iteration_limit.and.all(ieee_is_finite(report%residuals)).and.&
      all(report%residuals>=0d0).and.all(report%residuals<=stage_tolerances).and.ieee_is_finite(report%projector_overlap).and.&
      report%projector_overlap>=controller%controls%minimum_projector_overlap.and.report%electron_ok.and.&
      report%occupation_ok.and.report%hermitian_ok.and.report%symmetry_ok.and.report%real_space_ok.and.report%finite_ok
    local_integer=merge(1,0,local_accept)
    call MPI_Allreduce(local_integer,global_integer,1,MPI_INTEGER,MPI_MIN,icomm,ierr)
    if(ierr/=MPI_SUCCESS)then;accept=.false.;ok=.false.;message='stage decision reduction failed';return;endif
    call MPI_Allreduce(report%iteration,global_iteration,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    gap_integer=merge(1,0,report%gap_shrinking)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(gap_integer,global_gap_integer,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS)then;accept=.false.;ok=.false.;message='adaptive-stage canonicalization failed';return;endif
    accept=global_integer==1;ok=.true.;message=''
    if(.not.accept)return
    controller%accepted_state=state
    controller%accepted_lambda=controller%trial_lambda
    if(global_gap_integer/=0)then
      controller%step=max(controller%controls%minimum_step,controller%step*controller%controls%shrink_factor)
    else if(global_iteration<=controller%controls%iteration_limit/2.and..not.controller%rollback_since_accept)then
      controller%step=min(controller%controls%maximum_step,controller%step*controller%controls%growth_factor)
    endif
    controller%trial_active=.false.;controller%trace_valid=controller%accepted_state%trace_cache_valid
    controller%growth_streak=0
    controller%rollback_since_accept=.false.
#else
    accept=.false.;ok=.false.;message='continuation stage decision requires MPI'
#endif
  end subroutine decide_dg_hybrid_stage

  subroutine reject_dg_hybrid_trial(icomm,controller,state,reason,ok,message)
    integer,intent(in)::icomm
    type(s_dg_hybrid_controller),intent(inout)::controller
    type(s_dg_hybrid_trial_state),intent(inout)::state
    character(*),intent(in)::reason
    logical,intent(out)::ok;character(*),intent(out)::message
#ifdef USE_MPI
    integer::local_bad,global_bad,ierr
    local_bad=merge(0,1,controller%valid.and.controller%trial_active.and.len_trim(reason)>0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ok=.false.;message='invalid continuation rollback request';return;endif
    if(controller%rollback_count>=controller%controls%maximum_rollbacks)then
      ok=.false.;message='continuation rollback limit exhausted';return
    endif
    state=controller%accepted_state;controller%trial_lambda=controller%accepted_lambda
    controller%step=max(controller%controls%minimum_step,controller%step*controller%controls%shrink_factor)
    controller%rollback_count=controller%rollback_count+1;controller%trial_active=.false.
    controller%rollback_since_accept=.true.
    controller%trace_valid=state%trace_cache_valid;controller%inner_sample_count=0;controller%growth_streak=0
    controller%face_lambda=controller%accepted_lambda;ok=.true.;message=''
#else
    ok=.false.;message='continuation rollback requires MPI'
#endif
  end subroutine reject_dg_hybrid_trial

  logical function valid_controls(controls) result(valid)
    type(s_dg_hybrid_controller_controls),intent(in)::controls
    valid=all(ieee_is_finite([controls%initial_step,controls%minimum_step,controls%maximum_step,&
      controls%growth_factor,controls%shrink_factor,controls%residual_growth_limit,controls%density_damping,&
      controls%minimum_projector_overlap])).and.all(ieee_is_finite(controls%intermediate_tolerance)).and.&
      all(ieee_is_finite(controls%final_tolerance)).and.controls%minimum_step>0d0.and.&
      controls%initial_step>=controls%minimum_step.and.&
      controls%initial_step<=controls%maximum_step.and.controls%maximum_step<=1d0.and.controls%growth_factor>1d0.and.&
      controls%shrink_factor>0d0.and.controls%shrink_factor<1d0.and.controls%residual_growth_limit>1d0.and.&
      controls%density_damping>0d0.and.controls%density_damping<=1d0.and.controls%minimum_projector_overlap>=0d0.and.&
      controls%minimum_projector_overlap<=1d0.and.controls%maximum_rollbacks>=0.and.controls%iteration_limit>1.and.&
      all(controls%intermediate_tolerance>=controls%final_tolerance).and.all(controls%final_tolerance>0d0)
  end function valid_controls

  integer(int64) function controller_controls_fingerprint(controls) result(hash)
    type(s_dg_hybrid_controller_controls),intent(in)::controls
    integer::i
    hash=int(z'BB67AE8584CAA73B',int64)
    call mix(transfer(controls%initial_step,hash));call mix(transfer(controls%minimum_step,hash))
    call mix(transfer(controls%maximum_step,hash));call mix(transfer(controls%growth_factor,hash))
    call mix(transfer(controls%shrink_factor,hash));call mix(transfer(controls%residual_growth_limit,hash))
    call mix(transfer(controls%density_damping,hash));call mix(transfer(controls%minimum_projector_overlap,hash))
    call mix(int(controls%maximum_rollbacks,int64));call mix(int(controls%iteration_limit,int64))
    do i=1,residual_channel_count
      call mix(transfer(controls%intermediate_tolerance(i),hash));call mix(transfer(controls%final_tolerance(i),hash))
    enddo
  contains
    subroutine mix(value)
      integer(int64),intent(in)::value
      hash=ieor(ishftc(hash,9),value)
    end subroutine mix
  end function controller_controls_fingerprint

  logical function valid_state(state) result(valid)
    type(s_dg_hybrid_trial_state),intent(in)::state
    valid=allocated(state%density).and.allocated(state%potential).and.allocated(state%projector).and.&
      allocated(state%trace).and.allocated(state%occupations).and.allocated(state%eigenvalues).and.&
      allocated(state%mixing_history)
    if(.not.valid)return
    valid=size(state%projector)>0.and.size(state%occupations)>0.and.&
      size(state%eigenvalues)==size(state%occupations).and.&
      all(ieee_is_finite(state%density)).and.all(ieee_is_finite(state%potential)).and.&
      all(ieee_is_finite(real(state%projector))).and.all(ieee_is_finite(aimag(state%projector))).and.&
      all(ieee_is_finite(real(state%trace))).and.all(ieee_is_finite(aimag(state%trace))).and.&
      all(ieee_is_finite(state%occupations)).and.all(ieee_is_finite(state%eigenvalues)).and.&
      all(ieee_is_finite(state%mixing_history))
  end function valid_state

end module dg_hybrid_continuation_controller
