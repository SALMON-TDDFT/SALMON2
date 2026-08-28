#include "config.h"
module dg_hybrid_continuation_scf
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_continuation_state,only:s_dg_hybrid_continuation_state
  use dg_hybrid_continuation_controller,only:s_dg_hybrid_controller_controls,s_dg_hybrid_trial_state,&
    s_dg_hybrid_stage_report,s_dg_hybrid_controller,initialize_dg_hybrid_controller,&
    propose_dg_hybrid_trial,observe_dg_hybrid_inner_residuals,decide_dg_hybrid_stage,reject_dg_hybrid_trial,&
    dg_hybrid_stage_tolerances,validate_dg_hybrid_controller_contract
  use dg_hybrid_continuation_residuals,only:s_dg_hybrid_residuals
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  abstract interface
    subroutine volume_callback(lambda,density,state,ok)
      import real64,s_dg_hybrid_trial_state
      real(real64),intent(in)::lambda,density(:);type(s_dg_hybrid_trial_state),intent(inout)::state
      logical,intent(out)::ok
    end subroutine volume_callback
    subroutine solve_callback(lambda,iteration,state,ok)
      import real64,s_dg_hybrid_trial_state
      real(real64),intent(in)::lambda;integer,intent(in)::iteration
      type(s_dg_hybrid_trial_state),intent(inout)::state;logical,intent(out)::ok
    end subroutine solve_callback
    subroutine projector_callback(lambda,state,overlap,ok)
      import real64,s_dg_hybrid_trial_state
      real(real64),intent(in)::lambda;type(s_dg_hybrid_trial_state),intent(inout)::state
      real(real64),intent(out)::overlap;logical,intent(out)::ok
    end subroutine projector_callback
    subroutine density_trace_callback(lambda,input_density,state,output_density,ok)
      import real64,s_dg_hybrid_trial_state
      real(real64),intent(in)::lambda,input_density(:);type(s_dg_hybrid_trial_state),intent(inout)::state
      real(real64),intent(out)::output_density(:);logical,intent(out)::ok
    end subroutine density_trace_callback
    subroutine residual_callback(lambda,iteration,input_density,input_trace,state,residuals,projector_overlap,&
        electron_ok,occupation_ok,hermitian_ok,symmetry_ok,real_space_ok,gap_shrinking,ok)
      import real64,s_dg_hybrid_trial_state,s_dg_hybrid_residuals
      real(real64),intent(in)::lambda,input_density(:);integer,intent(in)::iteration
      complex(real64),intent(in)::input_trace(:,:);type(s_dg_hybrid_trial_state),intent(in)::state
      type(s_dg_hybrid_residuals),intent(out)::residuals;real(real64),intent(out)::projector_overlap
      logical,intent(out)::electron_ok,occupation_ok,hermitian_ok,symmetry_ok,real_space_ok,gap_shrinking,ok
    end subroutine residual_callback
    subroutine mix_callback(iteration,input_density,output_density,damping,mixed_density,ok)
      import real64
      integer,intent(in)::iteration;real(real64),intent(in)::input_density(:),output_density(:),damping
      real(real64),intent(out)::mixed_density(:);logical,intent(out)::ok
    end subroutine mix_callback
  end interface
  type,public::s_dg_hybrid_continuation_callbacks
    procedure(volume_callback),pointer,nopass::build_volume=>null()
    procedure(solve_callback),pointer,nopass::solve_full=>null()
    procedure(projector_callback),pointer,nopass::refresh_projector=>null()
    procedure(density_trace_callback),pointer,nopass::refresh_density_trace=>null()
    procedure(residual_callback),pointer,nopass::evaluate_residuals=>null()
    procedure(mix_callback),pointer,nopass::mix_density=>null()
    type(s_dg_hybrid_trial_state)::seed_state
    integer::face_count=0,maximum_inner_iterations=0,accepted_stages=0,rollback_count=0
    real(real64)::final_lambda=0d0
  end type s_dg_hybrid_continuation_callbacks
  public::run_dg_hybrid_continuation_scf
contains
  subroutine run_dg_hybrid_continuation_scf(comm,state,controls,callbacks,accepted_state,ok,message)
    integer,intent(in)::comm;type(s_dg_hybrid_continuation_state),intent(in)::state
    type(s_dg_hybrid_controller_controls),intent(in)::controls
    type(s_dg_hybrid_continuation_callbacks),intent(inout)::callbacks
    type(s_dg_hybrid_trial_state),intent(out)::accepted_state
    logical,intent(out)::ok;character(*),intent(out)::message
    logical::callbacks_complete
#ifdef USE_MPI
    integer::local_bad,global_bad,ierr
#endif
    callbacks_complete=associated(callbacks%build_volume).and.associated(callbacks%solve_full).and.&
      associated(callbacks%refresh_projector).and.associated(callbacks%refresh_density_trace).and.&
      associated(callbacks%evaluate_residuals).and.associated(callbacks%mix_density)
#ifdef USE_MPI
    local_bad=merge(0,1,callbacks_complete)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      ok=.false.;message='collectively incomplete mandatory continuation callback bundle';return
    endif
#else
    if(.not.callbacks_complete)then
      ok=.false.;message='incomplete mandatory continuation callback bundle';return
    endif
#endif
    call run_dg_hybrid_coupled_fixed_points(comm,state,controls,callbacks%seed_state,callbacks%face_count,&
      callbacks%build_volume,callbacks%solve_full,callbacks%refresh_projector,callbacks%refresh_density_trace,&
      callbacks%evaluate_residuals,callbacks%mix_density,callbacks%maximum_inner_iterations,accepted_state,&
      callbacks%final_lambda,callbacks%accepted_stages,callbacks%rollback_count,ok,message)
  end subroutine run_dg_hybrid_continuation_scf

  subroutine run_dg_hybrid_coupled_fixed_points(comm,continuation,controls,seed_state,face_count,&
      build_volume,solve_full,refresh_projector,refresh_density_trace,evaluate_residuals,mix_density,&
      maximum_inner_iterations,final_state,final_lambda,accepted_stages,rollback_count,ok,message)
    integer,intent(in)::comm,face_count,maximum_inner_iterations
    type(s_dg_hybrid_continuation_state),intent(in)::continuation
    type(s_dg_hybrid_controller_controls),intent(in)::controls
    type(s_dg_hybrid_trial_state),intent(in)::seed_state
    procedure(volume_callback)::build_volume
    procedure(solve_callback)::solve_full
    procedure(projector_callback)::refresh_projector
    procedure(density_trace_callback)::refresh_density_trace
    procedure(residual_callback)::evaluate_residuals
    procedure(mix_callback)::mix_density
    type(s_dg_hybrid_trial_state),intent(out)::final_state
    real(real64),intent(out)::final_lambda
    integer,intent(out)::accepted_stages,rollback_count
    logical,intent(out)::ok;character(*),intent(out)::message
#ifdef USE_MPI
    type(s_dg_hybrid_controller)::controller
    type(s_dg_hybrid_trial_state)::state
    type(s_dg_hybrid_trial_state)::input_gamma_state
    type(s_dg_hybrid_residuals)::residuals
    type(s_dg_hybrid_stage_report)::report
    real(real64),allocatable::input_density(:),output_density(:),mixed_density(:)
    complex(real64),allocatable::input_trace(:,:)
    real(real64)::projector_overlap,tolerances(4),lambda
    integer::iteration,ierr,local_bad,global_bad,minimum_iterations,maximum_iterations,previous_operator_epoch,&
      minimum_operator_epoch,maximum_operator_epoch
    integer(int64)::frozen_operator_structure,minimum_operator_structure,maximum_operator_structure,&
      minimum_operator_value,maximum_operator_value
    logical::callback_ok,stage_ok,reject_requested,decision_ok,rejected,fatal_error
    character(256)::controller_message
    ok=.false.;message='';final_lambda=0d0;accepted_stages=0;rollback_count=0
    frozen_operator_structure=0_int64
    call validate_dg_hybrid_controller_contract(comm,controls,decision_ok,controller_message)
    if(.not.decision_ok)then;message=trim(controller_message);return;endif
    call MPI_Allreduce(maximum_inner_iterations,minimum_iterations,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(maximum_inner_iterations,maximum_iterations,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_iterations/=maximum_iterations)then
      message='rank-disagreeing continuation inner iteration limit';return
    endif
    local_bad=merge(0,1,continuation%valid.and.continuation%lambda==0d0.and.&
      allocated(continuation%seed_density).and.allocated(continuation%seed_density_ids).and.&
      continuation%global_density_count>0.and.size(continuation%seed_density)==size(seed_state%density).and.&
      allocated(seed_state%trace).and.seed_state%trace_cache_valid.and.maximum_inner_iterations>0.and.face_count>=0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid coupled continuation SCF contract';return;endif
    allocate(input_density(size(continuation%seed_density)),output_density(size(continuation%seed_density)),&
      mixed_density(size(continuation%seed_density)),&
      input_trace(size(seed_state%trace,1),size(seed_state%trace,2)))
    state=seed_state;input_density=continuation%seed_density;input_trace=seed_state%trace

    ! Lambda zero is not assumed converged merely because the immutable DC seed exists.
    lambda=0d0
    call converge_current(.false.,stage_ok,rejected,fatal_error)
    if(fatal_error)return
    if(.not.stage_ok)then;message='lambda-zero coupled fixed point did not converge';return;endif
    accepted_stages=1
    call initialize_dg_hybrid_controller(comm,controls,0d0,state,face_count,controller,decision_ok,controller_message)
    if(.not.decision_ok)then;message=trim(controller_message);return;endif

    do while(controller%accepted_lambda<1d0)
      call propose_dg_hybrid_trial(comm,controller,state,decision_ok,controller_message)
      if(.not.decision_ok)then;message=trim(controller_message);return;endif
      lambda=controller%trial_lambda;input_density=controller%accepted_state%density
      input_trace=controller%accepted_state%trace;rejected=.false.
      call converge_current(.true.,stage_ok,rejected,fatal_error)
      if(fatal_error)return
      if(rejected)cycle
      if(.not.stage_ok)then
        call reject_dg_hybrid_trial(comm,controller,state,'inner iteration limit',decision_ok,controller_message)
        if(.not.decision_ok)then;message=trim(controller_message);return;endif
        cycle
      endif
      accepted_stages=accepted_stages+1
    enddo
    rollback_count=controller%rollback_count

    ! Rebuild the complete lambda-one operator and observables once more without mixing.
    state=controller%accepted_state;input_gamma_state=state
    lambda=1d0;input_density=state%density;input_trace=state%trace
    call execute_callbacks(1,callback_ok)
    if(.not.callback_ok)then;message='lambda-one final refresh callback failed';return;endif
    call fill_report(1)
    call stage_consensus(report,report%tolerances,stage_ok)
    if(.not.stage_ok)then;message='lambda-one fully refreshed residual gate failed';return;endif
    state%density=input_gamma_state%density;state%trace=input_gamma_state%trace
    state%density_epoch=input_gamma_state%density_epoch;state%trace_epoch=input_gamma_state%trace_epoch
    state%derived_epoch=input_gamma_state%derived_epoch;state%trace_cache_valid=.true.
    final_state=state;final_lambda=1d0
    ok=.true.;message=''
  contains
    subroutine converge_current(use_controller,converged,rejected_trial,fatal)
      logical,intent(in)::use_controller;logical,intent(out)::converged,rejected_trial,fatal
      converged=.false.;rejected_trial=.false.;fatal=.false.
      do iteration=1,maximum_inner_iterations
        call execute_callbacks(iteration,callback_ok)
        if(.not.callback_ok)then
          if(len_trim(message)==0)message='coupled continuation callback failed'
          fatal=.true.;return
        endif
        call fill_report(iteration)
        if(use_controller)then
          call observe_dg_hybrid_inner_residuals(comm,controller,report%residuals,reject_requested,decision_ok,controller_message)
          if(.not.decision_ok)then;message=trim(controller_message);fatal=.true.;return;endif
          if(reject_requested)then
            call reject_dg_hybrid_trial(comm,controller,state,'sustained residual growth',decision_ok,controller_message)
            if(.not.decision_ok)then;message=trim(controller_message);fatal=.true.;return;endif
            rejected_trial=.true.;return
          endif
          call decide_dg_hybrid_stage(comm,controller,state,report,converged,decision_ok,controller_message)
          if(.not.decision_ok)then;message=trim(controller_message);fatal=.true.;return;endif
        else
          call dg_hybrid_stage_tolerances(controls,lambda,tolerances)
          call stage_consensus(report,tolerances,converged)
        endif
        if(converged)return
        call mix_density(iteration,input_density,output_density,controls%density_damping,mixed_density,callback_ok)
        call callback_consensus(callback_ok.and.all(ieee_is_finite(mixed_density)),global_bad,ierr)
        if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='coupled density mixing failed';fatal=.true.;return;endif
        input_density=mixed_density;input_trace=state%trace;state%density=mixed_density
      enddo
    end subroutine converge_current
    subroutine execute_callbacks(inner_iteration,callbacks_ok)
      integer,intent(in)::inner_iteration;logical,intent(out)::callbacks_ok
      callbacks_ok=.false.;state%density=input_density
      previous_operator_epoch=state%operator_epoch
      call build_volume(lambda,input_density,state,callback_ok)
      callback_ok=callback_ok.and.state%operator_epoch>previous_operator_epoch.and.&
        state%operator_structure_fingerprint/=0.and.state%operator_value_fingerprint/=0
      call MPI_Allreduce(state%operator_epoch,minimum_operator_epoch,1,MPI_INTEGER,MPI_MIN,comm,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Allreduce(state%operator_epoch,maximum_operator_epoch,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Allreduce(state%operator_structure_fingerprint,minimum_operator_structure,1,&
        MPI_INTEGER8,MPI_MIN,comm,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Allreduce(state%operator_structure_fingerprint,maximum_operator_structure,1,&
        MPI_INTEGER8,MPI_MAX,comm,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Allreduce(state%operator_value_fingerprint,minimum_operator_value,1,&
        MPI_INTEGER8,MPI_MIN,comm,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Allreduce(state%operator_value_fingerprint,maximum_operator_value,1,&
        MPI_INTEGER8,MPI_MAX,comm,ierr)
      callback_ok=callback_ok.and.ierr==MPI_SUCCESS.and.minimum_operator_epoch==maximum_operator_epoch.and.&
        minimum_operator_structure==maximum_operator_structure.and.minimum_operator_value==maximum_operator_value
      if(frozen_operator_structure==0_int64.and.callback_ok)&
        frozen_operator_structure=state%operator_structure_fingerprint
      callback_ok=callback_ok.and.state%operator_structure_fingerprint==frozen_operator_structure
      if(.not.callback_ok)message='stale or inconsistent operator provenance'
      call callback_consensus(callback_ok,global_bad,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
      call solve_full(lambda,inner_iteration,state,callback_ok);call callback_consensus(callback_ok,global_bad,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
      call refresh_projector(lambda,state,projector_overlap,callback_ok)
      call callback_consensus(callback_ok.and.ieee_is_finite(projector_overlap),global_bad,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
      call refresh_density_trace(lambda,input_density,state,output_density,callback_ok)
      call callback_consensus(callback_ok.and.all(ieee_is_finite(output_density)).and.state%trace_cache_valid,global_bad,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
      call evaluate_residuals(lambda,inner_iteration,input_density,input_trace,state,residuals,projector_overlap,&
        report%electron_ok,report%occupation_ok,report%hermitian_ok,report%symmetry_ok,report%real_space_ok,&
        report%gap_shrinking,callback_ok)
      call callback_consensus(callback_ok,global_bad,ierr)
      callbacks_ok=ierr==MPI_SUCCESS.and.global_bad==0
    end subroutine execute_callbacks
    subroutine fill_report(inner_iteration)
      integer,intent(in)::inner_iteration
      report%residuals=[residuals%r_h,residuals%r_rho,residuals%r_t,residuals%r_s]
      call dg_hybrid_stage_tolerances(controls,lambda,report%tolerances)
      report%projector_overlap=projector_overlap;report%iteration=inner_iteration
      report%finite_ok=all(ieee_is_finite(report%residuals)).and.ieee_is_finite(projector_overlap)
    end subroutine fill_report
    subroutine stage_consensus(candidate,candidate_tolerances,passes)
      type(s_dg_hybrid_stage_report),intent(in)::candidate;real(real64),intent(in)::candidate_tolerances(4)
      logical,intent(out)::passes;integer::local
      passes=candidate%iteration>=1.and.candidate%iteration<=controls%iteration_limit.and.&
        all(candidate%residuals>=0d0).and.all(candidate%residuals<=candidate_tolerances).and.&
        candidate%projector_overlap>=controls%minimum_projector_overlap.and.candidate%electron_ok.and.&
        candidate%occupation_ok.and.candidate%hermitian_ok.and.candidate%symmetry_ok.and.&
        candidate%real_space_ok.and.candidate%finite_ok.and.state%trace_cache_valid
      local=merge(1,0,passes);call MPI_Allreduce(local,global_bad,1,MPI_INTEGER,MPI_MIN,comm,ierr)
      passes=ierr==MPI_SUCCESS.and.global_bad==1
    end subroutine stage_consensus
    subroutine callback_consensus(callback_result,bad,status)
      logical,intent(in)::callback_result;integer,intent(out)::bad,status;integer::local
      local=merge(0,1,callback_result);call MPI_Allreduce(local,bad,1,MPI_INTEGER,MPI_MAX,comm,status)
    end subroutine callback_consensus
#else
    ok=.false.;message='MPI is required for coupled DG continuation SCF';final_lambda=0d0
    accepted_stages=0;rollback_count=0
#endif
  end subroutine run_dg_hybrid_coupled_fixed_points
end module dg_hybrid_continuation_scf
