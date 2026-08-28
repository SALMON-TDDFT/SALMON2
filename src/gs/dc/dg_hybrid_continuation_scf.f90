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
  use dg_hybrid_continuation_acceptance,only:s_dg_hybrid_acceptance_result,validate_dg_hybrid_acceptance_receipt
#ifdef USE_MPI
  use mpi, only: MPI_Allreduce, MPI_BXOR, MPI_IN_PLACE, MPI_INTEGER, MPI_INTEGER8, MPI_MAX, MPI_MIN, MPI_SUCCESS, MPI_SUM
#endif
  implicit none
  private
  type,public::s_dg_hybrid_production_catalog
    logical::frozen=.false.
    integer::global_basis_count=0,global_face_count=0
    integer(int64)::analysis_fingerprint=0_int64,basis_fingerprint=0_int64
    integer(int64)::selection_fingerprint=0_int64,action_fingerprint=0_int64
    integer(int64)::metric_fingerprint=0_int64,face_topology_fingerprint=0_int64
    integer(int64),allocatable::row_ids(:)
    integer,allocatable::effective_wf_ids(:),effective_pw_ids(:)
    integer,allocatable::wf_action(:,:),pw_action(:,:)
    complex(real64),allocatable::metric_rows(:,:),interface_rows(:,:)
  end type s_dg_hybrid_production_catalog
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
    subroutine acceptance_callback(lambda,state,receipt)
      import real64,s_dg_hybrid_trial_state,s_dg_hybrid_acceptance_result
      real(real64),intent(in)::lambda
      type(s_dg_hybrid_trial_state),intent(in)::state
      type(s_dg_hybrid_acceptance_result),intent(out)::receipt
    end subroutine acceptance_callback
  end interface
  type,public::s_dg_hybrid_continuation_callbacks
    procedure(volume_callback),pointer,nopass::build_volume=>null()
    procedure(solve_callback),pointer,nopass::solve_full=>null()
    procedure(projector_callback),pointer,nopass::refresh_projector=>null()
    procedure(density_trace_callback),pointer,nopass::refresh_density_trace=>null()
    procedure(residual_callback),pointer,nopass::evaluate_residuals=>null()
    procedure(mix_callback),pointer,nopass::mix_density=>null()
    procedure(acceptance_callback),pointer,nopass::accept_candidate=>null()
    type(s_dg_hybrid_trial_state)::seed_state
    integer::face_count=0,global_face_count=-1,maximum_inner_iterations=0,accepted_stages=0,rollback_count=0
    integer(int64)::face_topology_fingerprint=0_int64
    real(real64)::final_lambda=0d0
  end type s_dg_hybrid_continuation_callbacks
  public::fingerprint_dg_hybrid_catalog_matrix,validate_dg_hybrid_production_catalog,run_dg_hybrid_continuation_scf
contains
  subroutine fingerprint_dg_hybrid_catalog_matrix(icomm,row_ids,rows,fingerprint,ok,message)
    integer,intent(in)::icomm
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::rows(:,:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,j,ierr,local_bad,global_bad
    integer(int64)::local_hash,entry_hash,real_bits,imaginary_bits
    local_bad=merge(0,1,size(rows,1)==size(row_ids).and.size(rows,2)>0.and.&
      all(row_ids>=1_int64).and.all(row_ids<=int(size(rows,2),int64)).and.&
      all(ieee_is_finite(real(rows))).and.all(ieee_is_finite(aimag(rows))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      fingerprint=0_int64;ok=.false.;message='invalid distributed matrix fingerprint payload';return
    endif
    local_hash=0_int64
    do j=1,size(rows,2)
      do i=1,size(rows,1)
        real_bits=transfer(real(rows(i,j),real64),real_bits)
        imaginary_bits=transfer(aimag(rows(i,j)),imaginary_bits)
        entry_hash=ieor(row_ids(i),ishftc(int(j,int64),11))
        entry_hash=ieor(entry_hash,ishftc(real_bits,23))
        entry_hash=ieor(entry_hash,ishftc(imaginary_bits,41))
        local_hash=ieor(local_hash,entry_hash)
      enddo
    enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,icomm,ierr)
    ok=ierr==MPI_SUCCESS
    if(.not.ok)then;fingerprint=0_int64;message='distributed matrix fingerprint reduction failed';return;endif
    if(fingerprint==0_int64)fingerprint=1907_int64
    message=''
#else
    fingerprint=0_int64;ok=.false.;message='distributed matrix fingerprint requires MPI'
#endif
  end subroutine fingerprint_dg_hybrid_catalog_matrix

  subroutine validate_dg_hybrid_production_catalog(icomm,catalog,ok,message)
    integer,intent(in)::icomm
    type(s_dg_hybrid_production_catalog),intent(in)::catalog
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,local_bad,global_bad,ierr
    integer,allocatable::ownership_count(:)
    integer(int64)::local_fingerprint,minimum_fingerprint,maximum_fingerprint,recomputed_metric_fingerprint
    logical::fingerprint_ok
    character(256)::fingerprint_message

    local_bad=0
    if(.not.catalog%frozen.or.catalog%global_basis_count<1.or.catalog%global_face_count<1)then
      local_bad=1
    elseif(.not.allocated(catalog%row_ids).or..not.allocated(catalog%effective_wf_ids).or.&
        .not.allocated(catalog%effective_pw_ids).or..not.allocated(catalog%wf_action).or.&
        .not.allocated(catalog%pw_action).or..not.allocated(catalog%metric_rows).or.&
        .not.allocated(catalog%interface_rows))then
      local_bad=1
    elseif(any(shape(catalog%metric_rows)/=[size(catalog%row_ids),catalog%global_basis_count]).or.&
        any(shape(catalog%interface_rows)/=shape(catalog%metric_rows)).or.&
        any(catalog%row_ids<1_int64).or.any(catalog%row_ids>int(catalog%global_basis_count,int64)).or.&
        size(catalog%effective_wf_ids)<1.or.size(catalog%effective_pw_ids)<1.or.&
        size(catalog%wf_action,1)/=size(catalog%effective_wf_ids).or.&
        size(catalog%pw_action,1)/=size(catalog%effective_pw_ids).or.&
        size(catalog%wf_action,2)<1.or.size(catalog%pw_action,2)<1.or.&
        catalog%analysis_fingerprint==0_int64.or.catalog%basis_fingerprint==0_int64.or.&
        catalog%selection_fingerprint==0_int64.or.catalog%action_fingerprint==0_int64.or.&
        catalog%metric_fingerprint==0_int64.or.catalog%face_topology_fingerprint==0_int64)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      ok=.false.;message='invalid frozen production continuation catalog';return
    endif
    allocate(ownership_count(catalog%global_basis_count));ownership_count=0
    do i=1,size(catalog%row_ids)
      ownership_count(int(catalog%row_ids(i)))=ownership_count(int(catalog%row_ids(i)))+1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,catalog%global_basis_count,MPI_INTEGER,MPI_SUM,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership_count/=1))then
      ok=.false.;message='production continuation rows are not owned exactly once';return
    endif
    call fingerprint_dg_hybrid_catalog_matrix(icomm,catalog%row_ids,catalog%metric_rows,&
      recomputed_metric_fingerprint,fingerprint_ok,fingerprint_message)
    if(.not.fingerprint_ok.or.recomputed_metric_fingerprint/=catalog%metric_fingerprint)then
      ok=.false.;message='production continuation metric fingerprint mismatch';return
    endif
    local_fingerprint=ieor(catalog%analysis_fingerprint,catalog%basis_fingerprint)
    local_fingerprint=ieor(local_fingerprint,catalog%selection_fingerprint)
    local_fingerprint=ieor(local_fingerprint,catalog%action_fingerprint)
    local_fingerprint=ieor(local_fingerprint,catalog%metric_fingerprint)
    local_fingerprint=ieor(local_fingerprint,catalog%face_topology_fingerprint)
    call MPI_Allreduce(local_fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,icomm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,icomm,ierr)
    ok=ierr==MPI_SUCCESS.and.minimum_fingerprint==maximum_fingerprint
    if(ok)then;message='';else;message='rank-disagreeing production continuation catalog';endif
#else
    ok=.false.;message='production continuation catalog validation requires MPI'
#endif
  end subroutine validate_dg_hybrid_production_catalog

  subroutine run_dg_hybrid_continuation_scf(icomm,state,catalog,controls,callbacks,accepted_state,ok,message)
    integer,intent(in)::icomm;type(s_dg_hybrid_continuation_state),intent(in)::state
    type(s_dg_hybrid_production_catalog),intent(in)::catalog
    type(s_dg_hybrid_controller_controls),intent(in)::controls
    type(s_dg_hybrid_continuation_callbacks),intent(inout)::callbacks
    type(s_dg_hybrid_trial_state),intent(out)::accepted_state
    logical,intent(out)::ok;character(*),intent(out)::message
    logical::callbacks_complete,catalog_ok
    character(256)::catalog_message
#ifdef USE_MPI
    integer::local_bad,global_bad,ierr
#endif
    call validate_dg_hybrid_production_catalog(icomm,catalog,catalog_ok,catalog_message)
    if(.not.catalog_ok)then;ok=.false.;message=trim(catalog_message);return;endif
    callbacks_complete=associated(callbacks%build_volume).and.associated(callbacks%solve_full).and.&
      associated(callbacks%refresh_projector).and.associated(callbacks%refresh_density_trace).and.&
      associated(callbacks%evaluate_residuals).and.associated(callbacks%mix_density).and.&
      associated(callbacks%accept_candidate)
#ifdef USE_MPI
    local_bad=merge(0,1,callbacks_complete)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      ok=.false.;message='collectively incomplete mandatory continuation callback bundle';return
    endif
#else
    if(.not.callbacks_complete)then
      ok=.false.;message='incomplete mandatory continuation callback bundle';return
    endif
#endif
    if(callbacks%face_count/=catalog%global_face_count.or.&
        callbacks%global_face_count/=catalog%global_face_count.or.&
        callbacks%face_topology_fingerprint/=catalog%face_topology_fingerprint)then
      ok=.false.;message='continuation callback metadata differs from frozen production catalog';return
    endif
    call run_dg_hybrid_coupled_fixed_points(icomm,state,controls,callbacks%seed_state,catalog%global_face_count,&
      callbacks%build_volume,callbacks%solve_full,callbacks%refresh_projector,callbacks%refresh_density_trace,&
      callbacks%evaluate_residuals,callbacks%mix_density,callbacks%accept_candidate,&
      callbacks%global_face_count,callbacks%face_topology_fingerprint,callbacks%maximum_inner_iterations,accepted_state,&
      callbacks%final_lambda,callbacks%accepted_stages,callbacks%rollback_count,ok,message)
    if(ok)then
      call validate_dg_hybrid_production_catalog(icomm,catalog,catalog_ok,catalog_message)
      if(.not.catalog_ok)then;ok=.false.;message='production catalog changed during continuation';endif
    endif
  end subroutine run_dg_hybrid_continuation_scf

  subroutine run_dg_hybrid_coupled_fixed_points(icomm,continuation,controls,seed_state,face_count,&
      build_volume,solve_full,refresh_projector,refresh_density_trace,evaluate_residuals,mix_density,&
      accept_candidate,global_face_count,face_topology_fingerprint,maximum_inner_iterations,final_state,final_lambda,&
      accepted_stages,&
      rollback_count,ok,message)
    integer,intent(in)::icomm,face_count,global_face_count,maximum_inner_iterations
    integer(int64),intent(in)::face_topology_fingerprint
    type(s_dg_hybrid_continuation_state),intent(in)::continuation
    type(s_dg_hybrid_controller_controls),intent(in)::controls
    type(s_dg_hybrid_trial_state),intent(in)::seed_state
    procedure(volume_callback)::build_volume
    procedure(solve_callback)::solve_full
    procedure(projector_callback)::refresh_projector
    procedure(density_trace_callback)::refresh_density_trace
    procedure(residual_callback)::evaluate_residuals
    procedure(mix_callback)::mix_density
    procedure(acceptance_callback)::accept_candidate
    type(s_dg_hybrid_trial_state),intent(out)::final_state
    real(real64),intent(out)::final_lambda
    integer,intent(out)::accepted_stages,rollback_count
    logical,intent(out)::ok;character(*),intent(out)::message
#ifdef USE_MPI
    type(s_dg_hybrid_controller)::controller
    type(s_dg_hybrid_trial_state)::state
    type(s_dg_hybrid_trial_state)::input_gamma_state
    type(s_dg_hybrid_residuals)::residuals
    type(s_dg_hybrid_acceptance_result)::acceptance_receipt
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
    call validate_dg_hybrid_controller_contract(icomm,controls,decision_ok,controller_message)
    if(.not.decision_ok)then;message=trim(controller_message);return;endif
    call MPI_Allreduce(maximum_inner_iterations,minimum_iterations,1,MPI_INTEGER,MPI_MIN,icomm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(maximum_inner_iterations,maximum_iterations,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_iterations/=maximum_iterations)then
      message='rank-disagreeing continuation inner iteration limit';return
    endif
    local_bad=merge(0,1,continuation%valid.and.continuation%lambda==0d0.and.&
      allocated(continuation%seed_density).and.allocated(continuation%seed_density_ids).and.&
      continuation%global_density_count>0.and.size(continuation%seed_density)==size(seed_state%density).and.&
      allocated(seed_state%trace).and.seed_state%trace_cache_valid.and.maximum_inner_iterations>0.and.face_count>=0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
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
    call initialize_dg_hybrid_controller(icomm,controls,0d0,state,face_count,controller,decision_ok,controller_message)
    if(.not.decision_ok)then;message=trim(controller_message);return;endif

    do while(controller%accepted_lambda<1d0)
      call propose_dg_hybrid_trial(icomm,controller,state,decision_ok,controller_message)
      if(.not.decision_ok)then;message=trim(controller_message);return;endif
      lambda=controller%trial_lambda;input_density=controller%accepted_state%density
      input_trace=controller%accepted_state%trace;rejected=.false.
      call converge_current(.true.,stage_ok,rejected,fatal_error)
      if(fatal_error)return
      if(rejected)cycle
      if(.not.stage_ok)then
        call reject_dg_hybrid_trial(icomm,controller,state,'inner iteration limit',decision_ok,controller_message)
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
    if(stage_ok)then
      call accept_candidate(lambda,state,acceptance_receipt)
      call validate_dg_hybrid_acceptance_receipt(icomm,acceptance_receipt,size(state%occupations),global_face_count,&
        state%operator_value_fingerprint,face_topology_fingerprint,stage_ok,controller_message)
    endif
    if(.not.stage_ok)then;message='lambda-one fully refreshed residual gate failed';return;endif
    state%density=input_gamma_state%density;state%trace=input_gamma_state%trace
    state%density_epoch=input_gamma_state%density_epoch;state%trace_epoch=input_gamma_state%trace_epoch
    state%derived_epoch=input_gamma_state%derived_epoch;state%trace_cache_valid=.true.
    call accept_candidate(lambda,state,acceptance_receipt)
    call validate_dg_hybrid_acceptance_receipt(icomm,acceptance_receipt,size(state%occupations),global_face_count,&
      state%operator_value_fingerprint,face_topology_fingerprint,stage_ok,controller_message)
    if(.not.stage_ok)then;message='published lambda-one state failed final acceptance oracle';return;endif
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
        call stage_consensus(report,report%tolerances,stage_ok)
        if(stage_ok)then
          call accept_candidate(lambda,state,acceptance_receipt)
          call validate_dg_hybrid_acceptance_receipt(icomm,acceptance_receipt,size(state%occupations),global_face_count,&
            state%operator_value_fingerprint,face_topology_fingerprint,stage_ok,controller_message)
        endif
        if(use_controller)then
          call observe_dg_hybrid_inner_residuals(icomm,controller,report%residuals,reject_requested,decision_ok,controller_message)
          if(.not.decision_ok)then;message=trim(controller_message);fatal=.true.;return;endif
          if(reject_requested)then
            call reject_dg_hybrid_trial(icomm,controller,state,'sustained residual growth',decision_ok,controller_message)
            if(.not.decision_ok)then;message=trim(controller_message);fatal=.true.;return;endif
            rejected_trial=.true.;return
          endif
          if(stage_ok)then
            call decide_dg_hybrid_stage(icomm,controller,state,report,converged,decision_ok,controller_message)
          else
            converged=.false.;decision_ok=.true.;controller_message=''
          endif
          if(.not.decision_ok)then;message=trim(controller_message);fatal=.true.;return;endif
        else
          converged=stage_ok
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
      call MPI_Allreduce(state%operator_epoch,minimum_operator_epoch,1,MPI_INTEGER,MPI_MIN,icomm,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Allreduce(state%operator_epoch,maximum_operator_epoch,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Allreduce(state%operator_structure_fingerprint,minimum_operator_structure,1,&
        MPI_INTEGER8,MPI_MIN,icomm,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Allreduce(state%operator_structure_fingerprint,maximum_operator_structure,1,&
        MPI_INTEGER8,MPI_MAX,icomm,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Allreduce(state%operator_value_fingerprint,minimum_operator_value,1,&
        MPI_INTEGER8,MPI_MIN,icomm,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Allreduce(state%operator_value_fingerprint,maximum_operator_value,1,&
        MPI_INTEGER8,MPI_MAX,icomm,ierr)
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
      local=merge(1,0,passes);call MPI_Allreduce(local,global_bad,1,MPI_INTEGER,MPI_MIN,icomm,ierr)
      passes=ierr==MPI_SUCCESS.and.global_bad==1
    end subroutine stage_consensus
    subroutine callback_consensus(callback_result,bad,status)
      logical,intent(in)::callback_result;integer,intent(out)::bad,status;integer::local
      local=merge(0,1,callback_result);call MPI_Allreduce(local,bad,1,MPI_INTEGER,MPI_MAX,icomm,status)
    end subroutine callback_consensus
#else
    ok=.false.;message='MPI is required for coupled DG continuation SCF';final_lambda=0d0
    accepted_stages=0;rollback_count=0
#endif
  end subroutine run_dg_hybrid_coupled_fixed_points
end module dg_hybrid_continuation_scf
