#include "config.h"
program test_dg_hybrid_continuation_controller_mpi
  use mpi, only: MPI_Allreduce, MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD, MPI_Finalize, MPI_Init, MPI_INTEGER, MPI_MAX, &
    MPI_SUCCESS
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_positive_inf
  use dg_hybrid_continuation_controller
  implicit none
  integer::icomm,id_rank,nproc,ierr,i,solve_count,gap_occupied_index,gap_unoccupied_index
  type(s_dg_hybrid_controller_controls)::controls
  type(s_dg_hybrid_controller)::controller
  type(s_dg_hybrid_controller)::limit_controller
  type(s_dg_hybrid_trial_state)::state,accepted,limit_state
  type(s_dg_hybrid_stage_report)::report
  type(s_dg_hybrid_candidate_acceptance)::candidate_acceptance
  real(real64)::t0(4),t1(4),lambda_before,step_before
  logical::ok,accept,meaningful_gap
  character(256)::message

  call MPI_Init(ierr);icomm=MPI_COMM_WORLD
  call MPI_Comm_rank(icomm,id_rank,ierr);call MPI_Comm_size(icomm,nproc,ierr)
  call fill_state(state,9)
  call reject_dg_hybrid_trial(icomm,controller,state,'early nonconvergence',ok,message)
  call require(.not.ok.and.trim(message)=='invalid continuation rollback request',&
    'uninitialized continuation rollback did not fail collectively with the named diagnostic')
  call default_dg_hybrid_controller_controls(controls)
  call require(controls%initial_step==0.125d0.and.controls%minimum_step==0.015625d0.and.&
    controls%maximum_step==0.5d0.and.controls%growth_factor==1.5d0.and.controls%shrink_factor==0.5d0,&
    'adaptive-step defaults changed')
  call require(controls%residual_growth_limit==4d0.and.controls%density_damping==0.5d0.and.&
    controls%minimum_projector_overlap==0.9d0.and.controls%maximum_rollbacks==8,&
    'continuation safety defaults changed')
  controls%iteration_limit=8
  controls%intermediate_tolerance=[1d-4,2d-4,3d-4,4d-4]
  controls%final_tolerance=[1d-8,2d-8,3d-8,4d-8]
  call dg_hybrid_stage_tolerances(controls,0d0,t0)
  call dg_hybrid_stage_tolerances(controls,1d0,t1)
  call require(all(t0==controls%intermediate_tolerance).and.all(t1==controls%final_tolerance),&
    'endpoint continuation tolerances are incorrect')
  call dg_hybrid_stage_tolerances(controls,0.25d0,t0);call dg_hybrid_stage_tolerances(controls,0.75d0,t1)
  call require(all(t1<=t0).and.all(t1>=controls%final_tolerance),'inexact tolerances are not monotone')
  call dg_hybrid_continuation_state_count([2d0,0d0],2,2,solve_count,meaningful_gap,&
    gap_occupied_index,gap_unoccupied_index,ok)
  call require(ok.and.solve_count==2.and.meaningful_gap.and.gap_occupied_index==1.and.gap_unoccupied_index==2,&
    'configured occupied-empty pair did not provide its available gap')
  call dg_hybrid_continuation_state_count([2d0,0.5d0],3,2,solve_count,meaningful_gap,&
    gap_occupied_index,gap_unoccupied_index,ok)
  call require(ok.and.solve_count==3.and.meaningful_gap.and.gap_occupied_index==2.and.gap_unoccupied_index==3,&
    'fractionally occupied boundary did not retain an available separation diagnostic')
  call dg_hybrid_continuation_state_count([1d0,1d0],2,2,solve_count,meaningful_gap,&
    gap_occupied_index,gap_unoccupied_index,ok)
  call require(ok.and.solve_count==2.and..not.meaningful_gap.and.gap_occupied_index==0.and.gap_unoccupied_index==0,&
    'degenerate fully retained occupied space was rejected without an extra state')
  call dg_hybrid_continuation_state_count([2d0,2d0,0d0],7,4,solve_count,meaningful_gap,&
    gap_occupied_index,gap_unoccupied_index,ok)
  call require(ok.and.solve_count==7.and.meaningful_gap.and.gap_occupied_index==2.and.&
    gap_unoccupied_index==3,'continuation did not request one complete LCFO eigensystem')
  call dg_hybrid_continuation_state_count([2d0,2d0,2d0],7,2,solve_count,meaningful_gap,&
    gap_occupied_index,gap_unoccupied_index,ok)
  call require(.not.ok,'symmetry target smaller than the occupied window was accepted')
  call dg_hybrid_continuation_state_count([2d0,2d0,2d0],11,7,solve_count,meaningful_gap,&
    gap_occupied_index,gap_unoccupied_index,ok)
  call require(ok.and.solve_count==11.and.gap_occupied_index==3.and.gap_unoccupied_index==4,&
    'complete solve was replaced by a requested-rank prefix')

  call initialize_dg_hybrid_candidate_acceptance(icomm,7,0.4d0,candidate_acceptance,ok,message)
  call require(ok.and.candidate_acceptance%construction_rank==7.and.&
    .not.candidate_acceptance%legacy_dynamic_rank,trim(message))
  if(nproc>1)then
    call record_dg_hybrid_complete_lcfo_solve(icomm,candidate_acceptance,7,1,&
      merge(1001_int64,1099_int64,id_rank==0),ok,message)
    call require(.not.ok.and.candidate_acceptance%phase==1,&
      'rank-disagreeing complete eigensystem fingerprint was accepted')
  endif
  call record_dg_hybrid_occupation_policy(icomm,candidate_acceptance,2,.true.,1002_int64,ok,message)
  call require(.not.ok,'occupation policy was recorded before the complete LCFO solve')
  call record_dg_hybrid_complete_lcfo_solve(icomm,candidate_acceptance,6,1,1001_int64,ok,message)
  call require(.not.ok,'a prefix eigensystem was accepted as the complete LCFO solve')
  call record_dg_hybrid_complete_lcfo_solve(icomm,candidate_acceptance,7,2,1001_int64,ok,message)
  call require(.not.ok,'two complete eigensolver calls were accepted for one final candidate')
  call record_dg_hybrid_complete_lcfo_solve(icomm,candidate_acceptance,7,1,1001_int64,ok,message)
  call require(ok.and.candidate_acceptance%solved_rank==7,trim(message))
  call record_dg_hybrid_occupation_policy(icomm,candidate_acceptance,2,.false.,1002_int64,ok,message)
  call require(.not.ok,'an electron-count failure was accepted as the Task 9 occupation policy')
  call record_dg_hybrid_occupation_policy(icomm,candidate_acceptance,2,.true.,1002_int64,ok,message)
  call require(ok.and.candidate_acceptance%occupied_rank==2,trim(message))
  call record_dg_hybrid_unconditional_gates(icomm,candidate_acceptance,.false.,.true.,1d-2,1d-12,1d-10,ok,message)
  call require(.not.ok,'occupied-projector failure was accepted before spectral extension')
  call record_dg_hybrid_spectral_certification(icomm,candidate_acceptance,3,4,4,.true.,.false.,.false.,&
    1003_int64,ok,message)
  call require(.not.ok,'empty-state extension repaired a failed occupied-projector gate')
  call record_dg_hybrid_unconditional_gates(icomm,candidate_acceptance,.true.,.false.,1d-12,1d-2,1d-10,ok,message)
  call require(.not.ok,'density failure was accepted before spectral extension')
  if(nproc>1)then
    call record_dg_hybrid_unconditional_gates(icomm,candidate_acceptance,.true.,.true.,&
      1d-12+real(id_rank,real64)*epsilon(1d0),2d-12,1d-10,ok,message)
    call require(.not.ok.and.candidate_acceptance%phase==3,&
      'rank-disagreeing occupied-projector defect was accepted')
    call record_dg_hybrid_unconditional_gates(icomm,candidate_acceptance,.true.,.true.,&
      1d-12,2d-12,1d-10+real(id_rank,real64)*epsilon(1d0),ok,message)
    call require(.not.ok.and.candidate_acceptance%phase==3,&
      'rank-disagreeing physical-gate tolerance was accepted')
  endif
  call record_dg_hybrid_unconditional_gates(icomm,candidate_acceptance,.true.,.true.,1d-12,2d-12,1d-10,ok,message)
  call require(ok.and.candidate_acceptance%occupied_gate.and.candidate_acceptance%density_gate,trim(message))
  call authorize_dg_hybrid_v4_publication(icomm,candidate_acceptance,4,4,.true.,ok,message)
  call require(.not.ok,'v4 publication was authorized before certification and localization')
  call record_dg_hybrid_spectral_certification(icomm,candidate_acceptance,3,4,4,.false.,.false.,.false.,&
    1003_int64,ok,message)
  call require(.not.ok,'explicit energy window was accepted without a proof state')
  call record_dg_hybrid_spectral_certification(icomm,candidate_acceptance,3,4,4,.true.,.false.,.false.,&
    1003_int64,ok,message)
  call require(ok.and.candidate_acceptance%certified_rank==4,trim(message))
  call record_dg_hybrid_certified_rt_basis(icomm,candidate_acceptance,7,1004_int64,1005_int64,ok,message)
  call require(ok.and.candidate_acceptance%rt_basis_rank==7.and.candidate_acceptance%certified_rank==4,&
    'localized v4 construction basis was not retained after low-energy certification')
  call record_dg_hybrid_certified_rt_basis(icomm,candidate_acceptance,4,1004_int64,1005_int64,ok,message)
  call require(.not.ok,'spectral prefix was accepted as the localized v4 propagation basis')
  call authorize_dg_hybrid_v4_publication(icomm,candidate_acceptance,3,7,.true.,ok,message)
  call require(.not.ok,'legacy checkpoint version 2 was authorized for certified RT publication')
  call authorize_dg_hybrid_v4_publication(icomm,candidate_acceptance,4,4,.true.,ok,message)
  call require(.not.ok,'spectral prefix was accepted as the localized v4 payload rank')
  call authorize_dg_hybrid_v4_publication(icomm,candidate_acceptance,4,7,.true.,ok,message)
  call require(ok.and.candidate_acceptance%published_rt_rank==7.and.&
    candidate_acceptance%certified_rank<candidate_acceptance%published_rt_rank,trim(message))

  call initialize_dg_hybrid_candidate_acceptance(icomm,5,-1d0,candidate_acceptance,ok,message)
  call require(ok.and.candidate_acceptance%legacy_dynamic_rank.and.&
    candidate_acceptance%legacy_warning_required,trim(message))
  call record_dg_hybrid_complete_lcfo_solve(icomm,candidate_acceptance,5,1,2001_int64,ok,message)
  call require(ok,trim(message))
  call record_dg_hybrid_occupation_policy(icomm,candidate_acceptance,2,.true.,2002_int64,ok,message)
  call require(ok,trim(message))
  call record_dg_hybrid_unconditional_gates(icomm,candidate_acceptance,.true.,.true.,1d-12,2d-12,1d-10,ok,message)
  call require(ok,trim(message))
  call record_dg_hybrid_spectral_certification(icomm,candidate_acceptance,3,5,5,.false.,.true.,.false.,&
    2003_int64,ok,message)
  call require(.not.ok,'legacy dynamic-rank certification omitted its explicit warning receipt')
  call record_dg_hybrid_spectral_certification(icomm,candidate_acceptance,3,5,5,.false.,.true.,.true.,&
    2003_int64,ok,message)
  call require(ok.and.candidate_acceptance%legacy_warning_observed,trim(message))
  call record_dg_hybrid_certified_rt_basis(icomm,candidate_acceptance,5,2004_int64,2005_int64,ok,message)
  call require(ok,trim(message))
  call authorize_dg_hybrid_v4_publication(icomm,candidate_acceptance,4,5,.true.,ok,message)
  call require(ok.and.candidate_acceptance%published_rt_rank==5,trim(message))
  call initialize_dg_hybrid_candidate_acceptance(icomm,5,-1d0-epsilon(1d0),candidate_acceptance,ok,message)
  call require(.not.ok,'a negative energy window other than exactly -1 was accepted')

  call fill_state(accepted,10)
  if(nproc>1)then
    if(id_rank==0)controls%growth_factor=1.6d0
    call initialize_dg_hybrid_controller(icomm,controls,0d0,accepted,3,controller,ok,message)
    call require(.not.ok,'rank-disagreeing controller controls were accepted')
    controls%growth_factor=1.5d0
    call initialize_dg_hybrid_controller(icomm,controls,0.01d0*id_rank,accepted,3+id_rank,controller,ok,message)
    call require(.not.ok,'rank-disagreeing accepted lambda was accepted')
    call initialize_dg_hybrid_controller(icomm,controls,0d0,accepted,3+id_rank,controller,ok,message)
    call require(ok.and.size(controller%face_lambda)==3+id_rank,&
      'rank-local canonical-face ownership was rejected')
    call initialize_dg_hybrid_controller(icomm,controls,0d0,accepted,merge(0,1,id_rank==0),controller,ok,message)
    call require(ok,'rank with no locally owned canonical faces was rejected')
    call require(size(controller%face_lambda)==merge(0,1,id_rank==0),&
      'rank-local zero-length face lambda allocation is incorrect')
  endif
  controls%growth_factor=ieee_value(1d0,ieee_positive_inf)
  call initialize_dg_hybrid_controller(icomm,controls,0d0,accepted,3,controller,ok,message)
  call require(.not.ok,'nonfinite controller controls were accepted')
  controls%growth_factor=1.5d0
  call initialize_dg_hybrid_controller(icomm,controls,0d0,accepted,3,controller,ok,message)
  call require(ok,trim(message))
  call require(controller%trace_valid.and.controller%accepted_state%trace_cache_valid,&
    'valid accepted trace checkpoint was invalidated at initialization')
  state=accepted
  call propose_dg_hybrid_trial(icomm,controller,state,ok,message)
  call require(ok.and.controller%trial_lambda==0.125d0,'initial lambda proposal is incorrect')
  call require(all(controller%face_lambda==controller%trial_lambda),'face-local lambda was proposed')
  call require(.not.controller%trace_valid,'trace cache survived a lambda state change')
  call require(.not.state%trace_cache_valid,'trial-state trace cache remained usable after lambda change')

  call passing_report(controller,report);report%residuals(3)=2d0*report%tolerances(3)
  call decide_dg_hybrid_stage(icomm,controller,state,report,accept,ok,message)
  call require(ok.and..not.accept,'stage passed while one residual channel failed')
  call passing_report(controller,report);report%residuals=2d0*report%tolerances
  report%tolerances=100d0*report%tolerances
  call decide_dg_hybrid_stage(icomm,controller,state,report,accept,ok,message)
  call require(ok.and..not.accept,'caller-supplied loose tolerances bypassed controller controls')
  call passing_report(controller,report);report%residuals(1)=-1d0
  call decide_dg_hybrid_stage(icomm,controller,state,report,accept,ok,message)
  call require(ok.and..not.accept,'negative residual was accepted')
  call passing_report(controller,report)
  call decide_dg_hybrid_stage(icomm,controller,state,report,accept,ok,message)
  call require(ok.and..not.accept,'stage accepted a stale interface trace cache')
  call passing_report(controller,report);report%iteration=3
  state%trace_cache_valid=.true.
  step_before=controller%step
  call decide_dg_hybrid_stage(icomm,controller,state,report,accept,ok,message)
  call require(ok.and.accept,'fully converged easy stage was rejected')
  call require(controller%step==min(controls%maximum_step,step_before*controls%growth_factor),&
    'easy accepted stage did not grow the next step')

  accepted=state;call propose_dg_hybrid_trial(icomm,controller,state,ok,message);call require(ok,trim(message))
  call mutate_state(state)
  step_before=controller%step;lambda_before=controller%accepted_lambda
  call reject_dg_hybrid_trial(icomm,controller,state,'forced rollback',ok,message)
  call require(ok,'forced rollback failed: '//trim(message))
  call require(equal_state(state,accepted),'rollback did not restore every accepted payload bit')
  call require(controller%accepted_lambda==lambda_before.and.controller%step==max(controls%minimum_step,&
    step_before*controls%shrink_factor),'rollback did not restore lambda and reduce its step')
  call require(controller%trace_valid.and.state%trace_cache_valid,&
    'rollback did not restore the valid accepted trace checkpoint')

  call propose_dg_hybrid_trial(icomm,controller,state,ok,message);call require(ok,trim(message))
  call passing_report(controller,report);report%iteration=1
  state%trace_cache_valid=.true.
  step_before=controller%step
  call decide_dg_hybrid_stage(icomm,controller,state,report,accept,ok,message)
  call require(ok.and.accept.and.controller%step==step_before,'fast retry grew the step after rollback')

  call propose_dg_hybrid_trial(icomm,controller,state,ok,message);call require(ok,trim(message))
  call observe_dg_hybrid_inner_residuals(icomm,controller,[1d-6,1d-6,1d-6,1d-6],accept,ok,message)
  call require(ok.and..not.accept,'first inner residual sample rejected a trial')
  if(nproc>1)then
    if(id_rank==0)then
      t0=[5d-6,5d-6,5d-6,5d-6]
    else;t0=[1d-6,1d-6,1d-6,1d-6]
    endif
  else;t0=[5d-6,5d-6,5d-6,5d-6]
  endif
  call observe_dg_hybrid_inner_residuals(icomm,controller,t0,accept,ok,message)
  call require(ok.and..not.accept,'one residual-growth event rejected a trial')
  if(nproc>1)then
    if(id_rank==1)then
      t0=[2.5d-5,2.5d-5,2.5d-5,2.5d-5]
    else;t0=[1d-6,1d-6,1d-6,1d-6]
    endif
  else;t0=[3d-5,3d-5,3d-5,3d-5]
  endif
  call observe_dg_hybrid_inner_residuals(icomm,controller,t0,accept,ok,message)
  call require(ok.and.accept,'two consecutive excessive growth events did not request rollback')

  call passing_report(controller,report);report%gap_shrinking=id_rank==0;report%iteration=min(4,2+id_rank)
  state%trace_cache_valid=.true.
  step_before=controller%step
  call decide_dg_hybrid_stage(icomm,controller,state,report,accept,ok,message)
  call require(ok.and.accept,'shrinking gap alone rejected an acceptable stage')
  call require(controller%step==max(controls%minimum_step,step_before*controls%shrink_factor),&
    'shrinking gap did not conservatively reduce the next step')
  call propose_dg_hybrid_trial(icomm,controller,state,ok,message);call require(ok,trim(message))
  state%occupations=[1.5d0,0.5d0];state%eigenvalues=[-0.25d0,-0.25d0];state%trace_cache_valid=.true.
  call passing_report(controller,report);report%gap_shrinking=.false.
  call decide_dg_hybrid_stage(icomm,controller,state,report,accept,ok,message)
  call require(ok.and.accept,'continuous fractionally occupied zero-gap crossing was rejected')
  call propose_dg_hybrid_trial(icomm,controller,state,ok,message);call require(ok,trim(message))
  call passing_report(controller,report);report%occupation_ok=.false.
  call decide_dg_hybrid_stage(icomm,controller,state,report,accept,ok,message)
  call require(ok.and..not.accept,'failed cluster-aware occupation was accepted')
  call passing_report(controller,report);report%projector_overlap=0.89d0
  call decide_dg_hybrid_stage(icomm,controller,state,report,accept,ok,message)
  call require(ok.and..not.accept,'occupied-projector discontinuity was accepted')
  call passing_report(controller,report);report%symmetry_ok=.false.
  call decide_dg_hybrid_stage(icomm,controller,state,report,accept,ok,message)
  call require(ok.and..not.accept,'symmetry failure was accepted')

  call initialize_dg_hybrid_controller(icomm,controls,0.875d0,accepted,3,limit_controller,ok,message)
  call require(ok,trim(message));limit_state=accepted
  call propose_dg_hybrid_trial(icomm,limit_controller,limit_state,ok,message)
  call require(ok.and.limit_controller%trial_lambda==1d0,'final lambda proposal is incorrect')
  limit_state%trace_cache_valid=.true.
  call passing_report(limit_controller,report);report%iteration=controls%iteration_limit+2
  call decide_dg_hybrid_stage(icomm,limit_controller,limit_state,report,accept,ok,message)
  call require(ok.and..not.accept,'more than one final refresh iteration was accepted')
  report%iteration=controls%iteration_limit+1
  call decide_dg_hybrid_stage(icomm,limit_controller,limit_state,report,accept,ok,message)
  call require(ok.and.accept,'the final refresh outside the ordinary iteration budget was rejected')

  call initialize_dg_hybrid_controller(icomm,controls,0d0,accepted,3,limit_controller,ok,message)
  call require(ok,trim(message));limit_state=accepted
  do i=1,controls%maximum_rollbacks
    call propose_dg_hybrid_trial(icomm,limit_controller,limit_state,ok,message);call require(ok,trim(message))
    call reject_dg_hybrid_trial(icomm,limit_controller,limit_state,'rollback-limit fixture',ok,message)
    call require(ok,'rollback was rejected before the configured limit')
  enddo
  call require(limit_controller%step==controls%minimum_step,'rollback step fell below or stopped above its minimum')
  call propose_dg_hybrid_trial(icomm,limit_controller,limit_state,ok,message);call require(ok,trim(message))
  call reject_dg_hybrid_trial(icomm,limit_controller,limit_state,'ninth rollback',ok,message)
  call require(.not.ok,'rollback beyond the configured limit was accepted')
  if(id_rank==0)write(*,'(a,i0,a)')'PASS hybrid continuation controller on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine fill_state(value,seed)
    type(s_dg_hybrid_trial_state),intent(out)::value;integer,intent(in)::seed
    allocate(value%density(3),value%potential(3),value%projector(2,2),value%trace(2,2),&
      value%occupations(2),value%eigenvalues(2),value%mixing_history(4))
    value%density=[(real(seed+i,real64),i=1,3)];value%potential=value%density+10d0
    value%projector=cmplx(reshape([(real(seed+i,real64),i=1,4)],[2,2]),0d0,real64)
    value%trace=2d0*value%projector;value%occupations=[2d0,2d0];value%eigenvalues=[-1d0,-0.5d0]
    value%mixing_history=[(real(seed+20+i,real64),i=1,4)]
    value%density_epoch=11;value%operator_epoch=12;value%projector_epoch=13;value%trace_epoch=14
    value%derived_epoch=15
    value%operator_structure_fingerprint=701_int64;value%operator_value_fingerprint=702_int64
    value%trace_cache_valid=.true.
  end subroutine fill_state
  subroutine mutate_state(value)
    type(s_dg_hybrid_trial_state),intent(inout)::value
    value%density=value%density+1d0;value%potential=value%potential+2d0
    value%projector=value%projector+(3d0,1d0);value%trace=value%trace+(4d0,-1d0)
    value%occupations=value%occupations/2d0;value%eigenvalues=value%eigenvalues+5d0
    value%mixing_history=-value%mixing_history
    value%density_epoch=101;value%operator_epoch=102;value%projector_epoch=103;value%trace_epoch=104
    value%derived_epoch=105
    value%operator_structure_fingerprint=801_int64;value%operator_value_fingerprint=802_int64
    value%trace_cache_valid=.true.
  end subroutine mutate_state
  logical function equal_state(a,b)
    type(s_dg_hybrid_trial_state),intent(in)::a,b
    equal_state=all(a%density==b%density).and.all(a%potential==b%potential).and.&
      all(a%projector==b%projector).and.all(a%trace==b%trace).and.all(a%occupations==b%occupations).and.&
      all(a%eigenvalues==b%eigenvalues).and.all(a%mixing_history==b%mixing_history).and.&
      a%density_epoch==b%density_epoch.and.a%operator_epoch==b%operator_epoch.and.&
      a%projector_epoch==b%projector_epoch.and.a%trace_epoch==b%trace_epoch.and.a%derived_epoch==b%derived_epoch.and.&
      a%operator_structure_fingerprint==b%operator_structure_fingerprint.and.&
      a%operator_value_fingerprint==b%operator_value_fingerprint.and.&
      (a%trace_cache_valid.eqv.b%trace_cache_valid)
  end function equal_state
  subroutine passing_report(ctrl,value)
    type(s_dg_hybrid_controller),intent(in)::ctrl;type(s_dg_hybrid_stage_report),intent(out)::value
    call dg_hybrid_stage_tolerances(ctrl%controls,ctrl%trial_lambda,value%tolerances)
    value%residuals=0.5d0*value%tolerances;value%projector_overlap=0.99d0;value%electron_ok=.true.
    value%occupation_ok=.true.;value%hermitian_ok=.true.;value%symmetry_ok=.true.;value%real_space_ok=.true.
    value%finite_ok=.true.;value%gap_shrinking=.false.;value%iteration=ctrl%controls%iteration_limit
  end subroutine passing_report
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label;integer::local_bad,global_bad
    local_bad=merge(0,1,condition);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;if(id_rank==0)write(0,'(a)')trim(label);error stop 1;endif
  end subroutine require
end program test_dg_hybrid_continuation_controller_mpi
