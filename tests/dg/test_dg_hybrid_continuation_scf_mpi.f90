#include "config.h"
program test_dg_hybrid_continuation_scf_mpi
  use mpi, only: MPI_Allreduce, MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_PRECISION, &
    MPI_Finalize, MPI_IN_PLACE, MPI_Init, MPI_INTEGER, MPI_MAX, MPI_SUCCESS, MPI_SUM
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_continuation_state,only:s_dg_hybrid_continuation_state
  use dg_hybrid_continuation_controller,only:s_dg_hybrid_controller_controls,s_dg_hybrid_trial_state,&
    default_dg_hybrid_controller_controls
  use dg_hybrid_continuation_residuals,only:s_dg_hybrid_residuals
  use dg_hybrid_continuation_acceptance,only:s_dg_hybrid_acceptance_result
  use dg_hybrid_continuation_scf,only:s_dg_hybrid_continuation_callbacks,run_dg_hybrid_continuation_scf
  implicit none
  integer,parameter::nglobal=4
  integer::icomm,id_rank,nproc,ierr,i,nlocal,position,phase,solve_count,accepted_stages,rollbacks,acceptance_count
  integer(int64),allocatable::ids(:)
  real(real64),allocatable::dc_density(:),last_built_density(:),last_accepted_density(:)
  real(real64)::final_lambda
  complex(real64)::last_built_trace
  integer::last_built_density_epoch,last_built_trace_epoch
  type(s_dg_hybrid_continuation_state)::continuation
  type(s_dg_hybrid_controller_controls)::controls
  type(s_dg_hybrid_trial_state)::seed,final_state
  type(s_dg_hybrid_continuation_callbacks)::callbacks
  logical::ok,first_volume,lambda_zero_passed,forced_growth_complete,poison_final,lambda_one_converged,fatal_positive,&
    lambda_zero_gate_delayed,stale_operator_positive,rank_divergent_operator,fail_symmetry_oracle,fail_grid_oracle
  character(256)::message

  call MPI_Init(ierr);icomm=MPI_COMM_WORLD
  call MPI_Comm_rank(icomm,id_rank,ierr);call MPI_Comm_size(icomm,nproc,ierr)
  nlocal=count([(mod(i-1,nproc)==id_rank,i=1,nglobal)])
  allocate(ids(nlocal),dc_density(nlocal),last_built_density(nlocal),last_accepted_density(nlocal));position=0
  do i=1,nglobal
    if(mod(i-1,nproc)/=id_rank)cycle
    position=position+1;ids(position)=i;dc_density(position)=0.15d0+0.01d0*i
  enddo
  continuation%valid=.true.;continuation%lambda=0d0;continuation%global_density_count=nglobal
  allocate(continuation%seed_density_ids(nlocal),continuation%seed_density(nlocal),continuation%mixed_density(nlocal))
  continuation%seed_density_ids=ids;continuation%seed_density=dc_density;continuation%mixed_density=dc_density
  call fill_state(seed)
  call default_dg_hybrid_controller_controls(controls)
  controls%intermediate_tolerance=[2d-5,2d-5,2d-5,2d-8]
  controls%final_tolerance=[2d-8,2d-8,2d-8,2d-10]
  controls%iteration_limit=80
  fail_symmetry_oracle=.false.;fail_grid_oracle=.false.;acceptance_count=0
  call configure_callbacks()
  callbacks%accept_candidate=>null();solve_count=0
  call run_dg_hybrid_continuation_scf(icomm,continuation,controls,callbacks,final_state,ok,message)
  call require(.not.ok.and.solve_count==0.and.index(message,'collective')>0,&
    'missing acceptance oracle was not rejected collectively before execution')
  callbacks%accept_candidate=>acceptance_oracle
  callbacks%build_volume=>null();solve_count=0
  call run_dg_hybrid_continuation_scf(icomm,continuation,controls,callbacks,final_state,ok,message)
  call require(.not.ok.and.solve_count==0.and.index(message,'collective')>0,&
    'incomplete callback bundle was not rejected collectively before execution')
  callbacks%build_volume=>volume_build
  if(nproc>1)then
    if(id_rank==0)controls%density_damping=0.4d0
    phase=0;solve_count=0;first_volume=.true.;lambda_zero_passed=.false.;forced_growth_complete=.true.
    poison_final=.false.;lambda_one_converged=.false.;fatal_positive=.false.;lambda_zero_gate_delayed=.false.
    stale_operator_positive=.false.
    callbacks%maximum_inner_iterations=80+id_rank
    call run_dg_hybrid_continuation_scf(icomm,continuation,controls,callbacks,final_state,ok,message)
    call require(.not.ok.and.solve_count==0,'rank-disagreeing controls reached the lambda-zero callbacks')
    controls%density_damping=0.5d0
    solve_count=0
    call run_dg_hybrid_continuation_scf(icomm,continuation,controls,callbacks,final_state,ok,message)
    call require(.not.ok.and.solve_count==0.and.index(message,'iteration limit')>0,&
      'rank-disagreeing inner iteration limit reached the lambda-zero callbacks')
  endif
  callbacks%maximum_inner_iterations=80
  phase=0;solve_count=0;first_volume=.true.;lambda_zero_passed=.false.;forced_growth_complete=.false.
  poison_final=.false.;lambda_one_converged=.false.;fatal_positive=.false.;lambda_zero_gate_delayed=.false.
  stale_operator_positive=.false.;rank_divergent_operator=.false.
  callbacks%seed_state=seed
  call run_dg_hybrid_continuation_scf(icomm,continuation,controls,callbacks,final_state,ok,message)
  final_lambda=callbacks%final_lambda;accepted_stages=callbacks%accepted_stages;rollbacks=callbacks%rollback_count
  call require(ok,trim(message))
  call require(abs(final_lambda-1d0)<1d-15.and.rollbacks>=1.and.accepted_stages>=2,&
    'continuation did not reach lambda one through accepted stages and rollback')
  call require(maxval(abs(final_state%density-fixed_density(1d0)))<5d-8,&
    'lambda-one density differs from the dense fixed-point reference')
  call require(final_state%trace_cache_valid.and.phase==5,'lambda-one state was not fully refreshed without mixing')
  call require(all(final_state%density==last_built_density),&
    'published final density does not match the final Hamiltonian build provenance')
  call require(final_state%trace(1,1)==last_built_trace.and.&
    final_state%density_epoch==last_built_density_epoch.and.final_state%trace_epoch==last_built_trace_epoch,&
    'published density and trace do not share the final operator-input Gamma provenance')
  call require(all(final_state%density==last_accepted_density),&
    'published final state differs from the state checked by the last acceptance oracle')
  call fill_state(seed);phase=0;solve_count=0;first_volume=.true.;lambda_zero_passed=.false.
  forced_growth_complete=.true.;poison_final=.true.;lambda_one_converged=.false.;lambda_zero_gate_delayed=.false.
  callbacks%seed_state=seed
  call run_dg_hybrid_continuation_scf(icomm,continuation,controls,callbacks,final_state,ok,message)
  call require(.not.ok.and.index(message,'lambda-one fully refreshed residual gate failed')>0,&
    'lambda-one refresh used stale intermediate tolerances')
  call fill_state(seed);phase=0;solve_count=0;first_volume=.true.;lambda_zero_passed=.false.
  forced_growth_complete=.true.;poison_final=.false.;lambda_one_converged=.false.;fatal_positive=.true.
  lambda_zero_gate_delayed=.false.
  callbacks%seed_state=seed
  call run_dg_hybrid_continuation_scf(icomm,continuation,controls,callbacks,final_state,ok,message)
  call require(.not.ok.and.index(message,'coupled continuation callback failed')>0,&
    'fatal positive-lambda callback failure was retried as an iteration rejection')
  call fill_state(seed);callbacks%seed_state=seed;phase=0;solve_count=0;first_volume=.true.;lambda_zero_passed=.false.
  forced_growth_complete=.true.;poison_final=.false.;lambda_one_converged=.false.;fatal_positive=.false.
  lambda_zero_gate_delayed=.false.;stale_operator_positive=.true.
  call run_dg_hybrid_continuation_scf(icomm,continuation,controls,callbacks,final_state,ok,message)
  call require(.not.ok.and.index(message,'operator provenance')>0,&
    'stale density-dependent operator provenance was accepted')
  if(nproc>1)then
    call fill_state(seed);callbacks%seed_state=seed;phase=0;solve_count=0;first_volume=.true.;lambda_zero_passed=.false.
    forced_growth_complete=.true.;poison_final=.false.;lambda_one_converged=.false.;fatal_positive=.false.
    lambda_zero_gate_delayed=.false.;stale_operator_positive=.false.;rank_divergent_operator=.true.
    call run_dg_hybrid_continuation_scf(icomm,continuation,controls,callbacks,final_state,ok,message)
    call require(.not.ok.and.index(message,'operator provenance')>0,&
      'rank-divergent operator provenance was accepted')
  endif
  call fill_state(seed);callbacks%seed_state=seed;phase=0;solve_count=0;acceptance_count=0
  first_volume=.true.;lambda_zero_passed=.false.;forced_growth_complete=.true.;poison_final=.false.
  lambda_one_converged=.false.;fatal_positive=.false.;lambda_zero_gate_delayed=.false.
  stale_operator_positive=.false.;rank_divergent_operator=.false.;fail_symmetry_oracle=.true.;fail_grid_oracle=.false.
  call run_dg_hybrid_continuation_scf(icomm,continuation,controls,callbacks,final_state,ok,message)
  call require(.not.ok.and.acceptance_count>0,'candidate stage bypassed the symmetry acceptance oracle')
  call fill_state(seed);callbacks%seed_state=seed;phase=0;solve_count=0;acceptance_count=0
  first_volume=.true.;lambda_zero_passed=.false.;forced_growth_complete=.true.;poison_final=.false.
  lambda_one_converged=.false.;fatal_positive=.false.;lambda_zero_gate_delayed=.false.
  stale_operator_positive=.false.;rank_divergent_operator=.false.;fail_symmetry_oracle=.false.;fail_grid_oracle=.true.
  call run_dg_hybrid_continuation_scf(icomm,continuation,controls,callbacks,final_state,ok,message)
  call require(.not.ok.and.acceptance_count>0,'candidate stage bypassed the reconstructed-grid acceptance oracle')
  if(id_rank==0)write(*,'(a,i0,a)')'PASS hybrid continuation SCF on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine configure_callbacks()
    callbacks%build_volume=>volume_build;callbacks%solve_full=>full_solve
    callbacks%refresh_projector=>projector_refresh;callbacks%refresh_density_trace=>density_trace_refresh
    callbacks%evaluate_residuals=>residual_evaluation;callbacks%mix_density=>density_mix
    callbacks%accept_candidate=>acceptance_oracle
    callbacks%seed_state=seed;callbacks%face_count=1;callbacks%global_face_count=1
    callbacks%face_topology_fingerprint=96_int64
    callbacks%maximum_inner_iterations=80
  end subroutine configure_callbacks
  subroutine acceptance_oracle(lambda,state,receipt)
    real(real64),intent(in)::lambda
    type(s_dg_hybrid_trial_state),intent(in)::state
    type(s_dg_hybrid_acceptance_result),intent(out)::receipt
    acceptance_count=acceptance_count+1
    last_accepted_density=state%density
    receipt=s_dg_hybrid_acceptance_result()
    receipt%valid=lambda>=0d0.and.state%trace_cache_valid
    receipt%symmetry_complete=.not.fail_symmetry_oracle;receipt%grid_complete=.not.fail_grid_oracle
    receipt%face_complete=.true.;receipt%expected_occupied_count=2;receipt%checked_occupied_count=2
    receipt%checked_face_count=1
    receipt%analysis_fingerprint=11_int64;receipt%basis_fingerprint=12_int64
    receipt%operator_fingerprint=state%operator_value_fingerprint;receipt%state_fingerprint=14_int64
    receipt%action_fingerprint=15_int64
    receipt%action_operator_fingerprint=state%operator_value_fingerprint
    receipt%face_topology_fingerprint=96_int64
  end subroutine acceptance_oracle
  subroutine fill_state(state)
    type(s_dg_hybrid_trial_state),intent(out)::state
    allocate(state%density(nlocal),state%potential(nlocal),state%occupations(2),state%eigenvalues(2),&
      state%mixing_history(2),state%projector(2,2),state%trace(1,1))
    state%density=9d0;state%potential=0d0;state%occupations=[1d0,1d0];state%eigenvalues=[-1d0,-1d0]
    state%mixing_history=0d0;state%projector=(0d0,0d0);state%projector(1,1)=1d0;state%projector(2,2)=1d0
    state%trace=(0d0,0d0)
    state%density_epoch=0;state%operator_epoch=0;state%projector_epoch=0;state%trace_epoch=0;state%derived_epoch=0
    state%trace_cache_valid=.true.
  end subroutine fill_state
  function fixed_density(lambda) result(value)
    real(real64),intent(in)::lambda;real(real64)::value(nlocal)
    do i=1,nlocal;value(i)=0.40d0+0.03d0*real(ids(i),real64)+0.20d0*lambda;enddo
  end function fixed_density
  subroutine volume_build(lambda,input_density,state,callback_ok)
    real(real64),intent(in)::lambda,input_density(:);logical,intent(out)::callback_ok
    type(s_dg_hybrid_trial_state),intent(inout)::state
    callback_ok=phase==0.or.phase==5;phase=1
    last_built_density=input_density
    last_built_trace=state%trace(1,1);last_built_density_epoch=state%density_epoch
    last_built_trace_epoch=state%trace_epoch
    if(.not.(lambda>0d0.and.stale_operator_positive))then
      state%operator_epoch=state%operator_epoch+1
      state%operator_structure_fingerprint=991_int64
      state%operator_value_fingerprint=1001_int64
      if(rank_divergent_operator)state%operator_value_fingerprint=state%operator_value_fingerprint+id_rank
    endif
    if(first_volume)then
      callback_ok=callback_ok.and.all(input_density==dc_density);first_volume=.false.
    endif
    if(lambda>0d0)callback_ok=callback_ok.and.lambda_zero_passed
  end subroutine volume_build
  subroutine full_solve(lambda,iteration,state,callback_ok)
    real(real64),intent(in)::lambda;integer,intent(in)::iteration
    type(s_dg_hybrid_trial_state),intent(inout)::state;logical,intent(out)::callback_ok
    real(real64)::angle
    callback_ok=phase==1.and.iteration>0
    if(solve_count==0)callback_ok=callback_ok.and.all(state%density==dc_density)
    if(lambda>0d0.and.fatal_positive)callback_ok=.false.
    phase=2;solve_count=solve_count+1
    angle=0.37d0*solve_count
    state%eigenvalues=[-1d0-0.1d0*lambda,-1d0-0.1d0*lambda]
    ! A dense unitary gauge rotation changes every solve; the projector callback must remove it.
    state%projector(1,1)=cmplx(cos(angle),0d0,real64);state%projector(1,2)=cmplx(-sin(angle),0d0,real64)
    state%projector(2,1)=cmplx(sin(angle),0d0,real64);state%projector(2,2)=cmplx(cos(angle),0d0,real64)
    state%operator_epoch=state%operator_epoch+1
  end subroutine full_solve
  subroutine projector_refresh(lambda,state,overlap,callback_ok)
    real(real64),intent(in)::lambda
    type(s_dg_hybrid_trial_state),intent(inout)::state
    real(real64),intent(out)::overlap;logical,intent(out)::callback_ok
    callback_ok=phase==2.and.lambda>=0d0;phase=3
    state%projector=(0d0,0d0);state%projector(1,1)=(1d0,0d0);state%projector(2,2)=(1d0,0d0)
    state%projector_epoch=state%operator_epoch;overlap=1d0
  end subroutine projector_refresh
  subroutine density_trace_refresh(lambda,input_density,state,output_density,callback_ok)
    real(real64),intent(in)::lambda,input_density(:)
    type(s_dg_hybrid_trial_state),intent(inout)::state
    real(real64),intent(out)::output_density(:);logical,intent(out)::callback_ok
    callback_ok=phase==3;phase=4
    output_density=fixed_density(lambda)+0.25d0*(input_density-fixed_density(lambda))
    state%density=output_density;state%trace(1,1)=cmplx(sum(output_density)/real(nglobal,real64),0d0,real64)
    call MPI_Allreduce(MPI_IN_PLACE,state%trace,1,MPI_DOUBLE_COMPLEX,MPI_SUM,icomm,ierr)
    state%density_epoch=state%density_epoch+1;state%trace_epoch=state%density_epoch
    state%derived_epoch=state%density_epoch;state%trace_cache_valid=.true.;callback_ok=callback_ok.and.ierr==MPI_SUCCESS
  end subroutine density_trace_refresh
  subroutine residual_evaluation(lambda,iteration,input_density,input_trace,state,residuals,projector_overlap,&
      electron_ok,occupation_ok,hermitian_ok,symmetry_ok,real_space_ok,gap_shrinking,callback_ok)
    real(real64),intent(in)::lambda,input_density(:);integer,intent(in)::iteration
    complex(real64),intent(in)::input_trace(:,:)
    type(s_dg_hybrid_trial_state),intent(in)::state;type(s_dg_hybrid_residuals),intent(out)::residuals
    real(real64),intent(out)::projector_overlap
    logical,intent(out)::electron_ok,occupation_ok,hermitian_ok,symmetry_ok,real_space_ok,gap_shrinking,callback_ok
    real(real64)::local_norm,global_norm
    callback_ok=phase==4.and.state%trace_cache_valid;phase=5
    local_norm=sum((state%density-input_density)**2)
    call MPI_Allreduce(local_norm,global_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierr)
    residuals%r_rho=sqrt(global_norm)/max(1d0,sqrt(sum_global_square(input_density)))
    residuals%r_t=abs(state%trace(1,1)-input_trace(1,1))/max(1d0,abs(input_trace(1,1)))
    residuals%r_h=1d-12;residuals%r_s=1d-12
    if(lambda==controls%initial_step.and..not.forced_growth_complete.and.iteration<=3)then
      select case(iteration)
      case(1);residuals%r_h=1d-6
      case(2);residuals%r_h=5d-6
      case(3);residuals%r_h=3d-5
      end select
      if(iteration==3)forced_growth_complete=.true.
    endif
    if(lambda==1d0.and.lambda_one_converged.and.poison_final)residuals%r_h=1d-6
    projector_overlap=1d0;electron_ok=.true.;occupation_ok=.true.;hermitian_ok=.true.;symmetry_ok=.true.
    real_space_ok=.true.;gap_shrinking=.false.;callback_ok=callback_ok.and.ierr==MPI_SUCCESS
    if(lambda==0d0.and.residuals%r_rho<=controls%intermediate_tolerance(2).and.&
        residuals%r_t<=controls%intermediate_tolerance(3))then
      if(.not.lambda_zero_gate_delayed)then
        symmetry_ok=.false.;lambda_zero_gate_delayed=.true.
      else;lambda_zero_passed=.true.
      endif
    endif
    if(lambda==1d0.and.residuals%r_rho<=controls%final_tolerance(2).and.&
        residuals%r_t<=controls%final_tolerance(3).and.residuals%r_h<=controls%final_tolerance(1))&
      lambda_one_converged=.true.
  end subroutine residual_evaluation
  subroutine density_mix(iteration,input_density,output_density,damping,mixed_density,callback_ok)
    integer,intent(in)::iteration;real(real64),intent(in)::input_density(:),output_density(:),damping
    real(real64),intent(out)::mixed_density(:);logical,intent(out)::callback_ok
    callback_ok=phase==5.and.iteration>0;mixed_density=input_density+damping*(output_density-input_density);phase=0
  end subroutine density_mix
  real(real64) function sum_global_square(values)
    real(real64),intent(in)::values(:);real(real64)::local_value
    local_value=sum(values**2);call MPI_Allreduce(local_value,sum_global_square,1,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierr)
  end function sum_global_square
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label;integer::local_bad,global_bad
    local_bad=merge(0,1,condition);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;if(id_rank==0)write(0,'(a)')trim(label);error stop 1;endif
  end subroutine require
end program test_dg_hybrid_continuation_scf_mpi
