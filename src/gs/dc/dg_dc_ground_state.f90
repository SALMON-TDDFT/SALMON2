!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
module dg_dc_ground_state
#ifdef USE_MPI
  use mpi
#endif
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dg_nodal_state, only: s_dg_nodal_common_state
  implicit none
  private

  integer, parameter, public :: DG_DC_PRECONVERGENCE=1
  integer, parameter, public :: DG_DC_CONTINUATION=2
  integer, parameter, public :: DG_DC_FULL_SCF=3
  integer, parameter, public :: DG_DC_UNMIXED_FIXED_POINT=4
  integer, parameter, public :: DG_DC_ACCEPTED=5
  integer, parameter, public :: DG_DC_FAILED=6

  type, public :: s_dg_dc_gs_controls
    real(8) :: intermediate_orbital_tolerance=1d-5
    real(8) :: intermediate_density_tolerance=1d-5
    real(8) :: final_orbital_tolerance=1d-7
    real(8) :: final_density_tolerance=1d-7
    real(8) :: subspace_tolerance=1d-7
    real(8) :: initial_lambda_step=0.125d0
    real(8) :: minimum_lambda_step=0.015625d0
    real(8) :: maximum_lambda_step=0.5d0
    real(8) :: allowed_residual_growth=4d0
    real(8) :: density_mix_rate=0.5d0
    real(8) :: hermiticity_tolerance=1d-10
    real(8) :: orthogonality_tolerance=1d-10
    real(8) :: face_balance_tolerance=1d-10
    real(8) :: electron_count_tolerance=1d-8
    real(8) :: minimum_projector_overlap=0.9d0
    integer :: maximum_scf_iterations=100
    integer :: maximum_eigensolver_iterations=500
    integer :: maximum_rollbacks=8
  end type s_dg_dc_gs_controls

  type, public :: s_dg_dc_gs_diagnostics
    real(8) :: orbital_residual=huge(1d0)
    real(8) :: density_residual=huge(1d0)
    real(8) :: subspace_residual=huge(1d0)
    real(8) :: projector_overlap=0d0
    real(8) :: hermiticity_defect=huge(1d0)
    real(8) :: orthogonality_defect=huge(1d0)
    real(8) :: face_balance_defect=huge(1d0)
    real(8) :: electron_number=0d0
    real(8) :: expected_electron_number=0d0
    real(8) :: interface_scale=-1d0
    integer :: eigensolver_iterations=0
    integer(8) :: hamiltonian_potential_epoch=0_8
    integer(8) :: updated_potential_epoch=0_8
    logical :: eigensolver_converged=.false.
    logical :: finite=.false.
  end type s_dg_dc_gs_diagnostics

  type, public :: s_dg_dc_gs_result
    integer :: phase=DG_DC_PRECONVERGENCE
    integer :: naccepted_steps=0,nrollbacks=0,mixing_reset_count=0,total_scf_iterations=0
    real(8) :: lambda=0d0,maximum_interface_scale=0d0,minimum_projector_overlap=1d0
    logical :: accepted=.false.,failed=.false.,unmixed_verified=.false.
    real(8), allocatable :: lambda_history(:),lambda_steps(:)
    type(s_dg_dc_gs_diagnostics) :: final_diagnostics
  end type s_dg_dc_gs_result

  abstract interface
    subroutine dg_dc_self_consistent_step(state,density,lambda,density_mix,reset_mixing,unmixed, &
                                          diagnostics,communicator,ok,message)
      import :: s_dg_nodal_common_state,s_dg_dc_gs_diagnostics
      type(s_dg_nodal_common_state), intent(inout) :: state
      real(8), intent(inout) :: density(:,:,:,:)
      real(8), intent(in) :: lambda,density_mix
      logical, intent(in) :: reset_mixing,unmixed
      type(s_dg_dc_gs_diagnostics), intent(out) :: diagnostics
      integer, intent(in) :: communicator
      logical, intent(out) :: ok
      character(*), intent(out) :: message
    end subroutine dg_dc_self_consistent_step
  end interface

  public :: default_dg_dc_gs_controls,run_dg_dc_ground_state
  public :: dg_dc_self_consistent_step

contains

  pure function default_dg_dc_gs_controls() result(controls)
    type(s_dg_dc_gs_controls) :: controls
  end function default_dg_dc_gs_controls

  subroutine run_dg_dc_ground_state(state,density,controls,scf_step,communicator,result,ok,message)
    type(s_dg_nodal_common_state), intent(inout) :: state
    real(8), intent(inout) :: density(:,:,:,:)
    type(s_dg_dc_gs_controls), intent(in) :: controls
    procedure(dg_dc_self_consistent_step) :: scf_step
    integer, intent(in) :: communicator
    type(s_dg_dc_gs_result), intent(out) :: result
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    type(s_dg_nodal_common_state), allocatable :: accepted_state,candidate_state
    real(8), allocatable :: accepted_density(:,:,:,:)
    type(s_dg_dc_gs_diagnostics) :: diagnostics
    real(8) :: lambda,trial_lambda,lambda_step,iteration_residual,previous_iteration_residual
    integer(8) :: geometry_fingerprint,operator_fingerprint
    integer :: max_history,scf_iteration,allocation_status
    logical :: stage_ok,converged,at_final,callback_ok,validation_ok,growth_ok,local_converged

    call validate_controls(controls,communicator,ok,message)
    if(.not.ok) then
      result%phase=DG_DC_FAILED; result%failed=.true.; return
    end if
    max_history=controls%maximum_rollbacks+ceiling(1d0/controls%minimum_lambda_step)+4
    allocate(result%lambda_history(max_history),result%lambda_steps(max_history),stat=allocation_status)
    call collective_and(allocation_status==0,communicator,stage_ok)
    if(.not.stage_ok) then
      result%phase=DG_DC_FAILED;result%failed=.true.;ok=.false.
      message='DG DC GS: continuation history allocation failed collectively';return
    end if
    result%lambda_history=0d0; result%lambda_steps=0d0
    geometry_fingerprint=state%geometry_fingerprint
    operator_fingerprint=state%operator_fingerprint
    allocate(accepted_state,source=state,stat=allocation_status)
    if(allocation_status==0) allocate(accepted_density,source=density,stat=allocation_status)
    call collective_and(allocation_status==0,communicator,stage_ok)
    if(.not.stage_ok) then
      result%phase=DG_DC_FAILED;result%failed=.true.;ok=.false.
      message='DG DC GS: continuation snapshot allocation failed collectively';return
    end if
    lambda=0d0
    lambda_step=controls%initial_lambda_step

    do
      if(lambda==0d0 .and. result%naccepted_steps==0) then
        trial_lambda=0d0
        result%phase=DG_DC_PRECONVERGENCE
      else
        trial_lambda=min(1d0,lambda+lambda_step)
        result%phase=merge(DG_DC_FULL_SCF,DG_DC_CONTINUATION,trial_lambda==1d0)
      end if
      at_final=trial_lambda==1d0
      converged=.false.
      previous_iteration_residual=huge(1d0)
      do scf_iteration=1,controls%maximum_scf_iterations
        call scf_step(state,density,trial_lambda,controls%density_mix_rate,scf_iteration==1,.false., &
          diagnostics,communicator,callback_ok,message)
        call collective_and(callback_ok,communicator,stage_ok)
        result%total_scf_iterations=result%total_scf_iterations+1
        if(scf_iteration==1) result%mixing_reset_count=result%mixing_reset_count+1
        if(.not.stage_ok) exit
        call occupied_projector_overlap(accepted_state,state,communicator,diagnostics%projector_overlap,stage_ok)
        if(.not.stage_ok) then
          message='DG DC GS: occupied-projector tracking failed collectively'
          exit
        end if
        call validate_step(state,diagnostics,trial_lambda,geometry_fingerprint,operator_fingerprint, &
          controls,at_final,communicator,validation_ok,message)
        stage_ok=validation_ok
        if(.not.stage_ok) exit
        result%maximum_interface_scale=max(result%maximum_interface_scale,diagnostics%interface_scale)
        result%minimum_projector_overlap=min(result%minimum_projector_overlap,diagnostics%projector_overlap)
        iteration_residual=max(diagnostics%orbital_residual,diagnostics%density_residual)
        growth_ok=.not.(scf_iteration>1 .and. iteration_residual> &
          controls%allowed_residual_growth*max(previous_iteration_residual,epsilon(1d0)))
        call collective_and(growth_ok,communicator,stage_ok)
        if(.not.stage_ok) then
          message='DG DC GS: residual growth exceeded continuation bound collectively'; exit
        end if
        previous_iteration_residual=iteration_residual
        local_converged=step_converged(diagnostics,controls,at_final)
        call collective_and(local_converged,communicator,converged)
        if(converged) then
          converged=.true.; exit
        end if
      end do

      if(.not.stage_ok .or. .not.converged) then
        call restore_nodal_state(state,accepted_state)
        density=accepted_density
        result%nrollbacks=result%nrollbacks+1
        lambda_step=0.5d0*lambda_step
        if(result%nrollbacks>controls%maximum_rollbacks .or. &
           lambda_step<controls%minimum_lambda_step*(1d0-10d0*epsilon(1d0))) then
          result%phase=DG_DC_FAILED; result%failed=.true.; ok=.false.
          if(len_trim(message)==0) message='DG DC GS: minimum continuation step or rollback limit reached'
          return
        end if
        cycle
      end if

      lambda=trial_lambda
      allocate(candidate_state,source=state,stat=allocation_status)
      call collective_and(allocation_status==0,communicator,stage_ok)
      if(.not.stage_ok) then
        call restore_nodal_state(state,accepted_state);density=accepted_density
        result%phase=DG_DC_FAILED;result%failed=.true.;ok=.false.
        message='DG DC GS: accepted-state snapshot allocation failed collectively';return
      end if
      call move_alloc(candidate_state,accepted_state)
      accepted_density=density
      result%naccepted_steps=result%naccepted_steps+1
      result%lambda_history(result%naccepted_steps)=lambda
      result%lambda_steps(result%naccepted_steps)=merge(0d0,lambda_step,result%naccepted_steps==1)
      result%lambda=lambda
      if(lambda==1d0) exit
      lambda_step=min(controls%maximum_lambda_step,1.5d0*lambda_step)
    end do

    result%phase=DG_DC_UNMIXED_FIXED_POINT
    call scf_step(state,density,1d0,1d0,.true.,.true.,diagnostics,communicator,callback_ok,message)
    call collective_and(callback_ok,communicator,stage_ok)
    if(stage_ok) call occupied_projector_overlap(accepted_state,state,communicator,diagnostics%projector_overlap,stage_ok)
    if(.not.stage_ok) then
      call restore_nodal_state(state,accepted_state); density=accepted_density
      result%phase=DG_DC_FAILED; result%failed=.true.; ok=.false.
      if(len_trim(message)==0) message='DG DC GS: unmixed callback/projector failure'
      return
    end if
    call validate_step(state,diagnostics,1d0,geometry_fingerprint,operator_fingerprint,controls,.true., &
      communicator,validation_ok,message)
    local_converged=validation_ok .and. step_converged(diagnostics,controls,.true.)
    call collective_and(local_converged,communicator,stage_ok)
    if(.not.stage_ok) then
      call restore_nodal_state(state,accepted_state); density=accepted_density
      result%phase=DG_DC_FAILED; result%failed=.true.; ok=.false.
      if(len_trim(message)==0) message='DG DC GS: unmixed fixed-point verification failed'
      return
    end if
    result%unmixed_verified=.true.
    result%accepted=.true.
    result%phase=DG_DC_ACCEPTED
    result%final_diagnostics=diagnostics
    state%ground_state_ready=.true.
    state%dg_ground_state_residual=diagnostics%orbital_residual
    ok=.true.; message=''
  end subroutine run_dg_dc_ground_state

  subroutine restore_nodal_state(state,snapshot)
    type(s_dg_nodal_common_state), intent(inout) :: state
    type(s_dg_nodal_common_state), intent(in) :: snapshot
    state%enabled=snapshot%enabled
    state%initialized=snapshot%initialized
    state%ground_state_ready=snapshot%ground_state_ready
    state%geometry_fingerprint=snapshot%geometry_fingerprint
    state%operator_fingerprint=snapshot%operator_fingerprint
    state%dg_ground_state_residual=snapshot%dg_ground_state_residual
    state%psi_core=snapshot%psi_core
    state%hpsi_core=snapshot%hpsi_core
    state%occupations=snapshot%occupations
  end subroutine restore_nodal_state

  logical function step_converged(diagnostics,controls,final_stage)
    type(s_dg_dc_gs_diagnostics), intent(in) :: diagnostics
    type(s_dg_dc_gs_controls), intent(in) :: controls
    logical, intent(in) :: final_stage
    real(8) :: orbital_tolerance,density_tolerance
    orbital_tolerance=merge(controls%final_orbital_tolerance,controls%intermediate_orbital_tolerance,final_stage)
    density_tolerance=merge(controls%final_density_tolerance,controls%intermediate_density_tolerance,final_stage)
    step_converged=diagnostics%eigensolver_converged .and. &
      diagnostics%orbital_residual<=orbital_tolerance .and. &
      diagnostics%density_residual<=density_tolerance .and. &
      diagnostics%subspace_residual<=controls%subspace_tolerance
  end function step_converged

  subroutine validate_step(state,diagnostics,lambda,geometry_fingerprint,operator_fingerprint,controls, &
                           final_stage,communicator,ok,message)
    type(s_dg_nodal_common_state), intent(in) :: state
    type(s_dg_dc_gs_diagnostics), intent(in) :: diagnostics
    real(8), intent(in) :: lambda
    integer(8), intent(in) :: geometry_fingerprint,operator_fingerprint
    type(s_dg_dc_gs_controls), intent(in) :: controls
    logical, intent(in) :: final_stage
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    logical :: local_ok
    real(8) :: values(10)
    values=[diagnostics%orbital_residual,diagnostics%density_residual,diagnostics%subspace_residual, &
      diagnostics%projector_overlap,diagnostics%hermiticity_defect,diagnostics%orthogonality_defect, &
      diagnostics%face_balance_defect,diagnostics%electron_number,diagnostics%expected_electron_number, &
      diagnostics%interface_scale]
    local_ok=diagnostics%finite .and. all(ieee_is_finite(values)) .and. &
      diagnostics%eigensolver_iterations>=0 .and. &
      diagnostics%eigensolver_iterations<=controls%maximum_eigensolver_iterations .and. &
      abs(diagnostics%interface_scale-lambda)<=10d0*epsilon(1d0) .and. &
      diagnostics%hermiticity_defect<=controls%hermiticity_tolerance .and. &
      diagnostics%orthogonality_defect<=controls%orthogonality_tolerance .and. &
      diagnostics%face_balance_defect<=controls%face_balance_tolerance .and. &
      diagnostics%projector_overlap>=controls%minimum_projector_overlap .and. diagnostics%projector_overlap<=1d0 .and. &
      abs(diagnostics%electron_number-diagnostics%expected_electron_number)<=controls%electron_count_tolerance .and. &
      ((diagnostics%hamiltonian_potential_epoch==0_8 .and. diagnostics%updated_potential_epoch==0_8) .or. &
       diagnostics%updated_potential_epoch==diagnostics%hamiltonian_potential_epoch+1_8) .and. &
      state%geometry_fingerprint==geometry_fingerprint .and. state%operator_fingerprint==operator_fingerprint
    call collective_and(local_ok,communicator,ok)
    if(ok) then
      message=''
    else
      message='DG DC GS: invalid, stale, or nonfinite collective diagnostics'
    end if
  end subroutine validate_step

  subroutine validate_controls(controls,communicator,ok,message)
    type(s_dg_dc_gs_controls), intent(in) :: controls
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    real(8) :: values(14)
    logical :: local_ok
    values=[controls%intermediate_orbital_tolerance,controls%intermediate_density_tolerance, &
      controls%final_orbital_tolerance,controls%final_density_tolerance,controls%subspace_tolerance, &
      controls%initial_lambda_step,controls%minimum_lambda_step,controls%maximum_lambda_step, &
      controls%allowed_residual_growth,controls%density_mix_rate,controls%hermiticity_tolerance, &
      controls%orthogonality_tolerance,controls%face_balance_tolerance,controls%minimum_projector_overlap]
    local_ok=all(ieee_is_finite(values)) .and. all(values>0d0) .and. &
      controls%density_mix_rate<=1d0 .and. controls%initial_lambda_step>=controls%minimum_lambda_step .and. &
      controls%maximum_lambda_step>=controls%initial_lambda_step .and. controls%maximum_lambda_step<=1d0 .and. &
      ieee_is_finite(controls%electron_count_tolerance) .and. controls%electron_count_tolerance>0d0 .and. &
      controls%maximum_scf_iterations>0 .and. controls%maximum_eigensolver_iterations>0 .and. &
      controls%maximum_rollbacks>=0 .and. controls%maximum_rollbacks<=100000 .and. &
      controls%minimum_lambda_step>=1d-6
    call collective_and(local_ok,communicator,ok)
    if(ok) then
      message=''
    else
      message='DG DC GS: invalid continuation controls'
    end if
  end subroutine validate_controls

  subroutine occupied_projector_overlap(reference,current,communicator,overlap,ok)
    type(s_dg_nodal_common_state), intent(in) :: reference,current
    integer, intent(in) :: communicator
    real(8), intent(out) :: overlap
    logical, intent(out) :: ok
    complex(8), allocatable :: local_overlap(:,:),global_overlap(:,:)
    integer :: io,jo,is,noccupied,current_noccupied,ierr
    integer, allocatable :: reference_index(:),current_index(:)
    real(8) :: local_value,global_value

    ok=reference%nstate==current%nstate .and. reference%nspin==current%nspin .and. &
      all(reference%core_size==current%core_size)
    call collective_and(ok,communicator,ok)
    if(.not.ok) return
    overlap=1d0
    do is=1,reference%nspin
      noccupied=count(reference%occupations(:,is)>0d0)
      if(noccupied==0) cycle
      current_noccupied=count(current%occupations(:,is)>0d0)
      if(current_noccupied/=noccupied) then
        ok=.false.
        call collective_and(ok,communicator,ok)
        return
      end if
      allocate(reference_index(noccupied),current_index(noccupied),stat=ierr)
      call collective_and(ierr==0,communicator,ok)
      if(.not.ok) return
      reference_index=pack([(io,io=1,reference%nstate)],reference%occupations(:,is)>0d0)
      current_index=pack([(io,io=1,current%nstate)],current%occupations(:,is)>0d0)
      allocate(local_overlap(noccupied,noccupied),global_overlap(noccupied,noccupied),stat=ierr)
      call collective_and(ierr==0,communicator,ok)
      if(.not.ok) return
      do jo=1,noccupied
      do io=1,noccupied
        local_overlap(io,jo)=sum(conjg(reference%psi_core(:,:,:,reference_index(io),is))* &
          current%psi_core(:,:,:,current_index(jo),is))
      end do
      end do
#ifdef USE_MPI
      call MPI_Allreduce(local_overlap,global_overlap,size(local_overlap),MPI_DOUBLE_COMPLEX,MPI_SUM,communicator,ierr)
      call collective_and(ierr==MPI_SUCCESS,communicator,ok)
      if(.not.ok) return
#else
      global_overlap=local_overlap
#endif
      local_value=sum(abs(global_overlap)**2)/dble(noccupied)
      overlap=min(overlap,max(0d0,min(1d0,local_value)))
      deallocate(local_overlap,global_overlap,reference_index,current_index)
    end do
    local_value=overlap
#ifdef USE_MPI
    call MPI_Allreduce(local_value,global_value,1,MPI_DOUBLE_PRECISION,MPI_MIN,communicator,ierr)
    call collective_and(ierr==MPI_SUCCESS,communicator,ok)
    if(ok) overlap=global_value
#else
    overlap=local_value
#endif
  end subroutine occupied_projector_overlap

  subroutine collective_and(local,communicator,global)
    logical, intent(in) :: local
    integer, intent(in) :: communicator
    logical, intent(out) :: global
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local,global,1,MPI_LOGICAL,MPI_LAND,communicator,ierr)
    if(ierr/=MPI_SUCCESS) global=.false.
#else
    global=local .and. communicator>=0
#endif
  end subroutine collective_and
end module dg_dc_ground_state
