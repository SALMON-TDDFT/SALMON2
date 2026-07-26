!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
module dg_dc_handoff
#ifdef USE_MPI
  use mpi
#endif
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dg_nodal_state, only: s_dg_nodal_common_state, initialize_dg_nodal_common_state, &
                            release_dg_nodal_common_state, validate_dg_nodal_common_state_mpi, &
                            validate_unique_core_ownership
  implicit none
  private

  type, public :: s_dg_dc_handoff_state
    logical :: enabled=.false., accepted=.false., materialized=.false.
    logical :: lcfo_started=.false.,wannier_started=.false.,window_truncated=.false.
    integer :: minimum_iteration=0, accepted_iteration=0,direct_dg_solve_count=0
    integer :: candidate_orbitals_per_atom=0, candidate_count=0, metric_rank=0
    real(8) :: tolerance=0.0d0, metric_rank_tolerance=0.0d0
    real(8) :: interface_scale=0.0d0
    real(8) :: accepted_interface_scale=0d0,interface_step=0d0
    real(8) :: minimum_interface_step=0d0,maximum_interface_step=0d0
    real(8) :: allowed_residual_growth=0d0,accepted_stage_residual=huge(1d0)
    real(8) :: trial_stage_first_residual=huge(1d0)
    real(8) :: direct_orbital_residual=huge(1d0),direct_orthogonality_defect=huge(1d0)
    real(8) :: direct_hermiticity_defect=huge(1d0),direct_face_balance_defect=huge(1d0)
    real(8) :: residual=huge(1.0d0), charge_error=huge(1.0d0)
    real(8) :: minimum_metric_pivot=0.0d0, span_escape_norm=0.0d0
    integer(8) :: geometry_fingerprint=0_8, operator_fingerprint=0_8
    integer(8) :: projector_fingerprint=0_8
    integer :: projector_generation=0,projector_retained_rank(6)=0,projector_required_rank(6)=0
    real(8) :: projector_projection_residual(6)=huge(1d0)
    real(8) :: projector_escape_norm(6)=huge(1d0)
    real(8) :: projector_residual_limit(6)=0d0,projector_escape_limit(6)=0d0
    real(8) :: face_action_norm(6)=huge(1d0),face_pair_balance(3)=huge(1d0)
    complex(8), allocatable :: dc_reference_core(:,:,:,:,:)
    real(8), allocatable :: density_snapshot(:)
    real(8), allocatable :: potential_snapshot(:,:)
    real(8), allocatable :: density_spin_snapshot(:,:),hartree_snapshot(:),xc_snapshot(:,:)
    integer(8), allocatable :: snapshot_owner(:)
    real(8), allocatable :: fragment_orbital_snapshot(:,:,:,:,:,:,:)
    real(8), allocatable :: occupation_snapshot(:,:,:),eigenvalue_snapshot(:,:,:)
    real(8) :: chemical_potential_snapshot=0d0
    real(8), allocatable :: continuation_scale_history(:),continuation_step_history(:)
    logical :: density_potential_preserved=.false.
    logical :: mixing_history_discarded=.false.
    logical :: continuation_initialized=.false.,direct_ground_state_complete=.false.,continuation_failed=.false.
    integer :: continuation_rollbacks=0,maximum_continuation_rollbacks=0
    integer :: stage_iteration_count=0
  end type s_dg_dc_handoff_state

  type(s_dg_dc_handoff_state), public, save :: dg_dc_handoff_runtime
  type(s_dg_nodal_common_state), public, save :: dg_dc_nodal_runtime

  public :: initialize_dg_dc_handoff, evaluate_dg_dc_handoff
  public :: materialize_dg_dc_candidates, update_dg_dc_nodal_state
  public :: preserve_dg_dc_density_potential, mark_dg_dc_mixing_discarded
  public :: initialize_dg_dc_direct_continuation,evaluate_dg_dc_direct_continuation
  public :: preserve_dg_dc_fragment_orbitals,restore_dg_dc_fragment_orbitals

contains

  subroutine initialize_dg_dc_handoff(state,enabled,minimum_iteration,tolerance,candidates_per_atom, &
                                       metric_tolerance,ok,message)
    type(s_dg_dc_handoff_state), intent(inout) :: state
    logical, intent(in) :: enabled
    integer, intent(in) :: minimum_iteration,candidates_per_atom
    real(8), intent(in) :: tolerance,metric_tolerance
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    type(s_dg_dc_handoff_state) :: candidate

    ok=minimum_iteration > 0 .and. candidates_per_atom > 0 .and. ieee_is_finite(tolerance) .and. &
       tolerance > 0.0d0 .and. ieee_is_finite(metric_tolerance) .and. metric_tolerance > 0.0d0
    if (.not.ok) then
      message='DG DC handoff: invalid controls'
      return
    end if
    candidate%enabled=enabled
    candidate%minimum_iteration=minimum_iteration
    candidate%tolerance=tolerance
    candidate%candidate_orbitals_per_atom=candidates_per_atom
    candidate%metric_rank_tolerance=metric_tolerance
    call move_handoff(candidate,state)
    message=''
  end subroutine initialize_dg_dc_handoff

  subroutine evaluate_dg_dc_handoff(state,iteration,residual,charge,expected_charge,fragment_orbitals_finite, &
                                     lcfo_started,wannier_started,window_truncated,communicator,accept,message)
    type(s_dg_dc_handoff_state), intent(inout) :: state
    integer, intent(in) :: iteration,communicator
    real(8), intent(in) :: residual,charge,expected_charge
    logical, intent(in) :: fragment_orbitals_finite,lcfo_started,wannier_started,window_truncated
    logical, intent(out) :: accept
    character(*), intent(out) :: message
    logical :: local_accept,metadata_ok,real_metadata_ok,integer_metadata_ok
    logical :: fragment_complete,route_clear
    real(8) :: charge_tolerance,minimum_values(3),maximum_values(3),local_values(3)
    integer :: minimum_iteration_seen,maximum_iteration_seen

    accept=.false.; message=''
    if (.not.state%enabled .or. state%accepted) return
    local_values=[residual,abs(charge-expected_charge),expected_charge]
    call collective_real_minmax(local_values,minimum_values,maximum_values,communicator,real_metadata_ok)
    call collective_integer_minmax(iteration,minimum_iteration_seen,maximum_iteration_seen,communicator,integer_metadata_ok)
    metadata_ok=real_metadata_ok .and. integer_metadata_ok
    if(.not.metadata_ok .or. minimum_iteration_seen/=maximum_iteration_seen) then
      accept=.false.
      message='DG DC handoff: invalid or rank-disagreeing iteration diagnostics'
      return
    end if
    charge_tolerance=max(state%tolerance,100.0d0*epsilon(1.0d0)* &
      max(1.0d0,maxval(abs([minimum_values(3),maximum_values(3)]))))
    call collective_and(fragment_orbitals_finite,communicator,fragment_complete)
    call collective_and(.not.lcfo_started .and. .not.wannier_started .and. .not.window_truncated, &
      communicator,route_clear)
    state%residual=maximum_values(1)
    state%charge_error=maximum_values(2)
    local_accept=minimum_iteration_seen >= state%minimum_iteration .and. &
      all(ieee_is_finite(minimum_values)) .and. all(ieee_is_finite(maximum_values)) .and. &
      maximum_values(1) <= state%tolerance .and. maximum_values(2) <= charge_tolerance .and. &
      fragment_complete .and. route_clear
    call collective_and(local_accept,communicator,accept)
    if (.not.accept) then
      if(minimum_iteration_seen<state%minimum_iteration)then
        message='DG DC handoff: minimum iteration not reached'
      else if(.not.all(ieee_is_finite(minimum_values)) .or. .not.all(ieee_is_finite(maximum_values)))then
        message='DG DC handoff: nonfinite residual or charge diagnostic'
      else if(maximum_values(1)>state%tolerance)then
        message='DG DC handoff: density residual exceeds tolerance'
      else if(maximum_values(2)>charge_tolerance)then
        message='DG DC handoff: electron-count error exceeds tolerance'
      else if(.not.fragment_complete)then
        message='DG DC handoff: fragment orbital payload is incomplete or nonfinite'
      else
        message='DG DC handoff: LCFO, Wannier, or window-truncation route already started'
      end if
      return
    end if
    state%accepted=.true.
    state%accepted_iteration=minimum_iteration_seen
  end subroutine evaluate_dg_dc_handoff

  subroutine materialize_dg_dc_candidates(state,nodal,candidate_box,fragment,core_size,buffer,owner,natom, &
                                           configured_count,geometry_fingerprint,operator_fingerprint, &
                                           communicator,ok,message)
    type(s_dg_dc_handoff_state), intent(inout) :: state
    type(s_dg_nodal_common_state), intent(inout) :: nodal
    complex(8), intent(in) :: candidate_box(:,:,:,:,:)
    integer, intent(in) :: fragment,core_size(3),buffer(3),natom,configured_count,communicator
    integer(8), intent(in) :: owner(:,:,:),geometry_fingerprint,operator_fingerprint
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    type(s_dg_dc_handoff_state) :: candidate_state
    type(s_dg_nodal_common_state) :: candidate_nodal
    integer :: nstate,nspin,rank_count
    real(8) :: minimum_pivot,global_minimum_pivot
    logical :: local_ok

    nstate=size(candidate_box,4); nspin=size(candidate_box,5)
    local_ok=state%accepted .and. fragment>0 .and. natom > 0 .and. configured_count>0 .and. &
      configured_count<=nstate .and. &
      all(shape(candidate_box(:,:,:,1,1))==core_size+2*buffer) .and. &
      all(ieee_is_finite(real(candidate_box,8))) .and. all(ieee_is_finite(aimag(candidate_box)))
    if(local_ok) then
      call candidate_metric_rank(candidate_box(:,:,:,1:configured_count,:),core_size,buffer, &
        state%metric_rank_tolerance,rank_count,minimum_pivot)
      local_ok=rank_count==configured_count
    else
      rank_count=0
      minimum_pivot=0d0
    end if
    call collective_and(local_ok,communicator,ok)
    if (.not.ok) then
      message='DG DC handoff: incomplete or rank-deficient candidate pool'
      return
    end if
    call collective_real_min(minimum_pivot,global_minimum_pivot,communicator,ok)
    if(.not.ok) then
      message='DG DC handoff: metric diagnostic reduction failed collectively'
      return
    end if
    call initialize_dg_nodal_common_state(candidate_nodal,fragment,core_size,core_size+2*buffer,buffer,1, &
      nstate,nspin,owner,geometry_fingerprint,operator_fingerprint,communicator,ok,message)
    if (.not.ok) return
    candidate_nodal%psi_core=candidate_box(buffer(1)+1:buffer(1)+core_size(1), &
      buffer(2)+1:buffer(2)+core_size(2),buffer(3)+1:buffer(3)+core_size(3),:,:)
    call validate_dg_nodal_common_state_mpi(candidate_nodal,communicator,ok,message)
    if (.not.ok) then
      call release_dg_nodal_common_state(candidate_nodal)
      return
    end if
    candidate_state=state
    if (allocated(candidate_state%dc_reference_core)) deallocate(candidate_state%dc_reference_core)
    allocate(candidate_state%dc_reference_core,source=candidate_nodal%psi_core)
    candidate_state%materialized=.true.
    candidate_state%candidate_count=configured_count
    candidate_state%metric_rank=rank_count
    candidate_state%minimum_metric_pivot=global_minimum_pivot
    candidate_state%geometry_fingerprint=geometry_fingerprint
    candidate_state%operator_fingerprint=operator_fingerprint
    call move_handoff(candidate_state,state)
    call move_nodal(candidate_nodal,nodal)
    message=''
  end subroutine materialize_dg_dc_candidates

  subroutine update_dg_dc_nodal_state(state,nodal,updated,communicator,ok,message)
    type(s_dg_dc_handoff_state), intent(inout) :: state
    type(s_dg_nodal_common_state), intent(inout) :: nodal
    complex(8), intent(in) :: updated(:,:,:,:,:)
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    complex(8), allocatable :: basis(:,:),work(:),q(:)
    integer :: ngrid,nstate,io,jo,nq
    real(8) :: norm_q,local_escape,global_escape
    logical :: local_ok

    local_ok=state%materialized .and. allocated(state%dc_reference_core) .and. &
      all(shape(updated)==shape(nodal%psi_core)) .and. all(ieee_is_finite(real(updated,8))) .and. &
      all(ieee_is_finite(aimag(updated)))
    call collective_and(local_ok,communicator,ok)
    if (.not.ok) then
      message='DG DC handoff: invalid nodal update'
      return
    end if
    ngrid=product(nodal%core_size)*nodal%nspin
    nstate=nodal%nstate
    allocate(basis(ngrid,nstate),work(ngrid),q(ngrid))
    do io=1,nstate
      basis(:,io)=reshape(state%dc_reference_core(:,:,:,io,:),[ngrid])
    end do
    nq=0
    do io=1,nstate
      q=basis(:,io)
      do jo=1,nq
        q=q-dot_product(basis(:,jo),q)*basis(:,jo)
      end do
      norm_q=sqrt(max(0.0d0,real(dot_product(q,q),8)))
      if (norm_q > state%metric_rank_tolerance) then
        nq=nq+1
        basis(:,nq)=q/norm_q
      end if
    end do
    local_escape=0.0d0
    do io=1,nstate
      work=reshape(updated(:,:,:,io,:),[ngrid])
      do jo=1,nq
        work=work-dot_product(basis(:,jo),work)*basis(:,jo)
      end do
      local_escape=local_escape+real(dot_product(work,work),8)
    end do
    call collective_sum(local_escape,global_escape,communicator,ok)
    if(.not.ok) then
      message='DG DC handoff: span diagnostic reduction failed collectively'
      return
    end if
    nodal%psi_core=updated
    state%span_escape_norm=sqrt(max(0.0d0,global_escape))
    message=''
  end subroutine update_dg_dc_nodal_state

  subroutine preserve_dg_dc_density_potential(state,density,potential,owner,communicator,ok,message, &
      density_spin,hartree,xc)
    type(s_dg_dc_handoff_state), intent(inout) :: state
    real(8), intent(in) :: density(:),potential(:,:)
    integer(8), intent(in) :: owner(:)
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    real(8),optional,intent(in)::density_spin(:,:),hartree(:),xc(:,:)
    logical :: local_ok,ownership_ok
    integer :: allocation_status
    integer(8) :: local_count,global_count,local_minimum,global_minimum,local_maximum,global_maximum
    real(8), allocatable :: density_candidate(:),potential_candidate(:,:)
    real(8),allocatable::density_spin_candidate(:,:),hartree_candidate(:),xc_candidate(:,:)
    integer(8), allocatable :: owner_candidate(:)
    local_ok=state%accepted .and. size(density)>0 .and. size(potential,1)==size(density) .and. &
      size(potential,2)>0 .and. all(ieee_is_finite(density)) .and. all(ieee_is_finite(potential))
    local_ok=local_ok .and. size(owner)==size(density) .and. all(owner>0_8)
    if(present(density_spin).or.present(hartree).or.present(xc))then
      local_ok=local_ok.and.present(density_spin).and.present(hartree).and.present(xc)
      if(local_ok)local_ok=size(density_spin,1)==size(density).and.size(hartree)==size(density).and. &
        all(shape(xc)==shape(density_spin)).and.all(ieee_is_finite(density_spin)).and. &
        all(ieee_is_finite(hartree)).and.all(ieee_is_finite(xc))
    endif
    call collective_and(local_ok,communicator,ok)
    if(.not.ok) then
      message='DG DC handoff: invalid density/potential snapshot'
      return
    end if
    local_count=int(size(owner),8)
    local_minimum=minval(owner)
    local_maximum=maxval(owner)
    call collective_integer8_sum(local_count,global_count,communicator,ownership_ok)
    if(ownership_ok) call collective_integer8_min(local_minimum,global_minimum,communicator,ownership_ok)
    if(ownership_ok) call collective_integer8_max(local_maximum,global_maximum,communicator,ownership_ok)
    ownership_ok=ownership_ok .and. global_minimum==1_8 .and. global_count==global_maximum
    call collective_and(ownership_ok,communicator,ok)
    if(.not.ok) then
      message='DG DC handoff: distributed snapshot ownership is not complete and unique'
      return
    end if
    call validate_unique_core_ownership(reshape(owner,[size(owner),1,1]),communicator,ok,message)
    if(.not.ok) return
    allocate(density_candidate,source=density,stat=allocation_status)
    local_ok=allocation_status==0
    if(local_ok) allocate(potential_candidate,source=potential,stat=allocation_status)
    local_ok=local_ok .and. allocation_status==0
    if(local_ok) allocate(owner_candidate,source=owner,stat=allocation_status)
    local_ok=local_ok .and. allocation_status==0
    if(local_ok.and.present(density_spin))then
      allocate(density_spin_candidate,source=density_spin,stat=allocation_status)
      local_ok=allocation_status==0
    endif
    if(local_ok.and.present(hartree))then
      allocate(hartree_candidate,source=hartree,stat=allocation_status)
      local_ok=allocation_status==0
    endif
    if(local_ok.and.present(xc))then
      allocate(xc_candidate,source=xc,stat=allocation_status)
      local_ok=allocation_status==0
    endif
    call collective_and(local_ok,communicator,ok)
    if(.not.ok) then
      if(allocated(density_candidate)) deallocate(density_candidate)
      if(allocated(potential_candidate)) deallocate(potential_candidate)
      if(allocated(owner_candidate)) deallocate(owner_candidate)
      if(allocated(density_spin_candidate))deallocate(density_spin_candidate)
      if(allocated(hartree_candidate))deallocate(hartree_candidate)
      if(allocated(xc_candidate))deallocate(xc_candidate)
      message='DG DC handoff: density/potential snapshot allocation failed collectively'
      return
    end if
    if(allocated(state%density_snapshot)) deallocate(state%density_snapshot)
    if(allocated(state%potential_snapshot)) deallocate(state%potential_snapshot)
    if(allocated(state%snapshot_owner)) deallocate(state%snapshot_owner)
    if(allocated(state%density_spin_snapshot))deallocate(state%density_spin_snapshot)
    if(allocated(state%hartree_snapshot))deallocate(state%hartree_snapshot)
    if(allocated(state%xc_snapshot))deallocate(state%xc_snapshot)
    call move_alloc(density_candidate,state%density_snapshot)
    call move_alloc(potential_candidate,state%potential_snapshot)
    call move_alloc(owner_candidate,state%snapshot_owner)
    if(allocated(density_spin_candidate))call move_alloc(density_spin_candidate,state%density_spin_snapshot)
    if(allocated(hartree_candidate))call move_alloc(hartree_candidate,state%hartree_snapshot)
    if(allocated(xc_candidate))call move_alloc(xc_candidate,state%xc_snapshot)
    state%density_potential_preserved=.true.
    message=''
  end subroutine preserve_dg_dc_density_potential

  subroutine mark_dg_dc_mixing_discarded(state)
    type(s_dg_dc_handoff_state), intent(inout) :: state
    state%mixing_history_discarded=.true.
  end subroutine mark_dg_dc_mixing_discarded

  subroutine preserve_dg_dc_fragment_orbitals(state,orbitals,communicator,ok,message, &
      occupations,chemical_potential,eigenvalues)
    type(s_dg_dc_handoff_state),intent(inout)::state
    real(8),intent(in)::orbitals(:,:,:,:,:,:,:)
    integer,intent(in)::communicator
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8),optional,intent(in)::occupations(:,:,:),chemical_potential,eigenvalues(:,:,:)
    real(8),allocatable::candidate(:,:,:,:,:,:,:)
    real(8),allocatable::occupation_candidate(:,:,:),eigenvalue_candidate(:,:,:)
    integer::allocation_status
    logical::local_ok
    local_ok=state%accepted.and.size(orbitals)>0.and.all(ieee_is_finite(orbitals))
    if(local_ok)allocate(candidate,source=orbitals,stat=allocation_status)
    if(.not.local_ok)allocation_status=1
    local_ok=local_ok.and.allocation_status==0
    if(local_ok.and.present(occupations))then
      allocate(occupation_candidate,source=occupations,stat=allocation_status)
      local_ok=allocation_status==0.and.all(ieee_is_finite(occupations)).and.all(occupations>=0d0)
    endif
    if(local_ok.and.present(eigenvalues))then
      allocate(eigenvalue_candidate,source=eigenvalues,stat=allocation_status)
      local_ok=allocation_status==0.and.all(ieee_is_finite(eigenvalues))
    endif
    if(local_ok.and.present(chemical_potential))local_ok=ieee_is_finite(chemical_potential)
    local_ok=local_ok.and.(present(occupations).eqv.present(chemical_potential)).and. &
      (present(occupations).eqv.present(eigenvalues))
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then
      if(allocated(candidate))deallocate(candidate)
      if(allocated(occupation_candidate))deallocate(occupation_candidate)
      if(allocated(eigenvalue_candidate))deallocate(eigenvalue_candidate)
      message='DG DC direct continuation: orbital snapshot failed collectively'
      return
    endif
    if(allocated(state%fragment_orbital_snapshot))deallocate(state%fragment_orbital_snapshot)
    if(allocated(state%occupation_snapshot))deallocate(state%occupation_snapshot)
    if(allocated(state%eigenvalue_snapshot))deallocate(state%eigenvalue_snapshot)
    call move_alloc(candidate,state%fragment_orbital_snapshot)
    if(allocated(occupation_candidate))call move_alloc(occupation_candidate,state%occupation_snapshot)
    if(allocated(eigenvalue_candidate))call move_alloc(eigenvalue_candidate,state%eigenvalue_snapshot)
    if(present(chemical_potential))state%chemical_potential_snapshot=chemical_potential
    message=''
  end subroutine preserve_dg_dc_fragment_orbitals

  subroutine restore_dg_dc_fragment_orbitals(state,orbitals,communicator,ok,message, &
      occupations,chemical_potential,eigenvalues)
    type(s_dg_dc_handoff_state),intent(in)::state
    real(8),intent(out)::orbitals(:,:,:,:,:,:,:)
    integer,intent(in)::communicator
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8),optional,intent(out)::occupations(:,:,:),chemical_potential,eigenvalues(:,:,:)
    logical::local_ok
    local_ok=allocated(state%fragment_orbital_snapshot)
    if(local_ok)local_ok=all(shape(orbitals)==shape(state%fragment_orbital_snapshot))
    if(local_ok.and.present(occupations))then
      local_ok=allocated(state%occupation_snapshot)
      if(local_ok)local_ok=all(shape(occupations)==shape(state%occupation_snapshot))
    endif
    if(local_ok.and.present(eigenvalues))then
      local_ok=allocated(state%eigenvalue_snapshot)
      if(local_ok)local_ok=all(shape(eigenvalues)==shape(state%eigenvalue_snapshot))
    endif
    local_ok=local_ok.and.(present(occupations).eqv.present(chemical_potential)).and. &
      (present(occupations).eqv.present(eigenvalues))
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)then
      orbitals=0d0
      message='DG DC direct continuation: orbital restore layout failed collectively'
      return
    endif
    orbitals=state%fragment_orbital_snapshot
    if(present(occupations))occupations=state%occupation_snapshot
    if(present(eigenvalues))eigenvalues=state%eigenvalue_snapshot
    if(present(chemical_potential))chemical_potential=state%chemical_potential_snapshot
    message=''
  end subroutine restore_dg_dc_fragment_orbitals

  subroutine initialize_dg_dc_direct_continuation(state,initial_step,minimum_step,maximum_step, &
      allowed_growth,maximum_rollbacks,ok,message)
    type(s_dg_dc_handoff_state),intent(inout)::state
    real(8),intent(in)::initial_step,minimum_step,maximum_step,allowed_growth
    integer,intent(in)::maximum_rollbacks
    logical,intent(out)::ok
    character(*),intent(out)::message
    ok=state%accepted.and.ieee_is_finite(initial_step).and.ieee_is_finite(minimum_step).and. &
      ieee_is_finite(maximum_step).and.ieee_is_finite(allowed_growth).and.minimum_step>0d0.and. &
      minimum_step<=initial_step.and.initial_step<=maximum_step.and.maximum_step<=1d0.and.allowed_growth>1d0.and. &
      maximum_rollbacks>=0
    if(.not.ok)then
      message='DG DC direct continuation: invalid initialization controls'
      return
    endif
    state%accepted_interface_scale=0d0
    state%interface_step=initial_step
    state%minimum_interface_step=minimum_step
    state%maximum_interface_step=maximum_step
    state%allowed_residual_growth=allowed_growth
    state%accepted_stage_residual=state%residual
    state%trial_stage_first_residual=huge(1d0)
    state%stage_iteration_count=0
    state%interface_scale=initial_step
    state%continuation_rollbacks=0
    state%maximum_continuation_rollbacks=maximum_rollbacks
    state%continuation_failed=.false.
    if(allocated(state%continuation_scale_history))deallocate(state%continuation_scale_history)
    if(allocated(state%continuation_step_history))deallocate(state%continuation_step_history)
    allocate(state%continuation_scale_history(1),state%continuation_step_history(1))
    state%continuation_scale_history=0d0
    state%continuation_step_history=0d0
    state%continuation_initialized=.true.
    state%direct_ground_state_complete=.false.
    message=''
  end subroutine initialize_dg_dc_direct_continuation

  subroutine evaluate_dg_dc_direct_continuation(state,residual,tolerance,unmixed,communicator, &
      advance,rollback,complete,message)
    type(s_dg_dc_handoff_state),intent(inout)::state
    real(8),intent(in)::residual,tolerance
    logical,intent(in)::unmixed
    integer,intent(in)::communicator
    logical,intent(out)::advance,rollback,complete
    character(*),intent(out)::message
    real(8)::local_residual(1),minimum_residual(1),maximum_residual(1),previous_accepted_scale
    logical::ok,collective_unmixed
    advance=.false.;rollback=.false.;complete=.false.
    local_residual(1)=residual
    call collective_real_minmax(local_residual,minimum_residual,maximum_residual,communicator,ok)
    call collective_and(unmixed,communicator,collective_unmixed)
    if(.not.ok.or..not.state%continuation_initialized.or..not.ieee_is_finite(maximum_residual(1)).or. &
       .not.ieee_is_finite(tolerance).or.tolerance<=0d0)then
      message='DG DC direct continuation: invalid collective stage diagnostics'
      return
    endif
    if(state%interface_scale>state%accepted_interface_scale)then
      state%stage_iteration_count=state%stage_iteration_count+1
      if(state%stage_iteration_count==1)state%trial_stage_first_residual=maximum_residual(1)
    endif
    if(state%stage_iteration_count>=1.and.maximum_residual(1)>state%allowed_residual_growth* &
       max(state%accepted_stage_residual,tiny(1d0)))then
      state%interface_step=0.5d0*state%interface_step
      state%interface_scale=state%accepted_interface_scale
      state%continuation_rollbacks=state%continuation_rollbacks+1
      state%stage_iteration_count=0
      rollback=.true.
      if(state%interface_step<state%minimum_interface_step.or. &
         state%continuation_rollbacks>state%maximum_continuation_rollbacks)then
        state%continuation_failed=.true.
        message='DG DC direct continuation: rollback limit or minimum interface step reached'
      else
        message=''
      endif
      return
    endif
    if(maximum_residual(1)>tolerance)return
    previous_accepted_scale=state%accepted_interface_scale
    state%accepted_interface_scale=state%interface_scale
    state%accepted_stage_residual=maximum_residual(1)
    call append_dg_dc_continuation_history(state,state%accepted_interface_scale)
    if(state%interface_scale<1d0)then
      state%interface_step=min(state%maximum_interface_step,1.5d0*state%interface_step)
      state%interface_scale=min(1d0,state%accepted_interface_scale+state%interface_step)
      state%stage_iteration_count=0
      state%trial_stage_first_residual=huge(1d0)
      advance=.true.
    else if(collective_unmixed)then
      state%direct_ground_state_complete=.true.
      complete=.true.
    else if(previous_accepted_scale<1d0)then
      state%stage_iteration_count=0
      state%trial_stage_first_residual=huge(1d0)
      advance=.true.
    endif
    message=''
  end subroutine evaluate_dg_dc_direct_continuation

  subroutine append_dg_dc_continuation_history(state,accepted_scale)
    type(s_dg_dc_handoff_state),intent(inout)::state
    real(8),intent(in)::accepted_scale
    real(8),allocatable::new_scale(:),new_step(:)
    integer::n
    n=size(state%continuation_scale_history)
    if(state%continuation_scale_history(n)==accepted_scale)return
    allocate(new_scale(n+1),new_step(n+1))
    new_scale(:n)=state%continuation_scale_history
    new_step(:n)=state%continuation_step_history
    new_scale(n+1)=accepted_scale
    new_step(n+1)=accepted_scale-new_scale(n)
    call move_alloc(new_scale,state%continuation_scale_history)
    call move_alloc(new_step,state%continuation_step_history)
  end subroutine append_dg_dc_continuation_history

  subroutine candidate_metric_rank(candidate_box,core_size,buffer,tolerance,rank_count,minimum_pivot)
    complex(8), intent(in) :: candidate_box(:,:,:,:,:)
    integer, intent(in) :: core_size(3),buffer(3)
    real(8), intent(in) :: tolerance
    integer, intent(out) :: rank_count
    real(8), intent(out) :: minimum_pivot
    complex(8), allocatable :: vectors(:,:),q(:)
    real(8), allocatable :: residual_norm2(:)
    integer :: ngrid,nstate,io,jo,selected
    real(8) :: pivot

    ngrid=product(core_size)*size(candidate_box,5)
    nstate=size(candidate_box,4)
    allocate(vectors(ngrid,nstate),q(ngrid),residual_norm2(nstate))
    do io=1,nstate
      vectors(:,io)=reshape(candidate_box(buffer(1)+1:buffer(1)+core_size(1), &
        buffer(2)+1:buffer(2)+core_size(2),buffer(3)+1:buffer(3)+core_size(3),io,:),[ngrid])
      residual_norm2(io)=max(0d0,real(dot_product(vectors(:,io),vectors(:,io)),8))
    end do
    rank_count=0
    minimum_pivot=huge(1d0)
    do
      selected=maxloc(residual_norm2,dim=1)
      pivot=residual_norm2(selected)
      if(.not.ieee_is_finite(pivot) .or. pivot<=tolerance) exit
      q=vectors(:,selected)/sqrt(pivot)
      rank_count=rank_count+1
      minimum_pivot=min(minimum_pivot,pivot)
      residual_norm2(selected)=-1d0
      do io=1,nstate
        if(residual_norm2(io)<0d0) cycle
        do jo=1,2
          vectors(:,io)=vectors(:,io)-dot_product(q,vectors(:,io))*q
        end do
        residual_norm2(io)=max(0d0,real(dot_product(vectors(:,io),vectors(:,io)),8))
      end do
    end do
    if(rank_count==0) minimum_pivot=0d0
  end subroutine candidate_metric_rank

  subroutine move_handoff(source,destination)
    type(s_dg_dc_handoff_state), intent(inout) :: source
    type(s_dg_dc_handoff_state), intent(inout) :: destination
    destination=source
    if (allocated(source%dc_reference_core)) then
      if (allocated(destination%dc_reference_core)) deallocate(destination%dc_reference_core)
      call move_alloc(source%dc_reference_core,destination%dc_reference_core)
    end if
  end subroutine move_handoff

  subroutine move_nodal(source,destination)
    type(s_dg_nodal_common_state), intent(inout) :: source,destination
    call release_dg_nodal_common_state(destination)
    destination=source
    call release_dg_nodal_common_state(source)
  end subroutine move_nodal

  subroutine collective_and(local,communicator,global)
    logical, intent(in) :: local
    integer, intent(in) :: communicator
    logical, intent(out) :: global
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local,global,1,MPI_LOGICAL,MPI_LAND,communicator,ierr)
    if (ierr/=MPI_SUCCESS) global=.false.
#else
    global=local .and. communicator>=0
#endif
  end subroutine collective_and

  subroutine collective_sum(local,global,communicator,ok)
    real(8), intent(in) :: local
    real(8), intent(out) :: global
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local,global,1,MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)
    call collective_and(ierr==MPI_SUCCESS,communicator,ok)
#else
    global=local
    ok=communicator>=0
#endif
  end subroutine collective_sum

  subroutine collective_real_min(local,global,communicator,ok)
    real(8), intent(in) :: local
    real(8), intent(out) :: global
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local,global,1,MPI_DOUBLE_PRECISION,MPI_MIN,communicator,ierr)
    call collective_and(ierr==MPI_SUCCESS,communicator,ok)
#else
    global=local; ok=communicator>=0
#endif
  end subroutine collective_real_min

  subroutine collective_integer_min(local,global,communicator,ok)
    integer, intent(in) :: local,communicator
    integer, intent(out) :: global
    logical, intent(out) :: ok
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local,global,1,MPI_INTEGER,MPI_MIN,communicator,ierr)
    call collective_and(ierr==MPI_SUCCESS,communicator,ok)
#else
    global=local; ok=communicator>=0
#endif
  end subroutine collective_integer_min

  subroutine collective_integer8_sum(local,global,communicator,ok)
    integer(8), intent(in) :: local
    integer(8), intent(out) :: global
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local,global,1,MPI_INTEGER8,MPI_SUM,communicator,ierr)
    call collective_and(ierr==MPI_SUCCESS,communicator,ok)
#else
    global=local; ok=communicator>=0
#endif
  end subroutine collective_integer8_sum

  subroutine collective_integer8_min(local,global,communicator,ok)
    integer(8), intent(in) :: local
    integer(8), intent(out) :: global
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local,global,1,MPI_INTEGER8,MPI_MIN,communicator,ierr)
    call collective_and(ierr==MPI_SUCCESS,communicator,ok)
#else
    global=local; ok=communicator>=0
#endif
  end subroutine collective_integer8_min

  subroutine collective_integer8_max(local,global,communicator,ok)
    integer(8), intent(in) :: local
    integer(8), intent(out) :: global
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local,global,1,MPI_INTEGER8,MPI_MAX,communicator,ierr)
    call collective_and(ierr==MPI_SUCCESS,communicator,ok)
#else
    global=local; ok=communicator>=0
#endif
  end subroutine collective_integer8_max

  subroutine collective_integer_minmax(local,minimum,maximum,communicator,ok)
    integer, intent(in) :: local,communicator
    integer, intent(out) :: minimum,maximum
    logical, intent(out) :: ok
#ifdef USE_MPI
    integer :: ierr1,ierr2
    logical :: stage_ok
    call MPI_Allreduce(local,minimum,1,MPI_INTEGER,MPI_MIN,communicator,ierr1)
    call collective_and(ierr1==MPI_SUCCESS,communicator,stage_ok)
    if(.not.stage_ok) then
      ok=.false.; return
    end if
    call MPI_Allreduce(local,maximum,1,MPI_INTEGER,MPI_MAX,communicator,ierr2)
    call collective_and(ierr2==MPI_SUCCESS,communicator,ok)
#else
    minimum=local; maximum=local; ok=communicator>=0
#endif
  end subroutine collective_integer_minmax

  subroutine collective_real_minmax(local,minimum,maximum,communicator,ok)
    real(8), intent(in) :: local(:)
    real(8), intent(out) :: minimum(:),maximum(:)
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
#ifdef USE_MPI
    integer :: ierr1,ierr2
    logical :: stage_ok
    call MPI_Allreduce(local,minimum,size(local),MPI_DOUBLE_PRECISION,MPI_MIN,communicator,ierr1)
    call collective_and(ierr1==MPI_SUCCESS,communicator,stage_ok)
    if(.not.stage_ok) then
      ok=.false.; return
    end if
    call MPI_Allreduce(local,maximum,size(local),MPI_DOUBLE_PRECISION,MPI_MAX,communicator,ierr2)
    call collective_and(ierr2==MPI_SUCCESS,communicator,ok)
#else
    minimum=local; maximum=local; ok=communicator>=0
#endif
  end subroutine collective_real_minmax
end module dg_dc_handoff
