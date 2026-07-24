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
    integer :: minimum_iteration=0, accepted_iteration=0
    integer :: candidate_orbitals_per_atom=0, candidate_count=0, metric_rank=0
    real(8) :: tolerance=0.0d0, metric_rank_tolerance=0.0d0
    real(8) :: residual=huge(1.0d0), charge_error=huge(1.0d0)
    real(8) :: minimum_metric_pivot=0.0d0, span_escape_norm=0.0d0
    integer(8) :: geometry_fingerprint=0_8, operator_fingerprint=0_8
    complex(8), allocatable :: dc_reference_core(:,:,:,:,:)
    real(8), allocatable :: density_snapshot(:)
    real(8), allocatable :: potential_snapshot(:,:)
    integer(8), allocatable :: snapshot_owner(:)
    logical :: density_potential_preserved=.false.
    logical :: mixing_history_discarded=.false.
  end type s_dg_dc_handoff_state

  type(s_dg_dc_handoff_state), public, save :: dg_dc_handoff_runtime
  type(s_dg_nodal_common_state), public, save :: dg_dc_nodal_runtime

  public :: initialize_dg_dc_handoff, evaluate_dg_dc_handoff
  public :: materialize_dg_dc_candidates, update_dg_dc_nodal_state
  public :: preserve_dg_dc_density_potential, mark_dg_dc_mixing_discarded

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

  subroutine evaluate_dg_dc_handoff(state,iteration,residual,charge,expected_charge,fragment_solved, &
                                     lcfo_started,wannier_started,window_truncated,communicator,accept,message)
    type(s_dg_dc_handoff_state), intent(inout) :: state
    integer, intent(in) :: iteration,communicator
    real(8), intent(in) :: residual,charge,expected_charge
    logical, intent(in) :: fragment_solved,lcfo_started,wannier_started,window_truncated
    logical, intent(out) :: accept
    character(*), intent(out) :: message
    logical :: local_accept,metadata_ok,real_metadata_ok,integer_metadata_ok
    real(8) :: charge_tolerance,minimum_values(3),maximum_values(3),local_values(3)
    integer :: minimum_iteration_seen,maximum_iteration_seen

    accept=.false.; message=''
    if (.not.state%enabled .or. state%accepted) return
    charge_tolerance=max(state%tolerance,100.0d0*epsilon(1.0d0)*max(1.0d0,abs(expected_charge)))
    local_accept=iteration >= state%minimum_iteration .and. ieee_is_finite(residual) .and. &
      ieee_is_finite(charge) .and. ieee_is_finite(expected_charge) .and. residual <= state%tolerance .and. &
      abs(charge-expected_charge) <= charge_tolerance .and. fragment_solved .and. &
      .not.lcfo_started .and. .not.wannier_started .and. .not.window_truncated
    call collective_and(local_accept,communicator,accept)
    if (.not.accept) then
      message='DG DC handoff: gate not satisfied collectively'
      return
    end if
    local_values=[residual,charge,expected_charge]
    call collective_real_minmax(local_values,minimum_values,maximum_values,communicator,real_metadata_ok)
    call collective_integer_minmax(iteration,minimum_iteration_seen,maximum_iteration_seen,communicator,integer_metadata_ok)
    metadata_ok=real_metadata_ok .and. integer_metadata_ok
    if(.not.metadata_ok .or. minimum_iteration_seen/=maximum_iteration_seen .or. &
       any(minimum_values/=maximum_values)) then
      accept=.false.
      message='DG DC handoff: rank-disagreeing gate diagnostics'
      return
    end if
    state%accepted=.true.
    state%accepted_iteration=minimum_iteration_seen
    state%residual=minimum_values(1)
    state%charge_error=abs(minimum_values(2)-minimum_values(3))
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
      configured_count==natom*state%candidate_orbitals_per_atom .and. &
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

  subroutine preserve_dg_dc_density_potential(state,density,potential,owner,communicator,ok,message)
    type(s_dg_dc_handoff_state), intent(inout) :: state
    real(8), intent(in) :: density(:),potential(:,:)
    integer(8), intent(in) :: owner(:)
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    logical :: local_ok,ownership_ok
    integer(8) :: local_count,global_count,local_minimum,global_minimum,local_maximum,global_maximum
    local_ok=state%accepted .and. size(density)>0 .and. size(potential,1)==size(density) .and. &
      size(potential,2)>0 .and. all(ieee_is_finite(density)) .and. all(ieee_is_finite(potential))
    local_ok=local_ok .and. size(owner)==size(density) .and. all(owner>0_8)
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
    if(allocated(state%density_snapshot)) deallocate(state%density_snapshot)
    if(allocated(state%potential_snapshot)) deallocate(state%potential_snapshot)
    if(allocated(state%snapshot_owner)) deallocate(state%snapshot_owner)
    allocate(state%density_snapshot,source=density)
    allocate(state%potential_snapshot,source=potential)
    allocate(state%snapshot_owner,source=owner)
    state%density_potential_preserved=.true.
    message=''
  end subroutine preserve_dg_dc_density_potential

  subroutine mark_dg_dc_mixing_discarded(state)
    type(s_dg_dc_handoff_state), intent(inout) :: state
    state%mixing_history_discarded=.true.
  end subroutine mark_dg_dc_mixing_discarded

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
