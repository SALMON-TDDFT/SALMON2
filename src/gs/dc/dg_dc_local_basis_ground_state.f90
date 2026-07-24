!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
module dg_dc_local_basis_ground_state
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private

  type, public :: s_dg_dc_local_basis_layout
    integer :: fragment_id=0
    integer :: local_basis_count=0
    integer :: global_basis_count=0
    integer :: global_band_count=0
    integer(8) :: geometry_fingerprint=0_8
    integer(8) :: operator_fingerprint=0_8
    integer, allocatable :: basis_offsets(:)
    integer, allocatable :: fragment_ids(:)
  end type s_dg_dc_local_basis_layout

  public :: initialize_dg_dc_local_basis_layout
  public :: release_dg_dc_local_basis_layout

contains

  subroutine initialize_dg_dc_local_basis_layout(layout,fragment_id,local_basis_count,global_band_count, &
      geometry_fingerprint,operator_fingerprint,communicator,ok,message)
    type(s_dg_dc_local_basis_layout), intent(inout) :: layout
    integer, intent(in) :: fragment_id,local_basis_count,global_band_count,communicator
    integer(8), intent(in) :: geometry_fingerprint,operator_fingerprint
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: nproc,i,ierr,allocation_status
    integer(8) :: running_basis_count
    integer, allocatable :: basis_counts(:),band_counts(:)
    integer(8), allocatable :: geometry_values(:),operator_values(:)
    logical, allocatable :: seen_fragment(:)
    logical :: valid,stage_ok,all_stage_ok

    call release_dg_dc_local_basis_layout(layout)
    nproc=1
#ifdef USE_MPI
    call MPI_Comm_size(communicator,nproc,ierr)
    if(ierr/=MPI_SUCCESS) then
      ok=.false.; message='DG DC local basis: communicator size query failed'; return
    end if
#else
    if(communicator<0) then
      ok=.false.; message='DG DC local basis: invalid serial communicator'; return
    end if
#endif
    allocate(basis_counts(nproc),band_counts(nproc),layout%fragment_ids(nproc),geometry_values(nproc), &
      operator_values(nproc),layout%basis_offsets(0:nproc),seen_fragment(nproc),stat=allocation_status)
    stage_ok=allocation_status==0
    call collective_logical_and(stage_ok,communicator,all_stage_ok)
    if(.not.all_stage_ok) then
      if(allocated(basis_counts)) deallocate(basis_counts)
      if(allocated(band_counts)) deallocate(band_counts)
      if(allocated(geometry_values)) deallocate(geometry_values)
      if(allocated(operator_values)) deallocate(operator_values)
      if(allocated(seen_fragment)) deallocate(seen_fragment)
      call release_dg_dc_local_basis_layout(layout)
      ok=.false.; message='DG DC local basis: layout allocation failed collectively'; return
    end if
#ifdef USE_MPI
    call MPI_Allgather(local_basis_count,1,MPI_INTEGER,basis_counts,1,MPI_INTEGER,communicator,ierr)
    call collective_logical_and(ierr==MPI_SUCCESS,communicator,all_stage_ok)
    if(all_stage_ok) call MPI_Allgather(global_band_count,1,MPI_INTEGER,band_counts,1,MPI_INTEGER, &
      communicator,ierr)
    call collective_logical_and(all_stage_ok .and. ierr==MPI_SUCCESS,communicator,all_stage_ok)
    if(all_stage_ok) call MPI_Allgather(fragment_id,1,MPI_INTEGER,layout%fragment_ids,1,MPI_INTEGER, &
      communicator,ierr)
    call collective_logical_and(all_stage_ok .and. ierr==MPI_SUCCESS,communicator,all_stage_ok)
    if(all_stage_ok) call MPI_Allgather(geometry_fingerprint,1,MPI_INTEGER8,geometry_values,1, &
      MPI_INTEGER8,communicator,ierr)
    call collective_logical_and(all_stage_ok .and. ierr==MPI_SUCCESS,communicator,all_stage_ok)
    if(all_stage_ok) call MPI_Allgather(operator_fingerprint,1,MPI_INTEGER8,operator_values,1, &
      MPI_INTEGER8,communicator,ierr)
    call collective_logical_and(all_stage_ok .and. ierr==MPI_SUCCESS,communicator,all_stage_ok)
    if(.not.all_stage_ok) then
      deallocate(basis_counts,band_counts,geometry_values,operator_values,seen_fragment)
      call release_dg_dc_local_basis_layout(layout)
      ok=.false.; message='DG DC local basis: layout gather failed collectively'; return
    end if
#else
    basis_counts=[local_basis_count]
    band_counts=[global_band_count]
    layout%fragment_ids=[fragment_id]
    geometry_values=[geometry_fingerprint]
    operator_values=[operator_fingerprint]
#endif
    valid=all(basis_counts>0) .and. all(band_counts==band_counts(1)) .and. band_counts(1)>0 .and. &
      all(layout%fragment_ids>0) .and. all(layout%fragment_ids<=nproc) .and. &
      all(geometry_values==geometry_values(1)) .and. &
      all(operator_values==operator_values(1))
    seen_fragment=.false.
    do i=1,nproc
      if(layout%fragment_ids(i)>=1 .and. layout%fragment_ids(i)<=nproc) then
        valid=valid .and. .not.seen_fragment(layout%fragment_ids(i))
        seen_fragment(layout%fragment_ids(i))=.true.
      end if
    end do
    valid=valid .and. all(seen_fragment)
    layout%basis_offsets(0)=0
    running_basis_count=0_8
    do i=1,nproc
      if(int(basis_counts(i),8)>int(huge(1),8)-running_basis_count) then
        valid=.false.
        layout%basis_offsets(i)=0
      else
        running_basis_count=running_basis_count+int(basis_counts(i),8)
        layout%basis_offsets(i)=int(running_basis_count)
      end if
    end do
    if(.not.valid) then
      call release_dg_dc_local_basis_layout(layout)
      ok=.false.; message='DG DC local basis: invalid or rank-disagreeing layout'; return
    end if
    layout%fragment_id=fragment_id
    layout%local_basis_count=local_basis_count
    layout%global_basis_count=layout%basis_offsets(nproc)
    layout%global_band_count=global_band_count
    layout%geometry_fingerprint=geometry_fingerprint
    layout%operator_fingerprint=operator_fingerprint
    ok=.true.
    message=''
    deallocate(basis_counts,band_counts,geometry_values,operator_values,seen_fragment)
  end subroutine initialize_dg_dc_local_basis_layout

  subroutine release_dg_dc_local_basis_layout(layout)
    type(s_dg_dc_local_basis_layout), intent(inout) :: layout
    if(allocated(layout%basis_offsets)) deallocate(layout%basis_offsets)
    if(allocated(layout%fragment_ids)) deallocate(layout%fragment_ids)
    layout%fragment_id=0
    layout%local_basis_count=0
    layout%global_basis_count=0
    layout%global_band_count=0
    layout%geometry_fingerprint=0_8
    layout%operator_fingerprint=0_8
  end subroutine release_dg_dc_local_basis_layout

  subroutine collective_logical_and(local_value,communicator,global_value)
    logical, intent(in) :: local_value
    integer, intent(in) :: communicator
    logical, intent(out) :: global_value
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local_value,global_value,1,MPI_LOGICAL,MPI_LAND,communicator,ierr)
    if(ierr/=MPI_SUCCESS) global_value=.false.
#else
    global_value=local_value
    if(communicator<0) global_value=.false.
#endif
  end subroutine collective_logical_and
end module dg_dc_local_basis_ground_state
