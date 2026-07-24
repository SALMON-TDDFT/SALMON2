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
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dg_nodal_sipg, only: s_dg_nodal_sipg_action,evaluate_dg_nodal_sipg_face
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
  public :: assemble_dg_dc_local_basis_sipg_pair
  public :: compose_dg_dc_local_basis_pair_hamiltonian

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

  subroutine assemble_dg_dc_local_basis_sipg_pair(local_value,local_outward_normal,neighbor_value, &
      neighbor_outward_normal,h_normal,face_weight,penalty_factor,physical_face,canonical_owner, &
      auxiliary_periodic_wrap,periodic_shift,expected_periodic_shift,neighbor_fragment,expected_neighbor_fragment, &
      h_local_local,h_local_neighbor,h_neighbor_local,h_neighbor_neighbor,ok,message)
    complex(8), intent(in) :: local_value(:,:),local_outward_normal(:,:)
    complex(8), intent(in) :: neighbor_value(:,:),neighbor_outward_normal(:,:)
    real(8), intent(in) :: h_normal,face_weight,penalty_factor
    logical, intent(in) :: physical_face,canonical_owner,auxiliary_periodic_wrap
    integer, intent(in) :: periodic_shift(3),expected_periodic_shift(3)
    integer, intent(in) :: neighbor_fragment,expected_neighbor_fragment
    complex(8), intent(out) :: h_local_local(:,:),h_local_neighbor(:,:)
    complex(8), intent(out) :: h_neighbor_local(:,:),h_neighbor_neighbor(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    type(s_dg_nodal_sipg_action) :: action
    integer :: ipoint,itrial,itest,info,nlocal,nneighbor
    complex(8) :: zero

    h_local_local=(0d0,0d0)
    h_local_neighbor=(0d0,0d0)
    h_neighbor_local=(0d0,0d0)
    h_neighbor_neighbor=(0d0,0d0)
    nlocal=size(local_value,2)
    nneighbor=size(neighbor_value,2)
    ok=size(local_value,1)>0 .and. nlocal>0 .and. nneighbor>0 .and. &
      all(shape(local_outward_normal)==shape(local_value)) .and. &
      size(neighbor_value,1)==size(local_value,1) .and. &
      all(shape(neighbor_outward_normal)==shape(neighbor_value)) .and. &
      all(shape(h_local_local)==[nlocal,nlocal]) .and. &
      all(shape(h_local_neighbor)==[nlocal,nneighbor]) .and. &
      all(shape(h_neighbor_local)==[nneighbor,nlocal]) .and. &
      all(shape(h_neighbor_neighbor)==[nneighbor,nneighbor]) .and. &
      h_normal>0d0 .and. face_weight>=0d0 .and. penalty_factor>0d0 .and. &
      ieee_is_finite(h_normal) .and. ieee_is_finite(face_weight) .and. ieee_is_finite(penalty_factor) .and. &
      all(abs(periodic_shift)<=1) .and. count(periodic_shift/=0)<=1 .and. &
      all(abs(expected_periodic_shift)<=1) .and. count(expected_periodic_shift/=0)<=1 .and. &
      neighbor_fragment>0 .and. expected_neighbor_fragment>0
    if(physical_face) ok=ok .and. neighbor_fragment==expected_neighbor_fragment .and. &
      all(periodic_shift==expected_periodic_shift)
    if(.not.ok) then
      message='DG DC local basis: invalid SIPG pair layout or controls'
      return
    end if
    if(.not.physical_face .or. .not.canonical_owner .or. auxiliary_periodic_wrap) then
      message=''
      return
    end if
    zero=(0d0,0d0)
    do ipoint=1,size(local_value,1)
      do itrial=1,nlocal
        call evaluate_dg_nodal_sipg_face(local_value(ipoint,itrial),zero, &
          local_outward_normal(ipoint,itrial),zero,h_normal,face_weight,penalty_factor,action,info)
        if(info/=0) ok=.false.
        do itest=1,nlocal
          h_local_local(itest,itrial)=h_local_local(itest,itrial)+ &
            conjg(local_value(ipoint,itest))*action%total_value(1)+ &
            conjg(local_outward_normal(ipoint,itest))*action%total_normal(1)
        end do
        do itest=1,nneighbor
          h_neighbor_local(itest,itrial)=h_neighbor_local(itest,itrial)+ &
            conjg(neighbor_value(ipoint,itest))*action%total_value(2)- &
            conjg(neighbor_outward_normal(ipoint,itest))*action%total_normal(2)
        end do
      end do
      do itrial=1,nneighbor
        call evaluate_dg_nodal_sipg_face(zero,neighbor_value(ipoint,itrial),zero, &
          -neighbor_outward_normal(ipoint,itrial),h_normal,face_weight,penalty_factor,action,info)
        if(info/=0) ok=.false.
        do itest=1,nlocal
          h_local_neighbor(itest,itrial)=h_local_neighbor(itest,itrial)+ &
            conjg(local_value(ipoint,itest))*action%total_value(1)+ &
            conjg(local_outward_normal(ipoint,itest))*action%total_normal(1)
        end do
        do itest=1,nneighbor
          h_neighbor_neighbor(itest,itrial)=h_neighbor_neighbor(itest,itrial)+ &
            conjg(neighbor_value(ipoint,itest))*action%total_value(2)- &
            conjg(neighbor_outward_normal(ipoint,itest))*action%total_normal(2)
        end do
      end do
    end do
    ok=ok .and. all(ieee_is_finite(real(h_local_local,8))) .and. &
      all(ieee_is_finite(aimag(h_local_local))) .and. &
      all(ieee_is_finite(real(h_local_neighbor,8))) .and. &
      all(ieee_is_finite(aimag(h_local_neighbor))) .and. &
      all(ieee_is_finite(real(h_neighbor_local,8))) .and. &
      all(ieee_is_finite(aimag(h_neighbor_local))) .and. &
      all(ieee_is_finite(real(h_neighbor_neighbor,8))) .and. &
      all(ieee_is_finite(aimag(h_neighbor_neighbor)))
    if(ok) then
      message=''
    else
      h_local_local=(0d0,0d0)
      h_local_neighbor=(0d0,0d0)
      h_neighbor_local=(0d0,0d0)
      h_neighbor_neighbor=(0d0,0d0)
      message='DG DC local basis: SIPG pair evaluation failed'
    end if
  end subroutine assemble_dg_dc_local_basis_sipg_pair

  subroutine compose_dg_dc_local_basis_pair_hamiltonian(volume_left,volume_right,h_local_local, &
      h_local_neighbor,h_neighbor_local,h_neighbor_neighbor,lambda,hamiltonian,ok,message)
    complex(8), intent(in) :: volume_left(:,:),volume_right(:,:),h_local_local(:,:),h_local_neighbor(:,:)
    complex(8), intent(in) :: h_neighbor_local(:,:),h_neighbor_neighbor(:,:)
    real(8), intent(in) :: lambda
    complex(8), intent(out) :: hamiltonian(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: nleft,nright
    nleft=size(volume_left,1)
    nright=size(volume_right,1)
    hamiltonian=(0d0,0d0)
    ok=nleft>0 .and. nright>0 .and. size(volume_left,2)==nleft .and. &
      size(volume_right,2)==nright .and. all(shape(h_local_local)==[nleft,nleft]) .and. &
      all(shape(h_local_neighbor)==[nleft,nright]) .and. &
      all(shape(h_neighbor_local)==[nright,nleft]) .and. &
      all(shape(h_neighbor_neighbor)==[nright,nright]) .and. &
      all(shape(hamiltonian)==[nleft+nright,nleft+nright]) .and. &
      ieee_is_finite(lambda) .and. lambda>=0d0 .and. lambda<=1d0
    if(.not.ok) then
      message='DG DC local basis: invalid pair Hamiltonian layout or lambda'
      return
    end if
    hamiltonian(1:nleft,1:nleft)=volume_left+lambda*h_local_local
    hamiltonian(1:nleft,nleft+1:nleft+nright)=lambda*h_local_neighbor
    hamiltonian(nleft+1:nleft+nright,1:nleft)=lambda*h_neighbor_local
    hamiltonian(nleft+1:nleft+nright,nleft+1:nleft+nright)=volume_right+lambda*h_neighbor_neighbor
    ok=all(ieee_is_finite(real(hamiltonian,8))) .and. all(ieee_is_finite(aimag(hamiltonian)))
    if(ok) then
      message=''
    else
      hamiltonian=(0d0,0d0)
      message='DG DC local basis: nonfinite pair Hamiltonian'
    end if
  end subroutine compose_dg_dc_local_basis_pair_hamiltonian

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
