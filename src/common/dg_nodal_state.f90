!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
module dg_nodal_state
#ifdef USE_MPI
  use mpi
#endif
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  integer, parameter, public :: nodal_face_count = 6

  type, public :: s_dg_nodal_face
    integer :: axis = 0
    integer :: side = 0
    integer :: neighbor_fragment = 0
    integer :: neighbor_rank = -1
    integer :: width = 0
    integer :: tangential_size(2) = 0
    complex(8), allocatable :: send_value(:,:,:,:,:)
    complex(8), allocatable :: recv_value(:,:,:,:,:)
    complex(8), allocatable :: send_normal(:,:,:,:,:)
    complex(8), allocatable :: recv_normal(:,:,:,:,:)
  end type s_dg_nodal_face

  type, public :: s_dg_nodal_common_state
    logical :: enabled = .false.
    logical :: initialized = .false.
    logical :: ground_state_ready = .false.
    integer :: fragment = 0
    integer :: nstate = 0
    integer :: nspin = 0
    integer :: halo_width = 0
    integer :: core_size(3) = 0
    integer :: box_size(3) = 0
    integer :: buffer(3) = 0
    integer(8) :: geometry_fingerprint = 0
    integer(8) :: operator_fingerprint = 0
    real(8) :: dg_ground_state_residual = huge(1.0d0)
    complex(8), allocatable :: psi_core(:,:,:,:,:)
    complex(8), allocatable :: hpsi_core(:,:,:,:,:)
    real(8), allocatable :: occupations(:,:)
    integer(8), allocatable :: core_owner(:,:,:)
    type(s_dg_nodal_face), allocatable :: faces(:)
  end type s_dg_nodal_common_state

  public :: initialize_dg_nodal_common_state
  public :: initialize_dg_nodal_common_state_local
  public :: release_dg_nodal_common_state
  public :: validate_dg_nodal_common_state_mpi
  public :: dg_nodal_face_slot, allocate_dg_nodal_face_buffers

contains

  pure integer function dg_nodal_face_slot(axis, side) result(iface)
    integer, intent(in) :: axis, side
    if (axis < 1 .or. axis > 3 .or. abs(side) /= 1) then
      iface = 0
    else
      iface = 2 * (axis - 1) + merge(1, 2, side < 0)
    end if
  end function dg_nodal_face_slot

  subroutine initialize_dg_nodal_common_state(state, fragment, core_size, box_size, buffer, halo_width, &
                                               nstate, nspin, core_owner, geometry_fingerprint, &
                                               operator_fingerprint, communicator, ok, message)
    type(s_dg_nodal_common_state), intent(inout) :: state
    integer, intent(in) :: fragment, core_size(3), box_size(3), buffer(3), halo_width, nstate, nspin
    integer(8), intent(in) :: core_owner(:,:,:), geometry_fingerprint, operator_fingerprint
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    type(s_dg_nodal_common_state) :: candidate
    logical :: local_ok

    call build_candidate(candidate, fragment, core_size, box_size, buffer, halo_width, nstate, nspin, &
                         core_owner, geometry_fingerprint, operator_fingerprint, local_ok, message)
    call collective_logical_and(local_ok, communicator, ok)
    if (.not.ok) then
      if (local_ok) message = 'nodal DG: another rank has invalid common-state metadata'
      call release_dg_nodal_common_state(candidate)
      return
    end if
    call validate_global_metadata(candidate, communicator, ok, message)
    if (.not.ok) then
      call release_dg_nodal_common_state(candidate)
      return
    end if
    call validate_unique_core_ownership(candidate%core_owner, communicator, ok, message)
    if (.not.ok) then
      call release_dg_nodal_common_state(candidate)
      return
    end if
    call move_state(candidate, state)
    message = ''
  end subroutine initialize_dg_nodal_common_state

  subroutine initialize_dg_nodal_common_state_local(state, fragment, core_size, box_size, buffer, halo_width, &
                                                     nstate, nspin, core_owner, geometry_fingerprint, &
                                                     operator_fingerprint, ok, message)
    type(s_dg_nodal_common_state), intent(inout) :: state
    integer, intent(in) :: fragment, core_size(3), box_size(3), buffer(3), halo_width, nstate, nspin
    integer(8), intent(in) :: core_owner(:,:,:), geometry_fingerprint, operator_fingerprint
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    type(s_dg_nodal_common_state) :: candidate

    call build_candidate(candidate, fragment, core_size, box_size, buffer, halo_width, nstate, nspin, &
                         core_owner, geometry_fingerprint, operator_fingerprint, ok, message)
    if (.not.ok) return
    call move_state(candidate, state)
    message = ''
  end subroutine initialize_dg_nodal_common_state_local

  subroutine build_candidate(candidate, fragment, core_size, box_size, buffer, halo_width, nstate, nspin, &
                             core_owner, geometry_fingerprint, operator_fingerprint, ok, message)
    type(s_dg_nodal_common_state), intent(out) :: candidate
    integer, intent(in) :: fragment, core_size(3), box_size(3), buffer(3), halo_width, nstate, nspin
    integer(8), intent(in) :: core_owner(:,:,:), geometry_fingerprint, operator_fingerprint
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: axis, side, iface

    ok = .false.
    message = ''
    if (fragment < 1 .or. any(core_size < 1) .or. any(buffer < 0) .or. &
        any(box_size /= core_size + 2*buffer)) then
      message = 'nodal DG: invalid core-buffer geometry'
      return
    end if
    if (halo_width < 1 .or. nstate < 1 .or. nspin < 1) then
      message = 'nodal DG: invalid common-state dimensions'
      return
    end if
    if (any(shape(core_owner) /= core_size) .or. any(core_owner < 1_8)) then
      message = 'nodal DG: invalid explicit core ownership'
      return
    end if
    if (geometry_fingerprint == 0_8 .or. operator_fingerprint == 0_8) then
      message = 'nodal DG: missing provenance fingerprint'
      return
    end if

    allocate(candidate%psi_core(core_size(1),core_size(2),core_size(3),nstate,nspin))
    allocate(candidate%hpsi_core(core_size(1),core_size(2),core_size(3),nstate,nspin))
    allocate(candidate%occupations(nstate,nspin))
    allocate(candidate%core_owner, source=core_owner)
    allocate(candidate%faces(nodal_face_count))
    candidate%psi_core = (0.0d0,0.0d0)
    candidate%hpsi_core = (0.0d0,0.0d0)
    candidate%occupations = 0.0d0
    candidate%fragment = fragment
    candidate%core_size = core_size
    candidate%box_size = box_size
    candidate%buffer = buffer
    candidate%halo_width = halo_width
    candidate%nstate = nstate
    candidate%nspin = nspin
    candidate%geometry_fingerprint = geometry_fingerprint
    candidate%operator_fingerprint = operator_fingerprint
    do axis=1,3
      do side=-1,1,2
        iface=dg_nodal_face_slot(axis,side)
        candidate%faces(iface)%axis=axis
        candidate%faces(iface)%side=side
      end do
    end do
    candidate%initialized = .true.
    ok = .true.
  end subroutine build_candidate

  subroutine validate_dg_nodal_common_state_mpi(state, communicator, ok, message)
    type(s_dg_nodal_common_state), intent(in) :: state
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    logical :: local_ok

    local_ok = state%initialized .and. allocated(state%psi_core) .and. allocated(state%hpsi_core) .and. &
      allocated(state%occupations) .and. allocated(state%core_owner) .and. allocated(state%faces) .and. &
      state%geometry_fingerprint /= 0_8 .and. state%operator_fingerprint /= 0_8
    if (local_ok) local_ok = all(shape(state%psi_core) == [state%core_size,state%nstate,state%nspin]) .and. &
      all(shape(state%hpsi_core) == [state%core_size,state%nstate,state%nspin]) .and. &
      all(shape(state%occupations) == [state%nstate,state%nspin]) .and. &
      all(shape(state%core_owner) == state%core_size) .and. size(state%faces) == nodal_face_count
    if (local_ok) local_ok = all(ieee_is_finite(real(state%psi_core,8))) .and. &
      all(ieee_is_finite(aimag(state%psi_core))) .and. all(ieee_is_finite(real(state%hpsi_core,8))) .and. &
      all(ieee_is_finite(aimag(state%hpsi_core))) .and. all(ieee_is_finite(state%occupations))
    call collective_logical_and(local_ok, communicator, ok)
    if (ok) then
      message = ''
    else
      message = 'nodal DG: common-state validation failed collectively'
    end if
  end subroutine validate_dg_nodal_common_state_mpi

  subroutine validate_global_metadata(state, communicator, ok, message)
    type(s_dg_nodal_common_state), intent(in) :: state
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
    character(*), intent(out) :: message
#ifdef USE_MPI
    integer(8) :: local_values(5), minimum(5), maximum(5)
    integer :: ierr
    logical :: stage_ok
    local_values=[int(state%nstate,8),int(state%nspin,8),int(state%halo_width,8), &
                  state%geometry_fingerprint,state%operator_fingerprint]
    call MPI_Allreduce(local_values,minimum,size(local_values),MPI_INTEGER8,MPI_MIN,communicator,ierr)
    call mpi_stage_agrees(ierr == MPI_SUCCESS,communicator,stage_ok)
    if (.not.stage_ok) then
      ok=.false.; message='nodal DG: metadata minimum reduction failed collectively'; return
    end if
    call MPI_Allreduce(local_values,maximum,size(local_values),MPI_INTEGER8,MPI_MAX,communicator,ierr)
    call mpi_stage_agrees(ierr == MPI_SUCCESS,communicator,stage_ok)
    if (.not.stage_ok) then
      ok=.false.; message='nodal DG: metadata maximum reduction failed collectively'; return
    end if
    ok=all(minimum == maximum)
#else
    ok=communicator >= 0
#endif
    if (ok) then
      message=''
    else
      message='nodal DG: rank-disagreeing common-state metadata'
    end if
  end subroutine validate_global_metadata

  subroutine validate_unique_core_ownership(owner, communicator, ok, message)
    integer(8), intent(in) :: owner(:,:,:)
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
    character(*), intent(out) :: message
#ifdef USE_MPI
    integer :: ierr, nproc, i, destination, nreceive
    integer, allocatable :: send_counts(:), receive_counts(:), send_displacements(:), receive_displacements(:)
    integer, allocatable :: cursor(:)
    integer(8), allocatable :: local_owner(:), send_owner(:), receive_owner(:)
    logical :: stage_ok, local_unique

    call MPI_Comm_size(communicator,nproc,ierr)
    call mpi_stage_agrees(ierr == MPI_SUCCESS,communicator,stage_ok)
    if (.not.stage_ok) then
      ok=.false.; message='nodal DG: communicator-size query failed collectively'; return
    end if
    allocate(send_counts(nproc),receive_counts(nproc),send_displacements(nproc),receive_displacements(nproc))
    send_counts=0
    allocate(local_owner(size(owner)))
    local_owner=reshape(owner,[size(owner)])
    do i=1,size(local_owner)
      destination=int(modulo(local_owner(i)-1_8,int(nproc,8)))+1
      send_counts(destination)=send_counts(destination)+1
    end do
    call MPI_Alltoall(send_counts,1,MPI_INTEGER,receive_counts,1,MPI_INTEGER,communicator,ierr)
    call mpi_stage_agrees(ierr == MPI_SUCCESS,communicator,stage_ok)
    if (.not.stage_ok) then
      ok=.false.; message='nodal DG: ownership-count exchange failed collectively'
      deallocate(local_owner,send_counts,receive_counts,send_displacements,receive_displacements)
      return
    end if
    send_displacements(1)=0
    receive_displacements(1)=0
    do i=2,nproc
      send_displacements(i)=send_displacements(i-1)+send_counts(i-1)
      receive_displacements(i)=receive_displacements(i-1)+receive_counts(i-1)
    end do
    allocate(send_owner(size(local_owner)),cursor(nproc))
    cursor=send_displacements+1
    do i=1,size(local_owner)
      destination=int(modulo(local_owner(i)-1_8,int(nproc,8)))+1
      send_owner(cursor(destination))=local_owner(i)
      cursor(destination)=cursor(destination)+1
    end do
    nreceive=sum(receive_counts)
    allocate(receive_owner(nreceive))
    call MPI_Alltoallv(send_owner,send_counts,send_displacements,MPI_INTEGER8,receive_owner,receive_counts, &
                       receive_displacements,MPI_INTEGER8,communicator,ierr)
    call mpi_stage_agrees(ierr == MPI_SUCCESS,communicator,stage_ok)
    if (stage_ok) then
      call sort_integer8(receive_owner)
      local_unique=.true.
      do i=2,nreceive
        if (receive_owner(i) == receive_owner(i-1)) local_unique=.false.
      end do
      call collective_logical_and(local_unique,communicator,ok)
    else
      ok=.false.
    end if
    deallocate(local_owner,send_owner,receive_owner,cursor,send_counts,receive_counts, &
               send_displacements,receive_displacements)
#else
    integer :: i, j
    integer(8), allocatable :: flat(:)
    allocate(flat(size(owner)))
    flat=reshape(owner,[size(owner)])
    ok=.true.
    do i=1,size(flat)-1
      do j=i+1,size(flat)
        if (flat(i) == flat(j)) ok=.false.
      end do
    end do
    deallocate(flat)
    if (communicator < 0) ok=.false.
#endif
    if (ok) then
      message=''
    else
      message='nodal DG: duplicate core ownership'
    end if
  end subroutine validate_unique_core_ownership

  subroutine sort_integer8(values)
    integer(8), intent(inout) :: values(:)
    integer :: root, tail, child
    integer(8) :: saved
    if (size(values) < 2) return
    do root=size(values)/2,1,-1
      saved=values(root)
      tail=root
      do while (2*tail <= size(values))
        child=2*tail
        if (child < size(values)) then
          if (values(child) < values(child+1)) child=child+1
        end if
        if (saved >= values(child)) exit
        values(tail)=values(child)
        tail=child
      end do
      values(tail)=saved
    end do
    do root=size(values),2,-1
      saved=values(root)
      values(root)=values(1)
      tail=1
      do while (2*tail < root)
        child=2*tail
        if (child < root-1) then
          if (values(child) < values(child+1)) child=child+1
        end if
        if (saved >= values(child)) exit
        values(tail)=values(child)
        tail=child
      end do
      values(tail)=saved
    end do
  end subroutine sort_integer8

#ifdef USE_MPI
  subroutine mpi_stage_agrees(local_ok,communicator,global_ok)
    logical, intent(in) :: local_ok
    integer, intent(in) :: communicator
    logical, intent(out) :: global_ok
    integer :: ierr
    call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,communicator,ierr)
    if (ierr /= MPI_SUCCESS) global_ok=.false.
  end subroutine mpi_stage_agrees
#endif

  subroutine collective_logical_and(local_value, communicator, global_value)
    logical, intent(in) :: local_value
    integer, intent(in) :: communicator
    logical, intent(out) :: global_value
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local_value,global_value,1,MPI_LOGICAL,MPI_LAND,communicator,ierr)
    if (ierr /= MPI_SUCCESS) global_value=.false.
#else
    global_value=local_value .and. communicator >= 0
#endif
  end subroutine collective_logical_and

  subroutine allocate_dg_nodal_face_buffers(face, core_size, halo_width, nstate, nspin, ok)
    type(s_dg_nodal_face), intent(inout) :: face
    integer, intent(in) :: core_size(3), halo_width, nstate, nspin
    logical, intent(out) :: ok
    integer :: tangential(2)

    ok=.false.
    if (face%axis < 1 .or. face%axis > 3 .or. abs(face%side) /= 1) return
    if (any(core_size < 1) .or. halo_width < 1 .or. nstate < 1 .or. nspin < 1) return
    select case(face%axis)
    case(1); tangential=[core_size(2),core_size(3)]
    case(2); tangential=[core_size(1),core_size(3)]
    case(3); tangential=[core_size(1),core_size(2)]
    end select
    if (allocated(face%send_value)) deallocate(face%send_value)
    if (allocated(face%recv_value)) deallocate(face%recv_value)
    if (allocated(face%send_normal)) deallocate(face%send_normal)
    if (allocated(face%recv_normal)) deallocate(face%recv_normal)
    allocate(face%send_value(halo_width,tangential(1),tangential(2),nstate,nspin))
    allocate(face%recv_value(halo_width,tangential(1),tangential(2),nstate,nspin))
    allocate(face%send_normal(halo_width,tangential(1),tangential(2),nstate,nspin))
    allocate(face%recv_normal(halo_width,tangential(1),tangential(2),nstate,nspin))
    face%send_value=(0.0d0,0.0d0); face%recv_value=(0.0d0,0.0d0)
    face%send_normal=(0.0d0,0.0d0); face%recv_normal=(0.0d0,0.0d0)
    face%width=halo_width
    face%tangential_size=tangential
    ok=.true.
  end subroutine allocate_dg_nodal_face_buffers

  subroutine release_dg_nodal_common_state(state)
    type(s_dg_nodal_common_state), intent(inout) :: state
    type(s_dg_nodal_common_state) :: empty
    call move_state(empty,state)
  end subroutine release_dg_nodal_common_state

  subroutine move_state(source, destination)
    type(s_dg_nodal_common_state), intent(inout) :: source
    type(s_dg_nodal_common_state), intent(inout) :: destination
    call release_allocations(destination)
    destination%enabled=source%enabled
    destination%initialized=source%initialized
    destination%ground_state_ready=source%ground_state_ready
    destination%fragment=source%fragment
    destination%nstate=source%nstate
    destination%nspin=source%nspin
    destination%halo_width=source%halo_width
    destination%core_size=source%core_size
    destination%box_size=source%box_size
    destination%buffer=source%buffer
    destination%geometry_fingerprint=source%geometry_fingerprint
    destination%operator_fingerprint=source%operator_fingerprint
    destination%dg_ground_state_residual=source%dg_ground_state_residual
    call move_alloc(source%psi_core,destination%psi_core)
    call move_alloc(source%hpsi_core,destination%hpsi_core)
    call move_alloc(source%occupations,destination%occupations)
    call move_alloc(source%core_owner,destination%core_owner)
    call move_alloc(source%faces,destination%faces)
    source%initialized=.false.
  end subroutine move_state

  subroutine release_allocations(state)
    type(s_dg_nodal_common_state), intent(inout) :: state
    if (allocated(state%psi_core)) deallocate(state%psi_core)
    if (allocated(state%hpsi_core)) deallocate(state%hpsi_core)
    if (allocated(state%occupations)) deallocate(state%occupations)
    if (allocated(state%core_owner)) deallocate(state%core_owner)
    if (allocated(state%faces)) deallocate(state%faces)
  end subroutine release_allocations

end module dg_nodal_state
