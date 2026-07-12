!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
!=======================================================================
! Six-face MPI exchange for core-owned nodal DG wavefunctions.
!=======================================================================
module rt_dg_nodal_mpi
#ifdef USE_MPI
  use mpi
#endif
  use rt_dg_nodal_types, only: s_dg_nodal_state, nodal_face_count, nodal_face_slot
  use rt_dg_nodal_halo, only: refresh_nodal_single_fragment_halos
  implicit none
  private
  public :: pack_nodal_face_values, exchange_nodal_face_halos

contains

  subroutine pack_nodal_face_values(state)
    type(s_dg_nodal_state), intent(inout) :: state
    integer :: axis, side, iface, layer, i1, i2, source_index

    if (.not. state%initialized) stop 'nodal DG: face pack before state initialization'
    do axis = 1, 3
      do side = -1, 1, 2
        iface = nodal_face_slot(axis, side)
        call validate_face_before_communication(state, iface, axis, side)
        do layer = 1, state%halo_width
          source_index = merge(layer, state%core_size(axis)-layer+1, side < 0)
          do i2 = 1, state%faces(iface)%tangential_size(2)
            do i1 = 1, state%faces(iface)%tangential_size(1)
              select case(axis)
              case(1)
                state%faces(iface)%send_value(layer,i1,i2,:,:) = &
                  state%psi_core(source_index,i1,i2,:,:)
              case(2)
                state%faces(iface)%send_value(layer,i1,i2,:,:) = &
                  state%psi_core(i1,source_index,i2,:,:)
              case(3)
                state%faces(iface)%send_value(layer,i1,i2,:,:) = &
                  state%psi_core(i1,i2,source_index,:,:)
              end select
            end do
          end do
        end do
      end do
    end do
  end subroutine pack_nodal_face_values

  subroutine exchange_nodal_face_halos(state, communicator)
    type(s_dg_nodal_state), intent(inout) :: state
    integer, intent(in) :: communicator
#ifdef USE_MPI
    integer :: iface, axis, side, peer, send_tag, recv_tag, ierr
    integer :: send_request(nodal_face_count), recv_request(nodal_face_count)

    call pack_nodal_face_values(state)
    do iface = 1, nodal_face_count
      axis = state%faces(iface)%axis
      side = state%faces(iface)%side
      call validate_face_before_communication(state, iface, axis, side)
      peer = state%faces(iface)%neighbor_rank
      recv_tag = nodal_message_tag(state%faces(iface)%neighbor_fragment, axis, -side)
      call MPI_Irecv(state%faces(iface)%recv_value, size(state%faces(iface)%recv_value), &
                     MPI_DOUBLE_COMPLEX, peer, recv_tag, communicator, recv_request(iface), ierr)
      if (ierr /= MPI_SUCCESS) stop 'nodal DG: MPI_Irecv failed'
    end do
    do iface = 1, nodal_face_count
      axis = state%faces(iface)%axis
      side = state%faces(iface)%side
      peer = state%faces(iface)%neighbor_rank
      send_tag = nodal_message_tag(state%fragment, axis, side)
      call MPI_Isend(state%faces(iface)%send_value, size(state%faces(iface)%send_value), &
                     MPI_DOUBLE_COMPLEX, peer, send_tag, communicator, send_request(iface), ierr)
      if (ierr /= MPI_SUCCESS) stop 'nodal DG: MPI_Isend failed'
    end do
    call MPI_Waitall(nodal_face_count, recv_request, MPI_STATUSES_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) stop 'nodal DG: MPI receive wait failed'
    call MPI_Waitall(nodal_face_count, send_request, MPI_STATUSES_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) stop 'nodal DG: MPI send wait failed'
#else
    if (any(state%faces(:)%neighbor_rank /= 0)) &
      stop 'nodal DG: remote face exchange requires an MPI build'
    call refresh_nodal_single_fragment_halos(state)
    if (communicator < 0) stop 'nodal DG: invalid serial communicator'
#endif
  end subroutine exchange_nodal_face_halos

  subroutine validate_face_before_communication(state, iface, axis, side)
    type(s_dg_nodal_state), intent(in) :: state
    integer, intent(in) :: iface, axis, side

    if (iface /= nodal_face_slot(axis, side)) stop 'nodal DG: face slot mismatch before communication'
    if (state%faces(iface)%axis /= axis .or. state%faces(iface)%side /= side) &
      stop 'nodal DG: face identity mismatch before communication'
    if (state%faces(iface)%neighbor_fragment < 1 .or. state%faces(iface)%neighbor_rank < 0) &
      stop 'nodal DG: invalid face neighbor before communication'
    if (.not. allocated(state%faces(iface)%send_value) .or. &
        .not. allocated(state%faces(iface)%recv_value)) &
      stop 'nodal DG: face buffers missing before communication'
    if (size(state%faces(iface)%send_value) /= size(state%faces(iface)%recv_value)) &
      stop 'nodal DG: face send/receive size mismatch'
  end subroutine validate_face_before_communication

  pure integer function nodal_message_tag(fragment, axis, side) result(tag)
    integer, intent(in) :: fragment, axis, side

    tag = (fragment - 1) * nodal_face_count + nodal_face_slot(axis, side)
  end function nodal_message_tag

end module rt_dg_nodal_mpi
