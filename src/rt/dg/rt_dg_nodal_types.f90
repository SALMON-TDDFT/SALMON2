!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
!=======================================================================
! Core-owned real-space state and face halos for matrix-free DG RT.
!=======================================================================
module rt_dg_nodal_types
  implicit none
  private

  integer, parameter, public :: nodal_face_count = 6

  type :: nodal_face_info
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
  end type nodal_face_info

  type :: s_dg_nodal_state
    logical :: enabled = .false.
    logical :: initialized = .false.
    integer :: fragment = 0
    integer :: core_size(3) = 0
    integer :: halo_width = 0
    integer :: nstate = 0
    integer :: nspin = 0
    logical :: dg_ground_state_ready = .false.
    real(8) :: dg_ground_state_residual = huge(1.0d0)
    ! Only core points are propagated. Face buffers are refreshed for every
    ! Hamiltonian action in the Taylor expansion.
    complex(8), allocatable :: psi_core(:,:,:,:,:)
    complex(8), allocatable :: hpsi_core(:,:,:,:,:)
    type(nodal_face_info) :: faces(nodal_face_count)
  end type s_dg_nodal_state

  public :: nodal_face_info, s_dg_nodal_state, nodal_face_slot
  public :: allocate_nodal_state, allocate_nodal_face_buffers
  public :: accept_nodal_dg_ground_state

contains

  pure integer function nodal_face_slot(axis, side) result(iface)
    integer, intent(in) :: axis, side

    if (axis < 1 .or. axis > 3 .or. abs(side) /= 1) then
      iface = 0
    else
      iface = 2 * (axis - 1) + merge(1, 2, side < 0)
    end if
  end function nodal_face_slot

  subroutine allocate_nodal_state(state, fragment, core_size, halo_width, nstate, nspin)
    type(s_dg_nodal_state), intent(inout) :: state
    integer, intent(in) :: fragment, core_size(3), halo_width, nstate, nspin
    integer :: axis, side, iface

    if (fragment < 1 .or. any(core_size < 1) .or. halo_width < 1 .or. &
        nstate < 1 .or. nspin < 1) stop 'nodal DG: invalid state dimensions'
    if (allocated(state%psi_core)) deallocate(state%psi_core)
    if (allocated(state%hpsi_core)) deallocate(state%hpsi_core)
    allocate(state%psi_core(core_size(1), core_size(2), core_size(3), nstate, nspin))
    allocate(state%hpsi_core(core_size(1), core_size(2), core_size(3), nstate, nspin))
    state%psi_core = (0.0d0, 0.0d0)
    state%hpsi_core = (0.0d0, 0.0d0)
    state%fragment = fragment
    state%core_size = core_size
    state%halo_width = halo_width
    state%nstate = nstate
    state%nspin = nspin
    state%dg_ground_state_ready = .false.
    state%dg_ground_state_residual = huge(1.0d0)
    do axis = 1, 3
      do side = -1, 1, 2
        iface = nodal_face_slot(axis, side)
        state%faces(iface)%axis = axis
        state%faces(iface)%side = side
      end do
    end do
    state%initialized = .true.
  end subroutine allocate_nodal_state

  subroutine allocate_nodal_face_buffers(face, core_size, halo_width, nstate, nspin)
    type(nodal_face_info), intent(inout) :: face
    integer, intent(in) :: core_size(3), halo_width, nstate, nspin
    integer :: tangential(2)

    if (face%axis < 1 .or. face%axis > 3 .or. abs(face%side) /= 1) &
      stop 'nodal DG: invalid face identity'
    select case(face%axis)
    case(1)
      tangential = [core_size(2), core_size(3)]
    case(2)
      tangential = [core_size(1), core_size(3)]
    case(3)
      tangential = [core_size(1), core_size(2)]
    end select
    face%width = halo_width
    face%tangential_size = tangential
    allocate(face%send_value(halo_width,tangential(1),tangential(2),nstate,nspin))
    allocate(face%recv_value(halo_width,tangential(1),tangential(2),nstate,nspin))
    allocate(face%send_normal(halo_width,tangential(1),tangential(2),nstate,nspin))
    allocate(face%recv_normal(halo_width,tangential(1),tangential(2),nstate,nspin))
    face%send_value = (0.0d0, 0.0d0)
    face%recv_value = (0.0d0, 0.0d0)
    face%send_normal = (0.0d0, 0.0d0)
    face%recv_normal = (0.0d0, 0.0d0)
  end subroutine allocate_nodal_face_buffers

  subroutine accept_nodal_dg_ground_state(state, residual, tolerance)
    type(s_dg_nodal_state), intent(inout) :: state
    real(8), intent(in) :: residual, tolerance

    if (residual < 0.0d0 .or. residual /= residual) stop 'nodal DG: invalid ground-state residual'
    if (tolerance <= 0.0d0) stop 'nodal DG: invalid ground-state tolerance'
    state%dg_ground_state_residual = residual
    if (residual > tolerance) stop 'nodal DG: ground state is not stationary under the DG Hamiltonian'
    state%dg_ground_state_ready = .true.
  end subroutine accept_nodal_dg_ground_state

end module rt_dg_nodal_types
