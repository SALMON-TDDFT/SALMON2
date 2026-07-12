!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
!=======================================================================
! Periodic six-face topology for core-owned real-space DG states.
!=======================================================================
module rt_dg_nodal_halo
  use rt_dg_nodal_types, only: nodal_face_info, nodal_face_count, nodal_face_slot
  implicit none
  private
  public :: initialize_nodal_face_topology, neighbor_fragment_from_coords
  public :: refresh_nodal_single_fragment_halos

contains

  subroutine initialize_nodal_face_topology(faces, fragment_coords, num_fragment, fragment_rank, halo_width)
    type(nodal_face_info), intent(out) :: faces(nodal_face_count)
    integer, intent(in) :: fragment_coords(3)
    integer, intent(in) :: num_fragment(3)
    integer, intent(in) :: fragment_rank(:)
    integer, intent(in) :: halo_width
    integer :: axis, side, iface, neighbor_coords(3), jfrag

    if (any(num_fragment < 1)) stop 'nodal DG: invalid fragment topology'
    if (any(fragment_coords < 0) .or. any(fragment_coords >= num_fragment)) &
      stop 'nodal DG: fragment coordinate out of bounds'
    if (halo_width < 1) stop 'nodal DG: halo width must be positive'
    if (size(fragment_rank) < product(num_fragment)) stop 'nodal DG: incomplete fragment-rank map'

    do axis = 1, 3
      do side = -1, 1, 2
        iface = nodal_face_slot(axis, side)
        neighbor_coords(:) = fragment_coords(:)
        neighbor_coords(axis) = modulo(fragment_coords(axis) + side, num_fragment(axis))
        jfrag = neighbor_fragment_from_coords(neighbor_coords, num_fragment)

        faces(iface)%axis = axis
        faces(iface)%side = side
        faces(iface)%neighbor_fragment = jfrag
        faces(iface)%neighbor_rank = fragment_rank(jfrag)
        faces(iface)%width = halo_width

        if (faces(iface)%axis /= axis .or. faces(iface)%side /= side) &
          stop 'nodal DG: face identity was corrupted'
      end do
    end do
  end subroutine initialize_nodal_face_topology

  pure integer function neighbor_fragment_from_coords(coords, num_fragment) result(ifrag)
    integer, intent(in) :: coords(3), num_fragment(3)

    ifrag = 1 + coords(3) + num_fragment(3) * &
      (coords(2) + num_fragment(2) * coords(1))
  end function neighbor_fragment_from_coords

  subroutine refresh_nodal_single_fragment_halos(state)
    use rt_dg_nodal_types, only: s_dg_nodal_state
    type(s_dg_nodal_state), intent(inout) :: state
    integer :: axis, side, iface, layer, i1, i2, source_index

    if (.not. state%initialized) stop 'nodal DG: halo refresh before state initialization'
    do axis = 1, 3
      do side = -1, 1, 2
        iface = nodal_face_slot(axis, side)
        if (.not. allocated(state%faces(iface)%recv_value)) stop 'nodal DG: halo buffer is not allocated'
        do layer = 1, state%halo_width
          source_index = merge(state%core_size(axis)-layer+1, layer, side < 0)
          do i2 = 1, state%faces(iface)%tangential_size(2)
            do i1 = 1, state%faces(iface)%tangential_size(1)
              select case(axis)
              case(1)
                state%faces(iface)%recv_value(layer,i1,i2,:,:) = &
                  state%psi_core(source_index,i1,i2,:,:)
              case(2)
                state%faces(iface)%recv_value(layer,i1,i2,:,:) = &
                  state%psi_core(i1,source_index,i2,:,:)
              case(3)
                state%faces(iface)%recv_value(layer,i1,i2,:,:) = &
                  state%psi_core(i1,i2,source_index,:,:)
              end select
            end do
          end do
        end do
      end do
    end do
  end subroutine refresh_nodal_single_fragment_halos

end module rt_dg_nodal_halo
