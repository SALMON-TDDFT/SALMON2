!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!

module rt_dg_fragment_layout
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  implicit none

  private
  public :: wrap_global_grid_index
  public :: get_fragment_grid_sender_rank
  public :: wrap_fragment_cartesian_index
  public :: cartesian_index_to_fragment
  public :: find_density_grid_owner
  public :: get_fragment_group_root_rank
  public :: fragment_from_rank_address
  public :: build_density_grid_owner_maps
  public :: build_fragment_global_boxes

contains

  subroutine build_density_grid_owner_maps(dg_frag)
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    integer :: ifrag, i_local, ix, iy, iz, ixg, iyg, izg, owner_rank, source_rank, source_root_rank, i
    integer :: subgroup_target_rank
    integer :: ifrag_count
    integer :: nx_max, ny_max, nz_max
    integer, allocatable :: recv_count(:), recv_cursor(:)
    integer, allocatable :: id_tmp(:)

    if (allocated(dg_frag%density_owner_map)) deallocate(dg_frag%density_owner_map)
    if (allocated(dg_frag%density_primary_local_map)) deallocate(dg_frag%density_primary_local_map)
    if (allocated(dg_frag%density_ixg_map)) deallocate(dg_frag%density_ixg_map)
    if (allocated(dg_frag%density_iyg_map)) deallocate(dg_frag%density_iyg_map)
    if (allocated(dg_frag%density_izg_map)) deallocate(dg_frag%density_izg_map)
    if (allocated(dg_frag%density_send_count)) deallocate(dg_frag%density_send_count)
    if (allocated(dg_frag%density_send_slot_map)) deallocate(dg_frag%density_send_slot_map)
    if (allocated(dg_frag%density_subgroup_send_count)) deallocate(dg_frag%density_subgroup_send_count)
    if (allocated(dg_frag%density_subgroup_send_slot_map)) deallocate(dg_frag%density_subgroup_send_slot_map)
    if (allocated(dg_frag%density_grid_points)) deallocate(dg_frag%density_grid_points)
    if (allocated(dg_frag%density_grid_point_count)) deallocate(dg_frag%density_grid_point_count)
    if (allocated(dg_frag%density_grid_bx)) deallocate(dg_frag%density_grid_bx)
    if (allocated(dg_frag%density_grid_by)) deallocate(dg_frag%density_grid_by)
    if (allocated(dg_frag%density_grid_bz)) deallocate(dg_frag%density_grid_bz)
    dg_frag%density_rhobf_box_cache_valid = .false.
    if (allocated(dg_frag%density_subgroup_self_ixg)) deallocate(dg_frag%density_subgroup_self_ixg)
    if (allocated(dg_frag%density_subgroup_self_iyg)) deallocate(dg_frag%density_subgroup_self_iyg)
    if (allocated(dg_frag%density_subgroup_self_izg)) deallocate(dg_frag%density_subgroup_self_izg)
    if (allocated(dg_frag%density_block_nblocks)) deallocate(dg_frag%density_block_nblocks)
    if (allocated(dg_frag%density_block_first_offset)) deallocate(dg_frag%density_block_first_offset)
    if (allocated(dg_frag%density_block_step)) deallocate(dg_frag%density_block_step)
    if (allocated(dg_frag%current_valid_grid_count)) deallocate(dg_frag%current_valid_grid_count)
    if (allocated(dg_frag%current_valid_ix)) deallocate(dg_frag%current_valid_ix)
    if (allocated(dg_frag%current_valid_iy)) deallocate(dg_frag%current_valid_iy)
    if (allocated(dg_frag%current_valid_iz)) deallocate(dg_frag%current_valid_iz)
    if (allocated(dg_frag%current_valid_ixg)) deallocate(dg_frag%current_valid_ixg)
    if (allocated(dg_frag%current_valid_iyg)) deallocate(dg_frag%current_valid_iyg)
    if (allocated(dg_frag%current_valid_izg)) deallocate(dg_frag%current_valid_izg)
    if (allocated(dg_frag%runtime_neighbor_pair_cache)) deallocate(dg_frag%runtime_neighbor_pair_cache)
    if (allocated(dg_frag%momentum_neighbor_pair_cache)) deallocate(dg_frag%momentum_neighbor_pair_cache)
    if (allocated(dg_frag%density_phi_block_cache)) deallocate(dg_frag%density_phi_block_cache)
    if (allocated(dg_frag%density_phi_block_count)) deallocate(dg_frag%density_phi_block_count)
    dg_frag%density_phi_block_size = 0
    dg_frag%density_phi_block_cache_valid = .false.
    if (allocated(dg_frag%density_phase_block_cache)) deallocate(dg_frag%density_phase_block_cache)
    dg_frag%density_phase_block_size = 0
    dg_frag%density_phase_block_npw = 0
    dg_frag%density_phase_block_cache_valid = .false.
    if (allocated(dg_frag%density_recv_map)) then
      do i = lbound(dg_frag%density_recv_map, 1), ubound(dg_frag%density_recv_map, 1)
        if (allocated(dg_frag%density_recv_map(i)%ixg)) deallocate(dg_frag%density_recv_map(i)%ixg)
        if (allocated(dg_frag%density_recv_map(i)%iyg)) deallocate(dg_frag%density_recv_map(i)%iyg)
        if (allocated(dg_frag%density_recv_map(i)%izg)) deallocate(dg_frag%density_recv_map(i)%izg)
        if (allocated(dg_frag%density_recv_map(i)%bx)) deallocate(dg_frag%density_recv_map(i)%bx)
        if (allocated(dg_frag%density_recv_map(i)%by)) deallocate(dg_frag%density_recv_map(i)%by)
        if (allocated(dg_frag%density_recv_map(i)%bz)) deallocate(dg_frag%density_recv_map(i)%bz)
      end do
      deallocate(dg_frag%density_recv_map)
    end if
    if (allocated(dg_frag%density_send_count)) deallocate(dg_frag%density_send_count)
    if (allocated(dg_frag%density_send_slot_map)) deallocate(dg_frag%density_send_slot_map)
    if (allocated(dg_frag%density_subgroup_send_count)) deallocate(dg_frag%density_subgroup_send_count)
    if (allocated(dg_frag%density_subgroup_send_slot_map)) deallocate(dg_frag%density_subgroup_send_slot_map)
    if (allocated(dg_frag%density_grid_points)) deallocate(dg_frag%density_grid_points)
    if (allocated(dg_frag%density_grid_point_count)) deallocate(dg_frag%density_grid_point_count)
    if (allocated(dg_frag%density_grid_bx)) deallocate(dg_frag%density_grid_bx)
    if (allocated(dg_frag%density_grid_by)) deallocate(dg_frag%density_grid_by)
    if (allocated(dg_frag%density_grid_bz)) deallocate(dg_frag%density_grid_bz)
    dg_frag%density_rhobf_box_cache_valid = .false.
    if (allocated(dg_frag%density_subgroup_self_ixg)) deallocate(dg_frag%density_subgroup_self_ixg)
    if (allocated(dg_frag%density_subgroup_self_iyg)) deallocate(dg_frag%density_subgroup_self_iyg)
    if (allocated(dg_frag%density_subgroup_self_izg)) deallocate(dg_frag%density_subgroup_self_izg)
    if (allocated(dg_frag%density_block_nblocks)) deallocate(dg_frag%density_block_nblocks)
    if (allocated(dg_frag%density_block_first_offset)) deallocate(dg_frag%density_block_first_offset)
    if (allocated(dg_frag%density_block_step)) deallocate(dg_frag%density_block_step)
    if (allocated(dg_frag%current_valid_grid_count)) deallocate(dg_frag%current_valid_grid_count)
    if (allocated(dg_frag%current_valid_ix)) deallocate(dg_frag%current_valid_ix)
    if (allocated(dg_frag%current_valid_iy)) deallocate(dg_frag%current_valid_iy)
    if (allocated(dg_frag%current_valid_iz)) deallocate(dg_frag%current_valid_iz)
    if (allocated(dg_frag%current_valid_ixg)) deallocate(dg_frag%current_valid_ixg)
    if (allocated(dg_frag%current_valid_iyg)) deallocate(dg_frag%current_valid_iyg)
    if (allocated(dg_frag%current_valid_izg)) deallocate(dg_frag%current_valid_izg)
    if (allocated(dg_frag%runtime_neighbor_pair_cache)) deallocate(dg_frag%runtime_neighbor_pair_cache)
    if (allocated(dg_frag%momentum_neighbor_pair_cache)) deallocate(dg_frag%momentum_neighbor_pair_cache)
    if (allocated(dg_frag%density_recv_map)) then
      do ifrag = lbound(dg_frag%density_recv_map, 1), ubound(dg_frag%density_recv_map, 1)
        if (allocated(dg_frag%density_recv_map(ifrag)%ixg)) deallocate(dg_frag%density_recv_map(ifrag)%ixg)
        if (allocated(dg_frag%density_recv_map(ifrag)%iyg)) deallocate(dg_frag%density_recv_map(ifrag)%iyg)
        if (allocated(dg_frag%density_recv_map(ifrag)%izg)) deallocate(dg_frag%density_recv_map(ifrag)%izg)
        if (allocated(dg_frag%density_recv_map(ifrag)%bx)) deallocate(dg_frag%density_recv_map(ifrag)%bx)
        if (allocated(dg_frag%density_recv_map(ifrag)%by)) deallocate(dg_frag%density_recv_map(ifrag)%by)
        if (allocated(dg_frag%density_recv_map(ifrag)%bz)) deallocate(dg_frag%density_recv_map(ifrag)%bz)
      end do
      deallocate(dg_frag%density_recv_map)
    end if
    if (.not. associated(dg_frag%mg)) return
    if (.not. allocated(dg_frag%nxyz_domain)) return
    if (.not. allocated(dg_frag%ixyz_frag)) return
    if (dg_frag%ifrag_end < dg_frag%ifrag_start) return

    call build_fragment_global_boxes(dg_frag)

    if (allocated(dg_frag%id_array)) then
      allocate(id_tmp(dg_frag%n_frag))
      id_tmp = 0
      if (dg_frag%is_frag_root) then
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          id_tmp(ifrag) = dg_frag%id + 1
        end do
      end if
      call comm_summation(id_tmp, dg_frag%id_array, dg_frag%n_frag, dg_frag%icomm)
      dg_frag%id_array = dg_frag%id_array - 1
      deallocate(id_tmp)
    end if

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    nx_max = max(1, maxval(dg_frag%nxyz_domain(1, dg_frag%ifrag_start:dg_frag%ifrag_end)))
    ny_max = max(1, maxval(dg_frag%nxyz_domain(2, dg_frag%ifrag_start:dg_frag%ifrag_end)))
    nz_max = max(1, maxval(dg_frag%nxyz_domain(3, dg_frag%ifrag_start:dg_frag%ifrag_end)))

    allocate(dg_frag%density_owner_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_primary_local_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_ixg_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_iyg_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_izg_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_send_count(0:dg_frag%isize-1))
    allocate(dg_frag%density_send_slot_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_subgroup_send_count(0:dg_frag%isize_frag-1))
    allocate(dg_frag%density_subgroup_send_slot_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_recv_map(0:dg_frag%isize-1))
    allocate(dg_frag%density_grid_points(nx_max * ny_max * nz_max, ifrag_count))
    allocate(dg_frag%density_grid_point_count(ifrag_count))
    dg_frag%density_owner_map = dg_frag%id
    dg_frag%density_primary_local_map = .false.
    dg_frag%density_ixg_map = 1
    dg_frag%density_iyg_map = 1
    dg_frag%density_izg_map = 1
    dg_frag%density_send_count = 0
    dg_frag%density_send_slot_map = 0
    dg_frag%density_subgroup_send_count = 0
    dg_frag%density_subgroup_send_slot_map = 0
    dg_frag%density_grid_point_count = 0

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      do iz = 1, dg_frag%nxyz_domain(3, ifrag)
        izg = wrap_global_grid_index(dg_frag%frag_core_lo(3, ifrag) + iz - 1, dg_frag%lgnum_total(3))
        do iy = 1, dg_frag%nxyz_domain(2, ifrag)
          iyg = wrap_global_grid_index(dg_frag%frag_core_lo(2, ifrag) + iy - 1, dg_frag%lgnum_total(2))
          do ix = 1, dg_frag%nxyz_domain(1, ifrag)
            ixg = wrap_global_grid_index(dg_frag%frag_core_lo(1, ifrag) + ix - 1, dg_frag%lgnum_total(1))
            dg_frag%density_ixg_map(ix, iy, iz, i_local) = ixg
            dg_frag%density_iyg_map(ix, iy, iz, i_local) = iyg
            dg_frag%density_izg_map(ix, iy, iz, i_local) = izg
            dg_frag%density_primary_local_map(ix, iy, iz, i_local) = .true.
            source_rank = dg_frag%id_array(ifrag)
            subgroup_target_rank = source_rank - dg_frag%id_array(ifrag)
            if (dg_frag%density_primary_local_map(ix, iy, iz, i_local)) then
              source_rank = get_fragment_grid_sender_rank(dg_frag%id_array(ifrag), dg_frag%nxyz_domain(:, ifrag), &
                                                           ix, iy, iz, dg_frag%parallel_mode_orbital)
              subgroup_target_rank = source_rank - dg_frag%id_array(ifrag)
              if (subgroup_target_rank < 0 .or. subgroup_target_rank > dg_frag%isize_frag - 1) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                     "[FATAL] density subgroup target out of range: rank=", dg_frag%id, &
                     " ifrag=", ifrag, " ix=", ix, " iy=", iy, " iz=", iz
                flush(6)
                stop "DG-Fragment RT: density subgroup target out of range"
              end if
              if (subgroup_target_rank /= dg_frag%id_frag) then
                dg_frag%density_subgroup_send_count(subgroup_target_rank) = &
                  dg_frag%density_subgroup_send_count(subgroup_target_rank) + 1
                dg_frag%density_subgroup_send_slot_map(ix, iy, iz, i_local) = &
                  dg_frag%density_subgroup_send_count(subgroup_target_rank)
              end if
            end if
            owner_rank = find_density_grid_owner(dg_frag, ixg, iyg, izg, dg_frag%id_array(ifrag))
            dg_frag%density_owner_map(ix, iy, iz, i_local) = owner_rank
            source_rank = dg_frag%id_array(ifrag)
            if (dg_frag%density_primary_local_map(ix, iy, iz, i_local) .and. owner_rank /= source_rank .and. &
                dg_frag%is_frag_root) then
              dg_frag%density_send_count(owner_rank) = dg_frag%density_send_count(owner_rank) + 1
              dg_frag%density_send_slot_map(ix, iy, iz, i_local) = dg_frag%density_send_count(owner_rank)
            end if
            dg_frag%density_grid_point_count(i_local) = dg_frag%density_grid_point_count(i_local) + 1
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%ix = ix
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%iy = iy
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%iz = iz
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%ixg = ixg
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%iyg = iyg
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%izg = izg
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%owner_rank = owner_rank
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%is_primary = &
              dg_frag%density_primary_local_map(ix, iy, iz, i_local)
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%send_slot = &
              dg_frag%density_send_slot_map(ix, iy, iz, i_local)
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%subgroup_send_slot = &
              dg_frag%density_subgroup_send_slot_map(ix, iy, iz, i_local)
          end do
        end do
      end do
    end do

    allocate(recv_count(0:dg_frag%isize-1), recv_cursor(0:dg_frag%isize-1))
    recv_count = 0
    do ifrag = 1, dg_frag%n_frag
      source_root_rank = dg_frag%id_array(ifrag)
      do iz = 1, dg_frag%nxyz_domain(3, ifrag)
        izg = wrap_global_grid_index(dg_frag%frag_core_lo(3, ifrag) + iz - 1, dg_frag%lgnum_total(3))
        do iy = 1, dg_frag%nxyz_domain(2, ifrag)
          iyg = wrap_global_grid_index(dg_frag%frag_core_lo(2, ifrag) + iy - 1, dg_frag%lgnum_total(2))
          do ix = 1, dg_frag%nxyz_domain(1, ifrag)
            ixg = wrap_global_grid_index(dg_frag%frag_core_lo(1, ifrag) + ix - 1, dg_frag%lgnum_total(1))
            source_rank = source_root_rank
            owner_rank = find_density_grid_owner(dg_frag, ixg, iyg, izg, source_root_rank)
            if (source_rank == dg_frag%id) cycle
            if (owner_rank == dg_frag%id) recv_count(source_rank) = recv_count(source_rank) + 1
          end do
        end do
      end do
    end do

    do source_rank = 0, dg_frag%isize - 1
      dg_frag%density_recv_map(source_rank)%npts = recv_count(source_rank)
      if (recv_count(source_rank) <= 0) cycle
      allocate(dg_frag%density_recv_map(source_rank)%ixg(recv_count(source_rank)))
      allocate(dg_frag%density_recv_map(source_rank)%iyg(recv_count(source_rank)))
      allocate(dg_frag%density_recv_map(source_rank)%izg(recv_count(source_rank)))
    end do

    recv_cursor = 0
    do ifrag = 1, dg_frag%n_frag
      source_root_rank = dg_frag%id_array(ifrag)
      do iz = 1, dg_frag%nxyz_domain(3, ifrag)
        izg = wrap_global_grid_index(dg_frag%frag_core_lo(3, ifrag) + iz - 1, dg_frag%lgnum_total(3))
        do iy = 1, dg_frag%nxyz_domain(2, ifrag)
          iyg = wrap_global_grid_index(dg_frag%frag_core_lo(2, ifrag) + iy - 1, dg_frag%lgnum_total(2))
          do ix = 1, dg_frag%nxyz_domain(1, ifrag)
            ixg = wrap_global_grid_index(dg_frag%frag_core_lo(1, ifrag) + ix - 1, dg_frag%lgnum_total(1))
            source_rank = source_root_rank
            if (source_rank == dg_frag%id) cycle
            owner_rank = find_density_grid_owner(dg_frag, ixg, iyg, izg, source_root_rank)
            if (owner_rank == dg_frag%id) then
              recv_cursor(source_rank) = recv_cursor(source_rank) + 1
              dg_frag%density_recv_map(source_rank)%ixg(recv_cursor(source_rank)) = ixg
              dg_frag%density_recv_map(source_rank)%iyg(recv_cursor(source_rank)) = iyg
              dg_frag%density_recv_map(source_rank)%izg(recv_cursor(source_rank)) = izg
            end if
          end do
        end do
      end do
    end do

    deallocate(recv_count, recv_cursor)
  end subroutine build_density_grid_owner_maps

  subroutine build_fragment_global_boxes(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: ifrag

    if (.not. allocated(dg_frag%nxyz_domain) .or. .not. allocated(dg_frag%ixyz_frag)) return

    if (allocated(dg_frag%frag_core_lo)) deallocate(dg_frag%frag_core_lo)
    if (allocated(dg_frag%frag_core_hi)) deallocate(dg_frag%frag_core_hi)
    if (allocated(dg_frag%frag_buf_lo)) deallocate(dg_frag%frag_buf_lo)
    if (allocated(dg_frag%frag_buf_hi)) deallocate(dg_frag%frag_buf_hi)
    allocate(dg_frag%frag_core_lo(3, dg_frag%n_frag), dg_frag%frag_core_hi(3, dg_frag%n_frag))
    allocate(dg_frag%frag_buf_lo(3, dg_frag%n_frag), dg_frag%frag_buf_hi(3, dg_frag%n_frag))

    do ifrag = 1, dg_frag%n_frag
      dg_frag%frag_core_lo(:, ifrag) = dg_frag%ixyz_frag(:, ifrag)
      dg_frag%frag_core_hi(:, ifrag) = dg_frag%ixyz_frag(:, ifrag) + dg_frag%nxyz_domain(:, ifrag) - 1
      dg_frag%frag_buf_lo(:, ifrag) = dg_frag%frag_core_lo(:, ifrag) - dg_frag%nxyz_buffer(:)
      dg_frag%frag_buf_hi(:, ifrag) = dg_frag%frag_core_hi(:, ifrag) + dg_frag%nxyz_buffer(:)
    end do

    if (dg_frag%parallel_mode_orbital .and. dg_frag%ifrag_group >= 1 .and. dg_frag%ifrag_group <= dg_frag%n_frag) then
      dg_frag%rank_core_lo(:) = dg_frag%frag_core_lo(:, dg_frag%ifrag_group)
      dg_frag%rank_core_hi(:) = dg_frag%frag_core_hi(:, dg_frag%ifrag_group)
      dg_frag%rank_buf_lo(:) = dg_frag%frag_buf_lo(:, dg_frag%ifrag_group)
      dg_frag%rank_buf_hi(:) = dg_frag%frag_buf_hi(:, dg_frag%ifrag_group)
    else
      dg_frag%rank_core_lo(:) = dg_frag%mg%is(:)
      dg_frag%rank_core_hi(:) = dg_frag%mg%ie(:)
      dg_frag%rank_buf_lo(:) = dg_frag%mg%is(:) - dg_frag%nxyz_buffer(:)
      dg_frag%rank_buf_hi(:) = dg_frag%mg%ie(:) + dg_frag%nxyz_buffer(:)
    end if
  end subroutine build_fragment_global_boxes

  integer function wrap_global_grid_index(ig_raw, ngrid) result(ig)
    implicit none
    integer, intent(in) :: ig_raw, ngrid

    ig = modulo(ig_raw - 1, ngrid) + 1
  end function wrap_global_grid_index

  integer function get_fragment_grid_sender_rank(root_rank, ndom, ix, iy, iz, parallel_mode_orbital) result(sender_rank)
    use salmon_global, only: nproc_rgrid
    implicit none
    integer, intent(in) :: root_rank, ndom(3), ix, iy, iz
    logical, intent(in) :: parallel_mode_orbital
    integer :: ipx, ipy, ipz, coords(3), nsize

    if (parallel_mode_orbital) then
      sender_rank = root_rank
      return
    end if

    ipx = max(1, nproc_rgrid(1))
    ipy = max(1, nproc_rgrid(2))
    ipz = max(1, nproc_rgrid(3))

    nsize = max(1, (ndom(1) + ipx - 1) / ipx)
    coords(1) = min(ipx - 1, max(0, (ix - 1) / nsize))
    nsize = max(1, (ndom(2) + ipy - 1) / ipy)
    coords(2) = min(ipy - 1, max(0, (iy - 1) / nsize))
    nsize = max(1, (ndom(3) + ipz - 1) / ipz)
    coords(3) = min(ipz - 1, max(0, (iz - 1) / nsize))

    sender_rank = root_rank + coords(1) + ipx * (coords(2) + ipy * coords(3))
  end function get_fragment_grid_sender_rank

  integer function wrap_fragment_cartesian_index(i, ndiv) result(iwrap)
    implicit none
    integer, intent(in) :: i, ndiv

    iwrap = modulo(i - 1, max(1, ndiv)) + 1
  end function wrap_fragment_cartesian_index

  integer function cartesian_index_to_fragment(dg_frag, idx) result(ifrag)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: idx(3)

    ifrag = ((idx(1) - 1) * dg_frag%num_fragment(2) + (idx(2) - 1)) * &
            dg_frag%num_fragment(3) + idx(3)
  end function cartesian_index_to_fragment

  integer function find_density_grid_owner(dg_frag, ixg, iyg, izg, hint_rank) result(owner)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ixg, iyg, izg
    integer, intent(in), optional :: hint_rank
    integer :: jrank, first_match, nfrag_ranks, hint_group, rank_group
    integer :: xlo, xhi, ylo, yhi, zlo, zhi

    owner = -1
    first_match = -1
    nfrag_ranks = max(1, dg_frag%isize_frag)

    hint_group = -1
    if (present(hint_rank) .and. .not. dg_frag%parallel_mode_orbital) hint_group = max(0, hint_rank) / nfrag_ranks
    do jrank = 0, dg_frag%isize - 1
      xlo = dg_frag%mg%is_all(1, jrank)
      xhi = dg_frag%mg%ie_all(1, jrank)
      if (ixg < xlo .or. ixg > xhi) cycle
      ylo = dg_frag%mg%is_all(2, jrank)
      yhi = dg_frag%mg%ie_all(2, jrank)
      if (iyg < ylo .or. iyg > yhi) cycle
      zlo = dg_frag%mg%is_all(3, jrank)
      zhi = dg_frag%mg%ie_all(3, jrank)
      if (izg < zlo .or. izg > zhi) cycle
      if (first_match < 0) first_match = jrank
      if (hint_group >= 0) then
        rank_group = jrank / nfrag_ranks
        if (rank_group == hint_group) then
          owner = jrank
          return
        end if
      end if
    end do
    if (first_match >= 0) then
      owner = first_match
    else
      owner = dg_frag%id
    end if
  end function find_density_grid_owner

  integer function get_fragment_group_root_rank(ifrag, nproc_frag) result(owner_rank)
    use salmon_global, only: num_fragment
    implicit none
    integer, intent(in) :: ifrag, nproc_frag
    integer :: ix_frag, iy_frag, iz_frag, rem, group_index

    if (ifrag < 1 .or. nproc_frag <= 0) then
      owner_rank = 0
    else
      ! DC fragment files are numbered with z as the fastest index:
      !   ifrag = ((ix-1) * ny + (iy-1)) * nz + iz.
      ! The orbital_sequential MPI layout numbers real-space rank groups with
      ! x as the fastest index.  Convert explicitly so each fragment subgroup
      ! reads the DC file matching its parent real-space box.
      ix_frag = (ifrag - 1) / max(1, num_fragment(2) * num_fragment(3)) + 1
      rem = modulo(ifrag - 1, max(1, num_fragment(2) * num_fragment(3)))
      iy_frag = rem / max(1, num_fragment(3)) + 1
      iz_frag = modulo(rem, max(1, num_fragment(3))) + 1
      group_index = ((iz_frag - 1) * max(1, num_fragment(2)) + (iy_frag - 1)) * &
                    max(1, num_fragment(1)) + (ix_frag - 1)
      owner_rank = group_index * nproc_frag
    end if
  end function get_fragment_group_root_rank

  integer function fragment_from_rank_address(iaddr, nfrag_axis) result(ifrag)
    implicit none
    integer, intent(in) :: iaddr(3), nfrag_axis(3)
    integer :: ix_frag, iy_frag, iz_frag

    ix_frag = iaddr(1) + 1
    iy_frag = iaddr(2) + 1
    iz_frag = iaddr(3) + 1
    ifrag = ((ix_frag - 1) * max(1, nfrag_axis(2)) + (iy_frag - 1)) * &
            max(1, nfrag_axis(3)) + iz_frag
  end function fragment_from_rank_address

end module rt_dg_fragment_layout
