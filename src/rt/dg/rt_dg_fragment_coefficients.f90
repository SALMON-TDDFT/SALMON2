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

module rt_dg_fragment_coefficients
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  implicit none

  private
  public :: rebuild_coef_owner_map
  public :: get_subgroup_block_owner_rank
  public :: get_fragment_coef_owner_rank

contains

  subroutine validate_coef_owner_map(dg_frag, stage_label)
    use communication, only: comm_get_max
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    character(*), intent(in) :: stage_label

    integer :: ispin, global_idx, probe_row, owner_rank
    integer :: owner_min, owner_max, invalid_local, invalid_global
    integer :: probe_count

    if (.not. allocated(dg_frag%coef_owner)) return
    if (size(dg_frag%coef_owner, 1) <= 0 .or. size(dg_frag%coef_owner, 2) <= 0) return

    do ispin = 1, min(dg_frag%nspin, size(dg_frag%coef_owner, 2))
      invalid_local = 0
      do global_idx = 1, size(dg_frag%coef_owner, 1)
        owner_rank = dg_frag%coef_owner(global_idx, ispin)
        if (owner_rank < -1 .or. owner_rank >= dg_frag%isize) invalid_local = invalid_local + 1
      end do
      invalid_global = invalid_local
      call comm_get_max(invalid_global, dg_frag%icomm)
      if (invalid_global > 0) then
        write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "[FATAL] invalid coef owner values after build: stage=", &
          trim(stage_label), " rank=", dg_frag%id, " ispin=", ispin, " invalid_local=", invalid_local
        stop "DG-Fragment RT: invalid coef owner values after build"
      end if

      probe_count = min(8, size(dg_frag%coef_owner, 1))
      do probe_row = 1, probe_count
        owner_rank = dg_frag%coef_owner(probe_row, ispin)
        owner_min = -owner_rank
        call comm_get_max(owner_min, dg_frag%icomm)
        owner_min = -owner_min
        owner_max = owner_rank
        call comm_get_max(owner_max, dg_frag%icomm)
        if (owner_min /= owner_max) then
          write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0,a,i0)') "[FATAL] inconsistent coef owner after build: stage=", &
            trim(stage_label), " rank=", dg_frag%id, " ispin=", ispin, " row=", probe_row, &
            " owner_min=", owner_min, " owner_max=", owner_max
          stop "DG-Fragment RT: inconsistent coef owner after build"
        end if
      end do
    end do
  end subroutine validate_coef_owner_map

  subroutine rebuild_local_coef_storage(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    integer :: ispin, global_idx, local_idx, old_local_idx
    integer :: local_coef_max, old_nstate, old_nspin
    integer, allocatable :: new_count(:)
    integer, allocatable :: new_global_ids(:,:)
    integer, allocatable :: new_global_to_local(:,:)
    complex(8), allocatable :: coef_new_rows(:,:,:)
    complex(8), allocatable :: coef_newbuf_rows(:,:,:)
    complex(8), allocatable :: coef_work_rows(:,:,:)
    logical :: has_coef, has_coef_new, has_coef_work

    if (.not. allocated(dg_frag%coef_owner)) return

    allocate(new_count(dg_frag%nspin))
    new_count(:) = 0
    do ispin = 1, dg_frag%nspin
      do global_idx = 1, dg_frag%n_mat_max
        if (dg_frag%coef_owner(global_idx, ispin) == dg_frag%id) then
          new_count(ispin) = new_count(ispin) + 1
        end if
      end do
    end do

    local_coef_max = max(1, maxval(new_count(1:dg_frag%nspin)))
    allocate(new_global_ids(local_coef_max, dg_frag%nspin))
    allocate(new_global_to_local(dg_frag%n_mat_max, dg_frag%nspin))
    new_global_ids(:, :) = 0
    new_global_to_local(:, :) = 0

    do ispin = 1, dg_frag%nspin
      local_idx = 0
      do global_idx = 1, dg_frag%n_mat_max
        if (dg_frag%coef_owner(global_idx, ispin) /= dg_frag%id) cycle
        local_idx = local_idx + 1
        new_global_ids(local_idx, ispin) = global_idx
        new_global_to_local(global_idx, ispin) = local_idx
      end do
    end do

    has_coef = allocated(dg_frag%coef)
    has_coef_new = allocated(dg_frag%coef_new)
    has_coef_work = allocated(dg_frag%coef_work)
    if (has_coef) then
      old_nstate = size(dg_frag%coef, 2)
      old_nspin = min(dg_frag%nspin, size(dg_frag%coef, 3))
      allocate(coef_new_rows(local_coef_max, old_nstate, size(dg_frag%coef, 3)))
      coef_new_rows(:, :, :) = (0.0d0, 0.0d0)
      do ispin = 1, old_nspin
        do local_idx = 1, new_count(ispin)
          global_idx = new_global_ids(local_idx, ispin)
          old_local_idx = 0
          if (allocated(dg_frag%coef_global_to_local)) then
            old_local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
          end if
          if (old_local_idx < 1 .or. old_local_idx > size(dg_frag%coef, 1)) cycle
          coef_new_rows(local_idx, 1:old_nstate, ispin) = &
            dg_frag%coef(old_local_idx, 1:old_nstate, ispin)
        end do
      end do
    end if

    if (has_coef_new) then
      old_nstate = size(dg_frag%coef_new, 2)
      old_nspin = min(dg_frag%nspin, size(dg_frag%coef_new, 3))
      allocate(coef_newbuf_rows(local_coef_max, old_nstate, size(dg_frag%coef_new, 3)))
      coef_newbuf_rows(:, :, :) = (0.0d0, 0.0d0)
      do ispin = 1, old_nspin
        do local_idx = 1, new_count(ispin)
          global_idx = new_global_ids(local_idx, ispin)
          old_local_idx = 0
          if (allocated(dg_frag%coef_global_to_local)) then
            old_local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
          end if
          if (old_local_idx < 1 .or. old_local_idx > size(dg_frag%coef_new, 1)) cycle
          coef_newbuf_rows(local_idx, 1:old_nstate, ispin) = &
            dg_frag%coef_new(old_local_idx, 1:old_nstate, ispin)
        end do
      end do
    end if

    if (has_coef_work) then
      old_nstate = size(dg_frag%coef_work, 2)
      old_nspin = min(dg_frag%nspin, size(dg_frag%coef_work, 3))
      allocate(coef_work_rows(local_coef_max, old_nstate, size(dg_frag%coef_work, 3)))
      coef_work_rows(:, :, :) = (0.0d0, 0.0d0)
      do ispin = 1, old_nspin
        do local_idx = 1, new_count(ispin)
          global_idx = new_global_ids(local_idx, ispin)
          old_local_idx = 0
          if (allocated(dg_frag%coef_global_to_local)) then
            old_local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
          end if
          if (old_local_idx < 1 .or. old_local_idx > size(dg_frag%coef_work, 1)) cycle
          coef_work_rows(local_idx, 1:old_nstate, ispin) = &
            dg_frag%coef_work(old_local_idx, 1:old_nstate, ispin)
        end do
      end do
    end if

    if (allocated(dg_frag%local_coef_count)) deallocate(dg_frag%local_coef_count)
    if (allocated(dg_frag%local_coef_global_ids)) deallocate(dg_frag%local_coef_global_ids)
    if (allocated(dg_frag%coef_global_to_local)) deallocate(dg_frag%coef_global_to_local)
    call move_alloc(new_count, dg_frag%local_coef_count)
    call move_alloc(new_global_ids, dg_frag%local_coef_global_ids)
    call move_alloc(new_global_to_local, dg_frag%coef_global_to_local)

    if (has_coef) then
      deallocate(dg_frag%coef)
      call move_alloc(coef_new_rows, dg_frag%coef)
    end if
    if (has_coef_new) then
      deallocate(dg_frag%coef_new)
      call move_alloc(coef_newbuf_rows, dg_frag%coef_new)
    end if
    if (has_coef_work) then
      deallocate(dg_frag%coef_work)
      call move_alloc(coef_work_rows, dg_frag%coef_work)
    end if
  end subroutine rebuild_local_coef_storage

  subroutine rebuild_coef_owner_map(dg_frag, stage_label)
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    character(*), intent(in) :: stage_label

    integer :: ispin, ifrag, io, global_idx, nbasis_iter, owner_new
    integer :: env_status, env_len
    integer :: conflict_local
    character(len=64) :: env_owner_io_split
    logical :: owner_root_only
    logical, save :: printed_owner_mode = .false.

    owner_root_only = .true.
    env_owner_io_split = ""
    env_status = 1
    env_len = 0
    call get_environment_variable("SALMON_DG_COEF_OWNER_IO_SPLIT", env_owner_io_split, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      select case(trim(adjustl(env_owner_io_split(1:env_len))))
      case('y','Y','yes','YES','true','TRUE','1')
        owner_root_only = .false.
      end select
    end if

    if (dg_frag%parallel_mode_orbital) then
      owner_root_only = .false.
    end if

    if (dg_frag%id == 0 .and. .not. printed_owner_mode) then
      if (dg_frag%parallel_mode_orbital) then
        write(*,'(1x,a)') "[INFO] DG coef owner mode: orbital row-split"
      else if (owner_root_only) then
        write(*,'(1x,a)') "[INFO] DG coef owner mode: root-only (io split disabled)"
      else
        write(*,'(1x,a)') "[INFO] DG coef owner mode: io-split (legacy)"
      end if
      flush(6)
      printed_owner_mode = .true.
    end if

    if (.not. allocated(dg_frag%id_array)) then
      write(*,'(1x,a,a,a,i0)') "[FATAL] id_array is not allocated before coef owner build: stage=", &
        trim(stage_label), " rank=", dg_frag%id
      stop "DG-Fragment RT: missing id_array before coef owner build"
    end if

    if (allocated(dg_frag%coef_owner)) deallocate(dg_frag%coef_owner)
    allocate(dg_frag%coef_owner(dg_frag%n_mat_max, dg_frag%nspin))
    dg_frag%coef_owner(:, :) = -1
    dg_frag%H_local_block_ids_valid = .false.
    conflict_local = 0
    do ispin = 1, dg_frag%nspin
      do ifrag = 1, dg_frag%n_frag
        nbasis_iter = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1))
        do io = 1, nbasis_iter
          global_idx = dg_frag%index_basis(io, ifrag, ispin)
          if (global_idx < 1 .or. global_idx > dg_frag%n_mat_max) cycle
          if (owner_root_only) then
            owner_new = dg_frag%id_array(ifrag)
          else
            owner_new = get_subgroup_block_owner_rank(dg_frag%id_array(ifrag), dg_frag%isize_frag, io, nbasis_iter)
          end if
          if (dg_frag%coef_owner(global_idx, ispin) == -1) then
            dg_frag%coef_owner(global_idx, ispin) = owner_new
          else if (dg_frag%coef_owner(global_idx, ispin) /= owner_new) then
            conflict_local = conflict_local + 1
          end if
        end do
      end do
    end do

    if (dg_frag%id == 0 .and. conflict_local > 0) then
      write(*,'(1x,a,a,a,i0)') "[WARN] coef_owner conflicts detected (kept first owner): stage=", &
        trim(stage_label), " conflicts=", conflict_local
    end if

    dg_frag%owned_coef_start = 0
    dg_frag%owned_coef_end = -1
    do global_idx = 1, dg_frag%n_mat_max
      if (any(dg_frag%coef_owner(global_idx, 1:dg_frag%nspin) == dg_frag%id)) then
        if (dg_frag%owned_coef_start == 0) dg_frag%owned_coef_start = global_idx
        dg_frag%owned_coef_end = global_idx
      end if
    end do

    call validate_coef_owner_map(dg_frag, stage_label)
    call rebuild_local_coef_storage(dg_frag)
  end subroutine rebuild_coef_owner_map

  integer function get_subgroup_block_owner_rank(root_rank, subgroup_size, local_row, n_local_rows) result(owner_rank)
    implicit none
    integer, intent(in) :: root_rank, subgroup_size, local_row, n_local_rows
    integer :: subgroup_root_rank, owner_offset
    integer :: block_base, block_rem, cutoff

    if (subgroup_size <= 1 .or. n_local_rows <= 0) then
      owner_rank = max(0, root_rank)
      return
    end if

    subgroup_root_rank = root_rank - mod(max(0, root_rank), subgroup_size)
    block_base = n_local_rows / subgroup_size
    block_rem = mod(n_local_rows, subgroup_size)
    cutoff = (block_base + 1) * block_rem

    if (local_row <= 0) then
      owner_offset = 0
    else if (block_base <= 0) then
      owner_offset = min(local_row - 1, subgroup_size - 1)
    else if (local_row <= cutoff) then
      owner_offset = (local_row - 1) / (block_base + 1)
    else
      owner_offset = block_rem + (local_row - cutoff - 1) / block_base
    end if

    owner_rank = subgroup_root_rank + min(owner_offset, subgroup_size - 1)
  end function get_subgroup_block_owner_rank

  integer function get_fragment_coef_owner_rank(dg_frag, ifrag, local_row, n_local_rows) result(owner_rank)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, local_row, n_local_rows

    if (ifrag < 1 .or. ifrag > size(dg_frag%id_array)) then
      owner_rank = 0
      return
    end if

    if (dg_frag%parallel_mode_orbital) then
      owner_rank = get_subgroup_block_owner_rank(dg_frag%id_array(ifrag), dg_frag%isize_frag, &
                                                 local_row, n_local_rows)
    else
      owner_rank = dg_frag%id_array(ifrag)
    end if
  end function get_fragment_coef_owner_rank

end module rt_dg_fragment_coefficients
