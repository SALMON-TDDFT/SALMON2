!=======================================================================
  ! Calculate Hamiltonian matrix in fragment basis
  !=======================================================================
  logical function is_momentum_neighbor_axis(lg, s1, n1, s2, n2) result(ok)
    implicit none
    integer, intent(in) :: lg, s1, n1, s2, n2
    integer :: e1, e2, s1_next, s2_next

    e1 = s1 + n1 - 1
    e2 = s2 + n2 - 1
    s1_next = modulo(e1, lg) + 1
    s2_next = modulo(e2, lg) + 1
    ok = ((s1 == s2) .and. (n1 == n2)) .or. (s1 == s2_next) .or. (s2 == s1_next)
  end function is_momentum_neighbor_axis

  logical function is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col) result(is_pair)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag_row, ifrag_col
    integer :: axis
    logical :: axis_ok(3)

    is_pair = .false.
    if (ifrag_row == ifrag_col) then
      is_pair = .true.
      return
    end if
    if (allocated(dg_frag%momentum_neighbor_pair_cache)) then
      if (ifrag_row >= 1 .and. ifrag_row <= size(dg_frag%momentum_neighbor_pair_cache, 1) .and. &
          ifrag_col >= 1 .and. ifrag_col <= size(dg_frag%momentum_neighbor_pair_cache, 2)) then
        is_pair = dg_frag%momentum_neighbor_pair_cache(ifrag_row, ifrag_col)
        return
      end if
    end if

    do axis = 1, 3
      axis_ok(axis) = is_momentum_neighbor_axis(dg_frag%lgnum_total(axis), &
        dg_frag%ixyz_frag(axis, ifrag_row), dg_frag%nxyz_domain(axis, ifrag_row), &
        dg_frag%ixyz_frag(axis, ifrag_col), dg_frag%nxyz_domain(axis, ifrag_col))
    end do

    is_pair = all(axis_ok)
  end function is_momentum_neighbor_pair

  subroutine ensure_momentum_neighbor_pair_cache(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: ifrag_row, ifrag_col, axis
    logical :: axis_ok(3)

    if (allocated(dg_frag%momentum_neighbor_pair_cache)) return
    allocate(dg_frag%momentum_neighbor_pair_cache(dg_frag%n_frag, dg_frag%n_frag))
    dg_frag%momentum_neighbor_pair_cache(:, :) = .false.
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (ifrag_row == ifrag_col) then
          dg_frag%momentum_neighbor_pair_cache(ifrag_row, ifrag_col) = .true.
        else
          do axis = 1, 3
            axis_ok(axis) = is_momentum_neighbor_axis(dg_frag%lgnum_total(axis), &
              dg_frag%ixyz_frag(axis, ifrag_row), dg_frag%nxyz_domain(axis, ifrag_row), &
              dg_frag%ixyz_frag(axis, ifrag_col), dg_frag%nxyz_domain(axis, ifrag_col))
          end do
          dg_frag%momentum_neighbor_pair_cache(ifrag_row, ifrag_col) = all(axis_ok)
        end if
      end do
    end do
  end subroutine ensure_momentum_neighbor_pair_cache

  integer function map_global_to_phi_box_coord_ham(ig, lb, ub, lgtot) result(iloc)
    implicit none
    integer, intent(in) :: ig, lb, ub, lgtot

    iloc = modulo(ig - 1, lgtot) + 1
    if (iloc < lb) then
      iloc = iloc + ((lb - iloc + lgtot - 1) / lgtot) * lgtot
    end if
    if (iloc > ub) then
      iloc = iloc - ((iloc - ub + lgtot - 1) / lgtot) * lgtot
    end if
    if (iloc < lb .or. iloc > ub) iloc = 0
  end function map_global_to_phi_box_coord_ham

  integer function map_global_to_phi_box_coord_ham_fragment(dg_frag, ifrag, axis, ig, lb, ub) result(iloc)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, axis, ig, lb, ub
    integer :: ig_wrap, support_lo, support_len

    if (dg_frag%use_buffered_basis) then
      support_lo = dg_frag%basis_support_lo(axis, ifrag)
      support_len = dg_frag%basis_support_hi(axis, ifrag) - dg_frag%basis_support_lo(axis, ifrag) + 1
      ig_wrap = support_lo + modulo(ig - support_lo, support_len)
      iloc = map_global_to_phi_box_coord_ham(ig_wrap, lb, ub, dg_frag%lgnum_total(axis))
    else
      iloc = map_global_to_phi_box_coord_ham(ig, lb, ub, dg_frag%lgnum_total(axis))
    end if
  end function map_global_to_phi_box_coord_ham_fragment

  integer function wrap_support_local_index(iloc, nloc) result(iwrap)
    implicit none
    integer, intent(in) :: iloc, nloc
    iwrap = modulo(iloc - 1, nloc) + 1
  end function wrap_support_local_index

  subroutine fill_buffered_support_local_box(dg_frag, ifrag, i_local, jo, iorg, nsup, phi_local)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, i_local, jo
    integer, intent(in) :: iorg(3), nsup(3)
    real(8), intent(out) :: phi_local(-4:,-4:,-4:)
    integer :: lx, ly, lz, sx, sy, sz
    integer :: gx, gy, gz, bx, by, bz
    integer :: phi_lb1, phi_lb2, phi_lb3, phi_ub1, phi_ub2, phi_ub3

    phi_lb1 = lbound(dg_frag%phi_frag, 1)
    phi_lb2 = lbound(dg_frag%phi_frag, 2)
    phi_lb3 = lbound(dg_frag%phi_frag, 3)
    phi_ub1 = ubound(dg_frag%phi_frag, 1)
    phi_ub2 = ubound(dg_frag%phi_frag, 2)
    phi_ub3 = ubound(dg_frag%phi_frag, 3)

    do lz = -4, nsup(3) + 4
      sz = wrap_support_local_index(lz, nsup(3))
      gz = iorg(3) + sz - 1
      bz = map_global_to_phi_box_coord_ham_fragment(dg_frag, ifrag, 3, gz, phi_lb3, phi_ub3)
      do ly = -4, nsup(2) + 4
        sy = wrap_support_local_index(ly, nsup(2))
        gy = iorg(2) + sy - 1
        by = map_global_to_phi_box_coord_ham_fragment(dg_frag, ifrag, 2, gy, phi_lb2, phi_ub2)
        do lx = -4, nsup(1) + 4
          sx = wrap_support_local_index(lx, nsup(1))
          gx = iorg(1) + sx - 1
          bx = map_global_to_phi_box_coord_ham_fragment(dg_frag, ifrag, 1, gx, phi_lb1, phi_ub1)
          if (bx == 0 .or. by == 0 .or. bz == 0) then
            phi_local(lx, ly, lz) = 0.0d0
          else
            phi_local(lx, ly, lz) = dg_frag%phi_frag(bx, by, bz, jo, i_local)
          end if
        end do
      end do
    end do
  end subroutine fill_buffered_support_local_box

  integer function find_momentum_block(dg_frag, ifrag_row, ifrag_col) result(iblk)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag_row, ifrag_col

    iblk = 0
    if (.not. allocated(dg_frag%momentum_block_map)) return
    if (ifrag_row < 1 .or. ifrag_row > size(dg_frag%momentum_block_map, 1)) return
    if (ifrag_col < 1 .or. ifrag_col > size(dg_frag%momentum_block_map, 2)) return
    iblk = dg_frag%momentum_block_map(ifrag_row, ifrag_col)
  end function find_momentum_block

  integer function find_matrix_block(block_map, ifrag_row, ifrag_col) result(iblk)
    implicit none
    integer, intent(in) :: block_map(:, :)
    integer, intent(in) :: ifrag_row, ifrag_col

    iblk = 0
    if (ifrag_row < 1 .or. ifrag_row > size(block_map, 1)) return
    if (ifrag_col < 1 .or. ifrag_col > size(block_map, 2)) return
    iblk = block_map(ifrag_row, ifrag_col)
  end function find_matrix_block

  logical function fragment_row_is_locally_owned(dg_frag, ifrag_row) result(is_local)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag_row
    integer :: ispin, io, global_idx

    is_local = .false.
    if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) return
    if (.not. allocated(dg_frag%coef_owner)) return
    if (.not. allocated(dg_frag%n_basis)) return
    if (.not. allocated(dg_frag%index_basis)) return

    do ispin = 1, dg_frag%nspin
      if (dg_frag%n_basis(ifrag_row, ispin) <= 0) cycle
      do io = 1, min(dg_frag%n_basis(ifrag_row, ispin), size(dg_frag%index_basis, 1))
        global_idx = dg_frag%index_basis(io, ifrag_row, ispin)
        if (global_idx < 1 .or. global_idx > size(dg_frag%coef_owner, 1)) cycle
        if (dg_frag%coef_owner(global_idx, ispin) /= dg_frag%id) cycle
        is_local = .true.
        return
      end do
    end do
  end function fragment_row_is_locally_owned

  logical function fragment_row_has_single_owner(dg_frag, ifrag_row, owner_rank) result(has_single_owner)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag_row
    integer, intent(out) :: owner_rank
    integer :: ispin, io, global_idx, owner_candidate

    has_single_owner = .false.
    owner_rank = -1
    if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) return
    if (.not. allocated(dg_frag%coef_owner)) return
    if (.not. allocated(dg_frag%n_basis)) return
    if (.not. allocated(dg_frag%index_basis)) return

    owner_candidate = -1
    do ispin = 1, dg_frag%nspin
      if (dg_frag%n_basis(ifrag_row, ispin) <= 0) cycle
      do io = 1, min(dg_frag%n_basis(ifrag_row, ispin), size(dg_frag%index_basis, 1))
        global_idx = dg_frag%index_basis(io, ifrag_row, ispin)
        if (global_idx < 1 .or. global_idx > size(dg_frag%coef_owner, 1)) cycle
        if (dg_frag%coef_owner(global_idx, ispin) < 0) cycle
        if (owner_candidate < 0) then
          owner_candidate = dg_frag%coef_owner(global_idx, ispin)
        else if (owner_candidate /= dg_frag%coef_owner(global_idx, ispin)) then
          owner_rank = -1
          return
        end if
      end do
    end do

    if (owner_candidate < 0) return
    owner_rank = owner_candidate
    has_single_owner = .true.
  end function fragment_row_has_single_owner

  logical function momentum_blocks_use_row_local_storage(dg_frag) result(use_row_local)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer :: ifrag, owner_rank
    logical :: has_basis

    use_row_local = .false.
    if (.not. allocated(dg_frag%coef_owner)) return
    if (.not. allocated(dg_frag%n_basis)) return

    do ifrag = 1, dg_frag%n_frag
      has_basis = any(dg_frag%n_basis(ifrag, 1:dg_frag%nspin) > 0)
      if (.not. has_basis) cycle
      if (.not. fragment_row_has_single_owner(dg_frag, ifrag, owner_rank)) return
    end do
    use_row_local = .true.
  end function momentum_blocks_use_row_local_storage

  subroutine init_matrix_blocks(dg_frag, blocks, block_map, n_blocks, row_local_only)
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(matrix_block_info), allocatable, intent(inout) :: blocks(:)
    integer, allocatable, intent(inout) :: block_map(:, :)
    integer, intent(out) :: n_blocks
    logical, intent(in), optional :: row_local_only
    integer :: ifrag_row, ifrag_col, iblk, nrow_max, ncol_max, owner_rank
    logical :: use_row_local

    if (allocated(blocks)) then
      do iblk = 1, size(blocks)
        if (allocated(blocks(iblk)%val)) deallocate(blocks(iblk)%val)
      end do
      deallocate(blocks)
    end if
    if (allocated(block_map)) deallocate(block_map)
    call ensure_momentum_neighbor_pair_cache(dg_frag)
    use_row_local = .false.
    if (present(row_local_only)) use_row_local = row_local_only

    n_blocks = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (use_row_local) then
          if (dg_frag%use_buffered_basis) then
            if (ifrag_row < dg_frag%ifrag_start .or. ifrag_row > dg_frag%ifrag_end) cycle
          else
            if (.not. fragment_row_has_single_owner(dg_frag, ifrag_row, owner_rank)) cycle
            if (owner_rank /= dg_frag%id) cycle
          end if
        end if
        if (is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) n_blocks = n_blocks + 1
      end do
    end do
    if (n_blocks <= 0) return

    allocate(blocks(n_blocks))
    allocate(block_map(dg_frag%n_frag, dg_frag%n_frag))
    block_map = 0

    iblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (use_row_local) then
          if (dg_frag%use_buffered_basis) then
            if (ifrag_row < dg_frag%ifrag_start .or. ifrag_row > dg_frag%ifrag_end) cycle
          else
            if (.not. fragment_row_has_single_owner(dg_frag, ifrag_row, owner_rank)) cycle
            if (owner_rank /= dg_frag%id) cycle
          end if
        end if
        if (.not. is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) cycle
        iblk = iblk + 1
        nrow_max = max(1, maxval(dg_frag%n_basis(ifrag_row, 1:dg_frag%nspin)))
        ncol_max = max(1, maxval(dg_frag%n_basis(ifrag_col, 1:dg_frag%nspin)))
        block_map(ifrag_row, ifrag_col) = iblk
        blocks(iblk)%ifrag_row = ifrag_row
        blocks(iblk)%ifrag_col = ifrag_col
        blocks(iblk)%nrow_max = nrow_max
        blocks(iblk)%ncol_max = ncol_max
        allocate(blocks(iblk)%val(nrow_max, ncol_max, dg_frag%nspin))
        blocks(iblk)%val = 0.0d0
      end do
    end do
  end subroutine init_matrix_blocks

  subroutine sync_dense_matrix_to_blocks(dg_frag, mat, blocks, block_map)
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    real(8), intent(in) :: mat(:, :, :)
    type(matrix_block_info), intent(inout) :: blocks(:)
    integer, intent(in) :: block_map(:, :)
    integer :: ifrag_row, ifrag_col, iblk, ispin, ii, jj, ig_i, ig_j
    integer :: nrow, ncol

    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        iblk = find_matrix_block(block_map, ifrag_row, ifrag_col)
        if (iblk <= 0) cycle
        blocks(iblk)%val(:, :, :) = 0.0d0
        do ispin = 1, dg_frag%nspin
          nrow = dg_frag%n_basis(ifrag_row, ispin)
          ncol = dg_frag%n_basis(ifrag_col, ispin)
          if (nrow <= 0 .or. ncol <= 0) cycle
          do jj = 1, ncol
            ig_j = dg_frag%index_basis(jj, ifrag_col, ispin)
            if (ig_j < 1 .or. ig_j > size(mat, 2)) cycle
            do ii = 1, nrow
              ig_i = dg_frag%index_basis(ii, ifrag_row, ispin)
              if (ig_i < 1 .or. ig_i > size(mat, 1)) cycle
              blocks(iblk)%val(ii, jj, ispin) = mat(ig_i, ig_j, ispin)
            end do
          end do
        end do
      end do
    end do
  end subroutine sync_dense_matrix_to_blocks

  subroutine sync_blocks_to_dense_matrix(dg_frag, blocks, block_map, mat)
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(in) :: blocks(:)
    integer, intent(in) :: block_map(:, :)
    real(8), intent(inout) :: mat(:, :, :)
    integer :: ifrag_row, ifrag_col, iblk, ispin, ii, jj, ig_i, ig_j
    integer :: nrow, ncol, idx_ii, idx_jj, valid_row_count, valid_col_count
    integer, allocatable :: row_gid(:), col_gid(:), valid_row_ids(:), valid_col_ids(:)

    mat(:, :, :) = 0.0d0
    allocate(row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1)))
    allocate(valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1)))
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        iblk = find_matrix_block(block_map, ifrag_row, ifrag_col)
        if (iblk <= 0) cycle
        do ispin = 1, dg_frag%nspin
          nrow = dg_frag%n_basis(ifrag_row, ispin)
          ncol = dg_frag%n_basis(ifrag_col, ispin)
          if (nrow <= 0 .or. ncol <= 0) cycle
          valid_row_count = 0
          do ii = 1, nrow
            row_gid(ii) = dg_frag%index_basis(ii, ifrag_row, ispin)
            if (row_gid(ii) < 1 .or. row_gid(ii) > size(mat, 1)) cycle
            valid_row_count = valid_row_count + 1
            valid_row_ids(valid_row_count) = ii
          end do
          valid_col_count = 0
          do jj = 1, ncol
            col_gid(jj) = dg_frag%index_basis(jj, ifrag_col, ispin)
            if (col_gid(jj) < 1 .or. col_gid(jj) > size(mat, 2)) cycle
            valid_col_count = valid_col_count + 1
            valid_col_ids(valid_col_count) = jj
          end do
          do idx_jj = 1, valid_col_count
            jj = valid_col_ids(idx_jj)
            ig_j = col_gid(jj)
            do idx_ii = 1, valid_row_count
              ii = valid_row_ids(idx_ii)
              ig_i = row_gid(ii)
              mat(ig_i, ig_j, ispin) = blocks(iblk)%val(ii, jj, ispin)
            end do
          end do
        end do
      end do
    end do
    deallocate(row_gid, col_gid, valid_row_ids, valid_col_ids)
  end subroutine sync_blocks_to_dense_matrix

  subroutine init_momentum_blocks(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: ifrag_row, ifrag_col, nblk, iblk
    integer :: nrow_max, ncol_max, owner_rank
    logical :: use_row_local

    if (allocated(dg_frag%momentum_blocks)) then
      do iblk = 1, size(dg_frag%momentum_blocks)
        if (allocated(dg_frag%momentum_blocks(iblk)%val)) deallocate(dg_frag%momentum_blocks(iblk)%val)
      end do
      deallocate(dg_frag%momentum_blocks)
    end if
    if (allocated(dg_frag%momentum_block_map)) deallocate(dg_frag%momentum_block_map)
    call ensure_momentum_neighbor_pair_cache(dg_frag)
    use_row_local = momentum_blocks_use_row_local_storage(dg_frag)
    if (dg_frag%use_buffered_basis) use_row_local = .true.

    nblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (use_row_local) then
          if (dg_frag%use_buffered_basis) then
            if (ifrag_row < dg_frag%ifrag_start .or. ifrag_row > dg_frag%ifrag_end) cycle
          else
            if (.not. fragment_row_has_single_owner(dg_frag, ifrag_row, owner_rank)) cycle
            if (owner_rank /= dg_frag%id) cycle
          end if
        end if
        if (is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) nblk = nblk + 1
      end do
    end do

    dg_frag%n_momentum_blocks = nblk
    if (nblk <= 0) return
    allocate(dg_frag%momentum_blocks(nblk))
    allocate(dg_frag%momentum_block_map(dg_frag%n_frag, dg_frag%n_frag))
    dg_frag%momentum_block_map = 0

    iblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (use_row_local) then
          if (dg_frag%use_buffered_basis) then
            if (ifrag_row < dg_frag%ifrag_start .or. ifrag_row > dg_frag%ifrag_end) cycle
          else
            if (.not. fragment_row_has_single_owner(dg_frag, ifrag_row, owner_rank)) cycle
            if (owner_rank /= dg_frag%id) cycle
          end if
        end if
        if (.not. is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) cycle
        iblk = iblk + 1
        dg_frag%momentum_block_map(ifrag_row, ifrag_col) = iblk
        dg_frag%momentum_blocks(iblk)%ifrag_row = ifrag_row
        dg_frag%momentum_blocks(iblk)%ifrag_col = ifrag_col
        nrow_max = max(1, maxval(dg_frag%n_basis(ifrag_row, 1:dg_frag%nspin)))
        ncol_max = max(1, maxval(dg_frag%n_basis(ifrag_col, 1:dg_frag%nspin)))
        dg_frag%momentum_blocks(iblk)%nrow_max = nrow_max
        dg_frag%momentum_blocks(iblk)%ncol_max = ncol_max
        allocate(dg_frag%momentum_blocks(iblk)%val(3, nrow_max, ncol_max, dg_frag%nspin))
        dg_frag%momentum_blocks(iblk)%val = 0.0d0
      end do
    end do
  end subroutine init_momentum_blocks

  subroutine reduce_matrix_blocks(dg_frag, blocks, label, icomm_reduce)
    use communication, only: comm_is_root, comm_summation, comm_get_max
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(inout) :: blocks(:)
    character(*), intent(in) :: label
    integer, intent(in) :: icomm_reduce
    integer, parameter :: reduce_chunk_size = 262144
    real(8), allocatable :: send_block(:), recv_block(:)
    integer :: iblk, ispin, ii, jj
    integer :: nrow, ncol, block_size, max_block_size, total_active_size
    integer :: total_active_min, total_active_max, max_block_size_global
    integer :: chunk_begin, chunk_count, offset_flat
    logical :: packed_diag
    logical :: enable_reduce_trace
    character(16) :: env_reduce_trace
    integer :: env_status

    enable_reduce_trace = .false.
    env_reduce_trace = ''
    call get_environment_variable('SALMON_DG_HMAT_REDUCE_TRACE', env_reduce_trace, status=env_status)
    if (env_status == 0) then
      select case(trim(adjustl(env_reduce_trace)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enable_reduce_trace = .true.
      end select
    end if

    max_block_size = 0
    total_active_size = 0
    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        packed_diag = (blocks(iblk)%ifrag_row == blocks(iblk)%ifrag_col) .and. (nrow == ncol)
        if (packed_diag) then
          block_size = nrow * (nrow + 1) / 2
        else
          block_size = nrow * ncol
        end if
        max_block_size = max(max_block_size, block_size)
        total_active_size = total_active_size + block_size
      end do
    end do

    max_block_size_global = max_block_size
    call comm_get_max(max_block_size_global, icomm_reduce)
    total_active_max = total_active_size
    call comm_get_max(total_active_max, icomm_reduce)
    total_active_min = -total_active_size
    call comm_get_max(total_active_min, icomm_reduce)
    total_active_min = -total_active_min

    if (total_active_min /= total_active_max) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] Hamiltonian block size mismatch: label=", &
        trim(label), " rank=", dg_frag%id, " local=", total_active_size, &
        " min=", total_active_min, " max=", total_active_max
      flush(6)
      stop 1
    end if

    if (enable_reduce_trace .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "        hamiltonian block reduce begin: label=", trim(label), &
        " total_active=", total_active_size, " max_block=", max_block_size_global, &
        " chunk_size=", reduce_chunk_size
      flush(6)
    end if

    if (total_active_size <= 0) return
    allocate(send_block(total_active_size), recv_block(total_active_size))

    offset_flat = 1
    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        packed_diag = (blocks(iblk)%ifrag_row == blocks(iblk)%ifrag_col) .and. (nrow == ncol)
        do jj = 1, ncol
          if (packed_diag) then
            do ii = 1, min(jj, nrow)
              send_block(offset_flat) = blocks(iblk)%val(ii, jj, ispin)
              offset_flat = offset_flat + 1
            end do
          else
            do ii = 1, nrow
              send_block(offset_flat) = blocks(iblk)%val(ii, jj, ispin)
              offset_flat = offset_flat + 1
            end do
          end if
        end do
      end do
    end do

    chunk_begin = 1
    do while (chunk_begin <= total_active_size)
      chunk_count = min(reduce_chunk_size, total_active_size - chunk_begin + 1)
      call comm_summation(send_block(chunk_begin:chunk_begin + chunk_count - 1), &
                          recv_block(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, icomm_reduce)
      chunk_begin = chunk_begin + chunk_count
    end do

    offset_flat = 1
    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        packed_diag = (blocks(iblk)%ifrag_row == blocks(iblk)%ifrag_col) .and. (nrow == ncol)
        do jj = 1, ncol
          if (packed_diag) then
            do ii = 1, min(jj, nrow)
              blocks(iblk)%val(ii, jj, ispin) = recv_block(offset_flat)
              if (ii /= jj) blocks(iblk)%val(jj, ii, ispin) = recv_block(offset_flat)
              offset_flat = offset_flat + 1
            end do
          else
            do ii = 1, nrow
              blocks(iblk)%val(ii, jj, ispin) = recv_block(offset_flat)
              offset_flat = offset_flat + 1
            end do
          end if
        end do
      end do
    end do

    deallocate(send_block, recv_block)
    if (enable_reduce_trace .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0)') "        hamiltonian block reduce done: label=", trim(label), &
        " total_active=", total_active_size
      flush(6)
    end if
  end subroutine reduce_matrix_blocks

  !=======================================================================
  ! Calculate initial Hamiltonian matrix from basis functions
  !
  ! Includes halo (ghost cell) exchange for accurate boundary treatment.
  ! System boundaries use PERIODIC boundary conditions (full system is periodic).
  ! Fragment boundaries are handled via MPI communication between neighboring fragments.
  ! The real-space fragment basis itself is shared across spin channels in the
  ! present non-SOI DG path; the nspin axis here labels spin-resolved projected
  ! matrices and basis indexing, not separate copies of phi_frag for each spin.
  !=======================================================================
  subroutine calculate_hamiltonian_matrix(dg_frag, system, lg, mg, stencil, &
                                         Vh, Vxc, Vpsl, pp, ppg)
    use structures
    use communication, only: comm_is_root, comm_summation, comm_get_max
    use parallelization, only: nproc_size_global
    use rt_dg_plane_wave, only: diagonalize_mixed_basis
    use rt_dg_fragment_ops, only: copy_matrix_blocks_metric_to_complex_dense, symmetrize_real_matrix_blocks
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: lg, mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_scalar),         intent(in)    :: Vh, Vxc(:), Vpsl
    type(s_pp_info),        intent(in)    :: pp
    type(s_pp_grid),        intent(in)    :: ppg
    
    integer :: ifrag, ispin, io, jo, i_local, nbf, nbf_raw, ig_i, ig_j
    integer :: ifrag_chk, i_local_chk, ix_chk, iy_chk, iz_chk, istate_chk
    integer :: nstate_chk_max, lg1_chk, lg2_chk, lg3_chk
    integer :: phi_lb1_chk, phi_ub1_chk, phi_lb2_chk, phi_ub2_chk, phi_lb3_chk, phi_ub3_chk
    integer :: iorg_chk(3), ndom_chk(3)
    integer :: gx, gy, gz, bx, by, bz
    integer :: ndom(3)
    integer :: i_halo, i_halo_active, jfrag, n_basis_halo, n_basis_halo_max, n_active_halo
    integer :: l(3), halo_send_idx(3), halo_recv_idx(3)
    integer :: face_axis, face_dir
    integer :: npts_halo, ipt
    integer :: face_idx_local, face_idx_remote
    integer :: ix_face, iy_face, iz_face
    real(8) :: hvol
    real(8) :: halo_integral_t, halo_integral_h, halo_integral_m, halo_integral_m_avgavg, halo_integral_m_avgjump
    real(8) :: halo_integral_t_sum, halo_integral_h_sum, halo_integral_m_sum, t_point, h_point, halo_val
    real(8) :: face_area, face_penalty_alpha, phi_local_face, phi_remote_face, phi_face_avg, phi_face_jump
    real(8) :: dphi_local_n, dphi_remote_n, dphi_flux_n, dphi_face_avg, dphi_face_jump
    real(8) :: sipg_ab_sum, sipg_ab_cap, sipg_ab_scale, sipg_ab_ref, sipg_ab_norm
    real(8) :: max_p
    real(8) :: Ac_zero(3)
    real(8) :: hmat_dense_mb, phi_frag_mb, halo_buf_mb, overlap_dense_mb, momentum_dense_mb
    integer :: is(3), ie(3)
    integer(8) :: byte_count
    integer(8) :: i_halo_chk
    real(8), allocatable :: T_phi(:,:,:)  ! Kinetic energy operator applied to basis (fragment-local)
    real(8), allocatable :: H_phi(:,:,:)  ! Hamiltonian-applied field H|phi_j> = T|phi_j> + V|phi_j> (fragment-local)
    real(8), allocatable :: V_total(:,:,:)  ! Total potential V = Vpsl + Vh + Vxc
    real(8), allocatable :: partial_t(:), partial_h(:), reduced_t(:), reduced_h(:)
    real(8), allocatable :: partial_th(:), reduced_th(:)
    real(8), allocatable :: halo_partial_t(:), halo_partial_h(:), halo_partial_m(:), halo_reduced_t(:), halo_reduced_h(:), halo_reduced_m(:)
    real(8), allocatable :: halo_partial_a(:), halo_partial_b(:), halo_partial_c(:), halo_reduced_a(:), halo_reduced_b(:), halo_reduced_c(:)
    real(8), allocatable :: halo_reduce_pair(:), halo_reduce_sum(:)
    real(8), allocatable :: halo_t_point_buf(:), halo_h_point_buf(:)
    integer, allocatable :: halo_active_list(:), halo_active_nbf(:), halo_active_iblk(:), halo_active_iblk_rev(:)
    integer, allocatable :: halo_active_face_axis(:), halo_active_face_dir(:), halo_active_face_idx_local(:), halo_active_face_idx_remote(:)
    type(matrix_block_info), allocatable :: H_diag_blocks(:), H_kin_diag_blocks(:)
    integer :: n_local_diag, nbf_max, i_diag, iblk, iblk_rev, nbf_diag, nbf_comm
    integer :: loc_s_dbg(3), loc_e_dbg(3)
    integer :: n_metric
    integer :: env_status, env_len
    logical :: release_dense_fragment_ops
    logical :: enable_startup_lowdin
    logical :: enable_startup_stationary_projection
    logical :: has_overlap_dbg
    logical :: need_halo_alloc
    complex(8), allocatable :: H_metric_ref(:,:)
    real(8) :: halo_integral_a_sum, halo_integral_b_sum, halo_integral_c_sum
    real(8) :: sipg_h_face_a, sipg_h_face_b, sipg_h_face_c
    character(len=32) :: env_sipg_disable_b_term, env_sipg_disable_c_term, env_sipg_enable_halo, env_sipg_use_avg_trace, env_sipg_use_flux_form
    character(len=32) :: env_sipg_use_weak_form, env_sipg_ab_norm_max
    character(len=32) :: env_startup_lowdin, env_startup_stationary_projection
    logical :: sipg_disable_b_term, sipg_disable_c_term, sipg_enable_fragment_halo_exchange, sipg_use_avg_trace, sipg_use_flux_form
    logical :: sipg_use_weak_form
    real(8) :: sipg_ab_norm_max
    
    release_dense_fragment_ops = (.not. dg_frag%yn_adaptive_basis) .and. &
      ((.not. dg_frag%use_plane_wave_basis) .or. dg_frag%n_plane_waves <= 0)
    env_sipg_disable_b_term = ''
    env_sipg_disable_c_term = ''
    env_sipg_enable_halo = ''
    env_sipg_use_avg_trace = ''
    env_sipg_use_flux_form = ''
    env_sipg_use_weak_form = ''
    env_sipg_ab_norm_max = ''
    sipg_disable_b_term = .false.
    sipg_disable_c_term = .false.
    sipg_enable_fragment_halo_exchange = .false.
    sipg_use_avg_trace = .false.
    sipg_use_flux_form = .false.
    sipg_use_weak_form = .false.
    sipg_ab_norm_max = -1.0d0
    enable_startup_lowdin = .false.
    enable_startup_stationary_projection = .false.
    env_startup_lowdin = ''
    env_startup_stationary_projection = ''
    call get_environment_variable('SALMON_DG_SIPG_DISABLE_B_TERM', env_sipg_disable_b_term, status=env_status)
    if (env_status == 0) then
      select case (trim(adjustl(env_sipg_disable_b_term)))
      case ('1', 'true', 'TRUE', 'on', 'ON', 'yes', 'YES')
        sipg_disable_b_term = .true.
      end select
    end if
    call get_environment_variable('SALMON_DG_SIPG_DISABLE_C_TERM', env_sipg_disable_c_term, status=env_status)
    if (env_status == 0) then
      select case (trim(adjustl(env_sipg_disable_c_term)))
      case ('1', 'true', 'TRUE', 'on', 'ON', 'yes', 'YES')
        sipg_disable_c_term = .true.
      end select
    end if
    call get_environment_variable('SALMON_DG_ENABLE_FRAGMENT_HALO_EXCHANGE', env_sipg_enable_halo, status=env_status)
    if (env_status == 0) then
      select case (trim(adjustl(env_sipg_enable_halo)))
      case ('1', 'true', 'TRUE', 'on', 'ON', 'yes', 'YES')
        sipg_enable_fragment_halo_exchange = .true.
      end select
    end if
    call get_environment_variable('SALMON_DG_SIPG_USE_AVG_TRACE', env_sipg_use_avg_trace, status=env_status)
    if (env_status == 0) then
      select case (trim(adjustl(env_sipg_use_avg_trace)))
      case ('1', 'true', 'TRUE', 'on', 'ON', 'yes', 'YES')
        sipg_use_avg_trace = .true.
      end select
    end if
    call get_environment_variable('SALMON_DG_SIPG_USE_FLUX_FORM', env_sipg_use_flux_form, status=env_status)
    if (env_status == 0) then
      select case (trim(adjustl(env_sipg_use_flux_form)))
      case ('1', 'true', 'TRUE', 'on', 'ON', 'yes', 'YES')
        sipg_use_flux_form = .true.
      end select
    end if
    call get_environment_variable('SALMON_DG_SIPG_USE_WEAK_FORM', env_sipg_use_weak_form, status=env_status)
    if (env_status == 0) then
      select case (trim(adjustl(env_sipg_use_weak_form)))
      case ('1', 'true', 'TRUE', 'on', 'ON', 'yes', 'YES')
        sipg_use_weak_form = .true.
      end select
    end if
    env_status = 1
    env_len = 0
    call get_environment_variable('SALMON_DG_SIPG_AB_NORM_MAX', env_sipg_ab_norm_max, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_sipg_ab_norm_max, *, iostat=env_status) sipg_ab_norm_max
      if (env_status /= 0) sipg_ab_norm_max = -1.0d0
    end if
    call get_environment_variable('SALMON_DG_STARTUP_LOWDIN', env_startup_lowdin, status=env_status)
    if (env_status == 0) then
      select case (trim(adjustl(env_startup_lowdin)))
      case ('1', 'true', 'TRUE', 'on', 'ON', 'yes', 'YES')
        enable_startup_lowdin = .true.
      end select
    end if
    call get_environment_variable('SALMON_DG_STARTUP_STATIONARY_PROJECTION', env_startup_stationary_projection, status=env_status)
    if (env_status == 0) then
      select case (trim(adjustl(env_startup_stationary_projection)))
      case ('1', 'true', 'TRUE', 'on', 'ON', 'yes', 'YES')
        enable_startup_stationary_projection = .true.
      end select
    end if
    if (.not. enable_startup_lowdin) enable_startup_stationary_projection = .false.
    if (.not. dg_frag%has_real_space_basis) then
      if (.not. allocated(dg_frag%H_mat)) then
        allocate(dg_frag%H_mat(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
      end if
      dg_frag%H_mat = 0.0d0
      return
    end if
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*)
      write(*,*) "=== Preparing Hamiltonian Matrix ==="
    end if

    nstate_chk_max = min(dg_frag%nstate_frag, size(dg_frag%phi_frag, 4))
    lg1_chk = dg_frag%lgnum_total(1)
    lg2_chk = dg_frag%lgnum_total(2)
    lg3_chk = dg_frag%lgnum_total(3)
    phi_lb1_chk = lbound(dg_frag%phi_frag, 1)
    phi_ub1_chk = ubound(dg_frag%phi_frag, 1)
    phi_lb2_chk = lbound(dg_frag%phi_frag, 2)
    phi_ub2_chk = ubound(dg_frag%phi_frag, 2)
    phi_lb3_chk = lbound(dg_frag%phi_frag, 3)
    phi_ub3_chk = ubound(dg_frag%phi_frag, 3)
    
    ! Step 1: Calculate momentum matrix elements (transition moments)
    ! Required for velocity gauge A·p coupling
    if (.not. allocated(dg_frag%momentum_blocks) .and. .not. allocated(dg_frag%momentum_mat)) then
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "  [1/3] Calculating momentum matrix elements (p_ij)..."
        write(*,*) "        Using 4th-order finite difference stencil"
      end if
      call calculate_momentum_matrix(dg_frag, system, mg, stencil)

      call calculate_overlap_matrix(dg_frag, system, mg)
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "        Momentum matrix calculated (for A·p coupling)"
        write(*,*) "        Overlap matrix S calculated (for generalized propagation)"
      end if
    else
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "  [1/3] Momentum matrix already available"
      end if
      if (.not. allocated(dg_frag%S_mat) .and. .not. allocated(dg_frag%S_mat_blocks) .and. &
          .not. allocated(dg_frag%S_mat_prop_blocks)) then
        call calculate_overlap_matrix(dg_frag, system, mg)
      end if
    end if

    ! Step 2: Allocate Hamiltonian matrix
    !
    ! Current Eq.(4) mapping in this routine:
    !   - volume kinetic term:    build_hpsi_for_basis() / apply_kinetic_to_basis()
    !   - volume local potential: build_hpsi_for_basis()
    !   - nonlocal PP term:       handled later in time evolution, not in H_mat here
    !   - DG face terms:          not yet assembled in the legacy path below
    !
    ! The next implementation stage adds the missing kinetic face terms onto the
    ! fragment-fragment blocks while preserving the current volume assembly.
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "  [2/3] Constructing Hamiltonian matrix H = T + V..."
    end if
    
    if (allocated(dg_frag%H_mat)) deallocate(dg_frag%H_mat)
    if (allocated(dg_frag%H_mat_kinetic)) deallocate(dg_frag%H_mat_kinetic)
    hmat_dense_mb = 0.0d0
    overlap_dense_mb = 0.0d0
    momentum_dense_mb = 0.0d0
    phi_frag_mb = 0.0d0
    halo_buf_mb = 0.0d0
    if (allocated(dg_frag%H_mat)) then
      byte_count = int(size(dg_frag%H_mat, 1), kind=8) * int(size(dg_frag%H_mat, 2), kind=8) * &
                   int(size(dg_frag%H_mat, 3), kind=8) * 8_8
      hmat_dense_mb = hmat_dense_mb + dble(byte_count) / (1024.0d0 * 1024.0d0)
    end if
    if (allocated(dg_frag%H_mat_kinetic)) then
      byte_count = int(size(dg_frag%H_mat_kinetic, 1), kind=8) * int(size(dg_frag%H_mat_kinetic, 2), kind=8) * &
                   int(size(dg_frag%H_mat_kinetic, 3), kind=8) * 8_8
      hmat_dense_mb = hmat_dense_mb + dble(byte_count) / (1024.0d0 * 1024.0d0)
    end if
    if (allocated(dg_frag%S_mat)) then
      byte_count = int(size(dg_frag%S_mat, 1), kind=8) * int(size(dg_frag%S_mat, 2), kind=8) * &
                   int(size(dg_frag%S_mat, 3), kind=8) * 8_8
      overlap_dense_mb = overlap_dense_mb + dble(byte_count) / (1024.0d0 * 1024.0d0)
    end if
    if (allocated(dg_frag%S_mat_prop)) then
      byte_count = int(size(dg_frag%S_mat_prop, 1), kind=8) * int(size(dg_frag%S_mat_prop, 2), kind=8) * &
                   int(size(dg_frag%S_mat_prop, 3), kind=8) * 8_8
      overlap_dense_mb = overlap_dense_mb + dble(byte_count) / (1024.0d0 * 1024.0d0)
    end if
    if (allocated(dg_frag%momentum_mat)) then
      byte_count = int(size(dg_frag%momentum_mat, 1), kind=8) * int(size(dg_frag%momentum_mat, 2), kind=8) * &
                   int(size(dg_frag%momentum_mat, 3), kind=8) * int(size(dg_frag%momentum_mat, 4), kind=8) * 8_8
      momentum_dense_mb = dble(byte_count) / (1024.0d0 * 1024.0d0)
    end if
    if (allocated(dg_frag%phi_frag)) then
      byte_count = int(size(dg_frag%phi_frag, 1), kind=8) * int(size(dg_frag%phi_frag, 2), kind=8) * &
                   int(size(dg_frag%phi_frag, 3), kind=8) * int(size(dg_frag%phi_frag, 4), kind=8) * &
                   int(size(dg_frag%phi_frag, 5), kind=8) * 8_8
      phi_frag_mb = dble(byte_count) / (1024.0d0 * 1024.0d0)
    end if
    if (allocated(dg_frag%halo)) then
      do i_halo_chk = 1, size(dg_frag%halo)
        if (allocated(dg_frag%halo(i_halo_chk)%buf_send)) then
          byte_count = int(size(dg_frag%halo(i_halo_chk)%buf_send, 1), kind=8) * &
                       int(size(dg_frag%halo(i_halo_chk)%buf_send, 2), kind=8) * &
                       int(size(dg_frag%halo(i_halo_chk)%buf_send, 3), kind=8) * &
                       int(size(dg_frag%halo(i_halo_chk)%buf_send, 4), kind=8) * &
                       int(size(dg_frag%halo(i_halo_chk)%buf_send, 5), kind=8) * 8_8
          halo_buf_mb = halo_buf_mb + dble(byte_count) / (1024.0d0 * 1024.0d0)
        end if
        if (allocated(dg_frag%halo(i_halo_chk)%buf_recv)) then
          byte_count = int(size(dg_frag%halo(i_halo_chk)%buf_recv, 1), kind=8) * &
                       int(size(dg_frag%halo(i_halo_chk)%buf_recv, 2), kind=8) * &
                       int(size(dg_frag%halo(i_halo_chk)%buf_recv, 3), kind=8) * &
                       int(size(dg_frag%halo(i_halo_chk)%buf_recv, 4), kind=8) * &
                       int(size(dg_frag%halo(i_halo_chk)%buf_recv, 5), kind=8) * 8_8
          halo_buf_mb = halo_buf_mb + dble(byte_count) / (1024.0d0 * 1024.0d0)
        end if
      end do
    end if

    n_local_diag = max(0, dg_frag%ifrag_end - dg_frag%ifrag_start + 1)
    if (n_local_diag > 0) then
      allocate(H_diag_blocks(n_local_diag), H_kin_diag_blocks(n_local_diag))
      do i_diag = 1, n_local_diag
        ifrag = dg_frag%ifrag_start + i_diag - 1
        nbf_max = max(1, maxval(dg_frag%n_basis(ifrag, 1:dg_frag%nspin)))
        H_diag_blocks(i_diag)%ifrag_row = ifrag
        H_diag_blocks(i_diag)%ifrag_col = ifrag
        H_diag_blocks(i_diag)%nrow_max = nbf_max
        H_diag_blocks(i_diag)%ncol_max = nbf_max
        allocate(H_diag_blocks(i_diag)%val(nbf_max, nbf_max, dg_frag%nspin))
        H_diag_blocks(i_diag)%val = 0.0d0
        H_kin_diag_blocks(i_diag)%ifrag_row = ifrag
        H_kin_diag_blocks(i_diag)%ifrag_col = ifrag
        H_kin_diag_blocks(i_diag)%nrow_max = nbf_max
        H_kin_diag_blocks(i_diag)%ncol_max = nbf_max
        allocate(H_kin_diag_blocks(i_diag)%val(nbf_max, nbf_max, dg_frag%nspin))
        H_kin_diag_blocks(i_diag)%val = 0.0d0
      end do
    end if
    
    ! Optionally enable flow-only face communication for Step2 SIPG coupling.
    ! Hamiltonian bulk terms stay fragment-local; only remote face phi/dphi are exchanged.
    if (sipg_enable_fragment_halo_exchange .and. dg_frag%n_halo <= 0) then
      call init_flow_halo_communication(dg_frag)
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') "Fragment flow-face exchange enabled for Step2 by SALMON_DG_ENABLE_FRAGMENT_HALO_EXCHANGE"
      end if
    end if

    if (sipg_enable_fragment_halo_exchange) then
      call exchange_flow_face_halo(dg_frag)
    end if
    
    hvol = system%hvol
    is = mg%is
    ie = mg%ie
    
    allocate(V_total(is(1):ie(1), is(2):ie(2), is(3):ie(3)))
    
    ! Construct total potential: V = Vpsl + Vh + Vxc
    ! Note: This is used for initial H_mat calculation
    do ispin = 1, system%nspin
      call build_total_potential_grid(mg, Vh, Vxc(ispin), Vpsl, V_total)
      
      ! Loop over fragments assigned to this rank
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        ! Calculate Hamiltonian matrix elements for this fragment
        ! H_ij = <φ_i | T + V | φ_j> = T_ij + V_ij
        nbf_raw = dg_frag%n_basis(ifrag, ispin)
        if (nbf_raw < 0) then
          write(*,*) "[FATAL] negative n_basis in Hamiltonian Step 2: rank=", dg_frag%id, &
            " ifrag=", ifrag, " ispin=", ispin, " n_basis=", nbf_raw
          stop 1
        end if
        nbf = min(nbf_raw, dg_frag%nstate_frag)
        if (nbf <= 0) cycle
        if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) then
          write(*,*) "[FATAL] hamiltonian step2 invalid i_local: rank=", dg_frag%id, &
            " ifrag=", ifrag, " i_local=", i_local, " phi_frag_dim5=", size(dg_frag%phi_frag, 5)
          stop 1
        end if
        ndom(:) = dg_frag%nxyz_domain(:, ifrag)
        if (any(ndom <= 0)) then
          write(*,*) "[FATAL] hamiltonian step2 non-positive ndom: rank=", dg_frag%id, &
            " ifrag=", ifrag, " ndom=", ndom
          stop 1
        end if
        allocate(T_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
        allocate(H_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
        if (nbf > size(dg_frag%index_basis, 1)) then
          write(*,*) "[FATAL] hamiltonian n_basis exceeds index_basis dim1: rank=", dg_frag%id, &
            " ifrag=", ifrag, " ispin=", ispin, " n_basis_eff=", nbf, " n_basis_raw=", nbf_raw, &
            " index_basis_dim1=", size(dg_frag%index_basis, 1)
          stop 1
        end if
        if (nbf > size(dg_frag%phi_frag, 4)) then
          write(*,*) "[FATAL] hamiltonian n_basis exceeds phi_frag dim4: rank=", dg_frag%id, &
            " ifrag=", ifrag, " ispin=", ispin, " n_basis_eff=", nbf, " n_basis_raw=", nbf_raw, &
            " phi_frag_dim4=", size(dg_frag%phi_frag, 4)
          stop 1
        end if
        ! Use communicator-invariant length for frag-subgroup reductions to avoid
        ! count mismatch when effective nbf differs by rank-local conditions.
        nbf_comm = dg_frag%nstate_frag
        if (nbf_comm < nbf) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') "[FATAL] step2 invalid nbf_comm: rank=", &
            dg_frag%id, " ifrag=", ifrag, " ispin=", ispin, " nbf=", nbf, " nbf_comm=", nbf_comm
          stop "DG-Fragment RT: nbf_comm smaller than local nbf"
        end if
        allocate(partial_t(nbf_comm), partial_h(nbf_comm), reduced_t(nbf_comm), reduced_h(nbf_comm))
        allocate(partial_th(2 * nbf_comm), reduced_th(2 * nbf_comm))
        call get_fragment_owned_range(dg_frag, ifrag, mg, loc_s_dbg, loc_e_dbg, has_overlap_dbg)
        do jo = 1, nbf
          ig_j = dg_frag%index_basis(jo, ifrag, ispin)
          if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle
          ! Eq.(4) volume contribution only:
          !   H_phi = T(phi_j) + V_local * phi_j
          ! on the fragment core/domain.
          !
          ! The DG kinetic face terms from Eq.(4) are not included here yet and
          ! must be added at the fragment-fragment block assembly stage below.
          call build_hpsi_for_basis(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi)

          ! Calculate matrix elements with all φ_i
          partial_t(:) = 0.0d0
          partial_h(:) = 0.0d0
          !$omp parallel do private(io, ig_i)
          do io = 1, nbf
            ig_i = dg_frag%index_basis(io, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle

            ! Kinetic energy matrix element: T_ij = ∫ φ_i (T|φ_j>) dr
            call integrate_basis_with_field(dg_frag, ifrag, i_local, io, mg, T_phi, hvol, partial_t(io))

            ! Store kinetic part
            call integrate_basis_with_field(dg_frag, ifrag, i_local, io, mg, H_phi, hvol, partial_h(io))

          end do
          !$omp end parallel do
          partial_th(1:nbf_comm) = partial_t(1:nbf_comm)
          partial_th(nbf_comm + 1:2 * nbf_comm) = partial_h(1:nbf_comm)
          call comm_summation(partial_th, reduced_th, 2 * nbf_comm, dg_frag%icomm_frag)
          reduced_t(1:nbf_comm) = reduced_th(1:nbf_comm)
          reduced_h(1:nbf_comm) = reduced_th(nbf_comm + 1:2 * nbf_comm)
          if (dg_frag%is_frag_root) then
            do io = 1, nbf
              ig_i = dg_frag%index_basis(io, ifrag, ispin)
              if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
              H_kin_diag_blocks(i_local)%val(io, jo, ispin) = reduced_t(io)
              H_diag_blocks(i_local)%val(io, jo, ispin) = reduced_h(io)
            end do
          end if

        end do  ! jo loop
        deallocate(partial_t, partial_h, reduced_t, reduced_h)
        deallocate(partial_th, reduced_th)
        deallocate(T_phi, H_phi)
        if (allocated(T_phi) .or. allocated(H_phi)) then
          write(*,*) "[FATAL] hamiltonian deallocate(T_phi,H_phi) failed: rank=", dg_frag%id, &
            " ifrag=", ifrag, " ispin=", ispin
          stop 1
        end if
          
        
      end do  ! ifrag loop
      
    end do  ! ispin loop
    
    call init_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, dg_frag%H_block_map, dg_frag%n_H_blocks)
    call init_matrix_blocks(dg_frag, dg_frag%H_mat_kinetic_blocks, dg_frag%H_block_map, dg_frag%n_H_blocks)
    call rebuild_local_h_block_ids(dg_frag)
    do i_diag = 1, n_local_diag
      iblk = find_matrix_block(dg_frag%H_block_map, H_diag_blocks(i_diag)%ifrag_row, H_diag_blocks(i_diag)%ifrag_col)
      if (iblk <= 0) cycle
      do ispin = 1, dg_frag%nspin
        nbf_diag = dg_frag%n_basis(H_diag_blocks(i_diag)%ifrag_row, ispin)
        if (nbf_diag <= 0) cycle
        dg_frag%H_mat_blocks(iblk)%val(1:nbf_diag, 1:nbf_diag, ispin) = H_diag_blocks(i_diag)%val(1:nbf_diag, 1:nbf_diag, ispin)
        dg_frag%H_mat_kinetic_blocks(iblk)%val(1:nbf_diag, 1:nbf_diag, ispin) = &
          H_kin_diag_blocks(i_diag)%val(1:nbf_diag, 1:nbf_diag, ispin)
      end do
    end do

    do ispin = 1, system%nspin
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        nbf_raw = dg_frag%n_basis(ifrag, ispin)
        nbf = min(nbf_raw, dg_frag%nstate_frag)
        if (nbf <= 0) cycle

        n_active_halo = 0
        n_basis_halo_max = 0
        do i_halo = 1, dg_frag%n_halo
          if (dg_frag%halo(i_halo)%ifrag_dst /= ifrag) cycle
          jfrag = dg_frag%halo(i_halo)%ifrag_src
          if (jfrag < 1) cycle
          n_basis_halo = min(dg_frag%n_basis(jfrag, ispin), dg_frag%nstate_frag)
          if (n_basis_halo <= 0) cycle
          if (.not. allocated(dg_frag%halo(i_halo)%buf_flow_recv)) cycle
          if (n_basis_halo > size(dg_frag%halo(i_halo)%buf_flow_recv, 4)) then
            n_basis_halo = size(dg_frag%halo(i_halo)%buf_flow_recv, 4)
            if (n_basis_halo <= 0) cycle
          end if
          if (.not. is_face_halo_neighbor(dg_frag%halo(i_halo))) cycle
          iblk = find_matrix_block(dg_frag%H_block_map, jfrag, ifrag)
          iblk_rev = find_matrix_block(dg_frag%H_block_map, ifrag, jfrag)
          if (iblk <= 0 .and. iblk_rev <= 0) cycle
          n_active_halo = n_active_halo + 1
          n_basis_halo_max = max(n_basis_halo_max, n_basis_halo)
        end do

        if (allocated(halo_active_list)) deallocate(halo_active_list)
        if (allocated(halo_active_nbf)) deallocate(halo_active_nbf)
        if (allocated(halo_active_iblk)) deallocate(halo_active_iblk)
        if (allocated(halo_active_iblk_rev)) deallocate(halo_active_iblk_rev)
        if (allocated(halo_active_face_axis)) deallocate(halo_active_face_axis)
        if (allocated(halo_active_face_dir)) deallocate(halo_active_face_dir)
        if (allocated(halo_active_face_idx_local)) deallocate(halo_active_face_idx_local)
        if (allocated(halo_active_face_idx_remote)) deallocate(halo_active_face_idx_remote)

        if (n_active_halo > 0) then
          allocate(halo_active_list(n_active_halo), halo_active_nbf(n_active_halo), halo_active_iblk(n_active_halo), halo_active_iblk_rev(n_active_halo))
          allocate(halo_active_face_axis(n_active_halo), halo_active_face_dir(n_active_halo), halo_active_face_idx_local(n_active_halo), halo_active_face_idx_remote(n_active_halo))
          i_halo_active = 0
          do i_halo = 1, dg_frag%n_halo
            if (dg_frag%halo(i_halo)%ifrag_dst /= ifrag) cycle
            jfrag = dg_frag%halo(i_halo)%ifrag_src
            if (jfrag < 1) cycle
            n_basis_halo = min(dg_frag%n_basis(jfrag, ispin), dg_frag%nstate_frag)
            if (n_basis_halo <= 0) cycle
            if (.not. allocated(dg_frag%halo(i_halo)%buf_flow_recv)) cycle
            if (n_basis_halo > size(dg_frag%halo(i_halo)%buf_flow_recv, 4)) then
              n_basis_halo = size(dg_frag%halo(i_halo)%buf_flow_recv, 4)
              if (n_basis_halo <= 0) cycle
            end if
            if (.not. is_face_halo_neighbor(dg_frag%halo(i_halo))) cycle
            iblk = find_matrix_block(dg_frag%H_block_map, jfrag, ifrag)
            iblk_rev = find_matrix_block(dg_frag%H_block_map, ifrag, jfrag)
            if (iblk <= 0 .and. iblk_rev <= 0) cycle
            i_halo_active = i_halo_active + 1
            halo_active_list(i_halo_active) = i_halo
            halo_active_nbf(i_halo_active) = n_basis_halo
            halo_active_iblk(i_halo_active) = iblk
            halo_active_iblk_rev(i_halo_active) = iblk_rev
            call get_halo_face_axis_dir(dg_frag%halo(i_halo), face_axis, face_dir)
            halo_active_face_axis(i_halo_active) = face_axis
            halo_active_face_dir(i_halo_active) = face_dir
            l = dg_frag%halo(i_halo)%length
            halo_active_face_idx_local(i_halo_active) = get_local_face_boundary_index(l(face_axis), face_dir)
            halo_active_face_idx_remote(i_halo_active) = get_remote_face_boundary_index(l(face_axis), face_dir)
          end do
        end if

        if (n_basis_halo_max > 0) then
          need_halo_alloc = .false.
          if (.not. allocated(halo_partial_t)) need_halo_alloc = .true.
          if (.not. allocated(halo_partial_h)) need_halo_alloc = .true.
          if (.not. allocated(halo_partial_m)) need_halo_alloc = .true.
          if (.not. allocated(halo_partial_a)) need_halo_alloc = .true.
          if (.not. allocated(halo_partial_b)) need_halo_alloc = .true.
          if (.not. allocated(halo_partial_c)) need_halo_alloc = .true.
          if (.not. allocated(halo_reduced_t)) need_halo_alloc = .true.
          if (.not. allocated(halo_reduced_h)) need_halo_alloc = .true.
          if (.not. allocated(halo_reduced_m)) need_halo_alloc = .true.
          if (.not. allocated(halo_reduced_a)) need_halo_alloc = .true.
          if (.not. allocated(halo_reduced_b)) need_halo_alloc = .true.
          if (.not. allocated(halo_reduced_c)) need_halo_alloc = .true.
          if (.not. allocated(halo_reduce_pair)) need_halo_alloc = .true.
          if (.not. allocated(halo_reduce_sum)) need_halo_alloc = .true.
          if (.not. need_halo_alloc) then
            if (size(halo_partial_t) < n_basis_halo_max) need_halo_alloc = .true.
            if (size(halo_partial_h) < n_basis_halo_max) need_halo_alloc = .true.
            if (size(halo_partial_m) < n_basis_halo_max) need_halo_alloc = .true.
            if (size(halo_partial_a) < n_basis_halo_max) need_halo_alloc = .true.
            if (size(halo_partial_b) < n_basis_halo_max) need_halo_alloc = .true.
            if (size(halo_partial_c) < n_basis_halo_max) need_halo_alloc = .true.
            if (size(halo_reduced_t) < n_basis_halo_max) need_halo_alloc = .true.
            if (size(halo_reduced_h) < n_basis_halo_max) need_halo_alloc = .true.
            if (size(halo_reduced_m) < n_basis_halo_max) need_halo_alloc = .true.
            if (size(halo_reduced_a) < n_basis_halo_max) need_halo_alloc = .true.
            if (size(halo_reduced_b) < n_basis_halo_max) need_halo_alloc = .true.
            if (size(halo_reduced_c) < n_basis_halo_max) need_halo_alloc = .true.
            if (size(halo_reduce_pair) < 6 * n_basis_halo_max) need_halo_alloc = .true.
            if (size(halo_reduce_sum) < 6 * n_basis_halo_max) need_halo_alloc = .true.
          end if
          if (need_halo_alloc) then
            if (allocated(halo_partial_t)) deallocate(halo_partial_t)
            if (allocated(halo_partial_h)) deallocate(halo_partial_h)
            if (allocated(halo_partial_m)) deallocate(halo_partial_m)
            if (allocated(halo_partial_a)) deallocate(halo_partial_a)
            if (allocated(halo_partial_b)) deallocate(halo_partial_b)
            if (allocated(halo_partial_c)) deallocate(halo_partial_c)
            if (allocated(halo_reduced_t)) deallocate(halo_reduced_t)
            if (allocated(halo_reduced_h)) deallocate(halo_reduced_h)
            if (allocated(halo_reduced_m)) deallocate(halo_reduced_m)
            if (allocated(halo_reduced_a)) deallocate(halo_reduced_a)
            if (allocated(halo_reduced_b)) deallocate(halo_reduced_b)
            if (allocated(halo_reduced_c)) deallocate(halo_reduced_c)
            if (allocated(halo_reduce_pair)) deallocate(halo_reduce_pair)
            if (allocated(halo_reduce_sum)) deallocate(halo_reduce_sum)
            allocate(halo_partial_t(n_basis_halo_max), halo_partial_h(n_basis_halo_max), halo_partial_m(n_basis_halo_max))
            allocate(halo_partial_a(n_basis_halo_max), halo_partial_b(n_basis_halo_max), halo_partial_c(n_basis_halo_max))
            allocate(halo_reduced_t(n_basis_halo_max), halo_reduced_h(n_basis_halo_max), halo_reduced_m(n_basis_halo_max))
            allocate(halo_reduced_a(n_basis_halo_max), halo_reduced_b(n_basis_halo_max), halo_reduced_c(n_basis_halo_max))
            allocate(halo_reduce_pair(6 * n_basis_halo_max), halo_reduce_sum(6 * n_basis_halo_max))
          end if
        end if

        do jo = 1, nbf
          do i_halo_active = 1, n_active_halo
            i_halo = halo_active_list(i_halo_active)
            jfrag = dg_frag%halo(i_halo)%ifrag_src
            n_basis_halo = halo_active_nbf(i_halo_active)
            iblk = halo_active_iblk(i_halo_active)
            iblk_rev = halo_active_iblk_rev(i_halo_active)
            face_axis = halo_active_face_axis(i_halo_active)
            face_dir = halo_active_face_dir(i_halo_active)
            face_idx_local = halo_active_face_idx_local(i_halo_active)
            face_idx_remote = halo_active_face_idx_remote(i_halo_active)
            l = dg_frag%halo(i_halo)%length

            ! Legacy off-diagonal coupling path.
            !
            ! This accumulates the SIPG face contribution for fragment-neighbor
            ! couplings. 
            !
            ! The remaining field dependence is carried by
            !   - the bulk A·∇ operator via momentum_blocks/momentum_mat
            !   - the A^2/2 term in the local potential path
            !   - the nonlocal projector path handled separately.
            !
            ! Only face-neighbor halos correspond to the surface integrals in
            ! Eq.(4). Edge/corner halos remain available for stencil support
            ! elsewhere, but they are excluded from the Hamiltonian block
            ! coupling path here.

            halo_partial_t(1:n_basis_halo) = 0.0d0
            halo_partial_h(1:n_basis_halo) = 0.0d0
            halo_partial_m(1:n_basis_halo) = 0.0d0
            halo_partial_a(1:n_basis_halo) = 0.0d0
            halo_partial_b(1:n_basis_halo) = 0.0d0
            halo_partial_c(1:n_basis_halo) = 0.0d0
            !$omp parallel do private(io,halo_integral_t,halo_integral_h,halo_integral_m,halo_integral_m_avgavg,halo_integral_m_avgjump,halo_integral_a_sum,halo_integral_b_sum,halo_integral_c_sum,ix_face,iy_face,iz_face,halo_send_idx,halo_recv_idx,phi_local_face,phi_remote_face,phi_face_avg,phi_face_jump,dphi_local_n,dphi_remote_n,dphi_flux_n,dphi_face_avg,dphi_face_jump,sipg_h_face_a,sipg_h_face_b,sipg_h_face_c,sipg_ab_sum,sipg_ab_cap,sipg_ab_scale,sipg_ab_ref,sipg_ab_norm) schedule(static)
            do io = 1, n_basis_halo
              halo_integral_t = 0.0d0
              halo_integral_h = 0.0d0
              halo_integral_m = 0.0d0
              halo_integral_a_sum = 0.0d0
              halo_integral_b_sum = 0.0d0
              halo_integral_c_sum = 0.0d0
              halo_integral_m_avgavg = 0.0d0
              halo_integral_m_avgjump = 0.0d0
              select case (face_axis)
              case (1)
                do iz_face = 1, l(3)
                  do iy_face = 1, l(2)
                    call get_halo_block_point_indices(dg_frag%halo(i_halo), face_idx_local, iy_face, iz_face, halo_send_idx, halo_recv_idx)
                    phi_local_face = face_trace_value_local(dg_frag, i_local, jo, face_axis, face_dir, &
                      halo_send_idx(1), halo_send_idx(2), halo_send_idx(3))
                    dphi_local_n = one_sided_face_derivative_local(dg_frag, i_local, jo, face_axis, face_dir, &
                      halo_send_idx(1), halo_send_idx(2), halo_send_idx(3))
                    phi_remote_face = face_trace_value_flow_halo(dg_frag%halo(i_halo)%buf_flow_recv, l, io, face_axis, &
                      face_idx_remote, iy_face, iz_face)
                    dphi_remote_n = face_derivative_value_flow_halo(dg_frag%halo(i_halo)%buf_flow_recv, l, io, face_axis, &
                      face_idx_remote, iy_face, iz_face, dg_frag%hgs)
                    ! SIPG off-diagonal face term using explicit average/jump roles:
                    !   - {{grad phi_local}} · [[phi_remote]]
                    !   - {{grad phi_remote}} · [[phi_local]]
                    !   + alpha [[phi_remote]] · [[phi_local]]
                    ! With basis support localized to each element, this reduces to
                    !   +0.5 * phi_remote * dphi_local_n
                    !   -0.5 * dphi_remote_n * phi_local
                    !   -alpha * phi_remote * phi_local
                    if (sipg_disable_b_term) then
                      ! A-only mode: always use the pure-A definition and bypass
                      ! weak/avg/flux mixed reconstructions.
                      sipg_h_face_a = 0.5d0 * phi_remote_face * dphi_local_n * face_area
                      sipg_h_face_b = 0.0d0
                    else if (sipg_use_weak_form) then
                      phi_face_avg = 0.5d0 * (phi_local_face + phi_remote_face)
                      phi_face_jump = phi_remote_face - phi_local_face
                      dphi_face_avg = 0.5d0 * (dphi_local_n + dphi_remote_n)
                      dphi_face_jump = dphi_local_n - dphi_remote_n
                      sipg_h_face_a = 0.5d0 * phi_face_avg * dphi_face_jump * face_area
                      sipg_h_face_b = 0.5d0 * phi_face_jump * dphi_face_avg * face_area
                    else if (sipg_use_flux_form) then
                      dphi_flux_n = 0.5d0 * (dphi_local_n + dphi_remote_n)
                      sipg_h_face_a = 0.5d0 * phi_remote_face * dphi_flux_n * face_area
                      sipg_h_face_b = -0.5d0 * phi_local_face * dphi_flux_n * face_area
                    else if (sipg_use_avg_trace) then
                      phi_face_avg = 0.5d0 * (phi_local_face + phi_remote_face)
                      sipg_h_face_a = 0.5d0 * phi_face_avg * dphi_local_n * face_area
                      sipg_h_face_b = -0.5d0 * dphi_remote_n * phi_face_avg * face_area
                    else
                      sipg_h_face_a = 0.5d0 * phi_remote_face * dphi_local_n * face_area
                      sipg_h_face_b = -0.5d0 * dphi_remote_n * phi_local_face * face_area
                    end if
                    sipg_h_face_c = -face_penalty_alpha * phi_remote_face * phi_local_face * face_area
                    if (sipg_ab_norm_max > 0.0d0) then
                      sipg_ab_norm = abs(sipg_h_face_a) + abs(sipg_h_face_b)
                      sipg_ab_ref = 0.5d0 * (abs(phi_remote_face) * abs(dphi_local_n) + &
                                            abs(phi_local_face) * abs(dphi_remote_n)) * face_area
                      sipg_ab_cap = sipg_ab_norm_max * max(sipg_ab_ref, 1.0d-30)
                      if (sipg_ab_norm > sipg_ab_cap) then
                        sipg_ab_scale = sipg_ab_cap / sipg_ab_norm
                        sipg_h_face_a = sipg_h_face_a * sipg_ab_scale
                        sipg_h_face_b = sipg_h_face_b * sipg_ab_scale
                      end if
                    end if
                    if (sipg_disable_c_term) sipg_h_face_c = 0.0d0
                    halo_integral_t = halo_integral_t + sipg_h_face_a + sipg_h_face_b + sipg_h_face_c
                    halo_integral_a_sum = halo_integral_a_sum + sipg_h_face_a
                    halo_integral_b_sum = halo_integral_b_sum + sipg_h_face_b
                    halo_integral_c_sum = halo_integral_c_sum + sipg_h_face_c
                    ! Intentionally split the face A-coupling into avg-avg and avg-jump halves
                    ! so the summed contribution keeps the original total strength.
                    halo_integral_m_avgavg = halo_integral_m_avgavg + 0.5d0 * dble(face_dir) * phi_remote_face * phi_local_face * face_area
                    halo_integral_m_avgjump = halo_integral_m_avgjump + 0.5d0 * dble(face_dir) * phi_remote_face * phi_local_face * face_area
                  end do
                end do
              case (2)
                do iz_face = 1, l(3)
                  do ix_face = 1, l(1)
                    call get_halo_block_point_indices(dg_frag%halo(i_halo), ix_face, face_idx_local, iz_face, halo_send_idx, halo_recv_idx)
                    phi_local_face = face_trace_value_local(dg_frag, i_local, jo, face_axis, face_dir, &
                      halo_send_idx(1), halo_send_idx(2), halo_send_idx(3))
                    dphi_local_n = one_sided_face_derivative_local(dg_frag, i_local, jo, face_axis, face_dir, &
                      halo_send_idx(1), halo_send_idx(2), halo_send_idx(3))
                    phi_remote_face = face_trace_value_flow_halo(dg_frag%halo(i_halo)%buf_flow_recv, l, io, face_axis, &
                      ix_face, face_idx_remote, iz_face)
                    dphi_remote_n = face_derivative_value_flow_halo(dg_frag%halo(i_halo)%buf_flow_recv, l, io, face_axis, &
                      ix_face, face_idx_remote, iz_face, dg_frag%hgs)
                    if (sipg_disable_b_term) then
                      sipg_h_face_a = 0.5d0 * phi_remote_face * dphi_local_n * face_area
                      sipg_h_face_b = 0.0d0
                    else if (sipg_use_weak_form) then
                      phi_face_avg = 0.5d0 * (phi_local_face + phi_remote_face)
                      phi_face_jump = phi_remote_face - phi_local_face
                      dphi_face_avg = 0.5d0 * (dphi_local_n + dphi_remote_n)
                      dphi_face_jump = dphi_local_n - dphi_remote_n
                      sipg_h_face_a = 0.5d0 * phi_face_avg * dphi_face_jump * face_area
                      sipg_h_face_b = 0.5d0 * phi_face_jump * dphi_face_avg * face_area
                    else if (sipg_use_flux_form) then
                      dphi_flux_n = 0.5d0 * (dphi_local_n + dphi_remote_n)
                      sipg_h_face_a = 0.5d0 * phi_remote_face * dphi_flux_n * face_area
                      sipg_h_face_b = -0.5d0 * phi_local_face * dphi_flux_n * face_area
                    else if (sipg_use_avg_trace) then
                      phi_face_avg = 0.5d0 * (phi_local_face + phi_remote_face)
                      sipg_h_face_a = 0.5d0 * phi_face_avg * dphi_local_n * face_area
                      sipg_h_face_b = -0.5d0 * dphi_remote_n * phi_face_avg * face_area
                    else
                      sipg_h_face_a = 0.5d0 * phi_remote_face * dphi_local_n * face_area
                      sipg_h_face_b = -0.5d0 * dphi_remote_n * phi_local_face * face_area
                    end if
                    sipg_h_face_c = -face_penalty_alpha * phi_remote_face * phi_local_face * face_area
                    if (sipg_ab_norm_max > 0.0d0) then
                      sipg_ab_norm = abs(sipg_h_face_a) + abs(sipg_h_face_b)
                      sipg_ab_ref = 0.5d0 * (abs(phi_remote_face) * abs(dphi_local_n) + &
                                            abs(phi_local_face) * abs(dphi_remote_n)) * face_area
                      sipg_ab_cap = sipg_ab_norm_max * max(sipg_ab_ref, 1.0d-30)
                      if (sipg_ab_norm > sipg_ab_cap) then
                        sipg_ab_scale = sipg_ab_cap / sipg_ab_norm
                        sipg_h_face_a = sipg_h_face_a * sipg_ab_scale
                        sipg_h_face_b = sipg_h_face_b * sipg_ab_scale
                      end if
                    end if
                    if (sipg_disable_c_term) sipg_h_face_c = 0.0d0
                    halo_integral_t = halo_integral_t + sipg_h_face_a + sipg_h_face_b + sipg_h_face_c
                    halo_integral_a_sum = halo_integral_a_sum + sipg_h_face_a
                    halo_integral_b_sum = halo_integral_b_sum + sipg_h_face_b
                    halo_integral_c_sum = halo_integral_c_sum + sipg_h_face_c
                    ! Intentionally split the face A-coupling into avg-avg and avg-jump halves
                    ! so the summed contribution keeps the original total strength.
                    halo_integral_m_avgavg = halo_integral_m_avgavg + 0.5d0 * dble(face_dir) * phi_remote_face * phi_local_face * face_area
                    halo_integral_m_avgjump = halo_integral_m_avgjump + 0.5d0 * dble(face_dir) * phi_remote_face * phi_local_face * face_area
                  end do
                end do
              case (3)
                do iy_face = 1, l(2)
                  do ix_face = 1, l(1)
                    call get_halo_block_point_indices(dg_frag%halo(i_halo), ix_face, iy_face, face_idx_local, halo_send_idx, halo_recv_idx)
                    phi_local_face = face_trace_value_local(dg_frag, i_local, jo, face_axis, face_dir, &
                      halo_send_idx(1), halo_send_idx(2), halo_send_idx(3))
                    dphi_local_n = one_sided_face_derivative_local(dg_frag, i_local, jo, face_axis, face_dir, &
                      halo_send_idx(1), halo_send_idx(2), halo_send_idx(3))
                    phi_remote_face = face_trace_value_flow_halo(dg_frag%halo(i_halo)%buf_flow_recv, l, io, face_axis, &
                      ix_face, iy_face, face_idx_remote)
                    dphi_remote_n = face_derivative_value_flow_halo(dg_frag%halo(i_halo)%buf_flow_recv, l, io, face_axis, &
                      ix_face, iy_face, face_idx_remote, dg_frag%hgs)
                    if (sipg_disable_b_term) then
                      sipg_h_face_a = 0.5d0 * phi_remote_face * dphi_local_n * face_area
                      sipg_h_face_b = 0.0d0
                    else if (sipg_use_weak_form) then
                      phi_face_avg = 0.5d0 * (phi_local_face + phi_remote_face)
                      phi_face_jump = phi_remote_face - phi_local_face
                      dphi_face_avg = 0.5d0 * (dphi_local_n + dphi_remote_n)
                      dphi_face_jump = dphi_local_n - dphi_remote_n
                      sipg_h_face_a = 0.5d0 * phi_face_avg * dphi_face_jump * face_area
                      sipg_h_face_b = 0.5d0 * phi_face_jump * dphi_face_avg * face_area
                    else if (sipg_use_flux_form) then
                      dphi_flux_n = 0.5d0 * (dphi_local_n + dphi_remote_n)
                      sipg_h_face_a = 0.5d0 * phi_remote_face * dphi_flux_n * face_area
                      sipg_h_face_b = -0.5d0 * phi_local_face * dphi_flux_n * face_area
                    else if (sipg_use_avg_trace) then
                      phi_face_avg = 0.5d0 * (phi_local_face + phi_remote_face)
                      sipg_h_face_a = 0.5d0 * phi_face_avg * dphi_local_n * face_area
                      sipg_h_face_b = -0.5d0 * dphi_remote_n * phi_face_avg * face_area
                    else
                      sipg_h_face_a = 0.5d0 * phi_remote_face * dphi_local_n * face_area
                      sipg_h_face_b = -0.5d0 * dphi_remote_n * phi_local_face * face_area
                    end if
                    sipg_h_face_c = -face_penalty_alpha * phi_remote_face * phi_local_face * face_area
                    if (sipg_ab_norm_max > 0.0d0) then
                      sipg_ab_norm = abs(sipg_h_face_a) + abs(sipg_h_face_b)
                      sipg_ab_ref = 0.5d0 * (abs(phi_remote_face) * abs(dphi_local_n) + &
                                            abs(phi_local_face) * abs(dphi_remote_n)) * face_area
                      sipg_ab_cap = sipg_ab_norm_max * max(sipg_ab_ref, 1.0d-30)
                      if (sipg_ab_norm > sipg_ab_cap) then
                        sipg_ab_scale = sipg_ab_cap / sipg_ab_norm
                        sipg_h_face_a = sipg_h_face_a * sipg_ab_scale
                        sipg_h_face_b = sipg_h_face_b * sipg_ab_scale
                      end if
                    end if
                    if (sipg_disable_c_term) sipg_h_face_c = 0.0d0
                    halo_integral_t = halo_integral_t + sipg_h_face_a + sipg_h_face_b + sipg_h_face_c
                    halo_integral_a_sum = halo_integral_a_sum + sipg_h_face_a
                    halo_integral_b_sum = halo_integral_b_sum + sipg_h_face_b
                    halo_integral_c_sum = halo_integral_c_sum + sipg_h_face_c
                    ! Intentionally split the face A-coupling into avg-avg and avg-jump halves
                    ! so the summed contribution keeps the original total strength.
                    halo_integral_m_avgavg = halo_integral_m_avgavg + 0.5d0 * dble(face_dir) * phi_remote_face * phi_local_face * face_area
                    halo_integral_m_avgjump = halo_integral_m_avgjump + 0.5d0 * dble(face_dir) * phi_remote_face * phi_local_face * face_area
                  end do
                end do
              end select
              ! dH/dA(face) = avg-avg term + avg-jump term on each interface.
              halo_integral_m = halo_integral_m_avgavg + halo_integral_m_avgjump
              ! Neighbor-face off-diagonal coupling is SIPG/kinetic only.
              ! Hartree/XC (Vh/Vxc) are volume-local terms and are not taken
              ! from neighbor blocks in this halo path.
              halo_integral_h = halo_integral_t
              halo_partial_t(io) = halo_integral_t
              halo_partial_h(io) = halo_integral_h
              halo_partial_m(io) = halo_integral_m
              halo_partial_a(io) = halo_integral_a_sum
              halo_partial_b(io) = halo_integral_b_sum
              halo_partial_c(io) = halo_integral_c_sum
            end do
            !$omp end parallel do

            halo_reduce_pair(1:n_basis_halo) = halo_partial_t(1:n_basis_halo)
            ! Keep H halo channel identical to kinetic SIPG channel by design.
            halo_reduce_pair(n_basis_halo + 1:2 * n_basis_halo) = halo_partial_t(1:n_basis_halo)
            halo_reduce_pair(2 * n_basis_halo + 1:3 * n_basis_halo) = halo_partial_m(1:n_basis_halo)
            halo_reduce_pair(3 * n_basis_halo + 1:4 * n_basis_halo) = halo_partial_a(1:n_basis_halo)
            halo_reduce_pair(4 * n_basis_halo + 1:5 * n_basis_halo) = halo_partial_b(1:n_basis_halo)
            halo_reduce_pair(5 * n_basis_halo + 1:6 * n_basis_halo) = halo_partial_c(1:n_basis_halo)
            call comm_summation(halo_reduce_pair, halo_reduce_sum, 6 * n_basis_halo, dg_frag%icomm_frag)
            halo_reduced_t(1:n_basis_halo) = halo_reduce_sum(1:n_basis_halo)
            halo_reduced_h(1:n_basis_halo) = halo_reduced_t(1:n_basis_halo)
            halo_reduced_m(1:n_basis_halo) = halo_reduce_sum(2 * n_basis_halo + 1:3 * n_basis_halo)
            halo_reduced_a(1:n_basis_halo) = halo_reduce_sum(3 * n_basis_halo + 1:4 * n_basis_halo)
            halo_reduced_b(1:n_basis_halo) = halo_reduce_sum(4 * n_basis_halo + 1:5 * n_basis_halo)
            halo_reduced_c(1:n_basis_halo) = halo_reduce_sum(5 * n_basis_halo + 1:6 * n_basis_halo)

            if (.not. dg_frag%is_frag_root) cycle

            if (iblk > 0 .and. iblk_rev > 0) then
              !$omp parallel do private(io,halo_integral_t_sum,halo_integral_h_sum,halo_integral_m_sum) schedule(static)
              do io = 1, n_basis_halo
                halo_integral_t_sum = halo_reduced_t(io)
                halo_integral_h_sum = halo_reduced_h(io)
                halo_integral_m_sum = halo_reduced_m(io)
                if (iblk <= size(dg_frag%H_mat_kinetic_blocks) .and. iblk <= size(dg_frag%H_mat_blocks) .and. &
                    iblk_rev <= size(dg_frag%H_mat_kinetic_blocks) .and. iblk_rev <= size(dg_frag%H_mat_blocks) .and. &
                    io <= size(dg_frag%H_mat_kinetic_blocks(iblk)%val, 1) .and. &
                    jo <= size(dg_frag%H_mat_kinetic_blocks(iblk)%val, 2) .and. &
                    io <= size(dg_frag%H_mat_blocks(iblk)%val, 1) .and. &
                    jo <= size(dg_frag%H_mat_blocks(iblk)%val, 2) .and. &
                    jo <= size(dg_frag%H_mat_kinetic_blocks(iblk_rev)%val, 1) .and. &
                    io <= size(dg_frag%H_mat_kinetic_blocks(iblk_rev)%val, 2) .and. &
                    jo <= size(dg_frag%H_mat_blocks(iblk_rev)%val, 1) .and. &
                    io <= size(dg_frag%H_mat_blocks(iblk_rev)%val, 2)) then
                  dg_frag%H_mat_kinetic_blocks(iblk)%val(io, jo, ispin) = &
                    dg_frag%H_mat_kinetic_blocks(iblk)%val(io, jo, ispin) + 0.5d0 * halo_integral_t_sum
                  dg_frag%H_mat_blocks(iblk)%val(io, jo, ispin) = &
                    dg_frag%H_mat_blocks(iblk)%val(io, jo, ispin) + 0.5d0 * halo_integral_h_sum
                  dg_frag%H_mat_kinetic_blocks(iblk_rev)%val(jo, io, ispin) = &
                    dg_frag%H_mat_kinetic_blocks(iblk_rev)%val(jo, io, ispin) + 0.5d0 * halo_integral_t_sum
                  dg_frag%H_mat_blocks(iblk_rev)%val(jo, io, ispin) = &
                    dg_frag%H_mat_blocks(iblk_rev)%val(jo, io, ispin) + 0.5d0 * halo_integral_h_sum
                end if
                if (allocated(dg_frag%momentum_blocks) .and. face_axis >= 1 .and. face_axis <= 3) then
                  if (iblk <= size(dg_frag%momentum_blocks) .and. iblk_rev <= size(dg_frag%momentum_blocks) .and. &
                      io <= size(dg_frag%momentum_blocks(iblk)%val, 2) .and. &
                      jo <= size(dg_frag%momentum_blocks(iblk)%val, 3) .and. &
                      jo <= size(dg_frag%momentum_blocks(iblk_rev)%val, 2) .and. &
                      io <= size(dg_frag%momentum_blocks(iblk_rev)%val, 3)) then
                    dg_frag%momentum_blocks(iblk)%val(face_axis, io, jo, ispin) = &
                      dg_frag%momentum_blocks(iblk)%val(face_axis, io, jo, ispin) + 0.5d0 * halo_integral_m_sum
                    dg_frag%momentum_blocks(iblk_rev)%val(face_axis, jo, io, ispin) = &
                      dg_frag%momentum_blocks(iblk_rev)%val(face_axis, jo, io, ispin) - 0.5d0 * halo_integral_m_sum
                  end if
                end if
              end do
              !$omp end parallel do
            else if (iblk > 0) then
              !$omp parallel do private(io,halo_integral_t_sum,halo_integral_h_sum,halo_integral_m_sum) schedule(static)
              do io = 1, n_basis_halo
                halo_integral_t_sum = halo_reduced_t(io)
                halo_integral_h_sum = halo_reduced_h(io)
                halo_integral_m_sum = halo_reduced_m(io)
                if (iblk <= size(dg_frag%H_mat_kinetic_blocks) .and. iblk <= size(dg_frag%H_mat_blocks) .and. &
                    io <= size(dg_frag%H_mat_kinetic_blocks(iblk)%val, 1) .and. &
                    jo <= size(dg_frag%H_mat_kinetic_blocks(iblk)%val, 2) .and. &
                    io <= size(dg_frag%H_mat_blocks(iblk)%val, 1) .and. &
                    jo <= size(dg_frag%H_mat_blocks(iblk)%val, 2)) then
                  dg_frag%H_mat_kinetic_blocks(iblk)%val(io, jo, ispin) = &
                    dg_frag%H_mat_kinetic_blocks(iblk)%val(io, jo, ispin) + 0.5d0 * halo_integral_t_sum
                  dg_frag%H_mat_blocks(iblk)%val(io, jo, ispin) = &
                    dg_frag%H_mat_blocks(iblk)%val(io, jo, ispin) + 0.5d0 * halo_integral_h_sum
                end if
                if (allocated(dg_frag%momentum_blocks) .and. face_axis >= 1 .and. face_axis <= 3) then
                  if (iblk <= size(dg_frag%momentum_blocks) .and. &
                      io <= size(dg_frag%momentum_blocks(iblk)%val, 2) .and. &
                      jo <= size(dg_frag%momentum_blocks(iblk)%val, 3)) then
                    dg_frag%momentum_blocks(iblk)%val(face_axis, io, jo, ispin) = &
                      dg_frag%momentum_blocks(iblk)%val(face_axis, io, jo, ispin) + 0.5d0 * halo_integral_m_sum
                  end if
                end if
              end do
              !$omp end parallel do
            else
              !$omp parallel do private(io,halo_integral_t_sum,halo_integral_h_sum,halo_integral_m_sum) schedule(static)
              do io = 1, n_basis_halo
                halo_integral_t_sum = halo_reduced_t(io)
                halo_integral_h_sum = halo_reduced_h(io)
                halo_integral_m_sum = halo_reduced_m(io)
                if (iblk_rev > 0 .and. iblk_rev <= size(dg_frag%H_mat_kinetic_blocks) .and. iblk_rev <= size(dg_frag%H_mat_blocks) .and. &
                    jo <= size(dg_frag%H_mat_kinetic_blocks(iblk_rev)%val, 1) .and. &
                    io <= size(dg_frag%H_mat_kinetic_blocks(iblk_rev)%val, 2) .and. &
                    jo <= size(dg_frag%H_mat_blocks(iblk_rev)%val, 1) .and. &
                    io <= size(dg_frag%H_mat_blocks(iblk_rev)%val, 2)) then
                  dg_frag%H_mat_kinetic_blocks(iblk_rev)%val(jo, io, ispin) = &
                    dg_frag%H_mat_kinetic_blocks(iblk_rev)%val(jo, io, ispin) + 0.5d0 * halo_integral_t_sum
                  dg_frag%H_mat_blocks(iblk_rev)%val(jo, io, ispin) = &
                    dg_frag%H_mat_blocks(iblk_rev)%val(jo, io, ispin) + 0.5d0 * halo_integral_h_sum
                end if
                if (allocated(dg_frag%momentum_blocks) .and. face_axis >= 1 .and. face_axis <= 3) then
                  if (iblk_rev > 0 .and. iblk_rev <= size(dg_frag%momentum_blocks) .and. &
                      jo <= size(dg_frag%momentum_blocks(iblk_rev)%val, 2) .and. &
                      io <= size(dg_frag%momentum_blocks(iblk_rev)%val, 3)) then
                    dg_frag%momentum_blocks(iblk_rev)%val(face_axis, jo, io, ispin) = &
                      dg_frag%momentum_blocks(iblk_rev)%val(face_axis, jo, io, ispin) - 0.5d0 * halo_integral_m_sum
                  end if
                end if
              end do
              !$omp end parallel do
            end if
          end do
        end do

        if (allocated(halo_active_list)) deallocate(halo_active_list)
        if (allocated(halo_active_nbf)) deallocate(halo_active_nbf)
        if (allocated(halo_active_iblk)) deallocate(halo_active_iblk)
        if (allocated(halo_active_iblk_rev)) deallocate(halo_active_iblk_rev)
        if (allocated(halo_active_face_axis)) deallocate(halo_active_face_axis)
        if (allocated(halo_active_face_dir)) deallocate(halo_active_face_dir)
        if (allocated(halo_active_face_idx_local)) deallocate(halo_active_face_idx_local)
        if (allocated(halo_active_face_idx_remote)) deallocate(halo_active_face_idx_remote)
      end do
    end do

    do i_diag = 1, n_local_diag
      if (allocated(H_diag_blocks(i_diag)%val)) deallocate(H_diag_blocks(i_diag)%val)
      if (allocated(H_kin_diag_blocks(i_diag)%val)) deallocate(H_kin_diag_blocks(i_diag)%val)
    end do
    if (allocated(halo_partial_t)) deallocate(halo_partial_t, halo_partial_h, halo_partial_m, halo_reduced_t, halo_reduced_h, halo_reduced_m)
    if (allocated(halo_partial_a)) deallocate(halo_partial_a, halo_partial_b, halo_partial_c, halo_reduced_a, halo_reduced_b, halo_reduced_c)
    if (allocated(halo_reduce_pair)) deallocate(halo_reduce_pair, halo_reduce_sum)
    if (allocated(halo_t_point_buf)) deallocate(halo_t_point_buf, halo_h_point_buf)
    if (allocated(H_diag_blocks)) deallocate(H_diag_blocks)
    if (allocated(H_kin_diag_blocks)) deallocate(H_kin_diag_blocks)
    ! CRITICAL: MPI aggregation of Hamiltonian matrix
    ! Each rank computed elements only for its assigned fragments.
    ! Reduce one fragment block at a time to avoid a single dense global allreduce.
    call reduce_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, "hmat", dg_frag%icomm)
    call reduce_matrix_blocks(dg_frag, dg_frag%H_mat_kinetic_blocks, "hmat-kinetic", dg_frag%icomm)
    call symmetrize_real_matrix_blocks(dg_frag, dg_frag%H_mat_blocks)
    call symmetrize_real_matrix_blocks(dg_frag, dg_frag%H_mat_kinetic_blocks)
    if (allocated(dg_frag%H_mat)) deallocate(dg_frag%H_mat)
    if (allocated(dg_frag%H_mat_kinetic)) deallocate(dg_frag%H_mat_kinetic)

    if (comm_is_root(dg_frag%id)) then
      write(*,*) "        Kinetic and potential terms computed"
    end if
    
    ! Step 3: Non-local pseudopotential is handled in time evolution
    ! with vector potential A(t), so it is not added to H_mat here.
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "  [3/3] Non-local PP handled in time evolution (A-dependent)"
    end if

    ! DG+PW must start from the orthogonalized mixed basis. Do not allow
    ! raw fragment/PW runtime propagation to begin from an uninitialized
    ! mixed representation.
    if (dg_frag%use_plane_wave_basis .and. dg_frag%n_plane_waves > 0) then
      Ac_zero(:) = 0.0d0
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "  [init] Diagonalizing mixed basis at startup (A=0)"
      end if
      call diagonalize_mixed_basis(dg_frag, system, Vh, Vxc, Vpsl, Ac_zero)
      if (.not. dg_frag%mixed_basis_ready) then
        stop "DG+PW startup requires mixed_basis_ready after diagonalize_mixed_basis"
      end if
      if (allocated(dg_frag%coef_new)) dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
    else if (.not. dg_frag%puredg_lowdin_applied) then
      if (enable_startup_lowdin) then
        call apply_startup_lowdin_puredg(dg_frag)
        if (enable_startup_stationary_projection) then
          call apply_startup_stationary_projection_puredg(dg_frag, mg, ppg, system)
        end if
      else if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') "[INFO] PureDG startup Lowdin/stationary projection disabled (set SALMON_DG_STARTUP_LOWDIN=1 to enable)"
        write(*,'(1x,a)') "[INFO] To additionally enable stationary projection, set SALMON_DG_STARTUP_STATIONARY_PROJECTION=1"
      end if
      if (allocated(dg_frag%coef_new)) dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
      dg_frag%puredg_lowdin_applied = .true.
    end if

    ! Initialize field-free reference Hamiltonian for adaptive-basis metric.
    if (allocated(dg_frag%H_mat_old)) then
      dg_frag%H_mat_old(:, :, :) = (0.0d0, 0.0d0)
      n_metric = min(dg_frag%nstate_frag, size(dg_frag%H_mat_old, 1), size(dg_frag%H_mat_old, 2))
      if (n_metric > 0) then
        allocate(H_metric_ref(n_metric, n_metric))
      end if
      do ispin = 1, min(dg_frag%nspin, size(dg_frag%H_mat_old,3))
        if (n_metric <= 0) cycle
        H_metric_ref(:, :) = (0.0d0, 0.0d0)
        call copy_matrix_blocks_metric_to_complex_dense(dg_frag, dg_frag%H_mat_blocks, ispin, n_metric, H_metric_ref)
        dg_frag%H_mat_old(1:n_metric, 1:n_metric, ispin) = H_metric_ref(:, :)
      end do
      if (allocated(H_metric_ref)) deallocate(H_metric_ref)
    end if
    
    deallocate(V_total)
    if (allocated(V_total)) then
      write(*,*) "[FATAL] V_total still allocated before return: rank=", dg_frag%id
      stop 1
    end if
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "=== Hamiltonian Matrix Ready ==="
      write(*,*)
    end if
    
  end subroutine calculate_hamiltonian_matrix

  !=======================================================================
  ! One-time startup Lowdin orthonormalization for PureDG coefficients
  !   C <- C (C^\dagger S C)^(-1/2)
  !=======================================================================
  subroutine apply_startup_lowdin_puredg(dg_frag)
    use communication, only: comm_is_root, comm_summation
    use rt_dg_fragment_ops, only: apply_overlap_operator, gather_full_coef_view, zero_nonowned_coefficients
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    integer :: ispin, n_frag, n_tot, nst, j, io
    integer :: info, lwork
    real(8), parameter :: eps_eval = 1.0d-12
    real(8), allocatable :: eval(:), rwork(:)
    complex(8), allocatable :: C(:,:), SC(:,:), M(:,:), M_sum(:,:), U(:,:), M_invhalf(:,:)
    complex(8), allocatable :: work(:), coef_frag_all(:,:), coef_pw_all(:,:)
    external :: zheev

    if (dg_frag%use_plane_wave_basis) return
    if (dg_frag%puredg_lowdin_applied) return

    n_frag = dg_frag%n_mat_max
    n_tot = n_frag
    nst = min(dg_frag%nstate_tot, n_tot)
    if (n_tot <= 0 .or. nst <= 0) then
      dg_frag%puredg_lowdin_applied = .true.
      return
    end if

    do ispin = 1, dg_frag%nspin
      call gather_full_coef_view(dg_frag, ispin, n_frag, dg_frag%nstate_tot, coef_frag_all, coef_pw_all)

      allocate(C(n_tot, nst), SC(n_tot, nst), M(nst, nst), M_sum(nst, nst), U(nst, nst), M_invhalf(nst, nst))
      allocate(eval(nst), rwork(max(1, 3 * nst - 2)))

      C(:, :) = coef_frag_all(1:n_frag, 1:nst)
      do io = 1, nst
        call apply_overlap_operator(dg_frag, ispin, C(:, io), SC(:, io), .true.)
      end do
      M(:, :) = matmul(conjg(transpose(C)), SC)
      call comm_summation(M, M_sum, nst * nst, dg_frag%icomm)
      M(:, :) = 0.5d0 * (M_sum(:, :) + conjg(transpose(M_sum(:, :))))

      U(:, :) = M(:, :)
      lwork = -1
      allocate(work(1))
      call ZHEEV('V', 'U', nst, U, nst, eval, work, lwork, rwork, info)
      lwork = int(real(work(1), kind=8)) + 1
      deallocate(work)
      allocate(work(lwork))
      call ZHEEV('V', 'U', nst, U, nst, eval, work, lwork, rwork, info)
      if (info /= 0) then
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a,i0)') "[WARN] PureDG startup Lowdin failed: ispin=", ispin, " info=", info
        end if
        deallocate(C, SC, M, M_sum, U, M_invhalf, eval, rwork, work)
        if (allocated(coef_frag_all)) deallocate(coef_frag_all)
        if (allocated(coef_pw_all)) deallocate(coef_pw_all)
        cycle
      end if

      M_invhalf(:, :) = (0.0d0, 0.0d0)
      do j = 1, nst
        if (eval(j) > eps_eval) then
          M_invhalf(:, j) = U(:, j) / sqrt(eval(j))
        end if
      end do
      M_invhalf(:, :) = matmul(M_invhalf, conjg(transpose(U)))
      C(:, :) = matmul(C, M_invhalf)

      dg_frag%coef(1:n_frag, 1:nst, ispin) = C(:, :)

      deallocate(C, SC, M, M_sum, U, M_invhalf, eval, rwork, work)
      if (allocated(coef_frag_all)) deallocate(coef_frag_all)
      if (allocated(coef_pw_all)) deallocate(coef_pw_all)
    end do

    call zero_nonowned_coefficients(dg_frag)
    if (allocated(dg_frag%coef_new)) dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
    dg_frag%puredg_lowdin_applied = .true.

    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') "[INFO] Applied one-time PureDG startup Lowdin orthonormalization"
    end if
  end subroutine apply_startup_lowdin_puredg

  !=======================================================================
  ! One-time startup stationary projection for PureDG coefficients
  !   1) Build dense H and S in DG coefficient space
  !   2) Solve generalized EVP H c = e S c via Lowdin transform
  !   3) Replace startup coefficients with S-orthonormal eigenvectors
  !=======================================================================
  subroutine apply_startup_stationary_projection_puredg(dg_frag, mg, ppg, system)
    use communication, only: comm_is_root, comm_summation
    use structures, only: s_rgrid, s_pp_grid, s_dft_system
    use rt_dg_fragment_ops, only: copy_matrix_blocks_metric_to_complex_dense, &
                                  copy_overlap_operator_to_dense, zero_nonowned_coefficients, &
                                  ensure_nonlocal_pp_matrix_A, copy_complex_matrix_blocks_metric_to_dense
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    type(s_dft_system),     intent(in)    :: system

    integer :: ispin, n_frag, nst, nocc, info, lwork, j, nev, env_status
    real(8), parameter :: eps_s = 1.0d-12
    real(8) :: eval_s_min, eval_s_max
    real(8) :: resid_max, resid_rms, resid_rel_max, resid_rel_rms
    real(8) :: resid_norm, hc_norm, s_norm
    logical :: allow_buffered_startup_stationary_projection
    logical :: buffered_wf_coef_mode
    character(len=32) :: env_allow_buffered_startup_stationary_projection
    character(len=32) :: env_buffered_coef_source
    real(8), allocatable :: eval_s(:), eval_h(:), rwork(:)
    complex(8), allocatable :: h_dense(:,:), h_dense_sum(:,:), s_dense(:,:), s_dense_sum(:,:), x_lowdin(:,:), h_ortho(:,:), u_h(:,:), c_eig(:,:)
    complex(8), allocatable :: hc(:), sc(:)
    complex(8), allocatable :: work(:)
    external :: zheev

    if (dg_frag%use_plane_wave_basis) return

    allow_buffered_startup_stationary_projection = .false.
    env_allow_buffered_startup_stationary_projection = ''
    call get_environment_variable('SALMON_DG_ALLOW_BUFFERED_STARTUP_STATIONARY_PROJECTION', &
      env_allow_buffered_startup_stationary_projection, status=env_status)
    if (env_status == 0) then
      select case (trim(adjustl(env_allow_buffered_startup_stationary_projection)))
      case ('1', 'true', 'TRUE', 'on', 'ON', 'yes', 'YES')
        allow_buffered_startup_stationary_projection = .true.
      end select
    end if

    buffered_wf_coef_mode = .false.
    env_buffered_coef_source = ''
    call get_environment_variable('SALMON_DG_BUFFERED_COEF_SOURCE', env_buffered_coef_source, status=env_status)
    if (env_status == 0) then
      select case (trim(adjustl(env_buffered_coef_source)))
      case ('wf', 'WF', 'unbuffered', 'UNBUFFERED', 'plain', 'PLAIN', '0')
        buffered_wf_coef_mode = .true.
      end select
    end if

    if (buffered_wf_coef_mode .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') "[INFO] PureDG startup stationary projection enabled in buffered+wf mode"
    end if

    if (dg_frag%use_buffered_basis .and. .not. allow_buffered_startup_stationary_projection) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') "[INFO] PureDG startup stationary projection skipped in buffered-basis mode (memory guard)"
        write(*,'(1x,a)') "[INFO] Set SALMON_DG_ALLOW_BUFFERED_STARTUP_STATIONARY_PROJECTION=1 to override"
      end if
      return
    end if
    if (.not. allocated(dg_frag%H_mat_blocks)) return
    if (.not. allocated(dg_frag%esp)) return

    n_frag = dg_frag%n_mat_max
    nst = min(dg_frag%nstate_tot, n_frag, size(dg_frag%esp, 1))
    if (n_frag <= 0 .or. nst <= 0) return
    if (nst < n_frag) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,i0,a,i0,a)') "[INFO] PureDG startup stationary projection skipped: reduced state mode is active (nst=", &
          nst, ", n_frag=", n_frag, ")"
        write(*,'(1x,a)') "[INFO] Use SALMON_DG_BUFFERED_FULL_STATE=1 if full-state startup stationary projection is required"
      end if
      return
    end if

    do ispin = 1, dg_frag%nspin
        nocc = nst
        if (allocated(dg_frag%nocc_spin) .and. size(dg_frag%nocc_spin) >= ispin) then
          nocc = min(nst, max(0, dg_frag%nocc_spin(ispin)))
        end if
        if (nocc <= 0) cycle

        nev = nocc
        allocate(h_dense(nst, nst), h_dense_sum(nst, nst), s_dense(nst, nst), s_dense_sum(nst, nst), x_lowdin(nst, nst), h_ortho(nst, nst), u_h(nst, nst), c_eig(nst, nst))
        allocate(hc(nst), sc(nst))
        allocate(eval_s(nst), eval_h(nst), rwork(max(1, 3 * nst - 2)))

        h_dense(:, :) = (0.0d0, 0.0d0)
        s_dense(:, :) = (0.0d0, 0.0d0)
        call copy_matrix_blocks_metric_to_complex_dense(dg_frag, dg_frag%H_mat_blocks, ispin, nst, h_dense)
        call copy_overlap_operator_to_dense(dg_frag, ispin, .true., s_dense)

        ! Add non-local pseudopotential contribution to H at A=0
        block
          real(8) :: Ac_zero(3)
          complex(8), allocatable :: h_nl_dense(:,:)
          Ac_zero(:) = 0.0d0
          call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_zero, .false.)
          if (allocated(dg_frag%H_nl_blocks)) then
            allocate(h_nl_dense(nst, nst))
            h_nl_dense(:, :) = (0.0d0, 0.0d0)
            call copy_complex_matrix_blocks_metric_to_dense(dg_frag, dg_frag%H_nl_blocks, ispin, nst, h_nl_dense)
            h_dense(:, :) = h_dense(:, :) + h_nl_dense(:, :)
            if (comm_is_root(dg_frag%id)) then
              write(*,'(1x,a,i0,a,1pe12.4)') &
                "[INFO] PureDG startup: V_NL added to H for ispin=", ispin, &
                " ||V_NL||_max=", maxval(abs(h_nl_dense(:, :)))
            end if
            deallocate(h_nl_dense)
          end if
        end block

        call comm_summation(h_dense, h_dense_sum, nst * nst, dg_frag%icomm)
        call comm_summation(s_dense, s_dense_sum, nst * nst, dg_frag%icomm)
        h_dense(:, :) = 0.5d0 * (h_dense_sum(:, :) + conjg(transpose(h_dense_sum(:, :))))
        s_dense(:, :) = 0.5d0 * (s_dense_sum(:, :) + conjg(transpose(s_dense_sum(:, :))))

        ! Lowdin transform matrix X = U * diag(s^{-1/2})
        x_lowdin(:, :) = s_dense(:, :)
        lwork = -1
        allocate(work(1))
        call ZHEEV('V', 'U', nst, x_lowdin, nst, eval_s, work, lwork, rwork, info)
        lwork = int(real(work(1), kind=8)) + 1
        deallocate(work)
        allocate(work(lwork))
        call ZHEEV('V', 'U', nst, x_lowdin, nst, eval_s, work, lwork, rwork, info)
        if (info /= 0) then
          if (comm_is_root(dg_frag%id)) then
            write(*,'(1x,a,i0,a,i0)') "[WARN] PureDG startup stationary projection skipped (S diagonalization failed): ispin=", ispin, " info=", info
          end if
          deallocate(h_dense, h_dense_sum, s_dense, s_dense_sum, x_lowdin, h_ortho, u_h, c_eig, eval_s, eval_h, rwork, work)
          cycle
        end if

        eval_s_min = huge(1.0d0)
        eval_s_max = 0.0d0
        do j = 1, nst
          eval_s_max = max(eval_s_max, eval_s(j))
          if (eval_s(j) > eps_s) then
            x_lowdin(:, j) = x_lowdin(:, j) / sqrt(eval_s(j))
            eval_s_min = min(eval_s_min, eval_s(j))
          else
            x_lowdin(:, j) = (0.0d0, 0.0d0)
          end if
        end do

        ! H_ortho = X^H H X and diagonalize
        h_ortho(:, :) = matmul(conjg(transpose(x_lowdin)), matmul(h_dense, x_lowdin))
        u_h(:, :) = h_ortho(:, :)
        deallocate(work)
        lwork = -1
        allocate(work(1))
        call ZHEEV('V', 'U', nst, u_h, nst, eval_h, work, lwork, rwork, info)
        lwork = int(real(work(1), kind=8)) + 1
        deallocate(work)
        allocate(work(lwork))
        call ZHEEV('V', 'U', nst, u_h, nst, eval_h, work, lwork, rwork, info)
        if (info /= 0) then
          if (comm_is_root(dg_frag%id)) then
            write(*,'(1x,a,i0,a,i0)') "[WARN] PureDG startup stationary projection skipped (H_ortho diagonalization failed): ispin=", ispin, " info=", info
          end if
          deallocate(h_dense, h_dense_sum, s_dense, s_dense_sum, x_lowdin, h_ortho, u_h, c_eig, eval_s, eval_h, rwork, work)
          cycle
        end if

        c_eig(:, :) = matmul(x_lowdin, u_h)
        block
          integer :: info_occ, lwork_occ
          real(8), allocatable :: eval_occ(:), rwork_occ(:)
          complex(8), allocatable :: gram_occ(:,:), gram_vecs(:,:), work_occ(:)

          allocate(eval_occ(nev), rwork_occ(max(1, 3 * nev - 2)))
          allocate(gram_occ(nev, nev), gram_vecs(nev, nev))
          gram_occ(:, :) = matmul(conjg(transpose(c_eig(:, 1:nev))), matmul(s_dense, c_eig(:, 1:nev)))
          gram_occ(:, :) = 0.5d0 * (gram_occ(:, :) + conjg(transpose(gram_occ(:, :))))
          gram_vecs(:, :) = gram_occ(:, :)

          lwork_occ = -1
          allocate(work_occ(1))
          call ZHEEV('V', 'U', nev, gram_vecs, nev, eval_occ, work_occ, lwork_occ, rwork_occ, info_occ)
          lwork_occ = int(real(work_occ(1), kind=8)) + 1
          deallocate(work_occ)
          allocate(work_occ(lwork_occ))
          call ZHEEV('V', 'U', nev, gram_vecs, nev, eval_occ, work_occ, lwork_occ, rwork_occ, info_occ)
          if (info_occ == 0) then
            gram_occ(:, :) = (0.0d0, 0.0d0)
            do j = 1, nev
              if (eval_occ(j) > eps_s) gram_occ(:, j) = gram_vecs(:, j) / sqrt(eval_occ(j))
            end do
            gram_occ(:, :) = matmul(gram_occ(:, :), conjg(transpose(gram_vecs(:, :))))
            c_eig(:, 1:nev) = matmul(c_eig(:, 1:nev), gram_occ(:, :))
          end if
          deallocate(eval_occ, rwork_occ, gram_occ, gram_vecs, work_occ)
        end block
        do j = 1, nev
          sc(:) = matmul(s_dense, c_eig(:, j))
          s_norm = real(sum(conjg(c_eig(:, j)) * sc(:)), kind=8)
          if (s_norm > eps_s) c_eig(:, j) = c_eig(:, j) / sqrt(s_norm)
        end do

        ! Quantify generalized eigen residual on occupied startup manifold:
        !   r = H c - S c e
        resid_max = 0.0d0
        resid_rms = 0.0d0
        resid_rel_max = 0.0d0
        resid_rel_rms = 0.0d0
        do j = 1, nev
          hc(:) = matmul(h_dense, c_eig(:, j))
          sc(:) = matmul(s_dense, c_eig(:, j))
          resid_norm = sqrt(sum(abs(hc(:) - sc(:) * eval_h(j))**2))
          hc_norm = sqrt(sum(abs(hc(:))**2))
          resid_max = max(resid_max, resid_norm)
          resid_rms = resid_rms + resid_norm * resid_norm
          resid_rel_max = max(resid_rel_max, resid_norm / max(hc_norm, 1.0d-30))
          resid_rel_rms = resid_rel_rms + (resid_norm / max(hc_norm, 1.0d-30))**2
        end do
        resid_rms = sqrt(resid_rms / max(1, nev))
        resid_rel_rms = sqrt(resid_rel_rms / max(1, nev))

        ! Keep startup projection on occupied manifold only.
        ! This avoids reshuffling non-occupied channels at t=0.
        dg_frag%coef(1:n_frag, 1:nev, ispin) = c_eig(1:n_frag, 1:nev)
        dg_frag%esp(1:nev, ispin) = eval_h(1:nev)

        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a,i0,a,1pe12.4,a,1pe12.4)') "[INFO] PureDG startup stationary projection applied: ispin=", ispin, &
            " nocc=", nev, " S-eig(min,max)=", eval_s_min, ", ", eval_s_max
          write(*,'(1x,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
            "[INFO] PureDG startup projection residual: ispin=", ispin, &
            " abs(max,rms)=", resid_max, ", ", resid_rms, &
            " rel(max,rms)=", resid_rel_max, ", ", resid_rel_rms
        end if

        deallocate(h_dense, h_dense_sum, s_dense, s_dense_sum, x_lowdin, h_ortho, u_h, c_eig, hc, sc, eval_s, eval_h, rwork, work)
    end do

    call zero_nonowned_coefficients(dg_frag)
    if (allocated(dg_frag%coef_new)) dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
  end subroutine apply_startup_stationary_projection_puredg

  subroutine reduce_matrix_fragment_blocks(dg_frag, mat, label, icomm_reduce)
    use communication, only: comm_is_root, comm_summation, comm_get_max
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    real(8), intent(inout) :: mat(:,:,:)
    character(*), intent(in) :: label
    integer, intent(in) :: icomm_reduce
    integer, parameter :: reduce_chunk_size = 262144
    real(8), allocatable :: send_block(:), recv_block(:)
    integer :: ifrag_row, ifrag_col, ispin, ii, jj, ig_i, ig_j
    integer :: idx_ii, idx_jj, valid_row_count, valid_col_count
    integer :: nrow, ncol, block_size, max_block_size, total_active_size
    integer :: total_active_min, total_active_max, max_block_size_global
    integer :: chunk_begin, chunk_count, offset_flat
    integer, allocatable :: row_gid(:), col_gid(:), valid_row_ids(:), valid_col_ids(:)
    logical :: enable_reduce_trace
    character(16) :: env_reduce_trace
    integer :: env_status

    enable_reduce_trace = .false.
    env_reduce_trace = ''
    call get_environment_variable('SALMON_DG_HMAT_REDUCE_TRACE', env_reduce_trace, status=env_status)
    if (env_status == 0) then
      select case(trim(adjustl(env_reduce_trace)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enable_reduce_trace = .true.
      end select
    end if

    max_block_size = 0
    total_active_size = 0
    do ispin = 1, dg_frag%nspin
      do ifrag_col = 1, dg_frag%n_frag
        ncol = dg_frag%n_basis(ifrag_col, ispin)
        if (ncol <= 0) cycle
        do ifrag_row = 1, dg_frag%n_frag
          nrow = dg_frag%n_basis(ifrag_row, ispin)
          if (nrow <= 0) cycle
          block_size = nrow * ncol
          max_block_size = max(max_block_size, block_size)
          total_active_size = total_active_size + block_size
        end do
      end do
    end do

    max_block_size_global = max_block_size
    call comm_get_max(max_block_size_global, icomm_reduce)
    total_active_max = total_active_size
    call comm_get_max(total_active_max, icomm_reduce)
    total_active_min = -total_active_size
    call comm_get_max(total_active_min, icomm_reduce)
    total_active_min = -total_active_min

    if (total_active_min /= total_active_max) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] Hamiltonian block size mismatch: label=", &
        trim(label), " rank=", dg_frag%id, " local=", total_active_size, &
        " min=", total_active_min, " max=", total_active_max
      flush(6)
      stop 1
    end if

    if (enable_reduce_trace .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "        hamiltonian block reduce begin: label=", trim(label), &
        " total_active=", total_active_size, " max_block=", max_block_size_global, &
        " chunk_size=", reduce_chunk_size
      flush(6)
    end if

    if (max_block_size_global <= 0) return
    allocate(send_block(max_block_size_global), recv_block(max_block_size_global))
    allocate(row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1)))
    allocate(valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1)))

    do ispin = 1, dg_frag%nspin
      do ifrag_col = 1, dg_frag%n_frag
        ncol = dg_frag%n_basis(ifrag_col, ispin)
        if (ncol <= 0) cycle
        do ifrag_row = 1, dg_frag%n_frag
          nrow = dg_frag%n_basis(ifrag_row, ispin)
          if (nrow <= 0) cycle
          valid_row_count = 0
          do ii = 1, nrow
            row_gid(ii) = dg_frag%index_basis(ii, ifrag_row, ispin)
            if (row_gid(ii) < 1 .or. row_gid(ii) > size(mat, 1)) cycle
            valid_row_count = valid_row_count + 1
            valid_row_ids(valid_row_count) = ii
          end do
          valid_col_count = 0
          do jj = 1, ncol
            col_gid(jj) = dg_frag%index_basis(jj, ifrag_col, ispin)
            if (col_gid(jj) < 1 .or. col_gid(jj) > size(mat, 2)) cycle
            valid_col_count = valid_col_count + 1
            valid_col_ids(valid_col_count) = jj
          end do
          block_size = valid_row_count * valid_col_count
          do idx_jj = 1, valid_col_count
            jj = valid_col_ids(idx_jj)
            ig_j = col_gid(jj)
            do idx_ii = 1, valid_row_count
              ii = valid_row_ids(idx_ii)
              ig_i = row_gid(ii)
              offset_flat = (idx_jj - 1) * valid_row_count + idx_ii
              send_block(offset_flat) = mat(ig_i, ig_j, ispin)
            end do
          end do

          chunk_begin = 1
          do while (chunk_begin <= block_size)
            chunk_count = min(reduce_chunk_size, block_size - chunk_begin + 1)
            call comm_summation(send_block(chunk_begin:chunk_begin + chunk_count - 1), &
                                recv_block(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, icomm_reduce)
            chunk_begin = chunk_begin + chunk_count
          end do

          do idx_jj = 1, valid_col_count
            jj = valid_col_ids(idx_jj)
            ig_j = col_gid(jj)
            do idx_ii = 1, valid_row_count
              ii = valid_row_ids(idx_ii)
              ig_i = row_gid(ii)
              offset_flat = (idx_jj - 1) * valid_row_count + idx_ii
              mat(ig_i, ig_j, ispin) = recv_block(offset_flat)
            end do
          end do
        end do
      end do
    end do

    deallocate(row_gid, col_gid, valid_row_ids, valid_col_ids)
    deallocate(send_block, recv_block)
    if (enable_reduce_trace .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0)') "        hamiltonian block reduce done: label=", trim(label), &
        " total_active=", total_active_size
      flush(6)
    end if
  end subroutine reduce_matrix_fragment_blocks

  !=======================================================================
  ! Build total local potential on the given grid:
  !   V_total = Vpsl + Vh + Vxc(ispin)
  !=======================================================================
  subroutine build_total_potential_grid(grid, Vh, Vxc_spin, Vpsl, V_total)
    use structures
    implicit none
    type(s_rgrid), intent(in) :: grid
    type(s_scalar), intent(in) :: Vh, Vxc_spin, Vpsl
    real(8), intent(out) :: V_total(grid%is(1):grid%ie(1), grid%is(2):grid%ie(2), grid%is(3):grid%ie(3))
    integer :: ix, iy, iz
    integer :: grid_x_lo, grid_x_hi, grid_y_lo, grid_y_hi, grid_z_lo, grid_z_hi

    grid_x_lo = grid%is(1)
    grid_x_hi = grid%ie(1)
    grid_y_lo = grid%is(2)
    grid_y_hi = grid%ie(2)
    grid_z_lo = grid%is(3)
    grid_z_hi = grid%ie(3)
    do iz = grid_z_lo, grid_z_hi
      do iy = grid_y_lo, grid_y_hi
        do ix = grid_x_lo, grid_x_hi
          V_total(ix, iy, iz) = Vpsl%f(ix, iy, iz) + Vh%f(ix, iy, iz) + Vxc_spin%f(ix, iy, iz)
        end do
      end do
    end do
  end subroutine build_total_potential_grid

  !=======================================================================
  ! Build T|phi_j> and H|phi_j>=T|phi_j>+V|phi_j> for one fragment/basis state
  !=======================================================================
  subroutine build_hpsi_for_basis(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ifrag, i_local, jo
    type(s_rgrid), intent(in) :: mg
    type(s_stencil), intent(in) :: stencil
    real(8), intent(in) :: V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8), intent(out) :: T_phi(:,:,:)
    real(8), intent(out) :: H_phi(:,:,:)
    integer :: gx, gy, gz, lx, ly, lz
    integer :: bx, by, bz
    integer :: gx0, gx1, gy0, gy1, gz0, gz1
    integer :: iorg(3), nsup(3)
    integer :: loc_s(3), loc_e(3)
    integer :: phi_lb1, phi_ub1, phi_lb2, phi_ub2, phi_lb3, phi_ub3
    integer :: v_lb1, v_ub1, v_lb2, v_ub2, v_lb3, v_ub3
    logical :: has_overlap

    if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) then
      write(*,*) "[FATAL] build_hpsi invalid i_local: rank=", dg_frag%id, &
        " ifrag=", ifrag, " i_local=", i_local, " phi_frag_dim5=", size(dg_frag%phi_frag, 5)
      stop 1
    end if
    if (jo < 1 .or. jo > size(dg_frag%phi_frag, 4)) then
      write(*,*) "[FATAL] build_hpsi invalid jo: rank=", dg_frag%id, &
        " ifrag=", ifrag, " jo=", jo, " phi_frag_dim4=", size(dg_frag%phi_frag, 4)
      stop 1
    end if
    call apply_kinetic_to_basis(dg_frag, i_local, jo, mg, stencil, T_phi)
    H_phi(:, :, :) = T_phi(:, :, :)

    call get_fragment_basis_owned_range(dg_frag, ifrag, mg, iorg, nsup, loc_s, loc_e, has_overlap)
    if (.not. has_overlap) return
    phi_lb1 = lbound(dg_frag%phi_frag, 1)
    phi_ub1 = ubound(dg_frag%phi_frag, 1)
    phi_lb2 = lbound(dg_frag%phi_frag, 2)
    phi_ub2 = ubound(dg_frag%phi_frag, 2)
    phi_lb3 = lbound(dg_frag%phi_frag, 3)
    phi_ub3 = ubound(dg_frag%phi_frag, 3)
    v_lb1 = lbound(V_total, 1)
    v_ub1 = ubound(V_total, 1)
    v_lb2 = lbound(V_total, 2)
    v_ub2 = ubound(V_total, 2)
    v_lb3 = lbound(V_total, 3)
    v_ub3 = ubound(V_total, 3)
    gx0 = iorg(1) + loc_s(1) - 1
    gx1 = iorg(1) + loc_e(1) - 1
    gy0 = iorg(2) + loc_s(2) - 1
    gy1 = iorg(2) + loc_e(2) - 1
    gz0 = iorg(3) + loc_s(3) - 1
    gz1 = iorg(3) + loc_e(3) - 1
    if (gx0 < v_lb1 .or. gx1 > v_ub1 .or. gy0 < v_lb2 .or. gy1 > v_ub2 .or. gz0 < v_lb3 .or. gz1 > v_ub3) then
      write(*,*) "[FATAL] build_hpsi V_total index out of bounds: rank=", dg_frag%id, &
        " ifrag=", ifrag, " g0=", gx0, gy0, gz0, " g1=", gx1, gy1, gz1, &
        " V_lb=", v_lb1, v_lb2, v_lb3, " V_ub=", v_ub1, v_ub2, v_ub3
      stop 1
    end if
!$omp parallel do private(lz, ly, lx, gz, gy, gx, bx, by, bz) schedule(static)
    do lz = loc_s(3), loc_e(3)
      gz = iorg(3) + lz - 1
      do ly = loc_s(2), loc_e(2)
        gy = iorg(2) + ly - 1
!$omp simd private(bx, by, bz)
        do lx = loc_s(1), loc_e(1)
          gx = iorg(1) + lx - 1
          bx = map_global_to_phi_box_coord_ham_fragment(dg_frag, ifrag, 1, gx, phi_lb1, phi_ub1)
          by = map_global_to_phi_box_coord_ham_fragment(dg_frag, ifrag, 2, gy, phi_lb2, phi_ub2)
          bz = map_global_to_phi_box_coord_ham_fragment(dg_frag, ifrag, 3, gz, phi_lb3, phi_ub3)
          if (bx == 0 .or. by == 0 .or. bz == 0) cycle
          H_phi(gx, gy, gz) = H_phi(gx, gy, gz) + V_total(gx, gy, gz) * dg_frag%phi_frag(bx, by, bz, jo, i_local)
        end do
      end do
    end do
!$omp end parallel do
  end subroutine build_hpsi_for_basis

  subroutine build_hpsi_for_basis_probe(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ifrag, i_local, jo
    type(s_rgrid), intent(in) :: mg
    type(s_stencil), intent(in) :: stencil
    real(8), intent(in) :: V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    complex(8), intent(out) :: T_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    complex(8), intent(out) :: H_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8), allocatable :: T_phi_re(:,:,:), H_phi_re(:,:,:)

    allocate(T_phi_re(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(H_phi_re(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    call build_hpsi_for_basis(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi_re, H_phi_re)
    T_phi(:, :, :) = cmplx(T_phi_re(:, :, :), 0.0d0, kind=8)
    H_phi(:, :, :) = cmplx(H_phi_re(:, :, :), 0.0d0, kind=8)
    deallocate(T_phi_re, H_phi_re)
  end subroutine build_hpsi_for_basis_probe

  subroutine get_phi_value_at_global_probe(dg_frag, ifrag, i_local, jo, gx, gy, gz, phi_val)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, i_local, jo, gx, gy, gz
    complex(8), intent(out) :: phi_val
    integer :: bx, by, bz

    bx = map_global_to_phi_box_coord_ham_fragment(dg_frag, ifrag, 1, gx, lbound(dg_frag%phi_frag, 1), ubound(dg_frag%phi_frag, 1))
    by = map_global_to_phi_box_coord_ham_fragment(dg_frag, ifrag, 2, gy, lbound(dg_frag%phi_frag, 2), ubound(dg_frag%phi_frag, 2))
    bz = map_global_to_phi_box_coord_ham_fragment(dg_frag, ifrag, 3, gz, lbound(dg_frag%phi_frag, 3), ubound(dg_frag%phi_frag, 3))
    if (bx == 0 .or. by == 0 .or. bz == 0) then
      phi_val = (0.0d0, 0.0d0)
      return
    end if
    phi_val = cmplx(dg_frag%phi_frag(bx, by, bz, jo, i_local), 0.0d0, kind=8)
  end subroutine get_phi_value_at_global_probe

  !=======================================================================
  ! Integrate one bra basis function against a real-space field
  !   integral = <phi_io | field>
  !=======================================================================
  subroutine integrate_basis_with_field(dg_frag, ifrag, i_local, io, mg, field, hvol, integral)
    use structures
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, i_local, io
    type(s_rgrid), intent(in) :: mg
    real(8), intent(in) :: field(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8), intent(in) :: hvol
    real(8), intent(out) :: integral
    real(8) :: partial
    integer :: gx, gy, gz, lx, ly, lz
    integer :: bx, by, bz
    integer :: iorg(3), nsup(3), loc_s(3), loc_e(3)
    integer :: gx0, gx1, gy0, gy1, gz0, gz1
    integer :: f_lb1, f_ub1, f_lb2, f_ub2, f_lb3, f_ub3
    logical :: has_overlap

    if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) then
      write(*,*) "[FATAL] integrate invalid i_local: rank=", dg_frag%id, &
        " ifrag=", ifrag, " i_local=", i_local, " phi_frag_dim5=", size(dg_frag%phi_frag, 5)
      stop 1
    end if
    if (io < 1 .or. io > size(dg_frag%phi_frag, 4)) then
      write(*,*) "[FATAL] integrate invalid io: rank=", dg_frag%id, &
        " ifrag=", ifrag, " io=", io, " phi_frag_dim4=", size(dg_frag%phi_frag, 4)
      stop 1
    end if
    if (dg_frag%use_buffered_basis) then
      iorg(:) = dg_frag%basis_support_lo(:, ifrag)
      nsup(:) = dg_frag%basis_support_hi(:, ifrag) - dg_frag%basis_support_lo(:, ifrag) + 1
    else
      iorg(:) = dg_frag%ixyz_frag(:, ifrag)
      nsup(:) = dg_frag%nxyz_domain(:, ifrag)
    end if
    call get_fragment_basis_owned_range(dg_frag, ifrag, mg, iorg, nsup, loc_s, loc_e, has_overlap)
    if (.not. has_overlap) then
      integral = 0.0d0
      return
    end if
    f_lb1 = lbound(dg_frag%phi_frag, 1)
    f_ub1 = ubound(dg_frag%phi_frag, 1)
    f_lb2 = lbound(dg_frag%phi_frag, 2)
    f_ub2 = ubound(dg_frag%phi_frag, 2)
    f_lb3 = lbound(dg_frag%phi_frag, 3)
    f_ub3 = ubound(dg_frag%phi_frag, 3)
    gx0 = iorg(1) + loc_s(1) - 1
    gx1 = iorg(1) + loc_e(1) - 1
    gy0 = iorg(2) + loc_s(2) - 1
    gy1 = iorg(2) + loc_e(2) - 1
    gz0 = iorg(3) + loc_s(3) - 1
    gz1 = iorg(3) + loc_e(3) - 1
    if (gx0 < lbound(field, 1) .or. gx1 > ubound(field, 1) .or. &
        gy0 < lbound(field, 2) .or. gy1 > ubound(field, 2) .or. &
        gz0 < lbound(field, 3) .or. gz1 > ubound(field, 3)) then
      write(*,*) "[FATAL] integrate global range exceeds field bounds: rank=", dg_frag%id, &
        " ifrag=", ifrag, " g0=", gx0, gy0, gz0, " g1=", gx1, gy1, gz1, " field_shape=", shape(field)
      stop 1
    end if
    partial = 0.0d0
    do lz = loc_s(3), loc_e(3)
      gz = iorg(3) + lz - 1
      do ly = loc_s(2), loc_e(2)
        gy = iorg(2) + ly - 1
        !$omp simd reduction(+:partial) private(bx, by, bz)
        do lx = loc_s(1), loc_e(1)
          gx = iorg(1) + lx - 1
          bx = map_global_to_phi_box_coord_ham_fragment(dg_frag, ifrag, 1, gx, f_lb1, f_ub1)
          by = map_global_to_phi_box_coord_ham_fragment(dg_frag, ifrag, 2, gy, f_lb2, f_ub2)
          bz = map_global_to_phi_box_coord_ham_fragment(dg_frag, ifrag, 3, gz, f_lb3, f_ub3)
          if (bx == 0 .or. by == 0 .or. bz == 0) cycle
          partial = partial + dg_frag%phi_frag(bx, by, bz, io, i_local) * field(gx, gy, gz) * hvol
        end do
      end do
    end do
    integral = partial
  end subroutine integrate_basis_with_field
  
  !=======================================================================
  ! Apply kinetic energy operator to a single basis function
  !
  ! T|φ> = -∇²/2 |φ> = -0.5 * Laplacian(φ)
  !
  ! Uses 4th-order finite difference stencil (requires ±4 grid points).
  ! With halo exchange, computation is valid over entire domain (1:nx, 1:ny, 1:nz).
  !
  ! System boundaries: PERIODIC (full system is periodic)
  ! Fragment boundaries: Halo exchange provides neighbor data via MPI
  !=======================================================================
  subroutine apply_kinetic_to_basis(dg_frag, i_local, jo, mg, stencil, T_phi)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer,                intent(in) :: i_local, jo
    type(s_rgrid),          intent(in) :: mg
    type(s_stencil),        intent(in) :: stencil
    real(8),                intent(out) :: T_phi(:,:,:)
    
    integer :: gx, gy, gz, ifrag, lx, ly, lz
    integer :: ix0, iy0, iz0
    integer :: lgx, lgy, lgz
    real(8) :: v, lap0
    real(8) :: lapt(4,3)
    integer :: iorg(3), nsup(3), loc_s(3), loc_e(3)
    integer :: phi_lb1, phi_lb2, phi_lb3, phi_ub1, phi_ub2, phi_ub3
    logical :: has_overlap
    real(8), allocatable :: phi_local(:,:,:)
    
    ! Extract stencil coefficients
    lap0 = stencil%coef_lap0
    lapt = stencil%coef_lap
    if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) then
      write(*,*) "[FATAL] kinetic invalid i_local: rank=", dg_frag%id, &
        " i_local=", i_local, " phi_frag_dim5=", size(dg_frag%phi_frag, 5)
      stop 1
    end if
    if (jo < 1 .or. jo > size(dg_frag%phi_frag, 4)) then
      write(*,*) "[FATAL] kinetic invalid jo: rank=", dg_frag%id, &
        " jo=", jo, " phi_frag_dim4=", size(dg_frag%phi_frag, 4)
      stop 1
    end if
    ifrag = dg_frag%ifrag_start + i_local - 1
    if (ifrag < dg_frag%ifrag_start .or. ifrag > dg_frag%ifrag_end) then
      write(*,*) "[FATAL] kinetic invalid ifrag from i_local: rank=", dg_frag%id, &
        " ifrag=", ifrag, " i_local=", i_local, " ifrag_start/end=", dg_frag%ifrag_start, dg_frag%ifrag_end
      stop 1
    end if
    call get_fragment_basis_owned_range(dg_frag, ifrag, mg, iorg, nsup, loc_s, loc_e, has_overlap)
    if (any(nsup <= 0)) then
      write(*,*) "[FATAL] kinetic non-positive basis support size: rank=", dg_frag%id, &
        " ifrag=", ifrag, " nsup=", nsup
      stop 1
    end if
    lgx = dg_frag%lgnum_total(1)
    lgy = dg_frag%lgnum_total(2)
    lgz = dg_frag%lgnum_total(3)
    phi_lb1 = lbound(dg_frag%phi_frag, 1)
    phi_lb2 = lbound(dg_frag%phi_frag, 2)
    phi_lb3 = lbound(dg_frag%phi_frag, 3)
    phi_ub1 = ubound(dg_frag%phi_frag, 1)
    phi_ub2 = ubound(dg_frag%phi_frag, 2)
    phi_ub3 = ubound(dg_frag%phi_frag, 3)
    T_phi = 0.0d0
    if (.not. has_overlap) return
    if (dg_frag%use_buffered_basis) then
      allocate(phi_local(-4:nsup(1)+4, -4:nsup(2)+4, -4:nsup(3)+4))
      call fill_buffered_support_local_box(dg_frag, ifrag, i_local, jo, iorg, nsup, phi_local)
    end if
    ! Note: phi_frag is allocated as (1-nb:nx+nb, 1-nb:ny+nb, 1-nb:nz+nb, ...)
    ! where nb = nxyz_buffer = 4 for 4th-order stencil
    ! The interior domain is (1:nx, 1:ny, 1:nz), and halo provides data for stencil
    ! operations near boundaries.
    !
    ! With halo exchange, stencil operations can access phi_frag(ix±4, iy±4, iz±4)
    ! for all interior points without boundary restrictions.
    
    ! Apply kinetic operator using finite difference stencil.
    ! With periodic halos already populated in phi_frag, convert each global
    ! coordinate to its periodic interior index once, then use ix0+/-n offsets
    ! directly for neighbor access.
    !
    ! Note: exchange_phi_frag_halo() must be called before this routine

!$omp parallel do private(lz, ly, lx, gz, gy, gx, ix0, iy0, iz0, v) schedule(static)
    do lz = loc_s(3), loc_e(3)
      gz = iorg(3) + lz - 1
      iz0 = modulo(gz - 1, lgz) + 1
      do ly = loc_s(2), loc_e(2)
        gy = iorg(2) + ly - 1
        iy0 = modulo(gy - 1, lgy) + 1
!$omp simd private(v, ix0)
        do lx = loc_s(1), loc_e(1)
          gx = iorg(1) + lx - 1
          ix0 = modulo(gx - 1, lgx) + 1
          if (dg_frag%use_buffered_basis) then
            v = lapt(1,1) * (phi_local(lx + 1, ly, lz) + phi_local(lx - 1, ly, lz)) + &
                lapt(2,1) * (phi_local(lx + 2, ly, lz) + phi_local(lx - 2, ly, lz)) + &
                lapt(3,1) * (phi_local(lx + 3, ly, lz) + phi_local(lx - 3, ly, lz)) + &
                lapt(4,1) * (phi_local(lx + 4, ly, lz) + phi_local(lx - 4, ly, lz))
            v = v + &
                lapt(1,2) * (phi_local(lx, ly + 1, lz) + phi_local(lx, ly - 1, lz)) + &
                lapt(2,2) * (phi_local(lx, ly + 2, lz) + phi_local(lx, ly - 2, lz)) + &
                lapt(3,2) * (phi_local(lx, ly + 3, lz) + phi_local(lx, ly - 3, lz)) + &
                lapt(4,2) * (phi_local(lx, ly + 4, lz) + phi_local(lx, ly - 4, lz))
            v = v + &
                lapt(1,3) * (phi_local(lx, ly, lz + 1) + phi_local(lx, ly, lz - 1)) + &
                lapt(2,3) * (phi_local(lx, ly, lz + 2) + phi_local(lx, ly, lz - 2)) + &
                lapt(3,3) * (phi_local(lx, ly, lz + 3) + phi_local(lx, ly, lz - 3)) + &
                lapt(4,3) * (phi_local(lx, ly, lz + 4) + phi_local(lx, ly, lz - 4))
          else
            ! Compute Laplacian using 4th-order finite difference
            ! Stencil accesses phi_frag(ix±4, iy±4, iz±4) which now includes halo
            v = lapt(1,1) * (dg_frag%phi_frag(ix0 + 1, iy0, iz0, jo, i_local) + &
                             dg_frag%phi_frag(ix0 - 1, iy0, iz0, jo, i_local)) + &
                lapt(2,1) * (dg_frag%phi_frag(ix0 + 2, iy0, iz0, jo, i_local) + &
                             dg_frag%phi_frag(ix0 - 2, iy0, iz0, jo, i_local)) + &
                lapt(3,1) * (dg_frag%phi_frag(ix0 + 3, iy0, iz0, jo, i_local) + &
                             dg_frag%phi_frag(ix0 - 3, iy0, iz0, jo, i_local)) + &
                lapt(4,1) * (dg_frag%phi_frag(ix0 + 4, iy0, iz0, jo, i_local) + &
                             dg_frag%phi_frag(ix0 - 4, iy0, iz0, jo, i_local))
            v = v + &
                lapt(1,2) * (dg_frag%phi_frag(ix0, iy0 + 1, iz0, jo, i_local) + &
                             dg_frag%phi_frag(ix0, iy0 - 1, iz0, jo, i_local)) + &
                lapt(2,2) * (dg_frag%phi_frag(ix0, iy0 + 2, iz0, jo, i_local) + &
                             dg_frag%phi_frag(ix0, iy0 - 2, iz0, jo, i_local)) + &
                lapt(3,2) * (dg_frag%phi_frag(ix0, iy0 + 3, iz0, jo, i_local) + &
                             dg_frag%phi_frag(ix0, iy0 - 3, iz0, jo, i_local)) + &
                lapt(4,2) * (dg_frag%phi_frag(ix0, iy0 + 4, iz0, jo, i_local) + &
                             dg_frag%phi_frag(ix0, iy0 - 4, iz0, jo, i_local))
            v = v + &
                lapt(1,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 1, jo, i_local) + &
                             dg_frag%phi_frag(ix0, iy0, iz0 - 1, jo, i_local)) + &
                lapt(2,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 2, jo, i_local) + &
                             dg_frag%phi_frag(ix0, iy0, iz0 - 2, jo, i_local)) + &
                lapt(3,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 3, jo, i_local) + &
                             dg_frag%phi_frag(ix0, iy0, iz0 - 3, jo, i_local)) + &
                lapt(4,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 4, jo, i_local) + &
                             dg_frag%phi_frag(ix0, iy0, iz0 - 4, jo, i_local))
          end if
          
          ! T|φ> = (-∇²/2)|φ> = lap0*φ - 0.5 * (sum of neighbor terms)
          if (dg_frag%use_buffered_basis) then
            T_phi(gx, gy, gz) = lap0 * phi_local(lx, ly, lz) - 0.5d0 * v
          else
            T_phi(gx, gy, gz) = lap0 * dg_frag%phi_frag(ix0, iy0, iz0, jo, i_local) - 0.5d0 * v
          end if
          
        end do
      end do
    end do
!$omp end parallel do
    if (allocated(phi_local)) deallocate(phi_local)
    
  end subroutine apply_kinetic_to_basis

  subroutine get_fragment_owned_range(dg_frag, ifrag, mg, loc_s, loc_e, has_overlap)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag
    type(s_rgrid), intent(in) :: mg
    integer, intent(out) :: loc_s(3), loc_e(3)
    logical, intent(out) :: has_overlap

    integer :: iorg(3), ndom(3), g_s(3), g_e(3), ov_s(3), ov_e(3)

    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    g_s(:) = iorg(:)
    g_e(:) = iorg(:) + ndom(:) - 1
    ov_s(:) = max(g_s(:), mg%is(:))
    ov_e(:) = min(g_e(:), mg%ie(:))

    has_overlap = all(ov_s(:) <= ov_e(:))
    if (.not. has_overlap) then
      loc_s(:) = 1
      loc_e(:) = 0
      return
    end if

    loc_s(:) = ov_s(:) - iorg(:) + 1
    loc_e(:) = ov_e(:) - iorg(:) + 1
  end subroutine get_fragment_owned_range

  subroutine get_fragment_basis_owned_range(dg_frag, ifrag, mg, iorg, nsup, loc_s, loc_e, has_overlap)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag
    type(s_rgrid), intent(in) :: mg
    integer, intent(out) :: iorg(3), nsup(3), loc_s(3), loc_e(3)
    logical, intent(out) :: has_overlap

    integer :: g_s(3), g_e(3), ov_s(3), ov_e(3)

    if (dg_frag%use_buffered_basis) then
      iorg(:) = dg_frag%basis_support_lo(:, ifrag)
      nsup(:) = dg_frag%basis_support_hi(:, ifrag) - dg_frag%basis_support_lo(:, ifrag) + 1
    else
      iorg(:) = dg_frag%ixyz_frag(:, ifrag)
      nsup(:) = dg_frag%nxyz_domain(:, ifrag)
    end if

    g_s(:) = iorg(:)
    g_e(:) = iorg(:) + nsup(:) - 1
    ov_s(:) = max(g_s(:), mg%is(:))
    ov_e(:) = min(g_e(:), mg%ie(:))

    has_overlap = all(ov_s(:) <= ov_e(:))
    if (.not. has_overlap) then
      loc_s(:) = 1
      loc_e(:) = 0
      return
    end if

    loc_s(:) = ov_s(:) - iorg(:) + 1
    loc_e(:) = ov_e(:) - iorg(:) + 1
  end subroutine get_fragment_basis_owned_range

  subroutine get_fragment_local_range(dg_frag, ndom, loc_s, loc_e)
    use salmon_global, only: nproc_rgrid
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ndom(3)
    integer, intent(out) :: loc_s(3), loc_e(3)

    integer :: ipx, ipy, ipz, coords(3), nsize

    ipx = max(1, nproc_rgrid(1))
    ipy = max(1, nproc_rgrid(2))
    ipz = max(1, nproc_rgrid(3))

    if (dg_frag%id_frag < 0 .or. dg_frag%id_frag >= ipx * ipy * ipz) then
      stop "DG-Fragment RT: invalid fragment-local MPI rank in get_fragment_local_range"
    end if

    coords(1) = mod(dg_frag%id_frag, ipx)
    coords(2) = mod(dg_frag%id_frag / ipx, ipy)
    coords(3) = dg_frag%id_frag / max(1, ipx * ipy)

    nsize = (ndom(1) + ipx - 1) / ipx
    loc_s(1) = 1 + nsize * coords(1)
    loc_e(1) = min(ndom(1), loc_s(1) + nsize - 1)

    nsize = (ndom(2) + ipy - 1) / ipy
    loc_s(2) = 1 + nsize * coords(2)
    loc_e(2) = min(ndom(2), loc_s(2) + nsize - 1)

    nsize = (ndom(3) + ipz - 1) / ipz
    loc_s(3) = 1 + nsize * coords(3)
    loc_e(3) = min(ndom(3), loc_s(3) + nsize - 1)
  end subroutine get_fragment_local_range

  subroutine get_halo_face_axis_dir(halo, axis, dir)
    implicit none
    type(halo_info), intent(in) :: halo
    integer, intent(out) :: axis, dir
    integer :: idir

    axis = 0
    dir = 0
    do idir = 1, 3
      if (halo%dvec(idir) == 0) cycle
      if (axis /= 0) stop "DG-Fragment RT: halo spans multiple face normals"
      axis = idir
      dir = halo%dvec(idir)
    end do
    if (axis == 0 .or. abs(dir) /= 1) stop "DG-Fragment RT: invalid halo face direction"
  end subroutine get_halo_face_axis_dir

  logical function is_face_halo_neighbor(halo) result(is_face)
    implicit none
    type(halo_info), intent(in) :: halo

    is_face = (count(halo%dvec /= 0) == 1)
  end function is_face_halo_neighbor

  real(8) function get_dg_face_area_element(hgs, axis) result(face_area)
    implicit none
    real(8), intent(in) :: hgs(3)
    integer, intent(in) :: axis

    select case (axis)
    case (1)
      face_area = hgs(2) * hgs(3)
    case (2)
      face_area = hgs(1) * hgs(3)
    case (3)
      face_area = hgs(1) * hgs(2)
    case default
      stop "DG-Fragment RT: invalid face axis"
    end select
  end function get_dg_face_area_element

  real(8) function get_dg_face_penalty_coefficient(hgs, axis) result(alpha_penalty)
    implicit none
    real(8), intent(in) :: hgs(3)
    integer, intent(in) :: axis

    if (axis < 1 .or. axis > 3) stop "DG-Fragment RT: invalid penalty axis"
    if (hgs(axis) <= 0.0d0) stop "DG-Fragment RT: non-positive penalty spacing"
    alpha_penalty = 1.0d0 / hgs(axis)
  end function get_dg_face_penalty_coefficient

  integer function get_local_face_boundary_index(len_axis, face_dir) result(idx_face)
    implicit none
    integer, intent(in) :: len_axis, face_dir

    if (face_dir > 0) then
      idx_face = len_axis
    else
      idx_face = 1
    end if
  end function get_local_face_boundary_index

  integer function get_remote_face_boundary_index(len_axis, face_dir) result(idx_face)
    implicit none
    integer, intent(in) :: len_axis, face_dir

    if (face_dir > 0) then
      idx_face = 1
    else
      idx_face = len_axis
    end if
  end function get_remote_face_boundary_index

  real(8) function face_trace_value_local(dg_frag, i_local, jo, axis, face_dir, idx1, idx2, idx3) result(phi_face)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: i_local, jo, axis, face_dir
    integer, intent(in) :: idx1, idx2, idx3
    if (axis < 1 .or. axis > 3) stop "DG-Fragment RT: invalid local face axis"
    if (abs(face_dir) /= 1) stop "DG-Fragment RT: invalid local face direction"
    phi_face = dg_frag%phi_frag(idx1, idx2, idx3, jo, i_local)
  end function face_trace_value_local

  real(8) function face_trace_value_halo(buf_recv, length, io, axis, face_dir, idx1, idx2, idx3) result(phi_face)
    implicit none
    real(8), intent(in) :: buf_recv(:,:,:,:,:)
    integer, intent(in) :: length(3), io, axis, face_dir
    integer, intent(in) :: idx1, idx2, idx3
    if (axis < 1 .or. axis > 3) stop "DG-Fragment RT: invalid halo face axis"
    if (abs(face_dir) /= 1) stop "DG-Fragment RT: invalid halo face direction"
    if (idx1 < 1 .or. idx1 > length(1)) stop "DG-Fragment RT: invalid halo face idx1"
    if (idx2 < 1 .or. idx2 > length(2)) stop "DG-Fragment RT: invalid halo face idx2"
    if (idx3 < 1 .or. idx3 > length(3)) stop "DG-Fragment RT: invalid halo face idx3"
    phi_face = buf_recv(idx1, idx2, idx3, io, 1)
  end function face_trace_value_halo

  real(8) function face_trace_value_flow_halo(buf_recv, length, io, axis, idx1, idx2, idx3) result(phi_face)
    implicit none
    real(8), intent(in) :: buf_recv(:,:,:,:,:)
    integer, intent(in) :: length(3), io, axis
    integer, intent(in) :: idx1, idx2, idx3

    select case (axis)
    case (1)
      phi_face = buf_recv(1, idx2, idx3, io, 1)
    case (2)
      phi_face = buf_recv(idx1, 1, idx3, io, 1)
    case (3)
      phi_face = buf_recv(idx1, idx2, 1, io, 1)
    case default
      stop "DG-Fragment RT: invalid flow halo face axis"
    end select
  end function face_trace_value_flow_halo

  real(8) function one_sided_face_derivative_local(dg_frag, i_local, jo, axis, face_dir, idx1, idx2, idx3) result(dphi_dn)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: i_local, jo, axis, face_dir
    integer, intent(in) :: idx1, idx2, idx3
    real(8) :: f0, f1, f2, h

    select case (axis)
    case (1)
      f0 = dg_frag%phi_frag(idx1, idx2, idx3, jo, i_local)
      f1 = dg_frag%phi_frag(idx1-face_dir, idx2, idx3, jo, i_local)
      f2 = dg_frag%phi_frag(idx1-2*face_dir, idx2, idx3, jo, i_local)
      h = dg_frag%hgs(1)
    case (2)
      f0 = dg_frag%phi_frag(idx1, idx2, idx3, jo, i_local)
      f1 = dg_frag%phi_frag(idx1, idx2-face_dir, idx3, jo, i_local)
      f2 = dg_frag%phi_frag(idx1, idx2-2*face_dir, idx3, jo, i_local)
      h = dg_frag%hgs(2)
    case (3)
      f0 = dg_frag%phi_frag(idx1, idx2, idx3, jo, i_local)
      f1 = dg_frag%phi_frag(idx1, idx2, idx3-face_dir, jo, i_local)
      f2 = dg_frag%phi_frag(idx1, idx2, idx3-2*face_dir, jo, i_local)
      h = dg_frag%hgs(3)
    case default
      stop "DG-Fragment RT: invalid local face axis"
    end select

    ! Local samples are taken into element interior (opposite to outward normal),
    ! so the outward-normal derivative is the sign-flipped inward one-sided value.
    dphi_dn = (3.0d0 * f0 - 4.0d0 * f1 + f2) / (2.0d0 * h)
  end function one_sided_face_derivative_local

  real(8) function one_sided_face_derivative_halo(buf_recv, length, io, axis, face_dir, idx1, idx2, idx3, hgs) result(dphi_dn)
    implicit none
    real(8), intent(in) :: buf_recv(:,:,:,:,:)
    integer, intent(in) :: length(3), io, axis, face_dir
    integer, intent(in) :: idx1, idx2, idx3
    real(8), intent(in) :: hgs(3)
    real(8) :: f0, f1, f2, h
    integer :: i0, i1, i2

    select case (axis)
    case (1)
      if (face_dir > 0) then
        i0 = 1; i1 = 2; i2 = 3
      else
        i0 = length(1); i1 = length(1)-1; i2 = length(1)-2
      end if
      f0 = buf_recv(i0, idx2, idx3, io, 1)
      f1 = buf_recv(i1, idx2, idx3, io, 1)
      f2 = buf_recv(i2, idx2, idx3, io, 1)
      h = hgs(1)
    case (2)
      if (face_dir > 0) then
        i0 = 1; i1 = 2; i2 = 3
      else
        i0 = length(2); i1 = length(2)-1; i2 = length(2)-2
      end if
      f0 = buf_recv(idx1, i0, idx3, io, 1)
      f1 = buf_recv(idx1, i1, idx3, io, 1)
      f2 = buf_recv(idx1, i2, idx3, io, 1)
      h = hgs(2)
    case (3)
      if (face_dir > 0) then
        i0 = 1; i1 = 2; i2 = 3
      else
        i0 = length(3); i1 = length(3)-1; i2 = length(3)-2
      end if
      f0 = buf_recv(idx1, idx2, i0, io, 1)
      f1 = buf_recv(idx1, idx2, i1, io, 1)
      f2 = buf_recv(idx1, idx2, i2, io, 1)
      h = hgs(3)
    case default
      stop "DG-Fragment RT: invalid halo face axis"
    end select

    ! Halo samples are ordered from interface into remote interior,
    ! aligned with local outward normal; use matching 2nd-order one-sided form.
    dphi_dn = (-3.0d0 * f0 + 4.0d0 * f1 - f2) / (2.0d0 * h)
  end function one_sided_face_derivative_halo

  real(8) function face_derivative_value_flow_halo(buf_recv, length, io, axis, idx1, idx2, idx3, hgs) result(dphi_dn)
    implicit none
    real(8), intent(in) :: buf_recv(:,:,:,:,:)
    integer, intent(in) :: length(3), io, axis
    integer, intent(in) :: idx1, idx2, idx3
    real(8), intent(in) :: hgs(3)

    select case (axis)
    case (1)
      dphi_dn = buf_recv(1, idx2, idx3, io, 2)
    case (2)
      dphi_dn = buf_recv(idx1, 1, idx3, io, 2)
    case (3)
      dphi_dn = buf_recv(idx1, idx2, 1, io, 2)
    case default
      stop "DG-Fragment RT: invalid flow halo derivative axis"
    end select
  end function face_derivative_value_flow_halo

  subroutine apply_gradient_to_basis_ops_local_2d(dg_frag, i_local, jo, mg, stencil, loc_s, loc_e, grad_phi, grad_local_2d)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer,                intent(in) :: i_local, jo
    type(s_rgrid),          intent(in) :: mg
    type(s_stencil),        intent(in) :: stencil
    integer,                intent(in) :: loc_s(3), loc_e(3)
    real(8),                intent(out) :: grad_phi(:,:,:,:)
    real(8),                intent(out) :: grad_local_2d(:,:)

    integer :: lx, ly, lz, ifrag, iorg(3), nsup(3), nloc1, nloc2, ipt
    integer :: gxg, gyg, gzg
    integer :: ix0, iy0, iz0
    integer :: p_lb1, p_ub1, p_lb2, p_ub2, p_lb3, p_ub3
    integer :: sxp1, sxm1, sxp2, sxm2, sxp3, sxm3, sxp4, sxm4
    integer :: syp1, sym1, syp2, sym2, syp3, sym3, syp4, sym4
    integer :: szp1, szm1, szp2, szm2, szp3, szm3, szp4, szm4
    real(8) :: nabt(4,3), gx, gy, gz
    logical :: has_overlap

    nabt = stencil%coef_nab
    ifrag = dg_frag%ifrag_start + i_local - 1
    if (dg_frag%use_buffered_basis) then
      iorg(:) = dg_frag%basis_support_lo(:, ifrag)
      nsup(:) = dg_frag%basis_support_hi(:, ifrag) - dg_frag%basis_support_lo(:, ifrag) + 1
    else
      iorg(:) = dg_frag%ixyz_frag(:, ifrag)
      nsup(:) = dg_frag%nxyz_domain(:, ifrag)
    end if
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub3 = ubound(dg_frag%phi_frag, 3)
    nloc1 = loc_e(1) - loc_s(1) + 1
    nloc2 = loc_e(2) - loc_s(2) + 1
    !$omp parallel do private(lx, ly, lz, ipt, gxg, gyg, gzg, gx, gy, gz, ix0, iy0, iz0) schedule(static)
    do lz = 1, nsup(3)
      do ly = 1, nsup(2)
        !$omp simd private(ipt, gxg, gyg, gzg, gx, gy, gz, ix0, iy0, iz0)
        do lx = 1, nsup(1)
          gxg = iorg(1) + lx - 1
          gyg = iorg(2) + ly - 1
          gzg = iorg(3) + lz - 1
          ix0 = map_global_to_phi_box_coord_ham(gxg, p_lb1, p_ub1, dg_frag%lgnum_total(1))
          iy0 = map_global_to_phi_box_coord_ham(gyg, p_lb2, p_ub2, dg_frag%lgnum_total(2))
          iz0 = map_global_to_phi_box_coord_ham(gzg, p_lb3, p_ub3, dg_frag%lgnum_total(3))
          if (ix0 == 0 .or. iy0 == 0 .or. iz0 == 0) then
            gx = 0.0d0
            gy = 0.0d0
            gz = 0.0d0
          else
            if (dg_frag%use_buffered_basis) then
              sxp1 = map_global_to_phi_box_coord_ham(iorg(1) + wrap_support_local_index(lx + 1, nsup(1)) - 1, p_lb1, p_ub1, dg_frag%lgnum_total(1))
              sxm1 = map_global_to_phi_box_coord_ham(iorg(1) + wrap_support_local_index(lx - 1, nsup(1)) - 1, p_lb1, p_ub1, dg_frag%lgnum_total(1))
              sxp2 = map_global_to_phi_box_coord_ham(iorg(1) + wrap_support_local_index(lx + 2, nsup(1)) - 1, p_lb1, p_ub1, dg_frag%lgnum_total(1))
              sxm2 = map_global_to_phi_box_coord_ham(iorg(1) + wrap_support_local_index(lx - 2, nsup(1)) - 1, p_lb1, p_ub1, dg_frag%lgnum_total(1))
              sxp3 = map_global_to_phi_box_coord_ham(iorg(1) + wrap_support_local_index(lx + 3, nsup(1)) - 1, p_lb1, p_ub1, dg_frag%lgnum_total(1))
              sxm3 = map_global_to_phi_box_coord_ham(iorg(1) + wrap_support_local_index(lx - 3, nsup(1)) - 1, p_lb1, p_ub1, dg_frag%lgnum_total(1))
              sxp4 = map_global_to_phi_box_coord_ham(iorg(1) + wrap_support_local_index(lx + 4, nsup(1)) - 1, p_lb1, p_ub1, dg_frag%lgnum_total(1))
              sxm4 = map_global_to_phi_box_coord_ham(iorg(1) + wrap_support_local_index(lx - 4, nsup(1)) - 1, p_lb1, p_ub1, dg_frag%lgnum_total(1))
              syp1 = map_global_to_phi_box_coord_ham(iorg(2) + wrap_support_local_index(ly + 1, nsup(2)) - 1, p_lb2, p_ub2, dg_frag%lgnum_total(2))
              sym1 = map_global_to_phi_box_coord_ham(iorg(2) + wrap_support_local_index(ly - 1, nsup(2)) - 1, p_lb2, p_ub2, dg_frag%lgnum_total(2))
              syp2 = map_global_to_phi_box_coord_ham(iorg(2) + wrap_support_local_index(ly + 2, nsup(2)) - 1, p_lb2, p_ub2, dg_frag%lgnum_total(2))
              sym2 = map_global_to_phi_box_coord_ham(iorg(2) + wrap_support_local_index(ly - 2, nsup(2)) - 1, p_lb2, p_ub2, dg_frag%lgnum_total(2))
              syp3 = map_global_to_phi_box_coord_ham(iorg(2) + wrap_support_local_index(ly + 3, nsup(2)) - 1, p_lb2, p_ub2, dg_frag%lgnum_total(2))
              sym3 = map_global_to_phi_box_coord_ham(iorg(2) + wrap_support_local_index(ly - 3, nsup(2)) - 1, p_lb2, p_ub2, dg_frag%lgnum_total(2))
              syp4 = map_global_to_phi_box_coord_ham(iorg(2) + wrap_support_local_index(ly + 4, nsup(2)) - 1, p_lb2, p_ub2, dg_frag%lgnum_total(2))
              sym4 = map_global_to_phi_box_coord_ham(iorg(2) + wrap_support_local_index(ly - 4, nsup(2)) - 1, p_lb2, p_ub2, dg_frag%lgnum_total(2))
              szp1 = map_global_to_phi_box_coord_ham(iorg(3) + wrap_support_local_index(lz + 1, nsup(3)) - 1, p_lb3, p_ub3, dg_frag%lgnum_total(3))
              szm1 = map_global_to_phi_box_coord_ham(iorg(3) + wrap_support_local_index(lz - 1, nsup(3)) - 1, p_lb3, p_ub3, dg_frag%lgnum_total(3))
              szp2 = map_global_to_phi_box_coord_ham(iorg(3) + wrap_support_local_index(lz + 2, nsup(3)) - 1, p_lb3, p_ub3, dg_frag%lgnum_total(3))
              szm2 = map_global_to_phi_box_coord_ham(iorg(3) + wrap_support_local_index(lz - 2, nsup(3)) - 1, p_lb3, p_ub3, dg_frag%lgnum_total(3))
              szp3 = map_global_to_phi_box_coord_ham(iorg(3) + wrap_support_local_index(lz + 3, nsup(3)) - 1, p_lb3, p_ub3, dg_frag%lgnum_total(3))
              szm3 = map_global_to_phi_box_coord_ham(iorg(3) + wrap_support_local_index(lz - 3, nsup(3)) - 1, p_lb3, p_ub3, dg_frag%lgnum_total(3))
              szp4 = map_global_to_phi_box_coord_ham(iorg(3) + wrap_support_local_index(lz + 4, nsup(3)) - 1, p_lb3, p_ub3, dg_frag%lgnum_total(3))
              szm4 = map_global_to_phi_box_coord_ham(iorg(3) + wrap_support_local_index(lz - 4, nsup(3)) - 1, p_lb3, p_ub3, dg_frag%lgnum_total(3))
              if (sxp1 < 1 .or. sxm1 < 1 .or. sxp2 < 1 .or. sxm2 < 1 .or. sxp3 < 1 .or. sxm3 < 1 .or. sxp4 < 1 .or. sxm4 < 1) then
                gx = 0.0d0
              else
                gx = nabt(1,1) * (dg_frag%phi_frag(sxp1, iy0, iz0, jo, i_local) - dg_frag%phi_frag(sxm1, iy0, iz0, jo, i_local)) + &
                     nabt(2,1) * (dg_frag%phi_frag(sxp2, iy0, iz0, jo, i_local) - dg_frag%phi_frag(sxm2, iy0, iz0, jo, i_local)) + &
                     nabt(3,1) * (dg_frag%phi_frag(sxp3, iy0, iz0, jo, i_local) - dg_frag%phi_frag(sxm3, iy0, iz0, jo, i_local)) + &
                     nabt(4,1) * (dg_frag%phi_frag(sxp4, iy0, iz0, jo, i_local) - dg_frag%phi_frag(sxm4, iy0, iz0, jo, i_local))
              end if
              if (syp1 < 1 .or. sym1 < 1 .or. syp2 < 1 .or. sym2 < 1 .or. syp3 < 1 .or. sym3 < 1 .or. syp4 < 1 .or. sym4 < 1) then
                gy = 0.0d0
              else
                gy = nabt(1,2) * (dg_frag%phi_frag(ix0, syp1, iz0, jo, i_local) - dg_frag%phi_frag(ix0, sym1, iz0, jo, i_local)) + &
                     nabt(2,2) * (dg_frag%phi_frag(ix0, syp2, iz0, jo, i_local) - dg_frag%phi_frag(ix0, sym2, iz0, jo, i_local)) + &
                     nabt(3,2) * (dg_frag%phi_frag(ix0, syp3, iz0, jo, i_local) - dg_frag%phi_frag(ix0, sym3, iz0, jo, i_local)) + &
                     nabt(4,2) * (dg_frag%phi_frag(ix0, syp4, iz0, jo, i_local) - dg_frag%phi_frag(ix0, sym4, iz0, jo, i_local))
              end if
              if (szp1 < 1 .or. szm1 < 1 .or. szp2 < 1 .or. szm2 < 1 .or. szp3 < 1 .or. szm3 < 1 .or. szp4 < 1 .or. szm4 < 1) then
                gz = 0.0d0
              else
                gz = nabt(1,3) * (dg_frag%phi_frag(ix0, iy0, szp1, jo, i_local) - dg_frag%phi_frag(ix0, iy0, szm1, jo, i_local)) + &
                     nabt(2,3) * (dg_frag%phi_frag(ix0, iy0, szp2, jo, i_local) - dg_frag%phi_frag(ix0, iy0, szm2, jo, i_local)) + &
                     nabt(3,3) * (dg_frag%phi_frag(ix0, iy0, szp3, jo, i_local) - dg_frag%phi_frag(ix0, iy0, szm3, jo, i_local)) + &
                     nabt(4,3) * (dg_frag%phi_frag(ix0, iy0, szp4, jo, i_local) - dg_frag%phi_frag(ix0, iy0, szm4, jo, i_local))
              end if
            else
              gx = nabt(1,1) * (dg_frag%phi_frag(ix0 + 1, iy0, iz0, jo, i_local) - &
                                dg_frag%phi_frag(ix0 - 1, iy0, iz0, jo, i_local)) + &
                   nabt(2,1) * (dg_frag%phi_frag(ix0 + 2, iy0, iz0, jo, i_local) - &
                                dg_frag%phi_frag(ix0 - 2, iy0, iz0, jo, i_local)) + &
                   nabt(3,1) * (dg_frag%phi_frag(ix0 + 3, iy0, iz0, jo, i_local) - &
                                dg_frag%phi_frag(ix0 - 3, iy0, iz0, jo, i_local)) + &
                   nabt(4,1) * (dg_frag%phi_frag(ix0 + 4, iy0, iz0, jo, i_local) - &
                                dg_frag%phi_frag(ix0 - 4, iy0, iz0, jo, i_local))
              gy = nabt(1,2) * (dg_frag%phi_frag(ix0, iy0 + 1, iz0, jo, i_local) - &
                                dg_frag%phi_frag(ix0, iy0 - 1, iz0, jo, i_local)) + &
                   nabt(2,2) * (dg_frag%phi_frag(ix0, iy0 + 2, iz0, jo, i_local) - &
                                dg_frag%phi_frag(ix0, iy0 - 2, iz0, jo, i_local)) + &
                   nabt(3,2) * (dg_frag%phi_frag(ix0, iy0 + 3, iz0, jo, i_local) - &
                                dg_frag%phi_frag(ix0, iy0 - 3, iz0, jo, i_local)) + &
                   nabt(4,2) * (dg_frag%phi_frag(ix0, iy0 + 4, iz0, jo, i_local) - &
                                dg_frag%phi_frag(ix0, iy0 - 4, iz0, jo, i_local))
              gz = nabt(1,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 1, jo, i_local) - &
                                dg_frag%phi_frag(ix0, iy0, iz0 - 1, jo, i_local)) + &
                   nabt(2,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 2, jo, i_local) - &
                                dg_frag%phi_frag(ix0, iy0, iz0 - 2, jo, i_local)) + &
                   nabt(3,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 3, jo, i_local) - &
                                dg_frag%phi_frag(ix0, iy0, iz0 - 3, jo, i_local)) + &
                   nabt(4,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 4, jo, i_local) - &
                                dg_frag%phi_frag(ix0, iy0, iz0 - 4, jo, i_local))
            end if
          end if

          grad_phi(lx, ly, lz, 1) = gx
          grad_phi(lx, ly, lz, 2) = gy
          grad_phi(lx, ly, lz, 3) = gz

          if (lx >= loc_s(1) .and. lx <= loc_e(1) .and. &
              ly >= loc_s(2) .and. ly <= loc_e(2) .and. &
              lz >= loc_s(3) .and. lz <= loc_e(3)) then
            ipt = ((lz - loc_s(3)) * nloc2 + (ly - loc_s(2))) * nloc1 + (lx - loc_s(1)) + 1
            grad_local_2d(ipt, 1) = gx
            grad_local_2d(ipt, 2) = gy
            grad_local_2d(ipt, 3) = gz
          end if
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine apply_gradient_to_basis_ops_local_2d

  subroutine apply_gradient_at_phi_box_point(dg_frag, i_local, jo, mg, stencil, phi_idx, grad_x, grad_y, grad_z)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer,                intent(in) :: i_local, jo
    type(s_rgrid),          intent(in) :: mg
    type(s_stencil),        intent(in) :: stencil
    integer,                intent(in) :: phi_idx(3)
    real(8),                intent(out) :: grad_x, grad_y, grad_z

    integer :: ix0, iy0, iz0
    integer :: p_lb1, p_lb2, p_lb3, p_ub1, p_ub2, p_ub3
    real(8) :: nabt(4,3)

    nabt = stencil%coef_nab
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub3 = ubound(dg_frag%phi_frag, 3)
    ix0 = map_global_to_phi_box_coord_ham(phi_idx(1), p_lb1, p_ub1, dg_frag%lgnum_total(1))
    iy0 = map_global_to_phi_box_coord_ham(phi_idx(2), p_lb2, p_ub2, dg_frag%lgnum_total(2))
    iz0 = map_global_to_phi_box_coord_ham(phi_idx(3), p_lb3, p_ub3, dg_frag%lgnum_total(3))
    if (ix0 == 0 .or. iy0 == 0 .or. iz0 == 0) then
      grad_x = 0.0d0
      grad_y = 0.0d0
      grad_z = 0.0d0
      return
    end if
    grad_x = nabt(1,1) * (dg_frag%phi_frag(ix0 + 1, iy0, iz0, jo, i_local) - &
                           dg_frag%phi_frag(ix0 - 1, iy0, iz0, jo, i_local)) + &
             nabt(2,1) * (dg_frag%phi_frag(ix0 + 2, iy0, iz0, jo, i_local) - &
                           dg_frag%phi_frag(ix0 - 2, iy0, iz0, jo, i_local)) + &
             nabt(3,1) * (dg_frag%phi_frag(ix0 + 3, iy0, iz0, jo, i_local) - &
                           dg_frag%phi_frag(ix0 - 3, iy0, iz0, jo, i_local)) + &
             nabt(4,1) * (dg_frag%phi_frag(ix0 + 4, iy0, iz0, jo, i_local) - &
                           dg_frag%phi_frag(ix0 - 4, iy0, iz0, jo, i_local))

    grad_y = nabt(1,2) * (dg_frag%phi_frag(ix0, iy0 + 1, iz0, jo, i_local) - &
                           dg_frag%phi_frag(ix0, iy0 - 1, iz0, jo, i_local)) + &
             nabt(2,2) * (dg_frag%phi_frag(ix0, iy0 + 2, iz0, jo, i_local) - &
                           dg_frag%phi_frag(ix0, iy0 - 2, iz0, jo, i_local)) + &
             nabt(3,2) * (dg_frag%phi_frag(ix0, iy0 + 3, iz0, jo, i_local) - &
                           dg_frag%phi_frag(ix0, iy0 - 3, iz0, jo, i_local)) + &
             nabt(4,2) * (dg_frag%phi_frag(ix0, iy0 + 4, iz0, jo, i_local) - &
                           dg_frag%phi_frag(ix0, iy0 - 4, iz0, jo, i_local))

    grad_z = nabt(1,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 1, jo, i_local) - &
                           dg_frag%phi_frag(ix0, iy0, iz0 - 1, jo, i_local)) + &
             nabt(2,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 2, jo, i_local) - &
                           dg_frag%phi_frag(ix0, iy0, iz0 - 2, jo, i_local)) + &
             nabt(3,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 3, jo, i_local) - &
                           dg_frag%phi_frag(ix0, iy0, iz0 - 3, jo, i_local)) + &
             nabt(4,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 4, jo, i_local) - &
                           dg_frag%phi_frag(ix0, iy0, iz0 - 4, jo, i_local))

  end subroutine apply_gradient_at_phi_box_point

  subroutine apply_kinetic_at_phi_box_point(dg_frag, i_local, jo, mg, stencil, phi_idx, t_val)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer,                intent(in) :: i_local, jo
    type(s_rgrid),          intent(in) :: mg
    type(s_stencil),        intent(in) :: stencil
    integer,                intent(in) :: phi_idx(3)
    real(8),                intent(out) :: t_val

    integer :: gx, gy, gz
    integer :: phi_lb1, phi_lb2, phi_lb3, phi_ub1, phi_ub2, phi_ub3
    integer :: ix0, iy0, iz0
    real(8) :: lap0, lapt(4,3), v

    lap0 = stencil%coef_lap0
    lapt = stencil%coef_lap
    phi_lb1 = lbound(dg_frag%phi_frag, 1)
    phi_lb2 = lbound(dg_frag%phi_frag, 2)
    phi_lb3 = lbound(dg_frag%phi_frag, 3)
    phi_ub1 = ubound(dg_frag%phi_frag, 1)
    phi_ub2 = ubound(dg_frag%phi_frag, 2)
    phi_ub3 = ubound(dg_frag%phi_frag, 3)

    gx = modulo(phi_idx(1) - 1, mg%num(1)) + 1
    gy = modulo(phi_idx(2) - 1, mg%num(2)) + 1
    gz = modulo(phi_idx(3) - 1, mg%num(3)) + 1
    ix0 = gx
    iy0 = gy
    iz0 = gz
    if (ix0 - 4 < phi_lb1 .or. ix0 + 4 > phi_ub1 .or. &
        iy0 - 4 < phi_lb2 .or. iy0 + 4 > phi_ub2 .or. &
        iz0 - 4 < phi_lb3 .or. iz0 + 4 > phi_ub3) then
      t_val = 0.0d0
      return
    end if

    v = lapt(1,1) * (dg_frag%phi_frag(ix0 + 1, iy0, iz0, jo, i_local) + dg_frag%phi_frag(ix0 - 1, iy0, iz0, jo, i_local)) + &
        lapt(2,1) * (dg_frag%phi_frag(ix0 + 2, iy0, iz0, jo, i_local) + dg_frag%phi_frag(ix0 - 2, iy0, iz0, jo, i_local)) + &
        lapt(3,1) * (dg_frag%phi_frag(ix0 + 3, iy0, iz0, jo, i_local) + dg_frag%phi_frag(ix0 - 3, iy0, iz0, jo, i_local)) + &
        lapt(4,1) * (dg_frag%phi_frag(ix0 + 4, iy0, iz0, jo, i_local) + dg_frag%phi_frag(ix0 - 4, iy0, iz0, jo, i_local))
    v = v + &
        lapt(1,2) * (dg_frag%phi_frag(ix0, iy0 + 1, iz0, jo, i_local) + dg_frag%phi_frag(ix0, iy0 - 1, iz0, jo, i_local)) + &
        lapt(2,2) * (dg_frag%phi_frag(ix0, iy0 + 2, iz0, jo, i_local) + dg_frag%phi_frag(ix0, iy0 - 2, iz0, jo, i_local)) + &
        lapt(3,2) * (dg_frag%phi_frag(ix0, iy0 + 3, iz0, jo, i_local) + dg_frag%phi_frag(ix0, iy0 - 3, iz0, jo, i_local)) + &
        lapt(4,2) * (dg_frag%phi_frag(ix0, iy0 + 4, iz0, jo, i_local) + dg_frag%phi_frag(ix0, iy0 - 4, iz0, jo, i_local))
    v = v + &
        lapt(1,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 1, jo, i_local) + dg_frag%phi_frag(ix0, iy0, iz0 - 1, jo, i_local)) + &
        lapt(2,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 2, jo, i_local) + dg_frag%phi_frag(ix0, iy0, iz0 - 2, jo, i_local)) + &
        lapt(3,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 3, jo, i_local) + dg_frag%phi_frag(ix0, iy0, iz0 - 3, jo, i_local)) + &
        lapt(4,3) * (dg_frag%phi_frag(ix0, iy0, iz0 + 4, jo, i_local) + dg_frag%phi_frag(ix0, iy0, iz0 - 4, jo, i_local))

    t_val = lap0 * dg_frag%phi_frag(ix0, iy0, iz0, jo, i_local) - 0.5d0 * v
  end subroutine apply_kinetic_at_phi_box_point

  subroutine apply_kinetic_and_hamiltonian_at_phi_box_point(dg_frag, i_local, jo, mg, stencil, V_total, phi_idx, t_val, h_val)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer,                intent(in) :: i_local, jo
    type(s_rgrid),          intent(in) :: mg
    type(s_stencil),        intent(in) :: stencil
    real(8),                intent(in) :: V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    integer,                intent(in) :: phi_idx(3)
    real(8),                intent(out) :: t_val, h_val

    integer :: gx, gy, gz

    call apply_kinetic_at_phi_box_point(dg_frag, i_local, jo, mg, stencil, phi_idx, t_val)
    h_val = t_val
    gx = modulo(phi_idx(1) - 1, mg%num(1)) + 1
    gy = modulo(phi_idx(2) - 1, mg%num(2)) + 1
    gz = modulo(phi_idx(3) - 1, mg%num(3)) + 1
    if (gx < mg%is(1) .or. gx > mg%ie(1) .or. &
        gy < mg%is(2) .or. gy > mg%ie(2) .or. &
        gz < mg%is(3) .or. gz > mg%ie(3)) return
    if (phi_idx(1) < lbound(dg_frag%phi_frag, 1) .or. phi_idx(1) > ubound(dg_frag%phi_frag, 1) .or. &
        phi_idx(2) < lbound(dg_frag%phi_frag, 2) .or. phi_idx(2) > ubound(dg_frag%phi_frag, 2) .or. &
        phi_idx(3) < lbound(dg_frag%phi_frag, 3) .or. phi_idx(3) > ubound(dg_frag%phi_frag, 3)) return
    h_val = h_val + V_total(gx, gy, gz) * dg_frag%phi_frag(phi_idx(1), phi_idx(2), phi_idx(3), jo, i_local)
  end subroutine apply_kinetic_and_hamiltonian_at_phi_box_point

  subroutine apply_hamiltonian_at_phi_box_point(dg_frag, i_local, jo, mg, stencil, V_total, phi_idx, h_val)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer,                intent(in) :: i_local, jo
    type(s_rgrid),          intent(in) :: mg
    type(s_stencil),        intent(in) :: stencil
    real(8),                intent(in) :: V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    integer,                intent(in) :: phi_idx(3)
    real(8),                intent(out) :: h_val

    real(8) :: t_dummy

    call apply_kinetic_and_hamiltonian_at_phi_box_point(dg_frag, i_local, jo, mg, stencil, V_total, phi_idx, t_dummy, h_val)
  end subroutine apply_hamiltonian_at_phi_box_point

  !=======================================================================
  ! Add non-local pseudopotential contribution to Hamiltonian matrix
  !
  ! Calculates <φ_i|V_NL|φ_j> = Σ_ilma <φ_i|proj_ilma> V_ilma <proj_ilma|φ_j>
  ! where proj_ilma are the pseudopotential projector functions
  !
  ! NUMERICAL ACCURACY: Store unnormalized overlaps, apply rinv_uvu once
  ! This prevents rinv_uvu^2 error amplification and follows SALMON convention
  !=======================================================================
  subroutine add_nonlocal_pp_matrix(dg_frag, mg, ppg, nspin, hvol)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_rgrid),          intent(in) :: mg
    type(s_pp_grid),        intent(in) :: ppg
    integer,                intent(in)    :: nspin
    real(8),                intent(in)    :: hvol
    
    stop "DG-Fragment RT: add_nonlocal_pp_matrix fallback path is retired; use ensure_nonlocal_pp_matrix_A"
    
  end subroutine add_nonlocal_pp_matrix

  !=======================================================================
  ! Calculate momentum matrix elements in fragment basis (velocity gauge)
  !=======================================================================
  subroutine calculate_momentum_matrix(dg_frag, system, mg, stencil)
    use structures
    use communication, only: comm_is_root, comm_summation, comm_get_groupinfo, comm_get_min, comm_get_max
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    
    integer :: ifrag, i_local, ispin, io, jo, idir, nbf, jo_progress_stride
    integer :: ix, iy, iz, is(3), ie(3), i_halo, jfrag, n_basis_halo, ig_row, ig_col, ig_i, ig_j, l(3), d(3)
    integer :: lx, ly, lz, gx, gy, gz, iorg(3), ndom(3), ndom_grad(3), loc_s(3), loc_e(3), phi_loc_s(3), phi_loc_e(3), halo_s(3), halo_e(3)
    integer :: lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi
    integer :: halo_send_idx(3), halo_recv_idx(3)
    integer :: phi_lb1, phi_ub1, phi_lb2, phi_ub2, phi_lb3, phi_ub3
    integer :: grad_lb1, grad_ub1, grad_lb2, grad_ub2, grad_lb3, grad_ub3
    integer :: iblk, iblk_rev, iblk_self, ii, jj, mat_size, ni, nj, ndiag
    integer :: npts_local, npts_halo, ipt, nx_local, ny_local
    integer :: n_basis_halo_max, npts_halo_max
    logical :: log_frag_progress, has_overlap
    logical :: has_halo_work, use_row_local_momentum
    logical :: enable_hmat_timing_trace
    logical, save :: hmat_timing_trace_initialized = .false.
    logical, save :: hmat_timing_trace_cached = .false.
    logical :: enable_momentum_probe
    real(8) :: hvol, integral
    real(8) :: momentum_gb
    real(8) :: max_p, pavg
    real(8) :: t0, t1
    real(8) :: time_halo_exchange, time_self_integral, time_halo_integral
    real(8) :: time_grad_total
    real(8) :: time_block_reduce, time_antisym
    real(8) :: time_reduce_pack, time_reduce_comm, time_reduce_unpack
    real(8) :: frag_grad_start, frag_self_start, frag_halo_start
    real(8), allocatable :: grad_phi(:,:,:,:)  ! gradient of basis function (x,y,z components, fragment-local)
    real(8), allocatable :: grad_local_2d(:,:), phi_local_2d(:,:), self_proj(:,:)
    real(8), allocatable :: grad_halo_2d(:,:), halo_buf_2d(:,:), halo_proj(:,:)
    character(len=32) :: env_hmat_trace, env_momentum_probe
    integer :: env_hmat_len, env_hmat_stat
    
    if (.not. dg_frag%has_real_space_basis) return
    enable_momentum_probe = .false.
    env_momentum_probe = ''
    call get_environment_variable('SALMON_DG_MOMENTUM_PROBE', env_momentum_probe, &
      length=env_hmat_len, status=env_hmat_stat)
    if (env_hmat_stat /= 0 .or. env_hmat_len <= 0) then
      call get_environment_variable('SALMON_DG_TRANSITION_PROBE', env_momentum_probe, &
        length=env_hmat_len, status=env_hmat_stat)
    end if
    if (env_hmat_stat == 0 .and. env_hmat_len > 0) then
      if (env_momentum_probe(1:1) == '1' .or. env_momentum_probe(1:1) == 'y' .or. &
          env_momentum_probe(1:1) == 'Y' .or. env_momentum_probe(1:1) == 't' .or. &
          env_momentum_probe(1:1) == 'T') then
        enable_momentum_probe = .true.
      end if
    end if

    time_halo_exchange = 0.0d0
    time_self_integral = 0.0d0
    time_halo_integral = 0.0d0
    time_grad_total = 0.0d0
    time_block_reduce = 0.0d0
    time_antisym = 0.0d0
    time_reduce_pack = 0.0d0
    time_reduce_comm = 0.0d0
    time_reduce_unpack = 0.0d0
    use_row_local_momentum = momentum_blocks_use_row_local_storage(dg_frag)
    if (dg_frag%use_buffered_basis) use_row_local_momentum = .true.

    if (.not. hmat_timing_trace_initialized) then
      env_hmat_trace = ''
      env_hmat_len = 0
      env_hmat_stat = 0
      call get_environment_variable('SALMON_DG_HMAT_TIMING_TRACE', env_hmat_trace, &
        length=env_hmat_len, status=env_hmat_stat)
      if (env_hmat_stat == 0 .and. env_hmat_len > 0) then
        if (env_hmat_trace(1:1) == '1' .or. env_hmat_trace(1:1) == 'y' .or. &
            env_hmat_trace(1:1) == 'Y' .or. env_hmat_trace(1:1) == 't' .or. &
            env_hmat_trace(1:1) == 'T') then
          hmat_timing_trace_cached = .true.
        end if
      end if
      hmat_timing_trace_initialized = .true.
    end if
    enable_hmat_timing_trace = hmat_timing_trace_cached
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "        Computing transition moments: <φ_i|∇|φ_j>"
      flush(6)
    end if
    if (enable_momentum_probe) then
      write(*,'(1x,a,i0,a)') "        momentum-probe: rank=", dg_frag%id, " stage=enter"
      flush(6)
    end if
    momentum_gb = real(3_8 * int(dg_frag%n_mat_max, kind=8) * int(dg_frag%n_mat_max, kind=8) * &
      int(dg_frag%nspin, kind=8) * 8_8, 8) / 1.0d9
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,f10.3,a)') "        n_mat_max=", dg_frag%n_mat_max, &
        " momentum_mat GB=", momentum_gb, " per rank"
      flush(6)
    end if
    
    ! Allocate momentum matrix: (3 directions, n_mat_max x n_mat_max, nspin)
    ! Momentum matrix elements for vector potential coupling: p_ij = <phi_i|p|phi_j>
    ! In velocity gauge: H(t) = H_0 - i*A(t)·∇ + A(t)^2/2
    ! The A·p term couples to momentum matrix elements
    ! The A^2/2 term is diagonal (diamagnetic contribution)
    if (allocated(dg_frag%momentum_mat)) deallocate(dg_frag%momentum_mat)
    call init_momentum_blocks(dg_frag)
    is = mg%is
    ie = mg%ie
    hvol = system%hvol
    
    ! Buffered basis already carries its own periodic support, so neighbor halo
    ! exchange is only needed for the raw fragment-domain path.
    if (.not. dg_frag%use_buffered_basis) then
      call cpu_time(t0)
      call exchange_phi_frag_halo(dg_frag)
      call cpu_time(t1)
      time_halo_exchange = time_halo_exchange + (t1 - t0)
    end if
    if (enable_momentum_probe) then
      write(*,'(1x,a,i0,a)') "        momentum-probe: rank=", dg_frag%id, " stage=after-halo"
      flush(6)
    end if
    
    ! Loop over spin
    do ispin = 1, system%nspin
      ! Loop over local fragments
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        log_frag_progress = dg_frag%is_frag_root .and. ifrag == dg_frag%ifrag_start
        frag_grad_start = time_grad_total
        frag_self_start = time_self_integral
        frag_halo_start = time_halo_integral
        call get_fragment_basis_owned_range(dg_frag, ifrag, mg, iorg, ndom, loc_s, loc_e, has_overlap)
        if (enable_momentum_probe) then
          write(*,'(1x,a,i0,a,i0,a,l1,a,3(i0,1x),a,3(i0,1x))') "        momentum-probe: rank=", dg_frag%id, &
            " ifrag=", ifrag, " has_overlap=", has_overlap, " loc_s=", loc_s(1), loc_s(2), loc_s(3), &
            " loc_e=", loc_e(1), loc_e(2), loc_e(3)
          flush(6)
        end if
        if (.not. has_overlap) cycle
        
        ! Cache the local basis matrix once per fragment; it does not depend on jo.
        if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) then
          write(*,'(1x,a,i0,a,i0)') "DG-Fragment RT invalid i_local=", i_local, " phi_frag dim5=", size(dg_frag%phi_frag, 5)
          stop "DG-Fragment RT: invalid fragment-local basis index"
        end if
        phi_lb1 = lbound(dg_frag%phi_frag, 1)
        phi_ub1 = ubound(dg_frag%phi_frag, 1)
        phi_lb2 = lbound(dg_frag%phi_frag, 2)
        phi_ub2 = ubound(dg_frag%phi_frag, 2)
        phi_lb3 = lbound(dg_frag%phi_frag, 3)
        phi_ub3 = ubound(dg_frag%phi_frag, 3)
        phi_loc_s(:) = iorg(:) + loc_s(:) - 1
        phi_loc_e(:) = iorg(:) + loc_e(:) - 1
        if (phi_loc_s(1) < phi_lb1 .or. phi_loc_e(1) > phi_ub1 .or. &
            phi_loc_s(2) < phi_lb2 .or. phi_loc_e(2) > phi_ub2 .or. &
            phi_loc_s(3) < phi_lb3 .or. phi_loc_e(3) > phi_ub3) then
          write(*,'(1x,a,1x,3(i0,1x),a,1x,3(i0,1x),a,1x,3(i0,1x),a,1x,3(i0,1x))') &
            "DG-Fragment RT momentum phi_frag local range out of bounds: phi_loc_s=", &
            phi_loc_s(1), phi_loc_s(2), phi_loc_s(3), "phi_loc_e=", phi_loc_e(1), phi_loc_e(2), phi_loc_e(3), &
            "lb=", phi_lb1, phi_lb2, phi_lb3, "ub=", phi_ub1, phi_ub2, phi_ub3
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
            "        momentum context: rank=", dg_frag%id, " id_frag=", dg_frag%id_frag, " ifrag=", ifrag, &
            " i_local=", i_local, " mg_is=", mg%is(1), mg%is(2), mg%is(3), " mg_ie=", mg%ie(1), mg%ie(2), mg%ie(3), &
            " loc_s=", loc_s(1), loc_s(2), loc_s(3), " loc_e=", loc_e(1), loc_e(2), loc_e(3)
          call flush(6)
          stop "DG-Fragment RT: momentum phi_frag local range out of bounds"
        end if
        nbf = dg_frag%n_basis(ifrag, ispin)
        jo_progress_stride = max(1, nbf / 4)
        npts_local = (loc_e(1) - loc_s(1) + 1) * (loc_e(2) - loc_s(2) + 1) * (loc_e(3) - loc_s(3) + 1)
        lx_lo = phi_loc_s(1)
        lx_hi = phi_loc_e(1)
        ly_lo = phi_loc_s(2)
        ly_hi = phi_loc_e(2)
        lz_lo = phi_loc_s(3)
        lz_hi = phi_loc_e(3)
        nx_local = lx_hi - lx_lo + 1
        ny_local = ly_hi - ly_lo + 1
        allocate(phi_local_2d(npts_local, nbf), grad_local_2d(npts_local, 3), self_proj(nbf, 3))
        ndom_grad(:) = ndom(:)
        if (dg_frag%use_buffered_basis) then
          ndom_grad(:) = dg_frag%basis_support_hi(:, ifrag) - dg_frag%basis_support_lo(:, ifrag) + 1
        end if
        allocate(grad_phi(1:ndom_grad(1), 1:ndom_grad(2), 1:ndom_grad(3), 3))
        if (enable_momentum_probe) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        momentum-probe: rank=", dg_frag%id, &
            " ifrag=", ifrag, " nbf=", nbf, " npts_local=", npts_local
          flush(6)
        end if

        if (nbf > size(dg_frag%phi_frag, 4)) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "DG-Fragment RT invalid n_basis=", nbf, &
            " phi_frag dim4=", size(dg_frag%phi_frag, 4), " ifrag=", ifrag, " ispin=", ispin
          stop "DG-Fragment RT: invalid basis-function count"
        end if
        !$omp parallel do private(io,lz,ly,lx,ipt) schedule(static)
        do io = 1, nbf
          do lz = lz_lo, lz_hi
            do ly = ly_lo, ly_hi
              !$omp simd private(ipt)
              do lx = lx_lo, lx_hi
                ipt = ((lz - lz_lo) * ny_local + (ly - ly_lo)) * nx_local + (lx - lx_lo) + 1
                phi_local_2d(ipt, io) = dg_frag%phi_frag(lx, ly, lz, io, i_local)
              end do
            end do
          end do
        end do
        !$omp end parallel do
        iblk_self = find_momentum_block(dg_frag, ifrag, ifrag)

        n_basis_halo_max = 0
        npts_halo_max = 0
        do i_halo = 1, dg_frag%n_halo
          if (dg_frag%halo(i_halo)%ifrag_dst /= ifrag) cycle
          jfrag = dg_frag%halo(i_halo)%ifrag_src
          if (jfrag < 1) cycle
          n_basis_halo_max = max(n_basis_halo_max, dg_frag%n_basis(jfrag, ispin))
          l = dg_frag%halo(i_halo)%length
          npts_halo_max = max(npts_halo_max, l(1) * l(2) * l(3))
        end do
        has_halo_work = (n_basis_halo_max > 0 .and. npts_halo_max > 0)
        if (has_halo_work) then
          allocate(grad_halo_2d(npts_halo_max, 3), halo_buf_2d(npts_halo_max, n_basis_halo_max), halo_proj(n_basis_halo_max, 3))
        end if

        ! Loop over basis functions in fragment j (ket side)
        ! Keep this loop serial to avoid per-thread duplication of large grad_phi buffers.
        ! Parallelism is still provided inside apply_gradient_to_basis and SIMD in accumulations.
        do jo = 1, nbf
          if (enable_momentum_probe .and. (jo == 1 .or. jo == nbf)) then
            write(*,'(1x,a,i0,a,i0,a,i0)') "        momentum-probe: rank=", dg_frag%id, &
              " ifrag=", ifrag, " jo=", jo
            flush(6)
          end if
          if (enable_momentum_probe .and. jo == 1) then
            write(*,'(1x,a,i0,a,i0,a)') "        momentum-probe: rank=", dg_frag%id, &
              " ifrag=", ifrag, " stage=before-apply-gradient"
            flush(6)
          end if
          call cpu_time(t0)
          call apply_gradient_to_basis_ops_local_2d(dg_frag, i_local, jo, mg, stencil, loc_s, loc_e, grad_phi, grad_local_2d)
          call cpu_time(t1)
          time_grad_total = time_grad_total + (t1 - t0)
          if (enable_momentum_probe .and. jo == 1) then
            write(*,'(1x,a,i0,a,i0,a)') "        momentum-probe: rank=", dg_frag%id, &
              " ifrag=", ifrag, " stage=after-apply-gradient"
            flush(6)
          end if

          grad_lb1 = lbound(grad_phi, 1)
          grad_ub1 = ubound(grad_phi, 1)
          grad_lb2 = lbound(grad_phi, 2)
          grad_ub2 = ubound(grad_phi, 2)
          grad_lb3 = lbound(grad_phi, 3)
          grad_ub3 = ubound(grad_phi, 3)
          if (loc_s(1) < grad_lb1 .or. loc_e(1) > grad_ub1 .or. &
              loc_s(2) < grad_lb2 .or. loc_e(2) > grad_ub2 .or. &
              loc_s(3) < grad_lb3 .or. loc_e(3) > grad_ub3) then
            write(*,'(1x,a,1x,3(i0,1x),a,1x,3(i0,1x),a,1x,3(i0,1x),a,1x,3(i0,1x))') &
              "DG-Fragment RT momentum grad_phi local range out of bounds: loc_s=", &
              loc_s(1), loc_s(2), loc_s(3), "loc_e=", loc_e(1), loc_e(2), loc_e(3), &
              "lb=", grad_lb1, grad_lb2, grad_lb3, "ub=", grad_ub1, grad_ub2, grad_ub3
            stop "DG-Fragment RT: momentum grad_phi local range out of bounds"
          end if

          ig_j = dg_frag%index_basis(jo, ifrag, ispin)
          if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) then
            cycle
          end if
          call cpu_time(t0)
          call dgemm('T', 'N', nbf, 3, npts_local, hvol, phi_local_2d, npts_local, &
            grad_local_2d, npts_local, 0.0d0, self_proj, nbf)
          call cpu_time(t1)
          time_self_integral = time_self_integral + (t1 - t0)
          if (enable_momentum_probe .and. jo == 1) then
            write(*,'(1x,a,i0,a,i0,a)') "        momentum-probe: rank=", dg_frag%id, &
              " ifrag=", ifrag, " stage=after-self-dgemm"
            flush(6)
          end if

          do io = 1, nbf
            ig_i = dg_frag%index_basis(io, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            do idir = 1, 3
              integral = self_proj(io, idir)
              if (iblk_self > 0) then
                if (io <= size(dg_frag%momentum_blocks(iblk_self)%val, 2) .and. &
                    jo <= size(dg_frag%momentum_blocks(iblk_self)%val, 3)) then
                  dg_frag%momentum_blocks(iblk_self)%val(idir, io, jo, ispin) = integral
                end if
              end if
            end do
          end do

          ig_col = ig_j
          if (ig_col >= 1 .and. ig_col <= dg_frag%n_mat_max) then
            if (enable_momentum_probe .and. jo == 1) then
              write(*,'(1x,a,i0,a,i0,a,i0)') "        momentum-probe: rank=", dg_frag%id, &
                " ifrag=", ifrag, " stage=before-halo-loop ig_col=", ig_col
              flush(6)
            end if
            do i_halo = 1, dg_frag%n_halo
              if (dg_frag%use_buffered_basis) cycle
              if (enable_momentum_probe .and. jo == 1) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,3(i0,1x))') "        momentum-probe: rank=", dg_frag%id, &
                  " ifrag=", ifrag, " i_halo=", i_halo, " src=", dg_frag%halo(i_halo)%ifrag_src, &
                  " dst=", dg_frag%halo(i_halo)%ifrag_dst, " len=", dg_frag%halo(i_halo)%length(1), &
                  dg_frag%halo(i_halo)%length(2), dg_frag%halo(i_halo)%length(3)
                flush(6)
              end if
              if (dg_frag%halo(i_halo)%ifrag_dst /= ifrag) cycle
              jfrag = dg_frag%halo(i_halo)%ifrag_src
              n_basis_halo = dg_frag%n_basis(jfrag, ispin)
              l = dg_frag%halo(i_halo)%length
              if (enable_momentum_probe .and. jo == 1) then
                write(*,'(1x,a,i0,a,i0,a,l1,a,l1,a,i0)') "        momentum-probe: rank=", dg_frag%id, &
                  " ifrag=", ifrag, " buf_recv_alloc=", allocated(dg_frag%halo(i_halo)%buf_recv), &
                  " buf_send_alloc=", allocated(dg_frag%halo(i_halo)%buf_send), " n_basis_halo=", n_basis_halo
                flush(6)
              end if
              if (size(dg_frag%halo(i_halo)%buf_recv, 5) < 1) then
                write(*,*) "[FATAL] momentum halo buf dim5 invalid: rank=", dg_frag%id, " i_halo=", i_halo
                stop 1
              end if
              if (enable_momentum_probe .and. jo == 1) then
                write(*,'(1x,a,i0,a,5(i0,1x))') "        momentum-probe: rank=", dg_frag%id, &
                  " buf_recv_shape=", size(dg_frag%halo(i_halo)%buf_recv, 1), size(dg_frag%halo(i_halo)%buf_recv, 2), &
                  size(dg_frag%halo(i_halo)%buf_recv, 3), size(dg_frag%halo(i_halo)%buf_recv, 4), &
                  size(dg_frag%halo(i_halo)%buf_recv, 5)
                flush(6)
              end if
              if (n_basis_halo > size(dg_frag%halo(i_halo)%buf_recv, 4)) then
                write(*,*) "[FATAL] momentum halo basis exceeds buf dim4: rank=", dg_frag%id, &
                  " i_halo=", i_halo, " n_basis_halo=", n_basis_halo, " buf_dim4=", size(dg_frag%halo(i_halo)%buf_recv, 4)
                stop 1
              end if
              if (n_basis_halo > size(dg_frag%index_basis, 1)) then
                write(*,*) "[FATAL] momentum halo basis exceeds index_basis dim1: rank=", dg_frag%id, &
                  " i_halo=", i_halo, " jfrag=", jfrag, " n_basis_halo=", n_basis_halo, &
                  " index_basis_dim1=", size(dg_frag%index_basis, 1)
                stop 1
              end if
              npts_halo = l(1) * l(2) * l(3)
              if (npts_halo <= 0 .or. n_basis_halo <= 0) cycle
              if (enable_momentum_probe .and. jo == 1) then
                write(*,'(1x,a,i0,a,i0)') "        momentum-probe: rank=", dg_frag%id, &
                  " npts_halo=", npts_halo
                flush(6)
              end if

              if (.not. has_halo_work) cycle
              if (enable_momentum_probe .and. jo == 1) then
                write(*,'(1x,a,i0,a)') "        momentum-probe: rank=", dg_frag%id, " stage=after-halo-alloc"
                flush(6)
              end if

              ipt = 0
              do iz = 1, l(3)
                do iy = 1, l(2)
                  do ix = 1, l(1)
                    if (enable_momentum_probe .and. jo == 1 .and. ix == 1 .and. iy == 1 .and. iz == 1) then
                      write(*,'(1x,a,i0,a)') "        momentum-probe: rank=", dg_frag%id, " stage=before-first-halo-index"
                      flush(6)
                    end if
                    call get_halo_block_point_indices(dg_frag%halo(i_halo), ix, iy, iz, halo_send_idx, halo_recv_idx)
                    if (enable_momentum_probe .and. jo == 1 .and. ix == 1 .and. iy == 1 .and. iz == 1) then
                      write(*,'(1x,a,i0,a,3(i0,1x),a,3(i0,1x))') "        momentum-probe: rank=", dg_frag%id, &
                        " first-halo send_idx=", halo_send_idx(1), halo_send_idx(2), halo_send_idx(3), &
                        " recv_idx=", halo_recv_idx(1), halo_recv_idx(2), halo_recv_idx(3)
                      flush(6)
                    end if
                    ipt = ipt + 1
                    ! For stencil evaluation use source-side indices; recv-side halo points can
                    ! sit at the outer ghost edge where +/-4 access may go out of bounds.
                    call apply_gradient_at_phi_box_point(dg_frag, i_local, jo, mg, stencil, halo_send_idx, &
                      grad_halo_2d(ipt, 1), grad_halo_2d(ipt, 2), grad_halo_2d(ipt, 3))
                  end do
                end do
              end do

              !$omp parallel do private(io,iz,iy,ix,ipt) schedule(static)
              do io = 1, n_basis_halo
                do iz = 1, l(3)
                  do iy = 1, l(2)
                    !$omp simd private(ipt)
                    do ix = 1, l(1)
                      ipt = ((iz - 1) * l(2) + (iy - 1)) * l(1) + ix
                      halo_buf_2d(ipt, io) = dg_frag%halo(i_halo)%buf_recv(ix, iy, iz, io, 1)
                    end do
                  end do
                end do
              end do
              !$omp end parallel do

              call cpu_time(t0)
              call dgemm('T', 'N', n_basis_halo, 3, npts_halo, hvol, halo_buf_2d, npts_halo_max, &
                grad_halo_2d, npts_halo_max, 0.0d0, halo_proj, n_basis_halo_max)
              call cpu_time(t1)
              time_halo_integral = time_halo_integral + (t1 - t0)

              iblk = find_momentum_block(dg_frag, jfrag, ifrag)
              iblk_rev = find_momentum_block(dg_frag, ifrag, jfrag)
              if (iblk > 0 .and. iblk_rev > 0) then
                !$omp parallel do private(io,ig_row,idir,integral) schedule(static)
                do io = 1, n_basis_halo
                  ig_row = dg_frag%index_basis(io, jfrag, ispin)
                  if (ig_row < 1 .or. ig_row > dg_frag%n_mat_max) cycle
                  if (io <= size(dg_frag%momentum_blocks(iblk)%val, 2) .and. &
                      jo <= size(dg_frag%momentum_blocks(iblk)%val, 3) .and. &
                      jo <= size(dg_frag%momentum_blocks(iblk_rev)%val, 2) .and. &
                      io <= size(dg_frag%momentum_blocks(iblk_rev)%val, 3)) then
                    do idir = 1, 3
                      integral = halo_proj(io, idir)
                      dg_frag%momentum_blocks(iblk)%val(idir, io, jo, ispin) = &
                        dg_frag%momentum_blocks(iblk)%val(idir, io, jo, ispin) + 0.5d0 * integral
                      dg_frag%momentum_blocks(iblk_rev)%val(idir, jo, io, ispin) = &
                        dg_frag%momentum_blocks(iblk_rev)%val(idir, jo, io, ispin) - 0.5d0 * integral
                    end do
                  end if
                end do
                !$omp end parallel do
              else if (iblk > 0) then
                !$omp parallel do private(io,ig_row,idir,integral) schedule(static)
                do io = 1, n_basis_halo
                  ig_row = dg_frag%index_basis(io, jfrag, ispin)
                  if (ig_row < 1 .or. ig_row > dg_frag%n_mat_max) cycle
                  if (io <= size(dg_frag%momentum_blocks(iblk)%val, 2) .and. &
                      jo <= size(dg_frag%momentum_blocks(iblk)%val, 3)) then
                    do idir = 1, 3
                      integral = halo_proj(io, idir)
                      dg_frag%momentum_blocks(iblk)%val(idir, io, jo, ispin) = &
                        dg_frag%momentum_blocks(iblk)%val(idir, io, jo, ispin) + 0.5d0 * integral
                    end do
                  end if
                end do
                !$omp end parallel do
              else if (iblk_rev > 0) then
                !$omp parallel do private(io,ig_row,idir,integral) schedule(static)
                do io = 1, n_basis_halo
                  ig_row = dg_frag%index_basis(io, jfrag, ispin)
                  if (ig_row < 1 .or. ig_row > dg_frag%n_mat_max) cycle
                  if (jo <= size(dg_frag%momentum_blocks(iblk_rev)%val, 2) .and. &
                      io <= size(dg_frag%momentum_blocks(iblk_rev)%val, 3)) then
                    do idir = 1, 3
                      integral = halo_proj(io, idir)
                      dg_frag%momentum_blocks(iblk_rev)%val(idir, jo, io, ispin) = &
                        dg_frag%momentum_blocks(iblk_rev)%val(idir, jo, io, ispin) - 0.5d0 * integral
                    end do
                  end if
                end do
                !$omp end parallel do
              end if

            end do
            if (enable_momentum_probe .and. jo == 1) then
              write(*,'(1x,a,i0,a,i0,a)') "        momentum-probe: rank=", dg_frag%id, &
                " ifrag=", ifrag, " stage=after-halo-loop"
              flush(6)
            end if
          end if

        end do  ! jo
        if (enable_momentum_probe) then
          write(*,'(1x,a,i0,a,i0,a)') "        momentum-probe: rank=", dg_frag%id, &
            " ifrag=", ifrag, " stage=after-jo-loop"
          flush(6)
        end if
        if (allocated(grad_phi)) deallocate(grad_phi)
        if (allocated(grad_halo_2d)) deallocate(grad_halo_2d)
        if (allocated(halo_buf_2d)) deallocate(halo_buf_2d)
        if (allocated(halo_proj)) deallocate(halo_proj)
        deallocate(phi_local_2d, grad_local_2d, self_proj)
        if (enable_momentum_probe) then
          write(*,'(1x,a,i0,a,i0,a)') "        momentum-probe: rank=", dg_frag%id, &
            " ifrag=", ifrag, " stage=after-frag-dealloc"
          flush(6)
        end if
      end do  ! ifrag
    end do  ! ispin
    if (enable_momentum_probe) then
      write(*,'(1x,a,i0,a,l1)') "        momentum-probe: rank=", dg_frag%id, &
        " stage=before-reduce use_row_local=", use_row_local_momentum
      flush(6)
    end if
    
    ! MPI aggregation of fragment-neighbor momentum blocks.
    if (.not. use_row_local_momentum) then
    block
      real(8), allocatable :: send_flat(:), recv_flat(:)
      real(8) :: t_reduce_start, t_reduce_end
      integer :: nrow, ncol
      integer :: total_size, offset_flat, offset_base, chunk_size, chunk_begin, chunk_count
      integer :: total_size_min, total_size_max, nblk_min, nblk_max, ifrag_chk
      real(8) :: meta_sig_blocks, meta_sig_basis
      real(8) :: meta_sig_blocks_min(1), meta_sig_blocks_max(1)
      real(8) :: meta_sig_basis_min(1), meta_sig_basis_max(1)
      call cpu_time(t_reduce_start)
      total_size = 0
      do iblk = 1, dg_frag%n_momentum_blocks
        nrow = dg_frag%momentum_blocks(iblk)%nrow_max
        ncol = dg_frag%momentum_blocks(iblk)%ncol_max
        total_size = total_size + 3 * nrow * ncol * dg_frag%nspin
      end do
      nblk_max = dg_frag%n_momentum_blocks
      call comm_get_max(nblk_max, dg_frag%icomm)
      nblk_min = -dg_frag%n_momentum_blocks
      call comm_get_max(nblk_min, dg_frag%icomm)
      nblk_min = -nblk_min
      if (nblk_min /= nblk_max) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] momentum reduce metadata mismatch: rank=", &
          dg_frag%id, " nblk=", dg_frag%n_momentum_blocks, " min=", nblk_min, " max=", nblk_max
        flush(6)
        stop "DG-Fragment RT: momentum reduce n_momentum_blocks mismatch across MPI ranks"
      end if

      meta_sig_blocks = 0.0d0
      do iblk = 1, dg_frag%n_momentum_blocks
        meta_sig_blocks = meta_sig_blocks + &
          dble(iblk) * 1.0d0 + &
          dble(dg_frag%momentum_blocks(iblk)%ifrag_row) * 1.0d3 + &
          dble(dg_frag%momentum_blocks(iblk)%ifrag_col) * 1.0d6 + &
          dble(dg_frag%momentum_blocks(iblk)%nrow_max) * 1.0d9 + &
          dble(dg_frag%momentum_blocks(iblk)%ncol_max) * 1.0d12
      end do
      meta_sig_blocks_min(1) = meta_sig_blocks
      meta_sig_blocks_max(1) = meta_sig_blocks
      call comm_get_min(meta_sig_blocks_min, meta_sig_blocks_min, 1, dg_frag%icomm)
      call comm_get_max(meta_sig_blocks_max, meta_sig_blocks_max, 1, dg_frag%icomm)

      meta_sig_basis = 0.0d0
      do ispin = 1, dg_frag%nspin
        do ifrag_chk = 1, dg_frag%n_frag
          meta_sig_basis = meta_sig_basis + &
            dble(ifrag_chk) * 1.0d0 + &
            dble(ispin) * 1.0d3 + &
            dble(dg_frag%n_basis(ifrag_chk, ispin)) * 1.0d6
        end do
      end do
      meta_sig_basis_min(1) = meta_sig_basis
      meta_sig_basis_max(1) = meta_sig_basis
      call comm_get_min(meta_sig_basis_min, meta_sig_basis_min, 1, dg_frag%icomm)
      call comm_get_max(meta_sig_basis_max, meta_sig_basis_max, 1, dg_frag%icomm)

      if (abs(meta_sig_blocks_max(1) - meta_sig_blocks_min(1)) > 1.0d-9 .or. &
          abs(meta_sig_basis_max(1) - meta_sig_basis_min(1)) > 1.0d-9) then
        write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
          "        [FATAL] momentum reduce metadata mismatch: rank=", dg_frag%id, &
          " block_min=", meta_sig_blocks_min(1), " block_max=", meta_sig_blocks_max(1), &
          " basis_min=", meta_sig_basis_min(1), " basis_max=", meta_sig_basis_max(1)
        flush(6)
        stop "DG-Fragment RT: momentum reduce metadata signature mismatch across MPI ranks"
      end if

      ! All ranks must call the same number of collective operations with the
      ! same count. If total_size diverges, chunked allreduce can deadlock.
      total_size_max = total_size
      call comm_get_max(total_size_max, dg_frag%icomm)
      total_size_min = -total_size
      call comm_get_max(total_size_min, dg_frag%icomm)
      total_size_min = -total_size_min
      if (total_size_min /= total_size_max) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] momentum reduce size mismatch: rank=", &
          dg_frag%id, " local=", total_size, " min=", total_size_min, " max=", total_size_max
        flush(6)
        stop "DG-Fragment RT: momentum reduce total_size mismatch across MPI ranks"
      end if

      if (total_size > 0) then
        allocate(send_flat(total_size), recv_flat(total_size))
        call cpu_time(t0)
        offset_flat = 1
        do iblk = 1, dg_frag%n_momentum_blocks
          nrow = dg_frag%momentum_blocks(iblk)%nrow_max
          ncol = dg_frag%momentum_blocks(iblk)%ncol_max
          do ispin = 1, dg_frag%nspin
            do jj = 1, ncol
              offset_base = offset_flat - 1
              !$omp simd
              do ii = 1, nrow
                send_flat(offset_base + (ii - 1) * 3 + 1) = dg_frag%momentum_blocks(iblk)%val(1, ii, jj, ispin)
                send_flat(offset_base + (ii - 1) * 3 + 2) = dg_frag%momentum_blocks(iblk)%val(2, ii, jj, ispin)
                send_flat(offset_base + (ii - 1) * 3 + 3) = dg_frag%momentum_blocks(iblk)%val(3, ii, jj, ispin)
              end do
              offset_flat = offset_base + 3 * nrow + 1
            end do
          end do
        end do
        call cpu_time(t1)
        time_reduce_pack = time_reduce_pack + (t1 - t0)
        chunk_size = 262144
        call cpu_time(t0)
        chunk_begin = 1
        do while (chunk_begin <= total_size)
          chunk_count = min(chunk_size, total_size - chunk_begin + 1)
          call comm_summation(send_flat(chunk_begin:chunk_begin + chunk_count - 1), &
                              recv_flat(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, dg_frag%icomm)
          chunk_begin = chunk_begin + chunk_count
        end do
        call cpu_time(t1)
        time_reduce_comm = time_reduce_comm + (t1 - t0)
        call cpu_time(t0)
        offset_flat = 1
        do iblk = 1, dg_frag%n_momentum_blocks
          nrow = dg_frag%momentum_blocks(iblk)%nrow_max
          ncol = dg_frag%momentum_blocks(iblk)%ncol_max
          offset_base = offset_flat - 1
          do ispin = 1, dg_frag%nspin
            do jj = 1, ncol
              !$omp simd
              do ii = 1, nrow
                dg_frag%momentum_blocks(iblk)%val(1, ii, jj, ispin) = recv_flat(offset_base + &
                  ((((ispin - 1) * ncol + (jj - 1)) * nrow + (ii - 1)) * 3 + 1))
                dg_frag%momentum_blocks(iblk)%val(2, ii, jj, ispin) = recv_flat(offset_base + &
                  ((((ispin - 1) * ncol + (jj - 1)) * nrow + (ii - 1)) * 3 + 2))
                dg_frag%momentum_blocks(iblk)%val(3, ii, jj, ispin) = recv_flat(offset_base + &
                  ((((ispin - 1) * ncol + (jj - 1)) * nrow + (ii - 1)) * 3 + 3))
              end do
            end do
          end do
          offset_flat = offset_base + 3 * nrow * ncol * dg_frag%nspin + 1
        end do
        call cpu_time(t1)
        time_reduce_unpack = time_reduce_unpack + (t1 - t0)
        deallocate(send_flat, recv_flat)
      end if
      call cpu_time(t_reduce_end)
      time_block_reduce = time_block_reduce + (t_reduce_end - t_reduce_start)
    end block
    end if
    if (enable_momentum_probe) then
      write(*,'(1x,a,i0,a)') "        momentum-probe: rank=", dg_frag%id, " stage=after-reduce"
      flush(6)
    end if

    ! Enforce anti-symmetry blockwise: self blocks against themselves, off-diagonal
    ! blocks against the reverse ordered fragment pair.
    if (enable_momentum_probe) then
      write(*,'(1x,a,i0,a)') "        momentum-probe: rank=", dg_frag%id, " stage=before-antisym"
      flush(6)
    end if
    call cpu_time(t0)
    do ispin = 1, system%nspin
      do iblk = 1, dg_frag%n_momentum_blocks
        ifrag = dg_frag%momentum_blocks(iblk)%ifrag_row
        jfrag = dg_frag%momentum_blocks(iblk)%ifrag_col
        if (ifrag == jfrag) then
          ndiag = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%momentum_blocks(iblk)%val, 2), &
                      size(dg_frag%momentum_blocks(iblk)%val, 3))
          do idir = 1, 3
            do ii = 1, ndiag
              dg_frag%momentum_blocks(iblk)%val(idir, ii, ii, ispin) = 0.0d0
              do jj = ii + 1, ndiag
                pavg = 0.5d0 * (dg_frag%momentum_blocks(iblk)%val(idir, ii, jj, ispin) - &
                                dg_frag%momentum_blocks(iblk)%val(idir, jj, ii, ispin))
                dg_frag%momentum_blocks(iblk)%val(idir, ii, jj, ispin) = pavg
                dg_frag%momentum_blocks(iblk)%val(idir, jj, ii, ispin) = -pavg
              end do
            end do
          end do
        else
          iblk_rev = find_momentum_block(dg_frag, jfrag, ifrag)
          if (iblk_rev <= 0 .or. iblk >= iblk_rev) cycle
          ni = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%momentum_blocks(iblk)%val, 2), &
                   size(dg_frag%momentum_blocks(iblk_rev)%val, 3))
          nj = min(dg_frag%n_basis(jfrag, ispin), size(dg_frag%momentum_blocks(iblk)%val, 3), &
                   size(dg_frag%momentum_blocks(iblk_rev)%val, 2))
          do idir = 1, 3
            do ii = 1, ni
              do jj = 1, nj
                pavg = 0.5d0 * (dg_frag%momentum_blocks(iblk)%val(idir, ii, jj, ispin) - &
                                dg_frag%momentum_blocks(iblk_rev)%val(idir, jj, ii, ispin))
                dg_frag%momentum_blocks(iblk)%val(idir, ii, jj, ispin) = pavg
                dg_frag%momentum_blocks(iblk_rev)%val(idir, jj, ii, ispin) = -pavg
              end do
            end do
          end do
        end if
      end do
    end do
    call cpu_time(t1)
    time_antisym = time_antisym + (t1 - t0)
    if (enable_momentum_probe) then
      write(*,'(1x,a,i0,a)') "        momentum-probe: rank=", dg_frag%id, " stage=after-antisym"
      flush(6)
    end if
    if (enable_hmat_timing_trace .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,1pe12.4)') "        momentum antisym done time=", time_antisym
      flush(6)
    end if

    max_p = 0.0d0
    do iblk = 1, dg_frag%n_momentum_blocks
      max_p = max(max_p, maxval(abs(dg_frag%momentum_blocks(iblk)%val)))
    end do
    if (enable_momentum_probe) then
      write(*,'(1x,a,i0,a,1pe12.4)') "        momentum-probe: rank=", dg_frag%id, &
        " stage=after-maxp max=", max_p
      flush(6)
    end if
    if (enable_hmat_timing_trace .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
        "        momentum timing: halo_exchange=", time_halo_exchange, &
        " grad=", time_grad_total, " self=", time_self_integral, &
        " halo=", time_halo_integral, &
        " reduce=", time_block_reduce, " antisym=", time_antisym
      write(*,'(a,1pe12.4)') "        Max |momentum_mat|: ", max_p
      write(*,'(a,i0,a,i0,a)') "        Total matrix elements: ", &
                               3 * dg_frag%n_mat_max * dg_frag%n_mat_max * system%nspin, &
                               " (", 3, " directions × basis states × spin)"
    end if
    if (enable_momentum_probe) then
      write(*,'(1x,a,i0,a)') "        momentum-probe: rank=", dg_frag%id, " stage=exit"
      flush(6)
    end if
    
  end subroutine calculate_momentum_matrix

  !=======================================================================
  ! Calculate overlap matrix in DG basis (S_ij = <phi_i|phi_j>)
  !=======================================================================
  subroutine calculate_overlap_matrix(dg_frag, system, mg)
    use structures
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg

    integer :: ifrag, i_local, ispin, io, jo, iblk, iblk_rev, nbf, jo_progress_stride
    integer :: ix, iy, iz, is(3), ie(3), i_halo, jfrag, n_basis_halo
    integer :: ig_row, ig_col, l(3), d(3), ii, jj, halo_send_idx(3), halo_recv_idx(3)
    integer :: lx, ly, lz, gx, gy, gz, iorg(3), ndom(3), loc_s(3), loc_e(3), halo_s(3), halo_e(3)
    integer :: npts_local, nx_local, ny_local, ipt
    integer :: phi_loc_s(3), phi_loc_e(3)
    integer :: lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi
    integer :: phi_lb1, phi_lb2, phi_lb3, phi_ub1, phi_ub2, phi_ub3
    integer :: buf_lb1, buf_lb2, buf_lb3, buf_ub1, buf_ub2, buf_ub3
    logical :: log_frag_progress, release_dense_overlap, has_overlap, use_row_local_overlap
    logical :: build_prop_overlap_blocks
    logical :: enable_overlap_timing_trace, enable_overlap_probe
    logical, save :: overlap_timing_trace_initialized = .false.
    logical, save :: overlap_timing_trace_cached = .false.
    real(8) :: hvol, integral, savg, s_min, s_max, cond_est
    real(8) :: t0, t1, time_self_integral, time_halo_integral, time_reduce_total
    real(8) :: frag_self_start, frag_halo_start
    real(8), allocatable :: phi_local_2d(:,:), self_overlap(:,:)
    character(len=32) :: env_overlap_trace
    integer :: env_overlap_len, env_overlap_stat

    if (.not. dg_frag%has_real_space_basis) return
    enable_overlap_probe = .false.
    call get_environment_variable('SALMON_DG_OVERLAP_PROBE', env_overlap_trace, &
      length=env_overlap_len, status=env_overlap_stat)
    if (env_overlap_stat /= 0 .or. env_overlap_len <= 0) then
      call get_environment_variable('SALMON_DG_TRANSITION_PROBE', env_overlap_trace, &
        length=env_overlap_len, status=env_overlap_stat)
    end if
    if (env_overlap_stat == 0 .and. env_overlap_len > 0) then
      if (env_overlap_trace(1:1) == '1' .or. env_overlap_trace(1:1) == 'y' .or. &
          env_overlap_trace(1:1) == 'Y' .or. env_overlap_trace(1:1) == 't' .or. &
          env_overlap_trace(1:1) == 'T') then
        enable_overlap_probe = .true.
      end if
    end if
    if (enable_overlap_probe) then
      write(*,'(1x,a,i0,a)') "        overlap-probe: rank=", dg_frag%id, " stage=enter"
      flush(6)
    end if
    if (.not. allocated(dg_frag%index_basis) .or. .not. allocated(dg_frag%n_mat)) return

    release_dense_overlap = (.not. dg_frag%yn_adaptive_basis) .and. &
      ((.not. dg_frag%use_plane_wave_basis) .or. dg_frag%n_plane_waves <= 0)

    if (.not. overlap_timing_trace_initialized) then
      env_overlap_trace = ''
      env_overlap_len = 0
      env_overlap_stat = 0
      call get_environment_variable('SALMON_DG_OVERLAP_TIMING_TRACE', env_overlap_trace, &
        length=env_overlap_len, status=env_overlap_stat)
      if (env_overlap_stat == 0 .and. env_overlap_len > 0) then
        if (env_overlap_trace(1:1) == '1' .or. env_overlap_trace(1:1) == 'y' .or. &
            env_overlap_trace(1:1) == 'Y' .or. env_overlap_trace(1:1) == 't' .or. &
            env_overlap_trace(1:1) == 'T') then
          overlap_timing_trace_cached = .true.
        end if
      end if
      overlap_timing_trace_initialized = .true.
    end if
    enable_overlap_timing_trace = overlap_timing_trace_cached

    if (allocated(dg_frag%S_mat)) deallocate(dg_frag%S_mat)
    if (allocated(dg_frag%S_mat_prop)) deallocate(dg_frag%S_mat_prop)
    dg_frag%has_global_overlap_copy = .true.
    dg_frag%overlap_prop_root_authoritative = .false.
    use_row_local_overlap = dg_frag%use_buffered_basis
    build_prop_overlap_blocks = .not. dg_frag%use_buffered_basis
    if (use_row_local_overlap) release_dense_overlap = .true.
    call init_matrix_blocks(dg_frag, dg_frag%S_mat_blocks, dg_frag%S_block_map, dg_frag%n_S_blocks, &
      row_local_only=use_row_local_overlap)
    if (allocated(dg_frag%S_mat_blocks)) then
      do i_halo = 1, size(dg_frag%S_mat_blocks)
        dg_frag%S_mat_blocks(i_halo)%val(:, :, :) = 0.0d0
      end do
    end if
    if (enable_overlap_probe) then
      write(*,'(1x,a,i0,a)') "        overlap-probe: rank=", dg_frag%id, " stage=after-init-blocks"
      flush(6)
    end if

    is = mg%is
    ie = mg%ie
    hvol = system%hvol
    time_self_integral = 0.0d0
    time_halo_integral = 0.0d0
    time_reduce_total = 0.0d0
    phi_lb1 = lbound(dg_frag%phi_frag, 1)
    phi_lb2 = lbound(dg_frag%phi_frag, 2)
    phi_lb3 = lbound(dg_frag%phi_frag, 3)
    phi_ub1 = ubound(dg_frag%phi_frag, 1)
    phi_ub2 = ubound(dg_frag%phi_frag, 2)
    phi_ub3 = ubound(dg_frag%phi_frag, 3)

    if (.not. dg_frag%use_buffered_basis) then
      call exchange_phi_frag_halo(dg_frag)
    end if

    do ispin = 1, system%nspin
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        log_frag_progress = .true.
        frag_self_start = time_self_integral
        frag_halo_start = time_halo_integral
        if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) then
          write(*,*) "[FATAL] overlap invalid i_local: rank=", dg_frag%id, " id_frag=", dg_frag%id_frag, &
            " ifrag=", ifrag, " i_local=", i_local, " phi_dim5=", size(dg_frag%phi_frag, 5)
          stop 1
        end if
        iorg(:) = dg_frag%ixyz_frag(:, ifrag)
        ndom(:) = dg_frag%nxyz_domain(:, ifrag)
        call get_fragment_owned_range(dg_frag, ifrag, mg, loc_s, loc_e, has_overlap)
        if (.not. has_overlap) cycle
        nbf = dg_frag%n_basis(ifrag, ispin)
        jo_progress_stride = max(1, nbf / 4)
        phi_loc_s(:) = iorg(:) + loc_s(:) - 1
        phi_loc_e(:) = iorg(:) + loc_e(:) - 1
        lx_lo = phi_loc_s(1)
        lx_hi = phi_loc_e(1)
        ly_lo = phi_loc_s(2)
        ly_hi = phi_loc_e(2)
        lz_lo = phi_loc_s(3)
        lz_hi = phi_loc_e(3)
        if (phi_loc_s(1) < phi_lb1 .or. phi_loc_e(1) > phi_ub1 .or. &
            phi_loc_s(2) < phi_lb2 .or. phi_loc_e(2) > phi_ub2 .or. &
            phi_loc_s(3) < phi_lb3 .or. phi_loc_e(3) > phi_ub3) then
          write(*,*) "[FATAL] overlap local range out of bounds: ifrag=", ifrag, "ispin=", ispin, &
            "loc_s=", phi_loc_s(1), phi_loc_s(2), phi_loc_s(3), &
            "loc_e=", phi_loc_e(1), phi_loc_e(2), phi_loc_e(3), &
            "phi_lb=", phi_lb1, phi_lb2, phi_lb3, "phi_ub=", phi_ub1, phi_ub2, phi_ub3
          stop 1
        end if
        if (nbf > size(dg_frag%index_basis, 1)) then
          write(*,*) "[FATAL] overlap n_basis exceeds index_basis dim1: rank=", dg_frag%id, &
            " ifrag=", ifrag, " ispin=", ispin, " n_basis=", nbf, &
            " index_basis_dim1=", size(dg_frag%index_basis, 1)
          stop 1
        end if
        if (nbf > size(dg_frag%phi_frag, 4)) then
          write(*,*) "[FATAL] overlap n_basis exceeds phi dim4: rank=", dg_frag%id, &
            " ifrag=", ifrag, " ispin=", ispin, " n_basis=", nbf, &
            " phi_dim4=", size(dg_frag%phi_frag, 4)
          stop 1
        end if

        if (nbf <= 0) cycle
        npts_local = (lx_hi - lx_lo + 1) * (ly_hi - ly_lo + 1) * (lz_hi - lz_lo + 1)
        nx_local = lx_hi - lx_lo + 1
        ny_local = ly_hi - ly_lo + 1
        allocate(phi_local_2d(npts_local, nbf), self_overlap(nbf, nbf))
        !$omp parallel do private(io,lz,ly,lx,ipt) schedule(static)
        do io = 1, nbf
          do lz = lz_lo, lz_hi
            do ly = ly_lo, ly_hi
              !$omp simd private(ipt)
              do lx = lx_lo, lx_hi
                ipt = ((lz - lz_lo) * ny_local + (ly - ly_lo)) * nx_local + (lx - lx_lo) + 1
                phi_local_2d(ipt, io) = dg_frag%phi_frag(lx, ly, lz, io, i_local)
              end do
            end do
          end do
        end do
        !$omp end parallel do
        call cpu_time(t0)
        call dgemm('T', 'N', nbf, nbf, npts_local, hvol, phi_local_2d, npts_local, &
                   phi_local_2d, npts_local, 0.0d0, self_overlap, nbf)
        call cpu_time(t1)
        time_self_integral = time_self_integral + (t1 - t0)
        iblk = find_matrix_block(dg_frag%S_block_map, ifrag, ifrag)

        do jo = 1, nbf
          ig_col = dg_frag%index_basis(jo, ifrag, ispin)
          if (ig_col < 1 .or. ig_col > dg_frag%n_mat_max) cycle
          if (iblk <= 0) cycle
          if (jo > size(dg_frag%S_mat_blocks(iblk)%val, 2)) cycle

          do io = 1, nbf
            ig_row = dg_frag%index_basis(io, ifrag, ispin)
            if (ig_row < 1 .or. ig_row > dg_frag%n_mat_max) cycle
            integral = self_overlap(io, jo)
            if (io <= size(dg_frag%S_mat_blocks(iblk)%val, 1)) then
              dg_frag%S_mat_blocks(iblk)%val(io, jo, ispin) = integral
            end if
          end do
          ! Overlap is a core-domain metric.
          ! Halo values are temporary stencil/projector references only and must
          ! not contribute directly to S_ij. Off-diagonal fragment blocks remain
          ! zero unless they are supported by core-domain ownership on some rank.

        end do
        deallocate(phi_local_2d, self_overlap)
      end do
    end do
    if (enable_overlap_probe) then
      write(*,'(1x,a,i0,a)') "        overlap-probe: rank=", dg_frag%id, " stage=after-local-loops"
      flush(6)
    end if

    if (.not. use_row_local_overlap) then
      call cpu_time(t0)
      call reduce_matrix_blocks(dg_frag, dg_frag%S_mat_blocks, "smat-frag", dg_frag%icomm_frag)
      if (.not. dg_frag%is_frag_root) then
        if (allocated(dg_frag%S_mat_blocks)) then
          do i_halo = 1, size(dg_frag%S_mat_blocks)
            if (fragment_row_is_locally_owned(dg_frag, dg_frag%S_mat_blocks(i_halo)%ifrag_row)) cycle
            dg_frag%S_mat_blocks(i_halo)%val(:, :, :) = 0.0d0
          end do
        end if
      end if
      call cpu_time(t1)
      time_reduce_total = time_reduce_total + (t1 - t0)
    end if
    if (enable_overlap_probe) then
      write(*,'(1x,a,i0,a)') "        overlap-probe: rank=", dg_frag%id, " stage=after-reduce"
      flush(6)
    end if
    if (enable_overlap_timing_trace .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,1pe12.4)') "        overlap reduce done time=", time_reduce_total
      flush(6)
    end if

    if (build_prop_overlap_blocks) then
      call init_matrix_blocks(dg_frag, dg_frag%S_mat_prop_blocks, dg_frag%S_block_map, dg_frag%n_S_blocks, &
        row_local_only=use_row_local_overlap)
      if (allocated(dg_frag%S_mat_blocks) .and. allocated(dg_frag%S_mat_prop_blocks)) then
        do i_halo = 1, size(dg_frag%S_mat_blocks)
          dg_frag%S_mat_prop_blocks(i_halo)%val(:, :, :) = dg_frag%S_mat_blocks(i_halo)%val(:, :, :)
        end do
      end if
    else
      if (allocated(dg_frag%S_mat_prop_blocks)) then
        do i_halo = 1, size(dg_frag%S_mat_prop_blocks)
          if (allocated(dg_frag%S_mat_prop_blocks(i_halo)%val)) deallocate(dg_frag%S_mat_prop_blocks(i_halo)%val)
        end do
        deallocate(dg_frag%S_mat_prop_blocks)
      end if
    end if

    do ispin = 1, dg_frag%nspin
      do ifrag = 1, dg_frag%n_frag
        iblk = find_matrix_block(dg_frag%S_block_map, ifrag, ifrag)
        if (iblk <= 0) cycle
        do ii = 1, min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%S_mat_blocks(iblk)%val, 1), &
                       size(dg_frag%S_mat_blocks(iblk)%val, 2))
          if (dg_frag%S_mat_blocks(iblk)%val(ii, ii, ispin) < 1.0d-12) then
            dg_frag%S_mat_blocks(iblk)%val(ii, ii, ispin) = 1.0d0
            if (allocated(dg_frag%S_mat_prop_blocks)) then
              dg_frag%S_mat_prop_blocks(iblk)%val(ii, ii, ispin) = 1.0d0
            end if
          end if
        end do
      end do
    end do

    if (.not. release_dense_overlap) then
      allocate(dg_frag%S_mat(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
      allocate(dg_frag%S_mat_prop(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
      call sync_blocks_to_dense_matrix(dg_frag, dg_frag%S_mat_blocks, dg_frag%S_block_map, dg_frag%S_mat)
      do ispin = 1, dg_frag%nspin
        do ii = 1, dg_frag%n_mat_max
          if (dg_frag%S_mat(ii, ii, ispin) < 1.0d-12) dg_frag%S_mat(ii, ii, ispin) = 1.0d0
          do jj = ii + 1, dg_frag%n_mat_max
            savg = 0.5d0 * (dg_frag%S_mat(ii, jj, ispin) + dg_frag%S_mat(jj, ii, ispin))
            dg_frag%S_mat(ii, jj, ispin) = savg
            dg_frag%S_mat(jj, ii, ispin) = savg
          end do
        end do
      end do
      dg_frag%S_mat_prop(:, :, :) = dg_frag%S_mat(:, :, :)
    end if
    dg_frag%has_global_overlap_copy = .false.
    dg_frag%overlap_prop_root_authoritative = .false.
    if (enable_overlap_probe) then
      write(*,'(1x,a,i0,a)') "        overlap-probe: rank=", dg_frag%id, " stage=exit"
      flush(6)
    end if

  end subroutine calculate_overlap_matrix
