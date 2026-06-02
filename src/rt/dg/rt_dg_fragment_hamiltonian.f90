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
    integer :: axis, n_same, n_adjacent
    logical :: same_axis, adjacent_axis

    is_pair = .false.
    if (ifrag_row == ifrag_col) then
      is_pair = .true.
      return
    end if
    n_same = 0
    n_adjacent = 0
    do axis = 1, 3
      same_axis = (dg_frag%ixyz_frag(axis, ifrag_row) == dg_frag%ixyz_frag(axis, ifrag_col) .and. &
                   dg_frag%nxyz_domain(axis, ifrag_row) == dg_frag%nxyz_domain(axis, ifrag_col))
      adjacent_axis = is_momentum_neighbor_axis(dg_frag%lgnum_total(axis), &
        dg_frag%ixyz_frag(axis, ifrag_row), dg_frag%nxyz_domain(axis, ifrag_row), &
        dg_frag%ixyz_frag(axis, ifrag_col), dg_frag%nxyz_domain(axis, ifrag_col)) .and. .not. same_axis
      if (same_axis) n_same = n_same + 1
      if (adjacent_axis) n_adjacent = n_adjacent + 1
    end do

    is_pair = (n_same == 2 .and. n_adjacent == 1)
  end function is_momentum_neighbor_pair

  subroutine ensure_momentum_neighbor_pair_cache(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: ifrag_row, ifrag_col

    if (allocated(dg_frag%momentum_neighbor_pair_cache)) return
    allocate(dg_frag%momentum_neighbor_pair_cache(dg_frag%n_frag, dg_frag%n_frag))
    dg_frag%momentum_neighbor_pair_cache(:, :) = .false.
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (ifrag_row == ifrag_col) then
          dg_frag%momentum_neighbor_pair_cache(ifrag_row, ifrag_col) = .true.
        else
          dg_frag%momentum_neighbor_pair_cache(ifrag_row, ifrag_col) = &
            is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)
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

  integer function map_global_to_periodic_box_coord_ham(ig, lb, ub) result(iloc)
    implicit none
    integer, intent(in) :: ig, lb, ub
    integer :: nbox

    nbox = ub - lb + 1
    if (nbox <= 0) then
      iloc = 0
      return
    end if
    iloc = lb + modulo(ig - lb, nbox)
  end function map_global_to_periodic_box_coord_ham

  subroutine copy_periodic_global_scalar_to_rank_buffer(dg_frag, grid, field_global, field_buffer)
    use structures, only: s_rgrid, s_scalar
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_rgrid), intent(in) :: grid
    type(s_scalar), intent(in) :: field_global
    real(8), intent(out) :: field_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                                         dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                                         dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3))
    integer :: ix, iy, iz
    integer :: gx, gy, gz

!$omp parallel do private(ix,iy,gx,gy,gz) schedule(static)
    do iz = dg_frag%rank_buf_lo(3), dg_frag%rank_buf_hi(3)
      gz = map_global_to_phi_box_coord_ham(iz, grid%is(3), grid%ie(3), dg_frag%lgnum_total(3))
      do iy = dg_frag%rank_buf_lo(2), dg_frag%rank_buf_hi(2)
        gy = map_global_to_phi_box_coord_ham(iy, grid%is(2), grid%ie(2), dg_frag%lgnum_total(2))
!$omp simd
        do ix = dg_frag%rank_buf_lo(1), dg_frag%rank_buf_hi(1)
          gx = map_global_to_phi_box_coord_ham(ix, grid%is(1), grid%ie(1), dg_frag%lgnum_total(1))
          if (gx == 0 .or. gy == 0 .or. gz == 0) then
            field_buffer(ix, iy, iz) = 0.0d0
          else
            field_buffer(ix, iy, iz) = field_global%f(gx, gy, gz)
          end if
        end do
      end do
    end do
!$omp end parallel do
  end subroutine copy_periodic_global_scalar_to_rank_buffer

  integer function map_global_to_fragment_phi_box_coord_ham(dg_frag, ifrag, axis, ig, lb, ub) result(iloc)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, axis, ig, lb, ub
    integer :: core_lo, ndom, ig_wrap

    core_lo = dg_frag%ixyz_frag(axis, ifrag)
    ndom = dg_frag%nxyz_domain(axis, ifrag)
    if (ndom <= 0) then
      iloc = 0
      return
    end if
    ig_wrap = core_lo + modulo(ig - core_lo, ndom)
    ig_wrap = modulo(ig_wrap - 1, dg_frag%lgnum_total(axis)) + 1
    iloc = map_global_to_phi_box_coord_ham(ig_wrap, lb, ub, dg_frag%lgnum_total(axis))
  end function map_global_to_fragment_phi_box_coord_ham

  subroutine enforce_fragment_periodic_buffer_for_state_ham(dg_frag, ifrag, i_local, jo)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ifrag, i_local, jo

    integer :: ix, iy, iz
    integer :: sx, sy, sz
    integer :: lb1, ub1, lb2, ub2, lb3, ub3

    lb1 = lbound(dg_frag%phi_frag, 1)
    ub1 = ubound(dg_frag%phi_frag, 1)
    lb2 = lbound(dg_frag%phi_frag, 2)
    ub2 = ubound(dg_frag%phi_frag, 2)
    lb3 = lbound(dg_frag%phi_frag, 3)
    ub3 = ubound(dg_frag%phi_frag, 3)

    do iz = lb3, ub3
      sz = map_global_to_fragment_phi_box_coord_ham(dg_frag, ifrag, 3, iz, lb3, ub3)
      if (sz == 0) cycle
      do iy = lb2, ub2
        sy = map_global_to_fragment_phi_box_coord_ham(dg_frag, ifrag, 2, iy, lb2, ub2)
        if (sy == 0) cycle
        do ix = lb1, ub1
          sx = map_global_to_fragment_phi_box_coord_ham(dg_frag, ifrag, 1, ix, lb1, ub1)
          if (sx == 0) cycle
          dg_frag%phi_frag(ix, iy, iz, jo, i_local) = dg_frag%phi_frag(sx, sy, sz, jo, i_local)
        end do
      end do
    end do
  end subroutine enforce_fragment_periodic_buffer_for_state_ham

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

  logical function basis_row_is_locally_owned(dg_frag, ifrag, ispin, io) result(is_local)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, ispin, io
    integer :: global_idx

    is_local = .false.
    if (ifrag < 1 .or. ifrag > dg_frag%n_frag) return
    if (ispin < 1 .or. ispin > dg_frag%nspin) return
    if (.not. allocated(dg_frag%coef_owner)) return
    if (.not. allocated(dg_frag%index_basis)) return
    if (io < 1 .or. io > size(dg_frag%index_basis, 1)) return

    global_idx = dg_frag%index_basis(io, ifrag, ispin)
    if (global_idx < 1 .or. global_idx > size(dg_frag%coef_owner, 1)) return
    is_local = (dg_frag%coef_owner(global_idx, ispin) == dg_frag%id)
  end function basis_row_is_locally_owned

  subroutine init_matrix_blocks(dg_frag, blocks, block_map, n_blocks, diagonal_only)
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(matrix_block_info), allocatable, intent(inout) :: blocks(:)
    integer, allocatable, intent(inout) :: block_map(:, :)
    integer, intent(out) :: n_blocks
    logical, intent(in), optional :: diagonal_only
    integer :: ifrag_row, ifrag_col, iblk, nrow_max, ncol_max
    logical :: diagonal_blocks_only, local_fragment_only, include_pair

    if (allocated(blocks)) then
      do iblk = 1, size(blocks)
        if (allocated(blocks(iblk)%val)) deallocate(blocks(iblk)%val)
      end do
      deallocate(blocks)
    end if
    if (allocated(block_map)) deallocate(block_map)
    call ensure_momentum_neighbor_pair_cache(dg_frag)
    ! Galerkin volume operators are fragment-local.  S and local-potential H
    ! must stay block diagonal; inter-fragment coupling belongs in an explicit
    ! boundary/flow term, not in these dense volume blocks.
    diagonal_blocks_only = .true.
    if (present(diagonal_only)) diagonal_blocks_only = diagonal_only
    local_fragment_only = dg_frag%parallel_mode_orbital

    n_blocks = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (local_fragment_only) then
          if (ifrag_row /= dg_frag%ifrag_group) cycle
          if (diagonal_blocks_only .and. ifrag_col /= dg_frag%ifrag_group) cycle
        end if
        if (diagonal_blocks_only) then
          include_pair = (ifrag_row == ifrag_col)
        else
          include_pair = is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)
        end if
        if (include_pair) n_blocks = n_blocks + 1
      end do
    end do
    if (n_blocks <= 0) return

    allocate(blocks(n_blocks))
    allocate(block_map(dg_frag%n_frag, dg_frag%n_frag))
    block_map = 0

    iblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (local_fragment_only) then
          if (ifrag_row /= dg_frag%ifrag_group) cycle
          if (diagonal_blocks_only .and. ifrag_col /= dg_frag%ifrag_group) cycle
        end if
        if (diagonal_blocks_only) then
          include_pair = (ifrag_row == ifrag_col)
        else
          include_pair = is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)
        end if
        if (.not. include_pair) cycle
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

  subroutine init_momentum_blocks(dg_frag, diagonal_only)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    logical, intent(in), optional :: diagonal_only
    integer :: ifrag_row, ifrag_col, nblk, iblk
    integer :: nrow_max, ncol_max
    logical :: diagonal_blocks_only, local_fragment_only, include_pair

    if (allocated(dg_frag%momentum_blocks)) then
      do iblk = 1, size(dg_frag%momentum_blocks)
        if (allocated(dg_frag%momentum_blocks(iblk)%val)) deallocate(dg_frag%momentum_blocks(iblk)%val)
      end do
      deallocate(dg_frag%momentum_blocks)
    end if
    if (allocated(dg_frag%momentum_blocks_c)) then
      do iblk = 1, size(dg_frag%momentum_blocks_c)
        if (allocated(dg_frag%momentum_blocks_c(iblk)%val)) deallocate(dg_frag%momentum_blocks_c(iblk)%val)
      end do
      deallocate(dg_frag%momentum_blocks_c)
    end if
    if (allocated(dg_frag%momentum_block_map)) deallocate(dg_frag%momentum_block_map)
    call ensure_momentum_neighbor_pair_cache(dg_frag)
    ! The volume momentum operator is also fragment-local.  Boundary current
    ! corrections must be added as an explicit flow term so that current across
    ! fragment faces is not confused with dense neighbor volume coupling.
    diagonal_blocks_only = .true.
    if (present(diagonal_only)) diagonal_blocks_only = diagonal_only
    local_fragment_only = dg_frag%parallel_mode_orbital

    nblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (local_fragment_only) then
          if (ifrag_row /= dg_frag%ifrag_group) cycle
          if (diagonal_blocks_only .and. ifrag_col /= dg_frag%ifrag_group) cycle
        end if
        if (diagonal_blocks_only) then
          include_pair = (ifrag_row == ifrag_col)
        else
          include_pair = is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)
        end if
        if (include_pair) nblk = nblk + 1
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
        if (local_fragment_only) then
          if (ifrag_row /= dg_frag%ifrag_group) cycle
          if (diagonal_blocks_only .and. ifrag_col /= dg_frag%ifrag_group) cycle
        end if
        if (diagonal_blocks_only) then
          include_pair = (ifrag_row == ifrag_col)
        else
          include_pair = is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)
        end if
        if (.not. include_pair) cycle
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

  integer function face_neighbor_fragment_ham(dg_frag, ifrag, axis, side) result(jfrag)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, axis, side
    integer :: coords(3), rem

    jfrag = 0
    if (ifrag < 1 .or. ifrag > dg_frag%n_frag) return
    if (axis < 1 .or. axis > 3) return
    if (dg_frag%num_fragment(axis) <= 1) return

    coords(1) = (ifrag - 1) / (dg_frag%num_fragment(2) * dg_frag%num_fragment(3)) + 1
    rem = modulo(ifrag - 1, dg_frag%num_fragment(2) * dg_frag%num_fragment(3))
    coords(2) = rem / dg_frag%num_fragment(3) + 1
    coords(3) = modulo(rem, dg_frag%num_fragment(3)) + 1
    coords(axis) = modulo(coords(axis) - 1 + side + dg_frag%num_fragment(axis), &
                          dg_frag%num_fragment(axis)) + 1
    jfrag = ((coords(1) - 1) * dg_frag%num_fragment(2) + coords(2) - 1) * &
            dg_frag%num_fragment(3) + coords(3)
  end function face_neighbor_fragment_ham

  logical function is_periodic_boundary_face_ham(dg_frag, ifrag, axis, side) result(is_pbc)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, axis, side
    integer :: frag_lo, frag_hi

    is_pbc = .false.
    if (ifrag < 1 .or. ifrag > dg_frag%n_frag) return
    if (axis < 1 .or. axis > 3) return
    if (dg_frag%num_fragment(axis) <= 1) return
    frag_lo = dg_frag%ixyz_frag(axis, ifrag)
    frag_hi = frag_lo + dg_frag%nxyz_domain(axis, ifrag) - 1
    if (side > 0) then
      is_pbc = (frag_hi == dg_frag%lgnum_total(axis))
    else if (side < 0) then
      is_pbc = (frag_lo == 1)
    end if
  end function is_periodic_boundary_face_ham

  subroutine read_fragment_buffer_basis_box_ham(dg_frag, ifrag, phi_box, box_lo, box_hi)
    use filesystem, only: get_filehandle
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag
    real(8), allocatable, intent(out) :: phi_box(:,:,:,:)
    integer, intent(out) :: box_lo(3), box_hi(3)

    character(32), parameter :: bdir_frag = './data_dcdft/fragments/'
    character(32), parameter :: binfile_bfb = 'basis_functions_buffer.bin'
    integer, parameter :: basis_buffer_magic = -22022212
    integer, parameter :: basis_buffer_version = 1
    character(256) :: filename
    integer :: iunit, magic_file, version_file
    integer :: nxyz_domain(3), nxyz_buffer_file(3), nxyz_box(3)
    integer :: nspin_file, nstate_frag_file, n_basis_keep
    integer :: ispin_file, istate, ix, ixg, iyg, izg
    integer :: ix_box, iy_box, iz_box, ix_src, iy_src, iz_src
    integer, allocatable :: n_basis_frag(:)
    real(8), allocatable :: phi_read(:,:,:)
    logical :: has_file

    if (allocated(phi_box)) deallocate(phi_box)
    box_lo(:) = 0
    box_hi(:) = -1
    write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_bfb
    inquire(file=filename, exist=has_file)
    if (.not. has_file) then
      write(*,'(1x,a,i0,a,a)') '[FATAL] missing neighbor DG buffer basis at ifrag=', ifrag, &
        ' file=', trim(filename)
      stop 'DG-Fragment RT: missing neighbor basis buffer file'
    end if

    iunit = get_filehandle()
    open(iunit, file=filename, form='unformatted', access='stream', status='old')
    read(iunit) magic_file, version_file
    if (magic_file /= basis_buffer_magic .or. version_file /= basis_buffer_version) then
      write(*,'(1x,a,i0,a,i0,a,i0)') '[FATAL] invalid neighbor basis buffer header at ifrag=', &
        ifrag, ' magic=', magic_file, ' version=', version_file
      stop 'DG-Fragment RT: invalid neighbor basis buffer file'
    end if
    read(iunit) nxyz_domain(1:3), nxyz_buffer_file(1:3), nspin_file, nstate_frag_file
    do ix = 1, 3
      if (dg_frag%num_fragment(ix) > 1 .and. nxyz_buffer_file(ix) < dg_frag%nxyz_buffer(ix)) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') '[FATAL] neighbor DG buffer too small at ifrag=', &
          ifrag, ' axis=', ix, ' seed_buffer=', nxyz_buffer_file(ix), &
          ' required=', dg_frag%nxyz_buffer(ix)
        stop 'DG-Fragment RT: insufficient neighbor basis buffer'
      end if
    end do

    allocate(n_basis_frag(max(1, nspin_file)))
    read(iunit) n_basis_frag(1:max(1, nspin_file))
    n_basis_keep = 0
    if (nspin_file >= 1) then
      n_basis_keep = min(maxval(dg_frag%n_basis(ifrag, 1:dg_frag%nspin)), nstate_frag_file, dg_frag%nstate_frag)
    end if
    nxyz_box(:) = nxyz_domain(:) + 2 * nxyz_buffer_file(:)
    box_lo(:) = dg_frag%ixyz_frag(:, ifrag) - nxyz_buffer_file(:)
    box_hi(:) = dg_frag%ixyz_frag(:, ifrag) + nxyz_domain(:) - 1 + nxyz_buffer_file(:)
    allocate(phi_box(1:nxyz_box(1), 1:nxyz_box(2), 1:nxyz_box(3), max(1, n_basis_keep)))
    allocate(phi_read(1:nxyz_box(1), 1:nxyz_box(2), 1:nxyz_box(3)))
    phi_box = 0.0d0

    do ispin_file = 1, nspin_file
      do istate = 1, nstate_frag_file
        read(iunit) phi_read(1:nxyz_box(1), 1:nxyz_box(2), 1:nxyz_box(3))
        if (ispin_file /= 1 .or. istate > n_basis_keep) cycle
        do izg = box_lo(3), box_hi(3)
          iz_box = izg - (dg_frag%ixyz_frag(3, ifrag) - nxyz_buffer_file(3)) + 1
          if (iz_box < 1 .or. iz_box > nxyz_box(3)) then
            iz_src = nxyz_buffer_file(3) + modulo(izg - dg_frag%ixyz_frag(3, ifrag), nxyz_domain(3)) + 1
          else
            iz_src = iz_box
          end if
          do iyg = box_lo(2), box_hi(2)
            iy_box = iyg - (dg_frag%ixyz_frag(2, ifrag) - nxyz_buffer_file(2)) + 1
            if (iy_box < 1 .or. iy_box > nxyz_box(2)) then
              iy_src = nxyz_buffer_file(2) + modulo(iyg - dg_frag%ixyz_frag(2, ifrag), nxyz_domain(2)) + 1
            else
              iy_src = iy_box
            end if
            do ixg = box_lo(1), box_hi(1)
              ix_box = ixg - (dg_frag%ixyz_frag(1, ifrag) - nxyz_buffer_file(1)) + 1
              if (ix_box < 1 .or. ix_box > nxyz_box(1)) then
                ix_src = nxyz_buffer_file(1) + modulo(ixg - dg_frag%ixyz_frag(1, ifrag), nxyz_domain(1)) + 1
              else
                ix_src = ix_box
              end if
              phi_box(ixg - box_lo(1) + 1, iyg - box_lo(2) + 1, izg - box_lo(3) + 1, istate) = &
                phi_read(ix_src, iy_src, iz_src)
            end do
          end do
        end do
      end do
    end do
    close(iunit)
    deallocate(phi_read, n_basis_frag)
  end subroutine read_fragment_buffer_basis_box_ham

  real(8) function phi_box_value_ham(phi_box, box_lo, box_hi, lgtot, istate, gidx) result(val)
    implicit none
    real(8), intent(in) :: phi_box(:,:,:,:)
    integer, intent(in) :: box_lo(3), box_hi(3), lgtot(3), istate, gidx(3)
    integer :: ix, iy, iz

    val = 0.0d0
    if (istate < 1 .or. istate > size(phi_box, 4)) return
    ix = map_global_to_phi_box_coord_ham(gidx(1), box_lo(1), box_hi(1), lgtot(1))
    iy = map_global_to_phi_box_coord_ham(gidx(2), box_lo(2), box_hi(2), lgtot(2))
    iz = map_global_to_phi_box_coord_ham(gidx(3), box_lo(3), box_hi(3), lgtot(3))
    if (ix == 0 .or. iy == 0 .or. iz == 0) return
    val = phi_box(ix - box_lo(1) + 1, iy - box_lo(2) + 1, iz - box_lo(3) + 1, istate)
  end function phi_box_value_ham

  real(8) function phi_local_value_ham(dg_frag, i_local, istate, gidx) result(val)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: i_local, istate, gidx(3)
    integer :: ix, iy, iz

    val = 0.0d0
    if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) return
    if (istate < 1 .or. istate > size(dg_frag%phi_frag, 4)) return
    ix = map_global_to_phi_box_coord_ham(gidx(1), lbound(dg_frag%phi_frag, 1), &
                                         ubound(dg_frag%phi_frag, 1), dg_frag%lgnum_total(1))
    iy = map_global_to_phi_box_coord_ham(gidx(2), lbound(dg_frag%phi_frag, 2), &
                                         ubound(dg_frag%phi_frag, 2), dg_frag%lgnum_total(2))
    iz = map_global_to_phi_box_coord_ham(gidx(3), lbound(dg_frag%phi_frag, 3), &
                                         ubound(dg_frag%phi_frag, 3), dg_frag%lgnum_total(3))
    if (ix == 0 .or. iy == 0 .or. iz == 0) return
    val = dg_frag%phi_frag(ix, iy, iz, istate, i_local)
  end function phi_local_value_ham

  real(8) function phi_box_dn_ham(phi_box, box_lo, box_hi, lgtot, stencil, istate, gidx, axis, side) result(dn)
    use structures, only: s_stencil
    implicit none
    real(8), intent(in) :: phi_box(:,:,:,:)
    integer, intent(in) :: box_lo(3), box_hi(3), lgtot(3), istate, gidx(3), axis, side
    type(s_stencil), intent(in) :: stencil
    integer :: ix, iy, iz, ixc, iyc, izc
    real(8) :: nabt(4,3), grad_axis

    dn = 0.0d0
    if (istate < 1 .or. istate > size(phi_box, 4)) return
    ix = map_global_to_phi_box_coord_ham(gidx(1), box_lo(1), box_hi(1), lgtot(1))
    iy = map_global_to_phi_box_coord_ham(gidx(2), box_lo(2), box_hi(2), lgtot(2))
    iz = map_global_to_phi_box_coord_ham(gidx(3), box_lo(3), box_hi(3), lgtot(3))
    if (ix == 0 .or. iy == 0 .or. iz == 0) return
    ixc = ix - box_lo(1) + 1
    iyc = iy - box_lo(2) + 1
    izc = iz - box_lo(3) + 1
    nabt = stencil%coef_nab
    select case(axis)
    case(1)
      if (ixc - 4 < 1 .or. ixc + 4 > size(phi_box, 1)) return
      grad_axis = nabt(1,1) * (phi_box(ixc + 1, iyc, izc, istate) - phi_box(ixc - 1, iyc, izc, istate)) + &
                  nabt(2,1) * (phi_box(ixc + 2, iyc, izc, istate) - phi_box(ixc - 2, iyc, izc, istate)) + &
                  nabt(3,1) * (phi_box(ixc + 3, iyc, izc, istate) - phi_box(ixc - 3, iyc, izc, istate)) + &
                  nabt(4,1) * (phi_box(ixc + 4, iyc, izc, istate) - phi_box(ixc - 4, iyc, izc, istate))
    case(2)
      if (iyc - 4 < 1 .or. iyc + 4 > size(phi_box, 2)) return
      grad_axis = nabt(1,2) * (phi_box(ixc, iyc + 1, izc, istate) - phi_box(ixc, iyc - 1, izc, istate)) + &
                  nabt(2,2) * (phi_box(ixc, iyc + 2, izc, istate) - phi_box(ixc, iyc - 2, izc, istate)) + &
                  nabt(3,2) * (phi_box(ixc, iyc + 3, izc, istate) - phi_box(ixc, iyc - 3, izc, istate)) + &
                  nabt(4,2) * (phi_box(ixc, iyc + 4, izc, istate) - phi_box(ixc, iyc - 4, izc, istate))
    case(3)
      if (izc - 4 < 1 .or. izc + 4 > size(phi_box, 3)) return
      grad_axis = nabt(1,3) * (phi_box(ixc, iyc, izc + 1, istate) - phi_box(ixc, iyc, izc - 1, istate)) + &
                  nabt(2,3) * (phi_box(ixc, iyc, izc + 2, istate) - phi_box(ixc, iyc, izc - 2, istate)) + &
                  nabt(3,3) * (phi_box(ixc, iyc, izc + 3, istate) - phi_box(ixc, iyc, izc - 3, istate)) + &
                  nabt(4,3) * (phi_box(ixc, iyc, izc + 4, istate) - phi_box(ixc, iyc, izc - 4, istate))
    case default
      return
    end select
    dn = real(side, 8) * grad_axis
  end function phi_box_dn_ham

  real(8) function phi_local_dn_ham(dg_frag, i_local, mg, stencil, istate, gidx, axis, side) result(dn)
    use structures, only: s_rgrid, s_stencil
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: i_local, istate, gidx(3), axis, side
    type(s_rgrid), intent(in) :: mg
    type(s_stencil), intent(in) :: stencil
    real(8) :: gx, gy, gz

    dn = 0.0d0
    call apply_gradient_at_phi_box_point(dg_frag, i_local, istate, mg, stencil, gidx, gx, gy, gz)
    select case(axis)
    case(1)
      dn = real(side, 8) * gx
    case(2)
      dn = real(side, 8) * gy
    case(3)
      dn = real(side, 8) * gz
    end select
  end function phi_local_dn_ham

  subroutine add_dg_surface_flux_blocks(dg_frag, system, mg, stencil, H_blocks, T_blocks, block_map)
    use structures, only: s_dft_system, s_rgrid, s_stencil
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_rgrid), intent(in) :: mg
    type(s_stencil), intent(in) :: stencil
    type(matrix_block_info), intent(inout) :: H_blocks(:)
    type(matrix_block_info), intent(inout) :: T_blocks(:)
    integer, intent(in) :: block_map(:, :)

    real(8), parameter :: surface_penalty_factor = 10.0d0
    integer :: ifrag, jfrag, i_local, axis, side, ispin, io, jo
    integer :: iblk_self, iblk_cross, nrow, ncol, ix, iy, iz
    integer :: g_row(3), g_col(3), loop_lo(3), loop_hi(3)
    integer :: col_lo(3), col_hi(3)
    real(8) :: area_weight, alpha, u_l, u_r, v_l, dnu_l, dnu_r, dnv_l
    real(8) :: term_self, term_cross
    real(8) :: face_self_sum, face_cross_sum, face_self_max, face_cross_max
    integer :: face_self_count, face_cross_count
    logical :: pbc_face
    real(8), allocatable :: phi_col(:,:,:,:)

    if (.not. allocated(dg_frag%phi_frag)) return
    if (.not. allocated(dg_frag%n_basis)) return
    if (.not. allocated(dg_frag%index_basis)) return
    if (.not. dg_frag%parallel_mode_orbital) return

    ! Equation (4) of the DG-ALB formulation: keep the volume kinetic and
    ! local-potential terms fragment-local, then add only the surface
    ! jump/average/penalty contribution as sparse face-neighbor blocks.
    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) cycle
      if (.not. dg_frag%parallel_mode_orbital .and. .not. dg_frag%is_frag_root) cycle

      do axis = 1, 3
        if (dg_frag%num_fragment(axis) <= 1) cycle
        if (system%hgs(axis) <= 0.0d0) cycle
        area_weight = system%hvol / system%hgs(axis)
        alpha = surface_penalty_factor / system%hgs(axis)
        do side = -1, 1, 2
          jfrag = face_neighbor_fragment_ham(dg_frag, ifrag, axis, side)
          if (jfrag <= 0 .or. jfrag == ifrag) cycle
          iblk_self = find_matrix_block(block_map, ifrag, ifrag)
          iblk_cross = find_matrix_block(block_map, ifrag, jfrag)
          if (iblk_self <= 0 .and. iblk_cross <= 0) cycle

          call read_fragment_buffer_basis_box_ham(dg_frag, jfrag, phi_col, col_lo, col_hi)

          loop_lo(:) = dg_frag%ixyz_frag(:, ifrag)
          loop_hi(:) = dg_frag%ixyz_frag(:, ifrag) + dg_frag%nxyz_domain(:, ifrag) - 1
          if (side > 0) then
            loop_lo(axis) = dg_frag%ixyz_frag(axis, ifrag) + dg_frag%nxyz_domain(axis, ifrag) - 1
            loop_hi(axis) = loop_lo(axis)
          else
            loop_lo(axis) = dg_frag%ixyz_frag(axis, ifrag)
            loop_hi(axis) = loop_lo(axis)
          end if
          pbc_face = is_periodic_boundary_face_ham(dg_frag, ifrag, axis, side)
          face_self_sum = 0.0d0
          face_cross_sum = 0.0d0
          face_self_max = 0.0d0
          face_cross_max = 0.0d0
          face_self_count = 0
          face_cross_count = 0

          do ispin = 1, dg_frag%nspin
            nrow = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), size(dg_frag%phi_frag, 4))
            ncol = min(dg_frag%n_basis(jfrag, ispin), size(dg_frag%index_basis, 1), size(phi_col, 4))
            if (nrow <= 0) cycle

            do iz = loop_lo(3), loop_hi(3)
              do iy = loop_lo(2), loop_hi(2)
                do ix = loop_lo(1), loop_hi(1)
                  g_row = [ix, iy, iz]
                  g_col = g_row
                  if (side > 0) then
                    g_col(axis) = dg_frag%ixyz_frag(axis, jfrag)
                  else
                    g_col(axis) = dg_frag%ixyz_frag(axis, jfrag) + dg_frag%nxyz_domain(axis, jfrag) - 1
                  end if

                  if (iblk_self > 0) then
                    do jo = 1, nrow
                      u_l = phi_local_value_ham(dg_frag, i_local, jo, g_row)
                      dnu_l = phi_local_dn_ham(dg_frag, i_local, mg, stencil, jo, g_row, axis, side)
                      do io = 1, nrow
                        if (dg_frag%parallel_mode_orbital) then
                          if (.not. basis_row_is_locally_owned(dg_frag, ifrag, ispin, io)) cycle
                        end if
                        v_l = phi_local_value_ham(dg_frag, i_local, io, g_row)
                        dnv_l = phi_local_dn_ham(dg_frag, i_local, mg, stencil, io, g_row, axis, side)
                        term_self = (-0.25d0 * v_l * dnu_l - 0.25d0 * dnv_l * u_l + alpha * v_l * u_l) * area_weight
                        H_blocks(iblk_self)%val(io, jo, ispin) = H_blocks(iblk_self)%val(io, jo, ispin) + term_self
                        T_blocks(iblk_self)%val(io, jo, ispin) = T_blocks(iblk_self)%val(io, jo, ispin) + term_self
                        if (pbc_face) then
                          face_self_sum = face_self_sum + term_self * term_self
                          face_self_max = max(face_self_max, abs(term_self))
                          face_self_count = face_self_count + 1
                        end if
                      end do
                    end do
                  end if

                  if (iblk_cross > 0 .and. ncol > 0) then
                    do jo = 1, ncol
                      u_r = phi_box_value_ham(phi_col, col_lo, col_hi, dg_frag%lgnum_total, jo, g_col)
                      dnu_r = phi_box_dn_ham(phi_col, col_lo, col_hi, dg_frag%lgnum_total, stencil, &
                                             jo, g_col, axis, side)
                      do io = 1, nrow
                        if (dg_frag%parallel_mode_orbital) then
                          if (.not. basis_row_is_locally_owned(dg_frag, ifrag, ispin, io)) cycle
                        end if
                        v_l = phi_local_value_ham(dg_frag, i_local, io, g_row)
                        dnv_l = phi_local_dn_ham(dg_frag, i_local, mg, stencil, io, g_row, axis, side)
                        term_cross = (-0.25d0 * v_l * dnu_r + 0.25d0 * dnv_l * u_r - alpha * v_l * u_r) * area_weight
                        H_blocks(iblk_cross)%val(io, jo, ispin) = H_blocks(iblk_cross)%val(io, jo, ispin) + term_cross
                        T_blocks(iblk_cross)%val(io, jo, ispin) = T_blocks(iblk_cross)%val(io, jo, ispin) + term_cross
                        if (pbc_face) then
                          face_cross_sum = face_cross_sum + term_cross * term_cross
                          face_cross_max = max(face_cross_max, abs(term_cross))
                          face_cross_count = face_cross_count + 1
                        end if
                      end do
                    end do
                  end if
                end do
              end do
            end do
          end do

          if (pbc_face) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x))') &
              '[DG-PBC-FACE] kind=H rank=', dg_frag%id, ' ifrag=', ifrag, &
              ' jfrag=', jfrag, ' axis=', axis, ' side=', side, &
              ' iblk_self=', iblk_self, ' iblk_cross=', iblk_cross, &
              ' row_lo=', dg_frag%ixyz_frag(:, ifrag), ' row_hi=', &
              dg_frag%ixyz_frag(:, ifrag) + dg_frag%nxyz_domain(:, ifrag) - 1
            write(*,'(1x,a,2(i0,1x),a,2(i0,1x),a,2(1x,1pe14.6),a,2(1x,1pe14.6),a,2(i0,1x))') &
              '[DG-PBC-FACE-NORM] kind=H col_box_x=', col_lo(1), col_hi(1), &
              ' col_box_y=', col_lo(2), col_hi(2), &
              ' self_frob cross_frob=', sqrt(face_self_sum), sqrt(face_cross_sum), &
              ' self_max cross_max=', face_self_max, face_cross_max, &
              ' counts=', face_self_count, face_cross_count
            flush(6)
          end if

          if (allocated(phi_col)) deallocate(phi_col)
        end do
      end do
    end do
  end subroutine add_dg_surface_flux_blocks

  subroutine add_dg_surface_momentum_blocks(dg_frag, system)
    use structures, only: s_dft_system
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system

    integer :: ifrag, jfrag, i_local, axis, side, ispin, io, jo
    integer :: iblk_self, iblk_cross, nrow, ncol, ix, iy, iz
    integer :: g_row(3), g_col(3), loop_lo(3), loop_hi(3)
    integer :: col_lo(3), col_hi(3)
    real(8) :: area_weight, normal_sign
    real(8) :: u_l, u_r, v_l, term_self, term_cross
    real(8) :: face_self_sum, face_cross_sum, face_self_max, face_cross_max
    integer :: face_self_count, face_cross_count
    logical :: pbc_face
    real(8), allocatable :: phi_col(:,:,:,:)

    if (.not. allocated(dg_frag%momentum_blocks)) return
    if (.not. allocated(dg_frag%momentum_block_map)) return
    if (.not. allocated(dg_frag%phi_frag)) return
    if (.not. allocated(dg_frag%n_basis)) return
    if (.not. allocated(dg_frag%index_basis)) return
    if (.not. dg_frag%parallel_mode_orbital) return

    ! Momentum/current uses the DG face correction for the normal gradient.
    ! Local volume gradients stay inside each fragment; only u_R-u_L on a
    ! shared face is written into off-diagonal flow blocks.
    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) cycle

      do axis = 1, 3
        if (dg_frag%num_fragment(axis) <= 1) cycle
        if (system%hgs(axis) <= 0.0d0) cycle
        area_weight = system%hvol / system%hgs(axis)
        do side = -1, 1, 2
          normal_sign = real(side, kind=8)
          jfrag = face_neighbor_fragment_ham(dg_frag, ifrag, axis, side)
          if (jfrag <= 0 .or. jfrag == ifrag) cycle
          iblk_self = find_momentum_block(dg_frag, ifrag, ifrag)
          iblk_cross = find_momentum_block(dg_frag, ifrag, jfrag)
          if (iblk_self <= 0 .and. iblk_cross <= 0) cycle

          call read_fragment_buffer_basis_box_ham(dg_frag, jfrag, phi_col, col_lo, col_hi)

          loop_lo(:) = dg_frag%ixyz_frag(:, ifrag)
          loop_hi(:) = dg_frag%ixyz_frag(:, ifrag) + dg_frag%nxyz_domain(:, ifrag) - 1
          if (side > 0) then
            loop_lo(axis) = dg_frag%ixyz_frag(axis, ifrag) + dg_frag%nxyz_domain(axis, ifrag) - 1
            loop_hi(axis) = loop_lo(axis)
          else
            loop_lo(axis) = dg_frag%ixyz_frag(axis, ifrag)
            loop_hi(axis) = loop_lo(axis)
          end if
          pbc_face = is_periodic_boundary_face_ham(dg_frag, ifrag, axis, side)
          face_self_sum = 0.0d0
          face_cross_sum = 0.0d0
          face_self_max = 0.0d0
          face_cross_max = 0.0d0
          face_self_count = 0
          face_cross_count = 0

          do ispin = 1, dg_frag%nspin
            nrow = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), size(dg_frag%phi_frag, 4))
            ncol = min(dg_frag%n_basis(jfrag, ispin), size(dg_frag%index_basis, 1), size(phi_col, 4))
            if (nrow <= 0) cycle

            do iz = loop_lo(3), loop_hi(3)
              do iy = loop_lo(2), loop_hi(2)
                do ix = loop_lo(1), loop_hi(1)
                  g_row = [ix, iy, iz]
                  g_col = g_row
                  if (side > 0) then
                    g_col(axis) = dg_frag%ixyz_frag(axis, jfrag)
                  else
                    g_col(axis) = dg_frag%ixyz_frag(axis, jfrag) + dg_frag%nxyz_domain(axis, jfrag) - 1
                  end if

                  if (iblk_self > 0) then
                    do jo = 1, nrow
                      u_l = phi_local_value_ham(dg_frag, i_local, jo, g_row)
                      do io = 1, nrow
                        if (.not. basis_row_is_locally_owned(dg_frag, ifrag, ispin, io)) cycle
                        v_l = phi_local_value_ham(dg_frag, i_local, io, g_row)
                        term_self = -0.5d0 * normal_sign * v_l * u_l * area_weight
                        dg_frag%momentum_blocks(iblk_self)%val(axis, io, jo, ispin) = &
                          dg_frag%momentum_blocks(iblk_self)%val(axis, io, jo, ispin) + term_self
                        if (pbc_face) then
                          face_self_sum = face_self_sum + term_self * term_self
                          face_self_max = max(face_self_max, abs(term_self))
                          face_self_count = face_self_count + 1
                        end if
                      end do
                    end do
                  end if

                  if (iblk_cross > 0 .and. ncol > 0) then
                    do jo = 1, ncol
                      u_r = phi_box_value_ham(phi_col, col_lo, col_hi, dg_frag%lgnum_total, jo, g_col)
                      do io = 1, nrow
                        if (.not. basis_row_is_locally_owned(dg_frag, ifrag, ispin, io)) cycle
                        v_l = phi_local_value_ham(dg_frag, i_local, io, g_row)
                        term_cross = 0.5d0 * normal_sign * v_l * u_r * area_weight
                        dg_frag%momentum_blocks(iblk_cross)%val(axis, io, jo, ispin) = &
                          dg_frag%momentum_blocks(iblk_cross)%val(axis, io, jo, ispin) + term_cross
                        if (pbc_face) then
                          face_cross_sum = face_cross_sum + term_cross * term_cross
                          face_cross_max = max(face_cross_max, abs(term_cross))
                          face_cross_count = face_cross_count + 1
                        end if
                      end do
                    end do
                  end if
                end do
              end do
            end do
          end do

          if (pbc_face) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x))') &
              '[DG-PBC-FACE] kind=M rank=', dg_frag%id, ' ifrag=', ifrag, &
              ' jfrag=', jfrag, ' axis=', axis, ' side=', side, &
              ' iblk_self=', iblk_self, ' iblk_cross=', iblk_cross, &
              ' row_lo=', dg_frag%ixyz_frag(:, ifrag), ' row_hi=', &
              dg_frag%ixyz_frag(:, ifrag) + dg_frag%nxyz_domain(:, ifrag) - 1
            write(*,'(1x,a,2(i0,1x),a,2(i0,1x),a,2(1x,1pe14.6),a,2(1x,1pe14.6),a,2(i0,1x))') &
              '[DG-PBC-FACE-NORM] kind=M col_box_x=', col_lo(1), col_hi(1), &
              ' col_box_y=', col_lo(2), col_hi(2), &
              ' self_frob cross_frob=', sqrt(face_self_sum), sqrt(face_cross_sum), &
              ' self_max cross_max=', face_self_max, face_cross_max, &
              ' counts=', face_self_count, face_cross_count
            flush(6)
          end if

          if (allocated(phi_col)) deallocate(phi_col)
        end do
      end do
    end do
  end subroutine add_dg_surface_momentum_blocks

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
    logical, save :: reduce_trace_initialized = .false.
    logical, save :: enable_reduce_trace = .false.
    character(16) :: env_reduce_trace
    integer :: env_status

    if (.not. reduce_trace_initialized) then
      env_reduce_trace = ''
      call get_environment_variable('SALMON_DG_HMAT_REDUCE_TRACE', env_reduce_trace, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(env_reduce_trace)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enable_reduce_trace = .true.
        end select
      end if
      reduce_trace_initialized = .true.
    end if

    max_block_size = 0
    total_active_size = 0
    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = nrow * ncol
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

    if (max_block_size_global <= 0) return
    allocate(send_block(max_block_size_global), recv_block(max_block_size_global))

    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = nrow * ncol
        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            send_block(offset_flat) = blocks(iblk)%val(ii, jj, ispin)
            offset_flat = offset_flat + 1
          end do
        end do

        chunk_begin = 1
        do while (chunk_begin <= block_size)
          chunk_count = min(reduce_chunk_size, block_size - chunk_begin + 1)
          call comm_summation(send_block(chunk_begin:chunk_begin + chunk_count - 1), &
                              recv_block(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, icomm_reduce)
          chunk_begin = chunk_begin + chunk_count
        end do

        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            blocks(iblk)%val(ii, jj, ispin) = recv_block(offset_flat)
            offset_flat = offset_flat + 1
          end do
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

  subroutine trace_matrix_blocks_if_enabled(dg_frag, blocks, label)
    use communication, only: comm_is_root
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(in) :: blocks(:)
    character(*), intent(in) :: label

    logical, save :: initialized = .false.
    logical, save :: enabled = .false.
    character(16) :: env_trace
    integer :: env_status
    integer :: iblk, ispin, nrow, ncol
    real(8) :: frob, maxabs

    if (.not. initialized) then
      env_trace = ''
      call get_environment_variable('SALMON_DG_HMAT_BLOCK_TRACE', env_trace, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(env_trace)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enabled = .true.
        end select
      end if
      initialized = .true.
    end if
    if (.not. enabled) return
    if (.not. comm_is_root(dg_frag%id)) return

    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = min(dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin), size(blocks(iblk)%val, 1))
        ncol = min(dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin), size(blocks(iblk)%val, 2))
        if (nrow <= 0 .or. ncol <= 0) cycle
        frob = sqrt(sum(blocks(iblk)%val(1:nrow, 1:ncol, ispin)**2))
        maxabs = maxval(abs(blocks(iblk)%val(1:nrow, 1:ncol, ispin)))
        write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0,a,1pe14.6,a,1pe14.6)') &
          '[HMAT-BLOCK] label=', trim(label), ' iblk=', iblk, &
          ' row=', blocks(iblk)%ifrag_row, ' col=', blocks(iblk)%ifrag_col, &
          ' ispin=', ispin, ' frob=', frob, ' max=', maxabs
      end do
    end do
    flush(6)
  end subroutine trace_matrix_blocks_if_enabled

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
    integer :: loop_ifrag_start, loop_ifrag_end
    integer :: ndom(3)
    integer :: i_halo
    real(8) :: hvol
    real(8) :: max_p
    real(8) :: Ac_zero(3)
    logical :: is_local_fragment
    integer :: is(3), ie(3)
    real(8), allocatable :: T_phi(:,:,:)  ! Kinetic energy operator applied to basis (fragment-local)
    real(8), allocatable :: H_phi(:,:,:)  ! Hamiltonian-applied field H|phi_j> = T|phi_j> + V|phi_j> (fragment-local)
    real(8), allocatable :: V_total(:,:,:)  ! Total potential V = Vpsl + Vh + Vxc
    real(8), allocatable :: partial_t(:), partial_h(:), reduced_t(:), reduced_h(:)
    real(8), allocatable :: partial_th(:), reduced_th(:)
    type(matrix_block_info), allocatable :: H_diag_blocks(:), H_kin_diag_blocks(:)
    integer :: n_local_diag, nbf_max, i_diag, iblk, iblk_rev, nbf_diag, nbf_comm
    integer :: jo_s, jo_e, jo_loc, ncol_local
    integer :: n_metric
    logical :: release_dense_fragment_ops
    complex(8), allocatable :: H_metric_ref(:,:)
    
    release_dense_fragment_ops = (.not. dg_frag%yn_adaptive_basis) .and. &
      ((.not. dg_frag%use_plane_wave_basis) .or. dg_frag%n_plane_waves <= 0)

    if (.not. dg_frag%has_real_space_basis) then
      return
    end if

    ! Enforce fragment-local stencil policy: no halo communication path.
    dg_frag%n_halo = 0
    dg_frag%has_halo_exchange = .false.
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*)
      write(*,*) "=== Preparing Hamiltonian Matrix ==="
    end if

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
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "  [2/3] Constructing Hamiltonian matrix H = T + V..."
    end if
    
    if (allocated(dg_frag%H_mat)) deallocate(dg_frag%H_mat)
    if (allocated(dg_frag%H_mat_kinetic)) deallocate(dg_frag%H_mat_kinetic)
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
    
    ! Halo exchange removed: stencil operations use local phi_frag with fragment PBC buffer only.
    
    hvol = system%hvol
    ! Construct total potential: V = Vpsl + Vh + Vxc
    ! Note: This is used for initial H_mat calculation
    do ispin = 1, system%nspin
      loop_ifrag_start = dg_frag%ifrag_start
      loop_ifrag_end = dg_frag%ifrag_end

      do ifrag = loop_ifrag_start, loop_ifrag_end
        ndom(:) = dg_frag%nxyz_domain(:, ifrag)
        if (any(ndom <= 0)) then
          write(*,*) "[FATAL] hamiltonian step2 non-positive ndom: rank=", dg_frag%id, &
            " ifrag=", ifrag, " ndom=", ndom
          stop 1
        end if
        allocate(V_total(1:ndom(1), 1:ndom(2), 1:ndom(3)))
        call build_fragment_total_potential_grid(dg_frag, ifrag, mg, Vh, Vxc(ispin), Vpsl, V_total)

        is_local_fragment = .true.
        i_local = ifrag - dg_frag%ifrag_start + 1
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
        if (dg_frag%parallel_mode_orbital) then
          ! Matrix blocks are kept local to the row-owning orbital rank.  Build
          ! all ket columns locally so propagation does not need a global block
          ! replication/reduction step.
          jo_s = 1
          jo_e = nbf
        else
          call get_fragment_column_range(dg_frag, nbf, jo_s, jo_e)
        end if
        if (jo_s > jo_e) cycle
        if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) then
          write(*,*) "[FATAL] hamiltonian step2 invalid i_local: rank=", dg_frag%id, &
            " ifrag=", ifrag, " i_local=", i_local, " phi_frag_dim5=", size(dg_frag%phi_frag, 5)
          stop 1
        end if
        ! T_phi/H_phi/V_total are indexed by fragment-local box coordinates
        ! (1:ndom).  Orbital mode gathers V_total in global fragment order so
        ! column-split ranks all see the same full fragment potential.
        allocate(T_phi(1:ndom(1), 1:ndom(2), 1:ndom(3)))
        allocate(H_phi(1:ndom(1), 1:ndom(2), 1:ndom(3)))
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
        ! Orbital mode stores row-owned coefficient slices; build only those
        ! bra-side rows while still using all ket columns of the local fragment.
        ! Legacy mode keeps the old full-column spatial-box reduction.
        do jo = jo_s, jo_e
          ig_j = dg_frag%index_basis(jo, ifrag, ispin)
          if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle
          call build_hpsi_for_basis(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi)

          partial_t(:) = 0.0d0
          partial_h(:) = 0.0d0
          !$omp parallel do private(ig_i)
          do io = 1, nbf
            ig_i = dg_frag%index_basis(io, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            if (dg_frag%parallel_mode_orbital) then
              if (.not. basis_row_is_locally_owned(dg_frag, ifrag, ispin, io)) cycle
            end if

            ! Kinetic energy matrix element: T_ij = ∫ φ_i (T|φ_j>) dr
            call integrate_basis_with_field(dg_frag, ifrag, i_local, io, mg, T_phi, hvol, partial_t(io))

            ! Store kinetic part
            call integrate_basis_with_field(dg_frag, ifrag, i_local, io, mg, H_phi, hvol, partial_h(io))

          end do
          !$omp end parallel do
          partial_th(1:nbf_comm) = partial_t(1:nbf_comm)
          partial_th(nbf_comm + 1:2 * nbf_comm) = partial_h(1:nbf_comm)
          if (dg_frag%parallel_mode_orbital) then
            ! Orbital mode integrates the full fragment box locally and keeps
            ! only the row-local self block, so no subgroup sum is needed here.
            reduced_t(1:nbf_comm) = partial_t(1:nbf_comm)
            reduced_h(1:nbf_comm) = partial_h(1:nbf_comm)
            do io = 1, nbf
              ig_i = dg_frag%index_basis(io, ifrag, ispin)
              if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
              if (.not. basis_row_is_locally_owned(dg_frag, ifrag, ispin, io)) cycle
              H_kin_diag_blocks(i_local)%val(io, jo, ispin) = reduced_t(io)
              H_diag_blocks(i_local)%val(io, jo, ispin) = reduced_h(io)
            end do
          else
            ! Legacy real-space mode reduces parent-grid box contributions.
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
          end if

        end do  ! jo loop
        deallocate(partial_t, partial_h, reduced_t, reduced_h)
        deallocate(partial_th, reduced_th)
        deallocate(V_total)
        deallocate(T_phi, H_phi)
        if (allocated(T_phi) .or. allocated(H_phi)) then
          write(*,*) "[FATAL] hamiltonian deallocate(T_phi,H_phi) failed: rank=", dg_frag%id, &
            " ifrag=", ifrag, " ispin=", ispin
          stop 1
        end if
          
        
      end do  ! ifrag loop
      
    end do  ! ispin loop
    
    call init_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, dg_frag%H_block_map, dg_frag%n_H_blocks, &
                            diagonal_only=.false.)
    call init_matrix_blocks(dg_frag, dg_frag%H_mat_core_blocks, dg_frag%H_block_map, dg_frag%n_H_blocks, &
                            diagonal_only=.false.)
    call init_matrix_blocks(dg_frag, dg_frag%H_mat_kinetic_blocks, dg_frag%H_block_map, dg_frag%n_H_blocks, &
                            diagonal_only=.false.)
    call rebuild_local_h_block_ids(dg_frag)
    do i_diag = 1, n_local_diag
      iblk = find_matrix_block(dg_frag%H_block_map, H_diag_blocks(i_diag)%ifrag_row, H_diag_blocks(i_diag)%ifrag_col)
      if (iblk <= 0) cycle
      do ispin = 1, dg_frag%nspin
        nbf_diag = dg_frag%n_basis(H_diag_blocks(i_diag)%ifrag_row, ispin)
        if (nbf_diag <= 0) cycle
        dg_frag%H_mat_blocks(iblk)%val(1:nbf_diag, 1:nbf_diag, ispin) = H_diag_blocks(i_diag)%val(1:nbf_diag, 1:nbf_diag, ispin)
        dg_frag%H_mat_core_blocks(iblk)%val(1:nbf_diag, 1:nbf_diag, ispin) = &
          H_diag_blocks(i_diag)%val(1:nbf_diag, 1:nbf_diag, ispin)
        dg_frag%H_mat_kinetic_blocks(iblk)%val(1:nbf_diag, 1:nbf_diag, ispin) = &
          H_kin_diag_blocks(i_diag)%val(1:nbf_diag, 1:nbf_diag, ispin)
      end do
    end do

    call add_dg_surface_flux_blocks(dg_frag, system, mg, stencil, &
                                    dg_frag%H_mat_blocks, dg_frag%H_mat_kinetic_blocks, &
                                    dg_frag%H_block_map)

    ! Non-buffered halo projection route has been removed.

    do i_diag = 1, n_local_diag
      if (allocated(H_diag_blocks(i_diag)%val)) deallocate(H_diag_blocks(i_diag)%val)
      if (allocated(H_kin_diag_blocks(i_diag)%val)) deallocate(H_kin_diag_blocks(i_diag)%val)
    end do
    if (allocated(H_diag_blocks)) deallocate(H_diag_blocks)
    if (allocated(H_kin_diag_blocks)) deallocate(H_kin_diag_blocks)
    ! CRITICAL: MPI aggregation of Hamiltonian matrix
    ! Each rank computed elements only for its assigned fragments.
    ! Reduce one fragment block at a time to avoid a single dense global allreduce.
    if (.not. dg_frag%parallel_mode_orbital) then
      call reduce_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, "hmat", dg_frag%icomm)
      call reduce_matrix_blocks(dg_frag, dg_frag%H_mat_core_blocks, "hmat-core", dg_frag%icomm)
    end if
    if (.not. dg_frag%is_frag_root) then
      if (dg_frag%parallel_mode_orbital) then
        ! Orbital row-split keeps the reduced matrix replicated on all
        ! subgroup ranks; coefficient ownership selects the rows to update.
        continue
      else
        ! Legacy real-space mode: zero out blocks not owned by this rank.
        if (allocated(dg_frag%H_mat_blocks)) then
          do i_halo = 1, size(dg_frag%H_mat_blocks)
            if (fragment_row_is_locally_owned(dg_frag, dg_frag%H_mat_blocks(i_halo)%ifrag_row)) cycle
            dg_frag%H_mat_blocks(i_halo)%val(:, :, :) = 0.0d0
          end do
        end if
        if (allocated(dg_frag%H_mat_core_blocks)) then
          do i_halo = 1, size(dg_frag%H_mat_core_blocks)
            if (fragment_row_is_locally_owned(dg_frag, dg_frag%H_mat_core_blocks(i_halo)%ifrag_row)) cycle
            dg_frag%H_mat_core_blocks(i_halo)%val(:, :, :) = 0.0d0
          end do
        end if
      end if
    end if
    if (.not. dg_frag%parallel_mode_orbital) then
      call reduce_matrix_blocks(dg_frag, dg_frag%H_mat_kinetic_blocks, "hmat-kinetic", dg_frag%icomm)
    end if
    if (.not. dg_frag%is_frag_root) then
      if (dg_frag%parallel_mode_orbital) then
        ! Keep the replicated kinetic block for row-owned coefficient updates.
        continue
      else
        if (allocated(dg_frag%H_mat_kinetic_blocks)) then
          do i_halo = 1, size(dg_frag%H_mat_kinetic_blocks)
            if (fragment_row_is_locally_owned(dg_frag, dg_frag%H_mat_kinetic_blocks(i_halo)%ifrag_row)) cycle
            dg_frag%H_mat_kinetic_blocks(i_halo)%val(:, :, :) = 0.0d0
          end do
        end if
      end if
    end if
    call symmetrize_real_matrix_blocks(dg_frag, dg_frag%H_mat_blocks)
    call symmetrize_real_matrix_blocks(dg_frag, dg_frag%H_mat_core_blocks)
    call symmetrize_real_matrix_blocks(dg_frag, dg_frag%H_mat_kinetic_blocks)
    call trace_matrix_blocks_if_enabled(dg_frag, dg_frag%H_mat_blocks, "H")
    call trace_matrix_blocks_if_enabled(dg_frag, dg_frag%H_mat_kinetic_blocks, "T")
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

    ! Build initial mixed basis once (fragment + PW) with A=0.
    ! This is the first rollout of the mixed-basis orthogonalization path:
    ! initialization uses the mixed diagonalization, while later RT updates
    ! still fall back to the existing mixed refresh path.
    if (dg_frag%use_plane_wave_basis .and. dg_frag%n_plane_waves > 0) then
      Ac_zero(:) = 0.0d0
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "  [init] Diagonalizing mixed basis at startup (A=0)"
      end if
      call diagonalize_mixed_basis(dg_frag, system, Vh, Vxc, Vpsl, Ac_zero)
      if (allocated(dg_frag%coef_new)) dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
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
    
    if (allocated(V_total)) deallocate(V_total)
    if (allocated(V_total)) then
      write(*,*) "[FATAL] V_total still allocated before return: rank=", dg_frag%id
      stop 1
    end if
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "=== Hamiltonian Matrix Ready ==="
      write(*,*)
    end if
    
  end subroutine calculate_hamiltonian_matrix

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

    if (comm_is_root(dg_frag%id)) then
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
    if (comm_is_root(dg_frag%id)) then
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
!$omp parallel do collapse(2) private(ix)
    do iz = grid_z_lo, grid_z_hi
      do iy = grid_y_lo, grid_y_hi
!$omp simd
        do ix = grid_x_lo, grid_x_hi
          V_total(ix, iy, iz) = Vpsl%f(ix, iy, iz) + Vh%f(ix, iy, iz) + Vxc_spin%f(ix, iy, iz)
        end do
      end do
    end do
!$omp end parallel do
  end subroutine build_total_potential_grid

  subroutine build_fragment_total_potential_grid(dg_frag, ifrag, mg, Vh, Vxc_spin, Vpsl, V_total)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag
    type(s_rgrid), intent(in) :: mg
    type(s_scalar), intent(in) :: Vh, Vxc_spin, Vpsl
    real(8), intent(inout) :: V_total(:,:,:)

    integer :: lx, ly, lz, gx, gy, gz, gwx, gwy, gwz
    integer :: iorg(3), ndom(3)
    real(8) :: vh_val
    logical :: full_box_local

    V_total(:, :, :) = 0.0d0
    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    if (dg_frag%parallel_mode_orbital) then
      full_box_local = .true.
      do lz = 1, ndom(3)
        gz = iorg(3) + lz - 1
        gwz = map_global_to_periodic_box_coord_ham(gz, 1, dg_frag%lgnum_total(3))
        do ly = 1, ndom(2)
          gy = iorg(2) + ly - 1
          gwy = map_global_to_periodic_box_coord_ham(gy, 1, dg_frag%lgnum_total(2))
          do lx = 1, ndom(1)
            gx = iorg(1) + lx - 1
            gwx = map_global_to_periodic_box_coord_ham(gx, 1, dg_frag%lgnum_total(1))
            if (gwx < mg%is(1) .or. gwx > mg%ie(1) .or. &
                gwy < mg%is(2) .or. gwy > mg%ie(2) .or. &
                gwz < mg%is(3) .or. gwz > mg%ie(3)) then
              full_box_local = .false.
              exit
            end if
          end do
          if (.not. full_box_local) exit
        end do
        if (.not. full_box_local) exit
      end do
      if (.not. full_box_local) then
        write(*,'(1x,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
          "[FATAL] DG orbital H build requires the parent rank to own the full fragment box: rank=", dg_frag%id, &
          " ifrag=", ifrag, " mg_is=", mg%is(1), mg%is(2), mg%is(3), &
          " mg_ie=", mg%ie(1), mg%ie(2), mg%ie(3), " ndom=", ndom(1), ndom(2), ndom(3)
        write(*,'(1x,a)') "        Set process_allocation='orbital_sequential' and nproc_rgrid(1:3)=num_fragment(1:3);"
        write(*,'(1x,a)') "        use nproc_ob for orbital parallelism."
        flush(6)
        stop "DG-Fragment RT orbital H build requires local full fragment box"
      end if
    end if
!$omp parallel do private(ly,lx,gz,gy,gx,gwz,gwy,gwx,vh_val) schedule(static)
    do lz = 1, ndom(3)
      gz = iorg(3) + lz - 1
      gwz = map_global_to_periodic_box_coord_ham(gz, 1, dg_frag%lgnum_total(3))
      if (gwz < mg%is(3) .or. gwz > mg%ie(3)) cycle
      do ly = 1, ndom(2)
        gy = iorg(2) + ly - 1
        gwy = map_global_to_periodic_box_coord_ham(gy, 1, dg_frag%lgnum_total(2))
        if (gwy < mg%is(2) .or. gwy > mg%ie(2)) cycle
!$omp simd private(gx,gwx,vh_val)
        do lx = 1, ndom(1)
          gx = iorg(1) + lx - 1
          gwx = map_global_to_periodic_box_coord_ham(gx, 1, dg_frag%lgnum_total(1))
          if (gwx >= mg%is(1) .and. gwx <= mg%ie(1)) then
            vh_val = Vh%f(gwx, gwy, gwz)
            V_total(lx, ly, lz) = Vpsl%f(gwx, gwy, gwz) + vh_val + Vxc_spin%f(gwx, gwy, gwz)
          end if
        end do
      end do
    end do
!$omp end parallel do
  end subroutine build_fragment_total_potential_grid

  subroutine build_total_potential_grid_with_buffered_hartree(grid, dg_frag, Vh_buffer, Vxc_spin, Vpsl, V_total)
    use structures
    implicit none
    type(s_rgrid), intent(in) :: grid
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    real(8), intent(in) :: Vh_buffer(:,:,:)
    type(s_scalar), intent(in) :: Vxc_spin, Vpsl
    real(8), intent(out) :: V_total(grid%is(1):grid%ie(1), grid%is(2):grid%ie(2), grid%is(3):grid%ie(3))
    integer :: ix, iy, iz
    integer :: bx, by, bz
    integer :: b_lo1, b_lo2, b_lo3
    integer :: b_hi1, b_hi2, b_hi3

    b_lo1 = lbound(Vh_buffer, 1)
    b_hi1 = ubound(Vh_buffer, 1)
    b_lo2 = lbound(Vh_buffer, 2)
    b_hi2 = ubound(Vh_buffer, 2)
    b_lo3 = lbound(Vh_buffer, 3)
    b_hi3 = ubound(Vh_buffer, 3)

!$omp parallel do private(ix,iy,bx,by,bz) schedule(static)
    do iz = grid%is(3), grid%ie(3)
      bz = map_global_to_periodic_box_coord_ham(iz, b_lo3, b_hi3)
      do iy = grid%is(2), grid%ie(2)
        by = map_global_to_periodic_box_coord_ham(iy, b_lo2, b_hi2)
!$omp simd
        do ix = grid%is(1), grid%ie(1)
          bx = map_global_to_periodic_box_coord_ham(ix, b_lo1, b_hi1)
          V_total(ix, iy, iz) = Vpsl%f(ix, iy, iz) + Vh_buffer(bx, by, bz) + Vxc_spin%f(ix, iy, iz)
        end do
      end do
    end do
!$omp end parallel do
  end subroutine build_total_potential_grid_with_buffered_hartree

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
    real(8), intent(in) :: V_total(:,:,:)
    real(8), intent(out) :: T_phi(:,:,:)
    real(8), intent(out) :: H_phi(:,:,:)
    integer :: gx, gy, gz, lx, ly, lz
    integer :: bx, by, bz
    integer :: gx0, gx1, gy0, gy1, gz0, gz1
    integer :: iorg(3), ndom(3)
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

    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    call get_fragment_owned_range(dg_frag, ifrag, mg, loc_s, loc_e, has_overlap)
    if (.not. has_overlap) return
    if (loc_s(1) < lbound(T_phi, 1) .or. loc_e(1) > ubound(T_phi, 1) .or. &
        loc_s(2) < lbound(T_phi, 2) .or. loc_e(2) > ubound(T_phi, 2) .or. &
        loc_s(3) < lbound(T_phi, 3) .or. loc_e(3) > ubound(T_phi, 3)) then
      write(*,*) "[FATAL] build_hpsi local range exceeds T_phi bounds: rank=", dg_frag%id, &
        " ifrag=", ifrag, " loc_s=", loc_s, " loc_e=", loc_e, " T_shape=", shape(T_phi)
      stop 1
    end if
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
    if (loc_s(1) < v_lb1 .or. loc_e(1) > v_ub1 .or. &
        loc_s(2) < v_lb2 .or. loc_e(2) > v_ub2 .or. &
        loc_s(3) < v_lb3 .or. loc_e(3) > v_ub3) then
      write(*,*) "[FATAL] build_hpsi V_total index out of bounds: rank=", dg_frag%id, &
        " ifrag=", ifrag, " loc_s=", loc_s, " loc_e=", loc_e, &
        " V_lb=", v_lb1, v_lb2, v_lb3, " V_ub=", v_ub1, v_ub2, v_ub3
      stop 1
    end if
!$omp parallel do private(ly, lx, gz, gy, gx, bx, by, bz) schedule(static)
    do lz = loc_s(3), loc_e(3)
      gz = iorg(3) + lz - 1
      do ly = loc_s(2), loc_e(2)
        gy = iorg(2) + ly - 1
!$omp simd private(bx, by, bz)
        do lx = loc_s(1), loc_e(1)
          gx = iorg(1) + lx - 1
          bx = map_global_to_phi_box_coord_ham(gx, phi_lb1, phi_ub1, dg_frag%lgnum_total(1))
          by = map_global_to_phi_box_coord_ham(gy, phi_lb2, phi_ub2, dg_frag%lgnum_total(2))
          bz = map_global_to_phi_box_coord_ham(gz, phi_lb3, phi_ub3, dg_frag%lgnum_total(3))
          if (bx == 0 .or. by == 0 .or. bz == 0) cycle
          H_phi(lx, ly, lz) = H_phi(lx, ly, lz) + V_total(lx, ly, lz) * dg_frag%phi_frag(bx, by, bz, jo, i_local)
        end do
      end do
    end do
!$omp end parallel do
  end subroutine build_hpsi_for_basis

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
    real(8), intent(in) :: field(:,:,:)
    real(8), intent(in) :: hvol
    real(8), intent(out) :: integral
    real(8) :: partial
    integer :: gx, gy, gz, lx, ly, lz
    integer :: bx, by, bz
    integer :: ndom(3), loc_s(3), loc_e(3)
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
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    call get_fragment_owned_range(dg_frag, ifrag, mg, loc_s, loc_e, has_overlap)
    if (.not. has_overlap) then
      integral = 0.0d0
      return
    end if
    if (loc_s(1) < lbound(field, 1) .or. loc_e(1) > ubound(field, 1) .or. &
        loc_s(2) < lbound(field, 2) .or. loc_e(2) > ubound(field, 2) .or. &
        loc_s(3) < lbound(field, 3) .or. loc_e(3) > ubound(field, 3)) then
      write(*,*) "[FATAL] integrate local range exceeds field bounds: rank=", dg_frag%id, &
        " ifrag=", ifrag, " loc_s=", loc_s, " loc_e=", loc_e, " field_shape=", shape(field)
      stop 1
    end if
    f_lb1 = lbound(dg_frag%phi_frag, 1)
    f_ub1 = ubound(dg_frag%phi_frag, 1)
    f_lb2 = lbound(dg_frag%phi_frag, 2)
    f_ub2 = ubound(dg_frag%phi_frag, 2)
    f_lb3 = lbound(dg_frag%phi_frag, 3)
    f_ub3 = ubound(dg_frag%phi_frag, 3)
    gx0 = dg_frag%ixyz_frag(1, ifrag) + loc_s(1) - 1
    gx1 = dg_frag%ixyz_frag(1, ifrag) + loc_e(1) - 1
    gy0 = dg_frag%ixyz_frag(2, ifrag) + loc_s(2) - 1
    gy1 = dg_frag%ixyz_frag(2, ifrag) + loc_e(2) - 1
    gz0 = dg_frag%ixyz_frag(3, ifrag) + loc_s(3) - 1
    gz1 = dg_frag%ixyz_frag(3, ifrag) + loc_e(3) - 1
    partial = 0.0d0
    do lz = loc_s(3), loc_e(3)
      gz = dg_frag%ixyz_frag(3, ifrag) + lz - 1
      do ly = loc_s(2), loc_e(2)
        gy = dg_frag%ixyz_frag(2, ifrag) + ly - 1
        !$omp simd reduction(+:partial) private(bx, by, bz)
        do lx = loc_s(1), loc_e(1)
          gx = dg_frag%ixyz_frag(1, ifrag) + lx - 1
          bx = map_global_to_phi_box_coord_ham(gx, f_lb1, f_ub1, dg_frag%lgnum_total(1))
          by = map_global_to_phi_box_coord_ham(gy, f_lb2, f_ub2, dg_frag%lgnum_total(2))
          bz = map_global_to_phi_box_coord_ham(gz, f_lb3, f_ub3, dg_frag%lgnum_total(3))
          if (bx == 0 .or. by == 0 .or. bz == 0) cycle
          partial = partial + dg_frag%phi_frag(bx, by, bz, io, i_local) * field(lx, ly, lz) * hvol
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
    integer :: p_lb1, p_ub1, p_lb2, p_ub2, p_lb3, p_ub3
    integer :: lgx, lgy, lgz
    real(8) :: v, lap0
    real(8) :: lapt(4,3)
    integer :: ndom(3), loc_s(3), loc_e(3)
    logical :: has_overlap
    
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
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    if (any(ndom <= 0)) then
      write(*,*) "[FATAL] kinetic non-positive domain size: rank=", dg_frag%id, &
        " ifrag=", ifrag, " ndom=", ndom
      stop 1
    end if
    if (.not. allocated(dg_frag%phi_frag_has_seed_buffer) .or. &
        .not. dg_frag%phi_frag_has_seed_buffer(i_local)) then
      call enforce_fragment_periodic_buffer_for_state_ham(dg_frag, ifrag, i_local, jo)
    end if
    call get_fragment_owned_range(dg_frag, ifrag, mg, loc_s, loc_e, has_overlap)
    lgx = dg_frag%lgnum_total(1)
    lgy = dg_frag%lgnum_total(2)
    lgz = dg_frag%lgnum_total(3)
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub3 = ubound(dg_frag%phi_frag, 3)
    T_phi = 0.0d0
    if (.not. has_overlap) return
    ! Note: phi_frag is allocated as (1-nb:nx+nb, 1-nb:ny+nb, 1-nb:nz+nb, ...)
    ! where nb = nxyz_buffer = 4 for 4th-order stencil
    ! The interior domain is (1:nx, 1:ny, 1:nz), and halo provides data for stencil
    ! operations near boundaries.
    !
    ! With halo exchange, stencil operations can access phi_frag(ix±4, iy±4, iz±4)
    ! for all interior points without boundary restrictions.
    
    ! Apply kinetic operator using finite difference stencil.
    ! Convert each global coordinate to its periodic interior index once,
    ! then use ix0+/-n offsets directly for neighbor access inside phi_frag.

!$omp parallel do private(ly, lx, gz, gy, gx, ix0, iy0, iz0, v) schedule(static)
    do lz = loc_s(3), loc_e(3)
      gz = dg_frag%ixyz_frag(3, ifrag) + lz - 1
      iz0 = map_global_to_phi_box_coord_ham(gz, p_lb3, p_ub3, lgz)
      do ly = loc_s(2), loc_e(2)
        gy = dg_frag%ixyz_frag(2, ifrag) + ly - 1
        iy0 = map_global_to_phi_box_coord_ham(gy, p_lb2, p_ub2, lgy)
        do lx = loc_s(1), loc_e(1)
          gx = dg_frag%ixyz_frag(1, ifrag) + lx - 1
          ix0 = map_global_to_phi_box_coord_ham(gx, p_lb1, p_ub1, lgx)
          if (ix0 == 0 .or. iy0 == 0 .or. iz0 == 0) then
            write(*,'(1x,a,i0,a,i0)') "[FATAL] kinetic buffer coordinate out of bounds: rank=", &
              dg_frag%id, " ifrag=", ifrag
            stop "DG-Fragment RT: kinetic buffer coordinate out of bounds"
          end if
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
          
          ! T|φ> = (-∇²/2)|φ> = lap0*φ - 0.5 * (sum of neighbor terms)
          T_phi(lx, ly, lz) = lap0 * dg_frag%phi_frag(ix0, iy0, iz0, jo, i_local) - 0.5d0 * v
          
        end do
      end do
    end do
!$omp end parallel do
    
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
    if (dg_frag%parallel_mode_orbital) then
      ! Orbital mode replicates the fragment real-space box on every subgroup
      ! rank and distributes matrix work by basis columns, not by grid boxes.
      loc_s(:) = 1
      loc_e(:) = ndom(:)
      has_overlap = all(ndom(:) > 0)
      return
    end if
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

  subroutine get_fragment_mg_overlap_range(dg_frag, ifrag, mg, loc_s, loc_e, has_overlap)
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
  end subroutine get_fragment_mg_overlap_range

  subroutine get_fragment_column_range(dg_frag, ncol, col_s, col_e)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ncol
    integer, intent(out) :: col_s, col_e

    integer :: base, extra, rank_in_frag, nworker

    if (ncol <= 0) then
      col_s = 1
      col_e = 0
      return
    end if
    if (.not. dg_frag%parallel_mode_orbital .or. dg_frag%isize_frag <= 1) then
      col_s = 1
      col_e = ncol
      return
    end if

    nworker = max(1, dg_frag%isize_frag)
    rank_in_frag = max(0, min(dg_frag%id_frag, nworker - 1))
    base = ncol / nworker
    extra = mod(ncol, nworker)
    if (rank_in_frag < extra) then
      col_s = rank_in_frag * (base + 1) + 1
      col_e = col_s + base
    else
      col_s = extra * (base + 1) + (rank_in_frag - extra) * base + 1
      col_e = col_s + base - 1
    end if
    if (col_s > ncol) then
      col_s = 1
      col_e = 0
    else
      col_e = min(col_e, ncol)
    end if
  end subroutine get_fragment_column_range

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

  subroutine apply_gradient_to_basis_ops_local_2d(dg_frag, i_local, jo, mg, stencil, loc_s, loc_e, grad_phi, grad_local_2d)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer,                intent(in) :: i_local, jo
    type(s_rgrid),          intent(in) :: mg
    type(s_stencil),        intent(in) :: stencil
    integer,                intent(in) :: loc_s(3), loc_e(3)
    real(8),                intent(out) :: grad_phi(:,:,:,:)
    real(8),                intent(out) :: grad_local_2d(:,:)

    integer :: lx, ly, lz, ifrag, ndom(3), nloc1, nloc2, ipt
    integer :: gxg, gyg, gzg
    integer :: ix0, iy0, iz0
    integer :: p_lb1, p_ub1, p_lb2, p_ub2, p_lb3, p_ub3
    real(8) :: nabt(4,3), gx, gy, gz

    nabt = stencil%coef_nab
    ifrag = dg_frag%ifrag_start + i_local - 1
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub3 = ubound(dg_frag%phi_frag, 3)
    nloc1 = loc_e(1) - loc_s(1) + 1
    nloc2 = loc_e(2) - loc_s(2) + 1
    if (.not. allocated(dg_frag%phi_frag_has_seed_buffer) .or. &
        .not. dg_frag%phi_frag_has_seed_buffer(i_local)) then
      call enforce_fragment_periodic_buffer_for_state_ham(dg_frag, ifrag, i_local, jo)
    end if
    !$omp parallel do private(lx, ly, ipt, gxg, gyg, gzg, gx, gy, gz, ix0, iy0, iz0) schedule(static)
    do lz = 1, ndom(3)
      do ly = 1, ndom(2)
        !$omp simd private(ipt, gxg, gyg, gzg, gx, gy, gz, ix0, iy0, iz0)
        do lx = 1, ndom(1)
          gxg = dg_frag%ixyz_frag(1, ifrag) + lx - 1
          gyg = dg_frag%ixyz_frag(2, ifrag) + ly - 1
          gzg = dg_frag%ixyz_frag(3, ifrag) + lz - 1
          ix0 = map_global_to_phi_box_coord_ham(gxg, p_lb1, p_ub1, dg_frag%lgnum_total(1))
          iy0 = map_global_to_phi_box_coord_ham(gyg, p_lb2, p_ub2, dg_frag%lgnum_total(2))
          iz0 = map_global_to_phi_box_coord_ham(gzg, p_lb3, p_ub3, dg_frag%lgnum_total(3))
          if (ix0 == 0 .or. iy0 == 0 .or. iz0 == 0) then
            gx = 0.0d0
            gy = 0.0d0
            gz = 0.0d0
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
  ! Calculate momentum matrix elements in fragment basis (velocity gauge)
  !=======================================================================
  subroutine calculate_momentum_matrix(dg_frag, system, mg, stencil)
    use structures
    use communication, only: comm_is_root, comm_summation, comm_get_groupinfo, comm_get_min, comm_get_max, COMM_GROUP_NULL
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    
    integer :: ifrag, i_local, ispin, io, jo, idir, nbf, jo_progress_stride
    integer :: jloc, ncol_mom
    integer :: jo_s, jo_e
    integer :: ix, iy, iz, is(3), ie(3), i_halo, jfrag, n_basis_halo, ig_row, ig_col, ig_i, ig_j, l(3), d(3)
    integer :: lx, ly, lz, gx, gy, gz, iorg(3), ndom(3), loc_s(3), loc_e(3), phi_loc_s(3), phi_loc_e(3), halo_s(3), halo_e(3)
    integer :: lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi
    integer :: halo_send_idx(3), halo_recv_idx(3)
    integer :: phi_lb1, phi_ub1, phi_lb2, phi_ub2, phi_lb3, phi_ub3
    integer :: grad_lb1, grad_ub1, grad_lb2, grad_ub2, grad_lb3, grad_ub3
    integer :: iblk, iblk_rev, iblk_self, ii, jj, mat_size, ni, nj, ndiag
    integer :: npts_local, ipt, nx_local, ny_local
    logical :: log_frag_progress, has_overlap
    real(8) :: hvol, integral
    real(8) :: momentum_gb, momentum_block_gb
    real(8) :: max_p, pavg
    real(8) :: t0, t1
    real(8) :: time_halo_exchange, time_self_integral, time_halo_integral
    real(8) :: time_grad_total
    real(8) :: time_block_reduce, time_antisym
    real(8) :: time_reduce_pack, time_reduce_comm, time_reduce_unpack
    real(8) :: frag_grad_start, frag_self_start, frag_halo_start
    real(8), allocatable :: grad_phi(:,:,:,:)  ! gradient of basis function (x,y,z components, fragment-local)
    real(8), allocatable :: grad_local_2d(:,:), grad_all_2d(:,:,:), phi_local_2d(:,:), self_proj(:,:)
    real(8), allocatable :: mom_full(:, :, :, :)

    if (.not. dg_frag%has_real_space_basis) return
    ! Enforce fragment-local stencil policy: no halo communication path.
    dg_frag%n_halo = 0
    dg_frag%has_halo_exchange = .false.
    time_halo_exchange = 0.0d0
    time_self_integral = 0.0d0
    time_halo_integral = 0.0d0
    time_grad_total = 0.0d0
    time_block_reduce = 0.0d0
    time_antisym = 0.0d0
    time_reduce_pack = 0.0d0
    time_reduce_comm = 0.0d0
    time_reduce_unpack = 0.0d0
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "        Computing transition moments: <φ_i|∇|φ_j>"
      flush(6)
    end if
    
    ! Momentum matrix elements for vector potential coupling: p_ij = <phi_i|p|phi_j>
    ! In velocity gauge: H(t) = H_0 - i*A(t)·∇ + A(t)^2/2
    ! The A·p term couples to momentum matrix elements
    ! The A^2/2 term is diagonal (diamagnetic contribution)
    if (allocated(dg_frag%momentum_mat)) deallocate(dg_frag%momentum_mat)
    call init_momentum_blocks(dg_frag, diagonal_only=(.not. dg_frag%parallel_mode_orbital))
    momentum_gb = real(3_8 * int(dg_frag%n_mat_max, kind=8) * int(dg_frag%n_mat_max, kind=8) * &
      int(dg_frag%nspin, kind=8) * 8_8, 8) / 1.0d9
    momentum_block_gb = 0.0d0
    do iblk = 1, dg_frag%n_momentum_blocks
      momentum_block_gb = momentum_block_gb + real(3_8 * int(dg_frag%momentum_blocks(iblk)%nrow_max, kind=8) * &
        int(dg_frag%momentum_blocks(iblk)%ncol_max, kind=8) * int(dg_frag%nspin, kind=8) * 8_8, 8) / 1.0d9
    end do
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,f10.3,a)') "        n_mat_max=", dg_frag%n_mat_max, &
        " dense momentum_mat GB=", momentum_gb, " (not allocated)"
      write(*,'(1x,a,i0,a,f10.6,a)') "        momentum_blocks=", dg_frag%n_momentum_blocks, &
        " allocated GB=", momentum_block_gb, " per rank"
      flush(6)
    end if
    is = mg%is
    ie = mg%ie
    hvol = system%hvol
    
    ! Halo exchange removed: stencil operations use local phi_frag with fragment PBC buffer only.
    
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
        iorg(:) = dg_frag%ixyz_frag(:, ifrag)
        ndom(:) = dg_frag%nxyz_domain(:, ifrag)
        call get_fragment_owned_range(dg_frag, ifrag, mg, loc_s, loc_e, has_overlap)
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
        if (dg_frag%parallel_mode_orbital) then
          ! Momentum blocks are row-local in orbital mode.  Build the full
          ! ket-column range locally and avoid the old all-rank block reduce.
          jo_s = 1
          jo_e = nbf
        else
          call get_fragment_column_range(dg_frag, nbf, jo_s, jo_e)
        end if
        if (jo_s > jo_e) cycle
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
        ncol_mom = jo_e - jo_s + 1
        allocate(phi_local_2d(npts_local, nbf), grad_local_2d(npts_local, 3), &
                 grad_all_2d(npts_local, ncol_mom, 3), self_proj(nbf, ncol_mom))
        allocate(grad_phi(1:ndom(1), 1:ndom(2), 1:ndom(3), 3))
        grad_all_2d(:, :, :) = 0.0d0

        if (nbf > size(dg_frag%phi_frag, 4)) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "DG-Fragment RT invalid n_basis=", nbf, &
            " phi_frag dim4=", size(dg_frag%phi_frag, 4), " ifrag=", ifrag, " ispin=", ispin
          stop "DG-Fragment RT: invalid basis-function count"
        end if
        !$omp parallel do private(lz,ly,lx,ipt) schedule(static)
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

        ! Build all ket gradients first, then project in three BLAS-3 calls.
        ! This avoids nbf tiny DGEMM(nbf x 3) calls, which were dominating
        ! momentum timing as "self" on large DG runs.
        do jo = jo_s, jo_e
          jloc = jo - jo_s + 1
          call cpu_time(t0)
          call apply_gradient_to_basis_ops_local_2d(dg_frag, i_local, jo, mg, stencil, loc_s, loc_e, grad_phi, grad_local_2d)
          call cpu_time(t1)
          time_grad_total = time_grad_total + (t1 - t0)

          ig_j = dg_frag%index_basis(jo, ifrag, ispin)
          if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) then
            cycle
          end if
          grad_all_2d(:, jloc, 1) = grad_local_2d(:, 1)
          grad_all_2d(:, jloc, 2) = grad_local_2d(:, 2)
          grad_all_2d(:, jloc, 3) = grad_local_2d(:, 3)

          ig_col = ig_j
          if (ig_col < 1 .or. ig_col > dg_frag%n_mat_max) cycle
        end do  ! jo

        do idir = 1, 3
          call cpu_time(t0)
          call dgemm('T', 'N', nbf, ncol_mom, npts_local, hvol, phi_local_2d, npts_local, &
            grad_all_2d(1, 1, idir), npts_local, 0.0d0, self_proj, nbf)
          call cpu_time(t1)
          time_self_integral = time_self_integral + (t1 - t0)
          do io = 1, nbf
            ig_i = dg_frag%index_basis(io, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            do jloc = 1, ncol_mom
              jo = jo_s + jloc - 1
              integral = self_proj(io, jloc)
              if (iblk_self > 0) then
                if (io <= size(dg_frag%momentum_blocks(iblk_self)%val, 2) .and. &
                    jo <= size(dg_frag%momentum_blocks(iblk_self)%val, 3)) then
                  dg_frag%momentum_blocks(iblk_self)%val(idir, io, jo, ispin) = integral
                end if
              end if
            end do
          end do
        end do  ! idir
        if (allocated(grad_phi)) deallocate(grad_phi)
        deallocate(phi_local_2d, grad_local_2d, grad_all_2d, self_proj)
      end do  ! ifrag
    end do  ! ispin
    call add_dg_surface_momentum_blocks(dg_frag, system)
    
    ! Legacy real-space mode needs a global block aggregation.  Orbital mode
    ! builds the self block directly on the row-owning subgroup ranks.
    if (.not. dg_frag%parallel_mode_orbital) then
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

    ! Enforce anti-symmetry blockwise.  In orbital row-split mode self blocks
    ! are first gathered inside the fragment subgroup, because each rank stores
    ! only its bra rows.
    if (dg_frag%parallel_mode_orbital) then
      call cpu_time(t0)
      do ispin = 1, system%nspin
        do iblk = 1, dg_frag%n_momentum_blocks
          ifrag = dg_frag%momentum_blocks(iblk)%ifrag_row
          jfrag = dg_frag%momentum_blocks(iblk)%ifrag_col
          if (ifrag /= jfrag) cycle
          ndiag = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%momentum_blocks(iblk)%val, 2), &
                      size(dg_frag%momentum_blocks(iblk)%val, 3))
          if (ndiag <= 0) cycle
          allocate(mom_full(size(dg_frag%momentum_blocks(iblk)%val, 1), &
                            size(dg_frag%momentum_blocks(iblk)%val, 2), &
                            size(dg_frag%momentum_blocks(iblk)%val, 3), &
                            size(dg_frag%momentum_blocks(iblk)%val, 4)))
          if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
            call comm_summation(dg_frag%momentum_blocks(iblk)%val, mom_full, &
                                size(dg_frag%momentum_blocks(iblk)%val), dg_frag%icomm_frag)
          else
            mom_full(:, :, :, :) = dg_frag%momentum_blocks(iblk)%val(:, :, :, :)
          end if
          do idir = 1, 3
            do ii = 1, ndiag
              mom_full(idir, ii, ii, ispin) = 0.0d0
              do jj = ii + 1, ndiag
                pavg = 0.5d0 * (mom_full(idir, ii, jj, ispin) - mom_full(idir, jj, ii, ispin))
                mom_full(idir, ii, jj, ispin) = pavg
                mom_full(idir, jj, ii, ispin) = -pavg
              end do
            end do
            do ii = 1, ndiag
              if (basis_row_is_locally_owned(dg_frag, ifrag, ispin, ii)) then
                dg_frag%momentum_blocks(iblk)%val(idir, ii, 1:ndiag, ispin) = mom_full(idir, ii, 1:ndiag, ispin)
              else
                dg_frag%momentum_blocks(iblk)%val(idir, ii, 1:ndiag, ispin) = 0.0d0
              end if
            end do
          end do
          deallocate(mom_full)
        end do
      end do
      call cpu_time(t1)
      time_antisym = time_antisym + (t1 - t0)
    else
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
    end if
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,1pe12.4)') "        momentum antisym done time=", time_antisym
      flush(6)
    end if

    max_p = 0.0d0
    do iblk = 1, dg_frag%n_momentum_blocks
      max_p = max(max_p, maxval(abs(dg_frag%momentum_blocks(iblk)%val)))
    end do
    if (comm_is_root(dg_frag%id)) then
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
    
  end subroutine calculate_momentum_matrix

  !=======================================================================
  ! Calculate overlap matrix in DG basis (S_ij = <phi_i|phi_j>)
  !=======================================================================
  subroutine calculate_overlap_matrix(dg_frag, system, mg)
    use structures
    use communication, only: comm_is_root, comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg

    integer :: ifrag, jfrag, i_local, ispin, io, jo, iblk, nbf, jo_progress_stride
    integer :: jo_s, jo_e, jo_loc, ncol_local
    integer :: ix, iy, iz, is(3), ie(3), i_halo
    integer :: ig_row, ig_col, d(3), ii, jj
    integer :: lx, ly, lz, gx, gy, gz, iorg(3), ndom(3), loc_s(3), loc_e(3), halo_s(3), halo_e(3)
    integer :: npts_local, nx_local, ny_local, ipt
    integer :: phi_loc_s(3), phi_loc_e(3)
    integer :: lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi
    integer :: phi_lb1, phi_lb2, phi_lb3, phi_ub1, phi_ub2, phi_ub3
    integer :: buf_lb1, buf_lb2, buf_lb3, buf_ub1, buf_ub2, buf_ub3
    logical :: log_frag_progress, release_dense_overlap, has_overlap
    real(8) :: hvol, integral, savg, s_min, s_max, cond_est
    real(8) :: t0, t1, time_self_integral, time_halo_integral, time_reduce_total
    real(8) :: frag_self_start, frag_halo_start
    real(8), allocatable :: phi_local_2d(:,:), self_overlap(:,:)

    if (.not. dg_frag%has_real_space_basis) return
    if (.not. allocated(dg_frag%index_basis) .or. .not. allocated(dg_frag%n_mat)) return

    ! Enforce fragment-local stencil policy: no halo communication path.
    dg_frag%n_halo = 0
    dg_frag%has_halo_exchange = .false.

    release_dense_overlap = (.not. dg_frag%yn_adaptive_basis) .and. &
      ((.not. dg_frag%use_plane_wave_basis) .or. dg_frag%n_plane_waves <= 0)

    if (allocated(dg_frag%S_mat)) deallocate(dg_frag%S_mat)
    if (allocated(dg_frag%S_mat_prop)) deallocate(dg_frag%S_mat_prop)
    dg_frag%has_global_overlap_copy = .true.
    dg_frag%overlap_prop_root_authoritative = .false.
    call init_matrix_blocks(dg_frag, dg_frag%S_mat_blocks, dg_frag%S_block_map, dg_frag%n_S_blocks, &
                            diagonal_only=.true.)
    if (allocated(dg_frag%S_mat_blocks)) then
      do i_halo = 1, size(dg_frag%S_mat_blocks)
        dg_frag%S_mat_blocks(i_halo)%val(:, :, :) = 0.0d0
      end do
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
        if (dg_frag%parallel_mode_orbital) then
          ! Row-local overlap blocks need all ket columns on the owning rank;
          ! this avoids replicating the full S block set through allreduce.
          jo_s = 1
          jo_e = nbf
        else
          call get_fragment_column_range(dg_frag, nbf, jo_s, jo_e)
        end if
        ncol_local = max(0, jo_e - jo_s + 1)
        if (ncol_local <= 0) cycle
        npts_local = (lx_hi - lx_lo + 1) * (ly_hi - ly_lo + 1) * (lz_hi - lz_lo + 1)
        nx_local = lx_hi - lx_lo + 1
        ny_local = ly_hi - ly_lo + 1
        allocate(phi_local_2d(npts_local, nbf), self_overlap(nbf, ncol_local))
        !$omp parallel do private(lz,ly,lx,ipt) schedule(static)
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
        call dgemm('T', 'N', nbf, ncol_local, npts_local, hvol, phi_local_2d, npts_local, &
                   phi_local_2d(1, jo_s), npts_local, 0.0d0, self_overlap, nbf)
        call cpu_time(t1)
        time_self_integral = time_self_integral + (t1 - t0)
        iblk = find_matrix_block(dg_frag%S_block_map, ifrag, ifrag)

        do jo = jo_s, jo_e
          jo_loc = jo - jo_s + 1
          ig_col = dg_frag%index_basis(jo, ifrag, ispin)
          if (ig_col < 1 .or. ig_col > dg_frag%n_mat_max) cycle
          if (iblk <= 0) cycle
          if (jo > size(dg_frag%S_mat_blocks(iblk)%val, 2)) cycle

          do io = 1, nbf
            ig_row = dg_frag%index_basis(io, ifrag, ispin)
            if (ig_row < 1 .or. ig_row > dg_frag%n_mat_max) cycle
            integral = self_overlap(io, jo_loc)
            if (io <= size(dg_frag%S_mat_blocks(iblk)%val, 1)) then
              dg_frag%S_mat_blocks(iblk)%val(io, jo, ispin) = integral
            end if
          end do

        end do
        deallocate(phi_local_2d, self_overlap)
      end do
    end do

    call cpu_time(t0)
    ! Legacy ranks can hold partial real-space boxes; orbital mode already
    ! integrates the full fragment box for its local row block.
    if (.not. dg_frag%parallel_mode_orbital) then
      call reduce_matrix_blocks(dg_frag, dg_frag%S_mat_blocks, "smat-frag", dg_frag%icomm_frag)
    end if
    if (.not. dg_frag%is_frag_root) then
      if (dg_frag%parallel_mode_orbital) then
        ! Keep the replicated overlap block for row-owned coefficient updates.
        continue
      else
        if (allocated(dg_frag%S_mat_blocks)) then
          do i_halo = 1, size(dg_frag%S_mat_blocks)
            if (fragment_row_is_locally_owned(dg_frag, dg_frag%S_mat_blocks(i_halo)%ifrag_row)) cycle
            dg_frag%S_mat_blocks(i_halo)%val(:, :, :) = 0.0d0
          end do
        end if
      end if
    end if
    call cpu_time(t1)
    time_reduce_total = time_reduce_total + (t1 - t0)
    write(*,'(1x,a,1pe12.4)') "        overlap reduce done time=", time_reduce_total
    flush(6)

    call init_matrix_blocks(dg_frag, dg_frag%S_mat_prop_blocks, dg_frag%S_block_map, dg_frag%n_S_blocks, &
                            diagonal_only=.true.)
    if (allocated(dg_frag%S_mat_blocks) .and. allocated(dg_frag%S_mat_prop_blocks)) then
      do i_halo = 1, size(dg_frag%S_mat_prop_blocks)
        ifrag = dg_frag%S_mat_prop_blocks(i_halo)%ifrag_row
        jfrag = dg_frag%S_mat_prop_blocks(i_halo)%ifrag_col
        iblk = 0
        if (allocated(dg_frag%S_block_map)) then
          if (ifrag >= 1 .and. ifrag <= size(dg_frag%S_block_map, 1) .and. &
              jfrag >= 1 .and. jfrag <= size(dg_frag%S_block_map, 2)) &
            iblk = dg_frag%S_block_map(ifrag, jfrag)
        end if
        if (iblk > 0 .and. iblk <= size(dg_frag%S_mat_blocks)) &
          dg_frag%S_mat_prop_blocks(i_halo)%val(:, :, :) = dg_frag%S_mat_blocks(iblk)%val(:, :, :)
      end do
    end if

    do ispin = 1, dg_frag%nspin
      if (.not. allocated(dg_frag%S_mat_blocks)) cycle
      do iblk = 1, size(dg_frag%S_mat_blocks)
        ifrag = dg_frag%S_mat_blocks(iblk)%ifrag_row
        if (ifrag < 1 .or. ifrag > dg_frag%n_frag) cycle
        if (dg_frag%S_mat_blocks(iblk)%ifrag_col /= ifrag) cycle
        if (.not. allocated(dg_frag%S_mat_blocks(iblk)%val)) cycle
        do ii = 1, min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%S_mat_blocks(iblk)%val, 1), &
                       size(dg_frag%S_mat_blocks(iblk)%val, 2))
          if (dg_frag%S_mat_blocks(iblk)%val(ii, ii, ispin) < 1.0d-12) then
            dg_frag%S_mat_blocks(iblk)%val(ii, ii, ispin) = 1.0d0
            if (allocated(dg_frag%S_mat_prop_blocks)) then
              i_halo = find_matrix_block(dg_frag%S_block_map, ifrag, ifrag)
              if (i_halo > 0 .and. i_halo <= size(dg_frag%S_mat_prop_blocks)) &
                dg_frag%S_mat_prop_blocks(i_halo)%val(ii, ii, ispin) = 1.0d0
            end if
          end if
        end do
      end do
    end do

    dg_frag%has_global_overlap_copy = .false.
    dg_frag%overlap_prop_root_authoritative = .false.

  end subroutine calculate_overlap_matrix
