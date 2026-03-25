module rt_dg_fragment_ops
  use communication, only: comm_bcast, comm_get_max, comm_is_root, comm_summation, COMM_GROUP_NULL
  use rt_dg_fragment_types, only: s_dg_fragment_rt, matrix_block_info
  implicit none

  private
  public :: ensure_nonlocal_pp_matrix_A
  public :: ensure_overlap_prop_available
  public :: calculate_microscopic_current_dg
  public :: build_spatial_A_coupling_matrices
  public :: apply_gradient_to_basis
  public :: apply_momentum_blocks
  public :: apply_matrix_blocks
  public :: apply_matrix_blocks_batch
  public :: copy_matrix_blocks_to_complex_dense
  public :: pack_owned_coef
  public :: fetch_remote_coef_rows
  public :: pack_owned_coef_pw
  public :: fetch_remote_coef_pw_rows
  public :: gather_full_coef_view
  public :: zero_nonowned_coefficients

contains

  !=======================================================================
  ! Build non-local pseudopotential matrix with vector potential A(t)
  ! V_NL(A) = e^{-i A·r} V_NL e^{i A·r}
  !
  ! Uses the SALMON approximation: A is nearly constant within PP cutoff
  !=======================================================================
  subroutine build_nonlocal_pp_matrix_A(dg_frag, mg, ppg, nspin, hvol, Ac_tot, use_micro_A, Ac_micro, H_nl)
    use structures
    use math_constants, only: zi
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    integer,                intent(in) :: nspin
    real(8),                intent(in) :: hvol
    real(8),                intent(in) :: Ac_tot(3)
    logical,                intent(in) :: use_micro_A
    real(8),                intent(in), optional :: Ac_micro(:, :, :, :)
    complex(8),             intent(out) :: H_nl(dg_frag%n_mat_max, dg_frag%n_mat_max, nspin)

    integer :: ifrag, ispin, io, jo, i_local, ilma, ia, j, ix, iy, iz, ig_i, ig_j, nbf
    integer :: iorg(3), ndom(3), lx, ly, lz
    integer :: is(3), ie(3)
    real(8) :: x, y, z, phase
    real(8) :: A_local(3)
    complex(8), allocatable :: uVpsi(:,:,:)  ! (nstate_frag, Nlma, nspin)
    complex(8) :: overlap_i, overlap_j, nlpp_contrib

    if (ppg%Nlma == 0) then
      H_nl = (0.0d0, 0.0d0)
      return
    end if

    is = mg%is
    ie = mg%ie
    H_nl = (0.0d0, 0.0d0)

    allocate(uVpsi(dg_frag%nstate_frag, ppg%Nlma, nspin))

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      iorg(:) = dg_frag%ixyz_frag(:, ifrag)
      ndom(:) = dg_frag%nxyz_domain(:, ifrag)
      uVpsi = (0.0d0, 0.0d0)

      do ispin = 1, nspin
        !$omp parallel do collapse(2) private(ilma, io, ia, j, ix, iy, iz, x, y, z, phase, overlap_i)
        do ilma = 1, ppg%Nlma
          do io = 1, dg_frag%nstate_frag
            ia = ppg%ia_tbl(ilma)
            overlap_i = (0.0d0, 0.0d0)
            do j = 1, ppg%mps(ia)
              ix = ppg%jxyz(1, j, ia)
              iy = ppg%jxyz(2, j, ia)
              iz = ppg%jxyz(3, j, ia)

              if (ix >= is(1) .and. ix <= ie(1) .and. &
                  iy >= is(2) .and. iy <= ie(2) .and. &
                  iz >= is(3) .and. iz <= ie(3)) then
                lx = mod(ix - iorg(1), dg_frag%lgnum_total(1)) + 1
                ly = mod(iy - iorg(2), dg_frag%lgnum_total(2)) + 1
                lz = mod(iz - iorg(3), dg_frag%lgnum_total(3)) + 1
                if (lx < 1 .or. lx > ndom(1)) cycle
                if (ly < 1 .or. ly > ndom(2)) cycle
                if (lz < 1 .or. lz > ndom(3)) cycle
                x = ppg%rxyz(1, j, ia)
                y = ppg%rxyz(2, j, ia)
                z = ppg%rxyz(3, j, ia)
                if (use_micro_A .and. present(Ac_micro)) then
                  A_local(1:3) = Ac_micro(1:3, ix, iy, iz)
                else
                  A_local(1:3) = Ac_tot(1:3)
                end if
                phase = A_local(1) * x + A_local(2) * y + A_local(3) * z
                overlap_i = overlap_i + dg_frag%phi_frag(lx, ly, lz, io, i_local) * &
                            ppg%uV(j, ilma) * exp(-zi * phase) * hvol
              end if
            end do
            uVpsi(io, ilma, ispin) = overlap_i
          end do
        end do
        !$omp end parallel do

        nbf = dg_frag%n_basis(ifrag, ispin)
        !$omp parallel do collapse(2) private(jo, io, ig_i, ig_j, ilma, nlpp_contrib, overlap_i, overlap_j)
        do jo = 1, nbf
          do io = 1, nbf
            ig_i = dg_frag%index_basis(io, ifrag, ispin)
            ig_j = dg_frag%index_basis(jo, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle
            nlpp_contrib = (0.0d0, 0.0d0)
            do ilma = 1, ppg%Nlma
              overlap_i = uVpsi(io, ilma, ispin)
              overlap_j = uVpsi(jo, ilma, ispin)
              nlpp_contrib = nlpp_contrib + overlap_i * overlap_j * ppg%rinv_uvu(ilma)
            end do
            H_nl(ig_i, ig_j, ispin) = H_nl(ig_i, ig_j, ispin) + nlpp_contrib
          end do
        end do
        !$omp end parallel do
      end do
    end do

    deallocate(uVpsi)

  end subroutine build_nonlocal_pp_matrix_A

  !=======================================================================
  ! Ensure cached non-local PP matrix for current A(t)
  !=======================================================================
  subroutine ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_tot)
    use structures
    use salmon_global, only: ae_shape1, ae_shape2, theory
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    type(s_dft_system),     intent(in)    :: system
    real(8),                intent(in)    :: Ac_tot(3)

    real(8) :: delta_A
    logical :: reuse_allowed
    logical :: use_micro_A

    if (ppg%Nlma == 0 .or. .not. allocated(ppg%uV)) then
      if (allocated(dg_frag%H_nl_cache)) deallocate(dg_frag%H_nl_cache)
      dg_frag%has_nl_cache = .false.
      return
    end if

    if (.not. allocated(dg_frag%H_nl_cache)) then
      allocate(dg_frag%H_nl_cache(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
      dg_frag%has_nl_cache = .false.
    end if

    use_micro_A = (trim(theory) == 'single_scale_maxwell_tddft' .and. allocated(system%Ac_micro%v))

    reuse_allowed = (.not. use_micro_A) .and. &
                    (trim(ae_shape1) == 'impulse' .or. trim(ae_shape1) == 'none') .and. &
                    (trim(ae_shape2) == 'impulse' .or. trim(ae_shape2) == 'none')

    delta_A = maxval(abs(Ac_tot - dg_frag%Ac_nl_cache))
    if (.not. dg_frag%has_nl_cache .or. (.not. reuse_allowed) .or. delta_A > dg_frag%Ac_nl_cache_tol) then
      if (use_micro_A) then
        call build_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system%nspin, system%hvol, Ac_tot, &
             .true., system%Ac_micro%v, dg_frag%H_nl_cache)
      else
        call build_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system%nspin, system%hvol, Ac_tot, &
             .false., H_nl=dg_frag%H_nl_cache)
      end if
      dg_frag%Ac_nl_cache = Ac_tot
      dg_frag%has_nl_cache = .true.
    end if

  end subroutine ensure_nonlocal_pp_matrix_A

  subroutine ensure_overlap_prop_available(dg_frag, n_use)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: n_use

    integer :: n_copy

    if (n_use <= 0) return
    if (.not. dg_frag%overlap_prop_root_authoritative) return
    if (dg_frag%icomm_frag == COMM_GROUP_NULL) return

    n_copy = min(n_use, dg_frag%n_mat_max)
    if (n_copy <= 0) return

    if (allocated(dg_frag%S_mat_prop_c)) then
      call comm_bcast(dg_frag%S_mat_prop_c(1:n_copy, 1:n_copy, 1:dg_frag%nspin), dg_frag%icomm_frag, 0)
    end if
    if (allocated(dg_frag%S_mat_prop)) then
      call comm_bcast(dg_frag%S_mat_prop(1:n_copy, 1:n_copy, 1:dg_frag%nspin), dg_frag%icomm_frag, 0)
    end if
    if (allocated(dg_frag%S_mat_c)) then
      call comm_bcast(dg_frag%S_mat_c(1:n_copy, 1:n_copy, 1:dg_frag%nspin), dg_frag%icomm_frag, 0)
    end if
    if (allocated(dg_frag%S_mat)) then
      call comm_bcast(dg_frag%S_mat(1:n_copy, 1:n_copy, 1:dg_frag%nspin), dg_frag%icomm_frag, 0)
    end if
    if (allocated(dg_frag%S_mat_prop_blocks) .and. allocated(dg_frag%S_mat_prop) .and. &
        allocated(dg_frag%S_block_map)) then
      call sync_dense_matrix_to_blocks_runtime(dg_frag, dg_frag%S_mat_prop, dg_frag%S_mat_prop_blocks, dg_frag%S_block_map)
    end if
    if (allocated(dg_frag%S_mat_blocks) .and. allocated(dg_frag%S_mat) .and. &
        allocated(dg_frag%S_block_map)) then
      call sync_dense_matrix_to_blocks_runtime(dg_frag, dg_frag%S_mat, dg_frag%S_mat_blocks, dg_frag%S_block_map)
    end if
  end subroutine ensure_overlap_prop_available

  subroutine sync_dense_matrix_to_blocks_runtime(dg_frag, mat, blocks, block_map)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    real(8), intent(in) :: mat(:, :, :)
    type(matrix_block_info), intent(inout) :: blocks(:)
    integer, intent(in) :: block_map(:, :)

    integer :: ifrag_row, ifrag_col, iblk, ispin, ii, jj, ig_i, ig_j
    integer :: nrow, ncol

    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        iblk = find_matrix_block_runtime(block_map, ifrag_row, ifrag_col)
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
  end subroutine sync_dense_matrix_to_blocks_runtime

  subroutine init_matrix_blocks_runtime(dg_frag, blocks, block_map)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), allocatable, intent(inout) :: blocks(:)
    integer, allocatable, intent(inout) :: block_map(:, :)

    integer :: ifrag_row, ifrag_col, iblk, n_blocks, nrow_max, ncol_max

    if (allocated(blocks)) then
      do iblk = 1, size(blocks)
        if (allocated(blocks(iblk)%val)) deallocate(blocks(iblk)%val)
      end do
      deallocate(blocks)
    end if
    if (allocated(block_map)) deallocate(block_map)

    n_blocks = dg_frag%n_frag * dg_frag%n_frag
    if (n_blocks <= 0) return

    allocate(blocks(n_blocks))
    allocate(block_map(dg_frag%n_frag, dg_frag%n_frag))
    block_map = 0

    iblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
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
  end subroutine init_matrix_blocks_runtime

  subroutine sync_blocks_to_dense_matrix_runtime(dg_frag, blocks, block_map, mat)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(in) :: blocks(:)
    integer, intent(in) :: block_map(:, :)
    real(8), intent(inout) :: mat(:, :, :)

    integer :: ifrag_row, ifrag_col, iblk, ispin, ii, jj, ig_i, ig_j
    integer :: nrow, ncol

    mat(:, :, :) = 0.0d0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        iblk = find_matrix_block_runtime(block_map, ifrag_row, ifrag_col)
        if (iblk <= 0) cycle
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
              mat(ig_i, ig_j, ispin) = blocks(iblk)%val(ii, jj, ispin)
            end do
          end do
        end do
      end do
    end do
  end subroutine sync_blocks_to_dense_matrix_runtime

  subroutine reduce_matrix_blocks_runtime(dg_frag, blocks, label, icomm_reduce)
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
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] Matrix block size mismatch: label=", &
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
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0)') "        hamiltonian block reduce done: label=", trim(label), &
        " total_active=", total_active_size
      flush(6)
    end if
  end subroutine reduce_matrix_blocks_runtime

  integer function find_matrix_block_runtime(block_map, ifrag_row, ifrag_col) result(iblk)
    implicit none
    integer, intent(in) :: block_map(:, :)
    integer, intent(in) :: ifrag_row, ifrag_col

    iblk = 0
    if (ifrag_row < 1 .or. ifrag_row > size(block_map, 1)) return
    if (ifrag_col < 1 .or. ifrag_col > size(block_map, 2)) return
    iblk = block_map(ifrag_row, ifrag_col)
  end function find_matrix_block_runtime

  subroutine pack_owned_coef(dg_frag, ispin, row_ids, packed)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    integer, intent(in) :: row_ids(:)
    complex(8), intent(out) :: packed(:, :)

    integer :: irow, global_row

    packed(:, :) = (0.0d0, 0.0d0)
    if (.not. allocated(dg_frag%coef_owner)) return
    if (ispin < 1 .or. ispin > dg_frag%nspin) return

    do irow = 1, min(size(row_ids), size(packed, 1))
      global_row = row_ids(irow)
      if (global_row < 1 .or. global_row > size(dg_frag%coef, 1)) cycle
      if (dg_frag%coef_owner(global_row, ispin) /= dg_frag%id) cycle
      packed(irow, 1:size(packed, 2)) = dg_frag%coef(global_row, 1:size(packed, 2), ispin)
    end do
  end subroutine pack_owned_coef

  subroutine fetch_remote_coef_rows(dg_frag, ispin, row_ids, fetched)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    integer, intent(in) :: row_ids(:)
    complex(8), intent(out) :: fetched(:, :)

    integer :: irow, global_row, owner_rank
    complex(8), allocatable :: row_buf(:)

    fetched(:, :) = (0.0d0, 0.0d0)
    if (.not. allocated(dg_frag%coef_owner)) return
    if (ispin < 1 .or. ispin > dg_frag%nspin) return

    allocate(row_buf(size(fetched, 2)))
    do irow = 1, min(size(row_ids), size(fetched, 1))
      global_row = row_ids(irow)
      row_buf(:) = (0.0d0, 0.0d0)
      if (global_row < 1 .or. global_row > size(dg_frag%coef, 1)) then
        fetched(irow, :) = row_buf(:)
        cycle
      end if
      owner_rank = dg_frag%coef_owner(global_row, ispin)
      if (owner_rank < 0) then
        fetched(irow, :) = row_buf(:)
        cycle
      end if
      if (owner_rank == dg_frag%id) row_buf(:) = dg_frag%coef(global_row, 1:size(fetched, 2), ispin)
      call comm_bcast(row_buf, dg_frag%icomm, owner_rank)
      fetched(irow, :) = row_buf(:)
    end do
    deallocate(row_buf)
  end subroutine fetch_remote_coef_rows

  subroutine pack_owned_coef_pw(dg_frag, row_ids, packed)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: row_ids(:)
    complex(8), intent(out) :: packed(:, :, :)

    integer :: irow, pw_row

    packed(:, :, :) = (0.0d0, 0.0d0)
    if (.not. allocated(dg_frag%coef_pw_owner)) return
    if (.not. allocated(dg_frag%coef_pw)) return

    do irow = 1, min(size(row_ids), size(packed, 1))
      pw_row = row_ids(irow)
      if (pw_row < 1 .or. pw_row > size(dg_frag%coef_pw, 1)) cycle
      if (dg_frag%coef_pw_owner(pw_row) /= dg_frag%id) cycle
      packed(irow, 1:size(packed, 2), 1:size(packed, 3)) = &
        dg_frag%coef_pw(pw_row, 1:size(packed, 2), 1:size(packed, 3))
    end do
  end subroutine pack_owned_coef_pw

  subroutine fetch_remote_coef_pw_rows(dg_frag, row_ids, fetched)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: row_ids(:)
    complex(8), intent(out) :: fetched(:, :, :)

    integer :: irow, pw_row, owner_rank
    complex(8), allocatable :: row_buf(:, :)

    fetched(:, :, :) = (0.0d0, 0.0d0)
    if (.not. allocated(dg_frag%coef_pw_owner)) return
    if (.not. allocated(dg_frag%coef_pw)) return

    allocate(row_buf(size(fetched, 2), size(fetched, 3)))
    do irow = 1, min(size(row_ids), size(fetched, 1))
      pw_row = row_ids(irow)
      row_buf(:, :) = (0.0d0, 0.0d0)
      if (pw_row < 1 .or. pw_row > size(dg_frag%coef_pw, 1)) then
        fetched(irow, :, :) = row_buf(:, :)
        cycle
      end if
      owner_rank = dg_frag%coef_pw_owner(pw_row)
      if (owner_rank < 0) then
        fetched(irow, :, :) = row_buf(:, :)
        cycle
      end if
      if (owner_rank == dg_frag%id) row_buf(:, :) = dg_frag%coef_pw(pw_row, 1:size(fetched, 2), 1:size(fetched, 3))
      call comm_bcast(row_buf, dg_frag%icomm, owner_rank)
      fetched(irow, :, :) = row_buf(:, :)
    end do
    deallocate(row_buf)
  end subroutine fetch_remote_coef_pw_rows

  subroutine gather_full_coef_view(dg_frag, ispin, n_frag_rows, nstate_use, coef_frag, coef_pw)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ispin
    integer, intent(in) :: n_frag_rows
    integer, intent(in) :: nstate_use
    complex(8), allocatable, intent(out) :: coef_frag(:,:)
    complex(8), allocatable, intent(out) :: coef_pw(:,:)

    integer :: i, n_pw
    integer, allocatable :: frag_row_ids(:), pw_row_ids(:)
    complex(8), allocatable :: coef_pw_all(:,:,:)

    allocate(coef_frag(max(0, n_frag_rows), max(0, nstate_use)))
    coef_frag(:, :) = (0.0d0, 0.0d0)
    if (n_frag_rows > 0 .and. nstate_use > 0) then
      allocate(frag_row_ids(n_frag_rows))
      do i = 1, n_frag_rows
        frag_row_ids(i) = i
      end do
      call fetch_remote_coef_rows(dg_frag, ispin, frag_row_ids, coef_frag)
      deallocate(frag_row_ids)
    end if

    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    allocate(coef_pw(max(0, n_pw), max(0, nstate_use)))
    coef_pw(:, :) = (0.0d0, 0.0d0)
    if (n_pw > 0 .and. nstate_use > 0) then
      allocate(pw_row_ids(n_pw))
      do i = 1, n_pw
        pw_row_ids(i) = i
      end do
      allocate(coef_pw_all(n_pw, nstate_use, dg_frag%nspin))
      coef_pw_all(:, :, :) = (0.0d0, 0.0d0)
      call fetch_remote_coef_pw_rows(dg_frag, pw_row_ids, coef_pw_all)
      coef_pw(:, :) = coef_pw_all(:, :, min(max(ispin, 1), dg_frag%nspin))
      deallocate(coef_pw_all)
      deallocate(pw_row_ids)
    end if
  end subroutine gather_full_coef_view

  subroutine zero_nonowned_coefficients(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    integer :: ispin, i

    if (allocated(dg_frag%coef_owner) .and. allocated(dg_frag%coef)) then
      do ispin = 1, min(size(dg_frag%coef_owner, 2), size(dg_frag%coef, 3))
        do i = 1, min(size(dg_frag%coef_owner, 1), size(dg_frag%coef, 1))
          if (dg_frag%coef_owner(i, ispin) == dg_frag%id) cycle
          dg_frag%coef(i, :, ispin) = (0.0d0, 0.0d0)
        end do
      end do
    end if

    if (allocated(dg_frag%coef_pw_owner) .and. allocated(dg_frag%coef_pw)) then
      do i = 1, min(size(dg_frag%coef_pw_owner), size(dg_frag%coef_pw, 1))
        if (dg_frag%coef_pw_owner(i) == dg_frag%id) cycle
        dg_frag%coef_pw(i, :, :) = (0.0d0, 0.0d0)
      end do
    end if
  end subroutine zero_nonowned_coefficients

  subroutine apply_matrix_blocks(dg_frag, blocks, ispin, x, y)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(in) :: blocks(:)
    integer, intent(in) :: ispin
    complex(8), intent(in) :: x(:)
    complex(8), intent(inout) :: y(:)

    integer :: iblk, ifrag_row, ifrag_col
    integer :: nrow, ncol, ii, jj, ig_i, ig_j
    complex(8) :: xj

    if (ispin < 1 .or. ispin > dg_frag%nspin) return
    if (.not. allocated(dg_frag%index_basis)) return

    do iblk = 1, size(blocks)
      ifrag_row = blocks(iblk)%ifrag_row
      ifrag_col = blocks(iblk)%ifrag_col
      if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) cycle
      if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
      nrow = dg_frag%n_basis(ifrag_row, ispin)
      ncol = dg_frag%n_basis(ifrag_col, ispin)
      if (nrow <= 0 .or. ncol <= 0) cycle

      do jj = 1, ncol
        ig_j = dg_frag%index_basis(jj, ifrag_col, ispin)
        if (ig_j < 1 .or. ig_j > size(x)) cycle
        xj = x(ig_j)
        if (abs(xj) == 0.0d0) cycle
        do ii = 1, nrow
          ig_i = dg_frag%index_basis(ii, ifrag_row, ispin)
          if (ig_i < 1 .or. ig_i > size(y)) cycle
          y(ig_i) = y(ig_i) + blocks(iblk)%val(ii, jj, ispin) * xj
        end do
      end do
    end do
  end subroutine apply_matrix_blocks

  subroutine apply_matrix_blocks_batch(dg_frag, blocks, ispin, x, y)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(in) :: blocks(:)
    integer, intent(in) :: ispin
    complex(8), intent(in) :: x(:, :)
    complex(8), intent(inout) :: y(:, :)

    integer :: iblk, ifrag_row, ifrag_col
    integer :: nrow, ncol, ii, jj, ist, ig_i, ig_j
    complex(8) :: xj

    if (ispin < 1 .or. ispin > dg_frag%nspin) return
    if (.not. allocated(dg_frag%index_basis)) return

    do iblk = 1, size(blocks)
      ifrag_row = blocks(iblk)%ifrag_row
      ifrag_col = blocks(iblk)%ifrag_col
      if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) cycle
      if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
      nrow = dg_frag%n_basis(ifrag_row, ispin)
      ncol = dg_frag%n_basis(ifrag_col, ispin)
      if (nrow <= 0 .or. ncol <= 0) cycle

      do jj = 1, ncol
        ig_j = dg_frag%index_basis(jj, ifrag_col, ispin)
        if (ig_j < 1 .or. ig_j > size(x, 1)) cycle
        do ist = 1, size(x, 2)
          xj = x(ig_j, ist)
          if (abs(xj) == 0.0d0) cycle
          do ii = 1, nrow
            ig_i = dg_frag%index_basis(ii, ifrag_row, ispin)
            if (ig_i < 1 .or. ig_i > size(y, 1)) cycle
            y(ig_i, ist) = y(ig_i, ist) + blocks(iblk)%val(ii, jj, ispin) * xj
          end do
        end do
      end do
    end do
  end subroutine apply_matrix_blocks_batch

  subroutine copy_matrix_blocks_to_complex_dense(dg_frag, blocks, ispin, mat)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(in) :: blocks(:)
    integer, intent(in) :: ispin
    complex(8), intent(inout) :: mat(:, :)

    integer :: iblk, ifrag_row, ifrag_col
    integer :: nrow, ncol, ii, jj, ig_i, ig_j

    if (ispin < 1 .or. ispin > dg_frag%nspin) return
    if (.not. allocated(dg_frag%index_basis)) return

    do iblk = 1, size(blocks)
      ifrag_row = blocks(iblk)%ifrag_row
      ifrag_col = blocks(iblk)%ifrag_col
      if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) cycle
      if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
      nrow = dg_frag%n_basis(ifrag_row, ispin)
      ncol = dg_frag%n_basis(ifrag_col, ispin)
      if (nrow <= 0 .or. ncol <= 0) cycle

      do jj = 1, ncol
        ig_j = dg_frag%index_basis(jj, ifrag_col, ispin)
        if (ig_j < 1 .or. ig_j > size(mat, 2)) cycle
        do ii = 1, nrow
          ig_i = dg_frag%index_basis(ii, ifrag_row, ispin)
          if (ig_i < 1 .or. ig_i > size(mat, 1)) cycle
          mat(ig_i, ig_j) = cmplx(blocks(iblk)%val(ii, jj, ispin), 0.0d0, kind=8)
        end do
      end do
    end do
  end subroutine copy_matrix_blocks_to_complex_dense

  !=======================================================================
  ! Apply gradient operator to a basis function using finite differences
  !=======================================================================
  subroutine apply_gradient_to_basis(dg_frag, i_local, jo, mg, stencil, grad_phi)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer,                intent(in) :: i_local, jo
    type(s_rgrid),          intent(in) :: mg
    type(s_stencil),        intent(in) :: stencil
    real(8),                intent(out) :: grad_phi(:,:,:,:)

    integer :: lx, ly, lz, ifrag
    real(8) :: nabt(4,3)
    integer :: ndom(3)

    nabt = stencil%coef_nab
    ifrag = dg_frag%ifrag_start + i_local - 1
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)

    grad_phi = 0.0d0

    !$omp parallel do collapse(3) private(lx, ly, lz) schedule(static)
    do lz = 1, ndom(3)
      do ly = 1, ndom(2)
        do lx = 1, ndom(1)
          grad_phi(lx, ly, lz, 1) = &
              nabt(1,1) * (dg_frag%phi_frag(lx+1, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-1, ly, lz, jo, i_local)) + &
              nabt(2,1) * (dg_frag%phi_frag(lx+2, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-2, ly, lz, jo, i_local)) + &
              nabt(3,1) * (dg_frag%phi_frag(lx+3, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-3, ly, lz, jo, i_local)) + &
              nabt(4,1) * (dg_frag%phi_frag(lx+4, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-4, ly, lz, jo, i_local))

          grad_phi(lx, ly, lz, 2) = &
              nabt(1,2) * (dg_frag%phi_frag(lx, ly+1, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-1, lz, jo, i_local)) + &
              nabt(2,2) * (dg_frag%phi_frag(lx, ly+2, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-2, lz, jo, i_local)) + &
              nabt(3,2) * (dg_frag%phi_frag(lx, ly+3, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-3, lz, jo, i_local)) + &
              nabt(4,2) * (dg_frag%phi_frag(lx, ly+4, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-4, lz, jo, i_local))

          grad_phi(lx, ly, lz, 3) = &
              nabt(1,3) * (dg_frag%phi_frag(lx, ly, lz+1, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-1, jo, i_local)) + &
              nabt(2,3) * (dg_frag%phi_frag(lx, ly, lz+2, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-2, jo, i_local)) + &
              nabt(3,3) * (dg_frag%phi_frag(lx, ly, lz+3, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-3, jo, i_local)) + &
              nabt(4,3) * (dg_frag%phi_frag(lx, ly, lz+4, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-4, jo, i_local))
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine apply_gradient_to_basis

  subroutine apply_momentum_blocks(dg_frag, ispin, scale_vec, x, y)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer,                intent(in) :: ispin
    real(8),                intent(in) :: scale_vec(3)
    complex(8),             intent(in) :: x(:,:)
    complex(8),             intent(inout) :: y(:,:)

    integer :: iblk, idir, ib, jb, row_idx, col_idx, nstate
    real(8) :: scale

    if (.not. allocated(dg_frag%momentum_blocks)) return
    nstate = min(size(x, 2), size(y, 2))

    do iblk = 1, dg_frag%n_momentum_blocks
      if (.not. allocated(dg_frag%momentum_blocks(iblk)%val)) cycle
      do idir = 1, 3
        scale = scale_vec(idir)
        if (abs(scale) < 1.0d-30) cycle
        do jb = 1, dg_frag%n_basis(dg_frag%momentum_blocks(iblk)%ifrag_col, ispin)
          col_idx = dg_frag%index_basis(jb, dg_frag%momentum_blocks(iblk)%ifrag_col, ispin)
          if (col_idx < 1 .or. col_idx > size(x, 1)) cycle
          do ib = 1, dg_frag%n_basis(dg_frag%momentum_blocks(iblk)%ifrag_row, ispin)
            row_idx = dg_frag%index_basis(ib, dg_frag%momentum_blocks(iblk)%ifrag_row, ispin)
            if (row_idx < 1 .or. row_idx > size(y, 1)) cycle
            y(row_idx, 1:nstate) = y(row_idx, 1:nstate) + &
              scale * dg_frag%momentum_blocks(iblk)%val(idir, ib, jb, ispin) * x(col_idx, 1:nstate)
          end do
        end do
      end do
    end do

  end subroutine apply_momentum_blocks

  subroutine calculate_microscopic_current_dg(dg_frag, system, mg, stencil, curr)
    ! NOTE: DG microscopic current here intentionally excludes non-local PP contribution.
    use structures
    use salmon_global, only: nelec
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_vector),         intent(inout) :: curr

    integer :: ifrag, i_local, ispin, io, istate_frag, jstate_frag
    integer :: ig_i, ig_j
    integer :: ix, iy, iz, ixg, iyg, izg
    integer :: nxyz(3), ixyz0(3)
    integer :: nocc_per_spin
    real(8) :: occ_factor
    complex(8) :: coef_pair
    real(8) :: phi_i
    real(8), allocatable :: grad_phi(:,:,:,:)
    real(8), allocatable :: curr_local(:,:,:,:), curr_sum(:,:,:,:)
    real(8), allocatable :: w_local(:,:,:), w_sum(:,:,:)

    curr%v = 0.0d0
    if (.not. allocated(dg_frag%phi_frag)) return

    nocc_per_spin = min(dg_frag%nstate_tot, int(nelec / 2.0d0 + 1.0d-12))
    occ_factor = 2.0d0 / real(system%nspin, 8)

    allocate(curr_local(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(curr_sum(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(w_local(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(w_sum(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    curr_local = 0.0d0
    curr_sum = 0.0d0
    w_local = 0.0d0
    w_sum = 0.0d0

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
      ixyz0(1:3) = dg_frag%ixyz_frag(1:3, ifrag)

      do iz = 1, nxyz(3)
        do iy = 1, nxyz(2)
          do ix = 1, nxyz(1)
            ixg = modulo(ixyz0(1) + ix - 2, mg%num(1)) + 1
            iyg = modulo(ixyz0(2) + iy - 2, mg%num(2)) + 1
            izg = modulo(ixyz0(3) + iz - 2, mg%num(3)) + 1
            if (ixg < mg%is(1) .or. ixg > mg%ie(1)) cycle
            if (iyg < mg%is(2) .or. iyg > mg%ie(2)) cycle
            if (izg < mg%is(3) .or. izg > mg%ie(3)) cycle
            w_local(ixg, iyg, izg) = w_local(ixg, iyg, izg) + 1.0d0
          end do
        end do
      end do

      do ispin = 1, system%nspin
        do jstate_frag = 1, dg_frag%n_basis(ifrag, ispin)
          ig_j = dg_frag%index_basis(jstate_frag, ifrag, ispin)
          if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle

          allocate(grad_phi(1:nxyz(1), 1:nxyz(2), 1:nxyz(3), 3))
          call apply_gradient_to_basis(dg_frag, i_local, jstate_frag, mg, stencil, grad_phi)

          do istate_frag = 1, dg_frag%n_basis(ifrag, ispin)
            ig_i = dg_frag%index_basis(istate_frag, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle

            coef_pair = (0.0d0, 0.0d0)
            do io = 1, nocc_per_spin
              coef_pair = coef_pair + occ_factor * conjg(dg_frag%coef(ig_i, io, ispin)) * dg_frag%coef(ig_j, io, ispin)
            end do
            if (abs(aimag(coef_pair)) < 1.0d-18) cycle

            do iz = 1, nxyz(3)
              do iy = 1, nxyz(2)
                do ix = 1, nxyz(1)
                  ixg = modulo(ixyz0(1) + ix - 2, mg%num(1)) + 1
                  iyg = modulo(ixyz0(2) + iy - 2, mg%num(2)) + 1
                  izg = modulo(ixyz0(3) + iz - 2, mg%num(3)) + 1
                  if (ixg < mg%is(1) .or. ixg > mg%ie(1)) cycle
                  if (iyg < mg%is(2) .or. iyg > mg%ie(2)) cycle
                  if (izg < mg%is(3) .or. izg > mg%ie(3)) cycle

                  phi_i = dg_frag%phi_frag(ix, iy, iz, istate_frag, i_local)
                  curr_local(1, ixg, iyg, izg) = curr_local(1, ixg, iyg, izg) + aimag(coef_pair) * phi_i * grad_phi(ix, iy, iz, 1)
                  curr_local(2, ixg, iyg, izg) = curr_local(2, ixg, iyg, izg) + aimag(coef_pair) * phi_i * grad_phi(ix, iy, iz, 2)
                  curr_local(3, ixg, iyg, izg) = curr_local(3, ixg, iyg, izg) + aimag(coef_pair) * phi_i * grad_phi(ix, iy, iz, 3)
                end do
              end do
            end do
          end do

          deallocate(grad_phi)
        end do
      end do
    end do

    call comm_summation(curr_local, curr_sum, size(curr_local), dg_frag%icomm)
    call comm_summation(w_local, w_sum, size(w_local), dg_frag%icomm)

    do izg = mg%is(3), mg%ie(3)
      do iyg = mg%is(2), mg%ie(2)
        do ixg = mg%is(1), mg%ie(1)
          if (w_sum(ixg, iyg, izg) > 1.0d-12) then
            curr%v(:, ixg, iyg, izg) = curr_sum(:, ixg, iyg, izg) / w_sum(ixg, iyg, izg)
          else
            curr%v(:, ixg, iyg, izg) = curr_sum(:, ixg, iyg, izg)
          end if
        end do
      end do
    end do

    deallocate(curr_local, curr_sum, w_local, w_sum)

  end subroutine calculate_microscopic_current_dg

  !=======================================================================
  ! Build spatially resolved A·∇ and A^2/2 matrices from Ac_micro(r)
  !=======================================================================
  subroutine build_spatial_A_coupling_matrices(dg_frag, system, mg, stencil, ispin, Ap_mat, A2_mat)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    integer,                intent(in)    :: ispin
    real(8),                intent(out)   :: Ap_mat(dg_frag%n_mat_max, dg_frag%n_mat_max)
    real(8),                intent(out)   :: A2_mat(dg_frag%n_mat_max, dg_frag%n_mat_max)

    integer :: ifrag, i_local, io, jo
    integer :: ig_i, ig_j
    integer :: ix, iy, iz, ixg, iyg, izg
    integer :: nxyz(3), ixyz0(3)
    real(8) :: phi_i, A2val, Ap_int, A2_int
    real(8), allocatable :: grad_phi(:,:,:,:)
    real(8), allocatable :: Ap_spin(:,:,:), A2_spin(:,:,:)
    type(matrix_block_info), allocatable :: Ap_blocks(:), A2_blocks(:)
    integer, allocatable :: A_block_map(:,:)

    Ap_mat = 0.0d0
    A2_mat = 0.0d0
    if (.not. allocated(system%Ac_micro%v)) return

    allocate(Ap_spin(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    allocate(A2_spin(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    Ap_spin(:, :, :) = 0.0d0
    A2_spin(:, :, :) = 0.0d0

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
      ixyz0(1:3) = dg_frag%ixyz_frag(1:3, ifrag)

      do jo = 1, dg_frag%n_basis(ifrag, ispin)
        ig_j = dg_frag%index_basis(jo, ifrag, ispin)
        if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle

        allocate(grad_phi(1:nxyz(1), 1:nxyz(2), 1:nxyz(3), 3))
        call apply_gradient_to_basis(dg_frag, i_local, jo, mg, stencil, grad_phi)

        !$omp parallel do private(io, ig_i, Ap_int, A2_int, iz, iy, ix, ixg, iyg, izg, phi_i, A2val) schedule(static)
        do io = 1, dg_frag%n_basis(ifrag, ispin)
          ig_i = dg_frag%index_basis(io, ifrag, ispin)
          if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle

          Ap_int = 0.0d0
          A2_int = 0.0d0
          do iz = 1, nxyz(3)
            do iy = 1, nxyz(2)
              do ix = 1, nxyz(1)
                ixg = modulo(ixyz0(1) + ix - 2, mg%num(1)) + 1
                iyg = modulo(ixyz0(2) + iy - 2, mg%num(2)) + 1
                izg = modulo(ixyz0(3) + iz - 2, mg%num(3)) + 1
                if (ixg < mg%is(1) .or. ixg > mg%ie(1)) cycle
                if (iyg < mg%is(2) .or. iyg > mg%ie(2)) cycle
                if (izg < mg%is(3) .or. izg > mg%ie(3)) cycle

                phi_i = dg_frag%phi_frag(ix, iy, iz, io, i_local)
                Ap_int = Ap_int + phi_i * ( &
                         system%Ac_micro%v(1, ixg, iyg, izg) * grad_phi(ix, iy, iz, 1) + &
                         system%Ac_micro%v(2, ixg, iyg, izg) * grad_phi(ix, iy, iz, 2) + &
                         system%Ac_micro%v(3, ixg, iyg, izg) * grad_phi(ix, iy, iz, 3) ) * system%hvol

                A2val = system%Ac_micro%v(1, ixg, iyg, izg)**2 + &
                        system%Ac_micro%v(2, ixg, iyg, izg)**2 + &
                        system%Ac_micro%v(3, ixg, iyg, izg)**2
                A2_int = A2_int + 0.5d0 * phi_i * A2val * dg_frag%phi_frag(ix, iy, iz, jo, i_local) * system%hvol
              end do
            end do
          end do

          Ap_spin(ig_i, ig_j, ispin) = Ap_spin(ig_i, ig_j, ispin) + Ap_int
          A2_spin(ig_i, ig_j, ispin) = A2_spin(ig_i, ig_j, ispin) + A2_int
        end do
        !$omp end parallel do

        deallocate(grad_phi)
      end do
    end do

    call init_matrix_blocks_runtime(dg_frag, Ap_blocks, A_block_map)
    call init_matrix_blocks_runtime(dg_frag, A2_blocks, A_block_map)
    call sync_dense_matrix_to_blocks_runtime(dg_frag, Ap_spin, Ap_blocks, A_block_map)
    call sync_dense_matrix_to_blocks_runtime(dg_frag, A2_spin, A2_blocks, A_block_map)
    call reduce_matrix_blocks_runtime(dg_frag, Ap_blocks, "spatial-Ap", dg_frag%icomm)
    call reduce_matrix_blocks_runtime(dg_frag, A2_blocks, "spatial-A2", dg_frag%icomm)
    call sync_blocks_to_dense_matrix_runtime(dg_frag, Ap_blocks, A_block_map, Ap_spin)
    call sync_blocks_to_dense_matrix_runtime(dg_frag, A2_blocks, A_block_map, A2_spin)

    Ap_mat(:, :) = Ap_spin(:, :, ispin)
    A2_mat(:, :) = A2_spin(:, :, ispin)

    if (allocated(Ap_blocks)) then
      do io = 1, size(Ap_blocks)
        if (allocated(Ap_blocks(io)%val)) deallocate(Ap_blocks(io)%val)
      end do
      deallocate(Ap_blocks)
    end if
    if (allocated(A2_blocks)) then
      do io = 1, size(A2_blocks)
        if (allocated(A2_blocks(io)%val)) deallocate(A2_blocks(io)%val)
      end do
      deallocate(A2_blocks)
    end if
    if (allocated(A_block_map)) deallocate(A_block_map)
    deallocate(Ap_spin, A2_spin)

  end subroutine build_spatial_A_coupling_matrices

end module rt_dg_fragment_ops
