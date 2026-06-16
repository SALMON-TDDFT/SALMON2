  subroutine calculate_time_derivative(dg_frag, system, mg, ppg, Ac_tot, dcoef_dt, state_start, state_end)
    use structures
    use rt_dg_fragment_types, only: s_dg_fragment_rt, matrix_block_info, complex_matrix_block_info
    use rt_dg_fragment_ops, only: rebuild_local_h_block_ids, fetch_remote_coef_rows, ensure_nonlocal_pp_matrix_A
    use misc_routines, only: get_wtime
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    real(8),                intent(in)    :: Ac_tot(3)
    complex(8),             intent(out)   :: dcoef_dt(:,:,:)
    integer, optional,      intent(in)    :: state_start, state_end

    integer :: ispin
    integer :: irow, local_idx, local_col
    integer :: n_frag
    integer :: state_first, state_last, state0, state_s, state_e
    integer :: state_out_s, state_out_e, nstate_blk, nstate_work
    integer, parameter :: state_work_target_mb = 512
    integer, parameter :: state_work_vectors = 4
    integer(8) :: target_bytes, bytes_per_state
    real(8) :: A_squared
    complex(8), parameter :: zi = (0.0d0, 1.0d0)
    logical :: has_nonlocal, use_nonlocal_blocks
    logical :: has_so_nonlocal, output_is_block, output_is_local_rows
    logical :: has_overlap_operator
    logical :: trace_derivative
    logical, save :: derivative_timing_initialized = .false.
    logical, save :: enable_derivative_timing = .false.
    logical, save :: derivative_real_blas_initialized = .false.
    logical, save :: derivative_block_diag_h_initialized = .false.
    logical, save :: enable_block_diag_h = .false.
    integer, save :: derivative_timing_call_count = 0
    integer :: trace_call_id, env_status
    character(16) :: env_derivative_timing
    character(16) :: env_block_diag_h
    character(32) :: env_derivative_real_blas_min_ops
    integer(8), save :: derivative_real_blas_min_ops = 131072_8
    real(8) :: t0, t1
    real(8) :: time_setup, time_fetch, time_h_apply, time_nl_apply, time_m_apply
    real(8) :: time_finalize, time_scatter, time_total
    complex(8), allocatable, save :: coef_work(:,:), rhs_work(:,:)
    complex(8), allocatable, save :: h_work(:,:), m_work(:,:)
    complex(8), allocatable, save :: cmat_pack(:,:), cx_pack(:,:), cy_pack(:,:)
    real(8), allocatable, save :: rmat_pack(:,:), rx_pack(:,:), ix_pack(:,:), ry_pack(:,:), iy_pack(:,:)
    integer, allocatable, save :: needed_row_ids(:), needed_row_pos(:)
    integer, allocatable, save :: row_lid_work(:), row_pos_work(:)
    integer, allocatable, save :: col_lid_work(:), col_pos_work(:)
    integer, allocatable, save :: frag_map_lid(:,:), frag_map_pos(:,:), frag_map_count(:)
    logical, allocatable, save :: output_row_needed(:)
    logical, allocatable, save :: overlap_identity_verified(:)
    logical, save :: compact_row_cache_valid = .false.
    logical, save :: compact_frag_map_cache_valid = .false.
    integer, save :: compact_row_cache_generation = 0
    integer, save :: compact_frag_map_generation = -1
    integer, save :: compact_row_cache_ispin = -1
    integer, save :: compact_row_cache_n_frag = -1
    integer, save :: compact_row_cache_index_dim1 = -1
    integer, save :: compact_row_cache_n_basis_dim1 = -1
    integer, save :: compact_row_cache_n_basis_dim2 = -1
    integer, save :: compact_row_cache_h_local_count = -1
    integer, save :: compact_row_cache_h_nl_count = -1
    integer, save :: compact_row_cache_momentum_count = -1
    integer, save :: compact_row_cache_local_row_count = -1
    integer, allocatable, save :: compact_row_cache_h_local_ids(:)
    integer, allocatable, save :: compact_row_cache_h_nl_ids(:)
    integer, allocatable, save :: compact_row_cache_local_rows(:)

    if (.not. derivative_timing_initialized) then
      env_derivative_timing = ''
      call get_environment_variable('SALMON_DG_DERIV_TIMING', env_derivative_timing, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(env_derivative_timing)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enable_derivative_timing = .true.
        end select
      end if
      derivative_timing_initialized = .true.
    end if
    if (.not. derivative_real_blas_initialized) then
      env_derivative_real_blas_min_ops = ''
      call get_environment_variable('SALMON_DG_DERIV_REAL_BLAS_MIN_OPS', env_derivative_real_blas_min_ops, status=env_status)
      if (env_status == 0 .and. len_trim(env_derivative_real_blas_min_ops) > 0) then
        read(env_derivative_real_blas_min_ops, *, iostat=env_status) derivative_real_blas_min_ops
        if (env_status /= 0) derivative_real_blas_min_ops = 131072_8
      end if
      derivative_real_blas_initialized = .true.
    end if
    if (.not. derivative_block_diag_h_initialized) then
      env_block_diag_h = ''
      call get_environment_variable('SALMON_DG_BLOCK_DIAG_H', env_block_diag_h, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(env_block_diag_h)))
        case('0','n','N','no','NO','false','FALSE','off','OFF')
          enable_block_diag_h = .false.
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enable_block_diag_h = .true.
        end select
      end if
      derivative_block_diag_h_initialized = .true.
    end if
    trace_derivative = .false.
    trace_call_id = 0
    time_setup = 0.0d0
    time_fetch = 0.0d0
    time_h_apply = 0.0d0
    time_nl_apply = 0.0d0
    time_m_apply = 0.0d0
    time_finalize = 0.0d0
    time_scatter = 0.0d0
    time_total = 0.0d0
    if (enable_derivative_timing .and. dg_frag%id == 0 .and. derivative_timing_call_count < 12) then
      derivative_timing_call_count = derivative_timing_call_count + 1
      trace_call_id = derivative_timing_call_count
      trace_derivative = .true.
      t0 = get_wtime()
    end if

    dcoef_dt = (0.0d0, 0.0d0)

    if (dg_frag%use_plane_wave_basis .or. allocated(dg_frag%coef_pw)) then
      stop "DG derivative supports the pure fragment block-sparse route only"
    end if
    if (.not. allocated(dg_frag%H_mat_blocks)) then
      stop "DG derivative requires block-sparse H blocks"
    end if
    if (.not. dg_frag%coef_state_block_mode .and. .not. allocated(dg_frag%momentum_blocks)) then
      stop "DG derivative requires block-sparse DG velocity blocks"
    end if
    if (.not. dg_frag%coef_state_block_mode .and. dg_frag%dc_lcfo_seed_basis_cleaned .and. &
        .not. dg_frag%momentum_blocks_include_dg_flux) then
      stop "DG derivative requires covariant Flux contribution in DG velocity blocks"
    end if

    call rebuild_local_h_block_ids(dg_frag)
    if (.not. dg_frag%H_local_block_ids_valid .or. .not. allocated(dg_frag%H_local_block_ids)) then
      write(*,'(1x,a,i0)') '[FATAL] DG derivative requires row-owner-local H block ids: rank=', dg_frag%id
      stop 'DG-Fragment RT: missing row-owner-local H block ids'
    end if

    n_frag = dg_frag%n_mat_max
    if (n_frag <= 0) return
    if (.not. dg_frag%coef_state_block_mode .and. dg_frag%nstate_tot > size(dg_frag%coef, 2)) then
      stop "DG derivative invalid fragment coefficient state count"
    end if

    output_is_local_rows = allocated(dg_frag%coef_owner) .and. allocated(dg_frag%coef_global_to_local) .and. &
                           size(dcoef_dt, 1) == size(dg_frag%coef, 1) .and. size(dcoef_dt, 1) < n_frag
    if (.not. output_is_local_rows .and. size(dcoef_dt, 1) < n_frag) then
      stop "DG derivative output row count is neither global nor local-owned"
    end if

    output_is_block = present(state_start) .or. present(state_end)
    if (output_is_block) then
      if (.not. (present(state_start) .and. present(state_end))) then
        stop "DG derivative state_start/state_end must be provided together"
      end if
      state_first = max(1, state_start)
      state_last = min(dg_frag%nstate_tot, state_end)
      if (state_first > state_last) return
      if (size(dcoef_dt, 2) < state_last - state_first + 1) then
        stop "DG derivative output block has too few state columns"
      end if
    else
      state_first = 1
      state_last = dg_frag%nstate_tot
      if (size(dcoef_dt, 2) < dg_frag%nstate_tot) then
        stop "DG derivative output has too few state columns"
      end if
    end if

    A_squared = Ac_tot(1)**2 + Ac_tot(2)**2 + Ac_tot(3)**2
    if (A_squared /= A_squared) stop "NaN in DG derivative vector potential"

    has_nonlocal = (ppg%Nlma > 0 .and. allocated(ppg%uV))
    has_so_nonlocal = (allocated(ppg%uv_so) .and. allocated(dg_frag%phi_frag_c) .and. dg_frag%nspin == 2)
    if (has_so_nonlocal) then
      stop 'DG-Fragment RT: SO nonlocal projector route is disabled in compact RT derivative'
    end if
    if (has_nonlocal .and. .not. dg_frag%coef_state_block_mode) then
      call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_tot, .false.)
    end if
    use_nonlocal_blocks = has_nonlocal .and. allocated(dg_frag%H_nl_blocks) .and. &
                          allocated(dg_frag%H_nl_local_block_ids)

    has_overlap_operator = allocated(dg_frag%S_mat_prop_blocks) .or. allocated(dg_frag%S_mat_prop_c) .or. &
                           allocated(dg_frag%S_mat_prop) .or. allocated(dg_frag%S_mat_c) .or. allocated(dg_frag%S_mat)
    if (has_overlap_operator) then
      if (.not. (dg_frag%dc_lcfo_seed_basis_cleaned .and. .not. dg_frag%identity_seed_coefficients)) then
        stop 'DG-Fragment RT: compact derivative requires a DC-exported S-orthonormal seed basis'
      end if
    end if

    call build_needed_coefficient_rows(1)
    target_bytes = int(state_work_target_mb, kind=8) * 1024_8 * 1024_8
    bytes_per_state = 16_8 * int(max(1, size(needed_row_ids)), kind=8) * int(state_work_vectors, kind=8)
    nstate_work = max(1, min(state_last - state_first + 1, int(max(1_8, target_bytes / bytes_per_state))))
    if (trace_derivative) then
      t1 = get_wtime()
      time_setup = time_setup + (t1 - t0)
    end if

    do ispin = 1, dg_frag%nspin
      if (trace_derivative) t0 = get_wtime()
      call build_needed_coefficient_rows(ispin)
      call build_fragment_compact_maps(ispin)
      if (has_overlap_operator) call verify_local_overlap_identity(ispin)
      call ensure_work_arrays(size(needed_row_ids), nstate_work)
      if (trace_derivative) then
        t1 = get_wtime()
        time_setup = time_setup + (t1 - t0)
      end if

      do state0 = state_first, state_last, nstate_work
        nstate_blk = min(nstate_work, state_last - state0 + 1)
        state_s = state0
        state_e = state0 + nstate_blk - 1
        state_out_s = merge(state_s - state_first + 1, state_s, output_is_block)
        state_out_e = merge(state_e - state_first + 1, state_e, output_is_block)

        if (trace_derivative) t0 = get_wtime()
        call fetch_remote_coef_rows(dg_frag, ispin, needed_row_ids, coef_work(:, 1:nstate_blk), state_s, state_e)
        if (trace_derivative) then
          t1 = get_wtime()
          time_fetch = time_fetch + (t1 - t0)
        end if

        h_work(:, 1:nstate_blk) = (0.0d0, 0.0d0)
        if (trace_derivative) t0 = get_wtime()
        call apply_real_blocks_compact(dg_frag%H_mat_blocks, dg_frag%H_local_block_ids, ispin, &
                                       coef_work(:, 1:nstate_blk), h_work(:, 1:nstate_blk))
        if (trace_derivative) then
          t1 = get_wtime()
          time_h_apply = time_h_apply + (t1 - t0)
        end if
        if (use_nonlocal_blocks) then
          if (trace_derivative) t0 = get_wtime()
          call apply_complex_blocks_compact(dg_frag%H_nl_blocks, dg_frag%H_nl_local_block_ids, ispin, &
                                            coef_work(:, 1:nstate_blk), h_work(:, 1:nstate_blk))
          if (trace_derivative) then
            t1 = get_wtime()
            time_nl_apply = time_nl_apply + (t1 - t0)
          end if
        end if
        m_work(:, 1:nstate_blk) = (0.0d0, 0.0d0)
        ! momentum_blocks store the real gradient/Flux-gradient matrix G.
        ! The velocity-gauge Hamiltonian contribution is -i A.G, so the
        ! time derivative contains -A.G after applying -i H.
        if (allocated(dg_frag%momentum_blocks)) then
          if (trace_derivative) t0 = get_wtime()
          call apply_momentum_blocks_compact(ispin, Ac_tot, coef_work(:, 1:nstate_blk), m_work(:, 1:nstate_blk))
          if (trace_derivative) then
            t1 = get_wtime()
            time_m_apply = time_m_apply + (t1 - t0)
          end if
        end if

        if (trace_derivative) t0 = get_wtime()
        call add_diamagnetic_term(nstate_blk)
        h_work(:, 1:nstate_blk) = -zi * h_work(:, 1:nstate_blk)
        h_work(:, 1:nstate_blk) = h_work(:, 1:nstate_blk) - m_work(:, 1:nstate_blk)

        rhs_work(:, 1:nstate_blk) = h_work(:, 1:nstate_blk)
        if (trace_derivative) then
          t1 = get_wtime()
          time_finalize = time_finalize + (t1 - t0)
          t0 = get_wtime()
        end if
        call scatter_owned_derivative(ispin, nstate_blk, state_out_s, state_out_e)
        if (trace_derivative) then
          t1 = get_wtime()
          time_scatter = time_scatter + (t1 - t0)
        end if
      end do
    end do
    if (trace_derivative) then
      time_total = time_setup + time_fetch + time_h_apply + time_nl_apply + time_m_apply + time_finalize + time_scatter
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,8(a,1pe12.4))') &
        '[DG-DERIV] call=', trace_call_id, ' states=', state_last - state_first + 1, &
        ' nstate_work=', nstate_work, ' nneeded=', size(needed_row_ids), &
        ' setup=', time_setup, ' fetch=', time_fetch, ' h=', time_h_apply, &
        ' nl=', time_nl_apply, ' mom=', time_m_apply, ' finalize=', time_finalize, &
        ' scatter=', time_scatter, ' total=', time_total
      flush(6)
    end if

  contains

    subroutine ensure_work_arrays(nrow, ncol)
      integer, intent(in) :: nrow, ncol

      if (.not. allocated(coef_work)) then
        allocate(coef_work(nrow, ncol), rhs_work(nrow, ncol), h_work(nrow, ncol), m_work(nrow, ncol))
      else if (size(coef_work, 1) /= nrow .or. size(coef_work, 2) /= ncol) then
        deallocate(coef_work, rhs_work, h_work, m_work)
        allocate(coef_work(nrow, ncol), rhs_work(nrow, ncol), h_work(nrow, ncol), m_work(nrow, ncol))
      end if
    end subroutine ensure_work_arrays

    subroutine ensure_index_work_arrays(nmax)
      integer, intent(in) :: nmax

      if (.not. allocated(row_lid_work)) then
        allocate(row_lid_work(nmax), row_pos_work(nmax), col_lid_work(nmax), col_pos_work(nmax))
      else if (size(row_lid_work) < nmax) then
        deallocate(row_lid_work, row_pos_work, col_lid_work, col_pos_work)
        allocate(row_lid_work(nmax), row_pos_work(nmax), col_lid_work(nmax), col_pos_work(nmax))
      end if
    end subroutine ensure_index_work_arrays

    subroutine ensure_fragment_map_arrays(nmax)
      integer, intent(in) :: nmax

      if (.not. allocated(frag_map_lid)) then
        allocate(frag_map_lid(nmax, max(1, dg_frag%n_frag)), &
                 frag_map_pos(nmax, max(1, dg_frag%n_frag)), &
                 frag_map_count(max(1, dg_frag%n_frag)))
      else if (size(frag_map_lid, 1) < nmax .or. size(frag_map_lid, 2) /= max(1, dg_frag%n_frag)) then
        deallocate(frag_map_lid, frag_map_pos, frag_map_count)
        allocate(frag_map_lid(nmax, max(1, dg_frag%n_frag)), &
                 frag_map_pos(nmax, max(1, dg_frag%n_frag)), &
                 frag_map_count(max(1, dg_frag%n_frag)))
      end if
    end subroutine ensure_fragment_map_arrays

    subroutine ensure_dgemm_pack_arrays(nrow, ncol, nstate)
      integer, intent(in) :: nrow, ncol, nstate

      if (.not. allocated(rmat_pack)) then
        allocate(rmat_pack(nrow, ncol), rx_pack(ncol, nstate), ix_pack(ncol, nstate), &
                 ry_pack(nrow, nstate), iy_pack(nrow, nstate))
      else if (size(rmat_pack, 1) < nrow .or. size(rmat_pack, 2) < ncol .or. &
               size(rx_pack, 1) < ncol .or. size(rx_pack, 2) < nstate .or. &
               size(ry_pack, 1) < nrow .or. size(ry_pack, 2) < nstate) then
        deallocate(rmat_pack, rx_pack, ix_pack, ry_pack, iy_pack)
        allocate(rmat_pack(nrow, ncol), rx_pack(ncol, nstate), ix_pack(ncol, nstate), &
                 ry_pack(nrow, nstate), iy_pack(nrow, nstate))
      end if
    end subroutine ensure_dgemm_pack_arrays

    subroutine ensure_zgemm_pack_arrays(nrow, ncol, nstate)
      integer, intent(in) :: nrow, ncol, nstate

      if (.not. allocated(cmat_pack)) then
        allocate(cmat_pack(nrow, ncol), cx_pack(ncol, nstate), cy_pack(nrow, nstate))
      else if (size(cmat_pack, 1) < nrow .or. size(cmat_pack, 2) < ncol .or. &
               size(cx_pack, 1) < ncol .or. size(cx_pack, 2) < nstate .or. &
               size(cy_pack, 1) < nrow .or. size(cy_pack, 2) < nstate) then
        deallocate(cmat_pack, cx_pack, cy_pack)
        allocate(cmat_pack(nrow, ncol), cx_pack(ncol, nstate), cy_pack(nrow, nstate))
      end if
    end subroutine ensure_zgemm_pack_arrays

    subroutine build_fragment_compact_maps(ispin_current)
      integer, intent(in) :: ispin_current
      integer :: ifrag, ib, nbf, max_nbf, global_idx, compact_pos

      if (compact_frag_map_cache_valid .and. compact_frag_map_generation == compact_row_cache_generation) return
      max_nbf = max(1, size(dg_frag%index_basis, 1))
      call ensure_fragment_map_arrays(max_nbf)
      frag_map_count(:) = 0
      frag_map_lid(:, :) = 0
      frag_map_pos(:, :) = 0
      do ifrag = 1, dg_frag%n_frag
        nbf = min(dg_frag%n_basis(ifrag, ispin_current), size(dg_frag%index_basis, 1))
        do ib = 1, nbf
          global_idx = dg_frag%index_basis(ib, ifrag, ispin_current)
          if (global_idx < 1 .or. global_idx > size(needed_row_pos)) cycle
          compact_pos = needed_row_pos(global_idx)
          if (compact_pos <= 0) cycle
          frag_map_count(ifrag) = frag_map_count(ifrag) + 1
          frag_map_lid(frag_map_count(ifrag), ifrag) = ib
          frag_map_pos(frag_map_count(ifrag), ifrag) = compact_pos
        end do
      end do
      compact_frag_map_generation = compact_row_cache_generation
      compact_frag_map_cache_valid = .true.
    end subroutine build_fragment_compact_maps

    subroutine build_compact_basis_map(ifrag, ispin_current, lid, pos, nvalid)
      integer, intent(in) :: ifrag, ispin_current
      integer, intent(out) :: lid(:), pos(:), nvalid
      integer :: ib, nbf, global_idx, compact_pos

      nvalid = 0
      if (ifrag < 1 .or. ifrag > dg_frag%n_frag) return
      nbf = min(dg_frag%n_basis(ifrag, ispin_current), size(dg_frag%index_basis, 1))
      nbf = min(nbf, size(lid))
      do ib = 1, nbf
        global_idx = dg_frag%index_basis(ib, ifrag, ispin_current)
        if (global_idx < 1 .or. global_idx > size(needed_row_pos)) cycle
        compact_pos = needed_row_pos(global_idx)
        if (compact_pos <= 0) cycle
        nvalid = nvalid + 1
        lid(nvalid) = ib
        pos(nvalid) = compact_pos
      end do
    end subroutine build_compact_basis_map

    subroutine build_needed_coefficient_rows(ispin_current)
      integer, intent(in) :: ispin_current
      logical, allocatable :: row_needed(:)
      integer :: iblk_idx, iblk, ifrag_row, ifrag_col, global_idx, n_needed
      integer :: h_local_count, h_nl_count, local_row_count
      integer :: index_dim1, n_basis_dim1, n_basis_dim2

      h_local_count = -1
      h_nl_count = -1
      local_row_count = -1
      if (allocated(dg_frag%H_local_block_ids)) h_local_count = size(dg_frag%H_local_block_ids)
      if (allocated(dg_frag%H_nl_local_block_ids)) h_nl_count = size(dg_frag%H_nl_local_block_ids)
      if (allocated(dg_frag%H_local_rows)) local_row_count = size(dg_frag%H_local_rows)
      index_dim1 = -1
      n_basis_dim1 = -1
      n_basis_dim2 = -1
      if (allocated(dg_frag%index_basis)) index_dim1 = size(dg_frag%index_basis, 1)
      if (allocated(dg_frag%n_basis)) then
        n_basis_dim1 = size(dg_frag%n_basis, 1)
        n_basis_dim2 = size(dg_frag%n_basis, 2)
      end if

      if (compact_row_cache_valid .and. compact_row_cache_ispin == ispin_current .and. &
          compact_row_cache_n_frag == n_frag .and. &
          compact_row_cache_index_dim1 == index_dim1 .and. &
          compact_row_cache_n_basis_dim1 == n_basis_dim1 .and. &
          compact_row_cache_n_basis_dim2 == n_basis_dim2 .and. &
          compact_row_cache_h_local_count == h_local_count .and. &
          compact_row_cache_h_nl_count == h_nl_count .and. &
          compact_row_cache_momentum_count == dg_frag%n_momentum_blocks .and. &
          compact_row_cache_local_row_count == local_row_count .and. &
          allocated(needed_row_ids) .and. allocated(needed_row_pos) .and. allocated(output_row_needed)) then
        if (compact_row_id_cache_matches()) return
      end if

      if (allocated(needed_row_ids)) deallocate(needed_row_ids)
      if (.not. allocated(needed_row_pos)) then
        allocate(needed_row_pos(max(1, n_frag)))
      else if (size(needed_row_pos) /= max(1, n_frag)) then
        deallocate(needed_row_pos)
        allocate(needed_row_pos(max(1, n_frag)))
      end if
      if (.not. allocated(output_row_needed)) then
        allocate(output_row_needed(max(1, n_frag)))
      else if (size(output_row_needed) /= max(1, n_frag)) then
        deallocate(output_row_needed)
        allocate(output_row_needed(max(1, n_frag)))
      end if
      needed_row_pos(:) = 0
      output_row_needed(:) = .false.

      allocate(row_needed(max(1, n_frag)))
      row_needed(:) = .false.
      do iblk_idx = 1, size(dg_frag%H_local_block_ids)
        iblk = dg_frag%H_local_block_ids(iblk_idx)
        if (iblk < 1 .or. iblk > size(dg_frag%H_mat_blocks)) cycle
        ifrag_row = dg_frag%H_mat_blocks(iblk)%ifrag_row
        ifrag_col = dg_frag%H_mat_blocks(iblk)%ifrag_col
        call mark_fragment_basis_rows(row_needed, ifrag_row, ispin_current)
        call mark_output_fragment_rows(ifrag_row, ispin_current)
        call mark_fragment_basis_rows(row_needed, ifrag_col, ispin_current)
      end do
      if (allocated(dg_frag%H_nl_blocks) .and. allocated(dg_frag%H_nl_local_block_ids)) then
        do iblk_idx = 1, size(dg_frag%H_nl_local_block_ids)
          iblk = dg_frag%H_nl_local_block_ids(iblk_idx)
          if (iblk < 1 .or. iblk > size(dg_frag%H_nl_blocks)) cycle
          ifrag_row = dg_frag%H_nl_blocks(iblk)%ifrag_row
          ifrag_col = dg_frag%H_nl_blocks(iblk)%ifrag_col
          call mark_fragment_basis_rows(row_needed, ifrag_row, ispin_current)
          call mark_output_fragment_rows(ifrag_row, ispin_current)
          call mark_fragment_basis_rows(row_needed, ifrag_col, ispin_current)
        end do
      end if
      if (allocated(dg_frag%momentum_blocks)) then
        do iblk = 1, dg_frag%n_momentum_blocks
          if (iblk < 1 .or. iblk > size(dg_frag%momentum_blocks)) cycle
          ifrag_row = dg_frag%momentum_blocks(iblk)%ifrag_row
          ifrag_col = dg_frag%momentum_blocks(iblk)%ifrag_col
          if (.not. fragment_row_is_local(ifrag_row)) cycle
          call mark_fragment_basis_rows(row_needed, ifrag_row, ispin_current)
          call mark_output_fragment_rows(ifrag_row, ispin_current)
          call mark_fragment_basis_rows(row_needed, ifrag_col, ispin_current)
        end do
      end if

      n_needed = count(row_needed(1:n_frag))
      if (n_needed <= 0) stop 'DG-Fragment RT: derivative has no local compact rows'
      allocate(needed_row_ids(n_needed))
      n_needed = 0
      do global_idx = 1, n_frag
        if (.not. row_needed(global_idx)) cycle
        n_needed = n_needed + 1
        needed_row_ids(n_needed) = global_idx
        needed_row_pos(global_idx) = n_needed
      end do
      deallocate(row_needed)
      compact_row_cache_ispin = ispin_current
      compact_row_cache_n_frag = n_frag
      compact_row_cache_index_dim1 = index_dim1
      compact_row_cache_n_basis_dim1 = n_basis_dim1
      compact_row_cache_n_basis_dim2 = n_basis_dim2
      compact_row_cache_h_local_count = h_local_count
      compact_row_cache_h_nl_count = h_nl_count
      compact_row_cache_momentum_count = dg_frag%n_momentum_blocks
      compact_row_cache_local_row_count = local_row_count
      call save_compact_row_id_cache()
      compact_row_cache_generation = compact_row_cache_generation + 1
      compact_row_cache_valid = .true.
      compact_frag_map_cache_valid = .false.
    end subroutine build_needed_coefficient_rows

    logical function compact_row_id_cache_matches() result(matches)
      matches = .false.
      if (allocated(dg_frag%H_local_block_ids)) then
        if (.not. allocated(compact_row_cache_h_local_ids)) return
        if (size(compact_row_cache_h_local_ids) /= size(dg_frag%H_local_block_ids)) return
        if (any(compact_row_cache_h_local_ids /= dg_frag%H_local_block_ids)) return
      else if (allocated(compact_row_cache_h_local_ids)) then
        return
      end if
      if (allocated(dg_frag%H_nl_local_block_ids)) then
        if (.not. allocated(compact_row_cache_h_nl_ids)) return
        if (size(compact_row_cache_h_nl_ids) /= size(dg_frag%H_nl_local_block_ids)) return
        if (any(compact_row_cache_h_nl_ids /= dg_frag%H_nl_local_block_ids)) return
      else if (allocated(compact_row_cache_h_nl_ids)) then
        return
      end if
      if (allocated(dg_frag%H_local_rows)) then
        if (.not. allocated(compact_row_cache_local_rows)) return
        if (size(compact_row_cache_local_rows) /= size(dg_frag%H_local_rows)) return
        if (any(compact_row_cache_local_rows /= dg_frag%H_local_rows)) return
      else if (allocated(compact_row_cache_local_rows)) then
        return
      end if
      matches = .true.
    end function compact_row_id_cache_matches

    subroutine save_compact_row_id_cache()
      if (allocated(compact_row_cache_h_local_ids)) deallocate(compact_row_cache_h_local_ids)
      if (allocated(compact_row_cache_h_nl_ids)) deallocate(compact_row_cache_h_nl_ids)
      if (allocated(compact_row_cache_local_rows)) deallocate(compact_row_cache_local_rows)
      if (allocated(dg_frag%H_local_block_ids)) then
        allocate(compact_row_cache_h_local_ids(size(dg_frag%H_local_block_ids)))
        compact_row_cache_h_local_ids(:) = dg_frag%H_local_block_ids(:)
      end if
      if (allocated(dg_frag%H_nl_local_block_ids)) then
        allocate(compact_row_cache_h_nl_ids(size(dg_frag%H_nl_local_block_ids)))
        compact_row_cache_h_nl_ids(:) = dg_frag%H_nl_local_block_ids(:)
      end if
      if (allocated(dg_frag%H_local_rows)) then
        allocate(compact_row_cache_local_rows(size(dg_frag%H_local_rows)))
        compact_row_cache_local_rows(:) = dg_frag%H_local_rows(:)
      end if
    end subroutine save_compact_row_id_cache

    subroutine mark_fragment_basis_rows(row_needed, ifrag, ispin_current)
      logical, intent(inout) :: row_needed(:)
      integer, intent(in) :: ifrag, ispin_current
      integer :: ib, global_idx

      if (ifrag < 1 .or. ifrag > dg_frag%n_frag) return
      do ib = 1, min(dg_frag%n_basis(ifrag, ispin_current), size(dg_frag%index_basis, 1))
        global_idx = dg_frag%index_basis(ib, ifrag, ispin_current)
        if (global_idx < 1 .or. global_idx > size(row_needed)) cycle
        row_needed(global_idx) = .true.
      end do
    end subroutine mark_fragment_basis_rows

    subroutine mark_output_fragment_rows(ifrag, ispin_current)
      integer, intent(in) :: ifrag, ispin_current
      integer :: ib, global_idx

      if (ifrag < 1 .or. ifrag > dg_frag%n_frag) return
      do ib = 1, min(dg_frag%n_basis(ifrag, ispin_current), size(dg_frag%index_basis, 1))
        global_idx = dg_frag%index_basis(ib, ifrag, ispin_current)
        if (global_idx < 1 .or. global_idx > size(output_row_needed)) cycle
        output_row_needed(global_idx) = .true.
      end do
    end subroutine mark_output_fragment_rows

    subroutine verify_local_overlap_identity(ispin_current)
      integer, intent(in) :: ispin_current
      real(8), parameter :: s_identity_tol = 1.0d-10
      real(8) :: diag_dev, offdiag_max

      if (.not. allocated(overlap_identity_verified)) then
        allocate(overlap_identity_verified(max(1, dg_frag%nspin)))
        overlap_identity_verified(:) = .false.
      else if (size(overlap_identity_verified) < max(1, dg_frag%nspin)) then
        deallocate(overlap_identity_verified)
        allocate(overlap_identity_verified(max(1, dg_frag%nspin)))
        overlap_identity_verified(:) = .false.
      end if
      if (ispin_current <= size(overlap_identity_verified)) then
        if (overlap_identity_verified(ispin_current)) return
      end if

      diag_dev = 0.0d0
      offdiag_max = 0.0d0
      if (allocated(dg_frag%S_mat_prop_blocks)) then
        call scan_real_overlap_blocks(dg_frag%S_mat_prop_blocks, ispin_current, diag_dev, offdiag_max)
      else if (allocated(dg_frag%S_mat_blocks)) then
        call scan_real_overlap_blocks(dg_frag%S_mat_blocks, ispin_current, diag_dev, offdiag_max)
      else if (allocated(dg_frag%S_mat_prop_c)) then
        call scan_complex_overlap_dense(dg_frag%S_mat_prop_c(:, :, ispin_current), diag_dev, offdiag_max)
      else if (allocated(dg_frag%S_mat_c)) then
        call scan_complex_overlap_dense(dg_frag%S_mat_c(:, :, ispin_current), diag_dev, offdiag_max)
      else if (allocated(dg_frag%S_mat_prop)) then
        call scan_real_overlap_dense(dg_frag%S_mat_prop(:, :, ispin_current), diag_dev, offdiag_max)
      else if (allocated(dg_frag%S_mat)) then
        call scan_real_overlap_dense(dg_frag%S_mat(:, :, ispin_current), diag_dev, offdiag_max)
      end if
      if (diag_dev > s_identity_tol .or. offdiag_max > s_identity_tol) then
        write(*,'(1x,a,i0,a,1pe12.5,a,1pe12.5,a,1pe12.5)') &
          '[FATAL] compact DG derivative cannot skip non-identity S: ispin=', ispin_current, &
          ' diag_dev=', diag_dev, ' offdiag_max=', offdiag_max, ' tol=', s_identity_tol
        stop 'DG-Fragment RT: compact derivative requires identity propagation overlap'
      end if
      if (ispin_current <= size(overlap_identity_verified)) overlap_identity_verified(ispin_current) = .true.
    end subroutine verify_local_overlap_identity

    subroutine scan_real_overlap_blocks(blocks, ispin_block, diag_out, offdiag_out)
      type(matrix_block_info), intent(in) :: blocks(:)
      integer, intent(in) :: ispin_block
      real(8), intent(inout) :: diag_out, offdiag_out
      integer :: iblk, ii, jj, ifrag_row, ifrag_col, nrow, ncol
      integer :: row_gid, col_gid
      real(8) :: val_re

      do iblk = 1, size(blocks)
        ifrag_row = blocks(iblk)%ifrag_row
        ifrag_col = blocks(iblk)%ifrag_col
        if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) cycle
        if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
        nrow = min(dg_frag%n_basis(ifrag_row, ispin_block), size(blocks(iblk)%val, 1))
        ncol = min(dg_frag%n_basis(ifrag_col, ispin_block), size(blocks(iblk)%val, 2))
        do ii = 1, nrow
          row_gid = dg_frag%index_basis(ii, ifrag_row, ispin_block)
          if (.not. row_is_output(row_gid)) cycle
          do jj = 1, ncol
            col_gid = dg_frag%index_basis(jj, ifrag_col, ispin_block)
            val_re = blocks(iblk)%val(ii, jj, ispin_block)
            if (row_gid == col_gid) then
              diag_out = max(diag_out, abs(val_re - 1.0d0))
            else
              offdiag_out = max(offdiag_out, abs(val_re))
            end if
          end do
        end do
      end do
    end subroutine scan_real_overlap_blocks

    subroutine scan_real_overlap_dense(mat, diag_out, offdiag_out)
      real(8), intent(in) :: mat(:, :)
      real(8), intent(inout) :: diag_out, offdiag_out
      integer :: row_idx, col_idx

      do row_idx = 1, min(size(mat, 1), size(output_row_needed))
        if (.not. output_row_needed(row_idx)) cycle
        do col_idx = 1, size(mat, 2)
          if (row_idx == col_idx) then
            diag_out = max(diag_out, abs(mat(row_idx, col_idx) - 1.0d0))
          else
            offdiag_out = max(offdiag_out, abs(mat(row_idx, col_idx)))
          end if
        end do
      end do
    end subroutine scan_real_overlap_dense

    subroutine scan_complex_overlap_dense(mat, diag_out, offdiag_out)
      complex(8), intent(in) :: mat(:, :)
      real(8), intent(inout) :: diag_out, offdiag_out
      integer :: row_idx, col_idx

      do row_idx = 1, min(size(mat, 1), size(output_row_needed))
        if (.not. output_row_needed(row_idx)) cycle
        do col_idx = 1, size(mat, 2)
          if (row_idx == col_idx) then
            diag_out = max(diag_out, abs(mat(row_idx, col_idx) - cmplx(1.0d0, 0.0d0, kind=8)))
          else
            offdiag_out = max(offdiag_out, abs(mat(row_idx, col_idx)))
          end if
        end do
      end do
    end subroutine scan_complex_overlap_dense

    logical function fragment_row_is_local(ifrag) result(is_local)
      integer, intent(in) :: ifrag
      integer :: ib, ispin_current, global_idx

      if (ifrag < 1 .or. ifrag > dg_frag%n_frag) then
        is_local = .false.
      else if (allocated(dg_frag%H_local_rows)) then
        is_local = any(dg_frag%H_local_rows == ifrag)
      else if (allocated(dg_frag%coef_owner) .and. allocated(dg_frag%index_basis)) then
        is_local = .false.
        do ispin_current = 1, min(dg_frag%nspin, size(dg_frag%coef_owner, 2), size(dg_frag%index_basis, 3))
          do ib = 1, min(dg_frag%n_basis(ifrag, ispin_current), size(dg_frag%index_basis, 1))
            global_idx = dg_frag%index_basis(ib, ifrag, ispin_current)
            if (global_idx < 1 .or. global_idx > size(dg_frag%coef_owner, 1)) cycle
            if (dg_frag%coef_owner(global_idx, ispin_current) /= dg_frag%id) cycle
            is_local = .true.
            return
          end do
        end do
      else
        is_local = .false.
      end if
    end function fragment_row_is_local

    subroutine apply_real_blocks_compact(blocks, block_ids, ispin_current, x, y)
      type(matrix_block_info), intent(in) :: blocks(:)
      integer, intent(in) :: block_ids(:), ispin_current
      complex(8), intent(in) :: x(:, :)
      complex(8), intent(inout) :: y(:, :)
      integer :: iblk_idx, iblk

      do iblk_idx = 1, size(block_ids)
        iblk = block_ids(iblk_idx)
        if (iblk < 1 .or. iblk > size(blocks)) cycle
        call apply_one_real_block(blocks(iblk), ispin_current, x, y)
      end do
    end subroutine apply_real_blocks_compact

    subroutine apply_complex_blocks_compact(blocks, block_ids, ispin_current, x, y)
      type(complex_matrix_block_info), intent(in) :: blocks(:)
      integer, intent(in) :: block_ids(:), ispin_current
      complex(8), intent(in) :: x(:, :)
      complex(8), intent(inout) :: y(:, :)
      integer :: iblk_idx, iblk

      do iblk_idx = 1, size(block_ids)
        iblk = block_ids(iblk_idx)
        if (iblk < 1 .or. iblk > size(blocks)) cycle
        call apply_one_complex_block(blocks(iblk), ispin_current, x, y)
      end do
    end subroutine apply_complex_blocks_compact

    subroutine apply_one_real_block(block, ispin_current, x, y)
      type(matrix_block_info), intent(in) :: block
      integer, intent(in) :: ispin_current
      complex(8), intent(in) :: x(:, :)
      complex(8), intent(inout) :: y(:, :)
      integer :: iv, jv, ii, jj, istate, ifrag_row, ifrag_col, nrow, ncol
      integer :: row_pos, col_pos, nstate, nmax
      integer(8) :: block_ops
      real(8), parameter :: one = 1.0d0, zero = 0.0d0

      nstate = min(size(x, 2), size(y, 2))
      ifrag_row = block%ifrag_row
      ifrag_col = block%ifrag_col
      if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) return
      if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) return
      if (allocated(frag_map_count)) then
        nrow = frag_map_count(ifrag_row)
        ncol = frag_map_count(ifrag_col)
      else
        nmax = max(dg_frag%n_basis(ifrag_row, ispin_current), dg_frag%n_basis(ifrag_col, ispin_current))
        call ensure_index_work_arrays(max(1, nmax))
        call build_compact_basis_map(ifrag_row, ispin_current, row_lid_work, row_pos_work, nrow)
        call build_compact_basis_map(ifrag_col, ispin_current, col_lid_work, col_pos_work, ncol)
      end if
      if (nrow <= 0 .or. ncol <= 0) return
      block_ops = int(nrow, kind=8) * int(ncol, kind=8) * int(nstate, kind=8)
      if (nrow < 64 .or. ncol < 64 .or. block_ops < derivative_real_blas_min_ops) then
        call apply_one_real_block_scalar(block, ispin_current, x, y)
        return
      end if
      call ensure_dgemm_pack_arrays(nrow, ncol, nstate)
      do jv = 1, ncol
        if (allocated(frag_map_count)) then
          jj = frag_map_lid(jv, ifrag_col)
          col_pos = frag_map_pos(jv, ifrag_col)
        else
          jj = col_lid_work(jv)
          col_pos = col_pos_work(jv)
        end if
        do iv = 1, nrow
          if (allocated(frag_map_count)) then
            ii = frag_map_lid(iv, ifrag_row)
          else
            ii = row_lid_work(iv)
          end if
          rmat_pack(iv, jv) = block%val(ii, jj, ispin_current)
        end do
        do istate = 1, nstate
          rx_pack(jv, istate) = real(x(col_pos, istate), kind=8)
          ix_pack(jv, istate) = aimag(x(col_pos, istate))
        end do
      end do
      call dgemm('N', 'N', nrow, nstate, ncol, one, rmat_pack, size(rmat_pack, 1), &
                 rx_pack, size(rx_pack, 1), zero, ry_pack, size(ry_pack, 1))
      call dgemm('N', 'N', nrow, nstate, ncol, one, rmat_pack, size(rmat_pack, 1), &
                 ix_pack, size(ix_pack, 1), zero, iy_pack, size(iy_pack, 1))
      do iv = 1, nrow
        if (allocated(frag_map_count)) then
          row_pos = frag_map_pos(iv, ifrag_row)
        else
          row_pos = row_pos_work(iv)
        end if
        do istate = 1, nstate
          y(row_pos, istate) = y(row_pos, istate) + cmplx(ry_pack(iv, istate), iy_pack(iv, istate), kind=8)
        end do
      end do
    end subroutine apply_one_real_block

    subroutine apply_one_real_block_scalar(block, ispin_current, x, y)
      type(matrix_block_info), intent(in) :: block
      integer, intent(in) :: ispin_current
      complex(8), intent(in) :: x(:, :)
      complex(8), intent(inout) :: y(:, :)
      integer :: iv, jv, ii, jj, istate, ifrag_row, ifrag_col, nrow, ncol
      integer :: row_pos, col_pos, nstate, nmax
      real(8) :: val

      ifrag_row = block%ifrag_row
      ifrag_col = block%ifrag_col
      if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) return
      if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) return
      if (allocated(frag_map_count)) then
        nrow = frag_map_count(ifrag_row)
        ncol = frag_map_count(ifrag_col)
      else
        nmax = max(dg_frag%n_basis(ifrag_row, ispin_current), dg_frag%n_basis(ifrag_col, ispin_current))
        call ensure_index_work_arrays(max(1, nmax))
        call build_compact_basis_map(ifrag_row, ispin_current, row_lid_work, row_pos_work, nrow)
        call build_compact_basis_map(ifrag_col, ispin_current, col_lid_work, col_pos_work, ncol)
      end if
      if (nrow <= 0 .or. ncol <= 0) return
      nstate = min(size(x, 2), size(y, 2))
!$omp parallel do if(nrow*ncol*nstate >= 8192) private(iv,jv,ii,jj,row_pos,col_pos,val,istate) schedule(static)
      do iv = 1, nrow
        if (allocated(frag_map_count)) then
          ii = frag_map_lid(iv, ifrag_row)
          row_pos = frag_map_pos(iv, ifrag_row)
        else
          ii = row_lid_work(iv)
          row_pos = row_pos_work(iv)
        end if
        do jv = 1, ncol
          if (allocated(frag_map_count)) then
            jj = frag_map_lid(jv, ifrag_col)
            col_pos = frag_map_pos(jv, ifrag_col)
          else
            jj = col_lid_work(jv)
            col_pos = col_pos_work(jv)
          end if
          val = block%val(ii, jj, ispin_current)
          if (val == 0.0d0) cycle
          !$omp simd
          do istate = 1, nstate
            y(row_pos, istate) = y(row_pos, istate) + val * x(col_pos, istate)
          end do
        end do
      end do
!$omp end parallel do
    end subroutine apply_one_real_block_scalar

    subroutine apply_one_complex_block(block, ispin_current, x, y)
      type(complex_matrix_block_info), intent(in) :: block
      integer, intent(in) :: ispin_current
      complex(8), intent(in) :: x(:, :)
      complex(8), intent(inout) :: y(:, :)
      integer :: iv, jv, ii, jj, istate, ifrag_row, ifrag_col, nrow, ncol
      integer :: row_pos, col_pos, nstate, nmax
      integer(8) :: block_ops
      complex(8) :: val
      complex(8), parameter :: cone = (1.0d0, 0.0d0), czero = (0.0d0, 0.0d0)

      ifrag_row = block%ifrag_row
      ifrag_col = block%ifrag_col
      if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) return
      if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) return
      if (allocated(frag_map_count)) then
        nrow = frag_map_count(ifrag_row)
        ncol = frag_map_count(ifrag_col)
      else
        nmax = max(dg_frag%n_basis(ifrag_row, ispin_current), dg_frag%n_basis(ifrag_col, ispin_current))
        call ensure_index_work_arrays(max(1, nmax))
        call build_compact_basis_map(ifrag_row, ispin_current, row_lid_work, row_pos_work, nrow)
        call build_compact_basis_map(ifrag_col, ispin_current, col_lid_work, col_pos_work, ncol)
      end if
      if (nrow <= 0 .or. ncol <= 0) return
      nstate = min(size(x, 2), size(y, 2))
      block_ops = int(nrow, kind=8) * int(ncol, kind=8) * int(nstate, kind=8)
      if (nrow >= 64 .and. ncol >= 64 .and. block_ops >= derivative_real_blas_min_ops) then
        call ensure_zgemm_pack_arrays(nrow, ncol, nstate)
        do jv = 1, ncol
          if (allocated(frag_map_count)) then
            jj = frag_map_lid(jv, ifrag_col)
            col_pos = frag_map_pos(jv, ifrag_col)
          else
            jj = col_lid_work(jv)
            col_pos = col_pos_work(jv)
          end if
          do iv = 1, nrow
            if (allocated(frag_map_count)) then
              ii = frag_map_lid(iv, ifrag_row)
            else
              ii = row_lid_work(iv)
            end if
            cmat_pack(iv, jv) = block%val(ii, jj, ispin_current)
          end do
          do istate = 1, nstate
            cx_pack(jv, istate) = x(col_pos, istate)
          end do
        end do
        call zgemm('N', 'N', nrow, nstate, ncol, cone, cmat_pack, size(cmat_pack, 1), &
                   cx_pack, size(cx_pack, 1), czero, cy_pack, size(cy_pack, 1))
        do iv = 1, nrow
          if (allocated(frag_map_count)) then
            row_pos = frag_map_pos(iv, ifrag_row)
          else
            row_pos = row_pos_work(iv)
          end if
          do istate = 1, nstate
            y(row_pos, istate) = y(row_pos, istate) + cy_pack(iv, istate)
          end do
        end do
        return
      end if
!$omp parallel do if(nrow*ncol*nstate >= 8192) private(iv,jv,ii,jj,row_pos,col_pos,val,istate) schedule(static)
      do iv = 1, nrow
        if (allocated(frag_map_count)) then
          ii = frag_map_lid(iv, ifrag_row)
          row_pos = frag_map_pos(iv, ifrag_row)
        else
          ii = row_lid_work(iv)
          row_pos = row_pos_work(iv)
        end if
        do jv = 1, ncol
          if (allocated(frag_map_count)) then
            jj = frag_map_lid(jv, ifrag_col)
            col_pos = frag_map_pos(jv, ifrag_col)
          else
            jj = col_lid_work(jv)
            col_pos = col_pos_work(jv)
          end if
          val = block%val(ii, jj, ispin_current)
          if (abs(val) == 0.0d0) cycle
          !$omp simd
          do istate = 1, nstate
            y(row_pos, istate) = y(row_pos, istate) + val * x(col_pos, istate)
          end do
        end do
      end do
!$omp end parallel do
    end subroutine apply_one_complex_block

    subroutine add_diamagnetic_term(nstate_blk_current)
      integer, intent(in) :: nstate_blk_current
      integer :: irow

      if (abs(A_squared) < 1.0d-30) return
      do irow = 1, size(needed_row_ids)
        if (.not. row_is_output(needed_row_ids(irow))) cycle
        if (.not. row_is_owned(needed_row_ids(irow), ispin)) cycle
        h_work(irow, 1:nstate_blk_current) = h_work(irow, 1:nstate_blk_current) + &
                                             0.5d0 * A_squared * coef_work(irow, 1:nstate_blk_current)
      end do
    end subroutine add_diamagnetic_term

    subroutine subtract_static_seed_phase(ispin_current, state_s_current, nstate_blk_current)
      integer, intent(in) :: ispin_current, state_s_current, nstate_blk_current
      integer :: istate, state_global, irow
      real(8) :: eps

      if (.not. dg_frag%dc_lcfo_seed_basis_cleaned) return
      if (dg_frag%identity_seed_coefficients) return
      if (.not. allocated(dg_frag%esp)) then
        stop "DG-Fragment RT: DC-LCFO seed requires exported eigenvalues"
      end if
      if (ispin_current < 1 .or. ispin_current > size(dg_frag%esp, 2)) return
      do istate = 1, nstate_blk_current
        state_global = state_s_current + istate - 1
        if (state_global < 1 .or. state_global > size(dg_frag%esp, 1)) cycle
        eps = dg_frag%esp(state_global, ispin_current)
        if (eps == 0.0d0) cycle
        do irow = 1, size(needed_row_ids)
          h_work(irow, istate) = h_work(irow, istate) - eps * coef_work(irow, istate)
        end do
      end do
    end subroutine subtract_static_seed_phase

    subroutine apply_momentum_blocks_compact(ispin_current, scale_vec, x, y)
      integer, intent(in) :: ispin_current
      real(8), intent(in) :: scale_vec(3)
      complex(8), intent(in) :: x(:, :)
      complex(8), intent(inout) :: y(:, :)
      integer :: iblk, idir, iactive, iv, jv, ii, jj, istate, ifrag_row, ifrag_col
      integer :: nrow, ncol, nstate, row_pos, col_pos, nmax, nactive
      integer :: active_dirs(3)
      integer(8) :: block_ops
      real(8) :: scale, val
      real(8), parameter :: one = 1.0d0, zero = 0.0d0

      nstate = min(size(x, 2), size(y, 2))
      nactive = 0
      do idir = 1, 3
        if (abs(scale_vec(idir)) < 1.0d-30) cycle
        nactive = nactive + 1
        active_dirs(nactive) = idir
      end do
      if (nactive <= 0) return
      do iblk = 1, dg_frag%n_momentum_blocks
        if (.not. allocated(dg_frag%momentum_blocks(iblk)%val)) cycle
        ifrag_row = dg_frag%momentum_blocks(iblk)%ifrag_row
        ifrag_col = dg_frag%momentum_blocks(iblk)%ifrag_col
        if (.not. fragment_row_is_local(ifrag_row)) cycle
        if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) cycle
        if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
        if (allocated(frag_map_count)) then
          nrow = frag_map_count(ifrag_row)
          ncol = frag_map_count(ifrag_col)
        else
          nmax = max(dg_frag%n_basis(ifrag_row, ispin_current), dg_frag%n_basis(ifrag_col, ispin_current))
          call ensure_index_work_arrays(max(1, nmax))
          call build_compact_basis_map(ifrag_row, ispin_current, row_lid_work, row_pos_work, nrow)
          call build_compact_basis_map(ifrag_col, ispin_current, col_lid_work, col_pos_work, ncol)
        end if
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_ops = int(nrow, kind=8) * int(ncol, kind=8) * int(nstate, kind=8)
        if (nrow >= 64 .and. ncol >= 64 .and. block_ops >= derivative_real_blas_min_ops) then
          call ensure_dgemm_pack_arrays(nrow, ncol, nstate)
          do jv = 1, ncol
            if (allocated(frag_map_count)) then
              jj = frag_map_lid(jv, ifrag_col)
              col_pos = frag_map_pos(jv, ifrag_col)
            else
              jj = col_lid_work(jv)
              col_pos = col_pos_work(jv)
            end if
            do iv = 1, nrow
              if (allocated(frag_map_count)) then
                ii = frag_map_lid(iv, ifrag_row)
              else
                ii = row_lid_work(iv)
              end if
              val = 0.0d0
              do iactive = 1, nactive
                idir = active_dirs(iactive)
                val = val + scale_vec(idir) * dg_frag%momentum_blocks(iblk)%val(idir, ii, jj, ispin_current)
              end do
              rmat_pack(iv, jv) = val
            end do
            do istate = 1, nstate
              rx_pack(jv, istate) = real(x(col_pos, istate), kind=8)
              ix_pack(jv, istate) = aimag(x(col_pos, istate))
            end do
          end do
          call dgemm('N', 'N', nrow, nstate, ncol, one, rmat_pack, size(rmat_pack, 1), &
                     rx_pack, size(rx_pack, 1), zero, ry_pack, size(ry_pack, 1))
          call dgemm('N', 'N', nrow, nstate, ncol, one, rmat_pack, size(rmat_pack, 1), &
                     ix_pack, size(ix_pack, 1), zero, iy_pack, size(iy_pack, 1))
          do iv = 1, nrow
            if (allocated(frag_map_count)) then
              row_pos = frag_map_pos(iv, ifrag_row)
            else
              row_pos = row_pos_work(iv)
            end if
            do istate = 1, nstate
              y(row_pos, istate) = y(row_pos, istate) + cmplx(ry_pack(iv, istate), iy_pack(iv, istate), kind=8)
            end do
          end do
          cycle
        end if
!$omp parallel do if(nrow*ncol*nstate*nactive >= 8192) private(iv,jv,iactive,idir,ii,jj,row_pos,col_pos,scale,val,istate) schedule(static)
        do iv = 1, nrow
          if (allocated(frag_map_count)) then
            ii = frag_map_lid(iv, ifrag_row)
            row_pos = frag_map_pos(iv, ifrag_row)
          else
            ii = row_lid_work(iv)
            row_pos = row_pos_work(iv)
          end if
          do iactive = 1, nactive
            idir = active_dirs(iactive)
            scale = scale_vec(idir)
            do jv = 1, ncol
              if (allocated(frag_map_count)) then
                jj = frag_map_lid(jv, ifrag_col)
                col_pos = frag_map_pos(jv, ifrag_col)
              else
                jj = col_lid_work(jv)
                col_pos = col_pos_work(jv)
              end if
              val = dg_frag%momentum_blocks(iblk)%val(idir, ii, jj, ispin_current)
              if (val == 0.0d0) cycle
              !$omp simd
              do istate = 1, nstate
                y(row_pos, istate) = y(row_pos, istate) + scale * val * x(col_pos, istate)
              end do
            end do
          end do
        end do
!$omp end parallel do
      end do
    end subroutine apply_momentum_blocks_compact

    logical function row_is_owned(global_idx, ispin_current) result(is_owned)
      integer, intent(in) :: global_idx, ispin_current

      if (dg_frag%coef_state_block_mode) then
        if (.not. allocated(dg_frag%coef_global_to_local)) then
          is_owned = .false.
        else if (global_idx < 1 .or. global_idx > size(dg_frag%coef_global_to_local, 1)) then
          is_owned = .false.
        else
          is_owned = (dg_frag%coef_global_to_local(global_idx, ispin_current) > 0)
        end if
      else if (allocated(dg_frag%coef_owner)) then
        if (global_idx < 1 .or. global_idx > size(dg_frag%coef_owner, 1)) then
          is_owned = .false.
        else
          is_owned = (dg_frag%coef_owner(global_idx, ispin_current) == dg_frag%id)
        end if
      else
        is_owned = .true.
      end if
    end function row_is_owned

    logical function row_is_output(global_idx) result(is_output)
      integer, intent(in) :: global_idx

      if (.not. allocated(output_row_needed)) then
        is_output = .false.
      else if (global_idx < 1 .or. global_idx > size(output_row_needed)) then
        is_output = .false.
      else
        is_output = output_row_needed(global_idx)
      end if
    end function row_is_output

    subroutine scatter_owned_derivative(ispin_current, nstate_blk_current, state_out_s_current, state_out_e_current)
      integer, intent(in) :: ispin_current, nstate_blk_current, state_out_s_current, state_out_e_current
      integer :: irow, global_idx, local_idx

      do irow = 1, size(needed_row_ids)
        global_idx = needed_row_ids(irow)
        if (.not. row_is_output(global_idx)) cycle
        if (.not. row_is_owned(global_idx, ispin_current)) cycle
        if (output_is_local_rows) then
          local_idx = dg_frag%coef_global_to_local(global_idx, ispin_current)
          if (local_idx < 1 .or. local_idx > size(dcoef_dt, 1)) cycle
          dcoef_dt(local_idx, state_out_s_current:state_out_e_current, ispin_current) = &
            rhs_work(irow, 1:nstate_blk_current)
        else
          if (global_idx < 1 .or. global_idx > size(dcoef_dt, 1)) cycle
          dcoef_dt(global_idx, state_out_s_current:state_out_e_current, ispin_current) = &
            rhs_work(irow, 1:nstate_blk_current)
        end if
      end do
    end subroutine scatter_owned_derivative

  end subroutine calculate_time_derivative
