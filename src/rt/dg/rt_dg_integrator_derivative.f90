  subroutine calculate_time_derivative(dg_frag, system, mg, ppg, Ac_tot, dcoef_dt, state_start, state_end)
    use structures
    use rt_dg_fragment_types, only: s_dg_fragment_rt
    use rt_dg_fragment_ops, only: apply_momentum_blocks, apply_matrix_blocks_batch, &
                                  ensure_nonlocal_projector_overlap_cache, apply_nonlocal_projector_overlap_batch, &
                                  solve_overlap_operator_batch, rebuild_local_h_block_ids, fetch_remote_coef_rows
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    real(8),                intent(in)    :: Ac_tot(3)
    complex(8),             intent(out)   :: dcoef_dt(:,:,:)
    integer, optional,      intent(in)    :: state_start, state_end

    integer :: io, ispin, ispin_other
    integer :: n_frag, n_tot
    integer :: n_s
    integer :: state_first, state_last, state0, state_s, state_e
    integer :: state_out_s, state_out_e, nstate_blk, nstate_work
    integer :: irow
    integer, parameter :: state_work_target_mb = 512
    integer, parameter :: state_work_vectors = 5
    integer(8) :: target_bytes, bytes_per_state
    real(8) :: A_squared
    complex(8), parameter :: zi = (0.0d0, 1.0d0)
    logical :: has_nonlocal, has_so_nonlocal, output_is_block, output_is_local_rows
    logical :: use_local_h_blocks
    logical, parameter :: trace_first_block = .false.
    complex(8), allocatable, save :: coef_all(:,:), rhs_all(:,:)
    complex(8), allocatable, save :: dcoef_dt_h0(:,:), dcoef_dt_m(:,:), rhs_in(:,:)
    complex(8), allocatable :: coef_frag_other(:,:)
    complex(8), allocatable :: coef_needed(:,:)
    integer, allocatable :: needed_row_ids(:)

    dcoef_dt = (0.0d0, 0.0d0)

    if (dg_frag%use_plane_wave_basis .or. allocated(dg_frag%coef_pw)) then
      stop "DG derivative now supports the pure fragment block-sparse route only"
    end if
    if (.not. allocated(dg_frag%H_mat_blocks)) then
      stop "DG derivative requires block-sparse H blocks"
    end if
    if (.not. allocated(dg_frag%momentum_blocks)) then
      stop "DG derivative requires block-sparse momentum blocks"
    end if
    call rebuild_local_h_block_ids(dg_frag)

    n_frag = dg_frag%n_mat_max
    if (n_frag <= 0) return
    n_tot = n_frag
    if (dg_frag%nstate_tot > size(dg_frag%coef, 2)) then
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

    target_bytes = int(state_work_target_mb, kind=8) * 1024_8 * 1024_8
    bytes_per_state = 16_8 * int(max(1, n_tot), kind=8) * int(state_work_vectors, kind=8)
    nstate_work = max(1, min(state_last - state_first + 1, int(max(1_8, target_bytes / max(1_8, bytes_per_state)))))

    if (.not. allocated(coef_all)) then
      allocate(coef_all(n_tot, nstate_work), rhs_all(n_tot, nstate_work))
    else if (size(coef_all, 1) /= n_tot .or. size(coef_all, 2) /= nstate_work) then
      deallocate(coef_all, rhs_all)
      allocate(coef_all(n_tot, nstate_work), rhs_all(n_tot, nstate_work))
    end if
    if (.not. allocated(dcoef_dt_h0)) then
      allocate(dcoef_dt_h0(n_tot, nstate_work), dcoef_dt_m(n_tot, nstate_work))
    else if (size(dcoef_dt_h0, 1) /= n_tot .or. size(dcoef_dt_h0, 2) /= nstate_work) then
      deallocate(dcoef_dt_h0, dcoef_dt_m)
      allocate(dcoef_dt_h0(n_tot, nstate_work), dcoef_dt_m(n_tot, nstate_work))
    end if

    A_squared = Ac_tot(1)**2 + Ac_tot(2)**2 + Ac_tot(3)**2
    if (A_squared /= A_squared) stop "NaN in DG derivative vector potential"
    has_nonlocal = (ppg%Nlma > 0 .and. allocated(ppg%uV))
    has_so_nonlocal = (allocated(ppg%uv_so) .and. allocated(dg_frag%phi_frag_c) .and. dg_frag%nspin == 2)
    use_local_h_blocks = dg_frag%H_local_block_ids_valid .and. allocated(dg_frag%H_local_block_ids)
    if (.not. use_local_h_blocks) then
      write(*,'(1x,a,i0)') '[FATAL] DG derivative requires row-owner-local H block ids: rank=', dg_frag%id
      stop 'DG-Fragment RT: missing row-owner-local H block ids'
    end if
    if (has_nonlocal) then
      call ensure_nonlocal_projector_overlap_cache(dg_frag, mg, ppg, system%nspin, system%hvol, Ac_tot)
      if (.not. dg_frag%has_nl_projector_cache .or. .not. allocated(dg_frag%nl_projector_overlap) .or. &
          .not. allocated(dg_frag%nl_projector_overlap_halo)) &
        stop 'DG-Fragment RT: missing nonlocal PP projector cache'
    end if
    if (has_so_nonlocal) then
      stop 'DG-Fragment RT: SO nonlocal projector route is disabled in block-sparse RT'
    end if
    call build_needed_coefficient_rows()
    if (has_so_nonlocal) then
      if (.not. allocated(coef_frag_other)) then
        allocate(coef_frag_other(n_frag, nstate_work))
      else if (size(coef_frag_other, 1) /= n_frag .or. size(coef_frag_other, 2) /= nstate_work) then
        deallocate(coef_frag_other)
        allocate(coef_frag_other(n_frag, nstate_work))
      end if
    end if

    n_s = 0
    if (allocated(dg_frag%S_mat_prop_blocks) .or. allocated(dg_frag%S_mat_prop_c) .or. &
        allocated(dg_frag%S_mat_prop) .or. allocated(dg_frag%S_mat_c) .or. allocated(dg_frag%S_mat)) then
      n_s = n_frag
      if (.not. allocated(rhs_in)) then
        allocate(rhs_in(n_s, nstate_work))
      else if (size(rhs_in, 1) /= n_s .or. size(rhs_in, 2) /= nstate_work) then
        deallocate(rhs_in)
        allocate(rhs_in(n_s, nstate_work))
      end if
    end if

    do ispin = 1, dg_frag%nspin
      do state0 = state_first, state_last, nstate_work
        nstate_blk = min(nstate_work, state_last - state0 + 1)
        state_s = state0
        state_e = state0 + nstate_blk - 1
        state_out_s = merge(state_s - state_first + 1, state_s, output_is_block)
        state_out_e = merge(state_e - state_first + 1, state_e, output_is_block)

        coef_all(:, :) = (0.0d0, 0.0d0)
        if (trace_first_block .and. ispin == 1 .and. state0 == state_first) then
          write(*,'(1x,a,i0,a,i0,a,i0)') '[DG-DERIV-FIRST] gather begin n_frag=', n_frag, &
            ' state_e=', state_e, ' nstate_blk=', nstate_blk
          flush(6)
        end if
        call gather_needed_fragment_coef(ispin, nstate_blk, state_s, state_e, coef_all(:, 1:nstate_blk))
        if (trace_first_block .and. ispin == 1 .and. state0 == state_first) then
          write(*,'(1x,a)') '[DG-DERIV-FIRST] gather done'
          flush(6)
        end if
        call gather_so_partner_if_needed(ispin, nstate_blk, state_s, state_e)

        dcoef_dt_h0(:, :) = (0.0d0, 0.0d0)
        if (trace_first_block .and. ispin == 1 .and. state0 == state_first) then
          write(*,'(1x,a)') '[DG-DERIV-FIRST] H apply begin'
          flush(6)
        end if
        call apply_h_blocks_for_route(ispin, nstate_blk)
        if (trace_first_block .and. ispin == 1 .and. state0 == state_first) then
          write(*,'(1x,a)') '[DG-DERIV-FIRST] H apply done'
          flush(6)
        end if
        if (trace_first_block .and. ispin == 1 .and. state0 == state_first) then
          write(*,'(1x,a)') '[DG-DERIV-FIRST] nonlocal apply begin'
          flush(6)
        end if
        call apply_nonlocal_for_route(ispin, nstate_blk)
        if (trace_first_block .and. ispin == 1 .and. state0 == state_first) then
          write(*,'(1x,a)') '[DG-DERIV-FIRST] nonlocal apply done'
          flush(6)
        end if

        dcoef_dt_h0(1:n_frag, 1:nstate_blk) = dcoef_dt_h0(1:n_frag, 1:nstate_blk) + &
                                              0.5d0 * A_squared * coef_all(1:n_frag, 1:nstate_blk)
        dcoef_dt_h0(1:n_frag, 1:nstate_blk) = -zi * dcoef_dt_h0(1:n_frag, 1:nstate_blk)

        dcoef_dt_m(:, :) = (0.0d0, 0.0d0)
        if (trace_first_block .and. ispin == 1 .and. state0 == state_first) then
          write(*,'(1x,a)') '[DG-DERIV-FIRST] momentum apply begin'
          flush(6)
        end if
        if (allocated(dg_frag%H_local_rows)) then
          call apply_momentum_blocks(dg_frag, ispin, Ac_tot, &
                                     coef_all(1:n_frag, 1:nstate_blk), dcoef_dt_m(1:n_frag, 1:nstate_blk), &
                                     dg_frag%H_local_rows)
        else
          call apply_momentum_blocks(dg_frag, ispin, Ac_tot, &
                                     coef_all(1:n_frag, 1:nstate_blk), dcoef_dt_m(1:n_frag, 1:nstate_blk))
        end if
        if (trace_first_block .and. ispin == 1 .and. state0 == state_first) then
          write(*,'(1x,a)') '[DG-DERIV-FIRST] momentum apply done'
          flush(6)
        end if
        rhs_all(:, 1:nstate_blk) = dcoef_dt_h0(:, 1:nstate_blk) - dcoef_dt_m(:, 1:nstate_blk)
        if (trace_first_block .and. ispin == 1 .and. state0 == state_first) then
          write(*,'(1x,a)') '[DG-DERIV-FIRST] overlap solve begin'
          flush(6)
        end if
        call solve_overlap_for_route(ispin, nstate_blk)
        if (trace_first_block .and. ispin == 1 .and. state0 == state_first) then
          write(*,'(1x,a)') '[DG-DERIV-FIRST] overlap solve done'
          flush(6)
        end if
        call scatter_owned_derivative(ispin, nstate_blk, state_out_s, state_out_e)
        if (trace_first_block .and. ispin == 1 .and. state0 == state_first) then
          write(*,'(1x,a)') '[DG-DERIV-FIRST] scatter done'
          flush(6)
        end if
      end do
    end do

  contains

    subroutine gather_so_partner_if_needed(ispin_current, nstate_blk_current, state_s_current, state_e_current)
      integer, intent(in) :: ispin_current, nstate_blk_current, state_s_current, state_e_current

      if (.not. has_so_nonlocal) return
      ispin_other = 3 - ispin_current
      call gather_needed_fragment_coef(ispin_other, nstate_blk_current, state_s_current, state_e_current, &
                                       coef_frag_other(1:n_frag, 1:nstate_blk_current))
    end subroutine gather_so_partner_if_needed

    subroutine build_needed_coefficient_rows()
      logical, allocatable :: row_needed(:)
      integer :: iblk_idx, iblk, ifrag_row, ifrag_col, global_idx, n_needed

      if (allocated(needed_row_ids)) deallocate(needed_row_ids)
      allocate(row_needed(max(1, n_frag)))
      row_needed(:) = .false.
      do iblk_idx = 1, size(dg_frag%H_local_block_ids)
        iblk = dg_frag%H_local_block_ids(iblk_idx)
        if (iblk < 1 .or. iblk > size(dg_frag%H_mat_blocks)) cycle
        ifrag_row = dg_frag%H_mat_blocks(iblk)%ifrag_row
        ifrag_col = dg_frag%H_mat_blocks(iblk)%ifrag_col
        call mark_fragment_basis_rows(row_needed, ifrag_row)
        call mark_fragment_basis_rows(row_needed, ifrag_col)
      end do
      if (allocated(dg_frag%momentum_blocks) .and. allocated(dg_frag%H_local_rows)) then
        do iblk = 1, dg_frag%n_momentum_blocks
          if (iblk < 1 .or. iblk > size(dg_frag%momentum_blocks)) cycle
          ifrag_row = dg_frag%momentum_blocks(iblk)%ifrag_row
          ifrag_col = dg_frag%momentum_blocks(iblk)%ifrag_col
          if (.not. any(dg_frag%H_local_rows == ifrag_row)) cycle
          call mark_fragment_basis_rows(row_needed, ifrag_row)
          call mark_fragment_basis_rows(row_needed, ifrag_col)
        end do
      end if
      n_needed = count(row_needed(1:n_frag))
      allocate(needed_row_ids(max(1, n_needed)))
      needed_row_ids(:) = 0
      if (n_needed > 0) then
        n_needed = 0
        do global_idx = 1, n_frag
          if (.not. row_needed(global_idx)) cycle
          n_needed = n_needed + 1
          needed_row_ids(n_needed) = global_idx
        end do
      end if
      deallocate(row_needed)
    end subroutine build_needed_coefficient_rows

    subroutine mark_fragment_basis_rows(row_needed, ifrag)
      logical, intent(inout) :: row_needed(:)
      integer, intent(in) :: ifrag
      integer :: ib, ispin_mark, global_idx

      if (ifrag < 1 .or. ifrag > dg_frag%n_frag) return
      do ispin_mark = 1, dg_frag%nspin
        do ib = 1, min(dg_frag%n_basis(ifrag, ispin_mark), size(dg_frag%index_basis, 1))
          global_idx = dg_frag%index_basis(ib, ifrag, ispin_mark)
          if (global_idx < 1 .or. global_idx > size(row_needed)) cycle
          row_needed(global_idx) = .true.
        end do
      end do
    end subroutine mark_fragment_basis_rows

    subroutine gather_needed_fragment_coef(ispin_current, nstate_blk_current, state_s_current, state_e_current, coef_out)
      integer, intent(in) :: ispin_current, nstate_blk_current, state_s_current, state_e_current
      complex(8), intent(inout) :: coef_out(:, :)

      coef_out(:, :) = (0.0d0, 0.0d0)
      if (.not. allocated(needed_row_ids)) then
        write(*,'(1x,a,i0)') '[FATAL] DG derivative missing block-sparse row dependency map: rank=', dg_frag%id
        stop 'DG-Fragment RT: missing block-sparse row dependency map'
      end if
      if (size(needed_row_ids) >= n_frag) then
        write(*,'(1x,a,i0,a,i0,a,i0)') '[FATAL] DG derivative would gather all coefficient rows: rank=', &
          dg_frag%id, ' n_needed=', size(needed_row_ids), ' n_frag=', n_frag
        stop 'DG-Fragment RT: all-row coefficient gather is disabled'
      end if
      if (.not. allocated(coef_needed)) then
        allocate(coef_needed(size(needed_row_ids), nstate_blk_current))
      else if (size(coef_needed, 1) /= size(needed_row_ids) .or. size(coef_needed, 2) /= nstate_blk_current) then
        deallocate(coef_needed)
        allocate(coef_needed(size(needed_row_ids), nstate_blk_current))
      end if
      call fetch_remote_coef_rows(dg_frag, ispin_current, needed_row_ids, coef_needed, state_s_current, state_e_current)
      do irow = 1, size(needed_row_ids)
        if (needed_row_ids(irow) < 1 .or. needed_row_ids(irow) > size(coef_out, 1)) cycle
        coef_out(needed_row_ids(irow), 1:nstate_blk_current) = coef_needed(irow, 1:nstate_blk_current)
      end do
    end subroutine gather_needed_fragment_coef

    subroutine apply_h_blocks_for_route(ispin_current, nstate_blk_current)
      integer, intent(in) :: ispin_current, nstate_blk_current

      if (use_local_h_blocks) then
        call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin_current, &
                                       coef_all(:, 1:nstate_blk_current), dcoef_dt_h0(:, 1:nstate_blk_current), &
                                       dg_frag%H_local_block_ids)
      else
        call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin_current, &
                                       coef_all(:, 1:nstate_blk_current), dcoef_dt_h0(:, 1:nstate_blk_current))
      end if
    end subroutine apply_h_blocks_for_route

    subroutine apply_nonlocal_for_route(ispin_current, nstate_blk_current)
      integer, intent(in) :: ispin_current, nstate_blk_current

      if (has_nonlocal) then
        if (allocated(dg_frag%H_local_rows)) then
          call apply_nonlocal_projector_overlap_batch(dg_frag, mg, ppg, system%hvol, Ac_tot, ispin_current, &
                                                      coef_all(:, 1:nstate_blk_current), &
                                                      dcoef_dt_h0(:, 1:nstate_blk_current), &
                                                      dg_frag%H_local_rows)
        else
          call apply_nonlocal_projector_overlap_batch(dg_frag, mg, ppg, system%hvol, Ac_tot, ispin_current, &
                                                      coef_all(:, 1:nstate_blk_current), &
                                                      dcoef_dt_h0(:, 1:nstate_blk_current))
        end if
      end if
    end subroutine apply_nonlocal_for_route

    subroutine solve_overlap_for_route(ispin_current, nstate_blk_current)
      integer, intent(in) :: ispin_current, nstate_blk_current

      if (n_s <= 0) return
      rhs_in(:, 1:nstate_blk_current) = rhs_all(1:n_s, 1:nstate_blk_current)
      call solve_overlap_operator_batch(dg_frag, ispin_current, &
                                        rhs_in(:, 1:nstate_blk_current), &
                                        rhs_all(1:n_s, 1:nstate_blk_current), .true.)
    end subroutine solve_overlap_for_route

    subroutine scatter_owned_derivative(ispin_current, nstate_blk_current, state_out_s_current, state_out_e_current)
      integer, intent(in) :: ispin_current, nstate_blk_current, state_out_s_current, state_out_e_current
      integer :: local_idx

      if (allocated(dg_frag%coef_owner)) then
        do io = 1, n_frag
          if (dg_frag%coef_owner(io, ispin_current) /= dg_frag%id) cycle
          if (output_is_local_rows) then
            local_idx = dg_frag%coef_global_to_local(io, ispin_current)
            if (local_idx < 1 .or. local_idx > size(dcoef_dt, 1)) cycle
            dcoef_dt(local_idx, state_out_s_current:state_out_e_current, ispin_current) = rhs_all(io, 1:nstate_blk_current)
          else
            dcoef_dt(io, state_out_s_current:state_out_e_current, ispin_current) = rhs_all(io, 1:nstate_blk_current)
          end if
        end do
      else
        dcoef_dt(1:n_frag, state_out_s_current:state_out_e_current, ispin_current) = &
          rhs_all(1:n_frag, 1:nstate_blk_current)
      end if
    end subroutine scatter_owned_derivative

  end subroutine calculate_time_derivative
