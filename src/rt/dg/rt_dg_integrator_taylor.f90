  subroutine time_evolution_taylor4pc(dg_frag, system, info, rt, itt, dt, &
                                      lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
    rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use salmon_global, only: yn_fix_func, yn_dg_length_gauge
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional
    use rt_dg_fragment_types, only: s_dg_fragment_rt
    use rt_dg_fragment_ops, only: ensure_nonlocal_pp_matrix_A
    use misc_routines, only: get_wtime
    use communication, only: comm_get_min
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(inout) :: system
    type(s_parallel_info),  intent(in)    :: info
    type(s_rt),             intent(inout) :: rt
    integer,                intent(in)    :: itt
    real(8),                intent(in)    :: dt
    type(s_rgrid),          intent(in)    :: lg, mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_xc_functional),  intent(in)    :: xc_func
    type(s_sendrecv_grid),  intent(inout) :: srg, srg_scalar
    type(s_reciprocal_grid),intent(in)    :: fg
    type(s_poisson),        intent(inout) :: poisson
    type(s_pp_info),        intent(in)    :: pp
    type(s_pp_grid),        intent(in)    :: ppg
    type(s_pp_nlcc),        intent(in)    :: ppn
    type(s_scalar),         intent(inout) :: rho, Vh, Vpsl
    type(s_scalar),         intent(inout) :: rho_s(system%nspin), Vxc(system%nspin)
    type(s_dft_energy),     intent(inout) :: energy

    complex(8), allocatable, save :: coef_base(:,:,:)
    complex(8), allocatable, save :: coef_next(:,:,:)
    complex(8), allocatable, save :: term_blk(:,:,:)
    complex(8), allocatable, save :: deriv_blk(:,:,:)
    complex(8), allocatable, save :: result_blk(:,:,:)
    complex(8), allocatable, save :: work_in(:,:,:)
    complex(8), allocatable, save :: work_out(:,:,:)
    integer, allocatable, save :: owned_pw_row_ids_taylor(:)
    integer :: n, nrow_taylor, nstate_prop, nstate_work, state_global_first, state_global_last
    integer :: state0, state_s, state_e, nstate_blk
    integer :: n_pw, n_pw_owned
    integer :: it0, it1, nsub_taylor
    integer, parameter :: state_work_target_mb = 512
    integer, parameter :: state_work_vectors = 8
    real(8), parameter :: taylor4pc_dt_substep_target = 1.0d-2
    integer, parameter :: taylor4pc_substeps_default = 2
    integer, save :: cfg_state_work_target_mb = state_work_target_mb
    integer, save :: trace_taylor_call_count = 0
    logical, save :: taylor_env_initialized = .false.
    logical, save :: enable_taylor_timing = .false.
    integer(8) :: target_bytes, bytes_per_state
    real(8) :: nstate_work_min_in(1), nstate_work_min_out(1)
    real(8) :: Ac_mid(3), Ac_ham(3), E_mid(3)
    real(8) :: t_taylor0, t_taylor1, time_taylor_apply
    logical :: trace_taylor
    logical :: has_state_work
    logical :: use_length_gauge
    logical :: use_pw_taylor
    integer :: trace_call_id

    n = size(dg_frag%coef, 1)
    nrow_taylor = n
    n_pw = 0
    n_pw_owned = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    use_pw_taylor = (n_pw > 0)
    if (use_pw_taylor) then
      if (dg_frag%coef_state_block_mode) then
        stop "DG Taylor4-PC PW path does not yet support state-block coefficient ownership"
      end if
      call build_owned_pw_rows_for_taylor()
      nrow_taylor = n + n_pw_owned
    end if

    nstate_prop = dg_frag%nstate_tot
    if (allocated(dg_frag%nocc_spin)) then
      nstate_prop = min(dg_frag%nstate_tot, max(1, maxval(dg_frag%nocc_spin(1:dg_frag%nspin))))
    end if
    state_global_first = 1
    state_global_last = nstate_prop
    if (dg_frag%coef_state_block_mode) then
      state_global_first = max(1, dg_frag%coef_state_start)
      state_global_last = min(nstate_prop, dg_frag%coef_state_end)
      nstate_prop = max(0, state_global_last - state_global_first + 1)
    end if
    has_state_work = (nstate_prop > 0)
    it0 = max(lbound(rt%Ac_tot, 2), itt - 1)
    it1 = min(ubound(rt%Ac_tot, 2), itt)
    Ac_mid(:) = 0.5d0 * (rt%Ac_tot(:, it0) + rt%Ac_tot(:, it1))
    E_mid(:) = 0.5d0 * (rt%E_tot(:, it0) + rt%E_tot(:, it1))
    use_length_gauge = (yn_dg_length_gauge == 'y')
    Ac_ham(:) = Ac_mid(:)
    if (use_length_gauge) Ac_ham(:) = 0.0d0
    nsub_taylor = max(taylor4pc_substeps_default, ceiling(abs(dt) / taylor4pc_dt_substep_target))
    call initialize_taylor_runtime_options()

    if (yn_fix_func == 'n') then
      call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_ham, &
                                            lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                            rho, rho_s, Vh, Vxc, Vpsl, energy, .false.)
    end if
    if (ppg%Nlma > 0 .and. allocated(ppg%uV)) then
      call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_ham, .false.)
    end if

    if (.not. has_state_work) then
      return
    end if

    target_bytes = int(cfg_state_work_target_mb, kind=8) * 1024_8 * 1024_8
    bytes_per_state = 16_8 * int(max(1, nrow_taylor), kind=8) * int(max(1, dg_frag%nspin), kind=8) * &
                      int(state_work_vectors, kind=8)
    nstate_work = max(1, min(nstate_prop, int(max(1_8, target_bytes / max(1_8, bytes_per_state)))))
    if (use_pw_taylor) then
      nstate_work_min_in(1) = dble(nstate_work)
      call comm_get_min(nstate_work_min_in, nstate_work_min_out, 1, dg_frag%icomm)
      nstate_work = max(1, min(nstate_prop, int(nstate_work_min_out(1))))
    end if
    time_taylor_apply = 0.0d0
    trace_taylor = .false.
    trace_call_id = 0
    if (enable_taylor_timing .and. itt == 1 .and. dg_frag%id == 0 .and. trace_taylor_call_count < 5) then
      trace_taylor_call_count = trace_taylor_call_count + 1
      trace_call_id = trace_taylor_call_count
      trace_taylor = .true.
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
        '[DG-TAYLOR] enter itt=', itt, ' call=', trace_call_id, &
        ' nstate_prop=', nstate_prop, ' nstate_work=', nstate_work, &
        ' nsub=', nsub_taylor, ' state_mb=', cfg_state_work_target_mb
      flush(6)
    end if

    t_taylor0 = get_wtime()
    call ensure_state_arrays(nrow_taylor, nstate_prop, dg_frag%nspin)
    call ensure_block_arrays(nrow_taylor, nstate_work, dg_frag%nspin)
    call pack_taylor_coefficients(coef_base, nstate_prop)
    call unpack_taylor_coefficients(coef_base, nstate_prop)
    call apply_taylor4_with_current_hamiltonian(coef_base, Ac_ham, E_mid, coef_next, nstate_prop, nstate_work, &
                                                state_global_first)
    call unpack_taylor_coefficients(coef_next, nstate_prop)
    t_taylor1 = get_wtime()
    time_taylor_apply = t_taylor1 - t_taylor0
    if (trace_taylor) then
      write(*,'(1x,a,i0,a,1pe13.5)') &
        '[DG-TAYLOR] timing call=', trace_call_id, &
        ' apply_total=', time_taylor_apply
      flush(6)
    end if

  contains

    subroutine initialize_taylor_runtime_options()
      character(32) :: env_value
      integer :: env_status, env_len, env_int, read_status

      if (taylor_env_initialized) return
      env_value = ''
      call get_environment_variable('SALMON_DG_TAYLOR_TIMING', env_value, length=env_len, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(env_value(1:env_len))))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enable_taylor_timing = .true.
        end select
      end if
      env_value = ''
      call get_environment_variable('SALMON_DG_TAYLOR_STATE_MB', env_value, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_value(1:env_len), *, iostat=read_status) env_int
        if (read_status == 0) cfg_state_work_target_mb = max(1, env_int)
      end if
      taylor_env_initialized = .true.
    end subroutine initialize_taylor_runtime_options

    subroutine ensure_state_arrays(nrow, ncol, nspin)
      integer, intent(in) :: nrow, ncol, nspin

      if (.not. allocated(coef_base)) then
        allocate(coef_base(nrow, ncol, nspin), coef_next(nrow, ncol, nspin), &
                 work_in(nrow, ncol, nspin), work_out(nrow, ncol, nspin))
      else if (size(coef_base, 1) /= nrow .or. size(coef_base, 2) /= ncol .or. size(coef_base, 3) /= nspin) then
        deallocate(coef_base, coef_next, work_in, work_out)
        allocate(coef_base(nrow, ncol, nspin), coef_next(nrow, ncol, nspin), &
                 work_in(nrow, ncol, nspin), work_out(nrow, ncol, nspin))
      end if
    end subroutine ensure_state_arrays

    subroutine ensure_block_arrays(nrow, ncol, nspin)
      integer, intent(in) :: nrow, ncol, nspin

      if (.not. allocated(term_blk)) then
        allocate(term_blk(nrow, ncol, nspin), deriv_blk(nrow, ncol, nspin), result_blk(nrow, ncol, nspin))
      else if (size(term_blk, 1) /= nrow .or. size(term_blk, 2) /= ncol .or. size(term_blk, 3) /= nspin) then
        deallocate(term_blk, deriv_blk, result_blk)
        allocate(term_blk(nrow, ncol, nspin), deriv_blk(nrow, ncol, nspin), result_blk(nrow, ncol, nspin))
      end if
    end subroutine ensure_block_arrays

    subroutine build_owned_pw_rows_for_taylor()
      integer :: ipw, nowned

      if (allocated(owned_pw_row_ids_taylor)) deallocate(owned_pw_row_ids_taylor)
      if (n_pw <= 0) return
      if (.not. allocated(dg_frag%coef_pw_owner)) then
        stop "DG Taylor4-PC PW path requires PW row owners"
      end if
      nowned = 0
      if (.not. allocated(dg_frag%fp_local_pw_ids)) then
        stop "DG Taylor4-PC PW path requires prepared local PW row ids"
      end if
      do ipw = 1, size(dg_frag%fp_local_pw_ids)
        if (dg_frag%fp_local_pw_ids(ipw) < 1 .or. dg_frag%fp_local_pw_ids(ipw) > size(dg_frag%coef_pw_owner)) cycle
        if (dg_frag%coef_pw_owner(dg_frag%fp_local_pw_ids(ipw)) == dg_frag%id) nowned = nowned + 1
      end do
      n_pw_owned = nowned
      allocate(owned_pw_row_ids_taylor(max(1, n_pw_owned)))
      nowned = 0
      do ipw = 1, size(dg_frag%fp_local_pw_ids)
        if (dg_frag%fp_local_pw_ids(ipw) < 1 .or. dg_frag%fp_local_pw_ids(ipw) > size(dg_frag%coef_pw_owner)) cycle
        if (dg_frag%coef_pw_owner(dg_frag%fp_local_pw_ids(ipw)) /= dg_frag%id) cycle
        nowned = nowned + 1
        owned_pw_row_ids_taylor(nowned) = dg_frag%fp_local_pw_ids(ipw)
      end do
    end subroutine build_owned_pw_rows_for_taylor

    subroutine pack_taylor_coefficients(buffer, nstate_current)
      complex(8), intent(out) :: buffer(:, :, :)
      integer, intent(in) :: nstate_current
      integer :: ispin, ipw_slot, pw_row, state_col0, state_col1

      buffer(:, :, :) = (0.0d0, 0.0d0)
      buffer(1:n, 1:nstate_current, :) = dg_frag%coef(1:n, 1:nstate_current, :)
      if (.not. use_pw_taylor) return
      if (.not. allocated(owned_pw_row_ids_taylor)) return
      state_col0 = state_global_first
      state_col1 = state_global_first + nstate_current - 1
      if (state_col1 > size(dg_frag%coef_pw, 2)) then
        stop "DG Taylor4-PC PW coefficient state range exceeds coef_pw columns"
      end if
      do ispin = 1, dg_frag%nspin
        do ipw_slot = 1, n_pw_owned
          pw_row = owned_pw_row_ids_taylor(ipw_slot)
          if (pw_row < 1 .or. pw_row > size(dg_frag%coef_pw, 1)) cycle
          buffer(n + ipw_slot, 1:nstate_current, ispin) = dg_frag%coef_pw(pw_row, state_col0:state_col1, ispin)
        end do
      end do
    end subroutine pack_taylor_coefficients

    subroutine unpack_taylor_coefficients(buffer, nstate_current)
      complex(8), intent(in) :: buffer(:, :, :)
      integer, intent(in) :: nstate_current
      integer :: ispin, ipw_slot, pw_row, state_col0, state_col1

      dg_frag%coef(1:n, 1:nstate_current, :) = buffer(1:n, 1:nstate_current, :)
      if (.not. use_pw_taylor) return
      if (.not. allocated(owned_pw_row_ids_taylor)) return
      state_col0 = state_global_first
      state_col1 = state_global_first + nstate_current - 1
      if (state_col1 > size(dg_frag%coef_pw, 2)) then
        stop "DG Taylor4-PC PW coefficient state range exceeds coef_pw columns"
      end if
      do ispin = 1, dg_frag%nspin
        do ipw_slot = 1, n_pw_owned
          pw_row = owned_pw_row_ids_taylor(ipw_slot)
          if (pw_row < 1 .or. pw_row > size(dg_frag%coef_pw, 1)) cycle
          dg_frag%coef_pw(pw_row, state_col0:state_col1, ispin) = buffer(n + ipw_slot, 1:nstate_current, ispin)
        end do
      end do
    end subroutine unpack_taylor_coefficients

    subroutine apply_taylor4_with_current_hamiltonian(coef_in, Ac_use, E_use, coef_out, nstate_current, nstate_block, &
                                                      state_global_offset)
      complex(8), intent(in) :: coef_in(:, :, :)
      real(8), intent(in) :: Ac_use(3)
      real(8), intent(in) :: E_use(3)
      complex(8), intent(out) :: coef_out(:, :, :)
      integer, intent(in) :: nstate_current, nstate_block, state_global_offset
      integer :: isub

      work_in(:, :, :) = coef_in(:, :, :)
      do isub = 1, nsub_taylor
        call apply_taylor4_single_step(work_in, Ac_use, E_use, work_out, nstate_current, nstate_block, &
                                       dt / dble(nsub_taylor), state_global_offset)
        work_in(:, :, :) = work_out(:, :, :)
      end do
      coef_out(:, :, :) = work_in(:, :, :)
      call unpack_taylor_coefficients(coef_out, nstate_current)
    end subroutine apply_taylor4_with_current_hamiltonian

    subroutine apply_taylor4_single_step(coef_in, Ac_use, E_use, coef_out, nstate_current, nstate_block, dt_step, &
                                         state_global_offset)
      complex(8), intent(in) :: coef_in(:, :, :)
      real(8), intent(in) :: Ac_use(3)
      real(8), intent(in) :: E_use(3)
      complex(8), intent(out) :: coef_out(:, :, :)
      integer, intent(in) :: nstate_current, nstate_block, state_global_offset
      real(8), intent(in) :: dt_step
      integer :: order
      real(8) :: coeff

      coef_out(:, :, :) = coef_in(:, :, :)
      do state0 = 1, nstate_current, nstate_block
        nstate_blk = min(nstate_block, nstate_current - state0 + 1)
        state_s = state_global_offset + state0 - 1
        state_e = state_s + nstate_blk - 1
        term_blk(:, 1:nstate_blk, :) = coef_in(:, state0:state0+nstate_blk-1, :)
        result_blk(:, 1:nstate_blk, :) = coef_in(:, state0:state0+nstate_blk-1, :)
        coeff = 1.0d0
        do order = 1, 4
          call unpack_taylor_state_block(term_blk(:, 1:nstate_blk, :), state0, nstate_blk)
          deriv_blk(:, 1:nstate_blk, :) = (0.0d0, 0.0d0)
          call calculate_time_derivative(dg_frag, system, mg, ppg, Ac_use, &
                                         deriv_blk(:, 1:nstate_blk, :), state_s, state_e, E_use, .true.)
          call check_finite_taylor_block(order, state_s, state_e, deriv_blk(:, 1:nstate_blk, :))
          coeff = coeff * dt_step / dble(order)
          result_blk(:, 1:nstate_blk, :) = result_blk(:, 1:nstate_blk, :) + &
                                           coeff * deriv_blk(:, 1:nstate_blk, :)
          term_blk(:, 1:nstate_blk, :) = deriv_blk(:, 1:nstate_blk, :)
        end do
        coef_out(:, state0:state0+nstate_blk-1, :) = result_blk(:, 1:nstate_blk, :)
      end do
    end subroutine apply_taylor4_single_step

    subroutine unpack_taylor_state_block(buffer, state_local_first, nstate_current)
      complex(8), intent(in) :: buffer(:, :, :)
      integer, intent(in) :: state_local_first, nstate_current
      integer :: ispin, ipw_slot, pw_row, state_local_last, state_col0, state_col1

      state_local_last = state_local_first + nstate_current - 1
      dg_frag%coef(1:n, state_local_first:state_local_last, :) = buffer(1:n, 1:nstate_current, :)
      if (.not. use_pw_taylor) return
      if (.not. allocated(owned_pw_row_ids_taylor)) return
      state_col0 = state_global_first + state_local_first - 1
      state_col1 = state_col0 + nstate_current - 1
      if (state_col1 > size(dg_frag%coef_pw, 2)) then
        stop "DG Taylor4-PC PW coefficient state block exceeds coef_pw columns"
      end if
      do ispin = 1, dg_frag%nspin
        do ipw_slot = 1, n_pw_owned
          pw_row = owned_pw_row_ids_taylor(ipw_slot)
          if (pw_row < 1 .or. pw_row > size(dg_frag%coef_pw, 1)) cycle
          dg_frag%coef_pw(pw_row, state_col0:state_col1, ispin) = &
            buffer(n + ipw_slot, 1:nstate_current, ispin)
        end do
      end do
    end subroutine unpack_taylor_state_block

    subroutine check_finite_taylor_block(order, state_s_in, state_e_in, block)
      integer, intent(in) :: order, state_s_in, state_e_in
      complex(8), intent(in) :: block(:, :, :)

      if (any(real(block, kind=8) /= real(block, kind=8)) .or. &
          any(aimag(block) /= aimag(block))) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') '[FATAL] NaN in DG Taylor derivative: rank=', &
          dg_frag%id, ' order=', order, ' state_s=', state_s_in, ' state_e=', state_e_in
        stop 'DG-Fragment RT: NaN Taylor derivative block'
      end if
    end subroutine check_finite_taylor_block

  end subroutine time_evolution_taylor4pc
