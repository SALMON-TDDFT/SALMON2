  subroutine time_evolution_taylor4pc(dg_frag, system, info, rt, itt, dt, &
                                      lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
    rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use salmon_global, only: yn_fix_func
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional
    use rt_dg_fragment_types, only: s_dg_fragment_rt
    use rt_dg_fragment_ops, only: ensure_nonlocal_pp_matrix_A
    use misc_routines, only: get_wtime
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
    integer :: n, nstate_prop, nstate_work, state_global_first, state_global_last
    integer :: state0, state_s, state_e, nstate_blk
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
    real(8) :: Ac_mid(3)
    real(8) :: t_taylor0, t_taylor1, time_taylor_apply
    logical :: trace_taylor
    logical :: has_state_work
    integer :: trace_call_id

    n = size(dg_frag%coef, 1)
    if (dg_frag%use_plane_wave_basis .or. allocated(dg_frag%coef_pw)) then
      stop "DG Taylor4-PC supports the pure fragment block-sparse route only"
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
    nsub_taylor = max(taylor4pc_substeps_default, ceiling(abs(dt) / taylor4pc_dt_substep_target))
    call initialize_taylor_runtime_options()

    if (.not. has_state_work) then
      if (yn_fix_func == 'n') then
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_mid, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy, .false.)
        if (dg_frag%coef_state_block_mode) call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_mid, .false.)
      end if
      return
    end if

    call ensure_state_arrays(n, nstate_prop, dg_frag%nspin)
    coef_base(:, :, :) = dg_frag%coef(:, 1:nstate_prop, :)

    target_bytes = int(cfg_state_work_target_mb, kind=8) * 1024_8 * 1024_8
    bytes_per_state = 16_8 * int(max(1, n), kind=8) * int(max(1, dg_frag%nspin), kind=8) * &
                      int(state_work_vectors, kind=8)
    nstate_work = max(1, min(nstate_prop, int(max(1_8, target_bytes / max(1_8, bytes_per_state)))))
    call ensure_block_arrays(n, nstate_work, dg_frag%nspin)
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

    dg_frag%coef(:, 1:nstate_prop, :) = coef_base(:, :, :)
    if (yn_fix_func == 'n') then
      call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_mid, &
                                            lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                            rho, rho_s, Vh, Vxc, Vpsl, energy, .false.)
      if (dg_frag%coef_state_block_mode) call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_mid, .false.)
    end if
    t_taylor0 = get_wtime()
    call apply_taylor4_with_current_hamiltonian(coef_base, Ac_mid, coef_next, nstate_prop, nstate_work, &
                                                state_global_first)
    t_taylor1 = get_wtime()
    time_taylor_apply = t_taylor1 - t_taylor0
    dg_frag%coef(:, 1:nstate_prop, :) = coef_next(:, :, :)
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

    subroutine apply_taylor4_with_current_hamiltonian(coef_in, Ac_use, coef_out, nstate_current, nstate_block, &
                                                      state_global_offset)
      complex(8), intent(in) :: coef_in(:, :, :)
      real(8), intent(in) :: Ac_use(3)
      complex(8), intent(out) :: coef_out(:, :, :)
      integer, intent(in) :: nstate_current, nstate_block, state_global_offset
      integer :: isub

      work_in(:, :, :) = coef_in(:, :, :)
      do isub = 1, nsub_taylor
        call apply_taylor4_single_step(work_in, Ac_use, work_out, nstate_current, nstate_block, &
                                       dt / dble(nsub_taylor), state_global_offset)
        work_in(:, :, :) = work_out(:, :, :)
      end do
      coef_out(:, :, :) = work_in(:, :, :)
      dg_frag%coef(:, 1:nstate_current, :) = coef_out(:, 1:nstate_current, :)
    end subroutine apply_taylor4_with_current_hamiltonian

    subroutine apply_taylor4_single_step(coef_in, Ac_use, coef_out, nstate_current, nstate_block, dt_step, &
                                         state_global_offset)
      complex(8), intent(in) :: coef_in(:, :, :)
      real(8), intent(in) :: Ac_use(3)
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
          dg_frag%coef(:, state0:state0+nstate_blk-1, :) = term_blk(:, 1:nstate_blk, :)
          deriv_blk(:, 1:nstate_blk, :) = (0.0d0, 0.0d0)
          call calculate_time_derivative(dg_frag, system, mg, ppg, Ac_use, &
                                         deriv_blk(:, 1:nstate_blk, :), state_s, state_e)
          call check_finite_taylor_block(order, state_s, state_e, deriv_blk(:, 1:nstate_blk, :))
          coeff = coeff * dt_step / dble(order)
          result_blk(:, 1:nstate_blk, :) = result_blk(:, 1:nstate_blk, :) + &
                                           coeff * deriv_blk(:, 1:nstate_blk, :)
          term_blk(:, 1:nstate_blk, :) = deriv_blk(:, 1:nstate_blk, :)
        end do
        coef_out(:, state0:state0+nstate_blk-1, :) = result_blk(:, 1:nstate_blk, :)
      end do
    end subroutine apply_taylor4_single_step

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
