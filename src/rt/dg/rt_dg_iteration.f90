  subroutine tddft_dg_fragment_iteration(dg_frag, system, info, rt, itt, dt, &
                                         lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                         rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use salmon_global, only: yn_fix_func, theory, nelec
    use communication, only: comm_summation, comm_is_root
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional
    use timer, only: timer_begin, timer_end, LOG_CALC_TIME_PROPAGATION, LOG_CALC_RHO
    use rt_dg_fragment_ops, only: zero_nonowned_coefficients, apply_overlap_operator
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
    logical :: enable_iteration_trace, enable_puredg_occ_unitarity
    logical :: startup_pre_relax_enable
    logical :: startup_cg_basis_enable
    logical :: startup_occ_unitarity_enable
    logical :: startup_cg_basis_allowed, startup_cg_basis_residual_enable
    logical :: startup_pre_relax_current_enable
    logical :: startup_pre_relax_charge_enable
    logical :: pre_relax_step_accepted
    logical :: pre_relax_coef_converged, pre_relax_current_converged
    logical :: pre_relax_charge_converged
    character(len=64) :: env_trace, env_puredg_occ_unitarity
    character(len=64) :: env_startup_pre_relax_steps, env_startup_pre_relax_tol
    character(len=64) :: env_startup_pre_relax_dtau
    character(len=64) :: env_startup_pre_relax_current_tol
    character(len=64) :: env_startup_pre_relax_charge_tol
    character(len=64) :: env_startup_pre_relax_backtrack_max
    character(len=64) :: env_startup_pre_relax_backtrack_factor
    character(len=64) :: env_startup_cg_basis_once
    character(len=64) :: env_startup_cg_basis_max_cycles, env_startup_cg_basis_res_tol
    character(len=64) :: env_startup_occ_unitarity
    integer :: env_trace_len, env_trace_status
    integer :: startup_pre_relax_steps, startup_pre_relax_iter, parse_status
    integer :: startup_pre_relax_backtrack_max, startup_pre_relax_backtrack_iter
    integer :: startup_cg_basis_max_cycles, startup_cg_basis_cycle
    real(8) :: startup_pre_relax_tol, pre_relax_rel_change, pre_relax_norm_prev, pre_relax_norm_diff
    real(8) :: pre_relax_norm_prev_local, pre_relax_norm_prev_global
    real(8) :: pre_relax_norm_diff_local, pre_relax_norm_diff_global
    real(8) :: pre_relax_residual_local, pre_relax_residual_global, pre_relax_residual_step
    real(8) :: startup_pre_relax_dtau
    real(8) :: startup_pre_relax_current_tol, startup_pre_relax_charge_tol
    real(8) :: pre_relax_current_rms, pre_relax_charge_residual
    real(8) :: pre_relax_charge_residual_prev, pre_relax_charge_accept_limit
    real(8) :: startup_pre_relax_backtrack_factor, startup_pre_relax_dtau_trial
    real(8) :: Ac_zero(3)
    real(8) :: startup_cg_basis_res_tol, startup_cg_basis_residual
    real(8) :: startup_cg_basis_norm_local, startup_cg_basis_norm_global
    real(8) :: startup_cg_basis_res_local, startup_cg_basis_res_global
    complex(8), allocatable :: coef_prev(:,:,:)
    complex(8), allocatable :: coef_prev_pw(:,:,:)
    complex(8), allocatable :: pre_relax_dcoef(:,:,:), pre_relax_dcoef_pw(:,:,:)
    complex(8), allocatable :: startup_cg_basis_dcoef(:,:,:), startup_cg_basis_dcoef_pw(:,:,:)
    complex(8), parameter :: zi = (0.0d0, 1.0d0)
    logical :: basis_functions_changed
    logical, save :: startup_cg_basis_done = .false.
    logical, save :: startup_pre_relax_done = .false.

    enable_iteration_trace = .false.
    enable_puredg_occ_unitarity = .false.
    env_trace = ''
    env_puredg_occ_unitarity = ''
    env_startup_pre_relax_steps = ''
    env_startup_pre_relax_tol = ''
    env_startup_pre_relax_dtau = ''
    env_startup_pre_relax_current_tol = ''
    env_startup_pre_relax_charge_tol = ''
    env_startup_pre_relax_backtrack_max = ''
    env_startup_pre_relax_backtrack_factor = ''
    env_startup_cg_basis_once = ''
    env_startup_cg_basis_max_cycles = ''
    env_startup_cg_basis_res_tol = ''
    env_startup_occ_unitarity = ''
    startup_pre_relax_steps = 0
    startup_pre_relax_tol = 1.0d-8
    startup_pre_relax_dtau = 1.0d-4
    startup_pre_relax_current_tol = -1.0d0
    startup_pre_relax_charge_tol = -1.0d0
    startup_pre_relax_backtrack_max = 8
    startup_pre_relax_backtrack_factor = 0.5d0
    startup_pre_relax_current_enable = .false.
    startup_pre_relax_charge_enable = .false.
    startup_pre_relax_enable = .false.
    startup_cg_basis_enable = .false.
    startup_occ_unitarity_enable = .false.
    startup_cg_basis_allowed = .false.
    startup_cg_basis_residual_enable = .false.
    startup_cg_basis_max_cycles = 1
    startup_cg_basis_res_tol = -1.0d0
    call get_environment_variable("SALMON_DG_ITER_TRACE", env_trace, length=env_trace_len, status=env_trace_status)
    if (env_trace_status == 0 .and. env_trace_len > 0) then
      if (env_trace(1:1) == '1' .or. env_trace(1:1) == 'y' .or. env_trace(1:1) == 'Y' .or. &
          env_trace(1:1) == 't' .or. env_trace(1:1) == 'T') then
        enable_iteration_trace = .true.
      end if
    end if
    call get_environment_variable("SALMON_DG_PUREDG_OCC_UNITARY", env_puredg_occ_unitarity, &
                                  length=env_trace_len, status=env_trace_status)
    if (env_trace_status == 0 .and. env_trace_len > 0) then
      if (env_puredg_occ_unitarity(1:1) == '1' .or. env_puredg_occ_unitarity(1:1) == 'y' .or. &
          env_puredg_occ_unitarity(1:1) == 'Y' .or. env_puredg_occ_unitarity(1:1) == 't' .or. &
          env_puredg_occ_unitarity(1:1) == 'T') then
        enable_puredg_occ_unitarity = .true.
      end if
    end if
    call get_environment_variable("SALMON_DG_STARTUP_PRE_RELAX_STEPS", env_startup_pre_relax_steps, &
                                  length=env_trace_len, status=env_trace_status)
    if (env_trace_status == 0 .and. env_trace_len > 0) then
      read(env_startup_pre_relax_steps(1:env_trace_len), *, iostat=parse_status) startup_pre_relax_steps
      if (parse_status /= 0) startup_pre_relax_steps = 0
      startup_pre_relax_steps = max(0, startup_pre_relax_steps)
    end if
    call get_environment_variable("SALMON_DG_STARTUP_PRE_RELAX_TOL", env_startup_pre_relax_tol, &
                                  length=env_trace_len, status=env_trace_status)
    if (env_trace_status == 0 .and. env_trace_len > 0) then
      read(env_startup_pre_relax_tol(1:env_trace_len), *, iostat=parse_status) startup_pre_relax_tol
      if (parse_status /= 0 .or. startup_pre_relax_tol <= 0.0d0) startup_pre_relax_tol = 1.0d-8
    end if
    call get_environment_variable("SALMON_DG_STARTUP_PRE_RELAX_DTAU", env_startup_pre_relax_dtau, &
                                  length=env_trace_len, status=env_trace_status)
    if (env_trace_status == 0 .and. env_trace_len > 0) then
      read(env_startup_pre_relax_dtau(1:env_trace_len), *, iostat=parse_status) startup_pre_relax_dtau
      if (parse_status /= 0 .or. startup_pre_relax_dtau <= 0.0d0) startup_pre_relax_dtau = 1.0d-4
    end if
    call get_environment_variable("SALMON_DG_STARTUP_PRE_RELAX_CURRENT_TOL", env_startup_pre_relax_current_tol, &
                                  length=env_trace_len, status=env_trace_status)
    if (env_trace_status == 0 .and. env_trace_len > 0) then
      read(env_startup_pre_relax_current_tol(1:env_trace_len), *, iostat=parse_status) startup_pre_relax_current_tol
      if (parse_status /= 0) startup_pre_relax_current_tol = -1.0d0
    end if
    call get_environment_variable("SALMON_DG_STARTUP_PRE_RELAX_CHARGE_TOL", env_startup_pre_relax_charge_tol, &
                                  length=env_trace_len, status=env_trace_status)
    if (env_trace_status == 0 .and. env_trace_len > 0) then
      read(env_startup_pre_relax_charge_tol(1:env_trace_len), *, iostat=parse_status) startup_pre_relax_charge_tol
      if (parse_status /= 0) startup_pre_relax_charge_tol = -1.0d0
    end if
    call get_environment_variable("SALMON_DG_STARTUP_PRE_RELAX_BACKTRACK_MAX", env_startup_pre_relax_backtrack_max, &
                                  length=env_trace_len, status=env_trace_status)
    if (env_trace_status == 0 .and. env_trace_len > 0) then
      read(env_startup_pre_relax_backtrack_max(1:env_trace_len), *, iostat=parse_status) startup_pre_relax_backtrack_max
      if (parse_status /= 0) startup_pre_relax_backtrack_max = 8
      startup_pre_relax_backtrack_max = max(1, startup_pre_relax_backtrack_max)
    end if
    call get_environment_variable("SALMON_DG_STARTUP_PRE_RELAX_BACKTRACK_FACTOR", env_startup_pre_relax_backtrack_factor, &
                                  length=env_trace_len, status=env_trace_status)
    if (env_trace_status == 0 .and. env_trace_len > 0) then
      read(env_startup_pre_relax_backtrack_factor(1:env_trace_len), *, iostat=parse_status) startup_pre_relax_backtrack_factor
      if (parse_status /= 0 .or. startup_pre_relax_backtrack_factor <= 0.0d0 .or. startup_pre_relax_backtrack_factor >= 1.0d0) then
        startup_pre_relax_backtrack_factor = 0.5d0
      end if
    end if
    startup_pre_relax_current_enable = (startup_pre_relax_current_tol > 0.0d0)
    startup_pre_relax_charge_enable = (startup_pre_relax_charge_tol > 0.0d0)
    startup_pre_relax_enable = (startup_pre_relax_steps > 0)
    call get_environment_variable("SALMON_DG_STARTUP_CG_BASIS_ONCE", env_startup_cg_basis_once, &
                                  length=env_trace_len, status=env_trace_status)
    if (env_trace_status == 0 .and. env_trace_len > 0) then
      if (env_startup_cg_basis_once(1:1) == '1' .or. env_startup_cg_basis_once(1:1) == 'y' .or. &
          env_startup_cg_basis_once(1:1) == 'Y' .or. env_startup_cg_basis_once(1:1) == 't' .or. &
          env_startup_cg_basis_once(1:1) == 'T') then
        startup_cg_basis_enable = .true.
      end if
    end if
    call get_environment_variable("SALMON_DG_STARTUP_CG_BASIS_MAX_CYCLES", env_startup_cg_basis_max_cycles, &
                                  length=env_trace_len, status=env_trace_status)
    if (env_trace_status == 0 .and. env_trace_len > 0) then
      read(env_startup_cg_basis_max_cycles(1:env_trace_len), *, iostat=parse_status) startup_cg_basis_max_cycles
      if (parse_status /= 0) startup_cg_basis_max_cycles = 1
      startup_cg_basis_max_cycles = max(1, startup_cg_basis_max_cycles)
    end if
    call get_environment_variable("SALMON_DG_STARTUP_CG_BASIS_RES_TOL", env_startup_cg_basis_res_tol, &
                                  length=env_trace_len, status=env_trace_status)
    if (env_trace_status == 0 .and. env_trace_len > 0) then
      read(env_startup_cg_basis_res_tol(1:env_trace_len), *, iostat=parse_status) startup_cg_basis_res_tol
      if (parse_status /= 0) startup_cg_basis_res_tol = -1.0d0
    end if
    startup_cg_basis_residual_enable = (startup_cg_basis_res_tol > 0.0d0)
    call get_environment_variable("SALMON_DG_STARTUP_OCC_UNITARY", env_startup_occ_unitarity, &
                                  length=env_trace_len, status=env_trace_status)
    if (env_trace_status == 0 .and. env_trace_len > 0) then
      if (env_startup_occ_unitarity(1:1) == '1' .or. env_startup_occ_unitarity(1:1) == 'y' .or. &
          env_startup_occ_unitarity(1:1) == 'Y' .or. env_startup_occ_unitarity(1:1) == 't' .or. &
          env_startup_occ_unitarity(1:1) == 'T') then
        startup_occ_unitarity_enable = .true.
      end if
    end if

    if (itt == 1 .and. startup_occ_unitarity_enable) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') "[INFO] DG startup occupied-subspace unitarity begin"
      end if
      call stabilize_coeff_unitarity(dg_frag, 0, occupied_only=.true.)
      if (allocated(dg_frag%coef_new)) dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
      if (trim(yn_fix_func) == 'n') then
        Ac_zero(:) = 0.0d0
        call update_density_and_hamiltonian(dg_frag, system, info, rt, itt, Ac_zero, &
                                            lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                            rho, rho_s, Vh, Vxc, Vpsl, energy)
      end if
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') "[INFO] DG startup occupied-subspace unitarity end"
      end if
    end if

    if (itt == 1 .and. startup_cg_basis_enable .and. .not. startup_cg_basis_done) then
      Ac_zero(:) = 0.0d0
      startup_cg_basis_allowed = (info%npk == 1 .and. info%nporbital == 1 .and. product(info%nprgrid(1:3)) == 1)
      startup_cg_basis_allowed = startup_cg_basis_allowed .and. all(mg%num(1:3) == dg_frag%lgnum_total(1:3))
      if (.not. startup_cg_basis_allowed) then
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a)') "[INFO] DG startup CG basis update skipped: MPI/layout precondition not met"
          write(*,'(1x,a,2(i0,1x),a,i0)') "       parent npk/nporbital=", info%npk, info%nporbital, &
                                          " product(nprgrid)=", product(info%nprgrid(1:3))
        end if
      else
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a,1pe12.4)') "[INFO] DG startup CG basis update begin: max_cycles=", startup_cg_basis_max_cycles, &
                                         " res_tol=", startup_cg_basis_res_tol
        end if
        if (trim(yn_fix_func) == 'n') then
          call update_density_and_hamiltonian(dg_frag, system, info, rt, itt, Ac_zero, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
        end if
        do startup_cg_basis_cycle = 1, startup_cg_basis_max_cycles
          call update_fragment_basis_via_cg(dg_frag, system, info, mg, stencil, srg, ppg, Vh, Vxc, Vpsl, basis_functions_changed)
          if (.not. basis_functions_changed) exit
          if (allocated(dg_frag%coef_new)) dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
          call zero_nonowned_coefficients(dg_frag)

          if (startup_cg_basis_residual_enable) then
            if (.not. allocated(startup_cg_basis_dcoef)) then
              allocate(startup_cg_basis_dcoef(size(dg_frag%coef,1), size(dg_frag%coef,2), size(dg_frag%coef,3)))
            end if
            if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) then
              if (.not. allocated(startup_cg_basis_dcoef_pw)) then
                allocate(startup_cg_basis_dcoef_pw(size(dg_frag%coef_pw,1), size(dg_frag%coef_pw,2), size(dg_frag%coef_pw,3)))
              end if
              call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_zero, 0, startup_cg_basis_dcoef, startup_cg_basis_dcoef_pw)
            else
              call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_zero, 0, startup_cg_basis_dcoef)
            end if

            startup_cg_basis_norm_local = sum(abs(dg_frag%coef(:, :, :))**2)
            startup_cg_basis_res_local = sum(abs(startup_cg_basis_dcoef(:, :, :))**2)
            if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw) .and. allocated(startup_cg_basis_dcoef_pw)) then
              startup_cg_basis_norm_local = startup_cg_basis_norm_local + sum(abs(dg_frag%coef_pw(:, :, :))**2)
              startup_cg_basis_res_local = startup_cg_basis_res_local + sum(abs(startup_cg_basis_dcoef_pw(:, :, :))**2)
            end if
            call comm_summation(startup_cg_basis_norm_local, startup_cg_basis_norm_global, dg_frag%icomm)
            call comm_summation(startup_cg_basis_res_local, startup_cg_basis_res_global, dg_frag%icomm)
            startup_cg_basis_residual = sqrt(max(startup_cg_basis_res_global, 0.0d0)) / &
                                        max(sqrt(max(startup_cg_basis_norm_global, 0.0d0)), 1.0d-30)

            if (comm_is_root(dg_frag%id)) then
              write(*,'(1x,a,i0,a,1pe12.4)') "[INFO] DG startup CG cycle=", startup_cg_basis_cycle, &
                                             " residual=", startup_cg_basis_residual
            end if
            if (startup_cg_basis_residual <= startup_cg_basis_res_tol) exit
          else
            if (comm_is_root(dg_frag%id)) then
              write(*,'(1x,a,i0)') "[INFO] DG startup CG cycle completed (residual check disabled): ", startup_cg_basis_cycle
            end if
          end if

          if (trim(yn_fix_func) == 'n') then
            call update_density_and_hamiltonian(dg_frag, system, info, rt, itt, Ac_zero, &
                                                lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                                rho, rho_s, Vh, Vxc, Vpsl, energy)
          end if
        end do
        if (allocated(startup_cg_basis_dcoef)) deallocate(startup_cg_basis_dcoef)
        if (allocated(startup_cg_basis_dcoef_pw)) deallocate(startup_cg_basis_dcoef_pw)
        if (trim(yn_fix_func) == 'n') then
          call update_density_and_hamiltonian(dg_frag, system, info, rt, itt, Ac_zero, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
        end if
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a)') "[INFO] DG startup CG basis update end"
        end if
      end if
      startup_cg_basis_done = .true.
    end if

    if (itt == 1 .and. startup_pre_relax_enable .and. .not. startup_pre_relax_done) then
      if (trim(yn_fix_func) /= 'n') then
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a)') "[INFO] DG startup pre-relax skipped because yn_fix_func /= 'n'"
        end if
      else
        if (.not. allocated(coef_prev)) allocate(coef_prev(size(dg_frag%coef,1), size(dg_frag%coef,2), size(dg_frag%coef,3)))
        if (.not. allocated(pre_relax_dcoef)) allocate(pre_relax_dcoef(size(dg_frag%coef,1), size(dg_frag%coef,2), size(dg_frag%coef,3)))
        if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) then
          if (.not. allocated(coef_prev_pw)) then
            allocate(coef_prev_pw(size(dg_frag%coef_pw,1), size(dg_frag%coef_pw,2), size(dg_frag%coef_pw,3)))
          end if
          if (.not. allocated(pre_relax_dcoef_pw)) then
            allocate(pre_relax_dcoef_pw(size(dg_frag%coef_pw,1), size(dg_frag%coef_pw,2), size(dg_frag%coef_pw,3)))
          end if
        end if
        Ac_zero(:) = 0.0d0
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,i0,a,f8.4)') "[INFO] DG startup pre-relax begin: steps=", startup_pre_relax_steps, &
                                         " tol=", startup_pre_relax_tol, " dtau=", startup_pre_relax_dtau, &
                                         " current_tol=", startup_pre_relax_current_tol, " charge_tol=", startup_pre_relax_charge_tol, &
                                         " backtrack_max=", startup_pre_relax_backtrack_max, " factor=", startup_pre_relax_backtrack_factor
        end if
        do startup_pre_relax_iter = 1, startup_pre_relax_steps
          coef_prev(:, :, :) = dg_frag%coef(:, :, :)
          if (allocated(coef_prev_pw)) coef_prev_pw(:, :, :) = dg_frag%coef_pw(:, :, :)
          call update_density_and_hamiltonian(dg_frag, system, info, rt, itt, Ac_zero, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
          pre_relax_charge_residual_prev = abs(dg_frag%elec_num_raw - dble(nelec)) / max(dble(nelec), 1.0d0)
          pre_relax_charge_accept_limit = pre_relax_charge_residual_prev
          if (startup_pre_relax_charge_enable) then
            pre_relax_charge_accept_limit = max(pre_relax_charge_accept_limit, startup_pre_relax_charge_tol)
          end if
          if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw) .and. allocated(pre_relax_dcoef_pw)) then
            call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_zero, 0, pre_relax_dcoef, pre_relax_dcoef_pw)
          else
            call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_zero, 0, pre_relax_dcoef)
          end if
          pre_relax_step_accepted = .false.
          startup_pre_relax_dtau_trial = startup_pre_relax_dtau
          do startup_pre_relax_backtrack_iter = 1, startup_pre_relax_backtrack_max
            dg_frag%coef(:, :, :) = coef_prev(:, :, :) - startup_pre_relax_dtau_trial * zi * pre_relax_dcoef(:, :, :)
            if (allocated(coef_prev_pw)) then
              dg_frag%coef_pw(:, :, :) = coef_prev_pw(:, :, :) - startup_pre_relax_dtau_trial * zi * pre_relax_dcoef_pw(:, :, :)
            end if
            call zero_nonowned_coefficients(dg_frag)
            call update_density_and_hamiltonian(dg_frag, system, info, rt, itt, Ac_zero, &
                                                lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                                rho, rho_s, Vh, Vxc, Vpsl, energy)
            pre_relax_norm_prev_local = sum(abs(coef_prev(:, :, :))**2)
            pre_relax_norm_diff_local = sum(abs(dg_frag%coef(:, :, :) - coef_prev(:, :, :))**2)
            pre_relax_residual_local = sum(abs(pre_relax_dcoef(:, :, :))**2)
            if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw) .and. allocated(pre_relax_dcoef_pw)) then
              pre_relax_norm_prev_local = pre_relax_norm_prev_local + sum(abs(coef_prev_pw(:, :, :))**2)
              pre_relax_residual_local = pre_relax_residual_local + sum(abs(pre_relax_dcoef_pw(:, :, :))**2)
            end if
            call comm_summation(pre_relax_norm_prev_local, pre_relax_norm_prev_global, dg_frag%icomm)
            call comm_summation(pre_relax_norm_diff_local, pre_relax_norm_diff_global, dg_frag%icomm)
            call comm_summation(pre_relax_residual_local, pre_relax_residual_global, dg_frag%icomm)
            pre_relax_norm_prev = sqrt(max(pre_relax_norm_prev_global, 0.0d0))
            pre_relax_norm_diff = sqrt(max(pre_relax_norm_diff_global, 0.0d0))
            pre_relax_rel_change = pre_relax_norm_diff / max(pre_relax_norm_prev, 1.0d-30)
            pre_relax_residual_step = startup_pre_relax_dtau_trial * sqrt(max(pre_relax_residual_global, 0.0d0)) / max(pre_relax_norm_prev, 1.0d-30)
            pre_relax_current_rms = 0.0d0
            if (startup_pre_relax_current_enable) then
              call calculate_observables(dg_frag, system, mg, stencil, ppg, rt, itt, Vh, Vxc, Vpsl, rho)
              pre_relax_current_rms = sqrt(sum(dg_frag%current(1:3)**2))
            end if
            pre_relax_charge_residual = abs(dg_frag%elec_num_raw - dble(nelec)) / max(dble(nelec), 1.0d0)
            if (startup_pre_relax_charge_enable .and. pre_relax_charge_residual > pre_relax_charge_accept_limit) then
              if (startup_pre_relax_backtrack_iter < startup_pre_relax_backtrack_max) then
                if (comm_is_root(dg_frag%id)) then
                  write(*,'(1x,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') "[INFO] DG startup pre-relax backtrack iter=", startup_pre_relax_iter, &
                                                 " trial=", startup_pre_relax_backtrack_iter, " dtau=", startup_pre_relax_dtau_trial, &
                                                 " charge_res=", pre_relax_charge_residual, " prev_charge_res=", pre_relax_charge_residual_prev, &
                                                 " accept_limit=", pre_relax_charge_accept_limit
                end if
                dg_frag%coef(:, :, :) = coef_prev(:, :, :)
                if (allocated(coef_prev_pw)) dg_frag%coef_pw(:, :, :) = coef_prev_pw(:, :, :)
                call zero_nonowned_coefficients(dg_frag)
                startup_pre_relax_dtau_trial = startup_pre_relax_dtau_trial * startup_pre_relax_backtrack_factor
                cycle
              else
                dg_frag%coef(:, :, :) = coef_prev(:, :, :)
                if (allocated(coef_prev_pw)) dg_frag%coef_pw(:, :, :) = coef_prev_pw(:, :, :)
                call zero_nonowned_coefficients(dg_frag)
                call update_density_and_hamiltonian(dg_frag, system, info, rt, itt, Ac_zero, &
                                                    lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                                    rho, rho_s, Vh, Vxc, Vpsl, energy)
                if (comm_is_root(dg_frag%id)) then
                  write(*,'(1x,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a)') "[INFO] DG startup pre-relax rejected at iter=", startup_pre_relax_iter, &
                                                 " charge_res=", pre_relax_charge_residual, &
                                                 " prev_charge_res=", pre_relax_charge_residual_prev, &
                                                 " accept_limit=", pre_relax_charge_accept_limit, " after max backtracks"
                end if
                exit
              end if
            end if
            pre_relax_step_accepted = .true.
            exit
          end do
          if (.not. pre_relax_step_accepted) exit
          if (comm_is_root(dg_frag%id)) then
            write(*,'(1x,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') "[INFO] DG startup pre-relax iter=", startup_pre_relax_iter, &
                                           " rel_change=", pre_relax_rel_change, " residual_step=", pre_relax_residual_step, &
                                           " current_rms=", pre_relax_current_rms, " charge_res=", pre_relax_charge_residual, &
                                           " dtau=", startup_pre_relax_dtau_trial
          end if
          pre_relax_coef_converged = (max(pre_relax_rel_change, pre_relax_residual_step) <= startup_pre_relax_tol)
          pre_relax_current_converged = (.not. startup_pre_relax_current_enable) .or. &
                                        (pre_relax_current_rms <= startup_pre_relax_current_tol)
          pre_relax_charge_converged = (.not. startup_pre_relax_charge_enable) .or. &
                                       (pre_relax_charge_residual <= startup_pre_relax_charge_tol)
          if (pre_relax_coef_converged .and. pre_relax_current_converged .and. pre_relax_charge_converged) exit
        end do
        if (allocated(dg_frag%coef_new)) dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
        deallocate(coef_prev)
        if (allocated(coef_prev_pw)) deallocate(coef_prev_pw)
        deallocate(pre_relax_dcoef)
        if (allocated(pre_relax_dcoef_pw)) deallocate(pre_relax_dcoef_pw)
      end if
      startup_pre_relax_done = .true.
    end if

    ! Time evolution in fragment basis coefficient space
    if (enable_iteration_trace) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        iteration trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
        " stage=", "entry"
      flush(6)
    end if
    dg_frag%current_iteration = itt
    call debug_coef_metric("entry")
    call timer_begin(LOG_CALC_TIME_PROPAGATION)
    select case(dg_frag%time_integrator)
    case(1, 3)  ! SSPRK3 or RK4
      call time_evolution_rk(dg_frag, system, info, rt, itt, dt, &
                             lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                             rho, rho_s, Vh, Vxc, Vpsl, energy)
    case(2)  ! AETRS
      call time_evolution_aetrs(dg_frag, system, mg, stencil, ppg, rt, itt, dt)
    case default
      stop "Unknown time integrator for DG-Fragment method"
    end select
    call timer_end(LOG_CALC_TIME_PROPAGATION)
    if (enable_iteration_trace) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        iteration trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
        " stage=", "after-time-evolution"
      flush(6)
    end if
    call debug_coef_metric("after-time-evolution")

    ! Keep coefficient normalization after each propagation step for mixed DG+PW runs.
    ! Pure-DG path keeps the legacy flow to avoid unnecessary overlap-side effects.
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw) .and. dg_frag%n_plane_waves > 0) then
      call stabilize_coeff_unitarity(dg_frag, itt)
    else if (enable_puredg_occ_unitarity) then
      call stabilize_coeff_unitarity(dg_frag, itt, occupied_only=.true.)
    end if
    call debug_coef_metric("after-unitarity")

    ! Self-consistent update of density and Hamiltonian (if enabled)
    ! Performance note:
    !   - Coefficient space evolution: O(n_basis²) - very fast
    !   - Density reconstruction: O(n_frag × n_basis² × n_occ × n_grid) - expensive!
    !   - Hartree/XC: O(N_grid log N_grid) - moderate (same as standard RT-TDDFT)
    ! Recommendation:
    !   - Linear response (weak field): yn_fix_func='y' (no update, fast)
    !   - Nonlinear response (strong field): yn_fix_func='n' (self-consistent, slower)
    !   - Future optimization: update every N steps instead of every step
    if (yn_fix_func == 'n') then
      if (dg_frag%time_integrator /= 3 .or. dg_frag%yn_adaptive_basis) then
        ! For adaptive-basis mode, keep post-step update active even in RK4
        ! so basis-update detection/trigger logic runs.
        if (enable_iteration_trace) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        iteration stage: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, " stage=", "before-update-density-hmat"
          flush(6)
        end if
        call timer_begin(LOG_CALC_RHO)
        call update_density_and_hamiltonian(dg_frag, system, info, rt, itt, rt%Ac_tot(:,itt), &
                                            lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                            rho, rho_s, Vh, Vxc, Vpsl, energy)
        call timer_end(LOG_CALC_RHO)
        if (enable_iteration_trace) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        iteration stage: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, " stage=", "after-update-density-hmat"
          flush(6)
        end if
      end if
    end if
    
    ! Calculate observables
    if (enable_iteration_trace) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        iteration stage: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, " stage=", "before-calc-observables"
      flush(6)
    end if
    call calculate_observables(dg_frag, system, mg, stencil, ppg, rt, itt, Vh, Vxc, Vpsl, rho)
    if (enable_iteration_trace) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        iteration stage: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, " stage=", "after-calc-observables"
      flush(6)
    end if

    if (trim(theory) == 'single_scale_maxwell_tddft') then
      call calculate_microscopic_current_dg(dg_frag, system, mg, stencil, rt%j_e)
    end if
    
  contains

  subroutine debug_coef_metric(stage_label)
    character(len=*), intent(in) :: stage_label

    integer :: ispin_probe, n_frag_rows, n_pw, n_tot, nstate_probe, io
    integer :: env_probe_len, env_probe_status, nocc_probe
    integer :: coef_metric_maxitt, parse_status
    complex(8), allocatable :: vec(:), svec(:)
    real(8) :: c2(3), cs2(3), cs2_occ_sum, cs2_occ_min, cs2_occ_max
    real(8) :: cs2_occ_rocc_sum, cs2_state, occ_factor, cs2_occ_rocc_sum_global
    logical :: enable_coef_metric_probe
    character(len=64) :: env_probe

    enable_coef_metric_probe = .false.
    env_probe = ''
    call get_environment_variable("SALMON_DG_COEF_METRIC_PROBE", env_probe, length=env_probe_len, status=env_probe_status)
    if (env_probe_status == 0 .and. env_probe_len > 0) then
      if (env_probe(1:1) == '1' .or. env_probe(1:1) == 'y' .or. env_probe(1:1) == 'Y' .or. &
          env_probe(1:1) == 't' .or. env_probe(1:1) == 'T') then
        enable_coef_metric_probe = .true.
      end if
    end if
    if (.not. enable_coef_metric_probe) return
    coef_metric_maxitt = 1
    env_probe = ''
    call get_environment_variable("SALMON_DG_COEF_METRIC_MAXITT", env_probe, length=env_probe_len, status=env_probe_status)
    if (env_probe_status == 0 .and. env_probe_len > 0) then
      read(env_probe(1:env_probe_len), *, iostat=parse_status) coef_metric_maxitt
      if (parse_status /= 0) coef_metric_maxitt = 1
    end if
    if (coef_metric_maxitt < 1) coef_metric_maxitt = 1
    if (itt > coef_metric_maxitt) return
    if (dg_frag%nspin <= 0 .or. dg_frag%nstate_tot <= 0) return

    ispin_probe = 1
    n_frag_rows = dg_frag%n_mat_max
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n_frag_rows + n_pw
    nstate_probe = min(3, dg_frag%nstate_tot)
    if (n_tot <= 0 .or. nstate_probe <= 0) return

    allocate(vec(n_tot), svec(n_tot))

    c2(:) = 0.0d0
    cs2(:) = 0.0d0
    cs2_occ_sum = 0.0d0
    cs2_occ_rocc_sum = 0.0d0
    cs2_occ_min = huge(1.0d0)
    cs2_occ_max = -huge(1.0d0)
    nocc_probe = min(dg_frag%nstate_tot, dg_frag%nocc_spin(ispin_probe))
    if (n_pw > 0) then
      do io = 1, nstate_probe
        vec(:) = (0.0d0, 0.0d0)
        if (n_frag_rows > 0) vec(1:n_frag_rows) = dg_frag%coef(1:n_frag_rows, io, ispin_probe)
        vec(n_frag_rows+1:n_tot) = dg_frag%coef_pw(1:n_pw, io, ispin_probe)
        call apply_overlap_operator(dg_frag, ispin_probe, vec, svec, .true.)
        c2(io) = real(sum(conjg(vec) * vec), kind=8)
        cs2(io) = real(sum(conjg(vec) * svec), kind=8)
      end do
      do io = 1, nocc_probe
        vec(:) = (0.0d0, 0.0d0)
        if (n_frag_rows > 0) vec(1:n_frag_rows) = dg_frag%coef(1:n_frag_rows, io, ispin_probe)
        vec(n_frag_rows+1:n_tot) = dg_frag%coef_pw(1:n_pw, io, ispin_probe)
        call apply_overlap_operator(dg_frag, ispin_probe, vec, svec, .true.)
        cs2_state = real(sum(conjg(vec) * svec), kind=8)
        cs2_occ_sum = cs2_occ_sum + cs2_state
        cs2_occ_min = min(cs2_occ_min, cs2_state)
        cs2_occ_max = max(cs2_occ_max, cs2_state)
        occ_factor = 1.0d0
        if (allocated(system%rocc)) then
          if (io <= size(system%rocc, 1) .and. 1 <= size(system%rocc, 3)) then
            occ_factor = max(0.0d0, system%rocc(io, 1, 1))
          end if
        end if
        cs2_occ_rocc_sum = cs2_occ_rocc_sum + occ_factor * cs2_state
      end do
    else
      do io = 1, nstate_probe
        vec(:) = (0.0d0, 0.0d0)
        if (n_frag_rows > 0) vec(1:n_frag_rows) = dg_frag%coef(1:n_frag_rows, io, ispin_probe)
        call apply_overlap_operator(dg_frag, ispin_probe, vec, svec, .true.)
        c2(io) = real(sum(conjg(vec) * vec), kind=8)
        cs2(io) = real(sum(conjg(vec) * svec), kind=8)
      end do
      do io = 1, nocc_probe
        vec(:) = (0.0d0, 0.0d0)
        if (n_frag_rows > 0) vec(1:n_frag_rows) = dg_frag%coef(1:n_frag_rows, io, ispin_probe)
        call apply_overlap_operator(dg_frag, ispin_probe, vec, svec, .true.)
        cs2_state = real(sum(conjg(vec) * svec), kind=8)
        cs2_occ_sum = cs2_occ_sum + cs2_state
        cs2_occ_min = min(cs2_occ_min, cs2_state)
        cs2_occ_max = max(cs2_occ_max, cs2_state)
        occ_factor = 1.0d0
        if (allocated(system%rocc)) then
          if (io <= size(system%rocc, 1) .and. 1 <= size(system%rocc, 3)) then
            occ_factor = max(0.0d0, system%rocc(io, 1, 1))
          end if
        end if
        cs2_occ_rocc_sum = cs2_occ_rocc_sum + occ_factor * cs2_state
      end do
    end if

    call comm_summation(cs2_occ_rocc_sum, cs2_occ_rocc_sum_global, dg_frag%icomm)

    if (dg_frag%id == 0) then
      write(*,'(1x,a,a,a,3(1x,es12.4),a,3(1x,es12.4))') "        coef stage probe: stage=", trim(stage_label), &
        " c2=", c2(1), c2(2), c2(3), " cs2=", cs2(1), cs2(2), cs2(3)
      flush(6)
      write(*,'(1x,a,a,a,i0,a,3(1x,es12.4))') "        coef occ probe: stage=", trim(stage_label), &
        " nocc=", nocc_probe, " cs2(sum,min,max)=", cs2_occ_sum, cs2_occ_min, cs2_occ_max
      flush(6)
      write(*,'(1x,a,a,a,i0,a,4(1x,es12.4))') "        coef-density comparable: stage=", trim(stage_label), &
        " itt=", itt, " cs2_rocc_sum / (cs2_rocc_sum-nelec) / Ne_raw / (Ne_raw-nelec)=", &
        cs2_occ_rocc_sum_global, cs2_occ_rocc_sum_global - dble(nelec), dg_frag%elec_num_raw, dg_frag%elec_num_raw - dble(nelec)
      flush(6)
      write(*,'(1x,a,a,a,i0,a,3(1x,es12.4))') "        coef-density probe: stage=", trim(stage_label), &
        " itt=", itt, " Ne_raw/rho_scale/raw-target=", dg_frag%elec_num_raw, dg_frag%rho_scale_factor, &
        dg_frag%elec_num_raw - dble(nelec)
      flush(6)
    end if

    deallocate(vec, svec)
  end subroutine debug_coef_metric

  end subroutine tddft_dg_fragment_iteration
