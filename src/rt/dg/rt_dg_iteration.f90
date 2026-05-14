  subroutine tddft_dg_fragment_iteration(dg_frag, system, info, rt, itt, dt, &
                                         lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                         rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use communication, only: comm_is_root
    use salmon_global, only: yn_fix_func, theory
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional
    use timer, only: timer_begin, timer_end, LOG_CALC_TIME_PROPAGATION, LOG_CALC_RHO
    use rt_dg_fragment_ops, only: zero_nonowned_coefficients, &
                    reorthonormalize_mixed_occupied_subspace
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
    logical, save :: startup_pre_relax_done = .false.
    logical, save :: reorth_mixed_occ_initialized = .false.
    logical, save :: enable_reorth_mixed_occ = .false.
    integer :: pre_relax_steps, pre_relax_iter, env_len, env_status, parse_status
    real(8) :: pre_relax_tol, pre_relax_dtau, rel_change, norm_prev, norm_diff
    real(8) :: Ac_zero(3)
    character(len=64) :: env_val
    complex(8), allocatable :: dcoef_dt(:,:,:), coef_prev(:,:,:)
    complex(8), allocatable :: dcoef_dt_pw(:,:,:), coef_prev_pw(:,:,:)
    logical, save :: enable_eigenstate_check = .false.
    logical, save :: eigenstate_check_initialized = .false.
    logical, save :: occvirt_ref_mode_initialized = .false.
    logical, save :: enable_occvirt_ref_diag = .false.
    logical, save :: occvirt_ref_legacy_mode = .false.
    real(8) :: rayleigh_max_rel_dev
    integer :: rayleigh_nonev_count, io_rayleigh
    integer :: n_frag_ref, n_pw_ref, n_tot_ref, nstate_ref
    ! Time evolution in fragment basis coefficient space

    ! Initialize eigenstate check on iter 1
    if (itt == 1 .and. .not. eigenstate_check_initialized) then
      eigenstate_check_initialized = .true.
      env_val = ''
      call get_environment_variable('SALMON_DG_INITIAL_EIGENSTATE_CHECK', env_val, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        if (trim(env_val) == '1' .or. trim(env_val) == 'yes' .or. trim(env_val) == 'true') then
          enable_eigenstate_check = .true.
          if (comm_is_root(dg_frag%id)) then
            write(*,'(1x,a)') "[DG-CONFIG] SALMON_DG_INITIAL_EIGENSTATE_CHECK=ON"
            flush(6)
          end if
        end if
      end if
    end if

    if (itt == 1 .and. .not. startup_pre_relax_done) then
      pre_relax_steps = 0
      pre_relax_tol = 1.0d-10
      pre_relax_dtau = 1.0d-5

      env_val = ''
      call get_environment_variable('SALMON_DG_STARTUP_PRE_RELAX_STEPS', env_val, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_val(1:env_len), *, iostat=parse_status) pre_relax_steps
        if (parse_status /= 0) pre_relax_steps = 0
      end if

      env_val = ''
      call get_environment_variable('SALMON_DG_STARTUP_PRE_RELAX_TOL', env_val, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_val(1:env_len), *, iostat=parse_status) pre_relax_tol
        if (parse_status /= 0) pre_relax_tol = 1.0d-10
      end if

      env_val = ''
      call get_environment_variable('SALMON_DG_STARTUP_PRE_RELAX_DTAU', env_val, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_val(1:env_len), *, iostat=parse_status) pre_relax_dtau
        if (parse_status /= 0) pre_relax_dtau = 1.0d-5
      end if

      if (pre_relax_steps > 0 .and. yn_fix_func == 'n') then
        Ac_zero(:) = 0.0d0
        allocate(coef_prev(size(dg_frag%coef,1), size(dg_frag%coef,2), size(dg_frag%coef,3)))
        allocate(dcoef_dt(size(dg_frag%coef,1), size(dg_frag%coef,2), size(dg_frag%coef,3)))
        if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) then
          allocate(coef_prev_pw(size(dg_frag%coef_pw,1), size(dg_frag%coef_pw,2), size(dg_frag%coef_pw,3)))
          allocate(dcoef_dt_pw(size(dg_frag%coef_pw,1), size(dg_frag%coef_pw,2), size(dg_frag%coef_pw,3)))
        end if

        if (dg_frag%id == 0) then
          write(*,'(1x,a,i0,a,es12.4,a,es12.4)') 'startup pre-relax: steps=', pre_relax_steps, &
            ' tol=', pre_relax_tol, ' dtau=', pre_relax_dtau
          flush(6)
        end if

        do pre_relax_iter = 1, pre_relax_steps
          coef_prev = dg_frag%coef
          if (allocated(coef_prev_pw)) coef_prev_pw = dg_frag%coef_pw

          call update_density_and_hamiltonian(dg_frag, system, info, rt, itt, Ac_zero, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)

          if (allocated(dcoef_dt_pw)) then
            call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_zero, itt, dcoef_dt, dcoef_dt_pw)
          else
            call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_zero, itt, dcoef_dt)
          end if

          dg_frag%coef = dg_frag%coef - pre_relax_dtau * (0.0d0, 1.0d0) * dcoef_dt
          if (allocated(dcoef_dt_pw)) dg_frag%coef_pw = dg_frag%coef_pw - pre_relax_dtau * (0.0d0, 1.0d0) * dcoef_dt_pw

          call zero_nonowned_coefficients(dg_frag)

          norm_prev = sum(abs(coef_prev)**2)
          norm_diff = sum(abs(dg_frag%coef - coef_prev)**2)
          if (allocated(coef_prev_pw)) then
            norm_prev = norm_prev + sum(abs(coef_prev_pw)**2)
            norm_diff = norm_diff + sum(abs(dg_frag%coef_pw - coef_prev_pw)**2)
          end if
          rel_change = sqrt(norm_diff / max(norm_prev, 1.0d-30))

          if (dg_frag%id == 0) then
            write(*,'(1x,a,i0,a,es12.4)') 'startup pre-relax iter=', pre_relax_iter, ' rel_change=', rel_change
            flush(6)
          end if

          if (rel_change < pre_relax_tol) exit
        end do

        call update_density_and_hamiltonian(dg_frag, system, info, rt, itt, Ac_zero, &
                                            lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                            rho, rho_s, Vh, Vxc, Vpsl, energy)

        deallocate(coef_prev, dcoef_dt)
        if (allocated(coef_prev_pw)) deallocate(coef_prev_pw)
        if (allocated(dcoef_dt_pw)) deallocate(dcoef_dt_pw)
      end if

      startup_pre_relax_done = .true.
    end if

    if (itt == 1 .and. enable_eigenstate_check) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') "[DG-EIGENSTATE] Quick diagnostics: checking initial coef vs frozen H diagonal"
        if (allocated(dg_frag%H_mat) .and. allocated(dg_frag%esp)) then
          rayleigh_max_rel_dev = 0.0d0
          rayleigh_nonev_count = 0
          do io_rayleigh = 1, min(24, size(dg_frag%H_mat, 1), size(dg_frag%esp, 1))
            if (abs(dg_frag%esp(io_rayleigh, 1)) > 1.0d-20) then
              rayleigh_max_rel_dev = max(rayleigh_max_rel_dev, &
                abs(dg_frag%H_mat(io_rayleigh, io_rayleigh, 1) - dg_frag%esp(io_rayleigh, 1)) / abs(dg_frag%esp(io_rayleigh, 1)))
            end if
          end do
          write(*,'(1x,a,1pe12.4)') "[DG-EIGENSTATE] max|H_ii - esp_i|/|esp_i| (first 24 states)=", rayleigh_max_rel_dev
          write(*,'(1x,a)') "[DG-EIGENSTATE] if >> 1e-6: initial coef NOT eigenstate of frozen H,S"
        end if
      end if
    end if

    if (.not. occvirt_ref_mode_initialized) then
      env_val = ''
      call get_environment_variable('SALMON_DG_OCCVIRT_DIAG', env_val, length=env_len, status=env_status)
      enable_occvirt_ref_diag = .false.
      if (env_status == 0 .and. env_len > 0) then
        if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
            env_val(1:1) == 't' .or. env_val(1:1) == 'T') enable_occvirt_ref_diag = .true.
      end if

      env_val = ''
      call get_environment_variable('SALMON_DG_OCCVIRT_REF_MODE', env_val, length=env_len, status=env_status)
      occvirt_ref_legacy_mode = .false.
      if (env_status == 0 .and. env_len > 0) then
        enable_occvirt_ref_diag = .true.
        select case (env_val(1:1))
        case ('l','L','0')
          occvirt_ref_legacy_mode = .true.
        case default
          occvirt_ref_legacy_mode = .false.
        end select
      end if
      occvirt_ref_mode_initialized = .true.
      if (comm_is_root(dg_frag%id)) then
        if (.not. enable_occvirt_ref_diag) then
          write(*,'(1x,a)') '[DG-CONFIG] SALMON_DG_OCCVIRT_DIAG=OFF (default)'
        else if (occvirt_ref_legacy_mode) then
          write(*,'(1x,a)') '[DG-CONFIG] SALMON_DG_OCCVIRT_REF_MODE=legacy'
        else
          write(*,'(1x,a)') '[DG-CONFIG] SALMON_DG_OCCVIRT_REF_MODE=t0'
        end if
        flush(6)
      end if
    end if

    if (enable_occvirt_ref_diag .and. itt == 1 .and. .not. occvirt_ref_legacy_mode .and. .not. dg_frag%coef_ref_ready) then
      ! The occ/virt reference is a diagnostic snapshot in the local coefficient
      ! layout.  In orbital row-split mode coef(:,:,:) only stores local rows;
      ! using n_mat_max here would allocate a dense global matrix on every rank.
      n_frag_ref = size(dg_frag%coef, 1)
      n_pw_ref = 0
      if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw_ref = size(dg_frag%coef_pw, 1)
      n_tot_ref = n_frag_ref + n_pw_ref
      nstate_ref = dg_frag%nstate_tot
      if (n_tot_ref > 0 .and. nstate_ref > 0) then
        if (allocated(dg_frag%coef_ref_all)) then
          if (size(dg_frag%coef_ref_all, 1) /= n_tot_ref .or. size(dg_frag%coef_ref_all, 2) /= nstate_ref .or. &
              size(dg_frag%coef_ref_all, 3) /= dg_frag%nspin) then
            deallocate(dg_frag%coef_ref_all)
          end if
        end if
        if (.not. allocated(dg_frag%coef_ref_all)) allocate(dg_frag%coef_ref_all(n_tot_ref, nstate_ref, dg_frag%nspin))
        dg_frag%coef_ref_all(:, :, :) = (0.0d0, 0.0d0)
        dg_frag%coef_ref_all(1:n_frag_ref, 1:nstate_ref, :) = dg_frag%coef(1:n_frag_ref, 1:nstate_ref, :)
        if (n_pw_ref > 0) dg_frag%coef_ref_all(n_frag_ref+1:n_tot_ref, 1:nstate_ref, :) = dg_frag%coef_pw(1:n_pw_ref, 1:nstate_ref, :)
        dg_frag%coef_ref_ready = .true.
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a)') '[DG-OBS] occvirt reference fixed at t=0 (pre-propagation)'
          flush(6)
        end if
      end if
    end if

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

    if (.not. reorth_mixed_occ_initialized) then
      env_val = ''
      call get_environment_variable('SALMON_DG_REORTH_MIXED_OCC', env_val, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
            env_val(1:1) == 't' .or. env_val(1:1) == 'T') enable_reorth_mixed_occ = .true.
      end if
      reorth_mixed_occ_initialized = .true.
      if (enable_reorth_mixed_occ .and. comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') '[DG-CONFIG] SALMON_DG_REORTH_MIXED_OCC=ON'
        flush(6)
      end if
    end if

    if (enable_reorth_mixed_occ) then
      call reorthonormalize_mixed_occupied_subspace(dg_frag)
    end if
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
        call timer_begin(LOG_CALC_RHO)
        call update_density_and_hamiltonian(dg_frag, system, info, rt, itt, rt%Ac_tot(:,itt), &
                                            lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                            rho, rho_s, Vh, Vxc, Vpsl, energy)
        call timer_end(LOG_CALC_RHO)
      end if
    end if
    
    ! Calculate observables
    call calculate_observables(dg_frag, system, mg, stencil, ppg, rt, itt, Vh, Vxc, Vpsl, rho)

    if (trim(theory) == 'single_scale_maxwell_tddft') then
      call calculate_microscopic_current_dg(dg_frag, system, mg, stencil, rt%j_e)
    end if

  end subroutine tddft_dg_fragment_iteration
