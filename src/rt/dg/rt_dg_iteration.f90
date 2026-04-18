  subroutine tddft_dg_fragment_iteration(dg_frag, system, info, rt, itt, dt, &
                                         lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                         rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use salmon_global, only: yn_fix_func, theory, nelec
    use communication, only: comm_summation
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
    character(len=64) :: env_trace, env_puredg_occ_unitarity
    integer :: env_trace_len, env_trace_status

    enable_iteration_trace = .false.
    enable_puredg_occ_unitarity = .false.
    env_trace = ''
    env_puredg_occ_unitarity = ''
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
    call calculate_observables(dg_frag, system, mg, stencil, ppg, rt, itt, Vh, Vxc, Vpsl)
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
