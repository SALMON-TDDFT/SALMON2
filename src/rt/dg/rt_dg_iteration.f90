  subroutine tddft_dg_fragment_iteration(dg_frag, system, info, rt, itt, dt, &
                                         lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                         rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use salmon_global, only: yn_fix_func, theory
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
    logical, parameter :: enable_iteration_trace = .false.
    ! Time evolution in fragment basis coefficient space
    if (enable_iteration_trace) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        iteration trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
        " stage=", "entry"
      flush(6)
    end if
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
          write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        iteration stage: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " stage=", "before-update-density-hmat"
          flush(6)
        end if
        call timer_begin(LOG_CALC_RHO)
        call update_density_and_hamiltonian(dg_frag, system, info, rt, itt, rt%Ac_tot(:,itt), &
                                            lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                            rho, rho_s, Vh, Vxc, Vpsl, energy)
        call timer_end(LOG_CALC_RHO)
        if (enable_iteration_trace) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        iteration stage: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " stage=", "after-update-density-hmat"
          flush(6)
        end if
      end if
    end if
    
    ! Calculate observables
    call calculate_observables(dg_frag, system, mg, stencil, ppg, rt, itt, Vh, Vxc, Vpsl)

    if (trim(theory) == 'single_scale_maxwell_tddft') then
      call calculate_microscopic_current_dg(dg_frag, system, mg, stencil, rt%j_e)
    end if
    
  contains

  subroutine debug_coef_metric(stage_label)
    character(len=*), intent(in) :: stage_label

    integer :: ispin_probe, n_frag_rows, n_pw, n_tot, nstate_probe, io
    complex(8), allocatable :: vec(:), svec(:)
    real(8) :: c2(3), cs2(3)
    logical, parameter :: enable_coef_metric_probe = .false.

    if (.not. enable_coef_metric_probe) return
    if (itt /= 1) return
    if (dg_frag%id /= 0) return
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
    if (n_pw > 0) then
      do io = 1, nstate_probe
        vec(:) = (0.0d0, 0.0d0)
        if (n_frag_rows > 0) vec(1:n_frag_rows) = dg_frag%coef(1:n_frag_rows, io, ispin_probe)
        vec(n_frag_rows+1:n_tot) = dg_frag%coef_pw(1:n_pw, io, ispin_probe)
        call apply_overlap_operator(dg_frag, ispin_probe, vec, svec, .true.)
        c2(io) = real(sum(conjg(vec) * vec), kind=8)
        cs2(io) = real(sum(conjg(vec) * svec), kind=8)
      end do
    else
      do io = 1, nstate_probe
        vec(:) = (0.0d0, 0.0d0)
        if (n_frag_rows > 0) vec(1:n_frag_rows) = dg_frag%coef(1:n_frag_rows, io, ispin_probe)
        call apply_overlap_operator(dg_frag, ispin_probe, vec, svec, .true.)
        c2(io) = real(sum(conjg(vec) * vec), kind=8)
        cs2(io) = real(sum(conjg(vec) * svec), kind=8)
      end do
    end if

    write(*,'(1x,a,a,a,3(1x,es12.4),a,3(1x,es12.4))') "        coef stage probe: stage=", trim(stage_label), &
      " c2=", c2(1), c2(2), c2(3), " cs2=", cs2(1), cs2(2), cs2(3)
    flush(6)

    deallocate(vec, svec)
  end subroutine debug_coef_metric

  end subroutine tddft_dg_fragment_iteration
