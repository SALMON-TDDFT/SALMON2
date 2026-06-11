  subroutine tddft_dg_fragment_iteration(dg_frag, system, info, rt, itt, dt, &
                                         lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                         rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use salmon_global, only: yn_fix_func, theory
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional
    use timer, only: timer_begin, timer_end, LOG_CALC_TIME_PROPAGATION, LOG_CALC_RHO
    use rt_dg_fragment_ops, only: zero_nonowned_coefficients
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
    logical, parameter :: trace_first_iteration = .false.

    ! Time evolution in fragment basis coefficient space
    if (trace_first_iteration) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') "[DG-ITER-FIRST] entry rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
        " integrator=", dg_frag%time_integrator
      flush(6)
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
    if (trace_first_iteration) then
      write(*,'(1x,a)') "[DG-ITER-FIRST] after time evolution"
      flush(6)
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
        if (trace_first_iteration) then
          write(*,'(1x,a)') "[DG-ITER-FIRST] before post-step density/H update"
          flush(6)
        end if
        call timer_begin(LOG_CALC_RHO)
        call update_density_and_hamiltonian(dg_frag, system, info, rt, itt, rt%Ac_tot(:,itt), &
                                            lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                            rho, rho_s, Vh, Vxc, Vpsl, energy)
        call timer_end(LOG_CALC_RHO)
        if (trace_first_iteration) then
          write(*,'(1x,a)') "[DG-ITER-FIRST] after post-step density/H update"
          flush(6)
        end if
      end if
    end if
    
    ! Calculate observables
    if (trace_first_iteration) then
      write(*,'(1x,a)') "[DG-ITER-FIRST] before observables"
      flush(6)
    end if
    call calculate_observables(dg_frag, system, mg, stencil, ppg, rt, itt, Vh, Vxc, Vpsl)
    if (trace_first_iteration) then
      write(*,'(1x,a)') "[DG-ITER-FIRST] after observables"
      flush(6)
    end if

    if (trim(theory) == 'single_scale_maxwell_tddft') then
      if (trace_first_iteration) then
        write(*,'(1x,a)') "[DG-ITER-FIRST] before microscopic current"
        flush(6)
      end if
      call calculate_microscopic_current_dg(dg_frag, system, mg, stencil, rt%j_e)
      if (trace_first_iteration) then
        write(*,'(1x,a)') "[DG-ITER-FIRST] after microscopic current"
        flush(6)
      end if
    end if
  end subroutine tddft_dg_fragment_iteration
