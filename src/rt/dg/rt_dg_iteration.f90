  subroutine tddft_dg_fragment_iteration(dg_frag, system, info, rt, itt, dt, &
                                         lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                         rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use communication, only: comm_summation
    use salmon_global, only: yn_fix_func, theory
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional
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
    integer :: coef_size, coef_pw_size
    complex(8), allocatable :: coef_pw_work(:,:,:)

    ! Time evolution in fragment basis coefficient space
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

    ! Synchronize propagated coefficients across ranks to remove rank-dependent drift.
    coef_size = dg_frag%n_mat_max * dg_frag%nstate_tot * dg_frag%nspin
    call comm_summation(dg_frag%coef, dg_frag%coef_work, coef_size, dg_frag%icomm)
    dg_frag%coef = dg_frag%coef_work / real(max(1, dg_frag%isize), 8)
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw) .and. dg_frag%n_plane_waves > 0) then
      coef_pw_size = dg_frag%n_plane_waves * dg_frag%nstate_tot * dg_frag%nspin
      allocate(coef_pw_work(dg_frag%n_plane_waves, dg_frag%nstate_tot, dg_frag%nspin))
      call comm_summation(dg_frag%coef_pw, coef_pw_work, coef_pw_size, dg_frag%icomm)
      dg_frag%coef_pw = coef_pw_work / real(max(1, dg_frag%isize), 8)
      deallocate(coef_pw_work)
    end if

    ! Enforce orthonormality after synchronization to preserve unitary evolution.
    call stabilize_coeff_unitarity(dg_frag, itt)

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
        write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        iteration stage: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " stage=", "before-update-density-hmat"
        flush(6)
        call update_density_and_hamiltonian(dg_frag, system, info, rt, itt, rt%Ac_tot(:,itt), &
                                            lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                            rho, rho_s, Vh, Vxc, Vpsl, energy)
        write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        iteration stage: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " stage=", "after-update-density-hmat"
        flush(6)
      end if
    end if
    
    ! Calculate observables
    call calculate_observables(dg_frag, system, mg, stencil, ppg, rt, itt)

    if (trim(theory) == 'single_scale_maxwell_tddft') then
      call calculate_microscopic_current_dg(dg_frag, system, mg, stencil, rt%j_e)
    end if
    
  end subroutine tddft_dg_fragment_iteration
