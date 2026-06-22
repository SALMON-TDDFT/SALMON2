  subroutine time_evolution_expdiag(dg_frag, system, info, rt, itt, dt, &
                                    lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                    rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use salmon_global, only: yn_fix_func, yn_dg_length_gauge, ae_shape1, e_impulse, epdir_re1, yn_restart, nt
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional
    use rt_dg_fragment_types, only: s_dg_fragment_rt
    use rt_dg_fragment_ops, only: ensure_nonlocal_pp_matrix_A
    use rt_dg_plane_wave, only: diagnose_wpw_reduced_density, diagnose_wpw_reduced_embed_local, &
      diagnose_wpw_reduced_embed_prodcoef, initialize_wpw_reduced_self_projection
    use eigen_subdiag_sub, only: eigen_zheev
    use communication, only: comm_summation
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

    integer :: it0, it1, ispin, ifrag, i_local, iblk
    integer :: nw, nbf, nstate_prop, state_first, state_last
    integer :: io, jo, iw, jw, istate, global_idx, local_idx, local_col
    integer :: ib, jb, ib_global, jb_global
    real(8) :: Ac_ham(3), E_mid(3), rdot_diag
    real(8) :: cphase, sphase
    complex(8), parameter :: zi = (0.0d0, 1.0d0)
    complex(8), allocatable :: c_w(:), tmp_w(:), next_w(:)
    complex(8), allocatable :: h_eff(:,:), evec(:,:)
    real(8), allocatable :: eval(:)
    logical :: use_bpw_wannier_h, use_formal_wannier_h
    logical :: use_formal_seed_phase
    logical :: use_global_flux_exp
    logical, save :: global_flux_env_checked = .false.
    logical, save :: global_flux_exp_enabled = .false.
    character(len=32) :: global_flux_env
    logical, save :: global_field_env_checked = .false.
    logical, save :: global_field_exp_enabled = .false.
    character(len=32) :: global_field_env
    logical, save :: formal_local_field_env_checked = .false.
    logical, save :: formal_local_field_enabled = .false.
    character(len=32) :: formal_local_field_env
    logical, save :: expdiag_warned = .false.
    logical, save :: xi_split_env_checked = .false.
    logical, save :: xi_split_enabled = .true.
    logical, save :: project_h_env_checked = .false.
    logical, save :: project_h_for_fixed_func = .false.
    logical, save :: mixed_z_env_checked = .false.
    logical, save :: mixed_z_enabled = .false.
    logical, save :: wpw_red_env_checked = .false.
    logical, save :: wpw_red_expdiag_enabled = .false.
    logical, save :: wpw_red_prodop_env_checked = .false.
    logical, save :: wpw_red_prodop_diag_enabled = .false.
    logical, save :: wpw_red_prodop_full_action_env_checked = .false.
    logical, save :: wpw_red_prodop_full_action_diag_enabled = .false.
    logical, save :: wpw_red_init_project_env_checked = .false.
    logical, save :: wpw_red_init_project_enabled = .false.
    logical, save :: timing_env_checked = .false.
    logical, save :: timing_enabled = .false.
    character(32) :: xi_split_env
    character(32) :: project_h_env
    character(32) :: mixed_z_env
    character(32) :: wpw_red_env
    character(32) :: wpw_red_prodop_env
    character(32) :: wpw_red_prodop_full_action_env
    character(32) :: wpw_red_init_project_env
    character(32) :: timing_env
    real(8) :: t_all0, t0, t1, t_route0
    real(8) :: time_update_h, time_nonlocal, time_mixed_field, time_mixed_phase
    real(8) :: time_global_flux, time_local_loop, time_local_diag
    real(8) :: time_gather, time_scatter, time_field_diag, time_field_matmul
    real(8) :: time_flux_project, time_flux_comm, time_flux_scatter
    logical :: wpw_red_bad_coef
    logical :: trace_wpw_order

    if (yn_dg_length_gauge /= 'y') stop "DG expdiag integrator currently requires length gauge"
    if (.not. timing_env_checked) then
      timing_env = ' '
      call get_environment_variable('SALMON_DG_EXPDIAG_TIMING', timing_env)
      select case (trim(adjustl(timing_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        timing_enabled = .true.
      case default
        timing_enabled = .false.
      end select
      timing_env_checked = .true.
    end if
    t_all0 = get_wtime()
    time_update_h = 0.0d0
    time_nonlocal = 0.0d0
    time_mixed_field = 0.0d0
    time_mixed_phase = 0.0d0
    time_global_flux = 0.0d0
    time_local_loop = 0.0d0
    time_local_diag = 0.0d0
    time_gather = 0.0d0
    time_scatter = 0.0d0
    time_field_diag = 0.0d0
    time_field_matmul = 0.0d0
    time_flux_project = 0.0d0
    time_flux_comm = 0.0d0
    time_flux_scatter = 0.0d0
    trace_wpw_order = wpw_reduced_step_trace_enabled()
    if (dg_frag%coef_state_block_mode) &
      stop "DG expdiag integrator does not yet support state-block coefficient ownership"
    use_bpw_wannier_h = dg_frag%buffer_wannier_flux_seed_applied .and. &
      dg_frag%has_buffer_periodic_wannier_basis .and. allocated(dg_frag%buffer_wannier_coef) .and. &
      allocated(dg_frag%buffer_wannier_h_flux) .and. allocated(dg_frag%buffer_wannier_v)
    use_formal_wannier_h = dg_frag%has_formal_dg_wannier_basis .and. &
      allocated(dg_frag%dg_wannier_basis_coef) .and. allocated(dg_frag%dg_wannier_h0_local) .and. &
      allocated(dg_frag%dg_wannier_xi_local) .and. allocated(dg_frag%dg_wannier_ref_center) .and. &
      allocated(dg_frag%dg_wannier_nkeep)
    if (.not. use_bpw_wannier_h .and. .not. use_formal_wannier_h) &
      stop "DG expdiag integrator requires formal DG-Wannier or buffer-periodic Wannier data"
    if (.not. global_flux_env_checked) then
      global_flux_env = ' '
      call get_environment_variable('SALMON_DG_EXPDIAG_GLOBAL_FLUX', global_flux_env)
      global_flux_exp_enabled = trim(global_flux_env) == '1' .or. &
        trim(global_flux_env) == 'y' .or. trim(global_flux_env) == 'Y'
      global_flux_env_checked = .true.
    end if
    if (.not. global_field_env_checked) then
      global_field_env = ' '
      call get_environment_variable('SALMON_DG_EXPDIAG_GLOBAL_FIELD', global_field_env)
      select case (trim(adjustl(global_field_env)))
      case ('0','n','N','no','NO','false','FALSE','off','OFF')
        global_field_exp_enabled = .false.
      case default
        global_field_exp_enabled = .true.
      end select
      global_field_env_checked = .true.
    end if
    if (.not. xi_split_env_checked) then
      xi_split_env = ' '
      call get_environment_variable('SALMON_DG_EXPDIAG_XI_SPLIT', xi_split_env)
      select case (trim(adjustl(xi_split_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        xi_split_enabled = .true.
      case ('0','n','N','no','NO','false','FALSE','off','OFF')
        xi_split_enabled = .false.
      case default
        xi_split_enabled = .true.
      end select
      xi_split_env_checked = .true.
    end if
    if (.not. project_h_env_checked) then
      project_h_env = ' '
      call get_environment_variable('SALMON_DG_EXPDIAG_PROJECT_H', project_h_env)
      select case (trim(adjustl(project_h_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        project_h_for_fixed_func = .true.
      case default
        project_h_for_fixed_func = .false.
      end select
      project_h_env_checked = .true.
    end if
    if (.not. mixed_z_env_checked) then
      mixed_z_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z', mixed_z_env)
      select case (trim(adjustl(mixed_z_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_enabled = .true.
      case default
        mixed_z_enabled = .false.
      end select
      mixed_z_env_checked = .true.
    end if
    if (.not. wpw_red_env_checked) then
      wpw_red_env = ' '
      call get_environment_variable('SALMON_DG_WPW_REDUCED_EXPDIAG', wpw_red_env)
      select case (trim(adjustl(wpw_red_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        wpw_red_expdiag_enabled = .true.
      case default
        wpw_red_expdiag_enabled = .false.
      end select
      wpw_red_env_checked = .true.
    end if
    if (.not. wpw_red_prodop_env_checked) then
      wpw_red_prodop_env = ' '
      call get_environment_variable('SALMON_DG_WPW_REDUCED_PRODOP_DIAG', wpw_red_prodop_env)
      select case (trim(adjustl(wpw_red_prodop_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        wpw_red_prodop_diag_enabled = .true.
      case default
        wpw_red_prodop_diag_enabled = .false.
      end select
      if (.not. wpw_red_prodop_diag_enabled) then
        wpw_red_prodop_env = ' '
        call get_environment_variable('SALMON_DG_WPW_REDUCED_PRODOP_ACTION_DIAG', wpw_red_prodop_env)
        select case (trim(adjustl(wpw_red_prodop_env)))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          wpw_red_prodop_diag_enabled = .true.
        case default
          wpw_red_prodop_diag_enabled = .false.
        end select
      end if
      wpw_red_prodop_env_checked = .true.
    end if
    if (.not. wpw_red_prodop_full_action_env_checked) then
      wpw_red_prodop_full_action_env = ' '
      call get_environment_variable('SALMON_DG_WPW_REDUCED_PRODOP_FULL_ACTION_DIAG', &
        wpw_red_prodop_full_action_env)
      select case (trim(adjustl(wpw_red_prodop_full_action_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        wpw_red_prodop_full_action_diag_enabled = .true.
      case default
        wpw_red_prodop_full_action_diag_enabled = .false.
      end select
      wpw_red_prodop_full_action_env_checked = .true.
    end if
    if (.not. wpw_red_init_project_env_checked) then
      wpw_red_init_project_env = ' '
      call get_environment_variable('SALMON_DG_WPW_REDUCED_INIT_PROJECT', wpw_red_init_project_env)
      select case (trim(adjustl(wpw_red_init_project_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        wpw_red_init_project_enabled = .true.
      case default
        wpw_red_init_project_enabled = .false.
      end select
      wpw_red_init_project_env_checked = .true.
    end if
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw) .and. .not. mixed_z_enabled) &
      stop "DG expdiag integrator does not yet support BPW+PW mixed propagation"

    it0 = max(lbound(rt%Ac_tot, 2), itt - 1)
    it1 = min(ubound(rt%Ac_tot, 2), itt)
    Ac_ham(:) = 0.0d0
    E_mid(:) = 0.5d0 * (rt%E_tot(:, it0) + rt%E_tot(:, it1))
    if (trim(ae_shape1) == 'impulse' .and. yn_restart == 'n' .and. itt == 1) then
      E_mid(:) = e_impulse * epdir_re1(:) / max(abs(dt), 1.0d-300)
      if (dg_frag%id == 0) then
        write(*,'(1x,a,3(1x,1pe13.5),a,1pe13.5)') &
          '[DG-LG-IMPULSE] using first-step rectangular E field=', E_mid(:), ' dt=', dt
        flush(6)
      end if
    end if
    if (trace_wpw_order .and. dg_frag%id == 0) then
      write(*,'(1x,a,a,i0,a,1pe12.4,3(a,1pe12.4),3(a,1pe12.4),4(a,l1))') &
        '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, ' dt=', dt, &
        ' E_x=', E_mid(1), ' E_y=', E_mid(2), ' E_z=', E_mid(3), &
        ' Ac_x=', Ac_ham(1), ' Ac_y=', Ac_ham(2), ' Ac_z=', Ac_ham(3), &
        ' production_impulse_rect=', trim(ae_shape1) == 'impulse' .and. yn_restart == 'n' .and. itt == 1, &
        ' density_hxc_update_next=', yn_fix_func == 'n', &
        ' aux_reduced_next=', wpw_red_expdiag_enabled, &
        ' production_mixed_z=', mixed_z_enabled
    end if

    if (yn_fix_func == 'n') then
      if (trace_wpw_order .and. dg_frag%id == 0) then
        write(*,'(1x,a,a,i0,a)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
          ' stage=production-density-hamiltonian-update-begin'
      end if
      t0 = get_wtime()
      call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_ham, &
                                            lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                            rho, rho_s, Vh, Vxc, Vpsl, energy, .false.)
      time_update_h = time_update_h + (get_wtime() - t0)
      if (trace_wpw_order .and. dg_frag%id == 0) then
        write(*,'(1x,a,a,i0,a,1pe12.4)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
          ' stage=production-density-hamiltonian-update-end elapsed=', time_update_h
      end if
    end if
    if (ppg%Nlma > 0 .and. allocated(ppg%uV)) then
      if (trace_wpw_order .and. dg_frag%id == 0) then
        write(*,'(1x,a,a,i0,a)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
          ' stage=production-nonlocal-update-begin'
      end if
      t0 = get_wtime()
      call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_ham, .false.)
      time_nonlocal = time_nonlocal + (get_wtime() - t0)
      if (trace_wpw_order .and. dg_frag%id == 0) then
        write(*,'(1x,a,a,i0,a,1pe12.4)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
          ' stage=production-nonlocal-update-end elapsed=', time_nonlocal
      end if
    end if

    nstate_prop = dg_frag%nstate_tot
    if (allocated(dg_frag%nocc_spin)) then
      nstate_prop = min(dg_frag%nstate_tot, max(1, maxval(dg_frag%nocc_spin(1:dg_frag%nspin))))
    end if
    state_first = 1
    state_last = min(nstate_prop, size(dg_frag%coef, 2))

    call diagnose_wpw_reduced_embed_local(dg_frag, itt)
    call diagnose_wpw_reduced_embed_prodcoef(dg_frag, itt)

    if (wpw_red_expdiag_enabled) then
      if (trace_wpw_order .and. dg_frag%id == 0) then
        write(*,'(1x,a,a,i0,a)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
          ' stage=aux-reduced-static-expdiag-begin'
      end if
      call dryrun_wpw_reduced_expdiag(state_first, state_last, dt, wpw_red_bad_coef, itt, E_mid, Ac_ham, &
        'before-production-propagation')
      if (trace_wpw_order .and. dg_frag%id == 0) then
        write(*,'(1x,a,a,i0,a,l1)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
          ' stage=aux-reduced-static-expdiag-end bad_coef=', wpw_red_bad_coef
      end if
      call diagnose_wpw_reduced_density(dg_frag, system, rho, itt, nt, wpw_red_bad_coef, dt, &
        reproject_stage='before-production')
    end if

    use_global_flux_exp = use_formal_wannier_h .and. global_flux_exp_enabled .and. &
      dg_frag%has_global_wannier_flux_eigen .and. allocated(dg_frag%global_wannier_flux_evec) .and. &
      allocated(dg_frag%global_wannier_flux_eval) .and. allocated(dg_frag%global_wannier_coef)
    if (use_global_flux_exp) then
      if (mixed_z_enabled) call ensure_mixed_wannier_bpw_position(dg_frag)
      if (mixed_z_enabled .and. dg_frag%has_mixed_wannier_bpw_position) then
        if (sum(abs(E_mid(1:3))) > 1.0d-30) then
          if (trace_wpw_order .and. dg_frag%id == 0) then
            write(*,'(1x,a,a,i0,a)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
              ' stage=production-mixed-z-field-kick-begin'
          end if
          t_route0 = get_wtime()
          call apply_global_mixed_position_field_exp(E_mid, state_first, state_last)
          time_mixed_field = time_mixed_field + (get_wtime() - t_route0)
          if (trace_wpw_order .and. dg_frag%id == 0) then
            write(*,'(1x,a,a,i0,a,1pe12.4)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
              ' stage=production-mixed-z-field-kick-end elapsed=', time_mixed_field
          end if
        end if
        if (trace_wpw_order .and. dg_frag%id == 0) then
          write(*,'(1x,a,a,i0,a)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
            ' stage=production-mixed-z-flux-phase-begin'
        end if
        t_route0 = get_wtime()
        call apply_global_mixed_phase_exp(state_first, state_last)
        time_mixed_phase = time_mixed_phase + (get_wtime() - t_route0)
        if (trace_wpw_order .and. dg_frag%id == 0) then
          write(*,'(1x,a,a,i0,a,1pe12.4)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
            ' stage=production-mixed-z-flux-phase-end elapsed=', time_mixed_phase
        end if
        if (wpw_red_prodop_diag_enabled) then
          call diagnose_wpw_reduced_prodop_action('mixed-z', itt, state_first, state_last, &
            sum(abs(E_mid(1:3))) > 1.0d-30, .true., .false.)
          call diagnose_wpw_reduced_density(dg_frag, system, rho, itt, nt, wpw_red_bad_coef, dt, &
            reproject_stage='after-production', prodop_route='mixed-z', &
            prodop_field_included=sum(abs(E_mid(1:3))) > 1.0d-30, &
            prodop_mixed_z_included=.true., prodop_global_flux_included=.false., &
            prodop_kick_applied=trim(ae_shape1) == 'impulse' .and. yn_restart == 'n' .and. itt == 1, &
            prodop_predictor_corrector_included=.false.)
        end if
        if (wpw_red_prodop_full_action_diag_enabled) then
          call diagnose_wpw_reduced_prodop_full_action('mixed-z', itt, state_first, state_last, &
            sum(abs(E_mid(1:3))) > 1.0d-30, &
            trim(ae_shape1) == 'impulse' .and. yn_restart == 'n' .and. itt == 1, &
            .true., .false., .false.)
        end if
        if (.not. expdiag_warned .and. dg_frag%id == 0) then
          write(*,'(1x,a)') "[DG-EXPDIAG] dense Wannier+BPW-perp mixed-Z integrator enabled."
          write(*,'(1x,a)') "[DG-EXPDIAG] kick and polarization use the same mixed AA_R/BPW-perp Z when enabled."
          flush(6)
          expdiag_warned = .true.
        end if
        call print_expdiag_timing('mixed-z')
        return
      end if
      if (sum(abs(E_mid(1:3))) > 1.0d-30) then
        if (trace_wpw_order .and. dg_frag%id == 0) then
          write(*,'(1x,a,a,i0,a,l1)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
            ' stage=production-global-flux-field-kick-begin global_field=', global_field_exp_enabled
        end if
        if (.not. global_field_exp_enabled) then
          call apply_formal_local_field_exp(E_mid, state_first, state_last)
        else if (dg_frag%has_global_wannier_position .and. allocated(dg_frag%global_wannier_position)) then
          call apply_global_wannier_position_field_exp(E_mid, state_first, state_last)
        else
          stop "DG expdiag global Flux length-gauge path requires Wannier AA_R position matrix"
        end if
        if (trace_wpw_order .and. dg_frag%id == 0) then
          write(*,'(1x,a,a,i0,a)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
            ' stage=production-global-flux-field-kick-end'
        end if
      end if
      if (trace_wpw_order .and. dg_frag%id == 0) then
        write(*,'(1x,a,a,i0,a)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
          ' stage=production-global-flux-phase-begin'
      end if
      t_route0 = get_wtime()
      call apply_global_flux_eigen_exp(state_first, state_last)
      time_global_flux = time_global_flux + (get_wtime() - t_route0)
      if (trace_wpw_order .and. dg_frag%id == 0) then
        write(*,'(1x,a,a,i0,a,1pe12.4)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
          ' stage=production-global-flux-phase-end elapsed=', time_global_flux
      end if
      if (wpw_red_prodop_diag_enabled) then
        call diagnose_wpw_reduced_prodop_action('global-flux', itt, state_first, state_last, &
          sum(abs(E_mid(1:3))) > 1.0d-30, .false., .true.)
        call diagnose_wpw_reduced_density(dg_frag, system, rho, itt, nt, wpw_red_bad_coef, dt, &
          reproject_stage='after-production', prodop_route='global-flux', &
          prodop_field_included=sum(abs(E_mid(1:3))) > 1.0d-30, &
          prodop_mixed_z_included=.false., prodop_global_flux_included=.true., &
          prodop_kick_applied=trim(ae_shape1) == 'impulse' .and. yn_restart == 'n' .and. itt == 1, &
          prodop_predictor_corrector_included=.false.)
      end if
      if (wpw_red_prodop_full_action_diag_enabled) then
        call diagnose_wpw_reduced_prodop_full_action('global-flux', itt, state_first, state_last, &
          sum(abs(E_mid(1:3))) > 1.0d-30, &
          trim(ae_shape1) == 'impulse' .and. yn_restart == 'n' .and. itt == 1, &
          .false., .true., .false.)
      end if
      if (.not. expdiag_warned .and. dg_frag%id == 0) then
        write(*,'(1x,a)') "[DG-EXPDIAG] global Flux-eigen projection integrator enabled."
        write(*,'(1x,a)') "[DG-EXPDIAG] this diagnostic path is dense/allreduce and is not weak-scaling."
        if (global_field_exp_enabled) then
          write(*,'(1x,a)') "[DG-EXPDIAG] nonzero field steps use a global Wannier AA_R field kick before Flux phase."
        else
          write(*,'(1x,a)') "[DG-EXPDIAG] nonzero field steps use a local DG-Wannier field kick before Flux phase."
        end if
        flush(6)
        expdiag_warned = .true.
      end if
      call print_expdiag_timing('global-flux')
      return
    end if

    use_formal_seed_phase = use_formal_wannier_h .and. allocated(dg_frag%esp) .and. &
      sum(abs(E_mid(1:3))) <= 1.0d-30
    if (use_formal_seed_phase) then
      do ispin = 1, dg_frag%nspin
        if (ispin > size(dg_frag%esp, 2)) cycle
        do istate = state_first, state_last
          if (istate > size(dg_frag%esp, 1)) cycle
          cphase = cos(dg_frag%esp(istate,ispin) * dt)
          sphase = sin(dg_frag%esp(istate,ispin) * dt)
          dg_frag%coef(1:size(dg_frag%coef,1),istate,ispin) = &
            cmplx(cphase, -sphase, kind=8) * dg_frag%coef(1:size(dg_frag%coef,1),istate,ispin)
        end do
      end do
      if (.not. expdiag_warned .and. dg_frag%id == 0) then
        write(*,'(1x,a)') "[DG-EXPDIAG] formal Flux-eigen phase integrator enabled for field-free steps."
        write(*,'(1x,a)') "[DG-EXPDIAG] nonzero length-gauge field steps require a Flux-eigen projection path."
        flush(6)
        expdiag_warned = .true.
      end if
      if (wpw_red_prodop_diag_enabled) then
        call diagnose_wpw_reduced_prodop_action('formal-seed-phase', itt, state_first, state_last, &
          .false., .false., .false.)
        call diagnose_wpw_reduced_density(dg_frag, system, rho, itt, nt, wpw_red_bad_coef, dt, &
          reproject_stage='after-production', prodop_route='formal-seed-phase', &
          prodop_field_included=.false., prodop_mixed_z_included=.false., &
          prodop_global_flux_included=.false., prodop_kick_applied=.false., &
          prodop_predictor_corrector_included=.false.)
      end if
      if (wpw_red_prodop_full_action_diag_enabled) then
        call diagnose_wpw_reduced_prodop_full_action('formal-seed-phase', itt, state_first, state_last, &
          .false., .false., .false., .false., .false.)
      end if
      call print_expdiag_timing('formal-seed-phase')
      return
    end if
    if (use_formal_wannier_h) then
      if (.not. formal_local_field_env_checked) then
        formal_local_field_env = ' '
        call get_environment_variable('SALMON_DG_EXPDIAG_FORMAL_LOCAL_FIELD', formal_local_field_env)
        formal_local_field_enabled = trim(formal_local_field_env) == '1' .or. &
          trim(formal_local_field_env) == 'y' .or. trim(formal_local_field_env) == 'Y'
        formal_local_field_env_checked = .true.
      end if
      if (.not. formal_local_field_enabled) then
        if (dg_frag%id == 0) then
          write(*,'(1x,a)') "[FATAL] formal DG-Wannier expdiag received a nonzero length-gauge field."
          write(*,'(1x,a)') "[FATAL] The local-only field path is disabled because it gives nonphysical P(t)."
          write(*,'(1x,a)') "[FATAL] Implement Flux-eigen projection or neighbor-extended exponential propagation."
        end if
        stop "DG expdiag formal length-gauge field path is not implemented"
      end if
    end if
    t_route0 = get_wtime()
    if (xi_split_enabled .and. use_bpw_wannier_h) call apply_xi_flux_split_half(E_mid, 0.5d0 * dt, state_first, state_last)
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = ifrag - dg_frag%ifrag_start + 1
      if (use_formal_wannier_h) then
        if (i_local < 1 .or. i_local > size(dg_frag%dg_wannier_nkeep)) cycle
        nw = dg_frag%dg_wannier_nkeep(i_local)
        nbf = min(dg_frag%n_basis(ifrag, 1), dg_frag%nstate_frag, size(dg_frag%dg_wannier_basis_coef, 1))
      else
        if (i_local < 1 .or. i_local > size(dg_frag%buffer_wannier_nkeep)) cycle
        nw = dg_frag%buffer_wannier_nkeep(i_local)
        nbf = min(dg_frag%n_basis(ifrag, 1), dg_frag%nstate_frag, size(dg_frag%buffer_wannier_coef, 1))
      end if
      if (nw <= 0 .or. nbf <= 0) cycle

      allocate(h_eff(nw,nw), eval(nw), evec(nw,nw), c_w(nw), tmp_w(nw), next_w(nw))
      do ispin = 1, dg_frag%nspin
        if (use_formal_wannier_h) then
          h_eff(1:nw,1:nw) = dg_frag%dg_wannier_h0_local(1:nw,1:nw,ispin,i_local)
        else
          h_eff(1:nw,1:nw) = cmplx(dg_frag%buffer_wannier_h_flux(1:nw,1:nw,i_local), 0.0d0, kind=8)
        end if
        iblk = 0
        if (allocated(dg_frag%H_block_map)) iblk = find_matrix_block(dg_frag%H_block_map, ifrag, ifrag)
        if ((yn_fix_func == 'n' .or. project_h_for_fixed_func) .and. iblk > 0 .and. allocated(dg_frag%H_mat_blocks)) then
          if (iblk <= size(dg_frag%H_mat_blocks)) then
            h_eff(1:nw,1:nw) = (0.0d0, 0.0d0)
            do jw = 1, nw
              do iw = 1, nw
                do jb = 1, min(nbf, size(dg_frag%H_mat_blocks(iblk)%val, 2))
                  do ib = 1, min(nbf, size(dg_frag%H_mat_blocks(iblk)%val, 1))
                    if (use_formal_wannier_h) then
                      h_eff(iw,jw) = h_eff(iw,jw) + &
                        conjg(dg_frag%dg_wannier_basis_coef(ib,iw,ispin,i_local)) * &
                        cmplx(dg_frag%H_mat_blocks(iblk)%val(ib,jb,ispin), 0.0d0, kind=8) * &
                        dg_frag%dg_wannier_basis_coef(jb,jw,ispin,i_local)
                    else
                      h_eff(iw,jw) = h_eff(iw,jw) + &
                        cmplx(dg_frag%buffer_wannier_coef(ib,iw,ispin,i_local), 0.0d0, kind=8) * &
                        cmplx(dg_frag%H_mat_blocks(iblk)%val(ib,jb,ispin), 0.0d0, kind=8) * &
                        cmplx(dg_frag%buffer_wannier_coef(jb,jw,ispin,i_local), 0.0d0, kind=8)
                    end if
                  end do
                end do
              end do
            end do
          end if
        end if
        do iw = 1, nw
          do jw = 1, nw
            if (use_formal_wannier_h) then
              h_eff(jw,iw) = h_eff(jw,iw) &
                - E_mid(1) * dg_frag%dg_wannier_xi_local(1,jw,iw,ispin,i_local) &
                - E_mid(2) * dg_frag%dg_wannier_xi_local(2,jw,iw,ispin,i_local) &
                - E_mid(3) * dg_frag%dg_wannier_xi_local(3,jw,iw,ispin,i_local)
            else
              h_eff(jw,iw) = h_eff(jw,iw) &
                - E_mid(1) * dg_frag%buffer_wannier_v(1,jw,iw,i_local) &
                - E_mid(2) * dg_frag%buffer_wannier_v(2,jw,iw,i_local) &
                - E_mid(3) * dg_frag%buffer_wannier_v(3,jw,iw,i_local)
            end if
          end do
          rdot_diag = 0.0d0
          if (use_formal_wannier_h) then
            if (i_local >= 1 .and. i_local <= size(dg_frag%dg_wannier_ref_center, 2)) then
              rdot_diag = E_mid(1) * dg_frag%dg_wannier_ref_center(1,i_local) &
                        + E_mid(2) * dg_frag%dg_wannier_ref_center(2,i_local) &
                        + E_mid(3) * dg_frag%dg_wannier_ref_center(3,i_local)
            end if
          else
            if (allocated(dg_frag%buffer_wannier_frag_center)) then
              if (i_local >= 1 .and. i_local <= size(dg_frag%buffer_wannier_frag_center, 2)) then
                rdot_diag = E_mid(1) * dg_frag%buffer_wannier_frag_center(1,i_local) &
                          + E_mid(2) * dg_frag%buffer_wannier_frag_center(2,i_local) &
                          + E_mid(3) * dg_frag%buffer_wannier_frag_center(3,i_local)
              end if
            end if
          end if
          h_eff(iw,iw) = h_eff(iw,iw) - rdot_diag
        end do
        evec(1:nw,1:nw) = h_eff(1:nw,1:nw)
        t1 = get_wtime()
        call eigen_zheev(evec, eval, h_eff)
        time_local_diag = time_local_diag + (get_wtime() - t1)

        do istate = state_first, state_last
          c_w(:) = (0.0d0, 0.0d0)
          do iw = 1, nw
            do io = 1, nbf
              global_idx = dg_frag%index_basis(io, ifrag, ispin)
              local_idx = 0
              if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) &
                local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
              if (local_idx <= 0 .or. local_idx > size(dg_frag%coef, 1)) cycle
              if (use_formal_wannier_h) then
                c_w(iw) = c_w(iw) + conjg(dg_frag%dg_wannier_basis_coef(io, iw, ispin, i_local)) * &
                                      dg_frag%coef(local_idx, istate, ispin)
              else
                c_w(iw) = c_w(iw) + cmplx(dg_frag%buffer_wannier_coef(io, iw, ispin, i_local), 0.0d0, kind=8) * &
                                      dg_frag%coef(local_idx, istate, ispin)
              end if
            end do
          end do

          tmp_w(:) = matmul(conjg(transpose(h_eff(1:nw,1:nw))), c_w(:))
          do iw = 1, nw
            cphase = cos(eval(iw) * dt)
            sphase = sin(eval(iw) * dt)
            tmp_w(iw) = cmplx(cphase, -sphase, kind=8) * tmp_w(iw)
          end do
          next_w(:) = matmul(h_eff(1:nw,1:nw), tmp_w(:))

          do io = 1, nbf
            global_idx = dg_frag%index_basis(io, ifrag, ispin)
            local_idx = 0
            if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) &
              local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
            if (local_idx <= 0 .or. local_idx > size(dg_frag%coef, 1)) cycle
            dg_frag%coef(local_idx, istate, ispin) = (0.0d0, 0.0d0)
            do iw = 1, nw
              if (use_formal_wannier_h) then
                dg_frag%coef(local_idx, istate, ispin) = dg_frag%coef(local_idx, istate, ispin) + &
                  dg_frag%dg_wannier_basis_coef(io, iw, ispin, i_local) * next_w(iw)
              else
                dg_frag%coef(local_idx, istate, ispin) = dg_frag%coef(local_idx, istate, ispin) + &
                  cmplx(dg_frag%buffer_wannier_coef(io, iw, ispin, i_local), 0.0d0, kind=8) * next_w(iw)
              end if
            end do
          end do
        end do
      end do
      deallocate(h_eff, eval, evec, c_w, tmp_w, next_w)
    end do
    if (xi_split_enabled .and. use_bpw_wannier_h) call apply_xi_flux_split_half(E_mid, 0.5d0 * dt, state_first, state_last)
    time_local_loop = time_local_loop + (get_wtime() - t_route0)

    if (.not. expdiag_warned .and. dg_frag%id == 0) then
      if (use_formal_wannier_h) then
        write(*,'(1x,a)') "[DG-EXPDIAG] formal DG-Wannier exponential integrator enabled."
      else
        write(*,'(1x,a)') "[DG-EXPDIAG] experimental local BPW exponential integrator enabled."
      end if
      if (xi_split_enabled .and. use_bpw_wannier_h) then
        write(*,'(1x,a)') "[DG-EXPDIAG] dynamic neighbor xi_flux is applied by a Strang half-step split."
      else if (use_formal_wannier_h) then
        write(*,'(1x,a)') "[DG-EXPDIAG] nonzero field step uses local DG-Wannier length-gauge blocks."
      else
        write(*,'(1x,a)') &
          "[DG-EXPDIAG] dynamic neighbor xi_flux split disabled; set SALMON_DG_EXPDIAG_XI_SPLIT=1 to test."
      end if
      if (yn_fix_func == 'y' .and. .not. project_h_for_fixed_func) then
        write(*,'(1x,a)') "[DG-EXPDIAG] fixed-function propagation uses the seed flux Hamiltonian."
      else
        write(*,'(1x,a)') "[DG-EXPDIAG] propagation projects the current DG Hamiltonian into the active Wannier basis."
      end if
      flush(6)
      expdiag_warned = .true.
    end if
    call print_expdiag_timing('local-block')

  contains

    subroutine diagnose_wpw_reduced_prodop_action(route, step_use, state_s, state_e, &
        field_terms_included, mixed_z_terms_included, global_flux_terms_included)
      character(*), intent(in) :: route
      integer, intent(in) :: step_use, state_s, state_e
      logical, intent(in) :: field_terms_included, mixed_z_terms_included, global_flux_terms_included
      logical, save :: warned_full_action = .false.

      if (dg_frag%id /= 0) return
      write(*,'(1x,a,3(a,i0),a,a,5(a,l1))') &
        '[DG-WPW-RED-DIAG-PRODOP]', &
        ' step=', step_use, ' state_first=', state_s, ' state_last=', state_e, &
        ' route=', trim(route), &
        ' sampled_raw_state_action=', .true., &
        ' full_reduced_basis_action=', .false., &
        ' field_terms_included=', field_terms_included, &
        ' mixed_z_terms_included=', mixed_z_terms_included, &
        ' global_flux_terms_included=', global_flux_terms_included
      if (.not. warned_full_action) then
        write(*,'(1x,a)') &
          '[DG-WPW-RED-DIAG-PRODOP] full reduced-basis column action is not formed here because production mixed/global coefficient space differs from fragment-local WPW reduced space.'
        write(*,'(1x,a)') &
          '[DG-WPW-RED-DIAG-PRODOP] following density/coef diagnostics compare static auxiliary prediction against the reprojected post-production raw state.'
        warned_full_action = .true.
      end if
    end subroutine diagnose_wpw_reduced_prodop_action

    subroutine diagnose_wpw_reduced_prodop_full_action(route, step_use, state_s, state_e, &
        field_terms_included, kick_applied, mixed_z_terms_included, global_flux_terms_included, &
        predictor_corrector_included)
      character(*), intent(in) :: route
      integer, intent(in) :: step_use, state_s, state_e
      logical, intent(in) :: field_terms_included, kick_applied
      logical, intent(in) :: mixed_z_terms_included, global_flux_terms_included
      logical, intent(in) :: predictor_corrector_included
      real(8) :: nred_local(1), nred_global(1)
      integer :: i_local
      logical :: embedding_available, bad_full_action
      logical, save :: warned_embedding = .false.

      nred_local(1) = 0.0d0
      nred_global(1) = 0.0d0
      if (dg_frag%wpw_reduced_ready .and. allocated(dg_frag%wpw_reduced_dim)) then
        do i_local = 1, size(dg_frag%wpw_reduced_dim)
          nred_local(1) = nred_local(1) + dble(max(0, dg_frag%wpw_reduced_dim(i_local)))
        end do
      end if
      call comm_summation(nred_local, nred_global, 1, dg_frag%icomm)

      embedding_available = .false.
      bad_full_action = .not. embedding_available
      if (dg_frag%id /= 0) return

      write(*,'(1x,a,2(a,i0),a,a,8(a,l1),a,1pe12.4)') &
        '[DG-WPW-RED-DIAG-PRODOP-FULL-ACTION]', &
        ' step=', step_use, ' nred=', nint(nred_global(1)), &
        ' route=', trim(route), &
        ' field_on=', field_terms_included, &
        ' kick_applied=', kick_applied, &
        ' mixed_z_included=', mixed_z_terms_included, &
        ' global_flux_included=', global_flux_terms_included, &
        ' predictor_corrector_included=', predictor_corrector_included, &
        ' sampled_state_action=', .false., &
        ' full_reduced_basis_action=', .false., &
        ' fixed_context=', .true., &
        ' rel_matrix_diff_Snorm=', 0.0d0
      write(*,'(1x,a,2(a,1pe12.4),a,1pe12.4,a,1pe12.4,a,1pe12.4,2(a,l1))') &
        '[DG-WPW-RED-DIAG-PRODOP-FULL-ACTION]', &
        ' rel_state_action_diff_Snorm=', 0.0d0, &
        ' rel_coef_diff_Snorm=', 0.0d0, &
        ' coef_diff_Snorm=', 0.0d0, &
        ' overlap_S_normed=', 0.0d0, &
        ' Pz_diff=', 0.0d0, &
        ' embedding_available=', embedding_available, &
        ' bad=', bad_full_action
      if (.not. warned_embedding) then
        write(*,'(1x,a)') &
          '[DG-WPW-RED-DIAG-PRODOP-FULL-ACTION] full column action is diagnostic-only and currently blocked: WPW-reduced -> production mixed/global coefficient embedding is not implemented.'
        write(*,'(1x,a)') &
          '[DG-WPW-RED-DIAG-PRODOP-FULL-ACTION] production P space is global BPW-perp while WPW reduced P space is fragment-local/windowed; snapshot/apply/restore will be enabled only after that bridge exists.'
        warned_embedding = .true.
      end if
    end subroutine diagnose_wpw_reduced_prodop_full_action

    logical function wpw_reduced_expdiag_trace_enabled() result(enabled)
      character(len=32) :: env_trace
      integer :: env_len, env_stat
      logical, save :: initialized = .false.
      logical, save :: cached_enabled = .false.

      if (.not. initialized) then
        env_trace = ''
        call get_environment_variable('SALMON_DG_WPW_REDUCED_EXPDIAG_TRACE', env_trace, length=env_len, status=env_stat)
        if (env_stat == 0 .and. env_len > 0) then
          select case (adjustl(trim(env_trace(1:env_len))))
          case ('1','y','Y','yes','YES','true','TRUE','on','ON')
            cached_enabled = .true.
          case default
            cached_enabled = .false.
          end select
        end if
        initialized = .true.
      end if
      enabled = cached_enabled
    end function wpw_reduced_expdiag_trace_enabled

    logical function wpw_reduced_step_trace_enabled() result(enabled)
      character(len=32) :: env_trace
      integer :: env_len, env_stat
      logical, save :: checked = .false.
      logical, save :: cached = .false.

      if (.not. checked) then
        env_trace = ''
        call get_environment_variable('SALMON_DG_WPW_REDUCED_STEP_TRACE', env_trace, length=env_len, status=env_stat)
        if (env_stat == 0 .and. env_len > 0) then
          select case (adjustl(trim(env_trace(1:env_len))))
          case ('1','y','Y','yes','YES','true','TRUE','on','ON')
            cached = .true.
          case default
            cached = .false.
          end select
        else
          cached = .false.
        end if
        checked = .true.
      end if
      enabled = cached
    end function wpw_reduced_step_trace_enabled

    subroutine dryrun_wpw_reduced_expdiag(state_s, state_e, dt_use, bad_coef_any, step_use, E_use, Ac_use, stage_label)
      integer, intent(in) :: state_s, state_e
      real(8), intent(in) :: dt_use
      logical, intent(out), optional :: bad_coef_any
      integer, intent(in), optional :: step_use
      real(8), intent(in), optional :: E_use(3), Ac_use(3)
      character(*), intent(in), optional :: stage_label

      integer :: iloc_frag, ifrag_use, ispin_use, nred, nself, nneigh, nstate_store
      integer :: ist, iw_use, io_use, n_w_use, nbf_use, global_idx_use, local_idx_use
      integer :: lwork_red, info_red, ired, jred
      integer :: trace_step
      real(8) :: eig_min, eig_max, s_min_red, s_max_red, h_herm_red
      real(8) :: norm_before, norm_after, norm_diff, norm_diff_max, norm_before_sum, norm_after_sum
      real(8) :: max_coef_abs, max_neighbor_abs
      real(8) :: cphase_red, sphase_red
      real(8) :: trace_E(3), trace_Ac(3)
      logical :: coef_reallocated, any_propagated, bad_coef, did_init_project, trace_expdiag, trace_step_enabled
      character(len=64) :: trace_stage
      complex(8), allocatable :: H_work_red(:,:), S_work_red(:,:), evec_red(:,:)
      complex(8), allocatable :: cvec(:), tmp_vec(:), avec(:), next_vec(:), work_red(:)
      real(8), allocatable :: eval_red(:), seval_red(:), rwork_red(:)
      external :: ZHEGV, ZHEEV

      if (present(bad_coef_any)) bad_coef_any = .false.
      trace_expdiag = wpw_reduced_expdiag_trace_enabled()
      trace_step_enabled = wpw_reduced_step_trace_enabled()
      trace_step = -1
      if (present(step_use)) trace_step = step_use
      trace_E(:) = 0.0d0
      if (present(E_use)) trace_E(:) = E_use(:)
      trace_Ac(:) = 0.0d0
      if (present(Ac_use)) trace_Ac(:) = Ac_use(:)
      trace_stage = 'unknown'
      if (present(stage_label)) trace_stage = stage_label
      if (trace_step_enabled .and. dg_frag%id == 0) then
        write(*,'(1x,a,a,a,a,i0,a,1pe12.4,3(a,1pe12.4),3(a,1pe12.4),5(a,l1))') &
          '[DG-WPW-RED-DIAG-STEP-TRACE]', &
          ' stage=', trim(trace_stage), ' step=', trace_step, ' dt=', dt_use, &
          ' E_x=', trace_E(1), ' E_y=', trace_E(2), ' E_z=', trace_E(3), &
          ' Ac_x=', trace_Ac(1), ' Ac_y=', trace_Ac(2), ' Ac_z=', trace_Ac(3), &
          ' aux_field_applied=', .false., &
          ' aux_kick_applied=', .false., &
          ' aux_predictor_corrector=', .false., &
          ' Hred_static=', .true., &
          ' density_updated_before_aux=', yn_fix_func == 'n'
      end if
      if (.not. dg_frag%wpw_reduced_ready) then
        if (trace_expdiag .and. dg_frag%id == 0) then
          write(*,'(1x,a)') '[DG-WPW-RED-EXPDIAG] requested but reduced H/S blocks are not ready; dry-run skipped.'
        end if
        return
      end if
      if (.not. allocated(dg_frag%wpw_reduced_H) .or. .not. allocated(dg_frag%wpw_reduced_S)) return
      if (.not. allocated(dg_frag%wpw_reduced_dim) .or. .not. allocated(dg_frag%wpw_reduced_nself)) return

      nstate_store = max(1, size(dg_frag%coef, 2))
      if (.not. allocated(dg_frag%wpw_reduced_eval)) then
        allocate(dg_frag%wpw_reduced_eval(dg_frag%wpw_reduced_max_dim, dg_frag%nspin, size(dg_frag%wpw_reduced_dim)))
      else if (size(dg_frag%wpw_reduced_eval,1) /= dg_frag%wpw_reduced_max_dim .or. &
               size(dg_frag%wpw_reduced_eval,2) /= dg_frag%nspin .or. &
               size(dg_frag%wpw_reduced_eval,3) /= size(dg_frag%wpw_reduced_dim)) then
        deallocate(dg_frag%wpw_reduced_eval)
        allocate(dg_frag%wpw_reduced_eval(dg_frag%wpw_reduced_max_dim, dg_frag%nspin, size(dg_frag%wpw_reduced_dim)))
      end if
      if (.not. allocated(dg_frag%wpw_reduced_evec)) then
        allocate(dg_frag%wpw_reduced_evec(dg_frag%wpw_reduced_max_dim, dg_frag%wpw_reduced_max_dim, &
          dg_frag%nspin, size(dg_frag%wpw_reduced_dim)))
      else if (size(dg_frag%wpw_reduced_evec,1) /= dg_frag%wpw_reduced_max_dim .or. &
               size(dg_frag%wpw_reduced_evec,2) /= dg_frag%wpw_reduced_max_dim .or. &
               size(dg_frag%wpw_reduced_evec,3) /= dg_frag%nspin .or. &
               size(dg_frag%wpw_reduced_evec,4) /= size(dg_frag%wpw_reduced_dim)) then
        deallocate(dg_frag%wpw_reduced_evec)
        allocate(dg_frag%wpw_reduced_evec(dg_frag%wpw_reduced_max_dim, dg_frag%wpw_reduced_max_dim, &
          dg_frag%nspin, size(dg_frag%wpw_reduced_dim)))
      end if
      coef_reallocated = .false.
      if (.not. allocated(dg_frag%coef_wpw_self)) then
        allocate(dg_frag%coef_wpw_self(dg_frag%wpw_reduced_max_dim, nstate_store, dg_frag%nspin, size(dg_frag%wpw_reduced_dim)))
        coef_reallocated = .true.
      else if (size(dg_frag%coef_wpw_self,1) /= dg_frag%wpw_reduced_max_dim .or. &
               size(dg_frag%coef_wpw_self,2) /= nstate_store .or. &
               size(dg_frag%coef_wpw_self,3) /= dg_frag%nspin .or. &
               size(dg_frag%coef_wpw_self,4) /= size(dg_frag%wpw_reduced_dim)) then
        deallocate(dg_frag%coef_wpw_self)
        allocate(dg_frag%coef_wpw_self(dg_frag%wpw_reduced_max_dim, nstate_store, dg_frag%nspin, size(dg_frag%wpw_reduced_dim)))
        coef_reallocated = .true.
      end if
      if (.not. allocated(dg_frag%coef_wpw_neighbor_reduced)) then
        allocate(dg_frag%coef_wpw_neighbor_reduced(dg_frag%wpw_reduced_max_dim, nstate_store, dg_frag%nspin, size(dg_frag%wpw_reduced_dim)))
        coef_reallocated = .true.
      else if (size(dg_frag%coef_wpw_neighbor_reduced,1) /= dg_frag%wpw_reduced_max_dim .or. &
               size(dg_frag%coef_wpw_neighbor_reduced,2) /= nstate_store .or. &
               size(dg_frag%coef_wpw_neighbor_reduced,3) /= dg_frag%nspin .or. &
               size(dg_frag%coef_wpw_neighbor_reduced,4) /= size(dg_frag%wpw_reduced_dim)) then
        deallocate(dg_frag%coef_wpw_neighbor_reduced)
        allocate(dg_frag%coef_wpw_neighbor_reduced(dg_frag%wpw_reduced_max_dim, nstate_store, dg_frag%nspin, size(dg_frag%wpw_reduced_dim)))
        coef_reallocated = .true.
      end if
      dg_frag%wpw_reduced_eval(:, :, :) = 0.0d0
      dg_frag%wpw_reduced_evec(:, :, :, :) = (0.0d0, 0.0d0)
      if (coef_reallocated .or. .not. dg_frag%wpw_reduced_coef_initialized) then
        dg_frag%coef_wpw_self(:, :, :, :) = (0.0d0, 0.0d0)
        dg_frag%coef_wpw_neighbor_reduced(:, :, :, :) = (0.0d0, 0.0d0)
        dg_frag%wpw_reduced_coef_initialized = .false.
      end if
      if (.not. dg_frag%wpw_reduced_coef_initialized .and. wpw_red_init_project_enabled) then
        call initialize_wpw_reduced_self_projection(dg_frag, state_s, state_e, did_init_project)
        if (did_init_project) dg_frag%wpw_reduced_coef_initialized = .true.
      end if

      any_propagated = .false.
      do iloc_frag = 1, size(dg_frag%wpw_reduced_dim)
        ifrag_use = dg_frag%ifrag_start + iloc_frag - 1
        nred = dg_frag%wpw_reduced_dim(iloc_frag)
        nself = dg_frag%wpw_reduced_nself(iloc_frag)
        if (nred <= 0 .or. nself <= 0) cycle
        nneigh = max(0, nred - nself)
        allocate(H_work_red(nred,nred), S_work_red(nred,nred), evec_red(nred,nred))
        allocate(cvec(nred), tmp_vec(nred), avec(nred), next_vec(nred))
        allocate(eval_red(nred), seval_red(nred), rwork_red(max(1, 3*nred - 2)), work_red(1))

        do ispin_use = 1, dg_frag%nspin
          H_work_red(:, :) = dg_frag%wpw_reduced_H(1:nred, 1:nred, ispin_use, iloc_frag)
          S_work_red(:, :) = dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin_use, iloc_frag)
          h_herm_red = 0.0d0
          do jred = 1, nred
            do ired = 1, nred
              h_herm_red = max(h_herm_red, abs(H_work_red(ired,jred) - conjg(H_work_red(jred,ired))))
            end do
          end do

          evec_red(:, :) = S_work_red(:, :)
          lwork_red = -1
          call ZHEEV('N', 'U', nred, evec_red, nred, seval_red, work_red, lwork_red, rwork_red, info_red)
          lwork_red = max(1, int(real(work_red(1), kind=8)))
          deallocate(work_red)
          allocate(work_red(lwork_red))
          evec_red(:, :) = S_work_red(:, :)
          call ZHEEV('N', 'U', nred, evec_red, nred, seval_red, work_red, lwork_red, rwork_red, info_red)
          if (info_red == 0) then
            s_min_red = seval_red(1)
            s_max_red = seval_red(nred)
          else
            s_min_red = 0.0d0
            s_max_red = 0.0d0
          end if
          deallocate(work_red)
          allocate(work_red(1))

          evec_red(:, :) = H_work_red(:, :)
          S_work_red(:, :) = dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin_use, iloc_frag)
          lwork_red = -1
          call ZHEGV(1, 'V', 'U', nred, evec_red, nred, S_work_red, nred, eval_red, work_red, lwork_red, rwork_red, info_red)
          lwork_red = max(1, int(real(work_red(1), kind=8)))
          deallocate(work_red)
          allocate(work_red(lwork_red))
          evec_red(:, :) = H_work_red(:, :)
          S_work_red(:, :) = dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin_use, iloc_frag)
          call ZHEGV(1, 'V', 'U', nred, evec_red, nred, S_work_red, nred, eval_red, work_red, lwork_red, rwork_red, info_red)
          if (info_red == 0) then
            eig_min = eval_red(1)
            eig_max = eval_red(nred)
            dg_frag%wpw_reduced_eval(1:nred, ispin_use, iloc_frag) = eval_red(1:nred)
            dg_frag%wpw_reduced_evec(1:nred, 1:nred, ispin_use, iloc_frag) = evec_red(1:nred, 1:nred)
          else
            eig_min = 0.0d0
            eig_max = 0.0d0
          end if

          norm_diff_max = 0.0d0
          norm_before_sum = 0.0d0
          norm_after_sum = 0.0d0
          max_coef_abs = 0.0d0
          max_neighbor_abs = 0.0d0
          bad_coef = .false.
          do ist = state_s, min(state_e, nstate_store)
            cvec(:) = (0.0d0, 0.0d0)
            if (dg_frag%wpw_reduced_coef_initialized) then
              cvec(1:nself) = dg_frag%coef_wpw_self(1:nself, ist, ispin_use, iloc_frag)
              if (nneigh > 0) cvec(nself+1:nred) = &
                dg_frag%coef_wpw_neighbor_reduced(1:nneigh, ist, ispin_use, iloc_frag)
            else if (allocated(dg_frag%global_wannier_local_coef) .and. allocated(dg_frag%global_wannier_local_nkeep)) then
              if (iloc_frag <= size(dg_frag%global_wannier_local_nkeep)) then
                n_w_use = min(dg_frag%global_wannier_local_nkeep(iloc_frag), nself, size(dg_frag%global_wannier_local_coef, 2))
                nbf_use = min(dg_frag%n_basis(ifrag_use, ispin_use), size(dg_frag%global_wannier_local_coef, 1))
                do iw_use = 1, n_w_use
                  do io_use = 1, nbf_use
                    global_idx_use = dg_frag%index_basis(io_use, ifrag_use, ispin_use)
                    local_idx_use = 0
                    if (global_idx_use > 0 .and. global_idx_use <= dg_frag%n_mat_max) &
                      local_idx_use = dg_frag%coef_global_to_local(global_idx_use, ispin_use)
                    if (local_idx_use <= 0 .or. local_idx_use > size(dg_frag%coef, 1)) cycle
                    cvec(iw_use) = cvec(iw_use) + &
                      conjg(dg_frag%global_wannier_local_coef(io_use, iw_use, ispin_use, iloc_frag)) * &
                      dg_frag%coef(local_idx_use, ist, ispin_use)
                  end do
                end do
              end if
            end if

            tmp_vec(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin_use, iloc_frag), cvec(:))
            norm_before = real(sum(conjg(cvec(:)) * tmp_vec(:)), kind=8)
            if (info_red == 0) then
              avec(:) = matmul(conjg(transpose(evec_red(:, :))), tmp_vec(:))
              do ired = 1, nred
                cphase_red = cos(eval_red(ired) * dt_use)
                sphase_red = sin(eval_red(ired) * dt_use)
                avec(ired) = cmplx(cphase_red, -sphase_red, kind=8) * avec(ired)
              end do
              next_vec(:) = matmul(evec_red(:, :), avec(:))
              tmp_vec(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin_use, iloc_frag), next_vec(:))
              norm_after = real(sum(conjg(next_vec(:)) * tmp_vec(:)), kind=8)
            else
              norm_after = norm_before
            end if
            norm_diff = abs(norm_after - norm_before)
            norm_diff_max = max(norm_diff_max, norm_diff)
            norm_before_sum = norm_before_sum + norm_before
            norm_after_sum = norm_after_sum + norm_after
            if (info_red == 0) then
              dg_frag%coef_wpw_self(1:nself, ist, ispin_use, iloc_frag) = next_vec(1:nself)
              if (nneigh > 0) dg_frag%coef_wpw_neighbor_reduced(1:nneigh, ist, ispin_use, iloc_frag) = &
                next_vec(nself+1:nred)
              max_coef_abs = max(max_coef_abs, maxval(abs(next_vec(1:nred))))
              if (nneigh > 0) max_neighbor_abs = max(max_neighbor_abs, maxval(abs(next_vec(nself+1:nred))))
              do ired = 1, nred
                if (real(next_vec(ired), kind=8) /= real(next_vec(ired), kind=8)) bad_coef = .true.
                if (aimag(next_vec(ired)) /= aimag(next_vec(ired))) bad_coef = .true.
                if (abs(next_vec(ired)) > huge(1.0d0) * 1.0d-100) bad_coef = .true.
              end do
              any_propagated = .true.
            end if
          end do

          if (trace_expdiag .and. dg_frag%id == 0) then
            write(*,'(1x,a,a,i0,5(a,i0),10(a,1pe12.4),a,l1)') &
              '[DG-WPW-RED-EXPDIAG] propagate:', &
              ' ifrag=', ifrag_use, ' dim=', nred, ' nself=', nself, &
              ' nkeep=', dg_frag%wpw_reduced_nkeep(iloc_frag), &
              ' ndrop=', dg_frag%wpw_reduced_ndrop(iloc_frag), ' info=', info_red, &
              ' eig_min=', eig_min, ' eig_max=', eig_max, &
              ' S_eval_min=', s_min_red, ' S_eval_max=', s_max_red, &
              ' H_herm=', h_herm_red, ' norm_before=', norm_before_sum, &
              ' norm_after=', norm_after_sum, ' norm_diff_max=', norm_diff_max, &
              ' max_coef=', max_coef_abs, ' max_neighbor_coef=', max_neighbor_abs, &
              ' bad_coef=', bad_coef
          end if
          if (present(bad_coef_any)) bad_coef_any = bad_coef_any .or. bad_coef
        end do

        deallocate(H_work_red, S_work_red, evec_red, cvec, tmp_vec, avec, next_vec)
        deallocate(eval_red, seval_red, rwork_red, work_red)
      end do
      if (any_propagated) dg_frag%wpw_reduced_coef_initialized = .true.
    end subroutine dryrun_wpw_reduced_expdiag

    subroutine print_expdiag_timing(route)
      character(*), intent(in) :: route
      real(8) :: total_time

      if (.not. timing_enabled) return
      if (dg_frag%id /= 0) return
      total_time = get_wtime() - t_all0
      write(*,'(1x,a,a,a,i0,14(a,1pe11.4))') '[DG-EXPDIAG-TIMING] route=', trim(route), &
        ' itt=', itt, &
        ' total=', total_time, &
        ' update_h=', time_update_h, &
        ' nonlocal=', time_nonlocal, &
        ' mixed_field=', time_mixed_field, &
        ' mixed_phase=', time_mixed_phase, &
        ' global_flux=', time_global_flux, &
        ' local_loop=', time_local_loop, &
        ' local_diag=', time_local_diag, &
        ' gather=', time_gather, &
        ' scatter=', time_scatter, &
        ' field_diag=', time_field_diag, &
        ' field_matmul=', time_field_matmul, &
        ' flux_project=', time_flux_project, &
        ' flux_comm=', time_flux_comm
      flush(6)
    end subroutine print_expdiag_timing

    subroutine apply_formal_local_field_exp(E_use, state_s, state_e)
      real(8), intent(in) :: E_use(3)
      integer, intent(in) :: state_s, state_e
      integer :: ispin_current, ifrag_use, i_local_use
      integer :: nw_use, nbf_use, io_use, iw_use, jw_use, istate_use
      integer :: global_idx_use, local_idx_use
      real(8) :: rdot_use, pc, ps
      complex(8), allocatable :: field_h(:,:), field_vec(:,:), field_tmp(:)
      complex(8), allocatable :: cw(:), tw(:), nwc(:)
      real(8), allocatable :: field_eval(:)

      if (.not. use_formal_wannier_h) return
      do ifrag_use = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local_use = ifrag_use - dg_frag%ifrag_start + 1
        if (i_local_use < 1 .or. i_local_use > size(dg_frag%dg_wannier_nkeep)) cycle
        nw_use = dg_frag%dg_wannier_nkeep(i_local_use)
        nbf_use = min(dg_frag%n_basis(ifrag_use, 1), dg_frag%nstate_frag, size(dg_frag%dg_wannier_basis_coef, 1))
        if (nw_use <= 0 .or. nbf_use <= 0) cycle
        allocate(field_h(nw_use,nw_use), field_vec(nw_use,nw_use), field_eval(nw_use), field_tmp(nw_use))
        allocate(cw(nw_use), tw(nw_use), nwc(nw_use))
        do ispin_current = 1, dg_frag%nspin
          field_h = (0.0d0, 0.0d0)
          do iw_use = 1, nw_use
            do jw_use = 1, nw_use
              field_h(jw_use,iw_use) = field_h(jw_use,iw_use) &
                - E_use(1) * dg_frag%dg_wannier_xi_local(1,jw_use,iw_use,ispin_current,i_local_use) &
                - E_use(2) * dg_frag%dg_wannier_xi_local(2,jw_use,iw_use,ispin_current,i_local_use) &
                - E_use(3) * dg_frag%dg_wannier_xi_local(3,jw_use,iw_use,ispin_current,i_local_use)
            end do
            rdot_use = E_use(1) * dg_frag%dg_wannier_ref_center(1,i_local_use) &
                     + E_use(2) * dg_frag%dg_wannier_ref_center(2,i_local_use) &
                     + E_use(3) * dg_frag%dg_wannier_ref_center(3,i_local_use)
            field_h(iw_use,iw_use) = field_h(iw_use,iw_use) - rdot_use
          end do
          field_vec = field_h
          call eigen_zheev(field_vec, field_eval, field_h)
          do istate_use = state_s, state_e
            if (istate_use < 1 .or. istate_use > size(dg_frag%coef, 2)) cycle
            cw = (0.0d0, 0.0d0)
            do iw_use = 1, nw_use
              do io_use = 1, nbf_use
                global_idx_use = dg_frag%index_basis(io_use, ifrag_use, ispin_current)
                if (global_idx_use < 1 .or. global_idx_use > dg_frag%n_mat_max) cycle
                local_idx_use = dg_frag%coef_global_to_local(global_idx_use, ispin_current)
                if (local_idx_use < 1 .or. local_idx_use > size(dg_frag%coef, 1)) cycle
                cw(iw_use) = cw(iw_use) + conjg(dg_frag%dg_wannier_basis_coef(io_use,iw_use,ispin_current,i_local_use)) * &
                  dg_frag%coef(local_idx_use,istate_use,ispin_current)
              end do
            end do
            tw = matmul(conjg(transpose(field_h)), cw)
            do iw_use = 1, nw_use
              pc = cos(field_eval(iw_use) * dt)
              ps = sin(field_eval(iw_use) * dt)
              tw(iw_use) = cmplx(pc, -ps, kind=8) * tw(iw_use)
            end do
            nwc = matmul(field_h, tw)
            do io_use = 1, nbf_use
              global_idx_use = dg_frag%index_basis(io_use, ifrag_use, ispin_current)
              if (global_idx_use < 1 .or. global_idx_use > dg_frag%n_mat_max) cycle
              local_idx_use = dg_frag%coef_global_to_local(global_idx_use, ispin_current)
              if (local_idx_use < 1 .or. local_idx_use > size(dg_frag%coef, 1)) cycle
              dg_frag%coef(local_idx_use,istate_use,ispin_current) = (0.0d0, 0.0d0)
              do iw_use = 1, nw_use
                dg_frag%coef(local_idx_use,istate_use,ispin_current) = &
                  dg_frag%coef(local_idx_use,istate_use,ispin_current) + &
                  dg_frag%dg_wannier_basis_coef(io_use,iw_use,ispin_current,i_local_use) * nwc(iw_use)
              end do
            end do
          end do
        end do
        deallocate(field_h, field_vec, field_eval, field_tmp, cw, tw, nwc)
      end do
    end subroutine apply_formal_local_field_exp

    subroutine apply_global_wannier_position_field_exp(E_use, state_s, state_e)
      real(8), intent(in) :: E_use(3)
      integer, intent(in) :: state_s, state_e
      integer :: ispin_current, ifrag_row, i_local_row
      integer :: nstate_blk, nwann, nbf
      integer :: istate, iwann, jwann, ibasis, global_row, local_row, state_col
      real(8) :: phase_c, phase_s
      complex(8), allocatable :: cw_local(:,:), cw_global(:,:), cw_next(:,:)
      complex(8), allocatable :: field_h(:,:), field_vec(:,:), tmp(:,:)
      complex(8), allocatable :: next_local(:,:)
      real(8), allocatable :: field_eval(:)

      nstate_blk = max(0, state_e - state_s + 1)
      if (nstate_blk <= 0) return
      if (.not. allocated(dg_frag%global_wannier_position)) return
      if (.not. allocated(dg_frag%global_wannier_coef)) return
      if (.not. allocated(dg_frag%coef_global_to_local)) return

      nwann = size(dg_frag%global_wannier_position, 2)
      if (nwann <= 0) return

      allocate(cw_local(nwann,nstate_blk), cw_global(nwann,nstate_blk), cw_next(nwann,nstate_blk))
      allocate(field_h(nwann,nwann), field_vec(nwann,nwann), tmp(nwann,nstate_blk), field_eval(nwann))
      do ispin_current = 1, dg_frag%nspin
        cw_local = (0.0d0, 0.0d0)
        do ifrag_row = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local_row = ifrag_row - dg_frag%ifrag_start + 1
          if (i_local_row < 1 .or. i_local_row > size(dg_frag%global_wannier_coef, 4)) cycle
          nbf = min(dg_frag%n_basis(ifrag_row, ispin_current), size(dg_frag%global_wannier_coef, 1))
          do ibasis = 1, nbf
            global_row = dg_frag%index_basis(ibasis, ifrag_row, ispin_current)
            if (global_row < 1 .or. global_row > dg_frag%n_mat_max) cycle
            local_row = dg_frag%coef_global_to_local(global_row, ispin_current)
            if (local_row < 1 .or. local_row > size(dg_frag%coef, 1)) cycle
            do iwann = 1, nwann
              do istate = 1, nstate_blk
                state_col = state_s + istate - 1
                if (state_col < 1 .or. state_col > size(dg_frag%coef, 2)) cycle
                cw_local(iwann,istate) = cw_local(iwann,istate) + &
                  conjg(dg_frag%global_wannier_coef(ibasis, iwann, ispin_current, i_local_row)) * &
                  dg_frag%coef(local_row,state_col,ispin_current)
              end do
            end do
          end do
        end do

        call comm_summation(cw_local, cw_global, nwann*nstate_blk, dg_frag%icomm)

        field_h = (0.0d0, 0.0d0)
        do iwann = 1, nwann
          do jwann = 1, nwann
            field_h(iwann,jwann) = - E_use(1) * dg_frag%global_wannier_position(1,iwann,jwann) &
                                  - E_use(2) * dg_frag%global_wannier_position(2,iwann,jwann) &
                                  - E_use(3) * dg_frag%global_wannier_position(3,iwann,jwann)
          end do
        end do
        field_vec = field_h
        call eigen_zheev(field_vec, field_eval, field_h)
        tmp = matmul(conjg(transpose(field_h)), cw_global)
        do iwann = 1, nwann
          phase_c = cos(field_eval(iwann) * dt)
          phase_s = sin(field_eval(iwann) * dt)
          tmp(iwann,1:nstate_blk) = cmplx(phase_c, -phase_s, kind=8) * tmp(iwann,1:nstate_blk)
        end do
        cw_next = matmul(field_h, tmp)

        allocate(next_local(size(dg_frag%coef,1),nstate_blk))
        next_local = (0.0d0, 0.0d0)
        do ifrag_row = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local_row = ifrag_row - dg_frag%ifrag_start + 1
          if (i_local_row < 1 .or. i_local_row > size(dg_frag%global_wannier_coef, 4)) cycle
          nbf = min(dg_frag%n_basis(ifrag_row, ispin_current), size(dg_frag%global_wannier_coef, 1))
          do ibasis = 1, nbf
            global_row = dg_frag%index_basis(ibasis, ifrag_row, ispin_current)
            if (global_row < 1 .or. global_row > dg_frag%n_mat_max) cycle
            local_row = dg_frag%coef_global_to_local(global_row, ispin_current)
            if (local_row < 1 .or. local_row > size(dg_frag%coef, 1)) cycle
            do istate = 1, nstate_blk
              do iwann = 1, nwann
                next_local(local_row,istate) = next_local(local_row,istate) + &
                  dg_frag%global_wannier_coef(ibasis, iwann, ispin_current, i_local_row) * cw_next(iwann,istate)
              end do
            end do
          end do
        end do

        do istate = 1, nstate_blk
          state_col = state_s + istate - 1
          if (state_col < 1 .or. state_col > size(dg_frag%coef, 2)) cycle
          dg_frag%coef(1:size(dg_frag%coef,1),state_col,ispin_current) = &
            next_local(1:size(dg_frag%coef,1),istate)
        end do
        deallocate(next_local)
      end do
      deallocate(cw_local, cw_global, cw_next, field_h, field_vec, tmp, field_eval)
    end subroutine apply_global_wannier_position_field_exp

    subroutine apply_global_mixed_position_field_exp(E_use, state_s, state_e)
      real(8), intent(in) :: E_use(3)
      integer, intent(in) :: state_s, state_e
      integer :: ispin_current, nstate_blk, nmix, nw, np, istate, imix
      real(8) :: phase_c, phase_s, tloc
      complex(8), allocatable :: cmix(:,:), field_h(:,:), field_vec(:,:), tmp(:,:)
      real(8), allocatable :: field_eval(:)

      if (.not. dg_frag%has_mixed_wannier_bpw_position) return
      if (.not. allocated(dg_frag%mixed_wannier_bpw_z)) return
      if (.not. allocated(dg_frag%mixed_wannier_bpw_pcoef)) return
      nstate_blk = max(0, state_e - state_s + 1)
      if (nstate_blk <= 0) return
      nw = dg_frag%mixed_wannier_bpw_nw
      np = dg_frag%mixed_wannier_bpw_np
      nmix = dg_frag%mixed_wannier_bpw_nmix
      if (nw <= 0 .or. np <= 0 .or. nmix /= nw + np) return

      allocate(cmix(nmix,nstate_blk), field_h(nmix,nmix), field_vec(nmix,nmix), tmp(nmix,nstate_blk), field_eval(nmix))
      do ispin_current = 1, dg_frag%nspin
        if (ispin_current > size(dg_frag%mixed_wannier_bpw_z, 4)) cycle
        tloc = get_wtime()
        call gather_global_mixed_coefficients(ispin_current, state_s, state_e, cmix)
        time_gather = time_gather + (get_wtime() - tloc)
        field_h(:, :) = (0.0d0, 0.0d0)
        field_h(:, :) = -E_use(1) * dg_frag%mixed_wannier_bpw_z(1,1:nmix,1:nmix,ispin_current) &
                       -E_use(2) * dg_frag%mixed_wannier_bpw_z(2,1:nmix,1:nmix,ispin_current) &
                       -E_use(3) * dg_frag%mixed_wannier_bpw_z(3,1:nmix,1:nmix,ispin_current)
        field_vec = field_h
        tloc = get_wtime()
        call eigen_zheev(field_vec, field_eval, field_h)
        time_field_diag = time_field_diag + (get_wtime() - tloc)
        tloc = get_wtime()
        tmp = matmul(conjg(transpose(field_h)), cmix)
        do imix = 1, nmix
          phase_c = cos(field_eval(imix) * dt)
          phase_s = sin(field_eval(imix) * dt)
          tmp(imix,1:nstate_blk) = cmplx(phase_c, -phase_s, kind=8) * tmp(imix,1:nstate_blk)
        end do
        cmix = matmul(field_h, tmp)
        time_field_matmul = time_field_matmul + (get_wtime() - tloc)
        tloc = get_wtime()
        call scatter_global_mixed_coefficients(ispin_current, state_s, state_e, cmix)
        time_scatter = time_scatter + (get_wtime() - tloc)
      end do
      deallocate(cmix, field_h, field_vec, tmp, field_eval)
    end subroutine apply_global_mixed_position_field_exp

    subroutine apply_global_mixed_phase_exp(state_s, state_e)
      integer, intent(in) :: state_s, state_e
      integer :: ispin_current, nstate_blk, nmix, imix, istate
      real(8) :: phase_c, phase_s, tloc
      complex(8), allocatable :: cmix(:,:)

      if (.not. dg_frag%has_mixed_wannier_bpw_position) return
      if (.not. allocated(dg_frag%mixed_wannier_bpw_eval)) return
      nstate_blk = max(0, state_e - state_s + 1)
      if (nstate_blk <= 0) return
      nmix = dg_frag%mixed_wannier_bpw_nmix
      if (nmix <= 0) return

      allocate(cmix(nmix,nstate_blk))
      do ispin_current = 1, dg_frag%nspin
        if (ispin_current > size(dg_frag%mixed_wannier_bpw_eval, 2)) cycle
        tloc = get_wtime()
        call gather_global_mixed_coefficients(ispin_current, state_s, state_e, cmix)
        time_gather = time_gather + (get_wtime() - tloc)
        do istate = 1, nstate_blk
          do imix = 1, nmix
            phase_c = cos(dg_frag%mixed_wannier_bpw_eval(imix,ispin_current) * dt)
            phase_s = sin(dg_frag%mixed_wannier_bpw_eval(imix,ispin_current) * dt)
            cmix(imix,istate) = cmplx(phase_c, -phase_s, kind=8) * cmix(imix,istate)
          end do
        end do
        tloc = get_wtime()
        call scatter_global_mixed_coefficients(ispin_current, state_s, state_e, cmix)
        time_scatter = time_scatter + (get_wtime() - tloc)
      end do
      deallocate(cmix)
    end subroutine apply_global_mixed_phase_exp

    subroutine gather_global_mixed_coefficients(ispin_current, state_s, state_e, cmix)
      integer, intent(in) :: ispin_current, state_s, state_e
      complex(8), intent(out) :: cmix(:,:)
      integer :: ifrag_row, i_local_row, nbf, ibasis, iwann, istate, state_col
      integer :: global_row, local_row, nstate_blk, nw, np
      complex(8), allocatable :: cw_local(:,:), cw_global(:,:)

      nstate_blk = max(0, state_e - state_s + 1)
      nw = dg_frag%mixed_wannier_bpw_nw
      np = dg_frag%mixed_wannier_bpw_np
      cmix(:, :) = (0.0d0, 0.0d0)
      if (nstate_blk <= 0 .or. nw <= 0) return
      allocate(cw_local(nw,nstate_blk), cw_global(nw,nstate_blk))
      cw_local(:, :) = (0.0d0, 0.0d0)
      do ifrag_row = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local_row = ifrag_row - dg_frag%ifrag_start + 1
        if (i_local_row < 1 .or. i_local_row > size(dg_frag%global_wannier_coef, 4)) cycle
        nbf = min(dg_frag%n_basis(ifrag_row, ispin_current), size(dg_frag%global_wannier_coef, 1))
        do ibasis = 1, nbf
          global_row = dg_frag%index_basis(ibasis, ifrag_row, ispin_current)
          if (global_row < 1 .or. global_row > dg_frag%n_mat_max) cycle
          local_row = dg_frag%coef_global_to_local(global_row, ispin_current)
          if (local_row < 1 .or. local_row > size(dg_frag%coef, 1)) cycle
          do iwann = 1, nw
            do istate = 1, nstate_blk
              state_col = state_s + istate - 1
              if (state_col < 1 .or. state_col > size(dg_frag%coef, 2)) cycle
              cw_local(iwann,istate) = cw_local(iwann,istate) + &
                conjg(dg_frag%global_wannier_coef(ibasis, iwann, ispin_current, i_local_row)) * &
                dg_frag%coef(local_row,state_col,ispin_current)
            end do
          end do
        end do
      end do
      call comm_summation(cw_local, cw_global, nw*nstate_blk, dg_frag%icomm)
      cmix(1:nw,1:nstate_blk) = cw_global(1:nw,1:nstate_blk)
      if (np > 0 .and. allocated(dg_frag%mixed_wannier_bpw_pcoef)) then
        cmix(nw+1:nw+np,1:nstate_blk) = dg_frag%mixed_wannier_bpw_pcoef(1:np,state_s:state_e,ispin_current)
      end if
      deallocate(cw_local, cw_global)
    end subroutine gather_global_mixed_coefficients

    subroutine scatter_global_mixed_coefficients(ispin_current, state_s, state_e, cmix)
      integer, intent(in) :: ispin_current, state_s, state_e
      complex(8), intent(in) :: cmix(:,:)
      integer :: ifrag_row, i_local_row, nbf, ibasis, iwann, istate, state_col
      integer :: global_row, local_row, nstate_blk, nw, np
      complex(8), allocatable :: next_local(:,:)

      nstate_blk = max(0, state_e - state_s + 1)
      nw = dg_frag%mixed_wannier_bpw_nw
      np = dg_frag%mixed_wannier_bpw_np
      if (nstate_blk <= 0 .or. nw <= 0) return
      allocate(next_local(size(dg_frag%coef,1),nstate_blk))
      next_local(:, :) = (0.0d0, 0.0d0)
      do ifrag_row = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local_row = ifrag_row - dg_frag%ifrag_start + 1
        if (i_local_row < 1 .or. i_local_row > size(dg_frag%global_wannier_coef, 4)) cycle
        nbf = min(dg_frag%n_basis(ifrag_row, ispin_current), size(dg_frag%global_wannier_coef, 1))
        do ibasis = 1, nbf
          global_row = dg_frag%index_basis(ibasis, ifrag_row, ispin_current)
          if (global_row < 1 .or. global_row > dg_frag%n_mat_max) cycle
          local_row = dg_frag%coef_global_to_local(global_row, ispin_current)
          if (local_row < 1 .or. local_row > size(dg_frag%coef, 1)) cycle
          do istate = 1, nstate_blk
            do iwann = 1, nw
              next_local(local_row,istate) = next_local(local_row,istate) + &
                dg_frag%global_wannier_coef(ibasis, iwann, ispin_current, i_local_row) * cmix(iwann,istate)
            end do
          end do
        end do
      end do
      do istate = 1, nstate_blk
        state_col = state_s + istate - 1
        if (state_col < 1 .or. state_col > size(dg_frag%coef, 2)) cycle
        dg_frag%coef(1:size(dg_frag%coef,1),state_col,ispin_current) = next_local(1:size(dg_frag%coef,1),istate)
      end do
      if (np > 0 .and. allocated(dg_frag%mixed_wannier_bpw_pcoef)) then
        dg_frag%mixed_wannier_bpw_pcoef(1:np,state_s:state_e,ispin_current) = cmix(nw+1:nw+np,1:nstate_blk)
      end if
      deallocate(next_local)
    end subroutine scatter_global_mixed_coefficients

    subroutine apply_global_flux_eigen_exp(state_s, state_e)
      integer, intent(in) :: state_s, state_e
      integer :: ispin_current, ifrag_row, i_local_row
      integer :: nstate_blk, neig, nwann, nbf
      integer :: istate, eig, iwann, ibasis, global_row, local_row, state_col
      real(8) :: phase_c, phase_s
      real(8) :: tloc
      complex(8), allocatable :: amp_local(:,:), amp_global(:,:), phi_row(:)
      complex(8), allocatable :: next_local(:,:)

      nstate_blk = max(0, state_e - state_s + 1)
      if (nstate_blk <= 0) return
      if (.not. allocated(dg_frag%global_wannier_flux_evec)) return
      if (.not. allocated(dg_frag%global_wannier_flux_eval)) return
      if (.not. allocated(dg_frag%global_wannier_coef)) return
      if (.not. allocated(dg_frag%coef_global_to_local)) return

      nwann = size(dg_frag%global_wannier_flux_evec, 1)
      neig = min(size(dg_frag%global_wannier_flux_evec, 2), &
                 size(dg_frag%global_wannier_flux_eval, 1))
      if (nwann <= 0 .or. neig <= 0) return

      allocate(amp_local(neig,nstate_blk), amp_global(neig,nstate_blk))
      allocate(phi_row(neig))
      do ispin_current = 1, dg_frag%nspin
        if (ispin_current > size(dg_frag%global_wannier_flux_eval, 2)) cycle
        amp_local = (0.0d0, 0.0d0)
        tloc = get_wtime()
        do ifrag_row = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local_row = ifrag_row - dg_frag%ifrag_start + 1
          if (i_local_row < 1 .or. i_local_row > size(dg_frag%global_wannier_coef, 4)) cycle
          nbf = min(dg_frag%n_basis(ifrag_row, ispin_current), size(dg_frag%global_wannier_coef, 1))
          do ibasis = 1, nbf
            global_row = dg_frag%index_basis(ibasis, ifrag_row, ispin_current)
            if (global_row < 1 .or. global_row > dg_frag%n_mat_max) cycle
            local_row = dg_frag%coef_global_to_local(global_row, ispin_current)
            if (local_row < 1 .or. local_row > size(dg_frag%coef, 1)) cycle
            phi_row = (0.0d0, 0.0d0)
            do eig = 1, neig
              do iwann = 1, nwann
                phi_row(eig) = phi_row(eig) + &
                  dg_frag%global_wannier_coef(ibasis, iwann, ispin_current, i_local_row) * &
                  dg_frag%global_wannier_flux_evec(iwann, eig)
              end do
            end do
            do istate = 1, nstate_blk
              state_col = state_s + istate - 1
              if (state_col < 1 .or. state_col > size(dg_frag%coef, 2)) cycle
              do eig = 1, neig
                amp_local(eig,istate) = amp_local(eig,istate) + &
                  conjg(phi_row(eig)) * dg_frag%coef(local_row,state_col,ispin_current)
              end do
            end do
          end do
        end do
        time_flux_project = time_flux_project + (get_wtime() - tloc)

        tloc = get_wtime()
        call comm_summation(amp_local, amp_global, neig*nstate_blk, dg_frag%icomm)
        time_flux_comm = time_flux_comm + (get_wtime() - tloc)
        do istate = 1, nstate_blk
          do eig = 1, neig
            phase_c = cos(dg_frag%global_wannier_flux_eval(eig,ispin_current) * dt)
            phase_s = sin(dg_frag%global_wannier_flux_eval(eig,ispin_current) * dt)
            amp_global(eig,istate) = cmplx(phase_c, -phase_s, kind=8) * amp_global(eig,istate)
          end do
        end do

        allocate(next_local(size(dg_frag%coef,1),nstate_blk))
        next_local = (0.0d0, 0.0d0)
        tloc = get_wtime()
        do ifrag_row = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local_row = ifrag_row - dg_frag%ifrag_start + 1
          if (i_local_row < 1 .or. i_local_row > size(dg_frag%global_wannier_coef, 4)) cycle
          nbf = min(dg_frag%n_basis(ifrag_row, ispin_current), size(dg_frag%global_wannier_coef, 1))
          do ibasis = 1, nbf
            global_row = dg_frag%index_basis(ibasis, ifrag_row, ispin_current)
            if (global_row < 1 .or. global_row > dg_frag%n_mat_max) cycle
            local_row = dg_frag%coef_global_to_local(global_row, ispin_current)
            if (local_row < 1 .or. local_row > size(dg_frag%coef, 1)) cycle
            phi_row = (0.0d0, 0.0d0)
            do eig = 1, neig
              do iwann = 1, nwann
                phi_row(eig) = phi_row(eig) + &
                  dg_frag%global_wannier_coef(ibasis, iwann, ispin_current, i_local_row) * &
                  dg_frag%global_wannier_flux_evec(iwann, eig)
              end do
            end do
            do istate = 1, nstate_blk
              do eig = 1, neig
                next_local(local_row,istate) = next_local(local_row,istate) + &
                  phi_row(eig) * amp_global(eig,istate)
              end do
            end do
          end do
        end do
        time_flux_scatter = time_flux_scatter + (get_wtime() - tloc)

        do istate = 1, nstate_blk
          state_col = state_s + istate - 1
          if (state_col < 1 .or. state_col > size(dg_frag%coef, 2)) cycle
          dg_frag%coef(1:size(dg_frag%coef,1),state_col,ispin_current) = &
            next_local(1:size(dg_frag%coef,1),istate)
        end do
        deallocate(next_local)
      end do
      deallocate(amp_local, amp_global, phi_row)
    end subroutine apply_global_flux_eigen_exp

    subroutine apply_xi_flux_split_half(E_use, dt_half, state_s, state_e)
      real(8), intent(in) :: E_use(3), dt_half
      integer, intent(in) :: state_s, state_e
      integer :: ispin_current, iblk_idx, iblk, ifrag_row, ifrag_col
      integer :: nstate_blk, state_col, io, jo, istate
      integer :: nrow, ncol, row_gid, col_gid, row_local, col_pos
      integer :: nfetched
      real(8) :: val
      complex(8), allocatable :: fetched(:,:), delta(:,:)
      integer, allocatable :: fetched_ids(:)

      if (sum(abs(E_use(1:3))) <= 1.0d-30) return
      if (abs(dt_half) <= 0.0d0) return
      nstate_blk = max(0, state_e - state_s + 1)
      if (nstate_blk <= 0) return

      if (.not. dg_frag%buffer_wannier_xi_flux_available .or. &
          .not. allocated(dg_frag%buffer_wannier_xi_flux_blocks) .or. &
          .not. allocated(dg_frag%buffer_wannier_xi_flux_local_block_ids)) then
        stop "DG expdiag requires neighbor xi_flux blocks"
      end if

      do ispin_current = 1, dg_frag%nspin
        call exchange_xi_flux_neighbor_coef(ispin_current, state_s, state_e, fetched_ids, fetched, nfetched)
        allocate(delta(size(dg_frag%coef,1),nstate_blk))
        delta = (0.0d0, 0.0d0)
        do iblk_idx = 1, size(dg_frag%buffer_wannier_xi_flux_local_block_ids)
          iblk = dg_frag%buffer_wannier_xi_flux_local_block_ids(iblk_idx)
          if (iblk < 1 .or. iblk > size(dg_frag%buffer_wannier_xi_flux_blocks)) cycle
          ifrag_row = dg_frag%buffer_wannier_xi_flux_blocks(iblk)%ifrag_row
          ifrag_col = dg_frag%buffer_wannier_xi_flux_blocks(iblk)%ifrag_col
          if (ifrag_row == ifrag_col) cycle
          if (.not. allocated(dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val)) cycle
          nrow = min(dg_frag%n_basis(ifrag_row, ispin_current), &
                     size(dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val, 2))
          ncol = min(dg_frag%n_basis(ifrag_col, ispin_current), &
                     size(dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val, 3))
          do io = 1, nrow
            row_gid = dg_frag%index_basis(io, ifrag_row, ispin_current)
            if (row_gid < 1 .or. row_gid > dg_frag%n_mat_max) cycle
            row_local = dg_frag%coef_global_to_local(row_gid, ispin_current)
            if (row_local <= 0 .or. row_local > size(dg_frag%coef, 1)) cycle
            do jo = 1, ncol
              col_gid = dg_frag%index_basis(jo, ifrag_col, ispin_current)
              col_pos = find_needed_row_pos(col_gid, fetched_ids, nfetched)
              if (col_pos <= 0) cycle
              val = E_use(1) * dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val(1,io,jo,ispin_current) &
                  + E_use(2) * dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val(2,io,jo,ispin_current) &
                  + E_use(3) * dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val(3,io,jo,ispin_current)
              if (abs(val) <= 0.0d0) cycle
              do istate = 1, nstate_blk
                delta(row_local,istate) = delta(row_local,istate) - zi * dt_half * val * fetched(col_pos,istate)
              end do
            end do
          end do
        end do
        do istate = 1, nstate_blk
          state_col = state_s + istate - 1
          if (state_col < 1 .or. state_col > size(dg_frag%coef, 2)) cycle
          dg_frag%coef(1:size(dg_frag%coef,1),state_col,ispin_current) = &
            dg_frag%coef(1:size(dg_frag%coef,1),state_col,ispin_current) + &
            delta(1:size(dg_frag%coef,1),istate)
        end do
        if (allocated(fetched_ids)) deallocate(fetched_ids)
        if (allocated(fetched)) deallocate(fetched)
        deallocate(delta)
      end do
    end subroutine apply_xi_flux_split_half

    subroutine exchange_xi_flux_neighbor_coef(ispin_use, state_s, state_e, fetched_ids, fetched_coef, nfetched)
      use mpi, only: MPI_Isend, MPI_Irecv, MPI_Waitall, MPI_INTEGER, MPI_DOUBLE_COMPLEX, &
                     MPI_SUCCESS, MPI_STATUSES_IGNORE
      integer, intent(in) :: ispin_use, state_s, state_e
      integer, allocatable, intent(out) :: fetched_ids(:)
      complex(8), allocatable, intent(out) :: fetched_coef(:,:)
      integer, intent(out) :: nfetched

      integer :: nstate_blk, peer, npeer, iblk, iblk_idx, ifrag_row, ifrag_col
      integer :: i, ist, local_idx, gid, ierr, nreq, total_recv, total_send
      integer :: data_pos, row_pos, peer_owner
      integer, allocatable :: peers(:), send_counts(:), recv_counts(:), send_displs(:), recv_displs(:)
      integer, allocatable :: requests(:), send_ids(:), recv_ids(:)
      complex(8), allocatable :: send_buf(:), recv_buf(:)

      nfetched = 0
      nstate_blk = max(0, state_e - state_s + 1)
      allocate(fetched_ids(1), fetched_coef(1,max(1,nstate_blk)))
      fetched_ids(:) = 0
      fetched_coef(:,:) = (0.0d0, 0.0d0)
      if (nstate_blk <= 0) return
      if (.not. allocated(dg_frag%coef)) return
      if (.not. allocated(dg_frag%coef_global_to_local)) return
      if (.not. allocated(dg_frag%coef_owner)) return
      if (.not. allocated(dg_frag%local_coef_count)) return
      if (.not. allocated(dg_frag%local_coef_global_ids)) return
      if (.not. allocated(dg_frag%id_array)) return

      allocate(peers(max(1, dg_frag%isize)))
      peers(:) = -1
      npeer = 0
      if (allocated(dg_frag%buffer_wannier_xi_flux_blocks) .and. &
          allocated(dg_frag%buffer_wannier_xi_flux_local_block_ids)) then
        do iblk_idx = 1, size(dg_frag%buffer_wannier_xi_flux_local_block_ids)
          iblk = dg_frag%buffer_wannier_xi_flux_local_block_ids(iblk_idx)
          if (iblk < 1 .or. iblk > size(dg_frag%buffer_wannier_xi_flux_blocks)) cycle
          ifrag_row = dg_frag%buffer_wannier_xi_flux_blocks(iblk)%ifrag_row
          ifrag_col = dg_frag%buffer_wannier_xi_flux_blocks(iblk)%ifrag_col
          if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) cycle
          if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
          if (ifrag_row == ifrag_col) cycle
          if (fragment_is_local(ifrag_row)) then
            peer_owner = dg_frag%id_array(ifrag_col)
            call add_exchange_peer(peer_owner, peers, npeer)
          end if
          if (fragment_is_local(ifrag_col)) then
            peer_owner = dg_frag%id_array(ifrag_row)
            call add_exchange_peer(peer_owner, peers, npeer)
          end if
        end do
      end if

      if (npeer <= 0) then
        deallocate(peers)
        return
      end if

      total_send = 0
      if (ispin_use >= 1 .and. ispin_use <= size(dg_frag%local_coef_count)) &
        total_send = max(0, dg_frag%local_coef_count(ispin_use))

      allocate(send_counts(npeer), recv_counts(npeer), send_displs(npeer), recv_displs(npeer))
      allocate(requests(max(1, 2*npeer)))
      send_counts(:) = total_send
      recv_counts(:) = 0
      nreq = 0
      do i = 1, npeer
        peer = peers(i)
        nreq = nreq + 1
        call MPI_Irecv(recv_counts(i), 1, MPI_INTEGER, peer, 8521, dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux count recv failed"
      end do
      do i = 1, npeer
        peer = peers(i)
        nreq = nreq + 1
        call MPI_Isend(send_counts(i), 1, MPI_INTEGER, peer, 8521, dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux count send failed"
      end do
      call MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE, ierr)
      if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux count wait failed"

      send_displs(1) = 0
      recv_displs(1) = 0
      do i = 2, npeer
        send_displs(i) = send_displs(i-1) + send_counts(i-1)
        recv_displs(i) = recv_displs(i-1) + recv_counts(i-1)
      end do
      total_recv = recv_displs(npeer) + recv_counts(npeer)

      allocate(send_ids(max(1,total_send)), recv_ids(max(1,total_recv)))
      send_ids(:) = 0
      recv_ids(:) = 0
      if (total_send > 0) send_ids(1:total_send) = dg_frag%local_coef_global_ids(1:total_send, ispin_use)

      nreq = 0
      do i = 1, npeer
        peer = peers(i)
        if (recv_counts(i) <= 0) cycle
        nreq = nreq + 1
        call MPI_Irecv(recv_ids(recv_displs(i)+1), recv_counts(i), MPI_INTEGER, peer, 8522, dg_frag%icomm, &
                       requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux row recv failed"
      end do
      do i = 1, npeer
        peer = peers(i)
        if (send_counts(i) <= 0) cycle
        nreq = nreq + 1
        call MPI_Isend(send_ids(1), send_counts(i), MPI_INTEGER, peer, 8522, dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux row send failed"
      end do
      if (nreq > 0) then
        call MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE, ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux row wait failed"
      end if

      allocate(send_buf(max(1,total_send*nstate_blk)), recv_buf(max(1,total_recv*nstate_blk)))
      send_buf(:) = (0.0d0, 0.0d0)
      recv_buf(:) = (0.0d0, 0.0d0)
      do ist = 1, nstate_blk
        do i = 1, total_send
          gid = send_ids(i)
          if (gid < 1 .or. gid > size(dg_frag%coef_global_to_local, 1)) cycle
          local_idx = dg_frag%coef_global_to_local(gid, ispin_use)
          if (local_idx < 1 .or. local_idx > size(dg_frag%coef, 1)) cycle
          if (state_s + ist - 1 < 1 .or. state_s + ist - 1 > size(dg_frag%coef, 2)) cycle
          send_buf((ist-1)*total_send + i) = dg_frag%coef(local_idx, state_s + ist - 1, ispin_use)
        end do
      end do

      nreq = 0
      do i = 1, npeer
        peer = peers(i)
        if (recv_counts(i) <= 0) cycle
        nreq = nreq + 1
        call MPI_Irecv(recv_buf(recv_displs(i)*nstate_blk + 1), recv_counts(i)*nstate_blk, &
                       MPI_DOUBLE_COMPLEX, peer, 8523, dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux data recv failed"
      end do
      do i = 1, npeer
        peer = peers(i)
        if (send_counts(i) <= 0) cycle
        nreq = nreq + 1
        call MPI_Isend(send_buf(1), send_counts(i)*nstate_blk, MPI_DOUBLE_COMPLEX, peer, 8523, &
                       dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux data send failed"
      end do
      if (nreq > 0) then
        call MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE, ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux data wait failed"
      end if

      deallocate(fetched_ids, fetched_coef)
      nfetched = total_recv
      allocate(fetched_ids(max(1,nfetched)), fetched_coef(max(1,nfetched),nstate_blk))
      fetched_ids(:) = 0
      fetched_coef(:,:) = (0.0d0, 0.0d0)
      if (nfetched > 0) then
        fetched_ids(1:nfetched) = recv_ids(1:nfetched)
        do peer = 1, npeer
          do ist = 1, nstate_blk
            do i = 1, recv_counts(peer)
              row_pos = recv_displs(peer) + i
              data_pos = recv_displs(peer) * nstate_blk + (ist - 1) * recv_counts(peer) + i
              fetched_coef(row_pos,ist) = recv_buf(data_pos)
            end do
          end do
        end do
      end if

      deallocate(peers, send_counts, recv_counts, send_displs, recv_displs, requests)
      deallocate(send_ids, recv_ids, send_buf, recv_buf)
    end subroutine exchange_xi_flux_neighbor_coef

    logical function fragment_is_local(ifrag_use) result(is_local)
      integer, intent(in) :: ifrag_use
      is_local = (ifrag_use >= dg_frag%ifrag_start .and. ifrag_use <= dg_frag%ifrag_end)
    end function fragment_is_local

    subroutine add_exchange_peer(peer, peers, npeer)
      integer, intent(in) :: peer
      integer, intent(inout) :: peers(:)
      integer, intent(inout) :: npeer
      integer :: i
      if (peer < 0 .or. peer >= dg_frag%isize .or. peer == dg_frag%id) return
      do i = 1, npeer
        if (peers(i) == peer) return
      end do
      if (npeer >= size(peers)) return
      npeer = npeer + 1
      peers(npeer) = peer
    end subroutine add_exchange_peer

    integer function find_needed_row_pos(gid, row_ids, nneeded) result(pos)
      integer, intent(in) :: gid, nneeded
      integer, intent(in) :: row_ids(:)
      integer :: i

      pos = 0
      if (gid <= 0) return
      do i = 1, nneeded
        if (row_ids(i) == gid) then
          pos = i
          return
        end if
      end do
    end function find_needed_row_pos
  end subroutine time_evolution_expdiag
