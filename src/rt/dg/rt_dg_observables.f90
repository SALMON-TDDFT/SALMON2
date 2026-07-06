  subroutine calculate_observables(dg_frag, system, mg, stencil, ppg, rt, itt, Vh, Vxc, Vpsl, rho)
    use structures
    use rt_dg_fragment_types, only: s_dg_fragment_rt
    use rt_dg_fragment_ops, only: calculate_macroscopic_current_dg, calculate_nonlocal_current_dg, &
                                  calculate_local_wannier_polarization_dg, ensure_gradient_basis_cache
    use rt_dg_plane_wave, only: apply_wpw_reduced_pz_to_production, map_global_to_phi_box_coord_pw
    use communication, only: comm_summation
    use unusedvar_mod, only: salmon_unusedvar
    use salmon_global, only: yn_dg_length_gauge, yn_dg_mixed_z_local_pz_writeback_total, &
      yn_dg_mixed_z_local_current_writeback_total
    use misc_routines, only: get_wtime
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_pp_grid),        intent(in)    :: ppg
    type(s_rt),             intent(inout) :: rt
    integer,                intent(in)    :: itt
    type(s_scalar),         intent(in)    :: Vh, Vpsl
    type(s_scalar),         intent(in)    :: Vxc(system%nspin)
    type(s_scalar),         intent(in), optional :: rho
    integer :: ispin
    real(8) :: current_raw(3), current_nl_raw(3), polarization_raw(3), polarization_density(3)
    real(8) :: nelec_ref, ne_density
    real(8) :: pz_prod_raw_before, pz_prod_output_before, pz_prod_density
    real(8) :: pz_prod_output_after, pz_prod_output_delta, pz_prod_replace_tol
    real(8) :: pz_prod_raw_tol, pz_total_candidate_raw
    real(8) :: wwfull_pz_diag, wwfull_pz_offdiag, wwfull_pz_total, wwfull_pz_diff
    real(8) :: pz_total_candidate, pz_total_wp_combined, pz_total_pp, pz_total_diff
    real(8) :: pz_total_density, pz_total_output_after, pz_total_tol
    real(8) :: current_total_diff_ww(3), current_total_diff_full(3)
    real(8) :: current_para_candidate(3), current_nl_candidate(3), current_dia_candidate(3)
    real(8) :: current_para_total_candidate(3)
    real(8) :: current_candidate_total(3), current_candidate_diff(3)
    real(8) :: current_before_writeback(3), current_after_minus_candidate(3)
    real(8) :: current_after_minus_before(3)
    real(8) :: current_total_residual_norm, current_total_full_norm
    real(8) :: current_total_wp_norm, current_total_pp_norm, current_total_p_ratio
    real(8) :: current_total_tol
    real(8) :: prop_obs_current_prod_norm, prop_obs_current_candidate_norm, prop_obs_current_diff_norm
    real(8) :: prop_obs_rho_prod_int, prop_obs_rho_candidate_int, prop_obs_rho_diff_int
    real(8) :: prop_obs_rho_max_abs_diff, prop_obs_rho_rms_diff
    real(8) :: prop_obs_pz_prod, prop_obs_pz_candidate, prop_obs_pz_diff
    real(8) :: prop_obs_pz_prod_prev, prop_obs_pz_diff_to_current, prop_obs_pz_diff_to_prev
    real(8) :: prop_obs_candidate_coef_norm, prop_obs_candidate_occ_trace
    real(8) :: prop_obs_payload_current(3), prop_obs_payload_mom_current(3), prop_obs_current_candidate_vec(3)
    real(8) :: prop_obs_payload_mom_self_current(3), prop_obs_payload_mom_cross_current(3)
    real(8) :: prop_obs_current_para_prod_norm
    real(8) :: prop_obs_current_nl_prod_norm, prop_obs_current_dia_prod_norm
    real(8) :: prop_obs_current_fragbasis_vec(3)
    real(8) :: prop_obs_current_mom_block_counts(3)
    real(8) :: prop_obs_payload_pz, prop_obs_occ
    complex(8) :: prop_obs_zterm
    logical :: enable_prop_obs_series, prop_obs_field_on, prop_obs_bad
    logical :: prop_obs_payload_valid
    logical :: prop_obs_payload_current_ready
    logical :: prop_obs_prev_pz_available
    logical :: prop_obs_pz_available, prop_obs_current_available
    character(24) :: prop_obs_field_kind, prop_obs_candidate_source, prop_obs_payload_source
    character(32) :: prop_obs_candidate_contraction_source
    character(32) :: prop_obs_payload_coef_time_tag
    real(8) :: zww_src_sum_diag_expdiag, zww_src_sum_diag_actual, zww_src_sum_center
    real(8) :: zww_src_sum_cell_shift, zww_src_contrib_expdiag, zww_src_contrib_actual
    real(8) :: zww_src_contrib_diff, zww_src_contrib_diag_actual, zww_src_contrib_offdiag_actual
    real(8) :: zww_src_diag_contrib_diff, zww_src_offdiag_residual
    logical :: wwfull_replacement_ready
    logical :: zww_src_same_gid_order, zww_src_same_owner_mapping, zww_src_same_center_source
    logical :: zww_src_same_cell_shift_source, zww_src_offdiag_explains_diff, zww_src_bad
    character(32) :: pz_prod_block_reason, pz_total_block_reason, current_total_block_reason
    character(48) :: current_para_candidate_source
    character(32), save :: pol_trace_env = ''
    character(32), save :: pw_obs_skip_env = ''
    character(32), save :: pz_obs_hook_env = ''
    character(32), save :: current_total_hook_env = ''
    character(32), save :: current_total_writeback_env = ''
    character(32), save :: prop_obs_series_env = ''
    character(32), save :: prop_obs_series_source_env = ''
    character(32), save :: prop_obs_series_local_mixedz_env = ''
    character(32), save :: prop_obs_series_current_detail_env = ''
    character(32), save :: perf_count_env = ''
    logical, save :: pol_trace_initialized = .false.
    logical, save :: enable_pol_trace = .false.
    logical, save :: pw_obs_skip_initialized = .false.
    logical, save :: enable_pw_obs_skip = .false.
    logical, save :: pz_obs_hook_initialized = .false.
    logical, save :: enable_pz_total_hook = .false.
    logical, save :: enable_pz_total_writeback = .false.
    logical, save :: current_total_hook_initialized = .false.
    logical, save :: enable_current_total_hook = .false.
    logical, save :: enable_current_total_writeback = .false.
    logical, save :: prop_obs_series_initialized = .false.
    logical, save :: prop_obs_series_enabled = .false.
    logical, save :: enable_prop_obs_series_current_detail = .false.
    logical, save :: perf_count_initialized = .false.
    logical, save :: enable_perf_count = .false.
    logical, save :: prop_obs_prev_pz_valid = .false.
    character(24), save :: prop_obs_series_source_value = 'prod_reference'
    real(8), save :: prop_obs_prev_pz_value = 0.0d0
    integer, save :: prop_obs_prev_pz_step = -1
    integer :: env_status, env_len
    integer :: prop_obs_iw, prop_obs_jw, prop_obs_ist, prop_obs_ispin
    integer :: prop_obs_nw, prop_obs_nstate, prop_obs_nspin
    real(8) :: perf_obs_t0, perf_series_t0, perf_current_t0, perf_section_t0
    real(8) :: perf_pz_t0, perf_reduce_t0, perf_payload_current_t0
    real(8) :: perf_pz_section_t0
    real(8) :: stability_rho, stability_pz, stability_current, stability_norm, stability_energy
    logical :: stability_bad
    logical :: perf_final_step

    call salmon_unusedvar(Vh)
    call salmon_unusedvar(Vxc)
    call salmon_unusedvar(Vpsl)
    if (.not. perf_count_initialized) then
      call get_environment_variable('SALMON_DG_MIXED_Z_PERF_COUNT', perf_count_env, &
        length=env_len, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(perf_count_env(1:env_len))))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enable_perf_count = .true.
        end select
      end if
      perf_count_initialized = .true.
    end if
    dg_frag%mixed_z_perf_count_enabled = enable_perf_count
    perf_obs_t0 = 0.0d0
    if (enable_perf_count) then
      perf_obs_t0 = get_wtime()
      dg_frag%mixed_z_perf_nstep = dg_frag%mixed_z_perf_nstep + 1_8
    end if
    pz_total_candidate = 0.0d0
    pz_total_wp_combined = 0.0d0
    pz_total_pp = 0.0d0
    pz_total_diff = 0.0d0
    pz_total_tol = huge(1.0d0)
    pz_prod_raw_before = 0.0d0
    prop_obs_pz_available = .false.
    prop_obs_current_available = .false.
    prop_obs_bad = .false.
    prop_obs_payload_current(:) = 0.0d0
    prop_obs_payload_mom_current(:) = 0.0d0
    prop_obs_payload_mom_self_current(:) = 0.0d0
    prop_obs_payload_mom_cross_current(:) = 0.0d0
    prop_obs_current_candidate_vec(:) = 0.0d0
    prop_obs_payload_current_ready = .false.
    prop_obs_current_para_prod_norm = 0.0d0
    prop_obs_current_nl_prod_norm = 0.0d0
    prop_obs_current_dia_prod_norm = 0.0d0
    prop_obs_current_fragbasis_vec(:) = 0.0d0
    prop_obs_current_mom_block_counts(:) = 0.0d0
    ! Fixed PROP-OBS candidate_source values:
    !   prod_reference, local_mixedz, local_mixedz_writeback
    prop_obs_candidate_source = prop_obs_series_source_value

    if ((dg_frag%use_plane_wave_basis .or. allocated(dg_frag%coef_pw)) .and. &
        .not. dg_frag%has_mixed_wannier_bpw_position) then
      call ensure_mixed_wannier_bpw_position(dg_frag)
    end if
    if ((dg_frag%use_plane_wave_basis .or. allocated(dg_frag%coef_pw)) .and. &
        .not. dg_frag%has_mixed_wannier_bpw_position) then
      if (.not. pw_obs_skip_initialized) then
        call get_environment_variable('SALMON_DG_SKIP_UNSUPPORTED_PW_OBSERVABLE', pw_obs_skip_env, &
          length=env_len, status=env_status)
        if (env_status == 0) then
          select case(trim(adjustl(pw_obs_skip_env(1:env_len))))
          case('1','y','Y','yes','YES','true','TRUE','on','ON')
            enable_pw_obs_skip = .true.
          end select
        end if
        pw_obs_skip_initialized = .true.
      end if
      if (enable_pw_obs_skip) then
        dg_frag%current(:) = 0.0d0
        dg_frag%current_nl(:) = 0.0d0
        dg_frag%current_para(:) = 0.0d0
        dg_frag%current_dia(:) = 0.0d0
        dg_frag%current_total(:) = 0.0d0
        dg_frag%dipole_lg_raw(:) = 0.0d0
        dg_frag%polarization_lg(:) = 0.0d0
        rt%curr(:, itt) = 0.0d0
        return
      end if
      stop "DG-Fragment RT: mixed/PW observable route was removed"
    end if
    if (.not. prop_obs_series_initialized) then
      prop_obs_series_source_value = 'prod_reference'
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PROP_OBS_SERIES_DIAG', prop_obs_series_env, &
        length=env_len, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(prop_obs_series_env(1:env_len))))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          prop_obs_series_enabled = .true.
        end select
      end if
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PROP_OBS_SERIES_SOURCE', &
        prop_obs_series_source_env, length=env_len, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(prop_obs_series_source_env(1:env_len))))
        case('prod_reference')
          prop_obs_series_source_value = 'prod_reference'
          prop_obs_series_enabled = .true.
        case('local_mixedz')
          prop_obs_series_source_value = 'local_mixedz'
          prop_obs_series_enabled = .true.
        case('local_mixedz_writeback')
          prop_obs_series_source_value = 'local_mixedz_writeback'
          prop_obs_series_enabled = .true.
        end select
      end if
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PROP_OBS_SERIES_LOCAL_MIXEDZ', &
        prop_obs_series_local_mixedz_env, length=env_len, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(prop_obs_series_local_mixedz_env(1:env_len))))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          prop_obs_series_source_value = 'local_mixedz'
          prop_obs_series_enabled = .true.
        end select
      end if
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PROP_OBS_SERIES_CURRENT_DETAIL', &
        prop_obs_series_current_detail_env, length=env_len, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(prop_obs_series_current_detail_env(1:env_len))))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enable_prop_obs_series_current_detail = .true.
        end select
      end if
      prop_obs_series_initialized = .true.
    end if
    prop_obs_candidate_source = prop_obs_series_source_value
    if (.not. pz_obs_hook_initialized) then
      ! Pz exposes only the accepted TOTAL dry-run/writeback path.
      enable_pz_total_writeback = (yn_dg_mixed_z_local_pz_writeback_total == 'y')
      if (enable_pz_total_writeback) enable_pz_total_hook = .true.
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PZ_REPLACE_TOTAL', pz_obs_hook_env, &
        length=env_len, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(pz_obs_hook_env(1:env_len))))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enable_pz_total_hook = .true.
        end select
      end if
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PZ_DRYRUN_TOTAL', pz_obs_hook_env, &
        length=env_len, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(pz_obs_hook_env(1:env_len))))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enable_pz_total_hook = .true.
        end select
      end if
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PZ_WRITEBACK_TOTAL', pz_obs_hook_env, &
        length=env_len, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(pz_obs_hook_env(1:env_len))))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enable_pz_total_hook = .true.
          enable_pz_total_writeback = .true.
        case('0','n','N','no','NO','false','FALSE','off','OFF')
          enable_pz_total_writeback = .false.
        end select
      end if
      pz_obs_hook_initialized = .true.
    end if
    if (.not. dg_frag%parallel_mode_orbital) then
      stop "DG-Fragment RT: non-orbital observable route was removed"
    end if

    if (enable_perf_count) dg_frag%mixed_z_perf_wall_obs_prepare = &
      dg_frag%mixed_z_perf_wall_obs_prepare + (get_wtime() - perf_obs_t0)
    perf_current_t0 = 0.0d0
    if (enable_perf_count) perf_current_t0 = get_wtime()
    perf_section_t0 = 0.0d0
    if (enable_perf_count) perf_section_t0 = get_wtime()
    call calculate_macroscopic_current_dg(dg_frag, system, mg, stencil, current_raw)
    if (enable_perf_count) dg_frag%mixed_z_perf_wall_current_para = &
      dg_frag%mixed_z_perf_wall_current_para + (get_wtime() - perf_section_t0)
    if (enable_perf_count) perf_section_t0 = get_wtime()
    call calculate_nonlocal_current_dg(dg_frag, system, mg, ppg, rt%Ac_tot(:, itt), current_nl_raw)
    if (enable_perf_count) dg_frag%mixed_z_perf_wall_current_nl = &
      dg_frag%mixed_z_perf_wall_current_nl + (get_wtime() - perf_section_t0)
    if (system%ngrid > 0) then
      current_para_candidate(:) = current_raw(:) / dble(system%ngrid)
      current_nl_candidate(:) = current_nl_raw(:) / dble(system%ngrid)
    else
      current_para_candidate(:) = 0.0d0
      current_nl_candidate(:) = 0.0d0
    end if
    current_para_total_candidate(:) = current_para_candidate(:)
    current_para_candidate_source = 'momentum_blocks_WW'
    if (dg_frag%current_momentum_decomp_ready) then
      current_para_total_candidate(:) = &
        (dg_frag%current_momentum_self_raw(:) + dg_frag%current_momentum_cross_raw(:)) / &
        dble(max(1, system%ngrid))
      current_para_candidate_source = 'fragment_momentum_self_cross'
    end if
    dg_frag%current(:) = current_para_candidate(:)
    dg_frag%current_nl(:) = current_nl_candidate(:)
    dg_frag%current_para(:) = dg_frag%current(:)
    if (enable_perf_count) perf_section_t0 = get_wtime()
    nelec_ref = 0.0d0
    if (allocated(system%rocc)) then
      do ispin = 1, min(dg_frag%nspin, size(system%rocc, 3))
        nelec_ref = nelec_ref + sum(max(0.0d0, system%rocc(1:min(dg_frag%nstate_tot, size(system%rocc, 1)), 1, ispin)))
      end do
    else if (dg_frag%nspin == 1) then
      nelec_ref = 2.0d0 * dble(max(0, dg_frag%nocc_spin(1)))
    else
      nelec_ref = dble(sum(max(0, dg_frag%nocc_spin(1:dg_frag%nspin))))
    end if
    if (system%ngrid > 0) then
      ne_density = nelec_ref / dble(system%ngrid)
    else
      ne_density = 0.0d0
    end if
    current_dia_candidate(:) = ne_density * rt%Ac_tot(:, itt)
    dg_frag%current_dia(:) = current_dia_candidate(:)
    if (enable_perf_count) dg_frag%mixed_z_perf_wall_current_dia = &
      dg_frag%mixed_z_perf_wall_current_dia + (get_wtime() - perf_section_t0)
    dg_frag%current_total(:) = current_para_candidate(:) + current_nl_candidate(:) + current_dia_candidate(:)
    current_candidate_total(:) = current_para_total_candidate(:) + current_nl_candidate(:) + current_dia_candidate(:)
    prop_obs_current_available = all(current_candidate_total(:) == current_candidate_total(:))
    prop_obs_current_prod_norm = sqrt(sum(dg_frag%current_total(:)**2))
    prop_obs_current_candidate_norm = sqrt(sum(current_candidate_total(:)**2))
    prop_obs_current_diff_norm = sqrt(sum((current_candidate_total(:) - dg_frag%current_total(:))**2))
    prop_obs_current_para_prod_norm = sqrt(sum(current_para_candidate(:)**2))
    prop_obs_current_nl_prod_norm = sqrt(sum(current_nl_candidate(:)**2))
    prop_obs_current_dia_prod_norm = sqrt(sum(current_dia_candidate(:)**2))
    if (.not. current_total_hook_initialized) then
      ! Current replacement is TOTAL-only: para_WW + nonlocal + diamagnetic.
      ! Para-only or WW-only current replacement paths are intentionally absent.
      enable_current_total_writeback = (yn_dg_mixed_z_local_current_writeback_total == 'y')
      if (enable_current_total_writeback) enable_current_total_hook = .true.
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_CURRENT_DRYRUN_TOTAL', current_total_hook_env, &
        length=env_len, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(current_total_hook_env(1:env_len))))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enable_current_total_hook = .true.
        end select
      end if
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_CURRENT_WRITEBACK_TOTAL', current_total_writeback_env, &
        length=env_len, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(current_total_writeback_env(1:env_len))))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enable_current_total_hook = .true.
          enable_current_total_writeback = .true.
        case('0','n','N','no','NO','false','FALSE','off','OFF')
          enable_current_total_writeback = .false.
        end select
      end if
	      current_total_hook_initialized = .true.
	    end if
    if (enable_current_total_hook .and. dg_frag%id == 0) then
      dg_frag%mixed_z_perf_current_writeback_calls = dg_frag%mixed_z_perf_current_writeback_calls + 1_8
      current_total_diff_ww(:) = dg_frag%current_total(:) - dg_frag%current_para(:)
      current_total_diff_full(:) = current_total_diff_ww(:)
      current_candidate_diff(:) = dg_frag%current_total(:) - current_candidate_total(:)
      current_before_writeback(:) = dg_frag%current_total(:)
      current_total_wp_norm = 0.0d0
      current_total_pp_norm = 0.0d0
      current_total_full_norm = sqrt(sum(current_candidate_total(:)**2))
      current_total_residual_norm = sqrt(sum(current_candidate_diff(:)**2))
      if (sqrt(sum(dg_frag%current_total(:)**2)) > 0.0d0) then
        current_total_p_ratio = sqrt(sum(current_total_diff_ww(:)**2)) / sqrt(sum(dg_frag%current_total(:)**2))
      else
        current_total_p_ratio = 0.0d0
      end if
      current_total_tol = 1.0d-12 * max(1.0d0, sqrt(sum(dg_frag%current_total(:)**2)), &
        sqrt(sum(dg_frag%current_para(:)**2)))
      current_total_block_reason = 'none'
      if (.not. all(current_candidate_total(:) == current_candidate_total(:))) then
        current_total_block_reason = 'nonfinite_candidate'
      else if (sqrt(sum(current_candidate_diff(:)**2)) > current_total_tol) then
        current_total_block_reason = 'diff_above_tol'
      else if (.not. enable_current_total_writeback) then
        current_total_block_reason = 'disabled_by_flag'
      end if
      if (enable_current_total_writeback .and. sqrt(sum(current_candidate_diff(:)**2)) <= current_total_tol .and. &
          all(current_candidate_total(:) == current_candidate_total(:))) then
        dg_frag%current_total(:) = current_candidate_total(:)
      end if
      current_after_minus_candidate(:) = dg_frag%current_total(:) - current_candidate_total(:)
      current_after_minus_before(:) = dg_frag%current_total(:) - current_before_writeback(:)
      if (.not. enable_perf_count) then
        write(*,'(1x,a,1(a,i0),29(a,1pe12.4),10(a,l1),5(a,a))') &
          '[DG-MIXEDZ-LOCAL-CURRENT-TOTAL-PATH-CMP]', &
          ' step=', itt, &
          ' current_prod_reference_norm=', sqrt(sum(dg_frag%current_total(:)**2)), &
          ' current_para_WW_norm=', sqrt(sum(current_para_candidate(:)**2)), &
          ' current_nl_candidate_norm=', sqrt(sum(current_nl_candidate(:)**2)), &
          ' current_dia_candidate_norm=', sqrt(sum(current_dia_candidate(:)**2)), &
          ' current_candidate_total_norm=', sqrt(sum(current_candidate_total(:)**2)), &
          ' current_WW_contract_norm=', sqrt(sum(current_para_candidate(:)**2)), &
          ' current_WP_cross_contract_norm=', current_total_wp_norm, &
          ' current_PP_contract_norm=', current_total_pp_norm, &
          ' current_full_contract_norm=', current_total_full_norm, &
          ' current_nonlocal_norm=', sqrt(sum(current_nl_candidate(:)**2)), &
          ' current_diamagnetic_norm=', sqrt(sum(current_dia_candidate(:)**2)), &
          ' diff_prod_minus_WW_norm=', sqrt(sum(current_total_diff_ww(:)**2)), &
          ' diff_prod_minus_full_norm=', current_total_residual_norm, &
          ' diff_prod_minus_candidate_norm=', sqrt(sum(current_candidate_diff(:)**2)), &
          ' p_contribution_ratio=', current_total_p_ratio, &
          ' nelec_ref=', nelec_ref, &
          ' ne_density=', ne_density, &
          ' current_total_x=', dg_frag%current_total(1), &
          ' current_total_y=', dg_frag%current_total(2), &
          ' current_total_z=', dg_frag%current_total(3), &
          ' current_para_x=', current_para_candidate(1), &
          ' current_para_y=', current_para_candidate(2), &
          ' current_para_z=', current_para_candidate(3), &
          ' current_nl_x=', current_nl_candidate(1), &
          ' current_nl_y=', current_nl_candidate(2), &
          ' current_nl_z=', current_nl_candidate(3), &
          ' current_dia_x=', current_dia_candidate(1), &
          ' current_dia_y=', current_dia_candidate(2), &
          ' current_dia_z=', current_dia_candidate(3), &
          ' current_available=', .true., &
          ' nonlocal_included=', sqrt(sum(current_nl_candidate(:)**2)) > current_total_tol, &
          ' diamagnetic_included=', sqrt(sum(current_dia_candidate(:)**2)) > current_total_tol, &
          ' para_candidate_available=', .true., &
          ' nl_candidate_available=', .true., &
          ' dia_candidate_available=', .true., &
          ' current_replaced=', .false., &
          ' density_replaced=', .false., &
          ' production_replacement=', .false., &
          ' bad=', any(dg_frag%current_total(:) /= dg_frag%current_total(:)) .or. &
            any(dg_frag%current_para(:) /= dg_frag%current_para(:)) .or. &
            sqrt(sum(current_candidate_diff(:)**2)) > current_total_tol, &
          ' prod_ref_kind=', 'total_para_nl_dia', &
          ' para_source=', trim(current_para_candidate_source), &
          ' nl_source=', 'production_nonlocal_path', &
          ' dia_density_source=', 'WW_density_target', &
          ' reference_scope=', 'production_total_current'
        write(*,'(1x,a,1(a,i0),6(a,1pe12.4),7(a,l1),2(a,a))') &
          '[DG-MIXEDZ-LOCAL-CURRENT-WRITEBACK-TOTAL-CMP]', &
          ' step=', itt, &
          ' current_before_norm=', sqrt(sum(current_before_writeback(:)**2)), &
          ' current_candidate_total_norm=', sqrt(sum(current_candidate_total(:)**2)), &
          ' current_after_norm=', sqrt(sum(dg_frag%current_total(:)**2)), &
          ' diff_candidate_minus_before_norm=', sqrt(sum(current_candidate_diff(:)**2)), &
          ' after_minus_candidate_norm=', sqrt(sum(current_after_minus_candidate(:)**2)), &
          ' after_minus_before_norm=', sqrt(sum(current_after_minus_before(:)**2)), &
          ' candidate_available=', .true., &
          ' replacement_ready=', sqrt(sum(current_candidate_diff(:)**2)) <= current_total_tol, &
          ' replacement_applied=', enable_current_total_writeback .and. &
            sqrt(sum(current_candidate_diff(:)**2)) <= current_total_tol .and. &
            all(current_candidate_total(:) == current_candidate_total(:)), &
          ' production_value_modified=', enable_current_total_writeback .and. &
            sqrt(sum(current_candidate_diff(:)**2)) <= current_total_tol .and. &
            all(current_candidate_total(:) == current_candidate_total(:)), &
          ' current_replaced=', enable_current_total_writeback .and. &
            sqrt(sum(current_candidate_diff(:)**2)) <= current_total_tol .and. &
            all(current_candidate_total(:) == current_candidate_total(:)), &
          ' density_replaced=', .false., &
          ' bad=', any(dg_frag%current_total(:) /= dg_frag%current_total(:)) .or. &
            sqrt(sum(current_after_minus_candidate(:)**2)) > current_total_tol, &
          ' target_kind=', 'total_para_nl_dia', &
          ' replacement_block_reason=', trim(current_total_block_reason)
      end if
    end if
    if (enable_perf_count) dg_frag%mixed_z_perf_wall_current_writeback = &
      dg_frag%mixed_z_perf_wall_current_writeback + (get_wtime() - perf_current_t0)
    dg_frag%dipole_lg_raw(:) = 0.0d0
    dg_frag%polarization_lg(:) = 0.0d0
    if (yn_dg_length_gauge == 'y') then
      if (enable_perf_count) perf_pz_t0 = get_wtime()
      if (.not. pol_trace_initialized) then
        call get_environment_variable('SALMON_DG_POL_TRACE', pol_trace_env, length=env_len, status=env_status)
        if (env_status == 0) then
          select case(trim(adjustl(pol_trace_env(1:env_len))))
          case('1','y','Y','yes','YES','true','TRUE','on','ON')
            enable_pol_trace = .true.
          end select
        end if
        pol_trace_initialized = .true.
      end if
      if (enable_perf_count) perf_pz_section_t0 = get_wtime()
      call calculate_local_wannier_polarization_dg(dg_frag, system, polarization_raw)
      if (enable_perf_count) dg_frag%mixed_z_perf_wall_pz_prod = &
        dg_frag%mixed_z_perf_wall_pz_prod + (get_wtime() - perf_pz_section_t0)
      dg_frag%dipole_lg_raw(:) = polarization_raw(:)
      polarization_density(:) = 0.0d0
      if (system%ngrid > 0 .and. system%hvol > 0.0d0) &
        polarization_density(:) = polarization_raw(:) / (dble(system%ngrid) * system%hvol)
      if (.not. dg_frag%polarization_lg_ref_ready) then
        dg_frag%polarization_lg_ref(:) = polarization_density(:)
        dg_frag%polarization_lg_ref_ready = .true.
      end if
      dg_frag%polarization_lg(:) = polarization_density(:) - dg_frag%polarization_lg_ref(:)
      if (enable_perf_count) perf_pz_section_t0 = get_wtime()
      call apply_wpw_reduced_pz_to_production(dg_frag, system, itt)
      if (enable_perf_count) dg_frag%mixed_z_perf_wall_pz_reduced = &
        dg_frag%mixed_z_perf_wall_pz_reduced + (get_wtime() - perf_pz_section_t0)
      if ((enable_pz_total_hook .or. prop_obs_series_enabled) .and. &
          dg_frag%mixed_z_prod_pz_decomp_ready) then
        if (enable_perf_count) perf_pz_section_t0 = get_wtime()
        pz_prod_raw_before = dg_frag%dipole_lg_raw(3)
        pz_prod_output_before = dg_frag%polarization_lg(3)
        wwfull_pz_diag = dg_frag%mixed_z_prod_pz_ww_diag_raw
        wwfull_pz_offdiag = dg_frag%mixed_z_prod_pz_ww_offdiag_raw
        wwfull_pz_total = dg_frag%mixed_z_prod_pz_ww_raw
        wwfull_pz_diff = wwfull_pz_total - pz_prod_raw_before
        pz_total_candidate_raw = wwfull_pz_total
        wwfull_replacement_ready = abs(wwfull_pz_diff) <= &
          1.0d-10 * max(1.0d0, abs(pz_prod_raw_before), abs(wwfull_pz_total))
        if (system%ngrid > 0) then
          if (system%hvol > 0.0d0) pz_prod_density = &
            pz_total_candidate_raw / (dble(system%ngrid) * system%hvol)
        else
          pz_prod_density = 0.0d0
        end if
        if (dg_frag%id == 0) then
          if (enable_pz_total_hook .or. prop_obs_series_enabled) then
            pz_total_wp_combined = dg_frag%mixed_z_prod_pz_wp_raw
            pz_total_pp = dg_frag%mixed_z_prod_pz_pp_raw
            pz_total_candidate = wwfull_pz_total + pz_total_wp_combined + pz_total_pp
            pz_total_diff = pz_total_candidate - pz_prod_raw_before
            pz_total_tol = 1.0d-10 * max(1.0d0, abs(pz_prod_raw_before), abs(pz_total_candidate))
            if (system%ngrid > 0) then
              if (system%hvol > 0.0d0) pz_total_density = &
                pz_total_candidate / (dble(system%ngrid) * system%hvol)
            else
              pz_total_density = 0.0d0
            end if
            pz_total_output_after = pz_total_density - dg_frag%polarization_lg_ref(3)
            pz_total_block_reason = 'none'
            if (pz_total_candidate /= pz_total_candidate .or. pz_prod_raw_before /= pz_prod_raw_before) then
              pz_total_block_reason = 'nonfinite_candidate'
            else if (abs(pz_total_diff) > pz_total_tol) then
              pz_total_block_reason = 'diff_above_tol'
            else if (.not. enable_pz_total_writeback) then
              pz_total_block_reason = 'disabled_by_flag'
            end if
            if (enable_pz_total_writeback .and. abs(pz_total_diff) <= pz_total_tol .and. &
                pz_total_candidate == pz_total_candidate .and. pz_prod_raw_before == pz_prod_raw_before) then
              dg_frag%dipole_lg_raw(3) = pz_total_candidate
              dg_frag%polarization_lg(3) = pz_total_output_after
            end if
            prop_obs_pz_available = .true.
            if (enable_pz_total_hook) then
              dg_frag%mixed_z_perf_pz_writeback_calls = dg_frag%mixed_z_perf_pz_writeback_calls + 1_8
              if (.not. enable_perf_count) then
                write(*,'(1x,a,1(a,i0),18(a,1pe23.15),6(a,l1),2(a,a))') &
                  '[DG-MIXEDZ-LOCAL-PZ-TOTAL-PATH-CMP]', &
                  ' step=', itt, &
                  ' Pz_before=', pz_prod_raw_before, &
                  ' Pz_candidate_total=', pz_total_candidate, &
                  ' Pz_after=', dg_frag%dipole_lg_raw(3), &
                  ' diff_candidate_minus_before=', pz_total_diff, &
                  ' WW_full=', wwfull_pz_total, &
                  ' WW_diag=', wwfull_pz_diag, &
                  ' WW_offdiag=', wwfull_pz_offdiag, &
                  ' WW_same_owner=', dg_frag%mixed_z_prod_pz_ww_same_owner_raw, &
                  ' WW_cross_owner=', dg_frag%mixed_z_prod_pz_ww_cross_owner_raw, &
                  ' WP_combined=', pz_total_wp_combined, &
                  ' PP=', pz_total_pp, &
                  ' replacement_tol=', pz_total_tol, &
                  ' zww_diag_mean=', dg_frag%mixed_z_prod_zww_diag_mean, &
                  ' center_z_mean=', dg_frag%mixed_z_prod_center_z_mean, &
                  ' diag_minus_center_rms=', dg_frag%mixed_z_prod_diag_minus_center_rms, &
                  ' weighted_zww=', dg_frag%mixed_z_prod_weighted_zww_diag_sum, &
                  ' weighted_center=', dg_frag%mixed_z_prod_weighted_center_sum, &
                  ' weighted_diff=', dg_frag%mixed_z_prod_weighted_diff_sum, &
                  ' replacement_ready=', abs(pz_total_diff) <= pz_total_tol, &
                  ' replacement_applied=', enable_pz_total_writeback .and. abs(pz_total_diff) <= pz_total_tol .and. &
                    pz_total_candidate == pz_total_candidate .and. pz_prod_raw_before == pz_prod_raw_before, &
                  ' production_value_modified=', enable_pz_total_writeback .and. abs(pz_total_diff) <= pz_total_tol .and. &
                    pz_total_candidate == pz_total_candidate .and. pz_prod_raw_before == pz_prod_raw_before, &
                  ' decomp_ready=', dg_frag%mixed_z_prod_pz_decomp_ready, &
                  ' field_on=', (abs(rt%E_tot(3, itt)) + abs(rt%Ac_tot(3, itt))) > 1.0d-14, &
                  ' bad=', pz_total_candidate /= pz_total_candidate .or. pz_prod_raw_before /= pz_prod_raw_before, &
                  ' target_kind=', 'total', &
                  ' replacement_block_reason=', trim(pz_total_block_reason)
                write(*,'(1x,a,1(a,i0),6(a,1pe12.4),7(a,l1),2(a,a))') &
                  '[DG-MIXEDZ-LOCAL-PZ-WRITEBACK-TOTAL-CMP]', &
                  ' step=', itt, &
                  ' Pz_before=', pz_prod_raw_before, &
                  ' Pz_candidate_total=', pz_total_candidate, &
                  ' Pz_after=', dg_frag%dipole_lg_raw(3), &
                  ' diff_candidate_minus_before=', pz_total_diff, &
                  ' after_minus_candidate=', dg_frag%dipole_lg_raw(3) - pz_total_candidate, &
                  ' after_minus_before=', dg_frag%dipole_lg_raw(3) - pz_prod_raw_before, &
                  ' candidate_available=', dg_frag%mixed_z_prod_pz_decomp_ready, &
                  ' replacement_ready=', abs(pz_total_diff) <= pz_total_tol .and. &
                    pz_total_candidate == pz_total_candidate .and. pz_prod_raw_before == pz_prod_raw_before, &
                  ' replacement_applied=', enable_pz_total_writeback .and. abs(pz_total_diff) <= pz_total_tol .and. &
                    pz_total_candidate == pz_total_candidate .and. pz_prod_raw_before == pz_prod_raw_before, &
                  ' production_value_modified=', enable_pz_total_writeback .and. abs(pz_total_diff) <= pz_total_tol .and. &
                    pz_total_candidate == pz_total_candidate .and. pz_prod_raw_before == pz_prod_raw_before, &
                  ' current_replaced=', .false., &
                  ' density_replaced=', .false., &
                  ' bad=', pz_total_candidate /= pz_total_candidate .or. pz_prod_raw_before /= pz_prod_raw_before, &
                  ' target_kind=', 'total', &
                  ' replacement_block_reason=', trim(pz_total_block_reason)
              end if
            end if
          end if
        end if
        if (enable_perf_count) dg_frag%mixed_z_perf_wall_pz_candidate = &
          dg_frag%mixed_z_perf_wall_pz_candidate + (get_wtime() - perf_pz_section_t0)
      end if
      if (enable_perf_count) dg_frag%mixed_z_perf_wall_pz_writeback = &
        dg_frag%mixed_z_perf_wall_pz_writeback + (get_wtime() - perf_pz_t0)
      if (prop_obs_series_enabled .and. &
          (trim(prop_obs_candidate_source) == 'local_mixedz' .or. &
           trim(prop_obs_candidate_source) == 'local_mixedz_writeback') .and. &
          dg_frag%mixed_z_local_prop_payload_ready .and. &
          dg_frag%mixed_z_local_prop_payload_step == itt .and. &
          .not. dg_frag%mixed_z_local_prop_payload_bad) then
        if (enable_prop_obs_series_current_detail .or. .not. dg_frag%current_momentum_decomp_ready) then
          if (enable_perf_count) then
            perf_payload_current_t0 = get_wtime()
            dg_frag%mixed_z_perf_payload_current_calls = dg_frag%mixed_z_perf_payload_current_calls + 1_8
          end if
          call compute_payload_ww_para_current(prop_obs_payload_current, prop_obs_payload_mom_current, &
            prop_obs_payload_mom_self_current, prop_obs_payload_mom_cross_current, &
            prop_obs_current_mom_block_counts, prop_obs_payload_current_ready)
          if (enable_perf_count) dg_frag%mixed_z_perf_wall_payload_current = &
            dg_frag%mixed_z_perf_wall_payload_current + (get_wtime() - perf_payload_current_t0)
        else
          prop_obs_payload_current_ready = .true.
        end if
      end if
      if (prop_obs_series_enabled .and. dg_frag%id == 0 .and. &
          dg_frag%mixed_z_local_prop_rho_step == itt) then
        perf_series_t0 = 0.0d0
        if (enable_perf_count) then
          perf_series_t0 = get_wtime()
          dg_frag%mixed_z_perf_series_validation_calls = &
            dg_frag%mixed_z_perf_series_validation_calls + 1_8
        end if
        prop_obs_field_on = (abs(rt%E_tot(3, itt)) + abs(rt%Ac_tot(3, itt))) > 1.0d-14
        if (prop_obs_field_on) then
          prop_obs_field_kind = 'field_on'
        else
          prop_obs_field_kind = 'field_off'
        end if
        prop_obs_payload_source = 'prod_reference_coef'
        prop_obs_payload_valid = .true.
        prop_obs_candidate_contraction_source = 'prod_observable_source'
        prop_obs_candidate_coef_norm = dg_frag%mixed_z_local_prop_payload_prod_coef_norm
        prop_obs_candidate_occ_trace = dg_frag%mixed_z_local_prop_rho_candidate_int
        prop_obs_rho_prod_int = dg_frag%mixed_z_local_prop_rho_prod_int
        prop_obs_rho_candidate_int = dg_frag%mixed_z_local_prop_rho_candidate_int
        prop_obs_rho_diff_int = dg_frag%mixed_z_local_prop_rho_diff_int
        prop_obs_rho_max_abs_diff = dg_frag%mixed_z_local_prop_rho_max_abs_diff
        prop_obs_rho_rms_diff = dg_frag%mixed_z_local_prop_rho_rms_diff
        prop_obs_pz_prod = pz_prod_raw_before
        prop_obs_pz_candidate = pz_total_candidate
        prop_obs_pz_diff = prop_obs_pz_candidate - prop_obs_pz_prod
        prop_obs_pz_prod_prev = prop_obs_prev_pz_value
        prop_obs_prev_pz_available = prop_obs_prev_pz_valid .and. prop_obs_prev_pz_step == itt - 1
        prop_obs_pz_diff_to_current = prop_obs_pz_candidate - prop_obs_pz_prod
        prop_obs_pz_diff_to_prev = 0.0d0
        if (prop_obs_prev_pz_available) prop_obs_pz_diff_to_prev = prop_obs_pz_candidate - prop_obs_pz_prod_prev
        prop_obs_payload_coef_time_tag = 'prod_observable_source'
        if (trim(prop_obs_candidate_source) == 'local_mixedz' .or. &
            trim(prop_obs_candidate_source) == 'local_mixedz_writeback') then
          if (dg_frag%mixed_z_local_prop_payload_ready .and. &
              dg_frag%mixed_z_local_prop_payload_step == itt) then
            prop_obs_payload_source = 'local_prop_coef'
            prop_obs_payload_valid = .not. dg_frag%mixed_z_local_prop_payload_bad
            if (prop_obs_payload_valid) then
              ! local_mixedz series uses the propagated payload for source identity
              ! and coefficient diagnostics, but the observable values stay on the
              ! accepted production-compatible targets: rho/Pz from the established
              ! W/TOTAL paths, current from fragment-basis momentum self+cross.
              prop_obs_candidate_contraction_source = 'local_prop_coef'
              if (trim(dg_frag%mixed_z_local_prop_payload_build_route) == &
                  'normal_series_after_propagation') then
                prop_obs_payload_coef_time_tag = 'after_propagation'
              else
                prop_obs_payload_coef_time_tag = 'before_observable'
              end if
              prop_obs_candidate_coef_norm = dg_frag%mixed_z_local_prop_payload_coef_norm
              prop_obs_candidate_occ_trace = dg_frag%mixed_z_local_prop_payload_occ_trace
              prop_obs_payload_pz = 0.0d0
              if (allocated(dg_frag%mixed_z_local_prop_payload_wcoef) .and. &
                  allocated(dg_frag%mixed_wannier_bpw_z)) then
                prop_obs_nw = min(size(dg_frag%mixed_z_local_prop_payload_wcoef, 1), &
                  size(dg_frag%mixed_wannier_bpw_z, 2), size(dg_frag%mixed_wannier_bpw_z, 3))
                prop_obs_nstate = min(size(dg_frag%mixed_z_local_prop_payload_wcoef, 2), &
                  size(system%rocc, 1))
                prop_obs_nspin = min(size(dg_frag%mixed_z_local_prop_payload_wcoef, 3), &
                  size(system%rocc, 3), size(dg_frag%mixed_wannier_bpw_z, 4))
                do prop_obs_ispin = 1, prop_obs_nspin
                  do prop_obs_ist = 1, prop_obs_nstate
                    prop_obs_occ = max(0.0d0, system%rocc(prop_obs_ist, 1, prop_obs_ispin))
                    if (prop_obs_occ <= 0.0d0) cycle
                    do prop_obs_iw = 1, prop_obs_nw
                      do prop_obs_jw = 1, prop_obs_nw
                        prop_obs_zterm = &
                          conjg(dg_frag%mixed_z_local_prop_payload_wcoef(prop_obs_iw,prop_obs_ist,prop_obs_ispin)) * &
                          dg_frag%mixed_wannier_bpw_z(3,prop_obs_iw,prop_obs_jw,prop_obs_ispin) * &
                          dg_frag%mixed_z_local_prop_payload_wcoef(prop_obs_jw,prop_obs_ist,prop_obs_ispin)
                        prop_obs_payload_pz = prop_obs_payload_pz - prop_obs_occ * real(prop_obs_zterm, kind=8)
                      end do
                    end do
                  end do
                end do
              end if
              prop_obs_pz_candidate = pz_total_candidate
              prop_obs_pz_diff = prop_obs_pz_candidate - prop_obs_pz_prod
              prop_obs_pz_diff_to_current = prop_obs_pz_candidate - prop_obs_pz_prod
              prop_obs_pz_diff_to_prev = 0.0d0
              if (prop_obs_prev_pz_available) prop_obs_pz_diff_to_prev = &
                prop_obs_pz_candidate - prop_obs_pz_prod_prev
              prop_obs_current_available = prop_obs_payload_current_ready
              if (prop_obs_current_available) then
                if (dg_frag%current_momentum_decomp_ready) then
                  prop_obs_current_fragbasis_vec(:) = &
                    (dg_frag%current_momentum_self_raw(:) + dg_frag%current_momentum_cross_raw(:)) / &
                    dble(max(1, system%ngrid))
                  prop_obs_current_candidate_vec(:) = prop_obs_current_fragbasis_vec(:) + &
                    current_nl_candidate(:) + current_dia_candidate(:)
                else
                  prop_obs_current_candidate_vec(:) = prop_obs_payload_mom_current(:) / dble(max(1, system%ngrid)) + &
                    current_nl_candidate(:) + current_dia_candidate(:)
                end if
                prop_obs_current_candidate_norm = sqrt(sum(prop_obs_current_candidate_vec(:)**2))
                prop_obs_current_diff_norm = sqrt(sum((prop_obs_current_candidate_vec(:) - &
                  dg_frag%current_total(:))**2))
              else
                prop_obs_current_candidate_norm = 0.0d0
                prop_obs_current_diff_norm = prop_obs_current_prod_norm
              end if
            end if
          else
            prop_obs_payload_source = 'missing_payload'
            prop_obs_payload_valid = .false.
          end if
        end if
        prop_obs_bad = dg_frag%mixed_z_local_prop_rho_bad .or. &
          .not. dg_frag%mixed_z_local_prop_rho_ready .or. &
          dg_frag%mixed_z_local_prop_rho_step /= itt .or. &
          .not. prop_obs_payload_valid .or. &
          .not. prop_obs_pz_available .or. .not. prop_obs_current_available .or. &
          prop_obs_pz_candidate /= prop_obs_pz_candidate .or. &
          pz_prod_raw_before /= pz_prod_raw_before .or. &
          prop_obs_current_prod_norm /= prop_obs_current_prod_norm .or. &
          prop_obs_current_candidate_norm /= prop_obs_current_candidate_norm .or. &
          prop_obs_current_diff_norm /= prop_obs_current_diff_norm
        if (.not. enable_perf_count) then
          write(*,'(1x,a,2(a,i0),25(a,1pe12.4),8(a,l1),8(a,a))') &
            '[DG-MIXEDZ-LOCAL-PROP-OBS-SERIES-CMP]', &
            ' step=', itt, &
            ' coef_source_step=', dg_frag%mixed_z_local_prop_payload_step, &
            ' rho_prod=', prop_obs_rho_prod_int, &
            ' rho_candidate=', prop_obs_rho_candidate_int, &
            ' rho_diff=', prop_obs_rho_diff_int, &
            ' rho_max_abs_diff=', prop_obs_rho_max_abs_diff, &
            ' rho_rms_diff=', prop_obs_rho_rms_diff, &
            ' Pz_prod=', prop_obs_pz_prod, &
            ' Pz_candidate=', prop_obs_pz_candidate, &
            ' Pz_diff=', prop_obs_pz_diff, &
            ' Pz_prod_prev=', prop_obs_pz_prod_prev, &
            ' diff_to_current=', prop_obs_pz_diff_to_current, &
            ' diff_to_prev=', prop_obs_pz_diff_to_prev, &
            ' current_prod_norm=', prop_obs_current_prod_norm, &
            ' current_candidate_norm=', prop_obs_current_candidate_norm, &
            ' current_diff_norm=', prop_obs_current_diff_norm, &
            ' current_para_prod_norm=', prop_obs_current_para_prod_norm, &
            ' current_nl_prod_norm=', prop_obs_current_nl_prod_norm, &
            ' current_dia_prod_norm=', prop_obs_current_dia_prod_norm, &
            ' candidate_coef_norm=', prop_obs_candidate_coef_norm, &
            ' candidate_occ_trace=', prop_obs_candidate_occ_trace, &
            ' coef_norm=', dg_frag%mixed_z_local_prop_payload_coef_norm, &
            ' prod_coef_norm=', dg_frag%mixed_z_local_prop_payload_prod_coef_norm, &
            ' coef_diff_Snorm=', dg_frag%mixed_z_local_prop_payload_coef_diff_snorm, &
            ' rel_coef_diff_Snorm=', dg_frag%mixed_z_local_prop_payload_rel_coef_diff_snorm, &
            ' payload_dim=', dg_frag%mixed_z_local_prop_payload_dim, &
            ' payload_occ_trace=', dg_frag%mixed_z_local_prop_payload_occ_trace, &
            ' rho_candidate_available=', dg_frag%mixed_z_local_prop_rho_ready, &
            ' Pz_candidate_available=', prop_obs_pz_available, &
            ' current_candidate_available=', prop_obs_current_available, &
            ' payload_valid=', prop_obs_payload_valid, &
            ' payload_step_match=', dg_frag%mixed_z_local_prop_payload_step == itt, &
            ' prev_Pz_available=', prop_obs_prev_pz_available, &
            ' field_on=', prop_obs_field_on, &
            ' bad=', prop_obs_bad, &
            ' field_kind=', trim(prop_obs_field_kind), &
            ' candidate_source=', trim(prop_obs_candidate_source), &
            ' candidate_contraction_source=', trim(prop_obs_candidate_contraction_source), &
            ' payload_coef_time_tag=', trim(prop_obs_payload_coef_time_tag), &
            ' payload_source=', trim(prop_obs_payload_source), &
            ' payload_basis_kind=', trim(dg_frag%mixed_z_local_prop_payload_basis_kind), &
            ' payload_build_route=', trim(dg_frag%mixed_z_local_prop_payload_build_route), &
            ' payload_block_reason=', trim(dg_frag%mixed_z_local_prop_payload_block_reason)
        end if
        if (enable_perf_count) dg_frag%mixed_z_perf_wall_series = &
          dg_frag%mixed_z_perf_wall_series + (get_wtime() - perf_series_t0)
        prop_obs_prev_pz_value = prop_obs_pz_prod
        prop_obs_prev_pz_step = itt
        prop_obs_prev_pz_valid = .true.
      end if
      if (enable_pol_trace .and. itt <= 5 .and. dg_frag%id == 0) then
        write(*,'(1x,a,i0,9(a,1pe13.5))') '[DG-POL-TRACE] itt=', itt, &
          ' raw_x=', polarization_raw(1), ' raw_y=', polarization_raw(2), ' raw_z=', polarization_raw(3), &
          ' den_x=', polarization_density(1), ' den_y=', polarization_density(2), ' den_z=', polarization_density(3), &
          ' dP_x=', dg_frag%polarization_lg(1), ' dP_y=', dg_frag%polarization_lg(2), &
          ' dP_z=', dg_frag%polarization_lg(3)
        flush(6)
      end if
    end if
    if (enable_perf_count) then
      if (dg_frag%mixed_z_local_prop_rho_ready .and. dg_frag%mixed_z_local_prop_rho_step == itt) then
        stability_rho = dg_frag%mixed_z_local_prop_rho_candidate_int
      else
        stability_rho = dg_frag%rho_global_scaled_elec
      end if
      stability_pz = dg_frag%polarization_lg(3)
      stability_current = sqrt(sum(dg_frag%current_total(:)**2))
      stability_norm = dg_frag%coef_occ_norm
      stability_energy = dg_frag%total_energy
      stability_bad = stability_rho /= stability_rho .or. stability_pz /= stability_pz .or. &
        stability_current /= stability_current .or. stability_norm /= stability_norm .or. &
        stability_energy /= stability_energy
      if (.not. dg_frag%mixed_z_frag_local_stability_baseline_ready .and. .not. stability_bad) then
        dg_frag%mixed_z_frag_local_rho_baseline = stability_rho
        dg_frag%mixed_z_frag_local_pz_baseline = stability_pz
        dg_frag%mixed_z_frag_local_current_baseline = stability_current
        dg_frag%mixed_z_frag_local_norm_baseline = stability_norm
        dg_frag%mixed_z_frag_local_energy_baseline = stability_energy
        dg_frag%mixed_z_frag_local_stability_baseline_ready = .true.
      end if
      if (dg_frag%mixed_z_frag_local_stability_baseline_ready .and. .not. stability_bad) then
        dg_frag%mixed_z_frag_local_rho_drift_max = max(dg_frag%mixed_z_frag_local_rho_drift_max, &
          abs(stability_rho - dg_frag%mixed_z_frag_local_rho_baseline))
        dg_frag%mixed_z_frag_local_pz_drift_max = max(dg_frag%mixed_z_frag_local_pz_drift_max, &
          abs(stability_pz - dg_frag%mixed_z_frag_local_pz_baseline))
        dg_frag%mixed_z_frag_local_current_drift_max = max(dg_frag%mixed_z_frag_local_current_drift_max, &
          abs(stability_current - dg_frag%mixed_z_frag_local_current_baseline))
        dg_frag%mixed_z_frag_local_norm_drift_max = max(dg_frag%mixed_z_frag_local_norm_drift_max, &
          abs(stability_norm - dg_frag%mixed_z_frag_local_norm_baseline), abs(dg_frag%coef_occ_norm_drift))
        dg_frag%mixed_z_frag_local_energy_drift_max = max(dg_frag%mixed_z_frag_local_energy_drift_max, &
          abs(stability_energy - dg_frag%mixed_z_frag_local_energy_baseline))
      end if
      dg_frag%mixed_z_perf_wall_obs = dg_frag%mixed_z_perf_wall_obs + (get_wtime() - perf_obs_t0)
      perf_final_step = .false.
      if (dg_frag%mixed_z_perf_final_itt >= 0) perf_final_step = itt == dg_frag%mixed_z_perf_final_itt
      if (perf_final_step .and. dg_frag%id == 0) then
        write(*,'(1x,a,22(a,i0),30(a,1pe12.4),a,l1)') &
          '[DG-MIXEDZ-PERF-COUNT]', &
          ' nstep=', dg_frag%mixed_z_perf_nstep, &
          ' prop_writeback_calls=', dg_frag%mixed_z_perf_prop_writeback_calls, &
          ' rho_writeback_calls=', dg_frag%mixed_z_perf_rho_writeback_calls, &
          ' pz_writeback_calls=', dg_frag%mixed_z_perf_pz_writeback_calls, &
          ' current_writeback_calls=', dg_frag%mixed_z_perf_current_writeback_calls, &
          ' series_validation_calls=', dg_frag%mixed_z_perf_series_validation_calls, &
          ' payload_current_calls=', dg_frag%mixed_z_perf_payload_current_calls, &
          ' zgemm_calls=', dg_frag%mixed_z_perf_zgemm_calls, &
          ' eigensolve_calls=', dg_frag%mixed_z_perf_eigensolve_calls, &
          ' expdiag_calls=', dg_frag%mixed_z_perf_expdiag_calls, &
          ' mpi_reduce_calls=', dg_frag%mixed_z_perf_mpi_reduce_calls, &
          ' obs_mpi_reduce_calls=', dg_frag%mixed_z_perf_obs_mpi_reduce_calls, &
          ' prop_pack_calls=', dg_frag%mixed_z_perf_prop_pack_calls, &
          ' prop_unpack_calls=', dg_frag%mixed_z_perf_prop_unpack_calls, &
          ' prop_pack_allocs=', dg_frag%mixed_z_perf_prop_pack_allocs, &
          ' prop_unpack_allocs=', dg_frag%mixed_z_perf_prop_unpack_allocs, &
          ' prop_pack_zero_bytes=', dg_frag%mixed_z_perf_prop_pack_zero_bytes, &
          ' prop_unpack_zero_bytes=', dg_frag%mixed_z_perf_prop_unpack_zero_bytes, &
          ' prop_pack_w_copy_bytes=', dg_frag%mixed_z_perf_prop_pack_w_copy_bytes, &
          ' prop_pack_p_copy_bytes=', dg_frag%mixed_z_perf_prop_pack_p_copy_bytes, &
          ' prop_unpack_w_copy_bytes=', dg_frag%mixed_z_perf_prop_unpack_w_copy_bytes, &
          ' prop_unpack_p_copy_bytes=', dg_frag%mixed_z_perf_prop_unpack_p_copy_bytes, &
          ' wall_prop=', dg_frag%mixed_z_perf_wall_prop, &
          ' wall_prop_pack=', dg_frag%mixed_z_perf_wall_prop_pack, &
          ' wall_prop_comm=', dg_frag%mixed_z_perf_wall_prop_comm, &
          ' wall_prop_field_exp=', dg_frag%mixed_z_perf_wall_prop_field_exp, &
          ' wall_prop_phase=', dg_frag%mixed_z_perf_wall_prop_phase, &
          ' wall_prop_unpack=', dg_frag%mixed_z_perf_wall_prop_unpack, &
          ' wall_obs=', dg_frag%mixed_z_perf_wall_obs, &
          ' wall_series=', dg_frag%mixed_z_perf_wall_series, &
          ' wall_obs_prepare=', dg_frag%mixed_z_perf_wall_obs_prepare, &
          ' wall_rho_writeback=', dg_frag%mixed_z_perf_wall_rho_writeback, &
          ' wall_pz_writeback=', dg_frag%mixed_z_perf_wall_pz_writeback, &
          ' wall_pz_prod=', dg_frag%mixed_z_perf_wall_pz_prod, &
          ' wall_pz_reduced=', dg_frag%mixed_z_perf_wall_pz_reduced, &
          ' wall_pz_candidate=', dg_frag%mixed_z_perf_wall_pz_candidate, &
          ' wall_pz_build_cw=', dg_frag%mixed_z_perf_wall_pz_build_cw, &
          ' wall_pz_contract_z=', dg_frag%mixed_z_perf_wall_pz_contract_z, &
          ' wall_pz_reduce=', dg_frag%mixed_z_perf_wall_pz_reduce, &
          ' wall_current_writeback=', dg_frag%mixed_z_perf_wall_current_writeback, &
          ' wall_payload_current=', dg_frag%mixed_z_perf_wall_payload_current, &
          ' wall_current_para=', dg_frag%mixed_z_perf_wall_current_para, &
          ' wall_current_nl=', dg_frag%mixed_z_perf_wall_current_nl, &
          ' wall_current_nl_cache=', dg_frag%mixed_z_perf_wall_current_nl_cache, &
          ' wall_current_nl_setup=', dg_frag%mixed_z_perf_wall_current_nl_setup, &
          ' wall_current_nl_fetch=', dg_frag%mixed_z_perf_wall_current_nl_fetch, &
          ' wall_current_nl_project=', dg_frag%mixed_z_perf_wall_current_nl_project, &
          ' wall_current_nl_contract=', dg_frag%mixed_z_perf_wall_current_nl_contract, &
          ' wall_current_nl_reduce=', dg_frag%mixed_z_perf_wall_current_nl_reduce, &
          ' wall_current_dia=', dg_frag%mixed_z_perf_wall_current_dia, &
          ' wall_obs_mpi_reduce=', dg_frag%mixed_z_perf_wall_obs_mpi_reduce, &
          ' wall_obs_unaccounted=', max(0.0d0, dg_frag%mixed_z_perf_wall_obs - &
            (dg_frag%mixed_z_perf_wall_obs_prepare + dg_frag%mixed_z_perf_wall_pz_writeback + &
             dg_frag%mixed_z_perf_wall_current_writeback + dg_frag%mixed_z_perf_wall_payload_current + &
             dg_frag%mixed_z_perf_wall_series)), &
          ' bad=', .false.
        write(*,'(1x,a,1(a,1pe12.4),1(a,a),1(a,i0),5(a,1pe12.4),1(a,l1))') &
          '[DG-MIXEDZ-FRAG-LOCAL-STABILITY-SUMMARY]', &
          ' field_abs=', dg_frag%mixed_z_frag_local_field_abs_max, &
          ' field_block_kind=', trim(dg_frag%mixed_z_frag_local_field_block_kind), &
          ' nt=', int(dg_frag%mixed_z_perf_prop_writeback_calls), &
          ' rho_drift_max=', dg_frag%mixed_z_frag_local_rho_drift_max, &
          ' Pz_drift_max=', dg_frag%mixed_z_frag_local_pz_drift_max, &
          ' current_drift_max=', dg_frag%mixed_z_frag_local_current_drift_max, &
          ' norm_drift_max=', dg_frag%mixed_z_frag_local_norm_drift_max, &
          ' energy_drift_max=', dg_frag%mixed_z_frag_local_energy_drift_max, &
          ' bad=', stability_bad .or. .not. dg_frag%mixed_z_frag_local_stability_baseline_ready
      end if
    end if
    rt%curr(:, itt) = -dg_frag%current_total(:)

  contains

    subroutine compute_payload_ww_para_current(current_ww, current_mom_ww, current_mom_self_ww, &
      current_mom_cross_ww, momentum_block_counts, available)
      real(8), intent(out) :: current_ww(3)
      real(8), intent(out) :: current_mom_ww(3)
      real(8), intent(out) :: current_mom_self_ww(3)
      real(8), intent(out) :: current_mom_cross_ww(3)
      real(8), intent(out) :: momentum_block_counts(3)
      logical, intent(out) :: available

      integer :: ispin_c, ifrag_c, ilocal_c, nspin_c, nstate_c, nw_payload
      integer :: n_w_local, nbf, ib, iw, jw, ist, idir
      integer :: iblk, io, jo, nrow_blk, ncol_blk, ifrag_col, ilocal_col
      integer :: n_w_col, nbf_col
      integer :: n_global_dense, ifrag_row_dense, ifrag_col_dense, gid_basis, gid_w
      integer :: ix, iy, iz, gx, gy, gz, bx, by, bz
      integer :: p_lb1, p_lb2, p_lb3, p_ub1, p_ub2, p_ub3
      integer :: gid_i, gid_j
      real(8) :: vol_weight, occ_c
      real(8) :: current_local(3), current_sum(3)
      real(8) :: current_mom_local(3), current_mom_sum(3)
      real(8) :: current_mom_self_local(3), current_mom_self_sum(3)
      real(8) :: current_mom_cross_local(3), current_mom_cross_sum(3)
      real(8) :: current_mom_self_dense_local(3), current_mom_cross_dense_local(3)
      real(8) :: momentum_block_counts_local(3), momentum_block_counts_sum(3)
      real(8) :: momentum_block_counts_dense_local(3)
      complex(8), parameter :: zi_local = (0.0d0, 1.0d0)
      complex(8) :: c_i, c_j
      complex(8), allocatable :: wval(:), wgrad(:,:), jww(:,:,:), jmom_block(:,:,:)
      complex(8), allocatable :: cbasis_local(:,:), cbasis_sum(:,:)

      current_ww(:) = 0.0d0
      current_mom_ww(:) = 0.0d0
      current_mom_self_ww(:) = 0.0d0
      current_mom_cross_ww(:) = 0.0d0
      momentum_block_counts(:) = 0.0d0
      current_local(:) = 0.0d0
      current_sum(:) = 0.0d0
      current_mom_local(:) = 0.0d0
      current_mom_sum(:) = 0.0d0
      current_mom_self_local(:) = 0.0d0
      current_mom_self_sum(:) = 0.0d0
      current_mom_cross_local(:) = 0.0d0
      current_mom_cross_sum(:) = 0.0d0
      current_mom_self_dense_local(:) = 0.0d0
      current_mom_cross_dense_local(:) = 0.0d0
      momentum_block_counts_local(:) = 0.0d0
      momentum_block_counts_sum(:) = 0.0d0
      momentum_block_counts_dense_local(:) = 0.0d0
      available = .false.

      if (.not. allocated(dg_frag%mixed_z_local_prop_payload_wcoef)) return
      if (.not. allocated(dg_frag%global_wannier_local_coef)) return
      if (.not. allocated(dg_frag%global_wannier_local_nkeep)) return
      if (.not. allocated(dg_frag%global_wannier_local_ids)) return
      if (.not. allocated(dg_frag%index_basis)) return
      if (.not. allocated(dg_frag%coef_owner)) return
      if (.not. allocated(system%rocc)) return
      if (.not. allocated(dg_frag%phi_frag) .and. .not. allocated(dg_frag%phi_frag_c)) return

      call ensure_gradient_basis_cache(dg_frag, mg, stencil)
      if (.not. allocated(dg_frag%gradient_basis_cache)) return

      nw_payload = size(dg_frag%mixed_z_local_prop_payload_wcoef, 1)
      nstate_c = min(size(dg_frag%mixed_z_local_prop_payload_wcoef, 2), size(system%rocc, 1))
      nspin_c = min(dg_frag%nspin, system%nspin, size(system%rocc, 3), &
        size(dg_frag%mixed_z_local_prop_payload_wcoef, 3))
      if (nw_payload <= 0 .or. nstate_c <= 0 .or. nspin_c <= 0) return

      vol_weight = product(dg_frag%hgs(1:3)) / &
        product(dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3)))

      n_global_dense = min(size(dg_frag%coef_owner, 1), dg_frag%n_mat_max)
      if (n_global_dense <= 0) return
      do ispin_c = 1, nspin_c
        allocate(cbasis_local(n_global_dense, nstate_c), cbasis_sum(n_global_dense, nstate_c))
        cbasis_local(:, :) = (0.0d0, 0.0d0)
        cbasis_sum(:, :) = (0.0d0, 0.0d0)
        do ifrag_c = dg_frag%ifrag_start, dg_frag%ifrag_end
          ilocal_c = ifrag_c - dg_frag%ifrag_start + 1
          if (ilocal_c < 1 .or. ilocal_c > size(dg_frag%global_wannier_local_nkeep)) cycle
          n_w_local = min(dg_frag%global_wannier_local_nkeep(ilocal_c), &
            size(dg_frag%global_wannier_local_coef, 2), size(dg_frag%global_wannier_local_ids, 1))
          nbf = min(dg_frag%n_basis(ifrag_c, ispin_c), size(dg_frag%global_wannier_local_coef, 1), &
            size(dg_frag%index_basis, 1))
          if (n_w_local <= 0 .or. nbf <= 0) cycle
          do iw = 1, n_w_local
            gid_w = dg_frag%global_wannier_local_ids(iw, ilocal_c)
            if (gid_w < 1 .or. gid_w > nw_payload) cycle
            do ib = 1, nbf
              gid_basis = dg_frag%index_basis(ib, ifrag_c, ispin_c)
              if (gid_basis < 1 .or. gid_basis > n_global_dense) cycle
              do ist = 1, nstate_c
                cbasis_local(gid_basis, ist) = cbasis_local(gid_basis, ist) + &
                  dg_frag%global_wannier_local_coef(ib, iw, ispin_c, ilocal_c) * &
                  dg_frag%mixed_z_local_prop_payload_wcoef(gid_w, ist, ispin_c)
              end do
            end do
          end do
        end do
        if (enable_perf_count) then
          dg_frag%mixed_z_perf_mpi_reduce_calls = dg_frag%mixed_z_perf_mpi_reduce_calls + 1_8
          dg_frag%mixed_z_perf_obs_mpi_reduce_calls = dg_frag%mixed_z_perf_obs_mpi_reduce_calls + 1_8
          perf_reduce_t0 = get_wtime()
        end if
        call comm_summation(cbasis_local, cbasis_sum, size(cbasis_local), dg_frag%icomm)
        if (enable_perf_count) dg_frag%mixed_z_perf_wall_obs_mpi_reduce = &
          dg_frag%mixed_z_perf_wall_obs_mpi_reduce + (get_wtime() - perf_reduce_t0)
        if (allocated(dg_frag%momentum_blocks)) then
          do iblk = 1, dg_frag%n_momentum_blocks
            if (.not. allocated(dg_frag%momentum_blocks(iblk)%val)) cycle
            ifrag_row_dense = dg_frag%momentum_blocks(iblk)%ifrag_row
            ifrag_col_dense = dg_frag%momentum_blocks(iblk)%ifrag_col
            if (ifrag_row_dense < 1 .or. ifrag_row_dense > dg_frag%n_frag) cycle
            if (ifrag_col_dense < 1 .or. ifrag_col_dense > dg_frag%n_frag) cycle
            nrow_blk = min(dg_frag%n_basis(ifrag_row_dense, ispin_c), &
              size(dg_frag%momentum_blocks(iblk)%val, 2), size(dg_frag%index_basis, 1))
            ncol_blk = min(dg_frag%n_basis(ifrag_col_dense, ispin_c), &
              size(dg_frag%momentum_blocks(iblk)%val, 3), size(dg_frag%index_basis, 1))
            if (nrow_blk <= 0 .or. ncol_blk <= 0) cycle
            do io = 1, nrow_blk
              gid_i = dg_frag%index_basis(io, ifrag_row_dense, ispin_c)
              if (gid_i < 1 .or. gid_i > n_global_dense) cycle
              if (dg_frag%coef_owner(gid_i, ispin_c) /= dg_frag%id) cycle
              do jo = 1, ncol_blk
                gid_j = dg_frag%index_basis(jo, ifrag_col_dense, ispin_c)
                if (gid_j < 1 .or. gid_j > n_global_dense) cycle
                do ist = 1, nstate_c
                  occ_c = max(0.0d0, system%rocc(ist, 1, ispin_c))
                  if (occ_c <= 0.0d0) cycle
                  do idir = 1, 3
                    if (ifrag_row_dense == ifrag_col_dense) then
                      current_mom_self_dense_local(idir) = current_mom_self_dense_local(idir) + &
                        occ_c * aimag(conjg(cbasis_sum(gid_i, ist)) * &
                        cmplx(dg_frag%momentum_blocks(iblk)%val(idir, io, jo, ispin_c), 0.0d0, kind=8) * &
                        cbasis_sum(gid_j, ist))
                    else
                      current_mom_cross_dense_local(idir) = current_mom_cross_dense_local(idir) + &
                        occ_c * aimag(conjg(cbasis_sum(gid_i, ist)) * &
                        cmplx(dg_frag%momentum_blocks(iblk)%val(idir, io, jo, ispin_c), 0.0d0, kind=8) * &
                        cbasis_sum(gid_j, ist))
                    end if
                  end do
                end do
              end do
            end do
            if (ifrag_row_dense == ifrag_col_dense) then
              momentum_block_counts_dense_local(1) = momentum_block_counts_dense_local(1) + 1.0d0
            else
              momentum_block_counts_dense_local(2) = momentum_block_counts_dense_local(2) + 1.0d0
            end if
          end do
        end if
        deallocate(cbasis_local, cbasis_sum)
      end do

      if (allocated(dg_frag%phi_frag_c)) then
        p_lb1 = lbound(dg_frag%phi_frag_c, 1); p_ub1 = ubound(dg_frag%phi_frag_c, 1)
        p_lb2 = lbound(dg_frag%phi_frag_c, 2); p_ub2 = ubound(dg_frag%phi_frag_c, 2)
        p_lb3 = lbound(dg_frag%phi_frag_c, 3); p_ub3 = ubound(dg_frag%phi_frag_c, 3)
      else
        p_lb1 = lbound(dg_frag%phi_frag, 1); p_ub1 = ubound(dg_frag%phi_frag, 1)
        p_lb2 = lbound(dg_frag%phi_frag, 2); p_ub2 = ubound(dg_frag%phi_frag, 2)
        p_lb3 = lbound(dg_frag%phi_frag, 3); p_ub3 = ubound(dg_frag%phi_frag, 3)
      end if

      do ispin_c = 1, nspin_c
        do ifrag_c = dg_frag%ifrag_start, dg_frag%ifrag_end
          ilocal_c = ifrag_c - dg_frag%ifrag_start + 1
          if (ilocal_c < 1 .or. ilocal_c > size(dg_frag%global_wannier_local_nkeep)) cycle
          n_w_local = min(dg_frag%global_wannier_local_nkeep(ilocal_c), &
            size(dg_frag%global_wannier_local_coef, 2), size(dg_frag%global_wannier_local_ids, 1))
          if (n_w_local <= 0) cycle
          nbf = min(dg_frag%n_basis(ifrag_c, ispin_c), size(dg_frag%global_wannier_local_coef, 1), &
            size(dg_frag%gradient_basis_cache, 5))
          if (allocated(dg_frag%phi_frag_c)) then
            nbf = min(nbf, size(dg_frag%phi_frag_c, 4))
          else
            nbf = min(nbf, size(dg_frag%phi_frag, 4))
          end if
          if (nbf <= 0) cycle

          allocate(wval(n_w_local), wgrad(3,n_w_local), jww(3,n_w_local,n_w_local))
          jww(:, :, :) = (0.0d0, 0.0d0)
          do iz = 1, dg_frag%nxyz_domain(3, ifrag_c)
            gz = dg_frag%ixyz_frag(3, ifrag_c) + iz - 1
            bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
            if (bz < p_lb3 .or. bz > p_ub3) cycle
            do iy = 1, dg_frag%nxyz_domain(2, ifrag_c)
              gy = dg_frag%ixyz_frag(2, ifrag_c) + iy - 1
              by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
              if (by < p_lb2 .or. by > p_ub2) cycle
              do ix = 1, dg_frag%nxyz_domain(1, ifrag_c)
                gx = dg_frag%ixyz_frag(1, ifrag_c) + ix - 1
                bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
                if (bx < p_lb1 .or. bx > p_ub1) cycle
                wval(:) = (0.0d0, 0.0d0)
                wgrad(:, :) = (0.0d0, 0.0d0)
                do iw = 1, n_w_local
                  do ib = 1, nbf
                    if (allocated(dg_frag%phi_frag_c)) then
                      wval(iw) = wval(iw) + &
                        dg_frag%global_wannier_local_coef(ib, iw, ispin_c, ilocal_c) * &
                        dg_frag%phi_frag_c(bx, by, bz, ib, ilocal_c)
                    else
                      wval(iw) = wval(iw) + &
                        dg_frag%global_wannier_local_coef(ib, iw, ispin_c, ilocal_c) * &
                        cmplx(dg_frag%phi_frag(bx, by, bz, ib, ilocal_c), 0.0d0, kind=8)
                    end if
                    do idir = 1, 3
                      wgrad(idir, iw) = wgrad(idir, iw) + &
                        dg_frag%global_wannier_local_coef(ib, iw, ispin_c, ilocal_c) * &
                        cmplx(dg_frag%gradient_basis_cache(ix, iy, iz, idir, ib, ilocal_c), 0.0d0, kind=8)
                    end do
                  end do
                end do
                do jw = 1, n_w_local
                  do iw = 1, n_w_local
                    do idir = 1, 3
                      jww(idir, iw, jw) = jww(idir, iw, jw) - &
                        zi_local * conjg(wval(iw)) * wgrad(idir, jw) * vol_weight
                    end do
                  end do
                end do
              end do
            end do
          end do

          do ist = 1, nstate_c
            occ_c = max(0.0d0, system%rocc(ist, 1, ispin_c))
            if (occ_c <= 0.0d0) cycle
            do iw = 1, n_w_local
              gid_i = dg_frag%global_wannier_local_ids(iw, ilocal_c)
              if (gid_i < 1 .or. gid_i > nw_payload) cycle
              c_i = dg_frag%mixed_z_local_prop_payload_wcoef(gid_i, ist, ispin_c)
              do jw = 1, n_w_local
                gid_j = dg_frag%global_wannier_local_ids(jw, ilocal_c)
                if (gid_j < 1 .or. gid_j > nw_payload) cycle
                c_j = dg_frag%mixed_z_local_prop_payload_wcoef(gid_j, ist, ispin_c)
                do idir = 1, 3
                  current_local(idir) = current_local(idir) + &
                    occ_c * real(conjg(c_i) * jww(idir, iw, jw) * c_j, kind=8)
                end do
              end do
            end do
          end do
          deallocate(wval, wgrad, jww)
        end do
      end do

      current_mom_self_local(:) = current_mom_self_dense_local(:)
      current_mom_cross_local(:) = current_mom_cross_dense_local(:)
      momentum_block_counts_local(:) = momentum_block_counts_dense_local(:)
      current_mom_local(:) = current_mom_self_local(:) + current_mom_cross_local(:)
      if (enable_perf_count) then
        dg_frag%mixed_z_perf_mpi_reduce_calls = dg_frag%mixed_z_perf_mpi_reduce_calls + 5_8
        dg_frag%mixed_z_perf_obs_mpi_reduce_calls = dg_frag%mixed_z_perf_obs_mpi_reduce_calls + 5_8
        perf_reduce_t0 = get_wtime()
      end if
      call comm_summation(current_local, current_sum, 3, dg_frag%icomm)
      call comm_summation(current_mom_local, current_mom_sum, 3, dg_frag%icomm)
      call comm_summation(current_mom_self_local, current_mom_self_sum, 3, dg_frag%icomm)
      call comm_summation(current_mom_cross_local, current_mom_cross_sum, 3, dg_frag%icomm)
      call comm_summation(momentum_block_counts_local, momentum_block_counts_sum, 3, dg_frag%icomm)
      if (enable_perf_count) dg_frag%mixed_z_perf_wall_obs_mpi_reduce = &
        dg_frag%mixed_z_perf_wall_obs_mpi_reduce + (get_wtime() - perf_reduce_t0)
      current_ww(:) = current_sum(:)
      current_mom_ww(:) = current_mom_sum(:)
      current_mom_self_ww(:) = current_mom_self_sum(:)
      current_mom_cross_ww(:) = current_mom_cross_sum(:)
      momentum_block_counts(:) = momentum_block_counts_sum(:)
      available = all(current_ww(:) == current_ww(:)) .and. all(current_mom_ww(:) == current_mom_ww(:)) .and. &
        all(current_mom_self_ww(:) == current_mom_self_ww(:)) .and. &
        all(current_mom_cross_ww(:) == current_mom_cross_ww(:))
    end subroutine compute_payload_ww_para_current
  end subroutine calculate_observables
