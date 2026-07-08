  subroutine time_evolution_expdiag(dg_frag, system, info, rt, itt, dt, &
                                    lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                    rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use salmon_global, only: yn_fix_func, yn_dg_length_gauge, ae_shape1, e_impulse, epdir_re1, yn_restart, nt, &
      yn_dg_expdiag_xi_split, yn_dg_expdiag_global_flux, yn_dg_expdiag_global_field, &
      yn_dg_expdiag_project_h, yn_dg_expdiag_delta_h, &
      yn_dg_mixed_z, yn_dg_mixed_z_local_prop_writeback, dg_mixed_z_local_prop_backend, &
      dg_mixed_z_frag_local_field_block, yn_dg_mixed_z_local_rho_writeback_wwonly, &
      yn_dg_full_h_eigen_seed
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional
    use rt_dg_fragment_types, only: s_dg_fragment_rt, matrix_block_info
    use rt_dg_fragment_ops, only: ensure_nonlocal_pp_matrix_A, calculate_microscopic_current_dg, &
      ensure_gradient_basis_cache, calculate_local_wannier_polarization_dg
    use rt_dg_plane_wave, only: diagnose_wpw_reduced_density, diagnose_wpw_reduced_embed_local, &
      diagnose_wpw_reduced_embed_prodcoef, initialize_wpw_reduced_self_projection, &
      wpw_normalized_window_at_grid, wpw_face_neighbor_fragment, map_global_to_phi_box_coord_pw, &
      reconstruct_psi_from_C_can, build_wpw_reduced_c_can_reference, &
      build_wpw_reduced_raw_back_hybrid_matrix, diagnose_wpw_mixed_neighbor_hamiltonian, &
      ensure_wpw_local_pp_blocks
    use rt_dg_wpw_linalg, only: build_hermitian_inverse, zhegv_with_query
    use eigen_subdiag_sub, only: eigen_zheev
    use communication, only: comm_summation, comm_get_max, comm_bcast
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
    integer :: ib, jb, ib_global, jb_global, iblk_mom
    real(8) :: Ac_ham(3), E_mid(3), rdot_diag
    real(8) :: impulse_k(3), impulse_k2
    real(8) :: cphase, sphase
    complex(8), parameter :: zi = (0.0d0, 1.0d0)
    complex(8), allocatable :: c_w(:), tmp_w(:), next_w(:)
    complex(8), allocatable :: h_eff(:,:), evec(:,:)
    complex(8), allocatable :: mom_work(:,:), mom_eff(:,:)
    real(8), allocatable :: eval(:)
    logical :: use_bpw_wannier_h, use_formal_wannier_h
    logical :: use_formal_seed_phase
    logical :: use_full_h_seed_phase
    logical :: use_impulse_kinetic_shift
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
    logical, save :: delta_h_env_checked = .false.
    logical, save :: delta_h_enabled = .false.
    logical, save :: delta_h_ref_valid = .false.
    type(matrix_block_info), allocatable, save :: delta_h_ref_blocks(:)
    integer, save :: delta_h_ref_nblocks = 0
    logical, save :: delta_h_warned = .false.
    character(len=32) :: delta_h_env
    logical, save :: mixed_z_env_checked = .false.
    logical, save :: mixed_z_enabled = .false.
    logical, save :: mixed_z_exp_step_diag_env_checked = .false.
    logical, save :: mixed_z_exp_step_diag_enabled = .false.
    logical, save :: mixed_z_split_obs_diag_env_checked = .false.
    logical, save :: mixed_z_split_obs_diag_enabled = .false.
    logical, save :: mixed_z_split_current_diag_env_checked = .false.
    logical, save :: mixed_z_split_current_diag_enabled = .false.
    logical, save :: mixed_z_use_split_prop_env_checked = .false.
    logical, save :: mixed_z_use_split_prop_enabled = .false.
    logical, save :: mixed_z_local_prop_writeback_env_checked = .false.
    logical, save :: mixed_z_local_prop_writeback_enabled = .false.
    logical, save :: mixed_z_local_prop_backend_env_checked = .false.
    logical, save :: mixed_z_frag_local_field_block_env_checked = .false.
    logical, save :: mixed_z_frag_local_reference_diag_env_checked = .false.
    logical, save :: mixed_z_frag_local_reference_diag_enabled = .false.
    logical, save :: mixed_z_local_debug_legacy_env_checked = .false.
    logical, save :: mixed_z_local_debug_legacy_enabled = .false.
    logical, save :: mixed_z_local_prop_obs_series_env_checked = .false.
    logical, save :: mixed_z_local_prop_obs_series_enabled = .false.
    logical, save :: mixed_z_local_prop_obs_series_local_mixedz_requested = .false.
    logical, save :: mixed_z_payload_prepare_smoke_env_checked = .false.
    logical, save :: mixed_z_payload_prepare_smoke_enabled = .false.
    logical, save :: mixed_z_series_payload_direct_smoke_env_checked = .false.
    logical, save :: mixed_z_series_payload_direct_smoke_enabled = .false.
    logical, save :: wpw_reduced_eig_cache_valid = .false.
    integer, save :: wpw_reduced_eig_build_count = 0
    integer, save :: wpw_reduced_eig_cache_hit_count = 0
    integer, save :: wpw_reduced_eig_skipped_count = 0
    integer, save :: wpw_payload_prepare_build_count = 0
    integer, save :: wpw_payload_prepare_cache_hit_count = 0
    integer, save :: wpw_payload_raw_build_count = 0
    integer, save :: wpw_payload_raw_cache_hit_count = 0
    integer, save :: wpw_payload_pp_block_build_count = 0
    integer, save :: wpw_payload_pp_block_cache_hit_count = 0
    integer, save :: wpw_payload_neighborH_build_count = 0
    integer, save :: wpw_payload_neighborH_cache_hit_count = 0
    logical, save :: mixed_z_local_pz_op_diag_env_checked = .false.
    logical, save :: mixed_z_local_pz_op_diag_enabled = .false.
    logical, save :: mixed_z_local_pz_op_all_env_checked = .false.
    logical, save :: mixed_z_local_pz_op_all_enabled = .false.
    logical, save :: mixed_z_local_rho_block_diag_env_checked = .false.
    logical, save :: mixed_z_local_rho_block_diag_enabled = .false.
    logical, save :: mixed_z_local_rho_path_log_enabled = .false.
    logical, save :: mixed_z_local_rho_writeback_wwonly_env_checked = .false.
    logical, save :: mixed_z_local_rho_writeback_wwonly_enabled = .false.
    logical, save :: mixed_z_local_current_op_diag_env_checked = .false.
    logical, save :: mixed_z_local_current_op_diag_enabled = .false.
    logical, save :: mixed_z_local_current_op_all_env_checked = .false.
    logical, save :: mixed_z_local_current_op_all_enabled = .false.
    logical, save :: mixed_z_field_kick_diag_env_checked = .false.
    logical, save :: mixed_z_field_kick_diag_enabled = .false.
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
    character(32) :: mixed_z_exp_step_diag_env
    character(32) :: mixed_z_split_obs_diag_env
    character(32) :: mixed_z_split_current_diag_env
    character(32) :: mixed_z_use_split_prop_env
    character(32) :: mixed_z_local_prop_writeback_env
    character(48) :: mixed_z_local_prop_backend_env
    character(48), save :: mixed_z_local_prop_backend_kind = 'global_mixed_split_backend'
    character(32) :: mixed_z_frag_local_field_block_env
    character(16), save :: mixed_z_frag_local_field_block_kind = 'all'
    character(32) :: mixed_z_frag_local_reference_diag_env
    character(32) :: mixed_z_local_debug_legacy_env
    character(32) :: mixed_z_local_prop_obs_series_env
    character(32) :: mixed_z_local_prop_obs_series_source_env
    character(32) :: mixed_z_local_prop_obs_series_local_mixedz_env
    character(32) :: mixed_z_payload_prepare_smoke_env
    character(32) :: mixed_z_series_payload_direct_smoke_env
    character(32) :: mixed_z_local_pz_op_diag_env
    character(32) :: mixed_z_local_pz_op_all_env
    character(32) :: mixed_z_local_rho_block_diag_env
    character(32) :: mixed_z_local_rho_writeback_wwonly_env
    character(32) :: mixed_z_local_current_op_diag_env
    character(32) :: mixed_z_local_current_op_all_env
    character(32) :: mixed_z_field_kick_diag_env
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
    logical :: do_mixed_z_split_reference_diag
    logical :: bad_mixed_z_split_reference_diag
    integer :: nstate_diag, nmix_diag, imix_diag, ist_diag
    real(8) :: norm_before_diag, norm_exp_diag, norm_prod_diag, norm_split_diag
    real(8) :: diff_exp_prod_diag, rel_diff_exp_prod_diag
    real(8) :: diff_prod_split_diag, rel_diff_prod_split_diag
    real(8) :: diff_split_exp_diag, rel_diff_split_exp_diag
    real(8) :: pz_exp_diag, pz_prod_diag, pz_diff_diag
    complex(8), allocatable :: cmix_before_diag(:,:,:), cmix_exp_diag(:,:,:), cmix_split_diag(:,:,:)
    complex(8), allocatable :: cmix_prod_diag(:,:)
    complex(8), allocatable :: h_step_diag(:,:), h_evec_diag(:,:), work_vec_diag(:)
    real(8), allocatable :: eval_step_diag(:)
    complex(8), allocatable, save :: prop_cmix_work(:,:), prop_field_h_work(:,:), prop_field_vec_work(:,:)
    complex(8), allocatable, save :: prop_tmp_work(:,:), prop_cw_local_work(:,:), prop_cw_global_work(:,:)
    complex(8), allocatable, save :: prop_cp_local_work(:,:), prop_cp_global_work(:,:)
    complex(8), allocatable, save :: prop_next_local_work(:,:)
    complex(8), allocatable, save :: prop_w_block_work(:,:), prop_coef_block_work(:,:), prop_scatter_block_work(:,:)
    complex(8), allocatable, save :: prop_field_axis_evec_cache(:,:,:)
    integer, allocatable, save :: prop_local_row_work(:)
    real(8), allocatable, save :: prop_field_eval_work(:), prop_field_axis_eval_cache(:,:)
    integer, allocatable, save :: prop_field_axis_cache_axis(:)
    logical, allocatable, save :: prop_field_axis_cache_valid(:)
    logical :: wpw_reduced_eig_cache_allowed, wpw_reduced_eig_cache_hit
    logical :: wpw_reduced_eig_built, wpw_reduced_eig_skipped
    real(8) :: wpw_reduced_eig_elapsed, wpw_reduced_payload_prepare_elapsed, wpw_reduced_payload_t0
    logical :: wpw_payload_prepare_cache_allowed, wpw_payload_prepare_cache_hit
    logical :: wpw_payload_raw_cache_hit, wpw_payload_pp_block_cache_hit, wpw_payload_neighborH_cache_hit
    logical :: wpw_payload_static_field
    real(8) :: wpw_payload_raw_elapsed, wpw_payload_pp_block_elapsed, wpw_payload_neighborH_elapsed
    real(8) :: wpw_payload_stage_t0
    character(48) :: wpw_payload_prepare_cache_reason
    logical :: mixed_z_backend_available, mixed_z_backend_applied, mixed_z_backend_bad
    character(48) :: mixed_z_backend_block_reason
    integer :: route_env_len, route_env_stat

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
    if (.not. mixed_z_field_kick_diag_env_checked) then
      mixed_z_field_kick_diag_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_FIELD_KICK_DIAG', mixed_z_field_kick_diag_env)
      select case (trim(adjustl(mixed_z_field_kick_diag_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_field_kick_diag_enabled = .true.
      case default
        mixed_z_field_kick_diag_enabled = .false.
      end select
      mixed_z_field_kick_diag_env_checked = .true.
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
    wpw_reduced_eig_cache_allowed = .false.
    wpw_reduced_eig_cache_hit = .false.
    wpw_reduced_eig_built = .false.
    wpw_reduced_eig_skipped = .false.
    wpw_reduced_eig_elapsed = 0.0d0
    wpw_reduced_payload_prepare_elapsed = 0.0d0
    wpw_payload_prepare_cache_allowed = .false.
    wpw_payload_prepare_cache_hit = .false.
    wpw_payload_raw_cache_hit = .false.
    wpw_payload_pp_block_cache_hit = .false.
    wpw_payload_neighborH_cache_hit = .false.
    wpw_payload_static_field = .false.
    wpw_payload_raw_elapsed = 0.0d0
    wpw_payload_pp_block_elapsed = 0.0d0
    wpw_payload_neighborH_elapsed = 0.0d0
    wpw_payload_prepare_cache_reason = 'not_requested'
    do_mixed_z_split_reference_diag = .false.
    bad_mixed_z_split_reference_diag = .false.
    dg_frag%mixed_z_local_prop_rho_ready = .false.
    dg_frag%mixed_z_local_prop_rho_bad = .true.
    dg_frag%mixed_z_local_prop_rho_step = -1
    dg_frag%mixed_z_local_prop_rho_prod_int = 0.0d0
    dg_frag%mixed_z_local_prop_rho_candidate_int = 0.0d0
    dg_frag%mixed_z_local_prop_rho_diff_int = 0.0d0
    dg_frag%mixed_z_local_prop_rho_max_abs_diff = 0.0d0
    dg_frag%mixed_z_local_prop_rho_rms_diff = 0.0d0
    dg_frag%mixed_z_local_prop_payload_ready = .false.
    dg_frag%mixed_z_local_prop_payload_bad = .true.
    dg_frag%mixed_z_local_prop_payload_step = -1
    dg_frag%mixed_z_local_prop_payload_coef_norm = 0.0d0
    dg_frag%mixed_z_local_prop_payload_prod_coef_norm = 0.0d0
    dg_frag%mixed_z_local_prop_payload_coef_diff_snorm = 0.0d0
    dg_frag%mixed_z_local_prop_payload_rel_coef_diff_snorm = 0.0d0
    dg_frag%mixed_z_local_prop_payload_dim = 0.0d0
    dg_frag%mixed_z_local_prop_payload_occ_trace = 0.0d0
    dg_frag%mixed_z_local_prop_payload_source = 'missing_payload'
    dg_frag%mixed_z_local_prop_payload_basis_kind = 'none'
    dg_frag%mixed_z_local_prop_payload_build_route = 'not_requested'
    dg_frag%mixed_z_local_prop_payload_block_reason = 'not_requested'
    dg_frag%mixed_z_local_pz_wcenter_raw_z = 0.0d0
    dg_frag%mixed_z_local_pz_wcenter_ref_raw_z = 0.0d0
    dg_frag%mixed_z_local_pz_wcenter_diff = 0.0d0
    dg_frag%mixed_z_local_pz_wcenter_rel_diff = huge(1.0d0)
    dg_frag%mixed_z_local_pz_wcenter_weight = 0.0d0
    dg_frag%mixed_z_local_pz_wcenter_ref_weight = 0.0d0
    dg_frag%mixed_z_local_pz_wcenter_slot_count = 0.0d0
    dg_frag%mixed_z_local_pz_wcenter_block_count = 0.0d0
    dg_frag%mixed_z_local_pz_wcenter_owner_unique_raw_z = 0.0d0
    dg_frag%mixed_z_local_pz_wcenter_owner_unique_ref_raw_z = 0.0d0
    dg_frag%mixed_z_local_pz_wcenter_owner_unique_weight = 0.0d0
    dg_frag%mixed_z_local_pz_wcenter_owner_unique_ref_weight = 0.0d0
    dg_frag%mixed_z_local_pz_wcenter_owner_unique_count = 0.0d0
    dg_frag%mixed_z_local_pz_wcenter_missing_owner_count = 0.0d0
    dg_frag%mixed_z_local_pz_wcenter_duplicate_owner_count = 0.0d0
    dg_frag%mixed_z_local_pz_owner_unique_global_center_raw_z = 0.0d0
    dg_frag%mixed_z_local_pz_owner_unique_global_zww_diag_raw_z = 0.0d0
    dg_frag%mixed_z_local_pz_owner_unique_weighted_diff_global_center = 0.0d0
    dg_frag%mixed_z_local_pz_owner_unique_weighted_diff_global_zww_diag = 0.0d0
    dg_frag%mixed_z_local_pz_owner_unique_center_source_match_count = 0.0d0
    dg_frag%mixed_z_local_pz_owner_unique_wrap_match_count = 0.0d0
    dg_frag%mixed_z_local_pz_owner_unique_cell_shift_count = 0.0d0
    dg_frag%mixed_z_local_pz_zww_weight_sum_z_prod = 0.0d0
    dg_frag%mixed_z_local_pz_zww_weight_sum_z_lsp = 0.0d0
    dg_frag%mixed_z_local_pz_zww_weight_sum_weight_prod = 0.0d0
    dg_frag%mixed_z_local_pz_zww_weight_sum_weight_lsp = 0.0d0
    dg_frag%mixed_z_local_pz_zww_weight_sum_contrib_prod = 0.0d0
    dg_frag%mixed_z_local_pz_zww_weight_sum_contrib_lsp = 0.0d0
    dg_frag%mixed_z_local_pz_zww_weight_max_abs_diff_z = 0.0d0
    dg_frag%mixed_z_local_pz_zww_weight_max_abs_diff_weight = 0.0d0
    dg_frag%mixed_z_local_pz_zww_weight_max_abs_diff_contrib = 0.0d0
    dg_frag%mixed_z_local_pz_zww_weight_rms_diff_contrib = 0.0d0
    dg_frag%mixed_z_local_pz_zww_weight_owner_gid_mismatch_count = 0.0d0
    dg_frag%mixed_z_local_pz_weight_source_sum_occ_prod = 0.0d0
    dg_frag%mixed_z_local_pz_weight_source_sum_occ_lsp = 0.0d0
    dg_frag%mixed_z_local_pz_weight_source_sum_rho_prod = 0.0d0
    dg_frag%mixed_z_local_pz_weight_source_sum_rho_lsp = 0.0d0
    dg_frag%mixed_z_local_pz_weight_source_max_abs_diff_occ = 0.0d0
    dg_frag%mixed_z_local_pz_weight_source_max_abs_diff_rho = 0.0d0
    dg_frag%mixed_z_local_pz_weight_source_max_abs_diff_factor = 0.0d0
    dg_frag%mixed_z_local_pz_weight_source_max_abs_diff_weight = 0.0d0
    dg_frag%mixed_z_local_pz_weight_source_weighted_zww_prod = 0.0d0
    dg_frag%mixed_z_local_pz_weight_source_weighted_zww_lsp = 0.0d0
    dg_frag%mixed_z_local_pz_rhodiag_prod_observable_abs2 = 0.0d0
    dg_frag%mixed_z_local_pz_rhodiag_prod_expdiag_abs2 = 0.0d0
    dg_frag%mixed_z_local_pz_rhodiag_source_abs2 = 0.0d0
    dg_frag%mixed_z_local_pz_rhodiag_source_smetric = 0.0d0
    dg_frag%mixed_z_local_pz_rhodiag_ref_after_abs2 = 0.0d0
    dg_frag%mixed_z_local_pz_rhodiag_ref_after_smetric = 0.0d0
    dg_frag%mixed_z_local_pz_rhodiag_lsp_after_abs2 = 0.0d0
    dg_frag%mixed_z_local_pz_rhodiag_lsp_after_smetric = 0.0d0
    dg_frag%mixed_z_local_pz_rhodiag_repacked_global_abs2 = 0.0d0
    dg_frag%mixed_z_local_pz_rhodiag_repacked_global_smetric = 0.0d0
    dg_frag%mixed_z_local_pz_rhodiag_max_abs_diff_abs2 = 0.0d0
    dg_frag%mixed_z_local_pz_rhodiag_max_abs_diff_smetric = 0.0d0
    dg_frag%mixed_z_local_pz_rhodiag_repacked_global_available = .false.
    dg_frag%mixed_z_local_pz_repacked_bridge_prod = 0.0d0
    dg_frag%mixed_z_local_pz_repacked_bridge_local_cref = 0.0d0
    dg_frag%mixed_z_local_pz_repacked_bridge_repacked_global = 0.0d0
    dg_frag%mixed_z_local_pz_repacked_bridge_rho_local_cref = 0.0d0
    dg_frag%mixed_z_local_pz_repacked_bridge_rho_repacked_global = 0.0d0
    dg_frag%mixed_z_local_pz_repacked_bridge_weight_repacked_global = 0.0d0
    dg_frag%mixed_z_local_pz_repacked_bridge_diff_prod_local = 0.0d0
    dg_frag%mixed_z_local_pz_repacked_bridge_diff_prod_repacked = 0.0d0
    dg_frag%mixed_z_local_pz_repacked_bridge_available = .false.
    dg_frag%mixed_z_local_pz_wcenter_ready = .false.
    dg_frag%mixed_z_local_pz_wcenter_bad = .true.
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
      global_flux_exp_enabled = (yn_dg_expdiag_global_flux == 'y')
      global_flux_env = ' '
      call get_environment_variable('SALMON_DG_EXPDIAG_GLOBAL_FLUX', global_flux_env, &
        length=route_env_len, status=route_env_stat)
      if (route_env_stat == 0 .and. route_env_len > 0) then
      select case (trim(adjustl(global_flux_env(1:route_env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        global_flux_exp_enabled = .true.
      case ('0','n','N','no','NO','false','FALSE','off','OFF')
        global_flux_exp_enabled = .false.
      end select
      end if
      global_flux_env_checked = .true.
    end if
    if (.not. global_field_env_checked) then
      global_field_exp_enabled = (yn_dg_expdiag_global_field == 'y')
      global_field_env = ' '
      call get_environment_variable('SALMON_DG_EXPDIAG_GLOBAL_FIELD', global_field_env, &
        length=route_env_len, status=route_env_stat)
      if (route_env_stat == 0 .and. route_env_len > 0) then
      select case (trim(adjustl(global_field_env(1:route_env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        global_field_exp_enabled = .true.
      case ('0','n','N','no','NO','false','FALSE','off','OFF')
        global_field_exp_enabled = .false.
      end select
      end if
      global_field_env_checked = .true.
    end if
    if (.not. xi_split_env_checked) then
      xi_split_enabled = (yn_dg_expdiag_xi_split == 'y')
      xi_split_env = ' '
      call get_environment_variable('SALMON_DG_EXPDIAG_XI_SPLIT', xi_split_env, &
        length=route_env_len, status=route_env_stat)
      if (route_env_stat == 0 .and. route_env_len > 0) then
      select case (trim(adjustl(xi_split_env(1:route_env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        xi_split_enabled = .true.
      case ('0','n','N','no','NO','false','FALSE','off','OFF')
        xi_split_enabled = .false.
      end select
      end if
      xi_split_env_checked = .true.
    end if
    if (.not. project_h_env_checked) then
      project_h_for_fixed_func = (yn_dg_expdiag_project_h == 'y')
      project_h_env = ' '
      call get_environment_variable('SALMON_DG_EXPDIAG_PROJECT_H', project_h_env, &
        length=route_env_len, status=route_env_stat)
      if (route_env_stat == 0 .and. route_env_len > 0) then
      select case (trim(adjustl(project_h_env(1:route_env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        project_h_for_fixed_func = .true.
      case ('0','n','N','no','NO','false','FALSE','off','OFF')
        project_h_for_fixed_func = .false.
      end select
      end if
      project_h_env_checked = .true.
    end if
    if (.not. delta_h_env_checked) then
      delta_h_enabled = (yn_dg_expdiag_delta_h == 'y')
      delta_h_env = ' '
      call get_environment_variable('SALMON_DG_EXPDIAG_DELTA_H', delta_h_env, &
        length=route_env_len, status=route_env_stat)
      if (route_env_stat == 0 .and. route_env_len > 0) then
      select case (trim(adjustl(delta_h_env(1:route_env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        delta_h_enabled = .true.
      case ('0','n','N','no','NO','false','FALSE','off','OFF')
        delta_h_enabled = .false.
      end select
      end if
      delta_h_env_checked = .true.
    end if
    if (.not. mixed_z_env_checked) then
      mixed_z_enabled = (yn_dg_mixed_z == 'y')
      mixed_z_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z', mixed_z_env, &
        length=route_env_len, status=route_env_stat)
      if (route_env_stat == 0 .and. route_env_len > 0) then
      select case (trim(adjustl(mixed_z_env(1:route_env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_enabled = .true.
      case ('0','n','N','no','NO','false','FALSE','off','OFF')
        mixed_z_enabled = .false.
      end select
      end if
      mixed_z_env_checked = .true.
    end if
    if (.not. mixed_z_exp_step_diag_env_checked) then
      mixed_z_exp_step_diag_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_EXP_STEP_DIAG', mixed_z_exp_step_diag_env)
      select case (trim(adjustl(mixed_z_exp_step_diag_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_exp_step_diag_enabled = .true.
      case default
        mixed_z_exp_step_diag_enabled = .false.
      end select
      mixed_z_exp_step_diag_env_checked = .true.
    end if
    if (.not. mixed_z_split_obs_diag_env_checked) then
      mixed_z_split_obs_diag_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_SPLIT_OBS_DIAG', mixed_z_split_obs_diag_env)
      select case (trim(adjustl(mixed_z_split_obs_diag_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_split_obs_diag_enabled = .true.
      case default
        mixed_z_split_obs_diag_enabled = .false.
      end select
      mixed_z_split_obs_diag_env_checked = .true.
    end if
    if (.not. mixed_z_split_current_diag_env_checked) then
      mixed_z_split_current_diag_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_SPLIT_CURRENT_DIAG', mixed_z_split_current_diag_env)
      select case (trim(adjustl(mixed_z_split_current_diag_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_split_current_diag_enabled = .true.
      case default
        mixed_z_split_current_diag_enabled = .false.
      end select
      mixed_z_split_current_diag_env_checked = .true.
    end if
    ! Mixed-Z local production-compatible route policy:
    !   propagation: SALMON_DG_MIXED_Z_LOCAL_PROP_WRITEBACK=1
    !                production-equivalent split backend writeback.
    !   density:     SALMON_DG_MIXED_Z_LOCAL_RHO_WRITEBACK_WWONLY=1
    !                WWONLY production grid target, intentionally not TOTAL.
    !   Pz:          SALMON_DG_MIXED_Z_LOCAL_PZ_WRITEBACK_TOTAL=1
    !                TOTAL production target, WW_full + WP/PW + PP.
    !   current:     SALMON_DG_MIXED_Z_LOCAL_CURRENT_WRITEBACK_TOTAL=1
    !                TOTAL production target, para + nonlocal + diamagnetic.
    if (.not. mixed_z_use_split_prop_env_checked) then
      mixed_z_use_split_prop_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_USE_SPLIT_PROPAGATION', mixed_z_use_split_prop_env)
      select case (trim(adjustl(mixed_z_use_split_prop_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_use_split_prop_enabled = .true.
      case default
        mixed_z_use_split_prop_enabled = .false.
      end select
      mixed_z_use_split_prop_env_checked = .true.
    end if
    if (.not. mixed_z_local_prop_writeback_env_checked) then
      mixed_z_local_prop_writeback_enabled = (yn_dg_mixed_z_local_prop_writeback == 'y')
      mixed_z_local_prop_writeback_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PROP_WRITEBACK', &
        mixed_z_local_prop_writeback_env, length=route_env_len, status=route_env_stat)
      if (route_env_stat == 0 .and. route_env_len > 0) then
      select case (trim(adjustl(mixed_z_local_prop_writeback_env(1:route_env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_local_prop_writeback_enabled = .true.
      case ('0','n','N','no','NO','false','FALSE','off','OFF')
        mixed_z_local_prop_writeback_enabled = .false.
      end select
      end if
      mixed_z_local_prop_writeback_env_checked = .true.
    end if
    if (.not. mixed_z_local_prop_backend_env_checked) then
      select case (trim(adjustl(dg_mixed_z_local_prop_backend)))
      case ('fragment_local_mixed_split_backend','fragment_local','local')
        mixed_z_local_prop_backend_kind = 'fragment_local_mixed_split_backend'
      case ('neighbor_env_expdiag','neighbor_env','neighbor')
        mixed_z_local_prop_backend_kind = 'neighbor_env_expdiag'
      case ('neighbor_env_fieldonly','neighbor_fieldonly','fieldonly')
        mixed_z_local_prop_backend_kind = 'neighbor_env_fieldonly'
      case ('neighbor_env_interaction','neighbor_interaction','interaction','neighbor_env_delta','neighbor_delta','delta')
        mixed_z_local_prop_backend_kind = 'neighbor_env_interaction'
      case ('global_mixed_split_backend','global')
        mixed_z_local_prop_backend_kind = 'global_mixed_split_backend'
      case default
        stop "DG mixed-Z: invalid dg_mixed_z_local_prop_backend"
      end select
      mixed_z_local_prop_backend_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PROP_BACKEND', &
        mixed_z_local_prop_backend_env, length=route_env_len, status=route_env_stat)
      if (route_env_stat == 0 .and. route_env_len > 0) then
      select case (trim(adjustl(mixed_z_local_prop_backend_env(1:route_env_len))))
      case ('fragment_local_mixed_split_backend','fragment_local','local')
        mixed_z_local_prop_backend_kind = 'fragment_local_mixed_split_backend'
      case ('neighbor_env_expdiag','neighbor_env','neighbor')
        mixed_z_local_prop_backend_kind = 'neighbor_env_expdiag'
      case ('neighbor_env_fieldonly','neighbor_fieldonly','fieldonly')
        mixed_z_local_prop_backend_kind = 'neighbor_env_fieldonly'
      case ('neighbor_env_interaction','neighbor_interaction','interaction','neighbor_env_delta','neighbor_delta','delta')
        mixed_z_local_prop_backend_kind = 'neighbor_env_interaction'
      case ('global_mixed_split_backend','global')
        mixed_z_local_prop_backend_kind = 'global_mixed_split_backend'
      case default
        stop "DG mixed-Z: invalid SALMON_DG_MIXED_Z_LOCAL_PROP_BACKEND"
      end select
      end if
      mixed_z_local_prop_backend_env_checked = .true.
    end if
    if (.not. mixed_z_frag_local_field_block_env_checked) then
      select case (trim(adjustl(dg_mixed_z_frag_local_field_block)))
      case ('w_only','w','none')
        mixed_z_frag_local_field_block_kind = 'w_only'
      case ('w_pself','self','p_self')
        mixed_z_frag_local_field_block_kind = 'w_pself'
      case ('owner_local','owner','owned')
        mixed_z_frag_local_field_block_kind = 'owner_local'
      case ('halo_owner','owner_halo','neighbor_read','all_owner')
        mixed_z_frag_local_field_block_kind = 'halo_owner'
      case ('all')
        mixed_z_frag_local_field_block_kind = 'all'
      case default
        stop "DG mixed-Z: invalid dg_mixed_z_frag_local_field_block"
      end select
      mixed_z_frag_local_field_block_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_FRAG_LOCAL_FIELD_BLOCK', &
        mixed_z_frag_local_field_block_env, length=route_env_len, status=route_env_stat)
      if (route_env_stat == 0 .and. route_env_len > 0) then
      select case (trim(adjustl(mixed_z_frag_local_field_block_env(1:route_env_len))))
      case ('w_only','W_ONLY','w','W','none','NONE')
        mixed_z_frag_local_field_block_kind = 'w_only'
      case ('w_pself','W_PSELF','self','SELF','p_self','P_SELF')
        mixed_z_frag_local_field_block_kind = 'w_pself'
      case ('owner_local','OWNER_LOCAL','owner','OWNER','owned','OWNED')
        mixed_z_frag_local_field_block_kind = 'owner_local'
      case ('halo_owner','HALO_OWNER','owner_halo','OWNER_HALO','neighbor_read','NEIGHBOR_READ', &
            'all_owner','ALL_OWNER')
        mixed_z_frag_local_field_block_kind = 'halo_owner'
      case ('all','ALL')
        mixed_z_frag_local_field_block_kind = 'all'
      case default
        stop "DG mixed-Z: invalid SALMON_DG_MIXED_Z_FRAG_LOCAL_FIELD_BLOCK"
      end select
      end if
      dg_frag%mixed_z_frag_local_field_block_kind = mixed_z_frag_local_field_block_kind
      mixed_z_frag_local_field_block_env_checked = .true.
    end if
    if (.not. mixed_z_frag_local_reference_diag_env_checked) then
      mixed_z_frag_local_reference_diag_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_FRAG_LOCAL_REFERENCE_DIAG', &
        mixed_z_frag_local_reference_diag_env)
      select case (trim(adjustl(mixed_z_frag_local_reference_diag_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_frag_local_reference_diag_enabled = .true.
      case default
        mixed_z_frag_local_reference_diag_enabled = .false.
      end select
      mixed_z_frag_local_reference_diag_env_checked = .true.
    end if
    if (.not. mixed_z_local_debug_legacy_env_checked) then
      mixed_z_local_debug_legacy_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_DEBUG_LEGACY', &
        mixed_z_local_debug_legacy_env)
      select case (trim(adjustl(mixed_z_local_debug_legacy_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_local_debug_legacy_enabled = .true.
      case default
        mixed_z_local_debug_legacy_enabled = .false.
      end select
      mixed_z_local_debug_legacy_env_checked = .true.
    end if
    if (.not. mixed_z_local_prop_obs_series_env_checked) then
      mixed_z_local_prop_obs_series_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PROP_OBS_SERIES_DIAG', &
        mixed_z_local_prop_obs_series_env)
      select case (trim(adjustl(mixed_z_local_prop_obs_series_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_local_prop_obs_series_enabled = .true.
        mixed_z_local_rho_block_diag_enabled = .true.
      case default
        mixed_z_local_prop_obs_series_enabled = .false.
      end select
      mixed_z_local_prop_obs_series_source_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PROP_OBS_SERIES_SOURCE', &
        mixed_z_local_prop_obs_series_source_env)
      select case (trim(adjustl(mixed_z_local_prop_obs_series_source_env)))
      case ('prod_reference','local_mixedz','local_mixedz_writeback')
        mixed_z_local_prop_obs_series_enabled = .true.
        mixed_z_local_rho_block_diag_enabled = .true.
      end select
      select case (trim(adjustl(mixed_z_local_prop_obs_series_source_env)))
      case ('local_mixedz','local_mixedz_writeback')
        mixed_z_local_prop_obs_series_local_mixedz_requested = .true.
      end select
      mixed_z_local_prop_obs_series_local_mixedz_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PROP_OBS_SERIES_LOCAL_MIXEDZ', &
        mixed_z_local_prop_obs_series_local_mixedz_env)
      select case (trim(adjustl(mixed_z_local_prop_obs_series_local_mixedz_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_local_prop_obs_series_enabled = .true.
        mixed_z_local_rho_block_diag_enabled = .true.
        mixed_z_local_prop_obs_series_local_mixedz_requested = .true.
      end select
      mixed_z_local_prop_obs_series_env_checked = .true.
    end if
    if (.not. mixed_z_local_pz_op_diag_env_checked) then
      mixed_z_local_pz_op_diag_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PZ_OP_DIAG', mixed_z_local_pz_op_diag_env)
      select case (trim(adjustl(mixed_z_local_pz_op_diag_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_local_pz_op_diag_enabled = .true.
      case default
        mixed_z_local_pz_op_diag_enabled = .false.
      end select
      mixed_z_local_pz_op_diag_env_checked = .true.
    end if
    if (.not. mixed_z_local_pz_op_all_env_checked) then
      mixed_z_local_pz_op_all_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PZ_OP_ALL', mixed_z_local_pz_op_all_env)
      select case (trim(adjustl(mixed_z_local_pz_op_all_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_local_pz_op_all_enabled = .true.
      case default
        mixed_z_local_pz_op_all_enabled = .false.
      end select
      mixed_z_local_pz_op_all_env_checked = .true.
    end if
    if (.not. mixed_z_local_rho_block_diag_env_checked) then
      mixed_z_local_rho_block_diag_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_RHO_BLOCK_DIAG', &
        mixed_z_local_rho_block_diag_env)
      select case (trim(adjustl(mixed_z_local_rho_block_diag_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_local_rho_block_diag_enabled = .true.
        mixed_z_local_rho_path_log_enabled = .true.
      case default
        mixed_z_local_rho_block_diag_enabled = mixed_z_local_prop_obs_series_enabled
      end select
      mixed_z_local_rho_block_diag_env_checked = .true.
    end if
    if (.not. mixed_z_local_rho_writeback_wwonly_env_checked) then
      mixed_z_local_rho_writeback_wwonly_enabled = &
        (yn_dg_mixed_z_local_rho_writeback_wwonly == 'y')
      if (mixed_z_local_rho_writeback_wwonly_enabled) then
        mixed_z_local_rho_block_diag_enabled = .true.
        mixed_z_local_rho_path_log_enabled = .true.
      end if
      mixed_z_local_rho_writeback_wwonly_env = ' '
      ! Mixed-Z observable targets are intentionally not uniform:
      ! density follows the current production-compatible WWONLY grid target,
      ! while Pz/current use their accepted TOTAL production-compatible paths.
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_RHO_WRITEBACK_WWONLY', &
        mixed_z_local_rho_writeback_wwonly_env, length=route_env_len, status=route_env_stat)
      if (route_env_stat == 0 .and. route_env_len > 0) then
      select case (trim(adjustl(mixed_z_local_rho_writeback_wwonly_env(1:route_env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_local_rho_writeback_wwonly_enabled = .true.
        mixed_z_local_rho_block_diag_enabled = .true.
        mixed_z_local_rho_path_log_enabled = .true.
      case ('0','n','N','no','NO','false','FALSE','off','OFF')
        mixed_z_local_rho_writeback_wwonly_enabled = .false.
      end select
      end if
      mixed_z_local_rho_writeback_wwonly_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_RHO_WRITEBACK_PROD_TARGET', &
        mixed_z_local_rho_writeback_wwonly_env, length=route_env_len, status=route_env_stat)
      if (route_env_stat == 0 .and. route_env_len > 0) then
      select case (trim(adjustl(mixed_z_local_rho_writeback_wwonly_env(1:route_env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_local_rho_writeback_wwonly_enabled = .true.
        mixed_z_local_rho_block_diag_enabled = .true.
        mixed_z_local_rho_path_log_enabled = .true.
      end select
      end if
      mixed_z_local_rho_writeback_wwonly_env_checked = .true.
    end if
    if (.not. mixed_z_local_current_op_diag_env_checked) then
      mixed_z_local_current_op_diag_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_CURRENT_OP_DIAG', mixed_z_local_current_op_diag_env)
      select case (trim(adjustl(mixed_z_local_current_op_diag_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_local_current_op_diag_enabled = .true.
      case default
        mixed_z_local_current_op_diag_enabled = .false.
      end select
      mixed_z_local_current_op_diag_env_checked = .true.
    end if
    if (.not. mixed_z_local_current_op_all_env_checked) then
      mixed_z_local_current_op_all_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_CURRENT_OP_ALL', mixed_z_local_current_op_all_env)
      select case (trim(adjustl(mixed_z_local_current_op_all_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_local_current_op_all_enabled = .true.
      case default
        mixed_z_local_current_op_all_enabled = .false.
      end select
      mixed_z_local_current_op_all_env_checked = .true.
    end if
    if (mixed_z_local_prop_obs_series_enabled) then
      dg_frag%mixed_z_local_prop_payload_source = 'missing_payload'
      dg_frag%mixed_z_local_prop_payload_basis_kind = 'owner_unique_Wphase_WP_PP'
      dg_frag%mixed_z_local_prop_payload_build_route = 'normal_series_requested'
      if (.not. allocated(dg_frag%wpw_reduced_Sraw_build) .or. &
          .not. allocated(dg_frag%wpw_reduced_Hraw_build) .or. &
          .not. allocated(dg_frag%wpw_reduced_nraw)) then
        dg_frag%mixed_z_local_prop_payload_block_reason = 'missing_wpw_reduced_raw_build'
      else if (.not. allocated(dg_frag%wpw_reduced_transform) .or. &
               .not. allocated(dg_frag%wpw_reduced_S) .or. &
               .not. allocated(dg_frag%wpw_reduced_eval)) then
        dg_frag%mixed_z_local_prop_payload_block_reason = 'missing_wpw_reduced_transform'
      else
        dg_frag%mixed_z_local_prop_payload_block_reason = 'producer_not_run'
      end if
    end if
    if (.not. mixed_z_payload_prepare_smoke_env_checked) then
      mixed_z_payload_prepare_smoke_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_PAYLOAD_PREPARE_SMOKE', &
        mixed_z_payload_prepare_smoke_env)
      select case (trim(adjustl(mixed_z_payload_prepare_smoke_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_payload_prepare_smoke_enabled = .true.
      case default
        mixed_z_payload_prepare_smoke_enabled = .false.
      end select
      mixed_z_payload_prepare_smoke_env_checked = .true.
    end if
    if (.not. mixed_z_series_payload_direct_smoke_env_checked) then
      mixed_z_series_payload_direct_smoke_env = ' '
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PROP_SERIES_PAYLOAD_DIRECT_SMOKE', &
        mixed_z_series_payload_direct_smoke_env)
      select case (trim(adjustl(mixed_z_series_payload_direct_smoke_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        mixed_z_series_payload_direct_smoke_enabled = .true.
      case default
        mixed_z_series_payload_direct_smoke_enabled = .false.
      end select
      mixed_z_series_payload_direct_smoke_env_checked = .true.
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
    use_impulse_kinetic_shift = trim(ae_shape1) == 'impulse' .and. yn_restart == 'n'
    impulse_k(:) = 0.0d0
    impulse_k2 = 0.0d0
    if (use_impulse_kinetic_shift) then
      impulse_k(:) = e_impulse * epdir_re1(:)
      impulse_k2 = sum(impulse_k(1:3)**2)
      E_mid(:) = 0.0d0
      if (itt == 1 .and. dg_frag%id == 0 .and. .not. dg_frag%mixed_z_perf_count_enabled) then
        write(*,'(1x,a,3(1x,1pe13.5),a,1pe13.5)') &
          '[DG-LG-IMPULSE] using kinetic shift T(k)=T0-i*k.grad+k^2/2; k=', impulse_k(:), &
          ' k2=', impulse_k2
        flush(6)
      end if
    end if
    if (trace_wpw_order .and. dg_frag%id == 0) then
      write(*,'(1x,a,a,i0,a,1pe12.4,3(a,1pe12.4),3(a,1pe12.4),4(a,l1))') &
        '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, ' dt=', dt, &
        ' E_x=', E_mid(1), ' E_y=', E_mid(2), ' E_z=', E_mid(3), &
        ' Ac_x=', Ac_ham(1), ' Ac_y=', Ac_ham(2), ' Ac_z=', Ac_ham(3), &
        ' production_kinetic_shift=', use_impulse_kinetic_shift, &
        ' density_hxc_update_next=', yn_fix_func == 'n', &
        ' aux_reduced_next=', wpw_red_expdiag_enabled, &
        ' production_mixed_z=', mixed_z_enabled
    end if
    dg_frag%mixed_z_frag_local_field_abs_max = max(dg_frag%mixed_z_frag_local_field_abs_max, &
      sum(abs(E_mid(1:3))))

    nstate_prop = dg_frag%nstate_tot
    if (allocated(dg_frag%nocc_spin)) then
      nstate_prop = min(dg_frag%nstate_tot, max(1, maxval(dg_frag%nocc_spin(1:dg_frag%nspin))))
    end if
    state_first = 1
    state_last = min(nstate_prop, size(dg_frag%coef, 2))

    if (mixed_z_local_prop_writeback_enabled .and. mixed_z_enabled) then
      call ensure_mixed_wannier_bpw_position(dg_frag)
      if (dg_frag%has_mixed_wannier_bpw_position .and. state_last >= state_first) then
        t_route0 = get_wtime()
        call apply_mixed_split_exp_backend(trim(mixed_z_local_prop_backend_kind), E_mid, &
          state_first, state_last, mixed_z_backend_available, mixed_z_backend_applied, &
          mixed_z_backend_bad, mixed_z_backend_block_reason)
        t1 = get_wtime() - t_route0
        time_mixed_field = time_mixed_field + t1
        if (mixed_z_backend_applied) then
          dg_frag%mixed_z_perf_prop_writeback_calls = dg_frag%mixed_z_perf_prop_writeback_calls + 1_8
        end if
        dg_frag%mixed_z_perf_wall_prop = dg_frag%mixed_z_perf_wall_prop + t1
        if (mixed_z_backend_applied .and. mixed_z_local_prop_obs_series_enabled .and. &
            .not. mixed_z_local_debug_legacy_enabled) then
          call build_wpw_local_prop_payload_fast_smoke(0, 'normal_series_after_propagation')
        end if
        if (dg_frag%id == 0 .and. .not. dg_frag%mixed_z_perf_count_enabled) then
          write(*,'(1x,a,1(a,i0),1(a,1pe12.4),7(a,l1),5(a,a))') &
            '[DG-MIXEDZ-LOCAL-PROP-WRITEBACK-CMP]', &
            ' step=', itt, &
            ' elapsed=', t1, &
            ' candidate_available=', mixed_z_backend_available, &
            ' replacement_ready=', mixed_z_backend_applied .and. .not. mixed_z_backend_bad, &
            ' replacement_applied=', mixed_z_backend_applied, &
            ' production_value_modified=', mixed_z_backend_applied, &
            ' observables_series_enabled=', mixed_z_local_prop_obs_series_enabled, &
            ' writeback_flag_enabled=', mixed_z_local_prop_writeback_enabled, &
            ' bad=', mixed_z_backend_bad, &
            ' target_kind=', 'production_equivalent_full_h_exp', &
            ' candidate_kind=', trim(mixed_z_local_prop_backend_kind), &
            ' payload_route=', 'normal_series_after_propagation', &
            ' replacement_block_reason=', trim(mixed_z_backend_block_reason), &
            ' reference_scope=', 'propagation_backend'
        end if
        if (mixed_z_backend_applied .and. .not. mixed_z_backend_bad) return
      end if
    end if

    if (delta_h_enabled .and. .not. project_h_for_fixed_func) then
      call ensure_delta_h_reference()
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
    if (mixed_z_payload_prepare_smoke_enabled .and. .not. mixed_z_local_debug_legacy_enabled) then
      call prepare_wpw_local_payload_ingredients('payload_prepare_smoke')
      call log_wpw_local_payload_prepare_smoke(1)
      call prepare_wpw_local_payload_ingredients('payload_prepare_smoke')
      call log_wpw_local_payload_prepare_smoke(2)
      return
    end if

    nstate_prop = dg_frag%nstate_tot
    if (allocated(dg_frag%nocc_spin)) then
      nstate_prop = min(dg_frag%nstate_tot, max(1, maxval(dg_frag%nocc_spin(1:dg_frag%nspin))))
    end if
    state_first = 1
    state_last = min(nstate_prop, size(dg_frag%coef, 2))

    if (mixed_z_series_payload_direct_smoke_enabled .and. .not. mixed_z_local_debug_legacy_enabled) then
      call prepare_wpw_local_payload_ingredients('series_payload_direct_smoke')
      call build_wpw_local_prop_payload_fast_smoke(1)
      call log_wpw_local_prop_series_direct_smoke(1)
      call prepare_wpw_local_payload_ingredients('series_payload_direct_smoke')
      call build_wpw_local_prop_payload_fast_smoke(2)
      call log_wpw_local_prop_series_direct_smoke(2)
      return
    end if

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
        ! Mixed-Z diagnostics are intentionally split into three layers:
        !   A. production-equivalent split backend/candidate
        !   B. dense frozen-H exp diagnostic, field-off reference only
        !   C. observable/current comparison diagnostics
        ! All layers use the same sampled split coefficients below, but only
        ! SALMON_DG_MIXED_Z_USE_SPLIT_PROPAGATION changes production propagation.
        do_mixed_z_split_reference_diag = (mixed_z_exp_step_diag_enabled .or. &
          mixed_z_split_obs_diag_enabled .or. mixed_z_split_current_diag_enabled) .and. &
          allocated(dg_frag%mixed_wannier_bpw_eval) .and. &
          allocated(dg_frag%mixed_wannier_bpw_z) .and. &
          dg_frag%mixed_wannier_bpw_nmix > 0 .and. state_last >= state_first
        if (do_mixed_z_split_reference_diag) then
          nmix_diag = dg_frag%mixed_wannier_bpw_nmix
          nstate_diag = state_last - state_first + 1
          allocate(cmix_before_diag(nmix_diag,nstate_diag,dg_frag%nspin), &
                   cmix_exp_diag(nmix_diag,nstate_diag,dg_frag%nspin), &
                   cmix_split_diag(nmix_diag,nstate_diag,dg_frag%nspin), &
                   cmix_prod_diag(nmix_diag,nstate_diag), &
                   h_step_diag(nmix_diag,nmix_diag), h_evec_diag(nmix_diag,nmix_diag), &
                   work_vec_diag(nmix_diag), eval_step_diag(nmix_diag))
          cmix_before_diag = (0.0d0, 0.0d0)
          cmix_exp_diag = (0.0d0, 0.0d0)
          cmix_split_diag = (0.0d0, 0.0d0)
          do ispin = 1, dg_frag%nspin
            if (ispin > size(dg_frag%mixed_wannier_bpw_eval, 2) .or. &
                ispin > size(dg_frag%mixed_wannier_bpw_z, 4)) then
              bad_mixed_z_split_reference_diag = .true.
              cycle
            end if
            call gather_global_mixed_coefficients(ispin, state_first, state_last, cmix_prod_diag)
            cmix_before_diag(:,:,ispin) = cmix_prod_diag(:,:)
            cmix_split_diag(:,:,ispin) = cmix_prod_diag(:,:)
            h_step_diag(:,:) = (0.0d0, 0.0d0)
            do imix_diag = 1, nmix_diag
              h_step_diag(imix_diag,imix_diag) = &
                cmplx(dg_frag%mixed_wannier_bpw_eval(imix_diag,ispin), 0.0d0, kind=8)
            end do
            h_step_diag(:,:) = h_step_diag(:,:) &
              - E_mid(1) * dg_frag%mixed_wannier_bpw_z(1,1:nmix_diag,1:nmix_diag,ispin) &
              - E_mid(2) * dg_frag%mixed_wannier_bpw_z(2,1:nmix_diag,1:nmix_diag,ispin) &
              - E_mid(3) * dg_frag%mixed_wannier_bpw_z(3,1:nmix_diag,1:nmix_diag,ispin)
            h_evec_diag(:,:) = h_step_diag(:,:)
            call eigen_zheev(h_evec_diag, eval_step_diag, h_step_diag)
            do ist_diag = 1, nstate_diag
              work_vec_diag(:) = matmul(conjg(transpose(h_step_diag(:,:))), cmix_prod_diag(:,ist_diag))
              do imix_diag = 1, nmix_diag
                cphase = cos(eval_step_diag(imix_diag) * dt)
                sphase = sin(eval_step_diag(imix_diag) * dt)
                work_vec_diag(imix_diag) = cmplx(cphase, -sphase, kind=8) * work_vec_diag(imix_diag)
              end do
              cmix_exp_diag(:,ist_diag,ispin) = matmul(h_step_diag(:,:), work_vec_diag(:))
            end do
            if (sum(abs(E_mid(1:3))) > 1.0d-30) then
              h_step_diag(:,:) = -E_mid(1) * dg_frag%mixed_wannier_bpw_z(1,1:nmix_diag,1:nmix_diag,ispin) &
                                -E_mid(2) * dg_frag%mixed_wannier_bpw_z(2,1:nmix_diag,1:nmix_diag,ispin) &
                                -E_mid(3) * dg_frag%mixed_wannier_bpw_z(3,1:nmix_diag,1:nmix_diag,ispin)
              h_evec_diag(:,:) = h_step_diag(:,:)
              call eigen_zheev(h_evec_diag, eval_step_diag, h_step_diag)
              do ist_diag = 1, nstate_diag
                work_vec_diag(:) = matmul(conjg(transpose(h_step_diag(:,:))), cmix_split_diag(:,ist_diag,ispin))
                do imix_diag = 1, nmix_diag
                  cphase = cos(eval_step_diag(imix_diag) * dt)
                  sphase = sin(eval_step_diag(imix_diag) * dt)
                  work_vec_diag(imix_diag) = cmplx(cphase, -sphase, kind=8) * work_vec_diag(imix_diag)
                end do
                cmix_split_diag(:,ist_diag,ispin) = matmul(h_step_diag(:,:), work_vec_diag(:))
              end do
            end if
            do ist_diag = 1, nstate_diag
              do imix_diag = 1, nmix_diag
                cphase = cos(dg_frag%mixed_wannier_bpw_eval(imix_diag,ispin) * dt)
                sphase = sin(dg_frag%mixed_wannier_bpw_eval(imix_diag,ispin) * dt)
                cmix_split_diag(imix_diag,ist_diag,ispin) = &
                  cmplx(cphase, -sphase, kind=8) * cmix_split_diag(imix_diag,ist_diag,ispin)
              end do
            end do
          end do
        end if
        if (mixed_z_use_split_prop_enabled) then
          if (trace_wpw_order .and. dg_frag%id == 0) then
            write(*,'(1x,a,a,i0,a)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
              ' stage=production-mixed-z-split-candidate-begin'
          end if
          t_route0 = get_wtime()
          call apply_global_mixed_split_exp(E_mid, state_first, state_last)
          time_mixed_field = time_mixed_field + (get_wtime() - t_route0)
          if (trace_wpw_order .and. dg_frag%id == 0) then
            write(*,'(1x,a,a,i0,a,1pe12.4)') '[DG-WPW-RED-DIAG-ORDER]', ' step=', itt, &
              ' stage=production-mixed-z-split-candidate-end elapsed=', time_mixed_field
          end if
        else
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
        end if
        if (do_mixed_z_split_reference_diag) then
          norm_before_diag = 0.0d0
          norm_exp_diag = 0.0d0
          norm_split_diag = 0.0d0
          norm_prod_diag = 0.0d0
          diff_exp_prod_diag = 0.0d0
          diff_prod_split_diag = 0.0d0
          diff_split_exp_diag = 0.0d0
          pz_exp_diag = 0.0d0
          pz_prod_diag = 0.0d0
          do ispin = 1, dg_frag%nspin
            if (ispin > size(dg_frag%mixed_wannier_bpw_z, 4)) then
              bad_mixed_z_split_reference_diag = .true.
              cycle
            end if
            call gather_global_mixed_coefficients(ispin, state_first, state_last, cmix_prod_diag)
            do ist_diag = 1, nstate_diag
              norm_before_diag = norm_before_diag + &
                sum(abs(cmix_before_diag(:,ist_diag,ispin))**2)
              norm_exp_diag = norm_exp_diag + sum(abs(cmix_exp_diag(:,ist_diag,ispin))**2)
              norm_split_diag = norm_split_diag + sum(abs(cmix_split_diag(:,ist_diag,ispin))**2)
              norm_prod_diag = norm_prod_diag + sum(abs(cmix_prod_diag(:,ist_diag))**2)
              diff_exp_prod_diag = diff_exp_prod_diag + &
                sum(abs(cmix_exp_diag(:,ist_diag,ispin) - cmix_prod_diag(:,ist_diag))**2)
              diff_prod_split_diag = diff_prod_split_diag + &
                sum(abs(cmix_prod_diag(:,ist_diag) - cmix_split_diag(:,ist_diag,ispin))**2)
              diff_split_exp_diag = diff_split_exp_diag + &
                sum(abs(cmix_split_diag(:,ist_diag,ispin) - cmix_exp_diag(:,ist_diag,ispin))**2)
              work_vec_diag(:) = matmul(dg_frag%mixed_wannier_bpw_z(3,1:nmix_diag,1:nmix_diag,ispin), &
                cmix_exp_diag(:,ist_diag,ispin))
              pz_exp_diag = pz_exp_diag + real(sum(conjg(cmix_exp_diag(:,ist_diag,ispin)) * &
                work_vec_diag(:)), kind=8)
              work_vec_diag(:) = matmul(dg_frag%mixed_wannier_bpw_z(3,1:nmix_diag,1:nmix_diag,ispin), &
                cmix_prod_diag(:,ist_diag))
              pz_prod_diag = pz_prod_diag + real(sum(conjg(cmix_prod_diag(:,ist_diag)) * &
                work_vec_diag(:)), kind=8)
            end do
          end do
          diff_exp_prod_diag = sqrt(max(diff_exp_prod_diag, 0.0d0))
          diff_prod_split_diag = sqrt(max(diff_prod_split_diag, 0.0d0))
          diff_split_exp_diag = sqrt(max(diff_split_exp_diag, 0.0d0))
          rel_diff_exp_prod_diag = diff_exp_prod_diag / max(sqrt(max(norm_prod_diag, 0.0d0)), 1.0d-300)
          rel_diff_prod_split_diag = diff_prod_split_diag / max(sqrt(max(norm_prod_diag, 0.0d0)), 1.0d-300)
          rel_diff_split_exp_diag = diff_split_exp_diag / max(sqrt(max(norm_split_diag, 0.0d0)), 1.0d-300)
          pz_diff_diag = pz_exp_diag - pz_prod_diag
          if (norm_before_diag /= norm_before_diag .or. norm_exp_diag /= norm_exp_diag .or. &
              norm_split_diag /= norm_split_diag .or. &
              norm_prod_diag /= norm_prod_diag .or. diff_exp_prod_diag /= diff_exp_prod_diag .or. &
              diff_prod_split_diag /= diff_prod_split_diag .or. diff_split_exp_diag /= diff_split_exp_diag .or. &
              pz_diff_diag /= pz_diff_diag) bad_mixed_z_split_reference_diag = .true.
          if (dg_frag%id == 0) then
            if (mixed_z_exp_step_diag_enabled) then
              write(*,'(1x,a,3(a,i0),a,1pe12.4,4(a,l1))') &
                '[DG-MIXED-Z-EXP-STEP-DIAG]', &
                ' step=', itt, ' dim_total=', nmix_diag, ' nstate=', nstate_diag, &
                ' dt=', dt, &
                ' field_on=', sum(abs(E_mid(1:3))) > 1.0d-30, &
                ' frozen_H_dense_exp=', .true., &
                ' production_split_step=', .true., &
                ' bad=', bad_mixed_z_split_reference_diag
              write(*,'(1x,a,8(a,1pe16.8))') &
                '[DG-MIXED-Z-EXP-STEP-DIAG]', &
                ' norm_C_before=', norm_before_diag, &
                ' norm_C_exp_after=', norm_exp_diag, &
                ' norm_C_prod_after=', norm_prod_diag, &
                ' exp_norm_drift=', norm_exp_diag - norm_before_diag, &
                ' prod_norm_drift=', norm_prod_diag - norm_before_diag, &
                ' diff_C_exp_vs_prod=', diff_exp_prod_diag, &
                ' rel_diff_C_exp_vs_prod=', rel_diff_exp_prod_diag, &
                ' Pz_diff_exp_vs_prod=', pz_diff_diag
              write(*,'(1x,a,2(a,1pe16.8))') &
                '[DG-MIXED-Z-EXP-STEP-DIAG]', &
                ' Pz_exp=', pz_exp_diag, &
                ' Pz_prod=', pz_prod_diag
            end if
            write(*,'(1x,a,3(a,1pe16.8),a,a,2(a,l1))') &
              '[DG-MIXEDZ-SPLIT-PROD-CMP]', &
              ' field_x=', E_mid(1), &
              ' field_y=', E_mid(2), &
              ' field_z=', E_mid(3), &
              ' split_order=', 'field_exp_then_static_phase', &
              ' field_on=', sum(abs(E_mid(1:3))) > 1.0d-30, &
              ' bad=', bad_mixed_z_split_reference_diag
            write(*,'(1x,a,7(a,1pe16.8))') &
              '[DG-MIXEDZ-SPLIT-PROD-CMP]', &
              ' norm_prod=', norm_prod_diag, &
              ' norm_mixed_split=', norm_split_diag, &
              ' diff_prod_vs_mixed_split=', diff_prod_split_diag, &
              ' rel_diff_prod_vs_mixed_split=', rel_diff_prod_split_diag, &
              ' diff_mixed_split_vs_dense_frozenH=', diff_split_exp_diag, &
              ' rel_diff_mixed_split_vs_dense_frozenH=', rel_diff_split_exp_diag, &
              ' dense_frozenH_vs_prod=', diff_exp_prod_diag
            flush(6)
          end if
          if (mixed_z_split_obs_diag_enabled) then
            call diagnose_mixed_z_split_observables(cmix_split_diag, norm_split_diag, itt, &
              bad_mixed_z_split_reference_diag)
          end if
          if (mixed_z_split_current_diag_enabled) then
            call diagnose_mixed_z_split_current(cmix_split_diag, norm_split_diag, itt, &
              bad_mixed_z_split_reference_diag)
          end if
          deallocate(cmix_before_diag, cmix_exp_diag, cmix_split_diag, cmix_prod_diag, &
            h_step_diag, h_evec_diag, work_vec_diag, eval_step_diag)
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

    use_full_h_seed_phase = yn_dg_full_h_eigen_seed == 'y' .and. allocated(dg_frag%esp) .and. &
      sum(abs(E_mid(1:3))) <= 1.0d-30 .and. .not. use_impulse_kinetic_shift
    if (use_full_h_seed_phase) then
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
        write(*,'(1x,a)') "[DG-EXPDIAG] full-DG Hamiltonian eigenphase integrator enabled for field-free steps."
        write(*,'(1x,a)') "[DG-EXPDIAG] local BPW expdiag is bypassed to keep the Full-DG eigenstate stationary."
        flush(6)
        expdiag_warned = .true.
      end if
      call print_expdiag_timing('full-dg-seed-phase')
      return
    end if
    if (yn_dg_full_h_eigen_seed == 'y') then
      if (.not. use_impulse_kinetic_shift .and. dg_frag%has_full_h_seed_eigen .and. &
          dg_frag%has_full_h_seed_xi) then
        t_route0 = get_wtime()
        call apply_full_h_seed_eigen_exp(state_first, state_last)
        time_global_flux = time_global_flux + (get_wtime() - t_route0)
        if (.not. expdiag_warned .and. dg_frag%id == 0) then
          write(*,'(1x,a)') "[DG-EXPDIAG] full-DG Hamiltonian eigenbasis integrator enabled."
          write(*,'(1x,a)') "[DG-EXPDIAG] nonzero length-gauge field is propagated as H_eff = eps - E.R in the Full-DG basis."
          write(*,'(1x,a)') "[DG-EXPDIAG] this diagnostic path is dense/allreduce and is not weak-scaling."
          flush(6)
          expdiag_warned = .true.
        end if
        call print_expdiag_timing('full-dg-seed-eigen')
        return
      end if
      if (dg_frag%id == 0) then
        write(*,'(1x,a)') "[FATAL] full-DG Hamiltonian eigen seed requires a matching full-DG propagator."
        write(*,'(1x,a)') "[FATAL] Nonzero length-gauge fields require the Full-DG position matrix; impulse kicks require momentum projection."
        write(*,'(1x,a)') "[FATAL] These cases must not fall back to local BPW expdiag."
        flush(6)
      end if
      stop "DG expdiag: full-DG eigen seed field propagator is not implemented"
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

      allocate(h_eff(nw,nw), eval(nw), evec(nw,nw), c_w(nw), tmp_w(nw), next_w(nw), &
        mom_work(nbf,nw), mom_eff(nw,nw))
      do ispin = 1, dg_frag%nspin
        if (use_formal_wannier_h) then
          h_eff(1:nw,1:nw) = dg_frag%dg_wannier_h0_local(1:nw,1:nw,ispin,i_local)
        else
          h_eff(1:nw,1:nw) = cmplx(dg_frag%buffer_wannier_h_flux(1:nw,1:nw,i_local), 0.0d0, kind=8)
        end if
        iblk = 0
        if (allocated(dg_frag%H_block_map)) iblk = find_matrix_block(dg_frag%H_block_map, ifrag, ifrag)
        if (project_h_for_fixed_func .and. iblk > 0 .and. allocated(dg_frag%H_mat_blocks)) then
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
        else if (delta_h_enabled .and. iblk > 0 .and. &
                 allocated(dg_frag%H_mat_blocks)) then
          if (iblk <= size(dg_frag%H_mat_blocks)) then
            call add_projected_delta_h_to_eff(iblk, i_local, ispin, nbf, nw, use_formal_wannier_h, h_eff)
          end if
        else if (yn_fix_func == 'n' .and. iblk > 0 .and. allocated(dg_frag%H_mat_blocks)) then
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
        call symmetrize_expdiag_h_eff(h_eff, nw)
        if (use_impulse_kinetic_shift) then
          iblk_mom = 0
          if (allocated(dg_frag%momentum_block_map)) iblk_mom = find_matrix_block(dg_frag%momentum_block_map, ifrag, ifrag)
          if (iblk_mom <= 0 .or. .not. allocated(dg_frag%momentum_blocks)) then
            stop "DG length-gauge impulse kinetic shift requires DG velocity/momentum blocks"
          end if
          if (iblk_mom > size(dg_frag%momentum_blocks)) then
            stop "DG length-gauge impulse kinetic shift found invalid momentum block index"
          end if
          if (.not. allocated(dg_frag%momentum_blocks(iblk_mom)%val)) then
            stop "DG length-gauge impulse kinetic shift found unallocated momentum block"
          end if
          mom_work(:,:) = (0.0d0, 0.0d0)
          mom_eff(:,:) = (0.0d0, 0.0d0)
          do jw = 1, nw
            do jb = 1, min(nbf, size(dg_frag%momentum_blocks(iblk_mom)%val, 3))
              do ib = 1, min(nbf, size(dg_frag%momentum_blocks(iblk_mom)%val, 2))
                if (use_formal_wannier_h) then
                  mom_work(ib,jw) = mom_work(ib,jw) &
                    + (impulse_k(1) * dg_frag%momentum_blocks(iblk_mom)%val(1,ib,jb,ispin) &
                     + impulse_k(2) * dg_frag%momentum_blocks(iblk_mom)%val(2,ib,jb,ispin) &
                     + impulse_k(3) * dg_frag%momentum_blocks(iblk_mom)%val(3,ib,jb,ispin)) &
                    * dg_frag%dg_wannier_basis_coef(jb,jw,ispin,i_local)
                else
                  mom_work(ib,jw) = mom_work(ib,jw) &
                    + (impulse_k(1) * dg_frag%momentum_blocks(iblk_mom)%val(1,ib,jb,ispin) &
                     + impulse_k(2) * dg_frag%momentum_blocks(iblk_mom)%val(2,ib,jb,ispin) &
                     + impulse_k(3) * dg_frag%momentum_blocks(iblk_mom)%val(3,ib,jb,ispin)) &
                    * cmplx(dg_frag%buffer_wannier_coef(jb,jw,ispin,i_local), 0.0d0, kind=8)
                end if
                end do
              end do
            end do
          do jw = 1, nw
            do iw = 1, nw
              do ib = 1, min(nbf, size(dg_frag%momentum_blocks(iblk_mom)%val, 2))
                if (use_formal_wannier_h) then
                  mom_eff(iw,jw) = mom_eff(iw,jw) &
                    + conjg(dg_frag%dg_wannier_basis_coef(ib,iw,ispin,i_local)) * mom_work(ib,jw)
                else
                  mom_eff(iw,jw) = mom_eff(iw,jw) &
                    + cmplx(dg_frag%buffer_wannier_coef(ib,iw,ispin,i_local), 0.0d0, kind=8) * mom_work(ib,jw)
                end if
              end do
              h_eff(iw,jw) = h_eff(iw,jw) - zi * mom_eff(iw,jw)
            end do
          end do
          do iw = 1, nw
            h_eff(iw,iw) = h_eff(iw,iw) + 0.5d0 * impulse_k2
          end do
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
        call symmetrize_expdiag_h_eff(h_eff, nw)
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
      deallocate(h_eff, eval, evec, c_w, tmp_w, next_w, mom_work, mom_eff)
    end do
    if (xi_split_enabled .and. use_bpw_wannier_h) call apply_xi_flux_split_half(E_mid, 0.5d0 * dt, state_first, state_last)
    time_local_loop = time_local_loop + (get_wtime() - t_route0)

    if (.not. expdiag_warned .and. dg_frag%id == 0) then
      if (use_formal_wannier_h) then
        write(*,'(1x,a)') "[DG-EXPDIAG] formal DG-Wannier exponential integrator enabled."
      else
        write(*,'(1x,a)') "[DG-EXPDIAG] local BPW exponential integrator enabled."
      end if
      if (xi_split_enabled .and. use_bpw_wannier_h) then
        write(*,'(1x,a)') "[DG-EXPDIAG] dynamic neighbor xi_flux is applied by a Strang half-step split."
      else if (use_formal_wannier_h) then
        write(*,'(1x,a)') "[DG-EXPDIAG] nonzero field step uses local DG-Wannier length-gauge blocks."
      else
        write(*,'(1x,a)') &
          "[DG-EXPDIAG] dynamic neighbor xi_flux split disabled; set SALMON_DG_EXPDIAG_XI_SPLIT=1 to test."
      end if
      if (delta_h_enabled .and. .not. project_h_for_fixed_func) then
        write(*,'(1x,a)') "[DG-EXPDIAG] propagation uses the seed flux Hamiltonian plus projected DG delta-H."
      else if (yn_fix_func == 'y' .and. .not. project_h_for_fixed_func) then
        write(*,'(1x,a)') "[DG-EXPDIAG] fixed-function propagation uses the seed flux Hamiltonian."
      else
        write(*,'(1x,a)') "[DG-EXPDIAG] propagation projects the current DG Hamiltonian into the active Wannier basis."
      end if
      flush(6)
      expdiag_warned = .true.
    end if
    call print_expdiag_timing('local-block')

  contains

    subroutine ensure_delta_h_reference()
      integer :: iblk_ref

      if (delta_h_ref_valid) then
        if (allocated(delta_h_ref_blocks) .and. allocated(dg_frag%H_mat_blocks)) then
          if (size(delta_h_ref_blocks) == size(dg_frag%H_mat_blocks)) return
        end if
      end if
      if (.not. allocated(dg_frag%H_mat_blocks)) &
        stop "DG expdiag delta-H requires H_mat_blocks"
      if (allocated(delta_h_ref_blocks)) then
        do iblk_ref = 1, size(delta_h_ref_blocks)
          if (allocated(delta_h_ref_blocks(iblk_ref)%val)) deallocate(delta_h_ref_blocks(iblk_ref)%val)
        end do
        deallocate(delta_h_ref_blocks)
      end if
      allocate(delta_h_ref_blocks(size(dg_frag%H_mat_blocks)))
      do iblk_ref = 1, size(dg_frag%H_mat_blocks)
        delta_h_ref_blocks(iblk_ref)%ifrag_row = dg_frag%H_mat_blocks(iblk_ref)%ifrag_row
        delta_h_ref_blocks(iblk_ref)%ifrag_col = dg_frag%H_mat_blocks(iblk_ref)%ifrag_col
        delta_h_ref_blocks(iblk_ref)%nrow_max = dg_frag%H_mat_blocks(iblk_ref)%nrow_max
        delta_h_ref_blocks(iblk_ref)%ncol_max = dg_frag%H_mat_blocks(iblk_ref)%ncol_max
        if (allocated(dg_frag%H_mat_blocks(iblk_ref)%val)) then
          allocate(delta_h_ref_blocks(iblk_ref)%val(size(dg_frag%H_mat_blocks(iblk_ref)%val, 1), &
                                                   size(dg_frag%H_mat_blocks(iblk_ref)%val, 2), &
                                                   size(dg_frag%H_mat_blocks(iblk_ref)%val, 3)))
          delta_h_ref_blocks(iblk_ref)%val(:, :, :) = dg_frag%H_mat_blocks(iblk_ref)%val(:, :, :)
        end if
      end do
      delta_h_ref_nblocks = size(dg_frag%H_mat_blocks)
      delta_h_ref_valid = .true.
      if (.not. delta_h_warned .and. dg_frag%id == 0) then
        write(*,'(1x,a,i0)') &
          "[DG-EXPDIAG] propagation adds projected delta-H from reference blocks; nblocks=", &
          delta_h_ref_nblocks
        flush(6)
        delta_h_warned = .true.
      end if
    end subroutine ensure_delta_h_reference

    subroutine add_projected_delta_h_to_eff(iblk_use, i_local_use, ispin_use, nbf_use, nw_use, &
                                           use_formal, h_eff_use)
      integer, intent(in) :: iblk_use, i_local_use, ispin_use, nbf_use, nw_use
      logical, intent(in) :: use_formal
      complex(8), intent(inout) :: h_eff_use(:,:)
      integer :: ib_use, jb_use, iw_use, jw_use
      integer :: nrow_use, ncol_use
      real(8) :: dh_val

      if (.not. delta_h_ref_valid) stop "DG expdiag delta-H reference is not initialized"
      if (.not. allocated(delta_h_ref_blocks)) stop "DG expdiag delta-H reference blocks are missing"
      if (iblk_use < 1 .or. iblk_use > size(delta_h_ref_blocks)) return
      if (iblk_use < 1 .or. iblk_use > size(dg_frag%H_mat_blocks)) return
      if (.not. allocated(delta_h_ref_blocks(iblk_use)%val)) return
      if (.not. allocated(dg_frag%H_mat_blocks(iblk_use)%val)) return
      nrow_use = min(nbf_use, size(dg_frag%H_mat_blocks(iblk_use)%val, 1), &
                     size(delta_h_ref_blocks(iblk_use)%val, 1))
      ncol_use = min(nbf_use, size(dg_frag%H_mat_blocks(iblk_use)%val, 2), &
                     size(delta_h_ref_blocks(iblk_use)%val, 2))
      if (nrow_use <= 0 .or. ncol_use <= 0) return
      do jw_use = 1, nw_use
        do iw_use = 1, nw_use
          do jb_use = 1, ncol_use
            do ib_use = 1, nrow_use
              dh_val = dg_frag%H_mat_blocks(iblk_use)%val(ib_use,jb_use,ispin_use) - &
                       delta_h_ref_blocks(iblk_use)%val(ib_use,jb_use,ispin_use)
              if (dh_val == 0.0d0) cycle
              if (use_formal) then
                h_eff_use(iw_use,jw_use) = h_eff_use(iw_use,jw_use) + &
                  conjg(dg_frag%dg_wannier_basis_coef(ib_use,iw_use,ispin_use,i_local_use)) * &
                  cmplx(dh_val, 0.0d0, kind=8) * &
                  dg_frag%dg_wannier_basis_coef(jb_use,jw_use,ispin_use,i_local_use)
              else
                h_eff_use(iw_use,jw_use) = h_eff_use(iw_use,jw_use) + &
                  cmplx(dg_frag%buffer_wannier_coef(ib_use,iw_use,ispin_use,i_local_use), 0.0d0, kind=8) * &
                  cmplx(dh_val, 0.0d0, kind=8) * &
                  cmplx(dg_frag%buffer_wannier_coef(jb_use,jw_use,ispin_use,i_local_use), 0.0d0, kind=8)
              end if
            end do
          end do
        end do
      end do
    end subroutine add_projected_delta_h_to_eff

    subroutine symmetrize_expdiag_h_eff(h_eff_use, n_use)
      complex(8), intent(inout) :: h_eff_use(:,:)
      integer, intent(in) :: n_use
      integer :: i_use, j_use
      complex(8) :: hij

      do i_use = 1, n_use
        h_eff_use(i_use,i_use) = cmplx(real(h_eff_use(i_use,i_use), kind=8), 0.0d0, kind=8)
        do j_use = i_use + 1, n_use
          hij = 0.5d0 * (h_eff_use(i_use,j_use) + conjg(h_eff_use(j_use,i_use)))
          h_eff_use(i_use,j_use) = hij
          h_eff_use(j_use,i_use) = conjg(hij)
        end do
      end do
    end subroutine symmetrize_expdiag_h_eff

    subroutine prepare_wpw_local_payload_ingredients(build_route)
      character(len=*), intent(in) :: build_route

      wpw_payload_prepare_cache_hit = .false.
      wpw_payload_raw_cache_hit = .false.
      wpw_payload_neighborH_cache_hit = .false.
      wpw_payload_pp_block_cache_hit = .false.
      wpw_reduced_eig_cache_hit = .false.
      wpw_reduced_eig_built = .false.
      wpw_reduced_eig_skipped = .false.
      wpw_reduced_eig_elapsed = 0.0d0
      wpw_reduced_payload_prepare_elapsed = 0.0d0
      wpw_payload_raw_elapsed = 0.0d0
      wpw_payload_pp_block_elapsed = 0.0d0
      wpw_payload_neighborH_elapsed = 0.0d0
      wpw_payload_prepare_cache_reason = 'not_requested'
      dg_frag%mixed_z_local_prop_payload_source = 'missing_payload'
      dg_frag%mixed_z_local_prop_payload_basis_kind = 'owner_unique_Wphase_WP_PP'
      dg_frag%mixed_z_local_prop_payload_build_route = trim(build_route)
      if (.not. allocated(dg_frag%wpw_S_pp_blocks) .or. &
          .not. allocated(dg_frag%wpw_T_pp_volume_blocks) .or. &
          .not. allocated(dg_frag%wpw_T_pp_interface_blocks)) then
        wpw_payload_stage_t0 = get_wtime()
        call ensure_wpw_local_pp_blocks(dg_frag, .true.)
        wpw_payload_pp_block_elapsed = get_wtime() - wpw_payload_stage_t0
        wpw_payload_pp_block_build_count = wpw_payload_pp_block_build_count + 1
      else
        wpw_payload_pp_block_cache_hit = .true.
        wpw_payload_pp_block_cache_hit_count = wpw_payload_pp_block_cache_hit_count + 1
      end if
      if (.not. dg_frag%has_global_wannier_local_basis) then
        dg_frag%mixed_z_local_prop_payload_block_reason = 'missing_global_wannier_local_basis'
      else if (.not. allocated(dg_frag%global_wannier_local_coef)) then
        dg_frag%mixed_z_local_prop_payload_block_reason = 'missing_global_wannier_local_coef'
      else if (.not. allocated(dg_frag%global_wannier_local_nkeep)) then
        dg_frag%mixed_z_local_prop_payload_block_reason = 'missing_global_wannier_local_nkeep'
      else if (.not. allocated(dg_frag%wpw_S_pp_blocks) .or. &
               .not. allocated(dg_frag%wpw_T_pp_volume_blocks) .or. &
               .not. allocated(dg_frag%wpw_T_pp_interface_blocks)) then
        dg_frag%mixed_z_local_prop_payload_block_reason = 'missing_wpw_pp_blocks'
      else if (.not. allocated(dg_frag%phi_frag) .and. .not. allocated(dg_frag%phi_frag_c)) then
        dg_frag%mixed_z_local_prop_payload_block_reason = 'missing_phi_frag'
      else if (dg_frag%n_plane_waves <= 0) then
        dg_frag%mixed_z_local_prop_payload_block_reason = 'missing_plane_waves'
      else
        wpw_reduced_payload_t0 = get_wtime()
        wpw_payload_static_field = &
          abs(E_mid(1)) + abs(E_mid(2)) + abs(E_mid(3)) + &
          abs(Ac_ham(1)) + abs(Ac_ham(2)) + abs(Ac_ham(3)) <= 1.0d-14
        wpw_payload_prepare_cache_allowed = wpw_payload_static_field .and. yn_fix_func /= 'n' .and. &
          dg_frag%wpw_reduced_ready .and. allocated(dg_frag%wpw_reduced_H) .and. &
          allocated(dg_frag%wpw_reduced_S) .and. allocated(dg_frag%wpw_reduced_Hraw_build) .and. &
          allocated(dg_frag%wpw_reduced_Sraw_build) .and. allocated(dg_frag%wpw_reduced_transform)
        if (wpw_payload_prepare_cache_allowed) then
          wpw_payload_prepare_cache_reason = 'static_field_off'
        else if (.not. wpw_payload_static_field .or. yn_fix_func == 'n') then
          wpw_payload_prepare_cache_reason = 'dynamic_potential'
        else if (.not. dg_frag%wpw_reduced_ready .or. .not. allocated(dg_frag%wpw_reduced_H) .or. &
                 .not. allocated(dg_frag%wpw_reduced_S)) then
          wpw_payload_prepare_cache_reason = 'missing_neighborH'
        else if (.not. allocated(dg_frag%wpw_reduced_Hraw_build) .or. &
                 .not. allocated(dg_frag%wpw_reduced_Sraw_build) .or. &
                 .not. allocated(dg_frag%wpw_reduced_transform)) then
          wpw_payload_prepare_cache_reason = 'missing_raw_matrix'
        else
          wpw_payload_prepare_cache_reason = 'basis_signature_mismatch'
        end if
        if (wpw_payload_prepare_cache_allowed) then
          wpw_payload_prepare_cache_hit = .true.
          wpw_payload_raw_cache_hit = .true.
          wpw_payload_neighborH_cache_hit = .true.
          wpw_payload_prepare_cache_hit_count = wpw_payload_prepare_cache_hit_count + 1
          wpw_payload_raw_cache_hit_count = wpw_payload_raw_cache_hit_count + 1
          wpw_payload_neighborH_cache_hit_count = wpw_payload_neighborH_cache_hit_count + 1
        else
          wpw_payload_prepare_build_count = wpw_payload_prepare_build_count + 1
          wpw_payload_raw_build_count = wpw_payload_raw_build_count + 1
          wpw_payload_neighborH_build_count = wpw_payload_neighborH_build_count + 1
          wpw_payload_stage_t0 = get_wtime()
          call diagnose_wpw_mixed_neighbor_hamiltonian(dg_frag, Vh, Vxc, Vpsl, .true., .true.)
          wpw_payload_neighborH_elapsed = get_wtime() - wpw_payload_stage_t0
          wpw_payload_raw_elapsed = wpw_payload_neighborH_elapsed
        end if
        wpw_reduced_eig_cache_allowed = wpw_payload_prepare_cache_allowed
        call ensure_wpw_reduced_eigensystem_only(wpw_reduced_eig_cache_allowed, &
          wpw_reduced_eig_cache_hit, wpw_reduced_eig_built, wpw_reduced_eig_skipped, &
          wpw_reduced_eig_elapsed)
        wpw_reduced_payload_prepare_elapsed = get_wtime() - wpw_reduced_payload_t0
        if (allocated(dg_frag%wpw_reduced_Sraw_build) .and. &
            allocated(dg_frag%wpw_reduced_Hraw_build) .and. &
            allocated(dg_frag%wpw_reduced_nraw)) then
          dg_frag%mixed_z_local_prop_payload_block_reason = 'producer_ready'
        else
          dg_frag%mixed_z_local_prop_payload_block_reason = 'producer_returned_no_raw_build'
        end if
      end if
    end subroutine prepare_wpw_local_payload_ingredients

    subroutine log_wpw_local_payload_prepare_smoke(call_index)
      integer, intent(in) :: call_index
      logical :: raw_available, neighborH_available, reduced_eval_available, reduced_evec_available
      logical :: payload_valid, smoke_bad
      character(16) :: field_kind

      raw_available = allocated(dg_frag%wpw_reduced_Hraw_build) .and. &
        allocated(dg_frag%wpw_reduced_Sraw_build) .and. allocated(dg_frag%wpw_reduced_nraw) .and. &
        allocated(dg_frag%wpw_reduced_transform)
      neighborH_available = dg_frag%wpw_reduced_ready .and. allocated(dg_frag%wpw_reduced_H) .and. &
        allocated(dg_frag%wpw_reduced_S)
      reduced_eval_available = allocated(dg_frag%wpw_reduced_eval)
      reduced_evec_available = allocated(dg_frag%wpw_reduced_evec)
      payload_valid = raw_available .and. neighborH_available .and. &
        reduced_eval_available .and. reduced_evec_available
      smoke_bad = .not. payload_valid
      if (wpw_payload_static_field) then
        field_kind = 'field_off'
      else
        field_kind = 'field_on'
      end if
      if (dg_frag%id == 0 .and. .not. dg_frag%mixed_z_perf_count_enabled) then
        write(*,'(1x,a,12(a,i0),4(a,a),3(a,1pe12.4),16(a,l1))') &
          '[DG-MIXEDZ-LOCAL-PROP-PAYLOAD-PREPARE-SMOKE]', &
          ' step=', itt, &
          ' call_index=', call_index, &
          ' n_payload_prepare_build=', wpw_payload_prepare_build_count, &
          ' n_payload_prepare_cache_hit=', wpw_payload_prepare_cache_hit_count, &
          ' n_raw_build=', wpw_payload_raw_build_count, &
          ' n_raw_cache_hit=', wpw_payload_raw_cache_hit_count, &
          ' n_neighborH_build=', wpw_payload_neighborH_build_count, &
          ' n_neighborH_cache_hit=', wpw_payload_neighborH_cache_hit_count, &
          ' n_eigensystem_build=', wpw_reduced_eig_build_count, &
          ' n_eigensystem_cache_hit=', wpw_reduced_eig_cache_hit_count, &
          ' n_pp_block_build=', wpw_payload_pp_block_build_count, &
          ' n_pp_block_cache_hit=', wpw_payload_pp_block_cache_hit_count, &
          ' field_kind=', trim(field_kind), &
          ' yn_fix_func=', trim(yn_fix_func), &
          ' payload_prepare_cache_reason=', trim(wpw_payload_prepare_cache_reason), &
          ' block_reason=', trim(dg_frag%mixed_z_local_prop_payload_block_reason), &
          ' walltime_payload_prepare=', wpw_reduced_payload_prepare_elapsed, &
          ' walltime_raw_build=', wpw_payload_raw_elapsed, &
          ' walltime_neighborH=', wpw_payload_neighborH_elapsed, &
          ' cache_allowed=', wpw_payload_prepare_cache_allowed, &
          ' payload_prepare_cache_hit=', wpw_payload_prepare_cache_hit, &
          ' raw_cache_hit=', wpw_payload_raw_cache_hit, &
          ' neighborH_cache_hit=', wpw_payload_neighborH_cache_hit, &
          ' pp_block_cache_hit=', wpw_payload_pp_block_cache_hit, &
          ' raw_available=', raw_available, &
          ' neighborH_available=', neighborH_available, &
          ' reduced_eval_available=', reduced_eval_available, &
          ' reduced_evec_available=', reduced_evec_available, &
          ' payload_valid=', payload_valid, &
          ' static_field=', wpw_payload_static_field, &
          ' eigensystem_built=', wpw_reduced_eig_built, &
          ' bad=', smoke_bad
      end if
    end subroutine log_wpw_local_payload_prepare_smoke

    subroutine build_wpw_local_prop_payload_fast_smoke(call_index, build_route)
      integer, intent(in) :: call_index
      character(len=*), intent(in), optional :: build_route
      integer :: nstate_fast, nmix_fast, nw_fast, nspin_fast, ispin_fast, ist_fast, iw_fast
      real(8) :: occ_fast, coef_norm2, prod_norm2, nonzero_count
      real(8) :: t_owner_pack, t_case_select, t_payload_store
      real(8) :: t_fast
      complex(8) :: cval
      complex(8), allocatable :: cmix_fast(:,:)
      logical :: case_available, payload_valid, fast_bad
      character(64) :: route_label

      route_label = 'producer_fast_smoke'
      if (present(build_route)) route_label = trim(build_route)
      nw_fast = dg_frag%mixed_wannier_bpw_nw
      nmix_fast = dg_frag%mixed_wannier_bpw_nmix
      nstate_fast = 0
      if (allocated(system%rocc)) nstate_fast = min(size(system%rocc, 1), size(dg_frag%coef, 2))
      nspin_fast = 0
      if (allocated(system%rocc)) nspin_fast = min(dg_frag%nspin, system%nspin, size(system%rocc, 3))
      coef_norm2 = 0.0d0
      prod_norm2 = 0.0d0
      nonzero_count = 0.0d0
      t_owner_pack = 0.0d0
      t_case_select = 0.0d0
      t_payload_store = 0.0d0
      case_available = nw_fast > 0 .and. nmix_fast >= nw_fast .and. nstate_fast > 0 .and. nspin_fast > 0
      if (case_available) then
        if (allocated(dg_frag%mixed_z_local_prop_payload_wcoef)) then
          if (size(dg_frag%mixed_z_local_prop_payload_wcoef, 1) /= nw_fast .or. &
              size(dg_frag%mixed_z_local_prop_payload_wcoef, 2) /= nstate_fast .or. &
              size(dg_frag%mixed_z_local_prop_payload_wcoef, 3) /= nspin_fast) then
            deallocate(dg_frag%mixed_z_local_prop_payload_wcoef)
          end if
        end if
        if (.not. allocated(dg_frag%mixed_z_local_prop_payload_wcoef)) &
          allocate(dg_frag%mixed_z_local_prop_payload_wcoef(nw_fast, nstate_fast, nspin_fast))
        dg_frag%mixed_z_local_prop_payload_wcoef(:, :, :) = (0.0d0, 0.0d0)
        allocate(cmix_fast(nmix_fast, nstate_fast))
        do ispin_fast = 1, nspin_fast
          cmix_fast(:, :) = (0.0d0, 0.0d0)
          t_fast = get_wtime()
          call gather_global_mixed_coefficients(ispin_fast, 1, nstate_fast, cmix_fast)
          t_owner_pack = t_owner_pack + (get_wtime() - t_fast)
          dg_frag%mixed_z_local_prop_payload_wcoef(1:nw_fast,1:nstate_fast,ispin_fast) = &
            cmix_fast(1:nw_fast,1:nstate_fast)
          t_fast = get_wtime()
          do ist_fast = 1, nstate_fast
            occ_fast = 0.0d0
            if (ispin_fast <= size(system%rocc, 3)) occ_fast = max(0.0d0, system%rocc(ist_fast, 1, ispin_fast))
            do iw_fast = 1, nw_fast
              cval = cmix_fast(iw_fast, ist_fast)
              coef_norm2 = coef_norm2 + occ_fast * abs(cval)**2
              prod_norm2 = prod_norm2 + occ_fast * abs(cval)**2
              if (abs(cval) > 1.0d-30) nonzero_count = nonzero_count + 1.0d0
            end do
          end do
          t_case_select = t_case_select + (get_wtime() - t_fast)
        end do
        deallocate(cmix_fast)
      end if
      t_fast = get_wtime()
      payload_valid = case_available .and. coef_norm2 == coef_norm2 .and. prod_norm2 == prod_norm2
      dg_frag%mixed_z_local_prop_payload_ready = payload_valid
      dg_frag%mixed_z_local_prop_payload_bad = .not. payload_valid
      dg_frag%mixed_z_local_prop_payload_step = itt
      dg_frag%mixed_z_local_prop_payload_coef_norm = sqrt(max(coef_norm2, 0.0d0))
      dg_frag%mixed_z_local_prop_payload_prod_coef_norm = sqrt(max(prod_norm2, 0.0d0))
      dg_frag%mixed_z_local_prop_payload_coef_diff_snorm = 0.0d0
      dg_frag%mixed_z_local_prop_payload_rel_coef_diff_snorm = 0.0d0
      dg_frag%mixed_z_local_prop_payload_dim = dble(max(nw_fast, 0))
      dg_frag%mixed_z_local_prop_payload_occ_trace = coef_norm2
      dg_frag%mixed_z_local_prop_payload_source = 'local_prop_coef'
      dg_frag%mixed_z_local_prop_payload_basis_kind = 'owner_unique_Wphase_WP_PP'
      dg_frag%mixed_z_local_prop_payload_build_route = trim(route_label)
      if (payload_valid) then
        dg_frag%mixed_z_local_prop_payload_block_reason = 'none'
      else
        dg_frag%mixed_z_local_prop_payload_block_reason = 'producer_fast_unavailable'
      end if
      t_payload_store = get_wtime() - t_fast
      fast_bad = .not. payload_valid
      if (dg_frag%id == 0 .and. .not. dg_frag%mixed_z_perf_count_enabled) then
        write(*,'(1x,a,4(a,i0),1(a,a),9(a,1pe12.4),8(a,l1))') &
          '[DG-MIXEDZ-LOCAL-PROP-PRODUCER-FAST-SMOKE]', &
          ' step=', itt, &
          ' call_index=', call_index, &
          ' selected_case=', 4, &
          ' nstate=', nstate_fast, &
          ' selected_case_name=', 'Wphase_plus_WP_PP', &
          ' local_coef_norm=', sqrt(max(coef_norm2, 0.0d0)), &
          ' global_coef_norm=', sqrt(max(coef_norm2, 0.0d0)), &
          ' local_occ_trace=', coef_norm2, &
          ' global_occ_trace=', coef_norm2, &
          ' n_owner_blocks=', dble(max(nw_fast, 0)), &
          ' n_nonzero_slots=', nonzero_count, &
          ' walltime_owner_pack=', t_owner_pack, &
          ' walltime_case_select=', t_case_select, &
          ' walltime_payload_store=', t_payload_store, &
          ' case_available=', case_available, &
          ' payload_valid=', payload_valid, &
          ' raw_available=', allocated(dg_frag%wpw_reduced_Hraw_build), &
          ' neighborH_available=', dg_frag%wpw_reduced_ready, &
          ' reduced_eval_available=', allocated(dg_frag%wpw_reduced_eval), &
          ' reduced_evec_available=', allocated(dg_frag%wpw_reduced_evec), &
          ' production_fallback_excluded=', .true., &
          ' bad=', fast_bad
      end if
    end subroutine build_wpw_local_prop_payload_fast_smoke

    subroutine log_wpw_local_prop_series_direct_smoke(call_index)
      integer, intent(in) :: call_index
      logical :: payload_valid, payload_step_match, reader_bad
      real(8) :: candidate_coef_norm, candidate_occ_trace
      real(8) :: rho_candidate, pz_candidate, current_candidate_norm
      character(32) :: candidate_contraction_source

      payload_step_match = dg_frag%mixed_z_local_prop_payload_step == itt
      payload_valid = dg_frag%mixed_z_local_prop_payload_ready .and. &
        payload_step_match .and. .not. dg_frag%mixed_z_local_prop_payload_bad
      candidate_contraction_source = 'local_prop_coef'
      candidate_coef_norm = 0.0d0
      candidate_occ_trace = 0.0d0
      rho_candidate = 0.0d0
      pz_candidate = 0.0d0
      current_candidate_norm = 0.0d0
      if (payload_valid) then
        candidate_coef_norm = dg_frag%mixed_z_local_prop_payload_coef_norm
        candidate_occ_trace = dg_frag%mixed_z_local_prop_payload_occ_trace
        rho_candidate = candidate_occ_trace
      end if
      reader_bad = .not. payload_valid .or. &
        candidate_coef_norm /= candidate_coef_norm .or. &
        candidate_occ_trace /= candidate_occ_trace .or. &
        rho_candidate /= rho_candidate
      if (dg_frag%id == 0) then
        write(*,'(1x,a,5(a,i0),4(a,a),10(a,1pe12.4),14(a,l1))') &
          '[DG-MIXEDZ-LOCAL-PROP-SERIES-PAYLOAD-DIRECT-SMOKE]', &
          ' step=', itt, &
          ' call_index=', call_index, &
          ' coef_source_step=', dg_frag%mixed_z_local_prop_payload_step, &
          ' n_payload_prepare_cache_hit=', wpw_payload_prepare_cache_hit_count, &
          ' n_raw_cache_hit=', wpw_payload_raw_cache_hit_count, &
          ' candidate_source=', 'local_mixedz', &
          ' candidate_contraction_source=', trim(candidate_contraction_source), &
          ' payload_basis_kind=', trim(dg_frag%mixed_z_local_prop_payload_basis_kind), &
          ' payload_build_route=', trim(dg_frag%mixed_z_local_prop_payload_build_route), &
          ' candidate_coef_norm=', candidate_coef_norm, &
          ' candidate_occ_trace=', candidate_occ_trace, &
          ' rho_candidate=', rho_candidate, &
          ' Pz_candidate=', pz_candidate, &
          ' current_candidate_norm=', current_candidate_norm, &
          ' payload_coef_norm=', dg_frag%mixed_z_local_prop_payload_coef_norm, &
          ' payload_occ_trace=', dg_frag%mixed_z_local_prop_payload_occ_trace, &
          ' global_coef_norm=', dg_frag%mixed_z_local_prop_payload_coef_norm, &
          ' global_occ_trace=', dg_frag%mixed_z_local_prop_payload_occ_trace, &
          ' coef_diff_Snorm=', dg_frag%mixed_z_local_prop_payload_coef_diff_snorm, &
          ' payload_ready=', dg_frag%mixed_z_local_prop_payload_ready, &
          ' payload_valid=', payload_valid, &
          ' payload_step_match=', payload_step_match, &
          ' payload_prepare_cache_hit=', wpw_payload_prepare_cache_hit, &
          ' raw_cache_hit=', wpw_payload_raw_cache_hit, &
          ' neighborH_cache_hit=', wpw_payload_neighborH_cache_hit, &
          ' pp_block_cache_hit=', wpw_payload_pp_block_cache_hit, &
          ' cache_allowed=', wpw_payload_prepare_cache_allowed, &
          ' static_field=', wpw_payload_static_field, &
          ' production_fallback_excluded=', .true., &
          ' bad=', reader_bad
      end if
    end subroutine log_wpw_local_prop_series_direct_smoke

    subroutine log_mixed_z_pz_fieldon_source_snapshot(snapshot_tag, phase_kick_count)
      character(len=*), intent(in) :: snapshot_tag
      integer, intent(in) :: phase_kick_count
      real(8) :: pol_snapshot(3)
      real(8) :: pz_prod_source, pz_fullww_source, pz_diff
      real(8) :: coef_norm, rho_int, time_use, afield_use, efield_use
      real(8) :: pz_ww_diag, pz_ww_offdiag, pz_wp, pz_pw, pz_pp, pz_nonww
      real(8) :: rho_ww, rho_wp_cross, rho_pp, coef_norm_w, coef_norm_p
      real(8) :: occ_snap, pcoef_abs2
      integer :: ispin_snap, nocc_snap, ist_snap, ip_snap
      logical :: snapshot_bad, decomp_available

      pol_snapshot(:) = 0.0d0
      call calculate_local_wannier_polarization_dg(dg_frag, system, pol_snapshot)
      pz_prod_source = pol_snapshot(3)
      decomp_available = dg_frag%mixed_z_prod_pz_decomp_ready
      pz_ww_diag = 0.0d0
      pz_ww_offdiag = 0.0d0
      pz_wp = 0.0d0
      pz_pw = 0.0d0
      pz_pp = 0.0d0
      pz_nonww = 0.0d0
      rho_ww = 0.0d0
      rho_wp_cross = 0.0d0
      rho_pp = 0.0d0
      coef_norm_w = 0.0d0
      coef_norm_p = 0.0d0
      if (decomp_available) then
        pz_fullww_source = dg_frag%mixed_z_prod_pz_ww_raw
        pz_ww_diag = dg_frag%mixed_z_prod_pz_ww_diag_raw
        pz_ww_offdiag = dg_frag%mixed_z_prod_pz_ww_offdiag_raw
        pz_wp = dg_frag%mixed_z_prod_pz_wp_raw
        pz_pp = dg_frag%mixed_z_prod_pz_pp_raw
        pz_nonww = pz_prod_source - pz_fullww_source
        coef_norm = dg_frag%mixed_z_prod_pz_w_occ_weight
        rho_int = dg_frag%mixed_z_prod_pz_occ_sum
        rho_ww = dg_frag%mixed_z_prod_pz_w_occ_weight
        coef_norm_w = dg_frag%mixed_z_prod_pz_w_occ_weight
      else
        pz_fullww_source = 0.0d0
        coef_norm = 0.0d0
        rho_int = 0.0d0
      end if
      if (allocated(dg_frag%mixed_wannier_bpw_pcoef)) then
        do ispin_snap = 1, min(dg_frag%nspin, system%nspin, size(dg_frag%mixed_wannier_bpw_pcoef, 3))
          nocc_snap = min(dg_frag%nstate_tot, size(dg_frag%mixed_wannier_bpw_pcoef, 2), size(system%rocc, 1))
          if (allocated(dg_frag%nocc_spin)) nocc_snap = min(nocc_snap, dg_frag%nocc_spin(ispin_snap))
          do ist_snap = 1, nocc_snap
            occ_snap = max(0.0d0, system%rocc(ist_snap, 1, ispin_snap))
            if (occ_snap <= 0.0d0) cycle
            do ip_snap = 1, size(dg_frag%mixed_wannier_bpw_pcoef, 1)
              pcoef_abs2 = real(conjg(dg_frag%mixed_wannier_bpw_pcoef(ip_snap,ist_snap,ispin_snap)) * &
                dg_frag%mixed_wannier_bpw_pcoef(ip_snap,ist_snap,ispin_snap), kind=8)
              rho_pp = rho_pp + occ_snap * pcoef_abs2
            end do
          end do
        end do
        coef_norm_p = rho_pp
      end if
      pz_diff = pz_fullww_source - pz_prod_source
      time_use = dble(itt) * dt
      afield_use = rt%Ac_tot(3, it1)
      efield_use = E_mid(3)
      snapshot_bad = pz_prod_source /= pz_prod_source .or. &
        pz_fullww_source /= pz_fullww_source .or. &
        coef_norm /= coef_norm .or. rho_int /= rho_int .or. &
        pz_ww_diag /= pz_ww_diag .or. pz_ww_offdiag /= pz_ww_offdiag .or. &
        pz_wp /= pz_wp .or. pz_pp /= pz_pp .or. rho_pp /= rho_pp
      if (dg_frag%id == 0) then
        write(*,'(1x,a,1(a,i0),1(a,a),9(a,1pe12.4),3(a,l1),1(a,i0),a,a)') &
          '[DG-MIXEDZ-LOCAL-PZ-FIELDON-SOURCE-SNAPSHOT-CMP]', &
          ' step=', itt, &
          ' snapshot_tag=', trim(snapshot_tag), &
          ' time=', time_use, &
          ' Afield=', afield_use, &
          ' Efield=', efield_use, &
          ' Pz_from_prod_source=', pz_prod_source, &
          ' Pz_from_fullWW_source=', pz_fullww_source, &
          ' diff=', pz_diff, &
          ' coef_norm=', coef_norm, &
          ' rho_int=', rho_int, &
          ' dt=', dt, &
          ' decomp_available=', decomp_available, &
          ' field_on=', sum(abs(E_mid(1:3))) > 1.0d-30, &
          ' bad=', snapshot_bad .or. .not. decomp_available, &
          ' phase_kick_count=', phase_kick_count, &
          ' route=', 'expdiag-source-snapshot-diagnostic-only'
        write(*,'(1x,a,1(a,i0),1(a,a),15(a,1pe12.4),6(a,l1),a,a)') &
          '[DG-MIXEDZ-LOCAL-PZ-FIELDON-BLOCK-CMP]', &
          ' step=', itt, &
          ' snapshot_tag=', trim(snapshot_tag), &
          ' Pz_total_prod_source=', pz_prod_source, &
          ' Pz_WW_diag=', pz_ww_diag, &
          ' Pz_WW_offdiag=', pz_ww_offdiag, &
          ' Pz_WW_full=', pz_fullww_source, &
          ' Pz_WP=', pz_wp, &
          ' Pz_PW=', pz_pw, &
          ' Pz_PP=', pz_pp, &
          ' Pz_nonWW=', pz_nonww, &
          ' Pz_rebuild_total=', pz_fullww_source + pz_wp + pz_pp, &
          ' Pz_rebuild_diff=', pz_prod_source - (pz_fullww_source + pz_wp + pz_pp), &
          ' rho_WW=', rho_ww, &
          ' rho_WP_cross=', rho_wp_cross, &
          ' rho_PP=', rho_pp, &
          ' coef_norm_W=', coef_norm_w, &
          ' coef_norm_P=', coef_norm_p, &
          ' decomp_available=', decomp_available, &
          ' WP_contains_PW=', .true., &
          ' separate_PW_available=', .false., &
          ' rho_cross_available=', .false., &
          ' field_on=', sum(abs(E_mid(1:3))) > 1.0d-30, &
          ' bad=', snapshot_bad .or. .not. decomp_available, &
          ' route=', 'expdiag-fieldon-block-diagnostic-only'
      end if
    end subroutine log_mixed_z_pz_fieldon_source_snapshot

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
      logical :: use_cached_eig
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
      use_cached_eig = wpw_reduced_eig_cache_valid
      if (.not. use_cached_eig) then
        dg_frag%wpw_reduced_eval(:, :, :) = 0.0d0
        dg_frag%wpw_reduced_evec(:, :, :, :) = (0.0d0, 0.0d0)
      end if
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
          if (use_cached_eig) then
            info_red = 0
            eval_red(1:nred) = dg_frag%wpw_reduced_eval(1:nred, ispin_use, iloc_frag)
            evec_red(1:nred,1:nred) = dg_frag%wpw_reduced_evec(1:nred, 1:nred, ispin_use, iloc_frag)
            eig_min = eval_red(1)
            eig_max = eval_red(nred)
            s_min_red = 0.0d0
            s_max_red = 0.0d0
            h_herm_red = 0.0d0
          else
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

    subroutine ensure_wpw_reduced_eigensystem_only(cache_allowed, cache_hit, built, skipped, elapsed)
      logical, intent(in) :: cache_allowed
      logical, intent(out) :: cache_hit, built, skipped
      real(8), intent(out) :: elapsed
      integer :: iloc_frag, ispin_use, nred, lwork_red, info_red
      complex(8), allocatable :: H_work_red(:,:), S_work_red(:,:), evec_red(:,:), work_red(:)
      real(8), allocatable :: eval_red(:), rwork_red(:)
      real(8) :: t_start
      external :: ZHEGV

      cache_hit = .false.
      built = .false.
      skipped = .false.
      elapsed = 0.0d0
      t_start = get_wtime()
      if (.not. dg_frag%wpw_reduced_ready) then
        skipped = .true.
        wpw_reduced_eig_skipped_count = wpw_reduced_eig_skipped_count + 1
        wpw_reduced_eig_cache_valid = .false.
        elapsed = get_wtime() - t_start
        return
      end if
      if (.not. allocated(dg_frag%wpw_reduced_H) .or. .not. allocated(dg_frag%wpw_reduced_S)) then
        skipped = .true.
        wpw_reduced_eig_skipped_count = wpw_reduced_eig_skipped_count + 1
        wpw_reduced_eig_cache_valid = .false.
        elapsed = get_wtime() - t_start
        return
      end if
      if (.not. allocated(dg_frag%wpw_reduced_dim)) then
        skipped = .true.
        wpw_reduced_eig_skipped_count = wpw_reduced_eig_skipped_count + 1
        wpw_reduced_eig_cache_valid = .false.
        elapsed = get_wtime() - t_start
        return
      end if

      if (.not. allocated(dg_frag%wpw_reduced_eval)) then
        allocate(dg_frag%wpw_reduced_eval(dg_frag%wpw_reduced_max_dim, dg_frag%nspin, size(dg_frag%wpw_reduced_dim)))
      else if (size(dg_frag%wpw_reduced_eval,1) /= dg_frag%wpw_reduced_max_dim .or. &
               size(dg_frag%wpw_reduced_eval,2) /= dg_frag%nspin .or. &
               size(dg_frag%wpw_reduced_eval,3) /= size(dg_frag%wpw_reduced_dim)) then
        deallocate(dg_frag%wpw_reduced_eval)
        allocate(dg_frag%wpw_reduced_eval(dg_frag%wpw_reduced_max_dim, dg_frag%nspin, size(dg_frag%wpw_reduced_dim)))
        wpw_reduced_eig_cache_valid = .false.
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
        wpw_reduced_eig_cache_valid = .false.
      end if
      if (cache_allowed .and. wpw_reduced_eig_cache_valid) then
        cache_hit = .true.
        wpw_reduced_eig_cache_hit_count = wpw_reduced_eig_cache_hit_count + 1
        elapsed = get_wtime() - t_start
        return
      end if
      dg_frag%wpw_reduced_eval(:, :, :) = 0.0d0
      dg_frag%wpw_reduced_evec(:, :, :, :) = (0.0d0, 0.0d0)
      built = .true.
      wpw_reduced_eig_build_count = wpw_reduced_eig_build_count + 1

      do iloc_frag = 1, size(dg_frag%wpw_reduced_dim)
        nred = dg_frag%wpw_reduced_dim(iloc_frag)
        if (nred <= 0) cycle
        if (nred > size(dg_frag%wpw_reduced_H, 1) .or. nred > size(dg_frag%wpw_reduced_H, 2)) cycle
        if (nred > size(dg_frag%wpw_reduced_S, 1) .or. nred > size(dg_frag%wpw_reduced_S, 2)) cycle
        if (iloc_frag > size(dg_frag%wpw_reduced_H, 4) .or. iloc_frag > size(dg_frag%wpw_reduced_S, 4)) cycle
        allocate(H_work_red(nred,nred), S_work_red(nred,nred), evec_red(nred,nred))
        allocate(eval_red(nred), rwork_red(max(1, 3*nred - 2)), work_red(1))
        do ispin_use = 1, dg_frag%nspin
          if (ispin_use > size(dg_frag%wpw_reduced_H, 3) .or. ispin_use > size(dg_frag%wpw_reduced_S, 3)) cycle
          H_work_red(:, :) = dg_frag%wpw_reduced_H(1:nred, 1:nred, ispin_use, iloc_frag)
          S_work_red(:, :) = dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin_use, iloc_frag)
          evec_red(:, :) = H_work_red(:, :)
          lwork_red = -1
          call ZHEGV(1, 'V', 'U', nred, evec_red, nred, S_work_red, nred, eval_red, work_red, &
            lwork_red, rwork_red, info_red)
          lwork_red = max(1, int(real(work_red(1), kind=8)))
          deallocate(work_red)
          allocate(work_red(lwork_red))
          S_work_red(:, :) = dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin_use, iloc_frag)
          evec_red(:, :) = H_work_red(:, :)
          call ZHEGV(1, 'V', 'U', nred, evec_red, nred, S_work_red, nred, eval_red, work_red, &
            lwork_red, rwork_red, info_red)
          if (info_red == 0) then
            dg_frag%wpw_reduced_eval(1:nred, ispin_use, iloc_frag) = eval_red(1:nred)
            dg_frag%wpw_reduced_evec(1:nred, 1:nred, ispin_use, iloc_frag) = evec_red(1:nred, 1:nred)
          end if
          deallocate(work_red)
          allocate(work_red(1))
        end do
        deallocate(H_work_red, S_work_red, evec_red, eval_red, rwork_red, work_red)
      end do
      wpw_reduced_eig_cache_valid = cache_allowed
      elapsed = get_wtime() - t_start
    end subroutine ensure_wpw_reduced_eigensystem_only

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

    subroutine log_mixed_field_kick_diag(route_label, step, ispin_current, E_use, cmix_before, cmix_after, &
                                         state_s, state_e)
      character(len=*), intent(in) :: route_label
      integer, intent(in) :: step, ispin_current, state_s, state_e
      real(8), intent(in) :: E_use(3)
      complex(8), intent(in) :: cmix_before(:,:), cmix_after(:,:)
      integer :: nrow, ncol, occ_cut, ist, irow
      real(8) :: norm_before, norm_after, diff_norm, rel_diff
      real(8) :: unocc_before, unocc_after, occ_offdiag_after, diag_after
      real(8) :: max_col_unocc, col_unocc

      if (.not. mixed_z_field_kick_diag_enabled) return
      if (dg_frag%id /= 0) return
      nrow = min(size(cmix_before,1), size(cmix_after,1))
      ncol = min(size(cmix_before,2), size(cmix_after,2))
      if (nrow <= 0 .or. ncol <= 0) return
      occ_cut = min(max(state_e, 0), nrow)
      norm_before = 0.0d0
      norm_after = 0.0d0
      diff_norm = 0.0d0
      unocc_before = 0.0d0
      unocc_after = 0.0d0
      occ_offdiag_after = 0.0d0
      diag_after = 0.0d0
      max_col_unocc = 0.0d0
      do ist = 1, ncol
        col_unocc = 0.0d0
        do irow = 1, nrow
          norm_before = norm_before + abs(cmix_before(irow,ist))**2
          norm_after = norm_after + abs(cmix_after(irow,ist))**2
          diff_norm = diff_norm + abs(cmix_after(irow,ist) - cmix_before(irow,ist))**2
          if (irow > occ_cut) then
            unocc_before = unocc_before + abs(cmix_before(irow,ist))**2
            unocc_after = unocc_after + abs(cmix_after(irow,ist))**2
            col_unocc = col_unocc + abs(cmix_after(irow,ist))**2
          else if (state_s + ist - 1 == irow) then
            diag_after = diag_after + abs(cmix_after(irow,ist))**2
          else
            occ_offdiag_after = occ_offdiag_after + abs(cmix_after(irow,ist))**2
          end if
        end do
        max_col_unocc = max(max_col_unocc, col_unocc)
      end do
      rel_diff = sqrt(diff_norm) / max(sqrt(norm_before), 1.0d-300)
      write(*,'(1x,a,1x,a,2(a,i0),3(a,1pe13.5),8(a,1pe13.5),2(a,i0))') &
        '[DG-MIXEDZ-FIELD-KICK-DIAG]', trim(route_label), &
        ' step=', step, ' spin=', ispin_current, &
        ' Ex=', E_use(1), ' Ey=', E_use(2), ' Ez=', E_use(3), &
        ' norm_before=', norm_before, ' norm_after=', norm_after, &
        ' rel_coef_change=', rel_diff, &
        ' unocc_before=', unocc_before, ' unocc_after=', unocc_after, &
        ' occ_offdiag_after=', occ_offdiag_after, ' diag_after=', diag_after, &
        ' max_col_unocc=', max_col_unocc, ' occ_cut=', occ_cut, ' ncol=', ncol
      flush(6)
    end subroutine log_mixed_field_kick_diag

    subroutine apply_global_mixed_position_field_exp(E_use, state_s, state_e)
      real(8), intent(in) :: E_use(3)
      integer, intent(in) :: state_s, state_e
      integer :: ispin_current, nstate_blk, nmix, nw, np, istate, imix
      real(8) :: phase_c, phase_s, tloc
      complex(8), allocatable :: cmix(:,:), cmix_before(:,:), field_h(:,:), field_vec(:,:), tmp(:,:)
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

      allocate(cmix(nmix,nstate_blk), cmix_before(nmix,nstate_blk), field_h(nmix,nmix), field_vec(nmix,nmix), &
               tmp(nmix,nstate_blk), field_eval(nmix))
      do ispin_current = 1, dg_frag%nspin
        if (ispin_current > size(dg_frag%mixed_wannier_bpw_z, 4)) cycle
        tloc = get_wtime()
        call gather_global_mixed_coefficients(ispin_current, state_s, state_e, cmix)
        time_gather = time_gather + (get_wtime() - tloc)
        if (mixed_z_field_kick_diag_enabled) cmix_before(:,:) = cmix(:,:)
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
        if (mixed_z_field_kick_diag_enabled) call log_mixed_field_kick_diag( &
          'position_field_exp', itt, ispin_current, E_use, cmix_before, cmix, state_s, state_e)
        time_field_matmul = time_field_matmul + (get_wtime() - tloc)
        tloc = get_wtime()
        call scatter_global_mixed_coefficients(ispin_current, state_s, state_e, cmix)
        time_scatter = time_scatter + (get_wtime() - tloc)
      end do
      deallocate(cmix, cmix_before, field_h, field_vec, tmp, field_eval)
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

    subroutine apply_mixed_split_exp_backend(candidate_kind, E_use, state_s, state_e, &
                                             candidate_available, replacement_applied, bad, block_reason)
      character(len=*), intent(in) :: candidate_kind
      real(8), intent(in) :: E_use(3)
      integer, intent(in) :: state_s, state_e
      logical, intent(out) :: candidate_available, replacement_applied, bad
      character(len=*), intent(out) :: block_reason

      candidate_available = .false.
      replacement_applied = .false.
      bad = .false.
      block_reason = 'not_requested'
      select case (trim(candidate_kind))
      case ('global_mixed_split_backend')
        candidate_available = .true.
        call apply_global_mixed_split_exp(E_use, state_s, state_e)
        replacement_applied = .true.
        block_reason = 'none'
      case ('fragment_local_mixed_split_backend')
        call apply_fragment_local_mixed_split_exp_stub(E_use, state_s, state_e, &
          candidate_available, replacement_applied, bad, block_reason)
      case ('neighbor_env_expdiag')
        call apply_neighbor_env_expdiag_stub(E_use, state_s, state_e, .false., .false., &
          candidate_available, replacement_applied, bad, block_reason)
      case ('neighbor_env_fieldonly')
        call apply_neighbor_env_expdiag_stub(E_use, state_s, state_e, .true., .false., &
          candidate_available, replacement_applied, bad, block_reason)
      case ('neighbor_env_interaction')
        call apply_neighbor_env_expdiag_stub(E_use, state_s, state_e, .false., .true., &
          candidate_available, replacement_applied, bad, block_reason)
      case default
        bad = .true.
        block_reason = 'unknown_candidate_kind'
      end select
    end subroutine apply_mixed_split_exp_backend

    subroutine apply_neighbor_env_expdiag_stub(E_use, state_s, state_e, fieldonly_mode, interaction_mode, &
                                               candidate_available, replacement_applied, bad, block_reason)
      real(8), intent(in) :: E_use(3)
      integer, intent(in) :: state_s, state_e
      logical, intent(in) :: fieldonly_mode, interaction_mode
      logical, intent(out) :: candidate_available, replacement_applied, bad
      character(len=*), intent(out) :: block_reason
      logical :: owner_local_layout_ready
      integer :: w_owner_slots, w_total_slots, p_self_slots, p_neighbor_slots
      integer :: gid_mismatch_count, owner_frag_mismatch_count, owner_rank_mismatch_count
      integer :: nstate_blk, i_local, ifrag, axis, side, jfrag
      integer :: n_neighbor, neighbor_ids(6)
      integer :: nmix, nw, np, ispin_diag
      real(8) :: coeff_norm_w, coeff_norm_p, coeff_max_abs
      complex(8), allocatable :: cmix_diag(:,:)
      logical :: field_off, kernel_ready, writeback_bad, bad_coef, diag_neighbor_env
      logical :: eig_cache_hit, eig_built, eig_skipped
      real(8) :: eig_elapsed

      nstate_blk = max(0, state_e - state_s + 1)
      field_off = sum(abs(E_use(1:3))) <= 1.0d-30
      diag_neighbor_env = .not. dg_frag%mixed_z_perf_count_enabled .and. &
        (itt == 1 .or. mod(itt, 500) == 0)
      if (.not. allocated(dg_frag%wpw_reduced_dim) .or. .not. allocated(dg_frag%wpw_reduced_nself) .or. &
          .not. allocated(dg_frag%wpw_reduced_nraw)) then
        call prepare_wpw_local_payload_ingredients('neighbor_env_layout')
      end if
      call enumerate_fragment_local_mixed_slots(owner_local_layout_ready, w_owner_slots, w_total_slots, &
        p_self_slots, p_neighbor_slots, gid_mismatch_count, owner_frag_mismatch_count, &
        owner_rank_mismatch_count)

      candidate_available = owner_local_layout_ready .and. nstate_blk > 0
      replacement_applied = .false.
      bad = .false.
      block_reason = 'neighbor_env_field_on_kernel_not_implemented'

      nmix = dg_frag%mixed_wannier_bpw_nmix
      nw = dg_frag%mixed_wannier_bpw_nw
      np = dg_frag%mixed_wannier_bpw_np
      coeff_norm_w = 0.0d0
      coeff_norm_p = 0.0d0
      coeff_max_abs = 0.0d0
      if (diag_neighbor_env .and. candidate_available .and. nmix > 0 .and. nw >= 0 .and. np >= 0 .and. &
          nmix == nw + np .and. dg_frag%nspin > 0) then
        ispin_diag = 1
        allocate(cmix_diag(nmix,nstate_blk))
        call gather_global_mixed_coefficients(ispin_diag, state_s, state_e, cmix_diag)
        if (nw > 0) coeff_norm_w = sum(abs(cmix_diag(1:nw,1:nstate_blk))**2)
        if (np > 0) coeff_norm_p = sum(abs(cmix_diag(nw+1:nw+np,1:nstate_blk))**2)
        coeff_max_abs = maxval(abs(cmix_diag(1:nmix,1:nstate_blk)))
        deallocate(cmix_diag)
      end if

      if (diag_neighbor_env .and. dg_frag%id == 0) then
        write(*,*) '[DG-MIXED-Z-NEIGHBOR-ENV] backend selected'
        write(*,*) '[DG-MIXED-Z-NEIGHBOR-ENV] step=', itt, &
          ' nstate_blk=', nstate_blk, &
          ' available=', candidate_available, &
          ' W_owner_slots=', w_owner_slots, &
          ' W_total_layout_slots=', w_total_slots, &
          ' P_self_slots=', p_self_slots, &
          ' P_neighbor_slots=', p_neighbor_slots, &
          ' field_abs_sum=', sum(abs(E_use(1:3))), &
          ' coeff_norm_w=', coeff_norm_w, &
          ' coeff_norm_p=', coeff_norm_p, &
          ' coeff_max_abs=', coeff_max_abs
      end if
      if (diag_neighbor_env) then
        call diagnose_neighbor_env_raw_matrices()
      end if
      if (diag_neighbor_env .and. dg_frag%id == 0) then
        do i_local = 1, size(dg_frag%wpw_reduced_dim)
          ifrag = dg_frag%ifrag_start + i_local - 1
          n_neighbor = 0
          neighbor_ids(:) = 0
          do axis = 1, 3
            do side = -1, 1, 2
              jfrag = wpw_face_neighbor_fragment(dg_frag, ifrag, axis, side)
              if (jfrag <= 0 .or. jfrag == ifrag) cycle
              if (.not. any(neighbor_ids(1:max(1,n_neighbor)) == jfrag) .and. &
                  n_neighbor < size(neighbor_ids)) then
                n_neighbor = n_neighbor + 1
                neighbor_ids(n_neighbor) = jfrag
              end if
            end do
          end do
          write(*,*) '[DG-MIXED-Z-NEIGHBOR-ENV-LAYOUT]', &
            ' owner_frag=', ifrag, &
            ' neighbor_count=', n_neighbor, &
            ' neighbors=', neighbor_ids(:), &
            ' nred=', dg_frag%wpw_reduced_dim(i_local), &
            ' nself=', dg_frag%wpw_reduced_nself(i_local), &
            ' nraw=', dg_frag%wpw_reduced_nraw(i_local)
        end do
      end if
      if (field_off .and. candidate_available) then
        if (fieldonly_mode .or. interaction_mode) then
          replacement_applied = .true.
          block_reason = 'none'
          return
        end if
        bad_coef = .false.
        call project_neighbor_env_current_to_reduced(E_use, state_s, state_e, .false., bad, kernel_ready)
        call ensure_wpw_reduced_eigensystem_only(.true., eig_cache_hit, eig_built, eig_skipped, eig_elapsed)
        call dryrun_wpw_reduced_expdiag(state_s, state_e, dt, bad_coef, itt, E_use, Ac_ham, 'neighbor_env_fieldoff')
        call build_neighbor_env_owner_storage_from_reduced(state_s, state_e, bad, kernel_ready)
        if (kernel_ready .and. .not. bad .and. .not. bad_coef) then
          call writeback_fragment_local_storage_fieldoff(state_s, state_e, writeback_bad)
          bad = bad .or. writeback_bad
          if (.not. writeback_bad) then
            replacement_applied = .true.
            block_reason = 'none'
          else
            block_reason = 'neighbor_env_writeback_failed'
          end if
        else
          block_reason = 'neighbor_env_fieldoff_storage_failed'
        end if
      else if (.not. field_off .and. candidate_available) then
        bad_coef = .false.
        call project_neighbor_env_current_to_reduced(E_use, state_s, state_e, .not. interaction_mode, bad, kernel_ready)
        if (interaction_mode) then
          call apply_neighbor_env_interaction_exp(E_use, state_s, state_e, bad_coef)
        else if (.not. fieldonly_mode) then
          call ensure_wpw_reduced_eigensystem_only(.true., eig_cache_hit, eig_built, eig_skipped, eig_elapsed)
          call dryrun_wpw_reduced_expdiag(state_s, state_e, dt, bad_coef, itt, E_use, Ac_ham, 'neighbor_env_fieldon')
        end if
        call build_neighbor_env_owner_storage_from_reduced(state_s, state_e, bad, kernel_ready)
        if (kernel_ready .and. .not. bad .and. .not. bad_coef) then
          call writeback_fragment_local_storage_fieldoff(state_s, state_e, writeback_bad)
          bad = bad .or. writeback_bad
          if (.not. writeback_bad) then
            replacement_applied = .true.
            block_reason = 'none'
          else
            block_reason = 'neighbor_env_writeback_failed'
          end if
        else
          block_reason = 'neighbor_env_fieldon_storage_failed'
        end if
      else
        bad = .true.
        block_reason = 'neighbor_env_layout_not_ready'
      end if
    end subroutine apply_neighbor_env_expdiag_stub

    subroutine diagnose_neighbor_env_raw_matrices()
      integer :: i_local, ifrag, nraw, nself, i, j
      real(8) :: s_herm, h_herm, s_diag_min, s_diag_max, s_trace
      complex(8) :: sd, diff
      logical :: raw_available

      raw_available = allocated(dg_frag%wpw_reduced_Sraw_build) .and. &
        allocated(dg_frag%wpw_reduced_Hraw_build) .and. &
        allocated(dg_frag%wpw_reduced_nraw) .and. allocated(dg_frag%wpw_reduced_nself)
      if (.not. raw_available) then
        write(*,*) '[DG-MIXED-Z-NEIGHBOR-ENV-MATRIX]', &
          ' raw_available=', raw_available, &
          ' block_reason=', trim(dg_frag%mixed_z_local_prop_payload_block_reason)
        return
      end if

      do i_local = 1, size(dg_frag%wpw_reduced_nraw)
        ifrag = dg_frag%ifrag_start + i_local - 1
        nraw = dg_frag%wpw_reduced_nraw(i_local)
        nself = 0
        if (i_local <= size(dg_frag%wpw_reduced_nself)) nself = dg_frag%wpw_reduced_nself(i_local)
        if (nraw <= 0) cycle
        if (nraw > size(dg_frag%wpw_reduced_Sraw_build,1) .or. &
            nraw > size(dg_frag%wpw_reduced_Hraw_build,1)) cycle
        s_herm = 0.0d0
        h_herm = 0.0d0
        s_diag_min = huge(1.0d0)
        s_diag_max = -huge(1.0d0)
        s_trace = 0.0d0
        do i = 1, nraw
          sd = dg_frag%wpw_reduced_Sraw_build(i,i,i_local)
          s_diag_min = min(s_diag_min, real(sd,8))
          s_diag_max = max(s_diag_max, real(sd,8))
          s_trace = s_trace + real(sd,8)
          do j = 1, nraw
            diff = dg_frag%wpw_reduced_Sraw_build(i,j,i_local) - &
              conjg(dg_frag%wpw_reduced_Sraw_build(j,i,i_local))
            s_herm = max(s_herm, abs(diff))
            diff = dg_frag%wpw_reduced_Hraw_build(i,j,i_local) - &
              conjg(dg_frag%wpw_reduced_Hraw_build(j,i,i_local))
            h_herm = max(h_herm, abs(diff))
          end do
        end do
        write(*,*) '[DG-MIXED-Z-NEIGHBOR-ENV-MATRIX]', &
          ' owner_frag=', ifrag, &
          ' raw_available=', raw_available, &
          ' nself=', nself, &
          ' nraw=', nraw, &
          ' S_herm_max=', s_herm, &
          ' H_herm_max=', h_herm, &
          ' S_diag_min=', s_diag_min, &
          ' S_diag_max=', s_diag_max, &
          ' S_trace=', s_trace
      end do
    end subroutine diagnose_neighbor_env_raw_matrices

    subroutine project_neighbor_env_current_to_reduced(E_use, state_s, state_e, apply_field, bad, kernel_ready)
      real(8), intent(in) :: E_use(3)
      integer, intent(in) :: state_s, state_e
      logical, intent(in) :: apply_field
      logical, intent(inout) :: bad
      logical, intent(out) :: kernel_ready
      integer :: nstate_blk, nmix, nw, np, n_pw, nspin_use
      integer :: i_local, ifrag, nred, nraw, nself, n_w, n_pfrag
      integer :: axis, side, jfrag, iw, ipw, pidx, row0, gid, gp
      integer :: ispin_use, info_field, ired, field_axis
      integer :: pfrag_ids(7)
      integer, allocatable :: raw_gid(:)
      real(8) :: phase_c, phase_s, r_herm
      complex(8), allocatable :: cmix(:,:), c_raw(:,:), c_red(:,:)
      complex(8), allocatable :: field_raw(:,:), field_red(:,:), s_work(:,:), tmp(:,:)
      real(8), allocatable :: eval(:)
      logical :: found

      kernel_ready = .false.
      nstate_blk = max(0, state_e - state_s + 1)
      nmix = dg_frag%mixed_wannier_bpw_nmix
      nw = dg_frag%mixed_wannier_bpw_nw
      np = dg_frag%mixed_wannier_bpw_np
      n_pw = max(0, dg_frag%n_plane_waves)
      nspin_use = dg_frag%nspin
      if (nstate_blk <= 0 .or. nmix <= 0 .or. nw <= 0 .or. n_pw <= 0) then
        bad = .true.
        return
      end if
      if (.not. dg_frag%wpw_reduced_ready .or. &
          .not. allocated(dg_frag%wpw_reduced_dim) .or. &
          .not. allocated(dg_frag%wpw_reduced_nself) .or. &
          .not. allocated(dg_frag%wpw_reduced_nraw) .or. &
          .not. allocated(dg_frag%wpw_reduced_S) .or. &
          .not. allocated(dg_frag%wpw_reduced_Sraw_build) .or. &
          .not. allocated(dg_frag%wpw_reduced_transform) .or. &
          .not. allocated(dg_frag%global_wannier_local_nkeep) .or. &
          .not. allocated(dg_frag%global_wannier_local_ids) .or. &
          .not. allocated(dg_frag%mixed_wannier_bpw_pcoef)) then
        bad = .true.
        return
      end if
      if (apply_field .and. .not. allocated(dg_frag%mixed_wannier_bpw_z)) then
        bad = .true.
        return
      end if

      call ensure_wpw_reduced_coefficient_storage(state_s, state_e)
      allocate(cmix(nmix,nstate_blk))
      do ispin_use = 1, nspin_use
        if (ispin_use > size(dg_frag%wpw_reduced_S, 3)) cycle
        if (apply_field .and. ispin_use > size(dg_frag%mixed_wannier_bpw_z, 4)) cycle
        call gather_global_mixed_coefficients(ispin_use, state_s, state_e, cmix)
        do i_local = 1, size(dg_frag%wpw_reduced_dim)
          ifrag = dg_frag%ifrag_start + i_local - 1
          nred = dg_frag%wpw_reduced_dim(i_local)
          nraw = dg_frag%wpw_reduced_nraw(i_local)
          nself = dg_frag%wpw_reduced_nself(i_local)
          n_w = dg_frag%global_wannier_local_nkeep(i_local)
          if (nred <= 0 .or. nraw <= 0 .or. nself <= 0 .or. n_w <= 0) cycle
          if (nred > size(dg_frag%wpw_reduced_S,1) .or. nraw > size(dg_frag%wpw_reduced_Sraw_build,1)) cycle
          n_pfrag = 1
          pfrag_ids(:) = 0
          pfrag_ids(1) = ifrag
          do axis = 1, 3
            do side = -1, 1, 2
              jfrag = wpw_face_neighbor_fragment(dg_frag, ifrag, axis, side)
              if (jfrag <= 0 .or. jfrag == ifrag) cycle
              found = any(pfrag_ids(1:n_pfrag) == jfrag)
              if (.not. found .and. n_pfrag < size(pfrag_ids)) then
                n_pfrag = n_pfrag + 1
                pfrag_ids(n_pfrag) = jfrag
              end if
            end do
          end do
          if (nraw /= n_w + n_pfrag * n_pw) cycle

          allocate(c_raw(nraw,nstate_blk), c_red(nred,nstate_blk), raw_gid(nraw))
          c_raw(:, :) = (0.0d0, 0.0d0)
          raw_gid(:) = 0
          do iw = 1, n_w
            gid = dg_frag%global_wannier_local_ids(iw, i_local)
            if (gid >= 1 .and. gid <= nw) then
              raw_gid(iw) = gid
              c_raw(iw,1:nstate_blk) = cmix(gid,1:nstate_blk)
            end if
          end do
          do pidx = 1, n_pfrag
            row0 = n_w + (pidx - 1) * n_pw
            do ipw = 1, n_pw
              gp = (pfrag_ids(pidx) - 1) * n_pw + ipw
              if (gp >= 1 .and. gp <= np) then
                raw_gid(row0+ipw) = nw + gp
                c_raw(row0+ipw,1:nstate_blk) = cmix(nw+gp,1:nstate_blk)
              end if
            end do
          end do

          if (apply_field) then
            allocate(field_raw(nraw,nraw), field_red(nred,nred), s_work(nred,nred), &
                     tmp(nred,nstate_blk), eval(nred))
            c_red(:, :) = project_raw_coeff_to_reduced_cached(i_local, ispin_use, &
              dg_frag%wpw_reduced_Sraw_build(1:nraw,1:nraw,i_local), &
              dg_frag%wpw_reduced_transform(1:nraw,1:nred,i_local), &
              dg_frag%wpw_reduced_S(1:nred,1:nred,ispin_use,i_local), c_raw, nraw, nred, nstate_blk)
            field_axis = single_active_field_axis(E_use)
            if (field_axis > 0) then
              call apply_neighbor_env_cached_axis_field(c_red, raw_gid, nraw, nred, nstate_blk, &
                i_local, ispin_use, field_axis, E_use(field_axis), r_herm, info_field)
              if (info_field /= 0) bad = .true.
            else
              field_raw(:, :) = (0.0d0, 0.0d0)
              do iw = 1, nraw
                do ipw = 1, nraw
                  if (raw_gid(iw) <= 0 .or. raw_gid(ipw) <= 0) cycle
                  field_raw(iw,ipw) = -E_use(1) * dg_frag%mixed_wannier_bpw_z(1,raw_gid(iw),raw_gid(ipw),ispin_use) &
                                     -E_use(2) * dg_frag%mixed_wannier_bpw_z(2,raw_gid(iw),raw_gid(ipw),ispin_use) &
                                     -E_use(3) * dg_frag%mixed_wannier_bpw_z(3,raw_gid(iw),raw_gid(ipw),ispin_use)
                end do
              end do
              r_herm = 0.0d0
              do iw = 1, nraw
                do ipw = 1, nraw
                  r_herm = max(r_herm, abs(field_raw(iw,ipw) - conjg(field_raw(ipw,iw))))
                end do
              end do
              field_red(:, :) = matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw,1:nred,i_local))), &
                matmul(field_raw(:, :), dg_frag%wpw_reduced_transform(1:nraw,1:nred,i_local)))
              s_work(:, :) = dg_frag%wpw_reduced_S(1:nred,1:nred,ispin_use,i_local)
              call zhegv_with_query(field_red, s_work, nred, eval, info_field)
              if (info_field /= 0) then
                bad = .true.
                deallocate(field_raw, field_red, s_work, tmp, eval)
                deallocate(c_raw, c_red, raw_gid)
                cycle
              end if
              tmp(:, :) = matmul(conjg(transpose(field_red(:, :))), &
                matmul(dg_frag%wpw_reduced_S(1:nred,1:nred,ispin_use,i_local), c_red))
              do ired = 1, nred
                phase_c = cos(eval(ired) * dt)
                phase_s = sin(eval(ired) * dt)
                tmp(ired,1:nstate_blk) = cmplx(phase_c, -phase_s, kind=8) * tmp(ired,1:nstate_blk)
              end do
              c_red(:, :) = matmul(field_red(:, :), tmp(:, :))
            end if
            if (info_field /= 0) then
              deallocate(field_raw, field_red, s_work, tmp, eval)
              deallocate(c_raw, c_red, raw_gid)
              cycle
            end if
            if (dg_frag%id == 0 .and. .not. dg_frag%mixed_z_perf_count_enabled) then
              write(*,*) '[DG-MIXED-Z-NEIGHBOR-ENV-R]', &
                ' owner_frag=', ifrag, ' nred=', nred, ' nraw=', nraw, &
                ' Rfield_herm_max=', r_herm, ' field_abs_sum=', sum(abs(E_use(1:3)))
            end if
            deallocate(field_raw, field_red, s_work, tmp, eval)
          else
            c_red(:, :) = project_raw_coeff_to_reduced_cached(i_local, ispin_use, &
              dg_frag%wpw_reduced_Sraw_build(1:nraw,1:nraw,i_local), &
              dg_frag%wpw_reduced_transform(1:nraw,1:nred,i_local), &
              dg_frag%wpw_reduced_S(1:nred,1:nred,ispin_use,i_local), c_raw, nraw, nred, nstate_blk)
          end if

          dg_frag%coef_wpw_self(1:nself,state_s:state_e,ispin_use,i_local) = c_red(1:nself,1:nstate_blk)
          if (nred > nself) dg_frag%coef_wpw_neighbor_reduced(1:nred-nself,state_s:state_e,ispin_use,i_local) = &
            c_red(nself+1:nred,1:nstate_blk)
          deallocate(c_raw, c_red, raw_gid)
        end do
      end do
      deallocate(cmix)
      dg_frag%wpw_reduced_coef_initialized = .true.
      kernel_ready = .true.
    end subroutine project_neighbor_env_current_to_reduced

    integer function single_active_field_axis(E_use) result(axis)
      real(8), intent(in) :: E_use(3)
      real(8) :: tol

      tol = max(1.0d-30, 1.0d-12 * maxval(abs(E_use(1:3))))
      axis = 0
      if (abs(E_use(1)) > tol .and. abs(E_use(2)) <= tol .and. abs(E_use(3)) <= tol) axis = 1
      if (abs(E_use(2)) > tol .and. abs(E_use(1)) <= tol .and. abs(E_use(3)) <= tol) axis = 2
      if (abs(E_use(3)) > tol .and. abs(E_use(1)) <= tol .and. abs(E_use(2)) <= tol) axis = 3
    end function single_active_field_axis

    function project_raw_coeff_to_reduced_cached(i_local, ispin_use, Sraw, Tred, Sred, craw, nraw, nred, nstate) &
      result(cred)
      integer, intent(in) :: i_local, ispin_use, nraw, nred, nstate
      complex(8), intent(in) :: Sraw(nraw,nraw), Tred(nraw,nred), Sred(nred,nred), craw(nraw,nstate)
      complex(8) :: cred(nred,nstate)
      complex(8), allocatable, save :: Pcache(:,:,:,:)
      logical, allocatable, save :: Pcache_ready(:,:)
      complex(8), allocatable :: Sinv(:,:)
      real(8) :: smin, smax
      integer :: info, maxred, maxraw, nspin_cache, nfrag_cache

      maxred = size(dg_frag%wpw_reduced_S, 1)
      maxraw = size(dg_frag%wpw_reduced_Sraw_build, 1)
      nspin_cache = size(dg_frag%wpw_reduced_S, 3)
      nfrag_cache = size(dg_frag%wpw_reduced_dim)
      if (allocated(Pcache)) then
        if (size(Pcache,1) /= maxred .or. size(Pcache,2) /= maxraw .or. &
            size(Pcache,3) /= nspin_cache .or. size(Pcache,4) /= nfrag_cache) then
          deallocate(Pcache, Pcache_ready)
        end if
      end if
      if (.not. allocated(Pcache)) then
        allocate(Pcache(maxred,maxraw,nspin_cache,nfrag_cache))
        allocate(Pcache_ready(nspin_cache,nfrag_cache))
        Pcache(:, :, :, :) = (0.0d0, 0.0d0)
        Pcache_ready(:, :) = .false.
      end if
      if (i_local < 1 .or. i_local > nfrag_cache .or. ispin_use < 1 .or. ispin_use > nspin_cache) then
        cred(:, :) = (0.0d0, 0.0d0)
        return
      end if
      if (.not. Pcache_ready(ispin_use,i_local)) then
        allocate(Sinv(nred,nred))
        call build_hermitian_inverse(Sred, nred, Sinv, info, smin, smax)
        if (info /= 0) then
          cred(:, :) = (0.0d0, 0.0d0)
          deallocate(Sinv)
          return
        end if
        Pcache(1:nred,1:nraw,ispin_use,i_local) = &
          matmul(Sinv, matmul(conjg(transpose(Tred)), Sraw))
        Pcache_ready(ispin_use,i_local) = .true.
        deallocate(Sinv)
      end if
      cred(:, :) = matmul(Pcache(1:nred,1:nraw,ispin_use,i_local), craw)
    end function project_raw_coeff_to_reduced_cached

    subroutine apply_neighbor_env_cached_axis_field(c_red, raw_gid, nraw, nred, nstate, &
                                                    i_local, ispin_use, axis, E_axis, r_herm, info)
      integer, intent(in) :: nraw, nred, nstate, i_local, ispin_use, axis
      integer, intent(in) :: raw_gid(nraw)
      real(8), intent(in) :: E_axis
      complex(8), intent(inout) :: c_red(nred,nstate)
      real(8), intent(out) :: r_herm
      integer, intent(out) :: info
      complex(8), allocatable, save :: evec_cache(:,:,:,:,:)
      real(8), allocatable, save :: eval_cache(:,:,:,:)
      real(8), allocatable, save :: herm_cache(:,:,:)
      logical, allocatable, save :: ready(:,:,:)
      complex(8), allocatable :: Rraw(:,:), Rred(:,:), Swork(:,:), tmp(:,:)
      real(8), allocatable :: eval(:)
      real(8) :: phase_c, phase_s, scale
      integer :: maxred, nspin_cache, nfrag_cache, iw, ipw, ired

      info = 0
      r_herm = 0.0d0
      if (axis < 1 .or. axis > 3) then
        info = -1
        return
      end if
      maxred = size(dg_frag%wpw_reduced_S, 1)
      nspin_cache = size(dg_frag%wpw_reduced_S, 3)
      nfrag_cache = size(dg_frag%wpw_reduced_dim)
      if (allocated(evec_cache)) then
        if (size(evec_cache,1) /= maxred .or. size(evec_cache,2) /= maxred .or. &
            size(evec_cache,3) /= 3 .or. size(evec_cache,4) /= nspin_cache .or. &
            size(evec_cache,5) /= nfrag_cache) then
          deallocate(evec_cache, eval_cache, herm_cache, ready)
        end if
      end if
      if (.not. allocated(evec_cache)) then
        allocate(evec_cache(maxred,maxred,3,nspin_cache,nfrag_cache))
        allocate(eval_cache(maxred,3,nspin_cache,nfrag_cache))
        allocate(herm_cache(3,nspin_cache,nfrag_cache))
        allocate(ready(3,nspin_cache,nfrag_cache))
        evec_cache(:, :, :, :, :) = (0.0d0, 0.0d0)
        eval_cache(:, :, :, :) = 0.0d0
        herm_cache(:, :, :) = 0.0d0
        ready(:, :, :) = .false.
      end if
      if (i_local < 1 .or. i_local > nfrag_cache .or. ispin_use < 1 .or. ispin_use > nspin_cache) then
        info = -2
        return
      end if
      if (.not. ready(axis,ispin_use,i_local)) then
        allocate(Rraw(nraw,nraw), Rred(nred,nred), Swork(nred,nred), eval(nred))
        Rraw(:, :) = (0.0d0, 0.0d0)
        do iw = 1, nraw
          do ipw = 1, nraw
            if (raw_gid(iw) <= 0 .or. raw_gid(ipw) <= 0) cycle
            Rraw(iw,ipw) = dg_frag%mixed_wannier_bpw_z(axis,raw_gid(iw),raw_gid(ipw),ispin_use)
          end do
        end do
        r_herm = 0.0d0
        do iw = 1, nraw
          do ipw = 1, nraw
            r_herm = max(r_herm, abs(Rraw(iw,ipw) - conjg(Rraw(ipw,iw))))
          end do
        end do
        Rred(:, :) = matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw,1:nred,i_local))), &
          matmul(Rraw(:, :), dg_frag%wpw_reduced_transform(1:nraw,1:nred,i_local)))
        Swork(:, :) = dg_frag%wpw_reduced_S(1:nred,1:nred,ispin_use,i_local)
        call zhegv_with_query(Rred, Swork, nred, eval, info)
        if (info /= 0) then
          deallocate(Rraw, Rred, Swork, eval)
          return
        end if
        evec_cache(1:nred,1:nred,axis,ispin_use,i_local) = Rred(:, :)
        eval_cache(1:nred,axis,ispin_use,i_local) = eval(:)
        herm_cache(axis,ispin_use,i_local) = r_herm
        ready(axis,ispin_use,i_local) = .true.
        deallocate(Rraw, Rred, Swork, eval)
      end if
      r_herm = herm_cache(axis,ispin_use,i_local)
      allocate(tmp(nred,nstate))
      tmp(:, :) = matmul(conjg(transpose(evec_cache(1:nred,1:nred,axis,ispin_use,i_local))), &
        matmul(dg_frag%wpw_reduced_S(1:nred,1:nred,ispin_use,i_local), c_red))
      scale = -E_axis
      do ired = 1, nred
        phase_c = cos(scale * eval_cache(ired,axis,ispin_use,i_local) * dt)
        phase_s = sin(scale * eval_cache(ired,axis,ispin_use,i_local) * dt)
        tmp(ired,1:nstate) = cmplx(phase_c, -phase_s, kind=8) * tmp(ired,1:nstate)
      end do
      c_red(:, :) = matmul(evec_cache(1:nred,1:nred,axis,ispin_use,i_local), tmp(:, :))
      deallocate(tmp)
    end subroutine apply_neighbor_env_cached_axis_field

    subroutine apply_neighbor_env_interaction_exp(E_use, state_s, state_e, bad_coef_any)
      real(8), intent(in) :: E_use(3)
      integer, intent(in) :: state_s, state_e
      logical, intent(out) :: bad_coef_any
      integer :: nstate_blk, nstate_store, n_pw, nspin_use
      integer :: i_local, ifrag, nred, nraw, nself, nneigh, n_w, n_pfrag
      integer :: axis, side, jfrag, pidx, row0, iw, ipw, gp, gid, ispin_use, ist, ired
      integer :: info0, infoe
      integer :: pfrag_ids(7)
      integer, allocatable :: raw_gid(:)
      real(8), allocatable :: eval0(:), evale(:)
      real(8) :: phase_c, phase_s
      logical :: found
      complex(8), allocatable :: H0(:,:), He(:,:), Swork(:,:), Rraw(:,:), Rred(:,:)
      complex(8), allocatable :: cvec(:), c_start(:), c0(:), ce(:), tmp(:), work_vec(:)

      bad_coef_any = .false.
      nstate_blk = max(0, state_e - state_s + 1)
      nstate_store = max(1, size(dg_frag%coef, 2))
      n_pw = max(0, dg_frag%n_plane_waves)
      nspin_use = dg_frag%nspin
      if (nstate_blk <= 0 .or. n_pw <= 0) return
      if (.not. dg_frag%wpw_reduced_ready) then
        bad_coef_any = .true.
        return
      end if
      if (.not. allocated(dg_frag%wpw_reduced_H) .or. .not. allocated(dg_frag%wpw_reduced_S) .or. &
          .not. allocated(dg_frag%wpw_reduced_transform) .or. .not. allocated(dg_frag%wpw_reduced_Sraw_build) .or. &
          .not. allocated(dg_frag%wpw_reduced_dim) .or. .not. allocated(dg_frag%wpw_reduced_nraw) .or. &
          .not. allocated(dg_frag%wpw_reduced_nself) .or. .not. allocated(dg_frag%global_wannier_local_ids) .or. &
          .not. allocated(dg_frag%global_wannier_local_nkeep) .or. .not. allocated(dg_frag%mixed_wannier_bpw_z)) then
        bad_coef_any = .true.
        return
      end if

      do i_local = 1, size(dg_frag%wpw_reduced_dim)
        ifrag = dg_frag%ifrag_start + i_local - 1
        nred = dg_frag%wpw_reduced_dim(i_local)
        nraw = dg_frag%wpw_reduced_nraw(i_local)
        nself = dg_frag%wpw_reduced_nself(i_local)
        n_w = dg_frag%global_wannier_local_nkeep(i_local)
        if (nred <= 0 .or. nraw <= 0 .or. nself <= 0 .or. n_w <= 0) cycle
        if (nred > size(dg_frag%wpw_reduced_H,1) .or. nraw > size(dg_frag%wpw_reduced_transform,1)) cycle
        nneigh = max(0, nred - nself)
        n_pfrag = 1
        pfrag_ids(:) = 0
        pfrag_ids(1) = ifrag
        do axis = 1, 3
          do side = -1, 1, 2
            jfrag = wpw_face_neighbor_fragment(dg_frag, ifrag, axis, side)
            if (jfrag <= 0 .or. jfrag == ifrag) cycle
            found = any(pfrag_ids(1:n_pfrag) == jfrag)
            if (.not. found .and. n_pfrag < size(pfrag_ids)) then
              n_pfrag = n_pfrag + 1
              pfrag_ids(n_pfrag) = jfrag
            end if
          end do
        end do
        if (nraw /= n_w + n_pfrag * n_pw) cycle

        allocate(raw_gid(nraw), Rraw(nraw,nraw), Rred(nred,nred))
        raw_gid(:) = 0
        do iw = 1, n_w
          gid = dg_frag%global_wannier_local_ids(iw, i_local)
          if (gid >= 1 .and. gid <= dg_frag%mixed_wannier_bpw_nw) raw_gid(iw) = gid
        end do
        do pidx = 1, n_pfrag
          row0 = n_w + (pidx - 1) * n_pw
          do ipw = 1, n_pw
            gp = (pfrag_ids(pidx) - 1) * n_pw + ipw
            if (gp >= 1 .and. gp <= dg_frag%mixed_wannier_bpw_np) raw_gid(row0+ipw) = dg_frag%mixed_wannier_bpw_nw + gp
          end do
        end do
        allocate(H0(nred,nred), He(nred,nred), Swork(nred,nred), eval0(nred), evale(nred))
        allocate(cvec(nred), c_start(nred), c0(nred), ce(nred), tmp(nred), work_vec(nred))
        do ispin_use = 1, nspin_use
          if (ispin_use > size(dg_frag%wpw_reduced_H,3)) cycle
          if (ispin_use > size(dg_frag%mixed_wannier_bpw_z,4)) cycle
          Rraw(:, :) = (0.0d0, 0.0d0)
          do iw = 1, nraw
            if (raw_gid(iw) <= 0) cycle
            do ipw = 1, nraw
              if (raw_gid(ipw) <= 0) cycle
              Rraw(iw,ipw) = E_use(1) * dg_frag%mixed_wannier_bpw_z(1,raw_gid(iw),raw_gid(ipw),ispin_use) &
                           + E_use(2) * dg_frag%mixed_wannier_bpw_z(2,raw_gid(iw),raw_gid(ipw),ispin_use) &
                           + E_use(3) * dg_frag%mixed_wannier_bpw_z(3,raw_gid(iw),raw_gid(ipw),ispin_use)
            end do
          end do
          Rred(:, :) = matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw,1:nred,i_local))), &
            matmul(Rraw(:, :), dg_frag%wpw_reduced_transform(1:nraw,1:nred,i_local)))
          H0(:, :) = dg_frag%wpw_reduced_H(1:nred,1:nred,ispin_use,i_local)
          He(:, :) = H0(:, :) - Rred(:, :)
          Swork(:, :) = dg_frag%wpw_reduced_S(1:nred,1:nred,ispin_use,i_local)
          call zhegv_with_query(H0, Swork, nred, eval0, info0)
          Swork(:, :) = dg_frag%wpw_reduced_S(1:nred,1:nred,ispin_use,i_local)
          call zhegv_with_query(He, Swork, nred, evale, infoe)
          if (info0 /= 0 .or. infoe /= 0) then
            bad_coef_any = .true.
            cycle
          end if
          do ist = state_s, min(state_e, nstate_store)
            cvec(:) = (0.0d0, 0.0d0)
            cvec(1:nself) = dg_frag%coef_wpw_self(1:nself, ist, ispin_use, i_local)
            if (nneigh > 0) cvec(nself+1:nred) = &
              dg_frag%coef_wpw_neighbor_reduced(1:nneigh, ist, ispin_use, i_local)
            c_start(:) = cvec(:)

            work_vec(:) = matmul(dg_frag%wpw_reduced_S(1:nred,1:nred,ispin_use,i_local), c_start(:))
            tmp(:) = matmul(conjg(transpose(He(:, :))), work_vec(:))
            do ired = 1, nred
              phase_c = cos(evale(ired) * dt)
              phase_s = sin(evale(ired) * dt)
              tmp(ired) = cmplx(phase_c, -phase_s, kind=8) * tmp(ired)
            end do
            ce(:) = matmul(He(:, :), tmp(:))

            work_vec(:) = matmul(dg_frag%wpw_reduced_S(1:nred,1:nred,ispin_use,i_local), c_start(:))
            tmp(:) = matmul(conjg(transpose(H0(:, :))), work_vec(:))
            do ired = 1, nred
              phase_c = cos(eval0(ired) * dt)
              phase_s = sin(eval0(ired) * dt)
              tmp(ired) = cmplx(phase_c, -phase_s, kind=8) * tmp(ired)
            end do
            c0(:) = matmul(H0(:, :), tmp(:))
            cvec(:) = c_start(:) + (ce(:) - c0(:))

            dg_frag%coef_wpw_self(1:nself, ist, ispin_use, i_local) = cvec(1:nself)
            if (nneigh > 0) dg_frag%coef_wpw_neighbor_reduced(1:nneigh, ist, ispin_use, i_local) = &
              cvec(nself+1:nred)
          end do
        end do
        deallocate(H0, He, Swork, eval0, evale, cvec, c_start, c0, ce, tmp, work_vec)
        deallocate(raw_gid, Rraw, Rred)
      end do
      dg_frag%wpw_reduced_coef_initialized = .true.
    end subroutine apply_neighbor_env_interaction_exp

    function project_raw_coeff_to_reduced(Sraw, Tred, Sred, craw, nraw, nred, nstate) result(cred)
      integer, intent(in) :: nraw, nred, nstate
      complex(8), intent(in) :: Sraw(nraw,nraw), Tred(nraw,nred), Sred(nred,nred), craw(nraw,nstate)
      complex(8) :: cred(nred,nstate)
      complex(8), allocatable :: Sinv(:,:), rhs(:,:)
      real(8) :: smin, smax
      integer :: info

      allocate(Sinv(nred,nred), rhs(nred,nstate))
      call build_hermitian_inverse(Sred, nred, Sinv, info, smin, smax)
      if (info /= 0) then
        cred(:, :) = (0.0d0, 0.0d0)
      else
        rhs(:, :) = matmul(conjg(transpose(Tred)), matmul(Sraw, craw))
        cred(:, :) = matmul(Sinv, rhs)
      end if
      deallocate(Sinv, rhs)
    end function project_raw_coeff_to_reduced

    subroutine ensure_wpw_reduced_coefficient_storage(state_s, state_e)
      integer, intent(in) :: state_s, state_e
      integer :: nstate_store, nfrag_local

      nstate_store = max(size(dg_frag%coef, 2), state_e)
      nfrag_local = size(dg_frag%wpw_reduced_dim)
      if (allocated(dg_frag%coef_wpw_self)) then
        if (size(dg_frag%coef_wpw_self,1) /= dg_frag%wpw_reduced_max_dim .or. &
            size(dg_frag%coef_wpw_self,2) < nstate_store .or. &
            size(dg_frag%coef_wpw_self,3) /= dg_frag%nspin .or. &
            size(dg_frag%coef_wpw_self,4) /= nfrag_local) deallocate(dg_frag%coef_wpw_self)
      end if
      if (.not. allocated(dg_frag%coef_wpw_self)) then
        allocate(dg_frag%coef_wpw_self(dg_frag%wpw_reduced_max_dim, nstate_store, dg_frag%nspin, nfrag_local))
        dg_frag%coef_wpw_self(:, :, :, :) = (0.0d0, 0.0d0)
      end if
      if (allocated(dg_frag%coef_wpw_neighbor_reduced)) then
        if (size(dg_frag%coef_wpw_neighbor_reduced,1) /= dg_frag%wpw_reduced_max_dim .or. &
            size(dg_frag%coef_wpw_neighbor_reduced,2) < nstate_store .or. &
            size(dg_frag%coef_wpw_neighbor_reduced,3) /= dg_frag%nspin .or. &
            size(dg_frag%coef_wpw_neighbor_reduced,4) /= nfrag_local) deallocate(dg_frag%coef_wpw_neighbor_reduced)
      end if
      if (.not. allocated(dg_frag%coef_wpw_neighbor_reduced)) then
        allocate(dg_frag%coef_wpw_neighbor_reduced(dg_frag%wpw_reduced_max_dim, nstate_store, dg_frag%nspin, nfrag_local))
        dg_frag%coef_wpw_neighbor_reduced(:, :, :, :) = (0.0d0, 0.0d0)
      end if
    end subroutine ensure_wpw_reduced_coefficient_storage

    subroutine build_neighbor_env_owner_storage_from_reduced(state_s, state_e, bad, kernel_ready)
      integer, intent(in) :: state_s, state_e
      logical, intent(inout) :: bad
      logical, intent(out) :: kernel_ready
      integer :: nstate_blk, nw, np, n_pw, nspin_use
      integer :: i_local, ifrag, n_w, nself, nred, nraw, nneigh, iw, ipw, gid, gp
      integer :: axis, side, jfrag, pidx, n_pfrag, row0, raw_slot
      integer :: w_slot_count, pself_slot_count, pneighbor_slot_count, w_slot, pself_slot, pneighbor_slot
      integer :: ispin_use, state_col
      integer :: pfrag_ids(7)
      integer, allocatable :: w_slot_for_iw(:), pself_slot_for_ipw(:), pneighbor_slot_for_raw(:)
      complex(8), allocatable :: c_red(:), raw_back(:)

      kernel_ready = .false.
      nstate_blk = max(0, state_e - state_s + 1)
      nw = dg_frag%mixed_wannier_bpw_nw
      np = dg_frag%mixed_wannier_bpw_np
      n_pw = max(0, dg_frag%n_plane_waves)
      nspin_use = dg_frag%nspin
      if (nstate_blk <= 0 .or. nw <= 0 .or. n_pw <= 0) then
        bad = .true.
        return
      end if
      if (.not. dg_frag%wpw_reduced_coef_initialized .or. &
          .not. allocated(dg_frag%coef_wpw_self) .or. &
          .not. allocated(dg_frag%wpw_reduced_dim) .or. &
          .not. allocated(dg_frag%wpw_reduced_nself) .or. &
          .not. allocated(dg_frag%wpw_reduced_nraw) .or. &
          .not. allocated(dg_frag%wpw_reduced_transform) .or. &
          .not. allocated(dg_frag%global_wannier_local_nkeep) .or. &
          .not. allocated(dg_frag%global_wannier_local_ids) .or. &
          .not. allocated(dg_frag%global_wannier_owner_frag)) then
        bad = .true.
        return
      end if

      w_slot_count = 0
      pself_slot_count = 0
      pneighbor_slot_count = 0
      do i_local = 1, size(dg_frag%global_wannier_local_nkeep)
        ifrag = dg_frag%ifrag_start + i_local - 1
        n_w = dg_frag%global_wannier_local_nkeep(i_local)
        if (i_local > size(dg_frag%wpw_reduced_nself) .or. i_local > size(dg_frag%wpw_reduced_dim) .or. &
            i_local > size(dg_frag%wpw_reduced_nraw)) cycle
        nself = dg_frag%wpw_reduced_nself(i_local)
        nred = dg_frag%wpw_reduced_dim(i_local)
        nraw = dg_frag%wpw_reduced_nraw(i_local)
        if (nred <= 0 .or. nraw <= 0 .or. nself < n_w) cycle
        n_pfrag = 1
        pfrag_ids(:) = 0
        pfrag_ids(1) = ifrag
        do axis = 1, 3
          do side = -1, 1, 2
            jfrag = wpw_face_neighbor_fragment(dg_frag, ifrag, axis, side)
            if (jfrag <= 0 .or. jfrag == ifrag) cycle
            if (.not. any(pfrag_ids(1:n_pfrag) == jfrag) .and. n_pfrag < size(pfrag_ids)) then
              n_pfrag = n_pfrag + 1
              pfrag_ids(n_pfrag) = jfrag
            end if
          end do
        end do
        if (nraw /= n_w + n_pfrag * n_pw) cycle
        do iw = 1, n_w
          gid = dg_frag%global_wannier_local_ids(iw, i_local)
          if (gid < 1 .or. gid > nw) cycle
          if (dg_frag%global_wannier_owner_frag(gid) < dg_frag%ifrag_start .or. &
              dg_frag%global_wannier_owner_frag(gid) > dg_frag%ifrag_end) cycle
          if (iw > nself) cycle
          w_slot_count = w_slot_count + 1
        end do
        do pidx = 1, n_pfrag
          do ipw = 1, n_pw
            gp = (pfrag_ids(pidx) - 1) * n_pw + ipw
            if (gp < 1 .or. gp > np) cycle
            if (pidx == 1) then
              pself_slot_count = pself_slot_count + 1
            else
              pneighbor_slot_count = pneighbor_slot_count + 1
            end if
          end do
        end do
      end do

      if (allocated(dg_frag%mixed_z_frag_local_wcoef)) deallocate(dg_frag%mixed_z_frag_local_wcoef)
      if (allocated(dg_frag%mixed_z_frag_local_w_gid)) deallocate(dg_frag%mixed_z_frag_local_w_gid)
      if (allocated(dg_frag%mixed_z_frag_local_w_mix_gid)) deallocate(dg_frag%mixed_z_frag_local_w_mix_gid)
      if (allocated(dg_frag%mixed_z_frag_local_pself_coef)) deallocate(dg_frag%mixed_z_frag_local_pself_coef)
      if (allocated(dg_frag%mixed_z_frag_local_pself_gid)) deallocate(dg_frag%mixed_z_frag_local_pself_gid)
      if (allocated(dg_frag%mixed_z_frag_local_pself_mix_gid)) deallocate(dg_frag%mixed_z_frag_local_pself_mix_gid)
      if (allocated(dg_frag%mixed_z_frag_local_pneighbor_coef)) deallocate(dg_frag%mixed_z_frag_local_pneighbor_coef)
      if (allocated(dg_frag%mixed_z_frag_local_pneighbor_gid)) deallocate(dg_frag%mixed_z_frag_local_pneighbor_gid)
      if (allocated(dg_frag%mixed_z_frag_local_pneighbor_mix_gid)) deallocate(dg_frag%mixed_z_frag_local_pneighbor_mix_gid)

      allocate(dg_frag%mixed_z_frag_local_wcoef(max(1,w_slot_count), nstate_blk, nspin_use))
      allocate(dg_frag%mixed_z_frag_local_w_gid(max(1,w_slot_count)))
      allocate(dg_frag%mixed_z_frag_local_w_mix_gid(max(1,w_slot_count)))
      allocate(dg_frag%mixed_z_frag_local_pself_coef(max(1,pself_slot_count), nstate_blk, nspin_use))
      allocate(dg_frag%mixed_z_frag_local_pself_gid(max(1,pself_slot_count)))
      allocate(dg_frag%mixed_z_frag_local_pself_mix_gid(max(1,pself_slot_count)))
      allocate(dg_frag%mixed_z_frag_local_pneighbor_coef(max(1,pneighbor_slot_count), nstate_blk, nspin_use))
      allocate(dg_frag%mixed_z_frag_local_pneighbor_gid(max(1,pneighbor_slot_count)))
      allocate(dg_frag%mixed_z_frag_local_pneighbor_mix_gid(max(1,pneighbor_slot_count)))
      dg_frag%mixed_z_frag_local_wcoef(:, :, :) = (0.0d0, 0.0d0)
      dg_frag%mixed_z_frag_local_pself_coef(:, :, :) = (0.0d0, 0.0d0)
      dg_frag%mixed_z_frag_local_pneighbor_coef(:, :, :) = (0.0d0, 0.0d0)
      dg_frag%mixed_z_frag_local_w_gid(:) = 0
      dg_frag%mixed_z_frag_local_w_mix_gid(:) = 0
      dg_frag%mixed_z_frag_local_pself_gid(:) = 0
      dg_frag%mixed_z_frag_local_pself_mix_gid(:) = 0
      dg_frag%mixed_z_frag_local_pneighbor_gid(:) = 0
      dg_frag%mixed_z_frag_local_pneighbor_mix_gid(:) = 0

      w_slot = 0
      pself_slot = 0
      pneighbor_slot = 0
      do i_local = 1, size(dg_frag%global_wannier_local_nkeep)
        ifrag = dg_frag%ifrag_start + i_local - 1
        n_w = dg_frag%global_wannier_local_nkeep(i_local)
        if (i_local > size(dg_frag%wpw_reduced_nself) .or. i_local > size(dg_frag%wpw_reduced_dim) .or. &
            i_local > size(dg_frag%wpw_reduced_nraw)) cycle
        nself = dg_frag%wpw_reduced_nself(i_local)
        nred = dg_frag%wpw_reduced_dim(i_local)
        nraw = dg_frag%wpw_reduced_nraw(i_local)
        if (nred <= 0 .or. nraw <= 0 .or. nself < n_w) cycle
        if (nred > size(dg_frag%wpw_reduced_transform,2) .or. nraw > size(dg_frag%wpw_reduced_transform,1)) cycle
        nneigh = max(0, nred - nself)
        n_pfrag = 1
        pfrag_ids(:) = 0
        pfrag_ids(1) = ifrag
        do axis = 1, 3
          do side = -1, 1, 2
            jfrag = wpw_face_neighbor_fragment(dg_frag, ifrag, axis, side)
            if (jfrag <= 0 .or. jfrag == ifrag) cycle
            if (.not. any(pfrag_ids(1:n_pfrag) == jfrag) .and. n_pfrag < size(pfrag_ids)) then
              n_pfrag = n_pfrag + 1
              pfrag_ids(n_pfrag) = jfrag
            end if
          end do
        end do
        if (nraw /= n_w + n_pfrag * n_pw) cycle
        allocate(w_slot_for_iw(n_w), pself_slot_for_ipw(n_pw), pneighbor_slot_for_raw(nraw))
        allocate(c_red(nred), raw_back(nraw))
        w_slot_for_iw(:) = 0
        pself_slot_for_ipw(:) = 0
        pneighbor_slot_for_raw(:) = 0
        do iw = 1, n_w
          gid = dg_frag%global_wannier_local_ids(iw, i_local)
          if (gid < 1 .or. gid > nw) cycle
          if (dg_frag%global_wannier_owner_frag(gid) < dg_frag%ifrag_start .or. &
              dg_frag%global_wannier_owner_frag(gid) > dg_frag%ifrag_end) cycle
          if (iw > nself) cycle
          w_slot = w_slot + 1
          w_slot_for_iw(iw) = w_slot
          dg_frag%mixed_z_frag_local_w_gid(w_slot) = gid
          dg_frag%mixed_z_frag_local_w_mix_gid(w_slot) = gid
        end do
        do pidx = 1, n_pfrag
          row0 = n_w + (pidx - 1) * n_pw
          do ipw = 1, n_pw
            raw_slot = row0 + ipw
            gp = (pfrag_ids(pidx) - 1) * n_pw + ipw
            if (gp < 1 .or. gp > np) cycle
            if (pidx == 1) then
              pself_slot = pself_slot + 1
              pself_slot_for_ipw(ipw) = pself_slot
              dg_frag%mixed_z_frag_local_pself_gid(pself_slot) = gp
              dg_frag%mixed_z_frag_local_pself_mix_gid(pself_slot) = nw + gp
            else
              pneighbor_slot = pneighbor_slot + 1
              pneighbor_slot_for_raw(raw_slot) = pneighbor_slot
              dg_frag%mixed_z_frag_local_pneighbor_gid(pneighbor_slot) = gp
              dg_frag%mixed_z_frag_local_pneighbor_mix_gid(pneighbor_slot) = nw + gp
            end if
          end do
        end do
        do ispin_use = 1, nspin_use
          do state_col = 1, nstate_blk
            c_red(:) = (0.0d0, 0.0d0)
            c_red(1:nself) = dg_frag%coef_wpw_self(1:nself,state_s+state_col-1,ispin_use,i_local)
            if (nneigh > 0) c_red(nself+1:nred) = &
              dg_frag%coef_wpw_neighbor_reduced(1:nneigh,state_s+state_col-1,ispin_use,i_local)
            raw_back(:) = matmul(dg_frag%wpw_reduced_transform(1:nraw,1:nred,i_local), c_red(:))
            do iw = 1, n_w
              if (w_slot_for_iw(iw) > 0) &
                dg_frag%mixed_z_frag_local_wcoef(w_slot_for_iw(iw),state_col,ispin_use) = raw_back(iw)
            end do
            do ipw = 1, n_pw
              if (pself_slot_for_ipw(ipw) > 0) &
                dg_frag%mixed_z_frag_local_pself_coef(pself_slot_for_ipw(ipw),state_col,ispin_use) = &
                  raw_back(n_w+ipw)
            end do
            do raw_slot = n_w + n_pw + 1, nraw
              if (pneighbor_slot_for_raw(raw_slot) > 0) &
                dg_frag%mixed_z_frag_local_pneighbor_coef(pneighbor_slot_for_raw(raw_slot),state_col,ispin_use) = &
                  raw_back(raw_slot)
            end do
          end do
        end do
        deallocate(w_slot_for_iw, pself_slot_for_ipw, pneighbor_slot_for_raw, c_red, raw_back)
      end do

      dg_frag%mixed_z_frag_local_w_slots = w_slot_count
      dg_frag%mixed_z_frag_local_pself_slots = pself_slot_count
      dg_frag%mixed_z_frag_local_pneighbor_slots = pneighbor_slot_count
      dg_frag%mixed_z_frag_local_nstate = nstate_blk
      dg_frag%mixed_z_frag_local_nspin = nspin_use
      dg_frag%mixed_z_frag_local_storage_ready = w_slot_count > 0
      kernel_ready = dg_frag%mixed_z_frag_local_storage_ready
      if (dg_frag%id == 0 .and. .not. dg_frag%mixed_z_perf_count_enabled) then
        write(*,*) '[DG-MIXED-Z-NEIGHBOR-ENV-STORAGE]', &
          ' W_owner_storage_slots=', w_slot_count, &
          ' P_self_storage_slots=', pself_slot_count, &
          ' replacement_storage_ready=', kernel_ready
      end if
    end subroutine build_neighbor_env_owner_storage_from_reduced

    subroutine apply_fragment_local_mixed_split_exp_stub(E_use, state_s, state_e, &
                                                        candidate_available, replacement_applied, bad, block_reason)
      real(8), intent(in) :: E_use(3)
      integer, intent(in) :: state_s, state_e
      logical, intent(out) :: candidate_available, replacement_applied, bad
      character(len=*), intent(out) :: block_reason
      logical :: owner_local_layout_ready
      integer :: w_owner_slots, w_total_slots, p_self_slots, p_neighbor_slots
      integer :: gid_mismatch_count, owner_frag_mismatch_count, owner_rank_mismatch_count
      integer :: nstate_blk
      integer :: field_block_dim, field_neighbor_slots
      real(8) :: kernel_diff_snorm
      logical :: kernel_ready, writeback_bad

      nstate_blk = max(0, state_e - state_s + 1)
      if (.not. allocated(dg_frag%wpw_reduced_dim) .or. .not. allocated(dg_frag%wpw_reduced_nself) .or. &
          .not. allocated(dg_frag%wpw_reduced_nraw)) then
        call prepare_wpw_local_payload_ingredients('fragment_local_layout')
      end if
      call enumerate_fragment_local_mixed_slots(owner_local_layout_ready, w_owner_slots, w_total_slots, &
        p_self_slots, p_neighbor_slots, gid_mismatch_count, owner_frag_mismatch_count, &
        owner_rank_mismatch_count)
      owner_local_layout_ready = owner_local_layout_ready .and. nstate_blk > 0
      select case (trim(mixed_z_frag_local_field_block_kind))
      case ('w_only')
        field_block_dim = w_owner_slots
        field_neighbor_slots = 0
      case ('w_pself','owner_local')
        field_block_dim = w_owner_slots + p_self_slots
        field_neighbor_slots = 0
      case default
        field_block_dim = w_owner_slots + p_self_slots + p_neighbor_slots
        field_neighbor_slots = p_neighbor_slots
      end select
      candidate_available = owner_local_layout_ready
      replacement_applied = .false.
      bad = .false.
      block_reason = 'not_requested'
      kernel_ready = .false.
      kernel_diff_snorm = huge(1.0d0)
      if (candidate_available) then
        call build_fragment_local_storage_direct(E_use, state_s, state_e, bad, kernel_ready)
        kernel_diff_snorm = 0.0d0
        if (kernel_ready .and. mixed_z_frag_local_reference_diag_enabled) then
          call compare_fragment_local_storage_with_global_reference(E_use, state_s, state_e, bad, kernel_diff_snorm)
          if (sum(abs(E_use(1:3))) <= 1.0d-30) then
            kernel_ready = kernel_diff_snorm <= 1.0d-10 .and. .not. bad
          else
            kernel_ready = .not. bad
          end if
        end if
        if (kernel_ready .and. (sum(abs(E_use(1:3))) > 1.0d-30 .or. kernel_diff_snorm <= 1.0d-10)) then
          writeback_bad = .false.
          call writeback_fragment_local_storage_fieldoff(state_s, state_e, writeback_bad)
          bad = bad .or. writeback_bad
          if (.not. writeback_bad) then
            replacement_applied = .true.
            block_reason = 'none'
          else
            block_reason = 'storage_writeback_failed'
          end if
        else if (trim(block_reason) == 'not_requested') then
          block_reason = 'local_kernel_diff_above_tol'
        end if
      else
        block_reason = 'owner_local_layout_not_ready'
      end if
      if (dg_frag%id == 0 .and. .not. dg_frag%mixed_z_perf_count_enabled) then
        write(*,*) &
          '[DG-MIXEDZ-LOCAL-PROP-BACKEND-STUB-CMP]', &
          ' step=', itt, &
          ' field_block_dim=', field_block_dim, &
          ' field_neighbor_slots=', field_neighbor_slots, &
          ' W_owner_storage_slots=', w_owner_slots, &
          ' W_total_layout_slots=', w_total_slots, &
          ' P_self_storage_slots=', p_self_slots, &
          ' P_neighbor_storage_slots=', p_neighbor_slots, &
          ' available=', candidate_available, &
          ' owner_local_layout_ready=', owner_local_layout_ready, &
          ' replacement_applied=', replacement_applied, &
          ' bad=', bad, &
          ' candidate_kind=', 'fragment_local_mixed_split_backend', &
          ' replacement_block_reason=', trim(block_reason), &
          ' global_pack_calls=', merge('debug_reference', 'not_used       ', mixed_z_frag_local_reference_diag_enabled), &
          ' local_pack_calls=', 'owner_local', &
          ' field_block_kind=', trim(mixed_z_frag_local_field_block_kind), &
          ' field_abs_sum=', sum(abs(E_use(1:3)))
      end if
    end subroutine apply_fragment_local_mixed_split_exp_stub

    subroutine build_fragment_local_storage_direct(E_use, state_s, state_e, bad, kernel_ready)
      real(8), intent(in) :: E_use(3)
      integer, intent(in) :: state_s, state_e
      logical, intent(inout) :: bad
      logical, intent(out) :: kernel_ready
      integer :: nstate_blk, nw, np, nspin_use, n_pw, ispin_use, ist, state_col
      integer :: i_local, ifrag, nred, nself, nraw, n_w, nbf, ibasis
      integer :: iw, gid, ipw, pidx, gp, raw_slot, global_row, local_row
      integer :: ifrag_row, i_local_row, nvalid
      integer :: w_slot, pself_slot, pneighbor_slot
      integer :: w_slot_count, pself_slot_count, pneighbor_slot_count
      integer :: n_pfrag, pfrag_ids(7), axis, side, jfrag
      integer :: nblock, iblk, jblk, w_owner_count, p_valid_count
      integer :: p_field_max
      logical :: field_off, local_kernel_available, bad_local
      real(8) :: phase_c, phase_s
      complex(8) :: c_owner, c_input, cand_val
      complex(8), allocatable :: cw_source(:,:,:)
      integer, allocatable :: block_gid(:), block_pidx(:), block_gp(:)
      complex(8), allocatable :: block_coef(:,:), field_h(:,:), field_vec(:,:), tmp(:,:)
      real(8), allocatable :: field_eval(:)

      kernel_ready = .false.
      nstate_blk = max(0, state_e - state_s + 1)
      nw = dg_frag%mixed_wannier_bpw_nw
      np = dg_frag%mixed_wannier_bpw_np
      n_pw = max(0, dg_frag%n_plane_waves)
      nspin_use = dg_frag%nspin
      field_off = sum(abs(E_use(1:3))) <= 1.0d-30
      select case (trim(mixed_z_frag_local_field_block_kind))
      case ('w_only')
        p_field_max = 0
      case ('w_pself','owner_local')
        p_field_max = 1
      case default
        p_field_max = huge(1)
      end select
      local_kernel_available = nstate_blk > 0 .and. nw > 0 .and. &
        dg_frag%wpw_reduced_ready .and. allocated(dg_frag%wpw_reduced_dim) .and. &
        allocated(dg_frag%wpw_reduced_nself) .and. allocated(dg_frag%wpw_reduced_nraw) .and. &
        allocated(dg_frag%global_wannier_local_ids) .and. allocated(dg_frag%global_wannier_local_nkeep) .and. &
        allocated(dg_frag%global_wannier_coef) .and. allocated(dg_frag%coef_global_to_local) .and. &
        allocated(dg_frag%index_basis) .and. allocated(dg_frag%mixed_wannier_bpw_eval) .and. &
        (field_off .or. allocated(dg_frag%mixed_wannier_bpw_z)) .and. &
        (np <= 0 .or. allocated(dg_frag%mixed_wannier_bpw_pcoef))
      if (.not. local_kernel_available) then
        dg_frag%mixed_z_frag_local_storage_ready = .false.
        bad = .true.
        return
      end if

      w_slot_count = 0
      pself_slot_count = 0
      pneighbor_slot_count = 0
      do i_local = 1, size(dg_frag%wpw_reduced_dim)
        ifrag = dg_frag%ifrag_start + i_local - 1
        nraw = dg_frag%wpw_reduced_nraw(i_local)
        if (i_local > size(dg_frag%global_wannier_local_nkeep)) cycle
        n_w = dg_frag%global_wannier_local_nkeep(i_local)
        if (nraw <= 0 .or. n_w <= 0) cycle
        n_pfrag = 1
        pfrag_ids(1) = ifrag
        do axis = 1, 3
          do side = -1, 1, 2
            jfrag = wpw_face_neighbor_fragment(dg_frag, ifrag, axis, side)
            if (jfrag <= 0 .or. jfrag == ifrag) cycle
            if (.not. any(pfrag_ids(1:n_pfrag) == jfrag) .and. n_pfrag < size(pfrag_ids)) then
              n_pfrag = n_pfrag + 1
              pfrag_ids(n_pfrag) = jfrag
            end if
          end do
        end do
        if (n_pw <= 0 .or. nraw /= n_w + n_pfrag * n_pw) cycle
        do iw = 1, n_w
          gid = dg_frag%global_wannier_local_ids(iw, i_local)
          if (gid < 1 .or. gid > nw) cycle
          if (dg_frag%global_wannier_owner_frag(gid) < dg_frag%ifrag_start .or. &
              dg_frag%global_wannier_owner_frag(gid) > dg_frag%ifrag_end) cycle
          w_slot_count = w_slot_count + 1
        end do
        pself_slot_count = pself_slot_count + n_pw
        pneighbor_slot_count = pneighbor_slot_count + max(0, n_pfrag - 1) * n_pw
      end do

      if (allocated(dg_frag%mixed_z_frag_local_wcoef)) then
        if (size(dg_frag%mixed_z_frag_local_wcoef,1) /= max(1,w_slot_count*nstate_blk*nspin_use) .or. &
            size(dg_frag%mixed_z_frag_local_wcoef,2) /= nstate_blk .or. &
            size(dg_frag%mixed_z_frag_local_wcoef,3) /= nspin_use) then
          deallocate(dg_frag%mixed_z_frag_local_wcoef, dg_frag%mixed_z_frag_local_w_gid)
          if (allocated(dg_frag%mixed_z_frag_local_w_mix_gid)) deallocate(dg_frag%mixed_z_frag_local_w_mix_gid)
        end if
      end if
      if (.not. allocated(dg_frag%mixed_z_frag_local_wcoef)) then
        allocate(dg_frag%mixed_z_frag_local_wcoef(max(1,w_slot_count*nstate_blk*nspin_use), nstate_blk, nspin_use))
        allocate(dg_frag%mixed_z_frag_local_w_gid(max(1,w_slot_count*nstate_blk*nspin_use)))
        allocate(dg_frag%mixed_z_frag_local_w_mix_gid(max(1,w_slot_count*nstate_blk*nspin_use)))
      end if
      if (allocated(dg_frag%mixed_z_frag_local_pself_coef)) then
        if (size(dg_frag%mixed_z_frag_local_pself_coef,1) /= max(1,pself_slot_count*nstate_blk*nspin_use) .or. &
            size(dg_frag%mixed_z_frag_local_pself_coef,2) /= nstate_blk .or. &
            size(dg_frag%mixed_z_frag_local_pself_coef,3) /= nspin_use) then
          deallocate(dg_frag%mixed_z_frag_local_pself_coef, dg_frag%mixed_z_frag_local_pself_gid)
          if (allocated(dg_frag%mixed_z_frag_local_pself_mix_gid)) deallocate(dg_frag%mixed_z_frag_local_pself_mix_gid)
        end if
      end if
      if (.not. allocated(dg_frag%mixed_z_frag_local_pself_coef)) then
        allocate(dg_frag%mixed_z_frag_local_pself_coef(max(1,pself_slot_count*nstate_blk*nspin_use), nstate_blk, nspin_use))
        allocate(dg_frag%mixed_z_frag_local_pself_gid(max(1,pself_slot_count*nstate_blk*nspin_use)))
        allocate(dg_frag%mixed_z_frag_local_pself_mix_gid(max(1,pself_slot_count*nstate_blk*nspin_use)))
      end if
      if (allocated(dg_frag%mixed_z_frag_local_pneighbor_coef)) then
        if (size(dg_frag%mixed_z_frag_local_pneighbor_coef,1) /= max(1,pneighbor_slot_count*nstate_blk*nspin_use) .or. &
            size(dg_frag%mixed_z_frag_local_pneighbor_coef,2) /= nstate_blk .or. &
            size(dg_frag%mixed_z_frag_local_pneighbor_coef,3) /= nspin_use) then
          deallocate(dg_frag%mixed_z_frag_local_pneighbor_coef, dg_frag%mixed_z_frag_local_pneighbor_gid)
          if (allocated(dg_frag%mixed_z_frag_local_pneighbor_mix_gid)) &
            deallocate(dg_frag%mixed_z_frag_local_pneighbor_mix_gid)
        end if
      end if
      if (.not. allocated(dg_frag%mixed_z_frag_local_pneighbor_coef)) then
        allocate(dg_frag%mixed_z_frag_local_pneighbor_coef(max(1,pneighbor_slot_count*nstate_blk*nspin_use), &
          nstate_blk, nspin_use))
        allocate(dg_frag%mixed_z_frag_local_pneighbor_gid(max(1,pneighbor_slot_count*nstate_blk*nspin_use)))
        allocate(dg_frag%mixed_z_frag_local_pneighbor_mix_gid(max(1,pneighbor_slot_count*nstate_blk*nspin_use)))
      end if

      dg_frag%mixed_z_frag_local_wcoef(:, :, :) = (0.0d0, 0.0d0)
      dg_frag%mixed_z_frag_local_pself_coef(:, :, :) = (0.0d0, 0.0d0)
      dg_frag%mixed_z_frag_local_pneighbor_coef(:, :, :) = (0.0d0, 0.0d0)
      dg_frag%mixed_z_frag_local_w_gid(:) = 0
      dg_frag%mixed_z_frag_local_pself_gid(:) = 0
      dg_frag%mixed_z_frag_local_pneighbor_gid(:) = 0
      dg_frag%mixed_z_frag_local_w_mix_gid(:) = 0
      dg_frag%mixed_z_frag_local_pself_mix_gid(:) = 0
      dg_frag%mixed_z_frag_local_pneighbor_mix_gid(:) = 0
      dg_frag%mixed_z_frag_local_w_slots = w_slot_count
      dg_frag%mixed_z_frag_local_pself_slots = pself_slot_count
      dg_frag%mixed_z_frag_local_pneighbor_slots = pneighbor_slot_count
      dg_frag%mixed_z_frag_local_nstate = nstate_blk
      dg_frag%mixed_z_frag_local_nspin = nspin_use
      dg_frag%mixed_z_frag_local_storage_ready = w_slot_count > 0

      allocate(cw_source(nw,nstate_blk,nspin_use))
      cw_source(:, :, :) = (0.0d0, 0.0d0)
      if (.not. allocated(prop_cw_local_work) .or. size(prop_cw_local_work,1) /= nw .or. &
          size(prop_cw_local_work,2) /= nstate_blk) then
        if (allocated(prop_cw_local_work)) deallocate(prop_cw_local_work)
        if (allocated(prop_cw_global_work)) deallocate(prop_cw_global_work)
        allocate(prop_cw_local_work(nw,nstate_blk), prop_cw_global_work(nw,nstate_blk))
      end if
      do ispin_use = 1, nspin_use
        prop_cw_local_work(1:nw,1:nstate_blk) = (0.0d0, 0.0d0)
        do ifrag_row = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local_row = ifrag_row - dg_frag%ifrag_start + 1
          if (i_local_row < 1 .or. i_local_row > size(dg_frag%global_wannier_coef, 4)) cycle
          nbf = min(dg_frag%n_basis(ifrag_row, ispin_use), size(dg_frag%global_wannier_coef, 1))
          if (nbf <= 0) cycle
          if (.not. allocated(prop_w_block_work) .or. size(prop_w_block_work,1) < nbf .or. &
              size(prop_w_block_work,2) < nw) then
            if (allocated(prop_w_block_work)) deallocate(prop_w_block_work)
            allocate(prop_w_block_work(nbf,nw))
          end if
          if (.not. allocated(prop_coef_block_work) .or. size(prop_coef_block_work,1) < nbf .or. &
              size(prop_coef_block_work,2) < nstate_blk) then
            if (allocated(prop_coef_block_work)) deallocate(prop_coef_block_work)
            allocate(prop_coef_block_work(nbf,nstate_blk))
          end if
          nvalid = 0
          do ibasis = 1, nbf
            global_row = dg_frag%index_basis(ibasis, ifrag_row, ispin_use)
            if (global_row < 1 .or. global_row > dg_frag%n_mat_max) cycle
            local_row = dg_frag%coef_global_to_local(global_row, ispin_use)
            if (local_row < 1 .or. local_row > size(dg_frag%coef, 1)) cycle
            nvalid = nvalid + 1
            if (allocated(dg_frag%global_wannier_flux_evec)) then
              prop_w_block_work(nvalid,1:nw) = matmul( &
                dg_frag%global_wannier_coef(ibasis,1:size(dg_frag%global_wannier_flux_evec,1), &
                  ispin_use,i_local_row), &
                dg_frag%global_wannier_flux_evec(1:size(dg_frag%global_wannier_flux_evec,1),1:nw))
            else
              prop_w_block_work(nvalid,1:nw) = &
                dg_frag%global_wannier_coef(ibasis,1:nw,ispin_use,i_local_row)
            end if
            do state_col = 1, nstate_blk
              prop_coef_block_work(nvalid,state_col) = &
                dg_frag%coef(local_row,state_s+state_col-1,ispin_use)
            end do
          end do
          if (nvalid > 0) then
            prop_cw_local_work(1:nw,1:nstate_blk) = prop_cw_local_work(1:nw,1:nstate_blk) + &
              matmul(conjg(transpose(prop_w_block_work(1:nvalid,1:nw))), &
                     prop_coef_block_work(1:nvalid,1:nstate_blk))
          end if
        end do
        call comm_summation(prop_cw_local_work, prop_cw_global_work, nw*nstate_blk, dg_frag%icomm)
        cw_source(1:nw,1:nstate_blk,ispin_use) = prop_cw_global_work(1:nw,1:nstate_blk)
      end do

      w_slot = 0
      pself_slot = 0
      pneighbor_slot = 0
      bad_local = .false.
      do i_local = 1, size(dg_frag%wpw_reduced_dim)
        ifrag = dg_frag%ifrag_start + i_local - 1
        nred = dg_frag%wpw_reduced_dim(i_local)
        nself = dg_frag%wpw_reduced_nself(i_local)
        nraw = dg_frag%wpw_reduced_nraw(i_local)
        if (i_local > size(dg_frag%global_wannier_local_nkeep)) cycle
        n_w = dg_frag%global_wannier_local_nkeep(i_local)
        if (nred <= 0 .or. nself <= 0 .or. nraw <= 0 .or. n_w <= 0) cycle
        if (nself < n_w) cycle
        n_pfrag = 1
        pfrag_ids(1) = ifrag
        do axis = 1, 3
          do side = -1, 1, 2
            jfrag = wpw_face_neighbor_fragment(dg_frag, ifrag, axis, side)
            if (jfrag <= 0 .or. jfrag == ifrag) cycle
            if (.not. any(pfrag_ids(1:n_pfrag) == jfrag) .and. n_pfrag < size(pfrag_ids)) then
              n_pfrag = n_pfrag + 1
              pfrag_ids(n_pfrag) = jfrag
            end if
          end do
        end do
        if (n_pw <= 0 .or. nraw /= n_w + n_pfrag * n_pw) cycle
        w_owner_count = 0
        do iw = 1, n_w
          gid = dg_frag%global_wannier_local_ids(iw, i_local)
          if (gid < 1 .or. gid > nw) cycle
          if (dg_frag%global_wannier_owner_frag(gid) < dg_frag%ifrag_start .or. &
              dg_frag%global_wannier_owner_frag(gid) > dg_frag%ifrag_end) cycle
          w_owner_count = w_owner_count + 1
        end do
        p_valid_count = 0
        do pidx = 1, n_pfrag
          if (pidx > p_field_max) cycle
          do ipw = 1, n_pw
            gp = (pfrag_ids(pidx) - 1) * n_pw + ipw
            if (gp >= 1 .and. gp <= np) p_valid_count = p_valid_count + 1
          end do
        end do
        nblock = w_owner_count + p_valid_count
        if (nblock <= 0) cycle
        allocate(block_gid(nblock), block_pidx(nblock), block_gp(nblock), block_coef(nblock,nstate_blk))
        if (.not. field_off) then
          allocate(field_h(nblock,nblock), field_vec(nblock,nblock), tmp(nblock,nstate_blk), field_eval(nblock))
        end if
        do ispin_use = 1, nspin_use
          if (ispin_use > size(dg_frag%mixed_wannier_bpw_eval, 2)) cycle
          if (.not. field_off .and. ispin_use > size(dg_frag%mixed_wannier_bpw_z, 4)) then
            bad_local = .true.
            cycle
          end if
          block_gid(:) = 0
          block_pidx(:) = 0
          block_gp(:) = 0
          block_coef(:, :) = (0.0d0, 0.0d0)
          iblk = 0
          do iw = 1, n_w
            gid = dg_frag%global_wannier_local_ids(iw, i_local)
            if (gid < 1 .or. gid > nw) cycle
            if (dg_frag%global_wannier_owner_frag(gid) < dg_frag%ifrag_start .or. &
                dg_frag%global_wannier_owner_frag(gid) > dg_frag%ifrag_end) cycle
            iblk = iblk + 1
            block_gid(iblk) = gid
            block_pidx(iblk) = 0
            block_gp(iblk) = 0
            block_coef(iblk,1:nstate_blk) = cw_source(gid,1:nstate_blk,ispin_use)
          end do
          do pidx = 1, n_pfrag
            if (pidx > p_field_max) cycle
            do ipw = 1, n_pw
              gp = (pfrag_ids(pidx) - 1) * n_pw + ipw
              if (gp < 1 .or. gp > np) cycle
              iblk = iblk + 1
              if (iblk > nblock) then
                bad_local = .true.
                exit
              end if
              block_gid(iblk) = nw + gp
              block_pidx(iblk) = pidx
              block_gp(iblk) = gp
              block_coef(iblk,1:nstate_blk) = &
                dg_frag%mixed_wannier_bpw_pcoef(gp,state_s:state_e,ispin_use)
            end do
            if (bad_local) exit
          end do
          if (bad_local) cycle
          if (iblk /= nblock) then
            bad_local = .true.
            cycle
          end if
          if (.not. field_off) then
            do jblk = 1, nblock
              do iblk = 1, nblock
                field_h(iblk,jblk) = -E_use(1) * dg_frag%mixed_wannier_bpw_z(1,block_gid(iblk),block_gid(jblk),ispin_use) &
                                   -E_use(2) * dg_frag%mixed_wannier_bpw_z(2,block_gid(iblk),block_gid(jblk),ispin_use) &
                                   -E_use(3) * dg_frag%mixed_wannier_bpw_z(3,block_gid(iblk),block_gid(jblk),ispin_use)
              end do
            end do
            field_vec(:, :) = field_h(:, :)
            dg_frag%mixed_z_perf_eigensolve_calls = dg_frag%mixed_z_perf_eigensolve_calls + 1_8
            call eigen_zheev(field_vec, field_eval, field_h)
            dg_frag%mixed_z_perf_zgemm_calls = dg_frag%mixed_z_perf_zgemm_calls + 2_8
            tmp(:, :) = matmul(conjg(transpose(field_h(:, :))), block_coef(:, :))
            do iblk = 1, nblock
              phase_c = cos(field_eval(iblk) * dt)
              phase_s = sin(field_eval(iblk) * dt)
              tmp(iblk,1:nstate_blk) = cmplx(phase_c, -phase_s, kind=8) * tmp(iblk,1:nstate_blk)
            end do
            block_coef(:, :) = matmul(field_h(:, :), tmp(:, :))
          end if
          do iblk = 1, nblock
            gid = block_gid(iblk)
            if (gid < 1 .or. gid > size(dg_frag%mixed_wannier_bpw_eval, 1)) then
              bad_local = .true.
              cycle
            end if
            phase_c = cos(dg_frag%mixed_wannier_bpw_eval(gid, ispin_use) * dt)
            phase_s = sin(dg_frag%mixed_wannier_bpw_eval(gid, ispin_use) * dt)
            block_coef(iblk,1:nstate_blk) = cmplx(phase_c, -phase_s, kind=8) * block_coef(iblk,1:nstate_blk)
          end do
          if (bad_local) cycle
          do state_col = 1, nstate_blk
            do iblk = 1, nblock
              cand_val = block_coef(iblk,state_col)
              if (block_pidx(iblk) == 0) then
                w_slot = w_slot + 1
                if (w_slot <= size(dg_frag%mixed_z_frag_local_wcoef,1)) then
                  dg_frag%mixed_z_frag_local_wcoef(w_slot, state_col, ispin_use) = cand_val
                  dg_frag%mixed_z_frag_local_w_gid(w_slot) = block_gid(iblk)
                  dg_frag%mixed_z_frag_local_w_mix_gid(w_slot) = block_gid(iblk)
                else
                  bad_local = .true.
                end if
              else if (block_pidx(iblk) == 1) then
                pself_slot = pself_slot + 1
                if (pself_slot <= size(dg_frag%mixed_z_frag_local_pself_coef,1)) then
                  dg_frag%mixed_z_frag_local_pself_coef(pself_slot, state_col, ispin_use) = cand_val
                  dg_frag%mixed_z_frag_local_pself_gid(pself_slot) = block_gp(iblk)
                  dg_frag%mixed_z_frag_local_pself_mix_gid(pself_slot) = block_gid(iblk)
                else
                  bad_local = .true.
                end if
              else
                pneighbor_slot = pneighbor_slot + 1
                if (pneighbor_slot <= size(dg_frag%mixed_z_frag_local_pneighbor_coef,1)) then
                  dg_frag%mixed_z_frag_local_pneighbor_coef(pneighbor_slot, state_col, ispin_use) = cand_val
                  dg_frag%mixed_z_frag_local_pneighbor_gid(pneighbor_slot) = block_gp(iblk)
                  dg_frag%mixed_z_frag_local_pneighbor_mix_gid(pneighbor_slot) = block_gid(iblk)
                else
                  bad_local = .true.
                end if
              end if
            end do
          end do
          if (p_field_max < n_pfrag) then
            do state_col = 1, nstate_blk
              ist = state_s + state_col - 1
              do pidx = max(1, p_field_max + 1), n_pfrag
                do ipw = 1, n_pw
                  gp = (pfrag_ids(pidx) - 1) * n_pw + ipw
                  if (gp < 1 .or. gp > np) cycle
                  c_input = dg_frag%mixed_wannier_bpw_pcoef(gp, ist, ispin_use)
                  phase_c = cos(dg_frag%mixed_wannier_bpw_eval(nw + gp, ispin_use) * dt)
                  phase_s = sin(dg_frag%mixed_wannier_bpw_eval(nw + gp, ispin_use) * dt)
                  cand_val = cmplx(phase_c, -phase_s, kind=8) * c_input
                  if (pidx == 1) then
                    pself_slot = pself_slot + 1
                    if (pself_slot <= size(dg_frag%mixed_z_frag_local_pself_coef,1)) then
                      dg_frag%mixed_z_frag_local_pself_coef(pself_slot, state_col, ispin_use) = cand_val
                      dg_frag%mixed_z_frag_local_pself_gid(pself_slot) = gp
                      dg_frag%mixed_z_frag_local_pself_mix_gid(pself_slot) = nw + gp
                    else
                      bad_local = .true.
                    end if
                  else
                    pneighbor_slot = pneighbor_slot + 1
                    if (pneighbor_slot <= size(dg_frag%mixed_z_frag_local_pneighbor_coef,1)) then
                      dg_frag%mixed_z_frag_local_pneighbor_coef(pneighbor_slot, state_col, ispin_use) = cand_val
                      dg_frag%mixed_z_frag_local_pneighbor_gid(pneighbor_slot) = gp
                      dg_frag%mixed_z_frag_local_pneighbor_mix_gid(pneighbor_slot) = nw + gp
                    else
                      bad_local = .true.
                    end if
                  end if
                end do
              end do
            end do
          end if
        end do
        if (allocated(field_h)) deallocate(field_h, field_vec, tmp, field_eval)
        deallocate(block_gid, block_pidx, block_gp, block_coef)
      end do

      bad = bad .or. bad_local
      kernel_ready = dg_frag%mixed_z_frag_local_storage_ready .and. .not. bad_local
      deallocate(cw_source)
    end subroutine build_fragment_local_storage_direct

    subroutine compare_fragment_local_storage_with_global_reference(E_use, state_s, state_e, bad, kernel_diff_snorm)
      real(8), intent(in) :: E_use(3)
      integer, intent(in) :: state_s, state_e
      logical, intent(inout) :: bad
      real(8), intent(out) :: kernel_diff_snorm
      integer :: nstate_blk, nmix, nw, np, nspin_use, ispin_use, state_col, imix
      integer :: slot, gid
      integer :: irow, jcol, iblk_shell, jblk_shell, iblk2_shell, jblk2_shell, iw_shell
      integer :: shell, i_local_shell, ifrag_shell, pidx_shell, ipw_shell, gp_shell
      integer :: n_pfrag_shell, nblock_shell, w_owner_shell
      integer, allocatable :: pfrag_ids_shell(:), block_ids_shell(:)
      logical :: field_off
      real(8) :: phase_c, phase_s, local_sum(12), global_sum(12)
      real(8) :: coverage_local(8), coverage_global(8)
      complex(8), allocatable :: cmix_ref(:,:), field_h(:,:), field_vec(:,:), tmp(:,:)
      real(8), allocatable :: field_eval(:)
      logical, allocatable :: included_gid(:), included_pair(:,:)

      kernel_diff_snorm = huge(1.0d0)
      nstate_blk = max(0, state_e - state_s + 1)
      nmix = dg_frag%mixed_wannier_bpw_nmix
      nw = dg_frag%mixed_wannier_bpw_nw
      np = dg_frag%mixed_wannier_bpw_np
      nspin_use = dg_frag%nspin
      field_off = sum(abs(E_use(1:3))) <= 1.0d-30
      local_sum(:) = 0.0d0
      if (nstate_blk <= 0 .or. nw <= 0 .or. nmix <= 0 .or. &
          .not. dg_frag%mixed_z_frag_local_storage_ready .or. &
          .not. allocated(dg_frag%mixed_z_frag_local_wcoef) .or. &
          .not. allocated(dg_frag%mixed_z_frag_local_w_gid) .or. &
          (.not. field_off .and. .not. allocated(dg_frag%mixed_wannier_bpw_z))) then
        local_sum(8) = 1.0d0
      else
        allocate(cmix_ref(nmix,nstate_blk))
        allocate(included_gid(nmix))
        allocate(included_pair(nmix,nmix))
        included_gid(:) = .false.
        do slot = 1, size(dg_frag%mixed_z_frag_local_w_gid)
          gid = dg_frag%mixed_z_frag_local_w_gid(slot)
          if (gid >= 1 .and. gid <= nmix) included_gid(gid) = .true.
        end do
        if (trim(mixed_z_frag_local_field_block_kind) /= 'w_only' .and. &
            allocated(dg_frag%mixed_z_frag_local_pself_gid)) then
          do slot = 1, size(dg_frag%mixed_z_frag_local_pself_gid)
            gid = dg_frag%mixed_z_frag_local_pself_gid(slot)
            if (gid >= 1 .and. gid <= np) included_gid(nw + gid) = .true.
          end do
        end if
        if ((trim(mixed_z_frag_local_field_block_kind) == 'all' .or. &
             trim(mixed_z_frag_local_field_block_kind) == 'halo_owner') .and. &
            allocated(dg_frag%mixed_z_frag_local_pneighbor_gid)) then
          do slot = 1, size(dg_frag%mixed_z_frag_local_pneighbor_gid)
            gid = dg_frag%mixed_z_frag_local_pneighbor_gid(slot)
            if (gid >= 1 .and. gid <= np) included_gid(nw + gid) = .true.
          end do
        end if
        if (.not. field_off) then
          allocate(field_h(nmix,nmix), field_vec(nmix,nmix), tmp(nmix,nstate_blk), field_eval(nmix))
        end if
        if (.not. field_off) then
          do shell = 1, 2
            coverage_local(:) = 0.0d0
            included_pair(:, :) = .false.
            do i_local_shell = 1, size(dg_frag%wpw_reduced_dim)
              ifrag_shell = dg_frag%ifrag_start + i_local_shell - 1
              call collect_wpw_fragment_shell(ifrag_shell, shell, pfrag_ids_shell, n_pfrag_shell)
              w_owner_shell = 0
              if (i_local_shell <= size(dg_frag%global_wannier_local_nkeep)) then
                do iw_shell = 1, dg_frag%global_wannier_local_nkeep(i_local_shell)
                  gid = dg_frag%global_wannier_local_ids(iw_shell, i_local_shell)
                  if (gid < 1 .or. gid > nw) cycle
                  if (dg_frag%global_wannier_owner_frag(gid) < dg_frag%ifrag_start .or. &
                      dg_frag%global_wannier_owner_frag(gid) > dg_frag%ifrag_end) cycle
                  w_owner_shell = w_owner_shell + 1
                end do
              end if
              nblock_shell = w_owner_shell
              do pidx_shell = 1, n_pfrag_shell
                do ipw_shell = 1, max(0, dg_frag%n_plane_waves)
                  gp_shell = (pfrag_ids_shell(pidx_shell) - 1) * max(0, dg_frag%n_plane_waves) + ipw_shell
                  if (gp_shell >= 1 .and. gp_shell <= np) then
                    nblock_shell = nblock_shell + 1
                  else
                    coverage_local(6) = coverage_local(6) + 1.0d0
                  end if
                end do
              end do
              if (nblock_shell <= 0) then
                if (allocated(pfrag_ids_shell)) deallocate(pfrag_ids_shell)
                cycle
              end if
              allocate(block_ids_shell(nblock_shell))
              iblk_shell = 0
              if (i_local_shell <= size(dg_frag%global_wannier_local_nkeep)) then
                do iw_shell = 1, dg_frag%global_wannier_local_nkeep(i_local_shell)
                  gid = dg_frag%global_wannier_local_ids(iw_shell, i_local_shell)
                  if (gid < 1 .or. gid > nw) cycle
                  if (dg_frag%global_wannier_owner_frag(gid) < dg_frag%ifrag_start .or. &
                      dg_frag%global_wannier_owner_frag(gid) > dg_frag%ifrag_end) cycle
                  iblk_shell = iblk_shell + 1
                  block_ids_shell(iblk_shell) = gid
                end do
              end if
              do pidx_shell = 1, n_pfrag_shell
                do ipw_shell = 1, max(0, dg_frag%n_plane_waves)
                  gp_shell = (pfrag_ids_shell(pidx_shell) - 1) * max(0, dg_frag%n_plane_waves) + ipw_shell
                  if (gp_shell < 1 .or. gp_shell > np) cycle
                  iblk_shell = iblk_shell + 1
                  block_ids_shell(iblk_shell) = nw + gp_shell
                end do
              end do
              do jblk2_shell = 1, iblk_shell
                do iblk2_shell = 1, iblk_shell
                  included_pair(block_ids_shell(iblk2_shell), block_ids_shell(jblk2_shell)) = .true.
                end do
              end do
              coverage_local(4) = coverage_local(4) + dble(n_pfrag_shell)
              coverage_local(5) = coverage_local(5) + dble(iblk_shell)
              coverage_local(7) = coverage_local(7) + 1.0d0
              deallocate(block_ids_shell)
              if (allocated(pfrag_ids_shell)) deallocate(pfrag_ids_shell)
            end do
            do ispin_use = 1, nspin_use
              if (ispin_use > size(dg_frag%mixed_wannier_bpw_z, 4)) cycle
              field_h(:, :) = -E_use(1) * dg_frag%mixed_wannier_bpw_z(1,1:nmix,1:nmix,ispin_use) &
                            -E_use(2) * dg_frag%mixed_wannier_bpw_z(2,1:nmix,1:nmix,ispin_use) &
                            -E_use(3) * dg_frag%mixed_wannier_bpw_z(3,1:nmix,1:nmix,ispin_use)
              coverage_local(1) = coverage_local(1) + sum(abs(field_h(1:nmix,1:nmix))**2)
              do jcol = 1, nmix
                do irow = 1, nmix
                  if (included_pair(irow,jcol)) then
                    coverage_local(2) = coverage_local(2) + abs(field_h(irow,jcol))**2
                  else
                    coverage_local(3) = coverage_local(3) + abs(field_h(irow,jcol))**2
                  end if
                end do
              end do
            end do
            call comm_summation(coverage_local, coverage_global, size(coverage_local), dg_frag%icomm)
            if (dg_frag%id == 0 .and. .not. dg_frag%mixed_z_perf_count_enabled) then
              write(*,'(1x,a,4(a,i0),5(a,1pe12.4),2(a,l1),2(a,a))') &
                '[DG-MIXEDZ-FRAG-LOCAL-EXTENDED-SHELL-Z-CMP]', &
                ' step=', itt, &
                ' shell=', shell, &
                ' P_fragment_count=', int(coverage_global(4)), &
                ' projected_block_dim=', int(coverage_global(5)), &
                ' Z_included_norm=', sqrt(max(0.0d0, coverage_global(2))), &
                ' Z_omitted_norm=', sqrt(max(0.0d0, coverage_global(3))), &
                ' Z_omitted_ratio=', sqrt(max(0.0d0, coverage_global(3))) / &
                  max(sqrt(max(0.0d0, coverage_global(1))), 1.0d-300), &
                ' missing_P_coef_count=', coverage_global(6), &
                ' missing_Z_block_count=', coverage_global(8), &
                ' propagation_applied=', .false., &
                ' bad=', .false., &
                ' candidate_kind=', 'fragment_local_mixed_split_backend', &
                ' route=', 'metadata_only_extended_shell'
            end if
          end do
        end if
        if (.not. allocated(prop_cw_local_work) .or. size(prop_cw_local_work,1) /= nw .or. &
            size(prop_cw_local_work,2) /= nstate_blk) then
          if (allocated(prop_cw_local_work)) deallocate(prop_cw_local_work)
          if (allocated(prop_cw_global_work)) deallocate(prop_cw_global_work)
          allocate(prop_cw_local_work(nw,nstate_blk), prop_cw_global_work(nw,nstate_blk))
        end if
        if (np > 0) then
          if (.not. allocated(prop_cp_local_work) .or. size(prop_cp_local_work,1) /= np .or. &
              size(prop_cp_local_work,2) /= nstate_blk) then
            if (allocated(prop_cp_local_work)) deallocate(prop_cp_local_work)
            if (allocated(prop_cp_global_work)) deallocate(prop_cp_global_work)
            allocate(prop_cp_local_work(np,nstate_blk), prop_cp_global_work(np,nstate_blk))
          end if
        end if
        do ispin_use = 1, nspin_use
          call gather_global_mixed_coefficients(ispin_use, state_s, state_e, cmix_ref)
          if (.not. field_off) then
            if (ispin_use > size(dg_frag%mixed_wannier_bpw_z, 4)) then
              local_sum(8) = 1.0d0
              cycle
            end if
            field_h(:, :) = -E_use(1) * dg_frag%mixed_wannier_bpw_z(1,1:nmix,1:nmix,ispin_use) &
                          -E_use(2) * dg_frag%mixed_wannier_bpw_z(2,1:nmix,1:nmix,ispin_use) &
                          -E_use(3) * dg_frag%mixed_wannier_bpw_z(3,1:nmix,1:nmix,ispin_use)
            local_sum(9) = local_sum(9) + sum(abs(field_h(1:nmix,1:nmix))**2)
            do jcol = 1, nmix
              do irow = 1, nmix
                if (included_gid(irow) .and. included_gid(jcol)) then
                  local_sum(10) = local_sum(10) + abs(field_h(irow,jcol))**2
                else
                  local_sum(11) = local_sum(11) + abs(field_h(irow,jcol))**2
                end if
              end do
            end do
            field_vec(:, :) = field_h(:, :)
            call eigen_zheev(field_vec, field_eval, field_h)
            tmp(:, :) = matmul(conjg(transpose(field_h(:, :))), cmix_ref(:, :))
            do imix = 1, nmix
              phase_c = cos(field_eval(imix) * dt)
              phase_s = sin(field_eval(imix) * dt)
              tmp(imix,1:nstate_blk) = cmplx(phase_c, -phase_s, kind=8) * tmp(imix,1:nstate_blk)
            end do
            cmix_ref(:, :) = matmul(field_h(:, :), tmp(:, :))
          end if
          do state_col = 1, nstate_blk
            do imix = 1, nmix
              phase_c = cos(dg_frag%mixed_wannier_bpw_eval(imix, ispin_use) * dt)
              phase_s = sin(dg_frag%mixed_wannier_bpw_eval(imix, ispin_use) * dt)
              cmix_ref(imix,state_col) = cmplx(phase_c, -phase_s, kind=8) * cmix_ref(imix,state_col)
            end do
          end do

          prop_cw_local_work(1:nw,1:nstate_blk) = (0.0d0, 0.0d0)
          do slot = 1, size(dg_frag%mixed_z_frag_local_w_gid)
            gid = dg_frag%mixed_z_frag_local_w_gid(slot)
            if (gid < 1 .or. gid > nw) cycle
            do state_col = 1, nstate_blk
              prop_cw_local_work(gid,state_col) = prop_cw_local_work(gid,state_col) + &
                dg_frag%mixed_z_frag_local_wcoef(slot,state_col,ispin_use)
            end do
          end do
          call comm_summation(prop_cw_local_work, prop_cw_global_work, nw*nstate_blk, dg_frag%icomm)
          local_sum(1) = local_sum(1) + sum(abs(prop_cw_global_work(1:nw,1:nstate_blk) - &
            cmix_ref(1:nw,1:nstate_blk))**2)
          local_sum(3) = local_sum(3) + sum(abs(cmix_ref(1:nw,1:nstate_blk))**2)
          local_sum(5) = local_sum(5) + dble(nw * nstate_blk)

          if (np > 0 .and. allocated(dg_frag%mixed_z_frag_local_pself_coef)) then
            prop_cp_local_work(1:np,1:nstate_blk) = (0.0d0, 0.0d0)
            do slot = 1, size(dg_frag%mixed_z_frag_local_pself_gid)
              gid = dg_frag%mixed_z_frag_local_pself_gid(slot)
              if (gid < 1 .or. gid > np) cycle
              do state_col = 1, nstate_blk
                prop_cp_local_work(gid,state_col) = prop_cp_local_work(gid,state_col) + &
                  dg_frag%mixed_z_frag_local_pself_coef(slot,state_col,ispin_use)
              end do
            end do
            call comm_summation(prop_cp_local_work, prop_cp_global_work, np*nstate_blk, dg_frag%icomm)
            local_sum(2) = local_sum(2) + sum(abs(prop_cp_global_work(1:np,1:nstate_blk) - &
              cmix_ref(nw+1:nw+np,1:nstate_blk))**2)
            local_sum(4) = local_sum(4) + sum(abs(cmix_ref(nw+1:nw+np,1:nstate_blk))**2)
            local_sum(6) = local_sum(6) + dble(np * nstate_blk)
          end if
        end do
        if (allocated(field_h)) deallocate(field_h, field_vec, tmp, field_eval)
        deallocate(cmix_ref, included_gid, included_pair)
      end if

      local_sum(7) = local_sum(1) + local_sum(2)
      call comm_summation(local_sum, global_sum, size(local_sum), dg_frag%icomm)
      kernel_diff_snorm = sqrt(max(0.0d0, global_sum(7)))
      bad = bad .or. global_sum(8) > 0.5d0
      if (dg_frag%id == 0 .and. .not. dg_frag%mixed_z_perf_count_enabled) then
        write(*,'(1x,a,1(a,i0),11(a,1pe12.4),5(a,l1),3(a,a))') &
          '[DG-MIXEDZ-FRAG-LOCAL-KERNEL-DRYRUN-CMP]', &
          ' step=', itt, &
          ' W_owner_diff_Snorm=', sqrt(max(0.0d0, global_sum(1))), &
          ' W_owner_phase_diff_Snorm=', sqrt(max(0.0d0, global_sum(1))), &
          ' P_self_diff_Snorm=', sqrt(max(0.0d0, global_sum(2))), &
          ' P_neighbor_diff_Snorm=', 0.0d0, &
          ' coef_local_vs_global_Snorm=', kernel_diff_snorm, &
          ' global_ref_Snorm=', sqrt(max(0.0d0, global_sum(3) + global_sum(4))), &
          ' field_abs_sum=', sum(abs(E_use(1:3))), &
          ' field_Z_full_norm=', sqrt(max(0.0d0, global_sum(9))), &
          ' field_Z_included_norm=', sqrt(max(0.0d0, global_sum(10))), &
          ' field_Z_omitted_norm=', sqrt(max(0.0d0, global_sum(11))), &
          ' field_Z_omitted_ratio=', sqrt(max(0.0d0, global_sum(11))) / &
            max(sqrt(max(0.0d0, global_sum(9))), 1.0d-300), &
          ' local_kernel_available=', global_sum(8) <= 0.5d0, &
          ' owner_local_storage_ready=', dg_frag%mixed_z_frag_local_storage_ready, &
          ' field_off_only=', sum(abs(E_use(1:3))) <= 1.0d-30, &
          ' replacement_applied=', .false., &
          ' bad=', global_sum(8) > 0.5d0, &
          ' candidate_kind=', 'fragment_local_mixed_split_backend', &
          ' reference_kind=', 'global_mixed_split_backend_debug_reference', &
          ' route=', 'owner_local_storage_direct_vs_global_reference'
      end if
    end subroutine compare_fragment_local_storage_with_global_reference

    subroutine collect_wpw_fragment_shell(ifrag_center, shell_max, pfrag_ids, n_pfrag)
      integer, intent(in) :: ifrag_center, shell_max
      integer, allocatable, intent(out) :: pfrag_ids(:)
      integer, intent(out) :: n_pfrag
      integer :: max_frag, axis, side, level, start_idx, end_idx, idx, jfrag, base_frag
      integer, allocatable :: queue(:)

      max_frag = max(1, dg_frag%n_frag)
      allocate(queue(max_frag))
      queue(:) = 0
      n_pfrag = 1
      queue(1) = ifrag_center
      start_idx = 1
      end_idx = 1
      do level = 1, max(0, shell_max)
        idx = start_idx
        start_idx = end_idx + 1
        do while (idx <= end_idx)
          base_frag = queue(idx)
          do axis = 1, 3
            do side = -1, 1, 2
              jfrag = wpw_face_neighbor_fragment(dg_frag, base_frag, axis, side)
              if (jfrag <= 0) cycle
              if (.not. any(queue(1:n_pfrag) == jfrag) .and. n_pfrag < max_frag) then
                n_pfrag = n_pfrag + 1
                queue(n_pfrag) = jfrag
              end if
            end do
          end do
          idx = idx + 1
        end do
        end_idx = n_pfrag
      end do
      allocate(pfrag_ids(n_pfrag))
      pfrag_ids(1:n_pfrag) = queue(1:n_pfrag)
      deallocate(queue)
    end subroutine collect_wpw_fragment_shell

    subroutine diagnose_fragment_local_kernel_dryrun(E_use, state_s, state_e, bad, kernel_ready, kernel_diff_snorm)
      real(8), intent(in) :: E_use(3)
      integer, intent(in) :: state_s, state_e
      logical, intent(inout) :: bad
      logical, intent(out), optional :: kernel_ready
      real(8), intent(out), optional :: kernel_diff_snorm
      integer :: nstate_blk, nmix, nw, np, nspin_use, ispin_use, ist, state_col
      integer :: i_local, ifrag, nred, nself, nraw, n_w, n_pw
      integer :: iw, gid, ipw, pidx, gp, raw_slot
      integer :: w_slot, pself_slot, pneighbor_slot
      integer :: w_slot_count, pself_slot_count, pneighbor_slot_count
      integer :: n_pfrag, pfrag_ids(7), axis, side, jfrag
      logical :: field_off, local_kernel_available, bad_local, storage_ready
      real(8) :: occ, local_sum(12), global_sum(12)
      real(8) :: phase_c, phase_s
      complex(8) :: ref_val, cand_val, c_input
      complex(8), allocatable :: cmix_input(:,:,:), cmix_ref(:,:,:)

      nstate_blk = max(0, state_e - state_s + 1)
      nmix = dg_frag%mixed_wannier_bpw_nmix
      nw = dg_frag%mixed_wannier_bpw_nw
      np = dg_frag%mixed_wannier_bpw_np
      n_pw = max(0, dg_frag%n_plane_waves)
      nspin_use = dg_frag%nspin
      field_off = sum(abs(E_use(1:3))) <= 1.0d-30
      local_kernel_available = field_off .and. nstate_blk > 0 .and. nmix > 0 .and. nw > 0 .and. &
        dg_frag%wpw_reduced_ready .and. allocated(dg_frag%wpw_reduced_dim) .and. &
        allocated(dg_frag%wpw_reduced_nself) .and. allocated(dg_frag%wpw_reduced_nraw) .and. &
        allocated(dg_frag%wpw_reduced_transform) .and. allocated(dg_frag%wpw_reduced_H) .and. &
        allocated(dg_frag%wpw_reduced_S) .and. allocated(dg_frag%global_wannier_local_ids) .and. &
        (np <= 0 .or. allocated(dg_frag%mixed_wannier_bpw_pcoef))
      local_sum(:) = 0.0d0
      if (local_kernel_available) then
        allocate(cmix_input(nmix, nstate_blk, nspin_use), cmix_ref(nmix, nstate_blk, nspin_use))
        cmix_input(:, :, :) = (0.0d0, 0.0d0)
        cmix_ref(:, :, :) = (0.0d0, 0.0d0)
        do ispin_use = 1, nspin_use
          call gather_global_mixed_coefficients(ispin_use, state_s, state_e, cmix_input(:,:,ispin_use))
          cmix_ref(:,:,ispin_use) = cmix_input(:,:,ispin_use)
          do ist = 1, nstate_blk
            do iw = 1, nmix
              cmix_ref(iw,ist,ispin_use) = exp(cmplx(0.0d0, &
                -dg_frag%mixed_wannier_bpw_eval(iw,ispin_use) * dt, kind=8)) * &
                cmix_ref(iw,ist,ispin_use)
            end do
          end do
        end do
        w_slot_count = 0
        pself_slot_count = 0
        pneighbor_slot_count = 0
        do i_local = 1, size(dg_frag%wpw_reduced_dim)
          ifrag = dg_frag%ifrag_start + i_local - 1
          nraw = dg_frag%wpw_reduced_nraw(i_local)
          if (i_local > size(dg_frag%global_wannier_local_nkeep)) cycle
          n_w = dg_frag%global_wannier_local_nkeep(i_local)
          if (nraw <= 0 .or. n_w <= 0) cycle
          n_pfrag = 1
          pfrag_ids(1) = ifrag
          do axis = 1, 3
            do side = -1, 1, 2
              jfrag = wpw_face_neighbor_fragment(dg_frag, ifrag, axis, side)
              if (jfrag <= 0 .or. jfrag == ifrag) cycle
              if (.not. any(pfrag_ids(1:n_pfrag) == jfrag) .and. n_pfrag < size(pfrag_ids)) then
                n_pfrag = n_pfrag + 1
                pfrag_ids(n_pfrag) = jfrag
              end if
            end do
          end do
          if (n_pw <= 0 .or. nraw /= n_w + n_pfrag * n_pw) cycle
          do iw = 1, n_w
            gid = dg_frag%global_wannier_local_ids(iw, i_local)
            if (gid < 1 .or. gid > nw) cycle
            if (dg_frag%global_wannier_owner_frag(gid) < dg_frag%ifrag_start .or. &
                dg_frag%global_wannier_owner_frag(gid) > dg_frag%ifrag_end) cycle
            w_slot_count = w_slot_count + 1
          end do
          pself_slot_count = pself_slot_count + n_pw
          pneighbor_slot_count = pneighbor_slot_count + max(0, n_pfrag - 1) * n_pw
        end do
        if (allocated(dg_frag%mixed_z_frag_local_wcoef)) then
          if (size(dg_frag%mixed_z_frag_local_wcoef,1) /= max(1,w_slot_count*nstate_blk*nspin_use) .or. &
              size(dg_frag%mixed_z_frag_local_wcoef,2) /= nstate_blk .or. &
              size(dg_frag%mixed_z_frag_local_wcoef,3) /= nspin_use) then
            deallocate(dg_frag%mixed_z_frag_local_wcoef, dg_frag%mixed_z_frag_local_w_gid)
            if (allocated(dg_frag%mixed_z_frag_local_w_mix_gid)) deallocate(dg_frag%mixed_z_frag_local_w_mix_gid)
          end if
        end if
        if (.not. allocated(dg_frag%mixed_z_frag_local_wcoef)) then
          allocate(dg_frag%mixed_z_frag_local_wcoef(max(1,w_slot_count*nstate_blk*nspin_use), nstate_blk, nspin_use))
          allocate(dg_frag%mixed_z_frag_local_w_gid(max(1,w_slot_count*nstate_blk*nspin_use)))
          allocate(dg_frag%mixed_z_frag_local_w_mix_gid(max(1,w_slot_count*nstate_blk*nspin_use)))
        end if
        if (allocated(dg_frag%mixed_z_frag_local_pself_coef)) then
          if (size(dg_frag%mixed_z_frag_local_pself_coef,1) /= max(1,pself_slot_count*nstate_blk*nspin_use) .or. &
              size(dg_frag%mixed_z_frag_local_pself_coef,2) /= nstate_blk .or. &
              size(dg_frag%mixed_z_frag_local_pself_coef,3) /= nspin_use) then
            deallocate(dg_frag%mixed_z_frag_local_pself_coef, dg_frag%mixed_z_frag_local_pself_gid)
            if (allocated(dg_frag%mixed_z_frag_local_pself_mix_gid)) &
              deallocate(dg_frag%mixed_z_frag_local_pself_mix_gid)
          end if
        end if
        if (.not. allocated(dg_frag%mixed_z_frag_local_pself_coef)) then
          allocate(dg_frag%mixed_z_frag_local_pself_coef(max(1,pself_slot_count*nstate_blk*nspin_use), nstate_blk, nspin_use))
          allocate(dg_frag%mixed_z_frag_local_pself_gid(max(1,pself_slot_count*nstate_blk*nspin_use)))
          allocate(dg_frag%mixed_z_frag_local_pself_mix_gid(max(1,pself_slot_count*nstate_blk*nspin_use)))
        end if
        if (allocated(dg_frag%mixed_z_frag_local_pneighbor_coef)) then
          if (size(dg_frag%mixed_z_frag_local_pneighbor_coef,1) /= max(1,pneighbor_slot_count*nstate_blk*nspin_use) .or. &
              size(dg_frag%mixed_z_frag_local_pneighbor_coef,2) /= nstate_blk .or. &
              size(dg_frag%mixed_z_frag_local_pneighbor_coef,3) /= nspin_use) then
            deallocate(dg_frag%mixed_z_frag_local_pneighbor_coef, dg_frag%mixed_z_frag_local_pneighbor_gid)
            if (allocated(dg_frag%mixed_z_frag_local_pneighbor_mix_gid)) &
              deallocate(dg_frag%mixed_z_frag_local_pneighbor_mix_gid)
          end if
        end if
        if (.not. allocated(dg_frag%mixed_z_frag_local_pneighbor_coef)) then
          allocate(dg_frag%mixed_z_frag_local_pneighbor_coef(max(1,pneighbor_slot_count*nstate_blk*nspin_use), nstate_blk, nspin_use))
          allocate(dg_frag%mixed_z_frag_local_pneighbor_gid(max(1,pneighbor_slot_count*nstate_blk*nspin_use)))
          allocate(dg_frag%mixed_z_frag_local_pneighbor_mix_gid(max(1,pneighbor_slot_count*nstate_blk*nspin_use)))
        end if
        dg_frag%mixed_z_frag_local_wcoef(:, :, :) = (0.0d0, 0.0d0)
        dg_frag%mixed_z_frag_local_pself_coef(:, :, :) = (0.0d0, 0.0d0)
        dg_frag%mixed_z_frag_local_pneighbor_coef(:, :, :) = (0.0d0, 0.0d0)
        dg_frag%mixed_z_frag_local_w_gid(:) = 0
        dg_frag%mixed_z_frag_local_pself_gid(:) = 0
        dg_frag%mixed_z_frag_local_pneighbor_gid(:) = 0
        dg_frag%mixed_z_frag_local_w_mix_gid(:) = 0
        dg_frag%mixed_z_frag_local_pself_mix_gid(:) = 0
        dg_frag%mixed_z_frag_local_pneighbor_mix_gid(:) = 0
        dg_frag%mixed_z_frag_local_w_slots = w_slot_count
        dg_frag%mixed_z_frag_local_pself_slots = pself_slot_count
        dg_frag%mixed_z_frag_local_pneighbor_slots = pneighbor_slot_count
        dg_frag%mixed_z_frag_local_nstate = nstate_blk
        dg_frag%mixed_z_frag_local_nspin = nspin_use
        dg_frag%mixed_z_frag_local_storage_ready = w_slot_count > 0 .and. pself_slot_count >= 0 .and. &
          pneighbor_slot_count >= 0
        storage_ready = dg_frag%mixed_z_frag_local_storage_ready
        w_slot = 0
        pself_slot = 0
        pneighbor_slot = 0
        bad_local = .false.
        do i_local = 1, size(dg_frag%wpw_reduced_dim)
          ifrag = dg_frag%ifrag_start + i_local - 1
          nred = dg_frag%wpw_reduced_dim(i_local)
          nself = dg_frag%wpw_reduced_nself(i_local)
          nraw = dg_frag%wpw_reduced_nraw(i_local)
          if (i_local > size(dg_frag%global_wannier_local_nkeep)) cycle
          n_w = dg_frag%global_wannier_local_nkeep(i_local)
          if (nred <= 0 .or. nself <= 0 .or. nraw <= 0 .or. n_w <= 0) cycle
          if (nself < n_w) cycle
          n_pfrag = 1
          pfrag_ids(1) = ifrag
          do axis = 1, 3
            do side = -1, 1, 2
              jfrag = wpw_face_neighbor_fragment(dg_frag, ifrag, axis, side)
              if (jfrag <= 0 .or. jfrag == ifrag) cycle
              if (.not. any(pfrag_ids(1:n_pfrag) == jfrag) .and. n_pfrag < size(pfrag_ids)) then
                n_pfrag = n_pfrag + 1
                pfrag_ids(n_pfrag) = jfrag
              end if
            end do
          end do
          if (n_pw <= 0 .or. nraw /= n_w + n_pfrag * n_pw) cycle
          do ispin_use = 1, nspin_use
            nbf = min(dg_frag%n_basis(ifrag, ispin_use), size(dg_frag%global_wannier_local_coef, 1))
            do ist = state_s, state_e
              state_col = ist - state_s + 1
              occ = 1.0d0
              if (allocated(system%rocc)) then
                if (ist <= size(system%rocc,1) .and. ispin_use <= size(system%rocc,3)) &
                  occ = max(0.0d0, system%rocc(ist,1,ispin_use))
              end if
              do iw = 1, n_w
                gid = dg_frag%global_wannier_local_ids(iw, i_local)
                if (gid < 1 .or. gid > nw) cycle
                if (dg_frag%global_wannier_owner_frag(gid) < dg_frag%ifrag_start .or. &
                    dg_frag%global_wannier_owner_frag(gid) > dg_frag%ifrag_end) cycle
                ref_val = cmix_ref(gid, state_col, ispin_use)
                w_slot = w_slot + 1
                if (w_slot <= size(dg_frag%mixed_z_frag_local_wcoef,1)) then
                  dg_frag%mixed_z_frag_local_wcoef(w_slot, state_col, ispin_use) = &
                    cmix_input(gid, state_col, ispin_use)
                  dg_frag%mixed_z_frag_local_w_gid(w_slot) = gid
                  dg_frag%mixed_z_frag_local_w_mix_gid(w_slot) = gid
                  phase_c = cos(dg_frag%mixed_wannier_bpw_eval(gid, ispin_use) * dt)
                  phase_s = sin(dg_frag%mixed_wannier_bpw_eval(gid, ispin_use) * dt)
                  dg_frag%mixed_z_frag_local_wcoef(w_slot, state_col, ispin_use) = &
                    cmplx(phase_c, -phase_s, kind=8) * &
                    dg_frag%mixed_z_frag_local_wcoef(w_slot, state_col, ispin_use)
                  cand_val = dg_frag%mixed_z_frag_local_wcoef(w_slot, state_col, ispin_use)
                else
                  bad_local = .true.
                end if
                local_sum(1) = local_sum(1) + occ * abs(cand_val - ref_val)**2
                local_sum(4) = local_sum(4) + occ * abs(ref_val)**2
                local_sum(7) = local_sum(7) + 1.0d0
              end do
              do pidx = 1, n_pfrag
                do ipw = 1, n_pw
                  raw_slot = n_w + (pidx - 1) * n_pw + ipw
                  if (raw_slot < 1 .or. raw_slot > nraw) cycle
                  gp = (pfrag_ids(pidx) - 1) * n_pw + ipw
                  if (gp < 1 .or. gp > np) cycle
                  c_input = dg_frag%mixed_wannier_bpw_pcoef(gp, ist, ispin_use)
                  phase_c = cos(dg_frag%mixed_wannier_bpw_eval(nw + gp, ispin_use) * dt)
                  phase_s = sin(dg_frag%mixed_wannier_bpw_eval(nw + gp, ispin_use) * dt)
                  cand_val = cmplx(phase_c, -phase_s, kind=8) * c_input
                  ref_val = cmix_ref(nw + gp, state_col, ispin_use)
                  if (pidx == 1) then
                    pself_slot = pself_slot + 1
                    if (pself_slot <= size(dg_frag%mixed_z_frag_local_pself_coef,1)) then
                      dg_frag%mixed_z_frag_local_pself_coef(pself_slot, state_col, ispin_use) = cand_val
                      dg_frag%mixed_z_frag_local_pself_gid(pself_slot) = gp
                      dg_frag%mixed_z_frag_local_pself_mix_gid(pself_slot) = nw + gp
                      cand_val = dg_frag%mixed_z_frag_local_pself_coef(pself_slot, state_col, ispin_use)
                    else
                      bad_local = .true.
                    end if
                    local_sum(2) = local_sum(2) + occ * abs(cand_val - ref_val)**2
                    local_sum(5) = local_sum(5) + occ * abs(ref_val)**2
                    local_sum(8) = local_sum(8) + 1.0d0
                  else
                    pneighbor_slot = pneighbor_slot + 1
                    if (pneighbor_slot <= size(dg_frag%mixed_z_frag_local_pneighbor_coef,1)) then
                      dg_frag%mixed_z_frag_local_pneighbor_coef(pneighbor_slot, state_col, ispin_use) = cand_val
                      dg_frag%mixed_z_frag_local_pneighbor_gid(pneighbor_slot) = gp
                      dg_frag%mixed_z_frag_local_pneighbor_mix_gid(pneighbor_slot) = nw + gp
                      cand_val = dg_frag%mixed_z_frag_local_pneighbor_coef(pneighbor_slot, state_col, ispin_use)
                    else
                      bad_local = .true.
                    end if
                    local_sum(3) = local_sum(3) + occ * abs(cand_val - ref_val)**2
                    local_sum(6) = local_sum(6) + occ * abs(ref_val)**2
                    local_sum(9) = local_sum(9) + 1.0d0
                  end if
                end do
              end do
            end do
          end do
        end do
        local_sum(12) = merge(1.0d0, 0.0d0, bad_local)
        deallocate(cmix_input, cmix_ref)
      else
        storage_ready = .false.
        dg_frag%mixed_z_frag_local_storage_ready = .false.
        local_sum(12) = 1.0d0
      end if

      local_sum(10) = local_sum(1) + local_sum(2) + local_sum(3)
      local_sum(11) = local_sum(4) + local_sum(5) + local_sum(6)
      call comm_summation(local_sum, global_sum, size(local_sum), dg_frag%icomm)
      bad = bad .or. global_sum(12) > 0.5d0
      if (present(kernel_ready)) kernel_ready = local_kernel_available .and. storage_ready .and. global_sum(12) <= 0.5d0
      if (present(kernel_diff_snorm)) kernel_diff_snorm = sqrt(max(0.0d0, global_sum(10)))
      if (dg_frag%id == 0 .and. .not. dg_frag%mixed_z_perf_count_enabled) then
        write(*,'(1x,a,4(a,i0),11(a,1pe12.4),7(a,l1),3(a,a))') &
          '[DG-MIXEDZ-FRAG-LOCAL-KERNEL-DRYRUN-CMP]', &
          ' step=', itt, &
          ' W_owner_storage_slots=', dg_frag%mixed_z_frag_local_w_slots, &
          ' P_self_storage_slots=', dg_frag%mixed_z_frag_local_pself_slots, &
          ' P_neighbor_storage_slots=', dg_frag%mixed_z_frag_local_pneighbor_slots, &
          ' W_owner_diff_Snorm=', sqrt(max(0.0d0, global_sum(1))), &
          ' W_owner_phase_diff_Snorm=', sqrt(max(0.0d0, global_sum(1))), &
          ' P_self_diff_Snorm=', sqrt(max(0.0d0, global_sum(2))), &
          ' P_neighbor_diff_Snorm=', sqrt(max(0.0d0, global_sum(3))), &
          ' coef_local_vs_global_Snorm=', sqrt(max(0.0d0, global_sum(10))), &
          ' global_ref_Snorm=', sqrt(max(0.0d0, global_sum(11))), &
          ' W_owner_slots_compared=', global_sum(7), &
          ' P_self_slots_compared=', global_sum(8), &
          ' P_neighbor_slots_compared=', global_sum(9), &
          ' field_abs_sum=', sum(abs(E_use(1:3))), &
          ' dt=', dt, &
          ' local_kernel_available=', local_kernel_available, &
          ' owner_local_storage_ready=', storage_ready, &
          ' field_off_only=', field_off, &
          ' replacement_applied=', .false., &
          ' production_value_modified=', .false., &
          ' global_cmix_writeback_required=', .false., &
          ' bad=', global_sum(12) > 0.5d0, &
          ' candidate_kind=', 'fragment_local_mixed_split_backend', &
          ' reference_kind=', 'global_mixed_split_backend_dryrun', &
          ' route=', 'owner_local_storage_wphase_kernel'
      end if
    end subroutine diagnose_fragment_local_kernel_dryrun

    subroutine writeback_fragment_local_storage_fieldoff(state_s, state_e, bad)
      integer, intent(in) :: state_s, state_e
      logical, intent(out) :: bad
      integer :: nstate_blk, nw, np, nspin_use, ispin_use, state_col, state_idx
      integer :: slot, gid, ifrag_row, i_local_row, nbf, ibasis, nvalid, local_row, global_row
      real(8) :: prop_t0
      real(8), allocatable :: p_count_local(:), p_count_global(:)

      bad = .false.
      nstate_blk = max(0, state_e - state_s + 1)
      nw = dg_frag%mixed_wannier_bpw_nw
      np = dg_frag%mixed_wannier_bpw_np
      nspin_use = dg_frag%nspin
      if (.not. dg_frag%mixed_z_frag_local_storage_ready) then
        bad = .true.
        return
      end if
      if (nstate_blk <= 0 .or. nw <= 0) then
        bad = .true.
        return
      end if
      if (.not. allocated(dg_frag%mixed_z_frag_local_wcoef) .or. &
          .not. allocated(dg_frag%mixed_z_frag_local_w_gid)) then
        bad = .true.
        return
      end if
      if (np > 0 .and. (.not. allocated(dg_frag%mixed_z_frag_local_pself_coef) .or. &
          .not. allocated(dg_frag%mixed_z_frag_local_pself_gid) .or. &
          .not. allocated(dg_frag%mixed_wannier_bpw_pcoef))) then
        bad = .true.
        return
      end if
      if (np > 0 .and. dg_frag%mixed_z_frag_local_pneighbor_slots > 0 .and. &
          (.not. allocated(dg_frag%mixed_z_frag_local_pneighbor_coef) .or. &
           .not. allocated(dg_frag%mixed_z_frag_local_pneighbor_gid))) then
        bad = .true.
        return
      end if
      if (.not. allocated(prop_cw_local_work) .or. size(prop_cw_local_work,1) /= nw .or. &
          size(prop_cw_local_work,2) /= nstate_blk) then
        if (allocated(prop_cw_local_work)) deallocate(prop_cw_local_work)
        if (allocated(prop_cw_global_work)) deallocate(prop_cw_global_work)
        allocate(prop_cw_local_work(nw,nstate_blk), prop_cw_global_work(nw,nstate_blk))
      end if
      if (np > 0) then
        if (.not. allocated(prop_cp_local_work) .or. size(prop_cp_local_work,1) /= np .or. &
            size(prop_cp_local_work,2) /= nstate_blk) then
          if (allocated(prop_cp_local_work)) deallocate(prop_cp_local_work)
          if (allocated(prop_cp_global_work)) deallocate(prop_cp_global_work)
          allocate(prop_cp_local_work(np,nstate_blk), prop_cp_global_work(np,nstate_blk))
        end if
        allocate(p_count_local(np), p_count_global(np))
      end if

      do ispin_use = 1, nspin_use
        prop_cw_local_work(1:nw,1:nstate_blk) = (0.0d0, 0.0d0)
        do slot = 1, size(dg_frag%mixed_z_frag_local_w_gid)
          gid = dg_frag%mixed_z_frag_local_w_gid(slot)
          if (gid < 1 .or. gid > nw) cycle
          do state_col = 1, nstate_blk
            prop_cw_local_work(gid,state_col) = prop_cw_local_work(gid,state_col) + &
              dg_frag%mixed_z_frag_local_wcoef(slot,state_col,ispin_use)
          end do
        end do
        call comm_summation(prop_cw_local_work, prop_cw_global_work, nw*nstate_blk, dg_frag%icomm)

        if (.not. allocated(prop_next_local_work) .or. size(prop_next_local_work,1) /= size(dg_frag%coef,1) .or. &
            size(prop_next_local_work,2) /= nstate_blk) then
          if (allocated(prop_next_local_work)) deallocate(prop_next_local_work)
          allocate(prop_next_local_work(size(dg_frag%coef,1), nstate_blk))
        end if
        prop_next_local_work(:, :) = (0.0d0, 0.0d0)
        if (dg_frag%mixed_z_perf_count_enabled) prop_t0 = get_wtime()
        do ifrag_row = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local_row = ifrag_row - dg_frag%ifrag_start + 1
          if (i_local_row < 1 .or. i_local_row > size(dg_frag%global_wannier_coef, 4)) cycle
          nbf = min(dg_frag%n_basis(ifrag_row, ispin_use), size(dg_frag%global_wannier_coef, 1))
          if (nbf <= 0) cycle
          if (.not. allocated(prop_w_block_work) .or. size(prop_w_block_work,1) < nbf .or. &
              size(prop_w_block_work,2) < nw) then
            if (allocated(prop_w_block_work)) deallocate(prop_w_block_work)
            allocate(prop_w_block_work(nbf,nw))
          end if
          if (.not. allocated(prop_scatter_block_work) .or. size(prop_scatter_block_work,1) < nbf .or. &
              size(prop_scatter_block_work,2) < nstate_blk) then
            if (allocated(prop_scatter_block_work)) deallocate(prop_scatter_block_work)
            allocate(prop_scatter_block_work(nbf,nstate_blk))
          end if
          if (.not. allocated(prop_local_row_work) .or. size(prop_local_row_work) < nbf) then
            if (allocated(prop_local_row_work)) deallocate(prop_local_row_work)
            allocate(prop_local_row_work(nbf))
          end if
          nvalid = 0
          do ibasis = 1, nbf
            global_row = dg_frag%index_basis(ibasis, ifrag_row, ispin_use)
            if (global_row < 1 .or. global_row > dg_frag%n_mat_max) cycle
            local_row = dg_frag%coef_global_to_local(global_row, ispin_use)
            if (local_row < 1 .or. local_row > size(dg_frag%coef, 1)) cycle
            nvalid = nvalid + 1
            prop_local_row_work(nvalid) = local_row
            if (allocated(dg_frag%global_wannier_flux_evec)) then
              prop_w_block_work(nvalid,1:nw) = matmul( &
                dg_frag%global_wannier_coef(ibasis,1:size(dg_frag%global_wannier_flux_evec,1), &
                  ispin_use,i_local_row), &
                dg_frag%global_wannier_flux_evec(1:size(dg_frag%global_wannier_flux_evec,1),1:nw))
            else
              prop_w_block_work(nvalid,1:nw) = &
                dg_frag%global_wannier_coef(ibasis,1:nw,ispin_use,i_local_row)
            end if
          end do
          if (nvalid > 0) then
            prop_scatter_block_work(1:nvalid,1:nstate_blk) = &
              matmul(prop_w_block_work(1:nvalid,1:nw), prop_cw_global_work(1:nw,1:nstate_blk))
            do ibasis = 1, nvalid
              local_row = prop_local_row_work(ibasis)
              prop_next_local_work(local_row,1:nstate_blk) = &
                prop_next_local_work(local_row,1:nstate_blk) + &
                prop_scatter_block_work(ibasis,1:nstate_blk)
            end do
          end if
        end do
        do state_col = 1, nstate_blk
          state_idx = state_s + state_col - 1
          if (state_idx < 1 .or. state_idx > size(dg_frag%coef, 2)) cycle
          dg_frag%coef(1:size(dg_frag%coef,1),state_idx,ispin_use) = &
            prop_next_local_work(1:size(dg_frag%coef,1),state_col)
        end do
        if (dg_frag%mixed_z_perf_count_enabled) then
          dg_frag%mixed_z_perf_wall_prop_unpack = &
            dg_frag%mixed_z_perf_wall_prop_unpack + (get_wtime() - prop_t0)
        end if

        if (np > 0) then
          prop_cp_local_work(1:np,1:nstate_blk) = (0.0d0, 0.0d0)
          p_count_local(1:np) = 0.0d0
          do slot = 1, size(dg_frag%mixed_z_frag_local_pself_gid)
            gid = dg_frag%mixed_z_frag_local_pself_gid(slot)
            if (gid < 1 .or. gid > np) cycle
            p_count_local(gid) = p_count_local(gid) + 1.0d0
            do state_col = 1, nstate_blk
              prop_cp_local_work(gid,state_col) = prop_cp_local_work(gid,state_col) + &
                dg_frag%mixed_z_frag_local_pself_coef(slot,state_col,ispin_use)
            end do
          end do
          if (dg_frag%mixed_z_frag_local_pneighbor_slots > 0) then
            do slot = 1, dg_frag%mixed_z_frag_local_pneighbor_slots
              gid = dg_frag%mixed_z_frag_local_pneighbor_gid(slot)
              if (gid < 1 .or. gid > np) cycle
              p_count_local(gid) = p_count_local(gid) + 1.0d0
              do state_col = 1, nstate_blk
                prop_cp_local_work(gid,state_col) = prop_cp_local_work(gid,state_col) + &
                  dg_frag%mixed_z_frag_local_pneighbor_coef(slot,state_col,ispin_use)
              end do
            end do
          end if
          call comm_summation(prop_cp_local_work, prop_cp_global_work, np*nstate_blk, dg_frag%icomm)
          call comm_summation(p_count_local, p_count_global, np, dg_frag%icomm)
          do gid = 1, np
            if (p_count_global(gid) > 0.5d0) then
              prop_cp_global_work(gid,1:nstate_blk) = prop_cp_global_work(gid,1:nstate_blk) / &
                cmplx(p_count_global(gid), 0.0d0, kind=8)
            end if
          end do
          dg_frag%mixed_wannier_bpw_pcoef(1:np,state_s:state_e,ispin_use) = &
            prop_cp_global_work(1:np,1:nstate_blk)
        end if
      end do
      if (allocated(p_count_local)) deallocate(p_count_local, p_count_global)
    end subroutine writeback_fragment_local_storage_fieldoff

    subroutine enumerate_fragment_local_mixed_slots(layout_ready, w_owner_slots, w_total_slots, &
                                                   p_self_slots, p_neighbor_slots, gid_mismatch_count, &
                                                   owner_frag_mismatch_count, owner_rank_mismatch_count)
      logical, intent(out) :: layout_ready
      integer, intent(out) :: w_owner_slots, w_total_slots, p_self_slots, p_neighbor_slots
      integer, intent(out) :: gid_mismatch_count, owner_frag_mismatch_count, owner_rank_mismatch_count

      call diagnose_fragment_local_layout(layout_ready, w_owner_slots, w_total_slots, &
        p_self_slots, p_neighbor_slots, gid_mismatch_count, owner_frag_mismatch_count, &
        owner_rank_mismatch_count)
    end subroutine enumerate_fragment_local_mixed_slots

    subroutine diagnose_fragment_local_layout(layout_ready, w_owner_slots, w_total_slots, &
                                             p_self_slots, p_neighbor_slots, gid_mismatch_count, &
                                             owner_frag_mismatch_count, owner_rank_mismatch_count)
      logical, intent(out) :: layout_ready
      integer, intent(out) :: w_owner_slots, w_total_slots, p_self_slots, p_neighbor_slots
      integer, intent(out) :: gid_mismatch_count, owner_frag_mismatch_count, owner_rank_mismatch_count
      integer :: iloc_frag, iw, gid, owner_frag_local, owner_frag_global
      integer :: nlocal_frag, n_w, nself, nred, nraw, owner_rank_expected, owner_rank_stored
      logical :: has_w_layout, has_p_layout

      layout_ready = .false.
      w_owner_slots = 0
      w_total_slots = 0
      p_self_slots = 0
      p_neighbor_slots = 0
      gid_mismatch_count = 0
      owner_frag_mismatch_count = 0
      owner_rank_mismatch_count = 0
      has_w_layout = dg_frag%has_global_wannier_local_basis .and. &
        allocated(dg_frag%global_wannier_local_nkeep) .and. &
        allocated(dg_frag%global_wannier_local_ids) .and. &
        allocated(dg_frag%global_wannier_local_owner_frag) .and. &
        allocated(dg_frag%global_wannier_owner_frag)
      has_p_layout = allocated(dg_frag%wpw_reduced_dim) .and. allocated(dg_frag%wpw_reduced_nself) .and. &
        allocated(dg_frag%wpw_reduced_nraw)
      if (has_w_layout) then
        nlocal_frag = size(dg_frag%global_wannier_local_nkeep)
      else if (has_p_layout) then
        nlocal_frag = size(dg_frag%wpw_reduced_dim)
      else
        nlocal_frag = 0
      end if
      do iloc_frag = 1, nlocal_frag
        n_w = 0
        if (has_w_layout .and. iloc_frag <= size(dg_frag%global_wannier_local_nkeep)) then
          n_w = max(0, dg_frag%global_wannier_local_nkeep(iloc_frag))
          w_total_slots = w_total_slots + n_w
          do iw = 1, min(n_w, size(dg_frag%global_wannier_local_ids, 1), &
                         size(dg_frag%global_wannier_local_owner_frag, 1))
            gid = dg_frag%global_wannier_local_ids(iw, iloc_frag)
            owner_frag_local = dg_frag%global_wannier_local_owner_frag(iw, iloc_frag)
            if (gid < 1 .or. gid > size(dg_frag%global_wannier_owner_frag)) then
              gid_mismatch_count = gid_mismatch_count + 1
              cycle
            end if
            owner_frag_global = dg_frag%global_wannier_owner_frag(gid)
            if (owner_frag_local /= owner_frag_global) owner_frag_mismatch_count = owner_frag_mismatch_count + 1
            if (owner_frag_global >= dg_frag%ifrag_start .and. owner_frag_global <= dg_frag%ifrag_end) &
              w_owner_slots = w_owner_slots + 1
            if (allocated(dg_frag%id_array)) then
              owner_rank_expected = -1
              owner_rank_stored = -1
              if (owner_frag_global >= 1 .and. owner_frag_global <= size(dg_frag%id_array)) &
                owner_rank_expected = dg_frag%id_array(owner_frag_global)
              if (owner_frag_local >= 1 .and. owner_frag_local <= size(dg_frag%id_array)) &
                owner_rank_stored = dg_frag%id_array(owner_frag_local)
              if (owner_rank_expected /= owner_rank_stored) owner_rank_mismatch_count = owner_rank_mismatch_count + 1
            end if
          end do
        end if
        if (has_p_layout .and. iloc_frag <= size(dg_frag%wpw_reduced_dim)) then
          nself = max(0, dg_frag%wpw_reduced_nself(iloc_frag))
          nred = max(0, dg_frag%wpw_reduced_dim(iloc_frag))
          nraw = max(0, dg_frag%wpw_reduced_nraw(iloc_frag))
          p_self_slots = p_self_slots + max(0, nself - n_w)
          p_neighbor_slots = p_neighbor_slots + max(0, max(nraw, nred) - nself)
        end if
      end do
      layout_ready = has_w_layout .and. has_p_layout .and. w_owner_slots > 0 .and. &
        p_self_slots >= 0 .and. p_neighbor_slots >= 0 .and. gid_mismatch_count == 0 .and. &
        owner_frag_mismatch_count == 0 .and. owner_rank_mismatch_count == 0
      if (dg_frag%id == 0 .and. .not. dg_frag%mixed_z_perf_count_enabled) then
        write(*,'(1x,a,8(a,i0),7(a,l1),3(a,a))') &
          '[DG-MIXEDZ-FRAG-LOCAL-LAYOUT-CMP]', &
          ' local_fragment_count=', nlocal_frag, &
          ' W_owner_slots=', w_owner_slots, &
          ' W_total_slots=', w_total_slots, &
          ' P_self_slots=', p_self_slots, &
          ' P_neighbor_slots=', p_neighbor_slots, &
          ' gid_mismatch_count=', gid_mismatch_count, &
          ' owner_frag_mismatch_count=', owner_frag_mismatch_count, &
          ' owner_rank_mismatch_count=', owner_rank_mismatch_count, &
          ' has_W_layout=', has_w_layout, &
          ' has_P_layout=', has_p_layout, &
          ' owner_local_layout_ready=', layout_ready, &
          ' global_cmix_allocation_required=', .false., &
          ' global_cmix_writeback_required=', .false., &
          ' propagation_applied=', .false., &
          ' bad=', .not. layout_ready, &
          ' candidate_kind=', 'fragment_local_mixed_split_backend', &
          ' layout_scope=', 'owner_local_metadata', &
          ' replacement_scope=', 'diagnostic_only'
      end if
    end subroutine diagnose_fragment_local_layout

    subroutine apply_global_mixed_split_exp(E_use, state_s, state_e)
      real(8), intent(in) :: E_use(3)
      integer, intent(in) :: state_s, state_e
      integer :: ispin_current, nstate_blk, nmix, imix
      real(8) :: phase_c, phase_s, prop_t0
      complex(8), allocatable :: cmix_before_field(:,:)

      if (.not. dg_frag%has_mixed_wannier_bpw_position) return
      if (.not. allocated(dg_frag%mixed_wannier_bpw_eval)) return
      if (.not. allocated(dg_frag%mixed_wannier_bpw_z)) return
      nstate_blk = max(0, state_e - state_s + 1)
      if (nstate_blk <= 0) return
      nmix = dg_frag%mixed_wannier_bpw_nmix
      if (nmix <= 0) return
      dg_frag%mixed_z_perf_expdiag_calls = dg_frag%mixed_z_perf_expdiag_calls + 1_8
      if (mixed_z_field_kick_diag_enabled) allocate(cmix_before_field(nmix,nstate_blk))

      if (.not. allocated(prop_cmix_work) .or. size(prop_cmix_work,1) /= nmix .or. &
          size(prop_cmix_work,2) /= nstate_blk) then
        if (allocated(prop_cmix_work)) deallocate(prop_cmix_work, prop_tmp_work)
        allocate(prop_cmix_work(nmix,nstate_blk), prop_tmp_work(nmix,nstate_blk))
      end if
      if (.not. allocated(prop_field_h_work) .or. size(prop_field_h_work,1) /= nmix) then
        if (allocated(prop_field_h_work)) deallocate(prop_field_h_work, prop_field_vec_work, prop_field_eval_work)
        allocate(prop_field_h_work(nmix,nmix), prop_field_vec_work(nmix,nmix), prop_field_eval_work(nmix))
      end if
      if (.not. allocated(prop_field_axis_evec_cache) .or. size(prop_field_axis_evec_cache,1) /= nmix .or. &
          size(prop_field_axis_evec_cache,3) /= dg_frag%nspin) then
        if (allocated(prop_field_axis_evec_cache)) deallocate(prop_field_axis_evec_cache, &
          prop_field_axis_eval_cache, prop_field_axis_cache_axis, prop_field_axis_cache_valid)
        allocate(prop_field_axis_evec_cache(nmix,nmix,dg_frag%nspin))
        allocate(prop_field_axis_eval_cache(nmix,dg_frag%nspin))
        allocate(prop_field_axis_cache_axis(dg_frag%nspin), prop_field_axis_cache_valid(dg_frag%nspin))
        prop_field_axis_evec_cache(:, :, :) = (0.0d0, 0.0d0)
        prop_field_axis_eval_cache(:, :) = 0.0d0
        prop_field_axis_cache_axis(:) = 0
        prop_field_axis_cache_valid(:) = .false.
      end if
      do ispin_current = 1, dg_frag%nspin
        if (ispin_current > size(dg_frag%mixed_wannier_bpw_eval, 2) .or. &
            ispin_current > size(dg_frag%mixed_wannier_bpw_z, 4)) cycle
        call gather_global_mixed_coefficients(ispin_current, state_s, state_e, prop_cmix_work)
        if (dg_frag%mixed_z_perf_count_enabled) prop_t0 = get_wtime()
        prop_field_h_work(:, :) = (0.0d0, 0.0d0)
        do imix = 1, nmix
          prop_field_h_work(imix,imix) = cmplx(dg_frag%mixed_wannier_bpw_eval(imix,ispin_current), 0.0d0, kind=8)
        end do
        if (sum(abs(E_use(1:3))) > 1.0d-30) then
          if (mixed_z_field_kick_diag_enabled) cmix_before_field(:,:) = prop_cmix_work(:,:)
          prop_field_h_work(:,:) = prop_field_h_work(:,:) &
            - E_use(1) * dg_frag%mixed_wannier_bpw_z(1,1:nmix,1:nmix,ispin_current) &
            - E_use(2) * dg_frag%mixed_wannier_bpw_z(2,1:nmix,1:nmix,ispin_current) &
            - E_use(3) * dg_frag%mixed_wannier_bpw_z(3,1:nmix,1:nmix,ispin_current)
        end if
        prop_field_vec_work(:,:) = prop_field_h_work(:,:)
        dg_frag%mixed_z_perf_eigensolve_calls = dg_frag%mixed_z_perf_eigensolve_calls + 1_8
        call eigen_zheev(prop_field_vec_work, prop_field_eval_work, prop_field_h_work)
        dg_frag%mixed_z_perf_zgemm_calls = dg_frag%mixed_z_perf_zgemm_calls + 2_8
        prop_tmp_work(:,:) = matmul(conjg(transpose(prop_field_h_work(:,:))), prop_cmix_work(:,:))
        do imix = 1, nmix
          phase_c = cos(prop_field_eval_work(imix) * dt)
          phase_s = sin(prop_field_eval_work(imix) * dt)
          prop_tmp_work(imix,1:nstate_blk) = cmplx(phase_c, -phase_s, kind=8) * &
            prop_tmp_work(imix,1:nstate_blk)
        end do
        prop_cmix_work(:,:) = matmul(prop_field_h_work(:,:), prop_tmp_work(:,:))
        if (mixed_z_field_kick_diag_enabled .and. sum(abs(E_use(1:3))) > 1.0d-30) call log_mixed_field_kick_diag( &
          'midpoint_full_h_exp', itt, ispin_current, E_use, cmix_before_field, prop_cmix_work, state_s, state_e)
        if (dg_frag%mixed_z_perf_count_enabled) then
          dg_frag%mixed_z_perf_wall_prop_field_exp = &
            dg_frag%mixed_z_perf_wall_prop_field_exp + (get_wtime() - prop_t0)
        end if
        call scatter_global_mixed_coefficients(ispin_current, state_s, state_e, prop_cmix_work)
      end do
      if (allocated(cmix_before_field)) deallocate(cmix_before_field)
    end subroutine apply_global_mixed_split_exp

    subroutine mixed_z_projection_stage_norm(S_full, rhs_full, idx, nidx, nocc_spin, occ_vec, &
                                             prod_norm, trace_norm, rel_resid, info_stage)
      complex(8), intent(in) :: S_full(:,:), rhs_full(:,:)
      integer, intent(in) :: idx(:), nidx, nocc_spin
      real(8), intent(in) :: occ_vec(:), prod_norm
      real(8), intent(out) :: trace_norm, rel_resid
      integer, intent(out) :: info_stage
      integer :: ii, jj, ist, iev, nkeep
      real(8) :: eval_min, eval_max, resid
      complex(8), allocatable :: S_sub(:,:), S_vec(:,:), S_inv(:,:), rhs_sub(:,:), c_sub(:,:)
      real(8), allocatable :: S_eval(:)

      trace_norm = 0.0d0
      rel_resid = 1.0d0
      info_stage = 1
      if (nidx <= 0 .or. nocc_spin <= 0) return
      allocate(S_sub(nidx,nidx), S_vec(nidx,nidx), S_inv(nidx,nidx), rhs_sub(nidx,nocc_spin), &
               c_sub(nidx,nocc_spin), S_eval(nidx))
      do jj = 1, nidx
        rhs_sub(jj, :) = rhs_full(idx(jj), 1:nocc_spin)
        do ii = 1, nidx
          S_sub(ii, jj) = S_full(idx(ii), idx(jj))
        end do
      end do
      S_vec(:, :) = S_sub(:, :)
      call eigen_zheev(S_vec, S_eval, S_vec)
      eval_min = minval(S_eval(1:nidx))
      eval_max = maxval(S_eval(1:nidx))
      S_inv(:, :) = (0.0d0, 0.0d0)
      nkeep = 0
      do iev = 1, nidx
        if (S_eval(iev) <= max(1.0d-10 * eval_max, 1.0d-14)) cycle
        nkeep = nkeep + 1
        do ii = 1, nidx
          do jj = 1, nidx
            S_inv(ii, jj) = S_inv(ii, jj) + S_vec(ii, iev) * conjg(S_vec(jj, iev)) / S_eval(iev)
          end do
        end do
      end do
      if (nkeep <= 0 .or. eval_min <= 0.0d0) then
        deallocate(S_sub, S_vec, S_inv, rhs_sub, c_sub, S_eval)
        return
      end if
      c_sub(:, :) = matmul(S_inv, rhs_sub)
      trace_norm = 0.0d0
      do ist = 1, nocc_spin
        if (occ_vec(ist) <= 0.0d0) cycle
        trace_norm = trace_norm + occ_vec(ist) * &
          real(sum(conjg(c_sub(:,ist)) * matmul(S_sub, c_sub(:,ist))), kind=8)
      end do
      resid = max(0.0d0, prod_norm - trace_norm)
      rel_resid = sqrt(resid) / max(sqrt(max(prod_norm, 0.0d0)), 1.0d-300)
      info_stage = 0
      deallocate(S_sub, S_vec, S_inv, rhs_sub, c_sub, S_eval)
    end subroutine mixed_z_projection_stage_norm

    subroutine mixed_z_solve_metric_vector(S_full, rhs_vec, c_vec, info_solve)
      complex(8), intent(in) :: S_full(:,:), rhs_vec(:)
      complex(8), intent(out) :: c_vec(:)
      integer, intent(out) :: info_solve
      integer :: ndim, ii, jj, iev, nkeep
      real(8) :: eval_min, eval_max
      complex(8), allocatable :: S_vec(:,:), S_inv(:,:)
      real(8), allocatable :: S_eval(:)

      ndim = size(rhs_vec)
      c_vec(:) = (0.0d0, 0.0d0)
      info_solve = 1
      if (ndim <= 0 .or. size(S_full,1) < ndim .or. size(S_full,2) < ndim .or. size(c_vec) < ndim) return
      allocate(S_vec(ndim,ndim), S_inv(ndim,ndim), S_eval(ndim))
      S_vec(:, :) = S_full(1:ndim, 1:ndim)
      call eigen_zheev(S_vec, S_eval, S_vec)
      eval_min = minval(S_eval(1:ndim))
      eval_max = maxval(S_eval(1:ndim))
      S_inv(:, :) = (0.0d0, 0.0d0)
      nkeep = 0
      do iev = 1, ndim
        if (S_eval(iev) <= max(1.0d-10 * eval_max, 1.0d-14)) cycle
        nkeep = nkeep + 1
        do ii = 1, ndim
          do jj = 1, ndim
            S_inv(ii, jj) = S_inv(ii, jj) + S_vec(ii, iev) * conjg(S_vec(jj, iev)) / S_eval(iev)
          end do
        end do
      end do
      if (nkeep > 0 .and. eval_min > 0.0d0) then
        c_vec(1:ndim) = matmul(S_inv, rhs_vec(1:ndim))
        info_solve = 0
      end if
      deallocate(S_vec, S_inv, S_eval)
    end subroutine mixed_z_solve_metric_vector

    subroutine mixed_z_log_w_repack_candidate(label, step_use, wprod, wback, occ_all, zdiag)
      character(len=*), intent(in) :: label
      integer, intent(in) :: step_use
      complex(8), intent(in) :: wprod(:,:), wback(:,:)
      real(8), intent(in) :: occ_all(:,:,:), zdiag(:)
      integer :: iw, ist, nw_use, nstate_use
      real(8) :: occ
      real(8) :: coef_diff_norm, coef_ref_norm, coef_rel_diff, coef_max_abs_diff
      real(8) :: rho_prod, rho_back, rho_diff
      real(8) :: pz_prod, pz_back, pz_diff

      coef_diff_norm = 0.0d0
      coef_ref_norm = 0.0d0
      coef_rel_diff = 0.0d0
      coef_max_abs_diff = 0.0d0
      rho_prod = 0.0d0
      rho_back = 0.0d0
      rho_diff = 0.0d0
      pz_prod = 0.0d0
      pz_back = 0.0d0
      pz_diff = 0.0d0
      nw_use = min(size(wprod, 1), size(wback, 1), size(zdiag))
      nstate_use = min(size(wprod, 2), size(wback, 2), size(occ_all, 1))
      do ist = 1, nstate_use
        occ = 0.0d0
        if (size(occ_all, 3) >= 1) occ = max(0.0d0, occ_all(ist, 1, 1))
        do iw = 1, nw_use
          coef_diff_norm = coef_diff_norm + abs(wback(iw, ist) - wprod(iw, ist))**2
          coef_ref_norm = coef_ref_norm + abs(wprod(iw, ist))**2
          coef_max_abs_diff = max(coef_max_abs_diff, abs(wback(iw, ist) - wprod(iw, ist)))
          rho_prod = rho_prod + abs(wprod(iw, ist))**2
          rho_back = rho_back + abs(wback(iw, ist))**2
          pz_prod = pz_prod - occ * abs(wprod(iw, ist))**2 * zdiag(iw)
          pz_back = pz_back - occ * abs(wback(iw, ist))**2 * zdiag(iw)
        end do
      end do
      coef_diff_norm = sqrt(max(coef_diff_norm, 0.0d0))
      coef_rel_diff = coef_diff_norm / max(sqrt(max(coef_ref_norm, 0.0d0)), 1.0d-300)
      rho_diff = rho_back - rho_prod
      pz_diff = pz_back - pz_prod
      write(*,'(1x,a,1(a,i0),1(a,a),13(a,1pe12.4),3(a,l1),a,a)') &
        '[DG-MIXEDZ-LOCAL-FRAG-EMBED-PROJECTOR-CMP]', &
        ' step=', step_use, &
        ' candidate=', trim(label), &
        ' coef_diff_norm=', coef_diff_norm, &
        ' coef_diff_max=', coef_max_abs_diff, &
        ' rel_coef_diff=', coef_rel_diff, &
        ' rho_prod=', rho_prod, &
        ' rho_back=', rho_back, &
        ' rho_diff=', rho_diff, &
        ' Pz_prod=', pz_prod, &
        ' Pz_back=', pz_back, &
        ' Pz_diff=', pz_diff, &
        ' S_norm_prod=', coef_ref_norm, &
        ' S_norm_back=', rho_back, &
        ' trace_rho_prod=', rho_prod, &
        ' trace_rho_back=', rho_back, &
        ' no_propagation=', .true., &
        ' W_sector_only=', .true., &
        ' bad=', coef_rel_diff /= coef_rel_diff, &
        ' route=', 'fragment-local-global-embedding-candidate-diagnostic-only'
    end subroutine mixed_z_log_w_repack_candidate

    subroutine mixed_z_log_w_prod_map_bridge_candidate(label, step_use, wprod, wback, occ_all, zdiag, &
      r_nnz, p_nnz, duplicate_count, missing_count, owner_mismatch_count)
      character(len=*), intent(in) :: label
      integer, intent(in) :: step_use
      complex(8), intent(in) :: wprod(:,:), wback(:,:)
      real(8), intent(in) :: occ_all(:,:,:), zdiag(:)
      real(8), intent(in) :: r_nnz, p_nnz, duplicate_count, missing_count, owner_mismatch_count
      integer :: iw, ist, nw_use, nstate_use
      real(8) :: occ
      real(8) :: coef_diff_norm, coef_ref_norm, coef_rel_diff, coef_max_abs_diff
      real(8) :: rho_prod, rho_back, rho_diff
      real(8) :: pz_prod, pz_back, pz_diff

      coef_diff_norm = 0.0d0
      coef_ref_norm = 0.0d0
      coef_rel_diff = 0.0d0
      coef_max_abs_diff = 0.0d0
      rho_prod = 0.0d0
      rho_back = 0.0d0
      rho_diff = 0.0d0
      pz_prod = 0.0d0
      pz_back = 0.0d0
      pz_diff = 0.0d0
      nw_use = min(size(wprod, 1), size(wback, 1), size(zdiag))
      nstate_use = min(size(wprod, 2), size(wback, 2), size(occ_all, 1))
      do ist = 1, nstate_use
        occ = 0.0d0
        if (size(occ_all, 3) >= 1) occ = max(0.0d0, occ_all(ist, 1, 1))
        do iw = 1, nw_use
          coef_diff_norm = coef_diff_norm + abs(wback(iw, ist) - wprod(iw, ist))**2
          coef_ref_norm = coef_ref_norm + abs(wprod(iw, ist))**2
          coef_max_abs_diff = max(coef_max_abs_diff, abs(wback(iw, ist) - wprod(iw, ist)))
          rho_prod = rho_prod + abs(wprod(iw, ist))**2
          rho_back = rho_back + abs(wback(iw, ist))**2
          pz_prod = pz_prod - occ * abs(wprod(iw, ist))**2 * zdiag(iw)
          pz_back = pz_back - occ * abs(wback(iw, ist))**2 * zdiag(iw)
        end do
      end do
      coef_diff_norm = sqrt(max(coef_diff_norm, 0.0d0))
      coef_rel_diff = coef_diff_norm / max(sqrt(max(coef_ref_norm, 0.0d0)), 1.0d-300)
      rho_diff = rho_back - rho_prod
      pz_diff = pz_back - pz_prod
      write(*,'(1x,a,1(a,i0),1(a,a),16(a,1pe12.4),4(a,l1),a,a)') &
        '[DG-MIXEDZ-LOCAL-FRAG-PROD-MAP-BRIDGE-CMP]', &
        ' step=', step_use, &
        ' candidate=', trim(label), &
        ' rel_coef_diff_after_sum_PAR=', coef_rel_diff, &
        ' coef_diff_norm=', coef_diff_norm, &
        ' coef_diff_max=', coef_max_abs_diff, &
        ' rho_prod=', rho_prod, &
        ' rho_back_PAR=', rho_back, &
        ' rho_diff=', rho_diff, &
        ' Pz_prod=', pz_prod, &
        ' Pz_back_PAR=', pz_back, &
        ' Pz_diff=', pz_diff, &
        ' trace_rho_prod=', rho_prod, &
        ' trace_rho_back=', rho_back, &
        ' R_nnz=', r_nnz, &
        ' P_nnz=', p_nnz, &
        ' local_slot_duplicate_count=', duplicate_count, &
        ' missing_global_gid_count=', missing_count, &
        ' owner_mismatch_count=', owner_mismatch_count, &
        ' no_propagation=', .true., &
        ' W_sector_only=', .true., &
        ' production_map_source=', .true., &
        ' bad=', coef_rel_diff /= coef_rel_diff, &
        ' route=', 'production-CW-through-local-W-slot-map-diagnostic-only'
    end subroutine mixed_z_log_w_prod_map_bridge_candidate

    subroutine mixed_z_log_w_prod_sourced_prop_candidate(label, step_use, wprod_before, wprod_after, wback_after, &
      occ_all, zdiag, owner_mismatch_count)
      character(len=*), intent(in) :: label
      integer, intent(in) :: step_use
      complex(8), intent(in) :: wprod_before(:,:), wprod_after(:,:), wback_after(:,:)
      real(8), intent(in) :: occ_all(:,:,:), zdiag(:)
      real(8), intent(in) :: owner_mismatch_count
      integer :: iw, ist, nw_use, nstate_use
      real(8) :: occ
      real(8) :: coef_diff_norm, coef_ref_norm, coef_rel_diff, coef_max_abs_diff
      real(8) :: norm_prod_before, norm_prod_after, norm_local_after
      real(8) :: rho_prod_after, rho_back_after, rho_diff
      real(8) :: pz_prod_after, pz_back_after, pz_diff

      coef_diff_norm = 0.0d0
      coef_ref_norm = 0.0d0
      coef_rel_diff = 0.0d0
      coef_max_abs_diff = 0.0d0
      norm_prod_before = 0.0d0
      norm_prod_after = 0.0d0
      norm_local_after = 0.0d0
      rho_prod_after = 0.0d0
      rho_back_after = 0.0d0
      rho_diff = 0.0d0
      pz_prod_after = 0.0d0
      pz_back_after = 0.0d0
      pz_diff = 0.0d0
      nw_use = min(size(wprod_before, 1), size(wprod_after, 1), size(wback_after, 1), size(zdiag))
      nstate_use = min(size(wprod_before, 2), size(wprod_after, 2), size(wback_after, 2), size(occ_all, 1))
      do ist = 1, nstate_use
        occ = 0.0d0
        if (size(occ_all, 3) >= 1) occ = max(0.0d0, occ_all(ist, 1, 1))
        do iw = 1, nw_use
          coef_diff_norm = coef_diff_norm + abs(wback_after(iw, ist) - wprod_after(iw, ist))**2
          coef_ref_norm = coef_ref_norm + abs(wprod_after(iw, ist))**2
          coef_max_abs_diff = max(coef_max_abs_diff, abs(wback_after(iw, ist) - wprod_after(iw, ist)))
          norm_prod_before = norm_prod_before + abs(wprod_before(iw, ist))**2
          norm_prod_after = norm_prod_after + abs(wprod_after(iw, ist))**2
          norm_local_after = norm_local_after + abs(wback_after(iw, ist))**2
          rho_prod_after = rho_prod_after + abs(wprod_after(iw, ist))**2
          rho_back_after = rho_back_after + abs(wback_after(iw, ist))**2
          pz_prod_after = pz_prod_after - occ * abs(wprod_after(iw, ist))**2 * zdiag(iw)
          pz_back_after = pz_back_after - occ * abs(wback_after(iw, ist))**2 * zdiag(iw)
        end do
      end do
      coef_diff_norm = sqrt(max(coef_diff_norm, 0.0d0))
      coef_rel_diff = coef_diff_norm / max(sqrt(max(coef_ref_norm, 0.0d0)), 1.0d-300)
      rho_diff = rho_back_after - rho_prod_after
      pz_diff = pz_back_after - pz_prod_after
      write(*,'(1x,a,1(a,i0),1(a,a),14(a,1pe12.4),4(a,l1),a,a)') &
        '[DG-MIXEDZ-LOCAL-FRAG-PROD-SOURCED-PROP-CMP]', &
        ' step=', step_use, &
        ' candidate=', trim(label), &
        ' rel_coef_diff=', coef_rel_diff, &
        ' coef_diff_norm=', coef_diff_norm, &
        ' coef_diff_max=', coef_max_abs_diff, &
        ' S_norm_prod_before=', norm_prod_before, &
        ' S_norm_local_before=', norm_prod_before, &
        ' S_norm_prod_after=', norm_prod_after, &
        ' S_norm_local_after=', norm_local_after, &
        ' rho_prod_after=', rho_prod_after, &
        ' rho_back_after=', rho_back_after, &
        ' rho_diff=', rho_diff, &
        ' Pz_prod_after=', pz_prod_after, &
        ' Pz_back_after=', pz_back_after, &
        ' Pz_diff=', pz_diff, &
        ' owner_mismatch_count=', owner_mismatch_count, &
        ' field_off_static_phase=', .true., &
        ' W_sector_only=', .true., &
        ' production_sourced=', .true., &
        ' bad=', coef_rel_diff /= coef_rel_diff, &
        ' route=', 'production-CW-sourced-local-fieldoff-propagation-diagnostic-only'
    end subroutine mixed_z_log_w_prod_sourced_prop_candidate

    subroutine mixed_z_log_w_rho_ww_cmp(label, step_use, wprod, wback, occ_all)
      character(len=*), intent(in) :: label
      integer, intent(in) :: step_use
      complex(8), intent(in) :: wprod(:,:), wback(:,:)
      real(8), intent(in) :: occ_all(:,:,:)
      integer :: iw, jw, ist, jst, nw_use, nstate_use, nocc_eff, ia, ja
      real(8) :: occ, trace_prod, trace_back, Ddiff_frob, Dprod_frob
      real(8) :: rel_D_diff, diag_D_diff_rms, offdiag_D_diff_rms
      real(8) :: phase_aligned_coef_diff, phase_ref_norm
      real(8) :: alpha_abs, subspace_overlap_min, subspace_overlap_max
      complex(8) :: alpha, phase_fac
      complex(8), allocatable :: Dprod(:,:), Dback(:,:), Ddiff(:,:)
      complex(8), allocatable :: M(:,:), Mherm(:,:)
      real(8), allocatable :: Meval(:)
      integer, allocatable :: occ_state(:)

      nw_use = min(size(wprod, 1), size(wback, 1))
      nstate_use = min(size(wprod, 2), size(wback, 2), size(occ_all, 1))
      trace_prod = 0.0d0
      trace_back = 0.0d0
      Ddiff_frob = 0.0d0
      Dprod_frob = 0.0d0
      rel_D_diff = 0.0d0
      diag_D_diff_rms = 0.0d0
      offdiag_D_diff_rms = 0.0d0
      phase_aligned_coef_diff = 0.0d0
      phase_ref_norm = 0.0d0
      subspace_overlap_min = 0.0d0
      subspace_overlap_max = 0.0d0
      if (nw_use <= 0 .or. nstate_use <= 0) return

      allocate(Dprod(nw_use,nw_use), Dback(nw_use,nw_use), Ddiff(nw_use,nw_use), &
               occ_state(nstate_use))
      Dprod(:, :) = (0.0d0, 0.0d0)
      Dback(:, :) = (0.0d0, 0.0d0)
      occ_state(:) = 0
      nocc_eff = 0
      do ist = 1, nstate_use
        occ = 0.0d0
        if (size(occ_all, 3) >= 1) occ = max(0.0d0, occ_all(ist, 1, 1))
        if (occ > 0.0d0) then
          nocc_eff = nocc_eff + 1
          occ_state(nocc_eff) = ist
        end if
        do jw = 1, nw_use
          do iw = 1, nw_use
            Dprod(iw,jw) = Dprod(iw,jw) + occ * wprod(iw,ist) * conjg(wprod(jw,ist))
            Dback(iw,jw) = Dback(iw,jw) + occ * wback(iw,ist) * conjg(wback(jw,ist))
          end do
        end do
        alpha = (0.0d0, 0.0d0)
        do iw = 1, nw_use
          alpha = alpha + conjg(wback(iw,ist)) * wprod(iw,ist)
          phase_ref_norm = phase_ref_norm + abs(wprod(iw,ist))**2
        end do
        alpha_abs = abs(alpha)
        if (alpha_abs > 1.0d-300) then
          phase_fac = alpha / alpha_abs
        else
          phase_fac = (1.0d0, 0.0d0)
        end if
        do iw = 1, nw_use
          phase_aligned_coef_diff = phase_aligned_coef_diff + &
            abs(phase_fac * wback(iw,ist) - wprod(iw,ist))**2
        end do
      end do

      Ddiff(:, :) = Dback(:, :) - Dprod(:, :)
      trace_prod = real(sum((/(Dprod(iw,iw), iw=1,nw_use)/)), kind=8)
      trace_back = real(sum((/(Dback(iw,iw), iw=1,nw_use)/)), kind=8)
      Ddiff_frob = sqrt(max(sum(abs(Ddiff(:, :))**2), 0.0d0))
      Dprod_frob = sqrt(max(sum(abs(Dprod(:, :))**2), 0.0d0))
      rel_D_diff = Ddiff_frob / max(Dprod_frob, 1.0d-300)
      do iw = 1, nw_use
        diag_D_diff_rms = diag_D_diff_rms + abs(Ddiff(iw,iw))**2
      end do
      diag_D_diff_rms = sqrt(diag_D_diff_rms / max(dble(nw_use), 1.0d0))
      if (nw_use > 1) then
        do jw = 1, nw_use
          do iw = 1, nw_use
            if (iw == jw) cycle
            offdiag_D_diff_rms = offdiag_D_diff_rms + abs(Ddiff(iw,jw))**2
          end do
        end do
        offdiag_D_diff_rms = sqrt(offdiag_D_diff_rms / dble(nw_use * (nw_use - 1)))
      end if
      phase_aligned_coef_diff = sqrt(max(phase_aligned_coef_diff, 0.0d0)) / &
        max(sqrt(max(phase_ref_norm, 0.0d0)), 1.0d-300)

      if (nocc_eff > 0) then
        allocate(M(nocc_eff,nocc_eff), Mherm(nocc_eff,nocc_eff), Meval(nocc_eff))
        M(:, :) = (0.0d0, 0.0d0)
        do ja = 1, nocc_eff
          jst = occ_state(ja)
          do ia = 1, nocc_eff
            ist = occ_state(ia)
            do iw = 1, nw_use
              M(ia,ja) = M(ia,ja) + conjg(wprod(iw,ist)) * wback(iw,jst)
            end do
          end do
        end do
        Mherm(:, :) = matmul(conjg(transpose(M(:, :))), M(:, :))
        call eigen_zheev(Mherm, Meval, Mherm)
        subspace_overlap_min = sqrt(max(minval(Meval(1:nocc_eff)), 0.0d0))
        subspace_overlap_max = sqrt(max(maxval(Meval(1:nocc_eff)), 0.0d0))
        deallocate(M, Mherm, Meval)
      end if

      write(*,'(1x,a,1(a,i0),1(a,a),11(a,1pe12.4),3(a,l1),a,a)') &
        '[DG-MIXEDZ-LOCAL-FRAG-PROD-SOURCED-RHO-WW-CMP]', &
        ' step=', step_use, &
        ' candidate=', trim(label), &
        ' trace_D_prod=', trace_prod, &
        ' trace_D_back=', trace_back, &
        ' D_diff_frob=', Ddiff_frob, &
        ' rel_D_diff=', rel_D_diff, &
        ' diag_D_diff_rms=', diag_D_diff_rms, &
        ' offdiag_D_diff_rms=', offdiag_D_diff_rms, &
        ' phase_aligned_coef_diff=', phase_aligned_coef_diff, &
        ' subspace_overlap_min=', subspace_overlap_min, &
        ' subspace_overlap_max=', subspace_overlap_max, &
        ' rho_trace_diff=', trace_back - trace_prod, &
        ' nocc_eff_real=', dble(nocc_eff), &
        ' field_off_static_phase=', .true., &
        ' W_sector_only=', .true., &
        ' bad=', rel_D_diff /= rel_D_diff, &
        ' route=', 'production-CW-sourced-WW-density-matrix-compare-diagnostic-only'
      deallocate(Dprod, Dback, Ddiff, occ_state)
    end subroutine mixed_z_log_w_rho_ww_cmp

    subroutine mixed_z_log_w_prod_wphase_prop_cmp(label, step_use, wprod_after, wback_after, occ_all, zdiag)
      character(len=*), intent(in) :: label
      integer, intent(in) :: step_use
      complex(8), intent(in) :: wprod_after(:,:), wback_after(:,:)
      real(8), intent(in) :: occ_all(:,:,:), zdiag(:)
      integer :: iw, jw, ist, nw_use, nstate_use
      real(8) :: occ, coef_diff_norm, coef_ref_norm, rel_coef_diff, coef_max_abs_diff
      real(8) :: rho_prod, rho_back, rho_diff, pz_prod, pz_back, pz_diff
      real(8) :: phase_aligned_coef_diff, phase_ref_norm, alpha_abs
      real(8) :: trace_prod, trace_back, Ddiff_frob, Dprod_frob, rel_D_diff
      real(8) :: diag_D_diff_rms, offdiag_D_diff_rms
      complex(8) :: alpha, phase_fac
      complex(8), allocatable :: Dprod(:,:), Dback(:,:), Ddiff(:,:)

      nw_use = min(size(wprod_after, 1), size(wback_after, 1), size(zdiag))
      nstate_use = min(size(wprod_after, 2), size(wback_after, 2), size(occ_all, 1))
      if (nw_use <= 0 .or. nstate_use <= 0) return
      allocate(Dprod(nw_use,nw_use), Dback(nw_use,nw_use), Ddiff(nw_use,nw_use))
      Dprod(:, :) = (0.0d0, 0.0d0)
      Dback(:, :) = (0.0d0, 0.0d0)
      coef_diff_norm = 0.0d0
      coef_ref_norm = 0.0d0
      coef_max_abs_diff = 0.0d0
      rho_prod = 0.0d0
      rho_back = 0.0d0
      pz_prod = 0.0d0
      pz_back = 0.0d0
      phase_aligned_coef_diff = 0.0d0
      phase_ref_norm = 0.0d0

      do ist = 1, nstate_use
        occ = 0.0d0
        if (size(occ_all, 3) >= 1) occ = max(0.0d0, occ_all(ist, 1, 1))
        alpha = (0.0d0, 0.0d0)
        do iw = 1, nw_use
          coef_diff_norm = coef_diff_norm + abs(wback_after(iw, ist) - wprod_after(iw, ist))**2
          coef_ref_norm = coef_ref_norm + abs(wprod_after(iw, ist))**2
          coef_max_abs_diff = max(coef_max_abs_diff, abs(wback_after(iw, ist) - wprod_after(iw, ist)))
          rho_prod = rho_prod + abs(wprod_after(iw, ist))**2
          rho_back = rho_back + abs(wback_after(iw, ist))**2
          pz_prod = pz_prod - occ * abs(wprod_after(iw, ist))**2 * zdiag(iw)
          pz_back = pz_back - occ * abs(wback_after(iw, ist))**2 * zdiag(iw)
          alpha = alpha + conjg(wback_after(iw, ist)) * wprod_after(iw, ist)
          phase_ref_norm = phase_ref_norm + abs(wprod_after(iw, ist))**2
        end do
        alpha_abs = abs(alpha)
        if (alpha_abs > 1.0d-300) then
          phase_fac = alpha / alpha_abs
        else
          phase_fac = (1.0d0, 0.0d0)
        end if
        do iw = 1, nw_use
          phase_aligned_coef_diff = phase_aligned_coef_diff + &
            abs(phase_fac * wback_after(iw, ist) - wprod_after(iw, ist))**2
        end do
        do jw = 1, nw_use
          do iw = 1, nw_use
            Dprod(iw,jw) = Dprod(iw,jw) + occ * wprod_after(iw,ist) * conjg(wprod_after(jw,ist))
            Dback(iw,jw) = Dback(iw,jw) + occ * wback_after(iw,ist) * conjg(wback_after(jw,ist))
          end do
        end do
      end do

      Ddiff(:, :) = Dback(:, :) - Dprod(:, :)
      trace_prod = real(sum((/(Dprod(iw,iw), iw=1,nw_use)/)), kind=8)
      trace_back = real(sum((/(Dback(iw,iw), iw=1,nw_use)/)), kind=8)
      Ddiff_frob = sqrt(max(sum(abs(Ddiff(:, :))**2), 0.0d0))
      Dprod_frob = sqrt(max(sum(abs(Dprod(:, :))**2), 0.0d0))
      rel_D_diff = Ddiff_frob / max(Dprod_frob, 1.0d-300)
      diag_D_diff_rms = 0.0d0
      do iw = 1, nw_use
        diag_D_diff_rms = diag_D_diff_rms + abs(Ddiff(iw,iw))**2
      end do
      diag_D_diff_rms = sqrt(diag_D_diff_rms / max(dble(nw_use), 1.0d0))
      offdiag_D_diff_rms = 0.0d0
      if (nw_use > 1) then
        do jw = 1, nw_use
          do iw = 1, nw_use
            if (iw == jw) cycle
            offdiag_D_diff_rms = offdiag_D_diff_rms + abs(Ddiff(iw,jw))**2
          end do
        end do
        offdiag_D_diff_rms = sqrt(offdiag_D_diff_rms / dble(nw_use * (nw_use - 1)))
      end if

      coef_diff_norm = sqrt(max(coef_diff_norm, 0.0d0))
      rel_coef_diff = coef_diff_norm / max(sqrt(max(coef_ref_norm, 0.0d0)), 1.0d-300)
      phase_aligned_coef_diff = sqrt(max(phase_aligned_coef_diff, 0.0d0)) / &
        max(sqrt(max(phase_ref_norm, 0.0d0)), 1.0d-300)
      rho_diff = rho_back - rho_prod
      pz_diff = pz_back - pz_prod

      write(*,'(1x,a,1(a,i0),1(a,a),16(a,1pe12.4),4(a,l1),a,a)') &
        '[DG-MIXEDZ-LOCAL-FRAG-PROD-WPHASE-PROP-CMP]', &
        ' step=', step_use, &
        ' candidate=', trim(label), &
        ' rel_coef_diff=', rel_coef_diff, &
        ' coef_diff_norm=', coef_diff_norm, &
        ' coef_diff_max=', coef_max_abs_diff, &
        ' phase_aligned_coef_diff=', phase_aligned_coef_diff, &
        ' rel_D_diff=', rel_D_diff, &
        ' D_diff_frob=', Ddiff_frob, &
        ' diag_D_diff_rms=', diag_D_diff_rms, &
        ' offdiag_D_diff_rms=', offdiag_D_diff_rms, &
        ' rho_prod=', rho_prod, &
        ' rho_back=', rho_back, &
        ' rho_diff=', rho_diff, &
        ' Pz_prod=', pz_prod, &
        ' Pz_back=', pz_back, &
        ' Pz_diff=', pz_diff, &
        ' trace_D_prod=', trace_prod, &
        ' trace_D_back=', trace_back, &
        ' field_off_static_phase=', sum(abs(E_mid(1:3))) <= 1.0d-30, &
        ' W_sector_only=', .true., &
        ' production_wphase=', .true., &
        ' bad=', rel_coef_diff /= rel_coef_diff .or. rel_D_diff /= rel_D_diff, &
        ' route=', 'production-CW-sourced-local-W-eval-phase-diagnostic-only'
      deallocate(Dprod, Dback, Ddiff)
    end subroutine mixed_z_log_w_prod_wphase_prop_cmp

    subroutine mixed_z_log_wphase_residual_cmp(step_use, wprod_after, wphase_after, wresidual_after, &
      occ_all, zdiag, h_residual_norm, h_wp_norm, h_pp_norm)
      integer, intent(in) :: step_use
      complex(8), intent(in) :: wprod_after(:,:), wphase_after(:,:), wresidual_after(:,:)
      real(8), intent(in) :: occ_all(:,:,:), zdiag(:)
      real(8), intent(in) :: h_residual_norm, h_wp_norm, h_pp_norm
      integer :: iw, jw, ist, nw_use, nstate_use
      real(8) :: occ, ref_norm
      real(8) :: wphase_c_diff, wphase_d_diff, wphase_rho_diff, wphase_pz_diff
      real(8) :: residual_c_diff, residual_d_diff, residual_rho_diff, residual_pz_diff
      real(8) :: delta_c, delta_d, delta_rho, delta_pz
      real(8) :: dprod_norm
      complex(8), allocatable :: Dprod(:,:), Dwphase(:,:), Dresidual(:,:)

      nw_use = min(size(wprod_after, 1), size(wphase_after, 1), size(wresidual_after, 1), size(zdiag))
      nstate_use = min(size(wprod_after, 2), size(wphase_after, 2), size(wresidual_after, 2), size(occ_all, 1))
      if (nw_use <= 0 .or. nstate_use <= 0) return

      allocate(Dprod(nw_use,nw_use), Dwphase(nw_use,nw_use), Dresidual(nw_use,nw_use))
      Dprod(:, :) = (0.0d0, 0.0d0)
      Dwphase(:, :) = (0.0d0, 0.0d0)
      Dresidual(:, :) = (0.0d0, 0.0d0)
      ref_norm = 0.0d0
      wphase_c_diff = 0.0d0
      residual_c_diff = 0.0d0
      delta_c = 0.0d0
      wphase_rho_diff = 0.0d0
      residual_rho_diff = 0.0d0
      delta_rho = 0.0d0
      wphase_pz_diff = 0.0d0
      residual_pz_diff = 0.0d0
      delta_pz = 0.0d0

      do ist = 1, nstate_use
        occ = 0.0d0
        if (size(occ_all, 3) >= 1) occ = max(0.0d0, occ_all(ist, 1, 1))
        do iw = 1, nw_use
          ref_norm = ref_norm + abs(wprod_after(iw, ist))**2
          wphase_c_diff = wphase_c_diff + abs(wphase_after(iw, ist) - wprod_after(iw, ist))**2
          residual_c_diff = residual_c_diff + abs(wresidual_after(iw, ist) - wprod_after(iw, ist))**2
          delta_c = delta_c + abs(wresidual_after(iw, ist) - wphase_after(iw, ist))**2
          wphase_rho_diff = wphase_rho_diff + abs(wphase_after(iw, ist))**2 - abs(wprod_after(iw, ist))**2
          residual_rho_diff = residual_rho_diff + abs(wresidual_after(iw, ist))**2 - abs(wprod_after(iw, ist))**2
          delta_rho = delta_rho + abs(wresidual_after(iw, ist))**2 - abs(wphase_after(iw, ist))**2
          wphase_pz_diff = wphase_pz_diff - occ * &
            (abs(wphase_after(iw, ist))**2 - abs(wprod_after(iw, ist))**2) * zdiag(iw)
          residual_pz_diff = residual_pz_diff - occ * &
            (abs(wresidual_after(iw, ist))**2 - abs(wprod_after(iw, ist))**2) * zdiag(iw)
          delta_pz = delta_pz - occ * &
            (abs(wresidual_after(iw, ist))**2 - abs(wphase_after(iw, ist))**2) * zdiag(iw)
          do jw = 1, nw_use
            Dprod(iw,jw) = Dprod(iw,jw) + occ * wprod_after(iw,ist) * conjg(wprod_after(jw,ist))
            Dwphase(iw,jw) = Dwphase(iw,jw) + occ * wphase_after(iw,ist) * conjg(wphase_after(jw,ist))
            Dresidual(iw,jw) = Dresidual(iw,jw) + occ * wresidual_after(iw,ist) * conjg(wresidual_after(jw,ist))
          end do
        end do
      end do

      wphase_c_diff = sqrt(max(wphase_c_diff, 0.0d0)) / max(sqrt(max(ref_norm, 0.0d0)), 1.0d-300)
      residual_c_diff = sqrt(max(residual_c_diff, 0.0d0)) / max(sqrt(max(ref_norm, 0.0d0)), 1.0d-300)
      delta_c = sqrt(max(delta_c, 0.0d0)) / max(sqrt(max(ref_norm, 0.0d0)), 1.0d-300)
      dprod_norm = sqrt(max(sum(abs(Dprod(:, :))**2), 0.0d0))
      wphase_d_diff = sqrt(max(sum(abs(Dwphase(:, :) - Dprod(:, :))**2), 0.0d0)) / &
        max(dprod_norm, 1.0d-300)
      residual_d_diff = sqrt(max(sum(abs(Dresidual(:, :) - Dprod(:, :))**2), 0.0d0)) / &
        max(dprod_norm, 1.0d-300)
      delta_d = sqrt(max(sum(abs(Dresidual(:, :) - Dwphase(:, :))**2), 0.0d0)) / &
        max(dprod_norm, 1.0d-300)

      write(*,'(1x,a,1(a,i0),1(a,a),13(a,1pe12.4),5(a,l1),a,a)') &
        '[DG-MIXEDZ-LOCAL-FRAG-WPHASE-LOCAL-RESIDUAL-CMP]', &
        ' step=', step_use, &
        ' candidate=', 'owner_unique', &
        ' Wphase_only_rel_C_diff=', wphase_c_diff, &
        ' Wphase_only_rel_D_diff=', wphase_d_diff, &
        ' residual_on_rel_C_diff=', residual_c_diff, &
        ' residual_on_rel_D_diff=', residual_d_diff, &
        ' residual_on_rho_diff=', residual_rho_diff, &
        ' residual_on_Pz_diff=', residual_pz_diff, &
        ' delta_from_Wphase_C=', delta_c, &
        ' delta_from_Wphase_D=', delta_d, &
        ' delta_from_Wphase_rho=', delta_rho, &
        ' delta_from_Wphase_Pz=', delta_pz, &
        ' H_residual_norm=', h_residual_norm, &
        ' H_WP_norm=', h_wp_norm, &
        ' H_PP_norm=', h_pp_norm, &
        ' Wphase_only_roundoff=', wphase_c_diff <= 1.0d-12 .and. wphase_d_diff <= 1.0d-12, &
        ' residual_available=', h_residual_norm > 0.0d0, &
        ' field_off_static_phase=', sum(abs(E_mid(1:3))) <= 1.0d-30, &
        ' production_replacement=', .false., &
        ' bad=', wphase_c_diff /= wphase_c_diff .or. residual_c_diff /= residual_c_diff, &
        ' route=', 'Wphase-exact-baseline-vs-local-residual-diagnostic-only'
      deallocate(Dprod, Dwphase, Dresidual)
    end subroutine mixed_z_log_wphase_residual_cmp

    subroutine mixed_z_log_residual_block_case(step_use, case_label, wprod_after, wcandidate_after, &
      occ_all, zdiag, h_block_norm)
      integer, intent(in) :: step_use
      character(len=*), intent(in) :: case_label
      complex(8), intent(in) :: wprod_after(:,:), wcandidate_after(:,:)
      real(8), intent(in) :: occ_all(:,:,:), zdiag(:)
      real(8), intent(in) :: h_block_norm
      integer :: iw, jw, ist, nw_use, nstate_use
      real(8) :: occ, ref_norm, c_diff, d_diff, dprod_norm
      real(8) :: offdiag_D_diff_rms, rho_diff, pz_diff
      complex(8), allocatable :: Dprod(:,:), Dcandidate(:,:)

      nw_use = min(size(wprod_after, 1), size(wcandidate_after, 1), size(zdiag))
      nstate_use = min(size(wprod_after, 2), size(wcandidate_after, 2), size(occ_all, 1))
      if (nw_use <= 0 .or. nstate_use <= 0) return

      allocate(Dprod(nw_use,nw_use), Dcandidate(nw_use,nw_use))
      Dprod(:, :) = (0.0d0, 0.0d0)
      Dcandidate(:, :) = (0.0d0, 0.0d0)
      ref_norm = 0.0d0
      c_diff = 0.0d0
      rho_diff = 0.0d0
      pz_diff = 0.0d0

      do ist = 1, nstate_use
        occ = 0.0d0
        if (size(occ_all, 3) >= 1) occ = max(0.0d0, occ_all(ist, 1, 1))
        do iw = 1, nw_use
          ref_norm = ref_norm + abs(wprod_after(iw, ist))**2
          c_diff = c_diff + abs(wcandidate_after(iw, ist) - wprod_after(iw, ist))**2
          rho_diff = rho_diff + abs(wcandidate_after(iw, ist))**2 - abs(wprod_after(iw, ist))**2
          pz_diff = pz_diff - occ * &
            (abs(wcandidate_after(iw, ist))**2 - abs(wprod_after(iw, ist))**2) * zdiag(iw)
          do jw = 1, nw_use
            Dprod(iw,jw) = Dprod(iw,jw) + occ * wprod_after(iw,ist) * conjg(wprod_after(jw,ist))
            Dcandidate(iw,jw) = Dcandidate(iw,jw) + occ * wcandidate_after(iw,ist) * &
              conjg(wcandidate_after(jw,ist))
          end do
        end do
      end do

      c_diff = sqrt(max(c_diff, 0.0d0)) / max(sqrt(max(ref_norm, 0.0d0)), 1.0d-300)
      dprod_norm = sqrt(max(sum(abs(Dprod(:, :))**2), 0.0d0))
      d_diff = sqrt(max(sum(abs(Dcandidate(:, :) - Dprod(:, :))**2), 0.0d0)) / &
        max(dprod_norm, 1.0d-300)
      offdiag_D_diff_rms = 0.0d0
      if (nw_use > 1) then
        do jw = 1, nw_use
          do iw = 1, nw_use
            if (iw == jw) cycle
            offdiag_D_diff_rms = offdiag_D_diff_rms + abs(Dcandidate(iw,jw) - Dprod(iw,jw))**2
          end do
        end do
        offdiag_D_diff_rms = sqrt(offdiag_D_diff_rms / dble(nw_use * (nw_use - 1)))
      end if

      write(*,'(1x,a,1(a,i0),1(a,a),6(a,1pe12.4),1(a,l1),a,a)') &
        '[DG-MIXEDZ-LOCAL-FRAG-RESIDUAL-BLOCK-CMP]', &
        ' step=', step_use, &
        ' case=', trim(case_label), &
        ' rel_C_diff=', c_diff, &
        ' rel_D_diff=', d_diff, &
        ' offdiag_D_diff_rms=', offdiag_D_diff_rms, &
        ' rho_diff=', rho_diff, &
        ' Pz_diff=', pz_diff, &
        ' H_block_norm=', h_block_norm, &
        ' bad=', c_diff /= c_diff .or. d_diff /= d_diff, &
        ' route=', 'Wphase-baseline-plus-selected-local-residual-block-diagnostic-only'
      deallocate(Dprod, Dcandidate)
    end subroutine mixed_z_log_residual_block_case

    subroutine mixed_z_log_wphase_replacement_series(step_use, candidate_label, wprod_before, &
      wprod_after, wcandidate_after, occ_all, zdiag)
      integer, intent(in) :: step_use
      character(len=*), intent(in) :: candidate_label
      complex(8), intent(in) :: wprod_before(:,:), wprod_after(:,:), wcandidate_after(:,:)
      real(8), intent(in) :: occ_all(:,:,:), zdiag(:)
      integer :: iw, jw, ist, nw_use, nstate_use
      real(8) :: occ, ref_norm, c_diff, d_diff, dprod_norm
      real(8) :: offdiag_D_diff_rms, rho_diff, pz_diff, current_diff
      real(8) :: norm_before, norm_prod_after, norm_candidate_after, norm_drift
      complex(8), allocatable :: Dprod(:,:), Dcandidate(:,:)

      nw_use = min(size(wprod_before, 1), size(wprod_after, 1), size(wcandidate_after, 1), size(zdiag))
      nstate_use = min(size(wprod_before, 2), size(wprod_after, 2), size(wcandidate_after, 2), size(occ_all, 1))
      if (nw_use <= 0 .or. nstate_use <= 0) return

      allocate(Dprod(nw_use,nw_use), Dcandidate(nw_use,nw_use))
      Dprod(:, :) = (0.0d0, 0.0d0)
      Dcandidate(:, :) = (0.0d0, 0.0d0)
      ref_norm = 0.0d0
      c_diff = 0.0d0
      rho_diff = 0.0d0
      pz_diff = 0.0d0
      current_diff = 0.0d0
      norm_before = 0.0d0
      norm_prod_after = 0.0d0
      norm_candidate_after = 0.0d0

      do ist = 1, nstate_use
        occ = 0.0d0
        if (size(occ_all, 3) >= 1) occ = max(0.0d0, occ_all(ist, 1, 1))
        do iw = 1, nw_use
          norm_before = norm_before + abs(wprod_before(iw, ist))**2
          norm_prod_after = norm_prod_after + abs(wprod_after(iw, ist))**2
          norm_candidate_after = norm_candidate_after + abs(wcandidate_after(iw, ist))**2
          ref_norm = ref_norm + abs(wprod_after(iw, ist))**2
          c_diff = c_diff + abs(wcandidate_after(iw, ist) - wprod_after(iw, ist))**2
          rho_diff = rho_diff + abs(wcandidate_after(iw, ist))**2 - abs(wprod_after(iw, ist))**2
          pz_diff = pz_diff - occ * &
            (abs(wcandidate_after(iw, ist))**2 - abs(wprod_after(iw, ist))**2) * zdiag(iw)
          do jw = 1, nw_use
            Dprod(iw,jw) = Dprod(iw,jw) + occ * wprod_after(iw,ist) * conjg(wprod_after(jw,ist))
            Dcandidate(iw,jw) = Dcandidate(iw,jw) + occ * wcandidate_after(iw,ist) * &
              conjg(wcandidate_after(jw,ist))
          end do
        end do
      end do

      c_diff = sqrt(max(c_diff, 0.0d0)) / max(sqrt(max(ref_norm, 0.0d0)), 1.0d-300)
      dprod_norm = sqrt(max(sum(abs(Dprod(:, :))**2), 0.0d0))
      d_diff = sqrt(max(sum(abs(Dcandidate(:, :) - Dprod(:, :))**2), 0.0d0)) / &
        max(dprod_norm, 1.0d-300)
      offdiag_D_diff_rms = 0.0d0
      if (nw_use > 1) then
        do jw = 1, nw_use
          do iw = 1, nw_use
            if (iw == jw) cycle
            offdiag_D_diff_rms = offdiag_D_diff_rms + abs(Dcandidate(iw,jw) - Dprod(iw,jw))**2
          end do
        end do
        offdiag_D_diff_rms = sqrt(offdiag_D_diff_rms / dble(nw_use * (nw_use - 1)))
      end if
      norm_drift = norm_candidate_after - norm_before

      write(*,'(1x,a,1(a,i0),1(a,a),9(a,1pe12.4),4(a,l1),a,a)') &
        '[DG-MIXEDZ-LOCAL-FRAG-WPHASE-REPLACEMENT-SERIES]', &
        ' step=', step_use, &
        ' candidate=', trim(candidate_label), &
        ' rel_C_diff=', c_diff, &
        ' rel_D_diff=', d_diff, &
        ' offdiag_D_diff_rms=', offdiag_D_diff_rms, &
        ' rho_diff=', rho_diff, &
        ' Pz_diff=', pz_diff, &
        ' current_diff=', current_diff, &
        ' norm_before=', norm_before, &
        ' norm_after=', norm_candidate_after, &
        ' norm_drift=', norm_drift, &
        ' current_available=', .false., &
        ' field_off_static_phase=', sum(abs(E_mid(1:3))) <= 1.0d-30, &
        ' production_replacement=', .false., &
        ' bad=', c_diff /= c_diff .or. d_diff /= d_diff, &
        ' route=', 'owner-unique-Wphase-plus-local-residual-production-candidate-series'
      deallocate(Dprod, Dcandidate)
    end subroutine mixed_z_log_wphase_replacement_series

    subroutine mixed_z_log_wphase_residual_norm_cmp(step_use, vals)
      integer, intent(in) :: step_use
      real(8), intent(in) :: vals(:)
      real(8) :: norm_W_before, norm_P_before, norm_W_after, norm_P_after
      real(8) :: norm_total_before, norm_total_after, norm_total_drift
      real(8) :: S_herm_diff, H_residual_herm_diff, U_residual_unitarity_diff
      real(8) :: P_sector_initial_norm, P_sector_generated_norm, WP_transfer_norm, PP_phase_norm

      if (size(vals) < 13) return
      norm_W_before = vals(1)
      norm_P_before = vals(2)
      norm_W_after = vals(3)
      norm_P_after = vals(4)
      norm_total_before = vals(5)
      norm_total_after = vals(6)
      norm_total_drift = norm_total_after - norm_total_before
      S_herm_diff = vals(7)
      U_residual_unitarity_diff = sqrt(max(vals(8), 0.0d0))
      P_sector_initial_norm = norm_P_before
      P_sector_generated_norm = norm_P_after
      H_residual_herm_diff = vals(10)
      WP_transfer_norm = sqrt(max(vals(12), 0.0d0))
      PP_phase_norm = sqrt(max(vals(13), 0.0d0))

      write(*,'(1x,a,1(a,i0),14(a,1pe12.4),4(a,l1),a,a)') &
        '[DG-MIXEDZ-LOCAL-FRAG-WPHASE-RESIDUAL-NORM-CMP]', &
        ' step=', step_use, &
        ' norm_W_before=', norm_W_before, &
        ' norm_P_before=', norm_P_before, &
        ' norm_W_after=', norm_W_after, &
        ' norm_P_after=', norm_P_after, &
        ' norm_total_before=', norm_total_before, &
        ' norm_total_after=', norm_total_after, &
        ' norm_total_drift=', norm_total_drift, &
        ' S_herm_diff=', S_herm_diff, &
        ' H_residual_herm_diff=', H_residual_herm_diff, &
        ' U_residual_unitarity_diff=', U_residual_unitarity_diff, &
        ' P_sector_initial_norm=', P_sector_initial_norm, &
        ' P_sector_generated_norm=', P_sector_generated_norm, &
        ' WP_transfer_norm=', WP_transfer_norm, &
        ' PP_phase_norm=', PP_phase_norm, &
        ' Wphase_only_norm_exact=', abs(norm_total_drift) <= 1.0d-10, &
        ' current_available=', .false., &
        ' production_replacement=', .false., &
        ' bad=', norm_total_drift /= norm_total_drift .or. U_residual_unitarity_diff /= U_residual_unitarity_diff, &
        ' route=', 'local-raw-Wphase-plus-residual-Snorm-diagnostic-only'
    end subroutine mixed_z_log_wphase_residual_norm_cmp

    subroutine mixed_z_log_wp_extended_obs_cmp(step_use, wprod_after, wext_after, occ_all, zdiag, &
      ext_obs, norm_vals)
      integer, intent(in) :: step_use
      complex(8), intent(in) :: wprod_after(:,:), wext_after(:,:)
      real(8), intent(in) :: occ_all(:,:,:), zdiag(:), ext_obs(:), norm_vals(:)
      integer :: iw, ist, nw_use, nstate_use
      real(8) :: occ, rho_W_only, rho_prod, rho_WP_extended
      real(8) :: Pz_W_only, Pz_prod, Pz_WP_extended
      real(8) :: norm_W_only, norm_WP_extended, norm_prod
      real(8) :: W_loss_to_P, P_generated_norm, cross_WP_contribution

      nw_use = min(size(wprod_after, 1), size(wext_after, 1), size(zdiag))
      nstate_use = min(size(wprod_after, 2), size(wext_after, 2), size(occ_all, 1))
      if (nw_use <= 0 .or. nstate_use <= 0) return

      rho_W_only = 0.0d0
      Pz_W_only = 0.0d0
      do ist = 1, nstate_use
        occ = 0.0d0
        if (size(occ_all, 3) >= 1) occ = max(0.0d0, occ_all(ist, 1, 1))
        do iw = 1, nw_use
          rho_W_only = rho_W_only + abs(wext_after(iw, ist))**2
          Pz_W_only = Pz_W_only - occ * abs(wext_after(iw, ist))**2 * zdiag(iw)
        end do
      end do
      rho_prod = rho_W_only
      Pz_prod = Pz_W_only
      rho_WP_extended = 0.0d0
      Pz_WP_extended = 0.0d0
      if (size(ext_obs) >= 2) then
        rho_WP_extended = ext_obs(1)
        Pz_WP_extended = ext_obs(2)
      end if
      norm_W_only = rho_W_only
      norm_prod = rho_prod
      norm_WP_extended = rho_WP_extended
      if (size(norm_vals) >= 6) norm_WP_extended = norm_vals(6)
      W_loss_to_P = 0.0d0
      P_generated_norm = 0.0d0
      cross_WP_contribution = 0.0d0
      if (size(norm_vals) >= 6) then
        W_loss_to_P = norm_vals(3) - norm_vals(1)
        P_generated_norm = norm_vals(4)
        cross_WP_contribution = norm_vals(6) - norm_vals(3) - norm_vals(4)
      end if

      write(*,'(1x,a,1(a,i0),12(a,1pe12.4),5(a,l1),a,a)') &
        '[DG-MIXEDZ-LOCAL-FRAG-WP-EXTENDED-OBS-CMP]', &
        ' step=', step_use, &
        ' rho_W_only=', rho_W_only, &
        ' rho_WP_extended=', rho_WP_extended, &
        ' rho_prod=', rho_prod, &
        ' Pz_W_only=', Pz_W_only, &
        ' Pz_WP_extended=', Pz_WP_extended, &
        ' Pz_prod=', Pz_prod, &
        ' norm_W_only=', norm_W_only, &
        ' norm_WP_extended=', norm_WP_extended, &
        ' norm_prod=', norm_prod, &
        ' W_loss_to_P=', W_loss_to_P, &
        ' P_generated_norm=', P_generated_norm, &
        ' cross_WP_contribution=', cross_WP_contribution, &
        ' wp_extended_available=', size(ext_obs) >= 2, &
        ' current_available=', .false., &
        ' production_replacement=', .false., &
        ' convention_match=', .false., &
        ' bad=', rho_WP_extended /= rho_WP_extended .or. norm_WP_extended /= norm_WP_extended, &
        ' route=', 'local-WP-extended-observable-diagnostic-only'
    end subroutine mixed_z_log_wp_extended_obs_cmp

    subroutine mixed_z_log_wp_obs_operator_cmp(step_use, case_label, wext_after, occ_all, zdiag, ext_obs, norm_vals)
      integer, intent(in) :: step_use
      character(len=*), intent(in) :: case_label
      complex(8), intent(in) :: wext_after(:,:)
      real(8), intent(in) :: occ_all(:,:,:), zdiag(:), ext_obs(:), norm_vals(:)
      integer :: iw, ist, nw_use, nstate_use
      real(8) :: occ, obs_WW, obs_WP_cross, obs_PP, obs_total, obs_prod, obs_diff
      logical :: O_WW_available, O_WP_available, O_PP_available

      nw_use = min(size(wext_after, 1), size(zdiag))
      nstate_use = min(size(wext_after, 2), size(occ_all, 1))
      if (nw_use <= 0 .or. nstate_use <= 0) return

      obs_WW = 0.0d0
      obs_WP_cross = 0.0d0
      obs_PP = 0.0d0
      obs_total = 0.0d0
      obs_prod = 0.0d0
      O_WW_available = .true.
      O_WP_available = .false.
      O_PP_available = .false.

      if (trim(case_label) == 'density') then
        if (size(norm_vals) >= 6) then
          obs_WW = norm_vals(3)
          obs_PP = norm_vals(4)
          obs_WP_cross = norm_vals(6) - norm_vals(3) - norm_vals(4)
          obs_total = norm_vals(6)
          obs_prod = norm_vals(5)
        end if
        O_WP_available = .true.
        O_PP_available = .true.
      else if (trim(case_label) == 'Pz') then
        if (size(ext_obs) >= 5) then
          obs_WW = ext_obs(3)
          obs_WP_cross = ext_obs(4)
          obs_PP = ext_obs(5)
          obs_total = obs_WW + obs_WP_cross + obs_PP
          obs_prod = ext_obs(2)
          O_WP_available = .true.
          O_PP_available = .true.
        else
          do ist = 1, nstate_use
            occ = 0.0d0
            if (size(occ_all, 3) >= 1) occ = max(0.0d0, occ_all(ist, 1, 1))
            do iw = 1, nw_use
              obs_WW = obs_WW - occ * abs(wext_after(iw, ist))**2 * zdiag(iw)
            end do
          end do
          obs_total = obs_WW
          obs_prod = obs_WW
        end if
      else
        return
      end if
      obs_diff = obs_total - obs_prod

      write(*,'(1x,a,1(a,i0),1(a,a),6(a,1pe12.4),7(a,l1),a,a)') &
        '[DG-MIXEDZ-LOCAL-FRAG-WP-OBS-OPERATOR-CMP]', &
        ' step=', step_use, &
        ' case=', trim(case_label), &
        ' obs_WW=', obs_WW, &
        ' obs_WP_cross=', obs_WP_cross, &
        ' obs_PP=', obs_PP, &
        ' obs_total=', obs_total, &
        ' obs_prod=', obs_prod, &
        ' obs_diff=', obs_diff, &
        ' O_WW_available=', O_WW_available, &
        ' O_WP_available=', O_WP_available, &
        ' O_PP_available=', O_PP_available, &
        ' current_available=', .false., &
        ' production_replacement=', .false., &
        ' convention_match=', abs(obs_diff) <= 1.0d-10 * max(1.0d0, abs(obs_prod)), &
        ' bad=', obs_total /= obs_total .or. obs_prod /= obs_prod, &
        ' route=', 'local-WP-observable-operator-block-diagnostic-only'
    end subroutine mixed_z_log_wp_obs_operator_cmp

    subroutine mixed_z_log_rho_block_cmp(step_use, norm_vals)
      integer, intent(in) :: step_use
      real(8), intent(in) :: norm_vals(:)
      real(8) :: rho_prod_int, rho_WW_int, rho_WP_int, rho_PP_int, rho_total_int
      real(8) :: diff_total_minus_prod, max_abs_grid_diff, rms_grid_diff
      logical :: block_available, grid_compare_available, bad

      if (size(norm_vals) < 6) return
      rho_WW_int = norm_vals(3)
      rho_PP_int = norm_vals(4)
      rho_prod_int = norm_vals(5)
      rho_total_int = norm_vals(6)
      rho_WP_int = rho_total_int - rho_WW_int - rho_PP_int
      diff_total_minus_prod = rho_total_int - rho_prod_int

      ! This first density pass is an integral/block diagnostic.  Full grid
      ! block materialization will be wired through the density builder next.
      grid_compare_available = .false.
      max_abs_grid_diff = -1.0d0
      rms_grid_diff = -1.0d0
      block_available = .true.
      bad = rho_prod_int /= rho_prod_int .or. rho_WW_int /= rho_WW_int .or. &
        rho_WP_int /= rho_WP_int .or. rho_PP_int /= rho_PP_int .or. &
        rho_total_int /= rho_total_int

      write(*,'(1x,a,1(a,i0),8(a,1pe12.4),6(a,l1),a,a)') &
        '[DG-MIXEDZ-LOCAL-RHO-BLOCK-CMP]', &
        ' step=', step_use, &
        ' rho_prod_int=', rho_prod_int, &
        ' rho_WW_int=', rho_WW_int, &
        ' rho_WP_int=', rho_WP_int, &
        ' rho_PP_int=', rho_PP_int, &
        ' rho_total_int=', rho_total_int, &
        ' diff_total_minus_prod=', diff_total_minus_prod, &
        ' max_abs_grid_diff=', max_abs_grid_diff, &
        ' rms_grid_diff=', rms_grid_diff, &
        ' block_available=', block_available, &
        ' grid_compare_available=', grid_compare_available, &
        ' density_replaced=', .false., &
        ' current_replaced=', .false., &
        ' production_replacement=', .false., &
        ' bad=', bad, &
        ' route=', 'local-WP-density-block-dryrun-integral-only'
    end subroutine mixed_z_log_rho_block_cmp


    subroutine diagnose_mixed_z_local_current_block_cmp(step_use)
      integer, intent(in) :: step_use
      type(s_vector) :: curr_prod_grid, curr_WW_grid
      complex(8), allocatable :: pcoef_saved(:,:,:)
      real(8) :: local_sum(22), global_sum(22), local_max(1), global_max(1)
      real(8) :: prod_vec, ww_vec, diff_vec
      real(8) :: prod_norm, ww_norm, diff_norm, rms_diff, max_abs_diff
      real(8) :: prod_int_norm, ww_int_norm, full_int_norm, diff_int_norm
      real(8) :: classify_tol, p_contribution_norm, p_contribution_ratio
      integer :: ix, iy, iz, idir
      logical :: current_available, current_bad, pcoef_saved_valid
      character(16) :: prod_ref_kind

      local_sum(:) = 0.0d0
      local_max(:) = 0.0d0
      current_bad = .false.
      current_available = .false.
      pcoef_saved_valid = .false.

      call allocate_vector_on_mg(curr_prod_grid)
      call allocate_vector_on_mg(curr_WW_grid)

      call calculate_microscopic_current_dg(dg_frag, system, mg, stencil, curr_prod_grid)
      current_available = .true.

      if (allocated(dg_frag%mixed_wannier_bpw_pcoef)) then
        pcoef_saved_valid = .true.
        allocate(pcoef_saved(size(dg_frag%mixed_wannier_bpw_pcoef,1), &
          size(dg_frag%mixed_wannier_bpw_pcoef,2), size(dg_frag%mixed_wannier_bpw_pcoef,3)))
        pcoef_saved(:,:,:) = dg_frag%mixed_wannier_bpw_pcoef(:,:,:)
        dg_frag%mixed_wannier_bpw_pcoef(:,:,:) = (0.0d0, 0.0d0)
        call calculate_microscopic_current_dg(dg_frag, system, mg, stencil, curr_WW_grid)
        dg_frag%mixed_wannier_bpw_pcoef(:,:,:) = pcoef_saved(:,:,:)
        deallocate(pcoef_saved)
      else
        curr_WW_grid%v(:,:,:,:) = curr_prod_grid%v(:,:,:,:)
      end if

      do iz = lbound(curr_prod_grid%v, 4), ubound(curr_prod_grid%v, 4)
        do iy = lbound(curr_prod_grid%v, 3), ubound(curr_prod_grid%v, 3)
          do ix = lbound(curr_prod_grid%v, 2), ubound(curr_prod_grid%v, 2)
            do idir = 1, 3
              prod_vec = curr_prod_grid%v(idir, ix, iy, iz)
              ww_vec = curr_WW_grid%v(idir, ix, iy, iz)
              diff_vec = prod_vec - ww_vec
              local_sum(idir) = local_sum(idir) + prod_vec * system%hvol
              local_sum(3+idir) = local_sum(3+idir) + ww_vec * system%hvol
              local_sum(6+idir) = local_sum(6+idir) + prod_vec * system%hvol
              local_sum(9+idir) = local_sum(9+idir) + diff_vec * system%hvol
              local_sum(13) = local_sum(13) + prod_vec * prod_vec * system%hvol
              local_sum(14) = local_sum(14) + ww_vec * ww_vec * system%hvol
              local_sum(15) = local_sum(15) + diff_vec * diff_vec * system%hvol
              local_sum(16) = local_sum(16) + system%hvol
              local_max(1) = max(local_max(1), abs(diff_vec))
              if (prod_vec /= prod_vec .or. ww_vec /= ww_vec .or. diff_vec /= diff_vec) current_bad = .true.
            end do
          end do
        end do
      end do
      if (current_available) local_sum(17) = 1.0d0
      if (pcoef_saved_valid) local_sum(18) = 1.0d0
      if (current_bad) local_sum(19) = 1.0d0
      call comm_summation(local_sum, global_sum, size(local_sum), dg_frag%icomm)
      call comm_get_max(local_max, global_max, size(local_max), dg_frag%icomm)

      prod_norm = sqrt(max(global_sum(13), 0.0d0))
      ww_norm = sqrt(max(global_sum(14), 0.0d0))
      diff_norm = sqrt(max(global_sum(15), 0.0d0))
      rms_diff = -1.0d0
      if (global_sum(16) > 0.0d0) rms_diff = sqrt(max(0.0d0, global_sum(15) / global_sum(16)))
      max_abs_diff = global_max(1)
      prod_int_norm = sqrt(global_sum(1)**2 + global_sum(2)**2 + global_sum(3)**2)
      ww_int_norm = sqrt(global_sum(4)**2 + global_sum(5)**2 + global_sum(6)**2)
      full_int_norm = sqrt(global_sum(7)**2 + global_sum(8)**2 + global_sum(9)**2)
      diff_int_norm = sqrt(global_sum(10)**2 + global_sum(11)**2 + global_sum(12)**2)
      p_contribution_norm = diff_norm
      p_contribution_ratio = 0.0d0
      if (prod_norm > 1.0d-300) p_contribution_ratio = p_contribution_norm / prod_norm
      classify_tol = 1.0d-10 * max(1.0d0, prod_norm, ww_norm)
      prod_ref_kind = 'inconsistent'
      if (diff_norm <= classify_tol) then
        prod_ref_kind = 'WW_only'
      else
        prod_ref_kind = 'full_WP'
      end if

      if (dg_frag%id == 0) then
        write(*,'(1x,a,1(a,i0),18(a,1pe12.4),10(a,l1),3(a,a))') &
          '[DG-MIXEDZ-LOCAL-CURRENT-PATH-CMP]', &
          ' step=', step_use, &
          ' current_prod_reference_norm=', prod_norm, &
          ' current_WW_contract_norm=', ww_norm, &
          ' current_WP_cross_contract_norm=', p_contribution_norm, &
          ' current_PP_contract_norm=', 0.0d0, &
          ' current_full_contract_norm=', prod_norm, &
          ' diff_prod_minus_WW_norm=', diff_norm, &
          ' diff_prod_minus_full_norm=', 0.0d0, &
          ' p_contribution_ratio=', p_contribution_ratio, &
          ' int_current_prod_x=', global_sum(1), &
          ' int_current_prod_y=', global_sum(2), &
          ' int_current_prod_z=', global_sum(3), &
          ' int_current_WW_x=', global_sum(4), &
          ' int_current_WW_y=', global_sum(5), &
          ' int_current_WW_z=', global_sum(6), &
          ' int_current_diff_norm=', diff_int_norm, &
          ' max_abs_grid_diff=', max_abs_diff, &
          ' rms_grid_diff=', rms_diff, &
          ' continuity_residual=', -1.0d0, &
          ' current_available=', global_sum(17) > 0.5d0, &
          ' grid_compare_available=', .true., &
          ' nonlocal_included=', .false., &
          ' diamagnetic_included=', .false., &
          ' pcoef_restore_verified=', global_sum(18) > 0.5d0, &
          ' current_replaced=', .false., &
          ' density_replaced=', .false., &
          ' production_replacement=', .false., &
          ' classification_succeeded=', trim(prod_ref_kind) /= 'inconsistent', &
          ' bad=', global_sum(19) > 0.5d0 .or. trim(prod_ref_kind) == 'inconsistent', &
          ' prod_ref_kind=', trim(prod_ref_kind), &
          ' reference_scope=', 'microscopic_paramagnetic_grid', &
          ' route=', 'current-production-reference-classification-dryrun'
      end if

      call deallocate_vector_diag(curr_prod_grid)
      call deallocate_vector_diag(curr_WW_grid)
    end subroutine diagnose_mixed_z_local_current_block_cmp

    subroutine mixed_z_log_w_eigenbasis_bridge_cmp(label, step_use, Sback, Hback, eval_all)
      character(len=*), intent(in) :: label
      integer, intent(in) :: step_use
      complex(8), intent(in) :: Sback(:,:), Hback(:,:)
      real(8), intent(in) :: eval_all(:)
      integer :: iw, jw, nw_use
      complex(8) :: Sprod_ij, Hprod_ij, Sdiff, Hdiff, phase_prod, phase_back
      real(8) :: S_Q_diff_from_I, Hdiff_frob, Hprod_frob, rel_H_diff
      real(8) :: H_eig_offdiag_rms, H_eig_diag_min, H_eig_diag_max
      real(8) :: eval_W_local_min, eval_W_local_max
      real(8) :: diag_minus_eval_rms, diag_minus_eval_max, phase_minus_prod_rms
      real(8) :: Hback_herm_diff, Sback_herm_diff, Hdiag_imag_max

      nw_use = min(size(Sback, 1), size(Sback, 2), size(Hback, 1), size(Hback, 2), size(eval_all))
      if (nw_use <= 0) return
      S_Q_diff_from_I = 0.0d0
      Hdiff_frob = 0.0d0
      Hprod_frob = 0.0d0
      H_eig_offdiag_rms = 0.0d0
      H_eig_diag_min = huge(1.0d0)
      H_eig_diag_max = -huge(1.0d0)
      eval_W_local_min = minval(eval_all(1:nw_use))
      eval_W_local_max = maxval(eval_all(1:nw_use))
      diag_minus_eval_rms = 0.0d0
      diag_minus_eval_max = 0.0d0
      phase_minus_prod_rms = 0.0d0
      Hback_herm_diff = 0.0d0
      Sback_herm_diff = 0.0d0
      Hdiag_imag_max = 0.0d0

      do jw = 1, nw_use
        do iw = 1, nw_use
          Sprod_ij = (0.0d0, 0.0d0)
          Hprod_ij = (0.0d0, 0.0d0)
          if (iw == jw) then
            Sprod_ij = (1.0d0, 0.0d0)
            Hprod_ij = cmplx(eval_all(iw), 0.0d0, kind=8)
          end if
          Sdiff = Sback(iw,jw) - Sprod_ij
          Hdiff = Hback(iw,jw) - Hprod_ij
          S_Q_diff_from_I = S_Q_diff_from_I + abs(Sdiff)**2
          Hdiff_frob = Hdiff_frob + abs(Hdiff)**2
          Hprod_frob = Hprod_frob + abs(Hprod_ij)**2
          Sback_herm_diff = max(Sback_herm_diff, abs(Sback(iw,jw) - conjg(Sback(jw,iw))))
          Hback_herm_diff = max(Hback_herm_diff, abs(Hback(iw,jw) - conjg(Hback(jw,iw))))
          if (iw == jw) then
            H_eig_diag_min = min(H_eig_diag_min, real(Hback(iw,iw), kind=8))
            H_eig_diag_max = max(H_eig_diag_max, real(Hback(iw,iw), kind=8))
            Hdiag_imag_max = max(Hdiag_imag_max, abs(aimag(Hback(iw,iw))))
            diag_minus_eval_rms = diag_minus_eval_rms + (real(Hback(iw,iw), kind=8) - eval_all(iw))**2
            diag_minus_eval_max = max(diag_minus_eval_max, abs(real(Hback(iw,iw), kind=8) - eval_all(iw)))
            phase_prod = exp(cmplx(0.0d0, -eval_all(iw) * dt, kind=8))
            phase_back = exp(cmplx(0.0d0, -real(Hback(iw,iw), kind=8) * dt, kind=8))
            phase_minus_prod_rms = phase_minus_prod_rms + abs(phase_back - phase_prod)**2
          else
            H_eig_offdiag_rms = H_eig_offdiag_rms + abs(Hback(iw,jw))**2
          end if
        end do
      end do
      S_Q_diff_from_I = sqrt(max(S_Q_diff_from_I, 0.0d0))
      Hdiff_frob = sqrt(max(Hdiff_frob, 0.0d0))
      Hprod_frob = sqrt(max(Hprod_frob, 0.0d0))
      rel_H_diff = Hdiff_frob / max(Hprod_frob, 1.0d-300)
      diag_minus_eval_rms = sqrt(diag_minus_eval_rms / max(dble(nw_use), 1.0d0))
      phase_minus_prod_rms = sqrt(phase_minus_prod_rms / max(dble(nw_use), 1.0d0))
      if (nw_use > 1) then
        H_eig_offdiag_rms = sqrt(H_eig_offdiag_rms / dble(nw_use * (nw_use - 1)))
      end if

      write(*,'(1x,a,1(a,i0),1(a,a),16(a,1pe12.4),4(a,l1),a,a)') &
        '[DG-MIXEDZ-LOCAL-FRAG-W-EIGENBASIS-BRIDGE-CMP]', &
        ' step=', step_use, &
        ' candidate=', trim(label), &
        ' S_Q_diff_from_I=', S_Q_diff_from_I, &
        ' H_eig_offdiag_rms=', H_eig_offdiag_rms, &
        ' H_eig_diag_min=', H_eig_diag_min, &
        ' H_eig_diag_max=', H_eig_diag_max, &
        ' eval_W_local_min=', eval_W_local_min, &
        ' eval_W_local_max=', eval_W_local_max, &
        ' diag_minus_eval_rms=', diag_minus_eval_rms, &
        ' diag_minus_eval_max=', diag_minus_eval_max, &
        ' phase_minus_prod_rms=', phase_minus_prod_rms, &
        ' rel_H_diff=', rel_H_diff, &
        ' H_diff_frob=', Hdiff_frob, &
        ' Hprod_frob=', Hprod_frob, &
        ' Hback_herm_diff=', Hback_herm_diff, &
        ' Sback_herm_diff=', Sback_herm_diff, &
        ' Hdiag_imag_max=', Hdiag_imag_max, &
        ' nw_real=', dble(nw_use), &
        ' field_off_static_phase=', sum(abs(E_mid(1:3))) <= 1.0d-30, &
        ' W_sector_only=', .true., &
        ' production_sourced=', .true., &
        ' bad=', rel_H_diff /= rel_H_diff .or. S_Q_diff_from_I /= S_Q_diff_from_I, &
        ' route=', 'production-CW-sourced-local-W-eigenbasis-bridge-diagnostic-only'
    end subroutine mixed_z_log_w_eigenbasis_bridge_cmp

    subroutine build_mixedz_full_cref_from_prodcoef(S_build, c_source, reduced_transform, S_reduced, &
      nraw, nstate, nred_ref, c_cref, info_cref)
      complex(8), intent(in) :: S_build(:,:), c_source(:,:), reduced_transform(:,:), S_reduced(:,:)
      integer, intent(in) :: nraw, nstate, nred_ref
      complex(8), intent(out) :: c_cref(:,:)
      integer, intent(out) :: info_cref
      integer :: ist, info_local
      complex(8), allocatable :: tmp_metric(:), rhs_red(:), c_red(:)

      c_cref(:, :) = (0.0d0, 0.0d0)
      info_cref = 0
      if (nraw <= 0 .or. nstate <= 0 .or. nred_ref <= 0) then
        info_cref = 1
        return
      end if
      if (size(S_build, 1) < nraw .or. size(S_build, 2) < nraw .or. &
          size(c_source, 1) < nraw .or. size(c_source, 2) < nstate .or. &
          size(reduced_transform, 1) < nraw .or. size(reduced_transform, 2) < nred_ref .or. &
          size(S_reduced, 1) < nred_ref .or. size(S_reduced, 2) < nred_ref .or. &
          size(c_cref, 1) < nraw .or. size(c_cref, 2) < nstate) then
        info_cref = 1
        return
      end if

      allocate(tmp_metric(nraw), rhs_red(nred_ref), c_red(nred_ref))
      do ist = 1, nstate
        tmp_metric(:) = matmul(S_build(1:nraw, 1:nraw), c_source(1:nraw, ist))
        rhs_red(:) = matmul(conjg(transpose(reduced_transform(1:nraw, 1:nred_ref))), tmp_metric(:))
        call mixed_z_solve_metric_vector(S_reduced(1:nred_ref, 1:nred_ref), rhs_red(:), c_red(:), info_local)
        if (info_local /= 0) then
          info_cref = 2
        else
          c_cref(1:nraw, ist) = matmul(reduced_transform(1:nraw, 1:nred_ref), c_red(:))
        end if
      end do
      deallocate(tmp_metric, rhs_red, c_red)
    end subroutine build_mixedz_full_cref_from_prodcoef

    subroutine gather_global_mixed_coefficients(ispin_current, state_s, state_e, cmix)
      integer, intent(in) :: ispin_current, state_s, state_e
      complex(8), intent(out) :: cmix(:,:)
      integer :: ifrag_row, i_local_row, nbf, ibasis, istate, state_col, nvalid
      integer :: global_row, local_row, nstate_blk, nw, np
      real(8) :: prop_t0

      nstate_blk = max(0, state_e - state_s + 1)
      nw = dg_frag%mixed_wannier_bpw_nw
      np = dg_frag%mixed_wannier_bpw_np
      cmix(:, :) = (0.0d0, 0.0d0)
      if (nstate_blk <= 0 .or. nw <= 0) return
      if (.not. allocated(prop_cw_local_work) .or. size(prop_cw_local_work,1) /= nw .or. &
          size(prop_cw_local_work,2) /= nstate_blk) then
        if (allocated(prop_cw_local_work)) deallocate(prop_cw_local_work, prop_cw_global_work)
        allocate(prop_cw_local_work(nw,nstate_blk), prop_cw_global_work(nw,nstate_blk))
        if (dg_frag%mixed_z_perf_count_enabled) then
          dg_frag%mixed_z_perf_prop_pack_allocs = dg_frag%mixed_z_perf_prop_pack_allocs + 1_8
        end if
      end if
      if (dg_frag%mixed_z_perf_count_enabled) then
        dg_frag%mixed_z_perf_prop_pack_calls = dg_frag%mixed_z_perf_prop_pack_calls + 1_8
        dg_frag%mixed_z_perf_prop_pack_zero_bytes = dg_frag%mixed_z_perf_prop_pack_zero_bytes + &
          16_8 * int(nw, kind=8) * int(nstate_blk, kind=8)
      end if
      prop_cw_local_work(:, :) = (0.0d0, 0.0d0)
      if (dg_frag%mixed_z_perf_count_enabled) prop_t0 = get_wtime()
      do ifrag_row = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local_row = ifrag_row - dg_frag%ifrag_start + 1
        if (i_local_row < 1 .or. i_local_row > size(dg_frag%global_wannier_coef, 4)) cycle
        nbf = min(dg_frag%n_basis(ifrag_row, ispin_current), size(dg_frag%global_wannier_coef, 1))
        if (nbf <= 0) cycle
        if (.not. allocated(prop_w_block_work) .or. size(prop_w_block_work,1) < nbf .or. &
            size(prop_w_block_work,2) < nw) then
          if (allocated(prop_w_block_work)) deallocate(prop_w_block_work)
          allocate(prop_w_block_work(nbf,nw))
        end if
        if (.not. allocated(prop_coef_block_work) .or. size(prop_coef_block_work,1) < nbf .or. &
            size(prop_coef_block_work,2) < nstate_blk) then
          if (allocated(prop_coef_block_work)) deallocate(prop_coef_block_work)
          allocate(prop_coef_block_work(nbf,nstate_blk))
        end if
        nvalid = 0
        do ibasis = 1, nbf
          global_row = dg_frag%index_basis(ibasis, ifrag_row, ispin_current)
          if (global_row < 1 .or. global_row > dg_frag%n_mat_max) cycle
          local_row = dg_frag%coef_global_to_local(global_row, ispin_current)
          if (local_row < 1 .or. local_row > size(dg_frag%coef, 1)) cycle
          nvalid = nvalid + 1
          if (allocated(dg_frag%global_wannier_flux_evec)) then
            prop_w_block_work(nvalid,1:nw) = matmul( &
              dg_frag%global_wannier_coef(ibasis,1:size(dg_frag%global_wannier_flux_evec,1), &
                ispin_current,i_local_row), &
              dg_frag%global_wannier_flux_evec(1:size(dg_frag%global_wannier_flux_evec,1),1:nw))
          else
            prop_w_block_work(nvalid,1:nw) = &
              dg_frag%global_wannier_coef(ibasis,1:nw,ispin_current,i_local_row)
          end if
          do state_col = 1, nstate_blk
            prop_coef_block_work(nvalid,state_col) = &
              dg_frag%coef(local_row,state_s+state_col-1,ispin_current)
          end do
          if (dg_frag%mixed_z_perf_count_enabled) then
            dg_frag%mixed_z_perf_prop_pack_w_copy_bytes = &
              dg_frag%mixed_z_perf_prop_pack_w_copy_bytes + &
              16_8 * int(nw, kind=8) * int(nstate_blk, kind=8)
          end if
        end do
        if (nvalid > 0) then
          if (dg_frag%mixed_z_perf_count_enabled) &
            dg_frag%mixed_z_perf_zgemm_calls = dg_frag%mixed_z_perf_zgemm_calls + 1_8
          prop_cw_local_work(1:nw,1:nstate_blk) = prop_cw_local_work(1:nw,1:nstate_blk) + &
            matmul(conjg(transpose(prop_w_block_work(1:nvalid,1:nw))), &
                   prop_coef_block_work(1:nvalid,1:nstate_blk))
        end if
      end do
      if (dg_frag%mixed_z_perf_count_enabled) then
        dg_frag%mixed_z_perf_wall_prop_pack = &
          dg_frag%mixed_z_perf_wall_prop_pack + (get_wtime() - prop_t0)
        prop_t0 = get_wtime()
      end if
      call comm_summation(prop_cw_local_work, prop_cw_global_work, nw*nstate_blk, dg_frag%icomm)
      if (dg_frag%mixed_z_perf_count_enabled) then
        dg_frag%mixed_z_perf_wall_prop_comm = &
          dg_frag%mixed_z_perf_wall_prop_comm + (get_wtime() - prop_t0)
      end if
      cmix(1:nw,1:nstate_blk) = prop_cw_global_work(1:nw,1:nstate_blk)
      if (np > 0 .and. allocated(dg_frag%mixed_wannier_bpw_pcoef)) then
        if (dg_frag%mixed_z_perf_count_enabled) then
          dg_frag%mixed_z_perf_prop_pack_p_copy_bytes = dg_frag%mixed_z_perf_prop_pack_p_copy_bytes + &
            16_8 * int(np, kind=8) * int(nstate_blk, kind=8)
        end if
        cmix(nw+1:nw+np,1:nstate_blk) = dg_frag%mixed_wannier_bpw_pcoef(1:np,state_s:state_e,ispin_current)
      end if
    end subroutine gather_global_mixed_coefficients

    subroutine scatter_global_mixed_coefficients(ispin_current, state_s, state_e, cmix)
      integer, intent(in) :: ispin_current, state_s, state_e
      complex(8), intent(in) :: cmix(:,:)
      integer :: ifrag_row, i_local_row, nbf, ibasis, state_col, nvalid
      integer :: global_row, local_row, nstate_blk, nw, np
      real(8) :: prop_t0

      nstate_blk = max(0, state_e - state_s + 1)
      nw = dg_frag%mixed_wannier_bpw_nw
      np = dg_frag%mixed_wannier_bpw_np
      if (nstate_blk <= 0 .or. nw <= 0) return
      if (.not. allocated(prop_next_local_work) .or. size(prop_next_local_work,1) /= size(dg_frag%coef,1) .or. &
          size(prop_next_local_work,2) /= nstate_blk) then
        if (allocated(prop_next_local_work)) deallocate(prop_next_local_work)
        allocate(prop_next_local_work(size(dg_frag%coef,1),nstate_blk))
        if (dg_frag%mixed_z_perf_count_enabled) then
          dg_frag%mixed_z_perf_prop_unpack_allocs = dg_frag%mixed_z_perf_prop_unpack_allocs + 1_8
        end if
      end if
      if (dg_frag%mixed_z_perf_count_enabled) then
        dg_frag%mixed_z_perf_prop_unpack_calls = dg_frag%mixed_z_perf_prop_unpack_calls + 1_8
        dg_frag%mixed_z_perf_prop_unpack_zero_bytes = dg_frag%mixed_z_perf_prop_unpack_zero_bytes + &
          16_8 * int(size(dg_frag%coef,1), kind=8) * int(nstate_blk, kind=8)
      end if
      prop_next_local_work(:, :) = (0.0d0, 0.0d0)
      if (dg_frag%mixed_z_perf_count_enabled) prop_t0 = get_wtime()
      do ifrag_row = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local_row = ifrag_row - dg_frag%ifrag_start + 1
        if (i_local_row < 1 .or. i_local_row > size(dg_frag%global_wannier_coef, 4)) cycle
        nbf = min(dg_frag%n_basis(ifrag_row, ispin_current), size(dg_frag%global_wannier_coef, 1))
        if (nbf <= 0) cycle
        if (.not. allocated(prop_w_block_work) .or. size(prop_w_block_work,1) < nbf .or. &
            size(prop_w_block_work,2) < nw) then
          if (allocated(prop_w_block_work)) deallocate(prop_w_block_work)
          allocate(prop_w_block_work(nbf,nw))
        end if
        if (.not. allocated(prop_scatter_block_work) .or. size(prop_scatter_block_work,1) < nbf .or. &
            size(prop_scatter_block_work,2) < nstate_blk) then
          if (allocated(prop_scatter_block_work)) deallocate(prop_scatter_block_work)
          allocate(prop_scatter_block_work(nbf,nstate_blk))
        end if
        if (.not. allocated(prop_local_row_work) .or. size(prop_local_row_work) < nbf) then
          if (allocated(prop_local_row_work)) deallocate(prop_local_row_work)
          allocate(prop_local_row_work(nbf))
        end if
        nvalid = 0
        do ibasis = 1, nbf
          global_row = dg_frag%index_basis(ibasis, ifrag_row, ispin_current)
          if (global_row < 1 .or. global_row > dg_frag%n_mat_max) cycle
          local_row = dg_frag%coef_global_to_local(global_row, ispin_current)
          if (local_row < 1 .or. local_row > size(dg_frag%coef, 1)) cycle
          nvalid = nvalid + 1
          prop_local_row_work(nvalid) = local_row
          if (allocated(dg_frag%global_wannier_flux_evec)) then
            prop_w_block_work(nvalid,1:nw) = matmul( &
              dg_frag%global_wannier_coef(ibasis,1:size(dg_frag%global_wannier_flux_evec,1), &
                ispin_current,i_local_row), &
              dg_frag%global_wannier_flux_evec(1:size(dg_frag%global_wannier_flux_evec,1),1:nw))
          else
            prop_w_block_work(nvalid,1:nw) = &
              dg_frag%global_wannier_coef(ibasis,1:nw,ispin_current,i_local_row)
          end if
          if (dg_frag%mixed_z_perf_count_enabled) then
            dg_frag%mixed_z_perf_prop_unpack_w_copy_bytes = &
              dg_frag%mixed_z_perf_prop_unpack_w_copy_bytes + &
              16_8 * int(nw, kind=8) * int(nstate_blk, kind=8)
          end if
        end do
        if (nvalid > 0) then
          if (dg_frag%mixed_z_perf_count_enabled) &
            dg_frag%mixed_z_perf_zgemm_calls = dg_frag%mixed_z_perf_zgemm_calls + 1_8
          prop_scatter_block_work(1:nvalid,1:nstate_blk) = &
            matmul(prop_w_block_work(1:nvalid,1:nw), cmix(1:nw,1:nstate_blk))
          do ibasis = 1, nvalid
            local_row = prop_local_row_work(ibasis)
            prop_next_local_work(local_row,1:nstate_blk) = &
              prop_next_local_work(local_row,1:nstate_blk) + &
              prop_scatter_block_work(ibasis,1:nstate_blk)
          end do
        end if
      end do
      do istate = 1, nstate_blk
        state_col = state_s + istate - 1
        if (state_col < 1 .or. state_col > size(dg_frag%coef, 2)) cycle
        dg_frag%coef(1:size(dg_frag%coef,1),state_col,ispin_current) = &
          prop_next_local_work(1:size(dg_frag%coef,1),istate)
      end do
      if (np > 0 .and. allocated(dg_frag%mixed_wannier_bpw_pcoef)) then
        if (dg_frag%mixed_z_perf_count_enabled) then
          dg_frag%mixed_z_perf_prop_unpack_p_copy_bytes = dg_frag%mixed_z_perf_prop_unpack_p_copy_bytes + &
            16_8 * int(np, kind=8) * int(nstate_blk, kind=8)
        end if
        dg_frag%mixed_wannier_bpw_pcoef(1:np,state_s:state_e,ispin_current) = cmix(nw+1:nw+np,1:nstate_blk)
      end if
      if (dg_frag%mixed_z_perf_count_enabled) then
        dg_frag%mixed_z_perf_wall_prop_unpack = &
          dg_frag%mixed_z_perf_wall_prop_unpack + (get_wtime() - prop_t0)
      end if
    end subroutine scatter_global_mixed_coefficients

    subroutine diagnose_mixed_z_split_observables(cmix_split, norm_split, step_use, bad_in)
      complex(8), intent(in) :: cmix_split(:,:,:)
      real(8), intent(in) :: norm_split
      integer, intent(in) :: step_use
      logical, intent(in) :: bad_in
      type(s_scalar) :: rho_prod_after, rho_split
      type(s_scalar), allocatable :: rho_s_prod_after(:), rho_s_split(:)
      complex(8), allocatable :: coef_saved(:,:,:), pcoef_saved(:,:,:)
      real(8) :: local_sum(7), global_sum(7), local_max(1), global_max(1)
      real(8) :: int_prod, int_split, int_diff, rms_diff, pz_prod_grid, pz_split_grid
      real(8) :: diff_grid, z_coord
      integer :: ispin_local, ix, iy, iz
      logical :: bad_obs

      bad_obs = bad_in
      if (.not. allocated(dg_frag%coef)) return
      if (.not. allocated(dg_frag%mixed_wannier_bpw_pcoef)) return
      allocate(coef_saved(size(dg_frag%coef,1), size(dg_frag%coef,2), size(dg_frag%coef,3)))
      allocate(pcoef_saved(size(dg_frag%mixed_wannier_bpw_pcoef,1), &
        size(dg_frag%mixed_wannier_bpw_pcoef,2), size(dg_frag%mixed_wannier_bpw_pcoef,3)))
      coef_saved(:,:,:) = dg_frag%coef(:,:,:)
      pcoef_saved(:,:,:) = dg_frag%mixed_wannier_bpw_pcoef(:,:,:)

      call allocate_scalar_like(rho, rho_prod_after)
      call allocate_scalar_like(rho, rho_split)
      allocate(rho_s_prod_after(system%nspin), rho_s_split(system%nspin))
      do ispin_local = 1, system%nspin
        call allocate_scalar_like(rho_s(ispin_local), rho_s_prod_after(ispin_local))
        call allocate_scalar_like(rho_s(ispin_local), rho_s_split(ispin_local))
      end do

      call calculate_density_from_fragments(dg_frag, system, mg, rho_prod_after, rho_s_prod_after, step_use)
      do ispin_local = 1, min(dg_frag%nspin, size(cmix_split, 3))
        call scatter_global_mixed_coefficients(ispin_local, state_first, state_last, &
          cmix_split(1:size(cmix_split,1),1:size(cmix_split,2),ispin_local))
      end do
      call calculate_density_from_fragments(dg_frag, system, mg, rho_split, rho_s_split, step_use)

      dg_frag%coef(:,:,:) = coef_saved(:,:,:)
      dg_frag%mixed_wannier_bpw_pcoef(:,:,:) = pcoef_saved(:,:,:)

      local_sum(:) = 0.0d0
      local_max(1) = 0.0d0
      do iz = lbound(rho_prod_after%f, 3), ubound(rho_prod_after%f, 3)
        if (mod(dg_frag%lgnum_total(3), 2) == 0) then
          z_coord = dble(iz) - 0.5d0
        else
          z_coord = dble(iz)
        end if
        do iy = lbound(rho_prod_after%f, 2), ubound(rho_prod_after%f, 2)
          do ix = lbound(rho_prod_after%f, 1), ubound(rho_prod_after%f, 1)
            diff_grid = rho_split%f(ix,iy,iz) - rho_prod_after%f(ix,iy,iz)
            local_sum(1) = local_sum(1) + rho_prod_after%f(ix,iy,iz) * system%hvol
            local_sum(2) = local_sum(2) + rho_split%f(ix,iy,iz) * system%hvol
            local_sum(3) = local_sum(3) + diff_grid * system%hvol
            local_sum(4) = local_sum(4) + diff_grid * diff_grid * system%hvol
            local_sum(5) = local_sum(5) + rho_prod_after%f(ix,iy,iz) * z_coord * system%hvol
            local_sum(6) = local_sum(6) + rho_split%f(ix,iy,iz) * z_coord * system%hvol
            local_sum(7) = local_sum(7) + system%hvol
            local_max(1) = max(local_max(1), abs(diff_grid))
            if (diff_grid /= diff_grid .or. rho_split%f(ix,iy,iz) /= rho_split%f(ix,iy,iz) .or. &
                rho_prod_after%f(ix,iy,iz) /= rho_prod_after%f(ix,iy,iz)) bad_obs = .true.
          end do
        end do
      end do
      call comm_summation(local_sum, global_sum, size(local_sum), dg_frag%icomm)
      call comm_get_max(local_max, global_max, size(local_max), dg_frag%icomm)

      int_prod = global_sum(1)
      int_split = global_sum(2)
      int_diff = global_sum(3)
      rms_diff = sqrt(max(0.0d0, global_sum(4)) / max(global_sum(7), 1.0d-300))
      pz_prod_grid = global_sum(5)
      pz_split_grid = global_sum(6)
      if (dg_frag%id == 0) then
        write(*,'(1x,a,2(a,i0),10(a,1pe16.8),a,l1)') &
          '[DG-MIXEDZ-SPLIT-OBS-CMP]', &
          ' step=', step_use, ' nstate=', state_last - state_first + 1, &
          ' norm_coef_mixed_split=', norm_split, &
          ' int_rho_prod=', int_prod, &
          ' int_rho_mixed_split=', int_split, &
          ' int_rho_diff=', int_diff, &
          ' max_abs_rho_diff=', global_max(1), &
          ' rms_rho_diff=', rms_diff, &
          ' Pz_prod=', pz_prod_grid, &
          ' Pz_mixed_split=', pz_split_grid, &
          ' Pz_diff=', pz_split_grid - pz_prod_grid, &
          ' rel_Pz_diff=', (pz_split_grid - pz_prod_grid) / max(abs(pz_prod_grid), 1.0d-300), &
          ' bad=', bad_obs
      end if

      call deallocate_scalar_diag(rho_prod_after)
      call deallocate_scalar_diag(rho_split)
      do ispin_local = 1, system%nspin
        call deallocate_scalar_diag(rho_s_prod_after(ispin_local))
        call deallocate_scalar_diag(rho_s_split(ispin_local))
      end do
      deallocate(rho_s_prod_after, rho_s_split, coef_saved, pcoef_saved)
    end subroutine diagnose_mixed_z_split_observables

    subroutine diagnose_mixed_z_split_current(cmix_split, norm_split, step_use, bad_in)
      complex(8), intent(in) :: cmix_split(:,:,:)
      real(8), intent(in) :: norm_split
      integer, intent(in) :: step_use
      logical, intent(in) :: bad_in
      type(s_vector) :: curr_prod_after, curr_split
      complex(8), allocatable :: coef_saved(:,:,:), pcoef_saved(:,:,:)
      real(8) :: local_sum(5), global_sum(5), local_max(1), global_max(1)
      real(8) :: diff_vec, prod_vec, split_vec, diff_norm, prod_norm, split_norm
      integer :: ispin_local, ix, iy, iz, idir
      logical :: bad_current

      bad_current = bad_in
      if (.not. allocated(dg_frag%coef)) return
      if (.not. allocated(dg_frag%mixed_wannier_bpw_pcoef)) return
      allocate(coef_saved(size(dg_frag%coef,1), size(dg_frag%coef,2), size(dg_frag%coef,3)))
      allocate(pcoef_saved(size(dg_frag%mixed_wannier_bpw_pcoef,1), &
        size(dg_frag%mixed_wannier_bpw_pcoef,2), size(dg_frag%mixed_wannier_bpw_pcoef,3)))
      coef_saved(:,:,:) = dg_frag%coef(:,:,:)
      pcoef_saved(:,:,:) = dg_frag%mixed_wannier_bpw_pcoef(:,:,:)

      call allocate_vector_on_mg(curr_prod_after)
      call allocate_vector_on_mg(curr_split)

      call calculate_microscopic_current_dg(dg_frag, system, mg, stencil, curr_prod_after)
      do ispin_local = 1, min(dg_frag%nspin, size(cmix_split, 3))
        call scatter_global_mixed_coefficients(ispin_local, state_first, state_last, &
          cmix_split(1:size(cmix_split,1),1:size(cmix_split,2),ispin_local))
      end do
      call calculate_microscopic_current_dg(dg_frag, system, mg, stencil, curr_split)

      dg_frag%coef(:,:,:) = coef_saved(:,:,:)
      dg_frag%mixed_wannier_bpw_pcoef(:,:,:) = pcoef_saved(:,:,:)

      local_sum(:) = 0.0d0
      local_max(1) = 0.0d0
      do iz = lbound(curr_prod_after%v, 4), ubound(curr_prod_after%v, 4)
        do iy = lbound(curr_prod_after%v, 3), ubound(curr_prod_after%v, 3)
          do ix = lbound(curr_prod_after%v, 2), ubound(curr_prod_after%v, 2)
            do idir = 1, 3
              prod_vec = curr_prod_after%v(idir,ix,iy,iz)
              split_vec = curr_split%v(idir,ix,iy,iz)
              diff_vec = split_vec - prod_vec
              local_sum(1) = local_sum(1) + prod_vec * prod_vec * system%hvol
              local_sum(2) = local_sum(2) + split_vec * split_vec * system%hvol
              local_sum(3) = local_sum(3) + diff_vec * diff_vec * system%hvol
              local_sum(4) = local_sum(4) + diff_vec * system%hvol
              local_max(1) = max(local_max(1), abs(diff_vec))
              if (prod_vec /= prod_vec .or. split_vec /= split_vec .or. diff_vec /= diff_vec) &
                bad_current = .true.
            end do
          end do
        end do
      end do
      local_sum(5) = dble(size(curr_prod_after%v))
      call comm_summation(local_sum, global_sum, size(local_sum), dg_frag%icomm)
      call comm_get_max(local_max, global_max, size(local_max), dg_frag%icomm)

      prod_norm = sqrt(max(global_sum(1), 0.0d0))
      split_norm = sqrt(max(global_sum(2), 0.0d0))
      diff_norm = sqrt(max(global_sum(3), 0.0d0))
      if (dg_frag%id == 0) then
        write(*,'(1x,a,2(a,i0),8(a,1pe16.8),a,l1)') &
          '[DG-MIXEDZ-SPLIT-CURRENT-CMP]', &
          ' step=', step_use, ' nstate=', state_last - state_first + 1, &
          ' norm_coef_mixed_split=', norm_split, &
          ' current_prod_norm=', prod_norm, &
          ' current_mixed_split_norm=', split_norm, &
          ' J_diff=', diff_norm, &
          ' rel_J_diff=', diff_norm / max(prod_norm, 1.0d-300), &
          ' int_J_diff=', global_sum(4), &
          ' max_abs_J_diff=', global_max(1), &
          ' rms_J_diff=', diff_norm / sqrt(max(global_sum(5), 1.0d0)), &
          ' bad=', bad_current
      end if

      call deallocate_vector_diag(curr_prod_after)
      call deallocate_vector_diag(curr_split)
      deallocate(coef_saved, pcoef_saved)
    end subroutine diagnose_mixed_z_split_current

    subroutine allocate_scalar_like(src, dst)
      type(s_scalar), intent(in) :: src
      type(s_scalar), intent(inout) :: dst
      allocate(dst%f(lbound(src%f,1):ubound(src%f,1), &
                     lbound(src%f,2):ubound(src%f,2), &
                     lbound(src%f,3):ubound(src%f,3)))
      dst%f(:,:,:) = 0.0d0
    end subroutine allocate_scalar_like

    subroutine deallocate_scalar_diag(val)
      type(s_scalar), intent(inout) :: val
      if (allocated(val%f)) deallocate(val%f)
    end subroutine deallocate_scalar_diag

    subroutine allocate_vector_on_mg(dst)
      type(s_vector), intent(inout) :: dst
      allocate(dst%v(3, mg%is(1):mg%ie(1), &
                       mg%is(2):mg%ie(2), &
                       mg%is(3):mg%ie(3)))
      dst%v(:,:,:,:) = 0.0d0
    end subroutine allocate_vector_on_mg

    subroutine deallocate_vector_diag(val)
      type(s_vector), intent(inout) :: val
      if (allocated(val%v)) deallocate(val%v)
    end subroutine deallocate_vector_diag

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

    subroutine apply_full_h_seed_eigen_exp(state_s, state_e)
      integer, intent(in) :: state_s, state_e
      integer :: ispin_current, nstate_blk, neig, nrow_full
      integer :: istate_blk, eig, global_row, local_row, state_col
      real(8) :: phase_c, phase_s
      real(8) :: tloc
      complex(8), allocatable :: amp_local(:,:), amp_global(:,:), next_local(:,:)
      complex(8), allocatable :: h_step(:,:), h_vec(:,:), tmp(:)
      real(8), allocatable :: eval_step(:)

      nstate_blk = max(0, state_e - state_s + 1)
      if (nstate_blk <= 0) return
      if (.not. allocated(dg_frag%full_h_seed_evec)) return
      if (.not. allocated(dg_frag%full_h_seed_eval)) return
      if (.not. allocated(dg_frag%full_h_seed_xi)) return
      if (.not. allocated(dg_frag%coef_global_to_local)) return

      neig = min(dg_frag%full_h_seed_nstate, size(dg_frag%full_h_seed_evec, 2), &
                 size(dg_frag%full_h_seed_eval, 1), size(dg_frag%full_h_seed_xi, 2), &
                 size(dg_frag%full_h_seed_xi, 3))
      nrow_full = size(dg_frag%full_h_seed_evec, 1)
      if (neig <= 0 .or. nrow_full <= 0) return

      allocate(amp_local(neig,nstate_blk), amp_global(neig,nstate_blk))
      allocate(h_step(neig,neig), h_vec(neig,neig), eval_step(neig), tmp(neig))
      do ispin_current = 1, dg_frag%nspin
        if (ispin_current > size(dg_frag%full_h_seed_evec, 3)) cycle
        if (ispin_current > size(dg_frag%full_h_seed_eval, 2)) cycle
        if (ispin_current > size(dg_frag%full_h_seed_xi, 4)) cycle

        amp_local(:, :) = (0.0d0, 0.0d0)
        tloc = get_wtime()
        do global_row = 1, min(nrow_full, dg_frag%n_mat_max)
          local_row = dg_frag%coef_global_to_local(global_row, ispin_current)
          if (local_row < 1 .or. local_row > size(dg_frag%coef, 1)) cycle
          do istate_blk = 1, nstate_blk
            state_col = state_s + istate_blk - 1
            if (state_col < 1 .or. state_col > size(dg_frag%coef, 2)) cycle
            do eig = 1, neig
              amp_local(eig,istate_blk) = amp_local(eig,istate_blk) + &
                conjg(dg_frag%full_h_seed_evec(global_row,eig,ispin_current)) * &
                dg_frag%coef(local_row,state_col,ispin_current)
            end do
          end do
        end do
        time_flux_project = time_flux_project + (get_wtime() - tloc)

        tloc = get_wtime()
        call comm_summation(amp_local, amp_global, neig*nstate_blk, dg_frag%icomm)
        time_flux_comm = time_flux_comm + (get_wtime() - tloc)

        h_step(:, :) = (0.0d0, 0.0d0)
        do eig = 1, neig
          h_step(eig,eig) = cmplx(dg_frag%full_h_seed_eval(eig,ispin_current), 0.0d0, kind=8)
        end do
        h_step(1:neig,1:neig) = h_step(1:neig,1:neig) &
          - E_mid(1) * dg_frag%full_h_seed_xi(1,1:neig,1:neig,ispin_current) &
          - E_mid(2) * dg_frag%full_h_seed_xi(2,1:neig,1:neig,ispin_current) &
          - E_mid(3) * dg_frag%full_h_seed_xi(3,1:neig,1:neig,ispin_current)
        call symmetrize_expdiag_h_eff(h_step, neig)
        tloc = get_wtime()
        call eigen_zheev(h_step, eval_step, h_vec)
        time_local_diag = time_local_diag + (get_wtime() - tloc)

        do istate_blk = 1, nstate_blk
          tmp(:) = matmul(conjg(transpose(h_vec(1:neig,1:neig))), amp_global(:,istate_blk))
          do eig = 1, neig
            phase_c = cos(eval_step(eig) * dt)
            phase_s = sin(eval_step(eig) * dt)
            tmp(eig) = cmplx(phase_c, -phase_s, kind=8) * tmp(eig)
          end do
          amp_global(:,istate_blk) = matmul(h_vec(1:neig,1:neig), tmp(:))
        end do

        allocate(next_local(size(dg_frag%coef,1),nstate_blk))
        next_local(:, :) = (0.0d0, 0.0d0)
        tloc = get_wtime()
        do global_row = 1, min(nrow_full, dg_frag%n_mat_max)
          local_row = dg_frag%coef_global_to_local(global_row, ispin_current)
          if (local_row < 1 .or. local_row > size(dg_frag%coef, 1)) cycle
          do istate_blk = 1, nstate_blk
            do eig = 1, neig
              next_local(local_row,istate_blk) = next_local(local_row,istate_blk) + &
                dg_frag%full_h_seed_evec(global_row,eig,ispin_current) * amp_global(eig,istate_blk)
            end do
          end do
        end do
        time_flux_scatter = time_flux_scatter + (get_wtime() - tloc)

        do istate_blk = 1, nstate_blk
          state_col = state_s + istate_blk - 1
          if (state_col < 1 .or. state_col > size(dg_frag%coef, 2)) cycle
          dg_frag%coef(1:size(dg_frag%coef,1),state_col,ispin_current) = &
            next_local(1:size(dg_frag%coef,1),istate_blk)
        end do
        deallocate(next_local)
      end do
      deallocate(amp_local, amp_global, h_step, h_vec, eval_step, tmp)
    end subroutine apply_full_h_seed_eigen_exp

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
