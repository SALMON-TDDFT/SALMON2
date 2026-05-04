  subroutine calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s, itt_debug)
    use structures
    use salmon_global, only: nelec, nelec_spin
    use communication, only: comm_summation, comm_bcast, comm_alltoallv, comm_send, comm_recv, COMM_GROUP_NULL
    use rt_dg_fragment_ops, only: refresh_pw_coef_cache, gather_full_coef_view, copy_overlap_operator_to_dense, &
      apply_overlap_operator_batch, capture_occmap_pair_snapshot
    use rt_dg_fragment_types, only: density_grid_point_info
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_scalar),         intent(inout) :: rho
    type(s_scalar),         intent(inout) :: rho_s(system%nspin)
    integer, intent(in), optional :: itt_debug

    integer :: ifrag, io, i_local, ispin, iocc
    integer :: istate_frag
    integer :: ix, iy, iz, ixg, iyg, izg, bx, by, bz, owner_rank
    integer :: ig_i, nbf, nbf_max, ipw, n_pw, n_frag, n_tot, n_basis_mix, max_mixed_basis
    integer :: nxyz(3), ifrag_count, ngrid_max
    integer :: nocc_spin, nocc_cache
    integer :: irank, slot, npts, idx_local, idx_remote
    integer :: local_grid_count, remote_grid_count, valid_remote_grid_count
    integer :: igrid0, igrid, ngrid, npt_blk, io0, io1, nbatch, nstate, ipw0, npw_blk
    integer :: total_send_pts, subgroup_root_rank, block_idx_global
    integer :: send_total_count, recv_total_count
    integer :: handler_send_peers, handler_recv_peers
    integer :: nblocks_ifrag, first_block_offset, block_step_blocks, block_offset
    integer :: valid_basis_count
    integer :: owner_valid_count, slot0_count, slotp_count, owner_true_count, owner_false_count
    integer :: basis_gid_probe(3)
    integer :: handler_id_frag
    integer :: spin_offset
    integer :: ix0_frag, ix1_frag, iy0_frag, iy1_frag, iz0_frag, iz1_frag
    integer :: cmp_slot
    integer :: ixg_min_probe, ixg_max_probe, owner_valid_probe, nprobe_cols, iprobe
    integer, parameter :: grid_block_size = 1024, state_block_size = 64, rho_state_block_size = 16, pw_block_size = 128
    integer, parameter :: mixed_io_block_size = 64
    complex(8), parameter :: zzero = (0.0d0, 0.0d0), zone = (1.0d0, 0.0d0)
    real(8) :: phi_i, rho_contrib, rho_raw_contrib, rho_accum, rho_mix_accum
    real(8) :: total_charge, total_charge_local, scale_rho, rho_sum_local
    real(8) :: occ_factor
    real(8) :: phi_sample_probe
    real(8) :: boxL(3), inv_sqrt_vol, theta, inv_lgnum1
    real(8) :: t_total0, t_total1, t_cache0, t_cache1
    real(8) :: t_project0, t_project1, t_comm0, t_comm1, t_norm0, t_norm1
    real(8) :: t_setup0, t_setup1, t_psi0, t_psi1, t_rho0, t_rho1
    real(8) :: time_cache, time_project, time_comm, time_norm
    real(8) :: time_comm_subgroup_reduce, time_comm_pack, time_comm_exchange, time_comm_unpack
    real(8) :: time_project_setup, time_project_psi, time_project_rho
    real(8) :: time_project_grid_prep, time_project_phi_pack, time_project_overhead
    real(8) :: time_project_dmat_build
    real(8) :: t_dmat0, t_dmat1
    real(8) :: time_cache_pw_refresh, time_cache_phi_block_refresh
    logical :: use_mixed_density
    logical :: enable_density_phi_block_cache
    logical :: enable_density_stage_contrib_trace
    logical :: enable_density_state_charge_trace
    logical :: rebuilt_pw_cache, rebuilt_phi_block_cache
    logical :: need_pw_cache_alloc, need_pw_cache_expand
    logical :: need_phi_cache_alloc, need_phi_count_alloc, need_phi_cache_invalid, need_phi_cache_resize
    integer, allocatable :: ix_buf(:), iy_buf(:), iz_buf(:), owner_buf(:), ixg_buf(:), iyg_buf(:), izg_buf(:)
    integer, allocatable :: slot_buf(:), local_grid_ids(:), remote_grid_ids(:), valid_remote_grid_ids(:)
    integer, allocatable :: basis_gid(:), valid_basis_ids(:)
    integer, allocatable :: basis_gid_spin(:,:), valid_basis_ids_spin(:,:), valid_basis_count_spin(:)
    real(8), allocatable :: basis_sdiag_probe(:,:,:)
    real(8), allocatable :: phi_col_metric_total(:,:,:)
    real(8), allocatable :: basis_smat_probe(:,:,:,:)
    real(8), allocatable :: phi_gram_total(:,:,:,:)
    real(8), allocatable :: phi_frag_metric_total(:,:,:,:)
    real(8), allocatable :: basis_frag_metric_total(:,:,:,:)
    type(s_scalar), allocatable :: rho_send(:), rho_recv(:)
    integer, allocatable :: send_counts(:), recv_counts(:), send_displs(:), recv_displs(:)
    integer, allocatable :: send_total_local_by_rank(:), recv_total_local_by_rank(:)
    integer, allocatable :: send_total_all_ranks(:), recv_total_all_ranks(:)
    real(8), allocatable :: send_flat(:), recv_flat(:)
    real(8), allocatable :: rho_bf(:,:,:), rho_s_bf(:,:,:,:)
    real(8), allocatable :: rho_rank_buf(:,:,:,:), send_rank_buf(:,:,:,:)
    real(8), allocatable :: phi_blk(:,:), rho_blk(:), rho_blk_accum(:), rho_blk_reduced(:), coef_blk_re(:,:), coef_blk_im(:,:), psi_blk_re(:,:), psi_blk_im(:,:)
    real(8), allocatable :: coef_blk_ri(:,:), psi_blk_ri(:,:)
    real(8), allocatable :: density_mat_re(:,:), density_tmp(:,:)
    real(8), allocatable :: D_frag_re(:,:,:)   ! (nbf_max, nbf_max, nspin) pre-computed D per fragment
    real(8), allocatable :: coef_re_frag(:,:)   ! (nbf_max, nocc_cache) real coef for current fragment
    real(8), allocatable :: coef_im_frag(:,:)   ! (nbf_max, nocc_cache) imag coef for current fragment
    real(8), allocatable :: D_partial_re(:,:)    ! (nbf_max, nbf_max) partial D per rank
    real(8), allocatable :: coef_re_full(:,:,:)  ! (nbf_max, nocc_cache, nspin) upfront bcast coef (n_pw>0)
    real(8), allocatable :: coef_im_full(:,:,:)  ! (nbf_max, nocc_cache, nspin)
    real(8), allocatable :: rho_blk_partial(:)   ! (grid_block_size) partial rho for state slice
    real(8), allocatable :: state_charge_local(:,:), state_charge_global(:,:)
    real(8), allocatable :: state_coeff_c2_local(:,:), state_coeff_c2_global(:,:)
    real(8), allocatable :: state_psi2_raw_local(:,:), state_psi2_raw_global(:,:)
    real(8), allocatable :: state_psi2_occ_local(:,:), state_psi2_occ_global(:,:)
    real(8), allocatable :: state_psi2_dv_local(:,:), state_psi2_dv_global(:,:)
    real(8), allocatable :: state_psi2_owned_local(:,:), state_psi2_owned_global(:,:)
    real(8), allocatable :: state_import_core_local(:,:), state_import_core_global(:,:)
    real(8), allocatable :: probe_state_owned_local(:,:), probe_state_owned_global(:,:)
    real(8), allocatable :: probe_state_import_local(:,:), probe_state_import_global(:,:)
    real(8), allocatable :: occ_cache(:), occ_sqrt_cache(:)
    complex(8), allocatable :: coef_c_full(:,:), coef_c_frag(:,:)
    real(8) :: time_project_rho_reduce, time_project_phi_block_build
    real(8) :: D_partial_trace
    integer :: io_s_frag, io_e_frag, nocc_loc, nocc_per_rank_loc
    integer :: nblocks_max, block_cache_idx, npt_cache, rem_xy
    integer :: phi_lb1, phi_lb2, phi_lb3, phi_ub1, phi_ub2, phi_ub3
    integer :: phi_lg1, phi_lg2, phi_lg3
    integer :: ibuf_x, ibuf_y, ibuf_z
    integer :: rho_x_lo, rho_x_hi, rho_y_lo, rho_y_hi, rho_z_lo, rho_z_hi
    integer :: rho_s_x_lo, rho_s_x_hi, rho_s_y_lo, rho_s_y_hi, rho_s_z_lo, rho_s_z_hi
    complex(8), allocatable :: phase_cache(:,:), coef_pw_blk(:,:), pw_tmp_z(:,:)
    complex(8), allocatable :: density_mix(:,:,:), basis_mix_blk(:,:), density_mix_tmp(:,:)
    complex(8), allocatable :: basis_mix_blk_t(:,:), density_mix_tmp_t(:,:)
    complex(8), allocatable :: transform_frag_spin(:,:,:), transform_pw_spin(:,:,:)
    complex(8), allocatable :: mix_transform_spin(:,:), mix_overlap_spin(:,:), s_mix(:,:), s_mix_work(:,:)
    complex(8), allocatable :: coef_mix_eff(:,:), coef_mix_metric(:,:), coef_mix_spin(:,:,:)
    complex(8), allocatable :: d_raw_ff(:,:), d_raw_fp(:,:), d_raw_pp(:,:)
    integer, allocatable :: ipiv_mix(:)
    integer, allocatable :: n_basis_mix_spin(:)
    complex(8), allocatable :: coef_probe_full(:,:), coef_probe_pw(:,:), overlap_probe(:,:), overlap_vec(:)
    integer, allocatable :: subgroup_self_ixg_tmp(:), subgroup_self_iyg_tmp(:), subgroup_self_izg_tmp(:)
    logical, parameter :: enable_density_reconstruct_trace = .false.
    real(8) :: coef_norm_probe, rho_probe_charge, phi_norm_probe, psi_norm_probe
    real(8) :: phi_col_norm_probe(3), phi_col_hvol_probe(3), phi_sdiag_probe(3)
    real(8) :: coef_map_local_probe(3,3), coef_map_global_probe(3,3), coef_map_diff_probe(3,3)
    real(8) :: orbital_norm_probe_local(3), orbital_norm_probe_total(3)
    real(8) :: orbital_norm_frag_local(2,3), orbital_norm_frag_total(2,3)
    real(8) :: overlap_state_probe(3), overlap_elec_probe, coef_state_probe(3), overlap_diag_probe(3)
    real(8) :: overlap_self_probe(2,3), overlap_cross_probe(3)
    real(8) :: frag_trace_probe, frag_state_trace_probe(3)
    real(8) :: frag_state_real_probe(3)
    real(8) :: charge_root_tmp_global, charge_root_sum_global
    real(8) :: charge_weighted_total_global, charge_weighted_total_pre_norm
    real(8) :: charge_owner_local_pre_comm, charge_owner_global_pre_comm
    real(8) :: charge_rho_bf_all_local_raw, charge_rho_bf_all_global_raw
    real(8) :: charge_rho_bf_core_local_raw, charge_rho_bf_core_global_raw
    real(8) :: charge_imported_local, charge_imported_global
    real(8) :: charge_imported_core_local, charge_imported_core_global
    real(8) :: charge_imported_buffer_local, charge_imported_buffer_global
    real(8) :: charge_contract_residual
    real(8) :: charge_blk_all, charge_blk_owner, charge_blk_handler, charge_blk_slot0
    real(8) :: g211_blk_all_re, g211_blk_all_im, g211_blk_owner_re, g211_blk_owner_im
    real(8) :: g211_blk_handler_re, g211_blk_handler_im
    real(8) :: g211_root_sum_re_local, g211_root_sum_im_local
    real(8) :: g211_rank_buf_re_local, g211_rank_buf_im_local
    real(8) :: g211_rank_buf_re_global, g211_rank_buf_im_global
    real(8) :: g211_pre_total_re, g211_pre_total_im
    real(8) :: phase_theta, phase_re, phase_im, g211_pred_re, g211_pred_im
    real(8) :: g211_re_line, g211_im_line, rho_val
    real(8), allocatable :: g211_cos_x(:), g211_sin_x(:)
    real(8), allocatable :: kpw_hx(:), kpw_hy(:), kpw_hz(:)
    character(32) :: env_phi_block_cache
    character(32) :: env_stage_trace, env_state_trace
    character(32) :: env_point_dup_audit
    character(32) :: env_weight_path_trace
    character(32) :: env_state_rhobf_trace
    character(32) :: env_state_rhobf_io
    character(32) :: env_state_rhobf_spin
    character(32) :: env_owned_path_mode
    character(32) :: env_owned_path_scale
    character(32) :: env_imported_unpack_mode
    character(32) :: env_imported_unpack_scale
    character(32) :: env_rho_mix_mode
    character(32) :: env_rho_mix_trace
    character(32) :: env_rho_mix_raw_trace
    character(32) :: env_rho_mix_grid_compare
    character(96) :: env_point_probe
    integer :: env_status
    integer :: env_probe_status
    integer :: rho_mix_mode_kind
    integer :: info_lapack
    integer :: itt_tag
    real(8) :: state_charge_sum_spin, state_charge_sum_all
    real(8) :: state_coeff_c2_sum_spin, state_coeff_c2_sum_all
    real(8) :: state_psi2_raw_sum_spin, state_psi2_occ_sum_spin
    real(8) :: state_psi2_dv_sum_spin, state_psi2_owned_sum_spin
    real(8) :: state_import_sum_spin
    real(8) :: state_ratio_c2
    real(8) :: state_ratio_occ_raw, state_ratio_dv_occ, state_ratio_owned_dv
    real(8) :: state_total_q, state_dev_q, state_ratio_total2
    real(8) :: state_ratio_import_dv, state_ratio_total_dv
    real(8) :: psi2_val
    real(8) :: ifrag_owned_point_count_local, ifrag_owned_point_count_global
    real(8) :: ifrag_import_point_count_local, ifrag_import_point_count_global
    real(8) :: probe_owned_pre_weight_local, probe_owned_pre_weight_global
    real(8) :: probe_owned_add_local, probe_owned_add_global
    real(8) :: probe_import_send_pre_weight_local, probe_import_send_pre_weight_global
    real(8) :: probe_import_send_add_local, probe_import_send_add_global
    real(8) :: probe_import_unpack_pre_weight_local, probe_import_unpack_pre_weight_global
    real(8) :: probe_import_unpack_add_local, probe_import_unpack_add_global
    real(8) :: probe_partition_weight, probe_overlap_weight, probe_norm_weight
    real(8) :: probe_owned_apply_count_local, probe_owned_apply_count_global
    real(8) :: probe_owned_weight_sum_local, probe_owned_weight_sum_global
    real(8) :: probe_import_send_apply_count_local, probe_import_send_apply_count_global
    real(8) :: probe_import_unpack_apply_count_local, probe_import_unpack_apply_count_global
    real(8) :: probe_import_unpack_weight_sum_local, probe_import_unpack_weight_sum_global
    real(8) :: imported_unpack_norm_trigger_local, imported_unpack_norm_trigger_global
    real(8) :: state_rhobf_psi2_q_local, state_rhobf_psi2_q_global
    real(8) :: state_rhobf_psi2_occ_q_local, state_rhobf_psi2_occ_q_global
    real(8) :: state_rhobf_psi2_dv_q_local, state_rhobf_psi2_dv_q_global
    real(8) :: state_rhobf_psi2_owned_q_local, state_rhobf_psi2_owned_q_global
    real(8) :: state_rhobf_psi2_after_partition_q_local, state_rhobf_psi2_after_partition_q_global
    real(8) :: state_rhobf_psi2_after_slot_q_local, state_rhobf_psi2_after_slot_q_global
    real(8) :: state_rhobf_psi2_after_any_norm_q_local, state_rhobf_psi2_after_any_norm_q_global
    real(8) :: state_rhobf_state_total_q
    logical :: enable_owned_path_normalized
    logical :: enable_state_rhobf_trace
    logical :: enable_density_point_dup_audit
    logical :: enable_density_weight_path_trace
    logical :: enable_imported_unpack_normalized
    logical :: has_density_point_probe
    logical :: found_duplicate_point_local
    integer :: probe_ixg, probe_iyg, probe_izg
    integer :: dup_ixg_local, dup_iyg_local, dup_izg_local
    integer :: dup_src_local, dup_tgt_local, dup_slot_local
    integer :: first_imp_ixg_local, first_imp_iyg_local, first_imp_izg_local
    integer :: first_imp_src_local, first_imp_tgt_local, first_imp_slot_local
    integer :: dup_ixg_global, dup_iyg_global, dup_izg_global
    integer :: dup_src_global, dup_tgt_global, dup_slot_global
    real(8) :: dup_local_contrib, dup_import_contrib
    real(8) :: imported_unpack_weight
    real(8) :: owned_path_scale
    real(8) :: imported_unpack_scale
    real(8) :: s_mix_dev_frob, s_mix_offdiag_frob, s_mix_diag_min, s_mix_diag_max
    real(8) :: tr_ff_probe, tr_fp_probe, tr_pp_probe, fp_frob_probe, fp_maxabs_probe
    real(8) :: trs_ff_probe, trs_fp_probe, trs_pp_probe, trs_total_probe, fp_weight_frac_probe
    real(8) :: phase_applied_total_block, send_pre_block_raw
    complex(8) :: tr_fp_complex
    real(8) :: rho_mix_grid_l2, rho_mix_grid_ref, rho_mix_grid_max, rho_mix_grid_rel
    real(8) :: coef_c2_probe(3), coef_metric_probe(3)
    real(8) :: density_mix_trace_probe, density_mix_diag_min, density_mix_diag_max
    real(8) :: psi2_occ_val, psi2_dv_val
    real(8) :: psi2_after_partition_val, psi2_after_slot_val, psi2_after_any_norm_val
    integer :: state_rhobf_trace_io, state_rhobf_trace_spin
    logical :: enable_rho_mix_trace
    logical :: enable_rho_mix_raw_trace
    logical :: enable_rho_mix_grid_compare
    logical :: enable_ifrag_compare_trace
    character(32) :: env_ifrag_compare_trace
    logical :: enable_fp_decomp_audit
    character(32) :: env_fp_decomp_audit
    logical :: enable_fp_phase_fix
    character(32) :: env_fp_phase_fix
    logical :: enable_tf_occmap_only
    character(32) :: env_tf_occmap_only
    character(256) :: env_tf_occmap_itts
    character(256) :: env_tf_occmap_itts_work
    integer, parameter :: max_occmap_trace_itt = 32
    integer :: occmap_trace_itts(max_occmap_trace_itt), n_occmap_trace_itt
    integer :: occmap_trace_pos, occmap_trace_next, occmap_trace_val, occmap_trace_ios, occmap_trace_idx
    character(32) :: occmap_trace_tok
    logical :: trace_occmap_itt
    logical :: ifrag_grid_seen(2), ifrag_basis_seen(2), ifrag_fp_seen(2)
    logical :: ifrag_decomp_seen(2)
    logical :: ifrag_ff_occ_seen(2)
    integer :: cmp_nxyz(3, 2), cmp_nbf(2), cmp_valid(2), cmp_basis_gid(3, 2)
    integer :: cmp_frag_lo(3, 2), cmp_frag_hi(3, 2)
    integer :: cmp_npt(2), cmp_local(2), cmp_remote(2), cmp_valid_remote(2)
    integer :: cmp_owner_valid(2), cmp_slot0(2), cmp_slotp(2), cmp_owner_true(2), cmp_owner_false(2)
    real(8) :: cmp_phase_total(2), cmp_send_pre(2), cmp_trs_fp(2), cmp_fp_frac(2), cmp_fp_frob(2), cmp_fp_max(2)
    real(8) :: cmp_trs_ff(2), cmp_trs_pp(2), cmp_trs_tot(2)
    real(8), allocatable :: cmp_ff_occ(:,:), cmp_ff_occ_global(:,:)
    real(8), allocatable :: cmp_ff_gid_pre(:,:,:), cmp_ff_gid_post(:,:,:)
    real(8), allocatable :: cmp_ff_gid_pre_global(:,:,:), cmp_ff_gid_post_global(:,:,:)
    real(8), allocatable :: cmp_tf_gid_pre(:,:,:), cmp_tf_gid_full(:,:,:), cmp_tf_gid_int(:,:,:)
    real(8), allocatable :: cmp_tf_gid_pre_global(:,:,:), cmp_tf_gid_full_global(:,:,:), cmp_tf_gid_int_global(:,:,:)
    real(8), allocatable :: cmp_tf_gid_mode_pre(:,:,:,:), cmp_tf_gid_mode_full(:,:,:,:)
    real(8), allocatable :: cmp_tf_gid_mode_pre_global(:,:,:,:), cmp_tf_gid_mode_full_global(:,:,:,:)
    complex(8), allocatable :: cmp_tf_gid_mode_ovl(:,:,:,:), cmp_tf_gid_mode_ovl_global(:,:,:,:)
    integer, parameter :: n_focus_gid = 2
    integer, parameter :: max_occmap_ref_occ = 32768
    integer, save :: cmp_tf_m2_ref_static(max_occmap_ref_occ, n_focus_gid) = 0
    integer, parameter :: max_occmap_track_step = 64
    integer, save :: cmp_tf_track_m2(max_occmap_ref_occ, n_focus_gid, max_occmap_track_step) = 0
    integer, save :: cmp_tf_track_itt(max_occmap_track_step) = 0
    integer, save :: cmp_tf_track_nstep = 0
    integer, allocatable :: cmp_ff_dom_gid(:,:)
    real(8), allocatable :: cmp_ff_dom_gid_local(:,:), cmp_ff_dom_gid_global(:,:)
    real(8) :: cmp_recv_post_raw(2)
    integer, parameter :: n_top_fp_pairs = 5
    integer, parameter :: n_top_ff_gid = 5
    integer, parameter :: focus_gid_ids(n_focus_gid) = (/14, 26/)
    integer, parameter :: n_top_tf_mode = 5
    real(8) :: top_fp_contrib_val(n_top_fp_pairs), fp_pair_contrib
    integer :: top_fp_io(n_top_fp_pairs), top_fp_ipw(n_top_fp_pairs), top_fp_gid(n_top_fp_pairs)
    integer :: k_top, k_min, io_top, ipw_top
    integer :: dom_m, ikx_pw, iky_pw, ikz_pw
    real(8) :: dom_m_abs, dom_sign, g_abs, m_term_abs
    complex(8) :: m_term
    real(8), allocatable :: ifrag_recv_post_raw_local(:,:), ifrag_recv_post_raw_global(:,:)

    itt_tag = -1
    if (present(itt_debug)) itt_tag = itt_debug
    call capture_occmap_pair_snapshot(dg_frag, itt_tag, 'density_reconstruct_pre')

    if (itt_tag == 1) then
      write(*,'(1x,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] ENTER_DENSITY_RECON rank=', dg_frag%id, ' itt=', itt_tag, &
        ' coef_ref_ready=', dg_frag%coef_ref_ready
      flush(6)
    end if

    rho%f = 0.0d0
    do ispin = 1, system%nspin
      rho_s(ispin)%f = 0.0d0
    end do
    call cpu_time(t_total0)
    time_cache = 0.0d0
    time_project = 0.0d0
    time_comm = 0.0d0
    time_norm = 0.0d0
    time_comm_subgroup_reduce = 0.0d0
    time_comm_pack = 0.0d0
    time_comm_exchange = 0.0d0
    time_comm_unpack = 0.0d0
    time_project_setup = 0.0d0
    time_project_psi = 0.0d0
    time_project_rho = 0.0d0
    time_project_grid_prep = 0.0d0
    time_project_phi_pack = 0.0d0
    time_project_overhead = 0.0d0
    time_project_dmat_build = 0.0d0
    time_project_rho_reduce = 0.0d0
    time_project_phi_block_build = 0.0d0
    time_cache_pw_refresh = 0.0d0
    time_cache_phi_block_refresh = 0.0d0
    orbital_norm_probe_local(:) = 0.0d0
    orbital_norm_probe_total(:) = 0.0d0
    orbital_norm_frag_local(:, :) = 0.0d0
    orbital_norm_frag_total(:, :) = 0.0d0
    overlap_state_probe(:) = 0.0d0
    overlap_elec_probe = 0.0d0
    coef_state_probe(:) = 0.0d0
    overlap_diag_probe(:) = 0.0d0
    overlap_self_probe(:, :) = 0.0d0
    overlap_cross_probe(:) = 0.0d0
    frag_trace_probe = 0.0d0
    frag_state_trace_probe(:) = 0.0d0
    charge_root_tmp_global = 0.0d0
    charge_root_sum_global = 0.0d0
    charge_weighted_total_global = 0.0d0
    charge_weighted_total_pre_norm = 0.0d0
    charge_owner_local_pre_comm = 0.0d0
    charge_owner_global_pre_comm = 0.0d0
    charge_rho_bf_all_local_raw = 0.0d0
    charge_rho_bf_all_global_raw = 0.0d0
    charge_rho_bf_core_local_raw = 0.0d0
    charge_rho_bf_core_global_raw = 0.0d0
    charge_imported_local = 0.0d0
    charge_imported_global = 0.0d0
    charge_imported_core_local = 0.0d0
    charge_imported_core_global = 0.0d0
    charge_imported_buffer_local = 0.0d0
    charge_imported_buffer_global = 0.0d0
    charge_contract_residual = 0.0d0
    rebuilt_pw_cache = .false.
    rebuilt_phi_block_cache = .false.
    enable_density_phi_block_cache = .false.
    enable_density_stage_contrib_trace = .false.
    enable_density_state_charge_trace = .false.
    enable_density_point_dup_audit = .false.
    enable_density_weight_path_trace = .false.
    enable_state_rhobf_trace = .false.
    enable_owned_path_normalized = .false.
    enable_imported_unpack_normalized = .false.
    has_density_point_probe = .false.
    env_phi_block_cache = ''
    env_stage_trace = ''
    env_state_trace = ''
    env_point_dup_audit = ''
    env_weight_path_trace = ''
    env_state_rhobf_trace = ''
    env_state_rhobf_io = ''
    env_state_rhobf_spin = ''
    env_owned_path_mode = ''
    env_owned_path_scale = ''
    env_imported_unpack_mode = ''
    env_imported_unpack_scale = ''
    env_rho_mix_mode = 'legacy'
    env_rho_mix_trace = ''
    env_rho_mix_raw_trace = ''
    env_rho_mix_grid_compare = ''
    env_ifrag_compare_trace = ''
    env_fp_decomp_audit = ''
    env_fp_phase_fix = ''
    env_tf_occmap_only = ''
    env_tf_occmap_itts = ''
    env_tf_occmap_itts_work = ''
    rho_mix_mode_kind = 0
    enable_rho_mix_trace = .false.
    enable_rho_mix_raw_trace = .false.
    enable_rho_mix_grid_compare = .false.
    enable_ifrag_compare_trace = .false.
    enable_fp_decomp_audit = .false.
    enable_fp_phase_fix = .false.
    enable_tf_occmap_only = .false.
    n_occmap_trace_itt = 0
    occmap_trace_itts(:) = 0
    trace_occmap_itt = .false.
    ifrag_grid_seen(:) = .false.
    ifrag_basis_seen(:) = .false.
    ifrag_fp_seen(:) = .false.
    ifrag_decomp_seen(:) = .false.
    ifrag_ff_occ_seen(:) = .false.
    cmp_nxyz(:, :) = 0
    cmp_nbf(:) = 0
    cmp_valid(:) = 0
    cmp_basis_gid(:, :) = 0
    cmp_frag_lo(:, :) = 0
    cmp_frag_hi(:, :) = 0
    cmp_npt(:) = 0
    cmp_local(:) = 0
    cmp_remote(:) = 0
    cmp_valid_remote(:) = 0
    cmp_owner_valid(:) = 0
    cmp_slot0(:) = 0
    cmp_slotp(:) = 0
    cmp_owner_true(:) = 0
    cmp_owner_false(:) = 0
    cmp_phase_total(:) = 0.0d0
    cmp_send_pre(:) = 0.0d0
    cmp_trs_fp(:) = 0.0d0
    cmp_fp_frac(:) = 0.0d0
    cmp_fp_frob(:) = 0.0d0
    cmp_fp_max(:) = 0.0d0
    cmp_trs_ff(:) = 0.0d0
    cmp_trs_pp(:) = 0.0d0
    cmp_trs_tot(:) = 0.0d0
    cmp_recv_post_raw(:) = 0.0d0
    env_point_probe = ''
    probe_ixg = 0
    probe_iyg = 0
    probe_izg = 0
    dup_ixg_local = -1
    dup_iyg_local = -1
    dup_izg_local = -1
    dup_src_local = -1
    dup_tgt_local = -1
    dup_slot_local = -1
    dup_ixg_global = -1
    dup_iyg_global = -1
    dup_izg_global = -1
    dup_src_global = -1
    dup_tgt_global = -1
    dup_slot_global = -1
    dup_local_contrib = 0.0d0
    dup_import_contrib = 0.0d0
    first_imp_ixg_local = -1
    first_imp_iyg_local = -1
    first_imp_izg_local = -1
    first_imp_src_local = -1
    first_imp_tgt_local = -1
    first_imp_slot_local = -1
    found_duplicate_point_local = .false.
    probe_partition_weight = 1.0d0
    probe_overlap_weight = 1.0d0
    probe_norm_weight = 1.0d0
    probe_owned_pre_weight_local = 0.0d0
    probe_owned_pre_weight_global = 0.0d0
    probe_owned_add_local = 0.0d0
    probe_owned_add_global = 0.0d0
    probe_import_send_pre_weight_local = 0.0d0
    probe_import_send_pre_weight_global = 0.0d0
    probe_import_send_add_local = 0.0d0
    probe_import_send_add_global = 0.0d0
    probe_import_unpack_pre_weight_local = 0.0d0
    probe_import_unpack_pre_weight_global = 0.0d0
    probe_import_unpack_add_local = 0.0d0
    probe_import_unpack_add_global = 0.0d0
    probe_owned_apply_count_local = 0.0d0
    probe_owned_weight_sum_local = 0.0d0
    probe_owned_apply_count_global = 0.0d0
    probe_owned_weight_sum_global = 0.0d0
    probe_import_send_apply_count_local = 0.0d0
    probe_import_send_apply_count_global = 0.0d0
    probe_import_unpack_apply_count_local = 0.0d0
    probe_import_unpack_apply_count_global = 0.0d0
    probe_import_unpack_weight_sum_local = 0.0d0
    probe_import_unpack_weight_sum_global = 0.0d0
    imported_unpack_norm_trigger_local = 0.0d0
    imported_unpack_norm_trigger_global = 0.0d0
    state_rhobf_psi2_q_local = 0.0d0
    state_rhobf_psi2_q_global = 0.0d0
    state_rhobf_psi2_occ_q_local = 0.0d0
    state_rhobf_psi2_occ_q_global = 0.0d0
    state_rhobf_psi2_dv_q_local = 0.0d0
    state_rhobf_psi2_dv_q_global = 0.0d0
    state_rhobf_psi2_owned_q_local = 0.0d0
    state_rhobf_psi2_owned_q_global = 0.0d0
    state_rhobf_psi2_after_partition_q_local = 0.0d0
    state_rhobf_psi2_after_partition_q_global = 0.0d0
    state_rhobf_psi2_after_slot_q_local = 0.0d0
    state_rhobf_psi2_after_slot_q_global = 0.0d0
    state_rhobf_psi2_after_any_norm_q_local = 0.0d0
    state_rhobf_psi2_after_any_norm_q_global = 0.0d0
    state_rhobf_state_total_q = 0.0d0
    state_rhobf_trace_io = 1
    state_rhobf_trace_spin = 1
    owned_path_scale = 0.5d0
    imported_unpack_scale = 0.5d0
    call get_environment_variable('SALMON_DG_DENSITY_PHI_BLOCK_CACHE', env_phi_block_cache, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_phi_block_cache)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enable_density_phi_block_cache = .true.
      end select
    end if
    call get_environment_variable('SALMON_DG_DENSITY_STAGE_CONTRIB_TRACE', env_stage_trace, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_stage_trace)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enable_density_stage_contrib_trace = .true.
      end select
    end if
    call get_environment_variable('SALMON_DG_DENSITY_STATE_CHARGE_TRACE', env_state_trace, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_state_trace)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enable_density_state_charge_trace = .true.
      end select
    end if
    call get_environment_variable('SALMON_DG_DENSITY_POINT_DUP_AUDIT', env_point_dup_audit, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_point_dup_audit)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enable_density_point_dup_audit = .true.
      end select
    end if
    call get_environment_variable('SALMON_DG_DENSITY_WEIGHT_PATH_TRACE', env_weight_path_trace, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_weight_path_trace)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enable_density_weight_path_trace = .true.
      end select
    end if
    call get_environment_variable('SALMON_DG_STATE_RHOBF_TRACE', env_state_rhobf_trace, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_state_rhobf_trace)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enable_state_rhobf_trace = .true.
      end select
    end if
    call get_environment_variable('SALMON_DG_STATE_RHOBF_IO', env_state_rhobf_io, status=env_status)
    if (env_status == 0) then
      read(env_state_rhobf_io, *, iostat=env_status) state_rhobf_trace_io
      if (env_status /= 0) state_rhobf_trace_io = 1
    end if
    call get_environment_variable('SALMON_DG_STATE_RHOBF_SPIN', env_state_rhobf_spin, status=env_status)
    if (env_status == 0) then
      read(env_state_rhobf_spin, *, iostat=env_status) state_rhobf_trace_spin
      if (env_status /= 0) state_rhobf_trace_spin = 1
    end if
    call get_environment_variable('SALMON_DG_OWNED_PATH_MODE', env_owned_path_mode, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_owned_path_mode)))
      case('normalized','NORMALIZED','norm','NORM')
        enable_owned_path_normalized = .true.
      case default
        enable_owned_path_normalized = .false.
      end select
    end if
    call get_environment_variable('SALMON_DG_OWNED_PATH_SCALE', env_owned_path_scale, status=env_status)
    if (env_status == 0) then
      read(env_owned_path_scale, *, iostat=env_status) owned_path_scale
      if (env_status /= 0) owned_path_scale = 0.5d0
    end if
    call get_environment_variable('SALMON_DG_IMPORTED_UNPACK_MODE', env_imported_unpack_mode, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_imported_unpack_mode)))
      case('normalized','NORMALIZED','norm','NORM')
        enable_imported_unpack_normalized = .true.
      case default
        enable_imported_unpack_normalized = .false.
      end select
    end if
    call get_environment_variable('SALMON_DG_IMPORTED_UNPACK_SCALE', env_imported_unpack_scale, status=env_status)
    if (env_status == 0) then
      read(env_imported_unpack_scale, *, iostat=env_status) imported_unpack_scale
      if (env_status /= 0) imported_unpack_scale = 0.5d0
    end if
    call get_environment_variable('SALMON_DG_RHO_MIX_MODE', env_rho_mix_mode, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_rho_mix_mode)))
      case('orthonormal_cc','ORTHONORMAL_CC')
        rho_mix_mode_kind = 1
      case('metric_consistent','METRIC_CONSISTENT','overlap_metric','OVERLAP_METRIC')
        rho_mix_mode_kind = 2
      case default
        rho_mix_mode_kind = 0
        env_rho_mix_mode = 'legacy'
      end select
    else
      env_rho_mix_mode = 'legacy'
    end if
    call get_environment_variable('SALMON_DG_RHO_MIX_TRACE', env_rho_mix_trace, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_rho_mix_trace)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enable_rho_mix_trace = .true.
      end select
    end if
    enable_rho_mix_raw_trace = enable_rho_mix_trace
    enable_rho_mix_grid_compare = enable_rho_mix_trace
    call get_environment_variable('SALMON_DG_RHO_MIX_RAW_TRACE', env_rho_mix_raw_trace, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_rho_mix_raw_trace)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enable_rho_mix_raw_trace = .true.
      case('0','n','N','no','NO','false','FALSE','off','OFF')
        enable_rho_mix_raw_trace = .false.
      end select
    end if
    call get_environment_variable('SALMON_DG_RHO_MIX_GRID_COMPARE', env_rho_mix_grid_compare, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_rho_mix_grid_compare)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enable_rho_mix_grid_compare = .true.
      case('0','n','N','no','NO','false','FALSE','off','OFF')
        enable_rho_mix_grid_compare = .false.
      end select
    end if
    call get_environment_variable('SALMON_DG_IFRAG_COMPARE_TRACE', env_ifrag_compare_trace, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_ifrag_compare_trace)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enable_ifrag_compare_trace = .true.
      case('0','n','N','no','NO','false','FALSE','off','OFF')
        enable_ifrag_compare_trace = .false.
      end select
    end if
    call get_environment_variable('SALMON_DG_FP_DECOMP_AUDIT', env_fp_decomp_audit, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_fp_decomp_audit)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enable_fp_decomp_audit = .true.
      case('0','n','N','no','NO','false','FALSE','off','OFF')
        enable_fp_decomp_audit = .false.
      end select
    end if
    if (enable_fp_decomp_audit) then
      enable_rho_mix_raw_trace = .true.
    end if
    call get_environment_variable('SALMON_DG_FP_PHASE_FIX', env_fp_phase_fix, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_fp_phase_fix)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enable_fp_phase_fix = .true.
      case('0','n','N','no','NO','false','FALSE','off','OFF')
        enable_fp_phase_fix = .false.
      end select
    end if
    call get_environment_variable('SALMON_DG_TF_OCCMAP_ONLY', env_tf_occmap_only, status=env_status)
    if (env_status == 0) then
      select case (adjustl(trim(env_tf_occmap_only)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enable_tf_occmap_only = .true.
      case('0','n','N','no','NO','false','FALSE','off','OFF')
        enable_tf_occmap_only = .false.
      end select
    end if
    call get_environment_variable('SALMON_DG_TF_OCCMAP_ITTS', env_tf_occmap_itts, status=env_status)
    if (env_status == 0) then
      env_tf_occmap_itts_work = trim(env_tf_occmap_itts)
      do occmap_trace_idx = 1, len_trim(env_tf_occmap_itts_work)
        if (env_tf_occmap_itts_work(occmap_trace_idx:occmap_trace_idx) == ',') env_tf_occmap_itts_work(occmap_trace_idx:occmap_trace_idx) = ' '
      end do
      occmap_trace_pos = 1
      do while (occmap_trace_pos <= len_trim(env_tf_occmap_itts_work))
        do while (occmap_trace_pos <= len_trim(env_tf_occmap_itts_work) .and. env_tf_occmap_itts_work(occmap_trace_pos:occmap_trace_pos) == ' ')
          occmap_trace_pos = occmap_trace_pos + 1
        end do
        if (occmap_trace_pos > len_trim(env_tf_occmap_itts_work)) exit
        occmap_trace_next = occmap_trace_pos
        do while (occmap_trace_next <= len_trim(env_tf_occmap_itts_work) .and. env_tf_occmap_itts_work(occmap_trace_next:occmap_trace_next) /= ' ')
          occmap_trace_next = occmap_trace_next + 1
        end do
        occmap_trace_tok = ''
        occmap_trace_tok = env_tf_occmap_itts_work(occmap_trace_pos:occmap_trace_next-1)
        read(occmap_trace_tok, *, iostat=occmap_trace_ios) occmap_trace_val
        if (occmap_trace_ios == 0) then
          if (n_occmap_trace_itt < max_occmap_trace_itt) then
            n_occmap_trace_itt = n_occmap_trace_itt + 1
            occmap_trace_itts(n_occmap_trace_itt) = occmap_trace_val
          end if
        end if
        occmap_trace_pos = occmap_trace_next + 1
      end do
    end if
    if (n_occmap_trace_itt > 0) then
      trace_occmap_itt = .false.
      do occmap_trace_idx = 1, n_occmap_trace_itt
        if (itt_tag == occmap_trace_itts(occmap_trace_idx)) then
          trace_occmap_itt = .true.
          exit
        end if
      end do
    else
      trace_occmap_itt = (itt_tag == 1 .or. itt_tag == 4 .or. itt_tag == 8)
    end if
    call get_environment_variable('SALMON_DG_DENSITY_POINT_PROBE', env_point_probe, status=env_probe_status)
    if (env_probe_status == 0) then
      read(env_point_probe, *, iostat=env_status) probe_ixg, probe_iyg, probe_izg
      if (env_status == 0) has_density_point_probe = .true.
    end if
    need_pw_cache_alloc = .false.
    need_pw_cache_expand = .false.
    need_phi_cache_alloc = .false.
    need_phi_count_alloc = .false.
    need_phi_cache_invalid = .false.
    need_phi_cache_resize = .false.
    if (enable_density_reconstruct_trace .and. dg_frag%is_frag_root) then
      write(*,'(1x,a,i0,a,i0,a,i0)') "        density trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " stage=entry_root root_world=", dg_frag%id - dg_frag%id_frag
      flush(6)
    end if

    if (.not. allocated(dg_frag%phi_frag)) return
    inv_lgnum1 = 1.0d0 / dble(max(1, dg_frag%lgnum_total(1)))

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    g211_root_sum_re_local = 0.0d0
    g211_root_sum_im_local = 0.0d0
    g211_rank_buf_re_local = 0.0d0
    g211_rank_buf_im_local = 0.0d0
    g211_rank_buf_re_global = 0.0d0
    g211_rank_buf_im_global = 0.0d0
    ngrid_max = 0
    if (ifrag_count > 0) then
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        ngrid_max = max(ngrid_max, product(dg_frag%nxyz_domain(:, ifrag)))
      end do
    end if
    allocate(ifrag_recv_post_raw_local(system%nspin, max(1, ifrag_count)))
    allocate(ifrag_recv_post_raw_global(system%nspin, max(1, ifrag_count)))
    ifrag_recv_post_raw_local(:, :) = 0.0d0
    ifrag_recv_post_raw_global(:, :) = 0.0d0
    allocate(cmp_ff_occ(max(1, nocc_cache), 2))
    allocate(cmp_ff_occ_global(max(1, nocc_cache), 2))
    allocate(cmp_ff_dom_gid(max(1, nocc_cache), 2))
    allocate(cmp_ff_dom_gid_local(max(1, nocc_cache), 2))
    allocate(cmp_ff_dom_gid_global(max(1, nocc_cache), 2))
    cmp_ff_occ(:, :) = 0.0d0
    cmp_ff_occ_global(:, :) = 0.0d0
    cmp_ff_dom_gid(:, :) = 0
    cmp_ff_dom_gid_local(:, :) = 0.0d0
    cmp_ff_dom_gid_global(:, :) = 0.0d0
    nbf_max = max(1, maxval(dg_frag%n_basis(:, 1:system%nspin)))

    allocate(ix_buf(grid_block_size), iy_buf(grid_block_size), iz_buf(grid_block_size))
    allocate(owner_buf(grid_block_size), ixg_buf(grid_block_size), iyg_buf(grid_block_size), izg_buf(grid_block_size))
    allocate(slot_buf(grid_block_size), local_grid_ids(grid_block_size), remote_grid_ids(grid_block_size))
    allocate(valid_remote_grid_ids(grid_block_size))
    allocate(basis_gid(dg_frag%nstate_frag), valid_basis_ids(dg_frag%nstate_frag))
    allocate(basis_gid_spin(nbf_max, system%nspin), valid_basis_ids_spin(nbf_max, system%nspin), valid_basis_count_spin(system%nspin))
    if (enable_density_reconstruct_trace) then
      allocate(basis_sdiag_probe(3, system%nspin, max(1, ifrag_count)))
      allocate(phi_col_metric_total(3, system%nspin, max(1, ifrag_count)))
      allocate(basis_smat_probe(3, 3, system%nspin, max(1, ifrag_count)))
      allocate(phi_gram_total(3, 3, system%nspin, max(1, ifrag_count)))
      allocate(phi_frag_metric_total(nbf_max, nbf_max, system%nspin, max(1, ifrag_count)))
      allocate(basis_frag_metric_total(nbf_max, nbf_max, system%nspin, max(1, ifrag_count)))
      basis_sdiag_probe(:, :, :) = 0.0d0
      phi_col_metric_total(:, :, :) = 0.0d0
      basis_smat_probe(:, :, :, :) = 0.0d0
      phi_gram_total(:, :, :, :) = 0.0d0
      phi_frag_metric_total(:, :, :, :) = 0.0d0
      basis_frag_metric_total(:, :, :, :) = 0.0d0
    end if
    allocate(phi_blk(grid_block_size, dg_frag%nstate_frag))
    allocate(rho_blk(grid_block_size))
    allocate(rho_blk_accum(grid_block_size))
    allocate(rho_blk_reduced(grid_block_size))
    nocc_cache = 0
    if (allocated(dg_frag%nocc_spin)) then
      nocc_cache = maxval(dg_frag%nocc_spin(1:system%nspin))
    end if
    nocc_cache = max(1, nocc_cache)
    allocate(coef_blk_re(nbf_max, state_block_size), coef_blk_im(nbf_max, state_block_size))
    allocate(psi_blk_re(grid_block_size, state_block_size), psi_blk_im(grid_block_size, state_block_size))
    allocate(state_charge_local(max(1, nocc_cache), system%nspin))
    allocate(state_charge_global(max(1, nocc_cache), system%nspin))
    allocate(state_coeff_c2_local(max(1, nocc_cache), system%nspin))
    allocate(state_coeff_c2_global(max(1, nocc_cache), system%nspin))
    allocate(state_psi2_raw_local(max(1, nocc_cache), system%nspin))
    allocate(state_psi2_raw_global(max(1, nocc_cache), system%nspin))
    allocate(state_psi2_occ_local(max(1, nocc_cache), system%nspin))
    allocate(state_psi2_occ_global(max(1, nocc_cache), system%nspin))
    allocate(state_psi2_dv_local(max(1, nocc_cache), system%nspin))
    allocate(state_psi2_dv_global(max(1, nocc_cache), system%nspin))
    allocate(state_psi2_owned_local(max(1, nocc_cache), system%nspin))
    allocate(state_psi2_owned_global(max(1, nocc_cache), system%nspin))
    allocate(state_import_core_local(max(1, nocc_cache), system%nspin))
    allocate(state_import_core_global(max(1, nocc_cache), system%nspin))
    allocate(probe_state_owned_local(max(1, nocc_cache), system%nspin))
    allocate(probe_state_owned_global(max(1, nocc_cache), system%nspin))
    allocate(probe_state_import_local(max(1, nocc_cache), system%nspin))
    allocate(probe_state_import_global(max(1, nocc_cache), system%nspin))
    state_charge_local(:, :) = 0.0d0
    state_charge_global(:, :) = 0.0d0
    state_coeff_c2_local(:, :) = 0.0d0
    state_coeff_c2_global(:, :) = 0.0d0
    state_psi2_raw_local(:, :) = 0.0d0
    state_psi2_raw_global(:, :) = 0.0d0
    state_psi2_occ_local(:, :) = 0.0d0
    state_psi2_occ_global(:, :) = 0.0d0
    state_psi2_dv_local(:, :) = 0.0d0
    state_psi2_dv_global(:, :) = 0.0d0
    state_psi2_owned_local(:, :) = 0.0d0
    state_psi2_owned_global(:, :) = 0.0d0
    state_import_core_local(:, :) = 0.0d0
    state_import_core_global(:, :) = 0.0d0
    probe_state_owned_local(:, :) = 0.0d0
    probe_state_owned_global(:, :) = 0.0d0
    probe_state_import_local(:, :) = 0.0d0
    probe_state_import_global(:, :) = 0.0d0
    allocate(coef_blk_ri(nbf_max, 2*state_block_size), psi_blk_ri(grid_block_size, 2*state_block_size))
    allocate(pw_tmp_z(grid_block_size, state_block_size))
    allocate(density_mat_re(nbf_max, nbf_max), density_tmp(grid_block_size, nbf_max))
    allocate(coef_pw_blk(pw_block_size, state_block_size))
    ibuf_x = dg_frag%nxyz_buffer(1)
    ibuf_y = dg_frag%nxyz_buffer(2)
    ibuf_z = dg_frag%nxyz_buffer(3)
    rho_x_lo = mg%is(1)
    rho_x_hi = mg%ie(1)
    rho_y_lo = mg%is(2)
    rho_y_hi = mg%ie(2)
    rho_z_lo = mg%is(3)
    rho_z_hi = mg%ie(3)
    rho_s_x_lo = rho_x_lo
    rho_s_x_hi = rho_x_hi
    rho_s_y_lo = rho_y_lo
    rho_s_y_hi = rho_y_hi
    rho_s_z_lo = rho_z_lo
    rho_s_z_hi = rho_z_hi
    allocate(rho_bf(rho_x_lo-ibuf_x:rho_x_hi+ibuf_x, &
                    rho_y_lo-ibuf_y:rho_y_hi+ibuf_y, &
                    rho_z_lo-ibuf_z:rho_z_hi+ibuf_z))
    allocate(rho_s_bf(rho_s_x_lo:rho_s_x_hi, &
                      rho_s_y_lo:rho_s_y_hi, &
                      rho_s_z_lo:rho_s_z_hi, &
      system%nspin))
    allocate(rho_send(0:dg_frag%isize-1))
    allocate(rho_recv(0:dg_frag%isize-1))

    rho%f = 0.0d0
    rho_bf(:, :, :) = 0.0d0
    rho_s_bf(:, :, :, :) = 0.0d0
    rho_blk_reduced(:) = 0.0d0
    n_pw = max(0, dg_frag%n_plane_waves)
    n_frag = dg_frag%n_mat_max
    allocate(cmp_ff_gid_pre(max(1, nocc_cache), max(1, n_frag), 2))
    allocate(cmp_ff_gid_post(max(1, nocc_cache), max(1, n_frag), 2))
    allocate(cmp_ff_gid_pre_global(max(1, nocc_cache), max(1, n_frag), 2))
    allocate(cmp_ff_gid_post_global(max(1, nocc_cache), max(1, n_frag), 2))
    allocate(cmp_tf_gid_pre(max(1, nocc_cache), n_focus_gid, 2))
    allocate(cmp_tf_gid_full(max(1, nocc_cache), n_focus_gid, 2))
    allocate(cmp_tf_gid_int(max(1, nocc_cache), n_focus_gid, 2))
    allocate(cmp_tf_gid_pre_global(max(1, nocc_cache), n_focus_gid, 2))
    allocate(cmp_tf_gid_full_global(max(1, nocc_cache), n_focus_gid, 2))
    allocate(cmp_tf_gid_int_global(max(1, nocc_cache), n_focus_gid, 2))
    allocate(cmp_tf_gid_mode_pre(max(1, nocc_cache), max(1, max_mixed_basis), n_focus_gid, 2))
    allocate(cmp_tf_gid_mode_full(max(1, nocc_cache), max(1, max_mixed_basis), n_focus_gid, 2))
    allocate(cmp_tf_gid_mode_pre_global(max(1, nocc_cache), max(1, max_mixed_basis), n_focus_gid, 2))
    allocate(cmp_tf_gid_mode_full_global(max(1, nocc_cache), max(1, max_mixed_basis), n_focus_gid, 2))
    allocate(cmp_tf_gid_mode_ovl(max(1, nocc_cache), max(1, max_mixed_basis), n_focus_gid, 2))
    allocate(cmp_tf_gid_mode_ovl_global(max(1, nocc_cache), max(1, max_mixed_basis), n_focus_gid, 2))
    if (itt_tag <= 1) then
      cmp_tf_m2_ref_static(:, :) = 0
      cmp_tf_track_m2(:, :, :) = 0
      cmp_tf_track_itt(:) = 0
      cmp_tf_track_nstep = 0
    end if
    cmp_ff_gid_pre(:, :, :) = 0.0d0
    cmp_ff_gid_post(:, :, :) = 0.0d0
    cmp_ff_gid_pre_global(:, :, :) = 0.0d0
    cmp_ff_gid_post_global(:, :, :) = 0.0d0
    cmp_tf_gid_pre(:, :, :) = 0.0d0
    cmp_tf_gid_full(:, :, :) = 0.0d0
    cmp_tf_gid_int(:, :, :) = 0.0d0
    cmp_tf_gid_pre_global(:, :, :) = 0.0d0
    cmp_tf_gid_full_global(:, :, :) = 0.0d0
    cmp_tf_gid_int_global(:, :, :) = 0.0d0
    cmp_tf_gid_mode_pre(:, :, :, :) = 0.0d0
    cmp_tf_gid_mode_full(:, :, :, :) = 0.0d0
    cmp_tf_gid_mode_pre_global(:, :, :, :) = 0.0d0
    cmp_tf_gid_mode_full_global(:, :, :, :) = 0.0d0
    cmp_tf_gid_mode_ovl(:, :, :, :) = (0.0d0, 0.0d0)
    cmp_tf_gid_mode_ovl_global(:, :, :, :) = (0.0d0, 0.0d0)
    n_tot = n_frag + n_pw
    if (n_pw > 0) then
      allocate(phase_cache(grid_block_size, n_pw))
      allocate(kpw_hx(n_pw), kpw_hy(n_pw), kpw_hz(n_pw))
      kpw_hx(1:n_pw) = dg_frag%k_pw(1, 1:n_pw) * dg_frag%hgs(1)
      kpw_hy(1:n_pw) = dg_frag%k_pw(2, 1:n_pw) * dg_frag%hgs(2)
      kpw_hz(1:n_pw) = dg_frag%k_pw(3, 1:n_pw) * dg_frag%hgs(3)
    end if

    ! Mixed-basis density reconstruction is only valid when PW channels exist.
    ! For n_pw==0, stay on the pure fragment-basis path.
    use_mixed_density = (n_pw > 0 .and. dg_frag%mixed_basis_ready .and. allocated(dg_frag%mixed_transform) .and. &
      allocated(dg_frag%coef_mix) .and. allocated(dg_frag%mixed_basis_dim))
    if ((enable_density_reconstruct_trace .or. enable_rho_mix_trace) .and. dg_frag%id == 0) then
      write(*,'(1x,a,a,a,i0)') '        density rho_mix mode: mode=', trim(adjustl(env_rho_mix_mode)), &
        ' kind=', rho_mix_mode_kind
      flush(6)
    end if
    ! Density reconstruction uses subgroup-distributed projection and collective reductions on icomm_frag.
    subgroup_root_rank = dg_frag%id - dg_frag%id_frag
    total_send_pts = 0
    do irank = 0, dg_frag%isize - 1
      npts = dg_frag%density_send_count(irank)
      if (npts <= 0) cycle
      if (.not. target_rank_owned_by_handler(irank)) cycle
      allocate(rho_send(irank)%f((system%nspin + 1) * npts, 1, 1))
      rho_send(irank)%f(:, :, :) = 0.0d0
    end do
    do irank = 0, dg_frag%isize - 1
      npts = dg_frag%density_recv_map(irank)%npts
      if (npts <= 0) cycle
      if (.not. target_rank_owned_by_handler(irank)) cycle
      allocate(rho_recv(irank)%f((system%nspin + 1) * npts, 1, 1))
      rho_recv(irank)%f(:, :, :) = 0.0d0
    end do
    if (.not. allocated(dg_frag%density_matrix_frag)) then
      if (allocated(dg_frag%density_matrix_frag)) deallocate(dg_frag%density_matrix_frag)
      allocate(dg_frag%density_matrix_frag(nbf_max, nbf_max, system%nspin, max(1, ifrag_count)))
    else if (size(dg_frag%density_matrix_frag, 1) /= nbf_max .or. size(dg_frag%density_matrix_frag, 2) /= nbf_max .or. &
             size(dg_frag%density_matrix_frag, 3) /= system%nspin .or. size(dg_frag%density_matrix_frag, 4) /= ifrag_count) then
      deallocate(dg_frag%density_matrix_frag)
      allocate(dg_frag%density_matrix_frag(nbf_max, nbf_max, system%nspin, max(1, ifrag_count)))
    end if
    if (.not. allocated(dg_frag%density_matrix_frag_valid)) then
      if (allocated(dg_frag%density_matrix_frag_valid)) deallocate(dg_frag%density_matrix_frag_valid)
      allocate(dg_frag%density_matrix_frag_valid(system%nspin, max(1, ifrag_count)))
    else if (size(dg_frag%density_matrix_frag_valid, 1) /= system%nspin .or. &
             size(dg_frag%density_matrix_frag_valid, 2) /= ifrag_count) then
      deallocate(dg_frag%density_matrix_frag_valid)
      allocate(dg_frag%density_matrix_frag_valid(system%nspin, max(1, ifrag_count)))
    end if
    dg_frag%density_matrix_frag(:, :, :, :) = (0.0d0, 0.0d0)
    dg_frag%density_matrix_frag_valid(:, :) = .false.
    max_mixed_basis = 0
    if (use_mixed_density) then
      max_mixed_basis = maxval(dg_frag%mixed_basis_dim(1:system%nspin))
      use_mixed_density = (max_mixed_basis > 0)
    end if
    if (use_mixed_density) then
      allocate(density_mix(max_mixed_basis, max_mixed_basis, system%nspin))
      density_mix(:, :, :) = (0.0d0, 0.0d0)
      allocate(coef_mix_spin(max_mixed_basis, max(1, nocc_cache), system%nspin))
      coef_mix_spin(:, :, :) = (0.0d0, 0.0d0)
      allocate(basis_mix_blk(grid_block_size, max_mixed_basis))
      allocate(density_mix_tmp(grid_block_size, max_mixed_basis))
      allocate(basis_mix_blk_t(max_mixed_basis, grid_block_size))
      allocate(density_mix_tmp_t(max_mixed_basis, grid_block_size))
      allocate(transform_frag_spin(dg_frag%nstate_frag, max_mixed_basis, system%nspin))
      allocate(n_basis_mix_spin(system%nspin))
      n_basis_mix_spin(:) = 0
      transform_frag_spin(:, :, :) = (0.0d0, 0.0d0)
      if (n_pw > 0) then
        allocate(transform_pw_spin(n_pw, max_mixed_basis, system%nspin))
        transform_pw_spin(:, :, :) = (0.0d0, 0.0d0)
      end if
      do ispin = 1, system%nspin
        nocc_spin = dg_frag%nocc_spin(ispin)
        n_basis_mix = min(dg_frag%mixed_basis_dim(ispin), max_mixed_basis, size(dg_frag%coef_mix, 1))
        if (n_basis_mix <= 0 .or. nocc_spin <= 0) cycle
        if (allocated(s_mix)) deallocate(s_mix)
        if (allocated(mix_transform_spin)) deallocate(mix_transform_spin)
        if (allocated(mix_overlap_spin)) deallocate(mix_overlap_spin)
        if (allocated(s_mix_work)) deallocate(s_mix_work)
        if (allocated(coef_mix_eff)) deallocate(coef_mix_eff)
        if (allocated(coef_mix_metric)) deallocate(coef_mix_metric)
        if (allocated(ipiv_mix)) deallocate(ipiv_mix)

        if (rho_mix_mode_kind == 2 .or. enable_rho_mix_trace) then
          allocate(mix_transform_spin(n_tot, n_basis_mix), mix_overlap_spin(n_tot, n_basis_mix), s_mix(n_basis_mix, n_basis_mix))
          mix_transform_spin(:, :) = dg_frag%mixed_transform(1:n_tot, 1:n_basis_mix, ispin)
          call apply_overlap_operator_batch(dg_frag, ispin, mix_transform_spin, mix_overlap_spin, .false.)
          s_mix(:, :) = matmul(conjg(transpose(mix_transform_spin)), mix_overlap_spin)
          if (enable_rho_mix_trace .and. dg_frag%id == 0) then
            s_mix_dev_frob = 0.0d0
            s_mix_offdiag_frob = 0.0d0
            s_mix_diag_min = huge(1.0d0)
            s_mix_diag_max = -huge(1.0d0)
            do io = 1, n_basis_mix
              s_mix_diag_min = min(s_mix_diag_min, real(s_mix(io, io), kind=8))
              s_mix_diag_max = max(s_mix_diag_max, real(s_mix(io, io), kind=8))
              s_mix_dev_frob = s_mix_dev_frob + abs(s_mix(io, io) - (1.0d0, 0.0d0))**2
              do i_local = 1, n_basis_mix
                if (i_local == io) cycle
                s_mix_offdiag_frob = s_mix_offdiag_frob + abs(s_mix(i_local, io))**2
              end do
            end do
            s_mix_dev_frob = sqrt(max(0.0d0, s_mix_dev_frob + s_mix_offdiag_frob))
            s_mix_offdiag_frob = sqrt(max(0.0d0, s_mix_offdiag_frob))
            write(*,'(1x,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') '        rho_mix metric: ispin=', ispin, &
              ' dev_frob=', s_mix_dev_frob, ' offdiag_frob=', s_mix_offdiag_frob, ' diag_min=', s_mix_diag_min, ' diag_max=', s_mix_diag_max
            flush(6)
          end if
        end if

        if (rho_mix_mode_kind == 2 .and. allocated(s_mix)) then
          allocate(coef_mix_metric(n_basis_mix, nocc_spin), coef_mix_eff(n_basis_mix, nocc_spin))
          coef_mix_metric(:, :) = dg_frag%coef_mix(1:n_basis_mix, 1:nocc_spin, ispin)
          allocate(s_mix_work(n_basis_mix, n_basis_mix), ipiv_mix(n_basis_mix))
          s_mix_work(:, :) = s_mix(:, :)
          call zgesv(n_basis_mix, nocc_spin, s_mix_work, n_basis_mix, ipiv_mix, coef_mix_metric, n_basis_mix, info_lapack)
          if (info_lapack /= 0) then
            coef_mix_eff(:, :) = dg_frag%coef_mix(1:n_basis_mix, 1:nocc_spin, ispin)
            if (dg_frag%id == 0) then
              write(*,'(1x,a,i0,a,i0)') ' [WARN] rho_mix metric_consistent fallback to legacy for ispin=', ispin, ' zgesv info=', info_lapack
              flush(6)
            end if
          else
            coef_mix_eff(:, :) = coef_mix_metric(:, :)
          end if
        else
          allocate(coef_mix_eff(n_basis_mix, nocc_spin))
          coef_mix_eff(:, :) = dg_frag%coef_mix(1:n_basis_mix, 1:nocc_spin, ispin)
        end if

        coef_mix_spin(1:n_basis_mix, 1:nocc_spin, ispin) = coef_mix_eff(1:n_basis_mix, 1:nocc_spin)

        density_mix(1:n_basis_mix, 1:n_basis_mix, ispin) = (0.0d0, 0.0d0)
        do io = 1, nocc_spin
          occ_factor = 1.0d0
          if (allocated(system%rocc)) then
            if (io <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
              occ_factor = max(0.0d0, system%rocc(io, 1, ispin))
            end if
          end if
          if (occ_factor <= 0.0d0) cycle
          density_mix(1:n_basis_mix, 1:n_basis_mix, ispin) = density_mix(1:n_basis_mix, 1:n_basis_mix, ispin) + &
            occ_factor * matmul(coef_mix_eff(1:n_basis_mix, io:io), conjg(transpose(coef_mix_eff(1:n_basis_mix, io:io))))
        end do

        if (enable_rho_mix_trace .and. dg_frag%id == 0) then
          coef_c2_probe(:) = 0.0d0
          coef_metric_probe(:) = 0.0d0
          do io = 1, min(3, nocc_spin)
            coef_c2_probe(io) = real(sum(conjg(coef_mix_eff(:, io)) * coef_mix_eff(:, io)), kind=8)
            if (allocated(s_mix)) then
              coef_metric_probe(io) = real(sum(conjg(coef_mix_eff(:, io)) * matmul(s_mix, coef_mix_eff(:, io))), kind=8)
            else
              coef_metric_probe(io) = coef_c2_probe(io)
            end if
          end do
          density_mix_trace_probe = 0.0d0
          density_mix_diag_min = huge(1.0d0)
          density_mix_diag_max = -huge(1.0d0)
          do io = 1, n_basis_mix
            density_mix_trace_probe = density_mix_trace_probe + real(density_mix(io, io, ispin), kind=8)
            density_mix_diag_min = min(density_mix_diag_min, real(density_mix(io, io, ispin), kind=8))
            density_mix_diag_max = max(density_mix_diag_max, real(density_mix(io, io, ispin), kind=8))
          end do
          write(*,'(1x,a,i0,a,3(1pe12.4,1x),a,3(1pe12.4,1x),a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
            '        rho_mix coef/diag: ispin=', ispin, ' c2=', coef_c2_probe(1), coef_c2_probe(2), coef_c2_probe(3), &
            ' cSc=', coef_metric_probe(1), coef_metric_probe(2), coef_metric_probe(3), ' tr=', density_mix_trace_probe, &
            ' diag_min=', density_mix_diag_min, ' diag_max=', density_mix_diag_max
          flush(6)
        end if
      end do
    end if
    boxL(1) = dg_frag%hgs(1) * real(mg%num(1), 8)
    boxL(2) = dg_frag%hgs(2) * real(mg%num(2), 8)
    boxL(3) = dg_frag%hgs(3) * real(mg%num(3), 8)
    inv_sqrt_vol = 1.0d0 / sqrt(max(1.0d-16, boxL(1) * boxL(2) * boxL(3)))

    call cpu_time(t_setup0)
    if (.not. allocated(dg_frag%density_block_nblocks)) then
      allocate(dg_frag%density_block_nblocks(ifrag_count), dg_frag%density_block_first_offset(ifrag_count), &
               dg_frag%density_block_step(ifrag_count))
      block_idx_global = 0
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        ifrag_owned_point_count_local = 0.0d0
        ifrag_import_point_count_local = 0.0d0
        nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
        ngrid = nxyz(1) * nxyz(2) * nxyz(3)
        nblocks_ifrag = (ngrid + grid_block_size - 1) / grid_block_size
        dg_frag%density_block_nblocks(i_local) = nblocks_ifrag
        dg_frag%density_block_first_offset(i_local) = 0
        dg_frag%density_block_step(i_local) = 1
        block_idx_global = block_idx_global + nblocks_ifrag
      end do
    end if
    call cpu_time(t_setup1)
    time_project_overhead = time_project_overhead + (t_setup1 - t_setup0)

    if (n_pw > 0) then
      call cpu_time(t_cache0)
      need_pw_cache_alloc = (.not. allocated(dg_frag%coef_pw_full_cache))
      need_pw_cache_expand = (.not. need_pw_cache_alloc) .and. dg_frag%coef_pw_full_cache_nstate < nocc_cache
      if (need_pw_cache_alloc .or. need_pw_cache_expand) then
        call refresh_pw_coef_cache(dg_frag, nocc_cache)
        rebuilt_pw_cache = .true.
      end if
      call cpu_time(t_cache1)
      time_cache = time_cache + (t_cache1 - t_cache0)
      if (rebuilt_pw_cache) time_cache_pw_refresh = time_cache_pw_refresh + (t_cache1 - t_cache0)
    end if
    if (enable_density_reconstruct_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a,a,1pe12.4)') "        density trace: stage=after-pw-cache dt=", "", time_cache
      write(*,'(1x,a,l1,a,l1,a,1pe12.4)') "        density cache: pw_refresh rebuilt=", rebuilt_pw_cache, &
        " expand_only=", need_pw_cache_expand, " dt=", time_cache_pw_refresh
      flush(6)
    end if

    if (enable_density_reconstruct_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density trace: stage=before-frag-loop"
      flush(6)
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        density trace: comm-state icomm_frag=", dg_frag%icomm_frag, &
        " id_frag=", dg_frag%id_frag, " isize_frag=", dg_frag%isize_frag, " ifrag_group=", dg_frag%ifrag_group
      write(*,'(1x,a,l1,a,a)') "        density trace: parallel-mode orbital=", dg_frag%parallel_mode_orbital, &
        " mode=", trim(dg_frag%parallel_mode)
      write(*,'(1x,a,l1,a,l1,a,l1,a,l1,a,l1)') "        density trace: mixed-state use=", use_mixed_density, &
        " ready=", dg_frag%mixed_basis_ready, " tr=", allocated(dg_frag%mixed_transform), &
        " coef_mix=", allocated(dg_frag%coef_mix), " dim=", allocated(dg_frag%mixed_basis_dim)
      flush(6)
    end if
    if (n_pw == 0 .and. enable_density_reconstruct_trace) then
      allocate(overlap_probe(dg_frag%n_mat_max, dg_frag%n_mat_max))
      allocate(overlap_vec(dg_frag%n_mat_max))
      overlap_probe(:, :) = (0.0d0, 0.0d0)
      overlap_vec(:) = (0.0d0, 0.0d0)
      basis_sdiag_probe(:, :, :) = 0.0d0
      call gather_full_coef_view(dg_frag, 1, dg_frag%n_mat_max, nocc_cache, coef_probe_full, coef_probe_pw)
      call copy_overlap_operator_to_dense(dg_frag, 1, .false., overlap_probe)
      do io = 1, min(3, dg_frag%n_mat_max)
        overlap_diag_probe(io) = real(overlap_probe(io, io), kind=8)
      end do
      do io = 1, min(3, dg_frag%nocc_spin(1))
        coef_state_probe(io) = real(sum(conjg(coef_probe_full(:, io)) * coef_probe_full(:, io)), kind=8)
        overlap_vec(:) = matmul(overlap_probe, coef_probe_full(:, io))
        overlap_state_probe(io) = real(sum(conjg(coef_probe_full(:, io)) * overlap_vec(:)), kind=8)
        if (dg_frag%n_frag >= 1) then
          do ix = 1, dg_frag%n_basis(1, 1)
            ig_i = dg_frag%index_basis(ix, 1, 1)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            do iy = 1, dg_frag%n_basis(1, 1)
              if (dg_frag%index_basis(iy, 1, 1) < 1 .or. dg_frag%index_basis(iy, 1, 1) > dg_frag%n_mat_max) cycle
              overlap_self_probe(1, io) = overlap_self_probe(1, io) + real(conjg(coef_probe_full(ig_i, io)) * &
                overlap_probe(ig_i, dg_frag%index_basis(iy, 1, 1)) * coef_probe_full(dg_frag%index_basis(iy, 1, 1), io), kind=8)
            end do
          end do
        end if
        if (dg_frag%n_frag >= 2) then
          do ix = 1, dg_frag%n_basis(2, 1)
            ig_i = dg_frag%index_basis(ix, 2, 1)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            do iy = 1, dg_frag%n_basis(2, 1)
              if (dg_frag%index_basis(iy, 2, 1) < 1 .or. dg_frag%index_basis(iy, 2, 1) > dg_frag%n_mat_max) cycle
              overlap_self_probe(2, io) = overlap_self_probe(2, io) + real(conjg(coef_probe_full(ig_i, io)) * &
                overlap_probe(ig_i, dg_frag%index_basis(iy, 2, 1)) * coef_probe_full(dg_frag%index_basis(iy, 2, 1), io), kind=8)
            end do
          end do
          do ix = 1, dg_frag%n_basis(1, 1)
            ig_i = dg_frag%index_basis(ix, 1, 1)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            do iy = 1, dg_frag%n_basis(2, 1)
              if (dg_frag%index_basis(iy, 2, 1) < 1 .or. dg_frag%index_basis(iy, 2, 1) > dg_frag%n_mat_max) cycle
              overlap_cross_probe(io) = overlap_cross_probe(io) + 2.0d0 * real(conjg(coef_probe_full(ig_i, io)) * &
                overlap_probe(ig_i, dg_frag%index_basis(iy, 2, 1)) * coef_probe_full(dg_frag%index_basis(iy, 2, 1), io), kind=8)
            end do
          end do
        end if
      end do
      do io = 1, dg_frag%nocc_spin(1)
        occ_factor = 1.0d0
        if (allocated(system%rocc)) then
          if (io <= size(system%rocc, 1) .and. 1 <= size(system%rocc, 3)) then
            occ_factor = max(0.0d0, system%rocc(io, 1, 1))
          end if
        end if
        if (occ_factor <= 0.0d0) cycle
        overlap_vec(:) = matmul(overlap_probe, coef_probe_full(:, io))
        overlap_elec_probe = overlap_elec_probe + occ_factor * real(sum(conjg(coef_probe_full(:, io)) * overlap_vec(:)), kind=8)
      end do
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        do ispin = 1, system%nspin
          nprobe_cols = min(3, dg_frag%n_basis(ifrag, ispin))
          do iprobe = 1, dg_frag%n_basis(ifrag, ispin)
            ig_i = dg_frag%index_basis(iprobe, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            do io = 1, dg_frag%n_basis(ifrag, ispin)
              if (dg_frag%index_basis(io, ifrag, ispin) < 1 .or. dg_frag%index_basis(io, ifrag, ispin) > dg_frag%n_mat_max) cycle
              basis_frag_metric_total(iprobe, io, ispin, i_local) = &
                real(overlap_probe(ig_i, dg_frag%index_basis(io, ifrag, ispin)), kind=8)
            end do
          end do
          do iprobe = 1, nprobe_cols
            ig_i = dg_frag%index_basis(iprobe, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            basis_sdiag_probe(iprobe, ispin, i_local) = real(overlap_probe(ig_i, ig_i), kind=8)
            do io = 1, nprobe_cols
              if (dg_frag%index_basis(io, ifrag, ispin) < 1 .or. dg_frag%index_basis(io, ifrag, ispin) > dg_frag%n_mat_max) cycle
              basis_smat_probe(iprobe, io, ispin, i_local) = real(overlap_probe(ig_i, dg_frag%index_basis(io, ifrag, ispin)), kind=8)
            end do
          end do
        end do
      end do
      if (enable_density_reconstruct_trace .and. dg_frag%id == 0) then
        write(*,'(1x,a,3(1pe12.4,1x),a,1pe12.4,a,3(1pe12.4,1x),a,3(1pe12.4,1x))') &
          "        density overlap probe: states=", &
          overlap_state_probe(1), overlap_state_probe(2), overlap_state_probe(3), " total=", overlap_elec_probe, &
          " c2=", coef_state_probe(1), coef_state_probe(2), coef_state_probe(3), &
          " sdiag=", overlap_diag_probe(1), overlap_diag_probe(2), overlap_diag_probe(3)
        flush(6)
        if (dg_frag%n_frag >= 2) then
          write(*,'(1x,a,3(1pe12.4,1x),a,3(1pe12.4,1x),a,3(1pe12.4,1x))') &
            "        density overlap split: self1=", overlap_self_probe(1,1), overlap_self_probe(1,2), overlap_self_probe(1,3), &
            " self2=", overlap_self_probe(2,1), overlap_self_probe(2,2), overlap_self_probe(2,3), &
            " cross=", overlap_cross_probe(1), overlap_cross_probe(2), overlap_cross_probe(3)
          flush(6)
        end if
      end if
      deallocate(coef_probe_full, coef_probe_pw, overlap_probe, overlap_vec)
    end if
    allocate(coef_re_frag(nbf_max, max(1, nocc_cache)))
    allocate(coef_im_frag(nbf_max, max(1, nocc_cache)))
    allocate(D_partial_re(nbf_max, nbf_max))
    allocate(D_frag_re(nbf_max, nbf_max, system%nspin))
    allocate(coef_re_full(nbf_max, max(1, nocc_cache), system%nspin))
    allocate(coef_im_full(nbf_max, max(1, nocc_cache), system%nspin))
    allocate(occ_cache(max(1, nocc_cache)), occ_sqrt_cache(max(1, nocc_cache)))
    allocate(coef_c_full(nbf_max, max(1, nocc_cache)), coef_c_frag(nbf_max, max(1, nocc_cache)))
    call cpu_time(t_setup0)
    phi_lb1 = lbound(dg_frag%phi_frag, 1)
    phi_lb2 = lbound(dg_frag%phi_frag, 2)
    phi_lb3 = lbound(dg_frag%phi_frag, 3)
    phi_ub1 = ubound(dg_frag%phi_frag, 1)
    phi_ub2 = ubound(dg_frag%phi_frag, 2)
    phi_ub3 = ubound(dg_frag%phi_frag, 3)
    phi_lg1 = dg_frag%lgnum_total(1)
    phi_lg2 = dg_frag%lgnum_total(2)
    phi_lg3 = dg_frag%lgnum_total(3)
    if (enable_density_phi_block_cache) then
      need_phi_cache_alloc = (.not. allocated(dg_frag%density_phi_block_cache))
      need_phi_count_alloc = (.not. allocated(dg_frag%density_phi_block_count))
      need_phi_cache_invalid = (.not. dg_frag%density_phi_block_cache_valid)
      need_phi_cache_resize = dg_frag%density_phi_block_size /= grid_block_size
      if (need_phi_cache_alloc .or. need_phi_count_alloc .or. need_phi_cache_invalid .or. need_phi_cache_resize) then
        if (allocated(dg_frag%density_phi_block_cache)) deallocate(dg_frag%density_phi_block_cache)
        if (allocated(dg_frag%density_phi_block_count)) deallocate(dg_frag%density_phi_block_count)
        nblocks_max = max(1, maxval(dg_frag%density_block_nblocks))
        allocate(dg_frag%density_phi_block_cache(grid_block_size, dg_frag%nstate_frag, nblocks_max, ifrag_count))
        allocate(dg_frag%density_phi_block_count(ifrag_count))
        dg_frag%density_phi_block_cache(:, :, :, :) = 0.0d0
        dg_frag%density_phi_block_count(:) = 0
        i_local = 0
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = i_local + 1
          local_grid_count = dg_frag%density_grid_point_count(i_local)
          dg_frag%density_phi_block_count(i_local) = dg_frag%density_block_nblocks(i_local)
          do block_cache_idx = 1, dg_frag%density_phi_block_count(i_local)
            igrid0 = 1 + (block_cache_idx - 1) * grid_block_size
            npt_cache = min(grid_block_size, local_grid_count - igrid0 + 1)
!$omp parallel do private(igrid, ixg, iyg, izg, bx, by, bz, istate_frag) schedule(static)
            do igrid = 1, npt_cache
              ixg = dg_frag%density_grid_points(igrid0 + igrid - 1, i_local)%ixg
              iyg = dg_frag%density_grid_points(igrid0 + igrid - 1, i_local)%iyg
              izg = dg_frag%density_grid_points(igrid0 + igrid - 1, i_local)%izg
              bx = map_global_to_phi_box_coord_ham(ixg, phi_lb1, phi_ub1, phi_lg1)
              by = map_global_to_phi_box_coord_ham(iyg, phi_lb2, phi_ub2, phi_lg2)
              bz = map_global_to_phi_box_coord_ham(izg, phi_lb3, phi_ub3, phi_lg3)
              if (bx == 0 .or. by == 0 .or. bz == 0) cycle
              do istate_frag = 1, dg_frag%nstate_frag
                dg_frag%density_phi_block_cache(igrid, istate_frag, block_cache_idx, i_local) = &
                  dg_frag%phi_frag(bx, by, bz, istate_frag, i_local)
              end do
            end do
!$omp end parallel do
          end do
        end do
        dg_frag%density_phi_block_size = grid_block_size
        dg_frag%density_phi_block_cache_valid = .true.
        rebuilt_phi_block_cache = .true.
      end if
    else
      if (allocated(dg_frag%density_phi_block_cache)) deallocate(dg_frag%density_phi_block_cache)
      if (allocated(dg_frag%density_phi_block_count)) deallocate(dg_frag%density_phi_block_count)
      dg_frag%density_phi_block_size = grid_block_size
      dg_frag%density_phi_block_cache_valid = .false.
    end if
    call cpu_time(t_setup1)
    time_project_phi_block_build = time_project_phi_block_build + (t_setup1 - t_setup0)
    if (rebuilt_phi_block_cache) time_cache_phi_block_refresh = time_cache_phi_block_refresh + (t_setup1 - t_setup0)
    if (enable_density_reconstruct_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a,l1,a,l1,a,l1,a,l1,a,l1,a,l1,a,1pe12.4)') "        density cache: phi_block enabled=", &
        enable_density_phi_block_cache, " rebuilt=", rebuilt_phi_block_cache, " alloc=", need_phi_cache_alloc, " count_alloc=", need_phi_count_alloc, &
        " invalid=", need_phi_cache_invalid, " resize=", need_phi_cache_resize, " dt=", time_cache_phi_block_refresh
      flush(6)
    end if
    call cpu_time(t_project0)
      i_local = 0
      block_idx_global = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1

        nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
        ngrid = nxyz(1) * nxyz(2) * nxyz(3)
        nblocks_ifrag = dg_frag%density_block_nblocks(i_local)
        first_block_offset = dg_frag%density_block_first_offset(i_local)
        block_step_blocks = dg_frag%density_block_step(i_local)
        valid_basis_count_spin(:) = 0
        basis_gid_spin(:, :) = 0
        valid_basis_ids_spin(:, :) = 0
        coef_re_full(:, :, :) = 0.0d0
        coef_im_full(:, :, :) = 0.0d0
        do ispin = 1, system%nspin
          nbf = dg_frag%n_basis(ifrag, ispin)
          if (nbf <= 0) cycle
          do istate_frag = 1, nbf
            basis_gid(istate_frag) = dg_frag%index_basis(istate_frag, ifrag, ispin)
            if (basis_gid(istate_frag) < 1 .or. basis_gid(istate_frag) > dg_frag%n_mat_max) cycle
            valid_basis_count_spin(ispin) = valid_basis_count_spin(ispin) + 1
            basis_gid_spin(istate_frag, ispin) = basis_gid(istate_frag)
            valid_basis_ids_spin(valid_basis_count_spin(ispin), ispin) = istate_frag
          end do
        end do
        if (use_mixed_density) then
          call cpu_time(t_setup0)
          n_basis_mix_spin(:) = 0
          do ispin = 1, system%nspin
            nbf = dg_frag%n_basis(ifrag, ispin)
            n_basis_mix = min(dg_frag%mixed_basis_dim(ispin), max_mixed_basis)
            n_basis_mix_spin(ispin) = n_basis_mix
            if (nbf <= 0 .or. n_basis_mix <= 0) cycle
            transform_frag_spin(1:nbf, 1:n_basis_mix, ispin) = (0.0d0, 0.0d0)
            do istate_frag = 1, nbf
              ig_i = dg_frag%index_basis(istate_frag, ifrag, ispin)
              if (ig_i < 1 .or. ig_i > n_frag) cycle
              transform_frag_spin(istate_frag, 1:n_basis_mix, ispin) = dg_frag%mixed_transform(ig_i, 1:n_basis_mix, ispin)
            end do
            if (n_pw > 0) then
              transform_pw_spin(1:n_pw, 1:n_basis_mix, ispin) = dg_frag%mixed_transform(n_frag+1:n_tot, 1:n_basis_mix, ispin)
              ! Phase-fixing test: rotate each mode m so max-abs PW component is real-positive
              if (enable_fp_phase_fix) then
                block
                  integer :: im_pf, ipw_pf, ipw_max
                  real(8) :: abs_max_pf
                  complex(8) :: phase_fix
                  do im_pf = 1, n_basis_mix
                    ipw_max = 1
                    abs_max_pf = abs(transform_pw_spin(1, im_pf, ispin))
                    do ipw_pf = 2, n_pw
                      if (abs(transform_pw_spin(ipw_pf, im_pf, ispin)) > abs_max_pf) then
                        abs_max_pf = abs(transform_pw_spin(ipw_pf, im_pf, ispin))
                        ipw_max = ipw_pf
                      end if
                    end do
                    if (abs_max_pf > 1.0d-30) then
                      phase_fix = conjg(transform_pw_spin(ipw_max, im_pf, ispin)) / abs_max_pf
                      transform_pw_spin(1:n_pw, im_pf, ispin) = transform_pw_spin(1:n_pw, im_pf, ispin) * phase_fix
                      transform_frag_spin(1:nbf, im_pf, ispin) = transform_frag_spin(1:nbf, im_pf, ispin) * phase_fix
                    end if
                  end do
                end block
              end if
            end if

            if (enable_tf_occmap_only .and. ispin == 1 .and. trace_occmap_itt .and. &
                nocc_spin > 0 .and. n_basis_mix > 0 .and. ifrag <= 2) then
              block
                complex(8), allocatable :: ff_occ_amp_min(:)
                integer :: iocc_min, io_ff_min, im_ff_min, gid_i_min
                integer :: gid_focus_slot_min, i_focus_min
                integer :: ff_audit_unit_min, ff_audit_ios_min
                integer :: m_mode_min, n_mode_max_min, m_peak_if2_min, m_second_if2_min
                integer :: m_ref_it1_min, ref_valid_it1_min, same_m_ref_it1_min
                integer :: n_ref_gt_min, ref_in_top5_min, ref_rank_min
                integer :: trans_ref_min, trans_cur_min
                integer :: i_step_min, track_step_idx_min
                integer, allocatable :: occmap_transition_counts_min(:,:,:)
                real(8), allocatable :: occmap_transition_ratio_sum_min(:,:,:), occmap_transition_margin_sum_min(:,:,:)
                real(8) :: occ_factor_min, if2_peak_val_min, if2_second_val_min
                real(8) :: ref_if2_val_min, cur_ref_ratio_min, ref_abs_min, margin_abs_min, margin_rel_min
                real(8) :: top12_ratio_min
                complex(8) :: tf_mode_amp_min

                cmp_slot = 0
                if (ifrag == 1) cmp_slot = 1
                if (ifrag == 2) cmp_slot = 2

                if (cmp_slot > 0) then
                  allocate(ff_occ_amp_min(nbf))
                  do iocc_min = 1, nocc_spin
                    occ_factor_min = 1.0d0
                    if (allocated(system%rocc)) then
                      if (iocc_min <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
                        occ_factor_min = max(0.0d0, system%rocc(iocc_min, 1, ispin))
                      end if
                    end if
                    if (occ_factor_min <= 0.0d0) cycle

                    ff_occ_amp_min(1:nbf) = matmul(transform_frag_spin(1:nbf, 1:n_basis_mix, ispin), coef_mix_spin(1:n_basis_mix, iocc_min, ispin))
                    do io_ff_min = 1, nbf
                      gid_i_min = dg_frag%index_basis(io_ff_min, ifrag, ispin)
                      if (gid_i_min < 1 .or. gid_i_min > n_frag) cycle
                      gid_focus_slot_min = 0
                      do i_focus_min = 1, n_focus_gid
                        if (gid_i_min == focus_gid_ids(i_focus_min)) then
                          gid_focus_slot_min = i_focus_min
                          exit
                        end if
                      end do
                      if (gid_focus_slot_min <= 0) cycle
                      do im_ff_min = 1, n_basis_mix
                        tf_mode_amp_min = transform_frag_spin(io_ff_min, im_ff_min, ispin) * coef_mix_spin(im_ff_min, iocc_min, ispin)
                        cmp_tf_gid_mode_full(iocc_min, im_ff_min, gid_focus_slot_min, cmp_slot) = cmp_tf_gid_mode_full(iocc_min, im_ff_min, gid_focus_slot_min, cmp_slot) + &
                          occ_factor_min * real(conjg(ff_occ_amp_min(io_ff_min)) * tf_mode_amp_min, kind=8)
                      end do
                    end do
                  end do
                  deallocate(ff_occ_amp_min)
                end if

                call comm_summation(cmp_tf_gid_mode_full, cmp_tf_gid_mode_full_global, size(cmp_tf_gid_mode_full), dg_frag%icomm)
                if (dg_frag%is_frag_root .and. cmp_slot == 2) then
                  ff_audit_unit_min = -1
                  open(newunit=ff_audit_unit_min, file='ff_tf_mode_interference.log', status='unknown', position='append', &
                    action='write', iostat=ff_audit_ios_min)
                  if (ff_audit_ios_min == 0) then
                    write(ff_audit_unit_min,*) 'ff-tf-header-lite: itt=', itt_tag, ' gidA=', focus_gid_ids(1), ' gidB=', focus_gid_ids(2)
                    n_mode_max_min = size(cmp_tf_gid_mode_full_global, 2)
                    allocate(occmap_transition_counts_min(n_mode_max_min, n_mode_max_min, n_focus_gid))
                    allocate(occmap_transition_ratio_sum_min(n_mode_max_min, n_mode_max_min, n_focus_gid))
                    allocate(occmap_transition_margin_sum_min(n_mode_max_min, n_mode_max_min, n_focus_gid))
                    occmap_transition_counts_min(:, :, :) = 0
                    occmap_transition_ratio_sum_min(:, :, :) = 0.0d0
                    occmap_transition_margin_sum_min(:, :, :) = 0.0d0
                    track_step_idx_min = 0
                    if (n_occmap_trace_itt > 0) then
                      do i_step_min = 1, min(n_occmap_trace_itt, max_occmap_track_step)
                        if (occmap_trace_itts(i_step_min) == itt_tag) then
                          track_step_idx_min = i_step_min
                          exit
                        end if
                      end do
                    else
                      if (itt_tag >= 1 .and. itt_tag <= max_occmap_track_step) track_step_idx_min = itt_tag
                    end if
                    if (track_step_idx_min > 0) then
                      cmp_tf_track_nstep = max(cmp_tf_track_nstep, track_step_idx_min)
                      cmp_tf_track_itt(track_step_idx_min) = itt_tag
                    end if
                    do iocc_min = 1, min(nocc_spin, size(cmp_tf_gid_mode_full_global, 1))
                      do i_focus_min = 1, n_focus_gid
                        m_peak_if2_min = 0
                        m_second_if2_min = 0
                        if2_peak_val_min = 0.0d0
                        if2_second_val_min = 0.0d0
                        do m_mode_min = 1, n_mode_max_min
                          if (abs(cmp_tf_gid_mode_full_global(iocc_min, m_mode_min, i_focus_min, 2)) > abs(if2_peak_val_min)) then
                            if2_second_val_min = if2_peak_val_min
                            m_second_if2_min = m_peak_if2_min
                            if2_peak_val_min = cmp_tf_gid_mode_full_global(iocc_min, m_mode_min, i_focus_min, 2)
                            m_peak_if2_min = m_mode_min
                          else if (abs(cmp_tf_gid_mode_full_global(iocc_min, m_mode_min, i_focus_min, 2)) > abs(if2_second_val_min)) then
                            if2_second_val_min = cmp_tf_gid_mode_full_global(iocc_min, m_mode_min, i_focus_min, 2)
                            m_second_if2_min = m_mode_min
                          end if
                        end do
                        if (iocc_min <= max_occmap_ref_occ) then
                          if (itt_tag == 1 .and. m_peak_if2_min > 0) &
                            cmp_tf_m2_ref_static(iocc_min, i_focus_min) = m_peak_if2_min
                          m_ref_it1_min = cmp_tf_m2_ref_static(iocc_min, i_focus_min)
                        else
                          m_ref_it1_min = 0
                        end if
                        ref_valid_it1_min = 0
                        if (m_ref_it1_min > 0) ref_valid_it1_min = 1
                        same_m_ref_it1_min = 0
                        if (m_ref_it1_min > 0 .and. m_ref_it1_min == m_peak_if2_min) same_m_ref_it1_min = 1
                        ref_if2_val_min = 0.0d0
                        if (m_ref_it1_min > 0 .and. m_ref_it1_min <= n_mode_max_min) then
                          ref_if2_val_min = cmp_tf_gid_mode_full_global(iocc_min, m_ref_it1_min, i_focus_min, 2)
                        end if
                        cur_ref_ratio_min = 0.0d0
                        if (abs(ref_if2_val_min) > 1.0d-30) cur_ref_ratio_min = if2_peak_val_min / ref_if2_val_min
                        margin_abs_min = if2_peak_val_min - ref_if2_val_min
                        margin_rel_min = margin_abs_min / max(abs(if2_peak_val_min), 1.0d-30)
                        ref_in_top5_min = 0
                        ref_rank_min = 0
                        if (m_ref_it1_min > 0 .and. m_ref_it1_min <= n_mode_max_min) then
                          ref_abs_min = abs(cmp_tf_gid_mode_full_global(iocc_min, m_ref_it1_min, i_focus_min, 2))
                          n_ref_gt_min = 0
                          do m_mode_min = 1, n_mode_max_min
                            if (abs(cmp_tf_gid_mode_full_global(iocc_min, m_mode_min, i_focus_min, 2)) > ref_abs_min) n_ref_gt_min = n_ref_gt_min + 1
                          end do
                          if (n_ref_gt_min < 5) then
                            ref_in_top5_min = 1
                            ref_rank_min = n_ref_gt_min + 1
                          end if
                        end if
                        top12_ratio_min = 0.0d0
                        if (abs(if2_second_val_min) > 1.0d-30) then
                          top12_ratio_min = if2_peak_val_min / if2_second_val_min
                        end if
                        if (ref_valid_it1_min == 1 .and. m_peak_if2_min > 0 .and. m_peak_if2_min <= n_mode_max_min) then
                          occmap_transition_counts_min(m_ref_it1_min, m_peak_if2_min, i_focus_min) = &
                            occmap_transition_counts_min(m_ref_it1_min, m_peak_if2_min, i_focus_min) + 1
                          occmap_transition_ratio_sum_min(m_ref_it1_min, m_peak_if2_min, i_focus_min) = &
                            occmap_transition_ratio_sum_min(m_ref_it1_min, m_peak_if2_min, i_focus_min) + cur_ref_ratio_min
                          occmap_transition_margin_sum_min(m_ref_it1_min, m_peak_if2_min, i_focus_min) = &
                            occmap_transition_margin_sum_min(m_ref_it1_min, m_peak_if2_min, i_focus_min) + margin_rel_min
                        end if
                        if (track_step_idx_min > 0 .and. iocc_min <= max_occmap_ref_occ) then
                          cmp_tf_track_m2(iocc_min, i_focus_min, track_step_idx_min) = m_peak_if2_min
                        end if
                        write(ff_audit_unit_min,*) 'ff-tf-occmap-lite: itt=', itt_tag, ' occ_id=', iocc_min, ' gid=', focus_gid_ids(i_focus_min), &
                          ' m2_peak=', m_peak_if2_min, ' if2_full_m2=', if2_peak_val_min, &
                          ' if1_full_at_m2=', merge(cmp_tf_gid_mode_full_global(iocc_min, m_peak_if2_min, i_focus_min, 1), 0.0d0, m_peak_if2_min > 0), &
                          ' ref1_m2=', m_ref_it1_min, ' ref_valid=', ref_valid_it1_min, ' same_ref1=', same_m_ref_it1_min, &
                          ' m2_second=', m_second_if2_min, ' if2_full_m2_second=', if2_second_val_min, ' top12_ratio=', top12_ratio_min
                        if (ref_valid_it1_min == 1 .and. same_m_ref_it1_min == 0) then
                          write(ff_audit_unit_min,*) 'ff-tf-occmap-lite-mismatch: itt=', itt_tag, ' occ_id=', iocc_min, ' gid=', focus_gid_ids(i_focus_min), &
                            ' ref1_m2=', m_ref_it1_min, ' m2_peak=', m_peak_if2_min, &
                            ' ref_rank=', ref_rank_min, ' cur_val=', if2_peak_val_min, ' ref_val=', ref_if2_val_min, &
                            ' cur_ref_ratio=', cur_ref_ratio_min, ' margin_abs=', margin_abs_min, ' margin_rel=', margin_rel_min, &
                            ' ref_in_top5=', ref_in_top5_min
                        end if
                      end do
                    end do
                    do i_focus_min = 1, n_focus_gid
                      do trans_ref_min = 1, n_mode_max_min
                        do trans_cur_min = 1, n_mode_max_min
                          if (occmap_transition_counts_min(trans_ref_min, trans_cur_min, i_focus_min) > 0) then
                            write(ff_audit_unit_min,*) 'ff-tf-occmap-transition: itt=', itt_tag, ' gid=', focus_gid_ids(i_focus_min), &
                              ' ref_m2=', trans_ref_min, ' cur_m2=', trans_cur_min, &
                              ' count=', occmap_transition_counts_min(trans_ref_min, trans_cur_min, i_focus_min), &
                              ' avg_cur_ref_ratio=', occmap_transition_ratio_sum_min(trans_ref_min, trans_cur_min, i_focus_min) / &
                                dble(max(1, occmap_transition_counts_min(trans_ref_min, trans_cur_min, i_focus_min))), &
                              ' avg_margin_rel=', occmap_transition_margin_sum_min(trans_ref_min, trans_cur_min, i_focus_min) / &
                                dble(max(1, occmap_transition_counts_min(trans_ref_min, trans_cur_min, i_focus_min)))
                          end if
                        end do
                      end do
                    end do
                    if (track_step_idx_min > 0) then
                      if ((n_occmap_trace_itt > 0 .and. track_step_idx_min == min(n_occmap_trace_itt, max_occmap_track_step)) .or. &
                          (n_occmap_trace_itt == 0 .and. itt_tag == 8)) then
                        do i_focus_min = 1, n_focus_gid
                          do iocc_min = 1, min(nocc_spin, max_occmap_ref_occ)
                            m_ref_it1_min = cmp_tf_m2_ref_static(iocc_min, i_focus_min)
                            if (m_ref_it1_min <= 0) cycle
                            write(ff_audit_unit_min,'(a,i0,a,i0,a,i0,a)',advance='no') 'ff-tf-occmap-track: gid=', focus_gid_ids(i_focus_min), &
                              ' occ_id=', iocc_min, ' ref_m2=', m_ref_it1_min, ' seq='
                            do i_step_min = 1, cmp_tf_track_nstep
                              if (cmp_tf_track_itt(i_step_min) <= 0) cycle
                              if (i_step_min > 1) write(ff_audit_unit_min,'(a)',advance='no') ','
                              write(ff_audit_unit_min,'(i0)',advance='no') cmp_tf_track_m2(iocc_min, i_focus_min, i_step_min)
                            end do
                            write(ff_audit_unit_min,*)
                          end do
                        end do
                      end if
                    end if
                    deallocate(occmap_transition_counts_min)
                    deallocate(occmap_transition_ratio_sum_min)
                    deallocate(occmap_transition_margin_sum_min)
                    flush(ff_audit_unit_min)
                    close(ff_audit_unit_min)
                  end if
                end if
              end block
            end if

            if ((enable_rho_mix_raw_trace .or. enable_fp_decomp_audit) .and. dg_frag%is_frag_root .and. trace_occmap_itt) then
              tr_ff_probe = 0.0d0
              tr_fp_probe = 0.0d0
              tr_pp_probe = 0.0d0
              fp_frob_probe = 0.0d0
              fp_maxabs_probe = 0.0d0
              trs_ff_probe = 0.0d0
              trs_fp_probe = 0.0d0
              trs_pp_probe = 0.0d0
              trs_total_probe = 0.0d0
              fp_weight_frac_probe = 0.0d0
              tr_fp_complex = (0.0d0, 0.0d0)
              if (allocated(d_raw_ff)) deallocate(d_raw_ff)
              allocate(d_raw_ff(nbf, nbf))
              d_raw_ff(:, :) = matmul(transform_frag_spin(1:nbf, 1:n_basis_mix, ispin), &
                matmul(density_mix(1:n_basis_mix, 1:n_basis_mix, ispin), &
                  conjg(transpose(transform_frag_spin(1:nbf, 1:n_basis_mix, ispin)))))
              do io = 1, nbf
                tr_ff_probe = tr_ff_probe + real(d_raw_ff(io, io), kind=8)
              end do
              if (n_pw > 0) then
                if (allocated(d_raw_fp)) deallocate(d_raw_fp)
                if (allocated(d_raw_pp)) deallocate(d_raw_pp)
                allocate(d_raw_fp(nbf, n_pw), d_raw_pp(n_pw, n_pw))
                d_raw_fp(:, :) = matmul(transform_frag_spin(1:nbf, 1:n_basis_mix, ispin), &
                  matmul(density_mix(1:n_basis_mix, 1:n_basis_mix, ispin), &
                    conjg(transpose(transform_pw_spin(1:n_pw, 1:n_basis_mix, ispin)))))
                d_raw_pp(:, :) = matmul(transform_pw_spin(1:n_pw, 1:n_basis_mix, ispin), &
                  matmul(density_mix(1:n_basis_mix, 1:n_basis_mix, ispin), &
                    conjg(transpose(transform_pw_spin(1:n_pw, 1:n_basis_mix, ispin)))))
                do io = 1, n_pw
                  tr_pp_probe = tr_pp_probe + real(d_raw_pp(io, io), kind=8)
                end do
                do io = 1, min(nbf, n_pw)
                  tr_fp_probe = tr_fp_probe + real(d_raw_fp(io, io), kind=8)
                end do
                fp_frob_probe = sqrt(max(0.0d0, sum(abs(d_raw_fp(:, :))**2)))
                fp_maxabs_probe = maxval(abs(d_raw_fp(:, :)))
              end if

              trs_ff_probe = 0.0d0
              do io = 1, nbf
                ig_i = dg_frag%index_basis(io, ifrag, ispin)
                if (ig_i < 1 .or. ig_i > n_frag) cycle
                do i_local = 1, nbf
                  istate_frag = dg_frag%index_basis(i_local, ifrag, ispin)
                  if (istate_frag < 1 .or. istate_frag > n_frag) cycle
                  if (allocated(dg_frag%S_mat_c)) then
                    tr_fp_complex = dg_frag%S_mat_c(ig_i, istate_frag, ispin)
                  else if (allocated(dg_frag%S_mat)) then
                    tr_fp_complex = cmplx(dg_frag%S_mat(ig_i, istate_frag, ispin), 0.0d0, kind=8)
                  else
                    tr_fp_complex = (0.0d0, 0.0d0)
                    if (ig_i == istate_frag) tr_fp_complex = (1.0d0, 0.0d0)
                  end if
                  trs_ff_probe = trs_ff_probe + real(tr_fp_complex * d_raw_ff(i_local, io), kind=8)
                end do
              end do
              if (n_pw > 0) then
                tr_fp_complex = (0.0d0, 0.0d0)
                top_fp_contrib_val(:) = 0.0d0
                top_fp_io(:) = 0
                top_fp_ipw(:) = 0
                top_fp_gid(:) = 0
                do io = 1, nbf
                  ig_i = dg_frag%index_basis(io, ifrag, ispin)
                  if (ig_i < 1 .or. ig_i > n_frag) cycle
                  do i_local = 1, n_pw
                    tr_fp_complex = tr_fp_complex + dg_frag%S_mat_frag_pw(ig_i, i_local, ispin) * conjg(d_raw_fp(io, i_local))
                    if ((enable_ifrag_compare_trace .or. enable_fp_decomp_audit) .and. ispin == 1) then
                      fp_pair_contrib = real(dg_frag%S_mat_frag_pw(ig_i, i_local, ispin) * conjg(d_raw_fp(io, i_local)), kind=8)
                      k_min = 1
                      do k_top = 2, n_top_fp_pairs
                        if (abs(top_fp_contrib_val(k_top)) < abs(top_fp_contrib_val(k_min))) k_min = k_top
                      end do
                      if (abs(fp_pair_contrib) > abs(top_fp_contrib_val(k_min))) then
                        top_fp_contrib_val(k_min) = fp_pair_contrib
                        top_fp_io(k_min) = io
                        top_fp_ipw(k_min) = i_local
                        top_fp_gid(k_min) = ig_i
                      end if
                    end if
                  end do
                end do
                trs_fp_probe = 2.0d0 * real(tr_fp_complex, kind=8)
                trs_pp_probe = tr_pp_probe
              end if
              trs_total_probe = trs_ff_probe + trs_fp_probe + trs_pp_probe
              if (abs(trs_total_probe) > 1.0d-20) fp_weight_frac_probe = trs_fp_probe / trs_total_probe

              write(*,'(1x,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
                '        rho_mix raw map: ifrag=', ifrag, ' ispin=', ispin, ' tr_ff=', tr_ff_probe, ' tr_fp_diag=', tr_fp_probe, &
                ' tr_pp=', tr_pp_probe, ' fp_frob=', fp_frob_probe, ' fp_max=', fp_maxabs_probe
              write(*,'(1x,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
                '        rho_mix fp-weighted: ifrag=', ifrag, ' ispin=', ispin, ' trs_ff=', trs_ff_probe, ' trs_fp=', trs_fp_probe, &
                ' trs_pp=', trs_pp_probe, ' trs_tot=', trs_total_probe, ' fp_frac=', fp_weight_frac_probe
              if (enable_tf_occmap_only .and. ispin == 1 .and. (itt_tag == 1 .or. itt_tag == 4 .or. itt_tag == 8) .and. &
                  nocc_spin > 0 .and. n_basis_mix > 0 .and. ifrag <= 2) then
                block
                  complex(8), allocatable :: ff_occ_amp(:)
                  integer :: iocc, io_ff, im_ff, gid_i
                  integer :: gid_focus_slot, i_focus
                  integer :: ff_audit_unit, ff_audit_ios
                  integer :: m_mode, n_mode_max, m_peak_if2
                  integer :: m_ref_it1, ref_valid_it1, same_m_ref_it1
                  real(8) :: occ_factor, if2_peak_val
                  complex(8) :: tf_mode_amp

                  cmp_slot = 0
                  if (ifrag == 1) cmp_slot = 1
                  if (ifrag == 2) cmp_slot = 2
                  if (cmp_slot > 0) then
                    allocate(ff_occ_amp(nbf))
                    do iocc = 1, nocc_spin
                      occ_factor = 1.0d0
                      if (allocated(system%rocc)) then
                        if (iocc <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
                          occ_factor = max(0.0d0, system%rocc(iocc, 1, ispin))
                        end if
                      end if
                      if (occ_factor <= 0.0d0) cycle

                      ff_occ_amp(1:nbf) = matmul(transform_frag_spin(1:nbf, 1:n_basis_mix, ispin), coef_mix_spin(1:n_basis_mix, iocc, ispin))
                      do io_ff = 1, nbf
                        gid_i = dg_frag%index_basis(io_ff, ifrag, ispin)
                        if (gid_i < 1 .or. gid_i > n_frag) cycle
                        gid_focus_slot = 0
                        do i_focus = 1, n_focus_gid
                          if (gid_i == focus_gid_ids(i_focus)) then
                            gid_focus_slot = i_focus
                            exit
                          end if
                        end do
                        if (gid_focus_slot <= 0) cycle
                        do im_ff = 1, n_basis_mix
                          tf_mode_amp = transform_frag_spin(io_ff, im_ff, ispin) * coef_mix_spin(im_ff, iocc, ispin)
                          cmp_tf_gid_mode_full(iocc, im_ff, gid_focus_slot, cmp_slot) = cmp_tf_gid_mode_full(iocc, im_ff, gid_focus_slot, cmp_slot) + &
                            occ_factor * real(conjg(ff_occ_amp(io_ff)) * tf_mode_amp, kind=8)
                        end do
                      end do
                    end do
                    deallocate(ff_occ_amp)
                  end if

                  call comm_summation(cmp_tf_gid_mode_full, cmp_tf_gid_mode_full_global, size(cmp_tf_gid_mode_full), dg_frag%icomm)
                  if (dg_frag%is_frag_root .and. cmp_slot == 2) then
                    ff_audit_unit = -1
                    open(newunit=ff_audit_unit, file='ff_tf_mode_interference.log', status='unknown', position='append', &
                      action='write', iostat=ff_audit_ios)
                    if (ff_audit_ios == 0) then
                      write(ff_audit_unit,*) 'ff-tf-header-lite: itt=', itt_tag, ' gidA=', focus_gid_ids(1), ' gidB=', focus_gid_ids(2)
                      n_mode_max = size(cmp_tf_gid_mode_full_global, 2)
                      do iocc = 1, min(nocc_spin, size(cmp_tf_gid_mode_full_global, 1))
                        do i_focus = 1, n_focus_gid
                          m_peak_if2 = 0
                          if2_peak_val = 0.0d0
                          do m_mode = 1, n_mode_max
                            if (abs(cmp_tf_gid_mode_full_global(iocc, m_mode, i_focus, 2)) > abs(if2_peak_val)) then
                              if2_peak_val = cmp_tf_gid_mode_full_global(iocc, m_mode, i_focus, 2)
                              m_peak_if2 = m_mode
                            end if
                          end do
                          if (iocc <= max_occmap_ref_occ) then
                            if (itt_tag == 1 .and. m_peak_if2 > 0) cmp_tf_m2_ref_static(iocc, i_focus) = m_peak_if2
                            m_ref_it1 = cmp_tf_m2_ref_static(iocc, i_focus)
                          else
                            m_ref_it1 = 0
                          end if
                          ref_valid_it1 = 0
                          if (m_ref_it1 > 0) ref_valid_it1 = 1
                          same_m_ref_it1 = 0
                          if (m_ref_it1 > 0 .and. m_ref_it1 == m_peak_if2) same_m_ref_it1 = 1
                          write(ff_audit_unit,*) 'ff-tf-occmap-lite: itt=', itt_tag, ' occ_id=', iocc, ' gid=', focus_gid_ids(i_focus), &
                            ' m2_peak=', m_peak_if2, ' if2_full_m2=', if2_peak_val, &
                            ' if1_full_at_m2=', merge(cmp_tf_gid_mode_full_global(iocc, m_peak_if2, i_focus, 1), 0.0d0, m_peak_if2 > 0), &
                            ' ref1_m2=', m_ref_it1, ' ref_valid=', ref_valid_it1, ' same_ref1=', same_m_ref_it1
                        end do
                      end do
                      flush(ff_audit_unit)
                      close(ff_audit_unit)
                    end if
                  end if
                end block
              end if

              if (enable_fp_decomp_audit .and. ispin == 1 .and. itt_tag <= 10 .and. nocc_spin > 0 .and. n_basis_mix > 0) then
                block
                  complex(8), allocatable :: ff_occ_amp(:), ff_occ_work(:)
                  real(8), allocatable :: ff_gid_part(:)
                  real(8), allocatable :: ff_gid_pre_by_gid(:), ff_gid_post_by_gid(:)
                  integer :: iocc, io_ff, il_ff, gid_i, gid_j, gid_dom, im_ff
                  integer :: gid_k, ff_audit_unit, ff_audit_ios
                  integer :: gid_focus_slot, i_focus
                  integer :: top_mode_if1_pre(n_top_tf_mode), top_mode_if2_pre(n_top_tf_mode)
                  integer :: top_mode_if1_full(n_top_tf_mode), top_mode_if2_full(n_top_tf_mode)
                  integer :: k_min_mode_if1_pre, k_min_mode_if2_pre
                  integer :: k_min_mode_if1_full, k_min_mode_if2_full
                  integer :: m_mode, n_mode_max
                  complex(8) :: tf_mode_amp
                  real(8) :: tf_full_io
                  real(8) :: top_mode_val_if1_pre(n_top_tf_mode), top_mode_val_if2_pre(n_top_tf_mode)
                  real(8) :: top_mode_val_if1_full(n_top_tf_mode), top_mode_val_if2_full(n_top_tf_mode)
                  real(8) :: top_mode_abs_if1(n_top_tf_mode), top_mode_abs_if2(n_top_tf_mode)
                  real(8) :: top_mode_phase_if1(n_top_tf_mode), top_mode_phase_if2(n_top_tf_mode)
                  real(8) :: phase_diff, overlap_count
                  real(8) :: cross_ph1, cross_ph2, cross_dph
                  real(8) :: occ_weight_m2
                  real(8) :: if2_peak_val
                  real(8) :: tf_gid_m2_if1, tf_gid_m2_if2
                  real(8) :: diag_est_m2_if1, diag_est_m2_if2
                  integer :: top_mode_if1(n_top_tf_mode), top_mode_if2(n_top_tf_mode)
                  integer :: same_mode_flag, overlap_i, overlap_j
                  integer :: m_peak_if2
                  integer :: m_ref_it1, ref_valid_it1, same_m_ref_it1
                  complex(8) :: mode_ovl_if1, mode_ovl_if2
                  complex(8) :: mode_ovl_if1_at_if2, mode_ovl_if2_ref
                  integer :: top_gid_if1(n_top_ff_gid), top_gid_if2(n_top_ff_gid)
                  integer :: top_gid_if1_post(n_top_ff_gid), top_gid_if2_post(n_top_ff_gid)
                  integer :: k_min_if1, k_min_if2, k_min_if1_post, k_min_if2_post
                  real(8) :: top_val_if1(n_top_ff_gid), top_val_if2(n_top_ff_gid)
                  real(8) :: top_val_if1_post(n_top_ff_gid), top_val_if2_post(n_top_ff_gid)
                  real(8) :: pre_part, post_part
                  real(8) :: ff_occ_trace, gid_dom_part, pre_io_part
                  allocate(ff_occ_amp(nbf), ff_occ_work(nbf), ff_gid_part(nbf))
                  allocate(ff_gid_pre_by_gid(max(1, n_frag)), ff_gid_post_by_gid(max(1, n_frag)))
                  ff_occ_amp(:) = (0.0d0, 0.0d0)
                  ff_occ_work(:) = (0.0d0, 0.0d0)
                  ff_gid_part(:) = 0.0d0
                  ff_gid_pre_by_gid(:) = 0.0d0
                  ff_gid_post_by_gid(:) = 0.0d0
                  cmp_slot = 0
                  if (ifrag == 1) cmp_slot = 1
                  if (ifrag == 2) cmp_slot = 2
                  do iocc = 1, nocc_spin
                    occ_factor = 1.0d0
                    if (allocated(system%rocc)) then
                      if (iocc <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
                        occ_factor = max(0.0d0, system%rocc(iocc, 1, ispin))
                      end if
                    end if
                    if (occ_factor <= 0.0d0) cycle
                    ff_occ_amp(1:nbf) = matmul(transform_frag_spin(1:nbf, 1:n_basis_mix, ispin), coef_mix_spin(1:n_basis_mix, iocc, ispin))
                    ff_occ_work(1:nbf) = (0.0d0, 0.0d0)
                    ff_gid_part(1:nbf) = 0.0d0
                    ff_gid_pre_by_gid(:) = 0.0d0
                    ff_gid_post_by_gid(:) = 0.0d0
                    do io_ff = 1, nbf
                      gid_i = dg_frag%index_basis(io_ff, ifrag, ispin)
                      if (gid_i < 1 .or. gid_i > n_frag) cycle
                      gid_focus_slot = 0
                      do i_focus = 1, n_focus_gid
                        if (gid_i == focus_gid_ids(i_focus)) then
                          gid_focus_slot = i_focus
                          exit
                        end if
                      end do
                      ! pre-transform_frag proxy on gid: diagonal-only propagation from mixed coefficients.
                      pre_io_part = 0.0d0
                      do im_ff = 1, n_basis_mix
                        pre_io_part = pre_io_part + abs(transform_frag_spin(io_ff, im_ff, ispin))**2 * &
                          abs(coef_mix_spin(im_ff, iocc, ispin))**2
                      end do
                      pre_part = occ_factor * pre_io_part
                      ff_gid_pre_by_gid(gid_i) = ff_gid_pre_by_gid(gid_i) + pre_part
                      if (cmp_slot > 0 .and. gid_focus_slot > 0 .and. (itt_tag == 1 .or. itt_tag == 4 .or. itt_tag == 8)) then
                        cmp_tf_gid_pre(iocc, gid_focus_slot, cmp_slot) = cmp_tf_gid_pre(iocc, gid_focus_slot, cmp_slot) + pre_part
                        tf_full_io = occ_factor * abs(ff_occ_amp(io_ff))**2
                        cmp_tf_gid_full(iocc, gid_focus_slot, cmp_slot) = cmp_tf_gid_full(iocc, gid_focus_slot, cmp_slot) + tf_full_io
                        cmp_tf_gid_int(iocc, gid_focus_slot, cmp_slot) = cmp_tf_gid_int(iocc, gid_focus_slot, cmp_slot) + (tf_full_io - pre_part)
                        do im_ff = 1, n_basis_mix
                          tf_mode_amp = transform_frag_spin(io_ff, im_ff, ispin) * coef_mix_spin(im_ff, iocc, ispin)
                          cmp_tf_gid_mode_pre(iocc, im_ff, gid_focus_slot, cmp_slot) = cmp_tf_gid_mode_pre(iocc, im_ff, gid_focus_slot, cmp_slot) + &
                            occ_factor * abs(transform_frag_spin(io_ff, im_ff, ispin))**2 * abs(coef_mix_spin(im_ff, iocc, ispin))**2
                          cmp_tf_gid_mode_full(iocc, im_ff, gid_focus_slot, cmp_slot) = cmp_tf_gid_mode_full(iocc, im_ff, gid_focus_slot, cmp_slot) + &
                            occ_factor * real(conjg(ff_occ_amp(io_ff)) * tf_mode_amp, kind=8)
                          cmp_tf_gid_mode_ovl(iocc, im_ff, gid_focus_slot, cmp_slot) = cmp_tf_gid_mode_ovl(iocc, im_ff, gid_focus_slot, cmp_slot) + &
                            occ_factor * conjg(ff_occ_amp(io_ff)) * tf_mode_amp
                        end do
                      end if
                      do il_ff = 1, nbf
                        gid_j = dg_frag%index_basis(il_ff, ifrag, ispin)
                        if (gid_j < 1 .or. gid_j > n_frag) cycle
                        if (allocated(dg_frag%S_mat_c)) then
                          tr_fp_complex = dg_frag%S_mat_c(gid_i, gid_j, ispin)
                        else if (allocated(dg_frag%S_mat)) then
                          tr_fp_complex = cmplx(dg_frag%S_mat(gid_i, gid_j, ispin), 0.0d0, kind=8)
                        else
                          tr_fp_complex = (0.0d0, 0.0d0)
                          if (gid_i == gid_j) tr_fp_complex = (1.0d0, 0.0d0)
                        end if
                        ff_occ_work(io_ff) = ff_occ_work(io_ff) + tr_fp_complex * ff_occ_amp(il_ff)
                      end do
                    end do
                    ff_occ_trace = 0.0d0
                    gid_dom = 0
                    gid_dom_part = 0.0d0
                    do io_ff = 1, nbf
                      gid_i = dg_frag%index_basis(io_ff, ifrag, ispin)
                      if (gid_i < 1 .or. gid_i > n_frag) cycle
                      ff_gid_part(io_ff) = occ_factor * real(conjg(ff_occ_amp(io_ff)) * ff_occ_work(io_ff), kind=8)
                      ff_gid_post_by_gid(gid_i) = ff_gid_post_by_gid(gid_i) + ff_gid_part(io_ff)
                      ff_occ_trace = ff_occ_trace + ff_gid_part(io_ff)
                      if (abs(ff_gid_part(io_ff)) > abs(gid_dom_part)) then
                        gid_dom_part = ff_gid_part(io_ff)
                        gid_dom = gid_i
                      end if
                    end do
                    if (cmp_slot > 0 .and. iocc <= size(cmp_ff_occ, 1)) then
                      cmp_ff_occ(iocc, cmp_slot) = ff_occ_trace
                      cmp_ff_dom_gid(iocc, cmp_slot) = gid_dom
                      cmp_ff_dom_gid_local(iocc, cmp_slot) = real(gid_dom, kind=8)
                      if (itt_tag == 1 .or. itt_tag == 4 .or. itt_tag == 8) then
                        cmp_ff_gid_pre(iocc, 1:n_frag, cmp_slot) = ff_gid_pre_by_gid(1:n_frag)
                        cmp_ff_gid_post(iocc, 1:n_frag, cmp_slot) = ff_gid_post_by_gid(1:n_frag)
                      end if
                    end if
                    write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4,a,i0,a,1pe12.4)') &
                      '        ff-decomp occ: itt=', itt_tag, ' ifrag=', ifrag, ' occ_id=', iocc, ' trs_ff_occ=', ff_occ_trace, &
                      ' dom_gid=', gid_dom, ' dom_part=', gid_dom_part
                  end do
                  if (cmp_slot > 0) ifrag_ff_occ_seen(cmp_slot) = .true.
                  if (itt_tag >= 8 .and. ifrag <= 2) then
                    call comm_summation(cmp_ff_occ, cmp_ff_occ_global, size(cmp_ff_occ), dg_frag%icomm)
                    call comm_summation(cmp_ff_dom_gid_local, cmp_ff_dom_gid_global, size(cmp_ff_dom_gid_local), dg_frag%icomm)
                    if (dg_frag%is_frag_root .and. cmp_slot == 2) then
                      do iocc = 1, min(nocc_cache, size(cmp_ff_occ_global, 1))
                        write(*,'(1x,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,i0,a,i0)') &
                          '        ff-decomp occ-diff: itt=', itt_tag, ' occ_id=', iocc, ' if1=', cmp_ff_occ_global(iocc, 1), &
                          ' if2=', cmp_ff_occ_global(iocc, 2), ' d=', cmp_ff_occ_global(iocc, 1) - cmp_ff_occ_global(iocc, 2), &
                          ' gid1=', nint(cmp_ff_dom_gid_global(iocc, 1)), ' gid2=', nint(cmp_ff_dom_gid_global(iocc, 2))
                      end do
                      flush(6)
                    end if
                  end if
                  if ((itt_tag == 1 .or. itt_tag == 4 .or. itt_tag == 8) .and. ifrag <= 2) then
                    call comm_summation(cmp_tf_gid_mode_full, cmp_tf_gid_mode_full_global, size(cmp_tf_gid_mode_full), dg_frag%icomm)
                    if (dg_frag%is_frag_root .and. cmp_slot == 2) then
                      ff_audit_unit = -1
                      open(newunit=ff_audit_unit, file='ff_tf_mode_interference.log', status='unknown', position='append', &
                        action='write', iostat=ff_audit_ios)
                      if (ff_audit_ios == 0) then
                        write(ff_audit_unit,*) 'ff-tf-header-lite: itt=', itt_tag, ' gidA=', focus_gid_ids(1), ' gidB=', focus_gid_ids(2)
                        n_mode_max = size(cmp_tf_gid_mode_full_global, 2)
                        do iocc = 1, min(nocc_spin, size(cmp_tf_gid_pre_global, 1))
                          do i_focus = 1, n_focus_gid
                            m_peak_if2 = 0
                            if2_peak_val = 0.0d0
                            do m_mode = 1, n_mode_max
                              if (abs(cmp_tf_gid_mode_full_global(iocc, m_mode, i_focus, 2)) > abs(if2_peak_val)) then
                                if2_peak_val = cmp_tf_gid_mode_full_global(iocc, m_mode, i_focus, 2)
                                m_peak_if2 = m_mode
                              end if
                            end do

                            if (iocc <= max_occmap_ref_occ) then
                              if (itt_tag == 1 .and. m_peak_if2 > 0) cmp_tf_m2_ref_static(iocc, i_focus) = m_peak_if2
                              m_ref_it1 = cmp_tf_m2_ref_static(iocc, i_focus)
                            else
                              m_ref_it1 = 0
                            end if
                            ref_valid_it1 = 0
                            if (m_ref_it1 > 0) ref_valid_it1 = 1
                            same_m_ref_it1 = 0
                            if (m_ref_it1 > 0 .and. m_ref_it1 == m_peak_if2) same_m_ref_it1 = 1

                            write(ff_audit_unit,*) 'ff-tf-occmap-lite: itt=', itt_tag, ' occ_id=', iocc, ' gid=', focus_gid_ids(i_focus), &
                              ' m2_peak=', m_peak_if2, ' if2_full_m2=', if2_peak_val, &
                              ' if1_full_at_m2=', merge(cmp_tf_gid_mode_full_global(iocc, m_peak_if2, i_focus, 1), 0.0d0, m_peak_if2 > 0), &
                              ' ref1_m2=', m_ref_it1, ' ref_valid=', ref_valid_it1, ' same_ref1=', same_m_ref_it1
                          end do
                        end do
                        flush(ff_audit_unit)
                        close(ff_audit_unit)
                      end if
                    end if
                  end if
                end block
              end if
              if (enable_fp_decomp_audit .and. ispin == 1 .and. itt_tag <= 10) then
                write(*,'(1x,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
                  '        fp-decomp total: itt=', itt_tag, ' ifrag=', ifrag, ' trs_ff=', trs_ff_probe, &
                  ' trs_fp=', trs_fp_probe, ' trs_pp=', trs_pp_probe, ' trs_tot=', trs_total_probe
              end if
              if (enable_ifrag_compare_trace .and. ispin == 1) then
                cmp_slot = 0
                if (ifrag == 1) cmp_slot = 1
                if (ifrag == 2) cmp_slot = 2
                if (cmp_slot > 0) then
                  cmp_trs_fp(cmp_slot) = trs_fp_probe
                  cmp_fp_frac(cmp_slot) = fp_weight_frac_probe
                  cmp_fp_frob(cmp_slot) = fp_frob_probe
                  cmp_fp_max(cmp_slot) = fp_maxabs_probe
                  cmp_send_pre(cmp_slot) = send_pre_block_raw
                  cmp_trs_ff(cmp_slot) = trs_ff_probe
                  cmp_trs_pp(cmp_slot) = trs_pp_probe
                  cmp_trs_tot(cmp_slot) = trs_total_probe
                  ifrag_fp_seen(cmp_slot) = .true.
                  ifrag_decomp_seen(cmp_slot) = .true.
                end if
                if (enable_fp_decomp_audit .and. cmp_slot > 0 .and. dg_frag%is_frag_root .and. itt_tag <= 10) then
                  write(*,'(1x,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
                    '        fp-decomp ifrag-total: itt=', itt_tag, ' ifrag=', ifrag, ' ff=', cmp_trs_ff(cmp_slot), &
                    ' fp=', cmp_trs_fp(cmp_slot), ' pp=', cmp_trs_pp(cmp_slot), ' total=', cmp_trs_tot(cmp_slot)
                  if (ifrag_decomp_seen(1) .and. ifrag_decomp_seen(2)) then
                    write(*,'(1x,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
                      '        fp-decomp ifrag-diff: itt=', itt_tag, ' d_ff=', cmp_trs_ff(1)-cmp_trs_ff(2), &
                      ' d_fp=', cmp_trs_fp(1)-cmp_trs_fp(2), ' d_pp=', cmp_trs_pp(1)-cmp_trs_pp(2), &
                      ' d_tot=', cmp_trs_tot(1)-cmp_trs_tot(2)
                  end if
                end if
                if (cmp_slot > 0 .and. dg_frag%is_frag_root .and. enable_ifrag_compare_trace) then
                  write(*,'(1x,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
                    '        ifrag-compare fp: itt=', itt_tag, ' ifrag=', ifrag, ' trs_fp=', cmp_trs_fp(cmp_slot), &
                    ' frac=', cmp_fp_frac(cmp_slot), ' frob=', cmp_fp_frob(cmp_slot), ' max=', cmp_fp_max(cmp_slot)
                  write(*,'(1x,a,i0,a,1pe12.4)') '        ifrag-compare fp-flow: ifrag=', ifrag, &
                    ' send_pre=', cmp_send_pre(cmp_slot)
                  ! --- top contributing (frag_state, pw) pairs sorted by |contribution| ---
                  do k_top = 1, n_top_fp_pairs - 1
                    do io_top = k_top + 1, n_top_fp_pairs
                      if (abs(top_fp_contrib_val(io_top)) > abs(top_fp_contrib_val(k_top))) then
                        fp_pair_contrib = top_fp_contrib_val(k_top)
                        top_fp_contrib_val(k_top) = top_fp_contrib_val(io_top)
                        top_fp_contrib_val(io_top) = fp_pair_contrib
                        ipw_top = top_fp_io(k_top); top_fp_io(k_top) = top_fp_io(io_top); top_fp_io(io_top) = ipw_top
                        ipw_top = top_fp_ipw(k_top); top_fp_ipw(k_top) = top_fp_ipw(io_top); top_fp_ipw(io_top) = ipw_top
                        ipw_top = top_fp_gid(k_top); top_fp_gid(k_top) = top_fp_gid(io_top); top_fp_gid(io_top) = ipw_top
                      end if
                    end do
                  end do
                  do k_top = 1, n_top_fp_pairs
                    if (top_fp_io(k_top) == 0) cycle
                    write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,1pe12.4)') &
                      '        ifrag-compare fp-top: ifrag=', ifrag, ' rank=', k_top, &
                      ' fstate=', top_fp_io(k_top), ' gid=', top_fp_gid(k_top), &
                      ' ipw=', top_fp_ipw(k_top), ' contrib2=', 2.0d0 * top_fp_contrib_val(k_top)
                    ! factor decomposition: S_mat_frag_pw and d_raw_fp separately
                    block
                      complex(8) :: s_fac, d_fac
                      s_fac = dg_frag%S_mat_frag_pw(top_fp_gid(k_top), top_fp_ipw(k_top), ispin)
                      d_fac = d_raw_fp(top_fp_io(k_top), top_fp_ipw(k_top))
                      write(*,'(1x,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
                        '          ifrag-compare fp-factors: ifrag=', ifrag, ' rank=', k_top, &
                        ' S_re=', real(s_fac,8), ' S_im=', aimag(s_fac), &
                        ' D_re=', real(d_fac,8), ' D_im=', aimag(d_fac)
                      if (enable_fp_decomp_audit .and. itt_tag <= 10) then
                        ikx_pw = nint(dg_frag%k_pw(1, top_fp_ipw(k_top)))
                        iky_pw = nint(dg_frag%k_pw(2, top_fp_ipw(k_top)))
                        ikz_pw = nint(dg_frag%k_pw(3, top_fp_ipw(k_top)))
                        g_abs = sqrt((dg_frag%k_pw(1, top_fp_ipw(k_top)) * dg_frag%hgs(1))**2 + &
                                     (dg_frag%k_pw(2, top_fp_ipw(k_top)) * dg_frag%hgs(2))**2 + &
                                     (dg_frag%k_pw(3, top_fp_ipw(k_top)) * dg_frag%hgs(3))**2)
                        dom_m = 1
                        dom_m_abs = -1.0d0
                        do io = 1, n_basis_mix
                          m_term = transform_frag_spin(top_fp_io(k_top), io, ispin) * density_mix(io, io, ispin) * &
                                   conjg(transform_pw_spin(top_fp_ipw(k_top), io, ispin))
                          m_term_abs = abs(m_term)
                          if (m_term_abs > dom_m_abs) then
                            dom_m_abs = m_term_abs
                            dom_m = io
                          end if
                        end do
                        m_term = transform_frag_spin(top_fp_io(k_top), dom_m, ispin) * density_mix(dom_m, dom_m, ispin) * &
                                 conjg(transform_pw_spin(top_fp_ipw(k_top), dom_m, ispin))
                        dom_sign = 0.0d0
                        if (real(m_term, kind=8) > 0.0d0) dom_sign = 1.0d0
                        if (real(m_term, kind=8) < 0.0d0) dom_sign = -1.0d0
                        write(*,'(1x,a,i0,a,i0,a,i0,a,3(i0,1x),a,1pe12.4,a,3(1pe11.3,1x),a,i0)') &
                          '          fp-decomp pw-id: ifrag=', ifrag, ' rank=', k_top, ' ipw=', top_fp_ipw(k_top), &
                          ' ik=', ikx_pw, iky_pw, ikz_pw, ' |G|=', g_abs, &
                          ' kpw=', dg_frag%k_pw(1, top_fp_ipw(k_top)), dg_frag%k_pw(2, top_fp_ipw(k_top)), &
                          dg_frag%k_pw(3, top_fp_ipw(k_top)), ' gpw=', n_frag + top_fp_ipw(k_top)
                        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
                          '          fp-decomp dom-map: ifrag=', ifrag, ' occ_id=', top_fp_io(k_top), ' gid=', top_fp_gid(k_top), &
                          ' ipw=', top_fp_ipw(k_top), ' m=', dom_m, ' sign=', dom_sign, ' mabs=', dom_m_abs, &
                          ' S_re=', real(s_fac,8), ' S_im=', aimag(s_fac), ' D_re=', real(d_fac,8), ' D_im=', aimag(d_fac)
                      end if
                      ! Drill into d_raw_fp = T_f * D_mix * T_pw^H:
                      ! for rank=1 only, print per-mixmode m: Tf(io,m), Tpw(ipw,m), Dmix(m,m)
                      if (k_top == 1) then
                        block
                          integer :: im, n_mix_pr, io_k, ipw_k
                          complex(8) :: tf_m, tpw_m, dm_m, partial
                          io_k  = top_fp_io(k_top)
                          ipw_k = top_fp_ipw(k_top)
                          n_mix_pr = min(n_basis_mix, 5)
                          if (io_k >= 1 .and. io_k <= size(transform_frag_spin,1) .and. &
                              ipw_k >= 1 .and. ipw_k <= size(transform_pw_spin,1)) then
                            do im = 1, n_mix_pr
                              if (im > size(transform_frag_spin,2) .or. im > size(transform_pw_spin,2) .or. &
                                  im > size(density_mix,1)) exit
                              tf_m   = transform_frag_spin(io_k,  im, ispin)
                              tpw_m  = transform_pw_spin(ipw_k, im, ispin)
                              dm_m   = density_mix(im, im, ispin)
                              partial = tf_m * dm_m * conjg(tpw_m)
                              write(*,'(1x,a,i0,a,i0,a,1pe11.3,a,1pe11.3,a,1pe11.3,a,1pe11.3,a,1pe11.3,a,1pe11.3,a,1pe11.3,a,1pe11.3)') &
                                '          ifrag-compare fp-dmix: ifrag=', ifrag, ' rank=1 m=', im, &
                                ' Tf_re=', real(tf_m,8),  ' Tf_im=', aimag(tf_m), &
                                ' Tpw_re=', real(tpw_m,8), ' Tpw_im=', aimag(tpw_m), &
                                ' Dm_re=', real(dm_m,8), ' Dm_im=', aimag(dm_m), &
                                ' p_re=', real(partial,8), ' p_im=', aimag(partial)
                            end do
                          end if
                        end block
                      end if
                    end block
                  end do
                  flush(6)
                end if
              end if
              flush(6)
            end if
          end do
          call cpu_time(t_setup1)
          time_project_setup = time_project_setup + (t_setup1 - t_setup0)
        end if
        ! --- D pre-pass: compute density matrix for all spins before block loop ---
        D_frag_re(:,:,:) = 0.0d0
        if (n_pw == 0) then
          do ispin = 1, system%nspin
            nbf = dg_frag%n_basis(ifrag, ispin)
            nocc_spin = dg_frag%nocc_spin(ispin)
            if (nbf <= 0 .or. nocc_spin <= 0) cycle
            state_charge_local(1:nocc_spin, ispin) = 0.0d0
            state_coeff_c2_local(1:nocc_spin, ispin) = 0.0d0
            state_psi2_raw_local(1:nocc_spin, ispin) = 0.0d0
            state_psi2_occ_local(1:nocc_spin, ispin) = 0.0d0
            state_psi2_dv_local(1:nocc_spin, ispin) = 0.0d0
            state_psi2_owned_local(1:nocc_spin, ispin) = 0.0d0
            state_import_core_local(1:nocc_spin, ispin) = 0.0d0
            probe_state_owned_local(1:nocc_spin, ispin) = 0.0d0
            probe_state_import_local(1:nocc_spin, ispin) = 0.0d0
            occ_cache(1:nocc_spin) = 1.0d0
            if (allocated(system%rocc)) then
              do io = 1, nocc_spin
                if (io <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
                  occ_cache(io) = max(0.0d0, system%rocc(io, 1, ispin))
                end if
              end do
            end if
            occ_sqrt_cache(1:nocc_spin) = sqrt(occ_cache(1:nocc_spin))
            if (enable_density_reconstruct_trace) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,3(1pe12.4,1x))') &
                "        density occupancy trace: rank=", dg_frag%id, " ifrag=", ifrag, &
                " ispin=", ispin, " nocc=", nocc_spin, " occ[1-3]=", occ_cache(1:min(3,nocc_spin))
              flush(6)
            end if
            valid_basis_count = valid_basis_count_spin(ispin)
            call cpu_time(t_dmat0)
            ! Step 3a: each rank contributes its owned rows, then icomm_frag assembles
            ! the complete coefficient matrix. After zero_nonowned_coefficients(),
            ! the fragment root no longer holds a full copy by itself.
            coef_re_frag(1:nbf_max, 1:nocc_spin) = 0.0d0
            coef_im_frag(1:nbf_max, 1:nocc_spin) = 0.0d0
            coef_c_full(1:nbf_max, 1:nocc_spin) = (0.0d0, 0.0d0)
!$omp parallel do private(io, idx_local, istate_frag) schedule(static)
            do io = 1, nocc_spin
              do idx_local = 1, valid_basis_count
                istate_frag = valid_basis_ids_spin(idx_local, ispin)
                coef_c_full(istate_frag, io) = dg_frag%coef(basis_gid_spin(istate_frag, ispin), io, ispin)
              end do
            end do
!$omp end parallel do
            if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
              coef_c_frag(1:nbf_max, 1:nocc_spin) = (0.0d0, 0.0d0)
              if (itt_tag == 1) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] BEFORE_COMM_SUM_DENSITY_COEF rank=', &
                  dg_frag%id, ' itt=', itt_tag, ' ifrag=', ifrag, ' ispin=', ispin, ' nbf=', nbf, &
                  ' coef_ref_ready=', dg_frag%coef_ref_ready
                flush(6)
              end if
              call comm_summation(coef_c_full(1:nbf, 1:nocc_spin), coef_c_frag(1:nbf, 1:nocc_spin), &
                                  nbf * nocc_spin, dg_frag%icomm_frag)
              if (itt_tag == 1) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] AFTER_COMM_SUM_DENSITY_COEF rank=', &
                  dg_frag%id, ' itt=', itt_tag, ' ifrag=', ifrag, ' ispin=', ispin, ' nbf=', nbf, &
                  ' coef_ref_ready=', dg_frag%coef_ref_ready
                flush(6)
              end if
              coef_c_full(1:nbf_max, 1:nocc_spin) = coef_c_frag(1:nbf_max, 1:nocc_spin)
            end if
            coef_re_full(1:nbf_max, 1:nocc_spin, ispin) = real(coef_c_full(1:nbf_max, 1:nocc_spin), kind=8)
            coef_im_full(1:nbf_max, 1:nocc_spin, ispin) = aimag(coef_c_full(1:nbf_max, 1:nocc_spin))
            if (enable_density_state_charge_trace) then
              do io = 1, nocc_spin
                occ_factor = occ_cache(io)
                if (occ_factor <= 0.0d0) cycle
                state_coeff_c2_local(io, ispin) = state_coeff_c2_local(io, ispin) + occ_factor * &
                  sum(coef_re_full(1:nbf, io, ispin) * coef_re_full(1:nbf, io, ispin) + &
                      coef_im_full(1:nbf, io, ispin) * coef_im_full(1:nbf, io, ispin))
              end do
            end if
            ! Step 3b: each rank computes weighted C C^H on its state slice
            nocc_per_rank_loc = (nocc_spin + dg_frag%isize_frag - 1) / dg_frag%isize_frag
            io_s_frag = dg_frag%id_frag * nocc_per_rank_loc + 1
            io_e_frag = min((dg_frag%id_frag + 1) * nocc_per_rank_loc, nocc_spin)
            coef_norm_probe = sum(coef_re_full(1:nbf, 1:nocc_spin, ispin)**2 + coef_im_full(1:nbf, 1:nocc_spin, ispin)**2)
            if (enable_density_reconstruct_trace) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,1pe12.4,a,i0,a,i0)') "        density frag probe: rank=", dg_frag%id, &
                " ifrag=", ifrag, " ispin=", ispin, " nocc=", nocc_spin, " coef_norm=", coef_norm_probe, &
                " io_s=", io_s_frag, " io_e=", io_e_frag
              flush(6)
            end if
            coef_map_local_probe(:, :) = 0.0d0
            coef_map_global_probe(:, :) = 0.0d0
            coef_map_diff_probe(:, :) = 0.0d0
            do iprobe = 1, min(3, nbf)
              if (basis_gid_spin(iprobe, ispin) <= 0) cycle
              do io = 1, min(3, nocc_spin)
                coef_map_local_probe(iprobe, io) = coef_re_full(iprobe, io, ispin)
                coef_map_global_probe(iprobe, io) = real(dg_frag%coef(basis_gid_spin(iprobe, ispin), io, ispin), kind=8)
                coef_map_diff_probe(iprobe, io) = coef_map_local_probe(iprobe, io) - coef_map_global_probe(iprobe, io)
              end do
            end do
            if (enable_density_reconstruct_trace) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,3(1pe12.4,1x),a,3(1pe12.4,1x),a,3(1pe12.4,1x))') &
                "        density coef-map probe: rank=", dg_frag%id, " ifrag=", ifrag, " ispin=", ispin, &
                " row1_local=", coef_map_local_probe(1,1), coef_map_local_probe(1,2), coef_map_local_probe(1,3), &
                " row1_global=", coef_map_global_probe(1,1), coef_map_global_probe(1,2), coef_map_global_probe(1,3), &
                " row1_diff=", coef_map_diff_probe(1,1), coef_map_diff_probe(1,2), coef_map_diff_probe(1,3)
              flush(6)
            end if
            nocc_loc = max(0, io_e_frag - io_s_frag + 1)

            D_partial_re(1:nbf_max, 1:nbf_max) = 0.0d0
            if (enable_density_reconstruct_trace) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                "        density D_partial_re before-calc: rank=", dg_frag%id, " ifrag=", ifrag, &
                " ispin=", ispin, " nbf=", nbf, " nocc_loc=", nocc_loc
              flush(6)
            end if
            if (nocc_loc > 0 .and. nbf > 0) then
              nbatch = 0
              do io = io_s_frag, io_e_frag
                nbatch = nbatch + 1
                coef_blk_re(1:nbf, nbatch) = occ_sqrt_cache(io) * coef_re_full(1:nbf, io, ispin)
                coef_blk_im(1:nbf, nbatch) = occ_sqrt_cache(io) * coef_im_full(1:nbf, io, ispin)
                if (nbatch == state_block_size .or. io == io_e_frag) then
                  call dsyrk('U', 'N', nbf, nbatch, 1.0d0, coef_blk_re, nbf_max, 1.0d0, D_partial_re, nbf_max)
                  call dsyrk('U', 'N', nbf, nbatch, 1.0d0, coef_blk_im, nbf_max, 1.0d0, D_partial_re, nbf_max)
                  nbatch = 0
                end if
              end do
              if (enable_density_reconstruct_trace) then
                D_partial_trace = sum(D_partial_re(1:nbf, 1:nbf))
                write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4)') &
                  "        density D_partial_re after-calc: rank=", dg_frag%id, " ifrag=", ifrag, &
                  " ispin=", ispin, " trace=", D_partial_trace
                flush(6)
              end if
            end if
            ! Step 3c: AllReduce partial D across icomm_frag
            if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
              if (itt_tag == 1) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] BEFORE_COMM_SUM_DENSITY_DMAT rank=', &
                  dg_frag%id, ' itt=', itt_tag, ' ifrag=', ifrag, ' ispin=', ispin, ' coef_ref_ready=', dg_frag%coef_ref_ready
                flush(6)
              end if
              call comm_summation(D_partial_re(1:nbf, 1:nbf), D_frag_re(1:nbf, 1:nbf, ispin), &
                                  nbf * nbf, dg_frag%icomm_frag)
              if (itt_tag == 1) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] AFTER_COMM_SUM_DENSITY_DMAT rank=', &
                  dg_frag%id, ' itt=', itt_tag, ' ifrag=', ifrag, ' ispin=', ispin, ' coef_ref_ready=', dg_frag%coef_ref_ready
                flush(6)
              end if
            else
              D_frag_re(1:nbf_max, 1:nbf_max, ispin) = D_partial_re(1:nbf_max, 1:nbf_max)
            end if

            if (nbf > nbf_max) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] density nbf overflow: rank=", dg_frag%id, &
                " ifrag=", ifrag, " ispin=", ispin, " nbf=", nbf
              flush(6)
              stop "DG-Fragment RT: density nbf exceeds nbf_max"
            end if
            if (valid_basis_count > nbf_max) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] density valid_basis_count overflow: rank=", dg_frag%id, &
                " ifrag=", ifrag, " ispin=", ispin, " valid=", valid_basis_count
              flush(6)
              stop "DG-Fragment RT: density valid_basis_count exceeds nbf_max"
            end if

            do idx_local = 1, valid_basis_count
              istate_frag = valid_basis_ids_spin(idx_local, ispin)
              if (istate_frag < 1 .or. istate_frag > nbf_max) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') "[FATAL] density basis index out of range: rank=", dg_frag%id, &
                  " ifrag=", ifrag, " ispin=", ispin, " idx_local=", idx_local, " istate=", istate_frag
                flush(6)
                stop "DG-Fragment RT: density basis index out of range"
              end if
              ig_i = basis_gid_spin(istate_frag, ispin)
              if (ig_i < 0 .or. ig_i > dg_frag%n_mat_max) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') "[FATAL] density global basis id out of range: rank=", dg_frag%id, &
                  " ifrag=", ifrag, " ispin=", ispin, " istate=", istate_frag, " gid=", ig_i
                flush(6)
                stop "DG-Fragment RT: density global basis id out of range"
              end if
            end do
            ! Step 3d: symmetrize (copy upper triangle to lower)
!$omp parallel do private(io, istate_frag) schedule(static)
            do io = 1, nbf
              do istate_frag = io + 1, nbf
                D_frag_re(istate_frag, io, ispin) = D_frag_re(io, istate_frag, ispin)
              end do
            end do
!$omp end parallel do
            ! n_pw=0: D_frag_re is authoritative; density_matrix_frag (complex) not updated in this path
            dg_frag%density_matrix_frag_valid(ispin, i_local) = .true.
            if (enable_density_reconstruct_trace .and. dg_frag%is_frag_root) then
              frag_trace_probe = 0.0d0
              frag_state_trace_probe(:) = 0.0d0
              do io = 1, nbf
                ig_i = basis_gid_spin(io, ispin)
                if (ig_i <= 0) cycle
                do istate_frag = 1, nbf
                  if (basis_gid_spin(istate_frag, ispin) <= 0) cycle
                  if (allocated(dg_frag%S_mat)) then
                    occ_factor = dg_frag%S_mat(ig_i, basis_gid_spin(istate_frag, ispin), ispin)
                  else if (allocated(dg_frag%S_mat_c)) then
                    occ_factor = real(dg_frag%S_mat_c(ig_i, basis_gid_spin(istate_frag, ispin), ispin), kind=8)
                  else
                    occ_factor = 0.0d0
                  end if
                  frag_trace_probe = frag_trace_probe + D_frag_re(io, istate_frag, ispin) * occ_factor
                end do
              end do
              do io = 1, min(3, nocc_spin)
                do idx_local = 1, valid_basis_count
                  ig_i = valid_basis_ids_spin(idx_local, ispin)
                  if (ig_i <= 0 .or. basis_gid_spin(ig_i, ispin) <= 0) cycle
                  do istate_frag = 1, valid_basis_count
                    if (valid_basis_ids_spin(istate_frag, ispin) <= 0) cycle
                    if (basis_gid_spin(valid_basis_ids_spin(istate_frag, ispin), ispin) <= 0) cycle
                    if (allocated(dg_frag%S_mat)) then
                      occ_factor = dg_frag%S_mat(basis_gid_spin(ig_i, ispin), &
                                                         basis_gid_spin(valid_basis_ids_spin(istate_frag, ispin), ispin), ispin)
                    else if (allocated(dg_frag%S_mat_c)) then
                      occ_factor = real(dg_frag%S_mat_c(basis_gid_spin(ig_i, ispin), &
                                                         basis_gid_spin(valid_basis_ids_spin(istate_frag, ispin), ispin), ispin), kind=8)
                    else
                      occ_factor = 0.0d0
                    end if
                    frag_state_trace_probe(io) = frag_state_trace_probe(io) + &
                      (coef_re_full(ig_i, io, ispin) * coef_re_full(valid_basis_ids_spin(istate_frag, ispin), io, ispin) + &
                       coef_im_full(ig_i, io, ispin) * coef_im_full(valid_basis_ids_spin(istate_frag, ispin), io, ispin)) * occ_factor
                  end do
                end do
              end do
              write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4,a,3(1pe12.4,1x))') "        density frag metric probe: rank=", dg_frag%id, &
                " ifrag=", ifrag, " ispin=", ispin, " Tr(SffDff)=", frag_trace_probe, &
                " state_metric=", frag_state_trace_probe(1), frag_state_trace_probe(2), frag_state_trace_probe(3)
              flush(6)
            end if
            call cpu_time(t_dmat1)
            time_project_dmat_build = time_project_dmat_build + (t_dmat1 - t_dmat0)
          end do
        else
          ! n_pw > 0: upfront subgroup assembly of coef_re/im on all icomm_frag ranks
          coef_re_full(1:nbf_max, 1:nocc_cache, 1:system%nspin) = 0.0d0
          coef_im_full(1:nbf_max, 1:nocc_cache, 1:system%nspin) = 0.0d0
          do ispin = 1, system%nspin
            nbf = dg_frag%n_basis(ifrag, ispin)
            nocc_spin = dg_frag%nocc_spin(ispin)
            if (nbf <= 0 .or. nocc_spin <= 0) cycle
            valid_basis_count = valid_basis_count_spin(ispin)
            coef_c_full(1:nbf_max, 1:nocc_spin) = (0.0d0, 0.0d0)
!$omp parallel do private(io, idx_local, istate_frag) schedule(static)
            do io = 1, nocc_spin
              do idx_local = 1, valid_basis_count
                istate_frag = valid_basis_ids_spin(idx_local, ispin)
                coef_c_full(istate_frag, io) = dg_frag%coef(basis_gid_spin(istate_frag, ispin), io, ispin)
              end do
            end do
!$omp end parallel do
            if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
              coef_c_frag(1:nbf_max, 1:nocc_spin) = (0.0d0, 0.0d0)
              call comm_summation(coef_c_full(1:nbf, 1:nocc_spin), coef_c_frag(1:nbf, 1:nocc_spin), &
                                  nbf * nocc_spin, dg_frag%icomm_frag)
              coef_c_full(1:nbf_max, 1:nocc_spin) = coef_c_frag(1:nbf_max, 1:nocc_spin)
            end if
            coef_re_full(1:nbf_max, 1:nocc_spin, ispin) = real(coef_c_full(1:nbf_max, 1:nocc_spin), kind=8)
            coef_im_full(1:nbf_max, 1:nocc_spin, ispin) = aimag(coef_c_full(1:nbf_max, 1:nocc_spin))
          end do
          ! NOTE: coef_pw_full_cache is already populated on all ranks by refresh_pw_coef_cache
          ! (called before the fragment loop), so no bcast is needed here.

          ! Compute state range for this rank
          nocc_per_rank_loc = (nocc_cache + dg_frag%isize_frag - 1) / dg_frag%isize_frag
          io_s_frag = dg_frag%id_frag * nocc_per_rank_loc + 1
          io_e_frag = min((dg_frag%id_frag + 1) * nocc_per_rank_loc, nocc_cache)
        end if
        ! --- end D pre-pass ---
        do block_offset = first_block_offset, nblocks_ifrag - 1, block_step_blocks
          igrid0 = 1 + block_offset * grid_block_size
          npt_blk = min(grid_block_size, ngrid - igrid0 + 1)
          if (npt_blk <= 0 .or. ngrid <= 0) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') "[FATAL] density block size invalid: rank=", dg_frag%id, &
              " ifrag=", ifrag, " i_local=", i_local, " block_offset=", block_offset, " ngrid=", ngrid, " npt_blk=", npt_blk
            flush(6)
            stop "DG-Fragment RT: density block size invalid"
          end if
          local_grid_count = 0
          remote_grid_count = 0
          valid_remote_grid_count = 0
          call cpu_time(t_setup0)
          if (allocated(dg_frag%density_send_slot_map)) then
            call prepare_grid_buffers_owner_map(i_local, igrid0, npt_blk, nxyz, n_pw == 0 .and. .not. dg_frag%parallel_mode_orbital)
          else
            slot_buf(1:npt_blk) = 0
            call prepare_grid_buffers_owner_map_no_slot(i_local, igrid0, npt_blk, nxyz)
          end if
          if (enable_density_reconstruct_trace .and. block_offset == first_block_offset) then
            ixg = ixg_buf(1)
            iyg = iyg_buf(1)
            izg = izg_buf(1)
            bx = map_global_to_phi_box_coord_ham(ixg, lbound(dg_frag%phi_frag, 1), ubound(dg_frag%phi_frag, 1), dg_frag%lgnum_total(1))
            by = map_global_to_phi_box_coord_ham(iyg, lbound(dg_frag%phi_frag, 2), ubound(dg_frag%phi_frag, 2), dg_frag%lgnum_total(2))
            bz = map_global_to_phi_box_coord_ham(izg, lbound(dg_frag%phi_frag, 3), ubound(dg_frag%phi_frag, 3), dg_frag%lgnum_total(3))
            phi_sample_probe = 0.0d0
            if (bx /= 0 .and. by /= 0 .and. bz /= 0) then
              phi_sample_probe = dg_frag%phi_frag(bx, by, bz, 1, i_local)
            end if
            write(*,'(1x,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,1pe12.4)') "        density map probe: rank=", dg_frag%id, &
              " ifrag=", ifrag, " global=", ixg, iyg, izg, " local=", bx, by, bz, " phi1=", phi_sample_probe
            flush(6)
            ixg_min_probe = huge(1)
            ixg_max_probe = -huge(1)
            owner_valid_probe = 0
            do igrid = 1, npt_blk
              if (owner_buf(igrid) < 0) cycle
              owner_valid_probe = owner_valid_probe + 1
              ixg_min_probe = min(ixg_min_probe, ixg_buf(igrid))
              ixg_max_probe = max(ixg_max_probe, ixg_buf(igrid))
            end do
            if (owner_valid_probe <= 0) then
              ixg_min_probe = -1
              ixg_max_probe = -1
            end if
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') "        density block probe: rank=", dg_frag%id, &
              " ifrag=", ifrag, " npt_blk=", npt_blk, " owner_valid=", owner_valid_probe, &
              " ixg_min=", ixg_min_probe, " ixg_max=", ixg_max_probe
            flush(6)
          end if
          do igrid = 1, npt_blk
            if (owner_buf(igrid) < 0) cycle
            local_grid_count = local_grid_count + 1
            if (local_grid_count > grid_block_size) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] density local_grid_count overflow: rank=", dg_frag%id, &
                " ifrag=", ifrag, " local_grid_count=", local_grid_count, " npt_blk=", npt_blk
              flush(6)
              stop "DG-Fragment RT: density local_grid_count overflow"
            end if
            local_grid_ids(local_grid_count) = igrid
            if (slot_buf(igrid) > 0) then
              remote_grid_count = remote_grid_count + 1
              if (remote_grid_count > grid_block_size) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] density remote_grid_count overflow: rank=", dg_frag%id, &
                  " ifrag=", ifrag, " remote_grid_count=", remote_grid_count, " npt_blk=", npt_blk
                flush(6)
                stop "DG-Fragment RT: density remote_grid_count overflow"
              end if
              remote_grid_ids(remote_grid_count) = igrid
              valid_remote_grid_count = valid_remote_grid_count + 1
              if (valid_remote_grid_count > grid_block_size) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] density valid_remote_grid_count overflow: rank=", dg_frag%id, &
                  " ifrag=", ifrag, " valid_remote=", valid_remote_grid_count, " npt_blk=", npt_blk
                flush(6)
                stop "DG-Fragment RT: density valid_remote_grid_count overflow"
              end if
              valid_remote_grid_ids(valid_remote_grid_count) = igrid
            end if
          end do
          if (enable_ifrag_compare_trace .and. itt_tag >= 8 .and. itt_tag <= 10 .and. block_offset == first_block_offset) then
            owner_valid_count = 0
            slot0_count = 0
            slotp_count = 0
            owner_true_count = 0
            owner_false_count = 0
            do igrid = 1, npt_blk
              if (owner_buf(igrid) >= 0) owner_valid_count = owner_valid_count + 1
              if (slot_buf(igrid) == 0) slot0_count = slot0_count + 1
              if (slot_buf(igrid) > 0) slotp_count = slotp_count + 1
              if (target_rank_owned_by_handler(owner_buf(igrid))) then
                owner_true_count = owner_true_count + 1
              else
                owner_false_count = owner_false_count + 1
              end if
            end do
            cmp_slot = 0
            if (ifrag == 1) cmp_slot = 1
            if (ifrag == 2) cmp_slot = 2
            if (cmp_slot > 0) then
              cmp_nxyz(:, cmp_slot) = nxyz(:)
              cmp_frag_lo(1, cmp_slot) = dg_frag%ixyz_frag(1, ifrag)
              cmp_frag_lo(2, cmp_slot) = dg_frag%ixyz_frag(2, ifrag)
              cmp_frag_lo(3, cmp_slot) = dg_frag%ixyz_frag(3, ifrag)
              cmp_frag_hi(1, cmp_slot) = dg_frag%ixyz_frag(1, ifrag) + nxyz(1) - 1
              cmp_frag_hi(2, cmp_slot) = dg_frag%ixyz_frag(2, ifrag) + nxyz(2) - 1
              cmp_frag_hi(3, cmp_slot) = dg_frag%ixyz_frag(3, ifrag) + nxyz(3) - 1
              cmp_npt(cmp_slot) = npt_blk
              cmp_local(cmp_slot) = local_grid_count
              cmp_remote(cmp_slot) = remote_grid_count
              cmp_valid_remote(cmp_slot) = valid_remote_grid_count
              cmp_owner_valid(cmp_slot) = owner_valid_count
              cmp_slot0(cmp_slot) = slot0_count
              cmp_slotp(cmp_slot) = slotp_count
              cmp_owner_true(cmp_slot) = owner_true_count
              cmp_owner_false(cmp_slot) = owner_false_count
              ifrag_grid_seen(cmp_slot) = .true.
            end if
            if (cmp_slot > 0 .and. dg_frag%is_frag_root) then
              write(*,*) '        ifrag-compare gridprep: itt=', itt_tag, ' ifrag=', ifrag, ' ispin=', 1, &
                ' nxyz=', cmp_nxyz(1,cmp_slot), cmp_nxyz(2,cmp_slot), cmp_nxyz(3,cmp_slot), &
                ' lo=', cmp_frag_lo(1,cmp_slot), cmp_frag_lo(2,cmp_slot), cmp_frag_lo(3,cmp_slot), &
                ' hi=', cmp_frag_hi(1,cmp_slot), cmp_frag_hi(2,cmp_slot), cmp_frag_hi(3,cmp_slot)
              write(*,*) '        ifrag-compare owner: npt=', cmp_npt(cmp_slot), ' local=', cmp_local(cmp_slot), &
                ' remote=', cmp_remote(cmp_slot), ' valid_remote=', cmp_valid_remote(cmp_slot), &
                ' ownerT=', cmp_owner_true(cmp_slot), ' ownerF=', cmp_owner_false(cmp_slot)
              flush(6)
            end if
          end if
          if (n_pw == 0 .and. enable_density_state_charge_trace) then
            do igrid = 1, npt_blk
              if (owner_buf(igrid) < 0) cycle
              if (.not. target_rank_owned_by_handler(owner_buf(igrid))) cycle
              if (slot_buf(igrid) > 0) then
                ifrag_import_point_count_local = ifrag_import_point_count_local + 1.0d0
              else
                ifrag_owned_point_count_local = ifrag_owned_point_count_local + 1.0d0
              end if
            end do
          end if
          call cpu_time(t_setup1)
          time_project_grid_prep = time_project_grid_prep + (t_setup1 - t_setup0)

          if (n_pw > 0) then
!$omp parallel do private(ixg, iyg, izg, ipw, theta) schedule(static)
            do igrid = 1, npt_blk
              ixg = ixg_buf(igrid)
              iyg = iyg_buf(igrid)
              izg = izg_buf(igrid)
!$omp simd
              do ipw = 1, n_pw
                theta = kpw_hx(ipw) * real(ixg, 8) + kpw_hy(ipw) * real(iyg, 8) + kpw_hz(ipw) * real(izg, 8)
                phase_cache(igrid, ipw) = cmplx(cos(theta), sin(theta), kind=8) * inv_sqrt_vol
              end do
            end do
!$omp end parallel do
            phase_applied_total_block = 0.0d0
            if (enable_ifrag_compare_trace .and. itt_tag >= 8 .and. itt_tag <= 10 .and. block_offset == first_block_offset) then
              phase_applied_total_block = sum(abs(phase_cache(1:npt_blk, 1:n_pw)))
            end if
          else
            phase_applied_total_block = 0.0d0
          end if

          do ispin = 1, system%nspin
            nocc_spin = dg_frag%nocc_spin(ispin)
            nbf = dg_frag%n_basis(ifrag, ispin)
            if (nbf <= 0 .or. nocc_spin <= 0) cycle
            valid_basis_count = valid_basis_count_spin(ispin)
            send_pre_block_raw = 0.0d0
            if (enable_ifrag_compare_trace .and. itt_tag >= 8 .and. itt_tag <= 10 .and. ispin == 1 .and. block_offset == first_block_offset) then
              basis_gid_probe(:) = 0
              do io = 1, min(3, nbf)
                basis_gid_probe(io) = dg_frag%index_basis(io, ifrag, ispin)
              end do
              cmp_slot = 0
              if (ifrag == 1) cmp_slot = 1
              if (ifrag == 2) cmp_slot = 2
              if (cmp_slot > 0) then
                cmp_nbf(cmp_slot) = nbf
                cmp_valid(cmp_slot) = valid_basis_count
                cmp_basis_gid(:, cmp_slot) = basis_gid_probe(:)
                cmp_phase_total(cmp_slot) = phase_applied_total_block
                ifrag_basis_seen(cmp_slot) = .true.
              end if
              if (cmp_slot > 0 .and. dg_frag%is_frag_root) then
                write(*,*) '        ifrag-compare basis: itt=', itt_tag, ' ifrag=', ifrag, ' ispin=', ispin, ' block=', block_offset, &
                  ' nbf=', cmp_nbf(cmp_slot), ' valid=', cmp_valid(cmp_slot), ' gid=', cmp_basis_gid(1,cmp_slot), &
                  cmp_basis_gid(2,cmp_slot), cmp_basis_gid(3,cmp_slot), ' phase=', cmp_phase_total(cmp_slot)
                flush(6)
              end if
            end if

          call cpu_time(t_setup0)
          if (enable_density_phi_block_cache) then
            phi_blk(1:npt_blk, 1:nbf) = dg_frag%density_phi_block_cache(1:npt_blk, 1:nbf, block_offset + 1, i_local)
          else
            phi_blk(1:npt_blk, 1:nbf) = 0.0d0
!$omp parallel do private(igrid, ixg, iyg, izg, bx, by, bz, istate_frag) schedule(static)
            do igrid = 1, npt_blk
              ixg = ixg_buf(igrid)
              iyg = iyg_buf(igrid)
              izg = izg_buf(igrid)
              bx = map_global_to_phi_box_coord_ham(ixg, phi_lb1, phi_ub1, phi_lg1)
              by = map_global_to_phi_box_coord_ham(iyg, phi_lb2, phi_ub2, phi_lg2)
              bz = map_global_to_phi_box_coord_ham(izg, phi_lb3, phi_ub3, phi_lg3)
              if (bx == 0 .or. by == 0 .or. bz == 0) cycle
              do istate_frag = 1, nbf
                phi_blk(igrid, istate_frag) = dg_frag%phi_frag(bx, by, bz, istate_frag, i_local)
              end do
            end do
!$omp end parallel do
          end if
          call cpu_time(t_setup1)
          time_project_phi_pack = time_project_phi_pack + (t_setup1 - t_setup0)

            if (use_mixed_density) then
              n_basis_mix = n_basis_mix_spin(ispin)
              if (n_basis_mix <= 0) cycle
              if (enable_density_reconstruct_trace .and. block_offset == first_block_offset .and. ispin == 1) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') "        density mixed trace: rank=", dg_frag%id, &
                  " ifrag=", ifrag, " block_offset=", block_offset, " npt_blk=", npt_blk, " n_basis_mix=", n_basis_mix
                flush(6)
              end if
              call cpu_time(t_psi0)
              if (enable_density_reconstruct_trace .and. block_offset == first_block_offset .and. ispin == 1) then
                write(*,'(1x,a,i0,a)') "        density mixed trace: rank=", dg_frag%id, " stage=before-frag-transform"
                flush(6)
              end if
              basis_mix_blk(1:npt_blk, 1:n_basis_mix) = matmul(phi_blk(1:npt_blk, 1:nbf), &
                transform_frag_spin(1:nbf, 1:n_basis_mix, ispin))
              if (n_pw > 0) then
                if (enable_density_reconstruct_trace .and. block_offset == first_block_offset .and. ispin == 1) then
                  write(*,'(1x,a,i0,a)') "        density mixed trace: rank=", dg_frag%id, " stage=before-zgemm-pw"
                  flush(6)
                end if
                call zgemm('N', 'N', npt_blk, n_basis_mix, n_pw, zone, phase_cache, grid_block_size, &
                  transform_pw_spin(1, 1, ispin), n_pw, zone, basis_mix_blk, grid_block_size)
              end if
              if (enable_density_reconstruct_trace .and. block_offset == first_block_offset .and. ispin == 1) then
                write(*,'(1x,a,i0,a)') "        density mixed trace: rank=", dg_frag%id, " stage=before-zgemm-density"
                flush(6)
              end if
              call zgemm('N', 'N', npt_blk, n_basis_mix, n_basis_mix, zone, basis_mix_blk, grid_block_size, &
                density_mix(1, 1, ispin), max_mixed_basis, zzero, density_mix_tmp, grid_block_size)
              if (n_pw > 0) then
                basis_mix_blk_t(1:n_basis_mix, 1:npt_blk) = transpose(basis_mix_blk(1:npt_blk, 1:n_basis_mix))
                density_mix_tmp_t(1:n_basis_mix, 1:npt_blk) = transpose(density_mix_tmp(1:npt_blk, 1:n_basis_mix))
              end if
              if (enable_density_reconstruct_trace .and. block_offset == first_block_offset .and. ispin == 1) then
                write(*,'(1x,a,i0,a)') "        density mixed trace: rank=", dg_frag%id, " stage=before-rho-accum-loop"
                flush(6)
              end if
              call cpu_time(t_psi1)
              time_project_psi = time_project_psi + (t_psi1 - t_psi0)
              call cpu_time(t_rho0)

!$omp parallel do private(io, io0, io1, rho_mix_accum) schedule(static)
              do igrid = 1, npt_blk
                rho_mix_accum = 0.0d0
                if (n_pw > 0) then
                  do io0 = 1, n_basis_mix, mixed_io_block_size
                    io1 = min(n_basis_mix, io0 + mixed_io_block_size - 1)
!$omp simd reduction(+:rho_mix_accum)
                    do io = io0, io1
                      rho_mix_accum = rho_mix_accum + real(conjg(basis_mix_blk_t(io, igrid)) * density_mix_tmp_t(io, igrid), kind=8)
                    end do
                  end do
                else
!$omp simd reduction(+:rho_mix_accum)
                  do io = 1, n_basis_mix
                    rho_mix_accum = rho_mix_accum + real(conjg(basis_mix_blk(igrid, io)) * density_mix_tmp(igrid, io), kind=8)
                  end do
                end if
                rho_blk(igrid) = rho_mix_accum
              end do
!$omp end parallel do
              if (enable_density_reconstruct_trace .and. block_offset == first_block_offset .and. ispin == 1) then
                write(*,'(1x,a,i0,a)') "        density mixed trace: rank=", dg_frag%id, " stage=after-rho-accum-loop"
                flush(6)
              end if
              if (enable_rho_mix_grid_compare .and. block_offset == first_block_offset .and. dg_frag%id == 0) then
                rho_blk_reduced(1:npt_blk) = 0.0d0
                do io0 = 1, nocc_spin, state_block_size
                  nbatch = min(state_block_size, nocc_spin - io0 + 1)
                  call zgemm('N', 'N', npt_blk, nbatch, n_basis_mix, zone, basis_mix_blk, grid_block_size, &
                    coef_mix_spin(1, io0, ispin), max_mixed_basis, zzero, pw_tmp_z, grid_block_size)
                  do io = 1, nbatch
                    occ_factor = 1.0d0
                    if (allocated(system%rocc)) then
                      if (io0 + io - 1 <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
                        occ_factor = max(0.0d0, system%rocc(io0 + io - 1, 1, ispin))
                      end if
                    end if
                    if (occ_factor <= 0.0d0) cycle
!$omp parallel do private(igrid) schedule(static)
                    do igrid = 1, npt_blk
                      rho_blk_reduced(igrid) = rho_blk_reduced(igrid) + occ_factor * &
                        (real(pw_tmp_z(igrid, io), kind=8)**2 + aimag(pw_tmp_z(igrid, io))**2)
                    end do
!$omp end parallel do
                  end do
                end do
                rho_mix_grid_l2 = 0.0d0
                rho_mix_grid_ref = 0.0d0
                rho_mix_grid_max = 0.0d0
                do igrid = 1, npt_blk
                  rho_mix_grid_l2 = rho_mix_grid_l2 + (rho_blk(igrid) - rho_blk_reduced(igrid))**2
                  rho_mix_grid_ref = rho_mix_grid_ref + rho_blk_reduced(igrid)**2
                  rho_mix_grid_max = max(rho_mix_grid_max, abs(rho_blk(igrid) - rho_blk_reduced(igrid)))
                end do
                rho_mix_grid_l2 = sqrt(max(0.0d0, rho_mix_grid_l2))
                rho_mix_grid_ref = sqrt(max(0.0d0, rho_mix_grid_ref))
                rho_mix_grid_rel = 0.0d0
                if (rho_mix_grid_ref > 1.0d-20) rho_mix_grid_rel = rho_mix_grid_l2 / rho_mix_grid_ref
                write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4)') '        rho_mix grid2path: ifrag=', ifrag, &
                  ' ispin=', ispin, ' npt=', npt_blk, ' l2=', rho_mix_grid_l2, ' rel=', rho_mix_grid_rel, ' max=', rho_mix_grid_max
                flush(6)
              end if
              if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
                call cpu_time(t_rho1)
                call comm_summation(rho_blk(1:npt_blk), rho_blk_reduced(1:npt_blk), npt_blk, dg_frag%icomm_frag)
                rho_blk(1:npt_blk) = rho_blk_reduced(1:npt_blk)
                call cpu_time(t_setup1)
                time_project_rho_reduce = time_project_rho_reduce + (t_setup1 - t_rho1)
              end if


              ! orbital mode: all orbital ranks share the same real-space domain and
              ! have identical rho_blk after icomm_frag reduce; all must materialize.
              ! Remote send (outside-fragment owners) is still frag_root-only.
              if (dg_frag%is_frag_root .or. dg_frag%parallel_mode_orbital) then
                do igrid = 1, npt_blk
                  ixg = ixg_buf(igrid)
                  iyg = iyg_buf(igrid)
                  izg = izg_buf(igrid)
                  if (slot_buf(igrid) > 0) cycle
                    if (.not. target_rank_owned_by_handler(owner_buf(igrid))) cycle
                    if (ixg < rho_s_x_lo .or. ixg > rho_s_x_hi .or. &
                      iyg < rho_s_y_lo .or. iyg > rho_s_y_hi .or. &
                      izg < rho_s_z_lo .or. izg > rho_s_z_hi) cycle
                  rho_raw_contrib = rho_blk(igrid)
                  rho_contrib = rho_raw_contrib * merge(owned_path_scale, 1.0d0, enable_owned_path_normalized)
                  bx = map_global_to_phi_box_coord_ham(ixg, lbound(rho_bf, 1), ubound(rho_bf, 1), dg_frag%lgnum_total(1))
                  by = map_global_to_phi_box_coord_ham(iyg, lbound(rho_bf, 2), ubound(rho_bf, 2), dg_frag%lgnum_total(2))
                  bz = map_global_to_phi_box_coord_ham(izg, lbound(rho_bf, 3), ubound(rho_bf, 3), dg_frag%lgnum_total(3))
                  if (bx == 0 .or. by == 0 .or. bz == 0) then
                    write(*,'(1x,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
                      "[FATAL] density local rho_bf map failed: rank=", dg_frag%id, &
                      " id_frag=", dg_frag%id_frag, " idx=", ixg, iyg, izg, &
                      " lb=", lbound(rho_bf, 1), lbound(rho_bf, 2), lbound(rho_bf, 3), &
                      " ub=", ubound(rho_bf, 1), ubound(rho_bf, 2), ubound(rho_bf, 3)
                    flush(6)
                    stop "DG-Fragment RT: density local rho_bf map failed"
                  end if
                  rho_bf(bx, by, bz) = rho_bf(bx, by, bz) + rho_contrib
                  rho_s_bf(ixg, iyg, izg, ispin) = rho_s_bf(ixg, iyg, izg, ispin) + rho_contrib
                end do
              end if
              if (dg_frag%is_frag_root) then
                do idx_remote = 1, valid_remote_grid_count
                  igrid = valid_remote_grid_ids(idx_remote)
                  owner_rank = owner_buf(igrid)
                  slot = slot_buf(igrid)
                  if (owner_rank < 0 .or. owner_rank > dg_frag%isize - 1) then
                    write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                      "DG density remote owner out of range: rank=", dg_frag%id, &
                      " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                      " igrid=", igrid, " owner=", owner_rank, " valid_remote=", valid_remote_grid_count
                    flush(6)
                    stop "DG-Fragment RT: density remote owner out of range"
                  end if
                  if (.not. allocated(rho_send(owner_rank)%f)) then
                    write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                      "DG density remote send buffer missing: rank=", dg_frag%id, &
                      " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                      " igrid=", igrid, " owner=", owner_rank, " nsend=", dg_frag%density_send_count(owner_rank)
                    flush(6)
                    stop "DG-Fragment RT: density remote send buffer missing"
                  end if
                  if (slot < 1 .or. slot > dg_frag%density_send_count(owner_rank)) then
                    write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                      "DG density remote slot out of range: rank=", dg_frag%id, &
                      " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                      " igrid=", igrid, " owner=", owner_rank, " slot=", slot, &
                      " nsend=", dg_frag%density_send_count(owner_rank)
                    flush(6)
                    stop "DG-Fragment RT: density remote slot out of range"
                  end if
                  rho_contrib = rho_blk(igrid)
                  send_pre_block_raw = send_pre_block_raw + rho_contrib * system%hvol
                  rho_send(owner_rank)%f(slot, 1, 1) = rho_send(owner_rank)%f(slot, 1, 1) + rho_contrib
                  spin_offset = ispin * dg_frag%density_send_count(owner_rank)
                  if (spin_offset + slot < 1 .or. spin_offset + slot > size(rho_send(owner_rank)%f, 1)) then
                    write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                      "DG density remote spin slot out of range: rank=", dg_frag%id, &
                      " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                      " igrid=", igrid, " owner=", owner_rank, " slot=", slot, &
                      " spin_slot=", spin_offset + slot, " send_size=", size(rho_send(owner_rank)%f, 1)
                    flush(6)
                    stop "DG-Fragment RT: density remote spin slot out of range"
                  end if
                  rho_send(owner_rank)%f(spin_offset + slot, 1, 1) = rho_send(owner_rank)%f(spin_offset + slot, 1, 1) + rho_contrib
                end do
              end if
              call cpu_time(t_rho1)
              time_project_rho = time_project_rho + (t_rho1 - t_rho0)
            else
              ! D already computed in pre-pass for n_pw == 0
              rho_blk_accum(1:npt_blk) = 0.0d0
              if (n_pw == 0) then
                if (.not. allocated(rho_blk_partial)) allocate(rho_blk_partial(grid_block_size))
                rho_blk_partial(1:npt_blk) = 0.0d0
                if (enable_density_reconstruct_trace) then
                  nprobe_cols = min(3, nbf)
                  do iprobe = 1, nprobe_cols
                    do igrid = 1, npt_blk
                      if (.not. target_rank_owned_by_handler(owner_buf(igrid))) cycle
                      if (slot_buf(igrid) > 0) cycle
                      phi_col_metric_total(iprobe, ispin, i_local) = phi_col_metric_total(iprobe, ispin, i_local) + &
                        phi_blk(igrid, iprobe) * phi_blk(igrid, iprobe) * system%hvol
                    end do
                    do io = 1, nprobe_cols
                      do igrid = 1, npt_blk
                        if (.not. target_rank_owned_by_handler(owner_buf(igrid))) cycle
                        if (slot_buf(igrid) > 0) cycle
                        phi_gram_total(iprobe, io, ispin, i_local) = phi_gram_total(iprobe, io, ispin, i_local) + &
                          phi_blk(igrid, iprobe) * phi_blk(igrid, io) * system%hvol
                      end do
                    end do
                  end do
                  do iprobe = 1, nbf
                    do io = 1, nbf
                      do igrid = 1, npt_blk
                        if (.not. target_rank_owned_by_handler(owner_buf(igrid))) cycle
                        if (slot_buf(igrid) > 0) cycle
                        phi_frag_metric_total(iprobe, io, ispin, i_local) = phi_frag_metric_total(iprobe, io, ispin, i_local) + &
                          phi_blk(igrid, iprobe) * phi_blk(igrid, io) * system%hvol
                      end do
                    end do
                  end do
                end if
                ! Distribute occupied states across icomm_frag ranks to avoid
                ! duplicate accumulation before subgroup reduction.
                nocc_per_rank_loc = (nocc_spin + dg_frag%isize_frag - 1) / dg_frag%isize_frag
                io_s_frag = dg_frag%id_frag * nocc_per_rank_loc + 1
                io_e_frag = min((dg_frag%id_frag + 1) * nocc_per_rank_loc, nocc_spin)
                do io0 = io_s_frag, io_e_frag, state_block_size
                  nbatch = min(state_block_size, io_e_frag - io0 + 1)
                  call cpu_time(t_psi0)
                  coef_blk_ri(1:nbf, 1:nbatch) = coef_re_full(1:nbf, io0:io0+nbatch-1, ispin)
                  coef_blk_ri(1:nbf, nbatch+1:2*nbatch) = coef_im_full(1:nbf, io0:io0+nbatch-1, ispin)
                  call dgemm('N', 'N', npt_blk, 2*nbatch, nbf, 1.0d0, phi_blk, grid_block_size, &
                             coef_blk_ri, nbf_max, 0.0d0, psi_blk_ri, grid_block_size)
                  psi_blk_re(1:npt_blk, 1:nbatch) = psi_blk_ri(1:npt_blk, 1:nbatch)
                  psi_blk_im(1:npt_blk, 1:nbatch) = psi_blk_ri(1:npt_blk, nbatch+1:2*nbatch)
                  if (block_offset == first_block_offset .and. io0 == io_s_frag) then
                    phi_norm_probe = 0.0d0
                    psi_norm_probe = 0.0d0
                    phi_col_norm_probe(:) = 0.0d0
                    nprobe_cols = min(3, nbf)
                    do idx_local = 1, local_grid_count
                      igrid = local_grid_ids(idx_local)
                      phi_norm_probe = phi_norm_probe + sum(phi_blk(igrid, 1:nbf) * phi_blk(igrid, 1:nbf))
                      psi_norm_probe = psi_norm_probe + sum(psi_blk_re(igrid, 1:nbatch) * psi_blk_re(igrid, 1:nbatch) + &
                                                           psi_blk_im(igrid, 1:nbatch) * psi_blk_im(igrid, 1:nbatch))
                      do iprobe = 1, nprobe_cols
                        phi_col_norm_probe(iprobe) = phi_col_norm_probe(iprobe) + phi_blk(igrid, iprobe) * phi_blk(igrid, iprobe)
                      end do
                    end do
                    if (enable_density_reconstruct_trace) then
                      write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4,a,1pe12.4)') "        density psi probe: rank=", dg_frag%id, &
                        " ifrag=", ifrag, " ispin=", ispin, " phi_norm=", phi_norm_probe, " psi_norm=", psi_norm_probe
                      flush(6)
                      write(*,'(1x,a,i0,a,i0,a,3(1pe12.4,1x))') "        density phi-col probe: rank=", dg_frag%id, &
                        " ifrag=", ifrag, " cols=", phi_col_norm_probe(1), phi_col_norm_probe(2), phi_col_norm_probe(3)
                      flush(6)
                    end if
                  end if
                  call cpu_time(t_psi1)
                  time_project_psi = time_project_psi + (t_psi1 - t_psi0)
                  call cpu_time(t_rho0)
!$omp parallel do private(io, igrid, rho_accum, occ_factor) schedule(static)
                    do igrid = 1, npt_blk
                      rho_accum = 0.0d0
!$omp simd reduction(+:rho_accum)
                      do io = 1, nbatch
                        occ_factor = occ_cache(io0 + io - 1)
                        rho_accum = rho_accum + occ_factor * &
                          (psi_blk_re(igrid, io) * psi_blk_re(igrid, io) + psi_blk_im(igrid, io) * psi_blk_im(igrid, io))
                      end do
                      rho_blk_partial(igrid) = rho_blk_partial(igrid) + rho_accum
                    end do
!$omp end parallel do
                  if (enable_density_state_charge_trace .or. enable_state_rhobf_trace) then
                    ! Each rank accumulates its state partition over owned grid points
                    ! (slot==0 && target_rank_owned_by_handler).  All ranks participate;
                    ! after global allreduce the per-state sum equals the total owned_local
                    ! charge.  To compare with final_integrated, add imported_halo separately.
                    do io = 1, nbatch
                      occ_factor = occ_cache(io0 + io - 1)
                      if (occ_factor <= 0.0d0) cycle
                      do igrid = 1, npt_blk
                        psi2_val = psi_blk_re(igrid, io) * psi_blk_re(igrid, io) + &
                                   psi_blk_im(igrid, io) * psi_blk_im(igrid, io)
                        state_psi2_raw_local(io0 + io - 1, ispin) = state_psi2_raw_local(io0 + io - 1, ispin) + psi2_val
                        state_psi2_occ_local(io0 + io - 1, ispin) = state_psi2_occ_local(io0 + io - 1, ispin) + occ_factor * psi2_val
                        state_psi2_dv_local(io0 + io - 1, ispin) = state_psi2_dv_local(io0 + io - 1, ispin) + occ_factor * psi2_val * system%hvol
                        if (has_density_point_probe) then
                          if (ixg_buf(igrid) == probe_ixg .and. iyg_buf(igrid) == probe_iyg .and. izg_buf(igrid) == probe_izg) then
                            if (slot_buf(igrid) > 0) then
                              if (target_rank_owned_by_handler(owner_buf(igrid))) then
                                probe_state_import_local(io0 + io - 1, ispin) = probe_state_import_local(io0 + io - 1, ispin) + &
                                  merge(imported_unpack_scale, 1.0d0, enable_imported_unpack_normalized) * occ_factor * psi2_val * system%hvol
                              end if
                            else
                              if (target_rank_owned_by_handler(owner_buf(igrid))) then
                                probe_state_owned_local(io0 + io - 1, ispin) = probe_state_owned_local(io0 + io - 1, ispin) + &
                                  merge(owned_path_scale, 1.0d0, enable_owned_path_normalized) * occ_factor * psi2_val * system%hvol
                              end if
                            end if
                          end if
                        end if
                        if (slot_buf(igrid) > 0) then
                          if (target_rank_owned_by_handler(owner_buf(igrid))) then
                            state_import_core_local(io0 + io - 1, ispin) = state_import_core_local(io0 + io - 1, ispin) + &
                              merge(imported_unpack_scale, 1.0d0, enable_imported_unpack_normalized) * occ_factor * psi2_val * system%hvol
                          end if
                          cycle
                        end if
                        if (.not. target_rank_owned_by_handler(owner_buf(igrid))) cycle
                        state_psi2_owned_local(io0 + io - 1, ispin) = state_psi2_owned_local(io0 + io - 1, ispin) + &
                          merge(owned_path_scale, 1.0d0, enable_owned_path_normalized) * occ_factor * psi2_val * system%hvol
                        state_charge_local(io0 + io - 1, ispin) = state_charge_local(io0 + io - 1, ispin) + &
                          merge(owned_path_scale, 1.0d0, enable_owned_path_normalized) * occ_factor * psi2_val * system%hvol
                        if (enable_state_rhobf_trace) then
                          if (ispin == state_rhobf_trace_spin .and. io0 + io - 1 == state_rhobf_trace_io) then
                            state_rhobf_psi2_q_local = state_rhobf_psi2_q_local + psi2_val
                            psi2_occ_val = occ_factor * psi2_val
                            psi2_dv_val = psi2_occ_val * system%hvol
                            state_rhobf_psi2_occ_q_local = state_rhobf_psi2_occ_q_local + psi2_occ_val
                            state_rhobf_psi2_dv_q_local = state_rhobf_psi2_dv_q_local + psi2_dv_val
                            if (target_rank_owned_by_handler(owner_buf(igrid))) then
                              if (slot_buf(igrid) > 0) then
                                psi2_after_partition_val = psi2_dv_val
                                psi2_after_slot_val = psi2_after_partition_val
                                psi2_after_any_norm_val = psi2_after_slot_val * &
                                  merge(imported_unpack_scale, 1.0d0, enable_imported_unpack_normalized)
                              else
                                state_rhobf_psi2_owned_q_local = state_rhobf_psi2_owned_q_local + psi2_dv_val
                                psi2_after_partition_val = psi2_dv_val * &
                                  merge(owned_path_scale, 1.0d0, enable_owned_path_normalized)
                                psi2_after_slot_val = psi2_after_partition_val
                                psi2_after_any_norm_val = psi2_after_slot_val
                              end if
                              state_rhobf_psi2_after_partition_q_local = state_rhobf_psi2_after_partition_q_local + psi2_after_partition_val
                              state_rhobf_psi2_after_slot_q_local = state_rhobf_psi2_after_slot_q_local + psi2_after_slot_val
                              state_rhobf_psi2_after_any_norm_q_local = state_rhobf_psi2_after_any_norm_q_local + psi2_after_any_norm_val
                            end if
                          end if
                        end if
                      end do
                    end do
                  end if
                  if (enable_density_reconstruct_trace .and. io0 <= size(orbital_norm_probe_local)) then
                    do io = 1, min(nbatch, size(orbital_norm_probe_local) - io0 + 1)
                      do idx_local = 1, local_grid_count
                        igrid = local_grid_ids(idx_local)
                        orbital_norm_probe_local(io0 + io - 1) = orbital_norm_probe_local(io0 + io - 1) + &
                          (psi_blk_re(igrid, io) * psi_blk_re(igrid, io) + &
                           psi_blk_im(igrid, io) * psi_blk_im(igrid, io)) * system%hvol
                        if (ifrag <= size(orbital_norm_frag_local, 1)) then
                          orbital_norm_frag_local(ifrag, io0 + io - 1) = orbital_norm_frag_local(ifrag, io0 + io - 1) + &
                            (psi_blk_re(igrid, io) * psi_blk_re(igrid, io) + &
                             psi_blk_im(igrid, io) * psi_blk_im(igrid, io)) * system%hvol
                        end if
                      end do
                    end do
                  end if
                  call cpu_time(t_rho1)
                  time_project_rho = time_project_rho + (t_rho1 - t_rho0)
                end do
                if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
                  call cpu_time(t_setup0)
                  if (itt_tag == 1) then
                    write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] BEFORE_COMM_SUM_DENSITY_RHOBLK rank=', &
                      dg_frag%id, ' itt=', itt_tag, ' ifrag=', ifrag, ' ispin=', ispin, ' coef_ref_ready=', dg_frag%coef_ref_ready
                    flush(6)
                  end if
                  call comm_summation(rho_blk_partial(1:npt_blk), rho_blk_reduced(1:npt_blk), npt_blk, dg_frag%icomm_frag)
                  if (itt_tag == 1) then
                    write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] AFTER_COMM_SUM_DENSITY_RHOBLK rank=', &
                      dg_frag%id, ' itt=', itt_tag, ' ifrag=', ifrag, ' ispin=', ispin, ' coef_ref_ready=', dg_frag%coef_ref_ready
                    flush(6)
                  end if
                  rho_blk_accum(1:npt_blk) = rho_blk_reduced(1:npt_blk)
                  call cpu_time(t_setup1)
                  time_project_rho_reduce = time_project_rho_reduce + (t_setup1 - t_setup0)
                else
                  rho_blk_accum(1:npt_blk) = rho_blk_partial(1:npt_blk)
                end if
                if (enable_density_reconstruct_trace .and. block_offset == first_block_offset) then
                  rho_probe_charge = sum(rho_blk_accum(1:npt_blk)) * system%hvol
                  write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4)') "        density frag probe: rank=", dg_frag%id, &
                    " ifrag=", ifrag, " ispin=", ispin, " first_block_charge=", rho_probe_charge
                  flush(6)
                end if
              else
                ! n_pw > 0: state-distributed loop, no per-batch bcast
                if (.not. allocated(rho_blk_partial)) allocate(rho_blk_partial(grid_block_size))
                rho_blk_partial(1:npt_blk) = 0.0d0
                occ_cache(1:nocc_cache) = 0.0d0
                occ_cache(1:nocc_spin) = 1.0d0
                if (allocated(system%rocc)) then
                  do io = 1, nocc_spin
                    if (io <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
                      occ_cache(io) = max(0.0d0, system%rocc(io, 1, ispin))
                    end if
                  end do
                end if

                do io0 = io_s_frag, min(io_e_frag, nocc_spin), state_block_size
                  nbatch = min(state_block_size, min(io_e_frag, nocc_spin) - io0 + 1)

                  ! copy coef from upfront buffer (no bcast)
                  coef_blk_re(1:nbf, 1:nbatch) = coef_re_full(1:nbf, io0:io0+nbatch-1, ispin)
                  coef_blk_im(1:nbf, 1:nbatch) = coef_im_full(1:nbf, io0:io0+nbatch-1, ispin)

                  call cpu_time(t_psi0)
                  call dgemm('N', 'N', npt_blk, nbatch, nbf, 1.0d0, phi_blk, grid_block_size, &
                             coef_blk_re, nbf_max, 0.0d0, psi_blk_re, grid_block_size)
                  call dgemm('N', 'N', npt_blk, nbatch, nbf, 1.0d0, phi_blk, grid_block_size, &
                             coef_blk_im, nbf_max, 0.0d0, psi_blk_im, grid_block_size)
                  do ipw0 = 1, n_pw, pw_block_size
                    npw_blk = min(pw_block_size, n_pw - ipw0 + 1)
                    ! direct access from full_cache on all ranks (no bcast)
                    coef_pw_blk(1:npw_blk, 1:nbatch) = &
                      dg_frag%coef_pw_full_cache(ipw0:ipw0+npw_blk-1, io0:io0+nbatch-1, ispin)

                    call zgemm('N', 'N', npt_blk, nbatch, npw_blk, zone, phase_cache(1, ipw0), grid_block_size, &
                               coef_pw_blk, pw_block_size, zzero, pw_tmp_z, grid_block_size)
                    psi_blk_re(1:npt_blk, 1:nbatch) = psi_blk_re(1:npt_blk, 1:nbatch) + real(pw_tmp_z(1:npt_blk, 1:nbatch), kind=8)
                    psi_blk_im(1:npt_blk, 1:nbatch) = psi_blk_im(1:npt_blk, 1:nbatch) + aimag(pw_tmp_z(1:npt_blk, 1:nbatch))
                  end do
                  call cpu_time(t_psi1)
                  time_project_psi = time_project_psi + (t_psi1 - t_psi0)

                  call cpu_time(t_rho0)
                  do io1 = 1, nbatch, rho_state_block_size
                    nstate = min(rho_state_block_size, nbatch - io1 + 1)
!$omp parallel do private(io, idx_local, igrid, rho_accum) schedule(static)
                    do idx_local = 1, local_grid_count
                      igrid = local_grid_ids(idx_local)
                      rho_accum = 0.0d0
!$omp simd reduction(+:rho_accum)
                      do io = io1, io1 + nstate - 1
                        rho_accum = rho_accum + occ_cache(io0 + io - 1) * &
                                    (psi_blk_re(igrid, io) * psi_blk_re(igrid, io) + &
                                    psi_blk_im(igrid, io) * psi_blk_im(igrid, io))
                      end do
                      rho_blk_partial(igrid) = rho_blk_partial(igrid) + rho_accum
                    end do
!$omp end parallel do
                  end do
                  call cpu_time(t_rho1)
                  time_project_rho = time_project_rho + (t_rho1 - t_rho0)
                end do  ! io0

                ! AllReduce rho_blk_partial across icomm_frag → rho_blk_accum
                ! Note: comm_summation overwrites rho_blk_accum (does not accumulate into it)
                call cpu_time(t_rho0)
                if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
                  if (itt_tag == 1) then
                    write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] BEFORE_COMM_SUM_DENSITY_RHOBLK_MIX rank=', &
                      dg_frag%id, ' itt=', itt_tag, ' ifrag=', ifrag, ' ispin=', ispin, ' coef_ref_ready=', dg_frag%coef_ref_ready
                    flush(6)
                  end if
                  call comm_summation(rho_blk_partial(1:npt_blk), rho_blk_accum(1:npt_blk), npt_blk, dg_frag%icomm_frag)
                  if (itt_tag == 1) then
                    write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] AFTER_COMM_SUM_DENSITY_RHOBLK_MIX rank=', &
                      dg_frag%id, ' itt=', itt_tag, ' ifrag=', ifrag, ' ispin=', ispin, ' coef_ref_ready=', dg_frag%coef_ref_ready
                    flush(6)
                  end if
                else
                  rho_blk_accum(1:npt_blk) = rho_blk_partial(1:npt_blk)
                end if
                call cpu_time(t_rho1)
                time_project_rho_reduce = time_project_rho_reduce + (t_rho1 - t_rho0)
              end if
              if (itt_tag == 1) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] AFTER_DENSITY_RHOBLK_PHASE rank=', &
                  dg_frag%id, ' itt=', itt_tag, ' i_local=', i_local, ' ifrag=', ifrag, ' ispin=', ispin, &
                  ' coef_ref_ready=', dg_frag%coef_ref_ready
                flush(6)
              end if
              ! rho_blk_accum: filled by dgemm-path (n_pw==0) or AllReduce (n_pw>0)
              call cpu_time(t_rho0)
!$omp parallel private(igrid, owner_rank, ixg, iyg, izg, bx, by, bz, rho_contrib, rho_raw_contrib, slot, theta)
!$omp do schedule(static)
                  do igrid = 1, npt_blk
                    ! orbital mode: rho_blk_accum is AllReduced across icomm_frag;
                    ! all orbital ranks hold identical data and must write to rho_bf.
                    if (.not. dg_frag%is_frag_root .and. .not. dg_frag%parallel_mode_orbital) cycle
                    ixg = ixg_buf(igrid)
                    iyg = iyg_buf(igrid)
                    izg = izg_buf(igrid)
                    if (slot_buf(igrid) > 0) cycle
                    if (.not. target_rank_owned_by_handler(owner_buf(igrid))) cycle
                    if (ixg < rho_s_x_lo .or. ixg > rho_s_x_hi .or. &
                      iyg < rho_s_y_lo .or. iyg > rho_s_y_hi .or. &
                      izg < rho_s_z_lo .or. izg > rho_s_z_hi) cycle
                    rho_raw_contrib = rho_blk_accum(igrid)
                    rho_contrib = rho_raw_contrib * merge(owned_path_scale, 1.0d0, enable_owned_path_normalized)
                    if (enable_density_weight_path_trace .and. has_density_point_probe) then
                      if (ixg == probe_ixg .and. iyg == probe_iyg .and. izg == probe_izg) then
                        probe_owned_pre_weight_local = probe_owned_pre_weight_local + rho_raw_contrib * system%hvol
                        probe_owned_add_local = probe_owned_add_local + rho_contrib * system%hvol
                        probe_owned_apply_count_local = probe_owned_apply_count_local + 1.0d0
                        probe_owned_weight_sum_local = probe_owned_weight_sum_local + merge(owned_path_scale, 1.0d0, enable_owned_path_normalized)
                      end if
                    end if
                    bx = map_global_to_phi_box_coord_ham(ixg, lbound(rho_bf, 1), ubound(rho_bf, 1), dg_frag%lgnum_total(1))
                    by = map_global_to_phi_box_coord_ham(iyg, lbound(rho_bf, 2), ubound(rho_bf, 2), dg_frag%lgnum_total(2))
                    bz = map_global_to_phi_box_coord_ham(izg, lbound(rho_bf, 3), ubound(rho_bf, 3), dg_frag%lgnum_total(3))
                    if (bx == 0 .or. by == 0 .or. bz == 0) then
                      write(*,'(1x,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
                        "[FATAL] density accum rho_bf map failed: rank=", dg_frag%id, &
                        " id_frag=", dg_frag%id_frag, " idx=", ixg, iyg, izg, &
                        " lb=", lbound(rho_bf, 1), lbound(rho_bf, 2), lbound(rho_bf, 3), &
                        " ub=", ubound(rho_bf, 1), ubound(rho_bf, 2), ubound(rho_bf, 3)
                      flush(6)
                      stop "DG-Fragment RT: density accum rho_bf map failed"
                    end if
                    rho_bf(bx, by, bz) = rho_bf(bx, by, bz) + rho_contrib
                    rho_s_bf(ixg, iyg, izg, ispin) = rho_s_bf(ixg, iyg, izg, ispin) + rho_contrib
                  end do
!$omp end do nowait
!$omp do schedule(static)
                  do idx_remote = 1, valid_remote_grid_count
                    if (.not. dg_frag%is_frag_root) cycle
                    igrid = valid_remote_grid_ids(idx_remote)
                    owner_rank = owner_buf(igrid)
                    slot = slot_buf(igrid)
                    if (owner_rank < 0 .or. owner_rank > dg_frag%isize - 1) then
                      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                        "DG density accum owner out of range: rank=", dg_frag%id, &
                        " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                        " igrid=", igrid, " owner=", owner_rank, " valid_remote=", valid_remote_grid_count
                      flush(6)
                      stop "DG-Fragment RT: density accum owner out of range"
                    end if
                    if (dg_frag%density_send_count(owner_rank) <= 0) then
                      if (owner_rank == dg_frag%id) cycle
                      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                        "DG density accum invalid send count: rank=", dg_frag%id, &
                        " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                        " igrid=", igrid, " owner=", owner_rank, " nsend=", dg_frag%density_send_count(owner_rank)
                      flush(6)
                      stop "DG-Fragment RT: density accum invalid send count"
                    end if
                    if (.not. allocated(rho_send(owner_rank)%f)) then
                      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                        "DG density accum send buffer missing: rank=", dg_frag%id, &
                        " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                        " igrid=", igrid, " owner=", owner_rank, " nsend=", dg_frag%density_send_count(owner_rank)
                      flush(6)
                      stop "DG-Fragment RT: density accum send buffer missing"
                    end if
                    if (slot < 1 .or. slot > dg_frag%density_send_count(owner_rank)) then
                      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                        "DG density accum slot out of range: rank=", dg_frag%id, &
                        " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                        " igrid=", igrid, " owner=", owner_rank, " slot=", slot, &
                        " nsend=", dg_frag%density_send_count(owner_rank)
                      flush(6)
                      stop "DG-Fragment RT: density accum slot out of range"
                    end if
                    rho_contrib = rho_blk_accum(igrid)
                    send_pre_block_raw = send_pre_block_raw + rho_contrib * system%hvol
                    if (enable_density_weight_path_trace .and. has_density_point_probe) then
                      ixg = ixg_buf(igrid)
                      iyg = iyg_buf(igrid)
                      izg = izg_buf(igrid)
                      if (ixg == probe_ixg .and. iyg == probe_iyg .and. izg == probe_izg) then
                        probe_import_send_pre_weight_local = probe_import_send_pre_weight_local + rho_contrib * system%hvol
                        probe_import_send_add_local = probe_import_send_add_local + rho_contrib * system%hvol
                        probe_import_send_apply_count_local = probe_import_send_apply_count_local + 1.0d0
                      end if
                    end if
!$omp atomic update
                    rho_send(owner_rank)%f(slot, 1, 1) = rho_send(owner_rank)%f(slot, 1, 1) + rho_contrib
                    spin_offset = ispin * dg_frag%density_send_count(owner_rank)
                    if (spin_offset + slot < 1 .or. spin_offset + slot > size(rho_send(owner_rank)%f, 1)) then
                      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                        "DG density accum spin slot out of range: rank=", dg_frag%id, &
                        " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                        " igrid=", igrid, " owner=", owner_rank, " slot=", slot, &
                        " spin_slot=", spin_offset + slot, " send_size=", size(rho_send(owner_rank)%f, 1)
                      flush(6)
                      stop "DG-Fragment RT: density accum spin slot out of range"
                    end if
!$omp atomic update
                    rho_send(owner_rank)%f(spin_offset + slot, 1, 1) = rho_send(owner_rank)%f(spin_offset + slot, 1, 1) + rho_contrib
                  end do
!$omp end do
!$omp end parallel
                call cpu_time(t_rho1)
                time_project_rho = time_project_rho + (t_rho1 - t_rho0)
            end if
          end do
        end do
        if (n_pw == 0 .and. enable_density_state_charge_trace) then
          call comm_summation(ifrag_owned_point_count_local, ifrag_owned_point_count_global, dg_frag%icomm)
          call comm_summation(ifrag_import_point_count_local, ifrag_import_point_count_global, dg_frag%icomm)
          do ispin = 1, system%nspin
            nocc_spin = dg_frag%nocc_spin(ispin)
            if (nocc_spin <= 0) cycle
            call comm_summation(state_charge_local(1:nocc_spin, ispin), state_charge_global(1:nocc_spin, ispin), nocc_spin, dg_frag%icomm)
            call comm_summation(state_coeff_c2_local(1:nocc_spin, ispin), state_coeff_c2_global(1:nocc_spin, ispin), nocc_spin, dg_frag%icomm)
            call comm_summation(state_psi2_raw_local(1:nocc_spin, ispin), state_psi2_raw_global(1:nocc_spin, ispin), nocc_spin, dg_frag%icomm)
            call comm_summation(state_psi2_occ_local(1:nocc_spin, ispin), state_psi2_occ_global(1:nocc_spin, ispin), nocc_spin, dg_frag%icomm)
            call comm_summation(state_psi2_dv_local(1:nocc_spin, ispin), state_psi2_dv_global(1:nocc_spin, ispin), nocc_spin, dg_frag%icomm)
            call comm_summation(state_psi2_owned_local(1:nocc_spin, ispin), state_psi2_owned_global(1:nocc_spin, ispin), nocc_spin, dg_frag%icomm)
            call comm_summation(state_import_core_local(1:nocc_spin, ispin), state_import_core_global(1:nocc_spin, ispin), nocc_spin, dg_frag%icomm)
            if (has_density_point_probe) then
              call comm_summation(probe_state_owned_local(1:nocc_spin, ispin), probe_state_owned_global(1:nocc_spin, ispin), nocc_spin, dg_frag%icomm)
              call comm_summation(probe_state_import_local(1:nocc_spin, ispin), probe_state_import_global(1:nocc_spin, ispin), nocc_spin, dg_frag%icomm)
            end if
          end do
          if (enable_state_rhobf_trace) then
            call comm_summation(state_rhobf_psi2_q_local, state_rhobf_psi2_q_global, dg_frag%icomm)
            call comm_summation(state_rhobf_psi2_occ_q_local, state_rhobf_psi2_occ_q_global, dg_frag%icomm)
            call comm_summation(state_rhobf_psi2_dv_q_local, state_rhobf_psi2_dv_q_global, dg_frag%icomm)
            call comm_summation(state_rhobf_psi2_owned_q_local, state_rhobf_psi2_owned_q_global, dg_frag%icomm)
            call comm_summation(state_rhobf_psi2_after_partition_q_local, state_rhobf_psi2_after_partition_q_global, dg_frag%icomm)
            call comm_summation(state_rhobf_psi2_after_slot_q_local, state_rhobf_psi2_after_slot_q_global, dg_frag%icomm)
            call comm_summation(state_rhobf_psi2_after_any_norm_q_local, state_rhobf_psi2_after_any_norm_q_global, dg_frag%icomm)
          end if
          if (dg_frag%id == 0) then
            state_charge_sum_all = 0.0d0
            state_coeff_c2_sum_all = 0.0d0
            do ispin = 1, system%nspin
              nocc_spin = dg_frag%nocc_spin(ispin)
              if (nocc_spin <= 0) cycle
              state_charge_sum_spin = sum(state_charge_global(1:nocc_spin, ispin))
              state_coeff_c2_sum_spin = sum(state_coeff_c2_global(1:nocc_spin, ispin))
              state_psi2_raw_sum_spin = sum(state_psi2_raw_global(1:nocc_spin, ispin))
              state_psi2_occ_sum_spin = sum(state_psi2_occ_global(1:nocc_spin, ispin))
              state_psi2_dv_sum_spin = sum(state_psi2_dv_global(1:nocc_spin, ispin))
              state_psi2_owned_sum_spin = sum(state_psi2_owned_global(1:nocc_spin, ispin))
              state_import_sum_spin = sum(state_import_core_global(1:nocc_spin, ispin))
              state_charge_sum_all = state_charge_sum_all + state_charge_sum_spin
              state_coeff_c2_sum_all = state_coeff_c2_sum_all + state_coeff_c2_sum_spin
              write(*,'(1x,a,i0,a,i0)') '        density state-grid-count: ifrag=', ifrag, ' ispin=', ispin
              write(*,'(1x,a,3(a,1pe12.4))') '        density state-grid-values:', &
                ' owned_pts=', ifrag_owned_point_count_global, ' imported_pts=', ifrag_import_point_count_global, &
                ' imported_frac=', merge(ifrag_import_point_count_global / (ifrag_owned_point_count_global + ifrag_import_point_count_global), 0.0d0, &
                  ifrag_owned_point_count_global + ifrag_import_point_count_global > 1.0d-14)
              write(*,'(1x,a,i0,a,i0,a,i0,4(a,1pe12.4))') '        density state-stage sum: ifrag=', ifrag, &
                ' ispin=', ispin, ' nocc=', nocc_spin, ' psi2_raw=', state_psi2_raw_sum_spin, &
                ' psi2_occ=', state_psi2_occ_sum_spin, ' psi2_dv=', state_psi2_dv_sum_spin, &
                ' psi2_owned=', state_psi2_owned_sum_spin
              write(*,'(1x,a,i0,a,i0,a,i0,3(a,1pe12.4))') '        density state-import sum: ifrag=', ifrag, &
                ' ispin=', ispin, ' nocc=', nocc_spin, ' imported_total=', state_import_sum_spin, &
                ' owned_plus_imported=', state_charge_sum_spin + state_import_sum_spin, &
                ' imported/(owned+imported)=', merge(state_import_sum_spin / (state_charge_sum_spin + state_import_sum_spin), 0.0d0, &
                  state_charge_sum_spin + state_import_sum_spin > 1.0d-14)
              write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4)') '        density state-charge sum: ifrag=', ifrag, &
                ' ispin=', ispin, ' nocc=', nocc_spin, ' owned_total=', state_charge_sum_spin
              write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4)') '        density state-coef sum: ifrag=', ifrag, &
                ' ispin=', ispin, ' nocc=', nocc_spin, ' c2_total=', state_coeff_c2_sum_spin
              do io = 1, min(5, nocc_spin)
                state_ratio_c2 = 0.0d0
                state_ratio_occ_raw = 0.0d0
                state_ratio_dv_occ = 0.0d0
                state_ratio_owned_dv = 0.0d0
                state_ratio_import_dv = 0.0d0
                state_ratio_total_dv = 0.0d0
                state_total_q = state_charge_global(io, ispin) + state_import_core_global(io, ispin)
                state_dev_q = state_total_q - 2.0d0
                state_ratio_total2 = state_total_q / 2.0d0
                if (state_coeff_c2_global(io, ispin) > 1.0d-14) state_ratio_c2 = state_charge_global(io, ispin) / state_coeff_c2_global(io, ispin)
                if (state_psi2_raw_global(io, ispin) > 1.0d-14) state_ratio_occ_raw = state_psi2_occ_global(io, ispin) / state_psi2_raw_global(io, ispin)
                if (state_psi2_occ_global(io, ispin) > 1.0d-14) state_ratio_dv_occ = state_psi2_dv_global(io, ispin) / state_psi2_occ_global(io, ispin)
                if (state_psi2_dv_global(io, ispin) > 1.0d-14) state_ratio_owned_dv = state_psi2_owned_global(io, ispin) / state_psi2_dv_global(io, ispin)
                if (state_psi2_dv_global(io, ispin) > 1.0d-14) state_ratio_import_dv = state_import_core_global(io, ispin) / state_psi2_dv_global(io, ispin)
                if (state_psi2_dv_global(io, ispin) > 1.0d-14) state_ratio_total_dv = state_total_q / state_psi2_dv_global(io, ispin)
                write(*,'(1x,a,i0,a,i0,4(a,1pe12.4))') '        density state-stage top: ifrag=', ifrag, &
                  ' io=', io, ' psi2_raw=', state_psi2_raw_global(io, ispin), ' psi2_occ=', state_psi2_occ_global(io, ispin), &
                  ' psi2_dv=', state_psi2_dv_global(io, ispin), ' psi2_owned=', state_psi2_owned_global(io, ispin)
                write(*,'(1x,a,i0,a,i0,4(a,1pe12.4))') '        density state-stage ratio: ifrag=', ifrag, &
                  ' io=', io, ' occ/raw=', state_ratio_occ_raw, ' dv/occ=', state_ratio_dv_occ, &
                  ' owned/dv=', state_ratio_owned_dv, ' owned/final=', merge(state_psi2_owned_global(io, ispin) / state_charge_global(io, ispin), 0.0d0, state_charge_global(io, ispin) > 1.0d-14)
                write(*,'(1x,a,i0,a,i0,3(a,1pe12.4))') '        density state-import core top: ifrag=', ifrag, &
                  ' io=', io, ' imported_q=', state_import_core_global(io, ispin), &
                  ' imported/owned=', merge(state_import_core_global(io, ispin) / state_charge_global(io, ispin), 0.0d0, state_charge_global(io, ispin) > 1.0d-14), &
                  ' imported/(owned+imported)=', merge(state_import_core_global(io, ispin) / (state_charge_global(io, ispin) + state_import_core_global(io, ispin)), 0.0d0, &
                    state_charge_global(io, ispin) + state_import_core_global(io, ispin) > 1.0d-14)
                write(*,'(1x,a,i0,a,i0,4(a,1pe12.4))') '        density state-total top: ifrag=', ifrag, &
                  ' io=', io, ' total_q=', state_total_q, ' total/2=', state_ratio_total2, &
                  ' dev_from_2=', state_dev_q, ' total/dv=', state_ratio_total_dv
                write(*,'(1x,a,i0,a,i0,3(a,1pe12.4))') '        density state-weight top: ifrag=', ifrag, &
                  ' io=', io, ' owned/dv=', state_ratio_owned_dv, ' imported/dv=', state_ratio_import_dv, &
                  ' imported_plus_owned/dv=', state_ratio_owned_dv + state_ratio_import_dv
                if (has_density_point_probe) then
                  write(*,'(1x,a,i0,a,i0,a,3(i0,1x),a,1pe12.4,a,1pe12.4,a,1pe12.4)') '        density point-probe state: ifrag=', ifrag, &
                    ' io=', io, ' idx=', probe_ixg, probe_iyg, probe_izg, &
                    ' owned_q=', probe_state_owned_global(io, ispin), ' imported_q=', probe_state_import_global(io, ispin), &
                    ' imported/owned=', merge(probe_state_import_global(io, ispin) / probe_state_owned_global(io, ispin), 0.0d0, probe_state_owned_global(io, ispin) > 1.0d-14)
                end if
                write(*,'(1x,a,i0,a,i0,a,1pe12.4)') '        density state-charge top: ifrag=', ifrag, &
                  ' io=', io, ' owned_q=', state_charge_global(io, ispin)
                write(*,'(1x,a,i0,a,i0,a,1pe12.4,a,1pe12.4)') '        density state-coef top: ifrag=', ifrag, &
                  ' io=', io, ' c2_q=', state_coeff_c2_global(io, ispin), ' ratio_q/c2=', state_ratio_c2
              end do
            end do
            write(*,'(1x,a,i0,a,1pe12.4)') '        density state-charge all-spin: ifrag=', ifrag, &
              ' owned_total=', state_charge_sum_all
            write(*,'(1x,a,i0,a,1pe12.4)') '        density state-coef all-spin: ifrag=', ifrag, &
              ' c2_total=', state_coeff_c2_sum_all
            if (enable_state_rhobf_trace) then
              if (state_rhobf_trace_spin >= 1 .and. state_rhobf_trace_spin <= system%nspin) then
                if (state_rhobf_trace_io >= 1 .and. state_rhobf_trace_io <= dg_frag%nocc_spin(state_rhobf_trace_spin)) then
                  state_rhobf_state_total_q = state_charge_global(state_rhobf_trace_io, state_rhobf_trace_spin) + &
                    state_import_core_global(state_rhobf_trace_io, state_rhobf_trace_spin)
                  write(*,'(1x,a,i0,a,i0,9(a,1pe12.4))') '        density state-rhobf compare: io=', state_rhobf_trace_io, &
                    ' ispin=', state_rhobf_trace_spin, &
                    ' state_total_q=', state_rhobf_state_total_q, &
                    ' psi2_q=', state_rhobf_psi2_q_global, &
                    ' psi2_occ_q=', state_rhobf_psi2_occ_q_global, &
                    ' psi2_dv_q=', state_rhobf_psi2_dv_q_global, &
                    ' psi2_owned_q=', state_rhobf_psi2_owned_q_global, &
                    ' psi2_after_partition_q=', state_rhobf_psi2_after_partition_q_global, &
                    ' psi2_after_slot_q=', state_rhobf_psi2_after_slot_q_global, &
                    ' psi2_after_any_norm_q=', state_rhobf_psi2_after_any_norm_q_global, &
                    ' path/state=', merge(state_rhobf_psi2_after_any_norm_q_global / state_rhobf_state_total_q, 0.0d0, state_rhobf_state_total_q > 1.0d-14)
                  write(*,'(1x,a,4(a,1pe12.4))') '        density state-rhobf stage-ratio:', &
                    ' occ/psi2=', merge(state_rhobf_psi2_occ_q_global / state_rhobf_psi2_q_global, 0.0d0, state_rhobf_psi2_q_global > 1.0d-14), &
                    ' dv/occ=', merge(state_rhobf_psi2_dv_q_global / state_rhobf_psi2_occ_q_global, 0.0d0, state_rhobf_psi2_occ_q_global > 1.0d-14), &
                    ' partition/owned=', merge(state_rhobf_psi2_after_partition_q_global / state_rhobf_psi2_owned_q_global, 0.0d0, state_rhobf_psi2_owned_q_global > 1.0d-14), &
                    ' anynorm/partition=', merge(state_rhobf_psi2_after_any_norm_q_global / state_rhobf_psi2_after_partition_q_global, 0.0d0, state_rhobf_psi2_after_partition_q_global > 1.0d-14)
                else
                  write(*,'(1x,a,i0,a,i0,a,i0)') '        density state-rhobf compare skipped: io=', state_rhobf_trace_io, &
                    ' ispin=', state_rhobf_trace_spin, ' nocc_spin=', dg_frag%nocc_spin(state_rhobf_trace_spin)
                end if
              end if
            end if
            flush(6)
          end if
        end if
        if (n_pw == 0 .and. enable_density_reconstruct_trace) then
          do ispin = 1, system%nspin
            phi_col_hvol_probe(:) = 0.0d0
            phi_sdiag_probe(:) = 0.0d0
            frag_state_trace_probe(:) = 0.0d0
            frag_state_real_probe(:) = 0.0d0
            nprobe_cols = min(3, dg_frag%n_basis(ifrag, ispin))
            do iprobe = 1, nprobe_cols
              phi_col_hvol_probe(iprobe) = phi_col_metric_total(iprobe, ispin, i_local)
              phi_sdiag_probe(iprobe) = basis_sdiag_probe(iprobe, ispin, i_local)
            end do
            do io = 1, min(3, dg_frag%nocc_spin(ispin))
              do iprobe = 1, dg_frag%n_basis(ifrag, ispin)
                do nocc_loc = 1, dg_frag%n_basis(ifrag, ispin)
                  frag_state_real_probe(io) = frag_state_real_probe(io) + &
                    (coef_re_full(iprobe, io, ispin) * coef_re_full(nocc_loc, io, ispin) + &
                     coef_im_full(iprobe, io, ispin) * coef_im_full(nocc_loc, io, ispin)) * &
                    phi_frag_metric_total(iprobe, nocc_loc, ispin, i_local)
                  if (basis_gid_spin(iprobe, ispin) > 0 .and. basis_gid_spin(nocc_loc, ispin) > 0) then
                    frag_state_trace_probe(io) = frag_state_trace_probe(io) + &
                      (coef_re_full(iprobe, io, ispin) * coef_re_full(nocc_loc, io, ispin) + &
                       coef_im_full(iprobe, io, ispin) * coef_im_full(nocc_loc, io, ispin)) * &
                      basis_frag_metric_total(iprobe, nocc_loc, ispin, i_local)
                  end if
                end do
              end do
            end do
            if (enable_density_reconstruct_trace) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,3(1pe12.4,1x),a,3(1pe12.4,1x))') "        density phi-metric probe: rank=", dg_frag%id, &
                " ifrag=", ifrag, " ispin=", ispin, " hvol_cols=", phi_col_hvol_probe(1), phi_col_hvol_probe(2), phi_col_hvol_probe(3), &
                " sdiag=", phi_sdiag_probe(1), phi_sdiag_probe(2), phi_sdiag_probe(3)
              flush(6)
              write(*,'(1x,a,i0,a,i0,a,i0,a,3(1pe12.4,1x),a,3(1pe12.4,1x))') "        density phi-gram probe: rank=", dg_frag%id, &
                " ifrag=", ifrag, " ispin=", ispin, " g12g13g23=", phi_gram_total(1,2,ispin,i_local), phi_gram_total(1,3,ispin,i_local), &
                phi_gram_total(2,3,ispin,i_local), " s12s13s23=", basis_smat_probe(1,2,ispin,i_local), basis_smat_probe(1,3,ispin,i_local), &
                basis_smat_probe(2,3,ispin,i_local)
              flush(6)
              write(*,'(1x,a,i0,a,i0,a,i0,a,3(1pe12.4,1x),a,3(1pe12.4,1x))') "        density frag-metric compare: rank=", dg_frag%id, &
                " ifrag=", ifrag, " ispin=", ispin, " real=", frag_state_real_probe(1), frag_state_real_probe(2), frag_state_real_probe(3), &
                " Smat=", frag_state_trace_probe(1), frag_state_trace_probe(2), frag_state_trace_probe(3)
              flush(6)
            end if
          end do
        end if
        block_idx_global = block_idx_global + nblocks_ifrag
      end do
    call cpu_time(t_project1)
    time_project = time_project + (t_project1 - t_project0)
    if (allocated(D_frag_re)) deallocate(D_frag_re)
    if (allocated(coef_re_frag)) deallocate(coef_re_frag)
    if (allocated(coef_im_frag)) deallocate(coef_im_frag)
    if (allocated(D_partial_re)) deallocate(D_partial_re)
    if (allocated(coef_re_full)) deallocate(coef_re_full)
    if (allocated(coef_im_full)) deallocate(coef_im_full)
    if (allocated(coef_c_full)) deallocate(coef_c_full)
    if (allocated(coef_c_frag)) deallocate(coef_c_frag)
    if (allocated(occ_cache)) deallocate(occ_cache)
    if (allocated(occ_sqrt_cache)) deallocate(occ_sqrt_cache)
    if (allocated(rho_blk_partial)) deallocate(rho_blk_partial)
    if (enable_density_reconstruct_trace .and. dg_frag%is_frag_root) then
      write(*,'(1x,a,i0,a,i0,a,1pe12.4)') "        density trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " stage=after-project dt=", time_project
      write(*,'(1x,a,i0,a,i0,8(a,1pe12.4))') "        density trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, &
        " setup=", time_project_setup, " psi=", time_project_psi, " rho=", time_project_rho, &
        " grid=", time_project_grid_prep, " phi=", time_project_phi_pack, &
        " over=", time_project_overhead, " dmat=", time_project_dmat_build, &
        " rho_red=", time_project_rho_reduce
      flush(6)
    end if

    if (enable_density_reconstruct_trace .and. dg_frag%is_frag_root) then
      write(*,'(1x,a,i0,a,i0,a)') "        density trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " stage=before-comm"
      flush(6)
    end if
    call cpu_time(t_comm0)
    allocate(send_counts(1:dg_frag%isize), recv_counts(1:dg_frag%isize))
    allocate(send_displs(1:dg_frag%isize), recv_displs(1:dg_frag%isize))
    send_counts = 0
    recv_counts = 0
    do irank = 0, dg_frag%isize - 1
      if (allocated(rho_send(irank)%f)) send_counts(irank + 1) = size(rho_send(irank)%f)
      if (allocated(rho_recv(irank)%f)) recv_counts(irank + 1) = size(rho_recv(irank)%f)
    end do
    send_displs(1) = 0
    recv_displs(1) = 0
    do irank = 2, dg_frag%isize
      send_displs(irank) = send_displs(irank - 1) + send_counts(irank - 1)
      recv_displs(irank) = recv_displs(irank - 1) + recv_counts(irank - 1)
    end do
    send_total_count = sum(send_counts)
    recv_total_count = sum(recv_counts)
    handler_send_peers = 0
    handler_recv_peers = 0
    do irank = 0, dg_frag%isize - 1
      if (send_counts(irank + 1) > 0) handler_send_peers = handler_send_peers + 1
      if (recv_counts(irank + 1) > 0) handler_recv_peers = handler_recv_peers + 1
      if (.not. dg_frag%parallel_mode_orbital) cycle
      if (send_counts(irank + 1) > 0 .and. .not. target_rank_owned_by_handler(irank)) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') '[FATAL] orbital send handler mismatch: rank=', dg_frag%id, &
          ' id_frag=', dg_frag%id_frag, ' target_rank=', irank, ' send_count=', send_counts(irank + 1)
        flush(6)
        stop 'DG-Fragment RT orbital mode: non-handler rank has send payload'
      end if
      if (recv_counts(irank + 1) > 0 .and. .not. target_rank_owned_by_handler(irank)) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') '[FATAL] orbital recv handler mismatch: rank=', dg_frag%id, &
          ' id_frag=', dg_frag%id_frag, ' source_rank=', irank, ' recv_count=', recv_counts(irank + 1)
        flush(6)
        stop 'DG-Fragment RT orbital mode: non-handler rank has recv payload'
      end if
    end do
    if (itt_tag == 1) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] DENSITY_ALLTOALLV_CONTRACT_SUM rank=', dg_frag%id, &
        ' itt=', itt_tag, ' send_total=', send_total_count, ' recv_total=', recv_total_count, &
        ' coef_ref_ready=', dg_frag%coef_ref_ready
      flush(6)
      do irank = 0, dg_frag%isize - 1
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] DENSITY_ALLTOALLV_CONTRACT_PEER rank=', &
          dg_frag%id, ' peer=', irank, ' send_count=', send_counts(irank + 1), ' recv_count=', recv_counts(irank + 1), &
          ' sdispl=', send_displs(irank + 1), ' rdispl=', recv_displs(irank + 1), ' coef_ref_ready=', dg_frag%coef_ref_ready
        flush(6)
      end do

      allocate(send_total_local_by_rank(1:dg_frag%isize), recv_total_local_by_rank(1:dg_frag%isize))
      allocate(send_total_all_ranks(1:dg_frag%isize), recv_total_all_ranks(1:dg_frag%isize))
      send_total_local_by_rank = 0
      recv_total_local_by_rank = 0
      send_total_local_by_rank(dg_frag%id + 1) = send_total_count
      recv_total_local_by_rank(dg_frag%id + 1) = recv_total_count
      call comm_summation(send_total_local_by_rank, send_total_all_ranks, dg_frag%isize, dg_frag%icomm)
      call comm_summation(recv_total_local_by_rank, recv_total_all_ranks, dg_frag%isize, dg_frag%icomm)
      if (dg_frag%id == 0) then
        write(*,'(1x,a)') '[DG-HANG-TRACE] DENSITY_ALLTOALLV_CONTRACT_GLOBAL_SUMMARY_BEGIN'
        do irank = 0, dg_frag%isize - 1
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') '[DG-HANG-TRACE] DENSITY_ALLTOALLV_CONTRACT_GLOBAL_SUM rank=', irank, &
            ' sum_sendcounts=', send_total_all_ranks(irank + 1), ' sum_recvcounts=', recv_total_all_ranks(irank + 1), &
            ' diff=', send_total_all_ranks(irank + 1) - recv_total_all_ranks(irank + 1)
        end do
        write(*,'(1x,a,i0,a,i0,a,i0)') '[DG-HANG-TRACE] DENSITY_ALLTOALLV_CONTRACT_GLOBAL_TOTAL sum_sendcounts=', &
          sum(send_total_all_ranks), ' sum_recvcounts=', sum(recv_total_all_ranks), &
          ' diff=', sum(send_total_all_ranks) - sum(recv_total_all_ranks)
        write(*,'(1x,a)') '[DG-HANG-TRACE] DENSITY_ALLTOALLV_CONTRACT_GLOBAL_SUMMARY_END'
        flush(6)
      end if
      deallocate(send_total_local_by_rank, recv_total_local_by_rank)
      deallocate(send_total_all_ranks, recv_total_all_ranks)
    end if
    if (enable_density_reconstruct_trace) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        density alltoallv summary: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " send_total=", send_total_count, " recv_total=", recv_total_count
      write(*,'(1x,a,i0,a,i0,a,l1,a,i0,a,i0)') "        density handler summary: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " orbital=", dg_frag%parallel_mode_orbital, &
        " send_peers=", handler_send_peers, " recv_peers=", handler_recv_peers
      flush(6)
      do irank = 0, dg_frag%isize - 1
        if (send_counts(irank + 1) > 0 .or. recv_counts(irank + 1) > 0) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
            "        density alltoallv peer: rank=", dg_frag%id, " id_frag=", dg_frag%id_frag, &
            " peer=", irank, " send_count=", send_counts(irank + 1), " recv_count=", recv_counts(irank + 1), &
            " send_pts=", dg_frag%density_send_count(irank), " recv_pts=", dg_frag%density_recv_map(irank)%npts
          flush(6)
        end if
      end do
    end if
    charge_owner_local_pre_comm = 0.0d0
!$omp parallel do collapse(2) reduction(+:charge_owner_local_pre_comm) private(ix,iy,iz) schedule(static)
    do iz = rho_z_lo, rho_z_hi
      do iy = rho_y_lo, rho_y_hi
!$omp simd reduction(+:charge_owner_local_pre_comm)
        do ix = rho_x_lo, rho_x_hi
          charge_owner_local_pre_comm = charge_owner_local_pre_comm + rho_bf(ix, iy, iz)
        end do
      end do
    end do
!$omp end parallel do
    charge_owner_local_pre_comm = charge_owner_local_pre_comm * system%hvol
    call comm_summation(charge_owner_local_pre_comm, charge_owner_global_pre_comm, dg_frag%icomm)
    dg_frag%rho_owned_elec = charge_owner_global_pre_comm
    dg_frag%rho_local_all_elec = charge_owner_global_pre_comm
    allocate(send_flat(max(1, send_total_count)), recv_flat(max(1, recv_total_count)))
    send_flat = 0.0d0
    recv_flat = 0.0d0
    call cpu_time(t_setup0)
    do irank = 0, dg_frag%isize - 1
      if (.not. allocated(rho_send(irank)%f)) cycle
      if (send_counts(irank + 1) <= 0) cycle
      send_flat(send_displs(irank + 1)+1:send_displs(irank + 1)+send_counts(irank + 1)) = rho_send(irank)%f(:, 1, 1)
      deallocate(rho_send(irank)%f)
    end do
    call cpu_time(t_setup1)
    time_comm_pack = time_comm_pack + (t_setup1 - t_setup0)
    if (enable_density_reconstruct_trace) then
      write(*,'(1x,a,i0,a)') "        density collective: rank=", dg_frag%id, " stage=before-alltoallv"
      flush(6)
    end if
    if (itt_tag == 1) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] BEFORE_DENSITY_ALLTOALLV rank=', dg_frag%id, &
        ' itt=', itt_tag, ' i_local=', i_local, ' coef_ref_ready=', dg_frag%coef_ref_ready
      flush(6)
    end if
    call cpu_time(t_setup0)
    call comm_alltoallv(send_flat, send_counts, send_displs, recv_flat, recv_counts, recv_displs, dg_frag%icomm)
    call cpu_time(t_setup1)
    time_comm_exchange = time_comm_exchange + (t_setup1 - t_setup0)
    if (enable_density_reconstruct_trace) then
      write(*,'(1x,a,i0,a)') "        density collective: rank=", dg_frag%id, " stage=after-alltoallv"
      flush(6)
    end if
    if (itt_tag == 1) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] AFTER_DENSITY_ALLTOALLV rank=', dg_frag%id, &
        ' itt=', itt_tag, ' i_local=', i_local, ' coef_ref_ready=', dg_frag%coef_ref_ready
      flush(6)
    end if
    call cpu_time(t_setup0)
    do irank = 0, dg_frag%isize - 1
      if (.not. allocated(rho_recv(irank)%f)) cycle
      if (recv_counts(irank + 1) > 0) then
        rho_recv(irank)%f(:, 1, 1) = recv_flat(recv_displs(irank + 1)+1:recv_displs(irank + 1)+recv_counts(irank + 1))
      end if
      npts = dg_frag%density_recv_map(irank)%npts
      if (recv_counts(irank + 1) /= (system%nspin + 1) * npts) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') &
          "[FATAL] density unpack recv size mismatch: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " peer=", irank, &
          " recv_count=", recv_counts(irank + 1), " npts=", npts
        flush(6)
        stop "DG-Fragment RT: density unpack recv size mismatch"
      end if
!$omp parallel do private(slot, ixg, iyg, izg, bx, by, bz, ispin, spin_offset, rho_contrib, rho_raw_contrib, imported_unpack_weight, i_local, ix0_frag, ix1_frag, iy0_frag, iy1_frag, iz0_frag, iz1_frag) schedule(static)
      do slot = 1, dg_frag%density_recv_map(irank)%npts
      ixg = dg_frag%density_recv_map(irank)%ixg(slot)
      iyg = dg_frag%density_recv_map(irank)%iyg(slot)
      izg = dg_frag%density_recv_map(irank)%izg(slot)
      bx = map_global_to_phi_box_coord_ham(ixg, lbound(rho_bf, 1), ubound(rho_bf, 1), dg_frag%lgnum_total(1))
      by = map_global_to_phi_box_coord_ham(iyg, lbound(rho_bf, 2), ubound(rho_bf, 2), dg_frag%lgnum_total(2))
      bz = map_global_to_phi_box_coord_ham(izg, lbound(rho_bf, 3), ubound(rho_bf, 3), dg_frag%lgnum_total(3))
      if (bx == 0 .or. by == 0 .or. bz == 0) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
          "[FATAL] density unpack rho_bf bounds: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " peer=", irank, " slot=", slot, &
          " idx=", ixg, iyg, izg, " lb=", lbound(rho_bf, 1), lbound(rho_bf, 2), lbound(rho_bf, 3), &
          " ub=", ubound(rho_bf, 1), ubound(rho_bf, 2), ubound(rho_bf, 3)
        flush(6)
        stop "DG-Fragment RT: density unpack rho_bf bounds"
      end if
      if (ixg < lbound(rho_s_bf, 1) .or. ixg > ubound(rho_s_bf, 1) .or. &
          iyg < lbound(rho_s_bf, 2) .or. iyg > ubound(rho_s_bf, 2) .or. &
          izg < lbound(rho_s_bf, 3) .or. izg > ubound(rho_s_bf, 3)) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
          "[FATAL] density unpack rho_s_bf bounds: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " peer=", irank, " slot=", slot, &
          " idx=", ixg, iyg, izg, " lb=", lbound(rho_s_bf, 1), lbound(rho_s_bf, 2), lbound(rho_s_bf, 3), &
          " ub=", ubound(rho_s_bf, 1), ubound(rho_s_bf, 2), ubound(rho_s_bf, 3)
        flush(6)
        stop "DG-Fragment RT: density unpack rho_s_bf bounds"
      end if
      rho_raw_contrib = rho_recv(irank)%f(slot, 1, 1)
      if (enable_ifrag_compare_trace .and. itt_tag >= 8 .and. itt_tag <= 10) then
        do i_local = 1, ifrag_count
          ix0_frag = dg_frag%ixyz_frag(1, dg_frag%ifrag_start + i_local - 1)
          iy0_frag = dg_frag%ixyz_frag(2, dg_frag%ifrag_start + i_local - 1)
          iz0_frag = dg_frag%ixyz_frag(3, dg_frag%ifrag_start + i_local - 1)
          ix1_frag = ix0_frag + dg_frag%nxyz_domain(1, dg_frag%ifrag_start + i_local - 1) - 1
          iy1_frag = iy0_frag + dg_frag%nxyz_domain(2, dg_frag%ifrag_start + i_local - 1) - 1
          iz1_frag = iz0_frag + dg_frag%nxyz_domain(3, dg_frag%ifrag_start + i_local - 1) - 1
          if (ixg < ix0_frag .or. ixg > ix1_frag) cycle
          if (iyg < iy0_frag .or. iyg > iy1_frag) cycle
          if (izg < iz0_frag .or. izg > iz1_frag) cycle
!$omp atomic update
          ifrag_recv_post_raw_local(1, i_local) = ifrag_recv_post_raw_local(1, i_local) + rho_raw_contrib * system%hvol
        end do
      end if
      imported_unpack_weight = 1.0d0
      if (enable_imported_unpack_normalized) then
        imported_unpack_weight = imported_unpack_scale
        if (imported_unpack_weight < 0.999999d0) imported_unpack_norm_trigger_local = imported_unpack_norm_trigger_local + 1.0d0
      end if
      rho_contrib = rho_raw_contrib * imported_unpack_weight
      if (enable_density_weight_path_trace .and. has_density_point_probe) then
        if (ixg == probe_ixg .and. iyg == probe_iyg .and. izg == probe_izg) then
          probe_import_unpack_pre_weight_local = probe_import_unpack_pre_weight_local + rho_raw_contrib * system%hvol
          probe_import_unpack_add_local = probe_import_unpack_add_local + rho_contrib * system%hvol
          probe_import_unpack_apply_count_local = probe_import_unpack_apply_count_local + 1.0d0
          probe_import_unpack_weight_sum_local = probe_import_unpack_weight_sum_local + imported_unpack_weight
        end if
      end if
      if (enable_density_point_dup_audit) then
        if (ixg >= rho_x_lo .and. ixg <= rho_x_hi .and. &
            iyg >= rho_y_lo .and. iyg <= rho_y_hi .and. &
            izg >= rho_z_lo .and. izg <= rho_z_hi) then
          if (first_imp_ixg_local < 0) then
            first_imp_ixg_local = ixg
            first_imp_iyg_local = iyg
            first_imp_izg_local = izg
            first_imp_src_local = irank
            first_imp_tgt_local = dg_frag%id
            first_imp_slot_local = slot
          end if
          if (.not. found_duplicate_point_local) then
            if (rho_bf(bx, by, bz) * system%hvol > 1.0d-14 .and. rho_raw_contrib * system%hvol > 1.0d-14) then
              dup_ixg_local = ixg
              dup_iyg_local = iyg
              dup_izg_local = izg
              dup_src_local = irank
              dup_tgt_local = dg_frag%id
              dup_slot_local = slot
              dup_local_contrib = rho_bf(bx, by, bz) * system%hvol
              dup_import_contrib = rho_raw_contrib * system%hvol
              found_duplicate_point_local = .true.
            end if
          end if
        end if
      end if
      rho_bf(bx, by, bz) = rho_bf(bx, by, bz) + rho_contrib
      if (ixg >= rho_x_lo .and. ixg <= rho_x_hi .and. &
          iyg >= rho_y_lo .and. iyg <= rho_y_hi .and. &
          izg >= rho_z_lo .and. izg <= rho_z_hi) then
!$omp atomic update
        charge_imported_core_local = charge_imported_core_local + rho_contrib * system%hvol
      else
!$omp atomic update
        charge_imported_buffer_local = charge_imported_buffer_local + rho_contrib * system%hvol
      end if
      do ispin = 1, system%nspin
        spin_offset = ispin * npts
        if (spin_offset + slot < 1 .or. spin_offset + slot > size(rho_recv(irank)%f, 1)) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
            "[FATAL] density unpack spin slot bounds: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " peer=", irank, " slot=", slot, &
            " spin_slot=", spin_offset + slot, " recv_size=", size(rho_recv(irank)%f, 1)
          flush(6)
          stop "DG-Fragment RT: density unpack spin slot bounds"
        end if
        rho_contrib = rho_recv(irank)%f(spin_offset + slot, 1, 1) * imported_unpack_weight
        rho_s_bf(ixg, iyg, izg, ispin) = rho_s_bf(ixg, iyg, izg, ispin) + rho_contrib
      end do
      end do
!$omp end parallel do
      deallocate(rho_recv(irank)%f)
    end do
    call cpu_time(t_setup1)
    time_comm_unpack = time_comm_unpack + (t_setup1 - t_setup0)
    if (itt_tag == 1) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] AFTER_DENSITY_UNPACK rank=', dg_frag%id, &
        ' itt=', itt_tag, ' i_local=', i_local, ' coef_ref_ready=', dg_frag%coef_ref_ready
      flush(6)
    end if
    if (enable_density_point_dup_audit) then
      if (found_duplicate_point_local) then
        write(*,'(1x,a,i0,a,3(i0,1x),a,1pe12.4,a,1pe12.4,a,i0,a,i0,a,i0)') '        density point-dup first: rank=', dg_frag%id, ' idx=', &
          dup_ixg_local, dup_iyg_local, dup_izg_local, ' pre_q=', dup_local_contrib, ' imported_q=', dup_import_contrib, &
          ' slot=', dup_slot_local, ' source_rank=', dup_src_local, ' target_rank=', dup_tgt_local
        flush(6)
      else if (first_imp_ixg_local > 0) then
        write(*,'(1x,a,i0,a,3(i0,1x),a,i0,a,i0,a,i0)') '        density point-import first: rank=', dg_frag%id, ' idx=', &
          first_imp_ixg_local, first_imp_iyg_local, first_imp_izg_local, ' slot=', first_imp_slot_local, &
          ' source_rank=', first_imp_src_local, ' target_rank=', first_imp_tgt_local
        flush(6)
      end if
    end if
    if (enable_ifrag_compare_trace .and. itt_tag >= 8 .and. itt_tag <= 10) then
      call comm_summation(ifrag_recv_post_raw_local, ifrag_recv_post_raw_global, size(ifrag_recv_post_raw_local), dg_frag%icomm)
      if (dg_frag%id == 0 .and. ifrag_count >= 2) then
        cmp_recv_post_raw(1) = ifrag_recv_post_raw_global(1, 1)
        cmp_recv_post_raw(2) = ifrag_recv_post_raw_global(1, 2)
        write(*,'(1x,a,i0,a,1pe12.4,a,1pe12.4)') '        ifrag-compare recv_post_raw: itt=', itt_tag, &
          ' f1=', cmp_recv_post_raw(1), ' f2=', cmp_recv_post_raw(2)
        flush(6)
      end if
    end if
    deallocate(send_flat, recv_flat)
    deallocate(send_counts, recv_counts, send_displs, recv_displs)
    call cpu_time(t_comm1)
    time_comm = time_comm + (t_comm1 - t_comm0)
    if (enable_density_reconstruct_trace .and. dg_frag%is_frag_root) then
      write(*,'(1x,a,i0,a,i0,a,1pe12.4)') "        density trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " stage=after-comm dt=", time_comm
      write(*,'(1x,a,i0,a,i0,4(a,1pe12.4))') "        density trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, &
        " subgroup_reduce=", time_comm_subgroup_reduce, " pack=", time_comm_pack, &
        " exchange=", time_comm_exchange, " unpack=", time_comm_unpack
      flush(6)
    end if

    if (enable_density_reconstruct_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density trace: stage=before-normalize"
      flush(6)
    end if
    charge_rho_bf_all_local_raw = 0.0d0
!$omp parallel do collapse(2) reduction(+:charge_rho_bf_all_local_raw) private(ix,iy,iz) schedule(static)
    do iz = lbound(rho_bf, 3), ubound(rho_bf, 3)
      do iy = lbound(rho_bf, 2), ubound(rho_bf, 2)
!$omp simd reduction(+:charge_rho_bf_all_local_raw)
        do ix = lbound(rho_bf, 1), ubound(rho_bf, 1)
          charge_rho_bf_all_local_raw = charge_rho_bf_all_local_raw + rho_bf(ix, iy, iz)
        end do
      end do
    end do
!$omp end parallel do
    call comm_summation(charge_rho_bf_all_local_raw, charge_rho_bf_all_global_raw, dg_frag%icomm)

    ! rho_bf -> rho_s boundary: only the authoritative mg-local density is
    ! materialized below for downstream Hartree/XC/reconstruct consumers.
    if (lbound(rho%f, 1) /= rho_x_lo .or. ubound(rho%f, 1) /= rho_x_hi .or. &
        lbound(rho%f, 2) /= rho_y_lo .or. ubound(rho%f, 2) /= rho_y_hi .or. &
        lbound(rho%f, 3) /= rho_z_lo .or. ubound(rho%f, 3) /= rho_z_hi) then
      write(*,'(1x,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
        "[FATAL] density rho core bounds mismatch: rank=", dg_frag%id, &
        " rho_lb=", lbound(rho%f, 1), lbound(rho%f, 2), lbound(rho%f, 3), &
        " rho_ub=", ubound(rho%f, 1), ubound(rho%f, 2), ubound(rho%f, 3), &
        " core_lo=", rho_x_lo, rho_y_lo, rho_z_lo, &
        " core_hi=", rho_x_hi, rho_y_hi, rho_z_hi
      flush(6)
      stop "DG-Fragment RT: density rho core bounds mismatch"
    end if
    rho%f(rho_x_lo:rho_x_hi, rho_y_lo:rho_y_hi, rho_z_lo:rho_z_hi) = &
      rho_bf(rho_x_lo:rho_x_hi, rho_y_lo:rho_y_hi, rho_z_lo:rho_z_hi)
    do ispin = 1, system%nspin
      if (lbound(rho_s(ispin)%f, 1) /= lbound(rho_s_bf, 1) .or. ubound(rho_s(ispin)%f, 1) /= ubound(rho_s_bf, 1) .or. &
          lbound(rho_s(ispin)%f, 2) /= lbound(rho_s_bf, 2) .or. ubound(rho_s(ispin)%f, 2) /= ubound(rho_s_bf, 2) .or. &
          lbound(rho_s(ispin)%f, 3) /= lbound(rho_s_bf, 3) .or. ubound(rho_s(ispin)%f, 3) /= ubound(rho_s_bf, 3)) then
        write(*,'(1x,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
          "[FATAL] density rho_s copy bounds mismatch: rank=", dg_frag%id, " ispin=", ispin, &
          " rho_lb=", lbound(rho_s(ispin)%f, 1), lbound(rho_s(ispin)%f, 2), lbound(rho_s(ispin)%f, 3), &
          " rho_ub=", ubound(rho_s(ispin)%f, 1), ubound(rho_s(ispin)%f, 2), ubound(rho_s(ispin)%f, 3), &
          " bf_lb=", lbound(rho_s_bf, 1), lbound(rho_s_bf, 2), lbound(rho_s_bf, 3), &
          " bf_ub=", ubound(rho_s_bf, 1), ubound(rho_s_bf, 2), ubound(rho_s_bf, 3)
        flush(6)
        stop "DG-Fragment RT: density rho_s copy bounds mismatch"
      end if
      rho_s(ispin)%f(:, :, :) = rho_s_bf(:, :, :, ispin)
    end do
    if (enable_density_reconstruct_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density trace: stage=after-rho-copy"
      flush(6)
    end if
    charge_rho_bf_core_local_raw = 0.0d0
!$omp parallel do collapse(2) reduction(+:charge_rho_bf_core_local_raw) private(ix,iy,iz) schedule(static)
    do iz = rho_z_lo, rho_z_hi
      do iy = rho_y_lo, rho_y_hi
!$omp simd reduction(+:charge_rho_bf_core_local_raw)
        do ix = rho_x_lo, rho_x_hi
          charge_rho_bf_core_local_raw = charge_rho_bf_core_local_raw + rho_bf(ix, iy, iz)
        end do
      end do
    end do
!$omp end parallel do
    call comm_summation(charge_rho_bf_core_local_raw, charge_rho_bf_core_global_raw, dg_frag%icomm)
    if (enable_density_reconstruct_trace) then
      call comm_summation(orbital_norm_probe_local, orbital_norm_probe_total, 3, dg_frag%icomm)
      call comm_summation(orbital_norm_frag_local, orbital_norm_frag_total, size(orbital_norm_frag_local), dg_frag%icomm)
    end if
    call cpu_time(t_norm0)
    total_charge_local = 0.0d0
!$omp parallel do collapse(2) reduction(+:total_charge_local) private(ix,iy,iz) schedule(static)
    do iz = rho_z_lo, rho_z_hi
      do iy = rho_y_lo, rho_y_hi
!$omp simd reduction(+:total_charge_local)
        do ix = rho_x_lo, rho_x_hi
          total_charge_local = total_charge_local + rho%f(ix, iy, iz)
        end do
      end do
    end do
!$omp end parallel do
    total_charge_local = total_charge_local * system%hvol
    charge_weighted_total_pre_norm = total_charge_local
    charge_imported_local = total_charge_local - charge_owner_local_pre_comm
    call comm_summation(charge_imported_local, charge_imported_global, dg_frag%icomm)
    call comm_summation(charge_imported_core_local, charge_imported_core_global, dg_frag%icomm)
    call comm_summation(charge_imported_buffer_local, charge_imported_buffer_global, dg_frag%icomm)
    if (enable_density_weight_path_trace .and. has_density_point_probe) then
      call comm_summation(probe_owned_pre_weight_local, probe_owned_pre_weight_global, dg_frag%icomm)
      call comm_summation(probe_owned_add_local, probe_owned_add_global, dg_frag%icomm)
      call comm_summation(probe_import_send_pre_weight_local, probe_import_send_pre_weight_global, dg_frag%icomm)
      call comm_summation(probe_import_send_add_local, probe_import_send_add_global, dg_frag%icomm)
      call comm_summation(probe_import_unpack_pre_weight_local, probe_import_unpack_pre_weight_global, dg_frag%icomm)
      call comm_summation(probe_import_unpack_add_local, probe_import_unpack_add_global, dg_frag%icomm)
      call comm_summation(probe_owned_apply_count_local, probe_owned_apply_count_global, dg_frag%icomm)
      call comm_summation(probe_owned_weight_sum_local, probe_owned_weight_sum_global, dg_frag%icomm)
      call comm_summation(probe_import_send_apply_count_local, probe_import_send_apply_count_global, dg_frag%icomm)
      call comm_summation(probe_import_unpack_apply_count_local, probe_import_unpack_apply_count_global, dg_frag%icomm)
      call comm_summation(probe_import_unpack_weight_sum_local, probe_import_unpack_weight_sum_global, dg_frag%icomm)
      call comm_summation(imported_unpack_norm_trigger_local, imported_unpack_norm_trigger_global, dg_frag%icomm)
    end if
    dg_frag%rho_imported_elec = charge_imported_global
    if (enable_density_reconstruct_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a,a,1pe12.4)') "        density charge budget:", &
        " weighted_pre_norm=", charge_weighted_total_pre_norm
      flush(6)
      write(*,'(1x,a,3(1pe12.4,1x))') "        density orbital self-norm probe: ", &
        orbital_norm_probe_total(1), orbital_norm_probe_total(2), orbital_norm_probe_total(3)
      flush(6)
      write(*,'(1x,a,3(1pe12.4,1x),a,3(1pe12.4,1x))') "        density orbital frag-norm probe: frag1=", &
        orbital_norm_frag_total(1,1), orbital_norm_frag_total(1,2), orbital_norm_frag_total(1,3), &
        " frag2=", orbital_norm_frag_total(2,1), orbital_norm_frag_total(2,2), orbital_norm_frag_total(2,3)
      flush(6)
    end if
    if (enable_density_reconstruct_trace) then
      write(*,'(1x,a,i0,a,i0,a,1pe12.4)') "        density collective: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " stage=before-total-charge-sum total_charge_local=", total_charge_local
      flush(6)
    end if
    call comm_summation(total_charge_local, total_charge, dg_frag%icomm)
    charge_weighted_total_global = total_charge
    dg_frag%rho_global_raw_elec = total_charge
    charge_contract_residual = total_charge - (charge_owner_global_pre_comm + charge_imported_global)
    dg_frag%rho_contract_residual_elec = charge_contract_residual
    if (enable_density_stage_contrib_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a,6(1x,a,1pe12.4))') '        density stage contrib:', &
        'owned_local=', charge_owner_global_pre_comm, &
        'imported_halo=', charge_imported_global, &
        'imported_core=', charge_imported_core_global, &
        'imported_buffer=', charge_imported_buffer_global, &
        'overlap_corrected=', charge_owner_global_pre_comm + charge_imported_global, &
        'final_integrated=', total_charge
      write(*,'(1x,a,10(1x,a,1pe12.4))') '        density integration audit:', &
        'rho_bf_all_raw=', charge_rho_bf_all_global_raw, &
        'rho_bf_all_q=', charge_rho_bf_all_global_raw * system%hvol, &
        'rho_core_raw=', charge_rho_bf_core_global_raw, &
        'rho_core_q=', charge_rho_bf_core_global_raw * system%hvol, &
        'owner_mask_q=', charge_owner_global_pre_comm, &
        'rho_sum_raw_core=', total_charge_local / system%hvol, &
        'rho_sum_q_core=', total_charge_local, &
        'global_final_q=', total_charge, &
        'hvol=', system%hvol, &
        'nspin=', dble(system%nspin)
      flush(6)
    end if
    if (enable_density_weight_path_trace .and. has_density_point_probe .and. dg_frag%id == 0) then
      if (enable_owned_path_normalized) then
        write(*,'(1x,a,a,a,1pe12.4)') '        density weight-probe owned-mode: ', 'normalized', ' scale=', owned_path_scale
      else
        write(*,'(1x,a,a)') '        density weight-probe owned-mode: ', 'legacy'
      end if
      if (enable_imported_unpack_normalized) then
        write(*,'(1x,a,a,a,1pe12.4)') '        density weight-probe unpack-mode: ', 'normalized', ' scale=', imported_unpack_scale
      else
        write(*,'(1x,a,a)') '        density weight-probe unpack-mode: ', 'legacy'
      end if
      write(*,'(1x,a,3(i0,1x),3(a,1pe12.4))') '        density weight-probe point: idx=', &
        probe_ixg, probe_iyg, probe_izg, ' partition_w=', probe_partition_weight, ' overlap_w=', probe_overlap_weight, &
        ' norm_w=', probe_norm_weight
      write(*,'(1x,a,2(a,1pe12.4),a,1pe12.4)') '        density weight-probe owned-path:', &
        ' pre_send_q=', probe_owned_pre_weight_global, ' add_q=', probe_owned_add_global, ' apply_count=', probe_owned_apply_count_global
      write(*,'(1x,a,a,1pe12.4)') '        density weight-probe owned-weight:', ' avg_w=', &
        merge(probe_owned_weight_sum_global / probe_owned_apply_count_global, 0.0d0, &
          probe_owned_apply_count_global > 1.0d-14)
      write(*,'(1x,a,2(a,1pe12.4),a,1pe12.4)') '        density weight-probe imported-send-path:', &
        ' pre_send_q=', probe_import_send_pre_weight_global, ' add_to_send_q=', probe_import_send_add_global, &
        ' apply_count=', probe_import_send_apply_count_global
      write(*,'(1x,a,2(a,1pe12.4),a,1pe12.4)') '        density weight-probe imported-unpack-path:', &
        ' recv_q=', probe_import_unpack_pre_weight_global, ' add_after_unpack_q=', probe_import_unpack_add_global, &
        ' apply_count=', probe_import_unpack_apply_count_global
      write(*,'(1x,a,a,1pe12.4)') '        density weight-probe imported-unpack-weight:', ' avg_w=', &
        merge(probe_import_unpack_weight_sum_global / probe_import_unpack_apply_count_global, 0.0d0, &
          probe_import_unpack_apply_count_global > 1.0d-14)
      write(*,'(1x,a,1pe12.4)') '        density weight-probe imported-unpack-normalized-trigger-count=', &
        imported_unpack_norm_trigger_global
      flush(6)
    end if
    if (enable_density_reconstruct_trace) then
      write(*,'(1x,a,i0,a,i0,a,1pe12.4)') "        density collective: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " stage=after-total-charge-sum total_charge=", total_charge
      flush(6)
    end if
    dg_frag%elec_num_raw = total_charge
    dg_frag%rho_scale_factor = 1.0d0
    if (total_charge > 1.0d-14 .and. total_charge == total_charge) then
      ! Normalize reconstructed real-space density to target total electrons (nelec).
      ! Therefore elec_num_scaled follows the chosen electron-count convention and can
      ! remain fixed (e.g., 40) even when elec_num_raw varies over time.
      scale_rho = nelec / total_charge
      dg_frag%rho_scale_factor = scale_rho
      if (abs(scale_rho - 1.0d0) > 1.0d-14) then
        if (system%nspin == 1) then
!$omp parallel do collapse(3) private(ix,iy,iz,rho_sum_local) schedule(static)
          do iz = rho_z_lo, rho_z_hi
            do iy = rho_y_lo, rho_y_hi
              do ix = rho_x_lo, rho_x_hi
                rho_sum_local = rho_s(1)%f(ix, iy, iz) * scale_rho
                rho_s(1)%f(ix, iy, iz) = rho_sum_local
                rho%f(ix, iy, iz) = rho_sum_local
              end do
            end do
          end do
!$omp end parallel do
        else
!$omp parallel do collapse(3) private(ix,iy,iz,ispin,rho_sum_local) schedule(static)
          do iz = rho_z_lo, rho_z_hi
            do iy = rho_y_lo, rho_y_hi
              do ix = rho_x_lo, rho_x_hi
                rho_sum_local = 0.0d0
                do ispin = 1, system%nspin
                  rho_s(ispin)%f(ix, iy, iz) = rho_s(ispin)%f(ix, iy, iz) * scale_rho
                  rho_sum_local = rho_sum_local + rho_s(ispin)%f(ix, iy, iz)
                end do
                rho%f(ix, iy, iz) = rho_sum_local
              end do
            end do
          end do
!$omp end parallel do
        end if
      end if
      dg_frag%elec_num_scaled = total_charge * scale_rho
    else
      dg_frag%elec_num_scaled = total_charge
    end if
    dg_frag%rho_global_scaled_elec = dg_frag%elec_num_scaled
    if (enable_density_reconstruct_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a,i0,a)') "        density collective: rank=", dg_frag%id, " stage=before-scaled-charge-sum"
      flush(6)
      write(*,'(1x,a,i0,a,1pe12.4)') "        density collective: rank=", dg_frag%id, &
        " stage=after-scaled-charge-sum elec_num_scaled=", dg_frag%elec_num_scaled
      flush(6)
    end if
    call cpu_time(t_norm1)
    time_norm = time_norm + (t_norm1 - t_norm0)
    if (enable_density_reconstruct_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a,a,1pe12.4)') "        density trace: stage=after-normalize dt=", "", time_norm
      flush(6)
    end if

    deallocate(ix_buf, iy_buf, iz_buf, owner_buf, ixg_buf, iyg_buf, izg_buf)
    deallocate(slot_buf, local_grid_ids, remote_grid_ids, valid_remote_grid_ids)
    deallocate(basis_gid, valid_basis_ids)
    deallocate(state_charge_local, state_charge_global)
    deallocate(state_coeff_c2_local, state_coeff_c2_global)
    deallocate(state_psi2_raw_local, state_psi2_raw_global)
    deallocate(state_psi2_occ_local, state_psi2_occ_global)
    deallocate(state_psi2_dv_local, state_psi2_dv_global)
    deallocate(state_psi2_owned_local, state_psi2_owned_global)
    deallocate(state_import_core_local, state_import_core_global)
    deallocate(probe_state_owned_local, probe_state_owned_global)
    deallocate(probe_state_import_local, probe_state_import_global)
    if (allocated(basis_sdiag_probe)) deallocate(basis_sdiag_probe)
    if (allocated(phi_col_metric_total)) deallocate(phi_col_metric_total)
    if (allocated(basis_smat_probe)) deallocate(basis_smat_probe)
    if (allocated(phi_gram_total)) deallocate(phi_gram_total)
    if (allocated(phi_frag_metric_total)) deallocate(phi_frag_metric_total)
    if (allocated(basis_frag_metric_total)) deallocate(basis_frag_metric_total)
    deallocate(phi_blk, rho_blk, rho_blk_accum, rho_blk_reduced, coef_blk_re, coef_blk_im, psi_blk_re, psi_blk_im, &
      density_mat_re, density_tmp)
    if (allocated(coef_blk_ri)) deallocate(coef_blk_ri)
    if (allocated(psi_blk_ri)) deallocate(psi_blk_ri)
    if (allocated(pw_tmp_z)) deallocate(pw_tmp_z)
    if (allocated(phase_cache)) deallocate(phase_cache)
    if (allocated(kpw_hx)) deallocate(kpw_hx, kpw_hy, kpw_hz)
    if (allocated(density_mix)) deallocate(density_mix)
    if (allocated(basis_mix_blk)) deallocate(basis_mix_blk)
    if (allocated(density_mix_tmp)) deallocate(density_mix_tmp)
    if (allocated(basis_mix_blk_t)) deallocate(basis_mix_blk_t)
    if (allocated(density_mix_tmp_t)) deallocate(density_mix_tmp_t)
    if (allocated(transform_frag_spin)) deallocate(transform_frag_spin)
    if (allocated(transform_pw_spin)) deallocate(transform_pw_spin)
    if (allocated(mix_transform_spin)) deallocate(mix_transform_spin)
    if (allocated(mix_overlap_spin)) deallocate(mix_overlap_spin)
    if (allocated(s_mix)) deallocate(s_mix)
    if (allocated(s_mix_work)) deallocate(s_mix_work)
    if (allocated(coef_mix_eff)) deallocate(coef_mix_eff)
    if (allocated(coef_mix_metric)) deallocate(coef_mix_metric)
    if (allocated(coef_mix_spin)) deallocate(coef_mix_spin)
    if (allocated(d_raw_ff)) deallocate(d_raw_ff)
    if (allocated(d_raw_fp)) deallocate(d_raw_fp)
    if (allocated(d_raw_pp)) deallocate(d_raw_pp)
    if (allocated(ipiv_mix)) deallocate(ipiv_mix)
    if (allocated(n_basis_mix_spin)) deallocate(n_basis_mix_spin)
    if (allocated(g211_cos_x)) deallocate(g211_cos_x, g211_sin_x)
    if (allocated(ifrag_recv_post_raw_local)) deallocate(ifrag_recv_post_raw_local)
    if (allocated(ifrag_recv_post_raw_global)) deallocate(ifrag_recv_post_raw_global)
    if (allocated(cmp_ff_occ)) deallocate(cmp_ff_occ)
    if (allocated(cmp_ff_occ_global)) deallocate(cmp_ff_occ_global)
    if (allocated(cmp_ff_dom_gid)) deallocate(cmp_ff_dom_gid)
    if (allocated(cmp_ff_dom_gid_local)) deallocate(cmp_ff_dom_gid_local)
    if (allocated(cmp_ff_dom_gid_global)) deallocate(cmp_ff_dom_gid_global)
    if (allocated(cmp_ff_gid_pre)) deallocate(cmp_ff_gid_pre)
    if (allocated(cmp_ff_gid_post)) deallocate(cmp_ff_gid_post)
    if (allocated(cmp_ff_gid_pre_global)) deallocate(cmp_ff_gid_pre_global)
    if (allocated(cmp_ff_gid_post_global)) deallocate(cmp_ff_gid_post_global)
    if (allocated(cmp_tf_gid_pre)) deallocate(cmp_tf_gid_pre)
    if (allocated(cmp_tf_gid_full)) deallocate(cmp_tf_gid_full)
    if (allocated(cmp_tf_gid_int)) deallocate(cmp_tf_gid_int)
    if (allocated(cmp_tf_gid_pre_global)) deallocate(cmp_tf_gid_pre_global)
    if (allocated(cmp_tf_gid_full_global)) deallocate(cmp_tf_gid_full_global)
    if (allocated(cmp_tf_gid_int_global)) deallocate(cmp_tf_gid_int_global)
    if (allocated(cmp_tf_gid_mode_pre)) deallocate(cmp_tf_gid_mode_pre)
    if (allocated(cmp_tf_gid_mode_full)) deallocate(cmp_tf_gid_mode_full)
    if (allocated(cmp_tf_gid_mode_pre_global)) deallocate(cmp_tf_gid_mode_pre_global)
    if (allocated(cmp_tf_gid_mode_full_global)) deallocate(cmp_tf_gid_mode_full_global)
    if (allocated(cmp_tf_gid_mode_ovl)) deallocate(cmp_tf_gid_mode_ovl)
    if (allocated(cmp_tf_gid_mode_ovl_global)) deallocate(cmp_tf_gid_mode_ovl_global)
    deallocate(rho_bf, rho_s_bf)
    deallocate(rho_send, rho_recv)
    call cpu_time(t_total1)
    if (enable_density_reconstruct_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
        "        density timing: total=", t_total1 - t_total0, " cache=", time_cache, &
        " project=", time_project, " comm=", time_comm, " norm=", time_norm
      write(*,'(1x,a,4(a,1pe12.4),2(a,l1))') "        density cache timing:", &
        " pw_refresh=", time_cache_pw_refresh, " phi_blk_build=", time_cache_phi_block_refresh, &
        " cache_total=", time_cache + time_cache_phi_block_refresh, " step_total=", time_cache_pw_refresh + time_cache_phi_block_refresh, &
        " pw_rebuilt=", rebuilt_pw_cache, " phi_rebuilt=", rebuilt_phi_block_cache
      write(*,'(1x,a,7(a,1pe12.4))') "        density timing detail: project", &
        " setup=", time_project_setup, " psi=", time_project_psi, " rho=", time_project_rho, &
        " grid=", time_project_grid_prep, " phi=", time_project_phi_pack, &
        " phi_blk_build=", time_project_phi_block_build, " over=", time_project_overhead
      flush(6)
    end if
    if (itt_tag == 1) then
      write(*,'(1x,a,i0,a,i0,a,l1)') '[DG-HANG-TRACE] EXIT_DENSITY_RECON rank=', dg_frag%id, ' itt=', itt_tag, &
        ' coef_ref_ready=', dg_frag%coef_ref_ready
      flush(6)
    end if
    call capture_occmap_pair_snapshot(dg_frag, itt_tag, 'density_reconstruct_post')

  contains

    subroutine prepare_grid_buffers_owner_map(i_local_grid, igrid0_grid, npt_blk_grid, nxyz_grid, use_subgroup_slot)
      implicit none
      integer, intent(in) :: i_local_grid, igrid0_grid, npt_blk_grid
      integer, intent(in) :: nxyz_grid(3)
      logical, intent(in) :: use_subgroup_slot
      type(density_grid_point_info) :: point

      if (use_subgroup_slot .and. dg_frag%parallel_mode_orbital) then
        write(*,'(1x,a,i0,a,i0,a)') '[FATAL] orbital mode entered subgroup-slot path: rank=', dg_frag%id, &
          ' id_frag=', dg_frag%id_frag, ' path=prepare_grid_buffers_owner_map'
        flush(6)
        stop 'DG-Fragment RT orbital mode: subgroup real-space slot path is not allowed'
      end if

!$omp parallel do private(igrid, point) schedule(static)
      do igrid = 1, npt_blk_grid
        point = dg_frag%density_grid_points(igrid0_grid + igrid - 1, i_local_grid)
        ix_buf(igrid) = point%ix
        iy_buf(igrid) = point%iy
        iz_buf(igrid) = point%iz
        if (.not. point%is_primary) then
          ixg_buf(igrid) = 1
          iyg_buf(igrid) = 1
          izg_buf(igrid) = 1
          owner_buf(igrid) = -1
          slot_buf(igrid) = 0
          cycle
        end if
        ixg_buf(igrid) = point%ixg
        iyg_buf(igrid) = point%iyg
        izg_buf(igrid) = point%izg
        owner_buf(igrid) = point%owner_rank
        if (use_subgroup_slot) then
          slot_buf(igrid) = point%subgroup_send_slot
        else
          slot_buf(igrid) = point%send_slot
        end if
      end do
!$omp end parallel do
    end subroutine prepare_grid_buffers_owner_map

    subroutine prepare_grid_buffers_owner_map_no_slot(i_local_grid, igrid0_grid, npt_blk_grid, nxyz_grid)
      implicit none
      integer, intent(in) :: i_local_grid, igrid0_grid, npt_blk_grid
      integer, intent(in) :: nxyz_grid(3)
      type(density_grid_point_info) :: point

!$omp parallel do private(igrid, point) schedule(static)
      do igrid = 1, npt_blk_grid
        point = dg_frag%density_grid_points(igrid0_grid + igrid - 1, i_local_grid)
        ix_buf(igrid) = point%ix
        iy_buf(igrid) = point%iy
        iz_buf(igrid) = point%iz
        if (.not. point%is_primary) then
          ixg_buf(igrid) = 1
          iyg_buf(igrid) = 1
          izg_buf(igrid) = 1
          owner_buf(igrid) = -1
          slot_buf(igrid) = 0
          cycle
        end if
        ixg_buf(igrid) = point%ixg
        iyg_buf(igrid) = point%iyg
        izg_buf(igrid) = point%izg
        owner_buf(igrid) = point%owner_rank
        slot_buf(igrid) = 0
      end do
!$omp end parallel do
    end subroutine prepare_grid_buffers_owner_map_no_slot

    integer function find_grid_owner(ixg, iyg, izg) result(owner)
      implicit none
    integer, intent(in) :: ixg, iyg, izg
    integer :: jrank

    owner = dg_frag%id
    do jrank = 0, dg_frag%isize - 1
      if (ixg < mg%is_all(1, jrank) .or. ixg > mg%ie_all(1, jrank)) cycle
      if (iyg < mg%is_all(2, jrank) .or. iyg > mg%ie_all(2, jrank)) cycle
      if (izg < mg%is_all(3, jrank) .or. izg > mg%ie_all(3, jrank)) cycle
      owner = jrank
      return
    end do
  end function find_grid_owner

  logical function target_rank_owned_by_handler(target_rank) result(is_handler)
    implicit none
    integer, intent(in) :: target_rank
    integer :: handler_id_frag
    integer :: frag_size

    frag_size = max(1, dg_frag%isize_frag)
    handler_id_frag = modulo(target_rank, frag_size)
    is_handler = (dg_frag%id_frag == handler_id_frag)
  end function target_rank_owned_by_handler

  logical function local_fragments_overlap_rank_box(dg_frag_local, mg_local, rank_box) result(overlap)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag_local
    type(s_rgrid), intent(in) :: mg_local
    integer, intent(in) :: rank_box
    integer :: jfrag

    overlap = .false.
    do jfrag = dg_frag_local%ifrag_start, dg_frag_local%ifrag_end
      if (fragment_overlaps_box(dg_frag_local%ixyz_frag(:, jfrag), dg_frag_local%nxyz_domain(:, jfrag), &
          mg_local%is_all(:, rank_box), mg_local%ie_all(:, rank_box), mg_local%num)) then
        overlap = .true.
        return
      end if
    end do
  end function local_fragments_overlap_rank_box

  logical function rank_fragments_overlap_local_box(dg_frag_local, mg_local, source_rank) result(overlap)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag_local
    type(s_rgrid), intent(in) :: mg_local
    integer, intent(in) :: source_rank
    integer :: jfrag

    overlap = .false.
    do jfrag = 1, dg_frag_local%n_frag
      if (dg_frag_local%id_array(jfrag) /= source_rank) cycle
      if (fragment_overlaps_box(dg_frag_local%ixyz_frag(:, jfrag), dg_frag_local%nxyz_domain(:, jfrag), &
          mg_local%is, mg_local%ie, mg_local%num)) then
        overlap = .true.
        return
      end if
    end do
  end function rank_fragments_overlap_local_box

  logical function fragment_overlaps_box(ixyz_frag, nxyz_domain, box_is, box_ie, num_grid) result(overlap)
    implicit none
    integer, intent(in) :: ixyz_frag(3), nxyz_domain(3), box_is(3), box_ie(3), num_grid(3)
    integer :: idir

    overlap = .true.
    do idir = 1, 3
      if (.not. periodic_range_overlaps(ixyz_frag(idir), nxyz_domain(idir), box_is(idir), box_ie(idir), num_grid(idir))) then
        overlap = .false.
        return
      end if
    end do
  end function fragment_overlaps_box

  logical function periodic_range_overlaps(start_idx, length_idx, box_is, box_ie, num_grid) result(overlap)
    implicit none
    integer, intent(in) :: start_idx, length_idx, box_is, box_ie, num_grid
    integer :: current_start, remaining_len, segment_len, seg_is, seg_ie

    overlap = .false.
    if (length_idx <= 0) return
    current_start = mod(start_idx - 1, num_grid) + 1
    remaining_len = length_idx
    do while (remaining_len > 0)
      segment_len = min(remaining_len, num_grid - current_start + 1)
      seg_is = current_start
      seg_ie = current_start + segment_len - 1
      if (max(seg_is, box_is) <= min(seg_ie, box_ie)) then
        overlap = .true.
        return
      end if
      remaining_len = remaining_len - segment_len
      current_start = 1
    end do
  end function periodic_range_overlaps

  integer function map_global_to_phi_box_coord_ham(ig, lb, ub, lgtot) result(iloc)
    implicit none
    integer, intent(in) :: ig, lb, ub, lgtot

    iloc = modulo(ig - 1, lgtot) + 1
    if (iloc < lb) then
      iloc = iloc + ((lb - iloc + lgtot - 1) / lgtot) * lgtot
    end if
    if (iloc > ub) then
      iloc = iloc - ((iloc - ub + lgtot - 1) / lgtot) * lgtot
    end if
    if (iloc < lb .or. iloc > ub) iloc = 0
  end function map_global_to_phi_box_coord_ham

  end subroutine calculate_density_from_fragments
