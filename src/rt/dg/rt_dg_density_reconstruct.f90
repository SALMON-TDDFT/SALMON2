  subroutine calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s, itt_debug)
    use structures
    use salmon_global, only: nelec, nelec_spin
    use communication, only: comm_summation, comm_bcast, comm_alltoallv, comm_send, comm_recv, COMM_GROUP_NULL
    use rt_dg_fragment_ops, only: refresh_pw_coef_cache, gather_full_coef_view, copy_overlap_operator_to_dense
    use rt_dg_fragment_types, only: density_grid_point_info
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_scalar),         intent(inout) :: rho
    type(s_scalar),         intent(inout) :: rho_s(system%nspin)
    integer, intent(in), optional :: itt_debug

    integer :: ifrag, io, i_local, ispin
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
    integer :: nblocks_ifrag, first_block_offset, block_step_blocks, block_offset
    integer :: valid_basis_count
    integer :: handler_id_frag
    integer :: spin_offset
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
    character(96) :: env_point_probe
    integer :: env_status
    integer :: env_probe_status
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
    real(8) :: psi2_occ_val, psi2_dv_val
    real(8) :: psi2_after_partition_val, psi2_after_slot_val, psi2_after_any_norm_val
    integer :: state_rhobf_trace_io, state_rhobf_trace_spin

    itt_tag = -1
    if (present(itt_debug)) itt_tag = itt_debug

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
    ! Density reconstruction uses subgroup-distributed projection and collective reductions on icomm_frag.
    subgroup_root_rank = dg_frag%id - dg_frag%id_frag
    total_send_pts = 0
    do irank = 0, dg_frag%isize - 1
      npts = dg_frag%density_send_count(irank)
      if (npts <= 0) cycle
      allocate(rho_send(irank)%f((system%nspin + 1) * npts, 1, 1))
      rho_send(irank)%f(:, :, :) = 0.0d0
    end do
    do irank = 0, dg_frag%isize - 1
      npts = dg_frag%density_recv_map(irank)%npts
      if (npts <= 0) cycle
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
            occ_factor * matmul(dg_frag%coef_mix(1:n_basis_mix, io:io, ispin), &
                                conjg(transpose(dg_frag%coef_mix(1:n_basis_mix, io:io, ispin))))
        end do
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
            call prepare_grid_buffers_owner_map(i_local, igrid0, npt_blk, nxyz, n_pw == 0)
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
                theta = kpw_hx(ipw) * real(ixg - 1, 8) + kpw_hy(ipw) * real(iyg - 1, 8) + kpw_hz(ipw) * real(izg - 1, 8)
                phase_cache(igrid, ipw) = cmplx(cos(theta), sin(theta), kind=8) * inv_sqrt_vol
              end do
            end do
!$omp end parallel do
          end if

          do ispin = 1, system%nspin
            nocc_spin = dg_frag%nocc_spin(ispin)
            nbf = dg_frag%n_basis(ifrag, ispin)
            if (nbf <= 0 .or. nocc_spin <= 0) cycle
            valid_basis_count = valid_basis_count_spin(ispin)

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
              if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
                call cpu_time(t_rho1)
                call comm_summation(rho_blk(1:npt_blk), rho_blk_reduced(1:npt_blk), npt_blk, dg_frag%icomm_frag)
                rho_blk(1:npt_blk) = rho_blk_reduced(1:npt_blk)
                call cpu_time(t_setup1)
                time_project_rho_reduce = time_project_rho_reduce + (t_setup1 - t_rho1)
              end if


              if (dg_frag%is_frag_root) then
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

                do io0 = io_s_frag, io_e_frag, state_block_size
                  nbatch = min(state_block_size, io_e_frag - io0 + 1)

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
                        rho_accum = rho_accum + psi_blk_re(igrid, io) * psi_blk_re(igrid, io) + &
                                    psi_blk_im(igrid, io) * psi_blk_im(igrid, io)
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
                    if (.not. dg_frag%is_frag_root) cycle
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
!$omp parallel do private(slot, ixg, iyg, izg, bx, by, bz, ispin, spin_offset, rho_contrib, rho_raw_contrib, imported_unpack_weight) schedule(static)
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
    if (allocated(n_basis_mix_spin)) deallocate(n_basis_mix_spin)
    if (allocated(g211_cos_x)) deallocate(g211_cos_x, g211_sin_x)
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

  contains

    subroutine prepare_grid_buffers_owner_map(i_local_grid, igrid0_grid, npt_blk_grid, nxyz_grid, use_subgroup_slot)
      implicit none
      integer, intent(in) :: i_local_grid, igrid0_grid, npt_blk_grid
      integer, intent(in) :: nxyz_grid(3)
      logical, intent(in) :: use_subgroup_slot
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

    handler_id_frag = modulo(target_rank, dg_frag%isize_frag)
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
