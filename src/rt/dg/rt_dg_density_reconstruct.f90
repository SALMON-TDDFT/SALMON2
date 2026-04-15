  subroutine calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s)
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

    integer :: ifrag, io, jo, i_local, ispin
    integer :: ifrag_basis
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
    real(8) :: total_charge, total_charge_local, target_charge, rocc_total, scale_rho, rho_sum_local
    real(8) :: occ_factor
    real(8) :: phi_sample_probe
    real(8) :: phi_cache_global_probe, phi_cache_local_probe
    real(8) :: boxL(3), inv_sqrt_vol, theta, inv_lgnum1
    real(8) :: t_total0, t_total1, t_cache0, t_cache1
    real(8) :: t_project0, t_project1, t_comm0, t_comm1, t_norm0, t_norm1
    real(8) :: t_setup0, t_setup1, t_psi0, t_psi1, t_rho0, t_rho1
    real(8) :: t_copy0, t_copy1, t_cleanup0, t_cleanup1
    real(8) :: t_post_split0, t_post_split1, t_post_dealloc0, t_post_dealloc1
    real(8) :: t_post_owner0, t_post_owner1, t_post_remote0, t_post_remote1
    real(8) :: time_cache, time_project, time_comm, time_norm
    real(8) :: time_copy, time_cleanup
    real(8) :: time_unaccounted
    real(8) :: time_preproject, time_postproject_precomm, time_postcomm_prenorm
    real(8) :: time_post_split_pw, time_post_dealloc, time_post_precomm_other
    real(8) :: time_post_owner_accum, time_post_remote_pack, time_post_precomm_residual
    real(8) :: time_comm_subgroup_reduce, time_comm_pack, time_comm_exchange, time_comm_unpack
    real(8) :: time_project_setup, time_project_psi, time_project_rho
    real(8) :: time_project_grid_prep, time_project_phi_pack, time_project_overhead
    real(8) :: time_project_dmat_build
    real(8) :: t_dmat0, t_dmat1
    real(8) :: time_cache_pw_refresh, time_cache_phi_block_refresh
    logical :: use_mixed_density
    logical :: split_pw_density
    logical :: rebuilt_pw_cache, rebuilt_phi_block_cache
    logical :: need_pw_cache_alloc, need_pw_cache_expand
    logical :: need_phi_cache_alloc, need_phi_count_alloc, need_phi_cache_invalid, need_phi_cache_resize
    integer, allocatable :: ix_buf(:), iy_buf(:), iz_buf(:), owner_buf(:), ixg_buf(:), iyg_buf(:), izg_buf(:)
    integer, allocatable :: slot_buf(:), local_grid_ids(:), remote_grid_ids(:), valid_remote_grid_ids(:)
    integer, allocatable :: basis_gid(:), valid_basis_ids(:)
    integer, allocatable :: basis_gid_spin(:,:), valid_basis_ids_spin(:,:), valid_basis_count_spin(:)
    integer, allocatable :: phi_map_x(:), phi_map_y(:), phi_map_z(:)
    real(8), allocatable :: basis_sdiag_probe(:,:,:)
    real(8), allocatable :: phi_col_metric_total(:,:,:)
    real(8), allocatable :: phi_col_full_metric_total(:,:,:)
    real(8), allocatable :: basis_smat_probe(:,:,:,:)
    real(8), allocatable :: phi_gram_total(:,:,:,:)
    real(8), allocatable :: phi_frag_metric_total(:,:,:,:)
    real(8), allocatable :: basis_frag_metric_total(:,:,:,:)
    type(s_scalar), allocatable :: rho_send(:), rho_recv(:)
    integer, allocatable :: send_counts(:), recv_counts(:), send_displs(:), recv_displs(:)
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
    real(8), allocatable :: rho_blk_cross_partial(:)
    real(8), allocatable :: occ_cache(:), occ_sqrt_cache(:)
    complex(8), allocatable :: coef_c_full(:,:), coef_c_frag(:,:)
    complex(8), allocatable :: coef_full_cache(:,:,:), coef_pw_view_cache(:,:,:)
    real(8) :: time_project_rho_reduce, time_project_phi_block_build
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
    complex(8), allocatable :: overlap_probe_trace(:,:)
    integer, allocatable :: subgroup_self_ixg_tmp(:), subgroup_self_iyg_tmp(:), subgroup_self_izg_tmp(:)
    logical :: enable_density_reconstruct_trace
    logical :: enable_density_split_charge_probe
    logical :: enable_density_renormalize
    integer :: env_len, env_status
    character(len=64) :: env_val
    real(8) :: coef_norm_probe, rho_probe_charge, phi_norm_probe, psi_norm_probe
    real(8) :: phi_col_norm_probe(3), phi_col_hvol_probe(3), phi_sdiag_probe(3)
    real(8) :: coef_map_local_probe(3,3), coef_map_global_probe(3,3), coef_map_diff_probe(3,3)
    real(8) :: orbital_norm_probe_local(3), orbital_norm_probe_total(3)
    real(8) :: orbital_norm_frag_local(2,3), orbital_norm_frag_total(2,3)
    real(8) :: overlap_state_probe(3), overlap_elec_probe, coef_state_probe(3), overlap_diag_probe(3)
    real(8) :: overlap_self_probe(2,3), overlap_cross_probe(3)
    real(8) :: overlap_ff_probe(3), overlap_fp_probe(3), overlap_pp_probe(3)
    real(8) :: overlap_ff_total_probe, overlap_fp_total_probe, overlap_pp_total_probe
    real(8) :: frag_trace_probe, frag_state_trace_probe(3)
    real(8) :: frag_trace_probe_local_s
    real(8) :: frag_trace_sumspin
    real(8) :: frag_charge_sumspin
    real(8) :: frag_coeff_occ_norm_sumspin
    real(8) :: frag_d_diag_sumspin
    real(8) :: frag_owner_charge_sumspin
    real(8) :: frag_remote_send_sumspin
    real(8) :: frag_state_real_probe(3)
    real(8) :: phi_col_full_local(3), phi_col_full_global(3)
    real(8) :: block_charge_direct, block_charge_psi
    real(8) :: rho_direct_probe, coeff_re_probe, coeff_im_probe
    real(8) :: charge_root_tmp_global, charge_root_sum_global
    real(8) :: charge_weighted_total_global, charge_weighted_total_pre_norm
    real(8) :: charge_fragment_local, charge_fragment_global
    real(8) :: charge_pw_local, charge_pw_global
    real(8) :: charge_cross_local, charge_cross_global
    real(8) :: charge_fragment_preowner_local, charge_fragment_preowner_global
    real(8) :: charge_fragment_postowner_local, charge_fragment_postowner_global
    real(8) :: coeff_occ_norm_probe, d_diag_probe
    real(8) :: charge_blk_all, charge_blk_owner, charge_blk_handler, charge_blk_slot0
    real(8) :: g211_blk_all_re, g211_blk_all_im, g211_blk_owner_re, g211_blk_owner_im
    real(8) :: g211_blk_handler_re, g211_blk_handler_im
    real(8) :: g211_root_sum_re_local, g211_root_sum_im_local
    real(8) :: g211_rank_buf_re_local, g211_rank_buf_im_local
    real(8) :: g211_rank_buf_re_global, g211_rank_buf_im_global
    real(8) :: g211_pre_total_re, g211_pre_total_im
    real(8) :: phase_theta, phase_re, phase_im, g211_pred_re, g211_pred_im
    real(8) :: g211_re_line, g211_im_line, rho_val
    real(8) :: charge_mismatch_abs, charge_mismatch_tol
    real(8), allocatable :: g211_cos_x(:), g211_sin_x(:)
    real(8), allocatable :: kpw_hx(:), kpw_hy(:), kpw_hz(:)

    enable_density_reconstruct_trace = .false.
    enable_density_split_charge_probe = .false.
    enable_density_renormalize = .false.
    call get_environment_variable("SALMON_DG_DENSITY_TRACE", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        enable_density_reconstruct_trace = .true.
      end if
    end if
    call get_environment_variable("SALMON_DG_DENSITY_SPLIT_CHARGE_TRACE", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        enable_density_split_charge_probe = .true.
      end if
    end if
    call get_environment_variable("SALMON_DG_ALLOW_DENSITY_RENORMALIZE", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        enable_density_renormalize = .true.
      end if
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
    time_copy = 0.0d0
    time_postcomm_prenorm = 0.0d0
    time_post_split_pw = 0.0d0
    time_post_dealloc = 0.0d0
    time_post_precomm_other = 0.0d0
    time_post_owner_accum = 0.0d0
    time_post_remote_pack = 0.0d0
    time_post_precomm_residual = 0.0d0
    time_comm_subgroup_reduce = 0.0d0
    time_comm_pack = 0.0d0
    time_comm_exchange = 0.0d0
    time_comm_unpack = 0.0d0
    time_project_setup = 0.0d0
    time_project_phi_pack = 0.0d0
    time_project_overhead = 0.0d0
    time_project_dmat_build = 0.0d0
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
    overlap_ff_probe(:) = 0.0d0
    overlap_fp_probe(:) = 0.0d0
    overlap_pp_probe(:) = 0.0d0
    overlap_ff_total_probe = 0.0d0
    overlap_fp_total_probe = 0.0d0
    overlap_pp_total_probe = 0.0d0
    frag_trace_probe = 0.0d0
    frag_state_trace_probe(:) = 0.0d0
    charge_root_tmp_global = 0.0d0
    charge_root_sum_global = 0.0d0
    charge_weighted_total_global = 0.0d0
    charge_weighted_total_pre_norm = 0.0d0
    charge_fragment_local = 0.0d0
    charge_fragment_global = 0.0d0
    charge_pw_local = 0.0d0
    charge_pw_global = 0.0d0
    charge_cross_local = 0.0d0
    charge_cross_global = 0.0d0
    charge_fragment_preowner_local = 0.0d0
    charge_fragment_preowner_global = 0.0d0
    charge_fragment_postowner_local = 0.0d0
    charge_fragment_postowner_global = 0.0d0
    rebuilt_pw_cache = .false.
    rebuilt_phi_block_cache = .false.
    need_pw_cache_alloc = .false.
    need_pw_cache_expand = .false.
    need_phi_cache_alloc = .false.
    need_phi_count_alloc = .false.
    need_phi_cache_invalid = .false.
    need_phi_cache_resize = .false.
    if (enable_density_reconstruct_trace .and. dg_frag%is_frag_root) then
      write(*,'(1x,a,i0,a,i0,a,i0)') "        density trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " stage=entry_root root_world=", dg_frag%id - dg_frag%id_frag
      write(*,'(1x,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') "        density mg-range: rank=", dg_frag%id, &
        " mg_is=", mg%is(1), mg%is(2), mg%is(3), " mg_ie=", mg%ie(1), mg%ie(2), mg%ie(3), &
        " mg_is_all=", mg%is_all(1, dg_frag%id), mg%is_all(2, dg_frag%id), mg%is_all(3, dg_frag%id), &
        " mg_ie_all=", mg%ie_all(1, dg_frag%id), mg%ie_all(2, dg_frag%id), mg%ie_all(3, dg_frag%id)
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
      allocate(phi_col_full_metric_total(3, system%nspin, max(1, ifrag_count)))
      allocate(basis_smat_probe(3, 3, system%nspin, max(1, ifrag_count)))
      allocate(phi_gram_total(3, 3, system%nspin, max(1, ifrag_count)))
      allocate(phi_frag_metric_total(nbf_max, nbf_max, system%nspin, max(1, ifrag_count)))
      allocate(basis_frag_metric_total(nbf_max, nbf_max, system%nspin, max(1, ifrag_count)))
      basis_sdiag_probe(:, :, :) = 0.0d0
      phi_col_metric_total(:, :, :) = 0.0d0
      phi_col_full_metric_total(:, :, :) = 0.0d0
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
    split_pw_density = (n_pw > 0)
    use_mixed_density = (n_pw > 0 .and. .not. split_pw_density .and. dg_frag%mixed_basis_ready .and. allocated(dg_frag%mixed_transform) .and. &
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
          if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
            if (io <= size(dg_frag%occ_state, 1) .and. ispin <= size(dg_frag%occ_state, 2)) then
              occ_factor = max(0.0d0, dg_frag%occ_state(io, ispin))
            end if
          else if (allocated(system%rocc)) then
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
    boxL(1) = dg_frag%hgs(1) * real(dg_frag%lgnum_total(1), 8)
    boxL(2) = dg_frag%hgs(2) * real(dg_frag%lgnum_total(2), 8)
    boxL(3) = dg_frag%hgs(3) * real(dg_frag%lgnum_total(3), 8)
    inv_sqrt_vol = 1.0d0 / sqrt(max(1.0d-16, boxL(1) * boxL(2) * boxL(3)))

    call cpu_time(t_setup0)
    if (.not. allocated(dg_frag%density_block_nblocks)) then
      allocate(dg_frag%density_block_nblocks(ifrag_count), dg_frag%density_block_first_offset(ifrag_count), &
               dg_frag%density_block_step(ifrag_count))
      block_idx_global = 0
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        frag_trace_sumspin = 0.0d0
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
      call refresh_pw_coef_cache(dg_frag, nocc_cache)
      rebuilt_pw_cache = (need_pw_cache_alloc .or. need_pw_cache_expand)
      call cpu_time(t_cache1)
      time_cache = time_cache + (t_cache1 - t_cache0)
      time_cache_pw_refresh = time_cache_pw_refresh + (t_cache1 - t_cache0)
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
    if (enable_density_reconstruct_trace .and. n_pw == 0) then
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
        if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
          if (io <= size(dg_frag%occ_state, 1) .and. 1 <= size(dg_frag%occ_state, 2)) then
            occ_factor = max(0.0d0, dg_frag%occ_state(io, 1))
          end if
        else if (allocated(system%rocc)) then
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
    if (enable_density_reconstruct_trace .and. split_pw_density) then
      n_tot = dg_frag%n_mat_max + n_pw
      allocate(overlap_probe(n_tot, n_tot), overlap_vec(n_tot))
      overlap_probe(:, :) = (0.0d0, 0.0d0)
      overlap_vec(:) = (0.0d0, 0.0d0)
      call gather_full_coef_view(dg_frag, 1, dg_frag%n_mat_max, nocc_cache, coef_probe_full, coef_probe_pw)
      call copy_overlap_operator_to_dense(dg_frag, 1, .true., overlap_probe)
      do io = 1, min(3, dg_frag%nocc_spin(1))
        overlap_ff_probe(io) = real(sum(conjg(coef_probe_full(:, io)) * matmul(overlap_probe(1:dg_frag%n_mat_max, 1:dg_frag%n_mat_max), &
          coef_probe_full(:, io))), kind=8)
        if (n_pw > 0) then
          overlap_fp_probe(io) = 2.0d0 * real(sum(conjg(coef_probe_full(:, io)) * matmul(overlap_probe(1:dg_frag%n_mat_max, dg_frag%n_mat_max+1:n_tot), &
            coef_probe_pw(:, io))), kind=8)
          overlap_pp_probe(io) = real(sum(conjg(coef_probe_pw(:, io)) * matmul(overlap_probe(dg_frag%n_mat_max+1:n_tot, dg_frag%n_mat_max+1:n_tot), &
            coef_probe_pw(:, io))), kind=8)
        end if
      end do
      do io = 1, dg_frag%nocc_spin(1)
        occ_factor = 1.0d0
        if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
          if (io <= size(dg_frag%occ_state, 1) .and. 1 <= size(dg_frag%occ_state, 2)) occ_factor = max(0.0d0, dg_frag%occ_state(io, 1))
        else if (allocated(system%rocc)) then
          if (io <= size(system%rocc, 1) .and. 1 <= size(system%rocc, 3)) occ_factor = max(0.0d0, system%rocc(io, 1, 1))
        end if
        if (occ_factor <= 0.0d0) cycle
        overlap_ff_total_probe = overlap_ff_total_probe + occ_factor * real(sum(conjg(coef_probe_full(:, io)) * &
          matmul(overlap_probe(1:dg_frag%n_mat_max, 1:dg_frag%n_mat_max), coef_probe_full(:, io))), kind=8)
        if (n_pw > 0) then
          overlap_fp_total_probe = overlap_fp_total_probe + occ_factor * 2.0d0 * real(sum(conjg(coef_probe_full(:, io)) * &
            matmul(overlap_probe(1:dg_frag%n_mat_max, dg_frag%n_mat_max+1:n_tot), coef_probe_pw(:, io))), kind=8)
          overlap_pp_total_probe = overlap_pp_total_probe + occ_factor * real(sum(conjg(coef_probe_pw(:, io)) * &
            matmul(overlap_probe(dg_frag%n_mat_max+1:n_tot, dg_frag%n_mat_max+1:n_tot), coef_probe_pw(:, io))), kind=8)
        end if
      end do
      if (dg_frag%id == 0) then
        write(*,'(1x,a,3(1pe12.4,1x),a,3(1pe12.4,1x),a,3(1pe12.4,1x))') &
          "        density overlap FFPP probe: ff=", overlap_ff_probe(1), overlap_ff_probe(2), overlap_ff_probe(3), &
          " fp=", overlap_fp_probe(1), overlap_fp_probe(2), overlap_fp_probe(3), &
          " pp=", overlap_pp_probe(1), overlap_pp_probe(2), overlap_pp_probe(3)
        write(*,'(1x,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
          "        density overlap FFPP total: ff=", overlap_ff_total_probe, " fp=", overlap_fp_total_probe, &
          " pp=", overlap_pp_total_probe, " total=", overlap_ff_total_probe + overlap_fp_total_probe + overlap_pp_total_probe
        flush(6)
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
    allocate(coef_full_cache(max(1, dg_frag%n_mat_max), max(1, nocc_cache), system%nspin))
    coef_full_cache(:, :, :) = (0.0d0, 0.0d0)
    if (n_pw > 0) then
      allocate(coef_pw_view_cache(max(1, n_pw), max(1, nocc_cache), system%nspin))
      coef_pw_view_cache(:, :, :) = (0.0d0, 0.0d0)
    end if
    do ispin = 1, system%nspin
      nocc_spin = dg_frag%nocc_spin(ispin)
      if (nocc_spin <= 0) cycle
      call gather_full_coef_view(dg_frag, ispin, dg_frag%n_mat_max, nocc_spin, coef_probe_full, coef_probe_pw, 1, nocc_spin)
      coef_full_cache(1:dg_frag%n_mat_max, 1:nocc_spin, ispin) = coef_probe_full(1:dg_frag%n_mat_max, 1:nocc_spin)
      if (n_pw > 0 .and. allocated(coef_pw_view_cache) .and. allocated(coef_probe_pw)) then
        coef_pw_view_cache(1:n_pw, 1:nocc_spin, ispin) = coef_probe_pw(1:n_pw, 1:nocc_spin)
      end if
    end do
    if (allocated(coef_probe_full)) deallocate(coef_probe_full)
    if (allocated(coef_probe_pw)) deallocate(coef_probe_pw)
    if (enable_density_reconstruct_trace) then
      allocate(overlap_probe_trace(max(1, dg_frag%n_mat_max), max(1, dg_frag%n_mat_max)))
      overlap_probe_trace(:, :) = (0.0d0, 0.0d0)
      call copy_overlap_operator_to_dense(dg_frag, 1, .true., overlap_probe_trace)
    end if
    call cpu_time(t_setup0)
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
      phi_lb1 = lbound(dg_frag%phi_frag, 1)
      phi_lb2 = lbound(dg_frag%phi_frag, 2)
      phi_lb3 = lbound(dg_frag%phi_frag, 3)
      phi_ub1 = ubound(dg_frag%phi_frag, 1)
      phi_ub2 = ubound(dg_frag%phi_frag, 2)
      phi_ub3 = ubound(dg_frag%phi_frag, 3)
      phi_lg1 = dg_frag%lgnum_total(1)
      phi_lg2 = dg_frag%lgnum_total(2)
      phi_lg3 = dg_frag%lgnum_total(3)
      allocate(phi_map_x(phi_lg1), phi_map_y(phi_lg2), phi_map_z(phi_lg3))
      do ixg = 1, phi_lg1
        phi_map_x(ixg) = map_global_to_phi_box_coord_ham(ixg, phi_lb1, phi_ub1, phi_lg1)
      end do
      do iyg = 1, phi_lg2
        phi_map_y(iyg) = map_global_to_phi_box_coord_ham(iyg, phi_lb2, phi_ub2, phi_lg2)
      end do
      do izg = 1, phi_lg3
        phi_map_z(izg) = map_global_to_phi_box_coord_ham(izg, phi_lb3, phi_ub3, phi_lg3)
      end do
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
            bx = phi_map_x(modulo(ixg - 1, phi_lg1) + 1)
            by = phi_map_y(modulo(iyg - 1, phi_lg2) + 1)
            bz = phi_map_z(modulo(izg - 1, phi_lg3) + 1)
            if (bx == 0 .or. by == 0 .or. bz == 0) cycle
            do istate_frag = 1, dg_frag%nstate_frag
              dg_frag%density_phi_block_cache(igrid, istate_frag, block_cache_idx, i_local) = &
                dg_frag%phi_frag(bx, by, bz, istate_frag, i_local)
            end do
          end do
!$omp end parallel do
          if (enable_density_reconstruct_trace .and. block_cache_idx == 1 .and. npt_cache > 0) then
            ixg = dg_frag%density_grid_points(igrid0, i_local)%ixg
            iyg = dg_frag%density_grid_points(igrid0, i_local)%iyg
            izg = dg_frag%density_grid_points(igrid0, i_local)%izg
            ix = dg_frag%density_grid_points(igrid0, i_local)%ix
            iy = dg_frag%density_grid_points(igrid0, i_local)%iy
            iz = dg_frag%density_grid_points(igrid0, i_local)%iz
            bx = phi_map_x(modulo(ixg - 1, phi_lg1) + 1)
            by = phi_map_y(modulo(iyg - 1, phi_lg2) + 1)
            bz = phi_map_z(modulo(izg - 1, phi_lg3) + 1)
            phi_cache_global_probe = 0.0d0
            phi_cache_local_probe = 0.0d0
            if (bx /= 0 .and. by /= 0 .and. bz /= 0) then
              phi_cache_global_probe = dg_frag%phi_frag(bx, by, bz, 1, i_local)
            end if
            if (ix >= lbound(dg_frag%phi_frag, 1) .and. ix <= ubound(dg_frag%phi_frag, 1) .and. &
                iy >= lbound(dg_frag%phi_frag, 2) .and. iy <= ubound(dg_frag%phi_frag, 2) .and. &
                iz >= lbound(dg_frag%phi_frag, 3) .and. iz <= ubound(dg_frag%phi_frag, 3)) then
              phi_cache_local_probe = dg_frag%phi_frag(ix, iy, iz, 1, i_local)
            end if
            write(*,'(1x,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x),a,1pe12.4,a,1pe12.4)') &
              "        density cache coord probe: rank=", dg_frag%id, " ifrag=", ifrag, &
              " ixiyiz=", ix, iy, iz, " ixgiygizg=", ixg, iyg, izg, " mapped=", bx, by, bz, &
              " phi(globalmap)=", phi_cache_global_probe, " phi(local)=", phi_cache_local_probe
            flush(6)
          end if
        end do
      end do
      deallocate(phi_map_x, phi_map_y, phi_map_z)
      dg_frag%density_phi_block_size = grid_block_size
      dg_frag%density_phi_block_cache_valid = .true.
      rebuilt_phi_block_cache = .true.
    end if
    call cpu_time(t_setup1)
    time_project_phi_block_build = time_project_phi_block_build + (t_setup1 - t_setup0)
    if (rebuilt_phi_block_cache) time_cache_phi_block_refresh = time_cache_phi_block_refresh + (t_setup1 - t_setup0)
    if (enable_density_reconstruct_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a,l1,a,l1,a,l1,a,l1,a,l1,a,1pe12.4)') "        density cache: phi_block rebuilt=", &
        rebuilt_phi_block_cache, " alloc=", need_phi_cache_alloc, " count_alloc=", need_phi_count_alloc, &
        " invalid=", need_phi_cache_invalid, " resize=", need_phi_cache_resize, " dt=", time_cache_phi_block_refresh
      flush(6)
    end if
    if (enable_density_reconstruct_trace .and. n_pw == 0) then
      do i_local = 1, ifrag_count
        ifrag = dg_frag%ifrag_start + i_local - 1
        do ispin = 1, system%nspin
          nprobe_cols = min(3, dg_frag%n_basis(ifrag, ispin))
          if (nprobe_cols <= 0) cycle
          phi_col_full_local(:) = 0.0d0
          do iz = max(dg_frag%frag_core_lo(3, ifrag), mg%is(3)), min(dg_frag%frag_core_hi(3, ifrag), mg%ie(3))
            bz = map_global_to_phi_box_coord_ham(iz, phi_lb3, phi_ub3, phi_lg3)
            if (bz == 0) cycle
            do iy = max(dg_frag%frag_core_lo(2, ifrag), mg%is(2)), min(dg_frag%frag_core_hi(2, ifrag), mg%ie(2))
              by = map_global_to_phi_box_coord_ham(iy, phi_lb2, phi_ub2, phi_lg2)
              if (by == 0) cycle
              do ix = max(dg_frag%frag_core_lo(1, ifrag), mg%is(1)), min(dg_frag%frag_core_hi(1, ifrag), mg%ie(1))
                bx = map_global_to_phi_box_coord_ham(ix, phi_lb1, phi_ub1, phi_lg1)
                if (bx == 0) cycle
                do iprobe = 1, nprobe_cols
                  phi_col_full_local(iprobe) = phi_col_full_local(iprobe) + &
                    dg_frag%phi_frag(bx, by, bz, iprobe, i_local) * dg_frag%phi_frag(bx, by, bz, iprobe, i_local) * system%hvol
                end do
              end do
            end do
          end do
          if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
            call comm_summation(phi_col_full_local(1:nprobe_cols), phi_col_full_global(1:nprobe_cols), nprobe_cols, dg_frag%icomm_frag)
            phi_col_full_metric_total(1:nprobe_cols, ispin, i_local) = phi_col_full_global(1:nprobe_cols)
          else
            phi_col_full_metric_total(1:nprobe_cols, ispin, i_local) = phi_col_full_local(1:nprobe_cols)
          end if
        end do
      end do
    end if
    call cpu_time(t_project0)
      i_local = 0
      block_idx_global = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        ifrag_basis = ifrag
        frag_trace_sumspin = 0.0d0
        frag_charge_sumspin = 0.0d0
        frag_coeff_occ_norm_sumspin = 0.0d0
        frag_d_diag_sumspin = 0.0d0
        frag_owner_charge_sumspin = 0.0d0
        frag_remote_send_sumspin = 0.0d0

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
          nbf = dg_frag%n_basis(ifrag_basis, ispin)
          if (nbf <= 0) cycle
          do istate_frag = 1, nbf
            basis_gid(istate_frag) = dg_frag%index_basis(istate_frag, ifrag_basis, ispin)
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
            nbf = dg_frag%n_basis(ifrag_basis, ispin)
            n_basis_mix = min(dg_frag%mixed_basis_dim(ispin), max_mixed_basis)
            n_basis_mix_spin(ispin) = n_basis_mix
            if (nbf <= 0 .or. n_basis_mix <= 0) cycle
            transform_frag_spin(1:nbf, 1:n_basis_mix, ispin) = (0.0d0, 0.0d0)
            do istate_frag = 1, nbf
              ig_i = dg_frag%index_basis(istate_frag, ifrag_basis, ispin)
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
        if (n_pw == 0 .or. split_pw_density) then
          do ispin = 1, system%nspin
            nbf = dg_frag%n_basis(ifrag_basis, ispin)
            nocc_spin = dg_frag%nocc_spin(ispin)
            if (nbf <= 0 .or. nocc_spin <= 0) cycle
            occ_cache(1:nocc_spin) = 1.0d0
            if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
              do io = 1, nocc_spin
                if (io <= size(dg_frag%occ_state, 1) .and. ispin <= size(dg_frag%occ_state, 2)) then
                  occ_cache(io) = max(0.0d0, dg_frag%occ_state(io, ispin))
                end if
              end do
            else if (allocated(system%rocc)) then
              do io = 1, nocc_spin
                if (io <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
                  occ_cache(io) = max(0.0d0, system%rocc(io, 1, ispin))
                end if
              end do
            end if
            occ_sqrt_cache(1:nocc_spin) = sqrt(occ_cache(1:nocc_spin))
            if (enable_density_reconstruct_trace) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,3(1pe12.4,1x),a,3(1pe12.4,1x),a,l1)') &
                "        density occ probe: rank=", dg_frag%id, " ifrag=", ifrag, " ispin=", ispin, &
                " occ=", occ_cache(1), occ_cache(min(2,nocc_spin)), occ_cache(min(3,nocc_spin)), &
                " rocc=", system%rocc(1,1,ispin), system%rocc(min(2,size(system%rocc,1)),1,ispin), &
                        system%rocc(min(3,size(system%rocc,1)),1,ispin), &
                " buffered=", dg_frag%use_buffered_basis
              flush(6)
            end if
            valid_basis_count = valid_basis_count_spin(ispin)
            call cpu_time(t_dmat0)
            ! Step 3a: build the fragment-local coefficient view from the full
            ! distributed coefficient matrix. Owner-distributed rows alone are
            ! insufficient here because D = C C^H needs cross terms between rows
            ! owned by different subgroup ranks.
            coef_re_frag(1:nbf_max, 1:nocc_spin) = 0.0d0
            coef_im_frag(1:nbf_max, 1:nocc_spin) = 0.0d0
            coef_c_full(1:nbf_max, 1:nocc_spin) = (0.0d0, 0.0d0)
!$omp parallel do private(io, idx_local, istate_frag, ig_i) schedule(static)
            do io = 1, nocc_spin
              do idx_local = 1, valid_basis_count
                istate_frag = valid_basis_ids_spin(idx_local, ispin)
                ig_i = basis_gid_spin(istate_frag, ispin)
                if (ig_i < 1 .or. ig_i > size(coef_full_cache, 1)) cycle
                coef_c_full(istate_frag, io) = coef_full_cache(ig_i, io, ispin)
              end do
            end do
!$omp end parallel do
            coef_re_full(1:nbf_max, 1:nocc_spin, ispin) = real(coef_c_full(1:nbf_max, 1:nocc_spin), kind=8)
            coef_im_full(1:nbf_max, 1:nocc_spin, ispin) = aimag(coef_c_full(1:nbf_max, 1:nocc_spin))
            ! Step 3b: each rank computes weighted C C^H on its state slice
            nocc_per_rank_loc = (nocc_spin + dg_frag%isize_frag - 1) / dg_frag%isize_frag
            io_s_frag = dg_frag%id_frag * nocc_per_rank_loc + 1
            io_e_frag = min((dg_frag%id_frag + 1) * nocc_per_rank_loc, nocc_spin)
            coef_norm_probe = sum(coef_re_full(1:nbf, 1:nocc_spin, ispin)**2 + coef_im_full(1:nbf, 1:nocc_spin, ispin)**2)
            if (enable_density_reconstruct_trace) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,1pe12.4,a,i0,a,i0,a,i0)') "        density frag probe: rank=", dg_frag%id, &
                " ifrag=", ifrag, " ispin=", ispin, " nocc=", nocc_spin, " coef_norm=", coef_norm_probe, &
                " valid_basis=", valid_basis_count, &
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
            end if
            ! Step 3c: AllReduce partial D across icomm_frag
            if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
              call comm_summation(D_partial_re(1:nbf, 1:nbf), D_frag_re(1:nbf, 1:nbf, ispin), &
                                  nbf * nbf, dg_frag%icomm_frag)
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
            coeff_occ_norm_probe = 0.0d0
            do io = 1, nocc_spin
              coeff_occ_norm_probe = coeff_occ_norm_probe + occ_cache(io) * &
                sum(coef_re_full(1:nbf, io, ispin)**2 + coef_im_full(1:nbf, io, ispin)**2)
            end do
            d_diag_probe = 0.0d0
            do io = 1, nbf
              d_diag_probe = d_diag_probe + D_frag_re(io, io, ispin)
            end do
            frag_coeff_occ_norm_sumspin = frag_coeff_occ_norm_sumspin + coeff_occ_norm_probe
            frag_d_diag_sumspin = frag_d_diag_sumspin + d_diag_probe
            if (enable_density_reconstruct_trace .and. dg_frag%is_frag_root) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4)') "        density coef-build check: rank=", dg_frag%id, &
                " ifrag=", ifrag, " ispin=", ispin, " occ_c2_sum=", coeff_occ_norm_probe, &
                " D_diag_sum=", d_diag_probe, " diff=", d_diag_probe - coeff_occ_norm_probe
              flush(6)
            end if
            if (enable_density_reconstruct_trace .and. dg_frag%is_frag_root) then
              frag_trace_probe = 0.0d0
              frag_trace_probe_local_s = 0.0d0
              frag_state_trace_probe(:) = 0.0d0
              do io = 1, nbf
                ig_i = basis_gid_spin(io, ispin)
                if (ig_i <= 0) cycle
                do istate_frag = 1, nbf
                  if (basis_gid_spin(istate_frag, ispin) <= 0) cycle
                  if (allocated(overlap_probe_trace)) then
                    occ_factor = real(overlap_probe_trace(ig_i, basis_gid_spin(istate_frag, ispin)), kind=8)
                  else if (allocated(dg_frag%S_mat)) then
                    occ_factor = dg_frag%S_mat(ig_i, basis_gid_spin(istate_frag, ispin), ispin)
                  else if (allocated(dg_frag%S_mat_c)) then
                    occ_factor = real(dg_frag%S_mat_c(ig_i, basis_gid_spin(istate_frag, ispin), ispin), kind=8)
                  else
                    occ_factor = 0.0d0
                  end if
                  frag_trace_probe = frag_trace_probe + D_frag_re(io, istate_frag, ispin) * occ_factor
                  if (allocated(dg_frag%S_mat)) then
                    occ_factor = dg_frag%S_mat(io, istate_frag, ispin)
                  else if (allocated(dg_frag%S_mat_c)) then
                    occ_factor = real(dg_frag%S_mat_c(io, istate_frag, ispin), kind=8)
                  else
                    occ_factor = 0.0d0
                  end if
                  frag_trace_probe_local_s = frag_trace_probe_local_s + D_frag_re(io, istate_frag, ispin) * occ_factor
                end do
              end do
              do io = 1, min(3, nocc_spin)
                do idx_local = 1, valid_basis_count
                  ig_i = valid_basis_ids_spin(idx_local, ispin)
                  if (ig_i <= 0 .or. basis_gid_spin(ig_i, ispin) <= 0) cycle
                  do istate_frag = 1, valid_basis_count
                    if (valid_basis_ids_spin(istate_frag, ispin) <= 0) cycle
                    if (basis_gid_spin(valid_basis_ids_spin(istate_frag, ispin), ispin) <= 0) cycle
                    if (allocated(overlap_probe_trace)) then
                      occ_factor = real(overlap_probe_trace(basis_gid_spin(ig_i, ispin), &
                                                           basis_gid_spin(valid_basis_ids_spin(istate_frag, ispin), ispin)), kind=8)
                    else if (allocated(dg_frag%S_mat)) then
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
              write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,3(1pe12.4,1x))') "        density frag metric probe: rank=", dg_frag%id, &
                " ifrag=", ifrag, " ispin=", ispin, " Tr(SffDff)=", frag_trace_probe, &
                " Tr(SlocalD)=", frag_trace_probe_local_s, &
                " state_metric=", frag_state_trace_probe(1), frag_state_trace_probe(2), frag_state_trace_probe(3)
              flush(6)
            end if
            frag_trace_sumspin = frag_trace_sumspin + frag_trace_probe
            call cpu_time(t_dmat1)
            time_project_dmat_build = time_project_dmat_build + (t_dmat1 - t_dmat0)
          end do
        else
          ! n_pw > 0: keep the fragment coefficient view identical to the
          ! pure-DG path so we can isolate PW-only effects downstream.
          coef_re_full(1:nbf_max, 1:nocc_cache, 1:system%nspin) = 0.0d0
          coef_im_full(1:nbf_max, 1:nocc_cache, 1:system%nspin) = 0.0d0
          do ispin = 1, system%nspin
            nbf = dg_frag%n_basis(ifrag_basis, ispin)
            nocc_spin = dg_frag%nocc_spin(ispin)
            if (nbf <= 0 .or. nocc_spin <= 0) cycle
            valid_basis_count = valid_basis_count_spin(ispin)
            coef_c_full(1:nbf_max, 1:nocc_spin) = (0.0d0, 0.0d0)
!$omp parallel do private(io, idx_local, istate_frag) schedule(static)
            do io = 1, nocc_spin
              do idx_local = 1, valid_basis_count
                istate_frag = valid_basis_ids_spin(idx_local, ispin)
                if (basis_gid_spin(istate_frag, ispin) < 1 .or. basis_gid_spin(istate_frag, ispin) > size(coef_full_cache, 1)) cycle
                coef_c_full(istate_frag, io) = coef_full_cache(basis_gid_spin(istate_frag, ispin), io, ispin)
              end do
            end do
!$omp end parallel do
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
        if (enable_density_reconstruct_trace .and. dg_frag%is_frag_root) then
          write(*,'(1x,a,i0,a,i0,a,1pe12.4)') "        density frag metric sumspin: rank=", dg_frag%id, &
            " ifrag=", ifrag, " Tr(SffDff)_sumspin=", frag_trace_sumspin
          flush(6)
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
            call prepare_grid_buffers_owner_map(i_local, igrid0, npt_blk, nxyz, .false.)
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
            if (split_pw_density .and. owner_buf(igrid) /= dg_frag%id) cycle
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
            nbf = dg_frag%n_basis(ifrag_basis, ispin)
            if (nbf <= 0 .or. nocc_spin <= 0) cycle
            valid_basis_count = valid_basis_count_spin(ispin)

          call cpu_time(t_setup0)
          phi_blk(1:npt_blk, 1:nbf) = dg_frag%density_phi_block_cache(1:npt_blk, 1:nbf, block_offset + 1, i_local)
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
                    if (owner_buf(igrid) /= dg_frag%id) cycle
                    if (ixg < rho_s_x_lo .or. ixg > rho_s_x_hi .or. &
                      iyg < rho_s_y_lo .or. iyg > rho_s_y_hi .or. &
                      izg < rho_s_z_lo .or. izg > rho_s_z_hi) cycle
                  rho_raw_contrib = rho_blk(igrid)
                  rho_contrib = rho_raw_contrib
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
              if (n_pw == 0 .or. split_pw_density) then
                if (.not. allocated(rho_blk_partial)) allocate(rho_blk_partial(grid_block_size))
                if (split_pw_density .and. .not. allocated(rho_blk_cross_partial)) allocate(rho_blk_cross_partial(grid_block_size))
                rho_blk_partial(1:npt_blk) = 0.0d0
                if (split_pw_density) rho_blk_cross_partial(1:npt_blk) = 0.0d0
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
                io_s_frag = 1
                io_e_frag = nocc_spin
                do io0 = io_s_frag, io_e_frag, state_block_size
                  nbatch = min(state_block_size, io_e_frag - io0 + 1)
                  call cpu_time(t_psi0)
                  coef_blk_ri(1:nbf, 1:nbatch) = coef_re_full(1:nbf, io0:io0+nbatch-1, ispin)
                  coef_blk_ri(1:nbf, nbatch+1:2*nbatch) = coef_im_full(1:nbf, io0:io0+nbatch-1, ispin)
                  call dgemm('N', 'N', npt_blk, 2*nbatch, nbf, 1.0d0, phi_blk, grid_block_size, &
                             coef_blk_ri, nbf_max, 0.0d0, psi_blk_ri, grid_block_size)
                  psi_blk_re(1:npt_blk, 1:nbatch) = psi_blk_ri(1:npt_blk, 1:nbatch)
                  psi_blk_im(1:npt_blk, 1:nbatch) = psi_blk_ri(1:npt_blk, nbatch+1:2*nbatch)
                  if (split_pw_density) then
                    pw_tmp_z(1:npt_blk, 1:nbatch) = zzero
                    do ipw0 = 1, n_pw, pw_block_size
                      npw_blk = min(pw_block_size, n_pw - ipw0 + 1)
                      coef_pw_blk(1:npw_blk, 1:nbatch) = &
                        dg_frag%coef_pw_full_cache(ipw0:ipw0+npw_blk-1, io0:io0+nbatch-1, ispin)
                      call zgemm('N', 'N', npt_blk, nbatch, npw_blk, zone, phase_cache(1, ipw0), grid_block_size, &
                                 coef_pw_blk, pw_block_size, zone, pw_tmp_z, grid_block_size)
                    end do
                  end if
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
!$omp parallel do private(io, igrid, rho_accum, rho_raw_contrib, occ_factor) schedule(static)
                    do igrid = 1, npt_blk
                      rho_accum = 0.0d0
                      rho_raw_contrib = 0.0d0
!$omp simd reduction(+:rho_accum)
                      do io = 1, nbatch
                        occ_factor = occ_cache(io0 + io - 1)
                        rho_accum = rho_accum + occ_factor * &
                          (psi_blk_re(igrid, io) * psi_blk_re(igrid, io) + psi_blk_im(igrid, io) * psi_blk_im(igrid, io))
                        if (split_pw_density) then
                          rho_raw_contrib = rho_raw_contrib + 2.0d0 * occ_factor * &
                            (psi_blk_re(igrid, io) * real(pw_tmp_z(igrid, io), kind=8) + &
                             psi_blk_im(igrid, io) * aimag(pw_tmp_z(igrid, io)))
                        end if
                      end do
                      rho_blk_partial(igrid) = rho_blk_partial(igrid) + rho_accum
                      if (split_pw_density) rho_blk_cross_partial(igrid) = rho_blk_cross_partial(igrid) + rho_raw_contrib
                    end do
!$omp end parallel do
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
                rho_blk_accum(1:npt_blk) = rho_blk_partial(1:npt_blk)
                if (split_pw_density .and. enable_density_split_charge_probe) then
                  charge_fragment_preowner_local = charge_fragment_preowner_local + sum(rho_blk_accum(1:npt_blk)) * system%hvol
                end if
                if (enable_density_reconstruct_trace .and. block_offset == first_block_offset) then
                  rho_probe_charge = sum(rho_blk_accum(1:npt_blk)) * system%hvol
                  write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4)') "        density frag probe: rank=", dg_frag%id, &
                    " ifrag=", ifrag, " ispin=", ispin, " first_block_charge=", rho_probe_charge
                  flush(6)
                  if (dg_frag%is_frag_root) then
                    block_charge_psi = 0.0d0
                    block_charge_direct = 0.0d0
                    do idx_local = 1, local_grid_count
                      igrid = local_grid_ids(idx_local)
                      block_charge_psi = block_charge_psi + rho_blk_accum(igrid)
                      rho_direct_probe = 0.0d0
                      do iprobe = 1, nbf
                        coeff_re_probe = 0.0d0
                        coeff_im_probe = 0.0d0
                        do io = 1, nbf
                          coeff_re_probe = coeff_re_probe + D_frag_re(iprobe, io, ispin) * phi_blk(igrid, io)
                        end do
                        rho_direct_probe = rho_direct_probe + phi_blk(igrid, iprobe) * coeff_re_probe
                      end do
                      block_charge_direct = block_charge_direct + rho_direct_probe
                    end do
                    block_charge_psi = block_charge_psi * system%hvol
                    block_charge_direct = block_charge_direct * system%hvol
                    write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4,a,1pe12.4)') "        density block-check: rank=", dg_frag%id, &
                      " ifrag=", ifrag, " ispin=", ispin, " psi=", block_charge_psi, " direct=", block_charge_direct
                    flush(6)
                  end if
                end if
              else
                ! n_pw > 0: state-distributed loop, no per-batch bcast
                if (.not. allocated(rho_blk_partial)) allocate(rho_blk_partial(grid_block_size))
                rho_blk_partial(1:npt_blk) = 0.0d0

                do io0 = io_s_frag, io_e_frag, state_block_size
                  nbatch = min(state_block_size, io_e_frag - io0 + 1)

                  call cpu_time(t_psi0)
                  coef_blk_ri(1:nbf, 1:nbatch) = coef_re_full(1:nbf, io0:io0+nbatch-1, ispin)
                  coef_blk_ri(1:nbf, nbatch+1:2*nbatch) = coef_im_full(1:nbf, io0:io0+nbatch-1, ispin)
                  call dgemm('N', 'N', npt_blk, 2*nbatch, nbf, 1.0d0, phi_blk, grid_block_size, &
                             coef_blk_ri, nbf_max, 0.0d0, psi_blk_ri, grid_block_size)
                  psi_blk_re(1:npt_blk, 1:nbatch) = psi_blk_ri(1:npt_blk, 1:nbatch)
                  psi_blk_im(1:npt_blk, 1:nbatch) = psi_blk_ri(1:npt_blk, nbatch+1:2*nbatch)
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
!$omp parallel do private(io, idx_local, igrid, rho_accum, occ_factor) schedule(static)
                  do idx_local = 1, local_grid_count
                    igrid = local_grid_ids(idx_local)
                    rho_accum = 0.0d0
!$omp simd reduction(+:rho_accum)
                    do io = 1, nbatch
                      occ_factor = occ_cache(io0 + io - 1)
                      rho_accum = rho_accum + occ_factor * (psi_blk_re(igrid, io) * psi_blk_re(igrid, io) + &
                                  psi_blk_im(igrid, io) * psi_blk_im(igrid, io))
                    end do
                    rho_blk_partial(igrid) = rho_blk_partial(igrid) + rho_accum
                  end do
!$omp end parallel do
                  call cpu_time(t_rho1)
                  time_project_rho = time_project_rho + (t_rho1 - t_rho0)
                end do  ! io0

                ! AllReduce rho_blk_partial across icomm_frag → rho_blk_accum
                ! Note: comm_summation overwrites rho_blk_accum (does not accumulate into it)
                call cpu_time(t_rho0)
                if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
                  call comm_summation(rho_blk_partial(1:npt_blk), rho_blk_accum(1:npt_blk), npt_blk, dg_frag%icomm_frag)
                else
                  rho_blk_accum(1:npt_blk) = rho_blk_partial(1:npt_blk)
                end if
                call cpu_time(t_rho1)
                time_project_rho_reduce = time_project_rho_reduce + (t_rho1 - t_rho0)
              end if
              ! rho_blk_accum: filled by dgemm-path (n_pw==0) or AllReduce (n_pw>0)
                frag_charge_sumspin = frag_charge_sumspin + sum(rho_blk_accum(1:npt_blk)) * system%hvol
              call cpu_time(t_rho0)
                  call cpu_time(t_post_owner0)
                  do igrid = 1, npt_blk
                    if ((n_pw > 0 .and. .not. split_pw_density) .and. .not. dg_frag%is_frag_root) cycle
                    ixg = ixg_buf(igrid)
                    iyg = iyg_buf(igrid)
                    izg = izg_buf(igrid)
                    if (owner_buf(igrid) /= dg_frag%id) cycle
                    if (ixg < rho_s_x_lo .or. ixg > rho_s_x_hi .or. &
                      iyg < rho_s_y_lo .or. iyg > rho_s_y_hi .or. &
                      izg < rho_s_z_lo .or. izg > rho_s_z_hi) cycle
                    rho_raw_contrib = rho_blk_accum(igrid)
                    if (split_pw_density .and. enable_density_split_charge_probe) charge_fragment_postowner_local = charge_fragment_postowner_local + rho_raw_contrib * system%hvol
                    rho_contrib = rho_raw_contrib
                    rho_mix_accum = 0.0d0
                    if (split_pw_density) then
                      rho_mix_accum = rho_blk_cross_partial(igrid)
                    end if
                    frag_owner_charge_sumspin = frag_owner_charge_sumspin + (rho_contrib + rho_mix_accum) * system%hvol
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
                    rho_bf(bx, by, bz) = rho_bf(bx, by, bz) + rho_contrib + rho_mix_accum
                    rho_s_bf(ixg, iyg, izg, ispin) = rho_s_bf(ixg, iyg, izg, ispin) + rho_contrib + rho_mix_accum
                    if (split_pw_density .and. enable_density_split_charge_probe) then
                      charge_fragment_local = charge_fragment_local + rho_contrib * system%hvol
                      charge_cross_local = charge_cross_local + rho_mix_accum * system%hvol
                    end if
                  end do
                  call cpu_time(t_post_owner1)
                  time_post_owner_accum = time_post_owner_accum + (t_post_owner1 - t_post_owner0)
                  call cpu_time(t_post_remote0)
                  do idx_remote = 1, valid_remote_grid_count
                    if (split_pw_density) cycle
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
                    frag_remote_send_sumspin = frag_remote_send_sumspin + rho_contrib * system%hvol
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
                    rho_send(owner_rank)%f(spin_offset + slot, 1, 1) = rho_send(owner_rank)%f(spin_offset + slot, 1, 1) + rho_contrib
                  end do
                  call cpu_time(t_post_remote1)
                  time_post_remote_pack = time_post_remote_pack + (t_post_remote1 - t_post_remote0)
                call cpu_time(t_rho1)
                time_project_rho = time_project_rho + (t_rho1 - t_rho0)
            end if
          end do
        end do
        if ((n_pw == 0 .or. split_pw_density) .and. enable_density_reconstruct_trace) then
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
              write(*,'(1x,a,i0,a,i0,a,i0,a,3(1pe12.4,1x),a,3(1pe12.4,1x),a,3(1pe12.4,1x))') "        density phi-metric probe: rank=", dg_frag%id, &
                " ifrag=", ifrag, " ispin=", ispin, " hvol_cols=", phi_col_hvol_probe(1), phi_col_hvol_probe(2), phi_col_hvol_probe(3), &
                " full_cols=", phi_col_full_metric_total(1,ispin,i_local), phi_col_full_metric_total(2,ispin,i_local), &
                phi_col_full_metric_total(3,ispin,i_local), " sdiag=", phi_sdiag_probe(1), phi_sdiag_probe(2), phi_sdiag_probe(3)
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
        if (enable_density_reconstruct_trace .and. dg_frag%is_frag_root) then
          write(*,'(1x,a,i0,a,i0,a,1pe12.4)') "        density frag charge sumspin: rank=", dg_frag%id, &
            " ifrag=", ifrag, " charge_sumspin=", frag_charge_sumspin
          write(*,'(1x,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4)') "        density frag coef consistency: rank=", dg_frag%id, &
            " ifrag=", ifrag, " occ_c2_sumspin=", frag_coeff_occ_norm_sumspin, &
            " D_diag_sumspin=", frag_d_diag_sumspin, " TrSD_sumspin=", frag_trace_sumspin
          write(*,'(1x,a,i0,a,i0,a,1pe20.12,a,1pe20.12,a,1pe20.12)') "        density frag comm split: rank=", dg_frag%id, &
            " ifrag=", ifrag, " owner_apply=", frag_owner_charge_sumspin, &
            " remote_send=", frag_remote_send_sumspin, " total=", frag_owner_charge_sumspin + frag_remote_send_sumspin
          flush(6)
        end if
        block_idx_global = block_idx_global + nblocks_ifrag
      end do
    call cpu_time(t_project1)
    time_project = time_project + (t_project1 - t_project0)
    call cpu_time(t_post_split0)
    if (split_pw_density) then
      do ispin = 1, system%nspin
        nocc_spin = dg_frag%nocc_spin(ispin)
        if (nocc_spin <= 0) cycle
        occ_cache(1:nocc_spin) = 1.0d0
        if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
          do io = 1, nocc_spin
            if (io <= size(dg_frag%occ_state, 1) .and. ispin <= size(dg_frag%occ_state, 2)) then
              occ_cache(io) = max(0.0d0, dg_frag%occ_state(io, ispin))
            end if
          end do
        else if (allocated(system%rocc)) then
          do io = 1, nocc_spin
            if (io <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
              occ_cache(io) = max(0.0d0, system%rocc(io, 1, ispin))
            end if
          end do
        end if
        igrid0 = 1
        ngrid = (rho_x_hi - rho_x_lo + 1) * (rho_y_hi - rho_y_lo + 1) * (rho_z_hi - rho_z_lo + 1)
        do while (igrid0 <= ngrid)
          npt_blk = min(grid_block_size, ngrid - igrid0 + 1)
          do igrid = 1, npt_blk
            block_idx_global = igrid0 + igrid - 2
            ixg = rho_x_lo + mod(block_idx_global, rho_x_hi - rho_x_lo + 1)
            rem_xy = block_idx_global / (rho_x_hi - rho_x_lo + 1)
            iyg = rho_y_lo + mod(rem_xy, rho_y_hi - rho_y_lo + 1)
            izg = rho_z_lo + rem_xy / (rho_y_hi - rho_y_lo + 1)
            ix_buf(igrid) = ixg
            iy_buf(igrid) = iyg
            iz_buf(igrid) = izg
          end do
!$omp parallel do private(ipw, theta) schedule(static)
          do igrid = 1, npt_blk
!$omp simd
            do ipw = 1, n_pw
              theta = kpw_hx(ipw) * real(ix_buf(igrid) - 1, 8) + kpw_hy(ipw) * real(iy_buf(igrid) - 1, 8) + &
                      kpw_hz(ipw) * real(iz_buf(igrid) - 1, 8)
              phase_cache(igrid, ipw) = cmplx(cos(theta), sin(theta), kind=8) * inv_sqrt_vol
            end do
          end do
!$omp end parallel do
          rho_blk_accum(1:npt_blk) = 0.0d0
          do io0 = 1, nocc_spin, state_block_size
            nbatch = min(state_block_size, nocc_spin - io0 + 1)
            psi_blk_re(1:npt_blk, 1:nbatch) = 0.0d0
            psi_blk_im(1:npt_blk, 1:nbatch) = 0.0d0
            do ipw0 = 1, n_pw, pw_block_size
              npw_blk = min(pw_block_size, n_pw - ipw0 + 1)
              coef_pw_blk(1:npw_blk, 1:nbatch) = dg_frag%coef_pw_full_cache(ipw0:ipw0+npw_blk-1, io0:io0+nbatch-1, ispin)
              call zgemm('N', 'N', npt_blk, nbatch, npw_blk, zone, phase_cache(1, ipw0), grid_block_size, &
                         coef_pw_blk, pw_block_size, zzero, pw_tmp_z, grid_block_size)
              psi_blk_re(1:npt_blk, 1:nbatch) = psi_blk_re(1:npt_blk, 1:nbatch) + real(pw_tmp_z(1:npt_blk, 1:nbatch), kind=8)
              psi_blk_im(1:npt_blk, 1:nbatch) = psi_blk_im(1:npt_blk, 1:nbatch) + aimag(pw_tmp_z(1:npt_blk, 1:nbatch))
            end do
            do io1 = 1, nbatch
              occ_factor = occ_cache(io0 + io1 - 1)
              if (occ_factor <= 0.0d0) cycle
!$omp parallel do private(rho_accum) schedule(static)
              do igrid = 1, npt_blk
                rho_accum = occ_factor * (psi_blk_re(igrid, io1) * psi_blk_re(igrid, io1) + &
                           psi_blk_im(igrid, io1) * psi_blk_im(igrid, io1))
                rho_blk_accum(igrid) = rho_blk_accum(igrid) + rho_accum
              end do
!$omp end parallel do
            end do
          end do
          do igrid = 1, npt_blk
            ixg = ix_buf(igrid)
            iyg = iy_buf(igrid)
            izg = iz_buf(igrid)
            rho_bf(ixg, iyg, izg) = rho_bf(ixg, iyg, izg) + rho_blk_accum(igrid)
            rho_s_bf(ixg, iyg, izg, ispin) = rho_s_bf(ixg, iyg, izg, ispin) + rho_blk_accum(igrid)
            if (enable_density_split_charge_probe) charge_pw_local = charge_pw_local + rho_blk_accum(igrid) * system%hvol
          end do
          igrid0 = igrid0 + npt_blk
        end do
      end do
    end if
    call cpu_time(t_post_split1)
    time_post_split_pw = time_post_split_pw + (t_post_split1 - t_post_split0)

    call cpu_time(t_post_dealloc0)
    if (allocated(D_frag_re)) deallocate(D_frag_re)
    if (allocated(coef_re_frag)) deallocate(coef_re_frag)
    if (allocated(coef_im_frag)) deallocate(coef_im_frag)
    if (allocated(D_partial_re)) deallocate(D_partial_re)
    if (allocated(coef_re_full)) deallocate(coef_re_full)
    if (allocated(coef_im_full)) deallocate(coef_im_full)
    if (allocated(coef_c_full)) deallocate(coef_c_full)
    if (allocated(coef_c_frag)) deallocate(coef_c_frag)
    if (allocated(coef_full_cache)) deallocate(coef_full_cache)
    if (allocated(coef_pw_view_cache)) deallocate(coef_pw_view_cache)
    if (allocated(overlap_probe_trace)) deallocate(overlap_probe_trace)
    if (allocated(occ_cache)) deallocate(occ_cache)
    if (allocated(occ_sqrt_cache)) deallocate(occ_sqrt_cache)
    if (allocated(rho_blk_partial)) deallocate(rho_blk_partial)
    if (allocated(rho_blk_cross_partial)) deallocate(rho_blk_cross_partial)
    call cpu_time(t_post_dealloc1)
    time_post_dealloc = time_post_dealloc + (t_post_dealloc1 - t_post_dealloc0)
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
    if (split_pw_density .and. enable_density_reconstruct_trace .and. enable_density_split_charge_probe) then
      call comm_summation(charge_fragment_local, charge_fragment_global, dg_frag%icomm)
      call comm_summation(charge_cross_local, charge_cross_global, dg_frag%icomm)
      call comm_summation(charge_pw_local, charge_pw_global, dg_frag%icomm)
      call comm_summation(charge_fragment_preowner_local, charge_fragment_preowner_global, dg_frag%icomm)
      call comm_summation(charge_fragment_postowner_local, charge_fragment_postowner_global, dg_frag%icomm)
      if (dg_frag%id == 0) then
        write(*,'(1x,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') "        density split-charge: frag_preowner=", &
          charge_fragment_preowner_global, " frag_postowner=", charge_fragment_postowner_global, &
          " frag=", charge_fragment_global, " cross=", charge_cross_global, " pw=", charge_pw_global, &
          " total=", charge_fragment_global + charge_cross_global + charge_pw_global
        flush(6)
      end if
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
    call cpu_time(t_setup0)
    call comm_alltoallv(send_flat, send_counts, send_displs, recv_flat, recv_counts, recv_displs, dg_frag%icomm)
    call cpu_time(t_setup1)
    time_comm_exchange = time_comm_exchange + (t_setup1 - t_setup0)
    if (enable_density_reconstruct_trace) then
      write(*,'(1x,a,i0,a)') "        density collective: rank=", dg_frag%id, " stage=after-alltoallv"
      flush(6)
    end if
    call cpu_time(t_setup0)
    do irank = 0, dg_frag%isize - 1
      if (.not. allocated(rho_recv(irank)%f)) cycle
      if (recv_counts(irank + 1) <= 0) then
        deallocate(rho_recv(irank)%f)
        cycle
      end if
      rho_recv(irank)%f(:, 1, 1) = recv_flat(recv_displs(irank + 1)+1:recv_displs(irank + 1)+recv_counts(irank + 1))

      npts = dg_frag%density_recv_map(irank)%npts
      if (npts < 0) npts = -1
      if ((system%nspin + 1) > 0) then
        if (npts < 0 .or. recv_counts(irank + 1) /= (system%nspin + 1) * npts) then
          if (mod(recv_counts(irank + 1), system%nspin + 1) == 0) then
            npts = recv_counts(irank + 1) / (system%nspin + 1)
          end if
        end if
      end if
      if (npts < 0 .or. recv_counts(irank + 1) /= (system%nspin + 1) * npts) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') &
          "[FATAL] density unpack recv size mismatch: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " peer=", irank, &
          " recv_count=", recv_counts(irank + 1), " npts=", npts
        flush(6)
        stop "DG-Fragment RT: density unpack recv size mismatch"
      end if
!$omp parallel do private(slot, ixg, iyg, izg, bx, by, bz, ispin, spin_offset, rho_contrib, rho_raw_contrib) schedule(static)
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
      rho_contrib = rho_raw_contrib
      rho_bf(bx, by, bz) = rho_bf(bx, by, bz) + rho_contrib
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
        rho_contrib = rho_recv(irank)%f(spin_offset + slot, 1, 1)
        rho_s_bf(ixg, iyg, izg, ispin) = rho_s_bf(ixg, iyg, izg, ispin) + rho_contrib
      end do
      end do
!$omp end parallel do
      deallocate(rho_recv(irank)%f)
    end do
    call cpu_time(t_setup1)
    time_comm_unpack = time_comm_unpack + (t_setup1 - t_setup0)
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
    call cpu_time(t_copy0)
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
    ! Rebuild rho from spin channels to keep charge accounting consistent
    ! with the authoritative spin-resolved density field.
    if (system%nspin == 1) then
!$omp parallel do collapse(3) private(ix,iy,iz) schedule(static)
      do iz = rho_z_lo, rho_z_hi
        do iy = rho_y_lo, rho_y_hi
          do ix = rho_x_lo, rho_x_hi
            rho%f(ix, iy, iz) = rho_s(1)%f(ix, iy, iz)
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
              rho_sum_local = rho_sum_local + rho_s(ispin)%f(ix, iy, iz)
            end do
            rho%f(ix, iy, iz) = rho_sum_local
          end do
        end do
      end do
!$omp end parallel do
    end if
    call cpu_time(t_copy1)
    time_copy = time_copy + (t_copy1 - t_copy0)
    if (enable_density_reconstruct_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density trace: stage=after-rho-copy"
      flush(6)
    end if
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
    target_charge = dble(nelec)
    if (allocated(system%rocc)) then
      rocc_total = 0.0d0
      do ispin = 1, min(system%nspin, size(system%rocc, 3))
        do io = 1, size(system%rocc, 1)
          rocc_total = rocc_total + max(0.0d0, system%rocc(io, 1, ispin))
        end do
      end do
      if (rocc_total > 1.0d-14) then
        if (dg_frag%n_frag > 0 .and. rocc_total > 1.5d0 * target_charge) then
          target_charge = rocc_total / dble(dg_frag%n_frag)
        else
          target_charge = rocc_total
        end if
      end if
    end if
    if (enable_density_reconstruct_trace) then
      write(*,'(1x,a,i0,a,i0,a,1pe12.4)') "        density collective: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " stage=after-total-charge-sum total_charge=", total_charge
      write(*,'(1x,a,2(1pe12.4,1x))') "        density normalization target: target_charge nelec=", &
        target_charge, dble(nelec)
      flush(6)
    end if
    dg_frag%elec_num_raw = total_charge
    dg_frag%rho_scale_factor = 1.0d0
    charge_mismatch_abs = abs(total_charge - target_charge)
    charge_mismatch_tol = max(1.0d-8, 1.0d-2 * max(1.0d0, target_charge))
    if (.not. enable_density_renormalize .and. total_charge > 1.0d-14 .and. total_charge == total_charge) then
      if (charge_mismatch_abs > charge_mismatch_tol) then
        if (dg_frag%id == 0) then
          write(*,'(1x,a)') "[FATAL] DG density charge mismatch detected before renormalization."
          write(*,'(1x,a,3(1pe14.6,1x),a,l1)') "        charge raw/target/diff=", total_charge, target_charge, &
            total_charge - target_charge, " allow_renorm=", enable_density_renormalize
          write(*,'(1x,a)') "        set SALMON_DG_ALLOW_DENSITY_RENORMALIZE=1 only for temporary comparison runs"
          flush(6)
        end if
        stop "DG-Fragment RT: density charge mismatch before renormalization"
      end if
    end if
    if (total_charge > 1.0d-14 .and. total_charge == total_charge) then
      scale_rho = target_charge / total_charge
      dg_frag%rho_scale_factor = scale_rho
      if (enable_density_renormalize .and. abs(scale_rho - 1.0d0) > 1.0d-14) then
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

    if (allocated(dg_frag%rho_frag)) dg_frag%rho_frag(:, :, :) = rho%f(:, :, :)
    if (allocated(dg_frag%rho_s_frag)) then
      do ispin = 1, system%nspin
        dg_frag%rho_s_frag(:, :, :, ispin) = rho_s(ispin)%f(:, :, :)
      end do
    end if
    if (allocated(dg_frag%rho_ud_frag)) dg_frag%rho_ud_frag(:, :, :) = (0.0d0, 0.0d0)

    call cpu_time(t_cleanup0)
    deallocate(ix_buf, iy_buf, iz_buf, owner_buf, ixg_buf, iyg_buf, izg_buf)
    deallocate(slot_buf, local_grid_ids, remote_grid_ids, valid_remote_grid_ids)
    deallocate(basis_gid, valid_basis_ids)
    if (allocated(basis_sdiag_probe)) deallocate(basis_sdiag_probe)
    if (allocated(phi_col_metric_total)) deallocate(phi_col_metric_total)
    if (allocated(phi_col_full_metric_total)) deallocate(phi_col_full_metric_total)
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
    call cpu_time(t_cleanup1)
    time_cleanup = time_cleanup + (t_cleanup1 - t_cleanup0)
    call cpu_time(t_total1)
    time_unaccounted = (t_total1 - t_total0) - (time_cache + time_project + time_comm + time_norm + time_copy + time_cleanup)
    time_preproject = max(0.0d0, t_project0 - t_total0)
    time_postproject_precomm = max(0.0d0, t_comm0 - t_project1)
    time_post_precomm_other = max(0.0d0, time_postproject_precomm - time_post_split_pw - time_post_dealloc)
    time_post_precomm_residual = max(0.0d0, time_post_precomm_other - time_post_owner_accum - time_post_remote_pack)
    time_postcomm_prenorm = max(0.0d0, t_norm0 - t_comm1)
    if (enable_density_reconstruct_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
        "        density timing: total=", t_total1 - t_total0, " cache=", time_cache, &
        " project=", time_project, " comm=", time_comm, " norm=", time_norm
      write(*,'(1x,a,3(a,1pe12.4))') "        density timing extra:", &
        " copy=", time_copy, " cleanup=", time_cleanup, " unaccounted=", time_unaccounted
      write(*,'(1x,a,3(a,1pe12.4))') "        density timing phase:", &
        " pre_project=", time_preproject, " post_project_pre_comm=", time_postproject_precomm, &
        " post_comm_pre_norm=", time_postcomm_prenorm
      write(*,'(1x,a,3(a,1pe12.4))') "        density timing phase2:", &
        " split_pw=", time_post_split_pw, " dealloc=", time_post_dealloc, &
        " precomm_other=", time_post_precomm_other
      write(*,'(1x,a,3(a,1pe12.4))') "        density timing phase3:", &
        " owner_accum=", time_post_owner_accum, " remote_pack=", time_post_remote_pack, &
        " residual=", time_post_precomm_residual
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
