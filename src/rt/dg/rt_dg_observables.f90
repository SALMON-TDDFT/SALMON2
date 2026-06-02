  subroutine calculate_observables(dg_frag, system, mg, stencil, ppg, rt, itt, Vh, Vxc, Vpsl, rho)
    use salmon_global, only: theory
    use structures
    use communication, only: comm_summation, comm_get_max, comm_is_root
    use timer, only: timer_begin, timer_end, LOG_CALC_CURRENT
    use rt_dg_fragment_ops, only: apply_momentum_blocks, apply_matrix_blocks_batch, apply_nonlocal_pp_projector_batch, &
                    apply_mixed_hamiltonian, mixed_fp_coupling_active, copy_matrix_blocks_to_complex_dense, &
                    apply_overlap_operator_batch
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
    
    integer :: io, jo, ispin, idir, n, nocc, n_pw, n_tot, max_nocc
    integer :: nstate_use_diag, nocc_diag, nvirt_diag
    integer :: ix, iy, iz
    integer :: ifrag, jfrag, ib, jb, i_idx, j_idx
    integer :: ifrag_obs, local_row_obs, nbasis_obs, global_idx_obs
    integer :: owner_rank_obs, subgroup_root_rank, owner_offset, block_base, block_rem, cutoff
    integer :: iblk, nrow_blk, ncol_blk, n_diag_block_ids, idb
    integer :: env_len, env_status, parse_status
    integer :: current_block_trace_stride, current_block_trace_maxblocks
    integer :: nblk_trace, iblk_trace, ifrag_row_trace, ifrag_col_trace
    integer :: nrow_trace, ncol_trace, ii_trace, jj_trace, ig_i_trace, ig_j_trace
    integer :: idx_trace, idx_base, n_nonzero_blocks
    integer :: mij_audit_stride, mij_audit_topk, mij_audit_max_occ, mij_audit_max_emp
    integer :: mij_audit_dir, nstate_use, nemp, nocc_use, nemp_use
    integer :: occ_start, emp_start, io_occ, je_emp, iocc_global, iemp_global
    integer :: k_top, k_replace, top_spin_count, top_occ_count, top_emp_count
    integer :: iprobe, iblk_probe, ifrag_row_probe, ifrag_col_probe
    integer :: ii_probe, jj_probe, ig_i_probe, ig_j_probe
    integer :: io_probe, je_probe, nrow_probe, ncol_probe, nval_row, nval_col
    integer :: itop, irep, io_state
    logical :: do_interface_check
    logical :: enable_current_block_trace, do_current_block_trace
    logical :: enable_energy_trace, do_energy_trace
    logical :: enable_mij_audit, do_mij_audit, has_esp, enable_mij_block_audit
    real(8), allocatable :: interface_flow(:,:), dndt_frag(:)
    real(8) :: pair_residual, max_pair_residual, charge_balance_residual
    real(8) :: current_tmp, energy_tmp, pw_weight_local, kpw_dir
    real(8) :: jpara_state, jpara_top_abs
    real(8) :: Ac_tot(3), A_squared
    real(8) :: current_local(3), energy_local
    real(8) :: energy_sum_raw
    real(8) :: energy_static_local, energy_kin_local, energy_nl_local, energy_ap_local, energy_a2_local
    real(8) :: energy_static_sum, energy_kin_sum, energy_nl_sum, energy_ap_sum, energy_a2_sum
    real(8) :: energy_kin_diag_local, energy_kin_offdiag_local
    real(8) :: energy_kin_diag_sum, energy_kin_offdiag_sum
    real(8) :: kinetic_diag_abs_local, kinetic_offdiag_abs_local
    real(8) :: kinetic_diag_abs_sum, kinetic_offdiag_abs_sum
    real(8) :: kinetic_apply_diff_local, kinetic_apply_diff_sum
    real(8) :: energy_static_avg, energy_kin_avg, energy_nl_avg, energy_ap_avg, energy_a2_avg
    real(8) :: frag_reduce_factor
    real(8) :: nelec_ref
    real(8) :: ne_density
    real(8) :: occ_factor
    real(8) :: coef_occ_norm_local, coef_occ_norm_global
    real(8) :: rho_ff_local, rho_fp_local, rho_pp_local
    real(8) :: rho_ff_global, rho_fp_global, rho_pp_global
    real(8) :: rho_ff_state, rho_fp_state, rho_pp_state
    real(8) :: csc_occ_identity_norm_local, csc_occ_identity_max_local
    real(8) :: csc_occ_identity_norm_global, csc_occ_identity_max_global
    real(8) :: csc_occ_identity_max_in(1), csc_occ_identity_max_out(1)
    real(8) :: scalar_reduce_in(1), scalar_reduce_out(1)
    real(8) :: occvirt_leakage_local, occvirt_leakage_global
    real(8) :: leak_state_abs2, leak_abs2_max, leak_pair_abs2
    real(8) :: block_tmp, block_norm, block_norm_min, block_norm_max
    real(8) :: mij_audit_ewin, de, abs2_val, abs_val
    real(8) :: mij_sum_abs2, mij_sum_abs2_window, mij_f_proxy
    real(8) :: mij_sum_abs2_total, mij_sum_abs2_window_total, mij_f_proxy_total
    real(8) :: top_abs2_min
    real(8) :: self_abs_sum, iface_abs_sum, blk_abs_sum, iface_frac
    real(8) :: top_blk_abs(3)
    complex(8) :: minus_i
    complex(8) :: blk_m, blk_recon, mij_probe_val
    character(len=64) :: env_value
    complex(8), allocatable :: op_mat(:,:), tmp_mat(:,:), coef_all(:,:), tmp_all(:,:)
    complex(8), allocatable :: coef_occ_diag(:,:), s_coef_occ(:,:), gram_occ(:,:), leak_proj(:,:)
    complex(8), allocatable :: coef_occ_frag(:,:), coef_occ_pw(:,:), s_coef_occ_frag(:,:), s_coef_occ_pw(:,:)
    complex(8), allocatable :: coef_frag_all(:,:), coef_pw_all(:,:)
    complex(8), allocatable :: coef_occ(:,:), coef_emp(:,:), tmp_emp(:,:), mij_mat(:,:)
    real(8), allocatable :: current_block_local(:), current_block_sum(:)
    real(8), allocatable :: top_abs2(:), top_de(:)
    integer, allocatable :: top_spin(:), top_occ(:), top_emp(:)
    logical :: has_nonlocal, use_hmat_complex, use_mixed_current
    logical :: have_occvirt_ref, enable_occvirt_diag
    logical, save :: occvirt_ref_mode_initialized = .false.
    logical, save :: occvirt_ref_legacy_mode = .false.
    logical, save :: observables_env_initialized = .false.
    logical, save :: cfg_enable_occvirt_diag = .false.
    logical, save :: cfg_enable_obs_charge_check = .false.
    integer, save :: cfg_obs_charge_check_stride = 10
    real(8), save :: cfg_obs_charge_check_tol = 1.0d-8
    logical, save :: cfg_enable_current_block_trace = .false.
    integer, save :: cfg_current_block_trace_stride = 20
    integer, save :: cfg_current_block_trace_maxblocks = 0
    logical, save :: cfg_enable_energy_trace = .false.
    integer, save :: cfg_energy_trace_stride = 10
    logical, save :: cfg_enable_mij_audit = .false.
    logical, save :: cfg_enable_mij_block_audit = .false.
    integer, save :: cfg_mij_audit_stride = 50
    integer, save :: cfg_mij_audit_topk = 10
    integer, save :: cfg_mij_audit_max_occ = 0
    integer, save :: cfg_mij_audit_max_emp = 0
    integer, save :: cfg_mij_audit_dir = 3
    real(8), save :: cfg_mij_audit_ewin = 0.0d0
    logical :: use_spatial_A
    logical :: enable_obs_charge_check
    logical :: trace_obs
    integer :: obs_charge_check_stride
    real(8) :: obs_charge_check_tol
    real(8) :: obs_t0, obs_t1
    real(8) :: obs_charge_local, obs_charge_global, obs_charge_diff
    integer :: obs_charge_check_env_len, obs_charge_check_env_status
    real(8), allocatable :: Ap_mat(:,:), A2_mat(:,:)
    integer, allocatable :: diag_block_ids(:)
    logical, allocatable :: fp_row_owned(:), pw_row_owned(:)
    real(8), allocatable :: occ_weight(:)
    complex(8) :: mfp
    real(8), parameter :: unit_dir(3,3) = reshape((/ &
      1.0d0, 0.0d0, 0.0d0, &
      0.0d0, 1.0d0, 0.0d0, &
      0.0d0, 0.0d0, 1.0d0 /), (/3, 3/))
    integer, parameter :: mij_block_probe_count = 3
    integer, parameter :: mij_block_probe_occ(mij_block_probe_count) = (/ 2, 2, 3 /)
    integer, parameter :: mij_block_probe_emp(mij_block_probe_count) = (/ 184, 192, 188 /)
    ! Calculate local observables (only for assigned fragments)
    ! MPI aggregation will sum across all ranks
    current_local = 0.0d0
    coef_occ_norm_local = 0.0d0
    coef_occ_norm_global = 0.0d0
    rho_ff_local = 0.0d0
    rho_fp_local = 0.0d0
    rho_pp_local = 0.0d0
    rho_ff_global = 0.0d0
    rho_fp_global = 0.0d0
    rho_pp_global = 0.0d0
    csc_occ_identity_norm_local = 0.0d0
    csc_occ_identity_max_local = 0.0d0
    csc_occ_identity_norm_global = 0.0d0
    csc_occ_identity_max_global = 0.0d0
    occvirt_leakage_local = 0.0d0
    occvirt_leakage_global = 0.0d0
    have_occvirt_ref = .false.
    jpara_top_abs = 0.0d0
    dg_frag%occvirt_top_occ = 0
    dg_frag%occvirt_top_virt = 0
    dg_frag%occvirt_top_abs2 = 0.0d0
    dg_frag%jpara_top_occ_state = 0
    dg_frag%jpara_top_occ_value = 0.0d0
    energy_local = 0.0d0
    pw_weight_local = 0.0d0
    energy_static_local = 0.0d0
    energy_kin_local = 0.0d0
    energy_nl_local = 0.0d0
    energy_ap_local = 0.0d0
    energy_a2_local = 0.0d0
    energy_kin_diag_local = 0.0d0
    energy_kin_offdiag_local = 0.0d0
    kinetic_diag_abs_local = 0.0d0
    kinetic_offdiag_abs_local = 0.0d0
    kinetic_apply_diff_local = 0.0d0

    n = dg_frag%n_mat_max
    trace_obs = (itt == 1 .and. dg_frag%id == 0)
    if (trace_obs) then
      write(*,'(1x,a,i0,a,i0,a,i0)') '[DG-OBS] enter itt=', itt, &
        ' n_mat=', n, ' nspin=', dg_frag%nspin
      flush(6)
    end if
    use_spatial_A = (trim(theory) == 'single_scale_maxwell_tddft' .and. allocated(system%Ac_micro%v) .and. dg_frag%has_real_space_basis)
    obs_charge_local = 0.0d0
    obs_charge_global = 0.0d0
    obs_charge_diff = 0.0d0

    ! Observables are evaluated every RT step; diagnostic switches are static
    ! process settings, so read them once and reuse the parsed values.
    if (.not. observables_env_initialized) then
      env_value = ''
      call get_environment_variable("SALMON_DG_OBS_CHARGE_CHECK", env_value, length=obs_charge_check_env_len, status=obs_charge_check_env_status)
      if (obs_charge_check_env_status == 0 .and. obs_charge_check_env_len > 0) then
        if (env_value(1:1) == '1' .or. env_value(1:1) == 'y' .or. env_value(1:1) == 'Y' .or. &
            env_value(1:1) == 't' .or. env_value(1:1) == 'T') then
          cfg_enable_obs_charge_check = .true.
        end if
      end if
      env_value = ''
      call get_environment_variable("SALMON_DG_OBS_CHARGE_CHECK_STRIDE", env_value, length=obs_charge_check_env_len, status=obs_charge_check_env_status)
      if (obs_charge_check_env_status == 0 .and. obs_charge_check_env_len > 0) then
        read(env_value(1:obs_charge_check_env_len), *, iostat=parse_status) cfg_obs_charge_check_stride
        if (parse_status /= 0) cfg_obs_charge_check_stride = 10
      end if
      if (cfg_obs_charge_check_stride < 1) cfg_obs_charge_check_stride = 1
      env_value = ''
      call get_environment_variable("SALMON_DG_OBS_CHARGE_CHECK_TOL", env_value, length=obs_charge_check_env_len, status=obs_charge_check_env_status)
      if (obs_charge_check_env_status == 0 .and. obs_charge_check_env_len > 0) then
        read(env_value(1:obs_charge_check_env_len), *, iostat=parse_status) cfg_obs_charge_check_tol
        if (parse_status /= 0) cfg_obs_charge_check_tol = 1.0d-8
      end if

      env_value = ''
      call get_environment_variable("SALMON_DG_CURRENT_BLOCK_TRACE", env_value, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        if (env_value(1:1) == '1' .or. env_value(1:1) == 'y' .or. env_value(1:1) == 'Y' .or. &
            env_value(1:1) == 't' .or. env_value(1:1) == 'T') then
          cfg_enable_current_block_trace = .true.
        end if
      end if
      env_value = ''
      call get_environment_variable("SALMON_DG_CURRENT_BLOCK_TRACE_STRIDE", env_value, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_value(1:env_len), *, iostat=parse_status) cfg_current_block_trace_stride
        if (parse_status /= 0) cfg_current_block_trace_stride = 20
      end if
      if (cfg_current_block_trace_stride < 1) cfg_current_block_trace_stride = 1
      env_value = ''
      call get_environment_variable("SALMON_DG_CURRENT_BLOCK_TRACE_MAXBLOCKS", env_value, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_value(1:env_len), *, iostat=parse_status) cfg_current_block_trace_maxblocks
        if (parse_status /= 0) cfg_current_block_trace_maxblocks = 0
      end if
      if (cfg_current_block_trace_maxblocks < 0) cfg_current_block_trace_maxblocks = 0

      env_value = ''
      call get_environment_variable("SALMON_DG_OBS_ENERGY_TRACE", env_value, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        if (env_value(1:1) == '1' .or. env_value(1:1) == 'y' .or. env_value(1:1) == 'Y' .or. &
            env_value(1:1) == 't' .or. env_value(1:1) == 'T') then
          cfg_enable_energy_trace = .true.
        end if
      end if
      env_value = ''
      call get_environment_variable("SALMON_DG_OBS_ENERGY_TRACE_STRIDE", env_value, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_value(1:env_len), *, iostat=parse_status) cfg_energy_trace_stride
        if (parse_status /= 0) cfg_energy_trace_stride = 10
      end if
      if (cfg_energy_trace_stride < 1) cfg_energy_trace_stride = 1

      env_value = ''
      call get_environment_variable("SALMON_DG_MIJ_AUDIT", env_value, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        if (env_value(1:1) == '1' .or. env_value(1:1) == 'y' .or. env_value(1:1) == 'Y' .or. &
            env_value(1:1) == 't' .or. env_value(1:1) == 'T') then
          cfg_enable_mij_audit = .true.
        end if
      end if
      env_value = ''
      call get_environment_variable("SALMON_DG_MIJ_BLOCK_AUDIT", env_value, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        if (env_value(1:1) == '1' .or. env_value(1:1) == 'y' .or. env_value(1:1) == 'Y' .or. &
            env_value(1:1) == 't' .or. env_value(1:1) == 'T') then
          cfg_enable_mij_block_audit = .true.
        end if
      end if
      env_value = ''
      call get_environment_variable("SALMON_DG_MIJ_AUDIT_STRIDE", env_value, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_value(1:env_len), *, iostat=parse_status) cfg_mij_audit_stride
        if (parse_status /= 0) cfg_mij_audit_stride = 50
      end if
      if (cfg_mij_audit_stride < 1) cfg_mij_audit_stride = 1
      env_value = ''
      call get_environment_variable("SALMON_DG_MIJ_AUDIT_TOPK", env_value, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_value(1:env_len), *, iostat=parse_status) cfg_mij_audit_topk
        if (parse_status /= 0) cfg_mij_audit_topk = 10
      end if
      if (cfg_mij_audit_topk < 1) cfg_mij_audit_topk = 1
      env_value = ''
      call get_environment_variable("SALMON_DG_MIJ_AUDIT_MAX_OCC", env_value, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_value(1:env_len), *, iostat=parse_status) cfg_mij_audit_max_occ
        if (parse_status /= 0) cfg_mij_audit_max_occ = 0
      end if
      if (cfg_mij_audit_max_occ < 0) cfg_mij_audit_max_occ = 0
      env_value = ''
      call get_environment_variable("SALMON_DG_MIJ_AUDIT_MAX_EMP", env_value, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_value(1:env_len), *, iostat=parse_status) cfg_mij_audit_max_emp
        if (parse_status /= 0) cfg_mij_audit_max_emp = 0
      end if
      if (cfg_mij_audit_max_emp < 0) cfg_mij_audit_max_emp = 0
      env_value = ''
      call get_environment_variable("SALMON_DG_MIJ_AUDIT_DIR", env_value, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_value(1:env_len), *, iostat=parse_status) cfg_mij_audit_dir
        if (parse_status /= 0) cfg_mij_audit_dir = 3
      end if
      if (cfg_mij_audit_dir < 1 .or. cfg_mij_audit_dir > 3) cfg_mij_audit_dir = 3
      env_value = ''
      call get_environment_variable("SALMON_DG_MIJ_AUDIT_EWIN", env_value, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_value(1:env_len), *, iostat=parse_status) cfg_mij_audit_ewin
        if (parse_status /= 0) cfg_mij_audit_ewin = 0.0d0
      end if
      if (cfg_mij_audit_ewin < 0.0d0) cfg_mij_audit_ewin = 0.0d0

      observables_env_initialized = .true.
    end if

    enable_obs_charge_check = cfg_enable_obs_charge_check
    obs_charge_check_stride = cfg_obs_charge_check_stride
    obs_charge_check_tol = cfg_obs_charge_check_tol
    enable_current_block_trace = cfg_enable_current_block_trace
    current_block_trace_stride = cfg_current_block_trace_stride
    current_block_trace_maxblocks = cfg_current_block_trace_maxblocks
    enable_energy_trace = cfg_enable_energy_trace
    enable_occvirt_diag = cfg_enable_occvirt_diag
    do_energy_trace = enable_energy_trace .and. (itt == 1 .or. mod(itt, cfg_energy_trace_stride) == 0)
    do_current_block_trace = enable_current_block_trace .and. allocated(dg_frag%momentum_blocks) .and. &
      (itt == 1 .or. mod(itt, current_block_trace_stride) == 0)

    if (.not. occvirt_ref_mode_initialized) then
      env_value = ''
      call get_environment_variable('SALMON_DG_OCCVIRT_DIAG', env_value, length=env_len, status=env_status)
      cfg_enable_occvirt_diag = .false.
      if (env_status == 0 .and. env_len > 0) then
        if (env_value(1:1) == '1' .or. env_value(1:1) == 'y' .or. env_value(1:1) == 'Y' .or. &
            env_value(1:1) == 't' .or. env_value(1:1) == 'T') cfg_enable_occvirt_diag = .true.
      end if

      env_value = ''
      call get_environment_variable('SALMON_DG_OCCVIRT_REF_MODE', env_value, length=env_len, status=env_status)
      occvirt_ref_legacy_mode = .false.
      if (env_status == 0 .and. env_len > 0) then
        cfg_enable_occvirt_diag = .true.
        select case (env_value(1:1))
        case ('l','L','0')
          occvirt_ref_legacy_mode = .true.
        case default
          occvirt_ref_legacy_mode = .false.
        end select
      end if
      occvirt_ref_mode_initialized = .true.
      enable_occvirt_diag = cfg_enable_occvirt_diag
    end if

    enable_mij_audit = cfg_enable_mij_audit
    enable_mij_block_audit = cfg_enable_mij_block_audit
    mij_audit_stride = cfg_mij_audit_stride
    mij_audit_topk = cfg_mij_audit_topk
    mij_audit_max_occ = cfg_mij_audit_max_occ
    mij_audit_max_emp = cfg_mij_audit_max_emp
    mij_audit_dir = cfg_mij_audit_dir
    mij_audit_ewin = cfg_mij_audit_ewin

    do_mij_audit = enable_mij_audit .and. (dg_frag%id == 0) .and. &
      (itt == 1 .or. mod(itt, mij_audit_stride) == 0)

    if (enable_obs_charge_check .and. present(rho)) then
      if (itt == 1 .or. mod(itt, obs_charge_check_stride) == 0) then
        obs_charge_local = 0.0d0
        do iz = mg%is(3), mg%ie(3)
          do iy = mg%is(2), mg%ie(2)
            do ix = mg%is(1), mg%ie(1)
              obs_charge_local = obs_charge_local + rho%f(ix, iy, iz)
            end do
          end do
        end do
        obs_charge_local = obs_charge_local * system%hvol
        scalar_reduce_in(1) = obs_charge_local
        scalar_reduce_out(1) = 0.0d0
        call comm_summation(scalar_reduce_in, scalar_reduce_out, 1, dg_frag%icomm)
        obs_charge_global = scalar_reduce_out(1)
        obs_charge_diff = obs_charge_global - dg_frag%elec_num_raw
        if (dg_frag%id == 0) then
          if (abs(obs_charge_diff) > obs_charge_check_tol) then
            write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
              "[WARN] obs-charge-check: itt=", itt, " rho_integral=", obs_charge_global, &
              " Ne_raw=", dg_frag%elec_num_raw, " diff=", obs_charge_diff
          else
            write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
              "[INFO] obs-charge-check: itt=", itt, " rho_integral=", obs_charge_global, &
              " Ne_raw=", dg_frag%elec_num_raw, " diff=", obs_charge_diff
          end if
          flush(6)
        end if
      end if
    end if
    if (do_current_block_trace) then
      nblk_trace = size(dg_frag%momentum_blocks)
      if (current_block_trace_maxblocks > 0) nblk_trace = min(nblk_trace, current_block_trace_maxblocks)
      if (nblk_trace > 0) then
        allocate(current_block_local(3 * nblk_trace), current_block_sum(3 * nblk_trace))
        current_block_local(:) = 0.0d0
        current_block_sum(:) = 0.0d0
      else
        do_current_block_trace = .false.
      end if
    end if

    do_interface_check = .false.
    if (do_interface_check) then
      allocate(interface_flow(dg_frag%n_frag, dg_frag%n_frag))
      interface_flow = 0.0d0
    end if
    if (n <= 0) then
      current_local = 0.0d0
      energy_local = 0.0d0
      goto 1000
    end if
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = size(dg_frag%coef_pw, 1)
    n_tot = n + n_pw
    if (trace_obs) then
      write(*,'(1x,a,i0,a,i0,a,l1)') '[DG-OBS] basis ready n_pw=', n_pw, &
        ' n_tot=', n_tot, ' mixed_current=', n_pw > 0
      flush(6)
    end if
    if (dg_frag%parallel_mode_orbital .and. n_pw == 0) then
      call calculate_observables_orbital_local(dg_frag, system, mg, ppg, rt, itt)
      if (trace_obs) then
        write(*,'(1x,a)') '[DG-OBS] leave orbital-local'
        flush(6)
      end if
      return
    end if
    nstate_use_diag = dg_frag%nstate_tot
    if (.not. enable_occvirt_diag .and. allocated(dg_frag%coef_ref_all)) then
      deallocate(dg_frag%coef_ref_all)
      dg_frag%coef_ref_ready = .false.
    end if
    if (enable_occvirt_diag .and. allocated(dg_frag%coef_ref_all)) then
      if (size(dg_frag%coef_ref_all, 1) /= n_tot .or. size(dg_frag%coef_ref_all, 2) /= nstate_use_diag .or. &
          size(dg_frag%coef_ref_all, 3) /= dg_frag%nspin) then
        deallocate(dg_frag%coef_ref_all)
        dg_frag%coef_ref_ready = .false.
      end if
    end if
    if (enable_occvirt_diag .and. .not. dg_frag%coef_ref_ready .and. (occvirt_ref_legacy_mode .or. itt == 1)) then
      if (.not. allocated(dg_frag%coef_ref_all)) allocate(dg_frag%coef_ref_all(n_tot, nstate_use_diag, dg_frag%nspin))
      dg_frag%coef_ref_all(:, :, :) = (0.0d0, 0.0d0)
      do ispin = 1, dg_frag%nspin
        dg_frag%coef_ref_all(1:n, 1:nstate_use_diag, ispin) = dg_frag%coef(1:n, 1:nstate_use_diag, ispin)
        if (n_pw > 0) then
          dg_frag%coef_ref_all(n+1:n_tot, 1:nstate_use_diag, ispin) = dg_frag%coef_pw(1:n_pw, 1:nstate_use_diag, ispin)
        end if
      end do
      dg_frag%coef_ref_ready = .true.
      if (.not. occvirt_ref_legacy_mode .and. comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') '[DG-OBS][WARN] t0 reference fallback initialized inside calculate_observables'
        flush(6)
      end if
    end if
    have_occvirt_ref = enable_occvirt_diag .and. dg_frag%coef_ref_ready
    max_nocc = max(1, maxval(dg_frag%nocc_spin(1:dg_frag%nspin)))

    allocate(tmp_mat(n, max_nocc))
    allocate(coef_frag_all(n, max_nocc))
    allocate(occ_weight(max_nocc))
    allocate(fp_row_owned(n))
    allocate(pw_row_owned(max(1, n_pw)))
    if (n_pw > 0) then
      allocate(coef_pw_all(n_pw, max_nocc))
      allocate(coef_all(n_tot, max_nocc), tmp_all(n_tot, max_nocc))
    end if
    minus_i = cmplx(0.0d0, -1.0d0, kind=8)

    ! Current calculation via momentum operator matrix (velocity gauge)
    ! Following conventional RT implementation in density_matrix.f90:
    !   - momentum_mat stores <φ|∇|φ> (gradient operator)
    !   - Current: j = Im[<ψ|∇|ψ>] with factor 2 and normalization by ngrid
    !   - Sign: Testing -2.0 to match conventional RT direction
    if (trace_obs) then
      write(*,'(1x,a)') '[DG-OBS] current start'
      flush(6)
      call cpu_time(obs_t0)
    end if
    call timer_begin(LOG_CALC_CURRENT)
    do ispin = 1, dg_frag%nspin
      nocc = min(dg_frag%nocc_spin(ispin), min(dg_frag%nstate_tot, n))
      if (nocc <= 0) cycle
      use_mixed_current = (n_pw > 0 .and. mixed_fp_coupling_active(dg_frag, ispin))
      ! Coefficients in this routine are fragment-local.  Build an ownership
      ! mask from the rank's assigned fragment and use it only for the final
      ! contraction; the operator application still sees the complete local
      ! fragment vector.
      fp_row_owned(:) = .false.
      do ifrag_obs = dg_frag%ifrag_start, dg_frag%ifrag_end
        nbasis_obs = min(dg_frag%n_basis(ifrag_obs, ispin), size(dg_frag%index_basis, 1))
        do local_row_obs = 1, nbasis_obs
          global_idx_obs = dg_frag%index_basis(local_row_obs, ifrag_obs, ispin)
          if (global_idx_obs < 1 .or. global_idx_obs > n) cycle
          owner_rank_obs = dg_frag%id
          if (dg_frag%parallel_mode_orbital .and. allocated(dg_frag%coef_owner) .and. &
              global_idx_obs <= size(dg_frag%coef_owner, 1) .and. ispin <= size(dg_frag%coef_owner, 2)) then
            ! In orbital mode coefficient rows are distributed by the canonical
            ! coef_owner map.  index_basis may be compacted/remapped, so using
            ! the fragment-local row number here makes observables depend on
            ! the fragment subgroup size.
            owner_rank_obs = dg_frag%coef_owner(global_idx_obs, ispin)
          else if (dg_frag%parallel_mode_orbital .and. dg_frag%isize_frag > 1) then
            subgroup_root_rank = dg_frag%id_array(ifrag_obs) - mod(max(0, dg_frag%id_array(ifrag_obs)), dg_frag%isize_frag)
            block_base = nbasis_obs / dg_frag%isize_frag
            block_rem = mod(nbasis_obs, dg_frag%isize_frag)
            cutoff = (block_base + 1) * block_rem
            if (local_row_obs <= 0) then
              owner_offset = 0
            else if (block_base <= 0) then
              owner_offset = min(local_row_obs - 1, dg_frag%isize_frag - 1)
            else if (local_row_obs <= cutoff) then
              owner_offset = (local_row_obs - 1) / (block_base + 1)
            else
              owner_offset = block_rem + (local_row_obs - cutoff - 1) / block_base
            end if
            owner_rank_obs = subgroup_root_rank + min(owner_offset, dg_frag%isize_frag - 1)
          end if
          if (owner_rank_obs == dg_frag%id) fp_row_owned(global_idx_obs) = .true.
        end do
      end do
      pw_row_owned(:) = .false.
      if (n_pw > 0) then
        do jo = 1, n_pw
          owner_rank_obs = dg_frag%id
          if (allocated(dg_frag%coef_pw_owner)) owner_rank_obs = dg_frag%coef_pw_owner(jo)
          if (owner_rank_obs == dg_frag%id) pw_row_owned(jo) = .true.
        end do
      end if
      coef_frag_all(1:n, 1:nocc) = dg_frag%coef(1:n, 1:nocc, ispin)
      if (n_pw > 0) coef_pw_all(1:n_pw, 1:nocc) = dg_frag%coef_pw(1:n_pw, 1:nocc, ispin)
      if (n_pw > 0) then
        coef_all(1:n_tot, 1:nocc) = (0.0d0, 0.0d0)
        coef_all(1:n, 1:nocc) = coef_frag_all(1:n, 1:nocc)
        coef_all(n+1:n_tot, 1:nocc) = coef_pw_all(1:n_pw, 1:nocc)
      end if
      do i_idx = 1, n
        if (.not. fp_row_owned(i_idx)) cycle
        coef_occ_norm_local = coef_occ_norm_local + sum(abs(coef_frag_all(i_idx, 1:nocc))**2)
      end do
      if (n_pw > 0) then
        do jo = 1, n_pw
          if (.not. pw_row_owned(jo)) cycle
          coef_occ_norm_local = coef_occ_norm_local + sum(abs(coef_pw_all(jo, 1:nocc))**2)
        end do
      end if
      do idir = 1, 3
        ! momentum_mat = <phi|grad|phi>; use the same occupation weights as conventional current.
        if (use_mixed_current) then
          tmp_all(1:n_tot, 1:nocc) = (0.0d0, 0.0d0)
          if (allocated(dg_frag%momentum_blocks)) then
            call apply_momentum_blocks(dg_frag, ispin, unit_dir(:, idir), coef_all(1:n, 1:nocc), tmp_all(1:n, 1:nocc))
          else if (allocated(dg_frag%momentum_mat_c)) then
            if (.not. allocated(op_mat)) allocate(op_mat(n, n))
            op_mat(:, :) = dg_frag%momentum_mat_c(idir, 1:n, 1:n, ispin)
            call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
                       coef_all(1:n, 1:nocc), n, (0.0d0, 0.0d0), tmp_all(1:n, 1:nocc), n)
          else
            if (.not. allocated(op_mat)) allocate(op_mat(n, n))
            op_mat(:, :) = cmplx(dg_frag%momentum_mat(idir, 1:n, 1:n, ispin), 0.0d0, kind=8)
            call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
                       coef_all(1:n, 1:nocc), n, (0.0d0, 0.0d0), tmp_all(1:n, 1:nocc), n)
          end if
          do jo = 1, n_pw
            kpw_dir = dg_frag%k_pw(idir, jo)
            mfp = cmplx(0.0d0, kpw_dir, kind=8)
            tmp_all(n+jo, 1:nocc) = tmp_all(n+jo, 1:nocc) + mfp * coef_all(n+jo, 1:nocc)
            do io = 1, n
              mfp = cmplx(0.0d0, kpw_dir, kind=8) * dg_frag%S_mat_frag_pw(io, jo, ispin)
              tmp_all(io, 1:nocc) = tmp_all(io, 1:nocc) + mfp * coef_all(n+jo, 1:nocc)
              tmp_all(n+jo, 1:nocc) = tmp_all(n+jo, 1:nocc) - conjg(mfp) * coef_all(io, 1:nocc)
            end do
          end do
          current_tmp = 0.0d0
          do io_state = 1, nocc
            do i_idx = 1, n
              if (.not. fp_row_owned(i_idx)) cycle
              current_tmp = current_tmp + occ_weight(io_state) * &
                aimag(conjg(coef_frag_all(i_idx, io_state)) * tmp_all(i_idx, io_state))
            end do
          end do
          if (n_pw > 0) then
            do io_state = 1, nocc
              do jo = 1, n_pw
                if (.not. pw_row_owned(jo)) cycle
                current_tmp = current_tmp + occ_weight(io_state) * &
                  aimag(conjg(coef_pw_all(jo, io_state)) * tmp_all(n+jo, io_state))
              end do
            end do
          end if
        else if (allocated(dg_frag%momentum_blocks)) then
          tmp_mat(:, :) = (0.0d0, 0.0d0)
          call apply_momentum_blocks(dg_frag, ispin, unit_dir(:, idir), coef_frag_all(1:n, 1:nocc), tmp_mat)
        else if (allocated(dg_frag%momentum_mat_c)) then
          if (.not. allocated(op_mat)) allocate(op_mat(n, n))
          op_mat(:, :) = dg_frag%momentum_mat_c(idir, 1:n, 1:n, ispin)
          call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
                     coef_frag_all(1:n, 1:nocc), n, (0.0d0, 0.0d0), tmp_mat, n)
        else
          if (.not. allocated(op_mat)) allocate(op_mat(n, n))
          op_mat(:, :) = cmplx(dg_frag%momentum_mat(idir, 1:n, 1:n, ispin), 0.0d0, kind=8)
          call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
                     coef_frag_all(1:n, 1:nocc), n, (0.0d0, 0.0d0), tmp_mat, n)
        end if
        
        if (.not. use_mixed_current) then
          ! Match conventional calc_current(): paramagnetic current is +Im(psi^* grad psi).
          current_tmp = 0.0d0
          do io_state = 1, nocc
            do i_idx = 1, n
              if (.not. fp_row_owned(i_idx)) cycle
              current_tmp = current_tmp + occ_weight(io_state) * &
                aimag(conjg(coef_frag_all(i_idx, io_state)) * tmp_mat(i_idx, io_state))
            end do
          end do
          if (n_pw > 0) then
            do jo = 1, n_pw
              if (.not. pw_row_owned(jo)) cycle
              kpw_dir = dg_frag%k_pw(idir, jo)
              if (abs(kpw_dir) < 1.0d-15) cycle
              current_tmp = current_tmp + kpw_dir * sum(occ_weight(1:nocc) * abs(coef_pw_all(jo, 1:nocc))**2)
            end do
          end if
        end if

        if (idir == 3) then
          if (use_mixed_current) then
            do io_state = 1, nocc
              jpara_state = 0.0d0
              do i_idx = 1, n
                if (.not. fp_row_owned(i_idx)) cycle
                jpara_state = jpara_state + occ_weight(io_state) * &
                  aimag(conjg(coef_frag_all(i_idx, io_state)) * tmp_all(i_idx, io_state))
              end do
              if (n_pw > 0) then
                do jo = 1, n_pw
                  if (.not. pw_row_owned(jo)) cycle
                  jpara_state = jpara_state + occ_weight(io_state) * &
                    aimag(conjg(coef_pw_all(jo, io_state)) * tmp_all(n+jo, io_state))
                end do
              end if
              if (abs(jpara_state) > jpara_top_abs) then
                jpara_top_abs = abs(jpara_state)
                dg_frag%jpara_top_occ_state = io_state
                dg_frag%jpara_top_occ_value = jpara_state
              end if
            end do
          else
            do io_state = 1, nocc
              jpara_state = 0.0d0
              do i_idx = 1, n
                if (.not. fp_row_owned(i_idx)) cycle
                jpara_state = jpara_state + occ_weight(io_state) * &
                  aimag(conjg(coef_frag_all(i_idx, io_state)) * tmp_mat(i_idx, io_state))
              end do
              if (abs(jpara_state) > jpara_top_abs) then
                jpara_top_abs = abs(jpara_state)
                dg_frag%jpara_top_occ_state = io_state
                dg_frag%jpara_top_occ_value = jpara_state
              end if
            end do
          end if
        end if
        current_local(idir) = current_local(idir) + current_tmp

        if (do_current_block_trace .and. (.not. use_mixed_current) .and. allocated(dg_frag%momentum_blocks)) then
          do iblk_trace = 1, nblk_trace
            ifrag_row_trace = dg_frag%momentum_blocks(iblk_trace)%ifrag_row
            ifrag_col_trace = dg_frag%momentum_blocks(iblk_trace)%ifrag_col
            nrow_trace = dg_frag%n_basis(ifrag_row_trace, ispin)
            ncol_trace = dg_frag%n_basis(ifrag_col_trace, ispin)
            if (nrow_trace <= 0 .or. ncol_trace <= 0) cycle
            block_tmp = 0.0d0
            do io = 1, nocc
              do ii_trace = 1, nrow_trace
                ig_i_trace = dg_frag%index_basis(ii_trace, ifrag_row_trace, ispin)
                if (ig_i_trace < 1 .or. ig_i_trace > n) cycle
                if (.not. fp_row_owned(ig_i_trace)) cycle
                mfp = (0.0d0, 0.0d0)
                do jj_trace = 1, ncol_trace
                  ig_j_trace = dg_frag%index_basis(jj_trace, ifrag_col_trace, ispin)
                  if (ig_j_trace < 1 .or. ig_j_trace > n) cycle
                  mfp = mfp + cmplx(dg_frag%momentum_blocks(iblk_trace)%val(idir, ii_trace, jj_trace, ispin), 0.0d0, kind=8) * &
                    coef_frag_all(ig_j_trace, io)
                end do
                block_tmp = block_tmp + occ_weight(io) * aimag(conjg(coef_frag_all(ig_i_trace, io)) * mfp)
              end do
            end do
            idx_trace = (idir - 1) * nblk_trace + iblk_trace
            current_block_local(idx_trace) = current_block_local(idx_trace) + block_tmp
          end do
        end if
      end do
    end do
    call timer_end(LOG_CALC_CURRENT)
    if (trace_obs) then
      call cpu_time(obs_t1)
      write(*,'(1x,a,1pe12.4)') '[DG-OBS] current done time=', obs_t1 - obs_t0
      flush(6)
      write(*,'(1x,a)') '[DG-OBS] overlap diagnostics start'
      flush(6)
      call cpu_time(obs_t0)
    end if

    do ispin = 1, dg_frag%nspin
      nocc_diag = min(dg_frag%nocc_spin(ispin), min(dg_frag%nstate_tot, n_tot))
      if (nocc_diag <= 0) cycle

      allocate(coef_occ_diag(n_tot, nocc_diag), s_coef_occ(n_tot, nocc_diag), gram_occ(nocc_diag, nocc_diag))
      allocate(coef_occ_frag(n_tot, nocc_diag), coef_occ_pw(n_tot, nocc_diag), &
               s_coef_occ_frag(n_tot, nocc_diag), s_coef_occ_pw(n_tot, nocc_diag))
      coef_occ_diag(:, :) = (0.0d0, 0.0d0)
      coef_occ_diag(1:n, 1:nocc_diag) = dg_frag%coef(1:n, 1:nocc_diag, ispin)
      if (n_pw > 0) coef_occ_diag(n+1:n_tot, 1:nocc_diag) = dg_frag%coef_pw(1:n_pw, 1:nocc_diag, ispin)

      coef_occ_frag(:, :) = (0.0d0, 0.0d0)
      coef_occ_pw(:, :) = (0.0d0, 0.0d0)
      coef_occ_frag(1:n, 1:nocc_diag) = coef_occ_diag(1:n, 1:nocc_diag)
      if (n_pw > 0) coef_occ_pw(n+1:n_tot, 1:nocc_diag) = coef_occ_diag(n+1:n_tot, 1:nocc_diag)

      call apply_overlap_operator_batch(dg_frag, ispin, coef_occ_diag, s_coef_occ, .false.)
      call apply_overlap_operator_batch(dg_frag, ispin, coef_occ_frag, s_coef_occ_frag, .false.)
      call apply_overlap_operator_batch(dg_frag, ispin, coef_occ_pw, s_coef_occ_pw, .false.)
      gram_occ(:, :) = matmul(conjg(transpose(coef_occ_diag(:, :))), s_coef_occ(:, :))
      do io_state = 1, nocc_diag
        gram_occ(io_state, io_state) = gram_occ(io_state, io_state) - (1.0d0, 0.0d0)
      end do
      csc_occ_identity_norm_local = csc_occ_identity_norm_local + sum(abs(gram_occ(:, :))**2)
      csc_occ_identity_max_local = max(csc_occ_identity_max_local, maxval(abs(gram_occ(:, :))))

      do io_state = 1, nocc_diag
        occ_factor = 1.0d0
        if (allocated(system%rocc)) then
          if (io_state <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
            occ_factor = max(0.0d0, system%rocc(io_state, 1, ispin))
          end if
        end if
        if (occ_factor <= 0.0d0) cycle

        rho_ff_state = real(sum(conjg(coef_occ_frag(:, io_state)) * s_coef_occ_frag(:, io_state)), kind=8)
        rho_fp_state = real(sum(conjg(coef_occ_frag(:, io_state)) * s_coef_occ_pw(:, io_state)) + &
                            sum(conjg(coef_occ_pw(:, io_state)) * s_coef_occ_frag(:, io_state)), kind=8)
        rho_pp_state = real(sum(conjg(coef_occ_pw(:, io_state)) * s_coef_occ_pw(:, io_state)), kind=8)

        rho_ff_local = rho_ff_local + occ_factor * rho_ff_state
        rho_fp_local = rho_fp_local + occ_factor * rho_fp_state
        rho_pp_local = rho_pp_local + occ_factor * rho_pp_state
      end do

      if (have_occvirt_ref) then
        nvirt_diag = max(0, nstate_use_diag - nocc_diag)
        if (nvirt_diag > 0) then
          allocate(leak_proj(nstate_use_diag, nocc_diag))
          leak_proj(:, :) = matmul(conjg(transpose(dg_frag%coef_ref_all(1:n_tot, 1:nstate_use_diag, ispin))), s_coef_occ(:, :))
          leak_abs2_max = 0.0d0
          occvirt_leakage_local = occvirt_leakage_local + sum(abs(leak_proj(nocc_diag+1:nstate_use_diag, 1:nocc_diag))**2)
          do io_state = 1, nocc_diag
            leak_state_abs2 = sum(abs(leak_proj(nocc_diag+1:nstate_use_diag, io_state))**2)
            if (leak_state_abs2 > dg_frag%occvirt_top_abs2) then
              dg_frag%occvirt_top_abs2 = leak_state_abs2
              dg_frag%occvirt_top_occ = io_state
            end if
            do jo = nocc_diag + 1, nstate_use_diag
              leak_pair_abs2 = abs(leak_proj(jo, io_state))**2
              if (leak_pair_abs2 > leak_abs2_max) then
                leak_abs2_max = leak_pair_abs2
                dg_frag%occvirt_top_occ = io_state
                dg_frag%occvirt_top_virt = jo
                dg_frag%occvirt_top_abs2 = leak_pair_abs2
              end if
            end do
          end do
          deallocate(leak_proj)
        end if
      end if

      deallocate(coef_occ_diag, s_coef_occ, gram_occ)
      deallocate(coef_occ_frag, coef_occ_pw, s_coef_occ_frag, s_coef_occ_pw)
    end do

    scalar_reduce_in(1) = csc_occ_identity_norm_local
    scalar_reduce_out(1) = 0.0d0
    call comm_summation(scalar_reduce_in, scalar_reduce_out, 1, dg_frag%icomm)
    csc_occ_identity_norm_global = scalar_reduce_out(1)
    csc_occ_identity_max_in(1) = csc_occ_identity_max_local
    call comm_get_max(csc_occ_identity_max_in, csc_occ_identity_max_out, 1, dg_frag%icomm)
    csc_occ_identity_max_global = csc_occ_identity_max_out(1)
    scalar_reduce_in(1) = occvirt_leakage_local
    scalar_reduce_out(1) = 0.0d0
    call comm_summation(scalar_reduce_in, scalar_reduce_out, 1, dg_frag%icomm)
    occvirt_leakage_global = scalar_reduce_out(1)
    scalar_reduce_in(1) = rho_ff_local
    scalar_reduce_out(1) = 0.0d0
    call comm_summation(scalar_reduce_in, scalar_reduce_out, 1, dg_frag%icomm)
    rho_ff_global = scalar_reduce_out(1)
    scalar_reduce_in(1) = rho_fp_local
    scalar_reduce_out(1) = 0.0d0
    call comm_summation(scalar_reduce_in, scalar_reduce_out, 1, dg_frag%icomm)
    rho_fp_global = scalar_reduce_out(1)
    scalar_reduce_in(1) = rho_pp_local
    scalar_reduce_out(1) = 0.0d0
    call comm_summation(scalar_reduce_in, scalar_reduce_out, 1, dg_frag%icomm)
    rho_pp_global = scalar_reduce_out(1)
    if (trace_obs) then
      call cpu_time(obs_t1)
      write(*,'(1x,a,1pe12.4)') '[DG-OBS] overlap diagnostics done time=', obs_t1 - obs_t0
      flush(6)
    end if

    if (do_mij_audit) then
      has_esp = allocated(dg_frag%esp) .and. size(dg_frag%esp, 2) >= dg_frag%nspin
      mij_sum_abs2_total = 0.0d0
      mij_sum_abs2_window_total = 0.0d0
      mij_f_proxy_total = 0.0d0
      top_spin_count = 0
      top_occ_count = 0
      top_emp_count = 0
      allocate(top_abs2(mij_audit_topk), top_de(mij_audit_topk), top_spin(mij_audit_topk), top_occ(mij_audit_topk), top_emp(mij_audit_topk))
      top_abs2(:) = -1.0d0
      top_de(:) = -1.0d0
      top_spin(:) = 0
      top_occ(:) = 0
      top_emp(:) = 0

      do ispin = 1, dg_frag%nspin
        nstate_use = min(dg_frag%nstate_tot, n)
        nocc = min(dg_frag%nocc_spin(ispin), nstate_use)
        nemp = max(0, nstate_use - nocc)
        if (nocc <= 0 .or. nemp <= 0) cycle

        nocc_use = nocc
        if (mij_audit_max_occ > 0) nocc_use = min(nocc_use, mij_audit_max_occ)
        nemp_use = nemp
        if (mij_audit_max_emp > 0) nemp_use = min(nemp_use, mij_audit_max_emp)
        if (nocc_use <= 0 .or. nemp_use <= 0) cycle

        ! Prefer frontier states: highest occupied and lowest empty states.
        occ_start = max(1, nocc - nocc_use + 1)
        emp_start = nocc + 1

        allocate(coef_occ(n, nocc_use), coef_emp(n, nemp_use), tmp_emp(n, nemp_use), mij_mat(nocc_use, nemp_use))
        coef_occ(:, :) = dg_frag%coef(1:n, occ_start:occ_start+nocc_use-1, ispin)
        coef_emp(:, :) = dg_frag%coef(1:n, emp_start:emp_start+nemp_use-1, ispin)
        tmp_emp(:, :) = (0.0d0, 0.0d0)

        if (allocated(dg_frag%momentum_blocks)) then
          call apply_momentum_blocks(dg_frag, ispin, unit_dir(:, mij_audit_dir), coef_emp, tmp_emp)
        else if (allocated(dg_frag%momentum_mat_c)) then
          if (.not. allocated(op_mat)) allocate(op_mat(n, n))
          op_mat(:, :) = dg_frag%momentum_mat_c(mij_audit_dir, 1:n, 1:n, ispin)
          call zgemm('N', 'N', n, nemp_use, n, (1.0d0, 0.0d0), op_mat, n, coef_emp, n, (0.0d0, 0.0d0), tmp_emp, n)
        else
          if (.not. allocated(op_mat)) allocate(op_mat(n, n))
          op_mat(:, :) = cmplx(dg_frag%momentum_mat(mij_audit_dir, 1:n, 1:n, ispin), 0.0d0, kind=8)
          call zgemm('N', 'N', n, nemp_use, n, (1.0d0, 0.0d0), op_mat, n, coef_emp, n, (0.0d0, 0.0d0), tmp_emp, n)
        end if

        mij_mat(:, :) = matmul(conjg(transpose(coef_occ(:, :))), tmp_emp(:, :))

        if (enable_mij_block_audit) then
          do iprobe = 1, mij_block_probe_count
            iocc_global = mij_block_probe_occ(iprobe)
            iemp_global = mij_block_probe_emp(iprobe)
            if (iocc_global < occ_start .or. iocc_global > occ_start + nocc_use - 1) cycle
            if (iemp_global < emp_start .or. iemp_global > emp_start + nemp_use - 1) cycle

            io_probe = iocc_global - occ_start + 1
            je_probe = iemp_global - emp_start + 1
            mij_probe_val = mij_mat(io_probe, je_probe)
            blk_recon = (0.0d0, 0.0d0)
            self_abs_sum = 0.0d0
            iface_abs_sum = 0.0d0
            top_blk_abs(:) = 0.0d0

            do itop = 1, 3
              top_spin(itop) = 0
              top_occ(itop) = 0
              top_emp(itop) = 0
            end do

            if (allocated(dg_frag%momentum_blocks)) then
              do iblk_probe = 1, size(dg_frag%momentum_blocks)
                ifrag_row_probe = dg_frag%momentum_blocks(iblk_probe)%ifrag_row
                ifrag_col_probe = dg_frag%momentum_blocks(iblk_probe)%ifrag_col
                nrow_probe = dg_frag%n_basis(ifrag_row_probe, ispin)
                ncol_probe = dg_frag%n_basis(ifrag_col_probe, ispin)
                if (nrow_probe <= 0 .or. ncol_probe <= 0) cycle

                nval_row = size(dg_frag%momentum_blocks(iblk_probe)%val, 2)
                nval_col = size(dg_frag%momentum_blocks(iblk_probe)%val, 3)
                blk_m = (0.0d0, 0.0d0)
                do ii_probe = 1, min(nrow_probe, nval_row)
                  ig_i_probe = dg_frag%index_basis(ii_probe, ifrag_row_probe, ispin)
                  if (ig_i_probe < 1 .or. ig_i_probe > n) cycle
                  do jj_probe = 1, min(ncol_probe, nval_col)
                    ig_j_probe = dg_frag%index_basis(jj_probe, ifrag_col_probe, ispin)
                    if (ig_j_probe < 1 .or. ig_j_probe > n) cycle
                    blk_m = blk_m + conjg(coef_occ(ig_i_probe, io_probe)) * &
                      cmplx(dg_frag%momentum_blocks(iblk_probe)%val(mij_audit_dir, ii_probe, jj_probe, ispin), 0.0d0, kind=8) * &
                      coef_emp(ig_j_probe, je_probe)
                  end do
                end do

                blk_abs_sum = abs(blk_m)
                blk_recon = blk_recon + blk_m
                if (ifrag_row_probe == ifrag_col_probe) then
                  self_abs_sum = self_abs_sum + blk_abs_sum
                else
                  iface_abs_sum = iface_abs_sum + blk_abs_sum
                end if

                irep = 1
                do itop = 2, 3
                  if (top_blk_abs(itop) < top_blk_abs(irep)) irep = itop
                end do
                if (blk_abs_sum > top_blk_abs(irep)) then
                  top_blk_abs(irep) = blk_abs_sum
                  top_spin(irep) = iblk_probe
                  top_occ(irep) = ifrag_row_probe
                  top_emp(irep) = ifrag_col_probe
                end if
              end do
            else
              do ifrag_row_probe = 1, dg_frag%n_frag
                nrow_probe = dg_frag%n_basis(ifrag_row_probe, ispin)
                if (nrow_probe <= 0) cycle
                do ifrag_col_probe = 1, dg_frag%n_frag
                  ncol_probe = dg_frag%n_basis(ifrag_col_probe, ispin)
                  if (ncol_probe <= 0) cycle
                  blk_m = (0.0d0, 0.0d0)
                  do ii_probe = 1, nrow_probe
                    ig_i_probe = dg_frag%index_basis(ii_probe, ifrag_row_probe, ispin)
                    if (ig_i_probe < 1 .or. ig_i_probe > n) cycle
                    do jj_probe = 1, ncol_probe
                      ig_j_probe = dg_frag%index_basis(jj_probe, ifrag_col_probe, ispin)
                      if (ig_j_probe < 1 .or. ig_j_probe > n) cycle
                      if (allocated(dg_frag%momentum_mat_c)) then
                        blk_m = blk_m + conjg(coef_occ(ig_i_probe, io_probe)) * &
                          dg_frag%momentum_mat_c(mij_audit_dir, ig_i_probe, ig_j_probe, ispin) * coef_emp(ig_j_probe, je_probe)
                      else
                        blk_m = blk_m + conjg(coef_occ(ig_i_probe, io_probe)) * &
                          cmplx(dg_frag%momentum_mat(mij_audit_dir, ig_i_probe, ig_j_probe, ispin), 0.0d0, kind=8) * &
                          coef_emp(ig_j_probe, je_probe)
                      end if
                    end do
                  end do

                  blk_abs_sum = abs(blk_m)
                  blk_recon = blk_recon + blk_m
                  if (ifrag_row_probe == ifrag_col_probe) then
                    self_abs_sum = self_abs_sum + blk_abs_sum
                  else
                    iface_abs_sum = iface_abs_sum + blk_abs_sum
                  end if

                  irep = 1
                  do itop = 2, 3
                    if (top_blk_abs(itop) < top_blk_abs(irep)) irep = itop
                  end do
                  if (blk_abs_sum > top_blk_abs(irep)) then
                    top_blk_abs(irep) = blk_abs_sum
                    top_spin(irep) = (ifrag_row_probe - 1) * dg_frag%n_frag + ifrag_col_probe
                    top_occ(irep) = ifrag_row_probe
                    top_emp(irep) = ifrag_col_probe
                  end if
                end do
              end do
            end if

            if (self_abs_sum + iface_abs_sum > 0.0d0) then
              iface_frac = iface_abs_sum / (self_abs_sum + iface_abs_sum)
            else
              iface_frac = 0.0d0
            end if

            de = -1.0d0
            if (has_esp .and. size(dg_frag%esp, 1) >= iemp_global) then
              de = dg_frag%esp(iemp_global, ispin) - dg_frag%esp(iocc_global, ispin)
            end if
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
              "[DG-MIJ-BLOCK-AUDIT] itt=", itt, " spin=", ispin, " iocc=", iocc_global, " iemp=", iemp_global, &
              " absM=", abs(mij_probe_val), " absRecon=", abs(blk_recon), " absDiff=", abs(mij_probe_val - blk_recon), &
              " selfAbs=", self_abs_sum, " ifaceAbs=", iface_abs_sum, " ifaceFrac=", iface_frac
            write(*,'(1x,a,1pe14.6)') "[DG-MIJ-BLOCK-AUDIT] de=", de

            do itop = 1, 3
              if (top_blk_abs(itop) <= 0.0d0) cycle
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,1pe14.6)') &
                "[DG-MIJ-BLOCK-AUDIT-TOPBLK] rank=", itop, " block=", top_spin(itop), " ifrag_row=", top_occ(itop), &
                " ifrag_col=", top_emp(itop), " absBlk=", top_blk_abs(itop)
            end do
          end do
        end if

        mij_sum_abs2 = sum(abs(mij_mat(:, :))**2)
        mij_sum_abs2_total = mij_sum_abs2_total + mij_sum_abs2

        mij_sum_abs2_window = 0.0d0
        mij_f_proxy = 0.0d0
        if (has_esp .and. size(dg_frag%esp, 1) >= nstate_use) then
          do io_occ = 1, nocc_use
            iocc_global = occ_start + io_occ - 1
            do je_emp = 1, nemp_use
              iemp_global = emp_start + je_emp - 1
              abs2_val = abs(mij_mat(io_occ, je_emp))**2
              de = dg_frag%esp(iemp_global, ispin) - dg_frag%esp(iocc_global, ispin)
              if (mij_audit_ewin > 0.0d0) then
                if (de > 0.0d0 .and. de <= mij_audit_ewin) mij_sum_abs2_window = mij_sum_abs2_window + abs2_val
              end if
              if (de > 1.0d-10) mij_f_proxy = mij_f_proxy + 2.0d0 * abs2_val / de
            end do
          end do
        end if
        mij_sum_abs2_window_total = mij_sum_abs2_window_total + mij_sum_abs2_window
        mij_f_proxy_total = mij_f_proxy_total + mij_f_proxy
        top_spin_count = top_spin_count + nocc_use
        top_occ_count = top_occ_count + nocc_use
        top_emp_count = top_emp_count + nemp_use

        do io_occ = 1, nocc_use
          iocc_global = occ_start + io_occ - 1
          do je_emp = 1, nemp_use
            iemp_global = emp_start + je_emp - 1
            abs2_val = abs(mij_mat(io_occ, je_emp))**2
            k_replace = 0
            top_abs2_min = huge(1.0d0)
            do k_top = 1, mij_audit_topk
              if (top_abs2(k_top) < 0.0d0) then
                k_replace = k_top
                exit
              end if
              if (top_abs2(k_top) < top_abs2_min) then
                top_abs2_min = top_abs2(k_top)
                k_replace = k_top
              end if
            end do
            if (k_replace > 0) then
              if (top_abs2(k_replace) < 0.0d0 .or. abs2_val > top_abs2(k_replace)) then
                top_abs2(k_replace) = abs2_val
                if (has_esp .and. size(dg_frag%esp, 1) >= iemp_global) then
                  top_de(k_replace) = dg_frag%esp(iemp_global, ispin) - dg_frag%esp(iocc_global, ispin)
                else
                  top_de(k_replace) = -1.0d0
                end if
                top_spin(k_replace) = ispin
                top_occ(k_replace) = iocc_global
                top_emp(k_replace) = iemp_global
              end if
            end if
          end do
        end do

        deallocate(coef_occ, coef_emp, tmp_emp, mij_mat)
      end do

      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
        "[DG-MIJ-AUDIT] itt=", itt, " dir=", mij_audit_dir, " nspin=", dg_frag%nspin, &
        " occ_used_total=", top_occ_count, " emp_used_total=", top_emp_count, &
        " sum_abs2=", mij_sum_abs2_total, " sum_abs2_ewin=", mij_sum_abs2_window_total, " f_proxy=", mij_f_proxy_total
      if (mij_audit_ewin > 0.0d0) then
        write(*,'(1x,a,1pe14.6)') "[DG-MIJ-AUDIT] ewin=", mij_audit_ewin
      end if
      do k_top = 1, mij_audit_topk
        if (top_abs2(k_top) < 0.0d0) cycle
        abs_val = sqrt(max(0.0d0, top_abs2(k_top)))
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
          "[DG-MIJ-AUDIT-TOPK] rank=", k_top, " spin=", top_spin(k_top), " iocc=", top_occ(k_top), &
          " iemp=", top_emp(k_top), " absM=", abs_val, " absM2=", top_abs2(k_top), " de=", top_de(k_top)
      end do
      flush(6)
      deallocate(top_abs2, top_de, top_spin, top_occ, top_emp)
    end if
    
      ! Get vector potential at current time for energy calculation
      if (trace_obs) then
        write(*,'(1x,a)') '[DG-OBS] energy start'
        flush(6)
        call cpu_time(obs_t0)
      end if
      Ac_tot = rt%Ac_tot(:, itt)
      A_squared = Ac_tot(1)**2 + Ac_tot(2)**2 + Ac_tot(3)**2
      
      call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_tot)
      has_nonlocal = dg_frag%has_nl_cache

      ! Calculate total energy: E = <ψ|H(t)|ψ>
      ! H(t) = H_0 - i*A(t)·∇ + A²(t)/2 + V_NL(A)
      do ispin = 1, dg_frag%nspin
      nocc = min(dg_frag%nocc_spin(ispin), min(dg_frag%nstate_tot, n))
      if (nocc <= 0) cycle
      occ_weight(:) = 0.0d0
      call get_dg_spin_occ_info(dg_frag, system, ispin, occ_weight, nocc)
      fp_row_owned(:) = .false.
      do ifrag_obs = dg_frag%ifrag_start, dg_frag%ifrag_end
        nbasis_obs = min(dg_frag%n_basis(ifrag_obs, ispin), size(dg_frag%index_basis, 1))
        do local_row_obs = 1, nbasis_obs
          global_idx_obs = dg_frag%index_basis(local_row_obs, ifrag_obs, ispin)
          if (global_idx_obs < 1 .or. global_idx_obs > n) cycle
          owner_rank_obs = dg_frag%id
          if (dg_frag%parallel_mode_orbital .and. allocated(dg_frag%coef_owner) .and. &
              global_idx_obs <= size(dg_frag%coef_owner, 1) .and. ispin <= size(dg_frag%coef_owner, 2)) then
            ! Use the same row ownership as propagation and block application.
            ! The compact global row is index_basis(local_row, fragment, spin),
            ! not the fragment-local row number.
            owner_rank_obs = dg_frag%coef_owner(global_idx_obs, ispin)
          else if (dg_frag%parallel_mode_orbital .and. dg_frag%isize_frag > 1) then
            subgroup_root_rank = dg_frag%id_array(ifrag_obs) - mod(max(0, dg_frag%id_array(ifrag_obs)), dg_frag%isize_frag)
            block_base = nbasis_obs / dg_frag%isize_frag
            block_rem = mod(nbasis_obs, dg_frag%isize_frag)
            cutoff = (block_base + 1) * block_rem
            if (local_row_obs <= 0) then
              owner_offset = 0
            else if (block_base <= 0) then
              owner_offset = min(local_row_obs - 1, dg_frag%isize_frag - 1)
            else if (local_row_obs <= cutoff) then
              owner_offset = (local_row_obs - 1) / (block_base + 1)
            else
              owner_offset = block_rem + (local_row_obs - cutoff - 1) / block_base
            end if
            owner_rank_obs = subgroup_root_rank + min(owner_offset, dg_frag%isize_frag - 1)
          end if
          if (owner_rank_obs == dg_frag%id) fp_row_owned(global_idx_obs) = .true.
        end do
      end do
      pw_row_owned(:) = .false.
      if (n_pw > 0) then
        do jo = 1, n_pw
          owner_rank_obs = dg_frag%id
          if (allocated(dg_frag%coef_pw_owner)) owner_rank_obs = dg_frag%coef_pw_owner(jo)
          if (owner_rank_obs == dg_frag%id) pw_row_owned(jo) = .true.
        end do
      end if
      coef_frag_all(1:n, 1:nocc) = dg_frag%coef(1:n, 1:nocc, ispin)
      if (n_pw > 0) coef_pw_all(1:n_pw, 1:nocc) = dg_frag%coef_pw(1:n_pw, 1:nocc, ispin)
      if (n_pw > 0) then
        coef_all(1:n_tot, 1:nocc) = (0.0d0, 0.0d0)
        tmp_all(1:n_tot, 1:nocc) = (0.0d0, 0.0d0)
        coef_all(1:n, 1:nocc) = coef_frag_all(1:n, 1:nocc)
        coef_all(n+1:n_tot, 1:nocc) = coef_pw_all(1:n_pw, 1:nocc)
      end if
      use_hmat_complex = allocated(dg_frag%H_mat_c) .and. allocated(dg_frag%phi_frag_c)
      if (allocated(op_mat)) op_mat(:, :) = (0.0d0, 0.0d0)
      if (use_hmat_complex .or. (.not. allocated(dg_frag%H_mat_blocks)) .or. use_spatial_A .or. do_interface_check) then
        if (.not. allocated(op_mat)) allocate(op_mat(n, n))
        op_mat(:, :) = (0.0d0, 0.0d0)
        if (use_hmat_complex) then
          op_mat(:, :) = dg_frag%H_mat_c(1:n, 1:n, ispin)
        else if (.not. allocated(dg_frag%H_mat_blocks)) then
          op_mat(:, :) = cmplx(dg_frag%H_mat(1:n, 1:n, ispin), 0.0d0, kind=8)
        end if
      end if
      if (use_spatial_A) then
        if (.not. allocated(Ap_mat)) allocate(Ap_mat(n, n), A2_mat(n, n))
        call build_spatial_A_coupling_matrices(dg_frag, system, mg, stencil, ispin, Ap_mat, A2_mat)
        op_mat(:, :) = op_mat(:, :) + cmplx(A2_mat(:, :), 0.0d0, kind=8)
        op_mat(:, :) = op_mat(:, :) + minus_i * cmplx(Ap_mat(:, :), 0.0d0, kind=8)
      else
        if (n_pw > 0 .and. allocated(dg_frag%H_mat_frag_pw) .and. mixed_fp_coupling_active(dg_frag, ispin)) then
          call apply_mixed_hamiltonian(dg_frag, ispin, coef_all(1:n_tot, 1:nocc), tmp_all(1:n_tot, 1:nocc))
          if (has_nonlocal) then
            call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all(1:n, 1:nocc), &
              tmp_all(1:n, 1:nocc))
          end if
          tmp_all(1:n_tot, 1:nocc) = tmp_all(1:n_tot, 1:nocc) + 0.5d0 * A_squared * coef_all(1:n_tot, 1:nocc)
          if (allocated(dg_frag%momentum_blocks)) then
            call apply_momentum_blocks(dg_frag, ispin, Ac_tot, coef_all(1:n, 1:nocc), tmp_all(1:n, 1:nocc))
          else
            do idir = 1, 3
              if (allocated(dg_frag%momentum_mat_c)) then
                if (.not. allocated(op_mat)) allocate(op_mat(n, n))
                op_mat(:, :) = dg_frag%momentum_mat_c(idir, 1:n, 1:n, ispin)
              else
                if (.not. allocated(op_mat)) allocate(op_mat(n, n))
                op_mat(:, :) = cmplx(dg_frag%momentum_mat(idir, 1:n, 1:n, ispin), 0.0d0, kind=8)
              end if
              call zgemm('N', 'N', n, nocc, n, minus_i * Ac_tot(idir), op_mat, n, &
                         coef_all(1:n, 1:nocc), n, (1.0d0, 0.0d0), tmp_all(1:n, 1:nocc), n)
            end do
          end if
          do jo = 1, n_pw
            kpw_dir = dot_product(Ac_tot(1:3), dg_frag%k_pw(1:3, jo))
            mfp = cmplx(0.0d0, kpw_dir, kind=8)
            tmp_all(n+jo, 1:nocc) = tmp_all(n+jo, 1:nocc) + minus_i * mfp * coef_all(n+jo, 1:nocc)
            do io = 1, n
              mfp = cmplx(0.0d0, kpw_dir, kind=8) * dg_frag%S_mat_frag_pw(io, jo, ispin)
              tmp_all(io, 1:nocc) = tmp_all(io, 1:nocc) + minus_i * mfp * coef_all(n+jo, 1:nocc)
              tmp_all(n+jo, 1:nocc) = tmp_all(n+jo, 1:nocc) + minus_i * (-conjg(mfp)) * coef_all(io, 1:nocc)
            end do
          end do
        else if (allocated(dg_frag%momentum_blocks)) then
          tmp_mat(:, :) = (0.0d0, 0.0d0)
          if (.not. use_hmat_complex .and. allocated(dg_frag%H_mat_blocks)) then
            call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, coef_frag_all(1:n, 1:nocc), tmp_mat)
            if (has_nonlocal) then
              call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_frag_all(1:n, 1:nocc), &
                tmp_mat)
            end if
            tmp_mat(1:n, 1:nocc) = tmp_mat(1:n, 1:nocc) + 0.5d0 * A_squared * coef_frag_all(1:n, 1:nocc)
          else
            if (.not. allocated(op_mat)) allocate(op_mat(n, n))
            do io = 1, n
              op_mat(io, io) = op_mat(io, io) + 0.5d0 * A_squared
            end do
            call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
                       coef_frag_all(1:n, 1:nocc), n, (0.0d0, 0.0d0), tmp_mat, n)
            if (has_nonlocal) then
              call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_frag_all(1:n, 1:nocc), &
                tmp_mat)
            end if
          end if
          if (.not. allocated(op_mat)) allocate(op_mat(n, n))
          op_mat(:, 1:nocc) = (0.0d0, 0.0d0)
          call apply_momentum_blocks(dg_frag, ispin, Ac_tot, coef_frag_all(1:n, 1:nocc), op_mat(:, 1:nocc))
          tmp_mat(:, :) = tmp_mat(:, :) + minus_i * op_mat(:, 1:nocc)
        else
          if (.not. allocated(op_mat)) allocate(op_mat(n, n))
          do io = 1, n
            op_mat(io, io) = op_mat(io, io) + 0.5d0 * A_squared
          end do
          do idir = 1, 3
            if (allocated(dg_frag%momentum_mat_c)) then
              op_mat(:, :) = op_mat(:, :) + minus_i * Ac_tot(idir) * dg_frag%momentum_mat_c(idir, 1:n, 1:n, ispin)
            else
              op_mat(:, :) = op_mat(:, :) + minus_i * Ac_tot(idir) * dg_frag%momentum_mat(idir, 1:n, 1:n, ispin)
            end if
          end do
        end if
      end if

      if (do_interface_check) then
        do ifrag = 1, dg_frag%n_frag
          do jfrag = 1, dg_frag%n_frag
            if (jfrag == ifrag) cycle
            do io = 1, nocc
              do ib = 1, dg_frag%n_basis(ifrag, ispin)
                i_idx = dg_frag%index_basis(ib, ifrag, ispin)
                if (i_idx < 1 .or. i_idx > n) cycle
                do jb = 1, dg_frag%n_basis(jfrag, ispin)
                  j_idx = dg_frag%index_basis(jb, jfrag, ispin)
                  if (j_idx < 1 .or. j_idx > n) cycle
                  interface_flow(ifrag, jfrag) = interface_flow(ifrag, jfrag) + &
                    2.0d0 * aimag(conjg(coef_frag_all(i_idx, io)) * op_mat(i_idx, j_idx) * &
                                  coef_frag_all(j_idx, io))
                end do
              end do
            end do
          end do
        end do
      end if

      if (n_pw > 0 .and. allocated(dg_frag%H_mat_frag_pw) .and. mixed_fp_coupling_active(dg_frag, ispin) .and. .not. use_spatial_A) then
        energy_tmp = 0.0d0
        do io = 1, nocc
          do i_idx = 1, n
            if (.not. fp_row_owned(i_idx)) cycle
            energy_tmp = energy_tmp + occ_weight(io) * real(conjg(coef_frag_all(i_idx, io)) * tmp_all(i_idx, io))
          end do
          if (n_pw > 0) then
            do jo = 1, n_pw
              if (.not. pw_row_owned(jo)) cycle
              energy_tmp = energy_tmp + occ_weight(io) * real(conjg(coef_pw_all(jo, io)) * tmp_all(n+jo, io))
            end do
          end if
        end do
      else
        if (.not. allocated(dg_frag%momentum_blocks) .or. use_spatial_A) then
          call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
                     coef_frag_all(1:n, 1:nocc), n, (0.0d0, 0.0d0), tmp_mat, n)
          if (has_nonlocal) then
            call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_frag_all(1:n, 1:nocc), &
              tmp_mat)
          end if
        end if
        energy_tmp = 0.0d0
        do io = 1, nocc
          do i_idx = 1, n
            if (.not. fp_row_owned(i_idx)) cycle
            energy_tmp = energy_tmp + occ_weight(io) * real(conjg(coef_frag_all(i_idx, io)) * tmp_mat(i_idx, io))
          end do
        end do
      end if
        energy_local = energy_local + energy_tmp
      end do

      if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) then
        do ispin = 1, dg_frag%nspin
          nocc = min(dg_frag%nocc_spin(ispin), min(dg_frag%nstate_tot, n))
          if (nocc <= 0) cycle
          if (n_pw > 0) coef_pw_all(1:n_pw, 1:nocc) = dg_frag%coef_pw(1:n_pw, 1:nocc, ispin)
          pw_row_owned(:) = .false.
          do jo = 1, n_pw
            owner_rank_obs = dg_frag%id
            if (allocated(dg_frag%coef_pw_owner)) owner_rank_obs = dg_frag%coef_pw_owner(jo)
            if (owner_rank_obs == dg_frag%id) pw_row_owned(jo) = .true.
          end do
          do jo = 1, n_pw
            if (.not. pw_row_owned(jo)) cycle
            pw_weight_local = pw_weight_local + sum(abs(coef_pw_all(jo, 1:nocc))**2)
          end do
        end do
      end if
    if (do_interface_check) then
      allocate(dndt_frag(dg_frag%n_frag))
      dndt_frag = 0.0d0
      max_pair_residual = 0.0d0
      do ifrag = 1, dg_frag%n_frag
        do jfrag = 1, dg_frag%n_frag
          if (jfrag == ifrag) cycle
          dndt_frag(ifrag) = dndt_frag(ifrag) - interface_flow(ifrag, jfrag)
        end do
      end do

      do ifrag = 1, dg_frag%n_frag - 1
        do jfrag = ifrag + 1, dg_frag%n_frag
          pair_residual = interface_flow(ifrag, jfrag) + interface_flow(jfrag, ifrag)
          max_pair_residual = max(max_pair_residual, abs(pair_residual))
        end do
      end do
      charge_balance_residual = abs(sum(dndt_frag))

      deallocate(dndt_frag, interface_flow)
    end if


    if (allocated(Ap_mat)) deallocate(Ap_mat)
    if (allocated(A2_mat)) deallocate(A2_mat)
    if (allocated(op_mat)) deallocate(op_mat)
    deallocate(tmp_mat)
    if (allocated(occ_weight)) deallocate(occ_weight)
    if (allocated(coef_frag_all)) deallocate(coef_frag_all)
    if (allocated(coef_pw_all)) deallocate(coef_pw_all)
    if (allocated(coef_all)) deallocate(coef_all)
    if (allocated(tmp_all)) deallocate(tmp_all)
    if (allocated(fp_row_owned)) deallocate(fp_row_owned)
    if (allocated(pw_row_owned)) deallocate(pw_row_owned)
    ! Cache retained for reuse
    if (trace_obs) then
      call cpu_time(obs_t1)
      write(*,'(1x,a,1pe12.4)') '[DG-OBS] energy local done time=', obs_t1 - obs_t0
      flush(6)
    end if

  1000 continue
    
    ! MPI aggregation: sum local contributions from all ranks
    if (trace_obs) then
      write(*,'(1x,a)') '[DG-OBS] reductions start'
      flush(6)
      call cpu_time(obs_t0)
    end if
    if (do_energy_trace) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,1pe16.8,a,1pe16.8)') &
        "[DG-ENERGY-LOCAL] itt=", itt, " rank=", dg_frag%id, " id_frag=", dg_frag%id_frag, &
        " ifrag_start=", dg_frag%ifrag_start, " ifrag_end=", dg_frag%ifrag_end, &
        " isize_frag=", dg_frag%isize_frag, " energy_local=", energy_local, &
        " pw_weight_local=", pw_weight_local
      flush(6)
    end if
    call comm_summation(current_local, dg_frag%current, 3, dg_frag%icomm)
    scalar_reduce_in(1) = coef_occ_norm_local
    scalar_reduce_out(1) = 0.0d0
    call comm_summation(scalar_reduce_in, scalar_reduce_out, 1, dg_frag%icomm)
    coef_occ_norm_global = scalar_reduce_out(1)
    if (do_current_block_trace) then
      call comm_summation(current_block_local, current_block_sum, 3 * nblk_trace, dg_frag%icomm)
    end if
    scalar_reduce_in(1) = energy_local
    scalar_reduce_out(1) = 0.0d0
    call comm_summation(scalar_reduce_in, scalar_reduce_out, 1, dg_frag%icomm)
    dg_frag%total_energy = scalar_reduce_out(1)
    scalar_reduce_in(1) = pw_weight_local
    scalar_reduce_out(1) = 0.0d0
    call comm_summation(scalar_reduce_in, scalar_reduce_out, 1, dg_frag%icomm)
    dg_frag%pw_weight_raw = scalar_reduce_out(1)
    if (dg_frag%parallel_mode_orbital) then
      ! In orbital mode the observable contraction is row-owned: each rank
      ! contributes only its coefficient rows, so the world sum is already the
      ! full value and must not be averaged over subgroup ranks.
      frag_reduce_factor = 1.0d0
    else
      frag_reduce_factor = real(max(1, dg_frag%isize_frag), 8)
    end if
    energy_sum_raw = dg_frag%total_energy
    dg_frag%total_energy = dg_frag%total_energy / frag_reduce_factor
    if (do_energy_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a,i0,a,1pe16.8,a,1pe16.8,a,1pe16.8)') &
        "[DG-ENERGY-GLOBAL] itt=", itt, " energy_sum=", energy_sum_raw, &
        " frag_reduce=", frag_reduce_factor, " total_energy=", dg_frag%total_energy
      flush(6)
    end if

    ! Legacy paths can still carry replicated observable contributions.  The
    ! orbital path above is row-owned and the MPI sum is already complete.
    if (.not. dg_frag%parallel_mode_orbital) then
      dg_frag%current(:) = dg_frag%current(:) / real(max(1, dg_frag%isize), 8)
      dg_frag%coef_occ_norm = coef_occ_norm_global / real(max(1, dg_frag%isize), 8)
    else
      dg_frag%coef_occ_norm = coef_occ_norm_global
    end if
    if (itt == 1 .or. dg_frag%coef_occ_norm0 <= 0.0d0) dg_frag%coef_occ_norm0 = dg_frag%coef_occ_norm
    dg_frag%coef_occ_norm_drift = dg_frag%coef_occ_norm - dg_frag%coef_occ_norm0
    dg_frag%csc_occ_identity_norm = sqrt(max(0.0d0, csc_occ_identity_norm_global / real(max(1, dg_frag%isize), 8)))
    dg_frag%csc_occ_identity_max = csc_occ_identity_max_global
    dg_frag%occvirt_leakage = sqrt(max(0.0d0, occvirt_leakage_global / real(max(1, dg_frag%isize), 8)))
    dg_frag%rho_ff_elec = rho_ff_global / real(max(1, dg_frag%isize), 8)
    dg_frag%rho_fp_elec = rho_fp_global / real(max(1, dg_frag%isize), 8)
    dg_frag%rho_pp_elec = rho_pp_global / real(max(1, dg_frag%isize), 8)
    if (do_current_block_trace) then
      current_block_sum(:) = current_block_sum(:) / real(max(1, dg_frag%isize), 8)
      current_block_sum(:) = current_block_sum(:) / dble(system%ngrid)
      if (dg_frag%id == 0) then
        block_norm_min = huge(1.0d0)
        block_norm_max = 0.0d0
        n_nonzero_blocks = 0
        write(*,'(1x,a,i0,a,i0)') "        current-block-trace: itt=", itt, " nblocks=", nblk_trace
        do iblk_trace = 1, nblk_trace
          idx_base = iblk_trace
          block_norm = sqrt(current_block_sum(idx_base)**2 + &
            current_block_sum(nblk_trace + idx_base)**2 + current_block_sum(2 * nblk_trace + idx_base)**2)
          if (block_norm > 1.0d-30) then
            block_norm_min = min(block_norm_min, block_norm)
            block_norm_max = max(block_norm_max, block_norm)
            n_nonzero_blocks = n_nonzero_blocks + 1
          end if
          write(*,'(1x,a,i0,a,i0,a,i0,a,3(1pe14.6,1x),a,1pe14.6)') &
            "        current-block-trace: iblk=", iblk_trace, " ifrag_row=", dg_frag%momentum_blocks(iblk_trace)%ifrag_row, &
            " ifrag_col=", dg_frag%momentum_blocks(iblk_trace)%ifrag_col, " jdens=", &
            current_block_sum(idx_base), current_block_sum(nblk_trace + idx_base), current_block_sum(2 * nblk_trace + idx_base), &
            " norm=", block_norm
        end do
        if (n_nonzero_blocks > 1 .and. block_norm_min > 0.0d0) then
          write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
            "        current-block-trace summary: itt=", itt, " norm_min=", block_norm_min, &
            " norm_max=", block_norm_max, " max/min=", block_norm_max / block_norm_min
        else
          write(*,'(1x,a,i0,a)') "        current-block-trace summary: itt=", itt, " insufficient nonzero blocks"
        end if
        flush(6)
      end if
      deallocate(current_block_local, current_block_sum)
    end if
    if (.not. dg_frag%parallel_mode_orbital) then
      dg_frag%pw_weight_raw = dg_frag%pw_weight_raw / real(max(1, dg_frag%isize), 8)
    end if
    if (trace_obs) then
      call cpu_time(obs_t1)
      write(*,'(1x,a,1pe12.4)') '[DG-OBS] reductions done time=', obs_t1 - obs_t0
      flush(6)
    end if
    dg_frag%energy_kinetic = 0.0d0
    dg_frag%energy_nonlocal = 0.0d0
    dg_frag%energy_Ap = 0.0d0
    dg_frag%energy_A2 = 0.0d0

    ! Normalize by global grid count exactly as conventional calc_current().
    ! This avoids decomposition-dependent scaling from local/grid-view differences.
    dg_frag%current(:) = dg_frag%current(:) / dble(system%ngrid)

    dg_frag%current_para(:) = dg_frag%current(:)
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
      ! conventional calc_current() adds A * sum(|psi|^2) and normalizes by ngrid.
      ! Do not divide by hvol here, otherwise the diamagnetic part is too large.
      ne_density = nelec_ref / dble(system%ngrid)
    else
      ne_density = 0.0d0
    end if
    dg_frag%current_dia(:) = ne_density * rt%Ac_tot(:, itt)
    dg_frag%current_total(:) = dg_frag%current_para(:) + dg_frag%current_dia(:)
    if (.not. dg_frag%elec_num_baseline_ready) then
      dg_frag%elec_num_raw_t0 = dg_frag%elec_num_raw
      dg_frag%elec_num_scaled_t0 = dg_frag%elec_num_scaled
      dg_frag%elec_num_baseline_ready = .true.
    end if
    dg_frag%rho_drift_indicator = dg_frag%elec_num_raw - dg_frag%elec_num_raw_t0
    
    ! Store in rt structure for output (J_total = J_para + J_dia to match standard periodic SALMON)
    rt%curr(:, itt) = dg_frag%current_total(:)
    if (trace_obs) then
      write(*,'(1x,a)') '[DG-OBS] leave'
      flush(6)
    end if
    
  end subroutine calculate_observables

  subroutine calculate_observables_orbital_local(dg_frag, system, mg, ppg, rt, itt)
    use structures
    use communication, only: comm_summation
    use rt_dg_fragment_ops, only: apply_momentum_blocks, apply_matrix_blocks_batch, &
                                  apply_complex_matrix_blocks_batch, apply_nonlocal_pp_projector_batch, &
                                  ensure_nonlocal_pp_matrix_A
    use misc_routines, only: get_wtime
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    type(s_rt),             intent(inout) :: rt
    integer,                intent(in)    :: itt

    integer :: ispin, idir, io, iloc, nlocal, nocc
    integer :: env_status, env_len
    real(8) :: Ac_tot(3), A_squared, current_tmp, energy_tmp
    real(8) :: current_local(3), current_nl_local(3), current_nl_tmp
    real(8) :: coef_occ_norm_local, coef_occ_norm_global
    real(8) :: energy_local, energy_sum_raw
    real(8) :: nelec_ref, ne_density
    real(8) :: t0, t1
    real(8) :: scalar_reduce_in(1), scalar_reduce_out(1)
    real(8) :: nl_fd_delta, Ac_plus(3), Ac_minus(3)
    real(8), allocatable :: occ_weight(:)
    complex(8), allocatable :: tmp_h(:,:), tmp_p(:,:), tmp_nlp(:,:), tmp_nlm(:,:)
    complex(8) :: minus_i
    logical :: has_nonlocal, trace_obs, enable_nl_current_fd
    character(len=32) :: env_value
    real(8), parameter :: unit_dir(3,3) = reshape((/ &
      1.0d0, 0.0d0, 0.0d0, &
      0.0d0, 1.0d0, 0.0d0, &
      0.0d0, 0.0d0, 1.0d0 /), (/3, 3/))

    trace_obs = (itt == 1 .and. dg_frag%id == 0)
    if (trace_obs) then
      write(*,'(1x,a)') '[DG-OBS-LOCAL] enter'
      flush(6)
      t0 = get_wtime()
    end if

    nlocal = size(dg_frag%coef, 1)
    current_local(:) = 0.0d0
    current_nl_local(:) = 0.0d0
    coef_occ_norm_local = 0.0d0
    energy_local = 0.0d0
    dg_frag%pw_weight_raw = 0.0d0
    minus_i = cmplx(0.0d0, -1.0d0, kind=8)
    Ac_tot = rt%Ac_tot(:, itt)
    A_squared = Ac_tot(1)**2 + Ac_tot(2)**2 + Ac_tot(3)**2
    enable_nl_current_fd = .false.
    env_value = ''
    call get_environment_variable("SALMON_DG_NL_CURRENT_FD", env_value, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_value(1:1) == '1' .or. env_value(1:1) == 'y' .or. env_value(1:1) == 'Y' .or. &
          env_value(1:1) == 't' .or. env_value(1:1) == 'T') then
        enable_nl_current_fd = .true.
      end if
    end if
    nl_fd_delta = 1.0d-5

    call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_tot)
    has_nonlocal = dg_frag%has_nl_cache

    do ispin = 1, dg_frag%nspin
      nocc = min(dg_frag%nocc_spin(ispin), min(dg_frag%nstate_tot, size(dg_frag%coef, 2)))
      if (nocc <= 0) cycle
      allocate(occ_weight(nocc))
      occ_weight(:) = 0.0d0
      call get_dg_spin_occ_info(dg_frag, system, ispin, occ_weight, nocc)
      allocate(tmp_h(nlocal, nocc), tmp_p(nlocal, nocc))
      if (has_nonlocal .and. enable_nl_current_fd) then
        allocate(tmp_nlp(nlocal, nocc), tmp_nlm(nlocal, nocc))
      end if

      tmp_h(:, :) = (0.0d0, 0.0d0)
      do idir = 1, 3
        tmp_p(:, :) = (0.0d0, 0.0d0)
        call apply_momentum_blocks(dg_frag, ispin, unit_dir(:, idir), &
                                   dg_frag%coef(1:nlocal, 1:nocc, ispin), tmp_p)
        current_tmp = 0.0d0
        do io = 1, nocc
          do iloc = 1, nlocal
            current_tmp = current_tmp + occ_weight(io) * &
              aimag(conjg(dg_frag%coef(iloc, io, ispin)) * tmp_p(iloc, io))
          end do
        end do
        current_local(idir) = current_local(idir) + current_tmp
        if (has_nonlocal .and. enable_nl_current_fd) then
          Ac_plus(:) = Ac_tot(:)
          Ac_minus(:) = Ac_tot(:)
          Ac_plus(idir) = Ac_plus(idir) + nl_fd_delta
          Ac_minus(idir) = Ac_minus(idir) - nl_fd_delta
          tmp_nlp(:, :) = (0.0d0, 0.0d0)
          tmp_nlm(:, :) = (0.0d0, 0.0d0)
          call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_plus)
          if (allocated(dg_frag%H_nl_local_block_ids)) then
            call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, &
                                                   dg_frag%coef(1:nlocal, 1:nocc, ispin), tmp_nlp, &
                                                   dg_frag%H_nl_local_block_ids)
          else
            call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, &
                                                   dg_frag%coef(1:nlocal, 1:nocc, ispin), tmp_nlp)
          end if
          call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_minus)
          if (allocated(dg_frag%H_nl_local_block_ids)) then
            call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, &
                                                   dg_frag%coef(1:nlocal, 1:nocc, ispin), tmp_nlm, &
                                                   dg_frag%H_nl_local_block_ids)
          else
            call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, &
                                                   dg_frag%coef(1:nlocal, 1:nocc, ispin), tmp_nlm)
          end if
          current_nl_tmp = 0.0d0
          do io = 1, nocc
            do iloc = 1, nlocal
              current_nl_tmp = current_nl_tmp + occ_weight(io) * real( &
                conjg(dg_frag%coef(iloc, io, ispin)) * (tmp_nlp(iloc, io) - tmp_nlm(iloc, io)), kind=8)
            end do
          end do
          current_nl_tmp = current_nl_tmp / (2.0d0 * nl_fd_delta)
          current_local(idir) = current_local(idir) + current_nl_tmp
          current_nl_local(idir) = current_nl_local(idir) + current_nl_tmp
          call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_tot)
        end if
        if (abs(Ac_tot(idir)) > 1.0d-30) then
          tmp_h(:, :) = tmp_h(:, :) + minus_i * Ac_tot(idir) * tmp_p(:, :)
        end if
      end do

      if (allocated(dg_frag%H_mat_blocks)) then
        if (allocated(dg_frag%H_local_block_ids)) then
          call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, &
                                         dg_frag%coef(1:nlocal, 1:nocc, ispin), tmp_h, dg_frag%H_local_block_ids)
        else
          call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, &
                                         dg_frag%coef(1:nlocal, 1:nocc, ispin), tmp_h)
        end if
      end if
      if (has_nonlocal) then
        if (allocated(dg_frag%H_nl_blocks)) then
          if (allocated(dg_frag%H_nl_local_block_ids)) then
            call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, &
                                                   dg_frag%coef(1:nlocal, 1:nocc, ispin), tmp_h, &
                                                   dg_frag%H_nl_local_block_ids)
          else
            call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, &
                                                   dg_frag%coef(1:nlocal, 1:nocc, ispin), tmp_h)
          end if
        else
          call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, &
                                                 dg_frag%coef(1:nlocal, 1:nocc, ispin), tmp_h)
        end if
      end if
      tmp_h(:, :) = tmp_h(:, :) + 0.5d0 * A_squared * dg_frag%coef(1:nlocal, 1:nocc, ispin)

      energy_tmp = 0.0d0
      do io = 1, nocc
        do iloc = 1, nlocal
          energy_tmp = energy_tmp + occ_weight(io) * real(conjg(dg_frag%coef(iloc, io, ispin)) * tmp_h(iloc, io), kind=8)
        end do
      end do
      energy_local = energy_local + energy_tmp
      coef_occ_norm_local = coef_occ_norm_local + sum(abs(dg_frag%coef(1:nlocal, 1:nocc, ispin))**2)

      if (allocated(tmp_nlp)) deallocate(tmp_nlp)
      if (allocated(tmp_nlm)) deallocate(tmp_nlm)
      deallocate(tmp_h, tmp_p, occ_weight)
    end do

    if (trace_obs) then
      t1 = get_wtime()
      write(*,'(1x,a,1pe12.4)') '[DG-OBS-LOCAL] local contractions done time=', t1 - t0
      flush(6)
      t0 = get_wtime()
      write(*,'(1x,a)') '[DG-OBS-LOCAL] reductions start'
      flush(6)
    end if
    call comm_summation(current_local, dg_frag%current, 3, dg_frag%icomm)
    call comm_summation(current_nl_local, dg_frag%current_nl, 3, dg_frag%icomm)
    scalar_reduce_in(1) = coef_occ_norm_local
    scalar_reduce_out(1) = 0.0d0
    call comm_summation(scalar_reduce_in, scalar_reduce_out, 1, dg_frag%icomm)
    coef_occ_norm_global = scalar_reduce_out(1)
    scalar_reduce_in(1) = energy_local
    scalar_reduce_out(1) = 0.0d0
    call comm_summation(scalar_reduce_in, scalar_reduce_out, 1, dg_frag%icomm)
    dg_frag%total_energy = scalar_reduce_out(1)
    if (trace_obs) then
      t1 = get_wtime()
      write(*,'(1x,a,1pe12.4)') '[DG-OBS-LOCAL] reductions done time=', t1 - t0
      flush(6)
    end if

    energy_sum_raw = dg_frag%total_energy
    dg_frag%total_energy = energy_sum_raw
    dg_frag%coef_occ_norm = coef_occ_norm_global
    if (itt == 1 .or. dg_frag%coef_occ_norm0 <= 0.0d0) dg_frag%coef_occ_norm0 = dg_frag%coef_occ_norm
    dg_frag%coef_occ_norm_drift = dg_frag%coef_occ_norm - dg_frag%coef_occ_norm0
    dg_frag%csc_occ_identity_norm = 0.0d0
    dg_frag%csc_occ_identity_max = 0.0d0
    dg_frag%occvirt_leakage = 0.0d0
    dg_frag%rho_ff_elec = 0.0d0
    dg_frag%rho_fp_elec = 0.0d0
    dg_frag%rho_pp_elec = 0.0d0
    dg_frag%energy_kinetic = 0.0d0
    dg_frag%energy_nonlocal = 0.0d0
    dg_frag%energy_Ap = 0.0d0
    dg_frag%energy_A2 = 0.0d0

    dg_frag%current(:) = dg_frag%current(:) / dble(system%ngrid)
    dg_frag%current_nl(:) = dg_frag%current_nl(:) / dble(system%ngrid)
    dg_frag%current_para(:) = dg_frag%current(:)
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
      ! conventional calc_current() adds A * sum(|psi|^2) and normalizes by ngrid.
      ! Do not divide by hvol here, otherwise the diamagnetic part is too large.
      ne_density = nelec_ref / dble(system%ngrid)
    else
      ne_density = 0.0d0
    end if
    dg_frag%current_dia(:) = ne_density * rt%Ac_tot(:, itt)
    dg_frag%current_total(:) = dg_frag%current_para(:) + dg_frag%current_dia(:)
    if (.not. dg_frag%elec_num_baseline_ready) then
      dg_frag%elec_num_raw_t0 = dg_frag%elec_num_raw
      dg_frag%elec_num_scaled_t0 = dg_frag%elec_num_scaled
      dg_frag%elec_num_baseline_ready = .true.
    end if
    dg_frag%rho_drift_indicator = dg_frag%elec_num_raw - dg_frag%elec_num_raw_t0
    rt%curr(:, itt) = dg_frag%current_total(:)

    if (trace_obs) then
      write(*,'(1x,a)') '[DG-OBS-LOCAL] leave'
      flush(6)
    end if
  end subroutine calculate_observables_orbital_local

  integer function find_matrix_block_local(block_map, ifrag_row, ifrag_col) result(iblk)
    implicit none
    integer, intent(in) :: block_map(:, :)
    integer, intent(in) :: ifrag_row, ifrag_col

    iblk = 0
    if (ifrag_row < 1 .or. ifrag_row > size(block_map, 1)) return
    if (ifrag_col < 1 .or. ifrag_col > size(block_map, 2)) return
    iblk = block_map(ifrag_row, ifrag_col)
  end function find_matrix_block_local

  subroutine build_total_potential_grid_local(grid, Vh, Vxc_spin, Vpsl, V_total)
    use structures
    implicit none
    type(s_rgrid), intent(in) :: grid
    type(s_scalar), intent(in) :: Vh, Vxc_spin, Vpsl
    real(8), intent(out) :: V_total(grid%is(1):grid%ie(1), grid%is(2):grid%ie(2), grid%is(3):grid%ie(3))
    integer :: ix, iy, iz

!$omp parallel do collapse(2) private(ix,iy,iz)
    do iz = grid%is(3), grid%ie(3)
      do iy = grid%is(2), grid%ie(2)
!$omp simd
        do ix = grid%is(1), grid%ie(1)
          V_total(ix, iy, iz) = Vpsl%f(ix, iy, iz) + Vh%f(ix, iy, iz) + Vxc_spin%f(ix, iy, iz)
        end do
      end do
    end do
!$omp end parallel do
  end subroutine build_total_potential_grid_local

  integer function map_global_to_phi_box_coord_obs(ig, lb, ub, lgtot) result(iloc)
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
  end function map_global_to_phi_box_coord_obs

  subroutine build_local_potential_applied_basis_local(dg_frag, ifrag, i_local, jo, mg, V_total, V_phi)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, i_local, jo
    type(s_rgrid), intent(in) :: mg
    real(8), intent(in) :: V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    complex(8), intent(out) :: V_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    integer :: iorg(3), ndom(3), loc_s(3), loc_e(3)
    integer :: lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi
    integer :: g_s(3), g_e(3), ov_s(3), ov_e(3)
    integer :: lx, ly, lz, gx, gy, gz, bx, by, bz
    integer :: p_lb1, p_ub1, p_lb2, p_ub2, p_lb3, p_ub3
    complex(8) :: phi0
    logical :: has_overlap

    V_phi(:, :, :) = (0.0d0, 0.0d0)
    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    g_s(:) = iorg(:)
    g_e(:) = iorg(:) + ndom(:) - 1
    ov_s(:) = max(g_s(:), mg%is(:))
    ov_e(:) = min(g_e(:), mg%ie(:))
    has_overlap = all(ov_s(:) <= ov_e(:))
    if (.not. has_overlap) return
    loc_s(:) = ov_s(:) - iorg(:) + 1
    loc_e(:) = ov_e(:) - iorg(:) + 1
    lx_lo = loc_s(1)
    lx_hi = loc_e(1)
    ly_lo = loc_s(2)
    ly_hi = loc_e(2)
    lz_lo = loc_s(3)
    lz_hi = loc_e(3)
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub3 = ubound(dg_frag%phi_frag, 3)

    if (allocated(dg_frag%phi_frag_c)) then
!$omp parallel do private(lz, ly, lx, gx, gy, gz, bx, by, bz, phi0) schedule(static)
      do lz = lz_lo, lz_hi
        gz = ov_s(3) + (lz - lz_lo)
        do ly = ly_lo, ly_hi
          gy = ov_s(2) + (ly - ly_lo)
          do lx = lx_lo, lx_hi
            gx = ov_s(1) + (lx - lx_lo)
            bx = map_global_to_phi_box_coord_obs(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
            by = map_global_to_phi_box_coord_obs(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
            bz = map_global_to_phi_box_coord_obs(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
            if (bx == 0 .or. by == 0 .or. bz == 0) cycle
            phi0 = dg_frag%phi_frag_c(bx, by, bz, jo, i_local)
            V_phi(gx, gy, gz) = V_total(gx, gy, gz) * phi0
          end do
        end do
      end do
!$omp end parallel do
    else
!$omp parallel do private(lz, ly, lx, gx, gy, gz, bx, by, bz, phi0) schedule(static)
      do lz = lz_lo, lz_hi
        gz = ov_s(3) + (lz - lz_lo)
        do ly = ly_lo, ly_hi
          gy = ov_s(2) + (ly - ly_lo)
          do lx = lx_lo, lx_hi
            gx = ov_s(1) + (lx - lx_lo)
            bx = map_global_to_phi_box_coord_obs(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
            by = map_global_to_phi_box_coord_obs(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
            bz = map_global_to_phi_box_coord_obs(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
            if (bx == 0 .or. by == 0 .or. bz == 0) cycle
            phi0 = cmplx(dg_frag%phi_frag(bx, by, bz, jo, i_local), 0.0d0, kind=8)
            V_phi(gx, gy, gz) = V_total(gx, gy, gz) * phi0
          end do
        end do
      end do
!$omp end parallel do
    end if
  end subroutine build_local_potential_applied_basis_local

  subroutine integrate_local_basis_with_field_local(dg_frag, ifrag, i_local, io, mg, field, hvol, integral)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, i_local, io
    type(s_rgrid), intent(in) :: mg
    complex(8), intent(in) :: field(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8), intent(in) :: hvol
    complex(8), intent(out) :: integral
    integer :: iorg(3), ndom(3), loc_s(3), loc_e(3)
    integer :: lx_lo, lx_hi, ly_lo, ly_hi, lz_lo, lz_hi
    integer :: g_s(3), g_e(3), ov_s(3), ov_e(3)
    integer :: lx, ly, lz, gx, gy, gz, bx, by, bz
    integer :: p_lb1, p_ub1, p_lb2, p_ub2, p_lb3, p_ub3
    logical :: has_overlap

    integral = (0.0d0, 0.0d0)
    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    g_s(:) = iorg(:)
    g_e(:) = iorg(:) + ndom(:) - 1
    ov_s(:) = max(g_s(:), mg%is(:))
    ov_e(:) = min(g_e(:), mg%ie(:))
    has_overlap = all(ov_s(:) <= ov_e(:))
    if (.not. has_overlap) return
    loc_s(:) = ov_s(:) - iorg(:) + 1
    loc_e(:) = ov_e(:) - iorg(:) + 1
    lx_lo = loc_s(1)
    lx_hi = loc_e(1)
    ly_lo = loc_s(2)
    ly_hi = loc_e(2)
    lz_lo = loc_s(3)
    lz_hi = loc_e(3)
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub3 = ubound(dg_frag%phi_frag, 3)

    if (allocated(dg_frag%phi_frag_c)) then
      do lz = lz_lo, lz_hi
        gz = ov_s(3) + (lz - lz_lo)
        do ly = ly_lo, ly_hi
          gy = ov_s(2) + (ly - ly_lo)
          do lx = lx_lo, lx_hi
            gx = ov_s(1) + (lx - lx_lo)
            bx = map_global_to_phi_box_coord_obs(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
            by = map_global_to_phi_box_coord_obs(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
            bz = map_global_to_phi_box_coord_obs(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
            if (bx == 0 .or. by == 0 .or. bz == 0) cycle
            integral = integral + conjg(dg_frag%phi_frag_c(bx, by, bz, io, i_local)) * field(gx, gy, gz) * hvol
          end do
        end do
      end do
    else
      do lz = lz_lo, lz_hi
        gz = ov_s(3) + (lz - lz_lo)
        do ly = ly_lo, ly_hi
          gy = ov_s(2) + (ly - ly_lo)
          do lx = lx_lo, lx_hi
            gx = ov_s(1) + (lx - lx_lo)
            bx = map_global_to_phi_box_coord_obs(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
            by = map_global_to_phi_box_coord_obs(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
            bz = map_global_to_phi_box_coord_obs(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
            if (bx == 0 .or. by == 0 .or. bz == 0) cycle
            integral = integral + cmplx(dg_frag%phi_frag(bx, by, bz, io, i_local), 0.0d0, kind=8) * field(gx, gy, gz) * hvol
          end do
        end do
      end do
    end if
  end subroutine integrate_local_basis_with_field_local
