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
    logical :: enable_mij_audit, do_mij_audit, has_esp, enable_mij_block_audit
    real(8), allocatable :: interface_flow(:,:), dndt_frag(:)
    real(8) :: pair_residual, max_pair_residual, charge_balance_residual
    real(8) :: current_tmp, energy_tmp, pw_weight_local, kpw_dir
    real(8) :: jpara_state, jpara_top_abs
    real(8) :: Ac_tot(3), A_squared
    real(8) :: current_local(3), energy_local
    real(8) :: energy_static_local, energy_kin_local, energy_nl_local, energy_ap_local, energy_a2_local
    real(8) :: energy_static_sum, energy_kin_sum, energy_nl_sum, energy_ap_sum, energy_a2_sum
    real(8) :: energy_kin_rs_sum, energy_one_rs_sum
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
    complex(8), allocatable :: tmp_probe(:,:), dense_probe_mat(:,:), dense_probe_out(:,:)
    complex(8), allocatable :: coef_occ(:,:), coef_emp(:,:), tmp_emp(:,:), mij_mat(:,:)
    real(8), allocatable :: current_block_local(:), current_block_sum(:)
    real(8), allocatable :: top_abs2(:), top_de(:)
    integer, allocatable :: top_spin(:), top_occ(:), top_emp(:)
    logical :: has_nonlocal, use_hmat_complex, use_mixed_current
    logical :: have_occvirt_ref
    logical, save :: occvirt_ref_mode_initialized = .false.
    logical, save :: occvirt_ref_legacy_mode = .false.
    logical :: require_dense_nl
    logical :: use_spatial_A
    logical :: enable_obs_charge_check
    integer :: obs_charge_check_stride
    real(8) :: obs_charge_check_tol
    real(8) :: obs_charge_local, obs_charge_global, obs_charge_diff
    integer :: obs_charge_check_env_len, obs_charge_check_env_status
    logical, parameter :: enable_energy_component_probe = .false.
    real(8), allocatable :: Ap_mat(:,:), A2_mat(:,:)
    integer, allocatable :: diag_block_ids(:)
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
    energy_kin_rs_sum = 0.0d0
    energy_one_rs_sum = 0.0d0

    n = dg_frag%n_mat_max
    use_spatial_A = (trim(theory) == 'single_scale_maxwell_tddft' .and. allocated(system%Ac_micro%v) .and. dg_frag%has_real_space_basis)
    enable_obs_charge_check = .false.
    obs_charge_check_stride = 10
    obs_charge_check_tol = 1.0d-8
    obs_charge_local = 0.0d0
    obs_charge_global = 0.0d0
    obs_charge_diff = 0.0d0
    env_value = ''
    call get_environment_variable("SALMON_DG_OBS_CHARGE_CHECK", env_value, length=obs_charge_check_env_len, status=obs_charge_check_env_status)
    if (obs_charge_check_env_status == 0 .and. obs_charge_check_env_len > 0) then
      if (env_value(1:1) == '1' .or. env_value(1:1) == 'y' .or. env_value(1:1) == 'Y' .or. &
          env_value(1:1) == 't' .or. env_value(1:1) == 'T') then
        enable_obs_charge_check = .true.
      end if
    end if
    env_value = ''
    call get_environment_variable("SALMON_DG_OBS_CHARGE_CHECK_STRIDE", env_value, length=obs_charge_check_env_len, status=obs_charge_check_env_status)
    if (obs_charge_check_env_status == 0 .and. obs_charge_check_env_len > 0) then
      read(env_value(1:obs_charge_check_env_len), *, iostat=parse_status) obs_charge_check_stride
      if (parse_status /= 0) obs_charge_check_stride = 10
    end if
    if (obs_charge_check_stride < 1) obs_charge_check_stride = 1
    env_value = ''
    call get_environment_variable("SALMON_DG_OBS_CHARGE_CHECK_TOL", env_value, length=obs_charge_check_env_len, status=obs_charge_check_env_status)
    if (obs_charge_check_env_status == 0 .and. obs_charge_check_env_len > 0) then
      read(env_value(1:obs_charge_check_env_len), *, iostat=parse_status) obs_charge_check_tol
      if (parse_status /= 0) obs_charge_check_tol = 1.0d-8
    end if
    enable_current_block_trace = .false.
    current_block_trace_stride = 20
    current_block_trace_maxblocks = 0
    env_value = ''
    call get_environment_variable("SALMON_DG_CURRENT_BLOCK_TRACE", env_value, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_value(1:1) == '1' .or. env_value(1:1) == 'y' .or. env_value(1:1) == 'Y' .or. &
          env_value(1:1) == 't' .or. env_value(1:1) == 'T') then
        enable_current_block_trace = .true.
      end if
    end if
    env_value = ''
    call get_environment_variable("SALMON_DG_CURRENT_BLOCK_TRACE_STRIDE", env_value, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_value(1:env_len), *, iostat=parse_status) current_block_trace_stride
      if (parse_status /= 0) current_block_trace_stride = 20
    end if
    if (current_block_trace_stride < 1) current_block_trace_stride = 1
    env_value = ''
    call get_environment_variable("SALMON_DG_CURRENT_BLOCK_TRACE_MAXBLOCKS", env_value, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_value(1:env_len), *, iostat=parse_status) current_block_trace_maxblocks
      if (parse_status /= 0) current_block_trace_maxblocks = 0
    end if
    if (current_block_trace_maxblocks < 0) current_block_trace_maxblocks = 0
    do_current_block_trace = enable_current_block_trace .and. allocated(dg_frag%momentum_blocks) .and. &
      (itt == 1 .or. mod(itt, current_block_trace_stride) == 0)

    if (.not. occvirt_ref_mode_initialized) then
      env_value = ''
      call get_environment_variable('SALMON_DG_OCCVIRT_REF_MODE', env_value, length=env_len, status=env_status)
      occvirt_ref_legacy_mode = .false.
      if (env_status == 0 .and. env_len > 0) then
        select case (env_value(1:1))
        case ('l','L','0')
          occvirt_ref_legacy_mode = .true.
        case default
          occvirt_ref_legacy_mode = .false.
        end select
      end if
      occvirt_ref_mode_initialized = .true.
    end if

    enable_mij_audit = .false.
    enable_mij_block_audit = .false.
    mij_audit_stride = 50
    mij_audit_topk = 10
    mij_audit_max_occ = 0
    mij_audit_max_emp = 0
    mij_audit_dir = 3
    mij_audit_ewin = 0.0d0
    env_value = ''
    call get_environment_variable("SALMON_DG_MIJ_AUDIT", env_value, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_value(1:1) == '1' .or. env_value(1:1) == 'y' .or. env_value(1:1) == 'Y' .or. &
          env_value(1:1) == 't' .or. env_value(1:1) == 'T') then
        enable_mij_audit = .true.
      end if
    end if
    env_value = ''
    call get_environment_variable("SALMON_DG_MIJ_BLOCK_AUDIT", env_value, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_value(1:1) == '1' .or. env_value(1:1) == 'y' .or. env_value(1:1) == 'Y' .or. &
          env_value(1:1) == 't' .or. env_value(1:1) == 'T') then
        enable_mij_block_audit = .true.
      end if
    end if
    env_value = ''
    call get_environment_variable("SALMON_DG_MIJ_AUDIT_STRIDE", env_value, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_value(1:env_len), *, iostat=parse_status) mij_audit_stride
      if (parse_status /= 0) mij_audit_stride = 50
    end if
    if (mij_audit_stride < 1) mij_audit_stride = 1
    env_value = ''
    call get_environment_variable("SALMON_DG_MIJ_AUDIT_TOPK", env_value, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_value(1:env_len), *, iostat=parse_status) mij_audit_topk
      if (parse_status /= 0) mij_audit_topk = 10
    end if
    if (mij_audit_topk < 1) mij_audit_topk = 1
    env_value = ''
    call get_environment_variable("SALMON_DG_MIJ_AUDIT_MAX_OCC", env_value, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_value(1:env_len), *, iostat=parse_status) mij_audit_max_occ
      if (parse_status /= 0) mij_audit_max_occ = 0
    end if
    if (mij_audit_max_occ < 0) mij_audit_max_occ = 0
    env_value = ''
    call get_environment_variable("SALMON_DG_MIJ_AUDIT_MAX_EMP", env_value, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_value(1:env_len), *, iostat=parse_status) mij_audit_max_emp
      if (parse_status /= 0) mij_audit_max_emp = 0
    end if
    if (mij_audit_max_emp < 0) mij_audit_max_emp = 0
    env_value = ''
    call get_environment_variable("SALMON_DG_MIJ_AUDIT_DIR", env_value, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_value(1:env_len), *, iostat=parse_status) mij_audit_dir
      if (parse_status /= 0) mij_audit_dir = 3
    end if
    if (mij_audit_dir < 1 .or. mij_audit_dir > 3) mij_audit_dir = 3
    env_value = ''
    call get_environment_variable("SALMON_DG_MIJ_AUDIT_EWIN", env_value, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_value(1:env_len), *, iostat=parse_status) mij_audit_ewin
      if (parse_status /= 0) mij_audit_ewin = 0.0d0
    end if
    if (mij_audit_ewin < 0.0d0) mij_audit_ewin = 0.0d0

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
        call comm_summation(obs_charge_local, obs_charge_global, dg_frag%icomm)
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
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n + n_pw
    nstate_use_diag = dg_frag%nstate_tot
    if (allocated(dg_frag%coef_ref_all)) then
      if (size(dg_frag%coef_ref_all, 1) /= n_tot .or. size(dg_frag%coef_ref_all, 2) /= nstate_use_diag .or. &
          size(dg_frag%coef_ref_all, 3) /= dg_frag%nspin) then
        deallocate(dg_frag%coef_ref_all)
        dg_frag%coef_ref_ready = .false.
      end if
    end if
    if (.not. dg_frag%coef_ref_ready .and. (occvirt_ref_legacy_mode .or. itt == 1)) then
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
    have_occvirt_ref = dg_frag%coef_ref_ready
    max_nocc = max(1, maxval(dg_frag%nocc_spin(1:dg_frag%nspin)))

    allocate(tmp_mat(n, max_nocc))
    if (enable_energy_component_probe) allocate(tmp_probe(n, max_nocc))
    allocate(coef_frag_all(n, max_nocc))
    allocate(occ_weight(max_nocc))
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
    call timer_begin(LOG_CALC_CURRENT)
    do ispin = 1, dg_frag%nspin
      nocc = min(dg_frag%nocc_spin(ispin), min(dg_frag%nstate_tot, n))
      if (nocc <= 0) cycle
      use_mixed_current = (n_pw > 0 .and. mixed_fp_coupling_active(dg_frag, ispin))
      coef_frag_all(1:n, 1:nocc) = dg_frag%coef(1:n, 1:nocc, ispin)
      if (n_pw > 0) then
        coef_pw_all(1:n_pw, 1:nocc) = dg_frag%coef_pw(1:n_pw, 1:nocc, ispin)
      end if
      if (n_pw > 0) then
        coef_all(1:n_tot, 1:nocc) = (0.0d0, 0.0d0)
        coef_all(1:n, 1:nocc) = coef_frag_all(1:n, 1:nocc)
        coef_all(n+1:n_tot, 1:nocc) = coef_pw_all(1:n_pw, 1:nocc)
      end if
      coef_occ_norm_local = coef_occ_norm_local + sum(abs(coef_frag_all(1:n, 1:nocc))**2)
      if (n_pw > 0) coef_occ_norm_local = coef_occ_norm_local + sum(abs(coef_pw_all(1:n_pw, 1:nocc))**2)
      do idir = 1, 3
        ! momentum_mat = <φ|∇|φ>, need to apply -i via aimag() and include factor 2
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
          current_tmp = sum(aimag(conjg(coef_all(1:n_tot, 1:nocc)) * tmp_all(1:n_tot, 1:nocc)))
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
          ! Factor -2.0: -1 for operator sign convention, 2 for Im[ψ*∇ψ] normalization
          current_tmp = sum(aimag(conjg(coef_frag_all(1:n, 1:nocc)) * tmp_mat(1:n, 1:nocc)))
          if (n_pw > 0) then
            do jo = 1, n_pw
              kpw_dir = dg_frag%k_pw(idir, jo)
              if (abs(kpw_dir) < 1.0d-15) cycle
              current_tmp = current_tmp + kpw_dir * sum(abs(coef_pw_all(jo, 1:nocc))**2)
            end do
          end if
        end if

        if (idir == 3) then
          if (use_mixed_current) then
            do io_state = 1, nocc
              jpara_state = -2.0d0 * sum(aimag(conjg(coef_all(1:n_tot, io_state)) * tmp_all(1:n_tot, io_state)))
              if (abs(jpara_state) > jpara_top_abs) then
                jpara_top_abs = abs(jpara_state)
                dg_frag%jpara_top_occ_state = io_state
                dg_frag%jpara_top_occ_value = jpara_state
              end if
            end do
          else
            do io_state = 1, nocc
              jpara_state = -2.0d0 * sum(aimag(conjg(coef_frag_all(1:n, io_state)) * tmp_mat(1:n, io_state)))
              if (abs(jpara_state) > jpara_top_abs) then
                jpara_top_abs = abs(jpara_state)
                dg_frag%jpara_top_occ_state = io_state
                dg_frag%jpara_top_occ_value = jpara_state
              end if
            end do
          end if
        end if
        current_local(idir) = current_local(idir) - 2.0d0 * current_tmp

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
                mfp = (0.0d0, 0.0d0)
                do jj_trace = 1, ncol_trace
                  ig_j_trace = dg_frag%index_basis(jj_trace, ifrag_col_trace, ispin)
                  if (ig_j_trace < 1 .or. ig_j_trace > n) cycle
                  mfp = mfp + cmplx(dg_frag%momentum_blocks(iblk_trace)%val(idir, ii_trace, jj_trace, ispin), 0.0d0, kind=8) * &
                    coef_frag_all(ig_j_trace, io)
                end do
                block_tmp = block_tmp + aimag(conjg(coef_frag_all(ig_i_trace, io)) * mfp)
              end do
            end do
            idx_trace = (idir - 1) * nblk_trace + iblk_trace
            current_block_local(idx_trace) = current_block_local(idx_trace) - 2.0d0 * block_tmp
          end do
        end if
      end do
    end do
    call timer_end(LOG_CALC_CURRENT)

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

    call comm_summation(csc_occ_identity_norm_local, csc_occ_identity_norm_global, dg_frag%icomm)
    csc_occ_identity_max_in(1) = csc_occ_identity_max_local
    call comm_get_max(csc_occ_identity_max_in, csc_occ_identity_max_out, 1, dg_frag%icomm)
    csc_occ_identity_max_global = csc_occ_identity_max_out(1)
    call comm_summation(occvirt_leakage_local, occvirt_leakage_global, dg_frag%icomm)
    call comm_summation(rho_ff_local, rho_ff_global, dg_frag%icomm)
    call comm_summation(rho_fp_local, rho_fp_global, dg_frag%icomm)
    call comm_summation(rho_pp_local, rho_pp_global, dg_frag%icomm)

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
      Ac_tot = rt%Ac_tot(:, itt)
      A_squared = Ac_tot(1)**2 + Ac_tot(2)**2 + Ac_tot(3)**2
      
      require_dense_nl = (.not. allocated(dg_frag%H_mat_blocks)) .or. &
                         (allocated(dg_frag%H_mat_c) .and. allocated(dg_frag%phi_frag_c)) .or. &
                         use_spatial_A .or. do_interface_check
      call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_tot, require_dense_nl)
      has_nonlocal = dg_frag%has_nl_cache

      ! Calculate total energy: E = <ψ|H(t)|ψ>
      ! H(t) = H_0 - i*A(t)·∇ + A²(t)/2 + V_NL(A)
      do ispin = 1, dg_frag%nspin
      nocc = min(dg_frag%nocc_spin(ispin), min(dg_frag%nstate_tot, n))
      if (nocc <= 0) cycle
      occ_weight(:) = 0.0d0
      call get_dg_spin_occ_info(dg_frag, system, ispin, occ_weight, nocc)
      coef_frag_all(1:n, 1:nocc) = dg_frag%coef(1:n, 1:nocc, ispin)
      if (n_pw > 0) then
        coef_pw_all(1:n_pw, 1:nocc) = dg_frag%coef_pw(1:n_pw, 1:nocc, ispin)
      end if
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
        if (has_nonlocal .and. allocated(dg_frag%H_nl_cache) .and. ((.not. allocated(dg_frag%H_mat_blocks)) .or. use_hmat_complex)) then
          op_mat(:, :) = op_mat(:, :) + dg_frag%H_nl_cache(1:n, 1:n, ispin)
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
            if (allocated(dg_frag%H_nl_cache) .and. ((.not. allocated(dg_frag%H_mat_blocks)) .or. use_hmat_complex)) then
              tmp_all(1:n, 1:nocc) = tmp_all(1:n, 1:nocc) + &
                matmul(dg_frag%H_nl_cache(1:n, 1:n, ispin), coef_all(1:n, 1:nocc))
            else
              call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all(1:n, 1:nocc), &
                tmp_all(1:n, 1:nocc))
            end if
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
            if (enable_energy_component_probe) then
              tmp_probe(:, :) = (0.0d0, 0.0d0)
              if (allocated(dg_frag%H_mat_kinetic_blocks)) then
                call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_kinetic_blocks, ispin, coef_frag_all(1:n, 1:nocc), tmp_probe)
                if (itt == 1 .or. itt == 40) then
                  n_diag_block_ids = 0
                  do iblk = 1, size(dg_frag%H_mat_kinetic_blocks)
                    if (dg_frag%H_mat_kinetic_blocks(iblk)%ifrag_row /= dg_frag%H_mat_kinetic_blocks(iblk)%ifrag_col) cycle
                    n_diag_block_ids = n_diag_block_ids + 1
                  end do
                  if (allocated(diag_block_ids)) deallocate(diag_block_ids)
                  if (n_diag_block_ids > 0) then
                    allocate(diag_block_ids(n_diag_block_ids))
                    idb = 0
                    do iblk = 1, size(dg_frag%H_mat_kinetic_blocks)
                      if (dg_frag%H_mat_kinetic_blocks(iblk)%ifrag_row /= dg_frag%H_mat_kinetic_blocks(iblk)%ifrag_col) cycle
                      idb = idb + 1
                      diag_block_ids(idb) = iblk
                    end do
                    if (.not. allocated(dense_probe_out)) allocate(dense_probe_out(n, nocc))
                    dense_probe_out(:, :) = (0.0d0, 0.0d0)
                    call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_kinetic_blocks, ispin, coef_frag_all(1:n, 1:nocc), &
                      dense_probe_out, diag_block_ids)
                    do io = 1, nocc
                      energy_kin_diag_local = energy_kin_diag_local + occ_weight(io) * &
                        sum(real(conjg(coef_frag_all(1:n, io)) * dense_probe_out(1:n, io)))
                      energy_kin_offdiag_local = energy_kin_offdiag_local + occ_weight(io) * &
                        sum(real(conjg(coef_frag_all(1:n, io)) * (tmp_probe(1:n, io) - dense_probe_out(1:n, io))))
                    end do
                    deallocate(diag_block_ids)
                    deallocate(dense_probe_out)
                  end if
                end if
                if (itt == 1 .or. itt == 40) then
                  allocate(dense_probe_mat(n, n), dense_probe_out(n, nocc))
                  dense_probe_mat(:, :) = (0.0d0, 0.0d0)
                  call copy_matrix_blocks_to_complex_dense(dg_frag, dg_frag%H_mat_kinetic_blocks, ispin, dense_probe_mat)
                  dense_probe_out(:, :) = matmul(dense_probe_mat(:, :), coef_frag_all(1:n, 1:nocc))
                  kinetic_apply_diff_local = kinetic_apply_diff_local + &
                    sum(abs(tmp_probe(1:n, 1:nocc) - dense_probe_out(1:n, 1:nocc)))
                  deallocate(dense_probe_mat, dense_probe_out)
                end if
                do io = 1, nocc
                  energy_kin_local = energy_kin_local + occ_weight(io) * &
                    sum(real(conjg(coef_frag_all(1:n, io)) * tmp_probe(1:n, io)))
                end do
              end if
              do io = 1, nocc
                energy_static_local = energy_static_local + occ_weight(io) * &
                  sum(real(conjg(coef_frag_all(1:n, io)) * tmp_mat(1:n, io)))
              end do
            end if
            if (has_nonlocal) then
              if (allocated(dg_frag%H_nl_cache)) then
                if (enable_energy_component_probe) then
                  tmp_probe(:, :) = matmul(dg_frag%H_nl_cache(1:n, 1:n, ispin), coef_frag_all(1:n, 1:nocc))
                  do io = 1, nocc
                    energy_nl_local = energy_nl_local + occ_weight(io) * &
                      sum(real(conjg(coef_frag_all(1:n, io)) * tmp_probe(1:n, io)))
                  end do
                end if
                tmp_mat(:, :) = tmp_mat(:, :) + &
                  matmul(dg_frag%H_nl_cache(1:n, 1:n, ispin), coef_frag_all(1:n, 1:nocc))
              else
                if (enable_energy_component_probe) then
                  tmp_probe(:, :) = (0.0d0, 0.0d0)
                  call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_frag_all(1:n, 1:nocc), &
                    tmp_probe)
                  do io = 1, nocc
                    energy_nl_local = energy_nl_local + occ_weight(io) * &
                      sum(real(conjg(coef_frag_all(1:n, io)) * tmp_probe(1:n, io)))
                  end do
                end if
                call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_frag_all(1:n, 1:nocc), &
                  tmp_mat)
              end if
            end if
            if (enable_energy_component_probe) then
              do io = 1, nocc
                energy_a2_local = energy_a2_local + occ_weight(io) * 0.5d0 * A_squared * &
                  sum(abs(coef_frag_all(1:n, io))**2)
              end do
            end if
            tmp_mat(1:n, 1:nocc) = tmp_mat(1:n, 1:nocc) + 0.5d0 * A_squared * coef_frag_all(1:n, 1:nocc)
          else
            if (.not. allocated(op_mat)) allocate(op_mat(n, n))
            do io = 1, n
              op_mat(io, io) = op_mat(io, io) + 0.5d0 * A_squared
            end do
            call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
                       coef_frag_all(1:n, 1:nocc), n, (0.0d0, 0.0d0), tmp_mat, n)
          end if
          if (.not. allocated(op_mat)) allocate(op_mat(n, n))
          op_mat(:, 1:nocc) = (0.0d0, 0.0d0)
          call apply_momentum_blocks(dg_frag, ispin, Ac_tot, coef_frag_all(1:n, 1:nocc), op_mat(:, 1:nocc))
          if (enable_energy_component_probe) then
            do io = 1, nocc
              energy_ap_local = energy_ap_local + occ_weight(io) * &
                sum(real(conjg(coef_frag_all(1:n, io)) * (minus_i * op_mat(:, io))))
            end do
          end if
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
          energy_tmp = energy_tmp + occ_weight(io) * sum(real(conjg(coef_all(1:n_tot, io)) * tmp_all(1:n_tot, io)))
        end do
      else
        if (.not. allocated(dg_frag%momentum_blocks) .or. use_spatial_A) then
        call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
                   coef_frag_all(1:n, 1:nocc), n, (0.0d0, 0.0d0), tmp_mat, n)
        end if
        energy_tmp = 0.0d0
        do io = 1, nocc
          energy_tmp = energy_tmp + occ_weight(io) * sum(real(conjg(coef_frag_all(1:n, io)) * tmp_mat(1:n, io)))
        end do
      end if
        energy_local = energy_local + energy_tmp
      end do

      if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) then
        do ispin = 1, dg_frag%nspin
          nocc = min(dg_frag%nocc_spin(ispin), min(dg_frag%nstate_tot, n))
          if (nocc <= 0) cycle
          coef_pw_all(1:n_pw, 1:nocc) = dg_frag%coef_pw(1:n_pw, 1:nocc, ispin)
          pw_weight_local = pw_weight_local + sum(abs(coef_pw_all(:, 1:nocc))**2)
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

    if (enable_energy_component_probe .and. allocated(dg_frag%H_mat_kinetic_blocks)) then
      do iblk = 1, size(dg_frag%H_mat_kinetic_blocks)
        do ispin = 1, dg_frag%nspin
          nrow_blk = dg_frag%n_basis(dg_frag%H_mat_kinetic_blocks(iblk)%ifrag_row, ispin)
          ncol_blk = dg_frag%n_basis(dg_frag%H_mat_kinetic_blocks(iblk)%ifrag_col, ispin)
          if (nrow_blk <= 0 .or. ncol_blk <= 0) cycle
          if (dg_frag%H_mat_kinetic_blocks(iblk)%ifrag_row == dg_frag%H_mat_kinetic_blocks(iblk)%ifrag_col) then
            kinetic_diag_abs_local = kinetic_diag_abs_local + &
              sum(abs(dg_frag%H_mat_kinetic_blocks(iblk)%val(1:nrow_blk, 1:ncol_blk, ispin)))
          else
            kinetic_offdiag_abs_local = kinetic_offdiag_abs_local + &
              sum(abs(dg_frag%H_mat_kinetic_blocks(iblk)%val(1:nrow_blk, 1:ncol_blk, ispin)))
          end if
        end do
      end do
    end if

    if (allocated(Ap_mat)) deallocate(Ap_mat)
    if (allocated(A2_mat)) deallocate(A2_mat)
    if (allocated(op_mat)) deallocate(op_mat)
    deallocate(tmp_mat)
    if (allocated(tmp_probe)) deallocate(tmp_probe)
    if (allocated(occ_weight)) deallocate(occ_weight)
    if (allocated(coef_frag_all)) deallocate(coef_frag_all)
    if (allocated(coef_pw_all)) deallocate(coef_pw_all)
    if (allocated(coef_all)) deallocate(coef_all)
    if (allocated(tmp_all)) deallocate(tmp_all)
    ! Cache retained for reuse

  1000 continue
    
    if (n_pw == 0) then
      call compute_realspace_energy_probe(dg_frag, system, mg, stencil, itt, Vh, Vxc, Vpsl, &
                                          energy_kin_local, energy_local, energy_kin_rs_sum, energy_one_rs_sum)
    end if

    if (enable_energy_component_probe .and. n_pw == 0 .and. (itt == 1 .or. itt == 40)) then
      write(*,'(1x,a,i0,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
        "        energy-local probe: rank=", dg_frag%id, " itt=", itt, &
        " local_total=", energy_local, " local_static=", energy_static_local, " local_kin=", energy_kin_local
      flush(6)
    end if

    ! MPI aggregation: sum local contributions from all ranks
    call comm_summation(current_local, dg_frag%current, 3, dg_frag%icomm)
    call comm_summation(coef_occ_norm_local, coef_occ_norm_global, dg_frag%icomm)
    if (do_current_block_trace) then
      call comm_summation(current_block_local, current_block_sum, 3 * nblk_trace, dg_frag%icomm)
    end if
    call comm_summation(energy_local, dg_frag%total_energy, dg_frag%icomm)
    call comm_summation(pw_weight_local, dg_frag%pw_weight_raw, dg_frag%icomm)
    frag_reduce_factor = real(max(1, dg_frag%isize_frag), 8)
    dg_frag%total_energy = dg_frag%total_energy / frag_reduce_factor
    if (enable_energy_component_probe) then
      call comm_summation(energy_static_local, energy_static_sum, dg_frag%icomm)
      call comm_summation(energy_kin_local, energy_kin_sum, dg_frag%icomm)
      call comm_summation(energy_nl_local, energy_nl_sum, dg_frag%icomm)
      call comm_summation(energy_ap_local, energy_ap_sum, dg_frag%icomm)
      call comm_summation(energy_a2_local, energy_a2_sum, dg_frag%icomm)
      call comm_summation(energy_kin_diag_local, energy_kin_diag_sum, dg_frag%icomm)
      call comm_summation(energy_kin_offdiag_local, energy_kin_offdiag_sum, dg_frag%icomm)
      call comm_summation(kinetic_diag_abs_local, kinetic_diag_abs_sum, dg_frag%icomm)
      call comm_summation(kinetic_offdiag_abs_local, kinetic_offdiag_abs_sum, dg_frag%icomm)
      call comm_summation(kinetic_apply_diff_local, kinetic_apply_diff_sum, dg_frag%icomm)
      energy_static_sum = energy_static_sum / frag_reduce_factor
      energy_kin_sum = energy_kin_sum / frag_reduce_factor
      energy_nl_sum = energy_nl_sum / frag_reduce_factor
      energy_ap_sum = energy_ap_sum / frag_reduce_factor
      energy_a2_sum = energy_a2_sum / frag_reduce_factor
      energy_kin_diag_sum = energy_kin_diag_sum / frag_reduce_factor
      energy_kin_offdiag_sum = energy_kin_offdiag_sum / frag_reduce_factor
      energy_static_avg = energy_static_sum
      energy_kin_avg = energy_kin_sum
      energy_nl_avg = energy_nl_sum
      energy_ap_avg = energy_ap_sum
      energy_a2_avg = energy_a2_sum
    end if

    ! Current and PW weight are replicated over all ranks, so these remain world-averaged.
    dg_frag%current(:) = dg_frag%current(:) / real(max(1, dg_frag%isize), 8)
    dg_frag%coef_occ_norm = coef_occ_norm_global / real(max(1, dg_frag%isize), 8)
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
    dg_frag%pw_weight_raw = dg_frag%pw_weight_raw / real(max(1, dg_frag%isize), 8)
    if (enable_energy_component_probe) then
      kinetic_diag_abs_sum = kinetic_diag_abs_sum / real(max(1, dg_frag%isize), 8)
      kinetic_offdiag_abs_sum = kinetic_offdiag_abs_sum / real(max(1, dg_frag%isize), 8)
      if (dg_frag%id == 0 .and. n_pw == 0 .and. (itt == 1 .or. mod(itt, 10) == 0)) then
        write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
          "        energy-components: itt=", itt, " total=", dg_frag%total_energy, &
          " static=", energy_static_sum, " kin=", energy_kin_sum, " nl=", energy_nl_sum, &
          " Ap=", energy_ap_sum, " A2=", energy_a2_sum
        write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
          "        energy-components(avg): itt=", itt, " static=", energy_static_avg, " kin=", energy_kin_avg, &
          " nl=", energy_nl_avg, " Ap=", energy_ap_avg, " A2=", energy_a2_avg
        if (itt == 1 .or. itt == 40) then
          write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6)') &
            "        kinetic-split: itt=", itt, " diag=", energy_kin_diag_sum, " offdiag=", energy_kin_offdiag_sum
        end if
        write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6)') &
          "        kinetic-block-norms: itt=", itt, " diag_abs=", kinetic_diag_abs_sum, &
          " offdiag_abs=", kinetic_offdiag_abs_sum
        if (itt == 1 .or. itt == 40) then
          write(*,'(1x,a,i0,a,1pe14.6)') &
            "        kinetic-apply-diff: itt=", itt, " abs_sum=", kinetic_apply_diff_sum
          write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
            "        energy-global compare: itt=", itt, &
            " static_mat=", energy_static_sum, " static_rs=", energy_one_rs_sum, &
            " kin_mat=", energy_kin_sum, " kin_rs=", energy_kin_rs_sum, &
            " vloc_mat=", energy_static_sum - energy_kin_sum, " vloc_rs=", energy_one_rs_sum - energy_kin_rs_sum
          write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
            "        energy-global compare(avg): itt=", itt, &
            " static_mat=", energy_static_avg, " static_rs=", energy_one_rs_sum, &
            " kin_mat=", energy_kin_avg, " kin_rs=", energy_kin_rs_sum, &
            " vloc_mat=", energy_static_avg - energy_kin_avg, " vloc_rs=", energy_one_rs_sum - energy_kin_rs_sum
        end if
        flush(6)
      end if
    end if
    dg_frag%energy_kinetic = 0.0d0
    dg_frag%energy_nonlocal = 0.0d0
    dg_frag%energy_Ap = 0.0d0
    dg_frag%energy_A2 = 0.0d0
    if (enable_energy_component_probe) then
      dg_frag%energy_kinetic = energy_kin_sum
      dg_frag%energy_nonlocal = energy_nl_sum
      dg_frag%energy_Ap = energy_ap_sum
      dg_frag%energy_A2 = energy_a2_sum
    end if

    if (n_pw == 0 .and. (itt == 1 .or. itt == 40)) then
      call debug_vloc_block_probe(dg_frag, system, mg, stencil, Vh, Vxc, Vpsl, itt)
    end if
    
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
      ne_density = nelec_ref / dble(system%ngrid)
    else
      ne_density = 0.0d0
    end if
    dg_frag%current_dia(:) = -ne_density * rt%Ac_tot(:, itt)
    dg_frag%current_total(:) = dg_frag%current_para(:) + dg_frag%current_dia(:)
    if (.not. dg_frag%elec_num_baseline_ready) then
      dg_frag%elec_num_raw_t0 = dg_frag%elec_num_raw
      dg_frag%elec_num_scaled_t0 = dg_frag%elec_num_scaled
      dg_frag%elec_num_baseline_ready = .true.
    end if
    dg_frag%rho_drift_indicator = dg_frag%elec_num_raw - dg_frag%elec_num_raw_t0
    
    ! Store in rt structure for output
    rt%curr(:, itt) = dg_frag%current(:)
    
  end subroutine calculate_observables

  subroutine debug_vloc_block_probe(dg_frag, system, mg, stencil, Vh, Vxc, Vpsl, itt)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_rgrid), intent(in) :: mg
    type(s_stencil), intent(in) :: stencil
    type(s_scalar), intent(in) :: Vh, Vpsl
    type(s_scalar), intent(in) :: Vxc(system%nspin)
    integer, intent(in) :: itt

    integer :: ifrag, i_local, ispin, iblk, nbf, jo, io, nprobe
    real(8), allocatable :: V_total(:,:,:)
    complex(8), allocatable :: V_phi(:,:,:), T_phi(:,:,:), H_phi(:,:,:)
    complex(8) :: integral_v, integral_t, integral_h
    real(8) :: vdiag_direct(3), vdiag_mat(3), tdiag_direct(3), tdiag_mat(3), hdiag_direct(3), hdiag_mat(3)
    real(8) :: voff12_direct, voff12_mat

    if (.not. dg_frag%is_frag_root) return
    if (.not. allocated(dg_frag%H_mat_blocks) .or. .not. allocated(dg_frag%H_mat_kinetic_blocks)) return
    if (dg_frag%ifrag_end < dg_frag%ifrag_start) return

    ispin = 1
    ifrag = dg_frag%ifrag_start
    i_local = 1
    nbf = min(3, dg_frag%n_basis(ifrag, ispin))
    if (nbf <= 0) return

    iblk = find_matrix_block_local(dg_frag%H_block_map, ifrag, ifrag)
    if (iblk <= 0) return

    allocate(V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(V_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(T_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(H_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    call build_total_potential_grid_local(mg, Vh, Vxc(ispin), Vpsl, V_total)

    vdiag_direct(:) = 0.0d0
    vdiag_mat(:) = 0.0d0
    tdiag_direct(:) = 0.0d0
    tdiag_mat(:) = 0.0d0
    hdiag_direct(:) = 0.0d0
    hdiag_mat(:) = 0.0d0
    voff12_direct = 0.0d0
    voff12_mat = 0.0d0

    do jo = 1, nbf
      call build_hpsi_for_basis_probe(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi)
      call integrate_local_basis_with_field_local(dg_frag, ifrag, i_local, jo, mg, T_phi, system%hvol, integral_t)
      call integrate_local_basis_with_field_local(dg_frag, ifrag, i_local, jo, mg, H_phi, system%hvol, integral_h)
      call build_local_potential_applied_basis_local(dg_frag, ifrag, i_local, jo, mg, V_total, V_phi)
      call integrate_local_basis_with_field_local(dg_frag, ifrag, i_local, jo, mg, V_phi, system%hvol, integral_v)
      tdiag_direct(jo) = real(integral_t, kind=8)
      hdiag_direct(jo) = real(integral_h, kind=8)
      vdiag_direct(jo) = real(integral_v, kind=8)
      tdiag_mat(jo) = dg_frag%H_mat_kinetic_blocks(iblk)%val(jo, jo, ispin)
      hdiag_mat(jo) = dg_frag%H_mat_blocks(iblk)%val(jo, jo, ispin)
      vdiag_mat(jo) = dg_frag%H_mat_blocks(iblk)%val(jo, jo, ispin) - dg_frag%H_mat_kinetic_blocks(iblk)%val(jo, jo, ispin)
    end do

    if (nbf >= 2) then
      call build_local_potential_applied_basis_local(dg_frag, ifrag, i_local, 2, mg, V_total, V_phi)
      call integrate_local_basis_with_field_local(dg_frag, ifrag, i_local, 1, mg, V_phi, system%hvol, integral_v)
      voff12_direct = real(integral_v, kind=8)
      voff12_mat = dg_frag%H_mat_blocks(iblk)%val(1, 2, ispin) - dg_frag%H_mat_kinetic_blocks(iblk)%val(1, 2, ispin)
    end if

    write(*,'(1x,a,i0,a,i0,a,3(1pe14.6,1x),a,3(1pe14.6,1x))') &
      "        static-diag probe: rank=", dg_frag%id, " itt=", itt, " h_mat=", &
      hdiag_mat(1), hdiag_mat(2), hdiag_mat(3), " h_rs=", hdiag_direct(1), hdiag_direct(2), hdiag_direct(3)
    write(*,'(1x,a,i0,a,i0,a,3(1pe14.6,1x),a,3(1pe14.6,1x))') &
      "        static-diag probe: rank=", dg_frag%id, " itt=", itt, " t_mat=", &
      tdiag_mat(1), tdiag_mat(2), tdiag_mat(3), " t_rs=", tdiag_direct(1), tdiag_direct(2), tdiag_direct(3)
    write(*,'(1x,a,i0,a,i0,a,3(1pe14.6,1x),a,3(1pe14.6,1x))') &
      "        static-diag probe: rank=", dg_frag%id, " itt=", itt, " v_mat=", &
      vdiag_mat(1), vdiag_mat(2), vdiag_mat(3), " diag_rs=", vdiag_direct(1), vdiag_direct(2), vdiag_direct(3)
    if (nbf >= 2) then
      write(*,'(1x,a,i0,a,i0,a,1pe14.6,a,1pe14.6)') &
        "        static-diag probe: rank=", dg_frag%id, " itt=", itt, " v12_mat=", voff12_mat, " v12_rs=", voff12_direct
    end if
    flush(6)

    deallocate(V_total, V_phi, T_phi, H_phi)
  end subroutine debug_vloc_block_probe

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

    do iz = grid%is(3), grid%ie(3)
      do iy = grid%is(2), grid%ie(2)
        do ix = grid%is(1), grid%ie(1)
          V_total(ix, iy, iz) = Vpsl%f(ix, iy, iz) + Vh%f(ix, iy, iz) + Vxc_spin%f(ix, iy, iz)
        end do
      end do
    end do
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

  subroutine compute_realspace_energy_probe(dg_frag, system, mg, stencil, itt, Vh, Vxc, Vpsl, energy_kin_mat, energy_one_mat, kin_sum_out, one_sum_out)
    use structures
    use communication, only: comm_summation, comm_is_root
    use parallelization, only: nproc_id_global
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    integer,                intent(in)    :: itt
    type(s_scalar),         intent(in)    :: Vh, Vpsl
    type(s_scalar),         intent(in)    :: Vxc(system%nspin)
    real(8),                intent(in)    :: energy_kin_mat, energy_one_mat
    real(8),                intent(out)   :: kin_sum_out, one_sum_out

    integer :: ispin, io, ifrag, i_local, jo, nbf, ig_j
    integer :: loc_s(3), loc_e(3), gx, gy, gz, ixg, iyg, izg
    logical :: has_overlap
    complex(8), allocatable :: psi(:,:,:), tpsi(:,:,:), hpsi(:,:,:)
    real(8), allocatable :: V_total(:,:,:)
    complex(8), allocatable :: T_phi(:,:,:), H_phi(:,:,:)
    real(8), allocatable :: occ_weight_spin(:)
    complex(8) :: coeff, ztmp
    real(8) :: kin_local, one_local, kin_sum, one_sum
    real(8) :: frag_reduce_factor
    integer :: nocc_spin

    kin_sum_out = 0.0d0
    one_sum_out = 0.0d0
    if (dg_frag%use_plane_wave_basis) return
    if (.not. dg_frag%has_real_space_basis) return
    if (dg_frag%n_mat_max <= 0) return

    kin_local = 0.0d0
    one_local = 0.0d0

    do ispin = 1, dg_frag%nspin
      if (dg_frag%nocc_spin(ispin) <= 0) cycle
      nocc_spin = min(dg_frag%nocc_spin(ispin), min(dg_frag%nstate_tot, dg_frag%n_mat_max))
      if (nocc_spin <= 0) cycle
      allocate(occ_weight_spin(nocc_spin))
      call get_dg_spin_occ_info(dg_frag, system, ispin, occ_weight_spin, nocc_spin)
      allocate(V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
      call build_total_potential_grid_local(mg, Vh, Vxc(ispin), Vpsl, V_total)
      allocate(psi(1:dg_frag%lgnum_total(1), 1:dg_frag%lgnum_total(2), 1:dg_frag%lgnum_total(3)))
      allocate(tpsi(1:dg_frag%lgnum_total(1), 1:dg_frag%lgnum_total(2), 1:dg_frag%lgnum_total(3)))
      allocate(hpsi(1:dg_frag%lgnum_total(1), 1:dg_frag%lgnum_total(2), 1:dg_frag%lgnum_total(3)))

      do io = 1, nocc_spin
        psi(:, :, :) = (0.0d0, 0.0d0)
        tpsi(:, :, :) = (0.0d0, 0.0d0)
        hpsi(:, :, :) = (0.0d0, 0.0d0)
        i_local = 0
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = i_local + 1
          nbf = min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag)
          if (nbf <= 0) cycle
          call get_fragment_owned_range(dg_frag, ifrag, mg, loc_s, loc_e, has_overlap)
          if (.not. has_overlap) cycle
          allocate(T_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
          allocate(H_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
          do jo = 1, nbf
            ig_j = dg_frag%index_basis(jo, ifrag, ispin)
            if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle
            coeff = dg_frag%coef(ig_j, io, ispin)
            if (abs(coeff) == 0.0d0) cycle
            call build_hpsi_for_basis_probe(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi)
            do gz = loc_s(3), loc_e(3)
              do gy = loc_s(2), loc_e(2)
                do gx = loc_s(1), loc_e(1)
                  ixg = dg_frag%ixyz_frag(1, ifrag) + gx - 1
                  iyg = dg_frag%ixyz_frag(2, ifrag) + gy - 1
                  izg = dg_frag%ixyz_frag(3, ifrag) + gz - 1
                  call get_phi_value_at_global_probe(dg_frag, ifrag, i_local, jo, &
                    ixg, iyg, izg, ztmp)
                  if (ztmp == (0.0d0, 0.0d0)) cycle
                  psi(ixg, iyg, izg) = psi(ixg, iyg, izg) + coeff * ztmp
                  tpsi(ixg, iyg, izg) = tpsi(ixg, iyg, izg) + coeff * T_phi(ixg, iyg, izg)
                  hpsi(ixg, iyg, izg) = hpsi(ixg, iyg, izg) + coeff * H_phi(ixg, iyg, izg)
                end do
              end do
            end do
          end do
          deallocate(T_phi, H_phi)
        end do
        ztmp = sum(conjg(psi) * tpsi)
        kin_local = kin_local + occ_weight_spin(io) * real(ztmp, kind=8) * system%hvol
        ztmp = sum(conjg(psi) * hpsi)
        one_local = one_local + occ_weight_spin(io) * real(ztmp, kind=8) * system%hvol
      end do

      deallocate(psi, tpsi, hpsi, V_total, occ_weight_spin)
    end do

    call comm_summation(kin_local, kin_sum, dg_frag%icomm)
    call comm_summation(one_local, one_sum, dg_frag%icomm)
    kin_sum_out = kin_sum
    one_sum_out = one_sum
  end subroutine compute_realspace_energy_probe
