  subroutine calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, dcoef_dt, dcoef_dt_pw)
    use structures
    use communication, only: comm_summation, comm_get_max, comm_maxloc_type
    use salmon_global, only: theory, yn_spinorbit
    use rt_dg_fragment_ops, only: apply_momentum_blocks, apply_matrix_blocks_batch, &
                                  apply_complex_matrix_blocks_batch, &
                                  apply_nonlocal_pp_projector_batch, &
                                  apply_nonlocal_pp_projector_batch_so, apply_mixed_hamiltonian, &
                                  solve_overlap_operator_batch, solve_overlap_operator_batch_local, mixed_fp_coupling_active, &
                                  copy_matrix_blocks_metric_to_complex_dense, copy_momentum_blocks_to_complex_dense, gather_full_coef_view, &
                                  apply_overlap_operator
    use misc_routines, only: get_wtime
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag  ! Changed to inout for cache updates
    type(s_dft_system),     intent(in) :: system
    type(s_rgrid),          intent(in) :: mg
    type(s_stencil),        intent(in) :: stencil
    type(s_pp_grid),        intent(in) :: ppg
    real(8),                intent(in) :: Ac_tot(3)  ! vector potential
    integer,                intent(in) :: itt
    complex(8),             intent(out) :: dcoef_dt(:,:,:)
    complex(8), optional,   intent(out) :: dcoef_dt_pw(:,:,:)
    
    integer :: io, jo, ispin, idir, iblk, nrow_blk, ncol_blk
    integer :: n, n_frag, n_pw, n_tot
    real(8) :: A_squared
    complex(8), parameter :: zi = (0.0d0, 1.0d0)  ! imaginary unit
    logical :: has_nonlocal, has_so_nonlocal, has_spin_mix, use_hmat_complex
    logical :: need_h0_dense, need_m_dense
    complex(8), allocatable, save :: H0c(:,:), M(:,:), dcoef_dt_h0(:,:), dcoef_dt_m(:,:), dcoef_dt_h0_core(:,:), &
                     dcoef_dt_nl(:,:), dcoef_dt_spinmix(:,:), rhs_pre_solve(:,:)
    complex(8), allocatable, save :: coef_all(:,:), rhs_all(:,:)
    complex(8), allocatable :: h_core_probe(:,:)
    complex(8), allocatable :: coef_frag_all(:,:), coef_pw_all(:,:)
    complex(8), allocatable :: coef_frag_other(:,:), coef_pw_other(:,:)
    complex(8), allocatable, save :: rhs_in(:,:), rhs_eig(:,:), Uc(:,:), Uc_keep(:,:)
    complex(8), allocatable, save :: coef_mix_all(:,:), rhs_mix(:,:), raw_rhs(:,:), op_mix(:,:), tmp_mix(:,:)
    complex(8), allocatable, save :: op_mix_sum(:,:), raw_rhs_sum(:,:)
    complex(8), allocatable, save :: mfp_vec(:), coef_pw_scaled(:,:), tmp_fp_fw(:,:), tmp_fp_bw(:,:), mfp_fp(:,:)
    complex(8) :: mfp
    real(8) :: huge_val
    real(8) :: t_coef_gather0, t_coef_gather1
    real(8) :: t_deriv0, t_deriv1, dt_deriv, dt_gather_local
    logical :: found_nan
    integer :: nan_io, nan_jo
    real(8) :: max_abs_h0, max_abs_m, op_mix_norm_h0, op_mix_norm_m
    integer :: n_s, n_basis, ispin_other
    integer :: state_s, state_e, state0, nstate_blk
    integer :: nocc_spin, io_occ
    complex(8), allocatable, save :: vec_occ(:), svec_occ(:), sdc_occ(:), vec_src(:), ssrc_vec(:)
    real(8), allocatable, save :: occ_metric_norm(:)
    real(8) :: occ_factor, csc, csc_local, c_s_dc_re, c_s_dc_re_local, alpha_real
    real(8) :: dndt_before, dndt_after, dndt_before_local, dndt_after_local
    logical :: use_spatial_A
    logical :: disable_mfp, disable_local_overlap_solve, bypass_overlap_solve, disable_block_h_apply, disable_nonlocal_pp, use_mixed_basis
    logical :: use_direct_pfp, use_fft_pfp
    logical :: use_prop_overlap, enforce_norm_tangent, enable_norm_deriv_check, enable_excitation_source_trace
    logical :: need_source_decomp, need_spinmix_work
    logical :: enable_overlap_path_trace, force_overlap_solve
    logical :: mfp_coupling_on, has_s_prop_blocks, has_s_blocks, has_s_prop_blocks_c, has_s_blocks_c
    logical :: has_s_prop_c, has_s_c
    character(len=32) :: overlap_gate_reason
    character(len=32) :: env_mfp
    character(len=32) :: env_mfp_mode
    integer :: env_len, env_stat
    real(8), allocatable, save :: Ap_mat(:,:), A2_mat(:,:)
    complex(8), allocatable, save :: ap_block_dense(:,:), ap_block_ref(:,:)
    integer, parameter :: state_block_size = 64
    integer, save :: derivative_step_trace_count = 0
    integer :: derivative_step_trace_id
    logical :: trace_deriv_step
    real(8) :: t_trace0, t_trace1
    real(8) :: t_h0_part0, t_h0_part1
    logical :: enable_deriv_trace
    logical :: enable_hermit_check
    logical :: enable_op_mix_trace
    logical :: enable_m_block_audit
    logical, save :: ap_block_check_initialized = .false.
    logical, save :: enable_ap_block_check = .false.
    logical, save :: derivative_env_initialized = .false.
    logical, save :: cfg_enable_deriv_trace = .false.
    logical, save :: cfg_enable_hermit_check = .false.
    logical, save :: cfg_enable_op_mix_trace = .false.
    logical, save :: cfg_enable_m_block_audit = .true.
    logical, save :: cfg_enable_overlap_path_trace = .true.
    logical, save :: cfg_force_overlap_solve = .false.
    logical, save :: cfg_disable_mfp = .false.
    logical, save :: cfg_disable_local_overlap_solve = .true.
    logical, save :: cfg_bypass_overlap_solve = .false.
    logical, save :: cfg_disable_block_h_apply = .false.
    logical, save :: cfg_disable_nonlocal_pp = .false.
    logical, save :: cfg_use_direct_pfp = .false.
    logical, save :: cfg_use_fft_pfp = .false.
    logical, save :: cfg_use_prop_overlap = .true.
    logical, save :: cfg_enforce_norm_tangent = .false.
    logical, save :: cfg_enable_norm_deriv_check = .false.
    logical, save :: cfg_enable_excitation_source_trace = .false.
    character(len=32) :: env_deriv_trace
    character(len=32) :: env_hermit_check
    character(len=32) :: env_op_mix_trace
    character(len=32) :: env_m_block_audit
    character(len=32), save :: env_ap_block_check
    integer :: env_trace_len, env_trace_stat
    real(8) :: ap_block_diff_max
    real(8) :: m_raw_ff, m_raw_fp, m_raw_pf, m_raw_pp, m_raw_all
    real(8) :: m_rhs_frag_norm, m_rhs_pw_norm, m_rhs_all_norm, m_rhs_max
    real(8) :: rhs_pre_norm, rhs_post_norm, rhs_pre_max, rhs_post_max
    real(8) :: rhs_norm_local_arr(4), rhs_norm_max_arr(4)
    real(8) :: h_split_local(4), h_split_global(4)
    type(comm_maxloc_type) :: rhs_pre_norm_loc, rhs_pre_norm_maxloc
    type(comm_maxloc_type) :: rhs_post_norm_loc, rhs_post_norm_maxloc
    real(8) :: excite_local(8), excite_global(8), excite_per_electron, eps_occ
    real(8) :: h_raw_ff, h_raw_fp, h_raw_pf, h_raw_pp, h_raw_all
    real(8) :: h_max_ff, h_max_fp, h_max_pf, h_max_pp, h_max_all
    real(8) :: m_mix_oo, m_mix_ov, m_mix_vo, m_mix_vv, m_mix_all
    integer :: nocc_basis, nvirt_basis
    integer :: jo_occ, source_idx
    complex(8) :: occ_overlap_local, occ_overlap
    real(8) :: source_norm(6), source_occ_proj(6), source_leak(6), source_norm_global
    
    ! Time derivative in velocity gauge:
    !   d/dt c_i = -i * (H_0_ij + A^2(t)/2 * delta_ij) * c_j - A(t)·<i|∇|j> * c_j
    ! H(t) = H_0 - i*A(t)·∇ + A^2(t)/2
    ! The A^2/2 term is the diamagnetic contribution
    ! Complex coefficients are ESSENTIAL for:
    ! - Phase evolution: exp(-iE_n*t)
    ! - Superposition states
    ! - Oscillatory responses (optical, currents)

    dcoef_dt = (0.0d0, 0.0d0)
    if (present(dcoef_dt_pw)) dcoef_dt_pw = (0.0d0, 0.0d0)
    t_deriv0 = get_wtime()
    dt_gather_local = 0.0d0
    derivative_step_trace_id = 0
    trace_deriv_step = .false.
    huge_val = huge(0.0d0) / 2.0d0

    if (.not. derivative_env_initialized) then
      env_deriv_trace = ''
      call get_environment_variable('SALMON_DG_DERIV_TRACE', env_deriv_trace, length=env_trace_len, status=env_trace_stat)
      if (env_trace_stat == 0 .and. env_trace_len > 0) then
        if (env_deriv_trace(1:1) == '1' .or. env_deriv_trace(1:1) == 'y' .or. env_deriv_trace(1:1) == 'Y' .or. &
            env_deriv_trace(1:1) == 't' .or. env_deriv_trace(1:1) == 'T') cfg_enable_deriv_trace = .true.
      end if

      env_hermit_check = ''
      call get_environment_variable('SALMON_DG_HERMIT_CHECK', env_hermit_check, length=env_trace_len, status=env_trace_stat)
      if (env_trace_stat == 0 .and. env_trace_len > 0) then
        if (env_hermit_check(1:1) == '1' .or. env_hermit_check(1:1) == 'y' .or. env_hermit_check(1:1) == 'Y' .or. &
            env_hermit_check(1:1) == 't' .or. env_hermit_check(1:1) == 'T') cfg_enable_hermit_check = .true.
      end if

      env_op_mix_trace = ''
      call get_environment_variable('SALMON_DG_OP_MIX_TRACE', env_op_mix_trace, length=env_trace_len, status=env_trace_stat)
      if (env_trace_stat == 0 .and. env_trace_len > 0) then
        if (env_op_mix_trace(1:1) == '1' .or. env_op_mix_trace(1:1) == 'y' .or. env_op_mix_trace(1:1) == 'Y' .or. &
            env_op_mix_trace(1:1) == 't' .or. env_op_mix_trace(1:1) == 'T') cfg_enable_op_mix_trace = .true.
      end if

      env_m_block_audit = ''
      call get_environment_variable('SALMON_DG_M_BLOCK_AUDIT', env_m_block_audit, length=env_trace_len, status=env_trace_stat)
      if (env_trace_stat == 0 .and. env_trace_len > 0) then
        if (env_m_block_audit(1:1) == '0' .or. env_m_block_audit(1:1) == 'n' .or. env_m_block_audit(1:1) == 'N' .or. &
            env_m_block_audit(1:1) == 'f' .or. env_m_block_audit(1:1) == 'F') cfg_enable_m_block_audit = .false.
        if (env_m_block_audit(1:1) == '1' .or. env_m_block_audit(1:1) == 'y' .or. env_m_block_audit(1:1) == 'Y' .or. &
            env_m_block_audit(1:1) == 't' .or. env_m_block_audit(1:1) == 'T') cfg_enable_m_block_audit = .true.
      end if

      env_mfp = ''
      call get_environment_variable('SALMON_DG_OVERLAP_SOLVE_TRACE', env_mfp, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        if (env_mfp(1:1) == '0' .or. env_mfp(1:1) == 'n' .or. env_mfp(1:1) == 'N' .or. &
            env_mfp(1:1) == 'f' .or. env_mfp(1:1) == 'F') cfg_enable_overlap_path_trace = .false.
        if (env_mfp(1:1) == '1' .or. env_mfp(1:1) == 'y' .or. env_mfp(1:1) == 'Y' .or. &
            env_mfp(1:1) == 't' .or. env_mfp(1:1) == 'T') cfg_enable_overlap_path_trace = .true.
      end if

      env_mfp = ''
      call get_environment_variable('SALMON_DG_OVERLAP_PATH_TRACE', env_mfp, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        if (env_mfp(1:1) == '0' .or. env_mfp(1:1) == 'n' .or. env_mfp(1:1) == 'N' .or. &
            env_mfp(1:1) == 'f' .or. env_mfp(1:1) == 'F') cfg_enable_overlap_path_trace = .false.
        if (env_mfp(1:1) == '1' .or. env_mfp(1:1) == 'y' .or. env_mfp(1:1) == 'Y' .or. &
            env_mfp(1:1) == 't' .or. env_mfp(1:1) == 'T') cfg_enable_overlap_path_trace = .true.
      end if

      env_mfp = ''
      call get_environment_variable('SALMON_DG_FORCE_OVERLAP_SOLVE', env_mfp, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        if (env_mfp(1:1) == '1' .or. env_mfp(1:1) == 'y' .or. env_mfp(1:1) == 'Y' .or. &
            env_mfp(1:1) == 't' .or. env_mfp(1:1) == 'T') cfg_force_overlap_solve = .true.
      end if

      env_mfp = ''
      call get_environment_variable('SALMON_DG_DISABLE_MFP', env_mfp, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        if (env_mfp(1:1) == '1' .or. env_mfp(1:1) == 'y' .or. env_mfp(1:1) == 'Y' .or. &
            env_mfp(1:1) == 't' .or. env_mfp(1:1) == 'T') cfg_disable_mfp = .true.
      end if

      env_mfp_mode = ''
      call get_environment_variable('SALMON_DG_MFP_MODE', env_mfp_mode, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        if (env_mfp_mode(1:1) == 'f' .or. env_mfp_mode(1:1) == 'F' .or. env_mfp_mode(1:1) == 'g' .or. env_mfp_mode(1:1) == 'G') then
          cfg_use_direct_pfp = .true.
          cfg_use_fft_pfp = .true.
        else if (env_mfp_mode(1:1) == 'p' .or. env_mfp_mode(1:1) == 'P' .or. &
                 env_mfp_mode(1:1) == 'd' .or. env_mfp_mode(1:1) == 'D' .or. &
                 env_mfp_mode(1:1) == 'n' .or. env_mfp_mode(1:1) == 'N' .or. &
                 env_mfp_mode(1:1) == '1') then
          cfg_use_direct_pfp = .true.
          cfg_use_fft_pfp = .false.
        end if
      end if

      env_mfp = ''
      call get_environment_variable('SALMON_DG_DISABLE_LOCAL_OVERLAP_SOLVE', env_mfp, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        if (env_mfp(1:1) == '0' .or. env_mfp(1:1) == 'n' .or. env_mfp(1:1) == 'N' .or. &
            env_mfp(1:1) == 'f' .or. env_mfp(1:1) == 'F') cfg_disable_local_overlap_solve = .false.
      end if

      env_mfp = ''
      call get_environment_variable('SALMON_DG_BYPASS_OVERLAP_SOLVE', env_mfp, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        if (env_mfp(1:1) == '1' .or. env_mfp(1:1) == 'y' .or. env_mfp(1:1) == 'Y' .or. &
            env_mfp(1:1) == 't' .or. env_mfp(1:1) == 'T') cfg_bypass_overlap_solve = .true.
      end if

      env_mfp = ''
      call get_environment_variable('SALMON_DG_DISABLE_BLOCK_H_APPLY', env_mfp, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        if (env_mfp(1:1) == '1' .or. env_mfp(1:1) == 'y' .or. env_mfp(1:1) == 'Y' .or. &
            env_mfp(1:1) == 't' .or. env_mfp(1:1) == 'T') cfg_disable_block_h_apply = .true.
      end if

      env_mfp = ''
      call get_environment_variable('SALMON_DG_DISABLE_NONLOCAL_PP', env_mfp, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        if (env_mfp(1:1) == '1' .or. env_mfp(1:1) == 'y' .or. env_mfp(1:1) == 'Y' .or. &
            env_mfp(1:1) == 't' .or. env_mfp(1:1) == 'T') cfg_disable_nonlocal_pp = .true.
      end if

      env_mfp = ''
      call get_environment_variable('SALMON_DG_DISABLE_PROP_OVERLAP', env_mfp, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        if (env_mfp(1:1) == '1' .or. env_mfp(1:1) == 'y' .or. env_mfp(1:1) == 'Y' .or. &
            env_mfp(1:1) == 't' .or. env_mfp(1:1) == 'T') cfg_use_prop_overlap = .false.
      end if

      env_mfp = ''
      call get_environment_variable('SALMON_DG_ENFORCE_NORM_TANGENT', env_mfp, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        if (env_mfp(1:1) == '1' .or. env_mfp(1:1) == 'y' .or. env_mfp(1:1) == 'Y' .or. &
            env_mfp(1:1) == 't' .or. env_mfp(1:1) == 'T') cfg_enforce_norm_tangent = .true.
      end if

      env_mfp = ''
      call get_environment_variable('SALMON_DG_NORM_DERIV_CHECK', env_mfp, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        if (env_mfp(1:1) == '1' .or. env_mfp(1:1) == 'y' .or. env_mfp(1:1) == 'Y' .or. &
            env_mfp(1:1) == 't' .or. env_mfp(1:1) == 'T') cfg_enable_norm_deriv_check = .true.
      end if

      env_mfp = ''
      call get_environment_variable('SALMON_DG_EXCITATION_SOURCE_TRACE', env_mfp, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        if (env_mfp(1:1) == '1' .or. env_mfp(1:1) == 'y' .or. env_mfp(1:1) == 'Y' .or. &
            env_mfp(1:1) == 't' .or. env_mfp(1:1) == 'T') cfg_enable_excitation_source_trace = .true.
      end if

      derivative_env_initialized = .true.
    end if

    enable_deriv_trace = cfg_enable_deriv_trace
    enable_hermit_check = cfg_enable_hermit_check
    enable_op_mix_trace = cfg_enable_op_mix_trace
    enable_m_block_audit = cfg_enable_m_block_audit
    enable_overlap_path_trace = cfg_enable_overlap_path_trace
    force_overlap_solve = cfg_force_overlap_solve
    disable_mfp = cfg_disable_mfp
    disable_local_overlap_solve = cfg_disable_local_overlap_solve
    bypass_overlap_solve = cfg_bypass_overlap_solve
    disable_block_h_apply = cfg_disable_block_h_apply
    disable_nonlocal_pp = cfg_disable_nonlocal_pp
    use_direct_pfp = cfg_use_direct_pfp
    use_fft_pfp = cfg_use_fft_pfp
    use_prop_overlap = cfg_use_prop_overlap
    enforce_norm_tangent = cfg_enforce_norm_tangent
    enable_norm_deriv_check = cfg_enable_norm_deriv_check
    enable_excitation_source_trace = cfg_enable_excitation_source_trace

    if (enable_deriv_trace .and. itt <= 2) then
      write(*,'(1x,a,i0,a,i0,a,3(1x,1pe12.4))') "        deriv-trace: rank=", dg_frag%id, " itt=", itt, &
        " Ac_tot=", Ac_tot(1), Ac_tot(2), Ac_tot(3)
      flush(6)
    end if
    ! Calculate A^2 (diamagnetic term)
    A_squared = Ac_tot(1)**2 + Ac_tot(2)**2 + Ac_tot(3)**2
    if (A_squared /= A_squared) then
      write(*,'(a,i0,a,i0,a,3es12.4)') "[NaN] A_squared: rank=", dg_frag%id, " itt=", itt, " Ac_tot=", &
        Ac_tot(1), Ac_tot(2), Ac_tot(3)
      stop "NaN in A_squared"
    end if

    has_nonlocal = (ppg%Nlma > 0 .and. allocated(ppg%uV))
    ! SOI projector path supports both real and complex fragment basis via cached complex views.
    has_so_nonlocal = (allocated(ppg%uv_so) .and. yn_spinorbit == 'y' .and. dg_frag%nspin == 2)
    ! Explicit spin-mixing block is not stored in current s_dg_fragment_rt.
    has_spin_mix = .false.
    
    n = size(dg_frag%coef, 1)
    if (n <= 0) return
    n_frag = size(dg_frag%coef, 1)
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n_frag + n_pw

    use_spatial_A = (trim(theory) == 'single_scale_maxwell_tddft' .and. allocated(system%Ac_micro%v) .and. dg_frag%has_real_space_basis)
    if (disable_nonlocal_pp) then
      has_nonlocal = .false.
      has_so_nonlocal = .false.
    end if
    if (itt == 1 .and. dg_frag%id == 0 .and. derivative_step_trace_count < 8) then
      derivative_step_trace_count = derivative_step_trace_count + 1
      derivative_step_trace_id = derivative_step_trace_count
      trace_deriv_step = .true.
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') '[DG-DERIV] enter itt=', itt, &
        ' call=', derivative_step_trace_id, ' n_frag=', n_frag, ' n_pw=', n_pw, &
        ' nstate_tot=', dg_frag%nstate_tot, ' nspin=', dg_frag%nspin
      write(*,'(1x,a,l1,a,l1,a,l1,a,l1,a,l1)') '[DG-DERIV] flags has_nonlocal=', has_nonlocal, &
        ' has_so_nonlocal=', has_so_nonlocal, ' use_spatial_A=', use_spatial_A, &
        ' mixed_ready=', dg_frag%mixed_basis_ready, ' orbital_mode=', dg_frag%parallel_mode_orbital
      flush(6)
    end if
    if (use_spatial_A) then
      if (.not. allocated(Ap_mat)) then
        allocate(Ap_mat(n, n), A2_mat(n, n))
      else if (size(Ap_mat, 1) /= n .or. size(Ap_mat, 2) /= n) then
        deallocate(Ap_mat, A2_mat)
        allocate(Ap_mat(n, n), A2_mat(n, n))
      end if
    end if

    need_source_decomp = enable_excitation_source_trace
    need_spinmix_work = need_source_decomp .or. has_spin_mix

    if (.not. allocated(dcoef_dt_h0)) then
      allocate(dcoef_dt_h0(n_tot, dg_frag%nstate_tot), dcoef_dt_m(n_tot, dg_frag%nstate_tot))
    else if (size(dcoef_dt_h0, 1) /= n_tot .or. size(dcoef_dt_h0, 2) /= dg_frag%nstate_tot) then
      deallocate(dcoef_dt_h0, dcoef_dt_m)
      allocate(dcoef_dt_h0(n_tot, dg_frag%nstate_tot), dcoef_dt_m(n_tot, dg_frag%nstate_tot))
    end if
    if (need_source_decomp) then
      if (.not. allocated(dcoef_dt_h0_core)) then
        allocate(dcoef_dt_h0_core(n_tot, dg_frag%nstate_tot), dcoef_dt_nl(n_tot, dg_frag%nstate_tot), &
                 rhs_pre_solve(n_tot, dg_frag%nstate_tot))
      else if (size(dcoef_dt_h0_core, 1) /= n_tot .or. size(dcoef_dt_h0_core, 2) /= dg_frag%nstate_tot) then
        deallocate(dcoef_dt_h0_core, dcoef_dt_nl, rhs_pre_solve)
        allocate(dcoef_dt_h0_core(n_tot, dg_frag%nstate_tot), dcoef_dt_nl(n_tot, dg_frag%nstate_tot), &
                 rhs_pre_solve(n_tot, dg_frag%nstate_tot))
      end if
    else
      if (allocated(dcoef_dt_h0_core)) deallocate(dcoef_dt_h0_core)
      if (allocated(dcoef_dt_nl)) deallocate(dcoef_dt_nl)
      if (allocated(rhs_pre_solve)) deallocate(rhs_pre_solve)
    end if
    if (need_spinmix_work) then
      if (.not. allocated(dcoef_dt_spinmix)) then
        allocate(dcoef_dt_spinmix(n_tot, dg_frag%nstate_tot))
      else if (size(dcoef_dt_spinmix, 1) /= n_tot .or. size(dcoef_dt_spinmix, 2) /= dg_frag%nstate_tot) then
        deallocate(dcoef_dt_spinmix)
        allocate(dcoef_dt_spinmix(n_tot, dg_frag%nstate_tot))
      end if
    else if (allocated(dcoef_dt_spinmix)) then
      deallocate(dcoef_dt_spinmix)
    end if
    if (.not. allocated(coef_all)) then
      allocate(coef_all(n_tot, dg_frag%nstate_tot), rhs_all(n_tot, dg_frag%nstate_tot))
    else if (size(coef_all, 1) /= n_tot .or. size(coef_all, 2) /= dg_frag%nstate_tot) then
      deallocate(coef_all, rhs_all)
      allocate(coef_all(n_tot, dg_frag%nstate_tot), rhs_all(n_tot, dg_frag%nstate_tot))
    end if

    do ispin = 1, dg_frag%nspin
      if (enable_deriv_trace .and. itt <= 2) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        deriv-trace: rank=", dg_frag%id, " itt=", itt, &
          " ispin=", ispin, " stage=spin-begin n_frag=", n_frag
        flush(6)
      end if
      ! Build H0c = H_0 + V_NL(A) + A^2/2
      use_hmat_complex = allocated(dg_frag%H_mat_c) .and. allocated(dg_frag%phi_frag_c)
      use_mixed_basis = (n_pw > 0 .and. dg_frag%mixed_basis_ready .and. allocated(dg_frag%mixed_transform) .and. &
                         allocated(dg_frag%mixed_basis_dim) .and. allocated(dg_frag%coef_mix))
      n_basis = 0
      if (use_mixed_basis) n_basis = min(dg_frag%mixed_basis_dim(ispin), n_tot, size(dg_frag%coef_mix, 1))
      need_h0_dense = use_spatial_A .or. use_hmat_complex .or. (.not. allocated(dg_frag%H_mat_blocks)) .or. &
              (n_pw == 0 .and. disable_block_h_apply) .or. use_mixed_basis
      need_m_dense = use_spatial_A .or. (.not. allocated(dg_frag%momentum_blocks)) .or. use_mixed_basis
      if (trace_deriv_step) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,l1,a,l1,a,l1,a,l1)') &
          '[DG-DERIV] spin begin call=', derivative_step_trace_id, ' ispin=', ispin, &
          ' n_basis=', n_basis, ' n_tot=', n_tot, ' use_mixed=', use_mixed_basis, &
          ' need_h0_dense=', need_h0_dense, ' need_m_dense=', need_m_dense, &
          ' has_blocks=', allocated(dg_frag%H_mat_blocks)
        flush(6)
      end if
      if (need_h0_dense) then
        write(*,'(1x,a,i0,a,i0,a,i0)') '[FATAL] derivative would require dense H0c(n_tot,n_tot); block/raw route only: rank=', &
             dg_frag%id, ' itt=', itt, ' ispin=', ispin
        stop 1
      end if
      if (need_m_dense) then
        write(*,'(1x,a,i0,a,i0,a,i0)') '[FATAL] derivative would require dense M(n_tot,n_tot); block/raw route only: rank=', &
             dg_frag%id, ' itt=', itt, ' ispin=', ispin
        stop 1
      end if
      if (allocated(H0c)) H0c(:, :) = (0.0d0, 0.0d0)
      if (allocated(M)) M(:, :) = (0.0d0, 0.0d0)

      if (need_h0_dense .and. use_hmat_complex) then
        H0c(1:n_frag, 1:n_frag) = dg_frag%H_mat_c(1:n_frag, 1:n_frag, ispin)
      else if (need_h0_dense .and. allocated(dg_frag%H_mat_blocks)) then
        call copy_matrix_blocks_metric_to_complex_dense(dg_frag, dg_frag%H_mat_blocks, ispin, n_frag, H0c(1:n_frag, 1:n_frag))
      else if (need_h0_dense .and. .not. allocated(dg_frag%H_mat_blocks)) then
        write(*,'(1x,a,i0,a,i0,a,i0)') '[FATAL] derivative requires H_mat_blocks (dense H_mat route removed): rank=', &
             dg_frag%id, ' itt=', itt, ' ispin=', ispin
        stop 1
      end if

      ! Fill PW-PW diagonal and FP off-diagonal blocks of H0c.
      ! These were missing, causing the mixed-basis propagation (T^H*H0c*T)
      ! to use incorrect energies for PW states and no static FP coupling.
      if (need_h0_dense .and. n_pw > 0) then
        if (allocated(dg_frag%H_mat_pw)) then
          H0c(n_frag+1:n_frag+n_pw, n_frag+1:n_frag+n_pw) = dg_frag%H_mat_pw(1:n_pw, 1:n_pw, ispin)
        else if (allocated(dg_frag%H_mat_pw_diag)) then
          do io = 1, n_pw
            H0c(n_frag+io, n_frag+io) = dg_frag%H_mat_pw_diag(io, ispin)
          end do
        else
          do io = 1, n_pw
            H0c(n_frag+io, n_frag+io) = cmplx(0.5d0 * sum(dg_frag%k_pw(1:3,io)**2), 0.0d0, kind=8)
          end do
        end if
        if (.not. allocated(dg_frag%H_mat_frag_pw)) then
          write(*,'(1x,a,i0,a,i0,a,i0)') '[FATAL] derivative requires H_mat_frag_pw when n_pw>0 (FP/PF off-diagonal is mandatory): rank=', &
               dg_frag%id, ' itt=', itt, ' ispin=', ispin
          stop 1
        end if
        H0c(1:n_frag, n_frag+1:n_frag+n_pw) = dg_frag%H_mat_frag_pw(1:n_frag, 1:n_pw, ispin)
        H0c(n_frag+1:n_frag+n_pw, 1:n_frag) = conjg(transpose(dg_frag%H_mat_frag_pw(1:n_frag, 1:n_pw, ispin)))
      end if

      if (use_spatial_A) then
        call build_spatial_A_coupling_matrices(dg_frag, system, mg, stencil, ispin, Ap_mat, A2_mat)
        H0c(:, :) = H0c(:, :) + cmplx(A2_mat(:, :), 0.0d0, kind=8)
        M(:, :) = cmplx(Ap_mat(:, :), 0.0d0, kind=8)
      else
        if (need_h0_dense) then
          do io = 1, n_tot
            H0c(io, io) = H0c(io, io) + 0.5d0 * A_squared
          end do
        end if

        ! Build M = A·<∇>
        if (need_m_dense .and. allocated(dg_frag%momentum_blocks)) then
          call copy_momentum_blocks_to_complex_dense(dg_frag, ispin, Ac_tot, M(1:n_frag, 1:n_frag))
        else if (.not. allocated(dg_frag%momentum_blocks)) then
          write(*,'(1x,a,i0,a,i0,a,i0)') '[FATAL] derivative requires momentum_blocks (dense momentum route removed): rank=', &
               dg_frag%id, ' itt=', itt, ' ispin=', ispin
          stop 1
        end if
        if (need_m_dense .and. n_pw > 0) then
          if (.not. allocated(mfp_vec)) then
            allocate(mfp_vec(n_pw))
          else if (size(mfp_vec, 1) /= n_pw) then
            deallocate(mfp_vec)
            allocate(mfp_vec(n_pw))
          end if
!$omp simd
          do io = 1, n_pw
            mfp_vec(io) = zi * dot_product(Ac_tot(1:3), dg_frag%k_pw(1:3, io))
            M(n_frag+io, n_frag+io) = mfp_vec(io)
          end do
          if (.not. disable_mfp .and. mixed_fp_coupling_active(dg_frag, ispin)) then
            if (use_direct_pfp .and. ((use_fft_pfp .and. allocated(dg_frag%P_mat_frag_pw_g)) .or. allocated(dg_frag%P_mat_frag_pw))) then
              do jo = 1, n_pw
                if (use_fft_pfp .and. allocated(dg_frag%P_mat_frag_pw_g)) then
                  M(1:n_frag, n_frag+jo) = Ac_tot(1) * dg_frag%P_mat_frag_pw_g(1, 1:n_frag, jo, ispin) + &
                                           Ac_tot(2) * dg_frag%P_mat_frag_pw_g(2, 1:n_frag, jo, ispin) + &
                                           Ac_tot(3) * dg_frag%P_mat_frag_pw_g(3, 1:n_frag, jo, ispin)
                else
                  M(1:n_frag, n_frag+jo) = Ac_tot(1) * dg_frag%P_mat_frag_pw(1, 1:n_frag, jo, ispin) + &
                                           Ac_tot(2) * dg_frag%P_mat_frag_pw(2, 1:n_frag, jo, ispin) + &
                                           Ac_tot(3) * dg_frag%P_mat_frag_pw(3, 1:n_frag, jo, ispin)
                end if
                M(n_frag+jo, 1:n_frag) = -conjg(M(1:n_frag, n_frag+jo))
              end do
            else
!$omp parallel do schedule(static) private(jo)
              do jo = 1, n_pw
                M(1:n_frag, n_frag+jo) = mfp_vec(jo) * dg_frag%S_mat_frag_pw(1:n_frag, jo, ispin)
                M(n_frag+jo, 1:n_frag) = -conjg(M(1:n_frag, n_frag+jo))
              end do
!$omp end parallel do
            end if
          end if
        end if
      end if
      if (need_m_dense) then
        if (any(real(M(:, :)) /= real(M(:, :))) .or. any(aimag(M(:, :)) /= aimag(M(:, :)))) then
            write(*,'(a,i0,a,i0,a,i0)') "[NaN] M: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
            stop "NaN in M"
        end if
      end if
      if (need_m_dense) then
        if (any(abs(M(:, :)) > huge_val)) then
          write(*,'(a,i0,a,i0,a,i0)') "[Inf] M: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
          stop "Inf in M"
        end if
      end if

      if (n_frag > size(dg_frag%coef, 1)) then
        write(*,'(1x,a,i0,a,i0,a,i0)') "[FATAL] derivative n_frag exceeds coef rows: rank=", dg_frag%id, &
          " ispin=", ispin, " n_frag=", n_frag
        stop "DG derivative invalid n_frag/coefficient shape"
      end if
      if (dg_frag%nstate_tot > size(dg_frag%coef, 2)) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] derivative nstate exceeds coef cols: rank=", dg_frag%id, &
          " ispin=", ispin, " nstate_tot=", dg_frag%nstate_tot, " coef_cols=", size(dg_frag%coef, 2)
        stop "DG derivative invalid nstate/coefficient shape"
      end if
      if (trace_deriv_step) then
        t_trace0 = get_wtime()
        write(*,'(1x,a,i0,a,i0,a,i0)') '[DG-DERIV] coef gather start call=', derivative_step_trace_id, &
          ' ispin=', ispin, ' block_size=', state_block_size
        flush(6)
      end if
      t_coef_gather0 = get_wtime()
      if (dg_frag%parallel_mode_orbital) then
        if (trace_deriv_step) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') '[DG-DERIV] coef local-view copy start call=', &
            derivative_step_trace_id, ' ispin=', ispin, ' n_rows=', n_frag, ' n_cols=', dg_frag%nstate_tot
          flush(6)
        end if
        ! Orbital mode stores only this rank's basis-row slice in coef(:,:,:).
        ! Rebuilding it with gather_full_coef_view would require a different
        ! collective sequence in each fragment subgroup, so the derivative path
        ! must consume the local row view directly.
        coef_all(1:n_frag, 1:dg_frag%nstate_tot) = dg_frag%coef(1:n_frag, 1:dg_frag%nstate_tot, ispin)
        if (n_pw > 0) coef_all(n_frag+1:n_tot, 1:dg_frag%nstate_tot) = &
          dg_frag%coef_pw(1:n_pw, 1:dg_frag%nstate_tot, ispin)
        if (trace_deriv_step) then
          write(*,'(1x,a,i0,a,i0)') '[DG-DERIV] coef local-view copy done call=', &
            derivative_step_trace_id, ' ispin=', ispin
          flush(6)
        end if
      else
        coef_all(:, :) = (0.0d0, 0.0d0)
        do state0 = 1, dg_frag%nstate_tot, state_block_size
          nstate_blk = min(state_block_size, dg_frag%nstate_tot - state0 + 1)
          state_s = state0
          state_e = state0 + nstate_blk - 1
          call gather_full_coef_view(dg_frag, ispin, n_frag, nstate_blk, coef_frag_all, coef_pw_all, state_s, state_e)
          coef_all(1:n_frag, state_s:state_e) = coef_frag_all(1:n_frag, 1:nstate_blk)
          if (n_pw > 0) coef_all(n_frag+1:n_tot, state_s:state_e) = coef_pw_all(1:n_pw, 1:nstate_blk)
        end do
      end if
      if (has_so_nonlocal .or. has_spin_mix) then
        ispin_other = 3 - ispin
        if (.not. allocated(coef_frag_other)) then
          allocate(coef_frag_other(max(0, n_frag), max(0, dg_frag%nstate_tot)))
        else if (size(coef_frag_other, 1) /= max(0, n_frag) .or. size(coef_frag_other, 2) /= max(0, dg_frag%nstate_tot)) then
          deallocate(coef_frag_other)
          allocate(coef_frag_other(max(0, n_frag), max(0, dg_frag%nstate_tot)))
        end if
        coef_frag_other(:, :) = (0.0d0, 0.0d0)
        if (n_pw > 0) then
          if (.not. allocated(coef_pw_other)) then
            allocate(coef_pw_other(max(0, n_pw), max(0, dg_frag%nstate_tot)))
          else if (size(coef_pw_other, 1) /= max(0, n_pw) .or. size(coef_pw_other, 2) /= max(0, dg_frag%nstate_tot)) then
            deallocate(coef_pw_other)
            allocate(coef_pw_other(max(0, n_pw), max(0, dg_frag%nstate_tot)))
          end if
          coef_pw_other(:, :) = (0.0d0, 0.0d0)
        end if
        if (dg_frag%parallel_mode_orbital) then
          coef_frag_other(1:n_frag, 1:dg_frag%nstate_tot) = &
            dg_frag%coef(1:n_frag, 1:dg_frag%nstate_tot, ispin_other)
          if (n_pw > 0) coef_pw_other(1:n_pw, 1:dg_frag%nstate_tot) = &
            dg_frag%coef_pw(1:n_pw, 1:dg_frag%nstate_tot, ispin_other)
        else
          do state0 = 1, dg_frag%nstate_tot, state_block_size
            nstate_blk = min(state_block_size, dg_frag%nstate_tot - state0 + 1)
            state_s = state0
            state_e = state0 + nstate_blk - 1
            call gather_full_coef_view(dg_frag, ispin_other, n_frag, nstate_blk, coef_frag_all, coef_pw_all, state_s, state_e)
            coef_frag_other(1:n_frag, state_s:state_e) = coef_frag_all(1:n_frag, 1:nstate_blk)
            if (n_pw > 0) coef_pw_other(1:n_pw, state_s:state_e) = coef_pw_all(1:n_pw, 1:nstate_blk)
          end do
        end if
      end if
      t_coef_gather1 = get_wtime()
      dt_gather_local = dt_gather_local + (t_coef_gather1 - t_coef_gather0)
      if (trace_deriv_step) then
        t_trace1 = get_wtime()
        write(*,'(1x,a,i0,a,i0,a,1pe12.4)') '[DG-DERIV] coef gather done call=', &
          derivative_step_trace_id, ' ispin=', ispin, ' time=', t_trace1 - t_trace0
        call trace_deriv_norm('coef', coef_all)
        flush(6)
      end if
      if (enable_deriv_trace .and. itt <= 2) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        deriv-trace: rank=", dg_frag%id, " itt=", itt, &
          " ispin=", ispin, " stage=after-coef-gather"
        flush(6)
      end if

      if (use_mixed_basis .and. n_basis > 0) then
        if (.not. allocated(coef_mix_all)) then
          allocate(coef_mix_all(n_basis, dg_frag%nstate_tot))
        else if (size(coef_mix_all, 1) /= n_basis .or. size(coef_mix_all, 2) /= dg_frag%nstate_tot) then
          deallocate(coef_mix_all)
          allocate(coef_mix_all(n_basis, dg_frag%nstate_tot))
        end if
        coef_mix_all(:, :) = dg_frag%coef_mix(1:n_basis, 1:dg_frag%nstate_tot, ispin)
      end if

      if (any(real(coef_all(:, :)) /= real(coef_all(:, :))) .or. &
          any(aimag(coef_all(:, :)) /= aimag(coef_all(:, :)))) then
        write(*,'(a,i0,a,i0,a,i0)') "[NaN] coef: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
        stop "NaN in coef before zgemm"
      end if
      if (any(abs(coef_all(:, :)) > huge_val)) then
        write(*,'(a,i0,a,i0,a,i0)') "[Inf] coef: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
        stop "Inf in coef before zgemm"
      end if


      if (need_h0_dense) then
        if (any(abs(H0c(:, :)) > huge_val)) then
          write(*,'(a,i0,a,i0,a,i0)') "[Inf] H0c: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
          stop "Inf in H0c"
        end if
      end if
      if (enable_m_block_audit .and. need_h0_dense .and. dg_frag%id == 0 .and. itt <= 5) then
        h_raw_ff = 0.0d0
        h_raw_fp = 0.0d0
        h_raw_pf = 0.0d0
        h_raw_pp = 0.0d0
        h_max_ff = 0.0d0
        h_max_fp = 0.0d0
        h_max_pf = 0.0d0
        h_max_pp = 0.0d0
        if (n_frag > 0) then
          h_raw_ff = sqrt(sum(abs(H0c(1:n_frag, 1:n_frag))**2))
          h_max_ff = maxval(abs(H0c(1:n_frag, 1:n_frag)))
        end if
        if (n_frag > 0 .and. n_pw > 0) then
          h_raw_fp = sqrt(sum(abs(H0c(1:n_frag, n_frag+1:n_tot))**2))
          h_raw_pf = sqrt(sum(abs(H0c(n_frag+1:n_tot, 1:n_frag))**2))
          h_max_fp = maxval(abs(H0c(1:n_frag, n_frag+1:n_tot)))
          h_max_pf = maxval(abs(H0c(n_frag+1:n_tot, 1:n_frag)))
        end if
        if (n_pw > 0) then
          h_raw_pp = sqrt(sum(abs(H0c(n_frag+1:n_tot, n_frag+1:n_tot))**2))
          h_max_pp = maxval(abs(H0c(n_frag+1:n_tot, n_frag+1:n_tot)))
        end if
        h_raw_all = sqrt(sum(abs(H0c(1:n_tot, 1:n_tot))**2))
        h_max_all = maxval(abs(H0c(1:n_tot, 1:n_tot)))
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') '[H0-BLOCK-AUDIT] itt=', itt, ' ispin=', ispin, &
             ' n_frag=', n_frag, ' n_pw=', n_pw, ' n_basis=', n_basis
        write(*,'(1x,a,5(1x,1pe14.6))') '[H0-BLOCK-AUDIT] frob ff fp pf pp all=', &
             h_raw_ff, h_raw_fp, h_raw_pf, h_raw_pp, h_raw_all
        write(*,'(1x,a,5(1x,1pe14.6))') '[H0-BLOCK-AUDIT] max  ff fp pf pp all=', &
             h_max_ff, h_max_fp, h_max_pf, h_max_pp, h_max_all
        flush(6)
      end if
      if (enable_hermit_check .and. need_h0_dense .and. itt <= 120 .and. dg_frag%id == 0) then
        max_abs_h0 = maxval(abs(H0c(1:n_tot, 1:n_tot) - transpose(conjg(H0c(1:n_tot, 1:n_tot)))))
        write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4)') '        hermit-check: rank=', dg_frag%id, ' itt=', itt, &
          ' ispin=', ispin, ' max|H-H^H|=', max_abs_h0
        flush(6)
      end if


      ! dcoef_dt = -i * H0c * coef - M * coef
      dcoef_dt_h0(:, :) = (0.0d0, 0.0d0)
      if (need_source_decomp) then
        dcoef_dt_h0_core(:, :) = (0.0d0, 0.0d0)
        dcoef_dt_nl(:, :) = (0.0d0, 0.0d0)
      end if
      if (need_spinmix_work) dcoef_dt_spinmix(:, :) = (0.0d0, 0.0d0)
      if (trace_deriv_step) then
        t_trace0 = get_wtime()
        write(*,'(1x,a,i0,a,i0,a,l1,a,l1,a,l1,a,l1)') '[DG-DERIV] H0 apply start call=', &
          derivative_step_trace_id, ' ispin=', ispin, ' use_mixed=', use_mixed_basis, &
          ' has_nl=', has_nonlocal, ' has_so_nl=', has_so_nonlocal, &
          ' block_h=', allocated(dg_frag%H_mat_blocks)
        flush(6)
      end if
      if (use_mixed_basis .and. n_basis > 0) then

        if (.not. allocated(rhs_mix)) then
          allocate(rhs_mix(n_basis, dg_frag%nstate_tot), raw_rhs(n_tot, dg_frag%nstate_tot), op_mix(n_basis, n_basis), &
                   tmp_mix(n_tot, n_basis), op_mix_sum(n_basis, n_basis), raw_rhs_sum(n_tot, dg_frag%nstate_tot))
        else if (size(rhs_mix, 1) /= n_basis .or. size(rhs_mix, 2) /= dg_frag%nstate_tot .or. &
                 size(raw_rhs, 1) /= n_tot .or. size(raw_rhs, 2) /= dg_frag%nstate_tot .or. &
                 size(op_mix, 1) /= n_basis .or. size(op_mix, 2) /= n_basis .or. &
                 size(tmp_mix, 1) /= n_tot .or. size(tmp_mix, 2) /= n_basis .or. &
                 size(op_mix_sum, 1) /= n_basis .or. size(op_mix_sum, 2) /= n_basis .or. &
                 size(raw_rhs_sum, 1) /= n_tot .or. size(raw_rhs_sum, 2) /= dg_frag%nstate_tot) then
          deallocate(rhs_mix, raw_rhs, op_mix, tmp_mix, op_mix_sum, raw_rhs_sum)
          allocate(rhs_mix(n_basis, dg_frag%nstate_tot), raw_rhs(n_tot, dg_frag%nstate_tot), op_mix(n_basis, n_basis), &
                   tmp_mix(n_tot, n_basis), op_mix_sum(n_basis, n_basis), raw_rhs_sum(n_tot, dg_frag%nstate_tot))
        end if
        rhs_mix(:, :) = (0.0d0, 0.0d0)
        raw_rhs(:, :) = (0.0d0, 0.0d0)
        call build_mixed_projected_operator(H0c, op_mix)
        if (enable_op_mix_trace .and. dg_frag%id == 0 .and. (itt <= 200 .or. mod(itt, 50) == 0)) then
          op_mix_norm_h0 = maxval(abs(op_mix(1:n_basis, 1:n_basis)))
          write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4)') '        op-mix-norm(H0): itt=', itt, ' ispin=', ispin, &
            ' n_basis=', n_basis, ' max|op_mix|=', op_mix_norm_h0
          flush(6)
        end if
        call zgemm('N', 'N', n_basis, dg_frag%nstate_tot, n_basis, (1.0d0, 0.0d0), op_mix, n_basis, &
                   coef_mix_all, n_basis, (0.0d0, 0.0d0), rhs_mix, n_basis)
        call expand_mixed_rhs_to_raw(rhs_mix, raw_rhs)
        dcoef_dt_h0(:, :) = raw_rhs(:, :)
        if (need_source_decomp) dcoef_dt_h0_core(:, :) = raw_rhs(:, :)
        if (has_nonlocal) then
          if (allocated(dg_frag%H_nl_blocks)) then
            if (need_source_decomp) then
              call apply_block_nonlocal_cache(dcoef_dt_nl)
              dcoef_dt_h0(1:n_frag, :) = dcoef_dt_h0(1:n_frag, :) + dcoef_dt_nl(1:n_frag, :)
            else
              call apply_block_nonlocal_cache(dcoef_dt_h0, .false.)
            end if
          else
            call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all(1:n_frag, :), dcoef_dt_h0(1:n_frag, :))
            if (need_source_decomp) &
              call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all(1:n_frag, :), dcoef_dt_nl(1:n_frag, :))
          end if
        end if
        if (has_so_nonlocal) then
          call apply_nonlocal_pp_projector_batch_so(dg_frag, mg, ppg, system, Ac_tot, ispin, &
               coef_all(1:n_frag, :), coef_frag_other(1:n_frag, 1:dg_frag%nstate_tot), dcoef_dt_h0(1:n_frag, :))
          if (need_source_decomp) &
            call apply_nonlocal_pp_projector_batch_so(dg_frag, mg, ppg, system, Ac_tot, ispin, &
                 coef_all(1:n_frag, :), coef_frag_other(1:n_frag, 1:dg_frag%nstate_tot), dcoef_dt_nl(1:n_frag, :))
        end if
        ! NOTE:
        ! In mixed-basis mode, H0c already includes the diamagnetic A^2/2 term
        ! (either from uniform A_squared on the diagonal or from A2_mat in the
        ! spatial-A path). Adding it again here double-counts A^2.
        if (need_source_decomp) then
          dcoef_dt_h0_core(1:n_tot, :) = -zi * dcoef_dt_h0_core(1:n_tot, :)
          dcoef_dt_nl(1:n_tot, :) = -zi * dcoef_dt_nl(1:n_tot, :)
        end if
        dcoef_dt_h0(1:n_tot, :) = -zi * dcoef_dt_h0(1:n_tot, :)
      else if (n_pw > 0 .and. allocated(dg_frag%H_mat_frag_pw)) then

        call apply_mixed_hamiltonian(dg_frag, ispin, coef_all(1:n_tot, :), dcoef_dt_h0(1:n_tot, :))
        if (need_source_decomp) dcoef_dt_h0_core(1:n_tot, :) = dcoef_dt_h0(1:n_tot, :)

        if (has_nonlocal) then
          if (allocated(dg_frag%H_nl_blocks)) then
            if (need_source_decomp) then
              call apply_block_nonlocal_cache(dcoef_dt_nl)
              dcoef_dt_h0(1:n_frag, :) = dcoef_dt_h0(1:n_frag, :) + dcoef_dt_nl(1:n_frag, :)
            else
              call apply_block_nonlocal_cache(dcoef_dt_h0, .false.)
            end if
          else
            call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all(1:n_frag, :), dcoef_dt_h0(1:n_frag, :))
            if (need_source_decomp) &
              call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all(1:n_frag, :), dcoef_dt_nl(1:n_frag, :))
          end if
        end if
        if (has_so_nonlocal) then
          call apply_nonlocal_pp_projector_batch_so(dg_frag, mg, ppg, system, Ac_tot, ispin, &
               coef_all(1:n_frag, :), coef_frag_other(1:n_frag, 1:dg_frag%nstate_tot), dcoef_dt_h0(1:n_frag, :))
          if (need_source_decomp) &
            call apply_nonlocal_pp_projector_batch_so(dg_frag, mg, ppg, system, Ac_tot, ispin, &
                 coef_all(1:n_frag, :), coef_frag_other(1:n_frag, 1:dg_frag%nstate_tot), dcoef_dt_nl(1:n_frag, :))
        end if
        if (need_source_decomp) then
          dcoef_dt_h0_core(1:n_tot, :) = -zi * dcoef_dt_h0_core(1:n_tot, :)
          dcoef_dt_nl(1:n_tot, :) = -zi * dcoef_dt_nl(1:n_tot, :)
        end if
        dcoef_dt_h0(1:n_tot, :) = -zi * dcoef_dt_h0(1:n_tot, :)
      else if (n_pw > 0 .and. allocated(dg_frag%H_mat_blocks) .and. .not. disable_block_h_apply) then

        if (allocated(dg_frag%H_local_block_ids)) then
          call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, coef_all(1:n_frag, :), &
                                         dcoef_dt_h0(1:n_frag, :), dg_frag%H_local_block_ids)
        else
          call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, coef_all(1:n_frag, :), dcoef_dt_h0(1:n_frag, :))
        end if
        if (allocated(dg_frag%H_mat_pw_diag)) then
          do io = 1, n_pw
            dcoef_dt_h0(n_frag+io, :) = dcoef_dt_h0(n_frag+io, :) + &
              dg_frag%H_mat_pw_diag(io, ispin) * coef_all(n_frag+io, :)
          end do
        else
          do io = 1, n_pw
            dcoef_dt_h0(n_frag+io, :) = dcoef_dt_h0(n_frag+io, :) + &
              0.5d0 * sum(dg_frag%k_pw(1:3, io)**2) * coef_all(n_frag+io, :)
          end do
        end if
        if (need_source_decomp) dcoef_dt_h0_core(1:n_tot, :) = dcoef_dt_h0(1:n_tot, :)

        if (has_nonlocal) then
          if (allocated(dg_frag%H_nl_blocks)) then
            if (need_source_decomp) then
              call apply_block_nonlocal_cache(dcoef_dt_nl)
              dcoef_dt_h0(1:n_frag, :) = dcoef_dt_h0(1:n_frag, :) + dcoef_dt_nl(1:n_frag, :)
            else
              call apply_block_nonlocal_cache(dcoef_dt_h0, .false.)
            end if
          else
            call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all(1:n_frag, :), dcoef_dt_h0(1:n_frag, :))
            if (need_source_decomp) &
              call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all(1:n_frag, :), dcoef_dt_nl(1:n_frag, :))
          end if
        end if
        if (has_so_nonlocal) then
          call apply_nonlocal_pp_projector_batch_so(dg_frag, mg, ppg, system, Ac_tot, ispin, &
               coef_all(1:n_frag, :), coef_frag_other(1:n_frag, 1:dg_frag%nstate_tot), dcoef_dt_h0(1:n_frag, :))
          if (need_source_decomp) &
            call apply_nonlocal_pp_projector_batch_so(dg_frag, mg, ppg, system, Ac_tot, ispin, &
                 coef_all(1:n_frag, :), coef_frag_other(1:n_frag, 1:dg_frag%nstate_tot), dcoef_dt_nl(1:n_frag, :))
        end if
        if (need_source_decomp) then
          dcoef_dt_h0_core(1:n_tot, :) = -zi * dcoef_dt_h0_core(1:n_tot, :)
          dcoef_dt_nl(1:n_tot, :) = -zi * dcoef_dt_nl(1:n_tot, :)
        end if
        dcoef_dt_h0(1:n_tot, :) = -zi * dcoef_dt_h0(1:n_tot, :)
      else if (n_pw == 0 .and. .not. disable_block_h_apply .and. .not. use_hmat_complex .and. allocated(dg_frag%H_mat_blocks)) then

        if (trace_deriv_step) t_h0_part0 = get_wtime()
        if (allocated(dg_frag%H_local_block_ids)) then
          call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, coef_all, dcoef_dt_h0, dg_frag%H_local_block_ids)
        else
          call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, coef_all, dcoef_dt_h0)
        end if
        if (allocated(dg_frag%H_mat_core_blocks) .and. enable_overlap_path_trace .and. &
            (itt <= 5 .or. mod(itt, 50) == 0)) then
          allocate(h_core_probe(n_tot, dg_frag%nstate_tot))
          h_core_probe(:, :) = (0.0d0, 0.0d0)
          if (allocated(dg_frag%H_local_block_ids)) then
            call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_core_blocks, ispin, coef_all, h_core_probe, &
                                           dg_frag%H_local_block_ids)
          else
            call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_core_blocks, ispin, coef_all, h_core_probe)
          end if
          h_split_local(1) = sum(abs(dcoef_dt_h0(:, :))**2)
          h_split_local(2) = sum(abs(h_core_probe(:, :))**2)
          h_split_local(3) = sum(abs(dcoef_dt_h0(:, :) - h_core_probe(:, :))**2)
          h_split_local(4) = sum(abs(coef_all(:, :))**2)
          h_split_global(:) = 0.0d0
          call comm_summation(h_split_local, h_split_global, 4, dg_frag%icomm)
          if (dg_frag%id == 0) then
            write(*,'(1x,a,i0,a,i0,a,4(1x,1pe14.6))') '[DG-H-SPLIT] itt=', itt, &
              ' ispin=', ispin, ' sqrtFull sqrtCoreRef sqrtFullMinusCoreRef sqrtCoef=', &
              sqrt(max(0.0d0, h_split_global(1))), sqrt(max(0.0d0, h_split_global(2))), &
              sqrt(max(0.0d0, h_split_global(3))), sqrt(max(0.0d0, h_split_global(4)))
            flush(6)
          end if
          deallocate(h_core_probe)
        end if
        if (trace_deriv_step) then
          t_h0_part1 = get_wtime()
          write(*,'(1x,a,i0,a,i0,a,1pe12.4)') '[DG-DERIV] H0 block apply done call=', &
            derivative_step_trace_id, ' ispin=', ispin, ' time=', t_h0_part1 - t_h0_part0
          call trace_deriv_norm('h0-block-raw', dcoef_dt_h0)
          call trace_h_blocks_summary()
          flush(6)
        end if
        if (need_source_decomp) dcoef_dt_h0_core(:, :) = dcoef_dt_h0(:, :)

        if (has_nonlocal) then
          if (allocated(dg_frag%H_nl_blocks)) then
            if (trace_deriv_step) then
              t_h0_part0 = get_wtime()
              write(*,'(1x,a,i0,a,i0,a)') '[DG-DERIV] H0 nonlocal start call=', &
                derivative_step_trace_id, ' ispin=', ispin, ' route=block-cache'
              flush(6)
            end if
            if (need_source_decomp) then
              call apply_block_nonlocal_cache(dcoef_dt_nl)
              dcoef_dt_h0(1:n_frag, :) = dcoef_dt_h0(1:n_frag, :) + dcoef_dt_nl(1:n_frag, :)
            else
              call apply_block_nonlocal_cache(dcoef_dt_h0, .false.)
            end if
            if (trace_deriv_step) then
              t_h0_part1 = get_wtime()
              write(*,'(1x,a,i0,a,i0,a,1pe12.4)') '[DG-DERIV] H0 nonlocal done call=', &
                derivative_step_trace_id, ' ispin=', ispin, ' time=', t_h0_part1 - t_h0_part0
              call trace_deriv_norm('h0-plus-nl-raw', dcoef_dt_h0)
              flush(6)
            end if
          else
            if (trace_deriv_step) then
              t_h0_part0 = get_wtime()
              write(*,'(1x,a,i0,a,i0,a)') '[DG-DERIV] H0 nonlocal start call=', &
                derivative_step_trace_id, ' ispin=', ispin, ' route=projector-direct'
              flush(6)
            end if
            call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all, dcoef_dt_h0)
            if (need_source_decomp) &
              call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all, dcoef_dt_nl)
            if (trace_deriv_step) then
              t_h0_part1 = get_wtime()
              write(*,'(1x,a,i0,a,i0,a,1pe12.4)') '[DG-DERIV] H0 nonlocal done call=', &
                derivative_step_trace_id, ' ispin=', ispin, ' time=', t_h0_part1 - t_h0_part0
              call trace_deriv_norm('h0-plus-nl-raw', dcoef_dt_h0)
              flush(6)
            end if
          end if
        end if
        if (has_so_nonlocal) then
          if (trace_deriv_step) t_h0_part0 = get_wtime()
          call apply_nonlocal_pp_projector_batch_so(dg_frag, mg, ppg, system, Ac_tot, ispin, &
               coef_all(1:n_frag, :), coef_frag_other(1:n_frag, 1:dg_frag%nstate_tot), dcoef_dt_h0(1:n_frag, :))
          if (need_source_decomp) &
            call apply_nonlocal_pp_projector_batch_so(dg_frag, mg, ppg, system, Ac_tot, ispin, &
                 coef_all(1:n_frag, :), coef_frag_other(1:n_frag, 1:dg_frag%nstate_tot), dcoef_dt_nl(1:n_frag, :))
          if (trace_deriv_step) then
            t_h0_part1 = get_wtime()
            write(*,'(1x,a,i0,a,i0,a,1pe12.4)') '[DG-DERIV] H0 so-nonlocal done call=', &
              derivative_step_trace_id, ' ispin=', ispin, ' time=', t_h0_part1 - t_h0_part0
            flush(6)
          end if
        end if
        if (need_source_decomp) then
          dcoef_dt_h0_core(1:n_frag, :) = -zi * dcoef_dt_h0_core(1:n_frag, :)
          dcoef_dt_nl(1:n_frag, :) = -zi * dcoef_dt_nl(1:n_frag, :)
        end if
        dcoef_dt_h0(1:n_frag, :) = -zi * dcoef_dt_h0(1:n_frag, :)
      else

        write(*,'(1x,a,i0,a,i0,a,i0)') '[FATAL] derivative reached removed dense H0 application path: rank=', &
             dg_frag%id, ' itt=', itt, ' ispin=', ispin
        stop 1

      end if

      if (has_spin_mix .and. n_frag > 0) then
        dcoef_dt_h0(1:n_frag, 1:dg_frag%nstate_tot) = dcoef_dt_h0(1:n_frag, 1:dg_frag%nstate_tot) + &
          dcoef_dt_spinmix(1:n_frag, 1:dg_frag%nstate_tot)
      end if
      if (trace_deriv_step) then
        t_trace1 = get_wtime()
        write(*,'(1x,a,i0,a,i0,a,1pe12.4)') '[DG-DERIV] H0 apply done call=', &
          derivative_step_trace_id, ' ispin=', ispin, ' time=', t_trace1 - t_trace0
        call trace_deriv_norm('h0', dcoef_dt_h0)
        flush(6)
      end if
      if (enable_deriv_trace .and. itt <= 2) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        deriv-trace: rank=", dg_frag%id, " itt=", itt, &
          " ispin=", ispin, " stage=after-h0-term"
        flush(6)
      end if


      if (use_mixed_basis .and. n_basis > 0) then
        if (trace_deriv_step) then
          t_trace0 = get_wtime()
          write(*,'(1x,a,i0,a,i0,a)') '[DG-DERIV] M apply start call=', &
            derivative_step_trace_id, ' ispin=', ispin, ' route=mixed'
          flush(6)
        end if
        rhs_mix(:, :) = (0.0d0, 0.0d0)
        call build_mixed_projected_operator(M, op_mix)
        if (enable_m_block_audit .and. dg_frag%id == 0 .and. itt <= 5) then
          m_raw_ff = 0.0d0
          m_raw_fp = 0.0d0
          m_raw_pf = 0.0d0
          m_raw_pp = 0.0d0
          if (n_frag > 0) m_raw_ff = sqrt(sum(abs(M(1:n_frag, 1:n_frag))**2))
          if (n_frag > 0 .and. n_pw > 0) then
            m_raw_fp = sqrt(sum(abs(M(1:n_frag, n_frag+1:n_tot))**2))
            m_raw_pf = sqrt(sum(abs(M(n_frag+1:n_tot, 1:n_frag))**2))
          end if
          if (n_pw > 0) m_raw_pp = sqrt(sum(abs(M(n_frag+1:n_tot, n_frag+1:n_tot))**2))
          m_raw_all = sqrt(sum(abs(M(1:n_tot, 1:n_tot))**2))

          nocc_basis = min(n_basis, dg_frag%nstate_tot)
          if (allocated(dg_frag%nocc_spin)) then
            if (ispin <= size(dg_frag%nocc_spin)) nocc_basis = min(nocc_basis, max(0, dg_frag%nocc_spin(ispin)))
          end if
          nvirt_basis = max(0, n_basis - nocc_basis)

          m_mix_oo = 0.0d0
          m_mix_ov = 0.0d0
          m_mix_vo = 0.0d0
          m_mix_vv = 0.0d0
          if (nocc_basis > 0) m_mix_oo = sqrt(sum(abs(op_mix(1:nocc_basis, 1:nocc_basis))**2))
          if (nocc_basis > 0 .and. nvirt_basis > 0) then
            m_mix_ov = sqrt(sum(abs(op_mix(1:nocc_basis, nocc_basis+1:n_basis))**2))
            m_mix_vo = sqrt(sum(abs(op_mix(nocc_basis+1:n_basis, 1:nocc_basis))**2))
          end if
          if (nvirt_basis > 0) m_mix_vv = sqrt(sum(abs(op_mix(nocc_basis+1:n_basis, nocc_basis+1:n_basis))**2))
          m_mix_all = sqrt(sum(abs(op_mix(1:n_basis, 1:n_basis))**2))

          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') '[M-BLOCK-AUDIT] itt=', itt, ' ispin=', ispin, &
               ' n_frag=', n_frag, ' n_pw=', n_pw, ' n_basis=', n_basis
          write(*,'(1x,a,i0,a,i0)') '[M-BLOCK-AUDIT] nocc_basis=', nocc_basis, ' nvirt_basis=', nvirt_basis
          write(*,'(1x,a,5(1x,1pe14.6))') '[M-BLOCK-AUDIT] raw_norms ff fp pf pp all=', &
               m_raw_ff, m_raw_fp, m_raw_pf, m_raw_pp, m_raw_all
          write(*,'(1x,a,5(1x,1pe14.6))') '[M-BLOCK-AUDIT] mix_norms oo ov vo vv all=', &
               m_mix_oo, m_mix_ov, m_mix_vo, m_mix_vv, m_mix_all
          flush(6)
        end if
        if (enable_op_mix_trace .and. dg_frag%id == 0 .and. (itt <= 200 .or. mod(itt, 50) == 0)) then
          op_mix_norm_m = maxval(abs(op_mix(1:n_basis, 1:n_basis)))
          write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4)') '        op-mix-norm(M): itt=', itt, ' ispin=', ispin, &
            ' n_basis=', n_basis, ' max|op_mix|=', op_mix_norm_m
          flush(6)
        end if
        if (enable_op_mix_trace .and. dg_frag%id == 0 .and. itt <= 5 .and. n_basis > 0) then
          write(*,'(1x,a,i0,a,i0,a,i0)') '[MIX-M-ROWNORM] itt=', itt, ' ispin=', ispin, ' n_basis=', n_basis
          write(*,'(1x,a)') '[MIX-M-ROWNORM] mode  row_norm_off_diag  row_max_off_diag'
          do io = 1, n_basis
            op_mix_norm_h0 = 0.0d0
            op_mix_norm_m = 0.0d0
            do jo = 1, n_basis
              if (io == jo) cycle
              op_mix_norm_h0 = op_mix_norm_h0 + abs(op_mix(io,jo))**2
              op_mix_norm_m  = max(op_mix_norm_m, abs(op_mix(io,jo)))
            end do
            write(*,'(1x,a,i5,2(1x,1pe14.6))') '[MIX-M-ROWNORM] ', io, sqrt(op_mix_norm_h0), op_mix_norm_m
          end do
          flush(6)
        end if
        call zgemm('N', 'N', n_basis, dg_frag%nstate_tot, n_basis, (1.0d0, 0.0d0), op_mix, n_basis, &
                   coef_mix_all, n_basis, (0.0d0, 0.0d0), rhs_mix, n_basis)
        call expand_mixed_rhs_to_raw(rhs_mix, raw_rhs)
        dcoef_dt_m(:, :) = raw_rhs(:, :)
      else if (allocated(dg_frag%momentum_blocks) .and. .not. use_spatial_A) then
        if (trace_deriv_step) then
          t_trace0 = get_wtime()
          write(*,'(1x,a,i0,a,i0,a,i0)') '[DG-DERIV] M apply start call=', &
            derivative_step_trace_id, ' ispin=', ispin, ' n_pw=', n_pw
          flush(6)
        end if
        dcoef_dt_m(:, :) = (0.0d0, 0.0d0)
        call apply_momentum_blocks(dg_frag, ispin, Ac_tot, coef_all(1:n_frag, :), dcoef_dt_m(1:n_frag, :))
        if (enable_ap_block_check .and. n_frag > 0 .and. dg_frag%nstate_tot > 0 .and. itt <= 2) then
          if (.not. allocated(ap_block_dense)) then
            allocate(ap_block_dense(n_frag, n_frag))
          else if (size(ap_block_dense, 1) /= n_frag .or. size(ap_block_dense, 2) /= n_frag) then
            deallocate(ap_block_dense)
            allocate(ap_block_dense(n_frag, n_frag))
          end if
          if (.not. allocated(ap_block_ref)) then
            allocate(ap_block_ref(n_frag, dg_frag%nstate_tot))
          else if (size(ap_block_ref, 1) /= n_frag .or. size(ap_block_ref, 2) /= dg_frag%nstate_tot) then
            deallocate(ap_block_ref)
            allocate(ap_block_ref(n_frag, dg_frag%nstate_tot))
          end if
          call copy_momentum_blocks_to_complex_dense(dg_frag, ispin, Ac_tot, ap_block_dense)
          call zgemm('N', 'N', n_frag, dg_frag%nstate_tot, n_frag, (1.0d0, 0.0d0), ap_block_dense, n_frag, &
                     coef_all, n_tot, (0.0d0, 0.0d0), ap_block_ref, n_frag)
          ap_block_diff_max = maxval(abs(ap_block_ref - dcoef_dt_m(1:n_frag, :)))
          if (dg_frag%id == 0) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,1pe14.6)') "        ap-block-check: rank=", dg_frag%id, " itt=", itt, &
              " ispin=", ispin, " max|block-dense|=", ap_block_diff_max
            flush(6)
          end if
        end if
        if (n_pw > 0) then
          if (.not. allocated(mfp_vec)) then
            allocate(mfp_vec(n_pw))
          else if (size(mfp_vec, 1) /= n_pw) then
            deallocate(mfp_vec)
            allocate(mfp_vec(n_pw))
          end if
!$omp parallel do schedule(static) private(io)
          do io = 1, n_pw
            mfp_vec(io) = zi * dot_product(Ac_tot(1:3), dg_frag%k_pw(1:3, io))
            dcoef_dt_m(n_frag+io, :) = dcoef_dt_m(n_frag+io, :) + mfp_vec(io) * coef_all(n_frag+io, :)
          end do
!$omp end parallel do
          if (.not. disable_mfp .and. mixed_fp_coupling_active(dg_frag, ispin)) then
            if (use_direct_pfp .and. ((use_fft_pfp .and. allocated(dg_frag%P_mat_frag_pw_g)) .or. allocated(dg_frag%P_mat_frag_pw))) then
              if (.not. allocated(tmp_fp_bw)) then
                allocate(tmp_fp_bw(n_pw, min(state_block_size, dg_frag%nstate_tot)))
              else if (size(tmp_fp_bw, 1) /= n_pw .or. size(tmp_fp_bw, 2) /= min(state_block_size, dg_frag%nstate_tot)) then
                deallocate(tmp_fp_bw)
                allocate(tmp_fp_bw(n_pw, min(state_block_size, dg_frag%nstate_tot)))
              end if
              if (.not. allocated(tmp_fp_fw)) then
                allocate(tmp_fp_fw(n_frag, min(state_block_size, dg_frag%nstate_tot)))
              else if (size(tmp_fp_fw, 1) /= n_frag .or. size(tmp_fp_fw, 2) /= min(state_block_size, dg_frag%nstate_tot)) then
                deallocate(tmp_fp_fw)
                allocate(tmp_fp_fw(n_frag, min(state_block_size, dg_frag%nstate_tot)))
              end if
              if (.not. allocated(mfp_fp)) then
                allocate(mfp_fp(n_frag, n_pw))
              else if (size(mfp_fp, 1) /= n_frag .or. size(mfp_fp, 2) /= n_pw) then
                deallocate(mfp_fp)
                allocate(mfp_fp(n_frag, n_pw))
              end if

              do jo = 1, n_pw
                if (use_fft_pfp .and. allocated(dg_frag%P_mat_frag_pw_g)) then
                  mfp_fp(1:n_frag, jo) = Ac_tot(1) * dg_frag%P_mat_frag_pw_g(1, 1:n_frag, jo, ispin) + &
                                         Ac_tot(2) * dg_frag%P_mat_frag_pw_g(2, 1:n_frag, jo, ispin) + &
                                         Ac_tot(3) * dg_frag%P_mat_frag_pw_g(3, 1:n_frag, jo, ispin)
                else
                  mfp_fp(1:n_frag, jo) = Ac_tot(1) * dg_frag%P_mat_frag_pw(1, 1:n_frag, jo, ispin) + &
                                         Ac_tot(2) * dg_frag%P_mat_frag_pw(2, 1:n_frag, jo, ispin) + &
                                         Ac_tot(3) * dg_frag%P_mat_frag_pw(3, 1:n_frag, jo, ispin)
                end if
              end do

              do state0 = 1, dg_frag%nstate_tot, state_block_size
                nstate_blk = min(state_block_size, dg_frag%nstate_tot - state0 + 1)
                state_s = state0
                state_e = state0 + nstate_blk - 1

                call zgemm('N', 'N', n_frag, nstate_blk, n_pw, (1.0d0, 0.0d0), &
                           mfp_fp, n_frag, coef_all(n_frag+1:n_tot, state_s:state_e), n_pw, &
                           (0.0d0, 0.0d0), tmp_fp_fw, n_frag)
                dcoef_dt_m(1:n_frag, state_s:state_e) = dcoef_dt_m(1:n_frag, state_s:state_e) + &
                  tmp_fp_fw(1:n_frag, 1:nstate_blk)

                call zgemm('C', 'N', n_pw, nstate_blk, n_frag, (1.0d0, 0.0d0), &
                           mfp_fp, n_frag, coef_all(1, state_s), n_tot, (0.0d0, 0.0d0), tmp_fp_bw, n_pw)
                dcoef_dt_m(n_frag+1:n_tot, state_s:state_e) = dcoef_dt_m(n_frag+1:n_tot, state_s:state_e) - &
                  tmp_fp_bw(1:n_pw, 1:nstate_blk)
              end do
            else
              if (.not. allocated(coef_pw_scaled)) then
                allocate(coef_pw_scaled(n_pw, min(state_block_size, dg_frag%nstate_tot)), &
                         tmp_fp_bw(n_pw, min(state_block_size, dg_frag%nstate_tot)))
              else if (size(coef_pw_scaled, 1) /= n_pw .or. size(coef_pw_scaled, 2) /= min(state_block_size, dg_frag%nstate_tot)) then
                deallocate(coef_pw_scaled, tmp_fp_bw)
                allocate(coef_pw_scaled(n_pw, min(state_block_size, dg_frag%nstate_tot)), &
                         tmp_fp_bw(n_pw, min(state_block_size, dg_frag%nstate_tot)))
              end if
              if (.not. allocated(tmp_fp_fw)) then
                allocate(tmp_fp_fw(n_frag, min(state_block_size, dg_frag%nstate_tot)))
              else if (size(tmp_fp_fw, 1) /= n_frag .or. size(tmp_fp_fw, 2) /= min(state_block_size, dg_frag%nstate_tot)) then
                deallocate(tmp_fp_fw)
                allocate(tmp_fp_fw(n_frag, min(state_block_size, dg_frag%nstate_tot)))
              end if

              do state0 = 1, dg_frag%nstate_tot, state_block_size
                nstate_blk = min(state_block_size, dg_frag%nstate_tot - state0 + 1)
                state_s = state0
                state_e = state0 + nstate_blk - 1

!$omp parallel do schedule(static) private(jo)
                do jo = 1, n_pw
                  coef_pw_scaled(jo, 1:nstate_blk) = mfp_vec(jo) * coef_all(n_frag+jo, state_s:state_e)
                end do
!$omp end parallel do

                call zgemm('N', 'N', n_frag, nstate_blk, n_pw, (1.0d0, 0.0d0), &
                           dg_frag%S_mat_frag_pw(1:n_frag, 1:n_pw, ispin), n_frag, &
                           coef_pw_scaled, n_pw, (0.0d0, 0.0d0), tmp_fp_fw, n_frag)
                dcoef_dt_m(1:n_frag, state_s:state_e) = dcoef_dt_m(1:n_frag, state_s:state_e) + &
                  tmp_fp_fw(1:n_frag, 1:nstate_blk)

                call zgemm('C', 'N', n_pw, nstate_blk, n_frag, (1.0d0, 0.0d0), &
                           dg_frag%S_mat_frag_pw(1:n_frag, 1:n_pw, ispin), n_frag, &
                           coef_all(1, state_s), n_tot, (0.0d0, 0.0d0), tmp_fp_bw, n_pw)
!$omp parallel do schedule(static) private(jo)
                do jo = 1, n_pw
                  dcoef_dt_m(n_frag+jo, state_s:state_e) = dcoef_dt_m(n_frag+jo, state_s:state_e) - &
                    conjg(mfp_vec(jo)) * tmp_fp_bw(jo, 1:nstate_blk)
                end do
!$omp end parallel do
              end do
            end if
          end if
        end if
      else
        write(*,'(1x,a,i0,a,i0,a,i0)') '[FATAL] derivative reached removed dense M application path: rank=', &
             dg_frag%id, ' itt=', itt, ' ispin=', ispin
        stop 1
      end if
      if (enable_m_block_audit .and. .not. (use_mixed_basis .and. n_basis > 0) .and. &
          allocated(dg_frag%momentum_blocks) .and. dg_frag%id == 0 .and. itt <= 5) then
        m_raw_ff = 0.0d0
        max_abs_m = 0.0d0
        do iblk = 1, size(dg_frag%momentum_blocks)
          if (.not. allocated(dg_frag%momentum_blocks(iblk)%val)) cycle
          if (ispin > size(dg_frag%momentum_blocks(iblk)%val, 4)) cycle
          nrow_blk = min(dg_frag%momentum_blocks(iblk)%nrow_max, size(dg_frag%momentum_blocks(iblk)%val, 2))
          ncol_blk = min(dg_frag%momentum_blocks(iblk)%ncol_max, size(dg_frag%momentum_blocks(iblk)%val, 3))
          if (nrow_blk <= 0 .or. ncol_blk <= 0) cycle
          m_raw_ff = m_raw_ff + sum(( &
            Ac_tot(1) * dg_frag%momentum_blocks(iblk)%val(1, 1:nrow_blk, 1:ncol_blk, ispin) + &
            Ac_tot(2) * dg_frag%momentum_blocks(iblk)%val(2, 1:nrow_blk, 1:ncol_blk, ispin) + &
            Ac_tot(3) * dg_frag%momentum_blocks(iblk)%val(3, 1:nrow_blk, 1:ncol_blk, ispin))**2)
          max_abs_m = max(max_abs_m, maxval(abs( &
            Ac_tot(1) * dg_frag%momentum_blocks(iblk)%val(1, 1:nrow_blk, 1:ncol_blk, ispin) + &
            Ac_tot(2) * dg_frag%momentum_blocks(iblk)%val(2, 1:nrow_blk, 1:ncol_blk, ispin) + &
            Ac_tot(3) * dg_frag%momentum_blocks(iblk)%val(3, 1:nrow_blk, 1:ncol_blk, ispin))))
        end do
        m_rhs_frag_norm = 0.0d0
        m_rhs_pw_norm = 0.0d0
        if (n_frag > 0) m_rhs_frag_norm = sqrt(sum(abs(dcoef_dt_m(1:n_frag, :))**2))
        if (n_pw > 0) m_rhs_pw_norm = sqrt(sum(abs(dcoef_dt_m(n_frag+1:n_tot, :))**2))
        m_rhs_all_norm = sqrt(sum(abs(dcoef_dt_m(1:n_tot, :))**2))
        m_rhs_max = maxval(abs(dcoef_dt_m(1:n_tot, :)))
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') '[M-BLOCK-AUDIT] itt=', itt, ' ispin=', ispin, &
             ' n_frag=', n_frag, ' n_pw=', n_pw, ' n_basis=', n_basis
        write(*,'(1x,a,2(1x,1pe14.6))') '[M-BLOCK-AUDIT] route=block raw_frob raw_max=', &
             sqrt(max(0.0d0, m_raw_ff)), max_abs_m
        write(*,'(1x,a,4(1x,1pe14.6))') '[M-BLOCK-AUDIT] applied_norm frag pw all max=', &
             m_rhs_frag_norm, m_rhs_pw_norm, m_rhs_all_norm, m_rhs_max
        flush(6)
      end if
      if (trace_deriv_step) then
        t_trace1 = get_wtime()
        write(*,'(1x,a,i0,a,i0,a,1pe12.4)') '[DG-DERIV] M apply done call=', &
          derivative_step_trace_id, ' ispin=', ispin, ' time=', t_trace1 - t_trace0
        call trace_deriv_norm('m', dcoef_dt_m)
        flush(6)
      end if
      if (enable_deriv_trace .and. itt <= 2) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        deriv-trace: rank=", dg_frag%id, " itt=", itt, &
          " ispin=", ispin, " stage=after-m-term"
        flush(6)
      end if


      rhs_all = dcoef_dt_h0 - dcoef_dt_m
      if (trace_deriv_step) call trace_deriv_norm('rhs-pre-s', rhs_all)
      if (need_source_decomp) rhs_pre_solve(:, :) = rhs_all(:, :)

      mfp_coupling_on = (n_pw > 0 .and. mixed_fp_coupling_active(dg_frag, ispin))
      has_s_prop_blocks_c = allocated(dg_frag%S_mat_prop_blocks_c)
      has_s_blocks_c = allocated(dg_frag%S_mat_blocks_c)
      has_s_prop_blocks = allocated(dg_frag%S_mat_prop_blocks)
      has_s_blocks = allocated(dg_frag%S_mat_blocks)
      has_s_prop_c = allocated(dg_frag%S_mat_prop_c)
      has_s_c = allocated(dg_frag%S_mat_c)
      overlap_gate_reason = 'none'
      n_s = 0
      ! Decide whether the RHS needs an S^{-1} application.  The mixed-basis
      ! route normally propagates in the orthonormalized mixed coordinates, so
      ! applying the fragment overlap again is skipped unless explicitly forced.
      ! Legacy fragment-only routes still use the available overlap storage.
      if (use_mixed_basis .and. n_basis > 0) then
        if (force_overlap_solve) then
          if (mfp_coupling_on) then
            n_s = n_tot
            overlap_gate_reason = 'force+mixed-fp-coupling'
          else if (has_s_prop_blocks_c .or. has_s_blocks_c .or. has_s_prop_blocks .or. has_s_blocks .or. has_s_prop_c .or. has_s_c) then
            n_s = n_frag
            overlap_gate_reason = 'force+fragment-overlap'
          else
            n_s = 0
            overlap_gate_reason = 'force+no-overlap-storage'
          end if
        else
          n_s = 0
          overlap_gate_reason = 'mixed-basis-priority'
        end if
      else if (mfp_coupling_on) then
        n_s = n_tot
        overlap_gate_reason = 'mixed-fp-coupling'
      else if (has_s_prop_blocks_c .or. has_s_blocks_c .or. has_s_prop_blocks .or. has_s_blocks .or. has_s_prop_c .or. has_s_c) then
        n_s = n_frag
        overlap_gate_reason = 'fragment-overlap'
      else
        overlap_gate_reason = 'no-overlap-storage'
      end if

      if (trace_deriv_step) then
        t_trace0 = get_wtime()
        write(*,'(1x,a,i0,a,i0,a,i0,a,a)') '[DG-DERIV] overlap start call=', &
          derivative_step_trace_id, ' ispin=', ispin, ' n_s=', n_s, &
          ' reason=', trim(overlap_gate_reason)
        flush(6)
      end if
      if (enable_overlap_path_trace .and. dg_frag%id == 0 .and. (itt <= 200 .or. mod(itt, 50) == 0)) then
        write(*,'(1x,a,i0,a,i0,a,5(i0,a),a,a)') '[OVERLAP-GATE] itt=', itt, ' ispin=', ispin, &
          ' n_frag=', n_frag, ' n_pw=', n_pw, ' n_tot=', n_tot, ' n_basis=', n_basis, ' n_s=', n_s, &
          ' reason=', trim(overlap_gate_reason)
        write(*,'(1x,a,8(a,l1),a,l1)') '[OVERLAP-GATE] flags:', &
          ' use_mixed=', use_mixed_basis, ' force=', force_overlap_solve, ' mixed_ready=', dg_frag%mixed_basis_ready, &
          ' mfp_coupling=', mfp_coupling_on, ' has_prop_blocks=', has_s_prop_blocks, ' has_blocks=', has_s_blocks, &
          ' has_prop_c=', has_s_prop_c, ' has_c=', has_s_c, ' use_prop=', use_prop_overlap
        flush(6)
      end if
      if (n_s > 0) then

        if (.not. allocated(rhs_in)) then
          allocate(rhs_in(n_s, dg_frag%nstate_tot))
        else if (size(rhs_in, 1) /= n_s .or. size(rhs_in, 2) /= dg_frag%nstate_tot) then
          deallocate(rhs_in)
          allocate(rhs_in(n_s, dg_frag%nstate_tot))
        end if
        rhs_in(:, :) = rhs_all(1:n_s, :)
        rhs_pre_norm = 0.0d0
        rhs_pre_max = 0.0d0
        if (enable_overlap_path_trace .and. (itt <= 5 .or. mod(itt, 50) == 0)) then
          rhs_pre_norm = sqrt(sum(abs(rhs_in(:, :))**2))
          rhs_pre_max = maxval(abs(rhs_in(:, :)))
        end if
        if (bypass_overlap_solve) then
          if (enable_overlap_path_trace .and. dg_frag%id == 0 .and. (itt <= 200 .or. mod(itt, 50) == 0)) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,l1,a,l1,a,l1)') '[OVERLAP-PATH] bypass: itt=', itt, ' ispin=', ispin, &
              ' n_s=', n_s, ' use_prop=', use_prop_overlap, ' force=', force_overlap_solve, ' mixed=', use_mixed_basis
            flush(6)
          end if
          rhs_all(1:n_s, :) = rhs_in(:, :)
        else
          ! Convert dC/dt from covariant RHS form to coefficient derivatives
          ! by solving S x = rhs over the active fragment/mixed subspace.
          if (enable_overlap_path_trace .and. dg_frag%id == 0 .and. (itt <= 200 .or. mod(itt, 50) == 0)) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,l1,a,l1,a,l1)') '[OVERLAP-SOLVE-ENTER] itt=', itt, ' ispin=', ispin, &
              ' n_s=', n_s, ' use_prop=', use_prop_overlap, ' force=', force_overlap_solve, ' mixed=', use_mixed_basis
            flush(6)
          end if
          call solve_overlap_operator_batch(dg_frag, ispin, rhs_in, rhs_all(1:n_s, :), use_prop_overlap)
        end if
        if (enable_overlap_path_trace .and. (itt <= 5 .or. mod(itt, 50) == 0)) then
          rhs_post_norm = sqrt(sum(abs(rhs_all(1:n_s, :))**2))
          rhs_post_max = maxval(abs(rhs_all(1:n_s, :)))
          rhs_norm_local_arr = (/ rhs_pre_norm, rhs_pre_max, rhs_post_norm, rhs_post_max /)
          rhs_norm_max_arr = 0.0d0
          call comm_get_max(rhs_norm_local_arr, rhs_norm_max_arr, 4, dg_frag%icomm)
          rhs_pre_norm_loc%val = rhs_pre_norm
          rhs_pre_norm_loc%rank = dg_frag%id
          rhs_post_norm_loc%val = rhs_post_norm
          rhs_post_norm_loc%rank = dg_frag%id
          call comm_get_max(rhs_pre_norm_loc, rhs_pre_norm_maxloc, dg_frag%icomm)
          call comm_get_max(rhs_post_norm_loc, rhs_post_norm_maxloc, dg_frag%icomm)
          if (dg_frag%id == 0) then
            write(*,'(1x,a,i0,a,i0,a,4(1x,1pe14.6),a,4(1x,1pe14.6))') '[OVERLAP-NORM] itt=', itt, ' ispin=', ispin, &
              ' local_pre premax post postmax=', rhs_pre_norm, rhs_pre_max, rhs_post_norm, rhs_post_max, &
              ' global_max=', rhs_norm_max_arr
            write(*,'(1x,a,i0,a,i0,a,2(1x,i0,1x,1pe14.6))') '[OVERLAP-NORM-MAXLOC] itt=', itt, ' ispin=', ispin, &
              ' pre_rank/pre post_rank/post=', rhs_pre_norm_maxloc%rank, rhs_pre_norm_maxloc%val, &
              rhs_post_norm_maxloc%rank, rhs_post_norm_maxloc%val
            flush(6)
          end if
        end if
        if (enable_overlap_path_trace .and. (itt <= 5 .or. mod(itt, 50) == 0)) then
          nocc_spin = dg_frag%nstate_tot
          if (allocated(dg_frag%nocc_spin)) then
            if (ispin <= size(dg_frag%nocc_spin)) nocc_spin = min(nocc_spin, max(0, dg_frag%nocc_spin(ispin)))
          end if
          excite_local(:) = 0.0d0
          if (n_s > 0) then
            if (.not. allocated(vec_occ)) then
              allocate(vec_occ(n_s), svec_occ(n_s), sdc_occ(n_s), vec_src(n_s), ssrc_vec(n_s))
            else if (size(vec_occ, 1) /= n_s) then
              deallocate(vec_occ, svec_occ, sdc_occ, vec_src, ssrc_vec)
              allocate(vec_occ(n_s), svec_occ(n_s), sdc_occ(n_s), vec_src(n_s), ssrc_vec(n_s))
            end if
          end if
          do io_occ = 1, nocc_spin
            occ_factor = 1.0d0
            if (allocated(system%rocc)) then
              if (io_occ <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
                occ_factor = max(0.0d0, system%rocc(io_occ, 1, ispin))
              end if
            end if
            excite_local(1) = excite_local(1) + occ_factor * sum(abs(dcoef_dt_h0(1:n_s, io_occ))**2)
            excite_local(2) = excite_local(2) + occ_factor * sum(abs(dcoef_dt_m(1:n_s, io_occ))**2)
            excite_local(3) = excite_local(3) + occ_factor * sum(abs(rhs_in(1:n_s, io_occ))**2)
            excite_local(4) = excite_local(4) + occ_factor * sum(abs(rhs_all(1:n_s, io_occ))**2)
            if (dg_frag%id == 0) excite_local(5) = excite_local(5) + occ_factor
            eps_occ = 0.0d0
            if (allocated(dg_frag%esp)) then
              if (io_occ <= size(dg_frag%esp, 1) .and. ispin <= size(dg_frag%esp, 2)) eps_occ = dg_frag%esp(io_occ, ispin)
            end if
            vec_occ(:) = coef_all(1:n_s, io_occ)
            call apply_overlap_operator(dg_frag, ispin, vec_occ, svec_occ, use_prop_overlap)
            vec_src(:) = dcoef_dt_h0(1:n_s, io_occ) + zi * eps_occ * svec_occ(:)
            excite_local(6) = excite_local(6) + occ_factor * sum(abs(vec_src(:))**2)
            vec_src(:) = rhs_in(1:n_s, io_occ) + zi * eps_occ * svec_occ(:)
            excite_local(7) = excite_local(7) + occ_factor * sum(abs(vec_src(:))**2)
            vec_src(:) = rhs_all(1:n_s, io_occ) + zi * eps_occ * vec_occ(:)
            excite_local(8) = excite_local(8) + occ_factor * sum(abs(vec_src(:))**2)
          end do
          call comm_summation(excite_local, excite_global, 8, dg_frag%icomm)
          excite_per_electron = 0.0d0
          if (excite_global(5) > 1.0d-30) excite_per_electron = sqrt(max(0.0d0, excite_global(8)) / excite_global(5))
          if (dg_frag%id == 0) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,6(1x,1pe14.6))') '[DG-EXCITE-AUDIT] itt=', itt, ' ispin=', ispin, &
              ' nocc=', nocc_spin, ' occsum sqrtH sqrtM sqrtPre sqrtPost postResPerOcc=', &
              excite_global(5), sqrt(max(0.0d0, excite_global(1))), sqrt(max(0.0d0, excite_global(2))), &
              sqrt(max(0.0d0, excite_global(3))), sqrt(max(0.0d0, excite_global(4))), excite_per_electron
            write(*,'(1x,a,i0,a,i0,a,3(1x,1pe14.6))') '[DG-PHASE-RESIDUAL] itt=', itt, ' ispin=', ispin, &
              ' sqrtHres sqrtPreRes sqrtPostRes=', sqrt(max(0.0d0, excite_global(6))), &
              sqrt(max(0.0d0, excite_global(7))), sqrt(max(0.0d0, excite_global(8)))
            flush(6)
          end if
        end if

      else if (enable_overlap_path_trace .and. dg_frag%id == 0 .and. (itt <= 200 .or. mod(itt, 50) == 0)) then
        write(*,'(1x,a,i0,a,i0,a,a)') '[OVERLAP-PATH] no-solve: itt=', itt, ' ispin=', ispin, &
          ' reason=', trim(overlap_gate_reason)
        flush(6)
      end if
      if (trace_deriv_step) then
        t_trace1 = get_wtime()
        write(*,'(1x,a,i0,a,i0,a,1pe12.4)') '[DG-DERIV] overlap done call=', &
          derivative_step_trace_id, ' ispin=', ispin, ' time=', t_trace1 - t_trace0
        call trace_deriv_norm('rhs-post-s', rhs_all)
        flush(6)
      end if

      if (enable_norm_deriv_check .and. dg_frag%id == 0 .and. itt <= 3) then
        write(*,'(1x,a,i0,a,i0,a,3(1x,i0),a,2(1x,l1))') '        norm-deriv-gate: itt=', itt, ' ispin=', ispin, &
          ' n_frag/n_pw/n_s=', n_frag, n_pw, n_s, ' enforce/check=', enforce_norm_tangent, enable_norm_deriv_check
        if (.not. (n_pw == 0 .and. n_frag > 0)) then
          write(*,'(1x,a)') '        norm-deriv-gate: skipped because (n_pw == 0 .and. n_frag > 0) is false'
        end if
        flush(6)
      end if

      if (n_pw == 0 .and. n_frag > 0 .and. (enforce_norm_tangent .or. enable_norm_deriv_check .or. enable_excitation_source_trace)) then
        nocc_spin = dg_frag%nstate_tot
        if (allocated(dg_frag%nocc_spin)) then
          if (ispin <= size(dg_frag%nocc_spin)) nocc_spin = min(nocc_spin, max(0, dg_frag%nocc_spin(ispin)))
        end if
        if (nocc_spin > 0) then
          if (n_s > 0) then
            if (.not. allocated(vec_occ)) then
              allocate(vec_occ(n_s), svec_occ(n_s), sdc_occ(n_s), vec_src(n_s), ssrc_vec(n_s))
            else if (size(vec_occ, 1) /= n_s) then
              deallocate(vec_occ, svec_occ, sdc_occ, vec_src, ssrc_vec)
              allocate(vec_occ(n_s), svec_occ(n_s), sdc_occ(n_s), vec_src(n_s), ssrc_vec(n_s))
            end if
          end if
          if (enable_excitation_source_trace) then
            if (.not. allocated(occ_metric_norm)) then
              allocate(occ_metric_norm(nocc_spin))
            else if (size(occ_metric_norm, 1) /= nocc_spin) then
              deallocate(occ_metric_norm)
              allocate(occ_metric_norm(nocc_spin))
            end if
            source_norm(:) = 0.0d0
            source_occ_proj(:) = 0.0d0
            do jo_occ = 1, nocc_spin
              if (n_s > 0) then
                vec_occ(:) = coef_all(1:n_s, jo_occ)
                call apply_overlap_operator(dg_frag, ispin, vec_occ, svec_occ, use_prop_overlap)
                csc_local = real(sum(conjg(vec_occ) * svec_occ), kind=8)
              else
                csc_local = real(sum(conjg(coef_all(1:n_frag, jo_occ)) * coef_all(1:n_frag, jo_occ)), kind=8)
              end if
              call comm_summation(csc_local, occ_metric_norm(jo_occ), dg_frag%icomm)
            end do
          end if
          dndt_before_local = 0.0d0
          dndt_after_local = 0.0d0
          do io_occ = 1, nocc_spin
            if (n_s > 0) then
              vec_occ(:) = coef_all(1:n_s, io_occ)
              call apply_overlap_operator(dg_frag, ispin, vec_occ, svec_occ, use_prop_overlap)
              csc_local = real(sum(conjg(vec_occ) * svec_occ), kind=8)
            else
              csc_local = real(sum(conjg(coef_all(1:n_frag, io_occ)) * coef_all(1:n_frag, io_occ)), kind=8)
            end if
            call comm_summation(csc_local, csc, dg_frag%icomm)
            if (csc <= 1.0d-20) cycle

            occ_factor = 1.0d0
            if (allocated(system%rocc)) then
              if (io_occ <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
                occ_factor = max(0.0d0, system%rocc(io_occ, 1, ispin))
              end if
            end if

            if (n_s > 0) then
              call apply_overlap_operator(dg_frag, ispin, rhs_all(1:n_s, io_occ), sdc_occ, use_prop_overlap)
              c_s_dc_re_local = real(sum(conjg(vec_occ) * sdc_occ), kind=8)
            else
              c_s_dc_re_local = real(sum(conjg(coef_all(1:n_frag, io_occ)) * rhs_all(1:n_frag, io_occ)), kind=8)
            end if
            call comm_summation(c_s_dc_re_local, c_s_dc_re, dg_frag%icomm)
            dndt_before_local = dndt_before_local + 2.0d0 * occ_factor * c_s_dc_re_local

            if (enable_excitation_source_trace) then
              do source_idx = 1, 6
                select case (source_idx)
                case (1)
                  if (n_s > 0) then
                    vec_src(:) = dcoef_dt_h0_core(1:n_s, io_occ)
                  end if
                case (2)
                  if (n_s > 0) then
                    vec_src(:) = dcoef_dt_nl(1:n_s, io_occ)
                  end if
                case (3)
                  if (n_s > 0) then
                    vec_src(:) = dcoef_dt_spinmix(1:n_s, io_occ)
                  end if
                case (4)
                  if (n_s > 0) then
                    vec_src(:) = dcoef_dt_m(1:n_s, io_occ)
                  end if
                case (5)
                  if (n_s > 0) then
                    vec_src(:) = rhs_pre_solve(1:n_s, io_occ)
                  end if
                case (6)
                  if (n_s > 0) then
                    vec_src(:) = rhs_all(1:n_s, io_occ)
                  end if
                end select

                if (n_s > 0) then
                  call apply_overlap_operator(dg_frag, ispin, vec_src, ssrc_vec, use_prop_overlap)
                  csc_local = real(sum(conjg(vec_src) * ssrc_vec), kind=8)
                else
                  select case (source_idx)
                  case (1)
                    csc_local = real(sum(conjg(dcoef_dt_h0_core(1:n_frag, io_occ)) * dcoef_dt_h0_core(1:n_frag, io_occ)), kind=8)
                  case (2)
                    csc_local = real(sum(conjg(dcoef_dt_nl(1:n_frag, io_occ)) * dcoef_dt_nl(1:n_frag, io_occ)), kind=8)
                  case (3)
                    csc_local = real(sum(conjg(dcoef_dt_spinmix(1:n_frag, io_occ)) * dcoef_dt_spinmix(1:n_frag, io_occ)), kind=8)
                  case (4)
                    csc_local = real(sum(conjg(dcoef_dt_m(1:n_frag, io_occ)) * dcoef_dt_m(1:n_frag, io_occ)), kind=8)
                  case (5)
                    csc_local = real(sum(conjg(rhs_pre_solve(1:n_frag, io_occ)) * rhs_pre_solve(1:n_frag, io_occ)), kind=8)
                  case default
                    csc_local = real(sum(conjg(rhs_all(1:n_frag, io_occ)) * rhs_all(1:n_frag, io_occ)), kind=8)
                  end select
                end if
                call comm_summation(csc_local, source_norm_global, dg_frag%icomm)
                source_norm(source_idx) = source_norm(source_idx) + occ_factor * max(0.0d0, source_norm_global)

                do jo_occ = 1, nocc_spin
                  if (occ_metric_norm(jo_occ) <= 1.0d-20) cycle
                  if (n_s > 0) then
                    occ_overlap_local = sum(conjg(coef_all(1:n_s, jo_occ)) * ssrc_vec)
                  else
                    select case (source_idx)
                    case (1)
                      occ_overlap_local = sum(conjg(coef_all(1:n_frag, jo_occ)) * dcoef_dt_h0_core(1:n_frag, io_occ))
                    case (2)
                      occ_overlap_local = sum(conjg(coef_all(1:n_frag, jo_occ)) * dcoef_dt_nl(1:n_frag, io_occ))
                    case (3)
                      occ_overlap_local = sum(conjg(coef_all(1:n_frag, jo_occ)) * dcoef_dt_spinmix(1:n_frag, io_occ))
                    case (4)
                      occ_overlap_local = sum(conjg(coef_all(1:n_frag, jo_occ)) * dcoef_dt_m(1:n_frag, io_occ))
                    case (5)
                      occ_overlap_local = sum(conjg(coef_all(1:n_frag, jo_occ)) * rhs_pre_solve(1:n_frag, io_occ))
                    case default
                      occ_overlap_local = sum(conjg(coef_all(1:n_frag, jo_occ)) * rhs_all(1:n_frag, io_occ))
                    end select
                  end if
                  call comm_summation(occ_overlap_local, occ_overlap, dg_frag%icomm)
                  source_occ_proj(source_idx) = source_occ_proj(source_idx) + occ_factor * (abs(occ_overlap)**2 / occ_metric_norm(jo_occ))
                end do
              end do
            end if

            if (enforce_norm_tangent) then
              alpha_real = c_s_dc_re / csc
              if (n_s > 0) then
                rhs_all(1:n_s, io_occ) = rhs_all(1:n_s, io_occ) - alpha_real * vec_occ
                call apply_overlap_operator(dg_frag, ispin, rhs_all(1:n_s, io_occ), sdc_occ, use_prop_overlap)
                c_s_dc_re_local = real(sum(conjg(vec_occ) * sdc_occ), kind=8)
              else
                rhs_all(1:n_frag, io_occ) = rhs_all(1:n_frag, io_occ) - alpha_real * coef_all(1:n_frag, io_occ)
                c_s_dc_re_local = real(sum(conjg(coef_all(1:n_frag, io_occ)) * rhs_all(1:n_frag, io_occ)), kind=8)
              end if
              call comm_summation(c_s_dc_re_local, c_s_dc_re, dg_frag%icomm)
            end if
            dndt_after_local = dndt_after_local + 2.0d0 * occ_factor * c_s_dc_re_local
          end do

          call comm_summation(dndt_before_local, dndt_before, dg_frag%icomm)
          call comm_summation(dndt_after_local, dndt_after, dg_frag%icomm)

          if (enable_norm_deriv_check .and. dg_frag%id == 0 .and. itt <= 120) then
            write(*,'(1x,a,i0,a,i0,a,2(1x,1pe12.4))') '        norm-deriv-check: itt=', itt, ' ispin=', ispin, &
              ' dNdt(before/after)=', dndt_before, dndt_after
            flush(6)
          end if
          if (enable_excitation_source_trace .and. dg_frag%id == 0 .and. itt <= 2) then
            source_leak(:) = max(0.0d0, source_norm(:) - source_occ_proj(:))
            write(*,'(1x,a,i0,a,i0,6(a,1pe12.4))') '        excitation-source-total: itt=', itt, ' ispin=', ispin, &
              ' h0_core=', source_norm(1), ' nl=', source_norm(2), ' spinmix=', source_norm(3), &
              ' m=', source_norm(4), ' rhs_preS=', source_norm(5), ' rhs_postS=', source_norm(6)
            write(*,'(1x,a,i0,a,i0,6(a,1pe12.4))') '        excitation-source-leak: itt=', itt, ' ispin=', ispin, &
              ' h0_core=', source_leak(1), ' nl=', source_leak(2), ' spinmix=', source_leak(3), &
              ' m=', source_leak(4), ' rhs_preS=', source_leak(5), ' rhs_postS=', source_leak(6)
            flush(6)
          end if
        end if
      end if
      if (enable_deriv_trace .and. itt <= 2) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        deriv-trace: rank=", dg_frag%id, " itt=", itt, &
          " ispin=", ispin, " stage=after-overlap-solve"
        flush(6)
      end if

      if (trace_deriv_step) then
        t_trace0 = get_wtime()
        write(*,'(1x,a,i0,a,i0)') '[DG-DERIV] copy-out start call=', derivative_step_trace_id, &
          ' ispin=', ispin
        flush(6)
      end if
      if (use_mixed_basis .and. n_basis > 0) then
        ! Mixed-basis RK stages keep a replicated canonical raw coefficient
        ! view on every orbital rank.  The derivative must therefore be
        ! replicated too; owner-only writes leave stale non-owner k rows that
        ! the RK update still reads.
        dcoef_dt(1:n_frag, 1:dg_frag%nstate_tot, ispin) = rhs_all(1:n_frag, :)
        if (present(dcoef_dt_pw) .and. n_pw > 0) then
          dcoef_dt_pw(1:n_pw, 1:dg_frag%nstate_tot, ispin) = rhs_all(n_frag+1:n_tot, :)
        end if
      else
        do io = 1, n_frag
          if (.not. allocated(dg_frag%coef_owner)) exit
          if (dg_frag%coef_owner(io, ispin) /= dg_frag%id) cycle
          dcoef_dt(io, 1:dg_frag%nstate_tot, ispin) = rhs_all(io, :)
        end do
        if (present(dcoef_dt_pw) .and. n_pw > 0) then
          do io = 1, n_pw
            if (.not. allocated(dg_frag%coef_pw_owner)) exit
            if (dg_frag%coef_pw_owner(io) /= dg_frag%id) cycle
            dcoef_dt_pw(io, 1:dg_frag%nstate_tot, ispin) = rhs_all(n_frag+io, :)
          end do
        end if
      end if
      if (trace_deriv_step) then
        t_trace1 = get_wtime()
        write(*,'(1x,a,i0,a,i0,a,1pe12.4)') '[DG-DERIV] copy-out done call=', &
          derivative_step_trace_id, ' ispin=', ispin, ' time=', t_trace1 - t_trace0
        flush(6)
      end if
      if (enable_deriv_trace .and. itt <= 2) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        deriv-trace: rank=", dg_frag%id, " itt=", itt, &
          " ispin=", ispin, " stage=spin-end"
        flush(6)
      end if
    end do

    t_deriv1 = get_wtime()
    dt_deriv = t_deriv1 - t_deriv0
    if (trace_deriv_step) then
      write(*,'(1x,a,i0,a,1pe12.4,a,1pe12.4)') '[DG-DERIV] leave call=', derivative_step_trace_id, &
        ' total=', dt_deriv, ' gather_total=', dt_gather_local
      flush(6)
    end if

    ! Cache retained for reuse

  contains

    logical function distributed_mixed_rows_active() result(active)
      implicit none
      active = dg_frag%parallel_mode_orbital .and. dg_frag%isize > 1 .and. allocated(dg_frag%coef_owner)
    end function distributed_mixed_rows_active

    logical function owns_fragment_row(row_idx) result(owned)
      implicit none
      integer, intent(in) :: row_idx
      owned = .false.
      if (row_idx < 1 .or. row_idx > n_frag) return
      if (.not. allocated(dg_frag%coef_owner)) return
      if (row_idx > size(dg_frag%coef_owner, 1) .or. ispin > size(dg_frag%coef_owner, 2)) return
      owned = (dg_frag%coef_owner(row_idx, ispin) == dg_frag%id)
    end function owns_fragment_row

    logical function owns_pw_row_global(pw_idx) result(owned)
      implicit none
      integer, intent(in) :: pw_idx
      owned = .false.
      if (pw_idx < 1 .or. pw_idx > n_pw) return
      if (.not. allocated(dg_frag%coef_pw_owner)) return
      if (pw_idx > size(dg_frag%coef_pw_owner, 1)) return
      ! The PW block is global in the mixed matrix.  coef_pw_owner is defined
      ! inside each fragment subgroup, so only the first subgroup contributes
      ! PW rows to world reductions to avoid counting the same PW row once per
      ! fragment.
      if (dg_frag%ifrag_group /= 1) return
      owned = (dg_frag%coef_pw_owner(pw_idx) == dg_frag%id)
    end function owns_pw_row_global

    subroutine next_owned_fragment_segment(start_row, row_s, row_e)
      implicit none
      integer, intent(inout) :: start_row
      integer, intent(out) :: row_s, row_e

      row_s = 0
      row_e = -1
      do while (start_row <= n_frag)
        if (owns_fragment_row(start_row)) exit
        start_row = start_row + 1
      end do
      if (start_row > n_frag) return
      row_s = start_row
      row_e = row_s
      do while (row_e < n_frag)
        if (.not. owns_fragment_row(row_e + 1)) exit
        row_e = row_e + 1
      end do
      start_row = row_e + 1
    end subroutine next_owned_fragment_segment

    subroutine next_owned_pw_segment(start_pw, pw_s, pw_e)
      implicit none
      integer, intent(inout) :: start_pw
      integer, intent(out) :: pw_s, pw_e

      pw_s = 0
      pw_e = -1
      do while (start_pw <= n_pw)
        if (owns_pw_row_global(start_pw)) exit
        start_pw = start_pw + 1
      end do
      if (start_pw > n_pw) return
      pw_s = start_pw
      pw_e = pw_s
      do while (pw_e < n_pw)
        if (.not. owns_pw_row_global(pw_e + 1)) exit
        pw_e = pw_e + 1
      end do
      start_pw = pw_e + 1
    end subroutine next_owned_pw_segment

    subroutine accumulate_projected_rows(raw_op, projected, row_s, row_e)
      implicit none
      complex(8), intent(in) :: raw_op(n_tot, n_tot)
      complex(8), intent(inout) :: projected(n_basis, n_basis)
      integer, intent(in) :: row_s, row_e
      integer :: nrow

      if (row_s < 1 .or. row_e < row_s) return
      nrow = row_e - row_s + 1
      call zgemm('N', 'N', nrow, n_basis, n_tot, (1.0d0, 0.0d0), &
                 raw_op(row_s, 1), n_tot, dg_frag%mixed_transform(1, 1, ispin), n_tot, &
                 (0.0d0, 0.0d0), tmp_mix(row_s, 1), n_tot)
      call zgemm('C', 'N', n_basis, n_basis, nrow, (1.0d0, 0.0d0), &
                 dg_frag%mixed_transform(row_s, 1, ispin), n_tot, tmp_mix(row_s, 1), n_tot, &
                 (1.0d0, 0.0d0), projected, n_basis)
    end subroutine accumulate_projected_rows

    subroutine build_mixed_projected_operator(raw_op, projected)
      implicit none
      complex(8), intent(in) :: raw_op(n_tot, n_tot)
      complex(8), intent(inout) :: projected(n_basis, n_basis)
      integer :: row_cursor, row_s, row_e
      integer :: pw_cursor, pw_s, pw_e

      if (distributed_mixed_rows_active()) then
        projected(:, :) = (0.0d0, 0.0d0)

        row_cursor = 1
        do
          call next_owned_fragment_segment(row_cursor, row_s, row_e)
          if (row_s <= 0) exit
          call accumulate_projected_rows(raw_op, projected, row_s, row_e)
        end do

        if (n_pw > 0) then
          pw_cursor = 1
          do
            call next_owned_pw_segment(pw_cursor, pw_s, pw_e)
            if (pw_s <= 0) exit
            call accumulate_projected_rows(raw_op, projected, n_frag + pw_s, n_frag + pw_e)
          end do
        end if

        call comm_summation(projected, op_mix_sum, n_basis * n_basis, dg_frag%icomm)
        projected(:, :) = op_mix_sum(:, :)
      else
        call zgemm('N', 'N', n_tot, n_basis, n_tot, (1.0d0, 0.0d0), raw_op, n_tot, &
                   dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin), n_tot, (0.0d0, 0.0d0), tmp_mix, n_tot)
        call zgemm('C', 'N', n_basis, n_basis, n_tot, (1.0d0, 0.0d0), &
                   dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin), n_tot, tmp_mix, n_tot, &
                   (0.0d0, 0.0d0), projected, n_basis)
      end if
    end subroutine build_mixed_projected_operator

    subroutine expand_raw_rows(rhs, raw_out, row_s, row_e)
      implicit none
      complex(8), intent(in) :: rhs(n_basis, dg_frag%nstate_tot)
      complex(8), intent(inout) :: raw_out(n_tot, dg_frag%nstate_tot)
      integer, intent(in) :: row_s, row_e
      integer :: nrow

      if (row_s < 1 .or. row_e < row_s) return
      nrow = row_e - row_s + 1
      call zgemm('N', 'N', nrow, dg_frag%nstate_tot, n_basis, (1.0d0, 0.0d0), &
                 dg_frag%mixed_transform(row_s, 1, ispin), n_tot, rhs, n_basis, &
                 (0.0d0, 0.0d0), raw_out(row_s, 1), size(raw_out, 1))
    end subroutine expand_raw_rows

    subroutine expand_mixed_rhs_to_raw(rhs, raw_out)
      implicit none
      complex(8), intent(in) :: rhs(n_basis, dg_frag%nstate_tot)
      complex(8), intent(inout) :: raw_out(n_tot, dg_frag%nstate_tot)
      integer :: row_cursor, row_s, row_e
      integer :: pw_cursor, pw_s, pw_e

      if (distributed_mixed_rows_active()) then
        raw_out(:, :) = (0.0d0, 0.0d0)

        row_cursor = 1
        do
          call next_owned_fragment_segment(row_cursor, row_s, row_e)
          if (row_s <= 0) exit
          call expand_raw_rows(rhs, raw_out, row_s, row_e)
        end do

        if (n_pw > 0) then
          pw_cursor = 1
          do
            call next_owned_pw_segment(pw_cursor, pw_s, pw_e)
            if (pw_s <= 0) exit
            call expand_raw_rows(rhs, raw_out, n_frag + pw_s, n_frag + pw_e)
          end do
        end if

        call comm_summation(raw_out, raw_rhs_sum, n_tot * dg_frag%nstate_tot, dg_frag%icomm)
        raw_out(:, :) = raw_rhs_sum(:, :)
      else
        call zgemm('N', 'N', n_tot, dg_frag%nstate_tot, n_basis, (1.0d0, 0.0d0), &
                   dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin), n_tot, rhs, n_basis, &
                   (0.0d0, 0.0d0), raw_out, n_tot)
      end if
    end subroutine expand_mixed_rhs_to_raw

    subroutine trace_deriv_norm(label, arr)
      implicit none
      character(len=*), intent(in) :: label
      complex(8), intent(in) :: arr(:, :)
      integer :: nr, nc
      real(8) :: coef_norm, arr_norm, ratio

      if (.not. trace_deriv_step) return
      nr = min(n_tot, size(arr, 1), size(coef_all, 1))
      nc = min(dg_frag%nstate_tot, size(arr, 2), size(coef_all, 2))
      if (nr <= 0 .or. nc <= 0) return
      coef_norm = sqrt(sum(abs(coef_all(1:nr, 1:nc))**2))
      arr_norm = sqrt(sum(abs(arr(1:nr, 1:nc))**2))
      ratio = 0.0d0
      if (coef_norm > 1.0d-300) ratio = arr_norm / coef_norm
      write(*,'(1x,a,i0,a,i0,a,a,a,3(1x,1pe14.6))') '[DG-DERIV-NORM] call=', derivative_step_trace_id, &
        ' ispin=', ispin, ' stage=', trim(label), ' norm ratio coef_norm=', arr_norm, ratio, coef_norm
    end subroutine trace_deriv_norm

    subroutine trace_h_blocks_summary()
      implicit none
      integer :: iblk, nrow, ncol
      real(8) :: hmax_self, hmax_off, hfrob_self, hfrob_off

      if (.not. trace_deriv_step) return
      if (.not. allocated(dg_frag%H_mat_blocks)) return
      hmax_self = 0.0d0
      hmax_off = 0.0d0
      hfrob_self = 0.0d0
      hfrob_off = 0.0d0
      do iblk = 1, size(dg_frag%H_mat_blocks)
        if (dg_frag%H_mat_blocks(iblk)%ifrag_row < 1 .or. dg_frag%H_mat_blocks(iblk)%ifrag_row > dg_frag%n_frag) cycle
        if (dg_frag%H_mat_blocks(iblk)%ifrag_col < 1 .or. dg_frag%H_mat_blocks(iblk)%ifrag_col > dg_frag%n_frag) cycle
        nrow = min(dg_frag%n_basis(dg_frag%H_mat_blocks(iblk)%ifrag_row, ispin), &
                   size(dg_frag%H_mat_blocks(iblk)%val, 1))
        ncol = min(dg_frag%n_basis(dg_frag%H_mat_blocks(iblk)%ifrag_col, ispin), &
                   size(dg_frag%H_mat_blocks(iblk)%val, 2))
        if (nrow <= 0 .or. ncol <= 0) cycle
        if (dg_frag%H_mat_blocks(iblk)%ifrag_row == dg_frag%H_mat_blocks(iblk)%ifrag_col) then
          hmax_self = max(hmax_self, maxval(abs(dg_frag%H_mat_blocks(iblk)%val(1:nrow, 1:ncol, ispin))))
          hfrob_self = hfrob_self + sum(dg_frag%H_mat_blocks(iblk)%val(1:nrow, 1:ncol, ispin)**2)
        else
          hmax_off = max(hmax_off, maxval(abs(dg_frag%H_mat_blocks(iblk)%val(1:nrow, 1:ncol, ispin))))
          hfrob_off = hfrob_off + sum(dg_frag%H_mat_blocks(iblk)%val(1:nrow, 1:ncol, ispin)**2)
        end if
      end do
      write(*,'(1x,a,i0,a,i0,a,4(1x,1pe14.6))') '[DG-H-BLOCK-NORM] call=', derivative_step_trace_id, &
        ' ispin=', ispin, ' max_self max_off frob_self frob_off=', hmax_self, hmax_off, sqrt(hfrob_self), sqrt(hfrob_off)
    end subroutine trace_h_blocks_summary

    subroutine apply_block_nonlocal_cache(nl_out, zero_output)
      implicit none
      complex(8), intent(inout) :: nl_out(n_tot, dg_frag%nstate_tot)
      logical, intent(in), optional :: zero_output
      logical :: do_zero

      do_zero = .true.
      if (present(zero_output)) do_zero = zero_output
      if (do_zero) nl_out(:, :) = (0.0d0, 0.0d0)
      if (.not. allocated(dg_frag%H_nl_blocks)) return
      if (allocated(dg_frag%H_nl_local_block_ids)) then
        call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, &
             coef_all(1:n_frag, :), nl_out(1:n_frag, :), dg_frag%H_nl_local_block_ids)
      else
        call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, &
             coef_all(1:n_frag, :), nl_out(1:n_frag, :))
      end if
    end subroutine apply_block_nonlocal_cache

  end subroutine calculate_time_derivative
