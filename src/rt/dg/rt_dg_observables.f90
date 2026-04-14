  subroutine calculate_observables(dg_frag, system, mg, stencil, ppg, rt, itt, Vh, Vxc, Vpsl)
    use salmon_global, only: theory
    use structures
    use communication, only: comm_summation
    use timer, only: timer_begin, timer_end, LOG_CALC_CURRENT
    use sym_vector_sub, only: sym_vector_xyz
    use rt_dg_fragment_ops, only: apply_momentum_blocks, apply_matrix_blocks_batch, apply_nonlocal_pp_projector_batch, &
                    apply_mixed_hamiltonian, mixed_fp_coupling_active, copy_matrix_blocks_to_complex_dense, &
                    gather_full_coef_view, apply_overlap_operator, copy_overlap_operator_to_dense, &
                    copy_momentum_blocks_to_complex_dense, solve_overlap_operator_batch
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
    
    integer :: io, jo, ispin, idir, n, nocc, n_pw, n_tot, max_nocc, n_mix
    integer :: active_state_cap
    integer :: ifrag, jfrag, ib, jb, i_idx, j_idx
    integer :: iblk, nrow_blk, ncol_blk, n_diag_block_ids, idb
    integer :: env_len, env_status, probe_stride, probe_nprint
    integer :: occ_trace_n
    integer :: current_norm_stride, current_component_stride, excitation_trace_stride, para_decomp_trace_stride
    integer :: transition_stride
    character(len=64) :: env_val
    logical :: do_interface_check
    logical :: enable_orbital_probe
    logical :: enable_transition_probe
    logical :: enable_electron_probe
    logical :: enable_obs_trace
    logical :: enable_occ_trace
    logical :: enable_realspace_probe
    logical :: enable_current_norm_trace
    logical :: use_s_metric_current
    logical :: use_s_adjoint_current
    logical :: enable_current_component_trace
    logical :: enable_excitation_trace, do_excitation_sample
    logical :: enable_para_decomp_trace, do_para_decomp_sample
    real(8), allocatable :: interface_flow(:,:), dndt_frag(:)
    real(8) :: pair_residual, max_pair_residual, charge_balance_residual
    real(8) :: current_tmp, energy_tmp, pw_weight_local, kpw_dir
    real(8) :: current_io, energy_io, energy_mix_io
    real(8) :: occ_i
    real(8) :: charge_trace_spin
    real(8) :: Ac_tot(3), A_squared
    real(8) :: current_local(3), current_para_local(3), current_dia_local(3), energy_local, energy_mix_local
    real(8) :: current_reduce_sum(3), current_after_isize(3)
    real(8) :: current_para_reduce_sum(3), current_para_after_isize(3)
    real(8) :: current_dia_reduce_sum(3), current_dia_after_isize(3)
    real(8) :: current_nl_reduce_sum(3), current_nl_after_isize(3)
    real(8) :: elec_coef_local, elec_plain_local
    real(8) :: elec_coef_sum, elec_plain_sum
    real(8) :: orth_def_fro_local, orth_offdiag_max_local, orth_diag_min_local, orth_diag_max_local
    real(8) :: orth_def_fro_sum, orth_offdiag_max_sum, orth_diag_min_sum, orth_diag_max_sum
    real(8) :: energy_static_local, energy_kin_local, energy_nl_local, energy_ap_local, energy_a2_local
    real(8) :: energy_static_sum, energy_kin_sum, energy_nl_sum, energy_ap_sum, energy_a2_sum
    real(8) :: energy_kin_rs_sum, energy_one_rs_sum, energy_nl_rs_sum
    real(8) :: energy_kin_diag_local, energy_kin_offdiag_local
    real(8) :: energy_kin_diag_sum, energy_kin_offdiag_sum
    real(8) :: kinetic_diag_abs_local, kinetic_offdiag_abs_local
    real(8) :: kinetic_diag_abs_sum, kinetic_offdiag_abs_sum
    real(8) :: kinetic_apply_diff_local, kinetic_apply_diff_sum
    real(8) :: energy_static_avg, energy_kin_avg, energy_nl_avg, energy_ap_avg, energy_a2_avg
    real(8) :: frag_reduce_factor
    real(8) :: dm_diag_norm_local, dm_offdiag_norm_local, dm_offdiag_max_local
    real(8) :: dm_diag_norm_sum, dm_offdiag_norm_sum, dm_offdiag_max_sum
    real(8) :: dm_offdiag_ratio, dm_diag_delta, dm_offdiag_delta, dm_offdiag_ratio_delta
    real(8), save :: dm_diag_norm_ref = -1.0d0, dm_offdiag_norm_ref = -1.0d0, dm_offdiag_ratio_ref = -1.0d0
    logical, save :: dm_excitation_ref_initialized = .false.
    real(8) :: para_cdiag_local(3), para_coff_local(3), para_pdiag_local(3), para_poff_local(3)
    real(8) :: para_cdiag_sum(3), para_coff_sum(3), para_pdiag_sum(3), para_poff_sum(3)
    real(8) :: para_cdiag_now(3), para_coff_now(3), para_pdiag_now(3), para_poff_now(3)
    real(8) :: para_coff_abs_local(3), para_poff_abs_local(3)
    real(8) :: para_coff_abs_sum(3), para_poff_abs_sum(3)
    real(8) :: para_coff_abs_now(3), para_poff_abs_now(3)
    complex(8) :: minus_i
    complex(8), allocatable :: op_mat(:,:), tmp_mat(:,:), coef_all(:,:), tmp_all(:,:)
    complex(8), allocatable :: coef_frag_all(:,:), coef_pw_all(:,:), coef_frag_view(:,:), coef_pw_view(:,:)
    complex(8), allocatable :: tmp_probe(:,:), dense_probe_mat(:,:), dense_probe_out(:,:)
    complex(8), allocatable :: overlap_vec(:), overlap_dense(:,:)
    complex(8), allocatable :: overlap_rhs(:,:), tmp_adj(:,:), tmp_solve(:,:)
    logical :: has_nonlocal, use_hmat_complex, use_mixed_current
    logical :: require_dense_nl
    logical :: use_spatial_A
    logical, parameter :: enable_energy_component_probe = .false.
    logical :: use_energy_components
    real(8), allocatable :: Ap_mat(:,:), A2_mat(:,:)
    integer, allocatable :: diag_block_ids(:)
    real(8), allocatable :: occ_weight(:)
    real(8), allocatable :: current_orb_local(:), current_orb_sum(:)
    real(8), allocatable :: energy_orb_local(:), energy_orb_sum(:)
    real(8), allocatable :: nl_state_local(:), nl_state_sum(:)
    real(8) :: current_diag_local(3), current_offdiag_local(3)
    real(8) :: current_diag_sum(3), current_offdiag_sum(3)
    real(8) :: current_nl_rs_sum(3)
    complex(8) :: current_blk_total, current_blk_diag, current_elem
    complex(8) :: mfp
    real(8), parameter :: unit_dir(3,3) = reshape((/ &
      1.0d0, 0.0d0, 0.0d0, &
      0.0d0, 1.0d0, 0.0d0, &
      0.0d0, 0.0d0, 1.0d0 /), (/3, 3/))

    Ac_tot = rt%Ac_tot(:, itt)

    ! Calculate local observables (only for assigned fragments)
    ! MPI aggregation will sum across all ranks
    current_local = 0.0d0
    current_para_local = 0.0d0
    current_dia_local = 0.0d0
    energy_local = 0.0d0
    energy_mix_local = 0.0d0
    pw_weight_local = 0.0d0
    elec_coef_local = 0.0d0
    elec_plain_local = 0.0d0
    orth_def_fro_local = 0.0d0
    orth_offdiag_max_local = 0.0d0
    orth_diag_min_local = huge(1.0d0)
    orth_diag_max_local = -huge(1.0d0)
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
    energy_nl_rs_sum = 0.0d0
    current_diag_local(:) = 0.0d0
    current_offdiag_local(:) = 0.0d0
    current_nl_rs_sum(:) = 0.0d0
    n = dg_frag%n_mat_max
    use_spatial_A = (trim(theory) == 'single_scale_maxwell_tddft' .and. allocated(system%Ac_micro%v) .and. dg_frag%has_real_space_basis)

    do_interface_check = .false.
    if (do_interface_check) then
      allocate(interface_flow(dg_frag%n_frag, dg_frag%n_frag))
      interface_flow = 0.0d0
    end if
    if (n <= 0) then
      current_local = 0.0d0
      energy_local = 0.0d0
      energy_mix_local = 0.0d0
      goto 1000
    end if
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n + n_pw
    use_energy_components = enable_energy_component_probe
    active_state_cap = max(1, min(dg_frag%nstate_tot, n))
    if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
      max_nocc = active_state_cap
    else
      max_nocc = max(1, maxval(dg_frag%nocc_spin(1:dg_frag%nspin)))
    end if

    enable_orbital_probe = .false.
    probe_stride = 1
    probe_nprint = 20
    enable_transition_probe = .false.
    enable_electron_probe = .false.
    enable_obs_trace = .false.
    enable_occ_trace = .false.
    ! Pure-DG (no PW rows) needs the real-space nonlocal current term for
    ! consistency with the conventional current decomposition.
    enable_realspace_probe = (n_pw == 0)
    enable_current_norm_trace = .false.
    enable_current_component_trace = .false.
    enable_excitation_trace = .false.
    enable_para_decomp_trace = .false.
    use_s_metric_current = .false.
    use_s_adjoint_current = .false.
    transition_stride = 1
    occ_trace_n = 12
    current_norm_stride = 1
    current_component_stride = 1
    excitation_trace_stride = 1
    para_decomp_trace_stride = 1
    call get_environment_variable("SALMON_DG_ORBITAL_PROBE", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        enable_orbital_probe = .true.
      end if
    end if
    if (enable_orbital_probe) then
      call get_environment_variable("SALMON_DG_ORBITAL_PROBE_STRIDE", env_val, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_val(1:env_len), *, iostat=env_status) probe_stride
        if (env_status /= 0 .or. probe_stride < 1) probe_stride = 1
      end if
      call get_environment_variable("SALMON_DG_ORBITAL_PROBE_N", env_val, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_val(1:env_len), *, iostat=env_status) probe_nprint
        if (env_status /= 0 .or. probe_nprint < 1) probe_nprint = 20
      end if
    end if
    call get_environment_variable("SALMON_DG_TRANSITION_PROBE", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        enable_transition_probe = .true.
      end if
    end if
    call get_environment_variable("SALMON_DG_ELECTRON_PROBE", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        enable_electron_probe = .true.
      end if
    end if
    call get_environment_variable("SALMON_DG_OBS_TRACE", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        enable_obs_trace = .true.
      end if
    end if
    call get_environment_variable("SALMON_DG_REALSPACE_PROBE", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        enable_realspace_probe = .true.
      end if
    end if
    call get_environment_variable("SALMON_DG_OCC_TRACE", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        enable_occ_trace = .true.
      end if
    end if
    if (enable_occ_trace) then
      call get_environment_variable("SALMON_DG_OCC_TRACE_N", env_val, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_val(1:env_len), *, iostat=env_status) occ_trace_n
        if (env_status /= 0 .or. occ_trace_n < 1) occ_trace_n = 12
      end if
    end if
    call get_environment_variable("SALMON_DG_CURRENT_NORM_TRACE", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        enable_current_norm_trace = .true.
      end if
    end if
    if (enable_current_norm_trace) then
      call get_environment_variable("SALMON_DG_CURRENT_NORM_TRACE_STRIDE", env_val, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_val(1:env_len), *, iostat=env_status) current_norm_stride
        if (env_status /= 0 .or. current_norm_stride < 1) current_norm_stride = 1
      end if
    end if
    call get_environment_variable("SALMON_DG_CURRENT_COMPONENT_TRACE", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        enable_current_component_trace = .true.
      end if
    end if
    if (enable_current_component_trace) then
      call get_environment_variable("SALMON_DG_CURRENT_COMPONENT_TRACE_STRIDE", env_val, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_val(1:env_len), *, iostat=env_status) current_component_stride
        if (env_status /= 0 .or. current_component_stride < 1) current_component_stride = 1
      end if
    end if
    call get_environment_variable("SALMON_DG_CURRENT_USE_S_METRIC", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        use_s_metric_current = .true.
      end if
    end if
    call get_environment_variable("SALMON_DG_CURRENT_USE_S_ADJOINT", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        use_s_adjoint_current = .true.
      end if
    end if
    call get_environment_variable("SALMON_DG_EXCITATION_TRACE", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        enable_excitation_trace = .true.
      end if
    end if
    if (enable_excitation_trace) then
      call get_environment_variable("SALMON_DG_EXCITATION_TRACE_STRIDE", env_val, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_val(1:env_len), *, iostat=env_status) excitation_trace_stride
        if (env_status /= 0 .or. excitation_trace_stride < 1) excitation_trace_stride = 1
      end if
    end if
    call get_environment_variable("SALMON_DG_PARA_DECOMP_TRACE", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        enable_para_decomp_trace = .true.
      end if
    end if
    if (enable_para_decomp_trace) then
      call get_environment_variable("SALMON_DG_PARA_DECOMP_TRACE_STRIDE", env_val, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_val(1:env_len), *, iostat=env_status) para_decomp_trace_stride
        if (env_status /= 0 .or. para_decomp_trace_stride < 1) para_decomp_trace_stride = 1
      end if
    end if
    if (enable_transition_probe) then
      call get_environment_variable("SALMON_DG_TRANSITION_PROBE_STRIDE", env_val, length=env_len, status=env_status)
      if (env_status == 0 .and. env_len > 0) then
        read(env_val(1:env_len), *, iostat=env_status) transition_stride
        if (env_status /= 0 .or. transition_stride < 1) transition_stride = 1
      end if
    end if

    allocate(tmp_mat(n, max_nocc))
    if (use_energy_components) allocate(tmp_probe(n, max_nocc))
    allocate(coef_frag_all(n, max_nocc))
    allocate(occ_weight(max_nocc))
    if (enable_orbital_probe) then
      allocate(current_orb_local(3 * max_nocc), current_orb_sum(3 * max_nocc))
      allocate(energy_orb_local(max_nocc), energy_orb_sum(max_nocc))
      current_orb_local(:) = 0.0d0
      energy_orb_local(:) = 0.0d0
    end if
    if (n_pw > 0) then
      allocate(coef_pw_all(n_pw, max_nocc))
      allocate(coef_all(n_tot, max_nocc), tmp_all(n_tot, max_nocc))
    end if
    if (enable_electron_probe .or. use_s_metric_current) then
      allocate(overlap_vec(n_tot))
      if (n_pw == 0) allocate(overlap_dense(n, n))
    end if
    if (use_s_adjoint_current .and. n_pw == 0) then
      if (.not. allocated(overlap_vec)) allocate(overlap_vec(n))
      allocate(overlap_rhs(n, max_nocc), tmp_adj(n, max_nocc), tmp_solve(n, max_nocc))
    end if
    minus_i = cmplx(0.0d0, -1.0d0, kind=8)
    dm_diag_norm_local = 0.0d0
    dm_offdiag_norm_local = 0.0d0
    dm_offdiag_max_local = 0.0d0
    para_cdiag_local(:) = 0.0d0
    para_coff_local(:) = 0.0d0
    para_pdiag_local(:) = 0.0d0
    para_poff_local(:) = 0.0d0
    para_coff_abs_local(:) = 0.0d0
    para_poff_abs_local(:) = 0.0d0

    if (enable_obs_trace) then
      write(*,'(1x,a,i0,a,i0,a,i0,a)') "        obs-trace: rank=", dg_frag%id, " itt=", itt, " n=", n, " stage=setup-done"
      flush(6)
    end if

    ! Current calculation via momentum operator matrix (velocity gauge)
    ! Following conventional RT implementation in density_matrix.f90:
    !   - momentum_mat stores <φ|∇|φ> (gradient operator)
    !   - Current: j = Im[<ψ|∇|ψ>] with factor 2 and final normalization by cell volume
    !   - Sign: Testing -2.0 to match conventional RT direction
    call timer_begin(LOG_CALC_CURRENT)
    do ispin = 1, dg_frag%nspin
      if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
        nocc = active_state_cap
      else
        nocc = min(dg_frag%nocc_spin(ispin), active_state_cap)
      end if
      if (nocc <= 0) cycle
      occ_weight(:) = 0.0d0
      do io = 1, nocc
        if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
          if (io <= size(dg_frag%occ_state, 1) .and. ispin <= size(dg_frag%occ_state, 2)) then
            occ_weight(io) = max(0.0d0, dg_frag%occ_state(io, ispin))
          end if
        else
          occ_weight(io) = system%rocc(io, 1, ispin)
        end if
      end do
      charge_trace_spin = sum(occ_weight(1:nocc))
      if (enable_occ_trace .and. dg_frag%id == 0 .and. (itt == 1 .or. mod(itt, 20) == 0)) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
          "        occ-trace: itt=", itt, " spin=", ispin, " nocc=", nocc, &
          " sum=", charge_trace_spin, " min=", minval(occ_weight(1:nocc)), " max=", maxval(occ_weight(1:nocc))
        write(*,'(1x,a,20(1x,1pe12.4))') "        occ-trace first:", &
          occ_weight(1:min(nocc, min(occ_trace_n, 20)))
      end if
      use_mixed_current = (n_pw > 0 .and. mixed_fp_coupling_active(dg_frag, ispin))
      call gather_full_coef_view(dg_frag, ispin, n, nocc, coef_frag_view, coef_pw_view, 1, nocc)
      coef_frag_all(1:n, 1:nocc) = coef_frag_view(1:n, 1:nocc)
      do_excitation_sample = enable_excitation_trace .and. (n_pw == 0) .and. &
                (itt == 0 .or. itt == 1 .or. mod(itt - 1, excitation_trace_stride) == 0)
      do_para_decomp_sample = enable_para_decomp_trace .and. (n_pw == 0) .and. (.not. use_mixed_current) .and. &
             (itt == 0 .or. itt == 1 .or. mod(itt - 1, para_decomp_trace_stride) == 0)
      if (do_excitation_sample) then
        do ib = 1, n
          do jb = 1, n
            current_elem = (0.0d0, 0.0d0)
            do io = 1, nocc
              if (occ_weight(io) <= 0.0d0) cycle
              current_elem = current_elem + occ_weight(io) * coef_frag_all(ib, io) * conjg(coef_frag_all(jb, io))
            end do
            if (ib == jb) then
              dm_diag_norm_local = dm_diag_norm_local + abs(current_elem)
            else
              dm_offdiag_norm_local = dm_offdiag_norm_local + abs(current_elem)
              dm_offdiag_max_local = max(dm_offdiag_max_local, abs(current_elem))
            end if
          end do
        end do
      end if
      if (n_pw > 0) coef_pw_all(1:n_pw, 1:nocc) = coef_pw_view(1:n_pw, 1:nocc)
      if (n_pw > 0) then
        coef_all(1:n_tot, 1:nocc) = (0.0d0, 0.0d0)
        coef_all(1:n, 1:nocc) = coef_frag_all(1:n, 1:nocc)
        coef_all(n+1:n_tot, 1:nocc) = coef_pw_all(1:n_pw, 1:nocc)
      end if
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
                       coef_all, n_tot, (0.0d0, 0.0d0), tmp_all, n_tot)
          else
            if (.not. allocated(op_mat)) allocate(op_mat(n, n))
            op_mat(:, :) = cmplx(dg_frag%momentum_mat(idir, 1:n, 1:n, ispin), 0.0d0, kind=8)
            call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
                       coef_all, n_tot, (0.0d0, 0.0d0), tmp_all, n_tot)
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
          do io = 1, nocc
            if (occ_weight(io) <= 0.0d0) cycle
            current_io = 0.0d0
            if (use_s_metric_current) then
              call apply_overlap_operator(dg_frag, ispin, tmp_all(1:n_tot, io), overlap_vec(1:n_tot), .true.)
              do ib = 1, n_tot
                current_io = current_io + aimag(conjg(coef_all(ib, io)) * overlap_vec(ib))
              end do
            else
              do ib = 1, n_tot
                current_io = current_io + aimag(conjg(coef_all(ib, io)) * tmp_all(ib, io))
              end do
            end if
            current_tmp = current_tmp + occ_weight(io) * current_io
            if (enable_orbital_probe) current_orb_local((idir - 1) * max_nocc + io) = &
              current_orb_local((idir - 1) * max_nocc + io) - 2.0d0 * occ_weight(io) * current_io
          end do
        else if (allocated(dg_frag%momentum_blocks)) then
          tmp_mat(:, :) = (0.0d0, 0.0d0)
          call apply_momentum_blocks(dg_frag, ispin, unit_dir(:, idir), coef_frag_all(1:n, 1:nocc), tmp_mat)
        else if (allocated(dg_frag%momentum_mat_c)) then
          if (.not. allocated(op_mat)) allocate(op_mat(n, n))
          op_mat(:, :) = dg_frag%momentum_mat_c(idir, 1:n, 1:n, ispin)
          call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
                     coef_frag_all, n, (0.0d0, 0.0d0), tmp_mat, n)
        else
          if (.not. allocated(op_mat)) allocate(op_mat(n, n))
          op_mat(:, :) = cmplx(dg_frag%momentum_mat(idir, 1:n, 1:n, ispin), 0.0d0, kind=8)
          call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
                     coef_frag_all(1:n, 1:nocc), n, (0.0d0, 0.0d0), tmp_mat, n)
        end if
        
        if (.not. use_mixed_current) then
          if (use_s_adjoint_current .and. n_pw == 0) then
            if (allocated(dg_frag%momentum_blocks)) then
              if (.not. allocated(op_mat)) allocate(op_mat(n, n))
              op_mat(:, :) = (0.0d0, 0.0d0)
              call copy_momentum_blocks_to_complex_dense(dg_frag, ispin, unit_dir(:, idir), op_mat)
            else if (allocated(dg_frag%momentum_mat_c)) then
              if (.not. allocated(op_mat)) allocate(op_mat(n, n))
              op_mat(:, :) = dg_frag%momentum_mat_c(idir, 1:n, 1:n, ispin)
            else
              if (.not. allocated(op_mat)) allocate(op_mat(n, n))
              op_mat(:, :) = cmplx(dg_frag%momentum_mat(idir, 1:n, 1:n, ispin), 0.0d0, kind=8)
            end if
            do io = 1, nocc
              call apply_overlap_operator(dg_frag, ispin, coef_frag_all(1:n, io), overlap_rhs(1:n, io), .true.)
            end do
            call zgemm('C', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, overlap_rhs, n, &
                       (0.0d0, 0.0d0), tmp_adj, n)
            call solve_overlap_operator_batch(dg_frag, ispin, tmp_adj(1:n, 1:nocc), tmp_solve(1:n, 1:nocc), .true.)
          end if
          ! Factor -2.0: -1 for operator sign convention, 2 for Im[psi*grad psi] normalization
          current_tmp = 0.0d0
          do io = 1, nocc
            if (occ_weight(io) <= 0.0d0) cycle
            current_io = 0.0d0
            if (use_s_adjoint_current .and. n_pw == 0) then
              do ib = 1, n
                current_io = current_io + aimag(conjg(coef_frag_all(ib, io)) * &
                  (0.5d0 * (tmp_mat(ib, io) - tmp_solve(ib, io))))
              end do
            else if (use_s_metric_current) then
              call apply_overlap_operator(dg_frag, ispin, tmp_mat(1:n, io), overlap_vec(1:n), .true.)
              do ib = 1, n
                current_io = current_io + aimag(conjg(coef_frag_all(ib, io)) * overlap_vec(ib))
              end do
            else
              do ib = 1, n
                current_io = current_io + aimag(conjg(coef_frag_all(ib, io)) * tmp_mat(ib, io))
              end do
            end if
            current_tmp = current_tmp + occ_weight(io) * current_io
            if (enable_orbital_probe) then
              current_orb_local((idir - 1) * max_nocc + io) = &
                current_orb_local((idir - 1) * max_nocc + io) - 2.0d0 * occ_weight(io) * current_io
            end if
          end do
        end if
        if (do_para_decomp_sample) then
          if (allocated(dg_frag%momentum_blocks)) then
            do iblk = 1, dg_frag%n_momentum_blocks
              ifrag = dg_frag%momentum_blocks(iblk)%ifrag_row
              jfrag = dg_frag%momentum_blocks(iblk)%ifrag_col
              if (ifrag <= 0 .or. ifrag > dg_frag%n_frag) cycle
              if (jfrag <= 0 .or. jfrag > dg_frag%n_frag) cycle
              nrow_blk = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), &
                             size(dg_frag%momentum_blocks(iblk)%val, 2))
              ncol_blk = min(dg_frag%n_basis(jfrag, ispin), size(dg_frag%index_basis, 1), &
                             size(dg_frag%momentum_blocks(iblk)%val, 3))
              if (nrow_blk <= 0 .or. ncol_blk <= 0) cycle
              do io = 1, nocc
                if (occ_weight(io) <= 0.0d0) cycle
                do ib = 1, nrow_blk
                  i_idx = dg_frag%index_basis(ib, ifrag, ispin)
                  if (i_idx < 1 .or. i_idx > n) cycle
                  do jb = 1, ncol_blk
                    j_idx = dg_frag%index_basis(jb, jfrag, ispin)
                    if (j_idx < 1 .or. j_idx > n) cycle
                    current_elem = occ_weight(io) * conjg(coef_frag_all(i_idx, io)) * &
                                   dg_frag%momentum_blocks(iblk)%val(idir, ib, jb, ispin) * coef_frag_all(j_idx, io)
                    if (i_idx == j_idx) then
                      para_cdiag_local(idir) = para_cdiag_local(idir) - 2.0d0 * aimag(current_elem)
                    else
                      para_coff_local(idir) = para_coff_local(idir) - 2.0d0 * aimag(current_elem)
                      para_coff_abs_local(idir) = para_coff_abs_local(idir) + 2.0d0 * abs(aimag(current_elem))
                    end if
                    if (ifrag == jfrag .and. ib == jb) then
                      para_pdiag_local(idir) = para_pdiag_local(idir) - 2.0d0 * aimag(current_elem)
                    else
                      para_poff_local(idir) = para_poff_local(idir) - 2.0d0 * aimag(current_elem)
                      para_poff_abs_local(idir) = para_poff_abs_local(idir) + 2.0d0 * abs(aimag(current_elem))
                    end if
                  end do
                end do
              end do
            end do
          end if
        end if
        if (enable_transition_probe .and. n_pw == 0 .and. (.not. use_mixed_current)) then
          if (allocated(dg_frag%momentum_blocks)) then
            current_blk_total = (0.0d0, 0.0d0)
            current_blk_diag = (0.0d0, 0.0d0)
            do iblk = 1, dg_frag%n_momentum_blocks
              ifrag = dg_frag%momentum_blocks(iblk)%ifrag_row
              jfrag = dg_frag%momentum_blocks(iblk)%ifrag_col
              if (ifrag <= 0 .or. ifrag > dg_frag%n_frag) cycle
              if (jfrag <= 0 .or. jfrag > dg_frag%n_frag) cycle
              nrow_blk = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), &
                             size(dg_frag%momentum_blocks(iblk)%val, 2))
              ncol_blk = min(dg_frag%n_basis(jfrag, ispin), size(dg_frag%index_basis, 1), &
                             size(dg_frag%momentum_blocks(iblk)%val, 3))
              if (nrow_blk <= 0 .or. ncol_blk <= 0) cycle
              do io = 1, nocc
                occ_i = 1.0d0
                if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
                  if (io <= size(dg_frag%occ_state, 1) .and. ispin <= size(dg_frag%occ_state, 2)) then
                    occ_i = max(0.0d0, dg_frag%occ_state(io, ispin))
                  end if
                else if (allocated(system%rocc)) then
                  if (io <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
                    occ_i = max(0.0d0, system%rocc(io, 1, ispin))
                  end if
                end if
                if (occ_i <= 0.0d0) cycle
                do ib = 1, nrow_blk
                  i_idx = dg_frag%index_basis(ib, ifrag, ispin)
                  if (i_idx < 1 .or. i_idx > n) cycle
                  do jb = 1, ncol_blk
                    j_idx = dg_frag%index_basis(jb, jfrag, ispin)
                    if (j_idx < 1 .or. j_idx > n) cycle
                    current_elem = occ_i * conjg(coef_frag_all(i_idx, io)) * &
                      dg_frag%momentum_blocks(iblk)%val(idir, ib, jb, ispin) * coef_frag_all(j_idx, io)
                    current_blk_total = current_blk_total + current_elem
                    if (ifrag == jfrag .and. ib == jb) then
                      current_blk_diag = current_blk_diag + current_elem
                    end if
                  end do
                end do
              end do
            end do
            current_diag_local(idir) = current_diag_local(idir) - 2.0d0 * aimag(current_blk_diag)
            current_offdiag_local(idir) = current_offdiag_local(idir) - 2.0d0 * aimag(current_blk_total - current_blk_diag)
          end if
        end if
        current_para_local(idir) = current_para_local(idir) - 2.0d0 * current_tmp
      end do
      current_dia_local(1:3) = current_dia_local(1:3) + Ac_tot(1:3) * charge_trace_spin
    end do
    current_local(1:3) = current_para_local(1:3) + current_dia_local(1:3)
    call timer_end(LOG_CALC_CURRENT)
    if (enable_obs_trace) then
      write(*,'(1x,a,i0,a,i0,a,3(1x,1pe12.4),a)') "        obs-trace: rank=", dg_frag%id, " itt=", itt, &
        " current_local=", current_local(1), current_local(2), current_local(3), " stage=current-done"
      flush(6)
    end if
    
      ! Get vector potential at current time for energy calculation
      A_squared = Ac_tot(1)**2 + Ac_tot(2)**2 + Ac_tot(3)**2
      
      require_dense_nl = (.not. allocated(dg_frag%H_mat_blocks)) .or. &
                         (allocated(dg_frag%H_mat_c) .and. allocated(dg_frag%phi_frag_c)) .or. &
                         use_spatial_A .or. do_interface_check
      call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_tot, require_dense_nl)
      has_nonlocal = dg_frag%has_nl_cache

      ! Calculate total energy: E = <ψ|H(t)|ψ>
      ! H(t) = H_0 - i*A(t)·∇ + A²(t)/2 + V_NL(A)
      do ispin = 1, dg_frag%nspin
      if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
        nocc = active_state_cap
      else
        nocc = min(dg_frag%nocc_spin(ispin), active_state_cap)
      end if
      if (nocc <= 0) cycle
      occ_weight(:) = 0.0d0
      do io = 1, nocc
        if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
          if (io <= size(dg_frag%occ_state, 1) .and. ispin <= size(dg_frag%occ_state, 2)) then
            occ_weight(io) = max(0.0d0, dg_frag%occ_state(io, ispin))
          end if
        else
          occ_weight(io) = system%rocc(io, 1, ispin)
        end if
      end do
      call gather_full_coef_view(dg_frag, ispin, n, nocc, coef_frag_view, coef_pw_view, 1, nocc)
      coef_frag_all(1:n, 1:nocc) = coef_frag_view(1:n, 1:nocc)
      if (n_pw > 0) coef_pw_all(1:n_pw, 1:nocc) = coef_pw_view(1:n_pw, 1:nocc)
      if (n_pw > 0) then
        coef_all(1:n_tot, 1:nocc) = (0.0d0, 0.0d0)
        tmp_all(1:n_tot, 1:nocc) = (0.0d0, 0.0d0)
        coef_all(1:n, 1:nocc) = coef_frag_all(1:n, 1:nocc)
        coef_all(n+1:n_tot, 1:nocc) = coef_pw_all(1:n_pw, 1:nocc)
      end if
      if (enable_electron_probe .and. dg_frag%id_frag == 0) then
        if (n_pw == 0) then
          overlap_dense(:, :) = (0.0d0, 0.0d0)
          call copy_overlap_operator_to_dense(dg_frag, ispin, .true., overlap_dense)
        end if
        do io = 1, nocc
          if (n_pw > 0) then
            call apply_overlap_operator(dg_frag, ispin, coef_all(1:n_tot, io), overlap_vec, .true.)
            elec_coef_local = elec_coef_local + occ_weight(io) * &
              real(sum(conjg(coef_all(1:n_tot, io)) * overlap_vec(1:n_tot)))
            elec_plain_local = elec_plain_local + occ_weight(io) * &
              real(sum(conjg(coef_all(1:n_tot, io)) * coef_all(1:n_tot, io)))
            energy_io = real(sum(conjg(coef_all(1:n_tot, io)) * overlap_vec(1:n_tot)), kind=8)
            orth_def_fro_local = orth_def_fro_local + (energy_io - 1.0d0) * (energy_io - 1.0d0)
            orth_diag_min_local = min(orth_diag_min_local, energy_io)
            orth_diag_max_local = max(orth_diag_max_local, energy_io)
            do jo = io + 1, nocc
              call apply_overlap_operator(dg_frag, ispin, coef_all(1:n_tot, jo), overlap_vec, .true.)
              current_io = abs(sum(conjg(coef_all(1:n_tot, io)) * overlap_vec(1:n_tot)))
              orth_def_fro_local = orth_def_fro_local + 2.0d0 * current_io * current_io
              orth_offdiag_max_local = max(orth_offdiag_max_local, current_io)
            end do
          else
            overlap_vec(1:n) = matmul(overlap_dense(1:n, 1:n), coef_frag_all(1:n, io))
            elec_coef_local = elec_coef_local + occ_weight(io) * &
              real(sum(conjg(coef_frag_all(1:n, io)) * overlap_vec(1:n)))
            elec_plain_local = elec_plain_local + occ_weight(io) * &
              real(sum(conjg(coef_frag_all(1:n, io)) * coef_frag_all(1:n, io)))
            energy_io = real(sum(conjg(coef_frag_all(1:n, io)) * overlap_vec(1:n)), kind=8)
            orth_def_fro_local = orth_def_fro_local + (energy_io - 1.0d0) * (energy_io - 1.0d0)
            orth_diag_min_local = min(orth_diag_min_local, energy_io)
            orth_diag_max_local = max(orth_diag_max_local, energy_io)
            do jo = io + 1, nocc
              overlap_vec(1:n) = matmul(overlap_dense(1:n, 1:n), coef_frag_all(1:n, jo))
              current_io = abs(sum(conjg(coef_frag_all(1:n, io)) * overlap_vec(1:n)))
              orth_def_fro_local = orth_def_fro_local + 2.0d0 * current_io * current_io
              orth_offdiag_max_local = max(orth_offdiag_max_local, current_io)
            end do
          end if
        end do
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
                         coef_all, n_tot, (1.0d0, 0.0d0), tmp_all, n_tot)
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
            if (use_energy_components) then
              tmp_probe(:, :) = (0.0d0, 0.0d0)
              if (allocated(dg_frag%H_mat_kinetic_blocks)) then
                call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_kinetic_blocks, ispin, coef_frag_all(1:n, 1:nocc), tmp_probe)
                if (enable_energy_component_probe .and. (itt == 1 .or. itt == 40)) then
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
                if (enable_energy_component_probe .and. (itt == 1 .or. itt == 40)) then
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
              if (allocated(dg_frag%H_nl_cache) .and. .not. dg_frag%use_buffered_basis) then
                if (use_energy_components) then
                  tmp_probe(:, :) = matmul(dg_frag%H_nl_cache(1:n, 1:n, ispin), coef_frag_all(1:n, 1:nocc))
                  do io = 1, nocc
                    energy_nl_local = energy_nl_local + occ_weight(io) * &
                      sum(real(conjg(coef_frag_all(1:n, io)) * tmp_probe(1:n, io)))
                    if (allocated(nl_state_local)) then
                      nl_state_local(io) = nl_state_local(io) + occ_weight(io) * &
                        sum(real(conjg(coef_frag_all(1:n, io)) * tmp_probe(1:n, io)))
                    end if
                  end do
                end if
                tmp_mat(:, :) = tmp_mat(:, :) + &
                  matmul(dg_frag%H_nl_cache(1:n, 1:n, ispin), coef_frag_all(1:n, 1:nocc))
              else
                if (use_energy_components) then
                  tmp_probe(:, :) = (0.0d0, 0.0d0)
                  call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_frag_all(1:n, 1:nocc), &
                    tmp_probe)
                  do io = 1, nocc
                    energy_nl_local = energy_nl_local + occ_weight(io) * &
                      sum(real(conjg(coef_frag_all(1:n, io)) * tmp_probe(1:n, io)))
                    if (allocated(nl_state_local)) then
                      nl_state_local(io) = nl_state_local(io) + occ_weight(io) * &
                        sum(real(conjg(coef_frag_all(1:n, io)) * tmp_probe(1:n, io)))
                    end if
                  end do
                end if
                call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_frag_all(1:n, 1:nocc), &
                  tmp_mat)
              end if
            end if
            if (use_energy_components) then
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
          if (use_energy_components) then
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
          energy_io = 0.0d0
          do ib = 1, n_tot
            energy_io = energy_io + occ_weight(io) * real(conjg(coef_all(ib, io)) * tmp_all(ib, io))
          end do
          energy_tmp = energy_tmp + energy_io
          if (enable_orbital_probe) energy_orb_local(io) = energy_orb_local(io) + energy_io
        end do
      else
        if (.not. allocated(dg_frag%momentum_blocks) .or. use_spatial_A) then
        call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
             coef_frag_all, n, (0.0d0, 0.0d0), tmp_mat, n)
        end if
        energy_tmp = 0.0d0
        do io = 1, nocc
          energy_io = 0.0d0
          do ib = 1, n
            energy_io = energy_io + occ_weight(io) * real(conjg(coef_frag_all(ib, io)) * tmp_mat(ib, io))
          end do
          energy_tmp = energy_tmp + energy_io
          if (enable_orbital_probe) energy_orb_local(io) = energy_orb_local(io) + energy_io
        end do
      end if
        energy_local = energy_local + energy_tmp

      if (n_pw > 0 .and. allocated(dg_frag%H_mat_frag_pw) .and. mixed_fp_coupling_active(dg_frag, ispin) .and. &
          dg_frag%mixed_basis_ready .and. allocated(dg_frag%coef_mix) .and. allocated(dg_frag%mixed_transform) .and. &
          .not. use_spatial_A) then
        n_mix = min(dg_frag%mixed_basis_dim(ispin), size(dg_frag%coef_mix, 1), size(dg_frag%mixed_transform, 2))
        if (n_mix > 0) then
          coef_all(1:n_tot, 1:nocc) = (0.0d0, 0.0d0)
          coef_all(1:n_tot, 1:nocc) = matmul(dg_frag%mixed_transform(1:n_tot, 1:n_mix, ispin), &
                                             dg_frag%coef_mix(1:n_mix, 1:nocc, ispin))
          tmp_all(1:n_tot, 1:nocc) = (0.0d0, 0.0d0)
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
                         coef_all, n_tot, (1.0d0, 0.0d0), tmp_all, n_tot)
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
          do io = 1, nocc
            energy_mix_io = 0.0d0
            do ib = 1, n_tot
              energy_mix_io = energy_mix_io + occ_weight(io) * real(conjg(coef_all(ib, io)) * tmp_all(ib, io))
            end do
            energy_mix_local = energy_mix_local + energy_mix_io
          end do
        end if
      end if
      end do

      if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) then
        do ispin = 1, dg_frag%nspin
          if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
            nocc = active_state_cap
          else
            nocc = min(dg_frag%nocc_spin(ispin), active_state_cap)
          end if
          if (nocc <= 0) cycle
          call gather_full_coef_view(dg_frag, ispin, n, nocc, coef_frag_view, coef_pw_view, 1, nocc)
          coef_pw_all(1:n_pw, 1:nocc) = coef_pw_view(1:n_pw, 1:nocc)
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

    if (use_energy_components .and. allocated(dg_frag%H_mat_kinetic_blocks)) then
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
    if (allocated(coef_frag_view)) deallocate(coef_frag_view)
    if (allocated(coef_pw_view)) deallocate(coef_pw_view)
    if (allocated(coef_all)) deallocate(coef_all)
    if (allocated(tmp_all)) deallocate(tmp_all)
    if (allocated(overlap_vec)) deallocate(overlap_vec)
    if (allocated(overlap_dense)) deallocate(overlap_dense)
    if (allocated(overlap_rhs)) deallocate(overlap_rhs)
    if (allocated(tmp_adj)) deallocate(tmp_adj)
    if (allocated(tmp_solve)) deallocate(tmp_solve)
    ! Cache retained for reuse

  1000 continue
    
    if (n_pw == 0 .and. enable_realspace_probe) then
      call compute_realspace_nonlocal_current_probe(dg_frag, system, mg, ppg, Ac_tot, current_nl_rs_sum)
      current_local(1:3) = current_local(1:3) + current_nl_rs_sum(1:3)
      call compute_realspace_energy_probe(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, Vh, Vxc, Vpsl, &
                                          energy_kin_local, energy_local, energy_kin_rs_sum, energy_one_rs_sum, energy_nl_rs_sum)
    end if

    if (enable_energy_component_probe .and. n_pw == 0 .and. (itt == 1 .or. itt == 40)) then
      write(*,'(1x,a,i0,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
        "        energy-local probe: rank=", dg_frag%id, " itt=", itt, &
        " local_total=", energy_local, " local_static=", energy_static_local, " local_kin=", energy_kin_local
      flush(6)
    end if

    ! MPI aggregation: sum local contributions from all ranks
    call comm_summation(current_local, dg_frag%current, 3, dg_frag%icomm)
    current_reduce_sum(:) = dg_frag%current(:)
    current_para_reduce_sum(:) = 0.0d0
    current_dia_reduce_sum(:) = 0.0d0
    current_nl_reduce_sum(:) = 0.0d0
    if (enable_current_component_trace) then
      call comm_summation(current_para_local, current_para_reduce_sum, 3, dg_frag%icomm)
      call comm_summation(current_dia_local, current_dia_reduce_sum, 3, dg_frag%icomm)
      if (n_pw == 0 .and. enable_realspace_probe) then
        call comm_summation(current_nl_rs_sum, current_nl_reduce_sum, 3, dg_frag%icomm)
      end if
    end if
    if (enable_transition_probe) then
      call comm_summation(current_diag_local, current_diag_sum, 3, dg_frag%icomm)
      call comm_summation(current_offdiag_local, current_offdiag_sum, 3, dg_frag%icomm)
    end if
    call comm_summation(energy_local, dg_frag%total_energy, dg_frag%icomm)
    call comm_summation(energy_mix_local, dg_frag%total_energy_mix, dg_frag%icomm)
    call comm_summation(pw_weight_local, dg_frag%pw_weight_raw, dg_frag%icomm)
    if (enable_electron_probe) then
      call comm_summation(elec_coef_local, elec_coef_sum, dg_frag%icomm)
      call comm_summation(elec_plain_local, elec_plain_sum, dg_frag%icomm)
      call comm_summation(orth_def_fro_local, orth_def_fro_sum, dg_frag%icomm)
      call comm_summation(orth_offdiag_max_local, orth_offdiag_max_sum, dg_frag%icomm)
      call comm_summation(orth_diag_min_local, orth_diag_min_sum, dg_frag%icomm)
      call comm_summation(orth_diag_max_local, orth_diag_max_sum, dg_frag%icomm)
    end if
    if (enable_orbital_probe) then
      call comm_summation(current_orb_local, current_orb_sum, 3 * max_nocc, dg_frag%icomm)
      call comm_summation(energy_orb_local, energy_orb_sum, max_nocc, dg_frag%icomm)
    end if
    if (enable_excitation_trace) then
      call comm_summation(dm_diag_norm_local, dm_diag_norm_sum, dg_frag%icomm)
      call comm_summation(dm_offdiag_norm_local, dm_offdiag_norm_sum, dg_frag%icomm)
      call comm_summation(dm_offdiag_max_local, dm_offdiag_max_sum, dg_frag%icomm)
      dm_diag_norm_sum = dm_diag_norm_sum / real(max(1, dg_frag%isize), 8)
      dm_offdiag_norm_sum = dm_offdiag_norm_sum / real(max(1, dg_frag%isize), 8)
      dm_offdiag_max_sum = dm_offdiag_max_sum / real(max(1, dg_frag%isize), 8)
      dm_offdiag_ratio = 0.0d0
      if (dm_diag_norm_sum > 0.0d0) dm_offdiag_ratio = dm_offdiag_norm_sum / dm_diag_norm_sum
      if ((.not. dm_excitation_ref_initialized) .or. itt == 0) then
        dm_diag_norm_ref = dm_diag_norm_sum
        dm_offdiag_norm_ref = dm_offdiag_norm_sum
        dm_offdiag_ratio_ref = dm_offdiag_ratio
        dm_excitation_ref_initialized = .true.
      end if
      dm_diag_delta = dm_diag_norm_sum - dm_diag_norm_ref
      dm_offdiag_delta = dm_offdiag_norm_sum - dm_offdiag_norm_ref
      dm_offdiag_ratio_delta = dm_offdiag_ratio - dm_offdiag_ratio_ref
      if (dg_frag%id == 0 .and. (itt == 0 .or. itt == 1 .or. mod(itt - 1, excitation_trace_stride) == 0)) then
        write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
          "        dg-excitation-trace: itt=", itt, &
          " dm_diag_norm=", dm_diag_norm_sum, " dm_offdiag_norm=", dm_offdiag_norm_sum, &
          " dm_offdiag_max=", dm_offdiag_max_sum, " offdiag/diag=", dm_offdiag_ratio, &
          " d_diag=", dm_diag_delta, " d_offdiag=", dm_offdiag_delta, " d_ratio=", dm_offdiag_ratio_delta
        flush(6)
      end if
    end if
    if (enable_para_decomp_trace) then
      call comm_summation(para_cdiag_local, para_cdiag_sum, 3, dg_frag%icomm)
      call comm_summation(para_coff_local, para_coff_sum, 3, dg_frag%icomm)
      call comm_summation(para_pdiag_local, para_pdiag_sum, 3, dg_frag%icomm)
      call comm_summation(para_poff_local, para_poff_sum, 3, dg_frag%icomm)
      call comm_summation(para_coff_abs_local, para_coff_abs_sum, 3, dg_frag%icomm)
      call comm_summation(para_poff_abs_local, para_poff_abs_sum, 3, dg_frag%icomm)
      para_cdiag_sum(:) = para_cdiag_sum(:) / real(max(1, dg_frag%isize), 8)
      para_coff_sum(:) = para_coff_sum(:) / real(max(1, dg_frag%isize), 8)
      para_pdiag_sum(:) = para_pdiag_sum(:) / real(max(1, dg_frag%isize), 8)
      para_poff_sum(:) = para_poff_sum(:) / real(max(1, dg_frag%isize), 8)
      para_coff_abs_sum(:) = para_coff_abs_sum(:) / real(max(1, dg_frag%isize), 8)
      para_poff_abs_sum(:) = para_poff_abs_sum(:) / real(max(1, dg_frag%isize), 8)
      para_cdiag_now(:) = para_cdiag_sum(:) / system%det_a
      para_coff_now(:) = para_coff_sum(:) / system%det_a
      para_pdiag_now(:) = para_pdiag_sum(:) / system%det_a
      para_poff_now(:) = para_poff_sum(:) / system%det_a
      para_coff_abs_now(:) = para_coff_abs_sum(:) / system%det_a
      para_poff_abs_now(:) = para_poff_abs_sum(:) / system%det_a
      call sym_vector_xyz(para_cdiag_now)
      call sym_vector_xyz(para_coff_now)
      call sym_vector_xyz(para_pdiag_now)
      call sym_vector_xyz(para_poff_now)
      call sym_vector_xyz(para_coff_abs_now)
      call sym_vector_xyz(para_poff_abs_now)
      if (dg_frag%id == 0 .and. (itt == 0 .or. itt == 1 .or. mod(itt - 1, para_decomp_trace_stride) == 0)) then
        write(*,'(1x,a,i0,a,3(1x,1pe14.6),a,3(1x,1pe14.6),a,3(1x,1pe14.6),a,3(1x,1pe14.6))') &
          "        dg-para-decomp-trace: itt=", itt, &
          " cdiag=", para_cdiag_now(1), para_cdiag_now(2), para_cdiag_now(3), &
          " coff=", para_coff_now(1), para_coff_now(2), para_coff_now(3), &
          " pdiag=", para_pdiag_now(1), para_pdiag_now(2), para_pdiag_now(3), &
          " poff=", para_poff_now(1), para_poff_now(2), para_poff_now(3)
        write(*,'(1x,a,i0,a,3(1x,1pe14.6),a,3(1x,1pe14.6))') &
          "        dg-para-decomp-abs: itt=", itt, &
          " coff_abs=", para_coff_abs_now(1), para_coff_abs_now(2), para_coff_abs_now(3), &
          " poff_abs=", para_poff_abs_now(1), para_poff_abs_now(2), para_poff_abs_now(3)
        flush(6)
      end if
    end if
    frag_reduce_factor = real(max(1, dg_frag%isize_frag), 8)
    dg_frag%total_energy = dg_frag%total_energy / frag_reduce_factor
    dg_frag%total_energy_mix = dg_frag%total_energy_mix / frag_reduce_factor
    if (enable_electron_probe) then
      elec_coef_sum = elec_coef_sum / frag_reduce_factor
      elec_plain_sum = elec_plain_sum / frag_reduce_factor
      orth_def_fro_sum = orth_def_fro_sum / frag_reduce_factor
      orth_diag_min_sum = orth_diag_min_sum / frag_reduce_factor
      orth_diag_max_sum = orth_diag_max_sum / frag_reduce_factor
    end if
    if (enable_orbital_probe) energy_orb_sum(:) = energy_orb_sum(:) / frag_reduce_factor
    ! Observables are currently evaluated from full gathered coefficient views on each rank.
    ! The world reduction therefore sums replicated current contributions and must be averaged back.
    dg_frag%current(:) = dg_frag%current(:) / real(max(1, dg_frag%isize), 8)
    current_after_isize(:) = dg_frag%current(:)
    if (enable_current_component_trace) then
      current_para_after_isize(:) = current_para_reduce_sum(:) / real(max(1, dg_frag%isize), 8)
      current_dia_after_isize(:) = current_dia_reduce_sum(:) / real(max(1, dg_frag%isize), 8)
      current_nl_after_isize(:) = current_nl_reduce_sum(:) / real(max(1, dg_frag%isize), 8)
    end if
    if (use_energy_components) then
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
      if (dg_frag%use_buffered_basis .and. n_pw == 0) then
        energy_kin_sum = energy_kin_rs_sum
        energy_kin_avg = energy_kin_rs_sum
        energy_nl_sum = energy_nl_rs_sum
        energy_nl_avg = energy_nl_rs_sum
      end if
    end if

    ! PW weight is replicated over all ranks after the world reduction.
    dg_frag%pw_weight_raw = dg_frag%pw_weight_raw / real(max(1, dg_frag%isize), 8)
    if (enable_transition_probe) then
      current_diag_sum(:) = current_diag_sum(:) / real(max(1, dg_frag%isize), 8)
      current_offdiag_sum(:) = current_offdiag_sum(:) / real(max(1, dg_frag%isize), 8)
    end if
    if (enable_orbital_probe) current_orb_sum(:) = current_orb_sum(:) / real(max(1, dg_frag%isize), 8)
    if (use_energy_components) then
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
    if (use_energy_components) then
      if (dg_frag%use_buffered_basis .and. n_pw == 0) then
        dg_frag%energy_kinetic = energy_kin_rs_sum
      else
        dg_frag%energy_kinetic = energy_kin_sum
      end if
      dg_frag%energy_nonlocal = energy_nl_sum
      dg_frag%energy_Ap = energy_ap_sum
      dg_frag%energy_A2 = energy_a2_sum
    end if
    if (dg_frag%use_buffered_basis .and. n_pw == 0) then
      dg_frag%total_energy = energy_one_rs_sum
    end if

    ! DG momentum matrices already include hvol in their matrix elements, so the
    ! coefficient contraction yields a cell-integrated current. Normalize by the
    ! simulation-cell volume to obtain the current density.
    dg_frag%current(:) = dg_frag%current(:) / system%det_a
    if (enable_current_component_trace) then
      current_para_reduce_sum(:) = current_para_after_isize(:) / system%det_a
      current_dia_reduce_sum(:) = current_dia_after_isize(:) / system%det_a
      current_nl_reduce_sum(:) = current_nl_after_isize(:) / system%det_a
      if (dg_frag%id == 0) then
        if (itt == 0 .or. itt == 1 .or. mod(itt - 1, current_component_stride) == 0) then
          write(*,'(1x,a,i0,a,3(1x,1pe14.6),a,3(1x,1pe14.6),a,3(1x,1pe14.6),a,3(1x,1pe14.6))') &
            "        current-component-trace: itt=", itt, &
            " para=", current_para_reduce_sum(1), current_para_reduce_sum(2), current_para_reduce_sum(3), &
            " dia=", current_dia_reduce_sum(1), current_dia_reduce_sum(2), current_dia_reduce_sum(3), &
            " nlrs=", current_nl_reduce_sum(1), current_nl_reduce_sum(2), current_nl_reduce_sum(3), &
            " total=", dg_frag%current(1), dg_frag%current(2), dg_frag%current(3)
          flush(6)
        end if
      end if
    end if
    if (enable_current_norm_trace .and. dg_frag%id == 0) then
      if (itt == 0 .or. itt == 1 .or. mod(itt - 1, current_norm_stride) == 0) then
        write(*,'(1x,a,i0,a,3(1x,1pe14.6),a,3(1x,1pe14.6),a,3(1x,1pe14.6),a,1pe14.6,a,i0,a,i0)') &
          "        current-norm-trace: itt=", itt, &
          " reduce=", current_reduce_sum(1), current_reduce_sum(2), current_reduce_sum(3), &
          " after_isize=", current_after_isize(1), current_after_isize(2), current_after_isize(3), &
          " after_det_a=", dg_frag%current(1), dg_frag%current(2), dg_frag%current(3), &
          " det_a=", system%det_a, " isize=", dg_frag%isize, " isize_frag=", dg_frag%isize_frag
        flush(6)
      end if
    end if
    call sym_vector_xyz(dg_frag%current)
    if (enable_orbital_probe) current_orb_sum(:) = current_orb_sum(:) / system%det_a
    if (enable_transition_probe) then
      current_diag_sum(:) = current_diag_sum(:) / system%det_a
      current_offdiag_sum(:) = current_offdiag_sum(:) / system%det_a
    end if

    if (enable_orbital_probe .and. dg_frag%id == 0) then
      if (mod(itt - 1, probe_stride) == 0) then
        write(*,'(1x,a,i0,a)') "        orbital-probe: itt=", itt, " (io, E, Jx, Jy, Jz)"
        do io = 1, min(max_nocc, probe_nprint)
          write(*,'(1x,i6,4(1x,1pe14.6))') io, energy_orb_sum(io), current_orb_sum(io), &
            current_orb_sum(max_nocc + io), current_orb_sum(2 * max_nocc + io)
        end do
        flush(6)
      end if
    end if
    if (enable_transition_probe .and. dg_frag%id == 0) then
      if (mod(itt - 1, transition_stride) == 0) then
        write(*,'(1x,a,i0,a,3(1x,1pe14.6),a,3(1x,1pe14.6),a,3(1x,1pe14.6))') &
          "        transition-probe: itt=", itt, " Jdiag=", current_diag_sum(1), current_diag_sum(2), current_diag_sum(3), &
          " Joffdiag=", current_offdiag_sum(1), current_offdiag_sum(2), current_offdiag_sum(3), &
          " Jtotal=", dg_frag%current(1), dg_frag%current(2), dg_frag%current(3)
        flush(6)
      end if
    end if
    if (enable_electron_probe .and. dg_frag%id == 0) then
      if (itt == 1 .or. mod(itt, 10) == 0) then
        if (dg_frag%use_buffered_basis) then
          write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
            "        electron-probe(buffered): itt=", itt, " coeff_metric_S=", elec_coef_sum, " coeff_metric_2=", elec_plain_sum, &
            " Ne_raw=", dg_frag%elec_num_raw, " rho_scale=", dg_frag%rho_scale_factor, &
            " ortho_fro=", sqrt(max(orth_def_fro_sum, 0.0d0)), " ortho_offmax=", orth_offdiag_max_sum, &
            " ortho_diag_min=", orth_diag_min_sum
          write(*,'(1x,a,1pe14.6)') "        electron-probe(buffered) ortho_diag_max=", orth_diag_max_sum
        else
          write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
            "        electron-probe: itt=", itt, " Ne_coef_S=", elec_coef_sum, " Ne_coef_2=", elec_plain_sum, &
            " Ne_raw=", dg_frag%elec_num_raw, " rho_scale=", dg_frag%rho_scale_factor, &
            " ortho_fro=", sqrt(max(orth_def_fro_sum, 0.0d0)), " ortho_offmax=", orth_offdiag_max_sum, &
            " ortho_diag_min=", orth_diag_min_sum
          write(*,'(1x,a,1pe14.6)') "        electron-probe ortho_diag_max=", orth_diag_max_sum
        end if
        flush(6)
      end if
    end if

    if (allocated(current_orb_local)) deallocate(current_orb_local)
    if (allocated(current_orb_sum)) deallocate(current_orb_sum)
    if (allocated(energy_orb_local)) deallocate(energy_orb_local)
    if (allocated(energy_orb_sum)) deallocate(energy_orb_sum)
    ! Store in rt structure for output
    rt%curr(:, itt) = dg_frag%current(:)
    
  end subroutine calculate_observables

  subroutine print_initial_electron_probe(dg_frag, system, mg, rho)
    use structures
    use communication, only: comm_summation
    use rt_dg_fragment_ops, only: apply_overlap_operator, gather_full_coef_view, copy_overlap_operator_to_dense
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_scalar),         intent(in)    :: rho

    integer :: ispin, io, n, n_pw, n_tot, nocc, nocc_report
    integer :: env_len, env_status
    character(len=64) :: env_val
    logical :: enable_electron_probe
    real(8) :: elec_coef_local, elec_plain_local, elec_rho_local
    real(8) :: elec_coef_sum, elec_plain_sum, elec_rho_sum, occ_i
    integer, allocatable :: occ_idx_report(:)
    real(8), allocatable :: occ_val_report(:), occ_sdiag_report(:), occ_c2_report(:)
    real(8), allocatable :: occ_frag_report(:), occ_pw_report(:)
    complex(8), allocatable :: coef_all(:,:), coef_frag_all(:,:), coef_pw_all(:,:), overlap_vec(:), overlap_dense(:,:)
    complex(8), allocatable :: coef_frag_view(:,:), coef_pw_view(:,:)

    enable_electron_probe = .false.
    call get_environment_variable("SALMON_DG_ELECTRON_PROBE", env_val, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      if (env_val(1:1) == '1' .or. env_val(1:1) == 'y' .or. env_val(1:1) == 'Y' .or. &
          env_val(1:1) == 't' .or. env_val(1:1) == 'T') then
        enable_electron_probe = .true.
      end if
    end if
    if (.not. enable_electron_probe) return

    n = dg_frag%n_mat_max
    if (n <= 0) return
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n + n_pw

    allocate(coef_frag_all(n, 1))
    if (n_pw > 0) then
      allocate(coef_pw_all(n_pw, 1), coef_all(n_tot, 1), overlap_vec(n_tot))
    else
      allocate(overlap_vec(n))
    end if

    elec_coef_local = 0.0d0
    elec_plain_local = 0.0d0
    elec_rho_local = sum(rho%f(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))) * system%hvol
    nocc_report = 0
    allocate(occ_idx_report(8), occ_val_report(8), occ_sdiag_report(8), occ_c2_report(8))
    allocate(occ_frag_report(8), occ_pw_report(8))
    occ_idx_report(:) = 0
    occ_val_report(:) = 0.0d0
    occ_sdiag_report(:) = 0.0d0
    occ_c2_report(:) = 0.0d0
    occ_frag_report(:) = 0.0d0
    occ_pw_report(:) = 0.0d0

    do ispin = 1, dg_frag%nspin
      if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
        nocc = dg_frag%nstate_tot
      else
        nocc = min(dg_frag%nocc_spin(ispin), dg_frag%nstate_tot)
      end if
      if (nocc <= 0) cycle
      call gather_full_coef_view(dg_frag, ispin, n, nocc, coef_frag_view, coef_pw_view, 1, nocc)
      if (dg_frag%id_frag == 0 .and. n_pw == 0) then
        if (.not. allocated(overlap_dense)) allocate(overlap_dense(n, n))
        overlap_dense(:, :) = (0.0d0, 0.0d0)
        call copy_overlap_operator_to_dense(dg_frag, ispin, .true., overlap_dense)
      end if
      do io = 1, nocc
        occ_i = 1.0d0
        if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
          if (io <= size(dg_frag%occ_state, 1) .and. ispin <= size(dg_frag%occ_state, 2)) then
            occ_i = max(0.0d0, dg_frag%occ_state(io, ispin))
          end if
        else if (allocated(system%rocc)) then
          if (io <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
            occ_i = max(0.0d0, system%rocc(io, 1, ispin))
          end if
        end if
        if (occ_i <= 0.0d0) cycle

        if (dg_frag%id_frag /= 0) cycle
        if (n_pw > 0) then
          coef_all(:, 1) = (0.0d0, 0.0d0)
          coef_all(1:n, 1) = coef_frag_view(1:n, io)
          coef_all(n+1:n_tot, 1) = coef_pw_view(1:n_pw, io)
          call apply_overlap_operator(dg_frag, ispin, coef_all(1:n_tot, 1), overlap_vec(1:n_tot), .true.)
          if (dg_frag%id == 0 .and. nocc_report < size(occ_idx_report)) then
            nocc_report = nocc_report + 1
            occ_idx_report(nocc_report) = io
            occ_val_report(nocc_report) = occ_i
            occ_sdiag_report(nocc_report) = real(sum(conjg(coef_all(1:n_tot, 1)) * overlap_vec(1:n_tot)))
            occ_c2_report(nocc_report) = real(sum(conjg(coef_all(1:n_tot, 1)) * coef_all(1:n_tot, 1)))
            occ_frag_report(nocc_report) = real(sum(conjg(coef_all(1:n, 1)) * overlap_vec(1:n)), kind=8)
            occ_pw_report(nocc_report) = real(sum(conjg(coef_all(n+1:n_tot, 1)) * overlap_vec(n+1:n_tot)), kind=8)
          end if
          elec_coef_local = elec_coef_local + occ_i * real(sum(conjg(coef_all(1:n_tot, 1)) * overlap_vec(1:n_tot)))
          elec_plain_local = elec_plain_local + occ_i * real(sum(conjg(coef_all(1:n_tot, 1)) * coef_all(1:n_tot, 1)))
        else
          coef_frag_all(1:n, 1) = coef_frag_view(1:n, io)
          if (dg_frag%id == 0 .and. nocc_report < size(occ_idx_report)) then
            nocc_report = nocc_report + 1
            occ_idx_report(nocc_report) = io
            occ_val_report(nocc_report) = occ_i
            occ_sdiag_report(nocc_report) = real(sum(conjg(coef_frag_all(1:n, 1)) * &
              matmul(overlap_dense(1:n, 1:n), coef_frag_all(1:n, 1))))
            occ_c2_report(nocc_report) = real(sum(conjg(coef_frag_all(1:n, 1)) * coef_frag_all(1:n, 1)))
          end if
          elec_coef_local = elec_coef_local + occ_i * real(sum(conjg(coef_frag_all(1:n, 1)) * &
            matmul(overlap_dense(1:n, 1:n), coef_frag_all(1:n, 1))))
          elec_plain_local = elec_plain_local + occ_i * real(sum(conjg(coef_frag_all(1:n, 1)) * coef_frag_all(1:n, 1)))
        end if
      end do
    end do

    call comm_summation(elec_coef_local, elec_coef_sum, dg_frag%icomm)
    call comm_summation(elec_plain_local, elec_plain_sum, dg_frag%icomm)
    call comm_summation(elec_rho_local, elec_rho_sum, dg_frag%icomm)
    elec_coef_sum = elec_coef_sum / real(max(1, dg_frag%isize_frag), 8)
    elec_plain_sum = elec_plain_sum / real(max(1, dg_frag%isize_frag), 8)
    elec_rho_sum = elec_rho_sum / real(max(1, dg_frag%isize_frag), 8)

    if (dg_frag%id == 0) then
      if (dg_frag%use_buffered_basis) then
        write(*,'(1x,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
          "        electron-probe-t0(buffered): coeff_metric_S=", elec_coef_sum, " coeff_metric_2=", elec_plain_sum, " Ne_rho=", elec_rho_sum
      else
        write(*,'(1x,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
          "        electron-probe-t0: Ne_coef_S=", elec_coef_sum, " Ne_coef_2=", elec_plain_sum, " Ne_rho=", elec_rho_sum
      end if
      do io = 1, nocc_report
        write(*,'(1x,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
          "        electron-probe-t0 occ-state=", occ_idx_report(io), " occ=", occ_val_report(io), &
          " sdiag=", occ_sdiag_report(io), " c2=", occ_c2_report(io), &
          " frag_w=", occ_frag_report(io), " pw_w=", occ_pw_report(io)
      end do
      flush(6)
    end if

    if (allocated(coef_all)) deallocate(coef_all)
    if (allocated(coef_pw_all)) deallocate(coef_pw_all)
    if (allocated(coef_frag_all)) deallocate(coef_frag_all)
    if (allocated(coef_frag_view)) deallocate(coef_frag_view)
    if (allocated(coef_pw_view)) deallocate(coef_pw_view)
    if (allocated(overlap_vec)) deallocate(overlap_vec)
    if (allocated(overlap_dense)) deallocate(overlap_dense)
    if (allocated(occ_idx_report)) deallocate(occ_idx_report)
    if (allocated(occ_val_report)) deallocate(occ_val_report)
    if (allocated(occ_sdiag_report)) deallocate(occ_sdiag_report)
    if (allocated(occ_c2_report)) deallocate(occ_c2_report)
    if (allocated(occ_frag_report)) deallocate(occ_frag_report)
    if (allocated(occ_pw_report)) deallocate(occ_pw_report)
  end subroutine print_initial_electron_probe

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

  integer function map_global_to_phi_box_coord_obs_fragment(dg_frag, ifrag, axis, ig, lb, ub) result(iloc)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, axis, ig, lb, ub
    integer :: ig_wrap, support_lo, support_len

    if (dg_frag%use_buffered_basis) then
      support_lo = dg_frag%basis_support_lo(axis, ifrag)
      support_len = dg_frag%basis_support_hi(axis, ifrag) - dg_frag%basis_support_lo(axis, ifrag) + 1
      ig_wrap = support_lo + modulo(ig - support_lo, support_len)
      iloc = map_global_to_phi_box_coord_obs(ig_wrap, lb, ub, dg_frag%lgnum_total(axis))
    else
      iloc = map_global_to_phi_box_coord_obs(ig, lb, ub, dg_frag%lgnum_total(axis))
    end if
  end function map_global_to_phi_box_coord_obs_fragment

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
    if (dg_frag%use_buffered_basis) then
      iorg(:) = dg_frag%basis_support_lo(:, ifrag)
      ndom(:) = dg_frag%basis_support_hi(:, ifrag) - dg_frag%basis_support_lo(:, ifrag) + 1
    else
      iorg(:) = dg_frag%ixyz_frag(:, ifrag)
      ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    end if
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
            bx = map_global_to_phi_box_coord_obs_fragment(dg_frag, ifrag, 1, gx, p_lb1, p_ub1)
            by = map_global_to_phi_box_coord_obs_fragment(dg_frag, ifrag, 2, gy, p_lb2, p_ub2)
            bz = map_global_to_phi_box_coord_obs_fragment(dg_frag, ifrag, 3, gz, p_lb3, p_ub3)
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
            bx = map_global_to_phi_box_coord_obs_fragment(dg_frag, ifrag, 1, gx, p_lb1, p_ub1)
            by = map_global_to_phi_box_coord_obs_fragment(dg_frag, ifrag, 2, gy, p_lb2, p_ub2)
            bz = map_global_to_phi_box_coord_obs_fragment(dg_frag, ifrag, 3, gz, p_lb3, p_ub3)
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
    if (dg_frag%use_buffered_basis) then
      iorg(:) = dg_frag%basis_support_lo(:, ifrag)
      ndom(:) = dg_frag%basis_support_hi(:, ifrag) - dg_frag%basis_support_lo(:, ifrag) + 1
    else
      iorg(:) = dg_frag%ixyz_frag(:, ifrag)
      ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    end if
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
            bx = map_global_to_phi_box_coord_obs_fragment(dg_frag, ifrag, 1, gx, p_lb1, p_ub1)
            by = map_global_to_phi_box_coord_obs_fragment(dg_frag, ifrag, 2, gy, p_lb2, p_ub2)
            bz = map_global_to_phi_box_coord_obs_fragment(dg_frag, ifrag, 3, gz, p_lb3, p_ub3)
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
            bx = map_global_to_phi_box_coord_obs_fragment(dg_frag, ifrag, 1, gx, p_lb1, p_ub1)
            by = map_global_to_phi_box_coord_obs_fragment(dg_frag, ifrag, 2, gy, p_lb2, p_ub2)
            bz = map_global_to_phi_box_coord_obs_fragment(dg_frag, ifrag, 3, gz, p_lb3, p_ub3)
            if (bx == 0 .or. by == 0 .or. bz == 0) cycle
            integral = integral + cmplx(dg_frag%phi_frag(bx, by, bz, io, i_local), 0.0d0, kind=8) * field(gx, gy, gz) * hvol
          end do
        end do
      end do
    end if
  end subroutine integrate_local_basis_with_field_local

  subroutine compute_realspace_energy_probe(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, Vh, Vxc, Vpsl, energy_kin_mat, energy_one_mat, kin_sum_out, one_sum_out, nl_sum_out)
    use structures
    use communication, only: comm_summation, comm_is_root, COMM_GROUP_NULL
    use parallelization, only: nproc_id_global
    use rt_dg_fragment_ops, only: apply_nonlocal_pp_projector_batch
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_pp_grid),        intent(in)    :: ppg
    real(8),                intent(in)    :: Ac_tot(3)
    integer,                intent(in)    :: itt
    type(s_scalar),         intent(in)    :: Vh, Vpsl
    type(s_scalar),         intent(in)    :: Vxc(system%nspin)
    real(8),                intent(in)    :: energy_kin_mat, energy_one_mat
    real(8),                intent(out)   :: kin_sum_out, one_sum_out, nl_sum_out

    integer :: ispin, io, ifrag, i_local, jo, nbf, ig_j, nocc
    integer :: core_s(3), core_e(3), ov_s(3), ov_e(3)
    integer :: gx, gy, gz, ixg, iyg, izg
    logical :: has_overlap
    real(8), allocatable :: V_total(:,:,:)
    complex(8), allocatable :: T_phi(:,:,:), H_phi(:,:,:)
    complex(8), allocatable :: psi_frag(:,:,:), tpsi_frag(:,:,:), hpsi_frag(:,:,:), vnlpsi_frag(:,:,:)
    complex(8), allocatable :: coef_probe(:,:), nl_probe(:,:)
    complex(8) :: coeff, ztmp, phi_val
    real(8) :: kin_local, one_local, kin_sum, one_sum, nl_sum
    real(8) :: kin_frag_sum, one_frag_sum
    real(8) :: frag_reduce_factor
    real(8) :: occ_probe
    logical :: use_buffered_direct_orbitals
    real(8), allocatable :: kin_state_local(:), kin_state_sum(:), one_state_local(:), one_state_sum(:)
    real(8), allocatable :: nl_state_local(:), nl_state_sum(:)

    kin_sum_out = 0.0d0
    one_sum_out = 0.0d0
    nl_sum_out = 0.0d0
    if (dg_frag%use_plane_wave_basis) return
    if (.not. dg_frag%has_real_space_basis) return
    if (dg_frag%n_mat_max <= 0) return

    kin_local = 0.0d0
    one_local = 0.0d0
    if (dg_frag%use_buffered_basis .and. (itt == 0 .or. itt == 1 .or. itt == 40)) then
      allocate(kin_state_local(dg_frag%nstate_tot), kin_state_sum(dg_frag%nstate_tot))
      allocate(one_state_local(dg_frag%nstate_tot), one_state_sum(dg_frag%nstate_tot))
      allocate(nl_state_local(dg_frag%nstate_tot), nl_state_sum(dg_frag%nstate_tot))
      kin_state_local(:) = 0.0d0
      one_state_local(:) = 0.0d0
      nl_state_local(:) = 0.0d0
    end if

    do ispin = 1, dg_frag%nspin
      use_buffered_direct_orbitals = dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)
      if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
        nocc = dg_frag%nstate_tot
      else
        nocc = dg_frag%nocc_spin(ispin)
      end if
      if (nocc <= 0) cycle
      allocate(V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
      call build_total_potential_grid_local(mg, Vh, Vxc(ispin), Vpsl, V_total)
      allocate(psi_frag(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
      allocate(tpsi_frag(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
      allocate(hpsi_frag(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
      allocate(vnlpsi_frag(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
      allocate(coef_probe(dg_frag%n_mat_max, 1), nl_probe(dg_frag%n_mat_max, 1))

      if (use_buffered_direct_orbitals) then
        i_local = 0
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = i_local + 1
          nbf = min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag)
          if (nbf <= 0) cycle
          core_s(:) = dg_frag%ixyz_frag(:, ifrag)
          core_e(:) = dg_frag%ixyz_frag(:, ifrag) + dg_frag%nxyz_domain(:, ifrag) - 1
          ov_s(:) = max(core_s(:), mg%is(:))
          ov_e(:) = min(core_e(:), mg%ie(:))
          has_overlap = all(ov_s(:) <= ov_e(:))
          allocate(T_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
          allocate(H_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
          do jo = 1, nbf
            ig_j = dg_frag%index_basis(jo, ifrag, ispin)
            if (ig_j < 1 .or. ig_j > dg_frag%nstate_tot) cycle
            occ_probe = dg_frag%occ_state(ig_j, ispin)
            if (occ_probe <= 0.0d0) cycle
            psi_frag(:, :, :) = (0.0d0, 0.0d0)
            call build_hpsi_for_basis_probe(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi)

            do gz = mg%is(3), mg%ie(3)
              do gy = mg%is(2), mg%ie(2)
                do gx = mg%is(1), mg%ie(1)
                  ixg = gx
                  iyg = gy
                  izg = gz
                  call get_phi_value_at_global_probe(dg_frag, ifrag, i_local, jo, ixg, iyg, izg, phi_val)
                  if (phi_val == (0.0d0, 0.0d0)) cycle
                  psi_frag(ixg, iyg, izg) = phi_val
                end do
              end do
            end do

            if (has_overlap) then
              ztmp = (0.0d0, 0.0d0)
              do gz = ov_s(3), ov_e(3)
                do gy = ov_s(2), ov_e(2)
                  do gx = ov_s(1), ov_e(1)
                    ixg = gx
                    iyg = gy
                    izg = gz
                    call get_phi_value_at_global_probe(dg_frag, ifrag, i_local, jo, ixg, iyg, izg, phi_val)
                    if (phi_val == (0.0d0, 0.0d0)) cycle
                    ztmp = ztmp + conjg(phi_val) * T_phi(ixg, iyg, izg)
                  end do
                end do
              end do
              kin_local = kin_local + occ_probe * real(ztmp, kind=8) * system%hvol
              if (allocated(kin_state_local)) kin_state_local(ig_j) = kin_state_local(ig_j) + occ_probe * real(ztmp, kind=8) * system%hvol

              ztmp = (0.0d0, 0.0d0)
              do gz = ov_s(3), ov_e(3)
                do gy = ov_s(2), ov_e(2)
                  do gx = ov_s(1), ov_e(1)
                    ixg = gx
                    iyg = gy
                    izg = gz
                    call get_phi_value_at_global_probe(dg_frag, ifrag, i_local, jo, ixg, iyg, izg, phi_val)
                    if (phi_val == (0.0d0, 0.0d0)) cycle
                    ztmp = ztmp + conjg(phi_val) * H_phi(ixg, iyg, izg)
                  end do
                end do
              end do
              one_local = one_local + occ_probe * real(ztmp, kind=8) * system%hvol
              if (allocated(one_state_local)) one_state_local(ig_j) = one_state_local(ig_j) + occ_probe * real(ztmp, kind=8) * system%hvol
            end if

            call apply_nonlocal_pp_realspace_probe(dg_frag, mg, ppg, system, Ac_tot, psi_frag, vnlpsi_frag)
            if (has_overlap) then
              ztmp = (0.0d0, 0.0d0)
              do gz = ov_s(3), ov_e(3)
                do gy = ov_s(2), ov_e(2)
                  do gx = ov_s(1), ov_e(1)
                    ztmp = ztmp + conjg(psi_frag(gx, gy, gz)) * vnlpsi_frag(gx, gy, gz)
                  end do
                end do
              end do
              if (allocated(nl_state_local)) nl_state_local(ig_j) = nl_state_local(ig_j) + occ_probe * real(ztmp, kind=8) * system%hvol
            end if
          end do
          deallocate(T_phi, H_phi)
        end do
      else
        do io = 1, nocc
          if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
            occ_probe = dg_frag%occ_state(io, ispin)
          else
            occ_probe = system%rocc(io, 1, ispin)
          end if
          if (occ_probe <= 0.0d0) cycle
          i_local = 0
          do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
            i_local = i_local + 1
            nbf = min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag)
            if (nbf <= 0) cycle
            core_s(:) = dg_frag%ixyz_frag(:, ifrag)
            core_e(:) = dg_frag%ixyz_frag(:, ifrag) + dg_frag%nxyz_domain(:, ifrag) - 1
            ov_s(:) = max(core_s(:), mg%is(:))
            ov_e(:) = min(core_e(:), mg%ie(:))
            has_overlap = all(ov_s(:) <= ov_e(:))
            if (.not. has_overlap) cycle
            psi_frag(:, :, :) = (0.0d0, 0.0d0)
            tpsi_frag(:, :, :) = (0.0d0, 0.0d0)
            hpsi_frag(:, :, :) = (0.0d0, 0.0d0)
            allocate(T_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
            allocate(H_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
            do jo = 1, nbf
              ig_j = dg_frag%index_basis(jo, ifrag, ispin)
              if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle
              coeff = dg_frag%coef(ig_j, io, ispin)
              if (abs(coeff) == 0.0d0) cycle
              call build_hpsi_for_basis_probe(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi)
              do gz = ov_s(3), ov_e(3)
                do gy = ov_s(2), ov_e(2)
                  do gx = ov_s(1), ov_e(1)
                    ixg = gx
                    iyg = gy
                    izg = gz
                    call get_phi_value_at_global_probe(dg_frag, ifrag, i_local, jo, ixg, iyg, izg, phi_val)
                    if (phi_val == (0.0d0, 0.0d0)) cycle
                    psi_frag(ixg, iyg, izg) = psi_frag(ixg, iyg, izg) + coeff * phi_val
                    tpsi_frag(ixg, iyg, izg) = tpsi_frag(ixg, iyg, izg) + coeff * T_phi(ixg, iyg, izg)
                    hpsi_frag(ixg, iyg, izg) = hpsi_frag(ixg, iyg, izg) + coeff * H_phi(ixg, iyg, izg)
                  end do
                end do
              end do
            end do
            deallocate(T_phi, H_phi)
            ztmp = (0.0d0, 0.0d0)
            do gz = ov_s(3), ov_e(3)
              do gy = ov_s(2), ov_e(2)
                do gx = ov_s(1), ov_e(1)
                  ztmp = ztmp + conjg(psi_frag(gx, gy, gz)) * tpsi_frag(gx, gy, gz)
                end do
              end do
            end do
            kin_local = kin_local + occ_probe * real(ztmp, kind=8) * system%hvol
            if (allocated(kin_state_local)) kin_state_local(io) = kin_state_local(io) + occ_probe * real(ztmp, kind=8) * system%hvol

            ztmp = (0.0d0, 0.0d0)
            do gz = ov_s(3), ov_e(3)
              do gy = ov_s(2), ov_e(2)
                do gx = ov_s(1), ov_e(1)
                  ztmp = ztmp + conjg(psi_frag(gx, gy, gz)) * hpsi_frag(gx, gy, gz)
                end do
              end do
            end do
            one_local = one_local + occ_probe * real(ztmp, kind=8) * system%hvol
            if (allocated(one_state_local)) one_state_local(io) = one_state_local(io) + occ_probe * real(ztmp, kind=8) * system%hvol
          end do
        end do
      end if

      deallocate(psi_frag, tpsi_frag, hpsi_frag, vnlpsi_frag, V_total, coef_probe, nl_probe)
    end do

    if (dg_frag%icomm_frag /= COMM_GROUP_NULL) then
      call comm_summation(kin_local, kin_frag_sum, dg_frag%icomm_frag)
      call comm_summation(one_local, one_frag_sum, dg_frag%icomm_frag)
      if (.not. dg_frag%is_frag_root) then
        kin_frag_sum = 0.0d0
        one_frag_sum = 0.0d0
      end if
      call comm_summation(kin_frag_sum, kin_sum, dg_frag%icomm)
      call comm_summation(one_frag_sum, one_sum, dg_frag%icomm)
      frag_reduce_factor = 1.0d0
    else
      call comm_summation(kin_local, kin_sum, dg_frag%icomm)
      call comm_summation(one_local, one_sum, dg_frag%icomm)
      frag_reduce_factor = real(max(1, dg_frag%isize_frag), 8)
      kin_sum = kin_sum / frag_reduce_factor
      one_sum = one_sum / frag_reduce_factor
    end if
    if (allocated(kin_state_local)) then
      if (dg_frag%icomm_frag /= COMM_GROUP_NULL) then
        call comm_summation(kin_state_local, kin_state_sum, dg_frag%nstate_tot, dg_frag%icomm_frag)
        call comm_summation(one_state_local, one_state_sum, dg_frag%nstate_tot, dg_frag%icomm_frag)
        call comm_summation(nl_state_local, nl_state_sum, dg_frag%nstate_tot, dg_frag%icomm_frag)
        if (.not. dg_frag%is_frag_root) then
          kin_state_sum(:) = 0.0d0
          one_state_sum(:) = 0.0d0
          nl_state_sum(:) = 0.0d0
        end if
        call comm_summation(kin_state_sum, kin_state_local, dg_frag%nstate_tot, dg_frag%icomm)
        call comm_summation(one_state_sum, one_state_local, dg_frag%nstate_tot, dg_frag%icomm)
        call comm_summation(nl_state_sum, nl_state_local, dg_frag%nstate_tot, dg_frag%icomm)
        kin_state_sum(:) = kin_state_local(:)
        one_state_sum(:) = one_state_local(:)
        nl_state_sum(:) = nl_state_local(:)
      else
        call comm_summation(kin_state_local, kin_state_sum, dg_frag%nstate_tot, dg_frag%icomm)
        call comm_summation(one_state_local, one_state_sum, dg_frag%nstate_tot, dg_frag%icomm)
        call comm_summation(nl_state_local, nl_state_sum, dg_frag%nstate_tot, dg_frag%icomm)
        kin_state_sum(:) = kin_state_sum(:) / frag_reduce_factor
        one_state_sum(:) = one_state_sum(:) / frag_reduce_factor
        nl_state_sum(:) = nl_state_sum(:) / frag_reduce_factor
      end if
      nl_sum = sum(nl_state_sum)
      if (dg_frag%id == 0) then
        do io = 1, dg_frag%nstate_tot
          if (io <= size(dg_frag%occ_state, 1)) then
            if (dg_frag%occ_state(io, 1) <= 0.0d0) cycle
            write(*,'(1x,a,i0,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
              "        buffered-rs-energy-state: itt=", itt, " io=", io, &
              " occ=", dg_frag%occ_state(io, 1), " kin=", kin_state_sum(io), " one=", one_state_sum(io), &
              " nloc=", nl_state_sum(io)
          end if
        end do
        flush(6)
      end if
      deallocate(kin_state_local, kin_state_sum, one_state_local, one_state_sum, nl_state_local, nl_state_sum)
    end if
    kin_sum_out = kin_sum
    one_sum_out = one_sum
    nl_sum_out = nl_sum
  end subroutine compute_realspace_energy_probe

  subroutine compute_realspace_nonlocal_current_probe(dg_frag, system, mg, ppg, Ac_tot, current_nl_out)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    real(8),                intent(in)    :: Ac_tot(3)
    real(8),                intent(out)   :: current_nl_out(3)

    integer :: ispin, io, ifrag, i_local, jo, nbf, ig_j, nocc
    integer :: gx, gy, gz, ixg, iyg, izg
    logical :: use_buffered_direct_orbitals
    complex(8), allocatable :: psi_frag(:,:,:)
    complex(8) :: coeff, phi_val
    real(8) :: occ_probe
    real(8) :: current_local(3), current_state(3)

    current_nl_out(:) = 0.0d0
    if (dg_frag%use_plane_wave_basis) return
    if (.not. dg_frag%has_real_space_basis) return
    if (dg_frag%n_mat_max <= 0) return

    do ispin = 1, dg_frag%nspin
      use_buffered_direct_orbitals = dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)
      if (use_buffered_direct_orbitals) then
        nocc = dg_frag%nstate_tot
      else
        nocc = dg_frag%nocc_spin(ispin)
      end if
      if (nocc <= 0) cycle

      allocate(psi_frag(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))

      if (use_buffered_direct_orbitals) then
        i_local = 0
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = i_local + 1
          nbf = min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag)
          if (nbf <= 0) cycle
          do jo = 1, nbf
            ig_j = dg_frag%index_basis(jo, ifrag, ispin)
            if (ig_j < 1 .or. ig_j > dg_frag%nstate_tot) cycle
            occ_probe = dg_frag%occ_state(ig_j, ispin)
            if (occ_probe <= 0.0d0) cycle
            psi_frag(:, :, :) = (0.0d0, 0.0d0)
            do gz = mg%is(3), mg%ie(3)
              do gy = mg%is(2), mg%ie(2)
                do gx = mg%is(1), mg%ie(1)
                  ixg = gx
                  iyg = gy
                  izg = gz
                  call get_phi_value_at_global_probe(dg_frag, ifrag, i_local, jo, ixg, iyg, izg, phi_val)
                  if (phi_val == (0.0d0, 0.0d0)) cycle
                  psi_frag(ixg, iyg, izg) = phi_val
                end do
              end do
            end do
            call calc_current_nonlocal_realspace_probe(dg_frag, mg, ppg, system, Ac_tot, psi_frag, current_state)
            current_nl_out(:) = current_nl_out(:) + occ_probe * current_state(:)
          end do
        end do
      else
        do io = 1, nocc
          occ_probe = system%rocc(io, 1, ispin)
          if (occ_probe <= 0.0d0) cycle
          psi_frag(:, :, :) = (0.0d0, 0.0d0)
          i_local = 0
          do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
            i_local = i_local + 1
            nbf = min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag)
            if (nbf <= 0) cycle
            do jo = 1, nbf
              ig_j = dg_frag%index_basis(jo, ifrag, ispin)
              if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle
              coeff = dg_frag%coef(ig_j, io, ispin)
              if (abs(coeff) == 0.0d0) cycle
              do gz = mg%is(3), mg%ie(3)
                do gy = mg%is(2), mg%ie(2)
                  do gx = mg%is(1), mg%ie(1)
                    ixg = gx
                    iyg = gy
                    izg = gz
                    call get_phi_value_at_global_probe(dg_frag, ifrag, i_local, jo, ixg, iyg, izg, phi_val)
                    if (phi_val == (0.0d0, 0.0d0)) cycle
                    psi_frag(ixg, iyg, izg) = psi_frag(ixg, iyg, izg) + coeff * phi_val
                  end do
                end do
              end do
            end do
          end do
          call calc_current_nonlocal_realspace_probe(dg_frag, mg, ppg, system, Ac_tot, psi_frag, current_state)
          current_nl_out(:) = current_nl_out(:) + occ_probe * current_state(:)
        end do
      end if

      deallocate(psi_frag)
    end do
  end subroutine compute_realspace_nonlocal_current_probe

  subroutine calc_current_nonlocal_realspace_probe(dg_frag, mg, ppg, system, Ac_tot, psi_frag, jw_frag)
    use structures
    use communication, only: comm_summation
    use salmon_global, only: theory
    use math_constants, only: zi
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_rgrid),          intent(in) :: mg
    type(s_pp_grid),        intent(in) :: ppg
    type(s_dft_system),     intent(in) :: system
    real(8),                intent(in) :: Ac_tot(3)
    complex(8),             intent(in) :: psi_frag(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8),                intent(out):: jw_frag(3)

    integer :: ilocal, ilma, ia, j, ix, iy, iz
    real(8) :: xcoord, ycoord, zcoord, phase
    real(8) :: A_local(3)
    logical :: use_micro_A
    complex(8) :: proj_amp, uVpsi, phase_factor
    complex(8) :: uVpsi_r(3)
    complex(8), allocatable :: proj_local(:), proj_global(:)

    jw_frag(:) = 0.0d0
    if (ppg%Nlma <= 0 .or. .not. allocated(ppg%uV)) return

    allocate(proj_local(ppg%Nlma), proj_global(ppg%Nlma))
    proj_local(:) = (0.0d0, 0.0d0)
    proj_global(:) = (0.0d0, 0.0d0)
    use_micro_A = (trim(theory) == 'single_scale_maxwell_tddft' .and. allocated(system%Ac_micro%v))

    do ilma = 1, ppg%Nlma
      ia = ppg%ia_tbl(ilma)
      proj_amp = (0.0d0, 0.0d0)
      do j = 1, ppg%mps(ia)
        ix = ppg%jxyz(1, j, ia)
        iy = ppg%jxyz(2, j, ia)
        iz = ppg%jxyz(3, j, ia)
        if (ix < mg%is(1) .or. ix > mg%ie(1) .or. iy < mg%is(2) .or. iy > mg%ie(2) .or. iz < mg%is(3) .or. iz > mg%ie(3)) cycle
        xcoord = ppg%rxyz(1, j, ia)
        ycoord = ppg%rxyz(2, j, ia)
        zcoord = ppg%rxyz(3, j, ia)
        if (use_micro_A) then
          A_local(1:3) = system%Ac_micro%v(1:3, ix, iy, iz)
        else
          A_local(1:3) = Ac_tot(1:3)
        end if
        phase = A_local(1) * xcoord + A_local(2) * ycoord + A_local(3) * zcoord
        phase_factor = ppg%uV(j, ilma) * exp(-zi * phase)
        proj_amp = proj_amp + conjg(phase_factor) * psi_frag(ix, iy, iz)
      end do
      proj_local(ilma) = proj_amp
    end do

    call comm_summation(proj_local, proj_global, ppg%Nlma, dg_frag%icomm_frag)

    do ilocal = 1, ppg%ilocal_nlma
      ilma = ppg%ilocal_nlma2ilma(ilocal)
      ia = ppg%ilocal_nlma2ia(ilocal)
      uVpsi_r(:) = (0.0d0, 0.0d0)
      do j = 1, ppg%mps(ia)
        ix = ppg%jxyz(1, j, ia)
        iy = ppg%jxyz(2, j, ia)
        iz = ppg%jxyz(3, j, ia)
        if (ix < mg%is(1) .or. ix > mg%ie(1) .or. iy < mg%is(2) .or. iy > mg%ie(2) .or. iz < mg%is(3) .or. iz > mg%ie(3)) cycle
        xcoord = ppg%rxyz(1, j, ia)
        ycoord = ppg%rxyz(2, j, ia)
        zcoord = ppg%rxyz(3, j, ia)
        if (use_micro_A) then
          A_local(1:3) = system%Ac_micro%v(1:3, ix, iy, iz)
        else
          A_local(1:3) = Ac_tot(1:3)
        end if
        phase = A_local(1) * xcoord + A_local(2) * ycoord + A_local(3) * zcoord
        phase_factor = ppg%uV(j, ilma) * exp(-zi * phase)
        uVpsi_r(1) = uVpsi_r(1) + conjg(phase_factor) * xcoord * psi_frag(ix, iy, iz)
        uVpsi_r(2) = uVpsi_r(2) + conjg(phase_factor) * ycoord * psi_frag(ix, iy, iz)
        uVpsi_r(3) = uVpsi_r(3) + conjg(phase_factor) * zcoord * psi_frag(ix, iy, iz)
      end do
      uVpsi = ppg%rinv_uvu(ilma) * proj_global(ilma)
      jw_frag(:) = jw_frag(:) + aimag(conjg(uVpsi_r(:)) * uVpsi)
    end do

    jw_frag(:) = 2.0d0 * jw_frag(:)
    deallocate(proj_local, proj_global)
  end subroutine calc_current_nonlocal_realspace_probe

  subroutine apply_nonlocal_pp_realspace_probe(dg_frag, mg, ppg, system, Ac_tot, psi_frag, vnlpsi_frag)
    use structures
    use communication, only: comm_summation
    use salmon_global, only: theory
    use math_constants, only: zi
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_rgrid),          intent(in) :: mg
    type(s_pp_grid),        intent(in) :: ppg
    type(s_dft_system),     intent(in) :: system
    real(8),                intent(in) :: Ac_tot(3)
    complex(8),             intent(in) :: psi_frag(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    complex(8),             intent(out):: vnlpsi_frag(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))

    integer :: ilma, ia, j, ix, iy, iz
    real(8) :: xcoord, ycoord, zcoord, phase
    real(8) :: A_local(3)
    logical :: use_micro_A
    complex(8) :: proj_amp, proj_val, phase_factor
    complex(8), allocatable :: proj_local(:), proj_global(:)

    vnlpsi_frag(:, :, :) = (0.0d0, 0.0d0)
    if (ppg%Nlma <= 0 .or. .not. allocated(ppg%uV)) return

    allocate(proj_local(ppg%Nlma), proj_global(ppg%Nlma))
    proj_local(:) = (0.0d0, 0.0d0)
    proj_global(:) = (0.0d0, 0.0d0)
    use_micro_A = (trim(theory) == 'single_scale_maxwell_tddft' .and. allocated(system%Ac_micro%v))

    do ilma = 1, ppg%Nlma
      ia = ppg%ia_tbl(ilma)
      proj_amp = (0.0d0, 0.0d0)
      do j = 1, ppg%mps(ia)
        ix = ppg%jxyz(1, j, ia)
        iy = ppg%jxyz(2, j, ia)
        iz = ppg%jxyz(3, j, ia)
        if (ix < mg%is(1) .or. ix > mg%ie(1) .or. iy < mg%is(2) .or. iy > mg%ie(2) .or. iz < mg%is(3) .or. iz > mg%ie(3)) cycle
        xcoord = ppg%rxyz(1, j, ia)
        ycoord = ppg%rxyz(2, j, ia)
        zcoord = ppg%rxyz(3, j, ia)
        if (use_micro_A) then
          A_local(1:3) = system%Ac_micro%v(1:3, ix, iy, iz)
        else
          A_local(1:3) = Ac_tot(1:3)
        end if
        phase = A_local(1) * xcoord + A_local(2) * ycoord + A_local(3) * zcoord
        phase_factor = ppg%uV(j, ilma) * exp(-zi * phase)
        proj_amp = proj_amp + conjg(phase_factor) * psi_frag(ix, iy, iz)
      end do
      proj_local(ilma) = proj_amp
    end do

    call comm_summation(proj_local, proj_global, ppg%Nlma, dg_frag%icomm_frag)

    do ilma = 1, ppg%Nlma
      ia = ppg%ia_tbl(ilma)
      proj_val = ppg%rinv_uvu(ilma) * proj_global(ilma)
      do j = 1, ppg%mps(ia)
        ix = ppg%jxyz(1, j, ia)
        iy = ppg%jxyz(2, j, ia)
        iz = ppg%jxyz(3, j, ia)
        if (ix < mg%is(1) .or. ix > mg%ie(1) .or. iy < mg%is(2) .or. iy > mg%ie(2) .or. iz < mg%is(3) .or. iz > mg%ie(3)) cycle
        xcoord = ppg%rxyz(1, j, ia)
        ycoord = ppg%rxyz(2, j, ia)
        zcoord = ppg%rxyz(3, j, ia)
        if (use_micro_A) then
          A_local(1:3) = system%Ac_micro%v(1:3, ix, iy, iz)
        else
          A_local(1:3) = Ac_tot(1:3)
        end if
        phase = A_local(1) * xcoord + A_local(2) * ycoord + A_local(3) * zcoord
        phase_factor = ppg%uV(j, ilma) * exp(-zi * phase)
        vnlpsi_frag(ix, iy, iz) = vnlpsi_frag(ix, iy, iz) + proj_val * phase_factor
      end do
    end do

    deallocate(proj_local, proj_global)
  end subroutine apply_nonlocal_pp_realspace_probe
