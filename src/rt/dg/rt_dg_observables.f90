  subroutine calculate_observables(dg_frag, system, mg, stencil, ppg, rt, itt, Vh, Vxc, Vpsl)
    use salmon_global, only: theory
    use structures
    use communication, only: comm_summation
    use timer, only: timer_begin, timer_end, LOG_CALC_CURRENT
    use rt_dg_fragment_ops, only: apply_momentum_blocks, apply_matrix_blocks_batch, apply_nonlocal_pp_projector_batch, &
                    apply_mixed_hamiltonian, mixed_fp_coupling_active, copy_matrix_blocks_to_complex_dense, &
                    gather_full_coef_view, apply_overlap_operator, copy_overlap_operator_to_dense
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
    
    integer :: io, jo, ispin, idir, n, nocc, n_pw, n_tot, max_nocc
    integer :: active_state_cap
    integer :: ifrag, jfrag, ib, jb, i_idx, j_idx
    integer :: iblk, nrow_blk, ncol_blk, n_diag_block_ids, idb
    integer :: env_len, env_status, probe_stride, probe_nprint
    integer :: transition_stride
    character(len=64) :: env_val
    logical :: do_interface_check
    logical :: enable_orbital_probe
    logical :: enable_transition_probe
    logical :: enable_electron_probe
    real(8), allocatable :: interface_flow(:,:), dndt_frag(:)
    real(8) :: pair_residual, max_pair_residual, charge_balance_residual
    real(8) :: current_tmp, energy_tmp, pw_weight_local, kpw_dir
    real(8) :: current_io, energy_io
    real(8) :: occ_i
    real(8) :: Ac_tot(3), A_squared
    real(8) :: current_local(3), energy_local
    real(8) :: elec_coef_local, elec_plain_local
    real(8) :: elec_coef_sum, elec_plain_sum
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
    complex(8) :: minus_i
    complex(8), allocatable :: op_mat(:,:), tmp_mat(:,:), coef_all(:,:), tmp_all(:,:)
    complex(8), allocatable :: coef_frag_all(:,:), coef_pw_all(:,:), coef_frag_view(:,:), coef_pw_view(:,:)
    complex(8), allocatable :: tmp_probe(:,:), dense_probe_mat(:,:), dense_probe_out(:,:)
    complex(8), allocatable :: overlap_vec(:), overlap_dense(:,:)
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
    complex(8) :: current_blk_total, current_blk_diag, current_elem
    complex(8) :: mfp
    real(8), parameter :: unit_dir(3,3) = reshape((/ &
      1.0d0, 0.0d0, 0.0d0, &
      0.0d0, 1.0d0, 0.0d0, &
      0.0d0, 0.0d0, 1.0d0 /), (/3, 3/))
    ! Calculate local observables (only for assigned fragments)
    ! MPI aggregation will sum across all ranks
    current_local = 0.0d0
    energy_local = 0.0d0
    pw_weight_local = 0.0d0
    elec_coef_local = 0.0d0
    elec_plain_local = 0.0d0
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
    current_diag_local(:) = 0.0d0
    current_offdiag_local(:) = 0.0d0
    use_energy_components = enable_energy_component_probe .or. dg_frag%use_buffered_basis

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
      goto 1000
    end if
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n + n_pw
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
    transition_stride = 1
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
    if (enable_electron_probe) then
      allocate(overlap_vec(n_tot))
      if (n_pw == 0) allocate(overlap_dense(n, n))
    end if
    minus_i = cmplx(0.0d0, -1.0d0, kind=8)

    ! Current calculation via momentum operator matrix (velocity gauge)
    ! Following conventional RT implementation in density_matrix.f90:
    !   - momentum_mat stores <φ|∇|φ> (gradient operator)
    !   - Current: j = Im[<ψ|∇|ψ>] with factor 2 and normalization by ngrid
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
        if (allocated(dg_frag%occ_state)) then
          if (io <= size(dg_frag%occ_state, 1) .and. ispin <= size(dg_frag%occ_state, 2)) then
            occ_weight(io) = max(0.0d0, dg_frag%occ_state(io, ispin))
          end if
        else
          occ_weight(io) = system%rocc(io, 1, ispin)
        end if
      end do
      use_mixed_current = (n_pw > 0 .and. mixed_fp_coupling_active(dg_frag, ispin))
      call gather_full_coef_view(dg_frag, ispin, n, nocc, coef_frag_view, coef_pw_view, 1, nocc)
      coef_frag_all(1:n, 1:nocc) = coef_frag_view(1:n, 1:nocc)
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
            do ib = 1, n_tot
              current_io = current_io + aimag(conjg(coef_all(ib, io)) * tmp_all(ib, io))
            end do
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
          ! Factor -2.0: -1 for operator sign convention, 2 for Im[ψ*∇ψ] normalization
          current_tmp = 0.0d0
          do io = 1, nocc
            if (occ_weight(io) <= 0.0d0) cycle
            current_io = 0.0d0
            do ib = 1, n
              current_io = current_io + aimag(conjg(coef_frag_all(ib, io)) * tmp_mat(ib, io))
            end do
            current_tmp = current_tmp + occ_weight(io) * current_io
            if (enable_orbital_probe) current_orb_local((idir - 1) * max_nocc + io) = &
              current_orb_local((idir - 1) * max_nocc + io) - 2.0d0 * occ_weight(io) * current_io
          end do
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
              nrow_blk = dg_frag%n_basis(ifrag, ispin)
              ncol_blk = dg_frag%n_basis(jfrag, ispin)
              if (nrow_blk <= 0 .or. ncol_blk <= 0) cycle
              do io = 1, nocc
                occ_i = 1.0d0
                if (allocated(dg_frag%occ_state)) then
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
        current_local(idir) = current_local(idir) - 2.0d0 * current_tmp
      end do
    end do
    call timer_end(LOG_CALC_CURRENT)
    
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
      if (dg_frag%use_buffered_basis .and. allocated(dg_frag%occ_state)) then
        nocc = active_state_cap
      else
        nocc = min(dg_frag%nocc_spin(ispin), active_state_cap)
      end if
      if (nocc <= 0) cycle
      occ_weight(:) = 0.0d0
      do io = 1, nocc
        if (allocated(dg_frag%occ_state)) then
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
          else
            elec_coef_local = elec_coef_local + occ_weight(io) * &
              real(sum(conjg(coef_frag_all(1:n, io)) * matmul(overlap_dense(1:n, 1:n), coef_frag_all(1:n, io))))
            elec_plain_local = elec_plain_local + occ_weight(io) * &
              real(sum(conjg(coef_frag_all(1:n, io)) * coef_frag_all(1:n, io)))
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
    ! Cache retained for reuse

  1000 continue
    
    if (n_pw == 0) then
      call compute_realspace_energy_probe(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, Vh, Vxc, Vpsl, &
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
    if (enable_transition_probe) then
      call comm_summation(current_diag_local, current_diag_sum, 3, dg_frag%icomm)
      call comm_summation(current_offdiag_local, current_offdiag_sum, 3, dg_frag%icomm)
    end if
    call comm_summation(energy_local, dg_frag%total_energy, dg_frag%icomm)
    call comm_summation(pw_weight_local, dg_frag%pw_weight_raw, dg_frag%icomm)
    if (enable_electron_probe) then
      call comm_summation(elec_coef_local, elec_coef_sum, dg_frag%icomm)
      call comm_summation(elec_plain_local, elec_plain_sum, dg_frag%icomm)
    end if
    if (enable_orbital_probe) then
      call comm_summation(current_orb_local, current_orb_sum, 3 * max_nocc, dg_frag%icomm)
      call comm_summation(energy_orb_local, energy_orb_sum, max_nocc, dg_frag%icomm)
    end if
    frag_reduce_factor = real(max(1, dg_frag%isize_frag), 8)
    dg_frag%total_energy = dg_frag%total_energy / frag_reduce_factor
    if (enable_electron_probe) then
      elec_coef_sum = elec_coef_sum / frag_reduce_factor
      elec_plain_sum = elec_plain_sum / frag_reduce_factor
    end if
    if (enable_orbital_probe) energy_orb_sum(:) = energy_orb_sum(:) / frag_reduce_factor
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
      end if
    end if

    ! Current and PW weight are replicated over all ranks, so these remain world-averaged.
    dg_frag%current(:) = dg_frag%current(:) / real(max(1, dg_frag%isize), 8)
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

    if (n_pw == 0 .and. (itt == 1 .or. itt == 40)) then
      call debug_vloc_block_probe(dg_frag, system, mg, stencil, Vh, Vxc, Vpsl, itt)
    end if
    
    ! Normalize by global grid count exactly as conventional calc_current().
    ! This avoids decomposition-dependent scaling from local/grid-view differences.
    dg_frag%current(:) = dg_frag%current(:) / dble(system%ngrid)
    if (enable_orbital_probe) current_orb_sum(:) = current_orb_sum(:) / dble(system%ngrid)
    if (enable_transition_probe) then
      current_diag_sum(:) = current_diag_sum(:) / dble(system%ngrid)
      current_offdiag_sum(:) = current_offdiag_sum(:) / dble(system%ngrid)
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
          write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
            "        electron-probe(buffered): itt=", itt, " coeff_metric_S=", elec_coef_sum, " coeff_metric_2=", elec_plain_sum, &
            " Ne_raw=", dg_frag%elec_num_raw, " rho_scale=", dg_frag%rho_scale_factor
        else
          write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
            "        electron-probe: itt=", itt, " Ne_coef_S=", elec_coef_sum, " Ne_coef_2=", elec_plain_sum, &
            " Ne_raw=", dg_frag%elec_num_raw, " rho_scale=", dg_frag%rho_scale_factor
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
    occ_idx_report(:) = 0
    occ_val_report(:) = 0.0d0
    occ_sdiag_report(:) = 0.0d0
    occ_c2_report(:) = 0.0d0

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
        if (allocated(dg_frag%occ_state)) then
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
        write(*,'(1x,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
          "        electron-probe-t0 occ-state=", occ_idx_report(io), " occ=", occ_val_report(io), &
          " sdiag=", occ_sdiag_report(io), " c2=", occ_c2_report(io)
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

  subroutine compute_realspace_energy_probe(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, Vh, Vxc, Vpsl, energy_kin_mat, energy_one_mat, kin_sum_out, one_sum_out)
    use structures
    use communication, only: comm_summation, comm_is_root
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
    real(8),                intent(out)   :: kin_sum_out, one_sum_out

    integer :: ispin, io, ifrag, i_local, jo, nbf, ig_j, nocc
    integer :: core_s(3), core_e(3), ov_s(3), ov_e(3)
    integer :: gx, gy, gz, ixg, iyg, izg
    logical :: has_overlap
    real(8), allocatable :: V_total(:,:,:)
    complex(8), allocatable :: T_phi(:,:,:), H_phi(:,:,:)
    complex(8), allocatable :: psi_frag(:,:,:), tpsi_frag(:,:,:), hpsi_frag(:,:,:)
    complex(8), allocatable :: coef_probe(:,:), nl_probe(:,:)
    complex(8) :: coeff, ztmp, phi_val
    real(8) :: kin_local, one_local, kin_sum, one_sum
    real(8) :: frag_reduce_factor
    real(8) :: occ_probe
    logical :: use_buffered_direct_orbitals
    real(8), allocatable :: kin_state_local(:), kin_state_sum(:), one_state_local(:), one_state_sum(:)
    real(8), allocatable :: nl_state_local(:), nl_state_sum(:)

    kin_sum_out = 0.0d0
    one_sum_out = 0.0d0
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
          if (.not. has_overlap) cycle
          allocate(T_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
          allocate(H_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
          do jo = 1, nbf
            ig_j = dg_frag%index_basis(jo, ifrag, ispin)
            if (ig_j < 1 .or. ig_j > dg_frag%nstate_tot) cycle
            occ_probe = dg_frag%occ_state(ig_j, ispin)
            if (occ_probe <= 0.0d0) cycle
            call build_hpsi_for_basis_probe(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi)

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
          end do
          deallocate(T_phi, H_phi)
        end do
        if (allocated(nl_state_local)) then
          do ig_j = 1, dg_frag%nstate_tot
            occ_probe = dg_frag%occ_state(ig_j, ispin)
            if (occ_probe <= 0.0d0) cycle
            coef_probe(:, :) = (0.0d0, 0.0d0)
            nl_probe(:, :) = (0.0d0, 0.0d0)
            coef_probe(ig_j, 1) = (1.0d0, 0.0d0)
            call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_probe, nl_probe)
            nl_state_local(ig_j) = occ_probe * real(nl_probe(ig_j, 1), kind=8)
          end do
        end if
      else
        do io = 1, nocc
          if (allocated(dg_frag%occ_state)) then
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

      deallocate(psi_frag, tpsi_frag, hpsi_frag, V_total, coef_probe, nl_probe)
    end do

    call comm_summation(kin_local, kin_sum, dg_frag%icomm)
    call comm_summation(one_local, one_sum, dg_frag%icomm)
    frag_reduce_factor = real(max(1, dg_frag%isize_frag), 8)
    kin_sum = kin_sum / frag_reduce_factor
    one_sum = one_sum / frag_reduce_factor
    if (allocated(kin_state_local)) then
      call comm_summation(kin_state_local, kin_state_sum, dg_frag%nstate_tot, dg_frag%icomm)
      call comm_summation(one_state_local, one_state_sum, dg_frag%nstate_tot, dg_frag%icomm)
      call comm_summation(nl_state_local, nl_state_sum, dg_frag%nstate_tot, dg_frag%icomm)
      kin_state_sum(:) = kin_state_sum(:) / frag_reduce_factor
      one_state_sum(:) = one_state_sum(:) / frag_reduce_factor
      nl_state_sum(:) = nl_state_sum(:) / frag_reduce_factor
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
  end subroutine compute_realspace_energy_probe
