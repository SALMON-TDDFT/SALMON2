  subroutine calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, dcoef_dt, dcoef_dt_pw)
    use structures
    use salmon_global, only: theory
    use rt_dg_fragment_ops, only: apply_momentum_blocks, apply_matrix_blocks_batch, apply_nonlocal_pp_projector_batch, apply_mixed_hamiltonian, &
                                  solve_overlap_operator_batch, solve_overlap_operator_batch_local, mixed_fp_coupling_active, &
                                  copy_matrix_blocks_metric_to_complex_dense, copy_momentum_blocks_to_complex_dense, gather_full_coef_view
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
    
    integer :: io, jo, ispin, idir
    integer :: n, n_frag, n_pw, n_tot
    real(8) :: A_squared
    complex(8), parameter :: zi = (0.0d0, 1.0d0)  ! imaginary unit
    logical :: has_nonlocal, use_hmat_complex
    logical :: need_h0_dense, need_m_dense
    complex(8), allocatable :: H0c(:,:), M(:,:), dcoef_dt_h0(:,:), dcoef_dt_m(:,:)
    complex(8), allocatable :: coef_all(:,:), rhs_all(:,:)
    complex(8), allocatable :: coef_frag_all(:,:), coef_pw_all(:,:)
    complex(8), allocatable :: rhs_in(:,:), rhs_eig(:,:), Uc(:,:), Uc_keep(:,:)
    complex(8), allocatable :: coef_mix_all(:,:), rhs_mix(:,:), raw_rhs(:,:), op_mix(:,:)
    complex(8) :: mfp
    real(8) :: huge_val
    logical :: found_nan
    integer :: nan_io, nan_jo
    real(8) :: max_abs_h0, max_abs_m
    integer :: n_s, n_basis
    logical :: use_spatial_A
    logical :: disable_mfp, use_mixed_basis
    character(len=32) :: env_mfp
    integer :: env_len, env_stat
    real(8), allocatable :: Ap_mat(:,:), A2_mat(:,:)
    logical, parameter :: enable_derivative_trace = .false.
    logical, parameter :: enable_derivative_progress = .true.
    
    ! Time derivative in velocity gauge:
    !   d/dt c_i = -i * (H_0_ij + A^2(t)/2 * delta_ij) * c_j - A(t)·<i|∇|j> * c_j
    ! H(t) = H_0 - i*A(t)·∇ + A^2(t)/2
    ! The A^2/2 term is the diamagnetic contribution
    ! Complex coefficients are ESSENTIAL for:
    ! - Phase evolution: exp(-iE_n*t)
    ! - Superposition states
    ! - Oscillatory responses (optical, currents)
    if (enable_derivative_trace) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
        " stage=", "entry"
      flush(6)
    end if
    
    dcoef_dt = (0.0d0, 0.0d0)
    if (present(dcoef_dt_pw)) dcoef_dt_pw = (0.0d0, 0.0d0)
    huge_val = huge(0.0d0) / 2.0d0
    
    ! Calculate A^2 (diamagnetic term)
    A_squared = Ac_tot(1)**2 + Ac_tot(2)**2 + Ac_tot(3)**2
    if (A_squared /= A_squared) then
      write(*,'(a,i0,a,i0,a,3es12.4)') "[NaN] A_squared: rank=", dg_frag%id, " itt=", itt, " Ac_tot=", &
        Ac_tot(1), Ac_tot(2), Ac_tot(3)
      stop "NaN in A_squared"
    end if

    if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        derivative trace: stage=before-ensure-nonlocal"
      flush(6)
    end if
    if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        derivative trace: stage=after-ensure-nonlocal"
      flush(6)
    end if
    has_nonlocal = (ppg%Nlma > 0 .and. allocated(ppg%uV))
    
    n = dg_frag%n_mat_max
    if (n <= 0) return
    n_frag = dg_frag%n_mat_max
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n_frag + n_pw
    if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        derivative trace: stage=after-dims"
      flush(6)
    end if

    use_spatial_A = (trim(theory) == 'single_scale_maxwell_tddft' .and. allocated(system%Ac_micro%v) .and. dg_frag%has_real_space_basis)
    disable_mfp = .false.
    env_mfp = ''
    call get_environment_variable('SALMON_DG_DISABLE_MFP', env_mfp, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      if (env_mfp(1:1) == '1' .or. env_mfp(1:1) == 'y' .or. env_mfp(1:1) == 'Y' .or. &
          env_mfp(1:1) == 't' .or. env_mfp(1:1) == 'T') then
        disable_mfp = .true.
      end if
    end if
    if (use_spatial_A) then
      allocate(Ap_mat(n, n), A2_mat(n, n))
    end if

    allocate(dcoef_dt_h0(n_tot, dg_frag%nstate_tot), dcoef_dt_m(n_tot, dg_frag%nstate_tot))
    allocate(coef_all(n_tot, dg_frag%nstate_tot), rhs_all(n_tot, dg_frag%nstate_tot))
    if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        derivative trace: stage=after-main-alloc"
      flush(6)
      write(*,'(1x,a)') "        derivative trace: stage=before-spin-loop"
      flush(6)
    end if

    do ispin = 1, dg_frag%nspin
      if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0 .and. ispin == 1) then
        write(*,'(1x,a)') "        derivative trace: stage=spin1-entry"
        flush(6)
      end if
      if (enable_derivative_trace) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
          " ispin=", ispin, " stage=", "spin-entry"
        flush(6)
      end if
      ! Build H0c = H_0 + V_NL(A) + A^2/2
      use_hmat_complex = allocated(dg_frag%H_mat_c) .and. allocated(dg_frag%phi_frag_c)
      use_mixed_basis = (n_pw > 0 .and. dg_frag%mixed_basis_ready .and. allocated(dg_frag%mixed_transform) .and. &
                         allocated(dg_frag%mixed_basis_dim) .and. allocated(dg_frag%coef_mix))
      n_basis = 0
      if (use_mixed_basis) n_basis = min(dg_frag%mixed_basis_dim(ispin), n_tot, size(dg_frag%coef_mix, 1))
      need_h0_dense = use_spatial_A .or. use_hmat_complex .or. (.not. allocated(dg_frag%H_mat_blocks)) .or. &
                      (n_pw > 0 .and. .not. allocated(dg_frag%H_mat_frag_pw)) .or. use_mixed_basis
      need_m_dense = use_spatial_A .or. (.not. allocated(dg_frag%momentum_blocks)) .or. use_mixed_basis
      if (need_h0_dense .and. .not. allocated(H0c)) allocate(H0c(n_tot, n_tot))
      if (need_m_dense .and. .not. allocated(M)) allocate(M(n_tot, n_tot))
      if (allocated(H0c)) H0c(:, :) = (0.0d0, 0.0d0)
      if (allocated(M)) M(:, :) = (0.0d0, 0.0d0)

      if (need_h0_dense .and. use_hmat_complex) then
        H0c(1:n_frag, 1:n_frag) = dg_frag%H_mat_c(1:n_frag, 1:n_frag, ispin)
      else if (need_h0_dense .and. allocated(dg_frag%H_mat_blocks)) then
        call copy_matrix_blocks_metric_to_complex_dense(dg_frag, dg_frag%H_mat_blocks, ispin, n_frag, H0c(1:n_frag, 1:n_frag))
      else if (need_h0_dense .and. .not. allocated(dg_frag%H_mat_blocks)) then
        H0c(1:n_frag, 1:n_frag) = cmplx(dg_frag%H_mat(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8)
      end if
      if (allocated(dg_frag%H_mat)) then
        if (any(dg_frag%H_mat(1:n, 1:n, ispin) /= dg_frag%H_mat(1:n, 1:n, ispin))) then
          write(*,'(a,i0,a,i0,a,i0)') "[NaN] H_mat: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
          stop "NaN in H_mat"
        end if
        if (any(abs(dg_frag%H_mat(1:n, 1:n, ispin)) > huge_val)) then
          write(*,'(a,i0,a,i0,a,i0)') "[Inf] H_mat: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
          stop "Inf in H_mat"
        end if
        if (maxval(abs(dg_frag%H_mat(1:n, 1:n, ispin))) > 1.0d150) then
          write(*,'(a,i0,a,i0,a,i0,a,es12.4)') "[WARN] |H_mat| huge: rank=", dg_frag%id, " itt=", itt, &
            " ispin=", ispin, " max=", maxval(abs(dg_frag%H_mat(1:n, 1:n, ispin)))
        end if
      end if
      
      if (has_nonlocal .and. allocated(dg_frag%H_nl_cache)) then
        if (any(real(dg_frag%H_nl_cache(1:n, 1:n, ispin)) /= real(dg_frag%H_nl_cache(1:n, 1:n, ispin))) .or. &
            any(aimag(dg_frag%H_nl_cache(1:n, 1:n, ispin)) /= aimag(dg_frag%H_nl_cache(1:n, 1:n, ispin)))) then
          write(*,'(a,i0,a,i0,a,i0)') "[NaN] H_nl_cache: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
          stop "NaN in H_nl_cache"
        end if
        if (any(abs(dg_frag%H_nl_cache(1:n, 1:n, ispin)) > huge_val)) then
          write(*,'(a,i0,a,i0,a,i0)') "[Inf] H_nl_cache: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
          stop "Inf in H_nl_cache"
        end if
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
          do idir = 1, 3
            if (allocated(dg_frag%momentum_mat_c)) then
              if (any(abs(dg_frag%momentum_mat_c(idir, 1:n_frag, 1:n_frag, ispin)) > huge_val)) then
                write(*,'(a,i0,a,i0,a,i0,a,i0)') "[Inf] momentum_mat_c: rank=", dg_frag%id, " itt=", itt, &
                  " ispin=", ispin, " idir=", idir
                stop "Inf in momentum_mat_c"
              end if
              M(1:n_frag, 1:n_frag) = M(1:n_frag, 1:n_frag) + Ac_tot(idir) * dg_frag%momentum_mat_c(idir, 1:n_frag, 1:n_frag, ispin)
            else
              if (any(abs(dg_frag%momentum_mat(idir, 1:n_frag, 1:n_frag, ispin)) > huge_val)) then
                write(*,'(a,i0,a,i0,a,i0,a,i0)') "[Inf] momentum_mat: rank=", dg_frag%id, " itt=", itt, &
                  " ispin=", ispin, " idir=", idir
                stop "Inf in momentum_mat"
              end if
              M(1:n_frag, 1:n_frag) = M(1:n_frag, 1:n_frag) + Ac_tot(idir) * dg_frag%momentum_mat(idir, 1:n_frag, 1:n_frag, ispin)
            end if
          end do
        end if
        if (need_m_dense .and. n_pw > 0) then
          do io = 1, n_pw
            M(n_frag+io, n_frag+io) = zi * dot_product(Ac_tot(1:3), dg_frag%k_pw(1:3, io))
          end do
          if (.not. disable_mfp .and. mixed_fp_coupling_active(dg_frag, ispin)) then
            do jo = 1, n_pw
              do io = 1, n_frag
                mfp = zi * dot_product(Ac_tot(1:3), dg_frag%k_pw(1:3, jo)) * dg_frag%S_mat_frag_pw(io, jo, ispin)
                M(io, n_frag+jo) = mfp
                M(n_frag+jo, io) = -conjg(mfp)
              end do
            end do
          end if
        end if
      end if
      if (need_m_dense .and. (any(real(M(:, :)) /= real(M(:, :))) .or. any(aimag(M(:, :)) /= aimag(M(:, :))))) then
          write(*,'(a,i0,a,i0,a,i0)') "[NaN] M: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
          stop "NaN in M"
      end if
      if (need_m_dense .and. any(abs(M(:, :)) > huge_val)) then
        write(*,'(a,i0,a,i0,a,i0)') "[Inf] M: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
        stop "Inf in M"
      end if
      if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0 .and. ispin == 1) then
        write(*,'(1x,a)') "        derivative trace: stage=spin1-after-operators"
        flush(6)
      end if

      if (allocated(coef_frag_all)) deallocate(coef_frag_all)
      if (allocated(coef_pw_all)) deallocate(coef_pw_all)
      call gather_full_coef_view(dg_frag, ispin, n_frag, dg_frag%nstate_tot, coef_frag_all, coef_pw_all)
      coef_all(:, :) = (0.0d0, 0.0d0)
      coef_all(1:n_frag, :) = coef_frag_all(1:n_frag, 1:dg_frag%nstate_tot)
      if (n_pw > 0) coef_all(n_frag+1:n_tot, :) = coef_pw_all(1:n_pw, 1:dg_frag%nstate_tot)
      if (use_mixed_basis .and. n_basis > 0) then
        allocate(coef_mix_all(n_basis, dg_frag%nstate_tot))
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
      if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0 .and. ispin == 1) then
        write(*,'(1x,a)') "        derivative trace: stage=spin1-after-coef-pack"
        flush(6)
      end if

      if (need_h0_dense .and. any(abs(H0c(:, :)) > huge_val)) then
        write(*,'(a,i0,a,i0,a,i0)') "[Inf] H0c: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
        stop "Inf in H0c"
      end if
      if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0 .and. ispin == 1) then
        write(*,'(1x,a)') "        derivative trace: stage=spin1-before-h0-apply"
        flush(6)
      end if

      ! dcoef_dt = -i * H0c * coef - M * coef
      dcoef_dt_h0(:, :) = (0.0d0, 0.0d0)
      if (enable_derivative_trace) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
          " ispin=", ispin, " stage=", "before-h0"
        flush(6)
      end if
      if (use_mixed_basis .and. n_basis > 0) then
        if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0 .and. ispin == 1) then
          write(*,'(1x,a)') "        derivative trace: stage=spin1-h0-branch-mixed-basis"
          flush(6)
        end if
        allocate(rhs_mix(n_basis, dg_frag%nstate_tot), raw_rhs(n_tot, dg_frag%nstate_tot), op_mix(n_basis, n_basis))
        rhs_mix(:, :) = (0.0d0, 0.0d0)
        raw_rhs(:, :) = (0.0d0, 0.0d0)
        op_mix(:, :) = matmul(conjg(transpose(dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin))), &
                              matmul(H0c(1:n_tot, 1:n_tot), dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin)))
        call zgemm('N', 'N', n_basis, dg_frag%nstate_tot, n_basis, -zi, op_mix, n_basis, &
                   coef_mix_all, n_basis, (0.0d0, 0.0d0), rhs_mix, n_basis)
        raw_rhs(:, :) = matmul(dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin), rhs_mix)
        dcoef_dt_h0(:, :) = raw_rhs(:, :)
        if (has_nonlocal) then
          call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all(1:n_frag, :), dcoef_dt_h0(1:n_frag, :))
        end if
        do io = 1, n_tot
          dcoef_dt_h0(io, :) = dcoef_dt_h0(io, :) + 0.5d0 * A_squared * coef_all(io, :)
        end do
        dcoef_dt_h0(1:n_tot, :) = -zi * dcoef_dt_h0(1:n_tot, :)
      else if (n_pw > 0 .and. allocated(dg_frag%H_mat_frag_pw)) then
        if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0 .and. ispin == 1) then
          write(*,'(1x,a)') "        derivative trace: stage=spin1-h0-branch-mixed-h"
          flush(6)
        end if
        if (enable_derivative_trace) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
            " ispin=", ispin, " stage=", "before-apply-mixed-h"
          flush(6)
        end if
        call apply_mixed_hamiltonian(dg_frag, ispin, coef_all(1:n_tot, :), dcoef_dt_h0(1:n_tot, :))
        if (enable_derivative_trace) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
            " ispin=", ispin, " stage=", "after-apply-mixed-h"
          flush(6)
        end if
        if (has_nonlocal) then
          call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all(1:n_frag, :), dcoef_dt_h0(1:n_frag, :))
        end if
        do io = 1, n_tot
          dcoef_dt_h0(io, :) = dcoef_dt_h0(io, :) + 0.5d0 * A_squared * coef_all(io, :)
        end do
        dcoef_dt_h0(1:n_tot, :) = -zi * dcoef_dt_h0(1:n_tot, :)
      else if (n_pw == 0 .and. .not. use_hmat_complex .and. allocated(dg_frag%H_mat_blocks)) then
        if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0 .and. ispin == 1) then
          write(*,'(1x,a)') "        derivative trace: stage=spin1-h0-branch-block"
          flush(6)
          write(*,'(1x,a)') "        derivative trace: stage=spin1-before-apply-block-h"
          flush(6)
        end if
        if (enable_derivative_trace) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
            " ispin=", ispin, " stage=", "before-apply-block-h"
          flush(6)
        end if
        if (allocated(dg_frag%H_local_block_ids)) then
          call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, coef_all(1:n_frag, :), &
            dcoef_dt_h0(1:n_frag, :), dg_frag%H_local_block_ids)
        else
          call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, coef_all(1:n_frag, :), dcoef_dt_h0(1:n_frag, :))
        end if
        if (enable_derivative_trace) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
            " ispin=", ispin, " stage=", "after-apply-block-h"
          flush(6)
        end if
        if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0 .and. ispin == 1) then
          write(*,'(1x,a)') "        derivative trace: stage=spin1-after-apply-block-h"
          flush(6)
        end if
        if (has_nonlocal) then
          if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0 .and. ispin == 1) then
            write(*,'(1x,a)') "        derivative trace: stage=spin1-before-apply-block-nl"
            flush(6)
          end if
          call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all(1:n_frag, :), dcoef_dt_h0(1:n_frag, :))
          if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0 .and. ispin == 1) then
            write(*,'(1x,a)') "        derivative trace: stage=spin1-after-apply-block-nl"
            flush(6)
          end if
        end if
        do io = 1, n_frag
          dcoef_dt_h0(io, :) = dcoef_dt_h0(io, :) + 0.5d0 * A_squared * coef_all(io, :)
        end do
        dcoef_dt_h0(1:n_frag, :) = -zi * dcoef_dt_h0(1:n_frag, :)
      else
        if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0 .and. ispin == 1) then
          write(*,'(1x,a)') "        derivative trace: stage=spin1-h0-branch-zgemm"
          flush(6)
        end if
        if (enable_derivative_trace) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
            " ispin=", ispin, " stage=", "before-zgemm-h"
          flush(6)
        end if
        call zgemm('N', 'N', n_tot, dg_frag%nstate_tot, n_tot, -zi, H0c, n_tot, &
                   coef_all, n_tot, (0.0d0, 0.0d0), dcoef_dt_h0, n_tot)
        if (enable_derivative_trace) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
            " ispin=", ispin, " stage=", "after-zgemm-h"
          flush(6)
        end if
      end if
      if ((enable_derivative_trace .or. enable_derivative_progress) .and. dg_frag%id == 0 .and. ispin == 1) then
        write(*,'(1x,a)') "        derivative trace: stage=spin1-after-h0-apply"
        flush(6)
      end if
      if (enable_derivative_trace) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
          " ispin=", ispin, " stage=", "after-h0"
          flush(6)
      end if
      max_abs_h0 = maxval(abs(dcoef_dt_h0))
      if (max_abs_h0 > 1.0d150) then
        write(*,'(a,i0,a,i0,a,i0,a,es12.4)') "[WARN] |dcoef_dt_h0| huge: rank=", dg_frag%id, &
          " itt=", itt, " ispin=", ispin, " max=", max_abs_h0
      end if
      found_nan = .false.
      nan_io = 0
      nan_jo = 0
      do io = 1, dg_frag%nstate_tot
        do nan_jo = 1, n_tot
          if (real(dcoef_dt_h0(nan_jo, io)) /= real(dcoef_dt_h0(nan_jo, io)) .or. &
              aimag(dcoef_dt_h0(nan_jo, io)) /= aimag(dcoef_dt_h0(nan_jo, io))) then
            found_nan = .true.
            nan_io = io
            exit
          end if
        end do
        if (found_nan) exit
      end do
      if (found_nan) then
        write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') "[NaN] dcoef_dt_h0: rank=", dg_frag%id, " itt=", itt, &
          " ispin=", ispin, " io=", nan_io, " jo=", nan_jo
        stop "NaN in H0c term"
      end if

      if (enable_derivative_trace) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
          " ispin=", ispin, " stage=", "before-m"
        flush(6)
      end if
      if (use_mixed_basis .and. n_basis > 0) then
        rhs_mix(:, :) = (0.0d0, 0.0d0)
        op_mix(:, :) = matmul(conjg(transpose(dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin))), &
                              matmul(M(1:n_tot, 1:n_tot), dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin)))
        call zgemm('N', 'N', n_basis, dg_frag%nstate_tot, n_basis, (1.0d0, 0.0d0), op_mix, n_basis, &
                   coef_mix_all, n_basis, (0.0d0, 0.0d0), rhs_mix, n_basis)
        raw_rhs(:, :) = matmul(dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin), rhs_mix)
        dcoef_dt_m(:, :) = raw_rhs(:, :)
      else if (allocated(dg_frag%momentum_blocks) .and. .not. use_spatial_A) then
        dcoef_dt_m(:, :) = (0.0d0, 0.0d0)
        call apply_momentum_blocks(dg_frag, ispin, Ac_tot, coef_all(1:n_frag, :), dcoef_dt_m(1:n_frag, :))
        if (n_pw > 0) then
          do io = 1, n_pw
            mfp = zi * dot_product(Ac_tot(1:3), dg_frag%k_pw(1:3, io))
            dcoef_dt_m(n_frag+io, :) = dcoef_dt_m(n_frag+io, :) + mfp * coef_all(n_frag+io, :)
          end do
          if (.not. disable_mfp .and. mixed_fp_coupling_active(dg_frag, ispin)) then
            do jo = 1, n_pw
              mfp = zi * dot_product(Ac_tot(1:3), dg_frag%k_pw(1:3, jo))
              do io = 1, n_frag
                dcoef_dt_m(io, :) = dcoef_dt_m(io, :) + mfp * dg_frag%S_mat_frag_pw(io, jo, ispin) * coef_all(n_frag+jo, :)
                dcoef_dt_m(n_frag+jo, :) = dcoef_dt_m(n_frag+jo, :) - conjg(mfp * dg_frag%S_mat_frag_pw(io, jo, ispin)) * coef_all(io, :)
              end do
            end do
          end if
        end if
      else
        call zgemm('N', 'N', n_tot, dg_frag%nstate_tot, n_tot, (1.0d0, 0.0d0), M, n_tot, &
                   coef_all, n_tot, (0.0d0, 0.0d0), dcoef_dt_m, n_tot)
      end if
      if (enable_derivative_trace) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
          " ispin=", ispin, " stage=", "after-m"
        flush(6)
      end if
      max_abs_m = maxval(abs(dcoef_dt_m))
      if (max_abs_m > 1.0d150) then
        write(*,'(a,i0,a,i0,a,i0,a,es12.4)') "[WARN] |dcoef_dt_m| huge: rank=", dg_frag%id, &
          " itt=", itt, " ispin=", ispin, " max=", max_abs_m
      end if
      found_nan = .false.
      nan_io = 0
      nan_jo = 0
      do io = 1, dg_frag%nstate_tot
        do nan_jo = 1, n_tot
          if (real(dcoef_dt_m(nan_jo, io)) /= real(dcoef_dt_m(nan_jo, io)) .or. &
              aimag(dcoef_dt_m(nan_jo, io)) /= aimag(dcoef_dt_m(nan_jo, io))) then
            found_nan = .true.
            nan_io = io
            exit
          end if
        end do
        if (found_nan) exit
      end do
      if (found_nan) then
        write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') "[NaN] dcoef_dt_m: rank=", dg_frag%id, " itt=", itt, &
          " ispin=", ispin, " io=", nan_io, " jo=", nan_jo
        stop "NaN in M term"
      end if

      rhs_all = dcoef_dt_h0 - dcoef_dt_m
      n_s = 0
      if (use_mixed_basis .and. n_basis > 0) then
        n_s = 0
      else if (n_pw > 0 .and. mixed_fp_coupling_active(dg_frag, ispin)) then
        n_s = n_tot
      else if (allocated(dg_frag%S_mat_prop_blocks) .or. allocated(dg_frag%S_mat_prop_c) .or. &
               allocated(dg_frag%S_mat_prop) .or. allocated(dg_frag%S_mat_c) .or. allocated(dg_frag%S_mat)) then
        n_s = n_frag
      end if
      if (n_s > 0) then
        if (enable_derivative_trace) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
            " ispin=", ispin, " n_s=", n_s, " stage=", "before-overlap-solve"
          flush(6)
        end if
        allocate(rhs_in(n_s, dg_frag%nstate_tot))
        rhs_in(:, :) = rhs_all(1:n_s, :)
        if (n_pw == 0 .and. allocated(dg_frag%H_local_rows) .and. size(dg_frag%H_local_rows) > 0) then
          call solve_overlap_operator_batch_local(dg_frag, ispin, rhs_in, rhs_all(1:n_s, :), .true.)
        else
          call solve_overlap_operator_batch(dg_frag, ispin, rhs_in, rhs_all(1:n_s, :), .true.)
        end if
        deallocate(rhs_in)
        if (enable_derivative_trace) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
            " ispin=", ispin, " n_s=", n_s, " stage=", "after-overlap-solve"
          flush(6)
        end if
      end if

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
      if (enable_derivative_trace) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
          " ispin=", ispin, " stage=", "spin-exit"
        flush(6)
      end if
      if (allocated(coef_frag_all)) deallocate(coef_frag_all)
      if (allocated(coef_pw_all)) deallocate(coef_pw_all)
      if (allocated(coef_mix_all)) deallocate(coef_mix_all)
      if (allocated(rhs_mix)) deallocate(rhs_mix)
      if (allocated(raw_rhs)) deallocate(raw_rhs)
      if (allocated(op_mix)) deallocate(op_mix)
    end do

    if (allocated(Ap_mat)) deallocate(Ap_mat)
    if (allocated(A2_mat)) deallocate(A2_mat)
    if (allocated(H0c)) deallocate(H0c)
   if (allocated(M)) deallocate(M)
    deallocate(dcoef_dt_h0, dcoef_dt_m, coef_all, rhs_all)
    if (enable_derivative_trace) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        derivative trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
        " stage=", "exit"
      flush(6)
    end if

    ! Cache retained for reuse
    
  end subroutine calculate_time_derivative
