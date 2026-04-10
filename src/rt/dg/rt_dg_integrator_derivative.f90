  subroutine calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, dcoef_dt, dcoef_dt_pw)
    use structures
    use salmon_global, only: theory
    use rt_dg_fragment_ops, only: apply_momentum_blocks, apply_matrix_blocks_batch, &
                                  apply_nonlocal_pp_projector_batch, &
                                  apply_nonlocal_pp_projector_batch_so, apply_mixed_hamiltonian, &
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
    logical :: has_nonlocal, has_so_nonlocal, use_hmat_complex
    logical :: need_h0_dense, need_m_dense
    complex(8), allocatable, save :: H0c(:,:), M(:,:), dcoef_dt_h0(:,:), dcoef_dt_m(:,:)
    complex(8), allocatable, save :: coef_all(:,:), rhs_all(:,:)
    complex(8), allocatable :: coef_frag_all(:,:), coef_pw_all(:,:)
    complex(8), allocatable :: coef_frag_other(:,:), coef_pw_other(:,:)
    complex(8), allocatable, save :: rhs_in(:,:), rhs_eig(:,:), Uc(:,:), Uc_keep(:,:)
    complex(8), allocatable, save :: coef_mix_all(:,:), rhs_mix(:,:), raw_rhs(:,:), op_mix(:,:), tmp_mix(:,:)
    complex(8), allocatable, save :: mfp_vec(:), coef_pw_scaled(:,:), tmp_fp_fw(:,:), tmp_fp_bw(:,:)
    complex(8) :: mfp
    real(8) :: huge_val
    real(8) :: t_coef_gather0, t_coef_gather1
    real(8) :: t_deriv0, t_deriv1, dt_deriv, dt_gather_local
    logical :: found_nan
    integer :: nan_io, nan_jo
    real(8) :: max_abs_h0, max_abs_m
    integer :: n_s, n_basis, ispin_other
    integer :: state_s, state_e, state0, nstate_blk
    logical :: use_spatial_A
    logical :: disable_mfp, use_mixed_basis
    character(len=32) :: env_mfp
    integer :: env_len, env_stat
    real(8), allocatable, save :: Ap_mat(:,:), A2_mat(:,:)
    integer, parameter :: state_block_size = 64
    
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
    call cpu_time(t_deriv0)
    dt_gather_local = 0.0d0
    huge_val = huge(0.0d0) / 2.0d0
    
    ! Calculate A^2 (diamagnetic term)
    A_squared = Ac_tot(1)**2 + Ac_tot(2)**2 + Ac_tot(3)**2
    if (A_squared /= A_squared) then
      write(*,'(a,i0,a,i0,a,3es12.4)') "[NaN] A_squared: rank=", dg_frag%id, " itt=", itt, " Ac_tot=", &
        Ac_tot(1), Ac_tot(2), Ac_tot(3)
      stop "NaN in A_squared"
    end if

    has_nonlocal = (ppg%Nlma > 0 .and. allocated(ppg%uV))
    has_so_nonlocal = (allocated(ppg%uv_so) .and. allocated(dg_frag%phi_frag_c) .and. dg_frag%nspin == 2)
    
    n = dg_frag%n_mat_max
    if (n <= 0) return
    n_frag = dg_frag%n_mat_max
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n_frag + n_pw

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
      if (.not. allocated(Ap_mat)) then
        allocate(Ap_mat(n, n), A2_mat(n, n))
      else if (size(Ap_mat, 1) /= n .or. size(Ap_mat, 2) /= n) then
        deallocate(Ap_mat, A2_mat)
        allocate(Ap_mat(n, n), A2_mat(n, n))
      end if
    end if

    if (.not. allocated(dcoef_dt_h0)) then
      allocate(dcoef_dt_h0(n_tot, dg_frag%nstate_tot), dcoef_dt_m(n_tot, dg_frag%nstate_tot))
    else if (size(dcoef_dt_h0, 1) /= n_tot .or. size(dcoef_dt_h0, 2) /= dg_frag%nstate_tot) then
      deallocate(dcoef_dt_h0, dcoef_dt_m)
      allocate(dcoef_dt_h0(n_tot, dg_frag%nstate_tot), dcoef_dt_m(n_tot, dg_frag%nstate_tot))
    end if
    if (.not. allocated(coef_all)) then
      allocate(coef_all(n_tot, dg_frag%nstate_tot), rhs_all(n_tot, dg_frag%nstate_tot))
    else if (size(coef_all, 1) /= n_tot .or. size(coef_all, 2) /= dg_frag%nstate_tot) then
      deallocate(coef_all, rhs_all)
      allocate(coef_all(n_tot, dg_frag%nstate_tot), rhs_all(n_tot, dg_frag%nstate_tot))
    end if

    do ispin = 1, dg_frag%nspin
      ! Build H0c = H_0 + V_NL(A) + A^2/2
      use_hmat_complex = allocated(dg_frag%H_mat_c) .and. allocated(dg_frag%phi_frag_c)
      use_mixed_basis = (n_pw > 0 .and. dg_frag%mixed_basis_ready .and. allocated(dg_frag%mixed_transform) .and. &
                         allocated(dg_frag%mixed_basis_dim) .and. allocated(dg_frag%coef_mix))
      n_basis = 0
      if (use_mixed_basis) n_basis = min(dg_frag%mixed_basis_dim(ispin), n_tot, size(dg_frag%coef_mix, 1))
      need_h0_dense = use_spatial_A .or. use_hmat_complex .or. (.not. allocated(dg_frag%H_mat_blocks)) .or. &
                      (n_pw > 0 .and. .not. allocated(dg_frag%H_mat_frag_pw)) .or. use_mixed_basis
      need_m_dense = use_spatial_A .or. (.not. allocated(dg_frag%momentum_blocks)) .or. use_mixed_basis
      if (need_h0_dense) then
        if (.not. allocated(H0c)) then
          allocate(H0c(n_tot, n_tot))
        else if (size(H0c, 1) /= n_tot .or. size(H0c, 2) /= n_tot) then
          deallocate(H0c)
          allocate(H0c(n_tot, n_tot))
        end if
      end if
      if (need_m_dense) then
        if (.not. allocated(M)) then
          allocate(M(n_tot, n_tot))
        else if (size(M, 1) /= n_tot .or. size(M, 2) /= n_tot) then
          deallocate(M)
          allocate(M(n_tot, n_tot))
        end if
      end if
      if (allocated(H0c)) H0c(:, :) = (0.0d0, 0.0d0)
      if (allocated(M)) M(:, :) = (0.0d0, 0.0d0)

      if (need_h0_dense .and. use_hmat_complex) then
        H0c(1:n_frag, 1:n_frag) = dg_frag%H_mat_c(1:n_frag, 1:n_frag, ispin)
      else if (need_h0_dense .and. allocated(dg_frag%H_mat_blocks)) then
        call copy_matrix_blocks_metric_to_complex_dense(dg_frag, dg_frag%H_mat_blocks, ispin, n_frag, H0c(1:n_frag, 1:n_frag))
      else if (need_h0_dense .and. .not. allocated(dg_frag%H_mat_blocks)) then
        H0c(1:n_frag, 1:n_frag) = cmplx(dg_frag%H_mat(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8)
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
!$omp parallel do schedule(static) private(jo)
            do jo = 1, n_pw
              M(1:n_frag, n_frag+jo) = mfp_vec(jo) * dg_frag%S_mat_frag_pw(1:n_frag, jo, ispin)
              M(n_frag+jo, 1:n_frag) = -conjg(M(1:n_frag, n_frag+jo))
            end do
!$omp end parallel do
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
      call cpu_time(t_coef_gather0)
      coef_all(:, :) = (0.0d0, 0.0d0)
      do state0 = 1, dg_frag%nstate_tot, state_block_size
        nstate_blk = min(state_block_size, dg_frag%nstate_tot - state0 + 1)
        state_s = state0
        state_e = state0 + nstate_blk - 1
        call gather_full_coef_view(dg_frag, ispin, n_frag, nstate_blk, coef_frag_all, coef_pw_all, state_s, state_e)
        coef_all(1:n_frag, state_s:state_e) = coef_frag_all(1:n_frag, 1:nstate_blk)
        if (n_pw > 0) coef_all(n_frag+1:n_tot, state_s:state_e) = coef_pw_all(1:n_pw, 1:nstate_blk)
      end do
      if (has_so_nonlocal) then
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
        do state0 = 1, dg_frag%nstate_tot, state_block_size
          nstate_blk = min(state_block_size, dg_frag%nstate_tot - state0 + 1)
          state_s = state0
          state_e = state0 + nstate_blk - 1
          call gather_full_coef_view(dg_frag, ispin_other, n_frag, nstate_blk, coef_frag_all, coef_pw_all, state_s, state_e)
          coef_frag_other(1:n_frag, state_s:state_e) = coef_frag_all(1:n_frag, 1:nstate_blk)
          if (n_pw > 0) coef_pw_other(1:n_pw, state_s:state_e) = coef_pw_all(1:n_pw, 1:nstate_blk)
        end do
      end if
      call cpu_time(t_coef_gather1)
      dt_gather_local = dt_gather_local + (t_coef_gather1 - t_coef_gather0)

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


      ! dcoef_dt = -i * H0c * coef - M * coef
      dcoef_dt_h0(:, :) = (0.0d0, 0.0d0)
      if (use_mixed_basis .and. n_basis > 0) then

        if (.not. allocated(rhs_mix)) then
          allocate(rhs_mix(n_basis, dg_frag%nstate_tot), raw_rhs(n_tot, dg_frag%nstate_tot), op_mix(n_basis, n_basis), &
                   tmp_mix(n_tot, n_basis))
        else if (size(rhs_mix, 1) /= n_basis .or. size(rhs_mix, 2) /= dg_frag%nstate_tot .or. &
                 size(raw_rhs, 1) /= n_tot .or. size(raw_rhs, 2) /= dg_frag%nstate_tot .or. &
                 size(op_mix, 1) /= n_basis .or. size(op_mix, 2) /= n_basis .or. &
                 size(tmp_mix, 1) /= n_tot .or. size(tmp_mix, 2) /= n_basis) then
          deallocate(rhs_mix, raw_rhs, op_mix, tmp_mix)
          allocate(rhs_mix(n_basis, dg_frag%nstate_tot), raw_rhs(n_tot, dg_frag%nstate_tot), op_mix(n_basis, n_basis), &
                   tmp_mix(n_tot, n_basis))
        end if
        rhs_mix(:, :) = (0.0d0, 0.0d0)
        raw_rhs(:, :) = (0.0d0, 0.0d0)
        call zgemm('N', 'N', n_tot, n_basis, n_tot, (1.0d0, 0.0d0), H0c, n_tot, &
                   dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin), n_tot, (0.0d0, 0.0d0), tmp_mix, n_tot)
        call zgemm('C', 'N', n_basis, n_basis, n_tot, (1.0d0, 0.0d0), &
                   dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin), n_tot, tmp_mix, n_tot, &
                   (0.0d0, 0.0d0), op_mix, n_basis)
        call zgemm('N', 'N', n_basis, dg_frag%nstate_tot, n_basis, -zi, op_mix, n_basis, &
                   coef_mix_all, n_basis, (0.0d0, 0.0d0), rhs_mix, n_basis)
        call zgemm('N', 'N', n_tot, dg_frag%nstate_tot, n_basis, (1.0d0, 0.0d0), &
                   dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin), n_tot, rhs_mix, n_basis, &
                   (0.0d0, 0.0d0), raw_rhs, n_tot)
        dcoef_dt_h0(:, :) = raw_rhs(:, :)
        if (has_nonlocal) then
          call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all(1:n_frag, :), dcoef_dt_h0(1:n_frag, :))
          if (has_so_nonlocal) then
            call apply_nonlocal_pp_projector_batch_so(dg_frag, mg, ppg, system, Ac_tot, ispin, &
                 coef_all(1:n_frag, :), coef_frag_other(1:n_frag, 1:dg_frag%nstate_tot), dcoef_dt_h0(1:n_frag, :))
          end if
        end if
        dcoef_dt_h0(1:n_tot, :) = dcoef_dt_h0(1:n_tot, :) + 0.5d0 * A_squared * coef_all(1:n_tot, :)
        dcoef_dt_h0(1:n_tot, :) = -zi * dcoef_dt_h0(1:n_tot, :)
      else if (n_pw > 0 .and. allocated(dg_frag%H_mat_frag_pw)) then

        call apply_mixed_hamiltonian(dg_frag, ispin, coef_all(1:n_tot, :), dcoef_dt_h0(1:n_tot, :))

        if (has_nonlocal) then
          call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all(1:n_frag, :), dcoef_dt_h0(1:n_frag, :))
          if (has_so_nonlocal) then
            call apply_nonlocal_pp_projector_batch_so(dg_frag, mg, ppg, system, Ac_tot, ispin, &
                 coef_all(1:n_frag, :), coef_frag_other(1:n_frag, 1:dg_frag%nstate_tot), dcoef_dt_h0(1:n_frag, :))
          end if
        end if
        dcoef_dt_h0(1:n_tot, :) = dcoef_dt_h0(1:n_tot, :) + 0.5d0 * A_squared * coef_all(1:n_tot, :)
        dcoef_dt_h0(1:n_tot, :) = -zi * dcoef_dt_h0(1:n_tot, :)
      else if (n_pw == 0 .and. .not. use_hmat_complex .and. allocated(dg_frag%H_mat_blocks)) then

        if (allocated(dg_frag%H_local_block_ids)) then
          call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, coef_all, dcoef_dt_h0, dg_frag%H_local_block_ids)
        else
          call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, coef_all, dcoef_dt_h0)
        end if

        if (has_nonlocal) then

          call apply_nonlocal_pp_projector_batch(dg_frag, mg, ppg, system, Ac_tot, ispin, coef_all, dcoef_dt_h0)
          if (has_so_nonlocal) then
            call apply_nonlocal_pp_projector_batch_so(dg_frag, mg, ppg, system, Ac_tot, ispin, &
                 coef_all(1:n_frag, :), coef_frag_other(1:n_frag, 1:dg_frag%nstate_tot), dcoef_dt_h0(1:n_frag, :))
          end if

        end if
        dcoef_dt_h0(1:n_frag, :) = dcoef_dt_h0(1:n_frag, :) + 0.5d0 * A_squared * coef_all(1:n_frag, :)
        dcoef_dt_h0(1:n_frag, :) = -zi * dcoef_dt_h0(1:n_frag, :)
      else

        call zgemm('N', 'N', n_tot, dg_frag%nstate_tot, n_tot, -zi, H0c, n_tot, &
                   coef_all, n_tot, (0.0d0, 0.0d0), dcoef_dt_h0, n_tot)

      end if


      if (use_mixed_basis .and. n_basis > 0) then
        rhs_mix(:, :) = (0.0d0, 0.0d0)
        call zgemm('N', 'N', n_tot, n_basis, n_tot, (1.0d0, 0.0d0), M, n_tot, &
                   dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin), n_tot, (0.0d0, 0.0d0), tmp_mix, n_tot)
        call zgemm('C', 'N', n_basis, n_basis, n_tot, (1.0d0, 0.0d0), &
                   dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin), n_tot, tmp_mix, n_tot, &
                   (0.0d0, 0.0d0), op_mix, n_basis)
        call zgemm('N', 'N', n_basis, dg_frag%nstate_tot, n_basis, (1.0d0, 0.0d0), op_mix, n_basis, &
                   coef_mix_all, n_basis, (0.0d0, 0.0d0), rhs_mix, n_basis)
        call zgemm('N', 'N', n_tot, dg_frag%nstate_tot, n_basis, (1.0d0, 0.0d0), &
                   dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin), n_tot, rhs_mix, n_basis, &
                   (0.0d0, 0.0d0), raw_rhs, n_tot)
        dcoef_dt_m(:, :) = raw_rhs(:, :)
      else if (allocated(dg_frag%momentum_blocks) .and. .not. use_spatial_A) then
        dcoef_dt_m(:, :) = (0.0d0, 0.0d0)
        call apply_momentum_blocks(dg_frag, ispin, Ac_tot, coef_all(1:n_frag, :), dcoef_dt_m(1:n_frag, :))
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
            if (.not. allocated(coef_pw_scaled)) then
              allocate(coef_pw_scaled(n_pw, dg_frag%nstate_tot), tmp_fp_bw(n_pw, dg_frag%nstate_tot))
            else if (size(coef_pw_scaled, 1) /= n_pw .or. size(coef_pw_scaled, 2) /= dg_frag%nstate_tot) then
              deallocate(coef_pw_scaled, tmp_fp_bw)
              allocate(coef_pw_scaled(n_pw, dg_frag%nstate_tot), tmp_fp_bw(n_pw, dg_frag%nstate_tot))
            end if
            if (.not. allocated(tmp_fp_fw)) then
              allocate(tmp_fp_fw(n_frag, dg_frag%nstate_tot))
            else if (size(tmp_fp_fw, 1) /= n_frag .or. size(tmp_fp_fw, 2) /= dg_frag%nstate_tot) then
              deallocate(tmp_fp_fw)
              allocate(tmp_fp_fw(n_frag, dg_frag%nstate_tot))
            end if

!$omp parallel do schedule(static) private(jo)
            do jo = 1, n_pw
              coef_pw_scaled(jo, :) = mfp_vec(jo) * coef_all(n_frag+jo, :)
            end do
!$omp end parallel do

            call zgemm('N', 'N', n_frag, dg_frag%nstate_tot, n_pw, (1.0d0, 0.0d0), &
                       dg_frag%S_mat_frag_pw(1:n_frag, 1:n_pw, ispin), n_frag, &
                       coef_pw_scaled, n_pw, (0.0d0, 0.0d0), tmp_fp_fw, n_frag)
            dcoef_dt_m(1:n_frag, :) = dcoef_dt_m(1:n_frag, :) + tmp_fp_fw(1:n_frag, :)

            call zgemm('C', 'N', n_pw, dg_frag%nstate_tot, n_frag, (1.0d0, 0.0d0), &
                       dg_frag%S_mat_frag_pw(1:n_frag, 1:n_pw, ispin), n_frag, &
                       coef_all, n_tot, (0.0d0, 0.0d0), tmp_fp_bw, n_pw)
!$omp parallel do schedule(static) private(jo)
            do jo = 1, n_pw
              dcoef_dt_m(n_frag+jo, :) = dcoef_dt_m(n_frag+jo, :) - conjg(mfp_vec(jo)) * tmp_fp_bw(jo, :)
            end do
!$omp end parallel do
          end if
        end if
      else
        call zgemm('N', 'N', n_tot, dg_frag%nstate_tot, n_tot, (1.0d0, 0.0d0), M, n_tot, &
                   coef_all, n_tot, (0.0d0, 0.0d0), dcoef_dt_m, n_tot)
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

        if (.not. allocated(rhs_in)) then
          allocate(rhs_in(n_s, dg_frag%nstate_tot))
        else if (size(rhs_in, 1) /= n_s .or. size(rhs_in, 2) /= dg_frag%nstate_tot) then
          deallocate(rhs_in)
          allocate(rhs_in(n_s, dg_frag%nstate_tot))
        end if
        rhs_in(:, :) = rhs_all(1:n_s, :)
        if (n_pw == 0 .and. allocated(dg_frag%H_local_rows) .and. size(dg_frag%H_local_rows) > 0) then
          call solve_overlap_operator_batch_local(dg_frag, ispin, rhs_in, rhs_all(1:n_s, :), .true.)
        else
          call solve_overlap_operator_batch(dg_frag, ispin, rhs_in, rhs_all(1:n_s, :), .true.)
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
    end do

    call cpu_time(t_deriv1)
    dt_deriv = t_deriv1 - t_deriv0

    ! Cache retained for reuse
    
  end subroutine calculate_time_derivative
