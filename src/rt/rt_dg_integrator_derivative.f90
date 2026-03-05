  subroutine calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, dcoef_dt, dcoef_dt_pw)
    use structures
    use salmon_global, only: theory
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
    complex(8), allocatable :: H0c(:,:), M(:,:), dcoef_dt_h0(:,:), dcoef_dt_m(:,:)
    complex(8), allocatable :: coef_all(:,:), rhs_all(:,:)
    complex(8), allocatable :: rhs_in(:,:), rhs_eig(:,:), Uc(:,:), Uc_keep(:,:)
    complex(8) :: mfp
    real(8) :: huge_val
    logical :: found_nan
    integer :: nan_io, nan_jo
    real(8) :: max_abs_h0, max_abs_m
    complex(8), allocatable :: S_eval(:,:), work_s(:)
    real(8), allocatable :: eval_s(:), rwork_s(:)
    integer :: info_eig, lwork_s, n_s, n_floor, n_keep, n_drop, i_keep_start
    real(8) :: eps_s, tau_drop, cond_s, amp_ratio, rhs_in_norm, rhs_out_norm, s_max
    real(8), parameter :: eps_s_abs = 1.0d-12, eps_s_rel = 1.0d-8
    logical :: use_spatial_A
    logical :: disable_mfp
    character(len=32) :: env_mfp
    integer :: env_len, env_stat
    real(8), allocatable :: Ap_mat(:,:), A2_mat(:,:)
    
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
    huge_val = huge(0.0d0) / 2.0d0
    
    ! Calculate A^2 (diamagnetic term)
    A_squared = Ac_tot(1)**2 + Ac_tot(2)**2 + Ac_tot(3)**2
    if (A_squared /= A_squared) then
      write(*,'(a,i0,a,i0,a,3es12.4)') "[NaN] A_squared: rank=", dg_frag%id, " itt=", itt, " Ac_tot=", &
        Ac_tot(1), Ac_tot(2), Ac_tot(3)
      stop "NaN in A_squared"
    end if

    call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_tot)
    has_nonlocal = dg_frag%has_nl_cache
    
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
      allocate(Ap_mat(n, n), A2_mat(n, n))
    end if

    allocate(H0c(n_tot, n_tot), M(n_tot, n_tot))
    allocate(dcoef_dt_h0(n_tot, dg_frag%nstate_tot), dcoef_dt_m(n_tot, dg_frag%nstate_tot))
    allocate(coef_all(n_tot, dg_frag%nstate_tot), rhs_all(n_tot, dg_frag%nstate_tot))

    do ispin = 1, dg_frag%nspin
      ! Build H0c = H_0 + V_NL(A) + A^2/2
      H0c(:, :) = (0.0d0, 0.0d0)
      use_hmat_complex = allocated(dg_frag%H_mat_c) .and. allocated(dg_frag%phi_frag_c)
      if (use_hmat_complex) then
        H0c(1:n_frag, 1:n_frag) = dg_frag%H_mat_c(1:n_frag, 1:n_frag, ispin)
      else
        H0c(1:n_frag, 1:n_frag) = cmplx(dg_frag%H_mat(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8)
      end if
      if (n_pw > 0 .and. allocated(dg_frag%H_mat_mixed)) then
        H0c(1:n_tot, 1:n_tot) = dg_frag%H_mat_mixed(1:n_tot, 1:n_tot, ispin)
      end if

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
      
      if (has_nonlocal) then
        H0c(1:n_frag, 1:n_frag) = H0c(1:n_frag, 1:n_frag) + dg_frag%H_nl_cache(1:n_frag, 1:n_frag, ispin)
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
        do io = 1, n_tot
          H0c(io, io) = H0c(io, io) + 0.5d0 * A_squared
        end do

        ! Build M = A·<∇>
        M(:, :) = (0.0d0, 0.0d0)
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
        if (n_pw > 0) then
          do io = 1, n_pw
            M(n_frag+io, n_frag+io) = zi * dot_product(Ac_tot(1:3), dg_frag%k_pw(1:3, io))
          end do
          if (.not. disable_mfp .and. allocated(dg_frag%S_mat_frag_pw)) then
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
      if (any(real(M(:, :)) /= real(M(:, :))) .or. any(aimag(M(:, :)) /= aimag(M(:, :)))) then
        write(*,'(a,i0,a,i0,a,i0)') "[NaN] M: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
        stop "NaN in M"
      end if
      if (any(abs(M(:, :)) > huge_val)) then
        write(*,'(a,i0,a,i0,a,i0)') "[Inf] M: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
        stop "Inf in M"
      end if

      coef_all(:, :) = (0.0d0, 0.0d0)
      coef_all(1:n_frag, :) = dg_frag%coef(1:n_frag, 1:dg_frag%nstate_tot, ispin)
      if (n_pw > 0) coef_all(n_frag+1:n_tot, :) = dg_frag%coef_pw(1:n_pw, 1:dg_frag%nstate_tot, ispin)

      if (any(real(coef_all(:, :)) /= real(coef_all(:, :))) .or. &
          any(aimag(coef_all(:, :)) /= aimag(coef_all(:, :)))) then
        write(*,'(a,i0,a,i0,a,i0)') "[NaN] coef: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
        stop "NaN in coef before zgemm"
      end if
      if (any(abs(coef_all(:, :)) > huge_val)) then
        write(*,'(a,i0,a,i0,a,i0)') "[Inf] coef: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
        stop "Inf in coef before zgemm"
      end if

      if (any(abs(H0c(:, :)) > huge_val)) then
        write(*,'(a,i0,a,i0,a,i0)') "[Inf] H0c: rank=", dg_frag%id, " itt=", itt, " ispin=", ispin
        stop "Inf in H0c"
      end if

      ! dcoef_dt = -i * H0c * coef - M * coef
      call zgemm('N', 'N', n_tot, dg_frag%nstate_tot, n_tot, -zi, H0c, n_tot, &
                 coef_all, n_tot, (0.0d0, 0.0d0), dcoef_dt_h0, n_tot)
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

      call zgemm('N', 'N', n_tot, dg_frag%nstate_tot, n_tot, (1.0d0, 0.0d0), M, n_tot, &
                 coef_all, n_tot, (0.0d0, 0.0d0), dcoef_dt_m, n_tot)
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
      if (n_pw > 0 .and. allocated(dg_frag%S_mat_mixed_prop)) then
        n_s = n_tot
        allocate(S_eval(n_s, n_s))
        S_eval(:, :) = dg_frag%S_mat_mixed_prop(1:n_s, 1:n_s, ispin)
      else if (allocated(dg_frag%S_mat_prop_c)) then
        n_s = n_frag
        allocate(S_eval(n_s, n_s))
        S_eval(:, :) = dg_frag%S_mat_prop_c(1:n_s, 1:n_s, ispin)
      else if (allocated(dg_frag%S_mat_prop)) then
        n_s = n_frag
        allocate(S_eval(n_s, n_s))
        S_eval(:, :) = cmplx(dg_frag%S_mat_prop(1:n_s, 1:n_s, ispin), 0.0d0, kind=8)
      else if (allocated(dg_frag%S_mat_c)) then
        n_s = n_frag
        allocate(S_eval(n_s, n_s))
        S_eval(:, :) = dg_frag%S_mat_c(1:n_s, 1:n_s, ispin)
      else if (allocated(dg_frag%S_mat)) then
        n_s = n_frag
        allocate(S_eval(n_s, n_s))
        S_eval(:, :) = cmplx(dg_frag%S_mat(1:n_s, 1:n_s, ispin), 0.0d0, kind=8)
      end if
      if (n_s > 0) then
        allocate(eval_s(n_s))
        lwork_s = max(1, 2*n_s)
        allocate(work_s(lwork_s), rwork_s(max(1, 3*n_s-2)))
        call zheev('V', 'U', n_s, S_eval, n_s, eval_s, work_s, lwork_s, rwork_s, info_eig)
        if (info_eig == 0) then
          s_max = max(eval_s(n_s), 1.0d0)
          eps_s = max(eps_s_abs, eps_s_rel * s_max)
          tau_drop = max(eps_s_abs, 1.0d-6 * s_max)
          i_keep_start = 1
          do while (i_keep_start <= n_s .and. eval_s(i_keep_start) < tau_drop)
            i_keep_start = i_keep_start + 1
          end do
          n_keep = max(0, n_s - i_keep_start + 1)
          n_drop = n_s - n_keep
          n_floor = 0
          if (n_keep > 0) then
            cond_s = eval_s(n_s) / max(eval_s(i_keep_start), tau_drop)
          else
            cond_s = huge(1.0d0)
          end if

          allocate(rhs_in(n_s, dg_frag%nstate_tot), Uc(n_s, n_s))
          rhs_in(:, :) = rhs_all(1:n_s, :)
          Uc(:, :) = S_eval(:, :)
          rhs_all(1:n_s, :) = (0.0d0, 0.0d0)
          if (n_keep > 0) then
            allocate(Uc_keep(n_s, n_keep), rhs_eig(n_keep, dg_frag%nstate_tot))
            Uc_keep(:, :) = Uc(:, i_keep_start:n_s)
            call zgemm('C', 'N', n_keep, dg_frag%nstate_tot, n_s, (1.0d0,0.0d0), Uc_keep, n_s, rhs_in, n_s, (0.0d0,0.0d0), rhs_eig, n_keep)
            do io = 1, n_keep
              if (eval_s(i_keep_start + io - 1) < eps_s) n_floor = n_floor + 1
              rhs_eig(io, :) = rhs_eig(io, :) / max(eval_s(i_keep_start + io - 1), eps_s)
            end do
            call zgemm('N', 'N', n_s, dg_frag%nstate_tot, n_keep, (1.0d0,0.0d0), Uc_keep, n_s, rhs_eig, n_keep, (0.0d0,0.0d0), rhs_all(1:n_s,:), n_s)
            deallocate(Uc_keep, rhs_eig)
          end if

          rhs_in_norm = sqrt(sum(abs(rhs_in(:, :))**2))
          rhs_out_norm = sqrt(sum(abs(rhs_all(1:n_s, :))**2))
          amp_ratio = rhs_out_norm / max(rhs_in_norm, 1.0d-300)
          if (dg_frag%id == 0 .and. mod(itt, 10) == 0) then
            write(*,'(1x,a,i0,a,i0,a,1pe11.3,a,1pe11.3,a,i0,a,i0,a,1pe11.3,a,1pe11.3)') &
              "[SINV] itt=", itt, " spin=", ispin, " cond=", cond_s, " amp=", amp_ratio, &
              " n_drop=", n_drop, " n_floor=", n_floor, " eps=", eps_s, " tau_drop=", tau_drop
          end if
          deallocate(rhs_in, Uc)
        else
          if (dg_frag%id == 0) then
            write(*,'(1x,a,i0,a,i0,a,i0)') '[WARN] dsyev failed in S^{-1} application: info=', info_eig, &
              ', itt=', itt, ', ispin=', ispin
          end if
        end if
        deallocate(work_s, rwork_s, eval_s, S_eval)
      end if

      dcoef_dt(1:n_frag, 1:dg_frag%nstate_tot, ispin) = rhs_all(1:n_frag, :)
      if (present(dcoef_dt_pw) .and. n_pw > 0) then
        dcoef_dt_pw(1:n_pw, 1:dg_frag%nstate_tot, ispin) = rhs_all(n_frag+1:n_tot, :)
      end if
    end do

    if (allocated(Ap_mat)) deallocate(Ap_mat)
    if (allocated(A2_mat)) deallocate(A2_mat)
    deallocate(H0c, M, dcoef_dt_h0, dcoef_dt_m, coef_all, rhs_all)

    ! Cache retained for reuse
    
  end subroutine calculate_time_derivative
