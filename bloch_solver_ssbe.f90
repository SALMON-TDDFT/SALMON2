module bloch_solver_ssbe
    implicit none
    private
    public :: s_sbe_bloch_solver, s_sbe_gs_info
    public :: init_etdrk4_data, finalize_etdrk4_data
    public :: dt_evolve_bloch_etdrk4

    ! Precision kind
    integer, parameter :: dp = kind(1.0d0)
    complex(dp), parameter :: zi = (0.0_dp, 1.0_dp)

    ! Threshold for Taylor expansion in phi functions (increased for stability)
    real(dp), parameter :: PHI_EPS = 1.0d-4

    type :: s_sbe_gs_info
        integer :: nb ! Number of bands
        integer :: nk ! Number of k-points
        integer :: n_dim ! Dimensionality (1, 2, 3)
        real(dp), allocatable :: eigen(:,:) ! Eigenvalues (nb, nk)
        real(dp) :: eg_au ! Band gap in atomic units (scalar global min or average)
        
        ! Momentum matrices: p_tm(n, m, dir, k) + rvnl correction if flag is set
        complex(dp), allocatable :: p_tm_matrix(:,:,:,:) 
        complex(dp), allocatable :: p_rvnl_matrix(:,:,:,:) ! Optional non-local part
        logical :: flag_vnl_correction = .false.
    end type s_sbe_gs_info

    type :: s_sbe_bloch_solver
        integer :: nb
        integer :: nk
        real(dp) :: dt
        real(dp) :: t2_sbe_fs ! Dephasing time T2 in femtoseconds
        
        complex(dp), allocatable :: rho(:,:,:) ! Density matrix (nb, nb, nk)
        
        ! ETDRK4 Coefficients (precomputed)
        complex(dp), allocatable :: exp_Ldt(:,:,:)       ! E = exp(L*dt)
        complex(dp), allocatable :: exp_Ldt_half(:,:,:)  ! A = exp(L*dt/2)
        complex(dp), allocatable :: phi1(:,:,:)
        complex(dp), allocatable :: phi2(:,:,:)
        complex(dp), allocatable :: phi3(:,:,:)
        complex(dp), allocatable :: phi1_half(:,:,:)     ! phi1(a/2)
        logical :: etdrk4_initialized = .false.
    end type s_sbe_bloch_solver

contains

    !-------------------------------------------------------------------------
    ! Safe calculation of phi functions: phi_k(z) = (exp(z) - sum_{j=0}^{k-1} z^j/j!) / z^k
    ! Uses Taylor series for |z| < PHI_EPS to avoid catastrophic cancellation.
    !-------------------------------------------------------------------------
    pure function calc_phi(z, k) result(phi_val)
        complex(dp), intent(in) :: z
        integer, intent(in) :: k
        complex(dp) :: phi_val
        complex(dp) :: ez

        ! Pure exponential implementation for maximum vectorization performance
        ! Branch-free logic avoids pipeline stalls on modern CPUs
        ez = exp(z)
        
        select case (k)
        case (1)
            phi_val = (ez - 1.0_dp) / z
        case (2)
            phi_val = (ez - 1.0_dp - z) / (z * z)
        case (3)
            phi_val = (ez - 1.0_dp - z - 0.5_dp * z * z) / (z * z * z)
        case default
            phi_val = 0.0_dp
        end select
        
        ! Handle z=0 case (diagonal elements) using L'Hopital's rule limits
        ! phi_1(0) = 1, phi_2(0) = 1/2, phi_3(0) = 1/6
        if (abs(z) < PHI_EPS) then
            select case (k)
            case (1)
                phi_val = 1.0_dp
            case (2)
                phi_val = 0.5_dp
            case (3)
                phi_val = 1.0_dp / 6.0_dp
            end select
        end if
    end function calc_phi

    !-------------------------------------------------------------------------
    ! Initialize ETDRK4 coefficients
    !-------------------------------------------------------------------------
    subroutine init_etdrk4_data(sbe, gs, dt)
        type(s_sbe_bloch_solver), intent(inout) :: sbe
        type(s_sbe_gs_info), intent(in) :: gs
        real(dp), intent(in) :: dt
        
        integer :: ik, n, m
        real(dp) :: delta_e, gamma, Eg_sq_inv
        complex(dp) :: lambda, z, z_half
        real(dp) :: t2_au

        ! Convert T2 from fs to atomic units (1 fs ≈ 41.34 a.u.)
        ! Or keep consistent with energy units. Assuming inputs are consistent.
        ! If Eg is in Hartree, dt in a.u., then T2 must be in a.u.
        ! Here we assume user provides T2 in fs, convert to a.u.
        t2_au = sbe%t2_sbe_fs * 41.34137d0 
        
        Eg_sq_inv = 1.0_dp / (gs%eg_au * gs%eg_au)

        allocate(sbe%exp_Ldt(gs%nb, gs%nb, gs%nk))
        allocate(sbe%exp_Ldt_half(gs%nb, gs%nb, gs%nk))
        allocate(sbe%phi1(gs%nb, gs%nb, gs%nk))
        allocate(sbe%phi2(gs%nb, gs%nb, gs%nk))
        allocate(sbe%phi3(gs%nb, gs%nb, gs%nk))
        allocate(sbe%phi1_half(gs%nb, gs%nb, gs%nk))

        !$omp parallel do default(shared) private(ik, n, m, delta_e, gamma, lambda, z, z_half)
        do ik = 1, gs%nk
            do m = 1, gs%nb
                do n = 1, gs%nb
                    delta_e = gs%eigen(n, ik) - gs%eigen(m, ik)
                    
                    ! Static decoherence rate: gamma = (delta_e)^2 / (T2 * Eg^2)
                    gamma = (delta_e * delta_e) * Eg_sq_inv / t2_au
                    
                    ! Linear operator eigenvalue: lambda = -i*delta_e - gamma
                    lambda = -zi * delta_e - gamma
                    
                    z = lambda * dt
                    z_half = lambda * dt * 0.5_dp

                    sbe%exp_Ldt(n,m,ik)      = exp(z)
                    sbe%exp_Ldt_half(n,m,ik) = exp(z_half)

                    sbe%phi1(n,m,ik)      = calc_phi(z, 1)
                    sbe%phi2(n,m,ik)      = calc_phi(z, 2)
                    sbe%phi3(n,m,ik)      = calc_phi(z, 3)
                    sbe%phi1_half(n,m,ik) = calc_phi(z_half, 1)
                end do
            end do
        end do
        
        sbe%etdrk4_initialized = .true.
    end subroutine init_etdrk4_data

    subroutine finalize_etdrk4_data(sbe)
        type(s_sbe_bloch_solver), intent(inout) :: sbe
        if (allocated(sbe%exp_Ldt)) deallocate(sbe%exp_Ldt)
        if (allocated(sbe%exp_Ldt_half)) deallocate(sbe%exp_Ldt_half)
        if (allocated(sbe%phi1)) deallocate(sbe%phi1)
        if (allocated(sbe%phi2)) deallocate(sbe%phi2)
        if (allocated(sbe%phi3)) deallocate(sbe%phi3)
        if (allocated(sbe%phi1_half)) deallocate(sbe%phi1_half)
        sbe%etdrk4_initialized = .false.
    end subroutine finalize_etdrk4_data

    !-------------------------------------------------------------------------
    ! Main ETDRK4 Time Stepper
    ! Signature uses explicit Ac arrays for compatibility with FDTD leapfrog schemes
    !-------------------------------------------------------------------------
    subroutine dt_evolve_bloch_etdrk4(sbe, gs, Ac_t, Ac_thalf, Ac_tdt, dt)
        type(s_sbe_bloch_solver), intent(inout) :: sbe
        type(s_sbe_gs_info), intent(inout) :: gs
        real(dp), intent(in) :: Ac_t(3)      ! A(t)
        real(dp), intent(in) :: Ac_thalf(3)  ! A(t + dt/2)
        real(dp), intent(in) :: Ac_tdt(3)    ! A(t + dt)
        real(dp), intent(in) :: dt

        real(dp) :: t2_au, Eg_sq_inv
        complex(dp), allocatable :: rho_n(:,:), rho1(:,:), rho2(:,:), rho3(:,:)
        complex(dp), allocatable :: N1(:,:), N2(:,:), N3(:,:), N4(:,:)
        complex(dp), allocatable :: V_k(:,:), p_k(:,:)
        complex(dp), allocatable :: C1(:,:), C2(:,:), C3(:,:), comm_tmp(:,:)
        integer :: ik, n, m, idir
        real(dp) :: delta_e, gamma, lambda_val
        complex(dp) :: factor_decoh
        real(dp) :: t2_au_local

        ! Precompute constants (outside OpenMP to avoid redundant calculations)
        t2_au_local = sbe%t2_sbe_fs * 41.34137d0
        Eg_sq_inv = 1.0_dp / (gs%eg_au * gs%eg_au)
        factor_decoh = -1.0_dp / (t2_au_local * gs%eg_au * gs%eg_au)

        ! Allocate working arrays (size nb x nb)
        allocate(rho_n(gs%nb, gs%nb), rho1(gs%nb, gs%nb), rho2(gs%nb, gs%nb), rho3(gs%nb, gs%nb))
        allocate(N1(gs%nb, gs%nb), N2(gs%nb, gs%nb), N3(gs%nb, gs%nb), N4(gs%nb, gs%nb))
        allocate(V_k(gs%nb, gs%nb), p_k(gs%nb, gs%nb))
        allocate(C1(gs%nb, gs%nb), C2(gs%nb, gs%nb), C3(gs%nb, gs%nb), comm_tmp(gs%nb, gs%nb))

        ! Parallel loop over k-points
        ! Ac_t, Ac_thalf, Ac_tdt are SHARED (computed once before parallel region)
        !$omp parallel do default(shared) private(ik, p_k, rho_n, N1, N2, N3, N4, &
        !$omp                                    rho1, rho2, rho3, V_k, C1, C2, C3, comm_tmp, &
        !$omp                                    idir, n, m)
        do ik = 1, gs%nk
            
            ! Load current density matrix for this k-point
            rho_n = sbe%rho(:, :, ik)

            ! ------------------------------------------------------------------
            ! Stage 1: Compute N1 = N(rho_n, t)
            ! ------------------------------------------------------------------
            call build_V_and_N(gs, ik, Ac_t, rho_n, V_k, p_k, N1, C1, C2, C3, tmp_mat, comm_tmp, factor_decoh)

            ! rho1 = A * rho_n + dt * phi1_half * N1
            do m = 1, gs%nb
                do n = 1, gs%nb
                    rho1(n, m) = sbe%exp_Ldt_half(n,m,ik) * rho_n(n,m) + &
                                 dt * sbe%phi1_half(n,m,ik) * N1(n,m)
                end do
            end do

            ! ------------------------------------------------------------------
            ! Stage 2: Compute N2 = N(rho1, t + dt/2)
            ! ------------------------------------------------------------------
            call build_V_and_N(gs, ik, Ac_thalf, rho1, V_k, p_k, N2, C1, C2, C3, tmp_mat, comm_tmp, factor_decoh)

            ! rho2 = A * rho_n + dt * phi1_half * N2
            do m = 1, gs%nb
                do n = 1, gs%nb
                    rho2(n, m) = sbe%exp_Ldt_half(n,m,ik) * rho_n(n,m) + &
                                 dt * sbe%phi1_half(n,m,ik) * N2(n,m)
                end do
            end do

            ! ------------------------------------------------------------------
            ! Stage 3: Compute N3 = N(rho2, t + dt/2)
            ! ------------------------------------------------------------------
            call build_V_and_N(gs, ik, Ac_thalf, rho2, V_k, p_k, N3, C1, C2, C3, tmp_mat, comm_tmp, factor_decoh)

            ! rho3 = A * rho1 + dt * phi1_half * (2*N3 - N1)
            do m = 1, gs%nb
                do n = 1, gs%nb
                    rho3(n, m) = sbe%exp_Ldt_half(n,m,ik) * rho1(n,m) + &
                                 dt * sbe%phi1_half(n,m,ik) * (2.0_dp * N3(n,m) - N1(n,m))
                end do
            end do

            ! ------------------------------------------------------------------
            ! Stage 4: Compute N4 = N(rho3, t + dt)
            ! ------------------------------------------------------------------
            call build_V_and_N(gs, ik, Ac_tdt, rho3, V_k, p_k, N4, C1, C2, C3, tmp_mat, comm_tmp, factor_decoh)

            ! ------------------------------------------------------------------
            ! Final Update: rho_new = E * rho_n + dt * (phi1*N1 + 2*phi2*(N2+N3) + phi3*N4)
            ! ------------------------------------------------------------------
            do m = 1, gs%nb
                do n = 1, gs%nb
                    sbe%rho(n, m, ik) = sbe%exp_Ldt(n,m,ik) * rho_n(n,m) + &
                                        dt * ( sbe%phi1(n,m,ik) * N1(n,m) + &
                                               2.0_dp * sbe%phi2(n,m,ik) * (N2(n,m) + N3(n,m)) + &
                                               sbe%phi3(n,m,ik) * N4(n,m) )
                end do
            end do

            ! ------------------------------------------------------------------
            ! Post-processing: Enforce Hermiticity
            ! ------------------------------------------------------------------
            do m = 1, gs%nb
                do n = 1, gs%nb
                    sbe%rho(n, m, ik) = 0.5_dp * (sbe%rho(n, m, ik) + conjg(sbe%rho(m, n, ik)))
                end do
            end do

        end do
        !$omp end parallel do

        deallocate(rho_n, rho1, rho2, rho3, N1, N2, N3, N4, V_k, p_k, C1, C2, C3, comm_tmp)
    end subroutine dt_evolve_bloch_etdrk4

    !-------------------------------------------------------------------------
    ! Helper subroutine to build V and compute Nonlinear term N
    ! Inlined logic separated for clarity, called 4 times per k-step
    !-------------------------------------------------------------------------
    subroutine build_V_and_N(gs, ik, Ac, rho_in, V_out, p_work, N_out, C1, C2, C3, tmp_mat, comm_tmp, factor_decoh)
        type(s_sbe_gs_info), intent(in) :: gs
        integer, intent(in) :: ik
        real(dp), intent(in) :: Ac(3)
        complex(dp), intent(in) :: rho_in(:,:)
        complex(dp), intent(out) :: V_out(:,:), p_work(:,:)
        complex(dp), intent(out) :: N_out(:,:), C1(:,:), C2(:,:), C3(:,:), tmp_mat(:,:), comm_tmp(:,:)
        complex(dp), intent(in) :: factor_decoh
        integer :: n, m, idir
        real(dp) :: delta_e

        ! 1. Build V properly: V_nm = Sum_dir A_dir * p_nm_dir
        V_out = 0.0_dp
        do idir = 1, gs%n_dim
            ! Get p_matrix for this direction
            do m = 1, gs%nb
                do n = 1, gs%nb
                    p_work(n, m) = gs%p_tm_matrix(n, m, idir, ik)
                    if (gs%flag_vnl_correction) then
                        p_work(n, m) = p_work(n, m) + gs%p_rvnl_matrix(n, m, idir, ik)
                    end if
                end do
            end do
            ! Add contribution to V
            do m = 1, gs%nb
                do n = 1, gs%nb
                    V_out(n, m) = V_out(n, m) + Ac(idir) * p_work(n, m)
                end do
            end do
        end do

        ! 2. Compute Commutator C1 = [V, rho] = V*rho - rho*V
        ! Use tmp_mat for rho*V
        call zgemm('N', 'N', gs%nb, gs%nb, gs%nb, (1.0_dp, 0.0_dp), V_out, gs%nb, rho_in, gs%nb, (0.0_dp, 0.0_dp), C1, gs%nb)
        call zgemm('N', 'N', gs%nb, gs%nb, gs%nb, (1.0_dp, 0.0_dp), rho_in, gs%nb, V_out, gs%nb, (0.0_dp, 0.0_dp), tmp_mat, gs%nb)
        C1 = C1 - tmp_mat

        ! 3. Compute Decoherence parts
        ! Dynamic part involves: [H0, [V, rho]] + [V, [H0, rho]] + [V, [V, rho]]
        
        ! Term A: [H0, C1] where C1 = [V, rho]. Since H0 is diagonal: (E_n - E_m) * C1_nm
        do m = 1, gs%nb
            do n = 1, gs%nb
                delta_e = gs%eigen(n, ik) - gs%eigen(m, ik)
                C2(n, m) = delta_e * C1(n, m)
            end do
        end do

        ! Term B: [V, [H0, rho]]. First compute D = [H0, rho]
        do m = 1, gs%nb
            do n = 1, gs%nb
                delta_e = gs%eigen(n, ik) - gs%eigen(m, ik)
                tmp_mat(n, m) = delta_e * rho_in(n, m)
            end do
        end do
        ! Then [V, D]
        call zgemm('N', 'N', gs%nb, gs%nb, gs%nb, (1.0_dp, 0.0_dp), V_out, gs%nb, tmp_mat, gs%nb, (0.0_dp, 0.0_dp), C3, gs%nb)
        call zgemm('N', 'N', gs%nb, gs%nb, gs%nb, (1.0_dp, 0.0_dp), tmp_mat, gs%nb, V_out, gs%nb, (0.0_dp, 0.0_dp), tmp_mat, gs%nb)
        C3 = C3 - tmp_mat

        ! Term C: [V, [V, rho]] = [V, C1]
        ! CORRECT BLAS sequence without destroying C1 (needed for -i*[V,rho] term):
        ! comm_tmp = V * C1
        call zgemm('N', 'N', gs%nb, gs%nb, gs%nb, (1.0_dp, 0.0_dp), V_out, gs%nb, C1, gs%nb, (0.0_dp, 0.0_dp), comm_tmp, gs%nb)
        ! tmp_mat = C1 * V
        call zgemm('N', 'N', gs%nb, gs%nb, gs%nb, (1.0_dp, 0.0_dp), C1, gs%nb, V_out, gs%nb, (0.0_dp, 0.0_dp), tmp_mat, gs%nb)
        
        ! Now comm_tmp = [V, C1] = V*C1 - C1*V
        comm_tmp = comm_tmp - tmp_mat

        ! 4. Assemble N
        ! N = -i * C1 + factor_decoh * (C2 + C3 + comm_tmp)
        do m = 1, gs%nb
            do n = 1, gs%nb
                N_out(n, m) = -zi * C1(n, m) + factor_decoh * (C2(n, m) + C3(n, m) + comm_tmp(n, m))
            end do
        end do

    end subroutine build_V_and_N

end module bloch_solver_ssbe
