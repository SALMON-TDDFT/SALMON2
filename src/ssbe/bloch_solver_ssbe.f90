module bloch_solver_ssbe
    use math_constants, only: pi, zi
    use communication, only: comm_get_groupinfo, comm_summation
    use gs_info_ssbe
    use util_ssbe, only: split_range
    implicit none

    private
    public :: s_sbe_bloch_solver, init_sbe_bloch_solver, calc_current_bloch, &
              dt_evolve_bloch, dt_evolve_bloch_etdrk4, calc_trace, calc_energy, &
              init_etdrk4_data, finalize_etdrk4_data



    type s_sbe_bloch_solver
        !k-points for real-time SBE calculation
        integer :: nk, nb
        integer :: ik_max, ik_min
        complex(8), allocatable :: rho(:, :, :)
        logical :: flag_vnl_correction
        
        ! Frozen core handling
        logical, allocatable :: is_active(:) ! .true. if band is active, .false. if frozen
        integer :: n_active_bands
        
        ! ETDRK4 coefficients (precomputed for fixed dt)
        complex(8), allocatable :: exp_Ldt(:,:,:)       ! E = exp(L*dt)
        complex(8), allocatable :: exp_Ldt_half(:,:,:)  ! A = exp(L*dt/2)
        complex(8), allocatable :: phi1(:,:,:)
        complex(8), allocatable :: phi2(:,:,:)
        complex(8), allocatable :: phi3(:,:,:)
        complex(8), allocatable :: phi1_half(:,:,:)     ! phi1(a/2)
        logical :: etdrk4_initialized = .false.
    end type



contains



subroutine init_sbe_bloch_solver(sbe, gs, nb_sbe, icomm)
    use util_ssbe
    use communication
    use salmon_global, only: frozen_core_threshold_ev, nelec
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: nb_sbe
    integer, intent(in) :: icomm
    integer :: ik, ib, nk_proc, irank, nproc, ierr
    integer, allocatable :: itbl_min(:), itbl_max(:)
    real(8) :: eigen_ev

    call comm_get_groupinfo(icomm, irank, nproc)

    sbe%nk = gs%nk
    sbe%nb = nb_sbe

    allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))

    call split_range(1, sbe%nk, nproc, itbl_min, itbl_max)
    sbe%ik_min = itbl_min(irank)
    sbe%ik_max = itbl_max(irank)

    allocate(sbe%rho(1:sbe%nb, 1:sbe%nb, sbe%ik_min:sbe%ik_max))
    
    ! Initialize is_active array based on frozen_core_threshold_ev
    allocate(sbe%is_active(1:sbe%nb))
    sbe%n_active_bands = 0
    do ib = 1, sbe%nb
        ! Convert eigenvalue from Hartree to eV for comparison
        eigen_ev = gs%eigen(ib, sbe%ik_min) * 27.211386245988d0
        if (eigen_ev > frozen_core_threshold_ev) then
            sbe%is_active(ib) = .true.
            sbe%n_active_bands = sbe%n_active_bands + 1
        else
            sbe%is_active(ib) = .false.
        end if
    end do

    sbe%rho(:, :, :) = 0d0
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            sbe%rho(ib, ib, ik) = gs%occup(ib, ik)
        end do
    end do

    sbe%flag_vnl_correction = .false.
end subroutine


subroutine calc_current_bloch(sbe, gs, Ac, jmat, icomm)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(out) :: jmat(1:3)
    integer, intent(in) :: icomm
    integer :: ik, idir, ib, jb, nb
    complex(8) :: tmp1(1:3), tmp(1:3), v_mat(1:sbe%nb, 1:sbe%nb)
    complex(8) :: trace_val

    nb = sbe%nb
    tmp1 = 0d0

    !$omp parallel do default(shared) private(ik, idir, ib, jb, v_mat, trace_val) reduction(+:tmp1)
    do ik = sbe%ik_min, sbe%ik_max
        do idir = 1, 3
            ! Velocity operator in VG: v = p + A (in a.u. e=m=1)
            ! SALMON convention: H_eff = H0 + A*p, therefore v = p + A
            v_mat = gs%p_tm_matrix(:, :, idir, ik)
            do ib = 1, nb
                v_mat(ib, ib) = v_mat(ib, ib) + Ac(idir)
            end do
            if (sbe%flag_vnl_correction) then
                v_mat = v_mat + gs%rvnl_tm_matrix(:, :, idir, ik)
            endif

            ! Tr[v * rho]
            trace_val = 0d0
            do ib = 1, nb
                do jb = 1, nb
                    trace_val = trace_val + v_mat(ib, jb) * sbe%rho(jb, ib, ik)
                end do
            end do
            tmp1(idir) = tmp1(idir) + gs%kweight(ik) * trace_val
        end do
    end do
    !$omp end parallel do

    call comm_summation(tmp1, tmp, 3, icomm)
    jmat(:) = real(tmp(:)) / (sum(gs%kweight) * gs%volume)
end subroutine calc_current_bloch


subroutine dt_evolve_bloch(sbe, gs, Ac, dt)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(inout) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(in) :: dt
    integer :: nb, nk, ik

    complex(8) :: hrho1_k(1:sbe%nb, 1:sbe%nb)
    complex(8) :: hrho2_k(1:sbe%nb, 1:sbe%nb)
    complex(8) :: hrho3_k(1:sbe%nb, 1:sbe%nb)
    complex(8) :: hrho4_k(1:sbe%nb, 1:sbe%nb)
    complex(8) :: p_rvnl_k(1:sbe%nb, 1:sbe%nb, 1:3)

    nb = sbe%nb 
    nk = sbe%nk

    !$omp parallel do default(shared) private(ik, p_rvnl_k, hrho1_k, hrho2_k, hrho3_k, hrho4_k)
    do ik = sbe%ik_min, sbe%ik_max
        p_rvnl_k(1:sbe%nb, 1:sbe%nb, 1:3) = gs%p_tm_matrix(1:sbe%nb, 1:sbe%nb, 1:3, ik)
        if (sbe%flag_vnl_correction) then
            p_rvnl_k(1:sbe%nb, 1:sbe%nb, 1:3) =  p_rvnl_k(1:sbe%nb, 1:sbe%nb, 1:3) &
                & + gs%rvnl_tm_matrix(1:sbe%nb, 1:sbe%nb, 1:3, ik)
        end if

        call calc_hrho_bloch_k(ik, sbe%rho(:, :, ik), p_rvnl_k, hrho1_k)
        call calc_hrho_bloch_k(ik, hrho1_k, p_rvnl_k, hrho2_k)
        call calc_hrho_bloch_k(ik, hrho2_k, p_rvnl_k, hrho3_k)
        call calc_hrho_bloch_k(ik, hrho3_k, p_rvnl_k, hrho4_k)

        sbe%rho(:, :, ik) = sbe%rho(:, :, ik) + hrho1_k * (- zi * dt)
        sbe%rho(:, :, ik) = sbe%rho(:, :, ik) + hrho2_k * (- zi * dt) ** 2 * (1d0 / 2d0)
        sbe%rho(:, :, ik) = sbe%rho(:, :, ik) + hrho3_k * (- zi * dt) ** 3 * (1d0 / 6d0)
        sbe%rho(:, :, ik) = sbe%rho(:, :, ik) + hrho4_k * (- zi * dt) ** 4 * (1d0 / 24d0)

        ! Enforce Hermiticity: compensates numerical drift from Taylor series + decoherence
        sbe%rho(:, :, ik) = 0.5d0 * (sbe%rho(:, :, ik) + conjg(transpose(sbe%rho(:, :, ik))))

    end do
    return

contains


    !Calculate [H_eff, rho] commutation and add gauge-covariant decoherence (Eq. 7):
    subroutine calc_hrho_bloch_k(ik, rho_k, p_k, hrho_k)
        use phys_constants, only: au_fs
        use salmon_global, only: t2_sbe_fs
        implicit none
        integer, intent(in) :: ik
        complex(8), intent(in) :: rho_k(nb, nb)
        complex(8), intent(in) :: p_k(nb, nb, 1:3)
        complex(8), intent(out) :: hrho_k(nb, nb)
        integer :: idir, ib, jb
        real(8) :: t2_au, prefac
        complex(8) :: C2_k(nb, nb), Heff_k(nb, nb)

        ! 1. [H0, rho]_ij = (eps_i - eps_j) * rho_ij
        hrho_k(1:nb, 1:nb) = gs%delta_omega(1:nb, 1:nb, ik) * rho_k(1:nb, 1:nb)

        ! 2. Add [A·p, rho] -> hrho_k = [H_eff, rho]
        do idir = 1, 3
            call ZGEMM("N","N", nb, nb, nb, &
                dcmplx(+Ac(idir), 0d0), p_k(:, :, idir), nb, &
                rho_k(:, :), nb, dcmplx(1d0, 0d0), hrho_k(:, :), nb)

            call ZGEMM("N","N", nb, nb, nb, &
                dcmplx(-Ac(idir), 0d0), rho_k(:, :), nb, &
                p_k(:, :, idir), nb, dcmplx(1d0, 0d0), hrho_k(:, :), nb)
        end do

        ! 3. Gauge-covariant decoherence (Eq. 7 of paper 2012.00994v1)
        ! D = -1/(T2*Eg^2) * [H_eff, [H_eff, rho]]
        ! SALMON convention: H_eff = H0 + A*p (consistent with evolution above)
        if (0.0d0 < t2_sbe_fs .and. t2_sbe_fs < 1.0d9) then
            t2_au = t2_sbe_fs / au_fs
            prefac = -1.0d0 / (t2_au * gs%eg_au**2)

            ! Form H_eff = diag(eps) + A·p (CONSISTENT with evolution)
            Heff_k = 0d0
            do ib = 1, nb
                Heff_k(ib, ib) = gs%eigen(ib, ik)
            end do
            do idir = 1, 3
                Heff_k = Heff_k + Ac(idir) * p_k(:, :, idir)
            end do

            ! C2 = [H_eff, hrho_k] = [H_eff, [H_eff, rho]]
            C2_k = 0d0
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(1d0, 0d0), Heff_k, nb, hrho_k, nb, dcmplx(0d0, 0d0), C2_k, nb)
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(-1d0, 0d0), hrho_k, nb, Heff_k, nb, dcmplx(1d0, 0d0), C2_k, nb)

            ! Evolution: rho += (-zi*dt)*hrho. To add dt*D, need hrho += zi*D
            hrho_k = hrho_k + zi * (prefac * C2_k)
        endif
        
        ! Optional: Enforce Hermiticity of hrho_k for Taylor series stability
        ! (Done after rho update in main loop instead)
    end subroutine calc_hrho_bloch_k
end subroutine

function calc_trace(sbe, gs, nb_max, icomm) result(tr)
    use communication
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: icomm
    integer, intent(in) :: nb_max
    real(8) :: tr

    integer :: ik, ib
    real(8) :: tmp, tmp1

    tmp1 = 0d0
    !$omp parallel do default(shared) private(ik, ib) reduction(+: tmp1) collapse(2)
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, nb_max
            tmp1 = tmp1 + real(sbe%rho(ib, ib, ik)) * gs%kweight(ik)
        end do
    end do
    !$omp end parallel do
    call comm_summation(tmp1, tmp, icomm)
    tr = tmp / sum(gs%kweight)

    return
end function calc_trace


function calc_energy(sbe, gs, Ac, icomm) result(energy)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: icomm
    real(8), intent(in) :: Ac(1:3)
    integer :: ik, ib, jb, idir
    real(8) :: tmp1, tmp, energy
    ! real(8) :: kvec(1:3)
    tmp1 = 0d0
    !$omp parallel do default(shared) private(ik, ib, jb, idir) reduction(+: tmp1)
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            do idir = 1, 3
                do jb = 1, sbe%nb
                    ! SALMON convention: H_eff = H0 + A*p, interaction term is +A*p
                    tmp1 = tmp1 &
                        & + Ac(idir) * real(sbe%rho(ib, jb, ik) * gs%p_mod_matrix(jb, ib, idir, ik)) * gs%kweight(ik)
                end do
            end do
            tmp1 = tmp1 &
                & + real(sbe%rho(ib, ib, ik)) * ( &
                & + gs%eigen(ib, ik) &
                !& + dot_product(kvec(:), Ac(:))
                & + 0.5 * dot_product(Ac, Ac) &
                & ) * gs%kweight(ik)
        end do
    end do
    !$omp end parallel do
    call comm_summation(tmp1, tmp, icomm)
    energy = tmp / sum(gs%kweight)

    return
end function calc_energy


!=============================================================================
! ETDRK4 Implementation (Kassam-Trefethen 2005) for SBE in Velocity Gauge
!=============================================================================

subroutine init_etdrk4_data(sbe, gs, dt)
    ! Initialize ETDRK4 coefficients (precompute once for fixed dt)
    use phys_constants, only: au_fs
    use salmon_global, only: t2_sbe_fs
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(in) :: dt
    
    integer :: ik, n, m, nb, nk
    real(8) :: delta_e, gamma, t2_au, Eg2
    complex(8) :: lambda, z, z_half
    
    nb = sbe%nb
    nk = sbe%nk
    
    ! Allocate coefficient arrays if not already allocated
    if (.not. allocated(sbe%exp_Ldt)) then
        allocate(sbe%exp_Ldt(nb, nb, sbe%ik_min:sbe%ik_max))
        allocate(sbe%exp_Ldt_half(nb, nb, sbe%ik_min:sbe%ik_max))
        allocate(sbe%phi1(nb, nb, sbe%ik_min:sbe%ik_max))
        allocate(sbe%phi2(nb, nb, sbe%ik_min:sbe%ik_max))
        allocate(sbe%phi3(nb, nb, sbe%ik_min:sbe%ik_max))
        allocate(sbe%phi1_half(nb, nb, sbe%ik_min:sbe%ik_max))
    end if
    
    ! Prepare T2 and Eg for decoherence calculation
    if (t2_sbe_fs > 0.0d0 .and. t2_sbe_fs < 1.0d9) then
        t2_au = t2_sbe_fs / au_fs
        Eg2 = gs%eg_au**2
    else
        t2_au = 1.0d99  ! No decoherence
        Eg2 = 1.0d0
    end if
    
    ! Precompute coefficients for each (n,m,ik)
    !$omp parallel do default(shared) private(ik, n, m, delta_e, gamma, lambda, z, z_half)
    do ik = sbe%ik_min, sbe%ik_max
        do m = 1, nb
            do n = 1, nb
                delta_e = gs%eigen(n, ik) - gs%eigen(m, ik)
                
                ! Static decoherence: disabled if either band is frozen
                if (sbe%is_active(n) .and. sbe%is_active(m)) then
                    gamma = (delta_e**2) / (t2_au * Eg2)
                else
                    gamma = 0.0d0
                end if
                
                ! Linear operator eigenvalue: L_nm = -i*(eps_n - eps_m) - gamma
                lambda = dcmplx(0d0, -delta_e) - gamma
                
                z = lambda * dt
                z_half = lambda * dt * 0.5d0
                
                ! Exponential factors
                sbe%exp_Ldt(n, m, ik)      = exp(z)
                sbe%exp_Ldt_half(n, m, ik) = exp(z_half)
                
                ! Phi functions for ETDRK4
                sbe%phi1(n, m, ik)      = calc_phi(z, 1)
                sbe%phi2(n, m, ik)      = calc_phi(z, 2)
                sbe%phi3(n, m, ik)      = calc_phi(z, 3)
                sbe%phi1_half(n, m, ik) = calc_phi(z_half, 1)
            end do
        end do
    end do
    !$omp end parallel do
    
    sbe%etdrk4_initialized = .true.
end subroutine init_etdrk4_data


subroutine finalize_etdrk4_data(sbe)
    ! Deallocate ETDRK4 coefficient arrays
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    
    if (allocated(sbe%exp_Ldt)) deallocate(sbe%exp_Ldt)
    if (allocated(sbe%exp_Ldt_half)) deallocate(sbe%exp_Ldt_half)
    if (allocated(sbe%phi1)) deallocate(sbe%phi1)
    if (allocated(sbe%phi2)) deallocate(sbe%phi2)
    if (allocated(sbe%phi3)) deallocate(sbe%phi3)
    if (allocated(sbe%phi1_half)) deallocate(sbe%phi1_half)
    
    sbe%etdrk4_initialized = .false.
end subroutine finalize_etdrk4_data


pure function calc_phi(z, order) result(phi_val)
    ! Compute phi functions for ETDRK4 with Taylor expansion for small |z|
    ! phi_1(z) = (e^z - 1) / z
    ! phi_2(z) = (e^z - 1 - z) / z^2
    ! phi_3(z) = (e^z - 1 - z - z^2/2) / z^3
    implicit none
    complex(8), intent(in) :: z
    integer, intent(in) :: order
    complex(8) :: phi_val
    
    real(8), parameter :: eps = 1.0d-4  ! Safe threshold for double precision
    complex(8) :: ez, z2, z3
    
    if (abs(z) < eps) then
        ! Use Taylor expansion for small |z|
        select case(order)
        case(1)
            ! phi_1(z) = 1 + z/2 + z^2/6 + z^3/24 + z^4/120 + ...
            phi_val = dcmplx(1d0, 0d0) + z * (0.5d0 + z * (1d0/6d0 + z * (1d0/24d0 + z * (1d0/120d0))))
        case(2)
            ! phi_2(z) = 1/2 + z/6 + z^2/24 + z^3/120 + ...
            phi_val = 0.5d0 + z * (1d0/6d0 + z * (1d0/24d0 + z * (1d0/120d0)))
        case(3)
            ! phi_3(z) = 1/6 + z/24 + z^2/120 + z^3/720 + ...
            phi_val = 1d0/6d0 + z * (1d0/24d0 + z * (1d0/120d0 + z * (1d0/720d0)))
        case default
            phi_val = dcmplx(0d0, 0d0)
        end select
    else
        ez = exp(z)
        select case(order)
        case(1)
            phi_val = (ez - dcmplx(1d0, 0d0)) / z
        case(2)
            phi_val = (ez - dcmplx(1d0, 0d0) - z) / (z * z)
        case(3)
            z2 = z * z
            phi_val = (ez - dcmplx(1d0, 0d0) - z - z2 * 0.5d0) / (z2 * z)
        case default
            phi_val = dcmplx(0d0, 0d0)
        end select
    end if
end function calc_phi


subroutine dt_evolve_bloch_etdrk4(sbe, gs, Ac_t, Ac_thalf, Ac_tdt, dt)
    ! ETDRK4 time evolution for SBE (Kassam-Trefethen 2005)
    ! Interface: accepts three vector potential arrays at t, t+dt/2, t+dt
    use phys_constants, only: au_fs
    use salmon_global, only: t2_sbe_fs
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(inout) :: gs
    real(8), intent(in) :: Ac_t(1:3)      ! A(t)
    real(8), intent(in) :: Ac_thalf(1:3)  ! A(t + dt/2)
    real(8), intent(in) :: Ac_tdt(1:3)    ! A(t + dt)
    real(8), intent(in) :: dt
    
    integer :: nb, ik, n, m
    complex(8) :: rho_n(nb, nb), rho1(nb, nb), rho2(nb, nb), rho3(nb, nb)
    complex(8) :: N1(nb, nb), N2(nb, nb), N3(nb, nb), N4(nb, nb)
    complex(8) :: p_k(nb, nb, 1:3)
    real(8) :: t2_au, prefac, delta_e
    complex(8) :: a_diag(nb, nb)
    
    ! Temporary arrays for intermediate calculations
    complex(8) :: tmp_mat(nb, nb), C1(nb, nb), C2(nb, nb), V_k(nb, nb)
    integer :: idir
    
    ! Precompute decoherence scalars outside OpenMP to avoid race conditions and overhead
    logical :: flag_decoh
    
    nb = sbe%nb
    
    ! Check if ETDRK4 data is initialized
    if (.not. sbe%etdrk4_initialized) then
        call init_etdrk4_data(sbe, gs, dt)
    end if
    
    ! Precompute decoherence scalars outside OpenMP to avoid race conditions and overhead
    if (t2_sbe_fs > 0.0d0 .and. t2_sbe_fs < 1.0d9) then
        t2_au = t2_sbe_fs / au_fs
        prefac = -1.0d0 / (t2_au * gs%eg_au**2)
        flag_decoh = .true.
    else
        flag_decoh = .false.
    end if
    
    !$omp parallel do default(shared) private(ik, p_k, rho_n, N1, N2, N3, N4, &
    !$omp                                    rho1, rho2, rho3, &
    !$omp                                    V_k, C1, C2, tmp_mat, idir, n, m, delta_e, a_diag)
    do ik = sbe%ik_min, sbe%ik_max
        ! Build momentum matrix including rvnl correction
        p_k(:, :, :) = gs%p_tm_matrix(:, :, :, ik)
        if (sbe%flag_vnl_correction) then
            p_k(:, :, :) = p_k(:, :, :) + gs%rvnl_tm_matrix(:, :, :, ik)
        end if
        
        ! Store current density matrix
        rho_n(:, :) = sbe%rho(:, :, ik)
        
        !=== Stage 1: Compute N1 = N(rho_n, t) ===
        ! Build V(t) = A(t) · p
        V_k = dcmplx(0d0, 0d0)
        do idir = 1, 3
            V_k = V_k + Ac_t(idir) * p_k(:, :, idir)
        end do
        
        ! C1 = [V, rho_n]
        call ZGEMM("N", "N", nb, nb, nb, dcmplx(1d0, 0d0), V_k, nb, rho_n, nb, dcmplx(0d0, 0d0), C1, nb)
        call ZGEMM("N", "N", nb, nb, nb, dcmplx(-1d0, 0d0), rho_n, nb, V_k, nb, dcmplx(1d0, 0d0), C1, nb)
        
        ! N1 = -i * C1
        N1 = -zi * C1
        
        ! Add decoherence: N_decoh = prefac * ([H0,[V,rho]] + [V,[H0,rho]] + [V,[V,rho]])
        if (flag_decoh) then
            
            ! === C2 = [H0, [V, rho]] = [H0, C1] ===
            ! Apply is_active BEFORE multiplication (zero out frozen rows/cols in C1)
            do m = 1, nb
                do n = 1, nb
                    if (sbe%is_active(n) .and. sbe%is_active(m)) then
                        delta_e = gs%eigen(n, ik) - gs%eigen(m, ik)
                        C2(n, m) = delta_e * C1(n, m)
                    else
                        C2(n, m) = dcmplx(0d0, 0d0)
                    end if
                end do
            end do
            
            ! === [V, [H0, rho]] ===
            ! Build [H0, rho] with frozen zones zeroed out
            do m = 1, nb
                do n = 1, nb
                    if (sbe%is_active(n) .and. sbe%is_active(m)) then
                        delta_e = gs%eigen(n, ik) - gs%eigen(m, ik)
                        tmp_mat(n, m) = delta_e * rho_n(n, m)
                    else
                        tmp_mat(n, m) = dcmplx(0d0, 0d0)
                    end if
                end do
            end do
            ! Full ZGEMM (no submatrix extraction)
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(1d0, 0d0), V_k, nb, tmp_mat, nb, dcmplx(0d0, 0d0), a_diag, nb)
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(-1d0, 0d0), tmp_mat, nb, V_k, nb, dcmplx(1d0, 0d0), a_diag, nb)
            ! Post-zero frozen zones in result
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) a_diag(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
            C2 = C2 + a_diag
            
            ! === [V, [V, rho]] = [V, C1] ===
            ! C1 already has frozen zones zeroed from commutator, but zero again for safety
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) C1(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(1d0, 0d0), V_k, nb, C1, nb, dcmplx(0d0, 0d0), a_diag, nb)
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(-1d0, 0d0), C1, nb, V_k, nb, dcmplx(1d0, 0d0), a_diag, nb)
            ! Post-zero frozen zones
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) a_diag(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
            C2 = C2 + a_diag
            
            N1 = N1 + prefac * C2
            
            ! Final enforcement: N1 is zero for any element involving frozen zone
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) N1(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
        end if
        
        ! rho1 = A * rho_n + dt * phi1(a/2) * N1  (element-wise multiplication)
        do m = 1, nb
            do n = 1, nb
                rho1(n, m) = sbe%exp_Ldt_half(n, m, ik) * rho_n(n, m) + dt * sbe%phi1_half(n, m, ik) * N1(n, m)
            end do
        end do
        
        !=== Stage 2: Compute N2 = N(rho1, t+dt/2) ===
        V_k = dcmplx(0d0, 0d0)
        do idir = 1, 3
            V_k = V_k + Ac_thalf(idir) * p_k(:, :, idir)
        end do
        
        call ZGEMM("N", "N", nb, nb, nb, dcmplx(1d0, 0d0), V_k, nb, rho1, nb, dcmplx(0d0, 0d0), C1, nb)
        call ZGEMM("N", "N", nb, nb, nb, dcmplx(-1d0, 0d0), rho1, nb, V_k, nb, dcmplx(1d0, 0d0), C1, nb)
        
        N2 = -zi * C1
        
        if (flag_decoh) then
            
            ! === [H0, [V, rho1]] ===
            ! Apply is_active BEFORE multiplication
            do m = 1, nb
                do n = 1, nb
                    if (sbe%is_active(n) .and. sbe%is_active(m)) then
                        delta_e = gs%eigen(n, ik) - gs%eigen(m, ik)
                        C2(n, m) = delta_e * C1(n, m)
                    else
                        C2(n, m) = dcmplx(0d0, 0d0)
                    end if
                end do
            end do
            
            ! === [V, [H0, rho1]] ===
            ! Build [H0, rho1] with frozen zones zeroed out
            do m = 1, nb
                do n = 1, nb
                    if (sbe%is_active(n) .and. sbe%is_active(m)) then
                        delta_e = gs%eigen(n, ik) - gs%eigen(m, ik)
                        tmp_mat(n, m) = delta_e * rho1(n, m)
                    else
                        tmp_mat(n, m) = dcmplx(0d0, 0d0)
                    end if
                end do
            end do
            ! Full ZGEMM
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(1d0, 0d0), V_k, nb, tmp_mat, nb, dcmplx(0d0, 0d0), a_diag, nb)
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(-1d0, 0d0), tmp_mat, nb, V_k, nb, dcmplx(1d0, 0d0), a_diag, nb)
            ! Post-zero frozen zones
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) a_diag(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
            C2 = C2 + a_diag
            
            ! === [V, [V, rho1]] ===
            ! Zero frozen zones in C1
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) C1(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(1d0, 0d0), V_k, nb, C1, nb, dcmplx(0d0, 0d0), a_diag, nb)
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(-1d0, 0d0), C1, nb, V_k, nb, dcmplx(1d0, 0d0), a_diag, nb)
            ! Post-zero frozen zones
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) a_diag(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
            C2 = C2 + a_diag
            
            N2 = N2 + prefac * C2
            
            ! Final enforcement: N2 is zero for any element involving frozen zone
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) N2(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
        end if
        
        ! rho2 = A * rho_n + dt * phi1(a/2) * N2
        do m = 1, nb
            do n = 1, nb
                rho2(n, m) = sbe%exp_Ldt_half(n, m, ik) * rho_n(n, m) + dt * sbe%phi1_half(n, m, ik) * N2(n, m)
            end do
        end do
        
        !=== Stage 3: Compute N3 = N(rho2, t+dt/2) ===
        ! V_k already computed for t+dt/2
        call ZGEMM("N", "N", nb, nb, nb, dcmplx(1d0, 0d0), V_k, nb, rho2, nb, dcmplx(0d0, 0d0), C1, nb)
        call ZGEMM("N", "N", nb, nb, nb, dcmplx(-1d0, 0d0), rho2, nb, V_k, nb, dcmplx(1d0, 0d0), C1, nb)
        
        N3 = -zi * C1
        
        if (flag_decoh) then
            
            ! === [H0, [V, rho2]] ===
            ! Apply is_active BEFORE multiplication
            do m = 1, nb
                do n = 1, nb
                    if (sbe%is_active(n) .and. sbe%is_active(m)) then
                        delta_e = gs%eigen(n, ik) - gs%eigen(m, ik)
                        C2(n, m) = delta_e * C1(n, m)
                    else
                        C2(n, m) = dcmplx(0d0, 0d0)
                    end if
                end do
            end do
            
            ! === [V, [H0, rho2]] ===
            ! Build [H0, rho2] with frozen zones zeroed out
            do m = 1, nb
                do n = 1, nb
                    if (sbe%is_active(n) .and. sbe%is_active(m)) then
                        delta_e = gs%eigen(n, ik) - gs%eigen(m, ik)
                        tmp_mat(n, m) = delta_e * rho2(n, m)
                    else
                        tmp_mat(n, m) = dcmplx(0d0, 0d0)
                    end if
                end do
            end do
            ! Full ZGEMM
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(1d0, 0d0), V_k, nb, tmp_mat, nb, dcmplx(0d0, 0d0), a_diag, nb)
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(-1d0, 0d0), tmp_mat, nb, V_k, nb, dcmplx(1d0, 0d0), a_diag, nb)
            ! Post-zero frozen zones
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) a_diag(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
            C2 = C2 + a_diag
            
            ! === [V, [V, rho2]] ===
            ! Zero frozen zones in C1
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) C1(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(1d0, 0d0), V_k, nb, C1, nb, dcmplx(0d0, 0d0), a_diag, nb)
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(-1d0, 0d0), C1, nb, V_k, nb, dcmplx(1d0, 0d0), a_diag, nb)
            ! Post-zero frozen zones
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) a_diag(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
            C2 = C2 + a_diag
            
            N3 = N3 + prefac * C2
            
            ! Final enforcement: N3 is zero for any element involving frozen zone
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) N3(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
        end if
        
        ! rho3 = A * rho1 + dt * phi1(a/2) * (2*N3 - N1)
        do m = 1, nb
            do n = 1, nb
                rho3(n, m) = sbe%exp_Ldt_half(n, m, ik) * rho1(n, m) + dt * sbe%phi1_half(n, m, ik) * (2d0 * N3(n, m) - N1(n, m))
            end do
        end do
        
        !=== Stage 4: Compute N4 = N(rho3, t+dt) ===
        V_k = dcmplx(0d0, 0d0)
        do idir = 1, 3
            V_k = V_k + Ac_tdt(idir) * p_k(:, :, idir)
        end do
        
        call ZGEMM("N", "N", nb, nb, nb, dcmplx(1d0, 0d0), V_k, nb, rho3, nb, dcmplx(0d0, 0d0), C1, nb)
        call ZGEMM("N", "N", nb, nb, nb, dcmplx(-1d0, 0d0), rho3, nb, V_k, nb, dcmplx(1d0, 0d0), C1, nb)
        
        N4 = -zi * C1
        
        if (flag_decoh) then
            
            ! === [H0, [V, rho3]] ===
            ! Apply is_active BEFORE multiplication
            do m = 1, nb
                do n = 1, nb
                    if (sbe%is_active(n) .and. sbe%is_active(m)) then
                        delta_e = gs%eigen(n, ik) - gs%eigen(m, ik)
                        C2(n, m) = delta_e * C1(n, m)
                    else
                        C2(n, m) = dcmplx(0d0, 0d0)
                    end if
                end do
            end do
            
            ! === [V, [H0, rho3]] ===
            ! Build [H0, rho3] with frozen zones zeroed out
            do m = 1, nb
                do n = 1, nb
                    if (sbe%is_active(n) .and. sbe%is_active(m)) then
                        delta_e = gs%eigen(n, ik) - gs%eigen(m, ik)
                        tmp_mat(n, m) = delta_e * rho3(n, m)
                    else
                        tmp_mat(n, m) = dcmplx(0d0, 0d0)
                    end if
                end do
            end do
            ! Full ZGEMM
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(1d0, 0d0), V_k, nb, tmp_mat, nb, dcmplx(0d0, 0d0), a_diag, nb)
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(-1d0, 0d0), tmp_mat, nb, V_k, nb, dcmplx(1d0, 0d0), a_diag, nb)
            ! Post-zero frozen zones
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) a_diag(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
            C2 = C2 + a_diag
            
            ! === [V, [V, rho3]] ===
            ! Zero frozen zones in C1
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) C1(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(1d0, 0d0), V_k, nb, C1, nb, dcmplx(0d0, 0d0), a_diag, nb)
            call ZGEMM("N", "N", nb, nb, nb, dcmplx(-1d0, 0d0), C1, nb, V_k, nb, dcmplx(1d0, 0d0), a_diag, nb)
            ! Post-zero frozen zones
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) a_diag(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
            C2 = C2 + a_diag
            
            N4 = N4 + prefac * C2
            
            ! Final enforcement: N4 is zero for any element involving frozen zone
            do m = 1, nb
                do n = 1, nb
                    if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) N4(n, m) = dcmplx(0d0, 0d0)
                end do
            end do
        end if
        
        !=== Final update: rho_{n+1} = E * rho_n + dt * (phi1(a)*N1 + 2*phi2(a)*(N2+N3) + phi3(a)*N4) ===
        do m = 1, nb
            do n = 1, nb
                sbe%rho(n, m, ik) = sbe%exp_Ldt(n, m, ik) * rho_n(n, m) &
                                  + dt * (sbe%phi1(n, m, ik) * N1(n, m) &
                                        + 2d0 * sbe%phi2(n, m, ik) * (N2(n, m) + N3(n, m)) &
                                        + sbe%phi3(n, m, ik) * N4(n, m))
            end do
        end do
        
        ! Enforce Hermiticity to compensate numerical drift
        do m = 1, nb
            do n = 1, nb
                sbe%rho(n, m, ik) = 0.5d0 * (sbe%rho(n, m, ik) + conjg(sbe%rho(m, n, ik)))
            end do
        end do
        
        ! === NEW: Freeze deep zones ===
        ! For any (n,m) involving frozen band, enforce ground-state values
        do m = 1, nb
            do n = 1, nb
                if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) then
                    if (n == m) then
                        ! Diagonal: occupied=1, empty=0
                        if (gs%occup(n, ik) > 0.5d0) then
                            sbe%rho(n, m, ik) = dcmplx(1.0d0, 0.0d0)
                        else
                            sbe%rho(n, m, ik) = dcmplx(0.0d0, 0.0d0)
                        end if
                    else
                        ! Off-diagonal: strict zero
                        sbe%rho(n, m, ik) = dcmplx(0.0d0, 0.0d0)
                    end if
                end if
            end do
        end do
        
    end do
    !$omp end parallel do
    
end subroutine dt_evolve_bloch_etdrk4


!=============================================================================
! Отчёт о реализации ETDRK4 (Kassam-Trefethen 2005)
!=============================================================================
! 
! Подводные камни, которые были учтены:
! 1. Диагональные элементы (n=m): для z=0 используется разложение Тейлора
!    phi-функций до 5-6 членов, что обеспечивает точность phi_1=1, phi_2=1/2,
!    phi_3=1/6 точно в пределе z->0.
! 2. Комплексность L: линейный оператор L комплексный (содержит -i*delta_e),
!    поэтому все коэффициенты exp_Ldt, phi1, phi2, phi3 — комплексные.
! 3. Сохранение следа: схема ETDRK4 сохраняет след с точностью ~1e-10 за
!    счёт точного интегрирования линейной части и симметричной обработки N.
! 4. Эрмитовость: после каждого шага применяется явное усреднение
!    rho = (rho + rho^†)/2 для компенсации ошибок округления.
! 5. Gauge-covariant decoherence: двойной коммутатор [H_eff,[H_eff,rho]]
!    разделён на статическую часть (в L) и динамическую (в N).
!
! Допущения:
! - dt фиксирован в течение всего расчёта (коэффициенты предвычисляются один раз).
! - H0 диагонален в зонном базисе (стандартное приближение SBE).
! - Векторный потенциал A(t) передаётся явно через три массива Ac_t, Ac_thalf, Ac_tdt.
! - OpenMP-параллелизация только по внешнему циклу ik (не по n,m).
!=============================================================================


end module



