module bloch_solver_ssbe
    use math_constants, only: pi, zi
    use phys_constants, only: au_ev
    use communication, only: comm_get_groupinfo, comm_summation, comm_bcast
    use gs_info_ssbe
    use util_ssbe, only: split_range
    implicit none

    private
    public :: s_sbe_bloch_solver, init_sbe_bloch_solver, calc_current_bloch, &
              dt_evolve_bloch_etdrk4, calc_trace, calc_energy, &
              init_etdrk4_data, finalize_etdrk4_data

    type s_sbe_bloch_solver
        !k-points for real-time SBE calculation
        integer :: nk, nb
        integer :: ik_max, ik_min
        complex(8), allocatable :: rho(:, :, :)
        logical :: flag_vnl_correction
        
        ! Frozen core handling
        logical, allocatable :: is_active(:) ! .true. if band is active, .false. if frozen
        integer :: n_active_bands = 0
        integer, allocatable :: active_idx(:)  ! Mapping: 1..n_active -> global band index
        
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
    use communication, only: comm_get_groupinfo, comm_summation, comm_bcast
    use salmon_global, only: frozen_core_threshold_ev, frozen_free_threshold_ev
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: nb_sbe
    integer, intent(in) :: icomm
    integer :: ik, ib, nk_proc, irank, nproc, ierr, count_active
    integer, allocatable :: itbl_min(:), itbl_max(:)
    real(8) :: eigen_ev, fermi_energy_ev
    integer, allocatable :: is_active_buf(:)

    call comm_get_groupinfo(icomm, irank, nproc)

    sbe%nk = gs%nk
    sbe%nb = nb_sbe

    allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))

    call split_range(1, sbe%nk, nproc, itbl_min, itbl_max)
    sbe%ik_min = itbl_min(irank)
    sbe%ik_max = itbl_max(irank)

    allocate(sbe%rho(1:sbe%nb, 1:sbe%nb, sbe%ik_min:sbe%ik_max))
    
    ! Calculate Fermi energy (average of HOMO and LUMO) in eV at Gamma point (ik=1)
    fermi_energy_ev = ((gs%eigen(gs%ne/2, 1) + gs%eigen(gs%ne/2 + 1, 1)) * 0.5d0) * au_ev
    
    ! Initialize is_active array based on thresholds relative to Fermi level
    allocate(sbe%is_active(1:sbe%nb))
    sbe%is_active = .false.  
    sbe%n_active_bands = 0   
    
    if (irank == 0) then
        do ib = 1, sbe%nb
            eigen_ev = gs%eigen(ib, 1) * au_ev
            if (eigen_ev > fermi_energy_ev + frozen_core_threshold_ev .and. &
                eigen_ev < fermi_energy_ev + frozen_free_threshold_ev) then
                sbe%is_active(ib) = .true.
                sbe%n_active_bands = sbe%n_active_bands + 1
            else
                sbe%is_active(ib) = .false.
            end if
        end do
    end if
    
    call comm_bcast(sbe%n_active_bands, icomm, 0)
    
    if (sbe%nb > 0) then
        allocate(is_active_buf(1:sbe%nb))
        is_active_buf(:) = 0
        
        if (irank == 0) then
            do ib = 1, sbe%nb
                if (sbe%is_active(ib)) then
                    is_active_buf(ib) = 1
                else
                    is_active_buf(ib) = 0
                end if
            end do
        end if
        
        call comm_bcast(is_active_buf, icomm, 0)
        
        do ib = 1, sbe%nb
            if (is_active_buf(ib) == 1) then
                sbe%is_active(ib) = .true.
            else
                sbe%is_active(ib) = .false.
            end if
        end do
        
        if (allocated(is_active_buf)) deallocate(is_active_buf)
    end if

    if (sbe%n_active_bands > 0) then
        allocate(sbe%active_idx(sbe%n_active_bands))
        count_active = 0
        do ib = 1, sbe%nb
            if (sbe%is_active(ib)) then
                count_active = count_active + 1
                sbe%active_idx(count_active) = ib
            end if
        end do
    else
        allocate(sbe%active_idx(1))  
        sbe%active_idx(1) = 1
    end if

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
            v_mat = gs%p_tm_matrix(:, :, idir, ik)
            do ib = 1, nb
                v_mat(ib, ib) = v_mat(ib, ib) + Ac(idir)
            end do
            if (sbe%flag_vnl_correction) then
                v_mat = v_mat + gs%rvnl_tm_matrix(:, :, idir, ik)
            endif

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


function calc_trace(sbe, gs, nb_max, icomm) result(tr)
    use communication, only: comm_get_groupinfo, comm_summation
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
    
    tmp1 = 0d0
    !$omp parallel do default(shared) private(ik, ib, jb, idir) reduction(+: tmp1)
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            do idir = 1, 3
                do jb = 1, sbe%nb
                    tmp1 = tmp1 &
                        & + Ac(idir) * real(sbe%rho(ib, jb, ik) * gs%p_mod_matrix(jb, ib, idir, ik)) * gs%kweight(ik)
                end do
            end do
            tmp1 = tmp1 &
                & + real(sbe%rho(ib, ib, ik)) * ( &
                & + gs%eigen(ib, ik) &
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
pure subroutine calc_etdrk4_coefficients(z, h, E, E2, Q, f1, f2, f3)
    ! Прямое вычисление коэффициентов ETDRK4.
    ! Использует ряд Тейлора для |z| < 0.1 (избегает деления на ноль)
    ! и точные формулы для |z| >= 0.1 (избегает потери точности).
    implicit none
    complex(8), intent(in) :: z
    real(8), intent(in) :: h
    complex(8), intent(out) :: E, E2, Q, f1, f2, f3
    
    real(8), parameter :: EPS = 0.1d0  ! Порог переключения
    complex(8) :: z2, z3, ez, ez2
    
    E  = exp(z)
    E2 = exp(z * 0.5d0)
    z2 = z * z
    z3 = z2 * z
    
    if (abs(z) < EPS) then
        ! Разложение в ряд Тейлора до O(z^4). 
        ! Точность при |z|<0.1 лучше 10^{-14}.
        ! Q = h * (e^{z/2} - 1) / z
        Q = h * (0.5d0 + z/8.0d0 + z2/48.0d0 + z3/384.0d0)
        
        ! f1 = h * (-4 - z + e^z(4 - 3z + z^2)) / z^3
        f1 = h * (1.0d0/6.0d0 + z/6.0d0 + 3.0d0*z2/40.0d0 + z3/45.0d0)
        
        ! f2 = h * (2 + z + e^z(-2 + z)) / z^3
        f2 = h * (1.0d0/6.0d0 + z/12.0d0 + z2/30.0d0 + z3/180.0d0)
        
        ! f3 = h * (-4 - 3z - z^2 + e^z(4 - z)) / z^3
        f3 = h * (1.0d0/6.0d0 + 0.0d0*z - z2/60.0d0 - z3/360.0d0)
    else
        ! Точные прямые формулы (вне зоны потери значимости)
        ez  = E
        ez2 = E2
        
        Q  = h * (ez2 - 1.0d0) / z
        f1 = h * (-4.0d0 - z + ez * (4.0d0 - 3.0d0*z + z2)) / z3
        f2 = h * (2.0d0 + z + ez * (-2.0d0 + z)) / z3
        f3 = h * (-4.0d0 - 3.0d0*z - z2 + ez * (4.0d0 - z)) / z3
    end if
end subroutine calc_etdrk4_coefficients

subroutine init_etdrk4_data(sbe, gs, dt)
    use phys_constants, only: au_fs
    use salmon_global, only: t2_sbe_fs
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(in) :: dt
    
    integer :: ik, n, m, nb
    real(8) :: delta_e, gamma, t2_au, Eg2
    complex(8) :: lambda, z
    
    nb = sbe%nb
    
    if (.not. allocated(sbe%exp_Ldt)) then
        allocate(sbe%exp_Ldt(nb, nb, sbe%ik_min:sbe%ik_max))
        allocate(sbe%exp_Ldt_half(nb, nb, sbe%ik_min:sbe%ik_max))
        allocate(sbe%phi1(nb, nb, sbe%ik_min:sbe%ik_max))      ! Будет хранить f1
        allocate(sbe%phi2(nb, nb, sbe%ik_min:sbe%ik_max))      ! Будет хранить f2
        allocate(sbe%phi3(nb, nb, sbe%ik_min:sbe%ik_max))      ! Будет хранить f3
        allocate(sbe%phi1_half(nb, nb, sbe%ik_min:sbe%ik_max)) ! Будет хранить Q
    end if
    
    if (t2_sbe_fs > 0.0d0 .and. t2_sbe_fs < 1.0d9) then
        t2_au = t2_sbe_fs / au_fs
        Eg2 = gs%eg_au**2
    else
        t2_au = 1.0d99
        Eg2 = 1.0d0
    end if
    
    !$omp parallel do default(shared) private(ik, n, m, delta_e, gamma, lambda, z)
    do ik = sbe%ik_min, sbe%ik_max
        do m = 1, nb
            do n = 1, nb
                delta_e = gs%eigen(n, ik) - gs%eigen(m, ik)
                
                if (sbe%is_active(n) .and. sbe%is_active(m)) then
                    gamma = (delta_e**2) / (t2_au * Eg2)
                else
                    gamma = 0.0d0
                end if
                
                lambda = dcmplx(0d0, -delta_e) - gamma
                z = lambda * dt
                
                call calc_etdrk4_coefficients(z, dt, &
                    sbe%exp_Ldt(n, m, ik), &
                    sbe%exp_Ldt_half(n, m, ik), &
                    sbe%phi1_half(n, m, ik), &  ! Это Q
                    sbe%phi1(n, m, ik), &       ! Это f1
                    sbe%phi2(n, m, ik), &       ! Это f2
                    sbe%phi3(n, m, ik))         ! Это f3
            end do
        end do
    end do
    !$omp end parallel do
    
    sbe%etdrk4_initialized = .true.
end subroutine init_etdrk4_data


subroutine finalize_etdrk4_data(sbe)
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


pure subroutine calc_phi_contour_all(z, phi_vals)
    implicit none
    complex(8), intent(in) :: z
    complex(8), intent(out) :: phi_vals(3)  
    
    integer, parameter :: M = 32  
    real(8), parameter :: R = 1.0d0  
    real(8), parameter :: TWOPI = 6.28318530717958647692d0
    real(8) :: theta
    complex(8) :: w, ez, w2, phi1_w, phi2_w, phi3_w
    complex(8) :: sum1, sum2, sum3
    integer :: j
    
    sum1 = dcmplx(0d0, 0d0)
    sum2 = dcmplx(0d0, 0d0)
    sum3 = dcmplx(0d0, 0d0)
    
    do j = 1, M
        theta = TWOPI * dble(j) / dble(M)
        w = z + R * exp(dcmplx(0d0, theta))
        
        ez = exp(w)
        w2 = w * w
        
        phi1_w = (ez - dcmplx(1d0, 0d0)) / w
        phi2_w = (ez - dcmplx(1d0, 0d0) - w) / w2
        phi3_w = (ez - dcmplx(1d0, 0d0) - w - 0.5d0 * w2) / (w2 * w)
        
        sum1 = sum1 + phi1_w
        sum2 = sum2 + phi2_w
        sum3 = sum3 + phi3_w
    end do
    
    phi_vals(1) = sum1 / dble(M)
    phi_vals(2) = sum2 / dble(M)
    phi_vals(3) = sum3 / dble(M)
end subroutine calc_phi_contour_all

subroutine dt_evolve_bloch_etdrk4(sbe, gs, Ac, dt)
    use phys_constants, only: au_fs
    use salmon_global, only: t2_sbe_fs
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(inout) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(in) :: dt
    
    integer :: ik, n, m, i, j, in, im, idir
    integer :: nb, nba
    real(8) :: t2_au, prefac, delta_e
    logical :: flag_decoh
    
    complex(8), allocatable :: rho_a(:, :), N_a(:, :), C1_a(:, :), C2_a(:, :), tmp_a(:, :), V_a(:, :)
    complex(8), allocatable :: p_k_full(:, :, :), rho_n_full(:, :)
    complex(8), allocatable :: rho1(:, :), rho2(:, :), rho3(:, :)
    complex(8), allocatable :: N1(:, :), N2(:, :), N3(:, :), N4(:, :)

    nb = sbe%nb
    nba = sbe%n_active_bands
    
    if (.not. sbe%etdrk4_initialized) call init_etdrk4_data(sbe, gs, dt)
    
    if (t2_sbe_fs > 0.0d0 .and. t2_sbe_fs < 1.0d9) then
        t2_au = t2_sbe_fs / au_fs
        prefac = -1.0d0 / (t2_au * gs%eg_au**2)
        flag_decoh = .true.
    else
        flag_decoh = .false.
    end if
    
    !$omp parallel private(rho_a, N_a, C1_a, C2_a, tmp_a, V_a) &
    !$omp            private(p_k_full, rho_n_full, rho1, rho2, rho3, N1, N2, N3, N4) &
    !$omp            shared(sbe, gs, Ac, dt, nb, nba, flag_decoh, prefac)
    
    if (nba > 0) then
        allocate(rho_a(nba, nba), N_a(nba, nba), C1_a(nba, nba), &
                 C2_a(nba, nba), tmp_a(nba, nba), V_a(nba, nba))
    else
        allocate(rho_a(1, 1), N_a(1, 1), C1_a(1, 1), &
                 C2_a(1, 1), tmp_a(1, 1), V_a(1, 1))
    end if
    
    allocate(p_k_full(nb, nb, 1:3))
    allocate(rho_n_full(nb, nb), rho1(nb, nb), rho2(nb, nb), rho3(nb, nb))
    allocate(N1(nb, nb), N2(nb, nb), N3(nb, nb), N4(nb, nb))
    
    !$omp do private(ik, i, j, idir, n, m, in, im, delta_e)
    do ik = sbe%ik_min, sbe%ik_max
        
        p_k_full(:, :, :) = gs%p_tm_matrix(:, :, :, ik)
        if (sbe%flag_vnl_correction) p_k_full(:, :, :) = p_k_full(:, :, :) + gs%rvnl_tm_matrix(:, :, :, ik)
        
        rho_n_full(:, :) = sbe%rho(:, :, ik)
        
        ! V_a = A·p (Константа на шаге для честного сравнения с Тейлором)
        V_a = dcmplx(0d0, 0d0)
        do idir = 1, 3
            do j = 1, nba
                im = sbe%active_idx(j)
                do i = 1, nba
                    in = sbe%active_idx(i)
                    V_a(i, j) = V_a(i, j) + Ac(idir) * p_k_full(in, im, idir)
                end do
            end do
        end do
        
        ! Извлекаем активный блок rho_n
        do j = 1, nba; do i = 1, nba
            in = sbe%active_idx(i); im = sbe%active_idx(j)
            rho_a(i, j) = rho_n_full(in, im)
        end do; end do
        
        ! =====================================================================
        ! STAGE 1: N1 = N(rho_n)  [В статье: Nv]
        ! =====================================================================
        call ZGEMM("N", "N", nba, nba, nba, dcmplx(1d0, 0d0), V_a, nba, rho_a, nba, dcmplx(0d0, 0d0), C1_a, nba)
        call ZGEMM("N", "N", nba, nba, nba, dcmplx(-1d0, 0d0), rho_a, nba, V_a, nba, dcmplx(1d0, 0d0), C1_a, nba)
        N_a = -zi * C1_a
        if (flag_decoh) then
            do j = 1, nba; do i = 1, nba
                in = sbe%active_idx(i); im = sbe%active_idx(j)
                delta_e = gs%eigen(in, ik) - gs%eigen(im, ik)
                C2_a(i, j) = delta_e * C1_a(i, j); tmp_a(i, j) = delta_e * rho_a(i, j)
            end do; end do
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(1d0, 0d0), V_a, nba, tmp_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(-1d0, 0d0), tmp_a, nba, V_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(1d0, 0d0), V_a, nba, C1_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(-1d0, 0d0), C1_a, nba, V_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            N_a = N_a + prefac * C2_a
        end if
        N1 = dcmplx(0d0, 0d0)
        do j = 1, nba; do i = 1, nba
            in = sbe%active_idx(i); im = sbe%active_idx(j); N1(in, im) = N_a(i, j)
        end do; end do
        
        ! a = E2 * rho_n + Q * N1  (rho1 соответствует 'a')
        do m = 1, nb; do n = 1, nb
            rho1(n, m) = sbe%exp_Ldt_half(n, m, ik) * rho_n_full(n, m) + sbe%phi1_half(n, m, ik) * N1(n, m)
        end do; end do
        
        do j = 1, nba; do i = 1, nba
            in = sbe%active_idx(i); im = sbe%active_idx(j); rho_a(i, j) = rho1(in, im)
        end do; end do

        ! =====================================================================
        ! STAGE 2: N2 = N(a)  [В статье: Na]
        ! =====================================================================
        call ZGEMM("N", "N", nba, nba, nba, dcmplx(1d0, 0d0), V_a, nba, rho_a, nba, dcmplx(0d0, 0d0), C1_a, nba)
        call ZGEMM("N", "N", nba, nba, nba, dcmplx(-1d0, 0d0), rho_a, nba, V_a, nba, dcmplx(1d0, 0d0), C1_a, nba)
        N_a = -zi * C1_a
        if (flag_decoh) then
            do j = 1, nba; do i = 1, nba
                in = sbe%active_idx(i); im = sbe%active_idx(j)
                delta_e = gs%eigen(in, ik) - gs%eigen(im, ik)
                C2_a(i, j) = delta_e * C1_a(i, j); tmp_a(i, j) = delta_e * rho_a(i, j)
            end do; end do
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(1d0, 0d0), V_a, nba, tmp_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(-1d0, 0d0), tmp_a, nba, V_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(1d0, 0d0), V_a, nba, C1_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(-1d0, 0d0), C1_a, nba, V_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            N_a = N_a + prefac * C2_a
        end if
        N2 = dcmplx(0d0, 0d0)
        do j = 1, nba; do i = 1, nba
            in = sbe%active_idx(i); im = sbe%active_idx(j); N2(in, im) = N_a(i, j)
        end do; end do
        
        ! b = E2 * rho_n + Q * N2  (rho2 соответствует 'b')
        do m = 1, nb; do n = 1, nb
            rho2(n, m) = sbe%exp_Ldt_half(n, m, ik) * rho_n_full(n, m) + sbe%phi1_half(n, m, ik) * N2(n, m)
        end do; end do
        
        do j = 1, nba; do i = 1, nba
            in = sbe%active_idx(i); im = sbe%active_idx(j); rho_a(i, j) = rho2(in, im)
        end do; end do

        ! =====================================================================
        ! STAGE 3: N3 = N(b)  [В статье: Nb]
        ! =====================================================================
        call ZGEMM("N", "N", nba, nba, nba, dcmplx(1d0, 0d0), V_a, nba, rho_a, nba, dcmplx(0d0, 0d0), C1_a, nba)
        call ZGEMM("N", "N", nba, nba, nba, dcmplx(-1d0, 0d0), rho_a, nba, V_a, nba, dcmplx(1d0, 0d0), C1_a, nba)
        N_a = -zi * C1_a
        if (flag_decoh) then
            do j = 1, nba; do i = 1, nba
                in = sbe%active_idx(i); im = sbe%active_idx(j)
                delta_e = gs%eigen(in, ik) - gs%eigen(im, ik)
                C2_a(i, j) = delta_e * C1_a(i, j); tmp_a(i, j) = delta_e * rho_a(i, j)
            end do; end do
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(1d0, 0d0), V_a, nba, tmp_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(-1d0, 0d0), tmp_a, nba, V_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(1d0, 0d0), V_a, nba, C1_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(-1d0, 0d0), C1_a, nba, V_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            N_a = N_a + prefac * C2_a
        end if
        N3 = dcmplx(0d0, 0d0)
        do j = 1, nba; do i = 1, nba
            in = sbe%active_idx(i); im = sbe%active_idx(j); N3(in, im) = N_a(i, j)
        end do; end do
        
        ! c = E2 * a + Q * (2*N3 - N1)  (rho3 соответствует 'c')
        do m = 1, nb; do n = 1, nb
            rho3(n, m) = sbe%exp_Ldt_half(n, m, ik) * rho1(n, m) + sbe%phi1_half(n, m, ik) * (2d0 * N3(n, m) - N1(n, m))
        end do; end do
        
        do j = 1, nba; do i = 1, nba
            in = sbe%active_idx(i); im = sbe%active_idx(j); rho_a(i, j) = rho3(in, im)
        end do; end do

        ! =====================================================================
        ! STAGE 4: N4 = N(c)  [В статье: Nc]
        ! =====================================================================
        call ZGEMM("N", "N", nba, nba, nba, dcmplx(1d0, 0d0), V_a, nba, rho_a, nba, dcmplx(0d0, 0d0), C1_a, nba)
        call ZGEMM("N", "N", nba, nba, nba, dcmplx(-1d0, 0d0), rho_a, nba, V_a, nba, dcmplx(1d0, 0d0), C1_a, nba)
        N_a = -zi * C1_a
        if (flag_decoh) then
            do j = 1, nba; do i = 1, nba
                in = sbe%active_idx(i); im = sbe%active_idx(j)
                delta_e = gs%eigen(in, ik) - gs%eigen(im, ik)
                C2_a(i, j) = delta_e * C1_a(i, j); tmp_a(i, j) = delta_e * rho_a(i, j)
            end do; end do
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(1d0, 0d0), V_a, nba, tmp_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(-1d0, 0d0), tmp_a, nba, V_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(1d0, 0d0), V_a, nba, C1_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            call ZGEMM("N", "N", nba, nba, nba, dcmplx(-1d0, 0d0), C1_a, nba, V_a, nba, dcmplx(1d0, 0d0), C2_a, nba)
            N_a = N_a + prefac * C2_a
        end if
        N4 = dcmplx(0d0, 0d0)
        do j = 1, nba; do i = 1, nba
            in = sbe%active_idx(i); im = sbe%active_idx(j); N4(in, im) = N_a(i, j)
        end do; end do
        
        ! =====================================================================
        ! FINAL UPDATE 
        ! =====================================================================
        do m = 1, nb; do n = 1, nb
        sbe%rho(n, m, ik) = sbe%exp_Ldt(n, m, ik) * rho_n_full(n, m) &
                      + sbe%phi1(n, m, ik) * N1(n, m) &
                      + 2.0d0 * sbe%phi2(n, m, ik) * (N2(n, m) + N3(n, m)) &
                      + sbe%phi3(n, m, ik) * N4(n, m)
        end do; end do
        
        ! Hermiticity
        do m = 1, nb; do n = 1, nb
            sbe%rho(n, m, ik) = 0.5d0 * (sbe%rho(n, m, ik) + conjg(sbe%rho(m, n, ik)))
        end do; end do
        
        ! Freeze deep zones
        do m = 1, nb; do n = 1, nb
            if (.not. (sbe%is_active(n) .and. sbe%is_active(m))) then
                if (n == m) then
                    if (gs%occup(n, ik) > 0.5d0) then
                        sbe%rho(n, m, ik) = dcmplx(1.0d0, 0.0d0)
                    else
                        sbe%rho(n, m, ik) = dcmplx(0.0d0, 0.0d0)
                    end if
                else
                    sbe%rho(n, m, ik) = dcmplx(0.0d0, 0.0d0)
                end if
            end if
        end do; end do
        
    end do
    !$omp end do
    
    deallocate(rho_a, N_a, C1_a, C2_a, tmp_a, V_a)
    deallocate(p_k_full, rho_n_full, rho1, rho2, rho3, N1, N2, N3, N4)
    !$omp end parallel
    
end subroutine dt_evolve_bloch_etdrk4

end module
