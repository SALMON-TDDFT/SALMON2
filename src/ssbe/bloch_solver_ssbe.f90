module bloch_solver_ssbe
    use math_constants, only: pi, zi
    use phys_constants, only: au_ev
    use communication, only: comm_get_groupinfo, comm_summation, comm_bcast
    use gs_info_ssbe
    use util_ssbe, only: split_range
    implicit none

    private
    public :: s_sbe_bloch_solver, init_sbe_bloch_solver, calc_current_bloch, &
              dt_evolve_bloch_cf4, calc_trace, calc_energy, calc_houston_population_k

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

        ! Houston-basis branch (wave-packet) positions X_a(k,t) used by the
        ! Kuhn-Zurek/Caldeira-Leggett dephasing kernel. By the origin-shift
        ! invariance of the dephasing map (only differences X_a-X_b enter),
        ! the choice X_a(0)=0 is physically irrelevant; it merely fixes a
        ! reproducible reference for restarts.
        real(8), allocatable :: X_branch(:, :)  ! (1:nb, ik_min:ik_max)

        ! Kuhn-Zurek/Caldeira-Leggett decoherence: lambda = kB*T / tau_m
        real(8) :: lambda_decoh = 0d0
        logical :: flag_decoh   = .false.
    end type

    !=========================================================================
    ! CF4 (commutator-free Magnus 4) + Suzuki-Yoshida composition constants
    !=========================================================================
    ! Two-point Gauss-Legendre nodes on [0,1]: c = 1/2 -+ sqrt(3)/6
    real(8), parameter :: cf4_c1 = 0.21132486540518713d0
    real(8), parameter :: cf4_c2 = 0.78867513459481287d0
    ! CF4 combination weights: alpha = 1/4 -+ sqrt(3)/6
    real(8), parameter :: cf4_alpha1 =  0.53867513459481287d0
    real(8), parameter :: cf4_alpha2 = -0.03867513459481287d0
    ! Suzuki-Yoshida triple-jump constants (4th-order composition of a
    ! 2nd-order base scheme): p1 + p2 + p1 = 1
    real(8), parameter :: yoshida_p1 =  1.35120719196d0
    real(8), parameter :: yoshida_p2 = -1.70241438392d0

contains

subroutine init_sbe_bloch_solver(sbe, gs, nb_sbe, icomm)
    use util_ssbe
    use communication, only: comm_get_groupinfo, comm_summation, comm_bcast
    use salmon_global, only: frozen_core_threshold_ev, frozen_free_threshold_ev, &
                             sbe_decoh_temperature_k, sbe_decoh_tau_m_fs
    use phys_constants, only: au_fs, kB_au
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: nb_sbe
    integer, intent(in) :: icomm
    integer :: ik, ib, nk_proc, irank, nproc, ierr, count_active
    integer, allocatable :: itbl_min(:), itbl_max(:)
    real(8) :: eigen_ev, fermi_energy_ev
    integer, allocatable :: is_active_buf(:)
    integer :: homo_idx, lumo_idx

    call comm_get_groupinfo(icomm, irank, nproc)

    sbe%nk = gs%nk
    sbe%nb = nb_sbe

    allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))
    call split_range(1, sbe%nk, nproc, itbl_min, itbl_max)
    sbe%ik_min = itbl_min(irank)
    sbe%ik_max = itbl_max(irank)

    allocate(sbe%rho(1:sbe%nb, 1:sbe%nb, sbe%ik_min:sbe%ik_max))
    
    ! =========================================================================
    ! ИНИЦИАЛИЗАЦИЯ rho: ИСПОЛЬЗУЕМ gs%occup (КАК В ОРИГИНАЛЕ!)
    ! =========================================================================
    sbe%rho(:, :, :) = 0d0
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            sbe%rho(ib, ib, ik) = gs%occup(ib, ik)  ! ← КРИТИЧЕСКОЕ ИСПРАВЛЕНИЕ
        end do
    end do
    
    sbe%flag_vnl_correction = .false.

    ! =========================================================================
    ! Houston-basis branch positions X_a(k,t): explicit zero-init.
    ! By the invariance theorem of the dephasing kernel under a global shift
    ! of the X_a origin (only X_a-X_b enters exp[-lambda(X_a-X_b)^2 tau]),
    ! the choice X_a(t=0) = 0 carries no physical content; it only fixes a
    ! reproducible convention for restarts.
    ! =========================================================================
    allocate(sbe%X_branch(1:sbe%nb, sbe%ik_min:sbe%ik_max))
    sbe%X_branch = 0d0

    ! =========================================================================
    ! Kuhn-Zurek/Caldeira-Leggett decoherence strength: lambda = kB*T / tau_m
    ! Both temperature and relaxation time must be positive to enable it;
    ! otherwise the dephasing map reduces identically to the identity (D=0),
    ! which is trivially CPTP.
    ! =========================================================================
    if (sbe_decoh_temperature_k > 0d0 .and. sbe_decoh_tau_m_fs > 0d0) then
        sbe%lambda_decoh = kB_au * sbe_decoh_temperature_k / (sbe_decoh_tau_m_fs / au_fs)
        sbe%flag_decoh   = .true.
    else
        sbe%lambda_decoh = 0d0
        sbe%flag_decoh   = .false.
    end if

    ! =========================================================================
    ! ЛОГИКА is_active (Frozen Core)
    ! =========================================================================
    
    ! 1. Calculate Fermi Energy (Assumes closed-shell / even number of electrons)
    
    homo_idx = gs%ne / 2
    lumo_idx = homo_idx + 1
    
    if (mod(gs%ne, 2) /= 0 .and. irank == 0) then
        write(*, '(a)') 'WARNING: Odd number of electrons. Fermi energy assumes closed-shell.'
    end if

    fermi_energy_ev = ((gs%eigen(homo_idx, 1) + gs%eigen(lumo_idx, 1)) * 0.5d0) 

    ! 2. Initialize active bands array
    allocate(sbe%is_active(1:sbe%nb))
    sbe%is_active = .false.  
    sbe%n_active_bands = 0   
    
    ! 3. Determine active bands on root rank
    if (irank == 0) then
        do ib = 1, sbe%nb
            eigen_ev = gs%eigen(ib, 1)
            ! Note: Ensure frozen_core_threshold_ev is negative if it represents a window below E_F
            if (eigen_ev > fermi_energy_ev + frozen_core_threshold_ev .and. &
                eigen_ev < fermi_energy_ev + frozen_free_threshold_ev) then
                sbe%is_active(ib) = .true.
                sbe%n_active_bands = sbe%n_active_bands + 1
            end if
        end do
    end if
    
    ! 4. Broadcast n_active_bands
    call comm_bcast(sbe%n_active_bands, icomm, 0)
    
    ! 5. Broadcast is_active logical array
    if (sbe%nb > 0) then
        allocate(is_active_buf(1:sbe%nb))
        
        if (irank == 0) then
            ! Modern Fortran: use merge() instead of a verbose do-loop
            is_active_buf = merge(1, 0, sbe%is_active)
        end if
        
        call comm_bcast(is_active_buf, icomm, 0)
        
        ! Element-wise comparison replaces the verbose do-loop
        sbe%is_active = (is_active_buf == 1)
        
        deallocate(is_active_buf)
    end if

    ! 6. Build active_idx array
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
        ! Modern Fortran handles zero-sized arrays natively. 
        ! If downstream legacy code crashes on size 0, revert to allocate(sbe%active_idx(1))
        allocate(sbe%active_idx(0))  
    end if

    ! 7. Diagnostic Print
    if (irank == 0) then
        write(*, '(a)') '=========================================='
        write(*, '(a)') 'DIAGNOSTIC: Frozen Core Check'
        write(*, '(a, f8.2, a)') '  frozen_core_threshold_ev = ', frozen_core_threshold_ev, ' eV'
        write(*, '(a, f8.2, a)') '  frozen_free_threshold_ev = ', frozen_free_threshold_ev, ' eV'
        write(*, '(a, f12.4, a)') '  Fermi energy (eV)      = ', fermi_energy_ev, ' eV'
        write(*, '(a, i4, a, i4)') '  n_active_bands         = ', sbe%n_active_bands, ' / ', sbe%nb
        write(*, '(a)') '----------------------------------------'
        write(*, '(a)') '  Band energies relative to Fermi level:'
        
        do ib = 1, min(sbe%nb, 100)  ! Print first 100 bands
            eigen_ev = gs%eigen(ib, 1) 
            write(*, '(a, i3, a, f10.4, a, f8.2, a, l1)') &
                '    Band ', ib, ': E = ', eigen_ev, ' eV, E-E_F = ', &
                (eigen_ev - fermi_energy_ev), ' eV, active = ', sbe%is_active(ib)
        end do
        
        if (sbe%nb > 100) write(*, '(a)') '    ... (more bands)'
        write(*, '(a)') '=========================================='
    end if

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
! CF4 (commutator-free Magnus, 4th order) propagator on Gauss-Legendre nodes,
! composed via the Suzuki-Yoshida triple-jump for the unitary part, combined
! with a strictly CPTP Kuhn-Zurek/Caldeira-Leggett dephasing map applied
! through Strang splitting with an exact Hadamard/Gaussian kernel in the
! instantaneous (Houston) eigenbasis:
!
!   rho(t+h) = D(h/2) o [ S2(p1 h) o S2(p2 h) o S2(p1 h) ] o D(h/2) [rho(t)]
!
! IMPORTANT: the Suzuki-Yoshida composition wraps ONLY the unitary CF4
! sub-steps S2(.), never the dephasing map D(.). A unitary step run with a
! negative sub-step (p2 h < 0) is simply a unitary rotation run backwards in
! time -- always valid. A dephasing step run with tau < 0 would replace the
! Hadamard/Gaussian kernel exp[-lambda (X_a-X_b)^2 tau] by its reciprocal,
! exp[+lambda (X_a-X_b)^2 |tau|], which is not positive semi-definite (it
! fails the Schoenberg/Bochner criterion for an RBF kernel and the Schur
! product theorem for the Hadamard map), and would break completely positive
! trace preservation. Hence D(h/2) is applied twice, with tau=+h/2>0 each
! time, by Strang splitting around the (always-safe) unitary composition.
!=============================================================================

subroutine dt_evolve_bloch_cf4(sbe, gs, t_start, dt, Ac_begin, Ac_end)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(inout) :: gs
    real(8), intent(in) :: t_start  ! time at the beginning of the step, rho(t_start) -> rho(t_start+dt)
    real(8), intent(in) :: dt
    real(8), intent(in) :: Ac_begin(1:3)  ! external vector potential A(t_start)
    real(8), intent(in) :: Ac_end(1:3)    ! external vector potential A(t_start+dt)

    real(8) :: tau_sub(3), t_sub(3)
    real(8) :: t_node(2, 3), s_node
    real(8) :: Ac_node(1:3, 2, 3)
    integer :: isub

    integer :: ik, nb, nba, i, j, idir, in, im

    complex(8), allocatable :: p_active(:, :, :)
    complex(8), allocatable :: rho_a(:, :)
    complex(8), allocatable :: H1(:, :), H2(:, :), HVG(:, :)
    real(8),    allocatable :: eigen_active(:)
    real(8),    allocatable :: V_begin(:), V_end(:), X_a(:)
    complex(8), allocatable :: p_k_full(:, :, :)
    complex(8), allocatable :: rho_n_full(:, :)

    nb  = sbe%nb
    nba = sbe%n_active_bands

    !-------------------------------------------------------------------------
    ! The external field is known only at the step endpoints (the analytic
    ! pulse in realtime_ssbe, or the macroscopic Maxwell field in
    ! multiscale_ssbe -- both callers supply Ac(t_start) and Ac(t_start+dt)).
    ! CF4(Gauss-Legendre)+Yoshida needs A at several intermediate sub-nodes;
    ! we obtain them by linear interpolation in time,
    !   A(t_start + s*dt) = (1-s) Ac_begin + s Ac_end,   s in [0,1],
    ! which is consistent with the existing multiscale convention (compare
    ! the "linear interpolation for A(t+dt/2)" used by the previous ETDRK4
    ! step) and strictly more accurate than the old approach of treating A as
    ! constant over the whole step. These nodes are identical for every
    ! k-point, so they are evaluated once before the OpenMP/k-point loop.
    !-------------------------------------------------------------------------
    tau_sub(1) = yoshida_p1 * dt
    tau_sub(2) = yoshida_p2 * dt
    tau_sub(3) = yoshida_p1 * dt

    t_sub(1) = t_start
    t_sub(2) = t_sub(1) + tau_sub(1)
    t_sub(3) = t_sub(2) + tau_sub(2)
    ! t_sub(3) + tau_sub(3) = t_start + dt, since p1 + p2 + p1 = 1

    do isub = 1, 3
        t_node(1, isub) = t_sub(isub) + cf4_c1 * tau_sub(isub)
        t_node(2, isub) = t_sub(isub) + cf4_c2 * tau_sub(isub)

        s_node = (t_node(1, isub) - t_start) / dt
        Ac_node(:, 1, isub) = (1d0 - s_node) * Ac_begin + s_node * Ac_end

        s_node = (t_node(2, isub) - t_start) / dt
        Ac_node(:, 2, isub) = (1d0 - s_node) * Ac_begin + s_node * Ac_end
    end do

    !$omp parallel default(shared) &
    !$omp    private(ik, i, j, idir, in, im, isub) &
    !$omp    private(p_active, rho_a, H1, H2, HVG, eigen_active, V_begin, V_end, X_a) &
    !$omp    private(p_k_full, rho_n_full)

    if (nba > 0) then
        allocate(p_active(nba, nba, 3), rho_a(nba, nba))
        allocate(H1(nba, nba), H2(nba, nba), HVG(nba, nba))
        allocate(eigen_active(nba), V_begin(nba), V_end(nba), X_a(nba))
    end if
    allocate(p_k_full(nb, nb, 1:3), rho_n_full(nb, nb))

    !$omp do
    do ik = sbe%ik_min, sbe%ik_max

        if (nba > 0) then
            p_k_full(:, :, :) = gs%p_tm_matrix(:, :, :, ik)
            if (sbe%flag_vnl_correction) p_k_full(:, :, :) = p_k_full(:, :, :) + gs%rvnl_tm_matrix(:, :, :, ik)
            rho_n_full(:, :) = sbe%rho(:, :, ik)

            ! Restrict to the active subspace (frozen core/free bands excluded)
            do idir = 1, 3
                do j = 1, nba
                    im = sbe%active_idx(j)
                    do i = 1, nba
                        in = sbe%active_idx(i)
                        p_active(i, j, idir) = p_k_full(in, im, idir)
                    end do
                end do
            end do
            do i = 1, nba
                eigen_active(i) = gs%eigen(sbe%active_idx(i), ik)
            end do
            do j = 1, nba; do i = 1, nba
                in = sbe%active_idx(i); im = sbe%active_idx(j)
                rho_a(i, j) = rho_n_full(in, im)
            end do; end do
            do i = 1, nba
                X_a(i) = sbe%X_branch(sbe%active_idx(i), ik)
            end do

            !-----------------------------------------------------------------
            ! Step 1: D(h/2) -- Strang/Hadamard dephasing, tau = +h/2 > 0
            !-----------------------------------------------------------------
            if (sbe%flag_decoh) then
                call build_HVG(nba, eigen_active, p_active, Ac_begin, HVG)
                call houston_dephase(nba, rho_a, HVG, p_active, Ac_begin, X_a, &
                                     sbe%lambda_decoh, 0.5d0 * dt, V_begin)
            else
                V_begin = 0d0
            end if

            !-----------------------------------------------------------------
            ! Step 2: S4_unitary = S2(p1 h) o S2(p2 h) o S2(p1 h)
            ! Each S2(tau) is a CF4 (two-exponential) commutator-free Magnus
            ! step on the two Gauss-Legendre nodes spanning that sub-interval.
            ! A negative tau (the middle Yoshida jump) is just a backward-time
            ! unitary rotation -- exact and unconditionally safe.
            !-----------------------------------------------------------------
            do isub = 1, 3
                call build_HVG(nba, eigen_active, p_active, Ac_node(:, 1, isub), H1)
                call build_HVG(nba, eigen_active, p_active, Ac_node(:, 2, isub), H2)
                call cf4_unitary_step(nba, rho_a, H1, H2, tau_sub(isub))
            end do

            !-----------------------------------------------------------------
            ! Step 3: D(h/2) -- Strang/Hadamard dephasing, tau = +h/2 > 0
            !-----------------------------------------------------------------
            if (sbe%flag_decoh) then
                call build_HVG(nba, eigen_active, p_active, Ac_end, HVG)
                call houston_dephase(nba, rho_a, HVG, p_active, Ac_end, X_a, &
                                     sbe%lambda_decoh, 0.5d0 * dt, V_end)
            else
                V_end = 0d0
            end if

            ! Branch-position update via the midpoint (average of endpoint)
            ! velocities -- consistent with the overall 4th-order accuracy of
            ! CF4 (a forward-Euler X_a += V_a(t_start)*h would degrade the
            ! scheme to 1st order in the branch coordinates).
            do i = 1, nba
                sbe%X_branch(sbe%active_idx(i), ik) = X_a(i) + 0.5d0 * (V_begin(i) + V_end(i)) * dt
            end do

            ! Scatter the evolved active block back into the full matrix
            do j = 1, nba; do i = 1, nba
                in = sbe%active_idx(i); im = sbe%active_idx(j)
                rho_n_full(in, im) = rho_a(i, j)
            end do; end do
            sbe%rho(:, :, ik) = rho_n_full(:, :)
        end if

        ! Hermiticity (numerical safeguard)
        do j = 1, nb; do i = 1, nb
            sbe%rho(i, j, ik) = 0.5d0 * (sbe%rho(i, j, ik) + conjg(sbe%rho(j, i, ik)))
        end do; end do

        ! Freeze deep core/high-energy zones
        do j = 1, nb; do i = 1, nb
            if (.not. (sbe%is_active(i) .and. sbe%is_active(j))) then
                if (i == j) then
                    if (gs%occup(i, ik) > 0.5d0) then
                        sbe%rho(i, j, ik) = dcmplx(2.0d0, 0.0d0)
                    else
                        sbe%rho(i, j, ik) = dcmplx(0.0d0, 0.0d0)
                    end if
                else
                    sbe%rho(i, j, ik) = dcmplx(0.0d0, 0.0d0)
                end if
            end if
        end do; end do

    end do
    !$omp end do

    if (nba > 0) then
        deallocate(p_active, rho_a, H1, H2, HVG, eigen_active, V_begin, V_end, X_a)
    end if
    deallocate(p_k_full, rho_n_full)
    !$omp end parallel

end subroutine dt_evolve_bloch_cf4


! Build the instantaneous velocity-gauge Hamiltonian in the active subspace:
!   H_VG(t) = diag(eigen) + A(t) . pi
subroutine build_HVG(nba, eigen_active, p_active, Ac, H)
    implicit none
    integer,    intent(in)  :: nba
    real(8),    intent(in)  :: eigen_active(nba)
    complex(8), intent(in)  :: p_active(nba, nba, 3)
    real(8),    intent(in)  :: Ac(3)
    complex(8), intent(out) :: H(nba, nba)
    integer :: i, idir

    H = Ac(1) * p_active(:, :, 1) + Ac(2) * p_active(:, :, 2) + Ac(3) * p_active(:, :, 3)
    do i = 1, nba
        H(i, i) = H(i, i) + eigen_active(i)
    end do
end subroutine build_HVG


! Single CF4 (commutator-free Magnus, 4th order) sub-step of length tau,
! evaluated on the two Gauss-Legendre Hamiltonians H1=H(t+c1*tau), H2=H(t+c2*tau):
!   Omega1 = tau (alpha1 H1 + alpha2 H2),  Omega2 = tau (alpha2 H1 + alpha1 H2)
!   rho <- exp(-i Omega2) exp(-i Omega1) rho exp(+i Omega1) exp(+i Omega2)
! Implemented as two successive exact unitary rotations (each built from an
! eigendecomposition of the Hermitian generator, so no Pade/Krylov truncation
! error is introduced -- the propagator is exactly unitary to machine precision).
subroutine cf4_unitary_step(nba, rho, H1, H2, tau)
    implicit none
    integer,    intent(in)    :: nba
    complex(8), intent(inout) :: rho(nba, nba)
    complex(8), intent(in)    :: H1(nba, nba), H2(nba, nba)
    real(8),    intent(in)    :: tau
    complex(8) :: Omega(nba, nba)

    Omega = tau * (cf4_alpha1 * H1 + cf4_alpha2 * H2)
    call apply_unitary_rotation(nba, rho, Omega)

    Omega = tau * (cf4_alpha2 * H1 + cf4_alpha1 * H2)
    call apply_unitary_rotation(nba, rho, Omega)
end subroutine cf4_unitary_step


! Apply rho -> U rho U^dagger with U = exp(-i*Omega) for Hermitian Omega,
! computed exactly via eigendecomposition Omega = W diag(lambda) W^dagger:
!   U rho U^dagger = W [ exp(-i lambda_i) (W^dagger rho W)_ij exp(+i lambda_j) ] W^dagger
subroutine apply_unitary_rotation(nba, rho, Omega)
    use eigen_lapack, only: eigen_zheev
    implicit none
    integer,    intent(in)    :: nba
    complex(8), intent(inout) :: rho(nba, nba)
    complex(8), intent(in)    :: Omega(nba, nba)

    real(8)    :: evals(nba)
    complex(8) :: W(nba, nba), t1(nba, nba), t2(nba, nba)
    integer :: i, j

    call eigen_zheev(Omega, evals, W)

    call ZGEMM('C', 'N', nba, nba, nba, dcmplx(1d0, 0d0), W,  nba, rho, nba, dcmplx(0d0, 0d0), t1, nba)
    call ZGEMM('N', 'N', nba, nba, nba, dcmplx(1d0, 0d0), t1, nba, W,   nba, dcmplx(0d0, 0d0), t2, nba)

    do j = 1, nba
        do i = 1, nba
            t2(i, j) = t2(i, j) * exp(dcmplx(0d0, -(evals(i) - evals(j))))
        end do
    end do

    call ZGEMM('N', 'N', nba, nba, nba, dcmplx(1d0, 0d0), W,  nba, t2, nba, dcmplx(0d0, 0d0), t1, nba)
    call ZGEMM('N', 'C', nba, nba, nba, dcmplx(1d0, 0d0), t1, nba, W,  nba, dcmplx(0d0, 0d0), rho, nba)
end subroutine apply_unitary_rotation


! Strang/Hadamard Kuhn-Zurek/Caldeira-Leggett dephasing step (exactly CPTP for
! any tau >= 0, by the Schoenberg/Bochner positive-definiteness of the
! Gaussian/RBF kernel combined with the Schur product theorem):
!   1) diagonalize the instantaneous H_VG(t) -> Houston (adiabatic) basis U, {E_a}
!   2) rotate rho~ = U^dagger rho U
!   3) rho~_ab <- exp[-lambda (X_a - X_b)^2 * tau] * rho~_ab     (Hadamard product
!      with a positive-semi-definite Gram/RBF matrix => exactly CPTP)
!   4) rotate back rho = U rho~ U^dagger
! Also returns the instantaneous branch (group) velocities in the field
! polarization direction, V_a = [(U^dagger pi U)_aa . e_hat] + (A . e_hat),
! i.e. the projection of v = p + A onto the unit vector of the external
! vector potential (or a fixed reference axis when A ~ 0). These feed the
! midpoint update of the branch positions X_a used by the next dephasing step.
subroutine houston_dephase(nba, rho, H, p_active, Ac, X, lambda, tau, V)
    use eigen_lapack, only: eigen_zheev
    implicit none
    integer,    intent(in)    :: nba
    complex(8), intent(inout) :: rho(nba, nba)
    complex(8), intent(in)    :: H(nba, nba)
    complex(8), intent(in)    :: p_active(nba, nba, 3)
    real(8),    intent(in)    :: Ac(3)
    real(8),    intent(in)    :: X(nba)
    real(8),    intent(in)    :: lambda, tau
    real(8),    intent(out)   :: V(nba)

    real(8)    :: evals(nba), ehat(3), Ac_norm
    complex(8) :: W(nba, nba), t1(nba, nba), t2(nba, nba)
    integer :: i, j, idir

    call eigen_zheev(H, evals, W)

    ! rho~ = U^dagger rho U
    call ZGEMM('C', 'N', nba, nba, nba, dcmplx(1d0, 0d0), W,  nba, rho, nba, dcmplx(0d0, 0d0), t1, nba)
    call ZGEMM('N', 'N', nba, nba, nba, dcmplx(1d0, 0d0), t1, nba, W,   nba, dcmplx(0d0, 0d0), t2, nba)

    ! Exact Hadamard/Gaussian kernel (PSD for tau >= 0)
    do j = 1, nba
        do i = 1, nba
            t2(i, j) = t2(i, j) * exp(-lambda * (X(i) - X(j))**2 * tau)
        end do
    end do

    ! rho = U rho~ U^dagger
    call ZGEMM('N', 'N', nba, nba, nba, dcmplx(1d0, 0d0), W,  nba, t2, nba, dcmplx(0d0, 0d0), t1, nba)
    call ZGEMM('N', 'C', nba, nba, nba, dcmplx(1d0, 0d0), t1, nba, W,  nba, dcmplx(0d0, 0d0), rho, nba)

    ! Branch velocities, projected on the polarization direction of A(t)
    Ac_norm = sqrt(dot_product(Ac, Ac))
    if (Ac_norm > 1.0d-12) then
        ehat = Ac / Ac_norm
    else
        ehat = (/ 1d0, 0d0, 0d0 /)
    end if

    V = 0d0
    do idir = 1, 3
        call ZGEMM('C', 'N', nba, nba, nba, dcmplx(1d0, 0d0), W,  nba, p_active(:, :, idir), nba, dcmplx(0d0, 0d0), t1, nba)
        call ZGEMM('N', 'N', nba, nba, nba, dcmplx(1d0, 0d0), t1, nba, W,                     nba, dcmplx(0d0, 0d0), t2, nba)
        do i = 1, nba
            V(i) = V(i) + ehat(idir) * (real(t2(i, i)) + Ac(idir))
        end do
    end do
end subroutine houston_dephase


! Population of band `ib_target` resolved per k-point, in the instantaneous
! Houston (adiabatic) eigenbasis of H_VG(t) = diag(eigen) + A(t).p -- i.e. the
! same basis used by the CPTP dephasing step. Used to monitor, e.g., the
! occupation of the lowest conduction band as a function of k during the
! real-time propagation:
!   pop_k(ik) = (W^dagger rho W)_{ib_target,ib_target},  H_VG(t) = W diag(E) W^dagger
! The result is summed over MPI ranks (each rank contributes only its own
! k-range, zero elsewhere) so that pop_k is identical and complete on every rank.
subroutine calc_houston_population_k(sbe, gs, Ac, ib_target, pop_k, icomm)
    use eigen_lapack, only: eigen_zheev
    implicit none
    type(s_sbe_bloch_solver), intent(in)  :: sbe
    type(s_sbe_gs_info),      intent(in)  :: gs
    real(8),                  intent(in)  :: Ac(3)
    integer,                  intent(in)  :: ib_target
    real(8),                  intent(out) :: pop_k(1:sbe%nk)
    integer,                  intent(in)  :: icomm

    real(8),    allocatable :: pop_local(:)
    real(8)    :: evals(sbe%nb)
    complex(8) :: H(sbe%nb, sbe%nb), W(sbe%nb, sbe%nb), t1(sbe%nb, sbe%nb), t2(sbe%nb, sbe%nb)
    integer :: ik, i, nb

    nb = sbe%nb
    allocate(pop_local(1:sbe%nk))
    pop_local = 0d0

    do ik = sbe%ik_min, sbe%ik_max
        H(:, :) = Ac(1) * gs%p_tm_matrix(:, :, 1, ik) &
                & + Ac(2) * gs%p_tm_matrix(:, :, 2, ik) &
                & + Ac(3) * gs%p_tm_matrix(:, :, 3, ik)
        if (sbe%flag_vnl_correction) then
            H(:, :) = H(:, :) &
                    & + Ac(1) * gs%rvnl_tm_matrix(:, :, 1, ik) &
                    & + Ac(2) * gs%rvnl_tm_matrix(:, :, 2, ik) &
                    & + Ac(3) * gs%rvnl_tm_matrix(:, :, 3, ik)
        end if
        do i = 1, nb
            H(i, i) = H(i, i) + gs%eigen(i, ik)
        end do

        call eigen_zheev(H, evals, W)

        ! rho~ = W^dagger rho W; diagonal element gives the Houston-basis population
        call ZGEMM('C', 'N', nb, nb, nb, dcmplx(1d0, 0d0), W,  nb, sbe%rho(:, :, ik), nb, dcmplx(0d0, 0d0), t1, nb)
        call ZGEMM('N', 'N', nb, nb, nb, dcmplx(1d0, 0d0), t1, nb, W,                 nb, dcmplx(0d0, 0d0), t2, nb)
        pop_local(ik) = real(t2(ib_target, ib_target))
    end do

    call comm_summation(pop_local, pop_k, sbe%nk, icomm)
    deallocate(pop_local)
end subroutine calc_houston_population_k


end module
