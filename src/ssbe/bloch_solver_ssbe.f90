module bloch_solver_ssbe
    use math_constants, only: pi, zi
    use communication, only: comm_get_groupinfo, comm_summation
    use gs_info_ssbe
    use util_ssbe, only: split_range
    implicit none



    type s_sbe_bloch_solver
        !k-points for real-time SBE calculation
        integer :: nk, nb
        integer :: ik_max, ik_min
        complex(8), allocatable :: rho(:, :, :)
        logical :: flag_vnl_correction
    end type



contains



subroutine init_sbe_bloch_solver(sbe, gs, nb_sbe, icomm)
    use util_ssbe
    use communication
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: nb_sbe
    integer, intent(in) :: icomm
    integer :: ik, ib, nk_proc, irank, nproc, ierr
    integer, allocatable :: itbl_min(:), itbl_max(:)

    call comm_get_groupinfo(icomm, irank, nproc)

    sbe%nk = gs%nk
    sbe%nb = nb_sbe

    allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))

    call split_range(1, sbe%nk, nproc, itbl_min, itbl_max)
    sbe%ik_min = itbl_min(irank)
    sbe%ik_max = itbl_max(irank)

    allocate(sbe%rho(1:sbe%nb, 1:sbe%nb, sbe%ik_min:sbe%ik_max))

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
            ! Velocity operator in VG: v = p - A*I (in a.u. e=m=1)
            v_mat = gs%p_tm_matrix(:, :, idir, ik)
            do ib = 1, nb
                v_mat(ib, ib) = v_mat(ib, ib) - Ac(idir)
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
    complex(8), parameter :: zi = dcmplx(0d0, 1d0)
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
        if (0.0d0 < t2_sbe_fs .and. t2_sbe_fs < 1.0d9) then
            t2_au = t2_sbe_fs / au_fs
            prefac = -1.0d0 / (t2_au * gs%eg_au**2)

            ! Form H_eff = diag(eps) + A·p
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


end module



