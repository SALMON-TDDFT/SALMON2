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
        complex(8), allocatable :: qnm(:, :, :)
        complex(8), allocatable :: grad_qnm(:, :, :, :)
        complex(8), allocatable :: qnm_new(:, :, :)
        complex(8), allocatable :: dqnm_stock(:, :, :, :)
        complex(8), allocatable :: dnm_i(:, :, :, :)
        real(8), allocatable    :: abs_dnm(:, :, :)
        complex(8), allocatable :: exp_iphi(:, :, :)
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
    use salmon_global, only: norder_correction
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(out) :: jmat(1:3)
    integer, intent(in) :: icomm
    integer :: ik, idir, ib, jb, kb, nb, ierr
    complex(8) :: tmp1(1:3), tmp(1:3)
    complex(8) :: sum1(1:3)
    complex(8) :: pnn(1:3),pni(1:3),pin(1:3)
    complex(8) :: pnj(1:3),pjn(1:3),pnk(1:3),pkn(1:3)
    complex(8) :: pij(1:3),pji(1:3),pik(1:3),pki(1:3),pjk(1:3),pkj(1:3)
    complex(8) :: pnn_Ac,pni_Ac,pin_Ac
    complex(8) :: pnj_Ac,pjn_Ac,pnk_Ac,pkn_Ac
    complex(8) :: pij_Ac,pji_Ac,pik_Ac,pki_Ac,pjk_Ac,pkj_Ac

    tmp1(1:3) = 0d0

    !$omp parallel do default(shared) private(ik,ib,jb,idir) reduction(+:tmp1)
    do ik = sbe%ik_min, sbe%ik_max
        do idir = 1, 3
            do ib = 1, sbe%nb
                do jb = 1, sbe%nb
                    tmp1(idir) = tmp1(idir) + gs%kweight(ik) * sbe%rho(jb, ib, ik) * ( &
                        & gs%p_tm_matrix(ib, jb, idir, ik) &
                    & )
                    if (sbe%flag_vnl_correction) then
                        tmp1(idir) = tmp1(idir) + gs%kweight(ik) * sbe%rho(jb, ib, ik) * ( &
                            & gs%rvnl_tm_matrix(ib, jb, idir, ik) &
                        & )
                    endif
                end do
            end do
        end do
    end do
    !$omp end parallel do

    if(norder_correction>=1)then
        do ik = sbe%ik_min, sbe%ik_max
            do idir = 1, 3
                do nb = 1, gs%ne/2
                    do ib = 1, sbe%nb
                        pin(idir) = gs%p_tm_matrix(ib, nb, idir, ik)
                        pni_Ac = gs%p_tm_matrix(nb, ib, 1, ik) * Ac(1) + &
                                 gs%p_tm_matrix(nb, ib, 2, ik) * Ac(2) + &
                                 gs%p_tm_matrix(nb, ib, 3, ik) * Ac(3)
                        if(nb /= ib) then
                            if(abs(gs%delta_omega(ib, nb, ik))> 1.d-3)then
                                tmp1(idir) = tmp1(idir) - gs%kweight(ik) * &
                                    2.d0 * sbe%rho(nb, nb, ik) * &
                                    dble(pni_Ac*pin(idir)) / gs%delta_omega(ib, nb, ik) 
                            end if
                        end if
                    end do
                end do
            end do
        end do
    end if

    if(norder_correction>=2)then
        do ik = sbe%ik_min, sbe%ik_max
            do idir = 1, 3
                do nb = 1, sbe%nb
                    pnn(idir) = gs%p_tm_matrix(nb, nb, idir, ik)
                    pnn_Ac = gs%p_tm_matrix(nb, nb, 1, ik) * Ac(1) + &
                             gs%p_tm_matrix(nb, nb, 2, ik) * Ac(2) + &
                             gs%p_tm_matrix(nb, nb, 3, ik) * Ac(3)
                    do ib = 1, sbe%nb
                        pni(idir) = gs%p_tm_matrix(nb, ib, idir, ik)
                        pni_Ac = gs%p_tm_matrix(nb, ib, 1, ik) * Ac(1) + &
                                 gs%p_tm_matrix(nb, ib, 2, ik) * Ac(2) + &
                                 gs%p_tm_matrix(nb, ib, 3, ik) * Ac(3)
                        pin_Ac = gs%p_tm_matrix(ib, nb, 1, ik) * Ac(1) + &
                                 gs%p_tm_matrix(ib, nb, 2, ik) * Ac(2) + &
                                 gs%p_tm_matrix(ib, nb, 3, ik) * Ac(3)
                        do jb = 1, sbe%nb
                            pjn(idir) = gs%p_tm_matrix(jb, nb, idir, ik)
                            pij(idir) = gs%p_tm_matrix(ib, jb, idir, ik)
                            pij_Ac = gs%p_tm_matrix(ib, jb, 1, ik) * Ac(1) + &
                                     gs%p_tm_matrix(ib, jb, 2, ik) * Ac(2) + &
                                     gs%p_tm_matrix(ib, jb, 3, ik) * Ac(3)
                            pjn_Ac = gs%p_tm_matrix(jb, nb, 1, ik) * Ac(1) + &
                                     gs%p_tm_matrix(jb, nb, 2, ik) * Ac(2) + &
                                     gs%p_tm_matrix(jb, nb, 3, ik) * Ac(3)
                            if(nb /= ib .and. nb /= jb) then
                                if(abs(gs%delta_omega(ib, nb, ik))> 1.d-3 .and. &
                                   abs(gs%delta_omega(jb, nb, ik))> 1.d-3)then
                                    tmp1(idir) = tmp1(idir) + gs%kweight(ik) * &
                                        sbe%rho(nb, nb, ik) / &
                                        gs%delta_omega(ib, nb, ik) / &
                                        gs%delta_omega(jb, nb, ik) * &
                                        (pjn(idir) * pij_Ac * pni_Ac &
                                       + pjn_Ac * (pij(idir) * pni_Ac  + pni(idir) * pij_Ac))
                                end if
                            end if
                        end do
                        if(nb /= ib) then
                            if(abs(gs%delta_omega(ib, nb, ik))> 1.d-3)then
                                tmp1(idir) = tmp1(idir) - gs%kweight(ik) * &
                                    sbe%rho(nb, nb, ik) * &
                                    (2.d0 * pnn_Ac * dble(pni(idir) * pin_Ac) + &
                                    pnn(idir)*abs(pni_Ac)**2) / &
                                    gs%delta_omega(ib, nb, ik)**2
                            end if
                        end if
                    end do
                end do
            end do
        end do
    end if

    if(norder_correction>=3)then
      do ik = sbe%ik_min, sbe%ik_max
        do idir = 1, 3
          do nb = 1, sbe%nb
            pnn(idir) = gs%p_tm_matrix(nb, nb, idir, ik)
            pnn_Ac = gs%p_tm_matrix(nb, nb, 1, ik) * Ac(1) + &
                     gs%p_tm_matrix(nb, nb, 2, ik) * Ac(2) + &
                     gs%p_tm_matrix(nb, nb, 3, ik) * Ac(3)
            do ib = 1, sbe%nb
              pni(idir) = gs%p_tm_matrix(nb, ib, idir, ik)
              pin(idir) = gs%p_tm_matrix(ib, nb, idir, ik)
              pni_Ac = gs%p_tm_matrix(nb, ib, 1, ik) * Ac(1) + &
                       gs%p_tm_matrix(nb, ib, 2, ik) * Ac(2) + &
                       gs%p_tm_matrix(nb, ib, 3, ik) * Ac(3)
              pin_Ac = gs%p_tm_matrix(ib, nb, 1, ik) * Ac(1) + &
                       gs%p_tm_matrix(ib, nb, 2, ik) * Ac(2) + &
                       gs%p_tm_matrix(ib, nb, 3, ik) * Ac(3)
              do jb = 1, sbe%nb
                pnj(idir) = gs%p_tm_matrix(nb, jb, idir, ik)
                pjn(idir) = gs%p_tm_matrix(jb, nb, idir, ik)
                pij(idir) = gs%p_tm_matrix(ib, jb, idir, ik)
                pji(idir) = gs%p_tm_matrix(jb, ib, idir, ik)
                pnj_Ac = gs%p_tm_matrix(nb, jb, 1, ik) * Ac(1) + &
                         gs%p_tm_matrix(nb, jb, 2, ik) * Ac(2) + &
                         gs%p_tm_matrix(nb, jb, 3, ik) * Ac(3)
                pjn_Ac = gs%p_tm_matrix(jb, nb, 1, ik) * Ac(1) + &
                         gs%p_tm_matrix(jb, nb, 2, ik) * Ac(2) + &
                         gs%p_tm_matrix(jb, nb, 3, ik) * Ac(3)
                pij_Ac = gs%p_tm_matrix(ib, jb, 1, ik) * Ac(1) + &
                         gs%p_tm_matrix(ib, jb, 2, ik) * Ac(2) + &
                         gs%p_tm_matrix(ib, jb, 3, ik) * Ac(3)
                pji_Ac = gs%p_tm_matrix(jb, ib, 1, ik) * Ac(1) + &
                         gs%p_tm_matrix(jb, ib, 2, ik) * Ac(2) + &
                         gs%p_tm_matrix(jb, ib, 3, ik) * Ac(3)
                do kb = 1, sbe%nb
                  pnk(idir) = gs%p_tm_matrix(nb, kb, idir, ik)
                  pkn(idir) = gs%p_tm_matrix(kb, nb, idir, ik)
                  pik(idir) = gs%p_tm_matrix(ib, kb, idir, ik)
                  pki(idir) = gs%p_tm_matrix(kb, ib, idir, ik)
                  pjk(idir) = gs%p_tm_matrix(jb, kb, idir, ik)
                  pkj(idir) = gs%p_tm_matrix(kb, jb, idir, ik)
                  pnk_Ac = gs%p_tm_matrix(nb, kb, 1, ik) * Ac(1) + &
                           gs%p_tm_matrix(nb, kb, 2, ik) * Ac(2) + &
                           gs%p_tm_matrix(nb, kb, 3, ik) * Ac(3)
                  pkn_Ac = gs%p_tm_matrix(kb, nb, 1, ik) * Ac(1) + &
                           gs%p_tm_matrix(kb, nb, 2, ik) * Ac(2) + &
                           gs%p_tm_matrix(kb, nb, 3, ik) * Ac(3)
                  pik_Ac = gs%p_tm_matrix(ib, kb, 1, ik) * Ac(1) + &
                           gs%p_tm_matrix(ib, kb, 2, ik) * Ac(2) + &
                           gs%p_tm_matrix(ib, kb, 3, ik) * Ac(3)
                  pki_Ac = gs%p_tm_matrix(kb, ib, 1, ik) * Ac(1) + &
                           gs%p_tm_matrix(kb, ib, 2, ik) * Ac(2) + &
                           gs%p_tm_matrix(kb, ib, 3, ik) * Ac(3)
                  pjk_Ac = gs%p_tm_matrix(jb, kb, 1, ik) * Ac(1) + &
                           gs%p_tm_matrix(jb, kb, 2, ik) * Ac(2) + &
                           gs%p_tm_matrix(jb, kb, 3, ik) * Ac(3)
                  pkj_Ac = gs%p_tm_matrix(kb, jb, 1, ik) * Ac(1) + &
                           gs%p_tm_matrix(kb, jb, 2, ik) * Ac(2) + &
                           gs%p_tm_matrix(kb, jb, 3, ik) * Ac(3)
                  if(nb /= ib .and. nb /= jb .and. nb /= kb) then
                    if(abs(gs%delta_omega(ib, nb, ik))> 1.d-3 .and. &
                       abs(gs%delta_omega(jb, nb, ik))> 1.d-3 .and. &
                       abs(gs%delta_omega(kb, nb, ik))> 1.d-3)then
                      tmp1(idir) = tmp1(idir) - gs%kweight(ik) * &
                                   sbe%rho(nb, nb, ik) / &
                                   gs%delta_omega(ib, nb, ik) / &
                                   gs%delta_omega(jb, nb, ik) / &
                                   gs%delta_omega(kb, nb, ik) * &
                                   dble(pji_Ac * ( pnk_Ac * ( pin(idir) * pkj_Ac + &
                                                              pkj(idir) * pin_Ac ) + &
                                                   pnk(idir) * pin_Ac * pkj_Ac ) + &
                                        pji(idir) * pin_Ac * pkj_Ac * pnk_Ac)
                    end if
                  end if
                end do
                if(nb /= ib .and. nb /= jb) then
                  if(abs(gs%delta_omega(ib, nb, ik))> 1.d-3 .and. &
                     abs(gs%delta_omega(jb, nb, ik))> 1.d-3)then
                    tmp1(idir) = tmp1(idir) + gs%kweight(ik) * &
                                 sbe%rho(nb, nb, ik) * &
                                 (gs%delta_omega(ib, nb, ik) + gs%delta_omega(jb, nb, ik)) / &
                                 2.d0 / &
                                 gs%delta_omega(ib, nb, ik)**2 / &
                                 gs%delta_omega(jb, nb, ik)**2 * &
                                 (pji(idir) * pnn_Ac * pin_Ac * pnj_Ac + &
                                  pjn(idir) * pni_Ac * (pin_Ac * pnj_Ac + &
                                                        pnn_Ac * pij_Ac) + &
                                  pji_Ac * (pnn_Ac * (pin(idir) * pnj_Ac + &
                                                      pnj(idir) * pin_Ac) + &
                                            pnn(idir) * pin_Ac * pnj_Ac) + &
                                  pjn_Ac * (pni(idir) * ( pin_Ac * pnj_Ac + &
                                                          pnn_Ac * pij_Ac ) + &
                                            pni_Ac * (pij(idir) * pnn_Ac + &
                                                      pin(idir) * pnj_Ac + &
                                                      pnj(idir) * pin_Ac + &
                                                      pnn(idir) * pij_Ac)))
                  end if
                end if
              end do
            end do
          end do
        end do
      end do
      do ik = sbe%ik_min, sbe%ik_max
        do idir = 1, 3
          do nb = 1, sbe%nb
            pnn(idir) = gs%p_tm_matrix(nb, nb, idir, ik)
            pnn_Ac = gs%p_tm_matrix(nb, nb, 1, ik) * Ac(1) + &
                     gs%p_tm_matrix(nb, nb, 2, ik) * Ac(2) + &
                     gs%p_tm_matrix(nb, nb, 3, ik) * Ac(3)
            sum1 = 0.d0
            do ib = 1, sbe%nb
              pni(idir) = gs%p_tm_matrix(nb, ib, idir, ik)
              pin(idir) = gs%p_tm_matrix(ib, nb, idir, ik)
              pni_Ac = gs%p_tm_matrix(nb, ib, 1, ik) * Ac(1) + &
                       gs%p_tm_matrix(nb, ib, 2, ik) * Ac(2) + &
                       gs%p_tm_matrix(nb, ib, 3, ik) * Ac(3)
              pin_Ac = gs%p_tm_matrix(ib, nb, 1, ik) * Ac(1) + &
                       gs%p_tm_matrix(ib, nb, 2, ik) * Ac(2) + &
                       gs%p_tm_matrix(ib, nb, 3, ik) * Ac(3)
              if(nb /= ib) then
                if(abs(gs%delta_omega(ib, nb, ik))> 1.d-3)then
                  sum1(idir) = sum1(idir) + &
                          (pni(idir) * pnn_Ac * pin_Ac + &
                           pni_Ac * (pin(idir) * pnn_Ac + &
                                    2.d0 * pnn(idir) * pin_Ac)) / &
                           gs%delta_omega(ib, nb, ik)**3
                end if
              end if
            end do
            tmp1(idir) = tmp1(idir) - gs%kweight(ik) * sbe%rho(nb, nb, ik) * pnn_Ac * sum1(idir)
          end do
        end do
      end do
    end if

    call comm_summation(tmp1, tmp, 3, icomm)

    jmat(:) = (real(tmp(1:3)) / sum(gs%kweight(:)) &
        & + Ac * calc_trace(sbe, gs, sbe%nb, icomm)) / gs%volume

    return
end subroutine calc_current_bloch


subroutine dt_evolve_bloch(sbe, gs, Ac, dt)
    use salmon_global, only: yn_sbe_gw_collision, sbe_deph_mode
    use sbe_collision_gw, only: add_collision_vg
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
    end do
    ! GW collision term (Phase 2), velocity gauge: explicit forward-Euler step
    ! applied to rho after the Taylor propagation.  OFF path untouched.
    if (yn_sbe_gw_collision == 'y') then
        call add_collision_vg(sbe%rho, gs%gamma_gw, gs%f0_ref, dt, &
          & sbe%nb, sbe%nk, sbe%ik_min, sbe%ik_max, sbe_deph_mode)
    end if
    return

contains


    !Calculate [H, rho] commutation:
    subroutine calc_hrho_bloch_k(ik, rho_k, p_k, hrho_k)
        implicit none
        integer, intent(in) :: ik
        complex(8), intent(in) :: rho_k(nb, nb)
        complex(8), intent(in) :: p_k(nb, nb, 1:3)
        complex(8), intent(out) :: hrho_k(nb, nb)
        integer :: idir
        !hrho = hrho + Ac(t) * (p * rho - rho * p)
        hrho_k(1:nb, 1:nb) = gs%delta_omega(1:nb, 1:nb, ik) * rho_k(1:nb, 1:nb)
        do idir=1, 3 !1:x, 2:y, 3:z
            ! hrho(1:nb, 1:nb, ik) = hrho(1:nb, 1:nb, ik) + Ac(idir) * (&
            ! & + matmul(gs%p_mod_matrix(1:nb, 1:nb, idir, ik), rho(1:nb, 1:nb, ik)) &
            ! & - matmul(rho(1:nb, 1:nb, ik), gs%p_mod_matrix(1:nb, 1:nb, idir, ik)) &
            ! & )

            call ZGEMM("N","N", nb, nb, nb, &
                dcmplx(+Ac(idir), 0d0), &
                p_k(:, :, idir),nb, &
                rho_k(:, :), nb, &
                dcmplx(1d0, 0d0), hrho_k(:, :),nb)

            call ZGEMM("N","N", nb, nb, nb, &
                dcmplx(-Ac(idir), 0d0), &
                rho_k(:, :), nb, &
                p_k(:, :, idir),nb, &
                dcmplx(1d0, 0d0), hrho_k(:, :), nb)

        end do !idir
        return
    end subroutine calc_hrho_bloch_k
end subroutine

function calc_trace(sbe, gs, nb_max, icomm) result(tr)
    use communication
    use salmon_global
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: icomm
    integer, intent(in) :: nb_max
    real(8) :: tr

    integer :: ik, ib
    real(8) :: tmp, tmp1

    tmp1 = 0d0
    select case(trim(gauge_sbe))
    case ("velocity_gauge")
        !$omp parallel do default(shared) private(ik, ib) reduction(+: tmp1) collapse(2)
        do ik = sbe%ik_min, sbe%ik_max
            do ib = 1, nb_max
                tmp1 = tmp1 + real(sbe%rho(ib, ib, ik)) * gs%kweight(ik)
            end do
        end do
        !$omp end parallel do
    case ("length_gauge")
        !$omp parallel do default(shared) private(ik, ib) reduction(+: tmp1) collapse(2)
        do ik = sbe%ik_min, sbe%ik_max
            do ib = 1, nb_max
                tmp1 = tmp1 + real(sbe%qnm_new(ib, ib, ik)) * gs%kweight(ik)
            end do
        end do
        !$omp end parallel do
    end select
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

subroutine prepare_qnm(sbe, gs, icomm)
  use salmon_global, only: epdir_re1, am_s, sbe_lg_diag, sbe_lg_degen
  implicit none
  type(s_sbe_bloch_solver), intent(inout) :: sbe
  type(s_sbe_gs_info), intent(in) :: gs
  integer, intent(in) :: icomm
  complex(8) :: dnm(sbe%nb, sbe%nb, sbe%nk)
  integer :: nb,nk
  integer :: ik,ib,jb,ii,jj

  nk = sbe%nk
  nb = sbe%nb

  allocate(sbe%qnm(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max))
  allocate(sbe%grad_qnm(sbe%nb, sbe%nb, 1:3, sbe%nk))
  allocate(sbe%qnm_new(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max))
  allocate(sbe%dqnm_stock(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max, am_s))
  allocate(sbe%dnm_i(sbe%nb, sbe%nb, 1:3, sbe%ik_min:sbe%ik_max))
  allocate(sbe%abs_dnm(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max))
  allocate(sbe%exp_iphi(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max))

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib=1,nb
  do jb=1,nb
    sbe%qnm(ib, jb, ik) = 0.d0
    sbe%qnm_new(ib, jb, ik) = 0.d0
    dnm(ib, jb, ik)=0.d0
    sbe%abs_dnm(ib, jb, ik)=0.d0
    sbe%exp_iphi(ib, jb, ik)=0.d0
  end do
  end do
  end do

  do ii=1,am_s
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib=1,nb
  do jb=1,nb
      sbe%dqnm_stock(ib, jb, ik, ii) = 0.d0
  end do
  end do
  end do
  end do

  !$omp parallel do default(shared) private(ik, ib, jb, jj) 
  do ik = sbe%ik_min, sbe%ik_max
    do jj=1,3
      do ib=1,nb
        do jb=1,nb
          dnm(ib, jb, ik) = dnm(ib, jb, ik) + epdir_re1(jj) * gs%d_matrix(ib, jb, jj, ik)
        end do
      end do
    end do
  end do

  !$omp parallel do default(shared) private(ik, ib, jb, jj) collapse(4)
  do ik = sbe%ik_min, sbe%ik_max
    do jj=1,3
      do ib=1,nb
        do jb=1,nb
          sbe%dnm_i(ib, jb, jj, ik) = epdir_re1(jj) * gs%d_matrix(ib, jb, jj, ik)
        end do
      end do
    end do
  end do

  sbe%abs_dnm=0.d0
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
    do ib=1,nb
      do jb=1,nb
        sbe%abs_dnm(ib, jb, ik) = abs(dnm(ib, jb, ik))
      end do
    end do
  end do

  ! LG-SBE diagnostic knockout (2): suppress near-degenerate dipole couplings.
  ! gs%delta_omega is (nb,nb,1:nk) on every rank; sbe%abs_dnm is the local
  ! k-slice (nb,nb,ik_min:ik_max), so slice delta_omega to make the WHERE
  ! conform (identical to the full array when nproc==1).
  if (sbe_lg_diag == 2 .or. sbe_lg_diag == 3) then
    where (abs(gs%delta_omega(:, :, sbe%ik_min:sbe%ik_max)) < 1.d-3) sbe%abs_dnm = 0.d0
  end if

  ! gicov X-full (R-1 representation): ALL off-diagonal coherences are carried
  ! in the PHYSICAL basis (qnm == rho).  Force exp_iphi=1 and abs_dnm=0 on
  ! every off-diagonal pair BEFORE the phase-fill below, so (i) the
  ! phase-fill's abs_dnm>=1.d-13 gate skips them and leaves exp_iphi=1, giving
  ! qnm_new == rho for every off-diagonal pair, and (ii) the dipole source
  ! terms (which multiply abs_dnm) vanish everywhere -- the interband coupling
  ! is supplied instead by the full-band gauge-covariant gradient (gicov_rhs),
  ! not a dipole source.  X-full has a SINGLE full-band block (block_id === 1),
  ! so "same-block" and "ib/=jb" now coincide; the loop body no longer
  ! references block_id.  gi/gifix/off are untouched (this whole block is
  ! guarded by sbe_lg_degen=='gicov').
  if (trim(sbe_lg_degen) == 'gicov') then
    !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
    do ik = sbe%ik_min, sbe%ik_max
      do ib=1,nb
        do jb=1,nb
          if (ib /= jb) then
            sbe%abs_dnm(ib, jb, ik) = 0.d0
            sbe%exp_iphi(ib, jb, ik) = (1.d0, 0.d0)
          end if
        end do
      end do
    end do
  end if

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
    do ib=1,nb
      do jb=1,nb
        if(sbe%abs_dnm(ib, jb, ik) >= 1.d-13) then
          sbe%exp_iphi(ib, jb, ik) = dnm(ib, jb, ik)/sbe%abs_dnm(ib, jb, ik)
        end if
        if(ib==jb)then
          sbe%qnm_new(ib, jb, ik) = sbe%rho(ib, jb, ik)
        else
          sbe%qnm_new(ib, jb, ik) = conjg(sbe%exp_iphi(ib, jb, ik)) *  &
                                 &  sbe%rho(ib, jb, ik)
        end if
      end do
    end do
  end do

end subroutine prepare_qnm

!===================================================================
! R-1 representation bridge (orientation-explicit).  qnm is the propagated
! variable; the covariant gradient needs PHYSICAL rho.
!   rho(i,j) =              qnm(i,j)                (i==j)
!            = exp_iphi(i,j) * qnm(i,j)             (i/=j)
!   qnm(i,j) =              rho(i,j)                (i==j)
!            = conjg(exp_iphi(i,j)) * rho(i,j)      (i/=j)
! In gicov (X-full), prepare_qnm sets exp_iphi=1 on EVERY off-diagonal pair
! (ib/=jb), not just same-block ones, so qnm == rho for all off-diagonal
! coherences -- there is no out-of-block dipole-derived phase left; the
! interband coupling is carried entirely by the full-band gauge-covariant
! gradient (gicov_rhs), not by a dipole-phase-carrying qnm/rho distinction.
!===================================================================
pure function rho_ij_from_q(sbe, ib, jb, ik) result(r)
  implicit none
  type(s_sbe_bloch_solver), intent(in) :: sbe
  integer, intent(in) :: ib, jb, ik
  complex(8) :: r
  if (ib == jb) then
    r = sbe%qnm(ib, jb, ik)
  else
    r = sbe%exp_iphi(ib, jb, ik) * sbe%qnm(ib, jb, ik)
  end if
end function rho_ij_from_q

pure function q_ij_from_rho(sbe, rho_ij, ib, jb, ik) result(q)
  implicit none
  type(s_sbe_bloch_solver), intent(in) :: sbe
  complex(8), intent(in) :: rho_ij
  integer, intent(in) :: ib, jb, ik
  complex(8) :: q
  if (ib == jb) then
    q = rho_ij
  else
    q = conjg(sbe%exp_iphi(ib, jb, ik)) * rho_ij
  end if
end function q_ij_from_rho

!===================================================================
! gicov RHS operator (Phase 3, X-full): instantaneous COHERENT drho/dt of the
! gauge-covariant length-gauge propagator, as a callable routine.  NO
! integrator / AB4 / dqnm_stock here (that is Task 5b) -- this returns the
! instantaneous physical drho only, so a RHS bug can be told apart from an
! integrator-stability issue.
!
!   drho = + sum_a E_a * D_cov[rho]_a         covariant transport (WHOLE field term) (1)
!          - i * delta_omega .* rho           band-energy coherent term       (2a; off-diag)
!          - rho / t_2                         legacy interband dephasing      (2b; off-diag)
!
! (1) covariant_grad_block returns D_cov rho = d_k rho - i[xi,rho] on the FULL
!     nb x nb density; it REPLACES the legacy grad_qnm term.  gs%u_transport is
!     fed AS-IS (orientation U, Task 2 -- NOT conjugate-transposed).  The
!     Cartesian k-spacing dk_a = b_matrix(a,a)/num_kgrid(a) reproduces the
!     legacy grad_k_array_nb2d_dcomplex spacing (nabt = bNmat(:,4)/(b(a,a)/N_a)),
!     so with U==I this reduces EXACTLY to the legacy gradient; the "+E*grad"
!     field prefactor/sign matches dt_evolve_bloch_lg (:663,:685).
!     X-full (block_id === 1, the full nb x nb overlap polar-factored as ONE
!     block): xi is now the FULL-BAND gauge-covariant connection, so its
!     off-diagonal entries (xi_inter) already ARE the physical interband
!     dipole -- the separate analytic dipole commutator that Phase-3-fixed-block
!     gicov used to add here would double-count it, so it is DELETED (proven
!     equivalent, not merely "unneeded", by src/ssbe/test/test_gicov_xfull.f90's
!     Test D: a REAL sigma_y rotation gives zero diagonal Berry connection, so
!     xi_full = xi_inter, and covariant_grad_block(full) is shown to equal
!     covariant_grad_block(bare) - i[xi_inter,rho] to stencil-matched tol).
!     gs%d_matrix is NO LONGER READ by this routine.
! (2) energy -i*delta_omega*rho (= -i[H0,rho]) and the 1/t_2 dephasing are kept
!     exactly as the legacy LG RHS (:680,:690); when GW dephasing is active the
!     t_2 term is suppressed (5b applies the GW collision as a separate substep).
!
! rho is reconstructed from sbe%qnm via rho_ij_from_q and gathered across ranks
! (comm_summation) because the covariant stencil is non-local in k.  Returns
! PHYSICAL drho (full nb x nb x nk).  5b maps drho -> dqnm elementwise via
! q_ij_from_rho if it steps qnm (diag: dqnm=drho; off-diag: dqnm=conjg(exp_iphi)*drho).
!
! NOTE (deviation from the brief signature): icomm is added as a trailing
! argument.  The covariant gradient must see neighbouring-k densities that may
! live on other ranks, so a gather is mandatory -- exactly as the legacy
! dt_evolve_bloch_lg gathers qnm via comm_summation before grad_k.  drho is
! returned for the full k range (nb,nb,nk); a caller slices its local k-range.
!===================================================================
subroutine gicov_rhs(sbe, gs, Efield, drho, icomm)
  use salmon_global, only: num_kgrid, t_2, yn_sbe_gw_collision, sbe_deph_mode
  use degenerate_block_ssbe, only: covariant_grad_block, theta_off
  implicit none
  type(s_sbe_bloch_solver), intent(in) :: sbe
  type(s_sbe_gs_info), intent(in) :: gs
  real(8), intent(in) :: Efield(1:3)
  complex(8), intent(out) :: drho(sbe%nb, sbe%nb, sbe%nk)
  integer, intent(in) :: icomm
  integer :: nb, nk, ik, ib, jb, axis
  real(8) :: dk(1:3)
  logical :: deph_by_gw
  complex(8), allocatable :: rho_loc(:,:,:), rho_full(:,:,:), Dq(:,:,:,:)
  complex(8) :: gterm

  nb = sbe%nb
  nk = sbe%nk

  allocate(rho_loc(nb,nb,nk), rho_full(nb,nb,nk), Dq(nb,nb,3,nk))

  ! ---- (0) physical rho on the local k-slice (zero elsewhere), gather to full k
  rho_loc = (0.d0, 0.d0)
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
    do ib = 1, nb
      do jb = 1, nb
        rho_loc(ib, jb, ik) = rho_ij_from_q(sbe, ib, jb, ik)
      end do
    end do
  end do
  call comm_summation(rho_loc, rho_full, nb*nb*nk, icomm)

  ! ---- (1) covariant transport (physical), WHOLE field term: + sum_a E_a D_cov[rho]_a
  do axis = 1, 3
    dk(axis) = gs%b_matrix(axis, axis) / dble(num_kgrid(axis))
  end do
  call covariant_grad_block(nb, nk, gs%nbvec, gs%bvec, num_kgrid, &
                            gs%u_transport, rho_full, dk, Dq, sbe%ik_min, sbe%ik_max)

  deph_by_gw = (yn_sbe_gw_collision == 'y' .and. trim(sbe_deph_mode) == 'gw')

  drho = (0.d0, 0.d0)
  ! only the local k-slice of drho is produced (the caller uses only ik_min:ik_max);
  ! Dq was likewise computed only on [ik_min:ik_max] by covariant_grad_block.
  !$omp parallel do default(shared) private(ik, ib, jb, gterm)
  do ik = sbe%ik_min, sbe%ik_max
    do ib = 1, nb
      do jb = 1, nb
        ! (1) covariant transport (intraband + interband, full-band)
        gterm = Efield(1) * Dq(ib, jb, 1, ik) + &
                Efield(2) * Dq(ib, jb, 2, ik) + &
                Efield(3) * Dq(ib, jb, 3, ik)
        drho(ib, jb, ik) = gterm
        ! (2) coherent band energy + phenomenological dephasing (off-diagonal only)
        if (ib /= jb) then
          drho(ib, jb, ik) = drho(ib, jb, ik) &
            & - zi * gs%delta_omega(ib, jb, ik) * rho_full(ib, jb, ik)
          ! Covariant T2: dephase only ENERGY-DISTINCT pairs. Inside a
          ! (near-)degenerate manifold (|delta_omega| <= theta_off, the same
          ! degeneracy scale the Wilson-line transport groups blocks by) the
          ! scalar -rho/t_2 relaxation is NOT U(N)-gauge-covariant -- it splits
          ! diagonal/off-diagonal, a partition U(N) rotations mix -- so it is
          ! skipped; intra-manifold relaxation belongs to the covariant (GW
          ! collision) channel, not to a phenomenological scalar T2.
          if (.not. deph_by_gw .and. &
            & abs(gs%delta_omega(ib, jb, ik)) > theta_off) then
            drho(ib, jb, ik) = drho(ib, jb, ik) - rho_full(ib, jb, ik) / t_2
          end if
        end if
      end do
    end do
  end do
  !$omp end parallel do

  deallocate(rho_loc, rho_full, Dq)
end subroutine gicov_rhs

subroutine dt_evolve_bloch_lg(sbe, gs, E, bj_am, dt, icomm)
  use salmon_global, only: am_s, num_kgrid, t_2, epdir_re1, &
    & yn_sbe_gw_collision, sbe_deph_mode, sbe_lg_diag, sbe_lg_degen
  use sbe_collision_gw, only: add_collision_diag, add_collision_offdiag
  use common_ssbe, only: grad_k_array_nb2d_dcomplex
  implicit none
  type(s_sbe_bloch_solver), intent(inout) :: sbe
  type(s_sbe_gs_info), intent(inout) :: gs
  real(8), intent(in) :: E(1:3)
  real(8), intent(in) :: bj_am(8,8)
  real(8), intent(in) :: dt
  integer, intent(in) :: icomm
  integer :: nb,nk
  integer :: ik,ib,jb
  integer :: ii,lb
  complex(8) :: shift_vector(3)
  real(8) :: abs_E, proj_E
  complex(8) :: qnm_tmp1(sbe%nb, sbe%nb, sbe%nk)
  complex(8) :: qnm_tmp(sbe%nb, sbe%nb, sbe%nk)

  nb = sbe%nb
  nk = sbe%nk

  ! gicov (R-3, Task 5b): self-starting RK4 on the PHYSICAL density via gicov_rhs
  ! (imaginary-axis stability 2*sqrt(2)~2.83), NO dqnm_stock / Adams-Moulton history,
  ! NO bare Euler; the GW collision is a separate VG-style substep.  gi/gifix/off/VG
  ! fall through to the AB4 path below (byte-unchanged).
  if (trim(sbe_lg_degen) == 'gicov') then
    call dt_evolve_bloch_lg_gicov(sbe, gs, E, dt, icomm)
    return
  end if

  shift_vector(:) = 0.d0 ! shift_vector is set to be 0

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    sbe%qnm(ib, jb, ik) = sbe%qnm_new(ib, jb, ik)
  end do
  end do
  end do

  do ii = 1,am_s-1
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    sbe%dqnm_stock(ib, jb, ik, ii) = sbe%dqnm_stock(ib, jb, ik, ii+1)
  end do
  end do
  end do
  end do

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    sbe%dqnm_stock(ib, jb, ik, am_s) = 0.d0
  end do
  end do
  end do

  qnm_tmp1=0.d0
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    qnm_tmp1(ib, jb, ik) = sbe%qnm(ib, jb, ik)
  end do
  end do
  end do

  call comm_summation(qnm_tmp1, qnm_tmp, nb*nb*nk, icomm)
  call grad_k_array_nb2d_dcomplex(nb,nk,gs%b_matrix,qnm_tmp,sbe%grad_qnm)

  ! LG-SBE diagnostic knockout (1): suppress the intraband k-gradient term
  if (sbe_lg_diag == 1 .or. sbe_lg_diag == 3) sbe%grad_qnm = 0.d0

  abs_E = sqrt(E(1)**2+E(2)**2+E(3)**2)
  proj_E = E(1)*epdir_re1(1) + E(2)*epdir_re1(2) + E(3)*epdir_re1(3) 

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    if(ib == jb) then
    ! qnn (diagonal part)
      if(abs_E >= 1.d-13) then
        sbe%dqnm_stock(ib, ib, ik, am_s) = E(1) * sbe%grad_qnm(ib, ib, 1, ik) + &
                                           E(2) * sbe%grad_qnm(ib, ib, 2, ik) + &
                                           E(3) * sbe%grad_qnm(ib, ib, 3, ik)
        do lb = 1,nb
          if(ib /= lb)then
            sbe%dqnm_stock(ib, ib, ik, am_s) = &
              & sbe%dqnm_stock(ib, ib, ik, am_s) + &
              & 2.d0 * proj_E * sbe%abs_dnm(ib, lb, ik) * &
              &        aimag(sbe%qnm(lb, ib, ik))
          end if
        end do
      end if
    else
    ! qnm (off-diagonal part)
      if (yn_sbe_gw_collision == 'y' .and. trim(sbe_deph_mode) == 'gw') then
        sbe%dqnm_stock(ib, jb, ik, am_s) = 0.d0      ! t_2 replaced by GW (added below)
      else
        sbe%dqnm_stock(ib, jb, ik, am_s) = - sbe%qnm(ib, jb, ik)/t_2
      end if
      if(abs_E >= 1.d-13) then
        sbe%dqnm_stock(ib, jb, ik, am_s) = &
          & sbe%dqnm_stock(ib, jb, ik, am_s) + &
          & E(1) * sbe%grad_qnm(ib, jb, 1, ik) + &
          & E(2) * sbe%grad_qnm(ib, jb, 2, ik) + &
          & E(3) * sbe%grad_qnm(ib, jb, 3, ik)
        sbe%dqnm_stock(ib, jb, ik, am_s) = &
          & sbe%dqnm_stock(ib, jb, ik, am_s) - &
          & zi * gs%delta_omega(ib, jb, ik) * sbe%qnm(ib, jb, ik)
        sbe%dqnm_stock(ib, jb, ik, am_s) = &
          & sbe%dqnm_stock(ib, jb, ik, am_s) - &
          & zi * (E(1) * shift_vector(1) + &
                  E(2) * shift_vector(2) + &
                  E(3) * shift_vector(3)) * sbe%qnm(ib, jb, ik)
        sbe%dqnm_stock(ib, jb, ik, am_s) = &
          & sbe%dqnm_stock(ib, jb, ik, am_s) + &
          & zi * proj_E * sbe%abs_dnm(ib, jb, ik) * &
          &     (sbe%qnm(ib, ib, ik) - sbe%qnm(jb, jb, ik))
        do lb = 1,nb
          if(lb /= ib .and. lb /= jb) then
            sbe%dqnm_stock(ib, jb, ik, am_s) = &
              & sbe%dqnm_stock(ib, jb, ik, am_s) - &
              & zi * proj_E * &
              &     (sbe%abs_dnm(ib, lb, ik) * &
              &      sbe%qnm(lb, jb, ik) - &
              &      sbe%abs_dnm(lb, jb, ik) * &
              &      sbe%qnm(ib, lb, ik)) * &
              &      sbe%exp_iphi(ib, lb, ik) * &
              &      sbe%exp_iphi(lb, jb, ik) * &
              &      conjg(sbe%exp_iphi(ib, jb, ik))
          end if
        end do
      end if
    end if
  end do
  end do
  end do

  ! GW collision term (Phase 2): diagonal QP-lifetime loss + off-diagonal
  ! dephasing into the same Adams-Moulton RHS slot.  OFF path untouched.
  if (yn_sbe_gw_collision == 'y') then
    call add_collision_diag(sbe%dqnm_stock(:, :, :, am_s), sbe%qnm, &
      & gs%gamma_gw, gs%f0_ref, sbe%nb, sbe%nk, sbe%ik_min, sbe%ik_max)
    call add_collision_offdiag(sbe%dqnm_stock(:, :, :, am_s), sbe%qnm, &
      & gs%gamma_gw, sbe%nb, sbe%nk, sbe%ik_min, sbe%ik_max, sbe_deph_mode)
  end if

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    sbe%qnm_new(ib, jb, ik) = sbe%qnm(ib, jb, ik)
  end do
  end do
  end do

  do ii = 1,am_s
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
      sbe%qnm_new(ib, jb, ik) = sbe%qnm_new(ib, jb, ik) + &
       &  bj_am(am_s+1-ii,am_s) * sbe%dqnm_stock(ib, jb, ik, ii) * dt
  end do
  end do
  end do
  end do

end subroutine dt_evolve_bloch_lg

!===================================================================
! gicov integrator (Phase 3, Task 5b).  Steps the PHYSICAL density rho one dt with
! a self-starting classical RK4 built on the (Task 5a machine-exact) gicov_rhs:
!
!   rho0 = rho_ij_from_q(qnm)                        (physical rho at step start)
!   k1 = f(rho0);  k2 = f(rho0 + dt/2 k1)
!   k3 = f(rho0 + dt/2 k2);  k4 = f(rho0 + dt k3)     f == gicov_rhs
!   rho_new = rho0 + dt/6 (k1 + 2 k2 + 2 k3 + k4)
!
! RK4's imaginary-axis stability limit 2*sqrt(2)~2.83 (vs the legacy AB4's 0.43) is
! the R-3 choice; RK4 is SELF-STARTING so there is NO dqnm_stock / Adams-Moulton
! history and NO bare Euler.  gicov_rhs reads sbe%qnm, reconstructs physical rho
! via rho_ij_from_q, and gathers neighbouring-k slices internally (comm_summation),
! so each stage writes its intermediate physical rho back into sbe%qnm through the
! R-1 bridge q_ij_from_rho before its gicov_rhs call.  Because prepare_qnm made
! exp_iphi a time-constant UNIT phase (same-block exp_iphi=1), the qnm<->rho bridge
! is lossless, so stepping rho through qnm is exact.
!
! The GW collision is applied as a SEPARATE substep on the physical rho AFTER the
! coherent RK4 (mirrors dt_evolve_bloch's post-Taylor add_collision_vg); it is a
! no-op when yn_sbe_gw_collision/='y' (e.g. the U3 gate: t_2=1e30, GW off).
!===================================================================
subroutine dt_evolve_bloch_lg_gicov(sbe, gs, E, dt, icomm)
  use salmon_global, only: yn_sbe_gw_collision, sbe_deph_mode
  use sbe_collision_gw, only: add_collision_vg
  implicit none
  type(s_sbe_bloch_solver), intent(inout) :: sbe
  type(s_sbe_gs_info), intent(inout) :: gs
  real(8), intent(in) :: E(1:3)
  real(8), intent(in) :: dt
  integer, intent(in) :: icomm
  integer :: nb, nk, ik, ib, jb
  complex(8), allocatable :: rho0(:,:,:), rho_new(:,:,:)
  complex(8), allocatable :: k1(:,:,:), k2(:,:,:), k3(:,:,:), k4(:,:,:)

  nb = sbe%nb
  nk = sbe%nk

  allocate(rho0   (nb, nb, sbe%ik_min:sbe%ik_max))
  allocate(rho_new(nb, nb, sbe%ik_min:sbe%ik_max))
  allocate(k1(nb, nb, nk), k2(nb, nb, nk), k3(nb, nb, nk), k4(nb, nb, nk))

  ! (0) advance the propagated state (previous step's qnm_new) into qnm; then
  !     reconstruct the physical rho0 at step start on the local k-slice.
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1, nb
  do jb = 1, nb
    sbe%qnm(ib, jb, ik) = sbe%qnm_new(ib, jb, ik)
  end do
  end do
  end do
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1, nb
  do jb = 1, nb
    rho0(ib, jb, ik) = rho_ij_from_q(sbe, ib, jb, ik)
  end do
  end do
  end do

  ! (1) RK4 stage 1 (sbe%qnm already holds rho0's representation)
  call gicov_rhs(sbe, gs, E, k1, icomm)
  ! (2) stage 2 at rho0 + (dt/2) k1
  call load_stage_qnm(sbe, rho0, k1, 0.5d0 * dt)
  call gicov_rhs(sbe, gs, E, k2, icomm)
  ! (3) stage 3 at rho0 + (dt/2) k2
  call load_stage_qnm(sbe, rho0, k2, 0.5d0 * dt)
  call gicov_rhs(sbe, gs, E, k3, icomm)
  ! (4) stage 4 at rho0 + dt k3
  call load_stage_qnm(sbe, rho0, k3, dt)
  call gicov_rhs(sbe, gs, E, k4, icomm)

  ! (5) combine on the local slice
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1, nb
  do jb = 1, nb
    rho_new(ib, jb, ik) = rho0(ib, jb, ik) + (dt / 6.d0) * &
      & (k1(ib, jb, ik) + 2.d0 * k2(ib, jb, ik) + 2.d0 * k3(ib, jb, ik) + k4(ib, jb, ik))
  end do
  end do
  end do

  ! (6) GW collision as a SEPARATE substep on the physical rho (VG-style).
  if (yn_sbe_gw_collision == 'y') then
    call add_collision_vg(rho_new, gs%gamma_gw, gs%f0_ref, dt, &
      & nb, nk, sbe%ik_min, sbe%ik_max, sbe_deph_mode)
  end if

  ! (7) write stepped physical rho into qnm_new; restore qnm to the step-start
  !     representation (AB4 leaves qnm holding the step-start state).
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1, nb
  do jb = 1, nb
    sbe%qnm_new(ib, jb, ik) = q_ij_from_rho(sbe, rho_new(ib, jb, ik), ib, jb, ik)
    sbe%qnm(ib, jb, ik)     = q_ij_from_rho(sbe, rho0(ib, jb, ik),    ib, jb, ik)
  end do
  end do
  end do

  deallocate(rho0, rho_new, k1, k2, k3, k4)

contains

  ! Write the RK stage density rho0 + a*kX into sbe%qnm (local slice) through the
  ! R-1 bridge so the next gicov_rhs sees the intermediate physical rho.
  subroutine load_stage_qnm(sbe, rho0, kX, a)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    complex(8), intent(in) :: rho0(nb, nb, sbe%ik_min:sbe%ik_max)
    complex(8), intent(in) :: kX(nb, nb, nk)
    real(8), intent(in) :: a
    integer :: ik, ib, jb
    !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
    do ik = sbe%ik_min, sbe%ik_max
    do ib = 1, nb
    do jb = 1, nb
      sbe%qnm(ib, jb, ik) = q_ij_from_rho(sbe, rho0(ib, jb, ik) + a * kX(ib, jb, ik), ib, jb, ik)
    end do
    end do
    end do
  end subroutine load_stage_qnm

end subroutine dt_evolve_bloch_lg_gicov

subroutine calc_current_bloch_lg(sbe, gs, jmat, icomm)
    use salmon_global, only: sbe_lg_degen
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(out) :: jmat(3)
    integer, intent(in) :: icomm
    integer :: ik, ib, jb, jj
    integer :: nk, nb
    complex(8) :: tmp1(3),tmp(3)
    complex(8) :: rho_ji

    nk = sbe%nk
    nb = sbe%nb

    tmp1(:) = 0.d0

    ! Pb4 (GI current): when sbe_lg_degen=='gi'/'gifix'/'gicov' the Pb3 xi substitution
    ! (gi/gifix) or the same-block exp_iphi=1 / covariant-transport representation (gicov)
    ! makes the dnm_i / exp_iphi phase machinery inconsistent inside degenerate blocks, so
    ! the delta_omega*|dnm|*aimag(qnm) current below is unreliable there.  Instead reconstruct
    ! the physical density matrix rho from the propagated qnm (rho = exp_iphi*qnm off the
    ! diagonal, rho = qnm on the diagonal; inverse of prepare_qnm's qnm=conjg(exp_iphi)*rho,
    ! == rho_ij_from_q) and contract it with the velocity v = p_tm_matrix (+ rvnl_tm_matrix
    ! when flag_vnl_correction), using the exact index/sign convention of the VG current
    ! calc_current_bloch (paramagnetic part; the length gauge carries no Ac*N term).  This
    ! reconstruction is Tr(v rho), gauge-invariant under a consistent GS gauge, so it is the
    ! minimal correct gicov current (subsumes Task 6's in-block current reconstruction).
    if (trim(sbe_lg_degen) == 'gi' .or. trim(sbe_lg_degen) == 'gifix' &
        & .or. trim(sbe_lg_degen) == 'gicov') then
        !$omp parallel do default(shared) private(ik,jj,ib,jb,rho_ji) reduction(+:tmp1)
        do ik = sbe%ik_min, sbe%ik_max
        do jj = 1,3
        do ib = 1,nb
        do jb = 1,nb
            if (ib == jb) then
                rho_ji = sbe%qnm_new(ib, ib, ik)
            else
                rho_ji = sbe%exp_iphi(jb, ib, ik) * sbe%qnm_new(jb, ib, ik)
            end if
            tmp1(jj) = tmp1(jj) + gs%kweight(ik) * rho_ji * gs%p_tm_matrix(ib, jb, jj, ik)
            if (sbe%flag_vnl_correction) then
                tmp1(jj) = tmp1(jj) + gs%kweight(ik) * rho_ji * gs%rvnl_tm_matrix(ib, jb, jj, ik)
            end if
        end do
        end do
        end do
        end do
        !$omp end parallel do

        call comm_summation(tmp1, tmp, 3, icomm)
        jmat(:) = dble(tmp(:)) / sum(gs%kweight(:)) / gs%volume
        return
    end if

    !$omp parallel do default(shared) private(ik,ib,jb,jj) reduction(+:tmp1)
    do ik = sbe%ik_min, sbe%ik_max
    do jj = 1,3
    do ib = 1,nb
    do jb = 1,nb
        if(ib == jb) then
            tmp1(jj) = tmp1(jj) + gs%grad_k_eigen(ib, jj, ik) *  &
                               & dble(sbe%qnm_new(ib, ib, ik)) * gs%kweight(ik)
        else
            tmp1(jj) = tmp1(jj) + gs%delta_omega(ib, jb, ik) * &
                               & abs(sbe%dnm_i(ib, jb, jj, ik)) * &
                               & aimag(sbe%qnm_new(jb, ib, ik)) * gs%kweight(ik)
        end if
    end do
    end do
    end do
    end do

    call comm_summation(tmp1, tmp, 3, icomm)
    jmat(:) = -dble(tmp(:)) / sum(gs%kweight(:)) / gs%volume

end subroutine calc_current_bloch_lg

subroutine adams_moulton_coefs(bj_am)
  implicit none
  real(8) :: bj_am(8,8)

  bj_am(1,1) = 1.d0

  !bj_am(1,2) = 3.d0/2.d0
  !bj_am(2,2) = -1.d0/2.d0

  !bj_am(1,3) = 23.d0/12.d0
  !bj_am(2,3) = -4.d0/3.d0
  !bj_am(3,3) = 5.d0/12.d0

  bj_am(1,4) = 55.d0/24.d0
  bj_am(2,4) = -59.d0/24.d0
  bj_am(3,4) = 37.d0/24.d0
  bj_am(4,4) = -3.d0/8.d0

  !bj_am(1,5) = 1901.d0/720.d0
  !bj_am(2,5) = -1387.d0/360.d0
  !bj_am(3,5) = 109.d0/30.d0
  !bj_am(4,5) = -637.d0/360.d0
  !bj_am(5,5) = 251.d0/720.d0

  !bj_am(1,6) = 4277.d0/1440.d0
  !bj_am(2,6) = -2641.d0/480.d0
  !bj_am(3,6) = 4991.d0/720.d0
  !bj_am(4,6) = -3649.d0/720.d0
  !bj_am(5,6) = 959.d0/480.d0
  !bj_am(6,6) = -95.d0/288.d0

  !bj_am(1,7) = 198721.d0/60480.d0
  !bj_am(2,7) = -18637.d0/2520.d0
  !bj_am(3,7) = 235183.d0/20160.d0
  !bj_am(4,7) = -10754.d0/945.d0
  !bj_am(5,7) = 135713.d0/20160.d0
  !bj_am(6,7) = -5603.d0/2520.d0
  !bj_am(7,7) = 19087.d0/60480.d0

  bj_am(1,8) = 16083.d0/4480.d0
  bj_am(2,8) = -1152169.d0/120960.d0
  bj_am(3,8) = 242653.d0/13440.d0
  bj_am(4,8) = -296053.d0/13440.d0
  bj_am(5,8) = 2102243.d0/120960.d0
  bj_am(6,8) = -115747.d0/13440.d0
  bj_am(7,8) = 32863.d0/13440.d0
  bj_am(8,8) = -5257.d0/17280.d0

end subroutine adams_moulton_coefs


end module

