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



