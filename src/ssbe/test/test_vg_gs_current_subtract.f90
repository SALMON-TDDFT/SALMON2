! src/ssbe/test/test_vg_gs_current_subtract.f90
! Velocity-gauge frozen-ground-state current subtraction
! (yn_sbe_gs_current_subtract): analytic-fixture gate.
!
! WHY THIS FILE EXISTS
!   Truncating the VG-SBE basis at nstate_sbe breaks the f-sum rule, so the
!   readout J = Tr[rho(t) v(k+A)] carries a spurious current proportional to
!   A(t) (a filled-band insulator in a COMPLETE basis carries none; observed as
!   the non-converging N_b-ladder of peak |J| in the gw campaign).  The
!   standard remedy (upstream SSBE practice) subtracts the frozen initial
!   density matrix evaluated with the SAME velocity operator at the same time:
!
!     J(t) = Tr[rho(t) v(k+A(t))] - Tr[rho0 v(k+A(t))],  rho0 = diag(occup).
!
!   calc_current_bloch implements this by calling the identical evaluation
!   routine (calc_current_bloch_core) a second time on rho0, so the definition
!   of v (bare p_tm, A*trace diamagnetic term, yn_vnl_correction,
!   norder_correction>=1 A-corrections) is shared by construction.  This test
!   pins the three defining properties on a small analytic fixture plus one
!   independent analytic evaluation of the legacy formula:
!
!     0. OFF path unchanged: flag='n' reproduces an INDEPENDENTLY coded
!        evaluation of the legacy norder=0 formula
!        J = (Re sum_k w tr[rho p] / sum w + A * tr[rho]_w) / V   (refactor
!        witness: the wrapper/core split did not move the legacy math), and
!        the frozen term itself is analytic: J(rho0, norder=0) = A*Ne/V
!        (diagonal p contribution cancels exactly on the TR-symmetric grid).
!     a. rho == rho0, flag ON  =>  J == 0 for EVERY A (all norder, vnl on/off):
!        the subtrahend is exactly the frozen-state current.
!     b. generic rho (perturbed populations + coherences):
!        J_on(rho) == J_off(rho) - J_off(rho0) for every A/norder/vnl combo,
!        where J_off(rho0) is obtained through the PUBLIC path
!        (init_sbe_bloch_solver => rho = diag(occup), flag off).  Non-vacuous:
!        the subtracted term is required to be well above roundoff for A/=0.
!     c. A == 0  =>  the subtrahend vanishes (time-reversal-symmetric k-grid:
!        p_bb(-k) = -p_bb(k) cancels the trace pairwise, and every norder
!        correction carries a factor of A), so J_on == J_off at A=0 to
!        accumulation roundoff (~1e-19 here; exact zero in exact arithmetic):
!        the post-pulse observables are untouched by the flag.
!
!   Fixture: nb=4, nk=2 with k2 the exact time-reversal partner of k1
!   (p(k2) = -conjg(p(k1)), same eigenvalues) -- the smallest grid on which the
!   physical A=0 statement (c) holds exactly, as it does for any real
!   materials' full-BZ k-grid.  Occupied bands 1..2 (occup=2, ne=4).  All gaps
!   >= 0.15 >> the 1e-3 delta_omega floor, so every norder>=1 correction term
!   is engaged (nothing is silently gated off).  Nonzero band-diagonal p AND
!   nonzero rvnl exercise the intraband + vnl branches.
!
! BUILD (already-built ninja tree at build_local/; single-process communication
! dummy).  Compile from a CLEAN dir (a stale .mod file in the repo root shadows
! the fresh build_local one -- use -I build_local ONLY).  Link the SAME objects
! the salmon executable built, minus main.f90.o, via an @objs.txt response
! file:
!
!   find <repo>/build_local/src/CMakeFiles/salmon.dir -name '*.o' ! -name 'main.f90.o' > objs.txt
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     -I<repo>/build_local -J<clean_dir> \
!     -c <repo>/src/ssbe/test/test_vg_gs_current_subtract.f90 -o <clean_dir>/test_vg_gs_current_subtract.o
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     @objs.txt <clean_dir>/test_vg_gs_current_subtract.o -o <clean_dir>/test_vg_gs_current_subtract \
!     -framework Accelerate -lm -ldl
!   <clean_dir>/test_vg_gs_current_subtract
!
program test_vg_gs_current_subtract
  use gs_info_ssbe,      only: s_sbe_gs_info
  use bloch_solver_ssbe, only: s_sbe_bloch_solver, init_sbe_bloch_solver, &
                               calc_current_bloch
  use salmon_global,     only: gauge_sbe, norder_correction, &
                               yn_sbe_gs_current_subtract
  implicit none

  integer, parameter :: nb = 4, nk = 2, ne = 4
  complex(8), parameter :: zi_ = (0d0, 1d0)
  real(8), parameter :: vol = 150.0d0
  ! A-sweep: zero, small/large along z, oblique, negative -- the subtraction
  ! must hold pointwise in A, not just for one field.
  integer, parameter :: nA = 5
  real(8), parameter :: Aset(3, nA) = reshape( (/ &
    &  0.00d0, 0.00d0, 0.00d0, &
    &  0.00d0, 0.00d0, 0.05d0, &
    &  0.00d0, 0.00d0, 0.80d0, &
    &  0.30d0, -0.20d0, 0.55d0, &
    & -0.15d0, 0.10d0, -0.40d0 /), (/3, nA/) )
  integer :: nfail

  nfail = 0
  call set_globals()

  call test_off_path_analytic(nfail)          ! check 0 (refactor witness)
  call test_rho0_gives_zero(nfail)            ! check a
  call test_difference_is_frozen_term(nfail)  ! check b (+ non-vacuous)
  call test_a0_identity(nfail)                ! check c

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_vg_gs_current_subtract checks passed."
  end if

contains

  !======================= salmon_global fixture ==============================
  subroutine set_globals()
    implicit none
    gauge_sbe = 'velocity_gauge'      ! calc_current_bloch is the VG readout
    norder_correction = 0             ! per-check loops override
    yn_sbe_gs_current_subtract = 'n'  ! per-check code overrides
  end subroutine set_globals

  !======================= assert helper (test_gicov_* style) =================
  subroutine check_true(cond, label, nfail)
    implicit none
    logical, intent(in) :: cond
    character(*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (cond) then
      write(*, '(a)') "ok    " // label
    else
      write(*, '(a)') "FAIL  " // label
      nfail = nfail + 1
    end if
  end subroutine check_true

  !======================= gs fixture =========================================
  ! nk=2 exact time-reversal pair: eigen(k2)=eigen(k1), p(k2)=-conjg(p(k1)),
  ! rvnl likewise (both Hermitian per k).  Bands 1..2 occupied (occup=2).
  subroutine build_gs(gs)
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    integer :: ib, jb, idir, ik
    real(8) :: re, im
    complex(8) :: h(nb, nb)

    gs%nk = nk
    gs%nb = nb
    gs%ne = ne
    gs%volume = vol
    allocate(gs%kweight(nk), gs%eigen(nb, nk), gs%occup(nb, nk))
    allocate(gs%delta_omega(nb, nb, nk))
    allocate(gs%p_tm_matrix(nb, nb, 3, nk), gs%rvnl_tm_matrix(nb, nb, 3, nk))
    gs%kweight(:) = 1.0d0

    gs%eigen(1, :) = -0.30d0
    gs%eigen(2, :) = -0.10d0
    gs%eigen(3, :) =  0.25d0
    gs%eigen(4, :) =  0.60d0
    do ik = 1, nk
      do jb = 1, nb
        do ib = 1, nb
          gs%delta_omega(ib, jb, ik) = gs%eigen(ib, ik) - gs%eigen(jb, ik)
        end do
      end do
    end do

    gs%occup(:, :) = 0d0
    gs%occup(1:ne/2, :) = 2d0

    ! Hermitian p at k1 from a deterministic dense seed (nonzero band-diagonal
    ! entries included -- the intraband channel and the TR cancellation at A=0
    ! must both be exercised); k2 = TR partner.
    do idir = 1, 3
      do jb = 1, nb
        do ib = 1, nb
          re = 0.11d0 * dble(ib) - 0.07d0 * dble(jb) + 0.05d0 * dble(idir) &
             & + 0.021d0 * dble(ib * jb)
          im = 0.013d0 * dble(ib - jb) * dble(idir)
          h(ib, jb) = dcmplx(re, im)
        end do
      end do
      gs%p_tm_matrix(:, :, idir, 1) = 0.5d0 * (h + transpose(conjg(h)))
      gs%p_tm_matrix(:, :, idir, 2) = -conjg(gs%p_tm_matrix(:, :, idir, 1))
    end do

    ! Hermitian rvnl (nonlocal-velocity correction), smaller magnitude,
    ! same TR structure.
    do idir = 1, 3
      do jb = 1, nb
        do ib = 1, nb
          re = 0.017d0 * dble(ib + jb) - 0.009d0 * dble(idir) &
             & + 0.004d0 * dble(ib * jb)
          im = 0.006d0 * dble(ib - jb) * dble(4 - idir)
          h(ib, jb) = dcmplx(re, im)
        end do
      end do
      gs%rvnl_tm_matrix(:, :, idir, 1) = 0.5d0 * (h + transpose(conjg(h)))
      gs%rvnl_tm_matrix(:, :, idir, 2) = -conjg(gs%rvnl_tm_matrix(:, :, idir, 1))
    end do
  end subroutine build_gs

  !======================= solver fixtures ====================================
  ! Ground-state solver through the PUBLIC path: init_sbe_bloch_solver sets
  ! rho = diag(gs%occup) -- the same rho0 the subtraction must reproduce.
  subroutine make_solver_gs(sbe, gs, vnl)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    logical, intent(in) :: vnl
    integer :: icomm
    icomm = 0
    call init_sbe_bloch_solver(sbe, gs, nb, icomm)
    sbe%flag_vnl_correction = vnl
  end subroutine make_solver_gs

  ! Generic excited-state density: Hermitian, populations perturbed away from
  ! occup (trace NOT preserved on purpose -- the A*trace terms of minuend and
  ! subtrahend must NOT be assumed to cancel) plus dense coherences, different
  ! at the two k (no TR symmetry imposed on rho: none is required).
  subroutine make_solver_excited(sbe, gs, vnl)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    logical, intent(in) :: vnl
    integer :: ik, ib, jb
    complex(8) :: c
    call make_solver_gs(sbe, gs, vnl)
    do ik = sbe%ik_min, sbe%ik_max
      sbe%rho(1, 1, ik) = 1.82d0 - 0.11d0 * dble(ik)
      sbe%rho(2, 2, ik) = 1.63d0 + 0.07d0 * dble(ik)
      sbe%rho(3, 3, ik) = 0.31d0 + 0.05d0 * dble(ik)
      sbe%rho(4, 4, ik) = 0.12d0 - 0.02d0 * dble(ik)
      do jb = 1, nb
        do ib = jb + 1, nb
          c = dcmplx(0.13d0 * dble(ib) - 0.05d0 * dble(jb) + 0.02d0 * dble(ik), &
                   & 0.09d0 * dble(ib - jb) - 0.03d0 * dble(ik))
          sbe%rho(ib, jb, ik) = c
          sbe%rho(jb, ib, ik) = conjg(c)
        end do
      end do
    end do
  end subroutine make_solver_excited

  !======================= independent legacy formula (norder=0) ==============
  ! Deliberately re-coded from the physics (NOT from calc_current_bloch_core):
  ! J = ( Re sum_k w tr[rho p] / sum_k w  +  A * (sum_k w tr[rho] / sum_k w) ) / V
  subroutine current_reference_norder0(sbe, gs, Ac, jref)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(in) :: Ac(3)
    real(8), intent(out) :: jref(3)
    integer :: ik, ib, jb, idir
    complex(8) :: s(3), pmat
    real(8) :: tr
    s(:) = (0d0, 0d0)
    tr = 0d0
    do ik = sbe%ik_min, sbe%ik_max
      do ib = 1, nb
        tr = tr + dble(sbe%rho(ib, ib, ik)) * gs%kweight(ik)
        do jb = 1, nb
          do idir = 1, 3
            pmat = gs%p_tm_matrix(ib, jb, idir, ik)
            if (sbe%flag_vnl_correction) pmat = pmat + gs%rvnl_tm_matrix(ib, jb, idir, ik)
            s(idir) = s(idir) + gs%kweight(ik) * sbe%rho(jb, ib, ik) * pmat
          end do
        end do
      end do
    end do
    jref(:) = (dble(s(:)) / sum(gs%kweight) + Ac(:) * (tr / sum(gs%kweight))) / vol
  end subroutine current_reference_norder0

  !======================= check 0: OFF path analytic =========================
  subroutine test_off_path_analytic(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe, sbe_gs
    real(8) :: jm(3), jref(3), errmax, err0
    integer :: ia, icomm
    logical :: vnl
    integer :: ivnl

    icomm = 0
    call build_gs(gs)
    errmax = 0d0
    norder_correction = 0
    yn_sbe_gs_current_subtract = 'n'
    do ivnl = 0, 1
      vnl = (ivnl == 1)
      call make_solver_excited(sbe, gs, vnl)
      do ia = 1, nA
        call calc_current_bloch(sbe, gs, Aset(:, ia), jm, icomm)
        call current_reference_norder0(sbe, gs, Aset(:, ia), jref)
        errmax = max(errmax, maxval(abs(jm - jref)))
      end do
      deallocate(sbe%rho)
    end do
    call check_true(errmax < 1d-13, "off-path == independent legacy formula (norder=0)", nfail)

    ! frozen-state current is analytic at norder=0: the band-diagonal p sums
    ! cancel exactly on the TR pair, leaving the pure diamagnetic A*Ne/V.
    call make_solver_gs(sbe_gs, gs, .false.)
    err0 = 0d0
    do ia = 1, nA
      call calc_current_bloch(sbe_gs, gs, Aset(:, ia), jm, icomm)
      err0 = max(err0, maxval(abs(jm - Aset(:, ia) * dble(ne) / vol)))
    end do
    deallocate(sbe_gs%rho)
    call check_true(err0 < 1d-14, "frozen-state current == A*Ne/V (norder=0, TR grid)", nfail)
  end subroutine test_off_path_analytic

  !======================= check a: rho==rho0, flag ON => J == 0 ==============
  subroutine test_rho0_gives_zero(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    real(8) :: jm(3), errmax
    integer :: ia, no, icomm, ivnl

    icomm = 0
    call build_gs(gs)
    errmax = 0d0
    yn_sbe_gs_current_subtract = 'y'
    do ivnl = 0, 1
      call make_solver_gs(sbe, gs, ivnl == 1)
      do no = 0, 3
        norder_correction = no
        do ia = 1, nA
          call calc_current_bloch(sbe, gs, Aset(:, ia), jm, icomm)
          errmax = max(errmax, maxval(abs(jm)))
        end do
      end do
      deallocate(sbe%rho)
    end do
    yn_sbe_gs_current_subtract = 'n'
    call check_true(errmax < 1d-14, "rho=rho0, flag on: J == 0 for all A/norder/vnl", nfail)
  end subroutine test_rho0_gives_zero

  !======================= check b: J_on == J_off - J_off(rho0) ===============
  subroutine test_difference_is_frozen_term(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe, sbe_gs
    real(8) :: j_on(3), j_off(3), j_frozen(3), errmax, submin
    integer :: ia, no, icomm, ivnl

    icomm = 0
    call build_gs(gs)
    errmax = 0d0
    submin = huge(1d0)
    do ivnl = 0, 1
      call make_solver_excited(sbe, gs, ivnl == 1)
      call make_solver_gs(sbe_gs, gs, ivnl == 1)
      do no = 0, 3
        norder_correction = no
        do ia = 1, nA
          yn_sbe_gs_current_subtract = 'n'
          call calc_current_bloch(sbe, gs, Aset(:, ia), j_off, icomm)
          call calc_current_bloch(sbe_gs, gs, Aset(:, ia), j_frozen, icomm)
          yn_sbe_gs_current_subtract = 'y'
          call calc_current_bloch(sbe, gs, Aset(:, ia), j_on, icomm)
          errmax = max(errmax, maxval(abs(j_on - (j_off - j_frozen))))
          ! non-vacuous: for A /= 0 the subtracted term must be well above
          ! roundoff (~A*Ne/V >= 0.05*4/150 ~ 1.3e-3), or the identity above
          ! would be the trivial 0 == 0 - 0.
          if (ia > 1) submin = min(submin, maxval(abs(j_frozen)))
        end do
      end do
      deallocate(sbe%rho, sbe_gs%rho)
    end do
    yn_sbe_gs_current_subtract = 'n'
    call check_true(errmax < 1d-14, "flag on/off difference == frozen-state term (public path)", nfail)
    call check_true(submin > 1d-4, "non-vacuous: frozen term well above roundoff for A/=0", nfail)
  end subroutine test_difference_is_frozen_term

  !======================= check c: A=0 => subtrahend == 0 ====================
  subroutine test_a0_identity(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    real(8) :: j_on(3), j_off(3), a0(3), errmax
    integer :: no, icomm, ivnl

    icomm = 0
    a0(:) = 0d0
    call build_gs(gs)
    errmax = 0d0
    do ivnl = 0, 1
      call make_solver_excited(sbe, gs, ivnl == 1)
      do no = 0, 3
        norder_correction = no
        yn_sbe_gs_current_subtract = 'n'
        call calc_current_bloch(sbe, gs, a0, j_off, icomm)
        yn_sbe_gs_current_subtract = 'y'
        call calc_current_bloch(sbe, gs, a0, j_on, icomm)
        errmax = max(errmax, maxval(abs(j_on - j_off)))
      end do
      deallocate(sbe%rho)
    end do
    yn_sbe_gs_current_subtract = 'n'
    ! TR pairwise cancellation is exact term-by-term, but the sequential
    ! accumulation order leaves O(ulp) residue (~1e-19 at this scale, J itself
    ! ~1e-3): assert roundoff-zero, far below any physical signal.
    call check_true(errmax < 1d-15, "A=0: subtrahend vanishes to roundoff (post-pulse untouched)", nfail)
  end subroutine test_a0_identity

end program test_vg_gs_current_subtract
