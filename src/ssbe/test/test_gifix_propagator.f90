! src/ssbe/test/test_gifix_propagator.f90
! LG-SBE Tier2 Task 4: multi-step propagator conservation test (production path).
!
! WHY THIS FILE EXISTS (guards the Pc-test false positive)
!   test_gisbe_covariance.f90 (Tier2 Pc, commit 11015404) asserts that the
!   non-Abelian connection xi is a gauge-COVARIANT tensor for a SINGLE
!   application of [xi,rho] -- 30/30 checks passed. Yet production still
!   diverged (see memory project_gw_lgsbe_divergence: GI-LG step20 blow-up).
!   A 1-step covariance check cannot see that failure because it never drives
!   the actual multistep Adams-Bashforth propagator; it only checks an
!   algebraic identity at t=0.
!
!   This file closes that gap: it builds a minimal fixed-block system (nb=4,
!   nk=4, doubly-degenerate block {1,2}) and runs it through the REAL
!   production propagator -- prepare_qnm + dt_evolve_bloch_lg in
!   bloch_solver_ssbe.f90, called in the exact order used by
!   realtime_ssbe.f90:44-105 (init_sbe_bloch_solver -> prepare_qnm ->
!   adams_moulton_coefs -> repeated dt_evolve_bloch_lg) -- for >=200 steps,
!   and asserts trace(qnm) and Hermiticity stay conserved throughout.
!
!   Positive (test_gifix_conserves): the {1,2} interband coupling is built
!   from the production xi kernel (degenerate_block_ssbe:xi_block_from_overlap,
!   already covariance-tested by Pc) on a smooth overlap -- an authentic,
!   bounded (O(1)) Hermitian non-Abelian connection, exactly what
!   sbe_lg_degen='gifix' places into gs%d_matrix for a fixed block. Must stay
!   bounded: trace conserved to 1e-10, Hermiticity norm <= 1e-8, every step.
!
!   Negative / teeth (test_gi_naive_diverges): the SAME system, but the
!   {1,2} coupling is instead the naive per-k dipole d = i*p/delta_omega --
!   what a per-k 'gi' block whose xi construction is rejected falls back to.
!   codex-r3: delta_omega must be NEAR-degenerate (2d-9: above the 1d-9
!   floor, below theta_on=5d-4), NOT exactly degenerate -- at exact
!   degeneracy prepare_qnm's own abs_dnm>=1d-13 gate (bloch_solver_ssbe.f90)
!   and the naive-dipole knockout would zero the coupling and the test would
!   spuriously pass. With a nonzero interband momentum p, |d| ~ |p|/domega
!   ~ 1.5e8 (matches the plan's estimate) is singular: dt*Omega_Rabi ~
!   0.05 * 2*(1e-6)*(1.5e8) ~ 16, far outside the explicit AB4 stability
!   region, so the SAME real propagator MUST blow up -- proving this test
!   has teeth (it would fail on the pre-xi-fix code path).
!
! BUILD (already-built CMake/ninja tree at build_local/; no MPI configured,
! so `communication` resolves to the single-process communication_dummy.f90
! stub -- comm_get_groupinfo always reports rank 0 of 1, comm_summation is a
! local copy). bloch_solver_ssbe.f90 is NOT standalone (unlike
! degenerate_block_ssbe.f90) -- it needs salmon_global, sbe_collision_gw,
! common_ssbe, gs_info_ssbe, so this test links against the SAME objects the
! `salmon` executable already built, minus main.f90.o (which supplies its own
! `program main` and would collide with this file's `program`):
!
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     -I<repo>/src/ssbe -I<repo>/build_local -J<scratch_dir> \
!     -c <repo>/src/ssbe/test/test_gifix_propagator.f90 -o <scratch_dir>/test_gifix_propagator.o
!
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     -L/opt/homebrew/opt/mysql-client/lib \
!     $(find <repo>/build_local/src/CMakeFiles/salmon.dir -name '*.o' ! -name 'main.f90.o') \
!     <scratch_dir>/test_gifix_propagator.o \
!     -o <scratch_dir>/test_gifix_propagator \
!     -framework Accelerate -lm -ldl
!
!   <scratch_dir>/test_gifix_propagator
!
program test_gifix_propagator
  use gs_info_ssbe,          only: s_sbe_gs_info
  use bloch_solver_ssbe,     only: s_sbe_bloch_solver, init_sbe_bloch_solver, &
                                    prepare_qnm, dt_evolve_bloch_lg, adams_moulton_coefs
  use degenerate_block_ssbe, only: xi_block_from_overlap, xi_sign
  use salmon_global,         only: epdir_re1, am_s, sbe_lg_diag, num_kgrid, t_2, &
                                    yn_sbe_gw_collision, sbe_deph_mode
  implicit none

  integer, parameter :: nb = 4, nk = 4
  complex(8), parameter :: zi_ = (0d0, 1d0)
  integer :: nfail

  nfail = 0
  call set_global_defaults()

  call test_gifix_conserves(nfail)
  call test_gi_naive_diverges(nfail)

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_gifix_propagator checks passed."
  end if

contains

  !======================= salmon_global fixture ==============================
  subroutine set_global_defaults()
    implicit none
    epdir_re1(1) = 1d0; epdir_re1(2) = 0d0; epdir_re1(3) = 0d0
    am_s = 4                          ! matches adams_moulton_coefs' populated column (:,4)
    sbe_lg_diag = 0                   ! no diagnostic knockouts: full production physics
    num_kgrid(1) = nk; num_kgrid(2) = 1; num_kgrid(3) = 1
    t_2 = 1d6                         ! long dephasing time -> negligible over the run
    yn_sbe_gw_collision = 'n'
    sbe_deph_mode = ''
  end subroutine set_global_defaults

  !======================= minimal fixed-block gs fixture ======================
  ! nb=4, nk=4; a doubly-degenerate pair {1,2} at delta_omega = +-domega12
  ! (k-independent: same physics at every k, so this is a genuinely FIXED
  ! block); bands {3,4} are decoupled non-degenerate spectators (d_matrix=0
  ! for every pair touching 3 or 4, so their populations/coherences are
  ! inert and cannot mask a {1,2} instability). coupling12 is placed at
  ! d_matrix(1,2,1,:) / d_matrix(2,1,1,:) = conjg(coupling12) (Hermitian),
  ! axis 1 only, matching epdir_re1=(1,0,0).
  !==============================================================================
  subroutine build_gs(gs, domega12, coupling12)
    implicit none
    type(s_sbe_gs_info), intent(out) :: gs
    real(8), intent(in) :: domega12
    complex(8), intent(in) :: coupling12
    real(8), parameter :: two_pi = 6.28318530717958647692d0
    integer :: ik, ib, jb
    real(8) :: eigen(nb)

    gs%nk = nk; gs%nb = nb; gs%ne = 4

    allocate(gs%eigen(nb, nk), gs%occup(nb, nk))
    allocate(gs%delta_omega(nb, nb, nk))
    allocate(gs%d_matrix(nb, nb, 3, nk))

    eigen(1) = 0.30d0 + 0.5d0 * domega12
    eigen(2) = 0.30d0 - 0.5d0 * domega12
    eigen(3) = 1.00d0
    eigen(4) = 1.60d0

    gs%b_matrix = 0d0
    gs%b_matrix(1, 1) = two_pi
    gs%b_matrix(2, 2) = two_pi
    gs%b_matrix(3, 3) = two_pi

    gs%d_matrix = (0d0, 0d0)
    do ik = 1, nk
      gs%eigen(:, ik) = eigen(:)
      ! NOTE: occupation across the {1,2} pair must be ASYMMETRIC (2,0), not
      ! (2,2). With equal occupation qnm11-qnm22==0 identically, which is the
      ! ONLY source term that drives the {1,2} coherence in dt_evolve_bloch_lg
      ! (the "+zi*proj_E*abs_dnm(ib,jb)*(qnm(ib,ib)-qnm(jb,jb))" term) -- an
      ! unbroken-symmetry trap that silently makes the huge naive-dipole
      ! coupling magnitude irrelevant (nothing ever excites it) and produced a
      ! false PASS for the "teeth" test below during development.
      gs%occup(1, ik) = 2d0; gs%occup(2, ik) = 0d0   ! near-degenerate pair, ASYMMETRIC occupation
      gs%occup(3, ik) = 0d0; gs%occup(4, ik) = 0d0   ! decoupled spectator conduction bands
      do jb = 1, nb
        do ib = 1, nb
          gs%delta_omega(ib, jb, ik) = eigen(ib) - eigen(jb)
        end do
      end do
      gs%d_matrix(1, 2, 1, ik) = coupling12
      gs%d_matrix(2, 1, 1, ik) = conjg(coupling12)
    end do
  end subroutine build_gs

  !======================= production xi fixture ===============================
  ! Calls the REAL Tier2 Pb3 kernel (degenerate_block_ssbe:xi_block_from_overlap,
  ! covariance-tested by test_gisbe_covariance.f90) on a smooth 2x2 degenerate
  ! link overlap, giving an authentic bounded Hermitian connection value --
  ! exactly what sbe_lg_degen='gifix' would place in gs%d_matrix for this block.
  !==============================================================================
  function xi_coupling_12() result(c)
    implicit none
    complex(8) :: c
    complex(8) :: Qk(2, 2), Qkb(2, 2), M(2, 2), xi_blk(2, 2)
    real(8) :: ang0, dk, resu
    integer :: info
    ang0 = 0.4d0; dk = 0.1d0
    Qk  = givens2(0.5d0 * ang0)
    Qkb = givens2(0.5d0 * (ang0 + dk))
    M   = matmul(transpose(conjg(Qk)), Qkb)
    call xi_block_from_overlap(M, dk, xi_sign, xi_blk, info, resu)
    if (info /= 0) then
      write(*, '(a,i0)') 'FATAL(test fixture): xi_block_from_overlap rejected the link, info=', info
      stop 1
    end if
    c = xi_blk(1, 2)
  end function xi_coupling_12

  function givens2(th) result(G)
    implicit none
    real(8), intent(in) :: th
    complex(8) :: G(2, 2)
    G(1, 1) = dcmplx(cos(th), 0d0); G(1, 2) = dcmplx(-sin(th), 0d0)
    G(2, 1) = dcmplx(sin(th), 0d0); G(2, 2) = dcmplx(cos(th), 0d0)
  end function givens2

  !======================= inline conservation probes ===========================
  ! trace_re_of / herm_norm_of are NOT existing production helpers (plan's
  ! codex-r1 note) -- inline computations over sbe%qnm_new using ordinary
  ! sbe/gs field access, in the spirit of the Pc test's small utility
  ! functions (hc/comm/mdiff, test_gisbe_covariance.f90:97-123).
  !==============================================================================
  real(8) function trace_re_of(sbe) result(tr)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    integer :: ik, ib
    tr = 0d0
    do ik = sbe%ik_min, sbe%ik_max
      do ib = 1, sbe%nb
        tr = tr + dble(sbe%qnm_new(ib, ib, ik))
      end do
    end do
  end function trace_re_of

  real(8) function herm_norm_of(sbe) result(hn)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    integer :: ik, ib, jb
    hn = 0d0
    do ik = sbe%ik_min, sbe%ik_max
      do jb = 1, sbe%nb
        do ib = 1, sbe%nb
          hn = max(hn, abs(sbe%qnm_new(ib, jb, ik) - conjg(sbe%qnm_new(jb, ib, ik))))
        end do
      end do
    end do
  end function herm_norm_of

  logical function is_finite(x) result(ok)
    implicit none
    real(8), intent(in) :: x
    ok = (x == x) .and. (abs(x) <= huge(1d0))     ! rejects NaN (x/=x) and +-Inf
  end function is_finite

  !======================= assert helper (ssbe style, mirrors ==================
  ! test_gisbe_covariance.f90:60-70's check_true) ================================
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

  !======================= (1) POSITIVE: gifix conserves =======================
  ! Real prepare_qnm + dt_evolve_bloch_lg, >=200 steps, dt=0.05, E just above
  ! the abs_E>=1d-13 gate (bloch_solver_ssbe.f90:662,682) so the xi/dipole
  ! path actually runs. trace must stay within 1e-10 of its t=0 value and
  ! herm_norm must stay <=1e-8 at EVERY step (not just the final one).
  !==============================================================================
  subroutine test_gifix_conserves(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    real(8) :: tr0, E(3), bj_am(8, 8), dt, tr, hn
    integer :: it, icomm
    logical :: trace_ok, herm_ok

    call build_gs(gs, 2d-9, xi_coupling_12())

    E(1) = 1d-6; E(2) = 0d0; E(3) = 0d0    ! small but >= 1e-13: engages the dipole path
    dt = 0.05d0
    icomm = 0                              ! single-process dummy communicator (no MPI build)

    call init_sbe_bloch_solver(sbe, gs, nb, icomm)
    call prepare_qnm(sbe, gs, icomm)
    call adams_moulton_coefs(bj_am)
    tr0 = trace_re_of(sbe)

    trace_ok = .true.; herm_ok = .true.
    do it = 1, 200
      call dt_evolve_bloch_lg(sbe, gs, E, bj_am, dt, icomm)
      tr = trace_re_of(sbe)
      hn = herm_norm_of(sbe)
      if (.not. is_finite(tr) .or. abs(tr - tr0) > 1d-10) trace_ok = .false.
      if (.not. is_finite(hn) .or. hn > 1d-8)             herm_ok  = .false.
    end do

    call check_true(trace_ok, "gifix: trace(qnm) conserved to 1e-10 over 200 real-propagator steps", nfail)
    call check_true(herm_ok,  "gifix: Hermiticity norm <= 1e-8 over 200 real-propagator steps", nfail)
  end subroutine test_gifix_conserves

  !======================= (2) NEGATIVE / TEETH: gi-style naive dipole =========
  ! SAME system, SAME real propagator, SAME E/dt/step count. Only the {1,2}
  ! coupling differs: naive d = i*p/domega with domega=2d-9 (near- but not
  ! exactly degenerate) and a nonzero interband p, giving |d| ~ 1.5e8. This
  ! is what a per-k 'gi' block falls back to when its xi construction is
  ! rejected. If this test PASSED (stayed bounded), it would mean the
  ! conservation checks above have no teeth. It must diverge.
  !==============================================================================
  subroutine test_gi_naive_diverges(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    real(8) :: E(3), bj_am(8, 8), dt, domega12, tr0, tr, hn
    complex(8) :: p12, d_naive
    integer :: it, icomm, it_diverged
    logical :: diverged

    domega12 = 2d-9                        ! above the 1d-9 floor, below theta_on=5d-4
    p12 = dcmplx(0.3d0, 0.1d0)             ! nonzero interband momentum matrix element
    d_naive = zi_ * p12 / domega12          ! naive dipole i*p/delta_omega
    call check_true(abs(d_naive) > 1d7, &
      "gi negative-control fixture: |i*p/domega| >> 1 (singular, ~1.5e8 as planned)", nfail)

    call build_gs(gs, domega12, d_naive)

    E(1) = 1d-6; E(2) = 0d0; E(3) = 0d0
    dt = 0.05d0
    icomm = 0

    call init_sbe_bloch_solver(sbe, gs, nb, icomm)
    call prepare_qnm(sbe, gs, icomm)
    call adams_moulton_coefs(bj_am)
    tr0 = trace_re_of(sbe)

    diverged = .false.; it_diverged = -1
    do it = 1, 200
      call dt_evolve_bloch_lg(sbe, gs, E, bj_am, dt, icomm)
      tr = trace_re_of(sbe)
      hn = herm_norm_of(sbe)
      if (.not. is_finite(tr) .or. .not. is_finite(hn) .or. &
          hn > 1d-8 .or. abs(tr - tr0) > 1d-6) then
        diverged = .true.
        it_diverged = it
        exit
      end if
    end do

    if (diverged) write(*, '(a,i0,a)') "      (diverged at real-propagator step ", it_diverged, ")"
    call check_true(diverged, &
      "TEETH: gi-style naive dipole (|d|~1.5e8) blows up the SAME real propagator gifix keeps bounded", nfail)
  end subroutine test_gi_naive_diverges

end program test_gifix_propagator
