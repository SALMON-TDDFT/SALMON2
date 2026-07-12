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
!   Positive A (test_gifix_conserves): the {1,2} interband coupling is built
!   from the production xi kernel (degenerate_block_ssbe:xi_block_from_overlap,
!   already covariance-tested by Pc) on a smooth overlap -- an authentic,
!   bounded (O(1)) Hermitian non-Abelian connection. This case proves only
!   that GIVEN an engaged xi the real propagator stays stable (trace conserved
!   to 1e-10, Hermiticity <= 1e-8, every step); it does NOT exercise Tasks
!   2-3's block-detection/xi-assembly code -- it fabricates xi directly.
!
!   Positive B (test_gifix_endtoend): closes that gap. It runs the ACTUAL
!   Tasks 2-3 path end-to-end: REAL prod_dk overlaps -> build_xi(...,
!   fixed_blocks=.true.) (= build_blocks_fixed's k-independent composite
!   partition + gap-isolation fail-closed + identity bypass) -> prepare_matrix's
!   exact blend into gs%d_matrix -> the real prepare_qnm/dt_evolve_bloch_lg
!   propagator. It asserts the identity bypass engaged (xi_ok=.true. on the
!   in-block pair), the isolated fixed block was NOT rejected (n_reject==0),
!   the connection is genuinely nonzero, AND the resulting xi keeps the
!   propagator conserving for 200 steps. THIS is what proves the T2/T3 fixed-
!   block path produces a propagator-stable xi (not just that some xi would be).
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
!   region, so the SAME real propagator MUST blow up. The naive per-k dipole
!   IS the pre-xi-fix fallback (what a per-k 'gi' block whose xi is rejected
!   uses), so its divergence is exactly the failure the gifix path removes --
!   giving the conservation checks above teeth (they are not vacuously
!   satisfied by a propagator that cannot be destabilised).
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
  use gs_info_ssbe,          only: s_sbe_gs_info, sbe_gs_set_replicated_kmap
  use bloch_solver_ssbe,     only: s_sbe_bloch_solver, init_sbe_bloch_solver, &
                                    prepare_qnm, dt_evolve_bloch_lg, adams_moulton_coefs
  use degenerate_block_ssbe, only: xi_block_from_overlap, xi_sign, &
                                    build_xi, blend, same_block, theta_on, theta_off
  use salmon_global,         only: epdir_re1, am_s, sbe_lg_diag, num_kgrid, t_2, &
                                    yn_sbe_gw_collision, sbe_deph_mode
  implicit none

  integer, parameter :: nb = 4, nk = 4
  complex(8), parameter :: zi_ = (0d0, 1d0)
  integer :: nfail

  nfail = 0
  call set_global_defaults()

  call test_gifix_conserves(nfail)    ! direct-rotation xi: propagator stable GIVEN an engaged xi
  call test_gifix_endtoend(nfail)     ! T2/T3 REAL path: build_xi(fixed_blocks=.true.) -> propagator
  call test_gi_naive_diverges(nfail)  ! teeth: naive per-k dipole fallback diverges the same propagator

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

    call sbe_gs_set_replicated_kmap(gs, nk)   ! replicated k layout (kmap = identity)

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

  !======================= nb x nb block frame for the e2e prod_dk =============
  ! Real 2x2 rotation on the {1,2} block, identity on the spectators {3,4}.
  ! A smooth k-sweep of `ang` gives the smooth Bloch frame whose neighbouring-k
  ! overlaps ARE the prod_dk that gifix's build_xi consumes.
  function frame12(ang) result(Q)
    implicit none
    real(8), intent(in) :: ang
    complex(8) :: Q(nb, nb)
    integer :: a
    real(8) :: c, s
    Q = (0d0, 0d0)
    do a = 1, nb
      Q(a, a) = (1d0, 0d0)
    end do
    c = cos(ang); s = sin(ang)
    Q(1, 1) = dcmplx(c, 0d0); Q(1, 2) = dcmplx(-s, 0d0)
    Q(2, 1) = dcmplx(s, 0d0); Q(2, 2) = dcmplx(c, 0d0)
  end function frame12

  !======================= END-TO-END gs via the REAL T2/T3 path ===============
  ! Builds gs entirely through Tasks 2-3's production code:
  !   (a) REAL prod_dk overlaps <u_a(k)|u_b(k+dk)> from a smooth {1,2}-block
  !       rotation (identity spectators), the way test_gisbe_covariance.f90's
  !       build_xi covariance test constructs prod_dk;
  !   (b) build_xi(..., fixed_blocks=.true.) -> build_blocks_fixed (k-independent
  !       composite partition, gap-isolation, fail-closed) + identity bypass ->
  !       the actual gs%xi/gs%xi_ok Tasks 2-3 emit;
  !   (c) prepare_matrix's exact blend (gs_info_ssbe.f90:490-521): d_matrix =
  !       w*xi + (1-w)*i*p/dw, with same_block/xi_ok/blend from the production
  !       module. For this tight fixed block (delta_omega=2d-9 << theta_on) w=1
  !       so d_matrix == the build_xi xi -- i.e. the propagator is driven by the
  !       xi that the real gifix GS path produced, NOT a hand-picked value.
  ! Returns n_reject / in-block xi_ok / in-block |xi| so the caller can assert
  ! the identity bypass engaged and the isolated block was not rejected.
  subroutine build_gs_endtoend(gs, domega12, n_reject_out, xi_ok_inblock, xi12_mag)
    implicit none
    type(s_sbe_gs_info), intent(out) :: gs
    real(8), intent(in) :: domega12
    integer, intent(out) :: n_reject_out
    logical, intent(out) :: xi_ok_inblock
    real(8), intent(out) :: xi12_mag
    real(8), parameter :: two_pi = 6.28318530717958647692d0
    integer, parameter :: nbvec = 3
    integer :: ik, ib, jb, ikx, a
    integer :: bvec(3, nbvec)
    real(8) :: eigen(nb), g, th0, resid, x, w, omega_eps
    complex(8) :: prod_dk(nb, nb, nbvec, nk)
    complex(8) :: Qk(nb, nb), Qn(nb, nb), dpdw(3)

    gs%nk = nk; gs%nb = nb; gs%ne = 4

    call sbe_gs_set_replicated_kmap(gs, nk)   ! replicated k layout (kmap = identity)
    allocate(gs%eigen(nb, nk), gs%occup(nb, nk))
    allocate(gs%delta_omega(nb, nb, nk))
    allocate(gs%d_matrix(nb, nb, 3, nk))
    allocate(gs%p_mod_matrix(nb, nb, 3, nk))
    allocate(gs%xi(nb, nb, 3, nk), gs%xi_ok(nb, nb, nk))

    eigen(1) = 0.30d0 + 0.5d0 * domega12
    eigen(2) = 0.30d0 - 0.5d0 * domega12
    eigen(3) = 1.00d0
    eigen(4) = 1.60d0

    gs%b_matrix = 0d0
    gs%b_matrix(1, 1) = two_pi; gs%b_matrix(2, 2) = two_pi; gs%b_matrix(3, 3) = two_pi
    gs%p_mod_matrix = (0d0, 0d0)   ! spectators inert: blend fallback branch -> 0

    do ik = 1, nk
      gs%eigen(:, ik) = eigen(:)
      gs%occup(1, ik) = 2d0; gs%occup(2, ik) = 0d0   ! ASYMMETRIC (drives the coherence)
      gs%occup(3, ik) = 0d0; gs%occup(4, ik) = 0d0
      do jb = 1, nb
        do ib = 1, nb
          gs%delta_omega(ib, jb, ik) = eigen(ib) - eigen(jb)
        end do
      end do
    end do

    ! (a) REAL prod_dk from a smooth k-dependent rotation of the {1,2} frame;
    !     y/z self-links = identity (mirrors test_gisbe_covariance.f90:447-456).
    bvec(:, 1) = (/ 1, 0, 0 /); bvec(:, 2) = (/ 0, 1, 0 /); bvec(:, 3) = (/ 0, 0, 1 /)
    g = 0.15d0; th0 = 0.20d0
    prod_dk = (0d0, 0d0)
    do ik = 1, nk
      ikx = mod(ik, nk) + 1                          ! periodic +x neighbour
      Qk = frame12(th0 + g * (ik  - 1))
      Qn = frame12(th0 + g * (ikx - 1))
      prod_dk(:, :, 1, ik) = matmul(transpose(conjg(Qk)), Qn)
      do a = 1, nb
        prod_dk(a, a, 2, ik) = (1d0, 0d0)            ! (0,1,0) self-link = I
        prod_dk(a, a, 3, ik) = (1d0, 0d0)            ! (0,0,1) self-link = I
      end do
    end do

    ! (b) THE T2/T3 PRODUCTION PATH: fixed composite block + identity bypass.
    call build_xi(nb, nk, nbvec, bvec, prod_dk, gs%eigen, gs%b_matrix, num_kgrid, &
                  gs%xi, gs%xi_ok, n_reject_out, resid, fixed_blocks=.true.)
    xi_ok_inblock = gs%xi_ok(1, 2, 1)
    xi12_mag      = maxval(abs(gs%xi(1, 2, :, 1)))

    ! (c) production blend (prepare_matrix, gs_info_ssbe.f90:490-521), verbatim.
    omega_eps = 1d-9
    gs%d_matrix = (0d0, 0d0)
    do ik = 1, nk
      do ib = 1, nb
        do jb = 1, nb
          x = abs(gs%delta_omega(ib, jb, ik))
          if (same_block(ib, jb, ik) .and. ib /= jb .and. gs%xi_ok(ib, jb, ik)) then
            w = blend(x, theta_on, theta_off)
            if (x > omega_eps) then
              dpdw(:) = zi_ * gs%p_mod_matrix(ib, jb, :, ik) / gs%delta_omega(ib, jb, ik)
            else
              dpdw(:) = (0d0, 0d0)
            end if
            gs%d_matrix(ib, jb, :, ik) = w * gs%xi(ib, jb, :, ik) + (1d0 - w) * dpdw(:)
          else if (omega_eps < x) then
            gs%d_matrix(ib, jb, :, ik) = zi_ * gs%p_mod_matrix(ib, jb, :, ik) / gs%delta_omega(ib, jb, ik)
          else
            gs%d_matrix(ib, jb, :, ik) = (0d0, 0d0)
          end if
        end do
      end do
    end do
  end subroutine build_gs_endtoend

  !======================= (1b) END-TO-END: T2/T3 path -> propagator ============
  ! Proves the ACTUAL Tasks 2-3 code path (build_xi with fixed_blocks=.true.,
  ! i.e. build_blocks_fixed + identity bypass) emits an xi that keeps the REAL
  ! prepare_qnm/dt_evolve_bloch_lg propagator conserving over 200 steps. This is
  ! the coverage the direct-rotation case (test_gifix_conserves) does NOT give:
  ! that case fabricates xi by calling xi_block_from_overlap directly, so it
  ! only shows "given an engaged xi the propagator is stable"; this case shows
  ! "the fixed-block path actually produces such an engaged, stable xi".
  subroutine test_gifix_endtoend(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    real(8) :: tr0, E(3), bj_am(8, 8), dt, tr, hn, xi12_mag
    integer :: it, icomm, n_reject
    logical :: trace_ok, herm_ok, xi_ok_inblock

    call build_gs_endtoend(gs, 2d-9, n_reject, xi_ok_inblock, xi12_mag)

    ! T2/T3 path health: identity bypass engaged xi for the in-block pair, and
    ! the gap-isolated fixed block survived build_xi's fail-closed check.
    call check_true(xi_ok_inblock, &
      "gifix e2e: build_xi(fixed_blocks=.true.) engaged xi_ok for in-block pair {1,2} (identity bypass)", nfail)
    call check_true(n_reject == 0, &
      "gifix e2e: build_xi n_reject==0 for the isolated fixed block", nfail)
    call check_true(xi12_mag > 1d-6, &
      "gifix e2e: build_xi produced a NONZERO in-block connection (xi actually engaged)", nfail)

    E(1) = 1d-6; E(2) = 0d0; E(3) = 0d0
    dt = 0.05d0
    icomm = 0

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

    call check_true(trace_ok, &
      "gifix e2e: trace conserved to 1e-10 (build_xi fixed_blocks xi -> real propagator, 200 steps)", nfail)
    call check_true(herm_ok, &
      "gifix e2e: Hermiticity <= 1e-8 (build_xi fixed_blocks xi -> real propagator, 200 steps)", nfail)
  end subroutine test_gifix_endtoend

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
