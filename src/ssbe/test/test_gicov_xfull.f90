! src/ssbe/test/test_gicov_xfull.f90
! LG-SBE gicov X-full (Task 1 of plans/2026-07-03-gicov-xfull.md): the two
! Task-1 validation tests for the RHS rewrite (drop the analytic interband
! dipole term from gicov_rhs because the full-band covariant transport already
! supplies it).
!
! WHY THIS FILE EXISTS
!   X-full replaces the fixed-block gicov (X-closed) with a SINGLE full-band
!   Wilson-line block (block_id === 1 for ALL bands): build_block_transport
!   polar-factors the WHOLE nb x nb overlap M, so covariant_grad_block's
!   D_cov rho = d_k rho - i[xi_full,rho], evaluated on the FULL density,
!   already contains the ENTIRE interband coupling (xi_full's off-diagonal
!   entries ARE the physical dipole), making the separate analytic
!   -i*sum_a E_a[d_matrix_a,rho] term in gicov_rhs pure double-counting. Two
!   tests validate this:
!
!   Test D (drop-dipole equivalence, pure kernel math): proves the underlying
!     CLAIM directly -- for a fixture with EXACTLY ZERO diagonal (intraband)
!     Berry connection (so xi_full == xi_inter, no separate intraband piece
!     to confuse the comparison), covariant_grad_block(full transport) equals
!     covariant_grad_block(bare, U=I) minus i[xi_inter,rho] to a
!     stencil-matched tolerance. Zero diagonal connection is obtained with a
!     REAL 2x2 rotation V(k) = exp(-i theta(k) sigma_y): sigma_y is pure
!     imaginary, so -i*theta*sigma_y is REAL antisymmetric and V(k) is a real
!     orthogonal matrix; its real columns u_n(k) satisfy <u_n|d_k u_n> = 0
!     exactly (real AND purely-imaginary => zero) -- unlike a generic complex
!     unitary transport, which can carry a nonzero intraband Berry phase that
!     would contaminate a naive "xi_full == xi_inter" comparison. This is
!     PURE KERNEL MATH: only degenerate_block_ssbe (covariant_grad_block +
!     build_block_transport), no gs_info_ssbe / bloch_solver_ssbe /
!     salmon_global -- standalone-compilable exactly like test_covariant_grad.f90
!     if extracted with its own zi_/two_pi/check_true/nbr/dotc helpers:
!       gfortran -O2 degenerate_block_ssbe.f90 test/test_gicov_xfull.f90 \
!                -o /tmp/t_xfull -framework Accelerate && /tmp/t_xfull
!     (the combined build below runs the WHOLE file, incl. Test U3-full,
!     against the SAME degenerate_block_ssbe.f90 source compiled into
!     build_local -- byte-identical kernel either way.)
!
!   Test U3-full (production path, full-band, WITH a band crossing): the
!     mirror-lesson gate (test_gicov_propagator.f90's own rationale) -- a
!     correct RHS rewrite is necessary but not sufficient; this drives the REAL
!     production propagator (prepare_qnm + dt_evolve_bloch_lg, sbe_lg_degen=
!     'gicov') for >=200 RK4 steps with block_id === 1 (ALL nb bands ONE
!     block, window option (A)) and a deliberate near-pi/2 "crossing" rotation
!     on one prod_dk link (bands 2,3) -- exactly the coarse-mesh band-crossing
!     scenario that broke the OLD fixed-block X-closed transport -- and asserts
!     boundedness + gauge invariance survive it.
!
!     CODEX REVIEW POINT (must, not optional): calc_current_bloch_lg's gicov
!     branch contracts the physical rho with gs%p_tm_matrix AND (when
!     sbe%flag_vnl_correction) gs%rvnl_tm_matrix (bloch_solver_ssbe.f90
!     :1109-1113). An "arbitrary per-k gauge W(k)" test MUST transform BOTH
!     of those current operators (p_tm_matrix, rvnl_tm_matrix) -> W(k)^H (.) W(k)
!     the SAME way it transforms prod_dk/u_transport/rho, or the gauged run's
!     current J(t) would spuriously differ from the ungauged run's even for a
!     CORRECT X-full implementation (a false test failure). This test enables
!     sbe%flag_vnl_correction = .true. and gauges rvnl_tm_matrix explicitly so
!     that omitting its transform would be caught.
!
! BUILD (combined; links build_local objects -- model on test_gicov_propagator.f90):
!   find <repo>/build_local/src/CMakeFiles/salmon.dir -name '*.o' ! -name 'main.f90.o' > objs.txt
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     -I<repo>/build_local -J<clean_dir> \
!     -c <repo>/src/ssbe/test/test_gicov_xfull.f90 -o <clean_dir>/test_gicov_xfull.o
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     @objs.txt <clean_dir>/test_gicov_xfull.o -o <clean_dir>/test_gicov_xfull \
!     -framework Accelerate -lm -ldl
!   <clean_dir>/test_gicov_xfull
!
! NOTE ON STRUCTURE: gfortran (even -std=f2008/f2018) rejects an internal
! procedure that itself CONTAINS further internal procedures (no recursive
! nesting), so every helper below is a flat, PROGRAM-level internal procedure
! (single level of CONTAINS) rather than nested inside test_D_.../test_U3full;
! Test D's (nb=2,nk=16) and Test U3-full's (nb=4,nk=8) own PARAMETER
! declarations are local to their own subroutines and shadow the identically
! named ones nowhere else declared at program scope (no conflict).
!
program test_gicov_xfull
  implicit none
  complex(8), parameter :: zi_ = (0d0, 1d0)
  real(8),    parameter :: two_pi = 6.28318530717958647692d0
  ! Test U3-full fixture size (shared by all its flat helper procedures below)
  integer, parameter :: nb = 4, nk = 8, nsteps = 200
  integer, parameter :: ik0 = 4     ! the single "band crossing" link (k=4 -> k=5)
  integer :: nfail

  nfail = 0
  call test_D_drop_dipole_equivalence(nfail)
  call test_U3full(nfail)

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_gicov_xfull checks passed."
  end if

contains

  !======================= assert / finite helpers (shared) ===================
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

  logical function is_finite(x) result(ok)
    implicit none
    real(8), intent(in) :: x
    ok = (x == x) .and. (abs(x) <= huge(1d0))     ! rejects NaN (x/=x) and +-Inf
  end function is_finite

  ! Periodic k-neighbour on a uniform grid num_kgrid, shift s along axis (same
  ! index ordering as build_ik_neighbor: i1 fastest). Generic; used by Test D.
  integer function nbr(ik, axis, s, num_kgrid) result(j)
    implicit none
    integer, intent(in) :: ik, axis, s, num_kgrid(3)
    integer :: i1, i2, i3, j1, j2, j3, n1, n2, n3
    n1 = num_kgrid(1); n2 = num_kgrid(2); n3 = num_kgrid(3)
    i1 = mod(ik - 1, n1); i2 = mod((ik - 1) / n1, n2); i3 = (ik - 1) / (n1 * n2)
    j1 = i1; j2 = i2; j3 = i3
    if (axis == 1) j1 = modulo(i1 + s, n1)
    if (axis == 2) j2 = modulo(i2 + s, n2)
    if (axis == 3) j3 = modulo(i3 + s, n3)
    j = j3 * n1 * n2 + j2 * n1 + j1 + 1
  end function nbr

  complex(8) function dotc(a, b) result(d)
    implicit none
    complex(8), intent(in) :: a(:), b(:)
    d = sum(conjg(a) * b)
  end function dotc

  !=============================================================================
  ! Test D -- drop-dipole equivalence (pure kernel math; see file header).
  !=============================================================================
  subroutine test_D_drop_dipole_equivalence(nfail)
    use degenerate_block_ssbe, only: covariant_grad_block, build_block_transport
    implicit none
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 2, nk = 16, nbvec = 3
    integer, parameter :: m_max = 4     ! (nk-1)/2=7 -> capped to 4 (full, uncapped regime)
    real(8), parameter :: cw(4) = (/ 4d0/5d0, -1d0/5d0, 4d0/105d0, -1d0/280d0 /)
    real(8), parameter :: theta_amp = 1.5d-4    ! small-rotation regime (Test-C-like scale, margin over the stencil-matched tol)
    integer :: num_kgrid(3), bvec(3, nbvec), block_id(nb, nk), n_reject
    complex(8) :: Vk(nb, nb, nk)
    complex(8) :: prod_dk(nb, nb, nbvec, nk)
    complex(8) :: U_full(nb, nb, 3, nk), U_bare(nb, nb, 3, nk)
    complex(8) :: rho(nb, nb, nk)
    complex(8) :: Dq_full(nb, nb, 3, nk), Dq_bare(nb, nb, 3, nk)
    complex(8) :: xi_full(nb, nb), xi_inter, diag1, diag2, com(nb, nb), tgt(nb, nb)
    complex(8) :: u1k(nb), u2k(nb), u1kp(nb), u1km(nb), u2kp(nb), u2km(nb)
    real(8) :: dk(3), theta(nk), th, resid, dqfull_scale, dqbare_scale, xi_mag, diag_mag
    integer :: ik, axis, n, kp, km, m

    num_kgrid = (/ nk, 1, 1 /)
    bvec(:, 1) = (/ 1, 0, 0 /); bvec(:, 2) = (/ 0, 1, 0 /); bvec(:, 3) = (/ 0, 0, 1 /)
    dk = (/ 0.4d0 / dble(nk), 1d0, 1d0 /)
    block_id(:, :) = 1     ! X-full: ONE full-band block (here the whole nb=2 system)

    ! theta(k): smooth, small-amplitude rotation angle (real sigma_y generator)
    do ik = 1, nk
      th = two_pi * dble(ik - 1) / dble(nk)
      theta(ik) = theta_amp * sin(th + 0.4d0)
    end do

    ! V(k) = exp(-i theta(k) sigma_y) = [[cos,-sin],[sin,cos]] -- REAL orthogonal
    ! (sigma_y is pure imaginary, so -i*theta*sigma_y is real antisymmetric).
    do ik = 1, nk
      Vk(1, 1, ik) = dcmplx(cos(theta(ik)), 0d0)
      Vk(1, 2, ik) = dcmplx(-sin(theta(ik)), 0d0)
      Vk(2, 1, ik) = dcmplx(sin(theta(ik)), 0d0)
      Vk(2, 2, ik) = dcmplx(cos(theta(ik)), 0d0)
    end do

    ! prod_dk(:,:,axis=1,ik) = V(k)^H V(k+e); axes 2,3 identity (num_kgrid=1
    ! there, but build_block_transport still probes every axis in bvec).
    prod_dk = (0d0, 0d0)
    do ik = 1, nk
      kp = nbr(ik, 1, 1, num_kgrid)
      prod_dk(:, :, 1, ik) = matmul(transpose(conjg(Vk(:, :, ik))), Vk(:, :, kp))
      do n = 1, nb
        prod_dk(n, n, 2, ik) = (1d0, 0d0)
        prod_dk(n, n, 3, ik) = (1d0, 0d0)
      end do
    end do

    ! U_full = build_block_transport(block_id===1) -- the actual X-full
    ! production path (polar factor of the FULL 2x2 overlap M). prod_dk is
    ! already exactly unitary here, so the polar factor reproduces it exactly.
    call build_block_transport(nb, nk, nbvec, bvec, prod_dk, block_id, num_kgrid, &
                                U_full, n_reject)
    call check_true(n_reject == 0, &
      "D: build_block_transport(block_id===1) accepts the full 2x2 overlap (no reject)", nfail)

    ! U_bare = identity everywhere (bare-stencil reference, matches Test A of
    ! test_covariant_grad.f90).
    U_bare = (0d0, 0d0)
    do ik = 1, nk
      do axis = 1, 3
        do n = 1, nb
          U_bare(n, n, axis, ik) = (1d0, 0d0)
        end do
      end do
    end do

    ! smooth Hermitian test density with BOTH diagonal and off-diagonal parts,
    ! so [xi_inter,rho] is genuinely nonzero.
    do ik = 1, nk
      th = two_pi * dble(ik - 1) / dble(nk)
      rho(1, 1, ik) = dcmplx(0.70d0 + 0.15d0 * cos(th), 0d0)
      rho(2, 2, ik) = dcmplx(0.55d0 + 0.10d0 * sin(th + 0.6d0), 0d0)
      ! NOTE: the phase factor's th-coefficient MUST be an integer (here 1) so
      ! rho is genuinely 2*pi-periodic in the k-index; a non-integer coefficient
      ! (e.g. 0.3*th) creates a real discontinuity at the periodic wraparound
      ! that the finite-difference stencil amplifies as O(1/dk) -- discovered
      ! empirically via a diverges-with-nk regression check during Task 1 TDD.
      rho(1, 2, ik) = (0.20d0 + 0.05d0 * cos(th + 1.1d0)) * exp(zi_ * th)
      rho(2, 1, ik) = conjg(rho(1, 2, ik))
    end do

    call covariant_grad_block(nb, nk, nbvec, bvec, num_kgrid, U_full, rho, dk, Dq_full, 1, nk, (/.true.,.true.,.true./))
    call covariant_grad_block(nb, nk, nbvec, bvec, num_kgrid, U_bare, rho, dk, Dq_bare, 1, nk, (/.true.,.true.,.true./))

    ! xi_inter(k) = i <u1(k)| stencil-D_k u2(k)>, using the SAME cw/m_max shells
    ! covariant_grad_block itself uses (not an analytic-continuum derivative) --
    ! this is what makes the comparison below "stencil-matched". diag1/diag2
    ! are the analogous intraband terms <u_n|stencil-D_k u_n>, asserted ~0.
    resid = 0d0; xi_mag = 0d0; diag_mag = 0d0
    do ik = 1, nk
      u1k = Vk(:, 1, ik); u2k = Vk(:, 2, ik)
      xi_inter = (0d0, 0d0); diag1 = (0d0, 0d0); diag2 = (0d0, 0d0)
      do m = 1, m_max
        kp = nbr(ik, 1,  m, num_kgrid)
        km = nbr(ik, 1, -m, num_kgrid)
        u1kp = Vk(:, 1, kp); u1km = Vk(:, 1, km)
        u2kp = Vk(:, 2, kp); u2km = Vk(:, 2, km)
        xi_inter = xi_inter + cw(m) * (dotc(u1k, u2kp) - dotc(u1k, u2km)) / dk(1)
        diag1    = diag1    + cw(m) * (dotc(u1k, u1kp) - dotc(u1k, u1km)) / dk(1)
        diag2    = diag2    + cw(m) * (dotc(u2k, u2kp) - dotc(u2k, u2km)) / dk(1)
      end do
      xi_inter = zi_ * xi_inter
      diag1    = zi_ * diag1
      diag2    = zi_ * diag2
      diag_mag = max(diag_mag, abs(diag1), abs(diag2))
      xi_mag   = max(xi_mag, abs(xi_inter))

      xi_full = (0d0, 0d0)
      xi_full(1, 2) = xi_inter
      xi_full(2, 1) = conjg(xi_inter)

      com = matmul(xi_full, rho(:, :, ik)) - matmul(rho(:, :, ik), xi_full)
      do axis = 1, 3
        tgt = Dq_bare(:, :, axis, ik)
        if (axis == 1) tgt = tgt - zi_ * com
        resid = max(resid, maxval(abs(Dq_full(:, :, axis, ik) - tgt)))
      end do
    end do

    dqfull_scale = maxval(abs(Dq_full(:, :, 1, :)))
    dqbare_scale = maxval(abs(Dq_bare(:, :, 1, :)))

    call check_true(diag_mag < 1d-6, &
      "D: diagonal (intraband) Berry connection ~0 for the REAL sigma_y rotation " // &
      "(confirms xi_full = xi_inter, no separate intraband piece)", nfail)
    call check_true(xi_mag > 1d-4, &
      "D: interband connection xi_inter is genuinely nonzero (not a vacuous check)", nfail)
    call check_true(dqfull_scale > 1d-3 .and. dqbare_scale > 1d-3, &
      "D: Dq_full and Dq_bare are both nonzero (nontrivial fixture)", nfail)
    call check_true(resid < 5d-6, &
      "D: covariant_grad_block(full) = covariant_grad_block(bare) - i[xi_inter,rho] " // &
      "to stencil-matched tol (drop-dipole equivalence -- the rewrite's core claim)", nfail)

    write(*, '(a,es10.2,a,es10.2,a,es10.2,a,es10.2,a,es10.2)') &
      "      D residuals  equivalence=", resid, "  diag-Berry=", diag_mag, &
      "  xi_inter=", xi_mag, "  Dq_full=", dqfull_scale, "  Dq_bare=", dqbare_scale
  end subroutine test_D_drop_dipole_equivalence

  !=============================================================================
  ! Test U3-full -- full-band (block_id===1), WITH a band crossing, arbitrary
  ! per-k gauge (incl. the current operators), >=200 RK4 steps. See file header.
  !=============================================================================
  subroutine test_U3full(nfail)
    use gs_info_ssbe, only: s_sbe_gs_info
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs, gsg
    complex(8) :: rho(nb, nb, nk), rhog(nb, nb, nk), W(nb, nb), Wh(nb, nb)
    real(8) :: E(3), dt, Nocc
    real(8) :: tr_u(0:nsteps), tr_g(0:nsteps), Ju(3, 0:nsteps), Jg(3, 0:nsteps)
    real(8) :: hn_u, hn_g, Jmax_u, Jmax_g
    real(8) :: rho_max_u, rho_max_g, coh_drift_u, coh_drift_g
    real(8) :: tr_drift, gauge_tr, gauge_J, rho_max_all, jmax_all
    logical :: ok_u, ok_g
    integer :: it, ik

    call set_globals_u3()

    ! ---- U3-full gate: ungauged vs gauged (incl. current operators) ----------
    E(1) = 0.10d0; E(2) = 0d0; E(3) = 0d0
    dt = 0.05d0

    call build_gs_u3(gs)
    call build_rho_init_u3(rho)
    call run_gicov_u3(gs, rho, E, dt, tr_u, Ju, hn_u, ok_u, Jmax_u, rho_max_u, coh_drift_u)
    Nocc = tr_u(0)

    call build_gs_u3(gsg)
    call gauge_gs_inplace_u3(gsg)
    do ik = 1, nk
      W  = make_W_u3(ik)
      Wh = hconj_u3(W)
      rhog(:, :, ik) = matmul(matmul(Wh, rho(:, :, ik)), W)
    end do
    call run_gicov_u3(gsg, rhog, E, dt, tr_g, Jg, hn_g, ok_g, Jmax_g, rho_max_g, coh_drift_g)

    ! (1) trace conserved every step (both runs)
    tr_drift = 0d0
    do it = 0, nsteps
      tr_drift = max(tr_drift, abs(tr_u(it) - Nocc), abs(tr_g(it) - tr_g(0)))
    end do
    call check_true(is_finite(tr_drift) .and. tr_drift < 1d-8, &
      "U3f.1: |Tr rho(t) - Nocc| < 1e-8 over 200 gicov RK4 steps (full-band, w/ crossing)", nfail)

    ! (2) Hermiticity bounded every step
    call check_true(is_finite(hn_u) .and. is_finite(hn_g) .and. hn_u <= 1d-8 .and. hn_g <= 1d-8, &
      "U3f.2: herm_norm(rho(t)) <= 1e-8 over 200 steps (no blow-up)", nfail)

    ! (3) no NaN/Inf anywhere
    call check_true(ok_u .and. ok_g, &
      "U3f.3: all trace/current/herm finite (no NaN) over 200 steps", nfail)

    ! (4) gauge invariance of Tr and current J(t) -- INCLUDING the p_tm_matrix/
    !     rvnl_tm_matrix transform (codex review point; see file header).
    gauge_tr = 0d0; gauge_J = 0d0
    do it = 0, nsteps
      gauge_tr = max(gauge_tr, abs(tr_g(it) - tr_u(it)))
      gauge_J  = max(gauge_J, abs(Jg(1, it) - Ju(1, it)), &
                              abs(Jg(2, it) - Ju(2, it)), abs(Jg(3, it) - Ju(3, it)))
    end do
    call check_true(is_finite(gauge_tr) .and. gauge_tr < 1d-8, &
      "U3f.4a: Tr(gauged) == Tr(ungauged) to 1e-8 every step", nfail)
    call check_true(is_finite(gauge_J) .and. gauge_J < 1d-8, &
      "U3f.4b: J(gauged) == J(ungauged) to 1e-8 every step (incl. gauged current operators)", nfail)

    ! (5) non-vacuous
    call check_true(Jmax_u > 1d-6, &
      "U3f.5: max|J(t)| >> 0 (dynamics non-trivial, boundedness not vacuous)", nfail)

    ! (6)/(7) direct |rho|/|J| ceilings (hardening; catches moderate blow-up
    ! that 1/2/4's roundoff-scale thresholds would miss)
    rho_max_all = max(rho_max_u, rho_max_g)
    call check_true(is_finite(rho_max_all) .and. rho_max_all < 10d0, &
      "U3f.6: |rho| bounded (max_t maxval(abs(rho(t))) < 10)", nfail)
    jmax_all = max(Jmax_u, Jmax_g)
    call check_true(is_finite(jmax_all) .and. jmax_all < 100d0, &
      "U3f.7: |J| bounded above (max_t |J(t)| < 100)", nfail)

    ! (8) some off-diagonal coherence genuinely driven (covariant term engaged,
    ! not just the diagonal energy term)
    call check_true(is_finite(coh_drift_u) .and. coh_drift_u > 1d-6, &
      "U3f.8: max_{t,k,i/=j} |rho_ij(t) - rho_ij(0)| > 1e-6 (coherences genuinely driven)", nfail)

    write(*, '(a,es10.2,a,es10.2,a,es10.2)') &
      "      U3f trace-drift=", tr_drift, "  herm(u/g)=", max(hn_u, hn_g), "  |J|max=", Jmax_u
    write(*, '(a,es10.2,a,es10.2)') &
      "      U3f gauge-resid  Tr=", gauge_tr, "  J=", gauge_J

    ! ---- TEETH: a non-covariant (broken-telescoping) transport corruption ----
    call run_teeth_corrupted_u3(nfail)
  end subroutine test_U3full

  !======================= salmon_global fixture (U3-full) =====================
  subroutine set_globals_u3()
    use salmon_global, only: epdir_re1, am_s, sbe_lg_diag, num_kgrid, t_2, &
                              yn_sbe_gw_collision, sbe_deph_mode, sbe_lg_degen
    implicit none
    epdir_re1(1) = 1d0; epdir_re1(2) = 0d0; epdir_re1(3) = 0d0
    am_s = 4
    num_kgrid(1) = nk; num_kgrid(2) = 1; num_kgrid(3) = 1
    t_2 = 1d30
    sbe_lg_diag = 0
    sbe_lg_degen = 'gicov'
    yn_sbe_gw_collision = 'n'
    sbe_deph_mode = ''
  end subroutine set_globals_u3

  !======================= trace/herm probes (over sbe%qnm_new) ===============
  real(8) function trace_re_of_u3(sbe) result(tr)
    use bloch_solver_ssbe, only: s_sbe_bloch_solver
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    integer :: ik, ib
    tr = 0d0
    do ik = sbe%ik_min, sbe%ik_max
      do ib = 1, sbe%nb
        tr = tr + dble(sbe%qnm_new(ib, ib, ik))
      end do
    end do
  end function trace_re_of_u3

  real(8) function herm_norm_of_u3(sbe) result(hn)
    use bloch_solver_ssbe, only: s_sbe_bloch_solver
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
  end function herm_norm_of_u3

  ! physical rho reconstructed AFTER a completed step directly from
  ! sbe%qnm_new (mirrors calc_current_bloch_lg's own gicov reconstruction);
  ! exp_iphi is time-constant (prepare_qnm sets it once), so this is exact.
  subroutine rho_from_qnm_new_u3(sbe, rho_out)
    use bloch_solver_ssbe, only: s_sbe_bloch_solver
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    complex(8), intent(out) :: rho_out(nb, nb, nk)
    integer :: ik, ib, jb
    rho_out = (0d0, 0d0)
    do ik = sbe%ik_min, sbe%ik_max
      do jb = 1, nb
        do ib = 1, nb
          if (ib == jb) then
            rho_out(ib, jb, ik) = sbe%qnm_new(ib, ib, ik)
          else
            rho_out(ib, jb, ik) = sbe%exp_iphi(ib, jb, ik) * sbe%qnm_new(ib, jb, ik)
          end if
        end do
      end do
    end do
  end subroutine rho_from_qnm_new_u3

  !======================= small matrix / gauge helpers (U3-full) =============
  function hconj_u3(A) result(Ah)
    implicit none
    complex(8), intent(in) :: A(nb, nb)
    complex(8) :: Ah(nb, nb)
    Ah = transpose(conjg(A))
  end function hconj_u3

  integer function knext_u3(ik, axis) result(kn)
    implicit none
    integer, intent(in) :: ik, axis
    if (axis == 1) then
      kn = mod(ik, nk) + 1
    else
      kn = ik            ! num_kgrid=1 on axes 2,3: neighbour is self
    end if
  end function knext_u3

  ! Per-k gauge W(k), restricted to mix ONLY the EXACTLY-DEGENERATE subspace
  ! {2,3} (independent phases on singletons 1,4). NOTE (found empirically
  ! during Task 1 TDD): gicov_rhs's energy term is "-i*delta_omega(ib,jb)*rho"
  ! -- a DIAGONAL-H0-in-whatever-basis-rho-is-expressed assumption, not a full
  ! matrix commutator -- so it is gauge-covariant ONLY under transforms that do
  ! NOT mix bands with DIFFERENT delta_omega (pure per-band phases, or mixing
  ! strictly within a degenerate subspace where delta_omega=0). A W(k) mixing
  ! bands with UNEQUAL eigenvalues breaks gauge invariance genuinely (an O(1)
  ! effect, not roundoff) -- this is a pre-existing structural property of
  ! gicov_rhs, unrelated to the Task 1 dipole-drop, and matches why
  ! test_gicov_rhs.f90/test_gicov_propagator.f90 always used an EXACTLY
  ! degenerate {2,3} block for their own gauge tests. build_gs_u3 below sets
  ! eigen(2)=eigen(3) for exactly this reason. The "crossing" transport
  ! (prod_dk) is a SEPARATE, unrelated structure and stays on bands {2,3}
  ! regardless (a transport/overlap property, not a gauge-freedom one).
  function make_W_u3(ik) result(W)
    implicit none
    integer, intent(in) :: ik
    complex(8) :: W(nb, nb)
    real(8) :: t, phi, chi, xi, psi, a1, a4, c, s
    t   = dble(ik)
    phi = 0.37d0 + 0.13d0 * t
    chi = 0.21d0 + 0.29d0 * t
    xi  = 0.53d0 - 0.17d0 * t
    psi = 0.11d0 + 0.19d0 * t
    a1  = 0.41d0 + 0.23d0 * t
    a4  = 0.61d0 - 0.31d0 * t
    c = cos(phi); s = sin(phi)
    W = (0d0, 0d0)
    W(1, 1) = exp(zi_ * a1)
    W(2, 2) = exp(zi_ * psi) * c * exp( zi_ * chi)
    W(2, 3) = exp(zi_ * psi) * s * exp( zi_ * xi )
    W(3, 2) = exp(zi_ * psi) * (-s) * exp(-zi_ * xi )
    W(3, 3) = exp(zi_ * psi) * c * exp(-zi_ * chi)
    W(4, 4) = exp(zi_ * a4)
  end function make_W_u3

  !======================= synthetic X-full gs fixture =========================
  ! nb=4, block_id === 1 (ALL bands ONE full-band block, window option A).
  ! prod_dk axis=1: per-band phase transport at every link EXCEPT the single
  ! "crossing" link ik0, where bands {2,3} undergo a near-pi/2 rotation (the
  ! coarse-mesh band-crossing scenario X-closed could not survive).
  ! u_transport is built by build_block_transport(block_id===1) from prod_dk
  ! (the actual production path, not hand-built). d_matrix = 0 (X-full: the
  ! covariant transport supplies the whole field term). p_tm_matrix and
  ! rvnl_tm_matrix are two INDEPENDENT Hermitian operators (different formulas)
  ! so a missed gauge-transform of either would show up distinctly in J(t).
  subroutine build_gs_u3(gs)
    use gs_info_ssbe, only: s_sbe_gs_info
    use degenerate_block_ssbe, only: build_block_transport
    use salmon_global, only: num_kgrid
    implicit none
    type(s_sbe_gs_info), intent(out) :: gs
    real(8) :: eigen(nb), dphi(nb), t, thetac, phia, phib
    integer :: ik, ib, jb, n_reject_local

    gs%nk = nk; gs%nb = nb; gs%ne = 6
    allocate(gs%eigen(nb, nk), gs%occup(nb, nk), gs%kweight(nk))
    allocate(gs%delta_omega(nb, nb, nk))
    allocate(gs%p_tm_matrix(nb, nb, 3, nk))
    allocate(gs%rvnl_tm_matrix(nb, nb, 3, nk))
    allocate(gs%d_matrix(nb, nb, 3, nk))
    allocate(gs%u_transport(nb, nb, 3, nk))
    allocate(gs%block_id(nb, nk))
    allocate(gs%bvec(3, 3))
    allocate(gs%prod_dk(nb, nb, 3, nk))

    ! bands {2,3} EXACTLY degenerate (0.90): the "crossing" transport link AND
    ! the make_W_u3 gauge both live on this pair -- see make_W_u3's comment for
    ! why (gicov_rhs's energy term needs delta_omega(2,3)=0 there to stay
    ! gauge-covariant under a mixing gauge). Singletons 1,4 have distinct
    ! energies (only ever given independent PHASES by make_W_u3, never mixed).
    eigen = (/ 0.30d0, 0.90d0, 0.90d0, 1.50d0 /)
    dphi  = (/ 0.11d0, 0.17d0, -0.13d0, 0.07d0 /)

    gs%volume = 1d0
    gs%b_matrix = 0d0
    gs%b_matrix(1, 1) = two_pi; gs%b_matrix(2, 2) = two_pi; gs%b_matrix(3, 3) = two_pi

    gs%nbvec = 3
    gs%bvec(:, 1) = (/ 1, 0, 0 /); gs%bvec(:, 2) = (/ 0, 1, 0 /); gs%bvec(:, 3) = (/ 0, 0, 1 /)

    gs%d_matrix = (0d0, 0d0)      ! X-full: no analytic dipole
    gs%prod_dk  = (0d0, 0d0)
    gs%block_id = 1               ! X-full: ONE full-band block covering ALL bands

    do ik = 1, nk
      t = dble(ik - 1)
      gs%kweight(ik) = 1d0
      gs%eigen(:, ik) = eigen(:)
      gs%occup(:, ik) = 1.5d0     ! equal, full occupation across ALL bands

      do jb = 1, nb
        do ib = 1, nb
          gs%delta_omega(ib, jb, ik) = eigen(ib) - eigen(jb)
        end do
      end do

      ! two INDEPENDENT Hermitian velocity-like operators, axis 1 only
      do ib = 1, nb
        gs%p_tm_matrix(ib, ib, 1, ik)    = dcmplx(0.50d0 + 0.10d0 * dble(ib), 0d0)
        gs%rvnl_tm_matrix(ib, ib, 1, ik) = dcmplx(0.20d0 + 0.05d0 * dble(ib), 0d0)
        do jb = ib + 1, nb
          gs%p_tm_matrix(ib, jb, 1, ik) = (0.40d0 + 0.10d0 * dble(ib + jb)) * &
                                          exp(zi_ * (0.3d0 * ib + 0.5d0 * jb + 0.2d0 * t))
          gs%p_tm_matrix(jb, ib, 1, ik) = conjg(gs%p_tm_matrix(ib, jb, 1, ik))
          gs%rvnl_tm_matrix(ib, jb, 1, ik) = (0.25d0 + 0.07d0 * dble(ib + jb)) * &
                                          exp(zi_ * (0.6d0 * ib - 0.4d0 * jb + 0.35d0 * t + 0.9d0))
          gs%rvnl_tm_matrix(jb, ib, 1, ik) = conjg(gs%rvnl_tm_matrix(ib, jb, 1, ik))
        end do
      end do
      gs%p_tm_matrix(:, :, 2, ik) = (0d0, 0d0); gs%p_tm_matrix(:, :, 3, ik) = (0d0, 0d0)
      gs%rvnl_tm_matrix(:, :, 2, ik) = (0d0, 0d0); gs%rvnl_tm_matrix(:, :, 3, ik) = (0d0, 0d0)

      ! prod_dk axis=1: per-band phase link, EXCEPT the crossing link ik0
      if (ik == ik0) then
        thetac = 1.45d0     ! near pi/2 (~83 deg) -- the "band crossing"
        phia = 0.20d0; phib = -0.30d0
        gs%prod_dk(1, 1, 1, ik) = exp(zi_ * dphi(1))
        gs%prod_dk(2, 2, 1, ik) =  cos(thetac) * exp( zi_ * phia)
        gs%prod_dk(2, 3, 1, ik) = -sin(thetac) * exp( zi_ * phib)
        gs%prod_dk(3, 2, 1, ik) =  sin(thetac) * exp(-zi_ * phib)
        gs%prod_dk(3, 3, 1, ik) =  cos(thetac) * exp(-zi_ * phia)
        gs%prod_dk(4, 4, 1, ik) = exp(zi_ * dphi(4))
      else
        gs%prod_dk(1, 1, 1, ik) = exp(zi_ * dphi(1))
        gs%prod_dk(2, 2, 1, ik) = exp(zi_ * dphi(2))
        gs%prod_dk(3, 3, 1, ik) = exp(zi_ * dphi(3))
        gs%prod_dk(4, 4, 1, ik) = exp(zi_ * dphi(4))
      end if
      do ib = 1, nb
        gs%prod_dk(ib, ib, 2, ik) = (1d0, 0d0)
        gs%prod_dk(ib, ib, 3, ik) = (1d0, 0d0)
      end do
    end do

    call build_block_transport(nb, nk, gs%nbvec, gs%bvec, gs%prod_dk, gs%block_id, &
                                num_kgrid, gs%u_transport, n_reject_local)
  end subroutine build_gs_u3

  ! Apply the arbitrary per-k gauge W(k) to EVERY gauge-covariant quantity:
  ! prod_dk (-> u_transport rebuilt via the production build_block_transport
  ! path) AND the current operators p_tm_matrix/rvnl_tm_matrix (codex review
  ! point -- omitting these would make J(t) spuriously gauge-dependent).
  subroutine gauge_gs_inplace_u3(gs)
    use gs_info_ssbe, only: s_sbe_gs_info
    use degenerate_block_ssbe, only: build_block_transport
    use salmon_global, only: num_kgrid
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    complex(8) :: Wi(nb, nb), Whi(nb, nb), Wkp(nb, nb)
    integer :: ik, kp, n_reject_local

    do ik = 1, nk
      Wi  = make_W_u3(ik)
      Whi = hconj_u3(Wi)
      kp  = knext_u3(ik, 1)
      Wkp = make_W_u3(kp)
      gs%prod_dk(:, :, 1, ik) = matmul(matmul(Whi, gs%prod_dk(:, :, 1, ik)), Wkp)
      gs%prod_dk(:, :, 2, ik) = matmul(matmul(Whi, gs%prod_dk(:, :, 2, ik)), Wi)
      gs%prod_dk(:, :, 3, ik) = matmul(matmul(Whi, gs%prod_dk(:, :, 3, ik)), Wi)
    end do
    call build_block_transport(nb, nk, gs%nbvec, gs%bvec, gs%prod_dk, gs%block_id, &
                                num_kgrid, gs%u_transport, n_reject_local)

    do ik = 1, nk
      Wi  = make_W_u3(ik)
      Whi = hconj_u3(Wi)
      gs%p_tm_matrix(:, :, 1, ik)    = matmul(matmul(Whi, gs%p_tm_matrix(:, :, 1, ik)),    Wi)
      gs%rvnl_tm_matrix(:, :, 1, ik) = matmul(matmul(Whi, gs%rvnl_tm_matrix(:, :, 1, ik)), Wi)
    end do
    ! d_matrix is 0 -- nothing to gauge there.
  end subroutine gauge_gs_inplace_u3

  ! Deliberately NON-covariant per-k corruption of u_transport: a per-k,
  ! per-band phase applied ONE-SIDED (not paired with the compensating W(k+e)
  ! on the other end of each link, unlike a genuine gauge transform), so the
  ! multi-shell Wilson-line composition no longer telescopes to anything
  ! physical.
  !
  ! NOTE (found empirically during Task 1 TDD): u_transport stays EXACTLY
  ! UNITARY at each individual link even after this corruption (Phi is a
  ! diagonal phase, and a unitary times a unitary is unitary), and
  ! covariant_grad_block's per-shell terms are conjugations U*rho*U^H, which
  ! preserve operator norm for ANY unitary U -- so |Dq| is bounded by a FIXED
  ! constant times max|rho| REGARDLESS of whether U is gauge-covariant. A
  ! phase-only corruption can therefore NOT cause literal numerical blow-up
  ! (verified: two very different corruption magnitudes/frequencies both left
  ! rho_max/Jmax essentially unchanged from the correct run). What a broken
  ! Wilson-line telescoping genuinely breaks is GAUGE INVARIANCE itself (the
  ! relation u_transport_gauged(k) = W(k)^H u_transport(k) W(k+e) no longer
  ! holds once corrupted one-sidedly) -- run_teeth_corrupted_u3 below tests
  ! THAT (comparing a corrupt-then-gauge run against a gauge-then-corrupt run),
  ! which is what the U3f.4a/4b gauge-invariance checks would actually have
  ! caught.
  subroutine corrupt_transport_u3(gs)
    use gs_info_ssbe, only: s_sbe_gs_info
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    complex(8) :: Phi(nb, nb)
    real(8) :: t
    integer :: ik
    do ik = 1, nk
      t = dble(ik)
      Phi = (0d0, 0d0)
      ! rapidly-varying (high "frequency" in k, large amplitude), deliberately
      ! NOT slowly-varying/smooth -- maximally breaks the Wilson-line
      ! composition Um(k)=U(k)U(k+e)...U(k+(m-1)e) that covariant_grad_block's
      ! telescoping cancellation needs.
      Phi(1, 1) = exp(zi_ * 4.7d0 * sin(3.1d0 * t))
      Phi(2, 2) = exp(zi_ * 5.3d0 * sin(1.7d0 * t + 0.9d0))
      Phi(3, 3) = exp(zi_ * 4.1d0 * sin(2.3d0 * t + 1.7d0))
      Phi(4, 4) = exp(zi_ * 6.1d0 * sin(0.9d0 * t + 2.3d0))
      gs%u_transport(:, :, 1, ik) = matmul(Phi, gs%u_transport(:, :, 1, ik))
    end do
  end subroutine corrupt_transport_u3

  !======================= smooth Hermitian initial density ====================
  ! EQUAL diagonal populations across ALL FOUR bands (same common wiggle every
  ! band, honouring the X-full "equal occupation" fixture premise) with
  ! nonzero coherences between EVERY band pair.
  subroutine build_rho_init_u3(rho)
    implicit none
    complex(8), intent(out) :: rho(nb, nb, nk)
    real(8) :: th, wig
    integer :: ik, ib, jb
    rho = (0d0, 0d0)
    do ik = 1, nk
      th  = two_pi * dble(ik - 1) / dble(nk)
      wig = 0.05d0 * cos(th + 0.3d0)
      do ib = 1, nb
        rho(ib, ib, ik) = dcmplx(1.50d0 + wig, 0d0)
      end do
      do ib = 1, nb
        do jb = ib + 1, nb
          rho(ib, jb, ik) = (0.05d0 + 0.015d0 * dble(ib + jb)) * &
                            exp(zi_ * (0.4d0 * ib + 0.7d0 * jb + th))
          rho(jb, ib, ik) = conjg(rho(ib, jb, ik))
        end do
      end do
    end do
  end subroutine build_rho_init_u3

  ! inject physical rho into sbe%qnm AND sbe%qnm_new via the q_ij_from_rho bridge
  subroutine set_state_from_rho_u3(sbe, rho)
    use bloch_solver_ssbe, only: s_sbe_bloch_solver, q_ij_from_rho
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    complex(8), intent(in) :: rho(nb, nb, nk)
    integer :: ik, ib, jb
    complex(8) :: q
    do ik = sbe%ik_min, sbe%ik_max
      do ib = 1, nb
        do jb = 1, nb
          q = q_ij_from_rho(sbe, rho(ib, jb, ik), ib, jb, ik)
          sbe%qnm(ib, jb, ik)     = q
          sbe%qnm_new(ib, jb, ik) = q
        end do
      end do
    end do
  end subroutine set_state_from_rho_u3

  !======================= drive the real gicov propagator =====================
  subroutine run_gicov_u3(gs, rho_init, E, dt, tr, Jout, hn_max, all_ok, Jmax, rho_max, coh_drift)
    use gs_info_ssbe, only: s_sbe_gs_info
    use bloch_solver_ssbe, only: s_sbe_bloch_solver, init_sbe_bloch_solver, &
                                  prepare_qnm, dt_evolve_bloch_lg, adams_moulton_coefs, &
                                  calc_current_bloch_lg
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    complex(8), intent(in) :: rho_init(nb, nb, nk)
    real(8), intent(in) :: E(3), dt
    real(8), intent(out) :: tr(0:nsteps), Jout(3, 0:nsteps), hn_max, Jmax
    real(8), intent(out) :: rho_max, coh_drift
    logical, intent(out) :: all_ok
    type(s_sbe_bloch_solver) :: sbe
    real(8) :: bj_am(8, 8), hn, jm(3)
    complex(8) :: rho_now(nb, nb, nk), rho0(nb, nb, nk)
    integer :: it, icomm, ik, ib, jb

    icomm = 0
    call init_sbe_bloch_solver(sbe, gs, nb, icomm)
    sbe%flag_vnl_correction = .true.     ! activates the rvnl_tm_matrix contraction
    call prepare_qnm(sbe, gs, icomm)
    call adams_moulton_coefs(bj_am)       ! gicov branch ignores it; interface parity
    call set_state_from_rho_u3(sbe, rho_init)

    hn_max = 0d0; Jmax = 0d0; all_ok = .true.; rho_max = 0d0; coh_drift = 0d0
    tr(0) = trace_re_of_u3(sbe)
    call calc_current_bloch_lg(sbe, gs, jm, icomm); Jout(:, 0) = jm(:)
    if (.not. (is_finite(tr(0)) .and. is_finite(jm(1)) .and. is_finite(jm(2)) .and. is_finite(jm(3)))) &
      all_ok = .false.

    call rho_from_qnm_new_u3(sbe, rho_now)
    rho_max = max(rho_max, maxval(abs(rho_now)))
    rho0 = rho_now

    do it = 1, nsteps
      call dt_evolve_bloch_lg(sbe, gs, E, bj_am, dt, icomm)
      tr(it) = trace_re_of_u3(sbe)
      hn = herm_norm_of_u3(sbe); hn_max = max(hn_max, hn)
      call calc_current_bloch_lg(sbe, gs, jm, icomm); Jout(:, it) = jm(:)
      Jmax = max(Jmax, maxval(abs(jm)))
      call rho_from_qnm_new_u3(sbe, rho_now)
      rho_max = max(rho_max, maxval(abs(rho_now)))
      do ik = 1, nk
        do jb = 1, nb
          do ib = 1, nb
            if (ib /= jb) coh_drift = max(coh_drift, abs(rho_now(ib, jb, ik) - rho0(ib, jb, ik)))
          end do
        end do
      end do
      if (.not. (is_finite(tr(it)) .and. is_finite(hn) .and. &
                 is_finite(jm(1)) .and. is_finite(jm(2)) .and. is_finite(jm(3)))) all_ok = .false.
    end do
  end subroutine run_gicov_u3

  !======================= TEETH: broken-telescoping transport =================
  ! A phase-only corruption of u_transport CANNOT cause literal numerical
  ! blow-up (see corrupt_transport_u3's note: conjugation by any unitary
  ! preserves norm, so |Dq| stays boundedly proportional to max|rho|
  ! regardless of covariance). What it DOES break is the defining property
  ! U3f.4a/4b actually check: gauge invariance. If the SAME corrupted physics
  ! were genuinely gauge-covariant, a "corrupt-then-run" and a
  ! "gauge-then-corrupt-then-run" of the identical corruption Phi(k) would
  ! give the SAME Tr(t)/J(t) (same physics, different gauge) -- exactly as
  ! U3f.4a/4b verify for the UNCORRUPTED transport to 1e-8. This asserts that
  ! relation genuinely BREAKS (a large, non-roundoff residual) once the
  ! one-sided Phi(k) corruption is present, proving the gauge-invariance gate
  ! has teeth (it would have caught this corruption).
  subroutine run_teeth_corrupted_u3(nfail)
    use gs_info_ssbe, only: s_sbe_gs_info
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gst_u, gst_g
    complex(8) :: rho0(nb, nb, nk), rho0g(nb, nb, nk), W(nb, nb), Wh(nb, nb)
    real(8) :: Et(3), dtt
    real(8) :: trt_u(0:nsteps), Jt_u(3, 0:nsteps), hn_max_u, Jmax_u, rho_max_u, coh_drift_u
    real(8) :: trt_g(0:nsteps), Jt_g(3, 0:nsteps), hn_max_g, Jmax_g, rho_max_g, coh_drift_g
    logical :: ok_u, ok_g
    real(8) :: tr_resid, j_resid
    logical :: broke
    integer :: ik, it

    Et(1) = 0.10d0; Et(2) = 0d0; Et(3) = 0d0
    dtt = 0.05d0

    ! "corrupt-then-run" (no gauge applied at all)
    call build_gs_u3(gst_u)
    call corrupt_transport_u3(gst_u)
    call build_rho_init_u3(rho0)
    call run_gicov_u3(gst_u, rho0, Et, dtt, trt_u, Jt_u, hn_max_u, ok_u, Jmax_u, rho_max_u, coh_drift_u)

    ! "gauge-then-corrupt-then-run" -- SAME Phi(k) corruption, applied AFTER
    ! the (otherwise correct) gauge transform of u_transport/current operators
    ! and the matching rho0 -> W^H rho0 W. If corruption were gauge-neutral,
    ! this would reproduce trt_u/Jt_u exactly (as U3f.4a/4b prove it does for
    ! the UNCORRUPTED transport).
    call build_gs_u3(gst_g)
    call gauge_gs_inplace_u3(gst_g)
    call corrupt_transport_u3(gst_g)
    do ik = 1, nk
      W  = make_W_u3(ik)
      Wh = hconj_u3(W)
      rho0g(:, :, ik) = matmul(matmul(Wh, rho0(:, :, ik)), W)
    end do
    call run_gicov_u3(gst_g, rho0g, Et, dtt, trt_g, Jt_g, hn_max_g, ok_g, Jmax_g, rho_max_g, coh_drift_g)

    tr_resid = 0d0; j_resid = 0d0
    do it = 0, nsteps
      tr_resid = max(tr_resid, abs(trt_g(it) - trt_u(it)))
      j_resid  = max(j_resid, abs(Jt_g(1, it) - Jt_u(1, it)), &
                              abs(Jt_g(2, it) - Jt_u(2, it)), abs(Jt_g(3, it) - Jt_u(3, it)))
    end do

    broke = (.not. is_finite(tr_resid)) .or. (.not. is_finite(j_resid)) &
            .or. (tr_resid > 1d-4) .or. (j_resid > 1d-4)

    write(*, '(a,es10.2,a,es10.2,a,l1,a,l1)') &
      "      TEETH corrupted-transport: Tr-gauge-resid=", tr_resid, "  J-gauge-resid=", j_resid, &
      "  finite(u/g)=", ok_u, "/", ok_g

    call check_true(broke, &
      "U3f-TEETH: a non-covariant per-k transport corruption (broken Wilson-line " // &
      "telescoping) breaks gauge invariance of Tr/J (>1e-4, vs <1e-8 for the correct " // &
      "transport) -- the U3f.4a/4b gauge-invariance gate has teeth", nfail)
  end subroutine run_teeth_corrupted_u3

end program test_gicov_xfull
