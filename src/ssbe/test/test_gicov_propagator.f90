! src/ssbe/test/test_gicov_propagator.f90
! LG-SBE gicov Phase 3, Task 5b: U3 multi-step trace-conservation gate (THE
! acceptance test for Approach B').
!
! WHY THIS FILE EXISTS
!   Task 5a proved the gicov right-hand side gicov_rhs machine-exact for a SINGLE
!   evaluation (trace-preserving / Hermitian / gauge-covariant / non-vacuous, all
!   ~1e-16).  A correct RHS is necessary but not sufficient: a one-shot property
!   test cannot see an INTEGRATOR-STABILITY failure (the mirror lesson of
!   test_gifix_propagator.f90 -- a 1-step covariance check passed yet production
!   AB4 diverged at step 20).  This file closes that gap for gicov: it runs the
!   REAL production propagator (prepare_qnm + dt_evolve_bloch_lg with
!   sbe_lg_degen='gicov') for >=200 steps and asserts the whole gicov path stays
!   bounded end-to-end AND that its physical observables (trace, current) are
!   invariant under an arbitrary per-k GS gauge W(k) -- the real-Si failure mode
!   the plan targets.
!
!   The gicov branch of dt_evolve_bloch_lg (Task 5b) steps the PHYSICAL density
!   rho with a self-starting RK4 built on gicov_rhs (imaginary-axis stability
!   2*sqrt(2)~2.83, vs AB4's 0.43) -- NO dqnm_stock / Adams-Moulton history, and
!   NO bare Euler.  The GW collision is a separate VG-style substep (a no-op here:
!   t_2=1e30, GW off).  The R-1 qnm<->rho bridge (rho_ij_from_q / q_ij_from_rho,
!   exp_iphi time-constant) carries the same-block coherences in the PHYSICAL
!   basis, so stepping rho through qnm is lossless.
!
!   Fixture (mirrors test_gicov_rhs's fixture so the covariant term is engaged):
!     nb=4, nk=8; EXACTLY-degenerate composite block {2,3} at eigen 0.90 isolated
!     by a 0.60 gap from singletons band 1 (0.30) and band 4 (1.50).  EQUAL, FULL
!     occupation on the block (gicov fail-closes on asymmetric occupations, unlike
!     test_gifix_propagator).  t_2=1e30 (in-block dephasing is non-covariant;
!     production uses the GW collision).  A nonzero block Wilson transport keeps
!     the covariant intraband term active; the 0.60 out-of-block gap keeps the
!     dipole commutator active.  The initial rho is a smooth Hermitian density
!     with EQUAL within-block populations (honours the gicov premise) plus nonzero
!     in-block AND out-of-block coherences, so the arbitrary gauge W(k) transforms
!     it NON-trivially (rho -> W^H rho W rotates the coherences).
!
!   Assertions (THE gate):
!     1. Trace conserved   |Tr rho(t) - Nocc| < 1e-8 every step (boundedness).
!     2. Hermiticity       herm_norm(rho(t)) <= 1e-8 every step (no blow-up).
!     3. No NaN            all trace/current/herm finite every step.
!     4. Gauge invariance  Tr and current J(t) of the gauged run match the
!                          ungauged run to < 1e-8 at every step.
!     5. Non-vacuous       the gicov run actually MOVES (max|J| well above 0), so
!                          boundedness is not the boundedness of a frozen state.
!     6. Direct |rho| ceiling: max_t maxval(abs(rho(t))) < 10 (hardening, Task
!        5b review I-1). Checks 1/2/4 only trip at hard overflow (~1e308) or via
!        roundoff-coupled 1e-8 thresholds -- a regression that grew |rho| to
!        O(1e4-1e6) over 200 steps would silently PASS them all. This is a direct
!        ceiling, comfortably above the observed O(1) but far below that blind spot.
!     7. Direct |J| ceiling: max_t |J(t)| < 100 (same hardening motivation, I-1).
!     8. In-block coherence driven: max_{t,k} |rho23(t) - rho23(0)| > 1e-6
!        (hardening, Task 5b review I-2). Check 5 only lower-bounds the GLOBAL
!        max|J|, which is carried mostly by the out-of-block dipole -- it does not
!        by itself prove the covariant intraband term INSIDE the degenerate block
!        {2,3} is doing anything. This asserts the in-block coherence itself moves.
!     9. TEETH (negative control): the OLD uncapped in-block dipole path (a
!        near-degenerate pair with the naive dipole i*p/delta_omega, |d|~1.5e8,
!        driven through the AB4 path, sbe_lg_degen/='gicov') MUST diverge within
!        the same window -- proving the conservation checks above have teeth.
!
! BUILD (already-built ninja tree at build_local/; single-process communication
! dummy).  Compile from a CLEAN dir (a stale ./degenerate_block_ssbe.mod in the
! repo root shadows the fresh build_local one -- use -I build_local ONLY).  Link
! the SAME objects the salmon executable built, minus main.f90.o, via an
! @objs.txt response file (avoids the argv-length ld errno=63):
!
!   find <repo>/build_local/src/CMakeFiles/salmon.dir -name '*.o' ! -name 'main.f90.o' > objs.txt
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     -I<repo>/build_local -J<clean_dir> \
!     -c <repo>/src/ssbe/test/test_gicov_propagator.f90 -o <clean_dir>/test_gicov_propagator.o
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     @objs.txt <clean_dir>/test_gicov_propagator.o -o <clean_dir>/test_gicov_propagator \
!     -framework Accelerate -lm -ldl
!   <clean_dir>/test_gicov_propagator
!
program test_gicov_propagator
  use gs_info_ssbe,          only: s_sbe_gs_info
  use bloch_solver_ssbe,     only: s_sbe_bloch_solver, init_sbe_bloch_solver, &
                                    prepare_qnm, dt_evolve_bloch_lg, adams_moulton_coefs, &
                                    q_ij_from_rho, calc_current_bloch_lg
  use salmon_global,         only: epdir_re1, am_s, sbe_lg_diag, num_kgrid, t_2, &
                                    yn_sbe_gw_collision, sbe_deph_mode, sbe_lg_degen
  implicit none

  integer, parameter :: nb = 4, nk = 8
  integer, parameter :: nsteps = 200
  complex(8), parameter :: zi_ = (0d0, 1d0)
  real(8), parameter :: two_pi = 6.28318530717958647692d0
  ! block partition: band1 -> block 1; bands 2,3 -> block 2; band4 -> block 3
  integer, parameter :: blk(nb) = (/ 1, 2, 2, 3 /)
  integer :: nfail

  nfail = 0
  call set_globals()

  call test_gicov_bounded_and_gauge_invariant(nfail)   ! U3 acceptance gate
  call test_teeth_naive_dipole_diverges(nfail)         ! negative control (teeth)

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_gicov_propagator checks passed."
  end if

contains

  !======================= salmon_global fixture ==============================
  subroutine set_globals()
    implicit none
    epdir_re1(1) = 1d0; epdir_re1(2) = 0d0; epdir_re1(3) = 0d0
    am_s = 4                     ! matches adams_moulton_coefs' populated column (gicov ignores it)
    num_kgrid(1) = nk; num_kgrid(2) = 1; num_kgrid(3) = 1
    t_2 = 1d30                   ! in-block dephasing negligible (non-covariant; GW replaces it)
    sbe_lg_diag = 0              ! no diagnostic knockouts: full production physics
    sbe_lg_degen = 'gicov'
    yn_sbe_gw_collision = 'n'
    sbe_deph_mode = ''
  end subroutine set_globals

  !======================= assert / finite helpers ============================
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

  !======================= conservation probes (over sbe%qnm_new) =============
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

  !======================= physical rho probe (from qnm_new) ===================
  ! Task 5b review I-1/I-2 hardening: reconstructs the physical density rho AFTER
  ! a completed step directly from sbe%qnm_new, mirroring calc_current_bloch_lg's
  ! own gicov reconstruction (bloch_solver_ssbe.f90 ~1099-1119): rho(i,j) ==
  ! qnm_new(i,j) on the diagonal, == exp_iphi(i,j)*qnm_new(i,j) off it. exp_iphi is
  ! a TIME-CONSTANT unit phase set once by prepare_qnm (same-block exp_iphi==1 in
  ! gicov), so this is exact at every step. NOTE: this deliberately does NOT reuse
  ! rho_ij_from_q(sbe,...) -- that helper reads sbe%qnm, which after
  ! dt_evolve_bloch_lg (gicov) still holds the PREVIOUS step's state (step (0) of
  ! dt_evolve_bloch_lg_gicov only advances qnm_new -> qnm at the START of the NEXT
  ! call), so rho_ij_from_q would silently read one step stale.
  subroutine rho_from_qnm_new(sbe, rho_out)
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
  end subroutine rho_from_qnm_new

  !======================= 2x2 U(2) gauge on the {2,3} block ==================
  ! Block-diagonal wrt blk() by construction: exactly-unitary U(2) on {2,3}, plus
  ! independent phases on singleton bands 1 and 4.  (Same W as test_gicov_rhs.)
  function make_W(ik) result(W)
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
  end function make_W

  function hconj(A) result(Ah)
    implicit none
    complex(8), intent(in) :: A(nb, nb)
    complex(8) :: Ah(nb, nb)
    Ah = transpose(conjg(A))
  end function hconj

  integer function knext(ik, axis) result(kn)
    implicit none
    integer, intent(in) :: ik, axis
    if (axis == 1) then
      kn = mod(ik, nk) + 1
    else
      kn = ik                       ! num_kgrid=1 on axes 2,3: neighbour is self
    end if
  end function knext

  !======================= d_matrix from p (gicov rule) =======================
  ! same-block -> 0; out-of-block with |delta_omega|>eps -> i*p/delta_omega.
  subroutine build_d_from_p(gs)
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    real(8), parameter :: eps = 1d-12
    integer :: ik, ib, jb
    gs%d_matrix = (0d0, 0d0)
    do ik = 1, nk
      do ib = 1, nb
        do jb = 1, nb
          if (blk(ib) == blk(jb)) then
            gs%d_matrix(ib, jb, 1:3, ik) = (0d0, 0d0)
          else if (abs(gs%delta_omega(ib, jb, ik)) > eps) then
            gs%d_matrix(ib, jb, 1:3, ik) = zi_ * gs%p_mod_matrix(ib, jb, 1:3, ik) &
                                         & / gs%delta_omega(ib, jb, ik)
          else
            gs%d_matrix(ib, jb, 1:3, ik) = (0d0, 0d0)
          end if
        end do
      end do
    end do
  end subroutine build_d_from_p

  !======================= synthetic gicov gs fixture =========================
  ! Exactly-degenerate composite block {2,3}; singletons {1},{4}; 0.60 gap.
  ! Hermitian momentum p on axis 1 only (epdir=(1,0,0)); nontrivial 2x2 block
  ! Wilson transport on axis 1.  p_tm_matrix := p_mod_matrix (same Hermitian
  ! velocity) so the current is Tr(p rho) with the SAME matrix the dipole is
  ! built from -> gauge-invariant under a consistent W.  Equal, FULL occupation
  ! on the block (gicov premise).
  subroutine build_gs(gs)
    implicit none
    type(s_sbe_gs_info), intent(out) :: gs
    real(8) :: eigen(nb), t, ang, c, s, phz
    integer :: ik, ib, jb

    gs%nk = nk; gs%nb = nb; gs%ne = 6
    allocate(gs%eigen(nb, nk), gs%occup(nb, nk), gs%kweight(nk))
    allocate(gs%delta_omega(nb, nb, nk))
    allocate(gs%p_mod_matrix(nb, nb, 3, nk))
    allocate(gs%p_tm_matrix(nb, nb, 3, nk))
    allocate(gs%d_matrix(nb, nb, 3, nk))
    allocate(gs%u_transport(nb, nb, 3, nk))
    allocate(gs%block_id(nb, nk))
    allocate(gs%bvec(3, 3))

    eigen(1) = 0.30d0
    eigen(2) = 0.90d0          ! block {2,3}: EXACTLY degenerate
    eigen(3) = 0.90d0
    eigen(4) = 1.50d0

    gs%volume = 1d0
    gs%b_matrix = 0d0
    gs%b_matrix(1, 1) = two_pi
    gs%b_matrix(2, 2) = two_pi
    gs%b_matrix(3, 3) = two_pi

    gs%nbvec = 3
    gs%bvec(:, 1) = (/ 1, 0, 0 /)
    gs%bvec(:, 2) = (/ 0, 1, 0 /)
    gs%bvec(:, 3) = (/ 0, 0, 1 /)

    gs%p_mod_matrix = (0d0, 0d0)
    gs%u_transport  = (0d0, 0d0)

    do ik = 1, nk
      t = dble(ik - 1)
      gs%kweight(ik) = 1d0
      gs%eigen(:, ik) = eigen(:)
      gs%block_id(:, ik) = blk(:)
      ! equal full occupation on the block {2,3}; band1 full, band4 empty
      gs%occup(1, ik) = 2d0; gs%occup(2, ik) = 2d0
      gs%occup(3, ik) = 2d0; gs%occup(4, ik) = 0d0
      do jb = 1, nb
        do ib = 1, nb
          gs%delta_omega(ib, jb, ik) = gs%eigen(ib, ik) - gs%eigen(jb, ik)
        end do
      end do

      ! Hermitian momentum on axis 1 (upper triangle set, lower = conjg)
      do ib = 1, nb
        do jb = ib + 1, nb
          phz = 0.3d0 * dble(ib) + 0.5d0 * dble(jb) + 0.2d0 * t
          gs%p_mod_matrix(ib, jb, 1, ik) = (0.4d0 + 0.1d0 * dble(ib + jb)) * exp(zi_ * phz)
          gs%p_mod_matrix(jb, ib, 1, ik) = conjg(gs%p_mod_matrix(ib, jb, 1, ik))
        end do
      end do

      ! block Wilson transport on axis 1: 2x2 U(2) on {2,3}, identity elsewhere
      gs%u_transport(1, 1, 1, ik) = (1d0, 0d0)
      gs%u_transport(4, 4, 1, ik) = (1d0, 0d0)
      ang = 0.30d0 + 0.15d0 * t
      c = cos(ang); s = sin(ang)
      phz = 0.10d0 + 0.05d0 * t
      gs%u_transport(2, 2, 1, ik) = c * exp( zi_ * phz)
      gs%u_transport(2, 3, 1, ik) = s
      gs%u_transport(3, 2, 1, ik) = -s
      gs%u_transport(3, 3, 1, ik) = c * exp(-zi_ * phz)
      do ib = 1, nb
        gs%u_transport(ib, ib, 2, ik) = (1d0, 0d0)
        gs%u_transport(ib, ib, 3, ik) = (1d0, 0d0)
      end do
    end do

    gs%p_tm_matrix = gs%p_mod_matrix    ! current velocity == dipole velocity (gauge-consistent)
    call build_d_from_p(gs)
  end subroutine build_gs

  !======================= apply the arbitrary per-k gauge W(k) ================
  ! p_mod, p_tm -> W(k)^H p W(k) (per axis); u_transport -> Wilson-line transform
  ! W(k)^H U(k,axis) W(k+e_axis); d_matrix rebuilt from the gauged p_mod.
  subroutine gauge_gs_inplace(gs)
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    complex(8) :: W(nb, nb), Wh(nb, nb)
    integer :: ik, axis, kk
    do ik = 1, nk
      W  = make_W(ik)
      Wh = hconj(W)
      do axis = 1, 3
        gs%p_mod_matrix(:, :, axis, ik) = matmul(matmul(Wh, gs%p_mod_matrix(:, :, axis, ik)), W)
        gs%p_tm_matrix (:, :, axis, ik) = matmul(matmul(Wh, gs%p_tm_matrix (:, :, axis, ik)), W)
      end do
    end do
    do ik = 1, nk
      Wh = hconj(make_W(ik))
      do axis = 1, 3
        kk = knext(ik, axis)
        gs%u_transport(:, :, axis, ik) = matmul(matmul(Wh, gs%u_transport(:, :, axis, ik)), make_W(kk))
      end do
    end do
    call build_d_from_p(gs)
  end subroutine gauge_gs_inplace

  !======================= smooth Hermitian initial density ===================
  ! EQUAL within-block populations rho22==rho33 (honours the gicov premise) with
  ! nonzero in-block coherence rho23 and out-of-block coherences, so W^H rho W is
  ! non-trivial (the coherences rotate under the gauge).
  subroutine build_rho_init(rho)
    implicit none
    complex(8), intent(out) :: rho(nb, nb, nk)
    real(8) :: th, phz, blkpop
    integer :: ik, ib, jb
    rho = (0d0, 0d0)
    do ik = 1, nk
      th = two_pi * dble(ik - 1) / dble(nk)
      blkpop = 0.70d0 + 0.04d0 * cos(th + 0.30d0)
      rho(1, 1, ik) = dcmplx(0.90d0 + 0.05d0 * cos(th),        0d0)
      rho(2, 2, ik) = dcmplx(blkpop,                            0d0)   ! block: EQUAL
      rho(3, 3, ik) = dcmplx(blkpop,                            0d0)   ! block: EQUAL
      rho(4, 4, ik) = dcmplx(0.15d0 + 0.05d0 * sin(th + 1.0d0), 0d0)
      do ib = 1, nb
        do jb = ib + 1, nb
          phz = 0.4d0 * dble(ib) + 0.7d0 * dble(jb) + th
          rho(ib, jb, ik) = (0.06d0 + 0.02d0 * dble(ib + jb)) * exp(zi_ * phz)
          rho(jb, ib, ik) = conjg(rho(ib, jb, ik))
        end do
      end do
    end do
  end subroutine build_rho_init

  ! inject physical rho into sbe%qnm AND sbe%qnm_new via the q_ij_from_rho bridge
  ! (qnm_new is the state dt_evolve_bloch_lg reads at step start; qnm mirrors it).
  subroutine set_state_from_rho(sbe, rho)
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
  end subroutine set_state_from_rho

  !======================= drive the real gicov propagator ====================
  ! init -> prepare_qnm -> inject rho -> nsteps of dt_evolve_bloch_lg (gicov RK4).
  ! Records trace(0:nsteps) and current(1:3, 0:nsteps); tracks max herm_norm,
  ! whether everything stayed finite, and the largest |current| seen.
  subroutine run_gicov(gs, rho_init, E, dt, tr, Jout, hn_max, all_ok, Jmax, rho_max, rho23_drift)
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    complex(8), intent(in) :: rho_init(nb, nb, nk)
    real(8), intent(in) :: E(3), dt
    real(8), intent(out) :: tr(0:nsteps), Jout(3, 0:nsteps), hn_max, Jmax
    real(8), intent(out) :: rho_max, rho23_drift    ! Task 5b review I-1/I-2 hardening
    logical, intent(out) :: all_ok
    type(s_sbe_bloch_solver) :: sbe
    real(8) :: bj_am(8, 8), hn, jm(3)
    complex(8) :: rho_now(nb, nb, nk), rho23_0(nk)
    integer :: it, icomm, ik

    icomm = 0
    call init_sbe_bloch_solver(sbe, gs, nb, icomm)
    call prepare_qnm(sbe, gs, icomm)
    call adams_moulton_coefs(bj_am)          ! gicov branch ignores bj_am; passed for the interface
    call set_state_from_rho(sbe, rho_init)

    hn_max = 0d0; Jmax = 0d0; all_ok = .true.; rho_max = 0d0; rho23_drift = 0d0
    tr(0) = trace_re_of(sbe)
    call calc_current_bloch_lg(sbe, gs, jm, icomm); Jout(:, 0) = jm(:)
    if (.not. (is_finite(tr(0)) .and. is_finite(jm(1)) .and. is_finite(jm(2)) .and. is_finite(jm(3)))) &
      all_ok = .false.

    call rho_from_qnm_new(sbe, rho_now)                ! t=0 physical rho (== rho_init, round-tripped)
    rho_max = max(rho_max, maxval(abs(rho_now)))
    do ik = 1, nk
      rho23_0(ik) = rho_now(2, 3, ik)                  ! in-block coherence baseline, all k
    end do

    do it = 1, nsteps
      call dt_evolve_bloch_lg(sbe, gs, E, bj_am, dt, icomm)
      tr(it) = trace_re_of(sbe)
      hn = herm_norm_of(sbe); hn_max = max(hn_max, hn)
      call calc_current_bloch_lg(sbe, gs, jm, icomm); Jout(:, it) = jm(:)
      Jmax = max(Jmax, maxval(abs(jm)))
      call rho_from_qnm_new(sbe, rho_now)
      rho_max = max(rho_max, maxval(abs(rho_now)))
      do ik = 1, nk
        rho23_drift = max(rho23_drift, abs(rho_now(2, 3, ik) - rho23_0(ik)))
      end do
      if (.not. (is_finite(tr(it)) .and. is_finite(hn) .and. &
                 is_finite(jm(1)) .and. is_finite(jm(2)) .and. is_finite(jm(3)))) all_ok = .false.
    end do
  end subroutine run_gicov

  !======================= U3 GATE: bounded + gauge invariant =================
  subroutine test_gicov_bounded_and_gauge_invariant(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs, gsg
    complex(8) :: rho(nb, nb, nk), rhog(nb, nb, nk), W(nb, nb), Wh(nb, nb)
    real(8) :: E(3), dt, Nocc
    real(8) :: tr_u(0:nsteps), tr_g(0:nsteps), Ju(3, 0:nsteps), Jg(3, 0:nsteps)
    real(8) :: hn_u, hn_g, Jmax_u, Jmax_g
    real(8) :: rho_max_u, rho_max_g, rho23_drift_u, rho23_drift_g   ! I-1/I-2 hardening
    real(8) :: tr_drift, gauge_tr, gauge_J, rho_max_all, jmax_all
    logical :: ok_u, ok_g
    integer :: it, ik

    E(1) = 0.10d0; E(2) = 0d0; E(3) = 0d0
    dt = 0.05d0

    ! ---- ungauged run ----
    call build_gs(gs)
    call build_rho_init(rho)
    call run_gicov(gs, rho, E, dt, tr_u, Ju, hn_u, ok_u, Jmax_u, rho_max_u, rho23_drift_u)
    Nocc = tr_u(0)

    ! ---- gauged run: same physics under an arbitrary per-k block gauge W(k) ----
    call build_gs(gsg)
    call gauge_gs_inplace(gsg)
    do ik = 1, nk
      W  = make_W(ik)
      Wh = hconj(W)
      rhog(:, :, ik) = matmul(matmul(Wh, rho(:, :, ik)), W)     ! rho -> W^H rho W
    end do
    call run_gicov(gsg, rhog, E, dt, tr_g, Jg, hn_g, ok_g, Jmax_g, rho_max_g, rho23_drift_g)

    ! ---- assertion 1: trace conserved to 1e-8 every step (ungauged & gauged) ---
    tr_drift = 0d0
    do it = 0, nsteps
      tr_drift = max(tr_drift, abs(tr_u(it) - Nocc), abs(tr_g(it) - tr_g(0)))
    end do
    call check_true(is_finite(tr_drift) .and. tr_drift < 1d-8, &
      "U3.1: |Tr rho(t) - Nocc| < 1e-8 over 200 gicov RK4 steps (trace conserved)", nfail)

    ! ---- assertion 2: Hermiticity bounded <= 1e-8 every step -------------------
    call check_true(is_finite(hn_u) .and. is_finite(hn_g) .and. hn_u <= 1d-8 .and. hn_g <= 1d-8, &
      "U3.2: herm_norm(rho(t)) <= 1e-8 over 200 steps (no blow-up)", nfail)

    ! ---- assertion 3: no NaN/Inf anywhere --------------------------------------
    call check_true(ok_u .and. ok_g, &
      "U3.3: all trace/current/herm finite (no NaN) over 200 steps", nfail)

    ! ---- assertion 4: gauge invariance of Tr and current J(t) ------------------
    gauge_tr = 0d0; gauge_J = 0d0
    do it = 0, nsteps
      gauge_tr = max(gauge_tr, abs(tr_g(it) - tr_u(it)))
      gauge_J  = max(gauge_J,  abs(Jg(1, it) - Ju(1, it)), &
                               abs(Jg(2, it) - Ju(2, it)), &
                               abs(Jg(3, it) - Ju(3, it)))
    end do
    call check_true(is_finite(gauge_tr) .and. gauge_tr < 1d-8, &
      "U3.4a: Tr(gauged) == Tr(ungauged) to 1e-8 every step (gauge-invariant trace)", nfail)
    call check_true(is_finite(gauge_J) .and. gauge_J < 1d-8, &
      "U3.4b: J(gauged) == J(ungauged) to 1e-8 every step (gauge-invariant current)", nfail)

    ! ---- assertion 5: non-vacuous (the propagator actually moves) --------------
    call check_true(Jmax_u > 1d-6, &
      "U3.5: max|J(t)| >> 0 (dynamics are non-trivial, boundedness not vacuous)", nfail)

    ! ---- assertion 6: DIRECT |rho| upper bound (hardening, Task 5b review I-1) -
    ! Checks 1/2/4 only trip at hard overflow (~1e308) or via roundoff-coupled
    ! 1e-8 thresholds -- a regression that grew |rho| to O(1e4-1e6) over the
    ! 200-step window would silently satisfy ALL of them. Ceiling 10 is chosen
    ! comfortably above the observed O(1) magnitude but far below that blind spot.
    rho_max_all = max(rho_max_u, rho_max_g)
    call check_true(is_finite(rho_max_all) .and. rho_max_all < 10d0, &
      "U3.6: |rho| bounded (max_t maxval(abs(rho(t))) < 10, no mild/moderate blow-up)", nfail)

    ! ---- assertion 7: DIRECT |J| upper bound (hardening, Task 5b review I-1) ---
    jmax_all = max(Jmax_u, Jmax_g)
    call check_true(is_finite(jmax_all) .and. jmax_all < 100d0, &
      "U3.7: |J| bounded above (max_t |J(t)| < 100, direct ceiling not just Tr/herm-coupled)", nfail)

    ! ---- assertion 8: in-block coherence rho(2,3,:) genuinely DRIVEN -----------
    ! (hardening, Task 5b review I-2). Assertion 5 only lower-bounds the GLOBAL
    ! max|J|, which is carried mostly by the out-of-block dipole (band1/4 <->
    ! block{2,3}) -- it does not by itself prove the covariant intraband term
    ! INSIDE the degenerate block is doing anything. This asserts the in-block
    ! coherence rho(2,3,ik) itself moves (max over k of |rho23(t)-rho23(0)|,
    ! well above roundoff), i.e. the covariant transport term is genuinely
    ! engaged, not vacuous.
    call check_true(is_finite(rho23_drift_u) .and. rho23_drift_u > 1d-6, &
      "U3.8: max_t,k |rho23(t) - rho23(0)| > 1e-6 (in-block coherence genuinely driven)", nfail)

    write(*, '(a,es10.2,a,es10.2,a,es10.2)') &
      "      trace-drift=", tr_drift, "  herm(u/g)=", max(hn_u, hn_g), "  |J|max=", Jmax_u
    write(*, '(a,es10.2,a,es10.2)') &
      "      gauge-resid  Tr=", gauge_tr, "  J=", gauge_J
    write(*, '(a,es10.2,a,es10.2)') &
      "      rho_max(u/g)=", rho_max_all, "  rho23-drift(u)=", rho23_drift_u
  end subroutine test_gicov_bounded_and_gauge_invariant

  !======================= near-degenerate gs for the teeth ===================
  ! A near- (not exactly) degenerate pair {2,3} at delta_omega=2e-9 with the
  ! NAIVE dipole d=i*p/delta_omega (|d|~1.5e8) and ASYMMETRIC occupation (drives
  ! the coherence).  Singletons {1},{4} inert (d=0).  This is the OLD uncapped
  ! in-block dipole path the AB4 propagator runs when sbe_lg_degen/='gicov'.
  subroutine build_gs_teeth(gs, domega, d_naive)
    implicit none
    type(s_sbe_gs_info), intent(out) :: gs
    real(8), intent(in) :: domega
    complex(8), intent(in) :: d_naive
    real(8) :: eigen(nb)
    integer :: ik, ib, jb

    gs%nk = nk; gs%nb = nb; gs%ne = 6
    allocate(gs%eigen(nb, nk), gs%occup(nb, nk), gs%kweight(nk))
    allocate(gs%delta_omega(nb, nb, nk))
    allocate(gs%d_matrix(nb, nb, 3, nk))

    eigen(1) = 0.30d0
    eigen(2) = 0.90d0 + 0.5d0 * domega
    eigen(3) = 0.90d0 - 0.5d0 * domega
    eigen(4) = 1.50d0

    gs%volume = 1d0
    gs%b_matrix = 0d0
    gs%b_matrix(1, 1) = two_pi; gs%b_matrix(2, 2) = two_pi; gs%b_matrix(3, 3) = two_pi

    gs%d_matrix = (0d0, 0d0)
    do ik = 1, nk
      gs%kweight(ik) = 1d0
      gs%eigen(:, ik) = eigen(:)
      gs%occup(1, ik) = 2d0; gs%occup(2, ik) = 2d0    ! ASYMMETRIC across {2,3}
      gs%occup(3, ik) = 0d0; gs%occup(4, ik) = 0d0
      do jb = 1, nb
        do ib = 1, nb
          gs%delta_omega(ib, jb, ik) = eigen(ib) - eigen(jb)
        end do
      end do
      gs%d_matrix(2, 3, 1, ik) = d_naive
      gs%d_matrix(3, 2, 1, ik) = conjg(d_naive)
    end do
  end subroutine build_gs_teeth

  !======================= TEETH: naive-dipole AB4 diverges ===================
  subroutine test_teeth_naive_dipole_diverges(nfail)
    implicit none
    integer, intent(inout) :: nfail
    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver) :: sbe
    real(8) :: E(3), bj_am(8, 8), dt, domega, tr0, tr, hn
    complex(8) :: p23, d_naive
    integer :: it, icomm, it_diverged
    logical :: diverged

    sbe_lg_degen = 'off'                    ! OLD AB4 path (NOT the gicov RK4 branch)
    domega  = 2d-9
    p23     = dcmplx(0.3d0, 0.1d0)
    d_naive = zi_ * p23 / domega            ! naive dipole i*p/delta_omega, |d|~1.5e8
    call check_true(abs(d_naive) > 1d7, &
      "teeth fixture: |i*p/domega| >> 1 (singular naive dipole, ~1.5e8)", nfail)

    call build_gs_teeth(gs, domega, d_naive)

    E(1) = 1d-6; E(2) = 0d0; E(3) = 0d0
    dt = 0.05d0
    icomm = 0

    call init_sbe_bloch_solver(sbe, gs, nb, icomm)
    call prepare_qnm(sbe, gs, icomm)
    call adams_moulton_coefs(bj_am)
    tr0 = trace_re_of(sbe)

    diverged = .false.; it_diverged = -1
    do it = 1, nsteps
      call dt_evolve_bloch_lg(sbe, gs, E, bj_am, dt, icomm)
      tr = trace_re_of(sbe)
      hn = herm_norm_of(sbe)
      if (.not. is_finite(tr) .or. .not. is_finite(hn) .or. &
          hn > 1d-8 .or. abs(tr - tr0) > 1d-6) then
        diverged = .true.; it_diverged = it; exit
      end if
    end do
    if (diverged) write(*, '(a,i0,a)') "      (teeth diverged at AB4 step ", it_diverged, ")"
    call check_true(diverged, &
      "TEETH: naive-dipole (|d|~1.5e8) blows up the SAME AB4 propagator gicov keeps bounded", nfail)

    sbe_lg_degen = 'gicov'                  ! restore
  end subroutine test_teeth_naive_dipole_diverges

end program test_gicov_propagator
