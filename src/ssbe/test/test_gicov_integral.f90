! src/ssbe/test/test_gicov_integral.f90
!
! Standalone unit test of the integral (covariant-Houston) transport kernels in
! gicov_integral_ssbe.f90, against a 2-band Dirac-cone analytic model -- the
! SAME decisive checks the numpy toy (toy_dirac_integral.py) validated, mapped
! onto the SALMON kernels.  This module has NO SALMON dependency (only LAPACK +
! the pure degenerate_block_ssbe), so it builds WITHOUT the full salmon object
! tree (unlike the propagator/checker tests that link build_local/):
!
!   gfortran -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     src/ssbe/sbe_lg_mode_ssbe.f90 src/ssbe/degenerate_block_ssbe.f90 \
!     src/ssbe/gicov_integral_ssbe.f90 \
!     src/ssbe/test/test_gicov_integral.f90 -o t -framework Accelerate && ./t
!
! (sbe_lg_mode_ssbe.f90 is a REQUIRED compile unit: degenerate_block_ssbe uses
! its predicates, and T8 below asserts the mode truth table directly.)
!
! Covers (task test list 1,2,4 + the SDD risk register):
!   T1 transport orientation + sign (risk #1): W O W^dagger built from an
!      ordered link chain reproduces the exact analytic co-moving Hamiltonian
!      OFF-AXIS; the wrong sandwich / wrong sign (x=kappa+a) does NOT.
!   T2 endpoint interpolation exactness (linear H) and floor() bracketing.
!   T3 step conserves trace / Hermiticity / populations (machine level).
!   T4 T2 covariance (risk #3/#7): exact-degenerate pair -> g(0)=0 -> NO
!      dephasing; a gapped pair decays by exp(-dt/T2).
!   T5 cache sizing (risk #5/#8): j_max=13 for graphene k99; int64 cache bytes
!      = 80(2j+1)nb^2 nk_local with no overflow at production scale (FIVE
!      complex matrices per shift: H~, the three v~, and the co-moving f0
!      reference F~0 -- see T12).
!   T6 single-axis / linear-polarization guard (risk #6): [110] and circular
!      are rejected on the whole trajectory; a single axis is accepted.
!   T7 floor() vs int() bracketing for a>0 (negative mesh shift).
!   T8 mode-predicate truth table (risk #2): the LEGACY modes (off/gi/gifix/
!      gicov) never select the integral machinery, and 'gicov_int' never
!      selects the finite-difference one -- the static half of the legacy
!      byte-regression guard (the full byte compare is the Fugaku gate).
!   T9 eigensolver-failure CONTRACT: a broken instantaneous eigenproblem is
!      REPORTED through the status flag, never absorbed.  Both kernels used to
!      swallow it (step_k froze the block = rho_out=rho; occupation_k fell back
!      to the FORBIDDEN diag(rho~)) -- and because both fallbacks conserve the
!      trace, the Ne/trace monitor could not see it.  Note zheev itself returns
!      info=0 on a NaN-poisoned matrix (verified on Accelerate), so the status
!      must ALSO reject non-finite eigenvalues, or the failure stays silent.
!  T10 bracket at the cached-span EDGE.  q = -j_max lands exactly on the node
!      +j_max (s = -q/dk), so the old floor()+"n_int+1" bracket asked for node
!      j_max+1: one past the cache -- an abort (or, without the guard, an
!      out-of-bounds read) at a perfectly legal endpoint of the pulse.  A zero
!      fractional weight must COLLAPSE the bracket onto the node instead.
!  T11 near-degenerate CHAIN closure.  "|eps_a-eps_b| <= floor" is not an
!      equivalence relation: for a chain a~b~c the pairwise gate leaves the
!      END pair (a,c) dephasing, which makes the result depend on the basis
!      chosen INSIDE the connected near-degenerate manifold -- the very thing
!      the floor exists to prevent.  The gate must act on the CONNECTED BLOCK.
!      T11 also pins the basis-invariant readout: within a degenerate block the
!      individual eigen-columns are arbitrary, so only the block total is
!      observable (gicov_int_occupation_k reports it shared equally).
!  T12 co-moving f0(x) REFERENCE.  _sbe_occ.data is an EXCITATION (the
!      non-gicov_int sibling reports rho_bb - gs%occup); the integral path was
!      reporting the ABSOLUTE instantaneous population, so an unexcited state
!      read out as "fully occupied" rather than as zero excitation.  The
!      reference must be the GS occupation TRANSPORTED into the same frame,
!      F~0 = W diag(f0(k_remote)) W^dagger, not a rigid "1 below nvb" block:
!      T12 shows the rigid reference manufactures a spurious excitation for a
!      METAL (fractional, k-varying f0), which is why F~0 is cached.
!
program test_gicov_integral
  use gicov_integral_ssbe
  implicit none
  integer :: nfail
  nfail = 0

  call t1_transport_sign(nfail)
  call t2_interp_floor(nfail)
  call t3_conservation(nfail)
  call t4_t2_covariance(nfail)
  call t5_cache_sizing(nfail)
  call t6_axis_guard(nfail)
  call t7_floor_shift(nfail)
  call t8_mode_predicates(nfail)
  call t9_eigensolver_status(nfail)
  call t10_bracket_edge(nfail)
  call t11_degenerate_closure(nfail)
  call t12_comoving_f0(nfail)

  if (nfail == 0) then
    write(*, '(a)') "ALL PASS (test_gicov_integral)"
  else
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    error stop 1
  end if

contains

  subroutine report(name, ok, nfail)
    character(*), intent(in)    :: name
    logical,      intent(in)    :: ok
    integer,      intent(inout) :: nfail
    if (ok) then
      write(*, '(a,a)') "  PASS  ", name
    else
      write(*, '(a,a)') "  FAIL  ", name
      nfail = nfail + 1
    end if
  end subroutine report

  ! 2-band Dirac eigendata: H(q) = vF(qx sx + qy sy); col1 lower (-vF r),
  ! col2 upper (+vF r).  Same convention as the toy's analytic_eig.
  subroutine eig2(qx, qy, w, V)
    real(8),    intent(in)  :: qx, qy
    real(8),    intent(out) :: w(2)
    complex(8), intent(out) :: V(2, 2)
    real(8), parameter :: vF = 1d0
    real(8)    :: r, th, s2
    complex(8) :: eph
    r = hypot(qx, qy)
    th = atan2(qy, qx)
    eph = exp(cmplx(0d0, th, 8))
    s2 = 1d0 / sqrt(2d0)
    V(1, 1) =  cmplx(s2, 0d0, 8);  V(2, 1) = -s2 * eph
    V(1, 2) =  cmplx(s2, 0d0, 8);  V(2, 2) =  s2 * eph
    w(1) = -vF * r;  w(2) = vF * r
  end subroutine eig2

  subroutine horb2(qx, qy, H)
    real(8),    intent(in)  :: qx, qy
    complex(8), intent(out) :: H(2, 2)
    real(8), parameter :: vF = 1d0
    H(1, 1) = (0d0, 0d0);  H(1, 2) = vF * cmplx(qx, -qy, 8)
    H(2, 1) = vF * cmplx(qx, qy, 8);  H(2, 2) = (0d0, 0d0)
  end subroutine horb2

  ! forward single-step link L(m) = polar(V(m)^dagger V(m+1)) between the mesh
  ! nodes m and m+1 (built via the SAME polar the kernel chain uses).
  subroutine flink(kx0, ky0, dk, m, L)
    use degenerate_block_ssbe, only: polar_unitary
    real(8),    intent(in)  :: kx0, ky0, dk
    integer,    intent(in)  :: m
    complex(8), intent(out) :: L(2, 2)
    complex(8) :: Va(2, 2), Vb(2, 2), Movl(2, 2)
    real(8)    :: wa(2), wb(2), sig
    integer    :: ierr
    call eig2(kx0 + m * dk, ky0, wa, Va)
    call eig2(kx0 + (m + 1) * dk, ky0, wb, Vb)
    Movl = matmul(conjg(transpose(Va)), Vb)
    call polar_unitary(Movl, 2, L, sig, ierr)
  end subroutine flink

  real(8) function maxdiff(A, B, n)
    integer,    intent(in) :: n
    complex(8), intent(in) :: A(n, n), B(n, n)
    integer :: i, j
    maxdiff = 0d0
    do j = 1, n
      do i = 1, n
        maxdiff = max(maxdiff, abs(A(i, j) - B(i, j)))
      end do
    end do
  end function maxdiff

  !----------------------------------------------------------------
  subroutine t1_transport_sign(nfail)
    integer, intent(inout) :: nfail
    real(8), parameter :: kx0 = 0.30d0, ky0 = 0.17d0, dk = 0.05d0
    complex(8) :: link(2, 2, 3), W(2, 2), Oloc(2, 2), Y(2, 2), Yref(2, 2)
    complex(8) :: Vk(2, 2), Vr(2, 2), Hr(2, 2)
    real(8)    :: wk(2), wr(2), err
    integer    :: n, s
    logical    :: ok
    call eig2(kx0, ky0, wk, Vk)                  ! frame at kappa = node 0
    ok = .true.
    ! forward telescoping n = 1,2,3
    do n = 1, 3
      do s = 1, n
        call flink(kx0, ky0, dk, s - 1, link(:, :, s))
      end do
      call gicov_int_chain(link(:, :, 1:n), 2, n, .false., W)
      ! band-frame operator at the remote node = diag(band energies there)
      call eig2(kx0 + n * dk, ky0, wr, Vr)
      Oloc = (0d0, 0d0);  Oloc(1, 1) = wr(1);  Oloc(2, 2) = wr(2)
      call gicov_int_transport_op(W, Oloc, 2, Y)
      ! reference co-moving operator V(kappa)^dagger H_orb(remote) V(kappa)
      call horb2(kx0 + n * dk, ky0, Hr)
      Yref = matmul(matmul(conjg(transpose(Vk)), Hr), Vk)
      err = maxdiff(Y, Yref, 2)
      if (err > 1d-11) ok = .false.
    end do
    call report("T1 forward transport orientation matches exact co-moving H (off-axis)", ok, nfail)

    ! backward n = -1: W(node0,node0-1) = adjoint of the (node-1 -> node0) link
    call flink(kx0, ky0, dk, -1, link(:, :, 1))
    call gicov_int_chain(link(:, :, 1:1), 2, 1, .true., W)
    call eig2(kx0 - dk, ky0, wr, Vr)
    Oloc = (0d0, 0d0);  Oloc(1, 1) = wr(1);  Oloc(2, 2) = wr(2)
    call gicov_int_transport_op(W, Oloc, 2, Y)
    call horb2(kx0 - dk, ky0, Hr)
    Yref = matmul(matmul(conjg(transpose(Vk)), Hr), Vk)
    call report("T1 backward (adjoint) transport matches exact co-moving H", &
              & maxdiff(Y, Yref, 2) < 1d-11, nfail)

    ! wrong sandwich W^dagger O W is O(dk) off-axis (decisive sign/orientation)
    do s = 1, 2
      call flink(kx0, ky0, dk, s - 1, link(:, :, s))
    end do
    call gicov_int_chain(link(:, :, 1:2), 2, 2, .false., W)
    call eig2(kx0 + 2 * dk, ky0, wr, Vr)
    Oloc = (0d0, 0d0);  Oloc(1, 1) = wr(1);  Oloc(2, 2) = wr(2)
    Y = matmul(matmul(conjg(transpose(W)), Oloc), W)     ! WRONG orientation
    call horb2(kx0 + 2 * dk, ky0, Hr)
    Yref = matmul(matmul(conjg(transpose(Vk)), Hr), Vk)
    call report("T1 wrong sandwich is detectably WRONG off-axis (guards orientation)", &
              & maxdiff(Y, Yref, 2) > 1d-3, nfail)
  end subroutine t1_transport_sign

  !----------------------------------------------------------------
  ! Endpoint interpolation of the bounded transported operator reproduces the
  ! exact co-moving Hamiltonian at x = kappa - a for arbitrary a (H_orb linear
  ! -> interpolation exact); the SIGN (x = kappa - a) is what makes it match.
  subroutine t2_interp_floor(nfail)
    integer, intent(inout) :: nfail
    real(8), parameter :: kx0 = 0.30d0, ky0 = 0.17d0, dk = 0.05d0
    complex(8) :: link(2, 2, 4), Wlo(2, 2), Whi(2, 2)
    complex(8) :: Ylo(2, 2), Yhi(2, 2), Y(2, 2), Ywrong(2, 2), Yref(2, 2), Yrefw(2, 2)
    complex(8) :: Vk(2, 2), Vr(2, 2), Hr(2, 2), Oloc(2, 2)
    real(8)    :: wk(2), wr(2), a, frac
    integer    :: n_int
    call eig2(kx0, ky0, wk, Vk)
    a = 0.6d0 * dk                                  ! x = kappa - a  falls between nodes -1 and 0
    call gicov_int_floor_shift(a, dk, n_int, frac)  ! s = -0.6 -> n_int=-1, frac=0.4
    ! bracket nodes n_int (=-1) and n_int+1 (=0)
    ! Y[-1]:
    call flink(kx0, ky0, dk, -1, link(:, :, 1))
    call gicov_int_chain(link(:, :, 1:1), 2, 1, .true., Wlo)
    call eig2(kx0 + n_int * dk, ky0, wr, Vr)
    Oloc = (0d0, 0d0); Oloc(1, 1) = wr(1); Oloc(2, 2) = wr(2)
    call gicov_int_transport_op(Wlo, Oloc, 2, Ylo)
    ! Y[0] = diag band energies at kappa (identity transport)
    call eig2(kx0, ky0, wr, Vr)
    Yhi = (0d0, 0d0); Yhi(1, 1) = wr(1); Yhi(2, 2) = wr(2)
    call gicov_int_interp(Ylo, Yhi, frac, 2, Y)
    ! exact reference at x = kappa - a
    call horb2(kx0 - a, ky0, Hr)
    Yref = matmul(matmul(conjg(transpose(Vk)), Hr), Vk)
    call report("T2 floor()+interp reproduces exact co-moving H at x=kappa-a", &
              & maxdiff(Y, Yref, 2) < 1d-11, nfail)
    ! WRONG sign x = kappa + a would match H_orb(kappa+a) instead -> mismatch
    call horb2(kx0 + a, ky0, Hr)
    Yrefw = matmul(matmul(conjg(transpose(Vk)), Hr), Vk)
    call report("T2 wrong-sign target (x=kappa+a) is a different operator off-axis", &
              & maxdiff(Yref, Yrefw, 2) > 1d-3, nfail)
    call report("T2 floor bracketing: a=+0.6dk -> n_int=-1, frac=0.4 (floor, not int)", &
              & (n_int == -1) .and. (abs(frac - 0.4d0) < 1d-12), nfail)
  end subroutine t2_interp_floor

  !----------------------------------------------------------------
  subroutine t3_conservation(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 3, nstep = 400
    complex(8) :: H(nb, nb), rho(nb, nb), rho2(nb, nb)
    complex(8) :: P(nb, nb), R(nb, nb), cwork(64)
    real(8)    :: eps(nb), rwork(64), tr0, herr, dt, gamma
    integer    :: a, b, it, ierr, blk(nb)
    ! a fixed generic Hermitian H and a Hermitian trace-1 initial rho
    H = (0d0, 0d0)
    H(1, 1) = 0.20d0; H(2, 2) = 0.55d0; H(3, 3) = 0.90d0
    H(1, 2) = cmplx(0.10d0, 0.03d0, 8); H(2, 1) = conjg(H(1, 2))
    H(1, 3) = cmplx(-0.04d0, 0.07d0, 8); H(3, 1) = conjg(H(1, 3))
    H(2, 3) = cmplx(0.06d0, -0.02d0, 8); H(3, 2) = conjg(H(2, 3))
    rho = (0d0, 0d0)
    rho(1, 1) = 0.6d0; rho(2, 2) = 0.3d0; rho(3, 3) = 0.1d0
    rho(1, 2) = cmplx(0.05d0, 0.02d0, 8); rho(2, 1) = conjg(rho(1, 2))
    rho(1, 3) = cmplx(0.01d0, -0.03d0, 8); rho(3, 1) = conjg(rho(1, 3))
    tr0 = real(rho(1, 1) + rho(2, 2) + rho(3, 3), 8)
    dt = 0.02d0; gamma = 1d0 / 50d0
    do it = 1, nstep
      call gicov_int_step_k(H, rho, nb, dt, gamma, 'step', 2d-3, 0d0, 1d-9, &
                          & eps, P, R, cwork, 64, rwork, blk, rho2, ierr)
      rho = rho2
    end do
    ! trace preserved (populations invariant in H's eigenbasis, exactly)
    call report("T3 trace conserved to machine level over 400 steps", &
              & abs(real(rho(1,1)+rho(2,2)+rho(3,3),8) - tr0) < 1d-12, nfail)
    herr = 0d0
    do b = 1, nb
      do a = 1, nb
        herr = max(herr, abs(rho(a, b) - conjg(rho(b, a))))
      end do
    end do
    call report("T3 Hermiticity preserved to machine level", herr < 1d-12, nfail)
  end subroutine t3_conservation

  !----------------------------------------------------------------
  subroutine t4_t2_covariance(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 2
    complex(8) :: H(nb, nb), rho(nb, nb), rho2(nb, nb)
    complex(8) :: P(nb, nb), R(nb, nb), cwork(64)
    real(8)    :: eps(nb), rwork(64), dt, gamma, c0, cN, expected
    integer    :: it, ierr, blk(nb)
    dt = 0.05d0; gamma = 1d0 / 10d0
    ! (a) EXACTLY degenerate pair: H = 0 (eps_1 = eps_2) -> g(0)=0 -> no decay
    H = (0d0, 0d0)
    rho = (0d0, 0d0)
    rho(1, 1) = 0.5d0; rho(2, 2) = 0.5d0
    rho(1, 2) = cmplx(0.3d0, 0d0, 8); rho(2, 1) = conjg(rho(1, 2))
    c0 = abs(rho(1, 2))
    do it = 1, 100
      call gicov_int_step_k(H, rho, nb, dt, gamma, 'step', 2d-3, 0d0, 1d-9, &
                          & eps, P, R, cwork, 64, rwork, blk, rho2, ierr)
      rho = rho2
    end do
    call report("T4 exact-degenerate pair is NOT dephased (g(0)=0 covariance)", &
              & abs(abs(rho(1, 2)) - c0) < 1d-12, nfail)
    ! (b) gapped pair: |eps_1 - eps_2| = 1.0 >> theta -> coherence decays exp(-gamma*dt*N)
    H = (0d0, 0d0); H(1, 1) = -0.5d0; H(2, 2) = 0.5d0
    rho = (0d0, 0d0)
    rho(1, 1) = 0.5d0; rho(2, 2) = 0.5d0
    rho(1, 2) = cmplx(0.3d0, 0d0, 8); rho(2, 1) = conjg(rho(1, 2))
    c0 = abs(rho(1, 2))
    do it = 1, 100
      call gicov_int_step_k(H, rho, nb, dt, gamma, 'step', 2d-3, 0d0, 1d-9, &
                          & eps, P, R, cwork, 64, rwork, blk, rho2, ierr)
      rho = rho2
    end do
    cN = abs(rho(1, 2))
    expected = c0 * exp(-gamma * dt * 100)
    call report("T4 gapped pair decays by exp(-dt/T2) (gate weight 1)", &
              & abs(cN - expected) < 1d-10, nfail)
  end subroutine t4_t2_covariance

  !----------------------------------------------------------------
  subroutine t5_cache_sizing(nfail)
    integer, intent(inout) :: nfail
    integer(8) :: nb_, want
    ! graphene k99 production: a_max=0.102424, dk=0.00788381 -> j_max = 13
    call report("T5 j_max = 13 for graphene k99 (ceil(a_max/dk))", &
              & gicov_int_jmax(0.102424d0, 0.00788381d0) == 13, nfail)
    ! int64 cache bytes = 80*(2j+1)*nb^2*nk_local, exact past 2**31.  80 = FIVE
    ! complex(8) matrices per shift: H~, v~_{1..3}, and the co-moving f0
    ! reference F~0 (T12) -- NOT 64: under-counting the cache is exactly the
    ! class of error the memory gate exists to block.
    call report("T5 cache bytes formula (j=13, nb=24, nk_local=1000): 5 matrices", &
              & gicov_int_cache_bytes(13, 24, 1000) == &
              & 80_8 * 27_8 * 24_8 * 24_8 * 1000_8, nfail)
    ! production-scale no-overflow (nb=192, nk_local=20000, j=13): > 2**31 bytes
    nb_ = 192_8
    want = 80_8 * 27_8 * nb_ * nb_ * 20000_8
    call report("T5 int64 cache bytes do not overflow at production scale", &
              & (gicov_int_cache_bytes(13, 192, 20000) == want) .and. (want > 2147483647_8), nfail)
  end subroutine t5_cache_sizing

  !----------------------------------------------------------------
  subroutine t6_axis_guard(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nt = 5
    real(8) :: q(3, nt), tol
    logical :: ok
    integer :: axis, it
    tol = 1d-8
    ! single axis (x): accepted, axis=1
    q = 0d0
    do it = 1, nt
      q(1, it) = 0.01d0 * it
    end do
    call gicov_int_axis_single(q, 3, nt, tol, ok, axis)
    call report("T6 single-axis (x) trajectory accepted, axis=1", ok .and. axis == 1, nfail)
    ! [110]: x and y both move -> rejected
    q = 0d0
    do it = 1, nt
      q(1, it) = 0.01d0 * it; q(2, it) = 0.01d0 * it
    end do
    call gicov_int_axis_single(q, 3, nt, tol, ok, axis)
    call report("T6 [110] non-axis linear polarization rejected", (.not. ok) .and. axis == 0, nfail)
    ! circular: x=cos, y=sin -> both nonzero over trajectory -> rejected
    q = 0d0
    do it = 1, nt
      q(1, it) = 0.02d0 * cos(0.5d0 * it); q(2, it) = 0.02d0 * sin(0.5d0 * it)
    end do
    call gicov_int_axis_single(q, 3, nt, tol, ok, axis)
    call report("T6 circular polarization rejected", .not. ok, nfail)
    ! single axis z accepted
    q = 0d0
    do it = 1, nt
      q(3, it) = -0.03d0 * it
    end do
    call gicov_int_axis_single(q, 3, nt, tol, ok, axis)
    call report("T6 single-axis (z) accepted, axis=3", ok .and. axis == 3, nfail)
  end subroutine t6_axis_guard

  !----------------------------------------------------------------
  subroutine t7_floor_shift(nfail)
    integer, intent(inout) :: nfail
    integer :: n_int
    real(8) :: frac
    ! a>0, dk=1: s=-a; a=0.5 -> s=-0.5 -> floor=-1 (int() would give 0 -> WRONG)
    call gicov_int_floor_shift(0.5d0, 1d0, n_int, frac)
    call report("T7 a=+0.5, dk=1 -> n_int=-1, frac=0.5 (floor, not int)", &
              & (n_int == -1) .and. (abs(frac - 0.5d0) < 1d-12), nfail)
    ! a<0: s=+; a=-1.3, dk=1 -> s=1.3 -> floor=1, frac=0.3
    call gicov_int_floor_shift(-1.3d0, 1d0, n_int, frac)
    call report("T7 a=-1.3, dk=1 -> n_int=1, frac=0.3", &
              & (n_int == 1) .and. (abs(frac - 0.3d0) < 1d-12), nfail)
  end subroutine t7_floor_shift

  !----------------------------------------------------------------
  ! Mode-predicate truth table (SDD risk #2: "enum added, but a scattered
  ! trim(sbe_lg_degen)=='gicov' branch was missed").  This is the LOCAL,
  ! statically-decidable half of the legacy byte-regression guard: every
  ! pre-existing deck names one of off/gi/gifix/gicov, and NONE of them may
  ! select the new integral machinery -- so the new code is unreachable for
  ! them by construction, whatever the shared call sites do.  (The full
  ! output byte-compare of a legacy 'gicov' + a default deck is the Fugaku
  ! gate; it cannot run in this standalone unit.)
  subroutine t8_mode_predicates(nfail)
    use sbe_lg_mode_ssbe, only: uses_prod_dk, uses_xfull_links, &
                              & uses_fd_gicov, uses_integral_gicov
    integer, intent(inout) :: nfail
    logical :: ok

    ! full truth table over every accepted sbe_lg_degen value
    ok = (.not. uses_prod_dk('off'))       .and. (.not. uses_xfull_links('off')) &
   .and. (.not. uses_fd_gicov('off'))      .and. (.not. uses_integral_gicov('off'))
    call report("T8 'off': selects no covariant machinery", ok, nfail)

    ok = uses_prod_dk('gi')                .and. (.not. uses_xfull_links('gi')) &
   .and. (.not. uses_fd_gicov('gi'))       .and. (.not. uses_integral_gicov('gi')) &
   .and. uses_prod_dk('gifix')             .and. (.not. uses_xfull_links('gifix')) &
   .and. (.not. uses_fd_gicov('gifix'))    .and. (.not. uses_integral_gicov('gifix'))
    call report("T8 'gi'/'gifix': prod_dk only (no X-full, no FD-gicov, no integral)", ok, nfail)

    ok = uses_prod_dk('gicov')             .and. uses_xfull_links('gicov') &
   .and. uses_fd_gicov('gicov')            .and. (.not. uses_integral_gicov('gicov'))
    call report("T8 'gicov': X-full + FD covariant, NEVER the integral path", ok, nfail)

    ok = uses_prod_dk('gicov_int')         .and. uses_xfull_links('gicov_int') &
   .and. (.not. uses_fd_gicov('gicov_int')) .and. uses_integral_gicov('gicov_int')
    call report("T8 'gicov_int': X-full + integral, NEVER the FD-gicov path", ok, nfail)

    ! the decisive legacy-regression assertion: no pre-existing mode reaches
    ! the integral propagator / the moving-gap T2 gate / the transport cache
    ok = (.not. uses_integral_gicov('off')) .and. (.not. uses_integral_gicov('gi')) &
   .and. (.not. uses_integral_gicov('gifix')) .and. (.not. uses_integral_gicov('gicov'))
    call report("T8 NO legacy mode (off/gi/gifix/gicov) enters the integral path", ok, nfail)
  end subroutine t8_mode_predicates

  !----------------------------------------------------------------
  ! Eigensolver-failure contract.  A corrupted instantaneous H~ must be
  ! REPORTED (ierr /= 0) so the caller can abort the whole communicator; it may
  ! never be silently absorbed, because both former fallbacks (freeze the block
  ! / read diag(rho~)) preserve the trace and therefore slip past the Ne monitor.
  subroutine t9_eigensolver_status(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 3
    complex(8) :: H(nb, nb), rho(nb, nb), rho2(nb, nb)
    complex(8) :: P(nb, nb), R(nb, nb), cwork(64)
    real(8)    :: eps(nb), rwork(64), nocc(nb), nan
    integer    :: ierr, ierr_ok, a, blk(nb)
    logical    :: ok

    rho = (0d0, 0d0)
    rho(1, 1) = 0.6d0; rho(2, 2) = 0.3d0; rho(3, 3) = 0.1d0
    rho(1, 2) = cmplx(0.20d0, 0.05d0, 8); rho(2, 1) = conjg(rho(1, 2))

    ! (a) clean, GAPPED, NON-diagonal H -> success, and the readout is the
    !     INSTANTANEOUS-eigenbasis projection, provably not diag(rho~):
    !     if the forbidden diag(rho~) fallback were ever reached it would
    !     return 0.6/0.3/0.1 -- so a mismatch here is the detector.
    H = (0d0, 0d0)
    H(1, 1) = 0.20d0; H(2, 2) = 0.55d0; H(3, 3) = 0.90d0
    H(1, 2) = cmplx(0.18d0, 0.04d0, 8); H(2, 1) = conjg(H(1, 2))
    H(1, 3) = cmplx(-0.09d0, 0.06d0, 8); H(3, 1) = conjg(H(1, 3))
    call gicov_int_occupation_k(H, rho, nb, 1d-9, eps, P, R, cwork, 64, rwork, nocc, blk, ierr_ok)
    ok = (ierr_ok == 0)
    ! the instantaneous populations must DIFFER from the frozen-basis diagonal
    ok = ok .and. (maxval(abs(nocc - (/ 0.6d0, 0.3d0, 0.1d0 /))) > 1d-3)
    ! ... while still summing to the (basis-invariant) trace
    ok = ok .and. (abs(sum(nocc) - 1.0d0) < 1d-12)
    call report("T9 occupation: ierr=0, instantaneous projection (NOT diag(rho~)), trace kept", &
              & ok, nfail)

    call gicov_int_step_k(H, rho, nb, 0.02d0, 1d0/50d0, 'step', 2d-3, 0d0, 1d-9, &
                        & eps, P, R, cwork, 64, rwork, blk, rho2, ierr_ok)
    ok = (ierr_ok == 0)
    ! and it actually EVOLVED (a frozen block would return rho unchanged)
    ok = ok .and. (maxdiff(rho2, rho, nb) > 1d-6)
    call report("T9 step: ierr=0 on a good H and the block genuinely evolves (no freeze)", &
              & ok, nfail)

    ! (b) non-finite H~ -> MUST be reported by BOTH kernels.  zheev returns
    !     info=0 here (Accelerate), so this exercises the finite-eigenvalue
    !     half of the status; without it the corruption would propagate.
    nan = 0d0; nan = nan / nan
    H = (0d0, 0d0)
    H(1, 1) = 0.20d0; H(2, 2) = 0.55d0
    H(3, 3) = cmplx(nan, 0d0, 8)
    call gicov_int_step_k(H, rho, nb, 0.02d0, 1d0/50d0, 'step', 2d-3, 0d0, 1d-9, &
                        & eps, P, R, cwork, 64, rwork, blk, rho2, ierr)
    call report("T9 step: non-finite H~ is REPORTED (ierr/=0), not silently frozen", &
              & ierr /= 0, nfail)
    call gicov_int_occupation_k(H, rho, nb, 1d-9, eps, P, R, cwork, 64, rwork, nocc, blk, ierr)
    ok = (ierr /= 0)
    ! and it must NOT have quietly returned the forbidden diag(rho~) fallback
    do a = 1, nb
      if (abs(nocc(a) - real(rho(a, a), 8)) > 1d-14) cycle
    end do
    ok = ok .and. (maxval(abs(nocc - (/ 0.6d0, 0.3d0, 0.1d0 /))) > 1d-14)
    call report("T9 occupation: non-finite H~ REPORTED, no diag(rho~) fallback value", &
              & ok, nfail)
  end subroutine t9_eigensolver_status

  !----------------------------------------------------------------
  ! Bracketing at the edge of the cached transport span.  gicov_int_bracket
  ! must never name a node outside [-jmax, jmax]: when the evaluation point
  ! falls exactly ON a node (frac = 0) the upper bracket must be the node
  ! itself, because the cache has nothing beyond the span and the upper
  ! endpoint would carry zero weight anyway.
  subroutine t10_bracket_edge(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: jmax = 13
    integer :: n_lo, n_hi, ierr
    real(8) :: frac
    complex(8) :: Ylo(2, 2), Yhi(2, 2), Y(2, 2)
    logical :: ok

    ! q = -jmax  ->  s = +jmax  ->  lands exactly on node +jmax.
    ! OLD behaviour: n_int=jmax, n_int+1 = jmax+1 -> outside the cache -> abort.
    call gicov_int_bracket(-dble(jmax), 1d0, jmax, n_lo, n_hi, frac, ierr)
    ok = (ierr == 0) .and. (n_lo == jmax) .and. (n_hi == jmax) &
   .and. (abs(frac) < 1d-15)
    call report("T10 q=-j_max: bracket collapses onto node +j_max (no j_max+1, no abort)", &
              & ok, nfail)

    ! the mirror endpoint q = +jmax -> s = -jmax -> node -jmax
    call gicov_int_bracket(dble(jmax), 1d0, jmax, n_lo, n_hi, frac, ierr)
    ok = (ierr == 0) .and. (n_lo == -jmax) .and. (n_hi == -jmax) &
   .and. (abs(frac) < 1d-15)
    call report("T10 q=+j_max: bracket collapses onto node -j_max", ok, nfail)

    ! ordinary interior point still brackets two DISTINCT nodes
    call gicov_int_bracket(-2.25d0, 1d0, jmax, n_lo, n_hi, frac, ierr)
    ok = (ierr == 0) .and. (n_lo == 2) .and. (n_hi == 3) &
   .and. (abs(frac - 0.25d0) < 1d-12)
    call report("T10 interior q=-2.25: n_lo=2, n_hi=3, frac=0.25 (normal bracket)", ok, nfail)

    ! genuinely outside the cached span -> reported, not silently clamped
    call gicov_int_bracket(-dble(jmax) - 0.5d0, 1d0, jmax, n_lo, n_hi, frac, ierr)
    call report("T10 q beyond the cached span is REPORTED (ierr/=0)", ierr /= 0, nfail)

    ! frac=0 interpolation must return the lower node EXACTLY (so collapsing the
    ! bracket is not merely safe, it is numerically identical to the old intent)
    Ylo = (0d0, 0d0);  Ylo(1, 1) = (0.7d0, 0d0);  Ylo(2, 2) = (0.2d0, 0d0)
    Yhi = (0d0, 0d0);  Yhi(1, 1) = (9.9d0, 0d0);  Yhi(2, 2) = (9.9d0, 0d0)
    call gicov_int_interp(Ylo, Ylo, 0d0, 2, Y)
    call report("T10 frac=0 interp on the collapsed bracket reproduces the node exactly", &
              & maxdiff(Y, Ylo, 2) < 1d-15, nfail)
  end subroutine t10_bracket_edge

  !----------------------------------------------------------------
  subroutine t11_degenerate_closure(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 3
    complex(8) :: H(nb, nb), rho(nb, nb), rho2(nb, nb)
    complex(8) :: P(nb, nb), R(nb, nb), cwork(64)
    real(8)    :: eps(nb), rwork(64), nocc(nb), dt, gamma, c0, cN, expected, fl
    integer    :: blk(nb), it, ierr
    logical    :: ok

    fl = 1d-3                 ! degeneracy floor
    dt = 0.05d0; gamma = 1d0 / 10d0

    ! ---- the block builder itself: a CHAIN 0, 0.75f, 1.5f is ONE block
    !      (each adjacent gap 0.75f <= floor) even though the end pair is
    !      1.5f > floor; 3.0f later is a genuinely separate block.
    eps = (/ 0d0, 0.75d0 * fl, 1.5d0 * fl /)
    call gicov_int_degen_blocks(eps, nb, fl, blk)
    call report("T11 connected closure: chain (0, .75f, 1.5f) is ONE block", &
              & (blk(1) == blk(2)) .and. (blk(2) == blk(3)), nfail)
    eps = (/ 0d0, 0.75d0 * fl, 3.0d0 * fl /)
    call gicov_int_degen_blocks(eps, nb, fl, blk)
    call report("T11 a gap > floor STARTS a new block (0, .75f | 3f)", &
              & (blk(1) == blk(2)) .and. (blk(3) /= blk(2)), nfail)

    ! ---- T2 gate on the chain.  theta=0 so the 'step' gate weight is 1 for
    !      ANY gap above the floor: the pairwise rule would then dephase the
    !      end pair (1,3) (gap 1.5f > floor), the closure must not.
    H = (0d0, 0d0)
    H(1, 1) = cmplx(0d0, 0d0, 8)
    H(2, 2) = cmplx(0.75d0 * fl, 0d0, 8)
    H(3, 3) = cmplx(1.5d0 * fl, 0d0, 8)
    rho = (0d0, 0d0)
    rho(1, 1) = 0.4d0; rho(2, 2) = 0.35d0; rho(3, 3) = 0.25d0
    rho(1, 3) = cmplx(0.3d0, 0d0, 8); rho(3, 1) = conjg(rho(1, 3))
    c0 = abs(rho(1, 3))
    do it = 1, 100
      call gicov_int_step_k(H, rho, nb, dt, gamma, 'step', 0d0, 0d0, fl, &
                          & eps, P, R, cwork, 64, rwork, blk, rho2, ierr)
      rho = rho2
    end do
    call report("T11 END pair of a near-degenerate CHAIN is NOT dephased (block closure)", &
              & (ierr == 0) .and. (abs(abs(rho(1, 3)) - c0) < 1d-12), nfail)

    ! ---- control: a genuinely separated pair still dephases (the gate has teeth)
    H = (0d0, 0d0)
    H(1, 1) = cmplx(0d0, 0d0, 8)
    H(2, 2) = cmplx(0.75d0 * fl, 0d0, 8)
    H(3, 3) = cmplx(50d0 * fl, 0d0, 8)          ! far above the floor -> own block
    rho = (0d0, 0d0)
    rho(1, 1) = 0.4d0; rho(2, 2) = 0.35d0; rho(3, 3) = 0.25d0
    rho(1, 3) = cmplx(0.3d0, 0d0, 8); rho(3, 1) = conjg(rho(1, 3))
    c0 = abs(rho(1, 3))
    do it = 1, 100
      call gicov_int_step_k(H, rho, nb, dt, gamma, 'step', 0d0, 0d0, fl, &
                          & eps, P, R, cwork, 64, rwork, blk, rho2, ierr)
      rho = rho2
    end do
    cN = abs(rho(1, 3))
    expected = c0 * exp(-gamma * dt * 100)
    call report("T11 control: a well-separated pair DOES dephase (gate still has teeth)", &
              & abs(cN - expected) < 1d-10, nfail)

    ! ---- basis-invariant readout: members of one degenerate block are not
    !      individually observable, so the block total is reported shared
    !      equally -- and the total is preserved.
    H = (0d0, 0d0)      ! exactly degenerate 1,2 ; 3 separate
    H(3, 3) = cmplx(50d0 * fl, 0d0, 8)
    rho = (0d0, 0d0)
    rho(1, 1) = 0.5d0; rho(2, 2) = 0.1d0; rho(3, 3) = 0.4d0
    call gicov_int_occupation_k(H, rho, nb, fl, eps, P, R, cwork, 64, rwork, &
                              & nocc, blk, ierr)
    ok = (ierr == 0)
    ok = ok .and. (abs(nocc(1) - nocc(2)) < 1d-12)        ! degenerate pair: equal
    ok = ok .and. (abs(nocc(1) + nocc(2) - 0.6d0) < 1d-12) ! block TOTAL invariant
    ok = ok .and. (abs(nocc(3) - 0.4d0) < 1d-12)          ! own block untouched
    ok = ok .and. (abs(sum(nocc) - 1.0d0) < 1d-12)        ! sum rule kept
    call report("T11 occupation is a degenerate-BLOCK sum (invariant, sum rule kept)", &
              & ok, nfail)
  end subroutine t11_degenerate_closure

  !----------------------------------------------------------------
  ! The co-moving f0(x) reference.  Set up a k-point whose transport W genuinely
  ! MIXES the bands (a rotation), a METAL-like fractional GS occupation, and an
  ! UNEXCITED co-moving density rho~ = F~0(x).  The excitation readout must then
  ! be exactly zero in every band.
  subroutine t12_comoving_f0(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 2
    complex(8) :: W(nb, nb), H(nb, nb), F0(nb, nb), rho(nb, nb), D(nb, nb)
    complex(8) :: Od(nb, nb), P(nb, nb), R(nb, nb), cwork(64)
    real(8)    :: eps(nb), rwork(64), nocc(nb), dn(nb), th, fl
    integer    :: blk(nb), ierr
    logical    :: ok

    fl = 1d-9
    th = 0.6d0                                   ! a transport that really mixes
    W = (0d0, 0d0)
    W(1, 1) = cmplx( cos(th), 0d0, 8);  W(1, 2) = cmplx(-sin(th), 0d0, 8)
    W(2, 1) = cmplx( sin(th), 0d0, 8);  W(2, 2) = cmplx( cos(th), 0d0, 8)

    ! co-moving Hamiltonian at the shifted point: H~ = W diag(eigen) W^dag
    Od = (0d0, 0d0);  Od(1, 1) = (-0.5d0, 0d0);  Od(2, 2) = (0.5d0, 0d0)
    call gicov_int_transport_op(W, Od, nb, H)

    ! METAL-like GS occupation (fractional, NOT a rigid 1/0 valence block),
    ! transported into the same frame: F~0 = W diag(f0) W^dag
    Od = (0d0, 0d0);  Od(1, 1) = (0.7d0, 0d0);  Od(2, 2) = (0.3d0, 0d0)
    call gicov_int_transport_op(W, Od, nb, F0)

    ! an UNEXCITED co-moving density is exactly that reference
    rho = F0

    ! (a) excitation against the TRANSPORTED reference: identically zero
    D = rho - F0
    call gicov_int_occupation_k(H, D, nb, fl, eps, P, R, cwork, 64, rwork, &
                              & dn, blk, ierr)
    call report("T12 unexcited co-moving state -> ZERO excitation (f0 subtracted)", &
              & (ierr == 0) .and. (maxval(abs(dn)) < 1d-12), nfail)

    ! (b) what the ABSOLUTE readout (the bug) reported instead: the occupation
    !     itself -- nonzero, and indistinguishable from a real excitation
    call gicov_int_occupation_k(H, rho, nb, fl, eps, P, R, cwork, 64, rwork, &
                              & nocc, blk, ierr)
    call report("T12 the old ABSOLUTE readout is NONzero for the same unexcited state", &
              & (ierr == 0) .and. (maxval(abs(nocc)) > 0.2d0), nfail)

    ! (c) a RIGID '1 below nvb' reference (no metal f0) manufactures a spurious
    !     excitation: this is why the GS occupation must be transported, not
    !     assumed.  diag(1,0) in the co-moving frame, subtracted from rho~.
    Od = (0d0, 0d0);  Od(1, 1) = (1d0, 0d0);  Od(2, 2) = (0d0, 0d0)
    D = rho - Od
    call gicov_int_occupation_k(H, D, nb, fl, eps, P, R, cwork, 64, rwork, &
                              & dn, blk, ierr)
    call report("T12 a RIGID 1/0 reference fakes excitation for a metal (so F~0 is needed)", &
              & (ierr == 0) .and. (maxval(abs(dn)) > 1d-2), nfail)

    ! (d) a GENUINE excitation is still reported, with the right sign and size:
    !     move 0.1 electrons from the lower to the upper instantaneous state.
    Od = (0d0, 0d0);  Od(1, 1) = (-0.1d0, 0d0);  Od(2, 2) = (0.1d0, 0d0)
    call gicov_int_transport_op(W, Od, nb, D)   ! same frame as H~/F~0
    rho = F0 + D
    D = rho - F0
    call gicov_int_occupation_k(H, D, nb, fl, eps, P, R, cwork, 64, rwork, &
                              & dn, blk, ierr)
    ok = (ierr == 0) .and. (abs(dn(1) + 0.1d0) < 1d-12) .and. (abs(dn(2) - 0.1d0) < 1d-12)
    call report("T12 a genuine excitation IS reported (-0.1 / +0.1, charge-conserving)", &
              & ok, nfail)
  end subroutine t12_comoving_f0

end program test_gicov_integral
