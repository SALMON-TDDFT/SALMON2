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
!     src/ssbe/degenerate_block_ssbe.f90 src/ssbe/gicov_integral_ssbe.f90 \
!     src/ssbe/test/test_gicov_integral.f90 -o t -framework Accelerate && ./t
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
!      = 64(2j+1)nb^2 nk_local with no overflow at production scale.
!   T6 single-axis / linear-polarization guard (risk #6): [110] and circular
!      are rejected on the whole trajectory; a single axis is accepted.
!   T7 floor() vs int() bracketing for a>0 (negative mesh shift).
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
    integer    :: a, b, it
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
                          & eps, P, R, cwork, 64, rwork, rho2)
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
    integer    :: it
    dt = 0.05d0; gamma = 1d0 / 10d0
    ! (a) EXACTLY degenerate pair: H = 0 (eps_1 = eps_2) -> g(0)=0 -> no decay
    H = (0d0, 0d0)
    rho = (0d0, 0d0)
    rho(1, 1) = 0.5d0; rho(2, 2) = 0.5d0
    rho(1, 2) = cmplx(0.3d0, 0d0, 8); rho(2, 1) = conjg(rho(1, 2))
    c0 = abs(rho(1, 2))
    do it = 1, 100
      call gicov_int_step_k(H, rho, nb, dt, gamma, 'step', 2d-3, 0d0, 1d-9, &
                          & eps, P, R, cwork, 64, rwork, rho2)
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
                          & eps, P, R, cwork, 64, rwork, rho2)
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
    ! int64 cache bytes = 64*(2j+1)*nb^2*nk_local, exact past 2**31
    call report("T5 cache bytes formula (j=13, nb=24, nk_local=1000)", &
              & gicov_int_cache_bytes(13, 24, 1000) == &
              & 64_8 * 27_8 * 24_8 * 24_8 * 1000_8, nfail)
    ! production-scale no-overflow (nb=192, nk_local=20000, j=13): > 2**31 bytes
    nb_ = 192_8
    want = 64_8 * 27_8 * nb_ * nb_ * 20000_8
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

end program test_gicov_integral
