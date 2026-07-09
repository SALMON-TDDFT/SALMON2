! src/ssbe/test/test_gisbe_covariance.f90
! LG-SBE Tier2  Pc = U(N) covariance gate for the non-Abelian Berry connection
! xi (degenerate_block_ssbe: xi_block_from_overlap / build_xi, Pb3) and for the
! Pb4 gauge-invariant (GI) current.
!
! WHAT THIS PROVES
!   Inside an *exactly* degenerate block (Delta_omega = 0) the eigenvectors carry
!   a U(N) gauge freedom: at each k they may be rotated |u_a> -> sum_b |u_b> U_ba(k),
!   i.e. u -> u.U(k).  A correct length-gauge coupling must therefore be
!     (i)  gauge COVARIANT for the xi-derived RHS commutator term  (E.d-type term),
!          RHS~  =  U^H (RHS) U   (transforms as a rank-2 U(N) tensor), and
!     (ii) gauge INVARIANT for the physical current  j = sum rho(jb,ib) p(ib,jb).
!   This file asserts both, and shows the *fix* (xi) passes the gate while the
!   naive dipole d = i p / Delta_omega (which is 0/0 in a degenerate block) fails.
!
! WHY THE GAUGE IS HELD CONSTANT ON THE FINITE-DIFFERENCE STENCIL (RHS test)
!   The discrete connection is  xi = s i logm(polar(M)) / dk  with the link overlap
!   M(a,b) = <u_a(k)|u_b(k+b)>.  Under u->u.U it becomes M~ = U(k)^H M U(k+b), so
!     xi~ = s i logm( U(k)^H polar(M) U(k+b) ) / dk .
!   When U is CONSTANT across the stencil (U(k)=U(k+b)=W) this collapses to the
!   exact tensor law  xi~ = W^H xi W  (similarity commutes with logm), hence
!   [xi~,rho~] = W^H [xi,rho] W.  When U VARIES across the stencil the extra piece
!   is precisely the pure-gauge inhomogeneous term  i U^H d_k U  --- a Berry-connection
!   transformation --- which is cancelled ONLY by the intraband k-gradient (grad)
!   term  E.d_k rho.  That grad term is deliberately OUT OF SCOPE here (it is also
!   the non-covariant one under block switching / Pe), so to isolate and test the
!   connection commutator alone we probe it with a stencil-constant gauge.  The
!   grad caveat is itself asserted (test_grad_caveat) so the boundary is explicit.
!   The current test (which never references xi) uses a FULLY k-dependent gauge.
!
! The kernel calls LAPACK (zheev/zgeev).  Standalone build from src/ssbe:
!     gfortran degenerate_block_ssbe.f90 test/test_gisbe_covariance.f90 \
!              -o t -llapack -lblas && ./t
! (CMake add_test is not the ssbe convention; run manually as above.)

program test_gisbe_covariance
  use degenerate_block_ssbe, only: xi_block_from_overlap, build_xi, xi_sign
  implicit none
  complex(8), parameter :: zi = (0d0, 1d0)
  integer :: nfail
  nfail = 0

  call test_current_invariance(nfail)   ! Pb4 GI current invariant, fully k-dependent U(2)
  call test_rhs_covariance(2, nfail)     ! [xi,rho] covariant, U(2) 2-fold-degenerate block
  call test_rhs_covariance(3, nfail)     ! [xi,rho] covariant, U(3) 3-fold-degenerate block
  call test_grad_caveat(nfail)           ! varying-stencil gauge: [xi,rho] is NOT a tensor (grad term)
  call test_negative_naive_dipole(nfail) ! d=i p/Delta_omega collapses -> gate FAILS; xi PASSES
  call test_build_xi_covariance(nfail)   ! end-to-end build_xi on a multi-k crossing spectrum

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_gisbe_covariance checks passed."
  end if

contains

  !======================= assert helpers (ssbe style) =======================
  subroutine check_true(cond, label, nfail)
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

  subroutine check_close_r(got, want, tol, label, nfail)
    real(8), intent(in) :: got, want, tol
    character(*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (abs(got - want) <= tol) then
      write(*, '(a)') "ok    " // label
    else
      write(*, '(a,es12.4)') "FAIL  " // label // "  |got-want|=", abs(got - want)
      nfail = nfail + 1
    end if
  end subroutine check_close_r

  subroutine check_int(got, want, label, nfail)
    integer, intent(in) :: got, want
    character(*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (got == want) then
      write(*, '(a)') "ok    " // label
    else
      write(*, '(a,2(1x,i0))') "FAIL  " // label // "  got/want=", got, want
      nfail = nfail + 1
    end if
  end subroutine check_int

  !======================= small matrix utilities ============================
  function hc(A) result(B)                       ! Hermitian conjugate A^H
    complex(8), intent(in) :: A(:, :)
    complex(8), allocatable :: B(:, :)
    B = conjg(transpose(A))
  end function hc

  function comm(A, B) result(C)                  ! commutator [A,B]
    complex(8), intent(in) :: A(:, :), B(:, :)
    complex(8), allocatable :: C(:, :)
    C = matmul(A, B) - matmul(B, A)
  end function comm

  function rot(W, A) result(R)                   ! covariant rotation W^H A W
    complex(8), intent(in) :: W(:, :), A(:, :)
    complex(8), allocatable :: R(:, :)
    R = matmul(matmul(conjg(transpose(W)), A), W)
  end function rot

  real(8) function mdiff(A, B)                    ! max_ij |A-B|
    complex(8), intent(in) :: A(:, :), B(:, :)
    mdiff = maxval(abs(A - B))
  end function mdiff

  real(8) function mmax(A)                        ! max_ij |A|
    complex(8), intent(in) :: A(:, :)
    mmax = maxval(abs(A))
  end function mmax

  ! Pb4 GI current contraction j = sum_{ib,jb} rho(jb,ib) p(ib,jb) = Tr(rho p).
  complex(8) function trace_rhop(rho, p) result(tr)
    complex(8), intent(in) :: rho(:, :), p(:, :)
    integer :: ib, jb, n
    n = size(rho, 1)
    tr = (0d0, 0d0)
    do jb = 1, n
      do ib = 1, n
        tr = tr + rho(jb, ib) * p(ib, jb)
      end do
    end do
  end function trace_rhop

  ! Fixed, reproducible, nontrivial complex-Hermitian matrix (shift varies entries).
  function herm_fixed(n, shift) result(H)
    integer, intent(in) :: n
    real(8), intent(in) :: shift
    complex(8), allocatable :: H(:, :)
    integer :: a, b
    allocate(H(n, n))
    do a = 1, n
      H(a, a) = dcmplx(1d0 + 0.3d0 * a + shift, 0d0)
      do b = a + 1, n
        H(a, b) = dcmplx(0.2d0 * a - 0.1d0 * b + shift, 0.15d0 * a + 0.05d0 * b)
        H(b, a) = conjg(H(a, b))
      end do
    end do
  end function herm_fixed

  ! n x n complex Givens (unitary) rotation in the (i,j) plane by angle th, phase ph.
  function givens_c(n, i, j, th, ph) result(G)
    integer, intent(in) :: n, i, j
    real(8), intent(in) :: th, ph
    complex(8), allocatable :: G(:, :)
    integer :: a
    real(8) :: c, s
    allocate(G(n, n)); G = (0d0, 0d0)
    do a = 1, n
      G(a, a) = (1d0, 0d0)
    end do
    c = cos(th); s = sin(th)
    G(i, i) = dcmplx(c, 0d0); G(j, j) = dcmplx(c, 0d0)
    G(i, j) = exp(zi * ph) * s
    G(j, i) = -exp(-zi * ph) * s
  end function givens_c

  ! General U(2) gauge  e^{i chi} [[e^{i al}c, e^{i be}s],[-e^{-i be}s, e^{-i al}c]].
  function u2_gauge(al, be, th, chi) result(U)
    real(8), intent(in) :: al, be, th, chi
    complex(8), allocatable :: U(:, :)
    real(8) :: c, s
    complex(8) :: g
    allocate(U(2, 2))
    c = cos(th); s = sin(th); g = exp(zi * chi)
    U(1, 1) =  g * exp( zi * al) * c
    U(1, 2) =  g * exp( zi * be) * s
    U(2, 1) = -g * exp(-zi * be) * s
    U(2, 2) =  g * exp(-zi * al) * c
  end function u2_gauge

  ! Nontrivial U(3) gauge as a product of three complex Givens rotations.
  function u3_gauge() result(U)
    complex(8), allocatable :: U(:, :)
    U = matmul(matmul(givens_c(3, 1, 2, 0.6d0,  0.4d0), &
                      givens_c(3, 1, 3, 0.5d0, -0.3d0)), &
                      givens_c(3, 2, 3, 0.7d0,  0.2d0))
  end function u3_gauge

  ! Smooth d x d Bloch-state frame u(ang) (columns) of an EXACTLY degenerate block;
  ! its k-derivative is nonzero so the link overlap yields a nonzero connection.
  function block_states(d, ang) result(Q)
    integer, intent(in) :: d
    real(8), intent(in) :: ang
    complex(8), allocatable :: Q(:, :)
    if (d == 2) then
      Q = givens_c(2, 1, 2, 0.5d0 * ang, 0d0)             ! real smooth rotation
    else
      Q = matmul(matmul(givens_c(3, 1, 2, 0.50d0 * ang,  0.20d0), &
                        givens_c(3, 1, 3, 0.30d0 * ang, -0.10d0)), &
                        givens_c(3, 2, 3, 0.40d0 * ang,  0.15d0))
    end if
  end function block_states

  ! Embed a d x d block gauge W at bands i0..i0+d-1 of an nb x nb identity
  ! (gauge acts inside the degenerate block, identity on all other bands).
  function embed(nb, i0, W) result(V)
    integer, intent(in) :: nb, i0
    complex(8), intent(in) :: W(:, :)
    complex(8), allocatable :: V(:, :)
    integer :: a, b, d
    d = size(W, 1)
    allocate(V(nb, nb)); V = (0d0, 0d0)
    do a = 1, nb
      V(a, a) = (1d0, 0d0)
    end do
    do b = 1, d
      do a = 1, d
        V(i0 + a - 1, i0 + b - 1) = W(a, b)
      end do
    end do
  end function embed

  ! nb x nb frame for the build_xi model: identity on bands 1 and nb, a real
  ! 2x2 rotation R(ang) on the exactly degenerate block {2,3}.
  function frame4(nb, ang) result(Q)
    integer, intent(in) :: nb
    real(8), intent(in) :: ang
    complex(8), allocatable :: Q(:, :)
    integer :: a
    real(8) :: c, s
    allocate(Q(nb, nb)); Q = (0d0, 0d0)
    do a = 1, nb
      Q(a, a) = (1d0, 0d0)
    end do
    c = cos(ang); s = sin(ang)
    Q(2, 2) = dcmplx(c, 0d0); Q(2, 3) = dcmplx(-s, 0d0)
    Q(3, 2) = dcmplx(s, 0d0); Q(3, 3) = dcmplx(c, 0d0)
  end function frame4

  !================================ tests ====================================

  ! (1) Pb4 GI current  j = sum rho(jb,ib) p(ib,jb)  is gauge INVARIANT under a
  !     FULLY k-dependent U(2) gauge on the degenerate block (a genuinely local
  !     gauge: a different U(2) at each k).  The current never references xi, so
  !     invariance needs no stencil constraint -- it is a pure tensor-trace fact:
  !       Tr( (V^H rho V)(V^H p V) ) = Tr(rho p).
  subroutine test_current_invariance(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 4, nk = 3, i0 = 2       ! block {2,3}
    complex(8) :: rho(nb, nb), p(nb, nb)
    complex(8), allocatable :: Wk(:, :), V(:, :)
    complex(8) :: j0, jt
    integer :: ik
    real(8) :: al, be, th, chi
    j0 = (0d0, 0d0); jt = (0d0, 0d0)
    do ik = 1, nk
      rho = herm_fixed(nb, 0.10d0 * ik)
      p   = 1.30d0 * herm_fixed(nb, 0.20d0 * ik + 0.30d0)   ! Hermitian velocity proxy
      j0  = j0 + trace_rhop(rho, p)
      ! k-DEPENDENT U(2) gauge inside the block, identity elsewhere:
      al = 0.30d0 * ik; be = -0.20d0 * ik; th = 0.40d0 + 0.15d0 * ik; chi = 0.10d0 * ik
      Wk = u2_gauge(al, be, th, chi)
      V  = embed(nb, i0, Wk)
      jt = jt + trace_rhop(rot(V, rho), rot(V, p))
    end do
    call check_close_r(abs(aimag(j0)), 0d0, 1d-12, "current: Tr(rho p) is real", nfail)
    call check_close_r(abs(jt - j0),   0d0, 1d-12, &
                       "Pb4 GI current invariant under k-dependent U(2) gauge", nfail)
    call check_true(abs(dble(j0)) > 1d-6, "current: reference current is nonzero", nfail)
  end subroutine test_current_invariance

  ! (2/3) The xi-derived RHS commutator term is gauge COVARIANT for a d-fold
  !       exactly-degenerate block under U(d) (d=2 and d=3).  The physical term is
  !       -i (E.eps) [xi, rho]; the scalar prefactor is irrelevant to the tensor
  !       law, so we test [xi, rho] directly.  Gauge held constant on the stencil
  !       (see header) so the connection is isolated from the grad term.
  subroutine test_rhs_covariance(d, nfail)
    integer, intent(in) :: d
    integer, intent(inout) :: nfail
    complex(8), allocatable :: W(:, :), Qk(:, :), Qkb(:, :), M(:, :), Mt(:, :)
    complex(8), allocatable :: rho(:, :), p(:, :)
    complex(8) :: xi(d, d), xit(d, d)
    real(8) :: dk, ang0, resu
    integer :: info, infot
    character(1) :: tag
    ang0 = 0.4d0; dk = 0.1d0
    if (d == 2) then
      W = u2_gauge(0.5d0, -0.4d0, 0.8d0, 0.3d0); tag = '2'
    else
      W = u3_gauge();                            tag = '3'
    end if
    ! link overlap M(a,b) = <u_a(k)|u_b(k+dk)> from the smooth degenerate frame
    Qk  = block_states(d, ang0)
    Qkb = block_states(d, ang0 + dk)
    M   = matmul(hc(Qk), Qkb)
    ! stencil-constant gauge: u -> u.W at BOTH k and k+dk  =>  M~ = W^H M W
    Mt  = rot(W, M)
    call xi_block_from_overlap(M,  dk, xi_sign, xi,  info,  resu)
    call xi_block_from_overlap(Mt, dk, xi_sign, xit, infot)
    call check_true(info == 0 .and. infot == 0, "U("//tag//"): xi kernel ok (info==0)", nfail)
    call check_close_r(mdiff(xi, hc(xi)), 0d0, 1d-10, "U("//tag//"): xi Hermitian", nfail)
    call check_true(mmax(xi) > 1d-6, "U("//tag//"): xi nonzero (real intra-block coupling)", nfail)
    ! underlying identity: xi is a U(N) tensor
    call check_close_r(mdiff(xit, rot(W, xi)), 0d0, 1d-10, &
                       "U("//tag//"): xi is a gauge tensor  xi~ = W^H xi W", nfail)
    ! RHS commutator covariance and current invariance in the same block
    rho = herm_fixed(d, 0.0d0)
    p   = 1.70d0 * herm_fixed(d, 0.5d0)
    call check_close_r(mdiff(comm(xit, rot(W, rho)), rot(W, comm(xi, rho))), 0d0, 1d-10, &
                       "U("//tag//"): RHS [xi,rho] covariant  RHS~ = W^H RHS W", nfail)
    call check_true(mmax(comm(xi, rho)) > 1d-6, "U("//tag//"): RHS coupling nonzero", nfail)
    call check_close_r(abs(trace_rhop(rot(W, rho), rot(W, p)) - trace_rhop(rho, p)), &
                       0d0, 1d-12, "U("//tag//"): block current Tr(rho p) invariant", nfail)
  end subroutine test_rhs_covariance

  ! (4) Boundary / caveat (documented, asserted): with a gauge that VARIES across
  !     the finite-difference stencil (U(k) /= U(k+dk)) the bare connection
  !     commutator [xi,rho] is NOT a tensor -- the discrepancy is the pure-gauge
  !     inhomogeneous term i U^H d_k U that the OUT-OF-SCOPE grad term cancels.
  !     We assert the discrepancy is O(1) large, i.e. covariance genuinely
  !     requires the stencil-constant probe used in test_rhs_covariance.
  subroutine test_grad_caveat(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: d = 2
    complex(8), allocatable :: Wk(:, :), Wkb(:, :), Qk(:, :), Qkb(:, :), M(:, :), Mt(:, :)
    complex(8), allocatable :: rho(:, :)
    complex(8) :: xi(d, d), xit(d, d)
    real(8) :: dk, ang0, resu
    integer :: info, infot
    ang0 = 0.4d0; dk = 0.1d0
    Qk  = block_states(d, ang0)
    Qkb = block_states(d, ang0 + dk)
    M   = matmul(hc(Qk), Qkb)
    ! two DIFFERENT gauges on the stencil: M~ = Wk^H M Wkb
    Wk  = u2_gauge(0.5d0, -0.4d0, 0.8d0,  0.3d0)
    Wkb = u2_gauge(1.1d0,  0.7d0, 0.3d0, -0.6d0)
    Mt  = matmul(matmul(hc(Wk), M), Wkb)
    call xi_block_from_overlap(M,  dk, xi_sign, xi,  info,  resu)
    call xi_block_from_overlap(Mt, dk, xi_sign, xit, infot)
    rho = herm_fixed(d, 0.0d0)
    call check_true(info == 0 .and. infot == 0, "grad-caveat: xi kernel ok", nfail)
    call check_true(mdiff(comm(xit, rot(Wk, rho)), rot(Wk, comm(xi, rho))) > 1d-1, &
                    "grad-caveat: varying-stencil gauge breaks bare [xi,rho] (grad term, out of scope)", &
                    nfail)
  end subroutine test_grad_caveat

  ! (5) NEGATIVE gate: the current fix uses xi; the legacy dipole d = i p/Delta_omega
  !     is 0/0 in an exactly degenerate block.  The code knocks out near-degenerate
  !     dipoles (|Delta_omega| < 1e-3 -> 0), so the naive block coupling COLLAPSES to
  !     zero (and is +-Inf if left unregularized).  Define the gate:
  !         GATE = (block RHS is nonzero)  .AND.  (block RHS is gauge-covariant).
  !     xi PASSES the gate; d = i p/Delta_omega FAILS it.  This proves the gate
  !     actually detects the bug the xi construction fixes.
  subroutine test_negative_naive_dipole(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: d = 2
    complex(8), allocatable :: W(:, :), Qk(:, :), Qkb(:, :), M(:, :), Mt(:, :)
    complex(8) :: xi(d, d), xit(d, d)
    complex(8) :: rho(d, d), p(d, d), d_naive(d, d), d_unreg(d, d)
    complex(8) :: RHS_xi(d, d), RHS_naive(d, d)
    real(8) :: dk, ang0, resu, dw, cov_xi, cov_naive
    integer :: info, infot, a, b
    logical :: gate_xi, gate_naive, finite_unreg
    ang0 = 0.4d0; dk = 0.1d0
    W   = u2_gauge(0.5d0, -0.4d0, 0.8d0, 0.3d0)
    Qk  = block_states(d, ang0)
    Qkb = block_states(d, ang0 + dk)
    M   = matmul(hc(Qk), Qkb)
    Mt  = rot(W, M)
    rho = herm_fixed(d, 0.0d0)
    p   = 1.70d0 * herm_fixed(d, 0.5d0)

    ! --- xi (the fix): nonzero, covariant ---
    call xi_block_from_overlap(M,  dk, xi_sign, xi,  info,  resu)
    call xi_block_from_overlap(Mt, dk, xi_sign, xit, infot)
    RHS_xi = comm(xi, rho)
    cov_xi = mdiff(comm(xit, rot(W, rho)), rot(W, RHS_xi))
    gate_xi = (mmax(RHS_xi) > 1d-6) .and. (cov_xi < 1d-10)

    ! --- naive d = i p / Delta_omega with exact degeneracy Delta_omega = 0 ---
    dw = 0.30d0 - 0.30d0                       ! runtime zero (avoids compile-time 1/0)
    if (abs(dw) < 1d-3) then                    ! matches the code's dipole knockout
      d_naive = (0d0, 0d0)
    else
      d_naive = zi * p / dw
    end if
    RHS_naive = comm(d_naive, rho)
    cov_naive = mdiff(comm(rot(W, d_naive), rot(W, rho)), rot(W, RHS_naive))  ! =0 (all zero)
    gate_naive = (mmax(RHS_naive) > 1d-6) .and. (cov_naive < 1d-10)

    ! --- unregularized 0/0 is non-finite (the raw math the knockout hides) ---
    d_unreg = zi * p / dw                       ! dw = 0 -> Inf / NaN
    finite_unreg = .true.
    do b = 1, d
      do a = 1, d
        if (.not. (d_unreg(a, b) == d_unreg(a, b) .and. abs(d_unreg(a, b)) <= huge(1d0))) &
          finite_unreg = .false.
      end do
    end do

    call check_true(gate_xi, "negative: GATE PASSES for xi (nonzero & covariant)", nfail)
    call check_true(.not. gate_naive, &
                    "negative: GATE FAILS for d=i p/Delta_omega (degenerate-block coupling collapses)", nfail)
    call check_true(mmax(RHS_xi) > 1d-6, "negative: xi gives NONZERO degenerate-block coupling", nfail)
    call check_close_r(mmax(RHS_naive), 0d0, 1d-30, &
                       "negative: d=i p/Delta_omega gives ZERO degenerate-block coupling", nfail)
    call check_true(.not. finite_unreg, &
                    "negative: unregularized d=i p/Delta_omega is non-finite (0/0 at exact degeneracy)", nfail)
  end subroutine test_negative_naive_dipole

  ! (6) End-to-end build_xi on a multi-k, dispersive spectrum with a genuine
  !     energy CROSSING of the two non-block bands (their order swaps between k=2
  !     and k=3), while bands {2,3} stay an exactly degenerate block at every k.
  !     Apply a (global-constant, hence stencil-constant on the periodic grid)
  !     U(2) gauge to the block via prod_dk~ = V^H prod_dk V and assert the whole
  !     xi field transforms as a tensor:  xi~(:,:,axis,ik) = V^H xi(:,:,axis,ik) V.
  subroutine test_build_xi_covariance(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 4, nk = 4, nbvec = 3
    integer :: num_kgrid(3), bvec(3, nbvec)
    complex(8) :: prod_dk(nb, nb, nbvec, nk), prod_t(nb, nb, nbvec, nk)
    real(8)    :: eigen(nb, nk), b_matrix(3, 3)
    complex(8) :: xi(nb, nb, 3, nk), xit(nb, nb, 3, nk)
    logical    :: xi_ok(nb, nb, nk), xi_okt(nb, nb, nk)
    complex(8), allocatable :: W(:, :), V(:, :), Qk(:, :), Qn(:, :)
    integer :: n_reject, n_reject_t, ik, ikx, axis, a
    real(8) :: g, th0, resmax, worst

    num_kgrid = (/ nk, 1, 1 /)
    bvec(:, 1) = (/ 1, 0, 0 /); bvec(:, 2) = (/ 0, 1, 0 /); bvec(:, 3) = (/ 0, 0, 1 /)
    b_matrix = 0d0; b_matrix(1, 1) = 0.4d0; b_matrix(2, 2) = 1d0; b_matrix(3, 3) = 1d0

    ! spectrum: block {2,3} exactly degenerate; bands 1 & 4 disperse and CROSS
    ! (gaps at grid points stay >> theta_off, so blocks are cleanly {1},{2,3},{4}).
    eigen(1, :) = (/ 0.36d0, 0.42d0, 0.48d0, 0.54d0 /)   ! rising
    eigen(2, :) = 0.30d0
    eigen(3, :) = 0.30d0
    eigen(4, :) = (/ 0.54d0, 0.48d0, 0.42d0, 0.36d0 /)   ! falling -> crosses band 1

    ! overlaps: +x link from smooth frame (block rotation), y/z self-links = I.
    g = 0.05d0; th0 = 0.30d0
    prod_dk = (0d0, 0d0)
    do ik = 1, nk
      ikx = mod(ik, nk) + 1                              ! periodic +x neighbour
      Qk  = frame4(nb, th0 + g * (ik  - 1))
      Qn  = frame4(nb, th0 + g * (ikx - 1))
      prod_dk(:, :, 1, ik) = matmul(hc(Qk), Qn)
      do a = 1, nb
        prod_dk(a, a, 2, ik) = (1d0, 0d0)                ! (0,1,0) self-link = identity
        prod_dk(a, a, 3, ik) = (1d0, 0d0)                ! (0,0,1) self-link = identity
      end do
    end do

    call build_xi(nb, nk, nbvec, bvec, prod_dk, eigen, b_matrix, num_kgrid, &
                  xi, xi_ok, n_reject, resmax)
    call check_int(n_reject, 0, "build_xi: no rejected links", nfail)
    call check_true(xi_ok(2, 3, 1), "build_xi: block {2,3} xi_ok", nfail)
    call check_true(mmax(xi(:, :, 1, 1)) > 1d-6, "build_xi: block xi nonzero on +x axis", nfail)

    ! global-constant U(2) gauge on the block; transform every prod_dk slice
    W = u2_gauge(0.7d0, 0.4d0, 0.6d0, -0.2d0)
    V = embed(nb, 2, W)
    do ik = 1, nk
      do axis = 1, nbvec
        prod_t(:, :, axis, ik) = rot(V, prod_dk(:, :, axis, ik))
      end do
    end do

    call build_xi(nb, nk, nbvec, bvec, prod_t, eigen, b_matrix, num_kgrid, &
                  xit, xi_okt, n_reject_t, resmax)
    call check_int(n_reject_t, 0, "build_xi(gauged): no rejected links", nfail)

    worst = 0d0
    do ik = 1, nk
      do axis = 1, 3
        worst = max(worst, mdiff(xit(:, :, axis, ik), rot(V, xi(:, :, axis, ik))))
      end do
    end do
    call check_close_r(worst, 0d0, 1d-10, &
                       "build_xi: xi field is a U(2) tensor  xi~ = V^H xi V  (all axes/k)", nfail)
  end subroutine test_build_xi_covariance

end program test_gisbe_covariance
