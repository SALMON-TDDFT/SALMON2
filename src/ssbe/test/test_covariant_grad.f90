! src/ssbe/test/test_covariant_grad.f90
! LG-SBE Tier2 Phase 2 -- unit tests for the pure gauge-covariant intraband
! k-derivative operator covariant_grad_block (degenerate_block_ssbe):
!
!     D_cov rho = d_k rho - i[xi, rho]
!
! computed on the FULL nb x nb density with the block-diagonal Wilson-line
! transport U_full (build_block_transport, Phase 1) parallel-transporting the
! neighbour densities into the k-frame before a 4-shell (+-1..+-4) central
! difference, weights c = (4/5,-1/5,4/105,-1/280) -- identical shells/weights
! to the production length-gauge k-gradient grad_k_array_nb2d_dcomplex
! (common_ssbe.f90:104-194, weights src/common/initialization.f90:1343).
!
! Test A -- U_full == I reduces the operator to the BARE 4-shell central
!           difference of rho (hand-coded reference, no common_ssbe link).
! Test B -- gauge COVARIANCE: under a per-k block gauge W(k) applied to both
!           rho -> W^H rho W and U_full -> W(k)^H U_full W(k+e), the operator
!           output transforms as a rank-2 tensor  Dq_g(k) = W(k)^H Dq(k) W(k).
! Test C -- R-5 SIGN gate: constant analytic xi = xi0*sigma_x on block {2,3},
!           U_full = exp(-i xi dk); the operator must reproduce
!           d_k rho - i[xi,rho] (d_k rho via the SAME stencil) to 1e-8, with
!           the "U" orientation (U_full fed as-is). The conjugate-transpose
!           (U^H) orientation is asserted to give the WRONG sign (>1e-4), so
!           the gate genuinely discriminates the commutator sign.
!
! Standalone build from src/ssbe (macOS: Accelerate provides zheev/zgemm):
!     gfortran degenerate_block_ssbe.f90 test/test_covariant_grad.f90 \
!              -o /tmp/t_cg -framework Accelerate && /tmp/t_cg
! (NOT -llapack -lblas; do not link common_ssbe.f90.)

program test_covariant_grad
  use degenerate_block_ssbe, only: covariant_grad_block
  implicit none
  complex(8), parameter :: zi = (0d0, 1d0)
  real(8),    parameter :: pi = 3.14159265358979323846d0
  integer :: nfail
  nfail = 0

  call test_A_bare_stencil(nfail)
  call test_B_covariance(nfail)
  call test_C_sign_gate(nfail)

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_covariant_grad checks passed."
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
      write(*, '(a,es12.4)') "ok    " // label // "  resid=", abs(got - want)
    else
      write(*, '(a,es12.4)') "FAIL  " // label // "  |got-want|=", abs(got - want)
      nfail = nfail + 1
    end if
  end subroutine check_close_r

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

  ! Embed a d x d block gauge W at bands i0..i0+d-1 of an nb x nb identity.
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

  ! Periodic k-neighbour on the uniform grid num_kgrid, shift s along axis.
  ! Same index ordering as build_ik_neighbor (i1 fastest).
  integer function nbr(ik, axis, s, num_kgrid) result(j)
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

  ! Smooth, reproducible, k-dependent Hermitian nb x nb field ("density").
  ! t in [0,1); off-diagonal entries are nonzero (so [xi,rho] does not vanish
  ! on any block).
  subroutine smooth_herm(nb, t, H)
    integer,    intent(in)  :: nb
    real(8),    intent(in)  :: t
    complex(8), intent(out) :: H(nb, nb)
    integer :: a, b
    real(8) :: tp
    tp = 2d0 * pi * t
    do a = 1, nb
      H(a, a) = dcmplx(1d0 + 0.5d0 * a + 0.3d0 * cos(tp + 0.2d0 * a), 0d0)
      do b = a + 1, nb
        H(a, b) = dcmplx(0.4d0 * cos(tp + 0.1d0 * (a + b)), &
                         0.3d0 * sin(tp + 0.15d0 * (a - b)))
        H(b, a) = conjg(H(a, b))
      end do
    end do
  end subroutine smooth_herm

  ! E = exp(-i H dk) for Hermitian H (n x n) via zheev (H = V w V^H).
  subroutine expm_mi_h_dk(H, n, dkv, E)
    integer,    intent(in)  :: n
    complex(8), intent(in)  :: H(n, n)
    real(8),    intent(in)  :: dkv
    complex(8), intent(out) :: E(n, n)
    complex(8) :: A(n, n), Vd(n, n)
    complex(8) :: cwork(4 * n)
    real(8)    :: w(n), rwork(4 * n)
    integer    :: info, k
    A = H
    call zheev('V', 'U', n, A, n, w, cwork, 4 * n, rwork, info)
    if (info /= 0) then
      write(*, '(a,i0)') "expm_mi_h_dk: zheev info=", info
      stop 2
    end if
    do k = 1, n
      Vd(:, k) = A(:, k) * exp(-zi * w(k) * dkv)
    end do
    E = matmul(Vd, conjg(transpose(A)))
  end subroutine expm_mi_h_dk

  real(8) function maxabs4(A)                     ! max |A| over a 4-D array
    complex(8), intent(in) :: A(:, :, :, :)
    maxabs4 = maxval(abs(A))
  end function maxabs4

  !================================ tests ====================================

  ! Test A -- U_full == I reduces the operator to the bare 4-shell central
  !           difference of rho (hand-coded reference).
  subroutine test_A_bare_stencil(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 4, nk = 12, nbvec = 3
    integer :: num_kgrid(3), bvec(3, nbvec)
    complex(8) :: U_full(nb, nb, 3, nk), rho(nb, nb, nk)
    complex(8) :: Dq(nb, nb, 3, nk), ref(nb, nb, 3, nk)
    real(8)    :: dk(3), c(4), bx
    integer    :: ik, axis, m, n, kp, km

    num_kgrid = (/ nk, 1, 1 /)
    bvec(:, 1) = (/ 1, 0, 0 /); bvec(:, 2) = (/ 0, 1, 0 /); bvec(:, 3) = (/ 0, 0, 1 /)
    bx = 0.4d0
    dk = (/ bx / nk, 1d0, 1d0 /)
    c  = (/ 4d0/5d0, -1d0/5d0, 4d0/105d0, -1d0/280d0 /)

    do ik = 1, nk
      call smooth_herm(nb, dble(ik - 1) / dble(nk), rho(:, :, ik))
    end do

    ! U_full = identity everywhere (all axes, all k)
    U_full = (0d0, 0d0)
    do ik = 1, nk
      do axis = 1, 3
        do n = 1, nb
          U_full(n, n, axis, ik) = (1d0, 0d0)
        end do
      end do
    end do

    call covariant_grad_block(nb, nk, nbvec, bvec, num_kgrid, U_full, rho, dk, Dq)

    ! hand-coded bare central difference (periodic), all axes
    ref = (0d0, 0d0)
    do axis = 1, 3
      do ik = 1, nk
        do m = 1, 4
          kp = nbr(ik, axis,  m, num_kgrid)
          km = nbr(ik, axis, -m, num_kgrid)
          ref(:, :, axis, ik) = ref(:, :, axis, ik) &
            + c(m) * (rho(:, :, kp) - rho(:, :, km)) / dk(axis)
        end do
      end do
    end do

    call check_close_r(maxabs4(Dq - ref), 0d0, 1d-12, &
      "A: U=I reduces to bare 4-shell central diff", nfail)
    call check_true(maxval(abs(Dq(:, :, 1, :))) > 1d-3, &
      "A: gradient nonzero on x-axis (nontrivial reference)", nfail)
  end subroutine test_A_bare_stencil

  ! Test B -- gauge covariance under a per-k block gauge W(k).
  subroutine test_B_covariance(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 4, nk = 12, nbvec = 3
    integer :: num_kgrid(3), bvec(3, nbvec)
    complex(8) :: U_full(nb, nb, 3, nk), Ug(nb, nb, 3, nk)
    complex(8) :: rho(nb, nb, nk), rhog(nb, nb, nk)
    complex(8) :: Dq(nb, nb, 3, nk), Dqg(nb, nb, 3, nk)
    complex(8) :: Wemb(nb, nb, nk)
    complex(8), allocatable :: Rk(:, :), Wk(:, :)
    real(8) :: dk(3), bx, worst
    integer :: ik, axis, n, ikp

    num_kgrid = (/ nk, 1, 1 /)
    bvec(:, 1) = (/ 1, 0, 0 /); bvec(:, 2) = (/ 0, 1, 0 /); bvec(:, 3) = (/ 0, 0, 1 /)
    bx = 0.4d0
    dk = (/ bx / nk, 1d0, 1d0 /)

    do ik = 1, nk
      call smooth_herm(nb, dble(ik - 1) / dble(nk), rho(:, :, ik))
    end do

    ! smooth block transport on {2,3} (axis 1); identity on axes 2,3
    U_full = (0d0, 0d0)
    do ik = 1, nk
      Rk = u2_gauge(0.20d0 + 0.13d0 * ik, -0.10d0 + 0.07d0 * ik, &
                    0.50d0 + 0.05d0 * ik,  0.15d0 * ik)
      U_full(:, :, 1, ik) = embed(nb, 2, Rk)
      do axis = 2, 3
        do n = 1, nb
          U_full(n, n, axis, ik) = (1d0, 0d0)
        end do
      end do
    end do

    call covariant_grad_block(nb, nk, nbvec, bvec, num_kgrid, U_full, rho, dk, Dq)

    ! genuinely k-dependent U(2) gauge on block {2,3}
    do ik = 1, nk
      Wk = u2_gauge(0.30d0 + 0.11d0 * ik, 0.20d0 - 0.09d0 * ik, &
                    0.40d0 + 0.06d0 * ik, -0.12d0 * ik)
      Wemb(:, :, ik) = embed(nb, 2, Wk)
    end do

    do ik = 1, nk
      rhog(:, :, ik) = rot(Wemb(:, :, ik), rho(:, :, ik))
    end do

    ! U_full(k) -> W(k)^H U_full(k) W(k+e_axis)
    Ug = (0d0, 0d0)
    do ik = 1, nk
      do axis = 1, 3
        ikp = nbr(ik, axis, 1, num_kgrid)
        Ug(:, :, axis, ik) = matmul(matmul(hc(Wemb(:, :, ik)), &
                                    U_full(:, :, axis, ik)), Wemb(:, :, ikp))
      end do
    end do

    call covariant_grad_block(nb, nk, nbvec, bvec, num_kgrid, Ug, rhog, dk, Dqg)

    worst = 0d0
    do ik = 1, nk
      do axis = 1, 3
        worst = max(worst, maxval(abs(Dqg(:, :, axis, ik) &
                                    - rot(Wemb(:, :, ik), Dq(:, :, axis, ik)))))
      end do
    end do
    call check_close_r(worst, 0d0, 1d-10, &
      "B: Dq covariant  Dq_g = W(k)^H Dq W(k)  (all k, all axes)", nfail)
    call check_true(maxval(abs(Dq(:, :, 1, :))) > 1d-3, &
      "B: ungauged gradient nonzero", nfail)
  end subroutine test_B_covariance

  ! Test C -- R-5 sign gate.  Constant analytic xi = xi0*sigma_x on block {2,3};
  ! U_full = exp(-i xi dk).  Operator must equal d_k rho - i[xi,rho] (d_k via
  ! the same stencil) to 1e-8 in the "U" orientation.  The U^H orientation
  ! must give the WRONG commutator sign (residual O(|[xi,rho]|)).
  subroutine test_C_sign_gate(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 4, nk = 32, nbvec = 3
    integer :: num_kgrid(3), bvec(3, nbvec)
    complex(8) :: U_full(nb, nb, 3, nk), U_H(nb, nb, 3, nk)
    complex(8) :: rho(nb, nb, nk), Dq(nb, nb, 3, nk), DqH(nb, nb, 3, nk)
    complex(8) :: tgt(nb, nb, 3, nk)
    complex(8) :: xi_blk(2, 2), U2(2, 2), xifull(nb, nb)
    real(8) :: dk(3), c(4), bx, xi0
    integer :: ik, axis, m, n, kp, km

    num_kgrid = (/ nk, 1, 1 /)
    bvec(:, 1) = (/ 1, 0, 0 /); bvec(:, 2) = (/ 0, 1, 0 /); bvec(:, 3) = (/ 0, 0, 1 /)
    bx = 0.4d0
    dk = (/ bx / nk, 1d0, 1d0 /)
    c  = (/ 4d0/5d0, -1d0/5d0, 4d0/105d0, -1d0/280d0 /)
    xi0 = 0.02d0

    ! xi = xi0 * sigma_x on block {2,3}
    xi_blk = (0d0, 0d0)
    xi_blk(1, 2) = dcmplx(xi0, 0d0)
    xi_blk(2, 1) = dcmplx(xi0, 0d0)
    ! U2 = exp(-i xi_blk dk(1))  (same constant link at every k, axis 1)
    call expm_mi_h_dk(xi_blk, 2, dk(1), U2)

    U_full = (0d0, 0d0)
    do ik = 1, nk
      U_full(:, :, 1, ik) = embed(nb, 2, U2)
      do axis = 2, 3
        do n = 1, nb
          U_full(n, n, axis, ik) = (1d0, 0d0)
        end do
      end do
    end do

    do ik = 1, nk
      call smooth_herm(nb, dble(ik - 1) / dble(nk), rho(:, :, ik))
    end do

    ! ---- "U" orientation (production) ----
    call covariant_grad_block(nb, nk, nbvec, bvec, num_kgrid, U_full, rho, dk, Dq)

    ! analytic target: d_k rho (same stencil) - i[xi,rho] on axis 1; axes 2,3 = 0.
    ! xi is an ADDITIVE connection -> ZERO outside the block (NOT embed, which
    ! identity-pads; that is only correct for a multiplicative unitary like U).
    xifull = (0d0, 0d0)
    xifull(2, 2) = xi_blk(1, 1); xifull(2, 3) = xi_blk(1, 2)
    xifull(3, 2) = xi_blk(2, 1); xifull(3, 3) = xi_blk(2, 2)
    tgt = (0d0, 0d0)
    do axis = 1, 3
      do ik = 1, nk
        do m = 1, 4
          kp = nbr(ik, axis,  m, num_kgrid)
          km = nbr(ik, axis, -m, num_kgrid)
          tgt(:, :, axis, ik) = tgt(:, :, axis, ik) &
            + c(m) * (rho(:, :, kp) - rho(:, :, km)) / dk(axis)
        end do
      end do
    end do
    do ik = 1, nk
      tgt(:, :, 1, ik) = tgt(:, :, 1, ik) - zi * comm(xifull, rho(:, :, ik))
    end do

    call check_close_r(maxabs4(Dq - tgt), 0d0, 1d-8, &
      "C: Dq = d_k rho - i[xi,rho]  (U orientation, SIGN gate)", nfail)
    call check_true(maxval(abs(comm(xifull, rho(:, :, 1)))) > 1d-4, &
      "C: commutator [xi,rho] nonzero (sign genuinely under test)", nfail)

    ! ---- U^H orientation must give the WRONG sign (gate discriminates) ----
    U_H = (0d0, 0d0)
    do ik = 1, nk
      do axis = 1, 3
        U_H(:, :, axis, ik) = hc(U_full(:, :, axis, ik))
      end do
    end do
    call covariant_grad_block(nb, nk, nbvec, bvec, num_kgrid, U_H, rho, dk, DqH)
    call check_true(maxabs4(DqH - tgt) > 1d-4, &
      "C: U^H orientation gives WRONG commutator sign (>1e-4) -> gate real", nfail)
  end subroutine test_C_sign_gate

end program test_covariant_grad
