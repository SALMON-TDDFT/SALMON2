! src/ssbe/test/test_xi_construction.f90
! Unit tests for the LG-SBE Tier2 Pb3 non-Abelian Berry-connection kernel
! (xi_block_from_overlap) and the cos^2 blend in degenerate_block_ssbe.
!
! The kernel calls LAPACK (zheev/zgeev), so this links against LAPACK/BLAS.
! Standalone build from src/ssbe:
!     gfortran degenerate_block_ssbe.f90 test/test_xi_construction.f90 \
!              -o t -llapack -lblas && ./t
!
! Model (analytic ground truth): a smooth 2-band rotation
!     |u_1(k)> = ( cos t,  sin t),  |u_2(k)> = (-sin t,  cos t),  t = t0 + g*k
! gives real orthonormal Bloch states with
!     <u_1|d_k u_2> = -g,   <u_2|d_k u_1> = +g   (diagonal 0),
! so the overlap M(a,b) = <u_a(k)|u_b(k+dk)> is the rotation R(g*dk), U=M, and
!     i*logm(U)/dk = i<u_n|d_k u_m>  (the natural Berry connection).
! SALMON's dipole is d_nm = i*p_nm/delta_omega with p_nm = (e_m-e_n)<u_n|d_k u_m>,
! hence d = -i<u_n|d_k u_m>.  Requiring xi = d (drop-in for the blend) pins the
! sign convention to xi_sign = -1, which these tests verify numerically.

program test_xi_construction
  use degenerate_block_ssbe, only: xi_block_from_overlap, xi_sign, blend
  implicit none
  integer :: nfail
  nfail = 0

  call test_sign_and_value(nfail)
  call test_hermitian_and_diag(nfail)
  call test_polar_rightstretch(nfail)
  call test_reject_singular(nfail)
  call test_reject_nearpi(nfail)
  call test_blend_ramp(nfail)

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_xi_construction checks passed."
  end if

contains

  !----------------------------- assert helpers ----------------------------
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

  subroutine check_close_c(got, want, tol, label, nfail)
    complex(8), intent(in) :: got, want
    real(8), intent(in) :: tol
    character(*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (abs(got - want) <= tol) then
      write(*, '(a)') "ok    " // label
    else
      write(*, '(a,es12.4)') "FAIL  " // label // "  |got-want|=", abs(got - want)
      write(*, '(6x,a,2es14.6,a,2es14.6)') "got=", got, "  want=", want
      nfail = nfail + 1
    end if
  end subroutine check_close_c

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

  ! Overlap M(a,b) = <u_a(k)|u_b(k+dk)> of the rotation model (real states).
  subroutine rotation_overlap(g, dk, t0, M)
    real(8), intent(in) :: g, dk, t0
    complex(8), intent(out) :: M(2, 2)
    real(8) :: ta, tb, u1a(2), u2a(2), u1b(2), u2b(2)
    ta = t0
    tb = t0 + g * dk
    u1a = (/  cos(ta),  sin(ta) /); u2a = (/ -sin(ta),  cos(ta) /)
    u1b = (/  cos(tb),  sin(tb) /); u2b = (/ -sin(tb),  cos(tb) /)
    M(1, 1) = dcmplx(dot_product(u1a, u1b), 0d0)
    M(1, 2) = dcmplx(dot_product(u1a, u2b), 0d0)
    M(2, 1) = dcmplx(dot_product(u2a, u1b), 0d0)
    M(2, 2) = dcmplx(dot_product(u2a, u2b), 0d0)
  end subroutine rotation_overlap

  !------------------------------- tests ------------------------------------

  ! (1) Sign + value: in the finite-gap limit xi (from overlaps, s=-1) equals the
  !     analytic dipole i*p/delta_omega; and s=+1 reproduces -that, so s=-1 is the
  !     unique sign making xi a drop-in dipole. This PINS xi_sign = -1.
  subroutine test_sign_and_value(nfail)
    integer, intent(inout) :: nfail
    complex(8), parameter :: zi = (0d0, 1d0)
    real(8) :: g, dk, t0, gap, p12
    complex(8) :: M(2, 2), xi(2, 2), xip(2, 2), d12_ref
    integer :: info
    real(8) :: resu

    g = 0.05d0; dk = 0.10d0; t0 = 0.37d0    ! small link phase g*dk = 5e-3
    gap = 1.0d-3                            ! finite intra-block gap (< theta_off)
    call rotation_overlap(g, dk, t0, M)

    ! analytic SALMON dipole d_12 = i*p_12/delta_omega, p_12 = (e2-e1)<u1|d_k u2>
    !   <u1|d_k u2> = -g,  delta_omega_12 = e1-e2 = gap  ->  p_12 = (-gap)*(-g)
    p12 = gap * g
    d12_ref = zi * p12 / gap                ! = zi*g

    call xi_block_from_overlap(M, dk, xi_sign, xi, info, resu)
    call check_true(info == 0,            "sign_and_value: info==0", nfail)
    call check_true(xi_sign == -1d0,      "sign_and_value: xi_sign pinned to -1", nfail)
    call check_close_c(xi(1, 2), d12_ref, 1d-9, "xi(1,2) == i*p/dw (finite-gap)", nfail)
    call check_close_c(xi(1, 2), zi * g,  1d-9, "xi(1,2) == i*g (analytic)",       nfail)
    call check_close_c(xi(2, 1), -zi * g, 1d-9, "xi(2,1) == -i*g (analytic)",      nfail)

    ! flipping the sign gives -d, confirming s=-1 (not +1) reproduces the dipole
    call xi_block_from_overlap(M, dk, +1d0, xip, info, resu)
    call check_close_c(xip(1, 2), -d12_ref, 1d-9, "s=+1 gives -d (wrong sign)", nfail)
    call check_true(abs(xip(1, 2) - d12_ref) > 0.5d0 * abs(d12_ref), &
                    "s=+1 does NOT match dipole -> s=-1 is unique", nfail)
  end subroutine test_sign_and_value

  ! (2) xi is Hermitian and its block-diagonal (intra-band) part vanishes here.
  subroutine test_hermitian_and_diag(nfail)
    integer, intent(inout) :: nfail
    real(8) :: g, dk, t0
    complex(8) :: M(2, 2), xi(2, 2)
    integer :: info
    g = 0.07d0; dk = 0.13d0; t0 = 1.1d0
    call rotation_overlap(g, dk, t0, M)
    call xi_block_from_overlap(M, dk, xi_sign, xi, info)
    call check_true(info == 0, "hermitian: info==0", nfail)
    call check_close_c(xi(2, 1), conjg(xi(1, 2)), 1d-10, "xi Hermitian (xi21=conj xi12)", nfail)
    call check_close_r(abs(xi(1, 1)), 0d0, 1d-10, "xi(1,1) ~ 0", nfail)
    call check_close_r(abs(xi(2, 2)), 0d0, 1d-10, "xi(2,2) ~ 0", nfail)
  end subroutine test_hermitian_and_diag

  ! (3) Polar step: a non-unitary overlap (right column stretch) must be projected
  !     back onto the same unitary U, so xi is unchanged and U^H U - I ~ 0.
  subroutine test_polar_rightstretch(nfail)
    integer, intent(inout) :: nfail
    complex(8), parameter :: zi = (0d0, 1d0)
    real(8) :: g, dk, t0, resu
    complex(8) :: M(2, 2), xi(2, 2)
    integer :: info
    g = 0.05d0; dk = 0.10d0; t0 = 0.37d0
    call rotation_overlap(g, dk, t0, M)
    M(:, 2) = 1.3d0 * M(:, 2)               ! right stretch -> removed by polar
    call xi_block_from_overlap(M, dk, xi_sign, xi, info, resu)
    call check_true(info == 0, "polar: info==0", nfail)
    call check_true(resu < 1d-8, "polar: |U^H U - I| ~ 0", nfail)
    call check_close_c(xi(1, 2), zi * g, 1d-9, "polar: xi(1,2) == i*g after stretch", nfail)
  end subroutine test_polar_rightstretch

  ! (4) Near-singular overlap (rank-deficient) is rejected with info==1.
  subroutine test_reject_singular(nfail)
    integer, intent(inout) :: nfail
    complex(8) :: M(2, 2), xi(2, 2)
    integer :: info
    M(1, 1) = (1d0, 0d0); M(1, 2) = (1d0, 0d0)
    M(2, 1) = (1d0, 0d0); M(2, 2) = (1d0, 0d0)   ! rank 1 -> sigma_min = 0
    call xi_block_from_overlap(M, 1d0, xi_sign, xi, info)
    call check_true(info == 1, "reject: near-singular -> info==1", nfail)
  end subroutine test_reject_singular

  ! (5) Eigenphase near +-pi (branch cut) is rejected with info==2.
  subroutine test_reject_nearpi(nfail)
    integer, intent(inout) :: nfail
    real(8), parameter :: pi = 3.14159265358979323846d0
    complex(8) :: M(2, 2), xi(2, 2)
    integer :: info
    ! rotation by 0.95*pi -> eigenphases +-0.95*pi > 0.9*pi
    call rotation_overlap(0.95d0 * pi, 1d0, 0.2d0, M)
    call xi_block_from_overlap(M, 1d0, xi_sign, xi, info)
    call check_true(info == 2, "reject: eigenphase near +-pi -> info==2", nfail)
  end subroutine test_reject_nearpi

  ! (6) blend() is the C^1 cos^2 ramp: 1 at/below theta_on, 0 at/above theta_off,
  !     0.5 at the midpoint, monotone in between (no discontinuity at the thetas).
  subroutine test_blend_ramp(nfail)
    integer, intent(inout) :: nfail
    call check_close_r(blend(0.0d0, 1d0, 2d0), 1d0, 1d-14, "blend: x<=on -> 1", nfail)
    call check_close_r(blend(1.0d0, 1d0, 2d0), 1d0, 1d-14, "blend: x==on -> 1", nfail)
    call check_close_r(blend(2.0d0, 1d0, 2d0), 0d0, 1d-14, "blend: x==off -> 0", nfail)
    call check_close_r(blend(3.0d0, 1d0, 2d0), 0d0, 1d-14, "blend: x>=off -> 0", nfail)
    call check_close_r(blend(1.5d0, 1d0, 2d0), 0.5d0, 1d-12, "blend: midpoint -> 0.5", nfail)
    call check_true(blend(1.25d0, 1d0, 2d0) > blend(1.75d0, 1d0, 2d0), &
                    "blend: monotone decreasing", nfail)
  end subroutine test_blend_ramp

end program test_xi_construction
