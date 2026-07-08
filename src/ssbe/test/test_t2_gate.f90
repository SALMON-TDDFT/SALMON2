! src/ssbe/test/test_t2_gate.f90
! T2 Delta-omega dephasing gate: runtime-configurable shape (hard 'step' /
! smooth Gaussian-notch 'gauss') + its own namelist parameters.
! Design: DFT/gw/gw_design/plans/2026-07-08-t2-gate-shape.md
!
! WHY THIS FILE EXISTS
!   The gicov length-gauge SBE's phenomenological T2 dephasing term was
!   hard-gated by the SHARED block-partition threshold theta_off (see
!   bloch_solver_ssbe.f90's gicov_rhs, degenerate_block_ssbe.f90).  This adds
!   THREE new &sbe keys (sbe_t2_gate_shape/theta/width) read ONLY by the T2
!   gate -- theta_off/theta_on/gap_margin (block partition / gi-blend /
!   VG-guard) are untouched.  Task 1 proves the plumbing defaults reproduce
!   the current behaviour bit-for-bit (shape='step', theta=2d-3 == the old
!   theta_off literal, width=0).  Task 2 unit-tests the pure weight helper
!   t2_gate_weight for both shapes + the exact-degeneracy floor clamp.
!   Task 3 unit-tests the checker predicate t2_gate_params_ok.  Task 4 proves
!   (a) shape='step' reproduces gicov_rhs's CURRENT hard-gate output bit-for-
!   bit (regression) and (b) shape='gauss' leaves an exactly-degenerate pair
!   undephased and keeps the dephasing increment Hermitian.
!
! BUILD (mirrors test_gicov_rhs.f90 -- already-built ninja tree at
! build_local/; single-process communication dummy).  Links the SAME objects
! the salmon executable built, minus main.f90.o:
!
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     -I<repo>/src/ssbe -I<repo>/build_local -J<scratch_dir> \
!     -c <repo>/src/ssbe/test/test_t2_gate.f90 -o <scratch_dir>/test_t2_gate.o
!
!   gfortran -fopenmp -cpp -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     $(find <repo>/build_local/src/CMakeFiles/salmon.dir -name '*.o' ! -name 'main.f90.o') \
!     <scratch_dir>/test_t2_gate.o -o <scratch_dir>/test_t2_gate \
!     -framework Accelerate -lm -ldl
!
!   cd <scratch_dir> && <scratch_dir>/test_t2_gate
!     (must run from a writable cwd: test_defaults creates/reads .namelist.tmp
!     via the REAL read_input_common, so it genuinely exercises the &sbe
!     defaults block in inputoutput.f90, not a hand-set fixture.)
!
program test_t2_gate
  use salmon_global, only: sbe_t2_gate_shape, sbe_t2_gate_theta, sbe_t2_gate_width
  use inputoutput,   only: read_input_common
  use gs_info_ssbe,  only: t2_gate_weight
  use input_checker_sbe, only: t2_gate_params_ok
  implicit none

  integer :: nfail
  nfail = 0

  call test_defaults(nfail)
  call test_gate_weight(nfail)
  call test_gate_params_ok(nfail)

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_t2_gate checks passed."
  end if

contains

  !======================= Task 1: namelist defaults ===========================
  ! Runs the REAL read_input_common with an empty .namelist.tmp (no &sbe group
  ! present at all) so this genuinely exercises the inputoutput.f90 defaults
  ! block, not a hand-set fixture (unlike the other ssbe standalone tests,
  ! which never run read_input -- see test_sbe_spinor.f90's fixture note).
  subroutine test_defaults(nfail)
    implicit none
    integer, intent(inout) :: nfail
    integer :: fh

    open(newunit=fh, file='.namelist.tmp', status='replace')
    close(fh)
    call read_input_common()

    call check_true(trim(sbe_t2_gate_shape) == 'step', 'default shape=step', nfail)
    call check_close_r(sbe_t2_gate_theta, 2d-3, 1d-15, 'default theta=2e-3', nfail)
    call check_close_r(sbe_t2_gate_width, 0d0,  1d-15, 'default width=0', nfail)
  end subroutine test_defaults

  !======================= Task 2: t2_gate_weight helper =======================
  ! Property tests for both shapes, per the "Gate weight definition" table
  ! (plan Sec. "Gate weight definition (single source of truth)"): step is a
  ! strict '>' hard gate; gauss is the Gaussian notch W=1-exp(-(dw/w)^2)
  ! (Thuemmler Eq.28, W(0)=0, quadratic onset, saturates to 1). Both shapes
  ! clamp to 0 for |delta_omega| <= floor (exact-degeneracy protection).
  subroutine test_gate_weight(nfail)
    implicit none
    integer, intent(inout) :: nfail

    ! step: strict > theta, exact-degeneracy 0
    call check_close_r(t2_gate_weight( 3d-3, 'step', 2d-3, 0d0, 1d-9), 1d0, 1d-15, &
      'step above', nfail)
    call check_close_r(t2_gate_weight( 1d-3, 'step', 2d-3, 0d0, 1d-9), 0d0, 1d-15, &
      'step below', nfail)
    call check_close_r(t2_gate_weight( 2d-3, 'step', 2d-3, 0d0, 1d-9), 0d0, 1d-15, &
      'step at-threshold=skip (strict >)', nfail)
    call check_close_r(t2_gate_weight( 0d0 , 'step', 2d-3, 0d0, 1d-9), 0d0, 1d-15, &
      'step at 0 = 0', nfail)

    ! gauss: W(0)=0, symmetric, saturates
    call check_close_r(t2_gate_weight( 0d0 , 'gauss', 0d0, 1d-3, 1d-9), 0d0, 1d-12, &
      'gauss at 0 = 0', nfail)
    call check_close_r(t2_gate_weight( 1d-3, 'gauss', 0d0, 1d-3, 1d-9), 1d0 - exp(-1d0), 1d-12, &
      'gauss at w', nfail)
    call check_close_r(t2_gate_weight( 1d-3, 'gauss', 0d0, 1d-3, 1d-9), &
                        t2_gate_weight(-1d-3, 'gauss', 0d0, 1d-3, 1d-9), 1d-15, &
      'gauss symmetric', nfail)
    call check_true(t2_gate_weight( 5d-3, 'gauss', 0d0, 1d-3, 1d-9) > 0.99d0, &
      'gauss saturates', nfail)

    ! floor clamp (both shapes): |dw|<=floor -> 0
    call check_close_r(t2_gate_weight( 5d-10, 'gauss', 0d0, 1d-3, 1d-9), 0d0, 1d-15, &
      'gauss floor clamp', nfail)
  end subroutine test_gate_weight

  !======================= Task 3: checker predicate t2_gate_params_ok ========
  ! shape must be 'step' or 'gauss' (case-sensitive: production never
  ! lowercases sbe_deph_mode-style string keys either -- see
  ! inputoutput.f90's "convert lowercase" block, which omits sbe_* strings);
  ! theta/width must be non-negative; gauss additionally requires width>0
  ! (width=0 would divide by zero in t2_gate_weight's gauss branch).
  subroutine test_gate_params_ok(nfail)
    implicit none
    integer, intent(inout) :: nfail

    call check_true(      t2_gate_params_ok('step' , 2d-3, 0d0 ), 'valid step ok', nfail)
    call check_true(      t2_gate_params_ok('gauss', 0d0 , 1d-3), 'valid gauss ok', nfail)
    call check_true(.not. t2_gate_params_ok('bogus', 2d-3, 0d0 ), 'bad shape rejected', nfail)
    call check_true(.not. t2_gate_params_ok('gauss', 0d0 , 0d0 ), 'gauss width=0 rejected', nfail)
    call check_true(.not. t2_gate_params_ok('step' ,-1d0 , 0d0 ), 'negative theta rejected', nfail)
  end subroutine test_gate_params_ok

  !======================= assert helpers (ssbe style, copied verbatim
  !======================= from test/test_block_transport.f90) ================
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

  subroutine check_close_r(got, want, tol, label, nfail)
    implicit none
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
    implicit none
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

end program test_t2_gate
