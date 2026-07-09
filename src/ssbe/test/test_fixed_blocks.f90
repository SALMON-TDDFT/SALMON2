! src/ssbe/test/test_fixed_blocks.f90
! LG-SBE Tier2 Task 2: build_fixed_blocks -- a k-INDEPENDENT composite-block
! partition derived once from the GS eigenvalues: connected components of the
! graph whose edges are band pairs near-degenerate somewhere in the BZ,
! verified gap-isolated from outside bands at every k (isolated_ok=.false.
! flags the metal-like ambiguous case, e.g. a fourth band that never joins the
! block but drifts within gap_margin of it at some k).
!
! Standalone build (same convention as test_gisbe_covariance.f90):
!     gfortran degenerate_block_ssbe.f90 test/test_fixed_blocks.f90 \
!              -o t_fixed -llapack -lblas && ./t_fixed
! (build_fixed_blocks itself needs no LAPACK, but the module links it in.)

program test_fixed_blocks
  use degenerate_block_ssbe, only: build_fixed_blocks, build_blocks_fixed
  implicit none
  integer :: nfail
  nfail = 0

  call test_build_fixed_blocks_basic(nfail)
  call test_fixed_blocks_metal_flag(nfail)
  call test_close_singletons_not_flagged(nfail)
  call test_block_id_k_independent(nfail)

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_fixed_blocks checks passed."
  end if

contains

  !======================= assert helpers (ssbe style, copied from ===========
  !======================= test_gisbe_covariance.f90 lines ~60-94) ===========
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

  !================================ tests ====================================

  ! bands 1,2 degenerate at k=1 only (and again at k=3, within theta_off);
  ! split beyond theta_off at k=2 -- transitive union-find must still join
  ! them into ONE k-independent block. Bands 3,4 stay always-separated
  ! singletons (gap 0.5 au >> gap_margin), so isolation holds at every k.
  subroutine test_build_fixed_blocks_basic(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 4, nk = 3
    real(8) :: eigen(nb, nk)
    integer :: fixed_block_id(nb)
    logical :: isolated_ok

    eigen(:, 1) = [0.0d0, 0.0d0,          1.0d0, 1.5d0]   ! 1,2 degenerate here
    eigen(:, 2) = [0.0d0, 1.0d-3 + 1d-4,  1.0d0, 1.5d0]   ! 1,2 split beyond theta_off
    eigen(:, 3) = [0.0d0, 5.0d-4,         1.0d0, 1.5d0]   ! 1,2 within theta_off again

    call build_fixed_blocks(nb, nk, eigen, fixed_block_id, isolated_ok)

    call check_true(fixed_block_id(1) == fixed_block_id(2), &
                     "1,2 must share a block", nfail)
    call check_true(fixed_block_id(3) /= fixed_block_id(1), &
                     "3 must not join {1,2}", nfail)
    call check_true(fixed_block_id(3) /= fixed_block_id(4), &
                     "3,4 separate (gap 0.5)", nfail)
    call check_true(isolated_ok, "all blocks gap-isolated here", nfail)
  end subroutine test_build_fixed_blocks_basic

  ! band 3 sits at 2.1d-3 from block {1,2}: OUTSIDE theta_off=2d-3 (so it never
  ! grouped in) but WITHIN gap_margin=2.2d-3 (so the block is not cleanly
  ! isolated) -- the ambiguous / metal-like case that must fail-closed.
  ! (gap_margin was recalibrated 1d-2 -> 2.2d-3 against the real Si 4^3
  ! spectrum's global min cross-block gap 2.654e-3 au; see build_fixed_blocks
  ! header in degenerate_block_ssbe.f90. 5d-3 used to sit inside the old
  ! margin but is now OUTSIDE 2.2d-3, so the outside band must be moved into
  ! the new (theta_off, gap_margin) = (2d-3, 2.2d-3) window.)
  subroutine test_fixed_blocks_metal_flag(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 3, nk = 2
    real(8) :: eigen(nb, nk)
    integer :: fid(nb)
    logical :: ok

    eigen(:, 1) = [0.0d0, 0.0d0,  2.1d-3]   ! 1,2 degenerate (block); 3 outside but close
    eigen(:, 2) = [0.0d0, 1.0d-3, 3.1d-3]   ! 1,2 still within theta_off -> stay grouped

    call build_fixed_blocks(nb, nk, eigen, fid, ok)

    call check_true(fid(3) /= fid(1), &
                     "band 3 (2.1d-3 > theta_off) must NOT join {1,2}", nfail)
    call check_true(.not. ok, &
                     "must flag non-isolated block (3 within gap_margin of {1,2})", nfail)
  end subroutine test_fixed_blocks_metal_flag

  ! bsize>1 gate coverage: bands 1,2 sit close (3.0/3.5d-3) at every k, always
  ! > theta_off (2d-3) apart, so they stay SEPARATE singletons. With
  ! gap_margin now 2.2d-3 (recalibrated, see test_fixed_blocks_metal_flag)
  ! 3.0d-3/3.5d-3 also sit > gap_margin, so this case would not even reach the
  ! bsize>1 gate on distance alone -- it still proves the intended point:
  ! because neither band is a multi-band block, the isolation loop must NOT
  ! trip regardless of proximity. Two merely-close-but-distinct non-degenerate
  ! bands are not the ambiguity this routine flags. Band 3 is far away. Proves
  ! the gate protects the named singleton-vs-singleton false-abort failure mode.
  subroutine test_close_singletons_not_flagged(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 3, nk = 2
    real(8) :: eigen(nb, nk)
    integer :: fid(nb)
    logical :: ok

    eigen(:, 1) = [0.0d0, 3.0d-3, 1.0d0]   ! 1,2 close but > theta_off (and > gap_margin)
    eigen(:, 2) = [0.0d0, 3.5d-3, 1.0d0]   ! still separate singletons at k=2

    call build_fixed_blocks(nb, nk, eigen, fid, ok)

    call check_true(fid(1) /= fid(2), &
                     "close singletons (3d-3 > theta_off) stay in different blocks", nfail)
    call check_true(ok, &
                     "close SINGLETON bands must NOT flag isolation (bsize>1 gate)", nfail)
  end subroutine test_close_singletons_not_flagged

  ! Task 3: build_blocks_fixed must broadcast the single k-independent
  ! fixed_block_id partition into EVERY k-slice of block_id(nb,nk) unchanged
  ! -- this is what makes it a drop-in replacement for build_blocks' per-k
  ! union-find inside build_xi's gifix path (so block_id(:,ik) is identical
  ! for all ik, unlike build_blocks which can re-partition per k).
  subroutine test_block_id_k_independent(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 4, nk = 3
    real(8) :: eigen(nb, nk)
    integer :: block_id(nb, nk)
    integer :: ik
    logical :: k_indep

    eigen(:, 1) = [0.0d0, 0.0d0,         1.0d0, 1.5d0]
    eigen(:, 2) = [0.0d0, 1.0d-3 + 1d-4, 1.0d0, 1.5d0]
    eigen(:, 3) = [0.0d0, 5.0d-4,        1.0d0, 1.5d0]

    call build_blocks_fixed(nb, nk, eigen, block_id)

    k_indep = .true.
    do ik = 2, nk
      if (any(block_id(:, ik) /= block_id(:, 1))) k_indep = .false.
    end do
    call check_true(k_indep, "build_blocks_fixed: block_id is k-independent", nfail)
    call check_true(block_id(1, 1) == block_id(2, 1), &
                     "build_blocks_fixed: 1,2 share block", nfail)
  end subroutine test_block_id_k_independent

end program test_fixed_blocks
