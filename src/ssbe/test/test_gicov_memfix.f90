! src/ssbe/test/test_gicov_memfix.f90
! LG-SBE gicov memory-scale bugfix: unit tests for the three pure helpers
! (degenerate_block_ssbe.f90) that gs_info_ssbe's read_prod_dk_data uses to
! (1) keep only the +x/+y/+z unit-shift columns of gs%prod_dk for
!     sbe_lg_degen='gicov' instead of the full (2*ndk+1)**3 shell -- the
!     dominant per-rank memory term (16*nbvec*nb_eff**2*nk bytes) that SIGBUS'd
!     at high nproc_k / large nk (k20^3, k24^3), because build_block_transport
!     (the ONLY prod_dk consumer on the gicov/X-full path) only ever reads the
!     3 +unit-shift columns (find_bvec(bvec,nbvec,1,0,0) etc.), never anything
!     else, at any ndk -- and
! (2) compute the file's expected record count in integer(8) so it stays exact
!     past 2**31-1, which a default-integer nk*nbvec*file_no**2 product
!     already overflows at a k20^3=8000 k-point mesh for a realistic exported
!     band count (silently wrapping negative, turning the reader's
!     `do irec=1,nrec` loop into a ZERO-iteration no-op with no error raised --
!     silent data corruption, not a crash).
!
! gi/gifix/off are untouched by design (build_xi's match_link_blocks scans
! every bvec column, so they still need the full shell) -- Test C below pins
! that gicov_prod_dk_nbvec is a no-op identity for every one of those modes.
!
! Pure module test, no SALMON dependency (same convention as
! test_block_transport.f90 / test_degenerate_block.f90):
!     gfortran degenerate_block_ssbe.f90 test/test_gicov_memfix.f90 \
!              -o t -framework Accelerate && ./t

program test_gicov_memfix
  use degenerate_block_ssbe, only: prod_dk_axis_slot, gicov_prod_dk_nbvec, expected_prod_dk_nrec
  implicit none
  integer :: nfail
  nfail = 0

  call test_axis_slot_kept_directions(nfail)
  call test_axis_slot_dropped_directions(nfail)
  call test_axis_slot_full_ndk1_shell_census(nfail)
  call test_nbvec_gicov_reduces(nfail)
  call test_nbvec_other_modes_unchanged(nfail)
  call test_nrec_small_case_matches_naive(nfail)
  call test_nrec_k20_overflow_class(nfail)
  call test_nrec_k24_overflow_class(nfail)
  call test_memory_reduction_estimate(nfail)

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_gicov_memfix checks passed."
  end if

contains

  !======================= assert helpers (ssbe style, copied verbatim
  !======================= from test/test_block_transport.f90) ===========
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

  subroutine check_int8(got, want, label, nfail)
    integer(8), intent(in) :: got, want
    character(*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (got == want) then
      write(*, '(a)') "ok    " // label
    else
      write(*, '(a,2(1x,i0))') "FAIL  " // label // "  got/want=", got, want
      nfail = nfail + 1
    end if
  end subroutine check_int8

  !======================= Test A: prod_dk_axis_slot, kept directions ======
  ! The 3 unit shifts build_block_transport actually reads must map to 3
  ! DISTINCT nonzero slots (no collision -- a collision would silently drop
  ! one axis's data by overwriting it with another's).
  subroutine test_axis_slot_kept_directions(nfail)
    integer, intent(inout) :: nfail
    integer :: sx, sy, sz
    sx = prod_dk_axis_slot(1, 0, 0)
    sy = prod_dk_axis_slot(0, 1, 0)
    sz = prod_dk_axis_slot(0, 0, 1)
    call check_true(sx > 0 .and. sy > 0 .and. sz > 0, &
      "A: +x/+y/+z all map to a kept (nonzero) slot", nfail)
    call check_true(sx /= sy .and. sy /= sz .and. sx /= sz, &
      "A: +x/+y/+z map to 3 DISTINCT slots (no collision)", nfail)
    call check_true(max(sx, sy, sz) <= 3 .and. min(sx, sy, sz) >= 1, &
      "A: kept slots are exactly {1,2,3}", nfail)
  end subroutine test_axis_slot_kept_directions

  !======================= Test B: prod_dk_axis_slot, dropped directions ===
  ! Everything build_block_transport never reads (negative shifts, diagonals,
  ! the zero shift, multi-shifts) must map to slot 0 (dropped).
  subroutine test_axis_slot_dropped_directions(nfail)
    integer, intent(inout) :: nfail
    call check_int(prod_dk_axis_slot(0, 0, 0), 0, "B: zero shift (0,0,0) dropped", nfail)
    call check_int(prod_dk_axis_slot(-1, 0, 0), 0, "B: -x (-1,0,0) dropped (never read by build_block_transport)", nfail)
    call check_int(prod_dk_axis_slot(0, -1, 0), 0, "B: -y (0,-1,0) dropped", nfail)
    call check_int(prod_dk_axis_slot(0, 0, -1), 0, "B: -z (0,0,-1) dropped", nfail)
    call check_int(prod_dk_axis_slot(1, 1, 0), 0, "B: diagonal (1,1,0) dropped", nfail)
    call check_int(prod_dk_axis_slot(1, 1, 1), 0, "B: diagonal (1,1,1) dropped", nfail)
    call check_int(prod_dk_axis_slot(2, 0, 0), 0, "B: multi-shift (2,0,0) dropped (ndk>=2 case)", nfail)
    call check_int(prod_dk_axis_slot(1, 0, -1), 0, "B: mixed (1,0,-1) dropped", nfail)
  end subroutine test_axis_slot_dropped_directions

  !======================= Test C: full ndk=1 27-shell census ==============
  ! Exhaustively walk the SAME (jdk1,jdk2,jdk3) in -1..1 enumeration
  ! read_prod_dk_data's writer/reader uses and confirm EXACTLY 3 of the 27
  ! combinations are kept, matching the "27 -> 3" reduction claim precisely.
  subroutine test_axis_slot_full_ndk1_shell_census(nfail)
    integer, intent(inout) :: nfail
    integer :: jdk1, jdk2, jdk3, islot, n_kept
    n_kept = 0
    do jdk3 = -1, 1
      do jdk2 = -1, 1
        do jdk1 = -1, 1
          islot = prod_dk_axis_slot(jdk1, jdk2, jdk3)
          if (islot > 0) n_kept = n_kept + 1
        end do
      end do
    end do
    call check_int(n_kept, 3, "C: exactly 3 of the full ndk=1 27-direction shell are kept", nfail)
  end subroutine test_axis_slot_full_ndk1_shell_census

  !======================= Test D: gicov_prod_dk_nbvec reduces for gicov ===
  subroutine test_nbvec_gicov_reduces(nfail)
    integer, intent(inout) :: nfail
    call check_int(gicov_prod_dk_nbvec('gicov', 1, 27), 3, &
      "D: gicov, ndk=1 (nbvec_full=27) -> stored nbvec=3", nfail)
    call check_int(gicov_prod_dk_nbvec('gicov', 2, 125), 3, &
      "D: gicov, ndk=2 (nbvec_full=125) -> stored nbvec=3 (generalizes beyond ndk=1)", nfail)
    call check_int(gicov_prod_dk_nbvec('gicov', 0, 1), 1, &
      "D: gicov, ndk=0 (nbvec_full=1, no +unit-shift data exists) -> unchanged (no-op)", nfail)
  end subroutine test_nbvec_gicov_reduces

  !======================= Test E: other modes are a no-op (unchanged) =====
  subroutine test_nbvec_other_modes_unchanged(nfail)
    integer, intent(inout) :: nfail
    call check_int(gicov_prod_dk_nbvec('gi', 1, 27), 27, "E: gi, ndk=1 -> full 27 kept (unchanged)", nfail)
    call check_int(gicov_prod_dk_nbvec('gifix', 1, 27), 27, "E: gifix, ndk=1 -> full 27 kept (unchanged)", nfail)
    call check_int(gicov_prod_dk_nbvec('off', 1, 27), 27, "E: off, ndk=1 -> full 27 kept (unchanged)", nfail)
  end subroutine test_nbvec_other_modes_unchanged

  !======================= Test F: nrec small case matches naive product ===
  subroutine test_nrec_small_case_matches_naive(nfail)
    integer, intent(inout) :: nfail
    integer(8) :: got
    got = expected_prod_dk_nrec(8, 27, 10)
    call check_int8(got, 21600_8, "F: small case nk=8,nbvec=27,no=10 -> nrec=21600 (sanity, still int32-safe)", nfail)
  end subroutine test_nrec_small_case_matches_naive

  !======================= Test G: k20-scale overflow class (the actual bug)
  ! nk=8000 (20^3), nbvec_full=27 (ndk=1), file_no=100: matches the task's own
  ! "k20 -> 2.16e9 > 2**31" note exactly. huge(0_4) = 2**31-1 = 2147483647;
  ! the true value 2,160,000,000 exceeds it by ~1.25e7 -- a default 4-byte
  ! INTEGER product here silently wraps negative (Fortran integer overflow is
  ! undefined behavior, so this test asserts the DECISIVE, portable fact --
  ! the true value doesn't fit in 32 bits -- via huge(), rather than
  ! performing the actual UB overflow arithmetic).
  subroutine test_nrec_k20_overflow_class(nfail)
    integer, intent(inout) :: nfail
    integer(8) :: got
    got = expected_prod_dk_nrec(8000, 27, 100)
    call check_int8(got, 2160000000_8, "G: k20^3 nk=8000,nbvec=27,no=100 -> nrec=2,160,000,000 exactly (int64)", nfail)
    call check_true(got > int(huge(0_4), 8), &
      "G: that value EXCEEDS huge(0_4)=2**31-1 -- default INTEGER(4) cannot hold it", nfail)
  end subroutine test_nrec_k20_overflow_class

  !======================= Test H: k24-scale overflow class ================
  subroutine test_nrec_k24_overflow_class(nfail)
    integer, intent(inout) :: nfail
    integer(8) :: got
    got = expected_prod_dk_nrec(13824, 27, 150)
    call check_int8(got, 8398080000_8, "H: k24^3 nk=13824,nbvec=27,no=150 -> nrec=8,398,080,000 exactly (int64)", nfail)
    call check_true(got > int(huge(0_4), 8), &
      "H: that value EXCEEDS huge(0_4)=2**31-1 -- default INTEGER(4) cannot hold it", nfail)
  end subroutine test_nrec_k24_overflow_class

  !======================= Test I: per-rank prod_dk memory reduction ========
  ! bytes = 16 (complex(8)) * nb_eff**2 * nbvec * nk. Representative
  ! nb_eff=30, nk=8000 (k20^3): the 27->3 column reduction must cut this
  ! EXACT array by a factor of 9 (independent of nb_eff/nk, since both share
  ! the same nb_eff**2*nk factor) -- this is the "largest term" memory claim.
  subroutine test_memory_reduction_estimate(nfail)
    integer, intent(inout) :: nfail
    integer(8) :: nb_eff, nk, bytes_before, bytes_after
    integer :: nbvec_before, nbvec_after
    nb_eff = 30_8; nk = 8000_8
    nbvec_before = gicov_prod_dk_nbvec('gi',    1, 27)  ! unchanged path -> 27
    nbvec_after  = gicov_prod_dk_nbvec('gicov', 1, 27)  ! reduced path   -> 3
    bytes_before = 16_8 * nb_eff * nb_eff * int(nbvec_before, 8) * nk
    bytes_after  = 16_8 * nb_eff * nb_eff * int(nbvec_after, 8) * nk
    call check_int8(bytes_before, 3110400000_8, &
      "I: prod_dk bytes/rank BEFORE (nb_eff=30,nk=8000,nbvec=27) = 3,110,400,000 (~2.9 GiB)", nfail)
    call check_int8(bytes_after, 345600000_8, &
      "I: prod_dk bytes/rank AFTER  (nb_eff=30,nk=8000,nbvec=3)  =   345,600,000 (~0.32 GiB)", nfail)
    call check_int8(bytes_before / bytes_after, 9_8, "I: reduction factor is EXACTLY 9x (27/3), independent of nb_eff/nk", nfail)
  end subroutine test_memory_reduction_estimate

end program test_gicov_memfix
