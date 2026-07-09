! src/ssbe/test/test_degenerate_block.f90
! Unit tests for degenerate_block_ssbe (LG-SBE Tier2 block identification).
! The module is dependency-free, so this compiles standalone, e.g. from
! src/ssbe:
!     gfortran degenerate_block_ssbe.f90 test/test_degenerate_block.f90 -o t && ./t
program test_degenerate_block
  use degenerate_block_ssbe
  implicit none
  integer :: nfail

  nfail = 0

  call test_two_fold(nfail)
  call test_three_fold(nfail)
  call test_non_degenerate(nfail)
  call test_transitive_chain(nfail)
  call test_tier0_fusion(nfail)
  call test_same_block_api(nfail)
  call test_match_ok(nfail)
  call test_match_split(nfail)

  if (nfail == 0) then
    write(*,*) "ALL DEGENERATE-BLOCK TESTS PASSED"
  else
    write(*,*) "FAILED:", nfail
    stop 1
  end if

contains

  subroutine check_int(got, want, label, nfail)
    integer, intent(in)      :: got, want
    character(*), intent(in) :: label
    integer, intent(inout)   :: nfail
    if (got /= want) then
      write(*,'(a,a,2(1x,i0))') "MISMATCH ", label, got, want
      nfail = nfail + 1
    else
      write(*,'(a,a)') "ok ", label
    end if
  end subroutine check_int

  subroutine check_true(cond, label, nfail)
    logical, intent(in)      :: cond
    character(*), intent(in) :: label
    integer, intent(inout)   :: nfail
    if (.not. cond) then
      write(*,'(a,a)') "MISMATCH ", label
      nfail = nfail + 1
    else
      write(*,'(a,a)') "ok ", label
    end if
  end subroutine check_true

  ! two nearly-degenerate bands + one isolated band -> [1,1,2]
  subroutine test_two_fold(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 3, nk = 1
    real(8) :: eigen(nb, nk)
    integer :: bid(nb, nk)
    eigen(:, 1) = [0d0, 1d-6, 1d0]
    call build_blocks(nb, nk, eigen, bid)
    call check_int(bid(1,1), 1, "two_fold b1", nfail)
    call check_int(bid(2,1), 1, "two_fold b2", nfail)
    call check_int(bid(3,1), 2, "two_fold b3", nfail)
  end subroutine test_two_fold

  ! three degenerate + one isolated -> [1,1,1,2]
  subroutine test_three_fold(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 4, nk = 1
    real(8) :: eigen(nb, nk)
    integer :: bid(nb, nk)
    eigen(:, 1) = [0d0, 1d-7, 2d-7, 0.5d0]
    call build_blocks(nb, nk, eigen, bid)
    call check_int(bid(1,1), 1, "three_fold b1", nfail)
    call check_int(bid(2,1), 1, "three_fold b2", nfail)
    call check_int(bid(3,1), 1, "three_fold b3", nfail)
    call check_int(bid(4,1), 2, "three_fold b4", nfail)
  end subroutine test_three_fold

  ! all gaps > theta_off -> distinct blocks [1,2,3]
  subroutine test_non_degenerate(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 3, nk = 1
    real(8) :: eigen(nb, nk)
    integer :: bid(nb, nk)
    eigen(:, 1) = [0d0, 0.1d0, 0.2d0]
    call build_blocks(nb, nk, eigen, bid)
    call check_int(bid(1,1), 1, "nondeg b1", nfail)
    call check_int(bid(2,1), 2, "nondeg b2", nfail)
    call check_int(bid(3,1), 3, "nondeg b3", nfail)
  end subroutine test_non_degenerate

  ! transitive closure: consecutive gaps 1.5e-3 (< theta_off) but end-to-end
  ! 4.5e-3 (> theta_off) must still coalesce into ONE block.
  subroutine test_transitive_chain(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 4, nk = 1
    real(8) :: eigen(nb, nk)
    integer :: bid(nb, nk)
    eigen(:, 1) = [0d0, 1.5d-3, 3.0d-3, 4.5d-3]
    call build_blocks(nb, nk, eigen, bid)
    call check_int(bid(1,1), 1, "chain b1", nfail)
    call check_int(bid(2,1), 1, "chain b2", nfail)
    call check_int(bid(3,1), 1, "chain b3", nfail)
    call check_int(bid(4,1), 1, "chain b4", nfail)
    call check_int(maxval(bid(:,1)), 1, "chain single block", nfail)
  end subroutine test_transitive_chain

  ! Tier0: bands {21,22} degenerate and {23,24,25} degenerate, separated by
  ! 10.87 meV = 3.995e-4 au < theta_off -> a single 5-band block.
  subroutine test_tier0_fusion(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 5, nk = 1
    real(8), parameter :: gap = 3.995d-4          ! 10.87 meV in atomic units
    real(8) :: eigen(nb, nk)
    integer :: bid(nb, nk)
    integer :: n
    ! module indices 1..5 <-> physical bands 21..25
    eigen(1, 1) = 0d0
    eigen(2, 1) = 1d-6
    eigen(3, 1) = gap
    eigen(4, 1) = gap + 1d-6
    eigen(5, 1) = gap + 2d-6
    call build_blocks(nb, nk, eigen, bid)
    do n = 1, nb
      call check_int(bid(n,1), 1, "tier0 member", nfail)
    end do
    call check_int(maxval(bid(:,1)), 1, "tier0 single block", nfail)
  end subroutine test_tier0_fusion

  ! same_block API after a build.
  subroutine test_same_block_api(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 3, nk = 1
    real(8) :: eigen(nb, nk)
    integer :: bid(nb, nk)
    eigen(:, 1) = [0d0, 1d-6, 1d0]
    call build_blocks(nb, nk, eigen, bid)
    call check_true(same_block(1, 2, 1),       "same_block deg pair",   nfail)
    call check_true(.not. same_block(1, 3, 1), "same_block across gap", nfail)
  end subroutine test_same_block_api

  ! match_link_blocks: clean block-diagonal link -> ok, n_fail=0.
  subroutine test_match_ok(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 3, nk = 2, nbvec = 1
    real(8)    :: eigen(nb, nk)
    integer    :: bid(nb, nk)
    integer    :: ik_neighbor(nbvec, nk)
    integer    :: link_member(nb, nbvec, nk)
    logical    :: link_ok(nbvec, nk)
    integer    :: n_fail_link
    complex(8) :: prod_dk(nb, nb, nbvec, nk)
    real(8)    :: c, s
    ! blocks {1,2}=block1, {3}=block2 at both k-points
    eigen(:, 1) = [0d0, 1d-6, 1d0]
    eigen(:, 2) = [0d0, 1d-6, 1d0]
    call build_blocks(nb, nk, eigen, bid)
    ! k1 -> k2 link only; k2 has no forward neighbour (0 => skip)
    ik_neighbor(1, 1) = 2
    ik_neighbor(1, 2) = 0
    ! block-diagonal overlap: 2x2 rotation inside {1,2}, identity on {3}
    c = cos(0.3d0); s = sin(0.3d0)
    prod_dk = dcmplx(0d0, 0d0)
    prod_dk(1, 1, 1, 1) = dcmplx( c, 0d0)
    prod_dk(1, 2, 1, 1) = dcmplx(-s, 0d0)
    prod_dk(2, 1, 1, 1) = dcmplx( s, 0d0)
    prod_dk(2, 2, 1, 1) = dcmplx( c, 0d0)
    prod_dk(3, 3, 1, 1) = dcmplx(1d0, 0d0)
    call match_link_blocks(nb, nk, nbvec, reshape([0,0,0],[3,1]), prod_dk, bid, &
                           ik_neighbor, link_member, link_ok, n_fail_link)
    call check_true(link_ok(1,1), "match_ok link ok",  nfail)
    call check_int(n_fail_link, 0, "match_ok n_fail=0", nfail)
  end subroutine test_match_ok

  ! match_link_blocks: block-1 leaks into block-2 across the link -> split
  ! detected (link fails, n_fail>=1).
  subroutine test_match_split(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 3, nk = 2, nbvec = 1
    real(8)    :: eigen(nb, nk)
    integer    :: bid(nb, nk)
    integer    :: ik_neighbor(nbvec, nk)
    integer    :: link_member(nb, nbvec, nk)
    logical    :: link_ok(nbvec, nk)
    integer    :: n_fail_link
    complex(8) :: prod_dk(nb, nb, nbvec, nk)
    eigen(:, 1) = [0d0, 1d-6, 1d0]
    eigen(:, 2) = [0d0, 1d-6, 1d0]
    call build_blocks(nb, nk, eigen, bid)
    ik_neighbor(1, 1) = 2
    ik_neighbor(1, 2) = 0
    ! permutation crossing the block boundary:
    !   band1(k)->col1 (blk1), band2(k)->col3 (blk2), band3(k)->col2 (blk1)
    prod_dk = dcmplx(0d0, 0d0)
    prod_dk(1, 1, 1, 1) = dcmplx(1d0, 0d0)
    prod_dk(2, 3, 1, 1) = dcmplx(1d0, 0d0)
    prod_dk(3, 2, 1, 1) = dcmplx(1d0, 0d0)
    call match_link_blocks(nb, nk, nbvec, reshape([0,0,0],[3,1]), prod_dk, bid, &
                           ik_neighbor, link_member, link_ok, n_fail_link)
    call check_true(.not. link_ok(1,1),  "match_split link fails", nfail)
    call check_true(n_fail_link >= 1,    "match_split n_fail>=1",  nfail)
  end subroutine test_match_split

end program test_degenerate_block
