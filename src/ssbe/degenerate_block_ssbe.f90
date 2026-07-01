! src/ssbe/degenerate_block_ssbe.f90
! LG-SBE Tier2: degenerate-block identification for the length-gauge
! semiconductor Bloch solver.
!
! Single responsibility -- "which bands form a degenerate block, and do those
! blocks stay dimensionally consistent across neighbouring k-points":
!
!   build_blocks       -- per k-point union-find over the eigen spectrum; bands
!                         whose energy separation is below theta_off are placed
!                         in the same connected component (block). Operates in
!                         the eigenbasis: band indices are NOT reordered or
!                         rotated.
!
!   same_block         -- 3-argument query same_block(ib,jb,ik) against the
!                         block assignment cached by the last build_blocks call.
!
!   match_link_blocks  -- for every k -> k+b link, use the prod_dk overlaps
!                         <u_{n,k}|u_{m,k+b}> to confirm that each block at k
!                         maps onto a single block at k+b of the SAME dimension.
!                         Links whose block dimension is not conserved (split /
!                         merge / cross-block leakage) are recorded so Pb3's xi
!                         construction can skip or diagnose them. This stage is
!                         the framework only (dimension check + fail record).
!
! Pure module (no SALMON dependencies) so it is unit-testable on its own, in the
! same spirit as sbe_collision_gw_core:
!     gfortran degenerate_block_ssbe.f90 test/test_degenerate_block.f90 -o t
!
! Thresholds (atomic units). theta_off is the block-membership (turn-"off")
! bound and is deliberately loose enough to swallow the Tier0 inter-sub-block
! gap (10.87 meV = 3.995e-4 au), so the {21,22} and {23,24,25} sub-blocks fuse
! into a single block. theta_on is the tighter core-degeneracy (turn-"on")
! bound, reserved for Pb3 dynamic hysteresis (theta_on < theta_off).
module degenerate_block_ssbe
  implicit none
  private
  public :: theta_on, theta_off
  public :: build_blocks
  public :: same_block
  public :: match_link_blocks

  real(8), parameter :: theta_on  = 5d-4   ! tight core-degeneracy bound  [au]
  real(8), parameter :: theta_off = 2d-3   ! block-membership bound       [au]

  ! Block assignment cached by the most recent build_blocks() call so that the
  ! 3-argument same_block(ib,jb,ik) query can be answered without threading
  ! block_id through every caller.
  integer, allocatable, save :: block_id_store(:, :)
  logical, save :: blocks_built = .false.

contains

  !-------------------------------------------------------------------
  ! Per-k union-find: group bands into connected components ("blocks")
  ! where any two members a,b satisfy |eigen(a)-eigen(b)| < theta_off.
  ! Connectivity is transitive, so a chain of near-degenerate states
  ! coalesces even when its end points are farther apart than theta_off.
  ! block_id(ib,ik) in 1..nblk(ik), numbered in ascending band order.
  !-------------------------------------------------------------------
  subroutine build_blocks(nb, nk, eigen, block_id)
    implicit none
    integer, intent(in)  :: nb, nk
    real(8), intent(in)  :: eigen(nb, nk)
    integer, intent(out) :: block_id(nb, nk)
    integer :: ik, ia, ib, ra, rb, r, nblk
    integer :: parent(nb)
    integer :: label(nb)

    do ik = 1, nk
      ! init: every band is its own root
      do ia = 1, nb
        parent(ia) = ia
      end do
      ! union all near-degenerate pairs (O(nb^2); nb is small)
      do ia = 1, nb
        do ib = ia + 1, nb
          if (abs(eigen(ia, ik) - eigen(ib, ik)) < theta_off) then
            ra = uf_find(parent, ia)
            rb = uf_find(parent, ib)
            if (ra /= rb) parent(rb) = ra
          end if
        end do
      end do
      ! relabel component roots to contiguous ids in ascending band order
      label(:) = 0
      nblk = 0
      do ia = 1, nb
        r = uf_find(parent, ia)
        if (label(r) == 0) then
          nblk = nblk + 1
          label(r) = nblk
        end if
        block_id(ia, ik) = label(r)
      end do
    end do

    ! cache the assignment for same_block()
    if (allocated(block_id_store)) deallocate(block_id_store)
    allocate(block_id_store(nb, nk))
    block_id_store(:, :) = block_id(:, :)
    blocks_built = .true.
  end subroutine build_blocks

  !-------------------------------------------------------------------
  ! Union-find "find" with iterative path compression (no recursion).
  !-------------------------------------------------------------------
  integer function uf_find(parent, i) result(root)
    implicit none
    integer, intent(inout) :: parent(:)
    integer, intent(in)    :: i
    integer :: j, nxt
    root = i
    do while (parent(root) /= root)
      root = parent(root)
    end do
    ! compress the path so future finds are O(1)
    j = i
    do while (parent(j) /= root)
      nxt = parent(j)
      parent(j) = root
      j = nxt
    end do
  end function uf_find

  !-------------------------------------------------------------------
  ! Are bands ib and jb in the same block at k-point ik? Uses the block
  ! assignment cached by the last build_blocks() call. Returns .false.
  ! (rather than aborting) if no build has run yet or indices are out of
  ! range, so callers can guard cheaply.
  !-------------------------------------------------------------------
  logical function same_block(ib, jb, ik) result(same)
    implicit none
    integer, intent(in) :: ib, jb, ik
    same = .false.
    if (.not. blocks_built) return
    if (.not. allocated(block_id_store)) return
    if (ik < 1 .or. ik > size(block_id_store, 2)) return
    if (ib < 1 .or. ib > size(block_id_store, 1)) return
    if (jb < 1 .or. jb > size(block_id_store, 1)) return
    same = (block_id_store(ib, ik) == block_id_store(jb, ik))
  end function same_block

  !-------------------------------------------------------------------
  ! Check block correspondence across each k -> k+b link and build the
  ! link-local member map used by Pb3's xi construction.
  !
  ! For each block B (dimension d) at k, the column weights
  !     w(m) = sum_{n in B} |<u_{n,k}|u_{m,k+b}>|^2
  ! identify the image columns (w > weight_thresh). The link is consistent
  ! for B iff those columns (i) all belong to ONE block at k+b and (ii) that
  ! target block has dimension d and (iii) there are exactly d image columns.
  ! Otherwise the block is a split/merge/leakage and is recorded (link_ok=F,
  ! n_fail incremented, link_member zeroed for its members).
  !
  ! ik_neighbor(iv,ik) gives the k index of k+b for link iv (<1 or >nk => the
  ! link is absent/boundary and is skipped). bvec is carried for diagnostics.
  !
  ! link_member(n,iv,ik) = the best-overlap target band index at k+b for band n
  ! at k within the matched target block (0 on failure). Within a degenerate
  ! block the individual band assignment is gauge-ambiguous; Pb3 uses the full
  ! block x block submatrix of prod_dk -- this greedy map is for
  ! initialisation / diagnostics only.
  !-------------------------------------------------------------------
  subroutine match_link_blocks(nb, nk, nbvec, bvec, prod_dk, block_id, ik_neighbor, &
                               link_member, link_ok, n_fail, weight_thresh)
    implicit none
    integer,    intent(in)  :: nb, nk, nbvec
    integer,    intent(in)  :: bvec(3, nbvec)            ! integer dk shifts (diagnostics)
    complex(8), intent(in)  :: prod_dk(nb, nb, nbvec, nk)
    integer,    intent(in)  :: block_id(nb, nk)
    integer,    intent(in)  :: ik_neighbor(nbvec, nk)    ! k index of k+b (<1 or >nk => skip)
    integer,    intent(out) :: link_member(nb, nbvec, nk)
    logical,    intent(out) :: link_ok(nbvec, nk)
    integer,    intent(out) :: n_fail
    real(8),    intent(in), optional :: weight_thresh
    integer :: ik, iv, ikpb, bk, nblk_k, n, m, d, dt, image_dim, tgt_block, best
    real(8) :: wthr, bestval
    real(8) :: w(nb)
    logical :: consistent

    wthr = 0.5d0
    if (present(weight_thresh)) wthr = weight_thresh

    link_member(:, :, :) = 0
    link_ok(:, :)        = .true.
    n_fail               = 0

    do ik = 1, nk
      do iv = 1, nbvec
        ikpb = ik_neighbor(iv, ik)
        if (ikpb < 1 .or. ikpb > nk) cycle    ! absent / boundary link

        nblk_k = maxval(block_id(:, ik))
        do bk = 1, nblk_k
          ! source-block dimension d
          d = 0
          do n = 1, nb
            if (block_id(n, ik) == bk) d = d + 1
          end do
          if (d == 0) cycle

          ! per-column overlap weight w(m) = sum_{n in bk} |<u_{n,k}|u_{m,k+b}>|^2
          do m = 1, nb
            w(m) = 0d0
            do n = 1, nb
              if (block_id(n, ik) == bk) w(m) = w(m) + abs(prod_dk(n, m, iv, ik))**2
            end do
          end do

          ! image detection: significant columns must all belong to one block
          ! at k+b, and there must be exactly d of them.
          tgt_block  = 0
          image_dim  = 0
          consistent = .true.
          do m = 1, nb
            if (w(m) > wthr) then
              image_dim = image_dim + 1
              if (tgt_block == 0) then
                tgt_block = block_id(m, ikpb)
              else if (block_id(m, ikpb) /= tgt_block) then
                consistent = .false.
              end if
            end if
          end do

          dt = 0
          if (tgt_block > 0) then
            do m = 1, nb
              if (block_id(m, ikpb) == tgt_block) dt = dt + 1
            end do
          end if

          if (consistent .and. tgt_block > 0 .and. image_dim == d .and. dt == d) then
            ! dimension conserved: record best-overlap band correspondence
            ! within the matched target block.
            do n = 1, nb
              if (block_id(n, ik) /= bk) cycle
              best    = 0
              bestval = -1d0
              do m = 1, nb
                if (block_id(m, ikpb) /= tgt_block) cycle
                if (abs(prod_dk(n, m, iv, ik)) > bestval) then
                  bestval = abs(prod_dk(n, m, iv, ik))
                  best    = m
                end if
              end do
              link_member(n, iv, ik) = best
            end do
          else
            ! split / merge / cross-block leakage: flag the link and this block.
            link_ok(iv, ik) = .false.
            n_fail = n_fail + 1
            do n = 1, nb
              if (block_id(n, ik) == bk) link_member(n, iv, ik) = 0
            end do
          end if
        end do
      end do
    end do
  end subroutine match_link_blocks

end module degenerate_block_ssbe
