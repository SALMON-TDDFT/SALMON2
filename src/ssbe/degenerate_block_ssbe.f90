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
  use sbe_lg_mode_ssbe, only: uses_xfull_links
  implicit none
  private
  public :: theta_on, theta_off
  public :: build_blocks
  public :: build_fixed_blocks
  public :: build_blocks_fixed
  public :: same_block
  public :: match_link_blocks
  ! Pb3 (non-Abelian xi + smooth blend):
  public :: blend                 ! cos^2 hysteresis ramp 1 -> 0
  public :: build_ik_neighbor     ! k -> k+b index table on the full uniform grid
  public :: find_bvec             ! locate the bvec column of an integer dk shift (0 = absent)
  public :: identity_kmap         ! kmap of the REPLICATED layout: kmap(ik) = ik (serial/full-nk callers)
  public :: xi_block_from_overlap ! polar + matrix-log kernel (LAPACK zheev/zgeev)
  public :: build_xi              ! orchestrator: blocks -> links -> xi(nb,nb,3,nk)
  public :: xi_sign               ! sign s in xi = s*i*logm(U)/dk (pinned by test)
  ! Phase 1 (gicov / Approach-B'): polar-only Wilson-line transport, no logm,
  ! so no eigenphase branch cut can occur (removes build_xi's info=2 reject
  ! mode by construction). Standalone kernel; consumed by later phases.
  public :: polar_unitary         ! U = M(M^H M)^{-1/2} polar factor (LAPACK zheev only)
  public :: build_block_transport ! orchestrator: fixed blocks -> block-diagonal U_full(nb,nb,3,nk)
  public :: build_fixed_blocks_closed
  public :: build_blocks_fixed_closed
  ! gicov memory-scale fix (LG-SBE Phase 1): prod_dk column-reduction helpers,
  ! consumed by gs_info_ssbe's read_prod_dk_data (see the header comment right
  ! before prod_dk_axis_slot below for why they live in this pure module).
  public :: prod_dk_axis_slot     ! which of the <=3 kept gicov columns a (jdk1,jdk2,jdk3) shift is, or 0
  public :: gicov_prod_dk_nbvec   ! stored gs%prod_dk 4th-dim extent (3 for gicov/ndk>=1, else nbvec_full)
  public :: expected_prod_dk_nrec ! int64-safe file_sbe_prod_dk record count nk*(2ndk+1)**3*no**2
  ! Phase 2 (gicov): pure gauge-covariant intraband k-derivative operator that
  ! consumes build_block_transport's U_full (Wilson-line transport, no logm).
  public :: covariant_grad_block  ! D_cov rho = d_k rho - i[xi,rho] (<=4-shell transported stencil, per-axis alias-capped)
  ! Halo support (pure, no MPI; consumed by bloch_solver's gicov halo exchange):
  public :: covariant_halo_needed ! mark every k whose rho the stencil READS for a local slice
  public :: build_halo_lists      ! derive deterministic per-rank send/recv k-lists from needed-sets
  public :: check_gicov_occupation ! metal detector: fractional / k-varying occupation

  real(8), parameter :: theta_on  = 5d-4   ! tight core-degeneracy bound  [au]
  real(8), parameter :: theta_off = 2d-3   ! block-membership bound       [au]

  ! Pb3 numerical conventions (see xi_block_from_overlap / build_xi):
  !   xi = s * i * logm(U) / dk_alpha, and the finite-gap limit must reproduce
  !   SALMON's dipole d = i * p_mod / delta_omega. The 2-band smooth-rotation
  !   model in test/test_xi_construction.f90 shows i*logm(U)/dk = i<u_n|d_k u_m>
  !   whereas d = -i<u_n|d_k u_m>, hence s = -1 (verified numerically there).
  real(8), parameter :: xi_sign       = -1d0   ! sign s (pinned by the unit test)
  real(8), parameter :: xi_sing_tol   = 1d-6   ! reject link if min singular value(M_blk) < this
  real(8), parameter :: xi_phi_reject = 0.9d0  ! reject link if any |phi_i| > xi_phi_reject*pi
  real(8), parameter :: pi_local = 3.14159265358979323846d0
  real(8), parameter :: wleak_default = 1d-1   ! leaked-weight union threshold |<u|u'>|^2, gated behind
                                                ! polar_unitary ierr==1 rank-deficiency (gicov X-closed)

  ! Block assignment cached by the most recent build_blocks() call so that the
  ! 3-argument same_block(ib,jb,ik) query can be answered without threading
  ! block_id through every caller.
  integer, allocatable, save :: block_id_store(:, :)
  logical, save :: blocks_built = .false.

  ! k-index neighbour tables of the covariant stencil, cached across calls (see
  ! get_ik_neighbor_cache).  Pure integer index maps derived from
  ! (num_kgrid, bvec) -- caching them is bit-neutral and removes an
  ! O(nbvec*nk) allocate+fill that was paid 4x per RK4 step.
  integer, allocatable, save :: ikn_store(:, :)   ! (nbvec, nk) forward k + bvec(:,iv)
  integer, allocatable, save :: bwd_store(:, :)   ! (nk, 3)    backward, per Cartesian axis
  integer, allocatable, save :: ikn_bvec(:, :)    ! (3, nbvec) cache key
  integer, save :: ikn_kgrid(3) = 0
  integer, save :: ikn_nk = 0, ikn_nbvec = 0
  logical, save :: ikn_built = .false.

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
  ! k-INDEPENDENT composite-block partition, determined once from the full
  ! eigen(nb,nk) spectrum: connected components of the graph whose edges are
  ! band pairs near-degenerate (< theta_off) at ANY k (transitive union-find,
  ! same threshold as build_blocks' per-k grouping, but unioned across all k
  ! so a pair that is close only at one k still fixes the two bands into one
  ! block for every k). This is the algorithmic core of the Tier2 gifix mode:
  ! blocks so defined are dimensionally stable across k by construction, so
  ! Pb3's per-link dimension-mismatch failure mode cannot occur.
  !
  ! isolated_ok is a COARSE energy pre-filter, NOT the primary metal guard --
  ! it is deliberately loose so it does not false-abort real insulator /
  ! semiconductor spectra that happen to carry a small but genuine cross-block
  ! gap. The physically-grounded guard is build_xi's per-block RUNTIME
  ! fail-closed check (xi_block_from_overlap info=1/2/3), which tests the
  ! actual computed overlaps at every link rather than a static energy
  ! heuristic; a spectrum that slips past this pre-filter but is genuinely
  ! metal-like is still caught there.
  !
  ! isolated_ok verifies the partition is physically meaningful: at every k, a
  ! MULTI-band block's members must stay farther than gap_margin from every
  ! band outside the block. gap_margin is deliberately LARGER than theta_off
  ! (2d-3 au, the grouping bound) -- otherwise a band that sits just outside
  ! theta_off (so it never joins the block) could still be arbitrarily close
  ! to it without ever tripping the isolation check, making the check
  ! tautological. isolated_ok=.false. signals the ambiguous / metal-like case
  ! (a would-be outside band drifting within gap_margin of a multi-band block
  ! at some k) and the caller fails closed rather than silently mis-treating
  ! it as an isolated insulating block. Singleton-vs-singleton proximity does
  ! NOT trip isolation (two merely-close-but-distinct non-degenerate bands are
  ! not an ambiguity this routine needs to flag).
  !
  ! gap_margin calibration (2026-07, real Si 4^3 pre-flight): the previous
  ! value 1d-2 au FALSE-ABORTED on the real Si 4^3 KS spectrum (96 isolation
  ! violations) even though that run is an ordinary semiconductor. The
  ! measured global minimum cross-block gap in that spectrum is 2.654e-3 au,
  ! at the block{17..22}/{23..28} boundary (bands 22/23); this gap is
  ! k-INDEPENDENT (symmetry-like -- reproduced at ik=2,15,17,45,53,62, i.e.
  ! not eigensolver noise). gap_margin=2.2d-3 keeps ~10% headroom above
  ! theta_off=2d-3 (so the check stays non-tautological) while sitting below
  ! the real flip point 2.654e-3, so this genuine semiconductor spectrum now
  ! passes.
  !-------------------------------------------------------------------
  subroutine build_fixed_blocks(nb, nk, eigen, fixed_block_id, isolated_ok)
    implicit none
    integer, intent(in)  :: nb, nk
    real(8), intent(in)  :: eigen(nb, nk)
    integer, intent(out) :: fixed_block_id(nb)
    logical, intent(out) :: isolated_ok
    real(8), parameter :: gap_margin = 2.2d-3 ! au; MUST be > theta_off so a band can be
                                               ! outside a block (>theta_off) yet fail isolation
                                               ! (<gap_margin). Calibrated to the real Si 4^3
                                               ! spectrum's global min cross-block gap 2.654e-3 au
                                               ! (k-independent, see header above); this is a
                                               ! COARSE pre-filter only -- build_xi's runtime
                                               ! overlap check (xi_block_from_overlap) is the
                                               ! physically-grounded guard.
    integer :: parent(nb), ia, ib, ik, root, nlab, lab(nb), ra, rb, bsize(nb)

    do ia = 1, nb
      parent(ia) = ia
    end do

    ! edge if bands ia,ib are within theta_off at ANY k (k-independent grouping)
    do ik = 1, nk
      do ia = 1, nb - 1
        do ib = ia + 1, nb
          if (abs(eigen(ia, ik) - eigen(ib, ik)) < theta_off) then
            ra = uf_find(parent, ia)
            rb = uf_find(parent, ib)
            if (ra /= rb) parent(rb) = ra
          end if
        end do
      end do
    end do

    ! contiguous ascending labels
    lab = 0; nlab = 0
    do ia = 1, nb
      root = uf_find(parent, ia)
      if (lab(root) == 0) then
        nlab = nlab + 1
        lab(root) = nlab
      end if
      fixed_block_id(ia) = lab(root)
    end do

    ! gap-isolation: only for MULTI-band blocks (two close singleton bands must
    ! not false-abort). At every k, a MULTI-band block's bands must be strictly
    ! farther than gap_margin from any band outside the block.
    bsize = 0
    do ia = 1, nb
      bsize(fixed_block_id(ia)) = bsize(fixed_block_id(ia)) + 1
    end do

    isolated_ok = .true.
    do ik = 1, nk
      do ia = 1, nb
        do ib = 1, nb
          if (fixed_block_id(ia) /= fixed_block_id(ib)) then
            if (bsize(fixed_block_id(ia)) > 1 .or. bsize(fixed_block_id(ib)) > 1) then
              if (abs(eigen(ia, ik) - eigen(ib, ik)) < gap_margin) isolated_ok = .false.
            end if
          end if
        end do
      end do
    end do
  end subroutine build_fixed_blocks

  !-------------------------------------------------------------------
  ! sbe_lg_degen='gifix' entry point: derive the k-independent composite
  ! partition (build_fixed_blocks) and broadcast it into the per-k
  ! block_id(nb,nk) shape build_xi/same_block expect, so the fixed-block
  ! mode is a drop-in replacement for build_blocks' per-k union-find.
  ! Fails closed (error stop) if the partition is not gap-isolated
  ! (metal-like disentanglement is out of scope for gifix -- use
  ! velocity_gauge instead).
  !-------------------------------------------------------------------
  subroutine build_blocks_fixed(nb, nk, eigen, block_id)
    implicit none
    integer, intent(in)  :: nb, nk
    real(8), intent(in)  :: eigen(nb, nk)
    integer, intent(out) :: block_id(nb, nk)
    integer :: fixed_block_id(nb), ik
    logical :: isolated_ok

    call build_fixed_blocks(nb, nk, eigen, fixed_block_id, isolated_ok)
    if (.not. isolated_ok) then
      write(*,*) 'ERROR sbe_lg_degen=gifix: a composite block is not gap-isolated ', &
                 '(metal-like disentanglement) - out of scope, use velocity_gauge'
      error stop 1
    end if

    do ik = 1, nk
      block_id(:, ik) = fixed_block_id(:)
    end do

    ! cache for same_block (mirror build_blocks:113-118, deallocate-first)
    if (allocated(block_id_store)) deallocate(block_id_store)
    allocate(block_id_store(nb, nk))
    block_id_store(:, :) = block_id(:, :)
    blocks_built = .true.
  end subroutine build_blocks_fixed

  !-------------------------------------------------------------------
  ! gicov (Option X-closed): k-independent composite partition that is
  ! closed under BOTH near-degeneracy (theta_off, as build_fixed_blocks)
  ! AND cross-link overlap. A degenerate block that is not closed across a
  ! BZ-boundary (umklapp) wrap link -- band identity scrambles, so its
  ! overlap sub-block is rank-deficient and polar_unitary would reject --
  ! is grown to include the bands its overlap image leaks onto, until the
  ! partition is closed. Then build_block_transport's polar factor is
  ! taken on a well-conditioned (closed) block.
  !
  ! Step 1: energy union-find (identical edges to build_fixed_blocks).
  ! Step 2: RANK-DEFICIENCY-GATED overlap closure, iterated to a fixed
  !   point over the +axis links. At each sweep, snapshot the current
  !   block labels; for every block B and every +axis link at every k,
  !   extract the |B|x|B| overlap sub-block M_B(a,b) = prod_dk(mem(a),
  !   mem(b),iv,ik) and test it with the EXISTING polar_unitary(M_B,d,U,
  !   sigma_min,ierr) -- the SAME sigma_min < xi_sing_tol gate
  !   build_block_transport itself fails a block on. Ordinary band
  !   rotation on a coarse k-mesh leaves a block's self-overlap healthy
  !   (sigma_min ~0.6-0.9: ierr==0) even though individual cross terms can
  !   be sizeable (~0.1-0.3) -- a scalar cross-overlap threshold cannot
  !   tell that apart from genuine umklapp scramble, which collapses
  !   sigma_min by ~9 orders of magnitude (~1e-9). So closure is triggered
  !   ONLY when polar_unitary reports ierr==1 (rank-deficient); ierr==3
  !   (LAPACK failure) fail-closes the whole build via error stop, exactly
  !   matching build_block_transport's ierr/=0 rejection. Only once a
  !   block IS singular do we union it with the outside bands its weight
  !   leaked onto (wl(m) = sum_{a in B} |prod_dk(mem(a),m,iv,ik)|^2 >
  !   wleak) -- wleak can never cascade a healthy block because it is
  !   gated behind the singularity test. Bounded: <= nb-1 unions total.
  ! Step 3: relabel + the same gap-isolation pre-filter as build_fixed_
  !   blocks (informational; the real fail-closed guard remains
  !   polar_unitary's near-singular check at transport-build time).
  !-------------------------------------------------------------------
  subroutine build_fixed_blocks_closed(nb, nk, nbvec, bvec, prod_dk, eigen, num_kgrid, &
                                       fixed_block_id, isolated_ok, wleak_in)
    implicit none
    integer,    intent(in)  :: nb, nk, nbvec
    integer,    intent(in)  :: bvec(3, nbvec), num_kgrid(3)
    complex(8), intent(in)  :: prod_dk(nb, nb, nbvec, nk)
    real(8),    intent(in)  :: eigen(nb, nk)
    integer,    intent(out) :: fixed_block_id(nb)
    logical,    intent(out) :: isolated_ok
    real(8),    intent(in), optional :: wleak_in
    real(8), parameter :: gap_margin = 2.2d-3
    integer :: parent(nb), label(nb), lab(nb), bsize(nb)
    integer :: iv_axis(3), ia, ib, ik, axis, iv, n, m, ra, rb, root, nlab
    integer :: d, a, bcol, mem(nb), ierr
    real(8) :: wleak, sig, wl
    complex(8), allocatable :: M_B(:,:), Udum(:,:)
    logical :: changed

    wleak = wleak_default
    if (present(wleak_in)) wleak = wleak_in

    ! ---- Step 1: energy union-find (edge if within theta_off at ANY k) ----
    do ia = 1, nb
      parent(ia) = ia
    end do
    do ik = 1, nk
      do ia = 1, nb - 1
        do ib = ia + 1, nb
          if (abs(eigen(ia, ik) - eigen(ib, ik)) < theta_off) then
            ra = uf_find(parent, ia); rb = uf_find(parent, ib)
            if (ra /= rb) parent(rb) = ra
          end if
        end do
      end do
    end do

    ! ---- Step 2: rank-deficiency-gated overlap closure, iterate to fixed point ----
    iv_axis(1) = find_bvec(bvec, nbvec, 1, 0, 0)
    iv_axis(2) = find_bvec(bvec, nbvec, 0, 1, 0)
    iv_axis(3) = find_bvec(bvec, nbvec, 0, 0, 1)

    changed = .true.
    do while (changed)
      changed = .false.
      ! snapshot current block label (root band) for every band
      do ia = 1, nb
        label(ia) = uf_find(parent, ia)
      end do
      do ik = 1, nk
        do axis = 1, 3
          iv = iv_axis(axis)
          if (iv == 0) cycle
          ! process each block once, keyed by its snapshot root band
          do root = 1, nb
            if (label(root) /= root) cycle
            ! snapshot block members
            d = 0
            do n = 1, nb
              if (label(n) == root) then; d = d + 1; mem(d) = n; end if
            end do
            ! overlap sub-block M_B and its smallest singular value (via polar_unitary)
            if (allocated(M_B)) deallocate(M_B, Udum)
            allocate(M_B(d, d), Udum(d, d))
            do bcol = 1, d
              do a = 1, d
                M_B(a, bcol) = prod_dk(mem(a), mem(bcol), iv, ik)
              end do
            end do
            call polar_unitary(M_B, d, Udum, sig, ierr)
            deallocate(M_B, Udum)
            if (ierr == 3) then            ! LAPACK failure -> fail-closed (build_block_transport rejects ierr/=0)
              write(*,'(a,i0,a,i0)') 'ERROR build_fixed_blocks_closed: polar_unitary LAPACK failure ik=', ik, ' axis=', axis
              error stop 1
            end if
            if (ierr /= 1) cycle           ! ierr==1 <=> sigma_min < xi_sing_tol: the ONLY close trigger
            ! block is genuinely rank-deficient on this link -> union with the
            ! outside bands its weight leaked onto
            do m = 1, nb
              if (label(m) == root) cycle
              wl = 0d0
              do a = 1, d
                wl = wl + abs(prod_dk(mem(a), m, iv, ik))**2
              end do
              if (wl > wleak) then
                ra = uf_find(parent, root); rb = uf_find(parent, m)
                if (ra /= rb) then; parent(rb) = ra; changed = .true.; end if
              end if
            end do
          end do
        end do
      end do
    end do

    ! ---- Step 3: contiguous ascending labels ----
    lab = 0; nlab = 0
    do ia = 1, nb
      root = uf_find(parent, ia)
      if (lab(root) == 0) then
        nlab = nlab + 1
        lab(root) = nlab
      end if
      fixed_block_id(ia) = lab(root)
    end do

    ! ---- gap-isolation pre-filter on the FINAL blocks (multi-band only) ----
    bsize = 0
    do ia = 1, nb
      bsize(fixed_block_id(ia)) = bsize(fixed_block_id(ia)) + 1
    end do
    isolated_ok = .true.
    do ik = 1, nk
      do ia = 1, nb
        do ib = 1, nb
          if (fixed_block_id(ia) /= fixed_block_id(ib)) then
            if (bsize(fixed_block_id(ia)) > 1 .or. bsize(fixed_block_id(ib)) > 1) then
              if (abs(eigen(ia, ik) - eigen(ib, ik)) < gap_margin) isolated_ok = .false.
            end if
          end if
        end do
      end do
    end do
  end subroutine build_fixed_blocks_closed

  !-------------------------------------------------------------------
  ! gicov single source of truth for the block partition: overlap-closed
  ! fixed blocks broadcast into per-k block_id AND cached for same_block
  ! (mirrors build_blocks_fixed, but with the overlap-closure builder).
  !
  ! The gap-isolation result is INFORMATIONAL for gicov X-closed (warn,
  ! do NOT error-stop): the overlap-closed block is closed by construction,
  ! and the REAL fail-closed guards are polar_unitary's near-singular check
  ! (in build_block_transport) + check_gicov_occupation (fractional/unequal
  ! occupation = metal). A hard abort here would false-abort a legitimately
  ! closed block that merely sits within gap_margin of a non-leaking outside
  ! band. The closure-size diagnostic (gate G0) is emitted by the CALLER
  ! (gs_info, under irank==0), since this pure module has no MPI rank info.
  !-------------------------------------------------------------------
  subroutine build_blocks_fixed_closed(nb, nk, nbvec, bvec, prod_dk, eigen, num_kgrid, block_id)
    implicit none
    integer,    intent(in)  :: nb, nk, nbvec
    integer,    intent(in)  :: bvec(3, nbvec), num_kgrid(3)
    complex(8), intent(in)  :: prod_dk(nb, nb, nbvec, nk)
    real(8),    intent(in)  :: eigen(nb, nk)
    integer,    intent(out) :: block_id(nb, nk)
    integer :: fixed_block_id(nb), ik
    logical :: isolated_ok

    call build_fixed_blocks_closed(nb, nk, nbvec, bvec, prod_dk, eigen, num_kgrid, &
                                   fixed_block_id, isolated_ok)
    if (.not. isolated_ok) then
      ! informational only (see header): the real guards run downstream.
      write(*,*) 'WARNING sbe_lg_degen=gicov: an overlap-closed block sits within ', &
                 'gap_margin of an outside band (informational; polar_unitary + ', &
                 'occupation guard are the real fail-closed checks)'
    end if

    do ik = 1, nk
      block_id(:, ik) = fixed_block_id(:)
    end do

    ! cache for same_block (mirror build_blocks_fixed: deallocate-first)
    if (allocated(block_id_store)) deallocate(block_id_store)
    allocate(block_id_store(nb, nk))
    block_id_store(:, :) = block_id(:, :)
    blocks_built = .true.
  end subroutine build_blocks_fixed_closed

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

  !===================================================================
  ! Pb3 : non-Abelian Berry connection xi + smooth blend
  !===================================================================

  !-------------------------------------------------------------------
  ! cos^2 hysteresis ramp, w(x): 1 for x<=xon, 0 for x>=xoff, C^1-smooth
  ! in between (zero slope at both ends, no discontinuity at the thetas).
  !-------------------------------------------------------------------
  real(8) function blend(x, xon, xoff) result(w)
    implicit none
    real(8), intent(in) :: x, xon, xoff
    real(8) :: t, c
    if (xoff <= xon) then                 ! degenerate window -> hard step
      if (x <= xon) then; w = 1d0; else; w = 0d0; end if
      return
    end if
    if (x <= xon) then
      w = 1d0
    else if (x >= xoff) then
      w = 0d0
    else
      t = (x - xon) / (xoff - xon)        ! 0 -> 1
      c = cos(0.5d0 * pi_local * t)       ! 1 -> 0
      w = c * c                           ! cos^2 : w(xon)=1, w(xoff)=0, C^1
    end if
  end function blend

  !-------------------------------------------------------------------
  ! k -> k+b index table for the FULL uniform grid nk = product(num_kgrid).
  ! Ordering matches common_ssbe / write_prod_dk:
  !     ik = i3*N2*N1 + i2*N1 + i1 + 1   (i1 fastest, 0-based i1,i2,i3).
  ! Neighbours are periodic (modulo), so every link exists on the full grid.
  !-------------------------------------------------------------------
  subroutine build_ik_neighbor(num_kgrid, bvec, nbvec, nk, ik_neighbor)
    implicit none
    integer, intent(in)  :: num_kgrid(3), nbvec, nk
    integer, intent(in)  :: bvec(3, nbvec)
    integer, intent(out) :: ik_neighbor(nbvec, nk)
    integer :: ik, iv, i1, i2, i3, j1, j2, j3, n1, n2, n3
    n1 = num_kgrid(1); n2 = num_kgrid(2); n3 = num_kgrid(3)
    do ik = 1, nk
      i1 = mod(ik - 1, n1)
      i2 = mod((ik - 1) / n1, n2)
      i3 = (ik - 1) / (n1 * n2)
      do iv = 1, nbvec
        j1 = modulo(i1 + bvec(1, iv), n1)
        j2 = modulo(i2 + bvec(2, iv), n2)
        j3 = modulo(i3 + bvec(3, iv), n3)
        ik_neighbor(iv, ik) = j3 * n1 * n2 + j2 * n1 + j1 + 1
      end do
    end do
  end subroutine build_ik_neighbor

  !-------------------------------------------------------------------
  ! Locate the bvec column equal to the integer shift (j1,j2,j3); 0 if absent.
  !-------------------------------------------------------------------
  integer function find_bvec(bvec, nbvec, j1, j2, j3) result(iv)
    implicit none
    integer, intent(in) :: nbvec, bvec(3, nbvec), j1, j2, j3
    integer :: k
    iv = 0
    do k = 1, nbvec
      if (bvec(1, k) == j1 .and. bvec(2, k) == j2 .and. bvec(3, k) == j3) then
        iv = k; return
      end if
    end do
  end function find_bvec

  !-------------------------------------------------------------------
  ! Non-Abelian Berry connection block from a k -> k+b overlap submatrix.
  !
  !   M_blk(a,b) = <u_{srcs(a),k} | u_{tgts(b),k+b}>   (d x d, both ascending)
  !   G  = M_blk^H M_blk                    (Hermitian, PD for a good link)
  !   U  = M_blk G^{-1/2}                    (polar factor: closest unitary)
  !   L  = logm(U) = V diag(i*phi) V^H       (phi in (-pi,pi], V unitary)
  !   xi_blk = s_sign * i * L / dk_alpha     (Hermitian)
  !
  ! info = 0 ok
  !      = 1 min singular value(M_blk) < xi_sing_tol           -> reject
  !      = 2 an eigenphase within xi_phi_reject*pi of +-pi     -> branch risk
  !      = 3 LAPACK failure
  ! resid_unitary (opt) = max|U^H U - I|, health check on the polar step.
  !
  ! LAPACK: zheev (Hermitian eig of G), zgeev (eig of the unitary U). No SALMON
  ! deps; the standalone test links with -llapack -lblas.
  !-------------------------------------------------------------------
  subroutine xi_block_from_overlap(M_blk, dk_alpha, s_sign, xi_blk, info, resid_unitary)
    implicit none
    complex(8), intent(in)  :: M_blk(:, :)
    real(8),    intent(in)  :: dk_alpha
    real(8),    intent(in)  :: s_sign
    complex(8), intent(out) :: xi_blk(:, :)
    integer,    intent(out) :: info
    real(8),    intent(out), optional :: resid_unitary
    integer :: d, i, j, k, ierr, lwork
    real(8) :: gmin, nrm
    complex(8) :: proj
    complex(8), allocatable :: A(:,:), Adg(:,:), Ginvsqrt(:,:), U(:,:), Uc(:,:)
    complex(8), allocatable :: VR(:,:), V(:,:), Vdg(:,:), Lmat(:,:), Idn(:,:)
    complex(8), allocatable :: lam(:), cwork(:), vldum(:,:)
    real(8),    allocatable :: w_eig(:), rwork(:), phi(:)
    complex(8), parameter :: zi = (0d0, 1d0)

    d    = size(M_blk, 1)
    info = 0
    xi_blk = (0d0, 0d0)

    allocate(A(d,d), Adg(d,d), Ginvsqrt(d,d), U(d,d), Uc(d,d))
    allocate(VR(d,d), V(d,d), Vdg(d,d), Lmat(d,d), Idn(d,d))
    allocate(lam(d), w_eig(d), phi(d))

    ! --- G = M^H M ; Hermitian eigendecomposition (zheev): A->vecs, w_eig ascending
    A = matmul(conjg(transpose(M_blk)), M_blk)
    lwork = max(1, 2*d - 1)
    allocate(cwork(lwork), rwork(max(1, 3*d - 2)))
    call zheev('V', 'U', d, A, d, w_eig, cwork, lwork, rwork, ierr)
    deallocate(cwork, rwork)
    if (ierr /= 0) then
      info = 3; call cleanup(); return
    end if

    ! --- near-singular reject: sigma_min(M) = sqrt(min eig of G) ---
    gmin = w_eig(1)
    if (gmin < xi_sing_tol * xi_sing_tol) then
      info = 1; call cleanup(); return
    end if

    ! --- G^{-1/2} = A diag(1/sqrt(w)) A^H ;  U = M G^{-1/2} ---
    do k = 1, d
      Adg(:, k) = A(:, k) / sqrt(w_eig(k))
    end do
    Ginvsqrt = matmul(Adg, conjg(transpose(A)))
    U = matmul(M_blk, Ginvsqrt)

    if (present(resid_unitary)) then
      Idn = (0d0, 0d0)
      do i = 1, d
        Idn(i, i) = (1d0, 0d0)
      end do
      resid_unitary = maxval(abs(matmul(conjg(transpose(U)), U) - Idn))
    end if

    ! --- eig of the unitary U (zgeev): lam(i) = e^{i phi_i}, VR = right vectors ---
    Uc = U
    lwork = max(1, 4*d)
    allocate(cwork(lwork), rwork(max(1, 2*d)), vldum(1,1))
    call zgeev('N', 'V', d, Uc, d, lam, vldum, 1, VR, d, cwork, lwork, rwork, ierr)
    deallocate(cwork, rwork, vldum)
    if (ierr /= 0) then
      info = 3; call cleanup(); return
    end if

    ! --- principal phases phi in (-pi,pi]; reject near +-pi (branch cut) ---
    do i = 1, d
      phi(i) = atan2(aimag(lam(i)), dble(lam(i)))
    end do
    if (maxval(abs(phi)) > xi_phi_reject * pi_local) then
      info = 2; call cleanup(); return
    end if

    ! --- U is normal (unitary): re-orthonormalise VR (modified Gram-Schmidt) so V
    !     is unitary and L = V diag(i phi) V^H is exact. Within a degenerate phase
    !     the basis choice inside the eigenspace does not change L.
    V = VR
    do i = 1, d
      do j = 1, i - 1
        proj = dot_product(V(:, j), V(:, i))   ! <V_j|V_i> (dot_product conjugates arg1)
        V(:, i) = V(:, i) - proj * V(:, j)
      end do
      nrm = sqrt(dble(dot_product(V(:, i), V(:, i))))
      V(:, i) = V(:, i) / nrm
    end do

    ! --- L = logm(U) = V diag(i phi) V^H (anti-Hermitian) ---
    do k = 1, d
      Vdg(:, k) = V(:, k) * (zi * phi(k))
    end do
    Lmat = matmul(Vdg, conjg(transpose(V)))

    ! --- xi = s * i * logm(U) / dk_alpha  (Hermitian) ---
    xi_blk = (s_sign * zi / dk_alpha) * Lmat

    call cleanup()

  contains
    subroutine cleanup()
      if (allocated(A)) deallocate(A, Adg, Ginvsqrt, U, Uc, VR, V, Vdg, Lmat, Idn, &
                                   lam, w_eig, phi)
    end subroutine cleanup
  end subroutine xi_block_from_overlap

  !-------------------------------------------------------------------
  ! Build xi(nb,nb,3,nk): for each k, each degenerate block (d>=2) and each
  ! Cartesian axis, take the +unit-shift overlap submatrix, polar-decompose and
  ! matrix-log it (xi_block_from_overlap). xi_ok(ib,jb,ik) flags pairs where the
  ! non-Abelian connection is trustworthy (all 3 axes succeeded); elsewhere xi=0
  ! and the caller falls back to i*p/delta_omega.
  !
  ! Axis alpha uses the reciprocal spacing dk_alpha = b_matrix(alpha,alpha)/N_alpha
  ! and the +e_alpha shift, matching the diagonal-reciprocal-lattice convention
  ! already used by the length-gauge k-gradient in common_ssbe.
  !
  ! Also calls build_blocks (so same_block() is cached for the caller's blend).
  !-------------------------------------------------------------------
  subroutine build_xi(nb, nk, nbvec, bvec, prod_dk, eigen, b_matrix, num_kgrid, &
                      xi, xi_ok, n_reject, resid_max, fixed_blocks)
    implicit none
    integer,    intent(in)  :: nb, nk, nbvec
    integer,    intent(in)  :: bvec(3, nbvec)
    complex(8), intent(in)  :: prod_dk(nb, nb, nbvec, nk)
    real(8),    intent(in)  :: eigen(nb, nk)
    real(8),    intent(in)  :: b_matrix(3, 3)
    integer,    intent(in)  :: num_kgrid(3)
    complex(8), intent(out) :: xi(nb, nb, 3, nk)
    logical,    intent(out) :: xi_ok(nb, nb, nk)
    integer,    intent(out) :: n_reject
    real(8),    intent(out), optional :: resid_max
    logical,    intent(in),  optional :: fixed_blocks
    integer, allocatable :: block_id(:, :), ik_neighbor(:, :), link_member(:, :, :)
    logical, allocatable :: link_ok_arr(:, :)
    complex(8), allocatable :: M_blk(:, :), xi_blk(:, :)
    integer :: iv_axis(3), srcs(nb), tgts(nb)
    integer :: ik, iv, axis, bk, nblk_k, d, dt, a, bcol, n, n0, tgt0, tgt_block
    integer :: ikpb, info_x, n_fail
    logical :: blk_ok, use_fixed
    integer :: reject_reason
    real(8) :: dk_alpha, ru, ru_max, dwmin
    character(96) :: reason_msg

    xi       = (0d0, 0d0)
    xi_ok    = .false.
    n_reject = 0
    ru_max   = 0d0

    use_fixed = .false.
    if (present(fixed_blocks)) use_fixed = fixed_blocks

    allocate(block_id(nb, nk), ik_neighbor(nbvec, nk))
    allocate(link_member(nb, nbvec, nk), link_ok_arr(nbvec, nk))

    if (use_fixed) then
      call build_blocks_fixed(nb, nk, eigen, block_id)
    else
      call build_blocks(nb, nk, eigen, block_id)
    end if
    call build_ik_neighbor(num_kgrid, bvec, nbvec, nk, ik_neighbor)

    if (use_fixed) then
      ! Identity match: band n at k maps to band n at k+b for every in-block
      ! band, by construction (the composite block is k-independent, so its
      ! member SET is identical at k and k+b -- only the individual band
      ! LABEL inside the block is gauge-ambiguous, same caveat as
      ! match_link_blocks' greedy map). This eliminates reason-3 (tgt0==0)
      ! and reason-4 (image-block dimension mismatch) by construction; only
      ! xi_block_from_overlap's own reason-5 rejection can still occur.
      link_member = 0; link_ok_arr = .true.; n_fail = 0
      do ik = 1, nk
        do iv = 1, nbvec
          ikpb = ik_neighbor(iv, ik)
          if (ikpb < 1 .or. ikpb > nk) cycle
          do n = 1, nb
            link_member(n, iv, ik) = n
          end do
        end do
      end do
    else
      call match_link_blocks(nb, nk, nbvec, bvec, prod_dk, block_id, ik_neighbor, &
                             link_member, link_ok_arr, n_fail)
    end if

    iv_axis(1) = find_bvec(bvec, nbvec, 1, 0, 0)
    iv_axis(2) = find_bvec(bvec, nbvec, 0, 1, 0)
    iv_axis(3) = find_bvec(bvec, nbvec, 0, 0, 1)

    do ik = 1, nk
      nblk_k = maxval(block_id(:, ik))
      do bk = 1, nblk_k
        ! ascending source-block bands
        d = 0
        do n = 1, nb
          if (block_id(n, ik) == bk) then
            d = d + 1
            srcs(d) = n
          end if
        end do
        if (d < 2) cycle    ! singleton block: no off-diagonal dipole to regularise

        allocate(M_blk(d, d), xi_blk(d, d))
        blk_ok = .true.
        reject_reason = 0

        do axis = 1, 3
          iv = iv_axis(axis)
          if (iv == 0) then
            reject_reason = 1               ! +axis shift absent in prod_dk
            blk_ok = .false.; exit
          end if
          ikpb = ik_neighbor(iv, ik)
          if (ikpb < 1 .or. ikpb > nk) then
            reject_reason = 2               ! neighbor k index invalid
            blk_ok = .false.; exit
          end if
          ! per-block link validity (match_link_blocks zeros members of bad blocks)
          n0   = srcs(1)
          tgt0 = link_member(n0, iv, ik)
          if (tgt0 == 0) then
            reject_reason = 3               ! link zeroed by match_link_blocks
            blk_ok = .false.; exit
          end if
          tgt_block = block_id(tgt0, ikpb)
          ! ascending image-block bands, verify dimension conserved
          dt = 0
          do n = 1, nb
            if (block_id(n, ikpb) == tgt_block) then
              dt = dt + 1
              if (dt <= d) tgts(dt) = n
            end if
          end do
          if (dt /= d) then
            reject_reason = 4               ! image-block dimension mismatch
            blk_ok = .false.; exit
          end if
          ! overlap submatrix, ascending <-> ascending (standard finite difference)
          do bcol = 1, d
            do a = 1, d
              M_blk(a, bcol) = prod_dk(srcs(a), tgts(bcol), iv, ik)
            end do
          end do
          dk_alpha = b_matrix(axis, axis) / dble(num_kgrid(axis))
          call xi_block_from_overlap(M_blk, dk_alpha, xi_sign, xi_blk, info_x, ru)
          if (info_x /= 0) then
            reject_reason = 5               ! xi_block_from_overlap LAPACK failure
            n_reject = n_reject + 1
            blk_ok = .false.; exit
          end if
          if (ru > ru_max) ru_max = ru
          do bcol = 1, d
            do a = 1, d
              xi(srcs(a), srcs(bcol), axis, ik) = xi_blk(a, bcol)
            end do
          end do
        end do

        if (.not. blk_ok .and. use_fixed) then
          ! reason=5 (xi_block_from_overlap rejection) conflates three distinct
          ! runtime causes under one code; unpack info_x so a GS-setup abort is
          ! diagnosable without re-instrumenting the build.
          if (reject_reason == 5) then
            select case (info_x)
            case (1)
              reason_msg = 'reason=5 info=1 (near-singular overlap block)'
            case (2)
              reason_msg = 'reason=5 info=2 (branch-cut |phi|>0.9pi -- coarse-mesh symptom, densify k)'
            case (3)
              reason_msg = 'reason=5 info=3 (LAPACK failure)'
            case default
              write(reason_msg, '(a,i0,a)') 'reason=5 info=', info_x, ' (unrecognised)'
            end select
          else
            write(reason_msg, '(a,i2)') 'reason=', reject_reason
          end if
          write(*, '(a,i6,a,a,a,32i4)') &
            'ERROR build_xi(gifix): rejected fixed block ik=', ik, &
            ' ', trim(reason_msg), ' bands=', srcs(1:d)
          error stop 1
        end if

        if (blk_ok) then
          do bcol = 1, d
            do a = 1, d
              xi_ok(srcs(a), srcs(bcol), ik) = .true.
            end do
          end do
        else
          ! rejected / partial: drop this block's xi, blend falls back to i*p/dw
          dwmin = minval(abs(eigen(srcs(2:d), ik) - eigen(srcs(1:d-1), ik)))
          write(*, '(a,i6,a,i2,a,es10.3,a,32i4)') &
            'build_xi REJECT: ik=', ik, ' reason=', reject_reason, &
            ' dwmin_au=', dwmin, ' bands=', srcs(1:d)
          do bcol = 1, d
            do a = 1, d
              xi(srcs(a), srcs(bcol), :, ik) = (0d0, 0d0)
            end do
          end do
        end if

        deallocate(M_blk, xi_blk)
      end do
    end do

    deallocate(block_id, ik_neighbor, link_member, link_ok_arr)
    if (present(resid_max)) resid_max = ru_max
  end subroutine build_xi

  !===================================================================
  ! Phase 1 (gicov / Approach-B'): polar-only Wilson-line transport.
  !===================================================================

  !-------------------------------------------------------------------
  ! Polar factor of a d x d overlap block: U = M (M^H M)^{-1/2}, the
  ! unique unitary closest to M in Frobenius norm -- NO logm, so no
  ! eigenphase branch cut can occur here by construction (this routine
  ! never calls zgeev). sigma_min = sqrt(min eig(M^H M)) is the smallest
  ! singular value of M.
  !
  ! This is the polar half of xi_block_from_overlap (the G=M^H M / zheev
  ! / G^{-1/2} / U=M.G^{-1/2} steps), lifted into a standalone reusable
  ! kernel for build_block_transport; xi_block_from_overlap itself is
  ! left untouched (duplicated, not refactored, to avoid any risk of
  ! regressing its existing covariance test suite).
  !
  ! ierr = 0 ok
  !      = 1 near-singular: sigma_min < xi_sing_tol (garbage U NOT
  !          returned as trustworthy -- caller must reject, not use it)
  !      = 3 LAPACK (zheev) failure
  !-------------------------------------------------------------------
  subroutine polar_unitary(M, d, U, sigma_min, ierr)
    implicit none
    integer,    intent(in)  :: d
    complex(8), intent(in)  :: M(d, d)
    complex(8), intent(out) :: U(d, d)
    real(8),    intent(out) :: sigma_min
    integer,    intent(out) :: ierr
    integer :: k, lapack_info, lwork
    complex(8), allocatable :: A(:,:), Adg(:,:), Ginvsqrt(:,:)
    complex(8), allocatable :: cwork(:)
    real(8),    allocatable :: w_eig(:), rwork(:)

    ierr      = 0
    sigma_min = 0d0
    U         = (0d0, 0d0)

    allocate(A(d,d), Adg(d,d), Ginvsqrt(d,d), w_eig(d))

    ! --- G = M^H M ; Hermitian eigendecomposition (zheev): A->vecs, w_eig ascending
    A = matmul(conjg(transpose(M)), M)
    lwork = max(1, 2*d - 1)
    allocate(cwork(lwork), rwork(max(1, 3*d - 2)))
    call zheev('V', 'U', d, A, d, w_eig, cwork, lwork, rwork, lapack_info)
    deallocate(cwork, rwork)
    if (lapack_info /= 0) then
      ierr = 3
      deallocate(A, Adg, Ginvsqrt, w_eig)
      return
    end if

    ! --- near-singular reject: sigma_min(M) = sqrt(min eig of G) ---
    sigma_min = sqrt(max(0d0, w_eig(1)))
    if (sigma_min < xi_sing_tol) then
      ierr = 1
      deallocate(A, Adg, Ginvsqrt, w_eig)
      return
    end if

    ! --- G^{-1/2} = A diag(1/sqrt(w)) A^H ;  U = M G^{-1/2} ---
    do k = 1, d
      Adg(:, k) = A(:, k) / sqrt(w_eig(k))
    end do
    Ginvsqrt = matmul(Adg, conjg(transpose(A)))
    U = matmul(M, Ginvsqrt)

    deallocate(A, Adg, Ginvsqrt, w_eig)
  end subroutine polar_unitary

  !-------------------------------------------------------------------
  ! Block-diagonal parallel-transport (Wilson-line) operator built ONLY
  ! from the polar factor of each k -> k+e_alpha overlap block (see
  ! polar_unitary) -- never logm. Umat(nb,nb,3,nk): identity on
  ! singleton bands / bands outside any >=2 block and on any axis whose
  ! +unit-shift column is absent from bvec; the polar-unitary factor of
  ! the overlap submatrix on each block of the CALLER-SUPPLIED fixed
  ! composite partition `block_id` (this routine no longer computes the
  ! partition itself -- it is built upstream by either build_blocks_fixed
  ! (energy-only) or build_blocks_fixed_closed (overlap-closed, gicov
  ! X-closed) and passed in).
  !
  ! PRECONDITION on block_id (caller's responsibility): it MUST be a
  ! k-INDEPENDENT fixed partition, i.e. the source band set at k and the
  ! target band set at k+e_alpha are IDENTICAL for every block (srcs ==
  ! tgts). The body below relies on this: M_blk(a,b) =
  ! prod_dk(srcs(a), srcs(b), iv, ik) uses a single srcs list with no
  ! separate tgts map. Both supported builders satisfy this by
  ! construction; a per-k-varying partition would violate the
  ! precondition and produce a nonsensical overlap block.
  !
  ! Restricted to fixed blocks only (gicov's usage): source and target
  ! band sets at k and k+e are the SAME block by construction, so
  ! M_blk(a,b) = prod_dk(srcs(a), srcs(b), iv, ik) with no separate tgts
  ! map (contrast build_xi's general per-k block correspondence via
  ! match_link_blocks, which fixed blocks do not need).
  !
  ! Fail-closed: a near-singular block (polar_unitary ierr/=0) aborts
  ! the WHOLE build via error stop rather than returning a garbage U for
  ! that block -- there is no "untrusted-but-returned" state (no U_ok
  ! output), matching the fixed-block gifix contract used by build_xi's
  ! use_fixed branch.
  !
  ! n_reject is always 0 on return: it exists only so callers can log
  ! this interface uniformly alongside build_xi's n_reject; a rejected
  ! block error-stops before this subroutine could ever return a
  ! nonzero count.
  !-------------------------------------------------------------------
  ! ik_lo:ik_hi is the k-range this rank OWNS -- prod_dk and Umat are indexed
  ! over exactly that range, and the body is k-LOCAL by construction: M_blk(a,b)
  ! = prod_dk(srcs(a), srcs(b), iv, ik) reads the SAME ik, and the polar factor
  ! of one k's overlap block depends on nothing else.  That is what lets the
  ! whole static-field group be rank-local without touching a single value:
  ! serial / replicated callers pass 1, nk and get the identical Umat.
  subroutine build_block_transport(nb, nk, nbvec, bvec, prod_dk, block_id, num_kgrid, Umat, n_reject, &
                                 & ik_lo, ik_hi)
    implicit none
    integer,    intent(in)  :: nb, nk, nbvec
    integer,    intent(in)  :: bvec(3, nbvec), num_kgrid(3)
    integer,    intent(in)  :: ik_lo, ik_hi
    complex(8), intent(in)  :: prod_dk(nb, nb, nbvec, ik_lo:ik_hi)
    integer,    intent(in)  :: block_id(nb, nk)
    complex(8), intent(out) :: Umat(nb, nb, 3, ik_lo:ik_hi)
    integer,    intent(out) :: n_reject
    complex(8), allocatable :: M_blk(:, :), Ublk(:, :)
    integer :: iv_axis(3), srcs(nb)
    integer :: ik, iv, axis, bk, nblk_k, d, a, bcol, n
    real(8) :: sigma_min
    integer :: ierr

    ! num_kgrid is accepted only for interface symmetry with build_xi (and any
    ! future non-fixed-block extension). It is NOT needed here: fixed blocks
    ! have tgts==srcs at every k (brief-pinned), and the full periodic grid
    ! guarantees every k+e link exists, so the target k-index never has to be
    ! looked up (contrast build_xi's use of build_ik_neighbor for general
    ! per-k block correspondence).

    n_reject = 0
    if (ik_hi < ik_lo) return    ! empty k-slice (nproc > nk)

    iv_axis(1) = find_bvec(bvec, nbvec, 1, 0, 0)
    iv_axis(2) = find_bvec(bvec, nbvec, 0, 1, 0)
    iv_axis(3) = find_bvec(bvec, nbvec, 0, 0, 1)

    ! default: identity everywhere (singletons / out-of-block bands / axes
    ! whose +unit-shift column is absent from bvec)
    Umat = (0d0, 0d0)
    do ik = ik_lo, ik_hi
      do axis = 1, 3
        do n = 1, nb
          Umat(n, n, axis, ik) = (1d0, 0d0)
        end do
      end do
    end do

    do ik = ik_lo, ik_hi
      nblk_k = maxval(block_id(:, ik))
      do bk = 1, nblk_k
        d = 0
        do n = 1, nb
          if (block_id(n, ik) == bk) then
            d = d + 1
            srcs(d) = n
          end if
        end do
        if (d < 2) cycle    ! singleton block: transport is trivially identity

        do axis = 1, 3
          iv = iv_axis(axis)
          if (iv == 0) cycle   ! +axis shift absent in prod_dk: leave identity

          allocate(M_blk(d, d), Ublk(d, d))
          do bcol = 1, d
            do a = 1, d
              M_blk(a, bcol) = prod_dk(srcs(a), srcs(bcol), iv, ik)
            end do
          end do

          call polar_unitary(M_blk, d, Ublk, sigma_min, ierr)
          if (ierr /= 0) then
            write(*, '(a,i6,a,32i4)') &
              'ERROR build_block_transport: near-singular overlap block ik=', ik, &
              ' bands=', srcs(1:d)
            error stop 1
          end if

          do bcol = 1, d
            do a = 1, d
              Umat(srcs(a), srcs(bcol), axis, ik) = Ublk(a, bcol)
            end do
          end do
          deallocate(M_blk, Ublk)
        end do
      end do
    end do

  end subroutine build_block_transport

  !===================================================================
  ! gicov memory-scale fix helpers (LG-SBE Phase 1, X-full).
  !
  ! Root cause: gs_info_ssbe's read_prod_dk_data (the GS-file reader) stores
  ! gs%prod_dk(nb_eff,nb_eff,nbvec,nk) with nbvec = (2*ndk+1)**3 on EVERY
  ! rank regardless of nproc_k (no k-decomposition at read time) -- at ndk=1,
  ! nbvec=27, so this single array dominates per-rank memory
  ! (16*27*nb_eff**2*nk bytes) and SIGBUS'd at high nproc_k/large nk (k20^3,
  ! k24^3). But build_block_transport just above -- the ONLY prod_dk consumer
  ! on the 'gicov' (X-full) path -- only ever looks up the THREE +unit-shift
  ! columns via find_bvec(bvec,nbvec,1,0,0)/(0,1,0)/(0,0,1) (iv_axis(1:3)
  ! above); it never touches a negative or multi-shift direction, at any ndk.
  ! So for 'gicov' the reader only needs to KEEP those <=3 columns, cutting
  ! the dominant term by up to 9x at ndk=1 (more at larger ndk) with no
  ! change to build_block_transport's result (bit-for-bit: it reads the same
  ! 3 physical columns either way). gi/gifix are UNCHANGED (build_xi's
  ! match_link_blocks, degenerate_block_ssbe.f90:542, scans every bvec
  ! column, so they keep the full shell).
  !
  ! These helpers are pure (no gs_info_ssbe / SALMON dependency) and live
  ! HERE rather than in gs_info_ssbe.f90 for the same reason this whole
  ! module is pure (see the file header): gs_info_ssbe.f90 pulls in
  ! communication/filesystem/salmon_global/inputoutput/util_ssbe/common_ssbe
  ! at module scope (inside init_sbe_gs_info's local `use`s), so a standalone
  ! test of anything living IN that file needs the full SALMON dependency
  ! tree built; a pure function living in THIS module needs only this file
  ! (gfortran degenerate_block_ssbe.f90 test/test_gicov_memfix.f90 -o t, same
  ! convention as test_block_transport.f90 / test_gicov_xfull.f90).
  ! gs_info_ssbe.f90's read_prod_dk_data calls these three directly.
  !===================================================================

  !-------------------------------------------------------------------
  ! Which of the <=3 KEPT gicov columns a raw (jdk1,jdk2,jdk3) dk-shift
  ! (as parsed from a file_sbe_prod_dk record) belongs to, if any:
  !   1 = +x (1,0,0), 2 = +y (0,1,0), 3 = +z (0,0,1), 0 = dropped.
  ! Fixed slot order (not file-dependent) so gs%bvec can be pre-populated
  ! once (see gicov_prod_dk_nbvec) instead of built by enumerating the file's
  ! own (2*ndk+1)**3 shift order.
  !-------------------------------------------------------------------
  pure integer function prod_dk_axis_slot(jdk1, jdk2, jdk3) result(islot)
    implicit none
    integer, intent(in) :: jdk1, jdk2, jdk3
    if (jdk1 == 1 .and. jdk2 == 0 .and. jdk3 == 0) then
      islot = 1
    else if (jdk1 == 0 .and. jdk2 == 1 .and. jdk3 == 0) then
      islot = 2
    else if (jdk1 == 0 .and. jdk2 == 0 .and. jdk3 == 1) then
      islot = 3
    else
      islot = 0
    end if
  end function prod_dk_axis_slot

  !-------------------------------------------------------------------
  ! STORED 4th-dim extent of gs%prod_dk (what read_prod_dk_data allocates
  ! and gs%nbvec is set to) -- NOT the file's own (2*ndk+1)**3 = nbvec_full,
  ! which the reader must still use unchanged to compute how many records
  ! the file contains (see expected_prod_dk_nrec: the writer's record count
  ! does not shrink just because the reader keeps fewer columns).
  !
  ! 'gicov' with ndk>=1 keeps only the 3 axis columns (the +unit-shift
  ! records exist in the file whenever ndk>=1, since the per-axis range is
  ! -ndk..ndk and 1 is inside it for any ndk>=1). Every other case --
  ! gi/gifix/off, OR 'gicov' with ndk==0 (no neighbour-shift data at all,
  ! nbvec_full is already 1) -- keeps the full shell UNCHANGED, so this
  ! degrades to a no-op identity for every pre-existing mode/edge-case.
  !-------------------------------------------------------------------
  pure integer function gicov_prod_dk_nbvec(sbe_lg_degen, ndk, nbvec_full) result(nbvec_out)
    implicit none
    character(*), intent(in) :: sbe_lg_degen
    integer, intent(in) :: ndk, nbvec_full
    if (uses_xfull_links(sbe_lg_degen) .and. ndk >= 1) then
      nbvec_out = 3
    else
      nbvec_out = nbvec_full
    end if
  end function gicov_prod_dk_nbvec

  !-------------------------------------------------------------------
  ! Writer's total record count for 'file_sbe_prod_dk' = nk*(2*ndk+1)**3*no**2
  ! (ik outermost/slowest; see read_prod_dk_data), computed in integer(8) so
  ! it stays exact past 2**31-1 -- hit already at a k20^3=8000 k-point mesh
  ! for a realistic export band count `no` (this was the previous default-
  ! integer overflow: nk*nbvec*file_no**2 silently wrapped negative/aliased,
  ! corrupting the "fewer/more records than expected" gate and the per-record
  ! ik_exp ordering check). Independent of gs_info state (no gs%/file I/O),
  ! so it is unit-testable without a real prod_dk.data file or MPI.
  !-------------------------------------------------------------------
  pure function expected_prod_dk_nrec(nk, nbvec_full, file_no) result(nrec)
    implicit none
    integer, intent(in) :: nk, nbvec_full, file_no
    integer(8) :: nrec
    nrec = int(nk, 8) * int(nbvec_full, 8) * int(file_no, 8) * int(file_no, 8)
  end function expected_prod_dk_nrec

  !===================================================================
  ! Phase 2 (gicov / Approach-B'): pure gauge-covariant intraband
  ! k-derivative operator.
  !===================================================================

  !-------------------------------------------------------------------
  ! Gauge-covariant intraband k-derivative
  !
  !     D_cov rho = d_k rho - i[xi, rho]
  !
  ! evaluated on the FULL nb x nb density using the block-diagonal Wilson-line
  ! transport U_full (build_block_transport: unitary on each fixed composite
  ! block, identity elsewhere; U_full(k) = <u(k)|u(k+e)> ~ I - i xi dk, polar,
  ! no logm) to parallel-transport the neighbour densities into the k-frame
  ! before an UP-TO-4-shell (+-1..+-m_max) CENTRAL difference, per-axis capped
  ! (see m_max below).  Weights
  !     c = (4/5, -1/5, 4/105, -1/280)  =  bNmat(:,4)
  ! are the SAME shells/weights as the production length-gauge k-gradient
  ! grad_k_array_nb2d_dcomplex (common_ssbe.f90:104-194; weights set by set_bN,
  ! src/common/initialization.f90:1343); a capped axis simply omits the
  ! highest shell(s) (truncated finite-difference order, weights NOT
  ! renormalised) rather than dropping the operator. No SALMON deps, no MPI
  ! (full-k arrays).
  !
  ! Per-axis adaptive shell cap (aliasing guard): shell m couples k to k+-m e,
  ! and on a periodic grid of size n = num_kgrid(alpha), k+m e and k-m e are
  ! the SAME index once 2m >= n. For the BARE stencil (U_full == I) the two
  ! terms are then literally the same rho value subtracted from itself and
  ! cancel exactly; but for the COVARIANT stencil the two composed Wilson
  ! lines Um(k) (forward, via k+e,k+2e,...) and Um(k-m e) (backward, via
  ! k-e,k-2e,...) reach that SAME point by two DIFFERENT paths around the
  ! periodic ring and are generally different unitaries, so the forward and
  ! backward transported terms do NOT cancel -- leaving a spurious O(1)
  ! contribution (this bit the production 8^3 mesh, where shell 4 aliases:
  ! k+4 == k-4 mod 8). To make this structurally impossible at ANY mesh size,
  ! each axis independently caps its shell count to
  !     m_max(alpha) = min(4, (num_kgrid(alpha) - 1) / 2)     (integer division)
  ! so that 2*m_max(alpha) < num_kgrid(alpha) always holds and no shell can
  ! alias. m_max=4 for num_kgrid(alpha)>=9 (unchanged full 4-shell, matches
  ! the legacy bare gradient's effective behaviour there); m_max=3 at the
  ! production nk_axis=8 (drops the aliasing shell 4 only); m_max=0 for a
  ! singleton axis (num_kgrid(alpha)=1), i.e. Dq(:,:,alpha,:) = 0.
  !
  ! For each Cartesian axis alpha, k-point k, shell m = 1..m_max(alpha):
  !   forward neighbour  k + m e  reached by composing the +e_alpha link m times;
  !   composed Wilson link (product of single FORWARD links)
  !     Um(k)     = U_full(k) U_full(k+e) ... U_full(k+(m-1)e)   [k -> k+m e]
  !     Um(k-m e) = U_full(k-m e) ... U_full(k-e)                [k-m e -> k]
  !   contribution
  !     c_m * [ Um(k) rho(k+m e) Um(k)^H  -  Um(k-m e)^H rho(k-m e) Um(k-m e) ] / dk(alpha)
  !   (forward neighbour transported by Um(k): U rho U^H; backward neighbour by
  !    Um(k-m e)^H: U^H rho U).  Dq(:,:,alpha,k) = sum_m of those.
  !
  ! With U_full == I this reduces EXACTLY to the (per-axis capped) bare central
  ! difference of rho (Test A, m_max=4 at its nk=12). It is gauge COVARIANT:
  ! under rho -> W(k)^H rho W(k) and U_full(k) -> W(k)^H U_full(k) W(k+e) the
  ! output transforms as a rank-2 tensor Dq(k) -> W(k)^H Dq(k) W(k) (Test B,
  ! telescoping W W^H = I). For constant xi with U_full = exp(-i xi dk) it
  ! equals d_k rho - i[xi,rho] (via the same capped stencil) to high order
  ! (Test C at nk=32/m_max=4, 1e-8 R-5 SIGN gate; Test D at the production
  ! nk=8/m_max=3, alias-free check).
  !
  ! Triple products use ZGEMM (zmm3).  ik_neighbor is built internally with
  ! build_ik_neighbor; the +axis columns are resolved with find_bvec; the
  ! backward neighbour map is the inverse permutation of the +axis forward map
  ! (well-defined on the full periodic grid).  An axis whose +unit-shift column
  ! is absent from bvec, OR whose m_max(axis)=0 (singleton axis), leaves
  ! Dq(:,:,axis,:) = 0.
  !-------------------------------------------------------------------
  ! kmap/nks: U_full and rho are stored on nks SLOTS, not on the global k-index
  ! -- kmap(ik) is where k lives (0 = not stored here).  A replicated caller
  ! passes nks = nk and the identity map, which reproduces the original
  ! `U_full(:,:,axis,ik)` indexing exactly; a k-local caller passes its owned
  ! slots plus the halo slots this very stencil's chain walk asked for
  ! (covariant_halo_needed), so every k read below is present by construction.
  subroutine covariant_grad_block(nb, nk, nbvec, bvec, num_kgrid, U_full, rho, dk, Dq, ik_lo, ik_hi, &
                                & axis_active, nks, kmap)
    implicit none
    integer,    intent(in)  :: nb, nk, nbvec, bvec(3, nbvec), num_kgrid(3)
    integer,    intent(in)  :: ik_lo, ik_hi   ! local k-range (MPI): only Dq(:,:,:,ik_lo:ik_hi) is computed, OMP-parallel over it. Serial/single-rank callers pass 1, nk. (Required, not optional: an interface mismatch then fails at compile time, not as a silent runtime present()-garbage SEGV.)
    integer,    intent(in)  :: nks             ! number of stored k slots (= nk when replicated)
    integer,    intent(in)  :: kmap(nk)        ! global k -> slot; 0 = not stored (never read here: the halo plan covers the stencil's own walk)
    complex(8), intent(in)  :: U_full(nb, nb, 3, nks)
    complex(8), intent(in)  :: rho(nb, nb, nks)
    real(8),    intent(in)  :: dk(3)
    complex(8), intent(out) :: Dq(nb, nb, 3, ik_lo:ik_hi)
    logical,    intent(in)  :: axis_active(3)  ! skip axes with E(axis)==0: their E*Dq contribution is exactly 0, so their Dq is left 0 (bit-for-bit). Linear pol => 1 axis, circular in-plane => 2. Callers with no field mask pass (/.true.,.true.,.true./).
    real(8), parameter :: cw(4) = (/ 4d0/5d0, -1d0/5d0, 4d0/105d0, -1d0/280d0 /)
    complex(8), allocatable :: Id(:,:)
    integer :: iv_axis(3), axis, iv, ik, n
    integer :: m_max(3), klo, khi

    allocate(Id(nb,nb))
    klo = ik_lo;  khi = ik_hi

    ! Neighbour tables are a pure function of (num_kgrid, bvec) -- the SAME
    ! tables every call.  Rebuilding them here cost an O(nbvec*nk) allocate +
    ! fill FOUR times per RK4 step; the module cache below builds them once and
    ! returns bit-identical indices (integer index maps, no arithmetic).
    call get_ik_neighbor_cache(num_kgrid, bvec, nbvec, nk)
    iv_axis(1) = find_bvec(bvec, nbvec, 1, 0, 0)
    iv_axis(2) = find_bvec(bvec, nbvec, 0, 1, 0)
    iv_axis(3) = find_bvec(bvec, nbvec, 0, 0, 1)

    ! Per-axis adaptive shell cap: 2*m_max(axis) < num_kgrid(axis) always, so
    ! shell m_max(axis) can never alias (k+m == k-m mod num_kgrid(axis)). See
    ! the subroutine header for the full rationale.
    do axis = 1, 3
      m_max(axis) = min(4, (num_kgrid(axis) - 1) / 2)
    end do

    Id = (0d0, 0d0)
    do n = 1, nb
      Id(n, n) = (1d0, 0d0)
    end do

    ! Zero only the local k-slice (all 3 axes). The caller contract is that
    ! Dq outside [ik_lo,ik_hi] is UNDEFINED (gicov_rhs's drho loop reads Dq
    ! at local ik only) -- zeroing the full nk array was a ~nb^2*3*nk memset
    ! per call (4x per RK4 step), the dominant fixed per-call floor observed
    ! in the strong-scaling scan. Skipped/absent axes at local ik stay 0 here,
    ! which the E*Dq sum in gicov_rhs relies on.
    Dq(:, :, :, klo:khi) = (0d0, 0d0)

    do axis = 1, 3
      if (.not. axis_active(axis)) cycle  ! field-inactive axis (E(axis)==0): Dq(:,:,axis,:) stays 0; E*Dq is exactly 0 => bit-for-bit skip (linear pol skips 2 axes, circular in-plane skips 1)
      iv = iv_axis(axis)
      if (iv == 0) cycle              ! +axis link absent in bvec: Dq(:,:,axis,:) stays 0
      if (m_max(axis) == 0) cycle     ! singleton axis (num_kgrid(axis)=1): Dq(:,:,axis,:) stays 0

      ! OMP over the local k-slice. The per-k body runs in cov_grad_one_k, whose
      ! nb x nb work matrices are automatic (stack) arrays = thread-local by
      ! construction (robust on gfortran and Fujitsu frt; avoids the flaky
      ! OMP-private-allocatable path). Dq(:,:,axis,ik) is a disjoint slice per
      ! ik, so the writes never race. Numerically identical to the old inline loop.
      !$omp parallel do default(shared) private(ik) schedule(static)
      do ik = klo, khi
        call cov_grad_one_k(nb, nbvec, nk, nks, ik, iv, axis, m_max(axis), cw, dk(axis), &
                            ikn_store, bwd_store(:, axis), kmap, U_full, rho, Id, Dq(:, :, axis, ik))
      end do
      !$omp end parallel do
    end do

    deallocate(Id)
  end subroutine covariant_grad_block


  !-------------------------------------------------------------------
  ! One-shot cache of the k-index neighbour tables used by the covariant
  ! stencil: the forward map ik_neighbor(iv, ik) = index of k + bvec(:,iv), and
  ! the per-axis backward map bwd(k + e_axis) = k (the inverse of the +axis
  ! permutation, well defined on the full periodic grid).  Both are pure integer
  ! functions of (num_kgrid, bvec), so caching them is bit-neutral; they were
  ! being rebuilt on every covariant_grad_block call = 4x per RK4 step.
  ! The cache key is (num_kgrid, bvec, nk): a different grid rebuilds it.
  !-------------------------------------------------------------------
  ! The slot map of a REPLICATED k-layout: every k is stored, at its own index.
  ! Serial and full-nk callers pass this, and every kmap-indexed expression then
  ! reduces to the original global-k one.
  pure function identity_kmap(nk) result(kmap)
    implicit none
    integer, intent(in) :: nk
    integer :: kmap(nk)
    integer :: ik
    do ik = 1, nk
      kmap(ik) = ik
    end do
  end function identity_kmap

  subroutine get_ik_neighbor_cache(num_kgrid, bvec, nbvec, nk)
    implicit none
    integer, intent(in) :: num_kgrid(3), nbvec, nk, bvec(3, nbvec)
    integer :: iv_axis(3), axis, iv, ik
    logical :: stale

    stale = .not. ikn_built
    if (.not. stale) then
      stale = (ikn_nk /= nk) .or. (ikn_nbvec /= nbvec) .or. &
            & any(ikn_kgrid(1:3) /= num_kgrid(1:3))
    end if
    if (.not. stale) then
      stale = any(ikn_bvec(1:3, 1:nbvec) /= bvec(1:3, 1:nbvec))
    end if
    if (.not. stale) return

    if (allocated(ikn_store))  deallocate(ikn_store)
    if (allocated(bwd_store))  deallocate(bwd_store)
    if (allocated(ikn_bvec))   deallocate(ikn_bvec)
    allocate(ikn_store(nbvec, nk), bwd_store(nk, 3), ikn_bvec(3, nbvec))

    call build_ik_neighbor(num_kgrid, bvec, nbvec, nk, ikn_store)
    iv_axis(1) = find_bvec(bvec, nbvec, 1, 0, 0)
    iv_axis(2) = find_bvec(bvec, nbvec, 0, 1, 0)
    iv_axis(3) = find_bvec(bvec, nbvec, 0, 0, 1)
    bwd_store(:, :) = 0
    do axis = 1, 3
      iv = iv_axis(axis)
      if (iv == 0) cycle
      do ik = 1, nk
        bwd_store(ikn_store(iv, ik), axis) = ik
      end do
    end do

    ikn_bvec(1:3, 1:nbvec) = bvec(1:3, 1:nbvec)
    ikn_kgrid(1:3) = num_kgrid(1:3)
    ikn_nk = nk;  ikn_nbvec = nbvec
    ikn_built = .true.
  end subroutine get_ik_neighbor_cache

  !-------------------------------------------------------------------
  ! Per-k covariant-gradient contribution for one axis (called under OMP).
  !   Dqk = sum_m cw(m) [ Um(k) rho(k+m e) Um(k)^H
  !                       - Um(k-m e)^H rho(k-m e) Um(k-m e) ] / dkax
  ! All nb x nb work matrices are automatic (per-call = thread-local), so this
  ! is safe to invoke from a parallel do over ik. Bit-for-bit identical to the
  ! original inlined loop body.
  !-------------------------------------------------------------------
  subroutine cov_grad_one_k(nb, nbvec, nk, nks, ik, iv, axis, mmax, cw, dkax, &
                            ik_neighbor, bwd, kmap, U_full, rho, Id, Dqk)
    implicit none
    integer,    intent(in)  :: nb, nbvec, nk, nks, ik, iv, axis, mmax
    real(8),    intent(in)  :: cw(4), dkax
    integer,    intent(in)  :: ik_neighbor(nbvec, nk), bwd(nk), kmap(nk)
    complex(8), intent(in)  :: U_full(nb, nb, 3, nks), rho(nb, nb, nks), Id(nb, nb)
    complex(8), intent(out) :: Dqk(nb, nb)
    complex(8) :: Umf(nb,nb), Umb(nb,nb), Utmp(nb,nb), tmp(nb,nb), fterm(nb,nb), bterm(nb,nb)
    integer :: m, kfwd, kbwd, krp

    Umf  = Id
    Umb  = Id
    kfwd = ik      ! position of the next forward link to append = k+(m-1)e
    kbwd = ik      ! stepped to k-m e inside the loop
    Dqk  = (0d0, 0d0)
    do m = 1, mmax
      ! ---- forward: Um(k) <- Um(k) * U_full(k+(m-1)e); neighbour = k+m e ----
      call zmm3('N', 'N', nb, Umf, U_full(:, :, axis, kmap(kfwd)), Utmp)
      Umf = Utmp
      krp = ik_neighbor(iv, kfwd)                     ! k+m e
      call zmm3('N', 'N', nb, Umf, rho(:, :, kmap(krp)), tmp)   ! Umf * rho(k+m e)
      call zmm3('N', 'C', nb, tmp, Umf, fterm)            ! ... * Umf^H

      ! ---- backward: base = k-m e; Um(k-m e) <- U_full(k-m e) * Um(k-(m-1)e) ----
      kbwd = bwd(kbwd)                                 ! k-m e
      call zmm3('N', 'N', nb, U_full(:, :, axis, kmap(kbwd)), Umb, Utmp)
      Umb = Utmp
      call zmm3('C', 'N', nb, Umb, rho(:, :, kmap(kbwd)), tmp) ! Umb^H * rho(k-m e)
      call zmm3('N', 'N', nb, tmp, Umb, bterm)           ! ... * Umb

      Dqk = Dqk + cw(m) * (fterm - bterm) / dkax
      kfwd = krp                                       ! next forward link at k+m e
    end do
  end subroutine cov_grad_one_k

  !-------------------------------------------------------------------
  ! Halo support (pure, no MPI): mark every k whose rho is READ when
  ! covariant_grad_block runs on the local slice [ik_lo, ik_hi]. Mirrors the
  ! chain-walk of cov_grad_one_k exactly (forward k+m e via ik_neighbor,
  ! backward k-m e via the inverse map), over ALL three axes = a superset of
  ! any axis_active mask, so the resulting halo plan is static even though
  ! the field direction varies in time. needed(ik_lo:ik_hi) is also set (local
  ! k included).  An empty slice (ik_lo > ik_hi, possible when nproc > nk)
  ! yields all-false.
  !
  ! The SAME set serves the STATIC Wilson links: cov_grad_one_k reads U_full at
  ! {k+(m-1)e} and {k-m e} for m = 1..m_max, both strict subsets of the {k+-m e}
  ! marked here, so one plan covers the rho halo AND the u_transport halo (that
  ! is why gs%kmap can be built from this routine alone).
  !-------------------------------------------------------------------
  subroutine covariant_halo_needed(nk, nbvec, bvec, num_kgrid, ik_lo, ik_hi, needed)
    implicit none
    integer, intent(in)  :: nk, nbvec, bvec(3, nbvec), num_kgrid(3)
    integer, intent(in)  :: ik_lo, ik_hi
    logical, intent(out) :: needed(nk)
    integer, allocatable :: ik_neighbor(:,:), bwd(:)
    integer :: iv_axis(3), m_max(3), axis, iv, ik, m, kfwd, kbwd

    needed = .false.
    if (ik_lo > ik_hi) return
    needed(ik_lo:ik_hi) = .true.

    allocate(ik_neighbor(nbvec, nk), bwd(nk))
    call build_ik_neighbor(num_kgrid, bvec, nbvec, nk, ik_neighbor)
    iv_axis(1) = find_bvec(bvec, nbvec, 1, 0, 0)
    iv_axis(2) = find_bvec(bvec, nbvec, 0, 1, 0)
    iv_axis(3) = find_bvec(bvec, nbvec, 0, 0, 1)
    do axis = 1, 3
      m_max(axis) = min(4, (num_kgrid(axis) - 1) / 2)
    end do

    do axis = 1, 3
      iv = iv_axis(axis)
      if (iv == 0) cycle
      if (m_max(axis) == 0) cycle
      do ik = 1, nk
        bwd(ik_neighbor(iv, ik)) = ik
      end do
      do ik = ik_lo, ik_hi
        kfwd = ik
        kbwd = ik
        do m = 1, m_max(axis)
          kfwd = ik_neighbor(iv, kfwd)   ! k+m e : rho(kfwd) is read
          needed(kfwd) = .true.
          kbwd = bwd(kbwd)               ! k-m e : rho(kbwd) is read
          needed(kbwd) = .true.
        end do
      end do
    end do
    deallocate(ik_neighbor, bwd)
  end subroutine covariant_halo_needed

  !-------------------------------------------------------------------
  ! Halo support (pure, no MPI): derive per-rank send/recv k-lists from the
  ! needed-sets of ALL ranks (needed_all(:,o) = covariant_halo_needed of rank
  ! o's slice). Deterministic ascending-k ordering on both sides, so sender
  ! and receiver agree on packing order with no metadata exchange:
  !   recv: from each owner o /= myrank, the k in o's range that I need
  !   send: to each requester r /= myrank, the k in MY range that r needs
  ! Flattened lists: partner p occupies k( ofs(p)+1 : ofs(p)+cnt(p) ).
  !-------------------------------------------------------------------
  subroutine build_halo_lists(nk, nproc, myrank, itbl_min, itbl_max, needed_all, &
                              nsrc, src_rank, src_cnt, src_ofs, src_k, &
                              ndst, dst_rank, dst_cnt, dst_ofs, dst_k)
    implicit none
    integer, intent(in)  :: nk, nproc, myrank
    integer, intent(in)  :: itbl_min(0:nproc-1), itbl_max(0:nproc-1)
    logical, intent(in)  :: needed_all(nk, 0:nproc-1)
    integer, intent(out) :: nsrc, ndst
    integer, allocatable, intent(out) :: src_rank(:), src_cnt(:), src_ofs(:), src_k(:)
    integer, allocatable, intent(out) :: dst_rank(:), dst_cnt(:), dst_ofs(:), dst_k(:)
    integer :: o, r, ik, c, tot, p

    ! ---- recv side: owners I need k from ----
    nsrc = 0;  tot = 0
    do o = 0, nproc - 1
      if (o == myrank) cycle
      c = 0
      do ik = itbl_min(o), itbl_max(o)
        if (needed_all(ik, myrank)) c = c + 1
      end do
      if (c > 0) then
        nsrc = nsrc + 1;  tot = tot + c
      end if
    end do
    allocate(src_rank(nsrc), src_cnt(nsrc), src_ofs(nsrc), src_k(tot))
    p = 0;  tot = 0
    do o = 0, nproc - 1
      if (o == myrank) cycle
      c = 0
      do ik = itbl_min(o), itbl_max(o)
        if (needed_all(ik, myrank)) c = c + 1
      end do
      if (c > 0) then
        p = p + 1
        src_rank(p) = o;  src_cnt(p) = c;  src_ofs(p) = tot
        do ik = itbl_min(o), itbl_max(o)                 ! ascending k
          if (needed_all(ik, myrank)) then
            tot = tot + 1;  src_k(tot) = ik
          end if
        end do
      end if
    end do

    ! ---- send side: requesters that need MY k ----
    ndst = 0;  tot = 0
    do r = 0, nproc - 1
      if (r == myrank) cycle
      c = 0
      do ik = itbl_min(myrank), itbl_max(myrank)
        if (needed_all(ik, r)) c = c + 1
      end do
      if (c > 0) then
        ndst = ndst + 1;  tot = tot + c
      end if
    end do
    allocate(dst_rank(ndst), dst_cnt(ndst), dst_ofs(ndst), dst_k(tot))
    p = 0;  tot = 0
    do r = 0, nproc - 1
      if (r == myrank) cycle
      c = 0
      do ik = itbl_min(myrank), itbl_max(myrank)
        if (needed_all(ik, r)) c = c + 1
      end do
      if (c > 0) then
        p = p + 1
        dst_rank(p) = r;  dst_cnt(p) = c;  dst_ofs(p) = tot
        do ik = itbl_min(myrank), itbl_max(myrank)       ! ascending k
          if (needed_all(ik, r)) then
            tot = tot + 1;  dst_k(tot) = ik
          end if
        end do
      end if
    end do
  end subroutine build_halo_lists

  !-------------------------------------------------------------------
  ! C = op(A) op(B), all n x n complex, via ZGEMM (BLAS).  ta/tb in {'N','C'}.
  !-------------------------------------------------------------------
  subroutine zmm3(ta, tb, n, A, B, C)
    implicit none
    character(1), intent(in)  :: ta, tb
    integer,      intent(in)  :: n
    complex(8),   intent(in)  :: A(n, n), B(n, n)
    complex(8),   intent(out) :: C(n, n)
    call zgemm(ta, tb, n, n, n, (1d0, 0d0), A, n, B, n, (0d0, 0d0), C, n)
  end subroutine zmm3

  !-------------------------------------------------------------------
  ! Metal detector for gs%occup(nb,nk).  is_metal=.true. means "the insulator
  ! T=0 rigid nvb-fill (gs%occup(1:nvb,:)=focc, same nvb bands full at every k)
  ! would be WRONG here -- keep the real per-(band,k) GS occupation instead."
  !
  ! Two independent signatures, either one sufficient:
  !   (1) genuine fractional occupation (finite-temperature smearing): any
  !       occup(ib,ik) strictly between (tol, focc-tol).
  !   (2) a k-VARYING occupied-band pattern (Fermi surface), which happens even
  !       at temperature==0 for a metal: either the per-k COUNT of occupied
  !       bands (occup>tol) is not the same at every k, or -- at fixed count --
  !       the occupied bands are not the CONTIGUOUS bottom block 1..n at some k
  !       (an empty band sits below a filled one: band crossing E_F away from
  !       k=Gamma). Both are impossible for an insulator by construction (the
  !       gap keeps exactly nvb bands full, contiguously, at every k).
  !
  ! Pure, no SALMON deps -- unit-testable the same way as this module's other
  ! standalone helpers (see the file header / test/test_degenerate_block.f90).
  !-------------------------------------------------------------------
  pure subroutine check_gicov_occupation(nb, nk, occup, focc, is_metal, frac_tol)
    implicit none
    integer,  intent(in)  :: nb, nk
    real(8),  intent(in)  :: occup(nb, nk)
    real(8),  intent(in)  :: focc
    logical,  intent(out) :: is_metal
    real(8),  intent(in), optional :: frac_tol
    real(8) :: tol
    integer :: ik, ib, ncount, ncount0

    tol = 1d-6
    if (present(frac_tol)) tol = frac_tol

    is_metal = .false.
    ncount0  = -1
    do ik = 1, nk
        ncount = 0
        do ib = 1, nb
            if (occup(ib, ik) > tol .and. occup(ib, ik) < focc - tol) is_metal = .true.  ! (1) fractional
            if (occup(ib, ik) > tol) then
                ncount = ncount + 1
            end if
            if (ib > 1) then
                if (occup(ib, ik) > tol .and. occup(ib - 1, ik) <= tol) is_metal = .true. ! (2) non-contiguous
            end if
        end do
        if (ncount0 < 0) then
            ncount0 = ncount
        else if (ncount /= ncount0) then
            is_metal = .true.                                                            ! (2) k-varying count
        end if
    end do
  end subroutine check_gicov_occupation

end module degenerate_block_ssbe
