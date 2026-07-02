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
  public :: build_fixed_blocks
  public :: build_blocks_fixed
  public :: same_block
  public :: match_link_blocks
  ! Pb3 (non-Abelian xi + smooth blend):
  public :: blend                 ! cos^2 hysteresis ramp 1 -> 0
  public :: build_ik_neighbor     ! k -> k+b index table on the full uniform grid
  public :: xi_block_from_overlap ! polar + matrix-log kernel (LAPACK zheev/zgeev)
  public :: build_xi              ! orchestrator: blocks -> links -> xi(nb,nb,3,nk)
  public :: xi_sign               ! sign s in xi = s*i*logm(U)/dk (pinned by test)

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
  ! k-INDEPENDENT composite-block partition, determined once from the full
  ! eigen(nb,nk) spectrum: connected components of the graph whose edges are
  ! band pairs near-degenerate (< theta_off) at ANY k (transitive union-find,
  ! same threshold as build_blocks' per-k grouping, but unioned across all k
  ! so a pair that is close only at one k still fixes the two bands into one
  ! block for every k). This is the algorithmic core of the Tier2 gifix mode:
  ! blocks so defined are dimensionally stable across k by construction, so
  ! Pb3's per-link dimension-mismatch failure mode cannot occur.
  !
  ! isolated_ok verifies the partition is physically meaningful: at every k, a
  ! MULTI-band block's members must stay farther than gap_margin from every
  ! band outside the block. gap_margin (1d-2 au) is deliberately LARGER than
  ! theta_off (2d-3 au, the grouping bound) -- otherwise a band that sits just
  ! outside theta_off (so it never joins the block) could still be arbitrarily
  ! close to it without ever tripping the isolation check, making the check
  ! tautological. isolated_ok=.false. signals the ambiguous / metal-like case
  ! (a would-be outside band drifting within gap_margin of a multi-band block
  ! at some k) and the caller fails closed rather than silently mis-treating
  ! it as an isolated insulating block. Singleton-vs-singleton proximity does
  ! NOT trip isolation (two merely-close-but-distinct non-degenerate bands are
  ! not an ambiguity this routine needs to flag).
  !-------------------------------------------------------------------
  subroutine build_fixed_blocks(nb, nk, eigen, fixed_block_id, isolated_ok)
    implicit none
    integer, intent(in)  :: nb, nk
    real(8), intent(in)  :: eigen(nb, nk)
    integer, intent(out) :: fixed_block_id(nb)
    logical, intent(out) :: isolated_ok
    real(8), parameter :: gap_margin = 1d-2   ! au; MUST be > theta_off so a band can be
                                               ! outside a block (>theta_off) yet fail isolation (<gap_margin)
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
          write(*, '(a,i6,a,i2,a,32i4)') &
            'ERROR build_xi(gifix): rejected fixed block ik=', ik, &
            ' reason=', reject_reason, ' bands=', srcs(1:d)
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

end module degenerate_block_ssbe
