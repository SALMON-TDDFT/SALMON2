! src/ssbe/test/test_gicov_klocal.f90
!
! Unit test for the RANK-LOCAL static-field layout (perf/sbe-static-fields-local-k).
!
! The k-decomposition of the static per-k matrices is a pure STORAGE change: the
! same numbers, held on fewer ranks, reached through a slot map instead of the
! global k index.  The two properties that make that claim true are exactly the
! two gated here, and both are gated as BIT-identity (not "close"):
!
!   A. build_block_transport is k-LOCAL.  Running it on a sub-range [lo,hi]
!      must reproduce, bit for bit, what the full-range run produced on those k.
!      This is what allows prod_dk / u_transport to be sliced at all: if the
!      polar factor at k depended on a neighbour, slicing would change values.
!
!   B. covariant_grad_block is INVARIANT under the slot map.  Feeding it a
!      SCRAMBLED slot layout (deliberately NOT the identity, and not even
!      contiguous) with the matching kmap must give bit-for-bit the same Dq as
!      the replicated identity layout.  This is what allows u_transport/rho to be
!      stored as "owned + halo" slots instead of full-nk arrays.
!
! A permutation is a stronger test than the real layout: the production kmap is
! monotone on the owned block, so an off-by-one in the slot arithmetic could
! cancel; under a scramble it cannot.
!
! PURE KERNEL MATH: only degenerate_block_ssbe -- no gs_info / bloch_solver / MPI.
!
! BUILD (standalone):
!   gfortran -O2 -ffree-line-length-none -fallow-argument-mismatch -w \
!     src/ssbe/sbe_lg_mode_ssbe.f90 src/ssbe/degenerate_block_ssbe.f90 \
!     src/ssbe/test/test_gicov_klocal.f90 -o t -framework Accelerate && ./t
program test_gicov_klocal
  use degenerate_block_ssbe, only: build_block_transport, covariant_grad_block, &
                                   covariant_halo_needed, identity_kmap
  implicit none
  integer :: nfail
  nfail = 0

  call test_transport_klocal(nfail)
  call test_grad_slot_invariance(nfail)

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  end if
  write(*, '(a)') "All test_gicov_klocal checks passed."

contains

  ! Smooth unitary frame V(k) = exp(i H(k)) -> genuinely unitary overlaps with a
  ! non-trivial non-Abelian connection (polar_unitary stays well conditioned).
  subroutine make_fixture(nb, nk, n1, n2, nbvec, bvec, num_kgrid, prod_dk, rho)
    integer, intent(in) :: nb, nk, n1, n2, nbvec
    integer, intent(out) :: bvec(3, nbvec), num_kgrid(3)
    complex(8), intent(out) :: prod_dk(nb, nb, nbvec, nk), rho(nb, nb, nk)
    integer :: ik, jk, iv, i, j, i1, i2, j1, j2
    complex(8) :: V(nb, nb, nk)
    real(8) :: s1, s2, tp
    tp = 8d0 * atan(1d0)
    num_kgrid = (/ n1, n2, 1 /)
    bvec(:, 1) = (/ 1, 0, 0 /)
    bvec(:, 2) = (/ 0, 1, 0 /)
    bvec(:, 3) = (/ 0, 0, 1 /)

    do ik = 1, nk
      i1 = mod(ik - 1, n1);  i2 = (ik - 1) / n1
      s1 = tp * dble(i1) / dble(n1);  s2 = tp * dble(i2) / dble(n2)
      ! a smooth SU(nb)-ish frame: a real rotation mixing (1,2) and (3,4) plus phases
      V(:, :, ik) = (0d0, 0d0)
      do i = 1, nb
        V(i, i, ik) = exp(cmplx(0d0, 0.31d0 * sin(s1) + 0.17d0 * dble(i) * cos(s2), 8))
      end do
      V(1, 2, ik) = V(1, 2, ik) + 0.28d0 * sin(s2)
      V(2, 1, ik) = V(2, 1, ik) - 0.28d0 * sin(s2)
      V(3, 4, ik) = V(3, 4, ik) + 0.22d0 * cos(s1)
      V(4, 3, ik) = V(4, 3, ik) - 0.22d0 * cos(s1)
      ! orthonormalize the columns (modified Gram-Schmidt) -> exactly unitary
      call mgs(V(:, :, ik), nb)
      do i = 1, nb
        do j = 1, nb
          rho(i, j, ik) = cmplx(0.3d0 * cos(s1 + 0.7d0 * dble(i)), &
                              & 0.2d0 * sin(s2 + 0.5d0 * dble(j)), 8)
        end do
      end do
    end do

    do ik = 1, nk
      i1 = mod(ik - 1, n1);  i2 = (ik - 1) / n1
      do iv = 1, nbvec
        j1 = modulo(i1 + bvec(1, iv), n1)
        j2 = modulo(i2 + bvec(2, iv), n2)
        jk = j2 * n1 + j1 + 1
        prod_dk(:, :, iv, ik) = matmul(conjg(transpose(V(:, :, ik))), V(:, :, jk))
      end do
    end do
  end subroutine make_fixture

  subroutine mgs(A, n)
    integer, intent(in) :: n
    complex(8), intent(inout) :: A(n, n)
    integer :: i, j
    complex(8) :: p
    real(8) :: nr
    do j = 1, n
      do i = 1, j - 1
        p = sum(conjg(A(:, i)) * A(:, j))
        A(:, j) = A(:, j) - p * A(:, i)
      end do
      nr = sqrt(dble(sum(conjg(A(:, j)) * A(:, j))))
      A(:, j) = A(:, j) / nr
    end do
  end subroutine mgs

  !===================================================================
  ! A. build_block_transport on a SUB-RANGE == the full-range result there.
  !===================================================================
  subroutine test_transport_klocal(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 4, n1 = 6, n2 = 5, nk = n1 * n2, nbvec = 3
    integer :: bvec(3, nbvec), num_kgrid(3), block_id(nb, nk)
    complex(8) :: prod_dk(nb, nb, nbvec, nk), rho(nb, nb, nk)
    complex(8) :: Ufull(nb, nb, 3, nk)
    complex(8), allocatable :: Uloc(:, :, :, :)
    integer :: nrej, lo, hi, ik, nbad

    call make_fixture(nb, nk, n1, n2, nbvec, bvec, num_kgrid, prod_dk, rho)
    block_id(:, :) = 1                                  ! X-full: one full-band block

    call build_block_transport(nb, nk, nbvec, bvec, prod_dk, block_id, num_kgrid, &
                             & Ufull, nrej, 1, nk)

    ! three disjoint slices that tile 1..nk, exactly as split_range would hand out
    nbad = 0
    do lo = 1, nk, 11
      hi = min(nk, lo + 10)
      allocate(Uloc(nb, nb, 3, lo:hi))
      call build_block_transport(nb, nk, nbvec, bvec, prod_dk(:, :, :, lo:hi), block_id, &
                               & num_kgrid, Uloc, nrej, lo, hi)
      do ik = lo, hi
        if (any(Uloc(:, :, :, ik) /= Ufull(:, :, :, ik))) nbad = nbad + 1
      end do
      deallocate(Uloc)
    end do
    call check(nbad == 0, "A: build_block_transport on a k-slice is BIT-identical to the full-range run", nfail)
  end subroutine test_transport_klocal

  !===================================================================
  ! B. covariant_grad_block is invariant under the slot map (kmap).
  !===================================================================
  subroutine test_grad_slot_invariance(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb = 4, n1 = 7, n2 = 6, nk = n1 * n2, nbvec = 3
    integer :: bvec(3, nbvec), num_kgrid(3), block_id(nb, nk)
    complex(8) :: prod_dk(nb, nb, nbvec, nk), rho(nb, nb, nk)
    complex(8) :: Ufull(nb, nb, 3, nk), Dq_ref(nb, nb, 3, nk)
    complex(8), allocatable :: Us(:, :, :, :), rs(:, :, :), Dq_s(:, :, :, :)
    integer :: kmap(nk), nrej, ik, lo, hi, ns, nbad, i
    logical :: needed(nk)
    real(8) :: dk(3)

    call make_fixture(nb, nk, n1, n2, nbvec, bvec, num_kgrid, prod_dk, rho)
    block_id(:, :) = 1
    dk = (/ 0.31d0, 0.27d0, 1d0 /)
    call build_block_transport(nb, nk, nbvec, bvec, prod_dk, block_id, num_kgrid, &
                             & Ufull, nrej, 1, nk)

    ! reference: replicated layout, identity slot map (= the pre-klocal code)
    call covariant_grad_block(nb, nk, nbvec, bvec, num_kgrid, Ufull, rho, dk, Dq_ref, &
                            & 1, nk, (/ .true., .true., .true. /), nk, identity_kmap(nk))

    ! k-local layout of a middle "rank": owned [lo,hi] + the stencil's own halo,
    ! packed into slots in a SCRAMBLED order (reverse) so a slot-arithmetic slip
    ! cannot cancel.
    lo = 15;  hi = 27
    call covariant_halo_needed(nk, nbvec, bvec, num_kgrid, lo, hi, needed)
    ns = count(needed)
    kmap(:) = 0
    i = ns
    do ik = 1, nk
      if (needed(ik)) then
        kmap(ik) = i        ! reverse packing: slot ns, ns-1, ...
        i = i - 1
      end if
    end do
    allocate(Us(nb, nb, 3, ns), rs(nb, nb, ns), Dq_s(nb, nb, 3, lo:hi))
    Us = (0d0, 0d0);  rs = (0d0, 0d0)
    do ik = 1, nk
      if (kmap(ik) > 0) then
        Us(:, :, :, kmap(ik)) = Ufull(:, :, :, ik)
        rs(:, :, kmap(ik))    = rho(:, :, ik)
      end if
    end do

    call covariant_grad_block(nb, nk, nbvec, bvec, num_kgrid, Us, rs, dk, Dq_s, &
                            & lo, hi, (/ .true., .true., .true. /), ns, kmap)

    nbad = 0
    do ik = lo, hi
      if (any(Dq_s(:, :, :, ik) /= Dq_ref(:, :, :, ik))) nbad = nbad + 1
    end do
    call check(nbad == 0, &
      & "B: covariant_grad_block over a scrambled owned+halo slot layout is BIT-identical", nfail)
    call check(ns < nk, "B: the k-local layout really does store fewer k than the full mesh", nfail)
    deallocate(Us, rs, Dq_s)
  end subroutine test_grad_slot_invariance

  subroutine check(ok, msg, nfail)
    logical, intent(in) :: ok
    character(*), intent(in) :: msg
    integer, intent(inout) :: nfail
    if (ok) then
      write(*, '(a,a)') "PASS ", msg
    else
      write(*, '(a,a)') "FAIL ", msg
      nfail = nfail + 1
    end if
  end subroutine check

end program test_gicov_klocal
