!
!  Copyright 2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!-----------------------------------------------------------------------------------------
! Joint (simultaneous) diagonalization of one or more mutually commuting
! unitary matrices, with a genuinely orthonormal eigenvector matrix even
! when eigenvalues are degenerate.
!
! A general (non-Hermitian) complex eigensolver does not guarantee that the
! eigenvectors it returns for a repeated eigenvalue are mutually orthogonal;
! it cannot be relied on to produce a unitary eigenvector matrix whenever the
! input has degenerate eigenvalues. A unitary matrix u is normal
! (u^dagger u = u u^dagger = I), so an orthonormal eigenbasis does exist in
! principle; the standard way to obtain it numerically, reusing only an
! ordinary Hermitian eigensolver (for which degenerate eigenvectors ARE
! guaranteed mutually orthogonal), is to split u into Hermitian "real" and
! "imaginary" parts,
!
!   a = (u + u^dagger)/2,   b = (u - u^dagger)/(2i),   u = a + i*b,
!
! which commute with each other (and with anything u itself commutes with,
! such as a Hamiltonian of which u represents a symmetry). Diagonalizing a
! first, and then, only within whatever degenerate block of a remains,
! diagonalizing the restriction of b to that block, recovers a genuine joint
! eigenbasis of a and b, hence of u = a + i*b. Because a's eigenvalue is
! cos(theta) and b's is sin(theta) for a phase e^{i theta}, the second step
! is exactly what is needed to distinguish phases that share the same real
! part (theta and -theta). Any degeneracy still remaining after both steps
! is a genuine (not merely apparent) multiple eigenvalue of u, for which any
! orthonormal basis of that eigenspace is an equally valid choice.
!
! diagonalize_commuting_unitary_family repeats this for a whole family of
! mutually commuting unitary matrices, refining only whichever blocks are
! still left degenerate by the matrices already used -- exactly what is
! needed when a single member of the family does not, by itself, separate
! every vector.
!-----------------------------------------------------------------------------------------
module eigen_unitary_sub
  use eigen_lapack, only: eigen_zheev
  use math_constants, only: zi
  implicit none

  private
  public :: diagonalize_unitary, diagonalize_commuting_unitary_family

contains

  !> Diagonalize a single n x n unitary matrix u.
  !>
  !> w(:,i)     : the i-th eigenvector (columns of w form a unitary matrix,
  !>              i.e. are exactly orthonormal, to machine precision, even
  !>              when u has degenerate eigenvalues).
  !> theta(i)   : the phase (radians, in (-pi,pi]) of the eigenvalue for
  !>              w(:,i), i.e. u.w(:,i) = exp(i*theta(i)) * w(:,i), read off
  !>              directly from w(:,i)^dagger . u . w(:,i). This alone is
  !>              not a check of how unitary u itself was: a vector that is
  !>              far from being an eigenvector can still give
  !>              |w^dagger u w| close to 1, and a genuine eigenvector of a
  !>              slightly non-unitary u can give it close to 1 too, so it
  !>              conflates whatever error u itself carries with whatever
  !>              residual this routine's diagonalization leaves. Checking
  !>              u's own unitarity (u^dagger u - I), the orthonormality of
  !>              w (w^dagger w - I), and, for a family, the off-diagonal
  !>              size of w^dagger u_k w, are the caller's responsibility
  !>              where needed.
  !> block_id(i): groups, for this particular u and this particular tol,
  !>              the columns that could not be separated within that
  !>              tolerance. For an (exactly or near-exactly) commuting
  !>              family, resolved at a tolerance fine enough to tell the
  !>              true eigenvalues apart, a shared block_id corresponds to
  !>              a genuine common eigenspace; with a coarser tol, distinct
  !>              eigenvalues closer together than tol can also end up
  !>              sharing a block_id -- that is tol doing its job, not a
  !>              bug. Any orthonormal basis of a given block is an equally
  !>              valid choice, and w's particular choice within it carries
  !>              no extra meaning. block_id is assigned in increasing
  !>              blocks of contiguous columns; it is not itself an
  !>              eigenvalue ordering.
  !> tol (optional): the *numerical* tolerance used to decide whether two
  !>              eigenvalues of the intermediate Hermitian matrices a, b
  !>              are "the same" for blocking purposes. This has nothing to
  !>              do with any physical energy scale; it defaults to a small
  !>              fixed multiple of machine epsilon.
  subroutine diagonalize_unitary(n, u, w, theta, block_id, tol)
    implicit none
    integer,    intent(in)  :: n
    complex(8), intent(in)  :: u(n,n)
    complex(8), intent(out) :: w(n,n)
    real(8),    intent(out) :: theta(n)
    integer,    intent(out) :: block_id(n)
    real(8),    intent(in), optional :: tol
    complex(8) :: a(n,n), b(n,n)
    real(8)    :: eps, ea(n)
    complex(8) :: va(n,n), uw
    complex(8), allocatable :: bblk(:,:), vb(:,:)
    real(8),    allocatable :: eb(:)
    integer :: i, j, i0, i1, j0, j1, nb, m

    eps = 1.0d2 * epsilon(1d0)
    if (present(tol)) eps = tol

    ! a, b: Hermitian real/imaginary parts of u (u = a + i*b); both commute
    ! with each other and with u itself.
    a = 0.5d0 * ( u + conjg(transpose(u)) )
    b = -0.5d0 * zi * ( u - conjg(transpose(u)) )

    call eigen_zheev(a, ea, va)   ! ea ascending; va exactly unitary

    w = va
    block_id = 0
    nb = 0
    i0 = 1
    do while (i0 <= n)
      i1 = i0
      do while (i1 < n)
        if ( abs(ea(i1+1)-ea(i0)) > eps ) exit
        i1 = i1 + 1
      end do
      ! columns i0..i1 share the same eigenvalue of a
      if (i1 > i0) then
        m = i1 - i0 + 1
        allocate( bblk(m,m), vb(m,m), eb(m) )
        bblk = matmul( conjg(transpose(va(:,i0:i1))), matmul(b,va(:,i0:i1)) )
        call eigen_zheev(bblk, eb, vb)   ! eb ascending; vb exactly unitary
        w(:,i0:i1) = matmul(va(:,i0:i1), vb)
        ! sub-group columns i0..i1 by (now resolved) eigenvalue of b;
        ! anything still tied here is a genuine degeneracy of u itself
        j0 = 1
        do while (j0 <= m)
          j1 = j0
          do while (j1 < m)
            if ( abs(eb(j1+1)-eb(j0)) > eps ) exit
            j1 = j1 + 1
          end do
          nb = nb + 1
          do j = j0, j1
            block_id(i0+j-1) = nb
          end do
          j0 = j1 + 1
        end do
        deallocate( bblk, vb, eb )
      else
        nb = nb + 1
        block_id(i0) = nb
      end if
      i0 = i1 + 1
    end do

    ! phase of each final eigenvector, read off directly from u
    do i = 1, n
      uw = dot_product( w(:,i), matmul(u,w(:,i)) )
      theta(i) = atan2( aimag(uw), real(uw,kind=8) )
    end do

  end subroutine diagonalize_unitary

  !> Joint (simultaneous) diagonalization of m mutually commuting n x n
  !> unitary matrices ulist(:,:,1:m). Refines the same unitary eigenvector
  !> matrix w, member by member, using diagonalize_unitary within whichever
  !> blocks are still unresolved by the members already used; blocks
  !> resolved down to a single column by an earlier member are left
  !> untouched by later ones. Any block still of size > 1 once all m
  !> members have been used shares the same block_id under every member of
  !> the family used so far, to within tol; for an (exactly or
  !> near-exactly) commuting family, resolved at a tolerance fine enough to
  !> tell the true eigenvalues apart, this is a genuine common eigenspace
  !> rather than an order-dependent artifact. Whether "same block_id" can
  !> further be read as "same sector" (or whatever the caller's blocks are
  !> meant to represent) is the caller's responsibility to establish, not a
  !> guarantee this routine makes on its own: it holds only if ulist is
  !> actually a complete enough set to separate every sector that needs
  !> separating, and only if the numerical separation in fact succeeds.
  !>
  !> w, block_id, tol: as in diagonalize_unitary, but for the whole family;
  !> block_id groups columns by their final, common resolution (again in
  !> increasing blocks of contiguous columns).
  subroutine diagonalize_commuting_unitary_family(n, m, ulist, w, block_id, tol)
    implicit none
    integer,    intent(in)  :: n, m
    complex(8), intent(in)  :: ulist(n,n,m)
    complex(8), intent(out) :: w(n,n)
    integer,    intent(out) :: block_id(n)
    real(8),    intent(in), optional :: tol
    real(8) :: eps
    ! bnd0(ib):bnd1(ib) is the (inclusive) column range of the ib-th current
    ! block, in the present column ordering of w; at most n blocks can ever
    ! exist, so fixed-size arrays of length n are always large enough.
    integer :: bnd0(n), bnd1(n), nblk
    integer :: bnd0_new(n), bnd1_new(n), nblk_new
    complex(8), allocatable :: ublk(:,:), wblk(:,:)
    real(8),    allocatable :: thblk(:)
    integer,    allocatable :: bidblk(:)
    integer :: k, ib, i0, i1, sz, j, jj

    eps = 1.0d2 * epsilon(1d0)
    if (present(tol)) eps = tol

    w = (0d0,0d0)
    do k = 1, n
      w(k,k) = (1d0,0d0)
    end do

    nblk = 1
    bnd0(1) = 1
    bnd1(1) = n

    do k = 1, m
      if (nblk == n) exit   ! every column already fully resolved

      nblk_new = 0
      do ib = 1, nblk
        i0 = bnd0(ib)
        i1 = bnd1(ib)
        sz = i1 - i0 + 1

        if (sz == 1) then
          nblk_new = nblk_new + 1
          bnd0_new(nblk_new) = i0
          bnd1_new(nblk_new) = i1
          cycle
        end if

        allocate( ublk(sz,sz), wblk(sz,sz), thblk(sz), bidblk(sz) )
        ublk = matmul( conjg(transpose(w(:,i0:i1))), matmul(ulist(:,:,k),w(:,i0:i1)) )
        call diagonalize_unitary(sz, ublk, wblk, thblk, bidblk, eps)
        w(:,i0:i1) = matmul(w(:,i0:i1), wblk)

        ! bidblk is assigned in contiguous, increasing groups (see
        ! diagonalize_unitary above), so this splits [i0,i1] into
        ! contiguous sub-ranges without needing to touch column order.
        j = 1
        do while (j <= sz)
          jj = j
          do while (jj < sz)
            if (bidblk(jj+1) /= bidblk(j)) exit
            jj = jj + 1
          end do
          nblk_new = nblk_new + 1
          bnd0_new(nblk_new) = i0 + j  - 1
          bnd1_new(nblk_new) = i0 + jj - 1
          j = jj + 1
        end do

        deallocate( ublk, wblk, thblk, bidblk )
      end do

      nblk = nblk_new
      bnd0(1:nblk) = bnd0_new(1:nblk)
      bnd1(1:nblk) = bnd1_new(1:nblk)
    end do

    block_id = 0
    do ib = 1, nblk
      do j = bnd0(ib), bnd1(ib)
        block_id(j) = ib
      end do
    end do

  end subroutine diagonalize_commuting_unitary_family

end module eigen_unitary_sub
