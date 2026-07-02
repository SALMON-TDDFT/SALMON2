! src/ssbe/test/test_block_transport.f90
! LG-SBE Tier2 Phase 1 -- unit tests for the polar (Wilson-line) parallel-
! transport kernel build_block_transport / polar_unitary
! (degenerate_block_ssbe). This is the branch-cut-free replacement for the
! logm half of xi_block_from_overlap: U = M (M^H M)^{-1/2} (the polar
! factor) ONLY -- never logm -- so no eigenphase branch cut
! (xi_block_from_overlap's info=2 case) can occur here by construction.
!
! Test A proves: (i) a link that would trip build_xi's 0.9*pi branch-cut
! guard (a near-pi block rotation) builds with ZERO rejects under the
! polar-only transport, and (ii) the resulting U is a gauge-covariant
! tensor  Ut(k) = W(k)^H U(k) W(k+e)  under a stencil-VARYING U(N) gauge
! (a stronger covariance claim than build_xi's stencil-constant test,
! possible because the polar decomposition of a well-conditioned block is
! UNIQUE: A^H M B = (A^H U B)(B^H P B) is itself a valid polar
! decomposition whenever A, B are unitary and M is invertible).
!
! Test B unit-tests polar_unitary directly: a near-singular block rejects
! with ierr=1 (no process abort -- only build_block_transport itself
! aborts, via error stop, on ierr/=0), and a well-conditioned unitary
! block gives ierr=0 and an exactly unitary U.
!
! Standalone build from src/ssbe (macOS: Accelerate provides zheev/zgemm):
!     gfortran degenerate_block_ssbe.f90 test/test_block_transport.f90 \
!              -o t -framework Accelerate && ./t

program test_block_transport
  use degenerate_block_ssbe, only: build_block_transport, polar_unitary
  implicit none
  complex(8), parameter :: zi = (0d0, 1d0)
  integer :: nfail
  nfail = 0

  call test_transport_unitary_and_covariant(nfail)
  call test_transport_singular_failclosed(nfail)

  if (nfail > 0) then
    write(*, '(a,i0,a)') "FAILED: ", nfail, " check(s)"
    stop 1
  else
    write(*, '(a)') "All test_block_transport checks passed."
  end if

contains

  !======================= assert helpers (ssbe style, copied verbatim
  !======================= from test/test_gisbe_covariance.f90) =========
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

  !======================= small matrix utilities (copied verbatim from
  !======================= test/test_gisbe_covariance.f90) ==============
  function hc(A) result(B)                       ! Hermitian conjugate A^H
    complex(8), intent(in) :: A(:, :)
    complex(8), allocatable :: B(:, :)
    B = conjg(transpose(A))
  end function hc

  function rot(W, A) result(R)                   ! covariant rotation W^H A W
    complex(8), intent(in) :: W(:, :), A(:, :)
    complex(8), allocatable :: R(:, :)
    R = matmul(matmul(conjg(transpose(W)), A), W)
  end function rot

  real(8) function mdiff(A, B)                    ! max_ij |A-B|
    complex(8), intent(in) :: A(:, :), B(:, :)
    mdiff = maxval(abs(A - B))
  end function mdiff

  real(8) function mmax(A)                        ! max_ij |A|
    complex(8), intent(in) :: A(:, :)
    mmax = maxval(abs(A))
  end function mmax

  ! General U(2) gauge  e^{i chi} [[e^{i al}c, e^{i be}s],[-e^{-i be}s, e^{-i al}c]].
  function u2_gauge(al, be, th, chi) result(U)
    real(8), intent(in) :: al, be, th, chi
    complex(8), allocatable :: U(:, :)
    real(8) :: c, s
    complex(8) :: g
    allocate(U(2, 2))
    c = cos(th); s = sin(th); g = exp(zi * chi)
    U(1, 1) =  g * exp( zi * al) * c
    U(1, 2) =  g * exp( zi * be) * s
    U(2, 1) = -g * exp(-zi * be) * s
    U(2, 2) =  g * exp(-zi * al) * c
  end function u2_gauge

  ! Embed a d x d block gauge W at bands i0..i0+d-1 of an nb x nb identity
  ! (gauge acts inside the degenerate block, identity on all other bands).
  function embed(nb, i0, W) result(V)
    integer, intent(in) :: nb, i0
    complex(8), intent(in) :: W(:, :)
    complex(8), allocatable :: V(:, :)
    integer :: a, b, d
    d = size(W, 1)
    allocate(V(nb, nb)); V = (0d0, 0d0)
    do a = 1, nb
      V(a, a) = (1d0, 0d0)
    end do
    do b = 1, d
      do a = 1, d
        V(i0 + a - 1, i0 + b - 1) = W(a, b)
      end do
    end do
  end function embed

  ! nb x nb frame for the build_xi/build_block_transport model: identity on
  ! bands 1 and nb, a real 2x2 rotation R(ang) on the exactly degenerate
  ! block {2,3}.
  function frame4(nb, ang) result(Q)
    integer, intent(in) :: nb
    real(8), intent(in) :: ang
    complex(8), allocatable :: Q(:, :)
    integer :: a
    real(8) :: c, s
    allocate(Q(nb, nb)); Q = (0d0, 0d0)
    do a = 1, nb
      Q(a, a) = (1d0, 0d0)
    end do
    c = cos(ang); s = sin(ang)
    Q(2, 2) = dcmplx(c, 0d0); Q(2, 3) = dcmplx(-s, 0d0)
    Q(3, 2) = dcmplx(s, 0d0); Q(3, 3) = dcmplx(c, 0d0)
  end function frame4

  ! n x n complex identity (not present in test_gisbe_covariance.f90).
  function eye(n) result(I_)
    integer, intent(in) :: n
    complex(8), allocatable :: I_(:, :)
    integer :: a
    allocate(I_(n, n)); I_ = (0d0, 0d0)
    do a = 1, n
      I_(a, a) = (1d0, 0d0)
    end do
  end function eye

  !================================ tests ====================================

  ! Test A -- the branch-cut-gone + covariance proof.
  subroutine test_transport_unitary_and_covariant(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nb=4, nk=4, nbvec=3
    integer :: num_kgrid(3), bvec(3,nbvec), ik, ikx, a, nrej, nrej_t
    complex(8) :: prod_dk(nb,nb,nbvec,nk), prod_t(nb,nb,nbvec,nk)
    complex(8) :: U(nb,nb,3,nk), Ut(nb,nb,3,nk), Wk(nb,nb,nk)
    real(8)    :: eigen(nb,nk), worst, uni
    complex(8), allocatable :: Wb(:,:)
    num_kgrid=(/nk,1,1/); bvec(:,1)=(/1,0,0/); bvec(:,2)=(/0,1,0/); bvec(:,3)=(/0,0,1/)
    eigen(1,:)=0.36d0; eigen(2,:)=0.30d0; eigen(3,:)=0.30d0; eigen(4,:)=0.54d0
    prod_dk=(0d0,0d0)
    do ik=1,nk
      ! near-pi (0.95*pi) real rotation on the {2,3} block along +x -> would trip build_xi's 0.9pi reject
      prod_dk(:,:,1,ik) = matmul(hc(frame4(nb,0d0)), frame4(nb,0.95d0*acos(-1d0)))
      do a=1,nb; prod_dk(a,a,2,ik)=(1d0,0d0); prod_dk(a,a,3,ik)=(1d0,0d0); end do
    end do
    call build_block_transport(nb,nk,nbvec,bvec,prod_dk,eigen,num_kgrid,U,nrej)
    call check_int(nrej,0,"transport: near-pi block builds with zero rejects (no branch-cut)",nfail)
    uni=0d0
    do ik=1,nk; uni=max(uni, mdiff(matmul(hc(U(:,:,1,ik)),U(:,:,1,ik)), eye(nb))); end do
    call check_close_r(uni,0d0,1d-12,"transport: U^H U = I on +x link",nfail)
    do ik=1,nk; Wb=u2_gauge(0.3d0*ik,0.2d0,0.7d0,-0.1d0*ik); Wk(:,:,ik)=embed(nb,2,Wb); end do
    do ik=1,nk; ikx=mod(ik,nk)+1
      prod_t(:,:,1,ik)=matmul(matmul(hc(Wk(:,:,ik)),prod_dk(:,:,1,ik)),Wk(:,:,ikx))
      do a=1,nb; prod_t(a,a,2,ik)=(1d0,0d0); prod_t(a,a,3,ik)=(1d0,0d0); end do
    end do
    call build_block_transport(nb,nk,nbvec,bvec,prod_t,eigen,num_kgrid,Ut,nrej_t)
    call check_int(nrej_t,0,"transport(gauged): zero rejects",nfail)
    worst=0d0
    do ik=1,nk; ikx=mod(ik,nk)+1
      worst=max(worst, mdiff(Ut(:,:,1,ik), matmul(matmul(hc(Wk(:,:,ik)),U(:,:,1,ik)),Wk(:,:,ikx))))
    end do
    call check_close_r(worst,0d0,1d-10,"transport: U covariant  Ut = W(k)^H U W(k+e)",nfail)
  end subroutine test_transport_unitary_and_covariant

  ! Test B -- polar_unitary fail-closed via ierr (no process abort): a
  ! near-singular block (two nearly-parallel columns) must reject with
  ! ierr==1, and a well-conditioned unitary block must give ierr==0 with an
  ! exactly unitary U. build_block_transport itself is NOT called here (it
  ! would error stop on a bad block).
  subroutine test_transport_singular_failclosed(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: d = 2
    complex(8) :: M(d,d), U(d,d)
    complex(8), allocatable :: Muni(:,:)
    real(8) :: sigma_min, sigma_ok
    integer :: ierr, ierr_ok

    ! near-singular: two columns nearly parallel
    M(1,1) = (1d0,0d0); M(1,2) = (1d0,0d0)
    M(2,1) = (1d0,0d0); M(2,2) = dcmplx(1d0+1d-14,0d0)
    call polar_unitary(M, d, U, sigma_min, ierr)
    call check_int(ierr, 1, "polar_unitary: near-singular block rejected (ierr=1)", nfail)

    ! well-conditioned unitary M -> ierr=0, U^H U = I
    Muni = u2_gauge(0.3d0, -0.5d0, 0.6d0, 0.2d0)
    call polar_unitary(Muni, d, U, sigma_ok, ierr_ok)
    call check_int(ierr_ok, 0, "polar_unitary: well-conditioned unitary M -> ierr=0", nfail)
    call check_close_r(mdiff(matmul(hc(U),U), eye(d)), 0d0, 1d-12, &
                       "polar_unitary: well-conditioned -> U^H U = I", nfail)
  end subroutine test_transport_singular_failclosed

end program test_block_transport
