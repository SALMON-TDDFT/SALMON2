!
!  Copyright 2019-2026 SALMON developers
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
!=======================================================================
! WPW reduced-basis linear algebra helpers
!=======================================================================

module rt_dg_wpw_linalg
  implicit none

  private
  public :: build_wpw_sorth_reduced_neighbor_block
  public :: apply_s_orthogonal_transform
  public :: build_symmetric_orthogonal_hamiltonian
  public :: build_hermitian_pseudoinverse
  public :: build_hermitian_inverse_sqrt
  public :: build_hermitian_inverse
  public :: build_wpw_reduced_c_can_reference
  public :: build_wpw_reduced_raw_back_hybrid_from_inverse
  public :: build_wpw_reduced_raw_back_hybrid_matrix
  public :: zheev_with_query
  public :: zhegv_with_query
  public :: rectangular_singular_minmax
  public :: wpw_local_herm_max

contains

  subroutine build_wpw_sorth_reduced_neighbor_block(H_in, S_in, n_self, n_ext, tol, H_red, S_red, &
      n_red, n_keep, n_drop, lambda_keep_min, lambda_keep_max, s_sn_after, s_nn_identity_err, &
      sss_min, sss_max, info, T_ext_red)
    implicit none
    integer, intent(in) :: n_self, n_ext
    complex(8), intent(in) :: H_in(n_ext,n_ext), S_in(n_ext,n_ext)
    real(8), intent(in) :: tol
    complex(8), allocatable, intent(out) :: H_red(:,:), S_red(:,:)
    complex(8), allocatable, intent(out), optional :: T_ext_red(:,:)
    integer, intent(out) :: n_red, n_keep, n_drop, info
    real(8), intent(out) :: lambda_keep_min, lambda_keep_max, s_sn_after, s_nn_identity_err, sss_min, sss_max

    integer :: i, j, n_neigh
    complex(8), allocatable :: H_sorth(:,:), S_sorth(:,:), sss_inv(:,:), xmat(:,:)
    complex(8), allocatable :: snn_vec(:,:), ymat(:,:), tmp(:,:)
    real(8), allocatable :: lambda(:)
    integer, allocatable :: keep_idx(:)

    n_red = n_self
    n_keep = 0
    n_drop = 0
    lambda_keep_min = 0.0d0
    lambda_keep_max = 0.0d0
    s_sn_after = 0.0d0
    s_nn_identity_err = 0.0d0
    sss_min = 0.0d0
    sss_max = 0.0d0
    info = 0
    if (allocated(H_red)) deallocate(H_red)
    if (allocated(S_red)) deallocate(S_red)
    if (present(T_ext_red)) then
      if (allocated(T_ext_red)) deallocate(T_ext_red)
    end if
    if (n_ext <= n_self .or. n_self <= 0) then
      info = -1
      return
    end if
    n_neigh = n_ext - n_self
    allocate(H_sorth(n_ext,n_ext), S_sorth(n_ext,n_ext), sss_inv(n_self,n_self), xmat(n_self,n_neigh))
    allocate(snn_vec(n_neigh,n_neigh), lambda(n_neigh), keep_idx(n_neigh))

    call build_hermitian_inverse(S_in(1:n_self, 1:n_self), n_self, sss_inv, info, sss_min, sss_max)
    if (info /= 0) then
      deallocate(H_sorth, S_sorth, sss_inv, xmat, snn_vec, lambda, keep_idx)
      return
    end if

    xmat(:, :) = matmul(sss_inv, S_in(1:n_self, n_self+1:n_ext))
    call apply_s_orthogonal_transform(H_in, xmat, n_self, n_ext, H_sorth)
    call apply_s_orthogonal_transform(S_in, xmat, n_self, n_ext, S_sorth)
    call hermitize_matrix(H_sorth, n_ext)
    call hermitize_matrix(S_sorth, n_ext)

    snn_vec(:, :) = S_sorth(n_self+1:n_ext, n_self+1:n_ext)
    call zheev_with_query(snn_vec, n_neigh, lambda, info)
    if (info /= 0) then
      deallocate(H_sorth, S_sorth, sss_inv, xmat, snn_vec, lambda, keep_idx)
      return
    end if

    n_keep = 0
    do i = 1, n_neigh
      if (lambda(i) > tol) then
        n_keep = n_keep + 1
        keep_idx(n_keep) = i
      end if
    end do
    n_drop = n_neigh - n_keep
    if (n_keep <= 0) then
      info = -2
      deallocate(H_sorth, S_sorth, sss_inv, xmat, snn_vec, lambda, keep_idx)
      return
    end if

    n_red = n_self + n_keep
    allocate(ymat(n_neigh,n_keep), tmp(n_neigh,n_keep), H_red(n_red,n_red), S_red(n_red,n_red))
    if (present(T_ext_red)) allocate(T_ext_red(n_ext, n_red))
    ymat(:, :) = (0.0d0, 0.0d0)
    do j = 1, n_keep
      i = keep_idx(j)
      ymat(:, j) = snn_vec(:, i) / sqrt(lambda(i))
    end do

    H_red(:, :) = (0.0d0, 0.0d0)
    S_red(:, :) = (0.0d0, 0.0d0)
    if (present(T_ext_red)) then
      T_ext_red(:, :) = (0.0d0, 0.0d0)
      do i = 1, n_self
        T_ext_red(i, i) = (1.0d0, 0.0d0)
      end do
      T_ext_red(1:n_self, n_self+1:n_red) = -matmul(xmat, ymat)
      T_ext_red(n_self+1:n_ext, n_self+1:n_red) = ymat(1:n_neigh, 1:n_keep)
    end if
    H_red(1:n_self, 1:n_self) = H_sorth(1:n_self, 1:n_self)
    S_red(1:n_self, 1:n_self) = S_sorth(1:n_self, 1:n_self)
    H_red(1:n_self, n_self+1:n_red) = matmul(H_sorth(1:n_self, n_self+1:n_ext), ymat)
    S_red(1:n_self, n_self+1:n_red) = matmul(S_sorth(1:n_self, n_self+1:n_ext), ymat)
    H_red(n_self+1:n_red, 1:n_self) = conjg(transpose(H_red(1:n_self, n_self+1:n_red)))
    S_red(n_self+1:n_red, 1:n_self) = conjg(transpose(S_red(1:n_self, n_self+1:n_red)))
    tmp(:, :) = matmul(H_sorth(n_self+1:n_ext, n_self+1:n_ext), ymat)
    H_red(n_self+1:n_red, n_self+1:n_red) = matmul(conjg(transpose(ymat)), tmp)
    tmp(:, :) = matmul(S_sorth(n_self+1:n_ext, n_self+1:n_ext), ymat)
    S_red(n_self+1:n_red, n_self+1:n_red) = matmul(conjg(transpose(ymat)), tmp)
    call hermitize_matrix(H_red, n_red)
    call hermitize_matrix(S_red, n_red)

    s_sn_after = sqrt(sum(abs(S_red(1:n_self, n_self+1:n_red))**2))
    s_nn_identity_err = 0.0d0
    do i = 1, n_keep
      do j = 1, n_keep
        if (i == j) then
          s_nn_identity_err = max(s_nn_identity_err, abs(S_red(n_self+i, n_self+j) - (1.0d0, 0.0d0)))
        else
          s_nn_identity_err = max(s_nn_identity_err, abs(S_red(n_self+i, n_self+j)))
        end if
      end do
    end do
    s_sn_after = sqrt(sum(abs(S_red(1:n_self, n_self+1:n_red))**2))
    lambda_keep_min = lambda(keep_idx(1))
    lambda_keep_max = lambda(keep_idx(n_keep))

    deallocate(H_sorth, S_sorth, sss_inv, xmat, snn_vec, lambda, keep_idx, ymat, tmp)
  end subroutine build_wpw_sorth_reduced_neighbor_block


  subroutine apply_s_orthogonal_transform(M_in, xmat, n_self, n_ext, M_out)
    implicit none
    integer, intent(in) :: n_self, n_ext
    complex(8), intent(in) :: M_in(n_ext,n_ext), xmat(n_self,n_ext-n_self)
    complex(8), intent(out) :: M_out(n_ext,n_ext)

    integer :: n_neigh
    complex(8), allocatable :: mss(:,:), msn(:,:), mns(:,:), mnn(:,:)

    n_neigh = n_ext - n_self
    allocate(mss(n_self,n_self), msn(n_self,n_neigh), mns(n_neigh,n_self), mnn(n_neigh,n_neigh))
    mss(:, :) = M_in(1:n_self, 1:n_self)
    msn(:, :) = M_in(1:n_self, n_self+1:n_ext)
    mns(:, :) = M_in(n_self+1:n_ext, 1:n_self)
    mnn(:, :) = M_in(n_self+1:n_ext, n_self+1:n_ext)

    M_out(:, :) = (0.0d0, 0.0d0)
    M_out(1:n_self, 1:n_self) = mss(:, :)
    M_out(1:n_self, n_self+1:n_ext) = msn(:, :) - matmul(mss, xmat)
    M_out(n_self+1:n_ext, 1:n_self) = mns(:, :) - matmul(conjg(transpose(xmat)), mss)
    M_out(n_self+1:n_ext, n_self+1:n_ext) = mnn(:, :) - matmul(mns, xmat) - &
      matmul(conjg(transpose(xmat)), msn) + matmul(matmul(conjg(transpose(xmat)), mss), xmat)
    deallocate(mss, msn, mns, mnn)
  end subroutine apply_s_orthogonal_transform


  subroutine build_symmetric_orthogonal_hamiltonian(H_in, S_in, n, H_orth, info, s_min, s_max)
    implicit none
    integer, intent(in) :: n
    complex(8), intent(in) :: H_in(n,n), S_in(n,n)
    complex(8), intent(out) :: H_orth(n,n)
    integer, intent(out) :: info
    real(8), intent(out) :: s_min, s_max

    integer :: i
    complex(8), allocatable :: svec(:,:), xmat(:,:), tmp(:,:)
    real(8), allocatable :: eval(:)

    allocate(svec(n,n), xmat(n,n), tmp(n,n), eval(n))
    svec(:, :) = S_in(:, :)
    call zheev_with_query(svec, n, eval, info)
    if (info /= 0) then
      H_orth(:, :) = (0.0d0, 0.0d0)
      s_min = 0.0d0
      s_max = 0.0d0
      deallocate(svec, xmat, tmp, eval)
      return
    end if
    s_min = eval(1)
    s_max = eval(n)
    if (s_min <= 1.0d-14) then
      info = -100
      H_orth(:, :) = (0.0d0, 0.0d0)
      deallocate(svec, xmat, tmp, eval)
      return
    end if
    xmat(:, :) = svec(:, :)
    do i = 1, n
      xmat(:, i) = xmat(:, i) / sqrt(eval(i))
    end do
    xmat(:, :) = matmul(xmat, conjg(transpose(svec)))
    tmp(:, :) = matmul(H_in, xmat)
    H_orth(:, :) = matmul(conjg(transpose(xmat)), tmp)
    call hermitize_matrix(H_orth, n)
    deallocate(svec, xmat, tmp, eval)
  end subroutine build_symmetric_orthogonal_hamiltonian


  subroutine build_hermitian_pseudoinverse(A_in, n, tol, A_inv, info, eval_min, eval_max, nkeep)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: tol
    complex(8), intent(in) :: A_in(n,n)
    complex(8), intent(out) :: A_inv(n,n)
    integer, intent(out) :: info, nkeep
    real(8), intent(out) :: eval_min, eval_max

    integer :: i
    complex(8), allocatable :: avec(:,:)
    real(8), allocatable :: eval(:)
    real(8) :: cutoff

    allocate(avec(n,n), eval(n))
    avec(:, :) = A_in(:, :)
    call zheev_with_query(avec, n, eval, info)
    if (info /= 0) then
      A_inv(:, :) = (0.0d0, 0.0d0)
      eval_min = 0.0d0
      eval_max = 0.0d0
      nkeep = 0
      deallocate(avec, eval)
      return
    end if
    eval_min = eval(1)
    eval_max = eval(n)
    cutoff = max(tol * max(eval_max, 0.0d0), 1.0d-14)
    A_inv(:, :) = (0.0d0, 0.0d0)
    nkeep = 0
    do i = 1, n
      if (eval(i) <= cutoff) cycle
      A_inv(:, :) = A_inv(:, :) + matmul(reshape(avec(:, i), [n,1]), reshape(conjg(avec(:, i)), [1,n])) / eval(i)
      nkeep = nkeep + 1
    end do
    if (nkeep <= 0) then
      info = -100
      A_inv(:, :) = (0.0d0, 0.0d0)
    else
      call hermitize_matrix(A_inv, n)
    end if
    deallocate(avec, eval)
  end subroutine build_hermitian_pseudoinverse


  subroutine build_hermitian_inverse_sqrt(A_in, n, tol, A_invsqrt, info, eval_min, eval_max, nkeep)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: tol
    complex(8), intent(in) :: A_in(n,n)
    complex(8), intent(out) :: A_invsqrt(n,n)
    integer, intent(out) :: info, nkeep
    real(8), intent(out) :: eval_min, eval_max

    integer :: i
    complex(8), allocatable :: avec(:,:)
    real(8), allocatable :: eval(:)
    real(8) :: cutoff

    allocate(avec(n,n), eval(n))
    avec(:, :) = A_in(:, :)
    call zheev_with_query(avec, n, eval, info)
    if (info /= 0) then
      A_invsqrt(:, :) = (0.0d0, 0.0d0)
      eval_min = 0.0d0
      eval_max = 0.0d0
      nkeep = 0
      deallocate(avec, eval)
      return
    end if
    eval_min = eval(1)
    eval_max = eval(n)
    cutoff = max(tol * max(eval_max, 0.0d0), 1.0d-14)
    A_invsqrt(:, :) = (0.0d0, 0.0d0)
    nkeep = 0
    do i = 1, n
      if (eval(i) <= cutoff) cycle
      A_invsqrt(:, :) = A_invsqrt(:, :) + &
        matmul(reshape(avec(:, i), [n,1]), reshape(conjg(avec(:, i)), [1,n])) / sqrt(eval(i))
      nkeep = nkeep + 1
    end do
    if (nkeep <= 0) then
      info = -100
      A_invsqrt(:, :) = (0.0d0, 0.0d0)
    else
      call hermitize_matrix(A_invsqrt, n)
    end if
    deallocate(avec, eval)
  end subroutine build_hermitian_inverse_sqrt


  subroutine build_hermitian_inverse(A_in, n, A_inv, info, eval_min, eval_max)
    implicit none
    integer, intent(in) :: n
    complex(8), intent(in) :: A_in(n,n)
    complex(8), intent(out) :: A_inv(n,n)
    integer, intent(out) :: info
    real(8), intent(out) :: eval_min, eval_max

    integer :: i
    complex(8), allocatable :: avec(:,:)
    real(8), allocatable :: eval(:)

    allocate(avec(n,n), eval(n))
    avec(:, :) = A_in(:, :)
    call zheev_with_query(avec, n, eval, info)
    if (info /= 0) then
      A_inv(:, :) = (0.0d0, 0.0d0)
      eval_min = 0.0d0
      eval_max = 0.0d0
      deallocate(avec, eval)
      return
    end if
    eval_min = eval(1)
    eval_max = eval(n)
    if (eval_min <= 1.0d-14) then
      info = -100
      A_inv(:, :) = (0.0d0, 0.0d0)
      deallocate(avec, eval)
      return
    end if
    A_inv(:, :) = (0.0d0, 0.0d0)
    do i = 1, n
      A_inv(:, :) = A_inv(:, :) + matmul(reshape(avec(:, i), [n,1]), reshape(conjg(avec(:, i)), [1,n])) / eval(i)
    end do
    call hermitize_matrix(A_inv, n)
    deallocate(avec, eval)
  end subroutine build_hermitian_inverse


  subroutine build_wpw_reduced_c_can_reference(S_build, S_red_inv, reduced_transform, c_source, &
      nraw, nred, c_can, c_red, info)
    implicit none
    integer, intent(in) :: nraw, nred
    complex(8), intent(in) :: S_build(:,:), S_red_inv(:,:), reduced_transform(:,:), c_source(:)
    complex(8), intent(out) :: c_can(:), c_red(:)
    integer, intent(out) :: info
    complex(8), allocatable :: tmp_raw(:), rhs_red(:)

    info = 0
    if (nraw <= 0 .or. nred <= 0) then
      info = 1
      return
    end if
    if (size(S_build,1) < nraw .or. size(S_build,2) < nraw .or. &
        size(S_red_inv,1) < nred .or. size(S_red_inv,2) < nred .or. &
        size(reduced_transform,1) < nraw .or. size(reduced_transform,2) < nred .or. &
        size(c_source) < nraw .or. size(c_can) < nraw .or. size(c_red) < nred) then
      info = 1
      return
    end if

    allocate(tmp_raw(nraw), rhs_red(nred))
    tmp_raw(:) = matmul(S_build(1:nraw,1:nraw), c_source(1:nraw))
    rhs_red(:) = matmul(conjg(transpose(reduced_transform(1:nraw,1:nred))), tmp_raw(:))
    c_red(1:nred) = matmul(S_red_inv(1:nred,1:nred), rhs_red(:))
    c_can(1:nraw) = matmul(reduced_transform(1:nraw,1:nred), c_red(1:nred))
    deallocate(tmp_raw, rhs_red)
  end subroutine build_wpw_reduced_c_can_reference


  subroutine build_wpw_reduced_raw_back_hybrid_from_inverse(S_source_inv, rhs_source, S_hybrid, &
      S_red_inv, reduced_transform, nraw, nred, c_can, c_red, c_source, info)
    implicit none
    integer, intent(in) :: nraw, nred
    complex(8), intent(in) :: S_source_inv(:,:), rhs_source(:), S_hybrid(:,:)
    complex(8), intent(in) :: S_red_inv(:,:), reduced_transform(:,:)
    complex(8), intent(out) :: c_can(:), c_red(:), c_source(:)
    integer, intent(out) :: info

    info = 0
    if (nraw <= 0 .or. nred <= 0) then
      info = 1
      return
    end if
    if (size(S_source_inv,1) < nraw .or. size(S_source_inv,2) < nraw .or. &
        size(rhs_source) < nraw .or. size(c_source) < nraw) then
      info = 1
      return
    end if

    c_source(1:nraw) = matmul(S_source_inv(1:nraw,1:nraw), rhs_source(1:nraw))
    call build_wpw_reduced_c_can_reference(S_hybrid, S_red_inv, reduced_transform, &
      c_source, nraw, nred, c_can, c_red, info)
  end subroutine build_wpw_reduced_raw_back_hybrid_from_inverse


  subroutine build_wpw_reduced_raw_back_hybrid_matrix(S_source, rhs_source, S_hybrid, S_red_inv, &
      reduced_transform, nraw, nred, nstate, c_can, c_source, info)
    implicit none
    integer, intent(in) :: nraw, nred, nstate
    complex(8), intent(in) :: S_source(:,:), rhs_source(:,:), S_hybrid(:,:)
    complex(8), intent(in) :: S_red_inv(:,:), reduced_transform(:,:)
    complex(8), intent(out) :: c_can(:,:), c_source(:,:)
    integer, intent(out) :: info

    integer :: ist, info_local, nkeep
    real(8) :: eval_min, eval_max
    complex(8), allocatable :: S_source_inv(:,:), c_red(:)

    info = 0
    if (nraw <= 0 .or. nred <= 0 .or. nstate <= 0) then
      info = 1
      return
    end if
    if (size(S_source,1) < nraw .or. size(S_source,2) < nraw .or. &
        size(rhs_source,1) < nraw .or. size(rhs_source,2) < nstate .or. &
        size(c_can,1) < nraw .or. size(c_can,2) < nstate .or. &
        size(c_source,1) < nraw .or. size(c_source,2) < nstate) then
      info = 1
      return
    end if

    allocate(S_source_inv(nraw,nraw), c_red(nred))
    call build_hermitian_pseudoinverse(S_source(1:nraw,1:nraw), nraw, 1.0d-8, &
      S_source_inv, info_local, eval_min, eval_max, nkeep)
    if (info_local /= 0) then
      c_can(1:nraw,1:nstate) = (0.0d0, 0.0d0)
      c_source(1:nraw,1:nstate) = (0.0d0, 0.0d0)
      info = 2
      deallocate(S_source_inv, c_red)
      return
    end if

    do ist = 1, nstate
      call build_wpw_reduced_raw_back_hybrid_from_inverse(S_source_inv, rhs_source(1:nraw,ist), &
        S_hybrid, S_red_inv, reduced_transform, nraw, nred, c_can(1:nraw,ist), c_red, &
        c_source(1:nraw,ist), info_local)
      if (info_local /= 0) info = 3
    end do
    deallocate(S_source_inv, c_red)
  end subroutine build_wpw_reduced_raw_back_hybrid_matrix


  subroutine zheev_with_query(A, n, eval, info)
    implicit none
    integer, intent(in) :: n
    complex(8), intent(inout) :: A(n,n)
    real(8), intent(out) :: eval(n)
    integer, intent(out) :: info

    integer :: lwork
    complex(8), allocatable :: work(:)
    real(8), allocatable :: rwork(:)
    external :: zheev

    allocate(work(1), rwork(max(1, 3*n - 2)))
    lwork = -1
    call ZHEEV('V', 'U', n, A, n, eval, work, lwork, rwork, info)
    lwork = max(1, int(real(work(1), kind=8)))
    deallocate(work)
    allocate(work(lwork))
    call ZHEEV('V', 'U', n, A, n, eval, work, lwork, rwork, info)
    deallocate(work, rwork)
  end subroutine zheev_with_query


  subroutine zhegv_with_query(H, S, n, eval, info)
    implicit none
    integer, intent(in) :: n
    complex(8), intent(inout) :: H(n,n), S(n,n)
    real(8), intent(out) :: eval(n)
    integer, intent(out) :: info

    integer :: lwork
    complex(8), allocatable :: work(:)
    real(8), allocatable :: rwork(:)
    external :: zhegv

    allocate(work(1), rwork(max(1, 3*n - 2)))
    lwork = -1
    call ZHEGV(1, 'V', 'U', n, H, n, S, n, eval, work, lwork, rwork, info)
    lwork = max(1, int(real(work(1), kind=8)))
    deallocate(work)
    allocate(work(lwork))
    call ZHEGV(1, 'V', 'U', n, H, n, S, n, eval, work, lwork, rwork, info)
    deallocate(work, rwork)
  end subroutine zhegv_with_query


  subroutine rectangular_singular_minmax(mat, nrow, ncol, sigma_min, sigma_max)
    implicit none
    integer, intent(in) :: nrow, ncol
    complex(8), intent(in) :: mat(nrow, ncol)
    real(8), intent(out) :: sigma_min, sigma_max

    integer :: ngram, info, lwork
    complex(8), allocatable :: gram(:,:), work(:)
    real(8), allocatable :: eval(:), rwork(:)
    external :: zheev

    sigma_min = 0.0d0
    sigma_max = 0.0d0
    if (nrow <= 0 .or. ncol <= 0) return
    ngram = min(nrow, ncol)
    allocate(gram(ngram, ngram), eval(ngram), rwork(max(1, 3*ngram - 2)), work(1))
    if (ncol <= nrow) then
      gram(:, :) = matmul(conjg(transpose(mat(1:nrow, 1:ncol))), mat(1:nrow, 1:ncol))
    else
      gram(:, :) = matmul(mat(1:nrow, 1:ncol), conjg(transpose(mat(1:nrow, 1:ncol))))
    end if
    call hermitize_matrix(gram, ngram)
    lwork = -1
    call ZHEEV('N', 'U', ngram, gram, ngram, eval, work, lwork, rwork, info)
    lwork = max(1, int(real(work(1), kind=8)))
    deallocate(work)
    allocate(work(lwork))
    if (ncol <= nrow) then
      gram(:, :) = matmul(conjg(transpose(mat(1:nrow, 1:ncol))), mat(1:nrow, 1:ncol))
    else
      gram(:, :) = matmul(mat(1:nrow, 1:ncol), conjg(transpose(mat(1:nrow, 1:ncol))))
    end if
    call hermitize_matrix(gram, ngram)
    call ZHEEV('N', 'U', ngram, gram, ngram, eval, work, lwork, rwork, info)
    if (info == 0) then
      sigma_min = sqrt(max(0.0d0, eval(1)))
      sigma_max = sqrt(max(0.0d0, eval(ngram)))
    end if
    deallocate(gram, eval, rwork, work)
  end subroutine rectangular_singular_minmax


  subroutine hermitize_matrix(mat, n)
    implicit none
    integer, intent(in) :: n
    complex(8), intent(inout) :: mat(n,n)
    integer :: i, j
    complex(8) :: avg

    do i = 1, n
      mat(i,i) = cmplx(real(mat(i,i), kind=8), 0.0d0, kind=8)
      do j = i + 1, n
        avg = 0.5d0 * (mat(i,j) + conjg(mat(j,i)))
        mat(i,j) = avg
        mat(j,i) = conjg(avg)
      end do
    end do
  end subroutine hermitize_matrix

  subroutine wpw_local_herm_max(mat, n, herm_max)
    implicit none
    integer, intent(in) :: n
    complex(8), intent(in) :: mat(n, n)
    real(8), intent(out) :: herm_max
    integer :: i, j

    herm_max = 0.0d0
    do j = 1, n
      do i = 1, n
        herm_max = max(herm_max, abs(mat(i, j) - conjg(mat(j, i))))
      end do
    end do
  end subroutine wpw_local_herm_max

end module rt_dg_wpw_linalg
