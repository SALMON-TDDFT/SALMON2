module rt_dg_wpw_sparse_blocks
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  type, public :: s_dg_wpw_sparse_blocks
    integer, allocatable :: wp_w_row_ids(:)
    integer, allocatable :: wp_pw_col_ids(:)
    complex(8), allocatable :: wp_values(:,:) ! (entry, spin)
    integer, allocatable :: pp_pw_row_ids(:)
    integer, allocatable :: pp_pw_col_ids(:)
    complex(8), allocatable :: pp_values(:,:) ! (entry, spin)
  end type s_dg_wpw_sparse_blocks

  public :: apply_wp_owned_columns
  public :: apply_pp_owned_rows

contains

  subroutine apply_wp_owned_columns(blocks, ispin, w_row_ids, pw_row_ids, xw, xp, &
                                    yw_partial, yp_owned, info)
    type(s_dg_wpw_sparse_blocks), intent(in) :: blocks
    integer, intent(in) :: ispin
    integer, intent(in) :: w_row_ids(:), pw_row_ids(:)
    complex(8), intent(in) :: xw(:,:), xp(:,:)
    complex(8), intent(inout) :: yw_partial(:,:), yp_owned(:,:)
    integer, intent(out) :: info

    integer :: ientry, iw, ip, ivec, nentry, nvec
    complex(8) :: value

    info = 0
    if (.not. allocated(blocks%wp_w_row_ids) .or. &
        .not. allocated(blocks%wp_pw_col_ids) .or. .not. allocated(blocks%wp_values)) then
      info = 1
      return
    end if
    nentry = size(blocks%wp_w_row_ids)
    if (size(blocks%wp_pw_col_ids) /= nentry .or. size(blocks%wp_values, 1) /= nentry) then
      info = 2
      return
    end if
    if (ispin < 1 .or. ispin > size(blocks%wp_values, 2)) then
      info = 3
      return
    end if
    nvec = size(xw, 2)
    if (size(xp, 2) /= nvec .or. size(yw_partial, 2) /= nvec .or. size(yp_owned, 2) /= nvec .or. &
        size(xw, 1) /= size(w_row_ids) .or. size(yw_partial, 1) /= size(w_row_ids) .or. &
        size(xp, 1) /= size(pw_row_ids) .or. size(yp_owned, 1) /= size(pw_row_ids)) then
      info = 4
      return
    end if
    if (.not. strictly_increasing(w_row_ids) .or. .not. strictly_increasing(pw_row_ids)) then
      info = 5
      return
    end if
    do ientry = 1, nentry
      iw = find_sorted_id(w_row_ids, blocks%wp_w_row_ids(ientry))
      ip = find_sorted_id(pw_row_ids, blocks%wp_pw_col_ids(ientry))
      if (iw <= 0 .or. ip <= 0) then
        info = 6
        return
      end if
      value = blocks%wp_values(ientry, ispin)
      if (.not. ieee_is_finite(real(value, 8)) .or. .not. ieee_is_finite(aimag(value))) then
        info = 7
        return
      end if
    end do

    do ientry = 1, nentry
      iw = find_sorted_id(w_row_ids, blocks%wp_w_row_ids(ientry))
      ip = find_sorted_id(pw_row_ids, blocks%wp_pw_col_ids(ientry))
      value = blocks%wp_values(ientry, ispin)
      do ivec = 1, nvec
        yw_partial(iw, ivec) = yw_partial(iw, ivec) + value * xp(ip, ivec)
        yp_owned(ip, ivec) = yp_owned(ip, ivec) + conjg(value) * xw(iw, ivec)
      end do
    end do
  end subroutine apply_wp_owned_columns

  subroutine apply_pp_owned_rows(blocks, ispin, pw_input_ids, pw_owned_ids, xp, yp_owned, info)
    type(s_dg_wpw_sparse_blocks), intent(in) :: blocks
    integer, intent(in) :: ispin
    integer, intent(in) :: pw_input_ids(:), pw_owned_ids(:)
    complex(8), intent(in) :: xp(:,:)
    complex(8), intent(inout) :: yp_owned(:,:)
    integer, intent(out) :: info

    integer :: ientry, irow, icol, ivec, nentry, nvec
    complex(8) :: value

    info = 0
    if (.not. allocated(blocks%pp_pw_row_ids) .or. &
        .not. allocated(blocks%pp_pw_col_ids) .or. .not. allocated(blocks%pp_values)) then
      info = 1
      return
    end if
    nentry = size(blocks%pp_pw_row_ids)
    if (size(blocks%pp_pw_col_ids) /= nentry .or. size(blocks%pp_values, 1) /= nentry) then
      info = 2
      return
    end if
    if (ispin < 1 .or. ispin > size(blocks%pp_values, 2)) then
      info = 3
      return
    end if
    nvec = size(xp, 2)
    if (size(xp, 1) /= size(pw_input_ids) .or. size(yp_owned, 1) /= size(pw_owned_ids) .or. &
        size(yp_owned, 2) /= nvec) then
      info = 4
      return
    end if
    if (.not. strictly_increasing(pw_input_ids) .or. .not. strictly_increasing(pw_owned_ids)) then
      info = 5
      return
    end if
    do ientry = 1, nentry
      irow = find_sorted_id(pw_owned_ids, blocks%pp_pw_row_ids(ientry))
      icol = find_sorted_id(pw_input_ids, blocks%pp_pw_col_ids(ientry))
      if (irow <= 0 .or. icol <= 0) then
        info = 6
        return
      end if
      value = blocks%pp_values(ientry, ispin)
      if (.not. ieee_is_finite(real(value, 8)) .or. .not. ieee_is_finite(aimag(value))) then
        info = 7
        return
      end if
    end do

    do ientry = 1, nentry
      irow = find_sorted_id(pw_owned_ids, blocks%pp_pw_row_ids(ientry))
      icol = find_sorted_id(pw_input_ids, blocks%pp_pw_col_ids(ientry))
      value = blocks%pp_values(ientry, ispin)
      do ivec = 1, nvec
        yp_owned(irow, ivec) = yp_owned(irow, ivec) + value * xp(icol, ivec)
      end do
    end do
  end subroutine apply_pp_owned_rows

  logical function strictly_increasing(ids) result(ok)
    integer, intent(in) :: ids(:)
    integer :: i

    ok = .true.
    do i = 2, size(ids)
      if (ids(i) <= ids(i-1)) then
        ok = .false.
        return
      end if
    end do
  end function strictly_increasing

  integer function find_sorted_id(ids, target) result(position)
    integer, intent(in) :: ids(:), target
    integer :: left, middle, right

    position = 0
    left = 1
    right = size(ids)
    do while (left <= right)
      middle = left + (right - left) / 2
      if (ids(middle) == target) then
        position = middle
        return
      else if (ids(middle) < target) then
        left = middle + 1
      else
        right = middle - 1
      end if
    end do
  end function find_sorted_id

end module rt_dg_wpw_sparse_blocks
