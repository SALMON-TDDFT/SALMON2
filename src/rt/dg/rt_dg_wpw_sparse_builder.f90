module rt_dg_wpw_sparse_builder
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use rt_dg_wpw_column_layout, only: s_dg_wpw_column_layout, wpw_column_pair, wpw_column_owner
  use rt_dg_wpw_sparse_blocks, only: s_dg_wpw_sparse_blocks
  implicit none
  private

  integer, parameter, public :: wpw_candidate_volume = 1
  integer, parameter, public :: wpw_candidate_face = 2

  type, public :: s_dg_wpw_sparse_candidates
    integer, allocatable :: wp_w_row_ids(:), wp_pw_col_ids(:)
    integer, allocatable :: wp_origin(:)
    complex(8), allocatable :: wp_h_values(:,:), wp_s_values(:,:)
    integer, allocatable :: pp_pw_row_ids(:), pp_pw_col_ids(:)
    integer, allocatable :: pp_origin(:)
    complex(8), allocatable :: pp_h_values(:,:), pp_s_values(:,:)
  end type s_dg_wpw_sparse_candidates

  public :: build_windowed_sparse_wpw_operators
  public :: wpw_normalized_window_at_grid
  public :: wpw_grad_chi

contains

  subroutine build_windowed_sparse_wpw_operators(layout, rank_id, nrank, candidates, h_blocks, s_blocks, info)
    type(s_dg_wpw_column_layout), intent(in) :: layout
    integer, intent(in) :: rank_id, nrank
    type(s_dg_wpw_sparse_candidates), intent(in) :: candidates
    type(s_dg_wpw_sparse_blocks), intent(out) :: h_blocks, s_blocks
    integer, intent(out) :: info

    integer :: i, nspin, nwp, npp, keep_wp, keep_pp, owner, owner_info
    integer :: fragment_id, g_id, pair_info

    info = 0
    if (trim(layout%basis_kind) /= 'windowed_kg' .or. layout%n_global_columns <= 0 .or. &
        layout%n_g_modes <= 0) then
      info = 1
      return
    end if
    if (rank_id < 0 .or. rank_id >= nrank .or. nrank <= 0) then
      info = 2
      return
    end if
    if (.not. allocated(layout%owned_column_ids) .or. .not. allocated(layout%pw_fragment_ids) .or. &
        .not. allocated(layout%pw_g_ids) .or. size(layout%pw_fragment_ids) /= size(layout%owned_column_ids) .or. &
        size(layout%pw_g_ids) /= size(layout%owned_column_ids)) then
      info = 10
      return
    end if
    do i = 1, size(layout%owned_column_ids)
      call wpw_column_pair(layout%owned_column_ids(i), layout%n_g_modes, fragment_id, g_id, pair_info)
      if (pair_info /= 0 .or. fragment_id /= layout%pw_fragment_ids(i) .or. g_id /= layout%pw_g_ids(i)) then
        info = 11
        return
      end if
    end do
    if (.not. candidates_ready(candidates, nspin, nwp, npp)) then
      info = 3
      return
    end if

    keep_wp = 0
    do i = 1, nwp
      call validate_column(candidates%wp_pw_col_ids(i), owner, fragment_id, g_id, owner_info)
      if (owner_info /= 0 .or. candidates%wp_w_row_ids(i) <= 0) then
        info = 4
        return
      end if
      if (.not. finite_row(candidates%wp_h_values(i,:)) .or. &
          .not. finite_row(candidates%wp_s_values(i,:))) then
        info = 5
        return
      end if
      if (owner /= rank_id) then
        info = 12
        return
      end if
      if (candidates%wp_origin(i) /= wpw_candidate_volume .and. &
          candidates%wp_origin(i) /= wpw_candidate_face) then
        info = 14
        return
      end if
      keep_wp = keep_wp + 1
    end do
    keep_pp = 0
    do i = 1, npp
      call validate_column(candidates%pp_pw_row_ids(i), owner, fragment_id, g_id, owner_info)
      if (owner_info /= 0) then
        info = 6
        return
      end if
      call wpw_column_pair(candidates%pp_pw_col_ids(i), layout%n_g_modes, fragment_id, g_id, pair_info)
      if (pair_info /= 0 .or. candidates%pp_pw_col_ids(i) > layout%n_global_columns) then
        info = 7
        return
      end if
      if (.not. finite_row(candidates%pp_h_values(i,:)) .or. &
          .not. finite_row(candidates%pp_s_values(i,:))) then
        info = 8
        return
      end if
      if (owner /= rank_id) then
        info = 12
        return
      end if
      if (candidates%pp_origin(i) /= wpw_candidate_volume) then
        info = 13
        return
      end if
      keep_pp = keep_pp + 1
    end do
    if (.not. pairs_strictly_increasing(candidates%wp_pw_col_ids, candidates%wp_w_row_ids) .or. &
        .not. pairs_strictly_increasing(candidates%pp_pw_row_ids, candidates%pp_pw_col_ids)) then
      info = 9
      return
    end if

    call allocate_blocks(h_blocks, keep_wp, keep_pp, nspin)
    call allocate_blocks(s_blocks, keep_wp, keep_pp, nspin)
    keep_wp = 0
    do i = 1, nwp
      keep_wp = keep_wp + 1
      h_blocks%wp_w_row_ids(keep_wp) = candidates%wp_w_row_ids(i)
      h_blocks%wp_pw_col_ids(keep_wp) = candidates%wp_pw_col_ids(i)
      h_blocks%wp_values(keep_wp,:) = candidates%wp_h_values(i,:)
      s_blocks%wp_w_row_ids(keep_wp) = candidates%wp_w_row_ids(i)
      s_blocks%wp_pw_col_ids(keep_wp) = candidates%wp_pw_col_ids(i)
      s_blocks%wp_values(keep_wp,:) = candidates%wp_s_values(i,:)
    end do
    keep_pp = 0
    do i = 1, npp
      keep_pp = keep_pp + 1
      h_blocks%pp_pw_row_ids(keep_pp) = candidates%pp_pw_row_ids(i)
      h_blocks%pp_pw_col_ids(keep_pp) = candidates%pp_pw_col_ids(i)
      h_blocks%pp_values(keep_pp,:) = candidates%pp_h_values(i,:)
      s_blocks%pp_pw_row_ids(keep_pp) = candidates%pp_pw_row_ids(i)
      s_blocks%pp_pw_col_ids(keep_pp) = candidates%pp_pw_col_ids(i)
      s_blocks%pp_values(keep_pp,:) = candidates%pp_s_values(i,:)
    end do
    info = 0

  contains
    subroutine validate_column(column_id, column_owner, k_id, local_g_id, validation_info)
      integer, intent(in) :: column_id
      integer, intent(out) :: column_owner, k_id, local_g_id, validation_info
      call wpw_column_pair(column_id, layout%n_g_modes, k_id, local_g_id, validation_info)
      if (validation_info /= 0 .or. column_id > layout%n_global_columns) return
      column_owner = wpw_column_owner(column_id, layout%n_global_columns, nrank, validation_info)
    end subroutine validate_column
  end subroutine build_windowed_sparse_wpw_operators

  complex(8) function wpw_normalized_window_at_grid(chi, plane_phase, omega_cell, info) result(value)
    real(8), intent(in) :: chi, omega_cell
    complex(8), intent(in) :: plane_phase
    integer, intent(out) :: info
    value = (0.0d0, 0.0d0)
    info = 0
    if (omega_cell <= 0.0d0 .or. .not. ieee_is_finite(chi) .or. &
        .not. ieee_is_finite(real(plane_phase,8)) .or. .not. ieee_is_finite(aimag(plane_phase))) then
      info = 1
      return
    end if
    value = chi * plane_phase / sqrt(omega_cell)
  end function wpw_normalized_window_at_grid

  subroutine wpw_grad_chi(chi, grad_chi, gvec, plane_phase, omega_cell, grad_p, info)
    real(8), intent(in) :: chi, grad_chi(3), gvec(3), omega_cell
    complex(8), intent(in) :: plane_phase
    complex(8), intent(out) :: grad_p(3)
    integer, intent(out) :: info
    integer :: idir
    grad_p = (0.0d0, 0.0d0)
    info = 0
    if (omega_cell <= 0.0d0 .or. .not. ieee_is_finite(chi) .or. &
        any(.not. ieee_is_finite(grad_chi)) .or. any(.not. ieee_is_finite(gvec)) .or. &
        .not. ieee_is_finite(real(plane_phase,8)) .or. .not. ieee_is_finite(aimag(plane_phase))) then
      info = 1
      return
    end if
    do idir = 1, 3
      grad_p(idir) = (grad_chi(idir) + cmplx(0.0d0, gvec(idir) * chi, kind=8)) * &
                     plane_phase / sqrt(omega_cell)
    end do
  end subroutine wpw_grad_chi

  logical function candidates_ready(candidates, nspin, nwp, npp) result(ok)
    type(s_dg_wpw_sparse_candidates), intent(in) :: candidates
    integer, intent(out) :: nspin, nwp, npp
    ok = .false.; nspin = 0; nwp = 0; npp = 0
    if (.not. allocated(candidates%wp_w_row_ids) .or. .not. allocated(candidates%wp_pw_col_ids) .or. &
        .not. allocated(candidates%wp_origin) .or. &
        .not. allocated(candidates%wp_h_values) .or. .not. allocated(candidates%wp_s_values) .or. &
        .not. allocated(candidates%pp_pw_row_ids) .or. .not. allocated(candidates%pp_pw_col_ids) .or. &
        .not. allocated(candidates%pp_origin) .or. &
        .not. allocated(candidates%pp_h_values) .or. .not. allocated(candidates%pp_s_values)) return
    nwp = size(candidates%wp_w_row_ids); npp = size(candidates%pp_pw_row_ids)
    nspin = size(candidates%wp_h_values,2)
    if (nspin <= 0) return
    ok = size(candidates%wp_pw_col_ids) == nwp .and. size(candidates%wp_origin) == nwp .and. &
         size(candidates%wp_h_values,1) == nwp .and. &
         all(shape(candidates%wp_s_values) == [nwp,nspin]) .and. &
         size(candidates%pp_pw_col_ids) == npp .and. size(candidates%pp_origin) == npp .and. &
         size(candidates%pp_h_values,1) == npp .and. &
         all(shape(candidates%pp_s_values) == [npp,nspin]) .and. size(candidates%pp_h_values,2) == nspin
  end function candidates_ready

  logical function finite_row(values) result(ok)
    complex(8), intent(in) :: values(:)
    integer :: i
    ok = .true.
    do i = 1, size(values)
      if (.not. ieee_is_finite(real(values(i),8)) .or. .not. ieee_is_finite(aimag(values(i)))) ok = .false.
    end do
  end function finite_row

  logical function pairs_strictly_increasing(primary_ids, secondary_ids) result(ok)
    integer, intent(in) :: primary_ids(:), secondary_ids(:)
    integer :: i
    ok = size(primary_ids) == size(secondary_ids)
    if (.not. ok) return
    do i = 2, size(primary_ids)
      if (primary_ids(i) < primary_ids(i-1) .or. &
          (primary_ids(i) == primary_ids(i-1) .and. secondary_ids(i) <= secondary_ids(i-1))) then
        ok = .false.
        return
      end if
    end do
  end function pairs_strictly_increasing

  subroutine allocate_blocks(blocks, nwp, npp, nspin)
    type(s_dg_wpw_sparse_blocks), intent(out) :: blocks
    integer, intent(in) :: nwp, npp, nspin
    allocate(blocks%wp_w_row_ids(nwp), blocks%wp_pw_col_ids(nwp), blocks%wp_values(nwp,nspin))
    allocate(blocks%pp_pw_row_ids(npp), blocks%pp_pw_col_ids(npp), blocks%pp_values(npp,nspin))
  end subroutine allocate_blocks
end module rt_dg_wpw_sparse_builder
