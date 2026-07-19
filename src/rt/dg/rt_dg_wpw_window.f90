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
! WPW window and periodic-grid geometry helpers
!=======================================================================

module rt_dg_wpw_window
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use dg_wpw_windows, only: evaluate_dg_wpw_normalized_windows
  use salmon_global, only: dg_wpw_window_buffer,dg_wpw_window_width
  implicit none

  private
  public :: wpw_window_buffer_axis
  public :: wpw_window_transition_width_axis
  public :: wpw_fragment_box_size
  public :: wpw_fragment_box_bounds
  public :: wpw_raw_window_at_grid
  public :: wpw_normalized_window_at_grid
  public :: wpw_face_neighbor_fragment
  public :: wpw_local_is_neighbor_pair
  public :: map_global_to_phi_box_coord_pw

contains

  integer function wpw_window_buffer_axis(dg_frag, axis) result(buf)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: axis
    buf = max(0, dg_wpw_window_buffer)
  end function wpw_window_buffer_axis


  integer function wpw_window_transition_width_axis(dg_frag, axis) result(width)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: axis
    width = dg_wpw_window_width
    width = max(0, min(width, max(0, wpw_window_buffer_axis(dg_frag, axis))))
  end function wpw_window_transition_width_axis


  subroutine wpw_fragment_box_size(dg_frag, ifrag, nx, ny, nz)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag
    integer, intent(out) :: nx, ny, nz
    integer :: lo(3), hi(3)

    call wpw_fragment_box_bounds(dg_frag, ifrag, lo, hi)
    nx = hi(1) - lo(1) + 1
    ny = hi(2) - lo(2) + 1
    nz = hi(3) - lo(3) + 1
  end subroutine wpw_fragment_box_size


  subroutine wpw_fragment_box_bounds(dg_frag, ifrag, lo, hi)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag
    integer, intent(out) :: lo(3), hi(3)

    if (allocated(dg_frag%frag_buf_lo) .and. allocated(dg_frag%frag_buf_hi)) then
      lo(1:3) = dg_frag%frag_buf_lo(1:3, ifrag)
      hi(1:3) = dg_frag%frag_buf_hi(1:3, ifrag)
    else
      lo(1) = dg_frag%ixyz_frag(1, ifrag) - wpw_window_buffer_axis(dg_frag, 1)
      lo(2) = dg_frag%ixyz_frag(2, ifrag) - wpw_window_buffer_axis(dg_frag, 2)
      lo(3) = dg_frag%ixyz_frag(3, ifrag) - wpw_window_buffer_axis(dg_frag, 3)
      hi(1) = dg_frag%ixyz_frag(1, ifrag) + dg_frag%nxyz_domain(1, ifrag) - 1 + wpw_window_buffer_axis(dg_frag, 1)
      hi(2) = dg_frag%ixyz_frag(2, ifrag) + dg_frag%nxyz_domain(2, ifrag) - 1 + wpw_window_buffer_axis(dg_frag, 2)
      hi(3) = dg_frag%ixyz_frag(3, ifrag) + dg_frag%nxyz_domain(3, ifrag) - 1 + wpw_window_buffer_axis(dg_frag, 3)
    end if
  end subroutine wpw_fragment_box_bounds


  subroutine wpw_raw_window_at_grid(dg_frag, ifrag, gx, gy, gz, q, grad_q)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, gx, gy, gz
    real(8), intent(out) :: q, grad_q(3)
    real(8) :: wx(3), dwx(3)

    call wpw_axis_window(dg_frag, 1, ifrag, gx, wx(1), dwx(1))
    call wpw_axis_window(dg_frag, 2, ifrag, gy, wx(2), dwx(2))
    call wpw_axis_window(dg_frag, 3, ifrag, gz, wx(3), dwx(3))

    q = wx(1) * wx(2) * wx(3)
    grad_q(1) = dwx(1) * wx(2) * wx(3)
    grad_q(2) = wx(1) * dwx(2) * wx(3)
    grad_q(3) = wx(1) * wx(2) * dwx(3)
  end subroutine wpw_raw_window_at_grid


  subroutine wpw_axis_window(dg_frag, axis, ifrag, g, w, dw)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: axis, ifrag, g
    real(8), intent(out) :: w, dw
    integer :: lo0, hi0, lo, hi, buf, width, ltot, shift0, shift
    real(8) :: t, arg, h
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0)

    w = 0.0d0
    dw = 0.0d0
    if (allocated(dg_frag%frag_core_lo) .and. allocated(dg_frag%frag_core_hi)) then
      lo0 = dg_frag%frag_core_lo(axis, ifrag)
      hi0 = dg_frag%frag_core_hi(axis, ifrag)
    else
      lo0 = dg_frag%ixyz_frag(axis, ifrag)
      hi0 = dg_frag%ixyz_frag(axis, ifrag) + dg_frag%nxyz_domain(axis, ifrag) - 1
    end if
    buf = wpw_window_buffer_axis(dg_frag, axis)
    width = wpw_window_transition_width_axis(dg_frag, axis)
    ltot = max(1, dg_frag%lgnum_total(axis))
    h = max(dg_frag%hgs(axis), 1.0d-30)
    shift0 = nint(dble(g - lo0) / dble(ltot))

    do shift = shift0 - 1, shift0 + 1
      lo = lo0 + shift * ltot
      hi = hi0 + shift * ltot
      if (g >= lo .and. g <= hi) then
        w = 1.0d0
        dw = 0.0d0
        return
      end if
      if (buf <= 0 .or. width <= 0) cycle
      if (g >= lo - width .and. g < lo) then
        t = dble(g - (lo - width)) / dble(width)
        t = max(0.0d0, min(1.0d0, t))
        arg = 0.5d0 * pi * t
        w = sin(arg)
        dw = 0.5d0 * pi * cos(arg) / (dble(width) * h)
        return
      end if
      if (g > hi .and. g <= hi + width) then
        t = dble(g - hi) / dble(width)
        t = max(0.0d0, min(1.0d0, t))
        arg = 0.5d0 * pi * t
        w = cos(arg)
        dw = -0.5d0 * pi * sin(arg) / (dble(width) * h)
        return
      end if
    end do
  end subroutine wpw_axis_window


  subroutine wpw_normalized_window_at_grid(dg_frag, ifrag, gx, gy, gz, chi, grad_chi)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, gx, gy, gz
    real(8), intent(out) :: chi, grad_chi(3)
    integer :: jfrag, nfrag, info
    integer :: total_grid(3), buffer(3), width(3)
    integer, allocatable :: core_lo(:,:), core_hi(:,:)
    real(8), allocatable :: chi_all(:), grad_chi_all(:,:)

    nfrag = dg_frag%n_frag
    chi = 0.0d0
    grad_chi(:) = 0.0d0
    if (ifrag < 1 .or. ifrag > nfrag) return
    allocate(core_lo(3, nfrag), core_hi(3, nfrag), chi_all(nfrag), grad_chi_all(3, nfrag))
    if (allocated(dg_frag%frag_core_lo) .and. allocated(dg_frag%frag_core_hi)) then
      core_lo = dg_frag%frag_core_lo(:, 1:nfrag)
      core_hi = dg_frag%frag_core_hi(:, 1:nfrag)
    else
      do jfrag = 1, nfrag
        core_lo(:, jfrag) = dg_frag%ixyz_frag(:, jfrag)
        core_hi(:, jfrag) = dg_frag%ixyz_frag(:, jfrag) + dg_frag%nxyz_domain(:, jfrag) - 1
      end do
    end if
    total_grid = max(1, dg_frag%lgnum_total)
    do jfrag = 1, 3
      buffer(jfrag) = wpw_window_buffer_axis(dg_frag, jfrag)
      width(jfrag) = wpw_window_transition_width_axis(dg_frag, jfrag)
    end do
    call evaluate_dg_wpw_normalized_windows(core_lo, core_hi, total_grid, dg_frag%hgs, buffer, width, &
      [gx, gy, gz], chi_all, grad_chi_all, info)
    if (info == 0) then
      chi = chi_all(ifrag)
      grad_chi = grad_chi_all(:, ifrag)
    end if
    deallocate(core_lo, core_hi, chi_all, grad_chi_all)
  end subroutine wpw_normalized_window_at_grid


  integer function wpw_face_neighbor_fragment(dg_frag, ifrag, axis, side) result(jfrag)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, axis, side
    integer :: cand, lo(3), hi(3), target, cand_face

    jfrag = 0
    target = modulo(dg_frag%ixyz_frag(axis, ifrag) + merge(dg_frag%nxyz_domain(axis, ifrag), -1, side > 0) - 1, &
      max(1, dg_frag%lgnum_total(axis))) + 1
    do cand = 1, dg_frag%n_frag
      if (cand == ifrag) cycle
      lo(:) = dg_frag%ixyz_frag(:, cand)
      hi(:) = dg_frag%ixyz_frag(:, cand) + dg_frag%nxyz_domain(:, cand) - 1
      if (side > 0) then
        cand_face = modulo(lo(axis) - 1, max(1, dg_frag%lgnum_total(axis))) + 1
      else
        cand_face = modulo(hi(axis) - 1, max(1, dg_frag%lgnum_total(axis))) + 1
      end if
      if (cand_face /= target) cycle
      if (.not. wpw_local_is_neighbor_pair(dg_frag, ifrag, cand)) cycle
      jfrag = cand
      return
    end do
  end function wpw_face_neighbor_fragment


  logical function wpw_local_is_neighbor_pair(dg_frag, ifrag, jfrag) result(is_neighbor)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, jfrag

    is_neighbor = .false.
    if (ifrag == jfrag) then
      is_neighbor = .true.
      return
    end if
    is_neighbor = wpw_local_axis_neighbor(dg_frag, 1, ifrag, jfrag) .and. &
      wpw_local_axis_neighbor(dg_frag, 2, ifrag, jfrag) .and. &
      wpw_local_axis_neighbor(dg_frag, 3, ifrag, jfrag)
  end function wpw_local_is_neighbor_pair


  logical function wpw_local_axis_neighbor(dg_frag, axis, ifrag, jfrag) result(ok)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: axis, ifrag, jfrag
    integer :: s1, s2, n1, n2, e1, e2, next1, next2, lg

    s1 = dg_frag%ixyz_frag(axis, ifrag)
    s2 = dg_frag%ixyz_frag(axis, jfrag)
    n1 = dg_frag%nxyz_domain(axis, ifrag)
    n2 = dg_frag%nxyz_domain(axis, jfrag)
    lg = max(1, dg_frag%lgnum_total(axis))
    e1 = s1 + n1 - 1
    e2 = s2 + n2 - 1
    next1 = modulo(e1, lg) + 1
    next2 = modulo(e2, lg) + 1
    ok = ((s1 == s2) .and. (n1 == n2)) .or. (s1 == next2) .or. (s2 == next1)
  end function wpw_local_axis_neighbor


  pure integer function map_global_to_phi_box_coord_pw(ig, phi_lo, phi_hi, lgnum) result(local_idx)
    implicit none
    integer, intent(in) :: ig, phi_lo, phi_hi, lgnum

    local_idx = modulo(ig - 1, lgnum) + 1
    if (local_idx < phi_lo) local_idx = local_idx + lgnum
    if (local_idx > phi_hi) local_idx = local_idx - lgnum
  end function map_global_to_phi_box_coord_pw

end module rt_dg_wpw_window
