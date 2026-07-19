!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!
module dg_wpw_windows
  implicit none
  private
  public :: evaluate_dg_wpw_normalized_windows

contains

  subroutine evaluate_dg_wpw_normalized_windows(core_lo, core_hi, total_grid, hgs, buffer, width, &
      grid, chi, grad_chi, info)
    integer, intent(in) :: core_lo(:,:), core_hi(:,:), total_grid(3), buffer(3), width(3), grid(3)
    real(8), intent(in) :: hgs(3)
    real(8), intent(out) :: chi(:), grad_chi(:,:)
    integer, intent(out) :: info
    integer :: k, nsupport
    real(8), allocatable :: q(:), grad_q(:,:)
    real(8) :: qsum, sqrt_qsum, q_grad_sum(3)
    real(8), parameter :: tiny_q = 1.0d-28

    info = 1
    chi = 0.0d0
    grad_chi = 0.0d0
    nsupport = size(core_lo, 2)
    if (nsupport < 1) return
    if (size(core_lo, 1) /= 3 .or. any(shape(core_hi) /= shape(core_lo))) return
    if (size(chi) /= nsupport) return
    if (size(grad_chi, 1) /= 3 .or. size(grad_chi, 2) /= nsupport) return
    if (any(total_grid < 1) .or. any(hgs <= 0.0d0)) return
    if (any(buffer < 0) .or. any(width < 0) .or. any(width > buffer)) return
    if (any(core_hi < core_lo)) return

    allocate(q(nsupport), grad_q(3, nsupport))
    do k = 1, nsupport
      call evaluate_raw_window(core_lo(:, k), core_hi(:, k), total_grid, hgs, width, grid, &
        q(k), grad_q(:, k))
    end do
    qsum = sum(q**2)
    if (qsum <= tiny_q) then
      deallocate(q, grad_q)
      return
    end if

    sqrt_qsum = sqrt(qsum)
    q_grad_sum = 0.0d0
    do k = 1, nsupport
      q_grad_sum = q_grad_sum + q(k) * grad_q(:, k)
    end do
    do k = 1, nsupport
      chi(k) = q(k) / sqrt_qsum
      grad_chi(:, k) = grad_q(:, k) / sqrt_qsum - q(k) * q_grad_sum / (qsum * sqrt_qsum)
    end do
    info = 0
    deallocate(q, grad_q)
  end subroutine evaluate_dg_wpw_normalized_windows


  subroutine evaluate_raw_window(core_lo, core_hi, total_grid, hgs, width, grid, q, grad_q)
    integer, intent(in) :: core_lo(3), core_hi(3), total_grid(3), width(3), grid(3)
    real(8), intent(in) :: hgs(3)
    real(8), intent(out) :: q, grad_q(3)
    real(8) :: w(3), dw(3)
    integer :: axis

    do axis = 1, 3
      call evaluate_axis_window(core_lo(axis), core_hi(axis), total_grid(axis), hgs(axis), &
        width(axis), grid(axis), w(axis), dw(axis))
    end do
    q = product(w)
    grad_q(1) = dw(1) * w(2) * w(3)
    grad_q(2) = w(1) * dw(2) * w(3)
    grad_q(3) = w(1) * w(2) * dw(3)
  end subroutine evaluate_raw_window


  subroutine evaluate_axis_window(core_lo, core_hi, total_grid, h, width, grid, w, dw)
    integer, intent(in) :: core_lo, core_hi, total_grid, width, grid
    real(8), intent(in) :: h
    real(8), intent(out) :: w, dw
    integer :: shift0, shift, lo, hi
    real(8) :: t, arg
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0)

    w = 0.0d0
    dw = 0.0d0
    shift0 = nint(dble(grid - core_lo) / dble(total_grid))
    do shift = shift0 - 1, shift0 + 1
      lo = core_lo + shift * total_grid
      hi = core_hi + shift * total_grid
      if (grid >= lo .and. grid <= hi) then
        w = 1.0d0
        return
      end if
      if (width <= 0) cycle
      if (grid >= lo - width .and. grid < lo) then
        t = dble(grid - (lo - width)) / dble(width)
        arg = 0.5d0 * pi * max(0.0d0, min(1.0d0, t))
        w = sin(arg)
        dw = 0.5d0 * pi * cos(arg) / (dble(width) * h)
        return
      end if
      if (grid > hi .and. grid <= hi + width) then
        t = dble(grid - hi) / dble(width)
        arg = 0.5d0 * pi * max(0.0d0, min(1.0d0, t))
        w = cos(arg)
        dw = -0.5d0 * pi * sin(arg) / (dble(width) * h)
        return
      end if
    end do
  end subroutine evaluate_axis_window

end module dg_wpw_windows
