program test_dg_wpw_windows
  use dg_wpw_windows, only: evaluate_dg_wpw_normalized_windows
  implicit none
  integer, parameter :: nsupport = 2
  integer :: core_lo(3, nsupport), core_hi(3, nsupport)
  integer :: total_grid(3), buffer(3), width(3), info, g, k
  real(8) :: hgs(3), chi(nsupport), grad_chi(3, nsupport)
  real(8) :: chi_periodic(nsupport), grad_periodic(3, nsupport)
  real(8) :: partition, derivative_identity(3)
  real(8), parameter :: tol = 2.0d-13

  core_lo(:, 1) = [1, 1, 1]
  core_hi(:, 1) = [4, 1, 1]
  core_lo(:, 2) = [5, 1, 1]
  core_hi(:, 2) = [8, 1, 1]
  total_grid = [8, 1, 1]
  hgs = [0.5d0, 1.0d0, 1.0d0]
  buffer = [2, 0, 0]
  width = [2, 0, 0]

  do g = 1, total_grid(1)
    call evaluate_dg_wpw_normalized_windows(core_lo, core_hi, total_grid, hgs, buffer, width, &
      [g, 1, 1], chi, grad_chi, info)
    call require(info == 0, 'valid support-local window evaluation failed')
    partition = sum(chi**2)
    derivative_identity = 0.0d0
    do k = 1, nsupport
      derivative_identity = derivative_identity + chi(k) * grad_chi(:, k)
    end do
    call require(abs(partition - 1.0d0) < tol, 'sum chi_K^2 is not one')
    call require(maxval(abs(derivative_identity)) < tol, 'sum chi_K grad chi_K is not zero')
  end do

  call evaluate_dg_wpw_normalized_windows(core_lo, core_hi, total_grid, hgs, buffer, width, &
    [1, 1, 1], chi, grad_chi, info)
  call require(info == 0, 'periodic reference evaluation failed')
  call evaluate_dg_wpw_normalized_windows(core_lo, core_hi, total_grid, hgs, buffer, width, &
    [9, 1, 1], chi_periodic, grad_periodic, info)
  call require(info == 0, 'periodic image evaluation failed')
  call require(maxval(abs(chi - chi_periodic)) < tol, 'periodic window values disagree')
  call require(maxval(abs(grad_chi - grad_periodic)) < tol, 'periodic window gradients disagree')

  call evaluate_dg_wpw_normalized_windows(core_lo(:, 1:1), core_hi(:, 1:1), total_grid, hgs, buffer, width, &
    [7, 1, 1], chi(1:1), grad_chi(:, 1:1), info)
  call require(info /= 0, 'missing support coverage did not fail closed')

  print *, 'PASS support-local normalized periodic WPW windows'

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // trim(message)
      error stop 1
    end if
  end subroutine require

end program test_dg_wpw_windows
