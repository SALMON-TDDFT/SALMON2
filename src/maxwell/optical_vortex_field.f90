module optical_vortex_field
  implicit none

contains

subroutine calc_optical_vortex_Ac_ext(t, x, y, nx, ny, hgs, A_ext)
  use math_constants, only: zi, pi
  use salmon_global, only: I_wcm2_1, E_amplitude1, ae_shape1, tw1, omega1, phi_cep1, t1_start, &
    & optical_vortex_center_x, optical_vortex_center_y, optical_vortex_radius, &
    & optical_vortex_charge, optical_vortex_polarization
  implicit none
  real(8), intent(in) :: t, x, y
  integer, intent(in) :: nx, ny
  real(8), intent(in) :: hgs(3)
  real(8), intent(out) :: A_ext(3)
  integer :: npower
  real(8) :: f0, tt, center_x, center_y, rxy, phi_xy, window, radial_factor, phase
  complex(8) :: zpol(3), zphase

  A_ext = 0d0

  if (I_wcm2_1 < 0d0) then
    f0 = E_amplitude1
  else
    f0 = 5.338d-9 * sqrt(I_wcm2_1)
  end if
  if (f0 == 0d0) return

  center_x = optical_vortex_center_x
  center_y = optical_vortex_center_y
  if (center_x < -1d20) center_x = 0.5d0 * nx * hgs(1)
  if (center_y < -1d20) center_y = 0.5d0 * ny * hgs(2)

  rxy = sqrt((x - center_x)**2 + (y - center_y)**2)
  if (rxy >= optical_vortex_radius) return

  phi_xy = atan2(y - center_y, x - center_x)
  tt = t - 0.5d0 * tw1 - t1_start

  select case(trim(ae_shape1))
  case('Acos2'); npower = 2
  case('Acos3'); npower = 3
  case('Acos4'); npower = 4
  case('Acos6'); npower = 6
  case('Acos8'); npower = 8
  case default
    stop 'unsupported ae_shape1 in calc_optical_vortex_Ac_ext'
  end select

  if (abs(tt) >= 0.5d0 * tw1) return

  window = sin(0.5d0 * pi * (1d0 - rxy / optical_vortex_radius))**2
  if (optical_vortex_charge == 0) then
    radial_factor = 1d0
  else
    radial_factor = (rxy / optical_vortex_radius)**abs(optical_vortex_charge)
  end if

  zpol = dcmplx(0d0, 0d0)
  select case(trim(optical_vortex_polarization))
  case('linear_x')
    zpol(1) = dcmplx(1d0, 0d0)
  case('linear_y')
    zpol(2) = dcmplx(1d0, 0d0)
  case('left_circular')
    zpol(1) = dcmplx(1d0 / sqrt(2d0), 0d0)
    zpol(2) = dcmplx(0d0, -1d0 / sqrt(2d0))
  case('right_circular')
    zpol(1) = dcmplx(1d0 / sqrt(2d0), 0d0)
    zpol(2) = dcmplx(0d0, 1d0 / sqrt(2d0))
  case default
    stop 'unsupported optical_vortex_polarization in calc_optical_vortex_Ac_ext'
  end select

  phase = omega1 * tt + phi_cep1 * 2d0 * pi + dble(optical_vortex_charge) * phi_xy
  zphase = exp(zi * phase)
  A_ext(:) = -f0 / omega1 * (cos(pi * tt / tw1))**npower * window * radial_factor * aimag(zpol(:) * zphase)
end subroutine calc_optical_vortex_Ac_ext

subroutine calc_optical_vortex_EB_ext(t, r, x, y, nx, ny, hgs, dt, E_ext, B_ext)
  use phys_constants, only: cspeed_au
  implicit none
  real(8), intent(in) :: t, r, x, y, hgs(3), dt
  integer, intent(in) :: nx, ny
  real(8), intent(out) :: E_ext(3), B_ext(3)
  real(8) :: Ac_t_plus(3), Ac_t_minus(3), Ac_x_plus(3), Ac_x_minus(3), Ac_y_plus(3), Ac_y_minus(3)
  real(8) :: dA_dt(3), dA_dx(3), dA_dy(3), tt

  tt = t - r / cspeed_au

  call calc_optical_vortex_Ac_ext(tt + 0.5d0 * dt, x, y, nx, ny, hgs, Ac_t_plus)
  call calc_optical_vortex_Ac_ext(tt - 0.5d0 * dt, x, y, nx, ny, hgs, Ac_t_minus)
  call calc_optical_vortex_Ac_ext(tt, x + 0.5d0 * hgs(1), y, nx, ny, hgs, Ac_x_plus)
  call calc_optical_vortex_Ac_ext(tt, x - 0.5d0 * hgs(1), y, nx, ny, hgs, Ac_x_minus)
  call calc_optical_vortex_Ac_ext(tt, x, y + 0.5d0 * hgs(2), nx, ny, hgs, Ac_y_plus)
  call calc_optical_vortex_Ac_ext(tt, x, y - 0.5d0 * hgs(2), nx, ny, hgs, Ac_y_minus)

  dA_dt(:) = (Ac_t_plus(:) - Ac_t_minus(:)) / dt
  dA_dx(:) = (Ac_x_plus(:) - Ac_x_minus(:)) / hgs(1)
  dA_dy(:) = (Ac_y_plus(:) - Ac_y_minus(:)) / hgs(2)

  E_ext(:) = -dA_dt(:)
  ! Plane-wave propagation along +z: dA/dz ≈ (1/c) dA/dt, so B_x ≈ -E_y and B_y ≈ E_x.
  ! For vortex fields with A_z=0, curl A gives B_x=B_y=0, which zeroes the Poynting
  ! vector transverse components and breaks the Lz_inc/ref flux decomposition.
  ! The plane-wave approximation is physically correct for the OAM flux calculation.
  B_ext(1) = -E_ext(2)
  B_ext(2) =  E_ext(1)
  B_ext(3) = cspeed_au * (dA_dx(2) - dA_dy(1))
end subroutine calc_optical_vortex_EB_ext

end module optical_vortex_field
