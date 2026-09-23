! Copyright 2026 SALMON developers
! Licensed under the Apache License, Version 2.0.
! See http://www.apache.org/licenses/LICENSE-2.0 for the license terms.
! Distributed on an AS IS basis, without warranties or conditions of any kind.
!
! Macroscopic xc vector potential, in SALMON A/c and electron-number-current conventions.
! See Sun et al., PRL 127, 077401 (2021), and Williams and Ullrich,
! JCTC 21, 4753 (2025), Eq. (26). All arguments use atomic units.
module tdcdft_lrc
  implicit none
  private
  public :: advance_xc_field,proca_coefficients,instant_screening
contains
  pure subroutine instant_screening(a,e,j,p,omega,reference,strength,field_floor,alpha0,response,alpha)
    implicit none
    real(8),intent(in) :: a(3),e(3),j(3),p(3),omega,reference,strength,field_floor,alpha0
    real(8),intent(inout) :: response,alpha
    real(8) :: norm2,increment
    ! Instantaneous vector least squares, not a time average or a dielectric inversion.
    norm2=sum(a*a)+sum((e/omega)**2)
    if(norm2<=field_floor**2) return
    response=(sum(a*j)-sum(e*p))/norm2
    increment=max(response-reference,0d0)
    alpha=alpha0/(1d0+strength*4d0*acos(-1d0)*increment/omega**2)
  end subroutine instant_screening

  subroutine proca_coefficients(a2,a0,alpha,restoring)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: a2,a0
    real(8),intent(out) :: alpha,restoring
    ! Paper uses (p-A/c)^2/2; SALMON uses (p+A/c)^2/2 and number-current density.
    if (.not.all(ieee_is_finite([a2,a0]))) error stop 'TDCDFT: a2 and a0 must be finite'
    if (abs(a2)<tiny(1d0)) error stop 'TDCDFT: a2 must be nonzero'
    alpha=-4d0*acos(-1d0)/a2
    restoring=a0/a2
    if (.not.all(ieee_is_finite([alpha,restoring]))) error stop 'TDCDFT: Proca coefficients overflow'
    if (restoring<0d0) error stop 'TDCDFT: require a0/a2 >= 0'
  end subroutine proca_coefficients

  pure subroutine advance_xc_field(dt,alpha,damping,restoring,current,a_old,a_now,a_next)
    implicit none
    real(8),intent(in) :: dt,alpha,damping,restoring,current(3),a_old(3),a_now(3)
    real(8),intent(out) :: a_next(3)
    ! Centered second-order discretization of A'' + damping*A' + restoring*A = alpha*j.
    a_next=((2d0-restoring*dt**2)*a_now-(1d0-0.5d0*damping*dt)*a_old &
           +alpha*dt**2*current)/(1d0+0.5d0*damping*dt)
  end subroutine advance_xc_field
end module tdcdft_lrc
