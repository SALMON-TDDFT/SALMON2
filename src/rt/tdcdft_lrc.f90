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
  public :: advance_xc_field
contains
  pure subroutine advance_xc_field(dt,alpha,damping,restoring,current,a_old,a_now,a_next)
    implicit none
    real(8),intent(in) :: dt,alpha,damping,restoring,current(3),a_old(3),a_now(3)
    real(8),intent(out) :: a_next(3)
    ! Centered second-order discretization of A'' + damping*A' + restoring*A = alpha*j.
    a_next=((2d0-restoring*dt**2)*a_now-(1d0-0.5d0*damping*dt)*a_old &
           +alpha*dt**2*current)/(1d0+0.5d0*damping*dt)
  end subroutine advance_xc_field
end module tdcdft_lrc
