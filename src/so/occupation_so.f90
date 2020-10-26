module occupation_so

  use spin_orbit_global, only: SPIN_ORBIT_ON

  implicit none

  private
  public :: init_occupation_so
  public :: SPIN_ORBIT_ON

contains

  subroutine init_occupation_so( occ, nelectrons )
    implicit none
    real(8),intent(inout) :: occ(:,:,:) ! (io,ik,is)
    integer,intent(in) :: nelectrons
    occ = 0.0d0
    occ(1:nelectrons,:,:)=1.0d0
  end subroutine init_occupation_so

end module occupation_so
