module spin_orbit_global

  implicit none

  private
  public :: read_sw_spin_orbit

  logical, public :: SPIN_ORBIT_ON = .false.

contains

  subroutine read_sw_spin_orbit( yn )
    implicit none
    character(1), intent(in) :: yn
    SPIN_ORBIT_ON = ( yn == 'y' )
!    logical :: flag
!    inquire( file='so.dat', exist=flag )
!    if( flag )then
!      open(100,file='so.dat')
!      read(100,*) SPIN_ORBIT_ON
!      close(100)
!    end if
  end subroutine read_sw_spin_orbit

end module spin_orbit_global
