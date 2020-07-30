module spin_orbit_global

  implicit none

  logical :: SPIN_ORBIT_ON = .false.

contains

  subroutine read_sw_spin_orbit
    implicit none
    logical :: flag
    inquire( file='so.dat', exist=flag )
    if( flag )then
      open(100,file='so.dat')
      read(100,*) SPIN_ORBIT_ON
      close(100)
    end if
  end subroutine read_sw_spin_orbit

end module spin_orbit_global
