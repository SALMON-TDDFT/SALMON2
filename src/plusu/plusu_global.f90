module plusU_global

  implicit none
  private
  public :: PLUS_U_ON
  public :: read_Hubbard_parameters
  public :: dm_plusU

  logical :: PLUS_U_ON = .true.

  real(8),allocatable :: U_Hubbard(:), J_Hubbard(:), U_eff(:)
  integer,allocatable :: orb_info_Hubbard(:,:)

  real(8),allocatable :: dm_plusU(:)

  logical :: read_parameters=.false.

contains

  subroutine read_Hubbard_parameters
    implicit none
    integer,parameter :: unit=100
    integer :: nlines,i
    if ( read_parameters ) return
    open(100,file="hubbard_param.dat",status="old")
    nlines=0
    do
      read(unit,*,END=9)
      nlines=nlines+1
    end do
9   continue
    allocate( U_Hubbard(nlines) ); U_Hubbard=0.0d0
    allocate( J_Hubbard(nlines) ); J_Hubbard=0.0d0
    allocate( orb_info_Hubbard(3,nlines) ); orb_info_Hubbard=0
    rewind unit
    do i=1,nlines
      read(unit,*) orb_info_Hubbard(1:3,i), U_Hubbard(i), J_Hubbard(i)
      write(*,'(1x,i2,2x,3i4,2f12.5)') i, orb_info_Hubbard(1:3,i), U_Hubbard(i), J_Hubbard(i)
    end do
    close(unit)
    read_parameters=.true.
  end subroutine read_Hubbard_parameters

end module plusU_global
