module plusU_global

  use salmon_communication, only: comm_get_globalinfo

  implicit none
  private
  public :: PLUS_U_ON
  public :: read_Hubbard_parameters
  public :: prep_Hubbard_parameters
  public :: dm_mms_nla
  public :: U_eff
  public :: V_eff

  logical :: PLUS_U_ON = .false.
  complex(8),allocatable :: dm_mms_nla(:,:,:,:,:,:) ! density matrix
  real(8),allocatable :: U_eff(:,:,:)
  real(8),allocatable :: V_eff(:,:)

  integer,allocatable :: orb_info_Hubbard(:,:)
  real(8),allocatable :: U_Hubbard(:)
  real(8),allocatable :: J_Hubbard(:)
  logical :: read_parameters=.false.
  integer :: npid

contains

  subroutine read_Hubbard_parameters
    implicit none
    integer,parameter :: unit=100
    integer :: nlines,i,ngid,nprocs
    character(2) :: unit_of_parameters
    real(8),parameter :: au2eV=27.2116d0
    if ( read_parameters ) return
    call comm_get_globalinfo( ngid, npid, nprocs )
    call check_plusU_ON_OFF
    if ( .not.PLUS_U_ON ) return 
    open(100,file="hubbard_param.dat",status="old")
    read(unit,*,END=9) ! the 1st line is assumed to be true/false switch which is already read
    read(unit,*,END=9) unit_of_parameters
    nlines=0
    do
      read(unit,*,END=9)
      nlines=nlines+1
    end do
9   continue
!
    allocate( U_Hubbard(nlines) ); U_Hubbard=0.0d0
    allocate( J_Hubbard(nlines) ); J_Hubbard=0.0d0
    allocate( orb_info_Hubbard(3,nlines) ); orb_info_Hubbard=0
!
    rewind unit
    read(unit,*)  ! skip
    read(unit,*)  ! skip
    if( npid == 0 )then
      write(*,'(1x,2x,2x,3a4,2a12,3x,"(unit is in ",a2,")")') "a","l","n","U","J",unit_of_parameters
    end if
    do i=1,nlines
      read(unit,*) orb_info_Hubbard(1:3,i), U_Hubbard(i), J_Hubbard(i) 
      if( npid == 0 )then
        write(*,'(1x,i2,2x,3i4,2f12.5)') i, orb_info_Hubbard(1:3,i), U_Hubbard(i), J_Hubbard(i)
      end if
    end do
    close(unit)
    if( unit_of_parameters=="eV" .or. unit_of_parameters=="ev" )then
      U_Hubbard(1:nlines)=U_Hubbard(1:nlines)/au2eV
      J_Hubbard(1:nlines)=J_Hubbard(1:nlines)/au2eV
    end if
    read_parameters=.true.
  end subroutine read_Hubbard_parameters

  subroutine check_plusU_ON_OFF
    implicit none
    logical :: flag
    inquire( FILE='hubbard_param.dat', EXIST=flag )
    if ( flag ) then
      open(100,file='hubbard_param.dat',status='old')
      read(100,*) PLUS_U_ON
      close(100)
      if( npid == 0 )then
        write(*,*) "hubbard_param.dat is found: PLUS_U_ON=",PLUS_U_ON
        if ( PLUS_U_ON ) write(*,*) "+U calculation is performed."
      end if
    end if
  end subroutine check_plusU_ON_OFF

  subroutine prep_Hubbard_parameters( U_eff )
    implicit none
    real(8),intent(out) :: U_eff(:,0:,:)
    integer :: i,a,l,n
    U_eff=0.0d0
    do i=1,size(J_Hubbard)
      a=orb_info_Hubbard(1,i)
      l=orb_info_Hubbard(2,i)
      n=orb_info_Hubbard(3,i)
      U_eff(a,l,n) = U_Hubbard(i) - J_Hubbard(i)
    end do
  end subroutine prep_Hubbard_parameters

end module plusU_global
