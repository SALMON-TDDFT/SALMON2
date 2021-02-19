module sym_sub

  use communication, only: comm_get_globalinfo, comm_is_root

  implicit none

  private
  public :: read_sw_symmetry
  public :: init_sym_sub

  logical,public :: DISPLAY     =.false.
  logical,public :: use_symmetry=.false.
  logical :: use_symmetry_dir(3) = .false.

  character(8)   :: sym_file    ='sym.dat'
  real(8),allocatable :: SymMatR(:,:,:)
  real(8),allocatable,public :: SymMatA(:,:,:)
  real(8),allocatable,public :: SymMatB(:,:,:)
  real(8),public :: Amat(3,3), Ainv(3,3) ! Each column of Amat (Bmat)
  real(8),public :: Bmat(3,3), Binv(3,3) !   is the (reciprocal) lattice vector
  logical :: flag_init=.false.

contains


  subroutine read_sw_symmetry( yn )
    implicit none
    character(*),intent(in) :: yn
    integer :: n,i
    if ( index(yn,'y') /= 0 ) use_symmetry_dir(:)=.true.
    n=len(trim(yn))
    do i = 1, n
      if ( yn(i:i) == 'n' ) use_symmetry_dir(i) = .false.
    end do
    use_symmetry = any( use_symmetry_dir )
  end subroutine read_sw_symmetry

  subroutine init_sym_sub( Amat_in, Bmat_in )
    implicit none
    real(8),intent(in) :: Amat_in(3,3), Bmat_in(3,3) ! Lattice vectors
    real(8) :: tmpmat(3,3), pi2
    real(8),allocatable :: work(:,:,:)
    integer :: ngid, npid, nprocs
    integer :: nsym, isym, n, j
    logical :: ok(3)

    if ( .not.use_symmetry ) return

    if ( flag_init ) return

    call comm_get_globalinfo( ngid, npid, nprocs )
    DISPLAY = comm_is_root(npid)

    if ( DISPLAY ) write(*,'(a60)') repeat("-",40)//" init_sym_sub(start)"

    call read_SymMat( use_symmetry )
    nsym=size(SymMatR,3)

    allocate( work(3,4,nsym) ); work=0.0d0

    n=0
    do isym=1,nsym
       ok=.true.
       if ( .not.use_symmetry_dir(1) ) then
          ok(1)=.false.
          if ( SymMatR(1,1,isym)==1.0d0 .and. SymMatR(2,1,isym)==0.0d0 .and. SymMatR(3,1,isym)==0.0d0 ) ok(1)=.true.
       end if
       if ( .not.use_symmetry_dir(2) ) then
          ok(2)=.false.
          if ( SymMatR(1,2,isym)==0.0d0 .and. SymMatR(2,2,isym)==1.0d0 .and. SymMatR(3,2,isym)==0.0d0 ) ok(2)=.true.
       end if
       if ( .not.use_symmetry_dir(3) ) then
          ok(3)=.false.
          if ( SymMatR(1,3,isym)==0.0d0 .and. SymMatR(2,3,isym)==0.0d0 .and. SymMatR(3,3,isym)==1.0d0 ) ok(3)=.true.
       end if
       if ( all(ok) ) then
          n=n+1
          work(:,:,n)=SymMatR(:,:,isym)
       end if
    end do

    nsym=n
    SymMatR=0.0d0
    SymMatR(:,:,1:nsym)=work(:,:,1:nsym)

    if ( DISPLAY ) then
       do isym=1,nsym
         write(*,'(1x,i4,3f10.5,2x,f10.5)') isym,(SymMatR(1,j,isym),j=1,4)
         write(*,'(1x,4x,3f10.5,2x,f10.5)')      (SymMatR(2,j,isym),j=1,4)
         write(*,'(1x,4x,3f10.5,2x,f10.5)')      (SymMatR(3,j,isym),j=1,4)
       end do
    end if

    deallocate( work )

! ---

    Amat=Amat_in
    Bmat=Bmat_in

    allocate( SymMatA(3,4,nsym) ); SymMatA=0.0d0
    allocate( SymMatB(3,4,nsym) ); SymMatB=0.0d0
    pi2=2.0d0*acos(-1.0d0)
    Ainv=transpose(Bmat)/pi2
    Binv=transpose(Amat)/pi2
    SymMatA(:,:,1:nsym)=SymMatR(:,:,1:nsym)
    do isym=1,nsym
       tmpmat=matmul( SymMatA(1:3,1:3,isym), Ainv )
       SymMatR(1:3,1:3,isym)=matmul( Amat, tmpmat )
       tmpmat=matmul( SymMatR(1:3,1:3,isym), Bmat )
       SymMatB(1:3,1:3,isym)=matmul( Binv, tmpmat )
       SymMatB(1:3,4,isym)=SymMatR(1:3,4,isym)
    end do
    flag_init=.true.
    if ( DISPLAY ) write(*,'(a60)') repeat("-",42)//" init_sym_sub(end)"
  end subroutine init_sym_sub


  subroutine read_SymMat( flag )
    implicit none
    logical,intent(out) :: flag
    integer,parameter :: unit=1001
    integer :: i, j, nsym

    inquire( FILE=sym_file, EXIST=flag )
    if ( .not.flag ) then
       if ( DISPLAY ) write(*,*) "symmetry-operation file ( "//sym_file//" ) can not be found."
       stop 'stop@read_SymMat'
    else
       if ( DISPLAY ) write(*,*) "symmetry-operation file is found ( "//sym_file//" )."
    end if

    open( unit, file='sym.dat', status='old' )

    i=0
    do
       read(unit,*,END=9)
       i=i+1
    end do
9   nsym=i/3

!--------------------------------------------------------------
! 3x3-marix data is stored in SymMatR(1:3,1:3,?)
! and the associated shift vector is stored in SymMatR(1:3,4,?)
!
    allocate( SymMatR(3,4,nsym) ); SymMatR=0.0d0
!
!--------------------------------------------------------------

    rewind unit
    do i=1,nsym
       read(unit,*) (SymMatR(1,j,i),j=1,4)
       read(unit,*) (SymMatR(2,j,i),j=1,4)
       read(unit,*) (SymMatR(3,j,i),j=1,4)
    end do

    close( unit )

  end subroutine read_SymMat


end module sym_sub
