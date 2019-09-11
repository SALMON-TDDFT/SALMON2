module sym_sub

  implicit none

  private
  public :: init_sym_sub

  logical,parameter,public :: use_symmetry=.false.
  !logical,parameter,public :: use_symmetry=.true.

  real(8),allocatable :: SymMatR(:,:,:)
  real(8),allocatable,public :: SymMatA(:,:,:)
  real(8),allocatable,public :: SymMatB(:,:,:)

contains


  subroutine init_sym_sub( Amat, Bmat, epdir )
    implicit none
    real(8),intent(in) :: Amat(3,3), Bmat(3,3) ! Lattice vectors
    real(8),intent(in) :: epdir(3)
    real(8) :: Ainv(3,3), Binv(3,3), tmpmat(3,3), pi2
    real(8),allocatable :: work(:,:,:)
    integer :: nsym, isym, n
    logical :: ok(3)

    write(*,'(a60)') repeat("-",40)//" init_sym_sub(start)"

    call read_SymMat
    nsym=size(SymMatR,3)

    allocate( work(3,4,nsym) ); work=0.0d0

    n=0
    do isym=1,nsym
       ok=.true.
       if ( epdir(1)/=0.0d0 ) then
          ok(1)=.false.
          if ( SymMatR(1,1,isym)==1.0d0 .and. SymMatR(2,1,isym)==0.0d0 .and. SymMatR(3,1,isym)==0.0d0 ) ok(1)=.true.
       end if
       if ( epdir(2)/=0.0d0 ) then
          ok(2)=.false.
          if ( SymMatR(1,2,isym)==0.0d0 .and. SymMatR(2,2,isym)==1.0d0 .and. SymMatR(3,2,isym)==0.0d0 ) ok(2)=.true.
       end if
       if ( epdir(3)/=0.0d0 ) then
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

    do isym=1,nsym
       write(*,'(1x,i4,3f10.5,2x,f10.5)') isym,SymMatR(1,:,isym)
       write(*,'(1x,4x,3f10.5,2x,f10.5)')      SymMatR(2,:,isym)
       write(*,'(1x,4x,3f10.5,2x,f10.5)')      SymMatR(3,:,isym)
    end do

    deallocate( work )

! ---

    allocate( SymMatA(3,4,nsym) ); SymMatA=0.0d0
    allocate( SymMatB(3,4,nsym) ); SymMatB=0.0d0
    pi2=2.0d0*acos(-1.0d0)
    Ainv=transpose(Bmat)/pi2
    Binv=transpose(Amat)/pi2
    do isym=1,nsym
       tmpmat=matmul( SymMatR(1:3,1:3,isym), Amat )
       SymMatA(1:3,1:3,isym)=matmul( Ainv, tmpmat )
       SymMatA(1:3,4,isym)=SymMatR(1:3,4,isym)
    end do
    do isym=1,nsym
       tmpmat=matmul( SymMatR(1:3,1:3,isym), Bmat )
       SymMatB(1:3,1:3,isym)=matmul( Binv, tmpmat )
       SymMatB(1:3,4,isym)=SymMatR(1:3,4,isym)
    end do
    write(*,'(a60)') repeat("-",42)//" init_sym_sub(end)"
  end subroutine init_sym_sub


  subroutine read_SymMat
    implicit none
    integer,parameter :: unit=1001
    integer :: i, j, nsym

    open( unit, file='cif_sym.dat', status='old' )

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
