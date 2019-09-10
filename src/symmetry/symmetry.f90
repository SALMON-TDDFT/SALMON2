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


  subroutine init_sym_sub( Amat, Bmat )
    implicit none
    real(8),intent(in) :: Amat(3,3), Bmat(3,3) ! Lattice vectors
    real(8) :: Ainv(3,3), Binv(3,3), tmpmat(3,3), pi2
    integer :: nsym, isym
    write(*,'(a60)') repeat("-",40)//" init_sym_sub(start)"
    call read_SymMat
    nsym=size(SymMatR,3)
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
