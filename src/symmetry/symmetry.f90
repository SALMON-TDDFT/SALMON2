module sym_sub

  use communication, only: comm_get_globalinfo, comm_is_root

  implicit none

  private
  public :: read_sw_symmetry
  public :: init_sym_sub, symmetry_validate_field, symmetry_validate_atoms_cartesian

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
    use salmon_global, only: xc, tdcdft, theory
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
    if (xc=='hse06'.or.tdcdft/='none') call symmetry_validate_group()
    flag_init=.true.
    ! Test the retained Cartesian group, not just the requested axis flags.
    ! The bounded HSE/TDCDFT RT implementation supports the z-field little group.
    if ((xc=='hse06'.or.tdcdft/='none').and. &
        (theory=='tddft_pulse'.or.theory=='tddft_response'.or.theory=='tddft')) &
      call symmetry_validate_field([0d0,0d0,1d0])
    if ( DISPLAY ) write(*,'(a60)') repeat("-",42)//" init_sym_sub(end)"
  end subroutine init_sym_sub


  subroutine symmetry_validate_atoms_cartesian(rion,kion)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(8),intent(in) :: rion(:,:)
    integer,intent(in) :: kion(:)
    real(8),allocatable :: fractional(:,:)
    logical,allocatable :: used(:)
    real(8) :: transformed(3),delta(3)
    integer :: a,b,s,na
    logical :: found
    if (.not.use_symmetry) return
    if (.not.flag_init) error stop 'Symmetry atom check before initialization'
    na=size(kion)
    if(size(rion,1)/=3.or.size(rion,2)/=na) error stop 'Symmetry: invalid atom layout'
    if(.not.all(ieee_is_finite(rion))) error stop 'Symmetry: nonfinite atom position'
    allocate(fractional(3,na),used(na))
    fractional=matmul(Ainv,rion)
    do s=1,size(SymMatA,3)
      used=.false.
      do a=1,na
        transformed=matmul(SymMatA(:,1:3,s),fractional(:,a))+SymMatA(:,4,s)
        found=.false.
        do b=1,na
          if(used(b).or.kion(a)/=kion(b))cycle
          delta=transformed-fractional(:,b)
          if(maxval(abs(delta-anint(delta)))<1d-8)then
            used(b)=.true.;found=.true.;exit
          endif
        enddo
        if(.not.found)error stop 'Symmetry: operation does not preserve atomic positions and species'
      enddo
    enddo
  end subroutine symmetry_validate_atoms_cartesian

  subroutine symmetry_validate_group()
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(8) :: rotation(3,3),translation(3),delta(3),cartesian(3,3),identity(3,3),determinant
    integer :: a,b,c,nsym
    logical :: found
    nsym=size(SymMatA,3)
    if (nsym==0) error stop 'Symmetry: empty retained group'
    if (.not.all(ieee_is_finite(SymMatA))) error stop 'Symmetry: nonfinite operation'
    identity=0d0
    do a=1,3
      identity(a,a)=1d0
    enddo
    do a=1,nsym
      rotation=SymMatA(:,1:3,a)
      if(maxval(abs(rotation-anint(rotation)))>1d-10) &
        error stop 'Symmetry: rotation does not preserve the fractional lattice'
      determinant=rotation(1,1)*(rotation(2,2)*rotation(3,3)-rotation(2,3)*rotation(3,2)) &
        -rotation(1,2)*(rotation(2,1)*rotation(3,3)-rotation(2,3)*rotation(3,1)) &
        +rotation(1,3)*(rotation(2,1)*rotation(3,2)-rotation(2,2)*rotation(3,1))
      if(abs(abs(determinant)-1d0)>1d-10) error stop 'Symmetry: nonunimodular lattice rotation'
      cartesian=matmul(Amat,matmul(rotation,Ainv))
      if(.not.all(ieee_is_finite(cartesian))) error stop 'Symmetry: nonfinite Cartesian rotation'
      if(maxval(abs(matmul(transpose(cartesian),cartesian)-identity))>1d-10) &
        error stop 'Symmetry: operation is not a Cartesian isometry'
    enddo
    ! Fractional translations are reduced modulo the simulation cell, not an
    ! assumed primitive cell. Conventional FCC cells require centering operations.
    do a=1,nsym
      do b=a+1,nsym
        delta=SymMatA(:,4,a)-SymMatA(:,4,b)
        if (maxval(abs(SymMatA(:,1:3,a)-SymMatA(:,1:3,b)))<1d-10.and. &
            maxval(abs(delta-anint(delta)))<1d-10) error stop 'Symmetry: duplicate operation'
      end do
      do b=1,nsym
        rotation=matmul(SymMatA(:,1:3,a),SymMatA(:,1:3,b))
        translation=matmul(SymMatA(:,1:3,a),SymMatA(:,4,b))+SymMatA(:,4,a)
        found=.false.
        do c=1,nsym
          delta=translation-SymMatA(:,4,c)
          if (maxval(abs(rotation-SymMatA(:,1:3,c)))<1d-10.and. &
              maxval(abs(delta-anint(delta)))<1d-10) then
            found=.true.
            exit
          end if
        end do
        if (.not.found) error stop 'Symmetry: operations do not form a closed group in the simulation cell'
      end do
    end do
  end subroutine symmetry_validate_group

  subroutine symmetry_validate_field(direction)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(8),intent(in) :: direction(3)
    real(8) :: rotated(3)
    integer :: isym
    if (.not.use_symmetry) return
    if (.not.flag_init.or..not.allocated(SymMatB)) error stop 'Symmetry field check before initialization'
    if (size(SymMatB,3)==0) error stop 'Symmetry: empty retained group'
    if (.not.all(ieee_is_finite(direction))) error stop 'Symmetry: nonfinite field direction'
    do isym=1,size(SymMatB,3)
      rotated=matmul(Bmat,matmul(SymMatB(:,1:3,isym),matmul(Binv,direction)))
      if (.not.all(ieee_is_finite(rotated))) error stop 'Symmetry: nonfinite field transformation'
      if (maxval(abs(rotated-direction))>1d-10*max(1d0,maxval(abs(direction)))) &
        error stop 'Symmetry: retained operation changes the RT field direction'
    end do
  end subroutine symmetry_validate_field

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
