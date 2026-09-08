module sym_sub

  use communication, only: comm_get_globalinfo, comm_is_root

  implicit none

  private
  public :: read_sw_symmetry
  public :: read_symmetry_file
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
    real(8), allocatable :: matrices(:,:,:)
    character(512) :: message

    inquire( FILE=sym_file, EXIST=flag )
    if ( .not.flag ) then
       if ( DISPLAY ) write(*,*) "symmetry-operation file ( "//sym_file//" ) can not be found."
       stop 'stop@read_SymMat'
    else
       if ( DISPLAY ) write(*,*) "symmetry-operation file is found ( "//sym_file//" )."
    end if

    call read_symmetry_file(sym_file, matrices, flag, message)
    if (.not. flag) then
       if (DISPLAY) write(*,*) trim(message)
       stop 'stop@read_SymMat'
    end if
    call move_alloc(matrices, SymMatR)

  end subroutine read_SymMat


  subroutine read_symmetry_file(filename, matrices, ok, message)
    implicit none
    character(*), intent(in) :: filename
    real(8), allocatable, intent(out) :: matrices(:,:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: unit, ios, nrow, irow, isym, row
    character(1024) :: line
    character(512) :: iomsg, row_message
    real(8) :: values(4)
    logical :: active, row_ok

    ok = .false.
    message = ''
    open(newunit=unit, file=trim(filename), status='old', action='read', &
         iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
       write(message,'(a,a,a,a)') 'cannot open symmetry file "', &
            trim(filename), '": ', trim(iomsg)
       return
    end if

    nrow = 0
    do
       read(unit,'(a)',iostat=ios,iomsg=iomsg) line
       if (ios < 0) exit
       if (ios > 0) then
          write(message,'(a,i0,a,a)') 'failed reading symmetry row ', &
               nrow + 1, ': ', trim(iomsg)
          close(unit)
          return
       end if
       call parse_symmetry_row(line, values, active, row_ok, row_message)
       if (.not. row_ok) then
          write(message,'(a,i0,a,a)') 'malformed symmetry row ', nrow + 1, &
               ': ', trim(row_message)
          close(unit)
          return
       end if
       if (.not. active) cycle
       nrow = nrow + 1
    end do

    if (nrow == 0) then
       message = 'symmetry file contains no operations'
       close(unit)
       return
    end if
    if (mod(nrow,3) /= 0) then
       write(message,'(a,i0,a)') 'incomplete symmetry operation: ', nrow, &
            ' nonblank rows is not a multiple of 3'
       close(unit)
       return
    end if

    allocate(matrices(3,4,nrow/3))
    matrices = 0.0d0
    rewind(unit)
    irow = 0
    do
       read(unit,'(a)',iostat=ios,iomsg=iomsg) line
       if (ios < 0) exit
       if (ios > 0) then
          write(message,'(a,i0,a,a)') 'failed rereading symmetry row ', &
               irow + 1, ': ', trim(iomsg)
          deallocate(matrices)
          close(unit)
          return
       end if
       call parse_symmetry_row(line, values, active, row_ok, row_message)
       if (.not. row_ok) then
          write(message,'(a,i0,a,a)') 'malformed symmetry row ', irow + 1, &
               ': ', trim(row_message)
          deallocate(matrices)
          close(unit)
          return
       end if
       if (.not. active) cycle
       irow = irow + 1
       isym = (irow - 1)/3 + 1
       row = mod(irow - 1,3) + 1
       matrices(row,:,isym) = values
    end do
    close(unit)
    ok = .true.
  end subroutine read_symmetry_file


  subroutine parse_symmetry_row(line, values, active, ok, message)
    implicit none
    character(*), intent(in) :: line
    real(8), intent(out) :: values(4)
    logical, intent(out) :: active, ok
    character(*), intent(out) :: message
    character(len(line)) :: data
    character(512) :: iomsg
    real(8) :: extra
    integer :: comment, hash_comment, bang_comment, ios

    data = line
    hash_comment = index(data, '#')
    bang_comment = index(data, '!')
    if (hash_comment == 0) then
       comment = bang_comment
    else if (bang_comment == 0) then
       comment = hash_comment
    else
       comment = min(hash_comment, bang_comment)
    end if
    if (comment == 1) then
       data = ''
    else if (comment > 1) then
       data = data(:comment-1)
    end if

    active = len_trim(data) > 0
    ok = .true.
    message = ''
    values = 0.0d0
    if (.not. active) return

    read(data,*,iostat=ios,iomsg=iomsg) values
    if (ios /= 0) then
       ok = .false.
       message = trim(iomsg)
       return
    end if
    read(data,*,iostat=ios,iomsg=iomsg) values, extra
    if (ios >= 0) then
       ok = .false.
       message = 'extra field after four matrix values'
    end if
  end subroutine parse_symmetry_row


end module sym_sub
