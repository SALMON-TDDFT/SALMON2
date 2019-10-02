module cif_format_module

  use atomic_number_module, only: get_atomic_number, get_element_name

  implicit none

  private
  public :: check_cif_format
  public :: read_atom_cif

  character(30),parameter :: keyword(11)=(/ &
       "_cell_length_a            ", &
       "_cell_length_b            ", &
       "_cell_length_c            ", &
       "_cell_angle_alpha         ", &
       "_cell_angle_beta          ", &
       "_cell_angle_gamma         ", &
       "_symmetry_equiv_pos_as_xyz", &
       "_atom_site_label          ", &
       "_atom_site_fract_x        ", &
       "_atom_site_fract_y        ", &
       "_atom_site_fract_z        " /)

contains


  subroutine check_cif_format( unit, ierr )
    implicit none
    integer,intent(in)  :: unit
    integer,intent(out) :: ierr
    character(30) :: cbuf
    ierr=-1
    rewind unit
    do
       read(unit,*,END=90) cbuf
       if ( any(cbuf==keyword) ) then
          ierr = 0
          exit
       end if
    end do
90  continue
  end subroutine check_cif_format


  subroutine read_atom_cif( unit )
    implicit none
    integer,intent(in) :: unit
    logical :: flag
    integer,allocatable :: ki_atom(:),md_atom(:),zn_atom(:)
    integer,parameter :: u1=10,u2=20,u3=30
    integer :: i,j,m,n,i0,i1,i2,i3,z,nsym,itmp(0:3)
    integer :: nbas,isym,natm
    integer,allocatable :: iatm(:),katm(:)
    character(40) :: cbuf, cbuf1, cbuf2
    character(40),allocatable :: cdummy(:)
    real(8),parameter :: bohr=0.529177210903d0
    real(8) :: alat(3,3),alatl(3),angle(3),deg2rad,rr
    real(8) :: R(3,4),asi(3),rsi(3),rtm(3),Rasi(3),ax
    real(8),allocatable :: rot(:,:,:),atm(:,:)
    real(8),allocatable :: aa_atom(:,:)

    do i=1,3
       call find_key( keyword(i), unit, flag )
       if ( flag ) read(unit,*) cbuf, cbuf2
       n=index(cbuf2,"(")-1
       if ( n == -1 ) n=len_trim(cbuf2)
       if ( n > 0 ) read(cbuf2(1:n),*) alatl(i)
       write(*,'(1x,3f15.10)') alatl(i)
    end do
    do i=1,3
       call find_key( keyword(i+3), unit, flag )
       if ( flag ) read(unit,*) cbuf, angle(i)
       write(*,'(1x,3f15.10)') angle(i)
    end do

    deg2rad = acos(-1.0d0)/180.0d0

    angle(:) = angle(:)*deg2rad
    alatl(:) = alatl(:)/bohr

    alat(:,1) = (/ 1.0d0, 0.0d0, 0.0d0 /)
    alat(1,2) = cos( angle(3) )
    alat(2,2) = sin( angle(3) )
    alat(1,3) = cos( angle(2) )
    alat(2,3) = cos( angle(1) )/sin( angle(3) )
    alat(3,3) = sqrt( 1.0d0 - alat(1,3)**2 - alat(2,3)**2 )

    alat(:,1) = alat(:,1)*alatl(1)
    alat(:,2) = alat(:,2)*alatl(2)
    alat(:,3) = alat(:,3)*alatl(3)

    where( abs(alat) < 1.d-5 )
       alat=0.0d0
    end where

! ---

    nsym=0

    call find_key( keyword(7), unit, flag )

    if ( flag ) then

       read(unit,*) cbuf

       open(u1,file="cif_sym.dat")

       i=0
       do
          read(unit,'(a)') cbuf

! There are two types of cif format
!   a)  1 'x, y, z'
!   b)  'x, y, z'
! following treatment absorbs the both types
!
          read(cbuf,*,END=22) cbuf1,cbuf2
          cbuf1=cbuf2
22        continue
          cbuf=cbuf1

          if ( cbuf /= "loop_" ) then
             i=i+1
             call chr_to_matrix( cbuf, R )
             write(u1,'(1x,3f20.15,2x,f20.15)') (R(1,j),j=1,4)
             write(u1,'(1x,3f20.15,2x,f20.15)') (R(2,j),j=1,4)
             write(u1,'(1x,3f20.15,2x,f20.15)') (R(3,j),j=1,4)
             cycle
          else
             exit
          end if

       end do
       nsym=i

       close(u1)

! ---

       if ( nsym == 0 ) then
          nsym=1
          allocate( rot(3,4,nsym) ); rot=0.0d0
          rot(1,1,1)=1.0d0
          rot(2,2,1)=1.0d0
          rot(3,3,1)=1.0d0
       else if ( nsym > 0 ) then
          allocate( rot(3,4,nsym) ); rot=0.0d0
          open(u1,file="cif_sym.dat",status="old")
          do i=1,nsym
             read(u1,*) (rot(1,j,i),j=1,4)
             read(u1,*) (rot(2,j,i),j=1,4)
             read(u1,*) (rot(3,j,i),j=1,4)
          end do
          close(u1)
       end if

       itmp=0
       i=0
       do
          read(unit,*,END=92) cbuf
          i=i+1
          if ( cbuf == keyword( 8) ) itmp(0)=i
          if ( cbuf == keyword( 9) ) itmp(1)=i
          if ( cbuf == keyword(10) ) itmp(2)=i
          if ( cbuf == keyword(11) ) itmp(3)=i
          if ( cbuf(1:5) /= "_atom" ) then
             backspace(unit)
             exit
          end if
       end do
92     continue

       n=maxval(itmp)

! ---

       allocate( cdummy(n) ); cdummy=""

       open(u2,file="cif_bas.dat")

       nbas=0
       do
          read(unit,*,END=94) cdummy(1:n)

          call get_atomic_number( cdummy(itmp(0)), z )
          if ( z > 0 ) then
             do i=1,3
                cbuf=cdummy(itmp(i))
                m=index(cbuf,"(")-1
                if ( m == -1 ) m=len_trim(cbuf)
                read(cbuf(1:m),*) asi(i)
             end do
             nbas=nbas+1
             !write(*,'(1x,i3,2x,a2,2x,i3,2x,3f15.10)') &
             !     nbas,cdummy(itmp(0)),z,(asi(i),i=1,3)
             do isym=1,nsym
                Rasi(:) = matmul( rot(:,1:3,isym), asi(:) )
                Rasi(:) = Rasi(:) + rot(:,4,isym)
                write(u2,'(1x,a3,2x,i3,2x,3f20.15)') cdummy(itmp(0)),z,Rasi(:)
             end do
          end if
       end do

94     close(u2)

       deallocate( cdummy )
       deallocate( rot )

! ---

       natm = nbas*nsym

       allocate( atm(3,natm) ); atm=0.0d0
       allocate( katm(natm)  ); katm=0
       allocate( iatm(118)   ); iatm=0

       open(u2,file="cif_bas.dat",status="old")

       n=0
       loop_i : do i=1,natm
          read(u2,*) cbuf,z,asi(:)
          call shift_aa_coordinates_atom( asi )
          rsi=matmul( alat, asi )
          do j=1,n
             rtm=matmul( alat, atm(:,j) )
             rr=sum( (rsi-rtm)**2 )
             if ( rr < 1.d-3 ) cycle loop_i
          end do
          if ( iatm(z) == 0 ) iatm(z)=maxval(iatm)+1
          n=n+1
          atm(:,n)=asi(:)
          katm(n)=iatm(z)
          !write(*,'(1x,2(i3,2x),3f10.5)') n,katm(n),(asi(j),j=1,3)
       end do loop_i

       close(u2)

! ---

       allocate( aa_atom(3,n) ); aa_atom=0.0d0
       allocate( ki_atom(n)   ); ki_atom=0
       allocate( md_atom(n)   ); md_atom=1

       aa_atom(:,:) = atm(:,1:n)
       ki_atom(:) = katm(1:n)

       m=maxval( iatm )
       allocate( zn_atom(m) ); zn_atom=1

       do z=1,size(iatm)
          i=iatm(z)
          if ( i > 0 ) zn_atom(i)=z
       end do

       deallocate( iatm )
       deallocate( katm )
       deallocate( atm  )

       open(u2,file="cif_bas_red.dat")
       write(u2,'("&atomic_red_coor")')

       open(u3,file="cif_bas_xyz.dat")
       write(u3,'("&atomic_coor")')

       do i=1,n
          call get_element_name( zn_atom(ki_atom(i)), cbuf(1:2) )
          write(u2,'(2x,a2,3f12.6,i3)') cbuf(1:2), aa_atom(:,i), ki_atom(i)
          rsi=matmul( alat, aa_atom(:,i) )
          write(u3,'(2x,a2,3f12.6,i3)') cbuf(1:2), rsi(:), ki_atom(i)
       end do

       write(u3,'("/")')
       close(u3)

       write(u2,'("/")')
       close(u2)


    end if ! flag

    deallocate( zn_atom )
    deallocate( md_atom )
    deallocate( ki_atom )
    deallocate( aa_atom )

  end subroutine read_atom_cif


  subroutine find_key( key, unit, flag )
    implicit none
    character(*),intent(in) :: key
    integer,intent(in)  :: unit
    logical,intent(out) :: flag
    character(30) :: cbuf
    flag=.false.
    rewind unit
    do
       read(unit,*,END=90) cbuf
       if ( cbuf == key ) then
          backspace(unit)
          flag=.true.
          return
       else
          cycle
       end if
    end do
90  continue
  end subroutine find_key


  subroutine chr_to_matrix( cbuf, R )
    implicit none
    character(*),intent(in) :: cbuf
    real(8),intent(out) :: R(3,4)
    character(10) :: str(3)
    integer :: i,j
    real(8) :: f
    R=0.0d0
    j=1
    f=1.0d0
    do i=1,len_trim(cbuf)
       if ( cbuf(i:i) == "-" ) then
          f=-1.0d0
       end if
       if ( cbuf(i:i) == "x" ) then
          R(j,1) = f
          f=1.0d0
       end if
       if ( cbuf(i:i) == "y" ) then
          R(j,2) = f
          f=1.0d0
       end if
       if ( cbuf(i:i) == "z" ) then
          R(j,3) = f
          f=1.0d0
       end if
       if ( 49 <= iachar(cbuf(i:i)) .and. iachar(cbuf(i:i)) <= 57 ) then
          if ( R(j,4) == 0.0d0 ) then
             R(j,4) = f*(iachar(cbuf(i:i))-48)
          else
             R(j,4) = R(j,4)/dble( iachar(cbuf(i:i))-48 )
          end if
          f=1.0d0
       end if
       if ( cbuf(i:i) == "," ) then
          j=j+1
          f=1.0d0
       end if
    end do
    !do i=1,3
    !   write(*,'(1x,3f15.10,3x,f15.10)') (R(i,j),j=1,4)
    !end do
  end subroutine chr_to_matrix


  subroutine shift_aa_coordinates_atom( aa_atom )
    implicit none
    real(8),intent(inout) :: aa_atom(:)
    integer :: j
    do j=1,3
       do
          if ( aa_atom(j) < 0.0d0 ) then
             aa_atom(j)=aa_atom(j)+1.0d0
             cycle
          else if ( aa_atom(j) >= 1.0d0 ) then
             aa_atom(j)=aa_atom(j)-1.0d0
             cycle
          else
             exit
          end if
       end do
    end do ! j
  end subroutine shift_aa_coordinates_atom


end module cif_format_module
