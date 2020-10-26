module sym_rho_sub

  use sym_sub, only: SymMatA, SymMatB, use_symmetry
  use communication

  implicit none

  private
  public :: init_sym_rho
  public :: sym_rho

  integer :: is(3),ie(3),num(3)
  integer :: icomm_r
  integer :: ng
  integer,allocatable :: nl(:)
  integer,allocatable :: grid_list(:,:,:)
  integer :: nsym

contains

  subroutine init_sym_rho( num_, is_, ie_, icomm_r_ )
    implicit none
    integer,intent(in) :: num_(3),is_(3),ie_(3),icomm_r_
    integer :: i1,i2,i3,isym,j1,j2,j3
    integer,allocatable :: icheck(:,:,:),itmp(:),icheck1(:,:,:)
    real(8) :: r1,r2,r3,t1,t2,t3,d1,d2,d3

    if ( .not.use_symmetry ) return

    num=num_
    is=is_-1
    ie=ie_-1
    icomm_r=icomm_r_
    nsym=size(SymMatA,3)

    allocate( icheck(0:num(1)-1,0:num(2)-1,0:num(3)-1) ); icheck=0
    allocate( icheck1(0:num(1)-1,0:num(2)-1,0:num(3)-1) ); icheck1=0
    do i3=0,num(3)-1
    do i2=0,num(2)-1
    do i1=0,num(1)-1
       do isym=1,nsym
          t1=i1-num(1)*SymMatA(1,4,isym)
          t2=i2-num(2)*SymMatA(2,4,isym)
          t3=i3-num(3)*SymMatA(3,4,isym)
          r1=t1*SymMatB(1,1,isym)+t2*SymMatB(2,1,isym)+t3*SymMatB(3,1,isym)
          r2=t1*SymMatB(1,2,isym)+t2*SymMatB(2,2,isym)+t3*SymMatB(3,2,isym)
          r3=t1*SymMatB(1,3,isym)+t2*SymMatB(2,3,isym)+t3*SymMatB(3,3,isym)
          if ( r1 < 0.0d0 ) r1=r1+num(1)
          if ( r2 < 0.0d0 ) r2=r2+num(2)
          if ( r3 < 0.0d0 ) r3=r3+num(3)
          if ( r1 > dble(num(1)) .or. r2 > dble(num(2)) .or. r3 > dble(num(3)) ) then
             write(*,*) "xxxxx ",r1,r2,r3
             stop "xxxxx"
          end if
!          write(*,'(1x,i4,2x,3i4,2x,3f8.3)') isym,i1,i2,i3,r1,r2,r3
          j1=nint(r1); d1=abs(r1-j1)
          j2=nint(r2); d2=abs(r2-j2)
          j3=nint(r3); d3=abs(r3-j3)
          if ( d1 > 1.d-3 .or. d2 > 1.d-3 .or. d3 > 1.d-3 ) then
             write(*,*) "yyyyy",r1,r2,r3
             stop "yyyyy"
          end if
          if ( icheck(j1,j2,j3) == 0 ) then
             icheck(j1,j2,j3) = 1 + i1 + i2*num(1) + i3*num(1)*num(2)
          end if
          icheck1(j1,j2,j3)=icheck1(j1,j2,j3)+1
       end do
    end do
    end do
    end do
    allocate( itmp(size(icheck)) ); itmp=0

    ng=0
    do i3=0,num(3)-1
    do i2=0,num(2)-1
    do i1=0,num(1)-1
       j1=icheck(i1,i2,i3)
       if ( ng==0 ) then
          ng=ng+1
          itmp(ng)=j1
       else if ( .not.any(itmp(1:ng)==j1) ) then
          ng=ng+1
          itmp(ng)=j1
       end if         
    end do
    end do
    end do
    !write(*,*) ng,itmp(1:ng)
    !write(*,*) count(icheck/=0),minval(icheck),maxval(icheck)

    allocate( nl(ng) ); nl=0

    do i1=1,ng
       nl(i1) = count( icheck == itmp(i1) )
    end do

    allocate( grid_list(4,maxval(nl),ng) ); grid_list=0

    nl=0
    do i3=0,num(3)-1
    do i2=0,num(2)-1
    do i1=0,num(1)-1
       j1=icheck(i1,i2,i3)
       do j2=1,ng
          if ( itmp(j2) == j1 ) exit
       end do
       if ( j2 > ng ) stop "zzzzz"
       nl(j2)=nl(j2)+1
       grid_list(1,nl(j2),j2)=i1
       grid_list(2,nl(j2),j2)=i2
       grid_list(3,nl(j2),j2)=i3
       grid_list(4,nl(j2),j2)=icheck1(i1,i2,i3)
    end do
    end do
    end do

    do i1=1,ng
    do i2=1,nl(i1)
       j1=grid_list(1,i2,i1)
       j2=grid_list(2,i2,i1)
       j3=grid_list(3,i2,i1)
       i3=grid_list(4,i2,i1)
!       write(*,'(1x,7i4)') i1,i2,j1,j2,j3,i3,icheck(j1,j2,j3)
    end do
    end do

    deallocate( itmp )
    deallocate( icheck1 )
    deallocate( icheck )

  end subroutine init_sym_rho

#ifdef test
  subroutine sym_rho( rho )
    implicit none
    real(8),intent(inout) :: rho(:,:,:)
    integer :: ig,il,i1,i2,i3
    real(8),allocatable :: sbuf(:), rbuf(:)
    if ( .not.use_symmetry ) return
    allocate( sbuf(ng) ); sbuf=0.0d0
    allocate( rbuf(ng) ); rbuf=0.0d0
    do ig=1,ng
       do il=1,nl(ig)
          i1=grid_list(1,il,ig)+1
          i2=grid_list(2,il,ig)+1
          i3=grid_list(3,il,ig)+1
          if ( is(1) <= i1 .and. i1 <= ie(1) .and. &
               is(2) <= i2 .and. i2 <= ie(2) .and. &
               is(3) <= i3 .and. i3 <= ie(3) ) then
             sbuf(ig)=sbuf(ig)+rho(i1-is(1)+1,i2-is(2)+1,i3-is(3)+1)
          end if
       end do !il
    end do ! ig

    call comm_summation( sbuf, rbuf, ng, icomm_r )

    !rbuf=rbuf/dble(nl)
    rbuf=rbuf/dble(nsym)
    do ig=1,ng
       do il=1,nl(ig)
          i1=grid_list(1,il,ig)+1
          i2=grid_list(2,il,ig)+1
          i3=grid_list(3,il,ig)+1
          if ( is(1) <= i1 .and. i1 <= ie(1) .and. &
               is(2) <= i2 .and. i2 <= ie(2) .and. &
               is(3) <= i3 .and. i3 <= ie(3) ) then
             rho(i1-is(1)+1,i2-is(2)+1,i3-is(3)+1)=rbuf(ig)
          end if
       end do !il
    end do ! ig
    deallocate( rbuf )
    deallocate( sbuf )
write(*,*) "sym_rho(1)"
  end subroutine sym_rho

#else

  subroutine sym_rho( rho )
    implicit none
    real(8),intent(inout) :: rho(:,:,:)
    integer :: i1,i2,i3,j1,j2,j3,isym,is1,is2,is3
    real(8) :: t1,t2,t3,r1,r2,r3,d1,d2,d3,fac
    real(8),allocatable :: work(:,:,:),wtmp(:,:,:)
    if ( .not.use_symmetry ) return
!
    allocate( work(0:num(1)-1,0:num(2)-1,0:num(3)-1) ); work=0.0d0
    allocate( wtmp(0:num(1)-1,0:num(2)-1,0:num(3)-1) ); wtmp=0.0d0
    wtmp(is(1):ie(1),is(2):ie(2),is(3):ie(3))=rho(:,:,:)
    call comm_summation( wtmp, work, size(work), icomm_r )
    rho=0.0d0
!
    is1=is(1)-1
    is2=is(2)-1
    is3=is(3)-1
    fac=1.0d0/nsym
!$omp parallel do collapse(3) private(i1,i2,i3,isym)  &
!$omp             private(t1,t2,t3,r1,r2,r3,j1,j2,j3,d1,d2,d3) &
!$omp             reduction(+:rho)
    do i3=is(3),ie(3)
    do i2=is(2),ie(2)
    do i1=is(1),ie(1)
       do isym=1,nsym
          t1=i1-num(1)*SymMatA(1,4,isym)
          t2=i2-num(2)*SymMatA(2,4,isym)
          t3=i3-num(3)*SymMatA(3,4,isym)
          r1=t1*SymMatB(1,1,isym)+t2*SymMatB(2,1,isym)+t3*SymMatB(3,1,isym)
          r2=t1*SymMatB(1,2,isym)+t2*SymMatB(2,2,isym)+t3*SymMatB(3,2,isym)
          r3=t1*SymMatB(1,3,isym)+t2*SymMatB(2,3,isym)+t3*SymMatB(3,3,isym)
          if ( r1 < 0.0d0 ) r1=r1+num(1)
          if ( r2 < 0.0d0 ) r2=r2+num(2)
          if ( r3 < 0.0d0 ) r3=r3+num(3)
          if ( r1 >= dble(num(1)) .or. r2 >= dble(num(2)) .or. r3 >= dble(num(3)) ) then
             write(*,*) "xxxxx ",r1,r2,r3
             stop "xxxxx"
          end if
          j1=nint(r1); d1=abs(r1-j1)
          j2=nint(r2); d2=abs(r2-j2)
          j3=nint(r3); d3=abs(r3-j3)
          if ( d1 > 1.d-3 .or. d2 > 1.d-3 .or. d3 > 1.d-3 ) then
             write(*,*) "yyyyy",r1,r2,r3
             stop "yyyyy"
          end if
          rho(i1-is1,i2-is2,i3-is3)=rho(i1-is1,i2-is2,i3-is3)+fac*work(j1,j2,j3)
       end do
    end do
    end do
    end do
!$omp end parallel do
    deallocate( work )
  end subroutine sym_rho
#endif

end module sym_rho_sub
