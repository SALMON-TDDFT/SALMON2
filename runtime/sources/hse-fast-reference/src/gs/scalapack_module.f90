!
!  Copyright 2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
module scalapack_module
  implicit none

  public :: create_gridmap
  public :: get_blocking_factor
  public :: init_blacs

private

contains
  subroutine create_gridmap(info)
    use structures, only: s_parallel_info
    implicit none
    type(s_parallel_info),intent(inout) :: info

    integer,parameter :: iclose_comm = 3 ! 1: rko, 2: r, 3: o
    integer :: ii,ix,iy,iz
    integer :: k1,k2
    integer :: px,py,pz,po,pk

    select case(iclose_comm)
    ! close on icomm_rko
    case(1)
      info%icomm_sl = info%icomm_rko

      info%npcol = int(sqrt(dble(info%isize_rko)))
      info%npcol = info%npcol + mod(info%npcol,2)
      do ii=1,100
        info%nprow=info%isize_rko/info%npcol
        if(info%nprow*info%npcol == info%isize_rko) exit
        info%npcol=info%npcol+2
      end do
      if (info%nprow > info%npcol) then
        k1          = info%nprow
        info%nprow = info%npcol
        info%npcol = k1
      end if
      if (info%nprow*info%npcol /= info%isize_rko) &
        stop 'scalapack_module: fatal error, please check nprow and npcol'

      allocate( info%gridmap(info%nprow,info%npcol) )
      info%gridmap = - 1
      ii = 0
      do k2=1,info%npcol
      do k1=1,info%nprow
        info%gridmap(k1,k2) = ii
        ii = ii + 1
      end do
      end do
      if (ii /= info%isize_rko) &
        stop 'scalapack_module: fatal error, please check gridmap'

    ! close on icomm_r
    case(2)
      info%icomm_sl = info%icomm_r

      info%npcol = int(sqrt(dble(info%isize_r)))
      info%npcol = info%npcol + mod(info%npcol,2)
      do ii=1,100
        info%nprow=info%isize_r/info%npcol
        if(info%nprow*info%npcol == info%isize_r) exit
        info%npcol=info%npcol+2
      end do
      if (info%nprow > info%npcol) then
        k1          = info%nprow
        info%nprow = info%npcol
        info%npcol = k1
      end if
      if (info%nprow*info%npcol /= info%isize_r) &
        stop 'scalapack_module: fatal error, please check nprow and npcol'

      allocate( info%gridmap(info%nprow,info%npcol) )
      info%gridmap = -1
      pk = info%iaddress(5)
      po = info%iaddress(4)
      ii = 0
      do iz=0,info%nprgrid(3)-1
      do iy=0,info%nprgrid(2)-1
      do ix=0,info%nprgrid(1)-1
        ii = ii + 1
        k1 = mod(ii-1,info%nprow)
        k2 = (ii-1)/info%nprow
        info%gridmap(k1+1,k2+1) = info%imap(ix,iy,iz,po,pk)
      end do
      end do
      end do
      if (ii /= info%isize_r) &
        stop 'scalapack_module: fatal error, please check gridmap'

    ! close on icomm_o
    case(3)
      info%icomm_sl = info%icomm_o

      info%npcol = int(sqrt(dble(info%isize_o)))
      info%npcol = info%npcol + mod(info%npcol,2)
      do ii=1,100
        info%nprow=info%isize_o/info%npcol
        if(info%nprow*info%npcol == info%isize_o) exit
        info%npcol=info%npcol+2
      end do
      if (info%nprow > info%npcol) then
        k1          = info%nprow
        info%nprow = info%npcol
        info%npcol = k1
      end if
      if (info%nprow*info%npcol /= info%isize_o) &
        stop 'scalapack_module: fatal error, please check nprow and npcol'

      allocate( info%gridmap(info%nprow,info%npcol) )
      info%gridmap = - 1
      pk = info%iaddress(5)
      pz = info%iaddress(3)
      py = info%iaddress(2)
      px = info%iaddress(1)
      ii = 0
      do k2=1,info%npcol
      do k1=1,info%nprow
        info%gridmap(k1,k2) = info%imap(px,py,pz,ii,pk)
        ii = ii + 1
      end do
      end do
      if (ii /= info%isize_o) &
        stop 'scalapack_module: fatal error, please check gridmap'
    end select
  end subroutine create_gridmap

  subroutine get_blocking_factor(info,n,mb,nb)
    use structures, only: s_parallel_info
    implicit none
    type(s_parallel_info),intent(in) :: info
    integer,intent(in)               :: n
    integer,intent(inout)            :: mb,nb

    if (.not. allocated(info%gridmap)) &
      stop 'scalapack_module: gridmap not constructed.'

    ! cyclic distribution
    mb = 1
    nb = 1
  end subroutine get_blocking_factor

  subroutine init_blacs(info,m)
    use structures, only: s_parallel_info
    use communication, only: comm_summation
    implicit none
    integer :: NUMROC
    type(s_parallel_info),intent(inout) :: info
    integer,intent(in) :: m

    integer :: n,mb,nb
    integer :: len_iwork
    integer :: ictxt,ierr

    integer :: npo, i, j, i_loc, j_loc, proc_row, proc_col, ip
    integer,allocatable :: icount(:)

    if (info%flag_blacs_gridinit) return

    if (.not. allocated(info%gridmap)) &
      stop 'scalapack_module: gridmap not constructed.'

    n = m
    call get_blocking_factor(info,n,mb,nb)

    call BLACS_PINFO( info%iam, info%nprocs )
    if (info%nprocs < 1) then
      !info%nprocs = info%isize_rko
      call BLACS_SETUP( info%iam, info%nprocs )
    end if

    call BLACS_GET( 0, 0, ictxt )
    call BLACS_GRIDMAP( ictxt, info%gridmap, info%nprow, info%nprow, info%npcol )
    call BLACS_GRIDINFO( ictxt, info%nprow, info%npcol, info%myrow, info%mycol )
    info%nrow_local = NUMROC( n, mb, info%myrow, 0, info%nprow )
    info%ncol_local = NUMROC( n, nb, info%mycol, 0, info%npcol )
    info%lda        = max(1, info%nrow_local)

    call DESCINIT( info%desca, n, n, mb, nb, 0, 0, ictxt, info%lda, ierr )
    call DESCINIT( info%descz, n, n, mb, nb, 0, 0, ictxt, info%lda, ierr )

    ! determine the working memory size of PDSYEVD and PZHEEVD
    len_iwork = 2 + 7*n + 8*info%npcol
    call set_from_pdsyevd
    call set_from_pzheevd

    info%flag_blacs_gridinit = .true.

!-----------

    npo = info%nporbital
    allocate( info%ndiv(0:npo-1), icount(0:npo-1) )

    icount = 0
    do i = 1, n
    do j = 1, n
      call INFOG2L( i, j, info%desca, info%nprow, info%npcol, info%myrow, info%mycol, i_loc, j_loc, proc_row, proc_col )
      if (info%myrow == proc_row .and. info%mycol == proc_col) then
        ip = info%irank_io(j)
        icount( ip ) = icount( ip ) + 1
      end if
    end do
    end do
    info%ndiv = icount

    allocate( info%i_tbl( maxval(info%ndiv), 0:npo-1), &
              info%j_tbl( maxval(info%ndiv), 0:npo-1), &
              info%iloc_tbl( maxval(info%ndiv), 0:npo-1), &
              info%jloc_tbl( maxval(info%ndiv), 0:npo-1) )

    icount(:) = 0
    info%i_tbl = 0
    info%j_tbl = 0
    info%iloc_tbl = 0
    info%jloc_tbl = 0

    do i = 1, n
    do j = 1, n
      call INFOG2L( i, j, info%desca, info%nprow, info%npcol, info%myrow, info%mycol, i_loc, j_loc, proc_row, proc_col )
      if (info%myrow == proc_row .and. info%mycol == proc_col) then
        ip = info%irank_io(j)
        icount( ip ) = icount( ip ) + 1
        info%i_tbl( icount(ip), ip ) = i
        info%j_tbl( icount(ip), ip ) = j
        info%iloc_tbl( icount(ip), ip ) = i_loc
        info%jloc_tbl( icount(ip), ip ) = j_loc
      end if
    end do
    end do

    deallocate(icount)

    return

  contains
    subroutine set_from_pdsyevd
      implicit none
      real(8), allocatable :: h_div(:,:), v_div(:,:), e(:)
      real(8) :: rtmp(1)
      integer :: itmp(1)
      integer :: len_work0,trilwmin

      allocate( h_div(info%nrow_local,info%ncol_local), &
                v_div(info%nrow_local,info%ncol_local), e(n) )

      trilwmin = 3*n + max( mb*(info%nrow_local+1), 3*mb )
      len_work0 = max( 1+6*n+2*info%nrow_local*info%ncol_local, trilwmin )

      call PDSYEVD( 'V', 'L', n, h_div, 1, 1, info%desca, e, v_div, 1, 1, info%descz, &
                    rtmp, -1, itmp, len_iwork, ierr )

      info%len_work = max( nint(rtmp(1))*10, len_work0*10 )

      deallocate( e, v_div, h_div )
    end subroutine

    subroutine set_from_pzheevd
      implicit none
      complex(8), allocatable :: h_div(:,:), v_div(:,:)
      real(8), allocatable    :: e(:)
      real(8)    :: rtmp(1)
      complex(8) :: ctmp(1)
      integer    :: itmp(1)

      allocate( h_div(info%nrow_local,info%ncol_local), &
                v_div(info%nrow_local,info%ncol_local), e(n) )

      call PZHEEVD( 'V', 'L', n, h_div, 1, 1, info%desca, e, v_div, 1, 1, info%descz, &
                    ctmp, -1, rtmp, -1, itmp, len_iwork, ierr )

      info%len_work = max( info%len_work, n+(info%nrow_local+info%ncol_local+mb)*mb )
      info%len_work = max( info%len_work, nint(real(ctmp(1))) )

      info%len_rwork = max( nint(rtmp(1)), (1+8*n+2*info%nrow_local*info%ncol_local)*2 )

      deallocate( e, v_div, h_div )
    end subroutine
  end subroutine init_blacs

end module scalapack_module
