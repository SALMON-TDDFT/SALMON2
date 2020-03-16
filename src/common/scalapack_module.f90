!
!  Copyright 2019 SALMON developers
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
  public :: get_blocking_factor, get_appropriate_matsize
  public :: init_blacs

private

contains
  subroutine create_gridmap(pinfo,info)
    use structures, only: s_process_info, s_parallel_info
    implicit none
    type(s_process_info),intent(inout)  :: pinfo
    type(s_parallel_info),intent(in) :: info

    integer,parameter :: iclose_comm = 3 ! 1: rko, 2: r, 3: o
    integer :: ii,ix,iy,iz
    integer :: k1,k2
    integer :: px,py,pz,po,pk

    select case(iclose_comm)
    ! close on icomm_rko
    case(1)
      pinfo%icomm_sl = info%icomm_rko

      pinfo%npcol = int(sqrt(dble(info%isize_rko)))
      pinfo%npcol = pinfo%npcol + mod(pinfo%npcol,2)
      do ii=1,100
        pinfo%nprow=info%isize_rko/pinfo%npcol
        if(pinfo%nprow*pinfo%npcol == info%isize_rko) exit
        pinfo%npcol=pinfo%npcol+2
      end do
      if (pinfo%nprow > pinfo%npcol) then
        k1          = pinfo%nprow
        pinfo%nprow = pinfo%npcol
        pinfo%npcol = k1
      end if
      if (pinfo%nprow*pinfo%npcol /= info%isize_rko) &
        stop 'scalapack_module: fatal error, please check nprow and npcol'

      allocate( pinfo%gridmap(pinfo%nprow,pinfo%npcol) )
      pinfo%gridmap = - 1
      ii = 0
      do k2=1,pinfo%npcol
      do k1=1,pinfo%nprow
        pinfo%gridmap(k1,k2) = ii
        ii = ii + 1
      end do
      end do
      if (ii /= info%isize_rko) &
        stop 'scalapack_module: fatal error, please check gridmap'

    ! close on icomm_r
    case(2)
      pinfo%icomm_sl = info%icomm_r

      pinfo%npcol = int(sqrt(dble(info%isize_r)))
      pinfo%npcol = pinfo%npcol + mod(pinfo%npcol,2)
      do ii=1,100
        pinfo%nprow=info%isize_r/pinfo%npcol
        if(pinfo%nprow*pinfo%npcol == info%isize_r) exit
        pinfo%npcol=pinfo%npcol+2
      end do
      if (pinfo%nprow > pinfo%npcol) then
        k1          = pinfo%nprow
        pinfo%nprow = pinfo%npcol
        pinfo%npcol = k1
      end if
      if (pinfo%nprow*pinfo%npcol /= info%isize_r) &
        stop 'scalapack_module: fatal error, please check nprow and npcol'

      allocate( pinfo%gridmap(pinfo%nprow,pinfo%npcol) )
      pinfo%gridmap = -1
      pk = info%iaddress(5)
      po = info%iaddress(4)
      ii = 0
      do iz=0,pinfo%nprgrid(3)-1
      do iy=0,pinfo%nprgrid(2)-1
      do ix=0,pinfo%nprgrid(1)-1
        ii = ii + 1
        k1 = mod(ii-1,pinfo%nprow)
        k2 = (ii-1)/pinfo%nprow
        pinfo%gridmap(k1+1,k2+1) = info%imap(ix,iy,iz,po,pk)
      end do
      end do
      end do
      if (ii /= info%isize_r) &
        stop 'scalapack_module: fatal error, please check gridmap'

    ! close on icomm_o
    case(3)
      pinfo%icomm_sl = info%icomm_o

      pinfo%npcol = int(sqrt(dble(info%isize_o)))
      pinfo%npcol = pinfo%npcol + mod(pinfo%npcol,2)
      do ii=1,100
        pinfo%nprow=info%isize_o/pinfo%npcol
        if(pinfo%nprow*pinfo%npcol == info%isize_o) exit
        pinfo%npcol=pinfo%npcol+2
      end do
      if (pinfo%nprow > pinfo%npcol) then
        k1          = pinfo%nprow
        pinfo%nprow = pinfo%npcol
        pinfo%npcol = k1
      end if
      if (pinfo%nprow*pinfo%npcol /= info%isize_o) &
        stop 'scalapack_module: fatal error, please check nprow and npcol'

      allocate( pinfo%gridmap(pinfo%nprow,pinfo%npcol) )
      pinfo%gridmap = - 1
      pk = info%iaddress(5)
      pz = info%iaddress(3)
      py = info%iaddress(2)
      px = info%iaddress(1)
      ii = 0
      do k2=1,pinfo%npcol
      do k1=1,pinfo%nprow
        pinfo%gridmap(k1,k2) = info%imap(px,py,pz,ii,pk)
        ii = ii + 1
      end do
      end do
      if (ii /= info%isize_o) &
        stop 'scalapack_module: fatal error, please check gridmap'
    end select
  end subroutine create_gridmap

  subroutine get_blocking_factor(pinfo,n,mb,nb)
    use structures, only: s_process_info
    implicit none
    type(s_process_info),intent(in) :: pinfo
    integer,intent(in)              :: n
    integer,intent(inout)           :: mb,nb
    integer :: k1,k2

    if (.not. allocated(pinfo%gridmap)) &
      stop 'scalapack_module: gridmap not constructed.'

    k1 = (n+pinfo%nprow-1)/pinfo%nprow
    k2 = (n+pinfo%npcol-1)/pinfo%npcol
    mb = min(k1,k2)
    nb = mb
  end subroutine get_blocking_factor

  subroutine get_appropriate_matsize(pinfo,info,n)
    use structures, only: s_process_info, s_parallel_info
    use communication, only: comm_is_root
    implicit none
    type(s_process_info),intent(in)     :: pinfo
    type(s_parallel_info),intent(in) :: info
    integer,intent(inout) :: n
    integer :: mb,nb  ! blocking factor
    integer :: k

    if (.not. allocated(pinfo%gridmap)) &
      stop 'scalapack_module: gridmap not constructed.'

    call get_blocking_factor(pinfo,n,mb,nb)

    if (nb * pinfo%npcol /= n) then
      k = max(nb*pinfo%npcol, n)
      k = min(k, (nb+1)*pinfo%npcol)
      if (comm_is_root(info%id_rko)) then
        print '(A)',       '[WARNING] scalapack_module'
        print '(A,2I6)' ,  '  nprow,npcol = ',pinfo%nprow,pinfo%npcol
        print '(A,2I6)' ,  '  mb,nb       = ',mb,nb
        print '(2(A,I6))', '  nb*npcol = ',nb*pinfo%npcol,' /= ',n
        print '(A,I6,A)',  '  appropriated value is ',k,' replaced it.'
      end if
      n = k
    end if
  end subroutine get_appropriate_matsize

  subroutine init_blacs(pinfo,info,m)
    use structures, only: s_process_info, s_parallel_info
    use communication, only: comm_summation
    implicit none
    integer :: NUMROC

    type(s_process_info),intent(inout)  :: pinfo
    type(s_parallel_info),intent(in) :: info
    integer,intent(in) :: m

    integer :: n,mb,nb
    integer :: len_iwork
    integer :: ictxt,ierr

    integer :: npo, i, j, i_loc, j_loc, proc_row, proc_col, ip
    integer,allocatable :: icount(:)

    if (pinfo%flag_blacs_gridinit) return

    if (.not. allocated(pinfo%gridmap)) &
      stop 'scalapack_module: gridmap not constructed.'

    n = m
    call get_appropriate_matsize(pinfo,info,n)
    call get_blocking_factor(pinfo,n,mb,nb)

    call BLACS_PINFO( pinfo%iam, pinfo%nprocs )
    if (pinfo%nprocs < 1) then
      !pinfo%nprocs = info%isize_rko
      call BLACS_SETUP( pinfo%iam, pinfo%nprocs )
    end if

    call BLACS_GET( 0, 0, ictxt )
    call BLACS_GRIDMAP( ictxt, pinfo%gridmap, pinfo%nprow, pinfo%nprow, pinfo%npcol )
    call BLACS_GRIDINFO( ictxt, pinfo%nprow, pinfo%npcol, pinfo%myrow, pinfo%mycol )
    pinfo%nrow_local = NUMROC( n, mb, pinfo%myrow, 0, pinfo%nprow )
    pinfo%ncol_local = NUMROC( n, nb, pinfo%mycol, 0, pinfo%npcol )
    pinfo%lda        = max(1, pinfo%nrow_local)

    call DESCINIT( pinfo%desca, n, n, mb, nb, 0, 0, ictxt, pinfo%lda, ierr )
    call DESCINIT( pinfo%descz, n, n, mb, nb, 0, 0, ictxt, pinfo%lda, ierr )

    ! determine the working memory size of PDSYEVD and PZHEEVD
    len_iwork = 2 + 7*n + 8*pinfo%npcol
    call set_from_pdsyevd
    call set_from_pzheevd

    pinfo%flag_blacs_gridinit = .true.

!-----------

    npo = pinfo%nporbital
    allocate( pinfo%ndiv(0:npo-1), icount(0:npo-1) )

    icount = 0
    do i = 1, n
    do j = 1, n
      call INFOG2L( i, j, pinfo%desca, pinfo%nprow, pinfo%npcol, pinfo%myrow, pinfo%mycol, i_loc, j_loc, proc_row, proc_col )
      if (pinfo%myrow == proc_row .and. pinfo%mycol == proc_col) then
        ip = info%irank_io(j)
        icount( ip ) = icount( ip ) + 1
      end if
    end do
    end do
    pinfo%ndiv = icount

    allocate( pinfo%i_tbl( max(1,pinfo%ndiv(info%id_o)), 0:npo-1), &
              pinfo%j_tbl( max(1,pinfo%ndiv(info%id_o)), 0:npo-1), &
              pinfo%iloc_tbl( max(1,pinfo%ndiv(info%id_o)), 0:npo-1), &
              pinfo%jloc_tbl( max(1,pinfo%ndiv(info%id_o)), 0:npo-1) )

    icount(:) = 0
    pinfo%i_tbl = 0
    pinfo%j_tbl = 0
    pinfo%iloc_tbl = 0
    pinfo%jloc_tbl = 0

    do i = 1, n
    do j = 1, n
      call INFOG2L( i, j, pinfo%desca, pinfo%nprow, pinfo%npcol, pinfo%myrow, pinfo%mycol, i_loc, j_loc, proc_row, proc_col )
      if (pinfo%myrow == proc_row .and. pinfo%mycol == proc_col) then
        ip = info%irank_io(j)
        icount( ip ) = icount( ip ) + 1
        pinfo%i_tbl( icount(ip), ip ) = i
        pinfo%j_tbl( icount(ip), ip ) = j
        pinfo%iloc_tbl( icount(ip), ip ) = i_loc
        pinfo%jloc_tbl( icount(ip), ip ) = j_loc
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

      allocate( h_div(pinfo%nrow_local,pinfo%ncol_local), &
                v_div(pinfo%nrow_local,pinfo%ncol_local), e(n) )

      trilwmin = 3*n + max( mb*(pinfo%nrow_local+1), 3*mb )
      len_work0 = max( 1+6*n+2*pinfo%nrow_local*pinfo%ncol_local, trilwmin )

      call PDSYEVD( 'V', 'L', n, h_div, 1, 1, pinfo%desca, e, v_div, 1, 1, pinfo%descz, &
                    rtmp, -1, itmp, len_iwork, ierr )

      pinfo%len_work = max( nint(rtmp(1))*10, len_work0*10 )

      deallocate( e, v_div, h_div )
    end subroutine

    subroutine set_from_pzheevd
      implicit none
      complex(8), allocatable :: h_div(:,:), v_div(:,:)
      real(8), allocatable    :: e(:)
      real(8)    :: rtmp(1)
      complex(8) :: ctmp(1)
      integer    :: itmp(1)

      allocate( h_div(pinfo%nrow_local,pinfo%ncol_local), &
                v_div(pinfo%nrow_local,pinfo%ncol_local), e(n) )

      call PZHEEVD( 'V', 'L', n, h_div, 1, 1, pinfo%desca, e, v_div, 1, 1, pinfo%descz, &
                    ctmp, -1, rtmp, -1, itmp, len_iwork, ierr )

      pinfo%len_work = max( pinfo%len_work, n+(pinfo%nrow_local+pinfo%ncol_local+mb)*mb )
      pinfo%len_work = max( pinfo%len_work, nint(real(ctmp(1))) )

      pinfo%len_rwork = max( nint(rtmp(1)), (1+8*n+2*pinfo%nrow_local*pinfo%ncol_local)*2 )

      deallocate( e, v_div, h_div )
    end subroutine
  end subroutine init_blacs

end module scalapack_module
