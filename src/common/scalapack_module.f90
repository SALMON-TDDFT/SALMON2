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
    use structures, only: s_process_info, s_orbital_parallel
    implicit none
    type(s_process_info),intent(inout)  :: pinfo
    type(s_orbital_parallel),intent(in) :: info

    integer,parameter :: iclose_comm = 2 ! 1: rko, 2: r
    integer :: ii,ix,iy,iz
    integer :: k1,k2

    select case(iclose_comm)
    ! close on icomm_rko
    case(1)
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
      ii = 0
      do iz=0,pinfo%nprgrid(3)-1
      do iy=0,pinfo%nprgrid(2)-1
      do ix=0,pinfo%nprgrid(1)-1
        ii = ii + 1
        k1 = mod(ii-1,pinfo%nprow)
        k2 = (ii-1)/pinfo%nprow
        pinfo%gridmap(k1+1,k2+1) = info%imap(ix,iy,iz,info%iaddress(4),info%iaddress(5))
      end do
      end do
      end do
      if (ii /= info%isize_r) &
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
    use structures, only: s_process_info, s_orbital_parallel
    use communication, only: comm_is_root
    implicit none
    type(s_process_info),intent(in)     :: pinfo
    type(s_orbital_parallel),intent(in) :: info
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
    use structures, only: s_process_info, s_orbital_parallel
    implicit none
    integer :: NUMROC

    type(s_process_info),intent(inout)  :: pinfo
    type(s_orbital_parallel),intent(in) :: info
    integer,intent(in) :: m

    integer :: n,mb,nb
    integer :: len_iwork,len_work0,trilwmin
    integer :: ictxt,ierr
    real(8) :: rtmp(1)
    integer, allocatable :: iwork(:)
    real(8), allocatable :: h_div(:,:), v_div(:,:), e(:)

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

    len_iwork = 2 + 7*n + 8*pinfo%npcol

    ! determine the working memory size from PDSYEVD
    allocate( h_div(pinfo%nrow_local,pinfo%ncol_local), &
              v_div(pinfo%nrow_local,pinfo%ncol_local), &
              e(n), iwork(len_iwork) )
    trilwmin = 3*n + max( mb*(pinfo%nrow_local+1), 3*mb )
    len_work0 = max( 1+6*n+2*pinfo%nrow_local*pinfo%ncol_local, trilwmin )
    call PDSYEVD( 'V', 'L', n, h_div, 1, 1, pinfo%desca, e, v_div, 1, 1, pinfo%descz, &
                  rtmp, -1, iwork, len_iwork, ierr )
    pinfo%len_work = max( nint(rtmp(1))*10, len_work0*10 )
    deallocate( iwork, e, v_div, h_div )

    pinfo%flag_blacs_gridinit = .true.
  end subroutine init_blacs

end module scalapack_module
