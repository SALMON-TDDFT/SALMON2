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


module sendrecv_grid
  use structures, only: s_rgrid

  implicit none

  public :: init_sendrecv_grid4d
  public :: alloc_cache_real8
  public :: update_overlap
  public :: s_sendrecv_grid4d


  integer, parameter :: iside_up   = 1
  integer, parameter :: iside_down = 2
  integer, parameter :: itype_send = 1
  integer, parameter :: itype_recv = 2

  ! TODO: Move type defination to "common/structures.f90"
  type s_pcomm_cache4d
    real(8), allocatable :: dbuf(:, :, :, :)
    complex(8), allocatable :: zbuf(:, :, :, :)
  end type s_pcomm_cache4d

  ! TODO: Move type defination to "common/structures.f90"
  type s_sendrecv_grid4d
    ! Size of grid system
    type(s_rgrid) :: rg
    ! Number of orbitals (4-th dimension of grid)
    integer :: nb
    ! Communicator
    integer :: icomm, myrank
    ! Neightboring MPI id (1:x,2:y,3:z, 1:upside,2:downside):
    integer :: neig(1:3, 1:2) 
    ! Communication requests (1:x,2:y,3:z, 1:upside,2:downside, 1:send,2:recv):
    integer :: ireq(1:3, 1:2, 1:2)
    ! PComm cache (1:x,2:y,3:z, 1:upside,2:downside, 1:src/2:dst)
    type(s_pcomm_cache4d) :: cache(1:3, 1:2, 1:2)
    ! Range (dim=1:x,2:y,3:z, dir=1:upside,2:downside, 1:src/2:dst, axis=1...3)
    integer :: is_block(1:3, 1:2, 1:2, 1:3)
    integer :: ie_block(1:3, 1:2, 1:2, 1:3)
    logical :: pcomm_initialized
  end type s_sendrecv_grid4d

  interface update_overlap
  module procedure update_overlap_array4d_real8
  module procedure update_overlap_array4d_complex8
  end interface


  contains

  ! Flip 1 to 2, 2 to 1:
  integer function flip(i)
    implicit none
    integer, intent(in) :: i
    flip = 3 - i
  end function flip

  ! Prepare unique MPI tag number (3~8) for send/recv based on
  ! the communication direction (idir, iside) from the sender:
  integer function tag(idir, iside)
    implicit none
    integer, intent(in) :: idir, iside
    tag = 2 * idir + iside
  end function tag

  ! Initializing s_sendrecv_grid4d structure:
  ! NOTE:
  ! * This subroutine is commonly used for real type and complex type.
  ! * The cache region MUST be allocated after this initialization 
  !   by using `alloc_cache_real8/complex8`.
  subroutine init_sendrecv_grid4d(srg, rg, nb, icomm, myrank, neig)
    implicit none
    type(s_sendrecv_grid4d), intent(inout) :: srg
    type(s_rgrid), intent(in) :: rg
    integer, intent(in) ::  nb, icomm, myrank, neig(1:3, 1:2)

    integer :: idir, iaxis
    integer :: is_block(1:3, 1:2, 1:2, 1:3)
    integer :: ie_block(1:3, 1:2, 1:2, 1:3)
    
    ! Calculate shape (upper and lower bounds) of overlapped region:
    ! NOTE:
    !
    !  nd nd         nd nd    U/D: up/down side of `idir`-th direction
    ! +==+--+-------+--+==+   S/R: send/recv region, representing the
    ! |DR|DS|       |US|UR|        inner and outer side of local grid.
    ! +==+--+-------+--+==+   
    !    is            ie     
    !    <------------->      
    !    Local grid size
    !
    ! is/ie_block(idir, iside, itype, iaxis): 
    !   lower/upper bounds of each blocks
    !   * idir: direction (1:x, 2:y, 3:z)
    !   * iside: kind of region (1:upside, 2:downside)
    !   * itype: kind of region (1:send, 2:recv)
    !   * iaxis: axis (1:x, 2:y, 3:z)
    do idir = 1, 3 ! 1:x,2:y,3:z
      do iaxis = 1, 3 ! 1:x,2:y,3:z
        if (idir == iaxis) then
          ! upside-send (US) block:
          is_block(idir, iside_up, itype_send, idir) = rg%ie(idir) - rg%nd + 1
          ie_block(idir, iside_up, itype_send, idir) = rg%ie(idir) 
          ! upside-recv (UR) block:
          is_block(idir, iside_up, itype_recv, idir) = rg%ie(idir) + 1
          ie_block(idir, iside_up, itype_recv, idir) = rg%ie(idir) + rg%nd
          ! downside-send (DS) block:
          is_block(idir, iside_down, itype_send, idir) = rg%is(idir)
          ie_block(idir, iside_down, itype_send, idir) = rg%is(idir) + rg%nd - 1
          ! downside-recv (DR) block:
          is_block(idir, iside_down, itype_recv, idir) = rg%is(idir) - rg%nd
          ie_block(idir, iside_down, itype_recv, idir) = rg%is(idir) - 1
        else
          is_block(idir, :, :, iaxis) = rg%is(iaxis)
          ie_block(idir, :, :, iaxis) = rg%ie(iaxis)
        end if
      end do
    end do

    ! Assign to s_sendrecv_grid4d structure:
    srg%rg = rg
    srg%nb = nb
    srg%neig = neig
    srg%is_block(:, :, :, :) = is_block
    srg%ie_block(:, :, :, :) = ie_block
    srg%icomm = icomm
    srg%myrank = myrank
    srg%ireq = -1
  end subroutine init_sendrecv_grid4d


  ! Allocate cache region for persistent communication:
  subroutine alloc_cache_real8(srg)
    implicit none
    type(s_sendrecv_grid4d), intent(inout) :: srg
    integer :: idir, iside, itype, is_b(3), ie_b(3)

    do idir = 1, 3 ! 1:x,2:y,3:z
      do iside = 1, 2 ! 1:up,2:down
        do itype = 1, 2 ! 1:send, 2:recv
            is_b(1:3) = srg%is_block(idir, iside, itype, 1:3)
            ie_b(1:3) = srg%ie_block(idir, iside, itype, 1:3)
            allocate(srg%cache(idir, iside, itype)%dbuf( &
              is_b(1):ie_b(1), is_b(2):ie_b(2), is_b(3):ie_b(3), 1:srg%nb))
        end do
      end do
    end do
    srg%pcomm_initialized = .false. ! Flag for persistent communication
  end subroutine


  ! Allocate cache region for persistent communication:
  subroutine alloc_cache_complex8(srg)
    implicit none
    type(s_sendrecv_grid4d), intent(inout) :: srg
    integer :: idir, iside, itype, is_b(3), ie_b(3)

    do idir = 1, 3 ! 1:x,2:y,3:z
      do iside = 1, 2 ! 1:up,2:down
        do itype = 1, 2 ! 1:send, 2:recv
            is_b(1:3) = srg%is_block(idir, iside, itype, 1:3)
            ie_b(1:3) = srg%ie_block(idir, iside, itype, 1:3)
            allocate(srg%cache(idir, iside, itype)%zbuf( &
              is_b(1):ie_b(1), is_b(2):ie_b(2), is_b(3):ie_b(3), 1:srg%nb))
        end do
      end do
    end do
    srg%pcomm_initialized = .false. ! Flag for persistent communication
  end subroutine


  subroutine update_overlap_array4d_real8(srg, data)
    use salmon_communication, only: comm_start_all, comm_wait_all, comm_proc_null
    implicit none
    type(s_sendrecv_grid4d), intent(inout) :: srg
    real(8), intent(inout) :: data( &
      srg%rg%is_array(1):srg%rg%ie_array(1), &
      srg%rg%is_array(2):srg%rg%ie_array(2), &
      srg%rg%is_array(3):srg%rg%ie_array(3), &
      1:srg%nb)
    integer :: idir, iside

    ! Exchange the overlap region with the neighboring node (or opposite side of itself).
    
    do idir = 1, 3 ! 1:x,2:y,3:z
      do iside = 1, 2 ! 1:up,2:down
        if (srg%neig(idir, iside) /= comm_proc_null) then
          if (srg%neig(idir, iside) /= srg%myrank) then
            ! Store the overlap reigion into the cache 
            call pack_cache(idir, iside) 
            ! In the first call of this subroutine, setup the persistent communication:
            if (.not. srg%pcomm_initialized) call init_pcomm(idir, iside)
            ! Start to communication
            call comm_start_all(srg%ireq(idir, iside, :))
          else
            ! NOTE: If neightboring nodes are itself (periodic with single proc),
            !       a simple side-to-side copy is used instead of the MPI comm.
            call copy_self(idir, iside)
          end if
        end if
      end do
    end do

    do idir = 1, 3 ! 1:x,2:y,3:z
      do iside = 1, 2 ! 1:up,2:down
        if (srg%neig(idir, iside) /= comm_proc_null) then
          if (srg%neig(idir, iside) /= srg%myrank) then
            ! Wait for recieving
            call comm_wait_all(srg%ireq(idir, iside, :))
            ! Write back the recieved cache
            call unpack_cache(idir, iside)
          end if
        end if
      end do
    end do

    srg%pcomm_initialized = .true. ! Update pcomm_initialized

    contains

    subroutine init_pcomm(jdir, jside)
      use salmon_communication, only: comm_send_init, comm_recv_init
      implicit none
      integer, intent(in) :: jdir, jside
      ! Send (and initialize persistent communication)
      srg%ireq(jdir, jside, itype_send) = comm_send_init( &
        srg%cache(jdir, jside, itype_send)%dbuf, &
        srg%neig(jdir, jside), &
        tag(jdir, jside), &
        srg%icomm)
      ! Recv (and initialize persistent communication)
      srg%ireq(jdir, jside, itype_recv) = comm_recv_init( &
        srg%cache(jdir, jside, itype_recv)%dbuf, &
        srg%neig(jdir, jside), &
        tag(jdir, flip(jside)), & ! `jside` in sender
        srg%icomm)
    end subroutine init_pcomm

    subroutine pack_cache(jdir, jside)
      use pack_unpack, only: copy_data
      implicit none
      integer, intent(in) :: jdir, jside
      integer :: is_s(1:3), ie_s(1:3) ! src region
      is_s(1:3) = srg%is_block(jdir, jside, itype_send, 1:3)
      ie_s(1:3) = srg%ie_block(jdir, jside, itype_send, 1:3)
      call copy_data( &
        data(is_s(1):ie_s(1), is_s(2):ie_s(2), is_s(3):ie_s(3), 1:srg%nb), &
        srg%cache(jdir, jside, itype_send)%dbuf)
    end subroutine pack_cache

    subroutine unpack_cache(jdir, jside)
      use pack_unpack, only: copy_data
      implicit none
      integer, intent(in) :: jdir, jside
      integer :: is_d(1:3), ie_d(1:3) ! dst region
      is_d(1:3) = srg%is_block(jdir, jside, itype_recv, 1:3)
      ie_d(1:3) = srg%ie_block(jdir, jside, itype_recv, 1:3)
      call copy_data( &
        srg%cache(jdir, jside, itype_recv)%dbuf, &
        data(is_d(1):ie_d(1), is_d(2):ie_d(2), is_d(3):ie_d(3), 1:srg%nb))
    end subroutine unpack_cache

    subroutine copy_self(jdir, jside)
      use pack_unpack, only: copy_data
      integer, intent(in) :: jdir, jside
      integer :: is_s(1:3), ie_s(1:3) ! src region
      integer :: is_d(1:3), ie_d(1:3) ! dst region
      is_s(1:3) = srg%is_block(jdir, flip(jside), itype_send, 1:3)
      ie_s(1:3) = srg%ie_block(jdir, flip(jside), itype_send, 1:3)
      is_d(1:3) = srg%is_block(jdir, jside, itype_recv, 1:3)
      ie_d(1:3) = srg%ie_block(jdir, jside, itype_recv, 1:3)
      call copy_data( &
        data(is_s(1):ie_s(1), is_s(2):ie_s(2), is_s(3):ie_s(3), 1:srg%nb), &
        data(is_d(1):ie_d(1), is_d(2):ie_d(2), is_d(3):ie_d(3), 1:srg%nb))
    end subroutine copy_self

  end subroutine update_overlap_array4d_real8


  subroutine update_overlap_array4d_complex8(srg, data)
    use salmon_communication, only: comm_start_all, comm_wait_all, comm_proc_null
    implicit none
    type(s_sendrecv_grid4d), intent(inout) :: srg
    complex(8), intent(inout) :: data( &
      srg%rg%is_array(1):srg%rg%ie_array(1), &
      srg%rg%is_array(2):srg%rg%ie_array(2), &
      srg%rg%is_array(3):srg%rg%ie_array(3), &
      1:srg%nb)
    integer :: idir, iside

    ! Exchange the overlap region with the neighboring node (or opposite side of itself).
    
    do idir = 1, 3 ! 1:x,2:y,3:z
      do iside = 1, 2 ! 1:up,2:down
        if (srg%neig(idir, iside) /= comm_proc_null) then
          if (srg%neig(idir, iside) /= srg%myrank) then
            ! Store the overlap reigion into the cache 
            call pack_cache(idir, iside) 
            ! In the first call of this subroutine, setup the persistent communication:
            if (.not. srg%pcomm_initialized) call init_pcomm(idir, iside)
            ! Start to communication
            call comm_start_all(srg%ireq(idir, iside, :))
          else
            ! NOTE: If neightboring nodes are itself (periodic with single proc),
            !       a simple side-to-side copy is used instead of the MPI comm.
            call copy_self(idir, iside)
          end if
        end if
      end do
    end do

    do idir = 1, 3 ! 1:x,2:y,3:z
      do iside = 1, 2 ! 1:up,2:down
        if (srg%neig(idir, iside) /= comm_proc_null) then
          if (srg%neig(idir, iside) /= srg%myrank) then
            ! Wait for recieving
            call comm_wait_all(srg%ireq(idir, iside, :))
            ! Write back the recieved cache
            call unpack_cache(idir, iside)
          end if
        end if
      end do
    end do

    srg%pcomm_initialized = .true. ! Update pcomm_initialized

    contains

    subroutine init_pcomm(jdir, jside)
      use salmon_communication, only: comm_send_init, comm_recv_init
      implicit none
      integer, intent(in) :: jdir, jside
      ! Send (and initialize persistent communication)
      srg%ireq(jdir, jside, itype_send) = comm_send_init( &
        srg%cache(jdir, jside, itype_send)%zbuf, &
        srg%neig(jdir, jside), &
        tag(jdir, jside), &
        srg%icomm)
      ! Recv (and initialize persistent communication)
      srg%ireq(jdir, jside, itype_recv) = comm_recv_init( &
        srg%cache(jdir, jside, itype_recv)%zbuf, &
        srg%neig(jdir, jside), &
        tag(jdir, flip(jside)), & ! `jside` in sender
        srg%icomm)
    end subroutine init_pcomm

    subroutine pack_cache(jdir, jside)
      use pack_unpack, only: copy_data
      implicit none
      integer, intent(in) :: jdir, jside
      integer :: is_s(1:3), ie_s(1:3) ! src region
      is_s(1:3) = srg%is_block(jdir, jside, itype_send, 1:3)
      ie_s(1:3) = srg%ie_block(jdir, jside, itype_send, 1:3)
      call copy_data( &
        data(is_s(1):ie_s(1), is_s(2):ie_s(2), is_s(3):ie_s(3), 1:srg%nb), &
        srg%cache(jdir, jside, itype_send)%zbuf)
    end subroutine pack_cache

    subroutine unpack_cache(jdir, jside)
      use pack_unpack, only: copy_data
      implicit none
      integer, intent(in) :: jdir, jside
      integer :: is_d(1:3), ie_d(1:3) ! dst region
      is_d(1:3) = srg%is_block(jdir, jside, itype_recv, 1:3)
      ie_d(1:3) = srg%ie_block(jdir, jside, itype_recv, 1:3)
      call copy_data( &
        srg%cache(jdir, jside, itype_recv)%zbuf, &
        data(is_d(1):ie_d(1), is_d(2):ie_d(2), is_d(3):ie_d(3), 1:srg%nb))
    end subroutine unpack_cache

    subroutine copy_self(jdir, jside)
      use pack_unpack, only: copy_data
      integer, intent(in) :: jdir, jside
      integer :: is_s(1:3), ie_s(1:3) ! src region
      integer :: is_d(1:3), ie_d(1:3) ! dst region
      is_s(1:3) = srg%is_block(jdir, flip(jside), itype_send, 1:3)
      ie_s(1:3) = srg%ie_block(jdir, flip(jside), itype_send, 1:3)
      is_d(1:3) = srg%is_block(jdir, jside, itype_recv, 1:3)
      ie_d(1:3) = srg%ie_block(jdir, jside, itype_recv, 1:3)
      call copy_data( &
        data(is_s(1):ie_s(1), is_s(2):ie_s(2), is_s(3):ie_s(3), 1:srg%nb), &
        data(is_d(1):ie_d(1), is_d(2):ie_d(2), is_d(3):ie_d(3), 1:srg%nb))
    end subroutine copy_self

  end subroutine update_overlap_array4d_complex8


end module sendrecv_grid


