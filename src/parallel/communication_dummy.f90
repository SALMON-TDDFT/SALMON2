!
!  Copyright 2019-2020 SALMON developers
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
! NOTE: this is a dummy for single node only application
module communication
  implicit none

  integer, private, parameter :: DEAD_BEEF       = int(z'7FFFDEAD')
  integer, private, parameter :: COMM_WORLD_ID   = int(z'7FFFFFFF')

  integer, public, parameter  :: COMM_GROUP_NULL = DEAD_BEEF
  integer, public, parameter  :: ROOT_PROCID     = 0
  integer, public, parameter  :: COMM_PROC_NULL  = DEAD_BEEF

  ! call once
  public :: comm_init
  public :: comm_finalize

  ! p2p communication
  ! application stops when a following routines is called in no-mpi environment
  public :: comm_send
  public :: comm_recv
  public :: comm_exchange

  ! p2p immediate communication
  ! application stops when a following routines is called in no-mpi environment
  public :: comm_isend
  public :: comm_irecv
  public :: comm_wait
  public :: comm_wait_all

  ! p2p persistent communication
  ! application stops when a following routines is called in no-mpi environment
  public :: comm_send_init
  public :: comm_recv_init
  public :: comm_start_all
  public :: comm_free_reqs

  ! collective communication
  public :: comm_sync_all
  public :: comm_summation
  public :: comm_bcast
  public :: comm_allgather
  public :: comm_allgatherv ! not implemented in no-mpi environment
  public :: comm_alltoall
  public :: comm_get_min
  public :: comm_get_max

  ! group (communicator)
  public :: comm_get_globalinfo
  public :: comm_get_groupinfo
  public :: comm_create_group
  public :: comm_create_group_byid
  public :: comm_free_group

  ! utils
  public :: comm_is_root
  public :: comm_show_error


  type, public :: comm_maxloc_type
    real(8) :: val
    integer :: rank
  end type

  interface comm_send
    ! 4-D array
    module procedure comm_send_array4d_double
    module procedure comm_send_array4d_dcomplex

    ! 5-D array
    module procedure comm_send_array5d_double
    module procedure comm_send_array5d_dcomplex
  end interface

  interface comm_recv
    ! 4-D array
    module procedure comm_recv_array4d_double
    module procedure comm_recv_array4d_dcomplex

    ! 5-D array
    module procedure comm_recv_array5d_double
    module procedure comm_recv_array5d_dcomplex
  end interface

  interface comm_exchange
    ! 3-D array
    module procedure comm_exchange_array3d_double
    module procedure comm_exchange_array3d_dcomplex

    ! 5-D array
    module procedure comm_exchange_array5d_double
    module procedure comm_exchange_array5d_dcomplex
  end interface

  interface comm_isend
    ! 3-D array
    module procedure comm_isend_array3d_double
    module procedure comm_isend_array3d_dcomplex

    ! 5-D array
    module procedure comm_isend_array5d_double
    module procedure comm_isend_array5d_dcomplex
  end interface

  interface comm_irecv
    ! 3-D array
    module procedure comm_irecv_array3d_double
    module procedure comm_irecv_array3d_dcomplex

    ! 5-D array
    module procedure comm_irecv_array5d_double
    module procedure comm_irecv_array5d_dcomplex
  end interface

  interface comm_send_init
    ! 3-D array
    module procedure comm_send_init_array3d_double
    module procedure comm_send_init_array3d_dcomplex

    ! 4-D array
    module procedure comm_send_init_array4d_double
    module procedure comm_send_init_array4d_dcomplex

    ! 5-D array
    module procedure comm_send_init_array5d_double
    module procedure comm_send_init_array5d_dcomplex
  end interface

  interface comm_recv_init
    ! 3-D array
    module procedure comm_recv_init_array3d_double
    module procedure comm_recv_init_array3d_dcomplex

    ! 4-D array
    module procedure comm_recv_init_array4d_double
    module procedure comm_recv_init_array4d_dcomplex

    ! 5-D array
    module procedure comm_recv_init_array5d_double
    module procedure comm_recv_init_array5d_dcomplex
  end interface

  interface comm_summation
    ! scalar
    module procedure comm_summation_integer
    module procedure comm_summation_double
    module procedure comm_summation_dcomplex

    ! 1-D array
    module procedure comm_summation_array1d_integer
    module procedure comm_summation_array1d_double
    module procedure comm_summation_array1d_dcomplex

    ! 2-D array
    module procedure comm_summation_array2d_integer
    module procedure comm_summation_array2d_double
    module procedure comm_summation_array2d_dcomplex

    ! 3-D array
    module procedure comm_summation_array3d_integer
    module procedure comm_summation_array3d_double
    module procedure comm_summation_array3d_dcomplex

    ! 4-D array
    module procedure comm_summation_array4d_double
    module procedure comm_summation_array4d_dcomplex

    ! 5-D array
    module procedure comm_summation_array5d_double
    module procedure comm_summation_array5d_dcomplex

    ! 6-D array
    module procedure comm_summation_array6d_double
    module procedure comm_summation_array6d_dcomplex

    ! in-place
    module procedure comm_sum_ip_array1d_integer
    module procedure comm_sum_ip_array2d_integer
    module procedure comm_sum_ip_array3d_integer
    module procedure comm_sum_ip_array3d_double
    module procedure comm_sum_ip_array5d_integer
  end interface

  interface comm_bcast
    ! scalar
    module procedure comm_bcast_integer
    module procedure comm_bcast_double
    module procedure comm_bcast_character
    module procedure comm_bcast_logical

    ! 1-D array
    module procedure comm_bcast_array1d_integer
    module procedure comm_bcast_array1d_double
    module procedure comm_bcast_array1d_character

    ! 2-D array
    module procedure comm_bcast_array2d_integer
    module procedure comm_bcast_array2d_double
    module procedure comm_bcast_array2d_character

    ! 3-D array
    module procedure comm_bcast_array3d_double
    module procedure comm_bcast_array3d_dcomplex

    ! 4-D array
    module procedure comm_bcast_array4d_double
    module procedure comm_bcast_array4d_dcomplex

    ! 5-D array
    module procedure comm_bcast_array5d_double
    module procedure comm_bcast_array5d_dcomplex
  end interface

  interface comm_allgather
    ! 1-D array
    module procedure comm_allgather_array1d_logical
  end interface

  interface comm_allgatherv
    ! 1-D array
    module procedure comm_allgatherv_array1d_double
  end interface

  interface comm_alltoall
    ! 1-D array
    module procedure comm_alltoall_array1d_complex
  end interface

  interface comm_get_min
    ! scalar
    module procedure comm_get_min_double

    ! 1-D array
    module procedure comm_get_min_array1d_double
  end interface

  interface comm_get_max
    ! scalar
    module procedure comm_get_maxloc
    module procedure comm_get_max_integer

    ! 1-D array
    module procedure comm_get_max_array1d_double
  end interface

  interface comm_logical_and
    ! scalar
    module procedure comm_logical_and_scalar
  end interface

  interface comm_logical_or
    ! 1-D array (in-place)
    module procedure comm_ip_logical_or_array1d
  end interface

  private :: get_rank, error_check, abort_show_message

#define ABORT_MESSAGE(target,msg) if(target/=COMM_PROC_NULL) call abort_show_message(msg)
#define UNUSED_VARIABLE(VAR)      if(.false.) call salmon_unusedvar(VAR)

contains
  subroutine comm_init
    implicit none
    ! no operation
  end subroutine

  subroutine comm_finalize
    implicit none
    ! no operation
  end subroutine

  subroutine comm_get_globalinfo(ngid, npid, nprocs)
    implicit none
    integer, intent(out) :: ngid, npid, nprocs
    ngid = COMM_WORLD_ID
    call get_rank(ngid, npid, nprocs)
  end subroutine

  subroutine comm_get_groupinfo(ngid, npid, nprocs)
    implicit none
    integer, intent(in)  :: ngid
    integer, intent(out) :: npid, nprocs
    call get_rank(ngid, npid, nprocs)
  end subroutine

  function comm_create_group(ngid, nprocs, key) result(ngid_dst)
    implicit none
    integer, intent(in) :: ngid, nprocs, key
    integer :: ngid_dst
    UNUSED_VARIABLE(key)
    UNUSED_VARIABLE(nprocs)
    ngid_dst = ngid
  end function

  function comm_create_group_byid(iparent, idlists) result(ichild)
    implicit none
    integer, intent(in) :: iparent    ! parent communicator
    integer, intent(in) :: idlists(:) ! include ranks in new communicator
    integer :: ichild
    UNUSED_VARIABLE(idlists)
    ichild = iparent
  end function

  subroutine comm_free_group(igroup)
    implicit none
    integer, intent(in) :: igroup
    UNUSED_VARIABLE(igroup)
    ! no operation
  end subroutine

  function comm_is_root(npid)
    implicit none
    integer, intent(in) :: npid
    logical :: comm_is_root
    comm_is_root = npid == ROOT_PROCID
  end function

  subroutine comm_show_error(errcode)
    implicit none
    integer, intent(in) :: errcode
    UNUSED_VARIABLE(errcode)
  end subroutine

  subroutine comm_sync_all(ngid)
    implicit none
    integer, intent(in), optional :: ngid
    UNUSED_VARIABLE(ngid)
    ! no operation
  end subroutine

  subroutine comm_send_array4d_double(invalue, ndest, ntag, ngroup)
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_array4d_double")
  end subroutine

  subroutine comm_send_array4d_dcomplex(invalue, ndest, ntag, ngroup)
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_array4d_dcomplex")
  end subroutine

  subroutine comm_recv_array4d_double(outvalue, nsrc, ntag, ngroup)
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_array4d_double")
  end subroutine

  subroutine comm_recv_array4d_dcomplex(outvalue, nsrc, ntag, ngroup)
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_array4d_dcomplex")
  end subroutine

  subroutine comm_send_array5d_double(invalue, ndest, ntag, ngroup)
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_array5d_double")
  end subroutine

  subroutine comm_send_array5d_dcomplex(invalue, ndest, ntag, ngroup)
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_array5d_dcomplex")
  end subroutine

  subroutine comm_recv_array5d_double(outvalue, nsrc, ntag, ngroup)
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_array5d_double")
  end subroutine

  subroutine comm_recv_array5d_dcomplex(outvalue, nsrc, ntag, ngroup)
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_array5d_dcomplex")
  end subroutine


  subroutine comm_exchange_array3d_double(invalue, ndest, outvalue, nsrc, ntag, ngroup)
    implicit none
    real(8), intent(in)  :: invalue(:,:,:)
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: nsrc, ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(nsrc)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_exchange_array3d_double")
  end subroutine

  subroutine comm_exchange_array3d_dcomplex(invalue, ndest, outvalue, nsrc, ntag, ngroup)
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: nsrc, ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(nsrc)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_exchange_array3d_dcomplex")
  end subroutine

  subroutine comm_exchange_array5d_double(invalue, ndest, outvalue, nsrc, ntag, ngroup)
    implicit none
    real(8), intent(in)  :: invalue(:,:,:,:,:)
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(nsrc)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_exchange_array5d_double")
  end subroutine

  subroutine comm_exchange_array5d_dcomplex(invalue, ndest, outvalue, nsrc, ntag, ngroup)
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(nsrc)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_exchange_array5d_dcomplex")
  end subroutine


  function comm_isend_array3d_double(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    real(8), intent(in) :: invalue(:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_isend_array3d_double")
    req = DEAD_BEEF
  end function

  function comm_isend_array3d_dcomplex(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(in) :: invalue(:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_isend_array3d_dcomplex")
    req = DEAD_BEEF
  end function

  function comm_isend_array5d_double(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_isend_array5d_double")
    req = DEAD_BEEF
  end function

  function comm_isend_array5d_dcomplex(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_isend_array5d_dcomplex")
    req = DEAD_BEEF
  end function

  function comm_irecv_array3d_double(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_irecv_array3d_double")
    req = DEAD_BEEF
  end function

  function comm_irecv_array3d_dcomplex(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_irecv_array3d_dcomplex")
    req = DEAD_BEEF
  end function

  function comm_irecv_array5d_double(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_irecv_array5d_double")
    req = DEAD_BEEF
  end function

  function comm_irecv_array5d_dcomplex(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_irecv_array5d_dcomplex")
    req = DEAD_BEEF
  end function


  subroutine comm_wait(req)
    implicit none
    integer, intent(in) :: req
    ABORT_MESSAGE(req,"comm_wait")
  end subroutine

  subroutine comm_wait_all(reqs)
    implicit none
    integer, intent(in) :: reqs(:)
    UNUSED_VARIABLE(reqs)
    ! do nothing
  end subroutine

  function comm_send_init_array3d_double(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    real(8), intent(in) :: invalue(:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_init_array3d_double")
    req = DEAD_BEEF
  end function

  function comm_send_init_array3d_dcomplex(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(in) :: invalue(:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_init_array3d_dcomplex")
    req = DEAD_BEEF
  end function

  function comm_send_init_array4d_double(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_init_array4d_double")
    req = DEAD_BEEF
  end function

  function comm_send_init_array4d_dcomplex(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_init_array4d_dcomplex")
    req = DEAD_BEEF
  end function

  function comm_send_init_array5d_double(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_init_array5d_double")
    req = DEAD_BEEF
  end function

  function comm_send_init_array5d_dcomplex(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_init_array5d_dcomplex")
    req = DEAD_BEEF
  end function

  function comm_recv_init_array3d_double(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_init_array3d_double")
    req = DEAD_BEEF
  end function

  function comm_recv_init_array3d_dcomplex(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_init_array3d_dcomplex")
    req = DEAD_BEEF
  end function

  function comm_recv_init_array4d_double(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_init_array4d_double")
    req = DEAD_BEEF
  end function

  function comm_recv_init_array4d_dcomplex(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_init_array4d_dcomplex")
    req = DEAD_BEEF
  end function

  function comm_recv_init_array5d_double(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_init_array5d_double")
    req = DEAD_BEEF
  end function

  function comm_recv_init_array5d_dcomplex(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_init_array5d_dcomplex")
    req = DEAD_BEEF
  end function

  subroutine comm_start_all(reqs)
    implicit none
    integer, intent(in) :: reqs(:)
    UNUSED_VARIABLE(reqs)
    ! do nothing
  end subroutine

  subroutine comm_free_reqs(reqs)
    implicit none
    integer, intent(in) :: reqs(:)
    UNUSED_VARIABLE(reqs)
    ! do nothing
  end subroutine

  subroutine comm_summation_integer(invalue, outvalue, ngroup, dest)
    implicit none
    integer, intent(in)  :: invalue
    integer, intent(out) :: outvalue
    integer, intent(in)  :: ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_integer")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_double(invalue, outvalue, ngroup, dest)
    implicit none
    real(8), intent(in)  :: invalue
    real(8), intent(out) :: outvalue
    integer, intent(in)  :: ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_double")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_dcomplex(invalue, outvalue, ngroup, dest)
    implicit none
    complex(8), intent(in)  :: invalue
    complex(8), intent(out) :: outvalue
    integer, intent(in)     :: ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_dcomplex")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array1d_integer(invalue, outvalue, N, ngroup, dest)
    implicit none
    integer, intent(in)  :: invalue(:)
    integer, intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_array1d_integer")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array1d_double(invalue, outvalue, N, ngroup, dest)
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_array1d_double")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array1d_dcomplex(invalue, outvalue, N, ngroup, dest)
    implicit none
    complex(8), intent(in)  :: invalue(:)
    complex(8), intent(out) :: outvalue(:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_array1d_dcomplex")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array2d_integer(invalue, outvalue, N, ngroup, dest)
    implicit none
    integer, intent(in)  :: invalue(:,:)
    integer, intent(out) :: outvalue(:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_array2d_integer")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array2d_double(invalue, outvalue, N, ngroup, dest)
    implicit none
    real(8), intent(in)  :: invalue(:,:)
    real(8), intent(out) :: outvalue(:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_array2d_double")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array2d_dcomplex(invalue, outvalue, N, ngroup, dest)
    implicit none
    complex(8), intent(in)  :: invalue(:,:)
    complex(8), intent(out) :: outvalue(:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_array2d_dcomplex")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array3d_integer(invalue, outvalue, N, ngroup, dest)
    implicit none
    integer, intent(in)  :: invalue(:,:,:)
    integer, intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_array3d_integer")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array3d_double(invalue, outvalue, N, ngroup, dest)
    implicit none
    real(8), intent(in)  :: invalue(:,:,:)
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_array3d_double")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array3d_dcomplex(invalue, outvalue, N, ngroup, dest)
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_array3d_dcomplex")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array4d_double(invalue, outvalue, N, ngroup, dest)
    implicit none
    real(8), intent(in)  :: invalue(:,:,:,:)
    real(8), intent(out) :: outvalue(:,:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_array4d_double")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array4d_dcomplex(invalue, outvalue, N, ngroup, dest)
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_array4d_dcomplex")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array5d_double(invalue, outvalue, N, ngroup, dest)
    implicit none
    real(8), intent(in)  :: invalue(:,:,:,:,:)
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_array5d_double")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array5d_dcomplex(invalue, outvalue, N, ngroup, dest)
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_array5d_dcomplex")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array6d_double(invalue, outvalue, N, ngroup, dest)
    implicit none
    real(8), intent(in)  :: invalue(:,:,:,:,:,:)
    real(8), intent(out) :: outvalue(:,:,:,:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_array6d_double")
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array6d_dcomplex(invalue, outvalue, N, ngroup, dest)
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:,:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:,:,:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    !NOP! ABORT_MESSAGE(ngroup,"comm_summation_array6d_dcomplex")
    outvalue = invalue
  end subroutine


  subroutine comm_sum_ip_array1d_integer(values, ngroup)
    implicit none
    integer, intent(inout) :: values(:)
    integer, intent(in)    :: ngroup
    UNUSED_VARIABLE(values)
    UNUSED_VARIABLE(ngroup)
    !NOP! ABORT_MESSAGE(ngroup,"comm_sum_ip_array1d_integer")
  end subroutine

  subroutine comm_sum_ip_array2d_integer(values, ngroup)
    implicit none
    integer, intent(inout) :: values(:,:)
    integer, intent(in)    :: ngroup
    UNUSED_VARIABLE(values)
    UNUSED_VARIABLE(ngroup)
    !NOP! ABORT_MESSAGE(ngroup,"comm_sum_ip_array2d_integer")
  end subroutine

  subroutine comm_sum_ip_array3d_double(values, ngroup)
    implicit none
    real(8), intent(inout) :: values(:,:,:)
    integer, intent(in)    :: ngroup
    UNUSED_VARIABLE(values)
    UNUSED_VARIABLE(ngroup)
    !NOP! ABORT_MESSAGE(ngroup,"comm_sum_ip_array3d_double")
  end subroutine

  subroutine comm_sum_ip_array3d_integer(values, ngroup)
    implicit none
    integer, intent(inout) :: values(:,:,:)
    integer, intent(in)    :: ngroup
    UNUSED_VARIABLE(values)
    UNUSED_VARIABLE(ngroup)
    !NOP! ABORT_MESSAGE(ngroup,"comm_sum_ip_array3d_integer")
  end subroutine

  subroutine comm_sum_ip_array5d_integer(values, ngroup)
    implicit none
    integer, intent(inout) :: values(:,:,:,:,:)
    integer, intent(in)    :: ngroup
    UNUSED_VARIABLE(values)
    UNUSED_VARIABLE(ngroup)
    !NOP! ABORT_MESSAGE(ngroup,"comm_sum_ip_array5d_integer")
  end subroutine



  subroutine comm_bcast_integer(val, ngroup, root)
    implicit none
    integer, intent(inout)        :: val
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_integer")
  end subroutine

  subroutine comm_bcast_double(val, ngroup, root)
    implicit none
    real(8), intent(inout)        :: val
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_double")
  end subroutine

  subroutine comm_bcast_character(val, ngroup, root)
    implicit none
    character(*), intent(inout)        :: val
    integer,      intent(in)           :: ngroup
    integer,      intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_character")
  end subroutine

  subroutine comm_bcast_logical(val, ngroup, root)
    implicit none
    logical, intent(inout)        :: val
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_logical")
  end subroutine

  subroutine comm_bcast_array1d_integer(val, ngroup, root)
    implicit none
    integer, intent(inout)        :: val(:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_array1d_integer")
  end subroutine

  subroutine comm_bcast_array1d_double(val, ngroup, root)
    implicit none
    real(8), intent(inout)        :: val(:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_array1d_double")
  end subroutine

  subroutine comm_bcast_array2d_integer(val, ngroup, root)
    implicit none
    integer, intent(inout)        :: val(:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_array2d_integer")
  end subroutine

  subroutine comm_bcast_array2d_double(val, ngroup, root)
    implicit none
    real(8), intent(inout)        :: val(:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_array2d_double")
  end subroutine

  subroutine comm_bcast_array3d_double(val, ngroup, root)
    implicit none
    real(8), intent(inout)        :: val(:,:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_array3d_double")
  end subroutine

  subroutine comm_bcast_array4d_double(val, ngroup, root)
    implicit none
    real(8), intent(inout)        :: val(:,:,:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_array4d_double")
  end subroutine

  subroutine comm_bcast_array5d_double(val, ngroup, root)
    implicit none
    real(8), intent(inout)        :: val(:,:,:,:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_array5d_double")
  end subroutine

  subroutine comm_bcast_array3d_dcomplex(val, ngroup, root)
    implicit none
    complex(8), intent(inout)     :: val(:,:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_array3d_dcomplex")
  end subroutine

  subroutine comm_bcast_array4d_dcomplex(val, ngroup, root)
    implicit none
    complex(8), intent(inout)     :: val(:,:,:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_array4d_dcomplex")
  end subroutine

  subroutine comm_bcast_array5d_dcomplex(val, ngroup, root)
    implicit none
    complex(8), intent(inout)     :: val(:,:,:,:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_array5d_dcomplex")
  end subroutine

  subroutine comm_bcast_array1d_character(val, ngroup, root)
    implicit none
    character(*), intent(inout)        :: val(:)
    integer,      intent(in)           :: ngroup
    integer,      intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_array1d_character")
  end subroutine

  subroutine comm_bcast_array2d_character(val, ngroup, root)
    implicit none
    character(*), intent(inout)        :: val(:,:)
    integer,      intent(in)           :: ngroup
    integer,      intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    !NOP! ABORT_MESSAGE(ngroup,"comm_bcast_array2d_character")
  end subroutine


  subroutine comm_allgather_array1d_logical(invalue, outvalue, ngroup)
    implicit none
    logical, intent(in)  :: invalue(:)
    logical, intent(out) :: outvalue(:,:)
    integer, intent(in)  :: ngroup
    UNUSED_VARIABLE(ngroup)
    !NOP! ABORT_MESSAGE(ngroup,"comm_allgather_array1d_logical")
    outvalue(:,1) = invalue(:)
  end subroutine


  subroutine comm_allgatherv_array1d_double(invalue, outvalue, ncounts, displs, ngroup)
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: ncounts(:)
    integer, intent(in)  :: displs(:)
    integer, intent(in)  :: ngroup
    UNUSED_VARIABLE(ngroup)
    !NOP! ABORT_MESSAGE(ngroup,"comm_allgatherv_array1d_double")
    outvalue(1+displs(1):1+displs(1)+ncounts(1)-1) = invalue(1+0:1+ncounts(1)-1)
  end subroutine


  subroutine comm_alltoall_array1d_complex(invalue, outvalue, ngroup, ncount)
    implicit none
    complex(8), intent(in)  :: invalue(:)
    complex(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: ngroup
    integer, intent(in)  :: ncount
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(ncount)
    !NOP! ABORT_MESSAGE(ngroup,"comm_alltoall_array1d_complex")
    outvalue = invalue
  end subroutine

  subroutine comm_get_min_double(svalue, ngroup)
    implicit none
    real(8), intent(inout) :: svalue
    integer, intent(in)    :: ngroup
    UNUSED_VARIABLE(svalue)
    UNUSED_VARIABLE(ngroup)
    !NOP! ABORT_MESSAGE(ngroup,"comm_get_min_double")
  end subroutine

  subroutine comm_get_min_array1d_double(invalue, outvalue, N, ngroup)
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    !NOP! ABORT_MESSAGE(ngroup,"comm_get_min_array1d_double")
    outvalue = invalue
  end subroutine

  subroutine comm_get_max_integer(svalue, ngroup)
    implicit none
    integer, intent(inout) :: svalue
    integer, intent(in)    :: ngroup
    UNUSED_VARIABLE(svalue)
    UNUSED_VARIABLE(ngroup)
    !NOP! ABORT_MESSAGE(ngroup,"comm_get_max_integer")
  end subroutine

  subroutine comm_get_max_array1d_double(invalue, outvalue, N, ngroup)
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    !NOP! ABORT_MESSAGE(ngroup,"comm_get_max_array1d_double")
    outvalue = invalue
  end subroutine

  subroutine comm_get_maxloc(invalue, outvalue, ngroup)
    implicit none
    type(comm_maxloc_type), intent(in)  :: invalue
    type(comm_maxloc_type), intent(out) :: outvalue
    integer, intent(in)                 :: ngroup
    UNUSED_VARIABLE(ngroup)
    !NOP! ABORT_MESSAGE(ngroup,"comm_get_maxloc")
    outvalue = invalue
  end subroutine

  subroutine comm_logical_and_scalar(invalue, outvalue, ngroup)
    implicit none
    logical, intent(in)  :: invalue
    logical, intent(out) :: outvalue
    integer, intent(in)  :: ngroup
    UNUSED_VARIABLE(ngroup)
    !NOP! ABORT_MESSAGE(ngroup,"comm_logical_and_scalar")
    outvalue = invalue
  end subroutine

  subroutine comm_ip_logical_or_array1d(values, ngroup)
    implicit none
    logical, intent(inout) :: values(:)
    integer, intent(in)    :: ngroup
    UNUSED_VARIABLE(values)
    UNUSED_VARIABLE(ngroup)
    !NOP! ABORT_MESSAGE(ngroup,"comm_ip_logical_or_array")
  end subroutine


  subroutine get_rank(comm, npid, nprocs)
    implicit none
    integer, intent(in)  :: comm
    integer, intent(out) :: npid, nprocs
    UNUSED_VARIABLE(comm)
    npid   = 0
    nprocs = 1
  end subroutine

  subroutine error_check(errcode)
    implicit none
    integer, intent(in) :: errcode
    UNUSED_VARIABLE(errcode)
    call abort_show_message('error_check')
  end subroutine

  subroutine abort_show_message(msg)
    implicit none
    character(*), intent(in) :: msg
    print '(A,A)', msg, ': this subroutine must not called (it takes MPI)'
#ifdef __INTEL_COMPILER
    call tracebackqq
#endif
    stop
  end subroutine
  
end module
