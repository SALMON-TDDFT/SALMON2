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
module parallelization
  implicit none

  !!! multi process
  ! global
  integer, public :: nproc_group_global
  integer, public :: nproc_id_global
  integer, public :: nproc_size_global

  ! TDKS eq.
  integer, public :: nproc_group_tdks
  integer, public :: nproc_id_tdks
  integer, public :: nproc_size_tdks

  ! call once
  public :: setup_parallel
  public :: end_parallel

  ! util
  public :: get_thread_id
  public :: get_nthreads
  public :: is_distributed_parallel

contains
  subroutine setup_parallel
    use communication
    implicit none
    call comm_init
    call comm_get_globalinfo(nproc_group_global, nproc_id_global, nproc_size_global)
  end subroutine

  subroutine end_parallel
    use communication
    implicit none
    call comm_finalize
  end subroutine

  function get_thread_id() result(nid)
#ifdef _OPENMP
    use omp_lib
#endif
    implicit none
    integer :: nid
#ifdef _OPENMP
    nid = omp_get_thread_num()
#else
    nid = 0
#endif
  end function

  function get_nthreads() result(nsize)
#ifdef _OPENMP
    use omp_lib
#endif
    implicit none
    integer :: nsize
#ifdef _OPENMP
    nsize = omp_get_max_threads()
#else
    nsize = 1
#endif
  end function

  function is_distributed_parallel() result(ret)
    implicit none
    logical :: ret
    ret = (nproc_size_global > 1)
  end function
end module
