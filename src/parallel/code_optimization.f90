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
! code optimization routines, variables, and others.

#include "config.h"

module code_optimization
  implicit none

  ! current/force OpenMP parallelization mode
  !   .TRUE.  : k,orbital parallel
  !   .FALSE. : rgrid parallel
  logical :: current_omp_mode
  logical :: force_omp_mode

  ! can we call optimized stencil code in hamiltonian?
  logical :: optimized_stencil_is_callable

  ! is a stencil computation in hamiltonian parallelized by OpenMP?
  logical :: stencil_is_parallelized_by_omp

  ! these tables store results of modulo
  integer,allocatable :: modx(:), mody(:), modz(:)

contains
  subroutine switch_stencil_optimization(nlens)
    use salmon_global, only: yn_want_stencil_hand_vectorization
    implicit none
    integer,intent(in) :: nlens(3) ! a length of x,y,z directions
#ifdef USE_OPT_EXPLICIT_VECTORIZATION
    ! unit-stride direction is should multiple of 4
    optimized_stencil_is_callable = (mod(nlens(1), 4) == 0) &
    &    .and. (yn_want_stencil_hand_vectorization == 'y')
#else
    optimized_stencil_is_callable = .false.
#endif
  end subroutine

  subroutine switch_openmp_parallelization(nlens)
    use salmon_global, only: stencil_openmp_mode,current_openmp_mode,force_openmp_mode
    implicit none
    integer,intent(in) :: nlens(3) ! a length of x,y,z directions

    select case(stencil_openmp_mode)
    case('orbital')
      stencil_is_parallelized_by_omp = .false.
    case('rgrid')
      stencil_is_parallelized_by_omp = .true.
    case default
      stencil_is_parallelized_by_omp = (nlens(1)*nlens(2)*nlens(3) >= 32*32*32)
    end select

    select case(current_openmp_mode)
    case('orbital')
      current_omp_mode = .true.
    case('rgrid')
      current_omp_mode = .false.
    case default
      current_omp_mode = .not. stencil_is_parallelized_by_omp
    end select

    select case(force_openmp_mode)
    case('orbital')
      force_omp_mode = .true.
    case('rgrid')
      force_omp_mode = .false.
    case default
      force_omp_mode = .not. stencil_is_parallelized_by_omp
    end select
  end subroutine

  subroutine set_modulo_tables(nlens)
    implicit none
    integer,intent(in) :: nlens(3) ! a length of x,y,z directions
    integer :: ix,iy,iz

    if (allocated(modx)) deallocate(modx)
    if (allocated(mody)) deallocate(mody)
    if (allocated(modz)) deallocate(modz)

    allocate(modx(0:nlens(1)*3-1))
    allocate(mody(0:nlens(2)*3-1))
    allocate(modz(0:nlens(3)*3-1))

    do ix=0,nlens(1)*3-1
      modx(ix) = mod(ix, nlens(1))
    end do
    do iy=0,nlens(2)*3-1
      mody(iy) = mod(iy, nlens(2))
    end do
    do iz=0,nlens(3)*3-1
      modz(iz) = mod(iz, nlens(3))
    end do
  end subroutine

  subroutine optimization_log(info)
    use structures, only: s_parallel_info
    use parallelization, only: is_distributed_parallel, get_nthreads
    implicit none
    type(s_parallel_info), intent(in) :: info
    print *, '======'
    if (is_distributed_parallel()) then
      print *, 'MPI distribution:'
      print *, '  nproc_k     :', info%npk
      print *, '  nproc_ob    :', info%nporbital
      print *, '  nproc_rgrid :', info%nprgrid
    end if
    print *, 'OpenMP parallelization:'
    print *, '  number of threads :', get_nthreads()
    print *, 'hpsi stencil:'
    print *, '  enables hand-coding vectorization :', optimized_stencil_is_callable
    print *, '  enables openmp parallelization    :', stencil_is_parallelized_by_omp
    print *, 'current:'
    print *, '  k,orbital parallelized of openmp  :', current_omp_mode
    print *, 'force:'
    print *, '  k,orbital parallelized of openmp  :', force_omp_mode
  end subroutine
end module
