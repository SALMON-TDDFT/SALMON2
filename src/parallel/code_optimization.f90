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
! code optimization routines, variables, and others.

#include "config.h"

module code_optimization
  implicit none

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
    use salmon_global, only: yn_want_stencil_openmp_parallelization, &
                             yn_force_stencil_openmp_parallelization, &
                             yn_force_stencil_sequential_computation
    implicit none
    integer,intent(in) :: nlens(3) ! a length of x,y,z directions
    logical :: ret
    ret = (yn_want_stencil_openmp_parallelization  == 'y') .and. (nlens(1)*nlens(2)*nlens(3) >= 32*32*32)

    if (yn_force_stencil_openmp_parallelization == 'y') ret = .true.  ! forcible option: use OpenMP parallelized code
    if (yn_force_stencil_sequential_computation == 'y') ret = .false. ! forcible option: use sequential code

    stencil_is_parallelized_by_omp = ret
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

  subroutine optimization_log(nproc_k, nproc_ob, nproc_domain_orbital, nproc_domain_general)
    use salmon_parallel, only: is_distributed_parallel, get_nthreads
    implicit none
    integer, intent(in) :: nproc_k, nproc_ob
    integer, intent(in) :: nproc_domain_orbital(3), nproc_domain_general(3)
    print *, '========== code optimization log =========='
    if (is_distributed_parallel()) then
      print *, 'MPI distribution:'
      print *, '  nproc_ob       :', nproc_ob
      print *, '  nproc_domain_orbital   :', nproc_domain_orbital
      print *, '  nproc_domain_general :', nproc_domain_general
    end if
    print *, 'OpenMP parallelization:'
    print *, '  number of threads :', get_nthreads()
    print *, 'hpsi stencil:'
    print *, '  enables hand-coding vectorization :', optimized_stencil_is_callable
    print *, '  enables openmp parallelization    :', stencil_is_parallelized_by_omp
  end subroutine
end module
