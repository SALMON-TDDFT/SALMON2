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
    implicit none
    integer,intent(in) :: nlens(3) ! a length of x,y,z directions
#ifdef SALMON_EXPLICIT_VECTORIZATION
    ! unit-stride direction is should multiple of 4
    optimized_stencil_is_callable = (mod(nlens(1), 4) == 0)
#else
    optimized_stencil_is_callable = .false.
#endif
  end subroutine

  subroutine switch_openmp_parallelization(nlens)
    implicit none
    integer,intent(in) :: nlens(3) ! a length of x,y,z directions
    stencil_is_parallelized_by_omp = (nlens(1)*nlens(2)*nlens(3) >= 32*32*32)
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

  subroutine optimization_log
    implicit none
    print *, '========== code optimization log =========='
    print *, 'hpsi stencil:'
    print *, '  enables hand-coding vectorization :', optimized_stencil_is_callable
    print *, '  enables openmp parallelization    :', stencil_is_parallelized_by_omp
  end subroutine
end module
