!
!  Copyright 2026 SALMON developers
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
! Host side of the cuFFT periodic Poisson solver, built with USE_CUDA.
! Kept out of poisson_periodic so its OpenACC kernels compile unchanged.
module poisson_cuda
  implicit none
  private
  public :: poisson_fft_cuda

  interface
    integer(c_int) function hartree_fft(rho, vh, coef, zrhog, work, nx, ny, nz, stream) &
                                        bind(C, name="hartree_fft")
      use iso_c_binding, only: c_int, c_long_long, c_double, c_double_complex
      ! Assumed-size dummies take the address without target, which would make the caller's
      ! OpenACC kernels assume aliasing.
      real(c_double), intent(in) :: rho(*), coef(*)
      real(c_double), intent(inout) :: vh(*)
      complex(c_double_complex), intent(inout) :: zrhog(*), work(*)
      integer(c_int), value :: nx, ny, nz
      integer(c_long_long), value :: stream
    end function hartree_fft
  end interface

  ! Status codes returned by the solver; keep in sync with cuda_stencil.h.
  integer, parameter :: st_ok = 0, st_device = 2, st_memory = 3, st_nofit = 4
  integer, parameter :: st_grid = 100    ! decided here: the grid is divided across ranks

contains

  ! Sets Vh and poisson%zrhoG_ele as poisson_ft does and returns .true., or returns .false.
  ! (both untouched) so poisson_ft computes them instead.
  logical function poisson_fft_cuda(lg,mg,fg,rho,Vh,poisson) result(done)
    use structures, only: s_rgrid, s_reciprocal_grid, s_scalar, s_poisson
    use iso_c_binding, only: c_long_long
    use openacc, only: acc_get_cuda_stream, acc_async_sync
    implicit none
    type(s_rgrid),intent(in) :: lg,mg
    type(s_reciprocal_grid),intent(in) :: fg
    type(s_scalar),intent(in) :: rho
    type(s_scalar) :: Vh
    type(s_poisson) :: poisson
    integer,save :: status = st_ok
    logical,save :: reported = .false.
    integer(c_long_long) :: stream

    done = .false.
    ! A failed first call (no plan, unsupported GPU) keeps poisson_ft for the rest of the run.
    if (status /= st_ok) return
    if (.not. (all(mg%is == lg%is) .and. all(mg%ie == lg%ie))) then
      call report(st_grid, reported)
      return
    end if
    stream = transfer(acc_get_cuda_stream(acc_async_sync), stream)
    status = hartree_fft(rho%f, Vh%f, fg%coef, poisson%zrhoG_ele, poisson%ff1z, &
                         lg%num(1), lg%num(2), lg%num(3), stream)
    done = (status == st_ok)
    call report(status, reported)
  end function poisson_fft_cuda

  ! One line on the root rank, the first time the solver is reached.
  subroutine report(status, reported)
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    implicit none
    integer,intent(in) :: status
    logical,intent(inout) :: reported
    character(64) :: why
    if (reported) return
    reported = .true.
    if (.not. comm_is_root(nproc_id_global)) return
    select case (status)
    case (st_ok)
      write(*,'(a)') ' CUDA Hartree FFT: active'
      return
    case (st_device); why = 'GPU older than sm_80'
    case (st_memory); why = 'arrays not in GPU-accessible memory'
    case (st_nofit);  why = 'cuFFT plan could not be created'
    case (st_grid);   why = 'grid divided across ranks'
    case default;     why = 'unknown status'
    end select
    write(*,'(a,a)') ' CUDA Hartree FFT: not used, OpenACC instead: ', trim(why)
  end subroutine report

end module poisson_cuda
