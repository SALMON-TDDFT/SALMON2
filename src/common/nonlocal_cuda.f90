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
! Host side of the tiled CUDA back-projection, built with USE_CUDA.
! Kept out of zpseudo so its OpenACC kernels compile unchanged.
module nonlocal_cuda
  implicit none
  private
  public :: back_projection_cuda

  interface
    integer(c_int) function back_projection_plan(max_vi, kmax, v2nlma, k2ilma, k2j, z2) &
                                                 bind(C, name="back_projection_plan")
      use iso_c_binding, only: c_int, c_long_long
      integer(c_int), value :: max_vi, kmax
      integer(c_int), intent(in) :: v2nlma(*), k2ilma(*), k2j(*)
      integer(c_long_long), value :: z2
    end function back_projection_plan

    integer(c_int) function back_projection(zekr, ubox, htpsi, v2j, nspin, nio, nik, nim, ub2, ub3, ub4, ub5, z3, &
                                            lb1, lb2, lb3, s2, s3, s4, s5, s6, s7, zoff0, hoff0, stream) &
                                            bind(C, name="back_projection")
      use iso_c_binding, only: c_int, c_long_long, c_double_complex
      ! Assumed-size dummies take the address without target, which would make the caller's
      ! OpenACC kernels assume aliasing.
      complex(c_double_complex), intent(in) :: zekr(*), ubox(*)
      complex(c_double_complex), intent(inout) :: htpsi(*)
      integer(c_int), intent(in) :: v2j(*)
      integer(c_int), value :: nspin, nio, nik, nim, lb1, lb2, lb3
      integer(c_long_long), value :: ub2, ub3, ub4, ub5, z3, s2, s3, s4, s5, s6, s7, zoff0, hoff0, stream
    end function back_projection
  end interface

  ! Status codes returned by the kernel; keep in sync with cuda_stencil.h.
  integer, parameter :: st_ok = 0, st_shape = 1, st_device = 2, st_memory = 3, st_nofit = 4

contains

  ! Adds phase 2 of zpseudo, htpsi += zekr_uV ubox, and returns .true., or returns .false.
  ! (htpsi untouched) so the caller runs its own kernel instead.
  logical function back_projection_cuda(ppg,ubox,htpsi,Nspin,io_s,io_e,ik_s,ik_e,im_s,im_e) result(done)
    use structures, only: s_pp_grid, s_orbital
    use iso_c_binding, only: c_long_long
    use openacc, only: acc_get_cuda_stream, acc_async_sync
    implicit none
    type(s_pp_grid),intent(in) :: ppg
    integer,intent(in) :: Nspin,io_s,io_e,ik_s,ik_e,im_s,im_e
    complex(8),intent(in) :: ubox(ppg%Nlma,Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e)
    type(s_orbital) :: htpsi
    integer,save :: plan_version = -1, plan_status = -1
    integer(c_long_long),save :: plan_z2 = -1
    logical,save :: reported = .false.
    integer :: ierr
    integer(c_long_long) :: z2,z3,ub2,ub3,ub4,ub5,s(7),stream

    z2 = size(ppg%zekr_uV, 1)
    if (ppg%vi_tables_version /= plan_version .or. z2 /= plan_z2) then
      plan_status = back_projection_plan(ppg%max_vi, size(ppg%k2ilma, 2), ppg%v2nlma, ppg%k2ilma, ppg%k2j, z2)
      plan_version = ppg%vi_tables_version
      plan_z2 = z2
    end if
    ierr = plan_status
    if (ierr == st_ok) then
      z3 = z2 * size(ppg%zekr_uV, 2)
      ub2 = size(ubox, 1)
      ub3 = ub2 * Nspin
      ub4 = ub3 * (io_e - io_s + 1)
      ub5 = ub4 * (ik_e - ik_s + 1)
      s(1) = 1
      s(2) = size(htpsi%zwf, 1)
      s(3) = s(2) * size(htpsi%zwf, 2)
      s(4) = s(3) * size(htpsi%zwf, 3)
      s(5) = s(4) * size(htpsi%zwf, 4)
      s(6) = s(5) * size(htpsi%zwf, 5)
      s(7) = s(6) * size(htpsi%zwf, 6)
      stream = transfer(acc_get_cuda_stream(acc_async_sync), stream)
      ierr = back_projection(ppg%zekr_uV, ubox, htpsi%zwf, ppg%v2j, Nspin, io_e - io_s + 1, ik_e - ik_s + 1, &
                             im_e - im_s + 1, ub2, ub3, ub4, ub5, z3, &
                             lbound(htpsi%zwf, 1), lbound(htpsi%zwf, 2), lbound(htpsi%zwf, 3), &
                             s(2), s(3), s(4), s(5), s(6), s(7), (ik_s - lbound(ppg%zekr_uV, 3)) * z3, &
                             (ik_s - lbound(htpsi%zwf, 6)) * s(6) + (im_s - lbound(htpsi%zwf, 7)) * s(7), stream)
    end if
    done = (ierr == st_ok)
    call report(ierr, reported)
  end function back_projection_cuda

  ! One line on the root rank, the first time the kernel is reached.
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
      write(*,'(a)') ' CUDA back-projection: active'
      return
    case (st_shape);  why = 'projector tables too large for the tiles'
    case (st_device); why = 'GPU older than sm_80'
    case (st_memory); why = 'arrays not in GPU-accessible memory'
    case (st_nofit);  why = 'shared memory does not fit this GPU'
    case default;     why = 'unknown status'
    end select
    write(*,'(a,a)') ' CUDA back-projection: not used, OpenACC instead: ', trim(why)
  end subroutine report

end module nonlocal_cuda
