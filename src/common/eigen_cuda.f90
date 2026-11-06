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
! Host side of the fused CUDA eigen-energy sums, built with USE_CUDA.
! Kept out of total_energy so its OpenACC kernels compile unchanged.
module eigen_cuda
  implicit none
  private
  public :: eigen_sums_cuda

  interface
    integer(c_int) function eigen_sums(t, h, s, v0, v1, npl, nz, nspin, nio, nik, pz, ps, po, pk, vz, &
                                       sums, stream) bind(C, name="eigen_sums")
      use iso_c_binding, only: c_int, c_long_long, c_double, c_double_complex
      ! Assumed-size dummies take the address without target, which would make the caller's
      ! OpenACC kernels assume aliasing.
      complex(c_double_complex), intent(in) :: t(*), h(*), s(*)
      real(c_double), intent(in) :: v0(*), v1(*)
      integer(c_int), value :: npl, nz, nspin, nio, nik
      integer(c_long_long), value :: pz, ps, po, pk, vz, stream
      real(c_double), intent(out) :: sums(*)
    end function eigen_sums
  end interface

  ! Status codes returned by the kernels; keep in sync with cuda_stencil.h.
  integer, parameter :: st_ok = 0, st_device = 2, st_memory = 3
  integer, parameter :: st_grid = 100    ! decided here: an x-y plane is not one contiguous run

contains

  ! Sets wrk1, E_local(1) and E_local(2) as calc_eigen_energy's complex loops do and returns .true.,
  ! or returns .false. (all untouched) so those loops run instead.
  logical function eigen_sums_cuda(tpsi,htpsi,ttpsi,V_local,system,info,mg,im,wrk1,E_local) result(done)
    use structures, only: s_orbital, s_scalar, s_dft_system, s_parallel_info, s_rgrid
    use iso_c_binding, only: c_long_long
    use openacc, only: acc_get_cuda_stream, acc_async_sync
    implicit none
    type(s_orbital),intent(in) :: tpsi,htpsi,ttpsi
    type(s_dft_system),intent(in) :: system
    type(s_scalar),intent(in) :: V_local(system%nspin)
    type(s_parallel_info),intent(in) :: info
    type(s_rgrid),intent(in) :: mg
    integer,intent(in) :: im
    real(8),intent(inout) :: wrk1(system%nspin,system%no,system%nk)
    real(8),intent(out) :: E_local(2)
    integer,save :: status = st_ok
    logical,save :: reported = .false.
    real(8),allocatable :: sums(:,:)
    integer :: ik,io,ispin,i,nspin,nio,nik
    integer(c_long_long) :: pz,ps,po,pk,stream

    done = .false.
    ! A failed first call keeps the OpenACC loops for the rest of the run.
    if (status /= st_ok) return
    nspin = system%nspin
    if (.not. (all(mg%is_array(1:2) == mg%is(1:2)) .and. mg%ie_array(1) == mg%ie(1) .and. &
               all(lbound(V_local(1)%f) == mg%is) .and. all(ubound(V_local(1)%f) == mg%ie))) then
      call report(st_grid, reported)
      return
    end if
    nio = info%io_e - info%io_s + 1
    nik = info%ik_e - info%ik_s + 1
    pz = size(tpsi%zwf, 1) * size(tpsi%zwf, 2)
    ps = pz * size(tpsi%zwf, 3)
    po = ps * size(tpsi%zwf, 4)
    pk = po * size(tpsi%zwf, 5)
    allocate(sums(3, nspin*nio*nik))
    stream = transfer(acc_get_cuda_stream(acc_async_sync), stream)
    status = eigen_sums(tpsi%zwf(mg%is(1),mg%is(2),mg%is(3),1,info%io_s,info%ik_s,im), &
                        htpsi%zwf(mg%is(1),mg%is(2),mg%is(3),1,info%io_s,info%ik_s,im), &
                        ttpsi%zwf(mg%is(1),mg%is(2),mg%is(3),1,info%io_s,info%ik_s,im), &
                        V_local(1)%f(mg%is(1),mg%is(2),mg%is(3)), V_local(nspin)%f(mg%is(1),mg%is(2),mg%is(3)), &
                        mg%num(1)*mg%num(2), mg%num(3), nspin, nio, nik, pz, ps, po, pk, &
                        int(mg%num(1)*mg%num(2), c_long_long), sums, stream)
    call report(status, reported)
    if (status /= st_ok) return
    wrk1 = 0d0
    E_local = 0d0
    ! The energies accumulate in the OpenACC loops' order: spin, then k, then orbital.
    do ispin = 1, nspin
    do ik = info%ik_s, info%ik_e
    do io = info%io_s, info%io_e
      i = ispin + nspin*((io - info%io_s) + nio*(ik - info%ik_s))
      wrk1(ispin,io,ik) = sums(1,i) * system%Hvol
      E_local(1) = E_local(1) + system%rocc(io,ik,ispin)*system%wtk(ik) * sums(2,i) * system%Hvol
      E_local(2) = E_local(2) + system%rocc(io,ik,ispin)*system%wtk(ik) * system%hvol * sums(3,i)
    end do
    end do
    end do
    done = .true.
  end function eigen_sums_cuda

  ! One line on the root rank, the first time the sums are reached.
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
      write(*,'(a)') ' CUDA eigen-energy sums: active'
      return
    case (st_device); why = 'GPU older than sm_80'
    case (st_memory); why = 'orbitals not in GPU-accessible memory'
    case (st_grid);   why = 'x-y plane divided across ranks'
    case default;     why = 'unknown status'
    end select
    write(*,'(a,a)') ' CUDA eigen-energy sums: not used, OpenACC instead: ', trim(why)
  end subroutine report

end module eigen_cuda
