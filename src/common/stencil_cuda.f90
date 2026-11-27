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
! Host side of the CUDA stencil kernels, built with USE_CUDA.
! Kept out of the callers so their OpenACC kernels compile unchanged.
module stencil_cuda
  implicit none
  private
  public :: hpsi_stencil_cuda, current_stencil_cuda

  interface
    integer(c_int) function hpsi_stencil(tp, hp, vp, nxi, nyi, nzi, norb, sx2, sx3, sx5, vx2, vx3, &
                                         lap0, lapt, nabt, stream, cfg) bind(C, name="hpsi_stencil")
      use iso_c_binding, only: c_int, c_long_long, c_double, c_double_complex
      ! Assumed-size dummies take the address without target, which would make the callers'
      ! OpenACC kernels assume aliasing.
      complex(c_double_complex), intent(in) :: tp(*)
      complex(c_double_complex), intent(inout) :: hp(*)
      real(c_double), intent(in) :: vp(*)
      integer(c_int), value :: nxi, nyi, nzi, norb
      integer(c_long_long), value :: sx2, sx3, sx5, vx2, vx3, stream
      real(c_double), value :: lap0
      real(c_double), intent(in) :: lapt(12), nabt(12)
      integer(c_int), intent(out) :: cfg(2)
    end function hpsi_stencil

    integer(c_int) function current_stencil(tp, nx, ny, nz, norb, sx2, sx3, sx5, nab, kac, bt, &
                                            rocc, wtk, out3, stream, cfg) bind(C, name="current_stencil")
      use iso_c_binding, only: c_int, c_long_long, c_double, c_double_complex
      complex(c_double_complex), intent(in) :: tp(*)
      integer(c_int), value :: nx, ny, nz, norb
      integer(c_long_long), value :: sx2, sx3, sx5, stream
      real(c_double), intent(in) :: nab(*), kac(*), bt(*), rocc(*)
      real(c_double), value :: wtk
      real(c_double), intent(inout) :: out3(*)
      integer(c_int), intent(out) :: cfg(2)
    end function current_stencil
  end interface

  ! Status codes returned by the kernels; keep in sync with cuda_stencil.h.
  integer, parameter :: st_ok = 0, st_shape = 1, st_device = 2, st_memory = 3, st_nofit = 4
  integer, parameter :: st_grid = 100    ! decided here: divided or non-periodic grid

contains

  ! Applies the orthogonal stencil plus V to all orbitals and spins and returns .true.,
  ! or returns .false. so the caller's OpenACC kernels compute htpsi instead.
  logical function hpsi_stencil_cuda(io_s,io_e,Nspin,is_array,ie_array,is,ie,idx,idy,idz,igs,ige, &
                                     tpsi,htpsi,V_local,lap0,lapt,nabt) result(done)
    use structures, only: s_scalar
    use iso_c_binding, only: c_long_long
    use openacc, only: acc_get_cuda_stream, acc_async_sync
    implicit none
    integer,intent(in) :: io_s,io_e,Nspin
    integer,intent(in) :: is_array(3),ie_array(3),is(3),ie(3)
    integer,intent(in) :: idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)
    integer,intent(in) :: igs(3),ige(3)
    complex(8),intent(in)  :: tpsi (is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3),1:Nspin,io_s:io_e)
    complex(8),intent(out) :: htpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3),1:Nspin,io_s:io_e)
    type(s_scalar),intent(in) :: V_local(1:Nspin)
    real(8),intent(in) :: lap0,lapt(12),nabt(12)
    logical,save :: reported = .false.
    integer :: ispin,ierr,cfg(2)
    integer(c_long_long) :: sx2,sx3,sx5,vx2,vx3,stream

    done = .false.
    cfg = 0
    if (.not. (all(igs == is) .and. all(ige == ie) .and. periodic(is,ie,idx,idy,idz))) then
      call report('hpsi stencil', st_grid, cfg, reported)
      return
    end if
    sx2 = ie_array(1) - is_array(1) + 1
    sx3 = sx2 * (ie_array(2) - is_array(2) + 1)
    sx5 = sx3 * (ie_array(3) - is_array(3) + 1) * Nspin
    stream = transfer(acc_get_cuda_stream(acc_async_sync), stream)
    do ispin = 1, Nspin
      vx2 = size(V_local(ispin)%f, 1)
      vx3 = vx2 * size(V_local(ispin)%f, 2)
      ierr = hpsi_stencil(tpsi(is(1),is(2),is(3),ispin,io_s), htpsi(is(1),is(2),is(3),ispin,io_s), &
                          V_local(ispin)%f(is(1),is(2),is(3)), &
                          ie(1)-is(1)+1, ie(2)-is(2)+1, ie(3)-is(3)+1, io_e-io_s+1, &
                          sx2, sx3, sx5, vx2, vx3, lap0, lapt, nabt, stream, cfg)
      if (ierr /= st_ok) exit
    end do
    done = (ierr == st_ok)
    call report('hpsi stencil', ierr, cfg, reported)
  end function hpsi_stencil_cuda

  ! Computes calc_current's stencil part for one spin and band group into jx,jy,jz and
  ! returns .true., or returns .false. (jx,jy,jz untouched) so the OpenACC kernel runs instead.
  logical function current_stencil_cuda(system,mg,stencil,info,psi,BT,ispin,im,jx,jy,jz) result(done)
    use structures, only: s_dft_system, s_rgrid, s_stencil, s_parallel_info, s_orbital
    use iso_c_binding, only: c_long_long
    use openacc, only: acc_get_cuda_stream, acc_async_sync
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_stencil),intent(in) :: stencil
    type(s_parallel_info),intent(in) :: info
    type(s_orbital) :: psi
    real(8),intent(in) :: BT(3,3)
    integer,intent(in) :: ispin,im
    real(8),intent(inout) :: jx,jy,jz
    logical,save :: reported = .false.
    integer :: ik,ierr,cfg(2)
    integer(c_long_long) :: sx2,sx3,sx5,stream
    real(8) :: out(3),kac(3)

    done = .false.
    cfg = 0
    if (info%if_divide_rspace .or. .not. periodic(mg%is,mg%ie,mg%idx,mg%idy,mg%idz)) then
      call report('current stencil', st_grid, cfg, reported)
      return
    end if
    sx2 = mg%ie_array(1) - mg%is_array(1) + 1
    sx3 = sx2 * (mg%ie_array(2) - mg%is_array(2) + 1)
    sx5 = sx3 * (mg%ie_array(3) - mg%is_array(3) + 1) * system%nspin
    stream = transfer(acc_get_cuda_stream(acc_async_sync), stream)
    out = 0d0
    ierr = st_ok
    do ik = info%ik_s, info%ik_e
      kac = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
      ierr = current_stencil(psi%zwf(mg%is(1),mg%is(2),mg%is(3),ispin,info%io_s,ik,im), &
                             mg%num(1), mg%num(2), mg%num(3), info%io_e - info%io_s + 1, &
                             sx2, sx3, sx5, stencil%coef_nab, kac, BT, &
                             system%rocc(info%io_s,ik,ispin), system%wtk(ik), out, stream, cfg)
      if (ierr /= st_ok) exit
    end do
    done = (ierr == st_ok)
    if (done) then
      jx = out(1)
      jy = out(2)
      jz = out(3)
    end if
    call report('current stencil', ierr, cfg, reported)
  end function current_stencil_cuda

  ! The kernels wrap neighbours themselves: idx/idy/idz must be the periodic wrap of [is,ie].
  logical function periodic(is,ie,idx,idy,idz)
    implicit none
    integer,intent(in) :: is(3),ie(3)
    integer,intent(in) :: idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)
    integer :: d,i,n(3)
    n = ie - is + 1
    periodic = .true.
    do d = -4, 4
      do i = is(1), ie(1)
        periodic = periodic .and. idx(i+d) == is(1) + modulo(i+d-is(1), n(1))
      end do
      do i = is(2), ie(2)
        periodic = periodic .and. idy(i+d) == is(2) + modulo(i+d-is(2), n(2))
      end do
      do i = is(3), ie(3)
        periodic = periodic .and. idz(i+d) == is(3) + modulo(i+d-is(3), n(3))
      end do
    end do
  end function periodic

  ! One line on the root rank, the first time each kernel is reached.
  subroutine report(kernel, status, cfg, reported)
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    implicit none
    character(*),intent(in) :: kernel
    integer,intent(in) :: status, cfg(2)
    logical,intent(inout) :: reported
    character(64) :: why
    if (reported) return
    reported = .true.
    if (.not. comm_is_root(nproc_id_global)) return
    select case (status)
    case (st_ok)
      write(*,'(a,a,a,i0,a,i0,a)') ' CUDA ', kernel, ': active (tile width ', cfg(1), ', lookahead ', cfg(2), ')'
      return
    case (st_shape);  why = 'grid too small'
    case (st_device); why = 'GPU older than sm_80'
    case (st_memory); why = 'orbitals not in GPU-accessible memory'
    case (st_nofit);  why = 'no configuration fits this GPU'
    case (st_grid);   why = 'grid divided across ranks or not periodic'
    case default;     why = 'unknown status'
    end select
    write(*,'(a,a,a,a)') ' CUDA ', kernel, ': not used, OpenACC instead: ', trim(why)
  end subroutine report

end module stencil_cuda
