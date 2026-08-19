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

#include "config.h"

module nonlocal_potential
  use nvtx_wrapper
  implicit none

! WARNING: We must not call these except for hpsi routine.

contains

subroutine dpseudo(tpsi,htpsi,info,nspin,ppg)
  use structures
  use communication, only: comm_summation
  use timer
  implicit none
  integer,intent(in) :: nspin
  type(s_parallel_info),intent(in) :: info
  type(s_pp_grid),intent(in) :: ppg
  type(s_orbital),intent(in) :: tpsi
  type(s_orbital) :: htpsi
  !
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb
  integer :: ilma,ia,j,ix,iy,iz,Nlma,ilocal
  real(8) :: uVpsi,wrk
  real(8),allocatable :: uVpsibox (:,:,:,:,:)
  real(8),allocatable :: uVpsibox2(:,:,:,:,:)

  call timer_begin(LOG_UHPSI_PSEUDO)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  norb = Nspin* info%numo * info%numk * info%numm

  Nlma = ppg%Nlma

  if(info%if_divide_rspace) then

    allocate(uVpsibox (Nlma,Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e))
    allocate(uVpsibox2(Nlma,Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e))

!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilma,ia,uVpsi,j,ix,iy,iz)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi = 0.d0
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          uVpsi = uVpsi + ppg%uV(j,ilma) * tpsi%rwf(ix,iy,iz,ispin,io,ik,im)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(ilma)
        uVpsibox(ilma,ispin,io,ik,im) = uVpsi
      end do

    end do
    end do
    end do
    end do
!$omp end parallel do

    call timer_end(LOG_UHPSI_PSEUDO)

    call timer_begin(LOG_UHPSI_PSEUDO_COMM)
    call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,info%icomm_r)
    call timer_end(LOG_UHPSI_PSEUDO_COMM)

    call timer_begin(LOG_UHPSI_PSEUDO)

!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz,wrk)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do ilocal=1,ppg%ilocal_nlma
        ilma=ppg%ilocal_nlma2ilma(ilocal)
        ia  =ppg%ilocal_nlma2ia  (ilocal)
        uVpsi = uVpsibox2(ilma,ispin,io,ik,im)
!OCL norecurrence
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          wrk = uVpsi * ppg%uV(j,ilma)
          htpsi%rwf(ix,iy,iz,ispin,io,ik,im) = htpsi%rwf(ix,iy,iz,ispin,io,ik,im) + wrk
        end do
      end do

    end do
    end do
    end do
    end do
!$omp end parallel do

    deallocate(uVpsibox,uVpsibox2)

  else

#ifdef USE_OPENACC
!$acc parallel loop collapse(4) private(im,ik,io,ispin,ilma,ia,uVpsi,j,ix,iy,iz,wrk)
#else
!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilma,ia,uVpsi,j,ix,iy,iz,wrk)
#endif
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi = 0.d0
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          uVpsi = uVpsi + ppg%uV(j,ilma) * tpsi%rwf(ix,iy,iz,ispin,io,ik,im)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(ilma)
!OCL norecurrence
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          wrk = uVpsi * ppg%uV(j,ilma)
          htpsi%rwf(ix,iy,iz,ispin,io,ik,im) = htpsi%rwf(ix,iy,iz,ispin,io,ik,im) + wrk
        end do
      end do

    end do
    end do
    end do
    end do
#ifdef USE_OPENACC
!$acc end parallel
#else
!$omp end parallel do
#endif

  end if

  call timer_end(LOG_UHPSI_PSEUDO)

  return
end subroutine dpseudo

!-----------------------------------------------------------------------------------------------------------------------------------

subroutine zpseudo(tpsi,htpsi,info,nspin,ppg)
  use structures
  use timer
  use iso_c_binding
#if defined(USE_OPENACC) && defined(USE_GEMM)
  use cublas
  use openacc
  use cudafor
  use mpi, only: MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_IN_PLACE
#ifdef USE_NCCL
  use nccl
  use mpi, only: MPI_CHARACTER
#endif
#endif
  implicit none
  intrinsic :: aimag
  integer,intent(in) :: nspin
  type(s_parallel_info),intent(in) :: info
  type(s_pp_grid),intent(in) :: ppg
  type(s_orbital),intent(in) :: tpsi
  type(s_orbital) :: htpsi
  !
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb
  integer :: ilma,ia,j,ix,iy,iz,Nlma,ilocal,vi,my_nlma,k
  complex(8) :: uVpsi,wrk
  complex(8),allocatable :: uVpsibox (:,:,:,:,:)
  complex(8),allocatable :: uVpsibox2(:,:,:,:,:)
#if defined(USE_OPENACC) && defined(USE_GEMM)
  complex(8), device, allocatable, save, target :: d_reduce(:,:,:,:,:)
  integer, save :: d_reduce_allocated = 0
#ifdef USE_NCCL
  ! NCCL reduces on the GPU; MPI_Allreduce does the arithmetic on the host,
  ! which costs ~30x on this buffer even though both take a device pointer.
  type(ncclComm), save :: nccl_comm
  integer(cuda_stream_kind), save :: nccl_stream
  logical, save :: nccl_ready = .false.
  type(ncclUniqueId) :: nccl_uid
  type(ncclResult) :: nccl_stat
  real(8), device, pointer :: d_reduce_re(:)
  integer :: nccl_rank, nccl_size, cuda_stat
#endif
#endif
#ifdef USE_OPENACC
  real(8),allocatable,save :: htpsi_zwf_r(:,:,:,:,:,:,:), htpsi_zwf_i(:,:,:,:,:,:,:)
  integer :: i1, i2, i3, i4, i5, i6, i7
  integer :: mps_max
  real(8) :: wrk_r
#endif
  complex(8) :: IMAGINARY_UNIT = (0, 1)
#if defined(USE_OPENACC) && defined(USE_CUDA)
  integer :: n,natom
  interface
    subroutine zpseudo_cuda(htpsi_zwf,n,im_s,im_e,ik_s,ik_e,io_s,io_e,Nspin,Nlma,ppg_nps,natom,mg_is_array_1,mg_ie_array_1,mg_is_array_2,mg_ie_array_2,mg_is_array_3,mg_ie_array_3,ppg_ia_tbl,ppg_mps,ppg_jxyz,ppg_zekr_uV,ppg_rinv_uvu,tpsi_zwf) bind(c)
      import
      ! Input integer
      integer(c_int), intent(in), value :: n
      integer(c_int), intent(in), value :: im_s
      integer(c_int), intent(in), value :: im_e
      integer(c_int), intent(in), value :: ik_s
      integer(c_int), intent(in), value :: ik_e
      integer(c_int), intent(in), value :: io_s
      integer(c_int), intent(in), value :: io_e
      integer(c_int), intent(in), value :: Nspin
      integer(c_int), intent(in), value :: Nlma
      integer(c_int), intent(in), value :: ppg_nps
      integer(c_int), intent(in), value :: natom
      integer(c_int), intent(in), value :: mg_is_array_1
      integer(c_int), intent(in), value :: mg_ie_array_1
      integer(c_int), intent(in), value :: mg_is_array_2
      integer(c_int), intent(in), value :: mg_ie_array_2
      integer(c_int), intent(in), value :: mg_is_array_3
      integer(c_int), intent(in), value :: mg_ie_array_3
      ! Output
      complex(c_double_complex), intent(inout) :: htpsi_zwf(mg_is_array_1:mg_ie_array_1, mg_is_array_2:mg_ie_array_2, mg_is_array_2:mg_ie_array_3, Nspin, io_e:io_s, ik_s:ik_e, im_s:im_e)
      ! Input
      integer(c_int)           , intent(in) :: ppg_ia_tbl(n * natom)
      integer(c_int)           , intent(in) :: ppg_mps(natom)
      integer(c_int)           , intent(in) :: ppg_jxyz(3, ppg_nps, natom)
      complex(c_double_complex), intent(in) :: ppg_zekr_uV(ppg_nps, ppg_Nlma, ik_s:ik_e)
      real(c_double)           , intent(in) :: ppg_rinv_uvu(n * natom)
      complex(c_double_complex), intent(in) :: tpsi_zwf(mg_is_array_1:mg_ie_array_1, mg_is_array_2:mg_ie_array_2,mg_is_array_2:mg_ie_array_3, Nspin, io_e:io_s, ik_s:ik_e, im_s:im_e)
    end subroutine zpseudo_cuda
  end interface
#endif
#if defined(USE_OPENACC) && defined(USE_GEMM)
  ! Batched-GEMM reformulation of phase 1 (projection); phase 2 (back-projection) is unchanged below.
  integer, parameter :: gemm_io_block = 64
  integer,save :: gemm_max_nproj = -1
  integer,allocatable,save :: gemm_nproj_atom(:), gemm_l2g(:,:)
  complex(8),allocatable,save :: gemm_zekr_packed(:,:,:)
  real(8),allocatable,save :: gemm_rinv_packed(:,:)
  complex(8),allocatable,save :: gemm_wf_packed(:,:,:), gemm_out_packed(:,:,:)
  type(cublasHandle),save :: gemm_handle
  integer :: gemm_ia, gemm_p, gemm_ilma, gemm_j, gemm_natom
  integer :: gemm_io_blk_s, gemm_this_block, gemm_io_local, gemm_stat
  integer :: mpi_ierr
#endif

#ifdef FORTRAN_COMPILER_HAS_2MB_ALIGNED_ALLOCATION
!dir$ attributes align : 2097152 :: uVpsibox, uVpsibox2
#endif

  call nvtxStartRange('zpseudo', __LINE__)
  call timer_begin(LOG_UHPSI_PSEUDO)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e

  Nlma = ppg%Nlma

#if defined(USE_OPENACC) && defined(USE_GEMM)
  ! GEMM path, both domain and non-domain decomposition: ppg%jxyz/mps are already local,
  ! so under domain decomposition the GEMM yields per-rank partial sums, reduced in phase 1.5.

  gemm_natom = size(ppg%mps)

  ! One-time setup: build zero-padded dense companion arrays, since ppg%zekr_uV/jxyz are only filled up to mps(ia) and reading beyond that is uninitialized.
  if (gemm_max_nproj < 0) then
    allocate(gemm_nproj_atom(gemm_natom))
    gemm_nproj_atom = 0
    do ilma=1,Nlma
      ia = ppg%ia_tbl(ilma)
      gemm_nproj_atom(ia) = gemm_nproj_atom(ia) + 1
    end do
    gemm_max_nproj = maxval(gemm_nproj_atom)

    allocate(gemm_l2g(gemm_max_nproj, gemm_natom))
    gemm_l2g = 0
    gemm_nproj_atom = 0   ! reused as a per-atom fill cursor below
    do ilma=1,Nlma
      ia = ppg%ia_tbl(ilma)
      gemm_nproj_atom(ia) = gemm_nproj_atom(ia) + 1
      gemm_l2g(gemm_nproj_atom(ia), ia) = ilma
    end do

    ! Kept 3D (no ik dimension): a slice through host_data/cublas triggers an nvfortran ICE.
    allocate(gemm_zekr_packed(ppg%nps, gemm_max_nproj, gemm_natom))
    allocate(gemm_rinv_packed(gemm_max_nproj, gemm_natom))
    gemm_zekr_packed = (0.d0,0.d0)
    gemm_rinv_packed = 0.d0
    do gemm_ia=1,gemm_natom
    do gemm_p=1,gemm_nproj_atom(gemm_ia)
      gemm_ilma = gemm_l2g(gemm_p, gemm_ia)
      gemm_rinv_packed(gemm_p, gemm_ia) = ppg%rinv_uvu(gemm_ilma)
    end do
    end do
!$acc enter data copyin(gemm_zekr_packed, gemm_rinv_packed, gemm_l2g, gemm_nproj_atom)

    allocate(gemm_wf_packed(ppg%nps, gemm_io_block, gemm_natom))
    allocate(gemm_out_packed(gemm_max_nproj, gemm_io_block, gemm_natom))
    gemm_wf_packed = (0.d0, 0.d0)
!$acc enter data copyin(gemm_wf_packed) create(gemm_out_packed)

    gemm_stat = cublasCreate(gemm_handle)
  end if

  ! cublas must share OpenACC's stream, or the scatter kernel may read the GEMM output early.
  gemm_stat = cublasSetStream(gemm_handle, acc_get_cuda_stream(acc_async_sync))

  ! Phase 1: batched GEMM projection, blocked over bands to bound gemm_wf_packed's memory.
  do im=im_s,im_e
  do ik=ik_s,ik_e
    ! Refresh per (call, k-point): zekr_uV carries the time-dependent A(t), so packing it
    ! only at setup would freeze the t=0 phases and silently drift the physics.
!$acc parallel loop collapse(3) present(ppg, gemm_zekr_packed, gemm_l2g, gemm_nproj_atom)
    do gemm_ia = 1, gemm_natom
    do gemm_p = 1, gemm_max_nproj
    do gemm_j = 1, ppg%nps
      if (gemm_p <= gemm_nproj_atom(gemm_ia) .and. gemm_j <= ppg%mps(gemm_ia)) then
        gemm_zekr_packed(gemm_j, gemm_p, gemm_ia) = &
            ppg%zekr_uV(gemm_j, gemm_l2g(gemm_p, gemm_ia), ik)
      end if
    end do
    end do
    end do

  do ispin=1,Nspin
    gemm_io_blk_s = io_s
    do while (gemm_io_blk_s <= io_e)
      gemm_this_block = min(gemm_io_block, io_e - gemm_io_blk_s + 1)

!$acc parallel loop collapse(3) present(tpsi,ppg,gemm_wf_packed)
      do gemm_ia = 1, gemm_natom
      do gemm_io_local = 1, gemm_this_block
      do gemm_j = 1, ppg%nps
        if (gemm_j <= ppg%mps(gemm_ia)) then
          gemm_wf_packed(gemm_j, gemm_io_local, gemm_ia) = tpsi%zwf( &
              ppg%jxyz(1,gemm_j,gemm_ia), ppg%jxyz(2,gemm_j,gemm_ia), ppg%jxyz(3,gemm_j,gemm_ia), &
              ispin, gemm_io_blk_s+gemm_io_local-1, ik, im)
        end if
      end do
      end do
      end do

!$acc host_data use_device(gemm_zekr_packed, gemm_wf_packed, gemm_out_packed)
      gemm_stat = cublasZgemmStridedBatched(gemm_handle, CUBLAS_OP_C, CUBLAS_OP_N, &
          gemm_max_nproj, gemm_this_block, ppg%nps, &
          (1.d0,0.d0), gemm_zekr_packed, ppg%nps, int(ppg%nps,8)*int(gemm_max_nproj,8), &
          gemm_wf_packed, ppg%nps, int(ppg%nps,8)*int(gemm_io_block,8), &
          (0.d0,0.d0), gemm_out_packed, gemm_max_nproj, int(gemm_max_nproj,8)*int(gemm_io_block,8), &
          gemm_natom)
!$acc end host_data

!$acc parallel loop collapse(3) present(ppg,gemm_out_packed,gemm_rinv_packed,gemm_l2g,gemm_nproj_atom)
      do gemm_ia = 1, gemm_natom
      do gemm_io_local = 1, gemm_this_block
      do gemm_p = 1, gemm_max_nproj
        if (gemm_p <= gemm_nproj_atom(gemm_ia)) then
          ppg%uVpsibox(gemm_l2g(gemm_p,gemm_ia), ispin, gemm_io_blk_s+gemm_io_local-1, ik, im) = &
              gemm_out_packed(gemm_p, gemm_io_local, gemm_ia) * gemm_rinv_packed(gemm_p, gemm_ia)
        end if
      end do
      end do
      end do

      gemm_io_blk_s = gemm_io_blk_s + gemm_this_block
    end do
  end do
  end do
  end do

  ! Phase 1.5: reduce each rank's partial uVpsi across rgrid ranks. Staged through a CUDA
  ! Fortran device buffer so CUDA-aware MPI sees a real device pointer, not managed memory.
  if(info%if_divide_rspace) then
    call timer_end(LOG_UHPSI_PSEUDO)
    call timer_begin(LOG_UHPSI_PSEUDO_COMM)

    norb = Nspin * info%numo * info%numk * info%numm

    if (d_reduce_allocated /= Nlma*norb) then
      if (d_reduce_allocated > 0) deallocate(d_reduce)
      allocate(d_reduce(Nlma, Nspin, io_e-io_s+1, ik_e-ik_s+1, im_e-im_s+1))
      d_reduce_allocated = Nlma*norb
    end if

    !$acc parallel loop collapse(5) deviceptr(d_reduce) present(ppg)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
    do ilma=1,Nlma
      d_reduce(ilma,ispin,io-io_s+1,ik-ik_s+1,im-im_s+1) = ppg%uVpsibox(ilma,ispin,io,ik,im)
    end do
    end do
    end do
    end do
    end do
    !$acc wait

#ifdef USE_NCCL
    if (.not. nccl_ready) then
      call MPI_Comm_rank(info%icomm_r, nccl_rank, mpi_ierr)
      call MPI_Comm_size(info%icomm_r, nccl_size, mpi_ierr)
      if (nccl_rank == 0) nccl_stat = ncclGetUniqueId(nccl_uid)
      call MPI_Bcast(nccl_uid%internal, 128, MPI_CHARACTER, 0, info%icomm_r, mpi_ierr)
      nccl_stat = ncclCommInitRank(nccl_comm, nccl_size, nccl_uid, nccl_rank)
      if (nccl_stat /= ncclSuccess) stop 'zpseudo: ncclCommInitRank failed'
      cuda_stat = cudaStreamCreate(nccl_stream)
      nccl_ready = .true.
    end if
    ! complex sum == double sum over twice as many elements
    call c_f_pointer(c_devloc(d_reduce), d_reduce_re, [2*Nlma*norb])
    nccl_stat = ncclAllReduce(d_reduce_re, d_reduce_re, 2*Nlma*norb, ncclDouble, ncclSum, nccl_comm, nccl_stream)
    cuda_stat = cudaStreamSynchronize(nccl_stream)
#else
    call MPI_Allreduce(MPI_IN_PLACE, d_reduce, int(Nlma*norb), MPI_DOUBLE_COMPLEX, MPI_SUM, info%icomm_r, mpi_ierr)
#endif

    !$acc parallel loop collapse(5) deviceptr(d_reduce) present(ppg)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
    do ilma=1,Nlma
      ppg%uVpsibox(ilma,ispin,io,ik,im) = d_reduce(ilma,ispin,io-io_s+1,ik-ik_s+1,im-im_s+1)
    end do
    end do
    end do
    end do
    end do
    !$acc wait

    call timer_end(LOG_UHPSI_PSEUDO_COMM)
    call timer_begin(LOG_UHPSI_PSEUDO)
  end if

  ! Phase 2 (back-projection): v2j-based scatter, no atomics.
!$acc kernels present(ppg,tpsi,htpsi)
!$acc loop collapse(5) independent private(ilocal,ilma,ia,uVpsi,vi,my_nlma,k,j,ix,iy,iz,wrk)
  do im=im_s,im_e
  do ik=ik_s,ik_e
  do io=io_s,io_e
  do ispin=1,Nspin
    do vi=0,ppg%max_vi-1
      my_nlma = ppg%v2nlma(vi)
      if (my_nlma < 1) cycle

      wrk = 0d0
!$acc loop seq
      do k=1,my_nlma
        ilma = ppg%k2ilma(vi,k)
        j    = ppg%k2j(vi,k)
        wrk  = wrk + ppg%uVpsibox(ilma,ispin,io,ik,im) * ppg%zekr_uV(j,ilma,ik)
      end do

      ix = ppg%v2j(1,vi)
      iy = ppg%v2j(2,vi)
      iz = ppg%v2j(3,vi)
      htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
    end do
  end do
  end do
  end do
  end do
!$acc end kernels

#else
  if(info%if_divide_rspace) then

    call timer_end(LOG_UHPSI_PSEUDO)

    call calc_uVpsi_rdivided(nspin,info,ppg,tpsi,uVpsibox,uVpsibox2)

    call timer_begin(LOG_UHPSI_PSEUDO)

#ifdef USE_OPENACC
! `acc atomic` does not support complex(8). Split it into the real and imag part.
    if (.not. allocated(htpsi_zwf_r)) then
      allocate( &
        htpsi_zwf_r( &
        lbound(htpsi%zwf,1):ubound(htpsi%zwf,1), &
        lbound(htpsi%zwf,2):ubound(htpsi%zwf,2), &
        lbound(htpsi%zwf,3):ubound(htpsi%zwf,3), &
        lbound(htpsi%zwf,4):ubound(htpsi%zwf,4), &
        lbound(htpsi%zwf,5):ubound(htpsi%zwf,5), &
        lbound(htpsi%zwf,6):ubound(htpsi%zwf,6), &
        lbound(htpsi%zwf,7):ubound(htpsi%zwf,7)  &
        ) &
        )
    end if
    if (.not. allocated(htpsi_zwf_i)) then
      allocate( &
        htpsi_zwf_i( &
        lbound(htpsi%zwf,1):ubound(htpsi%zwf,1), &
        lbound(htpsi%zwf,2):ubound(htpsi%zwf,2), &
        lbound(htpsi%zwf,3):ubound(htpsi%zwf,3), &
        lbound(htpsi%zwf,4):ubound(htpsi%zwf,4), &
        lbound(htpsi%zwf,5):ubound(htpsi%zwf,5), &
        lbound(htpsi%zwf,6):ubound(htpsi%zwf,6), &
        lbound(htpsi%zwf,7):ubound(htpsi%zwf,7)  &
        ) &
        )
    end if
    mps_max = size(ppg%jxyz, 2)
!$acc kernels
    htpsi_zwf_i = 0d0
    htpsi_zwf_r = 0d0

!$acc loop collapse(6) independent gang vector private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz,wrk,uVpsi)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
      do ilocal=1,ppg%ilocal_nlma
!OCL norecurrence
        do j=1,mps_max
          ilma=ppg%ilocal_nlma2ilma(ilocal)
          ia  =ppg%ilocal_nlma2ia  (ilocal)
          uVpsi = uVpsibox2(ispin,io,ik,im,ilma)
          if (j <= ppg%mps(ia)) then
            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)
            wrk = uVpsi * ppg%zekr_uV(j,ilma,ik)

            wrk_r = real(wrk)
            !$acc atomic update
            htpsi_zwf_r(ix,iy,iz,ispin,io,ik,im) = htpsi_zwf_r(ix,iy,iz,ispin,io,ik,im) + wrk_r
            !$acc end atomic

            wrk_r = aimag(wrk)
            !$acc atomic update
            htpsi_zwf_i(ix,iy,iz,ispin,io,ik,im) = htpsi_zwf_i(ix,iy,iz,ispin,io,ik,im) + wrk_r
            !$acc end atomic
          end if
        end do
      end do

    end do
    end do
    end do
    end do
!$acc loop collapse(7) independent gang vector
    do i1 = lbound(htpsi%zwf,1), ubound(htpsi%zwf,1)
    do i2 = lbound(htpsi%zwf,2), ubound(htpsi%zwf,2)
    do i3 = lbound(htpsi%zwf,3), ubound(htpsi%zwf,3)
    do i4 = lbound(htpsi%zwf,4), ubound(htpsi%zwf,4)
    do i5 = lbound(htpsi%zwf,5), ubound(htpsi%zwf,5)
    do i6 = lbound(htpsi%zwf,6), ubound(htpsi%zwf,6)
    do i7 = lbound(htpsi%zwf,7), ubound(htpsi%zwf,7)

      htpsi%zwf(i1,i2,i3,i4,i5,i6,i7) = htpsi%zwf(i1,i2,i3,i4,i5,i6,i7) + cmplx( &
          htpsi_zwf_r(i1,i2,i3,i4,i5,i6,i7), &
          htpsi_zwf_i(i1,i2,i3,i4,i5,i6,i7), &
          kind=8 &
      )

    end do
    end do
    end do
    end do
    end do
    end do
    end do
!$acc end kernels
#else
!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz,wrk)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do ilocal=1,ppg%ilocal_nlma
        ilma=ppg%ilocal_nlma2ilma(ilocal)
        ia  =ppg%ilocal_nlma2ia  (ilocal)
        uVpsi = uVpsibox2(ispin,io,ik,im,ilma)
!OCL norecurrence
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          wrk = uVpsi * ppg%zekr_uV(j,ilma,ik)

          htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
        end do
      end do

    end do
    end do
    end do
    end do
!$omp end parallel do
#endif

    deallocate(uVpsibox,uVpsibox2)

  else
#ifdef USE_OPENACC
#ifdef USE_CUDA
    natom=size(ppg%mps)
    n=size(ppg%rinv_uvu)/natom
    call zpseudo_cuda(htpsi%zwf,&
        n,&
        im_s,im_e,&
        ik_s,ik_e,&
        io_s,io_e,&
        Nspin,&
        Nlma,&
        ppg%nps,&
        natom,&
        lbound(htpsi%zwf,1),ubound(htpsi%zwf,1),&
        lbound(htpsi%zwf,2),ubound(htpsi%zwf,2),&
        lbound(htpsi%zwf,3),ubound(htpsi%zwf,3),&
        ppg%ia_tbl,&
        ppg%mps,&
        ppg%jxyz,&
        ppg%zekr_uV,&
        ppg%rinv_uvu,&
        tpsi%zwf)
#else
!$acc kernels present(ppg,tpsi,htpsi)
!$acc loop collapse(5) independent gang private(ilocal,ilma,ia,uVpsi,vi,my_nlma,k,j,ix,iy,iz,wrk)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi = 0.d0
!$acc loop vector(32) reduction(+:uVpsi)
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          uVpsi = uVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
        end do
        ppg%uVpsibox(ilma,ispin,io,ik,im) = uVpsi * ppg%rinv_uvu(ilma)
      end do
#ifdef USE_OPENACC
    end do
    end do
    end do
    end do
!$acc loop collapse(5) independent private(ilocal,ilma,ia,uVpsi,vi,my_nlma,k,j,ix,iy,iz,wrk)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
#endif
      do vi=0,ppg%max_vi-1
        my_nlma = ppg%v2nlma(vi)
        if (my_nlma < 1) cycle

        wrk = 0d0
!$acc loop seq
        do k=1,my_nlma
          ilma = ppg%k2ilma(vi,k)
          j    = ppg%k2j(vi,k)
          wrk  = wrk + ppg%uVpsibox(ilma,ispin,io,ik,im) * ppg%zekr_uV(j,ilma,ik)
        end do

        ix = ppg%v2j(1,vi)
        iy = ppg%v2j(2,vi)
        iz = ppg%v2j(3,vi)
        htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
      end do
    end do
    end do
    end do
    end do
!$acc end kernels
#endif
#else
!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilma,ia,uVpsi,j,ix,iy,iz,wrk)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi = 0.d0
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          uVpsi = uVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(ilma)
!OCL norecurrence
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          wrk = uVpsi * ppg%zekr_uV(j,ilma,ik)
          htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
        end do
      end do

    end do
    end do
    end do
    end do
!$omp end parallel do
#endif

  end if
#endif

  call nvtxEndRange()
  call timer_end(LOG_UHPSI_PSEUDO)

  return
end subroutine zpseudo

subroutine calc_uVpsi(nspin,info,ppg,tpsi,uVpsibox)
  use structures
  use timer
  implicit none
  integer        ,intent(in) :: nspin
  type(s_parallel_info),intent(in) :: info
  type(s_pp_grid),intent(in) :: ppg
  type(s_orbital),intent(in) :: tpsi
  complex(8)    ,allocatable :: uVpsibox (:,:,:,:,:)
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb
  integer :: ilma,ia,j,ix,iy,iz,Nlma
  complex(8) :: uVpsi
  integer :: ilocal

  call timer_begin(LOG_UHPSI_PSEUDO)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  norb = Nspin* info%numo * info%numk * info%numm

  Nlma = ppg%Nlma

  allocate(uVpsibox(Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e,Nlma))

#ifdef USE_OPENACC
!$acc parallel loop collapse(4) private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz)
#else
!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz)
#endif
  do im=im_s,im_e
  do ik=ik_s,ik_e
  do io=io_s,io_e
  do ispin=1,Nspin

    do ilocal=1,ppg%ilocal_nlma
      ilma=ppg%ilocal_nlma2ilma(ilocal)
      ia  =ppg%ilocal_nlma2ia  (ilocal)
      uVpsi = 0.d0
      do j=1,ppg%mps(ia)
        ix = ppg%jxyz(1,j,ia)
        iy = ppg%jxyz(2,j,ia)
        iz = ppg%jxyz(3,j,ia)
        uVpsi = uVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
      end do
      uVpsi = uVpsi * ppg%rinv_uvu(ilma)
      uVpsibox(ispin,io,ik,im,ilma) = uVpsi
    end do

  end do
  end do
  end do
  end do
#ifdef USE_OPENACC
!$acc end parallel
#else
!$omp end parallel do
#endif

  call timer_end(LOG_UHPSI_PSEUDO)

  return
end subroutine calc_uVpsi

subroutine calc_uVpsi_rdivided(nspin,info,ppg,tpsi,uVpsibox,uVpsibox2)
  use structures
  use timer
#ifdef FORTRAN_COMPILER_HAS_MPI_VERSION3
  use salmon_global, only: natom
  use communication, only: comm_wait_all,comm_show_error
  use mpi, only: MPI_SUM,MPI_DOUBLE_COMPLEX
#else
  use communication, only: comm_summation
#endif
  implicit none
  integer        ,intent(in) :: nspin
  type(s_parallel_info),intent(in) :: info
  type(s_pp_grid),intent(in) :: ppg
  type(s_orbital),intent(in) :: tpsi
  complex(8)    ,allocatable :: uVpsibox (:,:,:,:,:)
  complex(8)    ,allocatable :: uVpsibox2(:,:,:,:,:)
#ifdef USE_OPENACC
  complex(8)    ,allocatable, save, device :: dev_uVpsibox (:,:,:,:,:)
  complex(8)    ,allocatable, save, device :: dev_uVpsibox2(:,:,:,:,:)
#endif
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb
  integer :: ilma,ia,j,ix,iy,iz,Nlma,ilocal
  complex(8) :: uVpsi
#ifdef FORTRAN_COMPILER_HAS_MPI_VERSION3
  integer :: ireqs(natom),nreq,ierr,is,ie,ns
#endif

  call nvtxStartRange('calc_uVpsi_rdivided', __LINE__)
  call timer_begin(LOG_UHPSI_PSEUDO)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  norb = Nspin* info%numo * info%numk * info%numm

  Nlma = ppg%Nlma

  allocate(uVpsibox (Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e,Nlma))
  allocate(uVpsibox2(Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e,Nlma))

#ifdef USE_OPENACC
  if (.not. allocated(dev_uVpsibox)) then
    allocate(dev_uVpsibox (Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e,Nlma))
  end if

  if (.not. allocated(dev_uVpsibox2)) then
    allocate(dev_uVpsibox2(Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e,Nlma))
  end if

!$acc kernels
  dev_uVpsibox = 0d0
  dev_uVpsibox2 = 0d0
!$acc end kernels
#else
#ifdef FORTRAN_COMPILER_HAS_MPI_VERSION3
  uVpsibox2 = 0d0
#else
  uVpsibox  = 0d0
#endif
#endif

#ifdef USE_OPENACC
!$acc kernels
#else
!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz)
#endif
  do im=im_s,im_e
  do ik=ik_s,ik_e
  do io=io_s,io_e
  do ispin=1,Nspin

    do ilocal=1,ppg%ilocal_nlma
      ilma=ppg%ilocal_nlma2ilma(ilocal)
      ia  =ppg%ilocal_nlma2ia  (ilocal)
      uVpsi = 0.d0
      do j=1,ppg%mps(ia)
        ix = ppg%jxyz(1,j,ia)
        iy = ppg%jxyz(2,j,ia)
        iz = ppg%jxyz(3,j,ia)
        uVpsi = uVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
      end do
      uVpsi = uVpsi * ppg%rinv_uvu(ilma)
#ifdef USE_OPENACC
      dev_uVpsibox(ispin,io,ik,im,ilma) = uVpsi
#else
      uVpsibox(ispin,io,ik,im,ilma) = uVpsi
#endif
    end do

  end do
  end do
  end do
  end do
#ifdef USE_OPENACC
!$acc end kernels
#else
!$omp end parallel do
#endif

  call timer_end(LOG_UHPSI_PSEUDO)

  call timer_begin(LOG_UHPSI_PSEUDO_COMM)
#ifdef FORTRAN_COMPILER_HAS_MPI_VERSION3
! FIXME: This subroutine uses MPI functions directly...
  nreq = 0
  do ia=1,natom
    is = ppg%irange_atom(1,ia)
    ie = ppg%irange_atom(2,ia)
    ns = ie - is + 1

    if (ppg%ireferred_atom(ia)) then
      nreq = nreq + 1
#ifdef USE_OPENACC
      call MPI_Iallreduce( dev_uVpsibox (1,io_s,ik_s,im_s,is) &
                         , dev_uVpsibox2(1,io_s,ik_s,im_s,is) &
                         , ns*norb, MPI_DOUBLE_COMPLEX, MPI_SUM, ppg%icomm_atom(ia) &
                         , ireqs(nreq), ierr )
#else
      call MPI_Iallreduce( uVpsibox (1,io_s,ik_s,im_s,is) &
                         , uVpsibox2(1,io_s,ik_s,im_s,is) &
                         , ns*norb, MPI_DOUBLE_COMPLEX, MPI_SUM, ppg%icomm_atom(ia) &
                         , ireqs(nreq), ierr )
#endif
      call comm_show_error(ierr)
    !else
      ! uvpsibox2(:,:,:,:,ppg%irange_ia(1:2,ia)) does not use in this process...
      ! We can skip self copy, but zero clear required
    end if
  end do
  call comm_wait_all(ireqs(1:nreq))
#ifdef USE_OPENACC
!$acc kernels
uVpsibox2 = dev_uVpsibox2
uVpsibox = dev_uVpsibox
!$acc end kernels
#endif
#else
#ifdef USE_OPENACC
!$acc kernels
uVpsibox2 = dev_uVpsibox2
uVpsibox = dev_uVpsibox
!$acc end kernels
#endif
  call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,info%icomm_r)
#endif
  call timer_end(LOG_UHPSI_PSEUDO_COMM)

  call nvtxEndRange()
  return
end subroutine calc_uVpsi_rdivided

end module nonlocal_potential
