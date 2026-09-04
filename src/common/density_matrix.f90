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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

#include "config.h"

module density_matrix
  implicit none
  integer,private,parameter :: Nd = 4

contains

  subroutine calc_density(system,rho,psi,info,mg)
    use structures
    use communication, only: comm_summation
    use parallelization, only: get_thread_id,get_nthreads
    use misc_routines, only: ceiling_pow2
    use timer
    use sym_rho_sub, only: sym_rho
    use salmon_global, only: yn_spinorbit
    use noncollinear_module, only: calc_dm_noncollinear, rot_dm_noncollinear
    use nvtx_wrapper
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_orbital),intent(in) :: psi
    type(s_scalar) :: rho(system%nspin,info%im_s:info%im_e)
    !
    integer :: im,ispin,ik,io,is(3),ie(3),nsize,nspin,tid,ix,iy,iz,nthreads
    real(8) :: wrk2
    real(8),allocatable :: wrk(:,:,:,:)
    call nvtxStartRange('calc_density', __LINE__)
    call timer_begin(LOG_DENSITY_CALC)
    
    if(yn_spinorbit=='y') then
      call calc_dm_noncollinear( psi, system, info, mg )
      call rot_dm_noncollinear( rho, system, mg )
      call timer_end(LOG_DENSITY_CALC)
      return
    end if
    
    nspin = system%nspin

#ifdef FORTRAN_COMPILER_HAS_2MB_ALIGNED_ALLOCATION
!dir$ attributes align : 2097152 :: wrk
#endif

    is = mg%is
    ie = mg%ie
    nsize = mg%num(1) * mg%num(2) * mg%num(3)
    nthreads = get_nthreads()

    allocate(wrk(is(1):ie(1),is(2):ie(2),is(3):ie(3),0:ceiling_pow2(nthreads)-1))
    call timer_end(LOG_DENSITY_CALC)

    if(allocated(psi%rwf)) then

      do im=info%im_s,info%im_e
      do ispin=1,nspin
        call timer_begin(LOG_DENSITY_CALC)
        tid = 0
#ifdef USE_OPENACC
!$acc kernels copyin(is,ie)
        wrk(:,:,:,tid) = 0.d0

!$acc loop collapse(3) gang vector private(wrk2,ik,io,iz,iy,ix)
        do iz=is(3),ie(3)
        do iy=is(2),ie(2)
        do ix=is(1),ie(1)
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          wrk2 = abs( psi%rwf(ix,iy,iz,ispin,io,ik,im) )**2
          wrk(ix,iy,iz,tid) = wrk(ix,iy,iz,tid) + wrk2 * system%rocc(io,ik,ispin)*system%wtk(ik)
        end do
        end do
        end do
        end do
        end do
!$acc end kernels

#else
!$omp parallel private(ik,io,iz,iy,ix,wrk2) firstprivate(tid)
!$      tid = get_thread_id()
        wrk(:,:,:,tid) = 0.d0

!$omp do collapse(4)
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
        do iz=is(3),ie(3)
        do iy=is(2),ie(2)
        do ix=is(1),ie(1)
          wrk2 = abs( psi%rwf(ix,iy,iz,ispin,io,ik,im) )**2
          wrk(ix,iy,iz,tid) = wrk(ix,iy,iz,tid) + wrk2 * system%rocc(io,ik,ispin)*system%wtk(ik)
        end do
        end do
        end do
        end do
        end do
!$omp end do

        ix = size(wrk,4)/2
        do while(ix > 0)
          if(tid < ix .and. tid + ix < nthreads) then
            wrk(:,:,:,tid) = wrk(:,:,:,tid) + wrk(:,:,:,tid + ix)
          end if
          ix = ix/2
!$omp barrier
        end do

!$omp end parallel
#endif
        call timer_end(LOG_DENSITY_CALC)

        call timer_begin(LOG_DENSITY_COMM_COLL)
        call comm_summation(wrk(:,:,:,0),rho(ispin,im)%f(:,:,:),nsize,info%icomm_ko)
        call timer_end(LOG_DENSITY_COMM_COLL)
      end do
      end do

    else

      do im=info%im_s,info%im_e
      do ispin=1,nspin
        call timer_begin(LOG_DENSITY_CALC)
        tid = 0
#ifdef USE_OPENACC
!$acc kernels copyin(is,ie)
        wrk(:,:,:,tid) = 0.d0

!$acc loop collapse(3) gang vector private(wrk2,ik,io,iz,iy,ix)
        do iz=is(3),ie(3)
        do iy=is(2),ie(2)
        do ix=is(1),ie(1)
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          wrk2 = abs( psi%zwf(ix,iy,iz,ispin,io,ik,im) )**2
          wrk(ix,iy,iz,tid) = wrk(ix,iy,iz,tid) + wrk2 * system%rocc(io,ik,ispin)*system%wtk(ik)
        end do
        end do
        end do
        end do
        end do
!$acc end kernels

#else
!$omp parallel private(ik,io,iz,iy,ix,wrk2) firstprivate(tid)
!$      tid = get_thread_id()
        wrk(:,:,:,tid) = 0.d0

!$omp do collapse(4)
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
        do iz=is(3),ie(3)
        do iy=is(2),ie(2)
        do ix=is(1),ie(1)
          wrk2 = abs( psi%zwf(ix,iy,iz,ispin,io,ik,im) )**2
          wrk(ix,iy,iz,tid) = wrk(ix,iy,iz,tid) + wrk2 * system%rocc(io,ik,ispin)*system%wtk(ik)
        end do
        end do
        end do
        end do
        end do
!$omp end do

        ix = size(wrk,4)/2
        do while(ix > 0)
          if(tid < ix .and. tid + ix < nthreads) then
            wrk(:,:,:,tid) = wrk(:,:,:,tid) + wrk(:,:,:,tid + ix)
          end if
          ix = ix/2
!$omp barrier
        end do

!$omp end parallel
#endif
        call timer_end(LOG_DENSITY_CALC)

        call timer_begin(LOG_DENSITY_COMM_COLL)
        call comm_summation(wrk(:,:,:,0),rho(ispin,im)%f(:,:,:),nsize,info%icomm_ko)
        call timer_end(LOG_DENSITY_COMM_COLL)

        call timer_begin(LOG_DENSITY_CALC)
        call sym_rho( rho(ispin,im)%f(:,:,:) )
        call timer_end(LOG_DENSITY_CALC)
      end do
      end do

    end if

    deallocate(wrk)
    call nvtxEndRange
    return
  end subroutine calc_density

!===================================================================================================================================

  subroutine calc_current(system,mg,stencil,info,srg,psi,ppg,curr)
    use structures
    use salmon_global, only: yn_jm,yn_spinorbit
    use sendrecv_grid, only: update_overlap_complex8
    use communication, only: comm_summation
    use nonlocal_potential, only: calc_uVpsi_rdivided
    use pseudo_pt_current_so, only: calc_current_nonlocal_so &
                                  , calc_current_nonlocal_rdivided_so
    use sym_vector_sub, only: sym_vector_xyz
    use code_optimization, only: current_omp_mode
    use timer
    use iso_c_binding
    use nvtx_wrapper
#if defined(USE_OPENACC) && defined(USE_GEMM)
    use cublas
    use openacc
#endif
    implicit none
#if defined(USE_OPENACC)
    interface
      subroutine stencil_current_core_gpu(ik_s,ik_e,io_s,io_e,vec_k,vec_Ac,is_array,ie_array &
                                         ,is,ie,idx,idy,idz,nabt,ispin,im,spin_len,psi,BT,rocc,wtk,jx,jy,jz) bind(c)
        import
        ! Input
        integer(c_int), value :: ik_s
        integer(c_int), value :: ik_e
        integer(c_int), value :: io_s
        integer(c_int), value :: io_e
        integer(c_int), value :: ispin
        integer(c_int), value :: im
        integer(c_int), value :: spin_len
        ! Output
        real(c_double), intent(inout) :: jx
        real(c_double), intent(inout) :: jy
        real(c_double), intent(inout) :: jz
        ! Input (ptr)
        real(c_double), intent(in) :: vec_k(3:ik_e-ik_e+1)
        real(c_double), intent(in) :: vec_Ac(3)
        integer(c_int), intent(in) :: is_array(3)
        integer(c_int), intent(in) :: ie_array(3)
        integer(c_int), intent(in) :: is(3)
        integer(c_int), intent(in) :: ie(3)
        integer(c_int), intent(in) :: idx(is(1)-Nd:ie(1)+Nd)
        integer(c_int), intent(in) :: idy(is(2)-Nd:ie(2)+Nd)
        integer(c_int), intent(in) :: idz(is(3)-Nd:ie(3)+Nd)
        real(c_double), intent(in) :: nabt(Nd,3)
        complex(c_double_complex), intent(in) :: psi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_arary(3))
        real(c_double), intent(in) :: BT(3,3)
        real(c_double), intent(in) :: rocc(io_e-io_s+1,ik_e-ik_s+1)
        real(c_double), intent(in) :: wtk(ik_e-ik_s+1)
      end  subroutine stencil_current_core_gpu
    end interface
#endif
    type(s_dft_system),intent(in) :: system
    type(s_rgrid)  ,intent(in) :: mg
    type(s_stencil),intent(in) :: stencil
    type(s_parallel_info),intent(in) :: info
    type(s_sendrecv_grid)      :: srg
    type(s_orbital)            :: psi
    type(s_pp_grid),intent(in) :: ppg
    real(8) :: curr(3,system%nspin,info%im_s:info%im_e)
    !
    integer :: ispin,im,ik,io,nspin,ngrid
    real(8),dimension(3) :: wrk1,wrk2,wrk3,wrk4
    real(8) :: BT(3,3),kAc(3)
    complex(8),allocatable :: uVpsibox (:,:,:,:,:)
    complex(8),allocatable :: uVpsibox2(:,:,:,:,:)
    complex(8),allocatable :: uVpsi(:)
    real(8) :: jx,jy,jz
    ! Locals for stencil_current inlined below (OpenACC, non-CUDA branch).
    integer    :: ix,iy,iz
    real(8)    :: rtmp
    complex(8) :: cpsi
#if defined(USE_OPENACC) && defined(USE_GEMM)
    ! Batched-GEMM nonlocal current: zpseudo's phase-1 contraction, but four
    ! projections per (ilma,io) -- zekr_uV plus x/y/z-weighted copies, stacked
    ! into one packed matrix (4x zpseudo's) so a single GEMM returns all four.
    integer, parameter :: cur_io_block = 64
    integer,                     save :: cur_max_nproj = -1
    integer,    allocatable,     save :: cur_nproj_atom(:), cur_l2g(:,:)
    complex(8), allocatable,     save :: cur_zekr4(:,:,:), cur_wf(:,:,:), cur_out(:,:,:)
    real(8),    allocatable,     save :: cur_rinv(:,:)
    type(cublasHandle),          save :: cur_handle
    integer    :: cur_ia,cur_p,cur_j,cur_b,cur_natom,cur_ilma,cur_stat
    integer    :: cur_blk_s,cur_nb,cur_np4
    complex(8) :: cur_uv,cur_z
    real(8)    :: cur_w
#endif
    call nvtxStartRange('calc_current', __LINE__)
    call timer_begin(LOG_CURRENT_CALC)
#ifdef FORTRAN_COMPILER_HAS_2MB_ALIGNED_ALLOCATION
!dir$ attributes align : 2097152 :: uVpsibox, uVpsibox2
#endif

    nspin = system%nspin
    ngrid = system%ngrid

    BT = transpose(system%rmatrix_B)
    call timer_end(LOG_CURRENT_CALC)

    call timer_begin(LOG_CURRENT_CALC_UVPSI_RDIVIDED)
    if (info%if_divide_rspace .and. yn_jm=='n' .and. .not. yn_spinorbit=='y') then
#if !defined(USE_OPENACC)
      call calc_uVpsi_rdivided(nspin,info,ppg,psi,uVpsibox,uVpsibox2)
      allocate(uVpsi(ppg%Nlma))
#endif
    end if
    call timer_end(LOG_CURRENT_CALC_UVPSI_RDIVIDED)

  ! overlap region communication
    call timer_begin(LOG_CURRENT_COMM_HALO)
    if(info%if_divide_rspace) then
      call update_overlap_complex8(srg, mg, psi%zwf)
    end if
    call timer_end(LOG_CURRENT_COMM_HALO)

    do im=info%im_s,info%im_e
    do ispin=1,nspin
      call timer_begin(LOG_CURRENT_CALC)
      wrk4 = 0d0
#if defined(USE_OPENACC)
      jx = 0d0
      jy = 0d0
      jz = 0d0
      call timer_begin(LOG_CALC_STENCIL_CURRENT)
#if defined(USE_OPENACC) && defined(USE_CUDA)
      call stencil_current_core_gpu(info%ik_s,info%ik_e,info%io_s,info%io_e,system%vec_k,system%vec_Ac &
                                   ,mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx(1:),mg%idy(1:),mg%idz(1:) &
                                   ,stencil%coef_nab,ispin,im,nspin,psi%zwf,BT,system%rocc,system%wtk,jx,jy,jz)
!$acc enter data copyin(jx,jy,jz)
#else
!$acc enter data copyin(jx,jy,jz)
!$acc update device(system%vec_Ac)

! Grid-parallel stencil: collapse (ik,io,grid) so gang count is
! nk*numo*ngrid, not nk*numo; weighting is linear so it moves inside.
!$acc kernels copyin(BT,ispin,im)
!$acc loop gang private(ik,io,ix,iy,iz,kAc,wrk1,wrk2,wrk3,wrk4,rtmp,cpsi) &
!$acc&         reduction(+:jx,jy,jz) collapse(5) independent
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        kAc(1:3) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
        cpsi = conjg(psi%zwf(ix,iy,iz,ispin,io,ik,im))
        rtmp = abs(psi%zwf(ix,iy,iz,ispin,io,ik,im))**2
        wrk1(1) = kAc(1)*rtmp
        wrk1(2) = kAc(2)*rtmp
        wrk1(3) = kAc(3)*rtmp
        wrk2(1) = aimag(2d0 * ( stencil%coef_nab(1,1) * cpsi * psi%zwf(mg%idx(ix+1),iy,iz,ispin,io,ik,im) &
                              + stencil%coef_nab(2,1) * cpsi * psi%zwf(mg%idx(ix+2),iy,iz,ispin,io,ik,im) &
                              + stencil%coef_nab(3,1) * cpsi * psi%zwf(mg%idx(ix+3),iy,iz,ispin,io,ik,im) &
                              + stencil%coef_nab(4,1) * cpsi * psi%zwf(mg%idx(ix+4),iy,iz,ispin,io,ik,im) ))
        wrk2(2) = aimag(2d0 * ( stencil%coef_nab(1,2) * cpsi * psi%zwf(ix,mg%idy(iy+1),iz,ispin,io,ik,im) &
                              + stencil%coef_nab(2,2) * cpsi * psi%zwf(ix,mg%idy(iy+2),iz,ispin,io,ik,im) &
                              + stencil%coef_nab(3,2) * cpsi * psi%zwf(ix,mg%idy(iy+3),iz,ispin,io,ik,im) &
                              + stencil%coef_nab(4,2) * cpsi * psi%zwf(ix,mg%idy(iy+4),iz,ispin,io,ik,im) ))
        wrk2(3) = aimag(2d0 * ( stencil%coef_nab(1,3) * cpsi * psi%zwf(ix,iy,mg%idz(iz+1),ispin,io,ik,im) &
                              + stencil%coef_nab(2,3) * cpsi * psi%zwf(ix,iy,mg%idz(iz+2),ispin,io,ik,im) &
                              + stencil%coef_nab(3,3) * cpsi * psi%zwf(ix,iy,mg%idz(iz+3),ispin,io,ik,im) &
                              + stencil%coef_nab(4,3) * cpsi * psi%zwf(ix,iy,mg%idz(iz+4),ispin,io,ik,im) ))
        wrk3(1) = BT(1,1)*wrk2(1) + BT(1,2)*wrk2(2) + BT(1,3)*wrk2(3)
        wrk3(2) = BT(2,1)*wrk2(1) + BT(2,2)*wrk2(2) + BT(2,3)*wrk2(3)
        wrk3(3) = BT(3,1)*wrk2(1) + BT(3,2)*wrk2(2) + BT(3,3)*wrk2(3)
        wrk4 = (wrk1 + wrk3) * system%rocc(io,ik,ispin) * system%wtk(ik)
        jx = jx + wrk4(1)
        jy = jy + wrk4(2)
        jz = jz + wrk4(3)
      end do
      end do
      end do
      end do
      end do
!$acc end kernels
#endif
      call timer_end(LOG_CALC_STENCIL_CURRENT)
      
      if ( yn_jm == 'n' ) then
        if ( yn_spinorbit=='y' ) then
          call timer_begin(LOG_CURRENT_SO_NONLOCAL)
!$acc kernels copyin(ispin,im) copy(jx,jy,jz)
!$acc loop gang private(ik,io,wrk3,wrk4) reduction(+:jx,jy,jz) collapse(2) independent
          do ik=info%ik_s,info%ik_e
          do io=info%io_s,info%io_e
            call calc_current_nonlocal_so &
                 ( wrk3,psi%zwf(:,:,:,:,io,ik,im),ppg,mg%is_array,mg%ie_array,ik )
            wrk4 = wrk3 * system%rocc(io,ik,ispin) * system%wtk(ik)
            jx = jx + wrk4(1)
            jy = jy + wrk4(2)
            jz = jz + wrk4(3)
          end do
          end do
!$acc end kernels
!$acc exit data copyout(jx,jy,jz)
          call timer_end(LOG_CURRENT_SO_NONLOCAL)
        else ! yn_spinorbit=='y'
#if defined(USE_GEMM)
          cur_natom = size(ppg%mps)
          if (cur_max_nproj < 0) then
            allocate(cur_nproj_atom(cur_natom))
            cur_nproj_atom = 0
            do cur_ilma=1,ppg%Nlma
              cur_ia = ppg%ia_tbl(cur_ilma)
              cur_nproj_atom(cur_ia) = cur_nproj_atom(cur_ia) + 1
            end do
            cur_max_nproj = maxval(cur_nproj_atom)
            allocate(cur_l2g(cur_max_nproj, cur_natom))
            allocate(cur_rinv(cur_max_nproj, cur_natom))
            cur_l2g  = 0
            cur_rinv = 0d0
            cur_nproj_atom = 0   ! reused as a per-atom fill cursor
            do cur_ilma=1,ppg%Nlma
              cur_ia = ppg%ia_tbl(cur_ilma)
              cur_nproj_atom(cur_ia) = cur_nproj_atom(cur_ia) + 1
              cur_l2g (cur_nproj_atom(cur_ia), cur_ia) = cur_ilma
              cur_rinv(cur_nproj_atom(cur_ia), cur_ia) = ppg%rinv_uvu(cur_ilma)
            end do
            allocate(cur_zekr4(ppg%nps, 4*cur_max_nproj, cur_natom))
            allocate(cur_wf   (ppg%nps, cur_io_block,    cur_natom))
            allocate(cur_out  (4*cur_max_nproj, cur_io_block, cur_natom))
            cur_zekr4 = (0d0,0d0)
            cur_wf    = (0d0,0d0)
!$acc enter data copyin(cur_zekr4,cur_wf,cur_l2g,cur_nproj_atom,cur_rinv) create(cur_out)
            cur_stat = cublasCreate(cur_handle)
          end if
          cur_np4 = 4*cur_max_nproj
          ! cublas must share OpenACC's stream, or the reduction reads the GEMM early.
          cur_stat = cublasSetStream(cur_handle, acc_get_cuda_stream(acc_async_sync))

          do ik=info%ik_s,info%ik_e
            ! Refresh per (call, k-point): zekr_uV carries the time-dependent A(t).
!$acc parallel loop collapse(3) present(ppg,cur_zekr4,cur_l2g,cur_nproj_atom) private(cur_z)
            do cur_ia=1,cur_natom
            do cur_p=1,cur_max_nproj
            do cur_j=1,ppg%nps
              if (cur_p <= cur_nproj_atom(cur_ia) .and. cur_j <= ppg%mps(cur_ia)) then
                cur_z = ppg%zekr_uV(cur_j, cur_l2g(cur_p,cur_ia), ik)
                cur_zekr4(cur_j, cur_p,                   cur_ia) = cur_z
                cur_zekr4(cur_j, cur_p +   cur_max_nproj, cur_ia) = cur_z * ppg%Rxyz(1,cur_j,cur_ia)
                cur_zekr4(cur_j, cur_p + 2*cur_max_nproj, cur_ia) = cur_z * ppg%Rxyz(2,cur_j,cur_ia)
                cur_zekr4(cur_j, cur_p + 3*cur_max_nproj, cur_ia) = cur_z * ppg%Rxyz(3,cur_j,cur_ia)
              end if
            end do
            end do
            end do

            cur_blk_s = info%io_s
            do while (cur_blk_s <= info%io_e)
              cur_nb = min(cur_io_block, info%io_e - cur_blk_s + 1)
!$acc parallel loop collapse(3) present(psi,ppg,cur_wf)
              do cur_ia=1,cur_natom
              do cur_b=1,cur_nb
              do cur_j=1,ppg%nps
                if (cur_j <= ppg%mps(cur_ia)) then
                  cur_wf(cur_j,cur_b,cur_ia) = psi%zwf( &
                      ppg%jxyz(1,cur_j,cur_ia), ppg%jxyz(2,cur_j,cur_ia), ppg%jxyz(3,cur_j,cur_ia), &
                      ispin, cur_blk_s+cur_b-1, ik, im)
                end if
              end do
              end do
              end do

!$acc host_data use_device(cur_zekr4,cur_wf,cur_out)
              cur_stat = cublasZgemmStridedBatched(cur_handle, CUBLAS_OP_C, CUBLAS_OP_N, &
                  cur_np4, cur_nb, ppg%nps, &
                  (1d0,0d0), cur_zekr4, ppg%nps, int(ppg%nps,8)*int(cur_np4,8), &
                  cur_wf, ppg%nps, int(ppg%nps,8)*int(cur_io_block,8), &
                  (0d0,0d0), cur_out, cur_np4, int(cur_np4,8)*int(cur_io_block,8), &
                  cur_natom)
!$acc end host_data
              if (cur_stat /= CUBLAS_STATUS_SUCCESS) stop 'calc_current: cublasZgemmStridedBatched failed'

!$acc parallel loop collapse(2) present(cur_out,cur_rinv,cur_nproj_atom,system) &
!$acc&              private(cur_uv,cur_w,cur_p) reduction(+:jx,jy,jz)
              do cur_ia=1,cur_natom
              do cur_b=1,cur_nb
                cur_w = system%rocc(cur_blk_s+cur_b-1,ik,ispin) * system%wtk(ik)
!$acc loop seq
                do cur_p=1,cur_nproj_atom(cur_ia)
                  cur_uv = cur_out(cur_p,cur_b,cur_ia) * cur_rinv(cur_p,cur_ia)
                  jx = jx + 2d0*aimag(conjg(cur_out(cur_p +   cur_max_nproj,cur_b,cur_ia))*cur_uv)*cur_w
                  jy = jy + 2d0*aimag(conjg(cur_out(cur_p + 2*cur_max_nproj,cur_b,cur_ia))*cur_uv)*cur_w
                  jz = jz + 2d0*aimag(conjg(cur_out(cur_p + 3*cur_max_nproj,cur_b,cur_ia))*cur_uv)*cur_w
                end do
              end do
              end do

              cur_blk_s = cur_blk_s + cur_nb
            end do
          end do
#else
!$acc kernels copyin(ispin,im)
!$acc loop gang private(ik,io,wrk3,wrk4) reduction(+:jx,jy,jz) collapse(2) independent
          do ik=info%ik_s,info%ik_e
          do io=info%io_s,info%io_e
            call calc_current_nonlocal(wrk3,psi%zwf(:,:,:,ispin,io,ik,im),ppg,mg%is_array,mg%ie_array,ik)
            wrk4 = wrk3 * system%rocc(io,ik,ispin) * system%wtk(ik)
            jx = jx + wrk4(1)
            jy = jy + wrk4(2)
            jz = jz + wrk4(3)
          end do
          end do
!$acc end kernels
#endif
!$acc exit data copyout(jx,jy,jz)
        end if ! yn_spinorbit=='y'
      else ! yn_jm == 'n'
        wrk3=0.d0
      end if ! yn_jm == 'n'

      wrk4(1) = jx
      wrk4(2) = jy
      wrk4(3) = jz
#else
!$omp parallel do collapse(2) default(none) &
!$omp             private(ik,io,kAc,wrk1,wrk2,wrk3,uVpsi) &
!$omp             shared(info,system,mg,stencil,ppg,psi,uVpsibox2,BT,im,ispin,yn_jm) &
!$omp             shared(yn_spinorbit) &
!$omp             reduction(+:wrk4) if(current_omp_mode)
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e

        kAc(1:3) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
        call stencil_current(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz,stencil%coef_nab &
                            ,kAc,psi%zwf(:,:,:,ispin,io,ik,im),wrk1,wrk2)
        wrk2 = matmul(BT,wrk2)

        if ( yn_jm == 'n' ) then
          if ( yn_spinorbit=='y' ) then
            if ( info%if_divide_rspace ) then
              call calc_current_nonlocal_rdivided_so &
                   ( wrk3,psi%zwf(:,:,:,:,io,ik,im),ppg,mg%is_array,mg%ie_array,ik,info%icomm_r )
            else
              call calc_current_nonlocal_so &
                   ( wrk3,psi%zwf(:,:,:,:,io,ik,im),ppg,mg%is_array,mg%ie_array,ik )
            end if
          else
            if ( info%if_divide_rspace)then
              uVpsi(:) = uVpsibox2(ispin,io,ik,im,:)
              call calc_current_nonlocal_rdivided(wrk3,psi%zwf(:,:,:,ispin,io,ik,im),ppg,mg%is_array,mg%ie_array,ik,uVpsi)
            else
              call calc_current_nonlocal         (wrk3,psi%zwf(:,:,:,ispin,io,ik,im),ppg,mg%is_array,mg%ie_array,ik)
            end if
          end if
        else
          wrk3=0d0
        end if

        wrk4 = wrk4 + (wrk1 + wrk2 + wrk3) * system%rocc(io,ik,ispin)*system%wtk(ik)

      end do
      end do
!$omp end parallel do
#endif
      call timer_end(LOG_CURRENT_CALC)

      call timer_begin(LOG_CURRENT_COMM_COLL)
      call comm_summation(wrk4,wrk1,3,info%icomm_rko)

      curr(:,ispin,im) = wrk1 / dble(ngrid) ! ngrid = aLxyz/Hxyz
      call timer_end(LOG_CURRENT_COMM_COLL)

      call timer_begin(LOG_CURRENT_CALC)
      call sym_vector_xyz( curr(:,ispin,im) )
      call timer_end(LOG_CURRENT_CALC)
    end do
    end do

    if (info%if_divide_rspace .and. yn_jm=='n' .and. .not. yn_spinorbit=='y') then
#if !defined(USE_OPENACC)
      deallocate(uVpsibox,uVpsibox2,uVpsi)
#endif
    end if

    call nvtxEndRange
    return

  end subroutine calc_current

  subroutine stencil_current(is_array,ie_array,is,ie,idx,idy,idz,nabt,kAc,psi,j1,j2)
    !$acc routine vector
    integer   ,intent(in) :: is_array(3),ie_array(3),is(3),ie(3) &
                            ,idx(is(1)-Nd:ie(1)+Nd),idy(is(2)-Nd:ie(2)+Nd),idz(is(3)-Nd:ie(3)+Nd)
    real(8)   ,intent(in) :: nabt(Nd,3),kAc(3)
    complex(8),intent(in) :: psi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
    real(8)               :: j1(3),j2(3)
    !
    integer :: ix,iy,iz
    real(8) :: rtmp
    complex(8) :: cpsi,xtmp,ytmp,ztmp
    rtmp = 0d0
    xtmp = 0d0
    ytmp = 0d0
    ztmp = 0d0
#ifdef USE_OPENACC
!$acc loop vector collapse(3) private(iz,iy,ix,cpsi) reduction(+:rtmp,xtmp,ytmp,ztmp)
#else
!$omp parallel do collapse(2) private(iz,iy,ix,cpsi) reduction(+:rtmp,xtmp,ytmp,ztmp)
#endif
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)

!OCL swp
    do ix=is(1),ie(1)
      rtmp = rtmp + abs(psi(ix,iy,iz))**2
#ifndef USE_OPENACC
    end do

!OCL swp
    do ix=is(1),ie(1)
#endif
      cpsi = conjg(psi(ix,iy,iz))
      xtmp = xtmp + nabt(1,1) * cpsi * psi(idx(ix+1),iy,iz) &
                  + nabt(2,1) * cpsi * psi(idx(ix+2),iy,iz) &
                  + nabt(3,1) * cpsi * psi(idx(ix+3),iy,iz) &
                  + nabt(4,1) * cpsi * psi(idx(ix+4),iy,iz)
#ifndef USE_OPENACC
    end do

!OCL swp
    do ix=is(1),ie(1)
#endif
      cpsi = conjg(psi(ix,iy,iz))
      ytmp = ytmp + nabt(1,2) * cpsi * psi(ix,idy(iy+1),iz) &
                  + nabt(2,2) * cpsi * psi(ix,idy(iy+2),iz) &
                  + nabt(3,2) * cpsi * psi(ix,idy(iy+3),iz) &
                  + nabt(4,2) * cpsi * psi(ix,idy(iy+4),iz)
#ifndef USE_OPENACC
    end do

!OCL swp
    do ix=is(1),ie(1)
#endif
      cpsi = conjg(psi(ix,iy,iz))
      ztmp = ztmp + nabt(1,3) * cpsi * psi(ix,iy,idz(iz+1)) &
                  + nabt(2,3) * cpsi * psi(ix,iy,idz(iz+2)) &
                  + nabt(3,3) * cpsi * psi(ix,iy,idz(iz+3)) &
                  + nabt(4,3) * cpsi * psi(ix,iy,idz(iz+4))
    end do

    end do
    end do
#ifndef USE_OPENACC
!$omp end parallel do
#endif
    ! j1 = kAc(:) * rtmp
    j1(1) = kAc(1) * rtmp
    j1(2) = kAc(2) * rtmp
    j1(3) = kAc(3) * rtmp
    j2(1) = aimag(xtmp * 2d0)
    j2(2) = aimag(ytmp * 2d0)
    j2(3) = aimag(ztmp * 2d0)
    return
  end subroutine stencil_current

  subroutine calc_current_nonlocal(jw,psi,ppg,is_array,ie_array,ik)
    !$acc routine vector
    use structures
    implicit none
    integer   ,intent(in) :: is_array(3),ie_array(3),ik
    complex(8),intent(in) :: psi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
    type(s_pp_grid),intent(in) :: ppg
    real(8)               :: jw(3)
    real(8)               :: jw_1, jw_2, jw_3
    !
    integer    :: ilma,ia,j,ix,iy,iz
    real(8)    :: x,y,z
    complex(8) :: uVpsi,uVpsi_r(3)
#ifdef USE_OPENACC
    jw_1 = 0d0
    jw_2 = 0d0
    jw_3 = 0d0
!$acc loop vector private(ilma,ia,uVpsi,uVpsi_r,j,x,y,z,ix,iy,iz) reduction(+:jw_1, jw_2, jw_3)
#else
    jw = 0d0
!$omp parallel do private(ilma,ia,uVpsi,uVpsi_r,j,x,y,z,ix,iy,iz) reduction(+:jw)
#endif
    do ilma=1,ppg%Nlma
      ia=ppg%ia_tbl(ilma)
      uVpsi = 0d0
      uVpsi_r(1) = 0d0
      uVpsi_r(2) = 0d0
      uVpsi_r(3) = 0d0

      do j=1,ppg%Mps(ia)
        x = ppg%Rxyz(1,j,ia)
        y = ppg%Rxyz(2,j,ia)
        z = ppg%Rxyz(3,j,ia)
        ix = ppg%Jxyz(1,j,ia)
        iy = ppg%Jxyz(2,j,ia)
        iz = ppg%Jxyz(3,j,ia)
        uVpsi = uVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * psi(ix,iy,iz)
        uVpsi_r(1) = uVpsi_r(1) + conjg(ppg%zekr_uV(j,ilma,ik)) * x * psi(ix,iy,iz)
        uVpsi_r(2) = uVpsi_r(2) + conjg(ppg%zekr_uV(j,ilma,ik)) * y * psi(ix,iy,iz)
        uVpsi_r(3) = uVpsi_r(3) + conjg(ppg%zekr_uV(j,ilma,ik)) * z * psi(ix,iy,iz)
      end do
      uVpsi = uVpsi * ppg%rinv_uvu(ilma)
#ifdef USE_OPENACC
      jw_1 = jw_1 + aimag(conjg(uVpsi_r(1))*uVpsi)
      jw_2 = jw_2 + aimag(conjg(uVpsi_r(2))*uVpsi)
      jw_3 = jw_3 + aimag(conjg(uVpsi_r(3))*uVpsi)
#else
      jw = jw + aimag(conjg(uVpsi_r)*uVpsi)
#endif
    end do
#ifdef USE_OPENACC
    jw(1) = jw_1 * 2d0
    jw(2) = jw_2 * 2d0
    jw(3) = jw_3 * 2d0
#else
!$omp end parallel do
    jw = jw * 2d0
#endif
    return
  end subroutine calc_current_nonlocal

  subroutine calc_current_nonlocal_rdivided(jw,psi,ppg,is_array,ie_array,ik,uVpsibox)
    use structures
    implicit none
    integer   ,intent(in) :: is_array(3),ie_array(3),ik
    complex(8),intent(in) :: psi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
    type(s_pp_grid),intent(in) :: ppg
    complex(8),intent(in) :: uVpsibox(ppg%Nlma)
    real(8)               :: jw(3)
    !
    integer    :: ilocal,ilma,ia,j,ix,iy,iz
    real(8)    :: x,y,z
    complex(8) :: uVpsi,uVpsi_r(3)
    jw = 0d0
!$omp parallel do private(ilocal,ilma,ia,uVpsi,uVpsi_r,j,x,y,z,ix,iy,iz) reduction(+:jw)
    do ilocal=1,ppg%ilocal_nlma
      ilma=ppg%ilocal_nlma2ilma(ilocal)
      ia  =ppg%ilocal_nlma2ia  (ilocal)
      uVpsi_r = 0d0
      do j=1,ppg%Mps(ia)
        x = ppg%Rxyz(1,j,ia)
        y = ppg%Rxyz(2,j,ia)
        z = ppg%Rxyz(3,j,ia)
        ix = ppg%Jxyz(1,j,ia)
        iy = ppg%Jxyz(2,j,ia)
        iz = ppg%Jxyz(3,j,ia)
        uVpsi_r(1) = uVpsi_r(1) + conjg(ppg%zekr_uV(j,ilma,ik)) * x * psi(ix,iy,iz)
        uVpsi_r(2) = uVpsi_r(2) + conjg(ppg%zekr_uV(j,ilma,ik)) * y * psi(ix,iy,iz)
        uVpsi_r(3) = uVpsi_r(3) + conjg(ppg%zekr_uV(j,ilma,ik)) * z * psi(ix,iy,iz)
      end do
      uVpsi = uVpsibox(ilma)
      jw = jw + aimag(conjg(uVpsi_r)*uVpsi)
    end do
!$omp end parallel do
    jw = jw * 2d0
    return
  end subroutine calc_current_nonlocal_rdivided
  
!===================================================================================================================================
  
! curr(ispin) = \sum_{ik,io} [ system%rocc(io,ik,ispin)*system%wtk(ik)* curr_decomp(ispin,io,ik) ]
! curr: the current density
! curr_decomp: decomposition of the current density
  subroutine calc_current_decomposed(system,mg,stencil,info,srg,psi,ppg,curr_decomp)
    use structures
    use salmon_global, only: yn_jm,yn_spinorbit
    use sendrecv_grid, only: update_overlap_complex8
    use communication, only: comm_summation
    use nonlocal_potential, only: calc_uVpsi_rdivided
    use pseudo_pt_current_so, only: calc_current_nonlocal_so,calc_current_nonlocal_rdivided_so
    implicit none
    type(s_dft_system)   ,intent(in) :: system
    type(s_rgrid)        ,intent(in) :: mg
    type(s_stencil)      ,intent(in) :: stencil
    type(s_parallel_info),intent(in) :: info
    type(s_sendrecv_grid)            :: srg
    type(s_orbital)                  :: psi
    type(s_pp_grid)      ,intent(in) :: ppg
    real(8)                          :: curr_decomp(3,system%nspin,system%no,system%nk)
    !
    integer :: ispin,im,ik,io,nspin,ngrid
    real(8),dimension(3) :: wrk1,wrk2,wrk3
    real(8) :: BT(3,3),kAc(3)
    real(8) :: curr_wrk(3,system%nspin,system%no,system%nk)
    complex(8),allocatable :: uVpsibox (:,:,:,:,:)
    complex(8),allocatable :: uVpsibox2(:,:,:,:,:)
    complex(8),allocatable :: uVpsi(:)
    real(8) :: jx,jy,jz

    im = 1
    nspin = system%nspin
    ngrid = system%ngrid

    BT = transpose(system%rmatrix_B)

    if (info%if_divide_rspace .and. yn_jm=='n' .and. .not. yn_spinorbit=='y') then
      call calc_uVpsi_rdivided(nspin,info,ppg,psi,uVpsibox,uVpsibox2)
      allocate(uVpsi(ppg%Nlma))
    end if

  ! overlap region communication
    if(info%if_divide_rspace) then
      call update_overlap_complex8(srg, mg, psi%zwf)
    end if

    curr_wrk = 0d0
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    
      do ispin=1,nspin
        kAc(1:3) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
        call stencil_current(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz,stencil%coef_nab &
                            ,kAc,psi%zwf(:,:,:,ispin,io,ik,im),wrk1,wrk2)
        wrk2 = matmul(BT,wrk2)
        if ( yn_jm == 'n' ) then
          if ( yn_spinorbit=='y' ) then
            if ( info%if_divide_rspace ) then
              call calc_current_nonlocal_rdivided_so &
                   ( wrk3,psi%zwf(:,:,:,:,io,ik,im),ppg,mg%is_array,mg%ie_array,ik,info%icomm_r )
            else
              call calc_current_nonlocal_so &
                   ( wrk3,psi%zwf(:,:,:,:,io,ik,im),ppg,mg%is_array,mg%ie_array,ik )
            end if
          else
            if ( info%if_divide_rspace)then
              uVpsi(:) = uVpsibox2(ispin,io,ik,im,:)
              call calc_current_nonlocal_rdivided(wrk3,psi%zwf(:,:,:,ispin,io,ik,im),ppg,mg%is_array,mg%ie_array,ik,uVpsi)
            else
              call calc_current_nonlocal         (wrk3,psi%zwf(:,:,:,ispin,io,ik,im),ppg,mg%is_array,mg%ie_array,ik)
            end if
          end if
        else
          wrk3=0d0
        end if
        curr_wrk(:,ispin,io,ik) = (wrk1 + wrk2 + wrk3) / dble(ngrid) ! ngrid = aLxyz/Hxyz
      end do ! ispin
      
      if ( yn_spinorbit=='y' ) then
        curr_wrk(:,1,io,ik) = curr_wrk(:,1,io,ik) + curr_wrk(:,2,io,ik)
        curr_wrk(:,2,io,ik) = curr_wrk(:,1,io,ik)
      end if
      
    end do ! io
    end do ! ik
    
    call comm_summation(curr_wrk,curr_decomp,3*nspin*system%no*system%nk,info%icomm_rko)

    if (info%if_divide_rspace .and. yn_jm=='n' .and. .not. yn_spinorbit=='y') deallocate(uVpsibox,uVpsibox2,uVpsi)

    return
  end subroutine calc_current_decomposed

!===================================================================================================================================

  subroutine calc_microscopic_current(system,mg,stencil,info,psi,curr)
    use structures
    use communication, only: comm_summation
    use timer
    implicit none
    type(s_dft_system)   ,intent(in) :: system
    type(s_rgrid)        ,intent(in) :: mg
    type(s_stencil)      ,intent(in) :: stencil
    type(s_parallel_info),intent(in) :: info
    type(s_orbital)      ,intent(in) :: psi
    type(s_vector)                   :: curr ! electron number current density (without rho*A/c)
    !
    integer :: ispin,im,ik,io,is(3),ie(3),nsize,nspin,ix,iy,iz
    real(8) :: k(3),BT(3,3)
    real(8),allocatable :: wrk(:,:,:,:),wrk2(:,:,:,:)

    call timer_begin(LOG_MCURRENT_CALC)
    nspin = system%nspin

    if(.not. psi%update_zwf_overlap) stop "halo of wavefunction is not updated"
    if(info%im_s/=1 .or. info%im_e/=1) stop "error: im_s, im_e @ calc_microscopic_current"
    im = 1

    is = mg%is
    ie = mg%ie
    allocate(wrk(3,is(1):ie(1),is(2):ie(2),is(3):ie(3)),wrk2(3,is(1):ie(1),is(2):ie(2),is(3):ie(3)))
    nsize = 3* mg%num(1) * mg%num(2) * mg%num(3)

    ! transform from lattice-direction gradients to Cartesian (cf. calc_current);
    ! identity for an orthogonal, axis-aligned cell.
    BT = transpose(system%rmatrix_B)

    wrk2 = 0d0
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do ispin=1,nspin

      k(1:3) = system%vec_k(1:3,ik)
      call micro_current(mg%is_array,mg%ie_array,is,ie,mg%idx,mg%idy,mg%idz, &
      & stencil%coef_nab,k,BT,psi%zwf(:,:,:,ispin,io,ik,im),wrk)

!$omp parallel do collapse(2) private(ix,iy,iz)
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
        wrk2(1,ix,iy,iz) = wrk2(1,ix,iy,iz) + wrk(1,ix,iy,iz) * system%rocc(io,ik,ispin)*system%wtk(ik)
        wrk2(2,ix,iy,iz) = wrk2(2,ix,iy,iz) + wrk(2,ix,iy,iz) * system%rocc(io,ik,ispin)*system%wtk(ik)
        wrk2(3,ix,iy,iz) = wrk2(3,ix,iy,iz) + wrk(3,ix,iy,iz) * system%rocc(io,ik,ispin)*system%wtk(ik)
      end do
      end do
      end do

!     call nonlocal_part

    end do
    end do
    end do

    call timer_end(LOG_MCURRENT_CALC)

    call timer_begin(LOG_MCURRENT_COMM_COLL)
    call comm_summation(wrk2,curr%v,nsize,info%icomm_ko)
    call timer_end(LOG_MCURRENT_COMM_COLL)

    deallocate(wrk,wrk2)
    return

  contains
  
# define DX(dt) idx(ix+(dt)),iy,iz
# define DY(dt) ix,idy(iy+(dt)),iz
# define DZ(dt) ix,iy,idz(iz+(dt))

    subroutine micro_current(is_array,ie_array,is,ie,idx,idy,idz,nabt,k,BT,tpsi,jw)
      implicit none
      integer   ,intent(in) :: is_array(3),ie_array(3),is(3),ie(3), &
                             & idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)
      real(8)   ,intent(in) :: nabt(Nd,3),k(3),BT(3,3)
      complex(8),intent(in) :: tpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
      real(8)               :: jw(3,is(1):ie(1),is(2):ie(2),is(3):ie(3))
      !
      integer :: ix,iy,iz
      complex(8) :: xtmp,ytmp,ztmp,px,py,pz
      real(8) :: gx,gy,gz,rho_p
!$omp parallel do collapse(2) private(iz,iy,ix,xtmp,ytmp,ztmp,px,py,pz,gx,gy,gz,rho_p)
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)

!OCL swp
        do ix=is(1),ie(1)
          px = conjg(tpsi(ix,iy,iz))
          xtmp = nabt(1,1) * ( tpsi(DX(1)) - tpsi(DX(-1)) ) &
               + nabt(2,1) * ( tpsi(DX(2)) - tpsi(DX(-2)) ) &
               + nabt(3,1) * ( tpsi(DX(3)) - tpsi(DX(-3)) ) &
               + nabt(4,1) * ( tpsi(DX(4)) - tpsi(DX(-4)) )
          jw(1,ix,iy,iz) = aimag(px*xtmp)
        end do

!OCL swp
        do ix=is(1),ie(1)
          py = conjg(tpsi(ix,iy,iz))
          ytmp = nabt(1,2) * ( tpsi(DY(1)) - tpsi(DY(-1)) ) &
               + nabt(2,2) * ( tpsi(DY(2)) - tpsi(DY(-2)) ) &
               + nabt(3,2) * ( tpsi(DY(3)) - tpsi(DY(-3)) ) &
               + nabt(4,2) * ( tpsi(DY(4)) - tpsi(DY(-4)) )
          jw(2,ix,iy,iz) = aimag(py*ytmp)
        end do

!OCL swp
        do ix=is(1),ie(1)
          pz = conjg(tpsi(ix,iy,iz))
          ztmp = nabt(1,3) * ( tpsi(DZ(1)) - tpsi(DZ(-1)) ) &
               + nabt(2,3) * ( tpsi(DZ(2)) - tpsi(DZ(-2)) ) &
               + nabt(3,3) * ( tpsi(DZ(3)) - tpsi(DZ(-3)) ) &
               + nabt(4,3) * ( tpsi(DZ(4)) - tpsi(DZ(-4)) )
          jw(3,ix,iy,iz) = aimag(pz*ztmp)
        end do

        ! rotate the lattice-direction paramagnetic current to Cartesian via
        ! BT = transpose(rmatrix_B) (cf. calc_current), then add the Cartesian
        ! k-point term. BT is the identity for an orthogonal, axis-aligned cell.
!OCL swp
        do ix=is(1),ie(1)
          gx = jw(1,ix,iy,iz)
          gy = jw(2,ix,iy,iz)
          gz = jw(3,ix,iy,iz)
          rho_p = abs(tpsi(ix,iy,iz))**2
          jw(1,ix,iy,iz) = BT(1,1)*gx + BT(1,2)*gy + BT(1,3)*gz + k(1)*rho_p
          jw(2,ix,iy,iz) = BT(2,1)*gx + BT(2,2)*gy + BT(2,3)*gz + k(2)*rho_p
          jw(3,ix,iy,iz) = BT(3,1)*gx + BT(3,2)*gy + BT(3,3)*gz + k(3)*rho_p
        end do

      end do
      end do
!$omp end parallel do
      return
    end subroutine micro_current

  end subroutine calc_microscopic_current

end module
