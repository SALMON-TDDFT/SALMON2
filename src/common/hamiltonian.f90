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

module hamiltonian
  use nvtx_wrapper
  implicit none
  integer,private,parameter :: Nd = 4

contains

!===================================================================================================================================

SUBROUTINE hpsi(tpsi,htpsi,info,mg,V_local,system,stencil,srg,ppg,ttpsi)
  use structures
  use stencil_sub
  use nonlocal_potential
  use pseudo_pt_plusU_sub, only: pseudo_plusU, PLUS_U_ON
  use pseudo_pt_so_sub, only: pseudo_so
  use noncollinear_module, only: op_xc_noncollinear
  use sendrecv_grid, only: s_sendrecv_grid, update_overlap_real8, update_overlap_complex8
  use salmon_global, only: yn_want_communication_overlapping,yn_periodic,yn_jm,yn_symmetrized_stencil, &
          absorbing_boundary, yn_spinorbit
  use timer
  use code_optimization, only: stencil_is_parallelized_by_omp
  use communication, only: comm_summation
  implicit none

  external :: zstencil_typical_seq
  !$acc routine(zstencil_typical_seq) worker
  external :: zstencil_typical_gpu

  type(s_dft_system)   ,intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_rgrid)  ,intent(in) :: mg
  type(s_scalar) ,intent(in) :: V_local(system%Nspin)
  type(s_stencil),intent(in) :: stencil
  type(s_sendrecv_grid),intent(inout) :: srg
  type(s_pp_grid),intent(in) :: ppg
  type(s_orbital)            :: tpsi,htpsi
  type(s_orbital),optional   :: ttpsi
  !
  integer :: nspin,ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb,ix,iy,iz
  real(8) :: k_nabt(Nd,3),k_lap0,kAc(3)
  logical :: if_kAc,if_singlescale
  logical :: is_enable_overlapping
  !real(8) :: tmp,tmp1
  real(8) :: kAc0(3)

  call nvtxStartRange('hpsi', __LINE__)
  call timer_begin(LOG_UHPSI_ALL)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  nspin = system%nspin
  norb = Nspin* info%numo * info%numk * info%numm
  
  if_kAc = (yn_periodic=='y')
  if_singlescale = allocated(system%Ac_micro%v)

  ! check: can we execute computation/communication overlapping
  if (if_singlescale) then
    is_enable_overlapping = (yn_want_communication_overlapping == 'y') .and. &
                            stencil%if_orthogonal .and. &
                            info%if_divide_rspace
  else
    is_enable_overlapping = (yn_want_communication_overlapping == 'y') .and. &
                            stencil%if_orthogonal .and. &
                            info%if_divide_rspace .and. &
                            (im_e - im_s) <= 0    .and. &
                            (ik_e - ik_s) <= 0
  end if

  if(allocated(tpsi%rwf)) then

  ! overlap region communication
    call timer_begin(LOG_UHPSI_UPDATE_OVERLAP)
    if(info%if_divide_rspace) then
      call update_overlap_real8(srg, mg, tpsi%rwf)
    end if
    call timer_end(LOG_UHPSI_UPDATE_OVERLAP)

  ! stencil
    call nvtxStartRange('dstencil loop', __LINE__)
    call timer_begin(LOG_UHPSI_STENCIL)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
      call dstencil(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz &
                   ,tpsi%rwf(:,:,:,ispin,io,ik,im),htpsi%rwf(:,:,:,ispin,io,ik,im) &
                   ,V_local(ispin)%f,stencil%coef_lap0,stencil%coef_lap)
    end do
    end do
    end do
    end do
    call timer_end(LOG_UHPSI_STENCIL)
    call nvtxEndRange()

    ! nonlocal potential
    if ( yn_spinorbit=='y' ) then
      ! pseudopotential
      if(yn_jm=='n') call dpseudo(tpsi,htpsi,info,Nspin,ppg)
    else
      ! pseudopotential
      if(yn_jm=='n') call dpseudo(tpsi,htpsi,info,Nspin,ppg)
    end if

    ! DFT+U
    if ( PLUS_U_ON ) then
      call pseudo_plusU(tpsi,htpsi,system,info,ppg)
    end if

    ! meta-GGA: tau-operator part of Vxc
    if (system%xc_payload%use_tau_operator) then
      call add_xc_tau_operator(htpsi,tpsi,info,mg,system,stencil,srg)
    end if

  else

  ! overlap region communication
    call timer_begin(LOG_UHPSI_UPDATE_OVERLAP)
    if(info%if_divide_rspace .and. .not. is_enable_overlapping) then
      call update_overlap_complex8(srg, mg, tpsi%zwf)
    end if
    call timer_end(LOG_UHPSI_UPDATE_OVERLAP)

  ! stencil
    call timer_begin(LOG_UHPSI_STENCIL)
    if(stencil%if_orthogonal .and. .not.if_singlescale) then
    ! orthogonal lattice (general)
    
      if(stencil_is_parallelized_by_omp .or. is_enable_overlapping) then
      
        call nvtxStartRange('zstencil loop A', __LINE__)
        do im=im_s,im_e
        do ik=ik_s,ik_e
          if(if_kAc) then
            kAc(1:3) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
            k_lap0 = stencil%coef_lap0 + 0.5d0* sum(kAc(1:3)**2)
            k_nabt(:,1) = kAc(1) * stencil%coef_nab(:,1)
            k_nabt(:,2) = kAc(2) * stencil%coef_nab(:,2)
            k_nabt(:,3) = kAc(3) * stencil%coef_nab(:,3)
          else
            k_lap0 = stencil%coef_lap0
            k_nabt = 0d0
          end if
          if (is_enable_overlapping) then
            call zstencil_overlapped
          else
#ifdef USE_OPENACC
            call zstencil_typical_gpu(io_s, io_e, Nspin,mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz &
                          ,mg%is,mg%ie &
                          ,tpsi%zwf(:,:,:,:,:,ik,im),htpsi%zwf(:,:,:,:,:,ik,im) &
                          ,V_local(:),k_lap0,stencil%coef_lap,k_nabt &
                          )
#else
            do io=io_s,io_e
            do ispin=1,Nspin
              call zstencil(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz &
                            ,tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
                            ,V_local(ispin)%f,k_lap0,stencil%coef_lap,k_nabt)
            end do
            end do
#endif
          end if
        end do
        end do
      
        call nvtxEndRange()
      else
        call nvtxStartRange('zstencil loop B', __LINE__)
      ! OpenMP parallelization: k-point & orbital indices
      
#ifdef USE_OPENACC
        do im=im_s,im_e
        do ik=ik_s,ik_e
          if(if_kAc) then
            kAc(1:3) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
            k_lap0 = stencil%coef_lap0 + 0.5d0* sum(kAc(1:3)**2)
            k_nabt(:,1) = kAc(1) * stencil%coef_nab(:,1)
            k_nabt(:,2) = kAc(2) * stencil%coef_nab(:,2)
            k_nabt(:,3) = kAc(3) * stencil%coef_nab(:,3)
          else
            k_lap0 = stencil%coef_lap0
            k_nabt = 0d0
          end if
          call zstencil_typical_gpu(io_s, io_e, Nspin,mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz &
                          ,mg%is,mg%ie &
                          ,tpsi%zwf(:,:,:,:,:,ik,im),htpsi%zwf(:,:,:,:,:,ik,im) &
                          ,V_local(:),k_lap0,stencil%coef_lap,k_nabt &
                          )
        end do
        end do
#else
!$omp parallel do collapse(4) default(none) &
!$omp          private(im,ik,io,ispin,kAc,k_lap0,k_nabt) &
!$omp          shared(im_s,im_e,ik_s,ik_e,io_s,io_e,nspin,if_kac,system,stencil,mg,tpsi,htpsi,V_local)
        do im=im_s,im_e
        do ik=ik_s,ik_e
        do io=io_s,io_e
        do ispin=1,Nspin
          if(if_kAc) then
            kAc(1:3) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
            k_lap0 = stencil%coef_lap0 + 0.5d0* sum(kAc(1:3)**2)
            k_nabt(:,1) = kAc(1) * stencil%coef_nab(:,1)
            k_nabt(:,2) = kAc(2) * stencil%coef_nab(:,2)
            k_nabt(:,3) = kAc(3) * stencil%coef_nab(:,3)
          else
            k_lap0 = stencil%coef_lap0
            k_nabt = 0d0
          end if
          call zstencil(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz &
                            ,tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
                            ,V_local(ispin)%f,k_lap0,stencil%coef_lap,k_nabt)
        end do
        end do
        end do
        end do
!$omp end parallel do
#endif
        call nvtxEndRange()
        
      end if

      ! absorbing boundary condition
      if(absorbing_boundary=='z')then
         call add_imaginary_potential_for_absorbing_boundary_z(system,tpsi,htpsi)
      endif
      
    else if(stencil%if_orthogonal .and. if_singlescale) then
    ! orthogonal lattice, single-scale Maxwell-TDDFT
    
      if(yn_symmetrized_stencil=='y') then

!$omp parallel do collapse(4) default(none) &
!$omp          private(im,ik,io,ispin) &
!$omp          shared(im_s,im_e,ik_s,ik_e,io_s,io_e,nspin,mg,tpsi,htpsi,V_local,system,stencil)
        do im=im_s,im_e
        do ik=ik_s,ik_e
        do io=io_s,io_e
        do ispin=1,Nspin
          call zstencil_microAc_symmetrized(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz &
                        ,tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
                        ,V_local(ispin)%f,system%Ac_micro%v,stencil%coef_lap0 &
                        ,stencil%coef_lap,stencil%coef_nab,system%vec_k(1:3,ik))
        end do
        end do
        end do
        end do
!$omp end parallel do
        
      else if(stencil_is_parallelized_by_omp .or. is_enable_overlapping) then
        ! OpenMP parallelization: rgrid

        if (is_enable_overlapping) then
          call zstencil_microac_overlapped
        else
          do im=im_s,im_e
          do ik=ik_s,ik_e
          do io=io_s,io_e
          do ispin=1,Nspin
            call zstencil_microAc(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz &
                          ,tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
                          ,V_local(ispin)%f,system%Ac_micro%v,system%div_Ac%f,stencil%coef_lap0 &
                          ,stencil%coef_lap,stencil%coef_nab,system%vec_k(1:3,ik))
          end do
          end do
          end do
          end do
        end if

      else
        ! OpenMP parallelization: k-point & orbital indices

!$omp parallel do collapse(4) default(none) &
!$omp          private(im,ik,io,ispin) &
!$omp          shared(im_s,im_e,ik_s,ik_e,io_s,io_e,nspin,mg,tpsi,htpsi,V_local,system,stencil)
        do im=im_s,im_e
        do ik=ik_s,ik_e
        do io=io_s,io_e
        do ispin=1,Nspin
          call zstencil_microAc(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz &
                        ,tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
                        ,V_local(ispin)%f,system%Ac_micro%v,system%div_Ac%f,stencil%coef_lap0 &
                        ,stencil%coef_lap,stencil%coef_nab,system%vec_k(1:3,ik))
        end do
        end do
        end do
        end do
!$omp end parallel do

      end if

      ! absorbing boundary condition
      if(absorbing_boundary=='z')then
         call add_imaginary_potential_for_absorbing_boundary_z(system,tpsi,htpsi)
      endif


    else if(.not.stencil%if_orthogonal) then
    ! non-orthogonal lattice
    
#ifdef USE_OPENACC
!$acc update device(system%vec_Ac)
!$acc data copyin(tpsi,  &
!$acc             tpsi%zwf(mg%is_array(1):mg%ie_array(1),  &
!$acc                      mg%is_array(2):mg%ie_array(2),  &
!$acc                      mg%is_array(3):mg%ie_array(3),  &
!$acc                      1:Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e))  &
!$acc      copyout(htpsi,  &
!$acc              htpsi%zwf(mg%is_array(1):mg%ie_array(1),  &
!$acc                        mg%is_array(2):mg%ie_array(2),  &
!$acc                        mg%is_array(3):mg%ie_array(3),  &
!$acc                        1:Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e))
!$acc parallel present(system,mg,V_local,stencil,tpsi,htpsi)
!$acc loop collapse(3) private(kAc,kAc0,k_lap0) gang
#else
!$omp parallel do collapse(4) default(none) &
!$omp private(im,ik,io,ispin,kAc,k_lap0) &
!$omp shared(im_s,im_e,ik_s,ik_e,io_s,io_e,nspin,if_kac,system,stencil,mg,tpsi,htpsi,V_local)
#endif
      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin
        kAc = 0d0
        k_lap0 = 0d0
        if(if_kAc) then
          kAc(1:3) = system%vec_k(1:3,ik) + system%vec_Ac(1:3) ! Cartesian vector k+A/c
          k_lap0 = stencil%coef_lap0 + 0.5d0* sum(kAc(1:3)**2)
#ifdef USE_OPENACC
          kAc0 = kAc
          kAc(1) = system%rmatrix_B(1,1) * kAc0(1) + system%rmatrix_B(1,2) * kAc0(2) + system%rmatrix_B(1,3) * kAc0(3)
          kAc(2) = system%rmatrix_B(2,1) * kAc0(1) + system%rmatrix_B(2,2) * kAc0(2) + system%rmatrix_B(2,3) * kAc0(3)
          kAc(3) = system%rmatrix_B(3,1) * kAc0(1) + system%rmatrix_B(3,2) * kAc0(2) + system%rmatrix_B(3,3) * kAc0(3)
#else
          kAc(1:3) = matmul(system%rmatrix_B,kAc) ! B* (k+A/c)
#endif

        end if
          call zstencil_nonorthogonal(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz &
                                     ,tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
                                     ,V_local(ispin)%f,k_lap0,stencil%coef_lap,stencil%coef_nab,kAc,stencil%coef_F)
      end do
      end do
      end do
      end do
#ifdef USE_OPENACC
!$acc end parallel
!$acc end data
#else
!$omp end parallel do
#endif
      
    end if
    call timer_end(LOG_UHPSI_STENCIL)

  ! subtraction
    call nvtxStartRange('subtraction', __LINE__)
    call timer_begin(LOG_UHPSI_SUBTRACTION)
    if(present(ttpsi)) then
      if(allocated(tpsi%rwf)) then
#ifdef USE_OPENACC
        !$acc kernels loop private(im,ik,io,ispin,iz,iy,ix) collapse(7) independent
#else
        !$omp parallel do collapse(6) default(none) &
        !$omp          private(im,ik,io,ispin,iz,iy,ix) &
        !$omp          shared(im_s,im_e,ik_s,ik_e,io_s,io_e,nspin,mg) &
        !$omp          shared(ttpsi,htpsi,V_local,tpsi)
#endif
        do im=im_s,im_e
        do ik=ik_s,ik_e
        do io=io_s,io_e
        do ispin=1,Nspin
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            ttpsi%rwf(ix,iy,iz,ispin,io,ik,im) = htpsi%rwf(ix,iy,iz,ispin,io,ik,im) &
                                               - V_local(ispin)%f(ix,iy,iz) * tpsi%rwf(ix,iy,iz,ispin,io,ik,im)
          end do
          end do
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
      else
#ifdef USE_OPENACC
        !$acc kernels loop private(im,ik,io,ispin,iz,iy,ix) collapse(7) independent
#else
        !$omp parallel do collapse(6) default(none) &
        !$omp          private(im,ik,io,ispin,iz,iy,ix) &
        !$omp          shared(im_s,im_e,ik_s,ik_e,io_s,io_e,nspin,mg) &
        !$omp          shared(ttpsi,htpsi,V_local,tpsi)
#endif
        do im=im_s,im_e
        do ik=ik_s,ik_e
        do io=io_s,io_e
        do ispin=1,Nspin
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            ttpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) &
                                               - V_local(ispin)%f(ix,iy,iz) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
          end do
          end do
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
      end if
    end if
    call nvtxEndRange()
    call timer_end(LOG_UHPSI_SUBTRACTION)

  ! nonlocal potential
    call nvtxStartRange('nonlocal potential', __LINE__)
    if(yn_jm=='n') then
      if ( yn_spinorbit=='y' ) then
        call op_xc_noncollinear( tpsi, htpsi, info, mg )
        call pseudo_so(tpsi,htpsi,info,nspin,ppg,mg)
      else
      ! pseudopotential
        call zpseudo(tpsi,htpsi,info,nspin,ppg)
      end if
      if ( PLUS_U_ON ) then
        call pseudo_plusU(tpsi,htpsi,system,info,ppg)
      end if
    end if
    call nvtxEndRange()

    ! meta-GGA: tau-operator part of Vxc
    if (system%xc_payload%use_tau_operator) then
      call add_xc_tau_operator(htpsi,tpsi,info,mg,system,stencil,srg)
    end if

  end if

  call timer_end(LOG_UHPSI_ALL)
  call nvtxEndRange()

  return
contains
  subroutine zstencil_overlapped
    use sendrecv_grid, only: srg_pack, srg_communication, srg_unpack, &
                             update_overlap_complex8
    use code_optimization, only: modx,mody,modz,optimized_stencil_is_callable
    use communication, only: comm_proc_null, comm_get_groupinfo
    implicit none
    integer :: igs(3),ige(3)
    integer,parameter :: nyblk=8, nzblk=8
    integer :: ibx,iby,ibz
    integer :: iplane,ibs(3),ibe(3)
    integer :: is(3),ie(3)
    logical :: is_divided(3)
    integer :: myrank,nprocs

    call comm_get_groupinfo(srg%icomm, myrank, nprocs)
    do iplane=1,3
      is_divided(iplane) = srg%neig(1,iplane) /= comm_proc_null .and. &
                           srg%neig(1,iplane) /= myrank
    end do

    is(:) = mg%is(:)
    ie(:) = mg%ie(:)
    do iplane=1,3
      if (is_divided(iplane)) then
        is(iplane) = is(iplane) + 4
        ie(iplane) = ie(iplane) - 4
      end if
    end do

! phase 1. pack halo region
    call timer_begin(LOG_UHPSI_OVL_PHASE1)
    call update_overlap_complex8(srg, mg, tpsi%zwf, srg_pack)
    call timer_end  (LOG_UHPSI_OVL_PHASE1)

! phase 2. halo communication and computation without halo region
    call timer_begin(LOG_UHPSI_OVL_PHASE2)
!$omp parallel default(none) &
!$omp          private(io,ispin,igs,ige,ibx,iby,ibz) &
!$omp          shared(is,ie,ik,im,io_s,io_e,nspin,mg,tpsi,htpsi,V_local,k_lap0,stencil,k_nabt,srg,modx,mody,modz) &
!$omp          shared(optimized_stencil_is_callable)

! halo communication by master thread (tid = 0)
!$omp master
    call timer_begin(LOG_UHPSI_OVL_PHASE2_COMM)
    call update_overlap_complex8(srg, mg, tpsi%zwf, srg_communication)
    call timer_end  (LOG_UHPSI_OVL_PHASE2_COMM)
!$omp end master

! A computation with multi-thread except master thread,
! but master thread can join this loop if the communication completed before computation done.
!$omp do collapse(4) schedule(dynamic,1)
    do io=io_s,io_e
    do ispin=1,Nspin
    do ibz=is(3),ie(3),nzblk
    do iby=is(2),ie(2),nyblk
      igs(3) = ibz ; ige(3) = min(ibz + nzblk - 1, ie(3))
      igs(2) = iby ; ige(2) = min(iby + nyblk - 1, ie(2))
      igs(1) = is(1) ; ige(1) = ie(1)
#ifdef USE_OPT_EXPLICIT_VECTORIZATION
      if (optimized_stencil_is_callable) then
        call zstencil_tuned_seq(mg%is_array,mg%ie_array,mg%is,mg%ie,modx,mody,modz,igs,ige &
                               ,tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
                               ,V_local(ispin)%f,k_lap0,stencil%coef_lap,k_nabt)
      else
#endif
        call zstencil_typical_seq(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz,igs,ige &
                                 ,tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
                                 ,V_local(ispin)%f,k_lap0,stencil%coef_lap,k_nabt)
#ifdef USE_OPT_EXPLICIT_VECTORIZATION
      end if
#endif
    end do
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    call timer_end  (LOG_UHPSI_OVL_PHASE2)

! phase 3. unpack halo region
    call timer_begin(LOG_UHPSI_OVL_PHASE3)
    call update_overlap_complex8(srg, mg, tpsi%zwf, srg_unpack)
    call timer_end  (LOG_UHPSI_OVL_PHASE3)

    is(:) = mg%is(:)
    ie(:) = mg%ie(:)
! phase 4. computation with halo region
    call timer_begin(LOG_UHPSI_OVL_PHASE4)
!$omp parallel default(none) &
!$omp          firstprivate(is,ie) &
!$omp          private(io,ispin,iplane,igs,ige,ibx,iby,ibz,ibs,ibe) &
!$omp          shared(is_divided,ik,im,io_s,io_e,nspin,mg,tpsi,htpsi,V_local,k_lap0,stencil,k_nabt,srg,modx,mody,modz) &
!$omp          shared(optimized_stencil_is_callable)
    do iplane=1,6

      if (.not. is_divided((iplane+1)/2)) cycle

      ibs(:) = is(:)
      ibe(:) = ie(:)
      select case(iplane)
        case(1) ! update X (up)
          ibs(1) = mg%ie(1) - 4 + 1
          ibe(1) = mg%ie(1)
        case(2) ! update X (down)
          ibs(1) = mg%is(1)
          ibe(1) = mg%is(1) + 4 - 1
        case(3) ! update Y (up)
          ibs(2) = mg%ie(2) - 4 + 1
          ibe(2) = mg%ie(2)
        case(4) ! update Y (down)
          ibs(2) = mg%is(2)
          ibe(2) = mg%is(2) + 4 - 1
        case(5) ! update Z (up)
          ibs(3) = mg%ie(3) - 4 + 1
          ibe(3) = mg%ie(3)
        case(6) ! update Z (down)
          ibs(3) = mg%is(3)
          ibe(3) = mg%is(3) + 4 - 1
      end select

      select case(iplane)
        case(2)
          is(1) = is(1) + 4
          ie(1) = ie(1) - 4
        case(4)
          is(2) = is(2) + 4
          ie(2) = ie(2) - 4
      end select

#ifndef USE_OPENACC
!$omp do collapse(4) schedule(dynamic,1)
#endif
      do io=io_s,io_e
      do ispin=1,Nspin
      do ibz=ibs(3),ibe(3),nzblk
      do iby=ibs(2),ibe(2),nyblk
        igs(3) = ibz ; ige(3) = min(ibz + nzblk - 1, ibe(3))
        igs(2) = iby ; ige(2) = min(iby + nyblk - 1, ibe(2))
        igs(1) = ibs(1) ; ige(1) = ibe(1)
#ifdef USE_OPT_EXPLICIT_VECTORIZATION
        if (optimized_stencil_is_callable) then
          call zstencil_tuned_seq(mg%is_array,mg%ie_array,mg%is,mg%ie,modx,mody,modz,igs,ige &
                                 ,tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
                                 ,V_local(ispin)%f,k_lap0,stencil%coef_lap,k_nabt)
        else
#endif
          call zstencil_typical_seq(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz,igs,ige &
                                   ,tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
                                   ,V_local(ispin)%f,k_lap0,stencil%coef_lap,k_nabt)
#ifdef USE_OPT_EXPLICIT_VECTORIZATION
        end if
#endif
      end do
      end do
      end do
      end do
#ifndef USE_OPENACC
!$omp end do nowait
#endif
    end do
!$omp end parallel
    call timer_end  (LOG_UHPSI_OVL_PHASE4)
  end subroutine zstencil_overlapped

  subroutine zstencil_microac_overlapped
    use sendrecv_grid, only: srg_pack, srg_communication, srg_unpack, &
                             update_overlap_complex8
    use communication, only: comm_proc_null, comm_get_groupinfo
    implicit none
    integer :: igs(3),ige(3)
    integer,parameter :: nyblk=8, nzblk=8
    integer :: ibx,iby,ibz
    integer :: iplane,ibs(3),ibe(3)
    integer :: is(3),ie(3)
    logical :: is_divided(3)
    integer :: myrank,nprocs

    call comm_get_groupinfo(srg%icomm, myrank, nprocs)
    do iplane=1,3
      is_divided(iplane) = srg%neig(1,iplane) /= comm_proc_null .and. &
                           srg%neig(1,iplane) /= myrank
    end do

    is(:) = mg%is(:)
    ie(:) = mg%ie(:)
    do iplane=1,3
      if (is_divided(iplane)) then
        is(iplane) = is(iplane) + 4
        ie(iplane) = ie(iplane) - 4
      end if
    end do

! phase 1. pack halo region
    call timer_begin(LOG_UHPSI_OVL_PHASE1)
    call update_overlap_complex8(srg, mg, tpsi%zwf, srg_pack)
    call timer_end  (LOG_UHPSI_OVL_PHASE1)

! phase 2. halo communication and computation without halo region
    call timer_begin(LOG_UHPSI_OVL_PHASE2)
!$omp parallel default(none) &
!$omp          private(ik,im,io,ispin,igs,ige,ibx,iby,ibz) &
!$omp          shared(is,ie,im_s,im_e,ik_s,ik_e,io_s,io_e,nspin,mg,tpsi,htpsi,V_local,k_lap0,stencil,system,srg)

! halo communication by master thread (tid = 0)
!$omp master
    call timer_begin(LOG_UHPSI_OVL_PHASE2_COMM)
    call update_overlap_complex8(srg, mg, tpsi%zwf, srg_communication)
    call timer_end  (LOG_UHPSI_OVL_PHASE2_COMM)
!$omp end master

! A computation with multi-thread except master thread,
! but master thread can join this loop if the communication completed before computation done.
    do im=im_s,im_e
    do ik=ik_s,ik_e
!$omp do collapse(4) schedule(dynamic,1)
    do io=io_s,io_e
    do ispin=1,Nspin
    do ibz=is(3),ie(3),nzblk
    do iby=is(2),ie(2),nyblk
      igs(3) = ibz ; ige(3) = min(ibz + nzblk - 1, ie(3))
      igs(2) = iby ; ige(2) = min(iby + nyblk - 1, ie(2))
      igs(1) = is(1) ; ige(1) = ie(1)
      call zstencil_microAc_typical_seq( &
              mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz,igs,ige &
             ,tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
             ,V_local(ispin)%f,system%Ac_micro%v,system%div_Ac%f,stencil%coef_lap0 &
             ,stencil%coef_lap,stencil%coef_nab,system%vec_k(1:3,ik))
    end do
    end do
    end do
    end do
!$omp end do nowait
    end do
    end do
!$omp end parallel
    call timer_end  (LOG_UHPSI_OVL_PHASE2)

! phase 3. unpack halo region
    call timer_begin(LOG_UHPSI_OVL_PHASE3)
    call update_overlap_complex8(srg, mg, tpsi%zwf, srg_unpack)
    call timer_end  (LOG_UHPSI_OVL_PHASE3)

    is(:) = mg%is(:)
    ie(:) = mg%ie(:)
! phase 4. computation with halo region
    call timer_begin(LOG_UHPSI_OVL_PHASE4)
!$omp parallel default(none) &
!$omp          firstprivate(is,ie) &
!$omp          private(im,ik,io,ispin,iplane,igs,ige,ibx,iby,ibz,ibs,ibe) &
!$omp          shared(is_divided,im_s,im_e,ik_s,ik_e,io_s,io_e,nspin,mg,tpsi,htpsi,V_local,k_lap0,stencil,system,srg)
    do iplane=1,6

      if (.not. is_divided((iplane+1)/2)) cycle

      ibs(:) = is(:)
      ibe(:) = ie(:)
      select case(iplane)
        case(1) ! update X (up)
          ibs(1) = mg%ie(1) - 4 + 1
          ibe(1) = mg%ie(1)
        case(2) ! update X (down)
          ibs(1) = mg%is(1)
          ibe(1) = mg%is(1) + 4 - 1
        case(3) ! update Y (up)
          ibs(2) = mg%ie(2) - 4 + 1
          ibe(2) = mg%ie(2)
        case(4) ! update Y (down)
          ibs(2) = mg%is(2)
          ibe(2) = mg%is(2) + 4 - 1
        case(5) ! update Z (up)
          ibs(3) = mg%ie(3) - 4 + 1
          ibe(3) = mg%ie(3)
        case(6) ! update Z (down)
          ibs(3) = mg%is(3)
          ibe(3) = mg%is(3) + 4 - 1
      end select

      select case(iplane)
        case(2)
          is(1) = is(1) + 4
          ie(1) = ie(1) - 4
        case(4)
          is(2) = is(2) + 4
          ie(2) = ie(2) - 4
      end select

      do im=im_s,im_e
      do ik=ik_s,ik_e
#ifndef USE_OPENACC
!$omp do collapse(4) schedule(dynamic,1)
#endif
      do io=io_s,io_e
      do ispin=1,Nspin
      do ibz=ibs(3),ibe(3),nzblk
      do iby=ibs(2),ibe(2),nyblk
        igs(3) = ibz ; ige(3) = min(ibz + nzblk - 1, ibe(3))
        igs(2) = iby ; ige(2) = min(iby + nyblk - 1, ibe(2))
        igs(1) = ibs(1) ; ige(1) = ibe(1)
        call zstencil_microAc_typical_seq( &
                mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz,igs,ige &
               ,tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
               ,V_local(ispin)%f,system%Ac_micro%v,system%div_Ac%f,stencil%coef_lap0 &
               ,stencil%coef_lap,stencil%coef_nab,system%vec_k(1:3,ik))
      end do
      end do
      end do
      end do
#ifndef USE_OPENACC
!$omp end do nowait
#endif
      end do
      end do
    end do
!$omp end parallel
    call timer_end  (LOG_UHPSI_OVL_PHASE4)
  end subroutine zstencil_microac_overlapped

  subroutine add_imaginary_potential_for_absorbing_boundary_z(system,tpsi,htpsi)
    use structures, only: s_dft_system,s_orbital
    use salmon_global, only: al,imagnary_potential_w0,imagnary_potential_dr
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_orbital) :: tpsi,htpsi
    real(8) :: W0, dr, z,z0,z1,z2
    complex(8) :: W

    call nvtxStartRange('add_imaginary_potential_for_absorbing_boundary_z', __LINE__)
    dr = imagnary_potential_dr
    W0 = imagnary_potential_w0

!    z0 = lg%num(3) * system%hgs(3)   !cell length in z
    z0 = al(3)   !cell length in z
    z1 = dr      ! left boundary
    z2 = z0-dr   ! right boundary

    do im=im_s,im_e
!$omp parallel do collapse(4) &
!$omp          private(ik,io,ispin,ix,iy,iz,z,w)
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
       do iz = mg%is(3),mg%ie(3)
          z  = iz*system%hgs(3)
          if( z .le. z1 ) then
             w = dcmplx( 0d0, -w0*(z1-z)/dr )
          else if( z .ge. z2 ) then
             w = dcmplx( 0d0, -w0*(z-z2)/dr )
          else
             cycle
          endif

          do iy = mg%is(2),mg%ie(2)
          do ix = mg%is(1),mg%ie(1)
             htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = &
                   htpsi%zwf(ix,iy,iz,ispin,io,ik,im) &
                   + w * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
          end do
          end do
      end do
    end do
    end do
    end do
!$omp end parallel do
    end do
    call nvtxEndRange()

  end subroutine add_imaginary_potential_for_absorbing_boundary_z

end subroutine hpsi

!===================================================================================================================================

! Stop on a condition the tau operator cannot run under; every rank reaches it
! identically, so all finalize together.
subroutine fail_tau_operator(message)
  use communication, only: comm_is_root
  use parallelization, only: end_parallel, nproc_id_global
  implicit none
  character(*), intent(in) :: message

  if(comm_is_root(nproc_id_global)) then
    write(*,"(A)") 'Error in tau operator: '//trim(message)
  end if
  call end_parallel
  error stop 1
end subroutine fail_tau_operator

!===================================================================================================================================

! Generalized-Kohn-Sham meta-GGA kinetic-energy-density operator:
!   H_tau psi = -1/2 sum_c D_c(vtau * D_c psi), vtau = dE_xc/dtau,
!   D_c psi = (B^T nabla psi)_c + i*(k_c - u_c(r))*psi, u = j/n.
!
! tau_GI = 1/2 sum_occ f |D psi|^2 is stationary in u at u = j/n, so this D with u
! frozen is the exact derivative; D is anti-Hermitian, so the occupied-state sum of
! <psi|H_tau|psi> equals Integral vtau*tau_GI, which is payload%e_tau.
!
! GAUGE: the kinetic stencil carries k+A, and A shifts j/n by exactly A, so k-u and
! tau_GI keep their A-free values.  Building both from vec_k alone is then exact for a
! uniform A(t) and for the microscopic A(r,t) of single-scale Maxwell-TDDFT.
!
! In a ground state with time-reversal symmetry j = 0, hence u = 0 and every number here
! is what it was before u entered.
!
! Discrete divergence here matches the gradient used elsewhere to build tau; no effect on other paths.
!
! Two grid sweeps: pass 1 builds flux w_c = vtau*D_c psi; pass 2 applies the
! coef_nab divergence to the flux (exact adjoint of pass 1's gradient).
!
! Under real-space domain decomposition the divergence reads across the domain
! boundary, so the flux is halo-exchanged between the two passes, like the
! wavefunction's halo.
subroutine add_xc_tau_operator(htpsi,tpsi,info,mg,system,stencil,srg)
  use structures
  use salmon_global, only: yn_periodic
  use math_constants, only: zi
  use stencil_sub, only: calc_gradient_psi, calc_gradient_field
  use sendrecv_grid, only: update_overlap_real8, update_overlap_complex8
  implicit none
  type(s_orbital),            intent(inout) :: htpsi
  type(s_orbital),            intent(in)    :: tpsi
  type(s_parallel_info),      intent(in)    :: info
  type(s_rgrid),              intent(in)    :: mg
  type(s_dft_system),         intent(in)    :: system
  type(s_stencil),            intent(in)    :: stencil
  type(s_sendrecv_grid),      intent(inout) :: srg
  !
  integer :: im, ik, io, ispin, ix, iy, iz, c, d, n
  real(8) :: kvec(3), kvec_u(3), Bmat(3,3), vt, cnab(4,3), rwc(3), rdivg
  complex(8) :: psi0, wc(3), divg
  ! gpsi = (B^T nabla psi)_c ; the i*(k_c - u_c)*psi part of D_c psi is added in the
  ! grid loop below, where kvec_u = k - u(r).
  complex(8), allocatable :: gpsi(:,:,:,:)         ! (3, grid), one orbital at a time
  complex(8), allocatable :: uloc(:,:,:,:)         ! (grid, 3)                    no r-space division
  complex(8), allocatable :: uorb(:,:,:,:,:,:,:,:) ! (grid, nspin,io,ik,im, 3)    r-space division
  ! real (Gamma-only) path
  real(8),    allocatable :: rgpsi(:,:,:,:)
  real(8),    allocatable :: rgw(:,:,:,:)
  real(8),    allocatable :: rwvec(:,:,:,:)         ! (grid, 3)                    no r-space division
  real(8),    allocatable :: ruorb(:,:,:,:,:,:,:,:) ! (grid, nspin,io,ik,im, 3)    r-space division

#ifdef USE_OPENACC
  call fail_tau_operator("support is unavailable for OpenACC builds")
#endif
  if (.not. system%xc_payload%use_tau_operator) return
  if (.not. allocated(system%xc_payload%vtau%f)) then
    call fail_tau_operator("payload is enabled without a vtau field")
  end if
  if (.not. system%xc_payload%vtau_has_shadow_values) then
    call fail_tau_operator("requires vtau to be promoted to the halo grid layout")
  end if
  if (.not. allocated(system%xc_payload%uvel%v)) then
    call fail_tau_operator("payload is enabled without a velocity field")
  end if
  if (lbound(system%xc_payload%vtau%f,1) > mg%is_array(1) .or. &
      ubound(system%xc_payload%vtau%f,1) < mg%ie_array(1) .or. &
      lbound(system%xc_payload%vtau%f,2) > mg%is_array(2) .or. &
      ubound(system%xc_payload%vtau%f,2) < mg%ie_array(2) .or. &
      lbound(system%xc_payload%vtau%f,3) > mg%is_array(3) .or. &
      ubound(system%xc_payload%vtau%f,3) < mg%ie_array(3)) then
    call fail_tau_operator("requires vtau bounds compatible with mg%is_array:mg%ie_array")
  end if
  if (lbound(system%xc_payload%uvel%v,2) > mg%is(1) .or. &
      ubound(system%xc_payload%uvel%v,2) < mg%ie(1) .or. &
      lbound(system%xc_payload%uvel%v,3) > mg%is(2) .or. &
      ubound(system%xc_payload%uvel%v,3) < mg%ie(2) .or. &
      lbound(system%xc_payload%uvel%v,4) > mg%is(3) .or. &
      ubound(system%xc_payload%uvel%v,4) < mg%ie(3)) then
    call fail_tau_operator("requires uvel bounds compatible with mg%is:mg%ie")
  end if

  ! calc_gradient_psi and calc_gradient_field handle orthogonal and non-orthogonal
  ! lattices alike through rmatrix_B.
  Bmat = system%rmatrix_B
  cnab = stencil%coef_nab

  if (allocated(tpsi%rwf) .and. .not. info%if_divide_rspace) then
    ! Real orbitals: Gamma point only and j = 0 there, so u = 0 and D_c is the real
    ! Cartesian gradient.  calc_gradient_field carries B on both passes.
    allocate(rgpsi(3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(rgw  (3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(rwvec(  mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3),3))
    rwvec = 0d0
    do im=info%im_s,info%im_e
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do ispin=1,system%nspin
      call calc_gradient_field(mg,stencil%coef_nab,system%rmatrix_B, &
                               tpsi%rwf(:,:,:,ispin,io,ik,im),rgpsi)
!$omp parallel do collapse(2) private(ix,iy,iz,c)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        do c=1,3
          rwvec(ix,iy,iz,c) = system%xc_payload%vtau%f(mg%idx(ix),mg%idy(iy),mg%idz(iz)) * rgpsi(c,ix,iy,iz)
        end do
      end do
      end do
      end do
!$omp end parallel do
      do c=1,3
        call calc_gradient_field(mg,stencil%coef_nab,system%rmatrix_B,rwvec(:,:,:,c),rgw)
!$omp parallel do collapse(2) private(ix,iy,iz)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          ! accumulate -1/2 d/dx_c ( vtau d psi / d x_c )
          htpsi%rwf(ix,iy,iz,ispin,io,ik,im) = htpsi%rwf(ix,iy,iz,ispin,io,ik,im) - 0.5d0 * rgw(c,ix,iy,iz)
        end do
        end do
        end do
!$omp end parallel do
      end do
    end do
    end do
    end do
    end do
    deallocate(rgpsi,rgw,rwvec)
    return
  end if

  if (allocated(tpsi%rwf)) then
    ! Real orbitals under r-space division: same two passes as the complex branch,
    ! with u = 0 and the flux held for every local orbital so that one halo exchange
    ! per component feeds the divergence.
    allocate(rgpsi(3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(ruorb(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3), &
                   system%nspin,info%io_s:info%io_e,info%ik_s:info%ik_e,info%im_s:info%im_e,3))
    ruorb = 0d0
    do im=info%im_s,info%im_e
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do ispin=1,system%nspin
      call calc_gradient_field(mg,cnab,Bmat,tpsi%rwf(:,:,:,ispin,io,ik,im),rgpsi)
!$omp parallel do collapse(2) private(ix,iy,iz,c,d,vt,rwc)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        vt = system%xc_payload%vtau%f(mg%idx(ix),mg%idy(iy),mg%idz(iz))
        do c=1,3
          rwc(c) = vt * rgpsi(c,ix,iy,iz)
        end do
        do d=1,3
          ruorb(ix,iy,iz,ispin,io,ik,im,d) = Bmat(d,1)*rwc(1) + Bmat(d,2)*rwc(2) + Bmat(d,3)*rwc(3)
        end do
      end do
      end do
      end do
!$omp end parallel do
    end do
    end do
    end do
    end do
    do d=1,3
      call update_overlap_real8(srg, mg, ruorb(:,:,:,:,:,:,:,d))
    end do
    do im=info%im_s,info%im_e
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do ispin=1,system%nspin
!$omp parallel do collapse(2) private(ix,iy,iz,n,rdivg)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        rdivg = 0d0
        do n=1,4
          rdivg = rdivg + cnab(n,1)*( ruorb(mg%idx(ix+n),iy,iz,ispin,io,ik,im,1) - ruorb(mg%idx(ix-n),iy,iz,ispin,io,ik,im,1) ) &
                        + cnab(n,2)*( ruorb(ix,mg%idy(iy+n),iz,ispin,io,ik,im,2) - ruorb(ix,mg%idy(iy-n),iz,ispin,io,ik,im,2) ) &
                        + cnab(n,3)*( ruorb(ix,iy,mg%idz(iz+n),ispin,io,ik,im,3) - ruorb(ix,iy,mg%idz(iz-n),ispin,io,ik,im,3) )
        end do
        htpsi%rwf(ix,iy,iz,ispin,io,ik,im) = htpsi%rwf(ix,iy,iz,ispin,io,ik,im) - 0.5d0*rdivg
      end do
      end do
      end do
!$omp end parallel do
    end do
    end do
    end do
    end do
    deallocate(rgpsi,ruorb)
    return
  end if

  ! Complex orbitals.
  if (.not. info%if_divide_rspace) then
    ! No r-space division: build and consume the flux one orbital at a time,
    ! with no halo exchange (mg%idx/idy/idz already wrap periodically).
    allocate(gpsi(3,mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3)))
    allocate(uloc(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3),3))
    ! The halo layer is written by nothing below but IS read by the divergence.
    ! On a periodic grid mg%idx wraps back into the interior, but on an isolated
    ! one it addresses these cells for real, and they have to be zero.
    uloc = (0d0,0d0)
    do im=info%im_s,info%im_e
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do ispin=1,system%nspin
      ! vec_k alone; see the GAUGE note in this routine's header.
      kvec = 0d0
      if (yn_periodic == 'y') kvec(1:3) = system%vec_k(1:3,ik)
      call calc_gradient_psi(tpsi%zwf(:,:,:,ispin,io,ik,im),gpsi, &
           mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz, &
           cnab,Bmat)
!$omp parallel do collapse(2) private(ix,iy,iz,c,d,psi0,vt,wc,kvec_u)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        psi0 = tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
        vt   = system%xc_payload%vtau%f(mg%idx(ix),mg%idy(iy),mg%idz(iz))
        kvec_u(1:3) = kvec(1:3) - system%xc_payload%uvel%v(1:3,ix,iy,iz)
        do c=1,3
          wc(c) = vt * (gpsi(c,ix,iy,iz) + zi*kvec_u(c)*psi0)
        end do
        do d=1,3
          uloc(ix,iy,iz,d) = Bmat(d,1)*wc(1) + Bmat(d,2)*wc(2) + Bmat(d,3)*wc(3)
        end do
        htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) &
             - 0.5d0*zi*( kvec_u(1)*wc(1) + kvec_u(2)*wc(2) + kvec_u(3)*wc(3) )
      end do
      end do
      end do
!$omp end parallel do
!$omp parallel do collapse(2) private(ix,iy,iz,n,divg)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        divg = (0d0,0d0)
        do n=1,4
          divg = divg + cnab(n,1)*( uloc(mg%idx(ix+n),iy,iz,1) - uloc(mg%idx(ix-n),iy,iz,1) ) &
                      + cnab(n,2)*( uloc(ix,mg%idy(iy+n),iz,2) - uloc(ix,mg%idy(iy-n),iz,2) ) &
                      + cnab(n,3)*( uloc(ix,iy,mg%idz(iz+n),3) - uloc(ix,iy,mg%idz(iz-n),3) )
        end do
        htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) - 0.5d0*divg
      end do
      end do
      end do
!$omp end parallel do
    end do
    end do
    end do
    end do
    deallocate(gpsi,uloc)
    return
  end if

  ! r-space division: the divergence reads across the domain boundary, so the
  ! flux has to be halo-exchanged between the two passes.  It is built for every
  ! local orbital first so that the exchange is one batched call per component.
  allocate(gpsi(3,mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3)))
  allocate(uorb(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3), &
                system%nspin,info%io_s:info%io_e,info%ik_s:info%ik_e,info%im_s:info%im_e,3))
  uorb = (0d0,0d0)
  do im=info%im_s,info%im_e
  do ik=info%ik_s,info%ik_e
  do io=info%io_s,info%io_e
  do ispin=1,system%nspin
    ! vec_k alone; see the GAUGE note in this routine's header.
    kvec = 0d0
    if (yn_periodic == 'y') kvec(1:3) = system%vec_k(1:3,ik)
    call calc_gradient_psi(tpsi%zwf(:,:,:,ispin,io,ik,im),gpsi, &
         mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz, &
         cnab,Bmat)
!$omp parallel do collapse(2) private(ix,iy,iz,c,d,psi0,vt,wc,kvec_u)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      psi0 = tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
      vt   = system%xc_payload%vtau%f(mg%idx(ix),mg%idy(iy),mg%idz(iz))
      kvec_u(1:3) = kvec(1:3) - system%xc_payload%uvel%v(1:3,ix,iy,iz)
      do c=1,3
        wc(c) = vt * (gpsi(c,ix,iy,iz) + zi*kvec_u(c)*psi0)
      end do
      do d=1,3
        uorb(ix,iy,iz,ispin,io,ik,im,d) = Bmat(d,1)*wc(1) + Bmat(d,2)*wc(2) + Bmat(d,3)*wc(3)
      end do
      htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) &
           - 0.5d0*zi*( kvec_u(1)*wc(1) + kvec_u(2)*wc(2) + kvec_u(3)*wc(3) )
    end do
    end do
    end do
!$omp end parallel do
  end do
  end do
  end do
  end do
  do d=1,3
    call update_overlap_complex8(srg, mg, uorb(:,:,:,:,:,:,:,d))
  end do
  do im=info%im_s,info%im_e
  do ik=info%ik_s,info%ik_e
  do io=info%io_s,info%io_e
  do ispin=1,system%nspin
!$omp parallel do collapse(2) private(ix,iy,iz,n,divg)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      divg = (0d0,0d0)
      do n=1,4
        divg = divg + cnab(n,1)*( uorb(mg%idx(ix+n),iy,iz,ispin,io,ik,im,1) - uorb(mg%idx(ix-n),iy,iz,ispin,io,ik,im,1) ) &
                    + cnab(n,2)*( uorb(ix,mg%idy(iy+n),iz,ispin,io,ik,im,2) - uorb(ix,mg%idy(iy-n),iz,ispin,io,ik,im,2) ) &
                    + cnab(n,3)*( uorb(ix,iy,mg%idz(iz+n),ispin,io,ik,im,3) - uorb(ix,iy,mg%idz(iz-n),ispin,io,ik,im,3) )
      end do
      htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) - 0.5d0*divg
    end do
    end do
    end do
!$omp end parallel do
  end do
  end do
  end do
  end do
  deallocate(gpsi,uorb)

end subroutine add_xc_tau_operator

!===================================================================================================================================

subroutine update_vlocal(mg,nspin,Vh,Vpsl,Vxc,Vlocal)
  use structures
  use timer
  use nvtx_wrapper
  implicit none
  type(s_rgrid), intent(in) :: mg
  integer       ,intent(in) :: nspin
  type(s_scalar),intent(in) :: Vh,Vpsl,Vxc(nspin)
  type(s_scalar)            :: Vlocal(nspin)
  !
  integer :: is,ix,iy,iz
  call nvtxStartRange('update_vlocal', __LINE__)

  do is=1,nspin
#ifdef USE_OPENACC
!$acc parallel loop collapse(3) private(ix,iy,iz) copyin(vpsl, vh, vxc)
#else
!$omp parallel do collapse(2) private(ix,iy,iz)
#endif
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      Vlocal(is)%f(ix,iy,iz) = Vpsl%f(ix,iy,iz) + Vh%f(ix,iy,iz) + Vxc(is)%f(ix,iy,iz)
    end do
    end do
    end do
#ifdef USE_OPENACC
!$acc end parallel
#endif
  end do

  call nvtxEndRange
  return
end subroutine update_vlocal

!===================================================================================================================================

subroutine update_kvector_nonlocalpt(ik_s,ik_e,system,ppg)
  use math_constants,only : zi
  use structures
  use salmon_global, only: yn_spinorbit
  use update_kvector_so_sub, only: update_kvector_so
  use update_kvector_plusU_sub, only: update_kvector_plusU, PLUS_U_ON
  implicit none
  integer           ,intent(in) :: ik_s,ik_e !,n_max
  type(s_dft_system),intent(in) :: system
  type(s_pp_grid)               :: ppg
  !
  integer :: ilma,iatom,j,ik
  real(8) :: x,y,z
  complex(8) :: ekr
  real(8),allocatable :: kAc(:,:)
  
  allocate(kAc(3,ik_s:ik_e))
  do ik=ik_s,ik_e
    kAc(1:3,ik) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
  end do
  
  if ( yn_spinorbit=='y' ) then
    call update_kvector_so( ppg, kAc, ik_s, ik_e )
  end if
  if ( PLUS_U_ON ) then
    call update_kvector_plusU( ppg, kAc, ik_s, ik_e )
  end if
  
  if(.not.allocated(ppg%zekr_uV)) allocate(ppg%zekr_uV(ppg%nps,ppg%nlma,ik_s:ik_e))

#ifdef USE_OPENACC
!$acc kernels copyin(ppg)
!$acc loop collapse(2) gang private(ik,ilma,iatom,j,x,y,z,ekr)
#else
!$omp parallel do collapse(2) private(ik,ilma,iatom,j,x,y,z,ekr)
#endif
  do ik=ik_s,ik_e
    do ilma=1,ppg%nlma
      iatom = ppg%ia_tbl(ilma)
!$acc loop vector(128)
      do j=1,ppg%mps(iatom)
        x = ppg%rxyz(1,j,iatom)
        y = ppg%rxyz(2,j,iatom)
        z = ppg%rxyz(3,j,iatom)
        ekr = exp(zi*(kAc(1,ik)*x+kAc(2,ik)*y+kAc(3,ik)*z))
        ppg%zekr_uV(j,ilma,ik) = conjg(ekr) * ppg%uv(j,ilma)
      end do
    end do
  end do
#ifdef USE_OPENACC
!$acc end kernels
#else
!$omp end parallel do  
#endif

  deallocate(kAc)
  return
end subroutine update_kvector_nonlocalpt

subroutine update_kvector_nonlocalpt_microAc(ik_s,ik_e,system,ppg)
  use math_constants,only : zi
  use structures
  use timer
!  use fdtd_coulomb_gauge, only: line_integral
  implicit none
  integer           ,intent(in) :: ik_s,ik_e !,n_max
  type(s_dft_system),intent(in) :: system
!  type(s_rgrid)     ,intent(in) :: lg
!  real(8)           ,intent(in) :: vec_Ac(3,lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  type(s_pp_grid)               :: ppg
  !
  integer :: ilma,iatom,j,ik,ix,iy,iz
  real(8) :: kAc(3),x,y,z ! ,Hgs(3),r1_r0(3),r1(3),r0(3),integral
  complex(8) :: ekr
!  integer,allocatable :: index(:)
!  real(8),allocatable :: A_lerp(:,:),line(:,:),wrk(:)
!  Hgs = system%Hgs
!  allocate(A_lerp(3,n_max),line(3,n_max),wrk(n_max),index(n_max))

  call timer_begin(LOG_SS_UPDATE_NONLOCALPT_MICROAC)

  if(.not.allocated(ppg%zekr_uV)) allocate(ppg%zekr_uV(ppg%nps,ppg%nlma,ik_s:ik_e))
  do ilma=1,ppg%nlma
    iatom = ppg%ia_tbl(ilma)
    do j=1,ppg%mps(iatom)

!    ! C. Pickard & F. Mauri, PRL 91, 196401 (2003).
!      ix = ppg%jxyz(1,j,iatom) - lg%is(1) ! lg%is:lg%ie --> 0:ng%num-1
!      iy = ppg%jxyz(2,j,iatom) - lg%is(2)
!      iz = ppg%jxyz(3,j,iatom) - lg%is(3)
!      r1(1) = dble(ix)*hgs(1)
!      r1(2) = dble(iy)*hgs(2)
!      r1(3) = dble(iz)*hgs(3)
!      r1_r0 = ppg%rxyz(1:3,j,iatom) - system%Rion(1:3,iatom)
!      r0 = r1 - r1_r0
!      ! path: r0 --> r1 = (ix*hgs(1),iy*hgs(2),iz*hgs(3))
!      call line_integral(integral,r0,vec_Ac,lg%num(1),lg%num(2),lg%num(3),ix,iy,iz,Hgs(1),Hgs(2),Hgs(3) &
!            ,A_lerp,line,wrk,n_max,index)

      ix = ppg%jxyz(1,j,iatom)
      iy = ppg%jxyz(2,j,iatom)
      iz = ppg%jxyz(3,j,iatom)
      x = ppg%rxyz(1,j,iatom)
      y = ppg%rxyz(2,j,iatom)
      z = ppg%rxyz(3,j,iatom)
      do ik=ik_s,ik_e

!        k = system%vec_k(:,ik)
!        ekr = exp(zi*( k(1)*x+k(2)*y+k(3)*z + integral ))

      ! approximation: vector potential is almost constant in typical cutoff radius of pseudopotentials
        kAc = system%vec_k(:,ik) + system%Ac_micro%v(:,ix,iy,iz)
        ekr = exp(zi*( kAc(1)*x+kAc(2)*y+kAc(3)*z ))
        ppg%zekr_uV(j,ilma,ik) = conjg(ekr) * ppg%uv(j,ilma)
      end do
    end do
  end do
!  deallocate(A_lerp,line,wrk,index)

  call timer_end(LOG_SS_UPDATE_NONLOCALPT_MICROAC)

  return
end subroutine update_kvector_nonlocalpt_microAc

end module hamiltonian
