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
  implicit none
  integer,private,parameter :: Nd = 4

contains

subroutine fail_tau_operator(message)
  use communication, only: comm_is_root
  use parallelization, only: end_parallel, nproc_id_global
  implicit none
  character(*), intent(in) :: message

  if(comm_is_root(nproc_id_global)) then
    write(*,"(A)") 'Error in tau operator: '//trim(message)
  end if
  call end_parallel
  stop 1
end subroutine fail_tau_operator

!===================================================================================================================================

SUBROUTINE hpsi(tpsi,htpsi,info,mg,V_local,system,stencil,srg,ppg,ttpsi,xc_payload)
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
  type(s_xc_operator_payload), optional, intent(in) :: xc_payload
  !
  integer :: nspin,ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb,ix,iy,iz
  real(8) :: k_nabt(Nd,3),k_lap0,kAc(3)
  logical :: if_kAc,if_singlescale
  logical :: is_enable_overlapping
  !real(8) :: tmp,tmp1
  real(8) :: kAc0(3)

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

    if (present(xc_payload)) then
      if (xc_payload%use_tau_operator) then
        call add_xc_tau_operator(htpsi,tpsi,info,mg,system,stencil,srg,ppg,xc_payload)
      end if
    else if (system%xc_payload%use_tau_operator) then
      call add_xc_tau_operator(htpsi,tpsi,info,mg,system,stencil,srg,ppg,system%xc_payload)
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
      
      else
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
    call timer_begin(LOG_UHPSI_SUBTRACTION)
    if(present(ttpsi)) then
      if(allocated(tpsi%rwf)) then
        !$omp parallel do collapse(6) default(none) &
        !$omp          private(im,ik,io,ispin,iz,iy,ix) &
        !$omp          shared(im_s,im_e,ik_s,ik_e,io_s,io_e,nspin,mg) &
        !$omp          shared(ttpsi,htpsi,V_local,tpsi)
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
        !$omp end parallel do
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
    call timer_end(LOG_UHPSI_SUBTRACTION)

  ! nonlocal potential
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

    if (present(xc_payload)) then
      if (xc_payload%use_tau_operator) then
        call add_xc_tau_operator(htpsi,tpsi,info,mg,system,stencil,srg,ppg,xc_payload)
      end if
    else if (system%xc_payload%use_tau_operator) then
      call add_xc_tau_operator(htpsi,tpsi,info,mg,system,stencil,srg,ppg,system%xc_payload)
    end if

  end if

  call timer_end(LOG_UHPSI_ALL)

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

  end subroutine add_imaginary_potential_for_absorbing_boundary_z

end subroutine hpsi

!===================================================================================================================================

! The variational generalized-KS tau operator is -1/2 (nabla+ik).(vtau (nabla+ik) psi)
! for tau = 1/2 sum_occ |(nabla+ik)psi|^2 with vtau = dE/dtau (bare, as stored in the
! payload). The accumulated corr below builds the bare -(nabla+ik).(vtau (nabla+ik) psi),
! so the 1/2 is applied at the htpsi accumulation in every branch.
subroutine add_xc_tau_operator(htpsi,tpsi,info,mg,system,stencil,srg,ppg,xc_payload)
  use structures
  use salmon_global, only: yn_periodic
  use math_constants, only: zi
  implicit none
  type(s_orbital),            intent(inout) :: htpsi
  type(s_orbital),            intent(in)    :: tpsi
  type(s_parallel_info),      intent(in)    :: info
  type(s_rgrid),              intent(in)    :: mg
  type(s_dft_system),         intent(in)    :: system
  type(s_stencil),            intent(in)    :: stencil
  type(s_sendrecv_grid),      intent(inout) :: srg
  type(s_pp_grid),            intent(in)    :: ppg
  type(s_xc_operator_payload), optional, intent(in) :: xc_payload
  integer :: im, ik, io, ispin, ix, iy, iz, n
  integer :: ixp, ixm, iyp, iym, izp, izm
  real(8) :: a0, aplus, aminus, kvec(3), corr_r
  complex(8) :: psi0, corr, lap_axis, dprod_axis, dpsi_axis
  complex(8), allocatable :: wk_dpsi(:,:,:,:), wk_apsi(:,:,:), wk_dapsi(:,:,:,:), wk_adpsi(:,:,:,:)

#ifdef USE_OPENACC
  call fail_tau_operator("support is unavailable for OpenACC builds")
#endif
  if (present(xc_payload)) then
    if (.not. xc_payload%use_tau_operator) return
  else
    return
  end if
  if (.not. allocated(xc_payload%vtau%f)) then
    call fail_tau_operator("payload is enabled without vtau field")
  end if
  if (.not. xc_payload%vtau_has_shadow_values) then
    call fail_tau_operator("requires vtau shadow values to be prepared")
  end if
  if (allocated(system%Ac_micro%v)) then
    call fail_tau_operator("support is unavailable for single-scale Maxwell-TDDFT")
  end if
  if (lbound(xc_payload%vtau%f,1) > mg%is_array(1) .or. ubound(xc_payload%vtau%f,1) < mg%ie_array(1) .or. &
      lbound(xc_payload%vtau%f,2) > mg%is_array(2) .or. ubound(xc_payload%vtau%f,2) < mg%ie_array(2) .or. &
      lbound(xc_payload%vtau%f,3) > mg%is_array(3) .or. ubound(xc_payload%vtau%f,3) < mg%ie_array(3)) then
    call fail_tau_operator("requires vtau bounds compatible with mg%is_array:mg%ie_array")
  end if

  if (allocated(tpsi%rwf)) then
    if (.not. stencil%if_orthogonal) then
      call fail_tau_operator("nonorthogonal tau operator with real orbitals is unsupported")
    end if
!$omp parallel do collapse(4) default(none) &
!$omp          private(im,ik,io,ispin,iz,iy,ix,n,ixp,ixm,iyp,iym,izp,izm,a0,aplus,aminus,corr_r) &
!$omp          shared(info,mg,system,stencil,tpsi,htpsi,xc_payload)
    do im=info%im_s,info%im_e
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do ispin=1,system%nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        a0 = xc_payload%vtau%f(mg%idx(ix), mg%idy(iy), mg%idz(iz))
        corr_r = 0d0

        do n=1,4
          ixp = mg%idx(ix+n)
          ixm = mg%idx(ix-n)
          aplus = 0.5d0 * (a0 + xc_payload%vtau%f(ixp, mg%idy(iy), mg%idz(iz)))
          aminus = 0.5d0 * (a0 + xc_payload%vtau%f(ixm, mg%idy(iy), mg%idz(iz)))
          corr_r = corr_r - stencil%coef_lap(n,1) * ( &
                 aplus  * (tpsi%rwf(ixp,iy,iz,ispin,io,ik,im) - tpsi%rwf(ix,iy,iz,ispin,io,ik,im)) + &
                 aminus * (tpsi%rwf(ixm,iy,iz,ispin,io,ik,im) - tpsi%rwf(ix,iy,iz,ispin,io,ik,im)) )
        end do

        do n=1,4
          iyp = mg%idy(iy+n)
          iym = mg%idy(iy-n)
          aplus = 0.5d0 * (a0 + xc_payload%vtau%f(mg%idx(ix), iyp, mg%idz(iz)))
          aminus = 0.5d0 * (a0 + xc_payload%vtau%f(mg%idx(ix), iym, mg%idz(iz)))
          corr_r = corr_r - stencil%coef_lap(n,2) * ( &
                 aplus  * (tpsi%rwf(ix,iyp,iz,ispin,io,ik,im) - tpsi%rwf(ix,iy,iz,ispin,io,ik,im)) + &
                 aminus * (tpsi%rwf(ix,iym,iz,ispin,io,ik,im) - tpsi%rwf(ix,iy,iz,ispin,io,ik,im)) )
        end do

        do n=1,4
          izp = mg%idz(iz+n)
          izm = mg%idz(iz-n)
          aplus = 0.5d0 * (a0 + xc_payload%vtau%f(mg%idx(ix), mg%idy(iy), izp))
          aminus = 0.5d0 * (a0 + xc_payload%vtau%f(mg%idx(ix), mg%idy(iy), izm))
          corr_r = corr_r - stencil%coef_lap(n,3) * ( &
                 aplus  * (tpsi%rwf(ix,iy,izp,ispin,io,ik,im) - tpsi%rwf(ix,iy,iz,ispin,io,ik,im)) + &
                 aminus * (tpsi%rwf(ix,iy,izm,ispin,io,ik,im) - tpsi%rwf(ix,iy,iz,ispin,io,ik,im)) )
        end do

        htpsi%rwf(ix,iy,iz,ispin,io,ik,im) = htpsi%rwf(ix,iy,iz,ispin,io,ik,im) + 0.5d0 * corr_r
      end do
      end do
      end do
    end do
    end do
    end do
    end do
!$omp end parallel do
    return
  end if

  if (.not. stencil%if_orthogonal) then
    ! single work-array set reused across the (serial) orbital/k loop. The heavy
    ! grid loops inside add_xc_tau_operator_nonorth_complex are OMP-parallelized
    ! over the grid (mirrors calc_laplacian_field), so the orbital/k loop must NOT
    ! be collapsed for OMP here: doing so starves the threads when few orbitals/
    ! k-points are local per rank (e.g. SiO2 with ob/k decomposition, where
    ! ob/proc * k/proc < nthreads). Grid parallelism is unaffected by that.
    allocate(wk_dpsi (3,mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3)))
    allocate(wk_apsi (  mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3)))
    allocate(wk_dapsi(3,mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3)))
    allocate(wk_adpsi(3,mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3)))
    do im=info%im_s,info%im_e
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do ispin=1,system%nspin
      kvec = 0d0
      if (yn_periodic == 'y') kvec(1:3) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
      call add_xc_tau_operator_nonorth_complex( &
           htpsi%zwf(:,:,:,ispin,io,ik,im), tpsi%zwf(:,:,:,ispin,io,ik,im), &
           mg, stencil, system, xc_payload%vtau%f, kvec, &
           wk_dpsi, wk_apsi, wk_dapsi, wk_adpsi)
    end do
    end do
    end do
    end do
    deallocate(wk_dpsi,wk_apsi,wk_dapsi,wk_adpsi)
    return
  end if

!$omp parallel do collapse(4) default(none) &
!$omp          private(im,ik,io,ispin,iz,iy,ix,n,ixp,ixm,iyp,iym,izp,izm,a0,aplus,aminus,kvec,psi0,corr,lap_axis,dprod_axis,dpsi_axis) &
!$omp          shared(info,mg,system,stencil,tpsi,htpsi,xc_payload,yn_periodic)
  do im=info%im_s,info%im_e
  do ik=info%ik_s,info%ik_e
  do io=info%io_s,info%io_e
  do ispin=1,system%nspin
    kvec = 0d0
    if (yn_periodic == 'y') kvec(1:3) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      a0 = xc_payload%vtau%f(mg%idx(ix), mg%idy(iy), mg%idz(iz))
      psi0 = tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
      corr = 0d0

      lap_axis = 0d0
      dprod_axis = 0d0
      dpsi_axis = 0d0
      do n=1,4
        ixp = mg%idx(ix+n)
        ixm = mg%idx(ix-n)
        aplus = 0.5d0 * (a0 + xc_payload%vtau%f(ixp, mg%idy(iy), mg%idz(iz)))
        aminus = 0.5d0 * (a0 + xc_payload%vtau%f(ixm, mg%idy(iy), mg%idz(iz)))
        lap_axis = lap_axis + stencil%coef_lap(n,1) * ( &
                 aplus  * (tpsi%zwf(ixp,iy,iz,ispin,io,ik,im) - psi0) + &
                 aminus * (tpsi%zwf(ixm,iy,iz,ispin,io,ik,im) - psi0) )
        dprod_axis = dprod_axis + stencil%coef_nab(n,1) * ( &
                   xc_payload%vtau%f(ixp, mg%idy(iy), mg%idz(iz)) * tpsi%zwf(ixp,iy,iz,ispin,io,ik,im) - &
                   xc_payload%vtau%f(ixm, mg%idy(iy), mg%idz(iz)) * tpsi%zwf(ixm,iy,iz,ispin,io,ik,im) )
        dpsi_axis = dpsi_axis + stencil%coef_nab(n,1) * ( &
                  tpsi%zwf(ixp,iy,iz,ispin,io,ik,im) - tpsi%zwf(ixm,iy,iz,ispin,io,ik,im) )
      end do
      corr = corr - lap_axis - zi * kvec(1) * (dprod_axis + a0 * dpsi_axis) + a0 * kvec(1)**2 * psi0

      lap_axis = 0d0
      dprod_axis = 0d0
      dpsi_axis = 0d0
      do n=1,4
        iyp = mg%idy(iy+n)
        iym = mg%idy(iy-n)
        aplus = 0.5d0 * (a0 + xc_payload%vtau%f(mg%idx(ix), iyp, mg%idz(iz)))
        aminus = 0.5d0 * (a0 + xc_payload%vtau%f(mg%idx(ix), iym, mg%idz(iz)))
        lap_axis = lap_axis + stencil%coef_lap(n,2) * ( &
                 aplus  * (tpsi%zwf(ix,iyp,iz,ispin,io,ik,im) - psi0) + &
                 aminus * (tpsi%zwf(ix,iym,iz,ispin,io,ik,im) - psi0) )
        dprod_axis = dprod_axis + stencil%coef_nab(n,2) * ( &
                   xc_payload%vtau%f(mg%idx(ix), iyp, mg%idz(iz)) * tpsi%zwf(ix,iyp,iz,ispin,io,ik,im) - &
                   xc_payload%vtau%f(mg%idx(ix), iym, mg%idz(iz)) * tpsi%zwf(ix,iym,iz,ispin,io,ik,im) )
        dpsi_axis = dpsi_axis + stencil%coef_nab(n,2) * ( &
                  tpsi%zwf(ix,iyp,iz,ispin,io,ik,im) - tpsi%zwf(ix,iym,iz,ispin,io,ik,im) )
      end do
      corr = corr - lap_axis - zi * kvec(2) * (dprod_axis + a0 * dpsi_axis) + a0 * kvec(2)**2 * psi0

      lap_axis = 0d0
      dprod_axis = 0d0
      dpsi_axis = 0d0
      do n=1,4
        izp = mg%idz(iz+n)
        izm = mg%idz(iz-n)
        aplus = 0.5d0 * (a0 + xc_payload%vtau%f(mg%idx(ix), mg%idy(iy), izp))
        aminus = 0.5d0 * (a0 + xc_payload%vtau%f(mg%idx(ix), mg%idy(iy), izm))
        lap_axis = lap_axis + stencil%coef_lap(n,3) * ( &
                 aplus  * (tpsi%zwf(ix,iy,izp,ispin,io,ik,im) - psi0) + &
                 aminus * (tpsi%zwf(ix,iy,izm,ispin,io,ik,im) - psi0) )
        dprod_axis = dprod_axis + stencil%coef_nab(n,3) * ( &
                   xc_payload%vtau%f(mg%idx(ix), mg%idy(iy), izp) * tpsi%zwf(ix,iy,izp,ispin,io,ik,im) - &
                   xc_payload%vtau%f(mg%idx(ix), mg%idy(iy), izm) * tpsi%zwf(ix,iy,izm,ispin,io,ik,im) )
        dpsi_axis = dpsi_axis + stencil%coef_nab(n,3) * ( &
                  tpsi%zwf(ix,iy,izp,ispin,io,ik,im) - tpsi%zwf(ix,iy,izm,ispin,io,ik,im) )
      end do
      corr = corr - lap_axis - zi * kvec(3) * (dprod_axis + a0 * dpsi_axis) + a0 * kvec(3)**2 * psi0

      htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + 0.5d0 * corr
    end do
    end do
    end do
  end do
  end do
  end do
  end do
!$omp end parallel do
end subroutine add_xc_tau_operator

!===================================================================================================================================

subroutine add_xc_tau_operator_nonorth_complex(htpsi, tpsi, mg, stencil, system, vtau, kvec_cart, dpsi, apsi, dapsi, adpsi)
  use structures
  use math_constants, only: zi
  implicit none
  type(s_rgrid),      intent(in)    :: mg
  type(s_stencil),    intent(in)    :: stencil
  type(s_dft_system), intent(in)    :: system
  complex(8),         intent(inout) :: htpsi(mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), &
                                             mg%is_array(3):mg%ie_array(3))
  complex(8),         intent(in)    :: tpsi(mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), &
                                             mg%is_array(3):mg%ie_array(3))
  real(8),            intent(in)    :: vtau(mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), &
                                            mg%is_array(3):mg%ie_array(3))
  real(8),            intent(in)    :: kvec_cart(3)
  complex(8),         intent(inout) :: dpsi(3, mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), mg%is_array(3):mg%ie_array(3))
  complex(8),         intent(inout) :: apsi(mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), mg%is_array(3):mg%ie_array(3))
  complex(8),         intent(inout) :: dapsi(3, mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), mg%is_array(3):mg%ie_array(3))
  complex(8),         intent(inout) :: adpsi(3, mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), mg%is_array(3):mg%ie_array(3))
  integer :: ix, iy, iz, n, ixc, iyc, izc, jp, jm
  real(8) :: a0, k2, kvec_grid(3), ap, am, cf(6), lap(4,3), nab(4,3)
  complex(8) :: psi0, corr, d1, d2, d3, c21, c31, c12, c32, c13, c23
  integer :: idx(mg%is(1)-4:mg%ie(1)+4), idy(mg%is(2)-4:mg%ie(2)+4), idz(mg%is(3)-4:mg%ie(3)+4)

  ! local copies (idx as plain 1D arrays + stencil coefs) so the hot loop below
  ! vectorizes instead of repeatedly dereferencing derived-type components
  idx = mg%idx(mg%is(1)-4:mg%ie(1)+4)
  idy = mg%idy(mg%is(2)-4:mg%ie(2)+4)
  idz = mg%idz(mg%is(3)-4:mg%ie(3)+4)
  cf = stencil%coef_F; lap = stencil%coef_lap; nab = stencil%coef_nab

  call calc_tau_operator_axis_derivatives_complex(tpsi, mg, stencil%coef_nab, dpsi)

  apsi = 0d0
  adpsi = 0d0
!$omp parallel do collapse(2) default(none) private(ix,iy,iz,a0) shared(mg,vtau,tpsi,dpsi,apsi,adpsi)
  do iz = mg%is(3), mg%ie(3)
  do iy = mg%is(2), mg%ie(2)
  do ix = mg%is(1), mg%ie(1)
    a0 = vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz))
    apsi(ix,iy,iz) = a0 * tpsi(ix,iy,iz)
    adpsi(1,ix,iy,iz) = a0 * dpsi(1,ix,iy,iz)
    adpsi(2,ix,iy,iz) = a0 * dpsi(2,ix,iy,iz)
    adpsi(3,ix,iy,iz) = a0 * dpsi(3,ix,iy,iz)
  end do
  end do
  end do
!$omp end parallel do

  call calc_tau_operator_axis_derivatives_complex(apsi, mg, stencil%coef_nab, dapsi)

  kvec_grid = matmul(system%rmatrix_B, kvec_cart)
  k2 = sum(kvec_cart(1:3)**2)

  ! inlined direct(idir) + cross(iaxis,idir): the face-averaged vtau (ap/am) is
  ! computed once per face and reused for the diagonal (psi) and both off-diagonal
  ! (dpsi) terms. Same arithmetic/order as calc_tau_operator_{direct_axis,cross_component}.
!$omp parallel do collapse(2) default(none) &
!$omp   private(ix,iy,iz,ixc,iyc,izc,n,jp,jm,a0,psi0,ap,am,d1,d2,d3,c21,c31,c12,c32,c13,c23,corr) &
!$omp   shared(mg,vtau,tpsi,dpsi,dapsi,adpsi,htpsi,idx,idy,idz,cf,lap,nab,kvec_grid,k2)
  do iz = mg%is(3), mg%ie(3)
  do iy = mg%is(2), mg%ie(2)
    iyc = idy(iy); izc = idz(iz)
    do ix = mg%is(1), mg%ie(1)
      ixc = idx(ix)
      a0 = vtau(ixc, iyc, izc)
      psi0 = tpsi(ix,iy,iz)
      d1 = 0d0; d2 = 0d0; d3 = 0d0
      c21 = 0d0; c31 = 0d0; c12 = 0d0; c32 = 0d0; c13 = 0d0; c23 = 0d0
      do n = 1, 4
        jp = idx(ix+n); jm = idx(ix-n)
        ap = 0.5d0 * (a0 + vtau(jp, iyc, izc))
        am = 0.5d0 * (a0 + vtau(jm, iyc, izc))
        d1  = d1  + lap(n,1) * (ap * (tpsi(jp,iy,iz) - psi0) + am * (tpsi(jm,iy,iz) - psi0))
        c21 = c21 + nab(n,1) * (ap * dpsi(2,jp,iy,iz) - am * dpsi(2,jm,iy,iz))
        c31 = c31 + nab(n,1) * (ap * dpsi(3,jp,iy,iz) - am * dpsi(3,jm,iy,iz))
      end do
      do n = 1, 4
        jp = idy(iy+n); jm = idy(iy-n)
        ap = 0.5d0 * (a0 + vtau(ixc, jp, izc))
        am = 0.5d0 * (a0 + vtau(ixc, jm, izc))
        d2  = d2  + lap(n,2) * (ap * (tpsi(ix,jp,iz) - psi0) + am * (tpsi(ix,jm,iz) - psi0))
        c12 = c12 + nab(n,2) * (ap * dpsi(1,ix,jp,iz) - am * dpsi(1,ix,jm,iz))
        c32 = c32 + nab(n,2) * (ap * dpsi(3,ix,jp,iz) - am * dpsi(3,ix,jm,iz))
      end do
      do n = 1, 4
        jp = idz(iz+n); jm = idz(iz-n)
        ap = 0.5d0 * (a0 + vtau(ixc, iyc, jp))
        am = 0.5d0 * (a0 + vtau(ixc, iyc, jm))
        d3  = d3  + lap(n,3) * (ap * (tpsi(ix,iy,jp) - psi0) + am * (tpsi(ix,iy,jm) - psi0))
        c13 = c13 + nab(n,3) * (ap * dpsi(1,ix,iy,jp) - am * dpsi(1,ix,iy,jm))
        c23 = c23 + nab(n,3) * (ap * dpsi(2,ix,iy,jp) - am * dpsi(2,ix,iy,jm))
      end do
      corr = 0d0
      corr = corr - cf(1) * d1
      corr = corr - cf(2) * d2
      corr = corr - cf(3) * d3
      corr = corr - 0.5d0 * cf(4) * (c23 + c32)
      corr = corr - 0.5d0 * cf(5) * (c13 + c31)
      corr = corr - 0.5d0 * cf(6) * (c12 + c21)
      corr = corr - zi * ( kvec_grid(1) * (dapsi(1,ix,iy,iz) + adpsi(1,ix,iy,iz)) &
                         + kvec_grid(2) * (dapsi(2,ix,iy,iz) + adpsi(2,ix,iy,iz)) &
                         + kvec_grid(3) * (dapsi(3,ix,iy,iz) + adpsi(3,ix,iy,iz)) )
      corr = corr + a0 * k2 * psi0
      ! variational 1/2 of the generalized-KS tau operator (see add_xc_tau_operator)
      htpsi(ix,iy,iz) = htpsi(ix,iy,iz) + 0.5d0 * corr
    end do
  end do
  end do
!$omp end parallel do

end subroutine add_xc_tau_operator_nonorth_complex

subroutine calc_tau_operator_axis_derivatives_complex(box, mg, nabt, deriv)
  use structures
  implicit none
  type(s_rgrid), intent(in) :: mg
  real(8),       intent(in) :: nabt(4,3)
  complex(8),    intent(in) :: box(mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), &
                                   mg%is_array(3):mg%ie_array(3))
  complex(8), intent(out) :: deriv(3, mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), &
                                   mg%is_array(3):mg%ie_array(3))
  integer :: ix, iy, iz
  complex(8) :: w(3)

  deriv = 0d0
!$omp parallel do collapse(2) default(none) private(ix,iy,iz,w) shared(mg,nabt,box,deriv)
  do iz = mg%is(3), mg%ie(3)
  do iy = mg%is(2), mg%ie(2)
  do ix = mg%is(1), mg%ie(1)
    w(1) = nabt(1,1) * (box(mg%idx(ix+1), iy, iz) - box(mg%idx(ix-1), iy, iz)) &
         + nabt(2,1) * (box(mg%idx(ix+2), iy, iz) - box(mg%idx(ix-2), iy, iz)) &
         + nabt(3,1) * (box(mg%idx(ix+3), iy, iz) - box(mg%idx(ix-3), iy, iz)) &
         + nabt(4,1) * (box(mg%idx(ix+4), iy, iz) - box(mg%idx(ix-4), iy, iz))
    w(2) = nabt(1,2) * (box(ix, mg%idy(iy+1), iz) - box(ix, mg%idy(iy-1), iz)) &
         + nabt(2,2) * (box(ix, mg%idy(iy+2), iz) - box(ix, mg%idy(iy-2), iz)) &
         + nabt(3,2) * (box(ix, mg%idy(iy+3), iz) - box(ix, mg%idy(iy-3), iz)) &
         + nabt(4,2) * (box(ix, mg%idy(iy+4), iz) - box(ix, mg%idy(iy-4), iz))
    w(3) = nabt(1,3) * (box(ix, iy, mg%idz(iz+1)) - box(ix, iy, mg%idz(iz-1))) &
         + nabt(2,3) * (box(ix, iy, mg%idz(iz+2)) - box(ix, iy, mg%idz(iz-2))) &
         + nabt(3,3) * (box(ix, iy, mg%idz(iz+3)) - box(ix, iy, mg%idz(iz-3))) &
         + nabt(4,3) * (box(ix, iy, mg%idz(iz+4)) - box(ix, iy, mg%idz(iz-4)))
    deriv(:,ix,iy,iz) = w(:)
  end do
  end do
  end do
!$omp end parallel do
end subroutine calc_tau_operator_axis_derivatives_complex

pure function calc_tau_operator_direct_axis_complex(vtau, tpsi, mg, lapt, idir, ix, iy, iz) result(val)
  use structures
  implicit none
  type(s_rgrid), intent(in) :: mg
  real(8),       intent(in) :: vtau(mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), &
                                    mg%is_array(3):mg%ie_array(3))
  complex(8),    intent(in) :: tpsi(mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), &
                                    mg%is_array(3):mg%ie_array(3))
  real(8),       intent(in) :: lapt(4,3)
  integer,       intent(in) :: idir, ix, iy, iz
  complex(8) :: val
  complex(8) :: psi0
  real(8) :: a0, aplus, aminus

  val = 0d0
  psi0 = tpsi(ix,iy,iz)
  a0 = vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz))

  select case(idir)
  case(1)
    aplus = 0.5d0 * (a0 + vtau(mg%idx(ix+1), mg%idy(iy), mg%idz(iz)))
    aminus = 0.5d0 * (a0 + vtau(mg%idx(ix-1), mg%idy(iy), mg%idz(iz)))
    val = val + lapt(1,1) * (aplus * (tpsi(mg%idx(ix+1),iy,iz) - psi0) + aminus * (tpsi(mg%idx(ix-1),iy,iz) - psi0))
    aplus = 0.5d0 * (a0 + vtau(mg%idx(ix+2), mg%idy(iy), mg%idz(iz)))
    aminus = 0.5d0 * (a0 + vtau(mg%idx(ix-2), mg%idy(iy), mg%idz(iz)))
    val = val + lapt(2,1) * (aplus * (tpsi(mg%idx(ix+2),iy,iz) - psi0) + aminus * (tpsi(mg%idx(ix-2),iy,iz) - psi0))
    aplus = 0.5d0 * (a0 + vtau(mg%idx(ix+3), mg%idy(iy), mg%idz(iz)))
    aminus = 0.5d0 * (a0 + vtau(mg%idx(ix-3), mg%idy(iy), mg%idz(iz)))
    val = val + lapt(3,1) * (aplus * (tpsi(mg%idx(ix+3),iy,iz) - psi0) + aminus * (tpsi(mg%idx(ix-3),iy,iz) - psi0))
    aplus = 0.5d0 * (a0 + vtau(mg%idx(ix+4), mg%idy(iy), mg%idz(iz)))
    aminus = 0.5d0 * (a0 + vtau(mg%idx(ix-4), mg%idy(iy), mg%idz(iz)))
    val = val + lapt(4,1) * (aplus * (tpsi(mg%idx(ix+4),iy,iz) - psi0) + aminus * (tpsi(mg%idx(ix-4),iy,iz) - psi0))
  case(2)
    aplus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy+1), mg%idz(iz)))
    aminus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy-1), mg%idz(iz)))
    val = val + lapt(1,2) * (aplus * (tpsi(ix,mg%idy(iy+1),iz) - psi0) + aminus * (tpsi(ix,mg%idy(iy-1),iz) - psi0))
    aplus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy+2), mg%idz(iz)))
    aminus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy-2), mg%idz(iz)))
    val = val + lapt(2,2) * (aplus * (tpsi(ix,mg%idy(iy+2),iz) - psi0) + aminus * (tpsi(ix,mg%idy(iy-2),iz) - psi0))
    aplus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy+3), mg%idz(iz)))
    aminus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy-3), mg%idz(iz)))
    val = val + lapt(3,2) * (aplus * (tpsi(ix,mg%idy(iy+3),iz) - psi0) + aminus * (tpsi(ix,mg%idy(iy-3),iz) - psi0))
    aplus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy+4), mg%idz(iz)))
    aminus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy-4), mg%idz(iz)))
    val = val + lapt(4,2) * (aplus * (tpsi(ix,mg%idy(iy+4),iz) - psi0) + aminus * (tpsi(ix,mg%idy(iy-4),iz) - psi0))
  case(3)
    aplus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz+1)))
    aminus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz-1)))
    val = val + lapt(1,3) * (aplus * (tpsi(ix,iy,mg%idz(iz+1)) - psi0) + aminus * (tpsi(ix,iy,mg%idz(iz-1)) - psi0))
    aplus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz+2)))
    aminus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz-2)))
    val = val + lapt(2,3) * (aplus * (tpsi(ix,iy,mg%idz(iz+2)) - psi0) + aminus * (tpsi(ix,iy,mg%idz(iz-2)) - psi0))
    aplus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz+3)))
    aminus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz-3)))
    val = val + lapt(3,3) * (aplus * (tpsi(ix,iy,mg%idz(iz+3)) - psi0) + aminus * (tpsi(ix,iy,mg%idz(iz-3)) - psi0))
    aplus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz+4)))
    aminus = 0.5d0 * (a0 + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz-4)))
    val = val + lapt(4,3) * (aplus * (tpsi(ix,iy,mg%idz(iz+4)) - psi0) + aminus * (tpsi(ix,iy,mg%idz(iz-4)) - psi0))
  end select
end function calc_tau_operator_direct_axis_complex

pure function calc_tau_operator_cross_component_complex(vtau, dpsi, iaxis, mg, nabt, idir, ix, iy, iz) result(val)
  use structures
  implicit none
  type(s_rgrid), intent(in) :: mg
  real(8),       intent(in) :: vtau(mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), &
                                    mg%is_array(3):mg%ie_array(3))
  complex(8),    intent(in) :: dpsi(3, mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), &
                                    mg%is_array(3):mg%ie_array(3))
  real(8),       intent(in) :: nabt(4,3)
  integer,       intent(in) :: iaxis, idir, ix, iy, iz
  complex(8) :: val
  complex(8) :: flux_p, flux_m

  val = 0d0
  select case(idir)
  case(1)
    flux_p = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix+1), mg%idy(iy), mg%idz(iz))) * dpsi(iaxis, mg%idx(ix+1), iy, iz)
    flux_m = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix-1), mg%idy(iy), mg%idz(iz))) * dpsi(iaxis, mg%idx(ix-1), iy, iz)
    val = val + nabt(1,1) * (flux_p - flux_m)
    flux_p = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix+2), mg%idy(iy), mg%idz(iz))) * dpsi(iaxis, mg%idx(ix+2), iy, iz)
    flux_m = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix-2), mg%idy(iy), mg%idz(iz))) * dpsi(iaxis, mg%idx(ix-2), iy, iz)
    val = val + nabt(2,1) * (flux_p - flux_m)
    flux_p = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix+3), mg%idy(iy), mg%idz(iz))) * dpsi(iaxis, mg%idx(ix+3), iy, iz)
    flux_m = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix-3), mg%idy(iy), mg%idz(iz))) * dpsi(iaxis, mg%idx(ix-3), iy, iz)
    val = val + nabt(3,1) * (flux_p - flux_m)
    flux_p = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix+4), mg%idy(iy), mg%idz(iz))) * dpsi(iaxis, mg%idx(ix+4), iy, iz)
    flux_m = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix-4), mg%idy(iy), mg%idz(iz))) * dpsi(iaxis, mg%idx(ix-4), iy, iz)
    val = val + nabt(4,1) * (flux_p - flux_m)
  case(2)
    flux_p = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy+1), mg%idz(iz))) * dpsi(iaxis, ix, mg%idy(iy+1), iz)
    flux_m = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy-1), mg%idz(iz))) * dpsi(iaxis, ix, mg%idy(iy-1), iz)
    val = val + nabt(1,2) * (flux_p - flux_m)
    flux_p = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy+2), mg%idz(iz))) * dpsi(iaxis, ix, mg%idy(iy+2), iz)
    flux_m = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy-2), mg%idz(iz))) * dpsi(iaxis, ix, mg%idy(iy-2), iz)
    val = val + nabt(2,2) * (flux_p - flux_m)
    flux_p = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy+3), mg%idz(iz))) * dpsi(iaxis, ix, mg%idy(iy+3), iz)
    flux_m = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy-3), mg%idz(iz))) * dpsi(iaxis, ix, mg%idy(iy-3), iz)
    val = val + nabt(3,2) * (flux_p - flux_m)
    flux_p = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy+4), mg%idz(iz))) * dpsi(iaxis, ix, mg%idy(iy+4), iz)
    flux_m = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy-4), mg%idz(iz))) * dpsi(iaxis, ix, mg%idy(iy-4), iz)
    val = val + nabt(4,2) * (flux_p - flux_m)
  case(3)
    flux_p = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz+1))) * dpsi(iaxis, ix, iy, mg%idz(iz+1))
    flux_m = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz-1))) * dpsi(iaxis, ix, iy, mg%idz(iz-1))
    val = val + nabt(1,3) * (flux_p - flux_m)
    flux_p = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz+2))) * dpsi(iaxis, ix, iy, mg%idz(iz+2))
    flux_m = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz-2))) * dpsi(iaxis, ix, iy, mg%idz(iz-2))
    val = val + nabt(2,3) * (flux_p - flux_m)
    flux_p = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz+3))) * dpsi(iaxis, ix, iy, mg%idz(iz+3))
    flux_m = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz-3))) * dpsi(iaxis, ix, iy, mg%idz(iz-3))
    val = val + nabt(3,3) * (flux_p - flux_m)
    flux_p = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz+4))) * dpsi(iaxis, ix, iy, mg%idz(iz+4))
    flux_m = 0.5d0 * (vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz)) + vtau(mg%idx(ix), mg%idy(iy), mg%idz(iz-4))) * dpsi(iaxis, ix, iy, mg%idz(iz-4))
    val = val + nabt(4,3) * (flux_p - flux_m)
  end select
end function calc_tau_operator_cross_component_complex

subroutine update_vlocal(mg,nspin,Vh,Vpsl,Vxc,Vlocal)
  use structures
  use timer
  use nvtx
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
