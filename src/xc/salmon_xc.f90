!
!  Copyright 2018-2020 SALMON developers
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
!-----------------------------------------------------------------------------------------

#include "config.h"

#ifdef USE_LIBXC
#include "xc_version.h"
#endif

module salmon_xc
  use structures, only: s_xc_functional, s_xc_operator_payload, s_xc_aux_fields
  use builtin_pz, only: exc_cor_pz
  use builtin_pz_sp, only: exc_cor_pz_sp
  use builtin_pzm, only: exc_cor_pzm
  use builtin_r2scan, only: exc_cor_r2scan, grad_floor
  use builtin_tbmbj, only: exc_cor_tbmbj
  use builtin_pw, only: exc_cor_pw

#ifdef USE_LIBXC
#if XC_MAJOR_VERSION <= 4 
  use xc_f90_types_m
#endif
  use xc_f90_lib_m
#endif

  implicit none

! List of Exchange Correlation Functionals
  integer, parameter :: salmon_xctype_none  = 0
  integer, parameter :: salmon_xctype_pz    = 1
  integer, parameter :: salmon_xctype_pzm   = 2
  integer, parameter :: salmon_xctype_pbe   = 3
  integer, parameter :: salmon_xctype_tbmbj = 4
  integer, parameter :: salmon_xctype_tpss  = 5
  integer, parameter :: salmon_xctype_vs98  = 6
  integer, parameter :: salmon_xctype_pw    = 7
  integer, parameter :: salmon_xctype_r2scan = 8
#ifdef USE_LIBXC
  integer, parameter :: salmon_xctype_libxc = 101
#endif

  ! workspace used in exchange_correlation
  real(8),allocatable :: rho_tmp(:,:,:)
  real(8),allocatable :: rho_s_tmp(:,:,:,:)
  real(8),allocatable :: eexc_tmp(:,:,:)
  real(8),allocatable :: vxc_tmp(:,:,:)
  real(8),allocatable :: vxc_s_tmp(:,:,:,:)

  ! workspace used in exec_builtin_pz
  real(8),allocatable :: rho_s_1d(:)
  real(8),allocatable :: rho_s_sp_1d(:,:)
  real(8),allocatable :: exc_1d(:)
  real(8),allocatable :: eexc_1d(:)
  real(8),allocatable :: vexc_1d(:)
  real(8),allocatable :: vexc_sp_1d(:,:)

contains


! wrapper for calc_xc
  subroutine exchange_correlation(system, xc_func, mg, srg_scalar, srg, rho_s, pp, ppn, info, spsi, stencil, Vxc, E_xc, eexc, xc_payload, aux_override)
    use communication, only: comm_summation
    use structures
    use sendrecv_grid, only: update_overlap_real8
    use stencil_sub, only: calc_gradient_field, calc_laplacian_field
    use salmon_global, only: yn_spinorbit
    use noncollinear_module, only: rot_vxc_noncollinear
    use nvtx
    implicit none
    type(s_dft_system)      ,intent(inout), target :: system
    type(s_xc_functional)   ,intent(in) :: xc_func
    type(s_rgrid)           ,intent(in) :: mg
    type(s_sendrecv_grid)               :: srg_scalar, srg
    type(s_scalar)          ,intent(in) :: rho_s(system%nspin)
    type(s_pp_info)         ,intent(in) :: pp
    type(s_pp_nlcc)         ,intent(in) :: ppn
    type(s_parallel_info)   ,intent(in) :: info
    type(s_orbital)                     :: spsi
    type(s_stencil)         ,intent(in) :: stencil
    type(s_scalar)                      :: Vxc(system%nspin)
    real(8)                             :: E_xc
    type(s_scalar)          ,optional   :: eexc
    type(s_xc_operator_payload), intent(inout), optional, target :: xc_payload
    type(s_xc_aux_fields)   ,optional, intent(in) :: aux_override
    !
    integer :: ix,iy,iz,is,nspin,idir
    real(8) :: tot_exc
    ! real(8) :: rho_tmp(mg%num(1), mg%num(2), mg%num(3))
    ! real(8) :: rho_s_tmp(mg%num(1), mg%num(2), mg%num(3), 2)
    ! real(8) :: eexc_tmp(mg%num(1), mg%num(2), mg%num(3))
    ! real(8) :: vxc_tmp(mg%num(1), mg%num(2), mg%num(3))
    ! real(8) :: vxc_s_tmp(mg%num(1), mg%num(2), mg%num(3), 2)
    real(8),allocatable :: rhd(:,:,:), delr(:,:,:,:), grho(:,:,:,:), lrho(:,:,:), j(:,:,:,:), tau(:,:,:)
    real(8),allocatable :: rdedd_tmp(:,:,:,:),rdedd(:,:,:),drdedd_tmp(:,:,:,:),drdedd(:,:,:)
    real(8),allocatable :: delr_s(:,:,:,:,:),j_s(:,:,:,:,:),tau_s(:,:,:,:)
    real(8),allocatable :: rdedd_tmp_s(:,:,:,:,:),drdedd_s(:,:,:,:)
    type(s_xc_operator_payload), pointer :: payload
    logical :: use_payload
    
    call nvtxStartRange('exchange_correlation', __LINE__)

    use_payload = present(xc_payload) .or. (xc_func%xctype(1) == salmon_xctype_r2scan)
    nullify(payload)
    if (present(xc_payload)) then
      payload => xc_payload
    else if (xc_func%xctype(1) == salmon_xctype_r2scan) then
      payload => system%xc_payload
    end if

    if (use_payload) then
      call reset_xc_operator_payload(payload)
    else if (xc_operator_payload_is_active(system%xc_payload)) then
      call reset_xc_operator_payload(system%xc_payload)
    end if
    
    nspin = system%nspin

    if (nspin==1) then
      if (.not.allocated(rho_tmp)) allocate(rho_tmp(mg%num(1), mg%num(2), mg%num(3)))
      if (.not.allocated(vxc_tmp)) allocate(vxc_tmp(mg%num(1), mg%num(2), mg%num(3)))
    else if(nspin==2)then
      if (.not.allocated(rho_s_tmp)) allocate(rho_s_tmp(mg%num(1), mg%num(2), mg%num(3),2))
      if (.not.allocated(vxc_s_tmp)) allocate(vxc_s_tmp(mg%num(1), mg%num(2), mg%num(3),2))
    endif
    if (.not.allocated(eexc_tmp)) allocate(eexc_tmp(mg%num(1), mg%num(2), mg%num(3)))

    if(nspin==1)then
#ifdef USE_OPENACC
!$acc kernels loop collapse(3) private(iz,iy,ix)
#else
!$omp parallel do collapse(2) private(iz,iy,ix)
#endif
      do iz=1,mg%num(3)
      do iy=1,mg%num(2)
      do ix=1,mg%num(1)
        rho_tmp(ix,iy,iz)=rho_s(1)%f(mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)
      end do
      end do
      end do
#ifdef USE_OPENACC
!$acc end kernels
#else
!$omp end parallel do
#endif
    else if(nspin==2)then
#ifdef USE_OPENACC
!$acc kernels
!$acc loop collapse(4) private(is,iz,iy,ix)
#else
!$omp parallel private(is,iz,iy,ix)
#endif
      do is=1,2
#ifndef USE_OPENACC
!$omp do collapse(2)
#endif
      do iz=1,mg%num(3)
      do iy=1,mg%num(2)
      do ix=1,mg%num(1)
        rho_s_tmp(ix,iy,iz,is)=rho_s(is)%f(mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)
      end do
      end do
      end do
#ifndef USE_OPENACC
!$omp end do
#endif
      end do
#ifdef USE_OPENACC
!$acc end kernels
#else
!$omp end parallel
#endif
    end if

    if(xc_func%use_gradient) then ! meta GGA

      allocate (rhd (mg%is_array(1):mg%ie_array(1), &
                     mg%is_array(2):mg%ie_array(2), &
                     mg%is_array(3):mg%ie_array(3)))
      allocate (grho(3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)), &
              & lrho(mg%num(1), mg%num(2), mg%num(3)) )

       allocate (rdedd(mg%is_array(1):mg%ie_array(1), &
                         mg%is_array(2):mg%ie_array(2), &
                         mg%is_array(3):mg%ie_array(3)))
       allocate (drdedd_tmp(3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
       rhd = 0.d0
       grho = 0.d0
       lrho = 0.d0
       rdedd     =0.d0
       drdedd_tmp=0.d0
       if(nspin==1)then
         allocate (delr(mg%num(1), mg%num(2), mg%num(3) ,3))
         allocate (rdedd_tmp(mg%num(1), mg%num(2), mg%num(3),3), &
                   drdedd(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
         delr = 0.d0
         rdedd_tmp =0.d0
         drdedd    =0.d0
         if (xc_func%use_kinetic_energy) then
           allocate (j(mg%num(1), mg%num(2), mg%num(3) ,3), &
                     tau(mg%num(1), mg%num(2), mg%num(3)) )
           j = 0.d0
           tau = 0.d0
         end if
       elseif(nspin==2)then
         allocate (delr_s(mg%num(1), mg%num(2), mg%num(3),2,3))
         allocate (rdedd_tmp_s(mg%num(1), mg%num(2), mg%num(3),2,3), &
                   drdedd_s(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),2))
         delr_s = 0.d0 
         rdedd_tmp_s =0.d0
         drdedd_s    =0.d0
         if (xc_func%use_kinetic_energy) then
           allocate (j_s(mg%num(1), mg%num(2), mg%num(3),2,3), &
                     tau_s(mg%num(1), mg%num(2), mg%num(3),2) )
           j_s = 0.d0
           tau_s = 0.d0
         end if
       endif
    endif

    if(xc_func%use_gradient) then 
    if(nspin==1)then
!$omp parallel do collapse(2) private(ix,iy,iz)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        rhd(ix,iy,iz)=dble(rho_s(1)%f(ix,iy,iz))
      enddo
      enddo
      enddo
!$omp end parallel do

      if(info%if_divide_rspace) call update_overlap_real8(srg_scalar, mg, rhd)
      call calc_gradient_field(mg,stencil%coef_nab,system%rmatrix_B,rhd,grho)
      call calc_laplacian_field(mg,stencil,rhd,lrho(1:mg%num(1),1:mg%num(2),1:mg%num(3)))
      
!$omp parallel do collapse(2) private(iz,iy,ix)
      do iz=1,mg%num(3)
      do iy=1,mg%num(2)
      do ix=1,mg%num(1)
        delr(ix,iy,iz,1:3) = grho(1:3,mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)
      end do
      end do
      end do
!$omp end parallel do

      if (xc_func%use_kinetic_energy) then
        if (xc_func%xctype(1) == salmon_xctype_tbmbj) then
          if ((present(aux_override) .and. aux_override%use_tau) .and. (present(aux_override) .and. aux_override%use_j)) then
            continue
          else if (present(aux_override) .and. aux_override%use_tau) then
            call calc_tau_from_orbitals(system,mg,info,srg,stencil,spsi,rj=j)
          else if (present(aux_override) .and. aux_override%use_j) then
            call calc_tau_from_orbitals(system,mg,info,srg,stencil,spsi,tau=tau)
          else
            call calc_tau_from_orbitals(system,mg,info,srg,stencil,spsi,tau=tau,rj=j)
          end if

          if (present(aux_override) .and. aux_override%use_tau) then
!$omp parallel do collapse(2) private(iz,iy,ix)
            do iz=1,mg%num(3)
            do iy=1,mg%num(2)
            do ix=1,mg%num(1)
              tau(ix,iy,iz) = aux_override%tau%f(mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)
            end do
            end do
            end do
!$omp end parallel do
          end if

          if (present(aux_override) .and. aux_override%use_j) then
!$omp parallel do collapse(2) private(iz,iy,ix)
            do iz=1,mg%num(3)
            do iy=1,mg%num(2)
            do ix=1,mg%num(1)
              j(ix,iy,iz,1:3) = aux_override%j%v(1:3,mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)
            end do
            end do
            end do
!$omp end parallel do
          end if
        else if (present(aux_override) .and. aux_override%use_tau) then
!$omp parallel do collapse(2) private(iz,iy,ix)
          do iz=1,mg%num(3)
          do iy=1,mg%num(2)
          do ix=1,mg%num(1)
            tau(ix,iy,iz) = aux_override%tau%f(mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)
          end do
          end do
          end do
!$omp end parallel do
        else
          call calc_tau_from_orbitals(system,mg,info,srg,stencil,spsi,tau=tau)
        end if
      end if

    elseif(nspin==2)then
!$omp parallel do collapse(2) private(ix,iy,iz)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        rhd(ix,iy,iz)=dble(rho_s(1)%f(ix,iy,iz))
      enddo
      enddo
      enddo
!$omp end parallel do

      if(info%if_divide_rspace) call update_overlap_real8(srg_scalar, mg, rhd)
      call calc_gradient_field(mg,stencil%coef_nab,system%rmatrix_B,rhd,grho)
      call calc_laplacian_field(mg,stencil,rhd,lrho(1:mg%num(1),1:mg%num(2),1:mg%num(3)))

!$omp parallel do collapse(2) private(iz,iy,ix)
      do iz=1,mg%num(3)
      do iy=1,mg%num(2)
      do ix=1,mg%num(1)
        delr_s(ix,iy,iz,1,1:3) =grho(1:3,mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)
      end do
      end do
      end do
!$omp end parallel do  

!$omp parallel do collapse(2) private(ix,iy,iz)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        rhd (ix,iy,iz)=dble(rho_s(2)%f(ix,iy,iz))
      enddo
      enddo
      enddo
!$omp end parallel do

      if(info%if_divide_rspace) call update_overlap_real8(srg_scalar, mg, rhd)
      call calc_gradient_field(mg,stencil%coef_nab,system%rmatrix_B,rhd,grho)
      call calc_laplacian_field(mg,stencil,rhd,lrho(1:mg%num(1),1:mg%num(2),1:mg%num(3)))

!$omp parallel do collapse(2) private(iz,iy,ix)
      do iz=1,mg%num(3)
      do iy=1,mg%num(2)
      do ix=1,mg%num(1)
        delr_s(ix,iy,iz,2,1:3)=grho(1:3,mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)
      end do
      end do
      end do
!$omp end parallel do
  
    end if    
    end if

    if (xc_func%use_gradient) then
!      if(nspin==2) stop "error: GGA or metaGGA & spin/='unpolarized'"
      if(nspin==1)then
        if (xc_func%use_kinetic_energy) then
          if (xc_func%xctype(1) == salmon_xctype_tbmbj) then
            if (use_payload) then
              call calc_xc(xc_func, pp, rho=rho_tmp, eexc=eexc_tmp, vxc=vxc_tmp, rdedd=rdedd_tmp , grho=delr, &
                     &     rlrho=lrho, tau=tau, rj=j, rho_nlcc=ppn%rho_nlcc, payload=payload)
            else
              call calc_xc(xc_func, pp, rho=rho_tmp, eexc=eexc_tmp, vxc=vxc_tmp, rdedd=rdedd_tmp , grho=delr, &
                     &     rlrho=lrho, tau=tau, rj=j, rho_nlcc=ppn%rho_nlcc)
            end if
          else
            if (use_payload) then
              call calc_xc(xc_func, pp, rho=rho_tmp, eexc=eexc_tmp, vxc=vxc_tmp, rdedd=rdedd_tmp , grho=delr, &
                     &     rlrho=lrho, tau=tau, rho_nlcc=ppn%rho_nlcc, payload=payload)
            else
              call calc_xc(xc_func, pp, rho=rho_tmp, eexc=eexc_tmp, vxc=vxc_tmp, rdedd=rdedd_tmp , grho=delr, &
                     &     rlrho=lrho, tau=tau, rho_nlcc=ppn%rho_nlcc)
            end if
          end if
        else
          if (use_payload) then
            call calc_xc(xc_func, pp, rho=rho_tmp, eexc=eexc_tmp, vxc=vxc_tmp, rdedd=rdedd_tmp , grho=delr, &
                   &     rlrho=lrho, rho_nlcc=ppn%rho_nlcc, payload=payload)
          else
            call calc_xc(xc_func, pp, rho=rho_tmp, eexc=eexc_tmp, vxc=vxc_tmp, rdedd=rdedd_tmp , grho=delr, &
                   &     rlrho=lrho, rho_nlcc=ppn%rho_nlcc)
          end if
        end if
      elseif(nspin==2)then
!!!!!   Currently, only gga is working  !!!!!!!!!!!!!!!!!
        if (use_payload) then
          call calc_xc(xc_func, pp, rho_s=rho_s_tmp, grho_s=delr_s, &
           & eexc=eexc_tmp,vxc_s=vxc_s_tmp,rdedd_s=rdedd_tmp_s,rho_nlcc=ppn%rho_nlcc, payload=payload)
        else
          call calc_xc(xc_func, pp, rho_s=rho_s_tmp, grho_s=delr_s, &
           & eexc=eexc_tmp,vxc_s=vxc_s_tmp,rdedd_s=rdedd_tmp_s,rho_nlcc=ppn%rho_nlcc)
        end if
!               &     rlrho_s=lrho_s, tau_s=tau_s, rj_s=j_s, rho_nlcc=ppn%rho_nlcc)
      endif 
    else
      if(nspin==1)then
        if (use_payload) then
          call calc_xc(xc_func, pp, rho=rho_tmp, eexc=eexc_tmp, vxc=vxc_tmp, rho_nlcc=ppn%rho_nlcc, payload=payload)
        else
          call calc_xc(xc_func, pp, rho=rho_tmp, eexc=eexc_tmp, vxc=vxc_tmp, rho_nlcc=ppn%rho_nlcc)
        end if
      else if(nspin==2)then
        if (use_payload) then
          call calc_xc(xc_func, pp, rho_s=rho_s_tmp, eexc=eexc_tmp, vxc_s=vxc_s_tmp, rho_nlcc=ppn%rho_nlcc, payload=payload)
        else
          call calc_xc(xc_func, pp, rho_s=rho_s_tmp, eexc=eexc_tmp, vxc_s=vxc_s_tmp, rho_nlcc=ppn%rho_nlcc)
        end if
      end if
    end if

    if (use_payload) then
      call finalize_xc_payload(payload)
      if (xc_func%use_gradient .and. nspin == 1) then
        allocate(payload%grho%v(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
        payload%grho%v = grho
      end if
    end if

!!!!To include the sigma contribution to GGA Vxc potential !!!!!!!
    if (xc_func%use_gradient) then
    if(nspin==1)then
      drdedd=0.d0
      drdedd_tmp=0.d0
      if (use_payload) then
        allocate(payload%rdedd%v(3, mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), mg%is_array(3):mg%ie_array(3)))
        payload%rdedd%v = 0d0
      end if
      do idir=1,3 
        rdedd=0.d0
        do iz=1,mg%num(3)
        do iy=1,mg%num(2)
        do ix=1,mg%num(1)
           rdedd(mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)=rdedd_tmp(ix,iy,iz,idir)
        enddo
        enddo
        enddo
        if(info%if_divide_rspace) call update_overlap_real8(srg_scalar, mg, rdedd)
        if (use_payload) payload%rdedd%v(idir,:,:,:) = rdedd
        call calc_gradient_field(mg,stencil%coef_nab,system%rmatrix_B,rdedd,drdedd_tmp)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
           drdedd(ix,iy,iz)=drdedd(ix,iy,iz)+drdedd_tmp(idir,ix,iy,iz)
        enddo
        enddo
        enddo
      enddo
      if (use_payload) payload%rdedd_has_shadow_values = .true.
    elseif(nspin==2)then
      drdedd_s=0.d0
      drdedd_tmp=0.d0
      do is=1,2
        do idir=1,3
          rdedd=0.d0
          do iz=1,mg%num(3)
          do iy=1,mg%num(2)
          do ix=1,mg%num(1)
             rdedd(mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)=rdedd_tmp_s(ix,iy,iz,is,idir)
          enddo
          enddo
          enddo
          if(info%if_divide_rspace) call update_overlap_real8(srg_scalar, mg,rdedd)
          call calc_gradient_field(mg,stencil%coef_nab,system%rmatrix_B,rdedd,drdedd_tmp)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
             drdedd_s(ix,iy,iz,is)=drdedd_s(ix,iy,iz,is)+drdedd_tmp(idir,ix,iy,iz)
          enddo
          enddo
          enddo
        enddo
      enddo
    endif
    endif
!!! To include tau contribution to Vxc for MGGA functional  !!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(nspin==1)then
#ifdef USE_OPENACC
!$acc kernels loop collapse(3) private(iz,iy,ix)
#else
!$omp parallel do collapse(2) private(iz,iy,ix)
#endif
      do iz=1,mg%num(3)
      do iy=1,mg%num(2)
      do ix=1,mg%num(1)
        if (xc_func%use_gradient) then
          Vxc(1)%f(mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)=vxc_tmp(ix,iy,iz) & 
        &  +drdedd(mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)
        else
          Vxc(1)%f(mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)=vxc_tmp(ix,iy,iz)
        endif
      end do
      end do
      end do
#ifdef USE_OPENACC
!$acc end kernels
#else
!$omp end parallel do
#endif
    else if(nspin==2)then
#ifdef USE_OPENACC
!$acc kernels
!$acc loop collapse(4) private(is,iz,iy,ix)
#else
!$omp parallel private(is,iz,iy,ix)
#endif
      do is=1,2
#ifndef USE_OPENACC
!$omp do collapse(2)
#endif
      do iz=1,mg%num(3)
      do iy=1,mg%num(2)
      do ix=1,mg%num(1)
        if (xc_func%use_gradient) then
          Vxc(is)%f(mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)=vxc_s_tmp(ix,iy,iz,is) &
         &  +drdedd_s(mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1,is)
        else
          Vxc(is)%f(mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)=vxc_s_tmp(ix,iy,iz,is)
        endif
      end do
      end do
      end do
#ifndef USE_OPENACC
!$omp end do
#endif
      end do
#ifdef USE_OPENACC
!$acc end kernels
#else
!$omp end parallel
#endif
    end if

    tot_exc=0.d0
#ifdef USE_OPENACC
!$acc kernels loop collapse(3) reduction(+:tot_exc) private(iz,iy,ix)
#else
!$omp parallel do collapse(2) reduction(+:tot_exc) private(iz,iy,ix)
#endif
    do iz=1,mg%num(3)
    do iy=1,mg%num(2)
    do ix=1,mg%num(1)
      tot_exc=tot_exc+eexc_tmp(ix,iy,iz)
    end do
    end do
    end do
#ifdef USE_OPENACC
!$acc end kernels
#else
!$omp end parallel do
#endif
    tot_exc = tot_exc*system%hvol

    call comm_summation(tot_exc,E_xc,info%icomm_r)
    
    if(present(eexc)) then
      do iz=1,mg%num(3)
      do iy=1,mg%num(2)
      do ix=1,mg%num(1)
        eexc%f(mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1) = eexc_tmp(ix,iy,iz)
      end do
      end do
      end do    
    end if
    
    if(yn_spinorbit=='y') then
      call rot_vxc_noncollinear( Vxc, system, mg )
    end if

    call nvtxEndRange
    return
    
  contains

    subroutine finalize_xc_payload(xc_payload)
      use sendrecv_grid, only: update_overlap_real8
      implicit none
      type(s_xc_operator_payload), intent(inout) :: xc_payload
      real(8), allocatable :: vtau_local(:,:,:)
      real(8), allocatable :: vsigma_local(:,:,:)
      integer :: ix, iy, iz

      if (xc_payload%use_tau_operator) then
        if (.not. allocated(xc_payload%vtau%f)) then
          stop "error: tau operator payload is enabled without vtau field"
        end if
      end if
      if (.not. xc_payload%use_tau_operator .and. .not. allocated(xc_payload%vsigma%f)) return

      if (xc_payload%use_tau_operator) then
        allocate(vtau_local(mg%num(1), mg%num(2), mg%num(3)))
        vtau_local = xc_payload%vtau%f

        deallocate(xc_payload%vtau%f)
        allocate(xc_payload%vtau%f(mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), mg%is_array(3):mg%ie_array(3)))
        xc_payload%vtau%f = 0d0

        do iz = 1, mg%num(3)
        do iy = 1, mg%num(2)
        do ix = 1, mg%num(1)
          xc_payload%vtau%f(mg%is(1)+ix-1, mg%is(2)+iy-1, mg%is(3)+iz-1) = vtau_local(ix, iy, iz)
        end do
        end do
        end do

        if (info%if_divide_rspace) call update_overlap_real8(srg_scalar, mg, xc_payload%vtau%f)
        xc_payload%vtau_has_shadow_values = .true.
        deallocate(vtau_local)
      end if

      if (allocated(xc_payload%vsigma%f)) then
        allocate(vsigma_local(mg%num(1), mg%num(2), mg%num(3)))
        vsigma_local = xc_payload%vsigma%f

        deallocate(xc_payload%vsigma%f)
        allocate(xc_payload%vsigma%f(mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), mg%is_array(3):mg%ie_array(3)))
        xc_payload%vsigma%f = 0d0

        do iz = 1, mg%num(3)
        do iy = 1, mg%num(2)
        do ix = 1, mg%num(1)
          xc_payload%vsigma%f(mg%is(1)+ix-1, mg%is(2)+iy-1, mg%is(3)+iz-1) = vsigma_local(ix, iy, iz)
        end do
        end do
        end do

        if (info%if_divide_rspace) call update_overlap_real8(srg_scalar, mg, xc_payload%vsigma%f)
        xc_payload%vsigma_has_shadow_values = .true.
        deallocate(vsigma_local)
      end if
    end subroutine finalize_xc_payload

  end subroutine exchange_correlation

  subroutine calc_tau_from_orbitals(system, mg, info, srg, stencil, spsi, tau, rj)
    use communication, only: comm_summation
    use sendrecv_grid, only: update_overlap_real8, update_overlap_complex8
    use math_constants,only : zi
    use stencil_sub, only: calc_gradient_psi
    use structures, only: s_dft_system, s_rgrid, s_parallel_info, s_sendrecv_grid, s_stencil, s_orbital
    implicit none
    type(s_dft_system), intent(in) :: system
    type(s_rgrid), intent(in) :: mg
    type(s_parallel_info), intent(in) :: info
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_stencil), intent(in) :: stencil
    type(s_orbital), intent(inout) :: spsi
    real(8), intent(out), optional :: tau(mg%num(1),mg%num(2),mg%num(3))
    real(8), intent(out), optional :: rj(mg%num(1),mg%num(2),mg%num(3),3)
    integer :: im,ik,io,ispin,ix,iy,iz
    real(8) :: k(3),occ
    complex(8) :: zs(3),p
    logical :: need_tau, need_rj
    real(8), allocatable :: tau_tmp1(:,:,:), tau_tmp2(:,:,:)
    real(8), allocatable :: j_tmp1(:,:,:,:), j_tmp2(:,:,:,:)
    complex(8) :: gtpsi(3,mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3))

    if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ calc_tau_from_orbitals"
    im = 1
    need_tau = present(tau)
    need_rj = present(rj)
    if (.not. need_tau .and. .not. need_rj) stop "error: calc_tau_from_orbitals requires tau and/or rj"

    if (need_tau) then
      allocate(tau_tmp1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
      allocate(tau_tmp2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
      tau_tmp1 = 0d0
      tau_tmp2 = 0d0
    end if

    if (need_rj) then
      allocate(j_tmp1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),3))
      allocate(j_tmp2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),3))
      j_tmp1 = 0d0
      j_tmp2 = 0d0
    end if

    if(allocated(spsi%rwf)) then
       if(info%if_divide_rspace) call update_overlap_real8(srg, mg, spsi%rwf)
       if(.not.allocated(spsi%zwf)) &
            allocate(spsi%zwf(mg%is_array(1):mg%ie_array(1), &
                             mg%is_array(2):mg%ie_array(2), &
                             mg%is_array(3):mg%ie_array(3), &
                             system%nspin,info%io_s:info%io_e,info%ik_s:info%ik_e,info%im_s:info%im_e))
       spsi%zwf = dcmplx(spsi%rwf)
    else
       if(info%if_divide_rspace) call update_overlap_complex8(srg, mg, spsi%zwf)
    endif

    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do ispin=1,system%nspin
      call calc_gradient_psi(spsi%zwf(:,:,:,ispin,io,ik,im),gtpsi,mg%is_array,mg%ie_array,mg%is,mg%ie, &
           mg%idx,mg%idy,mg%idz,stencil%coef_nab,system%rmatrix_B)

      occ = system%rocc(io,ik,ispin)*system%wtk(ik)
      k(1:3) = system%vec_k(1:3,ik)
!$omp parallel do collapse(2) private(iz,iy,ix,zs,p)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        p = spsi%zwf(ix,iy,iz,ispin,io,ik,im)
        zs(1:3) = gtpsi(1:3,ix,iy,iz) + zi * k(1:3) * p
        if (need_tau) tau_tmp1(ix,iy,iz) = tau_tmp1(ix,iy,iz) + (abs(zs(1))**2+abs(zs(2))**2+abs(zs(3))**2)*occ*0.5d0
        if (need_rj) j_tmp1(ix,iy,iz,1:3) = j_tmp1(ix,iy,iz,1:3) + aimag(conjg(p)*zs(1:3))*occ
      end do
      end do
      end do
!$omp end parallel do
    end do
    end do
    end do

    if (need_tau) call comm_summation(tau_tmp1,tau_tmp2,mg%num(1)*mg%num(2)*mg%num(3),info%icomm_ko)
    if (need_rj) call comm_summation(j_tmp1,j_tmp2,mg%num(1)*mg%num(2)*mg%num(3)*3,info%icomm_ko)

!$omp parallel do collapse(2) private(iz,iy,ix)
    do iz=1,mg%num(3)
    do iy=1,mg%num(2)
    do ix=1,mg%num(1)
      if (need_tau) tau(ix,iy,iz) = tau_tmp2(mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1)
      if (need_rj) rj(ix,iy,iz,1:3) = j_tmp2(mg%is(1)+ix-1,mg%is(2)+iy-1,mg%is(3)+iz-1,1:3)
    end do
    end do
    end do
!$omp end parallel do

    if(allocated(spsi%rwf)) deallocate(spsi%zwf)
    if (need_rj) then
      deallocate(j_tmp2, j_tmp1)
    end if
    if (need_tau) deallocate(tau_tmp2, tau_tmp1)
  end subroutine calc_tau_from_orbitals


  subroutine reset_xc_operator_payload(payload)
    implicit none
    type(s_xc_operator_payload), intent(inout) :: payload

    payload%use_tau_operator = .false.
    payload%use_laplacian_operator = .false.
    payload%vtau_has_shadow_values = .false.
    payload%rdedd_has_shadow_values = .false.
    payload%vsigma_has_shadow_values = .false.
    if (allocated(payload%vtau%f)) deallocate(payload%vtau%f)
    if (allocated(payload%vsigma%f)) deallocate(payload%vsigma%f)
    if (allocated(payload%grho%v)) deallocate(payload%grho%v)
    if (allocated(payload%rdedd%v)) deallocate(payload%rdedd%v)
  end subroutine reset_xc_operator_payload

  logical function xc_operator_payload_is_active(payload)
    implicit none
    type(s_xc_operator_payload), intent(in) :: payload

    xc_operator_payload_is_active = payload%use_tau_operator .or. payload%use_laplacian_operator .or. &
         payload%vtau_has_shadow_values .or. payload%rdedd_has_shadow_values .or. payload%vsigma_has_shadow_values .or. &
         allocated(payload%vtau%f) .or. allocated(payload%vsigma%f) .or. &
         allocated(payload%grho%v) .or. allocated(payload%rdedd%v)
  end function xc_operator_payload_is_active


  subroutine print_xc_info()
    implicit none
    integer :: vmajor, vminor, vmicro

#ifdef USE_LIBXC
    call xc_f90_version(vmajor, vminor, vmicro)
    print '(2X,A,1X,I1,".",I1,".",I1)', "Libxc: [enabled]", vmajor, vminor, vmicro
#else
    print '(2X,A)', "Libxc: [disabled]"
#endif
    return
  end subroutine print_xc_info



  subroutine init_xc(xc, spin, cval, xcname, xname, cname)
    implicit none
    type(s_xc_functional), intent(inout) :: xc
    character(*), intent(in)           :: spin
    real(8), intent(in)                :: cval
    character(*), intent(in), optional :: xcname
    character(*), intent(in), optional :: xname
    character(*), intent(in), optional :: cname

    ! Initialization of xc variable
    xc%xctype(1:3) = salmon_xctype_none
    xc%cval = cval
    xc%use_gradient = .false.
    xc%use_laplacian = .false.
    xc%use_kinetic_energy = .false.
    xc%use_current = .false.
    
    if(spin=='unpolarized') then
      xc%ispin=0
    else !if(spin=='polarized') then
      xc%ispin=1
    end if

#ifdef USE_LIBXC
    if( xcname=="none" .and. xname=="none" .and. cname=="none" ) then
       print '(A, A)', "Error! Exchange nad Correction functionals are not specified!"
       stop
    endif
#endif         

    ! Exchange correlation
    if (present(xcname) .and. (len_trim(xcname) > 0)) then
      call setup_xcfunc(xcname)
    end if

    ! Exchange only
    if (present(xname) .and. (len_trim(xname) > 0)) then
      call setup_xfunc(xname)
    end if

    ! Correlation only
    if (present(cname) .and. (len_trim(cname) > 0)) then
      call setup_cfunc(cname)
    end if

    return

  contains



    subroutine setup_xcfunc(name)
      implicit none
      character(*), intent(in) :: name

      select case(lower(name))
      case('none')

#ifdef USE_LIBXC
#else
        print '(A, A)', "Error! Exchange functional is not specified!"
        stop
#endif         
        return
      
      case ('pz')
        xc%xctype(1) = salmon_xctype_pz
        return

      case ('pzm')
        xc%xctype(1) = salmon_xctype_pzm
        return
        
      case ('pw')
        xc%xctype(1) = salmon_xctype_pw
        return

      case ('pbe')
      
        xc%xctype(1) = salmon_xctype_pbe
        xc%use_gradient = .true.
        stop "Error: xc=pbe is not available. please use libxc_pbe."

      case ('tbmbj')

        xc%xctype(1) = salmon_xctype_tbmbj
        xc%use_gradient = .true.
        xc%use_laplacian = .true.
        xc%use_kinetic_energy = .true.
        xc%use_current = .true.
        return

      case ('bj_pw')

        xc%xctype(1) = salmon_xctype_tbmbj; xc%cval = 1d0
        xc%use_gradient = .true.
        xc%use_laplacian = .true.
        xc%use_kinetic_energy = .true.
        xc%use_current = .true.
        return

      case ('tpss')

        xc%xctype(1) = salmon_xctype_tpss
        xc%use_gradient = .true.
        xc%use_laplacian = .true.
        xc%use_kinetic_energy = .true.
        xc%use_current = .true.
        stop "TPSS functional is not implemented" ! future work

      case ('vs98')

        xc%xctype(1) = salmon_xctype_vs98
        xc%use_gradient = .true.
        xc%use_laplacian = .true.
        xc%use_kinetic_energy = .true.
        xc%use_current = .true.
        stop "VS98 functional is not implemented" ! future work

      case ('r2scan')

        xc%xctype(1) = salmon_xctype_r2scan
        xc%use_gradient = .true.
        xc%use_laplacian = .true.
        xc%use_kinetic_energy = .true.
        return

      ! Please insert additional functional here:
      ! e.g.
      ! case ('additional_functional')
      !   initialization_process_of_functional
      !   return

#ifdef USE_LIBXC
      case('libxc_pz')
        xc%xctype(2) = salmon_xctype_libxc
        xc%xctype(3) = salmon_xctype_libxc
        call init_libxc('LDA_X', 2)
        call init_libxc('LDA_C_PZ', 3)
        return

      case('libxc_pzm')
        xc%xctype(2) = salmon_xctype_libxc
        xc%xctype(3) = salmon_xctype_libxc
        call init_libxc('LDA_X', 2)
        call init_libxc('LDA_C_PZ_MOD', 3)
        return
        
      case('libxc_pw')
        xc%xctype(2) = salmon_xctype_libxc
        xc%xctype(3) = salmon_xctype_libxc
        call init_libxc('LDA_X', 2)
        call init_libxc('LDA_C_PW', 3)
        return

      case('libxc_pbe')
        xc%xctype(2) = salmon_xctype_libxc
        xc%xctype(3) = salmon_xctype_libxc
        call init_libxc('GGA_X_PBE', 2)
        call init_libxc('GGA_C_PBE', 3)
        return
#endif

      case default

#ifdef USE_LIBXC
        xc%xctype(1) = salmon_xctype_libxc
        call init_libxc(name, 1)
        return
#endif
        print '(A, A)', "Error! Undefined exchange functional:", trim(name)
        stop
      end select
      return
    end subroutine



    subroutine setup_xfunc(name)
      implicit none
      character(*), intent(in) :: name



      select case(name)
      case('none')
        ! xc%xctype(2) = salmon_xctype_none ! default
        return

      ! Please insert additional functional here:
      ! e.g.
      ! case ('additional_functional')
      !   initialization_process_of_functional
      !   return

      case default

#ifdef USE_LIBXC
        xc%xctype(2) = salmon_xctype_libxc
        call init_libxc(name, 2)
        return
#endif

        print '(A, A)', "Error! Undefined exchange functional:", trim(name)
        stop

      end select
      return
    end subroutine



    subroutine setup_cfunc(name)
      implicit none
      character(*), intent(in) :: name


      select case(name)
      case('none')
        ! xc%xctype(3) = salmon_xctype_none ! default
        return

      ! Please insert additional functional here:
      ! e.g.
      ! case ('additional_functional')
      !   initialization_process_of_functional
      !   return

      case default

#ifdef USE_LIBXC
        xc%xctype(3) = salmon_xctype_libxc
        call init_libxc(name, 3)
        return
#endif

        print '(A, A)', "Undefined correlation functional:", trim(name)
        stop

      end select
      return
    end subroutine



#ifdef USE_LIBXC
    subroutine init_libxc(libxc_name, ii)
      implicit none
      character(*), intent(in) :: libxc_name
      integer, intent(in) :: ii
      integer ::  ixc

      ixc = xc_f90_functional_get_number(trim(libxc_name))

      if (ixc <= 0) then
        print '(A, A)', "Undefined libxc functional:", trim(libxc_name)
        stop
      end if

!      if (spin /= 'unpolarized') then
!        print '(A)', "Spin polarized is not available"
!        stop
!      end if

      if (spin == 'unpolarized') then
#if XC_MAJOR_VERSION <= 4
        call xc_f90_func_init(xc%func(ii), xc%info(ii), ixc, XC_UNPOLARIZED)
#else
        call xc_f90_func_init(xc%func(ii), ixc, XC_UNPOLARIZED)
        xc%info(ii) = xc_f90_func_get_info(xc%func(ii))
#endif         
      else
#if XC_MAJOR_VERSION <= 4
        call xc_f90_func_init(xc%func(ii), xc%info(ii), ixc, XC_POLARIZED)
#else
        call xc_f90_func_init(xc%func(ii), ixc, XC_POLARIZED)
        xc%info(ii) = xc_f90_func_get_info(xc%func(ii))
#endif
      endif


#if XC_MAJOR_VERSION <= 4
      select case (xc_f90_info_family(xc%info(ii)))
#else
      select case (xc_f90_func_info_get_family(xc%info(ii)))
#endif
      case (XC_FAMILY_LDA)
      case (XC_FAMILY_GGA)
        xc%use_gradient = .true.
      case ( XC_FAMILY_MGGA)
        xc%use_gradient = .true.
        xc%use_laplacian = .true.
        xc%use_kinetic_energy = .true.
      case (XC_FAMILY_HYB_GGA, XC_FAMILY_HYB_MGGA)
        print '(A, A)', "Hybrid is not available:", trim(libxc_name)
        stop
      case default
        print '(A, A)', "Unknown Family:", trim(libxc_name)
        stop
      end select

      return
    end subroutine init_libxc
#endif

  end subroutine init_xc



  subroutine finalize_xc(xc)
    implicit none
    type(s_xc_functional), intent(inout) :: xc

#ifdef USE_LIBXC
    if (xc%xctype(1) == salmon_xctype_libxc) call xc_f90_func_end(xc%func(1))
    if (xc%xctype(2) == salmon_xctype_libxc) call xc_f90_func_end(xc%func(2))
    if (xc%xctype(3) == salmon_xctype_libxc) call xc_f90_func_end(xc%func(3))
#endif

    return
  end subroutine



  subroutine calc_xc(xc, pp, rho, rho_s, exc, eexc, vxc, vxc_s, rdedd, rdedd_s, &
      & grho, grho_s, rlrho, rlrho_s, tau, tau_s, rj, rj_s, &
      & rho_nlcc, payload, &
      & nd, ifdx, ifdy, ifdz, nabx, naby, nabz)
!      & nd, ifdx, ifdy, ifdz, nabx, naby, nabz, Hxyz, aLxyz)
    use structures, only: s_pp_info, s_xc_operator_payload
    use nvtx
    implicit none
    type(s_xc_functional), intent(in) :: xc
    type(s_pp_info),       intent(in) :: pp
    real(8), intent(in), optional :: rho(:, :, :) ! ispin = 0
    real(8), intent(in), optional :: rho_s(:, :, :, :) ! ispin = 1
    real(8), intent(out), optional :: exc(:, :, :) ! epsilon_xc[rho]
    real(8), intent(out), optional :: eexc(:, :, :) ! rho * epsilon_xc[rho]
    real(8), intent(out), optional :: vxc(:, :, :) ! v_xc[rho] for ispin=0
    real(8), intent(out), optional :: vxc_s(:, :, :, :) ! v_xc[rho] ispin=1
    !real(8), intent(out), optional :: gvxc(:, :, :) ! v_xc[rho] for ispin=0
    !real(8), intent(out), optional :: gvxc_s(:, :, :, :) ! v_xc[rho] ispin=1
    real(8), intent(in), optional :: grho(:, :, :, :)
    real(8), intent(in), optional :: grho_s(:, :, :, :, :) ! ispin = 1
    real(8), intent(in), optional :: rlrho(:, :, :)
    real(8), intent(in), optional :: rlrho_s(:, :, :, :) ! ispin = 1
    real(8), intent(in), optional :: rj(:, :, :, :)
    real(8), intent(in), optional :: rj_s(:, :, :, :) ! ispin = 1
    real(8), intent(in), optional :: tau(:, :, :)
    real(8), intent(in), optional :: tau_s(:, :, :, :) ! ispin = 1
!
    real(8), intent(in), optional :: rho_nlcc(:, :, :)
    real(8), intent(out),optional :: rdedd(:, :, :, :) 
    real(8), intent(out),optional :: rdedd_s(:, :, :, :, :)
    type(s_xc_operator_payload), intent(out), optional :: payload

    !===============================================================
    ! NOTE:
    !   The following section (finite difference table) is required
    !   in the GGA/mGGA functionals which is originally used by ARTED subroutine.
    ! TODO:
    !   Prepare more simplified solution to call the built-in GGA/mGGA functionals
    integer, intent(in), optional :: nd
    integer, intent(in), optional :: ifdx(:, :)
    integer, intent(in), optional :: ifdy(:, :)
    integer, intent(in), optional :: ifdz(:, :)
    real(8), intent(in), optional :: nabx(:)
    real(8), intent(in), optional :: naby(:)
    real(8), intent(in), optional :: nabz(:)
    ! real(8), intent(in), optional :: Hxyz
    ! real(8), intent(in), optional :: aLxyz
    !===============================================================

    integer :: nx, ny, nz, nl
    call nvtxStartRange('calc_xc', __LINE__)

    ! Detect size of 3-dimensional grid
    if (xc%ispin == 0) then
      nx = ubound(rho, 1) - lbound(rho, 1) + 1;
      ny = ubound(rho, 2) - lbound(rho, 2) + 1;
      nz = ubound(rho, 3) - lbound(rho, 3) + 1;
    else
      nx = ubound(rho_s, 1) - lbound(rho_s, 1) + 1;
      ny = ubound(rho_s, 2) - lbound(rho_s, 2) + 1;
      nz = ubound(rho_s, 3) - lbound(rho_s, 3) + 1;
    end if
    nl = nx * ny * nz

!    print *, nx,ny,nz,xc%ispin

    ! Initialize output variables
#ifdef USE_OPENACC
    if (present(exc)) then
!$acc kernels
      exc = 0d0
!$acc end kernels
    end if
    if (present(eexc)) then
!$acc kernels
      eexc = 0d0
!$acc end kernels
    end if
    if (present(vxc)) then
!$acc kernels
      vxc = 0d0
!$acc end kernels
    end if
    if (present(vxc_s)) then
!$acc kernels
      vxc_s = 0d0
!$acc end kernels
    end if
    if (present(rdedd)) then
!$acc kernels
      rdedd = 0.d0
!$acc end kernels
    end if
    if (present(rdedd_s)) then
!$acc kernels
      rdedd_s = 0.d0
!$acc end kernels
    end if
#else
    if (present(exc)) exc = 0d0
    if (present(eexc)) eexc = 0d0
    if (present(vxc)) vxc = 0d0
    if (present(vxc_s)) vxc_s = 0d0
    if (present(rdedd)) rdedd = 0.d0
    if (present(rdedd_s)) rdedd_s = 0.d0
#endif
    if (present(payload)) then
      payload%use_tau_operator = .false.
      payload%use_laplacian_operator = .false.
      payload%vtau_has_shadow_values = .false.
      if (allocated(payload%vtau%f)) deallocate(payload%vtau%f)
      if (allocated(payload%grho%v)) deallocate(payload%grho%v)
    end if

    ! Exchange-Correlation
    select case (xc%xctype(1))
    case(salmon_xctype_pz)
      call exec_builtin_pz()
    case(salmon_xctype_pzm)
      call exec_builtin_pzm()
    case(salmon_xctype_pw)
      call exec_builtin_pw()
    case(salmon_xctype_tbmbj)
      call exec_builtin_tbmbj()
    case(salmon_xctype_r2scan)
      call exec_builtin_r2scan()
#ifdef USE_LIBXC
    case(salmon_xctype_libxc)
     if (xc%ispin == 0) then 
       call exec_libxc(1)
     else
       call exec_libxc_spin(1)
     endif
#endif
    end select

    ! Exchange Only
    select case (xc%xctype(2))
#ifdef USE_LIBXC
    case(salmon_xctype_libxc)
      if (xc%ispin == 0) then
        call exec_libxc(2)
      else
        call exec_libxc_spin(2)
      endif
#endif
    end select

    ! Correlation Only
    select case (xc%xctype(3))
#ifdef USE_LIBXC
    case(salmon_xctype_libxc)
      if (xc%ispin == 0) then
        call exec_libxc(3)
      else
        call exec_libxc_spin(3)
      endif
#endif
    end select

    call nvtxEndRange
    return

  contains

    subroutine exec_builtin_calc_ax(y, alpha, x, n)
      implicit none
      integer                :: n
      real(8), intent(inout) :: y(n)
      real(8)                :: alpha
      real(8), intent(in)    :: x(n)
!$acc kernels
      y(:) = alpha * x(:)
!$acc end kernels
    end subroutine exec_builtin_calc_ax

    subroutine exec_builtin_calc_axpy(y, alpha, x, n)
      implicit none
      integer                :: n
      real(8), intent(inout) :: y(n)
      real(8)                :: alpha
      real(8), intent(in)    :: x(n)
!$acc kernels
      y(:) = alpha * x(:) + y(:)
!$acc end kernels
    end subroutine exec_builtin_calc_axpy

    subroutine exec_builtin_pz()
      use nvtx
      implicit none
      call nvtxStartRange('exec_builtin_pz', __LINE__)

      if (xc%ispin == 0) then
        if (.not.allocated(rho_s_1d)) allocate(rho_s_1d(nl))
        if (.not.allocated(vexc_1d)) allocate(vexc_1d(nl))
      else
        if (.not.allocated(rho_s_sp_1d)) allocate(rho_s_sp_1d(nl,2))
        if (.not.allocated(vexc_sp_1d)) allocate(vexc_sp_1d(nl,2))
      endif
      if (.not.allocated(exc_1d)) allocate(exc_1d(nl))
      if (.not.allocated(eexc_1d)) allocate(eexc_1d(nl))

#ifdef USE_OPENACC
      if (xc%ispin == 0) then
        call exec_builtin_calc_ax(rho_s_1d, 0.5d0, rho, nl)
      else if (xc%ispin == 1) then
        call exec_builtin_calc_ax(rho_s_sp_1d, 1.0d0, rho_s, nl*2)
      end if
#else
      if (xc%ispin == 0) then
        rho_s_1d = reshape(rho, (/nl/)) * 0.5
      else if (xc%ispin == 1) then
        rho_s_sp_1d = reshape(rho_s, (/nl,2/))
      end if
#endif

#ifndef SALMON_DEBUG_NEGLECT_NLCC
      if (present(rho_nlcc)) then
#ifdef USE_OPENACC
        if ( xc%ispin == 0 ) then
          call exec_builtin_calc_axpy(rho_s_1d, 0.5d0, rho_nlcc, nl)
        else if ( xc%ispin == 1 ) then
          call exec_builtin_calc_axpy(rho_s_sp_1d(:,1), 0.5d0, rho_nlcc, nl)
          call exec_builtin_calc_axpy(rho_s_sp_1d(:,2), 0.5d0, rho_nlcc, nl)
        end if
#else
        if ( xc%ispin == 0 ) then
          rho_s_1d = rho_s_1d + reshape(rho_nlcc, (/nl/)) * 0.5
        else if ( xc%ispin == 1 ) then
          rho_s_sp_1d(:,1) = rho_s_sp_1d(:,1) + reshape(rho_nlcc, (/nl/)) * 0.5
          rho_s_sp_1d(:,2) = rho_s_sp_1d(:,2) + reshape(rho_nlcc, (/nl/)) * 0.5
        end if
#endif
      endif
#endif

      if (xc%ispin == 0) then
        call exc_cor_pz(nl, rho_s_1d, exc_1d, eexc_1d, vexc_1d)
      else if (xc%ispin == 1) then
        call exc_cor_pz_sp(nl, rho_s_sp_1d, exc_1d, eexc_1d, vexc_sp_1d)
      end if

      if (xc%ispin == 0) then
        if (present(vxc)) then
#ifdef USE_OPENACC
          call exec_builtin_calc_axpy(vxc, 1.0d0, vexc_1d, nl)
#else
          vxc = vxc + reshape(vexc_1d, (/nx, ny, nz/))
#endif
        endif
      else if(xc%ispin == 1) then
        if (present(vxc_s)) then
#ifdef USE_OPENACC
          call exec_builtin_calc_axpy(vxc_s, 1.0d0, vexc_sp_1d, nl*2)
#else
          vxc_s = vxc_s + reshape(vexc_sp_1d, (/nx, ny, nz,2/))
#endif
        endif
      end if

      if (present(exc)) then
#ifdef USE_OPENACC
        call exec_builtin_calc_axpy(exc, 1.0d0, exc_1d, nl)
#else
        exc = exc + reshape(exc_1d, (/nx, ny, nz/))
#endif
      endif

      if (present(eexc)) then
#ifdef USE_OPENACC
        call exec_builtin_calc_axpy(eexc, 1.0d0, eexc_1d, nl)
#else
         eexc = eexc + reshape(eexc_1d, (/nx, ny, nz/))
#endif
      endif

      call nvtxEndRange
      return
    end subroutine exec_builtin_pz



    subroutine exec_builtin_pzm()
      use nvtx
      implicit none
      real(8) :: rho_s_1d(nl)
      real(8) :: exc_1d(nl)
      real(8) :: eexc_1d(nl)
      real(8) :: vexc_1d(nl)
      call nvtxStartRange('exec_builtin_pzm', __LINE__)

      rho_s_1d = reshape(rho, (/nl/)) * 0.5

#ifndef SALMON_DEBUG_NEGLECT_NLCC
      if (present(rho_nlcc)) then
        rho_s_1d = rho_s_1d + reshape(rho_nlcc, (/nl/)) * 0.5
      endif
#endif

      call exc_cor_pzm(nl, rho_s_1d, exc_1d, eexc_1d, vexc_1d)

      if (present(vxc)) then
         vxc = vxc + reshape(vexc_1d, (/nx, ny, nz/))
      endif

      if (present(exc)) then
         exc = exc + reshape(exc_1d, (/nx, ny, nz/))
      endif

      if (present(eexc)) then
         eexc = eexc + reshape(eexc_1d, (/nx, ny, nz/))
      endif

      call nvtxEndRange
      return
    end subroutine exec_builtin_pzm
    
    
    
    subroutine exec_builtin_pw()
      implicit none
      real(8) :: rho_s_1d(nl)
      real(8) :: rho_s_sp_1d(nl,2)
      real(8) :: exc_1d(nl)
      real(8) :: eexc_1d(nl)
      real(8) :: vexc_1d(nl)
      real(8) :: vexc_sp_1d(nl,2)

      if (xc%ispin == 0) then
        rho_s_1d = reshape(rho, (/nl/)) * 0.5
      else if (xc%ispin == 1) then
        rho_s_sp_1d = reshape(rho_s, (/nl,2/))
      end if

#ifndef SALMON_DEBUG_NEGLECT_NLCC
      if (present(rho_nlcc)) then
        if ( xc%ispin == 0 ) then
          rho_s_1d = rho_s_1d + reshape(rho_nlcc, (/nl/)) * 0.5
        else if ( xc%ispin == 1 ) then
          rho_s_sp_1d(:,1) = rho_s_sp_1d(:,1) + reshape(rho_nlcc, (/nl/)) * 0.5
          rho_s_sp_1d(:,2) = rho_s_sp_1d(:,2) + reshape(rho_nlcc, (/nl/)) * 0.5
        end if
      endif
#endif

      if (xc%ispin == 0) then
        call exc_cor_pw(nl, rho_s_1d, exc_1d, eexc_1d, vexc_1d)
      else if (xc%ispin == 1) then
        stop "PW for spin-polarized system is not supported."
      end if

      if (xc%ispin == 0) then
        if (present(vxc)) then
           vxc = vxc + reshape(vexc_1d, (/nx, ny, nz/))
        endif
      else if(xc%ispin == 1) then
        if (present(vxc_s)) then
           vxc_s = vxc_s + reshape(vexc_sp_1d, (/nx, ny, nz,2/))
        endif
      end if

      if (present(exc)) then
         exc = exc + reshape(exc_1d, (/nx, ny, nz/))
      endif

      if (present(eexc)) then
         eexc = eexc + reshape(eexc_1d, (/nx, ny, nz/))
      endif

      return
    end subroutine exec_builtin_pw
    

    subroutine exec_builtin_tbmbj()
      use nvtx
      implicit none
      real(8) :: rho_1d(nl)
      real(8) :: rho_s_1d(nl)
      real(8) :: grho_s_1d(nl, 3)
      real(8) :: rlrho_s_1d(nl)
      real(8) :: tau_s_1d(nl)
      real(8) :: j_s_1d(nl, 3)
      real(8) :: eexc_1d(nl)
      real(8) :: vexc_1d(nl)
      call nvtxStartRange('exec_builtin_tbmbj', __LINE__)

      rho_1d = reshape(rho, (/nl/))
      rho_s_1d = rho_1d * 0.5
#ifndef SALMON_DEBUG_NEGLECT_NLCC
      if (present(rho_nlcc)) then
        rho_s_1d = rho_s_1d + reshape(rho_nlcc, (/nl/)) * 0.5
      endif
#endif

      grho_s_1d = reshape(grho(:, :, :, :), (/nl, 3/)) * 0.5
      rlrho_s_1d = reshape(rlrho(:, :, :), (/nl/)) * 0.5
      tau_s_1d = reshape(tau(:, :, :), (/nl/)) * 0.5
      j_s_1d = reshape(rj(:, :, :, :), (/nl, 3/)) * 0.5

      !call exc_cor_tbmbj(nl, rho_1d, rho_s_1d,  grho_s_1d, rlrho_s_1d, tau_s_1d, j_s_1d, xc%cval, eexc_1d, vexc_1d, Hxyz, aLxyz)
      call exc_cor_tbmbj(nl, rho_1d, rho_s_1d,  grho_s_1d, rlrho_s_1d, tau_s_1d, j_s_1d, xc%cval, eexc_1d, vexc_1d)

      if (present(vxc)) then
        vxc = vxc + reshape(vexc_1d, (/nx, ny, nz/))
      endif

      if (present(exc)) then
        ! NOTE: Take care for "zero-division error"
        exc = exc + reshape(eexc_1d, (/nx, ny, nz/)) / rho
      endif

      if (present(eexc)) then
        eexc = eexc + reshape(eexc_1d, (/nx, ny, nz/))
      endif

      call nvtxEndRange
      return
    end subroutine exec_builtin_tbmbj


    subroutine exec_builtin_r2scan()
      implicit none
      real(8) :: rho_1d(nl)
      real(8) :: rho_s_1d(nl)
      real(8) :: tau_s_1d(nl)
      real(8) :: eexc_1d(nl)
      real(8) :: vexc_1d(nl)
      real(8) :: vtau_1d(nl)
      real(8) :: vgrad_1d(nl)
      real(8) :: vsigma_1d(nl)
      real(8) :: grho_norm_1d(nl)
      real(8) :: rho_safe
      integer :: i, ix, iy, iz, idir

      if (xc%ispin /= 0) stop "r2SCAN supports only nspin=1"

      rho_1d = reshape(rho, (/nl/))
#ifndef SALMON_DEBUG_NEGLECT_NLCC
      if (present(rho_nlcc)) then
        rho_1d = rho_1d + reshape(rho_nlcc, (/nl/))
      end if
#endif
      rho_s_1d = rho_1d * 0.5d0
      grho_norm_1d = reshape(sqrt(grho(:,:,:,1)**2 + grho(:,:,:,2)**2 + grho(:,:,:,3)**2), (/nl/))
      tau_s_1d = reshape(tau(:, :, :), (/nl/)) * 0.5d0

      call exc_cor_r2scan(nl, rho_1d, rho_s_1d, grho_norm_1d, tau_s_1d, eexc_1d, vexc_1d, vtau_1d, vgrad_1d)

      if (present(vxc)) then
        vxc = vxc + reshape(vexc_1d, (/nx, ny, nz/))
      end if

      if (present(exc)) then
!$omp parallel do collapse(2) default(none) private(ix, iy, iz, i, rho_safe) &
!$omp shared(nx, ny, nz, rho_1d, eexc_1d, exc)
        do iz = 1, nz
        do iy = 1, ny
        do ix = 1, nx
          i = (iz-1) * nx * ny + (iy-1) * nx + ix
          rho_safe = max(rho_1d(i), 1d-18)
          exc(ix,iy,iz) = exc(ix,iy,iz) + eexc_1d(i) / rho_safe
        end do
        end do
        end do
!$omp end parallel do
      end if

      if (present(eexc)) then
        eexc = eexc + reshape(eexc_1d, (/nx, ny, nz/))
      end if

      if (present(rdedd)) then
!$omp parallel do collapse(2) default(none) private(ix, iy, iz, idir, i) &
!$omp shared(nx, ny, nz, rdedd, vgrad_1d, grho_norm_1d, grho)
        do iz = 1, nz
        do iy = 1, ny
        do ix = 1, nx
          i = (iz-1) * nx * ny + (iy-1) * nx + ix
          if (grho_norm_1d(i) <= grad_floor) cycle
          do idir = 1, 3
            rdedd(ix,iy,iz,idir) = rdedd(ix,iy,iz,idir) - &
                 vgrad_1d(i) * grho(ix,iy,iz,idir) / grho_norm_1d(i)
          end do
        end do
        end do
        end do
!$omp end parallel do
      end if

      if (present(payload)) then
        payload%use_tau_operator = .true.
        payload%use_laplacian_operator = .false.
        payload%vtau_has_shadow_values = .false.
        payload%vsigma_has_shadow_values = .false.
        if (allocated(payload%vtau%f)) deallocate(payload%vtau%f)
        if (allocated(payload%vsigma%f)) deallocate(payload%vsigma%f)
        allocate(payload%vtau%f(nx, ny, nz))
        allocate(payload%vsigma%f(nx, ny, nz))
        payload%vtau%f = reshape(vtau_1d, (/nx, ny, nz/))
        vsigma_1d = 0d0
        where (grho_norm_1d > grad_floor)
          vsigma_1d = 0.5d0 * vgrad_1d / grho_norm_1d
        end where
        payload%vsigma%f = reshape(vsigma_1d, (/nx, ny, nz/))
      end if
    end subroutine exec_builtin_r2scan



#ifdef USE_LIBXC
    subroutine exec_libxc(ii)
#if XC_MAJOR_VERSION >= 5
      use, intrinsic :: iso_c_binding
#endif
      use nvtx
      implicit none
      integer, intent(in) :: ii
      integer :: ix,iy,iz,idir
      ! character(256) :: name

      real(8) :: rho_1d(nl), grho_1d(nl, 3)
      real(8) :: sigma_1d(nl), rlrho_1d(nl), tau_1d(nl)
      real(8) :: exc_tmp_1d(nl), vxc_tmp_1d(nl)
      real(8) :: gvxc_tmp_1d(nl), lvxc_tmp_1d(nl), tvxc_tmp_1d(nl)
      real(8),allocatable :: gvxc_tmp(:,:,:)
#if XC_MAJOR_VERSION <= 4
      integer :: np
#else
      integer(c_size_t) :: np
#endif

      call nvtxStartRange('exec_libxc', __LINE__)

      rho_1d = reshape(rho, (/nl/))
      if (present(rho_nlcc)) then
        rho_1d = rho_1d + reshape(rho_nlcc, (/nl/))
      end if

      if (xc%use_gradient) then
        sigma_1d = reshape( &
          & grho(:,:,:,1)**2+grho(:,:,:,2)**2+grho(:,:,:,3)**2, (/nl/) &
          & )
      end if

      if (xc%use_laplacian) then
        rlrho_1d = reshape(rlrho, (/nl/))
      end if

      if (xc%use_kinetic_energy) then
        tau_1d = reshape(tau, (/nl/))
      end if

      np = nl
#if XC_MAJOR_VERSION <= 4
      select case (xc_f90_info_family(xc%info(ii)))
#else
      select case (xc_f90_func_info_get_family(xc%info(ii)))
#endif

      case(XC_FAMILY_LDA)
        call xc_f90_lda_exc_vxc( &
          & xc%func(ii), np, rho_1d(1), &
          & exc_tmp_1d(1), vxc_tmp_1d(1) &
          & )

      case(XC_FAMILY_GGA)
        call xc_f90_gga_exc_vxc( &
          & xc%func(ii), np, rho_1d(1), sigma_1d(1), &
          & exc_tmp_1d(1), vxc_tmp_1d(1), gvxc_tmp_1d(1) &
          & )

      case(XC_FAMILY_MGGA)
        call xc_f90_mgga_exc_vxc( &
          & xc%func(ii), np, rho_1d(1), sigma_1d(1), rlrho_1d(1), tau_1d(1), &
          & exc_tmp_1d(1), vxc_tmp_1d(1), &
          & gvxc_tmp_1d(1), lvxc_tmp_1d(1), tvxc_tmp_1d(1))
      end select

      if (present(exc)) then
        exc = exc + reshape(exc_tmp_1d, (/nx, ny, nz/))
      endif

      if (present(eexc)) then
        eexc = eexc + reshape(exc_tmp_1d, (/nx, ny, nz/)) * rho
      endif

      if (present(vxc)) then
         vxc = vxc + reshape(vxc_tmp_1d, (/nx, ny, nz/))
      endif

      if (xc%use_gradient) then
        allocate(gvxc_tmp(nx, ny, nz))
        gvxc_tmp = reshape(gvxc_tmp_1d, (/nx, ny, nz/))
        do idir=1,3
          do iz=1,nz
          do iy=1,ny
          do ix=1,nx
            rdedd(ix,iy,iz,idir)=rdedd(ix,iy,iz,idir)-2.d0*gvxc_tmp(ix,iy,iz)*grho(ix,iy,iz,idir)
          enddo
          enddo
          enddo
        enddo
        deallocate(gvxc_tmp)
      endif
      
      call nvtxEndRange
      
      return
    end subroutine exec_libxc
#endif

#ifdef USE_LIBXC
    subroutine exec_libxc_spin(ii)
#if XC_MAJOR_VERSION >= 5
      use, intrinsic :: iso_c_binding
#endif
!      use salmon_global, only: rho_threshold,grho_threshold
      implicit none
      integer, intent(in) :: ii
      integer :: ix,iy,iz,idir,is
      ! character(256) :: name

      real(8) :: rho_1d(2) !, grho_1d(nl, 3)
      real(8) :: sigma_1d(3), rlrho_1d(2), tau_1d(2)
      real(8) :: exc_tmp_1d(1), vxc_tmp_1d(2)
      real(8) :: gvxc_tmp_1d(3) !, lvxc_tmp_1d(1), tvxc_tmp_1d(1)
!      real(8),parameter :: rho_threshold=1.d-16, grho_threshold=1.d-32
#if XC_MAJOR_VERSION <= 4
      integer :: np
#else
      integer(c_size_t) :: np
#endif

      if(pp%flag_nlcc)then
        stop "libxc for spin-polarised systems with nlcc is currently not implemented."
      end if

      do iz=1,nz
      do iy=1,ny
      do ix=1,nx
         exc_tmp_1d=0.d0
         vxc_tmp_1d=0.d0
         gvxc_tmp_1d=0.d0
         rho_1d(1) = rho_s(ix,iy,iz,1)
         rho_1d(2) = rho_s(ix,iy,iz,2)
         if (xc%use_gradient) then
           sigma_1d(1) = grho_s(ix,iy,iz,1,1)**2+grho_s(ix,iy,iz,1,2)**2+grho_s(ix,iy,iz,1,3)**2
           sigma_1d(2) = grho_s(ix,iy,iz,1,1)*grho_s(ix,iy,iz,2,1) &
                        +grho_s(ix,iy,iz,1,2)*grho_s(ix,iy,iz,2,2) &
                        +grho_s(ix,iy,iz,1,3)*grho_s(ix,iy,iz,2,3)
           sigma_1d(3) = grho_s(ix,iy,iz,2,1)**2+grho_s(ix,iy,iz,2,2)**2+grho_s(ix,iy,iz,2,3)**2
         endif
!  default 1.d-6 and 1.d-10 in xc_gga_drivers.f90 of QE
!         if(rho_1d(1)<=rho_threshold)    rho_1d(1)=0.d0
!         if(rho_1d(2)<=rho_threshold)    rho_1d(2)=0.d0
!         if(sigma_1d(1)<=grho_threshold) sigma_1d(1)=0.d0
!         if(sigma_1d(2)<=grho_threshold) sigma_1d(2)=0.d0
!         if(sigma_1d(3)<=grho_threshold) sigma_1d(3)=0.d0

!      rho_1d = reshape(rho, (/nl/))
!      if (present(rho_nlcc)) then
!        rho_1d = rho_1d + reshape(rho_nlcc, (/nl/))
!      end if
!      if (xc%use_gradient) then
!        sigma_1d = reshape( &
!          & grho(:,:,:,1)**2+grho(:,:,:,2)**2+grho(:,:,:,3)**2, (/nl/) &
!          & )
!      end if
!      if (xc%use_laplacian) then
!        rlrho_1d = reshape(rlrho, (/nl/))
!      end if
!      if (xc%use_kinetic_energy) then
!        tau_1d = reshape(tau, (/nl/))
!      end if
!
         np =  1 !nl
#if XC_MAJOR_VERSION <= 4
         select case (xc_f90_info_family(xc%info(ii)))
#else
         select case (xc_f90_func_info_get_family(xc%info(ii)))
#endif

         case(XC_FAMILY_LDA)
           call xc_f90_lda_exc_vxc( &
             & xc%func(ii), 1, rho_1d(1), &
             & exc_tmp_1d(1), vxc_tmp_1d(1) &
             & )

         case(XC_FAMILY_GGA)
           call xc_f90_gga_exc_vxc( &
             & xc%func(ii), 1, rho_1d(1), sigma_1d(1), &
             & exc_tmp_1d(1), vxc_tmp_1d(1), gvxc_tmp_1d(1) &
             & )

         case(XC_FAMILY_MGGA)
            stop "libxc for spin-polarised systems with mGGA is currently not implemented."
!           call xc_f90_mgga_exc_vxc( &
!             & xc%func(ii), np, rho_1d(1), sigma_1d(1), rlrho_1d(1), tau_1d(1), &
!             & exc_tmp_1d(1), vxc_tmp_1d(1), &
!             & gvxc_tmp_1d(1), lvxc_tmp_1d(1), tvxc_tmp_1d(1))
         end select

         if (present(eexc)) then
           eexc(ix,iy,iz)=eexc(ix,iy,iz)+exc_tmp_1d(1)*(rho_1d(1)+rho_1d(2))
         endif
         if (present(vxc_s)) then
           vxc_s(ix,iy,iz,1)=vxc_s(ix,iy,iz,1)+vxc_tmp_1d(1)
           vxc_s(ix,iy,iz,2)=vxc_s(ix,iy,iz,2)+vxc_tmp_1d(2)
         endif
!      if (present(exc)) then
!        exc = exc + reshape(exc_tmp_1d, (/nx, ny, nz/))
!      endif
!      if (present(eexc)) then
!        eexc = eexc + reshape(exc_tmp_1d, (/nx, ny, nz/)) * rho
!      endif
!      if (present(vxc)) then
!         vxc = vxc + reshape(vxc_tmp_1d, (/nx, ny, nz/))
!      endif
         if (xc%use_gradient) then
           do idir=1,3
             rdedd_s(ix,iy,iz,1,idir)=rdedd_s(ix,iy,iz,1,idir) & 
                 & -2.d0*gvxc_tmp_1d(1)*grho_s(ix,iy,iz,1,idir) & 
                 &      -gvxc_tmp_1d(2)*grho_s(ix,iy,iz,2,idir)  
             rdedd_s(ix,iy,iz,2,idir)=rdedd_s(ix,iy,iz,2,idir) &
                 & -2.d0*gvxc_tmp_1d(3)*grho_s(ix,iy,iz,2,idir) &
                 &      -gvxc_tmp_1d(2)*grho_s(ix,iy,iz,1,idir)  
           enddo
         endif
      enddo
      enddo
      enddo


      return
    end subroutine exec_libxc_spin
#endif

  end subroutine calc_xc


  
  function lower(s) result(lower_s)
     implicit none
     character(*), intent(in) :: s
     character(len(s)) :: lower_s
     integer :: i
     
     do i = 1, len(s)
       if (('A' <= s(i:i)) .and. (s(i:i) <= 'Z')) then
         lower_s(i:i) = char(ichar(s(i:i)) + 32)
       else
         lower_s(i:i) = s(i:i)
       end if 
     end do
     
     return
   end function lower

end module
