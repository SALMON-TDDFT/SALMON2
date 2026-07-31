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
module Conjugate_Gradient
  implicit none

  type halo_cg_info
    integer :: id_src,id_dst,ifrag_src,ifrag_dst,dvec(3),length(3),dsp_send(3),axis
  end type halo_cg_info

  type(halo_cg_info), allocatable, target, save :: flux_halo_cache(:)
  integer, save :: flux_halo_cache_n = -1
  integer, save :: flux_halo_cache_frag = -1
  integer, save :: flux_halo_cache_nfrag = -1
  integer, save :: flux_halo_cache_isize = -1
  integer, save :: flux_halo_cache_ncore(3) = -1
  integer, save :: flux_halo_cache_nbuf(3) = -1
  integer, save :: flux_halo_cache_radius(3) = -1
  integer, save :: flux_halo_cache_num_fragment(3) = -1

contains

subroutine gscg_rwf(ncg,mg,system,info,stencil,ppg,vlocal,srg,xk,hxk,gk,cg,dc,flux_weight,flux_iter)
  use structures
  use timer
  use hamiltonian, only: hpsi
  use communication, only: comm_summation, comm_is_root
  use parallelization, only: nproc_id_global
  use salmon_global, only: yn_preconditioning, yn_dc_lcfo_flux
  !$ use omp_lib
  implicit none
  integer           ,intent(in) :: ncg
  type(s_rgrid)     ,intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_stencil)   ,intent(in) :: stencil
  type(s_pp_grid)   ,intent(in) :: ppg
  type(s_scalar)    ,intent(in) :: vlocal(system%nspin)
  type(s_orbital)               :: xk,hxk,gk
  type(s_sendrecv_grid)         :: srg
  type(s_cg)                    :: cg
  type(s_dcdft),optional,intent(in) :: dc
  real(8),optional,intent(in) :: flux_weight
  integer,optional,intent(in) :: flux_iter
  type(s_orbital) :: flux_hxk,flux_hxk_raw,flux_hpk
  !
  integer,parameter :: nd=4
  integer :: nspin,io,ispin,io_s,io_e,is(3),ie(3),iy,iz
  integer :: ierr,icg,reject_local,reject_global
  real(8),parameter :: ep0=0.0d0
  real(8),parameter :: ep1=1.0d-15
  real(8),parameter :: c1  = 2.0d0
  real(8) :: rwork(9),W(2),c
  real(8),allocatable :: rb(:,:),E(:,:),E1(:,:),gkgk(:,:),bk(:,:),res(:,:)
  real(8),allocatable :: utmp3(:,:,:),wtmp2(:,:,:)
  real(8) :: utmp2(2,2),btmp2(2,2)
  logical :: use_fixed_flux,any_rejected,allow_cg_early_exit

  if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ gscg"

  call timer_begin(LOG_GSCG_ISOLATED_CALC)
  nspin = system%nspin
  is = mg%is
  ie = mg%ie
  io_s = info%io_s
  io_e = info%io_e
  use_fixed_flux = present(dc) .and. present(flux_weight) .and. yn_dc_lcfo_flux == 'y'
  if(use_fixed_flux) use_fixed_flux = flux_weight > 0d0
  allow_cg_early_exit = .not. use_fixed_flux

  if(.not. allocated(cg%pk%rwf)) then
    call allocate_orbital_real(nspin,mg,info,cg%pk)
    call allocate_orbital_real(nspin,mg,info,cg%pre_gk)
    call allocate_orbital_real(nspin,mg,info,cg%hpk)
    !$acc enter data copyin(cg)
  end if
  
  allocate(rb(system%nspin,system%no))
  allocate(E (system%nspin,system%no))
  allocate(E1(system%nspin,system%no))
  allocate(gkgk(system%nspin,system%no))
  allocate(bk(system%nspin,system%no))
  allocate(wtmp2(6,system%nspin,system%no))
  allocate(utmp3(2,system%nspin,system%no))
  allocate(res(system%nspin,system%no))
  res = 0.0d0
  if(use_fixed_flux) then
    call allocate_orbital_real(nspin,mg,info,flux_hxk_raw)
    call allocate_orbital_real(nspin,mg,info,flux_hxk)
    call allocate_orbital_real(nspin,mg,info,flux_hpk)
    flux_hxk_raw%rwf = 0d0
    flux_hxk%rwf = 0d0
    flux_hpk%rwf = 0d0
    call apply_dc_lcfo_flux_hpsi_rwf(mg,system,info,stencil,dc,xk,flux_hxk_raw)
    flux_hxk%rwf = flux_weight*flux_hxk_raw%rwf
  end if
  
  call timer_end(LOG_GSCG_ISOLATED_CALC)

  call timer_begin(LOG_GSCG_ISOLATED_HPSI)
  call hpsi(xk,hxk,info,mg,vlocal,system,stencil,srg,ppg)
  if(use_fixed_flux) then
    call print_fixed_flux_diagnostics
    call add_fixed_flux_operator(hxk,flux_hxk)
  end if
  call timer_end(LOG_GSCG_ISOLATED_HPSI)

  call timer_begin(LOG_GSCG_ISOLATED_CALC)

  E1=1.0d10

  call inner_product(mg,system,info,xk,hxk,E)

#ifdef USE_OPENACC
!$acc parallel loop private(io,ispin,iz,iy) collapse(4)
#else
!$omp parallel do private(io,ispin,iz,iy) collapse(4)
#endif
  do io=io_s,io_e
  do ispin=1,nspin
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
    gk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = -2.0d0*( hxk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
    & - E(ispin,io)* xk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) )
  end do
  end do
  end do
  end do
  call inner_product(mg,system,info,gk,gk,rb)

  do icg=1,Ncg+1

#ifdef USE_OPENACC
!$acc parallel loop private(io,ispin,iz,iy) collapse(4)
#else
!$omp parallel do private(io,ispin,iz,iy) collapse(4)
#endif
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      cg%pre_gk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = gk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) ! pre_gk==Pgk
    end do
    end do
    end do
    end do

    res = rb /c1**2

! --- Convergence check ---

    if ( allow_cg_early_exit .and. all(rb < ep0) ) exit
    if ( allow_cg_early_exit .and. all(abs(E-E1)<ep1) ) exit
    if ( icg==Ncg+1 ) exit

! --- Preconditioning ---

    if(yn_preconditioning=='y')then
      call preconditioning_rgk(mg,system,info,gk,cg%pre_gk)
    end if

! --- orthogonalization

! ---

    call inner_product(mg,system,info,cg%pre_gk,gk,rb) ! pre_gk==Pgk

    if ( icg==1 ) then
#ifdef USE_OPENACC
!$acc parallel loop private(io,ispin,iz,iy) collapse(4)
#else
!$omp parallel do private(io,ispin,iz,iy) collapse(4)
#endif
      do io=io_s,io_e
      do ispin=1,nspin
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
        cg%pk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = cg%pre_gk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) ! pre_gk==Pgk
      end do
      end do
      end do
      end do        
    else
      bk = rb/gkgk
#ifdef USE_OPENACC
!$acc parallel loop private(io,ispin,iz,iy) collapse(4)
#else
!$omp parallel do private(io,ispin,iz,iy) collapse(4)
#endif
      do io=io_s,io_e
      do ispin=1,nspin
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
        cg%pk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = cg%pre_gk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
        & + bk(ispin,io)*cg%pk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1)
      end do
      end do
      end do
      end do
    end if
    gkgk = rb
    call timer_end(LOG_GSCG_ISOLATED_CALC)

    call timer_begin(LOG_GSCG_ISOLATED_HPSI)
    call hpsi(cg%pk,cg%hpk,info,mg,vlocal,system,stencil,srg,ppg)
    if(use_fixed_flux) then
      flux_hpk%rwf = 0d0
      call apply_dc_lcfo_flux_hpsi_rwf(mg,system,info,stencil,dc,cg%pk,flux_hpk)
      flux_hpk%rwf = flux_weight*flux_hpk%rwf
      call add_fixed_flux_operator(cg%hpk,flux_hpk)
    end if
    call timer_end(LOG_GSCG_ISOLATED_HPSI)

    call timer_begin(LOG_GSCG_ISOLATED_CALC)
    call inner_product(mg,system,info,xk   ,xk    ,wtmp2(1,:,:))
    call inner_product(mg,system,info,cg%pk,xk    ,wtmp2(2,:,:))
    call inner_product(mg,system,info,cg%pk,cg%pk ,wtmp2(3,:,:))
    call inner_product(mg,system,info,xk   ,hxk   ,wtmp2(4,:,:))
    call inner_product(mg,system,info,cg%pk,hxk   ,wtmp2(5,:,:))
    call inner_product(mg,system,info,cg%pk,cg%hpk,wtmp2(6,:,:))

    do io=io_s,io_e
    do ispin=1,nspin
      btmp2(1,1)=wtmp2(1,ispin,io)
      btmp2(2,1)=wtmp2(2,ispin,io)
      btmp2(1,2)=wtmp2(2,ispin,io)
      btmp2(2,2)=wtmp2(3,ispin,io)
      utmp2(1,1)=wtmp2(4,ispin,io)
      utmp2(2,1)=wtmp2(5,ispin,io)
      utmp2(1,2)=wtmp2(5,ispin,io)
      utmp2(2,2)=wtmp2(6,ispin,io)
      call dsygv(1,'V','U',2,utmp2,2,btmp2,2,W,rwork,9,ierr)
      if ( abs(W(1)-E(ispin,io))>1.d-1 .and. abs(W(2)-E(ispin,io))<=1.d-1 ) then
        utmp2(1,1)=utmp2(1,2)
        utmp2(2,1)=utmp2(2,2)
        W(1)=W(2)
      end if

      !- Fix the phase -
      c=utmp2(1,1)
      if( c<0.0d0 ) then
        utmp2(1,1)=-utmp2(1,1)
        utmp2(2,1)=-utmp2(2,1)
      end if
      utmp3(1:2,ispin,io) = utmp2(1:2,1)
      E1(ispin,io) = E(ispin,io)
      E(ispin,io)  = W(1)
    end do ! ispin
    end do ! io
    
#ifdef USE_OPENACC
!$acc parallel loop private(io,ispin,iz,iy) collapse(4)
#else
!$omp parallel do private(io,ispin,iz,iy) collapse(4)
#endif
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      hxk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = &
        &   utmp3(1,ispin,io)* hxk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
        & + utmp3(2,ispin,io)* cg%hpk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1)
      gk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = -2.0d0*( &
        & hxk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
        & - E(ispin,io)*( utmp3(1,ispin,io) * xk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
        &               + utmp3(2,ispin,io) * cg%pk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) ) )
    end do
    end do
    end do
    end do

    call inner_product(mg,system,info,gk,gk,rb)
    bk = -1d0
    do io=io_s,io_e
    do ispin=1,nspin
      if ( rb(ispin,io)/res(ispin,io)>1.0d8 ) then
        E(ispin,io) = E1(ispin,io)
        bk(ispin,io) = 1d0
      end if
    end do
    end do
    any_rejected = any(bk(:,io_s:io_e) >= 0d0)
    reject_local = 0
    if(any_rejected) reject_local = 1
    if(use_fixed_flux) then
      call comm_summation(reject_local,reject_global,dc%icomm_tot)
    else
      call comm_summation(reject_local,reject_global,info%icomm_rko)
    end if
    any_rejected = reject_global > 0

#ifdef USE_OPENACC
!$acc parallel loop private(io,ispin,iz,iy) collapse(4)
#else
!$omp parallel do private(io,ispin,iz,iy) collapse(4)
#endif
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      if ( bk(ispin,io) < 0d0 ) then
        xk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = &
        &   utmp3(1,ispin,io) * xk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
        & + utmp3(2,ispin,io) * cg%pk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1)
      end if
    end do
    end do
    end do
    end do

    if(any_rejected) then
      call timer_end(LOG_GSCG_ISOLATED_CALC)
      call timer_begin(LOG_GSCG_ISOLATED_HPSI)
      call hpsi(xk,hxk,info,mg,vlocal,system,stencil,srg,ppg)
      if(use_fixed_flux) then
        flux_hxk_raw%rwf = 0d0
        call apply_dc_lcfo_flux_hpsi_rwf(mg,system,info,stencil,dc,xk,flux_hxk_raw)
        flux_hxk%rwf = flux_weight*flux_hxk_raw%rwf
        call add_fixed_flux_operator(hxk,flux_hxk)
      end if
      call timer_end(LOG_GSCG_ISOLATED_HPSI)
      call timer_begin(LOG_GSCG_ISOLATED_CALC)
#ifdef USE_OPENACC
!$acc parallel loop private(io,ispin,iz,iy) collapse(4)
#else
!$omp parallel do private(io,ispin,iz,iy) collapse(4)
#endif
      do io=io_s,io_e
      do ispin=1,nspin
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
        gk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = -2.0d0*( &
        & hxk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
        & - E(ispin,io)*xk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) )
      end do
      end do
      end do
      end do
    end if
    call inner_product(mg,system,info,gk,gk,rb)

  end do ! icg
  
  deallocate( utmp3,wtmp2 )
  deallocate( bk,gkgk,E1,E,rb,res )
  if(allocated(flux_hxk_raw%rwf)) deallocate(flux_hxk_raw%rwf)
  if(allocated(flux_hxk%rwf)) deallocate(flux_hxk%rwf)
  if(allocated(flux_hpk%rwf)) deallocate(flux_hpk%rwf)

  call timer_end(LOG_GSCG_ISOLATED_CALC)

  return
contains

subroutine print_fixed_flux_diagnostics
  implicit none
  integer :: io,ispin,iter_value
  real(8) :: e_core,e_flux,e_flux_raw,n_core,n_flux,n_flux_raw
  real(8) :: e_core_local,e_flux_local,e_flux_raw_local
  real(8) :: n_core_local,n_flux_local,n_flux_raw_local
  real(8) :: e_core_frag,e_flux_frag,e_flux_raw_frag
  real(8) :: n_core_frag,n_flux_frag,n_flux_raw_frag
  real(8) :: occ,ratio_e,ratio_n,ratio_e_raw,ratio_n_raw,weight_value
  logical :: do_print

  do_print = .true.
  iter_value = -1
  weight_value = 0d0
  if(present(flux_weight)) weight_value = flux_weight
  if(present(flux_iter)) then
    iter_value = flux_iter
    if(present(flux_weight)) then
      if(weight_value < 1d0 - 1d-12) then
        do_print = mod(iter_value,100) == 0
      else
        do_print = mod(iter_value,1000) == 0
      end if
    else
      do_print = mod(iter_value,1000) == 0
    end if
  end if
  if(.not. do_print) return

  call inner_product(mg,system,info,xk,hxk,wtmp2(1,:,:))
  call inner_product(mg,system,info,xk,flux_hxk,wtmp2(2,:,:))
  call inner_product(mg,system,info,hxk,hxk,wtmp2(3,:,:))
  call inner_product(mg,system,info,flux_hxk,flux_hxk,wtmp2(4,:,:))
  call inner_product(mg,system,info,xk,flux_hxk_raw,wtmp2(5,:,:))
  call inner_product(mg,system,info,flux_hxk_raw,flux_hxk_raw,wtmp2(6,:,:))

  e_core_local = 0d0
  e_flux_local = 0d0
  e_flux_raw_local = 0d0
  n_core_local = 0d0
  n_flux_local = 0d0
  n_flux_raw_local = 0d0
  do io=io_s,io_e
  do ispin=1,nspin
    occ = system%rocc(io,1,ispin)*system%wtk(1)
    e_core_local = e_core_local + occ*wtmp2(1,ispin,io)
    e_flux_local = e_flux_local + occ*wtmp2(2,ispin,io)
    e_flux_raw_local = e_flux_raw_local + occ*wtmp2(5,ispin,io)
    n_core_local = n_core_local + occ*wtmp2(3,ispin,io)
    n_flux_local = n_flux_local + occ*wtmp2(4,ispin,io)
    n_flux_raw_local = n_flux_raw_local + occ*wtmp2(6,ispin,io)
  end do
  end do
  call comm_summation(e_core_local,e_core,info%icomm_o)
  call comm_summation(e_flux_local,e_flux,info%icomm_o)
  call comm_summation(e_flux_raw_local,e_flux_raw,info%icomm_o)
  call comm_summation(n_core_local,n_core,info%icomm_o)
  call comm_summation(n_flux_local,n_flux,info%icomm_o)
  call comm_summation(n_flux_raw_local,n_flux_raw,info%icomm_o)
  e_core_frag = e_core
  e_flux_frag = e_flux
  e_flux_raw_frag = e_flux_raw
  n_core_frag = n_core
  n_flux_frag = n_flux
  n_flux_raw_frag = n_flux_raw
  if(present(dc)) then
    if(dc%id_frag /= 0) then
      e_core_frag = 0d0
      e_flux_frag = 0d0
      e_flux_raw_frag = 0d0
      n_core_frag = 0d0
      n_flux_frag = 0d0
      n_flux_raw_frag = 0d0
    end if
    call comm_summation(e_core_frag,e_core,dc%icomm_tot)
    call comm_summation(e_flux_frag,e_flux,dc%icomm_tot)
    call comm_summation(e_flux_raw_frag,e_flux_raw,dc%icomm_tot)
    call comm_summation(n_core_frag,n_core,dc%icomm_tot)
    call comm_summation(n_flux_frag,n_flux,dc%icomm_tot)
    call comm_summation(n_flux_raw_frag,n_flux_raw,dc%icomm_tot)
  end if

  n_core = sqrt(max(0d0,n_core))
  n_flux = sqrt(max(0d0,n_flux))
  n_flux_raw = sqrt(max(0d0,n_flux_raw))
  ratio_e = 0d0
  ratio_n = 0d0
  ratio_e_raw = 0d0
  ratio_n_raw = 0d0
  if(abs(e_core) > 1d-300) ratio_e = e_flux/e_core
  if(n_core > 1d-300) ratio_n = n_flux/n_core
  if(abs(e_core) > 1d-300) ratio_e_raw = e_flux_raw/e_core
  if(n_core > 1d-300) ratio_n_raw = n_flux_raw/n_core
  if(comm_is_root(nproc_id_global)) then
    write(*,'(1x,a,i8,a,es12.5,a,es14.6,a,es14.6,a,es12.5,a,es14.6,a,es12.5)') &
    & '[DC-LCFO-FLUX-DIAG] iter=',iter_value, &
    & ' weight=',weight_value, &
    & ' Ecore=',e_core, &
    & ' Eflux_raw=',e_flux_raw, &
    & ' raw/Ecore=',ratio_e_raw, &
    & ' Eflux=',e_flux, &
    & ' weighted/Ecore=',ratio_e
    write(*,'(1x,a,i8,a,es14.6,a,es14.6,a,es12.5,a,es14.6,a,es12.5)') &
    & '[DC-LCFO-FLUX-DIAG-NORM] iter=',iter_value, &
    & ' Hcore_norm=',n_core, &
    & ' Hflux_raw_norm=',n_flux_raw, &
    & ' raw/Hcore=',ratio_n_raw, &
    & ' Hflux_norm=',n_flux, &
    & ' weighted/Hcore=',ratio_n
  end if
end subroutine print_fixed_flux_diagnostics

subroutine add_fixed_flux_operator(hpsi_target,flux_source)
  implicit none
  type(s_orbital),intent(inout) :: hpsi_target
  type(s_orbital),intent(in) :: flux_source

#ifdef USE_OPENACC
!$acc parallel loop private(io,ispin,iz,iy) collapse(4)
#else
!$omp parallel do private(io,ispin,iz,iy) collapse(4)
#endif
  do io=io_s,io_e
  do ispin=1,nspin
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
    hpsi_target%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = &
    & hpsi_target%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
    & + flux_source%rwf(is(1):ie(1),iy,iz,ispin,io,1,1)
  end do
  end do
  end do
  end do
end subroutine add_fixed_flux_operator

subroutine inner_product(mg,system,info,psi1,psi2,rbox)
  !$ use omp_lib
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_orbital),intent(in) :: psi1,psi2
  real(8),intent(out) :: rbox(system%nspin,system%no)
  !
  integer :: io,ispin,nspin
  integer :: ix,iy,iz
  real(8) :: rbox2(system%nspin,system%no)
  real(8) :: sum0
  nspin = system%nspin

  rbox2 = 0.d0
#ifdef USE_OPENACC
!$acc parallel loop collapse(2) private(io,ispin,sum0,iz,iy,ix)
#else
!$OMP parallel do collapse(2) private(io,ispin,sum0,iz,iy,ix)
#endif
  do io=info%io_s,info%io_e
  do ispin=1,nspin
    sum0 = 0d0
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      sum0 = sum0 + psi1%rwf(ix,iy,iz,ispin,io,1,1) * psi2%rwf(ix,iy,iz,ispin,io,1,1)
    end do
    end do
    end do
    rbox2(ispin,io) = sum0 * system%hvol
  end do
  end do
  call timer_end(LOG_GSCG_ISOLATED_CALC)

  call timer_begin(LOG_GSCG_ISOLATED_COMM_COLL)
  call comm_summation(rbox2,rbox,nspin*system%no,info%icomm_r)
  call timer_end(LOG_GSCG_ISOLATED_COMM_COLL)

  call timer_begin(LOG_GSCG_ISOLATED_CALC)
end subroutine inner_product

subroutine preconditioning_rgk(mg,system,info,gk,pre_gk)
  !$ use omp_lib
  use preconditioning_sub, only: dstencil_preconditioning
  use structures
  use sendrecv_grid, only: update_overlap_real8
  use salmon_global, only: alpha_pre
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_orbital),intent(inout) :: gk
  type(s_orbital),intent(inout) :: pre_gk
  !
  integer :: io,ik,ispin,nspin
  integer :: ix,iy,iz
  real(8) :: alpha
  integer :: is(3),ie(3)
  logical :: is_enable_overlapping
  nspin = system%nspin

  alpha = alpha_pre

  if(info%if_divide_rspace) then
    call update_overlap_real8(srg, mg, gk%rwf)
  end if

  do ik=info%ik_s,info%ik_e
  do io=info%io_s,info%io_e
  do ispin=1,nspin
    call dstencil_preconditioning(mg%is_array,mg%ie_array,mg%is,  &
                                  mg%ie,mg%idx,mg%idy,mg%idz,system%hgs, &
                                  gk%rwf(:,:,:,ispin,io,ik,1), &
                                  pre_gk%rwf(:,:,:,ispin,io,ik,1),alpha)
  end do
  end do
  end do

end subroutine preconditioning_rgk

end subroutine gscg_rwf

!===================================================================================================================================

subroutine gscg_zwf(ncg,mg,system,info,stencil,ppg,vlocal,srg,xk,hxk,gk,cg,dc,flux_weight,flux_iter)
  use structures
  use timer
  use hamiltonian, only: hpsi
  use communication, only: comm_summation
  use salmon_global, only: yn_spinorbit,yn_preconditioning,yn_dc_lcfo_flux
  !$ use omp_lib
  implicit none
  integer           ,intent(in) :: ncg
  type(s_rgrid)     ,intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_stencil)   ,intent(in) :: stencil
  type(s_pp_grid)   ,intent(in) :: ppg
  type(s_scalar)    ,intent(in) :: vlocal(system%nspin)
  type(s_orbital)               :: xk,hxk,gk
  type(s_sendrecv_grid)         :: srg
  type(s_cg)                    :: cg
  type(s_dcdft),optional,intent(in) :: dc
  real(8),optional,intent(in) :: flux_weight
  integer,optional,intent(in) :: flux_iter
  !
  integer,parameter :: nd=4
  integer :: nspin,ik,io,ispin,ik_s,ik_e,io_s,io_e,is(3),ie(3),iy,iz
  integer :: ierr,icg
  real(8),parameter :: ep0=0.0d0
  real(8),parameter :: ep1=1.0d-15
  real(8),parameter :: c1  = 2.0d0
  real(8) :: rwork(9),W(2),r,c,d
  real(8),allocatable :: rb(:,:,:),E(:,:,:),E1(:,:,:),gkgk(:,:,:),bk(:,:,:),res(:,:,:)
  complex(8),allocatable :: utmp3(:,:,:,:),wtmp2(:,:,:,:),zb(:,:,:)
  complex(8) :: utmp2(2,2),btmp2(2,2)
  complex(8) :: work(9),zphase,ztmp

  if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ gscg"
  if(present(dc) .and. present(flux_weight) .and. yn_dc_lcfo_flux == 'y') then
    if(flux_weight > 0d0) stop "DC-LCFO-Flux CG is implemented for real orbitals only"
  end if

  call timer_begin(LOG_GSCG_PERIODIC_CALC)
  nspin = system%nspin
  is = mg%is
  ie = mg%ie
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e

  if(.not. allocated(cg%pk%zwf)) then
    call allocate_orbital_complex(nspin,mg,info,cg%pk)
    call allocate_orbital_complex(nspin,mg,info,cg%pre_gk)
    call allocate_orbital_complex(nspin,mg,info,cg%hpk)
    !$acc enter data copyin(cg)
  end if
  
  allocate(zb     (system%nspin,system%no,system%nk))
  allocate(rb     (system%nspin,system%no,system%nk))
  allocate(E      (system%nspin,system%no,system%nk))
  allocate(E1     (system%nspin,system%no,system%nk))
  allocate(res    (system%nspin,system%no,system%nk))
  allocate(gkgk   (system%nspin,system%no,system%nk))
  allocate(bk     (system%nspin,system%no,system%nk))
  allocate(wtmp2(6,system%nspin,system%no,system%nk))
  allocate(utmp3(2,system%nspin,system%no,system%nk))
  res = 0.0d0

  E1 = 1.0d10

  call timer_end(LOG_GSCG_PERIODIC_CALC)

  call timer_begin(LOG_GSCG_PERIODIC_HPSI)
  call hpsi(xk,hxk,info,mg,vlocal,system,stencil,srg,ppg)
  call timer_end(LOG_GSCG_PERIODIC_HPSI)

  call timer_begin(LOG_GSCG_PERIODIC_CALC)
  call inner_product(mg,system,info,xk,hxk,zb)
  E = dble(zb)
 
#ifdef USE_OPENACC
!$acc parallel loop private(ik,io,ispin,iz,iy) collapse(5)
#else
!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
#endif
  do ik=ik_s,ik_e
  do io=io_s,io_e
  do ispin=1,nspin
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
    gk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = -hxk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) &
    & + E(ispin,io,ik) * xk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
  end do
  end do
  end do
  end do
  end do

  call inner_product(mg,system,info,gk,gk,zb)
  rb = dble(zb)

  do icg=1,Ncg+1

#ifdef USE_OPENACC
!$acc parallel loop private(ik,io,ispin,iz,iy) collapse(5)
#else
!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
#endif
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      cg%pre_gk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = gk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
    end do
    end do
    end do
    end do
    end do
    
    res = rb

! --- Convergence check ---

    if ( all(rb < ep0) ) exit
    if ( all(abs(E-E1)<ep1) ) exit
    if ( icg==Ncg+1 ) exit

! --- Preconditioning ---

    if(yn_preconditioning=='y')then
      call preconditioning_zgk(mg,system,info,gk,cg%pre_gk)
    end if

! --- orthogonalization

! ---

    call inner_product(mg,system,info,cg%pre_gk,gk,zb)
    rb = dble(zb)

    if ( icg==1 ) then

#ifdef USE_OPENACC
!$acc parallel loop private(ik,io,ispin,iz,iy) collapse(5)
#else
!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
#endif
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,nspin
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
        cg%pk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = cg%pre_gk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
      end do
      end do
      end do
      end do
      end do

    else

      bk = rb/gkgk

#ifdef USE_OPENACC
!$acc parallel loop private(ik,io,ispin,iz,iy) collapse(5)
#else
!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
#endif
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,nspin
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
        cg%pk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = &
        & cg%pre_gk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) + &
        & bk(ispin,io,ik) * cg%pk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
      end do
      end do
      end do
      end do
      end do
      
    end if

    gkgk = rb
    call timer_end(LOG_GSCG_PERIODIC_CALC)

    call timer_begin(LOG_GSCG_PERIODIC_HPSI)
    call hpsi(cg%pk,cg%hpk,info,mg,vlocal,system,stencil,srg,ppg)
    call timer_end(LOG_GSCG_PERIODIC_HPSI)

    call timer_begin(LOG_GSCG_PERIODIC_CALC)

    call inner_product(mg,system,info,xk,   xk,    wtmp2(1,:,:,:))
    call inner_product(mg,system,info,cg%pk,xk,    wtmp2(2,:,:,:))
    call inner_product(mg,system,info,cg%pk,cg%pk, wtmp2(3,:,:,:))
    call inner_product(mg,system,info,xk,   hxk,   wtmp2(4,:,:,:))
    call inner_product(mg,system,info,cg%pk,hxk,   wtmp2(5,:,:,:))
    call inner_product(mg,system,info,cg%pk,cg%hpk,wtmp2(6,:,:,:))

    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
      btmp2(1,1)=wtmp2(1,ispin,io,ik)
      btmp2(2,1)=wtmp2(2,ispin,io,ik)
      btmp2(1,2)=wtmp2(2,ispin,io,ik)
      btmp2(2,2)=wtmp2(3,ispin,io,ik)
      utmp2(1,1)=wtmp2(4,ispin,io,ik)
      utmp2(2,1)=wtmp2(5,ispin,io,ik)
      utmp2(1,2)=wtmp2(5,ispin,io,ik)
      utmp2(2,2)=wtmp2(6,ispin,io,ik)
      ztmp=btmp2(1,2)
      ztmp=conjg(ztmp)
      btmp2(1,2)=ztmp
      ztmp=utmp2(1,2)
      ztmp=conjg(ztmp)
      utmp2(1,2)=ztmp
      call zhegv(1,'V','U',2,utmp2,2,btmp2,2,W,work,9,rwork,ierr)
      if ( abs(W(1)-E(ispin,io,ik))>1.d-1 .and. abs(W(2)-E(ispin,io,ik))<=1.d-1 ) then
        utmp2(1,1)=utmp2(1,2)
        utmp2(2,1)=utmp2(2,2)
        W(1)=W(2)
      end if
      !- Fix the phase -
      ztmp=utmp2(1,1)
      r=abs(ztmp)
      if (r > tiny(1.0d0)) then
        c=real(ztmp)/r
        d=aimag(ztmp)/r
        zphase=dcmplx(c,-d)
      else
        zphase=dcmplx(1.0d0,0.0d0)
      end if
      utmp2(1,1)=utmp2(1,1)*zphase
      utmp2(2,1)=utmp2(2,1)*zphase

      utmp3(1:2,ispin,io,ik) = utmp2(1:2,1)

      E1(ispin,io,ik) = E(ispin,io,ik)
      E (ispin,io,ik) = W(1)
    end do ! ispin
    end do ! io
    end do ! ik
    
#ifdef USE_OPENACC
!$acc parallel loop private(ik,io,ispin,iz,iy) collapse(5)
#else
!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
#endif
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      hxk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = &
      &   utmp3(1,ispin,io,ik) * hxk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) &
      & + utmp3(2,ispin,io,ik) * cg%hpk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
      gk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = - hxk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) &
      & + E(ispin,io,ik)*( utmp3(1,ispin,io,ik) * xk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) &
      &                  + utmp3(2,ispin,io,ik) * cg%pk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) )
    end do
    end do
    end do
    end do
    end do
    
    call inner_product(mg,system,info,gk,gk,zb)
    rb = dble(zb)
    
    bk = -1d0
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
      if ( rb(ispin,io,ik)/res(ispin,io,ik)>1.0d8 ) then
        E(ispin,io,ik) = E1(ispin,io,ik)
        bk(ispin,io,ik) = 1d0
      end if
    end do
    end do
    end do

#ifdef USE_OPENACC
!$acc parallel loop private(ik,io,ispin,iz,iy) collapse(5)
#else
!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
#endif
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      if ( bk(ispin,io,ik) < 0d0 ) then
        xk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = &
        & utmp3(1,ispin,io,ik) * xk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) &
        & + utmp3(2,ispin,io,ik) * cg%pk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
      end if
    end do
    end do
    end do
    end do
    end do

  end do ! icg

  deallocate( utmp3,wtmp2 )
  deallocate( bk,gkgk,E1,E,rb,res )
  deallocate( zb )

  call timer_end(LOG_GSCG_PERIODIC_CALC)

  return
contains

subroutine inner_product(mg,system,info,psi1,psi2,zbox)
  !$ use omp_lib
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_orbital),intent(in) :: psi1,psi2
  complex(8),intent(out) :: zbox(system%nspin,system%no,system%nk)
  !
  integer :: io,ik,ispin,nspin
  integer :: ix,iy,iz
  complex(8) :: zbox2(system%nspin,system%no,system%nk)
  complex(8) :: sum0
  nspin = system%nspin

  zbox2(:,:,:) = 0.d0
  
  if(yn_spinorbit=='n') then
#ifdef USE_OPENACC
!$acc parallel loop collapse(2) private(ik,io,ispin,sum0,iz,iy,ix)
#else
!$OMP parallel do collapse(2) private(ik,io,ispin,sum0,iz,iy,ix)
#endif
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do ispin=1,nspin
      sum0 = 0d0
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        sum0 = sum0 + conjg(psi1%zwf(ix,iy,iz,ispin,io,ik,1))*psi2%zwf(ix,iy,iz,ispin,io,ik,1)
      end do
      end do
      end do
      zbox2(ispin,io,ik) = sum0 * system%hvol
    end do
    end do
    end do
  else
#ifdef USE_OPENACC
!$acc parallel loop collapse(2) private(ik,io,ispin,sum0,iz,iy,ix)
#else
!$OMP parallel do collapse(2) private(ik,io,ispin,sum0,iz,iy,ix)
#endif
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
      sum0 = 0d0
      do ispin=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        sum0 = sum0 + conjg(psi1%zwf(ix,iy,iz,ispin,io,ik,1))*psi2%zwf(ix,iy,iz,ispin,io,ik,1)
      end do
      end do
      end do
      end do
      zbox2(1,io,ik) = sum0 * system%hvol
      zbox2(2,io,ik) = zbox2(1,io,ik)
    end do
    end do
  end if
  call timer_end(LOG_GSCG_PERIODIC_CALC)

  call timer_begin(LOG_GSCG_PERIODIC_COMM_COLL)
  call comm_summation(zbox2,zbox,nspin*system%no*system%nk,info%icomm_r)
  call timer_end(LOG_GSCG_PERIODIC_COMM_COLL)

  call timer_begin(LOG_GSCG_PERIODIC_CALC)
end subroutine inner_product

subroutine preconditioning_zgk(mg,system,info,gk,pre_gk)
  !$ use omp_lib
  use preconditioning_sub, only: zstencil_preconditioning,zstencil_nonorthogonal_preconditioning
  use structures
  use sendrecv_grid, only: update_overlap_complex8
  use salmon_global, only: yn_want_communication_overlapping,alpha_pre
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_orbital),intent(inout) :: gk
  type(s_orbital),intent(inout) :: pre_gk
  !
  integer :: io,ik,ispin,nspin
  real(8) :: alpha
  integer :: is(3),ie(3)
  logical :: is_enable_overlapping
  nspin = system%nspin

  alpha = alpha_pre

  is_enable_overlapping = (yn_want_communication_overlapping == 'y') .and. &
                          stencil%if_orthogonal .and. &
                          info%if_divide_rspace

  if(info%if_divide_rspace .and. .not. is_enable_overlapping) then
    call update_overlap_complex8(srg, mg, gk%zwf)
  end if

  if(stencil%if_orthogonal) then
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do ispin=1,nspin
      call zstencil_preconditioning(mg%is_array,mg%ie_array,mg%is,  &
                                    mg%ie,mg%idx,mg%idy,mg%idz,system%hgs, &
                                    gk%zwf(:,:,:,ispin,io,ik,1), &
                                    pre_gk%zwf(:,:,:,ispin,io,ik,1),alpha)
    end do
    end do
    end do
  else
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do ispin=1,nspin
      call zstencil_nonorthogonal_preconditioning(mg%is_array,mg%ie_array,mg%is,  &
                                    mg%ie,mg%idx,mg%idy,mg%idz, &
                                    gk%zwf(:,:,:,ispin,io,ik,1), &
                                    pre_gk%zwf(:,:,:,ispin,io,ik,1), &
                                    stencil%coef_lap0_nd1, &
                                    stencil%coef_lap_nd1,stencil%coef_nab_nd1, &
                                    stencil%coef_F,alpha)
    end do
    end do
    end do
  end if

end subroutine preconditioning_zgk

end subroutine gscg_zwf

subroutine prepare_dc_lcfo_flux_halo_cache(dc,ncore,stencil_radius,npg_base)
  use structures
  use salmon_global, only: num_fragment
  implicit none
  type(s_dcdft),intent(in) :: dc
  integer,intent(in) :: ncore(3),stencil_radius(3),npg_base
  integer :: nh(3),ir1(3),ir2(3),d(3)
  integer :: lx,ly,lz,ifrag,n,idir,nonzero_dirs,rank_same_subgroup
  logical :: cache_valid

  cache_valid = allocated(flux_halo_cache)
  cache_valid = cache_valid .and. flux_halo_cache_frag == dc%i_frag
  cache_valid = cache_valid .and. flux_halo_cache_nfrag == dc%n_frag
  cache_valid = cache_valid .and. flux_halo_cache_isize == dc%isize_tot
  cache_valid = cache_valid .and. all(flux_halo_cache_ncore == ncore)
  cache_valid = cache_valid .and. all(flux_halo_cache_nbuf == dc%nxyz_buffer)
  cache_valid = cache_valid .and. all(flux_halo_cache_radius == stencil_radius)
  cache_valid = cache_valid .and. all(flux_halo_cache_num_fragment == num_fragment)
  if(cache_valid) return

  nh = 0
  do n=1,3
    if(dc%nxyz_buffer(n) > ncore(n)) stop "DC-LCFO-Flux CG: buffer > domain"
    if(num_fragment(n) > 1 .and. dc%nxyz_buffer(n) < stencil_radius(n)) &
    & stop "DC-LCFO-Flux CG: buffer is smaller than the active stencil radius"
    if(num_fragment(n) > 1) nh(n) = 1
  end do

  if(allocated(flux_halo_cache)) deallocate(flux_halo_cache)
  allocate(flux_halo_cache(26))
  flux_halo_cache_n = 0
  do lx=-nh(1),nh(1)
  do ly=-nh(2),nh(2)
  do lz=-nh(3),nh(3)
    if(lx==0 .and. ly==0 .and. lz==0) cycle
    nonzero_dirs = count([lx,ly,lz] /= 0)
    if(nonzero_dirs /= 1) cycle

    flux_halo_cache_n = flux_halo_cache_n + 1
    flux_halo_cache(flux_halo_cache_n)%dvec = [lx,ly,lz]
    flux_halo_cache(flux_halo_cache_n)%axis = 0
    do idir=1,3
      if(flux_halo_cache(flux_halo_cache_n)%dvec(idir) /= 0) &
      & flux_halo_cache(flux_halo_cache_n)%axis = idir
    end do
    flux_halo_cache(flux_halo_cache_n)%id_dst = -1
    flux_halo_cache(flux_halo_cache_n)%id_src = -1
    flux_halo_cache(flux_halo_cache_n)%ifrag_dst = -1
    flux_halo_cache(flux_halo_cache_n)%ifrag_src = -1
    do ifrag=1,dc%n_frag
      ir1 = dc%ixyz_frag(:,ifrag)
      ir2 = dc%ixyz_frag(:,dc%i_frag) + flux_halo_cache(flux_halo_cache_n)%dvec*ncore
      d = mod(ir1 - ir2, dc%lg_tot%num(1:3))
      if(all(d == 0) .and. flux_halo_cache(flux_halo_cache_n)%id_dst < 0) then
        rank_same_subgroup = -1
        if(dc%id_frag < npg_base) rank_same_subgroup = (ifrag-1)*npg_base + dc%id_frag
        flux_halo_cache(flux_halo_cache_n)%id_dst = rank_same_subgroup
        flux_halo_cache(flux_halo_cache_n)%ifrag_dst = ifrag
      end if

      ir2 = dc%ixyz_frag(:,dc%i_frag) - flux_halo_cache(flux_halo_cache_n)%dvec*ncore
      d = mod(ir1 - ir2, dc%lg_tot%num(1:3))
      if(all(d == 0) .and. flux_halo_cache(flux_halo_cache_n)%id_src < 0) then
        rank_same_subgroup = -1
        if(dc%id_frag < npg_base) rank_same_subgroup = (ifrag-1)*npg_base + dc%id_frag
        flux_halo_cache(flux_halo_cache_n)%id_src = rank_same_subgroup
        flux_halo_cache(flux_halo_cache_n)%ifrag_src = ifrag
      end if
    end do
    if(flux_halo_cache(flux_halo_cache_n)%id_dst < 0 .or. &
    &  flux_halo_cache(flux_halo_cache_n)%id_src < 0) &
    & stop "DC-LCFO-Flux CG: missing halo rank"

    do n=1,3
      select case(flux_halo_cache(flux_halo_cache_n)%dvec(n))
      case(0)
        flux_halo_cache(flux_halo_cache_n)%length(n) = ncore(n)
        flux_halo_cache(flux_halo_cache_n)%dsp_send(n) = 0
      case(1)
        flux_halo_cache(flux_halo_cache_n)%length(n) = dc%nxyz_buffer(n)
        flux_halo_cache(flux_halo_cache_n)%dsp_send(n) = ncore(n) - dc%nxyz_buffer(n)
      case(-1)
        flux_halo_cache(flux_halo_cache_n)%length(n) = dc%nxyz_buffer(n)
        flux_halo_cache(flux_halo_cache_n)%dsp_send(n) = 0
      end select
    end do
  end do
  end do
  end do

  flux_halo_cache_frag = dc%i_frag
  flux_halo_cache_nfrag = dc%n_frag
  flux_halo_cache_isize = dc%isize_tot
  flux_halo_cache_ncore = ncore
  flux_halo_cache_nbuf = dc%nxyz_buffer
  flux_halo_cache_radius = stencil_radius
  flux_halo_cache_num_fragment = num_fragment
end subroutine prepare_dc_lcfo_flux_halo_cache

subroutine apply_dc_lcfo_flux_hpsi_rwf(mg,system,info,stencil,dc,psi,hpsi)
  use structures
  use communication, only: comm_isend, comm_irecv, comm_wait_all, comm_proc_null
  use dc_fragment_geometry, only: get_fragment_domain
  implicit none
  type(s_rgrid),          intent(in)    :: mg
  type(s_dft_system),     intent(in)    :: system
  type(s_parallel_info),  intent(in)    :: info
  type(s_stencil),        intent(in)    :: stencil
  type(s_dcdft),          intent(in)    :: dc
  type(s_orbital),        intent(in)    :: psi
  type(s_orbital),        intent(inout) :: hpsi
  !
  type(halo_cg_info), pointer :: halo(:)
  integer :: ncore(3),l(3),lo(3),hi(3)
  integer :: i_halo,n_halo,n,axis,dist
  integer :: ix,iy,iz,ixg,iyg,izg,il(3),ih(3)
  integer :: io,ilo,ispin,nio,itag_send,itag_recv,itag_dir
  integer :: npg_base,rank_send,rank_recv
  integer :: dist_max,stencil_radius(3)
  integer :: ireq_send(1),ireq_recv(1)
  real(8) :: flux_coef
  real(8),allocatable :: buf_send(:,:,:,:,:),buf_recv(:,:,:,:,:)
  logical :: owns_send_face,owns_target_face

  if(.not. allocated(psi%rwf)) return
  if(.not. allocated(hpsi%rwf)) return
  if(dc%id_frag < 0) return
  if(info%ik_s /= 1 .or. info%ik_e /= 1) stop "DC-LCFO-Flux CG: ik/=1"
  if(info%im_s /= 1 .or. info%im_e /= 1) stop "DC-LCFO-Flux CG: im/=1"
  if(.not. stencil%if_orthogonal) stop "DC-LCFO-Flux CG: nonorthogonal lattice is not supported"
  if(dc%optimized_fragment_geometry) stop "DC-LCFO-Flux CG: optimized fragment geometry is not supported"
  if(mod(dc%isize_tot,dc%n_frag) /= 0) stop "DC-LCFO-Flux CG: MPI size must be divisible by number of fragments"

  call get_fragment_domain(dc,dc%i_frag,ncore)
  do n=1,3
    stencil_radius(n) = active_laplacian_radius(stencil,n)
  end do
  nio = info%io_e - info%io_s + 1
  if(nio <= 0) return
  if(info%io_s < lbound(psi%rwf,5) .or. info%io_e > ubound(psi%rwf,5)) &
  & stop "DC-LCFO-Flux CG: psi orbital range mismatch"
  if(info%io_s < lbound(hpsi%rwf,5) .or. info%io_e > ubound(hpsi%rwf,5)) &
  & stop "DC-LCFO-Flux CG: hpsi orbital range mismatch"
  if(system%nspin > ubound(psi%rwf,4) .or. system%nspin > ubound(hpsi%rwf,4)) &
  & stop "DC-LCFO-Flux CG: spin range mismatch"

  npg_base = dc%isize_tot / dc%n_frag

  call prepare_dc_lcfo_flux_halo_cache(dc,ncore,stencil_radius,npg_base)
  n_halo = flux_halo_cache_n
  halo => flux_halo_cache

  do i_halo=1,n_halo
    l = halo(i_halo)%length
    allocate(buf_send(l(1),l(2),l(3),system%nspin,nio))
    allocate(buf_recv(l(1),l(2),l(3),system%nspin,nio))
    buf_send = 0d0
    buf_recv = 0d0
    lo(1) = max(1,lbound(psi%rwf,1) - halo(i_halo)%dsp_send(1))
    lo(2) = max(1,lbound(psi%rwf,2) - halo(i_halo)%dsp_send(2))
    lo(3) = max(1,lbound(psi%rwf,3) - halo(i_halo)%dsp_send(3))
    hi(1) = min(l(1),ubound(psi%rwf,1) - halo(i_halo)%dsp_send(1))
    hi(2) = min(l(2),ubound(psi%rwf,2) - halo(i_halo)%dsp_send(2))
    hi(3) = min(l(3),ubound(psi%rwf,3) - halo(i_halo)%dsp_send(3))

    if(.not. any(lo > hi)) then
!$omp parallel do private(io,ilo,ispin,iz,iy,ix,ixg,iyg,izg) schedule(static)
      do io=info%io_s,info%io_e
        ilo = io - info%io_s + 1
        do ispin=1,system%nspin
        do iz=lo(3),hi(3)
        do iy=lo(2),hi(2)
        do ix=lo(1),hi(1)
          ixg = halo(i_halo)%dsp_send(1) + ix
          iyg = halo(i_halo)%dsp_send(2) + iy
          izg = halo(i_halo)%dsp_send(3) + iz
          buf_send(ix,iy,iz,ispin,ilo) = &
          & psi%rwf(ixg,iyg,izg,ispin,io,1,1)
        end do
        end do
        end do
        end do
      end do
!$omp end parallel do
    end if

    owns_send_face = owns_flux_send_face(mg,ncore,halo(i_halo)%dvec,stencil_radius)
    owns_target_face = owns_flux_target_face(mg,ncore,halo(i_halo)%dvec,stencil_radius)

    ! Tag convention follows the receiver-side halo label:
    ! dvec points from the source fragment to the destination fragment.
    ! A receiver with the same dvec waits for source=current-dvec, so both
    ! endpoints use source_fragment + offset(dvec)*n_frag.
    itag_dir = halo_tag_offset(halo(i_halo)%dvec)
    itag_send = dc%i_frag + itag_dir*dc%n_frag
    itag_recv = halo(i_halo)%ifrag_src + itag_dir*dc%n_frag
    rank_send = comm_proc_null
    rank_recv = comm_proc_null
    if(owns_send_face) rank_send = flux_send_destination_rank(info,dc,halo(i_halo),npg_base)
    if(owns_target_face) rank_recv = flux_recv_source_rank(info,dc,halo(i_halo),npg_base)
    ireq_send(1) = comm_isend(buf_send,rank_send,itag_send,dc%icomm_tot)
    ireq_recv(1) = comm_irecv(buf_recv,rank_recv,itag_recv,dc%icomm_tot)
    call comm_wait_all(ireq_recv)
    call comm_wait_all(ireq_send)

    if(.not. owns_target_face) then
      deallocate(buf_send,buf_recv)
      cycle
    end if

    axis = halo(i_halo)%axis
    dist_max = min(stencil_radius(axis),min(ncore(axis),l(axis)))
    lo = 1
    hi = l
    do n=1,3
      if(n /= axis) then
        lo(n) = max(lo(n),mg%is(n))
        hi(n) = min(hi(n),mg%ie(n))
      end if
    end do
    if(halo(i_halo)%dvec(axis) > 0) then
      lo(axis) = max(1,mg%is(axis))
      hi(axis) = min(dist_max,mg%ie(axis))
    else
      lo(axis) = max(1,ncore(axis) - mg%ie(axis) + 1)
      hi(axis) = min(dist_max,ncore(axis) - mg%is(axis) + 1)
    end if
    if(any(lo > hi)) then
      deallocate(buf_send,buf_recv)
      cycle
    end if

    if(halo(i_halo)%dvec(axis) > 0) then
!$omp parallel do collapse(3) private(iz,iy,ix,ih,il,dist,flux_coef,io,ilo,ispin) schedule(static)
      do iz=lo(3),hi(3)
      do iy=lo(2),hi(2)
      do ix=lo(1),hi(1)
        ih = [ix,iy,iz]
        dist = ih(axis)
        il = ih
        il(axis) = dist
        flux_coef = -0.5d0 * stencil%coef_lap(dist,axis)
        do io=info%io_s,info%io_e
          ilo = io - info%io_s + 1
          do ispin=1,system%nspin
            hpsi%rwf(il(1),il(2),il(3),ispin,io,1,1) = &
            & hpsi%rwf(il(1),il(2),il(3),ispin,io,1,1) &
            & + flux_coef * buf_recv(ix,iy,iz,ispin,ilo)
          end do
        end do
      end do
      end do
      end do
!$omp end parallel do
    else
!$omp parallel do collapse(3) private(iz,iy,ix,ih,il,dist,flux_coef,io,ilo,ispin) schedule(static)
      do iz=lo(3),hi(3)
      do iy=lo(2),hi(2)
      do ix=lo(1),hi(1)
        ih = [ix,iy,iz]
        dist = ih(axis)
        il = ih
        il(axis) = ncore(axis) - dist + 1
        flux_coef = -0.5d0 * stencil%coef_lap(dist,axis)
        do io=info%io_s,info%io_e
          ilo = io - info%io_s + 1
          do ispin=1,system%nspin
            hpsi%rwf(il(1),il(2),il(3),ispin,io,1,1) = &
            & hpsi%rwf(il(1),il(2),il(3),ispin,io,1,1) &
            & + flux_coef * buf_recv(ix,iy,iz,ispin,ilo)
          end do
        end do
      end do
      end do
      end do
!$omp end parallel do
    end if

    deallocate(buf_send,buf_recv)
  end do

end subroutine apply_dc_lcfo_flux_hpsi_rwf

integer function flux_send_destination_rank(info,dc,halo,npg_base) result(rank)
  use structures
  implicit none
  type(s_parallel_info),intent(in) :: info
  type(s_dcdft),intent(in) :: dc
  type(halo_cg_info),intent(in) :: halo
  integer,intent(in) :: npg_base
  integer :: local_rank,face_side

  if(halo%dvec(halo%axis) > 0) then
    face_side = 0
  else
    face_side = info%nprgrid(halo%axis) - 1
  end if
  local_rank = rank_at_flux_face(info,halo%axis,face_side)
  rank = (halo%ifrag_dst - 1)*npg_base + local_rank
  if(rank < 0 .or. rank >= dc%isize_tot) stop "DC-LCFO-Flux CG: invalid send rank"
end function flux_send_destination_rank

integer function flux_recv_source_rank(info,dc,halo,npg_base) result(rank)
  use structures
  implicit none
  type(s_parallel_info),intent(in) :: info
  type(s_dcdft),intent(in) :: dc
  type(halo_cg_info),intent(in) :: halo
  integer,intent(in) :: npg_base
  integer :: local_rank,face_side

  if(halo%dvec(halo%axis) > 0) then
    face_side = info%nprgrid(halo%axis) - 1
  else
    face_side = 0
  end if
  local_rank = rank_at_flux_face(info,halo%axis,face_side)
  rank = (halo%ifrag_src - 1)*npg_base + local_rank
  if(rank < 0 .or. rank >= dc%isize_tot) stop "DC-LCFO-Flux CG: invalid recv rank"
end function flux_recv_source_rank

integer function rank_at_flux_face(info,axis,face_side) result(rank)
  use structures
  implicit none
  type(s_parallel_info),intent(in) :: info
  integer,intent(in) :: axis,face_side
  integer :: iaddr(5)

  iaddr = info%iaddress
  iaddr(axis) = face_side
  rank = info%imap(iaddr(1),iaddr(2),iaddr(3),iaddr(4),iaddr(5))
end function rank_at_flux_face

integer function active_laplacian_radius(stencil,axis) result(radius)
  use structures
  implicit none
  type(s_stencil),intent(in) :: stencil
  integer,intent(in) :: axis
  integer :: dist

  radius = 0
  do dist=1,size(stencil%coef_lap,1)
    if(abs(stencil%coef_lap(dist,axis)) > 0d0) radius = dist
  end do
end function active_laplacian_radius

logical function owns_flux_target_face(mg,ncore,dvec,stencil_radius) result(owns)
  use structures
  implicit none
  type(s_rgrid),intent(in) :: mg
  integer,intent(in) :: ncore(3),dvec(3),stencil_radius(3)
  integer :: axis,idir,dist,ix(3)

  owns = .false.
  axis = 0
  do idir=1,3
    if(dvec(idir) /= 0) axis = idir
  end do
  if(axis == 0) return

  do dist=1,min(stencil_radius(axis),ncore(axis))
    ix = mg%is
    if(dvec(axis) > 0) then
      ix(axis) = dist
    else
      ix(axis) = ncore(axis) - dist + 1
    end if
    if(all(ix >= mg%is) .and. all(ix <= mg%ie)) then
      owns = .true.
      return
    end if
  end do
end function owns_flux_target_face

logical function owns_flux_send_face(mg,ncore,dvec,stencil_radius) result(owns)
  use structures
  implicit none
  type(s_rgrid),intent(in) :: mg
  integer,intent(in) :: ncore(3),dvec(3),stencil_radius(3)
  integer :: axis,idir,dist,ix(3)

  owns = .false.
  axis = 0
  do idir=1,3
    if(dvec(idir) /= 0) axis = idir
  end do
  if(axis == 0) return

  do dist=1,min(stencil_radius(axis),ncore(axis))
    ix = mg%is
    if(dvec(axis) > 0) then
      ix(axis) = ncore(axis) - dist + 1
    else
      ix(axis) = dist
    end if
    if(all(ix >= mg%is) .and. all(ix <= mg%ie)) then
      owns = .true.
      return
    end if
  end do
end function owns_flux_send_face

integer function halo_tag_offset(dvec) result(offset)
  implicit none
  integer,intent(in) :: dvec(3)

  offset = (dvec(1) + 1)*9 + (dvec(2) + 1)*3 + (dvec(3) + 1)
end function halo_tag_offset

end module Conjugate_Gradient
