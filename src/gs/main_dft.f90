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
!=======================================================================

#include "config.h"

subroutine main_dft
use math_constants, only: pi, zi
use structures
use inputoutput
use parallelization, only: nproc_id_global,nproc_group_global,adjust_elapse_time,nproc_size_global
use communication, only: comm_is_root, comm_summation, comm_bcast, comm_sync_all,comm_get_max
use salmon_xc
use timer
use scf_iteration_sub
use density_matrix, only: calc_density
use writefield
use salmon_pp, only: calc_nlcc
use hartree_sub, only: hartree
use force_sub
use write_sub
use read_gs
use code_optimization
use initialization_sub
use occupation
use prep_pp_sub
use mixing_sub
use checkpoint_restart_sub
use hamiltonian
use structure_opt_sub
use total_energy
use band_dft_sub
use init_gs, only: init_wf
use initialization_dft
use jellium, only: check_condition_jm
use dcdft
use dc_force, only: report_dc_frozen_force
use lcfo
#ifdef USE_HSE
use exx_native, only: exx_export_snapshot,exx_eigen_diagnostic_enabled,exx_export_eigen_pair,exx_refresh
use exx_native, only: exx_enabled,exx_pack
use exx_stationarity,only:spinor_stationarity
use subspace_diagonalization,only:ssdg
#endif
implicit none
integer :: ix,iy,iz
integer :: Miter,iatom,jj,nspin
logical :: pre_scf_root
real(8) :: sum1
character(100) :: comment_line

type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_parallel_info) :: info
type(s_sendrecv_grid) :: srg, srg_scalar
type(s_orbital) :: spsi,shpsi,sttpsi
type(s_dft_system) :: system
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_xc_functional) :: xc_func
type(s_scalar) :: rho,rho_jm,Vh,Vpsl
type(s_scalar),allocatable :: V_local(:),rho_s(:),Vxc(:)
type(s_reciprocal_grid) :: fg
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_dft_energy) :: energy
type(s_ewald_ion_ion) :: ewald
type(s_cg)     :: cg
type(s_mixing) :: mixing
type(s_ofile)  :: ofl
type(s_band_dft) ::band
type(s_opt) :: opt
type(s_dcdft) :: dc
type(s_unfold) :: unfold

logical :: rion_update
logical :: flag_opt_conv
integer :: Miopt, iopt,nopt_max,i
integer :: iter_band_kpt, iter_band_kpt_end, iter_band_kpt_stride
logical :: is_checkpoint_iter, is_shutdown_time
integer :: ilevel_print

if(theory=='dft_band'.and.iperiodic/=3) return

if(yn_dc=='y') then
  call init_dcdft(dc,pp,mixing,ewald)
  ilevel_print = 0
else
  ilevel_print = 3
end if

!check condition for using jellium model
if(yn_jm=='y') call check_condition_jm

! Set the temporary functional before initialization evaluates any potential.
exx_pre_scf_active=exx_pre_scf_threshold>0d0
pre_scf_root=comm_is_root(nproc_id_global)
if(yn_dc=='y')pre_scf_root=comm_is_root(dc%id_tot)
if(exx_pre_scf_active.and.pre_scf_root) &
  write(*,'(a,es14.6,a,i0)') 'EXX_PRE_SCF PBE stage: threshold=',exx_pre_scf_threshold, &
                          ' consecutive steps=',exx_pre_scf_steps
call init_xc(xc_func, spin, cval, xcname=xc, xname=xname, cname=cname)

call timer_begin(LOG_TOTAL)
call timer_begin(LOG_INIT_GS)


! please move folloings into initialization_dft
call init_dft(nproc_group_global,info,lg,mg,system,stencil,fg,poisson,srg,srg_scalar,ofl,unfold)
allocate( rho_s(system%nspin),V_local(system%nspin),Vxc(system%nspin) )

call initialization1_dft( system, energy, stencil, fg, poisson,  &
                          lg, mg,   &
                          info,  &
                          srg, srg_scalar,  &
                          rho, rho_jm, rho_s, Vh, V_local, Vpsl, Vxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,  &
                          ofl )

call initialization2_dft( Miter, nspin, rion_update,  &
                          system, energy, ewald, stencil, fg, poisson,&
                          lg, mg, info,   &
                          srg, srg_scalar,  &
                          rho, rho_jm, rho_s, Vh,V_local, Vpsl, Vxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,   &
                          xc_func, mixing )

#ifdef USE_HSE
! Explicit opt-in, read-only diagnostic: evaluate the initial physical
! Hamiltonian, print scalar metrics, and return before any SCF iteration.
block
  character(8) :: diagnostic_env
  integer :: diagnostic_flag,diagnostic_status,diagnostic_bad,diagnostic_total,diagnostic_stage
  integer :: diagnostic_saved_ncg,diagnostic_remaining
  complex(8),allocatable :: diagnostic_psi(:,:,:),diagnostic_hpsi(:,:,:)
  real(8) :: diagnostic_metrics(5),diagnostic_global(5)
  diagnostic_flag=0
  if(comm_is_root(nproc_id_global))then
    call get_environment_variable('SALMON_HSE_STATIONARITY',diagnostic_env,status=diagnostic_status)
    if(diagnostic_status==0.and.trim(diagnostic_env)=='1')diagnostic_flag=1
    if(diagnostic_status==0.and.trim(diagnostic_env)=='solve')diagnostic_flag=2
    if(diagnostic_status==0.and.trim(diagnostic_env)=='blocked')diagnostic_flag=3
  endif
  call comm_bcast(diagnostic_flag,nproc_group_global)
  if(diagnostic_flag>0)then
    if(.not.exx_enabled().or.yn_spinorbit/='y'.or.yn_dc/='n'.or. &
       info%isize_r/=1.or.info%isize_o/=1) &
      error stop 'HSE stationarity diagnostic requires native SOC k-only layout'
    call exx_refresh(system,mg,info,spsi)
    if(diagnostic_flag>=2.and.propagator/='hse_taylor4_full') &
      error stop 'Fixed-H stationarity solve requires direct exchange'
    allocate(diagnostic_psi(2*product(mg%num),system%no,info%numk), &
             diagnostic_hpsi(2*product(mg%num),system%no,info%numk))
    do diagnostic_stage=1,merge(3,1,diagnostic_flag>=2)
    ! No density, potential, occupation or Fock-source refresh in these passes.
    if(diagnostic_stage==2)then
      diagnostic_saved_ncg=ncg
      diagnostic_remaining=ncg
      do
        if(diagnostic_flag==3)ncg=min(8,diagnostic_remaining)
        call solve_orbitals(mg,system,info,stencil,spsi,shpsi,sttpsi,srg,cg,ppg,V_local, &
          2,nscf_init_no_diagonal)
        if(diagnostic_flag/=3)exit
        call ssdg(mg,system,info,stencil,spsi,shpsi,ppg,V_local,srg)
        diagnostic_remaining=diagnostic_remaining-ncg
        if(diagnostic_remaining<=0)exit
      enddo
      ncg=diagnostic_saved_ncg
    endif
    if(diagnostic_stage==3)call ssdg(mg,system,info,stencil,spsi,shpsi,ppg,V_local,srg)
    call hpsi(spsi,shpsi,info,mg,V_local,system,stencil,srg,ppg)
    call exx_pack(spsi,mg,info,diagnostic_psi)
    call exx_pack(shpsi,mg,info,diagnostic_hpsi)
    call spinor_stationarity(diagnostic_psi,diagnostic_hpsi, &
      system%rocc(:,info%ik_s:info%ik_e,1),system%hvol, &
      diagnostic_metrics(1),diagnostic_metrics(2),diagnostic_metrics(3),diagnostic_bad, &
      diagnostic_metrics(4),diagnostic_metrics(5))
    call comm_summation(diagnostic_bad,diagnostic_total,info%icomm_rko)
    if(diagnostic_total/=0)error stop 'Invalid HSE stationarity diagnostic data'
    call comm_get_max(diagnostic_metrics,diagnostic_global,5,info%icomm_rko)
    if(comm_is_root(nproc_id_global))then
      if(diagnostic_flag>=2)write(*,'(a,i0)') 'HSE_FIXED_H_STAGE initial/solver/finalRitz: ',diagnostic_stage
      write(*,'(a,5es24.16)') 'HSE_STATIONARITY gram_error/hermitian_error_Ha/commutator_Ha/outside_Ha/within_Ha: ',diagnostic_global
    endif
    enddo
    if(comm_is_root(nproc_id_global)) &
      write(*,'(a)') 'HSE_STATIONARITY diagnostic only; no SCF performed or state certified.'
    return
  endif
end block
#endif

Miopt = 0
nopt_max = 1
if(yn_opt=='y') call initialization_opt(Miopt,opt,system,flag_opt_conv,nopt_max,ofl)

call timer_end(LOG_INIT_GS)

!---------------------------------------- Opt Iteration


#ifdef __FUJITSU
call fipp_start ! performance profiling
#endif

Structure_Optimization_Iteration : do iopt= Miopt+1, nopt_max

if(iopt>=2)then
  call timer_begin(LOG_INIT_GS)
  Miter = 0        ! Miter: Iteration counter set to zero
  rion_update = .true.
  call dealloc_init_ps(ppg)
  call init_ps(lg,mg,system,info,fg,poisson,pp,ppg,Vpsl)
  call calc_nlcc(pp, system, mg, ppn)
  if(yn_auto_mixing=='y') call reset_mixing_rate(mixing)
  call timer_end(LOG_INIT_GS)
end if

!---------------------------------------- Band Iteration

if(theory=='dft_band')then
   call init_band_dft(system,band) ! --> system%wtk=0.0
   iter_band_kpt_end    = band%num_band_kpt
   iter_band_kpt_stride = system%nk
else
   iter_band_kpt_end    = 1
   iter_band_kpt_stride = 1
end if

call comm_sync_all
call timer_enable_sub
Band_Iteration : do iter_band_kpt= 1, iter_band_kpt_end, iter_band_kpt_stride

if(theory=='dft_band')then
   call calc_band_write(iter_band_kpt,system,band,info)
end if


call timer_begin(LOG_INIT_GS_ITERATION)

call timer_end(LOG_INIT_GS_ITERATION)


call timer_begin(LOG_GS_ITERATION)
!------------------------------------ SCF Iteration
!Iteration loop for SCF (DFT_Iteration)
call scf_iteration_dft( Miter,rion_update,sum1,  &
                        system,energy,ewald,  &
                        lg,mg,  &
                        info,  &
                        poisson,fg,  &
                        cg,mixing,  &
                        stencil,  &
                        srg,srg_scalar,   &
                        spsi,shpsi,sttpsi,  &
                        rho,rho_jm,rho_s,  &
                        V_local,Vh,Vxc,Vpsl,xc_func,  &
                        pp,ppg,ppn,  &
                        band, ilevel_print, dc)


if(theory=='dft_band')then
   call write_band(system,energy)
end if

! output the wavefunctions for next GS calculations
if(write_gs_wfn_k == 'y') then !this input keyword is going to be removed....
   select case(iperiodic)
   case(3)
      call write_wfn(lg,mg,spsi,info,system)
      ! Experimental Implementation of Inner-Product Outputs:
      ! call write_prod_dk_data(lg, mg, system, info, spsi)
   case(0)
      write(*,*) "error: write_gs_wfn_k='y' & iperiodic=0"
   end select
end if

! output transition moment : --> want to put out of the optmization loop in future
if(yn_out_tm  == 'y'.or. yn_out_tm_bin == 'y'.or.yn_out_gs_sgm_eps=='y') then
   select case(iperiodic)
   case(3)
      call write_k_data(system,stencil)  !need? (probably remove later)
      call write_tm_data(spsi,system,info,mg,stencil,srg,ppg,energy)
   case(0)
     write(*,*) "error: yn_out_tm='y','yn_out_tm_bin='y',yn_out_gs_sgm_eps='y' & iperiodic=0"
  end select
end if

   if(yn_dc_force_diagnostic=='y')then
     if(.not.(sum1<threshold))error stop 'DC force diagnostic: SCF not converged'
     call report_dc_frozen_force(dc,system,info,mg,pp,ppg,spsi,ewald,energy)
   endif
   ! force
   if(yn_jm=='n' .and. yn_dc=="n")then
     call calc_force(system,pp,fg,info,mg,stencil,poisson,srg,ppg,spsi,ewald)
     if(comm_is_root(nproc_id_global))then
        write(*,*) "===== force ====="
        do iatom=1,natom
           select case(unit_system)
           case('au','a.u.'); write(*,300)iatom,(system%Force(ix,iatom),ix=1,3)
           case('A_eV_fs'  ); write(*,300)iatom,(system%Force(ix,iatom)*au_energy_ev/au_length_aa,ix=1,3)
           end select
        end do
300   format(i6,3e16.8)
     end if
   end if

call timer_end(LOG_GS_ITERATION)

end do Band_Iteration
call timer_disable_sub


call timer_begin(LOG_DEINIT_GS_ITERATION)
if(yn_opt=='y') then
   call structure_opt_check(iopt,flag_opt_conv,system%Force)
   if(.not.flag_opt_conv) call structure_opt(opt,iopt,system)
   !! Rion is old variables to be removed
   !! but currently it is used in many subroutines.
   !!Rion(:,:) = system%Rion(:,:)

   write(comment_line,10) iopt
   call write_xyz(comment_line,"add","r  ",system,ofl)
10 format("#opt iteration step=",i5)

   if(comm_is_root(nproc_id_global))then
      write(*,*) "atomic coordinate"
      do iatom=1,natom
         write(*,20) "'"//trim(atom_name(iatom))//"'",  &
                   (system%Rion(jj,iatom)*ulength_from_au,jj=1,3), &
                   Kion(iatom), flag_opt_atom(iatom)
      end do
20    format(a5,3f16.8,i3,a3)
   end if

   if(flag_opt_conv) then
      call structure_opt_fin(opt)
   else
      is_checkpoint_iter = (checkpoint_interval >= 1) .and. (mod(iopt,checkpoint_interval) == 0)
      is_shutdown_time   = (time_shutdown > 0d0) .and. (adjust_elapse_time(timer_now(LOG_TOTAL)) > time_shutdown)

      if(is_checkpoint_iter .or. is_shutdown_time) then
         if (is_shutdown_time .and. comm_is_root(info%id_rko)) then
           print *, 'shutdown the calculation, iopt =', iopt
         end if

         call checkpoint_gs(lg,mg,system,info,spsi,iopt,mixing)
         call comm_sync_all
         call checkpoint_opt(iopt,opt)
         if(comm_is_root(nproc_id_global))then
            write(*,'(a,i5)')"  checkpoint data is printed: iopt=", iopt
         endif
         call comm_sync_all

         if (is_shutdown_time) then
           exit Structure_Optimization_Iteration
         end if
      endif
   endif

end if
call timer_end(LOG_DEINIT_GS_ITERATION)


if(yn_opt=='y')then
  if(flag_opt_conv)then
  exit Structure_Optimization_Iteration
  end if
end if
end do Structure_Optimization_Iteration

#ifdef __FUJITSU
call fipp_stop ! performance profiling
#endif


#ifdef USE_HSE
call exx_export_snapshot(system,mg,info,spsi,Miter,sum1,sum1<threshold)
if(exx_eigen_diagnostic_enabled(info))then
  call exx_refresh(system,mg,info,spsi)
  call hpsi(spsi,shpsi,info,mg,V_local,system,stencil,srg,ppg,sttpsi)
  call exx_export_eigen_pair(system,mg,info,spsi,shpsi)
end if
#endif
!------------ Writing part -----------
call timer_begin(LOG_WRITE_GS_RESULTS)

if(yn_dc=='y') then
  if(yn_dc_lcfo == 'y') call dc_lcfo(lg,mg,system,info,stencil,ppg,energy,v_local,spsi,shpsi,sttpsi,srg,dc)
  call write_total_dcdft(system,dc)
end if

! write GS: basic data
if(yn_dc=='n') call write_band_information(system,energy)
call write_eigen(ofl,system,energy)
call write_info_data(Miter,system,energy,pp)
call write_k_data(system,stencil)
if(yn_spinorbit=='y') call write_mag_decomposed_gs(system,mg,info,spsi)

! write GS: analysis option
if(yn_out_psi =='y') call write_psi(lg,mg,system,info,spsi)
if(yn_out_dns =='y') call write_dns(lg,mg,system,info,rho_s)
if(yn_out_dos =='y') call write_dos(system,energy)
if(yn_out_pdos=='y') call write_pdos(lg,mg,system,info,pp,energy,spsi)
if(yn_out_elf =='y') call write_elf(0,lg,mg,system,info,stencil,rho,srg,srg_scalar,spsi)

call timer_end(LOG_WRITE_GS_RESULTS)

! write GS: binary data for restart
call timer_begin(LOG_WRITE_GS_DATA)
if(write_gs_restart_data=="no") then
   if(comm_is_root(nproc_id_global)) &
      write(*,'(a)')"  no restart data writing."
else if(write_gs_restart_data.ne."checkpoint_only") then
   if(comm_is_root(nproc_id_global)) write(*,'(a)')"  writing restart data..."
   call checkpoint_gs(lg,mg,system,info,spsi,Miter,mixing,ofl%dir_out_restart)
   call comm_sync_all
   if(yn_opt=='y') then
      if(.not.flag_opt_conv) then
         call comm_sync_all
         call checkpoint_opt(nopt_max,opt,ofl%dir_out_restart)
         call comm_sync_all
      endif
   endif
else
   if(yn_self_checkpoint=='n') then
      if(comm_is_root(nproc_id_global)) then
           write(*,'(a)')"  no restart data writing:"
           write(*,'(a)')"  check input keywords if you need restart data"
       endif
   endif
endif
if(yn_self_checkpoint=='y') then
   if(comm_is_root(nproc_id_global)) &
   write(*,'(a)')"  writing restart data in checkpoint format ..."
   call checkpoint_gs(lg,mg,system,info,spsi,Miter,mixing)
   call comm_sync_all
endif
if(comm_is_root(nproc_id_global)) write(*,'(a)')"  writing completed."
call timer_end(LOG_WRITE_GS_DATA)

!call timer_begin(LOG_WRITE_GS_INFO)  !if needed, please take back, sory: AY
!call timer_end(LOG_WRITE_GS_INFO)

call finalize_xc(xc_func)

if(yn_dc=='y') then
! override (restore)
  nproc_group_global = dc%icomm_tot
  nproc_id_global = dc%id_tot
  nproc_size_global = dc%isize_tot
  call comm_sync_all
end if

call timer_end(LOG_TOTAL)


end subroutine main_dft
