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

module salmonpy_gs_solver
use salmonpy_gs_storage, only: s_gs_storage
implicit none
private
public :: gs_solve
contains
subroutine gs_solve(ctx,release_requests,retain,import_only)
use math_constants, only: pi, zi
use structures
use inputoutput
use parallelization, only: nproc_id_global,nproc_group_global,adjust_elapse_time,nproc_size_global
use communication, only: comm_is_root, comm_summation, comm_bcast, comm_sync_all
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
use sendrecv_grid, only: dealloc_cache
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
use lcfo
use atomic_orbital, only: verify_atomic_orbital_consistency !!DEBUG
implicit none
logical,intent(in) :: release_requests,retain
logical,optional,intent(in) :: import_only
type(s_gs_storage),intent(inout) :: ctx
integer :: report_unit,report_status,report_flag
character(64) :: report_header
real(8) :: report_energy
integer :: ix,iy,iz
integer :: Miter,iatom,jj,nspin,force_status_unit
real(8) :: sum1
character(100) :: comment_line

logical :: rion_update
logical :: flag_opt_conv
integer :: Miopt, iopt,nopt_max,i
integer :: iter_band_kpt, iter_band_kpt_end, iter_band_kpt_stride
logical :: is_checkpoint_iter, is_shutdown_time
integer :: ilevel_print

if(theory=='dft_band'.and.iperiodic/=3) return

if(yn_dc=='y') then
  call init_dcdft(ctx%dc,ctx%pp,ctx%mixing,ctx%ewald)
  ilevel_print = 0
else
  ilevel_print = 3
end if

!check condition for using jellium model
if(yn_jm=='y') call check_condition_jm

call init_xc(ctx%xc_func, spin, cval, xcname=xc, xname=xname, cname=cname)
ctx%xc_live=.true.
ctx%active=.true.

call timer_begin(LOG_TOTAL)
call timer_begin(LOG_INIT_GS)


! please move folloings into initialization_dft
call init_dft(nproc_group_global,ctx%info,ctx%lg,ctx%mg,ctx%system,ctx%stencil,ctx%fg,ctx%poisson, &
  ctx%srg,ctx%srg_scalar,ctx%ofl,ctx%unfold)
allocate( ctx%rho_s(ctx%system%nspin),ctx%V_local(ctx%system%nspin),ctx%Vxc(ctx%system%nspin) )

call initialization1_dft( ctx%system, ctx%energy, ctx%stencil, ctx%fg, ctx%poisson,  &
                          ctx%lg, ctx%mg,   &
                          ctx%info,  &
                          ctx%srg, ctx%srg_scalar,  &
                          ctx%rho, ctx%rho_jm, ctx%rho_s, ctx%Vh, ctx%V_local, ctx%Vpsl, ctx%Vxc,  &
                          ctx%spsi, ctx%shpsi, ctx%sttpsi,  &
                          ctx%pp, ctx%ppg, ctx%ppn,  &
                          ctx%ofl )

call initialization2_dft( Miter, nspin, rion_update,  &
                          ctx%system, ctx%energy, ctx%ewald, ctx%stencil, ctx%fg, ctx%poisson,&
                          ctx%lg, ctx%mg, ctx%info,   &
                          ctx%srg, ctx%srg_scalar,  &
                          ctx%rho, ctx%rho_jm, ctx%rho_s, ctx%Vh,ctx%V_local, ctx%Vpsl, ctx%Vxc,  &
                          ctx%spsi, ctx%shpsi, ctx%sttpsi,  &
                          ctx%pp, ctx%ppg, ctx%ppn,   &
                          ctx%xc_func, ctx%mixing )

! Import existing checkpoint through the original initialization, without SCF or writes.
if(present(import_only))then
  if(import_only)then
    call timer_end(LOG_INIT_GS)
    call timer_end(LOG_TOTAL)
    if(release_requests)then
      call dealloc_cache(ctx%srg)
      call dealloc_cache(ctx%srg_scalar)
    endif
    ctx%iterations=Miter
    ctx%imported=.true.
    return
  endif
endif
Miopt = 0
nopt_max = 1
if(yn_opt=='y') call initialization_opt(Miopt,ctx%opt,ctx%system,flag_opt_conv,nopt_max,ctx%ofl)

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
  call dealloc_init_ps(ctx%ppg)
  call init_ps(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%fg,ctx%poisson,ctx%pp,ctx%ppg,ctx%Vpsl)
  call calc_nlcc(ctx%pp, ctx%system, ctx%mg, ctx%ppn)
  if(yn_auto_mixing=='y') call reset_mixing_rate(ctx%mixing)
  call timer_end(LOG_INIT_GS)
end if

!---------------------------------------- Band Iteration

if(theory=='dft_band')then
   call init_band_dft(ctx%system,ctx%band) ! --> system%wtk=0.0
   iter_band_kpt_end    = ctx%band%num_band_kpt
   iter_band_kpt_stride = ctx%system%nk
else
   iter_band_kpt_end    = 1
   iter_band_kpt_stride = 1
end if

call comm_sync_all
call timer_enable_sub
Band_Iteration : do iter_band_kpt= 1, iter_band_kpt_end, iter_band_kpt_stride

if(theory=='dft_band')then
   call calc_band_write(iter_band_kpt,ctx%system,ctx%band,ctx%info)
end if


call timer_begin(LOG_INIT_GS_ITERATION)

call timer_end(LOG_INIT_GS_ITERATION)


call timer_begin(LOG_GS_ITERATION)
!------------------------------------ SCF Iteration
!Iteration loop for SCF (DFT_Iteration)
call scf_iteration_dft( Miter,rion_update,sum1,  &
                        ctx%system,ctx%energy,ctx%ewald,  &
                        ctx%lg,ctx%mg,  &
                        ctx%info,  &
                        ctx%poisson,ctx%fg,  &
                        ctx%cg,ctx%mixing,  &
                        ctx%stencil,  &
                        ctx%srg,ctx%srg_scalar,   &
                        ctx%spsi,ctx%shpsi,ctx%sttpsi,  &
                        ctx%rho,ctx%rho_jm,ctx%rho_s,  &
                        ctx%V_local,ctx%Vh,ctx%Vxc,ctx%Vpsl,ctx%xc_func,  &
                        ctx%pp,ctx%ppg,ctx%ppn,  &
                        ctx%band, ilevel_print, ctx%dc)


if(theory=='dft_band')then
   call write_band(ctx%system,ctx%energy)
   call write_band_status(iter_band_kpt,ctx%system,ctx%band)
end if

! output the wavefunctions for next GS calculations
if(write_gs_wfn_k == 'y') then !this input keyword is going to be removed....
   select case(iperiodic)
   case(3)
      call write_wfn(ctx%lg,ctx%mg,ctx%spsi,ctx%info,ctx%system)
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
      call write_k_data(ctx%system,ctx%stencil)  !need? (probably remove later)
      call write_tm_data(ctx%spsi,ctx%system,ctx%info,ctx%mg,ctx%stencil,ctx%srg,ctx%ppg,ctx%energy)
   case(0)
     write(*,*) "error: yn_out_tm='y','yn_out_tm_bin='y',yn_out_gs_sgm_eps='y' & iperiodic=0"
  end select
end if

   ! force
   if(yn_jm=='n' .and. yn_dc=="n")then
     call calc_force(ctx%system,ctx%pp,ctx%fg,ctx%info,ctx%mg,ctx%stencil,ctx%poisson,ctx%srg,ctx%ppg,ctx%spsi,ctx%ewald,ctx%Vxc)
     if(comm_is_root(nproc_id_global))then
        if (theory=='dft' .and. yn_opt=='n' .and. yn_md=='n') then
           open(newunit=force_status_unit,file='gs_forces.dat',status='replace')
           write(force_status_unit,'(a)') 'SALMONPY_FORCES_V1'
           write(force_status_unit,*) natom
           do iatom=1,natom
              write(force_status_unit,'(3es26.17e3)') ctx%system%Force(:,iatom)
           end do
           close(force_status_unit)
        end if
        write(*,*) "===== force ====="
        do iatom=1,natom
           select case(unit_system)
           case('au','a.u.'); write(*,300)iatom,(ctx%system%Force(ix,iatom),ix=1,3)
           case('A_eV_fs'  ); write(*,300)iatom,(ctx%system%Force(ix,iatom)*au_energy_ev/au_length_aa,ix=1,3)
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
   call structure_opt_check(iopt,flag_opt_conv,ctx%system%Force)
   if(.not.flag_opt_conv) call structure_opt(ctx%opt,iopt,ctx%system)
   !! Rion is old variables to be removed
   !! but currently it is used in many subroutines.
   !!Rion(:,:) = system%Rion(:,:)

   write(comment_line,10) iopt
   call write_xyz(comment_line,"add","r  ",ctx%system,ctx%ofl)
10 format("#opt iteration step=",i5)

   if(comm_is_root(nproc_id_global))then
      write(*,*) "atomic coordinate"
      do iatom=1,natom
         write(*,20) "'"//trim(atom_name(iatom))//"'",  &
                   (ctx%system%Rion(jj,iatom)*ulength_from_au,jj=1,3), &
                   Kion(iatom), flag_opt_atom(iatom)
      end do
20    format(a5,3f16.8,i3,a3)
   end if

   if(flag_opt_conv) then
      call structure_opt_fin(ctx%opt)
   else
      is_checkpoint_iter = (checkpoint_interval >= 1) .and. (mod(iopt,checkpoint_interval) == 0)
      is_shutdown_time   = (time_shutdown > 0d0) .and. (adjust_elapse_time(timer_now(LOG_TOTAL)) > time_shutdown)

      if(is_checkpoint_iter .or. is_shutdown_time) then
         if (is_shutdown_time .and. comm_is_root(ctx%info%id_rko)) then
           print *, 'shutdown the calculation, iopt =', iopt
         end if

         call checkpoint_gs(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%spsi,iopt,ctx%mixing)
         call comm_sync_all
         call checkpoint_opt(iopt,ctx%opt)
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


!------------ Writing part -----------
call timer_begin(LOG_WRITE_GS_RESULTS)

if(yn_dc=='y') then
  if(yn_dc_lcfo == 'y') call dc_lcfo(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%stencil,ctx%ppg, &
    ctx%energy,ctx%v_local,ctx%spsi,ctx%shpsi,ctx%sttpsi,ctx%srg,ctx%dc)
  call write_total_dcdft(ctx%system,ctx%dc)
end if

! write GS: basic data
if(yn_dc=='n') call write_band_information(ctx%system,ctx%energy)
call write_eigen(ctx%ofl,ctx%system,ctx%energy)
call write_info_data(Miter,ctx%system,ctx%energy,ctx%pp)
call write_k_data(ctx%system,ctx%stencil)
if(yn_spinorbit=='y') call write_mag_decomposed_gs(ctx%system,ctx%mg,ctx%info,ctx%spsi)

! write GS: analysis option
if(yn_out_psi =='y') call write_psi(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%spsi)
if(yn_out_dns =='y') call write_dns(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%rho_s)
if(yn_out_dos =='y') call write_dos(ctx%system,ctx%energy)
if(yn_out_pdos=='y') call write_pdos(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%pp,ctx%energy,ctx%spsi)
if(yn_out_elf =='y') call write_elf(0,ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%stencil,ctx%rho,ctx%srg,ctx%srg_scalar,ctx%spsi)

! if(comm_is_root(nproc_id_global)) call verify_atomic_orbital_consistency(system,pp,lg,mg,nproc_id_global) !!DEBUG

call timer_end(LOG_WRITE_GS_RESULTS)

! write GS: binary data for restart
call timer_begin(LOG_WRITE_GS_DATA)
if(write_gs_restart_data=="no") then
   if(comm_is_root(nproc_id_global)) &
      write(*,'(a)')"  no restart data writing."
else if(write_gs_restart_data.ne."checkpoint_only") then
   if(comm_is_root(nproc_id_global)) write(*,'(a)')"  writing restart data..."
   call checkpoint_gs(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%spsi,Miter,ctx%mixing,ctx%ofl%dir_out_restart)
   call comm_sync_all
   if(yn_opt=='y') then
      if(.not.flag_opt_conv) then
         call comm_sync_all
         call checkpoint_opt(nopt_max,ctx%opt,ctx%ofl%dir_out_restart)
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
   call checkpoint_gs(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%spsi,Miter,ctx%mixing)
   call comm_sync_all
endif
if(comm_is_root(nproc_id_global)) write(*,'(a)')"  writing completed."
call timer_end(LOG_WRITE_GS_DATA)

!call timer_begin(LOG_WRITE_GS_INFO)  !if needed, please take back, sory: AY
!call timer_end(LOG_WRITE_GS_INFO)

if(.not.retain)then
  call finalize_xc(ctx%xc_func)
  ctx%xc_live=.false.
endif

if(yn_dc=='y') then
! override (restore)
  nproc_group_global = ctx%dc%icomm_tot
  nproc_id_global = ctx%dc%id_tot
  nproc_size_global = ctx%dc%isize_tot
  call comm_sync_all
end if

call timer_end(LOG_TOTAL)


! The model caller releases scoped communicators after these local buffers retire.
if(release_requests)then
  call dealloc_cache(ctx%srg)
  call dealloc_cache(ctx%srg_scalar)
endif
ctx%iterations=Miter
ctx%residual=sum1
ctx%threshold=threshold
ctx%report_status=0
if(retain)then
  ! Snapshot the existing native convergence flag; do not reinterpret its last-iteration policy.
  report_status=0;report_flag=0
  if(comm_is_root(nproc_id_global))then
    open(newunit=report_unit,file='gs_status.dat',status='old',iostat=report_status)
    if(report_status==0)then
      read(report_unit,'(a)',iostat=report_status)report_header
      if(report_status==0.and.trim(report_header)/='SALMONPY_GS_V1')report_status=1
      if(report_status==0)read(report_unit,*,iostat=report_status)report_flag,ctx%iterations
      if(report_status==0)read(report_unit,*,iostat=report_status)ctx%residual,ctx%threshold,report_energy
      close(report_unit)
    endif
  endif
  call comm_bcast(report_status,nproc_group_global)
  call comm_bcast(report_flag,nproc_group_global)
  ctx%report_status=report_status
  ctx%converged=report_flag==1
endif
end subroutine gs_solve
end module salmonpy_gs_solver

subroutine main_dft
  use salmonpy_gs_solver, only: gs_solve
  use salmonpy_gs_storage, only: s_gs_storage
  type(s_gs_storage) :: ctx
  call gs_solve(ctx,.false.,.false.)
end subroutine main_dft

subroutine main_dft_model
  use salmonpy_gs_solver, only: gs_solve
  use salmonpy_gs_storage, only: s_gs_storage
  type(s_gs_storage) :: ctx
  call gs_solve(ctx,.true.,.false.)
end subroutine main_dft_model
