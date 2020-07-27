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
use parallelization, only: nproc_id_global,nproc_group_global,adjust_elapse_time
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
use jellium, only: check_condition_jm, make_rho_jm
implicit none
integer :: ix,iy,iz
integer :: Miter,iatom,jj,nspin
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
type(s_scalar) :: rho,rho_jm,Vh,Vpsl,rho_old,Vlocal_old
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

logical :: rion_update
logical :: flag_opt_conv
integer :: Miopt, iopt,nopt_max
integer :: iter_band_kpt, iter_band_kpt_end, iter_band_kpt_stride
logical :: is_checkpoint_iter, is_shutdown_time

if(theory=='dft_band'.and.iperiodic/=3) return

!check condition for using jellium model
if(yn_jm=='y') call check_condition_jm

call init_xc(xc_func, spin, cval, xcname=xc, xname=xname, cname=cname)

call timer_begin(LOG_TOTAL)
call timer_begin(LOG_INIT_GS)

! please move folloings into initialization_dft
call init_dft(nproc_group_global,info,lg,mg,system,stencil,fg,poisson,srg,srg_scalar,ofl)
allocate( rho_s(system%nspin),V_local(system%nspin),Vxc(system%nspin) )

call initialization1_dft( system, energy, stencil, fg, poisson,  &
                          lg, mg,   &
                          info,  &
                          srg, srg_scalar,  &
                          rho, rho_s, Vh, V_local, Vpsl, Vxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,  &
                          ofl )

if(yn_jm=='n')then
  call initialization2_dft( Miter, nspin, rion_update,  &
                            system, energy, ewald, stencil, fg, poisson,&
                            lg, mg, info,   &
                            srg, srg_scalar,  &
                            rho, rho_s, Vh,V_local, Vpsl, Vxc,  &
                            spsi, shpsi, sttpsi,  &
                            pp, ppg, ppn,   &
                            xc_func, mixing )
else
  !make positive back ground charge density for using jellium model
  call allocate_scalar(mg,rho_jm)
  call make_rho_jm(lg,mg,system,info,rho_jm)
  
  call initialization2_dft( Miter, nspin, rion_update,  &
                            system, energy, ewald, stencil, fg, poisson,&
                            lg, mg, info,   &
                            srg, srg_scalar,  &
                            rho, rho_s, Vh,V_local, Vpsl, Vxc,  &
                            spsi, shpsi, sttpsi,  &
                            pp, ppg, ppn,   &
                            xc_func, mixing,  &
                            rho_jm )
end if

Miopt = 0
nopt_max = 1
if(yn_opt=='y') call initialization_opt(Miopt,opt,system,flag_opt_conv,nopt_max)

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
   call init_band_dft(system,band)
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

if(allocated(rho_old%f))    deallocate(rho_old%f)
if(allocated(Vlocal_old%f)) deallocate(Vlocal_old%f)
call allocate_scalar(mg,rho_old)
call allocate_scalar(mg,Vlocal_old)

!$OMP parallel do private(iz,iy,ix)
do iz=mg%is(3),mg%ie(3)
do iy=mg%is(2),mg%ie(2)
do ix=mg%is(1),mg%ie(1)
   rho_old%f(ix,iy,iz)   = rho%f(ix,iy,iz)
   Vlocal_old%f(ix,iy,iz)= V_local(1)%f(ix,iy,iz)
end do
end do
end do

call timer_end(LOG_INIT_GS_ITERATION)


call timer_begin(LOG_GS_ITERATION)
!------------------------------------ SCF Iteration
!Iteration loop for SCF (DFT_Iteration)
if(yn_jm=='n')then
  call scf_iteration_dft( Miter,rion_update,sum1,  &
                          system,energy,ewald,  &
                          lg,mg,  &
                          info,  &
                          poisson,fg,  &
                          cg,mixing,  &
                          stencil,  &
                          srg,srg_scalar,   &
                          spsi,shpsi,sttpsi,  &
                          rho,rho_s,  &
                          V_local,Vh,Vxc,Vpsl,xc_func,  &
                          pp,ppg,ppn,  &
                          rho_old,Vlocal_old,  &
                          band, 3)
else
  call scf_iteration_dft( Miter,rion_update,sum1,  &
                          system,energy,ewald,  &
                          lg,mg,  &
                          info,  &
                          poisson,fg,  &
                          cg,mixing,  &
                          stencil,  &
                          srg,srg_scalar,   &
                          spsi,shpsi,sttpsi,  &
                          rho,rho_s,  &
                          V_local,Vh,Vxc,Vpsl,xc_func,  &
                          pp,ppg,ppn,  &
                          rho_old,Vlocal_old,  &
                          band, 3,  &
                          rho_jm )
end if


if(theory=='dft_band')then
   call write_band(system,energy)
end if

! output the wavefunctions for next GS calculations
if(write_gs_wfn_k == 'y') then !this input keyword is going to be removed....
   select case(iperiodic)
   case(3)
      call write_wfn(lg,mg,spsi,info,system)
      ! Experimental Implementation of Inner-Product Outputs:
      call write_prod_dk_data(lg, mg, system, info, spsi)
   case(0)
      write(*,*) "error: write_gs_wfn_k='y' & iperiodic=0"
   end select
end if

! output transition moment : --> want to put out of the optmization loop in future
if(yn_out_tm  == 'y') then
   select case(iperiodic)
   case(3)
      call write_k_data(system,stencil)  !need? (probably remove later)
      call write_tm_data(spsi,system,info,mg,stencil,srg,ppg)
   case(0)
     write(*,*) "error: yn_out_tm='y' & iperiodic=0"
  end select
end if

   ! force
   if(yn_jm=='n')then
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
   call write_xyz(comment_line,"add","r  ",system)
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


!------------ Writing part -----------
call timer_begin(LOG_WRITE_GS_RESULTS)

! write GS: basic data
call write_band_information(system,energy)
call write_eigen(ofl,system,energy)
call write_info_data(Miter,system,energy,pp)
call write_k_data(system,stencil)

! write GS: analysis option
if(yn_out_psi =='y') call write_psi(lg,mg,system,info,spsi)
if(yn_out_dns =='y') call write_dns(lg,mg,system,rho%f)
if(yn_out_dos =='y') call write_dos(system,energy)
if(yn_out_pdos=='y') call write_pdos(lg,mg,system,info,pp,energy,spsi)
if(yn_out_elf =='y') call write_elf(0,lg,mg,system,info,stencil,rho,srg,srg_scalar,spsi)

call timer_end(LOG_WRITE_GS_RESULTS)

! write GS: binary data for restart
call timer_begin(LOG_WRITE_GS_DATA)
call checkpoint_gs(lg,mg,system,info,spsi,Miter,mixing,ofl%dir_out_restart)
if(yn_opt=='y') then
   if(.not.flag_opt_conv) then
      call comm_sync_all
      call checkpoint_opt(nopt_max,opt,ofl%dir_out_restart)
      call comm_sync_all
   endif
endif
call timer_end(LOG_WRITE_GS_DATA)

!call timer_begin(LOG_WRITE_GS_INFO)  !if needed, please take back, sory: AY
!call timer_end(LOG_WRITE_GS_INFO)

call finalize_xc(xc_func)

call timer_end(LOG_TOTAL)


end subroutine main_dft
