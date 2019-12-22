!
!  Copyright 2019 SALMON developers
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


subroutine main_dft_md
use structures
use inputoutput
use parallelization, only: nproc_id_global,nproc_group_global
use communication, only: comm_is_root, comm_summation, comm_bcast
use salmon_xc
use timer
use scf_iteration_sub
use density_matrix, only: calc_density
use writefield
use salmon_pp, only: calc_nlcc
use hartree_sub, only: hartree
use force_sub
use write_sub, only: write_dft_md_data,write_xyz,write_band_information,write_eigen,write_info_data,write_dos,write_pdos
use read_gs
use code_optimization
use initialization_sub
use occupation
use input_pp_sub
use prep_pp_sub
use mixing_sub
use checkpoint_restart_sub
use hamiltonian
use salmon_total_energy
use init_gs, only: init_wf
use md_sub, only: init_md, update_pseudo_rt, &
                  time_evolution_step_md_part1,time_evolution_step_md_part2
                  
implicit none
type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_rgrid) :: ng
type(s_process_info) :: pinfo
type(s_orbital_parallel) :: info
type(s_field_parallel) :: info_field
type(s_sendrecv_grid) :: srg, srg_ng
type(s_orbital) :: spsi,shpsi,sttpsi
type(s_dft_system) :: system
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_xc_functional) :: xc_func
type(s_scalar) :: srho,sVh,sVpsl,rho_old,Vlocal_old
type(s_scalar),allocatable :: V_local(:),srho_s(:),sVxc(:)
type(s_reciprocal_grid) :: fg
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_dft_energy) :: energy
type(s_cg)     :: cg
type(s_mixing) :: mixing
type(s_band_dft) ::band
type(s_ofile) :: ofl
type(s_md) :: md

logical :: rion_update
integer :: ix,iy,iz, nspin, it, Miter
real(8) :: dt_h, sum1
character(100) :: file_atoms_coo, comment_line


call init_xc(xc_func, ispin, cval, xcname=xc, xname=xname, cname=cname)

call timer_begin(LOG_TOTAL)
call timer_begin(LOG_INIT_GS)

it=0

! please move folloings into initialization_dft 
call init_dft(nproc_group_global,pinfo,info,info_field,lg,mg,ng,system,stencil,fg,poisson,srg,srg_ng,ofl)
allocate( srho_s(system%nspin),V_local(system%nspin),sVxc(system%nspin) )


call initialization1_dft( system, energy, stencil, fg, poisson,  &
                          lg, mg, ng,  &
                          pinfo, info, info_field,  &
                          srg, srg_ng,  &
                          srho, srho_s, sVh, V_local, sVpsl, sVxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,  &
                          ofl )

call initialization2_dft( it, nspin, rion_update,  &
                          system, energy, stencil, fg, poisson,  &
                          lg, mg, ng,  &
                          info, info_field,   &
                          srg, srg_ng,  &
                          srho, srho_s, sVh,V_local, sVpsl, sVxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,   &
                          xc_func, mixing )

call initialization_dft_md( it, rion_update,  &
                          system, md, energy, stencil, fg, poisson,  &
                          lg, mg, ng,  &
                          info, info_field,   &
                          srg, srg_ng,  &
                          srho, srho_s, sVh,V_local, sVpsl, sVxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,   &
                          xc_func, mixing )

call timer_end(LOG_INIT_GS)

!-------------  MD Loop ------------------

!Export initial step
call write_dft_md_data(0,ofl,md)

write(comment_line,10) 0
call write_xyz(comment_line,"new","r  ",system)
10 format("#Adiabatic-MD time step=",i5)

dt_h       = dt*0.5d0

if(comm_is_root(nproc_id_global)) then
   write(*,'(3a)')"# 1.time[au], 2.Tene(ion)[au], 3.Uene(gs)[au], 4.Ework[au], ",&
                    "5.Etot[au], 6.Etot-Etot0[au], 7.Enh[au], 8.Htot[au], ", &
                    "9.Temperature_ion[K], 10.convergency[au]"
endif

rion_update = .true.

call timer_begin(LOG_GS_ITERATION)

MD_Loop : do it=1,nt

   call time_evolution_step_md_part1(it,system,md)

   call update_pseudo_rt(it,info,info_field,system,lg,mg,ng,poisson,fg,pp,ppg,ppn,sVpsl)

   if(allocated(rho_old%f))    deallocate(rho_old%f)
   if(allocated(Vlocal_old%f)) deallocate(Vlocal_old%f)
   call allocate_scalar(ng,rho_old)
   call allocate_scalar(ng,Vlocal_old)

!$OMP parallel do private(iz,iy,ix)
   do iz=ng%is(3),ng%ie(3)
   do iy=ng%is(2),ng%ie(2)
   do ix=ng%is(1),ng%ie(1)
      rho_old%f(ix,iy,iz)   = srho%f(ix,iy,iz)
      Vlocal_old%f(ix,iy,iz)= V_local(1)%f(ix,iy,iz)
   end do
   end do
   end do

   !-------------- SCF Iteration ----------------
   !Iteration loop for SCF (DFT_Iteration)
   Miter=0
   call scf_iteration_dft( Miter,rion_update,sum1,  &
                           system,energy,  &
                           lg,mg,ng,  &
                           info,info_field,pinfo,  &
                           poisson,fg,  &
                           cg,mixing,  &
                           stencil,  &
                           srg,srg_ng,   &
                           spsi,shpsi,sttpsi,  &
                           srho,srho_s,  &
                           V_local,sVh,sVxc,sVpsl,xc_func,  &
                           pp,ppg,ppn,  &
                           rho_old,Vlocal_old,  &
                           band,1 )

   ! force
   call calc_force(system,pp,fg,info,mg,stencil,srg,ppg,spsi)

   call time_evolution_step_md_part2(system,md)

   md%Uene  = energy%E_tot
   md%E_tot = energy%E_tot + md%Tene
   md%Htot  = energy%E_tot + md%E_nh  !for NHC

   !---Write section---
   ! Export to standard log file
   if(comm_is_root(nproc_id_global)) then
      write(*,120) it*dt, md%Tene,md%Uene,md%E_work,md%E_tot,md%E_tot-md%E_tot0,&
                   md%E_nh, md%Htot, md%Temperature, sum1
120   format(1x,f10.4, 7e20.10E3,f12.3,e15.8)
   endif

   ! Export to file_dft_md_data
   call write_dft_md_data(it,ofl,md)

   ! Export to file_trj
   if (yn_out_rvf_rt=='y' .and. mod(it,out_rvf_rt_step)==0)then
      write(comment_line,110) it, it*dt
110   format("#md-gs  step=",i8,"   time",e16.6)
      if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
           write(comment_line,112) trim(comment_line), md%xi_nh
112   format(a,"  xi_nh=",e18.10)
      call write_xyz(comment_line,"add","rvf",system)
   endif

   ! Export electronic density (cube or vtk)
   if(yn_out_dns_rt=='y' .and. mod(it,out_dns_rt_step)==0) then
      call write_dns(lg,mg,ng,srho%f,system%hgs,srho%f,it)
   end if

end do MD_Loop

call timer_end(LOG_GS_ITERATION)

!------------ Writing part -----------
call timer_begin(LOG_WRITE_GS_RESULTS)

!xxxxxxxxxxxx
!  call print_restart_data_md_gs
!
!  call comm_sync_all
!  if(comm_is_root(nproc_id_global)) write(*,*) 'This calculation is shutdown successfully!'
!xxxxxxxxxxxx

! write GS: basic data
!call write_band_information(system,energy)
!call write_eigen(ofl_eigen,system,energy)
!call write_info_data(it,system,energy,pp)

! write GS: analysis option
!if(yn_out_psi =='y') call write_psi(lg,mg,system,info,spsi)
!if(yn_out_dns =='y') call write_dns(lg,mg,ng,srho%f,system%hgs)
!if(yn_out_dos =='y') call write_dos(system,energy)
!if(yn_out_pdos=='y') call write_pdos(lg,mg,system,info,pp,energy,spsi)
!if(yn_out_elf =='y') call write_elf(0,lg,mg,ng,system,info,stencil,srho,srg,srg_ng,spsi)

call timer_end(LOG_WRITE_GS_RESULTS)

! write GS: binary data for restart
call timer_begin(LOG_WRITE_GS_DATA)
call checkpoint_gs(lg,mg,ng,system,info,spsi,it,mixing,ofl%dir_out_restart)
call timer_end(LOG_WRITE_GS_DATA)

!call timer_begin(LOG_WRITE_GS_INFO)  !if needed, please take back, sory: AY
!call timer_end(LOG_WRITE_GS_INFO)

call finalize_xc(xc_func)

call timer_end(LOG_TOTAL)


end subroutine main_dft_md
