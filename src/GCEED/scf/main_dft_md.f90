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
use salmon_parallel, only: nproc_id_global,nproc_group_global
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use salmon_xc
use timer
use scf_iteration_sub
use density_matrix, only: calc_density
use writefield
use global_variables_scf
use salmon_pp, only: calc_nlcc
use hartree_sub, only: hartree
use force_sub
use write_sub
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
type(s_ofile)  :: ofile
type(s_band_dft) ::band
type(s_rt) :: rt
type(s_md) :: md

integer :: ix,iy,iz, iatom,nspin
real(8) :: sum1
character(100) :: file_atoms_coo, comment_line

logical :: rion_update
integer :: it, Miter
real(8),allocatable :: FionE(:,:)

real(8) :: dt_h, Uene0, E_tot0
real(8) :: Htot, Enh, Enh_gkTlns, gkT, Qnh  !NHC xxx 


call init_xc(xc_func, ispin, cval, xcname=xc, xname=xname, cname=cname)

iSCFRT=1
iblacsinit=0

call timer_begin(LOG_TOTAL)
call timer_begin(LOG_INIT_GS)

call convert_input_scf(file_atoms_coo)
mixing%num_rho_stock = 21


! please move folloings into initialization_dft 
call init_dft(iSCFRT,nproc_group_global,pinfo,info,info_field,lg,mg,ng,system,stencil,fg,poisson,srg,srg_ng,ofile)
allocate( srho_s(system%nspin),V_local(system%nspin),sVxc(system%nspin) )


call initialization1_dft( system, energy, stencil, fg, poisson,  &
                          lg, mg, ng,  &
                          pinfo, info, info_field,  &
                          srg, srg_ng,  &
                          srho, srho_s, sVh, V_local, sVpsl, sVxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,  &
                          ofile )

call initialization2_dft( it, nspin, rion_update,  &
                          system, energy, stencil, fg, poisson,  &
                          lg, mg, ng,  &
                          info, info_field,   &
                          srg, srg_ng,  &
                          srho, srho_s, sVh,V_local, sVpsl, sVxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,   &
                          xc_func, mixing )


it=0
call init_A(nt,it,rt)
!electric field:---- not yet
!if(iperiodic==3) call calc_Aext(Mit)
!for isolated system, use Vbox ...

call init_md(system,md)


call timer_end(LOG_INIT_GS)

!-------------  MD Loop ------------------

write(comment_line,10) 0
call write_xyz(comment_line,"new","r  ",system)
10 format("#Adiabatic-MD time step=",i5)

dt_h       = dt*0.5d0

!  Enh_gkTlns = 0d0
!  Enh        = 0d0
!  if(ensemble=="NVT" .and. thermostat=="nose-hoover") then
!     gkT = 3d0*NI * kB/hartree2J*temperature0_ion_k
!     Qnh = gkT * thermostat_tau**2d0
!  endif

allocate( FionE(3,system%nion) )

if(comm_is_root(nproc_id_global)) then
   write(*,'(3a)')"# 1.time[au], 2.Tene(ion)[au], 3.Uene(gs)[au], 4.Ework[au], ",&
                    "5.Etot[au], 6.Etot-Etot0[au], 7.Enh[au], 8.Htot[au], ", &
                    "9.Temperature_ion[K], 10.convergency[au]"
endif

rion_update = .true.

call timer_begin(LOG_GS_ITERATION)

MD_Loop : do it=1,nt

   call time_evolution_step_md_part1(it,system,md)

   call update_pseudo_rt(it,info,info_field,system,stencil,lg,mg,ng,poisson,fg,pp,ppg,ppn,sVpsl)


   poisson%iterVh=1000   ! what's this? necessary?
   iflag_diisjump=0

   if(.not.allocated(norm_diff_psi_stock)) then
      if(it==1) allocate(norm_diff_psi_stock(itotMST,1))
   end if
   norm_diff_psi_stock = 1d9

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
                           info,info_field,  &
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

!   ! output the wavefunctions for next GS calculations
!   if(write_gs_wfn_k == 'y') then !this input keyword is going to be removed....
!      select case(iperiodic)
!      case(3)
!         call write_wfn(lg,mg,spsi,info,system)
!         ! Experimental Implementation of Inner-Product Outputs:
!         call write_prod_dk_data(lg, mg, system, info, spsi) 
!      case(0)
!         write(*,*) "error: write_gs_wfn_k='y' & iperiodic=0"
!      end select
!   end if

   !! output transition moment (currently no support in dft_md)
   !if(yn_out_tm == 'y') then
   !   select case(iperiodic)
   !   case(3)
   !      call write_k_data(system,stencil)
   !      call write_tm_data(spsi,system,info,mg,stencil,srg,ppg)
   !   case(0)
   !      write(*,*) "error: yn_out_tm='y' & iperiodic=0"
   !   end select
   !end if

   ! force
   if(iperiodic == 3 .and. yn_ffte=='y') then
      !NOTE: calc_force_salmon hangs under this configuration due to ppg%vpsl_atom
      !      does not allocate.
   else
      call calc_force_salmon(system,pp,fg,info,mg,stencil,srg,ppg,spsi)
   end if

   !force on ion directly from field --- should put in calc_force_salmon?
   do iatom=1,system%nion
      FionE(:,iatom) = pp%Zps(Kion(iatom)) * rt%E_tot(:,it)
   enddo
   system%Force(:,:) = system%Force(:,:) + FionE(:,:)


   call time_evolution_step_md_part2(system,md)

   if(it==1) then
      Uene0 = energy%E_tot
      E_tot0= energy%E_tot + md%Tene
   endif
   md%Uene  = energy%E_tot
   md%E_tot = energy%E_tot + md%Tene  !??
   Htot     = energy%E_tot + Enh      !for NHC
   !how about field energy??


   !---Write section---
   ! Export to standard log file
   if(comm_is_root(nproc_id_global)) then
      write(*,120) it*dt, md%Tene, md%Uene, md%E_work, md%E_tot, md%E_tot-E_tot0,&
                   Enh, Htot, md%Temperature, sum1
120   format(1x,f10.4, 7e20.10E3,f12.3,e15.8)
   endif


!   ! Export to file_rt_data
!   call write_md_gs_data(it,Vion_gs,Tion,Temperature_ion,Eall,Eall00,Ework_integ_fdR_gs,Enh,Htot)


   ! Export to file_trj
   if (yn_out_rvf_rt=='y' .and. mod(it,out_rvf_rt_step)==0)then
      write(comment_line,110) it, it*dt
110   format("#md-gs  step=",i8,"   time",e16.6)
!      if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
!           write(comment_line,112) trim(comment_line), xi_nh
!112   format(a,"  xi_nh=",e18.10)
      call write_xyz(comment_line,"add","rvf",system)
   endif

   ! Export electronic density (cube or vtk)
   if(yn_out_dns_rt=='y' .and. mod(it,out_dns_rt_step)==0) then
      call write_dns(lg,mg,ng,srho%f,matbox_m,matbox_m2,system%hgs,iscfrt,srho%f,it)
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
call write_band_information(system,energy)
call write_eigen(file_eigen,system,energy)
call write_info_data(it,system,energy,pp)

! write GS: analysis option
if(yn_out_psi =='y') call write_psi(lg,mg,system,info,spsi)
if(yn_out_dns =='y') call write_dns(lg,mg,ng,srho%f,matbox_m,matbox_m2,system%hgs,iscfrt)
if(yn_out_dos =='y') call write_dos(system,energy)
if(yn_out_pdos=='y') call write_pdos(lg,mg,system,info,pp,energy,spsi)
if(yn_out_elf =='y') call write_elf(iscfrt,0,lg,mg,ng,system,info,stencil,srho,srg,srg_ng,spsi)

call timer_end(LOG_WRITE_GS_RESULTS)

! write GS: binary data for restart
call timer_begin(LOG_WRITE_GS_DATA)
call write_bin(ofile%dir_out_restart,lg,mg,ng,system,info,spsi,it,mixing=mixing)
call timer_end(LOG_WRITE_GS_DATA)

!call timer_begin(LOG_WRITE_GS_INFO)  !if needed, please take back, sory: AY
!call timer_end(LOG_WRITE_GS_INFO)

call finalize_xc(xc_func)

call timer_end(LOG_TOTAL)


end subroutine main_dft_md
