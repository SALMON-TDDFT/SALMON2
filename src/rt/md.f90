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
module md_sub
  implicit none

contains

subroutine init_md(system,md)
  use structures, only: s_system,s_md
  use salmon_global, only: MI,out_rvf_rt, ensemble, thermostat,set_ini_velocity,step_velocity_scaling
  use salmon_communication, only: comm_is_root
  use salmon_parallel, only: nproc_id_global
  implicit none    
  type(s_system) :: system
  type(s_md) :: md

  if(out_rvf_rt=='n') then
     if (comm_is_root(nproc_id_global)) &
          write(*,*)" out_rvf_rt --> y : changed for md option"
     out_rvf_rt='y'
  endif

  allocate(md%Rion_last(3,MI))
  allocate(md%force_last(3,MI))

  md%E_work = 0d0

  if(ensemble=="NVT" .and. thermostat=="nose-hoover") md%xi_nh=0d0
  if(set_ini_velocity=='y' .or. step_velocity_scaling>=1) &
       call set_initial_velocity(system,md)
  if(set_ini_velocity=='r') call read_initial_velocity(system,md)

  !if(use_ms_maxwell == 'y' .and. use_potential_model=='n') then
  !   if(nproc_size_global.lt.nmacro) then
  !      write(*,*) "Error: "
  !      write(*,*) "  number of parallelization nodes must be equal to or larger than  number of macro grids"
  !      write(*,*) "  in md option with multi-scale"
  !      call end_parallel
  !      stop
  !   endif
  !endif

end subroutine init_md

subroutine set_initial_velocity(system,md)
  use structures, only: s_system,s_md
  use salmon_global, only: MI,Kion,temperature0_ion
  use salmon_parallel, only: nproc_id_global,nproc_group_global
  use salmon_communication, only: comm_is_root,comm_bcast
  use math_constants, only: Pi
  use const, only: umass,hartree2J,kB
  implicit none
  type(s_system) :: system
  type(s_md) :: md
  integer :: ia,ixyz,iseed
  real(8) :: rnd1,rnd2,rnd, sqrt_kT_im, kB_au, mass_au
  real(8) :: Temperature_ion, scale_v, Tion
  
  if (comm_is_root(nproc_id_global)) then
     write(*,*) "  Initial velocities with maxwell-boltzmann distribution was set"
     write(*,*) "  Set temperature is ", real(temperature0_ion)
  endif

  kB_au = kB/hartree2J  ![au/K]

  iseed= 123
  do ia=1,MI
     mass_au = umass * system%Mass(Kion(ia))
     sqrt_kT_im = sqrt( kB_au * temperature0_ion / mass_au )

     do ixyz=1,3
        call quickrnd(iseed,rnd1)
        call quickrnd(iseed,rnd2)
        rnd = sqrt(-2d0*log(rnd1))*cos(2d0*Pi*rnd2)
        system%Velocity(ixyz,ia) = rnd * sqrt_kT_im
     enddo
  enddo
  
  !!(check temperature)
  !Tion=0d0
  !do ia=1,MI
  !   Tion = Tion + 0.5d0*umass*system%Mass(Kion(ia))*sum(system%Velocity(:,ia)**2d0)
  !enddo
  !Temperature_ion = Tion * 2d0 / (3d0*MI) / kB_au
  !write(*,*)"  Temperature: random-vel",real(Temperature_ion)
  
 
  !center of mass of system is removed
  call remove_system_momentum(1,system)
  
  !scaling: set temperature exactly to input value
  Tion=0d0
  do ia=1,MI
     Tion = Tion + 0.5d0 * umass*system%Mass(Kion(ia)) * sum(system%Velocity(:,ia)**2d0)
  enddo
  Temperature_ion = Tion * 2d0 / (3d0*MI) / kB_au
  !write(*,*)"    Temperature: befor-scaling",real(Temperature_ion)
 
  scale_v = sqrt(temperature0_ion/Temperature_ion)
  if(temperature0_ion==0d0) scale_v=0d0
  system%Velocity(:,:) = system%Velocity(:,:) * scale_v
 
  !(check)
  Tion=0d0
  do ia=1,MI
     Tion = Tion + 0.5d0 * umass*system%Mass(Kion(ia)) * sum(system%Velocity(:,ia)**2d0)
  enddo
  Temperature_ion = Tion * 2d0 / (3d0*MI) / kB_au
  if (comm_is_root(nproc_id_global)) &
       write(*,*)"    Initial Temperature: after-scaling",real(Temperature_ion)

  md%Tene = Tion
  md%Temperature = Temperature_ion
  
  call comm_bcast(system%Velocity ,nproc_group_global)
 
end subroutine set_initial_velocity
   
subroutine read_initial_velocity(system,md)
  ! initial velocity for md option can be given by external file 
  ! specified by file_ini_velocity option
  ! format is :
  ! do i=1,MI
  !    vx(i)  vy(i)  vz(i)
  ! enddo
  ! xi_nh  !only for nose-hoover thermostat option
  use structures, only: s_system,s_md
  use salmon_global, only: MI,file_ini_velocity, ensemble, thermostat
  use salmon_parallel, only: nproc_id_global,nproc_group_global,end_parallel
  use salmon_communication, only: comm_is_root,comm_bcast
  implicit none
  type(s_system) :: system
  type(s_md) :: md
  integer :: ia,ixyz

  if(comm_is_root(nproc_id_global)) then
     write(*,*) "Read initial velocity for MD"
     write(*,*) "file_ini_velocity=", trim(file_ini_velocity)
     if(file_ini_velocity(1:4)=="none") then
        write(*,*) "set file name in file_ini_velocity keyword"
        call end_parallel
        stop
     endif

     open(411, file=file_ini_velocity, status="old")
     do ia=1,MI
        read(411,*) (system%Velocity(ixyz,ia),ixyz=1,3)
     enddo
     if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
        read(411,*,err=100,end=100) md%xi_nh !if no value, skip reading (xi_nh=0)
     endif
100  close(411)
  endif
  call comm_bcast(system%Velocity,nproc_group_global)

end subroutine read_initial_velocity

subroutine remove_system_momentum(flag_print_check,system)
  ! remove center of mass and momentum of whole system
  use structures, only: s_system
  use salmon_global, only: MI,Kion
  use salmon_communication, only: comm_is_root
  use salmon_parallel, only: nproc_id_global
  use const, only: umass
!  use scf_data, only: ppg,ppg_all
  implicit none
  type(s_system) :: system
  integer :: ia, flag_print_check
  real(8) :: v_com(3), sum_mass, mass_au

  !velocity of center of mass is removed
  v_com(:)=0d0
  sum_mass=0d0
  do ia=1,MI
     mass_au = umass * system%Mass(Kion(ia))
     v_com(:) = v_com(:) + mass_au * system%Velocity(:,ia)
     sum_mass = sum_mass + mass_au
  enddo
  v_com(:) = v_com(:)/sum_mass
  do ia=1,MI
     system%Velocity(:,ia) = system%Velocity(:,ia) - v_com(:)
  enddo

  !rotation of system is removed
  !(this is only for isolated system)--- do nothing

  !(check velocity of center of mass)
  if(flag_print_check==1) then
     v_com(:)=0d0
     do ia=1,MI
        v_com(:) = v_com(:) + umass*system%Mass(Kion(ia)) * system%Velocity(:,ia)
     enddo
     v_com(:) = v_com(:) / sum_mass
     if(comm_is_root(nproc_id_global)) write(*,*)"    v_com =",real(v_com(:))
  endif

end subroutine remove_system_momentum

subroutine cal_Tion_Temperature_ion(Ene_ion,Temp_ion,system)
  use structures, only: s_system
  use salmon_global, only: MI,Kion
  use const, only: umass,hartree2J,kB
  implicit none
  type(s_system) :: system
  integer :: ia
  real(8) :: mass_au, Ene_ion,Temp_ion

  Ene_ion = 0.d0
  do ia=1,MI
     mass_au = umass * system%Mass(Kion(ia))
     Ene_ion = Ene_ion + 0.5d0 * mass_au * sum(system%Velocity(:,ia)**2d0)
  enddo
  Temp_ion = Ene_ion * 2d0 / (3d0*MI) / (kB/hartree2J)

  return
end subroutine cal_Tion_Temperature_ion


subroutine time_evolution_step_md_part1(system,force,md)
  use structures, only: s_system, s_force, s_md
  use salmon_global, only: MI,Kion,dt, Rion
  use const, only: umass,hartree2J,kB
  implicit none
  type(s_system) :: system
  type(s_force),intent(in) :: force
  type(s_md) :: md
  integer :: iatom
  real(8) :: mass_au, dt_h

     dt_h = dt*0.5d0

     !update ion velocity with dt/2
     do iatom=1,MI
        mass_au = umass * system%Mass(Kion(iatom))
        system%Velocity(:,iatom) = system%Velocity(:,iatom) + force%F(:,iatom)/mass_au * dt_h
     enddo

     md%Rion_last(:,:) = system%Rion(:,:)
     md%force_last(:,:)= force%F(:,:)

     !update ion coordinate with dt
     do iatom=1,MI
        system%Rion(:,iatom) = system%Rion(:,iatom) + system%Velocity(:,iatom) *dt
     enddo
     Rion(:,:) = system%Rion(:,:) !copy (old variable, Rion, is still used in somewhere)

end subroutine 

subroutine update_pseudo_rt(itt,system,stencil,lg,ng,fg,ppg,ppg_all,ppn)
  use structures, only: s_system,s_stencil,s_rgrid,s_pp_nlcc,s_pp_grid,s_fourier_grid
  use salmon_global, only: iperiodic,step_update_ps,step_update_ps2
  use const, only: umass,hartree2J,kB
  implicit none
  type(s_system) :: system
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: ng
  type(s_fourier_grid) :: fg
  type(s_stencil),intent(inout) :: stencil
  type(s_pp_nlcc), intent(in) :: ppn
  type(s_pp_grid) :: ppg,ppg_all
  integer :: itt

  if(iperiodic==3) call get_fourier_grid_G_rt(system,lg,ng,fg)

  !update pseudopotential
  if (mod(itt,step_update_ps)==0 ) then
     !xxxx call prep_ps_periodic('update_all       ')
     call dealloc_init_ps(ppg,ppg_all,ppn)
     call init_ps(system%al,system%brl,stencil%matrix_A)
  else if (mod(itt,step_update_ps2)==0 ) then
     !xxxx call prep_ps_periodic('update_wo_realloc')
     !xxxxxxx this option is not yet made xxxxxx
     call dealloc_init_ps(ppg,ppg_all,ppn)
     call init_ps(system%al,system%brl,stencil%matrix_A)
  endif

end subroutine 

subroutine time_evolution_step_md_part2(system,force,md)
  use structures, only: s_system, s_force, s_md
  use salmon_global, only: MI,Kion,dt,stop_system_momt
  use const, only: umass,hartree2J,kB
  implicit none
  type(s_system) :: system
  type(s_force),intent(in) :: force
  type(s_md) :: md
  integer :: iatom
  real(8) :: mass_au,dt_h, aforce(3,MI), dR(3,MI)

  dt_h = dt*0.5d0
  aforce(:,:) = 0.5d0*( md%force_last(:,:) + force%F(:,:) )

  !update ion velocity with dt/2
  dR(:,:) = system%Rion(:,:) - md%Rion_last(:,:)
  do iatom=1,MI
     mass_au = umass * system%Mass(Kion(iatom))
     system%Velocity(:,iatom) = system%Velocity(:,iatom) + force%F(:,iatom)/mass_au * dt_h
     md%E_work = md%E_work - sum(aforce(:,iatom)*dR(:,iatom))
  enddo


  if (stop_system_momt=='y') call remove_system_momentum(0,system)
  call cal_Tion_Temperature_ion(md%Tene,md%Temperature,system)

end subroutine 

end module md_sub
