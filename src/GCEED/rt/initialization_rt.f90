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
!===================================================================================================================================
module initialization_rt_sub

  implicit none

contains

subroutine initialization_rt( Mit, system, energy, rt, md, singlescale,  &
                     stencil, fg, poisson,  &
                     lg, mg, ng,  &
                     info, info_field,  &
                     xc_func, dmat, ofl, j_e,  &
                     srg, srg_ng,  &
                     spsi_in, spsi_out, tpsi, srho, srho_s,  &
                     V_local, Vbox, sVh, sVh_stock1, sVh_stock2, sVxc, sVpsl,&
                     pp, ppg, ppn )
use inputoutput
use scf_data
use allocate_mat_sub
use deallocate_mat_sub

  use math_constants, only: pi, zi
  use structures
  use salmon_parallel, only: nproc_id_global, nproc_group_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast, comm_sync_all
  use salmon_xc
  use timer
  use write_sub, only: write_xyz,write_rt_data_0d,write_rt_data_3d,write_rt_energy_data, &
                       write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d
  use code_optimization
  use initialization_sub
  use input_pp_sub
  use prep_pp_sub
  use density_matrix, only: calc_density
  use writefield
  use salmon_pp, only: calc_nlcc
  use force_sub, only: calc_force_salmon
  use hamiltonian
  use md_sub, only: init_md
  use fdtd_coulomb_gauge, only: ls_singlescale
  use taylor_sub, only: taylor_coe
  use checkpoint_restart_sub
  use hartree_sub, only: hartree
  use salmon_Total_Energy
  implicit none

  real(8)       :: debye2au   ! [D]  -> [a.u.] 
  integer       :: Ntime

  type(s_rgrid) :: lg
  type(s_rgrid) :: mg
  type(s_rgrid) :: ng
  type(s_dft_system)  :: system
  type(s_rt) :: rt
  type(s_process_info) :: pinfo
  type(s_orbital_parallel) :: info
  type(s_field_parallel) :: info_field
  type(s_poisson) :: poisson
  type(s_stencil) :: stencil
  type(s_xc_functional) :: xc_func
  type(s_reciprocal_grid) :: fg
  type(s_dft_energy) :: energy
  type(s_md) :: md
  type(s_ofile) :: ofl
  type(s_scalar) :: sVpsl
  type(s_scalar) :: srho,sVh,sVh_stock1,sVh_stock2,Vbox
  type(s_scalar),allocatable :: srho_s(:),V_local(:),sVxc(:)
  type(s_dmatrix) :: dmat
  type(s_orbital) :: spsi_in,spsi_out
  type(s_orbital) :: tpsi ! temporary wavefunctions
  type(s_sendrecv_grid) :: srg,srg_ng
  type(s_pp_info) :: pp
  type(s_pp_grid) :: ppg
  type(s_pp_nlcc) :: ppn
  type(s_vector)  :: j_e ! microscopic electron number current density
  type(ls_singlescale) :: singlescale
  type(s_ofile) :: ofile
  
  integer :: iob, i1,i2,i3, iik,jspin, Mit
  integer :: idensity, idiffDensity
  integer :: jj, ix,iy,iz
  real(8) :: rbox_array(10), rbox_array2(10)
  real(8),allocatable :: R1(:,:,:)
  character(10) :: fileLaser
  character(100):: comment_line
  real(8) :: curr_e_tmp(3,2), curr_i_tmp(3)
  
  call timer_begin(LOG_INIT_RT)

  call init_xc(xc_func, ispin, cval, xcname=xc, xname=xname, cname=cname)
  
  iSCFRT=2
  
  iwdenstep=30 
  
  call convert_input_rt(Ntime)
  
  call set_filename
  
  if(comm_is_root(nproc_id_global))then
    write(*,*)
    write(*,*) "Total time step      =",Ntime
    write(*,*) "Time step[fs]        =",dt*au_time_fs
    write(*,*) "Energy range         =",Nenergy
    write(*,*) "Energy resolution[eV]=",dE*au_energy_ev
    write(*,*) "ikind_eext is         ", ikind_eext
    write(*,*) "Step for writing dens=", iwdenstep
    select case (ikind_eext)
      case(0)
        write(*,*) "Field strength[a.u.] =",e_impulse
      case(1,6,7,8,15)
        write(*,20) "Laser frequency     =",omega1*au_energy_ev,"[eV]"
        write(*,21) "Pulse width of laser=",tw1*au_time_fs, "[fs]"
        write(*,22) "Laser intensity     =",I_wcm2_1,       "[W/cm^2]"
      case(4,12)
        write(*,23) "Laser frequency     =",omega1*au_energy_ev,omega2*au_energy_ev,"[eV]"
        write(*,24) "Pulse width of laser=",tw1*au_time_fs,tw2*au_time_fs, "[fs]"
        write(*,25) "Laser intensity     =",I_wcm2_1,I_wcm2_2,       "[W/cm^2]"
        write(*,21) "delay time          =",t1_t2*au_time_fs,     "[fs]"
    end select
  20 format(a21,f5.2, a4)
  21 format(a21,f16.8,a4)
  22 format(a21,e16.8,a8)
  23 format(a21,2f5.2, a4)
  24 format(a21,2f16.8,a4)
  25 format(a21,2e16.8,a8)
  
  end if
  
  debye2au = 0.393428d0
  
  select case (ikind_eext)
    case(1)
      if(I_wcm2_1>=1.d-12)then
        E_amplitude1=sqrt(I_wcm2_1)*1.0d2*2.74492d1/(5.14223d11)!I[W/cm^2]->E[a.u.]
      end if
      if(I_wcm2_2>=1.d-12)then
        E_amplitude2=sqrt(I_wcm2_2)*1.0d2*2.74492d1/(5.14223d11)!I[W/cm^2]->E[a.u.]
      else
        if(abs(E_amplitude2)<=1.d-12)then
          E_amplitude2=0.d0
        end if
      end if
  end select
  
  call timer_end(LOG_INIT_RT)
  
  call timer_begin(LOG_READ_GS_DATA)
  
  ! +----------------+
  ! | initialization |
  ! +----------------+
  
  call init_dft(iSCFRT,nproc_group_global,pinfo,info,info_field,lg,mg,ng,system,stencil,fg,poisson,srg,srg_ng,ofile)
  
  call init_code_optimization
  call allocate_mat(ng,mg,lg) ! future work: remove this line
  
  call allocate_scalar(mg,srho)
  call allocate_scalar(mg,sVh)
  call allocate_scalar(mg,sVh_stock1)
  call allocate_scalar(mg,sVh_stock2)
  call allocate_scalar_with_shadow(lg,Nd,Vbox)
  call allocate_scalar(mg,sVpsl)
  allocate(srho_s(system%nspin),V_local(system%nspin),sVxc(system%nspin))
  do jspin=1,system%nspin
    call allocate_scalar(mg,srho_s(jspin))
    call allocate_scalar(mg,V_local(jspin))
    call allocate_scalar(mg,sVxc(jspin))
  end do
  allocate(ppg%Vpsl_atom(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),natom))
  call read_pslfile(system,pp,ppg)
  call init_ps(lg,mg,ng,system,info,info_field,fg,poisson,pp,ppg,sVpsl)
  
  call allocate_orbital_complex(system%nspin,mg,info,spsi_in)
  call allocate_orbital_complex(system%nspin,mg,info,spsi_out)
  call allocate_orbital_complex(system%nspin,mg,info,tpsi)
  call allocate_dmatrix(system%nspin,mg,info,dmat)
  
  call timer_begin(LOG_RESTART_SYNC)
  call timer_begin(LOG_RESTART_SELF)
  call restart_rt(lg,mg,ng,system,info,spsi_in,Mit,sVh_stock1=sVh_stock1,sVh_stock2=sVh_stock2)
  call timer_end(LOG_RESTART_SELF)
  call comm_sync_all
  call timer_end(LOG_RESTART_SYNC)
  if(yn_restart=='n') Mit=0
  
  call calc_nlcc(pp, system, mg, ppn)
  if (comm_is_root(nproc_id_global)) then
    write(*, '(1x, a, es23.15e3)') "Maximal rho_NLCC=", maxval(ppn%rho_nlcc)
    write(*, '(1x, a, es23.15e3)') "Maximal tau_NLCC=", maxval(ppn%tau_nlcc)
  end if
  
  call calc_density(system,srho_s,spsi_in,info,mg)
  srho%f = 0d0
  do jspin=1,system%nspin
     srho%f = srho%f + srho_s(jspin)%f
  end do
  if(yn_restart=='y')then
    sVh%f = 2.d0*sVh_stock1%f - sVh_stock2%f
    sVh_stock2%f = sVh_stock1%f
  end if
  call hartree(lg,mg,ng,info_field,system,poisson,srg_ng,stencil,srho,sVh,fg)
  call exchange_correlation(system,xc_func,ng,srg_ng,srho_s,ppn,info_field%icomm_all,sVxc,energy%E_xc)
  call allgatherv_vlocal(ng,mg,info_field,system%nspin,sVh,sVpsl,sVxc,V_local)
  if(yn_restart=='y')then
    sVh_stock1%f=sVh%f
  else if(yn_restart=='n')then
    sVh_stock1%f=sVh%f
    sVh_stock2%f=sVh%f
  end if
  
  allocate(energy%esp(system%no,system%nk,system%nspin))
  
  call timer_end(LOG_READ_GS_DATA)
  
  ! calculation of GS total energy
  call calc_eigen_energy(energy,spsi_in,spsi_out,tpsi,system,info,mg,V_local,stencil,srg,ppg)
  select case(iperiodic)
  case(0)
     call calc_Total_Energy_isolated(energy,system,info,ng,pp,srho_s,sVh,sVxc)
  case(3)
     call calc_Total_Energy_periodic(energy,system,pp,fg,.true.)
  end select
  energy%E_tot0 = energy%E_tot
  
  ! +-------------+
  ! | old fashion |
  ! +-------------+
  
  call timer_begin(LOG_READ_RT_DATA)
  
  allocate( rt%rIe(     0:Ntime) )
  allocate( rt%dDp_e( 3,0:Ntime) )
  allocate( rt%Dp_e(  3,0:Ntime) )
  allocate( rt%Dp_i(  3,0:Ntime) )
  allocate( rt%Qp_e(3,3,0:Ntime) )
  call init_A(Ntime,Mit,rt)
  itotNtime = Ntime
  if (yn_restart /= 'y') Mit=0
  call timer_end(LOG_READ_RT_DATA)
  
  
  call timer_begin(LOG_INIT_RT)
  
  ntmg=1
  ! 'Hartree' parameter
  
  poisson%iterVh = 0        ! Iteration counter
  
  
  if(comm_is_root(nproc_id_global))then
    write(*, *) 
    write(*, *) "dip2boundary", dip2boundary(1), dip2boundary(2)
    write(*, *) "dip2center", dip2center(1), dip2center(2)
    write(*, *) "dip2boundary[A]", dip2boundary(1)*au_length_aa, dip2boundary(2)*au_length_aa
    write(*, *) "dip2center[A]", dip2center(1)*au_length_aa, dip2center(2)*au_length_aa
    write(*, *) 
  end if
  call timer_end(LOG_INIT_RT)
  
  !Open output files and print header parts (Please move and put this kinds here!)
  !(output file names)
  if(comm_is_root(nproc_id_global))then
     write(ofl%file_rt_data,"(2A,'_rt.data')") trim(base_directory),trim(SYSname)
     write(ofl%file_rt_energy_data,"(2A,'_rt_energy.data')") trim(base_directory),trim(SYSname)
     write(ofl%file_response_data,"(2A,'_response.data')") trim(base_directory),trim(SYSname)
     write(ofl%file_pulse_data,"(2A,'_pulse.data')") trim(base_directory),trim(SYSname)
  endif
  call comm_bcast(ofl%file_rt_data,       nproc_group_global)
  call comm_bcast(ofl%file_rt_energy_data,nproc_group_global)
  call comm_bcast(ofl%file_response_data, nproc_group_global)
  call comm_bcast(ofl%file_pulse_data,    nproc_group_global)
  
  !(write header)
  if(comm_is_root(nproc_id_global))then
    !(header of standard output)
    select case(iperiodic)
    case(0)
      write(*,'(1x,a10,a11,a48,a15,a18,a10)') &
                  "time-step ", "time[fs]",   &
                  "Dipole moment(xyz)[A]"     &
                 ,"electrons", "Total energy[eV]", "iterVh"
    case(3)
      write(*,'(1x,a10,a11,a48,a15,a18)')   &
                  "time-step", "time[fs] ", &
                  "Current(xyz)[a.u.]",     &
                  "electrons", "Total energy[eV] "
    end select
    write(*,'("#",7("----------"))')
  
    !(header of SYSname_rt.data)
    select case(iperiodic)
    case(0) ; call write_rt_data_0d(-1,ofl,dt,system,rt)
    case(3) ; call write_rt_data_3d(-1,ofl,dt,system,curr_e_tmp,curr_i_tmp)
    end select
  
    !(header of SYSname_rt_energy.data)
    call write_rt_energy_data(-1,ofl,dt,energy,md)
  
    !(header in SYSname_proj.data)
    if(iwrite_projection==1)then
      open(41,file=file_Projection)
      write(41,'("#",5X,"time[fs]",4("    projection"))')
      write(41,'("#",13x,4("  orbital",i5))') (iwrite_projection_ob(jj),jj=1,4)
      write(41,'("#",13x,4("        k",i5))') (iwrite_projection_k(jj), jj=1,4)
      write(41,'("#",7("----------"))')
    end if
  end if
  
  if(yn_md=='y' .or. yn_out_rvf_rt=='y')then
     call write_xyz(comment_line,"new","rvf",system)
  endif
  
  !---------------------------- time-evolution
  
  call timer_begin(LOG_INIT_TIME_PROPAGATION)
  
  idensity=0
  idiffDensity=1
  fileLaser= "laser.out"
  
  allocate( R1(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2), &
                                 lg%is(3):lg%ie(3)))
  
  ! External Field Direction
  
  select case(mod(lg%num(1),2))
    case(1)
      do i1=lg%is(1),lg%ie(1)
         vecR(1,i1,:,:)=dble(i1)
      end do
    case(0)
      do i1=lg%is(1),lg%ie(1)
         vecR(1,i1,:,:)=dble(i1)-0.5d0
      end do
  end select
  
  select case(mod(lg%num(2),2))
    case(1)
      do i2=lg%is(2),lg%ie(2)
         vecR(2,:,i2,:)=dble(i2)
      end do
    case(0)
      do i2=lg%is(2),lg%ie(2)
         vecR(2,:,i2,:)=dble(i2)-0.5d0
      end do
  end select
  
  select case(mod(lg%num(3),2))
    case(1)
      do i3=lg%is(3),lg%ie(3)
         vecR(3,:,:,i3)=dble(i3)
      end do
    case(0)
      do i3=lg%is(3),lg%ie(3)
         vecR(3,:,:,i3)=dble(i3)-0.5d0
      end do
  end select
  
  !$OMP parallel do collapse(2) private(iz,iy,ix)
  do iz=lg%is(3),lg%ie(3)
  do iy=lg%is(2),lg%ie(2)
  do ix=lg%is(1),lg%ie(1)
     R1(ix,iy,iz)=(epdir_re1(1)*lg%coordinate(ix,1)+   &
                   epdir_re1(2)*lg%coordinate(iy,2)+   &
                   epdir_re1(3)*lg%coordinate(iz,3))
  end do
  end do
  end do
  
  if(num_dipole_source>=1)then
    allocate(vonf_sd(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(eonf_sd(3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    vonf_sd=0.d0
    eonf_sd=0.d0
    call set_vonf_sd(mg,lg,system%Hgs)
  end if
  
  if(yn_restart /= 'y')then
    rbox_array=0d0
    do i1=1,3
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
         rbox_array(i1) = rbox_array(i1) + vecR(i1,ix,iy,iz) * srho%f(ix,iy,iz)
      end do
      end do
      end do
    end do
  
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
       rbox_array(4) = rbox_array(4) + srho%f(ix,iy,iz)
    end do
    end do
    end do
  
    call comm_summation(rbox_array,rbox_array2,4,nproc_group_global)
    rt%Dp0_e(1:3) = -rbox_array2(1:3) * system%Hgs(1:3) * system%Hvol
  
  end if
  if(comm_is_root(nproc_id_global))then
    write(*,'(a30)', advance="no") "Static dipole moment(xyz) ="
    write(*,'(3e20.10)') (rt%Dp0_e(i1)*au_length_aa, i1=1,3)
    write(*,*)
  endif
  
  ! Initial wave function
  if(iperiodic==0 .and. ikind_eext==0 .and. yn_restart /= 'y')then
    do iik=info%ik_s,info%ik_e
    do iob=info%io_s,info%io_e
    do jspin=1,system%nspin
          spsi_in%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),jspin,iob,iik,1) &
          = exp(zi*e_impulse*R1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),  &
             mg%is(3):mg%ie(3)))   &
             *  spsi_in%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),jspin,iob,iik,1)
    end do
    end do
    end do
  end if
  
  rt%rIe(0)     = rbox_array2(4) * system%Hvol
  rt%dDp_e(:,0) = 0d0
  rt%Dp_e(:,0)  = 0d0
  rt%Dp_i(:,0)  = 0d0
  rt%Qp_e(:,:,0)= 0d0
  
    do itt=0,0
      if(yn_out_dns_rt=='y')then
        call write_dns(lg,mg,ng,srho%f,matbox_m,matbox_m2,system%hgs,iscfrt,srho%f,itt)
      end if
      if(yn_out_elf_rt=='y')then
        call write_elf(iscfrt,itt,lg,mg,ng,system,info,stencil,srho,srg,srg_ng,spsi_in)
      end if
      if(yn_out_estatic_rt=='y')then
        call write_estatic(lg,ng,system%hgs,stencil,info_field,sVh,srg_ng,itt)
      end if
    end do
  
  if(iperiodic==3) call calc_Aext(Mit,rt)
  
  if(yn_md=='y') call init_md(system,md)
  
  ! single-scale Maxwell-TDDFT
  if(use_singlescale=='y') then
    if(comm_is_root(nproc_id_global)) write(*,*) "single-scale Maxwell-TDDFT method"
    call allocate_vector(mg,j_e)
    call allocate_scalar(mg,system%div_Ac)
    call allocate_vector(mg,system%Ac_micro)
  end if
  
  !-------------------------------------------------- Time evolution
  
  !(force at initial step)
  if(yn_md=='y' .or. yn_out_rvf_rt=='y')then
     call calc_force_salmon(system,pp,fg,info,mg,stencil,srg,ppg,spsi_in)
  
     !open trj file for coordinate, velocity, and force (rvf) in xyz format
     write(comment_line,10) -1, 0.0d0
     if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
          &  write(comment_line,12) trim(comment_line), md%xi_nh
     call write_xyz(comment_line,"add","rvf",system)
  10 format("#rt   step=",i8,"   time",e16.6)
  12 format(a,"  xi_nh=",e18.10)
  end if
  
  call timer_end(LOG_INIT_TIME_PROPAGATION)
  
  
  call timer_begin(LOG_INIT_RT)
  allocate(zc(N_hamil))
  call taylor_coe(N_hamil,dt,zc)

  deallocate (R1)

  call timer_end(LOG_INIT_RT)

contains

subroutine init_code_optimization
  implicit none
  integer :: ignum(3)

  call switch_stencil_optimization(mg%num)
  call switch_openmp_parallelization(mg%num)

  if(iperiodic==3 .and. product(pinfo%npdomain_orbital)==1) then
     ignum = mg%num
  else
     ignum = mg%num + (nd*2)
  end if
  call set_modulo_tables(ignum)

  if (comm_is_root(nproc_id_global)) then
     call optimization_log(pinfo)
  end if
end subroutine init_code_optimization

end subroutine initialization_rt

end module initialization_rt_sub
