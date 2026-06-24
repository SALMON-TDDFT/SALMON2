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
!===================================================================================================================================
module initialization_rt_sub
  use nvtx_wrapper
  implicit none

contains

subroutine initialization_rt( Mit, system, energy, ewald, rt, md, &
                     singlescale,  &
                     stencil, fg, poisson,  &
                     lg, mg, info,  &
                     xc_func, ofl,  &
                     srg, srg_scalar,  &
                     spsi_in, spsi_out, tpsi, rho, rho_jm, rho_s,  &
                     V_local, Vbox, Vh, Vh_stock1, Vh_stock2, Vxc, Vpsl,&
                     pp, ppg, ppn  )
  use inputoutput
  use math_constants, only: pi, zi
  use structures
  use parallelization, only: nproc_id_global, nproc_group_global, end_parallel
  use communication, only: comm_is_root, comm_summation, comm_bcast, comm_sync_all
  use salmon_xc
  use timer
  use write_sub, only: write_xyz,write_rt_data_0d,write_rt_data_3d,write_rt_energy_data, &
                       write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d,write_stress_rt,&
                       init_projection,write_rt_spin,write_current_decomposed
  use code_optimization
  use initialization_sub
  use prep_pp_sub
  use density_matrix, only: calc_density,calc_microscopic_current
  use writefield
  use salmon_pp, only: calc_nlcc, read_pslfile
  use force_sub, only: calc_force
  use hamiltonian
  use md_sub, only: init_md
  use fdtd_coulomb_gauge, only: init_singlescale
  use checkpoint_restart_sub
  use hartree_sub, only: hartree
  use Total_Energy
  use em_field, only: set_vonf,calc_Ac_ext_t
  use dip, only: calc_dip
  use sendrecv_grid
  use salmon_global, only: quiet, theory, yn_out_stress
  use gram_schmidt_orth, only: gram_schmidt
  use jellium, only: make_rho_jm
  use filesystem, only: open_filehandle
  use lcfo, only: init_conventional_from_dcdft
  use stress_sub, only: calc_stress, refresh_stress_output_state, s_stress_field_state
  use nstate_active_util, only: calc_nstate_active
  implicit none
  integer,parameter :: Nd = 4

  real(8)       :: debye2au   ! [D]  -> [a.u.]

  type(s_rgrid) :: lg
  type(s_rgrid) :: mg
  type(s_dft_system)  :: system
  type(s_rt) :: rt
  type(s_parallel_info) :: info
  type(s_poisson) :: poisson
  type(s_stencil) :: stencil
  type(s_xc_functional) :: xc_func
  type(s_reciprocal_grid) :: fg
  type(s_dft_energy) :: energy
  type(s_ewald_ion_ion) :: ewald
  type(s_md) :: md
  type(s_ofile) :: ofl
  type(s_scalar) :: Vpsl
  type(s_scalar) :: rho,rho_jm,Vh,Vh_stock1,Vh_stock2,Vbox
  type(s_scalar),allocatable :: rho_s(:),V_local(:),Vxc(:)
  type(s_orbital) :: spsi_in,spsi_out
  type(s_orbital) :: tpsi ! temporary wavefunctions
  type(s_sendrecv_grid) :: srg,srg_scalar
  type(s_pp_info) :: pp
  type(s_pp_grid) :: ppg
  type(s_pp_nlcc) :: ppn
  type(s_singlescale) :: singlescale
  type(s_ofile) :: ofile
  integer :: iob, i1,iik,jspin, Mit, m, n
  integer :: idensity, idiffDensity
  integer :: jj, ix,iy,iz
  real(8) :: rbox_array2(4) !,tt
  real(8),allocatable :: R1(:,:,:)
  character(10) :: fileLaser
  character(100):: comment_line
  real(8) :: curr_e_tmp(3,2), curr_i_tmp(3)
  integer :: itt
  integer :: Mit_restart_loaded, stress_label_iter, stress_state_iter
  logical :: rion_update
  type(s_stress_field_state) :: field_state
  integer :: nact

  call nvtxStartRange('initialization_rt', __LINE__)

  call timer_begin(LOG_INIT_RT)

  call init_xc(xc_func, spin, cval, xcname=xc, xname=xname, cname=cname)
  
  if((.not. quiet) .and. comm_is_root(nproc_id_global))then
    write(*,*)
    write(*,*) "Total time step      =",Nt
    write(*,*) "Time step[fs]        =",dt*au_time_fs
    write(*,*) "Energy range         =",Nenergy
    write(*,*) "Energy resolution[eV]=",dE*au_energy_ev
    if(ae_shape1 == 'impulse')then 
      write(*,*) "Field strength[a.u.] =",e_impulse
    else
      if(I_wcm2_2<1.d-12.and.E_amplitude2<=1.d-12)then
        write(*,20) "Laser frequency     =",omega1*au_energy_ev,"[eV]"
        write(*,21) "Pulse width of laser=",tw1*au_time_fs, "[fs]"
        write(*,22) "Laser intensity     =",I_wcm2_1,       "[W/cm^2]"
      else
        write(*,23) "Laser frequency     =",omega1*au_energy_ev,omega2*au_energy_ev,"[eV]"
        write(*,24) "Pulse width of laser=",tw1*au_time_fs,tw2*au_time_fs, "[fs]"
        write(*,25) "Laser intensity     =",I_wcm2_1,I_wcm2_2,       "[W/cm^2]"
        write(*,21) "delay time          =",t1_t2*au_time_fs,     "[fs]"
      end if
    end if
  20 format(a21,f5.2, a4)
  21 format(a21,f16.8,a4)
  22 format(a21,e16.8,a8)
  23 format(a21,2f5.2, a4)
  24 format(a21,2f16.8,a4)
  25 format(a21,2e16.8,a8)
  endif
  
  
  debye2au = 0.393428d0

  if(ae_shape1 /= 'impulse')then 
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
  end if
  
  call timer_end(LOG_INIT_RT)
  
  call timer_begin(LOG_READ_RT_DATA)
  
  allocate( rt%rIe(     0:nt) )
  allocate( rt%dDp_e( 3,0:nt) )
  allocate( rt%Dp_e(  3,0:nt) )
  allocate( rt%Dp_i(  3,0:nt) )
  allocate( rt%Qp_e(3,3,0:nt) )
  
  allocate( rt%curr  (3,0:nt+1) )
  allocate( rt%E_ext (3,0:nt+1) )
  allocate( rt%E_ind (3,0:nt+1) )
  allocate( rt%E_tot (3,0:nt+1) )
  allocate( rt%Ac_ext(3,0:nt+1) )
  allocate( rt%Ac_ind(3,0:nt+1) )
  allocate( rt%Ac_tot(3,0:nt+1) )
  
  rt%curr  = 0d0
  rt%E_ext = 0d0
  rt%E_ind = 0d0
  rt%E_tot = 0d0
  rt%Ac_ext = 0d0
  rt%Ac_ind = 0d0
  rt%Ac_tot = 0d0
  
  call calc_Ac_ext_t(0d0, dt, 0, nt+1, rt%Ac_ext(:, 0:nt+1))

  do itt = 0, nt
    rt%E_ext(:, itt) = -(rt%Ac_ext(:, itt+1) - rt%Ac_ext(:, itt)) / dt
  end do
  rt%E_ext(:, nt+1) = rt%E_ext(:, nt) 
  
  rt%Ac_tot(1:3, 0:nt+1) = rt%Ac_ext(1:3, 0:nt+1)
  rt%E_tot(1:3, 0:nt+1) = rt%E_ext(1:3, 0:nt+1) 
  
  Mit_restart_loaded = Mit
  if (yn_restart == 'n') Mit=0
  call timer_end(LOG_READ_RT_DATA)
  
  call timer_begin(LOG_READ_GS_DATA)
  call nvtxStartRange('READ_GS_DATA', __LINE__)
  
  call init_dft(nproc_group_global,info,lg,mg,system,stencil,fg,poisson,srg,srg_scalar,ofile)
  
  call init_code_optimization
  
  call allocate_scalar(mg,rho)
  call allocate_scalar(mg,Vh)
  call allocate_scalar(mg,Vh_stock1)
  call allocate_scalar(mg,Vh_stock2)
  call allocate_scalar_with_shadow(lg,Nd,Vbox)
  call allocate_scalar(mg,Vpsl)
  allocate(rho_s(system%nspin),V_local(system%nspin),Vxc(system%nspin),rt%rho0_s(system%nspin))
  do jspin=1,system%nspin
    call allocate_scalar(mg,rho_s(jspin))
    call allocate_scalar(mg,V_local(jspin))
    call allocate_scalar(mg,Vxc(jspin))
    call allocate_scalar(mg,rt%rho0_s(jspin))
  end do
  if(yn_jm=='n') then
    call read_pslfile(system,pp)
    call init_ps(lg,mg,system,info,fg,poisson,pp,ppg,Vpsl)
  else
    !make positive back ground charge density for using jellium model
    call allocate_scalar(mg,rho_jm)
    call make_rho_jm(lg,mg,info,system,rho_jm)
  end if

  ! --- nstate_active: reduce the RT propagation orbital set to the active (occupied / above-threshold)
  !     orbitals before the wavefunction arrays, parallel orbital-distribution and srg are finalized.
  !     The full nstate is retained for projection (init_projection, below, is unaffected).
  call apply_nstate_active(system, info, srg, nact)

  call allocate_orbital_complex(system%nspin,mg,info,spsi_in)
  call allocate_orbital_complex(system%nspin,mg,info,spsi_out)
  call allocate_orbital_complex(system%nspin,mg,info,tpsi)
  
  if(propagator=='aetrs')then
    allocate(rt%vloc_t(system%nspin),rt%vloc_new(system%nspin),rt%vloc_old(system%nspin,2))
    do jspin=1,system%nspin
      call allocate_scalar(mg,rt%vloc_t(jspin))
      call allocate_scalar(mg,rt%vloc_new(jspin))
      call allocate_scalar(mg,rt%vloc_old(jspin,1))
      call allocate_scalar(mg,rt%vloc_old(jspin,2))
    end do
  end if

  !$acc enter data copyin(system)
  !$acc enter data copyin(info)
  !$acc enter data copyin(mg)
  !$acc enter data copyin(stencil)
  !$acc enter data copyin(V_local)
  !$acc enter data copyin(spsi_in,spsi_out,tpsi) 
  !$acc enter data copyin(ppg)
  
  call timer_begin(LOG_RESTART_SYNC)
  call timer_begin(LOG_RESTART_SELF)
  if(yn_conventional_from_dcdft=='n') then
    call restart_rt(lg,mg,system,info,spsi_in,Mit,rt,Vh_stock1=Vh_stock1,Vh_stock2=Vh_stock2)
  else
  ! conventional TDDFT but wavefunctions are reconstructed from DC-LCFO data
    call init_conventional_from_dcdft(lg,mg,system,info,spsi_in)
  end if
  Mit_restart_loaded = Mit
  if(yn_reset_step_restart=='y' ) Mit=0
  call timer_end(LOG_RESTART_SELF)
  call comm_sync_all
  call timer_end(LOG_RESTART_SYNC)
  if(yn_restart=='n') Mit=0

  if((gram_schmidt_interval == 0) ) then
    call gram_schmidt(system, mg, info, spsi_in)
  end if

  if(yn_jm=='n') then
    call calc_nlcc(pp, system, mg, ppn)
    if ((.not. quiet) .and. comm_is_root(nproc_id_global)) then
      write(*, '(1x, a, es23.15e3)') "Maximal rho_NLCC=", maxval(ppn%rho_nlcc)
      write(*, '(1x, a, es23.15e3)') "Maximal tau_NLCC=", maxval(ppn%tau_nlcc)
    end if
  end if
  
  call calc_density(system,rho_s,spsi_in,info,mg)
  rho%f = 0d0
  do jspin=1,system%nspin
     rho%f = rho%f + rho_s(jspin)%f
     rt%rho0_s(jspin)%f = rho_s(jspin)%f ! electron density @ t=0 (GS)
  end do
  if(yn_restart=='y')then
    Vh%f = 2.d0*Vh_stock1%f - Vh_stock2%f
    Vh_stock2%f = Vh_stock1%f
  end if
  spsi_in%update_zwf_overlap  = .false.
  spsi_out%update_zwf_overlap = .false.

  if(yn_jm=='y') rho%f = rho%f + rho_jm%f

  call hartree(lg,mg,info,system,fg,poisson,srg_scalar,stencil,rho,Vh)
  call exchange_correlation(system,xc_func,mg,srg_scalar,srg,rho_s,pp,ppn,info,spsi_in,stencil,Vxc,energy%E_xc)
  call update_vlocal(mg,system%nspin,Vh,Vpsl,Vxc,V_local)
  if(yn_restart=='y')then
    Vh_stock1%f=Vh%f
  else if(yn_restart=='n')then
    Vh_stock1%f=Vh%f
    Vh_stock2%f=Vh%f
  end if

  if(propagator=='aetrs')then
    do jspin=1,system%nspin
      rt%vloc_old(jspin,1)%f = V_local(jspin)%f
      rt%vloc_old(jspin,2)%f = V_local(jspin)%f
    end do
  end if
  
  allocate(energy%esp(system%no,system%nk,system%nspin))
  
  if(projection_option/='no') then
    call init_projection(system,lg,mg,info,stencil,Vpsl,xc_func,pp,ppn,fg,poisson,srg_scalar,rt,energy,ofl)
  end if
  
  call nvtxEndRange
  call timer_end(LOG_READ_GS_DATA)



  select case(iperiodic)
  case(0) 
     ewald%yn_bookkeep='n'  !to be input keyword??
  case(3)
     ewald%yn_bookkeep='y'
     call  init_nion_div(system,lg,mg,info)
  end select
  if(ewald%yn_bookkeep=='y') call init_ewald(system,info,ewald)
  
  ! calculation of GS total energy
  call calc_eigen_energy(energy,spsi_in,spsi_out,tpsi,system,info,mg,V_local,stencil,srg,ppg)
  if(yn_jm=='n') then
    rion_update = .true. ! it's first calculation
  else
    rion_update = .false.
  end if
  select case(iperiodic)
  case(0)
     call calc_Total_Energy_isolated(system,info,lg,mg,pp,ppg,fg,poisson,rho_s,Vh,Vxc,rion_update,energy)
  case(3)
     call calc_Total_Energy_periodic(mg,ewald,system,info,pp,ppg,fg,poisson,rion_update,energy)
  end select
  energy%E_tot0 = energy%E_tot
  
  call timer_begin(LOG_INIT_RT)
  
  ! 'Hartree' parameter
  poisson%iterVh = 0        ! Iteration counter
  
  call timer_end(LOG_INIT_RT)
  
  !Open output files and print header parts (Please move and put this kinds here!)
  !(output file names)
  if(comm_is_root(nproc_id_global))then
     write(ofl%file_rt_data,"(2A,'_rt.data')") trim(base_directory),trim(SYSname)
     write(ofl%file_rt_energy_data,"(2A,'_rt_energy.data')") trim(base_directory),trim(SYSname)
     write(ofl%file_response_data,"(2A,'_response.data')") trim(base_directory),trim(SYSname)
     write(ofl%file_pulse_data,"(2A,'_pulse.data')") trim(base_directory),trim(SYSname)
     if(yn_out_stress == 'y') write(ofl%file_stress_data,"(2A,'_stress.data')") trim(base_directory),trim(SYSname)
  endif
  call comm_bcast(ofl%file_rt_data,       nproc_group_global)
  call comm_bcast(ofl%file_rt_energy_data,nproc_group_global)
  call comm_bcast(ofl%file_response_data, nproc_group_global)
  call comm_bcast(ofl%file_pulse_data,    nproc_group_global)
  if(yn_out_stress == 'y') call comm_bcast(ofl%file_stress_data, nproc_group_global)
  
 ! write file header
  
  !(header of SYSname_rt.data)
  select case(iperiodic)
  case(0) ; call write_rt_data_0d(-1,ofl,dt,system,rt)
  case(3) ; call write_rt_data_3d(-1,ofl,dt,system,curr_e_tmp,curr_i_tmp)
  end select

  !(header of SYSname_rt_energy.data)
  call write_rt_energy_data(-1,ofl,dt,energy,md)
  if(yn_out_stress == 'y') call write_stress_rt(-1, ofl, dt, system, energy, pp)
  
  if(yn_spinorbit=='y') then
  !(header in SYSname_rt_spin.data)
    call write_rt_spin(-1,ofl,system,lg,mg,info,stencil,ppg,spsi_in)
  end if
  
  if(yn_out_current_decomposed == 'y') then
  !(header in SYSname_current_decomposed.data)
    call write_current_decomposed(-1,ofl,mg,system,info,stencil,srg,spsi_in,ppg)
  end if
  
  if(yn_md=='y' .or. yn_out_rvf_rt=='y')then
     call write_xyz(comment_line,"new","rvf",system,ofl)
  endif
  
  !---------------------------- time-evolution
  
  call timer_begin(LOG_INIT_TIME_PROPAGATION)
  
  idensity=0
  idiffDensity=1
  fileLaser= "laser.out"
  
  allocate( R1(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2), &
                                 lg%is(3):lg%ie(3)))
  
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
    call allocate_scalar(mg,rt%vonf)
    call set_vonf(mg,lg,system%Hgs,rt)
  end if
  
  if(yn_restart == 'n')then
    call calc_dip(info%icomm_r,lg,mg,system,rho_s,rbox_array2)
    rt%Dp0_e(1:3) = -rbox_array2(1:3) * system%Hgs(1:3) * system%Hvol
  end if
  if ((.not. quiet) .and. comm_is_root(nproc_id_global))then
    write(*,'(a30)', advance="no") "Static dipole moment(xyz) ="
    write(*,'(3e20.10)') (rt%Dp0_e(i1)*au_length_aa, i1=1,3)
    write(*,*)
  endif
  
  ! Initial wave function
  if(iperiodic==0 .and. ae_shape1 == 'impulse' .and. yn_restart == 'n')then
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
        call write_dns(lg,mg,system,info,rho_s,rt%rho0_s,itt)
      end if
      if(yn_out_elf_rt=='y')then
        call write_elf(itt,lg,mg,system,info,stencil,rho,srg,srg_scalar,spsi_in)
      end if
      if(yn_out_estatic_rt=='y')then
        call write_estatic(lg,mg,system,stencil,info,Vh,srg_scalar,itt)
      end if
    end do
  
  if(yn_md=='y') call init_md(system,md)
  
  ! preparation for projection
  if(projection_option/='no') then
    rt%E_old = energy%E_kin
  end if
  
  ! single-scale Maxwell-TDDFT
  singlescale%flag_use=.false.
  if(theory=='single_scale_maxwell_tddft') singlescale%flag_use=.true.

  if(singlescale%flag_use) then
    if(comm_is_root(nproc_id_global)) write(*,*) "single-scale Maxwell-TDDFT method"
    if(.not. stencil%if_orthogonal) stop "error: single-scale Maxwell-TDDFT & non-orthogonal lattice"
    call allocate_vector(mg,rt%j_e)
    call init_singlescale(mg,lg,info,system%hgs,system%rmatrix_B,rho,Vh &
    & ,srg_scalar,singlescale,system%Ac_micro,system%div_Ac)

    if(yn_out_dns_ac_je=='y')then
       itt=Mit
       call write_dns_ac_je(info,mg,system,rho%f,singlescale,itt,"new")
       call write_dns_ac_je(info,mg,system,rho%f,singlescale,itt,"bin")
    end if
    
  else if(yn_out_micro_je=='y') then
  ! not single-scale but output microscopic current rt%j_e
    call allocate_vector(mg,rt%j_e)
  end if

  !-------------------------------------------------- Time evolution
  
  !(force at initial step)
  if(yn_md=='y' .or. yn_out_rvf_rt=='y')then
     call calc_force(system,pp,fg,info,mg,stencil,poisson,srg,ppg,spsi_in,ewald,Vxc,ppn)
  
     !open trj file for coordinate, velocity, and force (rvf) in xyz format
     write(comment_line,10) -1, 0.0d0
     if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
          &  write(comment_line,12) trim(comment_line), md%xi_nh
     call write_xyz(comment_line,"add","rvf",system,ofl)
  10 format("#rt   step=",i8,"   time",e16.6)
  12 format(a,"  xi_nh=",e18.10)
  end if

  if(yn_out_stress == 'y' .and. theory /= 'multi_scale_maxwell_tddft') then
    if(yn_restart == 'n') then
      stress_label_iter = 0
      stress_state_iter = 0
    else if(yn_reset_step_restart == 'y') then
      stress_label_iter = 0
      stress_state_iter = Mit_restart_loaded
    else
      stress_label_iter = Mit_restart_loaded
      stress_state_iter = Mit_restart_loaded
    end if
    if(.not. singlescale%flag_use) system%vec_Ac(:) = rt%Ac_tot(:,stress_state_iter)
    call refresh_stress_output_state(system, theory, info, mg, stencil, srg, ppg, &
                                     energy, V_local, spsi_in, spsi_out, tpsi, field_state)
    call calc_stress(system, pp, fg, info, mg, stencil, poisson, srg, ppg, ppn, &
                     spsi_in, ewald, energy, xc_func, rho_s, Vxc, field_state)
    call write_stress_rt(stress_label_iter, ofl, dt, system, energy, pp)
  end if
  
  call timer_end(LOG_INIT_TIME_PROPAGATION)
  
  
  call timer_begin(LOG_INIT_RT)

  allocate(rt%zc(N_hamil))

  if(propagator=='aetrs')then
    do n=1,N_hamil
      rt%zc(n)=(-zi*dt/2.d0)**n
      do m=1,n
        rt%zc(n)=rt%zc(n)/m
      end do
    end do
  else
    do n=1,N_hamil
      rt%zc(n)=(-zi*dt)**n
      do m=1,n
        rt%zc(n)=rt%zc(n)/m
      end do
    end do
  end if


  deallocate (R1)

  call timer_end(LOG_INIT_RT)

  call nvtxEndRange

contains

subroutine init_code_optimization
  implicit none
  integer :: ignum(3)
  call nvtxStartRange('init_code_optimization', __LINE__)
  
  call switch_stencil_optimization(mg%num)
  call switch_openmp_parallelization(mg%num)

  if(iperiodic==3 .and. product(info%nprgrid)==1) then
     ignum = mg%num
  else
     ignum = mg%num + (nd*2)
  end if
  call set_modulo_tables(ignum)

  if ((.not. quiet) .and. comm_is_root(nproc_id_global)) then
     call optimization_log(info)
  end if
  call nvtxEndRange
end subroutine init_code_optimization

!---------------------------------------------------------------------------------------------------
! apply_nstate_active
!   Reduce the RT propagation orbital set (system%no) to the "active" orbitals (occupied or above
!   occ_threshold_rt) before the wavefunction arrays / orbital parallel-distribution / srg are
!   finalized. The full nstate is kept for projection (handled separately by init_projection).
!
!   Backward compatibility (HARD requirement): when both nstate_active<=0 and occ_threshold_rt<=0
!   (defaults) this routine returns immediately, leaving system%no, info and srg exactly as set by
!   init_dft -> the downstream code path is byte-identical to the previous behavior.
!
!   The "active" count is derived from the GS (full nstate) occupations, which are read here from the
!   global restart files (info.bin / occupation.bin). The GS file is eigenvalue-ordered, so the first
!   nact orbitals (index 1..nact) are the lowest-energy = occupied / above-threshold set, and the
!   existing read_bin path (called later in restart_rt) reads exactly those first nact orbitals.
!---------------------------------------------------------------------------------------------------
subroutine apply_nstate_active(system, info, srg, nact_out)
  use structures, only: s_dft_system, s_parallel_info, s_sendrecv_grid
  use communication, only: comm_is_root, comm_bcast
  implicit none
  type(s_dft_system)    :: system
  type(s_parallel_info) :: info
  type(s_sendrecv_grid) :: srg
  integer, intent(out)  :: nact_out

  integer :: nstate_full
  integer :: nact, no_old, ncopy
  real(8) :: dropped, ne_active
  real(8),allocatable :: rocc_full(:,:,:), roccbox(:,:,:)
  character(256) :: gdir, wdir, dir_file_in
  integer :: mk, mo, itt, nprocs, iu
  logical :: if_real_orbital_dummy
  integer :: ik, is, io
  integer :: nact_stored, nact_input

  nact_out = system%no

  ! Scope (initial version): tddft_pulse only. LR / multiscale are future work (spec section 10).
  if (theory /= 'tddft_pulse') return

  ! Backward-compatibility fast path: nothing requested -> leave everything untouched.
  if (nstate_active <= 0 .and. occ_threshold_rt <= 0d0) return

  ! full GS orbital count (projection basis size, and the # of orbitals stored in the GS restart)
  nstate_full = nstate
  no_old = system%no

  ! ---- restart-meta authority (Task 4) ----
  ! On an RT restart (yn_restart=='y') the propagated orbital count MUST equal the value fixed at the
  ! FIRST RT init, so RT->RT (self-)checkpoint chains stay consistent even if the user edits
  ! nstate_active / occ_threshold_rt in the restart input. The first-init value was written into the RT
  ! restart header as system%no (= nact); read it back and make it authoritative. We deliberately skip
  ! the occupation-based re-derivation on restart: it is both unnecessary (the meta is authoritative)
  ! and fragile on self-checkpoint chains (occupation.bin then lives in the per-rank dir, not gdir).
  if (yn_restart == 'y') then
    call read_restart_orbital_count(nact_stored)
    if (nact_stored >= 1) then
      ! input-derived value for the warn comparison (only the nstate_active branch needs no file read;
      ! the occ_threshold_rt branch can't be evaluated here without the occupations, so we skip the
      ! comparison in that case and just adopt the meta value).
      if (nstate_active > 0) then
        nact_input = min(nstate_active, nstate_full)
        if (nact_stored /= nact_input .and. comm_is_root(nproc_id_global)) then
          write(*,'(1x,a,i0,a,i0,a)') &
            '[nstate_active] WARN: restart input nstate_active=', nact_input, &
            ' differs from RT restart meta nact=', nact_stored, &
            '; using restart meta to keep the checkpoint chain consistent.'
        end if
      end if
      nact = nact_stored
      nact_out = nact
      if ((.not. quiet) .and. comm_is_root(nproc_id_global)) then
        write(*,'(1x,a,i0,a,i0,a)') &
          '[nstate_active] restart: propagate ', nact, ' of nstate=', nstate_full, &
          ' (from RT restart meta)'
      end if
      if (nact == no_old) return
      call resize_propagation_system(system, info, srg, nact)
      return
    end if
    ! nact_stored < 1: header unreadable; fall through to derive from input below (e.g. first RT
    ! continuation from a GS-only restart that predates this feature). Robust default.
  end if

  ! ---- read the GS (full nstate) occupations from the global restart files ----
  ! Directory selection mirrors restart_rt (info.bin / occupation.bin are global, root-written files).
  call generate_restart_directory_name(directory_read_data, gdir, wdir)

  allocate(rocc_full(nstate_full, system%nk, system%nspin))
  rocc_full = 0d0

  if (comm_is_root(nproc_id_global)) then
    ! info.bin : nk, no (= # orbitals stored in the GS restart), iter, nprocs, if_real_orbital
    dir_file_in = trim(gdir)//"info.bin"
    iu = 95
    open(iu, file=dir_file_in, form='unformatted', status='old')
    read(iu) mk
    read(iu) mo
    read(iu) itt
    read(iu) nprocs
    read(iu) if_real_orbital_dummy
    close(iu)
    ! occupation.bin holds a single record: rocc(1:mo, 1:mk, 1:nspin).
    ! Read it into a (mo,mk,nspin) buffer to preserve the (band,k,spin) layout, then copy the
    ! first ncopy=min(mo,nstate_full) bands into rocc_full (matches the GS energy ordering).
    allocate(roccbox(mo, mk, system%nspin))
    dir_file_in = trim(gdir)//"occupation.bin"
    open(iu, file=dir_file_in, form='unformatted', status='old')
    read(iu) roccbox(1:mo, 1:mk, 1:system%nspin)
    close(iu)
    ncopy = min(mo, nstate_full)
    rocc_full(1:ncopy, 1:min(mk,system%nk), 1:system%nspin) = &
         roccbox(1:ncopy, 1:min(mk,system%nk), 1:system%nspin)
    deallocate(roccbox)
  end if
  call comm_bcast(rocc_full, nproc_group_global)

  ! ---- derive the active orbital count from the (current) input + occupations (fresh RT) ----
  call calc_nstate_active(rocc_full, nstate_full, system%nk, system%nspin, system%wtk, &
                          nstate_active, occ_threshold_rt, nact, dropped)

  ! ---- safety valve: too few active orbitals (significant occupied weight discarded) -> abort ----
  if (dropped > 1d-3 * dble(nelec)) then
    if (comm_is_root(nproc_id_global)) then
      write(*,*) 'ERROR [nstate_active]: too few active orbitals; discarded occupation =', dropped
      write(*,*) '       (tolerance =', 1d-3*dble(nelec), ' = 1d-3 * nelec). ', &
                 'Increase nstate_active or lower occ_threshold_rt.'
    end if
    deallocate(rocc_full)
    call end_parallel
    stop
  end if

  ! ---- Ne sanity: sum of active occupation should reproduce nelec ----
  ne_active = 0d0
  do is = 1, system%nspin
    do ik = 1, system%nk
      do io = 1, nact
        ne_active = ne_active + rocc_full(io,ik,is) * system%wtk(ik)
      end do
    end do
  end do

  if ((.not. quiet) .and. comm_is_root(nproc_id_global)) then
    write(*,'(1x,a,i0,a,i0)')   '[nstate_active] propagate ', nact, ' of nstate=', nstate_full
    write(*,'(1x,a,es12.4)')    '[nstate_active] discarded occupation (dropped) =', dropped
    write(*,'(1x,a,f12.6,a,i0)')'[nstate_active] active electron number Ne =', ne_active, &
                                ' (nelec=', nelec
  end if

  deallocate(rocc_full)

  nact_out = nact

  ! Nothing to change if the propagation set already equals the active count
  ! (e.g. semiconductor without smearing already has system%no = nelec/2).
  if (nact == no_old) return

  ! ---- shrink/resize the propagation system to nact orbitals ----
  call resize_propagation_system(system, info, srg, nact)

  ! Task 4 (DONE): nact persistence across RT->RT checkpoint chains is handled via the existing RT
  !   restart header. write_bin already stores system%no (= nact) in info.bin; on restart the stored
  !   value is read back early (read_restart_orbital_count) above and made authoritative over the
  !   input-derived value (warn on mismatch). No separate meta record is needed.

end subroutine apply_nstate_active

!---------------------------------------------------------------------------------------------------
! resize_propagation_system
!   Shrink the RT propagation system to nact orbitals: set system%no, reallocate system%rocc, and
!   re-build the orbital parallel distribution (init_parallel_dft) and the wavefunction sendrecv
!   buffers (init_sendrecv_grid). Shared by the fresh-RT and the restart paths of apply_nstate_active.
!   References the host's mg (same as apply_nstate_active did inline).
!---------------------------------------------------------------------------------------------------
subroutine resize_propagation_system(system, info, srg, nact)
  use structures, only: s_dft_system, s_parallel_info, s_sendrecv_grid
  implicit none
  type(s_dft_system)    :: system
  type(s_parallel_info) :: info
  type(s_sendrecv_grid) :: srg
  integer, intent(in)   :: nact
  integer,dimension(2,3) :: neig

  system%no = nact

  if (allocated(system%rocc)) deallocate(system%rocc)
  allocate(system%rocc(system%no, system%nk, system%nspin))
  system%rocc = 0d0
  ! rocc is reloaded from the restart file later (read_bin), but keep it consistent here too.

  ! ---- re-distribute orbitals over the ob-communicator for the new system%no ----
  if (allocated(info%io_s_all)) deallocate(info%io_s_all)
  if (allocated(info%io_e_all)) deallocate(info%io_e_all)
  if (allocated(info%numo_all)) deallocate(info%numo_all)
  ! irank_io is self-guarded (deallocated inside init_parallel_dft)
  call init_parallel_dft(system, info)

  ! ---- re-create the wavefunction sendrecv buffers for the new orbital count ----
  call create_sendrecv_neig(neig, info)
  call init_sendrecv_grid(srg, mg, info%numo*info%numk*system%nspin, info%icomm_rko, neig)

end subroutine resize_propagation_system

end subroutine initialization_rt

end module initialization_rt_sub
