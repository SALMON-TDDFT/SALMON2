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

subroutine main_tddft
use math_constants, only: pi
use salmon_global
use structures
use parallelization, only: adjust_elapse_time, nproc_group_global
use communication, only: comm_is_root, comm_sync_all, comm_bcast, comm_summation
use salmon_xc, only: finalize_xc
use timer
use write_sub, only: write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d
use initialization_rt_sub
use checkpoint_restart_sub
use jellium, only: check_condition_jm
use rt_angular_momentum, only: write_local_angular_momentum_xy, flush_local_angular_momentum_xy
use rt_local_chern_marker, only: compute_local_chern_marker_from_orbital
use nvtx
use parallelization, only: nproc_id_global
implicit none

type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_dft_system)  :: system
type(s_rt) :: rt
type(s_parallel_info) :: info
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_xc_functional) :: xc_func
type(s_reciprocal_grid) :: fg
type(s_ewald_ion_ion) :: ewald
type(s_dft_energy) :: energy
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

integer :: Mit, itt
logical :: is_checkpoint_iter, is_shutdown_time, is_checkpoint

!check condition for using jellium model
if(yn_jm=='y') call check_condition_jm

call timer_begin(LOG_TOTAL)

call initialization_rt( Mit, system, energy, ewald, rt, md, &
                        singlescale,  &
                        stencil, fg, poisson,  &
                        lg, mg,   &
                        info,  &
                        xc_func, ofl,  &
                        srg, srg_scalar,  &
                        spsi_in, spsi_out, tpsi, rho, rho_jm, rho_s,  &
                        V_local, Vbox, Vh, Vh_stock1, Vh_stock2, Vxc, Vpsl,&
                        pp, ppg, ppn )

#ifdef __FUJITSU
call fapp_start('time_evol',1,0) ! performance profiling
#endif

call print_header()

if (yn_out_lcm_rt == 'y') then
  if (yn_dg_fragment_rt == 'y') stop 'yn_out_lcm_rt=y is not supported for DG-Fragment RT'
  call write_local_chern_marker_xy(Mit, mg, system, info, spsi_in)
end if
if (yn_out_lz_rt == 'y') then
  if (.not. singlescale%flag_use) stop 'yn_out_lz_rt=y requires theory=single_scale_maxwell_tddft'
  call write_local_angular_momentum_xy(Mit, lg, mg, system, info, singlescale, spsi_in)
end if

if (iperiodic == 3) then
  call write_initial_density_probe(system, info, mg, rho, rho_s, Vh, Vxc, Vpsl, 'full-initial-density')
end if

#ifdef USE_OPENACC
!$acc enter data copyin(rt, rt%zc)
!$acc enter data copyin(mg, mg%is, mg%ie)
!$acc enter data copyin(poisson)
!$acc enter data copyin(fg)
!$acc enter data copyin(lg)
!$acc enter data copyin(Vh, Vxc, Vpsl)
!$acc enter data copyin(ewald, pp, ppg)
#endif

call comm_sync_all
call timer_enable_sub
call timer_begin(LOG_RT_ITERATION)

! Check if DG-Fragment RT method is enabled
if (yn_dg_fragment_rt == 'y') then
  ! === DG-Fragment RT time evolution ===
  call time_evolution_dg_fragment(Mit, system, rt, info, lg, mg, stencil, xc_func, &
                                   srg, srg_scalar, fg, poisson, pp, ppg, ppn, rho, rho_s, Vh, Vxc, Vpsl, energy, &
                                   ofl, md, singlescale)
else
  ! === Standard real-space RT time evolution ===
  TE : do itt=Mit+1,nt
    call nvtxStartRange('main loop', itt)

  if(mod(itt,2)==1)then
    call time_evolution_step(Mit,nt,itt,lg,mg,system,rt,info,stencil,xc_func &
     & ,srg,srg_scalar,pp,ppg,ppn,spsi_in,spsi_out,tpsi,rho,rho_jm,rho_s,V_local,Vbox,Vh,Vh_stock1,Vh_stock2,Vxc &
     & ,Vpsl,fg,energy,ewald,md,ofl,poisson,singlescale)
	  else
	    call time_evolution_step(Mit,nt,itt,lg,mg,system,rt,info,stencil,xc_func &
	     & ,srg,srg_scalar,pp,ppg,ppn,spsi_out,spsi_in,tpsi,rho,rho_jm,rho_s,V_local,Vbox,Vh,Vh_stock1,Vh_stock2,Vxc &
	     & ,Vpsl,fg,energy,ewald,md,ofl,poisson,singlescale)
	  end if

      if (yn_out_lcm_rt == 'y') then
        if (mod(itt, out_lcm_rt_step) == 0) then
          if (mod(itt,2) == 1) then
            call write_local_chern_marker_xy(itt, mg, system, info, spsi_out)
          else
            call write_local_chern_marker_xy(itt, mg, system, info, spsi_in)
          end if
        end if
      end if
      if (yn_out_lz_rt == 'y') then
        if (mod(itt, out_lz_rt_step) == 0) then
          if (mod(itt,2) == 1) then
            call write_local_angular_momentum_xy(itt, lg, mg, system, info, singlescale, spsi_out)
          else
            call write_local_angular_momentum_xy(itt, lg, mg, system, info, singlescale, spsi_in)
          end if
        end if
      end if

	  is_checkpoint_iter = (checkpoint_interval >= 1) .and. (mod(itt,checkpoint_interval) == 0)
  is_shutdown_time   = (time_shutdown > 0d0) .and. (adjust_elapse_time(timer_now(LOG_TOTAL)) > time_shutdown)

  is_checkpoint = is_checkpoint_iter .or. is_shutdown_time
  call nvtxStartRange('comm_bcast', __LINE__)
  call comm_bcast(is_checkpoint,nproc_group_global)
  call nvtxEndRange

  call nvtxEndRange
  if(is_checkpoint) then
    if (is_shutdown_time .and. comm_is_root(info%id_rko)) then
      print *, 'shutdown the calculation, iter =', itt
    end if

    call timer_begin(LOG_CHECKPOINT_SYNC)
    call timer_begin(LOG_CHECKPOINT_SELF)
    if (mod(itt,2)==1) then
      call checkpoint_rt(lg,mg,system,info,spsi_out,itt,rt,Vh_stock1,Vh_stock2,singlescale)
    else
      call checkpoint_rt(lg,mg,system,info,spsi_in, itt,rt,Vh_stock1,Vh_stock2,singlescale)
    endif
    call timer_end(LOG_CHECKPOINT_SELF)
    call comm_sync_all
    call timer_end(LOG_CHECKPOINT_SYNC)

    if (is_shutdown_time) then
      exit TE
    end if
  endif

end do TE
end if  ! yn_dg_fragment_rt

if (yn_out_lz_rt == 'y') then
  call flush_local_angular_momentum_xy(system)
end if

call timer_end(LOG_RT_ITERATION)
call timer_disable_sub

#ifdef __FUJITSU
call fapp_stop('time_evol',1,0) ! performance profiling
#endif

close(030) ! laser


!--------------------------------- end of time-evolution

!------------ Writing part -----------


call timer_begin(LOG_WRITE_RT_RESULTS)

!
select case(iperiodic)
case(0)
  if(theory=="tddft_response")then
    call write_response_0d(ofl,rt)
  else
    call write_pulse_0d(ofl,rt)
  end if
case(3)
  if(theory=="tddft_response")then
    call write_response_3d(ofl,rt)
  else
    call write_pulse_3d(ofl,rt)
  end if
end select

if(comm_is_root(nproc_id_global))then
  close(ofl%fh_rt)  ! Close _rt.data file
end if

call timer_end(LOG_WRITE_RT_RESULTS)
call timer_end(LOG_TOTAL)

if(write_rt_wfn_k=='y')then
  call checkpoint_rt(lg,mg,system,info,spsi_out,Mit,rt,Vh_stock1,Vh_stock2,singlescale,ofl%dir_out_restart)
end if

call finalize_xc(xc_func)

contains

subroutine time_evolution_dg_fragment(Mit, system, rt, info, lg, mg, stencil, xc_func, &
                                       srg, srg_scalar, fg, poisson, pp, ppg, ppn, rho, rho_s, Vh, Vxc, Vpsl, energy, ofl, md, singlescale)
  use structures
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_fragment, only: init_dg_fragment_rt_std => init_dg_fragment_rt, &
                            tddft_dg_fragment_iteration_std => tddft_dg_fragment_iteration, &
                            finalize_dg_fragment_rt_std => finalize_dg_fragment_rt, &
                            calculate_hamiltonian_matrix_std => calculate_hamiltonian_matrix, &
                            print_initial_electron_probe_std => print_initial_electron_probe, &
                            calculate_observables_std => calculate_observables
  use rt_dg_fragment_soi, only: init_dg_fragment_rt_soi => init_dg_fragment_rt, &
                                tddft_dg_fragment_iteration_soi => tddft_dg_fragment_iteration, &
                                finalize_dg_fragment_rt_soi => finalize_dg_fragment_rt, &
                                calculate_hamiltonian_matrix_soi => calculate_hamiltonian_matrix
  use communication, only: comm_is_root
  use parallelization, only: nproc_id_global
  use sendrecv_grid, only: s_sendrecv_grid
  use salmon_xc, only: s_xc_functional
  use write_sub
  use salmon_global, only: theory, method_singlescale, yn_ffte, yn_jm, yn_spinorbit
  use fdtd_coulomb_gauge, only: fdtd_singlescale, fourier_singlescale
  use hamiltonian, only: update_kvector_nonlocalpt_microAc
  implicit none
  integer,                intent(in)    :: Mit
  type(s_dft_system),     intent(inout) :: system
  type(s_rt),             intent(inout) :: rt
  type(s_parallel_info),  intent(in)    :: info
  type(s_rgrid),          intent(in)    :: lg, mg
  type(s_stencil),        intent(in)    :: stencil
  type(s_xc_functional),  intent(in)    :: xc_func
  type(s_sendrecv_grid),  intent(inout) :: srg, srg_scalar
  type(s_reciprocal_grid),intent(in)    :: fg
  type(s_poisson),        intent(inout) :: poisson
  type(s_pp_info),        intent(in)    :: pp
  type(s_pp_grid),        intent(in)    :: ppg
  type(s_pp_nlcc),        intent(in)    :: ppn
  type(s_scalar),         intent(inout) :: rho, Vh, Vpsl
  type(s_scalar),         intent(inout) :: rho_s(system%nspin), Vxc(system%nspin)
  type(s_dft_energy),     intent(inout) :: energy
  type(s_ofile),          intent(inout) :: ofl
  type(s_md),             intent(inout) :: md
  type(s_singlescale),    intent(inout) :: singlescale
  
  type(s_dg_fragment_rt) :: dg_frag
  integer :: itt
  logical :: enable_main_trace
  character(len=32) :: env_main_trace
  integer :: env_main_len, env_main_stat

  enable_main_trace = .false.
  env_main_trace = ''
  call get_environment_variable('SALMON_DG_MAIN_TRACE', env_main_trace, length=env_main_len, status=env_main_stat)
  if (env_main_stat == 0 .and. env_main_len > 0) then
    if (env_main_trace(1:1) == '1' .or. env_main_trace(1:1) == 'y' .or. env_main_trace(1:1) == 'Y' .or. &
        env_main_trace(1:1) == 't' .or. env_main_trace(1:1) == 'T') then
      enable_main_trace = .true.
    end if
  end if
  
  ! Initialize DG-Fragment RT
  if (yn_spinorbit == 'y') then
    call init_dg_fragment_rt_soi(dg_frag, system, rt, info, lg, mg, ppg)
  else
    call init_dg_fragment_rt_std(dg_frag, system, rt, info, lg, mg, ppg)
  end if
  
  ! Calculate Hamiltonian matrix with initial potentials
  ! Note: This must be done after init when stencil, pp, ppg, and potentials are available
  if (yn_spinorbit == 'y') then
    call calculate_hamiltonian_matrix_soi(dg_frag, system, lg, mg, stencil, Vh, Vxc, Vpsl, pp, ppg)
  else
    call calculate_hamiltonian_matrix_std(dg_frag, system, lg, mg, stencil, Vh, Vxc, Vpsl, pp, ppg)
  end if
  
  ! H_mat_kinetic is constructed inside calculate_hamiltonian_matrix
  
  if (comm_is_root(nproc_id_global)) then
    write(*,*)
    write(*,*) "=== Starting DG-Fragment RT time evolution ==="
    write(*,*)
  end if

  if (iperiodic == 3) then
    call write_initial_density_probe(system, info, mg, rho, rho_s, Vh, Vxc, Vpsl, 'dg-initial-density')
  end if
  if (yn_spinorbit /= 'y') then
    if (enable_main_trace) then
      write(*,'(1x,a,a)') "        main-trace: stage=", "before-initial-electron-probe"
      flush(6)
    end if
    call print_initial_electron_probe_std(dg_frag, system, mg, rho)
    if (enable_main_trace) then
      write(*,'(1x,a,a)') "        main-trace: stage=", "after-initial-electron-probe"
      flush(6)
      write(*,'(1x,a,a)') "        main-trace: stage=", "before-initial-observables"
      flush(6)
    end if
    call calculate_observables_std(dg_frag, system, mg, stencil, ppg, rt, Mit, Vh, Vxc, Vpsl)
    if (enable_main_trace) then
      write(*,'(1x,a,a)') "        main-trace: stage=", "after-initial-observables"
      flush(6)
      write(*,'(1x,a,a)') "        main-trace: stage=", "before-initial-update-dg-energy"
      flush(6)
    end if
    call update_dg_rt_total_energy(system, info, mg, fg, poisson, ppg, rho, rho_s, Vh, Vxc, Vpsl, dg_frag, energy, Mit)
    if (enable_main_trace) then
      write(*,'(1x,a,a)') "        main-trace: stage=", "after-initial-update-dg-energy"
      flush(6)
    end if
  end if

  ! Time evolution loop
  do itt = Mit+1, nt
    if (yn_spinorbit == 'y') then
      call tddft_dg_fragment_iteration_soi(dg_frag, system, info, rt, itt, dt, &
                                           lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                           rho, rho_s, Vh, Vxc, Vpsl, energy)
    else
      call tddft_dg_fragment_iteration_std(dg_frag, system, info, rt, itt, dt, &
                                           lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                           rho, rho_s, Vh, Vxc, Vpsl, energy)
    end if
    
    ! Reconstruct the DFT total energy from the DG one-body expectation.
    if (enable_main_trace .and. itt <= 2) then
      write(*,'(1x,a,i0,a,a)') "        main-trace: itt=", itt, " stage=", "before-update-dg-energy"
      flush(6)
    end if
    call update_dg_rt_total_energy(system, info, mg, fg, poisson, ppg, rho, rho_s, Vh, Vxc, Vpsl, dg_frag, energy, itt)
    if (enable_main_trace .and. itt <= 2) then
      write(*,'(1x,a,i0,a,a)') "        main-trace: itt=", itt, " stage=", "after-update-dg-energy"
      flush(6)
    end if
    energy%elec_num = dg_frag%elec_num_scaled
    energy%elec_num_raw = dg_frag%elec_num_raw
    energy%pw_weight_raw = dg_frag%pw_weight_raw
    energy%E_onebody_mix = dg_frag%total_energy_mix

    if (theory == 'single_scale_maxwell_tddft') then
      singlescale%E_electron = energy%E_tot
      if (method_singlescale == '1d_fourier' .and. yn_ffte == 'y') then
        call fourier_singlescale(lg, mg, info, fg, rho, rt%j_e, Vh, poisson, singlescale)
      else
        call fdtd_singlescale(itt, lg, mg, system, info, rho, Vh, rt%j_e, srg_scalar, system%Ac_micro, system%div_Ac, singlescale)
      end if
      if (yn_jm == 'n') call update_kvector_nonlocalpt_microAc(info%ik_s, info%ik_e, system, ppg)
      rt%curr(1:3, itt) = singlescale%curr_ave(1:3)
    end if
    
    ! Write observable data (current and energy)
    ! DG-Fragment uses coefficient-space evolution, so observables are already computed
    ! Keep output fields consistent with conventional RT
    system%vec_Ac(1:3) = rt%Ac_tot(1:3, itt)
    system%vec_Ac_ext(1:3) = rt%Ac_ext(1:3, itt)
    system%vec_E(1:3) = rt%E_tot(1:3, itt)
    system%vec_E_ext(1:3) = rt%E_ext(1:3, itt)
    energy%E_work_elec = energy%E_work_elec - dot_product(rt%E_tot(1:3, itt), dg_frag%current(1:3)) * system%det_a * dt
    select case(iperiodic)
    case(0)
      ! Isolated system: output dipole moment
      call write_rt_data_0d(itt, ofl, dt, system, rt)
    case(3)
      ! Periodic system: output current density
      ! DG-Fragment current from calculate_observables
      ! curr_e needs shape (3, 2) for two spins, use same current for both spins
      if (enable_main_trace .and. itt <= 2) then
        write(*,'(1x,a,i0,a,a)') "        main-trace: itt=", itt, " stage=", "before-write-rt-data-3d"
        flush(6)
      end if
      call write_rt_data_3d(itt, ofl, dt, system, &
                            reshape([dg_frag%current(1), dg_frag%current(2), dg_frag%current(3), &
                                     dg_frag%current(1), dg_frag%current(2), dg_frag%current(3)], [3,2]), &
                            dg_frag%current(1:3))
      if (enable_main_trace .and. itt <= 2) then
        write(*,'(1x,a,i0,a,a)') "        main-trace: itt=", itt, " stage=", "after-write-rt-data-3d"
        flush(6)
      end if
    end select
    
    ! Write energy data
    if (enable_main_trace .and. itt <= 2) then
      write(*,'(1x,a,i0,a,a)') "        main-trace: itt=", itt, " stage=", "before-write-rt-energy"
      flush(6)
    end if
    call write_rt_energy_data(itt, ofl, dt, energy, md)
    if (enable_main_trace .and. itt <= 2) then
      write(*,'(1x,a,i0,a,a)') "        main-trace: itt=", itt, " stage=", "after-write-rt-energy"
      flush(6)
    end if
    
    ! Output progress
    if (comm_is_root(nproc_id_global) .and. mod(itt, 10) == 0) then
      write(*,'(1x,i10,f12.3,3e16.6e3,e20.10e3)') itt, dble(itt)*dt, dg_frag%current(:), dg_frag%total_energy
    end if
  end do
  
  ! Write response spectra (dielectric function, HHG)
  ! call write_dg_fragment_response_data(rt, system, ofl)  ! Not implemented yet
  
  ! Finalize
  if (yn_spinorbit == 'y') then
    call finalize_dg_fragment_rt_soi(dg_frag)
  else
    call finalize_dg_fragment_rt_std(dg_frag)
  end if
  
  if (comm_is_root(nproc_id_global)) then
    write(*,*)
    write(*,*) "=== DG-Fragment RT time evolution completed ==="
    write(*,*)
  end if
  
end subroutine time_evolution_dg_fragment

subroutine update_dg_rt_total_energy(system, info, mg, fg, poisson, ppg, rho, rho_s, Vh, Vxc, Vpsl, dg_frag, energy, itt)
  use structures, only: s_dft_system, s_parallel_info, s_rgrid, s_reciprocal_grid, s_poisson, s_pp_grid, &
                        s_scalar, s_dft_energy
  use communication, only: comm_summation, comm_is_root, comm_get_max
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use math_constants, only: zi
  use salmon_global, only: yn_jm
  use parallelization, only: nproc_id_global
  implicit none
  type(s_dft_system),     intent(in)    :: system
  type(s_parallel_info),  intent(in)    :: info
  type(s_rgrid),          intent(in)    :: mg
  type(s_reciprocal_grid),intent(in)    :: fg
  type(s_poisson),        intent(in)    :: poisson
  type(s_pp_grid),        intent(in)    :: ppg
  type(s_scalar),         intent(in)    :: rho, Vh, Vpsl
  type(s_scalar),         intent(in)    :: rho_s(system%nspin), Vxc(system%nspin)
  type(s_dg_fragment_rt), intent(in)    :: dg_frag
  type(s_dft_energy),     intent(inout) :: energy
  integer,                intent(in)    :: itt
  integer, parameter :: n_gprobe = 3
  integer, parameter :: gprobe_idx(3, n_gprobe) = reshape((/ 1,1,1, 2,1,1, 3,1,1 /), (/3, n_gprobe/))
  integer :: ix, iy, iz, ispin, ia, iprobe
  real(8) :: rho_vh_local, rho_vh_sum, rho_vh_sum_r
  real(8) :: rho_vxc_local, rho_vxc_sum, rho_vxc_sum_r
  real(8) :: rho_vpsl_local, rho_vpsl_sum, rho_vpsl_sum_r
  real(8) :: rho2_local, rho2_sum, rho2_sum_r
  real(8) :: rho_max_local(1), rho_max_sum(1), rho_max_sum_r(1)
  real(8) :: rho_g2_local, rho_g2_sum, rho_g2_sum_r
  real(8) :: rho_gmax_local(1), rho_gmax_sum(1), rho_gmax_sum_r(1)
  real(8) :: rho_gprobe_re_local(n_gprobe), rho_gprobe_im_local(n_gprobe)
  real(8) :: rho_gprobe_re_sum(n_gprobe), rho_gprobe_im_sum(n_gprobe)
  real(8) :: rho_gprobe_re_sum_r(n_gprobe), rho_gprobe_im_sum_r(n_gprobe)
  real(8) :: E_wrk(3), E_sum(3), E_sum_r(3), etmp, sysvol, g(3), r(3), Gd
  complex(8) :: rho_e, rho_i
  logical, parameter :: enable_dg_energy_helper_trace = .false.

  rho_vh_local = 0.0d0
  rho_vxc_local = 0.0d0
  rho_vpsl_local = 0.0d0
  rho2_local = 0.0d0
  rho_max_local(1) = 0.0d0
  rho_g2_local = 0.0d0
  rho_gmax_local(1) = 0.0d0
  rho_gprobe_re_local(:) = 0.0d0
  rho_gprobe_im_local(:) = 0.0d0
  E_wrk(:) = 0.0d0
  E_sum(:) = 0.0d0
  E_sum_r(:) = 0.0d0
  etmp = 0.0d0

  do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
      do ix = mg%is(1), mg%ie(1)
        rho_max_local(1) = max(rho_max_local(1), rho%f(ix, iy, iz))
        rho2_local = rho2_local + rho%f(ix, iy, iz) * rho%f(ix, iy, iz)
        rho_vh_local = rho_vh_local + rho%f(ix, iy, iz) * Vh%f(ix, iy, iz)
        rho_vpsl_local = rho_vpsl_local + rho%f(ix, iy, iz) * Vpsl%f(ix, iy, iz)
        do ispin = 1, system%nspin
          rho_vxc_local = rho_vxc_local + rho_s(ispin)%f(ix, iy, iz) * Vxc(ispin)%f(ix, iy, iz)
        end do
      end do
    end do
  end do

  rho_vh_local = rho_vh_local * system%Hvol
  rho_vxc_local = rho_vxc_local * system%Hvol
  rho_vpsl_local = rho_vpsl_local * system%Hvol
  rho2_local = rho2_local * system%Hvol

  call comm_summation(rho_vh_local, rho_vh_sum, dg_frag%icomm)
  call comm_summation(rho_vh_local, rho_vh_sum_r, info%icomm_r)
  call comm_summation(rho_vxc_local, rho_vxc_sum, dg_frag%icomm)
  call comm_summation(rho_vxc_local, rho_vxc_sum_r, info%icomm_r)
  call comm_summation(rho_vpsl_local, rho_vpsl_sum, dg_frag%icomm)
  call comm_summation(rho_vpsl_local, rho_vpsl_sum_r, info%icomm_r)
  call comm_summation(rho2_local, rho2_sum, dg_frag%icomm)
  call comm_summation(rho2_local, rho2_sum_r, info%icomm_r)
  call comm_get_max(rho_max_local, rho_max_sum, 1, dg_frag%icomm)
  call comm_get_max(rho_max_local, rho_max_sum_r, 1, info%icomm_r)

  sysvol = system%det_a
  do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
      do ix = mg%is(1), mg%ie(1)
        g(1) = fg%vec_G(1, ix, iy, iz)
        g(2) = fg%vec_G(2, ix, iy, iz)
        g(3) = fg%vec_G(3, ix, iy, iz)
        rho_e = poisson%zrhoG_ele(ix, iy, iz)
        rho_g2_local = rho_g2_local + sysvol * fg%coef(ix, iy, iz) * abs(rho_e)**2
        rho_gmax_local(1) = max(rho_gmax_local(1), abs(rho_e))
        do iprobe = 1, n_gprobe
          if (ix == gprobe_idx(1, iprobe) .and. iy == gprobe_idx(2, iprobe) .and. iz == gprobe_idx(3, iprobe)) then
            rho_gprobe_re_local(iprobe) = real(rho_e)
            rho_gprobe_im_local(iprobe) = aimag(rho_e)
          end if
        end do
        E_wrk(1) = E_wrk(1) + sysvol * fg%coef(ix, iy, iz) * (abs(rho_e)**2 * 0.5d0)
        if (yn_jm == 'n') then
          rho_i = ppg%zrhoG_ion(ix, iy, iz)
          E_wrk(2) = E_wrk(2) + sysvol * fg%coef(ix, iy, iz) * (-rho_e * conjg(rho_i))
          do ia = info%ia_s, info%ia_e
            r(:) = system%Rion(:, ia)
            Gd = g(1) * r(1) + g(2) * r(2) + g(3) * r(3)
            etmp = etmp + conjg(rho_e) * ppg%zVG_ion(ix, iy, iz, system%kion(ia)) * exp(-zi * Gd)
          end do
        end if
      end do
    end do
  end do

  call comm_summation(etmp, E_wrk(3), info%icomm_ko)
  call comm_summation(E_wrk, E_sum, 3, dg_frag%icomm)
  call comm_summation(E_wrk, E_sum_r, 3, info%icomm_r)
  call comm_summation(rho_g2_local, rho_g2_sum, dg_frag%icomm)
  call comm_summation(rho_g2_local, rho_g2_sum_r, info%icomm_r)
  call comm_summation(rho_gprobe_re_local, rho_gprobe_re_sum, n_gprobe, dg_frag%icomm)
  call comm_summation(rho_gprobe_im_local, rho_gprobe_im_sum, n_gprobe, dg_frag%icomm)
  call comm_summation(rho_gprobe_re_local, rho_gprobe_re_sum_r, n_gprobe, info%icomm_r)
  call comm_summation(rho_gprobe_im_local, rho_gprobe_im_sum_r, n_gprobe, info%icomm_r)
  call comm_get_max(rho_gmax_local, rho_gmax_sum, 1, dg_frag%icomm)
  call comm_get_max(rho_gmax_local, rho_gmax_sum_r, 1, info%icomm_r)

  if (dg_frag%isize_frag > 1) then
    rho_vh_sum = rho_vh_sum_r
    rho_vxc_sum = rho_vxc_sum_r
    rho_vpsl_sum = rho_vpsl_sum_r
    rho2_sum = rho2_sum_r
    rho_g2_sum = rho_g2_sum_r
    rho_gprobe_re_sum(:) = rho_gprobe_re_sum_r(:)
    rho_gprobe_im_sum(:) = rho_gprobe_im_sum_r(:)
    rho_gmax_sum(1) = rho_gmax_sum_r(1)
    E_sum(:) = E_sum_r(:)
  end if

  energy%E_kin = dg_frag%energy_kinetic
  energy%E_ion_nloc = dg_frag%energy_nonlocal
  energy%E_h = E_sum(1)
  energy%E_ion_loc = E_sum(2) + E_sum(3)
  energy%E_tot = dg_frag%total_energy - rho_vh_sum - rho_vxc_sum + energy%E_h + energy%E_xc + energy%E_ion_ion
  if (enable_dg_energy_helper_trace .and. comm_is_root(nproc_id_global) .and. (itt == 1 .or. mod(itt, 10) == 0)) then
    write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
      "        dg-energy-helper: itt=", itt, " E_one=", dg_frag%total_energy, " rhoVh=", rho_vh_sum, " rhoVxc=", rho_vxc_sum, &
      " rhoVpsl=", rho_vpsl_sum, " rho2=", rho2_sum, " rhomax=", rho_max_sum(1), &
      " rhoG2=", rho_g2_sum, " rhoGmax=", rho_gmax_sum(1), " rhoVh_r=", rho_vh_sum_r, &
      " rhoVxc_r=", rho_vxc_sum_r, " rhoVpsl_r=", rho_vpsl_sum_r, " rho2_r=", rho2_sum_r, &
      " rhomax_r=", rho_max_sum_r(1), " rhoG2_r=", rho_g2_sum_r, " rhoGmax_r=", rho_gmax_sum_r(1), &
      " E_kin=", energy%E_kin, &
      " E_ion_nloc=", energy%E_ion_nloc, " E_h=", energy%E_h, " E_xc=", energy%E_xc, &
      " E_ion_loc=", energy%E_ion_loc, " E_ion_ion=", energy%E_ion_ion, " E_tot=", energy%E_tot
    write(*,'(1x,a,i0,6(a,1pe14.6),3(a,1pe14.6))') &
      "        dg-rhoG-probe: itt=", itt, &
      " g111_re=", rho_gprobe_re_sum(1), " g111_im=", rho_gprobe_im_sum(1), &
      " g211_re=", rho_gprobe_re_sum(2), " g211_im=", rho_gprobe_im_sum(2), &
      " g311_re=", rho_gprobe_re_sum(3), " g311_im=", rho_gprobe_im_sum(3), &
      " g111_abs=", sqrt(rho_gprobe_re_sum(1)**2 + rho_gprobe_im_sum(1)**2), &
      " g211_abs=", sqrt(rho_gprobe_re_sum(2)**2 + rho_gprobe_im_sum(2)**2), &
      " g311_abs=", sqrt(rho_gprobe_re_sum(3)**2 + rho_gprobe_im_sum(3)**2)
    write(*,'(1x,a,i0,6(a,1pe14.6),3(a,1pe14.6))') &
      "        dg-rhoG-probe-r: itt=", itt, &
      " g111_re=", rho_gprobe_re_sum_r(1), " g111_im=", rho_gprobe_im_sum_r(1), &
      " g211_re=", rho_gprobe_re_sum_r(2), " g211_im=", rho_gprobe_im_sum_r(2), &
      " g311_re=", rho_gprobe_re_sum_r(3), " g311_im=", rho_gprobe_im_sum_r(3), &
      " g111_abs=", sqrt(rho_gprobe_re_sum_r(1)**2 + rho_gprobe_im_sum_r(1)**2), &
      " g211_abs=", sqrt(rho_gprobe_re_sum_r(2)**2 + rho_gprobe_im_sum_r(2)**2), &
      " g311_abs=", sqrt(rho_gprobe_re_sum_r(3)**2 + rho_gprobe_im_sum_r(3)**2)
    write(*,'(1x,a,i0,3(a,1pe14.6))') "        dg-rhoG-probe-gvec: itt=", itt, &
      " g211_gx=", fg%vec_G(1, 2, 1, 1), " g211_gy=", fg%vec_G(2, 2, 1, 1), " g211_gz=", fg%vec_G(3, 2, 1, 1)
    flush(6)
  end if
end subroutine update_dg_rt_total_energy

subroutine write_initial_density_probe(system, info, mg, rho, rho_s, Vh, Vxc, Vpsl, label)
  use structures, only: s_dft_system, s_parallel_info, s_rgrid, s_scalar
  use communication, only: comm_summation, comm_get_max, comm_is_root
  implicit none
  type(s_dft_system),    intent(in) :: system
  type(s_parallel_info), intent(in) :: info
  type(s_rgrid),         intent(in) :: mg
  type(s_scalar),        intent(in) :: rho, Vh, Vpsl
  type(s_scalar),        intent(in) :: rho_s(system%nspin), Vxc(system%nspin)
  character(*),          intent(in) :: label
  integer :: ix, iy, iz, ispin
  real(8) :: rho_vh_local, rho_vh_sum
  real(8) :: rho_vxc_local, rho_vxc_sum
  real(8) :: rho_vpsl_local, rho_vpsl_sum
  real(8) :: rho2_local, rho2_sum
  real(8) :: rho_max_local(1), rho_max_sum(1)

  rho_vh_local = 0.0d0
  rho_vxc_local = 0.0d0
  rho_vpsl_local = 0.0d0
  rho2_local = 0.0d0
  rho_max_local(1) = 0.0d0

  do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
      do ix = mg%is(1), mg%ie(1)
        rho_max_local(1) = max(rho_max_local(1), rho%f(ix, iy, iz))
        rho2_local = rho2_local + rho%f(ix, iy, iz) * rho%f(ix, iy, iz)
        rho_vh_local = rho_vh_local + rho%f(ix, iy, iz) * Vh%f(ix, iy, iz)
        rho_vpsl_local = rho_vpsl_local + rho%f(ix, iy, iz) * Vpsl%f(ix, iy, iz)
        do ispin = 1, system%nspin
          rho_vxc_local = rho_vxc_local + rho_s(ispin)%f(ix, iy, iz) * Vxc(ispin)%f(ix, iy, iz)
        end do
      end do
    end do
  end do

  rho_vh_local = rho_vh_local * system%Hvol
  rho_vxc_local = rho_vxc_local * system%Hvol
  rho_vpsl_local = rho_vpsl_local * system%Hvol
  rho2_local = rho2_local * system%Hvol

  call comm_summation(rho_vh_local, rho_vh_sum, info%icomm_r)
  call comm_summation(rho_vxc_local, rho_vxc_sum, info%icomm_r)
  call comm_summation(rho_vpsl_local, rho_vpsl_sum, info%icomm_r)
  call comm_summation(rho2_local, rho2_sum, info%icomm_r)
  call comm_get_max(rho_max_local, rho_max_sum, 1, info%icomm_r)

  if (comm_is_root(nproc_id_global)) then
    write(*,'(1x,a,a,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
      '        ', trim(label), ': rhoVh=', rho_vh_sum, ' rhoVxc=', rho_vxc_sum, &
      ' rhoVpsl=', rho_vpsl_sum, ' rho2=', rho2_sum, ' rhomax=', rho_max_sum(1)
    flush(6)
  end if
end subroutine write_initial_density_probe

subroutine write_local_chern_marker_xy(itt, mg, system, info, psi_fin)
  use structures, only: s_rgrid, s_dft_system, s_parallel_info, s_orbital
  use communication, only: comm_is_root, comm_summation
  use rt_local_chern_marker, only: compute_local_chern_marker_from_orbital
  use rt_local_chern_marker_soi, only: compute_local_chern_marker_from_orbital_soi => compute_local_chern_marker_from_orbital
  use filesystem, only: create_directory
  use inputoutput, only: t_unit_length
  use salmon_global, only: yn_spinorbit
  implicit none
  integer, intent(in) :: itt
  type(s_rgrid), intent(in) :: mg
  type(s_dft_system), intent(in) :: system
  type(s_parallel_info), intent(in) :: info
  type(s_orbital), intent(in) :: psi_fin
  real(8), allocatable :: marker(:,:,:)
  real(8), allocatable :: marker_xy(:,:), marker_xy_full_local(:,:), marker_xy_full(:,:)
  character(256) :: filename, filenum, map_directory
  integer :: ix, iy, iz, iunit, nx, ny

  allocate(marker(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
  allocate(marker_xy(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2)))
  if (yn_spinorbit == 'y') then
    call compute_local_chern_marker_from_orbital_soi(mg, system, info, psi_fin, marker)
  else
    call compute_local_chern_marker_from_orbital(mg, system, info, psi_fin, marker)
  end if

  marker_xy(:,:) = 0.0d0
  do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
      do ix = mg%is(1), mg%ie(1)
        marker_xy(ix,iy) = marker_xy(ix,iy) + marker(ix,iy,iz) * system%hgs(3)
      end do
    end do
  end do

  nx = maxval(mg%ie_all(1,:))
  ny = maxval(mg%ie_all(2,:))
  allocate(marker_xy_full_local(nx, ny), marker_xy_full(nx, ny))
  marker_xy_full_local(:,:) = 0.0d0
  marker_xy_full(:,:) = 0.0d0
  do iy = mg%is(2), mg%ie(2)
    do ix = mg%is(1), mg%ie(1)
      marker_xy_full_local(ix,iy) = marker_xy(ix,iy)
    end do
  end do
  call comm_summation(marker_xy_full_local, marker_xy_full, nx*ny, info%icomm_r)

  if (comm_is_root(nproc_id_global)) then
    write(filenum, '(i6.6)') itt
    map_directory = trim(base_directory)//trim(sysname)//'_lcm_xy/'
    call create_directory(trim(map_directory))
    filename = trim(map_directory)//trim(sysname)//'_lcm_xy_'//trim(adjustl(filenum))//'.data'
    open(newunit=iunit, file=trim(filename), status='replace', action='write')
    write(iunit,'(a)') '# Local Chern marker integrated along z'
    write(iunit,'(a)') '# x: x coordinate'
    write(iunit,'(a)') '# y: y coordinate'
    write(iunit,'(a)') '# local_chern_marker_zint: local Chern marker integrated over z'
    write(iunit,'(a,a,a)') '# 1:x[', trim(t_unit_length%name), '] 2:y[', trim(t_unit_length%name), '] 3:local_chern_marker_zint[none]'
    do iy = 1, ny
      do ix = 1, nx
        write(iunit,'(3(1x,es24.16))') dble(ix-1) * system%hgs(1) * t_unit_length%conv, &
                                       dble(iy-1) * system%hgs(2) * t_unit_length%conv, marker_xy_full(ix,iy)
      end do
      write(iunit,*)
    end do
    close(iunit)
  end if
  deallocate(marker_xy_full, marker_xy_full_local, marker_xy, marker)
end subroutine write_local_chern_marker_xy

subroutine print_header()
  use parallelization, only: nproc_id_global
  implicit none
  !(header of standard output)
  if(comm_is_root(nproc_id_global))then
    write(*,*)
    select case(iperiodic)
    case(0)
      write(*,'(1x,a10,a11,a48,a15,a18,a10)') &
                  "time-step ", "time[fs]",   &
                  "Dipole moment(xyz)[A]"     &
                ,"electrons", "Total energy[eV]", "iterVh"
    case(3)
      if(yn_jm=='n')then
        write(*,'(1x,a10,a11,a48,a15,a18)')   &
                    "time-step", "time[fs] ", &
                    "Current(xyz)[a.u.]",     &
                    "electrons", "Total energy[eV] "
      else
        write(*,'(1x,a10,a11,a48,a15,a18)')   &
                    "time-step", "time[fs] ", &
                    "Current(xyz)[a.u.]",     &
                    "electrons"
      end if
    end select
    write(*,'("#",7("----------"))')
  endif
end subroutine print_header

end subroutine main_tddft
