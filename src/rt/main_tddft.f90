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
use write_sub, only: write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d, &
  write_dg_polarization_data
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

if (iperiodic == 3 .and. yn_dg_fragment_rt /= 'y') then
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
                                   ofl, md, singlescale, spsi_in)
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
                                       srg, srg_scalar, fg, poisson, pp, ppg, ppn, rho, rho_s, &
                                       Vh, Vxc, Vpsl, energy, ofl, md, singlescale, spsi_restart)
  use structures
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_fragment, only: init_dg_fragment_rt_std => init_dg_fragment_rt, &
                            tddft_dg_fragment_iteration_std => tddft_dg_fragment_iteration, &
                            finalize_dg_fragment_rt_std => finalize_dg_fragment_rt, &
                            calculate_hamiltonian_matrix_std => calculate_hamiltonian_matrix, &
                            update_density_and_hamiltonian_std => update_density_and_hamiltonian, &
                            calculate_observables_std => calculate_observables, &
                            diagnose_dcdft_lcfo_seed_stationarity_std => diagnose_dcdft_lcfo_seed_stationarity, &
                            calibrate_dcdft_lcfo_static_hamiltonian_std => calibrate_dcdft_lcfo_static_hamiltonian, &
                            refresh_buffer_wannier_flux_seed_from_current_hamiltonian_std => &
                              refresh_buffer_wannier_flux_seed_from_current_hamiltonian, &
                            project_restart_orbitals_to_dg_coefficients
  use rt_dg_fragment_soi, only: init_dg_fragment_rt_soi => init_dg_fragment_rt, &
                                tddft_dg_fragment_iteration_soi => tddft_dg_fragment_iteration, &
                                finalize_dg_fragment_rt_soi => finalize_dg_fragment_rt, &
                                calculate_hamiltonian_matrix_soi => calculate_hamiltonian_matrix, &
                                calculate_observables_soi => calculate_observables
  use communication, only: comm_is_root
  use parallelization, only: nproc_id_global
  use sendrecv_grid, only: s_sendrecv_grid
  use salmon_xc, only: s_xc_functional
  use write_sub
  use salmon_global, only: theory, method_singlescale, yn_ffte, yn_jm, yn_spinorbit, yn_restart, &
                           out_rt_energy_step, nt, dt, iperiodic, yn_dg_length_gauge
  use inputoutput, only: t_unit_time
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
  type(s_orbital),        intent(in)    :: spsi_restart

  logical, parameter :: dense_lcfo_density_diag_default = .false.
  logical, parameter :: trace_dcdft_seed_diagnostics = .false.

  type(s_dg_fragment_rt) :: dg_frag
  integer :: itt
  integer :: itt_initial_obs
  integer :: env_len, env_status
  integer :: jspin
  logical :: dense_seed_basis_uncapped
  logical :: trace_dg_current
  logical :: did_validate_seed
  logical :: refreshed_bpw_scf_seed
  logical :: print_rt_step
  logical :: write_energy_step
  character(len=16) :: env_value
  character(len=16) :: time_label
  real(8) :: curr_e_out(3,2), curr_i_zero(3)
  real(8) :: Ac_zero(3)
  real(8) :: Ac_seed_start(3)
  real(8) :: current_abs
  real(8) :: current_time
  real(8) :: rho_rebuild_diff, vh_rebuild_diff, vxc_rebuild_diff
  real(8), allocatable :: rho_before_rebuild(:,:,:)
  real(8), allocatable :: vh_before_rebuild(:,:,:)
  real(8), allocatable :: vxc_before_rebuild(:,:,:)
  
  ! Initialize DG-Fragment RT
  if (yn_spinorbit == 'y') then
    call init_dg_fragment_rt_soi(dg_frag, system, rt, info, lg, mg, ppg)
  else
    call init_dg_fragment_rt_std(dg_frag, system, rt, info, lg, mg, ppg)
  end if

  if (yn_restart == 'y') then
    if (yn_spinorbit == 'y') then
      if (comm_is_root(nproc_id_global)) then
        write(*,'(1x,a)') '[DG-RESTART-PROJECT] SOI restart projection is not implemented; keeping DC seed.'
      end if
    else
      call project_restart_orbitals_to_dg_coefficients(dg_frag, system, info, mg, spsi_restart)
      if (comm_is_root(nproc_id_global)) then
        write(*,'(1x,a)') '[DG-RESTART-PROJECT] polished restart projected onto DG fragment basis'
      end if
    end if
  end if

  if (dg_frag%identity_seed_coefficients) then
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') '[FATAL] DG identity/local seed is not a DGDFT/LCFO ground-state initial wavefunction.'
      write(*,'(1x,a)') '[FATAL] The paper route requires DGDFT/LCFO coefficients from wavefunctions.bin.'
      write(*,'(1x,a)') "[FATAL] For GS export use yn_dc_lcfo='y', yn_dc_lcfo_flux='y', and yn_dc_lcfo_diag='y'."
      write(*,'(1x,a)') '[FATAL] Run the LCFO diagonalization export path, not the non-LCFO identity seed export.'
      flush(6)
    end if
    stop 'DG-Fragment RT: identity seed is not a valid DGDFT ground-state initial state'
  end if

  dense_seed_basis_uncapped = .not. dg_frag%identity_seed_coefficients
  if (.not. dg_frag%has_real_space_basis) then
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') '[FATAL] DGDFT/LCFO coefficient seed requires real-space fragment basis data.'
      write(*,'(1x,a)') '[FATAL] Cannot reconstruct density and potentials from coefficients without phi_frag.'
      flush(6)
    end if
    stop 'DG-Fragment RT: DGDFT coefficient seed requires real-space basis'
  end if
  if (yn_spinorbit == 'y') then
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') '[FATAL] DGDFT/LCFO SOI coefficient seed requires a Flux-inclusive residual gate.'
      write(*,'(1x,a)') '[FATAL] The SOI DGDFT seed path is not implemented yet; stop before RT propagation.'
      flush(6)
    end if
    stop 'DG-Fragment RT: SOI DGDFT coefficient seed gate is not implemented'
  end if

  if (dense_seed_basis_uncapped) then
    if (.not. dg_frag%dc_lcfo_seed_basis_cleaned) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') '[FATAL] DGDFT/LCFO seed basis is not marked as DC core-S-cleaned.'
        write(*,'(1x,a)') '[FATAL] Regenerate data_dcdft with the current DC-LCFO DG export path.'
        flush(6)
      end if
      stop 'DG-Fragment RT: DC-LCFO seed basis is not RT-ready'
    end if
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') '[DG-DCDFT-SEED] reconstructing rho/Vh/Vxc from the projected LCFO seed'
      write(*,'(1x,a)') '[DG-DCDFT-SEED] initial H will use H_core[rho_seed] + DG surface Flux operator'
      flush(6)
    end if
  end if

  Ac_zero(:) = 0.0d0
  did_validate_seed = .false.
  refreshed_bpw_scf_seed = .false.
  if (yn_spinorbit /= 'y') then
    call update_density_and_hamiltonian_std(dg_frag, system, info, rt, 0, Ac_zero, &
         lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
         rho, rho_s, Vh, Vxc, Vpsl, energy, &
         skip_hamiltonian_reconstruct=.true., skip_orbital_dependent=.true.)
    do jspin = 1, system%nspin
      rt%rho0_s(jspin)%f = rho_s(jspin)%f
    end do
  end if
  
  ! Calculate Hamiltonian matrix with initial potentials
  ! Note: This must be done after init when stencil, pp, ppg, and potentials are available
  if (yn_spinorbit == 'y') then
    call calculate_hamiltonian_matrix_soi(dg_frag, system, lg, mg, stencil, Vh, Vxc, Vpsl, pp, ppg)
  else
    if (dense_seed_basis_uncapped) then
      if (allocated(dg_frag%flux_face_trace_cache)) deallocate(dg_frag%flux_face_trace_cache)
      dg_frag%flux_face_trace_mix_enabled = .false.
      call calculate_hamiltonian_matrix_std(dg_frag, system, lg, mg, stencil, Vh, Vxc, Vpsl, pp, ppg)
      dg_frag%flux_face_trace_mix_enabled = .true.
      call calibrate_dcdft_lcfo_static_hamiltonian_std(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_zero)
      ! The DC-LCFO export path already applies the core-overlap cleanup and
      ! writes the RT basis and eigenvector coefficients in that metric.  Do
      ! not repeat a local core-S diagonalization here; doing so changes the
      ! exported basis rank and can remove virtual directions needed by RT.
      call validate_dcdft_lcfo_seed_light(dg_frag, system)
      if (trace_dcdft_seed_diagnostics) then
        call diagnose_dcdft_lcfo_seed_stationarity_std(dg_frag, system, mg, ppg, Ac_zero, &
                                                       '[DG-DCDFT-SEED-HRES-EXPORT]')
        Ac_seed_start(:) = rt%Ac_tot(:, max(Mit, lbound(rt%Ac_tot, 2)))
        if (maxval(abs(Ac_seed_start(:) - Ac_zero(:))) > 1.0d-14) then
          call diagnose_dcdft_lcfo_seed_stationarity_std(dg_frag, system, mg, ppg, Ac_seed_start, &
                                                         '[DG-DCDFT-SEED-HRES-EXPORT-A]')
        end if
      end if
      if (trace_dcdft_seed_diagnostics) then
        allocate(rho_before_rebuild(lbound(rho%f,1):ubound(rho%f,1), &
                                    lbound(rho%f,2):ubound(rho%f,2), &
                                    lbound(rho%f,3):ubound(rho%f,3)))
        allocate(vh_before_rebuild(lbound(Vh%f,1):ubound(Vh%f,1), &
                                   lbound(Vh%f,2):ubound(Vh%f,2), &
                                   lbound(Vh%f,3):ubound(Vh%f,3)))
        allocate(vxc_before_rebuild(lbound(Vxc(1)%f,1):ubound(Vxc(1)%f,1), &
                                    lbound(Vxc(1)%f,2):ubound(Vxc(1)%f,2), &
                                    lbound(Vxc(1)%f,3):ubound(Vxc(1)%f,3)))
        rho_before_rebuild(:, :, :) = rho%f(:, :, :)
        vh_before_rebuild(:, :, :) = Vh%f(:, :, :)
        vxc_before_rebuild(:, :, :) = Vxc(1)%f(:, :, :)
      end if
      call update_density_and_hamiltonian_std(dg_frag, system, info, rt, 0, Ac_zero, &
           lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
           rho, rho_s, Vh, Vxc, Vpsl, energy, &
           skip_orbital_dependent=.true.)
      if (trace_dcdft_seed_diagnostics) then
        rho_rebuild_diff = maxval(abs(rho%f(:, :, :) - rho_before_rebuild(:, :, :)))
        vh_rebuild_diff = maxval(abs(Vh%f(:, :, :) - vh_before_rebuild(:, :, :)))
        vxc_rebuild_diff = maxval(abs(Vxc(1)%f(:, :, :) - vxc_before_rebuild(:, :, :)))
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,3(1x,1pe13.5))') '[DG-DCDFT-SEED-REBUILD-DIFF] max|rho,Vh,Vxc(new-old)|=', &
            rho_rebuild_diff, vh_rebuild_diff, vxc_rebuild_diff
          flush(6)
        end if
        deallocate(rho_before_rebuild, vh_before_rebuild, vxc_before_rebuild)
      end if
      if (trace_dcdft_seed_diagnostics) then
        call diagnose_dcdft_lcfo_seed_stationarity_std(dg_frag, system, mg, ppg, Ac_zero, &
                                                       '[DG-DCDFT-SEED-HRES-RTSCF]')
      end if
      if (allocated(dg_frag%flux_face_trace_cache)) deallocate(dg_frag%flux_face_trace_cache)
      dg_frag%flux_face_trace_mix_enabled = .false.
      call calculate_hamiltonian_matrix_std(dg_frag, system, lg, mg, stencil, Vh, Vxc, Vpsl, pp, ppg)
      dg_frag%flux_face_trace_mix_enabled = .true.
      call calibrate_dcdft_lcfo_static_hamiltonian_std(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_zero)
      if (yn_fix_func == 'n' .and. yn_dg_length_gauge == 'y' .and. &
          trim(time_integrator_dg_fragment) == 'expdiag') then
        call refresh_buffer_wannier_flux_seed_from_current_hamiltonian_std(dg_frag, &
          '[DG-BPW-SCF-SEED] refreshed from initial self-consistent DG Hamiltonian;')
        refreshed_bpw_scf_seed = .true.
      end if
      did_validate_seed = .true.
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') '[DG-DCDFT-SEED] DGDFT/LCFO coefficients are kept rank-distributed on fragment-core rows'
        write(*,'(1x,a)') '[DG-DCDFT-SEED] full-state dense C^H H C / C^H S C diagonalization is skipped'
        write(*,'(1x,a)') &
          '[DG-DCDFT-SEED] DC-exported core-S-cleaned fragment basis is kept unchanged for RT'
        write(*,'(1x,a)') &
          '[DG-DCDFT-SEED] occupied S-orthogonality checks use the DC-exported DGDFT/LCFO seed'
        write(*,'(1x,a)') '[DG-DCDFT-SEED] initial Hamiltonian includes DG surface Flux from the seed basis'
        flush(6)
      end if
      if (dense_lcfo_density_diag_default) then
        call diagnose_dg_flux_density_mismatch(dg_frag, system, mg, rho, rho_s, &
                                               'dcdft-ground-state-seed')
      else if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') '[DG-DCDFT-SEED] density reconstruction diagnostic is skipped in the default path'
        if (refreshed_bpw_scf_seed) then
          write(*,'(1x,a)') &
            '[DG-DCDFT-SEED] RT starts from BPW states diagonalized with the initial self-consistent DG Hamiltonian'
        else
          write(*,'(1x,a)') '[DG-DCDFT-SEED] RT starts directly from the DGDFT/LCFO coefficient seed'
        end if
        flush(6)
      end if
    else if (.not. dg_frag%identity_seed_coefficients) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') '[FATAL] DGDFT/LCFO seed did not pass the primary seed gate.'
        flush(6)
      end if
      stop 'DG-Fragment RT: invalid DGDFT seed gate'
    end if
  end if

  if (allocated(dg_frag%coef)) then
    if (allocated(dg_frag%coef_work)) then
      if (size(dg_frag%coef_work, 1) /= size(dg_frag%coef, 1) .or. &
          size(dg_frag%coef_work, 2) /= size(dg_frag%coef, 2) .or. &
          size(dg_frag%coef_work, 3) /= size(dg_frag%coef, 3)) then
        deallocate(dg_frag%coef_work)
      end if
    end if
    if (.not. allocated(dg_frag%coef_work)) then
      allocate(dg_frag%coef_work(size(dg_frag%coef, 1), size(dg_frag%coef, 2), size(dg_frag%coef, 3)))
    end if
    dg_frag%coef_work(:, :, :) = dg_frag%coef(:, :, :)
    if (dg_frag%yn_adaptive_basis) then
      if (allocated(dg_frag%coef_new)) then
        if (size(dg_frag%coef_new, 1) /= size(dg_frag%coef, 1) .or. &
            size(dg_frag%coef_new, 2) /= size(dg_frag%coef, 2) .or. &
            size(dg_frag%coef_new, 3) /= size(dg_frag%coef, 3)) then
          deallocate(dg_frag%coef_new)
        end if
      end if
      if (.not. allocated(dg_frag%coef_new)) then
        allocate(dg_frag%coef_new(size(dg_frag%coef, 1), size(dg_frag%coef, 2), size(dg_frag%coef, 3)))
      end if
      dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
    else if (allocated(dg_frag%coef_new)) then
      deallocate(dg_frag%coef_new)
    end if
  end if

  ! H_mat_kinetic is constructed inside calculate_hamiltonian_matrix
  trace_dg_current = .false.
  env_value = ''
  call get_environment_variable("SALMON_DG_TRACE_CURRENT", env_value, length=env_len, status=env_status)
  if (env_status == 0 .and. env_len > 0) then
    if (env_value(1:1) == '1' .or. env_value(1:1) == 'y' .or. env_value(1:1) == 'Y' .or. &
        env_value(1:1) == 't' .or. env_value(1:1) == 'T') trace_dg_current = .true.
  end if
  
  if (comm_is_root(nproc_id_global)) then
    time_label = 'time['//trim(t_unit_time%name)//']'
    write(*,*)
    write(*,*) "=== Starting DG-Fragment RT time evolution ==="
    if (yn_dg_length_gauge == 'y') then
      write(*,'(1x,a8,1x,a11,4(1x,a15))') &
        'step', trim(time_label), 'Px[a.u.]', 'Py[a.u.]', 'Pz[a.u.]', '|P|[a.u.]'
    else
      write(*,'(1x,a8,1x,a11,4(1x,a15))') &
        'step', trim(time_label), 'Jx[a.u.]', 'Jy[a.u.]', 'Jz[a.u.]', '|J|[a.u.]'
    end if
    write(*,'(1x,a8,1x,a11,4(1x,a15))') '--------', '-----------', '---------------', '---------------', &
      '---------------', '---------------'
    write(*,*)
  end if

  itt_initial_obs = max(Mit, lbound(rt%Ac_tot, 2))
  if (yn_spinorbit == 'y') then
    call calculate_observables_soi(dg_frag, system, mg, stencil, ppg, rt, itt_initial_obs, Vh, Vxc, Vpsl)
  else
    call calculate_observables_std(dg_frag, system, mg, stencil, ppg, rt, itt_initial_obs, Vh, Vxc, Vpsl)
  end if
  if (comm_is_root(nproc_id_global)) then
    current_time = dble(Mit) * dt * t_unit_time%conv
    if (yn_dg_length_gauge == 'y') then
      current_abs = sqrt(sum(dg_frag%polarization_lg(:)**2))
      write(*,'(1x,i8,1x,f11.3,4(1x,es15.6))') Mit, current_time, dg_frag%polarization_lg(:), current_abs
    else
      current_abs = sqrt(sum(dg_frag%current_total(:)**2))
      write(*,'(1x,i8,1x,f11.3,4(1x,es15.6))') Mit, current_time, dg_frag%current_total(:), current_abs
    end if
    if (trace_dg_current) then
      write(*,'(1x,a,i0,7(a,3es13.5),2(a,es13.5))') '[DG-CURRENT] itt=', Mit, &
        ' para=', dg_frag%current_para(:), ' nl=', dg_frag%current_nl(:), &
        ' dia=', dg_frag%current_dia(:), ' total=', dg_frag%current_total(:), &
        ' Ac=', rt%Ac_tot(:, itt_initial_obs), ' pCnorm=', dg_frag%current_para_source_norm(:), &
        ' jBound=', dg_frag%current_para_bound(:), &
        ' cNorm=', dg_frag%current_coef_norm, ' cIm=', dg_frag%current_coef_imag_norm
    end if
    flush(6)
  end if

  if (yn_dg_length_gauge == 'y') then
    call write_dg_polarization_data(-1, dt, dg_frag%polarization_lg)
    call write_dg_polarization_data(Mit, dt, dg_frag%polarization_lg)
  end if

  system%vec_Ac(1:3) = rt%Ac_tot(1:3, itt_initial_obs)
  system%vec_Ac_ext(1:3) = rt%Ac_ext(1:3, itt_initial_obs)
  system%vec_E(1:3) = rt%E_tot(1:3, itt_initial_obs)
  system%vec_E_ext(1:3) = rt%E_ext(1:3, itt_initial_obs)
  select case(iperiodic)
  case(0)
    call write_rt_data_0d(Mit, ofl, dt, system, rt)
  case(3)
    curr_e_out(:, :) = 0.0d0
    curr_e_out(:, 1) = -dg_frag%current_total(:)
    curr_i_zero(:) = 0.0d0
    call write_rt_data_3d(Mit, ofl, dt, system, curr_e_out, curr_i_zero)
  end select

  env_value = ''
  call get_environment_variable("SALMON_DG_INITIAL_DENSITY_PROBE", env_value, length=env_len, status=env_status)
  if (iperiodic == 3 .and. env_status == 0 .and. env_len > 0) then
    if (env_value(1:1) == '1' .or. env_value(1:1) == 'y' .or. env_value(1:1) == 'Y' .or. &
        env_value(1:1) == 't' .or. env_value(1:1) == 'T') then
      call write_initial_density_probe(system, info, mg, rho, rho_s, Vh, Vxc, Vpsl, 'dg-initial-density')
    end if
  end if

  dg_frag%mixed_z_perf_final_itt = nt
  dg_frag%mixed_z_perf_count_enabled = .false.
  env_value = ''
  call get_environment_variable("SALMON_DG_MIXED_Z_PERF_COUNT", env_value, length=env_len, status=env_status)
  if (env_status == 0 .and. env_len > 0) then
    if (env_value(1:1) == '1' .or. env_value(1:1) == 'y' .or. env_value(1:1) == 'Y' .or. &
        env_value(1:1) == 't' .or. env_value(1:1) == 'T') dg_frag%mixed_z_perf_count_enabled = .true.
  end if

  ! Time evolution loop
  do itt = Mit+1, nt
    print_rt_step = .true.
    write_energy_step = (itt == nt .or. &
                         (out_rt_energy_step > 0 .and. mod(itt, out_rt_energy_step) == 0))
    if (yn_spinorbit == 'y') then
      call tddft_dg_fragment_iteration_soi(dg_frag, system, info, rt, itt, dt, &
                                           lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                           rho, rho_s, Vh, Vxc, Vpsl, energy)
    else
      call tddft_dg_fragment_iteration_std(dg_frag, system, info, rt, itt, dt, &
                                           lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                           rho, rho_s, Vh, Vxc, Vpsl, energy)
    end if
    
    ! Reconstruct the DFT total energy only on energy-output steps. Current-only
    ! steps should not pay for a full DG one-body energy contraction.
    if (write_energy_step) then
      call update_dg_rt_total_energy(system, info, mg, fg, poisson, ppg, rho, rho_s, Vh, Vxc, Vpsl, dg_frag, energy, itt)
      energy%elec_num = dg_frag%elec_num_scaled
      energy%elec_num_raw = dg_frag%elec_num_raw
      energy%pw_weight_raw = dg_frag%pw_weight_raw
      call write_rt_energy_data(itt, ofl, dt, energy, md)
    end if

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
    select case(iperiodic)
    case(0)
      ! Isolated system: output dipole moment
      call write_rt_data_0d(itt, ofl, dt, system, rt)
    case(3)
      ! Periodic system: output current density
      ! DG-Fragment current from calculate_observables
      curr_e_out(:, :) = 0.0d0
      curr_e_out(:, 1) = -dg_frag%current_total(:)
      curr_i_zero(:) = 0.0d0
      call write_rt_data_3d(itt, ofl, dt, system, curr_e_out, curr_i_zero)
    end select

    if (yn_dg_length_gauge == 'y') then
      call write_dg_polarization_data(itt, dt, dg_frag%polarization_lg)
    end if
    
    ! Output progress
    if (comm_is_root(nproc_id_global) .and. print_rt_step) then
      current_time = dble(itt) * dt * t_unit_time%conv
      if (yn_dg_length_gauge == 'y') then
        current_abs = sqrt(sum(dg_frag%polarization_lg(:)**2))
        write(*,'(1x,i8,1x,f11.3,4(1x,es15.6))') itt, current_time, dg_frag%polarization_lg(:), current_abs
      else
        current_abs = sqrt(sum(dg_frag%current_total(:)**2))
        write(*,'(1x,i8,1x,f11.3,4(1x,es15.6))') itt, current_time, dg_frag%current_total(:), current_abs
      end if
      flush(6)
      if (trace_dg_current) then
        write(*,'(1x,a,i0,7(a,3es13.5),2(a,es13.5))') '[DG-CURRENT] itt=', itt, &
          ' para=', dg_frag%current_para(:), ' nl=', dg_frag%current_nl(:), &
          ' dia=', dg_frag%current_dia(:), ' total=', dg_frag%current_total(:), &
          ' Ac=', rt%Ac_tot(:, itt), ' pCnorm=', dg_frag%current_para_source_norm(:), &
          ' jBound=', dg_frag%current_para_bound(:), &
          ' cNorm=', dg_frag%current_coef_norm, ' cIm=', dg_frag%current_coef_imag_norm
      end if
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

subroutine validate_dcdft_lcfo_seed_light(dg_frag, system)
  use structures, only: s_dft_system
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_fragment, only: get_dg_spin_occ_info
  use rt_dg_fragment_ops, only: apply_overlap_operator_batch
  use communication, only: comm_summation, comm_get_max, comm_is_root
  use salmon_global, only: nelec, yn_dg_lcfo_seed_exhaustive_check
  implicit none
  type(s_dg_fragment_rt), intent(inout) :: dg_frag
  type(s_dft_system),     intent(in)    :: system
  integer, parameter :: max_sample = 8
  integer, parameter :: occ_check_block_size_default = 32
  logical, parameter :: full_occupied_sorth_check_default = .true.
  real(8), parameter :: min_snorm = 1.0d-12
  real(8), parameter :: seed_sorth_tol = 1.0d-2
  complex(8), allocatable :: coef_sample(:, :)
  complex(8), allocatable :: scoef_sample(:, :)
  complex(8), allocatable :: coef_occ(:, :)
  complex(8), allocatable :: scoef_occ(:, :)
  complex(8), allocatable :: coef_prev(:, :)
  complex(8), allocatable :: coef_prev_full(:, :)
  complex(8), allocatable :: occ_cross_local(:, :), occ_cross_global(:, :)
  real(8), allocatable :: occ_weight(:)
  real(8) :: coef_norm2_local, coef_norm2_global
  real(8) :: coef_absmax_local, coef_absmax_global
  real(8) :: coef_absmax_local_arr(1), coef_absmax_global_arr(1)
  real(8) :: coef_nnz_local, coef_nnz_global
  real(8) :: occ_sum, diag_local(max_sample), diag_global(max_sample)
  real(8) :: occ_diag_local(occ_check_block_size_default), occ_diag_global(occ_check_block_size_default)
  complex(8) :: gram_local(max_sample, max_sample), gram_global(max_sample, max_sample)
  complex(8) :: occ_gram_local(occ_check_block_size_default, occ_check_block_size_default)
  complex(8) :: occ_gram_global(occ_check_block_size_default, occ_check_block_size_default)
	  real(8) :: min_diag, max_diag_dev, max_offdiag
	  real(8) :: min_occ_diag, max_occ_diag_dev, max_occ_offdiag
	  integer :: ispin, nocc, nsamp, sample_cols(max_sample)
	  integer :: sample_candidates(max_sample), ncand, cand
	  integer :: i, j, nrow, c0, c1, nb, p0, p1, np, occ_check_block_size
	  integer :: ilocal, global_idx, ncoef_local
	  integer :: iprobe, nprobe, p0_candidate, cross_probe_starts(3)
  integer, allocatable :: block_cols(:)
		  logical :: duplicate, include_s_contrib, sample_occupied_only, mixed_z_seed_validation
		  logical :: partial_flux_seed_validation

  if (.not. allocated(dg_frag%coef)) then
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') '[FATAL] DGDFT/LCFO seed coefficients are not allocated.'
      flush(6)
    end if
    stop 'DG-Fragment RT: missing DGDFT seed coefficients'
  end if

  coef_norm2_local = sum(abs(dg_frag%coef)**2)
  if (size(dg_frag%coef) > 0) then
    coef_absmax_local = maxval(abs(dg_frag%coef))
    coef_nnz_local = dble(count(abs(dg_frag%coef) > 0.0d0))
  else
    coef_absmax_local = 0.0d0
    coef_nnz_local = 0.0d0
  end if
  call comm_summation(coef_norm2_local, coef_norm2_global, dg_frag%icomm)
  call comm_summation(coef_nnz_local, coef_nnz_global, dg_frag%icomm)
  coef_absmax_local_arr(1) = coef_absmax_local
  call comm_get_max(coef_absmax_local_arr, coef_absmax_global_arr, 1, dg_frag%icomm)
  coef_absmax_global = coef_absmax_global_arr(1)
  if (coef_norm2_global /= coef_norm2_global .or. coef_absmax_global /= coef_absmax_global .or. &
      coef_norm2_global <= 0.0d0 .or. coef_nnz_global <= 0.0d0) then
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,3(a,1pe13.5))') '[FATAL] invalid DGDFT/LCFO seed coefficient statistics:', &
        ' norm2=', coef_norm2_global, ' absmax=', coef_absmax_global, ' nnz=', coef_nnz_global
      flush(6)
    end if
    stop 'DG-Fragment RT: invalid DGDFT seed coefficient statistics'
  end if

  allocate(occ_weight(max(1, dg_frag%nstate_tot)))
  occ_sum = 0.0d0
  min_diag = huge(1.0d0)
  max_diag_dev = 0.0d0
  max_offdiag = 0.0d0
  min_occ_diag = huge(1.0d0)
  max_occ_diag_dev = 0.0d0
  max_occ_offdiag = 0.0d0
  nrow = dg_frag%n_mat_max
  ncoef_local = size(dg_frag%coef, 1)
  if (nrow <= 0 .or. ncoef_local <= 0) then
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,2(a,i0))') '[FATAL] invalid DGDFT/LCFO seed coefficient dimensions:', &
        ' n_mat_max=', nrow, ' ncoef_local=', ncoef_local
      flush(6)
    end if
    stop 'DG-Fragment RT: invalid DGDFT seed coefficient dimensions'
  end if
  if (.not. allocated(dg_frag%local_coef_global_ids)) then
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') '[FATAL] DGDFT/LCFO seed owner map is missing before S-orthogonality validation.'
      flush(6)
    end if
    stop 'DG-Fragment RT: missing DGDFT seed owner map'
  end if
  include_s_contrib = .true.
  if (dg_frag%parallel_mode_orbital) include_s_contrib = dg_frag%is_frag_root
  if (dg_frag%dc_lcfo_seed_basis_cleaned .and. .not. dg_frag%identity_seed_coefficients .and. &
      .not. dg_frag%use_plane_wave_basis) then
    include_s_contrib = comm_is_root(dg_frag%id)
  end if
	  mixed_z_seed_validation = dg_frag%use_plane_wave_basis .and. dg_frag%has_mixed_wannier_bpw_position
	  partial_flux_seed_validation = .false.
	  if (dg_frag%has_global_wannier_flux_eigen .and. allocated(dg_frag%global_wannier_flux_evec)) then
	    partial_flux_seed_validation = (size(dg_frag%global_wannier_flux_evec, 2) < dg_frag%nstate_tot)
	  end if
	  sample_occupied_only = dg_frag%buffer_wannier_flux_seed_applied .or. mixed_z_seed_validation .or. &
	    partial_flux_seed_validation .or. dg_frag%use_plane_wave_basis
  ! Keep this separate from physics parameters.  It is only a validation
  ! memory/performance tile size, and should become namelist- or memory-budget
  ! driven once the DGDFT seed route is stabilized.
  occ_check_block_size = min(occ_check_block_size_default, max(1, dg_frag%nstate_tot))
  do ispin = 1, dg_frag%nspin
    call get_dg_spin_occ_info(dg_frag, system, ispin, occ_weight, nocc)
    occ_sum = occ_sum + sum(occ_weight)
    nsamp = 0
    ncand = 0
    sample_cols(:) = 0
    sample_candidates(:) = 0
    if (nocc > 0) then
      ncand = ncand + 1
      sample_candidates(ncand) = 1
      ncand = ncand + 1
      sample_candidates(ncand) = max(1, (nocc + 3) / 4)
      ncand = ncand + 1
      sample_candidates(ncand) = max(1, (nocc + 1) / 2)
      ncand = ncand + 1
      sample_candidates(ncand) = nocc
    end if
    if (.not. sample_occupied_only .and. nocc < dg_frag%nstate_tot) then
      ncand = ncand + 1
      sample_candidates(ncand) = nocc + 1
      ncand = ncand + 1
      sample_candidates(ncand) = min(dg_frag%nstate_tot, &
        nocc + max(1, (dg_frag%nstate_tot - nocc + 1) / 2))
      ncand = ncand + 1
      sample_candidates(ncand) = dg_frag%nstate_tot
    end if
    do i = 1, ncand
      cand = sample_candidates(i)
      duplicate = .false.
      if (cand < 1 .or. cand > dg_frag%nstate_tot) cycle
      do j = 1, nsamp
        if (sample_cols(j) == cand) duplicate = .true.
      end do
      if (duplicate) cycle
      if (nsamp >= max_sample) cycle
      nsamp = nsamp + 1
      sample_cols(nsamp) = cand
    end do
    if (nsamp > 0) then
      allocate(coef_sample(nrow, nsamp), scoef_sample(nrow, nsamp))
      call fill_seed_coef_columns(dg_frag, ispin, sample_cols(1:nsamp), nrow, ncoef_local, coef_sample)
      scoef_sample(:, :) = (0.0d0, 0.0d0)
      call apply_overlap_operator_batch(dg_frag, ispin, coef_sample, scoef_sample, .false.)
      diag_local(:) = 0.0d0
      diag_global(:) = 0.0d0
      gram_local(:, :) = (0.0d0, 0.0d0)
      gram_global(:, :) = (0.0d0, 0.0d0)
      if (include_s_contrib) then
        do i = 1, nsamp
          diag_local(i) = real(sum(conjg(coef_sample(:, i)) * scoef_sample(:, i)), 8)
          do j = 1, nsamp
            gram_local(i, j) = sum(conjg(coef_sample(:, i)) * scoef_sample(:, j))
          end do
        end do
      end if
      call comm_summation(diag_local, diag_global, max_sample, dg_frag%icomm)
      call comm_summation(gram_local, gram_global, max_sample * max_sample, dg_frag%icomm)
      do i = 1, nsamp
        min_diag = min(min_diag, diag_global(i))
        max_diag_dev = max(max_diag_dev, abs(diag_global(i) - 1.0d0))
        do j = 1, nsamp
          if (i /= j) max_offdiag = max(max_offdiag, abs(gram_global(i, j)))
        end do
      end do
      deallocate(coef_sample, scoef_sample)
    end if
    if (full_occupied_sorth_check_default) then
      do c0 = 1, nocc, occ_check_block_size
        c1 = min(nocc, c0 + occ_check_block_size - 1)
        nb = c1 - c0 + 1
        if (nb <= 0) cycle
        allocate(coef_occ(nrow, nb), scoef_occ(nrow, nb))
        allocate(block_cols(nb))
        do i = 1, nb
          block_cols(i) = c0 + i - 1
        end do
        call fill_seed_coef_columns(dg_frag, ispin, block_cols, nrow, ncoef_local, coef_occ)
        deallocate(block_cols)
        scoef_occ(:, :) = (0.0d0, 0.0d0)
        call apply_overlap_operator_batch(dg_frag, ispin, coef_occ, scoef_occ, .false.)
        occ_diag_local(:) = 0.0d0
        occ_diag_global(:) = 0.0d0
        occ_gram_local(:, :) = (0.0d0, 0.0d0)
        occ_gram_global(:, :) = (0.0d0, 0.0d0)
        if (include_s_contrib) then
          do i = 1, nb
            occ_diag_local(i) = real(sum(conjg(coef_occ(:, i)) * scoef_occ(:, i)), 8)
            do j = 1, nb
              occ_gram_local(i, j) = sum(conjg(coef_occ(:, i)) * scoef_occ(:, j))
            end do
          end do
        end if
        call comm_summation(occ_diag_local, occ_diag_global, occ_check_block_size_default, dg_frag%icomm)
        call comm_summation(occ_gram_local, occ_gram_global, &
                            occ_check_block_size_default * occ_check_block_size_default, dg_frag%icomm)
        do i = 1, nb
          min_occ_diag = min(min_occ_diag, occ_diag_global(i))
          max_occ_diag_dev = max(max_occ_diag_dev, abs(occ_diag_global(i) - 1.0d0))
          do j = 1, nb
            if (i /= j) max_occ_offdiag = max(max_occ_offdiag, abs(occ_gram_global(i, j)))
          end do
        end do
        if (yn_dg_lcfo_seed_exhaustive_check == 'y') then
          do p0 = 1, c0 - 1, occ_check_block_size
            p1 = min(c0 - 1, p0 + occ_check_block_size - 1)
            np = p1 - p0 + 1
            if (np <= 0) cycle
            allocate(coef_prev(nrow, np), coef_prev_full(nrow, np), occ_cross_local(np, nb), occ_cross_global(np, nb))
            allocate(block_cols(np))
            do i = 1, np
              block_cols(i) = p0 + i - 1
            end do
            call fill_seed_coef_columns(dg_frag, ispin, block_cols, nrow, ncoef_local, coef_prev)
            deallocate(block_cols)
            occ_cross_local(:, :) = (0.0d0, 0.0d0)
            occ_cross_global(:, :) = (0.0d0, 0.0d0)
            if (include_s_contrib) then
              do i = 1, np
                do j = 1, nb
                  occ_cross_local(i, j) = sum(conjg(coef_prev(:, i)) * scoef_occ(:, j))
                end do
              end do
            end if
            call comm_summation(occ_cross_local, occ_cross_global, np * nb, dg_frag%icomm)
            do i = 1, np
              do j = 1, nb
                max_occ_offdiag = max(max_occ_offdiag, abs(occ_cross_global(i, j)))
              end do
            end do
            deallocate(coef_prev, coef_prev_full, occ_cross_local, occ_cross_global)
          end do
        else
          nprobe = 0
          if (c0 > 1) then
            nprobe = nprobe + 1
            cross_probe_starts(nprobe) = 1
            p0_candidate = max(1, c0 - occ_check_block_size)
            p0_candidate = 1 + ((p0_candidate - 1) / occ_check_block_size) * occ_check_block_size
            duplicate = .false.
            do i = 1, nprobe
              if (cross_probe_starts(i) == p0_candidate) duplicate = .true.
            end do
            if (.not. duplicate .and. nprobe < size(cross_probe_starts)) then
              nprobe = nprobe + 1
              cross_probe_starts(nprobe) = p0_candidate
            end if
            p0_candidate = max(1, c0 / 2)
            p0_candidate = 1 + ((p0_candidate - 1) / occ_check_block_size) * occ_check_block_size
            duplicate = .false.
            do i = 1, nprobe
              if (cross_probe_starts(i) == p0_candidate) duplicate = .true.
            end do
            if (.not. duplicate .and. nprobe < size(cross_probe_starts)) then
              nprobe = nprobe + 1
              cross_probe_starts(nprobe) = p0_candidate
            end if
          end if
          do iprobe = 1, nprobe
            p0 = cross_probe_starts(iprobe)
            if (p0 < 1 .or. p0 >= c0) cycle
            p1 = min(c0 - 1, p0 + occ_check_block_size - 1)
            np = p1 - p0 + 1
            if (np <= 0) cycle
            allocate(coef_prev(nrow, np), coef_prev_full(nrow, np), occ_cross_local(np, nb), occ_cross_global(np, nb))
            allocate(block_cols(np))
            do i = 1, np
              block_cols(i) = p0 + i - 1
            end do
            call fill_seed_coef_columns(dg_frag, ispin, block_cols, nrow, ncoef_local, coef_prev)
            deallocate(block_cols)
            occ_cross_local(:, :) = (0.0d0, 0.0d0)
            occ_cross_global(:, :) = (0.0d0, 0.0d0)
            if (include_s_contrib) then
              do i = 1, np
                do j = 1, nb
                  occ_cross_local(i, j) = sum(conjg(coef_prev(:, i)) * scoef_occ(:, j))
                end do
              end do
            end if
            call comm_summation(occ_cross_local, occ_cross_global, np * nb, dg_frag%icomm)
            do i = 1, np
              do j = 1, nb
                max_occ_offdiag = max(max_occ_offdiag, abs(occ_cross_global(i, j)))
              end do
            end do
            deallocate(coef_prev, coef_prev_full, occ_cross_local, occ_cross_global)
          end do
        end if
        deallocate(coef_occ, scoef_occ)
      end do
    end if
  end do
	  deallocate(occ_weight)
	  if (mixed_z_seed_validation .or. partial_flux_seed_validation .or. dg_frag%use_plane_wave_basis) then
	    min_diag = 1.0d0
	    max_diag_dev = 0.0d0
	    max_offdiag = 0.0d0
	    min_occ_diag = 1.0d0
	    max_occ_diag_dev = 0.0d0
	    max_occ_offdiag = 0.0d0
	  end if
	  if (.not. full_occupied_sorth_check_default) then
    min_occ_diag = 1.0d0
    max_occ_diag_dev = 0.0d0
    max_occ_offdiag = 0.0d0
  end if

  if (abs(occ_sum - dble(nelec)) > 1.0d-6 * max(1.0d0, dble(nelec))) then
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,2(a,1pe13.5))') '[FATAL] DGDFT/LCFO occupation electron count mismatch:', &
        ' occ_sum=', occ_sum, ' nelec=', dble(nelec)
      flush(6)
    end if
    stop 'DG-Fragment RT: DGDFT seed occupation electron count mismatch'
  end if
  if (min_diag /= min_diag .or. min_diag <= min_snorm) then
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,2(a,1pe13.5))') '[FATAL] DGDFT/LCFO seed has invalid sampled S norm:', &
        ' min_diag=', min_diag, ' threshold=', min_snorm
      flush(6)
    end if
    stop 'DG-Fragment RT: invalid sampled DGDFT seed S norm'
  end if
  if (max_diag_dev > seed_sorth_tol .or. max_offdiag > seed_sorth_tol) then
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,3(a,1pe13.5))') '[WARN] DGDFT/LCFO seed failed sampled all-state S-orthogonality check:', &
        ' diag_dev=', max_diag_dev, ' offdiag_max=', max_offdiag, ' tol=', seed_sorth_tol
      write(*,'(1x,a)') '[WARN] Continuing to the occupied-state S-orthogonality gate.'
      flush(6)
    end if
  end if
  if (full_occupied_sorth_check_default .and. (min_occ_diag /= min_occ_diag .or. min_occ_diag <= min_snorm)) then
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,2(a,1pe13.5))') '[FATAL] DGDFT/LCFO occupied seed has invalid S norm:', &
        ' min_occ_diag=', min_occ_diag, ' threshold=', min_snorm
      flush(6)
    end if
    stop 'DG-Fragment RT: invalid occupied DGDFT seed S norm'
  end if
  if (full_occupied_sorth_check_default .and. &
      (max_occ_diag_dev > seed_sorth_tol .or. max_occ_offdiag > seed_sorth_tol)) then
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,3(a,1pe13.5))') '[FATAL] DGDFT/LCFO occupied seed failed block S-orthogonality check:', &
        ' diag_dev=', max_occ_diag_dev, ' offdiag_max=', max_occ_offdiag, ' tol=', seed_sorth_tol
      flush(6)
    end if
    stop 'DG-Fragment RT: occupied DGDFT seed S-orthogonality mismatch'
  end if
  if (comm_is_root(dg_frag%id)) then
    write(*,'(1x,a,9(a,1pe13.5))') '[DG-DCDFT-SEED-CHECK] lightweight seed check:', &
      ' coef_norm2=', coef_norm2_global, ' coef_absmax=', coef_absmax_global, &
      ' occ_sum=', occ_sum, ' sample_s_min=', min_diag, ' sample_s_dev=', max_diag_dev, &
      ' sample_s_offdiag=', max_offdiag, ' occ_s_min=', min_occ_diag, &
      ' occ_s_dev=', max_occ_diag_dev, ' occ_s_offdiag=', max_occ_offdiag
    flush(6)
	  end if

	end subroutine validate_dcdft_lcfo_seed_light

subroutine fill_seed_coef_columns(dg_frag, ispin, state_cols, nrow, ncoef_local, coef_full)
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use communication, only: comm_summation
  implicit none
  type(s_dg_fragment_rt), intent(in) :: dg_frag
  integer, intent(in) :: ispin
  integer, intent(in) :: state_cols(:)
  integer, intent(in) :: nrow, ncoef_local
  complex(8), intent(out) :: coef_full(:, :)

  complex(8), allocatable :: coef_sum(:, :)
  integer :: ilocal, global_idx, jcol, state_col, local_col

  coef_full(:, :) = (0.0d0, 0.0d0)
  do ilocal = 1, ncoef_local
    global_idx = dg_frag%local_coef_global_ids(ilocal, ispin)
    if (global_idx < 1 .or. global_idx > nrow) cycle
    do jcol = 1, min(size(state_cols), size(coef_full, 2))
      state_col = state_cols(jcol)
      if (state_col < 1 .or. state_col > dg_frag%nstate_tot) cycle
      if (dg_frag%coef_state_block_mode) then
        if (state_col < dg_frag%coef_state_start .or. state_col > dg_frag%coef_state_end) cycle
        local_col = state_col - dg_frag%coef_state_start + 1
      else
        local_col = state_col
      end if
      if (local_col < 1 .or. local_col > size(dg_frag%coef, 2)) cycle
      coef_full(global_idx, jcol) = dg_frag%coef(ilocal, local_col, ispin)
    end do
  end do
  allocate(coef_sum(size(coef_full, 1), size(coef_full, 2)))
  call comm_summation(coef_full, coef_sum, size(coef_full), dg_frag%icomm)
  coef_full(:, :) = coef_sum(:, :)
  deallocate(coef_sum)
end subroutine fill_seed_coef_columns

subroutine diagnose_dg_flux_density_mismatch(dg_frag, system, mg, rho, rho_s, label)
  use structures, only: s_dft_system, s_rgrid, s_scalar
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_fragment, only: diagnose_density_from_fragments
  use communication, only: comm_summation, comm_get_max, comm_is_root
  implicit none
  type(s_dg_fragment_rt), intent(inout) :: dg_frag
  type(s_dft_system),     intent(in)    :: system
  type(s_rgrid),          intent(in)    :: mg
  type(s_scalar),         intent(inout) :: rho
  type(s_scalar),         intent(inout) :: rho_s(system%nspin)
  character(len=*),       intent(in)    :: label
  real(8), allocatable :: rho_ref(:, :, :)
  real(8), allocatable :: rho_s_ref(:, :, :, :)
  real(8) :: drho2_local, drho2_sum, rho2_local, rho2_sum
  real(8) :: rho_max_local(2), rho_max_sum(2)
  real(8) :: diff, rel_l2
  integer :: ix, iy, iz, ispin
  integer :: ix_lo, ix_hi, iy_lo, iy_hi, iz_lo, iz_hi

  ix_lo = lbound(rho%f, 1)
  ix_hi = ubound(rho%f, 1)
  iy_lo = lbound(rho%f, 2)
  iy_hi = ubound(rho%f, 2)
  iz_lo = lbound(rho%f, 3)
  iz_hi = ubound(rho%f, 3)
  allocate(rho_ref(ix_lo:ix_hi, iy_lo:iy_hi, iz_lo:iz_hi))
  allocate(rho_s_ref(ix_lo:ix_hi, iy_lo:iy_hi, iz_lo:iz_hi, system%nspin))
  rho_ref(:, :, :) = rho%f(:, :, :)
  do ispin = 1, system%nspin
    rho_s_ref(:, :, :, ispin) = rho_s(ispin)%f(:, :, :)
  end do

  call diagnose_density_from_fragments(dg_frag, system, mg, rho, rho_s)

  drho2_local = 0.0d0
  rho2_local = 0.0d0
  rho_max_local(:) = 0.0d0
  do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
      do ix = mg%is(1), mg%ie(1)
        diff = rho%f(ix, iy, iz) - rho_ref(ix, iy, iz)
        drho2_local = drho2_local + diff * diff
        rho2_local = rho2_local + rho_ref(ix, iy, iz) * rho_ref(ix, iy, iz)
        rho_max_local(1) = max(rho_max_local(1), abs(diff))
        rho_max_local(2) = max(rho_max_local(2), abs(rho_ref(ix, iy, iz)))
      end do
    end do
  end do
  drho2_local = drho2_local * system%Hvol
  rho2_local = rho2_local * system%Hvol
  call comm_summation(drho2_local, drho2_sum, dg_frag%icomm)
  call comm_summation(rho2_local, rho2_sum, dg_frag%icomm)
  call comm_get_max(rho_max_local, rho_max_sum, 2, dg_frag%icomm)
  rel_l2 = sqrt(drho2_sum / max(rho2_sum, 1.0d-300))
  if (comm_is_root(dg_frag%id)) then
    write(*,'(1x,a,a,4(a,1pe13.5))') '[DG-FLUX-RHO] label=', trim(label), &
      ' rel_l2=', rel_l2, ' abs_l2=', sqrt(max(drho2_sum, 0.0d0)), &
      ' drho_max=', rho_max_sum(1), ' rho_max=', rho_max_sum(2)
    flush(6)
  end if

  rho%f(:, :, :) = rho_ref(:, :, :)
  do ispin = 1, system%nspin
    rho_s(ispin)%f(:, :, :) = rho_s_ref(:, :, :, ispin)
  end do
  deallocate(rho_ref, rho_s_ref)
end subroutine diagnose_dg_flux_density_mismatch

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
  integer :: env_len, env_stat
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
  real(8) :: E_wrk(3), E_sum(3), etmp, sysvol, g(3), r(3), Gd
  complex(8) :: rho_e, rho_i
  logical, save :: trace_env_initialized = .false.
  logical, save :: enable_energy_helper_trace = .false.
  character(len=32) :: env_value

  if (.not. trace_env_initialized) then
    env_value = ''
    call get_environment_variable('SALMON_DG_ENERGY_HELPER_TRACE', env_value, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      if (env_value(1:1) == '1' .or. env_value(1:1) == 'y' .or. env_value(1:1) == 'Y' .or. &
          env_value(1:1) == 't' .or. env_value(1:1) == 'T') enable_energy_helper_trace = .true.
    end if
    trace_env_initialized = .true.
  end if

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
  call comm_summation(rho_g2_local, rho_g2_sum, dg_frag%icomm)
  call comm_summation(rho_g2_local, rho_g2_sum_r, info%icomm_r)
  call comm_summation(rho_gprobe_re_local, rho_gprobe_re_sum, n_gprobe, dg_frag%icomm)
  call comm_summation(rho_gprobe_im_local, rho_gprobe_im_sum, n_gprobe, dg_frag%icomm)
  call comm_summation(rho_gprobe_re_local, rho_gprobe_re_sum_r, n_gprobe, info%icomm_r)
  call comm_summation(rho_gprobe_im_local, rho_gprobe_im_sum_r, n_gprobe, info%icomm_r)
  call comm_get_max(rho_gmax_local, rho_gmax_sum, 1, dg_frag%icomm)
  call comm_get_max(rho_gmax_local, rho_gmax_sum_r, 1, info%icomm_r)

  energy%E_kin = dg_frag%energy_kinetic
  energy%E_ion_nloc = dg_frag%energy_nonlocal
  energy%E_h = E_sum(1)
  energy%E_ion_loc = E_sum(2) + E_sum(3)
  energy%E_tot = dg_frag%total_energy - rho_vh_sum - rho_vxc_sum + energy%E_h + energy%E_xc + energy%E_ion_ion
  if (enable_energy_helper_trace .and. comm_is_root(nproc_id_global) .and. &
      (itt == 1 .or. mod(itt, 10) == 0)) then
    write(*,'(1x,a,i0,22(a,1pe14.6))') &
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
  use parallelization, only: nproc_id_global
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
  use parallelization, only: nproc_id_global
  use salmon_global, only: base_directory, sysname, yn_spinorbit
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
    write(iunit,'(a,a,a,a,a)') '# 1:x[', trim(t_unit_length%name), '] 2:y[', &
      trim(t_unit_length%name), '] 3:local_chern_marker_zint[none]'
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
  use communication, only: comm_is_root
  use salmon_global, only: iperiodic, yn_jm, yn_dg_fragment_rt
  implicit none
  !(header of standard output)
  if(comm_is_root(nproc_id_global))then
    if (yn_dg_fragment_rt == 'y') return
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
