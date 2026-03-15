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
  if (info%isize_r > 1) stop 'yn_out_lcm_rt=y currently requires orbital-only decomposition (isize_r=1)'
  call write_local_chern_marker_xy(0, mg, system, info, spsi_in)
end if
if (yn_out_lz_rt == 'y') then
  if (.not. singlescale%flag_use) stop 'yn_out_lz_rt=y requires theory=single_scale_maxwell_tddft'
  call write_local_angular_momentum_xy(0, lg, mg, system, info, singlescale)
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
          call write_local_angular_momentum_xy(itt, lg, mg, system, info, singlescale)
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
                            calculate_hamiltonian_matrix_std => calculate_hamiltonian_matrix
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
  
  ! Initialize DG-Fragment RT
  if (yn_spinorbit == 'y') then
    call init_dg_fragment_rt_soi(dg_frag, system, rt, info, lg, mg)
  else
    call init_dg_fragment_rt_std(dg_frag, system, rt, info, lg, mg)
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
    
    ! Store observables from DG-Fragment calculation (already calculated in tddft_dg_fragment_iteration)
    energy%E_tot = dg_frag%total_energy
    energy%elec_num = dg_frag%elec_num_scaled
    energy%elec_num_raw = dg_frag%elec_num_raw
    energy%pw_weight_raw = dg_frag%pw_weight_raw

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
    if (mod(itt, 1) == 0) then  ! Output at every timestep
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
        ! curr_e needs shape (3, 2) for two spins, use same current for both spins
        call write_rt_data_3d(itt, ofl, dt, system, &
                              reshape([dg_frag%current(1), dg_frag%current(2), dg_frag%current(3), &
                                       dg_frag%current(1), dg_frag%current(2), dg_frag%current(3)], [3,2]), &
                              dg_frag%current(1:3))
      end select
      
      ! Write energy data
      call write_rt_energy_data(itt, ofl, dt, energy, md)
    endif
    
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

subroutine write_local_chern_marker_xy(itt, mg, system, info, psi_fin)
  use structures, only: s_rgrid, s_dft_system, s_parallel_info, s_orbital
  use communication, only: comm_is_root
  use rt_local_chern_marker, only: compute_local_chern_marker_from_orbital
  implicit none
  integer, intent(in) :: itt
  type(s_rgrid), intent(in) :: mg
  type(s_dft_system), intent(in) :: system
  type(s_parallel_info), intent(in) :: info
  type(s_orbital), intent(in) :: psi_fin
  real(8), allocatable :: marker(:,:,:)
  real(8), allocatable :: marker_xy(:,:)
  character(256) :: filename, filenum
  integer :: ix, iy, iz, iunit

  allocate(marker(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
  allocate(marker_xy(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2)))
  call compute_local_chern_marker_from_orbital(mg, system, info, psi_fin, marker)

  marker_xy(:,:) = 0.0d0
  do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
      do ix = mg%is(1), mg%ie(1)
        marker_xy(ix,iy) = marker_xy(ix,iy) + marker(ix,iy,iz) * system%hgs(3)
      end do
    end do
  end do

  if (comm_is_root(nproc_id_global)) then
    write(filenum, '(i6.6)') itt
    filename = trim(base_directory)//trim(sysname)//'_lcm_xy_'//trim(adjustl(filenum))//'.data'
    open(newunit=iunit, file=trim(filename), status='replace', action='write')
    write(iunit,'(a)') '# Local Chern marker integrated along z'
    write(iunit,'(a)') '# x: x coordinate'
    write(iunit,'(a)') '# y: y coordinate'
    write(iunit,'(a)') '# local_chern_marker_zint: local Chern marker integrated over z'
    write(iunit,'(a)') '# 1:x[a.u.] 2:y[a.u.] 3:local_chern_marker_zint[none]'
    do iy = mg%is(2), mg%ie(2)
      do ix = mg%is(1), mg%ie(1)
        write(iunit,'(3(1x,es24.16))') mg%coordinate(ix,1), mg%coordinate(iy,2), marker_xy(ix,iy)
      end do
      write(iunit,*)
    end do
    close(iunit)
  end if
  deallocate(marker_xy, marker)
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
