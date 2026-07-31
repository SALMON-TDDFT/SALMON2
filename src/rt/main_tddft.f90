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
#ifdef USE_MPI
use mpi, only: MPI_Comm_rank,MPI_Bcast,MPI_INTEGER
#endif
use salmon_global
use structures
use parallelization, only: adjust_elapse_time, nproc_group_global
use communication, only: comm_is_root, comm_sync_all, comm_bcast, comm_summation
use salmon_xc, only: finalize_xc
use timer
use write_sub, only: write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d, &
  write_dg_polarization_data, write_dg_polarization_response_3d
use initialization_rt_sub
use checkpoint_restart_sub
use jellium, only: check_condition_jm
use rt_angular_momentum, only: write_local_angular_momentum_xy, flush_local_angular_momentum_xy
use rt_local_chern_marker, only: compute_local_chern_marker_from_orbital
use dg_overlapping_wannier_checkpoint, only: s_dg_overlapping_wannier_checkpoint, &
  read_dg_overlapping_wannier_checkpoint
use rt_dg_overlapping_wannier, only: s_dg_overlapping_wannier_rt_state, &
  initialize_dg_overlapping_wannier_rt,advance_dg_overlapping_wannier_rt,&
  write_dg_overlapping_wannier_rt_restart,read_dg_overlapping_wannier_rt_restart,&
  evaluate_dg_overlapping_wannier_observables,&
  write_dg_overlapping_wannier_rt_observable_sample
use em_field, only: calc_Ac_ext_t
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

if(yn_dg_overlapping_wannier_rt=='y')then
  call run_dg_overlapping_wannier_coefficient_rt()
  return
endif

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
    if (yn_dg_length_gauge == 'y') call write_dg_polarization_response_3d(ofl)
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

subroutine run_dg_overlapping_wannier_coefficient_rt()
  type(s_dg_overlapping_wannier_checkpoint)::checkpoint
  type(s_dg_overlapping_wannier_rt_state)::state
  complex(8),allocatable::coefficients(:,:)
  real(8),allocatable::vector_potential_samples(:,:)
  real(8)::electric_field(3),vector_potential(3),acceptance_gates(6),&
    polarization(3),current(3),cell_volume
  integer,allocatable::row_ids(:)
  integer::step,rank,ierr
  logical::ok,reusable
  character(256)::message
#ifdef USE_MPI
  call MPI_Comm_rank(nproc_group_global,rank,ierr)
#else
  rank=0
#endif
  acceptance_gates=[dg_dc_gs_final_density_tolerance,dg_dc_gs_final_orbital_tolerance,&
    10d0*dg_dc_gs_final_orbital_tolerance,dg_dc_gs_electron_count_tolerance,&
    1d0/dg_dc_metric_rank_tolerance,dg_ow_symmetry_tolerance]
  call read_dg_overlapping_wannier_checkpoint(nproc_group_global,'./overlapping_wannier_gs',&
    0,0,0_8,0_8,acceptance_gates,checkpoint,reusable,ok,message)
  if(.not.ok.or..not.reusable)then
    if(rank==0)write(0,'(a)')trim(message)
    error stop 'accepted V3 overlapping-Wannier checkpoint is required'
  endif
  if(trim(checkpoint%field_coupling_convention)/='cell_wrapped_length_velocity')&
    error stop 'unsupported overlapping-Wannier field convention'
  cell_volume=product(al)
  if(cell_volume<=0d0)error stop 'invalid overlapping-Wannier RT cell volume'
  ! The V3 reader certifies that every retained tail covers each physical
  ! periodic-grid id at least once; overlapping buffers may repeat IDs.
  ! With basis updates forbidden, every
  ! coefficient combination remains in that closed periodic support, so
  ! a nonzero representational tail escape is structurally impossible.
  allocate(row_ids(size(checkpoint%overlap_row_ids)));row_ids=int(checkpoint%overlap_row_ids)
  allocate(coefficients,source=checkpoint%coefficients)
  call initialize_dg_overlapping_wannier_rt(nproc_group_global,row_ids,checkpoint%overlap,&
    checkpoint%hamiltonian0,checkpoint%position,checkpoint%velocity,&
    checkpoint%basis_generation,checkpoint%geometry_generation,checkpoint%basis_fingerprint,&
    checkpoint%operator_fingerprint,checkpoint%hamiltonian_fingerprint,&
    checkpoint%observable_fingerprint,checkpoint%field_coupling_convention,&
    checkpoint%basis_generation,checkpoint%geometry_generation,&
    checkpoint%basis_fingerprint,checkpoint%operator_fingerprint,coefficients,state,ok,message)
  if(.not.ok)then
    if(rank==0)write(0,'(a)')trim(message)
    error stop 'overlapping-Wannier coefficient RT initialization failed'
  endif
  if(yn_dg_overlapping_wannier_rt_restart=='y')then
    call read_dg_overlapping_wannier_rt_restart(nproc_group_global,&
      './overlapping_wannier_rt.restart',coefficients,state,ok,message)
    if(.not.ok)then
      if(rank==0)write(0,'(a)')trim(message)
      error stop 'overlapping-Wannier coefficient RT restart failed'
    endif
  endif
  electric_field=0d0
  call evaluate_dg_overlapping_wannier_observables(nproc_group_global,coefficients,&
    checkpoint%occupations,cell_volume,state,polarization,current,ok,message)
  if(.not.ok)then
    if(rank==0)write(0,'(a)')trim(message)
    error stop 'overlapping-Wannier coefficient RT observable evaluation failed'
  endif
  call write_dg_overlapping_wannier_rt_observable_sample(nproc_group_global,&
    './overlapping_wannier_rt_observables.dat',electric_field,polarization,current,&
    cell_volume,state,yn_dg_overlapping_wannier_rt_restart=='y',ok,message)
  if(.not.ok)then
    if(rank==0)write(0,'(a)')trim(message)
    error stop 'overlapping-Wannier coefficient RT observable publication failed'
  endif
  allocate(vector_potential_samples(3,0:nt+1))
  call calc_Ac_ext_t(0d0,dt,0,nt+1,vector_potential_samples)
  do step=state%step+1,nt
    if(step==1.and.state%step==0.and.trim(ae_shape1)=='impulse')then
      ! The SALMON impulse is a vector-potential jump at t=0.  The value just
      ! before the first coefficient interval is zero, not the already-jumped
      ! sample stored at index zero.
      electric_field=-vector_potential_samples(:,step)/dt
    else
      electric_field=-(vector_potential_samples(:,step)-vector_potential_samples(:,step-1))/dt
    endif
    vector_potential=0d0
    call advance_dg_overlapping_wannier_rt(nproc_group_global,dt,electric_field,&
      vector_potential,coefficients,state,ok,message)
    if(.not.ok)then
      if(rank==0)write(0,'(a,i0,2a)')'coefficient RT failed at step ',step,': ',trim(message)
      error stop 'overlapping-Wannier coefficient RT propagation failed'
    endif
    call evaluate_dg_overlapping_wannier_observables(nproc_group_global,coefficients,&
      checkpoint%occupations,cell_volume,state,polarization,current,ok,message)
    if(.not.ok)then
      if(rank==0)write(0,'(a,i0,2a)')'observable evaluation failed at step ',step,': ',trim(message)
      error stop 'overlapping-Wannier coefficient RT observable evaluation failed'
    endif
    call write_dg_overlapping_wannier_rt_observable_sample(nproc_group_global,&
      './overlapping_wannier_rt_observables.dat',electric_field,polarization,current,&
      cell_volume,state,yn_dg_overlapping_wannier_rt_restart=='y',ok,message)
    if(.not.ok)then
      if(rank==0)write(0,'(a,i0,2a)')'observable publication failed at step ',step,': ',trim(message)
      error stop 'overlapping-Wannier coefficient RT observable publication failed'
    endif
  enddo
  call write_dg_overlapping_wannier_rt_restart(nproc_group_global,&
    './overlapping_wannier_rt.restart',coefficients,state,ok,message)
  if(.not.ok)then
    if(rank==0)write(0,'(a)')trim(message)
    error stop 'cannot publish overlapping-Wannier coefficient RT restart'
  endif
end subroutine

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
  use salmon_global, only: iperiodic, yn_jm
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
