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
#ifdef USE_MPI
use mpi, only: MPI_Allreduce,MPI_Allgather,MPI_Allgatherv,MPI_INTEGER8,MPI_BXOR,MPI_SUCCESS,MPI_INTEGER,&
  MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_MIN,MPI_MAX
#endif
use structures
use inputoutput
use salmon_global, only: yn_dc_lcfo_flux, yn_dc_lcfo_wannier, yn_dg_wpw_production, &
  yn_dg_dc_local_periodic, yn_dg_dc_overlapping_wannier, ncg, base_directory, &
  dg_dc_handoff_min_iter, dg_dc_handoff_tolerance, &
  dg_dc_candidate_orbitals_per_atom, dg_dc_metric_rank_tolerance, &
  dg_dc_gs_intermediate_orbital_tolerance,dg_dc_gs_intermediate_density_tolerance, &
  dg_dc_gs_final_orbital_tolerance,dg_dc_gs_final_density_tolerance,dg_dc_gs_subspace_tolerance, &
  dg_dc_gs_initial_lambda_step,dg_dc_gs_minimum_lambda_step,dg_dc_gs_maximum_lambda_step, &
  dg_dc_gs_allowed_residual_growth,dg_dc_gs_density_mix_rate,dg_dc_gs_hermiticity_tolerance, &
  dg_dc_gs_orthogonality_tolerance,dg_dc_gs_face_balance_tolerance,dg_dc_gs_electron_count_tolerance, &
  dg_dc_gs_minimum_projector_overlap,dg_dc_gs_maximum_scf_iterations, &
  dg_dc_gs_maximum_eigensolver_iterations,dg_dc_gs_maximum_rollbacks, &
  dg_dc_gs_sipg_penalty_factor,dg_dc_gs_target_lambda,dg_ow_boundary_value_tolerance,&
  dg_ow_boundary_gradient_tolerance,dg_ow_symmetry_tolerance
use dg_dc_handoff, only: dg_dc_handoff_runtime, dg_dc_nodal_runtime, initialize_dg_dc_handoff, &
  materialize_dg_dc_candidates
use dg_dc_ground_state, only: s_dg_dc_gs_controls,s_dg_dc_gs_result,s_dg_dc_gs_diagnostics, &
  default_dg_dc_gs_controls,run_dg_dc_ground_state,DG_DC_ACCEPTED
use dg_dc_ground_state_adapter, only: expand_dg_dc_global_candidate_axis,reconstruct_dg_dc_core_density, &
  mix_dg_dc_density_transaction,initialize_dg_dc_physical_faces,apply_dg_dc_sipg_operator_mpi, &
  compose_dg_dc_hamiltonian,build_dg_dc_interior_volume_action,execute_dg_dc_production_iteration
use dg_nodal_state, only: s_dg_nodal_common_state
use dg_dc_local_basis_ground_state, only: s_dg_dc_local_basis_layout, &
  initialize_dg_dc_local_basis_layout,orthonormalize_dg_dc_fragment_core_basis, &
  transform_dg_dc_fragment_buffer_basis, &
  project_dg_dc_local_basis_volume,assemble_dg_dc_local_basis_interface_rows, &
  compose_dg_dc_distributed_hamiltonian_rows,initialize_dg_dc_local_basis_coefficients, &
  assign_dg_dc_local_basis_occupations,solve_dg_dc_local_basis_bands_cg, &
  reconstruct_dg_dc_local_basis_density,validate_dg_dc_local_basis_density, &
  diagnose_dg_dc_local_basis_continuation, &
  diagnose_dg_dc_six_face_balance,build_dg_dc_six_face_neighbors, &
  s_dg_dc_local_basis_production_state,dg_dc_local_basis_state
use dg_ground_state_checkpoint, only: s_dg_ground_state_checkpoint,s_dg_dc_direct_checkpoint_state, &
  populate_dg_ground_state_checkpoint,publish_dg_ground_state_checkpoint,publish_dg_dc_direct_checkpoint
use dg_overlapping_wannier_construction, only: s_dg_overlapping_wannier_construction, &
  construct_dg_overlapping_wannier_basis,verify_dg_overlapping_wannier_periodic_closure
use dg_overlapping_wannier_metric, only: assemble_dg_overlapping_wannier_metric
use dg_overlapping_wannier_operators, only: assemble_dg_overlapping_wannier_weak_operators
use dg_overlapping_wannier_nonlocal, only: assemble_dg_overlapping_wannier_nonlocal
use dg_overlapping_wannier_scf, only: s_dg_overlapping_wannier_scf_state, &
  s_dg_overlapping_wannier_scf_result,run_dg_overlapping_wannier_scf, &
  compute_dg_overlapping_wannier_scf_fingerprint
use dg_overlapping_wannier_checkpoint, only: s_dg_overlapping_wannier_checkpoint, &
  write_dg_overlapping_wannier_checkpoint,read_dg_overlapping_wannier_checkpoint
use rt_dg_nodal_cg, only: solve_nodal_ground_state_cg_mpi
use rt_dg_nodal_rayleigh_ritz, only: rayleigh_ritz_nodal_subspace_mpi
#ifdef USE_EIGENEXA
use eigenexa_module, only: finalize_eigenexa
#endif
use parallelization, only: nproc_id_global,nproc_group_global,adjust_elapse_time,nproc_size_global
use communication, only: comm_is_root, comm_summation, comm_bcast, comm_sync_all, comm_get_max, comm_logical_and
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
use jellium, only: check_condition_jm
use dcdft
use dcdft_soi
use lcfo
use lcfo_flux
use lcfo_soi
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
type(s_scalar) :: rho,rho_jm,Vh,Vpsl
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
type(s_dcdft) :: dc

logical :: rion_update
logical :: flag_opt_conv
logical :: local_basis_route_active
integer :: Miopt, iopt,nopt_max,i
integer :: iter_band_kpt, iter_band_kpt_end, iter_band_kpt_stride
logical :: is_checkpoint_iter, is_shutdown_time
logical :: dg_handoff_ok
character(256) :: dg_handoff_message
type(s_dg_dc_gs_controls) :: dg_gs_controls
type(s_dg_dc_gs_result) :: dg_gs_result
real(8), allocatable :: dg_gs_density(:,:,:,:),dg_gs_raw_density(:,:,:,:),dg_gs_core_density(:,:,:,:)
real(8), allocatable :: dg_gs_mixed_density(:,:,:,:)
real(8), allocatable :: dg_gs_eigenvalues(:,:)
real(8) :: dg_gs_current_lambda,dg_gs_face_hermiticity,dg_gs_face_balance
integer, allocatable :: dg_gs_fragment_origins(:,:),dg_gs_fragment_sizes(:,:)
integer(8) :: dg_gs_potential_epoch,dg_gs_hamiltonian_potential_epoch
logical :: dg_gs_operator_ok
character(256) :: dg_gs_operator_message
type(s_dg_overlapping_wannier_construction) :: ow_basis
type(s_dg_overlapping_wannier_scf_state) :: ow_state
type(s_dg_overlapping_wannier_scf_result) :: ow_result
type(s_dg_overlapping_wannier_checkpoint) :: ow_checkpoint
complex(8),allocatable :: ow_s(:,:),ow_core_values(:,:),ow_core_gradients(:,:,:),&
  ow_box_values(:,:),ow_box_gradients(:,:,:)
integer(8),allocatable :: ow_core_ids(:),ow_row_ids(:)
integer,allocatable :: ow_tail_generation(:,:)
integer,allocatable :: ow_core_box_positions(:)
real(8),allocatable :: ow_core_weights(:)
integer :: ow_box_size(3),ow_core_size(3),ow_buffer(3)
integer(8) :: ow_symmetry_fingerprint
integer(8) :: ow_potential_epoch_snapshot
integer(8) :: ow_global_grid_count
real(8),allocatable :: ow_density_snapshot(:,:,:,:)
real(8),allocatable :: ow_work_density(:,:,:,:)
logical :: ow_transaction_active
integer :: ilevel_print

if(theory=='dft_band'.and.iperiodic/=3) return

if(yn_dc=='y') then
  call initialize_dg_dc_handoff(dg_dc_handoff_runtime,yn_dg_dc_local_periodic=='y', &
    dg_dc_handoff_min_iter,dg_dc_handoff_tolerance,dg_dc_candidate_orbitals_per_atom, &
    dg_dc_metric_rank_tolerance,dg_handoff_ok,dg_handoff_message)
  if(.not.dg_handoff_ok) stop 'invalid DG DC handoff controls'
  if(yn_dg_wpw_production=='y'.and.yn_dc_lcfo_flux/='y') &
    stop 'DG WPW production requires yn_dc_lcfo_flux=y'
  if(yn_spinorbit=='y') then
    call init_dcdft_soi(dc,pp,mixing,ewald)
  else
    call init_dcdft(dc,pp,mixing,ewald)
  end if
  ilevel_print = 0
else
  ilevel_print = 3
end if

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
                          rho, rho_jm, rho_s, Vh, V_local, Vpsl, Vxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,  &
                          ofl )

call initialization2_dft( Miter, nspin, rion_update,  &
                          system, energy, ewald, stencil, fg, poisson,&
                          lg, mg, info,   &
                          srg, srg_scalar,  &
                          rho, rho_jm, rho_s, Vh,V_local, Vpsl, Vxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,   &
                          xc_func, mixing )

if(yn_dg_dc_local_periodic=='y')then
  call comm_logical_and(system%if_real_orbital.and.allocated(spsi%rwf),dg_handoff_ok,dc%icomm_tot)
  if(.not.dg_handoff_ok)stop 'DG DC local-periodic route requires allocated real orbitals collectively'
endif

Miopt = 0
nopt_max = 1
if(yn_opt=='y') call initialization_opt(Miopt,opt,system,flag_opt_conv,nopt_max,ofl)

call timer_end(LOG_INIT_GS)

if(yn_dc == 'y' .and. yn_dc_lcfo_wannier == 'y' .and. dc_lcfo_wannier_import_only_requested()) then
  if(comm_is_root(nproc_id_global)) &
    write(*,'(1x,a)') '[DC-LCFO-W90-IMPORT] import-only mode: skip SCF and reuse external Wannier90 outputs'
  call dc_lcfo_wannier_import_only(dc)
  call timer_end(LOG_TOTAL)
  return
end if

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
   call init_band_dft(system,band) ! --> system%wtk=0.0
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

call timer_end(LOG_INIT_GS_ITERATION)


call timer_begin(LOG_GS_ITERATION)
!------------------------------------ SCF Iteration
!Iteration loop for SCF (DFT_Iteration)
call scf_iteration_dft( Miter,rion_update,sum1,  &
                        system,energy,ewald,  &
                        lg,mg,  &
                        info,  &
                        poisson,fg,  &
                        cg,mixing,  &
                        stencil,  &
                        srg,srg_scalar,   &
                        spsi,shpsi,sttpsi,  &
                        rho,rho_jm,rho_s,  &
                        V_local,Vh,Vxc,Vpsl,xc_func,  &
                        pp,ppg,ppn,  &
                        band, ilevel_print, dc)


if(theory=='dft_band')then
   call write_band(system,energy)
end if

! output the wavefunctions for next GS calculations
if(write_gs_wfn_k == 'y') then !this input keyword is going to be removed....
   select case(iperiodic)
   case(3)
      call write_wfn(lg,mg,spsi,info,system)
      ! Experimental Implementation of Inner-Product Outputs:
      ! call write_prod_dk_data(lg, mg, system, info, spsi)
   case(0)
      write(*,*) "error: write_gs_wfn_k='y' & iperiodic=0"
   end select
end if

! output transition moment : --> want to put out of the optmization loop in future
if(yn_out_tm  == 'y'.or.yn_out_gs_sgm_eps=='y') then
   select case(iperiodic)
   case(3)
      call write_k_data(system,stencil)  !need? (probably remove later)
      call write_tm_data(spsi,system,info,mg,stencil,srg,ppg,energy)
   case(0)
     write(*,*) "error: yn_out_tm='y',yn_out_gs_sgm_eps='y' & iperiodic=0"
  end select
end if

   ! force
   if(yn_jm=='n' .and. yn_dc=="n")then
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
   call write_xyz(comment_line,"add","r  ",system,ofl)
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

      if(yn_dg_dc_local_periodic/='y' .and. yn_dg_dc_overlapping_wannier/='y' .and. &
         (is_checkpoint_iter .or. is_shutdown_time)) then
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
local_basis_route_active=.false.

if(yn_dc=='y') then
  if(yn_dg_dc_overlapping_wannier == 'y') then
    local_basis_route_active=.true.
    call run_dg_overlapping_wannier_ground_state_for_main
    return
  else if(yn_dg_dc_local_periodic == 'y') then
    local_basis_route_active=.true.
    if(.not.dg_dc_handoff_runtime%direct_ground_state_complete) &
      stop 'DG DC route reached publication without a full-scale fixed point'
    call publish_dg_dc_direct_ground_state_for_main(dg_handoff_ok,dg_handoff_message)
    if(.not.dg_handoff_ok)stop 'DG DC direct checkpoint publication failed collectively'
  else if(yn_dc_lcfo_flux == 'y') then
    if(yn_spinorbit == 'y') then
      stop "yn_dc_lcfo_flux=y is not implemented for spin-orbit mode"
    else
      if(comm_is_root(nproc_id_global)) &
      & write(*,'(1x,a)') '[DC-LCFO-FLUX] export phase: build Flux-LCFO basis and coefficients'
#ifdef USE_EIGENEXA
      call finalize_eigenexa(info)
#endif
      call dc_lcfo_flux(lg,mg,system,info,stencil,xc_func,pp,ppn,ppg,energy,rho_s,v_local,&
        spsi,shpsi,sttpsi,srg,dc)
    end if
  else if(yn_dc_lcfo == 'y') then
    if(yn_spinorbit == 'y') then
      call dc_lcfo_soi(lg,mg,system,info,stencil,ppg,energy,v_local,spsi,shpsi,sttpsi,srg,dc)
    else
      call dc_lcfo(lg,mg,system,info,stencil,ppg,energy,v_local,spsi,shpsi,sttpsi,srg,dc)
    end if
  end if
  if(.not.local_basis_route_active .and. yn_spinorbit == 'y') then
    call write_total_dcdft_soi(system,dc)
  else if(.not.local_basis_route_active) then
    call write_total_dcdft(system,dc)
  end if
end if

! write GS: basic data
if(.not.local_basis_route_active) then
if(yn_dc=='n') call write_band_information(system,energy)
call write_eigen(ofl,system,energy)
call write_info_data(Miter,system,energy,pp)
call write_k_data(system,stencil)
if(yn_spinorbit=='y') call write_mag_decomposed_gs(system,mg,info,spsi)

! write GS: analysis option
if(yn_out_psi =='y') call write_psi(lg,mg,system,info,spsi)
if(yn_out_dns =='y') call write_dns(lg,mg,system,info,rho_s)
if(yn_out_dos =='y') call write_dos(system,energy)
if(yn_out_pdos=='y') call write_pdos(lg,mg,system,info,pp,energy,spsi)
if(yn_out_elf =='y') call write_elf(0,lg,mg,system,info,stencil,rho,srg,srg_scalar,spsi)
end if

call timer_end(LOG_WRITE_GS_RESULTS)

! write GS: binary data for restart
call timer_begin(LOG_WRITE_GS_DATA)
if(.not.local_basis_route_active) then
if(write_gs_restart_data=="no") then
   if(comm_is_root(nproc_id_global)) &
      write(*,'(a)')"  no restart data writing."
else if(write_gs_restart_data.ne."checkpoint_only") then
   if(comm_is_root(nproc_id_global)) write(*,'(a)')"  writing restart data..."
   call checkpoint_gs(lg,mg,system,info,spsi,Miter,mixing,ofl%dir_out_restart)
   call comm_sync_all
   if(yn_opt=='y') then
      if(.not.flag_opt_conv) then
         call comm_sync_all
         call checkpoint_opt(nopt_max,opt,ofl%dir_out_restart)
         call comm_sync_all
      endif
   endif
else
   if(yn_self_checkpoint=='n') then
      if(comm_is_root(nproc_id_global)) then
           write(*,'(a)')"  no restart data writing:"
           write(*,'(a)')"  check input keywords if you need restart data"
       endif
   endif
endif
if(yn_self_checkpoint=='y') then
   if(comm_is_root(nproc_id_global)) &
   write(*,'(a)')"  writing restart data in checkpoint format ..."
   call checkpoint_gs(lg,mg,system,info,spsi,Miter,mixing)
   call comm_sync_all
endif
if(comm_is_root(nproc_id_global)) write(*,'(a)')"  writing completed."
else if(comm_is_root(nproc_id_global)) then
  write(*,'(a)')'  DG local-basis result retained in memory; standard restart publication skipped.'
end if
call timer_end(LOG_WRITE_GS_DATA)

!call timer_begin(LOG_WRITE_GS_INFO)  !if needed, please take back, sory: AY
!call timer_end(LOG_WRITE_GS_INFO)

call finalize_xc(xc_func)

if(yn_dc=='y') then
! override (restore)
  nproc_group_global = dc%icomm_tot
  nproc_id_global = dc%id_tot
  nproc_size_global = dc%isize_tot
  call comm_sync_all
  call finalize_dcdft(dc)
end if

call timer_end(LOG_TOTAL)

contains

  subroutine checked_ow_extent_product(extent,value,ok)
    integer,intent(in)::extent(3)
    integer(8),intent(out)::value
    logical,intent(out)::ok
    integer::axis
    value=1_8;ok=all(extent>0)
    if(.not.ok)return
    do axis=1,3
      if(value>huge(value)/int(extent(axis),8))then
        value=0_8;ok=.false.;return
      endif
      value=value*int(extent(axis),8)
    enddo
  end subroutine

  subroutine run_dg_overlapping_wannier_ground_state_for_main()
    complex(8),allocatable::candidate(:,:),candidate_gradient(:,:,:),occupied_coefficients(:,:),&
      periodic_phase(:,:),retained_vectors(:,:)
    real(8),allocatable::weights(:),coordinate(:),spectrum(:),occupations(:),gradient_rotation(:,:,:)
    integer(8),allocatable::physical_ids(:),box_ids(:),symmetry_map(:,:)
    integer,allocatable::fragments(:)
    logical,allocatable::boundary(:),core_mask(:),pairs(:,:)
    integer::ix,iy,iz,io,p,nbox,ncore,ncandidate,noccupied,ntarget,nsym,rank,nproc,&
      raw_ix,raw_iy,raw_iz,core_index,rejected_rank,ownership_count,ierr,allocation_status,&
      local_candidate_count,local_occupied_count,candidate_offset
    integer(8)::expected_core_count,expected_box_count,basis_fingerprint,nbox8,ncore8,product8,nxy8
    real(8)::minimum_eigenvalue,condition_number,closure_residual
    logical::ok,reusable
    character(256)::message,prefix

    call MPI_Comm_rank(dc%icomm_tot,rank,ierr);call MPI_Comm_size(dc%icomm_tot,nproc,ierr)
    ok=system%nspin==1.and.system%if_real_orbital.and.allocated(spsi%rwf)
    call comm_logical_and(ok,reusable,dc%icomm_tot)
    if(.not.reusable)error stop 'overlapping-Wannier production requires Gamma real DC candidates'
    ok=nproc==dc%n_frag.and..not.dc%optimized_fragment_geometry
    call comm_logical_and(ok,reusable,dc%icomm_tot)
    if(.not.reusable)&
      error stop 'overlapping-Wannier production requires one rank per valid DC fragment'
    ow_core_size=dc%nxyz_domain_frag(:,dc%i_frag);ow_buffer=dc%nxyz_buffer
    do ix=1,3
      if(ow_buffer(ix)>huge(ow_box_size(ix))/2)&
        error stop 'overlapping-Wannier buffer extent overflow'
      if(ow_core_size(ix)>huge(ow_box_size(ix))-2*ow_buffer(ix))&
        error stop 'overlapping-Wannier buffer extent overflow'
    enddo
    ow_box_size=ow_core_size+2*ow_buffer
    ok=all(ow_buffer>=size(stencil%coef_nab,1).or.ow_core_size==dc%lg_tot%num)
    call comm_logical_and(ok,reusable,dc%icomm_tot)
    if(.not.reusable)error stop 'overlapping-Wannier buffer is smaller than the SALMON gradient stencil'
    ok=.not.any([size(spsi%rwf,1),size(spsi%rwf,2),size(spsi%rwf,3)]<ow_box_size)
    call comm_logical_and(ok,reusable,dc%icomm_tot)
    if(.not.reusable)&
      error stop 'overlapping-Wannier candidate buffer box is incomplete'
    call checked_ow_extent_product(ow_box_size,nbox8,ok)
    call checked_ow_extent_product(ow_core_size,ncore8,reusable);ok=ok.and.reusable
    if(.not.ok.or.nbox8>int(huge(nbox),8).or.ncore8>int(huge(ncore),8))&
      error stop 'overlapping-Wannier grid extent overflow'
    nbox=int(nbox8);ncore=int(ncore8)
    local_candidate_count=system%no
    local_occupied_count=count(system%rocc(1:local_candidate_count,info%ik_s,1)>1d-12)
    if(local_candidate_count>huge(ncandidate)/nproc.or.local_occupied_count>huge(noccupied)/nproc)&
      error stop 'overlapping-Wannier candidate direct-sum extent overflow'
    ncandidate=local_candidate_count*nproc;noccupied=local_occupied_count*nproc
    if(nbox8>huge(product8)/int(ncandidate,8))&
      error stop 'overlapping-Wannier candidate payload overflow'
    product8=int(ncandidate,8)*nbox8
    if(product8>int(huge(nbox),8)/3_8)&
      error stop 'overlapping-Wannier candidate payload exceeds collective count'
    candidate_offset=(dc%i_frag-1)*local_candidate_count
    ntarget=ncandidate;nsym=nproc
    if(noccupied<1.or.noccupied>ntarget)error stop 'invalid overlapping-Wannier occupied window'
    call checked_ow_extent_product(dc%lg_tot%num,expected_core_count,ok)
    if(.not.ok.or.expected_core_count>int(huge(nbox),8))&
      error stop 'overlapping-Wannier global grid exceeds addressable extent'
    ow_global_grid_count=expected_core_count
    nxy8=int(dc%lg_tot%num(1),8)*int(dc%lg_tot%num(2),8)
    if(.not.ok.or.nbox8>huge(expected_box_count)/int(nproc,8))&
      error stop 'overlapping-Wannier global extent overflow'
    expected_box_count=nbox8*int(nproc,8)
    allocate(candidate(ncandidate,nbox),candidate_gradient(3,ncandidate,nbox),weights(nbox),&
      coordinate(nbox),periodic_phase(3,nbox),physical_ids(nbox),box_ids(nbox),&
      symmetry_map(nbox,nsym),fragments(nbox),boundary(nbox),core_mask(nbox),&
      occupied_coefficients(ncandidate,noccupied),gradient_rotation(3,3,nsym),stat=allocation_status)
    call comm_logical_and(allocation_status==0,reusable,dc%icomm_tot)
    if(.not.reusable)error stop 'overlapping-Wannier production allocation failed'
    candidate=(0d0,0d0);candidate_gradient=(0d0,0d0);weights=system%hvol
    occupied_coefficients=(0d0,0d0)
    do p=1,nproc
      do io=1,local_occupied_count
        occupied_coefficients((p-1)*local_candidate_count+io,(p-1)*local_occupied_count+io)=1d0
      enddo
    enddo
    gradient_rotation=0d0
    do io=1,nsym;do ix=1,3;gradient_rotation(ix,ix,io)=1d0;enddo;enddo
    p=0;core_index=0
    do iz=1,ow_box_size(3);do iy=1,ow_box_size(2);do ix=1,ow_box_size(1)
      p=p+1;box_ids(p)=int(dc%i_frag-1,8)*nbox8+int(p,8);fragments(p)=dc%i_frag
      core_mask(p)=ix>ow_buffer(1).and.ix<=ow_buffer(1)+ow_core_size(1).and.&
        iy>ow_buffer(2).and.iy<=ow_buffer(2)+ow_core_size(2).and.&
        iz>ow_buffer(3).and.iz<=ow_buffer(3)+ow_core_size(3)
      boundary(p)=ix==1.or.ix==ow_box_size(1).or.iy==1.or.iy==ow_box_size(2).or.&
        iz==1.or.iz==ow_box_size(3)
      raw_ix=canonical_to_dc_index(ix,ow_core_size(1),ow_buffer(1))
      raw_iy=canonical_to_dc_index(iy,ow_core_size(2),ow_buffer(2))
      raw_iz=canonical_to_dc_index(iz,ow_core_size(3),ow_buffer(3))
      do io=1,local_candidate_count
        candidate(candidate_offset+io,p)=cmplx(spsi%rwf(raw_ix,raw_iy,raw_iz,1,io,1,1),0d0,8)
      enddo
      physical_ids(p)=1_8+int(modulo(dc%ixyz_frag(1,dc%i_frag)-1+ix-ow_buffer(1)-1,dc%lg_tot%num(1)),8)+&
        int(dc%lg_tot%num(1),8)*(int(modulo(dc%ixyz_frag(2,dc%i_frag)-1+iy-ow_buffer(2)-1,&
        dc%lg_tot%num(2)),8)+int(dc%lg_tot%num(2),8)*int(modulo(dc%ixyz_frag(3,dc%i_frag)-1+&
        iz-ow_buffer(3)-1,dc%lg_tot%num(3)),8))
      coordinate(p)=real(physical_ids(p),8)
      periodic_phase(1,p)=exp(cmplx(0d0,2d0*pi*real(modulo(physical_ids(p)-1_8,&
        int(dc%lg_tot%num(1),8)),8)/real(dc%lg_tot%num(1),8),8))
      periodic_phase(2,p)=exp(cmplx(0d0,2d0*pi*real(modulo((physical_ids(p)-1_8)/&
        int(dc%lg_tot%num(1),8),int(dc%lg_tot%num(2),8)),8)/real(dc%lg_tot%num(2),8),8))
      periodic_phase(3,p)=exp(cmplx(0d0,2d0*pi*real((physical_ids(p)-1_8)/&
        nxy8,8)/real(dc%lg_tot%num(3),8),8))
      if(core_mask(p))then
        core_index=core_index+1
      endif
    enddo;enddo;enddo
    call build_dc_translation_symmetry_map(nbox,symmetry_map,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'invalid DC translation symmetry';endif
    call periodic_box_gradients(candidate,ow_box_size,stencil%coef_nab,candidate_gradient)
    call construct_dg_overlapping_wannier_basis(dc%icomm_tot,ncandidate,ntarget,noccupied,physical_ids,&
      fragments,weights,coordinate,boundary,candidate,candidate_gradient,occupied_coefficients,&
      expected_core_count,1,dg_ow_boundary_value_tolerance,dg_ow_boundary_gradient_tolerance,&
      dg_dc_metric_rank_tolerance,ow_basis,ok,message,core_mask,box_ids,symmetry_map,&
      expected_box_count,dg_ow_symmetry_tolerance,periodic_phase)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'overlapping-Wannier construction gate failed';endif
    ow_box_values=ow_basis%value;ow_box_gradients=ow_basis%gradient
    call verify_dg_overlapping_wannier_periodic_closure(dc%icomm_tot,box_ids,symmetry_map,ow_box_values,&
      ow_box_gradients,ow_basis%symmetry_representation,gradient_rotation,expected_box_count,&
      dg_ow_symmetry_tolerance,closure_residual,ow_symmetry_fingerprint,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'overlapping-Wannier symmetry gate failed';endif
    allocate(ow_core_values(ntarget,ncore),ow_core_gradients(3,ntarget,ncore),&
      ow_core_weights(ncore),ow_core_ids(ncore),&
      ow_core_box_positions(ncore));core_index=0
    do p=1,nbox
      if(.not.core_mask(p))cycle
      core_index=core_index+1;ow_core_values(:,core_index)=ow_box_values(:,p)
      ow_core_gradients(:,:,core_index)=ow_box_gradients(:,:,p)
      ow_core_weights(core_index)=weights(p);ow_core_ids(core_index)=physical_ids(p)
      ow_core_box_positions(core_index)=p
    enddo
    allocate(pairs(ntarget,ntarget));pairs=.true.
    call assemble_dg_overlapping_wannier_metric(dc%icomm_tot,ntarget,ow_core_ids,ow_core_weights,&
      ow_core_values,pairs,expected_core_count,dg_dc_metric_rank_tolerance,ow_s,retained_vectors,&
      spectrum,minimum_eigenvalue,condition_number,rejected_rank,ownership_count,ok,message)
    if(.not.ok.or.rejected_rank/=0)error stop 'overlapping-Wannier metric gate failed'
    allocate(ow_row_ids(count(ow_basis%center_owner_rank==rank)))
    io=0
    do p=1,ntarget
      if(ow_basis%center_owner_rank(p)/=rank)cycle
      io=io+1;ow_row_ids(io)=p
    enddo
    allocate(ow_tail_generation(ntarget,ncore));ow_tail_generation=ow_basis%generation
    call compute_dg_overlapping_wannier_scf_fingerprint(dc%icomm_tot,ow_row_ids,&
      ow_s(int(ow_row_ids),:),ow_core_ids,ow_core_weights,ow_core_values,ow_tail_generation,&
      basis_fingerprint,ow_symmetry_fingerprint)
    prefix=trim(base_directory)//'/overlapping_wannier_gs'
    call read_dg_overlapping_wannier_checkpoint(dc%icomm_tot,trim(prefix),ow_basis%generation,1,&
      basis_fingerprint,dg_dc_operator_fingerprint(.false.),[dg_dc_gs_final_density_tolerance,&
      dg_dc_gs_final_orbital_tolerance,10d0*dg_dc_gs_final_orbital_tolerance,&
      dg_dc_gs_electron_count_tolerance,1d0/dg_dc_metric_rank_tolerance,dg_ow_symmetry_tolerance],&
      ow_checkpoint,reusable,ok,message)
    if(ok.and.reusable)then
      call restore_ow_checkpoint_density(ow_checkpoint,ok,message)
      if(.not.ok)error stop 'overlapping-Wannier checkpoint density restore failed'
      if(comm_is_root(nproc_id_global))write(*,'(a)')'[OW-GS] reused accepted route checkpoint'
      return
    endif
    allocate(occupations(noccupied))
    do p=1,nproc
      occupations((p-1)*local_occupied_count+1:p*local_occupied_count)=&
        system%rocc(1:local_occupied_count,info%ik_s,1)
    enddo
    allocate(ow_state%density(ncore),ow_state%potential(ncore),ow_state%coefficients(ntarget,noccupied),&
      ow_state%eigenvalues(noccupied),ow_state%density_history(ncore,2))
    do p=1,ncore
      ow_state%density(p)=dc%rho_tot_s(1)%f(int(modulo(ow_core_ids(p)-1_8,int(dc%lg_tot%num(1),8)))+1,&
        int(modulo((ow_core_ids(p)-1_8)/int(dc%lg_tot%num(1),8),int(dc%lg_tot%num(2),8)))+1,&
        int((ow_core_ids(p)-1_8)/nxy8)+1)
      ow_state%potential(p)=0d0
    enddo
    ow_state%coefficients=(0d0,0d0);do io=1,noccupied;ow_state%coefficients(io,io)=1d0;enddo
    ow_state%eigenvalues=0d0;ow_state%density_history(:,1)=ow_state%density
    ow_state%density_history(:,2)=ow_state%density;ow_state%history_count=1
    ow_state%basis_generation=ow_basis%generation;ow_state%geometry_generation=1
    ow_state%basis_fingerprint=basis_fingerprint
    ow_state%operator_fingerprint=dg_dc_operator_fingerprint(.false.)
    call run_dg_overlapping_wannier_scf(dc%icomm_tot,ow_row_ids,ow_s(int(ow_row_ids),:),ow_core_ids,&
      ow_core_weights,ow_core_values,ow_tail_generation,ow_basis%generation,1,basis_fingerprint,&
      occupations,closure_residual,dg_ow_symmetry_tolerance,ow_symmetry_fingerprint,&
      expected_core_count,dg_dc_gs_density_mix_rate,dg_dc_gs_maximum_scf_iterations,&
      dg_dc_gs_maximum_eigensolver_iterations,dg_dc_gs_final_density_tolerance,&
      dg_dc_gs_final_orbital_tolerance,ow_build_hamiltonian,ow_mix_density,ow_transaction,&
      ow_commit_transaction,ow_state,ow_result,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'overlapping-Wannier SCF gate failed';endif
    call populate_ow_checkpoint(occupations,condition_number,closure_residual)
    call write_dg_overlapping_wannier_checkpoint(dc%icomm_tot,trim(prefix),ow_checkpoint,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'overlapping-Wannier checkpoint publication failed';endif
  end subroutine

  subroutine periodic_box_gradients(values,box_size,gradient_coefficients,gradients)
    complex(8),intent(in)::values(:,:)
    integer,intent(in)::box_size(3)
    real(8),intent(in)::gradient_coefficients(:,:)
    complex(8),intent(out)::gradients(:,:,:)
    integer::i,j,k,p,pm,pp,axis,nstate,index(3),minus(3),plus(3),dist,radius
    nstate=size(values,1)
    radius=size(gradient_coefficients,1)
    do k=1,box_size(3);do j=1,box_size(2);do i=1,box_size(1)
      index=[i,j,k];p=i+box_size(1)*((j-1)+box_size(2)*(k-1))
      do axis=1,3
        gradients(axis,:,p)=(0d0,0d0)
        if(index(axis)>radius.and.index(axis)<=box_size(axis)-radius)then
          do dist=1,radius
            minus=index;plus=index;minus(axis)=index(axis)-dist;plus(axis)=index(axis)+dist
            pm=minus(1)+box_size(1)*((minus(2)-1)+box_size(2)*(minus(3)-1))
            pp=plus(1)+box_size(1)*((plus(2)-1)+box_size(2)*(plus(3)-1))
            gradients(axis,:,p)=gradients(axis,:,p)+gradient_coefficients(dist,axis)*&
              (values(:,pp)-values(:,pm))
          enddo
        else if(ow_core_size(axis)==dc%lg_tot%num(axis).and.ow_buffer(axis)==0)then
          do dist=1,radius
            minus=index;plus=index
            minus(axis)=modulo(index(axis)-dist-1,box_size(axis))+1
            plus(axis)=modulo(index(axis)+dist-1,box_size(axis))+1
            pm=minus(1)+box_size(1)*((minus(2)-1)+box_size(2)*(minus(3)-1))
            pp=plus(1)+box_size(1)*((plus(2)-1)+box_size(2)*(plus(3)-1))
            gradients(axis,:,p)=gradients(axis,:,p)+gradient_coefficients(dist,axis)*&
              (values(:,pp)-values(:,pm))
          enddo
        else
          minus=index;plus=index
          minus(axis)=max(1,index(axis)-1);plus(axis)=min(box_size(axis),index(axis)+1)
          pm=minus(1)+box_size(1)*((minus(2)-1)+box_size(2)*(minus(3)-1))
          pp=plus(1)+box_size(1)*((plus(2)-1)+box_size(2)*(plus(3)-1))
          gradients(axis,:,p)=(values(:,pp)-values(:,pm))/&
            (real(plus(axis)-minus(axis),8)*system%hgs(axis))
        endif
      enddo
    enddo;enddo;enddo
  end subroutine

  subroutine build_dc_translation_symmetry_map(nbox,symmetry_map,ok,message)
    integer,intent(in)::nbox
    integer(8),intent(out)::symmetry_map(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::rank_fragments(:)
    integer::source_fragment,operation,target_fragment,candidate_fragment,axis,p,match_count,ierr
    integer::translation(3),target_origin(3)
    logical::local_ok,global_ok

    allocate(rank_fragments(dc%n_frag))
    call MPI_Allgather(dc%i_frag,1,MPI_INTEGER,rank_fragments,1,MPI_INTEGER,dc%icomm_tot,ierr)
    local_ok=ierr==MPI_SUCCESS.and.nbox>0.and.size(symmetry_map,1)==nbox.and.&
      size(symmetry_map,2)==dc%n_frag.and.dc%n_frag==size(rank_fragments)
    do p=1,dc%n_frag
      local_ok=local_ok.and.count(rank_fragments==p)==1
    enddo
    local_ok=local_ok.and.all(dc%nxyz_domain_frag==spread(dc%nxyz_domain_frag(:,1),2,dc%n_frag))
    source_fragment=dc%i_frag
    if(local_ok)then
      do operation=1,dc%n_frag
        translation=dc%ixyz_frag(:,operation)-dc%ixyz_frag(:,1)
        target_origin=1+modulo(dc%ixyz_frag(:,source_fragment)-1+translation,dc%lg_tot%num)
        match_count=0;target_fragment=0
        do candidate_fragment=1,dc%n_frag
          if(all(1+modulo(dc%ixyz_frag(:,candidate_fragment)-1,dc%lg_tot%num)==target_origin))then
            match_count=match_count+1;target_fragment=candidate_fragment
          endif
        enddo
        if(match_count/=1)then
          local_ok=.false.;exit
        endif
        do p=1,nbox
          symmetry_map(p,operation)=int((target_fragment-1)*nbox+p,8)
        enddo
      enddo
    endif
    call comm_logical_and(local_ok,global_ok,dc%icomm_tot)
    ok=global_ok
    if(ok)then
      message=''
    else
      message='DC fragment origins do not form a uniform periodic translation group'
    endif
  end subroutine

  subroutine restore_ow_checkpoint_density(checkpoint,ok,message)
    type(s_dg_overlapping_wannier_checkpoint),intent(in)::checkpoint
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8),allocatable::local_density(:),global_density(:),density4(:,:,:,:)
    integer::p
    allocate(local_density(int(ow_global_grid_count)),global_density(int(ow_global_grid_count)))
    local_density=0d0
    do p=1,size(checkpoint%core_physical_ids)
      local_density(int(checkpoint%core_physical_ids(p)))=checkpoint%density(p)
    enddo
    call comm_summation(local_density,global_density,size(local_density),dc%icomm_tot)
    allocate(density4(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),1))
    density4(:,:,:,1)=reshape(global_density,dc%lg_tot%num)
    call dg_dc_update_potential_from_density(density4,ok,message)
  end subroutine

  subroutine ow_build_hamiltonian(comm,density,potential,hrows,new_potential,fingerprint,ok,message)
    integer,intent(in)::comm
    real(8),intent(in)::density(:),potential(:)
    complex(8),intent(out)::hrows(:,:)
    real(8),intent(out)::new_potential(:)
    integer(8),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8),allocatable::global_density(:),summed_density(:),core_potential(:)
    complex(8),allocatable::kinetic(:,:),local_matrix(:,:),nonlocal_matrix(:,:),global_h(:,:)
    integer::p,ix,iy,iz,nwann,owned_core,owned_projectors
    allocate(global_density(int(ow_global_grid_count)),summed_density(int(ow_global_grid_count)))
    global_density=0d0
    do p=1,size(ow_core_ids);global_density(int(ow_core_ids(p)))=density(p);enddo
    call comm_summation(global_density,summed_density,size(global_density),comm)
    if(allocated(ow_work_density))deallocate(ow_work_density)
    allocate(ow_work_density(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),1))
    ow_work_density(:,:,:,1)=reshape(summed_density,dc%lg_tot%num)
    call dg_dc_update_potential_from_density(ow_work_density,ok,message)
    if(.not.ok)return
    nwann=size(ow_box_values,1);allocate(core_potential(size(ow_core_ids)))
    p=0
    do iz=1,ow_core_size(3);do iy=1,ow_core_size(2);do ix=1,ow_core_size(1)
      p=p+1;core_potential(p)=v_local(1)%f(ix,iy,iz)
    enddo;enddo;enddo
    call assemble_dg_overlapping_wannier_weak_operators(comm,nwann,ow_core_ids,ow_core_weights,&
      ow_core_values,ow_core_gradients,core_potential,ow_global_grid_count,&
      kinetic,local_matrix,owned_core,ok,message)
    if(.not.ok)return
    call assemble_ow_nonlocal(comm,nonlocal_matrix,owned_projectors,ok,message)
    if(.not.ok)return
    allocate(global_h(nwann,nwann));global_h=kinetic+local_matrix+nonlocal_matrix
    if(maxval(abs(global_h-conjg(transpose(global_h))))>&
       dg_dc_gs_hermiticity_tolerance*max(1d0,maxval(abs(global_h))))then
      ok=.false.;message='weak overlapping-Wannier Hamiltonian is not Hermitian';return
    endif
    do p=1,size(ow_row_ids);hrows(p,:)=global_h(int(ow_row_ids(p)),:);enddo
    new_potential=core_potential
    fingerprint=dg_dc_operator_fingerprint(.false.)
    ok=all(ieee_is_finite(real(global_h))).and.all(ieee_is_finite(aimag(global_h)))
    if(ok)then;message='';else;message='nonfinite projected DC Hamiltonian';endif
  end subroutine

  subroutine assemble_ow_nonlocal(comm,matrix,ownership_count,ok,message)
    integer,intent(in)::comm
    complex(8),allocatable,intent(out)::matrix(:,:)
    integer,intent(out)::ownership_count
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::local_overlap(:,:),owned_overlap(:,:)
    integer,allocatable::projector_counts(:),displacements(:),local_atom_ids(:),local_ordinals(:),&
      all_atom_ids(:),all_ordinals(:),all_unique_ids(:)
    integer(8),allocatable::projector_ids(:)
    real(8),allocatable::strength(:)
    logical,allocatable::complete(:,:)
    logical::global_ok
    integer::ilma,ia,j,ix,iy,iz,p,q,nwann,ierr,rank,nproc,total_records,total_projectors,&
      canonical_index(3),nowned,global_index,ordinal

    nwann=size(ow_core_values,1)
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    allocate(projector_counts(nproc),displacements(nproc))
    call MPI_Allgather(ppg%Nlma,1,MPI_INTEGER,projector_counts,1,MPI_INTEGER,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.all(projector_counts>=0)
    total_records=0;displacements(1)=0
    do p=1,nproc
      if(total_records>huge(total_records)-projector_counts(p))ok=.false.
      if(p>1)displacements(p)=total_records
      if(ok)total_records=total_records+projector_counts(p)
    enddo
    call comm_logical_and(ok,global_ok,comm);ok=global_ok
    if(.not.ok)then;message='invalid fragment-local nonlocal projector namespace';return;endif
    if(total_records==0)then
      allocate(matrix(nwann,nwann));matrix=(0d0,0d0);ownership_count=0
      ok=.true.;message='';return
    endif
    allocate(local_overlap(nwann,ppg%Nlma),local_atom_ids(ppg%Nlma),local_ordinals(ppg%Nlma),&
      all_atom_ids(total_records),all_ordinals(total_records),all_unique_ids(total_records))
    local_overlap=(0d0,0d0);local_atom_ids=0;local_ordinals=0
    do ilma=1,ppg%Nlma
      ia=ppg%ia_tbl(ilma)
      call map_dc_atom_to_physical_atom(ia,local_atom_ids(ilma),ok)
      ordinal=count(ppg%ia_tbl(1:ilma)==ia);local_ordinals(ilma)=ordinal
      do j=1,ppg%mps(ia)
        ix=ppg%jxyz(1,j,ia);iy=ppg%jxyz(2,j,ia);iz=ppg%jxyz(3,j,ia)
        canonical_index=[dc_to_canonical_index(ix,ow_core_size(1),ow_buffer(1)),&
          dc_to_canonical_index(iy,ow_core_size(2),ow_buffer(2)),&
          dc_to_canonical_index(iz,ow_core_size(3),ow_buffer(3))]
        if(any(canonical_index<1).or.any(canonical_index>ow_box_size))then
          ok=.false.;cycle
        endif
        p=canonical_index(1)+ow_box_size(1)*((canonical_index(2)-1)+&
          ow_box_size(2)*(canonical_index(3)-1))
        local_overlap(:,ilma)=local_overlap(:,ilma)+ppg%uV(j,ilma)*ow_box_values(:,p)
      enddo
    enddo
    call comm_logical_and(ok,global_ok,comm);ok=global_ok
    if(.not.ok)then;message='cannot canonicalize complete physical atom/projector support';return;endif
    call MPI_Allgatherv(local_atom_ids,ppg%Nlma,MPI_INTEGER,all_atom_ids,projector_counts,&
      displacements,MPI_INTEGER,comm,ierr)
    call MPI_Allgatherv(local_ordinals,ppg%Nlma,MPI_INTEGER,all_ordinals,projector_counts,&
      displacements,MPI_INTEGER,comm,ierr)
    total_projectors=0
    do p=1,total_records
      all_unique_ids(p)=0
      do q=1,p-1
        if(all_atom_ids(q)==all_atom_ids(p).and.all_ordinals(q)==all_ordinals(p))then
          all_unique_ids(p)=all_unique_ids(q);exit
        endif
      enddo
      if(all_unique_ids(p)==0)then
        total_projectors=total_projectors+1;all_unique_ids(p)=total_projectors
      endif
    enddo
    nowned=0
    do ilma=1,ppg%Nlma
      global_index=displacements(rank+1)+ilma
      if(count(all_unique_ids(1:global_index)==all_unique_ids(global_index))==1)nowned=nowned+1
    enddo
    allocate(owned_overlap(nwann,nowned),projector_ids(nowned),strength(nowned),complete(nwann,nowned))
    p=0
    do ilma=1,ppg%Nlma
      global_index=displacements(rank+1)+ilma
      if(count(all_unique_ids(1:global_index)==all_unique_ids(global_index))/=1)cycle
      p=p+1;projector_ids(p)=int(all_unique_ids(global_index),8)
      strength(p)=system%hvol*ppg%rinv_uvu(ilma);owned_overlap(:,p)=local_overlap(:,ilma)
    enddo
    complete=.true.
    call assemble_dg_overlapping_wannier_nonlocal(comm,nwann,projector_ids,strength,owned_overlap,&
      complete,int(total_projectors,8),matrix,ownership_count,ok,message)
  end subroutine

  subroutine map_dc_atom_to_physical_atom(local_atom,physical_atom,ok)
    integer,intent(in)::local_atom
    integer,intent(out)::physical_atom
    logical,intent(out)::ok
    real(8)::position(3),delta(3),best,tolerance
    integer::atom,sx,sy,sz
    position=system%Rion(:,local_atom)+dc%rxyz_frag(:,dc%i_frag)
    best=huge(best);physical_atom=0
    do atom=1,dc%system_tot%nion
      if(system%kion(local_atom)/=dc%system_tot%kion(atom))cycle
      do sz=-1,1;do sy=-1,1;do sx=-1,1
        delta=position-dc%system_tot%Rion(:,atom)-sx*dc%system_tot%primitive_a(:,1)-&
          sy*dc%system_tot%primitive_a(:,2)-sz*dc%system_tot%primitive_a(:,3)
        if(sum(delta*delta)<best)then;best=sum(delta*delta);physical_atom=atom;endif
      enddo;enddo;enddo
    enddo
    tolerance=1d-16*max(1d0,maxval(abs(dc%system_tot%primitive_a)))**2
    ok=physical_atom>0.and.best<=tolerance
  end subroutine

  integer function dc_to_canonical_index(index,core_count,buffer_count)
    integer,intent(in)::index,core_count,buffer_count
    if(index<=core_count+buffer_count)then
      dc_to_canonical_index=buffer_count+index
    else
      dc_to_canonical_index=index-(core_count+buffer_count)
    endif
  end function

  subroutine ow_mix_density(comm,mixing_rate,current_density,raw_density,history,history_count,&
      mixed_density,new_history,new_history_count,ok,message)
    integer,intent(in)::comm,history_count
    real(8),intent(in)::mixing_rate,current_density(:),raw_density(:),history(:,:)
    real(8),intent(out)::mixed_density(:),new_history(:,:)
    integer,intent(out)::new_history_count
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::keep
    real(8),allocatable::accepted4(:,:,:,:),raw4(:,:,:,:),mixed4(:,:,:,:)
    allocate(accepted4(size(current_density),1,1,1),raw4(size(raw_density),1,1,1),&
      mixed4(size(mixed_density),1,1,1))
    accepted4(:,1,1,1)=current_density;raw4(:,1,1,1)=raw_density
    call mix_dg_dc_density_transaction(accepted4,raw4,mixing_rate,.false.,mixed4,ok,message)
    if(.not.ok)return
    mixed_density=mixed4(:,1,1,1)
    new_history=history;keep=min(history_count,size(new_history,2)-1)
    if(keep>0)new_history(:,1:keep)=history(:,history_count-keep+1:history_count)
    new_history(:,keep+1)=raw_density;new_history_count=keep+1
    ok=all(ieee_is_finite(mixed_density));message=''
  end subroutine

  subroutine ow_transaction(action,ok,message)
    integer,intent(in)::action
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::is
    ok=.true.;message=''
    select case(action)
    case(0)
      if(allocated(ow_density_snapshot))deallocate(ow_density_snapshot)
      allocate(ow_density_snapshot(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),system%nspin))
      do is=1,system%nspin;ow_density_snapshot(:,:,:,is)=dc%rho_tot_s(is)%f;enddo
      ow_potential_epoch_snapshot=dg_gs_potential_epoch
      ow_transaction_active=.true.
    case(-1)
      if(.not.ow_transaction_active)then;ok=.false.;message='missing overlapping-Wannier rollback snapshot';return;endif
      call dg_dc_update_potential_from_density(ow_density_snapshot,ok,message)
      if(ok)dg_gs_potential_epoch=ow_potential_epoch_snapshot
      ow_transaction_active=.false.
    case(1)
      ok=ow_transaction_active
      if(.not.ok)message='missing overlapping-Wannier prepared transaction'
    end select
  end subroutine

  subroutine ow_commit_transaction()
    ow_transaction_active=.false.
    if(allocated(ow_density_snapshot))deallocate(ow_density_snapshot)
  end subroutine

  subroutine populate_ow_checkpoint(occupations,condition_number,closure_residual)
    real(8),intent(in)::occupations(:),condition_number,closure_residual
    integer::rank,i,j,nowned,nbox,nproc,ierr
    integer(8)::tail_count8
    integer(8),allocatable::all_tail_ids(:)
    call MPI_Comm_rank(dc%icomm_tot,rank,i)
    call MPI_Comm_size(dc%icomm_tot,nproc,ierr)
    nowned=count(ow_basis%center_owner_rank==rank)
    if(size(ow_basis%physical_grid_ids)>huge(nbox)/nproc)&
      error stop 'overlapping-Wannier checkpoint box extent overflow'
    nbox=size(ow_basis%physical_grid_ids)*nproc
    tail_count8=int(nowned,8)*int(nbox,8)
    if(tail_count8>int(huge(nbox),8))error stop 'overlapping-Wannier checkpoint tail extent overflow'
    allocate(all_tail_ids(nbox))
    call MPI_Allgather(ow_basis%physical_grid_ids,size(ow_basis%physical_grid_ids),MPI_INTEGER8,&
      all_tail_ids,size(ow_basis%physical_grid_ids),MPI_INTEGER8,dc%icomm_tot,ierr)
    ow_checkpoint=s_dg_overlapping_wannier_checkpoint()
    ow_checkpoint%basis_generation=ow_basis%generation;ow_checkpoint%geometry_generation=1
    ow_checkpoint%basis_fingerprint=ow_state%basis_fingerprint
    ow_checkpoint%operator_fingerprint=ow_state%operator_fingerprint
    allocate(ow_checkpoint%center_owner,source=ow_basis%center_owner_rank)
    allocate(ow_checkpoint%overlap,source=ow_s)
    allocate(ow_checkpoint%coefficients,source=ow_state%coefficients)
    allocate(ow_checkpoint%occupations,source=occupations)
    allocate(ow_checkpoint%core_physical_ids,source=ow_core_ids)
    allocate(ow_checkpoint%density,source=ow_state%density)
    allocate(ow_checkpoint%tail_center(nowned),ow_checkpoint%tail_generation(nowned),&
      ow_checkpoint%tail_offsets(nowned+1),ow_checkpoint%tail_physical_ids(int(tail_count8)))
    i=0;ow_checkpoint%tail_offsets(1)=1
    do j=1,size(ow_basis%center_owner_rank)
      if(ow_basis%center_owner_rank(j)/=rank)cycle
      i=i+1;ow_checkpoint%tail_center(i)=j;ow_checkpoint%tail_generation(i)=ow_basis%generation
      ow_checkpoint%tail_offsets(i+1)=i*nbox+1
      ow_checkpoint%tail_physical_ids((i-1)*nbox+1:i*nbox)=all_tail_ids
    enddo
    ow_checkpoint%density_residual=ow_result%density_residual
    ow_checkpoint%unmixed_density_residual=ow_result%unmixed_density_residual
    ow_checkpoint%coefficient_residual=ow_result%coefficient_residual
    ow_checkpoint%orthogonality_defect=ow_result%orthogonality_defect
    ow_checkpoint%metric_condition=condition_number
    ow_checkpoint%charge_error=ow_result%integrated_charge-ow_result%trace_charge
    ow_checkpoint%density_tolerance=dg_dc_gs_final_density_tolerance
    ow_checkpoint%coefficient_tolerance=dg_dc_gs_final_orbital_tolerance
    ow_checkpoint%orthogonality_tolerance=10d0*dg_dc_gs_final_orbital_tolerance
    ow_checkpoint%charge_tolerance=dg_dc_gs_electron_count_tolerance
    ow_checkpoint%condition_limit=1d0/dg_dc_metric_rank_tolerance
    ow_checkpoint%symmetry_closure_residual=closure_residual
    ow_checkpoint%symmetry_tolerance=dg_ow_symmetry_tolerance
    ow_checkpoint%accepted=ow_result%converged
  end subroutine

  subroutine run_dg_dc_local_basis_ground_state_for_main()
    type(s_dg_dc_local_basis_layout) :: layout
    type(s_dg_dc_local_basis_production_state) :: candidate_state
    complex(8), allocatable :: raw_basis(:,:),orthonormal_basis(:,:),basis_transform(:,:)
    complex(8), allocatable :: full_raw_basis(:,:),full_fragment_basis(:,:)
    complex(8), allocatable :: basis(:,:,:,:),hbasis(:,:),volume_block(:,:),interface_rows(:,:)
    complex(8), allocatable :: hamiltonian_rows(:,:),overlap_rows(:,:),coefficient_rows(:,:)
    complex(8), allocatable :: accepted_coefficient_rows(:,:)
    real(8), allocatable :: eigenvalues(:),occupations(:),local_density(:),quadrature_weights(:)
    real(8), allocatable :: accepted_eigenvalues(:)
    real(8), allocatable :: lambda_history(:),lambda_steps(:)
    real(8), allocatable :: previous_density(:,:,:,:),trial_density(:,:,:,:),mixed_density(:,:,:,:)
    real(8), allocatable :: accepted_density(:,:,:,:)
    integer :: core_size(3),npoint,nraw,neffective,ix,iy,iz,io,point,global_row
    integer :: spsi_shape(7),sttpsi_shape(7),shpsi_shape(7)
    integer :: iteration,local_rank,local_first,allocation_status,full_point_count,rollback_count
    integer :: accepted_step_count,stage_iteration,eigensolver_iterations
    integer(8) :: iteration_operator_fingerprint
    real(8) :: density_residual_local,density_residual_global,mix_rate
    real(8) :: accepted_lambda,trial_lambda,lambda_step,previous_stage_residual,stage_residual
    real(8) :: orbital_tolerance,density_tolerance,orbital_residual
    real(8) :: projector_overlap,minimum_projector_seen,subspace_residual
    real(8) :: orthogonality_defect,hermiticity_defect,face_balance_defect
    logical :: ok,stage_ok,stage_converged,unmixed_gate,first_stage
    character(256) :: message

    ok=system%nspin==1 .and. system%if_real_orbital .and. allocated(spsi%rwf) .and. &
      allocated(sttpsi%rwf) .and. allocated(shpsi%rwf) .and. dc%isize_tot==dc%n_frag .and. &
      .not.dc%optimized_fragment_geometry
    call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
    if(.not.dg_handoff_ok) stop 'DG DC local-basis production topology is unsupported'
    spsi_shape=shape(spsi%rwf);sttpsi_shape=shape(sttpsi%rwf);shpsi_shape=shape(shpsi%rwf)
    core_size=dc%nxyz_domain_frag(:,dc%i_frag)
    ok=all(spsi_shape(1:3)==sttpsi_shape(1:3)) .and. &
      all(spsi_shape(1:3)==shpsi_shape(1:3)) .and. all(spsi_shape(1:3)>=core_size) .and. &
      spsi_shape(4)>=1 .and. sttpsi_shape(4)>=1 .and. shpsi_shape(4)>=1 .and. &
      spsi_shape(5)>=system%no .and. sttpsi_shape(5)>=system%no .and. shpsi_shape(5)>=system%no .and. &
      spsi_shape(6)>=1 .and. sttpsi_shape(6)>=1 .and. shpsi_shape(6)>=1 .and. &
      spsi_shape(7)>=1 .and. sttpsi_shape(7)>=1 .and. shpsi_shape(7)>=1
    call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
    if(.not.dg_handoff_ok) stop 'DG DC local-basis orbital shape mismatch collectively'
    npoint=product(core_size)
    nraw=system%no
    full_point_count=size(spsi%rwf,1)*size(spsi%rwf,2)*size(spsi%rwf,3)
    allocate(raw_basis(npoint,nraw),orthonormal_basis(npoint,nraw),basis_transform(nraw,nraw), &
      full_raw_basis(full_point_count,nraw),stat=allocation_status)
    ok=allocation_status==0
    call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
    if(.not.dg_handoff_ok) stop 'DG DC local-basis basis allocation failed collectively'
    point=0
    do iz=1,core_size(3); do iy=1,core_size(2); do ix=1,core_size(1)
      point=point+1
      do io=1,nraw
        raw_basis(point,io)=cmplx(spsi%rwf(ix,iy,iz,1,io,1,1),0d0,8)
      end do
    end do; end do; end do
    do io=1,nraw
      full_raw_basis(:,io)=reshape(cmplx(spsi%rwf(:,:,:,1,io,1,1),0d0,8),[full_point_count])
    end do
    call orthonormalize_dg_dc_fragment_core_basis(raw_basis,system%hvol,dg_dc_metric_rank_tolerance, &
      orthonormal_basis,basis_transform,neffective,ok,message)
    call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
    if(.not.dg_handoff_ok) then
      write(*,'(a,i0,2a)')'[DG-DC-LOCAL-BASIS] rank ',local_rank,' orthonormalization: ',trim(message)
      stop 'DG DC local-basis core orthonormalization failed collectively'
    end if
    allocate(basis(core_size(1),core_size(2),core_size(3),neffective), &
      full_fragment_basis(full_point_count,neffective),stat=allocation_status)
    ok=allocation_status==0
    call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
    if(.not.dg_handoff_ok) stop 'DG DC local-basis transformed basis allocation failed collectively'
    call transform_dg_dc_fragment_buffer_basis(full_raw_basis,basis_transform,neffective, &
      full_fragment_basis,ok,message)
    ok=ok .and. maxval(abs(aimag(full_fragment_basis)))<=1d-12*max(1d0,maxval(abs(full_fragment_basis)))
    call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
    if(.not.dg_handoff_ok) stop 'DG DC local-basis full buffer transform failed collectively'
    do io=1,neffective
      basis(:,:,:,io)=reshape(orthonormal_basis(:,io),core_size)
    end do
    call initialize_dg_dc_local_basis_layout(layout,dc%i_frag,neffective,dc%nstate_tot, &
      dg_dc_geometry_fingerprint(),dg_dc_operator_fingerprint(),dc%icomm_tot,ok,message)
    if(.not.ok) stop 'DG DC local-basis distributed layout failed'
    allocate(interface_rows(neffective,layout%global_basis_count),stat=allocation_status)
    ok=allocation_status==0
    call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
    if(.not.dg_handoff_ok) stop 'DG DC local-basis interface allocation failed collectively'
    call assemble_dg_dc_local_basis_interface_rows(layout,basis,dc%ixyz_frag,dc%nxyz_domain_frag, &
      dc%lg_tot%num,system%hgs,dg_dc_gs_sipg_penalty_factor,dc%icomm_tot,interface_rows,ok,message)
    call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
    if(.not.dg_handoff_ok) stop 'DG DC local-basis SIPG interface assembly failed collectively'
    call diagnose_dg_dc_six_face_balance(dc%ixyz_frag,dc%nxyz_domain_frag,dc%lg_tot%num, &
      dc%icomm_tot,face_balance_defect,ok,message)
    if(.not.ok) stop 'DG DC local-basis canonical face balance failed collectively'
    allocate(volume_block(neffective,neffective),hamiltonian_rows(neffective,layout%global_basis_count), &
      overlap_rows(neffective,layout%global_basis_count),coefficient_rows(neffective,layout%global_band_count), &
      eigenvalues(layout%global_band_count),occupations(layout%global_band_count),local_density(npoint), &
      quadrature_weights(npoint),stat=allocation_status)
    ok=allocation_status==0
    call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
    if(.not.dg_handoff_ok) stop 'DG DC local-basis SCF allocation failed collectively'
    quadrature_weights=system%hvol
    overlap_rows=(0d0,0d0)
    local_rank=dc%id_tot
    local_first=layout%basis_offsets(local_rank)+1
    do io=1,neffective
      global_row=local_first+io-1
      overlap_rows(io,global_row)=1d0
    end do
    call initialize_dg_dc_local_basis_coefficients(layout,coefficient_rows,ok,message)
    call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
    if(.not.dg_handoff_ok) stop 'DG DC local-basis coefficient initialization failed collectively'
    call assign_dg_dc_local_basis_occupations(dc%elec_num_tot,occupations,ok,message)
    call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
    if(.not.dg_handoff_ok) stop 'DG DC local-basis occupation initialization failed collectively'
    allocate(previous_density(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),1), &
      trial_density(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),1), &
      mixed_density(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),1), &
      accepted_density(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),1), &
      accepted_coefficient_rows(neffective,layout%global_band_count), &
      accepted_eigenvalues(layout%global_band_count), &
      lambda_history(dg_dc_gs_maximum_scf_iterations+dg_dc_gs_maximum_rollbacks+2), &
      lambda_steps(dg_dc_gs_maximum_scf_iterations+dg_dc_gs_maximum_rollbacks+2),stat=allocation_status)
    ok=allocation_status==0
    call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
    if(.not.dg_handoff_ok) stop 'DG DC local-basis density allocation failed collectively'
    previous_density(:,:,:,1)=dc%rho_tot_s(1)%f
    accepted_density=previous_density
    accepted_coefficient_rows=coefficient_rows
    accepted_eigenvalues=0d0
    lambda_history=0d0
    lambda_steps=0d0
    mix_rate=dg_dc_gs_density_mix_rate
    accepted_lambda=0d0
    lambda_step=dg_dc_gs_initial_lambda_step
    rollback_count=0
    accepted_step_count=0
    minimum_projector_seen=1d0
    iteration=0
    first_stage=.true.
    unmixed_gate=.false.
    dg_gs_potential_epoch=0_8
Continuation_Stages: do
      if(unmixed_gate) then
        trial_lambda=dg_dc_gs_target_lambda
      else if(first_stage) then
        trial_lambda=0d0
      else
        trial_lambda=min(dg_dc_gs_target_lambda,accepted_lambda+lambda_step)
      end if
      orbital_tolerance=merge(dg_dc_gs_final_orbital_tolerance, &
        dg_dc_gs_intermediate_orbital_tolerance,trial_lambda==dg_dc_gs_target_lambda)
      density_tolerance=merge(dg_dc_gs_final_density_tolerance, &
        dg_dc_gs_intermediate_density_tolerance,trial_lambda==dg_dc_gs_target_lambda)
      stage_converged=.false.
      stage_ok=.true.
      previous_stage_residual=huge(1d0)
      do stage_iteration=1,merge(1,dg_dc_gs_maximum_scf_iterations,unmixed_gate)
      iteration=iteration+1
      iteration_operator_fingerprint=dg_dc_operator_fingerprint()
      dg_gs_hamiltonian_potential_epoch=dg_gs_potential_epoch
      sttpsi%rwf=0d0
      do io=1,neffective
        sttpsi%rwf(:,:,:,1,io,1,1)=reshape(real(full_fragment_basis(:,io),8),shape(sttpsi%rwf(:,:,:,1,io,1,1)))
      end do
      call hpsi(sttpsi,shpsi,info,mg,v_local,system,stencil,srg,ppg)
      allocate(hbasis(npoint,neffective),stat=allocation_status)
      ok=allocation_status==0
      call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
      if(.not.dg_handoff_ok) stop 'DG DC local-basis H-basis allocation failed collectively'
      point=0
      do iz=1,core_size(3); do iy=1,core_size(2); do ix=1,core_size(1)
        point=point+1
        do io=1,neffective
          hbasis(point,io)=cmplx(shpsi%rwf(ix,iy,iz,1,io,1,1),0d0,8)
        end do
      end do; end do; end do
      call project_dg_dc_local_basis_volume(orthonormal_basis(:,1:neffective),hbasis,system%hvol, &
        volume_block,ok,message)
      deallocate(hbasis)
      call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
      if(.not.dg_handoff_ok) then; stage_ok=.false.; exit; end if
      call compose_dg_dc_distributed_hamiltonian_rows(layout,volume_block,interface_rows, &
        trial_lambda, &
        hamiltonian_rows,ok,message)
      call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
      if(.not.dg_handoff_ok) then; stage_ok=.false.; exit; end if
      call solve_dg_dc_local_basis_bands_cg(layout,hamiltonian_rows,overlap_rows,dc%icomm_tot, &
        dg_dc_gs_maximum_eigensolver_iterations,orbital_tolerance, &
        coefficient_rows,eigenvalues,ok,message,eigensolver_iterations,orbital_residual)
      call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
      if(.not.dg_handoff_ok) then; stage_ok=.false.; exit; end if
      call diagnose_dg_dc_local_basis_continuation(layout,hamiltonian_rows,coefficient_rows, &
        accepted_coefficient_rows,occupations,dc%icomm_tot,projector_overlap,orthogonality_defect, &
        hermiticity_defect,ok,message)
      subspace_residual=sqrt(max(0d0,1d0-projector_overlap))
      ok=ok .and. projector_overlap>=dg_dc_gs_minimum_projector_overlap .and. &
        subspace_residual<=dg_dc_gs_subspace_tolerance .and. &
        orthogonality_defect<=dg_dc_gs_orthogonality_tolerance .and. &
        hermiticity_defect<=dg_dc_gs_hermiticity_tolerance .and. &
        face_balance_defect<=dg_dc_gs_face_balance_tolerance .and. &
        layout%geometry_fingerprint==dg_dc_geometry_fingerprint() .and. &
        iteration_operator_fingerprint==dg_dc_operator_fingerprint()
      call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
      if(.not.dg_handoff_ok) then; stage_ok=.false.; exit; end if
      call reconstruct_dg_dc_local_basis_density(orthonormal_basis(:,1:neffective),coefficient_rows, &
        occupations,local_density,ok,message)
      call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
      if(.not.dg_handoff_ok) then; stage_ok=.false.; exit; end if
      call validate_dg_dc_local_basis_density(occupations,2d0,dc%elec_num_tot,local_density, &
        quadrature_weights,dc%icomm_tot,ok,message)
      if(.not.ok) then; stage_ok=.false.; exit; end if
      rho_s(1)%f=0d0
      point=0
      do iz=1,core_size(3); do iy=1,core_size(2); do ix=1,core_size(1)
        point=point+1
        rho_s(1)%f(ix,iy,iz)=local_density(point)
      end do; end do; end do
      call calc_rho_total_dcdft(system%nspin,lg,mg,info,rho_s,dc)
      trial_density(:,:,:,1)=dc%rho_tot_s(1)%f
      density_residual_local=sum((trial_density(:,:,:,1)-previous_density(:,:,:,1))**2)* &
        dc%system_tot%hvol/real(dc%isize_tot,8)
      call comm_summation(density_residual_local,density_residual_global,dc%icomm_tot)
      density_residual_global=sqrt(max(0d0,density_residual_global))
      if(unmixed_gate) then
        mixed_density=trial_density
      else
        mixed_density=(1d0-mix_rate)*previous_density+mix_rate*trial_density
      end if
      previous_density=mixed_density
      call dg_dc_update_potential_from_density(previous_density,ok,message)
      call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
      if(.not.dg_handoff_ok) then; stage_ok=.false.; exit; end if
      stage_residual=density_residual_global
      if(stage_iteration>1 .and. stage_residual>dg_dc_gs_allowed_residual_growth* &
        max(previous_stage_residual,epsilon(1d0))) then
        stage_ok=.false.; exit
      end if
      previous_stage_residual=stage_residual
      if(density_residual_global<=density_tolerance) then
        stage_converged=.true.; exit
      end if
      end do
      if(unmixed_gate) then
        if(.not.stage_ok .or. .not.stage_converged) then
          coefficient_rows=accepted_coefficient_rows
          eigenvalues=accepted_eigenvalues
          previous_density=accepted_density
          call dg_dc_update_potential_from_density(previous_density,ok,message)
          call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
          if(.not.dg_handoff_ok) stop 'DG DC local-basis unmixed rollback potential restore failed'
          stop 'DG DC local-basis unmixed fixed-point gate failed'
        end if
        exit Continuation_Stages
      end if
      if(.not.stage_ok .or. .not.stage_converged) then
        coefficient_rows=accepted_coefficient_rows
        eigenvalues=accepted_eigenvalues
        previous_density=accepted_density
        call dg_dc_update_potential_from_density(previous_density,ok,message)
        call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
        if(.not.dg_handoff_ok) stop 'DG DC local-basis rollback potential restore failed collectively'
        rollback_count=rollback_count+1
        lambda_step=0.5d0*lambda_step
        if(first_stage .or. rollback_count>dg_dc_gs_maximum_rollbacks .or. &
          lambda_step<dg_dc_gs_minimum_lambda_step) &
          stop 'DG DC local-basis continuation rollback limit reached'
        cycle Continuation_Stages
      end if
      accepted_step_count=accepted_step_count+1
      minimum_projector_seen=min(minimum_projector_seen,projector_overlap)
      lambda_history(accepted_step_count)=trial_lambda
      lambda_steps(accepted_step_count)=trial_lambda-accepted_lambda
      accepted_lambda=trial_lambda
      accepted_coefficient_rows=coefficient_rows
      accepted_eigenvalues=eigenvalues
      accepted_density=previous_density
      first_stage=.false.
      if(accepted_lambda==dg_dc_gs_target_lambda) then
        unmixed_gate=.true.
      else
        lambda_step=min(dg_dc_gs_maximum_lambda_step,1.5d0*lambda_step)
      end if
    end do Continuation_Stages
    if(comm_is_root(nproc_id_global)) write(*,'(a,i0,2(a,es24.16))') &
      '[DG-DC-LOCAL-BASIS] converged iterations=',iteration,' density_residual=',density_residual_global, &
      ' lambda=',dg_dc_gs_target_lambda
    allocate(candidate_state%coefficient_rows,source=coefficient_rows,stat=allocation_status)
    if(allocation_status==0) allocate(candidate_state%full_fragment_basis,source=full_fragment_basis, &
      stat=allocation_status)
    if(allocation_status==0) allocate(candidate_state%basis_transform, &
      source=basis_transform(:,1:neffective),stat=allocation_status)
    if(allocation_status==0) allocate(candidate_state%eigenvalues,source=eigenvalues,stat=allocation_status)
    if(allocation_status==0) allocate(candidate_state%occupations,source=occupations,stat=allocation_status)
    if(allocation_status==0) allocate(candidate_state%basis_offsets,source=layout%basis_offsets, &
      stat=allocation_status)
    if(allocation_status==0) allocate(candidate_state%fragment_ids,source=layout%fragment_ids, &
      stat=allocation_status)
    ok=allocation_status==0
    call comm_logical_and(ok,dg_handoff_ok,dc%icomm_tot)
    if(.not.dg_handoff_ok) stop 'DG DC local-basis result state allocation failed collectively'
    candidate_state%scf_iterations=iteration
    candidate_state%geometry_fingerprint=layout%geometry_fingerprint
    candidate_state%hamiltonian_operator_fingerprint=iteration_operator_fingerprint
    candidate_state%operator_fingerprint=dg_dc_operator_fingerprint()
    candidate_state%density_residual=density_residual_global
    candidate_state%interface_scale=dg_dc_gs_target_lambda
    candidate_state%fragment_id=layout%fragment_id
    candidate_state%local_basis_count=layout%local_basis_count
    candidate_state%global_basis_count=layout%global_basis_count
    candidate_state%global_band_count=layout%global_band_count
    candidate_state%raw_basis_count=nraw
    candidate_state%core_size=core_size
    candidate_state%full_spatial_shape=spsi_shape(1:3)
    candidate_state%ready=.true.
    if(allocated(dg_dc_local_basis_state%coefficient_rows)) deallocate(dg_dc_local_basis_state%coefficient_rows)
    if(allocated(dg_dc_local_basis_state%full_fragment_basis)) &
      deallocate(dg_dc_local_basis_state%full_fragment_basis)
    if(allocated(dg_dc_local_basis_state%basis_transform)) deallocate(dg_dc_local_basis_state%basis_transform)
    if(allocated(dg_dc_local_basis_state%eigenvalues)) deallocate(dg_dc_local_basis_state%eigenvalues)
    if(allocated(dg_dc_local_basis_state%occupations)) deallocate(dg_dc_local_basis_state%occupations)
    if(allocated(dg_dc_local_basis_state%basis_offsets)) deallocate(dg_dc_local_basis_state%basis_offsets)
    if(allocated(dg_dc_local_basis_state%fragment_ids)) deallocate(dg_dc_local_basis_state%fragment_ids)
    call move_alloc(candidate_state%coefficient_rows,dg_dc_local_basis_state%coefficient_rows)
    call move_alloc(candidate_state%full_fragment_basis,dg_dc_local_basis_state%full_fragment_basis)
    call move_alloc(candidate_state%basis_transform,dg_dc_local_basis_state%basis_transform)
    call move_alloc(candidate_state%eigenvalues,dg_dc_local_basis_state%eigenvalues)
    call move_alloc(candidate_state%occupations,dg_dc_local_basis_state%occupations)
    call move_alloc(candidate_state%basis_offsets,dg_dc_local_basis_state%basis_offsets)
    call move_alloc(candidate_state%fragment_ids,dg_dc_local_basis_state%fragment_ids)
    dg_dc_local_basis_state%scf_iterations=candidate_state%scf_iterations
    dg_dc_local_basis_state%geometry_fingerprint=candidate_state%geometry_fingerprint
    dg_dc_local_basis_state%operator_fingerprint=candidate_state%operator_fingerprint
    dg_dc_local_basis_state%hamiltonian_operator_fingerprint= &
      candidate_state%hamiltonian_operator_fingerprint
    dg_dc_local_basis_state%density_residual=candidate_state%density_residual
    dg_dc_local_basis_state%interface_scale=candidate_state%interface_scale
    dg_dc_local_basis_state%fragment_id=candidate_state%fragment_id
    dg_dc_local_basis_state%local_basis_count=candidate_state%local_basis_count
    dg_dc_local_basis_state%global_basis_count=candidate_state%global_basis_count
    dg_dc_local_basis_state%global_band_count=candidate_state%global_band_count
    dg_dc_local_basis_state%raw_basis_count=candidate_state%raw_basis_count
    dg_dc_local_basis_state%core_size=candidate_state%core_size
    dg_dc_local_basis_state%full_spatial_shape=candidate_state%full_spatial_shape
    dg_dc_local_basis_state%ready=.true.
    dg_gs_controls=default_dg_dc_gs_controls()
    dg_gs_controls%intermediate_orbital_tolerance=dg_dc_gs_intermediate_orbital_tolerance
    dg_gs_controls%intermediate_density_tolerance=dg_dc_gs_intermediate_density_tolerance
    dg_gs_controls%final_orbital_tolerance=dg_dc_gs_final_orbital_tolerance
    dg_gs_controls%final_density_tolerance=dg_dc_gs_final_density_tolerance
    dg_gs_controls%subspace_tolerance=dg_dc_gs_subspace_tolerance
    dg_gs_controls%initial_lambda_step=dg_dc_gs_initial_lambda_step
    dg_gs_controls%minimum_lambda_step=dg_dc_gs_minimum_lambda_step
    dg_gs_controls%maximum_lambda_step=dg_dc_gs_maximum_lambda_step
    dg_gs_controls%allowed_residual_growth=dg_dc_gs_allowed_residual_growth
    dg_gs_controls%density_mix_rate=dg_dc_gs_density_mix_rate
    dg_gs_controls%hermiticity_tolerance=dg_dc_gs_hermiticity_tolerance
    dg_gs_controls%orthogonality_tolerance=dg_dc_gs_orthogonality_tolerance
    dg_gs_controls%face_balance_tolerance=dg_dc_gs_face_balance_tolerance
    dg_gs_controls%electron_count_tolerance=dg_dc_gs_electron_count_tolerance
    dg_gs_controls%minimum_projector_overlap=dg_dc_gs_minimum_projector_overlap
    dg_gs_controls%maximum_scf_iterations=dg_dc_gs_maximum_scf_iterations
    dg_gs_controls%maximum_eigensolver_iterations=dg_dc_gs_maximum_eigensolver_iterations
    dg_gs_controls%maximum_rollbacks=dg_dc_gs_maximum_rollbacks
    if(allocated(dg_gs_result%lambda_history))deallocate(dg_gs_result%lambda_history)
    if(allocated(dg_gs_result%lambda_steps))deallocate(dg_gs_result%lambda_steps)
    allocate(dg_gs_result%lambda_history(accepted_step_count),dg_gs_result%lambda_steps(accepted_step_count))
    dg_gs_result%lambda_history=lambda_history(:accepted_step_count)
    dg_gs_result%lambda_steps=lambda_steps(:accepted_step_count)
    dg_gs_result%phase=DG_DC_ACCEPTED
    dg_gs_result%naccepted_steps=accepted_step_count
    dg_gs_result%nrollbacks=rollback_count
    dg_gs_result%mixing_reset_count=1
    dg_gs_result%total_scf_iterations=iteration
    dg_gs_result%lambda=accepted_lambda
    dg_gs_result%maximum_interface_scale=accepted_lambda
    dg_gs_result%minimum_projector_overlap=min(minimum_projector_seen,projector_overlap)
    dg_gs_result%accepted=.true.
    dg_gs_result%failed=.false.
    dg_gs_result%unmixed_verified=.true.
    dg_gs_result%final_diagnostics%orbital_residual=orbital_residual
    dg_gs_result%final_diagnostics%density_residual=density_residual_global
    dg_gs_result%final_diagnostics%subspace_residual=subspace_residual
    dg_gs_result%final_diagnostics%projector_overlap=projector_overlap
    dg_gs_result%final_diagnostics%hermiticity_defect=hermiticity_defect
    dg_gs_result%final_diagnostics%orthogonality_defect=orthogonality_defect
    dg_gs_result%final_diagnostics%face_balance_defect=face_balance_defect
    dg_gs_result%final_diagnostics%electron_number=sum(occupations)
    dg_gs_result%final_diagnostics%expected_electron_number=dc%elec_num_tot
    dg_gs_result%final_diagnostics%interface_scale=accepted_lambda
    dg_gs_result%final_diagnostics%eigensolver_iterations=eigensolver_iterations
    dg_gs_result%final_diagnostics%hamiltonian_potential_epoch=dg_gs_hamiltonian_potential_epoch
    dg_gs_result%final_diagnostics%updated_potential_epoch=dg_gs_potential_epoch
    dg_gs_result%final_diagnostics%eigensolver_converged=.true.
    dg_gs_result%final_diagnostics%finite=.true.
    call publish_dg_dc_ground_state_for_main(dg_handoff_ok,dg_handoff_message)
    if(.not.dg_handoff_ok) stop 'DG DC local-basis checkpoint publication failed collectively'
  end subroutine run_dg_dc_local_basis_ground_state_for_main

  subroutine materialize_dg_dc_candidates_for_main()
    complex(8), allocatable :: candidate_box(:,:,:,:,:)
    integer(8), allocatable :: owner(:,:,:)
    integer :: box_size(3),orbital_storage_size(3),core_size(3),buffer(3),maximum_candidate_count
    integer :: ix,iy,iz,io,is,ix_tot,iy_tot,iz_tot,raw_ix,raw_iy,raw_iz
    integer(8) :: geometry_fingerprint,operator_fingerprint
    logical :: local_preflight,global_preflight

    local_preflight=dc%isize_tot==dc%n_frag .and. (allocated(spsi%rwf) .or. allocated(spsi%zwf))
    call comm_logical_and(local_preflight,global_preflight,dc%icomm_tot)
    if(.not.global_preflight) then
      if(comm_is_root(nproc_id_global)) write(*,'(a,2(a,i0),2(a,l1))')'[DG-DC-GS] preflight failure', &
        ' ranks=',dc%isize_tot,' fragments=',dc%n_frag,' rwf=',allocated(spsi%rwf),' zwf=',allocated(spsi%zwf)
      stop 'DG DC local-periodic handoff topology/orbital preflight failed collectively'
    end if
    buffer=dc%nxyz_buffer
    core_size=dc%nxyz_domain_frag(:,dc%i_frag)
    if(allocated(spsi%rwf)) then
      orbital_storage_size=[size(spsi%rwf,1),size(spsi%rwf,2),size(spsi%rwf,3)]
    else
      orbital_storage_size=[size(spsi%zwf,1),size(spsi%zwf,2),size(spsi%zwf,3)]
    end if
    box_size=core_size+2*buffer
    local_preflight=all(orbital_storage_size>=box_size) .and. system%no>0
    call comm_logical_and(local_preflight,global_preflight,dc%icomm_tot)
    if(.not.global_preflight) then
      if(comm_is_root(nproc_id_global)) write(*,'(a,4(a,3(i0,1x)),3(a,i0))') &
        '[DG-DC-GS] candidate preflight failure',' storage=',orbital_storage_size,' core=',core_size, &
        ' buffer=',buffer,' physical_box=',box_size,' local_states=',system%no,' atoms=',natom, &
        ' candidates_per_atom=',dg_dc_candidate_orbitals_per_atom
      stop 'DG DC handoff candidate/geometry preflight failed collectively'
    end if
    maximum_candidate_count=system%no
    call comm_get_max(maximum_candidate_count,dc%icomm_tot)
    allocate(candidate_box(box_size(1),box_size(2),box_size(3),maximum_candidate_count,system%nspin))
    candidate_box=(0d0,0d0)
    do is=1,system%nspin
    do io=1,system%no
    do iz=1,box_size(3)
      raw_iz=canonical_to_dc_index(iz,core_size(3),buffer(3))
    do iy=1,box_size(2)
      raw_iy=canonical_to_dc_index(iy,core_size(2),buffer(2))
    do ix=1,box_size(1)
      raw_ix=canonical_to_dc_index(ix,core_size(1),buffer(1))
      if(allocated(spsi%rwf)) then
        candidate_box(ix,iy,iz,io,is)=cmplx(spsi%rwf(raw_ix,raw_iy,raw_iz,is,io,1,1),0d0,8)
      else
        candidate_box(ix,iy,iz,io,is)=spsi%zwf(raw_ix,raw_iy,raw_iz,is,io,1,1)
      end if
    end do
    end do
    end do
    end do
    end do
    allocate(owner(core_size(1),core_size(2),core_size(3)))
    do iz=1,core_size(3)
      iz_tot=dc%jxyz_tot(iz,3)
    do iy=1,core_size(2)
      iy_tot=dc%jxyz_tot(iy,2)
    do ix=1,core_size(1)
      ix_tot=dc%jxyz_tot(ix,1)
      owner(ix,iy,iz)=int(ix_tot,8)+int(dc%lg_tot%num(1),8)* &
        (int(iy_tot-1,8)+int(dc%lg_tot%num(2),8)*int(iz_tot-1,8))
    end do
    end do
    end do
    geometry_fingerprint=dg_dc_geometry_fingerprint()
    operator_fingerprint=dg_dc_operator_fingerprint()
    call materialize_dg_dc_candidates(dg_dc_handoff_runtime,dg_dc_nodal_runtime,candidate_box,dc%i_frag, &
      core_size,buffer,owner, &
      natom,system%no,geometry_fingerprint,operator_fingerprint,dc%icomm_tot,dg_handoff_ok, &
      dg_handoff_message)
    if(.not.dg_handoff_ok) then
      if(comm_is_root(nproc_id_global)) write(*,'(a,1x,a)')'[DG-DC-GS] materialization failure', &
        trim(dg_handoff_message)
      stop 'DG DC candidate materialization failed'
    end if
  end subroutine materialize_dg_dc_candidates_for_main

  subroutine run_dg_dc_ground_state_for_main()
    integer :: is,local_candidate_count
    local_candidate_count=dg_dc_handoff_runtime%candidate_count
    if(system%nspin/=1 .or. .not.system%if_real_orbital) &
      stop 'DG DC local-periodic production requires non-spin-polarized Gamma real orbitals'
    do is=1,system%nspin
      dg_dc_nodal_runtime%occupations(1:local_candidate_count,is)= &
        system%rocc(1:local_candidate_count,info%ik_s,is)
    end do
    call expand_dg_dc_global_candidate_axis(dg_dc_nodal_runtime,local_candidate_count,dc%icomm_tot, &
      dg_handoff_ok,dg_handoff_message)
    if(.not.dg_handoff_ok) then
      if(comm_is_root(nproc_id_global)) write(*,'(a,1x,a)')'[DG-DC-GS] global layout failure',trim(dg_handoff_message)
      stop 'DG DC global candidate layout failed collectively'
    end if
    allocate(dg_gs_fragment_origins(3,dc%n_frag),dg_gs_fragment_sizes(3,dc%n_frag))
    dg_gs_fragment_origins=dc%ixyz_frag
    dg_gs_fragment_sizes=dc%nxyz_domain_frag
    call initialize_dg_dc_physical_faces(dg_dc_nodal_runtime,dg_gs_fragment_origins,dg_gs_fragment_sizes, &
      dc%lg_tot%num,dc%icomm_tot,dg_handoff_ok,dg_handoff_message)
    if(.not.dg_handoff_ok) then
      if(comm_is_root(nproc_id_global)) write(*,'(a,1x,a)')'[DG-DC-GS] face topology failure', &
        trim(dg_handoff_message)
      stop 'DG DC physical face topology failed collectively'
    end if
    allocate(dg_gs_density(size(dc%rho_tot_s(1)%f,1),size(dc%rho_tot_s(1)%f,2), &
      size(dc%rho_tot_s(1)%f,3),system%nspin))
    do is=1,system%nspin
      dg_gs_density(:,:,:,is)=dc%rho_tot_s(is)%f
    end do
    allocate(dg_gs_raw_density,mold=dg_gs_density)
    allocate(dg_gs_mixed_density,mold=dg_gs_density)
    allocate(dg_gs_core_density(dg_dc_nodal_runtime%core_size(1),dg_dc_nodal_runtime%core_size(2), &
      dg_dc_nodal_runtime%core_size(3),system%nspin))
    allocate(dg_gs_eigenvalues(dg_dc_nodal_runtime%nstate,system%nspin))
    dg_gs_controls=default_dg_dc_gs_controls()
    dg_gs_potential_epoch=0_8
    dg_gs_controls%intermediate_orbital_tolerance=dg_dc_gs_intermediate_orbital_tolerance
    dg_gs_controls%intermediate_density_tolerance=dg_dc_gs_intermediate_density_tolerance
    dg_gs_controls%final_orbital_tolerance=dg_dc_gs_final_orbital_tolerance
    dg_gs_controls%final_density_tolerance=dg_dc_gs_final_density_tolerance
    dg_gs_controls%subspace_tolerance=dg_dc_gs_subspace_tolerance
    dg_gs_controls%initial_lambda_step=dg_dc_gs_initial_lambda_step
    dg_gs_controls%minimum_lambda_step=dg_dc_gs_minimum_lambda_step
    dg_gs_controls%maximum_lambda_step=dg_dc_gs_maximum_lambda_step
    dg_gs_controls%allowed_residual_growth=dg_dc_gs_allowed_residual_growth
    dg_gs_controls%density_mix_rate=dg_dc_gs_density_mix_rate
    dg_gs_controls%hermiticity_tolerance=dg_dc_gs_hermiticity_tolerance
    dg_gs_controls%orthogonality_tolerance=dg_dc_gs_orthogonality_tolerance
    dg_gs_controls%face_balance_tolerance=dg_dc_gs_face_balance_tolerance
    dg_gs_controls%electron_count_tolerance=dg_dc_gs_electron_count_tolerance
    dg_gs_controls%minimum_projector_overlap=dg_dc_gs_minimum_projector_overlap
    dg_gs_controls%maximum_scf_iterations=dg_dc_gs_maximum_scf_iterations
    dg_gs_controls%maximum_eigensolver_iterations=dg_dc_gs_maximum_eigensolver_iterations
    dg_gs_controls%maximum_rollbacks=dg_dc_gs_maximum_rollbacks
    call run_dg_dc_ground_state(dg_dc_nodal_runtime,dg_gs_density,dg_gs_controls, &
      dg_dc_salmon_scf_step,dc%icomm_tot,dg_gs_result,dg_handoff_ok,dg_handoff_message)
    if(.not.dg_handoff_ok .or. .not.dg_gs_result%accepted .or. dg_gs_result%lambda/=1d0 .or. &
      .not.dg_gs_result%unmixed_verified) then
      call dg_dc_update_potential_from_density(dg_gs_density,dg_handoff_ok,dg_handoff_message)
      if(.not.dg_handoff_ok) stop 'DG DC accepted potential rollback failed collectively'
      stop 'DG DC ground-state continuation failed collectively'
    end if
    call publish_dg_dc_ground_state_for_main(dg_handoff_ok,dg_handoff_message)
    if(.not.dg_handoff_ok) stop 'DG DC verified ground-state checkpoint publication failed collectively'
  end subroutine run_dg_dc_ground_state_for_main

  subroutine publish_dg_dc_direct_ground_state_for_main(ok,message)
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_dg_dc_direct_checkpoint_state)::checkpoint
    integer::neighbors(3,2),periodic_shifts(3,2),axis,side_index,iface,is
    integer::orbital_core_start(3),orbital_core_end(3)
    character(512)::checkpoint_root
    call build_dg_dc_six_face_neighbors(dc%i_frag,dc%ixyz_frag,dc%nxyz_domain_frag,dc%lg_tot%num, &
      neighbors,periodic_shifts,ok,message)
    if(.not.ok)return
    checkpoint%ready=.true.
    checkpoint%accepted=dg_dc_handoff_runtime%direct_ground_state_complete
    checkpoint%unmixed_verified=dg_dc_handoff_runtime%direct_ground_state_complete
    checkpoint%fragment_id=dc%i_frag
    checkpoint%nstate=system%no
    checkpoint%state_start=info%io_s
    checkpoint%state_count=info%io_e-info%io_s+1
    checkpoint%nspin=system%nspin
    checkpoint%scf_iterations=Miter
    checkpoint%core_size=dc%nxyz_domain_frag(:,dc%i_frag)
    checkpoint%full_spatial_shape=[size(spsi%rwf,1),size(spsi%rwf,2),size(spsi%rwf,3)]
    checkpoint%orbital_spatial_lower_bound=[lbound(spsi%rwf,1),lbound(spsi%rwf,2),lbound(spsi%rwf,3)]
    checkpoint%global_size=dc%lg_tot%num
    checkpoint%fragment_origin=dc%ixyz_frag(:,dc%i_frag)
    checkpoint%fragment_size=checkpoint%core_size
    checkpoint%density_origin=dc%mg_tot%is-dc%lg_tot%is
    checkpoint%density_size=dc%mg_tot%ie-dc%mg_tot%is+1
    orbital_core_start=max(mg%is,[1,1,1])
    orbital_core_end=min(mg%ie,checkpoint%core_size)
    checkpoint%orbital_core_local_origin=orbital_core_start
    checkpoint%orbital_core_size=max(0,orbital_core_end-orbital_core_start+1)
    checkpoint%orbital_core_origin=checkpoint%fragment_origin+min(orbital_core_start-1,checkpoint%core_size)
    do axis=1,3
    do side_index=1,2
      iface=2*(axis-1)+side_index
      checkpoint%face_neighbors(iface)=neighbors(axis,side_index)
    enddo
    enddo
    checkpoint%continuation_rollbacks=dg_dc_handoff_runtime%continuation_rollbacks
    checkpoint%projector_generation=dg_dc_handoff_runtime%projector_generation
    checkpoint%projector_retained_rank=dg_dc_handoff_runtime%projector_retained_rank
    checkpoint%projector_required_rank=dg_dc_handoff_runtime%projector_required_rank
    checkpoint%geometry_fingerprint=dg_dc_geometry_fingerprint()
    checkpoint%operator_fingerprint=dg_dc_operator_fingerprint()
    checkpoint%hamiltonian_operator_fingerprint=ieor(checkpoint%operator_fingerprint, &
      dg_dc_frozen_operator_fingerprint)
    if(checkpoint%hamiltonian_operator_fingerprint==0_8)checkpoint%hamiltonian_operator_fingerprint=1_8
    checkpoint%projector_fingerprint=dg_dc_handoff_runtime%projector_fingerprint
    checkpoint%density_residual=dg_dc_handoff_runtime%accepted_stage_residual
    checkpoint%interface_scale=dg_dc_handoff_runtime%accepted_interface_scale
    checkpoint%orbital_residual=dg_dc_handoff_runtime%direct_orbital_residual
    checkpoint%orthogonality_defect=dg_dc_handoff_runtime%direct_orthogonality_defect
    checkpoint%hermiticity_defect=dg_dc_handoff_runtime%direct_hermiticity_defect
    checkpoint%charge_error=dg_dc_handoff_runtime%charge_error
    checkpoint%face_balance_defect=dg_dc_handoff_runtime%direct_face_balance_defect
    checkpoint%density_tolerance=dg_dc_gs_final_density_tolerance
    checkpoint%orbital_tolerance=dg_dc_gs_final_orbital_tolerance
    checkpoint%orthogonality_tolerance=dg_dc_gs_orthogonality_tolerance
    checkpoint%hermiticity_tolerance=dg_dc_gs_hermiticity_tolerance
    checkpoint%charge_tolerance=dg_dc_gs_electron_count_tolerance
    checkpoint%face_balance_tolerance=dg_dc_gs_face_balance_tolerance
    checkpoint%projector_projection_residual=dg_dc_handoff_runtime%projector_projection_residual
    checkpoint%projector_escape_norm=dg_dc_handoff_runtime%projector_escape_norm
    checkpoint%projector_residual_limit=dg_dc_handoff_runtime%projector_residual_limit
    checkpoint%projector_escape_limit=dg_dc_handoff_runtime%projector_escape_limit
    checkpoint%face_action_norm=dg_dc_handoff_runtime%face_action_norm
    checkpoint%face_pair_balance=dg_dc_handoff_runtime%face_pair_balance
    checkpoint%fragment_orbitals=spsi%rwf
    checkpoint%occupations=system%rocc(info%io_s:info%io_e,:,:)
    checkpoint%continuation_scale_history=dg_dc_handoff_runtime%continuation_scale_history
    checkpoint%continuation_step_history=dg_dc_handoff_runtime%continuation_step_history
    allocate(checkpoint%density(checkpoint%density_size(1),checkpoint%density_size(2), &
      checkpoint%density_size(3),system%nspin))
    allocate(checkpoint%hartree(checkpoint%density_size(1),checkpoint%density_size(2), &
      checkpoint%density_size(3)))
    allocate(checkpoint%vxc,mold=checkpoint%density)
    allocate(checkpoint%vlocal,mold=checkpoint%density)
    checkpoint%hartree=dc%Vh_tot%f
    do is=1,system%nspin
      checkpoint%density(:,:,:,is)=dc%rho_tot_s(is)%f
      checkpoint%vxc(:,:,:,is)=dc%Vxc_tot(is)%f
      checkpoint%vlocal(:,:,:,is)=dc%vloc_tot(is)%f
    enddo
    checkpoint_root=trim(sysname)//'.dg_dc_direct'
    call publish_dg_dc_direct_checkpoint(trim(checkpoint_root),checkpoint,dc%icomm_tot,ok,message)
    if(ok.and.comm_is_root(nproc_id_global))then
      write(*,'(a,1x,a)')'[DG-DC-DIRECT] checkpoint',trim(checkpoint_root)
      write(*,'(a,2(a,es24.16),a,i0)')'[DG-DC-DIRECT] diagnostics', &
        ' orbital_residual=',checkpoint%orbital_residual, &
        ' density_residual=',checkpoint%density_residual, &
        ' rollbacks=',checkpoint%continuation_rollbacks
    endif
  end subroutine publish_dg_dc_direct_ground_state_for_main

  subroutine publish_dg_dc_ground_state_for_main(ok,message)
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_dg_ground_state_checkpoint)::checkpoint
    real(8),allocatable::checkpoint_density(:,:,:,:),checkpoint_vxc(:,:,:,:),checkpoint_vlocal(:,:,:,:)
    integer::is,axis,side_index,iface
    integer::face_neighbors(6)
    integer::neighbors(3,2),periodic_shifts(3,2)
    character(512)::checkpoint_root
    allocate(checkpoint_density(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),system%nspin))
    allocate(checkpoint_vxc,mold=checkpoint_density)
    allocate(checkpoint_vlocal,mold=checkpoint_vxc)
    do is=1,system%nspin
      checkpoint_density(:,:,:,is)=dc%rho_tot_s(is)%f
      checkpoint_vxc(:,:,:,is)=dc%Vxc_tot(is)%f
      checkpoint_vlocal(:,:,:,is)=dc%vloc_tot(is)%f
    end do
    call build_dg_dc_six_face_neighbors(dc%i_frag,dc%ixyz_frag,dc%nxyz_domain_frag,dc%lg_tot%num, &
      neighbors,periodic_shifts,ok,message)
    if(.not.ok)return
    do axis=1,3
    do side_index=1,2
      iface=2*(axis-1)+side_index
      face_neighbors(iface)=neighbors(axis,side_index)-1
    end do
    end do
    call populate_dg_ground_state_checkpoint(checkpoint,dg_dc_local_basis_state,dg_gs_result,dg_gs_controls, &
      checkpoint_density,dc%Vh_tot%f,checkpoint_vxc,checkpoint_vlocal,dc%lg_tot%num,dc%ixyz_frag, &
      dc%nxyz_domain_frag,face_neighbors,'DG_DC_GS',ok,message)
    if(.not.ok)return
    checkpoint_root=trim(sysname)//'.dg_dc_gs'
    call publish_dg_ground_state_checkpoint(trim(checkpoint_root),checkpoint,dc%icomm_tot,ok,message)
    if(ok .and. comm_is_root(nproc_id_global)) call report_dg_dc_ground_state_for_main(checkpoint_root)
  end subroutine publish_dg_dc_ground_state_for_main

  subroutine report_dg_dc_ground_state_for_main(checkpoint_root)
    character(*),intent(in)::checkpoint_root
    integer::ihistory
    write(*,'(a,1x,a)')'[DG-DC-GS] checkpoint',trim(checkpoint_root)//'.manifest'
    write(*,'(a,5(a,i0),3(a,es24.16))')'[DG-DC-GS] continuation', &
      ' accepted_steps=',dg_gs_result%naccepted_steps,' rollbacks=',dg_gs_result%nrollbacks, &
      ' mixing_resets=',dg_gs_result%mixing_reset_count,' scf_iterations=',dg_gs_result%total_scf_iterations, &
      ' metric_rank=',dg_dc_handoff_runtime%metric_rank, &
      ' lambda=',dg_gs_result%lambda,' projector_min=',dg_gs_result%minimum_projector_overlap, &
      ' dc_handoff_energy=',energy%E_tot
    do ihistory=1,dg_gs_result%naccepted_steps
      write(*,'(a,i0,2(a,es24.16))')'[DG-DC-GS] lambda_history step=',ihistory, &
        ' lambda=',dg_gs_result%lambda_history(ihistory),' delta=',dg_gs_result%lambda_steps(ihistory)
    end do
    write(*,'(a,8(a,es24.16),2(a,i0))')'[DG-DC-GS] acceptance', &
      ' orbital=',dg_gs_result%final_diagnostics%orbital_residual, &
      ' density=',dg_gs_result%final_diagnostics%density_residual, &
      ' subspace=',dg_gs_result%final_diagnostics%subspace_residual, &
      ' charge=',dg_gs_result%final_diagnostics%electron_number, &
      ' orthogonality=',dg_gs_result%final_diagnostics%orthogonality_defect, &
      ' hermiticity=',dg_gs_result%final_diagnostics%hermiticity_defect, &
      ' face_balance=',dg_gs_result%final_diagnostics%face_balance_defect, &
      ' fixed_point_lambda=',dg_gs_result%final_diagnostics%interface_scale, &
      ' eigensolver_iterations=',dg_gs_result%final_diagnostics%eigensolver_iterations, &
      ' potential_epoch=',dg_gs_result%final_diagnostics%updated_potential_epoch
  end subroutine report_dg_dc_ground_state_for_main

  subroutine dg_dc_salmon_scf_step(state,density_arg,lambda,density_mix,reset_mixing,unmixed, &
      diagnostics,communicator,step_ok,step_message)
    type(s_dg_nodal_common_state), intent(inout) :: state
    real(8), intent(inout) :: density_arg(:,:,:,:)
    real(8), intent(in) :: lambda,density_mix
    logical, intent(in) :: reset_mixing,unmixed
    type(s_dg_dc_gs_diagnostics), intent(out) :: diagnostics
    integer, intent(in) :: communicator
    logical, intent(out) :: step_ok
    character(*), intent(out) :: step_message
    call execute_dg_dc_production_iteration(state,density_arg,lambda,density_mix,reset_mixing,unmixed, &
      dg_dc_salmon_restore_step,dg_dc_salmon_solve_step,dg_dc_salmon_update_step,communicator, &
      diagnostics,step_ok,step_message)
  end subroutine dg_dc_salmon_scf_step

  subroutine dg_dc_salmon_restore_step(state,density_arg,communicator,step_ok,step_message)
    type(s_dg_nodal_common_state), intent(inout) :: state
    real(8), intent(inout) :: density_arg(:,:,:,:)
    integer, intent(in) :: communicator
    logical, intent(out) :: step_ok
    character(*), intent(out) :: step_message
    call dg_dc_update_potential_from_density(density_arg,step_ok,step_message)
  end subroutine dg_dc_salmon_restore_step

  subroutine dg_dc_salmon_solve_step(state,lambda,communicator,diagnostics,step_ok,step_message)
    type(s_dg_nodal_common_state), intent(inout) :: state
    real(8), intent(in) :: lambda
    integer, intent(in) :: communicator
    type(s_dg_dc_gs_diagnostics), intent(inout) :: diagnostics
    logical, intent(out) :: step_ok
    character(*), intent(out) :: step_message
    real(8) :: orbital_residual
    integer :: eigensolver_iterations

    dg_gs_current_lambda=lambda
    dg_gs_hamiltonian_potential_epoch=dg_gs_potential_epoch
    dg_gs_operator_ok=.true.
    dg_gs_operator_message=''
    call solve_nodal_ground_state_cg_mpi(state,dg_dc_apply_complete_h_for_main, &
      ncg, &
      merge(dg_gs_controls%final_orbital_tolerance,dg_gs_controls%intermediate_orbital_tolerance,lambda==1d0), &
      communicator,dg_gs_eigenvalues,orbital_residual,eigensolver_iterations,dg_dc_rotate_subspace_for_main,.false.)
    if(.not.dg_gs_operator_ok) then
      step_ok=.false.; step_message=dg_gs_operator_message; return
    end if
    call dg_dc_assign_ground_state_occupations(state,step_ok,step_message)
    if(.not.step_ok) return
    diagnostics%orbital_residual=orbital_residual
    diagnostics%subspace_residual=orbital_residual
    diagnostics%projector_overlap=0d0
    diagnostics%hermiticity_defect=dg_gs_face_hermiticity
    diagnostics%face_balance_defect=dg_gs_face_balance
    diagnostics%interface_scale=lambda
    diagnostics%eigensolver_iterations=eigensolver_iterations
    diagnostics%hamiltonian_potential_epoch=dg_gs_hamiltonian_potential_epoch
    diagnostics%eigensolver_converged=orbital_residual<=merge(dg_gs_controls%final_orbital_tolerance, &
      dg_gs_controls%intermediate_orbital_tolerance,lambda==1d0)
    step_ok=.true.
    step_message=''
  end subroutine dg_dc_salmon_solve_step

  subroutine dg_dc_salmon_update_step(state,density_arg,density_mix,unmixed,communicator,diagnostics, &
      step_ok,step_message)
    type(s_dg_nodal_common_state), intent(inout) :: state
    real(8), intent(inout) :: density_arg(:,:,:,:)
    real(8), intent(in) :: density_mix
    logical, intent(in) :: unmixed
    integer, intent(in) :: communicator
    type(s_dg_dc_gs_diagnostics), intent(inout) :: diagnostics
    logical, intent(out) :: step_ok
    character(*), intent(out) :: step_message
    real(8) :: density_residual,electron_number,orthogonality_defect

    call dg_dc_update_density_potential(state,density_arg,density_mix,unmixed,density_residual, &
      electron_number,step_ok,step_message)
    if(.not.step_ok) return
    call dg_dc_orthogonality_defect(state,communicator,orthogonality_defect,step_ok)
    if(.not.step_ok) then
      step_message='DG DC GS: orthogonality diagnostic reduction failed'
      return
    end if
    diagnostics%density_residual=density_residual
    diagnostics%orthogonality_defect=orthogonality_defect
    diagnostics%electron_number=electron_number
    diagnostics%expected_electron_number=dc%elec_num_tot
    diagnostics%updated_potential_epoch=dg_gs_potential_epoch
    diagnostics%finite=all(ieee_is_finite([diagnostics%orbital_residual,density_residual,electron_number, &
      orthogonality_defect,dg_gs_face_hermiticity,dg_gs_face_balance]))
    step_ok=diagnostics%finite
    if(step_ok) then
      step_message=''
    else
      step_message='DG DC GS: nonfinite SALMON production diagnostics'
    end if
  end subroutine dg_dc_salmon_update_step

  subroutine dg_dc_apply_complete_h_for_main(state,hpsi_complete)
    type(s_dg_nodal_common_state), intent(inout) :: state
    complex(8), intent(out) :: hpsi_complete(:,:,:,:,:)
    complex(8), allocatable :: hpsi_volume_nonlocal(:,:,:,:,:),hpsi_sipg(:,:,:,:,:)
    logical :: action_ok
    character(256) :: action_message
    if(.not.dg_gs_operator_ok) then
      hpsi_complete=(0d0,0d0)
      return
    end if
    allocate(hpsi_volume_nonlocal,mold=state%psi_core)
    allocate(hpsi_sipg,mold=state%psi_core)
    call dg_dc_apply_volume_nonlocal_for_main(state,hpsi_volume_nonlocal,action_ok,action_message)
    action_ok=action_ok .and. dg_gs_operator_ok
    call comm_logical_and(action_ok,dg_gs_operator_ok,dc%icomm_tot)
    if(.not.dg_gs_operator_ok) then
      hpsi_complete=(0d0,0d0); dg_gs_operator_message='DG DC volume/nonlocal Hamiltonian action failed'; return
    end if
    call apply_dg_dc_sipg_operator_mpi(state,dg_gs_fragment_origins,dg_gs_fragment_sizes, &
      dc%lg_tot%num,dc%system_tot%hgs,81d0,dc%icomm_tot,hpsi_sipg,dg_gs_face_hermiticity, &
      dg_gs_face_balance,action_ok,action_message)
    action_ok=action_ok .and. dg_gs_operator_ok
    call comm_logical_and(action_ok,dg_gs_operator_ok,dc%icomm_tot)
    if(.not.dg_gs_operator_ok) then
      hpsi_complete=(0d0,0d0); dg_gs_operator_message='DG DC SIPG Hamiltonian action failed'; return
    end if
    call compose_dg_dc_hamiltonian(hpsi_volume_nonlocal,hpsi_sipg,dg_gs_current_lambda, &
      hpsi_complete,action_ok,action_message)
    action_ok=action_ok .and. dg_gs_operator_ok
    call comm_logical_and(action_ok,dg_gs_operator_ok,dc%icomm_tot)
    if(.not.dg_gs_operator_ok) then
      hpsi_complete=(0d0,0d0); dg_gs_operator_message='DG DC complete Hamiltonian composition failed'
    end if
  end subroutine dg_dc_apply_complete_h_for_main

  subroutine dg_dc_apply_volume_nonlocal_for_main(state,hpsi_volume_nonlocal,ok,message)
    type(s_dg_nodal_common_state), intent(in) :: state
    complex(8), intent(out) :: hpsi_volume_nonlocal(:,:,:,:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: first_state,batch_count,slot,global_state,is,ix,iy,iz
    real(8) :: scale
    real(8), allocatable :: vlocal_core(:,:,:,:)
    complex(8), allocatable :: zero_extended_local(:,:,:,:,:),interior_volume(:,:,:,:,:)
    hpsi_volume_nonlocal=(0d0,0d0)
    ok=system%if_real_orbital .and. all(shape(hpsi_volume_nonlocal)==shape(state%psi_core))
    if(.not.ok) then
      message='DG DC GS: unsupported volume/nonlocal orbital layout'
      return
    end if
    scale=sqrt(dc%system_tot%hvol)
    do first_state=1,state%nstate,system%no
      batch_count=min(system%no,state%nstate-first_state+1)
      spsi%rwf=0d0
      do slot=1,batch_count
        global_state=first_state+slot-1
        do is=1,system%nspin
          do iz=1,state%core_size(3); do iy=1,state%core_size(2); do ix=1,state%core_size(1)
            spsi%rwf(ix,iy,iz,is,slot,info%ik_s,info%im_s)= &
              real(state%psi_core(ix,iy,iz,global_state,is),8)/scale
          end do; end do; end do
        end do
      end do
      call hpsi(spsi,shpsi,info,mg,v_local,system,stencil,srg,ppg)
      do slot=1,batch_count
        global_state=first_state+slot-1
        do is=1,system%nspin
          do iz=1,state%core_size(3); do iy=1,state%core_size(2); do ix=1,state%core_size(1)
            hpsi_volume_nonlocal(ix,iy,iz,global_state,is)= &
              cmplx(scale*shpsi%rwf(ix,iy,iz,is,slot,info%ik_s,info%im_s),0d0,8)
          end do; end do; end do
        end do
      end do
    end do
    allocate(vlocal_core(state%core_size(1),state%core_size(2),state%core_size(3),state%nspin))
    do is=1,state%nspin
      vlocal_core(:,:,:,is)=v_local(is)%f(1:state%core_size(1),1:state%core_size(2),1:state%core_size(3))
    end do
    allocate(zero_extended_local,mold=state%psi_core)
    allocate(interior_volume,mold=state%psi_core)
    call build_dg_dc_interior_volume_action(state%psi_core,vlocal_core,stencil%coef_lap0, &
      stencil%coef_lap,zero_extended_local,interior_volume,ok,message)
    if(.not.ok) return
    hpsi_volume_nonlocal=hpsi_volume_nonlocal-zero_extended_local+interior_volume
    ok=all(ieee_is_finite(real(hpsi_volume_nonlocal,8))) .and. &
      all(ieee_is_finite(aimag(hpsi_volume_nonlocal)))
    if(ok) then
      message=''
    else
      message='DG DC GS: nonfinite volume/nonlocal action'
    end if
  end subroutine dg_dc_apply_volume_nonlocal_for_main

  subroutine dg_dc_rotate_subspace_for_main(state,hpsi_complete,eigenvalues,communicator)
    type(s_dg_nodal_common_state), intent(inout) :: state
    complex(8), intent(inout) :: hpsi_complete(:,:,:,:,:)
    real(8), intent(out) :: eigenvalues(state%nstate,state%nspin)
    integer, intent(in) :: communicator
    call rayleigh_ritz_nodal_subspace_mpi(state,hpsi_complete,eigenvalues,communicator)
  end subroutine dg_dc_rotate_subspace_for_main

  subroutine dg_dc_assign_ground_state_occupations(state,ok,message)
    type(s_dg_nodal_common_state), intent(inout) :: state
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: nfilled
    real(8) :: remaining
    state%occupations=0d0
    nfilled=min(state%nstate,int(dc%elec_num_tot/2d0))
    if(nfilled>0) state%occupations(1:nfilled,1)=2d0
    remaining=dc%elec_num_tot-2d0*dble(nfilled)
    if(remaining>10d0*epsilon(1d0) .and. nfilled<state%nstate) state%occupations(nfilled+1,1)=remaining
    ok=remaining<2d0+10d0*epsilon(1d0) .and. dc%elec_num_tot<=2d0*dble(state%nstate)
    if(ok) then
      message=''
    else
      message='DG DC GS: insufficient global candidate states for electron count'
    end if
  end subroutine dg_dc_assign_ground_state_occupations

  subroutine dg_dc_update_density_potential(state,density_arg,density_mix,unmixed,density_residual, &
      electron_number,ok,message)
    type(s_dg_nodal_common_state), intent(in) :: state
    real(8), intent(inout) :: density_arg(:,:,:,:)
    real(8), intent(in) :: density_mix
    logical, intent(in) :: unmixed
    real(8), intent(out) :: density_residual,electron_number
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: is,ix,iy,iz
    real(8) :: local_norm,global_norm,local_charge
    call reconstruct_dg_dc_core_density(state,state%nstate,dg_gs_core_density,ok,message)
    if(.not.ok) return
    dg_gs_core_density=dg_gs_core_density/dc%system_tot%hvol
    do is=1,system%nspin
      rho_s(is)%f=0d0
      do iz=1,state%core_size(3); do iy=1,state%core_size(2); do ix=1,state%core_size(1)
        rho_s(is)%f(ix,iy,iz)=dg_gs_core_density(ix,iy,iz,is)
      end do; end do; end do
    end do
    call calc_rho_total_dcdft(system%nspin,lg,mg,info,rho_s,dc)
    do is=1,system%nspin
      dg_gs_raw_density(:,:,:,is)=dc%rho_tot_s(is)%f
    end do
    local_norm=sum((dg_gs_raw_density-density_arg)**2)/real(dc%isize_tot,8)
    call comm_summation(local_norm,global_norm,dc%icomm_tot)
    density_residual=sqrt(global_norm*dc%system_tot%hvol)
    local_charge=sum(dg_gs_raw_density)*dc%system_tot%hvol/real(dc%isize_tot,8)
    call comm_summation(local_charge,electron_number,dc%icomm_tot)
    call mix_dg_dc_density_transaction(density_arg,dg_gs_raw_density,density_mix,unmixed, &
      dg_gs_mixed_density,ok,message)
    if(.not.ok) return
    density_arg=dg_gs_mixed_density
    call dg_dc_update_potential_from_density(density_arg,ok,message)
  end subroutine dg_dc_update_density_potential

  subroutine dg_dc_update_potential_from_density(density_arg,ok,message)
    real(8), intent(in) :: density_arg(:,:,:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: is
    do is=1,system%nspin
      dc%rho_tot_s(is)%f=density_arg(:,:,:,is)
    end do
    dc%rho_tot%f=0d0
    do is=1,system%nspin
      dc%rho_tot%f=dc%rho_tot%f+dc%rho_tot_s(is)%f
    end do
    call hartree(dc%lg_tot,dc%mg_tot,dc%info_tot,dc%system_tot,dc%fg_tot,dc%poisson_tot, &
      dc%srg_scalar_tot,stencil,dc%rho_tot,dc%Vh_tot)
    call exchange_correlation(dc%system_tot,xc_func,dc%mg_tot,srg_scalar,srg,dc%rho_tot_s, &
      pp,ppn,dc%info_tot,spsi,stencil,dc%Vxc_tot,energy%E_xc)
    call update_vlocal(dc%mg_tot,dc%system_tot%nspin,dc%Vh_tot,dc%Vpsl_tot,dc%Vxc_tot,dc%vloc_tot)
    call calc_vlocal_fragment_dcdft(system%nspin,mg,v_local,dc)
    dg_gs_potential_epoch=dg_gs_potential_epoch+1_8
    ok=all(ieee_is_finite(dc%rho_tot%f)) .and. all(ieee_is_finite(dc%Vh_tot%f))
    do is=1,system%nspin
      ok=ok .and. all(ieee_is_finite(dc%Vxc_tot(is)%f)) .and. all(ieee_is_finite(dc%vloc_tot(is)%f))
    end do
    if(ok) then
      message=''
    else
      message='DG DC GS: nonfinite DC Hartree/XC/vlocal update'
    end if
  end subroutine dg_dc_update_potential_from_density

  subroutine dg_dc_orthogonality_defect(state,communicator,defect,ok)
    type(s_dg_nodal_common_state), intent(in) :: state
    integer, intent(in) :: communicator
    real(8), intent(out) :: defect
    logical, intent(out) :: ok
    complex(8), allocatable :: local_overlap(:,:),global_overlap(:,:)
    integer :: i,j,ierr
    allocate(local_overlap(state%nstate,state%nstate),global_overlap(state%nstate,state%nstate))
    local_overlap=(0d0,0d0)
    do j=1,state%nstate
      do i=1,state%nstate
        local_overlap(i,j)=sum(conjg(state%psi_core(:,:,:,i,1))*state%psi_core(:,:,:,j,1))
      end do
    end do
    call comm_summation(local_overlap,global_overlap,size(local_overlap),communicator)
    do i=1,state%nstate
      global_overlap(i,i)=global_overlap(i,i)-(1d0,0d0)
    end do
    defect=maxval(abs(global_overlap))
    ok=ieee_is_finite(defect)
  end subroutine dg_dc_orthogonality_defect

  integer function canonical_to_dc_index(index,core_count,buffer_count)
    integer, intent(in) :: index,core_count,buffer_count
    if(index<=buffer_count) then
      canonical_to_dc_index=core_count+buffer_count+index
    else if(index<=buffer_count+core_count) then
      canonical_to_dc_index=index-buffer_count
    else
      canonical_to_dc_index=core_count+index-(buffer_count+core_count)
    end if
  end function canonical_to_dc_index

  integer(8) function dg_dc_geometry_fingerprint()
    integer(8) :: local_hash
    integer :: ii,jj
    local_hash=not(0_8)
    call hash_integer(local_hash,dc%lg_tot%num(1))
    call hash_integer(local_hash,dc%lg_tot%num(2))
    call hash_integer(local_hash,dc%lg_tot%num(3))
    call hash_real(local_hash,dc%system_tot%Hvol)
    do jj=1,3
    do ii=1,3
      call hash_real(local_hash,dc%system_tot%primitive_a(ii,jj))
      call hash_real(local_hash,dc%system_tot%primitive_b(ii,jj))
    end do
    end do
    do jj=1,dc%n_frag
    do ii=1,3
      call hash_integer(local_hash,dc%nxyz_domain_frag(ii,jj))
      call hash_integer(local_hash,dc%ixyz_frag(ii,jj))
    end do
    end do
    do jj=1,dc%system_tot%nion
      call hash_integer(local_hash,dc%system_tot%kion(jj))
      do ii=1,3
        call hash_real(local_hash,dc%system_tot%Rion(ii,jj))
      end do
    end do
    if(local_hash==0_8) local_hash=1_8
    dg_dc_geometry_fingerprint=local_hash
  end function dg_dc_geometry_fingerprint

  integer(8) function dg_dc_operator_fingerprint(include_potential)
    logical,intent(in),optional::include_potential
    integer(8) :: local_hash,point_hash,global_potential_hash,point_owner
    integer :: ix,iy,iz,is,ii,jj,kk,mpi_error
    logical::hash_potential
    hash_potential=.true.;if(present(include_potential))hash_potential=include_potential
    local_hash=0_8;global_potential_hash=0_8
    if(hash_potential)then
    do is=1,dc%system_tot%nspin
    do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)
    do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
    do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
      point_owner=int(ix-dc%lg_tot%is(1)+1,8)+int(dc%lg_tot%num(1),8)* &
        (int(iy-dc%lg_tot%is(2),8)+int(dc%lg_tot%num(2),8)* &
        (int(iz-dc%lg_tot%is(3),8)+int(dc%lg_tot%num(3),8)*int(is-1,8)))
      point_hash=not(0_8)
      call hash_integer8(point_hash,point_owner)
      call hash_real(point_hash,dc%vloc_tot(is)%f(ix,iy,iz))
      local_hash=ieor(local_hash,point_hash)
    end do
    end do
    end do
    end do
#ifdef USE_MPI
    call MPI_Allreduce(local_hash,global_potential_hash,1,MPI_INTEGER8,MPI_BXOR,dc%icomm_tot,mpi_error)
    if(mpi_error/=MPI_SUCCESS) stop 'DG DC operator fingerprint reduction failed'
#else
    global_potential_hash=local_hash
#endif
    endif
    local_hash=not(0_8)
    call hash_integer8(local_hash,global_potential_hash)
    call hash_integer(local_hash,dc%system_tot%nspin)
    call hash_real(local_hash,dg_dc_gs_sipg_penalty_factor)
    call hash_real(local_hash,dg_dc_gs_target_lambda)
    call hash_real(local_hash,stencil%coef_lap0)
    do jj=1,3
    do ii=1,4
      call hash_real(local_hash,stencil%coef_lap(ii,jj))
      call hash_real(local_hash,stencil%coef_nab(ii,jj))
    end do
    end do
    call hash_integer(local_hash,pp%lmax)
    call hash_integer(local_hash,pp%nrmax)
    if(allocated(pp%zps)) then
      do ii=1,size(pp%zps)
        call hash_integer(local_hash,pp%zps(ii))
        call hash_real(local_hash,pp%rloc(ii))
        call hash_real(local_hash,pp%rps(ii))
      end do
    end if
    if(allocated(pp%rad)) then
      do jj=1,size(pp%rad,2)
      do ii=1,size(pp%rad,1)
        call hash_real(local_hash,pp%rad(ii,jj))
      end do
      end do
    end if
    if(allocated(pp%radnl)) then
      do jj=1,size(pp%radnl,2)
      do ii=1,size(pp%radnl,1)
        call hash_real(local_hash,pp%radnl(ii,jj))
      end do
      end do
    end if
    if(allocated(pp%vloctbl)) then
      do jj=1,size(pp%vloctbl,2)
      do ii=1,size(pp%vloctbl,1)
        call hash_real(local_hash,pp%vloctbl(ii,jj))
      end do
      end do
    end if
    if(allocated(pp%udvtbl)) then
      do kk=1,size(pp%udvtbl,3)
      do jj=1,size(pp%udvtbl,2)
      do ii=1,size(pp%udvtbl,1)
        call hash_real(local_hash,pp%udvtbl(ii,jj,kk))
      end do
      end do
      end do
    end if
    call hash_character(local_hash,trim(xc))
    if(local_hash==0_8) local_hash=1_8
    dg_dc_operator_fingerprint=local_hash
  end function dg_dc_operator_fingerprint

  subroutine hash_integer(hash,value)
    integer(8), intent(inout) :: hash
    integer, intent(in) :: value
    integer(8) :: bits
    integer :: ibyte
    bits=int(value,8)
    do ibyte=0,7
      call hash_byte(hash,int(ibits(bits,8*ibyte,8)))
    end do
    if(hash==0_8) hash=1_8
  end subroutine hash_integer

  subroutine hash_integer8(hash,value)
    integer(8), intent(inout) :: hash
    integer(8), intent(in) :: value
    integer :: ibyte
    do ibyte=0,7
      call hash_byte(hash,int(ibits(value,8*ibyte,8)))
    end do
    if(hash==0_8) hash=1_8
  end subroutine hash_integer8

  subroutine hash_real(hash,value)
    integer(8), intent(inout) :: hash
    real(8), intent(in) :: value
    integer(8) :: bits
    integer :: ibyte
    bits=transfer(value,bits)
    do ibyte=0,7
      call hash_byte(hash,int(ibits(bits,8*ibyte,8)))
    end do
    if(hash==0_8) hash=1_8
  end subroutine hash_real

  subroutine hash_character(hash,value)
    integer(8), intent(inout) :: hash
    character(*), intent(in) :: value
    integer :: ii
    do ii=1,len(value)
      call hash_byte(hash,iachar(value(ii:ii)))
    end do
  end subroutine hash_character

  subroutine hash_byte(hash,value)
    integer(8), intent(inout) :: hash
    integer, intent(in) :: value
    integer(8), parameter :: polynomial=int(z'C96C5795D7870F42',8)
    integer :: ibit
    hash=ieor(hash,int(iand(value,255),8))
    do ibit=1,8
      if(btest(hash,0)) then
        hash=ieor(shiftr(hash,1),polynomial)
      else
        hash=shiftr(hash,1)
      end if
    end do
  end subroutine hash_byte

end subroutine main_dft
