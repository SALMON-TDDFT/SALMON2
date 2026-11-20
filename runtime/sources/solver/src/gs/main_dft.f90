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

module salmonpy_gs_context
  use filesystem,only:resolve_path
use iso_c_binding
use timer,only:s_timer_context,timer_context_swap,timer_context_get,timer_initialize, &
 LOG_TOTAL,LOG_INIT_GS,LOG_GS_ITERATION,LOG_WRITE_GS_RESULTS
use structures
use salmonpy_model_inputs,only:s_model_inputs
use salmonpy_gs_seed,only:s_gs_seed
use sym_sub,only:s_hse_symmetry_workspace
use sym_rho_sub,only:s_hse_density_workspace
use sendrecv_grid,only:dealloc_cache
#ifdef USE_MPI
use communication,only:comm_begin_children,comm_end_children,comm_release_children, &
 comm_is_attached,comm_get_globalinfo,comm_get_max,comm_bcast,comm_child_counts
#endif
use salmon_xc,only:finalize_xc
#ifdef USE_HSE
use exx_rt_state,only:s_exx_rt_state,rt_state_clear,rt_values_copy
use exx_spatial,only:spatial_exx_state
use exx_distributed_gauge,only:s_exx_gauge
use exx_sparse_orbitals,only:s_sparse_orbitals
use exx_ace,only:s_exx_ace
use exx_symmetry,only:s_exx_symmetry_map
use exx_factor_history,only:s_factor_history
use exx_native,only:exx_unbind_system
#endif
implicit none
private
public::s_gs_context,gs_run,gs_close,gs_probe_write,gs_input_config,gs_export_seed
#ifdef USE_MPI
public::gs_seed_from_handle
#endif
! Internal retained GS storage. Actual contexts must have TARGET and a stable address.
! Never copy a live context; global inputs remain shared.
! Numerical-only snapshots; static communicators/Libxc/FFT handles stay live.
type s_scf_payload
#ifdef USE_HSE
 type(s_exx_rt_state) :: hse
#endif
 type(s_model_inputs) :: continuation_inputs
 type(s_dft_system) :: system
 type(s_orbital) :: spsi,shpsi,sttpsi
 type(s_cg) :: cg
 type(s_mixing) :: mixing
 type(s_scalar) :: rho,rho_jm,Vh,Vpsl
 type(s_scalar),allocatable :: V_local(:),rho_s(:),Vxc(:)
 type(s_poisson) :: poisson
 type(s_reciprocal_grid) :: fg
 type(s_dft_energy) :: energy
 integer :: iterations=0
 logical :: converged=.false.
 real(8) :: residual=0d0,threshold=0d0,cval=0d0
end type
type s_scf_snapshot_slot
 integer(c_int64_t) :: token=0
 type(s_scf_payload),allocatable :: value
end type
integer(c_int64_t),save :: scf_snapshot_serial=0
type s_gs_context
 private
 type(s_timer_context)::timers
 logical::active=.false.,imported=.false.,scf_capable=.false.
 integer :: initializations=0,executed_iterations=0
 type(s_scf_snapshot_slot),allocatable :: scf_snapshots(:)
 type(s_scf_payload),allocatable :: scf_pending
 integer::mpi_child_owner=0
 type(s_hse_symmetry_workspace)::symmetry
 type(s_hse_density_workspace)::density_symmetry
 type(s_model_inputs)::inputs,continuation_inputs
 integer::iterations=0
 logical::converged=.false.
 real(8)::residual=0d0,threshold=0d0
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
type(s_unfold) :: unfold

#ifdef USE_HSE
 type(s_exx_rt_state)::hse
 logical::hse_bound=.false.
#endif
end type
#ifdef USE_MPI
type gs_handle_slot
 type(s_gs_context),pointer::context=>null()
 integer(c_int64_t)::token=0
end type
type(gs_handle_slot),allocatable,save::gs_handles(:)
integer(c_int64_t),save::gs_handle_serial=0
integer(c_int64_t),parameter::gs_handle_prefix=int(z'4A00000000000000',c_int64_t)
logical,save::gs_host_initialized=.false.
#endif
contains
subroutine gs_input_config(value)
 use salmon_global,only:nscf,threshold,exx_pre_scf_threshold,exx_pre_scf_steps,convergence,nelec,nelec_spin
 use mixing_sub,only:mixing_input_config
 type(s_gs_run_config),intent(out)::value
 call mixing_input_config(value%mixing)
 value%normalization%nelec=nelec;value%normalization%nelec_spin=nelec_spin
 value%convergence_metric=convergence
 value%max_iterations=nscf;value%threshold=threshold
 value%pre_scf_threshold=exx_pre_scf_threshold;value%pre_scf_steps=exx_pre_scf_steps
end subroutine
subroutine gs_run(ctx,status,config,import_only,disk_bootstrap)
use,intrinsic::ieee_arithmetic,only:ieee_is_finite
use exx_functional,only:is_hybrid
use math_constants, only: pi, zi
use structures
use inputoutput
use parallelization, only: nproc_id_global,nproc_group_global,adjust_elapse_time,nproc_size_global
use communication, only: comm_is_root, comm_summation, comm_bcast, comm_sync_all,comm_get_max
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
use lcfo
use atomic_orbital, only: verify_atomic_orbital_consistency !!DEBUG
#ifdef USE_HSE
use exx_native,only:exx_pre_scf_initialize,exx_pre_scf,exx_set_phase,exx_bind_gs_system,exx_unbind_system
use exx_rt_state,only:s_exx_rt_state,rt_state_clear
#endif
implicit none
type(s_gs_context),target,intent(inout)::ctx
integer,intent(out)::status
type(s_gs_run_config),optional,intent(in)::config
type(s_gs_run_config)::selected
logical,optional,intent(in)::import_only,disk_bootstrap
logical::reading_only
#ifdef USE_HSE
 integer :: gs_hse_status
#endif
integer :: phase_status,gs_status_unit,force_status_unit
integer :: ix,iy,iz
integer :: Miter,iatom,jj,nspin
logical :: pre_scf_root
real(8) :: sum1
character(100) :: comment_line


logical :: rion_update
logical :: flag_opt_conv
integer :: Miopt, iopt,nopt_max,i
integer :: iter_band_kpt, iter_band_kpt_end, iter_band_kpt_stride
logical :: is_checkpoint_iter, is_shutdown_time
integer :: ilevel_print

status=1
if(ctx%active)return
reading_only=.false.
if(present(import_only))reading_only=import_only
if(reading_only)then
 if(yn_restart/='y'.or.read_gs_restart_data/='wfn')return
 if(yn_dc/='n'.or.yn_opt/='n'.or.yn_md/='n')return
endif
if(present(config))then
 selected=config
else
 call gs_input_config(selected)
endif
! Legacy -1 skips the >1 no-mixing loop, just as zero does.
if(selected%mixing%initial_no_mix_steps==-1)selected%mixing%initial_no_mix_steps=0
! Imported orbitals do not run convergence iterations or PBE warmup.
if(reading_only)then
 if(selected%mixing%initial_no_mix_steps==-1)selected%mixing%initial_no_mix_steps=0
 selected%pre_scf_threshold=0d0;selected%pre_scf_steps=1
endif
if(.not.valid_mixing_input(selected%mixing))return
if(selected%mixing%method=='simple_dm'.and.yn_spinorbit/='y')return
if(.not.reading_only)then
select case(selected%convergence_metric)
case('rho_dne','norm_rho','norm_rho_dng','norm_pot','norm_pot_dng')
case default
 return
end select
if(selected%convergence_metric=='rho_dne'.and.yn_dc/='y')then
 if(gs_density_normalization(selected%normalization,merge(1,2,spin=='unpolarized'))<=0d0)return
endif
if(selected%max_iterations<1.or.selected%max_iterations>=huge(0))return
if(selected%pre_scf_steps<1)return
if(.not.ieee_is_finite(selected%threshold))return
if(.not.ieee_is_finite(selected%pre_scf_threshold))return
if(selected%threshold<0d0.or.selected%pre_scf_threshold<0d0)return
! Explicit values must satisfy the same semantic requirements as namelist inputs.
if(selected%pre_scf_threshold>0d0)then
 if(theory/='dft'.or.yn_restart/='n'.or.yn_md/='n'.or.yn_opt/='n'.or..not.is_hybrid(xc))return
 if(yn_dc=='y'.and.temperature<=0d0)return
 if(checkpoint_interval>0.or.time_shutdown>0d0)return
 if(selected%mixing%method=='simple_potential')return
 if(selected%convergence_metric/='rho_dne'.and.selected%convergence_metric/='norm_rho'.and. &
    selected%convergence_metric/='norm_rho_dng')return
endif
endif
if(theory=='dft_band'.and.iperiodic/=3) return
#ifdef USE_MPI
if(comm_is_attached())then
 call comm_begin_children(ctx%mpi_child_owner,status)
 if(status/=0)return
endif
#endif
ctx%active=.true.
call ctx%inputs%capture()
ctx%system%gs_input=selected;ctx%system%gs_input_owned=.true.
ctx%mixing%input=selected%mixing;ctx%mixing%input_owned=.true.

if(yn_dc=='y') then
  call init_dcdft(ctx%dc,ctx%pp,ctx%mixing,ctx%ewald)
  ilevel_print = 0
else
  ilevel_print = 3
end if

!check condition for using jellium model
if(yn_jm=='y') call check_condition_jm

! Set the temporary functional before initialization evaluates any potential.
#ifdef USE_HSE
call exx_pre_scf_initialize(ctx%system,selected%pre_scf_threshold,selected%pre_scf_steps,phase_status)
if(phase_status/=0)error stop 'Invalid pre-SCF configuration'
if(yn_dc=='y')call exx_set_phase(ctx%dc%system_tot,exx_pre_scf(ctx%system))
pre_scf_root=comm_is_root(nproc_id_global)
if(yn_dc=='y')pre_scf_root=comm_is_root(ctx%dc%id_tot)
if(exx_pre_scf(ctx%system).and.pre_scf_root) &
  write(*,'(a,es14.6,a,i0)') 'EXX_PRE_SCF PBE stage: threshold=',ctx%system%pre_scf%threshold, &
                          ' consecutive steps=',ctx%system%pre_scf%required_steps
#endif
call init_xc(ctx%xc_func, spin, cval, xcname=xc, xname=xname, cname=cname)

call timer_begin(LOG_TOTAL)
call timer_begin(LOG_INIT_GS)


! please move folloings into initialization_dft
call init_dft(nproc_group_global,ctx%info,ctx%lg,ctx%mg,ctx%system,ctx%stencil,ctx%fg,ctx%poisson, &
  ctx%srg,ctx%srg_scalar,ctx%ofl,ctx%unfold)
#ifdef USE_HSE
call exx_bind_gs_system(ctx%system,ctx%info,ctx%hse,gs_hse_status,allow_symmetry=ctx%mpi_child_owner>0)
if(gs_hse_status==1)error stop 'GS HSE owner binding conflict'
ctx%hse_bound=gs_hse_status==0
if(ctx%mixing%scf_progress%stop_after>0.and.trim(xc)=='hse06'.and..not.ctx%hse_bound)then
 call timer_end(LOG_INIT_GS)
 call timer_end(LOG_TOTAL)
#ifdef USE_MPI
 if(ctx%mpi_child_owner>0)call comm_end_children(ctx%mpi_child_owner)
#endif
 status=15
 return
endif
if(ctx%hse_bound.and.comm_is_root(nproc_id_global))write(*,'(a)')'EXX_GS_OWNER explicit'
#endif
allocate( ctx%rho_s(ctx%system%nspin),ctx%v_local(ctx%system%nspin),ctx%vxc(ctx%system%nspin) )

ctx%initializations=ctx%initializations+1
call initialization1_dft( ctx%system, ctx%energy, ctx%stencil, ctx%fg, ctx%poisson,  &
                          ctx%lg, ctx%mg,   &
                          ctx%info,  &
                          ctx%srg, ctx%srg_scalar,  &
                          ctx%rho, ctx%rho_jm, ctx%rho_s, ctx%vh, ctx%v_local, ctx%vpsl, ctx%vxc,  &
                          ctx%spsi, ctx%shpsi, ctx%sttpsi,  &
                          ctx%pp, ctx%ppg, ctx%ppn,  &
                          ctx%ofl )

! Disk restart recreates static resources only: no guessed orbitals/density,
! mixing initialization, eigenvalue/energy evaluation, or SCF replay.
if(present(disk_bootstrap))then
 if(disk_bootstrap)then
  ctx%ewald%yn_bookkeep='y'
  call init_nion_div(ctx%system,ctx%lg,ctx%mg,ctx%info)
  call init_ewald(ctx%system,ctx%info,ctx%ewald)
  call ctx%continuation_inputs%capture()
  ctx%mixing%scf_progress%paused=.true.
  call timer_end(LOG_INIT_GS)
  call timer_end(LOG_TOTAL)
#ifdef USE_MPI
  if(ctx%mpi_child_owner>0)call comm_end_children(ctx%mpi_child_owner)
#endif
  status=0
  return
 endif
endif

call initialization2_dft( Miter, nspin, rion_update,  &
                          ctx%system, ctx%energy, ctx%ewald, ctx%stencil, ctx%fg, ctx%poisson,&
                          ctx%lg, ctx%mg, ctx%info,   &
                          ctx%srg, ctx%srg_scalar,  &
                          ctx%rho, ctx%rho_jm, ctx%rho_s, ctx%vh,ctx%v_local, ctx%vpsl, ctx%vxc,  &
                          ctx%spsi, ctx%shpsi, ctx%sttpsi,  &
                          ctx%pp, ctx%ppg, ctx%ppn,   &
                          ctx%xc_func, ctx%mixing )


! Import through the original initialization without entering SCF or GS output.
if(reading_only)then
 call timer_end(LOG_INIT_GS)
 call timer_end(LOG_TOTAL)
#ifdef USE_MPI
 call comm_end_children(ctx%mpi_child_owner)
#endif
 ctx%iterations=Miter;ctx%imported=.true.
 status=0
 return
endif
Miopt = 0
nopt_max = 1
if(yn_opt=='y') call initialization_opt(Miopt,ctx%opt,ctx%system,flag_opt_conv,nopt_max,ctx%ofl)

call timer_end(LOG_INIT_GS)

!---------------------------------------- Opt Iteration


#ifdef __FUJITSU
call fipp_start ! performance profiling
#endif

Structure_Optimization_Iteration : do iopt= Miopt+1, nopt_max

if(iopt>=2)then
  call timer_begin(LOG_INIT_GS)
  Miter = 0        ! Miter: Iteration counter set to zero
  rion_update = .true.
  call dealloc_init_ps(ctx%ppg)
  call init_ps(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%fg,ctx%poisson,ctx%pp,ctx%ppg,ctx%vpsl)
  call calc_nlcc(ctx%pp, ctx%system, ctx%mg, ctx%ppn)
  if(ctx%mixing%input%auto_mixing=='y') call reset_mixing_rate(ctx%mixing)
  call timer_end(LOG_INIT_GS)
end if

!---------------------------------------- Band Iteration

if(theory=='dft_band')then
   call init_band_dft(ctx%system,ctx%band) ! --> system%wtk=0.0
   iter_band_kpt_end    = ctx%band%num_band_kpt
   iter_band_kpt_stride = ctx%system%nk
else
   iter_band_kpt_end    = 1
   iter_band_kpt_stride = 1
end if

call comm_sync_all
call timer_enable_sub
Band_Iteration : do iter_band_kpt= 1, iter_band_kpt_end, iter_band_kpt_stride

if(theory=='dft_band')then
   call calc_band_write(iter_band_kpt,ctx%system,ctx%band,ctx%info)
end if


call timer_begin(LOG_INIT_GS_ITERATION)

call timer_end(LOG_INIT_GS_ITERATION)


call timer_begin(LOG_GS_ITERATION)
!------------------------------------ SCF Iteration
!Iteration loop for SCF (DFT_Iteration)
call scf_iteration_dft( Miter,rion_update,sum1,  &
                        ctx%system,ctx%energy,ctx%ewald,  &
                        ctx%lg,ctx%mg,  &
                        ctx%info,  &
                        ctx%poisson,ctx%fg,  &
                        ctx%cg,ctx%mixing,  &
                        ctx%stencil,  &
                        ctx%srg,ctx%srg_scalar,   &
                        ctx%spsi,ctx%shpsi,ctx%sttpsi,  &
                        ctx%rho,ctx%rho_jm,ctx%rho_s,  &
                        ctx%v_local,ctx%vh,ctx%vxc,ctx%vpsl,ctx%xc_func,  &
                        ctx%pp,ctx%ppg,ctx%ppn,  &
                        ctx%band, ilevel_print, ctx%dc)
ctx%converged=ctx%mixing%scf_converged
call ctx%continuation_inputs%capture()
ctx%executed_iterations=ctx%executed_iterations+Miter
if(ctx%mixing%scf_progress%paused)then
 ctx%iterations=Miter;ctx%residual=sum1;ctx%threshold=selected%threshold
 call timer_end(LOG_GS_ITERATION)
 call timer_disable_sub
 call timer_end(LOG_TOTAL)
#ifdef USE_MPI
 if(ctx%mpi_child_owner>0)call comm_end_children(ctx%mpi_child_owner)
#endif
 status=0
 return
endif
! Scalar bridge reports mirror native values without changing SCF arithmetic.
if(theory=='dft'.and.yn_opt=='n'.and.yn_md=='n'.and.yn_dc=='n')then
 if(comm_is_root(nproc_id_global))then
  open(newunit=gs_status_unit,file=resolve_path('gs_status.dat'),status='replace')
  write(gs_status_unit,'(a)')'SALMONPY_GS_V1'
  write(gs_status_unit,*)merge(1,0,ctx%converged),Miter
  write(gs_status_unit,'(3es26.17e3)')sum1,selected%threshold,ctx%energy%E_tot
  write(gs_status_unit,'(a)')trim(selected%convergence_metric)
  close(gs_status_unit)
 endif
endif


if(theory=='dft_band')then
   call write_band(ctx%system,ctx%energy)
end if

! output the wavefunctions for next GS calculations
if(write_gs_wfn_k == 'y') then !this input keyword is going to be removed....
   select case(iperiodic)
   case(3)
      call write_wfn(ctx%lg,ctx%mg,ctx%spsi,ctx%info,ctx%system)
      ! Experimental Implementation of Inner-Product Outputs:
      ! call write_prod_dk_data(lg, mg, system, info, spsi)
   case(0)
      write(*,*) "error: write_gs_wfn_k='y' & iperiodic=0"
   end select
end if

! output transition moment : --> want to put out of the optmization loop in future
if(yn_out_tm  == 'y'.or. yn_out_tm_bin == 'y'.or.yn_out_gs_sgm_eps=='y') then
   select case(iperiodic)
   case(3)
      call write_k_data(ctx%system,ctx%stencil)  !need? (probably remove later)
      call write_tm_data(ctx%spsi,ctx%system,ctx%info,ctx%mg,ctx%stencil,ctx%srg,ctx%ppg,ctx%energy)
   case(0)
     write(*,*) "error: yn_out_tm='y','yn_out_tm_bin='y',yn_out_gs_sgm_eps='y' & iperiodic=0"
  end select
end if

   ! force
   if(yn_jm=='n' .and. yn_dc=="n")then
     call calc_force(ctx%system,ctx%pp,ctx%fg,ctx%info,ctx%mg,ctx%stencil,ctx%poisson,ctx%srg,ctx%ppg,ctx%spsi,ctx%ewald,ctx%Vxc)
     if(comm_is_root(nproc_id_global))then
        if(theory=='dft'.and.yn_opt=='n'.and.yn_md=='n')then
         open(newunit=force_status_unit,file=resolve_path('gs_forces.dat'),status='replace')
         write(force_status_unit,'(a)')'SALMONPY_FORCES_V1'
         write(force_status_unit,*)natom
         do iatom=1,natom
          write(force_status_unit,'(3es26.17e3)')ctx%system%Force(:,iatom)
         enddo
         close(force_status_unit)
        endif
        write(*,*) "===== force ====="
        do iatom=1,natom
           select case(unit_system)
           case('au','a.u.'); write(*,300)iatom,(ctx%system%Force(ix,iatom),ix=1,3)
           case('A_eV_fs'  ); write(*,300)iatom,(ctx%system%Force(ix,iatom)*au_energy_ev/au_length_aa,ix=1,3)
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
   call structure_opt_check(iopt,flag_opt_conv,ctx%system%Force)
   if(.not.flag_opt_conv) call structure_opt(ctx%opt,iopt,ctx%system)
   !! Rion is old variables to be removed
   !! but currently it is used in many subroutines.
   !!Rion(:,:) = system%Rion(:,:)

   write(comment_line,10) iopt
   call write_xyz(comment_line,"add","r  ",ctx%system,ctx%ofl)
10 format("#opt iteration step=",i5)

   if(comm_is_root(nproc_id_global))then
      write(*,*) "atomic coordinate"
      do iatom=1,natom
         write(*,20) "'"//trim(atom_name(iatom))//"'",  &
                   (ctx%system%Rion(jj,iatom)*ulength_from_au,jj=1,3), &
                   Kion(iatom), flag_opt_atom(iatom)
      end do
20    format(a5,3f16.8,i3,a3)
   end if

   if(flag_opt_conv) then
      call structure_opt_fin(ctx%opt)
   else
      is_checkpoint_iter = (checkpoint_interval >= 1) .and. (mod(iopt,checkpoint_interval) == 0)
      is_shutdown_time   = (time_shutdown > 0d0) .and. (adjust_elapse_time(timer_now(LOG_TOTAL)) > time_shutdown)

      if(is_checkpoint_iter .or. is_shutdown_time) then
         if (is_shutdown_time .and. comm_is_root(ctx%info%id_rko)) then
           print *, 'shutdown the calculation, iopt =', iopt
         end if

         call checkpoint_gs(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%spsi,iopt,ctx%mixing)
         call comm_sync_all
         call checkpoint_opt(iopt,ctx%opt)
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

if(yn_dc=='y') then
  if(yn_dc_lcfo == 'y') call dc_lcfo(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%stencil,ctx%ppg,ctx%energy, &
    ctx%v_local,ctx%spsi,ctx%shpsi,ctx%sttpsi,ctx%srg,ctx%dc)
  call write_total_dcdft(ctx%system,ctx%dc)
end if

! write GS: basic data
if(yn_dc=='n') call write_band_information(ctx%system,ctx%energy)
call write_eigen(ctx%ofl,ctx%system,ctx%energy)
call write_info_data(Miter,ctx%system,ctx%energy,ctx%pp)
call write_k_data(ctx%system,ctx%stencil)
if(yn_spinorbit=='y') call write_mag_decomposed_gs(ctx%system,ctx%mg,ctx%info,ctx%spsi)

! write GS: analysis option
if(yn_out_psi =='y') call write_psi(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%spsi)
if(yn_out_dns =='y') call write_dns(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%rho_s)
if(yn_out_dos =='y') call write_dos(ctx%system,ctx%energy)
if(yn_out_pdos=='y') call write_pdos(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%pp,ctx%energy,ctx%spsi)
if(yn_out_elf =='y') call write_elf(0,ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%stencil,ctx%rho,ctx%srg,ctx%srg_scalar,ctx%spsi)

! if(comm_is_root(nproc_id_global)) call verify_atomic_orbital_consistency(system,pp,lg,mg,nproc_id_global) !!DEBUG

call timer_end(LOG_WRITE_GS_RESULTS)

! write GS: binary data for restart
call timer_begin(LOG_WRITE_GS_DATA)
if(write_gs_restart_data=="no") then
   if(comm_is_root(nproc_id_global)) &
      write(*,'(a)')"  no restart data writing."
else if(write_gs_restart_data.ne."checkpoint_only") then
   if(comm_is_root(nproc_id_global)) write(*,'(a)')"  writing restart data..."
   call checkpoint_gs(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%spsi,Miter,ctx%mixing,ctx%ofl%dir_out_restart)
   call comm_sync_all
   if(yn_opt=='y') then
      if(.not.flag_opt_conv) then
         call comm_sync_all
         call checkpoint_opt(nopt_max,ctx%opt,ctx%ofl%dir_out_restart)
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
   call checkpoint_gs(ctx%lg,ctx%mg,ctx%system,ctx%info,ctx%spsi,Miter,ctx%mixing)
   call comm_sync_all
endif
if(comm_is_root(nproc_id_global)) write(*,'(a)')"  writing completed."
call timer_end(LOG_WRITE_GS_DATA)

!call timer_begin(LOG_WRITE_GS_INFO)  !if needed, please take back, sory: AY
!call timer_end(LOG_WRITE_GS_INFO)

! Owned XC/HSE resources remain live until gs_close.

if(yn_dc=='y') then
! override (restore)
  nproc_group_global = ctx%dc%icomm_tot
  nproc_id_global = ctx%dc%id_tot
  nproc_size_global = ctx%dc%isize_tot
  call comm_sync_all
end if

call timer_end(LOG_TOTAL)


ctx%iterations=Miter;ctx%residual=sum1;ctx%threshold=selected%threshold
#ifdef USE_MPI
if(ctx%mpi_child_owner>0)call comm_end_children(ctx%mpi_child_owner)
#endif
status=0
end subroutine gs_run

subroutine gs_close(ctx)
 type(s_gs_context),target,intent(inout)::ctx
 type(s_gs_context)::empty
 if(.not.ctx%active)return
 call dealloc_cache(ctx%srg)
 call dealloc_cache(ctx%srg_scalar)
 call finalize_xc(ctx%xc_func)
#ifdef USE_HSE
 if(ctx%hse_bound)then
  call exx_unbind_system(ctx%system)
  call rt_state_clear(ctx%hse)
 endif
#endif
#ifdef USE_MPI
 call comm_release_children(ctx%mpi_child_owner)
#endif
 ctx=empty
end subroutine
subroutine gs_export_seed(ctx,seed,status)
 type(s_gs_context),intent(in)::ctx
 type(s_gs_seed),intent(out)::seed
 integer,intent(out)::status
 status=1
 if(.not.ctx%active)return
 seed%inputs=ctx%inputs
 seed%psi=ctx%spsi
 seed%occupations=ctx%system%rocc
 seed%grid=ctx%lg%num;seed%core_start=ctx%mg%is;seed%core_end=ctx%mg%ie
 seed%decomposition=[ctx%info%npk,ctx%info%nporbital,ctx%info%nprgrid]
 seed%orbital_start=ctx%info%io_s;seed%orbital_end=ctx%info%io_e
 seed%k_start=ctx%info%ik_s;seed%k_end=ctx%info%ik_e
 status=0
end subroutine
subroutine gs_probe_write(ctx,unit)
 type(s_gs_context),intent(in)::ctx
 integer,intent(in)::unit
 integer::i,j
 if(.not.ctx%active)error stop 'GS probe requires a live context'
 if(ctx%imported)then
  ! Initial import has no converged energy/residual/force report.
  write(unit)ctx%iterations,ctx%spsi%zwf,ctx%rho%f,ctx%vh%f,ctx%vpsl%f,ctx%system%rocc
  do i=1,size(ctx%vxc)
   write(unit)ctx%vxc(i)%f,ctx%v_local(i)%f,ctx%rho_s(i)%f
  enddo
  return
 endif
 write(unit)ctx%iterations,ctx%energy%E_tot,ctx%energy%E_xc
 write(unit)ctx%spsi%zwf,ctx%rho%f,ctx%vh%f,ctx%vpsl%f
 write(unit)ctx%system%rocc,ctx%system%Force
 write(unit)ctx%system%occupation_input,ctx%system%occupation_input_ready
 write(unit)ctx%system%gs_input%normalization
 write(unit)ctx%mixing%input,ctx%mixing%input_owned,ctx%mixing%num_rho_stock,ctx%mixing%flag_mix_zero
 write(unit)ctx%mixing%mixrate,ctx%mixing%alpha_mb,ctx%mixing%beta_p,ctx%mixing%convergence_value_prev
 if(allocated(ctx%mixing%rho_in))then
  do i=1,size(ctx%mixing%rho_in)
   if(allocated(ctx%mixing%rho_in(i)%f))write(unit)ctx%mixing%rho_in(i)%f
  enddo
 endif
 if(allocated(ctx%mixing%rho_out))then
  do i=1,size(ctx%mixing%rho_out)
   if(allocated(ctx%mixing%rho_out(i)%f))write(unit)ctx%mixing%rho_out(i)%f
  enddo
 endif
 if(allocated(ctx%mixing%Vh_in))then
  do i=1,size(ctx%mixing%Vh_in)
   if(allocated(ctx%mixing%Vh_in(i)%f))write(unit)ctx%mixing%Vh_in(i)%f
  enddo
 endif
 if(allocated(ctx%mixing%Vh_out))then
  do i=1,size(ctx%mixing%Vh_out)
   if(allocated(ctx%mixing%Vh_out(i)%f))write(unit)ctx%mixing%Vh_out(i)%f
  enddo
 endif
 if(allocated(ctx%mixing%rho_s_in))then
  do j=1,size(ctx%mixing%rho_s_in,2)
   do i=1,size(ctx%mixing%rho_s_in,1)
    if(allocated(ctx%mixing%rho_s_in(i,j)%f))write(unit)ctx%mixing%rho_s_in(i,j)%f
   enddo
  enddo
 endif
 if(allocated(ctx%mixing%rho_s_out))then
  do j=1,size(ctx%mixing%rho_s_out,2)
   do i=1,size(ctx%mixing%rho_s_out,1)
    if(allocated(ctx%mixing%rho_s_out(i,j)%f))write(unit)ctx%mixing%rho_s_out(i,j)%f
   enddo
  enddo
 endif
 if(allocated(ctx%mixing%Vxc_in))then
  do j=1,size(ctx%mixing%Vxc_in,2)
   do i=1,size(ctx%mixing%Vxc_in,1)
    if(allocated(ctx%mixing%Vxc_in(i,j)%f))write(unit)ctx%mixing%Vxc_in(i,j)%f
   enddo
  enddo
 endif
 if(allocated(ctx%mixing%Vxc_out))then
  do j=1,size(ctx%mixing%Vxc_out,2)
   do i=1,size(ctx%mixing%Vxc_out,1)
    if(allocated(ctx%mixing%Vxc_out(i,j)%f))write(unit)ctx%mixing%Vxc_out(i,j)%f
   enddo
  enddo
 endif
 do i=1,size(ctx%vxc)
  write(unit)ctx%vxc(i)%f,ctx%v_local(i)%f,ctx%rho_s(i)%f
 enddo
#ifdef USE_HSE
 if(allocated(ctx%hse%cached_source))write(unit)ctx%hse%cached_source
 if(allocated(ctx%hse%cached_action))write(unit)ctx%hse%cached_action
 if(allocated(ctx%hse%ace%factors))write(unit)ctx%hse%ace%factors
#endif
end subroutine
#ifdef USE_MPI
subroutine gs_seed_from_handle(token,seed,status)
 integer(c_int64_t),intent(in)::token
 type(s_gs_seed),intent(out)::seed
 integer,intent(out)::status
 integer::index
 call gs_lookup(token,index,status)
 if(status==0)call gs_export_seed(gs_handles(index)%context,seed,status)
end subroutine
logical function scf_supported() result(ok)
 use salmon_global,only:theory,yn_periodic,spin,yn_spinorbit,yn_dc,yn_opt,yn_md,yn_symmetry,xc,yn_jm,xname,cname
 use density_matrix_and_energy_plusU_sub,only:PLUS_U_ON
 ok=.false.
#ifdef USE_OPENACC
 return
#endif
 if(theory/='dft'.or.yn_periodic/='y'.or.spin/='unpolarized')return
 if(yn_spinorbit/='n'.or.yn_dc/='n'.or.yn_opt/='n'.or.yn_md/='n')return
 if(PLUS_U_ON.or.yn_jm/='n'.or.xname/='none'.or.cname/='none')return
 select case(trim(xc))
 case('pz','libxc_pz','libxc_pbe')
  ok=.true.
#ifdef USE_HSE
 case('hse06')
  ok=.true.
#endif
 end select
end function
logical function scf_context_ok(ctx) result(ok)
 type(s_gs_context),intent(in)::ctx
 ok=.false.
 if(.not.ctx%active.or.ctx%imported.or..not.ctx%scf_capable)return
#ifdef USE_HSE
 if(trim(ctx%system%xc_input%xcname)=='hse06')then
  if(.not.ctx%hse_bound)return
 endif
#endif
 ok=.true.
end function
integer(c_int) function sp_gs_scf_protocol() bind(C)
 sp_gs_scf_protocol=1
end function
integer(c_int) function sp_gs_start(path,steps,token) bind(C)
 character(c_char),intent(in)::path(*)
 integer(c_int),value::steps
 integer(c_int64_t),intent(out)::token
 if(steps<1)then
  token=0;sp_gs_start=11;return
 endif
 sp_gs_start=create_gs(path,token,.false.,steps)
end function
integer(c_int) function sp_gs_scf_info(token,values) bind(C)
 integer(c_int64_t),value::token
 integer(c_int),intent(out)::values(6)
 integer::i,index,status,ns
 values=0
 call gs_lookup(token,index,status);sp_gs_scf_info=status
 if(status/=0)return
 associate(ctx=>gs_handles(index)%context)
  if(.not.scf_context_ok(ctx))then
   sp_gs_scf_info=15;return
  endif
  ns=0
  if(allocated(ctx%scf_snapshots))then
   do i=1,size(ctx%scf_snapshots)
    if(allocated(ctx%scf_snapshots(i)%value))ns=ns+1
   enddo
  endif
  values=[merge(1,0,ctx%mixing%scf_progress%paused),ctx%iterations, &
          ctx%system%gs_input%max_iterations,ctx%initializations,ctx%executed_iterations,ns]
 end associate
end function
integer(c_int) function sp_gs_continue(token,steps) bind(C)
 use salmon_global,only:calc_mode
 use timer
 use force_sub,only:calc_force
 use sym_sub,only:swap_hse_symmetry_workspace
 use sym_rho_sub,only:swap_hse_density_workspace
 integer(c_int64_t),value::token
 integer(c_int),value::steps
 type(s_model_inputs)::previous
 type(s_timer_context)::host_timers
 type(s_hse_symmetry_workspace)::previous_symmetry
 type(s_hse_density_workspace)::previous_density_symmetry
 integer::index,status,old_iter,group,rank,nproc,reference_steps,restore_status
 real(8)::residual
 character(len(calc_mode))::old_mode
 logical :: rion_update
 call gs_lookup(token,index,status);sp_gs_continue=status
 if(status/=0)return
 call comm_get_globalinfo(group,rank,nproc)
 reference_steps=steps;call comm_bcast(reference_steps,group)
 status=0
 if(steps<1.or.steps/=reference_steps)status=11
 if(.not.scf_context_ok(gs_handles(index)%context))status=15
 if(allocated(gs_handles(index)%context%scf_pending))status=11
 call comm_get_max(status,group)
 sp_gs_continue=status
 if(status/=0)return
 associate(ctx=>gs_handles(index)%context)
  if(ctx%converged.or.ctx%iterations>=ctx%system%gs_input%max_iterations)return
  call previous%capture()
  call ctx%continuation_inputs%activate(status)
  call comm_get_max(status,group)
  if(status/=0)then
   call previous%activate(restore_status)
   sp_gs_continue=status;return
  endif
  call swap_hse_symmetry_workspace(previous_symmetry)
  call swap_hse_density_workspace(previous_density_symmetry)
  call swap_hse_symmetry_workspace(ctx%symmetry)
  call swap_hse_density_workspace(ctx%density_symmetry)
  old_mode=calc_mode;calc_mode='GS'
  old_iter=ctx%iterations
  ctx%mixing%scf_progress%continuing=.true.
  ! Avoid overflow; an oversized budget simply reaches the original total cap.
  ctx%mixing%scf_progress%stop_after=ctx%iterations+min(steps, &
    ctx%system%gs_input%max_iterations-ctx%iterations)
  call timer_context_swap(host_timers)
  call timer_context_swap(ctx%timers)
  call timer_initialize
  call timer_begin(LOG_TOTAL)
  call timer_enable_sub
  call timer_begin(LOG_GS_ITERATION)
  rion_update=ctx%mixing%scf_progress%rion_update
  call scf_iteration_dft(ctx%iterations,rion_update,residual, &
    ctx%system,ctx%energy,ctx%ewald,ctx%lg,ctx%mg,ctx%info,ctx%poisson,ctx%fg, &
    ctx%cg,ctx%mixing,ctx%stencil,ctx%srg,ctx%srg_scalar,ctx%spsi,ctx%shpsi,ctx%sttpsi, &
    ctx%rho,ctx%rho_jm,ctx%rho_s,ctx%V_local,ctx%Vh,ctx%Vxc,ctx%Vpsl,ctx%xc_func, &
    ctx%pp,ctx%ppg,ctx%ppn,ctx%band,3,ctx%dc)
  ctx%residual=residual;ctx%converged=ctx%mixing%scf_converged
  call ctx%continuation_inputs%capture()
  ctx%executed_iterations=ctx%executed_iterations+ctx%iterations-old_iter
  if(.not.ctx%mixing%scf_progress%paused)then
   call calc_force(ctx%system,ctx%pp,ctx%fg,ctx%info,ctx%mg,ctx%stencil,ctx%poisson, &
                   ctx%srg,ctx%ppg,ctx%spsi,ctx%ewald,ctx%Vxc)
  endif
  call timer_end(LOG_GS_ITERATION)
  call timer_disable_sub
  call timer_end(LOG_TOTAL)
  call timer_context_swap(ctx%timers)
  call timer_context_swap(host_timers)
  calc_mode=old_mode
  call swap_hse_density_workspace(ctx%density_symmetry)
  call swap_hse_symmetry_workspace(ctx%symmetry)
  call swap_hse_density_workspace(previous_density_symmetry)
  call swap_hse_symmetry_workspace(previous_symmetry)
  call previous%activate(restore_status)
  call comm_get_max(restore_status,group)
  sp_gs_continue=restore_status
 end associate
end function
subroutine scf_copy(ctx,value,restore)
 type(s_gs_context),target,intent(inout)::ctx
 type(s_scf_payload),intent(inout)::value
 logical,intent(in)::restore
 if(restore)then
  ctx%continuation_inputs=value%continuation_inputs
  ctx%system=value%system;ctx%spsi=value%spsi;ctx%shpsi=value%shpsi;ctx%sttpsi=value%sttpsi
  ctx%cg=value%cg;ctx%mixing=value%mixing;ctx%rho=value%rho;ctx%rho_jm=value%rho_jm
  ctx%Vh=value%Vh;ctx%Vpsl=value%Vpsl;ctx%V_local=value%V_local;ctx%rho_s=value%rho_s;ctx%Vxc=value%Vxc
  ctx%poisson=value%poisson;ctx%fg=value%fg;ctx%energy=value%energy
  ctx%iterations=value%iterations;ctx%converged=value%converged
  ctx%residual=value%residual;ctx%threshold=value%threshold;ctx%xc_func%cval=value%cval
#ifdef USE_HSE
  if(ctx%hse_bound)then
   call rt_values_copy(ctx%hse,value%hse,.true.)
   ctx%system%salmonpy_hse_owner=>ctx%hse
  endif
#endif
 else
  value%continuation_inputs=ctx%continuation_inputs
  value%system=ctx%system;value%spsi=ctx%spsi;value%shpsi=ctx%shpsi;value%sttpsi=ctx%sttpsi
  value%cg=ctx%cg;value%mixing=ctx%mixing;value%rho=ctx%rho;value%rho_jm=ctx%rho_jm
  value%Vh=ctx%Vh;value%Vpsl=ctx%Vpsl;value%V_local=ctx%V_local;value%rho_s=ctx%rho_s;value%Vxc=ctx%Vxc
  value%poisson=ctx%poisson;value%fg=ctx%fg;value%energy=ctx%energy
  value%iterations=ctx%iterations;value%converged=ctx%converged
  value%residual=ctx%residual;value%threshold=ctx%threshold;value%cval=ctx%xc_func%cval
  ! Never retain an alias back into a live context in the numerical payload.
  nullify(value%system%salmonpy_hse_owner)
#ifdef USE_HSE
  if(ctx%hse_bound)call rt_values_copy(ctx%hse,value%hse,.false.)
#endif
 endif
end subroutine
integer(c_int) function sp_gs_scf_snapshot(token,snapshot) bind(C)
 integer(c_int64_t),value::token
 integer(c_int64_t),intent(out)::snapshot
 type(s_scf_snapshot_slot),allocatable::larger(:)
 integer::index,status,i,n,j
 type(s_scf_payload),allocatable::saved
 snapshot=0
 call gs_lookup(token,index,status);sp_gs_scf_snapshot=status
 if(status/=0)return
 associate(ctx=>gs_handles(index)%context)
  if(.not.scf_context_ok(ctx))then
   sp_gs_scf_snapshot=15;return
  endif
  if(allocated(ctx%scf_pending))then
   sp_gs_scf_snapshot=11;return
  endif
  n=0;i=0
  if(allocated(ctx%scf_snapshots))then
   n=size(ctx%scf_snapshots)
   do i=1,n
    if(.not.allocated(ctx%scf_snapshots(i)%value))exit
   enddo
   if(i>n)i=0
  endif
  if(i==0)then
   allocate(larger(n+1))
   do j=1,n
    larger(j)%token=ctx%scf_snapshots(j)%token
    if(allocated(ctx%scf_snapshots(j)%value))call move_alloc(ctx%scf_snapshots(j)%value,larger(j)%value)
   enddo
   call move_alloc(larger,ctx%scf_snapshots)
   i=n+1
  endif
  scf_snapshot_serial=scf_snapshot_serial+1
  snapshot=scf_snapshot_serial
  allocate(saved)
  call scf_copy(ctx,saved,.false.)
  ctx%scf_snapshots(i)%token=snapshot
  call move_alloc(saved,ctx%scf_snapshots(i)%value)
 end associate
end function
integer(c_int) function scf_snapshot_action(token,snapshot,restore) result(status)
 integer(c_int64_t),intent(in)::token,snapshot
 logical,intent(in)::restore
 integer::index,i,group,rank,nproc,parts(2),expected(2),slot,local
 type(s_scf_payload)::candidate
 call gs_lookup(token,index,status)
 if(status/=0)return
 call comm_get_globalinfo(group,rank,nproc)
 parts=transfer(snapshot,parts);expected=parts;call comm_bcast(expected,group)
 slot=0;local=0
 associate(ctx=>gs_handles(index)%context)
  if(all(parts==expected).and.allocated(ctx%scf_snapshots))then
   do i=1,size(ctx%scf_snapshots)
    if(ctx%scf_snapshots(i)%token/=snapshot)cycle
    if(.not.allocated(ctx%scf_snapshots(i)%value))cycle
    slot=i;exit
   enddo
  endif
  if(slot==0.or.allocated(ctx%scf_pending))local=11
  call comm_get_max(local,group)
  if(local/=0)then
   status=local;return
  endif
  if(restore)then
   ! Validate all selections before writes; allocations are not a general OOM recovery contract.
   candidate=ctx%scf_snapshots(slot)%value
   call scf_copy(ctx,candidate,.true.)
  else
   deallocate(ctx%scf_snapshots(slot)%value)
   ctx%scf_snapshots(slot)%token=0
  endif
 end associate
end function
integer(c_int) function sp_gs_scf_restore(token,snapshot) bind(C)
 integer(c_int64_t),value::token,snapshot
 sp_gs_scf_restore=scf_snapshot_action(token,snapshot,.true.)
end function
integer(c_int) function sp_gs_scf_snapshot_close(token,snapshot) bind(C)
 integer(c_int64_t),value::token,snapshot
 sp_gs_scf_snapshot_close=scf_snapshot_action(token,snapshot,.false.)
end function
integer(c_int) function sp_gs_scf_candidate_protocol() bind(C)
 sp_gs_scf_candidate_protocol=1
end function
integer(c_int) function sp_gs_scf_candidate(token,snapshot,mode) bind(C) result(status)
 ! Numerical candidates only: resource owners stay live and are never cloned.
 ! Mode2 stages, mode3 publishes after the host ready vote, mode4 discards.
 integer(c_int64_t),value::token,snapshot
 integer(c_int),value::mode
 integer::index,group,rank,nproc,expected_mode,parts(2),expected(2),slot,i,local
 type(s_scf_payload),allocatable::candidate
 call gs_lookup(token,index,status)
 if(status/=0)return
 call comm_get_globalinfo(group,rank,nproc)
 expected_mode=mode;call comm_bcast(expected_mode,group)
 parts=transfer(snapshot,parts);expected=parts;call comm_bcast(expected,group)
 local=0;slot=0
 if(mode/=expected_mode.or.any(parts/=expected))local=11
 associate(ctx=>gs_handles(index)%context)
  if(.not.scf_context_ok(ctx))local=15
  select case(mode)
  case(2)
   if(allocated(ctx%scf_pending))local=11
   if(allocated(ctx%scf_snapshots))then
    do i=1,size(ctx%scf_snapshots)
     if(ctx%scf_snapshots(i)%token/=snapshot)cycle
     if(.not.allocated(ctx%scf_snapshots(i)%value))cycle
     slot=i;exit
    enddo
   endif
   if(slot==0)local=11
  case(3)
   if(.not.allocated(ctx%scf_pending).or.snapshot/=0_c_int64_t)local=11
  case(4)
   if(snapshot/=0_c_int64_t)local=11
  case default
   local=11
  end select
  call comm_get_max(local,group)
  status=local
  if(status/=0)return
  select case(mode)
  case(2)
   allocate(candidate)
   candidate=ctx%scf_snapshots(slot)%value
   call move_alloc(candidate,ctx%scf_pending)
  case(3)
   ! No parsing, resource initialization, orbital solve or exchange evaluation.
   ! Allocation/process failure during assignments is outside handled recovery.
   call move_alloc(ctx%scf_pending,candidate)
   call scf_copy(ctx,candidate,.true.)
  case(4)
   if(allocated(ctx%scf_pending))deallocate(ctx%scf_pending)
  end select
 end associate
end function
integer(c_int) function sp_gs_protocol() bind(C)
 sp_gs_protocol=1
end function
integer(c_int) function sp_gs_compute_protocol() bind(C)
 sp_gs_compute_protocol=1
end function
integer(c_int) function sp_gs_create(path,token) bind(C)
 character(c_char),intent(in)::path(*)
 integer(c_int64_t),intent(out)::token
 token=0;sp_gs_create=2
 if(sp_gs_compute_protocol()/=1)return
 sp_gs_create=create_gs(path,token,.false.)
end function
subroutine gs_lookup(token,index,status)
 integer(c_int64_t),intent(in)::token
 integer,intent(out)::index,status
 integer::group,rank,nproc,i,parts(2),expected(2)
 index=0;status=10
 if(.not.comm_is_attached())return
 call comm_get_globalinfo(group,rank,nproc)
 parts=transfer(token,parts);expected=parts
 call comm_bcast(expected,group)
 if(all(expected==parts).and.allocated(gs_handles))then
  do i=1,size(gs_handles)
   if(gs_handles(i)%token/=token.or..not.associated(gs_handles(i)%context))cycle
   index=i;status=0;exit
  enddo
 endif
 call comm_get_max(status,group)
 if(status/=0)index=0
end subroutine
integer(c_int) function sp_gs_import(path,token) bind(C)
 character(c_char),intent(in)::path(*)
 integer(c_int64_t),intent(out)::token
 sp_gs_import=create_gs(path,token,.true.)
end function
integer(c_int) function create_gs(path,token,reading_only,steps)
 use inputoutput,only:read_input,s_input_parser
 use salmon_global,only:yn_symmetry,calc_mode
 use plusU_global,only:read_Hubbard_parameters
 use sym_sub,only:s_hse_symmetry_workspace,swap_hse_symmetry_workspace,read_sw_symmetry
 use sym_rho_sub,only:s_hse_density_workspace,swap_hse_density_workspace
 use math_constants,only:set_math_constants
 use parallelization,only:setup_parallel
 use timer,only:timer_initialize
 logical,intent(in)::reading_only
 integer,optional,intent(in)::steps
 character(c_char),intent(in)::path(*)
 integer(c_int64_t),intent(out)::token
 type(s_input_parser)::parser
 type(s_timer_context)::host_timers
 type(s_hse_symmetry_workspace)::old_symmetry
 type(s_hse_density_workspace)::old_density_symmetry
 type(s_gs_context),pointer::candidate
 type(gs_handle_slot),allocatable::enlarged(:)
 character(len(calc_mode))::previous_mode
 character(4096)::filename,expected
 integer::group,rank,nproc,i,index,status,selection(2),root_selection(2)
 token=0;create_gs=10
 if(.not.comm_is_attached())return
 call comm_get_globalinfo(group,rank,nproc)
 filename='';status=0
 do i=1,len(filename)
  if(path(i)==c_null_char)exit
  filename(i:i)=path(i)
 enddo
 if(i>len(filename).or.len_trim(filename)==0)status=11
 expected=filename;call comm_bcast(expected,group)
 if(filename/=expected)status=11
 selection=0
 if(present(steps))selection=[1,steps]
 root_selection=selection;call comm_bcast(root_selection,group)
 if(any(selection/=root_selection))status=11
 if(present(steps))then
  if(steps<0)status=11
 endif
 call comm_get_max(status,group)
 if(status/=0)then
  create_gs=status;return
 endif
 if(.not.gs_host_initialized)then
  call setup_parallel
  call set_math_constants
  call timer_initialize
  gs_host_initialized=.true.
 endif
 call read_input(trim(filename),parser)
 call read_Hubbard_parameters
 allocate(candidate)
 candidate%scf_capable=scf_supported()
 call swap_hse_symmetry_workspace(old_symmetry)
 call swap_hse_density_workspace(old_density_symmetry)
 call read_sw_symmetry(yn_symmetry)
 previous_mode=calc_mode;calc_mode='GS'
 if(present(steps))then
  if(steps<0.or..not.scf_supported())then
   calc_mode=previous_mode
   call swap_hse_density_workspace(old_density_symmetry)
   call swap_hse_symmetry_workspace(old_symmetry)
   deallocate(candidate);create_gs=15;return
  endif
  candidate%mixing%scf_progress%stop_after=steps
 endif
 call timer_context_swap(host_timers)
 call timer_context_swap(candidate%timers)
 call timer_initialize
 if(present(steps))then
  call gs_run(candidate,status,import_only=reading_only,disk_bootstrap=steps==0)
 else
  call gs_run(candidate,status,import_only=reading_only)
 endif
 call timer_context_swap(candidate%timers)
 call timer_context_swap(host_timers)
 call swap_hse_density_workspace(candidate%density_symmetry)
 call swap_hse_symmetry_workspace(candidate%symmetry)
 calc_mode=previous_mode
 call swap_hse_density_workspace(old_density_symmetry)
 call swap_hse_symmetry_workspace(old_symmetry)
 call comm_get_max(status,group)
 if(status/=0)then
  call gs_close(candidate);deallocate(candidate)
  create_gs=status;return
 endif
 index=0
 if(allocated(gs_handles))then
  do i=1,size(gs_handles)
   if(.not.associated(gs_handles(i)%context))then
    index=i;exit
   endif
  enddo
 endif
 if(index==0)then
  i=0
  if(allocated(gs_handles))i=size(gs_handles)
  allocate(enlarged(i+1))
  if(i>0)enlarged(:i)=gs_handles
  call move_alloc(enlarged,gs_handles)
  index=i+1
 endif
 gs_handle_serial=gs_handle_serial+1
 token=gs_handle_prefix+gs_handle_serial
 gs_handles(index)%context=>candidate
 gs_handles(index)%token=token
 create_gs=0
end function
integer(c_int) function sp_gs_timer_protocol() bind(C)
 sp_gs_timer_protocol=1
end function
integer(c_int) function sp_gs_timers(token,seconds) bind(C) result(status)
 integer(c_int64_t),value :: token
 real(c_double),intent(out) :: seconds(4)
 integer :: index
 seconds=0d0
 call gs_lookup(token,index,status)
 if(status/=0)return
 associate(owner=>gs_handles(index)%context%timers)
  seconds=[timer_context_get(owner,LOG_TOTAL),timer_context_get(owner,LOG_INIT_GS), &
           timer_context_get(owner,LOG_GS_ITERATION),timer_context_get(owner,LOG_WRITE_GS_RESULTS)]
 end associate
end function
integer(c_int) function sp_gs_info(token,stats,values) bind(C)
 integer(c_int64_t),value::token
 integer(c_int),intent(out)::stats(4)
 real(c_double),intent(out)::values(3)
 integer::index,status
 stats=0;values=0d0
 call gs_lookup(token,index,status)
 sp_gs_info=status
 if(status/=0)return
 associate(s=>gs_handles(index)%context)
  stats=[s%iterations,merge(1,0,s%converged),s%system%nion,s%system%no]
  if(s%imported)stats(2)=-1
  if(.not.s%imported)values=[s%energy%E_tot,s%residual,s%threshold]
 end associate
end function
subroutine array_layout(s,quantity,ndim,shape_out,lower,status)
  type(s_gs_context),intent(in) :: s
  integer,intent(in) :: quantity
  integer(c_int),intent(out) :: ndim
  integer(c_int64_t),intent(out) :: shape_out(7),lower(7)
  integer,intent(out) :: status
  ndim=0;shape_out=1;lower=1;status=0
  select case(quantity)
  case(1)
    ndim=3;shape_out(:3)=shape(s%rho%f,kind=c_int64_t);lower(:3)=lbound(s%rho%f,kind=c_int64_t)
  case(2)
    ndim=3;shape_out(:3)=shape(s%system%rocc,kind=c_int64_t);lower(:3)=lbound(s%system%rocc,kind=c_int64_t)
  case(3)
    if(s%imported)then
      status=11;return
    endif
    ndim=3;shape_out(:3)=shape(s%energy%esp,kind=c_int64_t);lower(:3)=lbound(s%energy%esp,kind=c_int64_t)
  case(4)
    if(s%imported.or.s%mixing%scf_progress%paused)then
      status=11;return
    endif
    ndim=2;shape_out(:2)=shape(s%system%Force,kind=c_int64_t)
  case(5)
    ndim=7
    if(allocated(s%spsi%rwf))then
      shape_out=shape(s%spsi%rwf,kind=c_int64_t);lower=lbound(s%spsi%rwf,kind=c_int64_t)
    else
      shape_out=shape(s%spsi%zwf,kind=c_int64_t);lower=lbound(s%spsi%zwf,kind=c_int64_t)
    endif
  case(6)
    if(.not.allocated(s%mixing%scf_progress%residual_history))then
      status=11;return
    endif
    ndim=2;shape_out(:2)=[4_c_int64_t,int(s%mixing%scf_progress%residual_count,c_int64_t)]
  case default
    status=11
  end select
end subroutine

integer(c_int) function sp_gs_array_info(token,quantity,ndim,shape_out,lower) bind(C)
  integer(c_int64_t),value :: token
  integer(c_int),value :: quantity
  integer(c_int),intent(out) :: ndim
  integer(c_int64_t),intent(out) :: shape_out(7),lower(7)
  integer :: index,status
  call gs_lookup(token,index,status)
  if(status==0)call array_layout(gs_handles(index)%context,quantity,ndim,shape_out,lower,status)
  sp_gs_array_info=status
end function

integer(c_int) function sp_gs_read(token,quantity,real_ptr,imag_ptr,length) bind(C)
  integer(c_int64_t),value :: token,length
  integer(c_int),value :: quantity
  type(c_ptr),value :: real_ptr,imag_ptr
  real(c_double),pointer :: re(:),im(:)
  integer(c_int64_t) :: shape_out(7),lower(7),k
  integer(c_int) :: ndim
  integer :: index,status,i1,i2,i3,i4,i5,i6,i7,group,rank,nproc
  call gs_lookup(token,index,status)
  if(status/=0)then
    sp_gs_read=status;return
  endif
  call array_layout(gs_handles(index)%context,quantity,ndim,shape_out,lower,status)
  if(status==0)then
    if(length/=product(shape_out).or..not.c_associated(real_ptr))status=11
    if(quantity==5.and..not.c_associated(imag_ptr))status=11
  endif
  call comm_get_globalinfo(group,rank,nproc)
  call comm_get_max(status,group)
  sp_gs_read=status
  if(status/=0)return
  call c_f_pointer(real_ptr,re,[length])
  associate(s=>gs_handles(index)%context)
    select case(quantity)
    case(1)
      re=reshape(s%rho%f,[length])
    case(2)
      re=reshape(s%system%rocc,[length])
    case(3)
      re=reshape(s%energy%esp,[length])
    case(4)
      re=reshape(s%system%Force,[length])
    case(6)
      re=reshape(s%mixing%scf_progress%residual_history(:,:s%mixing%scf_progress%residual_count),[length])
    case(5)
      call c_f_pointer(imag_ptr,im,[length])
      k=0
      do i7=int(lower(7)),int(lower(7)+shape_out(7)-1)
      do i6=int(lower(6)),int(lower(6)+shape_out(6)-1)
      do i5=int(lower(5)),int(lower(5)+shape_out(5)-1)
      do i4=int(lower(4)),int(lower(4)+shape_out(4)-1)
      do i3=int(lower(3)),int(lower(3)+shape_out(3)-1)
      do i2=int(lower(2)),int(lower(2)+shape_out(2)-1)
      do i1=int(lower(1)),int(lower(1)+shape_out(1)-1)
        k=k+1
        if(allocated(s%spsi%rwf))then
          re(k)=s%spsi%rwf(i1,i2,i3,i4,i5,i6,i7);im(k)=0d0
        else
          re(k)=real(s%spsi%zwf(i1,i2,i3,i4,i5,i6,i7),c_double)
          im(k)=aimag(s%spsi%zwf(i1,i2,i3,i4,i5,i6,i7))
        endif
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
    end select
  end associate
end function

integer(c_int) function sp_gs_close(token) bind(C)
 integer(c_int64_t),value::token
 integer::index,status
 call gs_lookup(token,index,status)
 sp_gs_close=status
 if(status/=0)return
 if(allocated(gs_handles(index)%context%scf_pending))then
  sp_gs_close=11;return
 endif
 call gs_close(gs_handles(index)%context)
 deallocate(gs_handles(index)%context)
 nullify(gs_handles(index)%context)
end function
integer(c_int) function sp_gs_resources(counts) bind(C)
 integer(c_int),intent(out)::counts(3)
 integer::i,live,nslots
 counts=0
 if(allocated(gs_handles))then
  do i=1,size(gs_handles)
   if(associated(gs_handles(i)%context))counts(1)=counts(1)+1
  enddo
 endif
 call comm_child_counts(live,nslots)
 counts(2:3)=[live,nslots]
 sp_gs_resources=0
end function
#endif
#ifdef USE_MPI
#include "scf_disk_native.inc"
#include "scf_disk_io.inc"
#endif
end module salmonpy_gs_context

subroutine main_dft
 use salmon_global,only:theory,iperiodic
 use salmonpy_gs_context,only:s_gs_context,gs_run,gs_close
 implicit none
 type(s_gs_context),target::ctx
 integer::status
 if(theory=='dft_band'.and.iperiodic/=3)return
 call gs_run(ctx,status)
 if(status/=0)error stop 'Invalid GS run configuration'
 call gs_close(ctx)
end subroutine
