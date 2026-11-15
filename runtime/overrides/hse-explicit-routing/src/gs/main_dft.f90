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
use structures
use salmon_xc,only:finalize_xc
#ifdef USE_HSE
use exx_rt_state,only:s_exx_rt_state,rt_state_clear
use exx_native,only:exx_unbind_system
#endif
implicit none
private
public::s_gs_context,gs_run,gs_close,gs_probe_write,gs_input_config
! Internal retained GS storage. Actual contexts must have TARGET and a stable address.
! Never copy a live context; global inputs remain shared.
type s_gs_context
 private
 logical::active=.false.
 integer::iterations=0
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
contains
subroutine gs_input_config(value)
 use salmon_global,only:nscf,threshold,exx_pre_scf_threshold,exx_pre_scf_steps,convergence
 use mixing_sub,only:mixing_input_config
 type(s_gs_run_config),intent(out)::value
 call mixing_input_config(value%mixing)
 value%convergence_metric=convergence
 value%max_iterations=nscf;value%threshold=threshold
 value%pre_scf_threshold=exx_pre_scf_threshold;value%pre_scf_steps=exx_pre_scf_steps
end subroutine
subroutine gs_run(ctx,status,config)
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
#ifdef USE_HSE
 integer :: gs_hse_status
#endif
integer :: phase_status
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
if(present(config))then
 selected=config
else
 call gs_input_config(selected)
endif
select case(selected%convergence_metric)
case('rho_dne','norm_rho','norm_rho_dng','norm_pot','norm_pot_dng')
case default
 return
end select
if(.not.valid_mixing_input(selected%mixing))return
if(selected%mixing%method=='simple_dm'.and.yn_spinorbit/='y')return
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
if(theory=='dft_band'.and.iperiodic/=3) return
ctx%active=.true.
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
call exx_bind_gs_system(ctx%system,ctx%info,ctx%hse,gs_hse_status)
if(gs_hse_status==1)error stop 'GS HSE owner binding conflict'
ctx%hse_bound=gs_hse_status==0
if(ctx%hse_bound.and.comm_is_root(nproc_id_global))write(*,'(a)')'EXX_GS_OWNER explicit'
#endif
allocate( ctx%rho_s(ctx%system%nspin),ctx%v_local(ctx%system%nspin),ctx%vxc(ctx%system%nspin) )

call initialization1_dft( ctx%system, ctx%energy, ctx%stencil, ctx%fg, ctx%poisson,  &
                          ctx%lg, ctx%mg,   &
                          ctx%info,  &
                          ctx%srg, ctx%srg_scalar,  &
                          ctx%rho, ctx%rho_jm, ctx%rho_s, ctx%vh, ctx%v_local, ctx%vpsl, ctx%vxc,  &
                          ctx%spsi, ctx%shpsi, ctx%sttpsi,  &
                          ctx%pp, ctx%ppg, ctx%ppn,  &
                          ctx%ofl )

call initialization2_dft( Miter, nspin, rion_update,  &
                          ctx%system, ctx%energy, ctx%ewald, ctx%stencil, ctx%fg, ctx%poisson,&
                          ctx%lg, ctx%mg, ctx%info,   &
                          ctx%srg, ctx%srg_scalar,  &
                          ctx%rho, ctx%rho_jm, ctx%rho_s, ctx%vh,ctx%v_local, ctx%vpsl, ctx%vxc,  &
                          ctx%spsi, ctx%shpsi, ctx%sttpsi,  &
                          ctx%pp, ctx%ppg, ctx%ppn,   &
                          ctx%xc_func, ctx%mixing )


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
     call calc_force(ctx%system,ctx%pp,ctx%fg,ctx%info,ctx%mg,ctx%stencil,ctx%poisson,ctx%srg,ctx%ppg,ctx%spsi,ctx%ewald)
     if(comm_is_root(nproc_id_global))then
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


ctx%iterations=Miter
status=0
end subroutine gs_run

subroutine gs_close(ctx)
 type(s_gs_context),target,intent(inout)::ctx
 type(s_gs_context)::empty
 if(.not.ctx%active)return
 call finalize_xc(ctx%xc_func)
#ifdef USE_HSE
 if(ctx%hse_bound)then
  call exx_unbind_system(ctx%system)
  call rt_state_clear(ctx%hse)
 endif
#endif
 ctx=empty
end subroutine
subroutine gs_probe_write(ctx,unit)
 type(s_gs_context),intent(in)::ctx
 integer,intent(in)::unit
 integer::i,j
 if(.not.ctx%active)error stop 'GS probe requires a live context'
 write(unit)ctx%iterations,ctx%energy%E_tot,ctx%energy%E_xc
 write(unit)ctx%spsi%zwf,ctx%rho%f,ctx%vh%f,ctx%vpsl%f
 write(unit)ctx%system%rocc,ctx%system%Force
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
