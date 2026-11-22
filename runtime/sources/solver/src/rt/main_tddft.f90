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
use tddft_worker
use salmonpy_step_output
use dip, only: subdip,calc_dip
use density_matrix, only: calc_current
use hamiltonian, only: update_kvector_nonlocalpt
#ifdef USE_HSE
use exx_rt_state,only:s_exx_rt_state,rt_state_clear
use exx_native,only:exx_bind_system,exx_unbind_system
#endif
use math_constants, only: pi
use salmon_global
use structures
use parallelization, only: nproc_group_global
use communication, only: comm_is_root, comm_sync_all, comm_bcast
use salmon_xc, only: finalize_xc
use timer
use write_sub, only: write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d
use initialization_rt_sub
use checkpoint_restart_sub
use jellium, only: check_condition_jm
use nvtx_wrapper
use parallelization, only: nproc_id_global
implicit none

#ifdef USE_HSE
type(s_exx_rt_state),target :: independent_hse
character(8) :: context_mode
integer :: context_status
#endif
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
type(s_unfold) :: unfold

real(8) :: worker_A(3),worker_current(3,2),worker_ne,worker_density_integral(4)
integer :: Mit, itt
logical :: is_checkpoint_iter, is_shutdown_time, is_checkpoint

! Live-worker numerical savepoint. No owning MPI, files or Libxc handles copied.
type :: isolated_savepoint
 type(s_orbital)::spsi_in,spsi_out,tpsi
 type(s_scalar)::rho,rho_jm,Vh,Vh_stock1,Vh_stock2,Vbox,Vpsl
 type(s_scalar),allocatable::rho_s(:),V_local(:),Vxc(:)
 type(s_rt)::rt
 type(s_poisson)::poisson
 type(s_dft_energy)::energy
 type(s_xc_operator_payload)::xc
 real(8)::fields(3,4),cval
 integer::step
end type
type(isolated_savepoint),allocatable::worker_saved
type(step_output_state)::worker_outputs
integer::worker_operation,worker_status,worker_updates
logical::worker_outputs_closed

worker_updates=0
worker_outputs_closed=.false.

!check condition for using jellium model
if(yn_jm=='y') call check_condition_jm

call timer_begin(LOG_TOTAL)

#ifdef USE_HSE
! Explicit development opt-in; full-solver regression exercises this route.
context_mode=""
call get_environment_variable("SALMONPY_HSE_CONTEXT",context_mode,status=context_status)
if(context_status==0.and.trim(context_mode)=="1")then
  call exx_bind_system(system,independent_hse,context_status)
  if(context_status/=0)error stop "HSE owner binding failed"
endif
#endif
call initialization_rt( Mit, system, energy, ewald, rt, md, &
                        singlescale,  &
                        stencil, fg, poisson,  &
                        lg, mg,   &
                        info,  &
                        xc_func, ofl,  &
                        srg, srg_scalar,  &
                        spsi_in, spsi_out, tpsi, rho, rho_jm, rho_s,  &
                        V_local, Vbox, Vh, Vh_stock1, Vh_stock2, Vxc, Vpsl,&
                        pp, ppg, ppn, unfold )
#ifdef USE_HSE
if(trim(context_mode)=="1")then
  if(.not.associated(system%salmonpy_hse_owner,independent_hse))error stop "HSE owner lost during initialization"
  write(*,'(a)') "SALMONPY_HSE_NATIVE_CONTEXT active"
endif
#endif

if(worker_requested()) then
  if(worker_E_requested()) then
    call subdip(info%icomm_r,0,rt,lg,mg,system,rho_s,worker_ne,poisson,energy%E_tot,pp)
    call worker_reply_E(0,rt%dDp_e(:,0),rt%Dp_e(:,0)+rt%Dp_i(:,0),worker_ne,energy%E_tot,.true.)
  else
  system%vec_Ac=worker_initial_A
  call update_kvector_nonlocalpt(info%ik_s,info%ik_e,system,ppg)
  call calc_current(system,mg,stencil,info,srg,spsi_in,ppg,worker_current(:,1:system%nspin))
  spsi_in%update_zwf_overlap=.true.
  if(worker_norm_requested())then
    call calc_dip(info%icomm_r,lg,mg,system,rho_s,worker_density_integral)
    worker_electron_number=worker_density_integral(4)*system%Hvol
  endif
  call worker_reply(0,sum(worker_current(:,1:system%nspin),dim=2))
  end if
end if

#ifdef __FUJITSU
call fapp_start('time_evol',1,0) ! performance profiling
#endif

call print_header()

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
itt=Mit
TE : do while(itt<nt.or.worker_E_requested())
  if(worker_requested()) then
    if(.not.worker_next_A(worker_A,worker_operation)) exit TE
    if(worker_operation/=0)then
      call isolated_state_control(worker_operation,worker_status)
      call worker_state_reply(worker_operation,worker_status,itt,worker_updates)
      cycle TE
    endif
    if(itt>=nt.or.worker_outputs_closed)then
      call worker_state_reply(0,14,itt,worker_updates)
      cycle TE
    endif
  endif
  itt=itt+1
  if(worker_requested())then
    if(worker_E_requested()) then
      rt%E_ext(:,itt)=worker_A
      rt%E_tot(:,itt)=worker_A
    else
      rt%Ac_ext(:,itt)=worker_A
      rt%Ac_tot(:,itt)=worker_A
    end if
  end if
  call nvtxStartRange('main loop', itt)

  ! Initialization always loads spsi_in, also for an odd restart step Mit.
  if(mod(itt-Mit,2)==1)then
    call time_evolution_step(Mit,nt,itt,lg,mg,system,rt,info,stencil,xc_func &
     & ,srg,srg_scalar,pp,ppg,ppn,spsi_in,spsi_out,tpsi,rho,rho_jm,rho_s,V_local,Vbox,Vh,Vh_stock1,Vh_stock2,Vxc &
     & ,Vpsl,fg,energy,ewald,md,ofl,poisson,singlescale,unfold)
  else
    call time_evolution_step(Mit,nt,itt,lg,mg,system,rt,info,stencil,xc_func &
     & ,srg,srg_scalar,pp,ppg,ppn,spsi_out,spsi_in,tpsi,rho,rho_jm,rho_s,V_local,Vbox,Vh,Vh_stock1,Vh_stock2,Vxc &
     & ,Vpsl,fg,energy,ewald,md,ofl,poisson,singlescale,unfold)
  end if

  if(worker_E_requested())worker_updates=worker_updates+1
  if(worker_requested()) then
    if(worker_E_requested()) then
      ! subdip already integrated the occupied density at this endpoint.
      worker_electron_number=rt%rIe(itt)
      call worker_reply_E(itt,rt%dDp_e(:,itt),rt%Dp_e(:,itt)+rt%Dp_i(:,itt),worker_electron_number, &
        energy%E_tot,itt==1.or.itt==nt.or.mod(itt,out_rt_energy_step)==0)
    else
      if(worker_norm_requested())then
        call calc_dip(info%icomm_r,lg,mg,system,rho_s,worker_density_integral)
        worker_electron_number=worker_density_integral(4)*system%Hvol
      endif
      call worker_reply(itt,rt%curr(:,itt))
    end if
  end if

  is_checkpoint_iter = (checkpoint_interval >= 1) .and. (mod(itt,checkpoint_interval) == 0)
  ! timer_now is a local elapsed duration; absolute cross-node clock offsets cancel.
  is_shutdown_time   = (time_shutdown > 0d0) .and. (timer_now(LOG_TOTAL) > time_shutdown)
  call comm_bcast(is_shutdown_time,nproc_group_global)

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
    if (mod(itt-Mit,2)==1) then
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
call timer_end(LOG_RT_ITERATION)
call timer_disable_sub

#ifdef __FUJITSU
call fapp_stop('time_evol',1,0) ! performance profiling
#endif

close(030) ! laser


!--------------------------------- end of time-evolution


!------------ Writing part -----------


! Worker outputs are consumed at each endpoint. Do not analyze uncomputed tail.
if(worker_requested()) then
  if(comm_is_root(nproc_id_global)) close(ofl%fh_rt)
  call finalize_xc(xc_func)
  return
end if

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

#ifdef USE_HSE
call exx_unbind_system(system)
call rt_state_clear(independent_hse)
#endif
contains

subroutine isolated_state_control(operation,status)
 integer,intent(in)::operation
 integer,intent(out)::status
 status=0
 select case(operation)
 case(1)
  if(allocated(worker_saved))then
   status=14;return
  endif
  if(yn_out_dns_rt/='n'.or.yn_out_dns_ac_je/='n'.or.yn_out_micro_je/='n'.or. &
     yn_out_elf_rt/='n'.or.yn_out_estatic_rt/='n'.or.yn_out_current_decomposed/='n'.or.yn_out_rvf_rt/='n')then
   status=15;return
  endif
  call step_output_capture(worker_outputs,ofl,status)
  if(status/=0)return
  allocate(worker_saved,stat=status)
  if(status/=0)return
    worker_saved%spsi_in=spsi_in
    worker_saved%spsi_out=spsi_out
    worker_saved%tpsi=tpsi
    worker_saved%rho=rho
    worker_saved%rho_jm=rho_jm
    worker_saved%Vh=Vh
    worker_saved%Vh_stock1=Vh_stock1
    worker_saved%Vh_stock2=Vh_stock2
    worker_saved%Vbox=Vbox
    worker_saved%Vpsl=Vpsl
    worker_saved%rho_s=rho_s
    worker_saved%V_local=V_local
    worker_saved%Vxc=Vxc
    worker_saved%rt=rt
    worker_saved%poisson=poisson
    worker_saved%energy=energy
  worker_saved%xc=system%xc_payload
  worker_saved%cval=xc_func%cval
  worker_saved%fields(:,1)=system%vec_Ac
  worker_saved%fields(:,2)=system%vec_Ac_ext
  worker_saved%fields(:,3)=system%vec_E
  worker_saved%fields(:,4)=system%vec_E_ext
  worker_saved%step=itt
 case(2)
  if(.not.allocated(worker_saved).or.worker_outputs_closed)then
   status=14;return
  endif
  ! A partial close must suspend further writes until explicit reconnection.
  worker_outputs_closed=.true.
  call step_output_close(worker_outputs,status)
  if(status/=0)return
    spsi_in=worker_saved%spsi_in
    spsi_out=worker_saved%spsi_out
    tpsi=worker_saved%tpsi
    rho=worker_saved%rho
    rho_jm=worker_saved%rho_jm
    Vh=worker_saved%Vh
    Vh_stock1=worker_saved%Vh_stock1
    Vh_stock2=worker_saved%Vh_stock2
    Vbox=worker_saved%Vbox
    Vpsl=worker_saved%Vpsl
    rho_s=worker_saved%rho_s
    V_local=worker_saved%V_local
    Vxc=worker_saved%Vxc
    rt=worker_saved%rt
    poisson=worker_saved%poisson
    energy=worker_saved%energy
  system%xc_payload=worker_saved%xc
  xc_func%cval=worker_saved%cval
  system%vec_Ac=worker_saved%fields(:,1)
  system%vec_Ac_ext=worker_saved%fields(:,2)
  system%vec_E=worker_saved%fields(:,3)
  system%vec_E_ext=worker_saved%fields(:,4)
  itt=worker_saved%step
  spsi_in%update_zwf_overlap=.true.
  spsi_out%update_zwf_overlap=.true.
  tpsi%update_zwf_overlap=.true.
 case(3,4)
  if(.not.allocated(worker_saved))then
   status=14;return
  endif
  call step_output_reopen(worker_outputs,status,ofl)
  if(status/=0)return
  worker_outputs_closed=.false.
  if(operation==3)deallocate(worker_saved)
 case(5)
  ! Read-only monotonic update counter; never part of restored numerical state.
 case default
  status=11
 end select
end subroutine

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
