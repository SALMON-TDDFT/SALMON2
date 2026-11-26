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
use em_field,only:apply_near_magnetic,activate_near_magnetic
use salmonpy_step_output
use sendrecv_grid,only:dealloc_cache
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
use parallelization, only: nproc_group_global,nproc_size_global
use communication, only: comm_is_root, comm_sync_all, comm_bcast,comm_get_max
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
integer(8)::worker_checkpoint_limit
type(step_output_state)::worker_outputs
integer::worker_operation,worker_status,worker_updates
logical::worker_outputs_closed

worker_checkpoint_limit=0_8
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

call apply_near_magnetic(0,mg,lg,system,rt,info,ppg)

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
    if(.not.worker_next_A(worker_A,worker_operation,rt,itt+1)) exit TE
    if(worker_operation==12)then
      call calc_dip(info%icomm_r,lg,mg,system,rho_s,worker_density_integral)
      worker_electron_number=worker_density_integral(4)*system%Hvol
      call worker_reply_E(itt,rt%dDp_e(:,itt),rt%Dp_e(:,itt)+rt%Dp_i(:,itt),worker_electron_number, &
        energy%E_tot,itt==0.or.itt==1.or.itt==nt.or.mod(itt,out_rt_energy_step)==0)
      cycle TE
    endif
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
  if(comm_is_root(nproc_id_global))then
    close(ofl%fh_rt)
    if(worker_host_requested())close(ofl%fh_rt_energy)
  endif
  if(worker_host_requested())then
    call dealloc_cache(srg)
    call dealloc_cache(srg_scalar)
  endif
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

recursive subroutine isolated_state_control(operation,status)
 integer,intent(in)::operation
 integer,intent(out)::status
 status=0
 ! Read-only density/magnetic exports still have single-rank assumptions.
 if(worker_host_requested().and.(operation==8.or.operation==10.or.operation==11))then
  status=15;return
 endif
 select case(operation)
 case(1)
  if(allocated(worker_saved))then
   status=14;return
  endif
  if(yn_out_dns_rt/='n'.or.yn_out_dns_ac_je/='n'.or.yn_out_micro_je/='n'.or. &
     yn_out_elf_rt/='n'.or.yn_out_estatic_rt/='n'.or.yn_out_current_decomposed/='n'.or.yn_out_rvf_rt/='n')then
   status=15;return
  endif
  if(comm_is_root(nproc_id_global))call step_output_capture(worker_outputs,ofl,status)
  if(worker_host_requested())call comm_get_max(status,nproc_group_global)
  if(status/=0)return
  allocate(worker_saved,stat=status)
  if(worker_host_requested())call comm_get_max(status,nproc_group_global)
  if(status/=0)then
    if(allocated(worker_saved))deallocate(worker_saved)
    return
  endif
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
  if(comm_is_root(nproc_id_global))call step_output_close(worker_outputs,status)
  if(worker_host_requested())call comm_get_max(status,nproc_group_global)
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
    call activate_near_magnetic(system,rt,ppg)
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
  if(comm_is_root(nproc_id_global))call step_output_reopen(worker_outputs,status,ofl)
  if(worker_host_requested())call comm_get_max(status,nproc_group_global)
  if(status/=0)return
  worker_outputs_closed=.false.
  if(operation==3)deallocate(worker_saved)
 case(5)
  ! Read-only monotonic update counter; never part of restored numerical state.
 case(6)
  ! Saving copies numerical state only; no propagation or shared process handles.
  call isolated_state_control(1,status)
  if(status/=0)return
  call isolated_disk_io(.false.,status)
  deallocate(worker_saved)
 case(7)
  ! Only a fresh worker may import. Failed reads never install partial arrays.
  if(itt/=0.or.worker_updates/=0.or.allocated(worker_saved))then
   status=14;return
  endif
  call isolated_state_control(1,status)
  if(status/=0)return
  call isolated_disk_io(.true.,status)
  if(worker_host_requested())call comm_get_max(status,nproc_group_global)
  if(status/=0)then
   deallocate(worker_saved);return
  endif
  call isolated_state_control(2,status)

 case(8)
  call isolated_density_samples(status)
 case(9)
  call isolated_grid_stage(status)
 case(10,11)
  call isolated_orbital_response(operation==11,status)
 case default
  status=11
 end select
end subroutine

! Install both live E/H frames only after complete validation. No evolution.
subroutine isolated_grid_stage(status)
 use iso_c_binding,only:c_int
 use salmonpy_grid_field,only:grid_validate,grid_values
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 integer,intent(out)::status
 integer(c_int)::header(6)
 real(8)::geometry(7),origins(3,3,2),spacing(3),point(3),mapped(3),ev(3),hv(3)
 real(8),allocatable::frames(:,:,:,:,:,:)
 integer::u,ios,close_ios,allocation_status,k,n,valid,corner,d,index(3),atom,j,ilma
 integer(8)::bytes,position_in_file
 status=34
 if(.not.allocated(rt%near_grid_values))return
 if(worker_host_requested())then
   call isolated_grid_stage_host(status);return
 endif
 n=itt+1
 if(n<1.or.n>rt%max_steps)return
 status=35
 open(newunit=u,file='salmonpy-grid-frames.bin',access='stream',form='unformatted', &
      status='old',action='read',iostat=ios)
 if(ios/=0)return
 read(u,iostat=ios)header
 if(ios/=0)goto 900
 if(any(header/=[16909060_c_int,1_c_int,int(n,c_int),int(rt%near_grid_counts,c_int)]))goto 900
 read(u,iostat=ios)geometry,origins,spacing
 if(ios/=0)goto 900
 if(any(geometry(1:3)/=rt%near_geometry(1:3)).or.any(geometry(4:6)/=rt%near_grid_center))goto 900
 if(.not.ieee_is_finite(geometry(7)))goto 900
 if(abs(geometry(7)-rt%time_step)>1d-12*rt%time_step)goto 900
 if(any(origins/=rt%near_grid_origins).or.any(spacing/=rt%near_grid_spacing))goto 900
 allocate(frames(3,rt%near_grid_counts(1),rt%near_grid_counts(2),rt%near_grid_counts(3),2,2), &
          stat=allocation_status)
 if(allocation_status/=0)then
  status=32;goto 900
 endif
 read(u,iostat=ios)frames
 if(ios/=0)goto 900
 inquire(unit=u,size=bytes,pos=position_in_file,iostat=ios)
 if(ios/=0)goto 900
 if(bytes/=position_in_file-1_8)goto 900
 do k=1,2
  call grid_validate(rt%near_grid_counts,rt%near_grid_origins,rt%near_grid_spacing,frames(:,:,:,:,:,k),valid)
  if(valid/=0)goto 900
 enddo
 ! Fixed ions and geometry: check extended magnetic query locations once,
 ! before enabling H for an initially electric-only model. Rectangular coverage
 ! of endpoints also covers all straight integration paths between them.
 if(.not.rt%near_magnetic.and.any(frames(:,:,:,:,2,:)/=0d0))then
  do corner=0,7
   do d=1,3
    index(d)=merge(mg%ie_array(d),mg%is_array(d),btest(corner,d-1))
    point(d)=lg%coordinate(lg%is(d),d)+(index(d)-lg%is(d))*system%Hgs(d)
   enddo
   mapped=rt%near_grid_center+point-rt%near_geometry(1:3)
   call grid_values(rt%near_grid_counts,origins,spacing,frames(:,:,:,:,:,1),mapped,ev,hv,valid)
   if(valid/=0)goto 900
  enddo
  do ilma=1,ppg%nlma
   atom=ppg%ia_tbl(ilma)
   mapped=rt%near_grid_center+system%Rion(:,atom)-rt%near_geometry(1:3)
   call grid_values(rt%near_grid_counts,origins,spacing,frames(:,:,:,:,:,1),mapped,ev,hv,valid)
   if(valid/=0)goto 900
   do j=1,ppg%mps(atom)
    mapped=rt%near_grid_center+system%Rion(:,atom)+ppg%rxyz(:,j,atom)-rt%near_geometry(1:3)
    call grid_values(rt%near_grid_counts,origins,spacing,frames(:,:,:,:,:,1),mapped,ev,hv,valid)
    if(valid/=0)goto 900
   enddo
  enddo
 endif
 close(u,iostat=close_ios)
 if(close_ios/=0)return
 rt%near_grid_values(:,:,:,:,:,n:n+1)=frames
 if(any(frames(:,:,:,:,2,:)/=0d0))rt%near_magnetic=.true.
 status=0;return
900 continue
 close(u,iostat=close_ios)
end subroutine

! Host root reads once; only the owning molecular group receives the frames.
subroutine isolated_grid_stage_host(status)
 use iso_c_binding,only:c_int
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 use salmonpy_grid_field,only:grid_validate,grid_values
#ifdef USE_MPI
 use mpi,only:MPI_Bcast,MPI_DOUBLE_PRECISION,MPI_SUCCESS
#endif
 integer,intent(out)::status
 integer(c_int)::header(6)
 integer::n,u,ios,close_ios,astat,k,valid,corner,d,index(3),atom,j,ilma,ierr
 integer(8)::bytes,position_in_file
 real(8)::geometry(7),origins(3,3,2),spacing(3),point(3),mapped(3),ev(3),hv(3)
 real(8),allocatable::frames(:,:,:,:,:,:)
 status=34;n=itt+1
 if(.not.allocated(rt%near_grid_values))return
 if(n<1.or.n>rt%max_steps)return
 if(6_8*product(int(rt%near_grid_counts,8))*2_8>huge(1))return
 allocate(frames(3,rt%near_grid_counts(1),rt%near_grid_counts(2),rt%near_grid_counts(3),2,2),stat=astat)
 status=0
 if(astat/=0)status=32
 call comm_get_max(status,nproc_group_global)
 if(status/=0)return
 if(comm_is_root(nproc_id_global))then
   status=35
   open(newunit=u,file='salmonpy-grid-frames.bin',access='stream',form='unformatted', &
        status='old',action='read',iostat=ios)
   if(ios==0)then
     read(u,iostat=ios)header
     if(ios==0)then
       if(any(header/=[16909060_c_int,1_c_int,int(n,c_int),int(rt%near_grid_counts,c_int)]))ios=1
     endif
     if(ios==0)read(u,iostat=ios)geometry,origins,spacing
     if(ios==0)read(u,iostat=ios)frames
     if(ios==0)then
       inquire(unit=u,size=bytes,pos=position_in_file,iostat=ios)
       if(ios==0.and.bytes/=position_in_file-1_8)ios=1
     endif
     close(u,iostat=close_ios)
     if(ios==0.and.close_ios==0)status=0
   endif
 endif
 call comm_bcast(status,nproc_group_global)
 if(status/=0)return
 call comm_bcast(geometry,nproc_group_global)
 call comm_bcast(origins,nproc_group_global)
 call comm_bcast(spacing,nproc_group_global)
#ifdef USE_MPI
 call MPI_Bcast(frames,size(frames),MPI_DOUBLE_PRECISION,0,nproc_group_global,ierr)
 if(ierr/=MPI_SUCCESS)status=35
#else
 status=15
#endif
 if(any(geometry(1:3)/=rt%near_geometry(1:3)).or.any(geometry(4:6)/=rt%near_grid_center))status=35
 if(.not.ieee_is_finite(geometry(7)))status=35
 if(abs(geometry(7)-rt%time_step)>1d-12*rt%time_step)status=35
 if(any(origins/=rt%near_grid_origins).or.any(spacing/=rt%near_grid_spacing))status=35
 do k=1,2
   call grid_validate(rt%near_grid_counts,origins,spacing,frames(:,:,:,:,:,k),valid)
   if(valid/=0)status=35
 enddo
 if(status==0.and..not.rt%near_magnetic.and.any(frames(:,:,:,:,2,:)/=0d0))then
   do corner=0,7
     do d=1,3
       index(d)=merge(mg%ie_array(d),mg%is_array(d),btest(corner,d-1))
       point(d)=lg%coordinate(lg%is(d),d)+(index(d)-lg%is(d))*system%Hgs(d)
     enddo
     mapped=rt%near_grid_center+point-rt%near_geometry(1:3)
     call grid_values(rt%near_grid_counts,origins,spacing,frames(:,:,:,:,:,1),mapped,ev,hv,valid)
     if(valid/=0)status=35
   enddo
   do ilma=1,ppg%nlma
     atom=ppg%ia_tbl(ilma)
     mapped=rt%near_grid_center+system%Rion(:,atom)-rt%near_geometry(1:3)
     call grid_values(rt%near_grid_counts,origins,spacing,frames(:,:,:,:,:,1),mapped,ev,hv,valid)
     if(valid/=0)status=35
     do j=1,ppg%mps(atom)
       mapped=rt%near_grid_center+system%Rion(:,atom)+ppg%rxyz(:,j,atom)-rt%near_geometry(1:3)
       call grid_values(rt%near_grid_counts,origins,spacing,frames(:,:,:,:,:,1),mapped,ev,hv,valid)
       if(valid/=0)status=35
     enddo
   enddo
 endif
 call comm_get_max(status,nproc_group_global)
 if(status/=0)return
 rt%near_grid_values(:,:,:,:,:,n:n+1)=frames
 if(any(frames(:,:,:,:,2,:)/=0d0))rt%near_magnetic=.true.
end subroutine

! Read-only -d<E>/dH query. All scratch/model copies are local to this operation.
subroutine isolated_orbital_response(audit,status)
 use iso_c_binding,only:c_int
 use salmonpy_orbital_response,only:orbital_grid_magnetization
 use hamiltonian,only:hpsi
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 logical,intent(in)::audit
 integer,intent(out)::status
 type(s_orbital)::psi,applied
 type(s_dft_system)::trial_system
 type(s_rt)::trial_rt
 type(s_pp_grid)::trial_ppg
 real(8),allocatable::M(:,:,:,:),direction(:,:,:,:)
 real(8)::energies(2,2),epsilon,weight,expectation,derivative,errors(2),predicted
 integer::u,ios,close_ios,d,ix,iy,iz,im,ik,io,k,side,code
 integer(c_int)::header(6)
 status=34
 if(yn_spinorbit/='n'.or.yn_periodic/='n'.or.yn_md/='n')return
 if(.not.allocated(rt%near_grid_values))return
 ! Match the occupied orbital at the accepted endpoint, including restored parity.
 if(mod(itt-Mit,2)==1)then
  psi=spsi_out
 else
  psi=spsi_in
 endif
 call orbital_grid_magnetization(mg,lg,system,info,stencil,ppg,rt,psi,M,status)
 if(status/=0)return
 status=33
 if(.not.audit)then
  open(newunit=u,file='salmonpy-orbital-response.bin',access='stream',form='unformatted', &
       status='replace',action='write',iostat=ios)
  if(ios/=0)return
  header=[16909060_c_int,1_c_int,int(itt,c_int),int(rt%near_grid_counts,c_int)]
  write(u,iostat=ios)header
  if(ios==0)write(u,iostat=ios)rt%near_grid_origins,rt%near_grid_spacing
  if(ios==0)write(u,iostat=ios)M
  close(u,iostat=close_ios)
  if(ios/=0.or.close_ios/=0)return
 else
  ! Explicit audit only: finite differences of the actual Hamiltonian, frozen psi.
  allocate(direction(3,rt%near_grid_counts(1),rt%near_grid_counts(2),rt%near_grid_counts(3)),stat=code)
  status=32;if(code/=0)return
  do iz=1,rt%near_grid_counts(3)
  do iy=1,rt%near_grid_counts(2)
  do ix=1,rt%near_grid_counts(1)
  do d=1,3
   direction(d,ix,iy,iz)=sin(real(ix+2*iy+3*iz+d,8))
  enddo
  enddo
  enddo
  enddo
  predicted=-sum(M*direction)*product(rt%near_grid_spacing)
  do k=1,2
   epsilon=10d0**(-k-2)
   do side=1,2
    trial_system=system;trial_rt=rt;trial_ppg=ppg;applied=psi
    trial_rt%near_magnetic=.true.
    trial_rt%near_grid_values(:,:,:,:,2,rt%near_index)= &
      rt%near_grid_values(:,:,:,:,2,rt%near_index)+(2*side-3)*epsilon*direction
    call apply_near_magnetic(rt%near_index,mg,lg,trial_system,trial_rt,info,trial_ppg)
    call hpsi(psi,applied,info,mg,V_local,trial_system,stencil,srg,trial_ppg)
    expectation=0d0
    do im=info%im_s,info%im_e
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
     weight=system%rocc(io,ik,1)*system%wtk(ik)*system%Hvol
     do iz=mg%is(3),mg%ie(3)
     do iy=mg%is(2),mg%ie(2)
     do ix=mg%is(1),mg%ie(1)
      expectation=expectation+weight*real(conjg(psi%zwf(ix,iy,iz,1,io,ik,im))* &
                                          applied%zwf(ix,iy,iz,1,io,ik,im),8)
     enddo
     enddo
     enddo
    enddo
    enddo
    enddo
    energies(side,k)=expectation
   enddo
   derivative=(energies(2,k)-energies(1,k))/(2d0*epsilon)
   errors(k)=abs(derivative-predicted)
  enddo
  do k=1,2
   if(.not.ieee_is_finite(errors(k)))return
  enddo
  open(newunit=u,file='salmonpy-orbital-audit.bin',access='stream',form='unformatted', &
       status='replace',action='write',iostat=ios)
  if(ios/=0)return
  write(u,iostat=ios)predicted,errors,energies
  close(u,iostat=close_ios)
  if(ios/=0.or.close_ios/=0)return
 endif
 status=0
end subroutine

! Explicit read-only density transfer. A private cache is consumed by the
! Python owner; it is not part of physical state or a permanent output stream.
subroutine isolated_density_samples(status)
 use iso_c_binding,only:c_int
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 integer,intent(out)::status
 integer::u,ios,close_ios,a,ix,iy,iz,spin,n,allocation_status
 integer(8)::count
 real(8),allocatable::samples(:,:)
 integer(c_int)::header(4)
 status=32
 if(any(mg%is/=lg%is).or.any(mg%ie/=lg%ie))return
 count=product(int(mg%ie-mg%is+1,8))
 if(count<1.or.count>huge(1_c_int))return
 n=int(count)
 allocate(samples(5,n),stat=allocation_status)
 if(allocation_status/=0)return
 a=0
 do iz=mg%is(3),mg%ie(3)
 do iy=mg%is(2),mg%ie(2)
 do ix=mg%is(1),mg%ie(1)
  a=a+1
  samples(1:3,a)=[lg%coordinate(ix,1),lg%coordinate(iy,2),lg%coordinate(iz,3)]
  samples(4,a)=0d0
  do spin=1,size(rho_s)
   samples(4,a)=samples(4,a)-rho_s(spin)%f(ix,iy,iz)
  enddo
  samples(5,a)=system%Hvol
 enddo
 enddo
 enddo
 do a=1,n
 do spin=1,5
  if(.not.ieee_is_finite(samples(spin,a)))return
 enddo
 enddo
 if(system%Hvol<=0d0)return
 status=33
 open(newunit=u,file='salmonpy-charge-samples.bin',access='stream',form='unformatted', &
      status='replace',action='write',iostat=ios)
 if(ios/=0)return
 header=[16909060_c_int,1_c_int,int(itt,c_int),int(n,c_int)]
 write(u,iostat=ios)header
 if(ios==0)write(u,iostat=ios)samples
 close(u,iostat=close_ios)
 if(ios/=0.or.close_ios/=0)return
 status=0
end subroutine

integer(8) function checkpoint_io_limit()
 checkpoint_io_limit=worker_checkpoint_limit
end function

subroutine isolated_disk_io(reading,status)
 use salmonpy_grid_field,only:grid_validate
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 logical,intent(in)::reading
 integer,intent(out)::status
 integer::u,ios,close_status,header(8),expected(8),extra,grid_step,grid_status,k
 real(8)::saved_dt,cell(3,3)
 character(32)::magic,expected_magic
 character(80)::filename
 integer::layout(7),saved_layout(7)
 worker_checkpoint_limit=max(64_8*(int(nt,8)+1_8),64_8*int(system%ngrid,8))
 if(allocated(spsi_in%zwf))worker_checkpoint_limit=max(worker_checkpoint_limit,size(spsi_in%zwf,kind=8))
 if(allocated(spsi_in%rwf))worker_checkpoint_limit=max(worker_checkpoint_limit,size(spsi_in%rwf,kind=8))
 if(allocated(rt%near_Ecoeff))worker_checkpoint_limit=max(worker_checkpoint_limit,size(rt%near_Ecoeff,kind=8))
 if(allocated(rt%near_grid_values)) &
  worker_checkpoint_limit=max(worker_checkpoint_limit,size(rt%near_grid_values,kind=8))
 expected=[1,nt,system%nspin,system%no,system%nk,lg%num]
 header=expected;saved_dt=dt;cell=system%primitive_a;magic='SALMONPY_ISOLATED_DIRECT_V7'
 filename='isolated-state.bin'
 layout=[nproc_size_global,nproc_id_global,nproc_k,nproc_ob,nproc_rgrid]
 if(worker_host_requested())then
   magic='SALMONPY_ISOLATED_MPI_V1'
   write(filename,'(a,i8.8,a)')'isolated-state-',nproc_id_global,'.bin'
 endif
 expected_magic=magic
 if(reading)then
  open(newunit=u,file=trim(filename),access='stream',form='unformatted',status='old', &
    action='read',iostat=ios)
 else
  open(newunit=u,file=trim(filename),access='stream',form='unformatted',status='new', &
    action='write',iostat=ios)
 endif
 status=ios
 if(ios/=0)return
 if(reading)then
  read(u,iostat=ios)magic,header,saved_dt,cell
  if(ios==0)then
   if(magic/=expected_magic.or.any(header/=expected))ios=1
   if(saved_dt/=dt.or.any(cell/=system%primitive_a))ios=1
  endif
 else
  write(u,iostat=ios)magic,header,saved_dt,cell
 endif
 if(ios==0.and.worker_host_requested())then
   if(reading)then
     read(u,iostat=ios)saved_layout
     if(ios==0.and.any(saved_layout/=layout))ios=1
   else
     write(u,iostat=ios)layout
   endif
 endif
 if(ios==0)call isolated_payload_io(u,reading,worker_saved,ios)
 if(ios==0.and.reading)then
  if(allocated(worker_saved%rt%near_grid_values))then
   do k=1,3
    if(.not.ieee_is_finite(worker_saved%rt%near_grid_center(k)).or. &
       .not.ieee_is_finite(worker_saved%rt%near_geometry(k)))ios=1
   enddo
   do grid_step=0,nt+1
    if(ios/=0)exit
    call grid_validate(worker_saved%rt%near_grid_counts,worker_saved%rt%near_grid_origins, &
      worker_saved%rt%near_grid_spacing,worker_saved%rt%near_grid_values(:,:,:,:,:,grid_step),grid_status)
    if(grid_status/=0)ios=1
   enddo
  endif
 endif
 if(ios==0)then
  magic='SALMONPY_ISOLATED_COMPLETE_V7'
  if(reading)then
   read(u,iostat=ios)magic
   if(ios==0.and.magic/='SALMONPY_ISOLATED_COMPLETE_V7')ios=1
   if(ios==0)then
    read(u,iostat=ios)extra
    if(is_iostat_end(ios))then
     ios=0
    else
     ios=1
    endif
   endif
   if(worker_saved%step<0.or.worker_saved%step>nt)ios=1
   if(.not.ieee_is_finite(worker_saved%cval))ios=1
  else
   write(u,iostat=ios)magic
  endif
 endif
 close(u,iostat=close_status)
 if(close_status/=0)ios=close_status
 status=ios
end subroutine

subroutine isolated_payload_io(unit,reading,value,ios)
  integer,intent(in)::unit
  logical,intent(in)::reading
  type(isolated_savepoint),intent(inout)::value
  integer,intent(inout)::ios
  logical::present_array
  integer::lo(7),hi(7),i1,i2,k,allocation_status
  integer(8)::elements,extent
  if(ios/=0)return
  call checkpoint_io_s_orbital(unit,reading,value%spsi_in,ios)
  if(ios/=0)return
  call checkpoint_io_s_orbital(unit,reading,value%spsi_out,ios)
  if(ios/=0)return
  call checkpoint_io_s_orbital(unit,reading,value%tpsi,ios)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%rho,ios)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%rho_jm,ios)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%Vh,ios)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%Vh_stock1,ios)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%Vh_stock2,ios)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%Vbox,ios)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%Vpsl,ios)
  if(ios/=0)return
  call checkpoint_io_s_rt(unit,reading,value%rt,ios)
  if(ios/=0)return
  call checkpoint_io_s_poisson(unit,reading,value%poisson,ios)
  if(ios/=0)return
  call checkpoint_io_s_dft_energy(unit,reading,value%energy,ios)
  if(ios/=0)return
  present_array=allocated(value%rho_s)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%rho_s))deallocate(value%rho_s)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%rho_s))then
        if(any(lo(1:1)/=lbound(value%rho_s)).or.any(hi(1:1)/=ubound(value%rho_s)))then
          ios=1;return
        endif
      else
        allocate(value%rho_s(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%rho_s);hi(1:1)=ubound(value%rho_s)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    do i1=lbound(value%rho_s,1),ubound(value%rho_s,1)
    call checkpoint_io_s_scalar(unit,reading,value%rho_s(i1),ios)
    if(ios/=0)return
    enddo
  endif
  present_array=allocated(value%V_local)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%V_local))deallocate(value%V_local)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%V_local))then
        if(any(lo(1:1)/=lbound(value%V_local)).or.any(hi(1:1)/=ubound(value%V_local)))then
          ios=1;return
        endif
      else
        allocate(value%V_local(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%V_local);hi(1:1)=ubound(value%V_local)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    do i1=lbound(value%V_local,1),ubound(value%V_local,1)
    call checkpoint_io_s_scalar(unit,reading,value%V_local(i1),ios)
    if(ios/=0)return
    enddo
  endif
  present_array=allocated(value%Vxc)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%Vxc))deallocate(value%Vxc)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%Vxc))then
        if(any(lo(1:1)/=lbound(value%Vxc)).or.any(hi(1:1)/=ubound(value%Vxc)))then
          ios=1;return
        endif
      else
        allocate(value%Vxc(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%Vxc);hi(1:1)=ubound(value%Vxc)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    do i1=lbound(value%Vxc,1),ubound(value%Vxc,1)
    call checkpoint_io_s_scalar(unit,reading,value%Vxc(i1),ios)
    if(ios/=0)return
    enddo
  endif
  if(reading)then
    read(unit,iostat=ios) value%fields
  else
    write(unit,iostat=ios) value%fields
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios)value%cval,value%step,value%xc%use_tau_operator, &
      value%xc%vtau_has_shadow_values,value%xc%e_tau
  else
    write(unit,iostat=ios)value%cval,value%step,value%xc%use_tau_operator, &
      value%xc%vtau_has_shadow_values,value%xc%e_tau
  endif
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%xc%vtau,ios)
  call checkpoint_io_s_vector(unit,reading,value%xc%uvel,ios)
end subroutine

#include "checkpoint_arrays.inc"

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
