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

recursive subroutine isolated_state_control(operation,status)
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
  call step_output_reopen(worker_outputs,status,ofl)
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
  if(status/=0)then
   deallocate(worker_saved);return
  endif
  call isolated_state_control(2,status)

 case default
  status=11
 end select
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
 character(32)::magic
 worker_checkpoint_limit=max(64_8*(int(nt,8)+1_8),64_8*int(system%ngrid,8))
 if(allocated(spsi_in%zwf))worker_checkpoint_limit=max(worker_checkpoint_limit,size(spsi_in%zwf,kind=8))
 if(allocated(spsi_in%rwf))worker_checkpoint_limit=max(worker_checkpoint_limit,size(spsi_in%rwf,kind=8))
 if(allocated(rt%near_Ecoeff))worker_checkpoint_limit=max(worker_checkpoint_limit,size(rt%near_Ecoeff,kind=8))
 if(allocated(rt%near_grid_values)) &
  worker_checkpoint_limit=max(worker_checkpoint_limit,size(rt%near_grid_values,kind=8))
 expected=[1,nt,system%nspin,system%no,system%nk,lg%num]
 header=expected;saved_dt=dt;cell=system%primitive_a;magic='SALMONPY_ISOLATED_DIRECT_V7'
 if(reading)then
  open(newunit=u,file='isolated-state.bin',access='stream',form='unformatted',status='old', &
    action='read',iostat=ios)
 else
  open(newunit=u,file='isolated-state.bin',access='stream',form='unformatted',status='new', &
    action='write',iostat=ios)
 endif
 status=ios
 if(ios/=0)return
 if(reading)then
  read(u,iostat=ios)magic,header,saved_dt,cell
  if(ios==0)then
   if(magic/='SALMONPY_ISOLATED_DIRECT_V7'.or.any(header/=expected))ios=1
   if(saved_dt/=dt.or.any(cell/=system%primitive_a))ios=1
  endif
 else
  write(u,iostat=ios)magic,header,saved_dt,cell
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
