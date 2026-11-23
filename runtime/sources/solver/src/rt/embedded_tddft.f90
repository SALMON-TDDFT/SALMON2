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

! SALMONpy one-point embedded TDDFT. Native physics uses existing RT kernels.
#include "config.h"
#ifdef SALMONPY_FAST_CORE
module salmonpy_standard_tddft_context
  use filesystem,only:resolve_path
#else
module salmonpy_tddft_context
  use filesystem,only:resolve_path
#endif
use tddft_worker
use noncollinear_module,only:s_noncollinear_state,swap_noncollinear_state,snapshot_noncollinear_state
use salmonpy_gs_seed, only: s_gs_seed
#ifdef SALMONPY_FAST_CORE
use salmonpy_gs_context, only: gs_seed_from_handle
#else
use salmonpy_embedded_gs, only: borrow_gs
use salmonpy_gs_storage, only: s_gs_storage
use salmonpy_gs_model, only: s_gs_model
#endif
use density_matrix, only: calc_current
use hamiltonian, only: update_kvector_nonlocalpt
use math_constants, only: pi
use salmon_global
use structures
use inputoutput, only: s_input_parser
use salmonpy_model_inputs, only: s_model_inputs
#ifndef SALMONPY_FAST_CORE
use code_optimization, only: s_model_optimization,park_model_optimization,activate_model_optimization
#endif
use salmonpy_step_output
use parallelization, only: adjust_elapse_time, nproc_group_global
use communication, only: comm_is_root, comm_sync_all, comm_bcast,comm_get_globalinfo, &
  comm_begin_children,comm_end_children,comm_release_children,comm_child_counts
use salmon_xc, only: finalize_xc
#ifndef SALMONPY_FAST_CORE
use salmon_xc, only: finalize_xc,s_model_xc_workspace,park_model_xc_workspace,activate_model_xc_workspace
#endif
use timer
use write_sub, only: write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d
use initialization_rt_sub
use checkpoint_restart_sub
use jellium, only: check_condition_jm
use nvtx_wrapper
use parallelization, only: nproc_id_global
use iso_c_binding
use sendrecv_grid, only: dealloc_cache
implicit none
private
public :: sp_tddft_init,sp_tddft_step,sp_tddft_close,sp_tddft_checkpoint,sp_tddft_step_performance
! Only numerical state: communicators, output units and library handles stay in the fresh context.
type :: checkpoint_payload
  type(s_noncollinear_state) :: noncollinear
  type(s_orbital) :: spsi_in,spsi_out,tpsi
  type(s_scalar) :: rho,rho_jm,Vh,Vh_stock1,Vh_stock2,Vbox,Vpsl
  type(s_scalar),allocatable :: rho_s(:),V_local(:),Vxc(:)
  type(s_rt) :: rt
  type(s_poisson) :: poisson
  type(s_dft_energy) :: energy
  real(8) :: fields(3,5),worker_A(3),worker_current(3,2),accepted_A(3)
  integer :: Mit,itt
end type
! A model owns immutable initialization values, distinct from RT numerical owners.
! Registry tokens own models; state contexts borrow their immutable model owner.
type :: rt_model
  type(s_input_parser) :: parser
  type(s_gs_seed),allocatable :: seed
  type(s_model_inputs) :: runtime_inputs
#ifndef SALMONPY_FAST_CORE
  type(s_model_optimization) :: optimization
  type(s_model_xc_workspace) :: xc_workspace
#endif
  logical :: runtime_ready=.false.,closed=.false.
  integer :: references=0
  integer(c_int64_t) :: token=0_c_int64_t
  real(8) :: initial_A(3)=0d0
end type
type(s_gs_seed),pointer,save :: pending_gs=>null()
logical,save :: gs_timers_suspended=.false.,gs_saved_worker=.false.
type :: rt_context
  type(rt_model),pointer :: model=>null()
type(checkpoint_payload),allocatable :: step_saved,pending
real(8) :: step_saved_cval=0d0,pending_xc_cval=0d0
logical :: step_restored=.false.
type(s_xc_operator_payload),allocatable :: step_saved_xc
type(step_output_state) :: step_outputs
real(8) :: accepted_A(3)=0d0
integer(8) :: checkpoint_array_limit=0_8
integer :: mpi_child_owner=0
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

real(8) :: worker_A(3),worker_current(3,2)
integer :: Mit, itt
logical :: is_checkpoint_iter, is_shutdown_time, is_checkpoint

end type
type(rt_context),pointer,save :: state=>null()
type(rt_model),pointer,save :: model=>null(),active_model=>null()
type model_slot
  type(rt_model),pointer :: owner=>null()
  integer(c_int64_t) :: generation=0_c_int64_t
end type
type(model_slot),allocatable,save :: models(:)
logical,save :: registry_mode=.false.,registering=.false.,model_spawning=.false.
logical,save :: used=.false.,handle_mode=.false.,model_timers_alive=.false.,handle_initializing=.false.
integer,save :: handle_group=0
integer(c_int64_t),parameter :: handle_stride=2147483648_c_int64_t
type(s_timer_context),save :: legacy_timers,legacy_host_timers
logical,save :: legacy_timers_alive=.false.,legacy_host_alive=.false.,handle_timer_active=.false.
type(s_noncollinear_state),save :: handle_host_noncollinear
type(s_timer_context),save :: handle_host_timers
logical,save :: handle_host_timer_alive=.false.
type handle_slot
  type(s_noncollinear_state) :: noncollinear
  type(s_timer_context) :: timers
  logical :: timers_alive=.false.
  type(rt_context),pointer :: context=>null()
  integer(c_int64_t) :: generation=0_c_int64_t
  character(len(sysname)) :: label=''
end type
type(handle_slot),allocatable,save :: handles(:)

contains

integer function live_contexts() result(n)
  integer :: i
  logical :: selected
  n=live_handles();selected=.false.
  if(.not.associated(state))return
  if(allocated(handles))then
    do i=1,size(handles)
      if(associated(handles(i)%context,state))selected=.true.
    enddo
  endif
  if(.not.selected)n=n+1
end function
subroutine park_active_model()
  if(.not.associated(active_model))return
#ifndef SALMONPY_FAST_CORE
  call park_model_optimization(active_model%optimization)
  call park_model_xc_workspace(active_model%xc_workspace)
#endif
  nullify(active_model)
end subroutine
subroutine select_model(owner,initializing,status,checked)
  use inputoutput,only:initialize_inputoutput_units
  use sym_sub,only:read_sw_symmetry
  type(rt_model),pointer,intent(in) :: owner
  logical,intent(in) :: initializing
  integer,intent(out) :: status
  logical,optional,intent(in)::checked
  logical::verify
  verify=.true.;if(present(checked))verify=checked
  status=12
  if(.not.associated(owner))return
  if(.not.initializing.and..not.owner%runtime_ready)return
  if(.not.initializing.and.associated(active_model,owner))then
    status=0;return
  endif
  call park_active_model()
  if(initializing)then
    call owner%parser%values%activate(status)
  else
    call owner%runtime_inputs%activate(status)
  endif
  if(verify)call handle_agree(status)
  if(status/=0)return
  call initialize_inputoutput_units
  calc_mode='RT';worker_embedded=.true.;worker_initial_A=owner%initial_A
  embedded_dt=dt;embedded_steps=nt
  embedded_decomposition=[nproc_k,nproc_ob,nproc_rgrid]
  call read_sw_symmetry(yn_symmetry)
#ifndef SALMONPY_FAST_CORE
  call activate_model_optimization(owner%optimization)
  call activate_model_xc_workspace(owner%xc_workspace)
#endif
  active_model=>owner
end subroutine
subroutine release_model(owner)
  type(rt_model),pointer,intent(inout) :: owner
  if(.not.associated(owner))return
  if(.not.owner%closed.or.owner%references/=0)return
  if(associated(active_model,owner))call park_active_model()
  if(associated(model,owner))nullify(model)
  deallocate(owner);nullify(owner)
end subroutine
integer(c_int) function sp_tddft_memory_gs_protocol() bind(C)
  sp_tddft_memory_gs_protocol=1
end function
integer(c_int) function sp_tddft_park_for_gs() bind(C)
  integer :: status
  status=0
  if(associated(state).or.(used.and..not.registry_mode))status=12
  call handle_agree(status)
  sp_tddft_park_for_gs=status
  if(status/=0)return
  call park_active_model()
  gs_saved_worker=worker_embedded;worker_embedded=.false.
  if(model_timers_alive)then
    call timer_end(LOG_RT_ITERATION)
    call timer_disable_sub
    call timer_end(LOG_TOTAL)
    model_timers_alive=.false.;gs_timers_suspended=.true.
  endif
end function
subroutine sp_tddft_resume_after_gs() bind(C)
  worker_embedded=gs_saved_worker
  if(gs_timers_suspended.and.live_contexts()>0)then
    call timer_begin(LOG_TOTAL)
    call timer_enable_sub
    call timer_begin(LOG_RT_ITERATION)
    model_timers_alive=.true.
  endif
  gs_timers_suspended=.false.
end subroutine
integer(c_int) function sp_tddft_model_from_gs(gs_token,timestep,capacity,decomposition,model_token,state_token,J) bind(C)
  integer(c_int64_t),value :: gs_token
  real(c_double),value :: timestep
  integer(c_int),value :: capacity
  integer(c_int),intent(in) :: decomposition(5)
  integer(c_int64_t),intent(out) :: model_token,state_token
  real(c_double),intent(out) :: J(3)
#ifndef SALMONPY_FAST_CORE
  type(s_gs_storage),pointer :: gs
  type(s_gs_model),pointer :: gs_model
#endif
  type(s_gs_seed),target :: seed
  integer :: status
  model_token=0;state_token=0;J=0d0
#ifdef SALMONPY_FAST_CORE
  call gs_seed_from_handle(gs_token,seed,status,require_zero_field=.true.)
  if(status==0)then
    if(any(decomposition/=seed%decomposition))status=11
    if(associated(pending_gs))status=12
  endif
  call handle_agree(status)
  sp_tddft_model_from_gs=status
  if(status/=0)return
#else
  call borrow_gs(gs_token,gs,gs_model,status)
  if(status==0)then
    seed%decomposition=[gs%info%npk,gs%info%nporbital,gs%info%nprgrid]
    if(any(decomposition/=seed%decomposition))status=11
    if(associated(pending_gs))status=12
  endif
  call handle_agree(status)
  sp_tddft_model_from_gs=status
  if(status/=0)return
  call gs_model%copy_inputs(seed%inputs)
  seed%psi=gs%spsi
  seed%occupations=gs%system%rocc
  seed%grid=gs%lg%num;seed%core_start=gs%mg%is;seed%core_end=gs%mg%ie
  seed%orbital_start=gs%info%io_s;seed%orbital_end=gs%info%io_e
  seed%k_start=gs%info%ik_s;seed%k_end=gs%info%ik_e
#endif
  pending_gs=>seed
  sp_tddft_model_from_gs=sp_tddft_model_register('memory-gs'//c_null_char,timestep,capacity, &
      [0d0,0d0,0d0],decomposition,model_token,state_token,J)
  nullify(pending_gs)
end function

integer(c_int) function sp_tddft_registry_protocol() bind(C)
  sp_tddft_registry_protocol=1
end function
subroutine model_reserve(index,token,status)
  integer,intent(out)::index,status
  integer(c_int64_t),intent(out)::token
  type(model_slot),allocatable :: grown(:)
  integer::old,i,allocation_status,group,rank,n,parts(2),reference(2)
  status=0;index=0;old=0;token=0
  if(allocated(models))then
    old=size(models)
    do i=1,old
      if(associated(models(i)%owner))cycle
      if(models(i)%generation>=int(huge(0),c_int64_t))cycle
      index=i;exit
    enddo
  endif
  if(index==0)then
    if(old>huge(0)/2)status=16
    call handle_agree(status)
    if(status/=0)return
    allocate(grown(max(4,2*old)),stat=allocation_status)
    if(allocation_status/=0)status=16
    call handle_agree(status)
    if(status/=0)return
    if(old>0)grown(:old)=models
    call move_alloc(grown,models);index=old+1
  endif
  parts=[index,int(models(index)%generation+1_c_int64_t)];reference=parts
  call comm_get_globalinfo(group,rank,n);call comm_bcast(reference,group)
  if(any(reference/=parts))status=11
  call handle_agree(status)
  if(status/=0)return
  models(index)%generation=models(index)%generation+1_c_int64_t
  token=-(models(index)%generation*handle_stride+int(index,c_int64_t))
end subroutine
subroutine model_lookup(token,index,status)
  integer(c_int64_t),intent(in)::token
  integer,intent(out)::index,status
  integer::parts(2),reference(2),group,rank,n
  index=0;status=12
  if(.not.handle_main())return
  call comm_get_globalinfo(group,rank,n)
  status=0
  if(.not.registry_mode.or.group/=handle_group.or.token>=0_c_int64_t)status=12
  if(token<-(int(huge(0),c_int64_t)*handle_stride+int(huge(0),c_int64_t)))status=12
  call handle_agree(status)
  if(status/=0)return
  parts=[int(mod(-token,handle_stride)),int((-token)/handle_stride)];reference=parts
  call comm_bcast(reference,group)
  if(any(reference/=parts))status=11
  call handle_agree(status)
  if(status/=0)return
  index=parts(1)
  if(.not.allocated(models))then
    status=12
  else if(index<1.or.index>size(models))then
    status=12
  else if(models(index)%generation/=int(parts(2),c_int64_t))then
    status=12
  else if(.not.associated(models(index)%owner))then
    status=12
  endif
  call handle_agree(status)
end subroutine
integer(c_int) function sp_tddft_model_register(path,timestep,capacity,A,decomposition,model_token,state_token,J) bind(C)
  character(c_char),intent(in)::path(*)
  real(c_double),value::timestep
  integer(c_int),value::capacity
  real(c_double),intent(in)::A(3)
  integer(c_int),intent(in)::decomposition(5)
  integer(c_int64_t),intent(out)::model_token,state_token
  real(c_double),intent(out)::J(3)
  type(rt_model),pointer :: previous,failed,old_default
  integer::status,index,group,rank,n,restore_status
  integer(c_int64_t)::candidate
  model_token=0;state_token=0;J=0d0;sp_tddft_model_register=12
  if(.not.handle_main())return
  call comm_get_globalinfo(group,rank,n)
  status=0
  if(associated(state).or.(used.and..not.registry_mode))status=12
  if(registry_mode.and.group/=handle_group)status=12
  call handle_agree(status)
  if(status/=0)return
  call model_reserve(index,candidate,status)
  if(status/=0)then
    sp_tddft_model_register=status;return
  endif
  previous=>active_model;old_default=>model
  call park_active_model()
  registering=.true.
  status=sp_tddft_handle_init(path,timestep,capacity,A,decomposition,state_token,J)
  registering=.false.
  if(status==0)then
    registry_mode=.true.;model%token=candidate
    models(index)%owner=>model;model_token=candidate
  else
    if(associated(model))then
      if(.not.associated(model,old_default))then
        failed=>model;failed%closed=.true.;call release_model(failed)
      endif
    endif
    model=>old_default
    if(associated(previous))call select_model(previous,.false.,restore_status)
  endif
  sp_tddft_model_register=status
end function
integer(c_int) function sp_tddft_model_spawn(token,state_token,J) bind(C)
  integer(c_int64_t),value::token
  integer(c_int64_t),intent(out)::state_token
  real(c_double),intent(out)::J(3)
  integer::index,status
  state_token=0;J=0d0
  call model_lookup(token,index,status);sp_tddft_model_spawn=status
  if(status/=0)return
  model=>models(index)%owner
  model_spawning=.true.
  sp_tddft_model_spawn=sp_tddft_model_new(state_token,J)
  model_spawning=.false.
end function
integer(c_int) function sp_tddft_model_close(token) bind(C)
  integer(c_int64_t),value::token
  integer::index,status
  type(rt_model),pointer::owner
  call model_lookup(token,index,status);sp_tddft_model_close=status
  if(status/=0)return
  owner=>models(index)%owner;nullify(models(index)%owner)
  owner%closed=.true.;call release_model(owner)
  sp_tddft_model_close=0
end function
! Same immutable full-k semilocal model, common dt/capacity/decomposition.
! Pointer selection never copies numerical/MPI/Libxc owners. Calls are synchronous.
integer(c_int) function sp_tddft_handle_protocol() bind(C)
  sp_tddft_handle_protocol=1
end function
integer(c_int) function sp_tddft_model_protocol() bind(C)
  sp_tddft_model_protocol=1
end function
logical function handle_main()
  use communication,only:comm_is_attached
  use mpi,only:MPI_Is_thread_main,MPI_SUCCESS
  integer::ierr
  logical::main
  handle_main=.false.
  if(.not.comm_is_attached())return
  call MPI_Is_thread_main(main,ierr)
  handle_main=ierr==MPI_SUCCESS.and.main
end function
subroutine handle_agree(status)
  use communication,only:comm_get_max,comm_get_globalinfo
  integer,intent(inout)::status
  integer::group,rank,n
  call comm_get_globalinfo(group,rank,n)
  call comm_get_max(status,group)
end subroutine
integer function live_handles() result(n)
  integer::i
  n=0
  if(.not.allocated(handles))return
  do i=1,size(handles)
    if(associated(handles(i)%context))n=n+1
  enddo
end function
subroutine handle_reserve(index,token,status)
  use communication,only:comm_get_globalinfo
  integer,intent(out)::index,status
  integer(c_int64_t),intent(out)::token
  type(handle_slot),allocatable::grown(:)
  type(s_noncollinear_state) :: empty_noncollinear
  type(s_timer_context) :: empty_timers
  integer::i,old,n,allocation_status,group,rank,parts(2),reference(2)
  index=0;status=0;token=0;old=0
  if(allocated(handles))then
    old=size(handles)
    do i=1,old
      if(associated(handles(i)%context))cycle
      if(handles(i)%generation>=int(huge(0),c_int64_t))cycle
      index=i;exit
    enddo
  endif
  if(index==0)then
    if(old>huge(0)/2)status=16
    call handle_agree(status)
    if(status/=0)return
    allocate(grown(max(4,2*old)),stat=allocation_status)
    if(allocation_status/=0)status=16
    call handle_agree(status)
    if(status/=0)return
    if(old>0)grown(:old)=handles
    call move_alloc(grown,handles);index=old+1
  endif
  parts=[index,int(handles(index)%generation+1_c_int64_t)];reference=parts
  call comm_get_globalinfo(group,rank,n);call comm_bcast(reference,group)
  if(any(reference/=parts))status=11
  call handle_agree(status)
  if(status/=0)return
  handles(index)%generation=handles(index)%generation+1_c_int64_t
  token=handles(index)%generation*handle_stride+int(index,c_int64_t)
  handles(index)%noncollinear=empty_noncollinear
  handles(index)%timers=empty_timers
  handles(index)%timers_alive=.false.
end subroutine
subroutine handle_lookup(token,index,status,checked)
  use communication,only:comm_get_globalinfo
  integer(c_int64_t),intent(in)::token
  integer,intent(out)::index,status
  logical,optional,intent(in)::checked
  logical::verify
  integer::group,rank,n,parts(2),reference(2)
  verify=.true.;if(present(checked))verify=checked
  index=0;status=12
  if(.not.handle_main())return
  call comm_get_globalinfo(group,rank,n)
  status=0
  if(.not.handle_mode.or.group/=handle_group.or.token<=0_c_int64_t)status=12
  if(token/handle_stride>int(huge(0),c_int64_t))status=12
  if(verify)call handle_agree(status)
  if(status/=0)return
  parts=[int(mod(token,handle_stride)),int(token/handle_stride)]
  if(verify)then
    reference=parts;call comm_bcast(reference,group)
    if(any(reference/=parts))status=11
  endif
  if(verify)call handle_agree(status)
  if(status/=0)return
  index=parts(1)
  if(.not.allocated(handles))then
    status=12
  else if(index<1.or.index>size(handles))then
    status=12
  else if(handles(index)%generation/=int(parts(2),c_int64_t))then
    status=12
  else if(.not.associated(handles(index)%context))then
    status=12
  endif
  if(verify)call handle_agree(status)
end subroutine
subroutine handle_select(token,index,status,checked)
  integer(c_int64_t),intent(in)::token
  integer,intent(out)::index,status
  logical,optional,intent(in)::checked
  call handle_lookup(token,index,status,checked)
  if(status==0)call select_model(handles(index)%context%model,.false.,status,checked)
  if(status==0)then
    state=>handles(index)%context
    call handle_timer_enter(index)
  endif
end subroutine
subroutine legacy_timer_enter()
  call timer_context_swap(legacy_host_timers)
  call timer_context_swap(legacy_timers)
  call timer_initialize
  legacy_host_alive=model_timers_alive
  model_timers_alive=legacy_timers_alive
end subroutine
subroutine legacy_timer_leave()
  legacy_timers_alive=model_timers_alive
  model_timers_alive=legacy_host_alive
  call timer_context_swap(legacy_timers)
  call timer_context_swap(legacy_host_timers)
end subroutine
integer(c_int) function sp_tddft_legacy_timer_protocol() bind(C)
  sp_tddft_legacy_timer_protocol=1
end function
integer(c_int) function sp_tddft_legacy_timers(values) bind(C)
  real(c_double),intent(out)::values(4)
  values=0d0;sp_tddft_legacy_timers=12
  if(handle_mode.or..not.associated(state))return
  values=[timer_context_get(legacy_timers,LOG_TOTAL),timer_context_get(legacy_timers,LOG_INIT_RT), &
          timer_context_get(legacy_timers,LOG_CALC_RHO),timer_context_get(legacy_timers,LOG_CALC_HARTREE)]
  sp_tddft_legacy_timers=0
end function
subroutine handle_timer_enter(index)
  integer,intent(in) :: index
  call swap_noncollinear_state(handle_host_noncollinear)
  call swap_noncollinear_state(handles(index)%noncollinear)
  call timer_context_swap(handle_host_timers)
  call timer_context_swap(handles(index)%timers)
  call timer_initialize
  handle_timer_active=.true.
  handle_host_timer_alive=model_timers_alive
  model_timers_alive=handles(index)%timers_alive
end subroutine
subroutine handle_timer_leave(index)
  integer,intent(in) :: index
  call swap_noncollinear_state(handles(index)%noncollinear)
  call swap_noncollinear_state(handle_host_noncollinear)
  handles(index)%timers_alive=model_timers_alive
  model_timers_alive=handle_host_timer_alive
  call timer_context_swap(handles(index)%timers)
  call timer_context_swap(handle_host_timers)
  handle_timer_active=.false.
end subroutine
subroutine handle_unselect(index)
  integer,intent(in)::index
  handles(index)%context=>state
  call handle_timer_leave(index)
  nullify(state)
end subroutine
integer(c_int) function sp_tddft_handle_init(path,timestep,capacity,A,decomposition,token,J) bind(C)
  use communication,only:comm_get_globalinfo
  use sym_sub,only:use_symmetry
  use pseudo_pt_plusU_sub,only:PLUS_U_ON
  character(c_char),intent(in)::path(*)
  real(c_double),value::timestep
  integer(c_int),value::capacity
  real(c_double),intent(in)::A(3)
  integer(c_int),intent(in)::decomposition(5)
  integer(c_int64_t),intent(out)::token
  real(c_double),intent(out)::J(3)
  integer::index,status,group,rank,n,ignored
  integer(c_int64_t)::candidate
  token=0;candidate=0;J=0d0;sp_tddft_handle_init=12
  if(.not.handle_main())return
  status=0
  if(associated(state))status=12
  if(.not.registering.and.(used.or.handle_mode))status=12
  call handle_agree(status)
  if(status/=0)return
  call handle_reserve(index,candidate,status)
  if(status==0)then
    call handle_timer_enter(index)
    handle_initializing=.true.
    status=sp_tddft_init(path,timestep,capacity,A,decomposition,J)
    handle_initializing=.false.
  endif
  if(status==0)then
    select case(trim(xc))
    case('pz','tbmbj','libxc_pz','libxc_pbe')
    case default
      status=15
    end select
    if(use_symmetry.or..not.worker_spin_supported().or.yn_md/='n')status=15
    if(iperiodic/=3.or.propagator/='middlepoint'.or.trans_longi/='tr')status=15
    if(xname/='none'.or.cname/='none'.or.yn_jm/='n'.or.PLUS_U_ON)status=15
    if(projection_option/='no'.or.dm_unfold_option/='no'.or.state%singlescale%flag_use)status=15
    if(yn_out_rvf_rt/='n'.or.gram_schmidt_interval>=1.or.state%Mit/=0)status=15
    if(yn_out_dns_rt/='n'.or.yn_out_dns_ac_je/='n'.or.yn_out_micro_je/='n')status=15
    if(yn_out_elf_rt/='n'.or.yn_out_estatic_rt/='n'.or.yn_out_current_decomposed/='n')status=15
    if(yn_out_intraband_current/='n'.or.yn_out_mag_decomposed_rt/='n'.or.yn_out_mag_micro_rt/='n')status=15
    if(yn_out_spin_current_decomposed/='n'.or.yn_out_spin_current_micro/='n')status=15
    if(yn_out_mom_distr_rt/='n'.or.yn_out_rt_energy_components/='n')status=15
    call handle_agree(status)
    if(status/=0)ignored=sp_tddft_close()
  endif
  if(status==0)then
    handles(index)%label=trim(sysname)
    call comm_get_globalinfo(group,rank,n);handle_group=group;handle_mode=.true.
    call handle_unselect(index);token=candidate;J=-J
  else if(candidate>0_c_int64_t)then
    call handle_timer_leave(index)
  endif
  sp_tddft_handle_init=status
end function
integer(c_int) function sp_tddft_model_new(token,J) bind(C)
  use communication,only:comm_get_globalinfo
  integer(c_int64_t),intent(out)::token
  real(c_double),intent(out)::J(3)
  integer::status,index,group,rank,n
  integer(c_int64_t)::candidate
  token=0;J=0d0;sp_tddft_model_new=12
  if(.not.handle_main())return
  call comm_get_globalinfo(group,rank,n)
  status=0
  if(.not.handle_mode.or.group/=handle_group.or.associated(state))status=12
  if(registry_mode.and..not.model_spawning)status=12
  call handle_agree(status)
  if(status/=0)return
  call handle_reserve(index,candidate,status)
  if(status/=0)then
    sp_tddft_model_new=status;return
  endif
  write(handles(index)%label,'(a,i0,a,i0)')'native_tddft_h',index,'_g',handles(index)%generation
  call handle_timer_enter(index)
  call timer_initialize
  call open_initial_state(model,J,status,trim(handles(index)%label))
  if(status==0)then
    call handle_unselect(index);token=candidate;J=-J
  else
    call handle_timer_leave(index)
  endif
  sp_tddft_model_new=status
end function
integer(c_int) function sp_tddft_handle_step(token,Aold,Anew,J) bind(C)
  integer(c_int64_t),value::token
  real(c_double),intent(in)::Aold(3),Anew(3)
  real(c_double),intent(out)::J(3)
  integer::index,status
  J=0d0
  call handle_select(token,index,status);sp_tddft_handle_step=status
  if(status/=0)return
  sp_tddft_handle_step=sp_tddft_step(Aold,Anew,J);J=-J
  call handle_unselect(index)
end function
! State-token route with only local generation/bounds checks; same physical update.
integer(c_int) function sp_tddft_handle_step_performance(token,Aold,Anew,J) bind(C) result(status)
  integer(c_int64_t),value::token
  real(c_double),intent(in)::Aold(3),Anew(3)
  real(c_double),intent(out)::J(3)
  integer::index
  J=0d0
  call handle_select(token,index,status,.false.)
  if(status/=0)return
  status=sp_tddft_step_performance(Aold,Anew,J);J=-J
  call handle_unselect(index)
end function
integer(c_int) function sp_tddft_timer_protocol() bind(C)
  sp_tddft_timer_protocol=1
end function
integer(c_int) function sp_tddft_handle_timers(token,values) bind(C)
  integer(c_int64_t),value :: token
  real(c_double),intent(out) :: values(4)
  integer :: index,status
  values=0d0
  call handle_lookup(token,index,status)
  sp_tddft_handle_timers=status
  if(status/=0)return
  values=[timer_context_get(handles(index)%timers,LOG_TOTAL), &
          timer_context_get(handles(index)%timers,LOG_INIT_RT), &
          timer_context_get(handles(index)%timers,LOG_CALC_RHO), &
          timer_context_get(handles(index)%timers,LOG_CALC_HARTREE)]
end function
integer(c_int) function sp_tddft_handle_electron_number(token,value) bind(C)
  integer(c_int64_t),value :: token
  real(c_double),intent(out) :: value
  integer :: index,status
  value=0d0
  call handle_select(token,index,status)
  sp_tddft_handle_electron_number=status
  if(status/=0)return
  sp_tddft_handle_electron_number=sp_tddft_electron_number(value)
  call handle_unselect(index)
end function
integer(c_int) function sp_tddft_electron_number(value) bind(C)
  use dip,only:calc_dip
  use, intrinsic :: ieee_arithmetic,only:ieee_is_finite
  real(c_double),intent(out) :: value
  real(8) :: integrals(4)
  value=0d0;sp_tddft_electron_number=12
  if(.not.associated(state))return
  call calc_dip(state%info%icomm_r,state%lg,state%mg,state%system,state%rho_s,integrals)
  value=integrals(4)*state%system%Hvol
  sp_tddft_electron_number=4
  if(.not.ieee_is_finite(value))return
  if(value<=0d0)return
  sp_tddft_electron_number=0
end function
integer(c_int) function sp_tddft_handle_info(token,step,A,J) bind(C)
  integer(c_int64_t),value::token
  integer(c_int),intent(out)::step
  real(c_double),intent(out)::A(3),J(3)
  integer::index,status
  step=-1;A=0d0;J=0d0
  call handle_select(token,index,status);sp_tddft_handle_info=status
  if(status/=0)return
  step=state%itt;A=state%accepted_A
  if(state%itt==state%Mit)then
    J=-sum(state%worker_current(:,1:state%system%nspin),dim=2)
  else
    J=-state%rt%curr(:,state%itt)
  endif
  call handle_unselect(index)
end function
integer(c_int) function sp_tddft_handle_checkpoint(token,path,mode) bind(C)
  integer(c_int64_t),value::token
  character(c_char),intent(in)::path(*)
  integer(c_int),value::mode
  integer::index,status
  call handle_select(token,index,status);sp_tddft_handle_checkpoint=status
  if(status/=0)return
  sp_tddft_handle_checkpoint=sp_tddft_checkpoint(path,mode)
  call handle_unselect(index)
end function
integer(c_int) function sp_tddft_handle_transaction(token,mode) bind(C)
  integer(c_int64_t),value::token
  integer(c_int),value::mode
  integer::index,status
  call handle_select(token,index,status);sp_tddft_handle_transaction=status
  if(status/=0)return
  sp_tddft_handle_transaction=sp_tddft_step_transaction(mode)
  call handle_unselect(index)
end function
integer(c_int) function sp_tddft_handle_close(token) bind(C)
  integer(c_int64_t),value::token
  type(s_noncollinear_state) :: empty_noncollinear
  integer::index,status
  call handle_select(token,index,status);sp_tddft_handle_close=status
  if(status/=0)return
  sp_tddft_handle_close=sp_tddft_close()
  call handle_unselect(index)
  if(sp_tddft_handle_close==0)handles(index)%noncollinear=empty_noncollinear
end function
integer(c_int) function sp_tddft_handle_label(token,label,capacity) bind(C)
  integer(c_int64_t),value::token
  character(c_char),intent(out)::label(*)
  integer(c_int),value::capacity
  integer::index,status,i,n
  call handle_lookup(token,index,status);sp_tddft_handle_label=status
  if(status/=0)return
  n=len_trim(handles(index)%label)
  if(capacity<=n)status=13
  call handle_agree(status)
  if(status/=0)then
    sp_tddft_handle_label=status;return
  endif
  do i=1,n
    label(i)=handles(index)%label(i:i)
  enddo
  label(n+1)=c_null_char
end function

integer(c_int) function sp_tddft_context_mpi_protocol() bind(C)
  sp_tddft_context_mpi_protocol=1
end function
! Read-only diagnostics: live tracked children, registry slots, state%pending payload,
! allocated state. This is not total MPI usage or process memory.
integer(c_int) function sp_tddft_resources(values) bind(C)
  integer(c_int),intent(out)::values(4)
  values=0
  call comm_child_counts(values(1),values(2))
  block
  integer::i
  if(associated(state))then
    values(4)=1
    if(allocated(state%pending))values(3)=1
  endif
  if(allocated(handles))then
    do i=1,size(handles)
      if(.not.associated(handles(i)%context))cycle
      values(4)=1
      if(allocated(handles(i)%context%pending))values(3)=1
    enddo
  endif
  end block
  sp_tddft_resources=0
end function
! Protocol1: statuses10..14 reject collectively before numerical/latch changes.
integer(c_int) function sp_tddft_step_transaction_protocol() bind(C)
  sp_tddft_step_transaction_protocol=1
end function
integer(c_int) function sp_tddft_step_transaction(mode) bind(C)
  use communication,only:comm_get_max,comm_get_globalinfo
  integer(c_int),value::mode
  integer::status,reference,group,rank,count
  call comm_get_globalinfo(group,rank,count)
  status=0;reference=mode
  call comm_bcast(reference,group)
  if(mode/=reference.or.mode<0.or.mode>3)status=11
  if(.not.associated(state))status=12
  call comm_get_max(status,group)
  if(status/=0)then
    sp_tddft_step_transaction=status;return
  endif
  if(mode==0)then
    if(allocated(state%step_saved).or.allocated(state%pending))status=14
    if(yn_out_dns_rt/='n'.or.yn_out_dns_ac_je/='n'.or.yn_out_micro_je/='n')status=15
    if(yn_out_elf_rt/='n'.or.yn_out_estatic_rt/='n'.or.yn_out_current_decomposed/='n')status=15
    if(projection_option/='no'.or.dm_unfold_option/='no'.or.yn_out_rvf_rt/='n')status=15
  else if(mode/=2)then
    if(.not.allocated(state%step_saved))status=14
    if(mode==1.and.state%step_restored)status=14
  endif
  call comm_get_max(status,group)
  if(status/=0)then
    sp_tddft_step_transaction=status;return
  endif
  select case(mode)
  case(0)
    allocate(state%step_saved,state%step_saved_xc)
    call checkpoint_copy(state%step_saved,.false.)
    state%step_saved_xc=state%system%xc_payload;state%step_saved_cval=state%xc_func%cval
    state%step_restored=.false.
    if(rank==0)call step_output_capture(state%step_outputs,state%ofl,status)
    call comm_get_max(status,group)
    if(status/=0)deallocate(state%step_saved,state%step_saved_xc)
  case(1)
    call checkpoint_copy(state%step_saved,.true.)
    state%system%xc_payload=state%step_saved_xc;state%xc_func%cval=state%step_saved_cval
    call update_kvector_nonlocalpt(state%info%ik_s,state%info%ik_e,state%system,state%ppg)
    state%step_restored=.true.
    if(rank==0)call step_output_close(state%step_outputs,status)
  case(2)
    if(allocated(state%step_saved))then
      if(rank==0)then
        call step_output_reopen(state%step_outputs,status,state%ofl)
        if(status==0)call step_output_flush(state%step_outputs,status)
      endif
      call comm_get_max(status,group)
      if(status==0)then
        deallocate(state%step_saved,state%step_saved_xc)
        state%step_restored=.false.
      endif
    endif
  case(3)
    if(rank==0)call step_output_reopen(state%step_outputs,status,state%ofl)
  end select
  call comm_get_max(status,group)
  sp_tddft_step_transaction=status
end function

! Local feature query; no MPI initialization or state selection required.
integer(c_int) function sp_tddft_soc_spatial_protocol() bind(C)
#ifdef USE_OPENACC
  sp_tddft_soc_spatial_protocol=0
#else
  sp_tddft_soc_spatial_protocol=1
#endif
end function

integer(c_int) function sp_tddft_error_protocol() bind(C)
  sp_tddft_error_protocol=1
end function

integer(c_int) function sp_tddft_recovery_protocol() bind(C)
  sp_tddft_recovery_protocol=1
end function

integer(c_int) function sp_tddft_checkpoint_protocol() bind(C)
  sp_tddft_checkpoint_protocol=1
end function

integer(c_int) function sp_tddft_init(path,timestep,capacity,A,decomposition,J) bind(C)
  character(c_char),intent(in)::path(*)
  real(c_double),value::timestep
  integer(c_int),value::capacity
  real(c_double),intent(in)::A(3)
  integer(c_int),intent(in)::decomposition(5)
  real(c_double),intent(out)::J(3)
  if(.not.handle_timer_active)call legacy_timer_enter()
  sp_tddft_init=tddft_init_impl(path,timestep,capacity,A,decomposition,J)
  if(.not.handle_timer_active)call legacy_timer_leave()
end function
integer(c_int) function tddft_init_impl(path, timestep, capacity, A, decomposition, J)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use communication, only: comm_get_max, comm_get_globalinfo
  use inputoutput, only: read_input
  use parallelization, only: setup_parallel
  use math_constants, only: set_math_constants
  use sym_sub, only: read_sw_symmetry
  character(c_char),intent(in) :: path(*)
  real(c_double),value :: timestep
  integer(c_int),value :: capacity
  real(c_double),intent(in) :: A(3)
  integer(c_int),intent(in) :: decomposition(5)
  real(c_double),intent(out) :: J(3)
  character(4096) :: filename
  integer :: i,status,reference_capacity,reference_decomposition(5),admission_group,admission_rank,admission_size,child_owner
  real(8) :: reference_dt,reference_A(3)
  character(4096) :: reference_filename
  logical :: exists,first_init
  call comm_get_globalinfo(admission_group,admission_rank,admission_size)
  J=0d0
  status=0
  if(used.and..not.registering)status=12
  call comm_get_max(status,admission_group)
  if(status/=0)then
    tddft_init_impl=status;return
  endif
  filename=''
  do i=1,len(filename)
    if(path(i)==c_null_char) exit
    filename(i:i)=path(i)
  end do
  if(i>len(filename).or.len_trim(filename)==0)status=10
  if(.not.ieee_is_finite(timestep).or.timestep<=0d0.or.capacity<1.or.capacity>huge(capacity)-2)status=10
  do i=1,3
    if(.not.ieee_is_finite(A(i)))status=10
  enddo
  if(any(decomposition<1))status=10
  if(product(real(decomposition,8))/=real(admission_size,8))status=10
  if(.not.associated(pending_gs))then
    inquire(file=resolve_path(trim(filename)),exist=exists)
    if(.not.exists)status=10
  endif
  call comm_get_max(status,admission_group)
  if(status/=0)then
    tddft_init_impl=status;return
  endif
  reference_dt=timestep;reference_capacity=capacity;reference_A=A
  reference_decomposition=decomposition;reference_filename=filename
  call comm_bcast(reference_dt,admission_group)
  call comm_bcast(reference_capacity,admission_group)
  call comm_bcast(reference_A,admission_group)
  call comm_bcast(reference_decomposition,admission_group)
  call comm_bcast(reference_filename,admission_group)
  if(reference_dt/=timestep.or.reference_capacity/=capacity.or.any(reference_A/=A).or. &
     any(reference_decomposition/=decomposition).or.reference_filename/=filename)status=11
  call comm_get_max(status,admission_group)
  if(status/=0)then
    tddft_init_impl=status;return
  endif
  ! Only actual initialization consumes the process-owned lifetime.
  first_init=.not.used
  used=.true.
  worker_embedded=.true.
  worker_initial_A=A; embedded_dt=timestep; embedded_steps=capacity
  embedded_decomposition=decomposition
  if(first_init)then
    call set_math_constants
    call setup_parallel
  endif
  allocate(model)
  if(associated(pending_gs))then
    allocate(model%seed)
    model%seed=pending_gs
    call model%seed%inputs%activate(status)
    theory='tddft_pulse';yn_restart='n';yn_self_checkpoint='n';yn_reset_step_restart='n'
    sysname='rt';base_directory='./';directory_read_data='./data_for_restart/'
    propagator='middlepoint';trans_longi='tr';yn_predictor_corrector='n'
    call worker_configure()
    call model%parser%values%capture()
  else
    call read_input(trim(filename),model%parser)
#ifdef SALMONPY_FAST_CORE
    call worker_configure()
    call model%parser%values%capture()
#endif
  endif
  calc_mode='RT'
  call read_sw_symmetry(yn_symmetry)
  if(first_init)call timer_initialize
  model%initial_A=worker_initial_A
  call open_initial_state(model,J,status,'')
  tddft_init_impl=status
end function
subroutine open_initial_state(owner,J,status,name)
  use communication,only:comm_get_max
  use inputoutput,only:initialize_inputoutput_units
  type(rt_model),pointer,intent(in)::owner
  real(8),intent(out)::J(3)
  integer,intent(out)::status
  character(*),intent(in)::name
  character(len(sysname))::saved_name
  logical::exists
  integer::child_owner,restore_status
  type(rt_model),pointer :: previous
  status=0;J=0d0
  if(.not.associated(owner))status=12
  call handle_agree(status)
  if(status/=0)return
  previous=>active_model
  call select_model(owner,.true.,status)
  if(status/=0)return
  saved_name=sysname
  if(len_trim(name)>0)saved_name=name
  if(len_trim(name)>0.or.handle_initializing)then
    inquire(file=resolve_path(trim(base_directory)//trim(saved_name)//'_rt.data'),exist=exists)
    if(exists)status=1
    inquire(file=resolve_path(trim(base_directory)//trim(saved_name)//'_rt_energy.data'),exist=exists)
    if(exists)status=1
    if(yn_spinorbit=='y')then
      inquire(file=resolve_path(trim(base_directory)//trim(saved_name)//'_rt_spin.data'),exist=exists)
      if(exists)status=1
    endif
  endif
  call handle_agree(status)
  if(status/=0)then
    if(associated(previous))call select_model(previous,.false.,restore_status)
    return
  endif
  quiet=.false.
  call comm_begin_children(child_owner,status)
  call comm_get_max(status,nproc_group_global)
  if(status/=0)then
    call comm_end_children(child_owner)
    call comm_release_children(child_owner)
    status=16;return
  endif
  allocate(state);state%mpi_child_owner=child_owner;state%model=>owner
  saved_name=sysname
  if(len_trim(name)>0)sysname=name
if(.not.model_timers_alive)call timer_begin(LOG_TOTAL)

call initialization_rt( state%Mit, state%system, state%energy, state%ewald, state%rt, state%md, &
                        state%singlescale,  &
                        state%stencil, state%fg, state%poisson,  &
                        state%lg, state%mg,   &
                        state%info,  &
                        state%xc_func, state%ofl,  &
                        state%srg, state%srg_scalar,  &
                        state%spsi_in, state%spsi_out, state%tpsi, state%rho, state%rho_jm, state%rho_s,  &
                        state%V_local, state%Vbox, state%Vh, state%Vh_stock1, state%Vh_stock2, state%Vxc, state%Vpsl,&
                        state%pp, state%ppg, state%ppn, state%unfold, memory_gs=owner%seed )

if(worker_requested()) then
  state%system%vec_Ac=worker_initial_A
  call update_kvector_nonlocalpt(state%info%ik_s,state%info%ik_e,state%system,state%ppg)
  call calc_current(state%system,state%mg,state%stencil,state%info,state%srg,state%spsi_in,state%ppg, &
     & state%worker_current(:,1:state%system%nspin))
  state%spsi_in%update_zwf_overlap=.true.
  J=sum(state%worker_current(:,1:state%system%nspin),dim=2)
end if

  sysname=saved_name
  call comm_end_children(child_owner)
  state%itt=state%Mit
  call comm_sync_all
  if(.not.model_timers_alive)then
    call timer_enable_sub
    call timer_begin(LOG_RT_ITERATION)
    model_timers_alive=.true.
  endif

  state%accepted_A=worker_initial_A
  call owner%runtime_inputs%capture()
  owner%runtime_ready=.true.
  owner%references=owner%references+1
  status=0
end subroutine
integer(c_int) function sp_tddft_step(Aold,Anew,J) bind(C)
  use communication,only:comm_resume_children,comm_end_children,comm_get_max
  real(c_double),intent(in) :: Aold(3),Anew(3)
  real(c_double),intent(out) :: J(3)
  integer :: status
  logical :: entered
  J=0d0;status=0
  if(.not.associated(state))status=12
  call comm_get_max(status,nproc_group_global)
  sp_tddft_step=status
  if(status/=0)return
  call comm_resume_children(state%mpi_child_owner,status)
  entered=status==0
  call comm_get_max(status,nproc_group_global)
  if(status==0)then
    if(.not.handle_timer_active)call legacy_timer_enter()
    sp_tddft_step=tddft_step_impl(Aold,Anew,J)
    if(.not.handle_timer_active)call legacy_timer_leave()
  else
    sp_tddft_step=16
  endif
  if(entered)call comm_end_children(state%mpi_child_owner)
end function
integer(c_int) function tddft_step_impl(Aold,Anew,J)
  use pseudo_pt_plusU_sub,only:PLUS_U_ON
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use communication, only: comm_get_max, comm_get_globalinfo
  real(c_double),intent(in) :: Aold(3),Anew(3)
  real(c_double),intent(out) :: J(3)
  integer :: status,reference_step,admission_group,admission_rank,admission_size,i
  real(8) :: reference_old(3),reference_new(3),saved_xc_cval
  type(checkpoint_payload) :: before_step
  type(s_xc_operator_payload) :: saved_xc_payload
  call comm_get_globalinfo(admission_group,admission_rank,admission_size)
  J=0d0
  status=0
  if(.not.associated(state))status=12
  call comm_get_max(status,admission_group)
  if(status/=0)then
    tddft_step_impl=status;return
  endif
  do i=1,3
    if(.not.ieee_is_finite(Aold(i)).or..not.ieee_is_finite(Anew(i)))status=10
  enddo
  if(any(Aold/=state%accepted_A))status=10
  if(state%itt>=nt)status=13
  if(allocated(state%pending))status=14
  call comm_get_max(status,admission_group)
  if(status/=0)then
    tddft_step_impl=status;return
  endif
  reference_old=Aold;reference_new=Anew;reference_step=state%itt
  call comm_bcast(reference_old,admission_group)
  call comm_bcast(reference_new,admission_group)
  call comm_bcast(reference_step,admission_group)
  if(any(reference_old/=Aold).or.any(reference_new/=Anew).or.reference_step/=state%itt)status=11
  call comm_get_max(status,admission_group)
  if(status/=0)then
    tddft_step_impl=status;return
  endif
  ! The transaction covers fixed-ion periodic scalar/SOC semilocal propagation.
  select case(trim(xc))
  case('pz','tbmbj','libxc_pz','libxc_pbe')
  case default
    status=15
  end select
  if(xname/='none'.or.cname/='none')status=15
  if(.not.worker_spin_supported().or.yn_md/='n'.or.yn_jm/='n')status=15
  if(iperiodic/=3.or.propagator/='middlepoint'.or.trans_longi/='tr')status=15
  if(projection_option/='no'.or.dm_unfold_option/='no'.or.state%singlescale%flag_use)status=15
  if(yn_out_rvf_rt/='n'.or.PLUS_U_ON.or.gram_schmidt_interval>=1)status=15
  call comm_get_max(status,admission_group)
  if(status/=0)then
    tddft_step_impl=status;return
  endif
  call checkpoint_copy(before_step,.false.)
  saved_xc_cval=state%xc_func%cval
  saved_xc_payload=state%system%xc_payload
  state%rt%recover_errors=.true.
  state%rt%step_error=0
  call evolve_endpoint(Anew)

  status=state%rt%step_error
  call comm_get_max(status,admission_group)
  state%rt%recover_errors=.false.
  state%rt%step_error=0
  if(status/=0)then
    call checkpoint_copy(before_step,.true.)
    state%xc_func%cval=saved_xc_cval
    call update_kvector_nonlocalpt(state%info%ik_s,state%info%ik_e,state%system,state%ppg)
    ! The complete in-memory wavefunctions include their accepted halo bytes
    ! and original overlap flags; disk LOAD's normalization is not appropriate.
    state%system%xc_payload=saved_xc_payload
    tddft_step_impl=status
    return
  endif
  J=state%rt%curr(:,state%itt)

  state%accepted_A=Anew
  tddft_step_impl=0
end function
! Explicit unchecked fixed-ion semilocal update. Caller owns rank order/inputs.
! No diagnostic agreement, full-state backup or automatic failure restoration.
integer(c_int) function sp_tddft_step_performance(Aold,Anew,J) bind(C) result(status)
  use communication,only:comm_resume_children,comm_end_children
  real(c_double),intent(in)::Aold(3),Anew(3)
  real(c_double),intent(out)::J(3)
  status=12;J=0d0
  if(.not.associated(state))return
  status=13
  if(state%itt>=nt)return
  status=14
  if(allocated(state%pending))return
  call comm_resume_children(state%mpi_child_owner,status)
  if(status/=0)return
  if(.not.handle_timer_active)call legacy_timer_enter()
  state%rt%recover_errors=.false.;state%rt%step_error=0
  call evolve_endpoint(Anew)
  J=state%rt%curr(:,state%itt);state%accepted_A=Anew;status=0
  if(.not.handle_timer_active)call legacy_timer_leave()
  call comm_end_children(state%mpi_child_owner)
end function
subroutine evolve_endpoint(Anew)
  real(c_double),intent(in)::Anew(3)
  state%itt=state%itt+1
  state%rt%Ac_ext(:,state%itt)=Anew
  state%rt%Ac_tot(:,state%itt)=Anew
  if(mod(state%itt,2)==1)then
    call time_evolution_step(state%Mit,nt,state%itt,state%lg,state%mg,state%system,state%rt,state%info,state%stencil,state%xc_func &
     & ,state%srg,state%srg_scalar,state%pp,state%ppg,state%ppn,state%spsi_in,state%spsi_out,state%tpsi &
     & ,state%rho,state%rho_jm,state%rho_s,state%V_local,state%Vbox,state%Vh,state%Vh_stock1,state%Vh_stock2,state%Vxc &
     & ,state%Vpsl,state%fg,state%energy,state%ewald,state%md,state%ofl,state%poisson,state%singlescale,state%unfold)
  else
    call time_evolution_step(state%Mit,nt,state%itt,state%lg,state%mg,state%system,state%rt,state%info,state%stencil,state%xc_func &
     & ,state%srg,state%srg_scalar,state%pp,state%ppg,state%ppn,state%spsi_out,state%spsi_in,state%tpsi &
     & ,state%rho,state%rho_jm,state%rho_s,state%V_local,state%Vbox,state%Vh,state%Vh_stock1,state%Vh_stock2,state%Vxc &
     & ,state%Vpsl,state%fg,state%energy,state%ewald,state%md,state%ofl,state%poisson,state%singlescale,state%unfold)
  end if

end subroutine
integer(c_int) function sp_tddft_close() bind(C)
  if(.not.handle_timer_active)call legacy_timer_enter()
  sp_tddft_close=tddft_close_impl()
  if(.not.handle_timer_active)call legacy_timer_leave()
end function
integer(c_int) function tddft_close_impl()
  integer::output_status
  type(rt_model),pointer :: owner
  tddft_close_impl=1
  if(.not.associated(state)) return
  if(comm_is_root(nproc_id_global))then
    if(allocated(state%step_saved))then
      call step_output_close(state%step_outputs,output_status)
    else
      close(state%ofl%fh_rt)
      close(state%ofl%fh_rt_energy)
      if(yn_spinorbit=='y')close(state%ofl%fh_rt_spin)
    endif
  endif
  if(allocated(state%step_saved))deallocate(state%step_saved,state%step_saved_xc)
  if(allocated(state%pending))deallocate(state%pending)
  if(model_timers_alive)then
      call timer_end(LOG_RT_ITERATION)
      call timer_disable_sub
      call timer_end(LOG_TOTAL)
      model_timers_alive=.false.
  endif
  ! The embedded path never opens laser unit30; it may belong to the host.
  call dealloc_cache(state%srg)
  call dealloc_cache(state%srg_scalar)
  call finalize_xc(state%xc_func)
  call comm_release_children(state%mpi_child_owner)
  owner=>state%model;owner%references=owner%references-1
  deallocate(state);nullify(state)
  call release_model(owner)
  ! Other SALMON module configuration remains process-owned. Do not reopen.
  tddft_close_impl=0
end function

! Native ABI checkpoint v3. Caller authenticates the file, binary, GS and configuration.
! XC identity and scalar cval are stored; Libxc resources are rebuilt by initialization.
! TB-mBJ's automatic c is evaluation-local, not a persistent mutable XC field.
! All ranks call collectively with DISTINCT paths and an unchanged process topology.
! restore=0 creates a new file; restore=1 loads into a fresh initialized context, without stepping.
integer(c_int) function sp_tddft_checkpoint(path,restore) bind(C)
  use parallelization,only:nproc_size_global
  use communication,only:comm_get_max
  use, intrinsic :: ieee_arithmetic,only:ieee_is_finite
  character(c_char),intent(in)::path(*)
  integer(c_int),value::restore
  type(checkpoint_payload),allocatable::copy
  character(4096)::filename
  character(32)::magic
  character(256)::xc_names(3),expected_xc_names(3)
  character(1)::extra
  integer::i,u,ios,close_status,status,header(33),expected(33),allocation_status,mode
  real(8)::saved_dt,cell(3,3),saved_input_cval,saved_xc_cval
  integer::xc_codes(4),expected_xc_codes(4)
  logical::reading,opened,xc_flags(4),expected_xc_flags(4)
  sp_tddft_checkpoint=1
  status=0;opened=.false.;reading=restore==1.or.restore==2
  if(.not.associated(state))status=2
  call comm_get_max(status,nproc_group_global)
  if(status/=0)then
    sp_tddft_checkpoint=status;return
  endif
  mode=restore
  call comm_bcast(mode,nproc_group_global)
  if(mode/=restore.or.restore<0.or.restore>4)status=2
  if(allocated(state%step_saved).and.restore/=0)status=14
  if(restore==3.and..not.allocated(state%pending))status=2
  call comm_get_max(status,nproc_group_global)
  if(status/=0)then
    sp_tddft_checkpoint=status;return
  endif
  if(restore==4)then
    if(allocated(state%pending))deallocate(state%pending)
    sp_tddft_checkpoint=0;return
  endif
  if(restore==3)then
    if(.not.allocated(state%pending))return
    call checkpoint_install(state%pending,state%pending_xc_cval)
    deallocate(state%pending);sp_tddft_checkpoint=0;return
  endif
  if(restore==2.and.allocated(state%pending))status=2
  if(restore<0.or.restore>4)status=2
  select case(trim(xc))
  case('pz','tbmbj')
#ifdef USE_LIBXC
  case('libxc_pz','libxc_pbe')
#endif
  case default
    status=2
  end select
  if(xname/='none'.or.cname/='none')status=2
  if(.not.worker_spin_supported().or.yn_md/='n')status=2
  if(.not.ieee_is_finite(cval).or..not.ieee_is_finite(state%xc_func%cval))status=2
  if(iperiodic/=3.or.propagator/='middlepoint'.or.trans_longi/='tr')status=2
  if(projection_option/='no'.or.dm_unfold_option/='no'.or.state%singlescale%flag_use)status=2
  if(restore==1.and.state%itt/=state%Mit)status=2
  filename=''
  do i=1,len(filename)
    if(path(i)==c_null_char)exit
    filename(i:i)=path(i)
  enddo
  if(i>len(filename).or.len_trim(filename)==0)status=3
  call comm_get_max(status,nproc_group_global)
  if(status/=0)then
    sp_tddft_checkpoint=status;return
  endif
  allocate(copy,stat=allocation_status)
  if(allocation_status/=0)status=4
  call comm_get_max(status,nproc_group_global)
  if(status/=0)then
    sp_tddft_checkpoint=status;return
  endif
  state%checkpoint_array_limit=max(64_8*(int(nt,8)+1_8),64_8*int(state%system%ngrid,8))
  if(allocated(state%spsi_in%zwf))state%checkpoint_array_limit=max(state%checkpoint_array_limit,size(state%spsi_in%zwf,kind=8))
  if(allocated(state%spsi_in%rwf))state%checkpoint_array_limit=max(state%checkpoint_array_limit,size(state%spsi_in%rwf,kind=8))
  call checkpoint_copy(copy,.false.)
  expected=[3,nproc_size_global,nproc_id_global,nt,state%system%nspin,state%system%no,state%system%nk, &
    state%lg%num,state%mg%is,state%mg%ie,state%mg%is_array,state%mg%ie_array, &
    state%info%io_s,state%info%io_e,state%info%ik_s,state%info%ik_e,state%info%im_s,state%info%im_e, &
    storage_size(1d0),storage_size(1),storage_size((0d0,0d0)),0,0]
  if(reading)then
    open(newunit=u,file=resolve_path(trim(filename)),status='old',access='stream',form='unformatted',action='read',iostat=ios)
  else
    open(newunit=u,file=resolve_path(trim(filename)),status='new',access='stream',form='unformatted',action='write',iostat=ios)
  endif
  opened=ios==0
  if(opened)then
    magic='SALMONPY_TDDFT_NATIVE_V3'
    header=expected;saved_dt=dt;cell=state%system%primitive_a
    expected_xc_names=[character(256)::xc,xname,cname]
    expected_xc_codes=[state%xc_func%xctype,state%xc_func%ispin]
    expected_xc_flags=[state%xc_func%use_gradient,state%xc_func%use_laplacian, &
      state%xc_func%use_kinetic_energy,state%xc_func%use_current]
    xc_names=expected_xc_names;xc_codes=expected_xc_codes;xc_flags=expected_xc_flags
    saved_input_cval=cval;saved_xc_cval=state%xc_func%cval
    if(reading)then
      read(u,iostat=ios)magic,header,saved_dt,cell
      if(ios==0)then
        if(magic/='SALMONPY_TDDFT_NATIVE_V3'.or.any(header/=expected))ios=1
        if(saved_dt/=dt.or.any(cell/=state%system%primitive_a))ios=1
      endif
    else
      write(u,iostat=ios)magic,header,saved_dt,cell
    endif
    if(ios==0)then
      if(reading)then
        read(u,iostat=ios)xc_names,xc_codes,xc_flags,saved_input_cval,saved_xc_cval
        if(ios==0)then
          if(any(xc_names/=expected_xc_names).or.any(xc_codes/=expected_xc_codes))ios=1
          if(any(xc_flags.neqv.expected_xc_flags))ios=1
          if(saved_input_cval/=cval.or..not.ieee_is_finite(saved_xc_cval))ios=1
        endif
      else
        write(u,iostat=ios)xc_names,xc_codes,xc_flags,saved_input_cval,saved_xc_cval
      endif
    endif
    if(ios==0)call checkpoint_io(u,reading,copy,ios)
    if(ios==0)then
      magic='SALMONPY_TDDFT_COMPLETE_V3'
      if(reading)then
        read(u,iostat=ios)magic
        if(ios==0.and.magic/='SALMONPY_TDDFT_COMPLETE_V3')ios=1
        if(ios==0)then
          read(u,iostat=ios)extra
          if(is_iostat_end(ios))then
            ios=0
          else
            ios=1
          endif
        endif
      else
        write(u,iostat=ios)magic
      endif
    endif
    close(u,iostat=close_status)
    if(close_status/=0)ios=close_status
  endif
  if(ios/=0)status=4
  if(reading.and.status==0)then
    if(copy%Mit/=state%Mit.or.copy%itt<copy%Mit.or.copy%itt>nt)status=5
    do i=1,3
      if(.not.ieee_is_finite(copy%accepted_A(i)))status=5
      if(.not.ieee_is_finite(copy%fields(i,1)))status=5
    enddo
  endif
  call comm_get_max(status,nproc_group_global)
  if(status==0.and.reading)then
    if(restore==2)then
      call move_alloc(copy,state%pending)
      state%pending_xc_cval=saved_xc_cval
    else
      call checkpoint_install(copy,saved_xc_cval)
    endif
  endif
  sp_tddft_checkpoint=status
end function

subroutine checkpoint_install(copy,saved_xc_cval)
  type(checkpoint_payload),intent(inout)::copy
  real(8),intent(in)::saved_xc_cval
    call checkpoint_copy(copy,.true.)
    ! Preserve initialized func/info Libxc handles; only the scalar coefficient is restored.
    state%xc_func%cval=saved_xc_cval
    ! PP phase factors and communication caches are reconstructed, not serialized handles.
    call update_kvector_nonlocalpt(state%info%ik_s,state%info%ik_e,state%system,state%ppg)
    state%spsi_in%update_zwf_overlap=.true.
    state%spsi_out%update_zwf_overlap=.true.
    state%tpsi%update_zwf_overlap=.true.
end subroutine

subroutine checkpoint_copy(copy,restore)
  type(checkpoint_payload),intent(inout)::copy
  logical,intent(in)::restore
  type(s_noncollinear_state) :: replacement
  if(restore)then
    replacement=copy%noncollinear
    call swap_noncollinear_state(replacement)
    state%spsi_in=copy%spsi_in
    state%spsi_out=copy%spsi_out
    state%tpsi=copy%tpsi
    state%rho=copy%rho
    state%rho_jm=copy%rho_jm
    state%Vh=copy%Vh
    state%Vh_stock1=copy%Vh_stock1
    state%Vh_stock2=copy%Vh_stock2
    state%Vbox=copy%Vbox
    state%Vpsl=copy%Vpsl
    if(allocated(copy%rho_s))then
      state%rho_s=copy%rho_s
    else
      if(allocated(state%rho_s))deallocate(state%rho_s)
    endif
    if(allocated(copy%V_local))then
      state%V_local=copy%V_local
    else
      if(allocated(state%V_local))deallocate(state%V_local)
    endif
    if(allocated(copy%Vxc))then
      state%Vxc=copy%Vxc
    else
      if(allocated(state%Vxc))deallocate(state%Vxc)
    endif
    state%energy=copy%energy
    state%poisson=copy%poisson
    state%Mit=copy%Mit
    state%itt=copy%itt
    state%worker_A=copy%worker_A
    state%worker_current=copy%worker_current
    state%rt%Dp0_e=copy%rt%Dp0_e
    if(allocated(copy%rt%dDp_e))then
      state%rt%dDp_e=copy%rt%dDp_e
    else
      if(allocated(state%rt%dDp_e))deallocate(state%rt%dDp_e)
    endif
    if(allocated(copy%rt%Dp_e))then
      state%rt%Dp_e=copy%rt%Dp_e
    else
      if(allocated(state%rt%Dp_e))deallocate(state%rt%Dp_e)
    endif
    if(allocated(copy%rt%Dp_i))then
      state%rt%Dp_i=copy%rt%Dp_i
    else
      if(allocated(state%rt%Dp_i))deallocate(state%rt%Dp_i)
    endif
    if(allocated(copy%rt%Qp_e))then
      state%rt%Qp_e=copy%rt%Qp_e
    else
      if(allocated(state%rt%Qp_e))deallocate(state%rt%Qp_e)
    endif
    if(allocated(copy%rt%rIe))then
      state%rt%rIe=copy%rt%rIe
    else
      if(allocated(state%rt%rIe))deallocate(state%rt%rIe)
    endif
    if(allocated(copy%rt%curr))then
      state%rt%curr=copy%rt%curr
    else
      if(allocated(state%rt%curr))deallocate(state%rt%curr)
    endif
    if(allocated(copy%rt%E_ext))then
      state%rt%E_ext=copy%rt%E_ext
    else
      if(allocated(state%rt%E_ext))deallocate(state%rt%E_ext)
    endif
    if(allocated(copy%rt%E_ind))then
      state%rt%E_ind=copy%rt%E_ind
    else
      if(allocated(state%rt%E_ind))deallocate(state%rt%E_ind)
    endif
    if(allocated(copy%rt%E_tot))then
      state%rt%E_tot=copy%rt%E_tot
    else
      if(allocated(state%rt%E_tot))deallocate(state%rt%E_tot)
    endif
    if(allocated(copy%rt%Ac_ext))then
      state%rt%Ac_ext=copy%rt%Ac_ext
    else
      if(allocated(state%rt%Ac_ext))deallocate(state%rt%Ac_ext)
    endif
    if(allocated(copy%rt%Ac_ind))then
      state%rt%Ac_ind=copy%rt%Ac_ind
    else
      if(allocated(state%rt%Ac_ind))deallocate(state%rt%Ac_ind)
    endif
    if(allocated(copy%rt%Ac_tot))then
      state%rt%Ac_tot=copy%rt%Ac_tot
    else
      if(allocated(state%rt%Ac_tot))deallocate(state%rt%Ac_tot)
    endif
    if(allocated(copy%rt%zc))then
      state%rt%zc=copy%rt%zc
    else
      if(allocated(state%rt%zc))deallocate(state%rt%zc)
    endif
    if(allocated(copy%rt%vloc_t))then
      state%rt%vloc_t=copy%rt%vloc_t
    else
      if(allocated(state%rt%vloc_t))deallocate(state%rt%vloc_t)
    endif
    if(allocated(copy%rt%vloc_new))then
      state%rt%vloc_new=copy%rt%vloc_new
    else
      if(allocated(state%rt%vloc_new))deallocate(state%rt%vloc_new)
    endif
    if(allocated(copy%rt%vloc_old))then
      state%rt%vloc_old=copy%rt%vloc_old
    else
      if(allocated(state%rt%vloc_old))deallocate(state%rt%vloc_old)
    endif
    if(allocated(copy%rt%rho0_s))then
      state%rt%rho0_s=copy%rt%rho0_s
    else
      if(allocated(state%rt%rho0_s))deallocate(state%rt%rho0_s)
    endif
    state%rt%vonf=copy%rt%vonf
    state%rt%j_e=copy%rt%j_e
    state%rt%E_old=copy%rt%E_old
    state%rt%rho0_loaded=copy%rt%rho0_loaded
    state%accepted_A=copy%accepted_A
    state%system%vec_Ac=copy%fields(:,1)
    state%system%vec_Ac_ext=copy%fields(:,2)
    state%system%vec_E=copy%fields(:,3)
    state%system%vec_E_ext=copy%fields(:,4)
  else
    call snapshot_noncollinear_state(copy%noncollinear)
    copy%spsi_in=state%spsi_in
    copy%spsi_out=state%spsi_out
    copy%tpsi=state%tpsi
    copy%rho=state%rho
    copy%rho_jm=state%rho_jm
    copy%Vh=state%Vh
    copy%Vh_stock1=state%Vh_stock1
    copy%Vh_stock2=state%Vh_stock2
    copy%Vbox=state%Vbox
    copy%Vpsl=state%Vpsl
    if(allocated(state%rho_s))then
      copy%rho_s=state%rho_s
    else
      if(allocated(copy%rho_s))deallocate(copy%rho_s)
    endif
    if(allocated(state%V_local))then
      copy%V_local=state%V_local
    else
      if(allocated(copy%V_local))deallocate(copy%V_local)
    endif
    if(allocated(state%Vxc))then
      copy%Vxc=state%Vxc
    else
      if(allocated(copy%Vxc))deallocate(copy%Vxc)
    endif
    copy%energy=state%energy
    copy%poisson=state%poisson
    copy%Mit=state%Mit
    copy%itt=state%itt
    copy%worker_A=state%worker_A
    copy%worker_current=state%worker_current
    copy%rt%Dp0_e=state%rt%Dp0_e
    if(allocated(state%rt%dDp_e))then
      copy%rt%dDp_e=state%rt%dDp_e
    else
      if(allocated(copy%rt%dDp_e))deallocate(copy%rt%dDp_e)
    endif
    if(allocated(state%rt%Dp_e))then
      copy%rt%Dp_e=state%rt%Dp_e
    else
      if(allocated(copy%rt%Dp_e))deallocate(copy%rt%Dp_e)
    endif
    if(allocated(state%rt%Dp_i))then
      copy%rt%Dp_i=state%rt%Dp_i
    else
      if(allocated(copy%rt%Dp_i))deallocate(copy%rt%Dp_i)
    endif
    if(allocated(state%rt%Qp_e))then
      copy%rt%Qp_e=state%rt%Qp_e
    else
      if(allocated(copy%rt%Qp_e))deallocate(copy%rt%Qp_e)
    endif
    if(allocated(state%rt%rIe))then
      copy%rt%rIe=state%rt%rIe
    else
      if(allocated(copy%rt%rIe))deallocate(copy%rt%rIe)
    endif
    if(allocated(state%rt%curr))then
      copy%rt%curr=state%rt%curr
    else
      if(allocated(copy%rt%curr))deallocate(copy%rt%curr)
    endif
    if(allocated(state%rt%E_ext))then
      copy%rt%E_ext=state%rt%E_ext
    else
      if(allocated(copy%rt%E_ext))deallocate(copy%rt%E_ext)
    endif
    if(allocated(state%rt%E_ind))then
      copy%rt%E_ind=state%rt%E_ind
    else
      if(allocated(copy%rt%E_ind))deallocate(copy%rt%E_ind)
    endif
    if(allocated(state%rt%E_tot))then
      copy%rt%E_tot=state%rt%E_tot
    else
      if(allocated(copy%rt%E_tot))deallocate(copy%rt%E_tot)
    endif
    if(allocated(state%rt%Ac_ext))then
      copy%rt%Ac_ext=state%rt%Ac_ext
    else
      if(allocated(copy%rt%Ac_ext))deallocate(copy%rt%Ac_ext)
    endif
    if(allocated(state%rt%Ac_ind))then
      copy%rt%Ac_ind=state%rt%Ac_ind
    else
      if(allocated(copy%rt%Ac_ind))deallocate(copy%rt%Ac_ind)
    endif
    if(allocated(state%rt%Ac_tot))then
      copy%rt%Ac_tot=state%rt%Ac_tot
    else
      if(allocated(copy%rt%Ac_tot))deallocate(copy%rt%Ac_tot)
    endif
    if(allocated(state%rt%zc))then
      copy%rt%zc=state%rt%zc
    else
      if(allocated(copy%rt%zc))deallocate(copy%rt%zc)
    endif
    if(allocated(state%rt%vloc_t))then
      copy%rt%vloc_t=state%rt%vloc_t
    else
      if(allocated(copy%rt%vloc_t))deallocate(copy%rt%vloc_t)
    endif
    if(allocated(state%rt%vloc_new))then
      copy%rt%vloc_new=state%rt%vloc_new
    else
      if(allocated(copy%rt%vloc_new))deallocate(copy%rt%vloc_new)
    endif
    if(allocated(state%rt%vloc_old))then
      copy%rt%vloc_old=state%rt%vloc_old
    else
      if(allocated(copy%rt%vloc_old))deallocate(copy%rt%vloc_old)
    endif
    if(allocated(state%rt%rho0_s))then
      copy%rt%rho0_s=state%rt%rho0_s
    else
      if(allocated(copy%rt%rho0_s))deallocate(copy%rt%rho0_s)
    endif
    copy%rt%vonf=state%rt%vonf
    copy%rt%j_e=state%rt%j_e
    copy%rt%E_old=state%rt%E_old
    copy%rt%rho0_loaded=state%rt%rho0_loaded
    copy%accepted_A=state%accepted_A
    copy%fields=0d0
    copy%fields(:,1)=state%system%vec_Ac
    copy%fields(:,2)=state%system%vec_Ac_ext
    copy%fields(:,3)=state%system%vec_E
    copy%fields(:,4)=state%system%vec_E_ext
  endif
end subroutine
subroutine checkpoint_io(unit,reading,value,ios)
  integer,intent(in)::unit
  logical,intent(in)::reading
  type(checkpoint_payload),intent(inout)::value
  integer,intent(inout)::ios
  logical::present_array
  integer::lo(7),hi(7),i1,i2,k,allocation_status
  integer(8)::elements,extent
  if(ios/=0)return
  call checkpoint_io_noncollinear(unit,reading,value%noncollinear,ios,checkpoint_io_limit())
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
        if(extent<0_8.or.extent>state%checkpoint_array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>state%checkpoint_array_limit/extent)then
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
        if(extent<0_8.or.extent>state%checkpoint_array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>state%checkpoint_array_limit/extent)then
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
        if(extent<0_8.or.extent>state%checkpoint_array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>state%checkpoint_array_limit/extent)then
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
    read(unit,iostat=ios) value%worker_A
  else
    write(unit,iostat=ios) value%worker_A
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%worker_current
  else
    write(unit,iostat=ios) value%worker_current
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%accepted_A
  else
    write(unit,iostat=ios) value%accepted_A
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%Mit
  else
    write(unit,iostat=ios) value%Mit
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%itt
  else
    write(unit,iostat=ios) value%itt
  endif
  if(ios/=0)return
end subroutine
integer(8) function checkpoint_io_limit()
 checkpoint_io_limit=state%checkpoint_array_limit
end function
#include "noncollinear_disk_io.inc"
#include "checkpoint_arrays.inc"
end module
