! Copyright 2026 SALMONpy contributors. Apache-2.0.
! Explicitly attached, k-only fast-PR HSE context. One lifetime per process.
! Host owns MPI initialization/finalization and calls every API on the main
! thread collectively within the attached group. No world service is entered.
#include "config.h"
module salmonpy_fast_hse_mpi
 use filesystem,only:resolve_path
 use salmonpy_hse_model_inputs,only:s_hse_model_inputs
 use sym_sub,only:s_hse_symmetry_workspace,swap_hse_symmetry_workspace
 use sym_rho_sub,only:s_hse_density_workspace,swap_hse_density_workspace
 use iso_c_binding
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 use communication,only:comm_is_attached,comm_get_globalinfo,comm_get_max,comm_bcast
 use parallelization,only:setup_parallel,nproc_group_global,nproc_size_global
 use salmonpy_tddft_context,only:tddft_context,rt_open,rt_step,rt_close,rt_input_config, &
   rt_history_info,rt_step_transaction,rt_checkpoint,rt_current,rt_accepted_step,rt_vector_potential,rt_model_settings
 use salmonpy_gs_seed,only:s_gs_seed
 use salmonpy_gs_context,only:gs_seed_from_handle
 use structures,only:s_rt_run_config
 implicit none
 private
 public::sp_fast_hse_init_symmetry,sp_fast_hse_symmetry_protocol,sp_fast_hse_symmetry_kmpi_protocol
 public::sp_fast_hse_init,sp_fast_hse_step,sp_fast_hse_checkpoint,sp_fast_hse_close,sp_fast_hse_info
 type(tddft_context),pointer,save::state=>null()
 logical,save::used=.false.,ready=.false.,failed=.false.
 real(8),save::accepted_A(3)=0d0
 ! Registry owns numerical contexts; the legacy selector is installed only
 ! during one synchronous call. Registry assignment copies pointers, not owners.
 type handle_slot
  integer::model_index=0
  type(tddft_context),pointer::context=>null()
  integer(c_int64_t)::generation=0_c_int64_t
  real(8)::field(3)=0d0
  logical::broken=.false.
  character(80)::label=''
 end type
 type(handle_slot),allocatable,save::handles(:)
 integer(c_int64_t),parameter::handle_stride=2147483648_c_int64_t
 logical,save::handle_mode=.false.,model_symmetry=.false.
 integer,save::handle_group=0
 type model_slot
   type(s_gs_seed),allocatable::memory_gs
   type(s_hse_model_inputs),allocatable::values
   type(s_hse_symmetry_workspace)::symmetry
   type(s_hse_density_workspace)::density_symmetry
   integer(c_int64_t)::generation=0
   integer::children=0
   logical::closed=.false.,symmetric=.false.
 end type
 type(model_slot),allocatable,save::models(:)
 integer,save::active_model=0
 integer,save::gs_parked_model=-1
 logical,save::registry_mode=.false.


contains
integer(c_int) function sp_fast_hse_history_info(values) bind(C)
 integer(c_int),intent(out)::values(8)
 integer::status
 values=0;sp_fast_hse_history_info=1
 if(.not.attached_main_thread())return
 if(.not.ready.or.failed.or..not.associated(state))return
 call rt_history_info(state,values,status)
 sp_fast_hse_history_info=status
end function
integer(c_int) function sp_fast_hse_step_transaction_protocol() bind(C)
  sp_fast_hse_step_transaction_protocol=1
end function
integer(c_int) function sp_fast_hse_step_transaction(mode) bind(C)
  integer(c_int),value::mode
  integer::status
  sp_fast_hse_step_transaction=12
  if(.not.attached_main_thread())return
  status=0
  if(.not.ready.or.failed.or..not.associated(state))status=12
  call agree(status)
  if(status/=0)return
  call rt_step_transaction(state,mode,status)
  if(status==0.and.mode==1)call rt_vector_potential(state,accepted_A,status)
  call agree(status)
  sp_fast_hse_step_transaction=status
end function
integer(c_int) function sp_fast_hse_recovery_protocol() bind(C)
  sp_fast_hse_recovery_protocol=1
end function

integer(c_int) function sp_fast_hse_checkpoint_protocol() bind(C)
  sp_fast_hse_checkpoint_protocol=1
end function

 logical function attached_main_thread() result(ok)
#ifdef USE_MPI
 use mpi,only:MPI_Is_thread_main,MPI_SUCCESS
 integer::ierr
 logical::main_thread
 ok=.false.
 if(.not.comm_is_attached())return
 call MPI_Is_thread_main(main_thread,ierr)
 ok=ierr==MPI_SUCCESS.and.main_thread
#else
 ok=.false.
#endif
 end function
 subroutine agree(status)
 integer,intent(inout)::status
 integer::group,rank,count
 call comm_get_globalinfo(group,rank,count)
 call comm_get_max(status,group)
 end subroutine
 subroutine c_path(path,filename,status)
 character(c_char),intent(in)::path(*)
 character(*),intent(out)::filename
 integer,intent(out)::status
 integer::i
 filename='';status=1
 do i=1,len(filename)
  if(path(i)==c_null_char)exit
  filename(i:i)=path(i)
 enddo
 if(i>len(filename).or.len_trim(filename)==0)return
 status=0
 end subroutine
 integer(c_int) function sp_fast_hse_symmetry_kmpi_protocol() bind(C)
 sp_fast_hse_symmetry_kmpi_protocol=1
 end function
 integer(c_int) function sp_fast_hse_symmetry_protocol() bind(C)
 sp_fast_hse_symmetry_protocol=1
 end function
 integer(c_int) function sp_fast_hse_init(path,timestep,capacity,J) bind(C)
 character(c_char),intent(in)::path(*)
 real(c_double),value::timestep
 integer(c_int),value::capacity
 real(c_double),intent(out)::J(3)
 sp_fast_hse_init=initialize(path,timestep,capacity,J,.false.)
 end function
 integer(c_int) function sp_fast_hse_init_symmetry(path,timestep,capacity,J) bind(C)
 character(c_char),intent(in)::path(*)
 real(c_double),value::timestep
 integer(c_int),value::capacity
 real(c_double),intent(out)::J(3)
 sp_fast_hse_init_symmetry=initialize(path,timestep,capacity,J,.true.)
 end function
 integer(c_int) function initialize(path,timestep,capacity,J,allow_symmetry,registered,output_name,memory_gs)
 use inputoutput,only:read_input,s_input_parser
 use math_constants,only:set_math_constants
 use plusU_global,only:read_Hubbard_parameters
 use sym_sub,only:read_sw_symmetry
 use timer,only:timer_initialize
 use salmon_global,only:theory,xc,calc_mode,yn_symmetry,nproc_k,nproc_ob,nproc_rgrid
#ifdef USE_MPI
 use exx_worker_service,only:service_enabled
#endif
 type(s_gs_seed),optional,intent(in)::memory_gs
 logical,intent(in)::allow_symmetry
 character(*),optional,intent(in)::output_name
 character(80)::selected_name
 logical,optional,intent(in)::registered
 logical::registration,first_init
 character(c_char),intent(in)::path(*)
 real(c_double),value::timestep
 integer(c_int),value::capacity
 real(c_double),intent(out)::J(3)
 type(s_rt_run_config)::config
 type(s_input_parser)::parser
 character(4096)::filename
 integer::status,saved_capacity
 real(8)::saved_dt
 logical::exists
 initialize=1;J=0d0
 if(.not.attached_main_thread())return
 call c_path(path,filename,status)
 registration=.false.
 if(present(registered))registration=registered
 first_init=.not.used
 if((used.and..not.registration).or.capacity<1.or.capacity>=huge(0))status=1
 if(.not.ieee_is_finite(timestep).or.timestep<=0d0)status=1
 inquire(file=resolve_path(trim(filename)),exist=exists)
 if(.not.exists)status=1
 call agree(status)
 if(status/=0)return
 used=.true.
 if(first_init)then
  call set_math_constants
  call setup_parallel
 endif
 saved_dt=timestep;saved_capacity=capacity
 call comm_bcast(saved_dt,nproc_group_global)
 call comm_bcast(saved_capacity,nproc_group_global)
 if(saved_dt/=timestep.or.saved_capacity/=capacity)status=1
 call agree(status)
 if(status/=0)return
 if(registration)then
  call read_input(trim(filename),parser)
 else
  call read_input(trim(filename))
 endif
 if(theory/='tddft_response'.or.xc/='hse06')status=1
 if(allow_symmetry.and.(nproc_size_global<1.or.nproc_size_global>2))status=1
 if(nproc_k/=nproc_size_global.or.nproc_ob/=1.or.any(nproc_rgrid/=1))status=1
#ifdef USE_MPI
 if(service_enabled())status=1
#endif
 call agree(status)
 if(status/=0)return
 calc_mode='RT'
 call read_Hubbard_parameters
 call read_sw_symmetry(yn_symmetry)
 if(first_init)call timer_initialize
 call rt_input_config(config)
 config%delta_t=timestep;config%max_steps=capacity
 selected_name='native_fast_hse'
 if(present(output_name))selected_name=output_name
 call rt_open(state,trim(selected_name),status,config=config,allow_symmetry=allow_symmetry,memory_gs=memory_gs)
 call agree(status)
 if(status/=0)then
  failed=.true.;return
 endif
 call rt_current(state,J,status)
 call agree(status)
 if(status==0)call rt_vector_potential(state,accepted_A,status)
 call agree(status)
 if(status/=0)then
  call rt_close(state)
  failed=.true.;return
 endif
 ready=.true.;initialize=0
 end function
 integer(c_int) function sp_fast_hse_step(Aold,Anew,J,energy) bind(C)
 real(c_double),intent(in)::Aold(3),Anew(3)
 real(c_double),intent(out)::J(3),energy
 real(8)::reference(6)
 integer::status,i
 sp_fast_hse_step=1;J=0d0;energy=0d0
 if(.not.attached_main_thread())return
 status=0
 if(.not.ready.or.failed.or..not.associated(state))status=1
 do i=1,3
  if(.not.ieee_is_finite(Aold(i)).or..not.ieee_is_finite(Anew(i)))status=1
 enddo
 if(any(Aold/=accepted_A))status=1
 call agree(status)
 if(status/=0)return
 reference=[Aold,Anew]
 call comm_bcast(reference,nproc_group_global)
 if(any(reference/=[Aold,Anew]))status=1
 call agree(status)
 if(status/=0)return
 call rt_step(state,Anew,J,energy,status)
 call agree(status)
 if(status/=0)then
  select case(status)
  case(20:25,30,32,33,34)
    sp_fast_hse_step=status
  case default
    failed=.true.
  end select
  return
 endif
 accepted_A=Anew
 sp_fast_hse_step=0
 end function
 integer(c_int) function sp_fast_hse_checkpoint(path,restore) bind(C)
 character(c_char),intent(in)::path(*)
 integer(c_int),value::restore
 character(4096)::filename
 integer::status,mode
 sp_fast_hse_checkpoint=1
 if(.not.attached_main_thread())return
 status=0;filename=""
 if(restore<3)call c_path(path,filename,status)
 if(.not.ready.or.failed.or..not.associated(state))status=1
 if(restore<0.or.restore>4)status=1
 call agree(status)
 if(status/=0)return
 mode=restore
 call comm_bcast(mode,nproc_group_global)
 if(mode/=restore)status=1
 call agree(status)
 if(status/=0)return
 ! Caller passes distinct per-rank paths, authenticated together outside native.
 if(restore>=2)then
  call rt_checkpoint(state,trim(filename),.true.,status,stage=restore)
 else
  call rt_checkpoint(state,trim(filename),restore==1,status)
 endif
 call agree(status)
 if(status/=0)then
  if(restore==1)failed=.true.
  return
 endif
 if(restore==1.or.restore==3)then
  call rt_vector_potential(state,accepted_A,status)
  call agree(status)
  if(status/=0)then
   failed=.true.;return
  endif
 endif
 sp_fast_hse_checkpoint=0
 end function
 integer(c_int) function sp_fast_hse_info(step,A,J) bind(C)
 integer(c_int),intent(out)::step
 real(c_double),intent(out)::A(3),J(3)
 integer::status
 sp_fast_hse_info=1;step=-1;A=0d0;J=0d0
 if(.not.attached_main_thread())return
 status=0
 if(.not.ready.or.failed.or..not.associated(state))status=1
 call agree(status)
 if(status/=0)return
 call rt_vector_potential(state,A,status)
 call agree(status)
 if(status==0)call rt_current(state,J,status)
 call agree(status)
 if(status/=0)return
 step=rt_accepted_step(state)
 if(step<0)status=1
 call agree(status)
 if(status/=0)return
 sp_fast_hse_info=0
 end function
 integer(c_int) function sp_fast_hse_close() bind(C)
 integer::status
 sp_fast_hse_close=1
 if(.not.attached_main_thread())return
 status=0
 if(.not.associated(state))status=1
 call agree(status)
 if(status/=0)return
 ! Output/native resources close before the host releases its attachment.
 ! Never call end_parallel here: MPI belongs to the host application.
 call rt_close(state)
 ready=.false.
 sp_fast_hse_close=0
 end function
! Versioned, synchronous handle entry. Numerical operations retain the legacy
! implementation; only the selected pointer/accepted field/failure flags vary.
integer(c_int) function sp_fast_hse_handle_protocol() bind(C)
 sp_fast_hse_handle_protocol=1
end function
subroutine handle_lookup(token,index,status)
 integer(c_int64_t),intent(in)::token
 integer,intent(out)::index,status
 integer::group,rank,count,parts(2),reference(2)
 index=0;status=12
 if(.not.attached_main_thread())return
 call comm_get_globalinfo(group,rank,count)
 status=0
 if(.not.handle_mode.or.group/=handle_group.or.token<=0_c_int64_t)status=12
 if(token/handle_stride>int(huge(0),c_int64_t))status=12
 call agree(status)
 if(status/=0)return
 parts=[int(mod(token,handle_stride)),int(token/handle_stride)]
 reference=parts;call comm_bcast(reference,group)
 if(any(reference/=parts))status=11
 call agree(status)
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
 call agree(status)
end subroutine
subroutine handle_reserve(index,token,status)
 integer,intent(out)::index,status
 integer(c_int64_t),intent(out)::token
 type(handle_slot),allocatable::grown(:)
 integer::i,old,count,allocation_status,group,rank,n,reference(2),parts(2)
 status=0;index=0;token=0;old=0
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
  call agree(status)
  if(status/=0)return
  count=max(4,2*old)
  allocate(grown(count),stat=allocation_status)
  if(allocation_status/=0)status=16
  call agree(status)
  if(status/=0)return
  if(old>0)grown(:old)=handles
  call move_alloc(grown,handles)
  index=old+1
 endif
 parts=[index,int(handles(index)%generation+1_c_int64_t)]
 reference=parts
 call comm_get_globalinfo(group,rank,n)
 call comm_bcast(reference,group)
 if(any(reference/=parts))status=11
 call agree(status)
 if(status/=0)return
 handles(index)%generation=handles(index)%generation+1_c_int64_t
 token=handles(index)%generation*handle_stride+int(index,c_int64_t)
 handles(index)%field=0d0;handles(index)%broken=.false.
end subroutine
subroutine handle_select(token,index,status)
 integer(c_int64_t),intent(in)::token
 integer,intent(out)::index,status
 call handle_lookup(token,index,status)
 if(status/=0)return
 if(handles(index)%model_index>0)call select_model(handles(index)%model_index)
 state=>handles(index)%context
 accepted_A=handles(index)%field;failed=handles(index)%broken;ready=.true.
end subroutine
subroutine handle_unselect(index)
 integer,intent(in)::index
 handles(index)%context=>state
 handles(index)%field=accepted_A;handles(index)%broken=failed
 nullify(state);accepted_A=0d0;failed=.false.;ready=.false.
end subroutine
integer(c_int) function sp_fast_hse_handle_init(path,timestep,capacity,symmetry,token,J) bind(C)
 character(c_char),intent(in)::path(*)
 real(c_double),value::timestep
 integer(c_int),value::capacity,symmetry
 integer(c_int64_t),intent(out)::token
 real(c_double),intent(out)::J(3)
 integer::status,index,group,rank,n,reference_symmetry,reference_capacity
 real(8)::reference_dt
 integer(c_int64_t)::candidate
 character(4096)::filename,reference_path
 token=0;J=0d0;sp_fast_hse_handle_init=12
 if(.not.attached_main_thread())return
 call c_path(path,filename,status)
 if(used.or.handle_mode.or.symmetry<0.or.symmetry>1)status=12
 if(capacity<1.or.capacity>=huge(0))status=1
 if(.not.ieee_is_finite(timestep).or.timestep<=0d0)status=1
 call agree(status)
 if(status/=0)then
  sp_fast_hse_handle_init=status;return
 endif
 call comm_get_globalinfo(group,rank,n)
 reference_path=filename;reference_symmetry=symmetry
 reference_dt=timestep;reference_capacity=capacity
 call comm_bcast(reference_path,group);call comm_bcast(reference_symmetry,group)
 call comm_bcast(reference_dt,group);call comm_bcast(reference_capacity,group)
 if(reference_path/=filename.or.reference_symmetry/=symmetry)status=11
 if(reference_dt/=timestep.or.reference_capacity/=capacity)status=11
 call agree(status)
 if(status/=0)then
  sp_fast_hse_handle_init=status;return
 endif
 call handle_reserve(index,candidate,status)
 if(status/=0)then
  sp_fast_hse_handle_init=status;return
 endif
 status=initialize(path,timestep,capacity,J,symmetry==1)
 if(status==0)then
  handles(index)%label='native_fast_hse'
  handle_mode=.true.;handle_group=group;model_symmetry=symmetry==1
  call handle_unselect(index)
  token=candidate
 endif
 sp_fast_hse_handle_init=status
end function
integer(c_int) function sp_fast_hse_context_mpi_protocol() bind(C)
 sp_fast_hse_context_mpi_protocol=1
end function
integer(c_int) function sp_fast_hse_mpi_resources(values) bind(C)
 use communication,only:comm_child_counts
 integer(c_int),intent(out)::values(2)
 values=0;sp_fast_hse_mpi_resources=12
 if(.not.attached_main_thread())return
 call comm_child_counts(values(1),values(2))
 sp_fast_hse_mpi_resources=0
end function
integer(c_int) function sp_fast_hse_shared_symmetry_protocol() bind(C)
 sp_fast_hse_shared_symmetry_protocol=1
end function
integer(c_int) function sp_fast_hse_model_protocol() bind(C)
 sp_fast_hse_model_protocol=1
end function
integer(c_int) function sp_fast_hse_handle_new(model,timestep,capacity,token,J) bind(C)
 integer(c_int64_t),value::model
 real(c_double),value::timestep
 integer(c_int),value::capacity
 integer(c_int64_t),intent(out)::token
 real(c_double),intent(out)::J(3)
 integer::parent,status
 token=0;J=0d0
 call handle_lookup(model,parent,status)
 sp_fast_hse_handle_new=status
 if(status/=0)return
 if(handles(parent)%model_index>0)then
  if(models(handles(parent)%model_index)%closed)then
   sp_fast_hse_handle_new=12;return
  endif
  call select_model(handles(parent)%model_index)
 endif
 sp_fast_hse_handle_new=create_model_state(timestep,capacity,token,J)
end function
integer(c_int) function sp_fast_hse_model_new(timestep,capacity,token,J) bind(C)
 real(c_double),value::timestep
 integer(c_int),value::capacity
 integer(c_int64_t),intent(out)::token
 real(c_double),intent(out)::J(3)
 integer::group,rank,n,status
 token=0;J=0d0;sp_fast_hse_model_new=12
 if(.not.attached_main_thread())return
 call comm_get_globalinfo(group,rank,n)
 status=0
 if(.not.handle_mode.or.group/=handle_group.or.registry_mode)status=12
 call agree(status)
 if(status/=0)return
 sp_fast_hse_model_new=create_model_state(timestep,capacity,token,J)
end function
integer(c_int) function create_model_state(timestep,capacity,token,J)
 real(c_double),intent(in)::timestep
 integer(c_int),intent(in)::capacity
 integer(c_int64_t),intent(out)::token
 real(c_double),intent(out)::J(3)
 type(tddft_context),pointer::fresh
 type(s_rt_run_config)::config
 integer::status,index,reference_capacity
 real(8)::reference_dt
 integer(c_int64_t)::candidate
 token=0;J=0d0;nullify(fresh)
 status=0
 if(capacity<1.or.capacity>=huge(0))status=1
 if(.not.ieee_is_finite(timestep).or.timestep<=0d0)status=1
 reference_dt=timestep;reference_capacity=capacity
 call comm_bcast(reference_dt,nproc_group_global)
 call comm_bcast(reference_capacity,nproc_group_global)
 if(reference_dt/=timestep.or.reference_capacity/=capacity)status=11
 call agree(status)
 if(status/=0)then
  create_model_state=status;return
 endif
 call handle_reserve(index,candidate,status)
 if(status/=0)then
  create_model_state=status;return
 endif
 write(handles(index)%label,'(a,i0,a,i0)')'native_hse_h',index,'_g',handles(index)%generation
 call rt_input_config(config);config%delta_t=timestep;config%max_steps=capacity
 if(active_model>0)then
  if(allocated(models(active_model)%memory_gs))then
   call rt_open(fresh,trim(handles(index)%label),status,config=config, &
     allow_symmetry=model_symmetry,memory_gs=models(active_model)%memory_gs)
  else
   call rt_open(fresh,trim(handles(index)%label),status,config=config,allow_symmetry=model_symmetry)
  endif
 else
  call rt_open(fresh,trim(handles(index)%label),status,config=config,allow_symmetry=model_symmetry)
 endif
 call agree(status)
 if(status==0)call rt_current(fresh,J,status)
 call agree(status)
 if(status==0)call rt_vector_potential(fresh,handles(index)%field,status)
 call agree(status)
 if(status==0)then
  handles(index)%context=>fresh;token=candidate
  handles(index)%model_index=active_model
  if(active_model>0)models(active_model)%children=models(active_model)%children+1
 else if(associated(fresh))then
  call rt_close(fresh)
 endif
 create_model_state=status
end function
integer(c_int) function sp_fast_hse_handle_step(token,Aold,Anew,J,energy) bind(C)
 integer(c_int64_t),value::token
 real(c_double),intent(in)::Aold(3),Anew(3)
 real(c_double),intent(out)::J(3),energy
 integer::index,status
 J=0d0;energy=0d0
 call handle_select(token,index,status)
 sp_fast_hse_handle_step=status
 if(status/=0)return
 sp_fast_hse_handle_step=sp_fast_hse_step(Aold,Anew,J,energy)
 call handle_unselect(index)
end function
integer(c_int) function sp_fast_hse_handle_info(token,step,A,J) bind(C)
 integer(c_int64_t),value::token
 integer(c_int),intent(out)::step
 real(c_double),intent(out)::A(3),J(3)
 integer::index,status
 step=-1;A=0d0;J=0d0
 call handle_select(token,index,status)
 sp_fast_hse_handle_info=status
 if(status/=0)return
 sp_fast_hse_handle_info=sp_fast_hse_info(step,A,J)
 call handle_unselect(index)
end function
integer(c_int) function sp_fast_hse_handle_history(token,values) bind(C)
 integer(c_int64_t),value::token
 integer(c_int),intent(out)::values(8)
 integer::index,status
 values=0
 call handle_select(token,index,status)
 sp_fast_hse_handle_history=status
 if(status/=0)return
 sp_fast_hse_handle_history=sp_fast_hse_history_info(values)
 call handle_unselect(index)
end function
integer(c_int) function sp_fast_hse_handle_transaction(token,mode) bind(C)
 integer(c_int64_t),value::token
 integer(c_int),value::mode
 integer::index,status
 call handle_select(token,index,status)
 sp_fast_hse_handle_transaction=status
 if(status/=0)return
 sp_fast_hse_handle_transaction=sp_fast_hse_step_transaction(mode)
 call handle_unselect(index)
end function
integer(c_int) function sp_fast_hse_handle_checkpoint(token,path,mode) bind(C)
 integer(c_int64_t),value::token
 character(c_char),intent(in)::path(*)
 integer(c_int),value::mode
 integer::index,status
 call handle_select(token,index,status)
 sp_fast_hse_handle_checkpoint=status
 if(status/=0)return
 sp_fast_hse_handle_checkpoint=sp_fast_hse_checkpoint(path,mode)
 call handle_unselect(index)
end function
integer(c_int) function sp_fast_hse_handle_label(token,buffer,capacity) bind(C)
 integer(c_int64_t),value::token
 integer(c_int),value::capacity
 character(c_char),intent(out)::buffer(capacity)
 integer::index,status,i,n
 call handle_lookup(token,index,status)
 sp_fast_hse_handle_label=status
 if(status/=0)return
 n=len_trim(handles(index)%label)
 if(capacity<=n)status=13
 call agree(status)
 sp_fast_hse_handle_label=status
 if(status/=0)return
 do i=1,n
  buffer(i)=handles(index)%label(i:i)
 enddo
 buffer(n+1)=c_null_char
end function
integer(c_int) function sp_fast_hse_handle_close(token) bind(C)
 integer(c_int64_t),value::token
 integer::index,status
 call handle_select(token,index,status)
 sp_fast_hse_handle_close=status
 if(status/=0)return
 sp_fast_hse_handle_close=sp_fast_hse_close()
 call handle_unselect(index)
 if(sp_fast_hse_handle_close==0.and.handles(index)%model_index>0)then
  block
  integer::owner
  owner=handles(index)%model_index
  models(owner)%children=models(owner)%children-1
  handles(index)%model_index=0
  call release_model(owner)
  end block
 endif
end function

subroutine park_model()
 integer::old
 if(active_model==0)return
 old=active_model
 call models(old)%values%capture()
 call swap_hse_symmetry_workspace(models(old)%symmetry)
 call swap_hse_density_workspace(models(old)%density_symmetry)
 active_model=0
end subroutine
subroutine select_model(index)
 integer,intent(in)::index
 if(index==active_model)return
 call park_model()
 call models(index)%values%activate()
 call swap_hse_symmetry_workspace(models(index)%symmetry)
 call swap_hse_density_workspace(models(index)%density_symmetry)
 model_symmetry=models(index)%symmetric
 call rt_model_settings(.false.)
 active_model=index
end subroutine
subroutine release_model(index)
 integer,intent(in)::index
 type(s_hse_symmetry_workspace)::empty
 type(s_hse_density_workspace)::empty_density
 if(.not.models(index)%closed.or.models(index)%children/=0)return
 if(index==active_model)call park_model()
 if(allocated(models(index)%values))deallocate(models(index)%values)
 if(allocated(models(index)%memory_gs))deallocate(models(index)%memory_gs)
 models(index)%symmetry=empty
 models(index)%density_symmetry=empty_density
end subroutine
subroutine model_lookup(token,index,status)
 integer(c_int64_t),intent(in)::token
 integer,intent(out)::index,status
 integer(c_int64_t)::number,generation
 integer::parts(2),reference(2)
 parts=transfer(token,parts);reference=parts
 call comm_bcast(reference,nproc_group_global)
 status=12;index=0
 if(any(reference/=parts))then
  status=11
 else if(token<0.and.token/=-huge(token)-1.and.allocated(models))then
  number=modulo(-token,handle_stride);generation=(-token)/handle_stride
  if(number>=1.and.number<=size(models))then
   index=int(number)
   if(allocated(models(index)%values).and.models(index)%generation==generation)status=0
  endif
 endif
 call agree(status)
 if(status/=0)index=0
end subroutine
integer(c_int) function sp_fast_hse_registry_protocol() bind(C)
 sp_fast_hse_registry_protocol=1
end function
integer(c_int) function sp_fast_hse_memory_gs_protocol() bind(C)
 sp_fast_hse_memory_gs_protocol=1
end function
integer(c_int) function sp_fast_hse_gs_boundary(mode) bind(C)
 integer(c_int),value::mode
 integer::status,reference
 sp_fast_hse_gs_boundary=12
 if(.not.attached_main_thread())return
 status=0;reference=mode
 call comm_bcast(reference,nproc_group_global)
 if(reference/=mode)status=11
 if(mode==0)then
  if(gs_parked_model/=-1.or.associated(state))status=12
 else if(mode==1)then
  if(gs_parked_model<0)status=12
 else
  status=11
 endif
 call agree(status)
 if(status/=0)then
  sp_fast_hse_gs_boundary=status;return
 endif
 if(mode==0)then
  gs_parked_model=active_model
  call park_model()
 else
  if(gs_parked_model>0)call select_model(gs_parked_model)
  gs_parked_model=-1
 endif
 sp_fast_hse_gs_boundary=0
end function
integer(c_int) function sp_fast_hse_model_register(path,timestep,capacity,symmetry,model_token,state_token,J) bind(C)
 character(c_char),intent(in)::path(*)
 real(c_double),value::timestep
 integer(c_int),value::capacity,symmetry
 integer(c_int64_t),intent(out)::model_token,state_token
 real(c_double),intent(out)::J(3)
 sp_fast_hse_model_register=register_model(path,timestep,capacity,symmetry,model_token,state_token,J)
end function
integer(c_int) function sp_fast_hse_model_from_gs(gs_token,path,timestep,capacity,symmetry,model_token,state_token,J) bind(C)
 integer(c_int64_t),value::gs_token
 character(c_char),intent(in)::path(*)
 real(c_double),value::timestep
 integer(c_int),value::capacity,symmetry
 integer(c_int64_t),intent(out)::model_token,state_token
 real(c_double),intent(out)::J(3)
 type(s_gs_seed)::seed
 integer::status
 model_token=0;state_token=0;J=0d0
 call gs_seed_from_handle(gs_token,seed,status)
 if(status==0)then
  if(any(seed%decomposition/=[nproc_size_global,1,1,1,1]))status=11
 endif
 call agree(status)
 sp_fast_hse_model_from_gs=status
 if(status/=0)return
 sp_fast_hse_model_from_gs=register_model(path,timestep,capacity,symmetry,model_token,state_token,J,seed)
end function
integer(c_int) function register_model(path,timestep,capacity,symmetry,model_token,state_token,J,memory_gs)
 character(c_char),intent(in)::path(*)
 real(c_double),value::timestep
 integer(c_int),value::capacity,symmetry
 integer(c_int64_t),intent(out)::model_token,state_token
 real(c_double),intent(out)::J(3)
 type(s_gs_seed),optional,intent(in)::memory_gs
 type(model_slot),allocatable::grown(:)
 character(4096)::filename,reference_path
 integer::status,index,slot,i,group,rank,n,previous,ref_int(2)
 real(8)::ref_dt
 integer(c_int64_t)::candidate
 logical::exists
 model_token=0;state_token=0;J=0d0;register_model=12
 if(.not.attached_main_thread())return
 call comm_get_globalinfo(group,rank,n)
 call c_path(path,filename,status)
 if(used.and..not.registry_mode)status=12
 if(capacity<1.or.capacity>=huge(0).or.symmetry<0.or.symmetry>1)status=1
 if(.not.ieee_is_finite(timestep).or.timestep<=0d0)status=1
 inquire(file=resolve_path(trim(filename)),exist=exists)
 if(.not.exists)status=1
 reference_path=filename;ref_int=[capacity,symmetry];ref_dt=timestep
 call comm_bcast(reference_path,group);call comm_bcast(ref_int,group);call comm_bcast(ref_dt,group)
 if(reference_path/=filename.or.any(ref_int/=[capacity,symmetry]).or.ref_dt/=timestep)status=11
 call agree(status)
 register_model=status
 if(status/=0)return
 if(.not.allocated(models))allocate(models(4))
 index=0
 do i=1,size(models)
  if(allocated(models(i)%values))cycle
  if(models(i)%generation>=int(huge(0),c_int64_t))cycle
  index=i;exit
 enddo
 if(index==0)then
  index=size(models)+1
  allocate(grown(size(models)*2));grown(:size(models))=models
  call move_alloc(grown,models)
 endif
 previous=active_model
 call park_model()
 allocate(models(index)%values)
 models(index)%generation=models(index)%generation+1
 if(present(memory_gs))then
  allocate(models(index)%memory_gs)
  models(index)%memory_gs=memory_gs
 endif
 models(index)%children=0;models(index)%closed=.false.;models(index)%symmetric=symmetry==1
 call rt_model_settings(.true.)
 call handle_reserve(slot,candidate,status)
 if(status==0)then
  write(handles(slot)%label,'(a,i0,a,i0)')'native_hse_m',index,'_g',models(index)%generation
  status=initialize(path,timestep,capacity,J,symmetry==1,.true.,trim(handles(slot)%label),memory_gs)
 endif
 call agree(status)
 if(status==0)then
  call models(index)%values%capture()
  active_model=index;model_symmetry=symmetry==1
  handles(slot)%model_index=index
  models(index)%children=1
  handle_mode=.true.;registry_mode=.true.;handle_group=group
  call handle_unselect(slot)
  model_token=-(models(index)%generation*handle_stride+int(index,c_int64_t))
  state_token=candidate
 else
  if(associated(state))call rt_close(state)
  nullify(state);ready=.false.;failed=.false.
  active_model=index;models(index)%closed=.true.
  call release_model(index)
  if(previous>0)call select_model(previous)
 endif
 register_model=status
end function
integer(c_int) function sp_fast_hse_model_spawn(model,timestep,capacity,token,J) bind(C)
 integer(c_int64_t),value::model
 real(c_double),value::timestep
 integer(c_int),value::capacity
 integer(c_int64_t),intent(out)::token
 real(c_double),intent(out)::J(3)
 integer::index,status
 token=0;J=0d0;sp_fast_hse_model_spawn=12
 if(.not.attached_main_thread())return
 call model_lookup(model,index,status)
 if(status==0)then
  if(models(index)%closed)status=12
 endif
 call agree(status)
 sp_fast_hse_model_spawn=status
 if(status/=0)return
 call select_model(index)
 sp_fast_hse_model_spawn=create_model_state(timestep,capacity,token,J)
end function
integer(c_int) function sp_fast_hse_model_close(model) bind(C)
 integer(c_int64_t),value::model
 integer::index,status
 sp_fast_hse_model_close=12
 if(.not.attached_main_thread())return
 call model_lookup(model,index,status)
 sp_fast_hse_model_close=status
 if(status/=0)return
 models(index)%closed=.true.
 call release_model(index)
end function

end module
