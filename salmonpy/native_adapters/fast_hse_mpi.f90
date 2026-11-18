! Copyright 2026 SALMONpy contributors. Apache-2.0.
! Explicitly attached, k-only fast-PR HSE context. One lifetime per process.
! Host owns MPI initialization/finalization and calls every API on the main
! thread collectively within the attached group. No world service is entered.
#include "config.h"
module salmonpy_fast_hse_mpi
 use iso_c_binding
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 use communication,only:comm_is_attached,comm_get_globalinfo,comm_get_max,comm_bcast
 use parallelization,only:setup_parallel,nproc_group_global,nproc_size_global
 use salmonpy_tddft_context,only:tddft_context,rt_open,rt_step,rt_close,rt_input_config, &
   rt_history_info,rt_step_transaction,rt_checkpoint,rt_current,rt_accepted_step,rt_vector_potential
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
  type(tddft_context),pointer::context=>null()
  integer(c_int64_t)::generation=0_c_int64_t
  real(8)::field(3)=0d0
  logical::broken=.false.
  character(80)::label=''
 end type
 type(handle_slot),allocatable,save::handles(:)
 integer(c_int64_t),parameter::handle_stride=2147483648_c_int64_t
 logical,save::handle_mode=.false.
 integer,save::handle_group=0

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
 integer(c_int) function initialize(path,timestep,capacity,J,allow_symmetry)
 use inputoutput,only:read_input
 use math_constants,only:set_math_constants
 use plusU_global,only:read_Hubbard_parameters
 use sym_sub,only:read_sw_symmetry
 use timer,only:timer_initialize
 use salmon_global,only:theory,xc,calc_mode,yn_symmetry,nproc_k,nproc_ob,nproc_rgrid
#ifdef USE_MPI
 use exx_worker_service,only:service_enabled
#endif
 logical,intent(in)::allow_symmetry
 character(c_char),intent(in)::path(*)
 real(c_double),value::timestep
 integer(c_int),value::capacity
 real(c_double),intent(out)::J(3)
 type(s_rt_run_config)::config
 character(4096)::filename
 integer::status,saved_capacity
 real(8)::saved_dt
 logical::exists
 initialize=1;J=0d0
 if(.not.attached_main_thread())return
 call c_path(path,filename,status)
 if(used.or.capacity<1.or.capacity>=huge(0))status=1
 if(.not.ieee_is_finite(timestep).or.timestep<=0d0)status=1
 inquire(file=trim(filename),exist=exists)
 if(.not.exists)status=1
 call agree(status)
 if(status/=0)return
 used=.true.
 call set_math_constants
 call setup_parallel
 saved_dt=timestep;saved_capacity=capacity
 call comm_bcast(saved_dt,nproc_group_global)
 call comm_bcast(saved_capacity,nproc_group_global)
 if(saved_dt/=timestep.or.saved_capacity/=capacity)status=1
 call agree(status)
 if(status/=0)return
 call read_input(trim(filename))
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
 call timer_initialize
 call rt_input_config(config)
 config%delta_t=timestep;config%max_steps=capacity
 call rt_open(state,'native_fast_hse',status,config=config,allow_symmetry=allow_symmetry)
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
 integer::status
 sp_fast_hse_step=1;J=0d0;energy=0d0
 if(.not.attached_main_thread())return
 status=0
 if(.not.ready.or.failed.or..not.associated(state))status=1
 if(.not.all(ieee_is_finite(Aold)).or..not.all(ieee_is_finite(Anew)))status=1
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
  handle_mode=.true.;handle_group=group
  call handle_unselect(index)
  token=candidate
 endif
 sp_fast_hse_handle_init=status
end function
integer(c_int) function sp_fast_hse_handle_new(model,timestep,capacity,token,J) bind(C)
 integer(c_int64_t),value::model
 real(c_double),value::timestep
 integer(c_int),value::capacity
 integer(c_int64_t),intent(out)::token
 real(c_double),intent(out)::J(3)
 type(tddft_context),pointer::fresh
 type(s_rt_run_config)::config
 integer::status,parent,index,reference_capacity
 real(8)::reference_dt
 integer(c_int64_t)::candidate
 token=0;J=0d0;nullify(fresh)
 call handle_lookup(model,parent,status)
 sp_fast_hse_handle_new=status
 if(status/=0)return
 if(capacity<1.or.capacity>=huge(0))status=1
 if(.not.ieee_is_finite(timestep).or.timestep<=0d0)status=1
 reference_dt=timestep;reference_capacity=capacity
 call comm_bcast(reference_dt,nproc_group_global)
 call comm_bcast(reference_capacity,nproc_group_global)
 if(reference_dt/=timestep.or.reference_capacity/=capacity)status=11
 call agree(status)
 if(status/=0)then
  sp_fast_hse_handle_new=status;return
 endif
 call handle_reserve(index,candidate,status)
 if(status/=0)then
  sp_fast_hse_handle_new=status;return
 endif
 write(handles(index)%label,'(a,i0,a,i0)')'native_hse_h',index,'_g',handles(index)%generation
 call rt_input_config(config);config%delta_t=timestep;config%max_steps=capacity
 call rt_open(fresh,trim(handles(index)%label),status,config=config)
 call agree(status)
 if(status==0)call rt_current(fresh,J,status)
 call agree(status)
 if(status==0)call rt_vector_potential(fresh,handles(index)%field,status)
 call agree(status)
 if(status==0)then
  handles(index)%context=>fresh;token=candidate
 else if(associated(fresh))then
  call rt_close(fresh)
 endif
 sp_fast_hse_handle_new=status
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
end function

end module
