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
   rt_checkpoint,rt_current,rt_accepted_step,rt_vector_potential
 use structures,only:s_rt_run_config
 implicit none
 private
 public::sp_fast_hse_init,sp_fast_hse_step,sp_fast_hse_checkpoint,sp_fast_hse_close,sp_fast_hse_info
 type(tddft_context),pointer,save::state=>null()
 logical,save::used=.false.,ready=.false.,failed=.false.
 real(8),save::accepted_A(3)=0d0
contains
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
 integer(c_int) function sp_fast_hse_init(path,timestep,capacity,J) bind(C)
 use inputoutput,only:read_input
 use math_constants,only:set_math_constants
 use plusU_global,only:read_Hubbard_parameters
 use sym_sub,only:read_sw_symmetry
 use timer,only:timer_initialize
 use salmon_global,only:theory,xc,calc_mode,yn_symmetry,nproc_k,nproc_ob,nproc_rgrid
#ifdef USE_MPI
 use exx_worker_service,only:service_enabled
#endif
 character(c_char),intent(in)::path(*)
 real(c_double),value::timestep
 integer(c_int),value::capacity
 real(c_double),intent(out)::J(3)
 type(s_rt_run_config)::config
 character(4096)::filename
 integer::status,saved_capacity
 real(8)::saved_dt
 logical::exists
 sp_fast_hse_init=1;J=0d0
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
 call rt_open(state,'native_fast_hse',status,config=config)
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
 ready=.true.;sp_fast_hse_init=0
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
  case(20:25,30,32,33)
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
end module
