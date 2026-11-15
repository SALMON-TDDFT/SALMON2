program probe
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 use salmon_global
 use parallelization
 use communication,only:comm_is_root
 use inputoutput
 use math_constants
 use timer
 use plusU_global,only:read_Hubbard_parameters
 use sym_sub,only:read_sw_symmetry
 use structures,only:s_rt_run_config
 use salmonpy_tddft_context,only:rt_input_config,tddft_context,rt_open,rt_step,rt_close,rt_probe_write
 implicit none
 type(tddft_context),pointer::a=>null(),b=>null()
 type(s_rt_run_config)::config_a,config_b
 character(16)::capacity_mode
 integer::capacity_a,capacity_b,saved_nt
 character(16)::mode,dt_mode,omp_mode,mode_a,mode_b
 character(16)::saved_stencil,saved_current,saved_force
 real(8)::ac(3),j(3),energy,saved_dt,tt,dt_a,dt_b,local_dt
 integer::step,code,ua,ub,local_step,api_unit
 character(32)::api_name
 call set_math_constants
 call setup_parallel
 call read_input
 if(theory/='tddft_response')stop 99
 calc_mode='RT'
 call read_Hubbard_parameters
 call read_sw_symmetry(yn_symmetry)
 call timer_initialize
 call get_command_argument(1,mode)
 call get_command_argument(2,dt_mode)
 call get_command_argument(3,omp_mode)
 mode_a='auto';mode_b='auto'
 if(omp_mode=='different_omp'.and.mode/='native')then
  mode_a='orbital';mode_b='rgrid'
 endif
 dt_a=dt;dt_b=dt
 if(dt_mode=='different_dt'.and.mode/='solo_b')dt_b=1.5d0*dt
 call get_command_argument(4,capacity_mode)
 call rt_input_config(config_a);call rt_input_config(config_b)
 config_a%delta_t=dt_a;config_b%delta_t=dt_b
 config_a%stencil_mode=mode_a;config_a%current_mode=mode_a;config_a%force_mode=mode_a
 config_b%stencil_mode=mode_b;config_b%current_mode=mode_b;config_b%force_mode=mode_b
 capacity_a=nt;capacity_b=nt
 if(capacity_mode=='different_cap'.and.mode/='native')capacity_b=6
 config_a%max_steps=capacity_a;config_b%max_steps=capacity_b
 if(mode/='solo_b')then
  config_a%max_steps=0
  call rt_open(a,'invalid_config',code,config=config_a);if(code==0.or.associated(a))stop 12
  config_a%max_steps=huge(0)
  call rt_open(a,'invalid_config',code,config=config_a);if(code==0.or.associated(a))stop 13
  config_a%max_steps=capacity_a
  call rt_open(a,'mixed_config',code,config=config_a,delta_t=dt_a)
  if(code==0.or.associated(a))stop 14
  call rt_open(a,'invalid_omp',code,openmp_mode='invalid');if(code==0.or.associated(a))stop 10
  call rt_open(a,'invalid_dt',code,delta_t=-1d0);if(code==0.or.associated(a))stop 8
  call rt_open(a,'invalid_dt',code,delta_t=ieee_value(0d0,ieee_quiet_nan))
  if(code==0.or.associated(a))stop 9
  call rt_open(a,'a',code,config=config_a);if(code/=0)stop 1
  ! Another legacy calculation entering warmup must not affect A, B or reopened RT states.
  if(mode/='solo_a'.and.mode/='native')exx_pre_scf_active=.true.
  open(newunit=ua,file='a_trace.bin',access='stream',form='unformatted',status='new')
 endif
 if(mode/='solo_a'.and.mode/='native')then
  call rt_open(b,'b',code,config=config_b);if(code/=0)stop 2
  open(newunit=ub,file='b_trace.bin',access='stream',form='unformatted',status='new')
 endif
 do step=1,nt
  if(step==nt/2+1.and.mode=='reopen_a')then
   call rt_close(a);close(ua)
   call rt_open(a,'a_reopened',code,config=config_a);if(code/=0)stop 4
   open(newunit=ua,file='a_reopened_trace.bin',access='stream',form='unformatted',status='new')
  endif
  if(step==nt/2+1.and.mode=='reopen_b')then
   call rt_close(b);close(ub)
   call rt_open(b,'b_reopened',code,config=config_b);if(code/=0)stop 5
   open(newunit=ub,file='b_reopened_trace.bin',access='stream',form='unformatted',status='new')
  endif
  if(mode=='reverse')then
   call advance(b,2,ub);call advance(a,1,ua)
  else
   if(associated(a))call advance(a,1,ua)
   if(associated(b))call advance(b,2,ub)
  endif
 enddo
 if(associated(a))then
  call rt_close(a);close(ua)
 endif
 if(associated(b))then
  call rt_close(b);close(ub)
 endif
 call end_parallel
 print *,'Full TDDFT context probe PASS'
contains
 subroutine advance(state,label,unit)
 type(tddft_context),pointer,intent(inout)::state
 integer,intent(in)::label,unit
 local_dt=dt_a
 if(label==2)local_dt=dt_b
 local_step=step
 if(step>nt/2.and.((mode=='reopen_a'.and.label==1).or.(mode=='reopen_b'.and.label==2)))local_step=step-nt/2
 if((label==1.and.local_step>capacity_a).or.(label==2.and.local_step>capacity_b))then
  ac=0d0
  call rt_step(state,ac,j,energy,code);if(code==0)stop 11
  return
 endif
 if(step==3.and.mode=='interleaved')then
  ac=ieee_value(0d0,ieee_quiet_nan)
  call rt_step(state,ac,j,energy,code);if(code==0)stop 6
  saved_dt=hse_omega;hse_omega=2d0*saved_dt
  ac=0d0
  call rt_step(state,ac,j,energy,code);if(code/=2)stop 7
  hse_omega=saved_dt
 endif
 ac=[.0001d0*sin(.3d0*local_step*label),.00007d0*sin(.2d0*local_step*(3-label)),.00005d0*label*sin(.4d0*local_step)]
 if(mode=='native')then
  tt=local_step*dt-.012d0
  ac=0d0
  if(abs(tt)<.012d0)ac(3)=-.01d0/20d0*(cos(pi*tt/.024d0))**2*sin(20d0*tt)
 endif
 ! Deliberately poison the old global after initialization, without selecting a state.
 saved_stencil=stencil_openmp_mode;saved_current=current_openmp_mode;saved_force=force_openmp_mode
 if(omp_mode=='different_omp'.and.mode/='native'.and.mode/='solo_a'.and.mode/='solo_b')then
  stencil_openmp_mode='rgrid';current_openmp_mode='rgrid';force_openmp_mode='rgrid'
 endif
 saved_dt=dt
 saved_nt=nt
 config_a%max_steps=1;config_b%max_steps=1
 if(capacity_mode=="different_cap")nt=1
 if(dt_mode=='different_dt'.and.mode/='native'.and.mode/='solo_a'.and.mode/='solo_b')dt=7d0*dt_a
 call rt_step(state,ac,j,energy,code)
 dt=saved_dt
 nt=saved_nt
 config_a%max_steps=capacity_a;config_b%max_steps=capacity_b
 stencil_openmp_mode=saved_stencil;current_openmp_mode=saved_current;force_openmp_mode=saved_force
 if(code/=0)stop 3
 api_name='a_api.dat'
 if(label==2)api_name='b_api.dat'
 if(step>nt/2.and.mode=='reopen_a'.and.label==1)api_name='a_reopened_api.dat'
 if(step>nt/2.and.mode=='reopen_b'.and.label==2)api_name='b_reopened_api.dat'
 open(newunit=api_unit,file=trim(api_name),status='unknown',position='append')
 write(api_unit,'(5es26.17)')local_step*local_dt,j,energy
 close(api_unit)
 call rt_probe_write(state,unit)
 end subroutine
end program
