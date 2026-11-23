! Copyright 2026 SALMONpy contributors. Apache-2.0.
program probe
 use salmon_global
 use parallelization
 use inputoutput
 use math_constants
 use timer
 use plusU_global,only:read_Hubbard_parameters
 use sym_sub,only:read_sw_symmetry
 use salmonpy_tddft_context,only:s_rt_model_settings, &
   tddft_context,rt_open,rt_step,rt_close,rt_probe_write, &
   rt_checkpoint,rt_accepted_step
 implicit none
 type(s_rt_model_settings)::owner_ctx
 type(tddft_context),pointer::ctx=>null()
 character(16)::mode,split_text
 character(4096)::path
 integer::split,first,last,step,code,unit,calls,meta
 real(8)::ac(3),j(3),energy
 call get_command_argument(1,mode)
 call get_command_argument(2,split_text);read(split_text,*)split
 call get_command_argument(3,path)
 call set_math_constants
 call setup_parallel
 call read_input
 if(theory/='tddft_response'.or.xc/='hse06')stop 99
 calc_mode='RT'
 call read_Hubbard_parameters
 call read_sw_symmetry(yn_symmetry)
 call timer_initialize
 call rt_open(ctx,'run',code, &
     model_settings=owner_ctx);if(code/=0)stop 1
 first=0;last=nt
 if(mode=='reject')then
   call rt_checkpoint(ctx,trim(path),.true.,code)
   if(code==0)stop 8
   ! Header, payload and footer are all validated before a numerical commit.
   if(rt_accepted_step(ctx)/=0)stop 9
   open(newunit=meta,file='calls.dat',status='new')
   write(meta,*)0,0,rt_accepted_step(ctx)
   close(meta)
   call rt_close(ctx)
   call end_parallel
   print *,'Fast HSE corrupt/mismatch checkpoint rejected'
   stop
 endif
 if(mode=='resume')then
   call rt_checkpoint(ctx,trim(path),.true.,code);if(code/=0)stop 2
   first=rt_accepted_step(ctx);if(first/=split)stop 3
 else if(mode=='save')then
   last=split
 else if(mode/='continuous')then
   stop 4
 endif
 calls=0
 open(newunit=unit,file='trace.bin',status='new',access='stream',form='unformatted')
 open(newunit=meta,file='api.dat',status='new')
 do step=first+1,last
   ac=[.0001d0*sin(.3d0*step),.00007d0*sin(.4d0*step),.00005d0*sin(.4d0*step)]
   call rt_step(ctx,ac,j,energy,code);if(code/=0)stop 5
   calls=calls+1
   write(meta,'(5es26.17)')step*dt,j,energy
   call rt_probe_write(ctx,unit)
 enddo
 close(unit);close(meta)
 if(mode=='save')then
   call rt_checkpoint(ctx,trim(path),.false.,code);if(code/=0)stop 6
 endif
 open(newunit=meta,file='calls.dat',status='new')
 write(meta,*)first,calls,rt_accepted_step(ctx)
 close(meta)
 call rt_close(ctx)
 call end_parallel
 print *,'Fast HSE direct checkpoint probe PASS'
end program
