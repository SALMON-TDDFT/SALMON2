! Copyright 2026 SALMONpy contributors. Apache-2.0.
program probe
 use mpi,only:MPI_Init,MPI_Finalize,MPI_Comm_rank,MPI_COMM_WORLD,MPI_COMM_SELF
 use communication,only:comm_attach,comm_detach
 use salmon_global
 use parallelization
 use inputoutput
 use math_constants
 use timer
 use plusU_global,only:read_Hubbard_parameters
 use sym_sub,only:read_sw_symmetry,use_symmetry,SymMatA,Ainv,Binv
 use salmonpy_tddft_context,only:s_rt_model_settings, &
   tddft_context,rt_open,rt_step,rt_close,rt_probe_write, &
   rt_checkpoint,rt_accepted_step,rt_current
 implicit none
 type(s_rt_model_settings)::owner_ctx,owner_second
 type(tddft_context),pointer::ctx=>null(),second=>null()
 character(16)::mode,split_text
 character(4096)::path
 integer::split,first,last,step,code,unit,calls,meta,world_rank,ierr
 character(32)::rank_directory
 real(8)::ac(3),j(3),energy
 call get_command_argument(1,mode)
 call get_command_argument(2,split_text);read(split_text,*)split
 call get_command_argument(3,path)
 call MPI_Init(ierr);if(ierr/=0)stop 20
 call MPI_Comm_rank(MPI_COMM_WORLD,world_rank,ierr)
 write(rank_directory,'(a,i0)')'rank-',world_rank
 call chdir(trim(rank_directory),ierr);if(ierr/=0)stop 21
 if(mode=='inner-reject')then
   call comm_attach(MPI_COMM_WORLD,ierr)
 else
   call comm_attach(MPI_COMM_SELF,ierr)
 endif
 if(ierr/=0)stop 22
 call set_math_constants
 call setup_parallel
 call read_input('inputfile')
 if(theory/='tddft_response'.or.xc/='hse06')stop 99
 calc_mode='RT'
 call read_Hubbard_parameters
 call read_sw_symmetry(yn_symmetry)
 call timer_initialize
 if(mode=='default')then
   call rt_open(ctx,'run',code, &
     model_settings=owner_ctx)
 else
   call rt_open(ctx,'run',code,allow_symmetry=.true., &
     model_settings=owner_ctx)
 endif
 if(mode=='default'.or.mode=='open-reject'.or.mode=='inner-reject')then
   if(code==0.or.associated(ctx))stop 23
   call end_parallel
   call comm_detach(ierr);call MPI_Finalize(ierr)
   print *,'Bounded symmetry TDDFT open correctly rejected'
   stop
 endif
 if(code/=0)stop 1
 if(mode=='lifetime')then
   call rt_open(second,'second',code,allow_symmetry=.true., &
     model_settings=owner_second)
   if(code==0.or.associated(second))stop 33
   mode='continuous'
 endif
 first=0;last=nt
 if(mode=='reject')then
   call rt_checkpoint(ctx,'before-invalid-load.bin',.false.,code);if(code/=0)stop 24
   call rt_checkpoint(ctx,trim(path),.true.,code)
   if(code==0)stop 8
   call rt_checkpoint(ctx,'after-invalid-load.bin',.false.,code);if(code/=0)stop 25
   ! Header, payload and footer are all validated before a numerical commit.
   if(rt_accepted_step(ctx)/=0)stop 9
   open(newunit=meta,file='calls.dat',status='new')
   write(meta,*)0,0,rt_accepted_step(ctx)
   close(meta)
   call rt_close(ctx)
   call end_parallel
   call comm_detach(ierr);call MPI_Finalize(ierr)
   print *,'Fast HSE corrupt/mismatch checkpoint rejected'
   stop
 endif
 if(mode=='resume'.or.mode=='staged')then
   if(mode=='staged')then
     open(newunit=unit,file='before-stage.bin',status='new',access='stream',form='unformatted')
     call rt_probe_write(ctx,unit);close(unit)
     call rt_checkpoint(ctx,trim(path),.true.,code,stage=2);if(code/=0)stop 26
     if(rt_accepted_step(ctx)/=0)stop 27
     if(use_symmetry)then
       use_symmetry=.false.
       call rt_checkpoint(ctx,'',.true.,code,stage=3)
       if(code==0.or.rt_accepted_step(ctx)/=0)stop 34
       use_symmetry=.true.
       Ainv=2d0*Ainv
       call rt_checkpoint(ctx,'',.true.,code,stage=3)
       if(code==0.or.rt_accepted_step(ctx)/=0)stop 38
       Ainv=Ainv/2d0
     endif
     open(newunit=unit,file='after-stage.bin',status='new',access='stream',form='unformatted')
     call rt_probe_write(ctx,unit);close(unit)
     call rt_checkpoint(ctx,'',.true.,code,stage=3);if(code/=0)stop 28
   else
     call rt_checkpoint(ctx,trim(path),.true.,code);if(code/=0)stop 2
   endif
   first=rt_accepted_step(ctx);if(first/=split)stop 3
 else if(mode=='save')then
   last=split
 else if(mode/='continuous'.and.mode/='field'.and.mode/='global')then
   stop 4
 endif
 calls=0
 open(newunit=unit,file='trace.bin',status='new',access='stream',form='unformatted')
 open(newunit=meta,file='api.dat',status='new')
 do step=first+1,last
   ac=[0d0,0d0,.0001d0*sin(.3d0*step)]
   if(mode=='global'.and.step==4)then
     call rt_checkpoint(ctx,'before-global.bin',.false.,code);if(code/=0)stop 29
     use_symmetry=.false.
     call rt_step(ctx,ac,j,energy,code)
     if(code/=2.or.rt_accepted_step(ctx)/=3)stop 30
     use_symmetry=.true.
     SymMatA(1,1,2)=-SymMatA(1,1,2)
     call rt_step(ctx,ac,j,energy,code)
     if(code/=2.or.rt_accepted_step(ctx)/=3)stop 31
     SymMatA(1,1,2)=-SymMatA(1,1,2)
     Ainv=2d0*Ainv
     call rt_step(ctx,ac,j,energy,code)
     if(code/=2.or.rt_accepted_step(ctx)/=3)stop 35
     call rt_current(ctx,j,code);if(code==0)stop 36
     Ainv=Ainv/2d0
     Binv=2d0*Binv
     call rt_step(ctx,ac,j,energy,code)
     if(code/=2.or.rt_accepted_step(ctx)/=3)stop 37
     Binv=Binv/2d0
     call rt_checkpoint(ctx,'after-global.bin' ,.false.,code);if(code/=0)stop 32
   endif
   if(mode=='field'.and.step==4)then
     call rt_checkpoint(ctx,'before-invalid.bin',.false.,code);if(code/=0)stop 11
     call rt_step(ctx,[1d-250,0d0,0d0],j,energy,code)
     if(code/=34.or.rt_accepted_step(ctx)/=3)stop 12
     call rt_checkpoint(ctx,'after-invalid.bin',.false.,code);if(code/=0)stop 13
   endif
   call rt_step(ctx,ac,j,energy,code);if(code/=0)stop 5
   calls=calls+1
   write(meta,'(5es26.17)')step*dt,j,energy
   call rt_probe_write(ctx,unit)
 enddo
 close(unit);close(meta)
 call rt_checkpoint(ctx,'final-state.bin',.false.,code);if(code/=0)stop 39
 if(mode=='save')then
   call rt_checkpoint(ctx,trim(path),.false.,code);if(code/=0)stop 6
 endif
 open(newunit=meta,file='calls.dat',status='new')
 write(meta,*)first,calls,rt_accepted_step(ctx)
 close(meta)
 call rt_close(ctx)
 call end_parallel
 call comm_detach(ierr);call MPI_Finalize(ierr)
 print *,'Fast HSE direct checkpoint probe PASS'
end program
