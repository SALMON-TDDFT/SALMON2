! Internal same-model RT context independence and overlapping transactions.
! Copyright 2026 SALMONpy contributors. Apache-2.0.
program context_owners
 use mpi,only:MPI_Init,MPI_Finalize,MPI_Comm_rank,MPI_Comm_split,MPI_Comm_free, &
   MPI_Barrier,MPI_COMM_WORLD
 use communication,only:comm_attach,comm_detach
 use salmon_global
 use parallelization,only:setup_parallel,end_parallel,nproc_id_global,nproc_size_global
 use inputoutput,only:read_input
 use math_constants,only:set_math_constants
 use timer,only:timer_initialize
 use plusU_global,only:read_Hubbard_parameters
 use sym_sub,only:read_sw_symmetry
 use structures,only:s_rt_run_config
 use salmonpy_tddft_context,only:tddft_context,rt_open,rt_input_config,rt_step,rt_close, &
   rt_step_transaction,rt_checkpoint,rt_current,rt_accepted_step,rt_probe_write
 implicit none
 type(tddft_context),pointer::a=>null(),b=>null()
 type(s_rt_run_config)::ca,cb
 character(32)::mode,arg,directory,path_a,path_b
 character(:),allocatable::a_rt,a_energy,b_rt,b_energy
 real(8)::j(3),energy_a=0d0,energy_b=0d0,held_b(3)
 integer::ierr,world_rank,group_size,comm,group,step,status,ua,ub
 logical::dual,protected
 call get_command_argument(1,mode);call get_command_argument(2,arg);read(arg,*)group_size
 dual=mode=='dual'.or.mode=='protected';protected=mode=='protected'
 call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,world_rank,ierr)
 group=world_rank/group_size
 call MPI_Comm_split(MPI_COMM_WORLD,group,world_rank,comm,ierr)
 write(directory,'(a,i0)')'group-',group
 call chdir(trim(directory),ierr);if(ierr/=0)stop 20
 call comm_attach(comm,ierr);if(ierr/=0)stop 21
 call set_math_constants;call setup_parallel;call read_input('inputfile')
 calc_mode='RT';nproc_k=group_size;nproc_ob=1;nproc_rgrid=1
 call read_Hubbard_parameters;call read_sw_symmetry(yn_symmetry);call timer_initialize
 call rt_input_config(ca);ca%delta_t=.001d0;ca%max_steps=32
 cb=ca;cb%delta_t=.002d0;cb%max_steps=34
 ! Independent controls set the legacy dt/nt as well. The interleaved run leaves
 ! them at A's values and must obtain B's dt/capacity from the explicit context.
 if(mode=='single-b')then
   dt=cb%delta_t;nt=cb%max_steps
 else
   dt=ca%delta_t;nt=ca%max_steps
 endif
 if(dual.or.mode=='single-a'.or.index(mode,'late-')==1)then
   call rt_open(a,'a',status,config=ca);call require(status,1)
   if(nproc_id_global==0)open(newunit=ua,file='a-trace.dat',status='new')
 endif
 if(dual.or.mode=='single-b')then
   call rt_open(b,'b',status,config=cb);call require(status,2)
   if(nproc_id_global==0)open(newunit=ub,file='b-trace.dat',status='new')
 endif
 write(path_a,'(a,i0,a)')'a-',nproc_id_global,'.bin'
 write(path_b,'(a,i0,a)')'b-',nproc_id_global,'.bin'
 ! Open B only while A's rollback streams are closed, forcing unit reuse.
 if(index(mode,'late-')==1)then
   call rt_step_transaction(a,0,status);call require(status,201)
   call read_logs('a',a_rt,a_energy)
   call advance(a,1,1)
   call rt_step_transaction(a,1,status);call require(status,202)
   if(nproc_id_global==0)then
     call bytes_write('a_rt.data',a_rt)
     call bytes_write('a_rt_energy.data',a_energy)
   endif
   call MPI_Barrier(comm,ierr)
   call rt_open(b,'b',status,config=cb);call require(status,203)
   if(nproc_id_global==0)open(newunit=ub,file='b-trace.dat',status='new')
   if(mode=='late-failure'.or.mode=='late-retry')then
     if(nproc_id_global==0)then
       open(newunit=status,file='a_rt.data',status='old')
       close(status,status='delete')
     endif
     call MPI_Barrier(comm,ierr)
     call rt_step_transaction(a,2,status)
     if(status==0)stop 206
     if(mode=='late-retry')then
       if(nproc_id_global==0)call bytes_write('a_rt.data',a_rt)
       call MPI_Barrier(comm,ierr)
       call rt_step_transaction(a,2,status);call require(status,207)
     endif
   endif
   if(mode=='late-reopen')then
     call rt_step_transaction(a,2,status);call require(status,204)
     call advance(a,1,1)
   endif
   call rt_close(a)
   do step=1,33
     call advance(b,2,step)
     if(nproc_id_global==0)then
       call rt_current(b,j,status);call require_local(status,205)
       write(ub,'(i4,4es26.17)')step,j,energy_b
     endif
   enddo
   call dump(b,'b-final-')
   call rt_close(b)
   if(nproc_id_global==0)then
     close(ua);close(ub)
   endif
   goto 900
 endif
 do step=1,32
   if(protected.and.step==6)then
     call rt_step_transaction(a,0,status);call require(status,31)
     call rt_step_transaction(b,0,status);call require(status,32)
     call read_logs('a',a_rt,a_energy);call read_logs('b',b_rt,b_energy)
   endif
   if(associated(a))call advance(a,1,step)
   if(associated(b))call advance(b,2,step)
   if(protected.and.step==6)then
     call rt_current(b,held_b,status);call require(status,33)
     call rollback(a,'a',a_rt,a_energy)
     call advance(a,1,step)
     call rt_current(b,j,status);call require(status,34)
     if(any(j/=held_b).or.rt_accepted_step(b)/=step)stop 35
     call rollback(b,'b',b_rt,b_energy)
     call advance(b,2,step)
   endif
   if(protected.and.step==9)then
     call rt_checkpoint(a,trim(path_a),.false.,status);call require(status,41)
     call rt_checkpoint(b,trim(path_b),.false.,status);call require(status,42)
   endif
   if(protected.and.step==10)then
     call rt_current(b,held_b,status);call require(status,43)
     call rt_checkpoint(a,trim(path_a),.true.,status,stage=2);call require(status,44)
     call rt_checkpoint(b,trim(path_b),.true.,status,stage=2);call require(status,45)
     call rt_checkpoint(a,'',.true.,status,stage=3);call require(status,46)
     if(rt_accepted_step(a)/=9.or.rt_accepted_step(b)/=10)stop 47
     call advance(a,1,step)
     call rt_checkpoint(b,'',.true.,status,stage=4);call require(status,48)
     call rt_current(b,j,status);call require(status,49)
     if(any(j/=held_b))stop 50
   endif
   if(protected.and.step==12)then
     call rt_step_transaction(a,0,status);call require(status,51)
     call read_logs('a',a_rt,a_energy)
     ! B may stage a restore while only A has a protected step in flight.
     call rt_checkpoint(b,'b-other-'//trim(path_b),.false.,status)
     ! Only root-independent rank paths are valid for MPI shards.
     call require(status,52)
     call rt_checkpoint(b,'b-other-'//trim(path_b),.true.,status,stage=2)
     call require(status,53)
     call rt_checkpoint(b,'',.true.,status,stage=4);call require(status,54)
     call rt_step_transaction(a,2,status);call require(status,55)
   endif
   if(nproc_id_global==0)then
     if(associated(a))then
       call rt_current(a,j,status);call require_local(status,61)
       write(ua,'(i4,4es26.17)')step,j,energy_a
     endif
     if(associated(b))then
       call rt_current(b,j,status);call require_local(status,62)
       write(ub,'(i4,4es26.17)')step,j,energy_b
     endif
   endif
 enddo
 if(associated(a))call dump(a,'a-final-')
 if(dual)then
   if(protected)then
     call rt_step_transaction(a,0,status);call require(status,71)
     call rt_checkpoint(b,'b-final-'//trim(path_b),.false.,status);call require(status,72)
     call rt_checkpoint(b,'b-final-'//trim(path_b),.true.,status,stage=2);call require(status,73)
   endif
   call rt_close(a)
   if(protected)then
     call rt_checkpoint(b,'',.true.,status,stage=3);call require(status,74)
   endif
 endif
 if(associated(b))then
   call advance(b,2,33)
   if(nproc_id_global==0)then
     call rt_current(b,j,status);call require_local(status,75)
     write(ub,'(i4,4es26.17)')33,j,energy_b
   endif
 endif
 if(associated(b))call dump(b,'b-final-')
 if(associated(a))call rt_close(a)
 if(associated(b))call rt_close(b)
 if(nproc_id_global==0)then
   if(dual.or.mode=='single-a')close(ua)
   if(dual.or.mode=='single-b')close(ub)
 endif
900 continue
 call end_parallel;call comm_detach(ierr);call MPI_Comm_free(comm,ierr);call MPI_Finalize(ierr)
 print *,'Internal RT contexts PASS ',trim(mode)
contains
 subroutine dump(ctx,prefix)
 type(tddft_context),pointer,intent(in)::ctx
 character(*),intent(in)::prefix
 character(80)::path
 integer::u
 write(path,'(a,i0,a)')prefix,nproc_id_global,'.probe'
 open(newunit=u,file=trim(path),status='new',access='stream',form='unformatted')
 call rt_probe_write(ctx,u)
 close(u)
 end subroutine
 subroutine require(code,label)
 integer,intent(in)::code,label
 if(code/=0)then
   write(*,*)'context test failure',label,code
   stop 90
 endif
 end subroutine
 subroutine require_local(code,label)
 integer,intent(in)::code,label
 if(code/=0)then
   write(*,*)'context current failure',label,code
   stop 91
 endif
 end subroutine
 subroutine advance(ctx,which,index)
 type(tddft_context),pointer,intent(inout)::ctx
 integer,intent(in)::which,index
 real(8)::value(3),current(3),total
 integer::code
 value=0d0
 if(which==1)then
   value(3)=.0001d0*sin(.3d0*index)
 else
   value(1)=.0001d0*(cos(.2d0*index)-1d0)
   value(2)=.0001d0*sin(.2d0*index)
 endif
 call rt_step(ctx,value,current,total,code);call require(code,100+index)
 if(rt_accepted_step(ctx)/=index)stop 92
 if(which==1)energy_a=total
 if(which==2)energy_b=total
 end subroutine
 subroutine read_logs(prefix,first,second)
 character(*),intent(in)::prefix
 character(:),allocatable,intent(out)::first,second
 if(nproc_id_global==0)then
   call bytes_read(prefix//'_rt.data',first)
   call bytes_read(prefix//'_rt_energy.data',second)
 else
   first='';second=''
 endif
 end subroutine
 subroutine bytes_read(path,bytes)
 character(*),intent(in)::path
 character(:),allocatable,intent(out)::bytes
 integer::u,n
 open(newunit=u,file=path,status='old',access='stream',form='unformatted')
 inquire(unit=u,size=n);allocate(character(n)::bytes);read(u)bytes;close(u)
 end subroutine
 subroutine bytes_write(path,bytes)
 character(*),intent(in)::path,bytes
 integer::u
 open(newunit=u,file=path,status='replace',access='stream',form='unformatted')
 write(u)bytes;close(u)
 end subroutine
 subroutine rollback(ctx,prefix,first,second)
 type(tddft_context),pointer,intent(inout)::ctx
 character(*),intent(in)::prefix
 character(:),allocatable,intent(in)::first,second
 integer::code
 call rt_step_transaction(ctx,1,code);call require(code,81)
 if(nproc_id_global==0)then
   call bytes_write(prefix//'_rt.data',first)
   call bytes_write(prefix//'_rt_energy.data',second)
 endif
 call MPI_Barrier(comm,ierr)
 call rt_step_transaction(ctx,2,code);call require(code,82)
 end subroutine
end program
