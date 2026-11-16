! Copyright 2026 SALMONpy contributors. Apache-2.0.
! Dedicated serial fast-PR HSE process. Configuration comes from cwd/inputfile;
! stdin is reserved for commands, stdout contains prefixed protocol records.
! Native diagnostic output is deliberately retained between protocol records.
#include "config.h"
program fast_hse_worker
 use iso_fortran_env,only:input_unit,output_unit
 use salmon_global
 use parallelization,only:setup_parallel,end_parallel
 use inputoutput,only:read_input
 use math_constants,only:set_math_constants
 use timer,only:timer_initialize
 use plusU_global,only:read_Hubbard_parameters
 use sym_sub,only:read_sw_symmetry
 use salmonpy_tddft_context,only:tddft_context,rt_open,rt_step,rt_close,rt_current,rt_checkpoint,rt_accepted_step
 implicit none
 type(tddft_context),pointer::ctx=>null()
 integer::commands,ios,code,accepted
 real(8)::ac(3),current(3),energy
 character(8192)::line
 character(16)::command
 character(4096)::path
#if defined(USE_MPI) || defined(USE_OPENACC) || defined(USE_FFTW) || !defined(USE_HSE)
 call fail('BUILD',1)
#endif
 ! The audited fast-PR read_input has no filename argument and consumes unit *
 ! through EOF. Keep a separate descriptor to the original pipe before rebinding
 ! input_unit; never consume command bytes as namelist data. Unix host only.
 open(newunit=commands,file='/dev/stdin',status='old',action='read',iostat=ios)
 if(ios/=0)call fail('STDIN',ios)
 open(unit=input_unit,file='inputfile',status='old',action='read',iostat=ios)
 if(ios/=0)call fail('INPUT',ios)
 call set_math_constants
 call setup_parallel
 call read_input
 close(input_unit)
 if(theory/='tddft_response'.or.xc/='hse06')call fail('CONFIG',1)
 calc_mode='RT'
 call read_Hubbard_parameters
 call read_sw_symmetry(yn_symmetry)
 call timer_initialize
 call rt_open(ctx,'worker',code)
 if(code/=0)call fail('OPEN',code)
 call rt_current(ctx,current,code)
 if(code/=0)call fail('CURRENT',code)
 write(output_unit,'(a,1x,i0,3(1x,es26.17e3))')'SALMONPY_HSE_READY',rt_accepted_step(ctx),current
 flush(output_unit)
 do
  read(commands,'(a)',iostat=ios)line
  if(ios/=0)call fail('COMMAND_READ',ios)
  command=''
  read(line,*,iostat=ios)command
  if(ios/=0)call fail('COMMAND',ios)
  select case(trim(command))
  case('STEP')
   read(line,*,iostat=ios)command,ac
   if(ios/=0)call fail('STEP_INPUT',ios)
   call rt_step(ctx,ac,current,energy,code)
   if(code/=0)call fail('STEP',code)
   write(output_unit,'(a,1x,i0,4(1x,es26.17e3))')'SALMONPY_HSE_STEP',rt_accepted_step(ctx),current,energy
  case('SAVE','LOAD')
   ! Everything after the command is the path, including embedded spaces.
   line=adjustl(line)
   if(len_trim(line)<6.or.len_trim(line(6:))>len(path))call fail('PATH',1)
   if(line(5:5)/=' ')call fail('PATH',1)
   path=adjustl(line(6:))
   if(path(1:1)/='/')call fail('PATH',1)
   call rt_checkpoint(ctx,trim(path),command=='LOAD',code)
   if(code/=0)call fail(trim(command),code)
   if(command=='LOAD')then
    call rt_current(ctx,current,code)
    if(code/=0)call fail('CURRENT',code)
    write(output_unit,'(a,1x,i0,3(1x,es26.17e3))')'SALMONPY_HSE_LOADED',rt_accepted_step(ctx),current
   else
    write(output_unit,'(a,1x,i0)')'SALMONPY_HSE_SAVED',rt_accepted_step(ctx)
   endif
  case('CLOSE')
   accepted=rt_accepted_step(ctx)
   call rt_close(ctx)
   call end_parallel
   close(commands)
   write(output_unit,'(a,1x,i0)')'SALMONPY_HSE_CLOSED',accepted
   flush(output_unit)
   exit
  case default
   call fail('UNKNOWN_COMMAND',1)
  end select
  flush(output_unit)
 enddo
contains
 subroutine fail(operation,status)
 character(*),intent(in)::operation
 integer,intent(in)::status
 write(output_unit,'(a,1x,a,1x,i0)')'SALMONPY_HSE_ERROR',trim(operation),status
 flush(output_unit)
 stop 1
 end subroutine
end program
