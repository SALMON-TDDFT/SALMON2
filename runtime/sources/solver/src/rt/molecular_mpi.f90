! Copyright 2026 SALMONpy contributors. Apache-2.0.
#include "config.h"
module salmonpy_molecular_mpi
 use iso_c_binding
 implicit none
 private
 public::sp_molecular_mpi_protocol,sp_molecular_mpi_run
 ! The blocking molecular entry reserves a dedicated process for one lifetime.
 logical::used=.false.
contains
 integer(c_int) function sp_molecular_mpi_protocol() bind(C) result(version)
  version=0
#ifdef USE_MPI
#ifndef USE_OPENACC
  version=3
#endif
#endif
 end function
 integer(c_int) function sp_molecular_mpi_run(group,transport,controller,path,decomposition) bind(C) result(status)
#ifdef USE_MPI
  use iso_fortran_env,only:output_unit
  use mpi,only:MPI_Initialized,MPI_Finalized,MPI_Comm_size,MPI_Comm_test_inter, &
    MPI_Allreduce,MPI_Bcast,MPI_INTEGER,MPI_MAX,MPI_CHARACTER,MPI_SUCCESS
  use communication,only:comm_attach,comm_detach
  use parallelization,only:setup_parallel
  use math_constants,only:set_math_constants
  use inputoutput,only:read_input
  use salmon_global,only:quiet,yn_symmetry,calc_mode
  use sym_sub,only:read_sw_symmetry
  use timer,only:timer_initialize
  use tddft_worker,only:worker_channel_context,worker_channel_begin,worker_channel_end
#endif
  integer(c_int),value::group,transport,controller
  character(c_char),intent(in)::path(*)
  integer(c_int),intent(in)::decomposition(5)
#ifdef USE_MPI
  type(worker_channel_context),target::channel
  character(4096)::filename,reference
  integer::i,n,size_control,ierr,bad,total,detach_status,reference_decomposition(5)
  logical::initialized,finalized,inter,exists,opened
#endif
  status=15
#if defined(USE_MPI) && !defined(USE_OPENACC)
  call MPI_Initialized(initialized,ierr)
  call MPI_Finalized(finalized,ierr)
  if(.not.initialized.or.finalized)return
  call MPI_Comm_size(group,n,ierr)
  if(ierr/=MPI_SUCCESS)return
  call MPI_Comm_size(transport,size_control,ierr)
  if(ierr/=MPI_SUCCESS)return
  bad=0
  if(used)bad=12
  if(controller<0.or.controller>=size_control)bad=10
  if(any(decomposition<1))bad=10
  if(product(real(decomposition,8))/=real(n,8).or.decomposition(1)/=1)bad=10
  call MPI_Comm_test_inter(group,inter,ierr)
  if(inter.or.ierr/=MPI_SUCCESS)bad=10
  filename=''
  do i=1,len(filename)
   if(path(i)==c_null_char)exit
   filename(i:i)=path(i)
  enddo
  if(i>len(filename).or.len_trim(filename)==0)bad=10
  inquire(file=trim(filename),exist=exists)
  if(.not.exists)bad=10
  inquire(unit=30,opened=opened)
  if(opened)bad=15
  reference=filename;reference_decomposition=decomposition
  call MPI_Bcast(reference,len(reference),MPI_CHARACTER,0,group,ierr)
  call MPI_Bcast(reference_decomposition,5,MPI_INTEGER,0,group,ierr)
  if(reference/=filename.or.any(reference_decomposition/=decomposition))bad=11
  call MPI_Allreduce(bad,total,1,MPI_INTEGER,MPI_MAX,group,ierr)
  status=total
  if(status/=0)return
  call comm_attach(group,ierr)
  if(ierr/=MPI_SUCCESS)then
   status=16;return
  endif
  used=.true.
  channel%transport=transport;channel%controller=controller;channel%decomposition=decomposition
  call worker_channel_begin(channel)
  call set_math_constants
  call setup_parallel
  call read_input(trim(filename))
  calc_mode='RT'
  call read_sw_symmetry(yn_symmetry)
  call timer_initialize
  quiet=.false.
  call main_tddft
  flush(output_unit)
  call worker_channel_end()
  call comm_detach(detach_status)
  status=0
  if(detach_status/=MPI_SUCCESS)status=16
#endif
 end function
end module
