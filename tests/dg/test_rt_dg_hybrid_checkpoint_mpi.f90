#include "config.h"
program test_rt_dg_hybrid_checkpoint_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64
  use rt_dg_hybrid_checkpoint,only:collective_rt_dg_hybrid_publication_precondition,&
    collective_rt_dg_hybrid_publication_mapping_precondition
  implicit none
  integer::rank,nproc,ierr,n,i,row_count
  integer,allocatable::owners(:)
  integer(int64),allocatable::rows(:),occupied_rows(:)
  logical::ok,local_valid
  character(256)::mode,message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call get_command_argument(1,mode)
  n=max(4,nproc);allocate(owners(n));owners=[(mod(i-1,nproc),i=1,n)]
  row_count=count(owners==rank);allocate(rows(row_count),occupied_rows(row_count))
  rows=pack([(int(i,int64),i=1,n)],owners==rank);occupied_rows=rows
  local_valid=.true.;if(trim(mode)=='bad_local'.and.rank==0)local_valid=.false.
  call collective_rt_dg_hybrid_publication_precondition(MPI_COMM_WORLD,local_valid,n,2,ok,message)
  if(trim(mode)=='bad_local')then
    if(.not.ok)error stop 'PASS named collective v5 publication precondition rejection'
    error stop 'one-rank malformed v5 publication was accepted'
  endif
  if(.not.ok)error stop trim(message)
  if(trim(mode)=='bad_mapping'.and.rank==0)owners(1)=mod(owners(1)+1,nproc)
  call collective_rt_dg_hybrid_publication_mapping_precondition(MPI_COMM_WORLD,n,rows,owners,&
    occupied_rows,.true.,ok,message)
  if(trim(mode)=='bad_mapping')then
    if(.not.ok)error stop 'PASS named collective v5 mapping rejection before publication'
    error stop 'rank-disagreeing v5 mapping was accepted'
  endif
  if(.not.ok)error stop trim(message)
  if(rank==0)write(*,'(a,i0,a)')'PASS v5 publication preconditions on ',nproc,' ranks'
  call MPI_Finalize(ierr)
end program test_rt_dg_hybrid_checkpoint_mpi
