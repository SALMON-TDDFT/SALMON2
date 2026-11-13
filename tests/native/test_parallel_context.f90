program test_parallel_context
  use mpi
  use communication, only: comm_attach,comm_detach,comm_get_globalinfo,comm_sync_all,comm_init,comm_finalize
  implicit none
  integer :: ierr,rank,group,local_rank,n,solver,total,i,provided
  logical :: done
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_split(MPI_COMM_WORLD,merge(0,1,rank<2),rank,group,ierr)
  call comm_attach(group,ierr)
  if(ierr/=MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD,10,ierr)
  call comm_init ! Must not initialize MPI a second time.
  call comm_get_globalinfo(solver,local_rank,n)
  if(n/=merge(2,3,rank<2)) call MPI_Abort(MPI_COMM_WORLD,11,ierr)
  call MPI_Allreduce(local_rank,total,1,MPI_INTEGER,MPI_SUM,solver,ierr)
  if(total/=n*(n-1)/2) call MPI_Abort(MPI_COMM_WORLD,12,ierr)
  call comm_attach(group,ierr) ! Duplicate attachment must fail without losing context.
  if(ierr==MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD,13,ierr)
  do i=1,n ! Unequal iteration counts expose accidental world barriers.
    call comm_sync_all
  end do
  call comm_finalize ! Detach only; owner still owns MPI.
  call MPI_Finalized(done,ierr)
  if(done) stop 14
  call MPI_Barrier(group,ierr)
  call comm_detach(ierr)
  if(ierr==MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD,15,ierr)
  call MPI_Comm_free(group,ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(rank==0) print *, 'native subgroup lifecycle PASS'
  call MPI_Finalize(ierr)
end program
