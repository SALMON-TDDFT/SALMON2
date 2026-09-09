program test_rt_dg_hybrid_checkpoint_v4_reject_mpi
  use mpi
  use rt_dg_hybrid_checkpoint_v4,only:s_rt_dg_hybrid_v4_shard,read_rt_dg_hybrid_checkpoint_v4
  implicit none
  type(s_rt_dg_hybrid_v4_shard)::payload
  integer::ierr,rank,nproc
  logical::ok
  character(512)::prefix,message
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call get_command_argument(1,prefix)
  call read_rt_dg_hybrid_checkpoint_v4(MPI_COMM_WORLD,trim(prefix),payload,ok,message)
  if(ok)call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  if(rank==0)write(*,'(a,i0,2a)')'PASS v4 collective rejection ranks=',nproc,' diagnostic=',trim(message)
  call MPI_Finalize(ierr)
end program test_rt_dg_hybrid_checkpoint_v4_reject_mpi
