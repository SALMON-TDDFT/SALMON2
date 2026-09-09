program test_rt_dg_hybrid_v3_rejection_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use rt_dg_hybrid_initialization,only:s_rt_dg_hybrid_state,initialize_rt_dg_hybrid_from_checkpoint
  implicit none
  type(s_rt_dg_hybrid_state)::state
  integer::ierr,rank,nproc
  logical::ok
  character(512)::message
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call initialize_rt_dg_hybrid_from_checkpoint(MPI_COMM_WORLD,'absent-dense-v3.chk','tddft_response',&
    .true.,1,.false.,.false.,.false.,.false.,.false.,[1],1_int64,1_int64,&
    [1d-8,1d-8,1d-8,1d-8],state,ok,message)
  if(ok.or.index(message,'dense Hybrid v3 checkpoint is unsupported')==0)then
    write(*,'(a)')trim(message);call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  endif
  if(rank==0)write(*,'(a,i0,a)')'PASS dense v3 early rejection on ',nproc,' ranks'
  call MPI_Finalize(ierr)
end program test_rt_dg_hybrid_v3_rejection_mpi
