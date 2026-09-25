! Bounded initialization-only diagnostic: no SCF or time propagation.
program fft_auto_probe
  use mpi
  use hse_exchange
  implicit none
  type(hse_kernel) :: op
  real(8),allocatable :: k(:,:)
  real(8) :: t
  integer :: rank,np,status,n,m,nk,i,x,y,z,start,count,ierr,provided
  character(32) :: arg
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,status)
  if(provided<MPI_THREAD_FUNNELED)call MPI_Abort(MPI_COMM_WORLD,1,status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,np,status)
  call get_command_argument(1,arg);read(arg,*)m
  n=12;nk=m**3;allocate(k(3,nk));i=0
  do z=0,m-1;do y=0,m-1;do x=0,m-1
    i=i+1;k(:,i)=2*acos(-1d0)*real([x,y,z],8)/real(m*n,8)
  enddo;enddo;enddo
  start=1+rank*nk/np;count=(rank+1)*nk/np-rank*nk/np
  t=MPI_Wtime()
  call hse_kernel_init(op,n,m,1d0,k,.11d0,8,ierr,start,count)
  if(ierr/=0)call MPI_Abort(MPI_COMM_WORLD,2,status)
  write(*,'(a,i0,a,l1,a,3es14.5)')'rank=',rank,' contiguous=',op%contiguous_fft, &
    ' trial strided contiguous init_seconds=',op%fft_trial_seconds,MPI_Wtime()-t
  call hse_kernel_destroy(op)
  call MPI_Finalize(status)
end program
