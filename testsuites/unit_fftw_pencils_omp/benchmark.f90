program benchmark
  use mpi
  use omp_lib, only: omp_set_num_threads,omp_set_dynamic
  use fftw_pencils, only: pencil_transform,pencil_clear,fftw_pencil_seconds
  implicit none
  integer :: ierr,provided,np,rank,n(3),dims(2),coords(2),comm(2),nt,i,q,t,j,status
  integer :: teams(3)=[1,2,4]
  complex(8),allocatable :: x(:,:),y(:,:),back(:,:)
  real(8) :: start,times(4),maximum(4)
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  if(provided<MPI_THREAD_FUNNELED)error stop 'MPI thread support'
  call MPI_Comm_size(MPI_COMM_WORLD,np,ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  if(np/=4)error stop 'use 4 MPI ranks'
  dims=[2,2];coords=[modulo(rank,2),rank/2]
  call MPI_Comm_split(MPI_COMM_WORLD,coords(2),coords(1),comm(1),ierr)
  call MPI_Comm_split(MPI_COMM_WORLD,coords(1),coords(2),comm(2),ierr)
  n=[128,64,64];nt=product(n)/np
  allocate(x(nt,4),y(nt,4),back(nt,4))
  do q=1,4;do i=1,nt
    x(i,q)=cmplx(sin(.017d0*i+q+rank),cos(.093d0*i-q+rank),8)
  enddo;enddo
  call omp_set_dynamic(.false.)
  do t=1,size(teams)
    call omp_set_num_threads(teams(t))
    call pencil_transform(n,dims,coords,comm,x,y,-1,status,spectral_z=.true.)
    if(status/=0)error stop 'warmup'
    fftw_pencil_seconds(2:4)=0d0
    call MPI_Barrier(MPI_COMM_WORLD,ierr);start=MPI_Wtime()
    do j=1,10
      call pencil_transform(n,dims,coords,comm,x,y,-1,status,spectral_z=.true.)
      if(status/=0)error stop 'forward'
      call pencil_transform(n,dims,coords,comm,y,back,1,status,spectral_z=.true.)
      if(status/=0)error stop 'inverse'
    enddo
    times=[MPI_Wtime()-start,fftw_pencil_seconds(2:4)]/10d0
    if(maxval(abs(back-x))>2d-12)error stop 'roundtrip'
    call MPI_Allreduce(times,maximum,4,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(rank==0)write(*,'(a,i3,4es18.8)')'OMP / roundtrip / FFT / MPI / packing: ',teams(t),maximum
  enddo
  call pencil_clear()
  call MPI_Comm_free(comm(1),ierr);call MPI_Comm_free(comm(2),ierr)
  call MPI_Finalize(ierr)
end program
