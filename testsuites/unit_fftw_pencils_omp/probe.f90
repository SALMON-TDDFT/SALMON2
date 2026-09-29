program probe
  use mpi
  use omp_lib, only: omp_set_num_threads,omp_set_dynamic
  use fftw_pencils, only: pencil_transform,pencil_clear,fftw_pencil_plans_created
  implicit none
  integer :: ierr,provided,np,rank,n(3),dims(2),coords(2),comm(2),nt,i,q,t,mode,status,before
  integer :: scale
  integer :: teams(4)=[1,3,4,1]
  complex(8),allocatable :: x(:,:),y(:,:),back(:,:),reference(:,:)
  logical :: spectral
  character(16) :: argument,scale_argument
  call get_command_argument(1,argument)
  call get_command_argument(2,scale_argument)
  scale=1
  if(len_trim(scale_argument)>0)read(scale_argument,*)scale
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  if(provided<MPI_THREAD_FUNNELED)error stop 'MPI thread support'
  call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  if(np/=1.and.np/=2.and.np/=4)error stop 'use 1,2,4 MPI ranks'
  dims=[1,np]
  if(np==4)dims=[2,2]
  coords=[modulo(rank,dims(1)),rank/dims(1)]
  call MPI_Comm_split(MPI_COMM_WORLD,coords(2),coords(1),comm(1),ierr)
  call MPI_Comm_split(MPI_COMM_WORLD,coords(1),coords(2),comm(2),ierr)
  n=[24,12,8]*scale;nt=product(n)/np
  allocate(x(nt,7),y(nt,7),back(nt,7),reference(nt,7))
  do q=1,7;do i=1,nt
    x(i,q)=cmplx(sin(.017d0*i+q+rank),cos(.093d0*i-q+rank),8)
  enddo;enddo
  call omp_set_dynamic(.false.)
  do mode=1,2
    spectral=mode==2
    call pencil_clear()
    do t=1,size(teams)
      call omp_set_num_threads(teams(t));before=fftw_pencil_plans_created
      call pencil_transform(n,dims,coords,comm,x,y,-1,status,spectral)
      if(status/=0)error stop 'forward'
      if(t==1)reference=y
      if(maxval(abs(y-reference))>2d-10)error stop 'threaded forward differs'
      call pencil_transform(n,dims,coords,comm,y,back,1,status,spectral)
      if(status/=0.or.maxval(abs(back-x))>2d-12)error stop 'roundtrip'
      if(argument/='serial'.and.scale==5.and.teams(t)>1.and.fftw_pencil_plans_created-before<=12) &
        error stop 'FFT still uses serial batch plans'
    enddo
  enddo
  call pencil_clear()
  call MPI_Comm_free(comm(1),ierr);call MPI_Comm_free(comm(2),ierr)
  if(rank==0)write(*,*)'PASS threaded FFT, team changes, partial batches and spectral layouts'
  call MPI_Finalize(ierr)
end program
