program probe
  use mpi
  use communication, only: comm_create_group,comm_free_group,comm_get_max
  use fftw_blocks, only: block_transform
  implicit none
  integer :: ierr,provided,rank,peers,dims(3),coords(3),groups(3),n(3),m(3),lo(3),p(3),mode(3)
  integer :: a,b,g,x,y,z,color,key,status,j
  complex(8),allocatable :: v(:,:),f(:,:),back(:,:)
  complex(8) :: expected
  real(8) :: err,pi,global_error(1)
  character(32) :: arg
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,peers,ierr)
  call get_command_argument(1,arg);read(arg,*)dims
  if(product(dims)/=peers)error stop 'rank count'
  coords=[mod(rank,dims(1)),mod(rank/dims(1),dims(2)),rank/(dims(1)*dims(2))]
  do a=1,3
    color=0
    do b=1,3
      if(b/=a)color=color*dims(b)+coords(b)
    enddo
    key=coords(a);groups(a)=comm_create_group(MPI_COMM_WORLD,color,key)
  enddo
  n=[16,8,8]
  call get_command_argument(2,arg)
  if(len_trim(arg)>0)read(arg,*)n
  m=n/dims;lo=coords*m;pi=acos(-1d0)
  allocate(v(product(m),3),f(product(m),3),back(product(m),3))
  do b=1,3
    mode=[b,mod(b,2),1]
    g=0
    do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
      g=g+1;p=lo+[x,y,z]
      v(g,b)=exp(cmplx(0d0,2*pi*sum(real(mode*p,8)/n),8))
    enddo;enddo;enddo
  enddo
  call block_transform(n,dims,coords,groups,v,f,-1,status)
  if(status/=0)error stop 'forward'
  err=0d0
  do b=1,3
    mode=[b,mod(b,2),1];g=0
    do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
      g=g+1;p=lo+[x,y,z];expected=0d0
      if(all(p==mode))expected=real(product(n),8)
      err=max(err,abs(f(g,b)-expected)/product(n))
    enddo;enddo;enddo
  enddo
  call block_transform(n,dims,coords,groups,f,back,1,status)
  if(status/=0)error stop 'inverse'
  err=max(err,maxval(abs(back-v)))
  call comm_get_max([err],global_error,1,MPI_COMM_WORLD)
  err=global_error(1)
  if(err>2d-12)error stop 'FFT mismatch'
  do b=1,3
    j=1
    if(b==2)j=0
    if(b==3)j=3
    call block_transform(n,dims,coords,groups,v(:,:j),f(:,:j),-1,status)
    if(status/=0)error stop 'cache forward'
    call block_transform(n,dims,coords,groups,f(:,:j),back(:,:j),1,status)
    if(status/=0)error stop 'cache inverse'
    if(j>0)then
      if(maxval(abs(back(:,:j)-v(:,:j)))>2d-12)error stop 'cache mismatch'
    endif
  enddo
  if(rank==0)write(*,*)'PASS block FFT',dims,err
  do a=1,3
    call comm_free_group(groups(a))
  enddo
  call MPI_Finalize(ierr)
end program
