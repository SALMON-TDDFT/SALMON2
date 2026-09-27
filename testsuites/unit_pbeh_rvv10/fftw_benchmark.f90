program benchmark
  use mpi
  use fftw_pencils
  use rvv10_distributed
  implicit none
  integer :: ierr,rank,np,dims(2),coords(2),comm(2),n(3),m(3),lo(3),nt,q,i,j,status,reps,scale
  complex(8),allocatable :: input(:,:),output(:,:),back(:,:),a(:),b(:)
  real(8),allocatable :: rho(:),sigma(:),e(:),v(:),w(:),er(:),vr(:),wr(:)
  real(8) :: start,first,warm(2),whole(2),totals(8),maximum(8),err
  character(16) :: argument
  logical :: used
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  if(np/=2.and.np/=4)error stop 'use 2 or 4 ranks'
  scale=1
  call get_command_argument(1,argument)
  if(len_trim(argument)>0)read(argument,*)scale
  if(scale<1.or.scale>4)error stop 'benchmark scale must be1..4'
  n=[32,24,16]*scale;dims=[np/2,2];coords=[modulo(rank,dims(1)),rank/dims(1)]
  call MPI_Comm_split(MPI_COMM_WORLD,coords(2),coords(1),comm(1),ierr)
  call MPI_Comm_split(MPI_COMM_WORLD,coords(1),coords(2),comm(2),ierr)
  m=[n(1),n(2)/dims(1),n(3)/dims(2)];lo=[1,coords(1)*m(2)+1,coords(2)*m(3)+1];nt=product(m)
  allocate(input(nt,32),output(nt,32),back(nt,32),a(nt),b(nt),rho(nt),sigma(nt),e(nt),v(nt),w(nt), &
    er(nt),vr(nt),wr(nt))
  do q=1,32;do i=1,nt
    input(i,q)=cmplx(sin(.017d0*i+q+rank),cos(.09d0*i-q+rank),8)
  enddo;enddo
  rho=.01d0+.002d0*real(input(:,1));sigma=.00002d0*(1+.4d0*real(input(:,2)))
  call MPI_Barrier(MPI_COMM_WORLD,ierr);start=MPI_Wtime()
  call pencil_transform(n,dims,coords,comm,input,output,-1,status)
  if(status/=0)error stop 'FFTW first transform failed'
  first=MPI_Wtime()-start
  reps=20;fftw_pencil_seconds(2:4)=0d0
  call MPI_Barrier(MPI_COMM_WORLD,ierr);start=MPI_Wtime()
  do j=1,reps
    call pencil_transform(n,dims,coords,comm,input,output,-1,status)
    if(status/=0)error stop 'FFTW forward failed'
    call pencil_transform(n,dims,coords,comm,output,back,1,status)
    if(status/=0)error stop 'FFTW inverse failed'
  enddo
  warm(2)=(MPI_Wtime()-start)/reps
  if(maxval(abs(back-input))>1d-11)error stop 'benchmark roundtrip'
  call MPI_Barrier(MPI_COMM_WORLD,ierr);start=MPI_Wtime()
  do j=1,reps;do q=1,32
    a=input(:,q)
    call pzfft3dv_rvv10(a,b,n(1),n(2),n(3),dims(1),dims(2),-1,comm(1),comm(2))
    call pzfft3dv_rvv10(b,a,n(1),n(2),n(3),dims(1),dims(2),1,comm(1),comm(2))
  enddo;enddo
  warm(1)=(MPI_Wtime()-start)/reps
  if(maxval(abs(a-input(:,32)))>1d-11)error stop 'FFTE roundtrip'
  totals=[first,fftw_pencil_seconds(1),warm,fftw_pencil_seconds(2:4)/reps,real(fftw_pencil_plans_created,8)]
  call MPI_Allreduce(totals,maximum,8,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  if(rank==0)write(*,'(a,8es20.10)')'FFT first/setup/FFTE/FFTW/localFFT/MPI/packing/plans: ',maximum
  fftw_pencil_seconds(2:4)=0d0
  call MPI_Barrier(MPI_COMM_WORLD,ierr);start=MPI_Wtime()
  do j=1,reps
    call pencil_transform(n,dims,coords,comm,input,output,-1,status,spectral_z=.true.)
    if(status/=0)error stop 'FFTW Z forward failed'
    call pencil_transform(n,dims,coords,comm,output,back,1,status,spectral_z=.true.)
    if(status/=0)error stop 'FFTW Z inverse failed'
  enddo
  totals(1:4)=[(MPI_Wtime()-start)/reps,fftw_pencil_seconds(2:4)/reps]
  if(maxval(abs(back-input))>1d-11)error stop 'Z benchmark roundtrip'
  call MPI_Allreduce(totals,maximum,4,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  if(rank==0)write(*,'(a,4es20.10)')'FFT Z pair/localFFT/MPI/packing: ',maximum(1:4)
  do q=1,2
    call MPI_Barrier(MPI_COMM_WORLD,ierr);start=MPI_Wtime()
    do j=1,3
      call rvv10_evaluate_distributed(n,lo,m,[1,dims],[0,coords],[MPI_COMM_SELF,comm],MPI_COMM_WORLD, &
        [.7d0,.9d0,1.1d0],rho,sigma,5.3d0,.0093d0,32,e,v,w,used,status,q==2)
      if(.not.used.or.status/=0)error stop 'benchmark functional failed'
    enddo
    whole(q)=(MPI_Wtime()-start)/3
    if(q==1)then
      er=e;vr=v;wr=w
    else
      err=max(maxval(abs(e-er)),maxval(abs(v-vr)),maxval(abs(w-wr)))
      if(err>1d-10)error stop 'backend functional mismatch'
    endif
  enddo
  call MPI_Allreduce(whole,maximum,2,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  if(rank==0)write(*,'(a,2es20.10)')'RVV10 whole FFTE/FFTW: ',maximum(1:2)
  call pencil_clear();call MPI_Finalize(ierr)
end program
