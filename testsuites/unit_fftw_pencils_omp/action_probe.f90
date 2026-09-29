program action_probe
  use mpi
  use omp_lib, only: omp_set_num_threads,omp_set_dynamic
  use exx_spatial, only: spatial_exx_state,spatial_exx_apply
  implicit none
  type(spatial_exx_state) :: op
  integer :: ierr,provided,np,rank,dims(2),coords(2),comm(2),n(3),ng,i,j,t,status
  complex(8),allocatable :: target(:,:,:),reference(:,:,:),value(:,:,:)
  real(8) :: omega
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  if(provided<MPI_THREAD_FUNNELED)error stop 'MPI thread support'
  call MPI_Comm_size(MPI_COMM_WORLD,np,ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  if(np/=4)error stop 'use 4 ranks'
  dims=[2,2];coords=[modulo(rank,2),rank/2];n=[128,64,64];ng=product(n)/np
  call MPI_Comm_split(MPI_COMM_WORLD,coords(2),coords(1),comm(1),ierr)
  call MPI_Comm_split(MPI_COMM_WORLD,coords(1),coords(2),comm(2),ierr)
  allocate(op%source(ng,4),target(ng,5,1),reference(ng,5,1),value(ng,5,1))
  do j=1,5;do i=1,ng
    target(i,j,1)=cmplx(sin(.017d0*i+j+rank),cos(.023d0*i-j+rank),8)/sqrt(real(product(n),8))
  enddo;enddo
  op%source=target(:,1:4,1)
  call omp_set_dynamic(.false.)
  do t=1,2
    omega=0d0
    if(t==2)omega=.11d0
    call omp_set_num_threads(1)
    call spatial_exx_apply(op,n,[.5d0,.5d0,.5d0],dims,coords,comm,MPI_COMM_WORLD,0d0, &
      target,reference,status,omega=omega)
    if(status/=0)error stop 'serial exchange'
    call omp_set_num_threads(4)
    call spatial_exx_apply(op,n,[.5d0,.5d0,.5d0],dims,coords,comm,MPI_COMM_WORLD,0d0, &
      target,value,status,omega=omega)
    if(status/=0.or.maxval(abs(value-reference))>2d-12)error stop 'threaded exchange differs'
    call spatial_exx_apply(op,n,[.5d0,.5d0,.5d0],dims,coords,comm,MPI_COMM_WORLD,0d0, &
      target(:,1:0,:),value(:,1:0,:),status,omega=omega)
    if(status/=0)error stop 'empty targets'
  enddo
  if(rank==0)print *, 'PASS large-point HSE/Coulomb action, padded tail and empty targets'
  call MPI_Finalize(ierr)
end program
