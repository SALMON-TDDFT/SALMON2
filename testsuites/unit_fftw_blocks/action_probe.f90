program action_probe
  use mpi
  use communication, only: comm_create_group,comm_free_group,comm_get_max
  use exx_spatial, only: spatial_exx_state,spatial_exx_apply
  use exx_local_fft, only: compact_kernel_bounds
  implicit none
  type(spatial_exx_state) :: op,ref
  integer :: ierr,provided,rank,peers,dims(3),coords(3),groups(3),n(3),m(3),lo(3),p(3),distance(3)
  integer :: a,b,g,x,y,z,color,status,q,j,t,compact,lower(3),upper(3)
  complex(8),allocatable :: target(:,:,:),reference(:,:,:),value(:,:,:),all_target(:,:,:)
  real(8) :: omega,err,errors(1)
  character(32) :: arg,mode
  call compact_kernel_bounds([128,16,16],[31,16,16],lower,upper)
  if(any(lower/=[-30,0,0]).or.any(upper/=[30,15,15]))error stop 'kernel bounds'
  if(product(upper-lower+1)/=15616)error stop 'duplicate periodic kernel entries'
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,peers,ierr)
  call get_command_argument(1,arg);read(arg,*)dims
  if(product(dims)/=peers)error stop 'rank count'
  coords=[mod(rank,dims(1)),mod(rank/dims(1),dims(2)),rank/(dims(1)*dims(2))]
  do a=1,3
    color=0
    do b=1,3
      if(b/=a)color=color*dims(b)+coords(b)
    enddo
    groups(a)=comm_create_group(MPI_COMM_WORLD,color,coords(a))
  enddo
  call get_command_argument(2,mode)
  n=[16,16,16];m=n/dims;lo=coords*m
  allocate(ref%source(product(n),2),all_target(product(n),3,1),reference(product(n),3,1))
  allocate(op%source(product(m),2),target(product(m),3,1),value(product(m),3,1))
  g=0
  do z=0,n(3)-1;do y=0,n(2)-1;do x=0,n(1)-1
    g=g+1;p=[x,y,z]
    do j=1,2
      distance=modulo(p-[0,8,8]-(j-1)*[7,0,0]+n/2,n)-n/2
      ref%source(g,j)=0d0
      if(sum(distance**2)<=2)ref%source(g,j)=cmplx(.1d0,.01d0*j,8)
      if(trim(mode)=='slab'.and.abs(distance(1))<=1)ref%source(g,j)=cmplx(.1d0,.01d0*j,8)
    enddo
    do j=1,3
      all_target(g,j,1)=cmplx(sin(.1d0*x+j),cos(.2d0*y+z-j),8)
    enddo
  enddo;enddo;enddo
  g=0
  do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
    g=g+1;p=lo+[x,y,z];q=1+p(1)+n(1)*(p(2)+n(2)*p(3))
    op%source(g,:)=ref%source(q,:);target(g,:,1)=all_target(q,:,1)
  enddo;enddo;enddo
  do t=1,2
    omega=0d0
    if(t==2)omega=.11d0
    call spatial_exx_apply(ref,n,[.5d0,.5d0,.5d0],[1,1],[0,0],[MPI_COMM_SELF,MPI_COMM_SELF], &
      MPI_COMM_SELF,0d0,all_target,reference,status,omega=omega)
    if(status/=0)error stop 'reference exchange'
    do compact=0,1
      op%compact=compact==1
      call spatial_exx_apply(op,n,[.5d0,.5d0,.5d0],dims,coords,groups,MPI_COMM_WORLD,0d0, &
        target,value,status,omega=omega)
      if(status/=0)error stop 'block exchange'
      g=0;err=0d0
      do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
        g=g+1;p=lo+[x,y,z];q=1+p(1)+n(1)*(p(2)+n(2)*p(3))
        err=max(err,maxval(abs(value(g,:,1)-reference(q,:,1))))
      enddo;enddo;enddo
      call comm_get_max([err],errors,1,MPI_COMM_WORLD)
      if(errors(1)>2d-12)error stop 'exchange mismatch'
      if(compact==1.and.op%local_pairs/=6)error stop 'compact pairs not executed'
      if(rank==0)write(*,*)'PASS exchange',dims,t,compact,errors(1)
    enddo
  enddo
  call spatial_exx_apply(op,n,[.5d0,.5d0,.5d0],dims,coords,groups,MPI_COMM_WORLD,0d0, &
    target(:,1:0,:),value(:,1:0,:),status,omega=omega)
  if(status/=0)error stop 'empty targets'
  do a=1,3
    call comm_free_group(groups(a))
  enddo
  call MPI_Finalize(ierr)
end program
