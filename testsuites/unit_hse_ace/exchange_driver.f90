program exchange_driver
  use mpi
  use hse_spatial
  use hse_wannier
  use hse_wannier_gauge, only: gauge_transport
  implicit none
  integer,parameter :: n(3)=[8,8,8],no=3,nt=5
  real(8),parameter :: h(3)=[.7d0,.8d0,.9d0]
  complex(8) :: psi(product(n),no,1),target(product(n),nt,1),reference(product(n),nt,1)
  complex(8),allocatable :: local(:,:,:),trial(:,:,:),action(:,:,:),previous_saved(:,:,:)
  type(spatial_exx_state) :: spatial,full_spatial
  type(s_hse_wannier) :: serial
  integer :: np,rank,err,status,dims(2),coords(2),comm(2),m(3),lo(3),g,l,x,y,z,j,k,stage
  real(8) :: error,global_error,dv,bad_dv,minimum,omega
  character(32) :: argument
  complex(8) :: transported(no,no,1)
  call get_command_argument(1,argument)
  read(argument,*)omega
  call MPI_Init(err)
  call MPI_Comm_size(MPI_COMM_WORLD,np,err)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,err)
  dims=[np,1]
  if(np==4)dims=[2,2]
  coords=[modulo(rank,dims(1)),rank/dims(1)]
  call MPI_Comm_split(MPI_COMM_WORLD,coords(2),coords(1),comm(1),err)
  call MPI_Comm_split(MPI_COMM_WORLD,coords(1),coords(2),comm(2),err)
  m=[n(1),n(2)/dims(1),n(3)/dims(2)];lo=[0,coords(1)*m(2),coords(2)*m(3)]
  dv=product(h)
  do j=1,no
    do g=1,product(n)
      psi(g,j,1)=cmplx(sin(.07d0*g*j)+cos(.03d0*g*(j+1)),sin(.05d0*g*(j+2)),8)
    enddo
    do k=1,j-1
      psi(:,j,1)=psi(:,j,1)-psi(:,k,1)*sum(conjg(psi(:,k,1))*psi(:,j,1))*dv
    enddo
    psi(:,j,1)=psi(:,j,1)/sqrt(sum(abs(psi(:,j,1))**2)*dv)
  enddo
  do j=1,nt
    do g=1,product(n)
      target(g,j,1)=cmplx(cos(.013d0*g*j),sin(.017d0*g*(j+1)),8)
    enddo
  enddo
  allocate(local(product(m),no,1),trial(product(m),nt,1),action(product(m),nt,1))
  call wannier_init(serial,n,[1,1,1],h,reshape([0d0,0d0,0d0],[3,1]),omega,status,2.5d0)
  if(status/=0)error stop 'serial init'
  do stage=1,2
    if(stage==2)psi=psi*cmplx(cos(.13d0),sin(.13d0),8)
    call wannier_refresh_source(serial,psi,reshape([2d0,2d0,2d0],[no,1]),3,1d-7,status)
    if(status/=0)error stop 'serial refresh'
    call wannier_apply(serial,target,reference,status)
    if(status/=0)error stop 'serial action'
    l=0
    do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
      l=l+1;g=1+x+n(1)*(y+lo(2)+n(2)*(z+lo(3)))
      local(l,:,1)=psi(g,:,1);trial(l,:,1)=target(g,:,1)
    enddo;enddo;enddo
    call spatial_exx_refresh(spatial,n,h,dims,coords,comm,MPI_COMM_WORLD,local, &
      merge(3,0,stage==1),1d-7,status)
    if(status/=0)error stop 'spatial refresh'
    call spatial_exx_refresh(full_spatial,n,h,[1,1],[0,0],[MPI_COMM_SELF,MPI_COMM_SELF], &
      MPI_COMM_SELF,psi,merge(3,0,stage==1),1d-7,status)
    if(status/=0)error stop 'full grid localization'
    l=0
    do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
      l=l+1;g=1+x+n(1)*(y+lo(2)+n(2)*(z+lo(3)))
      if(maxval(abs(spatial%source(l,:)-full_spatial%source(g,:)))>1d-10)error stop 'localization mismatch'
    enddo;enddo;enddo
    if(stage==1)then
      previous_saved=spatial%previous
    else
      if(abs(spatial%min_singular-1d0)>1d-10)error stop 'transport overlap mismatch'
      if(maxval(abs(spatial%previous-previous_saved))>1d-10)error stop 'transport gauge mismatch'
    endif
    call spatial_exx_apply(spatial,n,h,dims,coords,comm,MPI_COMM_WORLD,2.5d0,trial,action,status,omega=omega)
    if(status/=0)error stop 'spatial action'
    error=0d0;l=0
    do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
      l=l+1;g=1+x+n(1)*(y+lo(2)+n(2)*(z+lo(3)))
      error=max(error,maxval(abs(action(l,:,1)-reference(g,:,1))))
    enddo;enddo;enddo
    call MPI_Allreduce(error,global_error,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,err)
    if(global_error>1d-10)error stop 'exchange mismatch'
    if(rank==0)print *, 'PASS spatial exchange ranks/stage/error ',np,stage,global_error,omega
  enddo
  call spatial_exx_apply(spatial,n,h,dims,coords,comm,MPI_COMM_WORLD,2.5d0,trial,action,status,omega=-.1d0)
  if(status==0)error stop 'negative screening accepted'
  bad_dv=dv
  if(rank==0)bad_dv=-1d0
  call gauge_transport(local,spatial%previous,bad_dv,transported,minimum,status,sum_grid)
  if(status==0)error stop 'invalid local transport volume accepted'
  call wannier_destroy(serial)
  call MPI_Comm_free(comm(1),err);call MPI_Comm_free(comm(2),err)
  call MPI_Finalize(err)
contains
  subroutine sum_grid(a)
    complex(8),intent(inout) :: a(:,:)
    integer :: code
    call MPI_Allreduce(MPI_IN_PLACE,a,size(a),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,code)
  end subroutine
end program
