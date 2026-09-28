program probe
 use mpi
 use iso_fortran_env, only: int64
 use fftw_pencils
 use exx_spatial_local
 implicit none
 type(s_exx_spatial_local) :: plan
 integer :: n(3)=[16,12,8],dims(2),coords(2),comm(2),m(3),lo(3),p(3),rank,np,ierr,g,x,y,z,j,status,kind
 integer(int64) :: pairs,points
 complex(8),allocatable :: source(:),target(:,:),action(:,:),ref(:,:),work(:,:),spec(:,:)
 real(8),allocatable :: multiplier(:)
 complex(8) :: metric(7,7),total(7,7)
 real(8) :: q2,err,err_all
 logical :: used
 call MPI_Init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 dims=[np,1];if(np==4)dims=[2,2]
 coords=[mod(rank,dims(1)),rank/dims(1)]
 call MPI_Comm_split(MPI_COMM_WORLD,coords(2),coords(1),comm(1),ierr)
 call MPI_Comm_split(MPI_COMM_WORLD,coords(1),coords(2),comm(2),ierr)
 m=[n(1),n(2)/dims(1),n(3)/dims(2)];lo=[0,coords(1)*m(2),coords(2)*m(3)]
 allocate(source(product(m)),target(product(m),7),action(product(m),7),ref(product(m),7), &
 work(product(m),7),spec(product(m),7),multiplier(product(m)))
 g=0
 do y=0,n(2)/dims(2)-1;do x=0,n(1)/dims(1)-1;do z=0,n(3)-1
 g=g+1;p=[x+coords(1)*n(1)/dims(1),y+coords(2)*n(2)/dims(2),z]
 where(p>=(n+1)/2)p=p-n
 q2=sum((2*acos(-1d0)*p/n)**2);multiplier(g)=1/(.7d0+q2)
 enddo;enddo;enddo
 call spatial_local_init(plan,n,dims,coords,comm,multiplier,status)
 if(status/=0)error stop 'init'
 do kind=1,3
 g=0
 do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
 g=g+1;p=[x,y,z]+lo
 source(g)=0
 if(all(modulo(p+1,n)<=2))source(g)=cmplx(.3d0+.01d0*sum(p),.2d0,8)
 if(kind==2)source(g)=source(g)*exp(cmplx(0d0,.83d0,8))
 if(kind==3)source(g)=1
 do j=1,6
 target(g,j)=cmplx(sin(real(sum(p*[2,3,5])+j,8)),cos(real(sum(p*[5,1,2])+2*j,8)),8)
 enddo
 target(g,7)=0
 enddo;enddo;enddo
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target,action,used,status,pairs,points)
 if(status/=0)error stop 'apply'
 if(kind==3)then
 if(used.or.any(action/=(0d0,0d0)))error stop 'broad fallback'
 cycle
 endif
 if(.not.used.or.pairs/=6.or.points>=6*product(n))error stop 'compact counter'
 do j=1,7
 work(:,j)=conjg(source)*target(:,j)
 enddo
 call pencil_transform(n,dims,coords,comm,work(:,1:4),spec(:,1:4),-1,status,spectral_z=.true.)
 call pencil_transform(n,dims,coords,comm,work(:,5:7),spec(:,5:7),-1,status,spectral_z=.true.)
 if(status/=0)error stop 'forward'
 do j=1,7
 spec(:,j)=spec(:,j)*multiplier
 enddo
 call pencil_transform(n,dims,coords,comm,spec(:,1:4),work(:,1:4),1,status,spectral_z=.true.)
 call pencil_transform(n,dims,coords,comm,spec(:,5:7),work(:,5:7),1,status,spectral_z=.true.)
 if(status/=0)error stop 'inverse'
 do j=1,7
 ref(:,j)=-source*work(:,j)
 enddo
 err=maxval(abs(ref-action));call MPI_Allreduce(err,err_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
 if(err_all>1d-11)error stop 'global parity'
 metric=matmul(conjg(transpose(target)),action)
 call MPI_Allreduce(metric,total,49,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
 if(maxval(abs(total-conjg(transpose(total))))>1d-10)error stop 'Hermiticity'
 enddo
 source=0
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target,action,used,status,pairs,points)
 if(status/=0.or..not.used.or.pairs/=0.or.any(action/=(0d0,0d0)))error stop 'zero support'
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target(:,1:0),action(:,1:0),used,status,pairs,points)
 if(status/=0.or..not.used.or.pairs/=0)error stop 'empty targets'
 ! A changed displacement tile (single point) must replace the previous plan.
 source=0
 if(all(lo==0))source(1)=cmplx(.7d0,.4d0,8)
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target,action,used,status,pairs,points)
 if(status/=0.or..not.used.or.pairs/=6.or.points/=6)error stop 'single point'
 do j=1,7
 work(:,j)=conjg(source)*target(:,j)
 enddo
 call pencil_transform(n,dims,coords,comm,work(:,1:4),spec(:,1:4),-1,status,spectral_z=.true.)
 call pencil_transform(n,dims,coords,comm,work(:,5:7),spec(:,5:7),-1,status,spectral_z=.true.)
 do j=1,7
 spec(:,j)=spec(:,j)*multiplier
 enddo
 call pencil_transform(n,dims,coords,comm,spec(:,1:4),work(:,1:4),1,status,spectral_z=.true.)
 call pencil_transform(n,dims,coords,comm,spec(:,5:7),work(:,5:7),1,status,spectral_z=.true.)
 do j=1,7
 ref(:,j)=-source*work(:,j)
 enddo
 err=maxval(abs(ref-action));call MPI_Allreduce(err,err_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
 if(err_all>1d-11)error stop 'single point parity'
 call spatial_local_destroy(plan)
 if(rank==0)print *, 'PASS spatial local',np
 call MPI_Finalize(ierr)
end program
