program probe
 use mpi
 use exx_adaptive_support
 implicit none
 integer :: ierr,rank,np,n(3)=[24,16,12],m(3),lo(3),g,x,y,z,j,status
 real(8) :: refcenter(3),mom(6),momglobal(6),outer,globalouter,inside,globalinside
 real(8),allocatable :: d2(:)
 real(8) :: h(3)=[.4d0,.6d0,.8d0],length(3),delta(3),radii(4),loss(4),before(4),after(4),local(4)
 complex(8),allocatable :: source(:,:),saved(:,:),fixed(:,:)
 logical :: protected(4)
 call MPI_Init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 m=[n(1),n(2),n(3)/np];lo=[0,0,rank*m(3)];length=n*h
 allocate(source(product(m),4),saved(product(m),4));g=0
 do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
  g=g+1
  delta=modulo(([x,y,z]+lo)*h-[9.4d0,.3d0,9.3d0]+length/2,length)-length/2
  source(g,1)=exp(-sum((delta/[.65d0,.8d0,.7d0])**2))*cmplx(.6d0,.8d0,8)
  delta=modulo(([x,y,z]+lo)*h-[2.1d0,3.2d0,4.1d0]+length/2,length)-length/2
  source(g,2)=exp(-sum((delta/[.55d0,.9d0,1.1d0])**2))
  source(g,3)=1d0;source(g,4)=0d0
 enddo;enddo;enddo
 saved=source
 call adaptive_source_mask(n,h,lo,m,MPI_COMM_WORLD,source,0d0,radii,loss,protected,status)
 if(status==0.or.any(source/=saved))error stop 'invalid fraction accepted or modified data'
 call adaptive_source_mask(n,h,lo,m,MPI_COMM_WORLD,source,1d0,radii,loss,protected,status)
 if(status/=0.or.any(source/=saved).or.any(loss/=0d0))error stop 'fraction one must be exact'
 call adaptive_source_mask(n,h,lo,m,MPI_COMM_WORLD,source,.999d0,radii,loss,protected,status)
 if(status/=0)error stop 'valid call failed'
 if(any(protected(:2)).or..not.all(protected(3:)))error stop 'protection policy'
 if(any(source(:,3:)/=saved(:,3:)))error stop 'protected data changed'
 local=sum(abs(saved)**2,dim=1)
 call MPI_Allreduce(local,before,4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 local=sum(abs(source)**2,dim=1)
 call MPI_Allreduce(local,after,4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 if(any(after(:2)<.999d0*before(:2)-1d-12))error stop 'retained norm too small'
 if(any(after(:2)>=before(:2)))error stop 'no actual truncation'
 if(any(abs(loss(:3)-(1d0-after(:3)/before(:3)))>1d-12))error stop 'loss mismatch'
 do j=1,2
  if(any(source(:,j)/=saved(:,j).and.source(:,j)/=(0d0,0d0)))error stop 'renormalized source'
 enddo
 ! Independently reconstruct the center and confirm that the boundary shell
 ! cannot be removed while satisfying the requested retained norm.
 allocate(d2(product(m)))
 do j=1,2
  mom=0d0;g=0
  do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
   g=g+1
   delta=([x,y,z]+lo)*h
   mom(:3)=mom(:3)+abs(saved(g,j))**2*cos(2d0*acos(-1d0)*delta/length)
   mom(4:)=mom(4:)+abs(saved(g,j))**2*sin(2d0*acos(-1d0)*delta/length)
  enddo;enddo;enddo
  call MPI_Allreduce(mom,momglobal,6,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  refcenter=atan2(momglobal(4:),momglobal(:3))*length/(2d0*acos(-1d0))
  g=0
  do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
   g=g+1;delta=modulo(([x,y,z]+lo)*h-refcenter+length/2,length)-length/2
   d2(g)=sum(delta**2)
   if(d2(g)<radii(j)**2-1d-10.and.source(g,j)/=saved(g,j))error stop 'hole inside sphere'
   if(d2(g)>radii(j)**2+1d-10.and.source(g,j)/=(0d0,0d0))error stop 'tail outside sphere'
  enddo;enddo;enddo
  outer=maxval(d2,mask=abs(source(:,j))>0d0)
  call MPI_Allreduce(outer,globalouter,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  inside=sum(abs(saved(:,j))**2,mask=d2<globalouter-1d-10)
  call MPI_Allreduce(inside,globalinside,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(globalinside>=.999d0*before(j))error stop 'radius unnecessarily retains an extra shell'
 enddo
 ! Quantile radius should remain small even when the WF straddles cell boundaries.
 if(any(radii(:2)>3d0))error stop 'periodic center or metric incorrect'
 ! Explicit R controls geometry even if the norm target is not met.
 source=saved
 call adaptive_source_mask(n,h,lo,m,MPI_COMM_WORLD,source,.999d0,radii,loss,protected,status,fixed_radius=.7d0)
 if(status/=0.or.any(abs(radii(:2)-.7d0)>1d-12))error stop 'fixed radius not respected'
 if(any(loss(:2)<=.001d0))error stop 'fixture must miss the diagnostic target'
 fixed=source
 source=saved
 call adaptive_source_mask(n,h,lo,m,MPI_COMM_WORLD,source,.5d0,radii,loss,protected,status,fixed_radius=.7d0)
 if(status/=0.or.any(source/=fixed))error stop 'norm target changed fixed support'
 source=saved
 call adaptive_source_mask(n,h,lo,m,MPI_COMM_WORLD,source,1d0,radii,loss,protected,status,fixed_radius=.7d0)
 if(status/=0.or.any(source/=fixed))error stop 'fraction one overrode fixed R'
 if(any(source(:,3:)/=saved(:,3:)))error stop 'fixed R changed protected source'
 local=sum(abs(source)**2,dim=1)
 call MPI_Allreduce(local,after,4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 if(any(abs(loss(:2)-(1d0-after(:2)/before(:2)))>1d-12))error stop 'fixed loss mismatch'
 if(rank==0)write(*,'(8es25.16)')radii,loss
 call MPI_Finalize(ierr)
end program
