program spatial_gamma_probe
 use mpi
 use hse_spatial
 implicit none
 integer,parameter :: n(3)=[64,16,16],no=32
 real(8),parameter :: h(3)=.5d0,tol=1d-9
 complex(8),allocatable :: psi(:,:,:),local(:,:,:),base(:,:),rot(:,:),column(:)
 type(spatial_exx_state) :: spatial,partitioned
 integer :: rank,np,err,ng,io,jo,g,x,y,z,center(3),m(3),lo(3),dims(2),coords(2),comm(2),status
 integer :: spatial_rank,orb_rank,comm_r,comm_o,orb_size,first,last,spatial_size
 real(8) :: d(3),theta,spread,gram_error
 complex(8) :: sine,gram(no,no),total(no,no)
 character(20) :: arg
 call MPI_Init(err)
 call MPI_Comm_size(MPI_COMM_WORLD,np,err);call MPI_Comm_rank(MPI_COMM_WORLD,rank,err)
 call get_command_argument(1,arg);read(arg,*)orb_size
 spatial_size=np/orb_size;orb_rank=rank/spatial_size;spatial_rank=modulo(rank,spatial_size)
 call MPI_Comm_split(MPI_COMM_WORLD,orb_rank,spatial_rank,comm_r,err)
 call MPI_Comm_split(MPI_COMM_WORLD,spatial_rank,orb_rank,comm_o,err)
 dims=[spatial_size,1];coords=[spatial_rank,0]
 call MPI_Comm_split(comm_r,coords(2),coords(1),comm(1),err)
 call MPI_Comm_split(comm_r,coords(1),coords(2),comm(2),err)
 m=[n(1),n(2)/dims(1),n(3)];lo=[0,coords(1)*m(2),0];ng=product(n)
 allocate(base(ng,no),rot(no,no),column(ng));rot=0d0
 do io=1,no
   center=[4+8*modulo(io-1,8),4+8*modulo((io-1)/8,2),4+8*((io-1)/16)]
   g=0
   do z=0,n(3)-1;do y=0,n(2)-1;do x=0,n(1)-1
     g=g+1;d=real([x,y,z]-center,8);d=d-n*anint(d/n);d=d*h
     base(g,io)=exp(-sum(d*d)/1.28d0)
   enddo;enddo;enddo
   do jo=1,io-1
     base(:,io)=base(:,io)-sum(conjg(base(:,jo))*base(:,io))*product(h)*base(:,jo)
   enddo
   base(:,io)=base(:,io)/sqrt(sum(abs(base(:,io))**2)*product(h));rot(io,io)=1d0
 enddo
 do io=1,no-1;do jo=io+1,no
   theta=.37d0*sin(real(io*jo,8));sine=sin(theta)*exp(cmplx(0d0,.13d0*(io+jo),8))
   column(:no)=rot(:,io)
   rot(:,io)=cos(theta)*column(:no)+sine*rot(:,jo)
   rot(:,jo)=-conjg(sine)*column(:no)+cos(theta)*rot(:,jo)
 enddo;enddo
 base=matmul(base,rot)
 allocate(psi(product(m),no,1));g=0
 do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
   g=g+1;psi(g,:,1)=base(1+x+lo(1)+n(1)*(y+lo(2)+n(2)*(z+lo(3))),:)
 enddo;enddo;enddo
 call spatial_exx_refresh(spatial,n,h,dims,coords,comm,comm_r,psi,100,tol,status)
 if(rank==0)print *, 'spatial',status,spatial%last_localization_status,spatial%iterations,spatial%gradient,spatial%spread
 if(status/=0.or.spatial%last_localization_status/=0.or.spatial%gradient>=tol)call MPI_Abort(MPI_COMM_WORLD,1,err)
 gram=matmul(conjg(transpose(spatial%source)),spatial%source)*product(h)
 call MPI_Allreduce(gram,total,no*no,MPI_DOUBLE_COMPLEX,MPI_SUM,comm_r,err)
 do io=1,no
   total(io,io)=total(io,io)-1d0
 enddo
 if(maxval(abs(total))>1d-10)call MPI_Abort(MPI_COMM_WORLD,2,err)
 first=no*orb_rank/orb_size+1;last=no*(orb_rank+1)/orb_size
 local=psi(:,first:last,:)
 call spatial_exx_refresh(partitioned,n,h,dims,coords,comm,comm_r,local,100,tol,status,comm_o=comm_o)
 if(rank==0)print *, 'orbital',status,partitioned%last_localization_status,partitioned%iterations,partitioned%gradient,partitioned%spread
 if(status/=0.or.partitioned%last_localization_status/=0.or.partitioned%gradient>=tol)call MPI_Abort(MPI_COMM_WORLD,3,err)
 if(abs(spatial%spread-partitioned%spread)>1d-8)call MPI_Abort(MPI_COMM_WORLD,4,err)
 if(maxval(abs(partitioned%source-spatial%source(:,first:last)))>1d-7)call MPI_Abort(MPI_COMM_WORLD,5,err)
 if(rank==0)print *, 'PASS Gamma localization'
 call MPI_Finalize(err)
end program
