#include "config.h"
program spatial_gamma_probe
 use omp_lib, only: omp_get_max_threads
 use exx_blas_threads, only: scope_entries,scope_active
 use mpi
 use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
 use exx_spatial, only: spatial_exx_state,spatial_exx_refresh,spatial_exx_canonical_source
 use exx_distributed_gauge, only: gauge_tiles_rotate
 implicit none
 integer,parameter :: n(3)=[64,16,16]
 integer :: no=32,bytes,peak_bytes
 real(8),parameter :: h(3)=.5d0,tol=1d-9
 complex(8),allocatable :: psi(:,:,:),local(:,:,:),base(:,:),rot(:,:),column(:),back(:,:)
 type(spatial_exx_state) :: spatial,partitioned,transport_only,copy
 integer :: rank,np,err,ng,io,jo,g,x,y,z,center(3),m(3),lo(3),dims(2),coords(2),comm(2),status
 integer :: spatial_rank,orb_rank,comm_r,comm_o,orb_size,first,last,spatial_size
 real(8) :: d(3),theta,spread,gram_error
 real(8),allocatable :: occupation(:,:)
 complex(8) :: sine
 complex(8),allocatable :: gram(:,:),total(:,:)
 character(20) :: arg
 call MPI_Init(err)
 call MPI_Comm_size(MPI_COMM_WORLD,np,err);call MPI_Comm_rank(MPI_COMM_WORLD,rank,err)
 call get_command_argument(1,arg);read(arg,*)orb_size
 if(command_argument_count()>1)then
   call get_command_argument(2,arg);read(arg,*)no
 endif
 allocate(gram(no,no),total(no,no))
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
 call spatial_exx_refresh(partitioned,n,h,dims,coords,comm,comm_r,local,100,tol,status,comm_o=comm_o,comm_matrix=MPI_COMM_WORLD)
 if(rank==0)print *, 'orbital',status,partitioned%last_localization_status,partitioned%iterations, &
   partitioned%gradient,partitioned%spread
 if(status/=0.or.partitioned%last_localization_status/=0.or.partitioned%gradient>=tol)call MPI_Abort(MPI_COMM_WORLD,3,err)
 if(abs(spatial%spread-partitioned%spread)>1d-8)call MPI_Abort(MPI_COMM_WORLD,4,err)
 if(maxval(abs(partitioned%source-spatial%source(:,first:last)))>1d-7)call MPI_Abort(MPI_COMM_WORLD,5,err)
#ifdef USE_SCALAPACK
 if(np>1)then
   if(.not.allocated(partitioned%gauge_tiles%matrix))call MPI_Abort(MPI_COMM_WORLD,6,err)
   bytes=16*size(partitioned%gauge_tiles%matrix)
   call MPI_Allreduce(bytes,peak_bytes,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err)
   if(peak_bytes>=16*no*no)call MPI_Abort(MPI_COMM_WORLD,7,err)
   if(rank==0)print *, 'GAUGE tile bytes N P replicated/max',no,np,16*no*no,peak_bytes
   allocate(back(size(local,1),size(local,2)))
   call gauge_tiles_rotate(partitioned%gauge_tiles,partitioned%source,comm_r,comm_o,back,status,adjoint=.true.)
   if(status/=0.or.maxval(abs(back-local(:,:,1)))>1d-10)call MPI_Abort(MPI_COMM_WORLD,8,err)
 endif
#endif
 ! A common complex phase must be removed by overlap-polar transport with no sweeps.
 psi=psi*exp(cmplx(0d0,.19d0,8));local=psi(:,first:last,:)
 call spatial_exx_refresh(spatial,n,h,dims,coords,comm,comm_r,psi,0,tol,status)
 if(status/=0)call MPI_Abort(MPI_COMM_WORLD,9,err)
 call spatial_exx_refresh(partitioned,n,h,dims,coords,comm,comm_r,local,0,tol,status, &
   comm_o=comm_o,comm_matrix=MPI_COMM_WORLD)
 if(status/=0.or.abs(partitioned%min_singular-1d0)>1d-10)then
   print *, 'TRANSPORT failure rank/status/minimum',rank,status,partitioned%min_singular
   call MPI_Abort(MPI_COMM_WORLD,10,err)
 endif
 if(maxval(abs(partitioned%source-spatial%source(:,first:last)))>1d-7)call MPI_Abort(MPI_COMM_WORLD,11,err)
 ! A singular transport resets U and schedules a projected-position seed.
 spatial%previous=0d0;partitioned%previous=0d0
 spatial%seed_localized=.true.;partitioned%seed_localized=.true.
 call spatial_exx_refresh(spatial,n,h,dims,coords,comm,comm_r,psi,0,tol,status)
 call spatial_exx_refresh(partitioned,n,h,dims,coords,comm,comm_r,local,0,tol,status, &
   comm_o=comm_o,comm_matrix=MPI_COMM_WORLD)
 if(status/=0.or..not.partitioned%seed_needed)call MPI_Abort(MPI_COMM_WORLD,12,err)
 call spatial_exx_refresh(spatial,n,h,dims,coords,comm,comm_r,psi,100,tol,status)
 call spatial_exx_refresh(partitioned,n,h,dims,coords,comm,comm_r,local,100,tol,status, &
   comm_o=comm_o,comm_matrix=MPI_COMM_WORLD)
 if(status/=0.or.partitioned%seed_needed.or.partitioned%last_localization_status/=0) &
   call MPI_Abort(MPI_COMM_WORLD,13,err)
 if(abs(partitioned%spread-spatial%spread)>1d-8)call MPI_Abort(MPI_COMM_WORLD,14,err)
 if(maxval(abs(partitioned%source-spatial%source(:,first:last)))>1d-7)call MPI_Abort(MPI_COMM_WORLD,15,err)
 ! Failed minimization must retain the accepted transported gauge.
 partitioned%retain_accepted_gauge=.true.;transport_only=partitioned
 do io=1,size(local,2)
   do g=1,size(local,1)
     local(g,io,1)=local(g,io,1)+1d-5*cmplx(sin(dble(g+io)),cos(dble(2*g+io)),8)
   enddo
 enddo
 call spatial_exx_refresh(transport_only,n,h,dims,coords,comm,comm_r,local,0,tol,status, &
   comm_o=comm_o,comm_matrix=MPI_COMM_WORLD)
 if(status/=0)call MPI_Abort(MPI_COMM_WORLD,16,err)
 call spatial_exx_refresh(partitioned,n,h,dims,coords,comm,comm_r,local,1,1d-300,status, &
   comm_o=comm_o,comm_matrix=MPI_COMM_WORLD)
 if(status/=0.or..not.partitioned%retained_gauge.or.partitioned%last_localization_status/=0) &
   call MPI_Abort(MPI_COMM_WORLD,17,err)
 if(maxval(abs(partitioned%source-transport_only%source))>1d-9)call MPI_Abort(MPI_COMM_WORLD,18,err)
 copy=transport_only
 allocate(occupation(size(local,2),1));occupation=2d0
 call spatial_exx_canonical_source(partitioned,local,occupation,comm_r,status,comm_o=comm_o)
 if(status/=0.or.allocated(partitioned%gauge_tiles%matrix))call MPI_Abort(MPI_COMM_WORLD,19,err)
#ifdef USE_SCALAPACK
 if(np>1)then
   call gauge_tiles_rotate(copy%gauge_tiles,copy%source,comm_r,comm_o,back,status,adjoint=.true.)
   if(status/=0.or.maxval(abs(back-local(:,:,1)))>1d-9)call MPI_Abort(MPI_COMM_WORLD,20,err)
 endif
#endif
 if(rank==np-1.and.size(copy%previous)>0)copy%previous(1,1,1)=ieee_value(0d0,ieee_quiet_nan)
 call spatial_exx_refresh(copy,n,h,dims,coords,comm,comm_r,local,0,tol,status, &
   comm_o=comm_o,comm_matrix=MPI_COMM_WORLD)
 if(status==0)call MPI_Abort(MPI_COMM_WORLD,21,err)
#ifdef USE_SCALAPACK
 if(np>1.and.omp_get_max_threads()>1.and.scope_entries==0)call MPI_Abort(MPI_COMM_WORLD,22,err)
#endif
 if(scope_active)call MPI_Abort(MPI_COMM_WORLD,23,err)
 if(rank==0)print *, 'PASS Gamma localization, inverse, transport, reset/seed, retain, copy and invalid state'

 call MPI_Finalize(err)
end program
