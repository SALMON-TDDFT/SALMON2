program ownerprobe
 use mpi
 use ptcn_native_matrix
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 implicit none
 type(ptcn_owner_matrix)::p
 complex(8),allocatable::u(:,:),h(:,:),r(:,:),x(:,:),full(:,:),hh(:,:),gold(:,:),m(:,:),z(:,:),t(:,:),w(:)
 real(8),allocatable::e(:),rw(:)
 integer::np,rank,n,ng,first,last,i,j,err,s,nocc
 real(8)::diff,total,orth,dv
 character(20)::arg
 call MPI_Init(s);call MPI_Comm_rank(MPI_COMM_WORLD,rank,s);call MPI_Comm_size(MPI_COMM_WORLD,np,s)
 call get_command_argument(1,arg);read(arg,*)n
 ng=2*n+3;first=rank*n/np+1;last=(rank+1)*n/np;nocc=max(1,n-1)
 call ptcn_owner_create(p,ng,n,first,last,MPI_COMM_WORLD,err)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,1,s)
 allocate(u(ng,last-first+1),h(ng,last-first+1),r(ng,last-first+1),x(ng,last-first+1))
 allocate(full(ng,n),hh(ng,n),gold(ng,n),m(n,n),z(n,n),t(n,n),w(max(1,2*n)),e(n),rw(max(1,3*n-2)))
 do j=1,n;do i=1,ng
  full(i,j)=cmplx(sin(.3d0*i*j),cos(.7d0*i+j),8)*.1d0
  if(i==j)full(i,j)=full(i,j)+1
  hh(i,j)=cmplx(cos(.2d0*i*j),sin(.4d0*i+j),8)
 enddo;enddo
 u=full(:,first:last);h=hh(:,first:last)
 call ptcn_owner_project(p,u,h,1d0,nocc,r,err)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,2,s)
 gold=hh
 gold(:,1:nocc)=hh(:,1:nocc)-matmul(full(:,1:nocc),matmul(conjg(transpose(full(:,1:nocc))),hh(:,1:nocc)))
 if(nocc<n)gold(:,nocc+1:n)=hh(:,nocc+1:n)-matmul(full(:,nocc+1:n), &
   matmul(conjg(transpose(full(:,nocc+1:n))),hh(:,nocc+1:n)))
 diff=0
 if(last>=first)diff=maxval(abs(r-gold(:,first:last)))
 call MPI_Allreduce(diff,total,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,s)
 if(total>1d-11)call MPI_Abort(MPI_COMM_WORLD,3,s)
 x=u;call ptcn_owner_lowdin(p,x,1d0,orth,err)
 if(err/=0.or.orth>1d-11)call MPI_Abort(MPI_COMM_WORLD,4,s)
 m=matmul(conjg(transpose(full)),full);call zheev('V','U',n,m,n,e,w,size(w),rw,err)
 z=m;do j=1,n;z(:,j)=z(:,j)/sqrt(e(j));enddo
 t=matmul(z,conjg(transpose(m)));gold=matmul(full,t)
 diff=0;if(last>=first)diff=maxval(abs(x-gold(:,first:last)))
 call MPI_Allreduce(diff,total,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,s)
 if(total>1d-10)call MPI_Abort(MPI_COMM_WORLD,5,s)
 call ptcn_owner_gram(p,x,1d0,orth,err)
 if(err/=0.or.orth>1d-11)call MPI_Abort(MPI_COMM_WORLD,6,s)
 if(rank==0)print *,'NATIVE_OWNER_MATRIX_PASS',n,np,total,orth
 u=full(:,first:last)
 if(rank==np-1.and.last>=first)u(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
 call ptcn_owner_gram(p,u,1d0,orth,err)
 if(err==0)call MPI_Abort(MPI_COMM_WORLD,7,s)
 dv=1d0
 if(rank==0)dv=2d0
 if(np==1)dv=0d0
 call ptcn_owner_gram(p,x,dv,orth,err)
 if(err==0)call MPI_Abort(MPI_COMM_WORLD,8,s)
 if(rank==0)print *,'NATIVE_GRAM_COLLECTIVE_REJECT_PASS',np
 call ptcn_owner_destroy(p);call MPI_Finalize(s)
end program
