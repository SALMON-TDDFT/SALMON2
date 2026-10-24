program probe
use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
use mpi
use exx_symmetry
use exx_k_exchange
use exx_spinor
use exx_k_orbits,only:apply_spinor_k_orbits
implicit none
type(s_exx_symmetry_map)::m
type(exx_k_kernel)::op,refop
real(8)::ops(3,4,4),ks(3,2),ws(2),h(3),pi,errmax,t0
complex(8)::u(96,2,2),x(96,3,2),ref(96,3,3)
complex(8),allocatable::es(:,:,:),et(:,:,:),out(:,:,:),left(:,:),right(:,:)
integer,allocatable::starts(:),counts(:)
integer::i,j,k,q,e,rank,np,repfirst,nrep,first,nlocal,p,rep,mode
character(20)::arg
call MPI_Init(e);call MPI_Comm_rank(MPI_COMM_WORLD,rank,e);call MPI_Comm_size(MPI_COMM_WORLD,np,e)
if(np>3)call MPI_Abort(MPI_COMM_WORLD,1,e)
call get_command_argument(1,arg)
ops=0;pi=acos(-1d0);ks=0;h=[1d0,2d0,3d0];ws=[1d0/3,2d0/3]
ks(1,2)=-2*pi/(3*4*h(1))
do q=1,4;do i=1,3;ops(i,i,q)=1;enddo;enddo
ops(1,1,2)=-1;ops(1,1,4)=-1
ops(:,4,3)=[.5d0,.5d0,0d0];ops(:,4,4)=ops(:,4,3)
call symmetry_init(m,[4,6,2],h,ks,ws,ops,ops,[3,1,1],e)
if(e/=0)call MPI_Abort(MPI_COMM_WORLD,2,e)
do k=1,2;do j=1,2;do i=1,96
u(i,j,k)=cmplx(sin(real(i+3*j+2*k,8)),cos(real(2*i-j+k,8)),8)/7d0
enddo;enddo;enddo
do k=1,2;do j=1,3;do i=1,96
x(i,j,k)=cmplx(cos(real(i+2*j+k,8)),sin(real(3*i-j+k,8)),8)/8d0
enddo;enddo;enddo
call symmetry_expand(m,1,u,x,es,et,e)
if(e/=0)call MPI_Abort(MPI_COMM_WORLD,3,e)
call exx_k_kernel_init(refop,[4,6,2],[3,1,1],h,m%full_k,.11d0,2,e,fft_layout='strided')
if(e/=0)call MPI_Abort(MPI_COMM_WORLD,4,e)
call exx_spinor_apply_distributed(refop,es,et,ref,[1],[3],0,MPI_COMM_SELF,transpose_local,e)
if(e/=0)call MPI_Abort(MPI_COMM_WORLD,5,e)
allocate(starts(np),counts(np));repfirst=1;nrep=2
if(np>1)then;repfirst=min(rank+1,3);nrep=merge(1,0,rank<2);endif
do p=1,np
if(np==1)then;starts(p)=1;counts(p)=3
else;starts(p)=m%first(min(p,3));counts(p)=0
if(p<=2)counts(p)=m%first(p+1)-m%first(p)
endif
enddo
first=starts(rank+1);nlocal=counts(rank+1)
allocate(out(96,3,nrep),left(48,size(es,2)),right(48,size(es,2)))
call exx_k_kernel_init(op,[4,6,2],[3,1,1],h,m%full_k,.11d0,2,e,phase_start=first, &
phase_count=nlocal,fft_layout=merge('contiguous','strided   ',trim(arg)=='contiguous'))
if(e/=0)call MPI_Abort(MPI_COMM_WORLD,6,e)
if(trim(arg)=='nan'.and.rank==0)x(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
do rep=0,2
call MPI_Barrier(MPI_COMM_WORLD,e);t0=MPI_Wtime()
call apply_spinor_k_orbits(op,m,x(:,:,repfirst:repfirst+nrep-1),out,starts,counts, &
repfirst,rank,MPI_COMM_WORLD,transpose_tiles,fill,e)
if(trim(arg)=='nan')then
if(e==0)call MPI_Abort(MPI_COMM_WORLD,9,e)
print *,'SPINOR_ORBIT_INVALID_PASS',rank
call MPI_Finalize(e)
stop
endif
if(e/=0)call MPI_Abort(MPI_COMM_WORLD,7,e)
errmax=0
do j=1,nrep
q=m%first(repfirst+j-1);errmax=max(errmax,maxval(abs(out(:,:,j)-ref(:,:,q))))
enddo
if(errmax>2d-11)then
print *,'ERROR',rank,errmax
call MPI_Abort(MPI_COMM_WORLD,8,e)
endif
print *,'SPINOR_ORBIT_PASS',rank,rep,errmax,MPI_Wtime()-t0
enddo
call MPI_Finalize(e)
contains
subroutine fill(j,points,density,s,t,columns)
integer,intent(in)::j,points(:),s,t,columns(:)
complex(8),intent(out)::density(:,:)
integer::a,b,q
q=first+j-1
left=es((s-1)*48+1:s*48,:,q);right=es((t-1)*48+1:t*48,:,q)
do b=1,size(left,2)
left(:,b)=left(:,b)*op%phase(:,j);right(:,b)=right(:,b)*op%phase(:,j)
enddo
density=matmul(left(points,:),conjg(transpose(right(columns,:))))
end subroutine
subroutine transpose_tiles(send,recv,count,comm)
complex(8),intent(in)::send(:)
complex(8),intent(out)::recv(:)
integer,intent(in)::count,comm
integer::e
call MPI_Alltoall(send,count,MPI_DOUBLE_COMPLEX,recv,count,MPI_DOUBLE_COMPLEX,comm,e)
end subroutine
subroutine transpose_local(send,recv,count,comm)
complex(8),intent(in)::send(:)
complex(8),intent(out)::recv(:)
integer,intent(in)::count,comm
recv=send
end subroutine
end program
