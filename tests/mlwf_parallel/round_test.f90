program round_test
use mpi
use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
use communication,only:comm_init,comm_finalize
use exx_distributed_gauge,only:s_exx_gauge,jacobi_matching,gauge_jacobi_round,jacobi_rotation
implicit none
integer::rank,np,e,n,context,desc(9),pr,pc,mr,mc,nb,nr,nc,i,j,k,l,rd,a,b,q,st
integer,external::sys2blacs_handle,numroc
integer,allocatable::mate(:),slot(:),seen(:,:)
complex(8),allocatable::full(:,:,:),ref(:,:,:),u(:,:),ur(:,:),links(:,:,:),col(:),row(:),rot(:,:)
real(8)::err,allerr,orth,theta,weights(6),timings(4)
complex(8)::s
real(8)::c
integer(kind=8)::accepted
integer,parameter::cases(5)=[1,3,8,17,65]
integer::test_comm,world_rank,permutation
complex(8)::blk(2,2,6)
real(8),allocatable::rr(:,:)
type(s_exx_gauge)::g
call comm_init
call MPI_Comm_rank(MPI_COMM_WORLD,world_rank,e);call MPI_Comm_size(MPI_COMM_WORLD,np,e)
do permutation=0,1
 call MPI_Comm_split(MPI_COMM_WORLD,0,merge(world_rank,np-1-world_rank,permutation==0),test_comm,e)
 call MPI_Comm_rank(test_comm,rank,e)
do q=1,5
 n=cases(q);pr=int(sqrt(dble(np)))
 do while(mod(np,pr)/=0);pr=pr-1;enddo
 pc=np/pr;context=sys2blacs_handle(test_comm)
 call blacs_gridinit(context,'R',pr,pc);call blacs_gridinfo(context,pr,pc,mr,mc)
 nb=2;if(n==65)nb=32
 nr=max(1,numroc(n,nb,mr,0,pr));nc=max(1,numroc(n,nb,mc,0,pc))
 call descinit(desc,n,n,nb,nb,0,0,context,nr,st)
 if(st/=0)call MPI_Abort(MPI_COMM_WORLD,1,e)
 g%n=n;g%comm=test_comm
 allocate(g%rows(nr),g%cols(nc),g%matrix(nr,nc),links(nr,nc,6),full(n,n,6),ref(n,n,6), &
 u(n,n),ur(n,n),mate(n),slot(n),seen(n,n),rr(4,(n+1)/2),col(n),row(n))
 g%rows=0;g%cols=0;g%matrix=0;links=0;u=0;seen=0;weights=1d0
 do i=1,n
  if(mod((i-1)/nb,pr)==mr)g%rows(((i-1)/(nb*pr))*nb+mod(i-1,nb)+1)=i
  if(mod((i-1)/nb,pc)==mc)g%cols(((i-1)/(nb*pc))*nb+mod(i-1,nb)+1)=i
  u(i,i)=1
  do j=1,n;do l=1,3
   full(i,j,l)=cmplx(sin(dble(i+3*j+l)),cos(dble(2*i+j+2*l)),8)/n
  enddo;enddo
 enddo
 do l=1,3;full(:,:,l+3)=transpose(conjg(full(:,:,l)));enddo
 ref=full;ur=u
 do j=1,nc;do i=1,nr
  a=g%rows(i);b=g%cols(j)
  if(a>0.and.b>0)then;links(i,j,:)=full(a,b,:);g%matrix(i,j)=u(a,b);endif
 enddo;enddo
 allerr=0d0
 do rd=1,n+mod(n,2)-1
  call jacobi_matching(n,rd,mate,slot)
  rr=0
  do a=1,n
   b=mate(a)
   if(a>=b)cycle
   seen(a,b)=seen(a,b)+1
   do l=1,6
    blk(:,:,l)=reshape([ref(a,a,l),ref(b,a,l),ref(a,b,l),ref(b,b,l)],[2,2])
   enddo
   call jacobi_rotation(blk,weights,n,1d-8,rr(:,slot(a)))
  enddo
  do a=1,n
   b=mate(a);if(a>=b)cycle
   c=rr(1,slot(a));s=cmplx(rr(2,slot(a)),rr(3,slot(a)),8)
   do l=1,6
    col=ref(:,a,l);row=ref(:,b,l)
    ref(:,a,l)=c*col+s*row;ref(:,b,l)=-conjg(s)*col+c*row
    col=ref(a,:,l);row=ref(b,:,l)
    ref(a,:,l)=c*col+conjg(s)*row;ref(b,:,l)=-s*col+c*row
   enddo
   col=ur(:,a);row=ur(:,b);ur(:,a)=c*col+s*row;ur(:,b)=-conjg(s)*col+c*row
  enddo
  call gauge_jacobi_round(g,links,desc,weights,1d-8,rd,st,timings,accepted)
  if(st/=0)call MPI_Abort(MPI_COMM_WORLD,2,e)
  err=0
  do j=1,nc;do i=1,nr
   a=g%rows(i);b=g%cols(j)
   if(a>0.and.b>0)err=max(err,maxval(abs(links(i,j,:)-ref(a,b,:))),abs(g%matrix(i,j)-ur(a,b)))
  enddo;enddo
  call MPI_Allreduce(err,allerr,1,MPI_DOUBLE_PRECISION,MPI_MAX,test_comm,e)
  if(allerr>1d-11)then
   print *,'ROUND_ERROR',n,np,rd,allerr;call MPI_Abort(MPI_COMM_WORLD,3,e)
  endif
 enddo
 do a=1,n;do b=a+1,n
  if(seen(a,b)/=1)call MPI_Abort(MPI_COMM_WORLD,4,e)
 enddo;enddo
 full(:,:,1)=matmul(transpose(conjg(ur)),ur)
 do a=1,n;full(a,a,1)=full(a,a,1)-1;enddo
 orth=maxval(abs(full(:,:,1)))
 if(orth>1d-11)call MPI_Abort(MPI_COMM_WORLD,5,e)
 weights=0d0
 call gauge_jacobi_round(g,links,desc,weights,1d-8,1,st,timings,accepted)
 if(st/=0.or.accepted/=0)call MPI_Abort(MPI_COMM_WORLD,7,e)
 do j=1,nc;do i=1,nr
  a=g%rows(i);b=g%cols(j)
  if(a>0.and.b>0)then
   if(maxval(abs(links(i,j,:)-ref(a,b,:)))>1d-11.or.abs(g%matrix(i,j)-ur(a,b))>1d-11) &
    call MPI_Abort(MPI_COMM_WORLD,8,e)
  endif
 enddo;enddo
 if(world_rank==0)print *,'PARALLEL_ROUND_PASS',n,np,permutation,nb,allerr,orth
 if(n>1)then
  if(rank==np-1)links(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call gauge_jacobi_round(g,links,desc,weights,1d-8,1,st,timings,accepted)
  if(st==0)call MPI_Abort(MPI_COMM_WORLD,6,e)
 endif
 deallocate(g%rows,g%cols,g%matrix,links,full,ref,u,ur,mate,slot,seen,rr,col,row)
 call blacs_gridexit(context)
enddo
call MPI_Comm_free(test_comm,e)
enddo
call comm_finalize
end program
