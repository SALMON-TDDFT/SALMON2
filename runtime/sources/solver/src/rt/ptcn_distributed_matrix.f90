module ptcn_distributed_matrix
 use mpi
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public::ptcn_matrix_grid,ptcn_grid_create,ptcn_grid_destroy,ptcn_global,ptcn_project,ptcn_lowdin
 type ptcn_matrix_grid
  integer::comm=MPI_COMM_NULL,context=-1,ng=0,n=0,block=0,nprow=0,npcol=0,myrow=0,mycol=0
  integer::dw(9)=0,dm(9)=0,wr=0,wc=0,wnr=0,wnc=0,mr=0,mc=0,mnr=0,mnc=0
 end type
 complex(8),parameter::one=(1d0,0d0),zero=(0d0,0d0)
contains
 integer function ptcn_global(local,block,coord,procs) result(global)
  integer,intent(in)::local,block,coord,procs
  global=((local-1)/block*procs+coord)*block+mod(local-1,block)+1
 end function
 subroutine collective(g,err)
  type(ptcn_matrix_grid),intent(in)::g
  integer,intent(inout)::err
  integer::v,status
  call MPI_Allreduce(err,v,1,MPI_INTEGER,MPI_MAX,g%comm,status);err=v
 end subroutine
 subroutine ptcn_grid_create(g,ng,n,block,comm,err)
  type(ptcn_matrix_grid),intent(out)::g
  integer,intent(in)::ng,n,block,comm
  integer,intent(out)::err
  integer::np,status,params(3),lo(3),hi(3)
  integer,external::sys2blacs_handle,numroc
  g%comm=comm;err=0;params=[ng,n,block]
  call MPI_Allreduce(params,lo,3,MPI_INTEGER,MPI_MIN,comm,status)
  call MPI_Allreduce(params,hi,3,MPI_INTEGER,MPI_MAX,comm,status)
  if(any(lo/=hi).or.minval(lo)<1.or.ng<n.or.maxval(hi)>1000000)err=1
  if(err/=0)return
  call MPI_Comm_size(comm,np,status)
  g%ng=ng;g%n=n;g%block=block;g%nprow=int(sqrt(real(np,8)))
  do while(mod(np,g%nprow)/=0);g%nprow=g%nprow-1;enddo
  g%npcol=np/g%nprow
  g%context=sys2blacs_handle(comm)
  call blacs_gridinit(g%context,'R',g%nprow,g%npcol)
  call blacs_gridinfo(g%context,g%nprow,g%npcol,g%myrow,g%mycol)
  g%wnr=numroc(ng,block,g%myrow,0,g%nprow);g%wnc=numroc(n,block,g%mycol,0,g%npcol)
  g%mnr=numroc(n,block,g%myrow,0,g%nprow);g%mnc=g%wnc
  g%wr=max(1,g%wnr);g%wc=max(1,g%wnc);g%mr=max(1,g%mnr);g%mc=max(1,g%mnc)
  call descinit(g%dw,ng,n,block,block,0,0,g%context,g%wr,status);err=abs(status)
  call descinit(g%dm,n,n,block,block,0,0,g%context,g%mr,status);err=max(err,abs(status))
  call collective(g,err)
 end subroutine
 subroutine ptcn_grid_destroy(g)
  type(ptcn_matrix_grid),intent(inout)::g
  if(g%context>=0)call blacs_gridexit(g%context)
  g%context=-1
 end subroutine
 subroutine validate(g,u,dv,err)
  type(ptcn_matrix_grid),intent(in)::g
  complex(8),intent(in)::u(:,:)
  real(8),intent(in)::dv
  integer,intent(out)::err
  integer::i,j,status
  real(8)::lo,hi
  err=0
  if(any(shape(u)/=[g%wr,g%wc]))err=1
  if(.not.ieee_is_finite(dv).or.dv<=0)err=1
  do j=1,size(u,2);do i=1,size(u,1)
   if(.not.ieee_is_finite(real(u(i,j),8)).or..not.ieee_is_finite(aimag(u(i,j))))err=1
  enddo;enddo
  call collective(g,err)
  if(err/=0)return
  call MPI_Allreduce(dv,lo,1,MPI_DOUBLE_PRECISION,MPI_MIN,g%comm,status)
  call MPI_Allreduce(dv,hi,1,MPI_DOUBLE_PRECISION,MPI_MAX,g%comm,status)
  if(lo/=hi)err=1
 end subroutine
 subroutine matrix_finite(g,a,err)
  type(ptcn_matrix_grid),intent(in)::g
  complex(8),intent(in)::a(:,:)
  integer,intent(out)::err
  integer::i,j
  err=0
  do j=1,g%mnc;do i=1,g%mnr
   if(.not.ieee_is_finite(real(a(i,j),8)).or..not.ieee_is_finite(aimag(a(i,j))))err=1
  enddo;enddo
  call collective(g,err)
 end subroutine
 subroutine ptcn_project(g,u,hu,dv,nocc,r,err,occupations)
  type(ptcn_matrix_grid),intent(in)::g
  complex(8),intent(in)::u(:,:),hu(:,:)
  complex(8),intent(out)::r(:,:)
  real(8),intent(in)::dv
  integer,intent(in)::nocc
  integer,intent(out)::err
  real(8),intent(in),optional::occupations(:)
  real(8),allocatable::flo(:),fhi(:)
  complex(8),allocatable::m(:,:)
  integer::e,b,first,last,nb,lo,hi,status,i,j,gi,gj,has_f,min_f,max_f
  call validate(g,u,dv,err);call validate(g,hu,dv,e);err=max(err,e)
  if(any(shape(r)/=[g%wr,g%wc]))err=1
  call MPI_Allreduce(nocc,lo,1,MPI_INTEGER,MPI_MIN,g%comm,status)
  call MPI_Allreduce(nocc,hi,1,MPI_INTEGER,MPI_MAX,g%comm,status)
  if(lo/=hi.or.lo<1.or.hi>g%n)err=1
  call collective(g,err);if(err/=0)return
  has_f=0
  if(present(occupations))has_f=1
  call MPI_Allreduce(has_f,min_f,1,MPI_INTEGER,MPI_MIN,g%comm,status)
  call MPI_Allreduce(has_f,max_f,1,MPI_INTEGER,MPI_MAX,g%comm,status)
  if(min_f/=max_f)err=1
  call collective(g,err);if(err/=0)return
  if(present(occupations))then
   if(size(occupations)/=g%n)err=1
   do i=1,size(occupations)
    if(.not.ieee_is_finite(occupations(i)))then
     err=1
    else if(occupations(i)<0d0.or.occupations(i)>2d0)then
     err=1
    endif
   enddo
   call collective(g,err);if(err/=0)return
   allocate(flo(g%n),fhi(g%n),stat=e);err=abs(e)
   call collective(g,err);if(err/=0)return
   call MPI_Allreduce(occupations,flo,g%n,MPI_DOUBLE_PRECISION,MPI_MIN,g%comm,status)
   call MPI_Allreduce(occupations,fhi,g%n,MPI_DOUBLE_PRECISION,MPI_MAX,g%comm,status)
   if(any(flo/=fhi))err=1
   call collective(g,err);if(err/=0)return
  endif
  allocate(m(g%mr,g%mc),stat=e);err=abs(e);call collective(g,err);if(err/=0)return
  m=zero;r=hu
  if(present(occupations))then
   ! Flattened spinor rows are included in the Gram product. Only exactly
   ! equal occupation numbers may share a gauge; no near-degeneracy tolerance.
   call pzgemm('C','N',g%n,g%n,g%ng,one*dv,u,1,1,g%dw,hu,1,1,g%dw,zero,m,1,1,g%dm)
   call matrix_finite(g,m,err);if(err/=0)return
   do j=1,g%mnc;gj=ptcn_global(j,g%block,g%mycol,g%npcol)
   do i=1,g%mnr;gi=ptcn_global(i,g%block,g%myrow,g%nprow)
    if(occupations(gi)/=occupations(gj))m(i,j)=zero
   enddo;enddo
   call pzgemm('N','N',g%ng,g%n,g%n,-one,u,1,1,g%dw,m,1,1,g%dm,one,r,1,1,g%dw)
   call validate(g,r,dv,err)
   return
  endif
  do b=1,2
   first=1;last=nocc
   if(b==2)then
    first=nocc+1;last=g%n
   endif
   nb=last-first+1;if(nb==0)cycle
   call pzgemm('C','N',nb,nb,g%ng,one*dv,u,1,first,g%dw,hu,1,first,g%dw,zero,m,first,first,g%dm)
   call matrix_finite(g,m,err);if(err/=0)return
   call pzgemm('N','N',g%ng,nb,nb,-one,u,1,first,g%dw,m,first,first,g%dm,one,r,1,first,g%dw)
  enddo
  call validate(g,r,dv,err)
 end subroutine
 subroutine ptcn_lowdin(g,x,dv,orth,err)
  type(ptcn_matrix_grid),intent(in)::g
  complex(8),intent(inout)::x(:,:)
  real(8),intent(in)::dv
  real(8),intent(out)::orth
  integer,intent(out)::err
  complex(8),allocatable::a(:,:),z(:,:),t(:,:),next(:,:),work(:)
  real(8),allocatable::eigen(:),rwork(:)
  integer,allocatable::iwork(:)
  complex(8)::wq(1)
  real(8)::rq(1),local
  integer::iq(1),lw,lr,li,status,s,i,j,gi,gj
  orth=huge(1d0)
  call validate(g,x,dv,err);if(err/=0)return
  allocate(a(g%mr,g%mc),z(g%mr,g%mc),t(g%mr,g%mc),next(g%wr,g%wc),eigen(g%n),stat=s)
  err=abs(s);call collective(g,err);if(err/=0)return
  a=zero;z=zero;t=zero
  call pzgemm('C','N',g%n,g%n,g%ng,one*dv,x,1,1,g%dw,x,1,1,g%dw,zero,a,1,1,g%dm)
  call matrix_finite(g,a,err);if(err/=0)return
  call pzheevd('V','U',g%n,a,1,1,g%dm,eigen,z,1,1,g%dm,wq,-1,rq,-1,iq,-1,status)
  err=abs(status)
  if(.not.ieee_is_finite(real(wq(1),8)).or..not.ieee_is_finite(rq(1)))err=1
  if(err==0)then
   if(real(wq(1),8)<1.or.real(wq(1),8)>real(huge(1)-1,8).or.rq(1)<1.or.rq(1)>real(huge(1)-1,8).or.iq(1)<1)err=1
  endif
  call collective(g,err);if(err/=0)return
  lw=ceiling(real(wq(1),8));lr=ceiling(rq(1));li=iq(1)
  allocate(work(lw),rwork(lr),iwork(li),stat=s);err=abs(s);call collective(g,err);if(err/=0)return
  call pzheevd('V','U',g%n,a,1,1,g%dm,eigen,z,1,1,g%dm,work,lw,rwork,lr,iwork,li,status)
  err=abs(status)
  do i=1,g%n
   if(.not.ieee_is_finite(eigen(i)))err=1
  enddo
  if(err==0)then
   if(minval(eigen)<=epsilon(1d0)*maxval(eigen))err=1
  endif
  call collective(g,err);if(err/=0)return
  a=z
  do j=1,g%mnc
   gj=ptcn_global(j,g%block,g%mycol,g%npcol)
   a(:,j)=a(:,j)/sqrt(eigen(gj))
  enddo
  call pzgemm('N','C',g%n,g%n,g%n,one,a,1,1,g%dm,z,1,1,g%dm,zero,t,1,1,g%dm)
  next=zero
  call pzgemm('N','N',g%ng,g%n,g%n,one,x,1,1,g%dw,t,1,1,g%dm,zero,next,1,1,g%dw)
  call validate(g,next,dv,err);if(err/=0)return
  call pzgemm('C','N',g%n,g%n,g%ng,one*dv,next,1,1,g%dw,next,1,1,g%dw,zero,a,1,1,g%dm)
  local=0
  do j=1,g%mnc;gj=ptcn_global(j,g%block,g%mycol,g%npcol)
  do i=1,g%mnr;gi=ptcn_global(i,g%block,g%myrow,g%nprow)
   if(gi==gj)a(i,j)=a(i,j)-one
   if(.not.ieee_is_finite(real(a(i,j),8)).or..not.ieee_is_finite(aimag(a(i,j))))err=1
   local=max(local,abs(a(i,j)))
  enddo;enddo
  call collective(g,err);if(err/=0)return
  call MPI_Allreduce(local,orth,1,MPI_DOUBLE_PRECISION,MPI_MAX,g%comm,status)
  x=next
 end subroutine
end module
