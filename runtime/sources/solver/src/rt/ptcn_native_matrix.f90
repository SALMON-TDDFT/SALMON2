! Bounded redistribution from SALMON contiguous orbital owners to block-cyclic waves.
module ptcn_native_matrix
 use mpi
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 use ptcn_distributed_matrix
 implicit none
 private
 public::ptcn_owner_matrix,ptcn_owner_create,ptcn_owner_destroy,ptcn_owner_project,ptcn_owner_lowdin,ptcn_owner_gram
 type ptcn_owner_matrix
  type(ptcn_matrix_grid)::grid
  integer::comm=MPI_COMM_NULL,context=-1,desc(9)=0,np=0,rank=0,first=0,last=-1,ng=0,n=0,cyclic=0
  integer,allocatable::starts(:),ends(:),sc(:),sd(:),rc(:),rd(:)
 end type
contains
 subroutine agree(p,err)
  type(ptcn_owner_matrix),intent(in)::p
  integer,intent(inout)::err
  integer::e,s
  call MPI_Allreduce(err,e,1,MPI_INTEGER,MPI_MAX,p%comm,s);err=e
 end subroutine
 subroutine ptcn_owner_create(p,ng,n,first,last,comm,err)
  type(ptcn_owner_matrix),intent(out)::p
  integer,intent(in)::ng,n,first,last,comm
  integer,intent(out)::err
  integer::s,r,j,k,c,expected
  integer,external::sys2blacs_handle,numroc
  p%comm=comm;p%ng=ng;p%n=n;p%first=first;p%last=last
  call MPI_Comm_size(comm,p%np,s);call MPI_Comm_rank(comm,p%rank,s)
  call ptcn_grid_create(p%grid,ng,n,min(32,max(1,n/max(1,p%np))),comm,err);if(err/=0)return
  allocate(p%starts(0:p%np-1),p%ends(0:p%np-1),p%sc(0:p%np-1),p%sd(0:p%np-1), &
    p%rc(0:p%np-1),p%rd(0:p%np-1),stat=s);err=abs(s);call agree(p,err);if(err/=0)return
  call MPI_Allgather(first,1,MPI_INTEGER,p%starts,1,MPI_INTEGER,comm,s)
  call MPI_Allgather(last,1,MPI_INTEGER,p%ends,1,MPI_INTEGER,comm,s)
  expected=1
  do r=0,p%np-1
   if(p%starts(r)/=expected.or.p%ends(r)<p%starts(r)-1)err=1
   expected=p%ends(r)+1
  enddo
  if(expected/=n+1)err=1
  if(int(ng,8)*int(max(0,last-first+1),8)>huge(1))err=1
  if(int(ng,8)*int((n+p%np-1)/p%np,8)>huge(1))err=1
  call agree(p,err);if(err/=0)return
  p%cyclic=numroc(n,1,p%rank,0,p%np)
  p%context=sys2blacs_handle(comm);call blacs_gridinit(p%context,'R',1,p%np)
  call descinit(p%desc,ng,n,1,1,0,0,p%context,ng,s);err=abs(s)
  p%sc=0;p%rc=0
  do j=first,last
   r=mod(j-1,p%np);p%sc(r)=p%sc(r)+ng
  enddo
  do r=0,p%np-1;do j=p%starts(r),p%ends(r)
   if(mod(j-1,p%np)==p%rank)p%rc(r)=p%rc(r)+ng
  enddo;enddo
  p%sd(0)=0;p%rd(0)=0
  do r=1,p%np-1
   p%sd(r)=p%sd(r-1)+p%sc(r-1);p%rd(r)=p%rd(r-1)+p%rc(r-1)
  enddo
  call agree(p,err)
 end subroutine
 subroutine ptcn_owner_destroy(p)
  type(ptcn_owner_matrix),intent(inout)::p
  if(p%context>=0)call blacs_gridexit(p%context)
  call ptcn_grid_destroy(p%grid)
  p%context=-1
  if(allocated(p%starts))deallocate(p%starts,p%ends,p%sc,p%sd,p%rc,p%rd)
 end subroutine
 subroutine to_grid(p,native,w,err)
  type(ptcn_owner_matrix),intent(in)::p
  complex(8),intent(in)::native(:,:)
  complex(8),intent(out)::w(:,:)
  integer,intent(out)::err
  complex(8),allocatable::send(:),recv(:),cyc(:,:)
  integer::r,j,k,pos,s,e,i
  err=0
  if(any(shape(native)/=[p%ng,p%last-p%first+1]))err=1
  do j=1,size(native,2);do i=1,size(native,1)
   if(.not.ieee_is_finite(real(native(i,j),8)).or..not.ieee_is_finite(aimag(native(i,j))))err=1
  enddo;enddo
  call agree(p,err);if(err/=0)return
  allocate(send(max(1,sum(p%sc))),recv(max(1,sum(p%rc))),cyc(p%ng,max(1,p%cyclic)),stat=s)
  err=abs(s);call agree(p,err);if(err/=0)return
  do r=0,p%np-1
   pos=p%sd(r)
   do j=p%first,p%last
    if(mod(j-1,p%np)/=r)cycle
    send(pos+1:pos+p%ng)=native(:,j-p%first+1);pos=pos+p%ng
   enddo
  enddo
  call MPI_Alltoallv(send,p%sc,p%sd,MPI_DOUBLE_COMPLEX,recv,p%rc,p%rd,MPI_DOUBLE_COMPLEX,p%comm,s)
  cyc=0
  do r=0,p%np-1
   pos=p%rd(r)
   do j=p%starts(r),p%ends(r)
    if(mod(j-1,p%np)/=p%rank)cycle
    k=(j-1)/p%np+1;cyc(:,k)=recv(pos+1:pos+p%ng);pos=pos+p%ng
   enddo
  enddo
  w=0
  call pzgemr2d(p%ng,p%n,cyc,1,1,p%desc,w,1,1,p%grid%dw,p%grid%context)
 end subroutine
 subroutine from_grid(p,w,native,err)
  type(ptcn_owner_matrix),intent(in)::p
  complex(8),intent(in)::w(:,:)
  complex(8),intent(out)::native(:,:)
  integer,intent(out)::err
  complex(8),allocatable::send(:),recv(:),cyc(:,:)
  integer::r,j,k,pos,s
  err=0
  if(any(shape(native)/=[p%ng,p%last-p%first+1]))err=1
  call agree(p,err);if(err/=0)return
  allocate(send(max(1,sum(p%rc))),recv(max(1,sum(p%sc))),cyc(p%ng,max(1,p%cyclic)),stat=s)
  err=abs(s);call agree(p,err);if(err/=0)return
  cyc=0
  call pzgemr2d(p%ng,p%n,w,1,1,p%grid%dw,cyc,1,1,p%desc,p%grid%context)
  do r=0,p%np-1
   pos=p%rd(r)
   do j=p%starts(r),p%ends(r)
    if(mod(j-1,p%np)/=p%rank)cycle
    k=(j-1)/p%np+1;send(pos+1:pos+p%ng)=cyc(:,k);pos=pos+p%ng
   enddo
  enddo
  call MPI_Alltoallv(send,p%rc,p%rd,MPI_DOUBLE_COMPLEX,recv,p%sc,p%sd,MPI_DOUBLE_COMPLEX,p%comm,s)
  do r=0,p%np-1
   pos=p%sd(r)
   do j=p%first,p%last
    if(mod(j-1,p%np)/=r)cycle
    native(:,j-p%first+1)=recv(pos+1:pos+p%ng);pos=pos+p%ng
   enddo
  enddo
 end subroutine
 subroutine ptcn_owner_project(p,u,hu,dv,nocc,r,err,occupations)
  type(ptcn_owner_matrix),intent(in)::p
  complex(8),intent(in)::u(:,:),hu(:,:)
  complex(8),intent(out)::r(:,:)
  real(8),intent(in)::dv
  integer,intent(in)::nocc
  integer,intent(out)::err
  complex(8),allocatable::a(:,:),b(:,:),c(:,:)
  real(8),intent(in),optional::occupations(:)
  integer::s
  allocate(a(p%grid%wr,p%grid%wc),b(p%grid%wr,p%grid%wc),c(p%grid%wr,p%grid%wc),stat=s)
  err=abs(s);call agree(p,err);if(err/=0)return
  call to_grid(p,u,a,err);if(err/=0)return
  call to_grid(p,hu,b,err);if(err/=0)return
  call ptcn_project(p%grid,a,b,dv,nocc,c,err,occupations);if(err/=0)return
  call from_grid(p,c,r,err)
 end subroutine
 subroutine ptcn_owner_lowdin(p,x,dv,orth,err)
  type(ptcn_owner_matrix),intent(in)::p
  complex(8),intent(inout)::x(:,:)
  real(8),intent(in)::dv
  real(8),intent(out)::orth
  integer,intent(out)::err
  complex(8),allocatable::a(:,:)
  integer::s
  allocate(a(p%grid%wr,p%grid%wc),stat=s);err=abs(s);call agree(p,err);if(err/=0)return
  call to_grid(p,x,a,err);if(err/=0)return
  call ptcn_lowdin(p%grid,a,dv,orth,err);if(err/=0)return
  call from_grid(p,a,x,err)
 end subroutine
 subroutine ptcn_owner_gram(p,x,dv,orth,err)
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  type(ptcn_owner_matrix),intent(in)::p
  complex(8),intent(in)::x(:,:)
  real(8),intent(in)::dv
  real(8),intent(out)::orth
  integer,intent(out)::err
  complex(8),allocatable::a(:,:),m(:,:)
  integer::s,i,j,gi,gj
  real(8)::local,lo,hi
  allocate(a(p%grid%wr,p%grid%wc),m(p%grid%mr,p%grid%mc),stat=s)
  err=abs(s);call agree(p,err);if(err/=0)return
  if(.not.ieee_is_finite(dv).or.dv<=0)err=1
  call agree(p,err);if(err/=0)return
  call MPI_Allreduce(dv,lo,1,MPI_DOUBLE_PRECISION,MPI_MIN,p%comm,s)
  call MPI_Allreduce(dv,hi,1,MPI_DOUBLE_PRECISION,MPI_MAX,p%comm,s)
  if(lo/=hi)err=1
  call agree(p,err);if(err/=0)return
  call to_grid(p,x,a,err);if(err/=0)return
  m=0
  call pzgemm('C','N',p%n,p%n,p%ng,cmplx(dv,0d0,8),a,1,1,p%grid%dw,a,1,1,p%grid%dw, &
    (0d0,0d0),m,1,1,p%grid%dm)
  local=0
  do j=1,p%grid%mnc;gj=ptcn_global(j,p%grid%block,p%grid%mycol,p%grid%npcol)
  do i=1,p%grid%mnr;gi=ptcn_global(i,p%grid%block,p%grid%myrow,p%grid%nprow)
   if(gi==gj)m(i,j)=m(i,j)-1
   if(.not.ieee_is_finite(real(m(i,j),8)).or..not.ieee_is_finite(aimag(m(i,j))))err=1
   local=max(local,abs(m(i,j)))
  enddo;enddo
  call agree(p,err);if(err/=0)return
  call MPI_Allreduce(local,orth,1,MPI_DOUBLE_PRECISION,MPI_MAX,p%comm,s)
 end subroutine
end module
