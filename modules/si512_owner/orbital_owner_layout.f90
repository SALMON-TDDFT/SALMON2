module orbital_owner_layout
  use mpi
  use iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use wannier_pair_owner, only: orbital_owners
  implicit none
  private
  public :: owner_layout,layout_prepare,layout_forward,layout_reverse
  type owner_layout
    integer :: comm=MPI_COMM_NULL,ng=0,norb=0,nrow=0,ncol=0
    logical :: ready=.false.
    integer,allocatable :: ids(:),sc(:),rc(:),sd(:),rd(:),send_row(:),send_col(:),recv_row(:),recv_col(:)
  end type
contains
  subroutine layout_prepare(p,ng,norb,rows,cols,comm,status)
    type(owner_layout),intent(out) :: p
    integer,intent(in) :: ng,norb,rows(:),cols(:),comm
    integer,intent(out) :: status
    integer :: rank,peers,e,bad,i,j,r,k,nlocal,ns,nr,params(2),pmin(2),pmax(2)
    integer,allocatable :: owners(:),cursor(:),slot(:)
    integer(int64),allocatable :: skey(:),rkey(:)
    logical,allocatable :: filled(:,:)
    call MPI_Comm_rank(comm,rank,e);call MPI_Comm_size(comm,peers,e)
    params=[ng,norb]
    call MPI_Allreduce(params,pmin,2,MPI_INTEGER,MPI_MIN,comm,e)
    call MPI_Allreduce(params,pmax,2,MPI_INTEGER,MPI_MAX,comm,e)
    bad=0
    if(any(pmin/=pmax).or.ng<1.or.norb<1)bad=1
    if(int(size(rows),int64)*int(size(cols),int64)>int(huge(ns),int64))bad=1
    if(any(rows<1).or.any(rows>ng).or.any(cols<1).or.any(cols>norb))bad=1
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    if(status/=0)return
    allocate(owners(norb),slot(norb),cursor(peers))
    call orbital_owners(norb,peers,owners)
    nlocal=count(owners==rank)
    if(int(ng,int64)*int(nlocal,int64)>int(huge(nr),int64))bad=1
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    if(status/=0)return
    p%comm=comm;p%ng=ng;p%norb=norb;p%nrow=size(rows);p%ncol=size(cols)
    allocate(p%ids(nlocal),p%sc(peers),p%rc(peers),p%sd(peers),p%rd(peers))
    k=0;slot=0
    do i=1,norb
      if(owners(i)/=rank)cycle
      k=k+1;p%ids(k)=i;slot(i)=k
    enddo
    p%sc=0
    do j=1,size(cols)
      r=owners(cols(j))+1;p%sc(r)=p%sc(r)+size(rows)
    enddo
    call MPI_Alltoall(p%sc,1,MPI_INTEGER,p%rc,1,MPI_INTEGER,comm,e)
    if(sum(int(p%rc,int64))>int(huge(nr),int64))bad=1
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    if(status/=0)return
    p%sd(1)=0;p%rd(1)=0
    do r=2,peers
      p%sd(r)=p%sd(r-1)+p%sc(r-1);p%rd(r)=p%rd(r-1)+p%rc(r-1)
    enddo
    ns=sum(p%sc);nr=sum(p%rc)
    if(nr/=ng*nlocal)bad=1
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    if(status/=0)return
    allocate(p%send_row(ns),p%send_col(ns),p%recv_row(nr),p%recv_col(nr),skey(max(1,ns)),rkey(max(1,nr)))
    cursor=p%sd
    do j=1,size(cols)
      r=owners(cols(j))+1
      do i=1,size(rows)
        cursor(r)=cursor(r)+1;k=cursor(r)
        p%send_row(k)=i;p%send_col(k)=j
        skey(k)=int(rows(i)-1,int64)+int(ng,int64)*int(cols(j)-1,int64)
      enddo
    enddo
    call MPI_Alltoallv(skey,p%sc,p%sd,MPI_INTEGER8,rkey,p%rc,p%rd,MPI_INTEGER8,comm,e)
    allocate(filled(ng,nlocal));filled=.false.
    do k=1,nr
      i=int(modulo(rkey(k),int(ng,int64)))+1;j=int(rkey(k)/int(ng,int64))+1
      if(j<1.or.j>norb)then
        bad=1;cycle
      endif
      if(slot(j)==0)then
        bad=1;cycle
      endif
      j=slot(j)
      if(filled(i,j))bad=1
      filled(i,j)=.true.;p%recv_row(k)=i;p%recv_col(k)=j
    enddo
    if(.not.all(filled))bad=1
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    p%ready=status==0
  end subroutine

  subroutine layout_forward(p,local,owned,status)
    type(owner_layout),intent(in) :: p
    complex(8),intent(in) :: local(:,:)
    complex(8),intent(out) :: owned(:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: sb(:),rb(:)
    integer :: bad,k,i,j,e,ns,nr
    ! Only a successfully prepared collective plan may be passed here.
    bad=0;owned=0
    if(size(local,1)/=p%nrow.or.size(local,2)/=p%ncol)bad=1
    if(size(owned,1)/=p%ng.or.size(owned,2)/=size(p%ids))bad=1
    do j=1,size(local,2)
      do i=1,size(local,1)
        if(.not.ieee_is_finite(real(local(i,j),8)))bad=1
        if(.not.ieee_is_finite(aimag(local(i,j))))bad=1
      enddo
    enddo
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,p%comm,e)
    if(status/=0)return
    ns=sum(p%sc);nr=sum(p%rc);allocate(sb(max(1,ns)),rb(max(1,nr)))
    do k=1,ns
      sb(k)=local(p%send_row(k),p%send_col(k))
    enddo
    call MPI_Alltoallv(sb,p%sc,p%sd,MPI_DOUBLE_COMPLEX,rb,p%rc,p%rd,MPI_DOUBLE_COMPLEX,p%comm,e)
    do k=1,nr
      owned(p%recv_row(k),p%recv_col(k))=rb(k)
    enddo
  end subroutine

  subroutine layout_reverse(p,owned,local,status)
    type(owner_layout),intent(in) :: p
    complex(8),intent(in) :: owned(:,:)
    complex(8),intent(out) :: local(:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: sb(:),rb(:)
    integer :: bad,k,i,j,e,ns,nr
    bad=0;local=0
    if(size(local,1)/=p%nrow.or.size(local,2)/=p%ncol)bad=1
    if(size(owned,1)/=p%ng.or.size(owned,2)/=size(p%ids))bad=1
    do j=1,size(owned,2)
      do i=1,size(owned,1)
        if(.not.ieee_is_finite(real(owned(i,j),8)))bad=1
        if(.not.ieee_is_finite(aimag(owned(i,j))))bad=1
      enddo
    enddo
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,p%comm,e)
    if(status/=0)return
    ns=sum(p%rc);nr=sum(p%sc);allocate(sb(max(1,ns)),rb(max(1,nr)))
    do k=1,ns
      sb(k)=owned(p%recv_row(k),p%recv_col(k))
    enddo
    call MPI_Alltoallv(sb,p%rc,p%rd,MPI_DOUBLE_COMPLEX,rb,p%sc,p%sd,MPI_DOUBLE_COMPLEX,p%comm,e)
    do k=1,nr
      local(p%send_row(k),p%send_col(k))=rb(k)
    enddo
  end subroutine
end module
