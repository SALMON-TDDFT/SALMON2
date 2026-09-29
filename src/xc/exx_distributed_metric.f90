#include "config.h"
! ACE metric tiles are unique over the combined spatial/orbital communicator.
! BLACS contexts live only during construction; ACE copies own plain arrays.
module exx_distributed_metric
  use communication, only: comm_get_groupinfo,comm_summation,comm_bcast,comm_get_max
  use exx_ace, only: s_exx_ace
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: distributed_metric_available,distributed_metric_build,distributed_metric_rotate,distributed_metric_apply
contains
  logical function distributed_metric_available(comm) result(available)
    implicit none
    integer,intent(in) :: comm
    integer :: rank,np
    available=.false.
#ifdef USE_SCALAPACK
    call comm_get_groupinfo(comm,rank,np)
    available=np>1
#endif
  end function

  subroutine distributed_metric_build(ace,u,w,dv,comm_o,comm,counts,first,nonzero,status)
    implicit none
    type(s_exx_ace),intent(inout) :: ace
    complex(8),intent(in) :: u(:,:),w(:,:)
    real(8),intent(in) :: dv
    integer,intent(in) :: comm_o,comm,counts(0:),first,nonzero
    integer,intent(out) :: status
#ifdef USE_SCALAPACK
    integer :: context,nprow,npcol,myrow,mycol,block,desc(9),nr,nc,n,rank,np,rank_o,np_o
    integer :: i,j,k,owner,col,ir,jc,bad
    integer,external :: sys2blacs_handle,numroc
    complex(8),allocatable :: a(:,:),z(:,:),column(:),partial(:),total(:),work(:)
    complex(8) :: query(1)
    real(8),allocatable :: eigen(:),rwork(:)
    real(8) :: rquery(1),local_norm(2),global_norm(2)
    n=sum(counts);status=1
    call comm_get_groupinfo(comm,rank,np)
    call comm_get_groupinfo(comm_o,rank_o,np_o)
    nprow=int(sqrt(dble(np)))
    do while(mod(np,nprow)/=0)
      nprow=nprow-1
    enddo
    npcol=np/nprow
    context=sys2blacs_handle(comm)
    call blacs_gridinit(context,'R',nprow,npcol)
    call blacs_gridinfo(context,nprow,npcol,myrow,mycol)
    block=min(32,max(1,n/max(nprow,npcol)))
    nr=max(1,numroc(n,block,myrow,0,nprow));nc=max(1,numroc(n,block,mycol,0,npcol))
    call descinit(desc,n,n,block,block,0,0,context,nr,bad)
    bad=abs(bad);call comm_get_max(bad,comm)
    if(bad/=0)then
      call blacs_gridexit(context)
      return
    endif
    allocate(ace%metric_rows(nr),ace%metric_cols(nc))
    ace%metric_rows=0;ace%metric_cols=0
    do i=1,n
      if(mod((i-1)/block,nprow)==myrow)then
        ir=((i-1)/(block*nprow))*block+mod(i-1,block)+1
        ace%metric_rows(ir)=i
      endif
      if(mod((i-1)/block,npcol)==mycol)then
        jc=((i-1)/(block*npcol))*block+mod(i-1,block)+1
        ace%metric_cols(jc)=i
      endif
    enddo
    allocate(a(nr,nc));a=0d0
    if(nonzero/=0)then
      allocate(column(size(w,1)),partial(n),total(n))
      k=0
      do owner=0,np_o-1
        do col=1,counts(owner)
          k=k+1
          if(rank_o==owner)column=w(:,col)
          call comm_bcast(column,comm_o,owner)
          partial=0d0
          do i=1,size(u,2)
            partial(first+i-1)=-dv*dot_product(u(:,i),column)
          enddo
          call comm_summation(partial,total,n,comm)
          do j=1,nc
            if(ace%metric_cols(j)/=k)cycle
            do i=1,nr
              if(ace%metric_rows(i)>0)a(i,j)=total(ace%metric_rows(i))
            enddo
          enddo
        enddo
      enddo
      deallocate(column,partial,total)
      allocate(z(nr,nc));z=0d0
      call pztranc(n,n,(1d0,0d0),a,1,1,desc,(0d0,0d0),z,1,1,desc)
      local_norm=[sum(abs(a)**2),sum(abs(a-z)**2)]
      call comm_summation(local_norm,global_norm,2,comm)
      bad=0
      if(.not.finite_matrix(a).or..not.finite_matrix(z))bad=1
      if(.not.ieee_is_finite(global_norm(1)).or..not.ieee_is_finite(global_norm(2)))bad=1
      if(global_norm(1)<=0d0.or.sqrt(global_norm(2))>1d-10*sqrt(global_norm(1)))bad=1
      call comm_get_max(bad,comm)
      if(bad==0)then
        a=.5d0*(a+z)
        allocate(eigen(n))
        call pzheev('V','U',n,a,1,1,desc,eigen,z,1,1,desc,query,-1,rquery,-1,bad)
        bad=abs(bad);call comm_get_max(bad,comm)
        if(bad==0)then
          allocate(work(max(1,int(real(query(1))))),rwork(max(1,int(rquery(1)))))
          call pzheev('V','U',n,a,1,1,desc,eigen,z,1,1,desc,work,size(work),rwork,size(rwork),bad)
          bad=abs(bad);call comm_get_max(bad,comm)
          if(bad==0)then
            do i=1,n
              if(.not.ieee_is_finite(eigen(i)))bad=1
            enddo
            if(eigen(n)<=0d0.or.eigen(1)<=1d-12*eigen(n))bad=1
            if(.not.finite_matrix(z))bad=1
            call comm_get_max(bad,comm)
            if(bad==0)then
              ace%condition=eigen(n)/eigen(1)
              do j=1,nc
                if(ace%metric_cols(j)>0)z(:,j)=z(:,j)/sqrt(eigen(ace%metric_cols(j)))
              enddo
              call move_alloc(z,ace%metric_factor)
            endif
          endif
        endif
      endif
    else
      bad=0
      call move_alloc(a,ace%metric_factor)
    endif
    call blacs_gridexit(context)
    if(bad/=0)return
    ace%metric_distributed=.true.;ace%metric_comm=comm;ace%metric_order=n
    status=0
#else
    status=1
#endif
  end subroutine

  subroutine distributed_metric_rotate(ace,w,comm_o,counts,first,output)
    implicit none
    type(s_exx_ace),intent(in) :: ace
    complex(8),intent(in) :: w(:,:)
    integer,intent(in) :: comm_o,counts(0:),first
    complex(8),intent(out) :: output(:,:)
    complex(8),allocatable :: column(:),partial(:),total(:)
    integer :: rank,np,owner,k,col,i,j
    call comm_get_groupinfo(comm_o,rank,np)
    allocate(column(size(w,1)),partial(ace%metric_order),total(ace%metric_order))
    output=0d0;k=0
    do owner=0,np-1
      do col=1,counts(owner)
        k=k+1
        if(rank==owner)column=w(:,col)
        call comm_bcast(column,comm_o,owner)
        partial=0d0
        do i=1,size(ace%metric_rows)
          if(ace%metric_rows(i)/=k)cycle
          do j=1,size(ace%metric_cols)
            if(ace%metric_cols(j)>0)partial(ace%metric_cols(j))=ace%metric_factor(i,j)
          enddo
        enddo
        call comm_summation(partial,total,size(total),ace%metric_comm)
        do j=1,size(output,2)
          output(:,j)=output(:,j)+column*total(first+j-1)
        enddo
      enddo
    enddo
  end subroutine

  subroutine distributed_metric_apply(ace,input,rotated,output,scratch)
    ! Preserve A (A^H input), including its cancellation behavior near rank loss.
    implicit none
    type(s_exx_ace),intent(in) :: ace
    complex(8),intent(in) :: input(:)
    complex(8),intent(out) :: rotated(:),output(:),scratch(:)
    integer :: i,j,ir,jc,n
    n=ace%metric_order;scratch=0d0
    do j=1,size(ace%metric_cols)
      jc=ace%metric_cols(j)
      if(jc==0)cycle
      do i=1,size(ace%metric_rows)
        ir=ace%metric_rows(i)
        if(ir>0)scratch(jc)=scratch(jc)+conjg(ace%metric_factor(i,j))*input(ir)
      enddo
    enddo
    call comm_summation(scratch,rotated,n,ace%metric_comm)
    scratch=0d0
    do j=1,size(ace%metric_cols)
      jc=ace%metric_cols(j)
      if(jc==0)cycle
      do i=1,size(ace%metric_rows)
        ir=ace%metric_rows(i)
        if(ir>0)scratch(ir)=scratch(ir)+ace%metric_factor(i,j)*rotated(jc)
      enddo
    enddo
    call comm_summation(scratch,output,n,ace%metric_comm)
  end subroutine

  logical function finite_matrix(values) result(finite)
    implicit none
    complex(8),intent(in) :: values(:,:)
    real(8) :: component
    integer :: i,j
    finite=.false.
    do j=1,size(values,2)
      do i=1,size(values,1)
        component=real(values(i,j),8)
        if(.not.ieee_is_finite(component))return
        component=aimag(values(i,j))
        if(.not.ieee_is_finite(component))return
      enddo
    enddo
    finite=.true.
  end function
end module
