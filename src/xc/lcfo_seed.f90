#include "config.h"
! Gamma seed QR with bounded row recovery and optional distributed storage.
module lcfo_seed
 use communication, only: comm_get_groupinfo,comm_summation,comm_bcast,comm_get_max
 use lcfo_dist_rows, only: lcfo_gather_root
 use exx_wannier_gauge, only: gauge_seed_select,gauge_seed_finish
 implicit none
 private
 public :: lcfo_seed_gamma
contains
 subroutine lcfo_seed_gamma(local,counts,comm,u,status,snapshot_unit,distributed)
  implicit none
  complex(8),intent(in) :: local(:,:)
  integer,intent(in) :: counts(:),comm
  complex(8),intent(out) :: u(:,:)
  integer,intent(out) :: status
  integer,optional,intent(in) :: snapshot_unit
  complex(8),allocatable :: columns(:,:),overlap(:,:),send(:,:),receive(:,:)
  integer,allocatable :: chosen(:)
  logical,optional,intent(in) :: distributed
  integer :: rank,np,n,nb,lo,first,width,j,row,tile,backend
  call comm_get_groupinfo(comm,rank,np)
  if(size(counts)/=np.or.any(counts<0))error stop 'LCFO seed: incompatible row counts'
  n=size(local,2);nb=sum(counts)
  if(size(local,1)/=counts(rank+1).or.any(shape(u)/=[n,n])) &
   error stop 'LCFO seed: incompatible dimensions'
  allocate(chosen(n));chosen=0;status=1
  backend=0
  if(rank==0.and.present(distributed))then
   if(distributed)backend=1
  endif
  call comm_bcast(backend,comm,0)
  if(backend==1)then
#if defined(USE_MPI) && defined(USE_SCALAPACK)
   if(rank==0)write(*,'(a)')'LCFO seed: distributed pivoted QR'
   call distributed_select(local,counts,comm,chosen,status,snapshot_unit)
#else
   error stop 'LCFO seed: distributed QR requires MPI and ScaLAPACK'
#endif
  else
   call lcfo_gather_root(local,counts,comm,columns,adjoint=.true.,snapshot_unit=snapshot_unit)
   if(rank==0)call gauge_seed_select(columns,chosen,status)
   deallocate(columns)
  endif
  call comm_bcast(status,comm,0)
  if(status/=0)return
  call comm_bcast(chosen,comm,0)
  tile=min(64,n)
  if(n>huge(1)/tile)error stop 'LCFO seed: row tile exceeds MPI count'
  allocate(send(n,tile))
  if(rank==0)then
   allocate(receive(n,tile),overlap(n,n))
  else
   allocate(receive(1,1),overlap(0,0))
  endif
  lo=sum(counts(:rank))
  do first=1,n,tile
   width=min(tile,n-first+1);send=0d0
   do j=1,width
    row=chosen(first+j-1)-lo
    if(row>=1.and.row<=size(local,1))send(:,j)=conjg(local(row,:))
   enddo
   ! Exactly one owner contributes to each selected row; all others send zero.
   call comm_summation(send,receive,n*width,comm,0)
   if(rank==0)overlap(:,first:first+width-1)=receive(:,1:width)
  enddo
  deallocate(send,receive,chosen)
  if(rank==0)call gauge_seed_finish(overlap,nb,u,status)
  call comm_bcast(status,comm,0)
 end subroutine
#if defined(USE_MPI) && defined(USE_SCALAPACK)
 subroutine distributed_select(local,counts,comm,chosen,status,snapshot_unit)
  use mpi, only: MPI_Gatherv,MPI_Send,MPI_Recv,MPI_DOUBLE_COMPLEX,MPI_SUCCESS,MPI_STATUS_IGNORE
  implicit none
  complex(8),intent(in) :: local(:,:)
  integer,intent(in) :: counts(:),comm
  integer,intent(out) :: chosen(:),status
  integer,optional,intent(in) :: snapshot_unit
  integer,parameter :: block=32
  integer :: rank,np,ierr,n,nb,context,nc,desc(9),owner,dest,lo,row,gcol,width,j,jc
  integer :: nw,nrw
  integer,allocatable :: pivots(:),partial(:),displs(:)
  complex(8),allocatable :: a(:,:),buffer(:,:),tau(:),work(:),column(:)
  real(8),allocatable :: rwork(:)
  complex(8) :: query(1)
  real(8) :: rquery(1)
  integer,external :: sys2blacs_handle,numroc
  external :: pzgeqpf
  call comm_get_groupinfo(comm,rank,np)
  n=size(local,2);nb=sum(counts);chosen=0;status=1
  ! Preserve the existing column-major snapshot without allocating global coefficients.
  if(present(snapshot_unit))then
   allocate(displs(np));displs(1)=0
   do j=2,np;displs(j)=displs(j-1)+counts(j-1);enddo
   if(rank==0)then
    allocate(column(nb))
   else
    allocate(column(1))
   endif
   do j=1,n
    call MPI_Gatherv(local(:,j),counts(rank+1),MPI_DOUBLE_COMPLEX,column,counts,displs, &
                     MPI_DOUBLE_COMPLEX,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)error stop 'LCFO seed: snapshot gather failed'
    if(rank==0)write(snapshot_unit)column
   enddo
   deallocate(column,displs)
  endif
  if(n<1.or.n>nb)return
  if(n>huge(1)/block)error stop 'LCFO seed: QR transfer exceeds MPI count'
  context=sys2blacs_handle(comm)
  call blacs_gridinit(context,'R',1,np)
  nc=max(1,numroc(nb,block,rank,0,np))
  call descinit(desc,n,nb,block,block,0,0,context,n,ierr)
  if(ierr/=0)error stop 'LCFO seed: invalid QR descriptor'
  allocate(a(n,nc),buffer(n,block),pivots(nc),tau(nc));a=0d0
  ! Each transfer lies within one cyclic column block and one original owner.
  ! Only source and destination participate; no root or all-rank coefficient copy.
  lo=0
  do owner=0,np-1
   row=1
   do while(row<=counts(owner+1))
    gcol=lo+row;dest=mod((gcol-1)/block,np)
    width=min(block-mod(gcol-1,block),counts(owner+1)-row+1)
    jc=((gcol-1)/(block*np))*block+mod(gcol-1,block)+1
    if(rank==owner)then
     do j=1,width;buffer(:,j)=conjg(local(row+j-1,:));enddo
     if(dest/=owner)then
      call MPI_Send(buffer,n*width,MPI_DOUBLE_COMPLEX,dest,0,comm,ierr)
      if(ierr/=MPI_SUCCESS)error stop 'LCFO seed: QR coefficient send failed'
     endif
    endif
    if(rank==dest)then
     if(dest/=owner)then
      call MPI_Recv(buffer,n*width,MPI_DOUBLE_COMPLEX,owner,0,comm,MPI_STATUS_IGNORE,ierr)
      if(ierr/=MPI_SUCCESS)error stop 'LCFO seed: QR coefficient receive failed'
     endif
     a(:,jc:jc+width-1)=buffer(:,1:width)
    endif
    row=row+width
   enddo
   lo=lo+counts(owner+1)
  enddo
  deallocate(buffer)
  call pzgeqpf(n,nb,a,1,1,desc,pivots,tau,query,-1,rquery,-1,status)
  status=abs(status)
  call comm_get_max(status,comm)
  if(status==0)then
   nw=max(1,int(real(query(1))));nrw=max(1,int(rquery(1)))
   allocate(work(nw),rwork(nrw))
   call pzgeqpf(n,nb,a,1,1,desc,pivots,tau,work,nw,rwork,nrw,status)
   status=abs(status)
   call comm_get_max(status,comm)
  endif
  if(status==0)then
   allocate(partial(n));partial=0
   do j=1,n
    if(mod((j-1)/block,np)/=rank)cycle
    jc=((j-1)/(block*np))*block+mod(j-1,block)+1
    partial(j)=pivots(jc)
   enddo
   call comm_summation(partial,chosen,n,comm)
   if(any(chosen<1).or.any(chosen>nb))status=1
  endif
  call blacs_gridexit(context)
 end subroutine
#endif
end module
