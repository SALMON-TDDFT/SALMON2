#include "config.h"
! Root pivoted QR with one global coefficient array, plus bounded row recovery.
module lcfo_seed
#ifdef USE_MPI
 use mpi
#endif
 use lcfo_dist_rows, only: lcfo_gather_root
 use hse_wannier_gauge, only: gauge_seed_select,gauge_seed_finish
 implicit none
 private
 public :: lcfo_seed_gamma
contains
 subroutine lcfo_seed_gamma(local,counts,comm,u,status,snapshot_unit)
  implicit none
  complex(8),intent(in) :: local(:,:)
  integer,intent(in) :: counts(:),comm
  complex(8),intent(out) :: u(:,:)
  integer,intent(out) :: status
  integer,optional,intent(in) :: snapshot_unit
  complex(8),allocatable :: columns(:,:),overlap(:,:),send(:,:),receive(:,:)
  integer,allocatable :: chosen(:)
  integer :: rank,np,ierr,n,nb,lo,first,width,j,row,tile
  rank=0;np=1
#ifdef USE_MPI
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,np,ierr)
#endif
  if(size(counts)/=np.or.any(counts<0))error stop 'LCFO seed: incompatible row counts'
  n=size(local,2);nb=sum(counts)
  if(size(local,1)/=counts(rank+1).or.any(shape(u)/=[n,n])) &
   error stop 'LCFO seed: incompatible dimensions'
  allocate(chosen(n));chosen=0;status=1
  call lcfo_gather_root(local,counts,comm,columns,adjoint=.true.,snapshot_unit=snapshot_unit)
  if(rank==0)call gauge_seed_select(columns,chosen,status)
  deallocate(columns)
#ifdef USE_MPI
  call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
  if(ierr/=MPI_SUCCESS)error stop 'LCFO seed: QR status broadcast failed'
#endif
  if(status/=0)return
#ifdef USE_MPI
  call MPI_Bcast(chosen,n,MPI_INTEGER,0,comm,ierr)
  if(ierr/=MPI_SUCCESS)error stop 'LCFO seed: pivot broadcast failed'
#endif
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
#ifdef USE_MPI
   ! Exactly one owner contributes to each selected row; all others send zero.
   call MPI_Reduce(send,receive,n*width,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
   if(ierr/=MPI_SUCCESS)error stop 'LCFO seed: selected row reduction failed'
#else
   receive(:,1:width)=send(:,1:width)
#endif
   if(rank==0)overlap(:,first:first+width-1)=receive(:,1:width)
  enddo
  deallocate(send,receive,chosen)
  if(rank==0)call gauge_seed_finish(overlap,nb,u,status)
#ifdef USE_MPI
  call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
  if(ierr/=MPI_SUCCESS)error stop 'LCFO seed: SVD status broadcast failed'
#endif
 end subroutine
end module
