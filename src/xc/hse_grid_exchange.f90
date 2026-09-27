#include "config.h"
! Mesh-row exchange with cyclic source ownership and streamed target columns.
! Complete global FFT grids are replicated, but complete source sets are not.
! This correctness-first backend still requires global-grid FFTs and collectives.
module hse_grid_exchange
 use iso_fortran_env, only: int64
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef USE_MPI
 use mpi
#endif
 use hse_wannier, only: s_hse_wannier,wannier_init,wannier_destroy,wannier_apply
 implicit none
 private
 public :: s_hse_grid_exchange,grid_exchange_init,grid_exchange_set_source,grid_exchange_apply,grid_exchange_free
 type s_hse_grid_exchange
  type(s_hse_wannier) :: kernel
  integer :: comm=0,rank=0,np=1,ng=0
  integer,allocatable :: indices(:)
  logical :: active=.false.,source_ready=.false.
 end type
contains
 subroutine grid_exchange_free(op)
  implicit none
  type(s_hse_grid_exchange),intent(inout) :: op
  call wannier_destroy(op%kernel)
  if(allocated(op%indices))deallocate(op%indices)
  op%active=.false.;op%source_ready=.false.;op%ng=0
 end subroutine

 subroutine collective_status(op,bad,status)
  implicit none
  type(s_hse_grid_exchange),intent(in) :: op
  integer,intent(in) :: bad
  integer,intent(out) :: status
#ifdef USE_MPI
  integer :: ierr
  call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,op%comm,ierr)
#else
  status=bad
#endif
 end subroutine

 subroutine grid_exchange_init(op,n,h,omega,indices,comm,fft_batch,fft_measure,status)
  implicit none
  type(s_hse_grid_exchange),intent(inout) :: op
  integer,intent(in) :: n(3),indices(:),comm,fft_batch
  real(8),intent(in) :: h(3),omega
  logical,intent(in) :: fft_measure
  integer,intent(out) :: status
  integer :: bad,i,ierr,settings(5)
  integer(int64) :: ng64
  integer,allocatable :: coverage(:),total(:)
  real(8) :: real_settings(4)
  call grid_exchange_free(op)
  op%comm=comm;op%rank=0;op%np=1
#ifdef USE_MPI
  call MPI_Comm_rank(comm,op%rank,ierr);call MPI_Comm_size(comm,op%np,ierr)
#endif
  ! Check rank agreement before allocating or using any rank-dependent counts.
  settings=[n,fft_batch,merge(1,0,fft_measure)];real_settings=[h,omega]
#ifdef USE_MPI
  call MPI_Bcast(settings,5,MPI_INTEGER,0,comm,ierr)
  call MPI_Bcast(real_settings,4,MPI_DOUBLE_PRECISION,0,comm,ierr)
#endif
  bad=0
  if(any(settings/=[n,fft_batch,merge(1,0,fft_measure)]).or.any(real_settings/=[h,omega]))bad=1
  if(any(n<1).or.any(h<=0d0).or.omega<=0d0.or.fft_batch<1.or.fft_batch>32)bad=1
  do i=1,3
   if(.not.ieee_is_finite(h(i)))bad=1
  enddo
  if(.not.ieee_is_finite(omega))bad=1
  ng64=1_int64
  do i=1,3
   if(n(i)<1)then
    bad=1;exit
   endif
   if(ng64>int(huge(1),int64)/n(i))then
    bad=1;exit
   endif
   ng64=ng64*n(i)
  enddo
  call collective_status(op,bad,status)
  if(status/=0)return
  op%ng=int(ng64)
  if(any(indices<1).or.any(indices>op%ng))bad=1
  call collective_status(op,bad,status)
  if(status/=0)return
  allocate(coverage(op%ng),total(op%ng));coverage=0
  do i=1,size(indices)
   coverage(indices(i))=coverage(indices(i))+1
  enddo
#ifdef USE_MPI
  call MPI_Allreduce(coverage,total,op%ng,MPI_INTEGER,MPI_SUM,comm,ierr)
#else
  total=coverage
#endif
  if(any(total/=1))bad=1
  call collective_status(op,bad,status)
  if(status/=0)return
  deallocate(coverage,total)
  call wannier_init(op%kernel,n,[1,1,1],h,reshape([0d0,0d0,0d0],[3,1]),omega,bad)
  call collective_status(op,bad,status)
  if(status/=0)then
   call grid_exchange_free(op)
   return
  endif
  op%kernel%fft_batch_size=fft_batch;op%kernel%fft_measure=fft_measure
  op%indices=indices;op%active=.true.
 end subroutine

 subroutine grid_exchange_set_source(op,source,status)
  implicit none
  type(s_hse_grid_exchange),intent(inout) :: op
  complex(8),intent(in) :: source(:,:)
  integer,intent(out) :: status
  complex(8),allocatable :: send(:),received(:)
  integer :: ns,ns_root,bad,i,j,owner,nowned,ierr
  status=1
  if(.not.op%active)return
  op%source_ready=.false.
  ns=size(source,2);ns_root=ns
#ifdef USE_MPI
  call MPI_Bcast(ns_root,1,MPI_INTEGER,0,op%comm,ierr)
#endif
  bad=0
  if(ns/=ns_root.or.size(source,1)/=size(op%indices))bad=1
  do j=1,size(source,2);do i=1,size(source,1)
   if(.not.ieee_is_finite(real(source(i,j),8)))bad=1
   if(.not.ieee_is_finite(aimag(source(i,j))))bad=1
  enddo;enddo
  call collective_status(op,bad,status)
  if(status/=0)return
  nowned=ns/op%np
  if(op%rank<mod(ns,op%np))nowned=nowned+1
  if(allocated(op%kernel%source))deallocate(op%kernel%source)
  allocate(op%kernel%source(op%ng,nowned),send(op%ng),received(op%ng))
  j=0
  do i=1,ns
   send=0d0;send(op%indices)=source(:,i);owner=mod(i-1,op%np)
#ifdef USE_MPI
   call MPI_Reduce(send,received,op%ng,MPI_DOUBLE_COMPLEX,MPI_SUM,owner,op%comm,ierr)
#else
   received=send
#endif
   if(op%rank==owner)then
    j=j+1;op%kernel%source(:,j)=received
   endif
  enddo
  op%source_ready=.true.;status=0
 end subroutine

 subroutine grid_exchange_apply(op,target,action,status)
  implicit none
  type(s_hse_grid_exchange),intent(inout) :: op
  complex(8),intent(in) :: target(:,:)
  complex(8),intent(out) :: action(:,:)
  integer,intent(out) :: status
  complex(8),allocatable :: send(:),full(:,:,:),partial(:,:,:),total(:)
  integer :: nt,nt_root,bad,i,j,ierr
  status=1;action=0d0
  if(.not.op%active)return
  nt=size(target,2);nt_root=nt
#ifdef USE_MPI
  call MPI_Bcast(nt_root,1,MPI_INTEGER,0,op%comm,ierr)
#endif
  bad=0
  if(.not.op%source_ready.or.nt/=nt_root.or.size(target,1)/=size(op%indices))bad=1
  if(any(shape(target)/=shape(action)))bad=1
  do j=1,size(target,2);do i=1,size(target,1)
   if(.not.ieee_is_finite(real(target(i,j),8)))bad=1
   if(.not.ieee_is_finite(aimag(target(i,j))))bad=1
  enddo;enddo
  call collective_status(op,bad,status)
  if(status/=0)return
  allocate(send(op%ng),full(op%ng,1,1),partial(op%ng,1,1),total(op%ng))
  do j=1,nt
   send=0d0;send(op%indices)=target(:,j)
#ifdef USE_MPI
   call MPI_Allreduce(send,full,op%ng,MPI_DOUBLE_COMPLEX,MPI_SUM,op%comm,ierr)
#else
   full(:,1,1)=send
#endif
   ! Ranks without sources contribute exact zero and need no worker FFT plans.
   bad=0;partial=0d0
   if(size(op%kernel%source,2)>0)call wannier_apply(op%kernel,full,partial,bad)
   call collective_status(op,bad,status)
   if(status/=0)return
#ifdef USE_MPI
   call MPI_Allreduce(partial,total,op%ng,MPI_DOUBLE_COMPLEX,MPI_SUM,op%comm,ierr)
#else
   total=partial(:,1,1)
#endif
   action(:,j)=total(op%indices)
  enddo
  status=0
 end subroutine
end module
