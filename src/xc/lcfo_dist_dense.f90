#include "config.h"
! Distributed orbital metrics; local-row ACE factors and coefficient transport.
! Orbital-space transforms remain replicated. Large decompositions use ScaLAPACK.
module lcfo_dist_dense
 use communication, only: comm_get_groupinfo,comm_summation,comm_bcast,comm_get_max
 use exx_ace,only:s_exx_ace
 use,intrinsic :: ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public :: lcfo_distributed_polar,lcfo_distributed_ace_build
 integer,parameter :: parallel_threshold=128
#ifdef USE_SCALAPACK
 integer,save :: saved_comm=-1,context=-1,nprow,npcol,myrow,mycol
#endif
 logical,save :: reported_parallel=.false.,reported_serial=.false.
  ! Rank-specific loops keep IEEE inquiries scalar on Fujitsu compilers.
  private :: salmon_all_finite,finite_real_1d,finite_real_2d
  interface salmon_all_finite
    module procedure finite_real_1d,finite_real_2d
  end interface
contains
 subroutine group_info(comm,rank,np)
  implicit none
  integer,intent(in) :: comm
  integer,intent(out) :: rank,np
  call comm_get_groupinfo(comm,rank,np)
 end subroutine
 subroutine sum_matrix(local,total,comm)
  implicit none
  complex(8),intent(in) :: local(:,:)
  complex(8),intent(out) :: total(:,:)
  integer,intent(in) :: comm
  call comm_summation(local,total,size(local),comm)
 end subroutine
 subroutine share_solution(transform,e,condition,status,comm)
  implicit none
  complex(8),intent(inout) :: transform(:,:)
  real(8),intent(inout) :: e,condition
  integer,intent(inout) :: status
  integer,intent(in) :: comm
  call comm_bcast(status,comm,0)
  call comm_bcast(e,comm,0)
  call comm_bcast(condition,comm,0)
  if(status==0)call comm_bcast(transform,comm,0)
 end subroutine

 subroutine lcfo_distributed_polar(current,previous,comm,rotation,minimum,status)
  implicit none
  complex(8),intent(in) :: current(:,:),previous(:,:)
  integer,intent(in) :: comm
  complex(8),intent(out) :: rotation(:,:)
  real(8),intent(out) :: minimum
  integer,intent(out) :: status
  complex(8),allocatable :: local(:,:),metric(:,:)
  real(8) :: condition
  integer :: n
  n=size(current,2);status=1;minimum=0d0;rotation=0d0
  if(any(shape(previous)/=shape(current)).or.any(shape(rotation)/=[n,n]))return
  allocate(local(n,n),metric(n,n))
  local=matmul(conjg(transpose(current)),previous)
  call sum_matrix(local,metric,comm)
  if(.not.salmon_all_finite(real(metric)).or..not.salmon_all_finite(aimag(metric)))return
  call solve_matrix(metric,comm,.true.,rotation,minimum,condition,status)
  if(status==0.and.minimum<1d-8)status=1
 end subroutine

 subroutine lcfo_distributed_ace_build(ace,c,w,dv,comm,status)
  implicit none
  type(s_exx_ace),intent(inout) :: ace
  complex(8),intent(in) :: c(:,:),w(:,:)
  real(8),intent(in) :: dv
  integer,intent(in) :: comm
  integer,intent(out) :: status
  complex(8),allocatable :: local(:,:),metric(:,:),transform(:,:)
  real(8) :: scale,smallest,maximum(1),local_max(1)
  integer :: n
  status=1
  if(allocated(ace%factors))deallocate(ace%factors)
  if(any(shape(c)/=shape(w)).or.dv<=0d0.or..not.ieee_is_finite(dv))return
  n=size(c,2);if(n<1)return
  allocate(local(n,n),metric(n,n),transform(n,n))
  local=-matmul(conjg(transpose(c)),w)*dv
  call sum_matrix(local,metric,comm)
  if(.not.salmon_all_finite(real(metric)).or..not.salmon_all_finite(aimag(metric)))return
  local_max=0d0;if(size(w)>0)local_max=maxval(abs(w))
  call comm_get_max(local_max,maximum,1,comm)
  ace%dv=dv;ace%condition=0d0
  if(maximum(1)==0d0)then
   allocate(ace%factors(size(c,1),n,1));ace%factors=0d0;status=0;return
  endif
  scale=sqrt(sum(abs(metric)**2))
  if(scale==0d0.or.sqrt(sum(abs(metric-conjg(transpose(metric)))**2))>1d-10*scale)return
  metric=.5d0*(metric+conjg(transpose(metric)))
  call solve_matrix(metric,comm,.false.,transform,smallest,ace%condition,status)
  if(status/=0)return
  allocate(ace%factors(size(c,1),n,1));ace%factors(:,:,1)=matmul(w,transform)
 end subroutine

 subroutine solve_matrix(matrix,comm,polar,transform,minimum,condition,status)
  implicit none
  complex(8),intent(inout) :: matrix(:,:)
  integer,intent(in) :: comm
  logical,intent(in) :: polar
  complex(8),intent(out) :: transform(:,:)
  real(8),intent(out) :: minimum,condition
  integer,intent(out) :: status
  complex(8),allocatable :: left(:,:),right(:,:),work(:)
  complex(8) :: query(1)
  real(8),allocatable :: e(:),rwork(:)
  integer :: n,rank,np,j
  external :: zgesvd,zheev
  n=size(matrix,1);call group_info(comm,rank,np)
  status=0;minimum=0d0;condition=0d0;transform=0d0
#ifdef USE_SCALAPACK
  if(n>=parallel_threshold.and.np>1)then
   if(rank==0.and..not.reported_parallel)write(*,'(a,2i8)')'LCFO dense ScaLAPACK dimension/ranks:',n,np
   reported_parallel=.true.
   call parallel_solve(matrix,comm,polar,transform,minimum,condition,status)
   return
  endif
#endif
  if(rank==0)then
   if(.not.reported_serial)write(*,'(a,i8)')'LCFO dense root LAPACK dimension:',n
   reported_serial=.true.
   allocate(e(n),rwork(max(1,5*n)))
   if(polar)then
    allocate(left(n,n),right(n,n))
    call zgesvd('A','A',n,n,matrix,n,e,left,n,right,n,query,-1,rwork,status)
    if(status==0)then
     allocate(work(max(1,int(real(query(1))))))
     call zgesvd('A','A',n,n,matrix,n,e,left,n,right,n,work,size(work),rwork,status)
     if(status==0)then
      minimum=minval(e);transform=matmul(left,right)
     endif
    endif
   else
    call zheev('V','U',n,matrix,n,e,query,-1,rwork,status)
    if(status==0)then
     allocate(work(max(1,int(real(query(1))))))
     call zheev('V','U',n,matrix,n,e,work,size(work),rwork,status)
     if(status==0)then
      if(e(n)<=0d0.or.e(1)<=1d-12*e(n))then
       status=1
      else
       minimum=e(1);condition=e(n)/e(1)
       do j=1,n;transform(:,j)=matrix(:,j)/sqrt(e(j));enddo
      endif
     endif
    endif
   endif
  endif
  call share_solution(transform,minimum,condition,status,comm)
 end subroutine

#ifdef USE_SCALAPACK
 subroutine setup_grid(comm,n,desc,nr,nc,block)
  implicit none
  integer,intent(in) :: comm,n
  integer,intent(out) :: desc(9),nr,nc,block
  integer :: rank,np,ierr
  integer,external :: sys2blacs_handle,numroc
  call group_info(comm,rank,np)
  if(saved_comm/=comm)then
   if(context>=0)call blacs_gridexit(context)
   context=sys2blacs_handle(comm)
   nprow=int(sqrt(dble(np)))
   do while(mod(np,nprow)/=0);nprow=nprow-1;enddo
   npcol=np/nprow
   call blacs_gridinit(context,'R',nprow,npcol)
   call blacs_gridinfo(context,nprow,npcol,myrow,mycol)
   saved_comm=comm
  endif
  block=min(32,max(1,n/max(nprow,npcol)))
  nr=max(1,numroc(n,block,myrow,0,nprow));nc=max(1,numroc(n,block,mycol,0,npcol))
  call descinit(desc,n,n,block,block,0,0,context,nr,ierr)
  if(ierr/=0)error stop 'LCFO dense: invalid ScaLAPACK descriptor'
 end subroutine
 subroutine pack_matrix(global,local,block)
  implicit none
  complex(8),intent(in) :: global(:,:)
  complex(8),intent(out) :: local(:,:)
  integer,intent(in) :: block
  integer :: i,j,ir,jc
  local=0d0
  do j=1,size(global,2)
   if(mod((j-1)/block,npcol)/=mycol)cycle
   jc=((j-1)/(block*npcol))*block+mod(j-1,block)+1
   do i=1,size(global,1)
    if(mod((i-1)/block,nprow)/=myrow)cycle
    ir=((i-1)/(block*nprow))*block+mod(i-1,block)+1
    local(ir,jc)=global(i,j)
   enddo
  enddo
 end subroutine
 subroutine collect_matrix(local,global,block,comm)
  implicit none
  complex(8),intent(in) :: local(:,:)
  complex(8),intent(out) :: global(:,:)
  integer,intent(in) :: block,comm
  complex(8),allocatable :: contribution(:,:)
  integer :: i,j,ir,jc
  allocate(contribution(size(global,1),size(global,2)));contribution=0d0
  do j=1,size(global,2)
   if(mod((j-1)/block,npcol)/=mycol)cycle
   jc=((j-1)/(block*npcol))*block+mod(j-1,block)+1
   do i=1,size(global,1)
    if(mod((i-1)/block,nprow)/=myrow)cycle
    ir=((i-1)/(block*nprow))*block+mod(i-1,block)+1
    contribution(i,j)=local(ir,jc)
   enddo
  enddo
  call sum_matrix(contribution,global,comm)
 end subroutine
 subroutine agree_status(status,comm)
  implicit none
  integer,intent(inout) :: status
  integer,intent(in) :: comm
  status=abs(status)
  call comm_get_max(status,comm)
 end subroutine
 subroutine parallel_solve(matrix,comm,polar,transform,minimum,condition,status)
  implicit none
  complex(8),intent(in) :: matrix(:,:)
  integer,intent(in) :: comm
  logical,intent(in) :: polar
  complex(8),intent(out) :: transform(:,:)
  real(8),intent(out) :: minimum,condition
  integer,intent(out) :: status
  complex(8),allocatable :: a(:,:),left(:,:),right(:,:),product(:,:),work(:)
  complex(8) :: query(1)
  real(8),allocatable :: e(:),rwork(:)
  real(8) :: rquery(1)
  integer :: desc(9),nr,nc,block,n,j,jc
  external :: pzgesvd,pzheev,pzgemm
  n=size(matrix,1);status=0;minimum=0d0;condition=0d0;transform=0d0
  call setup_grid(comm,n,desc,nr,nc,block)
  allocate(a(nr,nc),left(nr,nc),e(n));call pack_matrix(matrix,a,block)
  if(polar)then
   allocate(right(nr,nc),product(nr,nc),rwork(1+4*n))
   call pzgesvd('V','V',n,n,a,1,1,desc,e,left,1,1,desc,right,1,1,desc,query,-1,rwork,status)
   call agree_status(status,comm)
   if(status==0)then
    allocate(work(max(1,int(real(query(1))))))
    call pzgesvd('V','V',n,n,a,1,1,desc,e,left,1,1,desc,right,1,1,desc,work,size(work),rwork,status)
   endif
  else
   call pzheev('V','U',n,a,1,1,desc,e,left,1,1,desc,query,-1,rquery,-1,status)
   call agree_status(status,comm)
   if(status==0)then
    allocate(work(max(1,int(real(query(1))))),rwork(max(1,int(rquery(1)))))
    call pzheev('V','U',n,a,1,1,desc,e,left,1,1,desc,work,size(work),rwork,size(rwork),status)
   endif
  endif
  call agree_status(status,comm)
  if(status/=0)return
  if(.not.salmon_all_finite(e))then
   status=1;return
  endif
  if(polar)then
   minimum=minval(e)
   call pzgemm('N','N',n,n,n,(1d0,0d0),left,1,1,desc,right,1,1,desc,(0d0,0d0),product,1,1,desc)
   call collect_matrix(product,transform,block,comm)
  else
   if(e(n)<=0d0.or.e(1)<=1d-12*e(n))then
    status=1;return
   endif
   minimum=e(1);condition=e(n)/e(1)
   do j=1,n
    if(mod((j-1)/block,npcol)/=mycol)cycle
    jc=((j-1)/(block*npcol))*block+mod(j-1,block)+1
    left(:,jc)=left(:,jc)/sqrt(e(j))
   enddo
   call collect_matrix(left,transform,block,comm)
  endif
 end subroutine
#endif

  pure logical function finite_real_1d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:)
    real(8) :: value
    integer :: i
    finite=.false.
    do i=1,size(values,1)
      value=values(i)
      if(.not.ieee_is_finite(value))return
    enddo
    finite=.true.
  end function

  pure logical function finite_real_2d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:,:)
    real(8) :: value
    integer :: i,j
    finite=.false.
    do j=1,size(values,2)
      do i=1,size(values,1)
        value=values(i,j)
        if(.not.ieee_is_finite(value))return
      enddo
    enddo
    finite=.true.
  end function

end module
