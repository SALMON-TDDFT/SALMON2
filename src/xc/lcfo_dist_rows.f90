#include "config.h"
! Sparse requests between core-row owners. No replicated global coefficient array.
module lcfo_dist_rows
#ifdef USE_MPI
 use mpi
#endif
 implicit none
 private
 public :: s_lcfo_halo,lcfo_halo_init,lcfo_halo_get,lcfo_halo_sum,lcfo_halo_free,lcfo_gather_root
 type :: s_lcfo_halo
  logical :: ready=.false.
  integer :: comm=0,rank=0,nproc=1,nlocal=0,nselected=0
  integer,allocatable :: need_count(:),give_count(:),need_disp(:),give_disp(:),slots(:),rows(:)
 end type
contains
 subroutine lcfo_halo_init(plan,counts,selected,comm)
  type(s_lcfo_halo),intent(inout) :: plan
  integer,intent(in) :: counts(:),selected(:),comm
  integer,allocatable :: offsets(:),requests(:),received(:),cursor(:)
  integer :: p,j,owner,k,ierr
  call lcfo_halo_free(plan)
  plan%comm=comm
#ifdef USE_MPI
  call MPI_Comm_rank(comm,plan%rank,ierr)
  call MPI_Comm_size(comm,plan%nproc,ierr)
#endif
  if(size(counts)/=plan%nproc.or.any(counts<0))error stop 'LCFO halo: invalid row counts'
  allocate(offsets(plan%nproc+1));offsets(1)=0
  do p=1,plan%nproc;offsets(p+1)=offsets(p)+counts(p);enddo
  if(any(selected<1).or.any(selected>offsets(plan%nproc+1)))error stop 'LCFO halo: invalid request'
  plan%nlocal=counts(plan%rank+1);plan%nselected=size(selected)
  allocate(plan%need_count(plan%nproc),plan%give_count(plan%nproc), &
    plan%need_disp(plan%nproc),plan%give_disp(plan%nproc),cursor(plan%nproc))
  plan%need_count=0
  do j=1,size(selected)
   owner=count(offsets(2:)<selected(j))+1
   plan%need_count(owner)=plan%need_count(owner)+1
  enddo
#ifdef USE_MPI
  call MPI_Alltoall(plan%need_count,1,MPI_INTEGER,plan%give_count,1,MPI_INTEGER,comm,ierr)
  if(ierr/=MPI_SUCCESS)error stop 'LCFO halo: request count exchange failed'
#else
  plan%give_count=plan%need_count
#endif
  plan%need_disp(1)=0;plan%give_disp(1)=0
  do p=2,plan%nproc
   plan%need_disp(p)=plan%need_disp(p-1)+plan%need_count(p-1)
   plan%give_disp(p)=plan%give_disp(p-1)+plan%give_count(p-1)
  enddo
  allocate(requests(max(1,size(selected))),received(max(1,sum(plan%give_count))),plan%slots(size(selected)))
  cursor=plan%need_disp
  do j=1,size(selected)
   owner=count(offsets(2:)<selected(j))+1;cursor(owner)=cursor(owner)+1;k=cursor(owner)
   requests(k)=selected(j);plan%slots(k)=j
  enddo
#ifdef USE_MPI
  call MPI_Alltoallv(requests,plan%need_count,plan%need_disp,MPI_INTEGER, &
    received,plan%give_count,plan%give_disp,MPI_INTEGER,comm,ierr)
  if(ierr/=MPI_SUCCESS)error stop 'LCFO halo: request exchange failed'
#else
  received=requests
#endif
  plan%rows=received(:sum(plan%give_count))-offsets(plan%rank+1)
  if(any(plan%rows<1).or.any(plan%rows>plan%nlocal))error stop 'LCFO halo: request sent to wrong owner'
  plan%ready=.true.
 end subroutine

 subroutine lcfo_halo_get(plan,local,selected)
  type(s_lcfo_halo),intent(in) :: plan
  complex(8),intent(in) :: local(:,:)
  complex(8),allocatable,intent(out) :: selected(:,:)
  complex(8),allocatable :: send(:),recv(:)
  integer :: nc,k,j,ierr
  if(.not.plan%ready.or.size(local,1)/=plan%nlocal)error stop 'LCFO halo: incompatible local rows'
  nc=size(local,2)
  allocate(send(max(1,size(plan%rows)*nc)),recv(max(1,plan%nselected*nc)),selected(plan%nselected,nc))
  do k=1,size(plan%rows);do j=1,nc;send((k-1)*nc+j)=local(plan%rows(k),j);enddo;enddo
#ifdef USE_MPI
  call MPI_Alltoallv(send,plan%give_count*nc,plan%give_disp*nc,MPI_DOUBLE_COMPLEX, &
    recv,plan%need_count*nc,plan%need_disp*nc,MPI_DOUBLE_COMPLEX,plan%comm,ierr)
  if(ierr/=MPI_SUCCESS)error stop 'LCFO halo: value exchange failed'
#else
  recv=send
#endif
  do k=1,plan%nselected;do j=1,nc;selected(plan%slots(k),j)=recv((k-1)*nc+j);enddo;enddo
 end subroutine

 subroutine lcfo_halo_sum(plan,selected,local)
  ! Adjoint of get: sum overlapping fragment contributions on their core owner.
  type(s_lcfo_halo),intent(in) :: plan
  complex(8),intent(in) :: selected(:,:)
  complex(8),allocatable,intent(out) :: local(:,:)
  complex(8),allocatable :: send(:),recv(:)
  integer :: nc,k,j,ierr
  if(.not.plan%ready.or.size(selected,1)/=plan%nselected)error stop 'LCFO halo: incompatible selected rows'
  nc=size(selected,2)
  allocate(send(max(1,plan%nselected*nc)),recv(max(1,size(plan%rows)*nc)),local(plan%nlocal,nc));local=0d0
  do k=1,plan%nselected;do j=1,nc;send((k-1)*nc+j)=selected(plan%slots(k),j);enddo;enddo
#ifdef USE_MPI
  call MPI_Alltoallv(send,plan%need_count*nc,plan%need_disp*nc,MPI_DOUBLE_COMPLEX, &
    recv,plan%give_count*nc,plan%give_disp*nc,MPI_DOUBLE_COMPLEX,plan%comm,ierr)
  if(ierr/=MPI_SUCCESS)error stop 'LCFO halo: contribution exchange failed'
#else
  recv=send
#endif
  do k=1,size(plan%rows);do j=1,nc
   local(plan%rows(k),j)=local(plan%rows(k),j)+recv((k-1)*nc+j)
  enddo;enddo
 end subroutine

 subroutine lcfo_halo_free(plan)
  type(s_lcfo_halo),intent(inout) :: plan
  type(s_lcfo_halo) :: empty
  plan=empty
 end subroutine

 subroutine lcfo_gather_root(local,counts,comm,global)
  ! One-time MLWF seeding/diagnostic output only. Non-root receives no global rows.
  complex(8),intent(in) :: local(:,:)
  integer,intent(in) :: counts(:),comm
  complex(8),allocatable,intent(out) :: global(:,:)
  complex(8),allocatable :: buffer(:)
  integer,allocatable :: disps(:)
  integer :: rank,np,ierr,p,nc,lo
  rank=0;np=1
#ifdef USE_MPI
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,np,ierr)
#endif
  if(size(counts)/=np.or.size(local,1)/=counts(rank+1))error stop 'LCFO gather: incompatible rows'
  nc=size(local,2);allocate(disps(np));disps(1)=0
  do p=2,np;disps(p)=disps(p-1)+counts(p-1)*nc;enddo
  if(rank==0)then
   allocate(buffer(max(1,sum(counts)*nc)),global(sum(counts),nc))
  else
   allocate(buffer(1),global(0,0))
  endif
#ifdef USE_MPI
  call MPI_Gatherv(local,size(local),MPI_DOUBLE_COMPLEX,buffer,counts*nc,disps,MPI_DOUBLE_COMPLEX,0,comm,ierr)
  if(ierr/=MPI_SUCCESS)error stop 'LCFO gather: initial gather failed'
#else
  buffer=reshape(local,[size(local)])
#endif
  if(rank==0)then
   lo=0
   do p=1,np
    global(lo+1:lo+counts(p),:)=reshape(buffer(disps(p)+1:disps(p)+counts(p)*nc),[counts(p),nc])
    lo=lo+counts(p)
   enddo
  endif
 end subroutine
end module
