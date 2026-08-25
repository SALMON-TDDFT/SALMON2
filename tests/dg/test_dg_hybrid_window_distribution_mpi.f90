#include "config.h"
program test_dg_hybrid_window_distribution_mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_window_distribution,only:prepare_dg_hybrid_window_distribution,&
    redistribute_dg_hybrid_fragment_windows
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  integer::comm,rank,nproc,ierr,nowned,i,p
  integer,allocatable::fragment_ids(:),core_fragment_ids(:),row_action(:,:),fragment_action(:,:)
  integer(int64),allocatable::box_ids(:),core_ids(:),request_ids(:)
  real(real64),allocatable::box_windows(:,:),raw_windows(:,:)
  integer(int64)::fingerprint,workspace
  logical::ok,values_ok
  character(256)::message
#ifdef USE_MPI
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
#else
  comm=0;rank=0;nproc=1
#endif
  nowned=count([(mod(i-1,nproc)==rank,i=1,2)])
  allocate(fragment_ids(nowned),box_ids(4),box_windows(nowned,4),row_action(4,2))
  p=0
  do i=1,2
    if(mod(i-1,nproc)/=rank)cycle
    p=p+1;fragment_ids(p)=i
  enddo
  box_ids=[1_int64,2_int64,3_int64,4_int64]
  do p=1,nowned;do i=1,4
    box_windows(p,i)=real(10*fragment_ids(p)+i,real64)
  enddo;enddo
  row_action(:,1)=[1,2,3,4]
  row_action(:,2)=[4,3,2,1]
  if(rank==0)then
    allocate(core_ids(2),source=[1_int64,2_int64])
    allocate(core_fragment_ids(2),source=[1,1])
  elseif(rank==1)then
    allocate(core_ids(2),source=[3_int64,4_int64])
    allocate(core_fragment_ids(2),source=[2,2])
  else
    allocate(core_ids(0))
    allocate(core_fragment_ids(0))
  endif
  if(nproc==1)then
    deallocate(core_ids,core_fragment_ids);allocate(core_ids(4),source=box_ids)
    allocate(core_fragment_ids(4),source=[1,1,2,2])
  endif
  call prepare_dg_hybrid_window_distribution(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,core_fragment_ids,&
    row_action,raw_windows,fragment_action,workspace,fingerprint,ok,message)
  call require(ok,'valid window distribution rejected: '//trim(message))
  values_ok=.true.
  do p=1,size(core_ids);do i=1,2
    values_ok=values_ok.and.abs(raw_windows(i,p)-real(10*i+core_ids(p),real64))<1d-12
  enddo;enddo
  call require(values_ok,'distributed raw window mismatch')
  call require(all(fragment_action(:,1)==[1,2]),'identity fragment action mismatch')
  call require(all(fragment_action(:,2)==[2,1]),'reversal fragment action mismatch')
  call require(workspace<=80_int64,'window distribution workspace is not bounded')
  allocate(request_ids(3),source=[4_int64,1_int64,3_int64])
  call redistribute_dg_hybrid_fragment_windows(comm,4,2,fragment_ids,box_ids,box_windows,request_ids,&
    raw_windows,workspace,fingerprint,ok,message)
  call require(ok,'buffer window redistribution failed: '//trim(message))
  values_ok=.true.
  do p=1,size(request_ids);do i=1,2
    values_ok=values_ok.and.abs(raw_windows(i,p)-real(10*i+request_ids(p),real64))<1d-12
  enddo;enddo
  call require(values_ok,'buffer request window mismatch')
  call require(workspace<=48_int64,'buffer window request workspace is not bounded')
  if(rank==0.and.size(fragment_ids)>0)fragment_ids(1)=2
  call prepare_dg_hybrid_window_distribution(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,core_fragment_ids,&
    row_action,raw_windows,fragment_action,workspace,fingerprint,ok,message)
  call require(.not.ok,'duplicate fragment ownership must fail collectively')
  if(rank==0)write(*,'(a,i0,a,i0)')'WINDOW_DISTRIBUTION ranks=',nproc,' fingerprint=',fingerprint
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid window distribution on ',nproc,' ranks'
#ifdef USE_MPI
  call MPI_Finalize(ierr)
#endif
contains
  subroutine require(condition,text)
    logical,intent(in)::condition
    character(*),intent(in)::text
    integer::bad,global_bad
    bad=merge(0,1,condition)
#ifdef USE_MPI
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
#else
    global_bad=bad
#endif
    if(global_bad/=0)then
      if(rank==0)write(0,'(a)')trim(text)
#ifdef USE_MPI
      call MPI_Abort(comm,1,ierr)
#else
      error stop 1
#endif
    endif
  end subroutine require
end program test_dg_hybrid_window_distribution_mpi
