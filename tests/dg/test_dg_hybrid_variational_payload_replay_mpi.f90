#include "config.h"
program test_dg_hybrid_variational_payload_replay_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_variational_payload,only:write_dg_hybrid_variational_payload_bundle,&
    read_dg_hybrid_variational_payload_bundle
  implicit none
  integer,parameter::n=4
  integer::comm,rank,nproc,ierr,i,p,nrow,global_count
  integer(int64),allocatable::row_ids(:),read_row_ids(:)
  complex(real64),allocatable::metric(:,:),kinetic(:,:),nonlocal(:,:),interface(:,:),&
    read_metric(:,:),read_kinetic(:,:),read_nonlocal(:,:),read_interface(:,:)
  integer(int64)::basis_fingerprint,metric_fingerprint,interface_fingerprint
  logical::ok
  character(512)::prefix
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call get_command_argument(1,prefix)
  nrow=count([(mod(i-1,nproc)==rank,i=1,n)])
  allocate(row_ids(nrow),metric(nrow,n),kinetic(nrow,n),nonlocal(nrow,n),interface(nrow,n))
  p=0
  do i=1,n
    if(mod(i-1,nproc)/=rank)cycle
    p=p+1;row_ids(p)=int(i,int64)
  enddo
  do p=1,nrow;do i=1,n
    metric(p,i)=cmplx(real(10*row_ids(p)+i,real64),-real(i,real64),real64)
    kinetic(p,i)=2d0*metric(p,i);nonlocal(p,i)=3d0*metric(p,i);interface(p,i)=4d0*metric(p,i)
  enddo;enddo
  call write_dg_hybrid_variational_payload_bundle(comm,trim(prefix),n,row_ids,metric,kinetic,nonlocal,&
    interface,91_int64,92_int64,93_int64,ok,message)
  call require(ok,'payload bundle write failed: '//trim(message))
  call read_dg_hybrid_variational_payload_bundle(comm,trim(prefix),global_count,read_row_ids,read_metric,&
    read_kinetic,read_nonlocal,read_interface,basis_fingerprint,metric_fingerprint,&
    interface_fingerprint,ok,message)
  call require(ok,'payload bundle read failed: '//trim(message))
  call require(global_count==n.and.all(read_row_ids==row_ids),'payload bundle metadata mismatch')
  call require(all(read_metric==metric).and.all(read_kinetic==kinetic).and.&
    all(read_nonlocal==nonlocal).and.all(read_interface==interface),'payload bundle matrix mismatch')
  call require(basis_fingerprint==91_int64.and.metric_fingerprint==92_int64.and.&
    interface_fingerprint==93_int64,'payload bundle fingerprint mismatch')
  call write_dg_hybrid_variational_payload_bundle(comm,trim(prefix),n,row_ids,metric,kinetic,nonlocal,&
    interface,91_int64,92_int64,93_int64,ok,message)
  call require(.not.ok.and.index(message,'exists')>0,'payload bundle silently overwrote existing evidence')
  if(rank==0)write(*,'(a,i0,a)')'PASS variational payload replay on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    integer::bad,global_bad
    bad=merge(0,1,condition);call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;if(rank==0)write(0,'(a)')trim(label);call MPI_Abort(comm,1,ierr);endif
  end subroutine require
end program test_dg_hybrid_variational_payload_replay_mpi
