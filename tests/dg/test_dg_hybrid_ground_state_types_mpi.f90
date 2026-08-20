#include "config.h"
program test_dg_hybrid_ground_state_types_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_ground_state_types,only:s_dg_hybrid_ground_state,validate_dg_hybrid_ground_state
  implicit none
  integer,parameter::n=4,noccupied=2
  integer::comm,rank,nproc,ierr,nowned,row,i,position
  integer(int64),allocatable::row_ids(:),duplicate_ids(:)
  complex(real64),allocatable::coefficients(:,:)
  real(real64)::occupations(noccupied),eigenvalues(noccupied)
  type(s_dg_hybrid_ground_state)::state
  integer(int64)::fingerprint,reference_fingerprint,workspace
  logical::ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nowned=count([(mod(row-1,nproc)==rank,row=1,n)])
  allocate(row_ids(nowned),coefficients(nowned,noccupied));position=0
  do row=n,1,-1
    if(mod(row-1,nproc)/=rank)cycle
    position=position+1;row_ids(position)=row
    coefficients(position,1)=cmplx(0.1d0*row,0.03d0*row,real64)
    coefficients(position,2)=cmplx(-0.04d0*row,0.07d0*row,real64)
  enddo
  occupations=[2d0,2d0];eigenvalues=[-0.7d0,-0.2d0]
  call validate_dg_hybrid_ground_state(comm,n,noccupied,row_ids,coefficients,occupations,eigenvalues,4d0,&
    1101_int64,2202_int64,3303_int64,4404_int64,1d-12,state,workspace,fingerprint,ok,message)
  call require(ok,trim(message));reference_fingerprint=fingerprint
  call require(state%valid.and.state%global_count==n.and.state%noccupied==noccupied,&
    'validated hybrid state metadata is invalid')
  call require(size(state%coefficients,1)==nowned.and.size(state%coefficients,2)==noccupied,&
    'validated hybrid coefficient shape is invalid')
  call require(abs(sum(state%occupations)-4d0)<1d-14.and.workspace>0_int64,&
    'validated hybrid state receipts are invalid')

  if(nproc>1)then
    call validate_dg_hybrid_ground_state(comm,n,noccupied,row_ids,coefficients,occupations,eigenvalues,4d0,&
      merge(1102_int64,1101_int64,rank==0),2202_int64,3303_int64,4404_int64,1d-12,&
      state,workspace,fingerprint,ok,message)
    call require(.not.ok,'rank-disagreeing hybrid provenance was accepted')
  endif
  occupations(1)=-1d0
  call validate_dg_hybrid_ground_state(comm,n,noccupied,row_ids,coefficients,occupations,eigenvalues,4d0,&
    1101_int64,2202_int64,3303_int64,4404_int64,1d-12,state,workspace,fingerprint,ok,message)
  call require(.not.ok,'negative occupation was accepted');occupations=[2d0,2d0]
  call validate_dg_hybrid_ground_state(comm,n,noccupied,row_ids,coefficients,occupations,eigenvalues,3d0,&
    1101_int64,2202_int64,3303_int64,4404_int64,1d-12,state,workspace,fingerprint,ok,message)
  call require(.not.ok,'wrong hybrid electron count was accepted')

  allocate(duplicate_ids(nowned+merge(1,0,rank==0.and.nowned>0)))
  duplicate_ids(1:nowned)=row_ids
  if(size(duplicate_ids)>nowned)duplicate_ids(size(duplicate_ids))=row_ids(1)
  deallocate(coefficients);allocate(coefficients(size(duplicate_ids),noccupied));coefficients=(0d0,0d0)
  call validate_dg_hybrid_ground_state(comm,n,noccupied,duplicate_ids,coefficients,occupations,eigenvalues,4d0,&
    1101_int64,2202_int64,3303_int64,4404_int64,1d-12,state,workspace,fingerprint,ok,message)
  call require(.not.ok,'duplicate hybrid row ownership was accepted')

  if(rank==0)then
    write(*,'(a,i0,a,i0)')'HYBRID_GROUND_STATE_TYPES ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid ground-state types on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)error stop label
  end subroutine require
end program test_dg_hybrid_ground_state_types_mpi
