#include "config.h"
program test_dg_hybrid_wannier_selection_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_hybrid_wannier_selection,only:select_dg_hybrid_wannier_blocks
  implicit none
  integer::comm,rank,nproc,ierr,i
  integer::block_ids(8),permuted_ids(8),conjugates(4),permutation(8)
  integer,allocatable::accepted(:),rejected(:),complement_rank(:)
  integer(int64)::fingerprint,reference_fingerprint
  real(real64)::localization(8),permuted_localization(8),threshold
  logical::ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  block_ids=[1,1,2,2,3,3,4,4]
  conjugates=[2,1,3,4]
  localization=[0.10d0,0.15d0,0.20d0,0.70d0,0.12d0,0.18d0,0.08d0,0.09d0]
  threshold=0.50d0
  call select_dg_hybrid_wannier_blocks(comm,4,block_ids,conjugates,localization,threshold,&
    accepted,rejected,complement_rank,fingerprint,ok,message)
  call require(ok,trim(message))
  call require(all(accepted==[3,4]),'complete accepted block set is incorrect')
  call require(all(rejected==[1,2]),'conjugate rejection is not atomic')
  call require(all(complement_rank==[2,2,0,0]),'complement rank does not match rejected block size')
  call require(fingerprint/=0_int64,'selection fingerprint is zero')
  reference_fingerprint=fingerprint

  permutation=[8,3,6,1,7,4,2,5]
  do i=1,8
    permuted_ids(i)=block_ids(permutation(i))
    permuted_localization(i)=localization(permutation(i))
  enddo
  call select_dg_hybrid_wannier_blocks(comm,4,permuted_ids,conjugates,permuted_localization,threshold,&
    accepted,rejected,complement_rank,fingerprint,ok,message)
  call require(ok,trim(message))
  call require(fingerprint==reference_fingerprint,'selection depends on incoming Wannier order')

  if(nproc>1)then
    if(rank==0)threshold=0.49d0
  else
    threshold=-1d0
  endif
  call select_dg_hybrid_wannier_blocks(comm,4,block_ids,conjugates,localization,threshold,&
    accepted,rejected,complement_rank,fingerprint,ok,message)
  call require(.not.ok,'rank-disagreeing selection threshold was accepted')
  threshold=0.50d0

  localization(1)=ieee_value(0d0,ieee_quiet_nan)
  call select_dg_hybrid_wannier_blocks(comm,4,block_ids,conjugates,localization,threshold,&
    accepted,rejected,complement_rank,fingerprint,ok,message)
  call require(.not.ok,'nonfinite localization receipt was accepted')
  localization(1)=0.10d0

  if(rank==0)conjugates(1)=1
  call select_dg_hybrid_wannier_blocks(comm,4,block_ids,conjugates,localization,threshold,&
    accepted,rejected,complement_rank,fingerprint,ok,message)
  call require(.not.ok,'rank-disagreeing conjugate metadata was accepted')
  conjugates=[2,1,3,4]

  if(rank==0)then
    write(*,'(a,i0,a,i0)')'HYBRID_SELECTION ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid Wannier selection on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_failure,global_failure
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_failure/=0)error stop label
  end subroutine require
end program test_dg_hybrid_wannier_selection_mpi
