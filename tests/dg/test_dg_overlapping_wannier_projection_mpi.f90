#include "config.h"
program test_dg_overlapping_wannier_projection_mpi
  use mpi
  use dg_overlapping_wannier_projection,only:t_dg_projection_channel,&
    build_dg_complete_sp_manifest,dg_complete_sp_target_count,validate_dg_complete_sp_manifest
  implicit none
  type(t_dg_projection_channel),allocatable::channels(:)
  integer::ierr,rank,target
  logical::ok
  character(256)::message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

  call build_dg_complete_sp_manifest([7,11],channels,ok,message)
  call require(ok,trim(message))
  call require(size(channels)==8,'two complete s+p shells')
  call require(all(channels%atom==[7,7,7,7,11,11,11,11]),'atom-major ownership')
  call require(all(channels%l==[0,1,1,1,0,1,1,1]),'complete s+p angular ordering')
  call require(all(channels%m==[1,1,2,3,1,1,2,3]),'complete real-harmonic ordering')
  call validate_dg_complete_sp_manifest(channels,[7,11],ok,message)
  call require(ok,trim(message))
  channels(4)%m=2
  call validate_dg_complete_sp_manifest(channels,[7,11],ok,message)
  call require(.not.ok,'missing p member rejected')
  channels(4)%m=3
  call validate_dg_complete_sp_manifest(channels(:7),[7,11],ok,message)
  call require(.not.ok,'truncated shell rejected')
  call build_dg_complete_sp_manifest([7,7],channels,ok,message)
  call require(.not.ok,'duplicate core-owned atom rejected')

  call dg_complete_sp_target_count(2,target,ok,message)
  call require(ok.and.target==8,'complete s+p target count')
  call dg_complete_sp_target_count(0,target,ok,message)
  call require(.not.ok,'empty core atom set rejected')

  if(rank==0)write(*,'(a)')'PASS complete-s+p projection manifest'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,text)
    logical,intent(in)::condition
    character(*),intent(in)::text
    if(.not.condition)then
      write(0,'(a,i0,2a)')'rank ',rank,': ',trim(text)
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
  end subroutine
end program
