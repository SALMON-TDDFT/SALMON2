#include "config.h"
program test_dg_hybrid_low_energy_symmetry_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:real64
  use dg_hybrid_low_energy_symmetry,only:select_dg_hybrid_symmetry_target,&
    evaluate_dg_hybrid_low_energy_symmetry
  implicit none
  integer::comm,rank,nproc,ierr,i,target_rank
  real(real64)::eigenvalues(3),occupied_defect,target_defect,energy_defect
  real(real64)::clustered(4),unresolved(3)
  complex(real64)::metric(3,3),representation(3,3,2),coefficients(3,3)
  logical::ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)

  metric=0d0;coefficients=(0d0,0d0);representation=(0d0,0d0)
  do i=1,3
    metric(i,i)=1d0;coefficients(i,i)=1d0;representation(i,i,1)=1d0
  enddo
  representation(1,1,2)=1d0;representation(2,2,2)=1d0;representation(3,3,2)=0.5d0
  eigenvalues=[-1d0,0d0,1d0]
  call evaluate_dg_hybrid_low_energy_symmetry(comm,metric,representation,coefficients,eigenvalues,&
    1,2,1d-12,occupied_defect,target_defect,energy_defect,ok,message)
  call require(ok,'closed low-energy space was rejected: '//trim(message))
  call require(max(occupied_defect,target_defect,energy_defect)<1d-13,&
    'unused high-energy nonunitarity contaminated target defects')

  representation(2,2,2)=0d0;representation(3,2,2)=1d0
  call evaluate_dg_hybrid_low_energy_symmetry(comm,metric,representation,coefficients,eigenvalues,&
    1,2,1d-12,occupied_defect,target_defect,energy_defect,ok,message)
  call require(.not.ok.and.target_defect>0.5d0,'target leakage was accepted')

  representation=(0d0,0d0)
  do i=1,3
    representation(i,i,1)=1d0
  enddo
  representation(1,2,2)=1d0;representation(2,1,2)=1d0;representation(3,3,2)=1d0
  call evaluate_dg_hybrid_low_energy_symmetry(comm,metric,representation,coefficients,eigenvalues,&
    1,2,1d-12,occupied_defect,target_defect,energy_defect,ok,message)
  call require(.not.ok.and.target_defect<1d-13.and.energy_defect>0.5d0,&
    'nondegenerate target energy mixing was accepted')

  clustered=[-1d0,0d0,0d0,1d0]
  call select_dg_hybrid_symmetry_target(comm,clustered,2,1d-12,target_rank,ok,message)
  call require(ok.and.target_rank==3,'degenerate target boundary was not extended')
  unresolved=[-1d0,0d0,0d0]
  call select_dg_hybrid_symmetry_target(comm,unresolved,2,1d-12,target_rank,ok,message)
  call require(.not.ok,'unresolved terminal degeneracy was accepted')
  call select_dg_hybrid_symmetry_target(comm,clustered,2+merge(1,0,rank>0),1d-12,target_rank,ok,message)
  if(nproc>1)call require(.not.ok,'rank-disagreeing target request was accepted')

  if(rank==0)write(*,'(a,i0,a)')'PASS low-energy Hybrid symmetry on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(rank==0)write(0,'(a)')trim(label)
      call MPI_Abort(comm,1,ierr)
    endif
  end subroutine require
end program test_dg_hybrid_low_energy_symmetry_mpi
