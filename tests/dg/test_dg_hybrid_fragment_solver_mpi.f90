#include "config.h"
program test_dg_hybrid_fragment_solver_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis,build_dg_hybrid_fragment_basis
  use dg_hybrid_fragment_solver,only:solve_dg_hybrid_fragment_basis
  implicit none
  integer,parameter::npoint=4,nbasis=4,nstate=2
  integer::comm,rank,nproc,ierr,nowned,j,k
  integer(int64),allocatable::wf_ids(:),pw_ids(:)
  complex(real64),allocatable::wf_values(:,:),pw_values(:,:),coefficients(:,:)
  complex(real64)::vectors(npoint,nbasis),hop(npoint,npoint),sop(npoint,npoint)
  real(real64)::occupations(nstate),eigenvalues(nstate),density(npoint),residual,orthogonality,electron_count
  logical::core_mask(npoint),ok
  type(s_dg_hybrid_fragment_basis)::basis
  integer(int64)::workspace,fingerprint
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  vectors=(0d0,0d0);do j=1,nbasis;vectors(j,j)=1d0;enddo
  hop=(0d0,0d0);sop=(0d0,0d0)
  do j=1,npoint;hop(j,j)=real(j,real64);sop(j,j)=1d0;enddo
  nowned=count([(mod(j-1,nproc)==rank,j=1,nbasis)])
  allocate(wf_ids(count([(mod(j-1,nproc)==rank.and.j<=2,j=1,nbasis)])),&
    pw_ids(count([(mod(j-1,nproc)==rank.and.j>2,j=1,nbasis)])))
  allocate(wf_values(npoint,size(wf_ids)),pw_values(npoint,size(pw_ids)));j=0;k=0
  if(size(wf_ids)>0)wf_values=(0d0,0d0);if(size(pw_ids)>0)pw_values=(0d0,0d0)
  do nowned=1,nbasis
    if(mod(nowned-1,nproc)/=rank)cycle
    if(nowned<=2)then;j=j+1;wf_ids(j)=100+nstated(nowned);wf_values(:,j)=vectors(:,nowned)
    else;k=k+1;pw_ids(k)=100+nstated(nowned);pw_values(:,k)=vectors(:,nowned);endif
  enddo
  call build_dg_hybrid_fragment_basis(comm,1,wf_ids,wf_values,pw_ids,pw_values,0,0,basis,ok,message)
  call require(ok,trim(message));occupations=[2d0,0d0];core_mask=[.true.,.true.,.false.,.false.]
  call solve_dg_hybrid_fragment_basis(comm,basis,nstate,occupations,core_mask,apply_h,apply_s,1d-12,&
    coefficients,eigenvalues,density,electron_count,residual,orthogonality,workspace,fingerprint,ok,message)
  call require(ok,trim(message));call require(maxval(abs(eigenvalues-[1d0,2d0]))<1d-12,'fragment eigenvalues mismatch')
  call require(residual<1d-12.and.orthogonality<1d-12,'fragment eigensystem receipts mismatch')
  call require(abs(electron_count-2d0)<1d-12.and.abs(density(1)-2d0)<1d-12.and.&
    maxval(abs(density(2:)))<1d-12,'fragment core density mismatch')
  call require(size(coefficients,1)==size(basis%global_ids).and.size(coefficients,2)==nstate,&
    'fragment coefficients are not basis-row distributed')
  if(rank==0)write(*,'(a,i0,a,i0)')'HYBRID_FRAGMENT_SOLVER ranks=',nproc,' fingerprint=',fingerprint
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid fragment solver on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  integer(int64) function nstated(value)
    integer,intent(in)::value;nstated=int(value,int64)
  end function nstated
  subroutine apply_h(input,output,callback_ok)
    complex(real64),intent(in)::input(:,:);complex(real64),intent(out)::output(:,:);logical,intent(out)::callback_ok
    output=matmul(hop,input);callback_ok=.true.
  end subroutine apply_h
  subroutine apply_s(input,output,callback_ok)
    complex(real64),intent(in)::input(:,:);complex(real64),intent(out)::output(:,:);logical,intent(out)::callback_ok
    output=matmul(sop,input);callback_ok=.true.
  end subroutine apply_s
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text;logical::global_condition
    call MPI_Allreduce(condition,global_condition,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(.not.global_condition)then;if(rank==0)write(0,'(a)')trim(text);call MPI_Abort(comm,1,ierr);endif
  end subroutine require
end program test_dg_hybrid_fragment_solver_mpi
