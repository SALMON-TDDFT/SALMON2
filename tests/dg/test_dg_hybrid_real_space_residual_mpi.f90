#include "config.h"
program test_dg_hybrid_real_space_residual_mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use mpi
  use dg_hybrid_real_space_residual,only:evaluate_dg_hybrid_real_space_residual
  implicit none
  integer,parameter::nbasis=3,nstate=2,npoint_global=8
  integer::comm,rank,nproc,ierr,nlocal,p,i,row
  integer(int64),allocatable::point_ids(:),row_ids(:)
  real(real64),allocatable::weights(:)
  complex(real64),allocatable::basis(:,:),action(:,:),coefficients(:,:)
  real(real64)::eigenvalues(nstate),residual
  logical::ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=count([(mod(i-1,nproc)==rank,i=1,npoint_global)])
  allocate(point_ids(nlocal),weights(nlocal),basis(nbasis,nlocal),action(nbasis,nlocal))
  p=0
  do i=1,npoint_global
    if(mod(i-1,nproc)/=rank)cycle
    p=p+1;point_ids(p)=i;weights(p)=0.25d0
    basis(:,p)=[cmplx(1d0+0.1d0*i,0.02d0*i,real64),&
      cmplx(-0.3d0+0.04d0*i,0.07d0,real64),cmplx(0.2d0,-0.03d0*i,real64)]
  enddo
  nlocal=count([(mod(i-1,nproc)==rank,i=1,nbasis)])
  allocate(row_ids(nlocal),coefficients(nlocal,nstate));p=0
  do row=1,nbasis
    if(mod(row-1,nproc)/=rank)cycle
    p=p+1;row_ids(p)=row;coefficients(p,:)=(0d0,0d0)
    if(row<=nstate)coefficients(p,row)=1d0
  enddo
  eigenvalues=[-0.7d0,0.4d0]
  do p=1,size(point_ids)
    action(1,p)=eigenvalues(1)*basis(1,p)
    action(2,p)=eigenvalues(2)*basis(2,p)
    action(3,p)=cmplx(0.11d0,-0.08d0,real64)
  enddo
  call evaluate_dg_hybrid_real_space_residual(comm,npoint_global,point_ids,weights,nbasis,row_ids,&
    basis,action,coefficients,eigenvalues,residual,ok,message)
  call require(ok.and.residual<1d-13,'exact reconstructed DG action did not give zero residual')
  action(2,:)=action(2,:)+cmplx(0.02d0,-0.01d0,real64)
  call evaluate_dg_hybrid_real_space_residual(comm,npoint_global,point_ids,weights,nbasis,row_ids,&
    basis,action,coefficients,eigenvalues,residual,ok,message)
  call require(ok.and.residual>1d-4,'off-subspace real-space defect was not detected')
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid real-space residual on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;if(rank==0)write(0,'(a)')trim(label);error stop 1;endif
  end subroutine require
end program test_dg_hybrid_real_space_residual_mpi
