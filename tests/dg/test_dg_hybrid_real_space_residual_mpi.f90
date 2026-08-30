#include "config.h"
program test_dg_hybrid_real_space_residual_mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use mpi
  use dg_hybrid_real_space_residual,only:evaluate_dg_hybrid_real_space_residual,&
    evaluate_dg_hybrid_face_action_residuals
  implicit none
  integer,parameter::nbasis=3,nstate=2,npoint_global=8
  integer::comm,rank,nproc,ierr,nlocal,p,i,row
  integer(int64),allocatable::point_ids(:),row_ids(:)
  real(real64),allocatable::weights(:)
  complex(real64),allocatable::basis(:,:),action(:,:),coefficients(:,:),face_rows(:,:,:),face_actions(:,:,:),&
    hamiltonian_rows(:,:),metric_eigen_action(:,:)
  real(real64)::eigenvalues(nstate),residual,face_residuals(3)
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
  action(2,:)=action(2,:)-cmplx(0.02d0,-0.01d0,real64)
  call evaluate_dg_hybrid_real_space_residual(comm,npoint_global,point_ids,weights,nbasis,row_ids,&
    basis,action,coefficients,eigenvalues,residual,ok,message)
  call require(ok.and.residual<1d-13,'boundary-only fixture did not have zero volume coefficient residual')
  allocate(face_rows(size(row_ids),nbasis,3),face_actions(size(row_ids),nstate,3),&
    hamiltonian_rows(size(row_ids),nbasis),metric_eigen_action(size(row_ids),nstate))
  face_rows=(0d0,0d0)
  do p=1,size(row_ids)
    face_rows(p,int(row_ids(p)),1)=cmplx(0.3d0,0.1d0,real64)
    face_rows(p,int(row_ids(p)),2)=cmplx(-0.2d0,0.15d0,real64)
    face_rows(p,int(row_ids(p)),3)=cmplx(0.4d0,-0.05d0,real64)
  enddo
  hamiltonian_rows=sum(face_rows,dim=3)
  do p=1,size(row_ids)
    if(row_ids(p)<=nstate)hamiltonian_rows(p,int(row_ids(p)))=&
      hamiltonian_rows(p,int(row_ids(p)))+eigenvalues(int(row_ids(p)))
  enddo
  metric_eigen_action=matmul(hamiltonian_rows,global_test_coefficients())
  face_actions=(0d0,0d0)
  call evaluate_dg_hybrid_face_action_residuals(comm,nbasis,row_ids,face_rows,hamiltonian_rows,coefficients,&
    metric_eigen_action,1d0,&
    face_actions,face_residuals,ok,message)
  call require(ok.and.all(face_residuals>1d-4),&
    'zero coefficient and trace residuals masked missing SIPG boundary action')
  do i=1,3
    face_actions(:,:,i)=matmul(face_rows(:,:,i),global_test_coefficients())
  enddo
  call evaluate_dg_hybrid_face_action_residuals(comm,nbasis,row_ids,face_rows,hamiltonian_rows,coefficients,&
    metric_eigen_action,1d0,&
    face_actions,face_residuals,ok,message)
  call require(ok.and.all(face_residuals<1d-13),'exact SIPG component actions did not give zero residuals')
  do i=1,3
    row=findloc(row_ids,1_int64,dim=1)
    if(row>0)face_actions(row,1,i)=face_actions(row,1,i)+cmplx(0.07d0,-0.03d0,real64)
    call evaluate_dg_hybrid_face_action_residuals(comm,nbasis,row_ids,face_rows,hamiltonian_rows,coefficients,&
      metric_eigen_action,1d0,&
      face_actions,face_residuals,ok,message)
    call require(ok.and.face_residuals(i)>1d-4.and.all(face_residuals(pack([(p,p=1,3)],[(p/=i,p=1,3)]))<1d-13),&
      'SIPG component residuals are not independent')
    if(row>0)face_actions(row,1,i)=face_actions(row,1,i)-cmplx(0.07d0,-0.03d0,real64)
  enddo
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid real-space residual on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  function global_test_coefficients() result(global_coefficients)
    complex(real64)::global_coefficients(nbasis,nstate)
    global_coefficients=(0d0,0d0)
    global_coefficients(1,1)=1d0;global_coefficients(2,2)=1d0
  end function global_test_coefficients
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;if(rank==0)write(0,'(a)')trim(label);error stop 1;endif
  end subroutine require
end program test_dg_hybrid_real_space_residual_mpi
