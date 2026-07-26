#include "config.h"
program test_dg_overlapping_wannier_operators_mpi
  use mpi
  use dg_overlapping_wannier_operators,only:assemble_dg_overlapping_wannier_weak_operators
  implicit none
  integer::comm,rank,nproc,ierr,p,i,j,nlocal,index,owned
  integer(8),allocatable::ids(:)
  real(8),allocatable::weights(:),vlocal(:)
  complex(8),allocatable::values(:,:),gradients(:,:,:),kinetic(:,:),potential(:,:),&
    reference_kinetic(:,:),reference_potential(:,:),rotated_values(:,:),rotated_gradients(:,:,:)
  complex(8)::gauge(3,3)
  logical::ok
  character(256)::message
  real(8)::x,k
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=count([(mod(p-1,nproc)==rank,p=1,6)])
  allocate(ids(nlocal),weights(nlocal),vlocal(nlocal),values(3,nlocal),gradients(3,3,nlocal))
  index=0;k=0.7d0
  do p=1,6
    if(mod(p-1,nproc)/=rank)cycle
    index=index+1;x=0.2d0*p;ids(index)=p;weights(index)=0.4d0+0.03d0*p
    vlocal(index)=-0.5d0+0.08d0*p
    values(1,index)=1d0
    values(2,index)=exp(cmplx(0d0,k*x,8))
    values(3,index)=0.4d0*values(1,index)+0.3d0*values(2,index)
    gradients(:,1,index)=(0d0,0d0)
    gradients(:,2,index)=[cmplx(0d0,k,8)*values(2,index),(0d0,0d0),(0d0,0d0)]
    gradients(:,3,index)=0.3d0*gradients(:,2,index)
  enddo
  call assemble_dg_overlapping_wannier_weak_operators(comm,3,ids,weights,values,gradients,&
    vlocal,6_8,kinetic,potential,owned,ok,message)
  call require(ok,trim(message));call require(owned==6,'unique core ownership')
  call require(maxval(abs(kinetic-conjg(transpose(kinetic))))<1d-13,'Hermitian kinetic')
  call require(maxval(abs(potential-conjg(transpose(potential))))<1d-13,'Hermitian potential')
  call require(abs(kinetic(1,1))<1d-14,'constant kinetic reference')
  call require(abs(real(kinetic(2,2))-0.5d0*k*k*sum_global(weights))<1d-13,&
    'plane-wave kinetic reference')
  call require(abs(kinetic(2,3))>1d-8.and.abs(potential(1,3))>1d-8,'off-fragment tail blocks')
  reference_kinetic=kinetic;reference_potential=potential

  gauge=(0d0,0d0);gauge(1,2)=1d0;gauge(2,1)=-1d0;gauge(3,3)=1d0
  rotated_values=matmul(gauge,values);allocate(rotated_gradients(3,3,nlocal))
  do p=1,nlocal;do i=1,3
    rotated_gradients(i,:,p)=matmul(gauge,gradients(i,:,p))
  enddo;enddo
  call assemble_dg_overlapping_wannier_weak_operators(comm,3,ids,weights,rotated_values,&
    rotated_gradients,vlocal,6_8,kinetic,potential,owned,ok,message)
  call require(ok,trim(message))
  call require(maxval(abs(kinetic-matmul(conjg(gauge),matmul(reference_kinetic,transpose(gauge)))))<1d-12,&
    'kinetic retained-space covariance')
  call require(maxval(abs(potential-matmul(conjg(gauge),matmul(reference_potential,transpose(gauge)))))<1d-12,&
    'potential retained-space covariance')

  gauge=(0d0,0d0);gauge(1,1)=sqrt(0.5d0);gauge(1,2)=sqrt(0.5d0)
  gauge(2,1)=-sqrt(0.5d0);gauge(2,2)=sqrt(0.5d0);gauge(3,3)=1d0
  rotated_values=matmul(gauge,values)
  do p=1,nlocal;do i=1,3
    rotated_gradients(i,:,p)=matmul(gauge,gradients(i,:,p))
  enddo;enddo
  call assemble_dg_overlapping_wannier_weak_operators(comm,3,ids,weights,rotated_values,&
    rotated_gradients,vlocal,6_8,kinetic,potential,owned,ok,message)
  call require(ok,trim(message))
  call require(maxval(abs(kinetic-matmul(conjg(gauge),matmul(reference_kinetic,transpose(gauge)))))<1d-12,&
    'kinetic retained-space rotation covariance')
  call require(maxval(abs(potential-matmul(conjg(gauge),matmul(reference_potential,transpose(gauge)))))<1d-12,&
    'potential retained-space rotation covariance')

  if(rank==0.and.size(ids)>0)then
    ids=[ids,ids(1)];weights=[weights,weights(1)];vlocal=[vlocal,vlocal(1)]
    values=reshape([values,values(:,1)],[3,size(ids)])
    gradients=reshape([gradients,gradients(:,:,1)],[3,3,size(ids)])
  endif
  call assemble_dg_overlapping_wannier_weak_operators(comm,3,ids,weights,values,gradients,&
    vlocal,6_8,kinetic,potential,owned,ok,message)
  call require(.not.ok,'duplicate core rejection')
  if(rank==0)then
    write(*,'(a,i0,a,4(es24.16,1x))')'OPERATORS ranks=',nproc,' values=',real(reference_kinetic(2,2)),&
      aimag(reference_kinetic(2,3)),real(reference_potential(1,3)),aimag(reference_potential(2,3))
    write(*,'(a,i0,a)')'PASS overlapping-Wannier weak operators on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  real(8) function sum_global(local)
    real(8),intent(in)::local(:)
    real(8)::partial
    partial=sum(local);call MPI_Allreduce(partial,sum_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  end function
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    integer::lf,gf
    lf=merge(0,1,condition);call MPI_Allreduce(lf,gf,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(gf/=0)error stop label
  end subroutine
end program
