#include "config.h"
program test_dg_overlapping_wannier_nonlocal_mpi
  use mpi
  use dg_overlapping_wannier_nonlocal,only:assemble_dg_overlapping_wannier_nonlocal,&
    assemble_dg_overlapping_wannier_nonlocal_rows
  implicit none
  integer::comm,rank,nproc,ierr,p,i,j,nlocal,index,owned
  integer(8),allocatable::ids(:),row_ids(:)
  real(8),allocatable::strength(:)
  complex(8),allocatable::overlap(:,:),matrix(:,:),reference(:,:),rotated(:,:),matrix_rows(:,:)
  complex(8)::direct_local(3,3),direct_global(3,3)
  logical,allocatable::complete(:,:)
  complex(8)::gauge(3,3)
  logical::ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=count([(mod(p-1,nproc)==rank,p=1,3)])
  allocate(ids(nlocal),strength(nlocal),overlap(3,nlocal),complete(3,nlocal))
  index=0
  do p=1,3
    if(mod(p-1,nproc)/=rank)cycle
    index=index+1;ids(index)=p;strength(index)=0.2d0+0.1d0*p
    overlap(:,index)=[cmplx(0.4d0*p,0.03d0*p,8),cmplx(-0.2d0*p,0.05d0,8),&
      cmplx(0.1d0,0.07d0*p,8)]
  enddo
  complete=.true.
  direct_local=(0d0,0d0)
  do p=1,nlocal;do j=1,3;do i=1,3
    direct_local(i,j)=direct_local(i,j)+strength(p)*conjg(overlap(i,p))*overlap(j,p)
  enddo;enddo;enddo
  call MPI_Allreduce(direct_local,direct_global,9,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call assemble_dg_overlapping_wannier_nonlocal(comm,3,ids,strength,overlap,complete,3_8,&
    matrix,owned,ok,message)
  call require(ok,trim(message));call require(owned==3,'unique projector ownership')
  call require(maxval(abs(matrix-conjg(transpose(matrix))))<1d-13,'Hermitian nonlocal matrix')
  call require(maxval(abs(matrix-direct_global))<1d-13,'direct nonlocal reference')
  call require(abs(matrix(1,2))>1d-8,'overlapping tail-projector block')
  reference=matrix
  allocate(row_ids(count([(mod(i-1,nproc)==rank,i=1,3)])))
  index=0
  do i=1,3
    if(mod(i-1,nproc)/=rank)cycle
    index=index+1;row_ids(index)=i
  enddo
  call assemble_dg_overlapping_wannier_nonlocal_rows(comm,3,row_ids,ids,strength,overlap,complete,&
    3_8,matrix_rows,owned,ok,message)
  call require(ok,trim(message))
  if(size(row_ids)>0)then
    direct_local=(0d0,0d0)
    direct_local(1:size(row_ids),:)=matrix_rows-reference(int(row_ids),:)
  else
    direct_local=(0d0,0d0)
  endif
  call require(maxval(abs(direct_local))<1d-13,'row-owned nonlocal reference')
  gauge=(0d0,0d0);gauge(1,1)=sqrt(0.5d0);gauge(1,2)=sqrt(0.5d0)
  gauge(2,1)=-sqrt(0.5d0);gauge(2,2)=sqrt(0.5d0);gauge(3,3)=1d0
  rotated=matmul(gauge,overlap)
  call assemble_dg_overlapping_wannier_nonlocal(comm,3,ids,strength,rotated,complete,3_8,&
    matrix,owned,ok,message)
  call require(ok,trim(message))
  call require(maxval(abs(matrix-matmul(conjg(gauge),matmul(reference,transpose(gauge)))))<1d-12,&
    'nonlocal gauge covariance')
  if(nproc>1)then
    call assemble_dg_overlapping_wannier_nonlocal_rows(comm,merge(4,3,rank==0),row_ids,ids,&
      strength,overlap,complete,3_8,matrix_rows,owned,ok,message)
    call require(.not.ok,'rank-inconsistent row-owned nonlocal nwann rejection')
    call assemble_dg_overlapping_wannier_nonlocal_rows(comm,3,row_ids,ids,strength,overlap,&
      complete,merge(4_8,3_8,rank==0),matrix_rows,owned,ok,message)
    call require(.not.ok,'rank-inconsistent row-owned projector-count rejection')
    call assemble_dg_overlapping_wannier_nonlocal(comm,merge(4,3,rank==0),ids,strength,overlap,&
      complete,3_8,matrix,owned,ok,message)
    call require(.not.ok,'rank-inconsistent nonlocal contract rejection')
  endif
  if(rank==0.and.size(ids)>0)complete(1,1)=.false.
  call assemble_dg_overlapping_wannier_nonlocal(comm,3,ids,strength,overlap,complete,3_8,&
    matrix,owned,ok,message)
  call require(.not.ok,'incomplete tail-projector rejection');complete=.true.
  if(rank==0.and.size(ids)>0)then
    ids=[ids,ids(1)];strength=[strength,strength(1)]
    overlap=reshape([overlap,overlap(:,1)],[3,size(ids)])
    complete=reshape([complete,complete(:,1)],[3,size(ids)])
  endif
  call assemble_dg_overlapping_wannier_nonlocal(comm,3,ids,strength,overlap,complete,3_8,&
    matrix,owned,ok,message)
  call require(.not.ok,'duplicate projector rejection')
  if(rank==0)then
    write(*,'(a,i0,a,4(es24.16,1x))')'NONLOCAL ranks=',nproc,' values=',real(reference(1,1)),&
      aimag(reference(1,2)),real(reference(2,3)),aimag(reference(2,3))
    write(*,'(a,i0,a)')'PASS overlapping-Wannier nonlocal on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    integer::lf,gf
    lf=merge(0,1,condition);call MPI_Allreduce(lf,gf,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(gf/=0)error stop label
  end subroutine
end program
