#include "config.h"
program test_dg_overlapping_wannier_construction_mpi
  use mpi
  use dg_overlapping_wannier_construction,only:s_dg_overlapping_wannier_construction,&
    construct_dg_overlapping_wannier_basis,release_dg_overlapping_wannier_construction
  implicit none
  integer::comm,rank,nproc,ierr,i,p,nlocal,index
  integer(8),allocatable::ids(:)
  integer,allocatable::fragment(:)
  real(8),allocatable::weight(:),coordinate(:)
  logical,allocatable::boundary(:)
  complex(8),allocatable::candidate(:,:),gradient(:,:,:),occupied(:,:),base_occupied(:,:),&
    rotated(:,:),rotated_gradient(:,:,:)
  complex(8)::gauge(4,4)
  complex(8),allocatable::reference_projector(:,:),projector(:,:)
  type(s_dg_overlapping_wannier_construction)::result
  integer,allocatable::reference_owner(:),reference_center_fragment(:)
  integer(8)::reference_fingerprint
  logical::ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=count([(mod(p-1,nproc)==rank,p=1,8)])
  allocate(ids(nlocal),fragment(nlocal),weight(nlocal),coordinate(nlocal),boundary(nlocal))
  allocate(candidate(4,nlocal),gradient(3,4,nlocal),occupied(4,2))
  index=0
  do p=1,8
    if(mod(p-1,nproc)/=rank)cycle
    index=index+1;ids(index)=p;fragment(index)=merge(1,2,p<=4)
    weight(index)=0.5d0+0.05d0*p;coordinate(index)=dble(p)
    boundary(index)=p>=7
    candidate(:,index)=[cmplx(1d0+0.1d0*p,0.02d0*p,8),&
      cmplx(sin(0.7d0*p),0.03d0*p,8),cmplx(cos(0.4d0*p),-0.01d0*p,8),&
      cmplx((-1d0)**p*0.2d0*p,0.015d0*p,8)]
    if(boundary(index))candidate(:,index)=1d-9*candidate(:,index)
    do i=1,4
      gradient(:,i,index)=[0.03d0*i*candidate(i,index),-0.02d0*p*candidate(i,index),&
        0.01d0*(i+p)*candidate(i,index)]
    enddo
  enddo
  occupied=(0d0,0d0);occupied(1,1)=1d0;occupied(2,2)=1d0
  base_occupied=occupied

  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    candidate,gradient,occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message)
  call require(ok,trim(message))
  call require(result%candidate_rank==4.and.result%target_rank==3.and.result%retained_rank==3,&
    'rank policy')
  call require(result%occupied_inclusion_residual<1d-9,'occupied inclusion')
  call require(result%boundary_value_max<1d-8.and.result%boundary_gradient_max<1d-7,&
    'buffer boundary diagnostics')
  call require(all(result%physical_grid_ids==ids),'complete physical tail ids')
  call require(size(result%gradient,1)==3.and.size(result%gradient,3)==nlocal,'complete gradient tails')
  call local_projector(result%value,reference_projector)
  reference_owner=result%center_owner_rank;reference_center_fragment=result%center_owner_fragment
  reference_fingerprint=result%transform_fingerprint
  call release_dg_overlapping_wannier_construction(result)

  gauge=(0d0,0d0);gauge(1,2)=1d0;gauge(2,1)=-1d0
  gauge(3,3)=sqrt(0.5d0);gauge(3,4)=sqrt(0.5d0)
  gauge(4,3)=-sqrt(0.5d0);gauge(4,4)=sqrt(0.5d0)
  rotated=matmul(gauge,candidate)
  allocate(rotated_gradient(3,4,nlocal))
  do p=1,nlocal;do i=1,3
    rotated_gradient(i,:,p)=matmul(gauge,gradient(i,:,p))
  enddo;enddo
  occupied=matmul(conjg(gauge),occupied)
  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    rotated,rotated_gradient,occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message)
  call require(ok,trim(message));call local_projector(result%value,projector)
  call require(maxval(abs(projector-reference_projector))<1d-9,'candidate-window gauge invariance')
  call require(all(result%center_owner_rank==reference_owner),'deterministic center ownership')
  call require(result%transform_fingerprint==reference_fingerprint,'deterministic transform fingerprint')
  call release_dg_overlapping_wannier_construction(result)

  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    candidate,gradient,base_occupied,8_8,41,1d-12,1d-12,1d-9,result,ok,message)
  call require(.not.ok,'buffer-boundary gate')

  rotated=candidate;rotated(4,:)=rotated(3,:)
  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    rotated,gradient,base_occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message)
  call require(.not.ok,'candidate rank-loss gate')

  occupied=base_occupied;occupied(:,2)=occupied(:,1)
  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    candidate,gradient,occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message)
  call require(.not.ok,'occupied rank-loss gate')

  if(rank==0)then
    write(*,'(a,i0,a,i0,a,*(i0,1x))')'CONSTRUCTION ranks=',nproc,' fingerprint=',&
      reference_fingerprint,' centers=',reference_center_fragment
    write(*,'(a,i0,a)')'PASS overlapping-Wannier construction on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine local_projector(values,projection)
    complex(8),intent(in)::values(:,:)
    complex(8),allocatable,intent(out)::projection(:,:)
    allocate(projection(size(values,2),size(values,2)))
    projection=matmul(transpose(conjg(values)),values)
  end subroutine
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_failure,global_failure
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_failure/=0)error stop label
  end subroutine
end program
