#include "config.h"
program test_dg_overlapping_wannier_construction_mpi
  use mpi
  use dg_overlapping_wannier_construction,only:s_dg_overlapping_wannier_construction,&
    construct_dg_overlapping_wannier_basis,release_dg_overlapping_wannier_construction,&
    verify_dg_overlapping_wannier_periodic_closure
  implicit none
  integer::comm,rank,nproc,ierr,i,p,nlocal,nclosure,index,ncore
  integer(8),allocatable::ids(:),box_ids(:),symmetry_map(:,:),broken_symmetry_map(:,:)
  integer,allocatable::fragment(:)
  real(8),allocatable::weight(:),coordinate(:)
  logical,allocatable::boundary(:),core_mask(:)
  complex(8),allocatable::candidate(:,:),gradient(:,:,:),occupied(:,:),base_occupied(:,:),&
    rotated(:,:),rotated_gradient(:,:,:),periodic_phase(:,:)
  complex(8)::gauge(4,4)
  complex(8),allocatable::reference_projector(:,:),projector(:,:)
  type(s_dg_overlapping_wannier_construction)::result
  integer,allocatable::reference_owner(:),reference_center_fragment(:)
  integer(8),allocatable::reference_center_box_ids(:)
  integer(8)::reference_fingerprint
  integer(8)::closure_fingerprint
  integer(8),allocatable::closure_ids(:),closure_map(:,:)
  complex(8),allocatable::closure_values(:,:),closure_gradients(:,:,:)
  complex(8)::closure_representation(2,2,1)
  real(8)::closure_rotation(3,3,1),closure_residual
  logical::ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=count([(mod(p-1,nproc)==rank,p=1,12)])
  allocate(ids(nlocal),box_ids(nlocal),symmetry_map(nlocal,2),broken_symmetry_map(nlocal,2),&
    fragment(nlocal),weight(nlocal),coordinate(nlocal),boundary(nlocal))
  allocate(core_mask(nlocal))
  allocate(candidate(4,nlocal),gradient(3,4,nlocal),occupied(4,2),periodic_phase(3,nlocal))
  index=0
  do p=1,12
    if(mod(p-1,nproc)/=rank)cycle
    index=index+1
    core_mask(index)=p>=3.and.p<=10
    if(core_mask(index))then
      ids(index)=int(p-2,8)
    else
      select case(p)
      case(1);ids(index)=7_8
      case(2);ids(index)=8_8
      case(11);ids(index)=1_8
      case default;ids(index)=2_8
      end select
    endif
    fragment(index)=rank+1;box_ids(index)=p
    symmetry_map(index,:)=[int(p,8),int(13-p,8)]
    broken_symmetry_map(index,:)=[int(p,8),int(mod(p,12)+1,8)]
    weight(index)=1d0
    coordinate(index)=dble(p)-6.5d0+0.01d0*(dble(p)-6.5d0)**2
    boundary(index)=p==1.or.p==12
    do i=1,3
      periodic_phase(i,index)=exp(cmplx(0d0,2d0*acos(-1d0)*i*(dble(p)-0.5d0)/12d0,8))
    enddo
    candidate(:,index)=[cmplx(1d0,0d0,8),&
      cmplx(cos(2d0*acos(-1d0)*(dble(p)-0.5d0)/12d0),0d0,8),&
      cmplx(sin(2d0*acos(-1d0)*(dble(p)-0.5d0)/12d0),0d0,8),&
      cmplx(cos(4d0*acos(-1d0)*(dble(p)-0.5d0)/12d0),0d0,8)]
    candidate(:,index)=candidate(:,index)*[cmplx(1d0,0d0,8),cmplx(0d0,1d0,8),&
      cmplx(sqrt(0.5d0),sqrt(0.5d0),8),cmplx(sqrt(0.5d0),-sqrt(0.5d0),8)]
    if(boundary(index))candidate(:,index)=1d-9*candidate(:,index)
    do i=1,4
      gradient(:,i,index)=[0.03d0*i*candidate(i,index),-0.02d0*p*candidate(i,index),&
        0.01d0*(i+p)*candidate(i,index)]
    enddo
  enddo
  occupied=(0d0,0d0);occupied(2,1)=1d0;occupied(3,2)=1d0
  base_occupied=occupied

  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    candidate,gradient,occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message,core_mask,&
    box_ids,symmetry_map,12_8,1d-9,periodic_phase)
  call require(ok,trim(message))
  call require(result%candidate_rank==4.and.result%target_rank==3.and.result%retained_rank==3,&
    'rank policy')
  call require(result%occupied_inclusion_residual<1d-9,'occupied inclusion')
  call require(result%symmetry_closure_residual<1d-9,'retained symmetry closure')
  call require(allocated(result%symmetry_representation),'retained symmetry representation')
  call require(maxval(abs(result%symmetry_representation(:,:,2)-identity3()))>1d-3,&
    'nontrivial periodic-box symmetry representation')
  call require(maxval(abs(matmul(conjg(transpose(result%symmetry_representation(:,:,2))),&
    result%symmetry_representation(:,:,2))-identity3()))<1d-9,'retained symmetry unitarity')
  call require(result%boundary_value_max<1d-8.and.result%boundary_gradient_max<1d-7,&
    'buffer boundary diagnostics')
  call require(allocated(result%center_box_point_ids),'periodic-box center ids')
  call require(all(result%center_box_point_ids>=3_8.and.result%center_box_point_ids<=10_8),&
    'all retained Wannier centers are core owned')
  call require(all(result%physical_grid_ids==ids),'complete physical tail ids')
  ncore=size(result%physical_grid_ids)-count(core_mask)
  call MPI_Allreduce(ncore,index,1,MPI_INTEGER,MPI_SUM,comm,ierr)
  call require(index==4,'buffer-only tail retention')
  call require(size(result%gradient,1)==3.and.size(result%gradient,3)==nlocal,'complete gradient tails')
  call local_projector(result%value,reference_projector)
  reference_owner=result%center_owner_rank;reference_center_fragment=result%center_owner_fragment
  reference_center_box_ids=result%center_box_point_ids
  reference_fingerprint=result%transform_fingerprint
  call release_dg_overlapping_wannier_construction(result)

  closure_representation=(0d0,0d0);closure_representation(1,1,1)=1d0;closure_representation(2,2,1)=1d0
  closure_rotation=0d0
  do i=1,3;closure_rotation(i,i,1)=1d0;enddo
  nclosure=count([(mod(p-1,nproc)==rank,p=1,2)])
  allocate(closure_ids(nclosure),closure_map(nclosure,1),closure_values(2,nclosure),closure_gradients(3,2,nclosure))
  index=0
  do p=1,2
    if(mod(p-1,nproc)/=rank)cycle
    index=index+1;closure_ids(index)=p;closure_map(index,1)=3-p
    closure_values(:,index)=[(1d0,0d0),(2d0,0d0)];closure_gradients(:,:,index)=1d0
  enddo
  call verify_dg_overlapping_wannier_periodic_closure(comm,closure_ids,closure_map,closure_values,&
    closure_gradients,closure_representation,closure_rotation,2_8,1d-12,closure_residual,&
    closure_fingerprint,ok,message)
  call require(ok.and.closure_residual<1d-12,'authoritative periodic value/gradient closure')
  if(rank==0)closure_gradients(1,1,1)=2d0
  call verify_dg_overlapping_wannier_periodic_closure(comm,closure_ids,closure_map,closure_values,&
    closure_gradients,closure_representation,closure_rotation,2_8,1d-12,closure_residual,&
    closure_fingerprint,ok,message)
  call require(.not.ok,'one-sided periodic gradient-tail corruption rejected collectively')

  coordinate=-1000d0*coordinate+37d0
  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    candidate,gradient,base_occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message,core_mask,&
    box_ids,symmetry_map,12_8,1d-9,periodic_phase)
  call require(ok,trim(message));call local_projector(result%value,projector)
  call require(maxval(abs(projector-reference_projector))<1d-9,&
    'periodic localization is independent of unbounded scalar coordinates')
  call release_dg_overlapping_wannier_construction(result)

  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    candidate,gradient,base_occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message,core_mask,&
    box_ids,symmetry_map,12_8,1d-9)
  call require(.not.ok,'periodic symmetry requires periodic localization phase links')

  periodic_phase(1,:)=2d0*periodic_phase(1,:)
  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    candidate,gradient,base_occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message,core_mask,&
    box_ids,symmetry_map,12_8,1d-9,periodic_phase)
  call require(.not.ok,'non-unit periodic localization phase rejection')
  periodic_phase(1,:)=0.5d0*periodic_phase(1,:)

  if(nproc>1)then
    call construct_dg_overlapping_wannier_basis(comm,merge(5,4,rank==0),3,2,ids,fragment,weight,&
      coordinate,boundary,candidate,gradient,base_occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,&
      message,core_mask,box_ids,symmetry_map,12_8,1d-9,periodic_phase)
    call require(.not.ok,'rank-inconsistent construction contract rejection')
    call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
      candidate,gradient,base_occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message,core_mask,&
      box_ids,symmetry_map(:,1:merge(1,2,rank==0)),12_8,1d-9,periodic_phase)
    call require(.not.ok,'rank-inconsistent symmetry-count rejection')
  endif

  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    candidate,gradient,base_occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message,core_mask,&
    box_ids,broken_symmetry_map,12_8,1d-9,periodic_phase)
  call require(.not.ok,'non-group-closed periodic-box symmetry rejection')

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
    rotated,rotated_gradient,occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message,core_mask,&
    box_ids,symmetry_map,12_8,1d-9,periodic_phase)
  call require(ok,trim(message));call local_projector(result%value,projector)
  call require(maxval(abs(projector-reference_projector))<1d-9,'candidate-window gauge invariance')
  call require(all(result%center_owner_rank==reference_owner),'deterministic center ownership')
  call require(all(result%center_owner_fragment>=1.and.result%center_owner_fragment<=nproc),&
    'centers retain their real DC fragment ownership')
  call require(all(result%center_box_point_ids==reference_center_box_ids),&
    'deterministic symmetry-compatible centers under candidate gauge rotation')
  call require(result%transform_fingerprint==reference_fingerprint,'deterministic transform fingerprint')
  call release_dg_overlapping_wannier_construction(result)

  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    candidate,gradient,base_occupied,8_8,41,1d-12,1d-12,1d-9,result,ok,message,core_mask,&
    box_ids,symmetry_map,12_8,1d-9,periodic_phase)
  call require(.not.ok,'buffer-boundary gate')

  rotated=candidate;rotated(4,:)=rotated(3,:)
  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    rotated,gradient,base_occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message,core_mask,&
    box_ids,symmetry_map,12_8,1d-9,periodic_phase)
  call require(.not.ok,'candidate rank-loss gate')

  occupied=base_occupied;occupied(:,2)=occupied(:,1)
  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    candidate,gradient,occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message,core_mask,&
    box_ids,symmetry_map,12_8,1d-9,periodic_phase)
  call require(.not.ok,'occupied rank-loss gate')

  if(rank==0)then
    write(*,'(a,i0,a,i0,a,*(i0,1x))')'CONSTRUCTION ranks=',nproc,' fingerprint=',&
      reference_fingerprint,' centers=',reference_center_box_ids
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
  function identity3() result(identity)
    complex(8)::identity(3,3)
    integer::j
    identity=(0d0,0d0)
    do j=1,3;identity(j,j)=1d0;enddo
  end function
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_failure,global_failure
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_failure/=0)error stop label
  end subroutine
end program
