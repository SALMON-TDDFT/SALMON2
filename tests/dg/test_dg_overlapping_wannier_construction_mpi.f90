#include "config.h"
program test_dg_overlapping_wannier_construction_mpi
  use mpi
  use dg_overlapping_wannier_construction,only:s_dg_overlapping_wannier_construction,&
    construct_dg_overlapping_wannier_basis,release_dg_overlapping_wannier_construction,&
    verify_dg_overlapping_wannier_periodic_closure,assemble_dg_distributed_candidate_symmetry,&
    align_dg_fragment_wannier_gauge,replicate_dg_fragment_wannier_representative,&
    verify_dg_fragment_wannier_streaming_closure,verify_dg_fragment_center_orbit,&
    verify_dg_uniform_fragment_target_rank,assign_dg_overlapping_wannier_occupations
  implicit none
  integer::comm,rank,nproc,ierr,i,p,nlocal,nclosure,index,ncore,fragment_id
  integer(8),allocatable::ids(:),box_ids(:),symmetry_map(:,:),broken_symmetry_map(:,:)
  integer,allocatable::fragment(:)
  real(8),allocatable::weight(:),coordinate(:)
  logical,allocatable::boundary(:),core_mask(:)
  complex(8),allocatable::candidate(:,:),gradient(:,:,:),occupied(:,:),base_occupied(:,:),&
    rotated(:,:),rotated_gradient(:,:,:),periodic_phase(:,:)
  complex(8),allocatable::direct_small(:,:),direct_small_gradient(:,:,:),&
    direct_small_occupied(:,:),direct_projector(:,:)
  complex(8),allocatable::stream_values(:,:),stream_gradients(:,:,:)
  complex(8),allocatable::mismatch_values(:,:),mismatch_gradients(:,:,:)
  complex(8)::gauge(4,4)
  complex(8),allocatable::reference_projector(:,:),projector(:,:),seed_projector(:,:)
  complex(8),allocatable::distributed_candidate(:,:),distributed_overlap(:,:,:)
  integer(8),allocatable::distributed_map(:,:)
  real(8),allocatable::distributed_weight(:)
  real(8),allocatable::seed_values(:,:)
  real(8),allocatable::occupied_seed_values(:,:)
  real(8),allocatable::raw_seed_values(:,:)
  real(8)::occupations(3)
  type(s_dg_overlapping_wannier_construction)::result
  integer,allocatable::reference_owner(:),reference_center_fragment(:)
  integer(8),allocatable::reference_center_box_ids(:)
  integer(8)::reference_fingerprint
  integer(8)::closure_fingerprint,rounded_closure_fingerprint
  integer(8),allocatable::closure_ids(:),closure_map(:,:)
  integer(8),allocatable::stream_ids(:),stream_map(:,:)
  integer(8)::local_centers(2),global_centers(2,2),center_orbit_map(4,2)
  complex(8),allocatable::closure_values(:,:),closure_gradients(:,:,:)
  complex(8)::closure_representation(2,2,1)
  real(8)::closure_rotation(3,3,1),closure_residual
  real(8)::gauge_residual,gauge_correction,theta
  complex(8)::gauge_values(2,2),gauge_gradients(3,2,2)
  real(8)::gauge_weights(2)
  logical::ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call assign_dg_overlapping_wannier_occupations(5d0,occupations,ok,message)
  call require(ok.and.maxval(abs(occupations-[2d0,2d0,1d0]))<1d-14,&
    'overlapping-Wannier fractional occupation assignment')
  call assign_dg_overlapping_wannier_occupations(7d0,occupations,ok,message)
  call require(.not.ok,'overlapping-Wannier occupation capacity gate')
  call verify_dg_uniform_fragment_target_rank(comm,4,ok,message)
  call require(ok,'uniform fragment target rank accepted')
  if(nproc>1)then
    call verify_dg_uniform_fragment_target_rank(comm,4+merge(1,0,rank==nproc-1),ok,message)
    call require(.not.ok,'nonuniform fragment target rank rejected collectively')
  endif
  theta=0.17d0*rank
  gauge_values=reshape([cmplx(cos(theta),0d0,8),cmplx(sin(theta),0d0,8),&
    cmplx(-sin(theta),0d0,8),cmplx(cos(theta),0d0,8)],[2,2])
  gauge_gradients=0d0
  do i=1,3;gauge_gradients(i,:,:)=i*gauge_values;enddo
  gauge_weights=1d0
  call align_dg_fragment_wannier_gauge(comm,gauge_weights,gauge_values,gauge_gradients,1d-12,&
    gauge_residual,gauge_correction,ok,message)
  call require(ok,trim(message))
  call require(gauge_residual<1d-12,'fragment arbitrary-gauge alignment')
  call require(maxval(abs(gauge_values-identity2()))<1d-12,'fragment aligned reference values')
  fragment_id=nproc-rank
  gauge_values=real(fragment_id,8)*identity2()
  do i=1,3;gauge_gradients(i,:,:)=real(i*fragment_id,8)*identity2();enddo
  call replicate_dg_fragment_wannier_representative(comm,fragment_id,gauge_values,gauge_gradients,&
    gauge_residual,gauge_correction,ok,message)
  call require(ok,trim(message))
  call require(gauge_residual<1d-12,'representative fragment replication closure')
  call require(maxval(abs(gauge_values-identity2()))<1d-12,'representative values replicated')
  do i=1,3
    call require(maxval(abs(gauge_gradients(i,:,:)-i*identity2()))<1d-12,&
      'representative gradients replicated')
  enddo
  local_centers=[1_8,2_8]
  global_centers(:,1)=[3_8,4_8];global_centers(:,2)=[1_8,2_8]
  center_orbit_map(:,1)=[1_8,2_8,3_8,4_8]
  center_orbit_map(:,2)=[3_8,4_8,1_8,2_8]
  call verify_dg_fragment_center_orbit(local_centers,reshape(global_centers,[4]),&
    center_orbit_map,ok,message)
  call require(ok,'reversed-rank translated center orbit closure')
  global_centers(1,1)=2_8
  call verify_dg_fragment_center_orbit(local_centers,reshape(global_centers,[4]),&
    center_orbit_map,ok,message)
  call require(.not.ok,'broken translated center orbit rejected')
  if(nproc>1)then
    allocate(mismatch_values(2+mod(rank,2),2),mismatch_gradients(3,2+mod(rank,2),2))
    mismatch_values=1d0;mismatch_gradients=1d0
    call replicate_dg_fragment_wannier_representative(comm,fragment_id,mismatch_values,&
      mismatch_gradients,gauge_residual,gauge_correction,ok,message)
    call require(.not.ok,'rank-inconsistent representative shape rejected collectively')
    deallocate(mismatch_values,mismatch_gradients)
  endif
  allocate(stream_values(2*nproc,2),stream_gradients(3,2*nproc,2),stream_ids(2),&
    stream_map(2,nproc))
  do i=1,2*nproc
    stream_values(i,:)=[cmplx(modulo(i-1,2)+1,0d0,8),cmplx(-modulo(i-1,2)-1,0d0,8)]
  enddo
  do i=1,3;stream_gradients(i,:,:)=i*stream_values;enddo
  stream_ids=[int(2*rank+1,8),int(2*rank+2,8)]
  do i=1,nproc
    stream_map(:,i)=[int(2*modulo(rank+i-1,nproc)+1,8),int(2*modulo(rank+i-1,nproc)+2,8)]
  enddo
  call verify_dg_fragment_wannier_streaming_closure(comm,rank+1,2,stream_ids,stream_map,&
    stream_values,stream_gradients,1d-12,closure_residual,closure_fingerprint,ok,message)
  call require(ok,trim(message))
  call require(closure_residual<1d-12,'streaming fragment value-gradient closure')
  stream_values=stream_values+cmplx(1d-13,-1d-13,8)
  call verify_dg_fragment_wannier_streaming_closure(comm,rank+1,2,stream_ids,stream_map,&
    stream_values,stream_gradients,1d-12,closure_residual,rounded_closure_fingerprint,ok,message)
  call require(ok,trim(message))
  call require(rounded_closure_fingerprint==closure_fingerprint,&
    'symmetry fingerprint ignores sub-tolerance floating-point roundoff')
  deallocate(stream_values,stream_gradients,stream_ids,stream_map)
  if(nproc==2)then
    allocate(distributed_candidate(1,2),distributed_weight(2),distributed_map(2,2))
    distributed_candidate(1,:)=[cmplx(1+2*rank,0d0,8),cmplx(2+2*rank,0d0,8)]
    distributed_weight=1d0
    distributed_map(:,1)=[int(2*rank+1,8),int(2*rank+2,8)]
    distributed_map(:,2)=[int(2*(1-rank)+1,8),int(2*(1-rank)+2,8)]
    call assemble_dg_distributed_candidate_symmetry(comm,distributed_candidate,distributed_weight,&
      distributed_map,distributed_overlap,ok,message)
    call require(ok,trim(message))
    call require(all(shape(distributed_overlap)==[2,2,2]),'distributed direct-sum symmetry shape')
    call require(abs(distributed_overlap(rank+1,rank+1,1)-sum(abs(distributed_candidate(1,:))**2))<1d-12,&
      'distributed identity symmetry block')
    call require(abs(distributed_overlap(2-rank,rank+1,2)-&
      dot_product(distributed_candidate(1,:),[cmplx(3-2*rank,0d0,8),cmplx(4-2*rank,0d0,8)]))<1d-12,&
      'distributed translated symmetry block')
    deallocate(distributed_candidate,distributed_weight,distributed_map,distributed_overlap)
  endif
  nlocal=count([(mod(p-1,nproc)==rank,p=1,12)])
  allocate(ids(nlocal),box_ids(nlocal),symmetry_map(nlocal,2),broken_symmetry_map(nlocal,2),&
    fragment(nlocal),weight(nlocal),coordinate(nlocal),boundary(nlocal))
  allocate(core_mask(nlocal))
  allocate(candidate(4,nlocal),gradient(3,4,nlocal),occupied(4,2),periodic_phase(3,nlocal))
  allocate(raw_seed_values(1,nlocal))
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
    raw_seed_values(1,index)=cos(6d0*acos(-1d0)*(dble(p)-0.5d0)/12d0)
    if(boundary(index))then
      raw_seed_values(1,index)=1d-9*raw_seed_values(1,index)
    endif
    do i=1,4
      gradient(:,i,index)=[0.03d0*i*candidate(i,index),-0.02d0*p*candidate(i,index),&
        0.01d0*(i+p)*candidate(i,index)]
    enddo
  enddo
  occupied=(0d0,0d0);occupied(2,1)=1d0;occupied(3,2)=1d0
  base_occupied=occupied
  allocate(direct_small(5,nlocal),direct_small_gradient(3,5,nlocal),&
    direct_small_occupied(5,2))
  direct_small(1:4,:)=candidate;direct_small(5,:)=raw_seed_values(1,:)
  direct_small_gradient(:,1:4,:)=gradient
  do i=1,3
    direct_small_gradient(i,5,:)=0.01d0*i*raw_seed_values(1,:)
  enddo
  direct_small_occupied=(0d0,0d0);direct_small_occupied(1:4,:)=occupied
  call construct_dg_overlapping_wannier_basis(comm,5,3,2,ids,fragment,weight,coordinate,boundary,&
    direct_small,direct_small_gradient,direct_small_occupied,8_8,41,1d-8,1d-7,1d-9,&
    result,ok,message,core_mask,projection_seed_values=raw_seed_values)
  call require(ok,trim(message));call local_projector(result%value,direct_projector)
  call require(result%target_rank==3.and.result%projection_inclusion_residual<1d-9,&
    'raw complete-shell direct sum is retained exactly')
  call release_dg_overlapping_wannier_construction(result)
  allocate(seed_values(3,nlocal))
  allocate(occupied_seed_values(2,nlocal))
  occupied_seed_values(1,:)=aimag(candidate(2,:))
  occupied_seed_values(2,:)=real(candidate(3,:))
  call construct_dg_overlapping_wannier_basis(comm,4,2,2,ids,fragment,weight,coordinate,boundary,&
    candidate,gradient,occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message,core_mask,&
    projection_seed_values=occupied_seed_values)
  call require(ok,trim(message))
  call require(result%projection_inclusion_residual<1d-9,'occupied-complete projector inclusion')
  call release_dg_overlapping_wannier_construction(result)
  seed_values(1,:)=aimag(candidate(2,:));seed_values(2,:)=real(candidate(3,:))
  seed_values(3,:)=real(candidate(4,:))

  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    candidate,gradient,occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message,core_mask,&
    projection_seed_values=seed_values)
  call require(ok,trim(message))
  call require(result%target_rank==3.and.result%occupied_inclusion_residual<1d-9,&
    'frozen occupied plus projector-seed target')
  call require(result%projection_inclusion_residual<1d-9,'complete projector-seed inclusion')
  call local_projector(result%value,seed_projector)
  call release_dg_overlapping_wannier_construction(result)
  seed_values=1d-8*seed_values
  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    candidate,gradient,occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message,core_mask,&
    projection_seed_values=seed_values)
  call require(ok,'projector-seed rank test is scale aware')
  call release_dg_overlapping_wannier_construction(result)
  seed_values=1d8*seed_values
  seed_values(3,:)=seed_values(2,:)
  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    candidate,gradient,occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message,core_mask,&
    projection_seed_values=seed_values)
  call require(.not.ok,'linearly dependent projector complement rejected')
  seed_values(3,:)=real(candidate(4,:))
  seed_values(1,:)=real(candidate(1,:))
  call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
    candidate,gradient,occupied,8_8,41,1d-8,1d-7,1d-9,result,ok,message,core_mask,&
    projection_seed_values=seed_values)
  call require(ok,trim(message))
  call require(result%target_rank==4,&
    'complete shell excess residual expands the occupied direct-sum target')
  call require(result%occupied_inclusion_residual<1d-9.and.&
    result%projection_inclusion_residual<1d-9,&
    'expanded direct-sum target exactly includes occupied and complete shell')
  call release_dg_overlapping_wannier_construction(result)
  if(nproc==1)then
    core_mask=.true.
    call construct_dg_overlapping_wannier_basis(comm,4,3,2,ids,fragment,weight,coordinate,boundary,&
      candidate,gradient,occupied,12_8,41,1d-8,1d-7,1d-9,result,ok,message,core_mask,&
      box_ids,symmetry_map(:,1:1),12_8,1d-9,periodic_phase,candidate_axis_offset=0,&
      projection_seed_values=seed_values)
    call require(ok,trim(message))
    call require(result%target_rank==4,&
      'full-candidate direct sum is trivially closed under distributed symmetry')
    call release_dg_overlapping_wannier_construction(result)
    core_mask=box_ids>=3_8.and.box_ids<=10_8
  endif
  seed_values(1,:)=aimag(candidate(2,:))

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
    projection_seed_values=seed_values)
  call require(ok,trim(message));call local_projector(result%value,projector)
  call require(maxval(abs(projector-seed_projector))<1d-9,&
    'projector-seed target is candidate-gauge invariant')
  call release_dg_overlapping_wannier_construction(result)
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
  function identity2() result(identity)
    complex(8)::identity(2,2)
    identity=(0d0,0d0);identity(1,1)=1d0;identity(2,2)=1d0
  end function
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
