#include "config.h"
program test_dg_hybrid_divided_operator_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_variational_payload,only:s_dg_hybrid_fixed_payload,freeze_dg_hybrid_variational_payload
  use dg_hybrid_wannier_complement,only:compute_dg_hybrid_union_to_complete_binding
  use dg_hybrid_divided_operator,only:extract_dg_hybrid_fragment_self_block,&
    compose_dg_hybrid_complete_rows,dg_hybrid_fragment_directory_fingerprint,freeze_dg_hybrid_single_owner_payload
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  implicit none
  integer,parameter::nunion=6,ncomplete=5
  real(real64),parameter::metric_tolerance=1d-8,discarded_metric_eigenvalue=1d-10,&
    above_cutoff_metric_eigenvalue=1d-7
  integer::comm,fragment_comm,rank,nproc,ierr,i,p,local_fragment
  integer::basis_fragment(nunion),basis_local_slot(nunion),basis_generation(nunion),mutated_basis_local_slot(nunion)
  integer(int64),allocatable::union_row_ids(:),output_row_ids(:),bad_row_ids(:)
  complex(real64),allocatable::metric_rows(:,:),kinetic_rows(:,:),nonlocal_rows(:,:),&
    interface_rows(:,:),local_rows(:,:),above_cutoff_metric_rows(:,:),hamiltonian_rows(:,:),&
    complete_metric_rows(:,:),full_rank_metric_rows(:,:)
  complex(real64)::dense_metric(nunion,nunion),dense_kinetic(nunion,nunion),&
    dense_nonlocal(nunion,nunion),dense_interface(nunion,nunion),dense_local(nunion,nunion),&
    dense_hamiltonian(nunion,nunion),dense_full_rank_metric(nunion,nunion),&
    swapped_hamiltonian(nunion,nunion),identity_transform(nunion,nunion),&
    rectangular_transform(nunion,ncomplete),expected_complete_h(ncomplete,ncomplete),&
    expected_complete_s(ncomplete,ncomplete),terminal_projector(nunion,nunion),&
    terminal_column_gram(ncomplete,ncomplete),terminal_identity(ncomplete,ncomplete),&
    metric_pseudoinverse(nunion,nunion),wrong_span_transform(nunion,ncomplete),&
    nonorthogonal_transform(nunion,ncomplete),phase_rotated_transform(nunion,ncomplete),&
    full_rank_nonidentity(nunion,nunion)
  complex(real64),allocatable::swapped_local_rows(:,:),oversized_transform(:,:),rank_mismatched_transform(:,:)
  real(real64)::root_two,metric_scale,metric_cutoff,roundoff_floor,discarded_span_defect
  integer(int64)::identity_fingerprint,rectangular_fingerprint,swapped_fingerprint,failure_fingerprint,&
    basis_directory_fingerprint,identity_binding_fingerprint,rectangular_binding_fingerprint,&
    candidate_binding_fingerprint
  type(s_dg_hybrid_fixed_payload)::fixed_payload,full_rank_payload,mutated_payload,above_cutoff_payload
  logical::ok,near_null_identity_rejected,stale_transform_binding_rejected
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call check_single_owner_publication
  call build_reference_payload(dense_metric,dense_kinetic,dense_nonlocal,dense_interface,dense_local)
  dense_hamiltonian=dense_kinetic+dense_nonlocal+dense_interface+dense_local
  basis_fragment=[1,2,1,2,1,2]
  basis_local_slot=[3,2,1,3,2,1]
  basis_generation=17
  basis_directory_fingerprint=dg_hybrid_fragment_directory_fingerprint(&
    basis_fragment,basis_local_slot,basis_generation,1101_int64)
  call distribute_source_rows(nunion,union_row_ids)
  allocate(metric_rows(size(union_row_ids),nunion),kinetic_rows(size(union_row_ids),nunion),&
    nonlocal_rows(size(union_row_ids),nunion),interface_rows(size(union_row_ids),nunion),&
    local_rows(size(union_row_ids),nunion))
  call select_rows(dense_metric,union_row_ids,metric_rows)
  call select_rows(dense_kinetic,union_row_ids,kinetic_rows)
  call select_rows(dense_nonlocal,union_row_ids,nonlocal_rows)
  call select_rows(dense_interface,union_row_ids,interface_rows)
  call select_rows(dense_local,union_row_ids,local_rows)
  call freeze_dg_hybrid_variational_payload(comm,nunion,union_row_ids,metric_rows,kinetic_rows,&
    nonlocal_rows,interface_rows,1101_int64,1102_int64,1103_int64,fixed_payload,ok,message,&
    basis_directory_fingerprint=basis_directory_fingerprint)
  call require(ok,'fixed payload setup failed: '//trim(message))
  dense_full_rank_metric=dense_metric
  dense_full_rank_metric(1,2)=cmplx(0d0,-0.2d0,real64)
  dense_full_rank_metric(2,1)=conjg(dense_full_rank_metric(1,2))
  allocate(full_rank_metric_rows(size(union_row_ids),nunion))
  call select_rows(dense_full_rank_metric,union_row_ids,full_rank_metric_rows)
  call require(1.1d0>metric_tolerance*1.7d0,&
    'identity oracle metric is not unambiguously full rank at the Task 5 cutoff')
  call freeze_dg_hybrid_variational_payload(comm,nunion,union_row_ids,full_rank_metric_rows,kinetic_rows,&
    nonlocal_rows,interface_rows,1101_int64,1302_int64,1103_int64,full_rank_payload,ok,message,&
    basis_directory_fingerprint=basis_directory_fingerprint)
  call require(ok,'full-rank identity-oracle payload setup failed: '//trim(message))

  if(nproc==1)then
    call check_fragment_self_block(comm,1,[3,5,1])
    call check_fragment_self_block(comm,2,[6,2,4])
  else
    local_fragment=mod(rank,2)+1
    call MPI_Comm_split(comm,local_fragment,rank,fragment_comm,ierr)
    call require(ierr==MPI_SUCCESS,'fragment communicator split failed')
    if(local_fragment==1)then
      call check_fragment_self_block(fragment_comm,1,[3,5,1])
    else
      call check_fragment_self_block(fragment_comm,2,[6,2,4])
    endif
  endif
  call require(abs(dense_nonlocal(1,6))>1d-12.and.abs(dense_interface(1,6))<1d-15.and.&
    abs(dense_interface(3,4))>1d-12,'fixture does not distinguish remote projector and face support')
  call require(abs(dense_hamiltonian(1,6)-dense_nonlocal(1,6))<1d-15.and.&
    abs(dense_hamiltonian(3,4)-dense_interface(3,4))<1d-15,&
    'remote projector or Cartesian-face term is not independently identifiable')

  mutated_basis_local_slot=basis_local_slot
  mutated_basis_local_slot(3)=2;mutated_basis_local_slot(5)=1
  call distribute_rows(3,bad_row_ids)
  call extract_dg_hybrid_fragment_self_block(comm,1,bad_row_ids,basis_fragment,fixed_payload,local_rows,&
    hamiltonian_rows,complete_metric_rows,ok,message,basis_local_slot=mutated_basis_local_slot,&
    basis_generation=basis_generation,fragment_catalog_fingerprint=1101_int64,&
    fragment_directory_fingerprint=basis_directory_fingerprint)
  call require(.not.ok.and.index(message,'directory fingerprint')>0.and.&
    .not.allocated(hamiltonian_rows).and..not.allocated(complete_metric_rows),&
    'stale fragment tuple provenance accepted a reordered local-slot directory')

  identity_transform=(0d0,0d0)
  do i=1,nunion;identity_transform(i,i)=1d0;enddo
  call make_transform_binding(identity_transform,2101_int64,identity_binding_fingerprint)
  call distribute_rows(nunion,output_row_ids)
  call compose_dg_hybrid_complete_rows(comm,output_row_ids,full_rank_payload,local_rows,identity_transform,&
    hamiltonian_rows,complete_metric_rows,identity_fingerprint,ok,message,&
    metric_tolerance=metric_tolerance,complete_map_fingerprint=2101_int64,complete_map_rank=nunion,&
    complete_transform_binding_fingerprint=identity_binding_fingerprint)
  call require(ok,'identity complete composition failed: '//trim(message))
  call compare_distributed_rows(output_row_ids,hamiltonian_rows,dense_hamiltonian,&
    'identity Hamiltonian did not reconstruct every union contribution exactly once')
  call compare_distributed_rows(output_row_ids,complete_metric_rows,dense_full_rank_metric,&
    'identity metric did not reconstruct the union metric')
  call check_global_pencil(output_row_ids,hamiltonian_rows,complete_metric_rows,nunion,.true.)

  call compose_dg_hybrid_complete_rows(comm,output_row_ids,fixed_payload,local_rows,identity_transform,&
    hamiltonian_rows,complete_metric_rows,failure_fingerprint,ok,message,&
    metric_tolerance=metric_tolerance,complete_map_fingerprint=2101_int64,complete_map_rank=nunion,&
    complete_transform_binding_fingerprint=identity_binding_fingerprint)
  near_null_identity_rejected=.not.ok.and.&
    (index(message,'rank')>0.or.index(message,'cutoff')>0).and.failure_fingerprint==0_int64.and.&
    .not.allocated(hamiltonian_rows).and..not.allocated(complete_metric_rows)

  allocate(swapped_local_rows,source=local_rows);swapped_hamiltonian=dense_hamiltonian
  do p=1,size(union_row_ids)
    if(union_row_ids(p)==3_int64)swapped_local_rows(p,3)=dense_local(4,4)
    if(union_row_ids(p)==4_int64)swapped_local_rows(p,4)=dense_local(3,3)
  enddo
  swapped_hamiltonian(3,3)=dense_hamiltonian(3,3)-dense_local(3,3)+dense_local(4,4)
  swapped_hamiltonian(4,4)=dense_hamiltonian(4,4)-dense_local(4,4)+dense_local(3,3)
  call compose_dg_hybrid_complete_rows(comm,output_row_ids,full_rank_payload,swapped_local_rows,identity_transform,&
    hamiltonian_rows,complete_metric_rows,swapped_fingerprint,ok,message,&
    metric_tolerance=metric_tolerance,complete_map_fingerprint=2101_int64,complete_map_rank=nunion,&
    complete_transform_binding_fingerprint=identity_binding_fingerprint)
  call require(ok,'swapped-value fingerprint fixture composition failed: '//trim(message))
  call compare_distributed_rows(output_row_ids,hamiltonian_rows,swapped_hamiltonian,&
    'swapped-value fingerprint fixture produced the wrong Hamiltonian')
  call require(swapped_fingerprint/=identity_fingerprint,&
    'operator fingerprint did not bind values to their matrix positions')
  call check_long_period_fingerprint_collision

  root_two=sqrt(2d0);rectangular_transform=(0d0,0d0)
  rectangular_transform(1,1)=1d0/root_two
  rectangular_transform(2,1)=cmplx(0d0,1d0/root_two,real64)
  rectangular_transform(3,2)=1d0;rectangular_transform(4,3)=1d0
  rectangular_transform(5,4)=1d0;rectangular_transform(6,5)=1d0
  call make_transform_binding(rectangular_transform,2102_int64,rectangular_binding_fingerprint)
  terminal_identity=(0d0,0d0)
  do i=1,ncomplete;terminal_identity(i,i)=1d0;enddo
  terminal_column_gram=matmul(conjg(transpose(rectangular_transform)),rectangular_transform)
  terminal_projector=matmul(rectangular_transform,conjg(transpose(rectangular_transform)))
  metric_scale=2.6d0-discarded_metric_eigenvalue
  metric_cutoff=metric_tolerance*metric_scale
  roundoff_floor=64d0*epsilon(1d0)*metric_scale*real(nunion,real64)
  call require(roundoff_floor<discarded_metric_eigenvalue.and.discarded_metric_eigenvalue<metric_cutoff,&
    'metric fixture does not contain an unambiguous discarded near-null mode')
  metric_pseudoinverse=(0d0,0d0)
  metric_pseudoinverse(1,1)=0.5d0/metric_scale
  metric_pseudoinverse(2,2)=0.5d0/metric_scale
  metric_pseudoinverse(1,2)=cmplx(0d0,-0.5d0/metric_scale,real64)
  metric_pseudoinverse(2,1)=conjg(metric_pseudoinverse(1,2))
  do i=3,nunion;metric_pseudoinverse(i,i)=1d0/real(11+i,real64)*10d0;enddo
  call require(maxval(abs(terminal_column_gram-terminal_identity))<2d-14,&
    'rectangular terminal map does not have orthonormal columns')
  call require(maxval(abs(terminal_projector-matmul(dense_metric,metric_pseudoinverse)))<2d-13,&
    'rectangular terminal map does not equal the retained metric-range projector')
  discarded_span_defect=maxval(abs(matmul(terminal_projector,dense_metric)-dense_metric))
  call require(discarded_span_defect>0.4d0*discarded_metric_eigenvalue.and.&
    discarded_span_defect<0.6d0*discarded_metric_eigenvalue,&
    'metric fixture discarded an exact null rather than a finite near-null mode')
  expected_complete_h=matmul(conjg(transpose(rectangular_transform)),&
    matmul(dense_hamiltonian,rectangular_transform))
  expected_complete_s=matmul(conjg(transpose(rectangular_transform)),&
    matmul(dense_metric,rectangular_transform))
  call distribute_rows(ncomplete,output_row_ids)
  call compose_dg_hybrid_complete_rows(comm,output_row_ids,fixed_payload,local_rows,rectangular_transform,&
    hamiltonian_rows,complete_metric_rows,rectangular_fingerprint,ok,message,&
    metric_tolerance=metric_tolerance,complete_map_fingerprint=2102_int64,complete_map_rank=ncomplete,&
    complete_transform_binding_fingerprint=rectangular_binding_fingerprint)
  call require(ok,'rectangular complete composition failed: '//trim(message))
  call check_exact_null_composition
  call compare_distributed_rows(output_row_ids,hamiltonian_rows,expected_complete_h,&
    'rectangular Hamiltonian congruence is incorrect')
  call compare_distributed_rows(output_row_ids,complete_metric_rows,expected_complete_s,&
    'rectangular metric congruence is incorrect')
  call check_global_pencil(output_row_ids,hamiltonian_rows,complete_metric_rows,ncomplete,.true.)
  call require(identity_fingerprint/=0_int64.and.rectangular_fingerprint/=0_int64.and.&
    identity_fingerprint/=rectangular_fingerprint,'complete operator fingerprints are incomplete')

  phase_rotated_transform=rectangular_transform
  phase_rotated_transform(:,1)=cmplx(0d0,1d0,real64)*phase_rotated_transform(:,1)
  call compose_dg_hybrid_complete_rows(comm,output_row_ids,fixed_payload,local_rows,phase_rotated_transform,&
    hamiltonian_rows,complete_metric_rows,failure_fingerprint,ok,message,&
    metric_tolerance=metric_tolerance,complete_map_fingerprint=2102_int64,complete_map_rank=ncomplete,&
    complete_transform_binding_fingerprint=rectangular_binding_fingerprint)
  stale_transform_binding_rejected=.not.ok.and.&
    (index(message,'fingerprint')>0.or.index(message,'binding')>0).and.failure_fingerprint==0_int64.and.&
    .not.allocated(hamiltonian_rows).and..not.allocated(complete_metric_rows)
  if(rank==0.and..not.near_null_identity_rejected)&
    write(0,'(a)')'RED: near-null identity map bypassed Task 5 cutoff-rank authentication'
  if(rank==0.and..not.stale_transform_binding_rejected)&
    write(0,'(a)')'RED: phase-rotated terminal map reused a stale transform binding receipt'
  call require(near_null_identity_rejected,&
    'terminal identity map accepted a mode below the certified Task 5 cutoff or published failure outputs')
  call require(stale_transform_binding_rejected,&
    'phase-rotated terminal map accepted the original transform binding receipt or published failure outputs')

  call require(above_cutoff_metric_eigenvalue>&
    metric_tolerance*(2.6d0-above_cutoff_metric_eigenvalue),&
    'above-cutoff metric fixture does not increase the certified Task 5 rank')
  allocate(above_cutoff_metric_rows,source=metric_rows)
  do p=1,size(union_row_ids)
    if(union_row_ids(p)==1_int64)above_cutoff_metric_rows(p,2)=&
      cmplx(0d0,-(1.3d0-above_cutoff_metric_eigenvalue),real64)
    if(union_row_ids(p)==2_int64)above_cutoff_metric_rows(p,1)=&
      cmplx(0d0,1.3d0-above_cutoff_metric_eigenvalue,real64)
  enddo
  call freeze_dg_hybrid_variational_payload(comm,nunion,union_row_ids,above_cutoff_metric_rows,kinetic_rows,&
    nonlocal_rows,interface_rows,1101_int64,1202_int64,1103_int64,above_cutoff_payload,ok,message,&
    basis_directory_fingerprint=basis_directory_fingerprint)
  call require(ok,'above-cutoff metric payload setup failed: '//trim(message))
  call make_transform_binding(rectangular_transform,2105_int64,candidate_binding_fingerprint)
  call compose_dg_hybrid_complete_rows(comm,output_row_ids,above_cutoff_payload,local_rows,rectangular_transform,&
    hamiltonian_rows,complete_metric_rows,failure_fingerprint,ok,message,&
    metric_tolerance=metric_tolerance,complete_map_fingerprint=2105_int64,complete_map_rank=ncomplete,&
    complete_transform_binding_fingerprint=candidate_binding_fingerprint)
  call require(.not.ok.and.index(message,'rank')>0.and.failure_fingerprint==0_int64.and.&
    .not.allocated(hamiltonian_rows).and..not.allocated(complete_metric_rows),&
    'terminal map accepted a metric mode above the certified Task 5 cutoff')

  ! The terminal rectangular map is never allowed to feed back into fragment ownership.
  if(nproc==1)then
    call check_fragment_self_block(comm,1,[3,5,1])
    call check_fragment_self_block(comm,2,[6,2,4])
  elseif(local_fragment==1)then
    call check_fragment_self_block(fragment_comm,1,[3,5,1])
  else
    call check_fragment_self_block(fragment_comm,2,[6,2,4])
  endif

  if(allocated(bad_row_ids))deallocate(bad_row_ids)
  if(rank==0)then
    allocate(bad_row_ids(2));bad_row_ids=[1_int64,1_int64]
  else
    allocate(bad_row_ids(0))
  endif
  call compose_dg_hybrid_complete_rows(comm,bad_row_ids,fixed_payload,local_rows,rectangular_transform,&
    hamiltonian_rows,complete_metric_rows,failure_fingerprint,ok,message,&
    metric_tolerance=metric_tolerance,complete_map_fingerprint=2102_int64,complete_map_rank=ncomplete,&
    complete_transform_binding_fingerprint=rectangular_binding_fingerprint)
  call require(.not.ok.and.index(message,'output rows')>0,&
    'duplicate/incomplete complete-row ownership was accepted')

  if(nproc>1)then
    if(rank==0)then
      allocate(rank_mismatched_transform(nunion,ncomplete))
    else
      allocate(rank_mismatched_transform(nunion,ncomplete-1))
    endif
    rank_mismatched_transform=(0d0,0d0)
    call distribute_rows(size(rank_mismatched_transform,2),output_row_ids)
    call compose_dg_hybrid_complete_rows(comm,output_row_ids,fixed_payload,local_rows,rank_mismatched_transform,&
      hamiltonian_rows,complete_metric_rows,failure_fingerprint,ok,message,&
      metric_tolerance=metric_tolerance,complete_map_fingerprint=2102_int64,complete_map_rank=ncomplete,&
      complete_transform_binding_fingerprint=rectangular_binding_fingerprint)
    call require(.not.ok.and.index(message,'shape differs')>0.and.failure_fingerprint==0_int64.and.&
      .not.allocated(hamiltonian_rows).and..not.allocated(complete_metric_rows),&
      'rank-disagreeing terminal transform shape was not rejected collectively')
  endif

  wrong_span_transform=(0d0,0d0)
  wrong_span_transform(1,1)=1d0
  do i=2,ncomplete;wrong_span_transform(i+1,i)=1d0;enddo
  call make_transform_binding(wrong_span_transform,2102_int64,candidate_binding_fingerprint)
  call distribute_rows(ncomplete,output_row_ids)
  call compose_dg_hybrid_complete_rows(comm,output_row_ids,fixed_payload,local_rows,wrong_span_transform,&
    hamiltonian_rows,complete_metric_rows,failure_fingerprint,ok,message,&
    metric_tolerance=metric_tolerance,complete_map_fingerprint=2102_int64,complete_map_rank=ncomplete,&
    complete_transform_binding_fingerprint=candidate_binding_fingerprint)
  call require(.not.ok.and.index(message,'span')>0.and.failure_fingerprint==0_int64.and.&
    .not.allocated(hamiltonian_rows).and..not.allocated(complete_metric_rows),&
    'orthonormal positive-metric transform with the wrong retained span was accepted')

  nonorthogonal_transform=rectangular_transform;nonorthogonal_transform(:,1)=2d0*nonorthogonal_transform(:,1)
  call make_transform_binding(nonorthogonal_transform,2102_int64,candidate_binding_fingerprint)
  call compose_dg_hybrid_complete_rows(comm,output_row_ids,fixed_payload,local_rows,nonorthogonal_transform,&
    hamiltonian_rows,complete_metric_rows,failure_fingerprint,ok,message,&
    metric_tolerance=metric_tolerance,complete_map_fingerprint=2102_int64,complete_map_rank=ncomplete,&
    complete_transform_binding_fingerprint=candidate_binding_fingerprint)
  call require(.not.ok.and.index(message,'orthonormal')>0.and.failure_fingerprint==0_int64.and.&
    .not.allocated(hamiltonian_rows).and..not.allocated(complete_metric_rows),&
    'nonorthonormal terminal transform was accepted')

  allocate(oversized_transform(nunion,nunion+1));oversized_transform=(0d0,0d0)
  do i=1,nunion;oversized_transform(i,i)=1d0;enddo
  call make_transform_binding(oversized_transform,2103_int64,candidate_binding_fingerprint)
  call distribute_rows(nunion+1,output_row_ids)
  call compose_dg_hybrid_complete_rows(comm,output_row_ids,fixed_payload,local_rows,oversized_transform,&
    hamiltonian_rows,complete_metric_rows,failure_fingerprint,ok,message,&
    metric_tolerance=metric_tolerance,complete_map_fingerprint=2103_int64,complete_map_rank=nunion+1,&
    complete_transform_binding_fingerprint=candidate_binding_fingerprint)
  call require(.not.ok.and.index(message,'transform')>0.and.failure_fingerprint==0_int64.and.&
    .not.allocated(hamiltonian_rows).and..not.allocated(complete_metric_rows),&
    'terminal transform wider than the union was accepted')

  full_rank_nonidentity=(0d0,0d0)
  do i=1,nunion;full_rank_nonidentity(i,i)=1d0;enddo
  full_rank_nonidentity(1,1)=0d0;full_rank_nonidentity(2,2)=0d0
  full_rank_nonidentity(1,2)=1d0;full_rank_nonidentity(2,1)=1d0
  call make_transform_binding(full_rank_nonidentity,2104_int64,candidate_binding_fingerprint)
  call distribute_rows(nunion,output_row_ids)
  call compose_dg_hybrid_complete_rows(comm,output_row_ids,fixed_payload,local_rows,full_rank_nonidentity,&
    hamiltonian_rows,complete_metric_rows,failure_fingerprint,ok,message,&
    metric_tolerance=metric_tolerance,complete_map_fingerprint=2104_int64,complete_map_rank=nunion,&
    complete_transform_binding_fingerprint=candidate_binding_fingerprint)
  call require(.not.ok.and.index(message,'full-rank')>0.and.failure_fingerprint==0_int64.and.&
    .not.allocated(hamiltonian_rows).and..not.allocated(complete_metric_rows),&
    'nonidentity full-rank terminal transform was accepted')

  mutated_payload=fixed_payload
  do p=1,size(union_row_ids)
    mutated_payload%kinetic_rows(p,int(union_row_ids(p)))=&
      mutated_payload%kinetic_rows(p,int(union_row_ids(p)))+1d-3
  enddo
  call distribute_rows(3,bad_row_ids)
  call extract_dg_hybrid_fragment_self_block(comm,1,bad_row_ids,basis_fragment,mutated_payload,local_rows,&
    hamiltonian_rows,complete_metric_rows,ok,message,basis_local_slot=basis_local_slot,&
    basis_generation=basis_generation,fragment_catalog_fingerprint=1101_int64,&
    fragment_directory_fingerprint=basis_directory_fingerprint)
  call require(.not.ok.and.index(message,'fingerprint')>0.and..not.allocated(hamiltonian_rows).and.&
    .not.allocated(complete_metric_rows),'fragment extraction accepted a mutated immutable payload')
  call distribute_rows(ncomplete,output_row_ids)
  call compose_dg_hybrid_complete_rows(comm,output_row_ids,mutated_payload,local_rows,rectangular_transform,&
    hamiltonian_rows,complete_metric_rows,failure_fingerprint,ok,message,&
    metric_tolerance=metric_tolerance,complete_map_fingerprint=2102_int64,complete_map_rank=ncomplete,&
    complete_transform_binding_fingerprint=rectangular_binding_fingerprint)
  call require(.not.ok.and.index(message,'fingerprint')>0,'mutated immutable payload was accepted')

  if(rank==0)then
    write(*,'(a,i0,a,i0,a,i0)')'HYBRID_DIVIDED_OPERATOR ranks=',nproc,&
      ' identity=',identity_fingerprint,' rectangular=',rectangular_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid divided operator on ',nproc,' ranks'
  endif
  if(nproc>1)call MPI_Comm_free(fragment_comm,ierr)
  call MPI_Finalize(ierr)
contains
  subroutine check_single_owner_publication
    type(s_dg_hybrid_fragment_basis)::basis,bad_basis
    type(s_dg_hybrid_fixed_payload)::published
    integer,allocatable::owner(:),fragment(:),slot(:),generation(:)
    complex(real64),allocatable::metric(:,:),kinetic(:,:),zero(:,:)
    integer::f,n,nb,first,a,k,scenario
    integer(int64)::directory_fp,input_fp
    logical::same
    ! One rank per fragment, reversed rank assignment and local column order.
    f=nproc-rank;n=f;nb=nproc*(nproc+1)/2;first=f*(f-1)/2
    basis%fragment_id=f;basis%generation=17;basis%provenance_fingerprint=301_int64+f
    allocate(basis%global_ids(n),basis%sector(n),basis%buffer_point_ids(1),basis%buffer_values(1,n))
    basis%global_ids=[(int(first+n-a+1,int64),a=1,n)]
    basis%sector=1;basis%sector(n)=2;basis%buffer_point_ids=int(f,int64);basis%buffer_values=1d0
    allocate(metric(n,nb),kinetic(n,nb),zero(n,nb));metric=0d0;zero=0d0
    do a=1,n;metric(a,basis%global_ids(a))=1d0;enddo
    kinetic=2d0*metric
    call freeze_dg_hybrid_single_owner_payload(comm,nproc,basis,metric,kinetic,zero,zero,&
      501_int64,503_int64,505_int64,published,owner,fragment,slot,generation,directory_fp,ok,message)
    call require(ok,'single-owner publication failed: '//trim(message))
    same=all(published%row_ids==basis%global_ids).and.all(published%metric_rows==metric).and.&
      all(published%kinetic_rows==kinetic).and.all(generation==17).and.&
      published%basis_directory_fingerprint==directory_fp
    do k=1,nproc
      first=k*(k-1)/2
      do a=1,k
        same=same.and.owner(first+a)==nproc-k.and.fragment(first+a)==k.and.slot(first+a)==k-a+1
      enddo
    enddo
    call require(same,'single-owner publication reordered columns or lost rank/fragment/generation bindings')
    do scenario=1,5
      bad_basis=basis;input_fp=501_int64
      if(rank==0)then
        if(scenario==1)bad_basis%global_ids(1)=0_int64
        if(scenario==2)bad_basis%generation=0
        if(scenario==3)bad_basis%fragment_id=1
        if(scenario==4)input_fp=511_int64
        if(scenario==5)bad_basis%global_ids(1)=1_int64
      endif
      if(nproc==1.and.scenario>=3)cycle
      call freeze_dg_hybrid_single_owner_payload(comm,nproc,bad_basis,metric,kinetic,zero,zero,&
        input_fp,503_int64,505_int64,published,owner,fragment,slot,generation,directory_fp,ok,message)
      call require(.not.ok.and..not.published%frozen.and..not.allocated(owner).and.&
        .not.allocated(fragment).and..not.allocated(slot).and..not.allocated(generation).and.directory_fp==0_int64,&
        'invalid single-owner metadata published a payload or directory')
    enddo
    call freeze_dg_hybrid_single_owner_payload(comm,nproc+1,basis,metric,kinetic,zero,zero,&
      501_int64,503_int64,505_int64,published,owner,fragment,slot,generation,directory_fp,ok,message)
    call require(.not.ok.and..not.published%frozen,'unequal fragment and rank counts accepted')
    call freeze_dg_hybrid_single_owner_payload(comm,nproc,basis,metric(:,:nb-1),kinetic,zero,zero,&
      501_int64,503_int64,505_int64,published,owner,fragment,slot,generation,directory_fp,ok,message)
    call require(.not.ok.and..not.published%frozen.and..not.allocated(owner),&
      'invalid matrix extent published single-owner directory')
  end subroutine check_single_owner_publication

  subroutine check_exact_null_composition
    type(s_dg_hybrid_fixed_payload)::null_payload
    complex(real64)::null_metric(nunion,nunion),expected_s(ncomplete,ncomplete)
    complex(real64),allocatable::null_rows(:,:),h(:,:),s(:,:)
    integer(int64)::binding,receipt
    null_metric=dense_metric
    null_metric(1,2)=cmplx(0d0,-1.3d0,real64);null_metric(2,1)=conjg(null_metric(1,2))
    allocate(null_rows(size(union_row_ids),nunion))
    call select_rows(null_metric,union_row_ids,null_rows)
    call freeze_dg_hybrid_variational_payload(comm,nunion,union_row_ids,null_rows,kinetic_rows,&
      nonlocal_rows,interface_rows,1101_int64,1402_int64,1103_int64,null_payload,ok,message,&
      basis_directory_fingerprint=basis_directory_fingerprint)
    call require(ok,'exact-null payload setup failed: '//trim(message))
    call compute_dg_hybrid_union_to_complete_binding(comm,rectangular_transform,2102_int64,1d-12,&
      binding,ok,message)
    call require(ok,'exact-null transform binding failed: '//trim(message))
    call compose_dg_hybrid_complete_rows(comm,output_row_ids,null_payload,local_rows,rectangular_transform,&
      h,s,receipt,ok,message,metric_tolerance=1d-12,complete_map_fingerprint=2102_int64,&
      complete_map_rank=ncomplete,complete_transform_binding_fingerprint=binding)
    call require(ok,'exact-null terminal metric rejected as ambiguous: '//trim(message))
    expected_s=matmul(conjg(transpose(rectangular_transform)),matmul(null_metric,rectangular_transform))
    call compare_distributed_rows(output_row_ids,s,expected_s,'exact-null complete metric differs from congruence')
    call compare_distributed_rows(output_row_ids,h,expected_complete_h,'exact-null composition changed Hamiltonian')
  end subroutine check_exact_null_composition

  subroutine make_transform_binding(transform,map_fingerprint,binding_fingerprint)
    complex(real64),intent(in)::transform(:,:)
    integer(int64),intent(in)::map_fingerprint
    integer(int64),intent(out)::binding_fingerprint
    logical::binding_ok
    character(256)::binding_message
    call compute_dg_hybrid_union_to_complete_binding(comm,transform,map_fingerprint,metric_tolerance,&
      binding_fingerprint,binding_ok,binding_message)
    call require(binding_ok,'terminal transform binding setup failed: '//trim(binding_message))
  end subroutine make_transform_binding

  subroutine build_reference_payload(metric,kinetic,nonlocal,interface,local)
    complex(real64),intent(out)::metric(:,:),kinetic(:,:),nonlocal(:,:),interface(:,:),local(:,:)
    integer::j
    metric=(0d0,0d0);kinetic=(0d0,0d0);nonlocal=(0d0,0d0)
    interface=(0d0,0d0);local=(0d0,0d0)
    metric(1,1)=1.3d0;metric(2,2)=1.3d0
    metric(1,2)=cmplx(0d0,-(1.3d0-discarded_metric_eigenvalue),real64)
    metric(2,1)=conjg(metric(1,2))
    do j=1,nunion
      if(j>=3)metric(j,j)=1.1d0+0.1d0*j
      kinetic(j,j)=0.4d0+0.07d0*j
      nonlocal(j,j)=-0.03d0+0.004d0*j
      interface(j,j)=0.08d0+0.006d0*j
      local(j,j)=-0.2d0+0.02d0*j
    enddo
    call set_pair(kinetic,1,3,cmplx(0.11d0,0.02d0,real64))
    call set_pair(kinetic,2,6,cmplx(-0.07d0,0.015d0,real64))
    call set_pair(nonlocal,3,5,cmplx(0.035d0,-0.012d0,real64))
    call set_pair(nonlocal,2,4,cmplx(-0.026d0,0.009d0,real64))
    call set_pair(nonlocal,1,6,cmplx(0.045d0,0.017d0,real64))
    call set_pair(interface,1,5,cmplx(-0.028d0,0.011d0,real64))
    call set_pair(interface,3,4,cmplx(-0.09d0,0.025d0,real64))
    call set_pair(interface,4,6,cmplx(0.031d0,-0.008d0,real64))
    call set_pair(local,1,3,cmplx(-0.023d0,0.004d0,real64))
    call set_pair(local,2,6,cmplx(0.019d0,0.007d0,real64))
  end subroutine build_reference_payload

  subroutine set_pair(matrix,irow,icol,value)
    complex(real64),intent(inout)::matrix(:,:)
    integer,intent(in)::irow,icol
    complex(real64),intent(in)::value
    matrix(irow,icol)=value;matrix(icol,irow)=conjg(value)
  end subroutine set_pair

  subroutine distribute_rows(global_count,row_ids)
    integer,intent(in)::global_count
    integer(int64),allocatable,intent(out)::row_ids(:)
    integer::j,q
    allocate(row_ids(count([(mod(j,nproc)==rank,j=1,global_count)])))
    q=0
    do j=1,global_count
      if(mod(j,nproc)/=rank)cycle
      q=q+1;row_ids(q)=int(j,int64)
    enddo
    if(size(row_ids)>1)row_ids=row_ids(size(row_ids):1:-1)
  end subroutine distribute_rows

  subroutine distribute_source_rows(global_count,row_ids)
    integer,intent(in)::global_count
    integer(int64),allocatable,intent(out)::row_ids(:)
    integer::j,q
    allocate(row_ids(count([(mod(j-1,nproc)==rank,j=1,global_count)])))
    q=0
    do j=1,global_count
      if(mod(j-1,nproc)/=rank)cycle
      q=q+1;row_ids(q)=int(j,int64)
    enddo
    if(size(row_ids)>1)row_ids=row_ids(size(row_ids):1:-1)
  end subroutine distribute_source_rows

  subroutine select_rows(dense,row_ids,rows)
    complex(real64),intent(in)::dense(:,:)
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(out)::rows(:,:)
    integer::j
    do j=1,size(row_ids);rows(j,:)=dense(int(row_ids(j)),:);enddo
  end subroutine select_rows

  subroutine check_fragment_self_block(fragment_communicator,fragment,global_ids)
    integer,intent(in)::fragment_communicator,fragment,global_ids(:)
    integer(int64),allocatable::fragment_rows(:),bad_fragment_rows(:)
    complex(real64),allocatable::hff(:,:),sff(:,:)
    real(real64)::hamiltonian_defect,metric_defect
    integer::j
    call distribute_rows_on_comm(fragment_communicator,size(global_ids),fragment_rows)
    call extract_dg_hybrid_fragment_self_block(fragment_communicator,fragment,fragment_rows,basis_fragment,&
      fixed_payload,local_rows,hff,sff,ok,message,basis_local_slot=basis_local_slot,&
      basis_generation=basis_generation,fragment_catalog_fingerprint=1101_int64,&
      fragment_directory_fingerprint=basis_directory_fingerprint)
    call require(ok,'fragment self-block extraction failed: '//trim(message))
    call require(any(shape(hff)/=[size(fragment_rows),size(global_ids)]) .eqv. .false.,&
      'fragment Hamiltonian has the wrong distributed extent')
    call require(any(shape(sff)/=shape(hff)) .eqv. .false.,&
      'fragment metric has the wrong distributed extent')
    hamiltonian_defect=0d0;metric_defect=0d0
    do j=1,size(fragment_rows)
      hamiltonian_defect=max(hamiltonian_defect,maxval(abs(hff(j,:)-dense_hamiltonian(&
        global_ids(int(fragment_rows(j))),global_ids))))
      metric_defect=max(metric_defect,maxval(abs(sff(j,:)-dense_metric(&
        global_ids(int(fragment_rows(j))),global_ids))))
    enddo
    call require(hamiltonian_defect<2d-13,&
      'fragment Hamiltonian omitted a fixed self contribution')
    call require(metric_defect<2d-13,'fragment metric self block is incorrect')
    call MPI_Comm_rank(fragment_communicator,j,ierr)
    if(j==0)then
      allocate(bad_fragment_rows(2));bad_fragment_rows=[1_int64,1_int64]
    else
      allocate(bad_fragment_rows(0))
    endif
    call extract_dg_hybrid_fragment_self_block(fragment_communicator,fragment,bad_fragment_rows,basis_fragment,&
      fixed_payload,local_rows,hff,sff,ok,message,basis_local_slot=basis_local_slot,&
      basis_generation=basis_generation,fragment_catalog_fingerprint=1101_int64,&
      fragment_directory_fingerprint=basis_directory_fingerprint)
    call require(.not.ok.and.index(message,'output rows')>0,&
      'duplicate/incomplete fragment-row ownership was accepted')
  end subroutine check_fragment_self_block

  subroutine check_long_period_fingerprint_collision
    integer,parameter::nlong=64
    integer(int64),allocatable::source_ids(:),result_ids(:)
    complex(real64),allocatable::metric(:,:),kinetic(:,:),nonlocal(:,:),interface(:,:),&
      local_a(:,:),local_b(:,:),transform(:,:),result_h(:,:),result_s(:,:)
    type(s_dg_hybrid_fixed_payload)::payload
    integer(int64)::fingerprint_a,fingerprint_b,transform_binding_fingerprint
    integer::q,row
    logical::stage_ok
    character(256)::stage_message
    call distribute_source_rows(nlong,source_ids)
    allocate(metric(size(source_ids),nlong),kinetic(size(source_ids),nlong),&
      nonlocal(size(source_ids),nlong),interface(size(source_ids),nlong),&
      local_a(size(source_ids),nlong),local_b(size(source_ids),nlong),transform(nlong,nlong))
    metric=(0d0,0d0);kinetic=(0d0,0d0);nonlocal=(0d0,0d0);interface=(0d0,0d0)
    local_a=(0d0,0d0);local_b=(0d0,0d0);transform=(0d0,0d0)
    do q=1,size(source_ids)
      row=int(source_ids(q));metric(q,row)=1d0
      select case(row)
      case(1);local_a(q,2)=0.11d0;local_b(q,2)=0.22d0
      case(2);local_a(q,1)=0.11d0;local_a(q,64)=0.22d0
              local_b(q,1)=0.22d0;local_b(q,64)=0.11d0
      case(64);local_a(q,2)=0.22d0;local_b(q,2)=0.11d0
      end select
    enddo
    do q=1,nlong;transform(q,q)=1d0;enddo
    call freeze_dg_hybrid_variational_payload(comm,nlong,source_ids,metric,kinetic,nonlocal,interface,&
      3101_int64,3102_int64,3103_int64,payload,stage_ok,stage_message)
    call require(stage_ok,'long-period fingerprint payload setup failed: '//trim(stage_message))
    call make_transform_binding(transform,3104_int64,transform_binding_fingerprint)
    call distribute_rows(nlong,result_ids)
    call compose_dg_hybrid_complete_rows(comm,result_ids,payload,local_a,transform,result_h,result_s,&
      fingerprint_a,stage_ok,stage_message,metric_tolerance=metric_tolerance,&
      complete_map_fingerprint=3104_int64,complete_map_rank=nlong,&
      complete_transform_binding_fingerprint=transform_binding_fingerprint)
    call require(stage_ok,'first long-period fingerprint composition failed: '//trim(stage_message))
    call compose_dg_hybrid_complete_rows(comm,result_ids,payload,local_b,transform,result_h,result_s,&
      fingerprint_b,stage_ok,stage_message,metric_tolerance=metric_tolerance,&
      complete_map_fingerprint=3104_int64,complete_map_rank=nlong,&
      complete_transform_binding_fingerprint=transform_binding_fingerprint)
    call require(stage_ok,'second long-period fingerprint composition failed: '//trim(stage_message))
    call require(fingerprint_a/=fingerprint_b,&
      'operator fingerprint has a periodic value-position collision')
  end subroutine check_long_period_fingerprint_collision

  subroutine distribute_rows_on_comm(communicator,global_count,row_ids)
    integer,intent(in)::communicator,global_count
    integer(int64),allocatable,intent(out)::row_ids(:)
    integer::local_rank,local_nproc,j,q
    call MPI_Comm_rank(communicator,local_rank,ierr)
    call MPI_Comm_size(communicator,local_nproc,ierr)
    allocate(row_ids(count([(mod(j,local_nproc)==local_rank,j=1,global_count)])))
    q=0
    do j=1,global_count
      if(mod(j,local_nproc)/=local_rank)cycle
      q=q+1;row_ids(q)=int(j,int64)
    enddo
    if(size(row_ids)>1)row_ids=row_ids(size(row_ids):1:-1)
  end subroutine distribute_rows_on_comm

  subroutine compare_distributed_rows(row_ids,rows,reference,label)
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::rows(:,:),reference(:,:)
    character(*),intent(in)::label
    real(real64)::local_defect,global_defect
    integer::j
    local_defect=0d0
    do j=1,size(row_ids)
      local_defect=max(local_defect,maxval(abs(rows(j,:)-reference(int(row_ids(j)),:))))
    enddo
    call MPI_Allreduce(local_defect,global_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call require(ierr==MPI_SUCCESS.and.global_defect<4d-13,label)
  end subroutine compare_distributed_rows

  subroutine check_global_pencil(row_ids,hrows,srows,n,require_positive)
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::hrows(:,:),srows(:,:)
    integer,intent(in)::n
    logical,intent(in)::require_positive
    complex(real64),allocatable::full_h(:,:),full_s(:,:),factor(:,:)
    complex(real64)::value
    real(real64)::pivot
    integer::j,k,q
    allocate(full_h(n,n),full_s(n,n),factor(n,n));full_h=(0d0,0d0);full_s=(0d0,0d0)
    do j=1,size(row_ids)
      full_h(int(row_ids(j)),:)=hrows(j,:);full_s(int(row_ids(j)),:)=srows(j,:)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,full_h,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call require(ierr==MPI_SUCCESS,'Hamiltonian row collection failed')
    call MPI_Allreduce(MPI_IN_PLACE,full_s,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call require(ierr==MPI_SUCCESS,'metric row collection failed')
    call require(maxval(abs(full_h-conjg(transpose(full_h))))<4d-13,&
      'complete Hamiltonian is not Hermitian')
    call require(maxval(abs(full_s-conjg(transpose(full_s))))<4d-13,&
      'complete metric is not Hermitian')
    if(.not.require_positive)return
    factor=(0d0,0d0)
    do j=1,n
      value=full_s(j,j)
      do q=1,j-1;value=value-factor(j,q)*conjg(factor(j,q));enddo
      pivot=real(value,real64)
      call require(abs(aimag(value))<4d-13.and.pivot>1d-10,&
        'complete metric lost positive retained rank')
      factor(j,j)=sqrt(pivot)
      do k=j+1,n
        value=full_s(k,j)
        do q=1,j-1;value=value-factor(k,q)*conjg(factor(j,q));enddo
        factor(k,j)=value/factor(j,j)
      enddo
    enddo
  end subroutine check_global_pencil

  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::bad,global_bad
    bad=merge(0,1,condition)
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(rank==0)write(0,'(a)')trim(label)
      error stop 1
    endif
  end subroutine require
end program test_dg_hybrid_divided_operator_mpi
