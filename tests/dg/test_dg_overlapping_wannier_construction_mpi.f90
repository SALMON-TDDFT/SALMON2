#include "config.h"
program test_dg_overlapping_wannier_construction_mpi
  use mpi
  use dg_overlapping_wannier_types,only:s_dg_ow_distributed_layout,&
    initialize_dg_ow_distributed_layout,reserve_dg_ow_workspace,&
    release_dg_ow_workspace,release_dg_ow_distributed_layout
  use dg_overlapping_wannier_construction,only:s_dg_overlapping_wannier_construction,&
    construct_dg_overlapping_wannier_basis,release_dg_overlapping_wannier_construction,&
    verify_dg_overlapping_wannier_periodic_closure,assemble_dg_distributed_candidate_symmetry,&
    assemble_dg_distributed_basis_symmetry_overlap,&
    assemble_dg_distributed_basis_symmetry_overlap_rows,&
    build_dg_pointwise_affine_owner_map,&
    find_dg_group_identity,&
    select_dg_fixed_rank_symmetry_closed_subspace,&
    build_dg_distributed_symmetry_closed_basis,&
    accumulate_dg_lcfo_buffer_contributions_to_core,&
    measure_dg_rank_fixed_symmetry_residuals,&
    accept_dg_boundary_calibrated_symmetry,&
    solve_dg_affine_common_fixed_point,&
    compute_dg_periodic_wannier_centers,&
    verify_dg_wannier_center_affine_orbits,&
    align_dg_fragment_wannier_gauge,replicate_dg_fragment_wannier_representative,&
    verify_dg_fragment_wannier_streaming_closure,verify_dg_fragment_center_orbit,&
    verify_dg_uniform_fragment_target_rank,assign_dg_overlapping_wannier_occupations,&
    build_dg_balanced_orbital_ownership,transpose_dg_spatial_cores_to_orbital_owners,&
    exchange_dg_point_permuted_orbital_rows,&
    redistribute_dg_owned_orbitals_to_center_fragments,&
    assign_dg_periodic_centers_to_fragments,&
    verify_dg_fragment_subspace_density_covariance,build_dg_core_owned_occupied_subspace
  implicit none
  integer::comm,rank,nproc,ierr,i,p,point,nlocal,nclosure,index,ncore,fragment_id
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
  complex(8),allocatable::distributed_candidate(:,:),distributed_overlap(:,:,:),&
    distributed_basis(:,:),distributed_basis_overlap(:,:,:),orbit_seed(:,:),required_orbit_seed(:,:),&
    orbit_basis(:,:),&
    orbit_gram(:,:)
  complex(8),allocatable::distributed_basis_overlap_rows(:,:,:)
  complex(8)::lcfo_buffer_contribution(2,2),lcfo_core_value(2,1)
  integer(8)::lcfo_buffer_ids(2),lcfo_core_ids(1)
  integer::orbit_rank,required_orbit_rank,identity_operation
  integer(8),allocatable::distributed_map(:,:),invalid_orbit_map(:,:)
  integer(8),allocatable::affine_local_ids(:),affine_all_ids(:,:),affine_target_ids(:),&
    affine_second_ids(:)
  integer,allocatable::affine_target_owner(:),affine_target_local(:),affine_wrap(:,:),&
    affine_second_owner(:),affine_second_local(:),affine_second_wrap(:,:)
  integer::affine_rotation(3,3)
  real(8)::affine_translation(3)
  real(8),allocatable::distributed_weight(:)
  real(8),allocatable::seed_values(:,:)
  real(8),allocatable::occupied_seed_values(:,:)
  real(8),allocatable::raw_seed_values(:,:)
  real(8)::occupations(3)
  type(s_dg_overlapping_wannier_construction)::result
  type(s_dg_ow_distributed_layout)::distributed_layout
  integer,allocatable::reference_owner(:),reference_center_fragment(:)
  integer(8),allocatable::reference_center_box_ids(:)
  integer(8)::reference_fingerprint
  integer(8)::symmetry_workspace_peak
  integer(8)::row_overlap_workspace_peak
  integer(8),allocatable::distributed_overlap_row_ids(:)
  integer(8)::closure_fingerprint,rounded_closure_fingerprint
  integer(8),allocatable::closure_ids(:),closure_map(:,:)
  integer(8),allocatable::stream_ids(:),stream_map(:,:)
  integer(8)::local_centers(2),global_centers(2,2),center_orbit_map(4,2)
  complex(8),allocatable::closure_values(:,:),closure_gradients(:,:,:)
  complex(8)::closure_representation(2,2,1)
  real(8)::closure_rotation(3,3,1),closure_residual
  real(8)::gauge_residual,gauge_correction,theta
  complex(8)::gauge_values(2,2),gauge_gradients(3,2,2)
  complex(8)::mixed_wannier(2,2)
  complex(8)::closure_metric(6,6),closure_localizer(6,6),closure_group(6,6,2),&
    closure_occupied(6,2),scaled_closure_metric(6,6)
  complex(8),allocatable::closure_transform(:,:)
  integer::closure_product(2,2)
  real(8)::subspace_leakage,occupied_inclusion
  complex(8)::calibrated_basis(1,4),calibrated_representation(1,1,2)
  integer(8)::calibrated_map(4,2)
  logical::calibrated_boundary(4)
  real(8)::calibrated_total(2),calibrated_boundary_residual(2),calibrated_interior_residual(2)
  real(8)::calibrated_allowance,calibrated_strict_tolerance
  integer::affine_rotations(3,3,2)
  real(8)::affine_translations(3,2),affine_center(3),affine_residual
  logical::has_affine_center
  complex(8)::periodic_center_values(1,2),periodic_center_phases(3,2)
  real(8)::periodic_centers(3,1),periodic_center_magnitudes(3,1)
  complex(8),allocatable::transpose_local(:,:),transpose_owned(:,:)
  complex(8),allocatable::permuted_image(:,:)
  complex(8),allocatable::mismatched_owned(:,:)
  complex(8),allocatable::center_local_values(:,:)
  integer(8),allocatable::transpose_local_ids(:),transpose_global_ids(:)
  integer(8),allocatable::mismatched_global_ids(:)
  integer(8),allocatable::redistribution_buffer_ids(:)
  integer(8),allocatable::center_all_core_ids(:,:),assigned_center_ids(:)
  integer,allocatable::orbital_counts(:),orbital_displacements(:),orbital_owners(:),owned_orbitals(:)
  integer,allocatable::invalid_owned_orbitals(:)
  integer,allocatable::mismatched_owned_orbitals(:)
  integer,allocatable::center_owners(:),center_local_orbitals(:)
  integer,allocatable::center_fragments(:),assigned_center_owners(:),assigned_center_fragments(:)
  real(8)::assignment_centers(3,4)
  real(8)::orbit_centers(3,2)
  complex(8)::fractional_core_candidates(2,2)
  complex(8),allocatable::core_occupied_coefficients(:,:)
  integer(8)::mixed_map(2,1)
  real(8)::fractional_core_electrons
  real(8)::gauge_weights(2)
  logical::ok,transpose_values_ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call initialize_dg_ow_distributed_layout(comm,0,2*nproc,distributed_layout,ok,message)
  call require(.not.ok,'distributed layout rejects an empty global orbital space')
  call initialize_dg_ow_distributed_layout(comm,2*nproc+1,2*nproc,&
    distributed_layout,ok,message)
  call require(ok.and.distributed_layout%global_orbital_count==2*nproc+1.and.&
    distributed_layout%global_core_point_count==2*nproc.and.&
    distributed_layout%owned_orbital_count<=2*nproc+1.and.&
    distributed_layout%owned_core_point_count==2,&
    'two-dimensional orbital/core layout has bounded unique ownership')
  call require(count(distributed_layout%orbital_owner==rank)==&
    distributed_layout%owned_orbital_count.and.&
    count(distributed_layout%core_point_owner==rank)==&
    distributed_layout%owned_core_point_count.and.&
    all(distributed_layout%orbital_owner>=0).and.&
    all(distributed_layout%orbital_owner<nproc).and.&
    all(distributed_layout%core_point_owner>=0).and.&
    all(distributed_layout%core_point_owner<nproc),&
    'distributed layout assigns every orbital and core point exactly once')
  if(nproc>1)call require(distributed_layout%owned_orbital_count<2*nproc+1,&
    'distributed layout does not replicate every orbital')
  call reserve_dg_ow_workspace(distributed_layout,4096_8,ok,message)
  call require(ok.and.distributed_layout%current_workspace_bytes==4096_8.and.&
    distributed_layout%peak_workspace_bytes==4096_8,&
    'workspace accounting records current and peak bytes')
  call release_dg_ow_workspace(distributed_layout,4096_8,ok,message)
  call require(ok.and.distributed_layout%current_workspace_bytes==0_8.and.&
    distributed_layout%peak_workspace_bytes==4096_8,&
    'operation-local workspace release preserves peak receipt')
  call reserve_dg_ow_workspace(distributed_layout,huge(0_8),ok,message)
  call require(ok,'workspace accounting accepts the largest representable extent')
  call reserve_dg_ow_workspace(distributed_layout,1_8,ok,message)
  call require(.not.ok,'workspace accounting rejects extent overflow')
  call release_dg_ow_distributed_layout(distributed_layout)
  affine_rotations=0
  affine_rotations(:,:,1)=reshape([1,0,0,0,1,0,0,0,1],[3,3])
  affine_rotations(:,:,2)=reshape([-1,0,0,0,-1,0,0,0,1],[3,3])
  affine_translations=0d0;affine_translations(:,2)=[0.5d0,0.25d0,0d0]
  call solve_dg_affine_common_fixed_point(affine_rotations,affine_translations,1d-12,&
    has_affine_center,affine_center,affine_residual,ok,message)
  call require(ok.and.has_affine_center.and.affine_residual<1d-12.and.&
    maxval(abs(modulo(affine_center(1:2)-[0.25d0,0.125d0]+0.5d0,1d0)-0.5d0))<1d-12,&
    'non-origin noncentrosymmetric rotation center is solved from full affine operations')
  affine_translations(:,2)=[0d0,0d0,0.5d0]
  call solve_dg_affine_common_fixed_point(affine_rotations,affine_translations,1d-12,&
    has_affine_center,affine_center,affine_residual,ok,message)
  call require(ok.and..not.has_affine_center,&
    'valid screw affine action does not invent a common fixed point')
  call solve_dg_affine_common_fixed_point(affine_rotations(:,:,1:1),affine_translations(:,1:1),1d-12,&
    has_affine_center,affine_center,affine_residual,ok,message)
  call require(ok.and.has_affine_center.and.maxval(abs(affine_center))<1d-12,&
    'C1 has a deterministic canonical center without restricting its affine action')
  periodic_center_values=1d0
  periodic_center_phases(:,1)=exp(cmplx(0d0,2d0*acos(-1d0)*0.9d0,8))
  periodic_center_phases(:,2)=exp(cmplx(0d0,2d0*acos(-1d0)*0.1d0,8))
  call compute_dg_periodic_wannier_centers(comm,periodic_center_values,[0.5d0,0.5d0],&
    periodic_center_phases,periodic_centers,periodic_center_magnitudes,ok,message)
  call require(ok.and.maxval(min(periodic_centers,1d0-periodic_centers))<1d-12.and.&
    minval(periodic_center_magnitudes)>0.8d0,'periodic Wannier center crosses a cell face continuously')

  call build_dg_balanced_orbital_ownership(2*nproc+1,nproc,orbital_counts,&
    orbital_displacements,orbital_owners,ok,message)
  call require(ok.and.maxval(orbital_counts)-minval(orbital_counts)<=1,&
    'non-divisible orbital ownership is balanced')
  call require(sum(orbital_counts)==2*nproc+1.and.all(orbital_owners>=0).and.&
    all(orbital_owners<nproc),'balanced ownership covers every orbital exactly once')
  if(nproc>1)call require(maxval(orbital_counts)<2*nproc+1,&
    'no distributed orbital owner holds the complete full-system orbital set')
  allocate(transpose_local(2*nproc+1,2),transpose_local_ids(2))
  transpose_local_ids=[2_8*rank+1_8,2_8*rank+2_8]
  do point=1,2
    do i=1,2*nproc+1
      transpose_local(i,point)=cmplx(1000*i+transpose_local_ids(point),0d0,8)
    end do
  end do
  call transpose_dg_spatial_cores_to_orbital_owners(comm,transpose_local,transpose_local_ids,2,&
    owned_orbitals,transpose_global_ids,transpose_owned,ok,message)
  call require(ok.and.size(owned_orbitals)==orbital_counts(rank+1),&
    'spatial-to-orbital transpose returns only balanced owned orbitals')
  call require(size(transpose_global_ids)==2*nproc.and.&
    all(transpose_global_ids==[(int(i,8),i=1,2*nproc)]),&
    'orbital owner receives the complete unique-core physical grid')
  transpose_values_ok=.true.
  do point=1,size(transpose_global_ids)
    do i=1,size(owned_orbitals)
      transpose_values_ok=transpose_values_ok.and.abs(transpose_owned(i,point)-cmplx(&
        1000*owned_orbitals(i)+transpose_global_ids(point),0d0,8))<1d-14
    end do
  end do
  call require(transpose_values_ok,'MPI Alltoallv preserves orbital/core values')
  allocate(distributed_map(2,1),permuted_image(2*nproc+1,2))
  distributed_map(1,1)=int(2*modulo(rank+1,nproc)+1,8)
  distributed_map(2,1)=int(2*modulo(rank-1+nproc,nproc)+2,8)
  call exchange_dg_point_permuted_orbital_rows(comm,transpose_local,distributed_map(:,1),&
    permuted_image,ok,message)
  call require(ok.and.all(permuted_image(:,1)==cmplx(&
    [(1000*i+distributed_map(1,1),i=1,2*nproc+1)],0d0,8)).and.&
    all(permuted_image(:,2)==cmplx(&
    [(1000*i+distributed_map(2,1),i=1,2*nproc+1)],0d0,8)),&
    'sparse point exchange applies a cross-rank symmetry permutation exactly')
  distributed_map(1,1)=int(2*nproc+1,8)
  call exchange_dg_point_permuted_orbital_rows(comm,transpose_local,distributed_map(:,1),&
    permuted_image,ok,message)
  call require(.not.ok,'sparse point exchange collectively rejects a missing target')
  distributed_map(:,1)=int(2*rank+1,8)
  call exchange_dg_point_permuted_orbital_rows(comm,transpose_local,distributed_map(:,1),&
    permuted_image,ok,message)
  call require(.not.ok,'sparse point exchange collectively rejects duplicate targets')
  deallocate(distributed_map,permuted_image)
  allocate(center_owners(2*nproc+1))
  do i=1,size(center_owners);center_owners(i)=modulo(i-1,nproc);end do
  allocate(redistribution_buffer_ids(3))
  redistribution_buffer_ids=[transpose_local_ids,1_8]
  call redistribute_dg_owned_orbitals_to_center_fragments(comm,owned_orbitals,transpose_global_ids,&
    transpose_owned,center_owners,redistribution_buffer_ids,center_local_orbitals,&
    center_local_values,ok,message)
  call require(ok.and.all(center_local_orbitals==pack([(i,i=1,2*nproc+1)],center_owners==rank)),&
    'center-fragment redistribution receives exactly its centered orbitals')
  transpose_values_ok=size(center_local_values,2)==size(redistribution_buffer_ids)
  do point=1,size(redistribution_buffer_ids);do i=1,size(center_local_orbitals)
    transpose_values_ok=transpose_values_ok.and.abs(center_local_values(i,point)-cmplx(&
      1000*center_local_orbitals(i)+redistribution_buffer_ids(point),0d0,8))<1d-14
  end do;end do
  call require(transpose_values_ok,'center-fragment redistribution preserves core-buffer values')
  transpose_global_ids(size(transpose_global_ids))=transpose_global_ids(1)
  call redistribute_dg_owned_orbitals_to_center_fragments(comm,owned_orbitals,transpose_global_ids,&
    transpose_owned,center_owners,redistribution_buffer_ids,center_local_orbitals,&
    center_local_values,ok,message)
  call require(.not.ok,'center-fragment redistribution rejects duplicate global core IDs')
  transpose_global_ids=[(int(i,8),i=1,2*nproc)]
  call redistribute_dg_owned_orbitals_to_center_fragments(comm,owned_orbitals,transpose_global_ids,&
    transpose_owned,center_owners,redistribution_buffer_ids,center_local_orbitals,&
    center_local_values,ok,message)
  call require(ok,'valid core IDs remain accepted after duplicate-ID rejection')
  invalid_owned_orbitals=owned_orbitals
  if(rank==0.and.size(invalid_owned_orbitals)>0)invalid_owned_orbitals(1)=0
  call redistribute_dg_owned_orbitals_to_center_fragments(comm,invalid_owned_orbitals,&
    transpose_global_ids,transpose_owned,center_owners,redistribution_buffer_ids,&
    mismatched_owned_orbitals,mismatched_owned,ok,message)
  call require(.not.ok,'center-fragment redistribution rejects out-of-range owned orbital IDs')
  allocate(center_all_core_ids(8,nproc),center_fragments(nproc))
  center_all_core_ids=reshape([(int(i,8),i=1,8*nproc)],[8,nproc])
  center_fragments=[(i,i=1,nproc)]
  assignment_centers=0d0
  assignment_centers(1,1)=0.5d0/real(2*nproc,8)
  assignment_centers(1:2,2)=[0.5d0/real(2*nproc,8),0.25d0]
  assignment_centers(:,3)=[0.5d0/real(2*nproc,8),0.25d0,0.25d0]
  assignment_centers(:,4)=[real(2*nproc-1,8)/real(2*nproc,8),0.5d0,0.5d0]
  call assign_dg_periodic_centers_to_fragments([2*nproc,2,2],assignment_centers,&
    center_all_core_ids,center_fragments,1d-12,assigned_center_ids,assigned_center_owners,&
    assigned_center_fragments,ok,message)
  call require(ok.and.all(assigned_center_ids(1:3)==1_8).and.&
    all(assigned_center_owners(1:3)==0).and.assigned_center_ids(4)==int(8*nproc,8).and.&
    assigned_center_owners(4)==nproc-1,&
    'periodic face, edge, and corner ownership has deterministic lower-ID tie-breaking')
  if(nproc>1)then
    if(rank==nproc-1)then
      call transpose_dg_spatial_cores_to_orbital_owners(comm,transpose_local(1:2*nproc,:),&
        transpose_local_ids,2,mismatched_owned_orbitals,mismatched_global_ids,mismatched_owned,ok,message)
    else
      call transpose_dg_spatial_cores_to_orbital_owners(comm,transpose_local,transpose_local_ids,2,&
        mismatched_owned_orbitals,mismatched_global_ids,mismatched_owned,ok,message)
    end if
    call require(.not.ok,'spatial-to-orbital transpose rejects rank-dependent orbital counts')
  end if
  deallocate(transpose_local,transpose_local_ids,transpose_global_ids,transpose_owned,&
    owned_orbitals,orbital_counts,orbital_displacements,orbital_owners,center_owners,&
    center_local_orbitals,center_local_values,center_all_core_ids,center_fragments,&
    assigned_center_ids,assigned_center_owners,assigned_center_fragments,redistribution_buffer_ids)
  orbit_centers(:,1)=[0.1d0,0.2d0,0.3d0];orbit_centers(:,2)=[0.4d0,0.2d0,0.3d0]
  affine_translations=0d0;affine_translations(:,2)=[0.5d0,0.4d0,0d0]
  call verify_dg_wannier_center_affine_orbits(orbit_centers,affine_rotations,&
    affine_translations,1d-12,ok,message)
  call require(ok,'Wannier centers form a closed full affine orbit independently of owner fragment')
  orbit_centers(1,2)=0.45d0
  call verify_dg_wannier_center_affine_orbits(orbit_centers,affine_rotations,&
    affine_translations,1d-12,ok,message)
  call require(.not.ok,'broken full affine Wannier center orbit is rejected')
  calibrated_map(:,1)=int(rank*4,8)+[1_8,2_8,3_8,4_8]
  calibrated_map(:,2)=int(rank*4,8)+[2_8,1_8,4_8,3_8]
  calibrated_boundary=[.true.,.true.,.false.,.false.]
  calibrated_basis=1d0;calibrated_basis(1,1)=1.1d0
  call measure_dg_rank_fixed_symmetry_residuals(comm,calibrated_basis,[1d0,1d0,1d0,1d0],&
    calibrated_map,calibrated_boundary,calibrated_representation,calibrated_total,&
    calibrated_boundary_residual,calibrated_interior_residual,ok,message,&
    workspace_peak_bytes=symmetry_workspace_peak)
  call require(symmetry_workspace_peak>0_8,&
    'streamed symmetry residual measurement publishes a workspace receipt')
  call require(ok.and.calibrated_boundary_residual(2)>10d0*calibrated_interior_residual(2),&
    'rank-fixed residual identifies boundary stitching error')
  calibrated_allowance=1.01d0*calibrated_boundary_residual(2)
  calibrated_strict_tolerance=max(1d-12,2d0*calibrated_interior_residual(2))
  call accept_dg_boundary_calibrated_symmetry(calibrated_boundary_residual,&
    calibrated_interior_residual,calibrated_allowance,calibrated_strict_tolerance,ok,message)
  call require(ok,'measured stitching baseline accepts boundary-localized residual')
  calibrated_basis=1d0;calibrated_basis(1,3)=1.1d0
  call measure_dg_rank_fixed_symmetry_residuals(comm,calibrated_basis,[1d0,1d0,1d0,1d0],&
    calibrated_map,calibrated_boundary,calibrated_representation,calibrated_total,&
    calibrated_boundary_residual,calibrated_interior_residual,ok,message)
  call require(ok.and.calibrated_interior_residual(2)>10d0*calibrated_boundary_residual(2),&
    'rank-fixed residual rejects equivalent interior breaking')
  call accept_dg_boundary_calibrated_symmetry(calibrated_boundary_residual,&
    calibrated_interior_residual,calibrated_allowance,calibrated_strict_tolerance,ok,message)
  call require(.not.ok,'boundary allowance cannot excuse interior symmetry breaking')
  lcfo_core_ids(1)=int(rank+1,8)
  lcfo_buffer_ids=[int(rank+1,8),int(modulo(rank+1,nproc)+1,8)]
  lcfo_buffer_contribution(:,1)=[cmplx(rank+1d0,0d0,8),cmplx(10d0*(rank+1),0d0,8)]
  lcfo_buffer_contribution(:,2)=[cmplx(100d0*(rank+1),0d0,8),cmplx(1000d0*(rank+1),0d0,8)]
  call accumulate_dg_lcfo_buffer_contributions_to_core(comm,lcfo_buffer_ids,&
    lcfo_buffer_contribution,lcfo_core_ids,lcfo_core_value,ok,message)
  call require(ok,trim(message))
  if(nproc==1)then
    call require(maxval(abs(lcfo_core_value(:,1)-[cmplx(101d0,0d0,8),cmplx(1010d0,0d0,8)]))<1d-12,&
      'LCFO contributions sharing one physical point are accumulated once per basis fragment')
  else
    i=modulo(rank-1,nproc)
    call require(maxval(abs(lcfo_core_value(:,1)-[cmplx(rank+1d0+100d0*(i+1),0d0,8),&
      cmplx(10d0*(rank+1)+1000d0*(i+1),0d0,8)]))<1d-12,&
      'LCFO buffer tails from every covering fragment accumulate on the unique core owner')
  endif
  closure_metric=(0d0,0d0);closure_localizer=(0d0,0d0);closure_group=(0d0,0d0)
  closure_occupied=(0d0,0d0)
  do i=1,6
    closure_metric(i,i)=1d0;closure_group(i,i,1)=1d0
  enddo
  closure_group(2,1,2)=1d0;closure_group(1,2,2)=1d0
  closure_group(4,3,2)=1d0;closure_group(3,4,2)=1d0
  closure_group(6,5,2)=1d0;closure_group(5,6,2)=1d0
  closure_localizer(1,1)=0d0;closure_localizer(2,2)=0d0
  closure_localizer(3,3)=1d0;closure_localizer(4,4)=1d0
  closure_localizer(5,5)=2d0;closure_localizer(6,6)=2d0
  closure_occupied(1,1)=1d0;closure_occupied(2,2)=1d0
  closure_product=reshape([1,2,2,1],[2,2])
  call find_dg_group_identity(reshape([2,1,1,2],[2,2]),identity_operation,ok,message)
  call require(ok.and.identity_operation==2,'group identity is derived from the product table')
  call select_dg_fixed_rank_symmetry_closed_subspace(closure_metric,closure_occupied,&
    closure_localizer,closure_group,closure_product,4,1d-12,closure_transform,&
    occupied_inclusion,subspace_leakage,ok,message)
  call require(ok.and.all(shape(closure_transform)==[6,4]),trim(message))
  projector=matmul(closure_transform,conjg(transpose(closure_transform)))
  call require(maxval(abs(projector(1:4,1:4)-closure_metric(1:4,1:4)))<1d-12.and.&
    maxval(abs(projector(5:6,:)))<1d-12,'fixed-rank selector keeps complete symmetry blocks')
  call require(occupied_inclusion<1d-12.and.subspace_leakage<1d-12,&
    'fixed-rank selector contains occupied space and closes under the group')
  scaled_closure_metric=2d0*closure_metric
  call select_dg_fixed_rank_symmetry_closed_subspace(scaled_closure_metric,closure_occupied,&
    closure_localizer,closure_group,closure_product,4,1d-12,closure_transform,&
    occupied_inclusion,subspace_leakage,ok,message)
  call require(ok.and.occupied_inclusion<1d-12.and.subspace_leakage<1d-12,&
    'fixed-rank selector metric-orthonormalizes the occupied coefficient block')
  call select_dg_fixed_rank_symmetry_closed_subspace(closure_metric,closure_occupied,&
    closure_localizer,closure_group,closure_product,3,1d-12,closure_transform,&
    occupied_inclusion,subspace_leakage,ok,message)
  call require(.not.ok.and.trim(message)=='target rank cuts a symmetry-degenerate block',&
    'fixed-rank selector rejects an incomplete symmetry block')
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
  mixed_wannier=reshape([1d0,1d0,1d0,-1d0],[2,2])/sqrt(2d0)
  mixed_map(:,1)=[2_8,1_8]
  call verify_dg_fragment_subspace_density_covariance(mixed_wannier,mixed_map,1d-12,ok,message)
  call require(ok,'unitary-mixed Wannier subspace density is symmetry covariant')
  mixed_wannier(1,1)=2d0*mixed_wannier(1,1)
  call verify_dg_fragment_subspace_density_covariance(mixed_wannier,mixed_map,1d-12,ok,message)
  call require(.not.ok,'broken Wannier subspace density covariance rejected')
  fractional_core_candidates=0d0
  fractional_core_candidates(1,1)=sqrt(0.6d0)
  fractional_core_candidates(2,1)=sqrt(0.3d0)
  fractional_core_candidates(1,2)=sqrt(0.4d0)
  fractional_core_candidates(2,2)=sqrt(0.7d0)
  call build_dg_core_owned_occupied_subspace(fractional_core_candidates,[.true.,.false.],&
    [1d0,1d0],[2d0,2d0],2d0,core_occupied_coefficients,fractional_core_electrons,ok,message)
  call require(ok.and.size(core_occupied_coefficients,2)==1,&
    'ionic electron ownership fixes rank despite fractional instantaneous core charge')
  call require(abs(fractional_core_electrons-1.8d0)<1d-12,&
    'fractional instantaneous core charge remains diagnostic evidence')
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
    allocate(affine_local_ids(2),affine_all_ids(2,2))
    affine_local_ids=[int(2*rank+1,8),int(2*rank+2,8)]
    affine_all_ids=reshape([1_8,2_8,3_8,4_8],[2,2])
    affine_rotation=0
    affine_rotation(1,1)=-1;affine_rotation(2,2)=1;affine_rotation(3,3)=1
    affine_translation=[0.5d0,0d0,0d0]
    call build_dg_pointwise_affine_owner_map([4,1,1],affine_local_ids,affine_all_ids,&
      affine_rotation,affine_translation,1d-12,affine_target_ids,affine_target_owner,&
      affine_target_local,affine_wrap,ok,message)
    call require(ok,trim(message))
    call build_dg_pointwise_affine_owner_map([4,1,1],affine_target_ids,affine_all_ids,&
      affine_rotation,affine_translation,1d-12,affine_second_ids,affine_second_owner,&
      affine_second_local,affine_second_wrap,ok,message)
    call require(ok.and.all(affine_second_ids==affine_local_ids),&
      'pointwise affine inversion closes on global physical IDs')
    do i=1,2
      call require(all(matmul(affine_rotation,affine_wrap(:,i))+affine_second_wrap(:,i)==0),&
        'lattice-wrap cocycle makes inversion phase square to identity')
    enddo
    affine_rotation=0
    call build_dg_pointwise_affine_owner_map([4,1,1],affine_local_ids,affine_all_ids,&
      affine_rotation,affine_translation,1d-12,affine_second_ids,affine_second_owner,&
      affine_second_local,affine_second_wrap,ok,message)
    call require(.not.ok.and.trim(message)=='affine rotation must be unimodular',&
      'noninvertible affine rotation is rejected before owner mapping')
    if(rank==0)then
      call require(all(affine_target_ids==[3_8,2_8]),&
        'external-center inversion splits one source core across owners')
      call require(all(affine_target_owner==[1,0]).and.all(affine_target_local==[1,2]),&
        'split inversion resolves each target owner and local index')
      call require(all(affine_wrap==0),'rank-zero inversion images require no lattice wrap')
    else
      call require(all(affine_target_ids==[1_8,4_8]),&
        'pointwise inversion preserves global owner IDs independently of fragments')
      call require(all(affine_target_owner==[0,1]).and.all(affine_target_local==[1,2]),&
        'pointwise inversion owner resolution is rank independent')
      call require(affine_wrap(1,1)==0.and.affine_wrap(1,2)==-1,&
        'pointwise inversion records the negative periodic lattice wrap')
    endif
    deallocate(affine_local_ids,affine_all_ids,affine_target_ids,affine_target_owner,&
      affine_target_local,affine_wrap)
    if(allocated(affine_second_ids))deallocate(affine_second_ids)
    if(allocated(affine_second_owner))deallocate(affine_second_owner)
    if(allocated(affine_second_local))deallocate(affine_second_local)
    if(allocated(affine_second_wrap))deallocate(affine_second_wrap)

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

    allocate(distributed_basis(2,2),distributed_weight(2),distributed_map(2,3))
    distributed_basis=(0d0,0d0)
    distributed_basis(rank+1,:)=[(1d0,0d0),(2d0,0d0)]
    distributed_weight=1d0
    distributed_map(:,1)=[int(2*rank+1,8),int(2*rank+2,8)]
    distributed_map(:,2)=[int(2*(1-rank)+1,8),int(2*(1-rank)+2,8)]
    distributed_map(:,3)=[int(2*rank+1,8),int(2*(1-rank)+2,8)]
    call assemble_dg_distributed_basis_symmetry_overlap(comm,distributed_basis,distributed_weight,&
      distributed_map,distributed_basis_overlap,ok,message)
    call require(ok,trim(message))
    call require(maxval(abs(distributed_basis_overlap(:,:,1)-&
      reshape([(5d0,0d0),(0d0,0d0),(0d0,0d0),(5d0,0d0)],[2,2])))<1d-12,&
      'distributed full-basis identity overlap')
    call require(maxval(abs(distributed_basis_overlap(:,:,2)-&
      reshape([(0d0,0d0),(5d0,0d0),(5d0,0d0),(0d0,0d0)],[2,2])))<1d-12,&
      'distributed full-basis fragment-swap overlap')
    call require(maxval(abs(distributed_basis_overlap(:,:,3)-&
      reshape([(1d0,0d0),(4d0,0d0),(4d0,0d0),(1d0,0d0)],[2,2])))<1d-12,&
      'distributed full-basis operation may split one core across owners')
    call assemble_dg_distributed_basis_symmetry_overlap_rows(comm,distributed_basis,distributed_weight,&
      distributed_map,distributed_overlap_row_ids,distributed_basis_overlap_rows,&
      row_overlap_workspace_peak,ok,message)
    call require(ok.and.row_overlap_workspace_peak>0_8,trim(message))
    call require(maxval(abs(distributed_basis_overlap_rows-&
      distributed_basis_overlap(int(distributed_overlap_row_ids),:,:)))<1d-12,&
      'row-owned symmetry overlaps match the dense reference')
    call require(size(distributed_basis_overlap_rows,1)==1,&
      'two-rank symmetry overlap owns only one global basis row per rank')

    allocate(orbit_seed(1,2));orbit_seed=(0d0,0d0)
    if(rank==0)orbit_seed(1,1)=1d0
    call build_dg_distributed_symmetry_closed_basis(comm,orbit_seed,distributed_weight,&
      distributed_map(:,1:2),closure_product,1,2,1d-12,orbit_basis,orbit_rank,ok,message)
    call require(ok.and.orbit_rank==2.and.all(shape(orbit_basis)==[2,2]),trim(message))
    orbit_gram=matmul(orbit_basis,conjg(transpose(orbit_basis)))
    call MPI_Allreduce(MPI_IN_PLACE,orbit_gram,4,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call require(maxval(abs(orbit_gram-reshape([(1d0,0d0),(0d0,0d0),&
      (0d0,0d0),(1d0,0d0)],[2,2])))<1d-12,&
      'streamed symmetry orbit is globally orthonormal')
    call build_dg_distributed_symmetry_closed_basis(comm,orbit_seed,distributed_weight,&
      distributed_map(:,1:2),closure_product,1,1,1d-12,orbit_basis,orbit_rank,ok,message)
    call require(.not.ok.and.trim(message)=='required symmetry orbit exceeds target rank',&
      'required occupied orbit cannot be truncated to the target rank')
    call build_dg_distributed_symmetry_closed_basis(comm,orbit_seed,distributed_weight,&
      distributed_map(:,1:2),closure_product,1,2,1d-12,orbit_basis,orbit_rank,ok,message,&
      minimum_rank=1,required_retained_rank=required_orbit_rank)
    call require(ok.and.orbit_rank==2.and.required_orbit_rank==2,&
      'candidate builder may cross a minimum rank only after completing the symmetry orbit')
    allocate(required_orbit_seed(2,2));required_orbit_seed=(0d0,0d0)
    if(rank==0)then
      required_orbit_seed(1,1)=1d0
      required_orbit_seed(2,2)=1d0
    end if
    call build_dg_distributed_symmetry_closed_basis(comm,required_orbit_seed,distributed_weight,&
      distributed_map(:,1:2),closure_product,2,4,1d-12,orbit_basis,orbit_rank,ok,message,&
      minimum_rank=1,required_retained_rank=required_orbit_rank)
    call require(ok.and.orbit_rank==4.and.required_orbit_rank==4,&
      'minimum-rank candidate construction still processes every required seed')
    deallocate(required_orbit_seed)
    allocate(invalid_orbit_map,source=distributed_map(:,1:2))
    if(rank==0)invalid_orbit_map(:,2)=[3_8,3_8]
    call build_dg_distributed_symmetry_closed_basis(comm,orbit_seed,distributed_weight,&
      invalid_orbit_map,closure_product,1,2,1d-12,orbit_basis,orbit_rank,ok,message)
    call require(.not.ok.and.trim(message)=='point maps are not a closed permutation group',&
      'streamed symmetry closure rejects a nonbijective global point action')
    deallocate(invalid_orbit_map)
    deallocate(orbit_seed,orbit_gram)
    if(allocated(orbit_basis))deallocate(orbit_basis)
    deallocate(distributed_basis,distributed_weight,distributed_map,distributed_basis_overlap)
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
  call require(allocated(result%center_box_point_ids).and.&
    size(result%center_box_point_ids)==result%target_rank,&
    'symmetry-free local construction still publishes Wannier centers')
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
