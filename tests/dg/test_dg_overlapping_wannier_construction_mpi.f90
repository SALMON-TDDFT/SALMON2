#include "config.h"
program test_dg_overlapping_wannier_construction_mpi
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use mpi
  use dg_overlapping_wannier_types,only:s_dg_ow_distributed_layout,&
    initialize_dg_ow_distributed_layout,reserve_dg_ow_workspace,&
    release_dg_ow_workspace,release_dg_ow_distributed_layout
  use dg_overlapping_wannier_construction,only:s_dg_overlapping_wannier_construction,&
    s_dg_translation_orbit_accumulator,&
    s_dg_prepared_translation_action,&
    construct_dg_overlapping_wannier_basis,release_dg_overlapping_wannier_construction,&
    verify_dg_overlapping_wannier_periodic_closure,assemble_dg_distributed_candidate_symmetry,&
    assemble_dg_distributed_basis_symmetry_overlap,&
    assemble_dg_distributed_basis_symmetry_overlap_rows,&
    validate_dg_row_owned_group_representation,&
    validate_dg_streamed_affine_representation,&
    select_dg_group_generators,&
    gather_dg_single_symmetry_representation,&
    build_dg_pointwise_affine_owner_map,&
    find_dg_group_identity,&
    select_dg_fixed_rank_symmetry_closed_subspace,&
    build_dg_distributed_symmetry_closed_basis,&
    build_dg_group_averaged_occupied_candidates_reference,&
    orthonormalize_dg_distributed_seed_space,&
    accumulate_dg_lcfo_buffer_contributions_to_core,&
    measure_dg_rank_fixed_symmetry_residuals,&
    accept_dg_boundary_calibrated_symmetry,&
    solve_dg_affine_common_fixed_point,&
    compute_dg_periodic_wannier_centers,&
    verify_dg_wannier_center_affine_orbits,&
    diagnose_dg_point_center_gauge,&
    build_dg_finite_abelian_character_table,&
    inverse_dg_translation_character_orbits,&
    accumulate_dg_translation_character_orbit_sector,&
    accumulate_dg_translation_character_orbit_sector_values,&
    apply_dg_row_owned_orbital_transform_streamed,&
    materialize_dg_row_owned_sector_on_spatial_grid,&
    build_dg_translation_character_intertwining_phase,&
    prepare_dg_translation_character_action,&
    build_dg_translation_character_intertwining_phase_prepared,&
    release_dg_prepared_translation_action,&
    validate_dg_factored_point_cogroup_gauge,&
    align_dg_fragment_wannier_gauge,replicate_dg_fragment_wannier_representative,&
    verify_dg_fragment_wannier_streaming_closure,verify_dg_fragment_center_orbit,&
    verify_dg_uniform_fragment_target_rank,assign_dg_overlapping_wannier_occupations,&
    build_dg_balanced_orbital_ownership,transpose_dg_spatial_cores_to_orbital_owners,&
    exchange_dg_point_permuted_orbital_rows,&
    redistribute_dg_buffer_orbitals_to_center_fragments,&
    redistribute_dg_owned_orbitals_to_center_fragments,&
    assign_dg_periodic_centers_to_fragments,&
    verify_dg_fragment_subspace_density_covariance,build_dg_core_owned_occupied_subspace,&
    build_dg_smooth_partition_of_unity,compose_dg_buffered_orbital_tile_to_physical_grid,&
    build_dg_equal_count_spectral_windows,build_dg_spectral_density_descriptors
  implicit none
  type(s_dg_translation_orbit_accumulator)::inverse_accumulator
  type(s_dg_prepared_translation_action)::prepared_translation_action
  integer::comm,rank,nproc,ierr,i,j,b,p,point,nlocal,nclosure,index,ncore,fragment_id
  integer(8),allocatable::ids(:),box_ids(:),symmetry_map(:,:),broken_symmetry_map(:,:)
  integer,allocatable::fragment(:)
  real(8),allocatable::weight(:),coordinate(:)
  logical,allocatable::boundary(:),core_mask(:)
  complex(8),allocatable::candidate(:,:),gradient(:,:,:),occupied(:,:),base_occupied(:,:),&
    rotated(:,:),rotated_gradient(:,:,:),periodic_phase(:,:)
  complex(8),allocatable::direct_small(:,:),direct_small_gradient(:,:,:),&
    direct_small_occupied(:,:),direct_projector(:,:)
  complex(8),allocatable::stream_values(:,:),stream_gradients(:,:,:)
  complex(8),allocatable::spatial_basis(:,:),sector_coefficients(:,:),materialized_sector(:,:)
  integer(8),allocatable::sector_coefficient_ids(:)
  complex(8),allocatable::mismatch_values(:,:),mismatch_gradients(:,:,:)
  complex(8)::gauge(4,4)
  complex(8),allocatable::reference_projector(:,:),projector(:,:),seed_projector(:,:)
  complex(8),allocatable::distributed_candidate(:,:),distributed_overlap(:,:,:),&
    distributed_basis(:,:),distributed_basis_overlap(:,:,:),orbit_seed(:,:),required_orbit_seed(:,:),&
    orbit_basis(:,:),orthonormal_seed_basis(:,:),averaged_candidates(:,:),&
    orbit_gram(:,:)
  complex(8),allocatable::distributed_basis_overlap_rows(:,:,:)
  complex(8),allocatable::single_symmetry_representation(:,:)
  complex(8)::lcfo_buffer_contribution(2,2),lcfo_core_value(2,1)
  integer(8)::lcfo_buffer_ids(2),lcfo_core_ids(1)
  integer::orbit_rank,required_orbit_rank,averaged_rank,identity_operation
  integer(8),allocatable::distributed_map(:,:),invalid_orbit_map(:,:)
  integer(8),allocatable::affine_local_ids(:),affine_all_ids(:,:),affine_target_ids(:),&
    affine_second_ids(:)
  integer,allocatable::affine_target_owner(:),affine_target_local(:),affine_wrap(:,:),&
    affine_second_owner(:),affine_second_local(:),affine_second_wrap(:,:)
  integer::affine_rotation(3,3)
  real(8)::affine_translation(3)
  real(8),allocatable::distributed_weight(:)
  real(8),allocatable::averaged_spectrum(:)
  real(8),allocatable::seed_values(:,:)
  real(8),allocatable::occupied_seed_values(:,:)
  real(8),allocatable::raw_seed_values(:,:)
  real(8),allocatable::spectral_window_weights(:,:)
  real(8),allocatable::spectral_occupied_density(:),spectral_unoccupied_density(:,:),&
    spectral_shared_density(:,:),spectral_reference_descriptors(:,:)
  complex(8),allocatable::spectral_state_values(:,:)
  integer(8),allocatable::spectral_row_ids(:)
  real(8)::spectral_eigenvalues(10),spectral_occupations(10)
  real(8)::occupations(3)
  type(s_dg_overlapping_wannier_construction)::result
  type(s_dg_ow_distributed_layout)::distributed_layout
  integer,allocatable::reference_owner(:),reference_center_fragment(:)
  integer(8),allocatable::reference_center_box_ids(:)
  integer(8)::reference_fingerprint
  integer(8)::symmetry_workspace_peak
  integer(8)::row_overlap_workspace_peak
  real(8)::row_identity_defect,row_unitarity_defect,row_closure_defect
  real(8)::dense_identity_defect,dense_unitarity_defect,dense_closure_defect
  integer(8),allocatable::distributed_overlap_row_ids(:)
  integer(8)::closure_fingerprint,rounded_closure_fingerprint
  integer(8)::spectral_window_fingerprint,spectral_window_workspace
  integer(8)::spectral_density_fingerprint,spectral_density_workspace
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
  integer::cyclic_product(4,4)
  integer::character_canonical_operations(4),character_inverses(4),character_conjugates(4),&
    character_generator_count
  integer,allocatable::character_generators(:),character_words(:,:)
  integer::character_product(4,4)
  integer::character_permutation(4),character_inverse_permutation(4),character_permuted_product(4,4)
  real(8)::character_translations(3,4),character_permuted_translations(3,4)
  complex(8)::character_table(4,4),character_gram(4,4)
  integer(8)::character_fingerprint,character_reference_fingerprint
  integer(8)::spatial_sector_fingerprint
  integer::z6_product(6,6),z6_operations(6),z6_inverses(6),z6_conjugates(6),z6_exponent(6),&
    z6_exponent_to_operation(0:5),z6_generator_count
  integer,allocatable::z6_generators(:),z6_words(:,:)
  real(8)::z6_translations(3,6)
  complex(8)::z6_characters(6,6)
  integer::z2_product(2,2),z2_operations(2),z2_inverses(2),z2_conjugates(2),z2_generator_count
  integer,allocatable::z2_generators(:),z2_words(:,:)
  real(8)::z2_translations(3,2)
  complex(8)::z2_characters(2,2)
  integer,allocatable::group_generators(:)
  real(8)::subspace_leakage,occupied_inclusion
  real(8)::averaged_trace,averaged_closure
  integer(8)::averaged_workspace_peak
  integer(8)::partition_ids(2)
  real(8)::raw_partition_weight(2),raw_partition_gradient(3,2),partition_weight(2),&
    partition_gradient(3,2),partition_sum_defect,partition_gradient_defect
  integer(8),allocatable::variable_partition_ids(:)
  real(8),allocatable::variable_raw_partition_weight(:),variable_raw_partition_gradient(:,:),&
    variable_partition_weight(:),variable_partition_gradient(:,:)
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
  complex(8),allocatable::composed_buffer_values(:,:),composed_owned_values(:,:)
  integer(8),allocatable::transpose_local_ids(:),transpose_global_ids(:)
  integer(8),allocatable::mismatched_global_ids(:)
  integer(8),allocatable::redistribution_buffer_ids(:)
  integer(8),allocatable::composed_owned_ids(:)
  integer(8),allocatable::center_all_core_ids(:,:),assigned_center_ids(:)
  integer,allocatable::orbital_counts(:),orbital_displacements(:),orbital_owners(:),owned_orbitals(:)
  integer,allocatable::invalid_owned_orbitals(:)
  integer,allocatable::mismatched_owned_orbitals(:)
  integer,allocatable::center_owners(:),center_local_orbitals(:)
  integer,allocatable::center_fragments(:),assigned_center_owners(:),assigned_center_fragments(:)
  real(8)::assignment_centers(3,4)
  real(8)::orbit_centers(3,2),orbit_center_magnitudes(3,2)
  integer::failed_center_operation
  complex(8)::center_gauge_basis(2,2)
  real(8)::center_gauge_weights(2),center_gauge_centers(3,2),center_gauge_tau(3)
  integer::center_gauge_rotation(3,3)
  integer(8)::center_gauge_map(2),center_gauge_workspace
  real(8)::center_gauge_monomial,center_gauge_leakage,center_gauge_unitarity
  complex(8)::fractional_core_candidates(2,2)
  complex(8),allocatable::core_occupied_coefficients(:,:)
  integer(8)::mixed_map(2,1)
  integer(8)::composition_fingerprint,composition_workspace_peak
  real(8)::fractional_core_electrons
  complex(8),allocatable::inverse_sector_values(:,:,:),inverse_sector_gradients(:,:,:,:),&
    inverse_orbit_values(:,:,:),inverse_orbit_gradients(:,:,:,:),inverse_streamed_values(:,:,:),&
    inverse_streamed_gradients(:,:,:,:),inverse_values_only(:,:),inverse_rotated_values(:,:,:),&
    inverse_rotated_gradients(:,:,:,:),inverse_translated_values(:,:,:),inverse_translated_gradients(:,:,:,:),&
    streamed_transform_rows(:,:),streamed_input_values(:,:),streamed_input_gradients(:,:,:),&
    streamed_output_values(:,:),streamed_output_gradients(:,:,:),intertwining_phase(:),inverse_saved_sector(:,:)
  complex(8),allocatable::prepared_intertwining_phase(:)
  real(8),allocatable::inverse_density(:),inverse_rotated_density(:)
  complex(8)::inverse_internal_gauge(2,2)
  real(8)::inverse_density_defect,inverse_orthogonality_defect,inverse_gamma_defect
  integer(8)::inverse_fingerprint,inverse_workspace,intertwining_payload_fingerprint
  integer(8)::prepared_intertwining_fingerprint,prepared_intertwining_payload_fingerprint,&
    prepared_intertwining_workspace
  integer(8)::inverse_reference_fingerprint
  integer(8),allocatable::inverse_row_ids(:)
  integer(8)::intertwining_maps(12,4),intertwining_generator_maps(12,2)
  integer(8),allocatable::factored_point_maps(:,:),factored_translation_maps(:,:)
  integer::factored_point_product(2,2),factored_cocycle(2,2),factored_global_index
  integer(8),allocatable::semidirect_point_maps(:,:),semidirect_translation_maps(:,:)
  integer,allocatable::semidirect_product(:,:),semidirect_cocycle(:,:)
  integer::semidirect_units(4),semidirect_offsets(4),semidirect_left,semidirect_right,&
    semidirect_product_unit,semidirect_product_index,semidirect_shift
  integer::factored_generator_count,factored_checked_pair_count,factored_prepared_operation_count
  integer::factored_receipt_min,factored_receipt_max
  real(8),allocatable::factored_weights(:)
  real(8)::factored_identity_defect,factored_unitarity_defect,factored_closure_defect
  real(8)::gauge_weights(2)
  logical::ok,transpose_values_ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  character_translations=0d0;character_translations(1,:)=[0d0,0.5d0,0.25d0,0.75d0]
  character_product=reshape([1,2,3,4,2,1,4,3,3,4,2,1,4,3,1,2],[4,4])
  call build_dg_finite_abelian_character_table(character_translations,character_product,1,1d-12,&
    character_canonical_operations,character_inverses,character_generator_count,character_generators,&
    character_words,character_table,character_conjugates,character_fingerprint,ok,message)
  call require(ok,trim(message))
  character_gram=matmul(character_table,conjg(transpose(character_table)))
  call require(character_generator_count==1.and.size(character_generators)==1.and.&
    all(character_canonical_operations==[1,3,2,4]).and.all(character_inverses==[1,4,3,2]).and.&
    all(character_conjugates==[1,4,3,2]).and.all(character_words(:,1)==[0,1,2,3]).and.&
    maxval(abs(character_gram-&
    4d0*reshape([(1d0,0d0),(0d0,0d0),(0d0,0d0),(0d0,0d0),&
      (0d0,0d0),(1d0,0d0),(0d0,0d0),(0d0,0d0),&
      (0d0,0d0),(0d0,0d0),(1d0,0d0),(0d0,0d0),&
      (0d0,0d0),(0d0,0d0),(0d0,0d0),(1d0,0d0)],[4,4])))<1d-12,&
    'canonical Z4 generators, words, and characters are exact')
  character_reference_fingerprint=character_fingerprint
  character_permutation=[3,1,4,2]
  do i=1,4;character_inverse_permutation(character_permutation(i))=i;enddo
  character_permuted_translations=character_translations(:,character_permutation)
  do i=1,4;do j=1,4
    character_permuted_product(i,j)=character_inverse_permutation(&
      character_product(character_permutation(i),character_permutation(j)))
  enddo;enddo
  call build_dg_finite_abelian_character_table(character_permuted_translations,&
    character_permuted_product,character_inverse_permutation(1),1d-12,&
    character_canonical_operations,character_inverses,character_generator_count,character_generators,&
    character_words,character_table,character_conjugates,character_fingerprint,ok,message)
  call require(ok.and.character_fingerprint==character_reference_fingerprint,&
    'character table fingerprint is invariant under operation numbering')
  character_translations=0d0
  character_translations(1,:)=[0d0,0d0,0.5d0,0.5d0]
  character_translations(2,:)=[0d0,0.5d0,0d0,0.5d0]
  character_product=reshape([1,2,3,4,2,1,4,3,3,4,1,2,4,3,2,1],[4,4])
  call build_dg_finite_abelian_character_table(character_translations,character_product,1,1d-12,&
    character_canonical_operations,character_inverses,character_generator_count,character_generators,&
    character_words,character_table,character_conjugates,character_fingerprint,ok,message)
  character_gram=matmul(character_table,conjg(transpose(character_table)))
  call require(ok.and.character_generator_count==2.and.&
    maxval(abs(character_gram-4d0*reshape([(1d0,0d0),(0d0,0d0),(0d0,0d0),(0d0,0d0),&
      (0d0,0d0),(1d0,0d0),(0d0,0d0),(0d0,0d0),(0d0,0d0),(0d0,0d0),&
      (1d0,0d0),(0d0,0d0),(0d0,0d0),(0d0,0d0),(0d0,0d0),(1d0,0d0)],[4,4])))<1d-12,&
    'Z2xZ2 character table has two generators and four orthogonal characters')
  character_product(2,3)=2
  call build_dg_finite_abelian_character_table(character_translations,character_product,1,1d-12,&
    character_canonical_operations,character_inverses,character_generator_count,character_generators,&
    character_words,character_table,character_conjugates,character_fingerprint,ok,message)
  call require(.not.ok,'corrupt translation product table is rejected')
  character_product=reshape([1,2,3,4,2,1,4,3,3,4,1,2,4,3,2,1],[4,4])
  character_translations(:,4)=character_translations(:,1)
  character_translations(1,4)=-0.5d-12
  call build_dg_finite_abelian_character_table(character_translations,character_product,1,1d-12,&
    character_canonical_operations,character_inverses,character_generator_count,character_generators,&
    character_words,character_table,character_conjugates,character_fingerprint,ok,message)
  call require(.not.ok.and.trim(message)=='finite translation catalog is nonfaithful',&
    'periodically coincident translations are rejected across the zero boundary')
  character_translations=0d0
  character_translations(1,:)=[0d0,0d0,0.5d0,0.5d0]
  character_translations(2,:)=[0d0,0.5d0,0d0,0.5d0]
  call build_dg_finite_abelian_character_table(character_translations,character_product,2,1d-12,&
    character_canonical_operations,character_inverses,character_generator_count,character_generators,&
    character_words,character_table,character_conjugates,character_fingerprint,ok,message)
  call require(.not.ok.and.trim(message)=='designated identity is not the zero translation',&
    'incorrect designated identity is rejected explicitly')
  z6_exponent=[0,4,2,3,1,5]
  do i=1,6;z6_exponent_to_operation(z6_exponent(i))=i;enddo
  z6_translations=0d0
  do i=1,6
    z6_translations(1,i)=modulo(0.5d0*real(z6_exponent(i),8),1d0)
    z6_translations(2,i)=modulo(real(z6_exponent(i),8)/3d0,1d0)
    do j=1,6
      z6_product(i,j)=z6_exponent_to_operation(modulo(z6_exponent(i)+z6_exponent(j),6))
    enddo
  enddo
  call build_dg_finite_abelian_character_table(z6_translations,z6_product,1,1d-12,z6_operations,&
    z6_inverses,z6_generator_count,z6_generators,z6_words,z6_characters,z6_conjugates,&
    character_fingerprint,ok,message)
  call require(ok.and.z6_generator_count==1,&
    'Z6 uses the minimum one generator instead of irredundant order-two/order-three generators')
  z2_translations=0d0;z2_translations(1,2)=0.5d0
  z2_product=reshape([1,2,2,1],[2,2])
  call build_dg_finite_abelian_character_table(z2_translations,z2_product,1,1d-12,z2_operations,&
    z2_inverses,z2_generator_count,z2_generators,z2_words,z2_characters,z2_conjugates,&
    character_fingerprint,ok,message)
  call require(ok.and.z2_generator_count==1.and.all(z2_inverses==[1,2]).and.&
    all(z2_conjugates==[1,2]),'Z2 characters are self-conjugate with exact inverses')
  call build_dg_finite_abelian_character_table(z2_translations,z2_product,1,&
    8d0*acos(-1d0)/real(huge(0_8),8),z2_operations,z2_inverses,z2_generator_count,&
    z2_generators,z2_words,z2_characters,z2_conjugates,character_fingerprint,ok,message)
  call require(.not.ok.and.trim(message)==&
    'finite translation tolerance is too small for canonical int64 keys',&
    'phase fingerprint int64 overflow is rejected before quantization')
  character_translations=ieee_value(0d0,ieee_quiet_nan)
  call build_dg_finite_abelian_character_table(character_translations,character_product,1,1d-12,&
    character_canonical_operations,character_inverses,character_generator_count,character_generators,&
    character_words,character_table,character_conjugates,character_fingerprint,ok,message)
  call require(.not.ok.and.trim(message)=='invalid finite-abelian character-table contract',&
    'nonfinite translation catalog is rejected')
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
  call redistribute_dg_buffer_orbitals_to_center_fragments(comm,transpose_local,[1,2],&
    transpose_local_ids,center_owners,redistribution_buffer_ids,mismatched_owned_orbitals,&
    mismatched_owned,ok,message)
  call require(ok.and.all(mismatched_owned_orbitals==center_local_orbitals).and.&
    all(abs(mismatched_owned-center_local_values)<1d-14),&
    'direct buffer redistribution matches the established two-stage result')
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
  orbit_center_magnitudes=0.9d0;orbit_center_magnitudes(:,2)=[0.8d0,0.7d0,0.6d0]
  call verify_dg_wannier_center_affine_orbits(orbit_centers,affine_rotations,&
    affine_translations,1d-12,ok,message,moment_magnitudes=orbit_center_magnitudes,&
    failed_operation=failed_center_operation)
  call require(.not.ok.and.has_text(message,'operation=').and.has_text(message,'source=').and.&
    has_text(message,'nearest_residual=').and.has_text(message,'moment_min=').and.failed_center_operation==2,&
    'broken full affine Wannier center orbit reports actionable mismatch diagnostics')
  center_gauge_basis=(0d0,0d0);center_gauge_weights=1d0
  if(rank==0)then
    center_gauge_basis(1,1)=(1d0,0d0);center_gauge_basis(2,2)=(1d0,0d0)
  endif
  center_gauge_map=int(rank*2,8)+[1_8,2_8]
  if(rank==0)center_gauge_map=[2_8,1_8]
  center_gauge_rotation=0
  do i=1,3;center_gauge_rotation(i,i)=1;enddo
  center_gauge_tau=[0.5d0,0d0,0d0]
  center_gauge_centers=0d0;center_gauge_centers(1,:)=[0.25d0,0.75d0]
  call diagnose_dg_point_center_gauge(comm,center_gauge_basis,center_gauge_weights,center_gauge_map,&
    center_gauge_rotation,center_gauge_tau,center_gauge_centers,1d-12,center_gauge_monomial,&
    center_gauge_leakage,center_gauge_unitarity,center_gauge_workspace,ok,message)
  call require(ok.and.center_gauge_monomial<1d-12.and.center_gauge_leakage<1d-12.and.&
    center_gauge_unitarity<1d-12.and.center_gauge_workspace>0_8,&
    'exact center permutation has zero center-gauge leakage')
  if(rank==0)then
    theta=acos(-1d0)/8d0
    center_gauge_basis(:,1)=[cmplx(cos(theta),0d0,8),cmplx(-sin(theta),0d0,8)]
    center_gauge_basis(:,2)=[cmplx(sin(theta),0d0,8),cmplx(cos(theta),0d0,8)]
  endif
  call diagnose_dg_point_center_gauge(comm,center_gauge_basis,center_gauge_weights,center_gauge_map,&
    center_gauge_rotation,center_gauge_tau,center_gauge_centers,1d-12,center_gauge_monomial,&
    center_gauge_leakage,center_gauge_unitarity,center_gauge_workspace,ok,message)
  call require(ok.and.center_gauge_monomial>1d-3.and.center_gauge_leakage>1d-3.and.&
    center_gauge_unitarity<1d-12,'cross-center unitary mixing is diagnosed as leakage')
  center_gauge_tau=0d0;center_gauge_centers=0.25d0
  call diagnose_dg_point_center_gauge(comm,center_gauge_basis,center_gauge_weights,center_gauge_map,&
    center_gauge_rotation,center_gauge_tau,center_gauge_centers,1d-12,center_gauge_monomial,&
    center_gauge_leakage,center_gauge_unitarity,center_gauge_workspace,ok,message)
  call require(ok.and.center_gauge_monomial>1d-3.and.center_gauge_leakage<1d-12,&
    'repeated-center internal rotation is not center-block leakage')
  if(rank==0)write(*,'(a,3(es16.8,1x),i0)')'POINT_CENTER_GAUGE ',center_gauge_monomial,&
    center_gauge_leakage,center_gauge_unitarity,center_gauge_workspace
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
  partition_ids=[int(rank+1,8),int(modulo(rank+1,nproc)+1,8)]
  if(nproc==1)partition_ids=[1_8,2_8]
  raw_partition_weight=[1d0,0.5d0]
  raw_partition_gradient=0d0;raw_partition_gradient(1,:)=[0.2d0,-0.2d0]
  call build_dg_smooth_partition_of_unity(comm,partition_ids,raw_partition_weight,&
    raw_partition_gradient,partition_weight,partition_gradient,partition_sum_defect,&
    partition_gradient_defect,ok,message)
  call require(ok.and.partition_sum_defect<1d-12.and.partition_gradient_defect<1d-12,trim(message))
  if(nproc==1)then
    call require(maxval(abs(partition_weight-1d0))<1d-12,'unique partition points retain unit weight')
  else
    call require(maxval(abs(partition_weight-[2d0/3d0,1d0/3d0]))<1d-12,&
      'partition normalization is independent of fragment rank')
  endif
  allocate(composed_buffer_values(2,2))
  do p=1,2
    composed_buffer_values(1,p)=cmplx(real(partition_ids(p),8),0d0,8)
    composed_buffer_values(2,p)=cmplx(0d0,-real(partition_ids(p),8),8)
  enddo
  call compose_dg_buffered_orbital_tile_to_physical_grid(comm,partition_ids,partition_weight,&
    composed_buffer_values,composed_owned_ids,composed_owned_values,composition_fingerprint,&
    composition_workspace_peak,ok,message)
  call require(ok,trim(message))
  call require(composition_workspace_peak>0_8,'buffer composition publishes nonzero measured workspace')
  call require(size(composed_owned_ids)==size(composed_owned_values,2),&
    'buffer composition returns one value column per owned physical-grid ID')
  do p=1,size(composed_owned_ids)
    call require(maxval(abs(composed_owned_values(:,p)-&
      [cmplx(real(composed_owned_ids(p),8),0d0,8),&
       cmplx(0d0,-real(composed_owned_ids(p),8),8)]))<1d-12,&
      'buffer-first composition recovers the smooth full-system orbital across a fragment face')
  enddo
  deallocate(composed_buffer_values,composed_owned_ids,composed_owned_values)
  partition_ids=[int(rank+1,8),int(rank+nproc+1,8)]
  partition_weight=1d0
  allocate(composed_buffer_values(1,2));composed_buffer_values(1,:)=cmplx(real(partition_ids,8),0d0,8)
  call compose_dg_buffered_orbital_tile_to_physical_grid(comm,partition_ids,partition_weight,&
    composed_buffer_values,composed_owned_ids,composed_owned_values,composition_fingerprint,&
    composition_workspace_peak,ok,message)
  call require(ok,trim(message))
  call require(all(composed_owned_ids==[int(2*rank+1,8),int(2*rank+2,8)]),&
    'buffer composition uses contiguous physical-ID owners required by affine point exchange')
  deallocate(composed_buffer_values,composed_owned_ids,composed_owned_values)
  partition_ids=[int(rank+1,8),int(modulo(rank+1,nproc)+1,8)]
  if(nproc==1)partition_ids=[1_8,2_8]
  if(nproc==1)then
    partition_weight=1d0
  else
    partition_weight=[2d0/3d0,1d0/3d0]
  endif
  partition_weight(1)=0.5d0*partition_weight(1)
  allocate(composed_buffer_values(2,2));composed_buffer_values=(1d0,0d0)
  call compose_dg_buffered_orbital_tile_to_physical_grid(comm,partition_ids,partition_weight,&
    composed_buffer_values,composed_owned_ids,composed_owned_values,composition_fingerprint,&
    composition_workspace_peak,ok,message)
  call require(.not.ok,'buffer composition rejects incomplete partition coverage')
  if(allocated(composed_owned_ids))deallocate(composed_owned_ids)
  if(allocated(composed_owned_values))deallocate(composed_owned_values)
  partition_weight(1)=2d0*partition_weight(1)
  composed_buffer_values(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call compose_dg_buffered_orbital_tile_to_physical_grid(comm,partition_ids,partition_weight,&
    composed_buffer_values,composed_owned_ids,composed_owned_values,composition_fingerprint,&
    composition_workspace_peak,ok,message)
  call require(.not.ok,'buffer composition rejects nonfinite orbital values')
  deallocate(composed_buffer_values)
  allocate(composed_buffer_values(2,2));composed_buffer_values=(1d0,0d0)
  partition_ids(2)=partition_ids(1)
  call compose_dg_buffered_orbital_tile_to_physical_grid(comm,partition_ids,partition_weight,&
    composed_buffer_values,composed_owned_ids,composed_owned_values,composition_fingerprint,&
    composition_workspace_peak,ok,message)
  call require(.not.ok,'buffer composition rejects duplicate physical IDs within one fragment')
  deallocate(composed_buffer_values)
  partition_ids=[int(rank+1,8),int(modulo(rank+1,nproc)+1,8)]
  if(nproc==1)partition_ids=[1_8,2_8]
  raw_partition_weight(2)=-0.5d0
  call build_dg_smooth_partition_of_unity(comm,partition_ids,raw_partition_weight,&
    raw_partition_gradient,partition_weight,partition_gradient,partition_sum_defect,&
    partition_gradient_defect,ok,message)
  call require(.not.ok,'smooth partition rejects a negative fragment window')
  allocate(variable_partition_ids(rank+2),variable_raw_partition_weight(rank+2),&
    variable_raw_partition_gradient(3,rank+2),variable_partition_weight(rank+2),&
    variable_partition_gradient(3,rank+2))
  variable_partition_ids(1)=1_8
  do i=2,rank+2
    variable_partition_ids(i)=1000_8*int(rank+1,8)+int(i,8)
  enddo
  variable_raw_partition_weight=1d0;variable_raw_partition_weight(1)=real(rank+1,8)
  variable_raw_partition_gradient=0d0;variable_raw_partition_gradient(1,1)=real(rank,8)-0.5d0
  call build_dg_smooth_partition_of_unity(comm,variable_partition_ids,variable_raw_partition_weight,&
    variable_raw_partition_gradient,variable_partition_weight,variable_partition_gradient,&
    partition_sum_defect,partition_gradient_defect,ok,message)
  call require(ok.and.partition_sum_defect<1d-12.and.partition_gradient_defect<1d-12,trim(message))
  call require(abs(variable_partition_weight(1)-real(rank+1,8)/&
    real(nproc*(nproc+1)/2,8))<1d-12,'partition supports unequal fragment coverage sizes')
  call require(all(abs(variable_partition_weight(2:)-1d0)<1d-12),&
    'unique points retain unit partition weight')
  deallocate(variable_partition_ids,variable_raw_partition_weight,variable_raw_partition_gradient,&
    variable_partition_weight,variable_partition_gradient)
  do j=1,4;do i=1,4
    cyclic_product(i,j)=mod(i+j-2,4)+1
  enddo;enddo
  call select_dg_group_generators(cyclic_product,1,group_generators,ok,message)
  call require(ok.and.size(group_generators)==1.and.group_generators(1)==2,&
    'deterministic generator selection closes the complete cyclic group')
  cyclic_product(4,4)=5
  call select_dg_group_generators(cyclic_product,1,group_generators,ok,message)
  call require(.not.ok,'generator selection rejects an invalid product table')
  cyclic_product(4,4)=3
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
    call gather_dg_single_symmetry_representation(comm,distributed_overlap_row_ids,&
      distributed_basis_overlap_rows,2,0,single_symmetry_representation,&
      row_overlap_workspace_peak,ok,message)
    call require(ok.and.row_overlap_workspace_peak>0_8,trim(message))
    if(rank==0)then
      call require(maxval(abs(single_symmetry_representation-distributed_basis_overlap(:,:,2)))<1d-12,&
        'one-operation writer gather matches the dense affine reference')
    else
      call require(size(single_symmetry_representation)==0,&
        'nonwriter ranks retain no dense fixed-center representation')
    endif
    call gather_dg_single_symmetry_representation(comm,distributed_overlap_row_ids,&
      distributed_basis_overlap_rows,0,0,single_symmetry_representation,&
      row_overlap_workspace_peak,ok,message)
    call require(.not.ok,'invalid streamed symmetry operation index rejected collectively')
    call gather_dg_single_symmetry_representation(comm,distributed_overlap_row_ids,&
      distributed_basis_overlap_rows,2,nproc,single_symmetry_representation,&
      row_overlap_workspace_peak,ok,message)
    call require(.not.ok,'invalid streamed symmetry writer rank rejected collectively')
    distributed_overlap_row_ids(1)=distributed_overlap_row_ids(1)+1_8
    call gather_dg_single_symmetry_representation(comm,distributed_overlap_row_ids,&
      distributed_basis_overlap_rows,2,0,single_symmetry_representation,&
      row_overlap_workspace_peak,ok,message)
    call require(.not.ok,'noncontiguous streamed symmetry row ownership rejected collectively')
    distributed_overlap_row_ids(1)=distributed_overlap_row_ids(1)-1_8
    if(rank==0)distributed_basis_overlap_rows(1,1,2)=&
      cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
    call gather_dg_single_symmetry_representation(comm,distributed_overlap_row_ids,&
      distributed_basis_overlap_rows,2,0,single_symmetry_representation,&
      row_overlap_workspace_peak,ok,message)
    call require(.not.ok,'nonfinite streamed symmetry input rejected collectively')
    if(rank==0)distributed_basis_overlap_rows(1,1,2)=distributed_basis_overlap(1,1,2)
    call gather_dg_single_symmetry_representation(comm,distributed_overlap_row_ids,&
      distributed_basis_overlap_rows(:,:,1:merge(2,3,rank==0)),2,0,&
      single_symmetry_representation,row_overlap_workspace_peak,ok,message)
    call require(.not.ok,'rank-inconsistent streamed symmetry operation count rejected collectively')
    distributed_basis_overlap_rows(:,:,1:2)=distributed_basis_overlap_rows(:,:,1:2)/5d0
    call validate_dg_row_owned_group_representation(comm,distributed_overlap_row_ids,&
      distributed_basis_overlap_rows(:,:,1:2),closure_product,1,1d-12,row_identity_defect,&
      row_unitarity_defect,row_closure_defect,row_overlap_workspace_peak,ok,message)
    call require(ok.and.max(row_identity_defect,max(row_unitarity_defect,row_closure_defect))<1d-12,&
      trim(message))
    dense_identity_defect=row_identity_defect;dense_unitarity_defect=row_unitarity_defect
    dense_closure_defect=row_closure_defect
    call validate_dg_streamed_affine_representation(comm,distributed_basis/sqrt(5d0),&
      distributed_weight,distributed_map(:,1:2),1,0d0,1d-12,row_identity_defect,&
      row_unitarity_defect,row_closure_defect,row_overlap_workspace_peak,ok,message)
    call require(ok,trim(message))
    call require(abs(row_identity_defect-dense_identity_defect)<1d-12.and.&
      abs(row_unitarity_defect-dense_unitarity_defect)<1d-12.and.&
      abs(row_closure_defect-dense_closure_defect)<1d-12,&
      'streamed affine proof matches the dense two-operation reference')

    allocate(orbit_seed(1,2));orbit_seed=(0d0,0d0)
    if(rank==0)orbit_seed(1,1)=1d0
    call build_dg_group_averaged_occupied_candidates_reference(comm,orbit_seed,distributed_weight,&
      distributed_map(:,1:1),reshape([1],[1,1]),1,1d-12,averaged_candidates,averaged_spectrum,&
      averaged_rank,averaged_trace,averaged_closure,averaged_workspace_peak,ok,message)
    call require(ok.and.averaged_rank==1.and.abs(averaged_spectrum(1)-1d0)<1d-12.and.&
      abs(averaged_trace-1d0)<1d-12.and.averaged_closure<1d-12,&
      'identity-only group average preserves the occupied projector')
    deallocate(averaged_candidates,averaged_spectrum)
    call build_dg_group_averaged_occupied_candidates_reference(comm,orbit_seed,distributed_weight,&
      distributed_map(:,1:2),closure_product,1,1d-12,averaged_candidates,averaged_spectrum,&
      averaged_rank,averaged_trace,averaged_closure,averaged_workspace_peak,ok,message)
    call require(ok.and.averaged_rank==2.and.abs(averaged_trace-1d0)<1d-12,trim(message))
    call require(maxval(abs(averaged_spectrum-[0.5d0,0.5d0]))<1d-12.and.&
      averaged_closure<1d-12.and.averaged_workspace_peak>0_8,&
      'streamed group-averaged occupied projector matches dense two-point reference')
    orbit_gram=matmul(averaged_candidates,conjg(transpose(averaged_candidates)))
    call MPI_Allreduce(MPI_IN_PLACE,orbit_gram,4,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call require(maxval(abs(orbit_gram-reshape([(1d0,0d0),(0d0,0d0),&
      (0d0,0d0),(1d0,0d0)],[2,2])))<1d-12,&
      'group-averaged occupied candidates are globally orthonormal')
    deallocate(averaged_candidates,averaged_spectrum)
    allocate(invalid_orbit_map,source=distributed_map(:,1:2))
    if(rank==0)invalid_orbit_map(:,2)=[3_8,3_8]
    call build_dg_group_averaged_occupied_candidates_reference(comm,orbit_seed,distributed_weight,&
      invalid_orbit_map,closure_product,1,1d-12,averaged_candidates,averaged_spectrum,&
      averaged_rank,averaged_trace,averaged_closure,averaged_workspace_peak,ok,message)
    call require(.not.ok,'group-averaged occupied projector rejects a non-group point action')
    deallocate(invalid_orbit_map)
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
    call orthonormalize_dg_distributed_seed_space(comm,required_orbit_seed,distributed_weight,1d-12,&
      orthonormal_seed_basis,required_orbit_rank,ok,message)
    call require(ok.and.required_orbit_rank==2,trim(message))
    orbit_gram=matmul(orthonormal_seed_basis,conjg(transpose(orthonormal_seed_basis)))
    call MPI_Allreduce(MPI_IN_PLACE,orbit_gram,4,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call require(maxval(abs(orbit_gram-reshape([(1d0,0d0),(0d0,0d0),&
      (0d0,0d0),(1d0,0d0)],[2,2])))<1d-12,&
      'already affine-closed seeds are orthonormalized without orbit expansion')
    deallocate(orthonormal_seed_basis)
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

  character_translations=0d0
  character_translations(1,:)=[0d0,0d0,0.5d0,0.5d0]
  character_translations(2,:)=[0d0,0.5d0,0d0,0.5d0]
  character_product=reshape([1,2,3,4,2,1,4,3,3,4,1,2,4,3,2,1],[4,4])
  call build_dg_finite_abelian_character_table(character_translations,character_product,1,1d-12,&
    character_canonical_operations,character_inverses,character_generator_count,character_generators,&
    character_words,character_table,character_conjugates,character_fingerprint,ok,message)
  call require(ok,trim(message))
  do i=1,12
    do j=1,4
      intertwining_maps(i,j)=int(4*((i-1)/4)+ieor(mod(i-1,4),j-1)+1,8)
    enddo
  enddo
  intertwining_generator_maps=intertwining_maps(:,character_generators)
  call build_dg_translation_character_intertwining_phase(comm,box_ids,12,intertwining_generator_maps(int(box_ids),:),&
    [2,2],character_words,character_product,1,character_table(1,:),character_table(2,:),character_fingerprint,1d-12,&
    intertwining_phase,inverse_fingerprint,intertwining_payload_fingerprint,inverse_workspace,ok,message)
  call require(ok.and.maxval(abs(abs(intertwining_phase)-1d0))<1d-12,&
    'finite translation action builds a unit-modulus character intertwining phase')
  call prepare_dg_translation_character_action(comm,box_ids,12,intertwining_generator_maps(int(box_ids),:),&
    [2,2],character_words,character_product,1,character_fingerprint,1d-12,prepared_translation_action,ok,message)
  call require(ok.and.prepared_translation_action%workspace_peak_bytes>0_8.and.&
    prepared_translation_action%construction_collective_count>0,&
    'finite translation action is prepared once with memory and collective receipts')
  call build_dg_translation_character_intertwining_phase_prepared(comm,prepared_translation_action,&
    character_table(1,:),character_table(2,:),1d-12,prepared_intertwining_phase,&
    prepared_intertwining_fingerprint,prepared_intertwining_payload_fingerprint,prepared_intertwining_workspace,ok,message)
  call require(ok.and.maxval(abs(prepared_intertwining_phase-intertwining_phase))<1d-12.and.&
    prepared_intertwining_fingerprint==inverse_fingerprint.and.&
    prepared_intertwining_payload_fingerprint==intertwining_payload_fingerprint,&
    'prepared translation action reproduces the one-shot phase and provenance exactly')
  call release_dg_prepared_translation_action(prepared_translation_action)
  intertwining_generator_maps(1,1)=intertwining_generator_maps(1,2)
  call build_dg_translation_character_intertwining_phase(comm,box_ids,12,intertwining_generator_maps(int(box_ids),:),&
    [2,2],character_words,character_product,1,character_table(1,:),character_table(2,:),character_fingerprint,1d-12,&
    intertwining_phase,inverse_fingerprint,intertwining_payload_fingerprint,inverse_workspace,ok,message)
  call require(.not.ok,'intertwining phase rejects a non-permutation translation action')
  allocate(inverse_sector_values(nlocal,2,4),inverse_sector_gradients(3,nlocal,2,4))
  inverse_sector_values=(0d0,0d0);inverse_sector_gradients=(0d0,0d0)
  do p=1,nlocal
    do j=1,4;do i=1,2
      if(box_ids(p)==int(2*(j-1)+i,8))inverse_sector_values(p,i,j)=1d0
      inverse_sector_gradients(:,p,i,j)=real(i+2*j,8)*inverse_sector_values(p,i,j)
    enddo;enddo
  enddo
  call inverse_dg_translation_character_orbits(comm,box_ids,12,character_table,character_product,1,&
    character_fingerprint,&
    inverse_sector_values,inverse_sector_gradients,1d-12,inverse_orbit_values,inverse_orbit_gradients,&
    inverse_density_defect,inverse_orthogonality_defect,inverse_gamma_defect,inverse_fingerprint,&
    inverse_workspace,ok,message)
  call require(ok.and.inverse_density_defect<1d-12.and.inverse_orthogonality_defect<1d-12.and.&
    inverse_gamma_defect<1d-12.and.inverse_workspace>0_8,trim(message))
  do j=1,4
    call accumulate_dg_translation_character_orbit_sector(comm,inverse_accumulator,j,character_table,&
      character_fingerprint,inverse_sector_values(:,:,j),inverse_sector_gradients(:,:,:,j),j==1,j==4,1d-12,&
      inverse_streamed_values,inverse_streamed_gradients,inverse_workspace,ok,message)
    call require(ok,trim(message))
  enddo
  call require(maxval(abs(inverse_streamed_values-inverse_orbit_values))<1d-12.and.&
    maxval(abs(inverse_streamed_gradients-inverse_orbit_gradients))<1d-12,&
    'streamed character accumulation matches the complete inverse transform')
  if(nproc>1)then
    call accumulate_dg_translation_character_orbit_sector(comm,inverse_accumulator,1,character_table,&
      character_fingerprint,inverse_sector_values(:,:,1),inverse_sector_gradients(:,:,:,1),rank==0,.false.,1d-12,&
      inverse_streamed_values,inverse_streamed_gradients,inverse_workspace,ok,message)
    call require(.not.ok,'streamed inverse rejects rank-disagreeing initialization')
  endif
  call accumulate_dg_translation_character_orbit_sector(comm,inverse_accumulator,2,character_table,&
    character_fingerprint,inverse_sector_values(:,:,2),inverse_sector_gradients(:,:,:,2),.true.,.false.,1d-12,&
    inverse_streamed_values,inverse_streamed_gradients,inverse_workspace,ok,message)
  call require(ok,trim(message))
  call accumulate_dg_translation_character_orbit_sector(comm,inverse_accumulator,2,character_table,&
    character_fingerprint,inverse_sector_values(:,:,2),inverse_sector_gradients(:,:,:,2),.false.,.false.,1d-12,&
    inverse_streamed_values,inverse_streamed_gradients,inverse_workspace,ok,message)
  call require(.not.ok,'streamed inverse rejects a duplicate character sector')
  call accumulate_dg_translation_character_orbit_sector(comm,inverse_accumulator,3,character_table,&
    character_fingerprint,inverse_sector_values(:,:,3),inverse_sector_gradients(:,:,:,3),.false.,.true.,1d-12,&
    inverse_streamed_values,inverse_streamed_gradients,inverse_workspace,ok,message)
  call require(.not.ok,'streamed inverse finalize rejects missing character sectors')
  call accumulate_dg_translation_character_orbit_sector(comm,inverse_accumulator,1,character_table,&
    character_fingerprint,inverse_sector_values(:,:,1),inverse_sector_gradients(:,:,:,1),.false.,.false.,1d-12,&
    inverse_streamed_values,inverse_streamed_gradients,inverse_workspace,ok,message)
  call require(ok,trim(message))
  call accumulate_dg_translation_character_orbit_sector(comm,inverse_accumulator,4,character_table,&
    character_fingerprint,inverse_sector_values(:,:,4),inverse_sector_gradients(:,:,:,4),.false.,.true.,1d-12,&
    inverse_streamed_values,inverse_streamed_gradients,inverse_workspace,ok,message)
  call require(ok,'streamed inverse completes an out-of-order exactly-once sequence')
  do j=4,1,-1
    call accumulate_dg_translation_character_orbit_sector_values(comm,inverse_accumulator,j,character_table,&
      character_fingerprint,inverse_sector_values(:,:,j),j==4,j==1,1d-12,inverse_values_only,&
      inverse_workspace,ok,message)
    call require(ok,trim(message))
  enddo
  call require(maxval(abs(inverse_values_only-transpose(reshape(inverse_orbit_values,&
    [size(inverse_orbit_values,1),size(inverse_orbit_values,2)*size(inverse_orbit_values,3)]))))<1d-12,&
    'values-only streamed inverse matches the complete transform without dummy gradients')
  allocate(inverse_saved_sector,source=inverse_sector_values(:,:,1))
  inverse_sector_values(:,:,1)=cmplx(0.5d0*huge(1d0),0d0,8)
  call accumulate_dg_translation_character_orbit_sector_values(comm,inverse_accumulator,1,character_table,&
    character_fingerprint,inverse_sector_values(:,:,1),.true.,.false.,1d-12,inverse_values_only,&
    inverse_workspace,ok,message)
  call require(.not.ok,'values-only streamed inverse rejects finite-huge sector values before accumulation')
  inverse_sector_values(:,:,1)=inverse_saved_sector;deallocate(inverse_saved_sector)
  allocate(streamed_transform_rows(nlocal,12),streamed_input_values(12,2),streamed_input_gradients(3,12,2))
  streamed_transform_rows=(0d0,0d0)
  do p=1,nlocal;streamed_transform_rows(p,int(box_ids(p)))=1d0;enddo
  streamed_input_values=reshape([(cmplx(i,0d0,8),i=1,24)],[12,2])
  do i=1,3;streamed_input_gradients(i,:,:)=real(i,8)*streamed_input_values;enddo
  call apply_dg_row_owned_orbital_transform_streamed(comm,box_ids,12,streamed_transform_rows,&
    streamed_input_values,streamed_input_gradients,streamed_output_values,streamed_output_gradients,&
    inverse_workspace,ok,message)
  call require(ok.and.inverse_workspace>0_8.and.maxval(abs(streamed_output_values-streamed_input_values))<1d-12.and.&
    maxval(abs(streamed_output_gradients-streamed_input_gradients))<1d-12,&
    'row-owned transform streams values and gradients without a dense replicated matrix')
  if(mod(12,nproc)==0)then
    allocate(factored_point_maps(nlocal,2),factored_translation_maps(nlocal,2),factored_weights(nlocal))
    do p=1,nlocal
      factored_global_index=rank*nlocal+p
      factored_point_maps(p,1)=int(factored_global_index,8)
      factored_point_maps(p,2)=int(4*((factored_global_index-1)/4)+mod(factored_global_index,4)+1,8)
      factored_translation_maps(p,1)=factored_point_maps(p,1)
      factored_translation_maps(p,2)=int(4*((factored_global_index-1)/4)+mod(factored_global_index+1,4)+1,8)
    enddo
    factored_point_product=reshape([1,2,2,1],[2,2]);factored_cocycle=1
    factored_cocycle(2,2)=2;factored_weights=1d0
    call validate_dg_factored_point_cogroup_gauge(comm,transpose(streamed_transform_rows),factored_weights,&
      factored_point_maps,factored_translation_maps,factored_point_product,1,factored_cocycle,2,0d0,1d-12,&
      factored_identity_defect,factored_unitarity_defect,factored_closure_defect,inverse_workspace,ok,message)
    call require(ok.and.factored_closure_defect<1d-12,&
      'nontrivial factored point-cogroup cocycle proof: '//trim(message))
    factored_cocycle(2,2)=1
    call validate_dg_factored_point_cogroup_gauge(comm,transpose(streamed_transform_rows),factored_weights,&
      factored_point_maps,factored_translation_maps,factored_point_product,1,factored_cocycle,2,0d0,1d-12,&
      factored_identity_defect,factored_unitarity_defect,factored_closure_defect,inverse_workspace,ok,message)
    call require(.not.ok,'factored point-cogroup proof rejects an in-range corrupt translation cocycle')
    factored_cocycle(2,2)=2
    factored_point_maps(1,2)=factored_point_maps(2,2)
    call validate_dg_factored_point_cogroup_gauge(comm,transpose(streamed_transform_rows),factored_weights,&
      factored_point_maps,factored_translation_maps,factored_point_product,1,factored_cocycle,2,0d0,1d-12,&
      factored_identity_defect,factored_unitarity_defect,factored_closure_defect,inverse_workspace,ok,message)
    call require(.not.ok,'factored point-cogroup proof rejects a non-permutation point map')
    factored_global_index=rank*nlocal+1
    factored_point_maps(1,2)=int(4*((factored_global_index-1)/4)+mod(factored_global_index,4)+1,8)
    if(nproc>1)then
      if(rank==0)factored_point_product(2,2)=2
      call validate_dg_factored_point_cogroup_gauge(comm,transpose(streamed_transform_rows),factored_weights,&
        factored_point_maps,factored_translation_maps,factored_point_product,1,factored_cocycle,2,0d0,1d-12,&
        factored_identity_defect,factored_unitarity_defect,factored_closure_defect,inverse_workspace,ok,message)
      call require(.not.ok,'factored point-cogroup proof rejects rank-disagreeing product metadata')
      if(rank==0)factored_point_product(2,2)=1
    endif
    semidirect_units=[1,5,7,11];semidirect_offsets=[0,1,0,0]
    allocate(semidirect_point_maps(nlocal,4),semidirect_translation_maps(nlocal,12),&
      semidirect_product(4,4),semidirect_cocycle(4,4))
    do p=1,nlocal
      factored_global_index=rank*nlocal+p-1
      do i=1,4
        semidirect_point_maps(p,i)=int(mod(semidirect_units(i)*factored_global_index+semidirect_offsets(i),12)+1,8)
      enddo
      do i=1,12
        semidirect_translation_maps(p,i)=int(mod(factored_global_index+i-1,12)+1,8)
      enddo
    enddo
    do semidirect_left=1,4;do semidirect_right=1,4
      semidirect_product_unit=mod(semidirect_units(semidirect_left)*semidirect_units(semidirect_right),12)
      semidirect_product_index=findloc(semidirect_units,semidirect_product_unit,dim=1)
      semidirect_product(semidirect_right,semidirect_left)=semidirect_product_index
      semidirect_shift=modulo(semidirect_offsets(semidirect_right)+semidirect_units(semidirect_right)*&
        semidirect_offsets(semidirect_left)-semidirect_offsets(semidirect_product_index),12)
      semidirect_cocycle(semidirect_right,semidirect_left)=semidirect_shift+1
    enddo;enddo
    call validate_dg_factored_point_cogroup_gauge(comm,transpose(streamed_transform_rows),factored_weights,&
      semidirect_point_maps,semidirect_translation_maps,semidirect_product,1,semidirect_cocycle,12,0d0,1d-12,&
      factored_identity_defect,factored_unitarity_defect,factored_closure_defect,inverse_workspace,ok,message,&
      factored_generator_count,factored_checked_pair_count,factored_prepared_operation_count)
    call require(ok.and.factored_generator_count==2.and.factored_checked_pair_count<16.and.&
      factored_prepared_operation_count==4,&
      'factored cocycle uses a complete generator proof in a noncommuting affine action: '//trim(message))
    call MPI_Allreduce(factored_checked_pair_count,factored_receipt_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(factored_checked_pair_count,factored_receipt_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call require(factored_receipt_min==12.and.factored_receipt_max==12,&
      'generator dense-pair receipt is rank independent and equals 2*g*n-g^2')
    semidirect_product(2,3)=2
    call validate_dg_factored_point_cogroup_gauge(comm,transpose(streamed_transform_rows),factored_weights,&
      semidirect_point_maps,semidirect_translation_maps,semidirect_product,1,semidirect_cocycle,12,0d0,1d-12,&
      factored_identity_defect,factored_unitarity_defect,factored_closure_defect,inverse_workspace,ok,message,&
      factored_generator_count,factored_checked_pair_count)
    call require(.not.ok,'generator proof rejects a nonassociative point product before dense work')
    semidirect_product(2,3)=findloc(semidirect_units,&
      mod(semidirect_units(3)*semidirect_units(2),12),dim=1)
    do p=1,nlocal
      semidirect_point_maps(p,4)=int(mod(int(semidirect_point_maps(p,4)),12)+1,8)
    enddo
    call validate_dg_factored_point_cogroup_gauge(comm,transpose(streamed_transform_rows),factored_weights,&
      semidirect_point_maps,semidirect_translation_maps,semidirect_product,1,semidirect_cocycle,12,0d0,1d-12,&
      factored_identity_defect,factored_unitarity_defect,factored_closure_defect,inverse_workspace,ok,message,&
      factored_generator_count,factored_checked_pair_count)
    call require(.not.ok,'all-pair integer gate rejects a corrupt non-generator representative')
    deallocate(semidirect_point_maps,semidirect_translation_maps,semidirect_product,semidirect_cocycle)
    deallocate(factored_point_maps,factored_translation_maps,factored_weights)
  endif
  if(nproc>1)then
    call inverse_dg_translation_character_orbits(comm,box_ids,12,character_table,character_product,1,&
      character_fingerprint,inverse_sector_values,inverse_sector_gradients,merge(1d-11,1d-12,rank==0),&
      inverse_translated_values,inverse_translated_gradients,inverse_density_defect,inverse_orthogonality_defect,&
      inverse_gamma_defect,inverse_fingerprint,inverse_workspace,ok,message)
    call require(.not.ok,&
      'inverse character transform rejects rank-disagreeing metadata collectively')
  endif
  allocate(inverse_row_ids,source=box_ids)
  if(rank==0.and.size(inverse_row_ids)>1)inverse_row_ids(1)=inverse_row_ids(2)
  call inverse_dg_translation_character_orbits(comm,inverse_row_ids,12,character_table,character_product,1,&
    character_fingerprint,inverse_sector_values,inverse_sector_gradients,1d-12,inverse_translated_values,&
    inverse_translated_gradients,inverse_density_defect,inverse_orthogonality_defect,inverse_gamma_defect,&
    inverse_fingerprint,inverse_workspace,ok,message)
  call require(.not.ok,&
    'inverse character transform rejects duplicate or missing row ownership')
  inverse_reference_fingerprint=inverse_fingerprint
  allocate(inverse_density(nlocal),inverse_rotated_density(nlocal))
  inverse_density=sum(sum(abs(inverse_orbit_values)**2,dim=3),dim=2)
  inverse_internal_gauge=reshape([cmplx(cos(0.37d0),0d0,8),cmplx(-sin(0.37d0),0d0,8),&
    cmplx(sin(0.37d0),0d0,8),cmplx(cos(0.37d0),0d0,8)],[2,2])
  allocate(inverse_rotated_values(nlocal,2,4),inverse_rotated_gradients(3,nlocal,2,4))
  do j=1,4
    inverse_rotated_values(:,:,j)=matmul(inverse_sector_values(:,:,j),inverse_internal_gauge)
    do i=1,3
      inverse_rotated_gradients(i,:,:,j)=matmul(inverse_sector_gradients(i,:,:,j),inverse_internal_gauge)
    enddo
  enddo
  call inverse_dg_translation_character_orbits(comm,box_ids,12,character_table,character_product,1,&
    character_fingerprint,&
    inverse_rotated_values,inverse_rotated_gradients,1d-12,inverse_translated_values,&
    inverse_translated_gradients,inverse_density_defect,inverse_orthogonality_defect,inverse_gamma_defect,&
    inverse_fingerprint,inverse_workspace,ok,message)
  inverse_rotated_density=sum(sum(abs(inverse_translated_values)**2,dim=3),dim=2)
  call require(ok.and.maxval(abs(inverse_rotated_density-inverse_density))<1d-12,&
    'inverse character orbit density is invariant under internal-channel gauge rotation')
  do j=1,4
    inverse_rotated_values(:,:,j)=character_table(j,2)*inverse_sector_values(:,:,j)
    inverse_rotated_gradients(:,:,:,j)=character_table(j,2)*inverse_sector_gradients(:,:,:,j)
  enddo
  call inverse_dg_translation_character_orbits(comm,box_ids,12,character_table,character_product,1,&
    character_fingerprint,&
    inverse_rotated_values,inverse_rotated_gradients,1d-12,inverse_translated_values,&
    inverse_translated_gradients,inverse_density_defect,inverse_orthogonality_defect,inverse_gamma_defect,&
    inverse_fingerprint,inverse_workspace,ok,message)
  do j=1,4
    i=character_product(j,2)
    call require(maxval(abs(inverse_translated_values(:,:,j)-inverse_orbit_values(:,:,i)))<1d-12.and.&
      maxval(abs(inverse_translated_gradients(:,:,:,j)-inverse_orbit_gradients(:,:,:,i)))<1d-12,&
      'inverse character transform gives the known Z2xZ2 translation permutation')
  enddo
  call require(ok,'translation permutation preserves the inverse-character orbit contracts')
  character_table(2,2)=character_table(2,2)*exp(cmplx(0d0,0.1d0,8))
  call inverse_dg_translation_character_orbits(comm,box_ids,12,character_table,character_product,1,&
    character_fingerprint,&
    inverse_sector_values,inverse_sector_gradients,1d-12,inverse_translated_values,&
    inverse_translated_gradients,inverse_density_defect,inverse_orthogonality_defect,inverse_gamma_defect,&
    inverse_fingerprint,inverse_workspace,ok,message)
  call require(.not.ok,'inverse character transform rejects a corrupt character phase')
  character_table(2,2)=cmplx(1d300,0d0,8)
  call inverse_dg_translation_character_orbits(comm,box_ids,12,character_table,character_product,1,&
    character_fingerprint,inverse_sector_values,inverse_sector_gradients,1d-12,inverse_translated_values,&
    inverse_translated_gradients,inverse_density_defect,inverse_orthogonality_defect,inverse_gamma_defect,&
    inverse_fingerprint,inverse_workspace,ok,message)
  call require(.not.ok,'inverse character transform rejects a finite huge character collectively')
  character_product(2,2)=5
  call inverse_dg_translation_character_orbits(comm,box_ids,12,character_table,character_product,1,&
    character_fingerprint,inverse_sector_values,inverse_sector_gradients,1d-12,inverse_translated_values,&
    inverse_translated_gradients,inverse_density_defect,inverse_orthogonality_defect,inverse_gamma_defect,&
    inverse_fingerprint,inverse_workspace,ok,message)
  call require(.not.ok,'inverse character transform safely rejects an out-of-range product entry')
  character_product=reshape([1,2,3,4,2,1,4,3,3,4,1,2,4,3,2,1],[4,4])
  call build_dg_finite_abelian_character_table(character_translations,character_product,1,1d-12,&
    character_canonical_operations,character_inverses,character_generator_count,character_generators,&
    character_words,character_table,character_conjugates,character_fingerprint,ok,message)
  inverse_sector_gradients(1,1,1,1)=cmplx(1d300,0d0,8)
  call inverse_dg_translation_character_orbits(comm,box_ids,12,character_table,character_product,1,&
    character_fingerprint,inverse_sector_values,inverse_sector_gradients,1d-12,inverse_translated_values,&
    inverse_translated_gradients,inverse_density_defect,inverse_orthogonality_defect,inverse_gamma_defect,&
    inverse_fingerprint,inverse_workspace,ok,message)
  call require(.not.ok,'inverse character transform rejects a finite huge sector gradient before squaring')

  character_translations=0d0;character_translations(1,:)=[0d0,0.5d0,0.25d0,0.75d0]
  character_product=reshape([1,2,3,4,2,1,4,3,3,4,2,1,4,3,1,2],[4,4])
  call build_dg_finite_abelian_character_table(character_translations,character_product,1,1d-12,&
    character_canonical_operations,character_inverses,character_generator_count,character_generators,&
    character_words,character_table,character_conjugates,character_fingerprint,ok,message)
  call require(ok,trim(message))
  do i=1,4;character_inverse_permutation(character_canonical_operations(i))=i;enddo
  do i=1,4;do j=1,4
    character_permuted_product(i,j)=character_inverse_permutation(&
      character_product(character_canonical_operations(i),character_canonical_operations(j)))
  enddo;enddo
  inverse_sector_values=(0d0,0d0);inverse_sector_gradients=(0d0,0d0)
  do p=1,nlocal
    if(box_ids(p)==1_8)inverse_sector_values(p,1,1)=1d0
    if(box_ids(p)==2_8)inverse_sector_values(p,1,3)=1d0
    if(box_ids(p)==3_8)then
      inverse_sector_values(p,1,2)=sqrt(0.5d0);inverse_sector_values(p,1,4)=sqrt(0.5d0)
    endif
    if(box_ids(p)==4_8)then
      inverse_sector_values(p,1,2)=cmplx(0d0,sqrt(0.5d0),8)
      inverse_sector_values(p,1,4)=cmplx(0d0,-sqrt(0.5d0),8)
    endif
    do j=1,4;do i=1,3
      inverse_sector_gradients(i,p,1,j)=real(i,8)*inverse_sector_values(p,1,j)
    enddo;enddo
  enddo
  call inverse_dg_translation_character_orbits(comm,box_ids,12,character_table,character_permuted_product,1,&
    character_fingerprint,&
    inverse_sector_values(:,1:1,:),inverse_sector_gradients(:,:,1:1,:),1d-12,inverse_translated_values,&
    inverse_translated_gradients,inverse_density_defect,inverse_orthogonality_defect,inverse_gamma_defect,&
    inverse_fingerprint,inverse_workspace,ok,message)
  call require(ok.and.inverse_gamma_defect<1d-12.and.inverse_orthogonality_defect<1d-12,&
    'Z4 inverse character transform produces real orthonormal Gamma orbits')
  inverse_reference_fingerprint=inverse_fingerprint
  inverse_orbit_values=inverse_translated_values;inverse_orbit_gradients=inverse_translated_gradients
  do j=1,4
    inverse_rotated_values(:,1,j)=character_table(j,2)*inverse_sector_values(:,1,j)
    inverse_rotated_gradients(:,:,1,j)=character_table(j,2)*inverse_sector_gradients(:,:,1,j)
  enddo
  call inverse_dg_translation_character_orbits(comm,box_ids,12,character_table,character_permuted_product,1,&
    character_fingerprint,inverse_rotated_values(:,1:1,:),inverse_rotated_gradients(:,:,1:1,:),1d-12,&
    inverse_translated_values,inverse_translated_gradients,inverse_density_defect,inverse_orthogonality_defect,&
    inverse_gamma_defect,inverse_fingerprint,inverse_workspace,ok,message)
  do j=1,4
    i=character_permuted_product(j,character_inverses(2))
    call require(maxval(abs(inverse_translated_values(:,:,j)-inverse_orbit_values(:,:,i)))<1d-12.and.&
      maxval(abs(inverse_translated_gradients(:,:,:,j)-inverse_orbit_gradients(:,:,:,i)))<1d-12,&
      'Z4 generator acts by the documented inverse translation permutation')
  enddo
  inverse_internal_gauge=reshape([cmplx(sqrt(0.5d0),0d0,8),cmplx(0d0,sqrt(0.5d0),8),&
    cmplx(0d0,sqrt(0.5d0),8),cmplx(sqrt(0.5d0),0d0,8)],[2,2])
  inverse_sector_values=(0d0,0d0);inverse_sector_gradients=(0d0,0d0)
  do p=1,nlocal
    if(box_ids(p)<=8_8)then
      j=(int(box_ids(p))-1)/2+1;i=mod(int(box_ids(p))-1,2)+1
      do b=1,4
        inverse_sector_values(p,i,b)=0.5d0*character_table(b,j)
        inverse_sector_gradients(:,p,i,b)=real(i+j,8)*inverse_sector_values(p,i,b)
      enddo
    endif
  enddo
  call inverse_dg_translation_character_orbits(comm,box_ids,12,character_table,character_permuted_product,1,&
    character_fingerprint,inverse_sector_values,inverse_sector_gradients,1d-12,inverse_translated_values,&
    inverse_translated_gradients,inverse_density_defect,inverse_orthogonality_defect,inverse_gamma_defect,&
    inverse_reference_fingerprint,inverse_workspace,ok,message)
  call require(ok,'ungauged Z4 repeated-channel inverse transform')
  inverse_rotated_values=inverse_sector_values;inverse_rotated_gradients=inverse_sector_gradients
  inverse_rotated_values(:,:,2)=matmul(inverse_sector_values(:,:,2),inverse_internal_gauge)
  inverse_rotated_values(:,:,4)=matmul(inverse_sector_values(:,:,4),conjg(inverse_internal_gauge))
  do i=1,3
    inverse_rotated_gradients(i,:,:,2)=matmul(inverse_sector_gradients(i,:,:,2),inverse_internal_gauge)
    inverse_rotated_gradients(i,:,:,4)=matmul(inverse_sector_gradients(i,:,:,4),conjg(inverse_internal_gauge))
  enddo
  call inverse_dg_translation_character_orbits(comm,box_ids,12,character_table,character_permuted_product,1,&
    character_fingerprint,inverse_rotated_values,inverse_rotated_gradients,1d-12,inverse_translated_values,&
    inverse_translated_gradients,inverse_density_defect,inverse_orthogonality_defect,inverse_gamma_defect,&
    inverse_fingerprint,inverse_workspace,ok,message)
  call require(ok.and.inverse_gamma_defect<1d-12.and.inverse_density_defect<1d-12,&
    'Z4 conjugate sectors preserve Gamma reality under a complex internal gauge')
  call require(inverse_fingerprint==inverse_reference_fingerprint,&
    'inverse character fingerprint is invariant under conjugate complex internal gauges')

  b=count([(mod(i-1,nproc)==rank,i=1,4)])
  allocate(spatial_basis(4,nlocal),sector_coefficient_ids(b),sector_coefficients(b,2))
  do p=1,nlocal
    do i=1,4;spatial_basis(i,p)=cmplx(real(10*i+p,8),real(i-p,8),8);enddo
  enddo
  b=0
  do i=1,4
    if(mod(i-1,nproc)/=rank)cycle
    b=b+1;sector_coefficient_ids(b)=i
    sector_coefficients(b,:)=[cmplx(real(i,8),0d0,8),cmplx(0d0,real(i,8),8)]
  enddo
  call materialize_dg_row_owned_sector_on_spatial_grid(comm,sector_coefficient_ids,4,&
    sector_coefficients,spatial_basis,987_8,materialized_sector,spatial_sector_fingerprint,inverse_workspace,ok,message)
  call require(ok.and.maxval(abs(materialized_sector(:,1)-sum(spatial_basis*&
    spread([(real(i,8),i=1,4)],2,nlocal),dim=1)))<1d-12.and.&
    maxval(abs(materialized_sector(:,2)-cmplx(0d0,1d0,8)*sum(spatial_basis*&
    spread([(real(i,8),i=1,4)],2,nlocal),dim=1)))<1d-12,&
    'row-owned character sector materializes on the local spatial grid without dense coefficient replication')

  spectral_eigenvalues=[-2d0,-1d0,0.1d0,0.2d0,0.2d0,0.2d0,0.8d0,1.0d0,1.2d0,1.4d0]
  spectral_occupations=[2d0,2d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0]
  call build_dg_equal_count_spectral_windows(comm,spectral_eigenvalues,spectral_occupations,3,1d-10,&
    spectral_window_weights,spectral_window_fingerprint,spectral_window_workspace,ok,message)
  call require(ok.and.all(spectral_window_weights>=0d0).and.&
    maxval(abs(sum(spectral_window_weights(3:,:),dim=2)-1d0))<1d-12.and.&
    any(spectral_window_weights(3:,:)>0d0.and.spectral_window_weights(3:,:)<1d0).and.&
    maxval(abs(spectral_window_weights(4,:)-spectral_window_weights(5,:)))<1d-14.and.&
    maxval(abs(spectral_window_weights(5,:)-spectral_window_weights(6,:)))<1d-14,&
    'equal-count spectral windows preserve a boundary degeneracy and partition unoccupied states')
  if(nproc>1)then
    spectral_eigenvalues(6)=merge(0.2d0,0.21d0,rank==0)
    call build_dg_equal_count_spectral_windows(comm,spectral_eigenvalues,spectral_occupations,3,1d-10,&
      spectral_window_weights,spectral_window_fingerprint,spectral_window_workspace,ok,message)
    call require(.not.ok,'rank-disagreeing spectral eigenvalues are collectively rejected')
  endif
  spectral_eigenvalues(6)=0.2d0
  spectral_eigenvalues(5)=ieee_value(0d0,ieee_quiet_nan)
  call build_dg_equal_count_spectral_windows(comm,spectral_eigenvalues,spectral_occupations,3,1d-10,&
    spectral_window_weights,spectral_window_fingerprint,spectral_window_workspace,ok,message)
  call require(.not.ok,'nonfinite spectral eigenvalues are collectively rejected')
  spectral_eigenvalues(5)=0.2d0;spectral_eigenvalues(7)=0.15d0
  call build_dg_equal_count_spectral_windows(comm,spectral_eigenvalues,spectral_occupations,3,1d-10,&
    spectral_window_weights,spectral_window_fingerprint,spectral_window_workspace,ok,message)
  call require(.not.ok,'unsorted spectral eigenvalues are collectively rejected')
  spectral_eigenvalues(7)=0.8d0
  b=count([(mod(i-1,nproc)==rank,i=1,4)])
  if(allocated(spectral_window_weights))deallocate(spectral_window_weights)
  allocate(spectral_row_ids(b),spectral_state_values(4,b),spectral_window_weights(4,1))
  b=0
  do i=1,4
    if(mod(i-1,nproc)/=rank)cycle
    b=b+1;spectral_row_ids(b)=i;spectral_state_values(:,b)=(0d0,0d0)
    spectral_state_values(1,b)=merge((1d0,0d0),(0d0,0d0),i==1)
    spectral_state_values(2,b)=merge((1d0,0d0),(0d0,0d0),i==2)
    spectral_state_values(3,b)=merge((1d0,0d0),(0d0,0d0),i==3)
    spectral_state_values(4,b)=merge((1d0,0d0),(0d0,0d0),i==4)
  enddo
  spectral_occupations(1:4)=[1d0,1d0,0d0,0d0];spectral_window_weights=0d0
  spectral_window_weights(3:4,1)=1d0
  call build_dg_spectral_density_descriptors(comm,spectral_row_ids,4,spectral_state_values,&
    spectral_occupations(1:4),spectral_window_weights,1d-12,spectral_occupied_density,&
    spectral_unoccupied_density,spectral_shared_density,spectral_density_fingerprint,&
    spectral_density_workspace,ok,message)
  call require(ok.and.maxloc(spectral_occupied_density,dim=1)<=size(spectral_occupied_density).and.&
    maxval(spectral_shared_density)<1d-12,'ionic occupied and electron descriptors remain spatially separated')
  allocate(spectral_reference_descriptors(size(spectral_occupied_density),3))
  spectral_reference_descriptors(:,1)=spectral_occupied_density
  spectral_reference_descriptors(:,2)=spectral_unoccupied_density(:,1)
  spectral_reference_descriptors(:,3)=spectral_shared_density(:,1)
  spectral_state_values(1:2,:)=matmul(reshape([cmplx(sqrt(0.5d0),0d0,8),cmplx(0d0,sqrt(0.5d0),8),&
    cmplx(0d0,sqrt(0.5d0),8),cmplx(sqrt(0.5d0),0d0,8)],[2,2]),spectral_state_values(1:2,:))
  spectral_state_values(3:4,:)=matmul(reshape([cmplx(sqrt(0.5d0),0d0,8),cmplx(0d0,sqrt(0.5d0),8),&
    cmplx(0d0,sqrt(0.5d0),8),cmplx(sqrt(0.5d0),0d0,8)],[2,2]),spectral_state_values(3:4,:))
  call build_dg_spectral_density_descriptors(comm,spectral_row_ids,4,spectral_state_values,&
    spectral_occupations(1:4),spectral_window_weights,1d-12,spectral_occupied_density,&
    spectral_unoccupied_density,spectral_shared_density,spectral_density_fingerprint,&
    spectral_density_workspace,ok,message)
  call require(ok.and.maxval(abs(spectral_occupied_density-spectral_reference_descriptors(:,1)))<1d-12.and.&
    maxval(abs(spectral_unoccupied_density(:,1)-spectral_reference_descriptors(:,2)))<1d-12,&
    'spectral densities are invariant under unitary rotations inside equal-weight state blocks')

  if(rank==0)then
    write(*,'(a,i0)')'INVERSE_CHARACTER_FINGERPRINT ',inverse_fingerprint
    write(*,'(a,i0)')'SPATIAL_SECTOR_FINGERPRINT ',spatial_sector_fingerprint
    write(*,'(a,i0)')'SPECTRAL_WINDOW_FINGERPRINT ',spectral_window_fingerprint
    write(*,'(a,i0)')'SPECTRAL_DENSITY_FINGERPRINT ',spectral_density_fingerprint
    write(*,'(a,i0,a,i0,a,*(i0,1x))')'CONSTRUCTION ranks=',nproc,' fingerprint=',&
      reference_fingerprint,' centers=',reference_center_box_ids
    write(*,'(a,i0,a)')'PASS overlapping-Wannier construction on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  logical function has_text(text,pattern) result(found)
    intrinsic::index
    character(*),intent(in)::text,pattern
    found=index(text,pattern)>0
  end function
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
