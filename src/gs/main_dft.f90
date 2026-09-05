!
!  Copyright 2019-2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!=======================================================================

#include "config.h"

subroutine main_dft
use iso_fortran_env,only:int64,error_unit
use,intrinsic::ieee_arithmetic,only:ieee_is_finite
use math_constants, only: pi, zi
#ifdef USE_MPI
use mpi, only: MPI_Allreduce,MPI_Allgather,MPI_Allgatherv,MPI_Bcast,MPI_IN_PLACE,MPI_INTEGER8,&
  MPI_BXOR,MPI_SUCCESS,MPI_INTEGER,MPI_DOUBLE_COMPLEX,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_MIN,&
  MPI_MAX,MPI_COMM_SELF,MPI_2INTEGER,MPI_MINLOC,MPI_LOGICAL
use omp_lib, only: omp_get_max_threads
#endif
use structures
use inputoutput
use salmon_global, only: yn_dc_lcfo_flux, yn_dc_lcfo_wannier, yn_dg_hybrid_scf, &
  yn_dg_hybrid_continuation_scf, &
  yn_dg_hybrid_divided_scf, &
  yn_dg_dc_overlapping_wannier, ncg, base_directory, num_fragment, &
  dg_dc_metric_rank_tolerance, &
  dg_dc_gs_intermediate_orbital_tolerance,dg_dc_gs_intermediate_density_tolerance, &
  dg_dc_gs_final_orbital_tolerance,dg_dc_gs_final_density_tolerance,dg_dc_gs_subspace_tolerance, &
  dg_dc_gs_initial_lambda_step,dg_dc_gs_minimum_lambda_step,dg_dc_gs_maximum_lambda_step, &
  dg_dc_gs_allowed_residual_growth,dg_dc_gs_density_mix_rate,dg_dc_gs_hermiticity_tolerance, &
  dg_dc_gs_orthogonality_tolerance,dg_dc_gs_face_balance_tolerance,dg_dc_gs_electron_count_tolerance, &
  dg_dc_gs_minimum_projector_overlap,dg_dc_gs_maximum_scf_iterations, &
  dg_dc_gs_maximum_eigensolver_iterations,dg_dc_gs_maximum_rollbacks, &
  dg_dc_gs_sipg_penalty_factor,dg_dc_gs_target_lambda,dg_ow_boundary_value_tolerance,&
  dg_ow_boundary_gradient_tolerance,dg_ow_symmetry_tolerance,&
  dg_ow_localization_support_tolerance,dg_ow_localization_spread_tolerance,&
  dg_ow_localization_gradient_tolerance,dg_ow_localization_max_iterations,&
  dg_ow_candidate_states_per_fragment,dg_ow_target_wanniers_per_fragment,wannier_num_iter,&
  dg_ow_w90_initial_projection,wannier_pw_cutoff,wannier_pw_max,nscf,method_mixing,&
  dg_dc_seed_mode,dg_dc_seed_directory,dg_hybrid_symmetry_energy_window,temperature,&
  dg_hybrid_fragment_cg_steps,dg_hybrid_divided_mixing
use dg_dc_seed_checkpoint,only:s_dg_dc_seed_contract,s_dg_dc_seed_payload,&
  DG_DC_SEED_ABSENT,DG_DC_SEED_VALID,build_dg_dc_seed_contract,probe_dg_dc_seed,&
  read_dg_dc_seed,write_dg_dc_seed,restore_dg_dc_seed_payload,resolve_dg_dc_seed_mode
use dg_overlapping_wannier_construction, only: s_dg_overlapping_wannier_construction, &
  construct_dg_overlapping_wannier_basis,verify_dg_overlapping_wannier_periodic_closure,&
  replicate_dg_fragment_wannier_representative,verify_dg_fragment_wannier_streaming_closure,&
  verify_dg_fragment_center_orbit,verify_dg_uniform_fragment_target_rank
use dg_overlapping_wannier_construction, only: build_dg_core_owned_occupied_subspace
use dg_overlapping_wannier_construction, only: orthonormalize_dg_distributed_seed_space
use dg_overlapping_wannier_construction, only: compose_dg_buffered_orbital_tile_to_physical_grid
use dg_overlapping_wannier_construction, only: measure_dg_rank_fixed_symmetry_residuals
#ifdef USE_EIGENEXA
use dg_overlapping_wannier_construction, only: measure_dg_rank_fixed_symmetry_residuals_eigenexa
use dg_overlapping_wannier_construction, only: build_dg_group_averaged_occupied_candidates_eigenexa
use dg_overlapping_wannier_construction, only: build_dg_cocycle_averaged_occupied_candidates_eigenexa
use dg_overlapping_wannier_construction, only: split_dg_translation_character_sector_eigenexa
#endif
use dg_overlapping_wannier_construction, only: select_dg_fixed_rank_symmetry_closed_subspace
use dg_overlapping_wannier_construction, only: find_dg_group_identity
use dg_overlapping_wannier_construction, only: select_dg_group_generators
use dg_overlapping_wannier_construction, only: build_dg_finite_abelian_character_table
use dg_overlapping_wannier_construction, only: s_dg_translation_orbit_accumulator,&
  s_dg_prepared_translation_action,s_dg_prepared_spectral_basins,&
  prepare_dg_translation_character_action,&
  build_dg_translation_character_intertwining_phase_prepared,release_dg_prepared_translation_action,&
  materialize_dg_row_owned_sector_on_spatial_grid,&
  accumulate_dg_translation_character_orbit_sector_values,&
  apply_dg_row_owned_orbital_transform_streamed,&
  validate_dg_factored_point_cogroup_gauge
use dg_overlapping_wannier_construction, only: build_dg_occupied_empty_moment_descriptors,&
  build_dg_periodic_spectral_basins,prepare_dg_spectral_basin_operators,&
  project_dg_prepared_spectral_basin_operator,release_dg_prepared_spectral_basins,&
  diagonalize_dg_spectral_basin_operator,select_dg_spectral_basin_channel_ranks,&
  propagate_dg_spectral_basin_orbit_channels,build_dg_spectral_channel_generator_actions,&
  compose_dg_occupied_complement_trial_rows,prepare_dg_direct_retained_wannier_frame
use dg_overlapping_wannier_construction, only: build_dg_smooth_partition_of_unity,&
  redistribute_dg_row_owned_real_field_to_requests
use dg_overlapping_wannier_construction, only: assemble_dg_distributed_basis_symmetry_overlap
use dg_overlapping_wannier_construction, only: assemble_dg_distributed_basis_symmetry_overlap_rows,&
  gather_dg_single_symmetry_representation
use dg_overlapping_wannier_construction, only: build_dg_pointwise_affine_owner_map
use dg_overlapping_wannier_construction, only: solve_dg_affine_common_fixed_point
use dg_overlapping_wannier_construction, only: compute_dg_periodic_wannier_centers
use dg_overlapping_wannier_construction, only: verify_dg_wannier_center_affine_orbits,&
  diagnose_dg_point_center_gauge
use dg_overlapping_wannier_construction, only: transpose_dg_spatial_cores_to_orbital_owners,&
  exchange_dg_point_permuted_orbital_rows,reindex_dg_point_maps_between_row_layouts,&
  measure_dg_spatial_basis_covariance,&
  measure_dg_spatial_gradient_covariance,&
  measure_dg_grid_map_stencil_defect,&
  redistribute_dg_owned_orbitals_to_center_fragments,&
  assign_dg_periodic_centers_to_fragments
use dg_overlapping_wannier_projection, only: t_dg_projection_channel,&
  build_dg_complete_sp_manifest,evaluate_dg_periodic_sp_projectors,&
  dg_periodic_grid_point_owned,select_dg_sp_atomic_orbital_ordinals
use dg_overlapping_wannier_metric, only: assemble_dg_stitched_overlap_density_rows
use dg_overlapping_wannier_operators, only: assemble_dg_stitched_weak_operator_rows
use dg_overlapping_wannier_nonlocal, only: assemble_dg_overlapping_wannier_nonlocal,&
  assemble_dg_overlapping_wannier_nonlocal_rows,collect_dg_overlapping_wannier_projector_overlaps,&
  apply_dg_overlapping_wannier_nonlocal_action
use dg_overlapping_wannier_scf, only: s_dg_overlapping_wannier_scf_state, &
  s_dg_overlapping_wannier_scf_result, &
  compute_dg_overlapping_wannier_scf_fingerprint,mix_dg_overlapping_wannier_density_history
use dg_hybrid_scf,only:run_dg_hybrid_self_consistent_ground_state
use dg_hybrid_windowed_pw_types,only:s_dg_hybrid_basis_catalog,s_dg_hybrid_production_selection
use dg_hybrid_production_pw_basis,only:build_dg_hybrid_production_pw_basis
use dg_hybrid_window_distribution,only:redistribute_dg_hybrid_fragment_windows
use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
use dg_hybrid_fragment_wannier,only:s_dg_hybrid_fragment_wannier_cache,&
  build_dg_hybrid_fragment_wannier_from_dc_seed
use dg_hybrid_fragment_selection,only:s_dg_hybrid_core_selection,s_dg_hybrid_selected_catalog,&
  s_dg_hybrid_dc_reference,select_dg_hybrid_core_wannier,prepare_dg_hybrid_selected_catalog,&
  export_dg_hybrid_dc_reference
use dg_hybrid_fragment_admission,only:s_dg_hybrid_support_operator,s_dg_hybrid_admission_report,&
  prepare_dg_hybrid_selected_trial,export_dg_hybrid_selected_basis_frame
use dg_hybrid_projected_fragment_pipeline,only:s_dg_hybrid_projection_factorization_receipt,&
  s_dg_hybrid_core_projection_report,project_dg_hybrid_core_seeds
use dg_hybrid_production_support,only:s_dg_hybrid_production_support_receipt,&
  prepare_dg_hybrid_production_support
use dg_hybrid_fragment_subspace,only:s_dg_hybrid_fragment_subspace_state,&
  s_dg_hybrid_fragment_candidate_catalog,s_dg_hybrid_fragment_epoch_budget,fragment_seed,fragment_pw
use dg_hybrid_fragment_thermal,only:s_dg_hybrid_thermal_state,advance_dg_hybrid_thermal_state
use dg_hybrid_fragment_preconditioner,only:s_dg_hybrid_fragment_preconditioner,&
  s_dg_hybrid_preconditioner_key,prepare_dg_hybrid_frame_preconditioner,&
  apply_dg_hybrid_fragment_preconditioner
use dg_hybrid_divided_operator,only:freeze_dg_hybrid_single_owner_payload,&
  extract_dg_hybrid_fragment_self_block
use dg_hybrid_production_face_traces,only:s_dg_hybrid_production_face_trace,&
  freeze_dg_hybrid_basis_directory,materialize_dg_hybrid_production_face_collection,&
  assemble_dg_hybrid_production_interface_component_rows,materialize_dg_hybrid_production_interior,&
  reconstruct_dg_hybrid_production_interface_state,reconstruct_dg_hybrid_production_interface_actions
use dg_hybrid_broken_volume,only:assemble_dg_hybrid_broken_volume_rows,assemble_dg_hybrid_local_potential_rows
use dg_hybrid_variational_payload,only:s_dg_hybrid_fixed_payload,s_dg_hybrid_variational_iterate,&
  freeze_dg_hybrid_variational_payload,compose_dg_hybrid_variational_hamiltonian,&
  write_dg_hybrid_variational_payload_bundle
use dg_hybrid_continuation_residuals,only:s_dg_hybrid_residuals,evaluate_dg_hybrid_residuals
use dg_hybrid_real_space_residual,only:evaluate_dg_hybrid_real_space_residual,&
  evaluate_dg_hybrid_face_action_residuals
use dg_hybrid_continuation_controller,only:s_dg_hybrid_controller_controls,s_dg_hybrid_trial_state,&
  s_dg_hybrid_stage_schedule,&
  s_dg_hybrid_stage_report,s_dg_hybrid_controller,s_dg_hybrid_candidate_acceptance,&
  default_dg_hybrid_controller_controls,&
  dg_hybrid_stage_tolerances,initialize_dg_hybrid_controller,propose_dg_hybrid_trial,&
  observe_dg_hybrid_inner_residuals,decide_dg_hybrid_stage,reject_dg_hybrid_trial,&
  initialize_dg_hybrid_stage_schedule,begin_dg_hybrid_stage_solve,&
  schedule_dg_hybrid_candidate_checks,complete_dg_hybrid_stage_solve,&
  initialize_dg_hybrid_candidate_acceptance,&
  record_dg_hybrid_complete_lcfo_solve,record_dg_hybrid_occupation_policy,&
  record_dg_hybrid_unconditional_gates,record_dg_hybrid_spectral_certification,&
  record_dg_hybrid_certified_rt_basis,authorize_dg_hybrid_v3_publication
use dg_hybrid_low_energy_symmetry,only:evaluate_dg_hybrid_low_energy_symmetry,&
  certify_dg_hybrid_energy_window
use dg_hybrid_localization_first,only:s_dg_hybrid_localization_receipt,&
  prepare_dg_hybrid_localization_first_seed,build_dg_hybrid_localization_receipt
use dg_hybrid_continuation_state,only:s_dg_hybrid_scope_receipt,build_dg_hybrid_scope_receipt,&
  close_dg_hybrid_selection
use plusU_global,only:PLUS_U_ON
use dg_hybrid_projected_fragment_pipeline,only:build_dg_hybrid_projected_fragment_basis,&
  build_dg_hybrid_projected_local_fragment_basis
use dg_hybrid_fragment_solver,only:solve_dg_hybrid_fragment_spectrum,&
  reconstruct_dg_hybrid_fragment_density
use dc_fragment_occupation,only:determine_dc_fragment_occupations
use dg_hybrid_divided_scf,only:run_dg_hybrid_divided_scf
use dg_hybrid_divided_mixing,only:s_dg_hybrid_divided_mixing_state,&
  prepare_dg_hybrid_divided_mixing,accept_dg_hybrid_divided_mixing
use dg_hybrid_schwarz_state,only:s_dg_hybrid_schwarz_state,initialize_dg_hybrid_schwarz_state
use dg_hybrid_schwarz_operator,only:s_dg_hybrid_schwarz_schedule,&
  build_dg_hybrid_schwarz_schedule,apply_dg_hybrid_schwarz_hamiltonian,apply_dg_hybrid_schwarz_rows
use dg_hybrid_schwarz_solver,only:advance_dg_hybrid_schwarz_epoch,&
  assign_dg_hybrid_schwarz_occupations
use dg_hybrid_lcfo,only:assemble_dg_hybrid_lcfo_rows
use dg_nonlocal_projector_range,only:s_dg_nonlocal_range_receipt,analyze_dg_nonlocal_projector_range
use dg_hybrid_generalized_eigensystem,only:s_dg_hybrid_complete_eigensystem,&
  solve_dg_hybrid_generalized_scalapack,solve_dg_hybrid_generalized_once_and_publish,&
  solve_dg_hybrid_generalized_complete_once
use dg_hybrid_occupation_policy,only:s_dg_hybrid_occupation_result,derive_dg_hybrid_occupation_policy
use dg_hybrid_density,only:reconstruct_dg_hybrid_density,reconstruct_dg_hybrid_occupied_state
use dg_hybrid_ground_state_types,only:s_dg_hybrid_ground_state,s_dg_hybrid_spectral_certification,&
  validate_dg_hybrid_ground_state
use dg_hybrid_certified_rt_basis,only:s_dg_hybrid_certified_rt_basis,build_dg_hybrid_certified_rt_basis
use rt_dg_hybrid_checkpoint,only:write_rt_dg_hybrid_occupied_checkpoint,&
  s_rt_dg_hybrid_ground_state_payload,write_rt_dg_hybrid_ground_state_checkpoint,&
  fingerprint_rt_dg_hybrid_component,fingerprint_rt_dg_hybrid_ground_state_payload,&
  rt_dg_hybrid_ground_state_checkpoint_version,rt_dg_hybrid_energy_window_explicit,&
  rt_dg_hybrid_energy_window_legacy_dynamic,rt_dg_hybrid_vector_canonical_momentum
use dg_overlapping_wannier_solver, only: solve_dg_overlapping_wannier_coefficients
#ifdef USE_EIGENEXA
use dg_overlapping_wannier_solver, only: solve_dg_overlapping_wannier_generalized_eigenexa
#endif
use dg_overlapping_wannier_density,only:reconstruct_dg_overlapping_wannier_density
use dg_overlapping_wannier_checkpoint, only: s_dg_overlapping_wannier_checkpoint, &
  write_dg_overlapping_wannier_checkpoint,read_dg_overlapping_wannier_checkpoint,&
  compute_dg_overlapping_wannier_matrix_fingerprints
use dg_overlapping_wannier_observables, only: assemble_dg_overlapping_wannier_observables,&
  assemble_dg_cell_wrapped_position
use dg_overlapping_wannier_full_cell, only: project_dg_full_cell_hamiltonian_tiles,&
  s_dg_full_cell_redistribution_schedule,initialize_dg_full_cell_redistribution,&
  apply_dg_full_cell_redistribution_forward,apply_dg_full_cell_redistribution_reverse
use hamiltonian, only: hpsi
use nonlocal_potential,only:calc_uVpsi_rdivided
use sendrecv_grid, only: init_sendrecv_grid,dealloc_cache
use dg_overlapping_wannier_symmetry, only: select_dg_exact_fragment_subgroup,&
  build_dg_fragment_site_stabilizer,build_dg_fragment_group_representation,&
  promote_dg_exact_global_subgroup,project_dg_fragment_covariant_operators,&
  evaluate_dg_covariance_residuals_by_operation,fingerprint_dg_exact_fragment_symmetry
use dg_overlapping_wannier_symmetry, only: build_dg_fragment_permuted_representation,&
  build_dg_fragment_symmetry_orbits,factor_dg_affine_translation_cocycle,&
  symmetrize_dg_distributed_pencil_rows
use dg_overlapping_wannier_w90,only:setup_dg_w90_gamma_library,&
  assemble_dg_w90_gamma_a_matrix,assemble_dg_w90_gamma_matrices,run_dg_w90_gamma_library,&
  apply_dg_w90_gamma_transform,&
  inherit_dg_w90_affine_receipts,validate_dg_w90_generator_covariance,&
  project_dg_w90_reference_sector_operators,&
  anchor_dg_w90_reference_character_sector,align_dg_w90_character_sector_gauge,&
  align_dg_w90_character_sectors_by_periodic_phase,sew_dg_w90_periodic_phase_conjugate_sector,&
  export_dg_w90_replay_bundle,DG_W90_CONSTRAINED,DG_W90_UNCONSTRAINED
use lcfo_wannier_sawf, only: t_sawf_crystallographic_catalog,t_sawf_symop,&
  load_sawf_crystallographic_catalog_auto
use lcfo_wannier_sawf_dmn,only:t_sawf_dmn_writer,t_sawf_operation_index,&
  begin_sawf_dmn,append_sawf_dmn_operation,finish_sawf_dmn,abort_sawf_dmn,&
  build_sawf_operation_index,lookup_sawf_operation_product,&
  convert_sawf_pullback_to_active_representation
use lcfo_wannier_sawf_band, only: validate_sawf_fragment_symmetry_map,&
  build_sawf_fragment_buffer_point_map
#ifdef USE_EIGENEXA
use eigenexa_module, only: init_eigenexa_mod=>init_eigenexa,finalize_eigenexa
#endif
use parallelization, only: nproc_id_global,nproc_group_global,adjust_elapse_time,nproc_size_global
use communication, only: comm_is_root, comm_summation, comm_bcast, comm_sync_all, comm_get_max, comm_logical_and
use salmon_xc
use timer
use scf_iteration_sub
use density_matrix, only: calc_density
use writefield
use salmon_pp, only: calc_nlcc
use hartree_sub, only: hartree
use force_sub
use write_sub
use read_gs
use filesystem,only:atomic_create_directory
use code_optimization
use initialization_sub
use occupation
use prep_pp_sub
use mixing_sub
use checkpoint_restart_sub
use hamiltonian
use structure_opt_sub
use total_energy
use band_dft_sub
use init_gs, only: init_wf
use initialization_dft
use jellium, only: check_condition_jm
use dcdft
use dcdft_soi
use lcfo
use lcfo_flux
use lcfo_soi
implicit none
integer :: ix,iy,iz
integer :: Miter,iatom,jj,nspin
integer(8) :: dg_gs_potential_epoch
real(8) :: sum1
character(100) :: comment_line
character(1024) :: variational_payload_capture_prefix
integer :: variational_payload_capture_length,variational_payload_capture_status

type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_parallel_info) :: info
type(s_sendrecv_grid) :: srg, srg_scalar
type(s_orbital) :: spsi,shpsi,sttpsi
type(s_dft_system) :: system
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_xc_functional) :: xc_func
type(s_scalar) :: rho,rho_jm,Vh,Vpsl
type(s_scalar),allocatable :: V_local(:),rho_s(:),Vxc(:)
type(s_reciprocal_grid) :: fg
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_dft_energy) :: energy
type(s_ewald_ion_ion) :: ewald
type(s_cg)     :: cg
type(s_mixing) :: mixing
type(s_ofile)  :: ofl
type(s_band_dft) ::band
type(s_opt) :: opt
type(s_dcdft) :: dc
type(s_dg_dc_seed_contract) :: dg_dc_seed_contract
type(s_dg_dc_seed_payload) :: dg_dc_seed_payload

logical :: rion_update
logical :: flag_opt_conv
logical :: local_basis_route_active
logical :: dg_dc_seed_run_scf,dg_dc_seed_load,dg_dc_seed_publish,dg_dc_seed_fatal
logical :: dg_dc_seed_scf_skipped,dg_dc_seed_ok,dg_dc_seed_collective_ok
integer :: Miopt, iopt,nopt_max,i
integer :: dg_dc_seed_status
integer :: dg_dc_seed_rwf_bounds(14),dg_dc_seed_rho_bounds(6),dg_dc_seed_vloc_bounds(6)
integer(int64) :: dg_dc_seed_publication_id
integer(int64) :: dg_dc_seed_immutable_inputs(6),dg_dc_seed_ownership_map(8)
real(8) :: dg_dc_seed_electron_tolerance
character(512) :: dg_dc_seed_message
integer :: iter_band_kpt, iter_band_kpt_end, iter_band_kpt_stride
logical :: is_checkpoint_iter, is_shutdown_time
type(s_dg_overlapping_wannier_construction) :: ow_basis
type(s_dg_overlapping_wannier_scf_state) :: ow_state
type(s_dg_overlapping_wannier_scf_result) :: ow_result
type(s_dg_overlapping_wannier_checkpoint) :: ow_checkpoint
complex(8),allocatable :: ow_srows(:,:),ow_rhorows(:,:),ow_core_values(:,:),ow_core_gradients(:,:,:),&
  ow_box_values(:,:),ow_box_gradients(:,:,:),ow_last_kinetic_rows(:,:),&
  ow_last_local_rows(:,:),ow_last_nonlocal_rows(:,:),ow_published_hrows(:,:),ow_direct_nonlocal_rows(:,:)
integer(8),allocatable :: ow_core_ids(:),ow_row_ids(:)
integer(8),allocatable :: ow_box_physical_ids(:)
integer,allocatable :: ow_tail_generation(:,:)
integer,allocatable :: ow_core_box_positions(:)
real(8),allocatable :: ow_core_weights(:)
complex(8),allocatable :: dg_certified_localizer_values(:,:)
real(8),allocatable :: dg_certified_localizer_weights(:)
integer(8),allocatable :: dg_certified_localizer_grid_ids(:)
real(8),allocatable :: ow_partition_weight(:),ow_partition_gradient(:,:)
complex(8),allocatable :: ow_pencil_generator_representation(:,:,:)
integer,allocatable :: ow_pencil_generator_operations(:),ow_pencil_affine_product(:,:),&
  ow_pencil_translation_subgroup(:),ow_pencil_coset_representatives(:)
integer :: ow_box_size(3),ow_core_size(3),ow_buffer(3)
integer(8) :: ow_symmetry_fingerprint
integer(8) :: ow_potential_epoch_snapshot
integer(8) :: ow_global_grid_count
type(s_dg_hybrid_ground_state) :: ow_hybrid_ground_state
type(s_dg_hybrid_fragment_basis) :: divided_fragment_basis
type(s_dg_hybrid_divided_mixing_state) :: divided_mixing_state
type(s_dg_hybrid_schwarz_state) :: bounded_schwarz_state
type(s_dg_hybrid_schwarz_schedule) :: bounded_schwarz_schedule
complex(8),allocatable :: ow_hybrid_hrows(:,:),ow_hybrid_coefficients(:,:)
complex(8),allocatable :: divided_fragment_coefficients(:,:)
real(8),allocatable :: ow_hybrid_occupations(:),ow_hybrid_eigenvalues(:),ow_hybrid_potential(:),ow_hybrid_density(:),&
  ow_hybrid_density_history(:,:),ow_hybrid_new_history(:,:)
real(8),allocatable :: ow_hybrid_divided_total_density(:,:,:)
real(8),allocatable :: divided_fragment_eigenvalues(:),divided_fragment_density(:)
logical,allocatable :: ow_divided_core_mask(:)
character(16) :: ow_hybrid_divided_convergence
real(8) :: ow_hybrid_divided_threshold
integer(8) :: ow_hybrid_operator_fingerprint=0_8,ow_hybrid_metric_fingerprint=0_8
integer(8) :: divided_solver_workspace=0_8,divided_solver_fingerprint=0_8
real(8) :: divided_fragment_electron_count=0d0,divided_fragment_residual=huge(1d0),&
  divided_fragment_orthogonality=huge(1d0)
type(s_dg_hybrid_projection_factorization_receipt) :: bounded_projection_receipt
type(s_dg_hybrid_admission_report) :: bounded_admission_report
type(s_dg_hybrid_fixed_payload) :: bounded_fixed_payload
type(s_dg_hybrid_fragment_preconditioner) :: bounded_fragment_preconditioner
type(s_dg_hybrid_preconditioner_key) :: bounded_preconditioner_key
type(s_dg_hybrid_fragment_epoch_budget) :: bounded_epoch_budget
type(s_dg_hybrid_thermal_state),allocatable :: bounded_thermal_state
complex(8),allocatable :: bounded_interior_values(:,:),bounded_local_potential_rows(:,:),&
  bounded_fragment_h(:,:),bounded_fragment_s(:,:),bounded_reference_frame(:,:)
complex(8),allocatable :: bounded_schwarz_candidate_vectors(:,:)
real(8),allocatable :: bounded_schwarz_candidate_energies(:)
integer(8),allocatable :: bounded_schwarz_candidate_ids(:)
real(8),allocatable :: bounded_core_weights(:),bounded_point_weights(:)
integer,allocatable :: bounded_basis_fragment(:),bounded_basis_local_slot(:),&
  bounded_basis_generation(:),bounded_interior_fragment(:)
integer(8),allocatable :: bounded_core_ids(:)
integer(8) :: bounded_global_basis_fingerprint=0_8,bounded_directory_fingerprint=0_8,&
  bounded_frame_fingerprint=0_8,bounded_selection_fingerprint=0_8
integer(8) :: bounded_face_fingerprint=0_8,bounded_mapping_fingerprint=0_8,&
  bounded_candidate_fingerprint=0_8
integer :: bounded_last_peer_exchange_count=0
integer(8) :: divided_mixing_inventory_fingerprint=0_8
integer :: divided_mixing_basis_generation=0
integer :: ow_hybrid_history_count=0
real(8) :: ow_hybrid_mixing_rate=0d0
real(8) :: ow_hybrid_eigensystem_residual=huge(1d0),ow_hybrid_orthogonality=huge(1d0),&
  ow_hybrid_symmetry_defect=huge(1d0)
integer(8) :: ow_diag_h_local_bytes=0_8
integer :: ow_full_cell_component_mode=0
type(s_dg_full_cell_redistribution_schedule) :: ow_hpsi_redistribution
integer(8),allocatable :: ow_hpsi_grid_ids(:)
integer(8) :: ow_hpsi_redistribution_workspace=0_8
real(8),allocatable :: ow_density_snapshot(:,:,:,:)
real(8),allocatable :: ow_work_density(:,:,:,:)
real(8) :: ow_diag_t_hermiticity,ow_diag_vlocal_hermiticity,ow_diag_vnl_hermiticity,&
  ow_diag_h_hermiticity
logical :: ow_transaction_active
logical :: divided_mixing_rollback_pending=.false.
logical :: ow_direct_nonlocal_compared=.false.
logical :: ow_projector_stage_diagnosed=.false.
integer :: ilevel_print

interface
  subroutine build_dg_hybrid_retained_basis_representation(comm_arg,global_count_arg,core_ids_arg,&
      core_weights_arg,pencil_maps_arg,fragment_basis_arg,row_ids_arg,s_rows_arg,tolerance_arg,&
      representation_arg,closure_defect_arg,closure_ok_arg,callback_ok,callback_message)
    import::s_dg_hybrid_fragment_basis
    integer,intent(in)::comm_arg,global_count_arg
    integer(8),intent(in)::core_ids_arg(:),pencil_maps_arg(:,:),row_ids_arg(:)
    real(8),intent(in)::core_weights_arg(:),tolerance_arg
    complex(8),intent(in)::s_rows_arg(:,:)
    type(s_dg_hybrid_fragment_basis),intent(in)::fragment_basis_arg
    complex(8),allocatable,intent(out)::representation_arg(:,:,:)
    real(8),intent(out)::closure_defect_arg
    logical,intent(out)::closure_ok_arg,callback_ok
    character(*),intent(out)::callback_message
  end subroutine build_dg_hybrid_retained_basis_representation
end interface

dg_dc_seed_ok=.true.;dg_dc_seed_collective_ok=.true.
if(trim(dg_dc_seed_mode)/='off')then
  dg_dc_seed_ok=yn_dc=='y'.and.yn_dg_dc_overlapping_wannier=='y'.and.&
    trim(theory)=='dft'.and.iperiodic==3.and.yn_spinorbit=='n'.and.yn_opt/='y'.and.&
    .not.PLUS_U_ON.and.yn_hse/='y'.and.yn_fix_func/='y'.and.yn_jm/='y'
  call comm_logical_and(dg_dc_seed_ok,dg_dc_seed_collective_ok,nproc_group_global)
  if(.not.dg_dc_seed_collective_ok)&
    error stop 'DG DC seed mode is enabled outside its supported conventional-DC scope'
endif

if(theory=='dft_band'.and.iperiodic/=3) return

if(yn_dc=='y') then
  if(yn_spinorbit=='y') then
    call init_dcdft_soi(dc,pp,mixing,ewald)
  else
    call init_dcdft(dc,pp,mixing,ewald)
  end if
  ilevel_print = 0
else
  ilevel_print = 3
end if

!check condition for using jellium model
if(yn_jm=='y') call check_condition_jm

call init_xc(xc_func, spin, cval, xcname=xc, xname=xname, cname=cname)

call timer_begin(LOG_TOTAL)
call timer_begin(LOG_INIT_GS)


! please move folloings into initialization_dft
call init_dft(nproc_group_global,info,lg,mg,system,stencil,fg,poisson,srg,srg_scalar,ofl)
allocate( rho_s(system%nspin),V_local(system%nspin),Vxc(system%nspin) )

call initialization1_dft( system, energy, stencil, fg, poisson,  &
                          lg, mg,   &
                          info,  &
                          srg, srg_scalar,  &
                          rho, rho_jm, rho_s, Vh, V_local, Vpsl, Vxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,  &
                          ofl )

call initialization2_dft( Miter, nspin, rion_update,  &
                          system, energy, ewald, stencil, fg, poisson,&
                          lg, mg, info,   &
                          srg, srg_scalar,  &
                          rho, rho_jm, rho_s, Vh,V_local, Vpsl, Vxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,   &
                          xc_func, mixing )

dg_dc_seed_run_scf=.true.;dg_dc_seed_load=.false.;dg_dc_seed_publish=.false.
dg_dc_seed_fatal=.false.;dg_dc_seed_scf_skipped=.false.;dg_dc_seed_ok=.true.
dg_dc_seed_collective_ok=.true.
dg_dc_seed_status=DG_DC_SEED_ABSENT;dg_dc_seed_publication_id=0_int64
dg_dc_seed_message='';dg_dc_seed_contract=s_dg_dc_seed_contract()
dg_dc_seed_electron_tolerance=dg_dc_gs_electron_count_tolerance
if(trim(dg_dc_seed_mode)/='off')then
  dg_dc_seed_ok=.not.(yn_dc/='y'.or.yn_dg_dc_overlapping_wannier/='y'.or.&
    yn_spinorbit/='n'.or.system%nspin/=1.or.system%nk/=1.or.&
    .not.system%if_real_orbital.or..not.allocated(spsi%rwf).or.yn_opt=='y'.or.&
    theory=='dft_band'.or.PLUS_U_ON.or.yn_hse=='y'.or.yn_fix_func=='y'.or.yn_jm=='y')
  call comm_logical_and(dg_dc_seed_ok,dg_dc_seed_collective_ok,dc%icomm_tot)
  if(.not.dg_dc_seed_collective_ok)&
    error stop 'DG DC seed reuse requires supported one-shot real-Gamma overlapping-Wannier DC'
  call prepare_dg_dc_seed_contract_inputs(dg_dc_seed_rwf_bounds,dg_dc_seed_rho_bounds,&
    dg_dc_seed_vloc_bounds,dg_dc_seed_immutable_inputs,dg_dc_seed_ownership_map,&
    dg_dc_seed_ok,dg_dc_seed_message)
  call comm_logical_and(dg_dc_seed_ok,dg_dc_seed_collective_ok,dc%icomm_tot)
  if(.not.dg_dc_seed_collective_ok)error stop 'cannot prepare DG DC seed contract inputs'
  call build_dg_dc_seed_contract(dc%icomm_tot,dc%i_frag,dg_dc_seed_rwf_bounds,&
    dg_dc_seed_rho_bounds,dg_dc_seed_vloc_bounds,dg_dc_seed_immutable_inputs,&
    dg_dc_seed_ownership_map,dg_dc_seed_contract,dg_dc_seed_ok,dg_dc_seed_message)
  if(.not.dg_dc_seed_ok)error stop 'cannot build DG DC seed contract'
  dg_dc_seed_ok=dc%id_tot>=0.and.dc%id_tot<dc%isize_tot.and.&
    dc%id_tot==dg_dc_seed_contract%rank.and.&
    dc%isize_tot==dg_dc_seed_contract%mpi_size
  call comm_logical_and(dg_dc_seed_ok,dg_dc_seed_collective_ok,dc%icomm_tot)
  if(.not.dg_dc_seed_collective_ok)&
    error stop 'invalid preserved total-system topology for DG DC seed'
endif

select case(trim(dg_dc_seed_mode))
case('off')
  call resolve_dg_dc_seed_mode(dg_dc_seed_mode,DG_DC_SEED_ABSENT,&
    dg_dc_seed_run_scf,dg_dc_seed_load,dg_dc_seed_publish,dg_dc_seed_fatal,&
    dg_dc_seed_scf_skipped,dg_dc_seed_ok,dg_dc_seed_message)
case('write')
  call atomic_create_directory(trim(dg_dc_seed_directory),dc%icomm_tot,dc%id_tot)
  call resolve_dg_dc_seed_mode(dg_dc_seed_mode,DG_DC_SEED_ABSENT,&
    dg_dc_seed_run_scf,dg_dc_seed_load,dg_dc_seed_publish,dg_dc_seed_fatal,&
    dg_dc_seed_scf_skipped,dg_dc_seed_ok,dg_dc_seed_message)
case('read')
  call probe_dg_dc_seed(dc%icomm_tot,trim(dg_dc_seed_directory),dg_dc_seed_contract,&
    dc%system_tot%hvol,dc%elec_num_tot,dg_dc_seed_electron_tolerance,threshold,&
    dg_dc_seed_status,dg_dc_seed_publication_id,dg_dc_seed_message)
  call resolve_dg_dc_seed_mode(dg_dc_seed_mode,dg_dc_seed_status,&
    dg_dc_seed_run_scf,dg_dc_seed_load,dg_dc_seed_publish,dg_dc_seed_fatal,&
    dg_dc_seed_scf_skipped,dg_dc_seed_ok,dg_dc_seed_message)
case('auto')
  call atomic_create_directory(trim(dg_dc_seed_directory),dc%icomm_tot,dc%id_tot)
  call probe_dg_dc_seed(dc%icomm_tot,trim(dg_dc_seed_directory),dg_dc_seed_contract,&
    dc%system_tot%hvol,dc%elec_num_tot,dg_dc_seed_electron_tolerance,threshold,&
    dg_dc_seed_status,dg_dc_seed_publication_id,dg_dc_seed_message)
  call resolve_dg_dc_seed_mode(dg_dc_seed_mode,dg_dc_seed_status,&
    dg_dc_seed_run_scf,dg_dc_seed_load,dg_dc_seed_publish,dg_dc_seed_fatal,&
    dg_dc_seed_scf_skipped,dg_dc_seed_ok,dg_dc_seed_message)
case default
  dg_dc_seed_fatal=.true.;dg_dc_seed_message='unknown DG DC seed mode'
end select
if(dg_dc_seed_fatal)error stop 'DG DC seed is absent, invalid, or incompatible'
if(dg_dc_seed_load)then
  call read_dg_dc_seed(dc%icomm_tot,trim(dg_dc_seed_directory),dg_dc_seed_contract,&
    dc%system_tot%hvol,dc%elec_num_tot,dg_dc_seed_electron_tolerance,threshold,&
    dg_dc_seed_payload,dg_dc_seed_publication_id,dg_dc_seed_ok,dg_dc_seed_message)
  if(.not.dg_dc_seed_ok)error stop 'failed to read compatible DG DC seed'
  call restore_dg_dc_seed_payload(dg_dc_seed_payload,spsi%rwf,dc%rho_tot_s(1)%f,&
    dc%vloc_tot(1)%f,energy%esp,system%rocc,system%mu,sum1,Miter,&
    dg_dc_seed_ok,dg_dc_seed_message)
  call comm_logical_and(dg_dc_seed_ok,dg_dc_seed_collective_ok,dc%icomm_tot)
  if(.not.dg_dc_seed_collective_ok)error stop 'failed to restore compatible DG DC seed'
  call validate_dg_dc_seed_state(dg_dc_seed_ok,dg_dc_seed_message)
  if(.not.dg_dc_seed_ok)error stop 'restored DG DC seed state is inconsistent'
  call rebuild_dg_dc_seed_derived_state_dcdft(mg,info,system,spsi,rho,rho_s,&
    V_local,dc,dg_dc_seed_ok,dg_dc_seed_message)
  if(.not.dg_dc_seed_ok)error stop 'failed to rebuild DG DC derived state'
endif

Miopt = 0
nopt_max = 1
if(yn_opt=='y') call initialization_opt(Miopt,opt,system,flag_opt_conv,nopt_max,ofl)

call timer_end(LOG_INIT_GS)

if(yn_dc == 'y' .and. yn_dc_lcfo_wannier == 'y' .and. dc_lcfo_wannier_import_only_requested()) then
  if(comm_is_root(nproc_id_global)) &
    write(*,'(1x,a)') '[DC-LCFO-W90-IMPORT] import-only mode: skip SCF and reuse external Wannier90 outputs'
  call dc_lcfo_wannier_import_only(dc)
  call timer_end(LOG_TOTAL)
  return
end if

!---------------------------------------- Opt Iteration


#ifdef __FUJITSU
call fipp_start ! performance profiling
#endif

Structure_Optimization_Iteration : do iopt= Miopt+1, nopt_max

if(iopt>=2)then
  call timer_begin(LOG_INIT_GS)
  Miter = 0        ! Miter: Iteration counter set to zero
  rion_update = .true.
  call dealloc_init_ps(ppg)
  call init_ps(lg,mg,system,info,fg,poisson,pp,ppg,Vpsl)
  call calc_nlcc(pp, system, mg, ppn)
  if(yn_auto_mixing=='y') call reset_mixing_rate(mixing)
  call timer_end(LOG_INIT_GS)
end if

!---------------------------------------- Band Iteration

if(theory=='dft_band')then
   call init_band_dft(system,band) ! --> system%wtk=0.0
   iter_band_kpt_end    = band%num_band_kpt
   iter_band_kpt_stride = system%nk
else
   iter_band_kpt_end    = 1
   iter_band_kpt_stride = 1
end if

call comm_sync_all
call timer_enable_sub
Band_Iteration : do iter_band_kpt= 1, iter_band_kpt_end, iter_band_kpt_stride

if(theory=='dft_band')then
   call calc_band_write(iter_band_kpt,system,band,info)
end if


call timer_begin(LOG_INIT_GS_ITERATION)

call timer_end(LOG_INIT_GS_ITERATION)


call timer_begin(LOG_GS_ITERATION)
!------------------------------------ SCF Iteration
!Iteration loop for SCF (DFT_Iteration)
if(dg_dc_seed_run_scf) then
call scf_iteration_dft( Miter,rion_update,sum1,  &
                        system,energy,ewald,  &
                        lg,mg,  &
                        info,  &
                        poisson,fg,  &
                        cg,mixing,  &
                        stencil,  &
                        srg,srg_scalar,   &
                        spsi,shpsi,sttpsi,  &
                        rho,rho_jm,rho_s,  &
                        V_local,Vh,Vxc,Vpsl,xc_func,  &
                        pp,ppg,ppn,  &
                        band, ilevel_print, dc)
endif


if(theory=='dft_band')then
   call write_band(system,energy)
end if

! output the wavefunctions for next GS calculations
if(write_gs_wfn_k == 'y') then !this input keyword is going to be removed....
   select case(iperiodic)
   case(3)
      call write_wfn(lg,mg,spsi,info,system)
      ! Experimental Implementation of Inner-Product Outputs:
      ! call write_prod_dk_data(lg, mg, system, info, spsi)
   case(0)
      write(*,*) "error: write_gs_wfn_k='y' & iperiodic=0"
   end select
end if

! output transition moment : --> want to put out of the optmization loop in future
if(yn_out_tm  == 'y'.or.yn_out_gs_sgm_eps=='y') then
   select case(iperiodic)
   case(3)
      call write_k_data(system,stencil)  !need? (probably remove later)
      call write_tm_data(spsi,system,info,mg,stencil,srg,ppg,energy)
   case(0)
     write(*,*) "error: yn_out_tm='y',yn_out_gs_sgm_eps='y' & iperiodic=0"
  end select
end if

   ! force
   if(yn_jm=='n' .and. yn_dc=="n")then
     call calc_force(system,pp,fg,info,mg,stencil,poisson,srg,ppg,spsi,ewald)
     if(comm_is_root(nproc_id_global))then
        write(*,*) "===== force ====="
        do iatom=1,natom
           select case(unit_system)
           case('au','a.u.'); write(*,300)iatom,(system%Force(ix,iatom),ix=1,3)
           case('A_eV_fs'  ); write(*,300)iatom,(system%Force(ix,iatom)*au_energy_ev/au_length_aa,ix=1,3)
           end select
        end do
300   format(i6,3e16.8)
     end if
   end if

call timer_end(LOG_GS_ITERATION)

end do Band_Iteration
call timer_disable_sub


call timer_begin(LOG_DEINIT_GS_ITERATION)
if(yn_opt=='y') then
   call structure_opt_check(iopt,flag_opt_conv,system%Force)
   if(.not.flag_opt_conv) call structure_opt(opt,iopt,system)
   !! Rion is old variables to be removed
   !! but currently it is used in many subroutines.
   !!Rion(:,:) = system%Rion(:,:)

   write(comment_line,10) iopt
   call write_xyz(comment_line,"add","r  ",system,ofl)
10 format("#opt iteration step=",i5)

   if(comm_is_root(nproc_id_global))then
      write(*,*) "atomic coordinate"
      do iatom=1,natom
         write(*,20) "'"//trim(atom_name(iatom))//"'",  &
                   (system%Rion(jj,iatom)*ulength_from_au,jj=1,3), &
                   Kion(iatom), flag_opt_atom(iatom)
      end do
20    format(a5,3f16.8,i3,a3)
   end if

   if(flag_opt_conv) then
      call structure_opt_fin(opt)
   else
      is_checkpoint_iter = (checkpoint_interval >= 1) .and. (mod(iopt,checkpoint_interval) == 0)
      is_shutdown_time   = (time_shutdown > 0d0) .and. (adjust_elapse_time(timer_now(LOG_TOTAL)) > time_shutdown)

      if(yn_dg_dc_overlapping_wannier/='y' .and. &
         (is_checkpoint_iter .or. is_shutdown_time)) then
         if (is_shutdown_time .and. comm_is_root(info%id_rko)) then
           print *, 'shutdown the calculation, iopt =', iopt
         end if

         call checkpoint_gs(lg,mg,system,info,spsi,iopt,mixing)
         call comm_sync_all
         call checkpoint_opt(iopt,opt)
         if(comm_is_root(nproc_id_global))then
            write(*,'(a,i5)')"  checkpoint data is printed: iopt=", iopt
         endif
         call comm_sync_all

         if (is_shutdown_time) then
           exit Structure_Optimization_Iteration
         end if
      endif
   endif

end if
call timer_end(LOG_DEINIT_GS_ITERATION)


if(yn_opt=='y')then
  if(flag_opt_conv)then
  exit Structure_Optimization_Iteration
  end if
end if
end do Structure_Optimization_Iteration

#ifdef __FUJITSU
call fipp_stop ! performance profiling
#endif


!------------ Writing part -----------
call timer_begin(LOG_WRITE_GS_RESULTS)
local_basis_route_active=.false.

if(yn_dc=='y') then
  if(yn_dg_dc_overlapping_wannier == 'y') then
    local_basis_route_active=.true.
    if(.not.(sum1<threshold))&
      error stop 'overlapping-Wannier route requires a converged conventional DC state'
    if(dg_dc_seed_publish)then
      call validate_dg_dc_seed_state(dg_dc_seed_ok,dg_dc_seed_message)
      if(.not.dg_dc_seed_ok)error stop 'converged DG DC seed state is inconsistent'
      call capture_dg_dc_seed_payload_dcdft(system,energy,spsi,dc,sum1,Miter,&
        dg_dc_seed_payload,dg_dc_seed_ok,dg_dc_seed_message)
      call comm_logical_and(dg_dc_seed_ok,dg_dc_seed_collective_ok,dc%icomm_tot)
      if(.not.dg_dc_seed_collective_ok)error stop 'failed to capture converged DG DC seed state'
      call write_dg_dc_seed(dc%icomm_tot,trim(dg_dc_seed_directory),dg_dc_seed_contract,&
        dg_dc_seed_payload,dc%system_tot%hvol,dc%elec_num_tot,&
        dg_dc_seed_electron_tolerance,threshold,dg_dc_seed_publication_id,&
        dg_dc_seed_ok,dg_dc_seed_message)
      if(.not.dg_dc_seed_ok)error stop 'failed to publish converged DG DC seed state'
    endif
    if(dc%id_tot==0)write(*,'(a,a,a,i0,a,l1,a,i0,a,i0)')&
      '[DG-DC-SEED] mode=',trim(dg_dc_seed_mode),&
      ' publication_id=',dg_dc_seed_publication_id,&
      ' scf_skipped=',dg_dc_seed_scf_skipped,&
      ' mpi_size=',dc%isize_tot,&
      ' mapping_fingerprint=',dg_dc_seed_contract%ownership_fingerprint
    if(yn_dg_hybrid_divided_scf == 'y') then
      call run_dg_hybrid_divided_ground_state_for_main
    else if(yn_dg_hybrid_continuation_scf == 'y') then
      call run_dg_hybrid_continuation_ground_state_for_main
    else
      call run_dg_overlapping_wannier_ground_state_for_main
    end if
    return
  else if(yn_dc_lcfo_flux == 'y') then
    if(yn_spinorbit == 'y') then
      stop "yn_dc_lcfo_flux=y is not implemented for spin-orbit mode"
    else
      if(comm_is_root(nproc_id_global)) &
      & write(*,'(1x,a)') '[DC-LCFO-FLUX] export phase: build Flux-LCFO basis and coefficients'
#ifdef USE_EIGENEXA
      call finalize_eigenexa(info)
#endif
      call dc_lcfo_flux(lg,mg,system,info,stencil,ppg,energy,rho_s,v_local,&
        spsi,shpsi,sttpsi,srg,dc)
    end if
  else if(yn_dc_lcfo == 'y') then
    if(yn_spinorbit == 'y') then
      call dc_lcfo_soi(lg,mg,system,info,stencil,ppg,energy,v_local,spsi,shpsi,sttpsi,srg,dc)
    else
      call dc_lcfo(lg,mg,system,info,stencil,ppg,energy,v_local,spsi,shpsi,sttpsi,srg,dc)
    end if
  end if
  if(.not.local_basis_route_active .and. yn_spinorbit == 'y') then
    call write_total_dcdft_soi(system,dc)
  else if(.not.local_basis_route_active) then
    call write_total_dcdft(system,dc)
  end if
end if

! write GS: basic data
if(.not.local_basis_route_active) then
if(yn_dc=='n') call write_band_information(system,energy)
call write_eigen(ofl,system,energy)
call write_info_data(Miter,system,energy,pp)
call write_k_data(system,stencil)
if(yn_spinorbit=='y') call write_mag_decomposed_gs(system,mg,info,spsi)

! write GS: analysis option
if(yn_out_psi =='y') call write_psi(lg,mg,system,info,spsi)
if(yn_out_dns =='y') call write_dns(lg,mg,system,info,rho_s)
if(yn_out_dos =='y') call write_dos(system,energy)
if(yn_out_pdos=='y') call write_pdos(lg,mg,system,info,pp,energy,spsi)
if(yn_out_elf =='y') call write_elf(0,lg,mg,system,info,stencil,rho,srg,srg_scalar,spsi)
end if

call timer_end(LOG_WRITE_GS_RESULTS)

! write GS: binary data for restart
call timer_begin(LOG_WRITE_GS_DATA)
if(.not.local_basis_route_active) then
if(write_gs_restart_data=="no") then
   if(comm_is_root(nproc_id_global)) &
      write(*,'(a)')"  no restart data writing."
else if(write_gs_restart_data.ne."checkpoint_only") then
   if(comm_is_root(nproc_id_global)) write(*,'(a)')"  writing restart data..."
   call checkpoint_gs(lg,mg,system,info,spsi,Miter,mixing,ofl%dir_out_restart)
   call comm_sync_all
   if(yn_opt=='y') then
      if(.not.flag_opt_conv) then
         call comm_sync_all
         call checkpoint_opt(nopt_max,opt,ofl%dir_out_restart)
         call comm_sync_all
      endif
   endif
else
   if(yn_self_checkpoint=='n') then
      if(comm_is_root(nproc_id_global)) then
           write(*,'(a)')"  no restart data writing:"
           write(*,'(a)')"  check input keywords if you need restart data"
       endif
   endif
endif
if(yn_self_checkpoint=='y') then
   if(comm_is_root(nproc_id_global)) &
   write(*,'(a)')"  writing restart data in checkpoint format ..."
   call checkpoint_gs(lg,mg,system,info,spsi,Miter,mixing)
   call comm_sync_all
endif
if(comm_is_root(nproc_id_global)) write(*,'(a)')"  writing completed."
else if(comm_is_root(nproc_id_global)) then
  write(*,'(a)')'  DG local-basis result retained in memory; standard restart publication skipped.'
end if
call timer_end(LOG_WRITE_GS_DATA)

!call timer_begin(LOG_WRITE_GS_INFO)  !if needed, please take back, sory: AY
!call timer_end(LOG_WRITE_GS_INFO)

call finalize_xc(xc_func)

if(yn_dc=='y') then
! override (restore)
  nproc_group_global = dc%icomm_tot
  nproc_id_global = dc%id_tot
  nproc_size_global = dc%isize_tot
  call comm_sync_all
  call finalize_dcdft(dc)
end if

call timer_end(LOG_TOTAL)

contains

  subroutine prepare_dg_dc_seed_contract_inputs(rwf_bounds,rho_bounds,vloc_bounds,&
      immutable_inputs,ownership_map,ok,message)
    integer,intent(out)::rwf_bounds(14),rho_bounds(6),vloc_bounds(6)
    integer(int64),intent(out)::immutable_inputs(6),ownership_map(8)
    logical,intent(out)::ok
    character(*),intent(out)::message

    rwf_bounds=0;rho_bounds=0;vloc_bounds=0
    immutable_inputs=0_int64;ownership_map=0_int64;ok=.false.;message=''
    if(.not.allocated(spsi%rwf).or..not.allocated(dc%rho_tot_s).or.&
       .not.allocated(dc%vloc_tot).or..not.allocated(dc%rho_tot_s(1)%f).or.&
       .not.allocated(dc%vloc_tot(1)%f).or..not.allocated(dc%system_tot%Rion).or.&
       .not.allocated(dc%system_tot%kion).or..not.allocated(system%Rion).or.&
       .not.allocated(system%kion).or..not.allocated(dc%nxyz_domain_frag).or.&
       .not.allocated(dc%ixyz_frag).or..not.allocated(dc%rxyz_frag).or.&
       .not.allocated(dc%jxyz_tot))then
      message='DG DC seed contract arrays are not allocated';return
    endif
    rwf_bounds=[lbound(spsi%rwf),ubound(spsi%rwf)]
    rho_bounds=[lbound(dc%rho_tot_s(1)%f),ubound(dc%rho_tot_s(1)%f)]
    vloc_bounds=[lbound(dc%vloc_tot(1)%f),ubound(dc%vloc_tot(1)%f)]
    immutable_inputs(1)=int(z'4443445345454431',int64)
    immutable_inputs(2)=dg_dc_seed_cell_atom_fingerprint()
    immutable_inputs(3)=dg_dc_seed_fragment_topology_fingerprint()
    immutable_inputs(4)=dg_dc_seed_physics_fingerprint()
    immutable_inputs(5)=dg_dc_seed_operator_input_fingerprint()
    immutable_inputs(6)=dg_dc_seed_convergence_fingerprint()
    ownership_map(1)=int(z'4F574E4552534831',int64)
    ownership_map(2)=dg_dc_seed_grid_ownership_fingerprint(dc%mg_tot,dc%info_tot)
    ownership_map(3)=dg_dc_seed_grid_ownership_fingerprint(mg,info)
    ownership_map(4)=dg_dc_seed_orbital_ownership_fingerprint(info)
    ownership_map(5)=dg_dc_seed_fragment_map_fingerprint()
    ownership_map(6)=dg_dc_seed_ppg_ownership_fingerprint(dc%ppg_tot)
    ownership_map(7)=dg_dc_seed_ppg_ownership_fingerprint(ppg)
    ownership_map(8)=int(dc%id_tot+1,int64)
    if(any(immutable_inputs==0_int64).or.any(ownership_map==0_int64))then
      message='DG DC seed contract fingerprint is zero';return
    endif
    ok=.true.
  end subroutine prepare_dg_dc_seed_contract_inputs

  subroutine validate_dg_dc_seed_state(ok,message)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::local_bad,global_bad,ierr,iteration_min,iteration_max
    real(8)::common_values(2),minimum_values(2),maximum_values(2)
    integer::expected_state_bounds(6)

    ok=.false.;message='';local_bad=0
    expected_state_bounds=[1,1,1,system%no,system%nk,system%nspin]
    if(.not.allocated(spsi%rwf).or..not.allocated(dc%rho_tot_s(1)%f).or.&
       .not.allocated(dc%vloc_tot(1)%f).or..not.allocated(energy%esp).or.&
       .not.allocated(system%rocc))local_bad=1
    if(local_bad==0)then
      if(any([lbound(spsi%rwf),ubound(spsi%rwf)]/=dg_dc_seed_rwf_bounds).or.&
         any([lbound(dc%rho_tot_s(1)%f),ubound(dc%rho_tot_s(1)%f)]/=&
           dg_dc_seed_rho_bounds).or.&
         any([lbound(dc%vloc_tot(1)%f),ubound(dc%vloc_tot(1)%f)]/=&
           dg_dc_seed_vloc_bounds).or.&
         any([lbound(energy%esp),ubound(energy%esp)]/=expected_state_bounds).or.&
         any([lbound(system%rocc),ubound(system%rocc)]/=expected_state_bounds))local_bad=1
      if(.not.all(ieee_is_finite(energy%esp)).or.&
         .not.all(ieee_is_finite(system%rocc)).or.any(system%rocc<0d0).or.&
         any(system%rocc>2d0+100d0*epsilon(1d0)).or.&
         .not.ieee_is_finite(system%mu).or..not.ieee_is_finite(sum1).or.&
         sum1<0d0.or..not.(sum1<threshold).or.Miter<0)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid restored DG DC seed bounds or values';return
    endif
    common_values=[system%mu,sum1]
    call MPI_Allreduce(common_values,minimum_values,2,MPI_DOUBLE_PRECISION,MPI_MIN,&
      dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG DC seed scalar minimum failed';return;endif
    call MPI_Allreduce(common_values,maximum_values,2,MPI_DOUBLE_PRECISION,MPI_MAX,&
      dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG DC seed scalar maximum failed';return;endif
    call MPI_Allreduce(Miter,iteration_min,1,MPI_INTEGER,MPI_MIN,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG DC seed iteration minimum failed';return;endif
    call MPI_Allreduce(Miter,iteration_max,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum_values/=maximum_values).or.&
       iteration_min/=iteration_max)then
      message='rank-inconsistent DG DC seed scalar provenance';return
    endif
    ok=.true.
#else
    ok=.false.;message='DG DC seed restore requires MPI'
#endif
  end subroutine validate_dg_dc_seed_state

  integer(int64) function dg_dc_seed_cell_atom_fingerprint()result(hash)
    integer::ii,jj
    hash=int(z'6A09E667F3BCC909',int64)
    call hash_integer(hash,dc%lg_tot%num(1));call hash_integer(hash,dc%lg_tot%num(2))
    call hash_integer(hash,dc%lg_tot%num(3));call hash_real(hash,dc%system_tot%hvol)
    do jj=1,3
      call hash_real(hash,dc%system_tot%hgs(jj))
      do ii=1,3
        call hash_real(hash,dc%system_tot%primitive_a(ii,jj))
        call hash_real(hash,dc%system_tot%primitive_b(ii,jj))
        call hash_real(hash,dc%system_tot%rmatrix_a(ii,jj))
        call hash_real(hash,dc%system_tot%rmatrix_b(ii,jj))
      enddo
    enddo
    call hash_real(hash,dc%system_tot%det_a)
    call hash_integer(hash,dc%system_tot%nion)
    do jj=1,dc%system_tot%nion
      call hash_integer(hash,dc%system_tot%kion(jj))
      do ii=1,3;call hash_real(hash,dc%system_tot%Rion(ii,jj));enddo
    enddo
    call hash_integer(hash,system%nion)
    do jj=1,system%nion
      call hash_integer(hash,system%kion(jj))
      do ii=1,3;call hash_real(hash,system%Rion(ii,jj));enddo
    enddo
    if(hash==0_int64)hash=1_int64
  end function dg_dc_seed_cell_atom_fingerprint

  integer(int64) function dg_dc_seed_fragment_topology_fingerprint()result(hash)
    integer::axis,fragment
    hash=int(z'BB67AE8584CAA73B',int64)
    call hash_integer(hash,dc%n_frag);call hash_integer(hash,dc%i_frag)
    call hash_integer(hash,merge(1,0,dc%optimized_fragment_geometry))
    do axis=1,3
      call hash_integer(hash,dc%nxyz_domain(axis))
      call hash_integer(hash,dc%nxyz_buffer(axis))
    enddo
    do fragment=1,dc%n_frag;do axis=1,3
      call hash_integer(hash,dc%nxyz_domain_frag(axis,fragment))
      call hash_integer(hash,dc%ixyz_frag(axis,fragment))
      call hash_real(hash,dc%rxyz_frag(axis,fragment))
    enddo;enddo
    if(hash==0_int64)hash=1_int64
  end function dg_dc_seed_fragment_topology_fingerprint

  integer(int64) function dg_dc_seed_physics_fingerprint()result(hash)
    integer::i,j
    hash=int(z'3C6EF372FE94F82B',int64)
    call hash_integer(hash,dc%nstate_tot);call hash_integer(hash,dc%nstate_frag)
    call hash_real(hash,dc%elec_num_tot)
    call hash_integer(hash,dc%system_tot%nspin);call hash_integer(hash,dc%system_tot%no)
    call hash_integer(hash,dc%system_tot%nk)
    call hash_integer(hash,merge(1,0,dc%system_tot%if_real_orbital))
    call hash_integer(hash,system%nspin);call hash_integer(hash,system%no)
    call hash_integer(hash,system%nk);call hash_integer(hash,merge(1,0,system%if_real_orbital))
    call hash_real(hash,temperature)
    if(allocated(dc%system_tot%vec_k))then
      do j=1,size(dc%system_tot%vec_k,2);do i=1,size(dc%system_tot%vec_k,1)
        call hash_real(hash,dc%system_tot%vec_k(i,j))
      enddo;enddo
    endif
    if(allocated(dc%system_tot%wtk))then
      do i=1,size(dc%system_tot%wtk);call hash_real(hash,dc%system_tot%wtk(i));enddo
    endif
    call hash_character(hash,trim(calc_mode));call hash_character(hash,trim(theory))
    call hash_character(hash,trim(yn_spinorbit))
    if(hash==0_int64)hash=1_int64
  end function dg_dc_seed_physics_fingerprint

  integer(int64) function dg_dc_seed_operator_input_fingerprint()result(hash)
    integer::i,j
    hash=int(z'A54FF53A5F1D36F1',int64)
    call hash_integer(hash,merge(1,0,stencil%if_orthogonal))
    call hash_real(hash,stencil%coef_lap0);call hash_real(hash,stencil%coef_lap0_nd1)
    do j=1,3
      do i=1,4
        call hash_real(hash,stencil%coef_lap(i,j));call hash_real(hash,stencil%coef_nab(i,j))
      enddo
      call hash_real(hash,stencil%coef_lap_nd1(1,j));call hash_real(hash,stencil%coef_nab_nd1(1,j))
    enddo
    do i=1,6;call hash_real(hash,stencil%coef_f(i));enddo
    call hash_character(hash,trim(xc));call hash_character(hash,trim(xname))
    call hash_character(hash,trim(cname));call hash_character(hash,trim(alibxc))
    do i=1,3;call hash_integer(hash,xc_func%xctype(i));enddo
    call hash_integer(hash,xc_func%ispin);call hash_real(hash,xc_func%cval)
    call hash_integer(hash,merge(1,0,xc_func%use_gradient))
    call hash_integer(hash,merge(1,0,xc_func%use_laplacian))
    call hash_integer(hash,merge(1,0,xc_func%use_kinetic_energy))
    call hash_integer(hash,merge(1,0,xc_func%use_current))
    call hash_pp_info_for_dg_dc_seed(hash)
    if(hash==0_int64)hash=1_int64
  end function dg_dc_seed_operator_input_fingerprint

  integer(int64) function dg_dc_seed_convergence_fingerprint()result(hash)
    hash=int(z'510E527FADE682D1',int64)
    call hash_character(hash,trim(convergence));call hash_real(hash,threshold)
    call hash_character(hash,trim(method_mixing));call hash_real(hash,mixing%mixrate)
    call hash_real(hash,mixing%alpha_mb);call hash_real(hash,mixing%beta_p)
    call hash_integer(hash,nscf);call hash_integer(hash,nscf_init_redistribution)
    call hash_integer(hash,nscf_init_no_diagonal);call hash_integer(hash,nscf_init_mix_zero)
    if(hash==0_int64)hash=1_int64
  end function dg_dc_seed_convergence_fingerprint

  integer(int64) function dg_dc_seed_grid_ownership_fingerprint(grid,parallel)result(hash)
    type(s_rgrid),intent(in)::grid
    type(s_parallel_info),intent(in)::parallel
    integer::i
    hash=int(z'9B05688C2B3E6C1F',int64)
    do i=1,3
      call hash_integer(hash,grid%is(i));call hash_integer(hash,grid%ie(i))
      call hash_integer(hash,grid%num(i));call hash_integer(hash,grid%is_array(i))
      call hash_integer(hash,grid%ie_array(i));call hash_integer(hash,parallel%nprgrid(i))
      call hash_integer(hash,parallel%iaddress(i))
    enddo
    call hash_alloc_integer_rank2(hash,grid%is_all)
    call hash_alloc_integer_rank2(hash,grid%ie_all)
    call hash_alloc_integer_rank1(hash,grid%idx)
    call hash_alloc_integer_rank1(hash,grid%idy)
    call hash_alloc_integer_rank1(hash,grid%idz)
    if(hash==0_int64)hash=1_int64
  end function dg_dc_seed_grid_ownership_fingerprint

  integer(int64) function dg_dc_seed_orbital_ownership_fingerprint(parallel)result(hash)
    type(s_parallel_info),intent(in)::parallel
    integer::i
    hash=int(z'1F83D9ABFB41BD6B',int64)
    call hash_integer(hash,parallel%npk);call hash_integer(hash,parallel%nporbital)
    do i=1,5;call hash_integer(hash,parallel%iaddress(i));enddo
    call hash_integer(hash,parallel%im_s);call hash_integer(hash,parallel%im_e)
    call hash_integer(hash,parallel%numm);call hash_integer(hash,parallel%ik_s)
    call hash_integer(hash,parallel%ik_e);call hash_integer(hash,parallel%numk)
    call hash_integer(hash,parallel%io_s);call hash_integer(hash,parallel%io_e)
    call hash_integer(hash,parallel%numo)
    call hash_alloc_integer_rank5(hash,parallel%imap)
    call hash_alloc_integer_rank1(hash,parallel%irank_io)
    call hash_alloc_integer_rank1(hash,parallel%io_s_all)
    call hash_alloc_integer_rank1(hash,parallel%io_e_all)
    call hash_alloc_integer_rank1(hash,parallel%numo_all)
    if(hash==0_int64)hash=1_int64
  end function dg_dc_seed_orbital_ownership_fingerprint

  integer(int64) function dg_dc_seed_fragment_map_fingerprint()result(hash)
    hash=int(z'5BE0CD19137E2179',int64)
    call hash_integer(hash,dc%id_tot);call hash_integer(hash,dc%isize_tot)
    call hash_integer(hash,dc%i_frag);call hash_integer(hash,dc%id_frag)
    call hash_integer(hash,dc%isize_frag)
    call hash_alloc_integer_rank2(hash,dc%jxyz_tot)
    if(hash==0_int64)hash=1_int64
  end function dg_dc_seed_fragment_map_fingerprint

  integer(int64) function dg_dc_seed_ppg_ownership_fingerprint(grid)result(hash)
    type(s_pp_grid),intent(in)::grid
    hash=int(z'CBBB9D5DC1059ED8',int64)
    call hash_integer(hash,grid%nps);call hash_integer(hash,grid%nlma)
    call hash_integer(hash,grid%ilocal_nlma)
    call hash_alloc_integer_rank1(hash,grid%mps)
    call hash_alloc_integer_rank3(hash,grid%jxyz)
    call hash_alloc_integer_rank2(hash,grid%lma_tbl)
    call hash_alloc_integer_rank1(hash,grid%ia_tbl)
    call hash_alloc_integer_rank2(hash,grid%irange_atom)
    call hash_alloc_integer_rank1(hash,grid%ilocal_nlma2ilma)
    call hash_alloc_integer_rank1(hash,grid%ilocal_nlma2ia)
    call hash_alloc_integer_rank2(hash,grid%jxyz_min)
    call hash_alloc_integer_rank2(hash,grid%jxyz_max)
    if(hash==0_int64)hash=1_int64
  end function dg_dc_seed_ppg_ownership_fingerprint

  subroutine hash_pp_info_for_dg_dc_seed(hash)
    integer(int64),intent(inout)::hash
    call hash_real(hash,pp%zion);call hash_integer(hash,pp%lmax)
    call hash_integer(hash,pp%lmax0);call hash_integer(hash,pp%nrmax)
    call hash_integer(hash,pp%nrmax0);call hash_integer(hash,merge(1,0,pp%flag_nlcc))
    call hash_alloc_character_rank1(hash,pp%atom_symbol)
    call hash_alloc_real_rank1(hash,pp%rmass)
    call hash_alloc_integer_rank1(hash,pp%mr)
    call hash_alloc_integer_rank1(hash,pp%lref)
    call hash_alloc_integer_rank1(hash,pp%nrps)
    call hash_alloc_integer_rank1(hash,pp%mlps)
    call hash_alloc_integer_rank2(hash,pp%nproj)
    call hash_alloc_integer_rank1(hash,pp%num_orb)
    call hash_alloc_integer_rank1(hash,pp%zps)
    call hash_alloc_integer_rank1(hash,pp%nrloc)
    call hash_alloc_real_rank1(hash,pp%rloc)
    call hash_alloc_real_rank1(hash,pp%rps)
    call hash_alloc_real_rank2(hash,pp%anorm)
    call hash_alloc_integer_rank2(hash,pp%inorm)
    call hash_alloc_real_rank2(hash,pp%anorm_so)
    call hash_alloc_integer_rank2(hash,pp%inorm_so)
    call hash_alloc_real_rank2(hash,pp%rad)
    call hash_alloc_real_rank2(hash,pp%radnl)
    call hash_alloc_real_rank2(hash,pp%vloctbl)
    call hash_alloc_real_rank2(hash,pp%dvloctbl)
    call hash_alloc_real_rank3(hash,pp%udvtbl)
    call hash_alloc_real_rank3(hash,pp%dudvtbl)
    call hash_alloc_real_rank2(hash,pp%rho_pp_tbl)
    call hash_alloc_real_rank2(hash,pp%rho_nlcc_tbl)
    call hash_alloc_real_rank2(hash,pp%tau_nlcc_tbl)
    call hash_alloc_real_rank3(hash,pp%upp_f)
    call hash_alloc_real_rank3(hash,pp%vpp_f)
    call hash_alloc_real_rank3(hash,pp%vpp_f_so)
    call hash_alloc_real_rank2(hash,pp%upp)
    call hash_alloc_real_rank2(hash,pp%dupp)
    call hash_alloc_real_rank2(hash,pp%vpp)
    call hash_alloc_real_rank2(hash,pp%dvpp)
    call hash_alloc_real_rank2(hash,pp%vpp_so)
    call hash_alloc_real_rank2(hash,pp%dvpp_so)
    call hash_alloc_real_rank3(hash,pp%udvtbl_so)
    call hash_alloc_real_rank3(hash,pp%dudvtbl_so)
    call hash_alloc_real_rank1(hash,pp%rps_ao)
    call hash_alloc_integer_rank1(hash,pp%nrps_ao)
    call hash_alloc_real_rank3(hash,pp%upptbl_ao)
    call hash_alloc_real_rank3(hash,pp%dupptbl_ao)
  end subroutine hash_pp_info_for_dg_dc_seed

  subroutine hash_alloc_character_rank1(hash,values)
    integer(int64),intent(inout)::hash
    character(2),allocatable,intent(in)::values(:)
    integer::i
    call hash_integer(hash,merge(1,0,allocated(values)))
    if(.not.allocated(values))return
    call hash_integer(hash,lbound(values,1));call hash_integer(hash,ubound(values,1))
    do i=lbound(values,1),ubound(values,1);call hash_character(hash,values(i));enddo
  end subroutine hash_alloc_character_rank1

  subroutine hash_alloc_integer_rank1(hash,values)
    integer(int64),intent(inout)::hash
    integer,allocatable,intent(in)::values(:)
    integer::i
    call hash_integer(hash,merge(1,0,allocated(values)))
    if(.not.allocated(values))return
    call hash_integer(hash,lbound(values,1));call hash_integer(hash,ubound(values,1))
    do i=lbound(values,1),ubound(values,1);call hash_integer(hash,values(i));enddo
  end subroutine hash_alloc_integer_rank1

  subroutine hash_alloc_integer_rank2(hash,values)
    integer(int64),intent(inout)::hash
    integer,allocatable,intent(in)::values(:,:)
    integer::i,j
    call hash_integer(hash,merge(1,0,allocated(values)))
    if(.not.allocated(values))return
    do i=1,2
      call hash_integer(hash,lbound(values,i));call hash_integer(hash,ubound(values,i))
    enddo
    do j=lbound(values,2),ubound(values,2);do i=lbound(values,1),ubound(values,1)
      call hash_integer(hash,values(i,j))
    enddo;enddo
  end subroutine hash_alloc_integer_rank2

  subroutine hash_alloc_integer_rank3(hash,values)
    integer(int64),intent(inout)::hash
    integer,allocatable,intent(in)::values(:,:,:)
    integer::i,j,k,axis
    call hash_integer(hash,merge(1,0,allocated(values)))
    if(.not.allocated(values))return
    do axis=1,3
      call hash_integer(hash,lbound(values,axis));call hash_integer(hash,ubound(values,axis))
    enddo
    do k=lbound(values,3),ubound(values,3);do j=lbound(values,2),ubound(values,2)
      do i=lbound(values,1),ubound(values,1);call hash_integer(hash,values(i,j,k));enddo
    enddo;enddo
  end subroutine hash_alloc_integer_rank3

  subroutine hash_alloc_integer_rank5(hash,values)
    integer(int64),intent(inout)::hash
    integer,allocatable,intent(in)::values(:,:,:,:,:)
    integer::i1,i2,i3,i4,i5,axis
    call hash_integer(hash,merge(1,0,allocated(values)))
    if(.not.allocated(values))return
    do axis=1,5
      call hash_integer(hash,lbound(values,axis));call hash_integer(hash,ubound(values,axis))
    enddo
    do i5=lbound(values,5),ubound(values,5);do i4=lbound(values,4),ubound(values,4)
      do i3=lbound(values,3),ubound(values,3);do i2=lbound(values,2),ubound(values,2)
        do i1=lbound(values,1),ubound(values,1);call hash_integer(hash,values(i1,i2,i3,i4,i5));enddo
      enddo;enddo
    enddo;enddo
  end subroutine hash_alloc_integer_rank5

  subroutine hash_alloc_real_rank1(hash,values)
    integer(int64),intent(inout)::hash
    real(8),allocatable,intent(in)::values(:)
    integer::i
    call hash_integer(hash,merge(1,0,allocated(values)))
    if(.not.allocated(values))return
    call hash_integer(hash,lbound(values,1));call hash_integer(hash,ubound(values,1))
    do i=lbound(values,1),ubound(values,1);call hash_real(hash,values(i));enddo
  end subroutine hash_alloc_real_rank1

  subroutine hash_alloc_real_rank2(hash,values)
    integer(int64),intent(inout)::hash
    real(8),allocatable,intent(in)::values(:,:)
    integer::i,j,axis
    call hash_integer(hash,merge(1,0,allocated(values)))
    if(.not.allocated(values))return
    do axis=1,2
      call hash_integer(hash,lbound(values,axis));call hash_integer(hash,ubound(values,axis))
    enddo
    do j=lbound(values,2),ubound(values,2);do i=lbound(values,1),ubound(values,1)
      call hash_real(hash,values(i,j))
    enddo;enddo
  end subroutine hash_alloc_real_rank2

  subroutine hash_alloc_real_rank3(hash,values)
    integer(int64),intent(inout)::hash
    real(8),allocatable,intent(in)::values(:,:,:)
    integer::i,j,k,axis
    call hash_integer(hash,merge(1,0,allocated(values)))
    if(.not.allocated(values))return
    do axis=1,3
      call hash_integer(hash,lbound(values,axis));call hash_integer(hash,ubound(values,axis))
    enddo
    do k=lbound(values,3),ubound(values,3);do j=lbound(values,2),ubound(values,2)
      do i=lbound(values,1),ubound(values,1);call hash_real(hash,values(i,j,k));enddo
    enddo;enddo
  end subroutine hash_alloc_real_rank3

  subroutine run_dg_hybrid_continuation_ground_state_for_main
    call run_dg_overlapping_wannier_ground_state_for_main
  end subroutine run_dg_hybrid_continuation_ground_state_for_main

  subroutine run_dg_hybrid_divided_ground_state_for_main
    type(s_dg_hybrid_fragment_wannier_cache)::fragment_cache
    type(s_dg_hybrid_core_selection)::core_selection
    type(s_dg_hybrid_selected_catalog)::selected_catalog
    type(s_dg_hybrid_basis_catalog)::pw_catalog
    type(s_dg_hybrid_fragment_basis)::projected_basis
    type(s_dg_hybrid_fragment_basis),allocatable::fragment_bases(:)
    type(s_dg_hybrid_projection_factorization_receipt)::projection_receipt
    type(s_dg_hybrid_production_face_trace),allocatable::production_faces(:)
    type(s_dg_hybrid_support_operator)::support_operators(3)
    type(s_dg_hybrid_production_support_receipt)::support_receipt
    type(s_dg_hybrid_admission_report)::admission_report
    type(s_dg_hybrid_fragment_subspace_state)::fragment_state
    type(s_dg_hybrid_fragment_candidate_catalog)::candidate_catalog
    type(s_dg_hybrid_dc_reference)::dc_reference
    type(s_dg_hybrid_core_projection_report)::core_projection_report
    type(s_dg_hybrid_fixed_payload)::fixed_payload
    type(s_dg_hybrid_fragment_preconditioner)::fragment_preconditioner
    type(s_dg_hybrid_preconditioner_key)::preconditioner_key
    integer::nproc,rank,ierr,status,p,q,axis,index3(3),raw_grid(3),core_grid(3),global_point_count,&
      local_basis_count,total_basis_count,face_count,initial_count,guard_count,candidate_count,&
      pw_candidate_count,global_column
    integer,allocatable::fragment_ids(:),core_fragment_ids(:),row_action(:,:),basis_counts(:),&
      basis_displacements(:),effective_basis_ids(:),basis_owner(:),basis_fragment(:),&
      fragment_origins(:,:),fragment_sizes(:,:),projector_offsets(:),selected_seeds(:)
    integer(int64)::raw_count,byte_limit,pw_workspace,window_workspace,basis_workspace,&
      pw_fingerprint,window_fingerprint,basis_fingerprint,frame_fingerprint,face_fingerprint,&
      global_basis_fingerprint,metric_fingerprint,interface_fingerprint,directory_fingerprint,&
      local_potential_fingerprint,preconditioner_fingerprint
    integer(int64)::support_fingerprints(3)
    integer(int64),allocatable::candidate_grid_ids(:),core_ids(:),gathered_basis_ids(:),projector_grid_ids(:)
    complex(8),allocatable::buffer_candidates(:,:),projector_candidates(:,:),reference_frame(:,:),&
      interior_values(:,:),interior_gradients(:,:,:),interior_kinetic_action(:,:),kinetic_rows(:,:),&
      metric_rows(:,:),local_potential_rows(:,:),nonlocal_rows(:,:),interface_components(:,:,:),&
      interface_rows(:,:),fragment_h(:,:),fragment_s(:,:),seed_coefficients(:,:)
    real(8),allocatable::core_lower(:,:),core_extent(:,:),atom_positions(:,:),raw_weight(:),&
      raw_gradient(:,:),partition_weight(:),partition_gradient(:,:),box_windows(:,:),&
      core_coordinates(:,:),buffer_coordinates(:,:),core_windows(:,:),buffer_windows(:,:),&
      g_vectors(:,:),core_weights(:),projector_weights(:),unit_potential(:),local_potential(:),&
      initial_density(:),converged_density(:)
    complex(8),allocatable::projector_support_values(:)
    real(8)::axis_weight(3),axis_gradient(3),coordinate,sum_defect,gradient_defect,&
      denominator,convergence_value,electron_defect
    real(8)::volume_diagnostics(4),local_potential_diagnostics(2)
    real(8)::reciprocal_rotation(3,3,1),fragment_lattice(3,3),fragment_reciprocal_lattice(3,3)
    integer,allocatable::payload_owner(:),payload_fragment(:),payload_local_slot(:),payload_generation(:),&
      interior_fragment(:)
    character(8),allocatable::atom_symbols(:)
    integer::scf_iterations
    logical::ok,collective_ok
    character(512)::message

    call MPI_Comm_rank(dc%icomm_tot,rank,ierr);call MPI_Comm_size(dc%icomm_tot,nproc,ierr)
    ok=ierr==MPI_SUCCESS.and.nproc==dc%n_frag.and.dc%isize_frag==1.and.&
      dc%i_frag==rank+1.and..not.dc%optimized_fragment_geometry
    call comm_logical_and(ok,collective_ok,dc%icomm_tot)
    if(.not.collective_ok)error stop 'divided Hybrid production requires one rank per uniform DC fragment'
    ok=system%nspin==1.and.system%if_real_orbital.and.allocated(spsi%rwf).and.&
      allocated(energy%esp).and.allocated(system%rocc).and.allocated(system%kion).and.&
      allocated(system%Rion)
    call comm_logical_and(ok,collective_ok,dc%icomm_tot)
    if(.not.collective_ok)error stop 'divided Hybrid production requires a saved real Gamma DC seed'

    raw_grid=lg%num;core_grid=dc%nxyz_domain_frag(:,dc%i_frag)
    raw_count=product(int(raw_grid,int64))
    ok=all(raw_grid>0).and.all(core_grid>0).and.all(core_grid<=raw_grid).and.&
      raw_count>0_int64.and.raw_count<=int(huge(0),int64)
    call comm_logical_and(ok,collective_ok,dc%icomm_tot)
    if(.not.collective_ok)error stop 'divided Hybrid fragment grid extent is invalid'
    fragment_lattice=0d0;fragment_reciprocal_lattice=0d0
    do axis=1,3
      fragment_lattice(axis,axis)=dc%system_tot%hgs(axis)*real(raw_grid(axis),8)
      fragment_reciprocal_lattice(axis,axis)=2d0*acos(-1d0)/fragment_lattice(axis,axis)
    enddo
    allocate(candidate_grid_ids(int(raw_count)),buffer_candidates(0,int(raw_count)),&
      projector_candidates(0,int(raw_count)),core_lower(3,dc%n_frag),&
      core_extent(3,dc%n_frag),atom_symbols(system%nion),atom_positions(3,system%nion),stat=status)
    call comm_logical_and(status==0,ok,dc%icomm_tot)
    if(.not.ok)error stop 'divided Hybrid DC-to-Wannier staging allocation failed'
    candidate_grid_ids=[(int(p,int64),p=1,int(raw_count))]
    do p=1,system%nion
      atom_symbols(p)=pp%atom_symbol(system%kion(p))
    enddo
    atom_positions=system%Rion
    do p=1,dc%n_frag
      core_lower(:,p)=dc%rxyz_frag(:,p)
      core_extent(:,p)=real(dc%nxyz_domain_frag(:,p),8)*dc%system_tot%hgs
    enddo
    byte_limit=8_int64*1024_int64*1024_int64*1024_int64
    call build_dg_hybrid_fragment_wannier_from_dc_seed(dc%icomm_tot,dc%icomm_frag,info%icomm_o,&
      dc%i_frag,1,'dgfw',raw_grid,[1,1,1],raw_grid,&
      spsi%rwf,energy%esp,system%rocc,system%hvol,candidate_grid_ids,buffer_candidates,&
      projector_candidates,dg_dc_metric_rank_tolerance,fragment_lattice,fragment_reciprocal_lattice,&
      atom_symbols,atom_positions,wannier_num_iter,dg_ow_localization_gradient_tolerance,&
      byte_limit,fragment_cache,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local DC-to-Wannier construction failed'
    endif
    call select_dg_hybrid_core_wannier(dc%icomm_tot,dc%i_frag,fragment_cache,&
      fragment_lattice,dc%rxyz_frag(:,dc%i_frag),dc%system_tot%primitive_a,[0d0,0d0,0d0],&
      core_lower,core_extent,raw_grid,core_grid,dc%lg_tot%num,dc%jxyz_tot,core_selection,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local core-center selection failed'
    endif
    call prepare_dg_hybrid_selected_catalog(dc%icomm_tot,dc%i_frag,fragment_cache,&
      core_selection,selected_catalog,ok,message)
    if(.not.ok.or.core_selection%selected_count<1)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local selected Wannier catalog is empty or invalid'
    endif
    global_point_count=product(dc%lg_tot%num)
    allocate(raw_weight(int(raw_count)),raw_gradient(3,int(raw_count)),&
      partition_weight(int(raw_count)),partition_gradient(3,int(raw_count)),&
      box_windows(1,int(raw_count)),buffer_coordinates(3,int(raw_count)),&
      fragment_ids(1),core_ids(size(core_selection%core_row_slots)),&
      core_fragment_ids(size(core_selection%core_row_slots)),&
      core_coordinates(3,size(core_selection%core_row_slots)),&
      core_weights(size(core_selection%core_row_slots)),&
      row_action(size(core_selection%core_row_slots),1),stat=status)
    call comm_logical_and(status==0,collective_ok,dc%icomm_tot)
    if(.not.collective_ok)error stop 'divided Hybrid selected-basis staging allocation failed'
    do p=1,int(raw_count)
      index3(1)=modulo(p-1,raw_grid(1))+1
      index3(2)=modulo((p-1)/raw_grid(1),raw_grid(2))+1
      index3(3)=(p-1)/(raw_grid(1)*raw_grid(2))+1
      do axis=1,3
        if(index3(axis)<=core_grid(axis).or.dc%nxyz_buffer(axis)==0)then
          axis_weight(axis)=1d0;axis_gradient(axis)=0d0
        else if(index3(axis)<=core_grid(axis)+dc%nxyz_buffer(axis))then
          coordinate=real(core_grid(axis)+dc%nxyz_buffer(axis)+1-index3(axis),8)/&
            real(dc%nxyz_buffer(axis)+1,8)
          axis_weight(axis)=coordinate**2*(3d0-2d0*coordinate)
          axis_gradient(axis)=-6d0*coordinate*(1d0-coordinate)/&
            (real(dc%nxyz_buffer(axis)+1,8)*system%hgs(axis))
        else
          coordinate=real(index3(axis)-core_grid(axis)-dc%nxyz_buffer(axis),8)/&
            real(dc%nxyz_buffer(axis)+1,8)
          axis_weight(axis)=coordinate**2*(3d0-2d0*coordinate)
          axis_gradient(axis)=6d0*coordinate*(1d0-coordinate)/&
            (real(dc%nxyz_buffer(axis)+1,8)*system%hgs(axis))
        endif
      enddo
      raw_weight(p)=product(axis_weight)
      raw_gradient(1,p)=axis_gradient(1)*axis_weight(2)*axis_weight(3)
      raw_gradient(2,p)=axis_weight(1)*axis_gradient(2)*axis_weight(3)
      raw_gradient(3,p)=axis_weight(1)*axis_weight(2)*axis_gradient(3)
    enddo
    call build_dg_smooth_partition_of_unity(dc%icomm_tot,core_selection%physical_grid_ids,&
      raw_weight,raw_gradient,partition_weight,partition_gradient,sum_defect,gradient_defect,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local partition of unity failed'
    endif
    box_windows(1,:)=partition_weight;fragment_ids(1)=dc%i_frag
    core_ids=core_selection%physical_grid_ids(core_selection%core_row_slots)
    core_fragment_ids=dc%i_frag;core_weights=system%hvol;row_action(:,1)=int(core_ids)
    reciprocal_rotation=0d0
    do axis=1,3;reciprocal_rotation(axis,axis,1)=1d0;enddo
    do p=1,size(core_ids)
      core_coordinates(1,p)=real(modulo(core_ids(p)-1_int64,int(dc%lg_tot%num(1),int64)),8)*&
        dc%system_tot%hgs(1)
      core_coordinates(2,p)=real(modulo((core_ids(p)-1_int64)/int(dc%lg_tot%num(1),int64),&
        int(dc%lg_tot%num(2),int64)),8)*dc%system_tot%hgs(2)
      core_coordinates(3,p)=real((core_ids(p)-1_int64)/&
        (int(dc%lg_tot%num(1),int64)*int(dc%lg_tot%num(2),int64)),8)*dc%system_tot%hgs(3)
    enddo
    do p=1,int(raw_count)
      buffer_coordinates(1,p)=real(modulo(core_selection%physical_grid_ids(p)-1_int64,&
        int(dc%lg_tot%num(1),int64)),8)*dc%system_tot%hgs(1)
      buffer_coordinates(2,p)=real(modulo((core_selection%physical_grid_ids(p)-1_int64)/&
        int(dc%lg_tot%num(1),int64),int(dc%lg_tot%num(2),int64)),8)*dc%system_tot%hgs(2)
      buffer_coordinates(3,p)=real((core_selection%physical_grid_ids(p)-1_int64)/&
        (int(dc%lg_tot%num(1),int64)*int(dc%lg_tot%num(2),int64)),8)*dc%system_tot%hgs(3)
    enddo
    call build_dg_hybrid_production_pw_basis(dc%icomm_tot,global_point_count,dc%n_frag,&
      fragment_ids,core_selection%physical_grid_ids,box_windows,core_ids,core_fragment_ids,&
      core_coordinates,row_action,dc%system_tot%primitive_b,reciprocal_rotation,wannier_pw_cutoff,&
      16,dg_dc_metric_rank_tolerance,core_windows,g_vectors,pw_catalog,pw_workspace,&
      pw_fingerprint,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local PW catalog construction failed'
    endif
    call redistribute_dg_hybrid_fragment_windows(dc%icomm_tot,global_point_count,dc%n_frag,&
      fragment_ids,core_selection%physical_grid_ids,box_windows,core_selection%physical_grid_ids,&
      buffer_windows,window_workspace,window_fingerprint,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local buffer-window distribution failed'
    endif
    call build_dg_hybrid_projected_local_fragment_basis(dc%icomm_tot,global_point_count,dc%n_frag,&
      dc%i_frag,core_ids,core_weights,core_coordinates,core_windows,&
      core_selection%physical_grid_ids,selected_catalog%local_values,buffer_coordinates,buffer_windows,&
      pw_catalog,g_vectors,selected_catalog%wannier_owner,16,dg_dc_metric_rank_tolerance,&
      selected_catalog%fingerprint,projected_basis,basis_workspace,basis_fingerprint,ok,message,&
      basis_generation=core_selection%basis_generation,projection_receipt=projection_receipt)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local selected WF+PW projection failed'
    endif
    local_basis_count=size(projected_basis%global_ids)
    allocate(basis_counts(nproc),basis_displacements(nproc),fragment_bases(dc%n_frag),&
      fragment_origins(3,dc%n_frag),fragment_sizes(3,dc%n_frag),stat=status)
    call comm_logical_and(status==0,collective_ok,dc%icomm_tot)
    if(.not.collective_ok)error stop 'divided Hybrid basis-directory allocation failed'
    call MPI_Allgather(local_basis_count,1,MPI_INTEGER,basis_counts,1,MPI_INTEGER,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)error stop 'divided Hybrid basis-count exchange failed'
    basis_displacements(1)=0
    do p=2,nproc;basis_displacements(p)=basis_displacements(p-1)+basis_counts(p-1);enddo
    total_basis_count=sum(basis_counts)
    allocate(gathered_basis_ids(total_basis_count),effective_basis_ids(total_basis_count),stat=status)
    call comm_logical_and(status==0,collective_ok,dc%icomm_tot)
    if(.not.collective_ok)error stop 'divided Hybrid global basis inventory allocation failed'
    call MPI_Allgatherv(projected_basis%global_ids,local_basis_count,MPI_INTEGER8,gathered_basis_ids,&
      basis_counts,basis_displacements,MPI_INTEGER8,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.any([(count(gathered_basis_ids==int(p,int64))/=1,p=1,total_basis_count)]))&
      error stop 'divided Hybrid basis IDs are not a complete single-owner inventory'
    effective_basis_ids=[(p,p=1,total_basis_count)]
    do p=1,dc%n_frag
      if(p==dc%i_frag)then
        fragment_bases(p)=projected_basis
      else
        fragment_bases(p)%fragment_id=0;fragment_bases(p)%generation=projected_basis%generation
        allocate(fragment_bases(p)%global_ids(0),fragment_bases(p)%sector(0),&
          fragment_bases(p)%buffer_point_ids(0),fragment_bases(p)%buffer_values(0,0))
      endif
    enddo
    fragment_origins=dc%ixyz_frag;fragment_sizes=dc%nxyz_domain_frag
    call freeze_dg_hybrid_basis_directory(dc%icomm_tot,fragment_bases,effective_basis_ids,&
      basis_owner,basis_fragment,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local basis directory failed'
    endif
    call materialize_dg_hybrid_production_face_collection(dc%icomm_tot,fragment_origins,fragment_sizes,&
      dc%lg_tot%num,dc%system_tot%hgs,stencil%coef_nab,fragment_bases,basis_owner,basis_fragment,&
      effective_basis_ids,production_faces,ok,message,face_count,face_fingerprint)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local production face materialization failed'
    endif
    allocate(projector_offsets(ppg%nlma+1),projector_weights(ppg%nlma),stat=status)
    call comm_logical_and(status==0,collective_ok,dc%icomm_tot)
    if(.not.collective_ok)error stop 'divided Hybrid projector directory allocation failed'
    projector_offsets(1)=1
    do p=1,ppg%nlma
      projector_offsets(p+1)=projector_offsets(p)+ppg%mps(ppg%ia_tbl(p))
      projector_weights(p)=abs(system%hvol*ppg%rinv_uvu(p))
    enddo
    allocate(projector_grid_ids(projector_offsets(ppg%nlma+1)-1),&
      projector_support_values(projector_offsets(ppg%nlma+1)-1),stat=status)
    call comm_logical_and(status==0,collective_ok,dc%icomm_tot)
    if(.not.collective_ok)error stop 'divided Hybrid projector support allocation failed'
    q=0
    do p=1,ppg%nlma
      do axis=1,ppg%mps(ppg%ia_tbl(p))
        q=q+1
        index3=ppg%jxyz(:,axis,ppg%ia_tbl(p))
        projector_grid_ids(q)=int(dc%jxyz_tot(index3(1),1),int64)+&
          int(dc%lg_tot%num(1),int64)*(int(dc%jxyz_tot(index3(2),2)-1,int64)+&
          int(dc%lg_tot%num(2),int64)*int(dc%jxyz_tot(index3(3),3)-1,int64))
        projector_support_values(q)=ppg%uV(axis,p)
      enddo
    enddo
    call prepare_dg_hybrid_production_support(dc%icomm_tot,dc%i_frag,projected_basis%generation,&
      dc%lg_tot%num,stencil%coef_nab,projected_basis,production_faces,face_count,face_fingerprint,&
      projector_offsets,projector_grid_ids,projector_support_values,projector_weights,&
      support_operators,support_fingerprints,support_receipt,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local production support admission failed'
    endif
    initial_count=max(1,count(fragment_cache%physical_dc_seed_occupations>1d-10));guard_count=1
    call prepare_dg_hybrid_selected_trial(dc%icomm_tot,dc%i_frag,fragment_cache,core_selection,&
      projected_basis,projection_receipt,support_operators,support_fingerprints,core_weights,&
      [dg_dc_metric_rank_tolerance,dg_dc_gs_subspace_tolerance,dg_dc_gs_electron_count_tolerance,&
       dg_dc_gs_electron_count_tolerance],&
      [dg_ow_boundary_value_tolerance,dg_ow_boundary_gradient_tolerance,dg_dc_gs_subspace_tolerance],&
      wannier_pw_cutoff,initial_count,guard_count,dg_dc_gs_subspace_tolerance,&
      dg_dc_gs_orthogonality_tolerance,fragment_state,selected_seeds,admission_report,ok,message,&
      require_seed_reproduction=.false.)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local selected trial admission failed'
    endif
    if(rank==0)write(*,'(a,l1,a,i0,a,3es12.4)')&
      '[DG-HYBRID-DIVIDED-SEED] reproduction_required=',admission_report%seed_reproduction_required,&
      ' selected=',size(selected_seeds),' orbital/density/electron=',admission_report%core%orbital_residual,&
      admission_report%core%density_defect,admission_report%core%electron_defect
    call export_dg_hybrid_selected_basis_frame(dc%icomm_tot,dc%i_frag,fragment_cache,&
      core_selection,projected_basis,projection_receipt,reference_frame,frame_fingerprint,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local fixed reference frame export failed'
    endif
    allocate(interior_fragment(size(core_ids)),unit_potential(size(core_ids)),&
      local_potential(size(core_ids)),stat=status)
    call comm_logical_and(status==0,collective_ok,dc%icomm_tot)
    if(.not.collective_ok)error stop 'divided Hybrid volume staging allocation failed'
    interior_fragment=dc%i_frag;unit_potential=1d0
    do p=1,size(core_ids)
      q=core_selection%core_row_slots(p)-1
      index3(1)=modulo(q,raw_grid(1))+1;q=q/raw_grid(1)
      index3(2)=modulo(q,raw_grid(2))+1;index3(3)=q/raw_grid(2)+1
      local_potential(p)=v_local(1)%f(index3(1),index3(2),index3(3))
    enddo
    call materialize_dg_hybrid_production_interior(dc%icomm_tot,dc%lg_tot%num,stencil%coef_nab,&
      stencil%coef_lap0,stencil%coef_lap,fragment_bases,basis_owner,basis_fragment,effective_basis_ids,&
      core_ids,interior_fragment,interior_values,interior_gradients,interior_kinetic_action,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local production interior materialization failed'
    endif
    call assemble_dg_hybrid_broken_volume_rows(dc%icomm_tot,total_basis_count,&
      projected_basis%global_ids,basis_fragment,core_ids,interior_fragment,core_weights,interior_values,&
      interior_gradients,unit_potential,kinetic_rows,metric_rows,volume_diagnostics,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local kinetic/metric assembly failed'
    endif
    call assemble_dg_hybrid_local_potential_rows(dc%icomm_tot,total_basis_count,&
      projected_basis%global_ids,basis_fragment,core_ids,interior_fragment,core_weights,interior_values,&
      local_potential,local_potential_rows,local_potential_diagnostics,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local potential projection failed'
    endif
    call assemble_dg_hybrid_selected_nonlocal_rows(projected_basis,total_basis_count,nonlocal_rows,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local nonlocal projection failed'
    endif
    call assemble_dg_hybrid_production_interface_component_rows(dc%icomm_tot,total_basis_count,&
      projected_basis%global_ids,production_faces,dg_dc_gs_sipg_penalty_factor,&
      interface_components,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local SIPG interface assembly failed'
    endif
    allocate(interface_rows(size(projected_basis%global_ids),total_basis_count))
    interface_rows=sum(interface_components,dim=3)
    call ow_fingerprint_distributed_matrix(dc%icomm_tot,projected_basis%global_ids,metric_rows,&
      metric_fingerprint,ok)
    if(.not.ok)error stop 'fragment-local metric fingerprint failed'
    call ow_fingerprint_distributed_matrix(dc%icomm_tot,projected_basis%global_ids,interface_rows,&
      interface_fingerprint,ok)
    if(.not.ok)error stop 'fragment-local interface fingerprint failed'
    call MPI_Allreduce(basis_fingerprint,global_basis_fingerprint,1,MPI_INTEGER8,MPI_BXOR,&
      dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)error stop 'fragment-local basis fingerprint reduction failed'
    global_basis_fingerprint=ieor(global_basis_fingerprint,int(z'6A09E667F3BCC909',int64))
    if(global_basis_fingerprint==0_int64)global_basis_fingerprint=1_int64
    call freeze_dg_hybrid_single_owner_payload(dc%icomm_tot,dc%n_frag,projected_basis,metric_rows,&
      kinetic_rows,nonlocal_rows,interface_rows,global_basis_fingerprint,metric_fingerprint,&
      interface_fingerprint,fixed_payload,payload_owner,payload_fragment,payload_local_slot,&
      payload_generation,directory_fingerprint,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local variational payload freeze failed'
    endif
    preconditioner_key%fragment_id=dc%i_frag
    preconditioner_key%basis_generation=projected_basis%generation
    preconditioner_key%operator_epoch=1
    preconditioner_key%basis_fingerprint=admission_report%basis_fingerprint
    preconditioner_key%metric_fingerprint=admission_report%metric_fingerprint
    preconditioner_key%operator_fingerprint=fixed_payload%fingerprint
    preconditioner_key%reference_fingerprint=frame_fingerprint

    ! The extension inventory is expressed in the final projected coordinates.
    ! Reproject every immutable physical DC seed; never reinterpret a named WF
    ! column as a physical eigenstate and never pad a seed with zero PW entries.
    call export_dg_hybrid_dc_reference(dc%icomm_tot,dc%i_frag,fragment_cache,&
      core_selection,dc_reference,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local DC reference export failed'
    endif
    call project_dg_hybrid_core_seeds(dc%icomm_tot,&
      projected_basis%buffer_values(core_selection%core_row_slots,:),core_weights,&
      dc_reference%core_orbitals,dc_reference%occupations,&
      [dg_dc_metric_rank_tolerance,huge(1d0)/100d0,huge(1d0)/100d0,huge(1d0)/100d0],&
      core_selection%selected_count,wannier_pw_cutoff,seed_coefficients,&
      core_projection_report,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local DC seed reprojection failed'
    endif
    pw_candidate_count=count(projected_basis%sector==2)
    candidate_count=size(seed_coefficients,2)+pw_candidate_count
    allocate(candidate_catalog%coefficients(local_basis_count,candidate_count),&
      candidate_catalog%energies(candidate_count),candidate_catalog%ids(candidate_count),&
      candidate_catalog%source_kind(candidate_count),candidate_catalog%used(candidate_count),stat=status)
    call comm_logical_and(status==0,collective_ok,dc%icomm_tot)
    if(.not.collective_ok)error stop 'divided Hybrid extension catalog allocation failed'
    candidate_catalog%coefficients=0d0
    candidate_catalog%coefficients(:,:size(seed_coefficients,2))=seed_coefficients
    candidate_catalog%energies(:size(seed_coefficients,2))=dc_reference%energies
    candidate_catalog%source_kind(:size(seed_coefficients,2))=fragment_seed
    candidate_catalog%used=.false.
    candidate_catalog%used(selected_seeds)=.true.
    q=size(seed_coefficients,2)
    do p=1,local_basis_count
      if(projected_basis%sector(p)/=2)cycle
      q=q+1;candidate_catalog%coefficients(p,q)=1d0
      global_column=0
      do axis=1,size(payload_fragment)
        if(payload_fragment(axis)==dc%i_frag.and.payload_local_slot(axis)==p)then
          global_column=axis;exit
        endif
      enddo
      if(global_column<1)error stop 'divided Hybrid PW kinetic column is missing'
      denominator=real(metric_rows(p,global_column),8)
      if(.not.ieee_is_finite(denominator).or.denominator<=dg_dc_metric_rank_tolerance)&
        error stop 'divided Hybrid PW candidate has an unresolved metric norm'
      candidate_catalog%energies(q)=real(kinetic_rows(p,global_column),8)/denominator
      candidate_catalog%source_kind(q)=fragment_pw
    enddo
    candidate_catalog%ids=[(int(p,int64),p=1,candidate_count)]
    candidate_catalog%fragment_id=dc%i_frag
    candidate_catalog%basis_generation=projected_basis%generation
    candidate_catalog%basis_fingerprint=admission_report%basis_fingerprint
    candidate_catalog%metric_fingerprint=admission_report%metric_fingerprint

    call prepare_dg_hybrid_schwarz_candidate_inventory(dc%icomm_tot,candidate_catalog,&
      bounded_schwarz_candidate_ids,bounded_schwarz_candidate_energies,&
      bounded_schwarz_candidate_vectors,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'common Schwarz candidate inventory preparation failed'
    endif
    bounded_mapping_fingerprint=directory_fingerprint
    bounded_candidate_fingerprint=ieor(fixed_payload%fingerprint,ishftc(directory_fingerprint,23))
    if(bounded_candidate_fingerprint==0_int64)bounded_candidate_fingerprint=global_basis_fingerprint
    call initialize_dg_hybrid_schwarz_state(dc%icomm_tot,dc%i_frag,dc%n_frag,&
      projected_basis%generation,dc%elec_num_tot,300d0,2d0,guard_count,1d-10,&
      dg_dc_gs_subspace_tolerance,bounded_mapping_fingerprint,bounded_candidate_fingerprint,&
      bounded_schwarz_candidate_ids,bounded_schwarz_candidate_energies,&
      bounded_schwarz_candidate_vectors,bounded_schwarz_state,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'common Schwarz state initialization failed'
    endif
    call build_dg_hybrid_schwarz_schedule(dc%icomm_tot,dc%i_frag,projected_basis%generation,&
      projected_basis%global_ids,payload_owner,payload_fragment,payload_local_slot,payload_generation,&
      interface_rows,directory_fingerprint,face_fingerprint,bounded_mapping_fingerprint,&
      bounded_schwarz_schedule,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'immutable Schwarz neighbor schedule construction failed'
    endif

    ! Publish the immutable frame/operator context used by sibling callbacks.
    ! Each callback is MPI_COMM_SELF inside its fragment; only the common
    ! occupation and outer density loop communicate over dc%icomm_tot.
    divided_fragment_basis=projected_basis
    bounded_projection_receipt=projection_receipt
    bounded_admission_report=admission_report
    bounded_fixed_payload=fixed_payload
    bounded_preconditioner_key=preconditioner_key
    bounded_interior_values=interior_values
    bounded_local_potential_rows=local_potential_rows
    bounded_reference_frame=reference_frame
    bounded_core_weights=core_weights
    allocate(bounded_point_weights(size(projected_basis%buffer_point_ids)),source=system%hvol)
    bounded_basis_fragment=payload_fragment
    bounded_basis_local_slot=payload_local_slot
    bounded_basis_generation=payload_generation
    bounded_interior_fragment=interior_fragment
    bounded_core_ids=core_ids
    bounded_global_basis_fingerprint=global_basis_fingerprint
    bounded_directory_fingerprint=directory_fingerprint
    bounded_face_fingerprint=face_fingerprint
    bounded_frame_fingerprint=frame_fingerprint
    bounded_selection_fingerprint=core_selection%fingerprint
    divided_mixing_basis_generation=projected_basis%generation
    divided_mixing_inventory_fingerprint=bounded_schwarz_state%fingerprint
    if(allocated(ow_core_ids))deallocate(ow_core_ids)
    if(allocated(ow_core_weights))deallocate(ow_core_weights)
    allocate(ow_core_ids,source=core_ids);allocate(ow_core_weights,source=core_weights)
    ow_global_grid_count=int(global_point_count,int64)
    if(allocated(ow_divided_core_mask))deallocate(ow_divided_core_mask)
    allocate(ow_divided_core_mask(size(projected_basis%buffer_point_ids)),source=.false.)
    ow_divided_core_mask(core_selection%core_row_slots)=.true.

    call prepare_dg_hybrid_divided_dc_controls(dc,ow_hybrid_divided_convergence,&
      ow_hybrid_divided_threshold,ow_hybrid_divided_total_density,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local divided DC control preparation failed'
    endif
    allocate(initial_density(size(core_ids)))
    do p=1,size(core_ids)
      q=int(core_ids(p)-1_int64)
      index3(1)=modulo(q,dc%lg_tot%num(1))+1;q=q/dc%lg_tot%num(1)
      index3(2)=modulo(q,dc%lg_tot%num(2))+1;index3(3)=q/dc%lg_tot%num(2)+1
      initial_density(p)=ow_hybrid_divided_total_density(index3(1),index3(2),index3(3))
    enddo
    call run_dg_hybrid_divided_scf(dc%icomm_tot,global_point_count,core_ids,initial_density,&
      dc%system_tot%hvol,ow_hybrid_divided_convergence,ow_hybrid_divided_threshold,&
      update_dg_hybrid_divided_potential,solve_dg_hybrid_schwarz_fragments,&
      assemble_dg_hybrid_schwarz_core_density,mix_dg_hybrid_divided_density,nscf,core_weights,&
      dc%elec_num_tot,dg_dc_gs_electron_count_tolerance,converged_density,scf_iterations,&
      convergence_value,electron_defect,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local divided WF+PW SCF failed'
    endif
    if(rank==0)write(*,'(a,2(a,es12.4),a,i0)')'[DG-HYBRID-DIVIDED] selected basis prepared',&
      ' partition_sum_defect=',sum_defect,' partition_gradient_defect=',gradient_defect,&
      ' global_basis_count=',size(projected_basis%global_ids)
    if(rank==0)write(*,'(a,i0,2(a,es12.4))')'[DG-HYBRID-DIVIDED] local SCF converged iterations=',&
      scf_iterations,' density=',convergence_value,' electron_defect=',electron_defect
  end subroutine run_dg_hybrid_divided_ground_state_for_main

  subroutine solve_dg_hybrid_schwarz_fragments(iteration,callback_ok)
    integer,intent(in)::iteration
    logical,intent(out)::callback_ok
    real(8),allocatable::core_potential(:)
    real(8),allocatable::occupations(:),energies(:)
    real(8)::diagnostics(2)
    integer(8)::potential_fingerprint
    integer::local_iterations,extensions,rank_local,ierr_local
    logical::converged,rolled_back,extended
    character(512)::solver_message

    callback_ok=.false.
    if(.not.bounded_schwarz_state%valid.or..not.bounded_schwarz_schedule%valid.or.&
      .not.allocated(bounded_core_ids))return
    allocate(core_potential(size(bounded_core_ids)))
    call extract_dg_hybrid_core_local_potential(bounded_core_ids,core_potential,callback_ok)
    if(.not.callback_ok)return
    call assemble_dg_hybrid_local_potential_rows(dc%icomm_tot,bounded_fixed_payload%global_basis_count,&
      divided_fragment_basis%global_ids,bounded_basis_fragment,bounded_core_ids,&
      bounded_interior_fragment,bounded_core_weights,bounded_interior_values,core_potential,&
      bounded_local_potential_rows,diagnostics,callback_ok,solver_message)
    if(.not.callback_ok)then
      write(error_unit,'(a,a)')'bounded fragment potential projection: ',trim(solver_message);return
    endif
    call ow_fingerprint_distributed_matrix(dc%icomm_tot,divided_fragment_basis%global_ids,&
      bounded_local_potential_rows,potential_fingerprint,callback_ok)
    if(.not.callback_ok)return
    call assemble_dg_hybrid_schwarz_local_preconditioner_blocks(callback_ok)
    if(.not.callback_ok)then
      write(error_unit,'(a)')'Schwarz local preconditioner block assembly failed';return
    endif
    bounded_last_peer_exchange_count=0
    call advance_dg_hybrid_schwarz_epoch(dc%icomm_tot,bounded_schwarz_state%basis_generation,&
      dg_hybrid_fragment_cg_steps,dg_dc_gs_intermediate_orbital_tolerance,&
      dg_dc_gs_orthogonality_tolerance,dg_dc_gs_allowed_residual_growth,&
      apply_dg_hybrid_schwarz_h,apply_dg_hybrid_schwarz_s,&
      apply_dg_hybrid_schwarz_preconditioner,bounded_schwarz_state,local_iterations,&
      divided_fragment_residual,divided_fragment_orthogonality,converged,rolled_back,&
      callback_ok,solver_message)
    if(.not.callback_ok)then
      divided_mixing_rollback_pending=rolled_back
      write(error_unit,'(a,a)')'bounded Schwarz update: ',trim(solver_message);return
    endif
    if(local_iterations>dg_hybrid_fragment_cg_steps)then
      callback_ok=.false.;write(error_unit,'(a)')'Schwarz CG step cap was exceeded';return
    endif
    extensions=0
    do
      call assign_dg_hybrid_schwarz_occupations(dc%icomm_tot,bounded_schwarz_state%basis_generation,&
        300d0,2d0,dc%elec_num_tot,1d-10,dg_dc_gs_subspace_tolerance,&
        bounded_schwarz_candidate_ids,bounded_schwarz_candidate_energies,&
        bounded_schwarz_candidate_vectors,apply_dg_hybrid_schwarz_h,apply_dg_hybrid_schwarz_s,&
        bounded_schwarz_state,occupations,energies,extended,callback_ok,solver_message)
      if(.not.callback_ok)then
        write(error_unit,'(a,a)')'global Schwarz occupations: ',trim(solver_message);return
      endif
      if(.not.extended)exit
      extensions=extensions+1
    enddo
    divided_mixing_inventory_fingerprint=bounded_schwarz_state%fingerprint
    divided_mixing_rollback_pending=.false.
    system%mu=bounded_schwarz_state%chemical_potential
    call MPI_Comm_rank(dc%icomm_tot,rank_local,ierr_local)
    if(ierr_local/=MPI_SUCCESS)then;callback_ok=.false.;return;endif
    if(rank_local==0)write(*,'(a,i0,3(a,i0),3(a,es12.4))')'[DG-HYBRID-SCHWARZ] epoch=',iteration,&
      ' neighbor_exchanges=',bounded_last_peer_exchange_count,' accepted_cg_steps=',local_iterations,&
      ' common_extensions=',extensions,' residual=',divided_fragment_residual,&
      ' electron_defect=',bounded_schwarz_state%electron_defect,&
      ' temperature=',bounded_schwarz_state%temperature
    callback_ok=local_iterations<=dg_hybrid_fragment_cg_steps
  end subroutine solve_dg_hybrid_schwarz_fragments

  subroutine assemble_dg_hybrid_schwarz_local_preconditioner_blocks(callback_ok)
    logical,intent(out)::callback_ok
    integer::global_column,local_slot,local_count
    local_count=size(divided_fragment_basis%global_ids);callback_ok=.false.
    if(allocated(bounded_fragment_h))deallocate(bounded_fragment_h)
    if(allocated(bounded_fragment_s))deallocate(bounded_fragment_s)
    allocate(bounded_fragment_h(local_count,local_count),bounded_fragment_s(local_count,local_count))
    bounded_fragment_h=0d0;bounded_fragment_s=0d0
    do global_column=1,bounded_fixed_payload%global_basis_count
      if(bounded_basis_fragment(global_column)/=dc%i_frag)cycle
      local_slot=bounded_basis_local_slot(global_column)
      bounded_fragment_h(:,local_slot)=bounded_fixed_payload%kinetic_rows(:,global_column)+&
        bounded_fixed_payload%nonlocal_rows(:,global_column)+&
        bounded_fixed_payload%interface_rows(:,global_column)+&
        bounded_local_potential_rows(:,global_column)
      bounded_fragment_s(:,local_slot)=bounded_fixed_payload%metric_rows(:,global_column)
    enddo
    callback_ok=all(ieee_is_finite(real(bounded_fragment_h))).and.&
      all(ieee_is_finite(aimag(bounded_fragment_h))).and.&
      all(ieee_is_finite(real(bounded_fragment_s))).and.all(ieee_is_finite(aimag(bounded_fragment_s)))
  end subroutine assemble_dg_hybrid_schwarz_local_preconditioner_blocks

  subroutine apply_dg_hybrid_schwarz_h(input,output,callback_ok)
    complex(8),intent(in)::input(:,:)
    complex(8),intent(out)::output(:,:)
    logical,intent(out)::callback_ok
    complex(8),allocatable::candidate(:,:)
    integer::peer_exchanges
    character(512)::operator_message
    call apply_dg_hybrid_schwarz_hamiltonian(dc%icomm_tot,bounded_schwarz_schedule,&
      bounded_schwarz_state%basis_generation,bounded_directory_fingerprint,bounded_face_fingerprint,&
      bounded_mapping_fingerprint,divided_fragment_basis%global_ids,bounded_basis_fragment,&
      bounded_basis_local_slot,bounded_fixed_payload%kinetic_rows,bounded_fixed_payload%nonlocal_rows,&
      bounded_fixed_payload%interface_rows,bounded_local_potential_rows,input,candidate,&
      peer_exchanges,callback_ok,operator_message)
    if(.not.callback_ok)then
      write(error_unit,'(a,a)')'Schwarz H application: ',trim(operator_message);return
    endif
    callback_ok=all(shape(candidate)==shape(output))
    if(callback_ok)output=candidate
    bounded_last_peer_exchange_count=max(bounded_last_peer_exchange_count,peer_exchanges)
  end subroutine apply_dg_hybrid_schwarz_h

  subroutine apply_dg_hybrid_schwarz_s(input,output,callback_ok)
    complex(8),intent(in)::input(:,:)
    complex(8),intent(out)::output(:,:)
    logical,intent(out)::callback_ok
    complex(8),allocatable::candidate(:,:)
    integer::peer_exchanges
    character(512)::operator_message
    call apply_dg_hybrid_schwarz_rows(dc%icomm_tot,bounded_schwarz_schedule,&
      bounded_schwarz_state%basis_generation,bounded_directory_fingerprint,bounded_face_fingerprint,&
      bounded_mapping_fingerprint,divided_fragment_basis%global_ids,bounded_basis_fragment,&
      bounded_basis_local_slot,bounded_fixed_payload%metric_rows,input,candidate,peer_exchanges,&
      callback_ok,operator_message)
    if(.not.callback_ok)then
      write(error_unit,'(a,a)')'Schwarz S application: ',trim(operator_message);return
    endif
    callback_ok=all(shape(candidate)==shape(output))
    if(callback_ok)output=candidate
  end subroutine apply_dg_hybrid_schwarz_s

  subroutine apply_dg_hybrid_schwarz_preconditioner(input,output,callback_ok)
    complex(8),intent(in)::input(:,:)
    complex(8),intent(out)::output(:,:)
    logical,intent(out)::callback_ok
    real(8)::scale
    integer::p
    callback_ok=all(shape(input)==shape(output)).and.&
      size(input,1)==size(bounded_fragment_h,1)
    if(.not.callback_ok)return
    do p=1,size(input,1)
      scale=max(1d0,abs(bounded_fragment_h(p,p)))
      output(p,:)=input(p,:)/scale
    enddo
    callback_ok=all(ieee_is_finite(real(output))).and.all(ieee_is_finite(aimag(output)))
  end subroutine apply_dg_hybrid_schwarz_preconditioner

  subroutine assemble_dg_hybrid_schwarz_core_density(core_density,electron_count,callback_ok)
    real(8),intent(out)::core_density(:),electron_count
    logical,intent(out)::callback_ok
    complex(8),allocatable::orbital_values(:)
    real(8)::local_electron_count
    integer::p,slot,ierr_local
    callback_ok=.false.;core_density=0d0;electron_count=0d0
    if(.not.bounded_schwarz_state%valid.or..not.allocated(bounded_schwarz_state%occupations))return
    if(size(core_density)/=size(bounded_core_ids))return
    allocate(orbital_values(bounded_schwarz_state%trial_count))
    do p=1,size(bounded_core_ids)
      slot=findloc(divided_fragment_basis%buffer_point_ids,bounded_core_ids(p),dim=1)
      if(slot<1)return
      orbital_values=matmul(divided_fragment_basis%buffer_values(slot,:),&
        bounded_schwarz_state%coefficients)
      core_density(p)=bounded_schwarz_state%wspin*&
        sum(bounded_schwarz_state%occupations*abs(orbital_values)**2)
    enddo
    local_electron_count=sum(bounded_core_weights*core_density)
    call MPI_Allreduce(local_electron_count,electron_count,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
      dc%icomm_tot,ierr_local)
    callback_ok=ierr_local==MPI_SUCCESS.and.all(ieee_is_finite(core_density)).and.&
      ieee_is_finite(electron_count)
  end subroutine assemble_dg_hybrid_schwarz_core_density

  subroutine prepare_dg_hybrid_schwarz_candidate_inventory(comm_arg,catalog,ids,energies,vectors,&
      callback_ok,callback_message)
    integer,intent(in)::comm_arg
    type(s_dg_hybrid_fragment_candidate_catalog),intent(in)::catalog
    integer(8),allocatable,intent(out)::ids(:)
    real(8),allocatable,intent(out)::energies(:)
    complex(8),allocatable,intent(out)::vectors(:,:)
    logical,intent(out)::callback_ok
    character(*),intent(out)::callback_message
    integer,allocatable::order(:)
    real(8),allocatable::local_energies(:)
    integer::local_count,common_count,nproc_local,ierr_local,p,q,selected,temp
    logical::local_ok,global_ok

    callback_ok=.false.;callback_message='invalid fragment Schwarz candidate catalog'
    local_ok=allocated(catalog%ids).and.allocated(catalog%energies).and.&
      allocated(catalog%coefficients)
    if(local_ok)local_ok=size(catalog%ids)>0.and.size(catalog%energies)==size(catalog%ids).and.&
      size(catalog%coefficients,2)==size(catalog%ids).and.size(catalog%coefficients,1)>0.and.&
      all(catalog%ids>0_int64).and.all(ieee_is_finite(catalog%energies)).and.&
      all(ieee_is_finite(real(catalog%coefficients))).and.&
      all(ieee_is_finite(aimag(catalog%coefficients)))
    call comm_logical_and(local_ok,global_ok,comm_arg)
    if(.not.global_ok)return
    local_count=size(catalog%ids);common_count=local_count
    call MPI_Allreduce(MPI_IN_PLACE,common_count,1,MPI_INTEGER,MPI_MIN,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.common_count<1)then
      callback_message='common Schwarz candidate capacity reduction failed';return
    endif
    allocate(order(local_count));order=[(p,p=1,local_count)]
    do p=1,local_count-1
      selected=p
      do q=p+1,local_count
        if(catalog%energies(order(q))<catalog%energies(order(selected)).or.&
          (catalog%energies(order(q))==catalog%energies(order(selected)).and.&
           catalog%ids(order(q))<catalog%ids(order(selected))))selected=q
      enddo
      if(selected/=p)then;temp=order(p);order(p)=order(selected);order(selected)=temp;endif
    enddo
    allocate(ids(common_count),energies(common_count),local_energies(common_count),&
      vectors(size(catalog%coefficients,1),common_count))
    ids=catalog%ids(order(:common_count));local_energies=catalog%energies(order(:common_count))
    vectors=catalog%coefficients(:,order(:common_count));energies=local_energies
    call MPI_Comm_size(comm_arg,nproc_local,ierr_local)
    if(ierr_local/=MPI_SUCCESS)then;callback_message='Schwarz candidate size query failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,energies,common_count,MPI_DOUBLE_PRECISION,MPI_SUM,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS)then;callback_message='Schwarz candidate energy reduction failed';return;endif
    energies=energies/real(nproc_local,8)
    local_ok=all(energies(2:)>=energies(:common_count-1)).and.&
      all([(count(ids==ids(p))==1,p=1,common_count)])
    call comm_logical_and(local_ok,global_ok,comm_arg)
    if(.not.global_ok)then;callback_message='common Schwarz candidate ordering is invalid';return;endif
    callback_ok=.true.;callback_message=''
  end subroutine prepare_dg_hybrid_schwarz_candidate_inventory

  subroutine assemble_dg_hybrid_selected_nonlocal_rows(fragment_basis,global_count,matrix_rows,ok,message)
    type(s_dg_hybrid_fragment_basis),intent(in)::fragment_basis
    integer,intent(in)::global_count
    complex(8),allocatable,intent(out)::matrix_rows(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::local_overlap(:,:),owned_overlap(:,:)
    integer,allocatable::local_atom_ids(:),local_ordinals(:)
    integer(int64),allocatable::projector_ids(:)
    real(8),allocatable::local_matrix_strength(:),local_action_strength(:),owned_matrix_strength(:)
    logical,allocatable::complete(:,:)
    integer::ilma,ia,j,ix,iy,iz,ix_tot,iy_tot,iz_tot,basis,position,ordinal,&
      total_projectors,ownership_count,local_bad,global_bad,ierr
    integer(int64)::point_id

    ok=.false.;message='';local_bad=0
    allocate(local_overlap(global_count,ppg%nlma),local_atom_ids(ppg%nlma),&
      local_ordinals(ppg%nlma),local_matrix_strength(ppg%nlma),local_action_strength(ppg%nlma))
    local_overlap=0d0;local_atom_ids=0;local_ordinals=0
    local_matrix_strength=0d0;local_action_strength=0d0
    do ilma=1,ppg%nlma
      ia=ppg%ia_tbl(ilma)
      call map_dc_atom_to_physical_atom(ia,local_atom_ids(ilma),ok)
      if(.not.ok)then;local_bad=1;cycle;endif
      ordinal=count(ppg%ia_tbl(1:ilma)==ia);local_ordinals(ilma)=ordinal
      local_matrix_strength(ilma)=system%hvol*ppg%rinv_uvu(ilma)
      local_action_strength(ilma)=ppg%rinv_uvu(ilma)
      do j=1,ppg%mps(ia)
        ix=ppg%jxyz(1,j,ia);iy=ppg%jxyz(2,j,ia);iz=ppg%jxyz(3,j,ia)
        ix_tot=dc%jxyz_tot(ix,1);iy_tot=dc%jxyz_tot(iy,2);iz_tot=dc%jxyz_tot(iz,3)
        point_id=int(ix_tot,int64)+int(dc%lg_tot%num(1),int64)*(&
          int(iy_tot-1,int64)+int(dc%lg_tot%num(2),int64)*int(iz_tot-1,int64))
        position=findloc(fragment_basis%buffer_point_ids,point_id,dim=1)
        if(position<=0)then;local_bad=1;cycle;endif
        do basis=1,size(fragment_basis%global_ids)
          local_overlap(int(fragment_basis%global_ids(basis)),ilma)=&
            local_overlap(int(fragment_basis%global_ids(basis)),ilma)+&
            ppg%uV(j,ilma)*fragment_basis%buffer_values(position,basis)
        enddo
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='selected basis omits nonlocal projector support';return
    endif
    call collect_dg_overlapping_wannier_projector_overlaps(dc%icomm_tot,global_count,local_atom_ids,&
      local_ordinals,local_matrix_strength,local_action_strength,local_overlap,projector_ids,&
      owned_matrix_strength,owned_overlap,total_projectors,ok,message)
    if(.not.ok)return
    allocate(complete(global_count,size(projector_ids)));complete=.true.
    call assemble_dg_overlapping_wannier_nonlocal_rows(dc%icomm_tot,global_count,&
      fragment_basis%global_ids,projector_ids,owned_matrix_strength,owned_overlap,complete,&
      int(total_projectors,int64),matrix_rows,ownership_count,ok,message)
    if(.not.ok)return
    if(ownership_count/=total_projectors)then
      ok=.false.;message='selected nonlocal projector ownership is incomplete'
    endif
  end subroutine assemble_dg_hybrid_selected_nonlocal_rows

  subroutine checked_ow_extent_product(extent,value,ok)
    integer,intent(in)::extent(3)
    integer(8),intent(out)::value
    logical,intent(out)::ok
    integer::axis
    value=1_8;ok=all(extent>0)
    if(.not.ok)return
    do axis=1,3
      if(value>huge(value)/int(extent(axis),8))then
        value=0_8;ok=.false.;return
      endif
      value=value*int(extent(axis),8)
    enddo
  end subroutine

  subroutine run_dg_overlapping_wannier_ground_state_for_main()
    complex(8),allocatable::global_seed_values(:,:),global_closed_core(:,:),&
      local_occupied_values(:,:),orbital_owned_full_values(:,:),center_local_buffer_values(:,:),&
      adapted_occupied_candidates(:,:),translation_adapted_occupied(:,:)
    complex(8),allocatable::orthonormal_lcfo_occupied(:,:)
    complex(8),allocatable::occupied_overlap_local(:,:),occupied_overlap_global(:,:)
    complex(8),allocatable::translation_hamiltonian_overlap_local(:,:),&
      translation_hamiltonian_overlap(:,:),translation_occupied_hamiltonian(:,:)
    complex(8),allocatable::w90_anchors(:,:),w90_m_matrix(:,:,:),w90_a_matrix(:,:),w90_seed_a_matrix(:,:),&
      w90_seed_representation(:,:),w90_transform(:,:)
    complex(8),allocatable::fixed_center_rows(:,:,:),fixed_center_representation(:,:),fixed_center_identity(:,:)
    complex(8),allocatable::translation_generator_rows(:,:,:),translation_gamma_rows(:,:),&
      translation_spatial_gamma_rows(:,:),&
      translation_sector_rows(:,:),translation_characters(:,:),translation_generator_characters(:,:),&
      translation_gamma_local_row(:),translation_gamma_global_row(:),&
      translation_lcfo_local_row(:),translation_lcfo_global_row(:)
    complex(8),allocatable::translation_w90_rows(:,:),translation_lcfo_rows(:,:),translation_w90_operator(:,:),&
      translation_lcfo_operator(:,:),translation_reference_rows(:,:),translation_reference_spatial(:,:),&
      translation_target_spatial(:,:),translation_aligned_spatial(:,:),translation_conjugate_spatial(:,:),&
      translation_phase(:),translation_orbit_rows(:,:),translation_transform_rows(:,:),&
      transformed_box_values(:,:),transformed_box_gradients(:,:,:)
    complex(8),allocatable::lcfo_fragment_contribution(:,:),lcfo_occupied_core(:,:),lcfo_reference_core(:,:)
    complex(8),allocatable::composed_tile_values(:,:),projector_buffer_tile(:,:)
    complex(8),allocatable::ow_direct_core_gradients(:,:,:),ow_neighbor_plus_values(:,:),&
      ow_neighbor_minus_values(:,:),ow_map_probe_values(:,:),ow_map_probe_gradients(:,:,:)
    complex(8),allocatable::ow_scalar_probe(:,:),ow_vector_probe(:,:,:),ow_scalar_representation(:,:,:)
    complex(8),allocatable::spectral_complement_generator_rows(:,:,:),spectral_complement_trial_rows(:,:),&
      spectral_trial_rows(:,:),&
      spectral_representative_vectors(:,:),spectral_basin_operator(:,:),spectral_wannier_action_rows(:,:,:),&
      spectral_wannier_representation(:,:),spectral_spatial_trials(:,:),spectral_amn(:,:)
    complex(8),allocatable::one_shot_hrows(:,:)
    complex(8),allocatable::divided_lcfo_hrows(:,:),divided_lcfo_srows(:,:)
    complex(8),allocatable::dg_hybrid_interface_rows(:,:),dg_hybrid_interface_component_rows(:,:,:)
    complex(8),allocatable::dg_hybrid_final_trace(:,:),dg_hybrid_final_hamiltonian_rows(:,:)
    complex(8),allocatable::dg_hybrid_interior_values(:,:),dg_hybrid_interior_gradients(:,:,:),&
      dg_hybrid_interior_kinetic_action(:,:),dg_hybrid_interior_nonlocal_action(:,:),&
      dg_hybrid_kinetic_rows(:,:),dg_hybrid_zero_local_rows(:,:),dg_hybrid_nonlocal_rows(:,:)
    complex(8),allocatable::divided_basis_representation(:,:,:)
    real(8),allocatable::weights(:),spectrum(:),occupations(:),lcfo_retained_occupations(:),&
      lcfo_retained_eigenvalues(:),local_point_rotations(:,:,:)
    real(8),allocatable::ow_total_density_values(:)
    real(8),allocatable::projector_buffer_real(:,:)
    real(8),allocatable::one_shot_density(:),one_shot_potential(:)
    real(8),allocatable::hybrid_converged_density(:),ow_initial_occupied_density(:)
    real(8),allocatable::divided_initial_density(:),divided_converged_density(:)
    real(8),allocatable::divided_lcfo_point_weights(:)
    real(8),allocatable::dg_hybrid_interior_weights(:),dg_hybrid_unit_local_potential(:)
    real(8),allocatable::dg_hybrid_final_density(:)
    real(8),allocatable::spectral_occupied_density(:),spectral_empty_moments(:,:),&
      spectral_shared_density(:,:),spectral_basin_spectra(:,:),spectral_descriptor_eigenvalues(:),&
      spectral_descriptor_occupations(:)
    real(8),allocatable::ow_raw_partition_weight(:),ow_raw_partition_gradient(:,:),&
      ow_box_density(:)
    real(8),allocatable::localized_centers(:,:),localized_center_magnitudes(:,:)
    real(8),allocatable::w90_fractional(:,:),w90_spreads(:),w90_eigenvalues(:),w90_atoms_cart(:,:),&
      fixed_center_eigenvalues(:),adapted_occupied_spectrum(:)
    real(8),allocatable::translation_w90_values(:),translation_lcfo_values(:),translation_singular_values(:)
    real(8),allocatable::translation_adapted_spectrum(:)
    real(8),allocatable::occupied_density_before(:),occupied_density_after(:),occupied_density_difference(:),&
      occupied_pre_total_residual(:),occupied_pre_boundary_residual(:),occupied_pre_interior_residual(:)
    real(8),allocatable::ow_core_spatial_covariance_residual(:)
    real(8),allocatable::ow_scalar_probe_weights(:),ow_scalar_probe_residual(:)
    real(8),allocatable::ow_gradient_covariance_left(:),ow_gradient_covariance_transpose(:),&
      ow_gradient_covariance_candidates(:,:),ow_gradient_map_commutator(:)
    real(8),allocatable::ow_grid_stencil_defect(:)
    type(t_dg_projection_channel),allocatable::manifest_channels(:)
    type(t_dg_projection_channel),allocatable::projector_tile_channels(:)
    type(s_dg_overlapping_wannier_construction)::symmetry_basis
    type(s_dg_prepared_translation_action)::translation_prepared_action
    type(s_dg_prepared_spectral_basins)::spectral_prepared_basins
    integer(8),allocatable::physical_ids(:),local_symmetry_map(:,:),ow_pencil_generator_maps(:,:),&
      ow_reciprocal_operation_maps(:,:),&
      exact_fragment_symmetry_fingerprints(:),global_symmetry_map(:,:)
    integer(8),allocatable::lcfo_core_ids(:),initial_core_ids(:),ow_total_density_ids(:),&
      ow_neighbor_plus_ids(:),ow_neighbor_minus_ids(:),ow_gradient_identity_map(:,:)
    integer(8),allocatable::all_core_ids(:,:),localized_center_ids(:),orbital_owned_full_ids(:)
    integer(8),allocatable::fixed_center_symmetry_map(:,:),fixed_center_row_ids(:)
    integer(8),allocatable::reindexed_global_symmetry_map(:,:),reindexed_fixed_center_symmetry_map(:,:)
    integer(8),allocatable::translation_row_ids(:),translation_stream_row_ids(:)
    integer(8),allocatable::divided_lcfo_row_ids(:)
    integer,allocatable::divided_effective_ids(:),divided_requested_ids(:),divided_selection_effective_ids(:),&
      divided_added_ids(:),&
      divided_closure_parent(:),divided_closure_reason(:),divided_closure_action(:),divided_scope_selectors(:),&
      divided_basis_owner(:),divided_basis_fragment(:),divided_basis_local_slot(:),divided_basis_generation(:),&
      divided_metric_offsets(:),divided_metric_columns(:),divided_operator_offsets(:),divided_operator_columns(:)
    type(s_dg_hybrid_scope_receipt)::divided_scope_receipt
    type(s_dg_hybrid_production_selection)::divided_production_selection
    integer,allocatable::dg_hybrid_interior_fragment(:)
    real(8)::dg_hybrid_broken_diagnostics(4)
    integer(8),allocatable::translation_spatial_ids(:),translation_generator_maps(:,:)
    integer(8),allocatable::spectral_row_ids(:),spectral_stream_row_ids(:),spectral_complement_row_ids(:)
    integer,allocatable::local_point_product(:,:),local_point_integer_rotations(:,:,:),&
      translation_product(:,:),global_point_product(:,:),global_point_integer_rotations(:,:,:),&
      global_translation_subgroup(:),global_point_representatives(:),global_point_cogroup_product(:,:),&
      global_translation_cocycle(:,:),translation_canonical_product(:,:)
    integer,allocatable::rank_fragments(:)
    integer,allocatable::projector_atom_ids(:)
    integer,allocatable::fixed_center_product(:,:)
    integer,allocatable::global_affine_generators(:)
    integer,allocatable::spectral_basin_labels(:),spectral_basin_generator_maps(:,:),&
      spectral_selected_ranks(:),spectral_block_offsets(:),spectral_orbit_id(:),&
      spectral_orbit_representatives(:),spectral_global_basin_labels(:),spectral_single_basin_map(:,:)
    logical,allocatable::spectral_block_ends(:,:)
    integer,allocatable::translation_canonical_operations(:),translation_inverse_operations(:),&
      translation_character_generators(:),translation_element_words(:,:),translation_character_conjugates(:),&
      translation_generator_orders(:)
    type(t_sawf_symop),allocatable::fixed_center_operations(:)
    type(s_dg_hybrid_fragment_basis),allocatable::divided_fragment_bases(:)
    type(s_dg_hybrid_production_face_trace),allocatable::divided_production_faces(:)
    type(s_dg_hybrid_fixed_payload)::dg_hybrid_fixed_payload
    type(s_dg_hybrid_localization_receipt)::hybrid_localization_receipt
    type(t_sawf_dmn_writer)::fixed_center_dmn_writer
    integer,allocatable::center_owner_candidate(:),center_box_candidate(:),center_fragment_candidate(:)
    integer,allocatable::orbital_owned_ids(:),center_local_orbital_ids(:)
    integer,allocatable::w90_nncell(:,:)
    logical,allocatable::core_mask(:)
    logical,allocatable::translation_character_done(:)
    logical,allocatable::lcfo_boundary_mask(:)
    integer::ix,iy,iz,io,p,nbox,ncore,noccupied,nstate,ntarget,rank,nproc,gradient_distance,&
      raw_ix,raw_iy,raw_iz,core_index,ierr,allocation_status,&
      local_target_count,w90_nntot,projector_tile_first,projector_tile_last,projector_tile_count
    integer::translation_allocation_status
    integer::global_seed_count,global_retained_rank,global_occupied_count,global_projection_count,&
      global_required_retained_rank
    integer::global_identity_operation
    integer::fixed_center_group_order,fixed_center_operation,fixed_center_identity_operation,&
      adapted_occupied_rank,adapted_occupied_selected_block_dimension,orthonormal_lcfo_rank,&
      translation_identity_operation,translation_adapted_rank
    integer::global_point_cogroup_identity_operation
    integer::translation_character_generator_count,translation_sector_rank
    integer::translation_character,translation_partner,translation_global_core_count,translation_processed_count
    integer::translation_point_generator_count,translation_point_checked_pair_count
    integer::spectral_basin_count,spectral_orbit_count,spectral_representative_count,&
      spectral_representative_column,spectral_basin,spectral_orbit,spectral_target_basin,&
      spectral_complement_rank,spectral_complement_local_row
    integer::failed_operation
    integer::w90_replay_environment_status,w90_replay_environment_length,&
      w90_replay_enabled,w90_replay_enabled_min,w90_replay_enabled_max
    character(1024)::w90_replay_directory
    integer::lcfo_symmetry_worst_operation,lcfo_symmetry_worst_generator_index
    real(8),allocatable::lcfo_total_symmetry_residual(:),lcfo_boundary_symmetry_residual(:),&
      lcfo_interior_symmetry_residual(:)
    real(8),allocatable::occupied_affine_total_residual(:),occupied_affine_boundary_residual(:),&
      occupied_affine_interior_residual(:),projection_affine_total_residual(:),&
      projection_affine_boundary_residual(:),projection_affine_interior_residual(:)
    integer::representative_pair(2),local_pair(2)
    integer::complete_sp_core_atom_count
    integer(8)::expected_core_count,expected_box_count,basis_fingerprint,operator_fingerprint,&
      pseudopotential_fingerprint,nbox8,ncore8,product8,nxy8,local_exact_symmetry_fingerprint,&
      lcfo_symmetry_workspace_peak,adapted_occupied_workspace_peak,occupied_pre_closure_workspace_peak
    integer(8)::translation_adapted_workspace_peak
    integer(8)::adapted_occupied_hamiltonian_fingerprint
    integer(8)::translation_character_fingerprint,translation_sector_fingerprint,&
      translation_sector_workspace_peak
    integer(8)::translation_phase_fingerprint,translation_phase_payload_fingerprint,&
      translation_phase_workspace,translation_materialize_fingerprint,translation_materialize_workspace,&
      translation_anchor_fingerprint,translation_anchor_workspace,translation_alignment_fingerprint,&
      translation_alignment_workspace,translation_gamma_fingerprint,translation_gamma_workspace,&
      translation_inverse_workspace,translation_transform_workspace,translation_operator_fingerprint,&
      translation_operator_workspace
    integer(8)::translation_lcfo_fingerprint,translation_post_gauge_fingerprint,translation_global_core_count8
    integer(8)::composition_fingerprint,composition_workspace_peak,occupied_composition_peak,&
      occupied_composition_fingerprint,projector_composition_peak,projector_composition_fingerprint
    integer(8)::w90_coordinator_bytes,w90_workspace_peak,w90_seed_a_workspace,w90_byte_limit
    integer(8)::one_shot_workspace_peak,one_shot_operator_fingerprint
    integer(8)::hybrid_state_workspace,hybrid_state_fingerprint,hybrid_checkpoint_fingerprint,&
      hybrid_scf_fingerprint,hybrid_provenance(6)
    integer(8)::occupied_affine_workspace_peak,projection_affine_workspace_peak
    integer(8)::w90_symmetry_workspace_peak,w90_covariance_workspace
    integer(8)::center_gauge_workspace_peak
    integer(8)::fixed_center_group_fingerprint,fixed_center_operation_workspace,&
      fixed_center_dmn_workspace_peak
    integer(8)::w90_input_fingerprint,w90_transform_fingerprint,hybrid_localization_seed_fingerprint
    integer(8)::spectral_frame_fingerprint,spectral_density_fingerprint,spectral_basin_fingerprint,&
      spectral_operator_fingerprint,spectral_eigensystem_fingerprint,spectral_catalog_fingerprint,&
      spectral_channel_fingerprint,spectral_complement_channel_fingerprint,&
      spectral_action_fingerprint,spectral_action_aggregate_fingerprint,&
      spectral_workspace_peak,&
      spectral_operation_workspace
    integer(8)::ow_stitched_peak_elements,ow_density_redistribution_workspace
    integer(8)::divided_pw_fingerprint,divided_buffer_window_fingerprint,divided_fragment_fingerprint,&
      divided_lcfo_peak_elements,divided_lcfo_operator_fingerprint,divided_state_workspace,&
      divided_state_fingerprint,divided_final_solver_workspace,divided_final_solver_fingerprint,&
      divided_selection_fingerprint,divided_fixed_payload_fingerprint,divided_basis_directory_fingerprint
    real(8)::condition_number,closure_residual,spread_max,gauge_correction
    real(8)::adapted_occupied_trace,adapted_occupied_closure,adapted_occupied_gamma_defect,&
      translation_adapted_trace,translation_adapted_closure,translation_adapted_gamma_defect,&
      adapted_occupied_selected_edge,adapted_occupied_rejected_edge,adapted_occupied_cluster_gap
    real(8)::adapted_occupied_secondary_selected_edge,adapted_occupied_secondary_rejected_edge,&
      adapted_occupied_secondary_cluster_gap,adapted_occupied_secondary_residual
    real(8)::adapted_occupied_subspace_distance,adapted_occupied_density_interior_difference,&
      adapted_occupied_density_boundary_difference,local_occupied_density_interior_difference,&
      local_occupied_density_boundary_difference,adapted_occupied_closure_before,&
      local_occupied_density_interior_norm,local_occupied_density_boundary_norm,&
      global_occupied_density_interior_norm,global_occupied_density_boundary_norm
    real(8)::ow_partition_sum_defect,ow_partition_gradient_defect,window_axis(3),&
      window_axis_derivative(3),window_coordinate
    real(8)::ow_stitched_electron_count,ow_stitched_s_hermiticity,ow_stitched_rho_hermiticity
    real(8)::ow_stitched_minimum_pivot,ow_stitched_pivot_condition
    real(8)::one_shot_residual,one_shot_orthogonality,one_shot_condition,&
      one_shot_gamma_defect,one_shot_charge,one_shot_trace_charge,&
      one_shot_local_difference,one_shot_global_difference,one_shot_local_norm,one_shot_global_norm
    real(8)::hybrid_density_residual,hybrid_energy_residual,hybrid_eigensystem_residual,&
      hybrid_electron_defect,hybrid_symmetry_defect,hybrid_scf_receipts(5)
    real(8)::initial_occupied_charge_local,initial_occupied_charge
    integer::hybrid_iterations
    integer::divided_iterations,dg_hybrid_nonlocal_ownership_count,divided_global_basis_count
    real(8)::divided_convergence_value,divided_electron_defect
    real(8)::divided_final_residual,divided_final_orthogonality,divided_final_projector_defect
    real(8)::divided_full_basis_closure_defect
    logical::ok,reusable,localization_converged,global_inversion_present,center_diagnostic_ok,diagnostic_ok
    logical::fixed_center_inversion_present,writer_ok
    logical::translation_self_conjugate
    logical::divided_full_basis_closure_ok
    real(8)::fixed_center_fractional(3)
    complex(8),allocatable::core_periodic_phase(:,:),localization_transform(:,:),retained_identity(:,:)
    complex(8),allocatable::synchronized_local_representation(:,:,:)
    real(8)::localization_initial_spread,localization_final_spread,localization_maximum_gradient,&
      retained_raw_unitarity_defect,retained_unitarity_defect,retained_group_closure_defect,&
      global_retained_group_closure_defect,retained_closure_search_tolerance
    real(8),allocatable::global_point_rotations(:,:,:)
    real(8)::ow_gradient_identity_rotation(3,3,1)
    real(8),allocatable::global_point_fractional_translations(:,:)
    real(8)::w90_reciprocal_lattice(3,3),w90_lattice_inverse(3,3),w90_determinant,w90_spread(3)
    real(8)::w90_identity_defect,w90_unitarity_defect,w90_closure_defect,w90_covariance_defect,&
      w90_generator_covariance_defect
    real(8)::translation_identity_defect,translation_unitarity_defect,translation_commutator_defect,&
      translation_order_defect,translation_gamma_pairing_defect
    real(8)::translation_operator_defect,translation_anchor_defect,translation_alignment_defect,&
      translation_gamma_defect,translation_closure_defect
    real(8)::spectral_frame_defect,spectral_operator_hermiticity,spectral_operator_trace,&
      spectral_eigensystem_residual,spectral_channel_gram_defect,spectral_action_unitarity,&
      spectral_action_block_defect
    real(8)::translation_alignment_max_defect,translation_gamma_max_defect
    real(8)::ow_gradient_path_local(2),ow_gradient_path_global(2),&
      ow_spatial_covariance_relative,ow_spatial_covariance_absolute,&
      ow_gradient_covariance_absolute,ow_gradient_stencil_norm_bound,ow_stencil_axis_defect
    real(8)::ow_map_probe_angle,ow_map_probe_symbol
    real(8)::monomial_defect,center_block_leakage,center_representation_unitarity_defect
    type(s_dg_translation_orbit_accumulator)::translation_inverse_state
    integer::localization_iterations,localization_spread_evaluations
    integer::ow_saved_eigenexa_comm
    character(256)::message,prefix,center_failure_message,center_diagnostic_message,diagnostic_message
    character(8),allocatable::w90_atom_symbols(:)

    call MPI_Comm_rank(dc%icomm_tot,rank,ierr);call MPI_Comm_size(dc%icomm_tot,nproc,ierr)
    ok=system%nspin==1.and.system%if_real_orbital.and.allocated(spsi%rwf)
    call comm_logical_and(ok,reusable,dc%icomm_tot)
    if(.not.reusable)error stop 'overlapping-Wannier production requires Gamma real DC candidates'
    ok=nproc==dc%n_frag.and..not.dc%optimized_fragment_geometry
    call comm_logical_and(ok,reusable,dc%icomm_tot)
    if(.not.reusable)&
      error stop 'overlapping-Wannier production requires one rank per valid DC fragment'
    ow_core_size=dc%nxyz_domain_frag(:,dc%i_frag);ow_buffer=dc%nxyz_buffer
    do ix=1,3
      if(ow_buffer(ix)>huge(ow_box_size(ix))/2)&
        error stop 'overlapping-Wannier buffer extent overflow'
      if(ow_core_size(ix)>huge(ow_box_size(ix))-2*ow_buffer(ix))&
        error stop 'overlapping-Wannier buffer extent overflow'
    enddo
    ow_box_size=ow_core_size+2*ow_buffer
    ok=all(ow_buffer>=size(stencil%coef_nab,1).or.ow_core_size==dc%lg_tot%num)
    call comm_logical_and(ok,reusable,dc%icomm_tot)
    if(.not.reusable)error stop 'overlapping-Wannier buffer is smaller than the SALMON gradient stencil'
    ok=.not.any([size(spsi%rwf,1),size(spsi%rwf,2),size(spsi%rwf,3)]<ow_box_size)
    call comm_logical_and(ok,reusable,dc%icomm_tot)
    if(.not.reusable)&
      error stop 'overlapping-Wannier candidate buffer box is incomplete'
    call checked_ow_extent_product(ow_box_size,nbox8,ok)
    call checked_ow_extent_product(ow_core_size,ncore8,reusable);ok=ok.and.reusable
    if(.not.ok.or.nbox8>int(huge(nbox),8).or.ncore8>int(huge(ncore),8))&
      error stop 'overlapping-Wannier grid extent overflow'
    nbox=int(nbox8);ncore=int(ncore8)
    nstate=ceiling(0.5d0*dc%elec_num_tot)
    if(mod(nstate,nproc)/=0)error stop 'LCFO occupied rank is not rank balanced'
    noccupied=nstate/nproc
    call checked_ow_extent_product(dc%lg_tot%num,expected_core_count,ok)
    if(.not.ok.or.expected_core_count>int(huge(nbox),8))&
      error stop 'overlapping-Wannier global grid exceeds addressable extent'
    ow_global_grid_count=expected_core_count
    nxy8=int(dc%lg_tot%num(1),8)*int(dc%lg_tot%num(2),8)
    if(.not.ok.or.nbox8>huge(expected_box_count)/int(nproc,8))&
      error stop 'overlapping-Wannier global extent overflow'
    expected_box_count=nbox8*int(nproc,8)
    allocate(weights(nbox),ow_box_density(nbox),physical_ids(nbox),core_mask(nbox),stat=allocation_status)
    call comm_logical_and(allocation_status==0,reusable,dc%icomm_tot)
    if(.not.reusable)error stop 'overlapping-Wannier production allocation failed'
    ok=allocated(rho_s)
    if(ok)ok=size(rho_s)>=1
    if(ok)ok=all(lbound(rho_s(1)%f)<=[1,1,1]).and.all(ubound(rho_s(1)%f)>=ow_box_size)
    call comm_logical_and(ok,reusable,dc%icomm_tot)
    if(.not.reusable)error stop 'overlapping-Wannier fragment density buffer is incomplete'
    weights=system%hvol
    p=0;core_index=0
    do iz=1,ow_box_size(3);do iy=1,ow_box_size(2);do ix=1,ow_box_size(1)
      p=p+1
      core_mask(p)=ix>ow_buffer(1).and.ix<=ow_buffer(1)+ow_core_size(1).and.&
        iy>ow_buffer(2).and.iy<=ow_buffer(2)+ow_core_size(2).and.&
        iz>ow_buffer(3).and.iz<=ow_buffer(3)+ow_core_size(3)
      physical_ids(p)=1_8+int(modulo(dc%ixyz_frag(1,dc%i_frag)-1+ix-ow_buffer(1)-1,dc%lg_tot%num(1)),8)+&
        int(dc%lg_tot%num(1),8)*(int(modulo(dc%ixyz_frag(2,dc%i_frag)-1+iy-ow_buffer(2)-1,&
        dc%lg_tot%num(2)),8)+int(dc%lg_tot%num(2),8)*int(modulo(dc%ixyz_frag(3,dc%i_frag)-1+&
        iz-ow_buffer(3)-1,dc%lg_tot%num(3)),8))
      if(core_mask(p))then
        core_index=core_index+1
      endif
    enddo;enddo;enddo
    allocate(ow_box_physical_ids,source=physical_ids)
    allocate(ow_raw_partition_weight(nbox),ow_raw_partition_gradient(3,nbox),&
      ow_partition_weight(nbox),ow_partition_gradient(3,nbox))
    do p=1,nbox
      raw_ix=modulo(p-1,ow_box_size(1))+1
      raw_iy=modulo((p-1)/ow_box_size(1),ow_box_size(2))+1
      raw_iz=(p-1)/(ow_box_size(1)*ow_box_size(2))+1
      do ix=1,3
        select case(ix)
        case(1);iy=raw_ix
        case(2);iy=raw_iy
        case default;iy=raw_iz
        end select
        if(ow_buffer(ix)==0.or.(iy>ow_buffer(ix).and.iy<=ow_buffer(ix)+ow_core_size(ix)))then
          window_axis(ix)=1d0;window_axis_derivative(ix)=0d0
        elseif(iy<=ow_buffer(ix))then
          window_coordinate=real(iy,8)/real(ow_buffer(ix)+1,8)
          window_axis(ix)=window_coordinate**2*(3d0-2d0*window_coordinate)
          window_axis_derivative(ix)=6d0*window_coordinate*(1d0-window_coordinate)/&
            (real(ow_buffer(ix)+1,8)*system%hgs(ix))
        else
          window_coordinate=real(ow_box_size(ix)+1-iy,8)/real(ow_buffer(ix)+1,8)
          window_axis(ix)=window_coordinate**2*(3d0-2d0*window_coordinate)
          window_axis_derivative(ix)=-6d0*window_coordinate*(1d0-window_coordinate)/&
            (real(ow_buffer(ix)+1,8)*system%hgs(ix))
        endif
      enddo
      ow_raw_partition_weight(p)=product(window_axis)
      ow_raw_partition_gradient(1,p)=window_axis_derivative(1)*window_axis(2)*window_axis(3)
      ow_raw_partition_gradient(2,p)=window_axis(1)*window_axis_derivative(2)*window_axis(3)
      ow_raw_partition_gradient(3,p)=window_axis(1)*window_axis(2)*window_axis_derivative(3)
    enddo
    call build_dg_smooth_partition_of_unity(dc%icomm_tot,physical_ids,ow_raw_partition_weight,&
      ow_raw_partition_gradient,ow_partition_weight,ow_partition_gradient,ow_partition_sum_defect,&
      ow_partition_gradient_defect,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'overlapping-Wannier smooth partition failed';endif
    if(rank==0)write(*,'(a,2(a,es16.8))')'[OW-GS-DIAGNOSTIC] smooth_partition',&
      ' sum_defect=',ow_partition_sum_defect,' gradient_defect=',ow_partition_gradient_defect
    ncore8=int(dc%mg_tot%ie(1)-dc%mg_tot%is(1)+1,8)*&
      int(dc%mg_tot%ie(2)-dc%mg_tot%is(2)+1,8)*int(dc%mg_tot%ie(3)-dc%mg_tot%is(3)+1,8)
    if(ncore8<1_8.or.ncore8>int(huge(ncore),8))&
      error stop 'distributed total-density slab extent overflow'
    allocate(ow_total_density_ids(int(ncore8)),ow_total_density_values(int(ncore8)),&
      stat=allocation_status)
    call comm_logical_and(allocation_status==0,reusable,dc%icomm_tot)
    if(.not.reusable)error stop 'distributed total-density slab allocation failed'
    core_index=0
    do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)
    do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
    do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
      core_index=core_index+1
      ow_total_density_ids(core_index)=int(ix,8)+int(dc%lg_tot%num(1),8)*&
        (int(iy-1,8)+int(dc%lg_tot%num(2),8)*int(iz-1,8))
      ow_total_density_values(core_index)=dc%rho_tot_s(1)%f(ix,iy,iz)
    enddo
    enddo
    enddo
    call redistribute_dg_row_owned_real_field_to_requests(dc%icomm_tot,expected_core_count,&
      ow_total_density_ids,ow_total_density_values,physical_ids,ow_box_density,&
      ow_density_redistribution_workspace,ok,message)
    if(.not.ok)then
      write(0,'(a)')trim(message);error stop 'distributed total-density buffer materialization failed'
    endif
    if(rank==0)write(*,'(a,i0)')'[OW-GS-DIAGNOSTIC] total_density_buffer_workspace_peak_bytes=',&
      ow_density_redistribution_workspace
    pseudopotential_fingerprint=ow_collective_operator_fingerprint(dc%icomm_tot)
    allocate(lcfo_core_ids(ncore),lcfo_boundary_mask(ncore));core_index=0
    do p=1,nbox
      if(.not.core_mask(p))cycle
      core_index=core_index+1;lcfo_core_ids(core_index)=physical_ids(p)
      ix=modulo(p-1,ow_box_size(1))+1
      iy=modulo((p-1)/ow_box_size(1),ow_box_size(2))+1
      iz=(p-1)/(ow_box_size(1)*ow_box_size(2))+1
      lcfo_boundary_mask(core_index)=ix-ow_buffer(1)<=size(stencil%coef_nab,1).or.&
        ix-ow_buffer(1)>ow_core_size(1)-size(stencil%coef_nab,1).or.&
        iy-ow_buffer(2)<=size(stencil%coef_nab,1).or.&
        iy-ow_buffer(2)>ow_core_size(2)-size(stencil%coef_nab,1).or.&
        iz-ow_buffer(3)<=size(stencil%coef_nab,1).or.&
        iz-ow_buffer(3)>ow_core_size(3)-size(stencil%coef_nab,1)
    end do
    allocate(projector_atom_ids(dc%system_tot%nion));projector_atom_ids=[(i,i=1,dc%system_tot%nion)]
    call build_dg_complete_sp_manifest(projector_atom_ids,manifest_channels,ok,message)
    deallocate(projector_atom_ids)
    if(.not.ok)then;write(0,'(a)')trim(message)
      error stop 'overlapping-Wannier complete-s+p projector catalog failed';endif
    global_projection_count=size(manifest_channels)
    if(mod(global_projection_count,nproc)/=0.or.mod(dc%system_tot%nion,nproc)/=0)&
      error stop 'global complete-s+p projector catalog is not rank balanced'
    local_target_count=global_projection_count/nproc
    if(nstate>huge(ntarget)-global_projection_count)error stop 'LCFO Wannier target rank overflow'
    ntarget=nstate+global_projection_count
    if(ntarget<1)error stop 'LCFO Wannier target rank is invalid'
    call dc_lcfo(lg,mg,system,info,stencil,ppg,energy,v_local,spsi,shpsi,sttpsi,srg,dc,&
      retained_count=ntarget,retained_box_count=nstate,&
      retained_box_contribution=lcfo_fragment_contribution,&
      retained_occupations=lcfo_retained_occupations,retained_eigenvalues=lcfo_retained_eigenvalues,&
      write_files=.false.)
    if(size(lcfo_retained_occupations)/=ntarget.or.&
        abs(sum(lcfo_retained_occupations)-dc%elec_num_tot)>&
        1d3*epsilon(1d0)*max(1d0,dc%elec_num_tot)) &
      error stop 'LCFO retained occupations do not match the target space'
    call compose_dg_buffered_orbital_tile_to_physical_grid(dc%icomm_tot,physical_ids,&
      ow_partition_weight,lcfo_fragment_contribution,ow_core_ids,global_seed_values,&
      composition_fingerprint,composition_workspace_peak,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message)
      error stop 'LCFO occupied buffer composition failed';endif
    occupied_composition_peak=composition_workspace_peak
    occupied_composition_fingerprint=composition_fingerprint
    deallocate(lcfo_fragment_contribution)
    if(size(ow_core_ids)/=ncore.or.any(shape(global_seed_values)/=[nstate,ncore]))&
      error stop 'LCFO occupied buffer composition has unexpected ownership shape'
    allocate(lcfo_occupied_core(nstate,ncore),source=global_seed_values)
    if(int(ncore,8)>huge(0_8)/int(nproc,8))error stop 'spectral global core extent overflows int64'
    translation_global_core_count8=int(ncore,8)*int(nproc,8)
    if(translation_global_core_count8>int(huge(0),8))&
      error stop 'spectral global core extent exceeds default integer'
    translation_global_core_count=int(translation_global_core_count8)
    complete_sp_core_atom_count=dc%system_tot%nion/nproc
    if(complete_sp_core_atom_count<1.or.&
        local_target_count/=4*complete_sp_core_atom_count)&
      error stop 'complete-s+p target is not four channels per core-owned atom'
    if(rank==0)write(*,'(a,i0,a,i0)')&
      '[OW-GS-DIAGNOSTIC] complete_sp_core_atom_count=',complete_sp_core_atom_count,&
      ' complete_sp_shell_channels=',local_target_count
    global_occupied_count=nstate
    retained_closure_search_tolerance=sqrt(sqrt(dg_ow_symmetry_tolerance))
    if(rank==0)write(*,'(a,i0,a,i0,a,i0)')&
      '[OW-GS-DIAGNOSTIC] lcfo_occupied_rank=',nstate,&
      ' lcfo_projection_localizer_rank=',global_projection_count,' lcfo_target_rank=',ntarget
    allocate(ow_basis%center_box_point_ids(ntarget),ow_basis%center_owner_rank(ntarget),&
      ow_basis%center_owner_fragment(ntarget))
    ow_basis%center_box_point_ids=0_8;ow_basis%center_owner_rank=-1
    ow_basis%center_owner_fragment=-1
    global_seed_count=ntarget
    call move_alloc(global_seed_values,composed_tile_values)
    allocate(global_seed_values(global_seed_count,ncore),ow_core_weights(ncore),&
      ow_core_box_positions(ncore),core_periodic_phase(3,ncore));global_seed_values=(0d0,0d0)
    global_seed_values(1:nstate,:)=composed_tile_values(1:nstate,:);deallocate(composed_tile_values)
    ow_core_weights=system%hvol;ow_core_box_positions=0;core_periodic_phase=(0d0,0d0)
    projector_composition_peak=0_8;projector_composition_fingerprint=0_8
    do projector_tile_first=1,global_projection_count,32
      projector_tile_last=min(global_projection_count,projector_tile_first+31)
      projector_tile_count=projector_tile_last-projector_tile_first+1
      call build_ow_complete_sp_projectors(physical_ids,pseudopotential_fingerprint,&
        projector_tile_channels,projector_buffer_real,ok,message,&
        manifest_channels(projector_tile_first:projector_tile_last))
      if(.not.ok)then;write(0,'(a)')trim(message)
        error stop 'overlapping-Wannier complete-s+p projector tile failed';endif
      allocate(projector_buffer_tile(projector_tile_count,nbox))
      projector_buffer_tile=cmplx(projector_buffer_real,0d0,8)
      deallocate(projector_buffer_real,projector_tile_channels)
      call compose_dg_buffered_orbital_tile_to_physical_grid(dc%icomm_tot,physical_ids,&
        ow_partition_weight,projector_buffer_tile,lcfo_core_ids,composed_tile_values,&
        composition_fingerprint,composition_workspace_peak,ok,message)
      deallocate(projector_buffer_tile)
      if(.not.ok.or.any(lcfo_core_ids/=ow_core_ids).or.&
          any(shape(composed_tile_values)/=[projector_tile_count,ncore]))then
        write(0,'(a)')trim(message);error stop 'complete-s+p buffer composition failed'
      endif
      global_seed_values(nstate+projector_tile_first:nstate+projector_tile_last,:)=composed_tile_values
      projector_composition_peak=max(projector_composition_peak,composition_workspace_peak)
      projector_composition_fingerprint=ieor(projector_composition_fingerprint,&
        ishftc(composition_fingerprint,mod(projector_tile_first,63)))
      deallocate(composed_tile_values)
    enddo
    ! This is the physical source frame for the LCFO discriminator: actual retained
    ! LCFO eigenfunctions for the occupied block and the actual radial manifest
    ! projectors for the complement.  It need not be orthonormal; projecting the
    ! weighted rank-one sum is the operator definition.
    allocate(lcfo_reference_core(ntarget,ncore))
    lcfo_reference_core(1:nstate,:)=lcfo_occupied_core
    lcfo_reference_core(nstate+1:ntarget,:)=global_seed_values(nstate+1:ntarget,:)
    do p=1,ncore
      raw_ix=int(modulo(ow_core_ids(p)-1_8,int(dc%lg_tot%num(1),8)))+1
      raw_iy=int(modulo((ow_core_ids(p)-1_8)/int(dc%lg_tot%num(1),8),&
        int(dc%lg_tot%num(2),8)))+1
      raw_iz=int((ow_core_ids(p)-1_8)/nxy8)+1
      lcfo_boundary_mask(p)=modulo(raw_ix-1,ow_core_size(1))<size(stencil%coef_nab,1).or.&
        modulo(raw_ix-1,ow_core_size(1))>=ow_core_size(1)-size(stencil%coef_nab,1).or.&
        modulo(raw_iy-1,ow_core_size(2))<size(stencil%coef_nab,1).or.&
        modulo(raw_iy-1,ow_core_size(2))>=ow_core_size(2)-size(stencil%coef_nab,1).or.&
        modulo(raw_iz-1,ow_core_size(3))<size(stencil%coef_nab,1).or.&
        modulo(raw_iz-1,ow_core_size(3))>=ow_core_size(3)-size(stencil%coef_nab,1)
    enddo
    if(rank==0)write(*,'(a,4(a,i0))')'[OW-GS-DIAGNOSTIC] buffer_composition',&
      ' occupied_workspace_peak_bytes=',occupied_composition_peak,&
      ' projector_workspace_peak_bytes=',projector_composition_peak,&
      ' occupied_fingerprint=',occupied_composition_fingerprint,&
      ' projector_fingerprint=',projector_composition_fingerprint
    call prepare_ow_global_point_action(ow_core_ids,global_symmetry_map,global_point_integer_rotations,&
      global_point_rotations,global_point_fractional_translations,global_point_product,&
      global_translation_subgroup,global_point_representatives,global_point_cogroup_product,&
      global_translation_cocycle,&
      global_inversion_present,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'global point-action construction failed';end if
    ! HYBRID_LOCALIZATION_FIRST_BRANCH_BEGIN
    if(yn_dg_hybrid_continuation_scf=='y')then
    ! HYBRID_LOCALIZATION_FIRST_ARM_BEGIN
    call find_dg_group_identity(global_point_product,global_identity_operation,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'global group identity construction failed';endif
    call find_dg_group_identity(global_point_cogroup_product,global_point_cogroup_identity_operation,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'point cogroup identity construction failed';endif
    call select_dg_group_generators(global_point_product,global_identity_operation,&
      global_affine_generators,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'global affine generator selection failed';endif
    if(rank==0)write(*,'(a,2(a,i0))')'[OW-GS-DIAGNOSTIC] affine_generator_proof',&
      ' group_order=',size(global_point_product,1),' generator_count=',size(global_affine_generators)
    if(size(global_affine_generators)>0)then
      call measure_dg_grid_map_stencil_defect(dc%icomm_tot,ow_core_ids,&
        global_symmetry_map(:,global_affine_generators),dc%lg_tot%num,system%hgs,&
        global_point_rotations(:,:,global_affine_generators),ow_grid_stencil_defect,ok,message)
      if(.not.ok)then;write(0,'(a)')trim(message);error stop 'affine grid-stencil diagnostic failed';endif
      if(rank==0)write(*,'(a,es16.8,a,i0)')'[OW-GS-DIAGNOSTIC] affine grid-stencil defect max=',&
        maxval(ow_grid_stencil_defect),' operation=',maxloc(ow_grid_stencil_defect,dim=1)
    else
      allocate(ow_grid_stencil_defect(0))
      if(rank==0)write(*,'(a,es16.8,a,i0)')'[OW-GS-DIAGNOSTIC] affine grid-stencil defect max=',&
        0d0,' operation=',global_identity_operation
    endif
    deallocate(ow_grid_stencil_defect)

    call prepare_dg_hybrid_localization_first_seed(dc%icomm_tot,ow_core_ids,global_seed_values,&
      ow_core_weights,dg_dc_metric_rank_tolerance,global_closed_core,global_retained_rank,&
      hybrid_localization_seed_fingerprint,ok,message)
    if(.not.ok.or.global_retained_rank/=ntarget)then
      if(.not.ok)write(0,'(a)')trim(message)
      error stop 'Hybrid localization-first occupied+s+p seed preparation failed'
    endif
    global_required_retained_rank=ntarget
    allocate(lcfo_total_symmetry_residual(size(global_affine_generators)),&
      lcfo_boundary_symmetry_residual(size(global_affine_generators)),&
      lcfo_interior_symmetry_residual(size(global_affine_generators)))
    if(size(global_affine_generators)>0)then
      call measure_dg_rank_fixed_symmetry_residuals(dc%icomm_tot,global_closed_core,ow_core_weights,&
        global_symmetry_map(:,global_affine_generators),lcfo_boundary_mask,&
        total_residual=lcfo_total_symmetry_residual,boundary_residual=lcfo_boundary_symmetry_residual,&
        interior_residual=lcfo_interior_symmetry_residual,ok=ok,message=message,&
        workspace_peak_bytes=lcfo_symmetry_workspace_peak)
      if(.not.ok)then
        write(0,'(a)')trim(message);error stop 'Hybrid construction-basis symmetry diagnostic failed'
      endif
      lcfo_symmetry_worst_generator_index=maxloc(lcfo_total_symmetry_residual,dim=1)
      lcfo_symmetry_worst_operation=global_affine_generators(lcfo_symmetry_worst_generator_index)
      global_retained_group_closure_defect=maxval(lcfo_total_symmetry_residual)
      if(rank==0)write(*,'(a,i0,3(a,es16.8))')&
        '[HYBRID-WF-SYMMETRY-DIAGNOSTIC] complete_basis_worst_operation=',lcfo_symmetry_worst_operation,&
        ' total=',lcfo_total_symmetry_residual(lcfo_symmetry_worst_generator_index),&
        ' boundary=',lcfo_boundary_symmetry_residual(lcfo_symmetry_worst_generator_index),&
        ' interior=',lcfo_interior_symmetry_residual(lcfo_symmetry_worst_generator_index)
    else
      lcfo_symmetry_worst_generator_index=0;lcfo_symmetry_worst_operation=global_identity_operation
      global_retained_group_closure_defect=0d0;lcfo_symmetry_workspace_peak=0_8
      if(rank==0)write(*,'(a,i0,3(a,es16.8))')&
        '[HYBRID-WF-SYMMETRY-DIAGNOSTIC] complete_basis_worst_operation=',lcfo_symmetry_worst_operation,&
        ' total=',0d0,' boundary=',0d0,' interior=',0d0
    endif

    allocate(w90_anchors(ntarget,ncore),source=global_seed_values,stat=allocation_status)
    call MPI_Allreduce(MPI_IN_PLACE,allocation_status,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_status/=0)&
      error stop 'Hybrid localization-first seed-anchor retention failed collectively'
    w90_byte_limit=8_8*1024_8*1024_8*1024_8
    call assemble_dg_w90_gamma_a_matrix(dc%icomm_tot,global_closed_core,w90_anchors,&
      ow_core_weights,dg_ow_symmetry_tolerance,w90_byte_limit,w90_seed_a_matrix,&
      w90_seed_a_workspace,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'Hybrid localization-first A assembly failed';endif
    deallocate(global_seed_values)
    if(allocated(lcfo_occupied_core))deallocate(lcfo_occupied_core)

    allocate(initial_core_ids(ncore))
    core_index=0
    do p=1,nbox
      if(.not.core_mask(p))cycle
      core_index=core_index+1
      if(core_index>ncore)error stop 'Hybrid localization-first core row extent overflow'
      initial_core_ids(core_index)=physical_ids(p)
      ow_core_weights(core_index)=weights(p)
      ow_core_box_positions(core_index)=p
      core_periodic_phase(1,core_index)=exp(cmplx(0d0,2d0*pi*real(modulo(physical_ids(p)-1_8,&
        int(dc%lg_tot%num(1),8)),8)/real(dc%lg_tot%num(1),8),8))
      core_periodic_phase(2,core_index)=exp(cmplx(0d0,2d0*pi*real(modulo((physical_ids(p)-1_8)/&
        int(dc%lg_tot%num(1),8),int(dc%lg_tot%num(2),8)),8)/real(dc%lg_tot%num(2),8),8))
      core_periodic_phase(3,core_index)=exp(cmplx(0d0,2d0*pi*real((physical_ids(p)-1_8)/&
        nxy8,8)/real(dc%lg_tot%num(3),8),8))
    enddo
    if(core_index/=ncore.or.any(ow_core_box_positions<1))&
      error stop 'Hybrid localization-first core box map is incomplete'
    call materialize_ow_distributed_core_to_buffer(dc%icomm_tot,global_closed_core,ow_core_ids,&
      initial_core_ids,ow_core_values,ok,message)
    if(.not.ok)then
      write(0,'(a)')trim(message)
      error stop 'Hybrid localization-first retained core redistribution failed'
    endif
    call reindex_dg_point_maps_between_row_layouts(dc%icomm_tot,ow_core_ids,initial_core_ids,&
      global_symmetry_map,reindexed_global_symmetry_map,ok,message)
    if(.not.ok)then
      write(0,'(a)')trim(message)
      error stop 'Hybrid localization-first symmetry-map row-layout transition failed'
    endif
    call move_alloc(reindexed_global_symmetry_map,global_symmetry_map)
    global_closed_core=ow_core_values
    ow_core_ids=initial_core_ids
    deallocate(initial_core_ids)
    call invert_ow_lattice(dc%system_tot%primitive_a,w90_lattice_inverse,w90_determinant,ok)
    if(.not.ok)error stop 'Wannier90 lattice is singular'
    w90_reciprocal_lattice=2d0*pi*transpose(w90_lattice_inverse)
    allocate(w90_atom_symbols(dc%system_tot%nion),w90_atoms_cart(3,dc%system_tot%nion))
    w90_atoms_cart=dc%system_tot%Rion
    do io=1,dc%system_tot%nion
      if(dc%system_tot%kion(io)<1.or.dc%system_tot%kion(io)>size(pp%atom_symbol))&
        error stop 'Wannier90 atom species is outside the pseudopotential table'
      w90_atom_symbols(io)=pp%atom_symbol(dc%system_tot%kion(io))
    enddo
    call setup_dg_w90_gamma_library(dc%icomm_tot,'overlapping_wannier_mlwf',&
      dc%system_tot%primitive_a,w90_reciprocal_lattice,w90_atom_symbols,w90_atoms_cart,&
      ntarget,ntarget,wannier_num_iter,dg_ow_w90_initial_projection,DG_W90_UNCONSTRAINED,&
      w90_nntot,w90_nncell,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'Wannier90 Gamma setup failed';endif
    allocate(w90_fractional(3,ncore),w90_eigenvalues(ntarget))
    do p=1,ncore
      w90_fractional(:,p)=[real(modulo(ow_core_ids(p)-1_8,int(dc%lg_tot%num(1),8)),8)/&
        real(dc%lg_tot%num(1),8),real(modulo((ow_core_ids(p)-1_8)/int(dc%lg_tot%num(1),8),&
        int(dc%lg_tot%num(2),8)),8)/real(dc%lg_tot%num(2),8),real((ow_core_ids(p)-1_8)/nxy8,8)/&
        real(dc%lg_tot%num(3),8)]
    enddo
    w90_eigenvalues=0d0
    call assemble_dg_w90_gamma_matrices(dc%icomm_tot,global_closed_core,w90_anchors,&
      ow_core_weights,w90_fractional,w90_nncell,w90_byte_limit,w90_m_matrix,w90_a_matrix,&
      w90_coordinator_bytes,w90_workspace_peak,ok,message,precomputed_a_matrix=w90_seed_a_matrix)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'Wannier90 M/A assembly failed';endif
    deallocate(w90_seed_a_matrix)
    call fingerprint_ow_w90_matrices(dc%icomm_tot,w90_m_matrix,w90_a_matrix,&
      w90_input_fingerprint,ok)
    if(.not.ok)error stop 'Wannier90 M/A fingerprint failed'
    w90_input_fingerprint=ieor(w90_input_fingerprint,hybrid_localization_seed_fingerprint)
    if(w90_input_fingerprint==0_8)w90_input_fingerprint=1_8
    w90_replay_directory=''
    call get_environment_variable('SALMON_DG_W90_REPLAY_DIRECTORY',w90_replay_directory,&
      length=w90_replay_environment_length,status=w90_replay_environment_status,trim_name=.true.)
    w90_replay_enabled=merge(1,0,w90_replay_environment_status==0.and.w90_replay_environment_length>0)
    call MPI_Allreduce(w90_replay_enabled,w90_replay_enabled_min,1,MPI_INTEGER,MPI_MIN,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)error stop 'Wannier90 replay enable MIN agreement failed'
    call MPI_Allreduce(w90_replay_enabled,w90_replay_enabled_max,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.w90_replay_enabled_min/=w90_replay_enabled_max)&
      error stop 'Wannier90 replay enable state disagrees across ranks'
    if(w90_replay_enabled==1)then
      call export_dg_w90_replay_bundle(dc%icomm_tot,'.','overlapping_wannier_mlwf',&
        trim(w90_replay_directory),'overlapping_wannier_mlwf',DG_W90_UNCONSTRAINED,&
        w90_eigenvalues,w90_a_matrix,w90_m_matrix,w90_nncell,ok,message)
      if(.not.ok)then;write(0,'(a)')trim(message);error stop 'Wannier90 replay export failed';endif
    endif
    deallocate(w90_anchors,w90_fractional)
    call run_dg_w90_gamma_library(dc%icomm_tot,'overlapping_wannier_mlwf',&
      dc%system_tot%primitive_a,w90_reciprocal_lattice,w90_atom_symbols,w90_atoms_cart,&
      w90_m_matrix,w90_a_matrix,w90_eigenvalues,huge(1d0)/4d0,dg_ow_symmetry_tolerance,&
      wannier_num_iter,w90_transform,localized_centers,w90_spreads,w90_spread,ok,message,&
      localization_iterations)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'Wannier90 MLWF optimization failed';endif
    deallocate(w90_m_matrix,w90_a_matrix,w90_eigenvalues,w90_atom_symbols,w90_atoms_cart,w90_nncell)
    localized_centers=matmul(w90_lattice_inverse,localized_centers)
    call apply_dg_w90_gamma_transform(dc%icomm_tot,ow_core_ids,ow_core_values,&
      transform=w90_transform,centers=localized_centers,tolerance=dg_ow_symmetry_tolerance,&
      ok=ok,message=message,spreads=w90_spreads)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'Wannier90 MLWF gauge canonicalization failed';endif
    call fingerprint_ow_w90_transform(dc%icomm_tot,w90_transform,w90_transform_fingerprint,ok)
    if(.not.ok)error stop 'Wannier90 canonical transform fingerprint failed'
    call build_dg_hybrid_localization_receipt(dc%icomm_tot,ntarget,global_retained_rank,&
      hybrid_localization_seed_fingerprint,w90_transform,localized_centers,w90_spreads,.true.,&
      localization_iterations,dg_ow_symmetry_tolerance,hybrid_localization_receipt,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'Hybrid localization-first receipt failed';endif
    if(rank==0)write(*,'(a,3(a,i0),5(a,es16.8))')&
      '[HYBRID-WF-LOCALIZATION] symmetry_constraint=off',&
      ' raw_rank=',hybrid_localization_receipt%raw_rank,&
      ' retained_rank=',hybrid_localization_receipt%retained_rank,&
      ' iterations=',hybrid_localization_receipt%iterations,&
      ' spread_min=',hybrid_localization_receipt%spread_min,&
      ' spread_max=',hybrid_localization_receipt%spread_max,&
      ' spread_mean=',hybrid_localization_receipt%spread_mean,&
      ' spread_total=',hybrid_localization_receipt%spread_total,&
      ' unitarity=',hybrid_localization_receipt%transform_unitarity_defect
    deallocate(w90_spreads)
    w90_identity_defect=0d0
    w90_unitarity_defect=hybrid_localization_receipt%transform_unitarity_defect
    w90_closure_defect=global_retained_group_closure_defect
    w90_covariance_defect=global_retained_group_closure_defect
    w90_symmetry_workspace_peak=lcfo_symmetry_workspace_peak;w90_covariance_workspace=0_8
    localization_initial_spread=w90_spread(1);localization_final_spread=w90_spread(1)
    localization_maximum_gradient=0d0;localization_spread_evaluations=0;localization_converged=.true.
    translation_post_gauge_fingerprint=hybrid_localization_receipt%transform_fingerprint

    allocate(localized_center_magnitudes(3,ntarget))
    call compute_dg_periodic_wannier_centers(dc%icomm_tot,ow_core_values,ow_core_weights,&
      core_periodic_phase,localized_centers,localized_center_magnitudes,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'localized Wannier center measurement failed';endif
    call diagnose_dg_point_center_gauge(dc%icomm_tot,ow_core_values,ow_core_weights,&
      global_symmetry_map(:,lcfo_symmetry_worst_operation),&
      global_point_integer_rotations(:,:,lcfo_symmetry_worst_operation),&
      global_point_fractional_translations(:,lcfo_symmetry_worst_operation),localized_centers,&
      retained_closure_search_tolerance,monomial_defect,center_block_leakage,&
      center_representation_unitarity_defect,center_gauge_workspace_peak,center_diagnostic_ok,&
      center_diagnostic_message)
    if(.not.all(ieee_is_finite([monomial_defect,center_block_leakage,&
        center_representation_unitarity_defect])))&
      error stop 'Hybrid WF symmetry diagnostic is nonfinite'
    if(rank==0)write(*,'(a,i0,3(a,es16.8),a,l1)')&
      '[HYBRID-WF-SYMMETRY-DIAGNOSTIC] operation=',lcfo_symmetry_worst_operation,&
      ' individual_wf_monomial=',monomial_defect,' center_block_leakage=',center_block_leakage,&
      ' projected_action_unitarity=',center_representation_unitarity_defect,&
      ' projected_action_closed=',center_diagnostic_ok

    call materialize_ow_distributed_core_to_buffer(dc%icomm_tot,ow_core_values,ow_core_ids,&
      physical_ids,ow_box_values,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'localized core-to-buffer streaming failed';endif
    allocate(ow_direct_core_gradients(3,ntarget,ncore),ow_neighbor_plus_ids(ncore),&
      ow_neighbor_minus_ids(ncore));ow_direct_core_gradients=(0d0,0d0)
    do ix=1,3
      do gradient_distance=1,size(stencil%coef_nab,1)
        do core_index=1,ncore
          raw_ix=int(modulo(ow_core_ids(core_index)-1_8,int(dc%lg_tot%num(1),8)))
          raw_iy=int(modulo((ow_core_ids(core_index)-1_8)/int(dc%lg_tot%num(1),8),&
            int(dc%lg_tot%num(2),8)))
          raw_iz=int((ow_core_ids(core_index)-1_8)/nxy8)
          select case(ix)
          case(1)
            ow_neighbor_minus_ids(core_index)=1_8+int(modulo(raw_ix-gradient_distance,&
              dc%lg_tot%num(1)),8)+int(dc%lg_tot%num(1),8)*(int(raw_iy,8)+int(dc%lg_tot%num(2),8)*int(raw_iz,8))
            ow_neighbor_plus_ids(core_index)=1_8+int(modulo(raw_ix+gradient_distance,&
              dc%lg_tot%num(1)),8)+int(dc%lg_tot%num(1),8)*(int(raw_iy,8)+int(dc%lg_tot%num(2),8)*int(raw_iz,8))
          case(2)
            ow_neighbor_minus_ids(core_index)=1_8+int(raw_ix,8)+int(dc%lg_tot%num(1),8)*(&
              int(modulo(raw_iy-gradient_distance,dc%lg_tot%num(2)),8)+&
              int(dc%lg_tot%num(2),8)*int(raw_iz,8))
            ow_neighbor_plus_ids(core_index)=1_8+int(raw_ix,8)+int(dc%lg_tot%num(1),8)*(&
              int(modulo(raw_iy+gradient_distance,dc%lg_tot%num(2)),8)+&
              int(dc%lg_tot%num(2),8)*int(raw_iz,8))
          case default
            ow_neighbor_minus_ids(core_index)=1_8+int(raw_ix,8)+int(dc%lg_tot%num(1),8)*(&
              int(raw_iy,8)+int(dc%lg_tot%num(2),8)*&
              int(modulo(raw_iz-gradient_distance,dc%lg_tot%num(3)),8))
            ow_neighbor_plus_ids(core_index)=1_8+int(raw_ix,8)+int(dc%lg_tot%num(1),8)*(&
              int(raw_iy,8)+int(dc%lg_tot%num(2),8)*&
              int(modulo(raw_iz+gradient_distance,dc%lg_tot%num(3)),8))
          end select
        enddo
        call materialize_ow_distributed_core_to_buffer(dc%icomm_tot,ow_core_values,ow_core_ids,&
          ow_neighbor_plus_ids,ow_neighbor_plus_values,ok,message)
        if(ok)call materialize_ow_distributed_core_to_buffer(dc%icomm_tot,ow_core_values,ow_core_ids,&
          ow_neighbor_minus_ids,ow_neighbor_minus_values,ok,message)
        if(.not.ok)then;write(0,'(a)')trim(message);error stop 'direct core gradient neighbor stream failed';endif
        ow_direct_core_gradients(ix,:,:)=ow_direct_core_gradients(ix,:,:)+&
          stencil%coef_nab(gradient_distance,ix)*(ow_neighbor_plus_values-ow_neighbor_minus_values)
        deallocate(ow_neighbor_plus_values,ow_neighbor_minus_values)
      enddo
    enddo
    deallocate(ow_neighbor_plus_ids,ow_neighbor_minus_ids,ow_core_values)
    allocate(ow_box_gradients(3,ntarget,nbox));call periodic_box_gradients(ow_box_values,ow_box_size,&
      stencil%coef_nab,ow_box_gradients)
    ! HYBRID_LOCALIZATION_FIRST_ARM_END
    else
    ! HYBRID_CONSTRAINED_LEGACY_ARM_BEGIN
    call prepare_ow_fixed_center_group(ow_core_ids,fixed_center_operations,&
      fixed_center_symmetry_map,fixed_center_product,fixed_center_fractional,&
      fixed_center_inversion_present,fixed_center_group_fingerprint,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'fixed-center point-group construction failed';end if
    fixed_center_group_order=size(fixed_center_operations)
    if(fixed_center_group_order>48)error stop 'fixed_center_group_order>48'
    if(.not.fixed_center_inversion_present)error stop 'fixed-center point group lacks inversion'
    call find_dg_group_identity(fixed_center_product,fixed_center_identity_operation,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'fixed-center group identity construction failed';end if
    call find_dg_group_identity(global_point_product,global_identity_operation,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'global group identity construction failed';end if
    call find_dg_group_identity(global_point_cogroup_product,global_point_cogroup_identity_operation,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'point cogroup identity construction failed';end if
    allocate(translation_product(size(global_translation_subgroup),size(global_translation_subgroup)))
    do io=1,size(global_translation_subgroup);do i=1,size(global_translation_subgroup)
      translation_product(i,io)=findloc(global_translation_subgroup,&
        global_point_product(global_translation_subgroup(i),global_translation_subgroup(io)),dim=1)
    enddo;enddo
    if(any(translation_product<1))error stop 'translation subgroup product is not closed'
    translation_identity_operation=findloc(global_translation_subgroup,global_identity_operation,dim=1)
    if(translation_identity_operation<1)error stop 'translation subgroup lacks the affine identity'
    allocate(translation_canonical_operations(size(global_translation_subgroup)),&
      translation_inverse_operations(size(global_translation_subgroup)),&
      translation_characters(size(global_translation_subgroup),size(global_translation_subgroup)),&
      translation_character_conjugates(size(global_translation_subgroup)))
    call build_dg_finite_abelian_character_table(&
      global_point_fractional_translations(:,global_translation_subgroup),translation_product,&
      translation_identity_operation,dg_ow_symmetry_tolerance,translation_canonical_operations,&
      translation_inverse_operations,translation_character_generator_count,&
      translation_character_generators,translation_element_words,translation_characters,&
      translation_character_conjugates,translation_character_fingerprint,ok,message)
    if(.not.ok)then
      write(0,'(a)')trim(message);error stop 'translation character catalog construction failed'
    endif
    allocate(translation_canonical_product(size(translation_product,1),size(translation_product,2)))
    do io=1,size(translation_product,1);do i=1,size(translation_product,2)
      translation_canonical_product(i,io)=findloc(translation_canonical_operations,&
        translation_product(translation_canonical_operations(i),translation_canonical_operations(io)),dim=1)
    enddo;enddo
    if(any(translation_canonical_product<1))error stop 'canonical translation product is not closed'
    allocate(translation_generator_orders(translation_character_generator_count),&
      translation_generator_characters(size(global_translation_subgroup),&
      translation_character_generator_count))
    do i=1,translation_character_generator_count
      translation_generator_orders(i)=1
      do while(translation_generator_orders(i)<=size(global_translation_subgroup))
        if(maxval(abs(translation_characters(:,translation_character_generators(i))**&
            translation_generator_orders(i)-1d0))<=100d0*dg_ow_symmetry_tolerance)exit
        translation_generator_orders(i)=translation_generator_orders(i)+1
      enddo
      if(translation_generator_orders(i)>size(global_translation_subgroup))&
        error stop 'translation character generator order is not finite'
      translation_generator_characters(:,i)=&
        translation_characters(:,translation_character_generators(i))
    enddo
    call select_dg_group_generators(global_point_product,global_identity_operation,&
      global_affine_generators,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'global affine generator selection failed';end if
    if(rank==0)write(*,'(a,2(a,i0))')'[OW-GS-DIAGNOSTIC] affine_generator_proof',&
      ' group_order=',size(global_point_product,1),' generator_count=',size(global_affine_generators)
    if(rank==0)then
      do i=1,size(global_affine_generators)
        ow_stencil_axis_defect=0d0
        do ix=1,3
          do gradient_distance=1,size(stencil%coef_nab,1)
            ow_stencil_axis_defect=max(ow_stencil_axis_defect,abs(stencil%coef_nab(gradient_distance,ix)-&
              sum(abs(global_point_rotations(:,ix,global_affine_generators(i)))*&
              stencil%coef_nab(gradient_distance,:))))
          enddo
        enddo
        write(*,'(2(a,i0),a,9(i3,1x),a,9(f7.3,1x),a,es16.8)')&
          '[OW-GS-DIAGNOSTIC] affine_generator index=',i,' operation=',global_affine_generators(i),&
          ' integer_R=',global_point_integer_rotations(:,:,global_affine_generators(i)),&
          ' cartesian_R=',global_point_rotations(:,:,global_affine_generators(i)),&
          ' stencil_axis_defect=',ow_stencil_axis_defect
      enddo
    endif
    call measure_dg_grid_map_stencil_defect(dc%icomm_tot,ow_core_ids,&
      global_symmetry_map(:,global_affine_generators),dc%lg_tot%num,system%hgs,&
      global_point_rotations(:,:,global_affine_generators),ow_grid_stencil_defect,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'affine grid-stencil diagnostic failed';endif
    if(rank==0)write(*,'(a,es16.8,a,i0)')'[OW-GS-DIAGNOSTIC] affine grid-stencil defect max=',&
      maxval(ow_grid_stencil_defect),' operation=',maxloc(ow_grid_stencil_defect,dim=1)
    deallocate(ow_grid_stencil_defect)
    allocate(ow_map_probe_values(1,ncore),ow_map_probe_gradients(3,1,ncore))
    do p=1,ncore
      raw_ix=int(modulo(ow_core_ids(p)-1_8,int(dc%lg_tot%num(1),8)))
      raw_iy=int(modulo((ow_core_ids(p)-1_8)/int(dc%lg_tot%num(1),8),&
        int(dc%lg_tot%num(2),8)))
      raw_iz=int((ow_core_ids(p)-1_8)/nxy8)
      ow_map_probe_angle=2d0*acos(-1d0)*(real(raw_ix,8)/real(dc%lg_tot%num(1),8)+&
        2d0*real(raw_iy,8)/real(dc%lg_tot%num(2),8)+&
        3d0*real(raw_iz,8)/real(dc%lg_tot%num(3),8))
      ow_map_probe_values(1,p)=exp(cmplx(0d0,ow_map_probe_angle,8))
    enddo
    do ix=1,3
      ow_map_probe_symbol=0d0
      do gradient_distance=1,size(stencil%coef_nab,1)
        ow_map_probe_symbol=ow_map_probe_symbol+2d0*stencil%coef_nab(gradient_distance,ix)*&
          sin(2d0*acos(-1d0)*real(ix*gradient_distance,8)/real(dc%lg_tot%num(ix),8))
      enddo
      ow_map_probe_gradients(ix,1,:)=cmplx(0d0,ow_map_probe_symbol,8)*ow_map_probe_values(1,:)
    enddo
    call measure_ow_discrete_gradient_map_commutator(dc%icomm_tot,ow_map_probe_values,&
      ow_map_probe_gradients,ow_core_weights,ow_core_ids,&
      global_symmetry_map(:,global_affine_generators),dc%lg_tot%num,stencil%coef_nab,&
      global_point_rotations(:,:,global_affine_generators),ow_gradient_map_commutator,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'affine plane-wave commutator diagnostic failed';endif
    if(rank==0)write(*,'(a,*(es16.8,1x))')&
      '[OW-GS-DIAGNOSTIC] affine plane-wave commutator=',ow_gradient_map_commutator
    deallocate(ow_map_probe_values,ow_map_probe_gradients,ow_gradient_map_commutator)
    allocate(lcfo_total_symmetry_residual(size(global_affine_generators)),&
      lcfo_boundary_symmetry_residual(size(global_affine_generators)),&
      lcfo_interior_symmetry_residual(size(global_affine_generators)))
#ifdef USE_EIGENEXA
    call orthonormalize_dg_distributed_seed_space(dc%icomm_tot,lcfo_occupied_core(1:nstate,:),&
      ow_core_weights,dg_dc_metric_rank_tolerance,orthonormal_lcfo_occupied,&
      orthonormal_lcfo_rank,ok,message)
    if(.not.ok.or.orthonormal_lcfo_rank/=nstate)then
      write(0,'(a)')trim(message);error stop 'LCFO occupied subspace lost rank before symmetry adaptation'
    endif
    ow_saved_eigenexa_comm=info%icomm_o
    call finalize_eigenexa(info)
    info%icomm_o=dc%icomm_tot
    call init_eigenexa_mod(info,nstate,direct_block_only=.true.)
    allocate(occupied_pre_total_residual(fixed_center_group_order),&
      occupied_pre_boundary_residual(fixed_center_group_order),&
      occupied_pre_interior_residual(fixed_center_group_order))
    call measure_dg_rank_fixed_symmetry_residuals_eigenexa(info,dc%icomm_tot,orthonormal_lcfo_occupied,&
      ow_core_weights,fixed_center_symmetry_map,lcfo_boundary_mask,&
      total_residual=occupied_pre_total_residual,boundary_residual=occupied_pre_boundary_residual,&
      interior_residual=occupied_pre_interior_residual,ok=ok,message=message,&
      workspace_peak_bytes=occupied_pre_closure_workspace_peak)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'pre-adaptation occupied closure measurement failed';endif
    adapted_occupied_closure_before=maxval(occupied_pre_total_residual)
    deallocate(occupied_pre_total_residual,occupied_pre_boundary_residual,occupied_pre_interior_residual)
    call finalize_eigenexa(info)
    if(nstate>huge(nstate)/size(global_translation_subgroup))&
      error stop 'translation-subgroup occupied orbit rank overflow'
    call init_eigenexa_mod(info,nstate*size(global_translation_subgroup),direct_block_only=.true.)
    call build_dg_group_averaged_occupied_candidates_eigenexa(info,dc%icomm_tot,&
      orthonormal_lcfo_occupied,ow_core_weights,&
      global_symmetry_map(:,global_translation_subgroup),translation_product,&
      translation_identity_operation,nstate,dg_ow_symmetry_tolerance,&
      translation_adapted_occupied,translation_adapted_spectrum,translation_adapted_rank,&
      translation_adapted_trace,translation_adapted_closure,translation_adapted_gamma_defect,&
      translation_adapted_workspace_peak,ok,message)
    call finalize_eigenexa(info)
    if(.not.ok.or.translation_adapted_rank/=nstate)then
      write(0,'(a)')trim(message);error stop 'translation-subgroup occupied adaptation failed'
    endif
    if(rank==0)write(*,'(a,2(a,i0),4(a,es16.8),a,i0)')&
      '[OW-GS-DIAGNOSTIC] translation_adapted_occupied',&
      ' input_rank=',nstate,' selected_rank=',translation_adapted_rank,&
      ' trace=',translation_adapted_trace,' closure=',translation_adapted_closure,&
      ' gamma_defect=',translation_adapted_gamma_defect,&
      ' electron_count_drift=',abs(translation_adapted_trace-real(nstate,8)),&
      ' workspace_peak_bytes=',translation_adapted_workspace_peak
    if(nstate>huge(nstate)/nstate)error stop 'occupied Hamiltonian MPI count overflow'
    allocate(translation_hamiltonian_overlap_local(nstate,nstate),&
      translation_hamiltonian_overlap(nstate,nstate),translation_occupied_hamiltonian(nstate,nstate),&
      stat=allocation_status)
    call MPI_Allreduce(allocation_status,translation_allocation_status,1,MPI_INTEGER,MPI_MAX,&
      dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.translation_allocation_status/=0)&
      error stop 'occupied Hamiltonian workspace allocation failed collectively'
    translation_hamiltonian_overlap_local=(0d0,0d0)
    do io=1,nstate;do i=1,nstate
      translation_hamiltonian_overlap_local(i,io)=sum(ow_core_weights*&
        conjg(lcfo_occupied_core(i,:))*translation_adapted_occupied(io,:))
    enddo;enddo
    call MPI_Allreduce(translation_hamiltonian_overlap_local,translation_hamiltonian_overlap,&
      nstate*nstate,MPI_DOUBLE_COMPLEX,MPI_SUM,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)error stop 'occupied Hamiltonian overlap reduction failed'
    translation_occupied_hamiltonian=(0d0,0d0)
    do io=1,nstate;do i=1,nstate
      translation_occupied_hamiltonian(i,io)=sum(conjg(translation_hamiltonian_overlap(:,i))*&
        lcfo_retained_eigenvalues(1:nstate)*translation_hamiltonian_overlap(:,io))
    enddo;enddo
    if(.not.all(ieee_is_finite(real(translation_occupied_hamiltonian))).or.&
        .not.all(ieee_is_finite(aimag(translation_occupied_hamiltonian))))&
      error stop 'occupied Hamiltonian projection is nonfinite'
    if(maxval(abs(translation_occupied_hamiltonian-&
        conjg(transpose(translation_occupied_hamiltonian))))>&
        dg_ow_symmetry_tolerance*max(1d0,maxval(abs(translation_occupied_hamiltonian))))&
      error stop 'occupied Hamiltonian projection is not Hermitian'
    if(nstate>huge(nstate)/size(global_point_representatives))&
      error stop 'point-cogroup occupied orbit rank overflow'
    call init_eigenexa_mod(info,nstate*size(global_point_representatives),direct_block_only=.true.)
    call build_dg_cocycle_averaged_occupied_candidates_eigenexa(info,dc%icomm_tot,&
      translation_adapted_occupied,ow_core_weights,&
      global_symmetry_map(:,global_translation_subgroup),&
      global_symmetry_map(:,global_point_representatives),global_point_cogroup_product,&
      global_translation_cocycle,global_point_cogroup_identity_operation,nstate,dg_ow_symmetry_tolerance,&
      adapted_occupied_candidates,adapted_occupied_spectrum,adapted_occupied_rank,&
      adapted_occupied_trace,adapted_occupied_closure,adapted_occupied_gamma_defect,&
      adapted_occupied_workspace_peak,ok,message,adapted_occupied_selected_edge,&
      adapted_occupied_rejected_edge,adapted_occupied_cluster_gap,&
      occupied_hamiltonian=translation_occupied_hamiltonian,&
      secondary_selected_edge=adapted_occupied_secondary_selected_edge,&
      secondary_rejected_edge=adapted_occupied_secondary_rejected_edge,&
      secondary_cluster_gap=adapted_occupied_secondary_cluster_gap,&
      primary_boundary_dimension=adapted_occupied_selected_block_dimension,&
      hamiltonian_fingerprint=adapted_occupied_hamiltonian_fingerprint,&
      secondary_eigensystem_residual=adapted_occupied_secondary_residual)
    call finalize_eigenexa(info)
    if(.not.ok.or.adapted_occupied_rank/=nstate)then
      if(rank==0)write(0,'(a,3(a,es24.16))')&
        '[OW-GS-DIAGNOSTIC] point_cogroup_adaptation_rejected',&
        ' selected_edge=',adapted_occupied_selected_edge,&
        ' rejected_edge=',adapted_occupied_rejected_edge,&
        ' cluster_gap=',adapted_occupied_cluster_gap
      write(0,'(a)')trim(message);error stop 'point-cogroup occupied-subspace adaptation failed'
    endif
    if(rank==0)write(*,'(a,a,i0,4(a,es24.16),a,i0)')&
      '[OW-GS-DIAGNOSTIC] point_cogroup_hamiltonian_tiebreak',&
      ' primary_boundary_dimension=',adapted_occupied_selected_block_dimension,&
      ' selected_edge=',adapted_occupied_secondary_selected_edge,&
      ' rejected_edge=',adapted_occupied_secondary_rejected_edge,&
      ' cluster_gap=',adapted_occupied_secondary_cluster_gap,&
      ' residual=',adapted_occupied_secondary_residual,&
      ' fingerprint=',adapted_occupied_hamiltonian_fingerprint
    occupied_composition_fingerprint=ieor(ishftc(occupied_composition_fingerprint,11),&
      adapted_occupied_hamiltonian_fingerprint)
    deallocate(translation_hamiltonian_overlap_local,translation_hamiltonian_overlap,&
      translation_occupied_hamiltonian)
    adapted_occupied_workspace_peak=max(adapted_occupied_workspace_peak,translation_adapted_workspace_peak,&
      occupied_pre_closure_workspace_peak,occupied_composition_peak,projector_composition_peak)
    allocate(occupied_overlap_local(nstate,nstate),occupied_overlap_global(nstate,nstate))
    do io=1,nstate;do i=1,nstate
      occupied_overlap_local(i,io)=sum(ow_core_weights*conjg(orthonormal_lcfo_occupied(i,:))*&
        adapted_occupied_candidates(io,:))
    enddo;enddo
    call MPI_Allreduce(occupied_overlap_local,occupied_overlap_global,nstate*nstate,&
      MPI_DOUBLE_COMPLEX,MPI_SUM,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)error stop 'occupied subspace-distance reduction failed'
    adapted_occupied_subspace_distance=sqrt(max(0d0,1d0-sum(abs(occupied_overlap_global)**2)/real(nstate,8)))
    deallocate(occupied_overlap_local,occupied_overlap_global)
    allocate(occupied_density_before(ncore),occupied_density_after(ncore),occupied_density_difference(ncore))
    occupied_density_before=sum(abs(lcfo_occupied_core(1:nstate,:))**2,dim=1)
    occupied_density_after=sum(abs(adapted_occupied_candidates)**2,dim=1)
    allocate(spectral_occupied_density(ncore),source=occupied_density_after)
    occupied_density_difference=abs(occupied_density_after-occupied_density_before)
    local_occupied_density_interior_difference=sum(ow_core_weights*occupied_density_difference**2,&
      mask=.not.lcfo_boundary_mask)
    local_occupied_density_boundary_difference=sum(ow_core_weights*occupied_density_difference**2,&
      mask=lcfo_boundary_mask)
    local_occupied_density_interior_norm=sum(ow_core_weights*occupied_density_before**2,&
      mask=.not.lcfo_boundary_mask)
    local_occupied_density_boundary_norm=sum(ow_core_weights*occupied_density_before**2,&
      mask=lcfo_boundary_mask)
    call MPI_Allreduce(local_occupied_density_interior_difference,&
      adapted_occupied_density_interior_difference,1,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
    call MPI_Allreduce(local_occupied_density_boundary_difference,&
      adapted_occupied_density_boundary_difference,1,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
    call MPI_Allreduce(local_occupied_density_interior_norm,global_occupied_density_interior_norm,&
      1,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
    call MPI_Allreduce(local_occupied_density_boundary_norm,global_occupied_density_boundary_norm,&
      1,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)error stop 'occupied density-difference reduction failed'
    adapted_occupied_density_interior_difference=sqrt(adapted_occupied_density_interior_difference/&
      max(tiny(1d0),global_occupied_density_interior_norm))
    adapted_occupied_density_boundary_difference=sqrt(adapted_occupied_density_boundary_difference/&
      max(tiny(1d0),global_occupied_density_boundary_norm))
    deallocate(occupied_density_before,occupied_density_after,occupied_density_difference)
    deallocate(lcfo_occupied_core)
    adapted_occupied_selected_block_dimension=count(&
      abs(adapted_occupied_spectrum-adapted_occupied_selected_edge)<=&
      dg_ow_symmetry_tolerance*max(1d0,abs(adapted_occupied_selected_edge)))
    global_seed_values(1:nstate,:)=adapted_occupied_candidates
    deallocate(adapted_occupied_candidates,adapted_occupied_spectrum,translation_adapted_occupied,&
      translation_adapted_spectrum,orthonormal_lcfo_occupied)
    if(rank==0)write(*,'(a,2(a,i0),7(a,es16.8),a,i0)')&
      '[OW-GS-DIAGNOSTIC] point_cogroup_adapted_occupied',&
      ' input_rank=',nstate,' selected_rank=',adapted_occupied_rank,&
      ' trace=',adapted_occupied_trace,' closure=',adapted_occupied_closure,&
      ' gamma_defect=',adapted_occupied_gamma_defect,&
      ' electron_count_drift=',abs(adapted_occupied_trace-real(nstate,8)),&
      ' selected_edge=',adapted_occupied_selected_edge,&
      ' rejected_edge=',adapted_occupied_rejected_edge,&
      ' cluster_gap=',adapted_occupied_cluster_gap,&
      ' workspace_peak_bytes=',adapted_occupied_workspace_peak
    allocate(occupied_affine_total_residual(size(global_affine_generators)),&
      occupied_affine_boundary_residual(size(global_affine_generators)),&
      occupied_affine_interior_residual(size(global_affine_generators)),&
      projection_affine_total_residual(size(global_affine_generators)),&
      projection_affine_boundary_residual(size(global_affine_generators)),&
      projection_affine_interior_residual(size(global_affine_generators)))
    call init_eigenexa_mod(info,nstate,direct_block_only=.true.)
    call measure_dg_rank_fixed_symmetry_residuals_eigenexa(info,dc%icomm_tot,&
      global_seed_values(1:nstate,:),ow_core_weights,&
      global_symmetry_map(:,global_affine_generators),lcfo_boundary_mask,&
      occupied_affine_total_residual,occupied_affine_boundary_residual,&
      occupied_affine_interior_residual,ok,message,occupied_affine_workspace_peak)
    call finalize_eigenexa(info)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'occupied affine diagnostic failed';endif
    call init_eigenexa_mod(info,global_projection_count,direct_block_only=.true.)
    call measure_dg_rank_fixed_symmetry_residuals_eigenexa(info,dc%icomm_tot,&
      global_seed_values(nstate+1:ntarget,:),ow_core_weights,&
      global_symmetry_map(:,global_affine_generators),lcfo_boundary_mask,&
      projection_affine_total_residual,projection_affine_boundary_residual,&
      projection_affine_interior_residual,ok,message,projection_affine_workspace_peak)
    call finalize_eigenexa(info)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'projection affine diagnostic failed';endif
    if(rank==0)write(*,'(a,6(a,es16.8))')'[OW-GS-DIAGNOSTIC] affine_block_closure',&
      ' occupied_total=',maxval(occupied_affine_total_residual),&
      ' occupied_boundary=',maxval(occupied_affine_boundary_residual),&
      ' occupied_interior=',maxval(occupied_affine_interior_residual),&
      ' projection_total=',maxval(projection_affine_total_residual),&
      ' projection_boundary=',maxval(projection_affine_boundary_residual),&
      ' projection_interior=',maxval(projection_affine_interior_residual)
    deallocate(occupied_affine_total_residual,occupied_affine_boundary_residual,&
      occupied_affine_interior_residual,projection_affine_total_residual,&
      projection_affine_boundary_residual,projection_affine_interior_residual)
    call orthonormalize_dg_distributed_seed_space(dc%icomm_tot,global_seed_values,ow_core_weights,&
      dg_dc_metric_rank_tolerance,global_closed_core,global_retained_rank,ok,message)
    global_required_retained_rank=global_retained_rank
    if(.not.ok.or.global_retained_rank/=ntarget)then
      write(0,'(a)')trim(message);error stop 'adapted occupied plus complete-s+p seed lost target rank'
    endif
    call init_eigenexa_mod(info,size(global_seed_values,1),direct_block_only=.true.)
    call measure_dg_rank_fixed_symmetry_residuals_eigenexa(info,dc%icomm_tot,global_closed_core,ow_core_weights,&
      global_symmetry_map(:,global_affine_generators),lcfo_boundary_mask,&
      total_residual=lcfo_total_symmetry_residual,&
      boundary_residual=lcfo_boundary_symmetry_residual,&
      interior_residual=lcfo_interior_symmetry_residual,ok=ok,message=message,&
      workspace_peak_bytes=lcfo_symmetry_workspace_peak)
    call finalize_eigenexa(info)
    info%icomm_o=ow_saved_eigenexa_comm
    call init_eigenexa_mod(info,system%no)
#else
    ok=.false.;message='production overlapping-Wannier symmetry requires EigenExa'
#endif
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'LCFO occupied symmetry measurement failed';end if
    if(rank==0)write(*,'(a,i0)')'[OW-GS-DIAGNOSTIC] lcfo_symmetry_workspace_peak_bytes=',&
      lcfo_symmetry_workspace_peak
    lcfo_symmetry_worst_generator_index=maxloc(lcfo_total_symmetry_residual,dim=1)
    lcfo_symmetry_worst_operation=global_affine_generators(lcfo_symmetry_worst_generator_index)
    if(rank==0)write(*,'(a,i0,3(a,es16.8))')&
      '[OW-GS-DIAGNOSTIC] lcfo_symmetry_worst_operation=',lcfo_symmetry_worst_operation,&
      ' total_residual=',lcfo_total_symmetry_residual(lcfo_symmetry_worst_generator_index),&
      ' boundary_residual=',lcfo_boundary_symmetry_residual(lcfo_symmetry_worst_generator_index),&
      ' interior_residual=',lcfo_interior_symmetry_residual(lcfo_symmetry_worst_generator_index)
    if(maxval(lcfo_total_symmetry_residual)>dg_ow_symmetry_tolerance)&
      error stop 'LCFO seed space is not closed under the full affine action'
    if(rank==0)write(*,'(a,3(a,i0))')'[OW-GS-DIAGNOSTIC] LCFO_occupied_closure',&
      ' input_rank=',ntarget,' closure_rank=',global_required_retained_rank,&
      ' candidate_rank=',global_retained_rank
    if(.not.ok.or.global_retained_rank<ntarget.or.&
        global_required_retained_rank/=ntarget)then
      write(0,'(a)')trim(message);error stop 'global symmetry-closed Wannier construction failed'
    end if
    allocate(spectral_empty_moments(ncore,1),spectral_shared_density(ncore,0),&
      spectral_basin_generator_maps(ncore,size(global_affine_generators)))
    spectral_empty_moments(:,1)=max(0d0,sum(abs(global_closed_core)**2,dim=1)-spectral_occupied_density)
    do i=1,size(global_affine_generators)
      spectral_basin_generator_maps(:,i)=int(global_symmetry_map(:,global_affine_generators(i)))
    enddo
    call build_dg_periodic_spectral_basins(dc%icomm_tot,ow_core_ids,dc%lg_tot%num,&
      spectral_occupied_density,spectral_empty_moments,spectral_shared_density,&
      spectral_basin_generator_maps,dg_ow_symmetry_tolerance,spectral_basin_labels,&
      spectral_basin_count,spectral_single_basin_map,spectral_basin_fingerprint,&
      spectral_workspace_peak,ok,message)
    if(.not.ok)then
      if(rank==0)write(0,'(a)')trim(message)
      error stop 'pre-Wannier occupied/empty spectral basin construction failed'
    endif
    if(rank==0)write(*,'(a,2(a,i0),a,i0)')'[OW-GS-DIAGNOSTIC] pre_wannier_spectral_basins',&
      ' basin_count=',spectral_basin_count,' generator_count=',size(global_affine_generators),&
      ' workspace_peak_bytes=',spectral_workspace_peak
    deallocate(spectral_basin_generator_maps,spectral_occupied_density,spectral_empty_moments,&
      spectral_shared_density)
    allocate(w90_anchors(ntarget,ncore),source=global_seed_values,stat=allocation_status)
    call MPI_Allreduce(MPI_IN_PLACE,allocation_status,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_status/=0)&
      error stop 'established Wannier90 seed-anchor retention failed collectively'
    w90_byte_limit=8_8*1024_8*1024_8*1024_8
    call assemble_dg_w90_gamma_a_matrix(dc%icomm_tot,global_closed_core,w90_anchors,&
      ow_core_weights,dg_ow_symmetry_tolerance,w90_byte_limit,w90_seed_a_matrix,&
      w90_seed_a_workspace,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'pre-DMN Wannier90 A assembly failed';endif
    if(allocated(global_seed_values))deallocate(global_seed_values)
    ! Diagonalize the spatial-basin projectors in the retained frame.  Only
    ! spectra are retained for every basin; representative eigenvectors are
    ! regenerated after the symmetry-consistent ranks have been selected.
    call fingerprint_ow_spatial_frame(dc%icomm_tot,ow_core_ids,global_closed_core,ow_core_weights,&
      dg_ow_symmetry_tolerance,spectral_frame_fingerprint,spectral_frame_defect,ok)
    if(.not.ok)error stop 'spectral retained-frame fingerprint failed'
    call prepare_dg_spectral_basin_operators(dc%icomm_tot,ow_core_ids,int(expected_core_count),&
      global_closed_core,ow_core_weights,spectral_basin_labels,spectral_basin_count,&
      spectral_frame_fingerprint,spectral_frame_defect,spectral_basin_fingerprint,&
      dg_ow_symmetry_tolerance,spectral_prepared_basins,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'spectral basin operator preparation failed';endif
    allocate(spectral_basin_spectra(ntarget,spectral_basin_count),&
      spectral_block_ends(ntarget,spectral_basin_count),spectral_selected_ranks(spectral_basin_count))
    spectral_block_ends=.false.
    do spectral_basin=1,spectral_basin_count
      call project_dg_prepared_spectral_basin_operator(dc%icomm_tot,spectral_prepared_basins,&
        global_closed_core,ow_core_weights,spectral_basin,spectral_basin_operator,&
        spectral_operator_hermiticity,spectral_operator_trace,spectral_operator_fingerprint,&
        spectral_operation_workspace,ok,message)
      if(.not.ok)then;write(0,'(a)')trim(message);error stop 'spectral basin operator projection failed';endif
      spectral_workspace_peak=max(spectral_workspace_peak,spectral_operation_workspace)
      call diagonalize_dg_spectral_basin_operator(dc%icomm_tot,spectral_basin_operator,&
        spectral_operator_fingerprint,dg_ow_symmetry_tolerance,spectrum,spectral_block_offsets,&
        spectral_eigensystem_residual,spectral_eigensystem_fingerprint,spectral_operation_workspace,ok,message)
      if(.not.ok)then;write(0,'(a)')trim(message);error stop 'spectral basin eigensystem failed';endif
      spectral_workspace_peak=max(spectral_workspace_peak,spectral_operation_workspace)
      spectral_basin_spectra(:,spectral_basin)=spectrum
      do i=1,size(spectral_block_offsets)-1
        spectral_block_ends(spectral_block_offsets(i+1)-1,spectral_basin)=.true.
      enddo
      deallocate(spectrum,spectral_block_offsets,spectral_basin_operator)
    enddo
    call select_dg_spectral_basin_channel_ranks(dc%icomm_tot,spectral_basin_spectra,&
      spectral_block_ends,spectral_single_basin_map,ntarget,dg_ow_symmetry_tolerance,&
      spectral_selected_ranks,spectral_catalog_fingerprint,spectral_operation_workspace,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'spectral basin channel-rank selection failed';endif
    spectral_workspace_peak=max(spectral_workspace_peak,spectral_operation_workspace)

    ! Find one deterministic representative per basin orbit and retain only
    ! the selected eigenvectors for those representatives.
    allocate(spectral_orbit_id(spectral_basin_count),spectral_orbit_representatives(spectral_basin_count))
    spectral_orbit_id=0;spectral_orbit_count=0
    do spectral_basin=1,spectral_basin_count
      if(spectral_orbit_id(spectral_basin)/=0)cycle
      spectral_orbit_count=spectral_orbit_count+1
      spectral_orbit_representatives(spectral_orbit_count)=spectral_basin
      spectral_orbit_id(spectral_basin)=spectral_orbit_count
      do
        spectral_representative_count=count(spectral_orbit_id==spectral_orbit_count)
        do spectral_target_basin=1,spectral_basin_count
          if(spectral_orbit_id(spectral_target_basin)/=spectral_orbit_count)cycle
          do i=1,size(spectral_single_basin_map,2)
            spectral_orbit=spectral_single_basin_map(spectral_target_basin,i)
            if(spectral_orbit_id(spectral_orbit)==0)spectral_orbit_id(spectral_orbit)=spectral_orbit_count
          enddo
        enddo
        if(count(spectral_orbit_id==spectral_orbit_count)==spectral_representative_count)exit
      enddo
    enddo
    spectral_representative_count=0
    do spectral_orbit=1,spectral_orbit_count
      spectral_representative_count=spectral_representative_count+&
        spectral_selected_ranks(spectral_orbit_representatives(spectral_orbit))
    enddo
    allocate(spectral_representative_vectors(ntarget,spectral_representative_count))
    spectral_representative_column=1
    do spectral_orbit=1,spectral_orbit_count
      spectral_basin=spectral_orbit_representatives(spectral_orbit)
      call project_dg_prepared_spectral_basin_operator(dc%icomm_tot,spectral_prepared_basins,&
        global_closed_core,ow_core_weights,spectral_basin,spectral_basin_operator,&
        spectral_operator_hermiticity,spectral_operator_trace,spectral_operator_fingerprint,&
        spectral_operation_workspace,ok,message)
      if(.not.ok)then;write(0,'(a)')trim(message);error stop 'representative basin projection failed';endif
      call diagonalize_dg_spectral_basin_operator(dc%icomm_tot,spectral_basin_operator,&
        spectral_operator_fingerprint,dg_ow_symmetry_tolerance,spectrum,spectral_block_offsets,&
        spectral_eigensystem_residual,spectral_eigensystem_fingerprint,spectral_operation_workspace,ok,message)
      if(.not.ok)then;write(0,'(a)')trim(message);error stop 'representative basin eigensystem failed';endif
      i=spectral_selected_ranks(spectral_basin)
      if(i>0)then
        spectral_representative_vectors(:,spectral_representative_column:spectral_representative_column+i-1)=&
          spectral_basin_operator(:,1:i)
        spectral_representative_column=spectral_representative_column+i
      endif
      deallocate(spectrum,spectral_block_offsets,spectral_basin_operator)
    enddo
    call release_dg_prepared_spectral_basins(spectral_prepared_basins)

    call assemble_dg_distributed_basis_symmetry_overlap_rows(dc%icomm_tot,global_closed_core,&
      ow_core_weights,global_symmetry_map(:,global_affine_generators),spectral_row_ids,fixed_center_rows,&
      spectral_operation_workspace,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'spectral generator action assembly failed';endif
    call propagate_dg_spectral_basin_orbit_channels(dc%icomm_tot,spectral_row_ids,fixed_center_rows,&
      spectral_single_basin_map,spectral_selected_ranks,spectral_representative_vectors,&
      fixed_center_group_fingerprint,spectral_catalog_fingerprint,dg_ow_symmetry_tolerance,&
      spectral_trial_rows,spectral_channel_gram_defect,spectral_channel_fingerprint,&
      spectral_operation_workspace,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'spectral basin channel propagation failed';endif
    call build_dg_spectral_channel_generator_actions(dc%icomm_tot,spectral_row_ids,fixed_center_rows,&
      spectral_trial_rows,spectral_single_basin_map,spectral_selected_ranks,fixed_center_group_fingerprint,&
      spectral_channel_fingerprint,dg_ow_symmetry_tolerance,spectral_wannier_action_rows,&
      spectral_action_unitarity,spectral_action_block_defect,spectral_action_fingerprint,&
      spectral_operation_workspace,ok,message)
    deallocate(fixed_center_rows,spectral_representative_vectors,spectral_orbit_id,&
      spectral_orbit_representatives,spectral_basin_spectra,spectral_block_ends)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'spectral Wannier generator action failed';endif
    spectral_workspace_peak=max(spectral_workspace_peak,spectral_operation_workspace)
    if(rank==0)write(*,'(a,3(a,i0),3(a,es16.8),a,i0)')'[OW-GS-DIAGNOSTIC] spectral_wannier_frame',&
      ' basin_count=',spectral_basin_count,' orbit_count=',spectral_orbit_count,&
      ' representative_columns=',spectral_representative_count,' frame_defect=',spectral_frame_defect,&
      ' gram_defect=',spectral_channel_gram_defect,' action_defect=',spectral_action_block_defect,&
      ' workspace_peak_bytes=',spectral_workspace_peak
    deallocate(spectral_wannier_action_rows)
    deallocate(spectral_selected_ranks,spectral_single_basin_map,spectral_basin_labels)
#ifdef USE_EIGENEXA
    if(translation_character_generator_count==0)then
      call assemble_dg_distributed_basis_symmetry_overlap_rows(dc%icomm_tot,global_closed_core,&
        ow_core_weights,global_symmetry_map(:,global_identity_operation:global_identity_operation),&
        translation_row_ids,fixed_center_rows,fixed_center_operation_workspace,ok,message)
      if(.not.ok)then
        write(0,'(a)')trim(message);error stop 'trivial translation identity row assembly failed'
      endif
      allocate(translation_generator_rows(size(translation_row_ids),ntarget,0),&
        translation_gamma_rows(size(translation_row_ids),ntarget),stat=allocation_status)
      call MPI_Allreduce(allocation_status,translation_allocation_status,1,MPI_INTEGER,MPI_MAX,&
        dc%icomm_tot,ierr)
      if(ierr/=MPI_SUCCESS.or.translation_allocation_status/=0)&
        error stop 'trivial translation row allocation failed collectively'
      deallocate(fixed_center_rows)
    endif
    do i=1,translation_character_generator_count
      io=translation_canonical_operations(translation_character_generators(i))
      io=global_translation_subgroup(io)
      call assemble_dg_distributed_basis_symmetry_overlap_rows(dc%icomm_tot,global_closed_core,&
        ow_core_weights,global_symmetry_map(:,io:io),translation_stream_row_ids,fixed_center_rows,&
        fixed_center_operation_workspace,ok,message)
      if(.not.ok)then
        write(0,'(a)')trim(message);error stop 'translation generator row assembly failed'
      endif
      if(i==1)then
        allocate(translation_row_ids(size(translation_stream_row_ids)),source=translation_stream_row_ids,&
          stat=allocation_status)
        call MPI_Allreduce(allocation_status,translation_allocation_status,1,MPI_INTEGER,MPI_MAX,&
          dc%icomm_tot,ierr)
        if(ierr/=MPI_SUCCESS.or.translation_allocation_status/=0)&
          error stop 'translation row ownership allocation failed collectively'
        allocate(translation_generator_rows(size(translation_row_ids),ntarget,&
          translation_character_generator_count),translation_gamma_rows(size(translation_row_ids),ntarget),&
          stat=allocation_status)
        call MPI_Allreduce(allocation_status,translation_allocation_status,1,MPI_INTEGER,MPI_MAX,&
          dc%icomm_tot,ierr)
        if(ierr/=MPI_SUCCESS.or.translation_allocation_status/=0)&
          error stop 'translation generator allocation failed collectively'
      elseif(any(translation_stream_row_ids/=translation_row_ids))then
        error stop 'translation generator row ownership changed during streaming'
      endif
      translation_generator_rows(:,:,i)=fixed_center_rows(:,:,1)
      deallocate(translation_stream_row_ids,fixed_center_rows)
    enddo
    allocate(translation_gamma_local_row(ntarget),translation_gamma_global_row(ntarget),stat=allocation_status)
    call MPI_Allreduce(allocation_status,translation_allocation_status,1,MPI_INTEGER,MPI_MAX,&
      dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.translation_allocation_status/=0)&
      error stop 'translation Gamma row allocation failed collectively'
    translation_gamma_rows=(0d0,0d0)
    do io=1,ntarget
      do i=1,ntarget
        translation_gamma_local_row(i)=sum(ow_core_weights*conjg(global_closed_core(io,:))*&
          conjg(global_closed_core(i,:)))
      enddo
      call MPI_Allreduce(translation_gamma_local_row,translation_gamma_global_row,ntarget,&
        MPI_DOUBLE_COMPLEX,MPI_SUM,dc%icomm_tot,ierr)
      if(ierr/=MPI_SUCCESS)error stop 'translation Gamma row reduction failed'
      p=findloc(translation_row_ids,int(io,8),dim=1)
      if(p>0)translation_gamma_rows(p,:)=translation_gamma_global_row
    enddo
    deallocate(translation_gamma_local_row,translation_gamma_global_row)
    ow_saved_eigenexa_comm=info%icomm_o
    call finalize_eigenexa(info)
    if(ntarget>huge(ntarget)/2)error stop 'translation-sector EigenExa extent overflow'
    info%icomm_o=dc%icomm_tot
    call init_eigenexa_mod(info,2*ntarget,direct_block_only=.true.)
    call split_dg_translation_character_sector_eigenexa(info,dc%icomm_tot,translation_row_ids,&
      translation_generator_rows,translation_gamma_rows,translation_generator_characters,&
      translation_generator_orders,translation_element_words,translation_character_conjugates,1,&
      dg_ow_symmetry_tolerance,translation_character_fingerprint,translation_sector_rows,translation_sector_rank,&
      translation_identity_defect,translation_unitarity_defect,translation_commutator_defect,&
      translation_order_defect,translation_gamma_pairing_defect,translation_sector_fingerprint,&
      translation_sector_workspace_peak,ok,message)
    call finalize_eigenexa(info)
    info%icomm_o=ow_saved_eigenexa_comm
    call init_eigenexa_mod(info,system%no)
    if(.not.ok)then
      write(0,'(a)')trim(message);error stop 'translation character-sector split failed'
    endif
    if(rank==0)write(*,'(a,2(a,i0),5(a,es16.8),a,i0)')&
      '[OW-GS-DIAGNOSTIC] translation_character_sector',&
      ' character_count=',size(translation_characters,1),' multiplicity=',translation_sector_rank,&
      ' identity_defect=',translation_identity_defect,' unitarity_defect=',translation_unitarity_defect,&
      ' commutator_defect=',translation_commutator_defect,' order_defect=',translation_order_defect,&
      ' gamma_pairing_defect=',translation_gamma_pairing_defect,&
      ' workspace_peak_bytes=',translation_sector_workspace_peak
    deallocate(translation_sector_rows)
#endif
    deallocate(lcfo_core_ids)
    global_retained_rank=ntarget
    allocate(fixed_center_identity(ntarget,ntarget),fixed_center_eigenvalues(ntarget))
    fixed_center_identity=(0d0,0d0);fixed_center_eigenvalues=0d0
    do io=1,ntarget;fixed_center_identity(io,io)=1d0;enddo
    fixed_center_dmn_workspace_peak=0_8
    spectral_action_aggregate_fingerprint=fixed_center_group_fingerprint;writer_ok=.true.
    if(rank==0)call begin_sawf_dmn(fixed_center_dmn_writer,'overlapping_wannier_mlwf.dmn',&
      ntarget,ntarget,fixed_center_group_order,dg_ow_symmetry_tolerance,writer_ok,message)
    call MPI_Bcast(writer_ok,1,MPI_LOGICAL,0,dc%icomm_tot,ierr)
    if(.not.writer_ok)error stop 'fixed-center DMN transaction could not begin'
    do fixed_center_operation=1,fixed_center_group_order
      call assemble_dg_distributed_basis_symmetry_overlap_rows(dc%icomm_tot,global_closed_core,&
        ow_core_weights,fixed_center_symmetry_map(:,fixed_center_operation:fixed_center_operation),&
        fixed_center_row_ids,fixed_center_rows,fixed_center_operation_workspace,ok,message)
      if(.not.ok)then
        if(rank==0)call abort_sawf_dmn(fixed_center_dmn_writer)
        write(0,'(a)')trim(message);error stop 'fixed-center row representation assembly failed'
      endif
      fixed_center_dmn_workspace_peak=max(fixed_center_dmn_workspace_peak,fixed_center_operation_workspace)
      call gather_dg_single_symmetry_representation(dc%icomm_tot,fixed_center_row_ids,&
        fixed_center_rows,1,0,fixed_center_representation,fixed_center_operation_workspace,ok,message)
      if(.not.ok)then
        if(rank==0)call abort_sawf_dmn(fixed_center_dmn_writer)
        write(0,'(a)')trim(message);error stop 'fixed-center representation gather failed'
      endif
      fixed_center_dmn_workspace_peak=max(fixed_center_dmn_workspace_peak,fixed_center_operation_workspace)
      if(rank==0)call convert_sawf_pullback_to_active_representation(&
        fixed_center_representation,ok,message)
      call MPI_Bcast(ok,1,MPI_LOGICAL,0,dc%icomm_tot,ierr)
      if(.not.ok)then
        if(rank==0)write(0,'(a)')trim(message)
        if(rank==0)call abort_sawf_dmn(fixed_center_dmn_writer)
        error stop 'fixed-center pullback representation conversion failed'
      endif
      allocation_status=0
      if(rank==0)then
        allocate(w90_seed_representation(ntarget,ntarget),stat=allocation_status)
        if(allocation_status==0)w90_seed_representation=matmul(conjg(transpose(w90_seed_a_matrix)),&
          matmul(fixed_center_representation,w90_seed_a_matrix))
      endif
      call MPI_Bcast(allocation_status,1,MPI_INTEGER,0,dc%icomm_tot,ierr)
      if(ierr/=MPI_SUCCESS.or.allocation_status/=0)then
        if(rank==0)call abort_sawf_dmn(fixed_center_dmn_writer)
        error stop 'fixed-center seed representation allocation failed'
      endif
      spectral_action_aggregate_fingerprint=ieor(spectral_action_aggregate_fingerprint,&
        ishftc(int(fixed_center_operation,8),modulo(fixed_center_operation,63)))
      writer_ok=.true.
      if(rank==0)call append_sawf_dmn_operation(fixed_center_dmn_writer,fixed_center_operation,&
        w90_seed_representation,fixed_center_representation,fixed_center_eigenvalues,&
        w90_seed_a_matrix,fixed_center_operation==fixed_center_identity_operation,writer_ok,message)
      if(rank==0.and..not.writer_ok)write(0,'(a,i0,2a)')&
        '[OW-GS-DIAGNOSTIC] fixed-center DMN append operation=',fixed_center_operation,&
        ' rejected: ',trim(message)
      call MPI_Bcast(writer_ok,1,MPI_LOGICAL,0,dc%icomm_tot,ierr)
      if(.not.writer_ok)then
        if(rank==0)call abort_sawf_dmn(fixed_center_dmn_writer)
        error stop 'fixed-center DMN operation append failed'
      endif
      deallocate(fixed_center_row_ids,fixed_center_rows,fixed_center_representation)
      if(rank==0)deallocate(w90_seed_representation)
    enddo
    writer_ok=.true.
    if(rank==0)call finish_sawf_dmn(fixed_center_dmn_writer,&
      fixed_center_operations,writer_ok,message)
    call MPI_Bcast(writer_ok,1,MPI_LOGICAL,0,dc%icomm_tot,ierr)
    if(.not.writer_ok)then
      if(rank==0)write(0,'(a)')trim(message)
      if(rank==0)call abort_sawf_dmn(fixed_center_dmn_writer)
      error stop 'fixed-center DMN transaction could not finish'
    endif
    deallocate(fixed_center_identity,fixed_center_eigenvalues)
    deallocate(spectral_trial_rows,spectral_row_ids)
    if(allocated(spectral_wannier_action_rows))deallocate(spectral_wannier_action_rows)
    allocate(initial_core_ids(ncore))
    core_index=0
    do p=1,nbox
      if(.not.core_mask(p))cycle
      core_index=core_index+1
      initial_core_ids(core_index)=physical_ids(p)
      ow_core_weights(core_index)=weights(p)
      ow_core_box_positions(core_index)=p
      core_periodic_phase(1,core_index)=exp(cmplx(0d0,2d0*pi*real(modulo(physical_ids(p)-1_8,&
        int(dc%lg_tot%num(1),8)),8)/real(dc%lg_tot%num(1),8),8))
      core_periodic_phase(2,core_index)=exp(cmplx(0d0,2d0*pi*real(modulo((physical_ids(p)-1_8)/&
        int(dc%lg_tot%num(1),8),int(dc%lg_tot%num(2),8)),8)/real(dc%lg_tot%num(2),8),8))
      core_periodic_phase(3,core_index)=exp(cmplx(0d0,2d0*pi*real((physical_ids(p)-1_8)/&
        nxy8,8)/real(dc%lg_tot%num(3),8),8))
    enddo
    if(core_index/=ncore)error stop 'initial retained core extent is incomplete'
    call materialize_ow_distributed_core_to_buffer(dc%icomm_tot,global_closed_core,ow_core_ids,&
      initial_core_ids,ow_core_values,ok,message)
    if(.not.ok)then
      write(0,'(a)')trim(message)
      error stop 'initial retained core redistribution failed'
    endif
    call reindex_dg_point_maps_between_row_layouts(dc%icomm_tot,ow_core_ids,initial_core_ids,&
      global_symmetry_map,reindexed_global_symmetry_map,ok,message)
    if(.not.ok)then
      write(0,'(a)')trim(message)
      error stop 'global symmetry-map row-layout transition failed'
    endif
    call reindex_dg_point_maps_between_row_layouts(dc%icomm_tot,ow_core_ids,initial_core_ids,&
      fixed_center_symmetry_map,reindexed_fixed_center_symmetry_map,ok,message)
    if(.not.ok)then
      write(0,'(a)')trim(message)
      error stop 'fixed-center symmetry-map row-layout transition failed'
    endif
    call move_alloc(reindexed_global_symmetry_map,global_symmetry_map)
    call move_alloc(reindexed_fixed_center_symmetry_map,fixed_center_symmetry_map)
    global_closed_core=ow_core_values
    ow_core_ids=initial_core_ids
    deallocate(initial_core_ids)
    call invert_ow_lattice(dc%system_tot%primitive_a,w90_lattice_inverse,w90_determinant,ok)
    if(.not.ok)error stop 'Wannier90 lattice is singular'
    w90_reciprocal_lattice=2d0*pi*transpose(w90_lattice_inverse)
    allocate(w90_atom_symbols(dc%system_tot%nion),w90_atoms_cart(3,dc%system_tot%nion))
    w90_atoms_cart=dc%system_tot%Rion
    do io=1,dc%system_tot%nion
      if(dc%system_tot%kion(io)<1.or.dc%system_tot%kion(io)>size(pp%atom_symbol))&
        error stop 'Wannier90 atom species is outside the pseudopotential table'
      w90_atom_symbols(io)=pp%atom_symbol(dc%system_tot%kion(io))
    enddo
    call setup_dg_w90_gamma_library(dc%icomm_tot,'overlapping_wannier_mlwf',&
      dc%system_tot%primitive_a,w90_reciprocal_lattice,w90_atom_symbols,w90_atoms_cart,&
      ntarget,ntarget,wannier_num_iter,dg_ow_w90_initial_projection,DG_W90_CONSTRAINED,&
      w90_nntot,w90_nncell,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'Wannier90 Gamma setup failed';endif
    allocate(w90_fractional(3,ncore),w90_eigenvalues(ntarget))
    do p=1,ncore
      w90_fractional(:,p)=[real(modulo(ow_core_ids(p)-1_8,int(dc%lg_tot%num(1),8)),8)/&
        real(dc%lg_tot%num(1),8),real(modulo((ow_core_ids(p)-1_8)/int(dc%lg_tot%num(1),8),&
        int(dc%lg_tot%num(2),8)),8)/real(dc%lg_tot%num(2),8),real((ow_core_ids(p)-1_8)/nxy8,8)/&
        real(dc%lg_tot%num(3),8)]
    enddo
    w90_eigenvalues=0d0
    call assemble_dg_w90_gamma_matrices(dc%icomm_tot,global_closed_core,w90_anchors,&
      ow_core_weights,w90_fractional,w90_nncell,w90_byte_limit,w90_m_matrix,w90_a_matrix,&
      w90_coordinator_bytes,w90_workspace_peak,ok,message,precomputed_a_matrix=w90_seed_a_matrix)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'Wannier90 M/A assembly failed';endif
    call fingerprint_ow_w90_matrices(dc%icomm_tot,w90_m_matrix,w90_a_matrix,&
      w90_input_fingerprint,ok)
    if(.not.ok)error stop 'Wannier90 M/A fingerprint failed'
    w90_input_fingerprint=ieor(w90_input_fingerprint,spectral_action_aggregate_fingerprint)
    if(w90_input_fingerprint==0_8)w90_input_fingerprint=1_8
    w90_replay_directory=''
    call get_environment_variable('SALMON_DG_W90_REPLAY_DIRECTORY',w90_replay_directory,&
      length=w90_replay_environment_length,status=w90_replay_environment_status,trim_name=.true.)
    w90_replay_enabled=merge(1,0,w90_replay_environment_status==0.and.w90_replay_environment_length>0)
    call MPI_Allreduce(w90_replay_enabled,w90_replay_enabled_min,1,MPI_INTEGER,MPI_MIN,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)error stop 'Wannier90 replay enable MIN agreement failed'
    call MPI_Allreduce(w90_replay_enabled,w90_replay_enabled_max,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.w90_replay_enabled_min/=w90_replay_enabled_max)&
      error stop 'Wannier90 replay enable state disagrees across ranks'
    if(w90_replay_enabled==1)then
      call export_dg_w90_replay_bundle(dc%icomm_tot,'.','overlapping_wannier_mlwf',&
        trim(w90_replay_directory),'overlapping_wannier_mlwf',DG_W90_CONSTRAINED,&
        w90_eigenvalues,w90_a_matrix,w90_m_matrix,w90_nncell,ok,message)
      if(.not.ok)then;write(0,'(a)')trim(message);error stop 'Wannier90 replay export failed';endif
    endif
    deallocate(w90_anchors,w90_fractional)
    call run_dg_w90_gamma_library(dc%icomm_tot,'overlapping_wannier_mlwf',&
      dc%system_tot%primitive_a,w90_reciprocal_lattice,w90_atom_symbols,w90_atoms_cart,&
      w90_m_matrix,w90_a_matrix,w90_eigenvalues,huge(1d0)/4d0,dg_ow_symmetry_tolerance,&
      wannier_num_iter,w90_transform,localized_centers,w90_spreads,w90_spread,ok,message,&
      localization_iterations)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'Wannier90 MLWF optimization failed';endif
    deallocate(w90_m_matrix,w90_a_matrix,w90_eigenvalues)
    deallocate(w90_atom_symbols,w90_atoms_cart,w90_nncell)
    w90_covariance_defect=0d0;w90_covariance_workspace=0_8
    do fixed_center_operation=1,fixed_center_group_order
      call assemble_dg_distributed_basis_symmetry_overlap_rows(dc%icomm_tot,global_closed_core,&
        ow_core_weights,fixed_center_symmetry_map(:,fixed_center_operation:fixed_center_operation),&
        fixed_center_row_ids,fixed_center_rows,&
        fixed_center_operation_workspace,ok,message)
      if(.not.ok)then
        write(0,'(a)')trim(message);error stop 'post-Wannier fixed-center assembly failed'
      endif
      call gather_dg_single_symmetry_representation(dc%icomm_tot,fixed_center_row_ids,fixed_center_rows,&
        1,0,fixed_center_representation,fixed_center_operation_workspace,ok,message)
      if(.not.ok)then;write(0,'(a)')trim(message);error stop 'post-Wannier generator gather failed';endif
      allocation_status=0
      if(rank==0)then
        call convert_sawf_pullback_to_active_representation(fixed_center_representation,ok,message)
        allocate(w90_seed_representation(ntarget,ntarget),stat=allocation_status)
        if(allocation_status==0.and.ok)w90_seed_representation=&
          matmul(conjg(transpose(w90_seed_a_matrix)),matmul(fixed_center_representation,w90_seed_a_matrix))
      endif
      call MPI_Bcast(ok,1,MPI_LOGICAL,0,dc%icomm_tot,ierr)
      call MPI_Bcast(allocation_status,1,MPI_INTEGER,0,dc%icomm_tot,ierr)
      if(ierr/=MPI_SUCCESS.or..not.ok.or.allocation_status/=0)&
        error stop 'post-Wannier target representation construction failed'
      call validate_dg_w90_generator_covariance(dc%icomm_tot,w90_transform,fixed_center_representation,&
        w90_seed_representation,dg_ow_symmetry_tolerance,w90_generator_covariance_defect,&
        fixed_center_operation_workspace,ok,message)
      if(.not.ok)then
        write(0,'(a)')trim(message);error stop 'post-Wannier fixed-center covariance failed'
      endif
      w90_covariance_defect=max(w90_covariance_defect,w90_generator_covariance_defect)
      w90_covariance_workspace=max(w90_covariance_workspace,fixed_center_operation_workspace)
      deallocate(fixed_center_row_ids,fixed_center_rows,fixed_center_representation)
      if(rank==0)deallocate(w90_seed_representation)
    enddo
    deallocate(w90_seed_a_matrix)
    if(rank==0)write(*,'(a,a,es16.8,a,i0)')'[OW-GS-DIAGNOSTIC] Wannier90 fixed-center covariance passed',&
      ' defect=',w90_covariance_defect,' workspace_peak_bytes=',w90_covariance_workspace
    localized_centers=matmul(w90_lattice_inverse,localized_centers)
    call apply_dg_w90_gamma_transform(dc%icomm_tot,ow_core_ids,ow_core_values,&
      transform=w90_transform,centers=localized_centers,tolerance=dg_ow_symmetry_tolerance,&
      ok=ok,message=message,spreads=w90_spreads)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'Wannier90 MLWF gauge canonicalization failed';endif
    deallocate(w90_spreads)
    call fingerprint_ow_w90_transform(dc%icomm_tot,w90_transform,w90_transform_fingerprint,ok)
    if(.not.ok)error stop 'Wannier90 canonical transform fingerprint failed'
    ! Column permutation and diagonal pivot phases are a unitary gauge P.  The
    ! corresponding symmetry representation is P^H D P and its covariance
    ! residual is R P, so the accepted max-entry defect above is unchanged.
    ! Downstream character-sector construction rebuilds representations in the
    ! canonicalized spatial basis rather than retaining the pre-gauge matrices.
    call inherit_dg_w90_affine_receipts(w90_transform,&
      max(maxval(lcfo_total_symmetry_residual),w90_covariance_defect),&
      dg_ow_symmetry_tolerance,&
      w90_identity_defect,w90_unitarity_defect,w90_closure_defect,w90_symmetry_workspace_peak,&
      ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'post-MLWF row-owned group validation failed';endif
    global_retained_group_closure_defect=w90_closure_defect
    localization_initial_spread=w90_spread(1);localization_final_spread=w90_spread(1)
    localization_maximum_gradient=0d0
    localization_spread_evaluations=0;localization_converged=.true.

#ifdef USE_EIGENEXA
    allocate(translation_w90_rows(size(translation_row_ids),ntarget),&
      translation_lcfo_rows(size(translation_row_ids),ntarget),&
      translation_w90_values(ntarget),translation_lcfo_values(ntarget),stat=allocation_status)
    call MPI_Allreduce(allocation_status,translation_allocation_status,1,MPI_INTEGER,MPI_MAX,&
      dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.translation_allocation_status/=0)&
      error stop 'translation reference operator allocation failed collectively'
    do p=1,size(translation_row_ids)
      translation_w90_rows(p,:)=w90_transform(int(translation_row_ids(p)),:)
    enddo
    deallocate(w90_transform)
    translation_lcfo_rows=(0d0,0d0)
    allocate(translation_lcfo_local_row(ntarget),translation_lcfo_global_row(ntarget))
    do i=1,ntarget
      do io=1,ntarget
        translation_lcfo_local_row(io)=sum(ow_core_weights*&
          conjg(global_closed_core(i,:))*lcfo_reference_core(io,:))
      enddo
      call MPI_Allreduce(translation_lcfo_local_row,translation_lcfo_global_row,ntarget,&
        MPI_DOUBLE_COMPLEX,MPI_SUM,dc%icomm_tot,ierr)
      if(ierr/=MPI_SUCCESS)error stop 'LCFO source overlap reduction failed'
      p=findloc(translation_row_ids,int(i,8),dim=1)
      if(p>0)translation_lcfo_rows(p,:)=translation_lcfo_global_row
    enddo
    deallocate(translation_lcfo_local_row,translation_lcfo_global_row,lcfo_reference_core)
    translation_lcfo_fingerprint=ieor(occupied_composition_fingerprint,projector_composition_fingerprint)
    if(translation_lcfo_fingerprint==0_8)translation_lcfo_fingerprint=1_8
    do io=1,ntarget
      translation_w90_values(io)=modulo(localized_centers(1,io),1d0)+&
        sqrt(2d0)*modulo(localized_centers(2,io),1d0)+sqrt(3d0)*modulo(localized_centers(3,io),1d0)
      if(io<=nstate)then
        translation_lcfo_values(io)=lcfo_retained_eigenvalues(io)
      else
        i=io-nstate
        translation_lcfo_values(io)=modulo(dot_product(matmul(w90_lattice_inverse,&
          dc%system_tot%Rion(:,manifest_channels(i)%atom)),[1d0,sqrt(2d0),sqrt(3d0)]),1d0)+&
          0.01d0*real(manifest_channels(i)%l,8)+0.001d0*real(manifest_channels(i)%m,8)+&
          0.0001d0*real(manifest_channels(i)%radial,8)
      endif
    enddo
    deallocate(localized_centers)
    ow_saved_eigenexa_comm=info%icomm_o
    call finalize_eigenexa(info);info%icomm_o=dc%icomm_tot
    call init_eigenexa_mod(info,2*ntarget,direct_block_only=.true.)
    call split_dg_translation_character_sector_eigenexa(info,dc%icomm_tot,translation_row_ids,&
      translation_generator_rows,translation_gamma_rows,translation_generator_characters,&
      translation_generator_orders,translation_element_words,translation_character_conjugates,1,&
      dg_ow_symmetry_tolerance,translation_character_fingerprint,translation_sector_rows,translation_sector_rank,&
      translation_identity_defect,translation_unitarity_defect,translation_commutator_defect,&
      translation_order_defect,translation_gamma_pairing_defect,translation_sector_fingerprint,&
      translation_sector_workspace_peak,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'post-W90 reference character split failed';endif
    call project_dg_w90_reference_sector_operators(dc%icomm_tot,translation_row_ids,translation_sector_rows,&
      translation_w90_rows,translation_w90_values,translation_lcfo_rows,translation_lcfo_values,ntarget,&
      w90_transform_fingerprint,translation_lcfo_fingerprint,w90_unitarity_defect,0d0,&
      dg_ow_symmetry_tolerance,translation_w90_operator,translation_lcfo_operator,translation_operator_defect,&
      translation_operator_fingerprint,translation_operator_workspace,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'reference-sector physical operator projection failed';endif
    deallocate(translation_w90_rows,translation_lcfo_rows,translation_w90_values,translation_lcfo_values)
    call anchor_dg_w90_reference_character_sector(dc%icomm_tot,translation_row_ids,translation_sector_rows,&
      translation_w90_operator,translation_lcfo_operator,ntarget,w90_transform_fingerprint,translation_lcfo_fingerprint,&
      w90_unitarity_defect,translation_operator_defect,dg_ow_symmetry_tolerance,translation_reference_rows,&
      translation_anchor_defect,translation_anchor_fingerprint,translation_anchor_workspace,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'reference translation character W90 anchor failed';endif
    allocate(translation_spatial_ids(ncore),stat=allocation_status)
    call MPI_Allreduce(allocation_status,translation_allocation_status,1,MPI_INTEGER,MPI_MAX,&
      dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.translation_allocation_status/=0)&
      error stop 'translation periodic-position allocation failed collectively'
    do p=1,ncore
      translation_spatial_ids(p)=int(rank*ncore+p,8)
    enddo
    call materialize_dg_row_owned_sector_on_spatial_grid(dc%icomm_tot,translation_row_ids,ntarget,&
      translation_reference_rows,global_closed_core,w90_input_fingerprint,translation_reference_spatial,&
      translation_materialize_fingerprint,translation_materialize_workspace,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'reference character spatial materialization failed';endif
    deallocate(translation_sector_rows,translation_w90_operator,translation_lcfo_operator,translation_reference_rows)

    allocate(translation_generator_maps(ncore,&
      translation_character_generator_count),translation_character_done(size(translation_characters,1)),stat=allocation_status)
    call MPI_Allreduce(allocation_status,translation_allocation_status,1,MPI_INTEGER,MPI_MAX,&
      dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.translation_allocation_status/=0)&
      error stop 'translation spatial action allocation failed collectively'
    do i=1,translation_character_generator_count
      io=global_translation_subgroup(translation_canonical_operations(translation_character_generators(i)))
      translation_generator_maps(:,i)=global_symmetry_map(:,io)
    enddo
    call prepare_dg_translation_character_action(dc%icomm_tot,translation_spatial_ids,&
      translation_global_core_count,translation_generator_maps,translation_generator_orders,&
      translation_element_words,translation_canonical_product,1,translation_character_fingerprint,&
      dg_ow_symmetry_tolerance,translation_prepared_action,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'translation action preparation failed';endif
    translation_gamma_fingerprint=int(z'6A09E667F3BCC909',8)
    translation_gamma_fingerprint=ieor(ishftc(translation_gamma_fingerprint,9),translation_global_core_count8)
    if(translation_gamma_fingerprint==0_8)translation_gamma_fingerprint=1_8
    translation_character_done=.false.;translation_processed_count=0
    translation_alignment_max_defect=0d0;translation_gamma_max_defect=0d0
    translation_post_gauge_fingerprint=ieor(translation_operator_fingerprint,translation_anchor_fingerprint)
    translation_post_gauge_fingerprint=ieor(ishftc(translation_post_gauge_fingerprint,11),&
      translation_materialize_fingerprint)
    do translation_character=1,size(translation_characters,1)
      if(translation_character_done(translation_character))cycle
      if(translation_character==1)then
        allocate(translation_aligned_spatial,source=translation_reference_spatial)
      else
        call split_dg_translation_character_sector_eigenexa(info,dc%icomm_tot,translation_row_ids,&
          translation_generator_rows,translation_gamma_rows,translation_generator_characters,&
          translation_generator_orders,translation_element_words,translation_character_conjugates,&
          translation_character,dg_ow_symmetry_tolerance,translation_character_fingerprint,&
          translation_sector_rows,translation_sector_rank,translation_identity_defect,translation_unitarity_defect,&
          translation_commutator_defect,translation_order_defect,translation_gamma_pairing_defect,&
          translation_sector_fingerprint,translation_sector_workspace_peak,ok,message)
        if(.not.ok)then;write(0,'(a)')trim(message);error stop 'target translation character split failed';endif
        call materialize_dg_row_owned_sector_on_spatial_grid(dc%icomm_tot,translation_row_ids,ntarget,&
          translation_sector_rows,global_closed_core,w90_input_fingerprint,translation_target_spatial,&
          translation_materialize_fingerprint,translation_materialize_workspace,ok,message)
        deallocate(translation_sector_rows)
        if(.not.ok)then;write(0,'(a)')trim(message);error stop 'target character spatial materialization failed';endif
        call build_dg_translation_character_intertwining_phase_prepared(dc%icomm_tot,translation_prepared_action,&
          translation_characters(1,:),&
          translation_characters(translation_character,:),dg_ow_symmetry_tolerance,translation_phase,translation_phase_fingerprint,&
          translation_phase_payload_fingerprint,translation_phase_workspace,ok,message)
        if(.not.ok)then;write(0,'(a)')trim(message);error stop 'translation character intertwining phase failed';endif
        call align_dg_w90_character_sectors_by_periodic_phase(dc%icomm_tot,translation_spatial_ids,&
          translation_reference_spatial,translation_target_spatial,translation_phase,translation_global_core_count,&
          translation_phase_fingerprint,translation_phase_payload_fingerprint,dg_ow_symmetry_tolerance,&
          translation_aligned_spatial,translation_singular_values,translation_alignment_defect,&
          translation_alignment_fingerprint,translation_alignment_workspace,ok,message,ow_core_weights)
        deallocate(translation_target_spatial,translation_phase,translation_singular_values)
        if(.not.ok)then;write(0,'(a)')trim(message);error stop 'periodic-phase character alignment failed';endif
        translation_alignment_max_defect=max(translation_alignment_max_defect,translation_alignment_defect)
        translation_post_gauge_fingerprint=ieor(ishftc(translation_post_gauge_fingerprint,11),&
          translation_alignment_fingerprint)
      endif
      translation_partner=translation_character_conjugates(translation_character)
      translation_self_conjugate=translation_partner==translation_character
      if(translation_self_conjugate)then
        allocate(translation_conjugate_spatial,source=translation_aligned_spatial)
      else
        call split_dg_translation_character_sector_eigenexa(info,dc%icomm_tot,translation_row_ids,&
          translation_generator_rows,translation_gamma_rows,translation_generator_characters,&
          translation_generator_orders,translation_element_words,translation_character_conjugates,&
          translation_partner,dg_ow_symmetry_tolerance,translation_character_fingerprint,&
          translation_sector_rows,translation_sector_rank,translation_identity_defect,translation_unitarity_defect,&
          translation_commutator_defect,translation_order_defect,translation_gamma_pairing_defect,&
          translation_sector_fingerprint,translation_sector_workspace_peak,ok,message)
        if(.not.ok)then;write(0,'(a)')trim(message);error stop 'conjugate translation character split failed';endif
        call materialize_dg_row_owned_sector_on_spatial_grid(dc%icomm_tot,translation_row_ids,ntarget,&
          translation_sector_rows,global_closed_core,w90_input_fingerprint,translation_conjugate_spatial,&
          translation_materialize_fingerprint,translation_materialize_workspace,ok,message)
        deallocate(translation_sector_rows)
        if(.not.ok)then;write(0,'(a)')trim(message);error stop 'conjugate character spatial materialization failed';endif
      endif
      allocate(translation_spatial_gamma_rows(ncore,0),stat=allocation_status)
      call MPI_Allreduce(allocation_status,translation_allocation_status,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
      if(ierr/=MPI_SUCCESS.or.translation_allocation_status/=0)error stop 'implicit spatial Gamma allocation failed'
      call sew_dg_w90_periodic_phase_conjugate_sector(dc%icomm_tot,translation_spatial_ids,&
        translation_aligned_spatial,translation_spatial_gamma_rows,translation_conjugate_spatial,&
        translation_global_core_count,translation_gamma_fingerprint,translation_self_conjugate,0d0,&
        dg_ow_symmetry_tolerance,translation_target_spatial,translation_gamma_defect,&
        translation_gamma_workspace,ok,message,.true.,ow_core_weights)
      deallocate(translation_spatial_gamma_rows,translation_conjugate_spatial)
      if(.not.ok)then;write(0,'(a)')trim(message);error stop 'spatial Gamma character sewing failed';endif
      translation_gamma_max_defect=max(translation_gamma_max_defect,translation_gamma_defect)
      translation_post_gauge_fingerprint=ieor(ishftc(translation_post_gauge_fingerprint,11),&
        ieor(translation_gamma_fingerprint,int(translation_character,8)))
      if(translation_self_conjugate)then
        deallocate(translation_aligned_spatial);call move_alloc(translation_target_spatial,translation_aligned_spatial)
      endif
      call accumulate_dg_translation_character_orbit_sector_values(dc%icomm_tot,translation_inverse_state,&
        translation_character,translation_characters,translation_character_fingerprint,translation_aligned_spatial,&
        translation_processed_count==0,translation_processed_count+1==size(translation_characters,1),&
        dg_ow_symmetry_tolerance,translation_orbit_rows,translation_inverse_workspace,ok,message)
      if(.not.ok)then;write(0,'(a)')trim(message);error stop 'translation inverse accumulation failed';endif
      translation_character_done(translation_character)=.true.;translation_processed_count=translation_processed_count+1
      deallocate(translation_aligned_spatial)
      if(.not.translation_self_conjugate)then
        call accumulate_dg_translation_character_orbit_sector_values(dc%icomm_tot,translation_inverse_state,&
          translation_partner,translation_characters,translation_character_fingerprint,translation_target_spatial,&
          .false.,translation_processed_count+1==size(translation_characters,1),dg_ow_symmetry_tolerance,&
          translation_orbit_rows,translation_inverse_workspace,ok,message)
        if(.not.ok)then;write(0,'(a)')trim(message);error stop 'conjugate inverse accumulation failed';endif
        translation_character_done(translation_partner)=.true.;translation_processed_count=translation_processed_count+1
        deallocate(translation_target_spatial)
      endif
    enddo
    call release_dg_prepared_translation_action(translation_prepared_action)
    call finalize_eigenexa(info);info%icomm_o=ow_saved_eigenexa_comm;call init_eigenexa_mod(info,system%no)
    if(translation_processed_count/=size(translation_characters,1))error stop 'translation character schedule incomplete'
    deallocate(ow_core_values);call move_alloc(translation_orbit_rows,ow_core_values)
    deallocate(global_closed_core)
    deallocate(translation_reference_spatial,translation_generator_maps,translation_spatial_ids,&
      translation_character_done,translation_generator_rows,translation_row_ids,translation_gamma_rows,&
      translation_canonical_product,translation_product)
    call validate_dg_factored_point_cogroup_gauge(dc%icomm_tot,ow_core_values,ow_core_weights,&
      global_symmetry_map(:,global_point_representatives),global_symmetry_map(:,global_translation_subgroup),&
      global_point_cogroup_product,global_point_cogroup_identity_operation,&
      global_translation_cocycle,size(translation_characters,1),max(translation_alignment_max_defect,&
      translation_gamma_max_defect),dg_ow_symmetry_tolerance,translation_identity_defect,&
      translation_unitarity_defect,translation_closure_defect,translation_transform_workspace,ok,message,&
      translation_point_generator_count,translation_point_checked_pair_count)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'post-character point-cogroup cocycle proof failed';endif
    if(rank==0)write(*,'(a,2(a,i0),2(a,es16.8),a,i0)')'[OW-GS-DIAGNOSTIC] factored_point_cogroup_proof',&
      ' generator_count=',translation_point_generator_count,' checked_dense_pairs=',translation_point_checked_pair_count,&
      ' identity_defect=',translation_identity_defect,' closure_defect=',translation_closure_defect,&
      ' workspace_peak_bytes=',translation_transform_workspace
    global_retained_group_closure_defect=max(global_retained_group_closure_defect,translation_closure_defect)
    allocate(localized_centers(3,ntarget),localized_center_magnitudes(3,ntarget))
    call compute_dg_periodic_wannier_centers(dc%icomm_tot,ow_core_values,ow_core_weights,&
      core_periodic_phase,localized_centers,localized_center_magnitudes,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'localized Wannier center measurement failed';end if
    if(rank==0)write(*,'(a,2(a,es12.4))')'[OW-GS-DIAGNOSTIC] periodic_center_magnitude',&
      ' minimum=',minval(localized_center_magnitudes),' maximum=',maxval(localized_center_magnitudes)
    call verify_dg_wannier_center_affine_orbits(localized_centers,global_point_integer_rotations,&
      global_point_fractional_translations,retained_closure_search_tolerance,ok,message,&
      moment_magnitudes=localized_center_magnitudes,failed_operation=failed_operation)
    if(.not.ok)then
      center_failure_message=message
      call diagnose_dg_point_center_gauge(dc%icomm_tot,ow_core_values,ow_core_weights,&
        global_symmetry_map(:,failed_operation),global_point_integer_rotations(:,:,failed_operation),&
        global_point_fractional_translations(:,failed_operation),localized_centers,&
        retained_closure_search_tolerance,monomial_defect,center_block_leakage,&
        center_representation_unitarity_defect,center_gauge_workspace_peak,center_diagnostic_ok,&
        center_diagnostic_message)
      if(rank==0.and.center_diagnostic_ok)write(*,'(a,i0,3(a,es16.8),a,i0)')&
        '[OW-GS-DIAGNOSTIC] point_center_gauge failed_operation=',failed_operation,&
        ' monomial_defect=',monomial_defect,' center_block_leakage=',center_block_leakage,&
        ' representation_unitarity_defect=',center_representation_unitarity_defect,&
        ' workspace_peak_bytes=',center_gauge_workspace_peak
      if(.not.center_diagnostic_ok)then
        if(rank==0)write(0,'(2a)')'point center-gauge diagnostic failed: ',trim(center_diagnostic_message)
        error stop 'point center-gauge diagnostic failed'
      endif
      if(rank==0)write(*,'(2a)')&
        '[OW-GS-DIAGNOSTIC] nonmonomial center action retained: ',trim(center_failure_message)
    end if
    allocate(ow_reciprocal_operation_maps,source=global_symmetry_map(:,global_point_representatives),&
      stat=allocation_status)
    call comm_logical_and(allocation_status==0,reusable,dc%icomm_tot)
    if(.not.reusable)error stop 'physical reciprocal-operation map allocation failed'
    allocate(ow_pencil_generator_maps,source=global_symmetry_map(:,global_affine_generators),&
      stat=allocation_status)
    call comm_logical_and(allocation_status==0,reusable,dc%icomm_tot)
    if(.not.reusable)error stop 'compact pencil-generator map allocation failed'
    deallocate(global_symmetry_map)
    call materialize_ow_distributed_core_to_buffer(dc%icomm_tot,ow_core_values,ow_core_ids,&
      physical_ids,ow_box_values,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'post-character core-to-buffer streaming failed';endif
    allocate(ow_direct_core_gradients(3,ntarget,ncore),ow_neighbor_plus_ids(ncore),&
      ow_neighbor_minus_ids(ncore));ow_direct_core_gradients=(0d0,0d0)
    do ix=1,3
      do gradient_distance=1,size(stencil%coef_nab,1)
        do core_index=1,ncore
          raw_ix=int(modulo(ow_core_ids(core_index)-1_8,int(dc%lg_tot%num(1),8)))
          raw_iy=int(modulo((ow_core_ids(core_index)-1_8)/int(dc%lg_tot%num(1),8),&
            int(dc%lg_tot%num(2),8)))
          raw_iz=int((ow_core_ids(core_index)-1_8)/nxy8)
          select case(ix)
          case(1)
            ow_neighbor_minus_ids(core_index)=1_8+int(modulo(raw_ix-gradient_distance,&
              dc%lg_tot%num(1)),8)+int(dc%lg_tot%num(1),8)*(int(raw_iy,8)+int(dc%lg_tot%num(2),8)*int(raw_iz,8))
            ow_neighbor_plus_ids(core_index)=1_8+int(modulo(raw_ix+gradient_distance,&
              dc%lg_tot%num(1)),8)+int(dc%lg_tot%num(1),8)*(int(raw_iy,8)+int(dc%lg_tot%num(2),8)*int(raw_iz,8))
          case(2)
            ow_neighbor_minus_ids(core_index)=1_8+int(raw_ix,8)+int(dc%lg_tot%num(1),8)*(&
              int(modulo(raw_iy-gradient_distance,dc%lg_tot%num(2)),8)+&
              int(dc%lg_tot%num(2),8)*int(raw_iz,8))
            ow_neighbor_plus_ids(core_index)=1_8+int(raw_ix,8)+int(dc%lg_tot%num(1),8)*(&
              int(modulo(raw_iy+gradient_distance,dc%lg_tot%num(2)),8)+&
              int(dc%lg_tot%num(2),8)*int(raw_iz,8))
          case default
            ow_neighbor_minus_ids(core_index)=1_8+int(raw_ix,8)+int(dc%lg_tot%num(1),8)*(&
              int(raw_iy,8)+int(dc%lg_tot%num(2),8)*&
              int(modulo(raw_iz-gradient_distance,dc%lg_tot%num(3)),8))
            ow_neighbor_plus_ids(core_index)=1_8+int(raw_ix,8)+int(dc%lg_tot%num(1),8)*(&
              int(raw_iy,8)+int(dc%lg_tot%num(2),8)*&
              int(modulo(raw_iz+gradient_distance,dc%lg_tot%num(3)),8))
          end select
        enddo
        call materialize_ow_distributed_core_to_buffer(dc%icomm_tot,ow_core_values,ow_core_ids,&
          ow_neighbor_plus_ids,ow_neighbor_plus_values,ok,message)
        if(ok)call materialize_ow_distributed_core_to_buffer(dc%icomm_tot,ow_core_values,ow_core_ids,&
          ow_neighbor_minus_ids,ow_neighbor_minus_values,ok,message)
        if(.not.ok)then;write(0,'(a)')trim(message);error stop 'direct core gradient neighbor stream failed';endif
        ow_direct_core_gradients(ix,:,:)=ow_direct_core_gradients(ix,:,:)+&
          stencil%coef_nab(gradient_distance,ix)*(ow_neighbor_plus_values-ow_neighbor_minus_values)
        deallocate(ow_neighbor_plus_values,ow_neighbor_minus_values)
      enddo
    enddo
    deallocate(ow_neighbor_plus_ids,ow_neighbor_minus_ids)
    deallocate(ow_core_values)
    allocate(ow_box_gradients(3,ntarget,nbox));call periodic_box_gradients(ow_box_values,ow_box_size,&
      stencil%coef_nab,ow_box_gradients)
#endif
    ! HYBRID_CONSTRAINED_LEGACY_ARM_END
    endif
    ! HYBRID_LOCALIZATION_FIRST_BRANCH_END
    if(.not.allocated(ow_pencil_generator_maps))then
      allocate(ow_reciprocal_operation_maps,source=global_symmetry_map(:,global_point_representatives),&
        stat=allocation_status)
      call comm_logical_and(allocation_status==0,reusable,dc%icomm_tot)
      if(.not.reusable)error stop 'physical reciprocal-operation map allocation failed'
      allocate(ow_pencil_generator_maps,source=global_symmetry_map(:,global_affine_generators),&
        stat=allocation_status)
      call comm_logical_and(allocation_status==0,reusable,dc%icomm_tot)
      if(.not.reusable)error stop 'compact pencil-generator map allocation failed'
      deallocate(global_symmetry_map)
    endif
    if(allocated(global_closed_core))deallocate(global_closed_core)
    if(allocated(lcfo_reference_core))deallocate(lcfo_reference_core)
    if(allocated(w90_transform))deallocate(w90_transform)
    if(rank==0)write(*,'(a,5(a,es12.4),3(a,i0))')'[OW-GS-DIAGNOSTIC] Wannier90_MLWF',&
      ' gauge_spread=',w90_spread(3),' total_spread=',w90_spread(1),&
      ' identity_defect=',w90_identity_defect,' unitarity_defect=',w90_unitarity_defect,&
      ' closure_defect=',w90_closure_defect,' coordinator_bytes=',w90_coordinator_bytes,&
      ' workspace_peak_bytes=',w90_workspace_peak,&
      ' symmetry_workspace_peak_bytes=',w90_symmetry_workspace_peak
    closure_residual=global_retained_group_closure_defect
    local_exact_symmetry_fingerprint=fingerprint_dg_exact_fragment_symmetry(&
      global_point_integer_rotations,global_point_product,dg_ow_symmetry_tolerance,&
      global_point_fractional_translations)
    if(local_exact_symmetry_fingerprint==0_8)&
      error stop 'invalid exact fragment symmetry checkpoint evidence'
    allocate(exact_fragment_symmetry_fingerprints(nproc))
    call MPI_Allgather(local_exact_symmetry_fingerprint,1,MPI_INTEGER8,&
      exact_fragment_symmetry_fingerprints,1,MPI_INTEGER8,dc%icomm_tot,ierr)
    ow_symmetry_fingerprint=0_8
    do p=1,nproc
      ow_symmetry_fingerprint=ieor(ow_symmetry_fingerprint,ishftc(&
        ieor(exact_fragment_symmetry_fingerprints(p),int(p,8)),modulo(13*p,63)))
    end do
    deallocate(exact_fragment_symmetry_fingerprints)
    ow_symmetry_fingerprint=ieor(ishftc(ow_symmetry_fingerprint,17),translation_post_gauge_fingerprint)
    if(allocated(ow_core_values))deallocate(ow_core_values)
    allocate(ow_core_values(ntarget,ncore),ow_core_gradients(3,ntarget,ncore))
    do core_index=1,ncore
      p=ow_core_box_positions(core_index)
      ow_core_values(:,core_index)=ow_box_values(:,p)
      ow_core_gradients(:,:,core_index)=ow_box_gradients(:,:,p)
    end do
    ow_gradient_path_local=0d0
    do core_index=1,ncore
      ow_gradient_path_local(1)=ow_gradient_path_local(1)+ow_core_weights(core_index)*&
        sum(abs(ow_core_gradients(:,:,core_index)-ow_direct_core_gradients(:,:,core_index))**2)
      ow_gradient_path_local(2)=ow_gradient_path_local(2)+ow_core_weights(core_index)*&
        sum(abs(ow_direct_core_gradients(:,:,core_index))**2)
    enddo
    call MPI_Allreduce(ow_gradient_path_local,ow_gradient_path_global,2,MPI_DOUBLE_PRECISION,&
      MPI_SUM,dc%icomm_tot,ierr)
    if(rank==0)write(*,'(2(a,es16.8))')'[OW-GS-DIAGNOSTIC] buffer/direct gradient relative defect=',&
      sqrt(max(0d0,ow_gradient_path_global(1))/max(tiny(1d0),ow_gradient_path_global(2))),&
      ' direct gradient norm=',sqrt(max(0d0,ow_gradient_path_global(2)))
    allocate(ow_scalar_probe(1,ncore),ow_vector_probe(3,1,ncore),&
      ow_scalar_representation(1,1,size(ow_pencil_generator_maps,2)),ow_scalar_probe_weights(ncore))
    ow_scalar_representation=(1d0,0d0);ow_scalar_probe_weights=1d0
    ow_scalar_probe(1,:)=cmplx(ow_partition_weight(ow_core_box_positions),0d0,8)
    ow_vector_probe(:,1,:)=cmplx(ow_partition_gradient(:,ow_core_box_positions),0d0,8)
    if(size(ow_pencil_generator_maps,2)>0)then
    call assemble_dg_distributed_basis_symmetry_overlap(dc%icomm_tot,ow_core_values,ow_core_weights,&
      ow_pencil_generator_maps,ow_pencil_generator_representation,ok,message)
    if(ok)call measure_dg_spatial_basis_covariance(dc%icomm_tot,ow_core_values,ow_core_weights,&
      ow_pencil_generator_maps,ow_pencil_generator_representation,&
      ow_core_spatial_covariance_residual,ok,message)
    if(ok.and.rank==0)write(*,'(a,es16.8,a,i0)')&
      '[OW-GS-DIAGNOSTIC] core spatial generator covariance max=',&
      maxval(ow_core_spatial_covariance_residual),' operation=',&
      maxloc(ow_core_spatial_covariance_residual,dim=1)
    if(ok)then
      ow_spatial_covariance_relative=maxval(ow_core_spatial_covariance_residual)
      ow_spatial_covariance_absolute=ow_spatial_covariance_relative*sqrt(real(ntarget,8))
    endif
    if(allocated(ow_core_spatial_covariance_residual))deallocate(ow_core_spatial_covariance_residual)
    call measure_dg_spatial_basis_covariance(dc%icomm_tot,ow_scalar_probe,ow_scalar_probe_weights,&
      ow_pencil_generator_maps,ow_scalar_representation,ow_scalar_probe_residual,ok,message)
    if(ok.and.rank==0)write(*,'(a,es16.8,a,i0)')&
      '[OW-GS-DIAGNOSTIC] core partition-weight covariance max=',maxval(ow_scalar_probe_residual),&
      ' operation=',maxloc(ow_scalar_probe_residual,dim=1)
    if(allocated(ow_scalar_probe_residual))deallocate(ow_scalar_probe_residual)
    call measure_dg_spatial_gradient_covariance(dc%icomm_tot,ow_vector_probe,ow_scalar_probe_weights,&
      ow_pencil_generator_maps,ow_scalar_representation,&
      global_point_rotations(:,:,global_affine_generators),&
      ow_gradient_covariance_left,ow_gradient_covariance_transpose,ok,message)
    if(ok.and.rank==0)write(*,'(2(a,es16.8,a,i0))')&
      '[OW-GS-DIAGNOSTIC] core partition-gradient covariance R max=',&
      maxval(ow_gradient_covariance_left),' operation=',maxloc(ow_gradient_covariance_left,dim=1),&
      ' RT max=',maxval(ow_gradient_covariance_transpose),&
      ' operation=',maxloc(ow_gradient_covariance_transpose,dim=1)
    if(allocated(ow_gradient_covariance_left))deallocate(ow_gradient_covariance_left)
    if(allocated(ow_gradient_covariance_transpose))deallocate(ow_gradient_covariance_transpose)
    if(ok)call measure_dg_spatial_gradient_covariance(dc%icomm_tot,ow_core_gradients,ow_core_weights,&
      ow_pencil_generator_maps,ow_pencil_generator_representation,&
      global_point_rotations(:,:,global_affine_generators),&
      ow_gradient_covariance_left,ow_gradient_covariance_transpose,ok,message,&
      orbital_action_residual=ow_gradient_covariance_candidates)
    if(ok.and.rank==0)write(*,'(2(a,es16.8,a,i0))')&
      '[OW-GS-DIAGNOSTIC] core gradient covariance R max=',maxval(ow_gradient_covariance_left),&
      ' operation=',maxloc(ow_gradient_covariance_left,dim=1),&
      ' RT max=',maxval(ow_gradient_covariance_transpose),&
      ' operation=',maxloc(ow_gradient_covariance_transpose,dim=1)
    if(ok.and.rank==0)then
      ow_gradient_stencil_norm_bound=2d0*maxval(sum(abs(stencil%coef_nab),dim=1))
      ow_gradient_covariance_absolute=max(maxval(ow_gradient_covariance_left),&
        maxval(ow_gradient_covariance_transpose))*sqrt(max(0d0,ow_gradient_path_global(2)))
      write(*,'(4(a,es16.8))')'[OW-GS-DIAGNOSTIC] gradient covariance absolute=',&
        ow_gradient_covariance_absolute,' differentiated basis-error bound=',&
        ow_gradient_stencil_norm_bound*ow_spatial_covariance_absolute,&
        ' stencil operator norm bound=',ow_gradient_stencil_norm_bound,&
        ' spatial covariance absolute=',ow_spatial_covariance_absolute
    endif
    if(ok.and.rank==0)then
      write(*,'(2(a,es16.8,a,i0))')'[OW-GS-DIAGNOSTIC] gradient C^T: R=',&
        maxval(ow_gradient_covariance_candidates(1,:)),' op=',&
        maxloc(ow_gradient_covariance_candidates(1,:),dim=1),' RT=',&
        maxval(ow_gradient_covariance_candidates(2,:)),' op=',&
        maxloc(ow_gradient_covariance_candidates(2,:),dim=1)
      write(*,'(2(a,es16.8,a,i0))')'[OW-GS-DIAGNOSTIC] gradient C: R=',&
        maxval(ow_gradient_covariance_candidates(3,:)),' op=',&
        maxloc(ow_gradient_covariance_candidates(3,:),dim=1),' RT=',&
        maxval(ow_gradient_covariance_candidates(4,:)),' op=',&
        maxloc(ow_gradient_covariance_candidates(4,:),dim=1)
      write(*,'(2(a,es16.8,a,i0))')'[OW-GS-DIAGNOSTIC] gradient C*: R=',&
        maxval(ow_gradient_covariance_candidates(5,:)),' op=',&
        maxloc(ow_gradient_covariance_candidates(5,:),dim=1),' RT=',&
        maxval(ow_gradient_covariance_candidates(6,:)),' op=',&
        maxloc(ow_gradient_covariance_candidates(6,:),dim=1)
      write(*,'(2(a,es16.8,a,i0))')'[OW-GS-DIAGNOSTIC] gradient C^H: R=',&
        maxval(ow_gradient_covariance_candidates(7,:)),' op=',&
        maxloc(ow_gradient_covariance_candidates(7,:),dim=1),' RT=',&
        maxval(ow_gradient_covariance_candidates(8,:)),' op=',&
        maxloc(ow_gradient_covariance_candidates(8,:),dim=1)
    endif
    else
      allocate(ow_pencil_generator_representation(ntarget,ntarget,0))
      ow_spatial_covariance_relative=0d0;ow_spatial_covariance_absolute=0d0
      ow_gradient_covariance_absolute=0d0;ow_gradient_stencil_norm_bound=0d0
      ok=.true.;message=''
      if(rank==0)then
        write(*,'(a,es16.8,a,i0)')'[OW-GS-DIAGNOSTIC] core spatial generator covariance max=',&
          0d0,' operation=',global_identity_operation
        write(*,'(a,es16.8,a,i0)')'[OW-GS-DIAGNOSTIC] core partition-weight covariance max=',&
          0d0,' operation=',global_identity_operation
        write(*,'(2(a,es16.8,a,i0))')'[OW-GS-DIAGNOSTIC] core partition-gradient covariance R max=',&
          0d0,' operation=',global_identity_operation,' RT max=',0d0,' operation=',global_identity_operation
        write(*,'(2(a,es16.8,a,i0))')'[OW-GS-DIAGNOSTIC] core gradient covariance R max=',&
          0d0,' operation=',global_identity_operation,' RT max=',0d0,' operation=',global_identity_operation
      endif
    endif
    allocate(ow_gradient_identity_map(ncore,1))
    ow_gradient_identity_map(:,1)=[(int(rank,8)*int(ncore,8)+int(p,8),p=1,ncore)]
    ow_gradient_identity_rotation=0d0
    do i=1,3;ow_gradient_identity_rotation(i,i,1)=1d0;enddo
    call measure_ow_discrete_gradient_map_commutator(dc%icomm_tot,ow_scalar_probe,&
      ow_vector_probe,ow_scalar_probe_weights,ow_core_ids,ow_gradient_identity_map,&
      dc%lg_tot%num,stencil%coef_nab,ow_gradient_identity_rotation,&
      ow_gradient_map_commutator,ok,message)
    if(ok.and.rank==0)write(*,'(a,es16.8)')&
      '[OW-GS-DIAGNOSTIC] partition analytic/finite-difference gradient defect=',&
      ow_gradient_map_commutator(1)
    if(allocated(ow_gradient_map_commutator))deallocate(ow_gradient_map_commutator)
    if(ok)call measure_ow_discrete_gradient_map_commutator(dc%icomm_tot,ow_core_values,&
      ow_direct_core_gradients,ow_core_weights,ow_core_ids,ow_gradient_identity_map,&
      dc%lg_tot%num,stencil%coef_nab,ow_gradient_identity_rotation,&
      ow_gradient_map_commutator,ok,message)
    if(ok.and.rank==0)write(*,'(a,es16.8)')&
      '[OW-GS-DIAGNOSTIC] finite-difference/value reconstruction defect=',&
      ow_gradient_map_commutator(1)
    if(allocated(ow_gradient_map_commutator))deallocate(ow_gradient_map_commutator)
    deallocate(ow_gradient_identity_map)
    if(size(ow_pencil_generator_maps,2)>0.and.ok)then
    call measure_ow_discrete_gradient_map_commutator(dc%icomm_tot,ow_core_values,&
      ow_direct_core_gradients,ow_core_weights,ow_core_ids,ow_pencil_generator_maps,&
      dc%lg_tot%num,stencil%coef_nab,global_point_rotations(:,:,global_affine_generators),&
      ow_gradient_map_commutator,ok,message)
    if(ok.and.rank==0)write(*,'(a,es16.8,a,i0)')&
      '[OW-GS-DIAGNOSTIC] finite-difference/map commutator max=',&
      maxval(ow_gradient_map_commutator),' operation=',maxloc(ow_gradient_map_commutator,dim=1)
    if(allocated(ow_gradient_map_commutator))deallocate(ow_gradient_map_commutator)
    if(allocated(ow_gradient_covariance_left))deallocate(ow_gradient_covariance_left)
    if(allocated(ow_gradient_covariance_transpose))deallocate(ow_gradient_covariance_transpose)
    if(allocated(ow_gradient_covariance_candidates))deallocate(ow_gradient_covariance_candidates)
    if(ok)call measure_dg_spatial_gradient_covariance(dc%icomm_tot,ow_direct_core_gradients,&
      ow_core_weights,ow_pencil_generator_maps,ow_pencil_generator_representation,&
      global_point_rotations(:,:,global_affine_generators),&
      ow_gradient_covariance_left,ow_gradient_covariance_transpose,ok,message)
    if(ok.and.rank==0)write(*,'(2(a,es16.8,a,i0))')&
      '[OW-GS-DIAGNOSTIC] direct core gradient covariance R max=',maxval(ow_gradient_covariance_left),&
      ' operation=',maxloc(ow_gradient_covariance_left,dim=1),&
      ' RT max=',maxval(ow_gradient_covariance_transpose),&
      ' operation=',maxloc(ow_gradient_covariance_transpose,dim=1)
    else if(size(ow_pencil_generator_maps,2)==0.and.ok)then
      if(rank==0)then
        write(*,'(a,es16.8,a,i0)')'[OW-GS-DIAGNOSTIC] finite-difference/map commutator max=',&
          0d0,' operation=',global_identity_operation
        write(*,'(2(a,es16.8,a,i0))')&
          '[OW-GS-DIAGNOSTIC] direct core gradient covariance R max=',0d0,&
          ' operation=',global_identity_operation,' RT max=',0d0,' operation=',global_identity_operation
      endif
    endif
    if(allocated(ow_gradient_covariance_left))deallocate(ow_gradient_covariance_left)
    if(allocated(ow_gradient_covariance_transpose))deallocate(ow_gradient_covariance_transpose)
    deallocate(ow_direct_core_gradients)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'pencil generator representation failed';endif
    allocate(ow_pencil_generator_operations,source=global_affine_generators)
    allocate(ow_pencil_affine_product,source=global_point_product)
    allocate(ow_pencil_translation_subgroup,source=global_translation_subgroup)
    allocate(ow_pencil_coset_representatives,source=global_point_representatives)
    allocate(all_core_ids(ncore,nproc),rank_fragments(nproc))
    call MPI_Allgather(ow_core_ids,ncore,MPI_INTEGER8,all_core_ids,ncore,MPI_INTEGER8,dc%icomm_tot,ierr)
    call MPI_Allgather(dc%i_frag,1,MPI_INTEGER,rank_fragments,1,MPI_INTEGER,dc%icomm_tot,ierr)
    call assign_dg_periodic_centers_to_fragments(dc%lg_tot%num,localized_centers,all_core_ids,&
      rank_fragments,retained_closure_search_tolerance,localized_center_ids,center_owner_candidate,&
      center_fragment_candidate,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'localized Wannier center ownership failed';end if
    allocate(center_box_candidate(ntarget))
    do io=1,ntarget
      if(center_owner_candidate(io)==rank)then
        p=findloc(ow_core_ids,localized_center_ids(io),dim=1)
        if(p<1)error stop 'localized center owner does not contain its core ID'
        center_box_candidate(io)=ow_core_box_positions(p)
      else
        center_box_candidate(io)=0
      end if
    end do
    call MPI_Allreduce(MPI_IN_PLACE,center_box_candidate,ntarget,MPI_INTEGER,MPI_SUM,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.any(center_box_candidate<1))error stop 'localized center box reduction failed'
    ow_basis%center_owner_rank=center_owner_candidate
    ow_basis%center_owner_fragment=center_fragment_candidate
    do io=1,ntarget
      ow_basis%center_box_point_ids(io)=int(center_owner_candidate(io),8)*int(nbox,8)+&
        int(center_box_candidate(io),8)
    end do
    ow_basis%generation=1
    allocate(ow_row_ids(count(ow_basis%center_owner_rank==rank)))
    io=0
    do p=1,ntarget
      if(ow_basis%center_owner_rank(p)/=rank)cycle
      io=io+1;ow_row_ids(io)=p
    enddo
    i=findloc(global_affine_generators,5,dim=1)
    if(i>0)then
      call diagnose_ow_total_nonlocal_projector_range(localized_centers,&
        global_point_integer_rotations(:,:,5),global_point_rotations(:,:,5),&
        global_point_fractional_translations(:,5),ow_pencil_generator_representation(:,:,i),&
        num_fragment,ok,message)
      if(.not.ok.and.rank==0)write(*,'(2a)')&
        '[OW-GS-DIAGNOSTIC] total nonlocal projector range unavailable: ',trim(message)
    else if(rank==0)then
      write(*,'(a)')'[OW-GS-DIAGNOSTIC] operation 5 is not an affine generator; projector range unavailable'
    endif
    call transpose_dg_spatial_cores_to_orbital_owners(dc%icomm_tot,ow_core_values,ow_core_ids,32,&
      orbital_owned_ids,orbital_owned_full_ids,orbital_owned_full_values,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'localized orbital ownership transpose failed';end if
    call redistribute_dg_owned_orbitals_to_center_fragments(dc%icomm_tot,orbital_owned_ids,&
      orbital_owned_full_ids,orbital_owned_full_values,center_owner_candidate,physical_ids,&
      center_local_orbital_ids,center_local_buffer_values,ok,message)
    if(.not.ok.or.any(center_local_orbital_ids/=&
        pack([(io,io=1,ntarget)],center_owner_candidate==rank)))then
      write(0,'(a)')trim(message);error stop 'localized center-fragment redistribution failed'
    end if
    deallocate(orbital_owned_ids,orbital_owned_full_ids,orbital_owned_full_values,&
      center_local_orbital_ids,center_local_buffer_values,all_core_ids,rank_fragments,localized_center_ids)
    if(rank==0)write(*,'(a,2(es12.4,1x))')'[OW-GS-DIAGNOSTIC] localized_center_moment_minmax=',&
      minval(localized_center_magnitudes),maxval(localized_center_magnitudes)
    deallocate(localized_centers,localized_center_magnitudes,center_owner_candidate,&
      center_box_candidate,center_fragment_candidate)
    call assemble_dg_stitched_overlap_density_rows(dc%icomm_tot,ntarget,ow_row_ids,physical_ids,&
      ow_partition_weight,ow_box_values,ow_box_density,system%hvol,expected_core_count,&
      dc%elec_num_tot,dg_dc_metric_rank_tolerance,dg_dc_gs_electron_count_tolerance,&
      ow_srows,ow_rhorows,ow_stitched_electron_count,&
      ow_stitched_s_hermiticity,ow_stitched_rho_hermiticity,ow_stitched_minimum_pivot,&
      ow_stitched_pivot_condition,ow_stitched_peak_elements,ok,message)
    if(rank==0)write(*,'(a,7(a,es24.16),a,i0)')'[OW-GS-DIAGNOSTIC] stitched_overlap_density',&
      ' electrons=',ow_stitched_electron_count,' s_hermiticity=',ow_stitched_s_hermiticity,&
      ' rho_hermiticity=',ow_stitched_rho_hermiticity,' minimum_pivot=',ow_stitched_minimum_pivot,&
      ' pivot_condition=',ow_stitched_pivot_condition,' expected_electrons=',dc%elec_num_tot,&
      ' electron_drift=',ow_stitched_electron_count-dc%elec_num_tot,&
      ' workspace_peak_elements=',ow_stitched_peak_elements
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'stitched overlap-density gate failed';endif
    condition_number=ow_stitched_pivot_condition
    allocate(spectrum(2));spectrum=[ow_stitched_minimum_pivot,&
      ow_stitched_minimum_pivot*ow_stitched_pivot_condition]
    allocate(ow_tail_generation(ntarget,ncore));ow_tail_generation=ow_basis%generation
    call compute_dg_overlapping_wannier_scf_fingerprint(dc%icomm_tot,ow_row_ids,&
      ow_srows,ow_core_ids,ow_core_weights,ow_core_values,ow_tail_generation,&
      basis_fingerprint,ow_symmetry_fingerprint)
    operator_fingerprint=ow_collective_operator_fingerprint(dc%icomm_tot)
    prefix='./overlapping_wannier_gs'
    call read_dg_overlapping_wannier_checkpoint(dc%icomm_tot,trim(prefix),ow_basis%generation,1,&
      basis_fingerprint,operator_fingerprint,[dg_dc_gs_final_density_tolerance,&
      dg_dc_gs_final_orbital_tolerance,10d0*dg_dc_gs_final_orbital_tolerance,&
      dg_dc_gs_electron_count_tolerance,1d0/dg_dc_metric_rank_tolerance,dg_ow_symmetry_tolerance],&
      ow_checkpoint,reusable,ok,message)
    if(ok.and.reusable.and.yn_dg_hybrid_scf/='y'.and.yn_dg_hybrid_continuation_scf/='y'.and.&
        yn_dg_hybrid_divided_scf/='y')then
      call restore_ow_checkpoint_density(ow_checkpoint,ok,message)
      if(.not.ok)error stop 'overlapping-Wannier checkpoint density restore failed'
      if(rank==0)write(*,'(a)')'[OW-GS] reused accepted route checkpoint'
      return
    endif
    if(rank==0)then
      write(*,'(3a)')'[OW-GS-DIAGNOSTIC] checkpoint read rejected: ',trim(message),''
      if(ok)write(*,'(a,4(z16.16,1x))')&
        '[OW-GS-DIAGNOSTIC] stored/expected basis/operator fingerprints=',&
        ow_checkpoint%basis_fingerprint,basis_fingerprint,&
        ow_checkpoint%operator_fingerprint,operator_fingerprint
    endif
    allocate(occupations(nstate))
    occupations=lcfo_retained_occupations(1:nstate)
    if(any(occupations<0d0).or.any(occupations>2d0).or.&
        abs(sum(occupations)-dc%elec_num_tot)>dg_dc_gs_electron_count_tolerance)&
      error stop 'retained LCFO occupation spectrum violates the electron-count gate'
    allocate(ow_initial_occupied_density(ncore))
    call redistribute_dg_row_owned_real_field_to_requests(dc%icomm_tot,expected_core_count,&
      ow_total_density_ids,ow_total_density_values,ow_core_ids,ow_initial_occupied_density,&
      ow_density_redistribution_workspace,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message)
      error stop 'converged DC+LCFO density redistribution failed';endif
    deallocate(ow_total_density_ids,ow_total_density_values)
    initial_occupied_charge_local=sum(ow_core_weights*ow_initial_occupied_density)
    call MPI_Allreduce(initial_occupied_charge_local,initial_occupied_charge,1,&
      MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
    if(rank==0)write(*,'(a,3(a,es24.16))')'[OW-GS-DIAGNOSTIC] converged DC+LCFO initial density',&
      ' electrons=',initial_occupied_charge,' expected=',dc%elec_num_tot,&
      ' difference=',initial_occupied_charge-dc%elec_num_tot
    if(ierr/=MPI_SUCCESS.or..not.all(ieee_is_finite(ow_initial_occupied_density)).or.&
        abs(initial_occupied_charge-dc%elec_num_tot)>&
        dg_dc_gs_electron_count_tolerance*max(1d0,dc%elec_num_tot))&
      error stop 'converged DC+LCFO initial density violates electron-count contract'
    if(yn_dg_hybrid_divided_scf=='y'.or.yn_dg_hybrid_continuation_scf=='y')then
      call prepare_dg_hybrid_divided_production_basis(dc%icomm_tot,dc%n_frag,dc%i_frag,dc%lg_tot%num,&
        dc%ixyz_frag,dc%nxyz_domain_frag,dc%system_tot%hgs,ncore,nbox,nxy8,expected_core_count,&
        physical_ids,ow_core_ids,ow_core_weights,&
        ow_core_values,ow_box_values,ow_basis%center_owner_fragment,divided_fragment_basis,&
        size(ow_basis%center_owner_fragment),size(ow_reciprocal_operation_maps,2),&
        ow_reciprocal_operation_maps,ow_raw_partition_weight,w90_reciprocal_lattice,&
        global_point_rotations(:,:,global_point_representatives),wannier_pw_cutoff,wannier_pw_max,&
        dg_ow_symmetry_tolerance,&
        basis_fingerprint,divided_pw_fingerprint,&
        divided_buffer_window_fingerprint,divided_fragment_fingerprint,divided_production_selection,ok,message)
      if(.not.ok)write(0,'(a)')trim(message)
      if(.not.ok)error stop 'divided Hybrid production basis preparation failed'
      if(rank==0)write(*,'(a,2(a,es24.16),2(a,i0))')'[HYBRID-PW-CUTOFF]',&
        ' requested=',divided_production_selection%requested_cutoff,&
        ' effective=',divided_production_selection%effective_cutoff,&
        ' shell_added=',divided_production_selection%shell_added,&
        ' orbit_added=',divided_production_selection%orbit_added
      if(basis_fingerprint==0_8.or.divided_pw_fingerprint==0_8)&
        error stop 'divided Hybrid LCFO symmetry handoff fingerprints are missing'
      if(rank==0)write(*,'(a,i0,a,i0)')&
        '[HYBRID-LCFO-SYMMETRY-HANDOFF] wannier_fingerprint=',basis_fingerprint,&
        ' production_fingerprint=',divided_pw_fingerprint
      allocate(divided_requested_ids,source=divided_production_selection%requested_packet_ids)
      call close_dg_hybrid_selection(dc%icomm_tot,divided_requested_ids,&
        divided_production_selection%packet_ids,divided_production_selection%packet_action,&
        divided_selection_effective_ids,divided_closure_parent,divided_closure_action,&
        divided_selection_fingerprint,ok,message)
      if(.not.ok)write(0,'(a)')trim(message)
      if(.not.ok)error stop 'divided Hybrid authoritative selection closure failed'
      if(allocated(ow_divided_core_mask))deallocate(ow_divided_core_mask)
      allocate(ow_divided_core_mask,source=core_mask)
      call prepare_dg_hybrid_divided_dc_controls(dc,ow_hybrid_divided_convergence,&
        ow_hybrid_divided_threshold,ow_hybrid_divided_total_density,ok,message)
      if(.not.ok)write(0,'(a)')trim(message)
      if(.not.ok)error stop 'divided Hybrid DC control preparation failed'
      allocate(divided_initial_density(ncore))
      do p=1,ncore
        ix=int(modulo(ow_core_ids(p)-1_8,int(dc%lg_tot%num(1),8)))+1
        iy=int(modulo((ow_core_ids(p)-1_8)/int(dc%lg_tot%num(1),8),int(dc%lg_tot%num(2),8)))+1
        iz=int((ow_core_ids(p)-1_8)/nxy8)+1
        divided_initial_density(p)=ow_hybrid_divided_total_density(ix,iy,iz)
      enddo
      allocate(divided_lcfo_row_ids,source=divided_fragment_basis%global_ids)
        allocate(divided_fragment_bases(dc%n_frag))
        do p=1,dc%n_frag
          allocate(divided_fragment_bases(p)%global_ids(0),divided_fragment_bases(p)%sector(0),&
            divided_fragment_bases(p)%buffer_point_ids(0),divided_fragment_bases(p)%buffer_values(0,0))
        enddo
        divided_fragment_bases(dc%i_frag)=divided_fragment_basis
        divided_global_basis_count=0
        if(size(divided_fragment_basis%global_ids)>0)&
          divided_global_basis_count=int(maxval(divided_fragment_basis%global_ids))
        call MPI_Allreduce(MPI_IN_PLACE,divided_global_basis_count,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
        if(ierr/=MPI_SUCCESS.or.divided_global_basis_count<1)&
          error stop 'divided Hybrid global basis count failed'
        allocate(divided_effective_ids(divided_global_basis_count))
        divided_effective_ids=[(p,p=1,size(divided_effective_ids))]
        call freeze_dg_hybrid_basis_directory(dc%icomm_tot,divided_fragment_bases,divided_effective_ids,&
          divided_basis_owner,divided_basis_fragment,ok,message)
        if(.not.ok)write(0,'(a)')trim(message)
        if(.not.ok)error stop 'divided Hybrid fixed-basis directory construction failed'
        allocate(dg_hybrid_interior_fragment(size(ow_core_ids)),&
          dg_hybrid_interior_weights(size(ow_core_ids)),dg_hybrid_unit_local_potential(size(ow_core_ids)))
        dg_hybrid_interior_fragment=dc%i_frag
        dg_hybrid_interior_weights=system%hvol
        dg_hybrid_unit_local_potential=1d0
        call materialize_dg_hybrid_production_interior(dc%icomm_tot,dc%lg_tot%num,stencil%coef_nab,&
          stencil%coef_lap0,stencil%coef_lap,&
          divided_fragment_bases,divided_basis_owner,divided_basis_fragment,divided_effective_ids,&
          ow_core_ids,dg_hybrid_interior_fragment,dg_hybrid_interior_values,dg_hybrid_interior_gradients,&
          dg_hybrid_interior_kinetic_action,ok,message)
        if(.not.ok)write(0,'(a)')trim(message)
        if(.not.ok)error stop 'divided Hybrid production interior materialization failed'
        call assemble_dg_hybrid_broken_volume_rows(dc%icomm_tot,size(divided_effective_ids),&
          divided_lcfo_row_ids,divided_basis_fragment,ow_core_ids,dg_hybrid_interior_fragment,&
          dg_hybrid_interior_weights,dg_hybrid_interior_values,dg_hybrid_interior_gradients,&
          dg_hybrid_unit_local_potential,dg_hybrid_kinetic_rows,divided_lcfo_srows,&
          dg_hybrid_broken_diagnostics,ok,message)
        if(.not.ok)write(0,'(a)')trim(message)
        if(.not.ok)error stop 'divided Hybrid broken-volume kinetic assembly failed'
        call ow_fingerprint_distributed_matrix(dc%icomm_tot,divided_lcfo_row_ids,&
          divided_lcfo_srows,divided_lcfo_operator_fingerprint,ok)
        if(.not.ok)error stop 'divided Hybrid variational metric fingerprint failed'
        allocate(dg_hybrid_nonlocal_rows(size(divided_lcfo_row_ids),size(divided_effective_ids)))
        call assemble_dg_hybrid_divided_nonlocal_rows(divided_fragment_basis,divided_lcfo_row_ids,&
          size(divided_effective_ids),dg_hybrid_nonlocal_rows,dg_hybrid_interior_nonlocal_action,&
          dg_hybrid_nonlocal_ownership_count,ok,message)
        if(.not.ok)write(0,'(a)')trim(message)
        if(.not.ok)error stop 'divided Hybrid complete nonlocal assembly failed'
        call materialize_dg_hybrid_production_face_collection(dc%icomm_tot,dc%ixyz_frag,&
          dc%nxyz_domain_frag,dc%lg_tot%num,dc%system_tot%hgs,stencil%coef_nab,&
          divided_fragment_bases,divided_basis_owner,divided_basis_fragment,divided_effective_ids,&
          divided_production_faces,ok,message)
        if(.not.ok)write(0,'(a)')trim(message)
        if(.not.ok)error stop 'divided Hybrid production face materialization failed'
        call assemble_dg_hybrid_production_interface_component_rows(dc%icomm_tot,size(divided_effective_ids),&
          divided_lcfo_row_ids,divided_production_faces,dg_dc_gs_sipg_penalty_factor,&
          dg_hybrid_interface_component_rows,ok,message)
        if(.not.ok)write(0,'(a)')trim(message)
        if(.not.ok)error stop 'divided Hybrid SIPG interface assembly failed'
        allocate(dg_hybrid_interface_rows(size(divided_lcfo_row_ids),size(divided_effective_ids)))
        dg_hybrid_interface_rows=sum(dg_hybrid_interface_component_rows,dim=3)
        variational_payload_capture_prefix=''
        call get_environment_variable('SALMON_DG_VARIATIONAL_PAYLOAD_CAPTURE',&
          variational_payload_capture_prefix,length=variational_payload_capture_length,&
          status=variational_payload_capture_status,trim_name=.true.)
        if(variational_payload_capture_status==0.and.variational_payload_capture_length>0)then
          call write_dg_hybrid_variational_payload_bundle(dc%icomm_tot,&
            trim(variational_payload_capture_prefix),size(divided_effective_ids),divided_lcfo_row_ids,&
            divided_lcfo_srows,dg_hybrid_kinetic_rows,dg_hybrid_nonlocal_rows,dg_hybrid_interface_rows,&
            divided_fragment_fingerprint,divided_lcfo_operator_fingerprint,&
            divided_buffer_window_fingerprint,ok,message)
          if(.not.ok)write(0,'(a)')trim(message)
          if(.not.ok)error stop 'divided Hybrid variational payload capture failed'
          if(nproc_id_global==0)write(*,'(a,a)')&
            '[HYBRID-VARIATIONAL-PAYLOAD-CAPTURE] prefix=',trim(variational_payload_capture_prefix)
        endif
        if(yn_dg_hybrid_divided_scf=='y')then
          call freeze_dg_hybrid_single_owner_payload(dc%icomm_tot,dc%n_frag,divided_fragment_basis,&
            divided_lcfo_srows,dg_hybrid_kinetic_rows,dg_hybrid_nonlocal_rows,dg_hybrid_interface_rows,&
            divided_fragment_fingerprint,divided_lcfo_operator_fingerprint,divided_buffer_window_fingerprint,&
            dg_hybrid_fixed_payload,divided_basis_owner,divided_basis_fragment,divided_basis_local_slot,&
            divided_basis_generation,divided_basis_directory_fingerprint,ok,message)
        else
          call freeze_dg_hybrid_variational_payload(dc%icomm_tot,size(divided_effective_ids),&
            divided_lcfo_row_ids,divided_lcfo_srows,dg_hybrid_kinetic_rows,dg_hybrid_nonlocal_rows,&
            dg_hybrid_interface_rows,divided_fragment_fingerprint,divided_lcfo_operator_fingerprint,&
            divided_buffer_window_fingerprint,dg_hybrid_fixed_payload,ok,message)
        endif
        if(.not.ok)write(0,'(a)')trim(message)
        if(.not.ok)error stop 'divided Hybrid fixed variational payload freeze failed'
        if(yn_dg_hybrid_divided_scf=='y')then
          divided_mixing_basis_generation=maxval(divided_basis_generation)
          divided_mixing_inventory_fingerprint=divided_basis_directory_fingerprint
        endif
        if(.not.dg_hybrid_fixed_payload%frozen.or.dg_hybrid_fixed_payload%fingerprint==0_8)&
          error stop 'divided Hybrid shared variational payload fingerprint is invalid'
        divided_fixed_payload_fingerprint=dg_hybrid_fixed_payload%fingerprint
        if(rank==0)write(*,'(a,i0)')'[HYBRID-SHARED-VARIATIONAL-PAYLOAD] fingerprint=',&
          divided_fixed_payload_fingerprint
      if(yn_dg_hybrid_continuation_scf=='y')then
        if(divided_fixed_payload_fingerprint==0_8.or.&
          divided_fixed_payload_fingerprint/=dg_hybrid_fixed_payload%fingerprint)&
          error stop 'DG continuation shared variational payload fingerprint changed'
        call build_dg_hybrid_retained_basis_representation(dc%icomm_tot,int(expected_core_count),&
          ow_core_ids,ow_core_weights,ow_pencil_generator_maps,divided_fragment_basis,&
          divided_lcfo_row_ids,divided_lcfo_srows,dg_ow_symmetry_tolerance,&
          divided_basis_representation,divided_full_basis_closure_defect,&
          divided_full_basis_closure_ok,ok,message)
        if(.not.ok)write(0,'(a)')trim(message)
        if(.not.ok)error stop 'DG continuation retained-basis symmetry representation failed'
        if(rank==0)write(*,'(a,i0,a,es24.16)')&
          '[HYBRID-RETAINED-BASIS-SYMMETRY] closed=',merge(1,0,divided_full_basis_closure_ok),&
          ' defect=',divided_full_basis_closure_defect
        call selection_added_members(divided_requested_ids,divided_selection_effective_ids,divided_added_ids)
        allocate(divided_closure_reason(size(divided_closure_parent)),source=1)
        call build_dg_hybrid_scope_receipt(dc%icomm_tot,merge(1,0,theory=='dft'),iperiodic==3,&
          dc%system_tot%nspin,yn_spinorbit=='y',PLUS_U_ON,yn_hse=='y',yn_fix_func=='y',yn_jm=='y',&
          xc_func%xctype,divided_scope_receipt,ok,message)
        if(.not.ok)write(0,'(a)')trim(message)
        if(.not.ok)error stop 'DG continuation supported-scope receipt failed'
        allocate(divided_scope_selectors(8));divided_scope_selectors=[divided_scope_receipt%theory_code,&
          merge(1,0,divided_scope_receipt%periodic),divided_scope_receipt%nspin,&
          merge(1,0,divided_scope_receipt%spinorbit),merge(1,0,divided_scope_receipt%plus_u),&
          merge(1,0,divided_scope_receipt%hse),merge(1,0,divided_scope_receipt%fix_func),&
          merge(1,0,divided_scope_receipt%jm)]
        call build_checkpoint_topology_graphs(divided_effective_ids,divided_lcfo_row_ids,divided_basis_fragment,&
          divided_production_faces,divided_metric_offsets,divided_metric_columns,&
          divided_operator_offsets,divided_operator_columns)
        call run_dg_hybrid_concrete_continuation(divided_initial_density,divided_effective_ids,&
          divided_lcfo_row_ids,divided_basis_fragment,dg_hybrid_interior_fragment,&
          dg_hybrid_interior_weights,dg_hybrid_interior_values,dg_hybrid_interior_gradients,&
          dg_hybrid_interior_kinetic_action,dg_hybrid_interior_nonlocal_action,dg_hybrid_fixed_payload,&
          dg_hybrid_interface_component_rows,divided_production_faces,ntarget,&
          ow_pencil_generator_maps,divided_basis_representation,&
          global_point_rotations(:,:,global_affine_generators),&
          divided_requested_ids,divided_selection_effective_ids,divided_added_ids,&
          divided_closure_parent,divided_closure_reason,&
          divided_closure_action,divided_scope_selectors,pseudopotential_fingerprint,&
          divided_scope_receipt%fingerprint,divided_selection_fingerprint,divided_metric_offsets,divided_metric_columns,&
          divided_operator_offsets,divided_operator_columns,&
          ow_hybrid_ground_state,dg_hybrid_final_density,dg_hybrid_final_trace,&
          dg_hybrid_final_hamiltonian_rows)
        return
      else
        if(divided_fixed_payload_fingerprint==0_8.or.&
          divided_fixed_payload_fingerprint/=dg_hybrid_fixed_payload%fingerprint)&
          error stop 'divided Hybrid shared variational payload fingerprint changed'
        call run_dg_hybrid_divided_scf(dc%icomm_tot,int(expected_core_count),ow_core_ids,&
          divided_initial_density,dc%system_tot%hvol,&
          ow_hybrid_divided_convergence,ow_hybrid_divided_threshold,&
          update_dg_hybrid_divided_potential,solve_dg_hybrid_divided_fragments,&
          assemble_dg_hybrid_divided_core_density,mix_dg_hybrid_divided_density,nscf,&
          ow_core_weights,dc%elec_num_tot,dg_dc_gs_electron_count_tolerance,&
          divided_converged_density,divided_iterations,divided_convergence_value,&
          divided_electron_defect,ok,message)
        if(.not.ok)write(0,'(a)')trim(message)
        if(.not.ok)error stop 'divided Hybrid SCF failed'
        if(rank==0)write(*,'(a,i0,a,es16.8,a,es16.8)')&
          '[OW-GS] divided WF+PW SCF converged iterations=',divided_iterations,&
          ' density=',divided_convergence_value,' electron_defect=',divided_electron_defect
        allocate(divided_lcfo_point_weights(size(divided_fragment_basis%buffer_point_ids)),source=system%hvol)
        call assemble_dg_hybrid_lcfo_rows(dc%icomm_tot,divided_fragment_basis,divided_lcfo_row_ids,&
          divided_lcfo_point_weights,apply_dg_hybrid_divided_fragment_hpsi,&
          apply_dg_hybrid_divided_fragment_metric,divided_lcfo_hrows,divided_lcfo_srows,&
          divided_lcfo_peak_elements,divided_lcfo_operator_fingerprint,ok,message)
        if(.not.ok)write(0,'(a)')trim(message)
        if(.not.ok)error stop 'divided Hybrid LCFO assembly failed'
      endif
      call build_dg_hybrid_retained_basis_representation(dc%icomm_tot,int(expected_core_count),&
        ow_core_ids,ow_core_weights,ow_pencil_generator_maps,divided_fragment_basis,&
        divided_lcfo_row_ids,divided_lcfo_srows,dg_ow_symmetry_tolerance,&
        divided_basis_representation,divided_full_basis_closure_defect,&
        divided_full_basis_closure_ok,ok,message)
      if(.not.ok)write(0,'(a)')trim(message)
      if(.not.ok)error stop 'divided Hybrid retained-basis symmetry representation failed'
      if(rank==0)write(*,'(a,i0,a,es24.16)')&
        '[HYBRID-RETAINED-BASIS-SYMMETRY] closed=',merge(1,0,divided_full_basis_closure_ok),&
        ' defect=',divided_full_basis_closure_defect
      call solve_dg_hybrid_generalized_once_and_publish(dc%icomm_tot,size(divided_lcfo_hrows,2),nstate,&
        divided_lcfo_row_ids,divided_lcfo_hrows,divided_lcfo_srows,dg_dc_gs_final_orbital_tolerance,&
        occupations,dc%elec_num_tot,divided_fragment_fingerprint,divided_lcfo_operator_fingerprint,&
        divided_lcfo_operator_fingerprint,divided_pw_fingerprint,solve_final_dg_hybrid_divided_lcfo,&
        ow_hybrid_ground_state,divided_state_workspace,divided_state_fingerprint,divided_final_residual,&
        divided_final_orthogonality,divided_final_projector_defect,divided_final_solver_workspace,&
        divided_final_solver_fingerprint,ok,message)
      if(.not.ok)write(0,'(a)')trim(message)
      if(.not.ok)error stop 'divided Hybrid final LCFO solve failed'
      if(rank==0)write(*,'(a,3(a,es16.8))')'[OW-GS] divided WF+PW LCFO solved once',&
        ' residual=',divided_final_residual,' orthogonality=',divided_final_orthogonality,&
        ' projector=',divided_final_projector_defect
      hybrid_provenance=[divided_pw_fingerprint,divided_buffer_window_fingerprint,&
        divided_fragment_fingerprint,divided_solver_fingerprint,divided_lcfo_operator_fingerprint,&
        divided_final_solver_fingerprint]
      hybrid_scf_receipts=[divided_convergence_value,divided_final_residual,divided_final_orthogonality,&
        divided_final_projector_defect,divided_electron_defect]
      call write_rt_dg_hybrid_occupied_checkpoint(dc%icomm_tot,'./overlapping_wannier_occupied.chk',&
        ow_hybrid_ground_state%global_count,ow_hybrid_ground_state%owned_row_ids,&
        ow_hybrid_ground_state%coefficients,ow_hybrid_ground_state%occupations,&
        ow_hybrid_ground_state%eigenvalues,divided_pw_fingerprint,divided_fragment_fingerprint,&
        hybrid_provenance,divided_lcfo_operator_fingerprint,divided_state_fingerprint,&
        hybrid_scf_receipts,max(dg_dc_gs_final_orbital_tolerance,dg_dc_gs_electron_count_tolerance,&
        dg_ow_symmetry_tolerance,maxval(hybrid_scf_receipts)),&
        hybrid_checkpoint_fingerprint,ok,message)
      if(.not.ok)write(0,'(a)')trim(message)
      if(.not.ok)error stop 'divided Hybrid occupied checkpoint failed'
      return
    endif
    if(yn_dg_hybrid_scf=='y')then
      if(rank==0)write(*,'(a)')'[OW-GS] starting distributed fixed-basis complex ScaLAPACK SCF'
      if(dc%system_tot%nspin/=1)error stop 'distributed hybrid SCF currently requires nspin=1'
      call ow_fingerprint_distributed_matrix(dc%icomm_tot,ow_row_ids,ow_srows,&
        ow_hybrid_metric_fingerprint,ok)
      if(.not.ok)error stop 'distributed hybrid metric fingerprint failed'
      allocate(ow_hybrid_hrows(size(ow_row_ids),ntarget),ow_hybrid_occupations(nstate),&
        ow_hybrid_eigenvalues(nstate),ow_hybrid_potential(ncore),ow_hybrid_density(ncore),&
        ow_hybrid_density_history(ncore,2),ow_hybrid_new_history(ncore,2))
      ow_hybrid_occupations=occupations;ow_hybrid_eigenvalues=0d0;ow_hybrid_potential=0d0
      ow_hybrid_density=ow_initial_occupied_density
      ow_hybrid_density_history(:,1)=ow_hybrid_density
      ow_hybrid_density_history(:,2)=ow_hybrid_density
      ow_hybrid_new_history=ow_hybrid_density_history;ow_hybrid_history_count=0
      ow_hybrid_mixing_rate=dg_dc_gs_density_mix_rate
      call run_dg_hybrid_self_consistent_ground_state(dc%icomm_tot,int(expected_core_count),ow_core_ids,&
        ow_hybrid_density,basis_fingerprint,ow_hybrid_metric_fingerprint,ow_hybrid_update_potential,&
        ow_hybrid_assemble_hamiltonian,ow_hybrid_solve_occupied,ow_hybrid_reconstruct_density,&
        ow_hybrid_density_mix,dg_dc_gs_maximum_scf_iterations,dg_dc_gs_final_density_tolerance,&
        dg_dc_gs_final_orbital_tolerance,dg_dc_gs_final_orbital_tolerance,&
        min(dg_dc_gs_electron_count_tolerance,dg_ow_symmetry_tolerance),hybrid_converged_density,&
        hybrid_iterations,hybrid_density_residual,hybrid_energy_residual,hybrid_eigensystem_residual,&
        hybrid_electron_defect,hybrid_symmetry_defect,hybrid_scf_fingerprint,ok,message)
      if(.not.ok)then;write(0,'(a)')trim(message);error stop 'distributed hybrid SCF failed';endif
      call validate_dg_hybrid_ground_state(dc%icomm_tot,ntarget,nstate,ow_row_ids,ow_hybrid_coefficients,&
        occupations,ow_hybrid_eigenvalues,dc%elec_num_tot,basis_fingerprint,ow_hybrid_metric_fingerprint,&
        ow_hybrid_operator_fingerprint,ow_symmetry_fingerprint,dg_dc_gs_electron_count_tolerance,&
        ow_hybrid_ground_state,hybrid_state_workspace,hybrid_state_fingerprint,ok,message)
      if(.not.ok)then;write(0,'(a)')trim(message);error stop 'distributed hybrid state validation failed';endif
      ow_hybrid_ground_state%converged=.true.
      hybrid_provenance=[ow_hybrid_metric_fingerprint,spectral_basin_fingerprint,&
        spectral_action_aggregate_fingerprint,translation_post_gauge_fingerprint,&
        w90_input_fingerprint,w90_transform_fingerprint]
      hybrid_scf_receipts=[hybrid_density_residual,hybrid_energy_residual,hybrid_eigensystem_residual,&
        hybrid_electron_defect,hybrid_symmetry_defect]
      call write_rt_dg_hybrid_occupied_checkpoint(dc%icomm_tot,'./overlapping_wannier_occupied.chk',&
        ntarget,ow_row_ids,ow_hybrid_coefficients,occupations,ow_hybrid_eigenvalues,spectral_catalog_fingerprint,&
        basis_fingerprint,hybrid_provenance,ow_hybrid_operator_fingerprint,hybrid_state_fingerprint,&
        hybrid_scf_receipts,max(dg_dc_gs_final_density_tolerance,dg_dc_gs_final_orbital_tolerance,&
        dg_dc_gs_electron_count_tolerance,dg_ow_symmetry_tolerance),hybrid_checkpoint_fingerprint,ok,message)
      if(.not.ok)then;write(0,'(a)')trim(message);error stop 'distributed hybrid checkpoint failed';endif
      if(rank==0)write(*,'(a,i0,5(a,es16.8))')'[OW-GS] hybrid SCF converged iterations=',hybrid_iterations,&
        ' density=',hybrid_density_residual,' band_energy_change=',hybrid_energy_residual,&
        ' eigensystem=',hybrid_eigensystem_residual,' electrons=',hybrid_electron_defect,&
        ' symmetry=',hybrid_symmetry_defect
      return
    endif
    allocate(ow_state%density(ncore),ow_state%potential(ncore),ow_state%coefficients(ntarget,nstate),&
      ow_state%eigenvalues(nstate),ow_state%density_history(ncore,2))
    ow_state%coefficients=(0d0,0d0);do io=1,nstate;ow_state%coefficients(io,io)=1d0;enddo
    allocate(one_shot_density(ncore))
    one_shot_density=ow_initial_occupied_density
    ow_state%density=one_shot_density;ow_state%potential=0d0
    ow_state%eigenvalues=0d0;ow_state%density_history(:,1)=ow_state%density
    ow_state%density_history(:,2)=ow_state%density;ow_state%history_count=1
    ow_state%basis_generation=ow_basis%generation;ow_state%geometry_generation=1
    ow_state%basis_fingerprint=basis_fingerprint
    ow_state%operator_fingerprint=operator_fingerprint
    allocate(one_shot_hrows(size(ow_row_ids),ntarget),one_shot_potential(ncore))
    ow_scalar_probe(1,:)=cmplx(ow_state%density,0d0,8)
    call measure_dg_spatial_basis_covariance(dc%icomm_tot,ow_scalar_probe,ow_scalar_probe_weights,&
      ow_pencil_generator_maps,ow_scalar_representation,ow_scalar_probe_residual,ok,message)
    if(ok.and.rank==0)write(*,'(a,es16.8,a,i0)')&
      '[OW-GS-DIAGNOSTIC] core input-density covariance max=',maxval(ow_scalar_probe_residual),&
      ' operation=',maxloc(ow_scalar_probe_residual,dim=1)
    if(allocated(ow_scalar_probe_residual))deallocate(ow_scalar_probe_residual)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'one-shot input density covariance failed';endif
    call ow_build_hamiltonian(dc%icomm_tot,ow_state%density,one_shot_hrows,&
      one_shot_potential,one_shot_operator_fingerprint,ok,message)
    ow_scalar_probe(1,:)=cmplx(one_shot_potential,0d0,8)
    call measure_dg_spatial_basis_covariance(dc%icomm_tot,ow_scalar_probe,ow_scalar_probe_weights,&
      ow_pencil_generator_maps,ow_scalar_representation,ow_scalar_probe_residual,diagnostic_ok,diagnostic_message)
    if(diagnostic_ok.and.rank==0)write(*,'(a,es16.8,a,i0)')&
      '[OW-GS-DIAGNOSTIC] core updated-potential covariance max=',maxval(ow_scalar_probe_residual),&
      ' operation=',maxloc(ow_scalar_probe_residual,dim=1)
    if(allocated(ow_scalar_probe_residual))deallocate(ow_scalar_probe_residual)
    deallocate(ow_scalar_probe,ow_vector_probe,ow_scalar_representation,ow_scalar_probe_weights)
    if(allocated(ow_pencil_generator_maps))deallocate(ow_pencil_generator_maps)
    if(allocated(ow_reciprocal_operation_maps))deallocate(ow_reciprocal_operation_maps)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'one-shot stitched Hamiltonian build failed';endif
    if(one_shot_operator_fingerprint/=operator_fingerprint)&
      error stop 'one-shot stitched Hamiltonian operator fingerprint mismatch'
#ifdef USE_EIGENEXA
    ow_saved_eigenexa_comm=info%icomm_o
    call finalize_eigenexa(info);info%icomm_o=dc%icomm_tot
    call init_eigenexa_mod(info,ntarget,direct_block_only=.true.)
    call solve_dg_overlapping_wannier_generalized_eigenexa(info,dc%icomm_tot,ow_row_ids,&
      one_shot_hrows,ow_srows,nstate,dg_dc_gs_final_orbital_tolerance,dg_dc_metric_rank_tolerance,&
      dg_dc_gs_final_orbital_tolerance,ow_state%coefficients,ow_state%eigenvalues,&
      one_shot_residual,one_shot_orthogonality,one_shot_condition,one_shot_gamma_defect,&
      one_shot_workspace_peak,ok,message)
    call finalize_eigenexa(info);info%icomm_o=ow_saved_eigenexa_comm
    call init_eigenexa_mod(info,system%no)
#else
    ok=.false.;message='one-shot overlapping-Wannier generalized solve requires EigenExa'
#endif
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'one-shot generalized EigenExa solve failed';endif
    call reconstruct_dg_overlapping_wannier_density(dc%icomm_tot,ow_core_ids,ow_core_weights,&
      ow_core_values,ow_tail_generation,ow_basis%generation,ow_state%coefficients,occupations,&
      expected_core_count,ow_row_ids,ow_srows,dg_dc_gs_final_density_tolerance,one_shot_density,&
      one_shot_charge,one_shot_trace_charge,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'one-shot occupied density reconstruction failed';endif
    one_shot_local_difference=sum((one_shot_density-ow_state%density)**2)
    one_shot_local_norm=sum(ow_state%density**2)
    call MPI_Allreduce(one_shot_local_difference,one_shot_global_difference,1,MPI_DOUBLE_PRECISION,&
      MPI_SUM,dc%icomm_tot,ierr)
    call MPI_Allreduce(one_shot_local_norm,one_shot_global_norm,1,MPI_DOUBLE_PRECISION,&
      MPI_SUM,dc%icomm_tot,ierr)
    ow_state%density=one_shot_density;ow_state%potential=one_shot_potential
    ow_state%operator_fingerprint=operator_fingerprint;ow_state%accepted=.true.
    ow_result=s_dg_overlapping_wannier_scf_result();ow_result%converged=.true.
    ow_result%iterations=1;ow_result%hamiltonian_rebuilds=1
    ow_result%density_residual=sqrt(one_shot_global_difference/max(tiny(1d0),one_shot_global_norm))
    ow_result%unmixed_density_residual=ow_result%density_residual
    ow_result%coefficient_residual=one_shot_residual
    ow_result%orthogonality_defect=one_shot_orthogonality
    ow_result%integrated_charge=one_shot_charge;ow_result%trace_charge=one_shot_trace_charge
    ow_result%symmetry_closure_residual=global_retained_group_closure_defect
    condition_number=one_shot_condition
    if(allocated(ow_published_hrows))deallocate(ow_published_hrows)
    allocate(ow_published_hrows,source=one_shot_hrows)
    if(rank==0)write(*,'(a,5(a,es16.8),a,i0)')'[OW-GS-DIAGNOSTIC] one_shot_generalized_eigenexa',&
      ' density_change=',ow_result%density_residual,' residual=',one_shot_residual,&
      ' s_orthogonality=',one_shot_orthogonality,' metric_condition=',one_shot_condition,&
      ' gamma_real_defect=',one_shot_gamma_defect,' workspace_peak_bytes=',one_shot_workspace_peak
    call populate_ow_checkpoint(occupations,condition_number,global_retained_group_closure_defect,&
      operator_fingerprint,&
      localization_initial_spread,localization_final_spread,&
      localization_maximum_gradient,localization_iterations,localization_converged,&
      w90_input_fingerprint,w90_transform_fingerprint,w90_spread,w90_coordinator_bytes,&
      w90_workspace_peak,w90_byte_limit,[w90_identity_defect,w90_unitarity_defect,w90_closure_defect],&
      .true.,w90_symmetry_workspace_peak,global_point_integer_rotations,&
      global_point_fractional_translations,adapted_occupied_subspace_distance,&
      abs(adapted_occupied_trace-real(nstate,8)),adapted_occupied_density_interior_difference,&
      adapted_occupied_density_boundary_difference,adapted_occupied_closure_before,&
      adapted_occupied_closure,adapted_occupied_selected_edge,adapted_occupied_rejected_edge,&
      adapted_occupied_cluster_gap,adapted_occupied_selected_block_dimension,&
      adapted_occupied_workspace_peak)
    call write_dg_overlapping_wannier_checkpoint(dc%icomm_tot,trim(prefix),ow_checkpoint,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'overlapping-Wannier checkpoint publication failed';endif
    call compute_ow_periodic_spread(dc%icomm_tot,spread_max)
    call write_ow_ground_state_evidence(spectrum,noccupied,nproc,rank,spread_max,&
      local_target_count,complete_sp_core_atom_count)
    deallocate(ow_box_gradients)
  end subroutine

  subroutine build_ow_complete_sp_projectors(physical_ids,pseudopotential_fingerprint,&
      channels,values,ok,message,requested_channels)
    integer(8),intent(in)::physical_ids(:),pseudopotential_fingerprint
    type(t_dg_projection_channel),allocatable,intent(out)::channels(:)
    real(8),allocatable,intent(out)::values(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(t_dg_projection_channel),intent(in),optional::requested_channels(:)
    integer,allocatable::core_atom_ids(:),radial_count(:,:),atomic_orbital_ordinals(:,:)
    real(8),allocatable::positions(:,:),radial_grid(:,:,:),radial_projector(:,:,:)
    real(8)::lattice_inverse(3,3),fractional(3),determinant
    integer::atom,axis,grid_index(3),relative_index(3),core_atom_count,species,ll,&
      point,nxy,allocation_status

    ok=.false.;message=''
    if(size(physical_ids)<1.or.dc%system_tot%nion<1.or.pseudopotential_fingerprint==0_8)then
      message='invalid production complete-s+p projector request';return
    endif
    if(.not.allocated(pp%nrps_ao).or..not.allocated(pp%mlps).or..not.allocated(pp%nproj).or.&
        .not.allocated(pp%rad).or..not.allocated(pp%upptbl_ao))then
      message='complete-s+p projector pseudopotential tables are unavailable';return
    endif
    call invert_ow_lattice(dc%system_tot%primitive_a,lattice_inverse,determinant,ok)
    if(.not.ok)then;message='complete-s+p projector lattice is singular';return;endif

    if(present(requested_channels))then
      if(size(requested_channels)<1)then;message='empty complete-s+p requested channel tile';return;endif
      allocate(channels,source=requested_channels)
    else
      allocate(core_atom_ids(dc%system_tot%nion));core_atom_count=0
      do atom=1,dc%system_tot%nion
        fractional=modulo(matmul(lattice_inverse,dc%system_tot%Rion(:,atom)),1d0)
        grid_index=modulo(floor(fractional*real(dc%lg_tot%num,8)),dc%lg_tot%num)
        if(dg_periodic_grid_point_owned(grid_index,dc%ixyz_frag(:,dc%i_frag),&
            dc%nxyz_domain_frag(:,dc%i_frag),dc%lg_tot%num))then
          core_atom_count=core_atom_count+1;core_atom_ids(core_atom_count)=atom
        endif
      enddo
      if(core_atom_count<1)then;message='complete-s+p fragment owns no core atom';return;endif
      call build_dg_complete_sp_manifest(core_atom_ids(1:core_atom_count),channels,ok,message)
      if(.not.ok)return
    endif

    allocate(positions(3,size(physical_ids)),radial_grid(pp%nrmax,2,dc%system_tot%nion),&
      radial_projector(pp%nrmax,2,dc%system_tot%nion),radial_count(2,dc%system_tot%nion),&
      atomic_orbital_ordinals(2,size(pp%mlps)),&
      stat=allocation_status)
    if(allocation_status/=0)then
      message='cannot allocate complete-s+p pseudo-atomic orbital workspace';return
    endif
    radial_grid=0d0;radial_projector=0d0;radial_count=0
    call select_dg_sp_atomic_orbital_ordinals(pp%mlps,pp%nproj,&
      atomic_orbital_ordinals,ok,message)
    if(.not.ok)return
    nxy=dc%lg_tot%num(1)*dc%lg_tot%num(2)
    do point=1,size(physical_ids)
      grid_index(1)=int(modulo(physical_ids(point)-1_8,int(dc%lg_tot%num(1),8)))
      grid_index(2)=int(modulo((physical_ids(point)-1_8)/int(dc%lg_tot%num(1),8),&
        int(dc%lg_tot%num(2),8)))
      grid_index(3)=int((physical_ids(point)-1_8)/int(nxy,8))
      positions(:,point)=matmul(dc%system_tot%primitive_a,&
        real(grid_index,8)/real(dc%lg_tot%num,8))
    enddo
    do atom=1,dc%system_tot%nion
      species=dc%system_tot%kion(atom)
      if(species<1.or.species>size(pp%mlps).or.species>size(pp%rad,2).or.&
          species>size(pp%nrps_ao).or.species>size(pp%upptbl_ao,3))then
        message='complete-s+p atom species is outside pseudopotential tables';return
      endif
      if(pp%mlps(species)<1.or.pp%nrps_ao(species)<2.or.pp%nrps_ao(species)>pp%nrmax.or.&
          pp%nrps_ao(species)>size(pp%rad,1).or.pp%nrps_ao(species)>size(pp%upptbl_ao,1).or.&
          lbound(pp%upptbl_ao,2)>0.or.ubound(pp%upptbl_ao,2)<1)then
        message='complete-s+p pseudo-atomic orbital table is incomplete';return
      endif
      do ll=0,1
        if(atomic_orbital_ordinals(ll+1,species)<lbound(pp%upptbl_ao,2).or.&
            atomic_orbital_ordinals(ll+1,species)>ubound(pp%upptbl_ao,2))then
          message='complete-s+p pseudo-atomic orbital ordinal is invalid';return
        endif
        radial_count(ll+1,atom)=pp%nrps_ao(species)
        radial_grid(1:pp%nrps_ao(species),ll+1,atom)=&
          pp%rad(1:pp%nrps_ao(species),species)
        radial_projector(1:pp%nrps_ao(species),ll+1,atom)=&
          pp%upptbl_ao(1:pp%nrps_ao(species),atomic_orbital_ordinals(ll+1,species),species)
      enddo
    enddo
    call evaluate_dg_periodic_sp_projectors(dc%system_tot%primitive_a,lattice_inverse,positions,&
      dc%system_tot%Rion,radial_grid,radial_projector,radial_count,pseudopotential_fingerprint,&
      channels,values,ok,message)
  end subroutine

  subroutine invert_ow_lattice(lattice,inverse,determinant,ok)
    real(8),intent(in)::lattice(3,3)
    real(8),intent(out)::inverse(3,3),determinant
    logical,intent(out)::ok
    determinant=lattice(1,1)*(lattice(2,2)*lattice(3,3)-lattice(2,3)*lattice(3,2))-&
      lattice(1,2)*(lattice(2,1)*lattice(3,3)-lattice(2,3)*lattice(3,1))+&
      lattice(1,3)*(lattice(2,1)*lattice(3,2)-lattice(2,2)*lattice(3,1))
    ok=abs(determinant)>1d-14
    if(.not.ok)then;inverse=0d0;return;endif
    inverse(1,:)=[lattice(2,2)*lattice(3,3)-lattice(2,3)*lattice(3,2),&
      lattice(1,3)*lattice(3,2)-lattice(1,2)*lattice(3,3),&
      lattice(1,2)*lattice(2,3)-lattice(1,3)*lattice(2,2)]/determinant
    inverse(2,:)=[lattice(2,3)*lattice(3,1)-lattice(2,1)*lattice(3,3),&
      lattice(1,1)*lattice(3,3)-lattice(1,3)*lattice(3,1),&
      lattice(1,3)*lattice(2,1)-lattice(1,1)*lattice(2,3)]/determinant
    inverse(3,:)=[lattice(2,1)*lattice(3,2)-lattice(2,2)*lattice(3,1),&
      lattice(1,2)*lattice(3,1)-lattice(1,1)*lattice(3,2),&
      lattice(1,1)*lattice(2,2)-lattice(1,2)*lattice(2,1)]/determinant
  end subroutine

  subroutine build_ow_fragment_permutation_representation(local_target_count,nbox,symmetry_map,&
      product_table,ok,message)
    integer,intent(in)::local_target_count,nbox
    integer(8),intent(in)::symmetry_map(:,:)
    integer,allocatable,intent(out)::product_table(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::rank_fragment(:),target_fragment_local(:),target_fragment_all(:,:)
    integer::nproc,rank,ierr,source_rank,target_rank,operation,iw,target_fragment,&
      left,right,product,middle_fragment,middle_rank,result_fragment

    call MPI_Comm_size(dc%icomm_tot,nproc,ierr);call MPI_Comm_rank(dc%icomm_tot,rank,ierr)
    allocate(rank_fragment(nproc),target_fragment_local(nproc),target_fragment_all(nproc,nproc))
    call MPI_Allgather(dc%i_frag,1,MPI_INTEGER,rank_fragment,1,MPI_INTEGER,dc%icomm_tot,ierr)
    ok=ierr==MPI_SUCCESS.and.size(symmetry_map,2)==nproc
    do operation=1,nproc
      target_fragment_local(operation)=int((symmetry_map(1,operation)-1_8)/int(nbox,8))+1
    enddo
    call MPI_Allgather(target_fragment_local,nproc,MPI_INTEGER,target_fragment_all,nproc,&
      MPI_INTEGER,dc%icomm_tot,ierr)
    ok=ok.and.ierr==MPI_SUCCESS
    allocate(product_table(nproc,nproc));product_table=0
    do left=1,nproc;do right=1,nproc
      do product=1,nproc
        do source_rank=0,nproc-1
          middle_fragment=target_fragment_all(right,source_rank+1)
          middle_rank=findloc(rank_fragment,middle_fragment,dim=1)-1
          if(middle_rank<0)exit
          result_fragment=target_fragment_all(left,middle_rank+1)
          if(result_fragment/=target_fragment_all(product,source_rank+1))exit
        end do
        if(source_rank==nproc)then;product_table(left,right)=product;exit;end if
      end do
      if(product_table(left,right)==0)ok=.false.
    end do;end do
    if(allocated(ow_basis%symmetry_representation))deallocate(ow_basis%symmetry_representation)
    allocate(ow_basis%symmetry_representation(local_target_count*nproc,local_target_count*nproc,nproc))
    ow_basis%symmetry_representation=(0d0,0d0)
    do operation=1,nproc
      do source_rank=0,nproc-1
        target_fragment=target_fragment_all(operation,source_rank+1)
        target_rank=findloc(rank_fragment,target_fragment,dim=1)-1
        if(target_rank<0)then;ok=.false.;cycle;endif
        do iw=1,local_target_count
          ow_basis%symmetry_representation(target_rank*local_target_count+iw,&
            source_rank*local_target_count+iw,operation)=1d0
        enddo
      enddo
    enddo
    if(ok)then;message='';else;message='fragment translation is not a rank permutation';endif
  end subroutine

  subroutine materialize_ow_distributed_core_to_buffer(comm,core_values,core_ids,buffer_ids,&
      buffer_values,ok,message)
    integer,intent(in)::comm
    complex(8),intent(in)::core_values(:,:)
    integer(8),intent(in)::core_ids(:),buffer_ids(:)
    complex(8),allocatable,intent(out)::buffer_values(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::owner_values(:,:)
    integer(8),allocatable::owner_ids(:),sorted_ids(:)
    integer,allocatable::sorted_positions(:),buffer_source(:)
    logical,allocatable::filled(:)
    integer,parameter::core_stream_tile=64
    integer::rank,nproc,ierr,owner,point,source_point,nwann,ncore,tile_first,tile_count

    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    nwann=size(core_values,1);ncore=size(core_ids)
    ok=nwann>0.and.ncore>0.and.size(core_values,2)==ncore.and.size(buffer_ids)>0
    if(.not.ok)then;message='invalid distributed core-to-buffer contract';return;end if
    allocate(buffer_values(nwann,size(buffer_ids)),&
      owner_values(nwann,min(ncore,core_stream_tile)),owner_ids(ncore),&
      sorted_ids(ncore),sorted_positions(ncore),buffer_source(size(buffer_ids)),filled(size(buffer_ids)))
    buffer_values=(0d0,0d0);filled=.false.
    do owner=0,nproc-1
      if(rank==owner)owner_ids=core_ids
      call MPI_Bcast(owner_ids,ncore,MPI_INTEGER8,owner,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;ok=.false.;message='core ID stream failed';return;end if
      sorted_ids=owner_ids;sorted_positions=[(point,point=1,ncore)]
      call sort_ow_id_positions(sorted_ids,sorted_positions)
      do point=1,size(buffer_ids)
        if(filled(point))then
          buffer_source(point)=0
        else
          buffer_source(point)=find_sorted_ow_id(sorted_ids,sorted_positions,buffer_ids(point))
        endif
      end do
      do tile_first=1,ncore,core_stream_tile
        tile_count=min(core_stream_tile,ncore-tile_first+1)
        if(rank==owner)owner_values(:,1:tile_count)=core_values(:,tile_first:tile_first+tile_count-1)
        call MPI_Bcast(owner_values,nwann*tile_count,MPI_DOUBLE_COMPLEX,owner,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;ok=.false.;message='core value stream failed';return;end if
        do point=1,size(buffer_ids)
          source_point=buffer_source(point)
          if(source_point<tile_first.or.source_point>=tile_first+tile_count)cycle
          buffer_values(:,point)=owner_values(:,source_point-tile_first+1);filled(point)=.true.
        end do
      end do
    end do
    ok=all(filled).and.all(ieee_is_finite(real(buffer_values))).and.&
      all(ieee_is_finite(aimag(buffer_values)))
    if(ok)then;message='';else;message='distributed core stream did not cover the local buffer';end if
  end subroutine materialize_ow_distributed_core_to_buffer

  subroutine materialize_ow_global_tails(comm,local_values,local_gradients,local_physical_ids,&
      local_target_count,nproc,ok,message)
    integer,intent(in)::comm,local_target_count,nproc
    complex(8),intent(in)::local_values(:,:),local_gradients(:,:,:)
    integer(8),intent(in)::local_physical_ids(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(8),allocatable::all_physical_ids(:,:)
    integer,allocatable::rank_fragment(:),source_position(:,:),sorted_position(:)
    integer(8),allocatable::representative_centers(:),sorted_physical_ids(:)
    complex(8),allocatable::all_values(:,:,:),local_gradient_axis(:,:),all_gradient_axis(:,:,:)
    integer::nbox,ntarget,source,other,point,source_point,iw,axis,ierr,allocation_status
    logical::global_ok

    nbox=size(local_physical_ids);ntarget=local_target_count*nproc
    ok=size(local_values,1)==local_target_count.and.size(local_values,2)==nbox.and.&
      all(shape(local_gradients)==[3,local_target_count,nbox]).and.&
      ow_global_grid_count>0_8.and.ow_global_grid_count<=int(huge(point),8)
    if(.not.ok)then;message='invalid fragment-local Wannier tail contract';return;endif
    allocate(all_physical_ids(nbox,nproc),rank_fragment(nproc),source_position(nbox,nproc),&
      sorted_physical_ids(nbox),sorted_position(nbox),&
      stat=allocation_status)
    if(allocation_status/=0)then;ok=.false.;message='cannot allocate Wannier tail physical-ID map';return;endif
    call MPI_Allgather(local_physical_ids,nbox,MPI_INTEGER8,all_physical_ids,nbox,MPI_INTEGER8,comm,ierr)
    ok=ierr==MPI_SUCCESS
    call MPI_Allgather(dc%i_frag,1,MPI_INTEGER,rank_fragment,1,MPI_INTEGER,comm,ierr)
    ok=ok.and.ierr==MPI_SUCCESS.and.all(rank_fragment>=1).and.all(rank_fragment<=dc%n_frag)
    do source=1,nproc
      do other=source+1,nproc
        if(rank_fragment(source)==rank_fragment(other))ok=.false.
      enddo
    enddo
    do source=1,nproc
      sorted_physical_ids=all_physical_ids(:,source)
      sorted_position=[(point,point=1,nbox)]
      call sort_ow_id_positions(sorted_physical_ids,sorted_position)
      do point=1,nbox
        if(sorted_physical_ids(point)<1_8.or.sorted_physical_ids(point)>ow_global_grid_count)then
          ok=.false.;cycle
        endif
        if(point>1)then
          if(sorted_physical_ids(point)==sorted_physical_ids(point-1))ok=.false.
        endif
        source_position(point,source)=find_sorted_ow_id(sorted_physical_ids,sorted_position,&
          local_physical_ids(point))
      enddo
    enddo
    call comm_logical_and(ok,global_ok,comm);ok=global_ok
    if(.not.ok)then;message='duplicate or invalid physical ID in fragment Wannier tail';return;endif

    allocate(ow_box_values(ntarget,nbox),ow_box_gradients(3,ntarget,nbox),&
      all_values(local_target_count,nbox,nproc),stat=allocation_status)
    if(allocation_status/=0)then;ok=.false.;message='cannot allocate retained global Wannier tails';return;endif
    ow_box_values=(0d0,0d0);ow_box_gradients=(0d0,0d0)
    call MPI_Allgather(local_values,local_target_count*nbox,MPI_DOUBLE_COMPLEX,all_values,&
      local_target_count*nbox,MPI_DOUBLE_COMPLEX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='cannot gather retained Wannier values';return;endif
    do source=1,nproc
      do point=1,nbox
        source_point=source_position(point,source)
        if(source_point==0)cycle
        do iw=1,local_target_count
          ow_box_values((source-1)*local_target_count+iw,point)=all_values(iw,source_point,source)
        enddo
      enddo
    enddo
    deallocate(all_values)

    allocate(local_gradient_axis(local_target_count,nbox),&
      all_gradient_axis(local_target_count,nbox,nproc),stat=allocation_status)
    if(allocation_status/=0)then;ok=.false.;message='cannot allocate retained Wannier gradients';return;endif
    do axis=1,3
      local_gradient_axis=local_gradients(axis,:,:)
      call MPI_Allgather(local_gradient_axis,local_target_count*nbox,MPI_DOUBLE_COMPLEX,all_gradient_axis,&
        local_target_count*nbox,MPI_DOUBLE_COMPLEX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;ok=.false.;message='cannot gather retained Wannier gradients';return;endif
      do source=1,nproc
        do point=1,nbox
          source_point=source_position(point,source)
          if(source_point==0)cycle
          do iw=1,local_target_count
            ow_box_gradients(axis,(source-1)*local_target_count+iw,point)=&
              all_gradient_axis(iw,source_point,source)
          enddo
        enddo
      enddo
    enddo

    if(allocated(ow_basis%center_owner_rank))deallocate(ow_basis%center_owner_rank)
    if(allocated(ow_basis%center_owner_fragment))deallocate(ow_basis%center_owner_fragment)
    if(.not.allocated(ow_basis%center_box_point_ids).or.&
        size(ow_basis%center_box_point_ids)/=local_target_count)then
      ok=.false.;message='representative Wannier center payload is incomplete';return
    endif
    allocate(representative_centers,source=ow_basis%center_box_point_ids)
    deallocate(ow_basis%center_box_point_ids)
    allocate(ow_basis%center_box_point_ids(ntarget))
    allocate(ow_basis%center_owner_rank(ntarget),ow_basis%center_owner_fragment(ntarget))
    do source=1,nproc
      do iw=1,local_target_count
        ow_basis%center_box_point_ids((source-1)*local_target_count+iw)=&
          int(rank_fragment(source)-1,8)*int(nbox,8)+representative_centers(iw)
      enddo
      ow_basis%center_owner_rank((source-1)*local_target_count+1:source*local_target_count)=source-1
      ow_basis%center_owner_fragment((source-1)*local_target_count+1:source*local_target_count)=&
        rank_fragment(source)
    enddo
    ok=all(ieee_is_finite(real(ow_box_values))).and.all(ieee_is_finite(aimag(ow_box_values))).and.&
      all(ieee_is_finite(real(ow_box_gradients))).and.all(ieee_is_finite(aimag(ow_box_gradients)))
    if(ok)then;message='';else;message='nonfinite retained global Wannier tail';endif
  end subroutine

  subroutine sort_ow_id_positions(ids,positions)
    integer(8),intent(inout)::ids(:)
    integer,intent(inout)::positions(:)
    integer(8),allocatable::work_ids(:)
    integer,allocatable::work_positions(:)
    integer::width,left,middle,right,i,j,k,n
    logical::choose_left
    n=size(ids)
    if(size(positions)/=n)error stop 'invalid Wannier physical-ID sort payload'
    allocate(work_ids(n),work_positions(n));width=1
    do while(width<n)
      do left=1,n,2*width
        middle=min(left+width,n+1);right=min(left+2*width,n+1)
        i=left;j=middle
        do k=left,right-1
          if(i>=middle)then
            choose_left=.false.
          else if(j>=right)then
            choose_left=.true.
          else
            choose_left=ids(i)<=ids(j)
          endif
          if(choose_left)then
            work_ids(k)=ids(i);work_positions(k)=positions(i);i=i+1
          else
            work_ids(k)=ids(j);work_positions(k)=positions(j);j=j+1
          endif
        enddo
      enddo
      ids=work_ids;positions=work_positions
      if(width>n/2)then;width=n;else;width=2*width;endif
    enddo
  end subroutine

  integer function find_sorted_ow_id(sorted_ids,sorted_positions,target) result(position)
    integer(8),intent(in)::sorted_ids(:),target
    integer,intent(in)::sorted_positions(:)
    integer::lower,upper,middle
    position=0
    if(size(sorted_positions)/=size(sorted_ids))return
    lower=1;upper=size(sorted_ids)
    do while(lower<=upper)
      middle=lower+(upper-lower)/2
      if(sorted_ids(middle)<target)then
        lower=middle+1
      else if(sorted_ids(middle)>target)then
        upper=middle-1
      else
        position=sorted_positions(middle);return
      endif
    enddo
  end function

  subroutine periodic_box_gradients(values,box_size,gradient_coefficients,gradients)
    complex(8),intent(in)::values(:,:)
    integer,intent(in)::box_size(3)
    real(8),intent(in)::gradient_coefficients(:,:)
    complex(8),intent(out)::gradients(:,:,:)
    integer::i,j,k,p,pm,pp,axis,nstate,index(3),minus(3),plus(3),dist,radius
    nstate=size(values,1)
    radius=size(gradient_coefficients,1)
    do k=1,box_size(3);do j=1,box_size(2);do i=1,box_size(1)
      index=[i,j,k];p=i+box_size(1)*((j-1)+box_size(2)*(k-1))
      do axis=1,3
        gradients(axis,:,p)=(0d0,0d0)
        if(index(axis)>radius.and.index(axis)<=box_size(axis)-radius)then
          do dist=1,radius
            minus=index;plus=index;minus(axis)=index(axis)-dist;plus(axis)=index(axis)+dist
            pm=minus(1)+box_size(1)*((minus(2)-1)+box_size(2)*(minus(3)-1))
            pp=plus(1)+box_size(1)*((plus(2)-1)+box_size(2)*(plus(3)-1))
            gradients(axis,:,p)=gradients(axis,:,p)+gradient_coefficients(dist,axis)*&
              (values(:,pp)-values(:,pm))
          enddo
        else if(ow_core_size(axis)==dc%lg_tot%num(axis).and.ow_buffer(axis)==0)then
          do dist=1,radius
            minus=index;plus=index
            minus(axis)=modulo(index(axis)-dist-1,box_size(axis))+1
            plus(axis)=modulo(index(axis)+dist-1,box_size(axis))+1
            pm=minus(1)+box_size(1)*((minus(2)-1)+box_size(2)*(minus(3)-1))
            pp=plus(1)+box_size(1)*((plus(2)-1)+box_size(2)*(plus(3)-1))
            gradients(axis,:,p)=gradients(axis,:,p)+gradient_coefficients(dist,axis)*&
              (values(:,pp)-values(:,pm))
          enddo
        else
          minus=index;plus=index
          minus(axis)=max(1,index(axis)-1);plus(axis)=min(box_size(axis),index(axis)+1)
          pm=minus(1)+box_size(1)*((minus(2)-1)+box_size(2)*(minus(3)-1))
          pp=plus(1)+box_size(1)*((plus(2)-1)+box_size(2)*(plus(3)-1))
          gradients(axis,:,p)=(values(:,pp)-values(:,pm))/&
            (real(plus(axis)-minus(axis),8)*system%hgs(axis))
        endif
      enddo
    enddo;enddo;enddo
  end subroutine

  subroutine measure_ow_discrete_gradient_map_commutator(comm,values,gradients,weights,&
      physical_ids,target_rows,grid_size,gradient_coefficients,rotations,residual,ok,message)
    integer,intent(in)::comm,grid_size(3)
    complex(8),intent(in)::values(:,:),gradients(:,:,:)
    real(8),intent(in)::weights(:),gradient_coefficients(:,:),rotations(:,:,:)
    integer(8),intent(in)::physical_ids(:),target_rows(:,:)
    real(8),allocatable,intent(out)::residual(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::image(:,:),mapped_gradient(:,:),expected(:,:),differentiated_image(:,:),&
      plus_values(:,:),minus_values(:,:)
    integer(8),allocatable::plus_ids(:),minus_ids(:)
    real(8)::local_norms(2),global_norms(2)
    integer::nstate,nlocal,noperation,operation,axis,source_axis,distance,point,status,ierr
    integer::index(3)
    logical::collective_ok

    nstate=size(values,1);nlocal=size(values,2);noperation=size(target_rows,2)
    ok=nstate>0.and.nlocal>0.and.noperation>0.and.size(gradients,1)==3.and.&
      size(gradients,2)==nstate.and.size(gradients,3)==nlocal.and.size(weights)==nlocal.and.&
      size(physical_ids)==nlocal.and.size(target_rows,1)==nlocal.and.&
      all(shape(rotations)==[3,3,noperation]).and.size(gradient_coefficients,2)==3
    call comm_logical_and(ok,collective_ok,comm)
    if(.not.collective_ok)then;ok=.false.;message='invalid finite-difference/map commutator contract';return;endif
    allocate(residual(noperation),image(nstate,nlocal),mapped_gradient(nstate,nlocal),&
      expected(nstate,nlocal),differentiated_image(nstate,nlocal),plus_ids(nlocal),minus_ids(nlocal),&
      stat=status)
    call comm_logical_and(status==0,collective_ok,comm)
    if(.not.collective_ok)then
      if(allocated(residual))deallocate(residual)
      if(allocated(image))deallocate(image)
      if(allocated(mapped_gradient))deallocate(mapped_gradient)
      if(allocated(expected))deallocate(expected)
      if(allocated(differentiated_image))deallocate(differentiated_image)
      if(allocated(plus_ids))deallocate(plus_ids)
      if(allocated(minus_ids))deallocate(minus_ids)
      ok=.false.;message='finite-difference/map commutator allocation failed';return
    endif
    do operation=1,noperation
      call exchange_dg_point_permuted_orbital_rows(comm,values,target_rows(:,operation),image,ok,message)
      if(.not.ok)return
      local_norms=0d0
      do axis=1,3
        expected=(0d0,0d0)
        do source_axis=1,3
          call exchange_dg_point_permuted_orbital_rows(comm,gradients(source_axis,:,:),&
            target_rows(:,operation),mapped_gradient,ok,message)
          if(.not.ok)return
          expected=expected+rotations(source_axis,axis,operation)*mapped_gradient
        enddo
        differentiated_image=(0d0,0d0)
        do distance=1,size(gradient_coefficients,1)
          do point=1,nlocal
            index(1)=int(modulo(physical_ids(point)-1_8,int(grid_size(1),8)))
            index(2)=int(modulo((physical_ids(point)-1_8)/int(grid_size(1),8),int(grid_size(2),8)))
            index(3)=int((physical_ids(point)-1_8)/int(grid_size(1)*grid_size(2),8))
            index(axis)=modulo(index(axis)+distance,grid_size(axis))
            plus_ids(point)=1_8+int(index(1),8)+int(grid_size(1),8)*(&
              int(index(2),8)+int(grid_size(2),8)*int(index(3),8))
            index(axis)=modulo(index(axis)-2*distance,grid_size(axis))
            minus_ids(point)=1_8+int(index(1),8)+int(grid_size(1),8)*(&
              int(index(2),8)+int(grid_size(2),8)*int(index(3),8))
          enddo
          call materialize_ow_distributed_core_to_buffer(comm,image,physical_ids,plus_ids,&
            plus_values,ok,message)
          if(ok)call materialize_ow_distributed_core_to_buffer(comm,image,physical_ids,minus_ids,&
            minus_values,ok,message)
          if(.not.ok)return
          differentiated_image=differentiated_image+gradient_coefficients(distance,axis)*&
            (plus_values-minus_values)
          deallocate(plus_values,minus_values)
        enddo
        do point=1,nlocal
          local_norms(1)=local_norms(1)+weights(point)*&
            sum(abs(differentiated_image(:,point)-expected(:,point))**2)
          local_norms(2)=local_norms(2)+weights(point)*sum(abs(expected(:,point))**2)
        enddo
      enddo
      call MPI_Allreduce(local_norms,global_norms,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;ok=.false.;message='finite-difference/map reduction failed';return;endif
      residual(operation)=sqrt(max(0d0,global_norms(1))/max(tiny(1d0),global_norms(2)))
    enddo
    ok=all(residual>=0d0.and.residual<huge(1d0))
    if(ok)then;message='';else;message='finite-difference/map residual is not finite';endif
  end subroutine measure_ow_discrete_gradient_map_commutator

  subroutine build_dc_translation_symmetry_map(nbox,symmetry_map,ok,message)
    integer,intent(in)::nbox
    integer(8),intent(out)::symmetry_map(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::rank_fragments(:)
    integer::source_fragment,operation,target_fragment,candidate_fragment,axis,p,match_count,ierr
    integer::translation(3),target_origin(3)
    logical::local_ok,global_ok

    allocate(rank_fragments(dc%n_frag))
    call MPI_Allgather(dc%i_frag,1,MPI_INTEGER,rank_fragments,1,MPI_INTEGER,dc%icomm_tot,ierr)
    local_ok=ierr==MPI_SUCCESS.and.nbox>0.and.size(symmetry_map,1)==nbox.and.&
      size(symmetry_map,2)==dc%n_frag.and.dc%n_frag==size(rank_fragments)
    do p=1,dc%n_frag
      local_ok=local_ok.and.count(rank_fragments==p)==1
    enddo
    local_ok=local_ok.and.all(dc%nxyz_domain_frag==spread(dc%nxyz_domain_frag(:,1),2,dc%n_frag))
    source_fragment=dc%i_frag
    if(local_ok)then
      do operation=1,dc%n_frag
        translation=dc%ixyz_frag(:,operation)-dc%ixyz_frag(:,1)
        target_origin=1+modulo(dc%ixyz_frag(:,source_fragment)-1+translation,dc%lg_tot%num)
        match_count=0;target_fragment=0
        do candidate_fragment=1,dc%n_frag
          if(all(1+modulo(dc%ixyz_frag(:,candidate_fragment)-1,dc%lg_tot%num)==target_origin))then
            match_count=match_count+1;target_fragment=candidate_fragment
          endif
        enddo
        if(match_count/=1)then
          local_ok=.false.;exit
        endif
        do p=1,nbox
          symmetry_map(p,operation)=int((target_fragment-1)*nbox+p,8)
        enddo
      enddo
    endif
    call comm_logical_and(local_ok,global_ok,dc%icomm_tot)
    ok=global_ok
    if(ok)then
      message=''
    else
      message='DC fragment origins do not form a uniform periodic translation group'
    endif
  end subroutine

  subroutine replicate_ow_global_symmetry_orbit(values,gradients,center_box_ids,&
      residual,correction,ok,message)
    complex(8),intent(inout)::values(:,:),gradients(:,:,:)
    integer(8),intent(inout)::center_box_ids(:)
    real(8),intent(out)::residual,correction
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(t_sawf_crystallographic_catalog)::catalog
    complex(8),allocatable::reference_values(:,:),reference_gradients(:,:,:),mapped_values(:,:),&
      mapped_gradients(:,:,:)
    real(8),allocatable::positions(:,:)
    integer,allocatable::species(:),source_to_target(:),point_map(:),fragment_maps(:,:),&
      fragment_orbit(:),orbit_representative(:),rank_fragment(:)
    integer(8),allocatable::reference_centers(:)
    integer::rank,ierr,nwann,nbox,reference_fragment,&
      operation,point,axis,target_axis,atom,max_targets,valid_operation_count,reference_rank
    real(8)::lattice_inverse(3,3),determinant,grid_residual,center_grid(3),local_correction
    logical::inverse_ok,catalog_ok,grid_ok,fragment_ok,center_available,map_ok,found
    character(256)::detail

    ok=.false.;message='';residual=huge(1d0);correction=huge(1d0)
    call MPI_Comm_rank(dc%icomm_tot,rank,ierr);nwann=size(values,1);nbox=size(values,2)
    if(nwann<1.or.nbox/=product(ow_box_size).or.size(center_box_ids)/=nwann.or.&
        any(shape(gradients)/=[3,nwann,nbox]))then
      message='invalid full-system symmetry-orbit replication contract';return
    end if
    allocate(reference_values(nwann,nbox),reference_gradients(3,nwann,nbox),&
      mapped_values(nwann,nbox),mapped_gradients(3,nwann,nbox),reference_centers(nwann))
    call invert_ow_lattice(dc%system_tot%primitive_a,lattice_inverse,determinant,inverse_ok)
    if(.not.inverse_ok)then;message='full-system orbit lattice is singular';return;end if
    allocate(positions(3,dc%system_tot%nion),species(dc%system_tot%nion))
    do atom=1,dc%system_tot%nion
      positions(:,atom)=modulo(matmul(lattice_inverse,dc%system_tot%Rion(:,atom)),1d0)
      species(atom)=dc%system_tot%kion(atom)
    end do
    call load_sawf_crystallographic_catalog_auto(dc%system_tot%primitive_a,positions,species,&
      dg_ow_symmetry_tolerance,catalog,catalog_ok,detail)
    if(.not.catalog_ok)then;message='full-system orbit catalog: '//trim(detail);return;end if
    allocate(fragment_maps(dc%n_frag,size(catalog%operations)),rank_fragment(dc%n_frag))
    valid_operation_count=0
    do operation=1,size(catalog%operations)
      call validate_sawf_fragment_symmetry_map(catalog%operations(operation),dc%lg_tot%num,&
        dc%ixyz_frag-1,dc%nxyz_domain_frag,[0,0,0],dg_ow_symmetry_tolerance,grid_ok,&
        fragment_ok,max_targets,source_to_target,grid_residual,center_available,center_grid,detail)
      if(.not.(grid_ok.and.fragment_ok))cycle
      valid_operation_count=valid_operation_count+1
      fragment_maps(:,valid_operation_count)=source_to_target
    end do
    if(valid_operation_count<1)then;message='no full-system operation preserves the fragment partition';return;end if
    call build_dg_fragment_symmetry_orbits(fragment_maps(:,1:valid_operation_count),&
      fragment_orbit,orbit_representative,catalog_ok,detail)
    if(.not.catalog_ok)then;message='full-system fragment orbit: '//trim(detail);return;end if
    if(rank==0)write(*,'(a,i0,a,*(i0,1x))')&
      '[OW-GS-DIAGNOSTIC] fragment_symmetry_orbit_count=',maxval(fragment_orbit),&
      ' representatives=',pack([(operation,operation=1,dc%n_frag)],&
      orbit_representative==[(operation,operation=1,dc%n_frag)])
    reference_fragment=orbit_representative(dc%i_frag)
    call MPI_Allgather(dc%i_frag,1,MPI_INTEGER,rank_fragment,1,MPI_INTEGER,dc%icomm_tot,ierr)
    reference_rank=findloc(rank_fragment,reference_fragment,dim=1)-1
    if(reference_rank<0)then;message='fragment-orbit representative has no MPI owner';return;end if
    if(rank==reference_rank)then
      reference_values=values;reference_gradients=gradients;reference_centers=center_box_ids
    end if
    call MPI_Bcast(reference_values,size(reference_values),MPI_DOUBLE_COMPLEX,&
      reference_rank,dc%icomm_tot,ierr)
    call MPI_Bcast(reference_gradients,size(reference_gradients),MPI_DOUBLE_COMPLEX,&
      reference_rank,dc%icomm_tot,ierr)
    call MPI_Bcast(reference_centers,size(reference_centers),MPI_INTEGER8,&
      reference_rank,dc%icomm_tot,ierr)
    if(any(reference_centers<1_8).or.any(reference_centers>int(nbox,8)))then
      write(message,'(a,2(i0,1x),a,i0)')'representative center IDs outside local box min/max=',&
        minval(reference_centers),maxval(reference_centers),' nbox=',nbox;return
    end if
    found=.false.
    do operation=1,size(catalog%operations)
      call validate_sawf_fragment_symmetry_map(catalog%operations(operation),dc%lg_tot%num,&
        dc%ixyz_frag-1,dc%nxyz_domain_frag,[0,0,0],dg_ow_symmetry_tolerance,grid_ok,&
        fragment_ok,max_targets,source_to_target,grid_residual,center_available,center_grid,detail)
      if(.not.(grid_ok.and.fragment_ok))cycle
      if(source_to_target(reference_fragment)/=dc%i_frag)cycle
      call build_sawf_fragment_buffer_point_map(catalog%operations(operation),dc%lg_tot%num,&
        dc%ixyz_frag(:,reference_fragment)-1,dc%nxyz_domain_frag(:,reference_fragment),&
        dc%ixyz_frag(:,dc%i_frag)-1,dc%nxyz_domain_frag(:,dc%i_frag),ow_buffer,&
        dg_ow_symmetry_tolerance,point_map,map_ok,detail)
      if(map_ok)then;found=.true.;exit;end if
    end do
    if(.not.found)then;message='no full-system affine operation connects representative fragment';return;end if
    if(any(point_map<1).or.any(point_map>nbox))then
      write(message,'(a,2(i0,1x),a,i0)')'full-system point map outside local box min/max=',&
        minval(point_map),maxval(point_map),' nbox=',nbox;return
    end if
    mapped_values=(0d0,0d0);mapped_gradients=(0d0,0d0)
    do point=1,nbox
      mapped_values(:,point_map(point))=reference_values(:,point)
      do axis=1,3;do target_axis=1,3
        mapped_gradients(target_axis,:,point_map(point))=mapped_gradients(target_axis,:,point_map(point))+&
          catalog%operations(operation)%R(target_axis,axis)*reference_gradients(axis,:,point)
      end do;end do
    end do
    local_correction=max(maxval(abs(values-mapped_values)),maxval(abs(gradients-mapped_gradients)))
    call MPI_Allreduce(local_correction,correction,1,MPI_DOUBLE_PRECISION,MPI_MAX,dc%icomm_tot,ierr)
    values=mapped_values;gradients=mapped_gradients
    do point=1,nwann;center_box_ids(point)=int(point_map(int(reference_centers(point))),8);end do
    residual=0d0;ok=.true.;message=''
  end subroutine replicate_ow_global_symmetry_orbit

  subroutine prepare_ow_exact_fragment_symmetry(local_symmetry_map,point_integer_rotations,point_rotations, &
      point_product,ok,message)
    integer(8),allocatable,intent(out)::local_symmetry_map(:,:)
    integer,allocatable,intent(out)::point_integer_rotations(:,:,:)
    real(8),allocatable,intent(out)::point_rotations(:,:,:)
    integer,allocatable,intent(out)::point_product(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(t_sawf_crystallographic_catalog)::catalog
    real(8),allocatable::fragment_positions(:,:)
    integer,allocatable::fragment_species(:),fragment_atom_index(:),selected(:),source_to_target(:),point_map(:)
    logical,allocatable::fragment_atom_mask(:),operation_allowed(:)
    real(8)::lattice_inverse(3,3),determinant,fragment_center(3),grid_residual,center_grid(3),site_residual
    integer::atom,axis,iop,ifrag,nfragment_atoms,max_targets,relative_index(3),grid_index(3)
    logical::inverse_ok,grid_ok,fragment_ok,center_available,map_ok
    character(256)::detail

    ok=.false.;message='';ifrag=dc%i_frag
    call invert_ow_lattice(dc%system_tot%primitive_a,lattice_inverse,determinant,inverse_ok)
    if(.not.inverse_ok)then;message='exact fragment symmetry lattice is singular';return;end if
    allocate(fragment_atom_mask(dc%system_tot%nion),fragment_atom_index(dc%system_tot%nion))
    fragment_atom_mask=.false.;fragment_atom_index=0;nfragment_atoms=0
    do atom=1,dc%system_tot%nion
      grid_index=modulo(floor(modulo(matmul(lattice_inverse,dc%system_tot%Rion(:,atom)),1d0)* &
        real(dc%lg_tot%num,8)),dc%lg_tot%num)
      do axis=1,3
        relative_index(axis)=modulo(grid_index(axis)-(dc%ixyz_frag(axis,ifrag)-1)+ow_buffer(axis), &
          dc%lg_tot%num(axis))-ow_buffer(axis)
      end do
      if(all(relative_index>=-ow_buffer).and.all(relative_index<ow_core_size+ow_buffer))then
        nfragment_atoms=nfragment_atoms+1;fragment_atom_mask(atom)=.true.
        fragment_atom_index(nfragment_atoms)=atom
      end if
    end do
    if(nfragment_atoms<1)then;message='exact buffered fragment contains no instantaneous atom';return;end if
    allocate(fragment_positions(3,nfragment_atoms),fragment_species(nfragment_atoms))
    do atom=1,nfragment_atoms
      fragment_positions(:,atom)=modulo(matmul(lattice_inverse,&
        dc%system_tot%Rion(:,fragment_atom_index(atom))),1d0)
      fragment_species(atom)=dc%system_tot%kion(fragment_atom_index(atom))
    end do
    call load_sawf_crystallographic_catalog_auto(dc%system_tot%primitive_a,fragment_positions, &
      fragment_species,dg_ow_symmetry_tolerance,catalog,ok,detail)
    if(.not.ok)then;message='fragment crystallographic catalog: '//trim(detail);return;end if
    allocate(operation_allowed(size(catalog%operations)));operation_allowed=.false.
    do iop=1,size(catalog%operations)
      call validate_sawf_fragment_symmetry_map(catalog%operations(iop),dc%lg_tot%num,&
        dc%ixyz_frag-1,dc%nxyz_domain_frag,[0,0,0],dg_ow_symmetry_tolerance,grid_ok,&
        fragment_ok,max_targets,source_to_target,grid_residual,center_available,center_grid,detail)
      operation_allowed(iop)=grid_ok.and.fragment_ok
      if(operation_allowed(iop))operation_allowed(iop)=source_to_target(ifrag)==ifrag
    end do
    fragment_center=modulo((real(dc%ixyz_frag(:,ifrag)-1,8)+0.5d0*real(ow_core_size,8))/ &
      real(dc%lg_tot%num,8),1d0)
    call build_dg_fragment_site_stabilizer(catalog%integer_rotation,catalog%fractional_translation,&
      fragment_center,operation_allowed,dg_ow_symmetry_tolerance,selected,point_product,&
      site_residual,ok,detail)
    if(.not.ok)then;message='fragment site stabilizer: '//trim(detail);return;end if
    allocate(local_symmetry_map(product(ow_box_size),size(selected)),&
      point_integer_rotations(3,3,size(selected)),point_rotations(3,3,size(selected)))
    do iop=1,size(selected)
      call build_sawf_fragment_buffer_point_map(catalog%operations(selected(iop)),dc%lg_tot%num,&
        dc%ixyz_frag(:,ifrag)-1,dc%nxyz_domain_frag(:,ifrag),&
        dc%ixyz_frag(:,ifrag)-1,dc%nxyz_domain_frag(:,ifrag),ow_buffer,&
        dg_ow_symmetry_tolerance,point_map,map_ok,detail)
      if(.not.map_ok)then;message='fragment point map: '//trim(detail);ok=.false.;return;end if
      local_symmetry_map(:,iop)=int(point_map,8)
      point_integer_rotations(:,:,iop)=catalog%integer_rotation(:,:,selected(iop))
      point_rotations(:,:,iop)=catalog%operations(selected(iop))%R
    end do
    write(*,'(a,i0,a,i0,a,i0,a,i0,2a,a,es12.4)')&
      '[OW-GS-DIAGNOSTIC] fragment=',ifrag,' exact_site_group_order=',size(selected),&
      ' space_group_number=',catalog%space_group_number,' hall_number=',catalog%hall_number,&
      ' point_group_symbol=',trim(catalog%point_group_symbol),' site_residual=',site_residual
    ok=.true.;message=''
  end subroutine prepare_ow_exact_fragment_symmetry

  subroutine restore_ow_checkpoint_density(checkpoint,ok,message)
    type(s_dg_overlapping_wannier_checkpoint),intent(in)::checkpoint
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8),allocatable::local_density(:),global_density(:),density4(:,:,:,:)
    integer::p
    allocate(local_density(int(ow_global_grid_count)),global_density(int(ow_global_grid_count)))
    local_density=0d0
    do p=1,size(checkpoint%core_physical_ids)
      local_density(int(checkpoint%core_physical_ids(p)))=checkpoint%density(p)
    enddo
    call comm_summation(local_density,global_density,size(local_density),dc%icomm_tot)
    allocate(density4(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),1))
    density4(:,:,:,1)=reshape(global_density,dc%lg_tot%num)
    call dg_dc_update_potential_from_density(density4,ok,message)
  end subroutine

  subroutine diagnose_ow_total_nonlocal_projector_range(centers,integer_rotation,cartesian_rotation,&
      translation,wannier_representation,fragment_shape,ok,message)
    real(8),intent(in)::centers(:,:),cartesian_rotation(3,3),translation(3)
    integer,intent(in)::integer_rotation(3,3),fragment_shape(3)
    complex(8),intent(in)::wannier_representation(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_dg_nonlocal_range_receipt)::receipt
    integer::rank,ierr,nx,ny,owned_count,ix,iy,iz,g,p,q,ia,ik,ll,l,l0,m,ilma,lm,radial,&
      nnz,target_atom,target_channel,source_first,target_first,allocation_status
    integer,allocatable::offsets(:),projector_atom(:),channel_l(:),channel_m(:),channel_radial(:),atom_map(:)
    integer(8),allocatable::support_ids(:)
    complex(8),allocatable::support_values(:),mg_wannier(:,:),projector_representation(:,:)
    real(8),allocatable::strength(:),fractional_atoms(:,:)
    real(8)::inverse(3,3),determinant,delta(3),cart_delta(3),distance,best,pblock(3,3),basis(3,3)
    logical::redistribution_ok
    character(256)::redistribution_message
    integer(8)::redistribution_workspace

    ok=.false.;message='';call MPI_Comm_rank(dc%icomm_tot,rank,ierr)
    if(ierr/=MPI_SUCCESS.or.size(centers,1)/=3.or.size(centers,2)/=size(ow_core_values,1).or.&
      any(shape(wannier_representation)/=[size(ow_core_values,1),size(ow_core_values,1)]).or.&
      dc%ppg_tot%Nlma<1.or.any(fragment_shape<1))then
      message='invalid total nonlocal projector diagnostic contract';return
    endif
    call invert_ow_lattice(dc%system_tot%primitive_a,inverse,determinant,redistribution_ok)
    if(.not.redistribution_ok)return
    nx=dc%lg_tot%num(1);ny=dc%lg_tot%num(2)
    owned_count=product(dc%mg_tot%ie-dc%mg_tot%is+1)
    if(.not.allocated(ow_hpsi_grid_ids))then
      allocate(ow_hpsi_grid_ids(owned_count),stat=allocation_status)
      if(allocation_status/=0)then;message='cannot allocate total nonlocal grid IDs';return;endif
      g=0
      do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3);do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
      do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
        g=g+1;ow_hpsi_grid_ids(g)=int(ix,8)+int(nx,8)*(int(iy-1,8)+int(ny,8)*int(iz-1,8))
      enddo;enddo;enddo
    endif
    if(.not.ow_hpsi_redistribution%initialized)then
      call initialize_dg_full_cell_redistribution(dc%icomm_tot,int(ow_global_grid_count),ow_core_ids,&
        ow_hpsi_grid_ids,ow_hpsi_redistribution,redistribution_workspace,redistribution_ok,&
        redistribution_message)
      if(.not.redistribution_ok)then;message=trim(redistribution_message);return;endif
      ow_hpsi_redistribution_workspace=redistribution_workspace
    endif
    allocate(mg_wannier(size(ow_core_values,1),owned_count),stat=allocation_status)
    if(allocation_status/=0)then;message='cannot allocate redistributed diagnostic Wannier tile';return;endif
    call apply_dg_full_cell_redistribution_forward(ow_hpsi_redistribution,ow_core_values,mg_wannier,&
      redistribution_ok,redistribution_message)
    if(.not.redistribution_ok)then;message=trim(redistribution_message);return;endif

    allocate(offsets(dc%ppg_tot%Nlma+1),projector_atom(dc%ppg_tot%Nlma),&
      channel_l(dc%ppg_tot%Nlma),channel_m(dc%ppg_tot%Nlma),channel_radial(dc%ppg_tot%Nlma),&
      strength(dc%ppg_tot%Nlma),fractional_atoms(3,dc%system_tot%nion),atom_map(dc%system_tot%nion))
    channel_l=-1;channel_m=0;channel_radial=0
    do ia=1,dc%system_tot%nion
      fractional_atoms(:,ia)=modulo(matmul(inverse,dc%system_tot%Rion(:,ia)),1d0)
      ik=dc%system_tot%kion(ia);lm=0;l0=0
      do ll=0,pp%mlps(ik)
        radial=0
        do l=l0,l0+pp%nproj(ll,ik)-1
          if(pp%inorm(l,ik)==0)cycle
          radial=radial+1
          do m=-ll,ll
            lm=lm+1;ilma=dc%ppg_tot%lma_tbl(lm,ia)
            channel_l(ilma)=ll;channel_m(ilma)=m;channel_radial(ilma)=radial
          enddo
        enddo
        l0=l
      enddo
    enddo
    if(any(channel_l<0).or.any(channel_l>1))then
      message='total nonlocal diagnostic currently requires complete s/p pseudopotential shells';return
    endif
    do ia=1,dc%system_tot%nion
      best=huge(1d0);target_atom=0
      do ik=1,dc%system_tot%nion
        if(dc%system_tot%kion(ik)/=dc%system_tot%kion(ia))cycle
        delta=matmul(real(integer_rotation,8),fractional_atoms(:,ia))+translation-fractional_atoms(:,ik)
        delta=delta-anint(delta);cart_delta=matmul(dc%system_tot%primitive_a,delta)
        distance=sqrt(sum(cart_delta*cart_delta))
        if(distance<best)then;best=distance;target_atom=ik;endif
      enddo
      if(target_atom<1.or.best>max(1d-8,dg_ow_symmetry_tolerance))then
        message='operation 5 has no exact same-species total-system atom partner';return
      endif
      atom_map(ia)=target_atom
    enddo
    do ia=1,size(atom_map)-1
      do q=ia+1,size(atom_map)
        if(atom_map(ia)==atom_map(q))then
          message='operation 5 total-system atom map is not one-to-one';return
        endif
      enddo
    enddo

    nnz=0;offsets(1)=1
    do ilma=1,dc%ppg_tot%Nlma
      ia=dc%ppg_tot%ia_tbl(ilma)
      do q=1,dc%ppg_tot%mps(ia)
        ix=dc%ppg_tot%jxyz(1,q,ia);iy=dc%ppg_tot%jxyz(2,q,ia);iz=dc%ppg_tot%jxyz(3,q,ia)
        if(ix<dc%mg_tot%is(1).or.ix>dc%mg_tot%ie(1).or.iy<dc%mg_tot%is(2).or.&
          iy>dc%mg_tot%ie(2).or.iz<dc%mg_tot%is(3).or.iz>dc%mg_tot%ie(3))cycle
        nnz=nnz+1
      enddo
      offsets(ilma+1)=nnz+1
    enddo
    allocate(support_ids(nnz),support_values(nnz),projector_representation(dc%ppg_tot%Nlma,dc%ppg_tot%Nlma))
    nnz=0;projector_representation=0d0
    do ilma=1,dc%ppg_tot%Nlma
      ia=dc%ppg_tot%ia_tbl(ilma);projector_atom(ilma)=ia
      strength(ilma)=system%hvol*dc%ppg_tot%rinv_uvu(ilma)
      do q=1,dc%ppg_tot%mps(ia)
        ix=dc%ppg_tot%jxyz(1,q,ia);iy=dc%ppg_tot%jxyz(2,q,ia);iz=dc%ppg_tot%jxyz(3,q,ia)
        if(ix<dc%mg_tot%is(1).or.ix>dc%mg_tot%ie(1).or.iy<dc%mg_tot%is(2).or.&
          iy>dc%mg_tot%ie(2).or.iz<dc%mg_tot%is(3).or.iz>dc%mg_tot%ie(3))cycle
        nnz=nnz+1;support_ids(nnz)=int(ix,8)+int(nx,8)*(int(iy-1,8)+int(ny,8)*int(iz-1,8))
        support_values(nnz)=dc%ppg_tot%zekr_uV(q,ilma,1)
      enddo
    enddo
    basis=0d0;basis(2,1)=-1d0;basis(3,2)=1d0;basis(1,3)=-1d0
    pblock=matmul(transpose(basis),matmul(cartesian_rotation,basis))
    do ilma=1,dc%ppg_tot%Nlma
      target_atom=atom_map(projector_atom(ilma));target_channel=0
      do q=1,dc%ppg_tot%Nlma
        if(dc%ppg_tot%ia_tbl(q)==target_atom.and.channel_l(q)==channel_l(ilma).and.&
          channel_radial(q)==channel_radial(ilma))then
          if(channel_l(ilma)==0.and.channel_m(q)==0)target_channel=q
          if(channel_l(ilma)==1)then
            projector_representation(q,ilma)=pblock(channel_m(q)+2,channel_m(ilma)+2)
            target_channel=q
          endif
        endif
      enddo
      if(channel_l(ilma)==0.and.target_channel>0)projector_representation(target_channel,ilma)=1d0
      if(target_channel==0)then;message='operation 5 has an unmatched radial projector channel';return;endif
    enddo
    call analyze_dg_nonlocal_projector_range(dc%icomm_tot,ow_hpsi_grid_ids,mg_wannier,centers,&
      dc%system_tot%primitive_a,dc%system_tot%Rion,dc%system_tot%kion,projector_atom,strength,&
      offsets,support_ids,support_values,integer_rotation,translation,fragment_shape,16,&
      wannier_representation,projector_representation,receipt,ok,message,ow_row_ids,ow_direct_nonlocal_rows)
    if(ok.and.rank==0)write(*,'(a,5(a,es16.8),2(a,i0))')&
      '[OW-GS-DIAGNOSTIC] total_nonlocal_projector_range_operation5',&
      ' local=',receipt%local_contribution,' adjacent=',receipt%adjacent_contribution,&
      ' remote=',receipt%remote_contribution,' max_remote_fraction=',receipt%maximum_remote_fraction,&
      ' covariance=',receipt%symmetry_pair_defect,' unmatched=',receipt%unmatched_channel_count,&
      ' workspace_peak_bytes=',receipt%workspace_peak_bytes
  end subroutine diagnose_ow_total_nonlocal_projector_range

  subroutine ow_build_hamiltonian(comm,density,hrows,new_potential,fingerprint,ok,message,update_auxiliary_pencil)
    integer,intent(in)::comm
    real(8),intent(in)::density(:)
    complex(8),intent(out)::hrows(:,:)
    real(8),intent(out)::new_potential(:)
    integer(8),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,intent(in),optional::update_auxiliary_pencil
    real(8),allocatable::global_density(:),summed_density(:),core_potential(:),box_potential(:)
    complex(8),allocatable::kinetic_rows(:,:),local_rows(:,:),nonlocal_rows(:,:),boundary_rows(:,:),&
      full_rows(:,:),&
      core_potential_values(:,:),box_potential_values(:,:),sym_h_rows(:,:),sym_s_rows(:,:),sym_rho_rows(:,:)
    complex(8),allocatable::component_rows(:,:,:)
    real(8)::kinetic_scale,local_scale,nonlocal_scale,hamiltonian_scale,&
      stitched_t_hermiticity,stitched_v_hermiticity,weight_gradient_trace,&
      local_direct_difference,global_direct_difference,direct_nonlocal_scale
    real(8)::pencil_before(3),pencil_after(3),boundary_artifact_change,boundary_artifact_magnitude
    real(8)::component_covariance(3)
    logical::finite_t,finite_local,finite_nonlocal,finite_h,update_auxiliary
    integer::p,ix,iy,iz,nwann,owned_projectors,rank,ierr
    integer(8)::stitched_operator_peak_elements
    integer(8)::full_cell_workspace_peak
    integer(8)::pencil_symmetry_workspace_peak
    call MPI_Comm_rank(comm,rank,ierr)
    update_auxiliary=.true.;if(present(update_auxiliary_pencil))update_auxiliary=update_auxiliary_pencil
    allocate(global_density(int(ow_global_grid_count)),summed_density(int(ow_global_grid_count)))
    global_density=0d0
    do p=1,size(ow_core_ids);global_density(int(ow_core_ids(p)))=density(p);enddo
    call comm_summation(global_density,summed_density,size(global_density),comm)
    if(allocated(ow_work_density))deallocate(ow_work_density)
    allocate(ow_work_density(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),1))
    ow_work_density(:,:,:,1)=reshape(summed_density,dc%lg_tot%num)
    call dg_dc_update_potential_from_density(ow_work_density,ok,message)
    if(.not.ok)return
    nwann=size(ow_core_values,1);allocate(core_potential(size(ow_core_ids)))
    do p=1,size(ow_core_ids)
      ix=int(modulo(ow_core_ids(p)-1_8,int(dc%lg_tot%num(1),8)))+1
      iy=int(modulo((ow_core_ids(p)-1_8)/int(dc%lg_tot%num(1),8),&
        int(dc%lg_tot%num(2),8)))+1
      iz=int((ow_core_ids(p)-1_8)/(int(dc%lg_tot%num(1),8)*int(dc%lg_tot%num(2),8)))+1
      core_potential(p)=dc%vloc_tot(1)%f(ix,iy,iz)
    enddo
    new_potential=core_potential
    ow_full_cell_component_mode=0
    call project_dg_full_cell_hamiltonian_tiles(comm,int(ow_global_grid_count),ow_core_ids,&
      ow_core_weights,ow_core_values,ow_row_ids,min(16,nwann),apply_ow_full_cell_hpsi_tile,&
      full_rows,full_cell_workspace_peak,ok,message)
    if(.not.ok)return
    ow_full_cell_component_mode=1
    call project_dg_full_cell_hamiltonian_tiles(comm,int(ow_global_grid_count),ow_core_ids,&
      ow_core_weights,ow_core_values,ow_row_ids,min(16,nwann),apply_ow_full_cell_hpsi_tile,&
      kinetic_rows,stitched_operator_peak_elements,ok,message)
    if(.not.ok)return
    ow_full_cell_component_mode=2
    call project_dg_full_cell_hamiltonian_tiles(comm,int(ow_global_grid_count),ow_core_ids,&
      ow_core_weights,ow_core_values,ow_row_ids,min(16,nwann),apply_ow_full_cell_hpsi_tile,&
      local_rows,stitched_operator_peak_elements,ok,message)
    ow_full_cell_component_mode=0
    if(.not.ok)return
    allocate(nonlocal_rows(size(ow_row_ids),nwann),boundary_rows(size(ow_row_ids),nwann))
    nonlocal_rows=full_rows-kinetic_rows-local_rows;boundary_rows=(0d0,0d0)
    hrows=full_rows
    if(rank==0)write(*,'(a,i0)')'[OW-GS-DIAGNOSTIC] full_cell_hpsi_workspace_peak_bytes=',&
      full_cell_workspace_peak
    allocate(component_rows(size(ow_row_ids),nwann,3))
    component_rows(:,:,1)=kinetic_rows;component_rows(:,:,2)=local_rows
    component_rows(:,:,3)=nonlocal_rows
    ! Measure the raw projected components before the strict covariance gate.
    ! A covariance rejection must not hide whether the defect was already a
    ! failure of the underlying Hermitian operator projection.
    call ow_distributed_hermiticity(comm,ow_row_ids,kinetic_rows,ow_diag_t_hermiticity,&
      kinetic_scale,finite_t)
    call ow_distributed_hermiticity(comm,ow_row_ids,local_rows,ow_diag_vlocal_hermiticity,&
      local_scale,finite_local)
    call ow_distributed_hermiticity(comm,ow_row_ids,nonlocal_rows,ow_diag_vnl_hermiticity,&
      nonlocal_scale,finite_nonlocal)
    if(allocated(ow_direct_nonlocal_rows).and..not.ow_direct_nonlocal_compared)then
      local_direct_difference=0d0;direct_nonlocal_scale=0d0
      if(size(nonlocal_rows)>0)then
        local_direct_difference=maxval(abs(nonlocal_rows-ow_direct_nonlocal_rows))
        direct_nonlocal_scale=maxval(abs(ow_direct_nonlocal_rows))
      endif
      call MPI_Allreduce(local_direct_difference,global_direct_difference,1,&
        MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      call MPI_Allreduce(MPI_IN_PLACE,direct_nonlocal_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(rank==0)write(*,'(a,3(a,es16.8))')&
        '[OW-GS-DIAGNOSTIC] total_projector_direct/hpsi_nonlocal_difference',&
        ' difference=',global_direct_difference,' direct_scale=',direct_nonlocal_scale,&
        ' hpsi_scale=',nonlocal_scale
      ow_direct_nonlocal_compared=.true.
    endif
    if(rank==0)write(*,'(a,6(a,es16.8))')'[OW-GS-DIAGNOSTIC] raw projected component Hermiticity',&
      ' kinetic_defect=',ow_diag_t_hermiticity,' kinetic_scale=',kinetic_scale,&
      ' local_defect=',ow_diag_vlocal_hermiticity,' local_scale=',local_scale,&
      ' nonlocal_defect=',ow_diag_vnl_hermiticity,' nonlocal_scale=',nonlocal_scale
    call symmetrize_dg_distributed_pencil_rows(comm,ow_row_ids,hrows,ow_srows,ow_rhorows,&
      boundary_rows,ow_pencil_generator_representation,ow_pencil_generator_operations,&
      ow_pencil_affine_product,ow_pencil_translation_subgroup,ow_pencil_coset_representatives,&
      dg_ow_symmetry_tolerance,sym_h_rows,sym_s_rows,sym_rho_rows,pencil_before,pencil_after,&
      boundary_artifact_change,boundary_artifact_magnitude,pencil_symmetry_workspace_peak,ok,message,&
      component_rows,component_covariance,update_auxiliary)
    deallocate(component_rows)
    if(.not.ok)return
    hrows=sym_h_rows
    if(update_auxiliary)then
      ow_srows=sym_s_rows;ow_rhorows=sym_rho_rows
      ow_hybrid_symmetry_defect=maxval(pencil_after)
    else
      ow_hybrid_symmetry_defect=pencil_after(1)
    endif
    if(rank==0)write(*,'(a,8(a,es16.8),a,i0)')'[OW-GS-DIAGNOSTIC] stitched_pencil_symmetry',&
      ' h_before=',pencil_before(1),' s_before=',pencil_before(2),' rho_before=',pencil_before(3),&
      ' h_after=',pencil_after(1),' s_after=',pencil_after(2),' rho_after=',pencil_after(3),&
      ' artifact_change=',boundary_artifact_change,' artifact_magnitude=',boundary_artifact_magnitude,&
      ' workspace_peak_elements=',pencil_symmetry_workspace_peak
    if(rank==0)write(*,'(a,3(a,es16.8))')'[OW-GS-DIAGNOSTIC] stitched_component_covariance',&
      ' kinetic=',component_covariance(1),' local=',component_covariance(2),&
      ' nonlocal=',component_covariance(3)
    ow_diag_h_local_bytes=max(ow_diag_h_local_bytes,int(size(hrows),8)*16_8)
    call ow_distributed_hermiticity(comm,ow_row_ids,hrows,ow_diag_h_hermiticity,&
      hamiltonian_scale,finite_h)
    if(.not.finite_t.or..not.finite_local.or..not.finite_nonlocal.or..not.finite_h)then
      ok=.false.;message='nonfinite projected DC Hamiltonian component';return
    endif
    if(ow_diag_t_hermiticity>dg_dc_gs_hermiticity_tolerance*max(1d0,kinetic_scale))then
      ok=.false.;message='weak overlapping-Wannier kinetic matrix is not Hermitian';return
    endif
    if(ow_diag_vlocal_hermiticity>dg_dc_gs_hermiticity_tolerance*max(1d0,local_scale))then
      ok=.false.;message='weak overlapping-Wannier local matrix is not Hermitian';return
    endif
    if(ow_diag_vnl_hermiticity>dg_dc_gs_hermiticity_tolerance*max(1d0,nonlocal_scale))then
      ok=.false.;message='overlapping-Wannier nonlocal matrix is not Hermitian';return
    endif
    if(ow_diag_h_hermiticity>dg_dc_gs_hermiticity_tolerance*max(1d0,hamiltonian_scale))then
      ok=.false.;message='weak overlapping-Wannier Hamiltonian is not Hermitian';return
    endif
    if(allocated(ow_last_kinetic_rows))deallocate(ow_last_kinetic_rows)
    if(allocated(ow_last_local_rows))deallocate(ow_last_local_rows)
    if(allocated(ow_last_nonlocal_rows))deallocate(ow_last_nonlocal_rows)
    allocate(ow_last_kinetic_rows,source=kinetic_rows)
    allocate(ow_last_local_rows,source=local_rows)
    allocate(ow_last_nonlocal_rows,source=nonlocal_rows)
    new_potential=core_potential
    fingerprint=ow_collective_operator_fingerprint(comm)
    ok=.true.;message=''
  end subroutine

  subroutine apply_dg_hybrid_divided_fragment_hpsi(tile_in,tile_out,callback_ok)
    complex(8),intent(in)::tile_in(:,:)
    complex(8),intent(out)::tile_out(:,:)
    logical,intent(out)::callback_ok
    type(s_parallel_info)::tile_info
    type(s_orbital)::tile_psi,tile_hpsi
    type(s_sendrecv_grid)::tile_srg
    integer::width,p,ix,iy,iz,io,allocation_status

    callback_ok=.false.;tile_out=(0d0,0d0);width=size(tile_in,2)
    if(size(tile_in,1)/=product(ow_box_size).or.any(shape(tile_out)/=shape(tile_in)).or.width<1)return
    if(any(mg%is_array>[1,1,1]).or.any(mg%ie_array<ow_box_size).or.system%nspin/=1)return
    tile_info=info;tile_info%im_s=1;tile_info%im_e=1;tile_info%numm=1
    tile_info%ik_s=1;tile_info%ik_e=1;tile_info%numk=1
    tile_info%io_s=1;tile_info%io_e=width;tile_info%numo=width;tile_info%if_divide_orbit=.false.
    allocate(tile_psi%zwf(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),&
      mg%is_array(3):mg%ie_array(3),1,1:width,1,1),&
      tile_hpsi%zwf(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),&
      mg%is_array(3):mg%ie_array(3),1,1:width,1,1),stat=allocation_status)
    if(allocation_status/=0)return
    tile_psi%zwf=(0d0,0d0);tile_hpsi%zwf=(0d0,0d0);p=0
    do iz=1,ow_box_size(3);do iy=1,ow_box_size(2);do ix=1,ow_box_size(1)
      p=p+1
      do io=1,width;tile_psi%zwf(ix,iy,iz,1,io,1,1)=tile_in(p,io);enddo
    enddo;enddo;enddo
    call init_sendrecv_grid(tile_srg,mg,width,info%icomm_rko,srg%neig)
    call hpsi(tile_psi,tile_hpsi,tile_info,mg,v_local,system,stencil,tile_srg,ppg)
    p=0
    do iz=1,ow_box_size(3);do iy=1,ow_box_size(2);do ix=1,ow_box_size(1)
      p=p+1
      do io=1,width;tile_out(p,io)=tile_hpsi%zwf(ix,iy,iz,1,io,1,1);enddo
    enddo;enddo;enddo
    callback_ok=all(ieee_is_finite(real(tile_out))).and.all(ieee_is_finite(aimag(tile_out)))
    call dealloc_cache(tile_srg)
    deallocate(tile_psi%zwf,tile_hpsi%zwf)
  end subroutine apply_dg_hybrid_divided_fragment_hpsi

  subroutine apply_dg_hybrid_divided_fragment_metric(tile_in,tile_out,callback_ok)
    complex(8),intent(in)::tile_in(:,:)
    complex(8),intent(out)::tile_out(:,:)
    logical,intent(out)::callback_ok
    tile_out=tile_in;callback_ok=all(ieee_is_finite(real(tile_in))).and.all(ieee_is_finite(aimag(tile_in)))
  end subroutine apply_dg_hybrid_divided_fragment_metric

  subroutine solve_dg_hybrid_divided_fragments(iteration,callback_ok)
    integer,intent(in)::iteration
    logical,intent(out)::callback_ok
    integer::fragment_basis_count,fragment_state_count,maximum_fragment_state_count,&
      fragment_point_count,ierr_local,&
      local_bad,global_bad,allocation_status
    real(8),allocatable::fragment_occupations(:),fragment_point_weights(:),fragment_core_norms(:),&
      representative_energies(:,:),representative_core_norms(:,:),all_fragment_occupations(:,:)
    real(8)::chemical_potential,common_electron_count
    logical,allocatable::representative_mask(:)
    character(256)::solver_message

    callback_ok=.false.;fragment_basis_count=size(divided_fragment_basis%global_ids)
    call MPI_Allreduce(MPI_IN_PLACE,fragment_basis_count,1,MPI_INTEGER,MPI_SUM,dc%icomm_frag,ierr_local)
    if(ierr_local/=MPI_SUCCESS)return
    fragment_state_count=min(system%no,fragment_basis_count)
    maximum_fragment_state_count=fragment_state_count
    call MPI_Allreduce(MPI_IN_PLACE,maximum_fragment_state_count,1,MPI_INTEGER,MPI_MAX,&
      dc%icomm_tot,ierr_local)
    if(ierr_local/=MPI_SUCCESS)return
    fragment_point_count=size(divided_fragment_basis%buffer_point_ids)
    local_bad=0
    if(iteration<1.or.fragment_state_count<1.or.size(system%rocc,1)<fragment_state_count)local_bad=1
    if(.not.allocated(ow_divided_core_mask))then
      local_bad=1
    elseif(size(ow_divided_core_mask)/=fragment_point_count)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr_local)
    if(ierr_local/=MPI_SUCCESS)return
    if(global_bad/=0)return
    if(allocated(divided_fragment_coefficients))deallocate(divided_fragment_coefficients)
    if(allocated(divided_fragment_eigenvalues))deallocate(divided_fragment_eigenvalues)
    if(allocated(divided_fragment_density))deallocate(divided_fragment_density)
    allocate(fragment_occupations(fragment_state_count),fragment_point_weights(fragment_point_count),&
      fragment_core_norms(fragment_state_count),divided_fragment_eigenvalues(fragment_state_count),&
      divided_fragment_density(fragment_point_count),&
      representative_energies(maximum_fragment_state_count,dc%n_frag),&
      representative_core_norms(maximum_fragment_state_count,dc%n_frag),&
      representative_mask(dc%n_frag),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr_local)
    if(ierr_local/=MPI_SUCCESS)return
    if(global_bad/=0)return
    fragment_point_weights=system%hvol
    call solve_dg_hybrid_fragment_spectrum(dc%icomm_frag,divided_fragment_basis,fragment_state_count,&
      ow_divided_core_mask,fragment_point_weights,apply_dg_hybrid_divided_fragment_hpsi,&
      apply_dg_hybrid_divided_fragment_metric,1d-12,divided_fragment_coefficients,&
      divided_fragment_eigenvalues,fragment_core_norms,divided_fragment_residual,&
      divided_fragment_orthogonality,divided_solver_workspace,divided_solver_fingerprint,&
      callback_ok,solver_message)
    if(.not.callback_ok)write(0,'(2a)')'divided fragment spectrum: ',trim(solver_message)
    local_bad=merge(0,1,callback_ok)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr_local)
    if(ierr_local/=MPI_SUCCESS)then;callback_ok=.false.;return;endif
    if(global_bad/=0)then;callback_ok=.false.;return;endif

    representative_energies=0d0;representative_core_norms=0d0;representative_mask=.false.
    if(dc%id_frag==0)then
      representative_mask(dc%i_frag)=.true.
      representative_energies(:,dc%i_frag)=divided_fragment_eigenvalues(fragment_state_count)
      representative_energies(1:fragment_state_count,dc%i_frag)=divided_fragment_eigenvalues
      representative_core_norms(1:fragment_state_count,dc%i_frag)=fragment_core_norms
    endif
    call determine_dc_fragment_occupations(dc%icomm_tot,representative_energies,&
      representative_core_norms,representative_mask,max(0d0,temperature),2d0,&
      dc%elec_num_tot,dg_dc_gs_electron_count_tolerance,chemical_potential,&
      all_fragment_occupations,common_electron_count,callback_ok,solver_message,allow_unordered=.true.)
    if(.not.callback_ok)then
      write(0,'(2a)')'divided common occupation: ',trim(solver_message);return
    endif
    if(abs(common_electron_count-dc%elec_num_tot)>dg_dc_gs_electron_count_tolerance)then
      callback_ok=.false.;write(0,'(a)')'divided common occupation electron count mismatch';return
    endif
    fragment_occupations=all_fragment_occupations(1:fragment_state_count,dc%i_frag)
    system%mu=chemical_potential
    system%rocc(:,1,1)=0d0
    system%rocc(1:fragment_state_count,1,1)=fragment_occupations
    call reconstruct_dg_hybrid_fragment_density(dc%icomm_frag,divided_fragment_basis,&
      divided_fragment_coefficients,fragment_occupations,2d0,ow_divided_core_mask,&
      fragment_point_weights,divided_fragment_density,divided_fragment_electron_count,&
      callback_ok,solver_message)
    if(.not.callback_ok)write(0,'(2a)')'divided fragment density: ',trim(solver_message)
    local_bad=merge(0,1,callback_ok)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr_local)
    if(ierr_local/=MPI_SUCCESS)then;callback_ok=.false.;return;endif
    callback_ok=global_bad==0
  end subroutine solve_dg_hybrid_divided_fragments

  subroutine gather_dg_hybrid_divided_core_density(core_density,total_density,callback_ok)
    real(8),intent(in)::core_density(:)
    real(8),allocatable,intent(out)::total_density(:,:,:)
    logical,intent(out)::callback_ok
    real(8),allocatable::local_density(:),global_density(:)
    integer::p,ierr_local

    callback_ok=.false.
    if(size(core_density)/=size(ow_core_ids).or.ow_global_grid_count<1_8)return
    allocate(local_density(int(ow_global_grid_count)),global_density(int(ow_global_grid_count)))
    local_density=0d0
    do p=1,size(ow_core_ids);local_density(int(ow_core_ids(p)))=core_density(p);enddo
    call MPI_Allreduce(local_density,global_density,size(global_density),MPI_DOUBLE_PRECISION,MPI_SUM,&
      dc%icomm_tot,ierr_local)
    if(ierr_local/=MPI_SUCCESS)return
    allocate(total_density(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3)))
    total_density=reshape(global_density,dc%lg_tot%num);callback_ok=all(ieee_is_finite(total_density))
  end subroutine gather_dg_hybrid_divided_core_density

  subroutine update_dg_hybrid_divided_potential(core_density,callback_ok)
    real(8),intent(in)::core_density(:)
    logical,intent(out)::callback_ok
    character(256)::potential_message

    call dg_dc_update_potential_from_distributed_density(ow_core_ids,core_density,callback_ok,potential_message)
    if(.not.callback_ok)write(0,'(2a)')'divided potential update: ',trim(potential_message)
  end subroutine update_dg_hybrid_divided_potential

  subroutine run_dg_hybrid_concrete_continuation(dc_seed_density,effective_ids,row_ids,basis_fragment,&
      interior_fragment,interior_weights,interior_values,interior_gradients,interior_kinetic_action,&
      interior_nonlocal_action,fixed_payload,interface_component_rows,production_faces,&
      symmetry_target_rank_arg,core_symmetry_maps_arg,basis_representation,&
      cartesian_rotations_arg,&
      requested_ids_arg,selection_effective_ids_arg,added_ids_arg,closure_parent_arg,&
      closure_reason_arg,closure_action_arg,scope_selectors_arg,pseudopotential_fingerprint_arg,&
      scope_fingerprint_arg,selection_fingerprint_arg,metric_offsets_arg,metric_columns_arg,&
      operator_offsets_arg,operator_columns_arg,&
      final_ground_state,final_density,final_trace,final_hamiltonian_rows)
    real(8),intent(in)::dc_seed_density(:),interior_weights(:)
    integer,intent(in)::effective_ids(:),basis_fragment(:),interior_fragment(:),symmetry_target_rank_arg,&
      requested_ids_arg(:),&
      selection_effective_ids_arg(:),added_ids_arg(:),&
      closure_parent_arg(:),closure_reason_arg(:),closure_action_arg(:),scope_selectors_arg(:)
    integer,intent(in)::metric_offsets_arg(:),metric_columns_arg(:),operator_offsets_arg(:),operator_columns_arg(:)
    integer(8),intent(in)::row_ids(:),core_symmetry_maps_arg(:,:)
    integer(8),intent(in)::pseudopotential_fingerprint_arg,scope_fingerprint_arg,selection_fingerprint_arg
    complex(8),intent(in)::interior_values(:,:),interior_gradients(:,:,:),&
      interior_kinetic_action(:,:),interior_nonlocal_action(:,:)
    complex(8),intent(in)::interface_component_rows(:,:,:)
    complex(8),intent(in)::basis_representation(:,:,:)
    real(8),intent(in)::cartesian_rotations_arg(:,:,:)
    type(s_dg_hybrid_fixed_payload),intent(in)::fixed_payload
    type(s_dg_hybrid_production_face_trace),intent(in)::production_faces(:)
    type(s_dg_hybrid_ground_state),intent(out)::final_ground_state
    real(8),allocatable,intent(out)::final_density(:)
    complex(8),allocatable,intent(out)::final_trace(:,:),final_hamiltonian_rows(:,:)
    real(8)::accepted_lambda,trial_lambda,electron_count,max_residual,&
      projector_change,local_norm,global_norm,&
      local_projector_scale,global_projector_scale,symmetry_residual
    real(8)::projector_symmetry_residual,low_energy_occupied_symmetry_defect,&
      occupied_symmetry_defect,target_symmetry_defect,target_energy_symmetry_defect,&
      density_symmetry_defect,input_density_symmetry_defect,output_density_symmetry_defect,&
      physical_symmetry_defect
    real(8)::broken_diagnostics(4)
    real(8),allocatable::rho_in(:),rho_out(:),local_potential(:),checkpoint_coordinates(:,:),&
      eigenvalues(:),solver_eigenvalues(:),checkpoint_cartesian_rotations(:,:,:)
    complex(8),allocatable::local_rows(:,:),coefficients(:,:),solver_coefficients(:,:),checkpoint_position(:,:,:),&
      gamma_rows(:,:),projector_rows(:,:),s_coefficients(:,:),interface_state(:,:),&
      previous_interface_state(:,:),hc(:,:),sc_epsilon(:,:),full_action_values(:,:),interface_component_actions(:,:,:),&
      full_metric(:,:),full_coefficients(:,:),full_b_rt(:,:),certified_representation(:,:,:),&
      checkpoint_basis_representation(:,:,:),checkpoint_momentum(:,:,:),&
      c_cert_rows(:,:),&
      certified_scalar_operators(:,:,:),certified_vector_operators(:,:,:,:),certified_tensor_operators(:,:,:,:,:),&
      rt_metric(:,:),rt_kinetic(:,:),rt_nonlocal(:,:),rt_local(:,:),rt_sipg(:,:),rt_hamiltonian(:,:),&
      full_component(:,:),rt_basis_values(:,:),projected_position(:,:)
    type(s_dg_hybrid_variational_iterate)::iterate
    type(s_dg_hybrid_residuals)::residuals
    type(s_dg_hybrid_controller_controls)::continuation_controls
    type(s_dg_hybrid_trial_state)::trial_state
    type(s_dg_hybrid_stage_report)::stage_report
    type(s_dg_hybrid_controller)::continuation_controller
    type(s_dg_hybrid_candidate_acceptance)::candidate_acceptance
    type(s_dg_hybrid_complete_eigensystem)::complete_eigensystem
    type(s_dg_hybrid_occupation_result)::occupation_result
    type(s_dg_hybrid_spectral_certification)::spectral_certification
    type(s_dg_hybrid_certified_rt_basis)::certified_rt_basis
    type(s_rt_dg_hybrid_ground_state_payload)::checkpoint_payload
    type(s_dg_hybrid_stage_schedule)::stage_schedule
    integer(8)::solver_fingerprint,final_operator_fingerprint,&
      final_state_workspace,final_state_fingerprint,kinetic_fingerprint,nonlocal_fingerprint,&
      local_fingerprint,sipg_fingerprint,checkpoint_fingerprint,seed_fingerprint,&
      rt_metric_fingerprint,rt_kinetic_fingerprint,rt_nonlocal_fingerprint,rt_local_fingerprint,&
      rt_sipg_fingerprint,rt_hamiltonian_fingerprint,rt_basis_fingerprint,rt_density_fingerprint,&
      rt_ownership_fingerprint,rt_grid_ownership_fingerprint,rt_scalar_fingerprint,&
      rt_vector_fingerprint,rt_tensor_fingerprint,rt_representation_fingerprint,&
      rt_rotation_fingerprint,rt_payload_fingerprint,initial_state_fingerprint,&
      catalog_ids_fingerprint,catalog_generation_fingerprint,catalog_ordering_fingerprint,&
      catalog_ownership_fingerprint,catalog_provenance_fingerprint,catalog_fingerprint,&
      symmetry_receipt_fingerprint,handoff_fingerprint
    integer::iteration,ierr_local,rank_local,gap_occupied_index,gap_unoccupied_index,p,&
      occupied_symmetry_rank,extended_target_rank,&
      face_value_count,face_value_position,face_point_count,face_point_position,face_basis_count,face_basis_position,&
      face_weight_count,face_weight_position,face_observable_count,face_observable_position,&
      interface_cursor,owned_face_count,face_slot,rt_row_count,rt_row_position,nproc_local,certified_rank,&
      construction_index,dynamic_requested_rank
    integer,allocatable::checkpoint_face_owner(:),construction_ownership(:),construction_generations(:),&
      construction_ordering(:),rt_row_owner_keys(:)
    integer(8),allocatable::construction_ids(:),rt_row_ids(:)
    logical::stage_converged,reject_trial,local_ok,accept_stage,final_refresh_performed,cheap_candidate,&
      run_solve,run_expensive,refresh_scheduled,meaningful_gap,occupation_kernel_ok,&
      seed_identity_accepted,lambda_zero_accepted,symmetry_evaluation_ok,&
      electron_gate_ok,occupied_gate_ok,density_gate_ok,payload_ready
    character(256)::continuation_message
    character(16)::window_mode
    real(8)::occupied_unoccupied_gap,accepted_gap,cluster_tolerance
    real(8)::hamiltonian_hermiticity,hamiltonian_scale
    real(8)::real_space_residual,interface_action_residuals(3),local_energy_parts(3),global_energy_parts(3),&
      final_energy_receipt(7)
    type(s_dft_energy)::checkpoint_energy
    complex(8),allocatable::energy_local_coefficients(:,:),energy_global_coefficients(:,:)
    integer::energy_row,energy_state,energy_gx,energy_gy,energy_gz,energy_ix,energy_iy,energy_iz
    logical::hamiltonian_finite

    call MPI_Comm_rank(dc%icomm_tot,rank_local,ierr_local)
    call MPI_Comm_size(dc%icomm_tot,nproc_local,ierr_local)
    if(ierr_local/=MPI_SUCCESS)error stop 'DG continuation communicator query failed'
    allocate(rho_in,source=dc_seed_density)
    seed_identity_accepted=all(rho_in==dc_seed_density);lambda_zero_accepted=.false.
    occupied_symmetry_rank=0
    if(size(core_symmetry_maps_arg,1)/=size(dc_seed_density).or.&
        size(core_symmetry_maps_arg,2)/=size(basis_representation,3))&
      error stop 'DG continuation core symmetry map shape mismatch'
    if(any(shape(cartesian_rotations_arg)/=[3,3,size(basis_representation,3)]))&
      error stop 'DG continuation Cartesian rotation shape mismatch'
    if(any(shape(interior_gradients)/=[3,size(effective_ids),size(ow_core_ids)]))&
      error stop 'DG continuation construction-gradient shape mismatch'
    meaningful_gap=.false.
    gap_occupied_index=0;gap_unoccupied_index=0
    allocate(local_potential(size(ow_core_ids)))
    allocate(previous_interface_state(0,3))
    accepted_lambda=0d0;trial_lambda=0d0;accepted_gap=huge(1d0)
    call default_dg_hybrid_controller_controls(continuation_controls)
    continuation_controls%iteration_limit=nscf
    continuation_controller%controls=continuation_controls
    allocate(full_action_values(size(effective_ids),size(ow_core_ids)))
    do while(accepted_lambda<1d0.or.trial_lambda==0d0)
      stage_converged=.false.;reject_trial=.false.;final_refresh_performed=.false.
      call initialize_dg_hybrid_stage_schedule(nscf,stage_schedule)
stage_pass: do
        call begin_dg_hybrid_stage_solve(stage_schedule,run_solve,iteration)
        if(.not.run_solve)exit stage_pass
        call dg_dc_update_potential_from_distributed_density(ow_core_ids,rho_in,local_ok,continuation_message)
        if(.not.local_ok)then;write(0,'(a)')trim(continuation_message);error stop 'DG continuation potential update failed';endif
        call extract_dg_hybrid_core_local_potential(ow_core_ids,local_potential,local_ok)
        if(.not.local_ok)error stop 'DG continuation local-potential extraction failed'
        call assemble_dg_hybrid_local_potential_rows(dc%icomm_tot,size(effective_ids),&
          row_ids,basis_fragment,ow_core_ids,interior_fragment,interior_weights,interior_values,&
          local_potential,local_rows,broken_diagnostics(1:2),local_ok,continuation_message)
        if(.not.local_ok)then;write(0,'(a)')trim(continuation_message);error stop 'DG continuation local projection failed';endif
        call compose_dg_hybrid_variational_hamiltonian(dc%icomm_tot,fixed_payload,&
          local_rows,trial_lambda,iteration,iterate,local_ok,continuation_message)
        if(.not.local_ok)then;write(0,'(a)')trim(continuation_message);error stop 'DG continuation composition failed';endif
        call initialize_dg_hybrid_candidate_acceptance(dc%icomm_tot,size(effective_ids),&
          dg_hybrid_symmetry_energy_window,candidate_acceptance,local_ok,continuation_message)
        if(.not.local_ok)error stop 'DG candidate acceptance initialization failed'
        call solve_dg_hybrid_generalized_complete_once(dc%icomm_tot,size(effective_ids),row_ids,&
          iterate%hamiltonian_rows,fixed_payload%metric_rows,dg_dc_gs_final_orbital_tolerance,&
          solve_dg_hybrid_generalized_scalapack,complete_eigensystem,local_ok,continuation_message)
        if(.not.local_ok)then
          write(0,'(a)')trim(continuation_message);error stop 'DG continuation eigensolve failed'
        endif
        solver_coefficients=complete_eigensystem%coefficients
        solver_eigenvalues=complete_eigensystem%eigenvalues
        max_residual=complete_eigensystem%maximum_residual
        solver_fingerprint=complete_eigensystem%fingerprint
        call record_dg_hybrid_complete_lcfo_solve(dc%icomm_tot,candidate_acceptance,&
          complete_eigensystem%global_count,complete_eigensystem%eigensolve_count,&
          complete_eigensystem%fingerprint,local_ok,continuation_message)
        if(.not.local_ok)error stop 'DG complete LCFO solve receipt failed'
        call derive_dg_hybrid_occupation_policy(dc%icomm_tot,complete_eigensystem%eigenvalues,&
          dc%elec_num_tot,max(0d0,temperature),dg_dc_gs_electron_count_tolerance,&
          occupation_result,local_ok,continuation_message)
        if(.not.local_ok)then
          write(0,'(a)')trim(continuation_message);error stop 'DG continuation occupation policy failed'
        endif
        electron_gate_ok=occupation_result%valid.and.&
          abs(occupation_result%electron_count-dc%elec_num_tot)<=dg_dc_gs_electron_count_tolerance
        call record_dg_hybrid_occupation_policy(dc%icomm_tot,candidate_acceptance,&
          occupation_result%noccupied,electron_gate_ok,occupation_result%fingerprint,local_ok,continuation_message)
        if(.not.local_ok)error stop 'DG candidate occupation receipt failed'
        occupied_symmetry_rank=occupation_result%noccupied
        extended_target_rank=occupation_result%noccupied
        if(dg_hybrid_symmetry_energy_window==-1d0)then
          extended_target_rank=max(occupation_result%noccupied,&
            min(symmetry_target_rank_arg,size(effective_ids)))
          do while(extended_target_rank<size(solver_eigenvalues))
            cluster_tolerance=max(dg_ow_symmetry_tolerance,64d0*epsilon(1d0))*&
              max(1d0,abs(solver_eigenvalues(extended_target_rank)),&
              abs(solver_eigenvalues(extended_target_rank+1)))
            if(solver_eigenvalues(extended_target_rank+1)-solver_eigenvalues(extended_target_rank)>&
                cluster_tolerance)exit
            extended_target_rank=extended_target_rank+1
          enddo
        endif
        if(allocated(coefficients))deallocate(coefficients)
        allocate(coefficients,source=solver_coefficients(:,1:occupation_result%noccupied))
        if(allocated(eigenvalues))deallocate(eigenvalues)
        allocate(eigenvalues,source=solver_eigenvalues(1:occupation_result%noccupied))
        meaningful_gap=occupation_result%noccupied<size(solver_eigenvalues)
        if(meaningful_gap)then
          gap_occupied_index=occupation_result%noccupied
          gap_unoccupied_index=occupation_result%noccupied+1
          occupied_unoccupied_gap=solver_eigenvalues(gap_unoccupied_index)-solver_eigenvalues(gap_occupied_index)
        else
          occupied_unoccupied_gap=huge(1d0)
        endif
        call reconstruct_dg_hybrid_occupied_state(dc%icomm_tot,size(effective_ids),row_ids,fixed_payload%metric_rows,&
          interior_values,interior_weights,coefficients,&
          occupation_result%occupations(1:occupation_result%noccupied),rho_out,gamma_rows,projector_rows,&
          s_coefficients,electron_count,local_ok,continuation_message)
        if(.not.local_ok)then;write(0,'(a)')trim(continuation_message);error stop 'DG continuation state reconstruction failed';endif
        call reconstruct_dg_hybrid_production_interface_state(dc%icomm_tot,size(effective_ids),&
          row_ids,coefficients,occupation_result%occupations(1:occupation_result%noccupied),&
          production_faces,interface_state,local_ok,continuation_message)
        if(.not.local_ok)then;write(0,'(a)')trim(continuation_message);error stop 'DG continuation trace reconstruction failed';endif
        if(size(previous_interface_state,1)==0)then
          deallocate(previous_interface_state);allocate(previous_interface_state,source=interface_state)
        endif
        call form_dg_hybrid_coefficient_actions(dc%icomm_tot,row_ids,&
          iterate%hamiltonian_rows,fixed_payload%metric_rows,coefficients,eigenvalues,hc,sc_epsilon,local_ok)
        if(.not.local_ok)error stop 'DG continuation coefficient residual action failed'
        call evaluate_dg_hybrid_residuals(dc%icomm_tot,hc,sc_epsilon,coefficients,s_coefficients,&
          rho_out,rho_in,interface_state,previous_interface_state,residuals,local_ok,continuation_message)
        if(.not.local_ok)then;write(0,'(a)')trim(continuation_message);error stop 'DG continuation residual evaluation failed';endif
        if(continuation_controller%valid)then
          local_norm=sum(abs(projector_rows-continuation_controller%accepted_state%projector)**2)
          local_projector_scale=sum(abs(continuation_controller%accepted_state%projector)**2)
          call MPI_Allreduce(local_norm,global_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr_local)
          if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(local_projector_scale,global_projector_scale,1,&
            MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr_local)
          if(ierr_local/=MPI_SUCCESS)error stop 'DG continuation projector-change reduction failed'
          projector_change=sqrt(global_norm)/max(1d0,sqrt(global_projector_scale))
        else
          projector_change=0d0
        endif
        if(continuation_controller%valid.and.continuation_controller%trial_active)then
          call observe_dg_hybrid_inner_residuals(dc%icomm_tot,continuation_controller,&
            [residuals%r_h,residuals%r_rho,residuals%r_t,residuals%r_s],reject_trial,local_ok,&
            continuation_message)
          if(.not.local_ok)error stop 'DG continuation residual-growth observation failed'
          if(reject_trial)exit
          call dg_hybrid_stage_tolerances(continuation_controls,trial_lambda,stage_report%tolerances)
        else
          stage_report%tolerances=continuation_controls%intermediate_tolerance
        endif
        cheap_candidate=all([residuals%r_h,residuals%r_rho,residuals%r_t,residuals%r_s]<=&
          stage_report%tolerances).and.abs(electron_count-dc%elec_num_tot)<=&
          dg_dc_gs_electron_count_tolerance.and.projector_change<=0.1d0
        occupation_kernel_ok=occupation_result%valid.and.&
          all(ieee_is_finite(occupation_result%occupations)).and.&
          all(occupation_result%occupations>=0d0).and.all(occupation_result%occupations<=2d0).and.&
          abs(occupation_result%electron_count-dc%elec_num_tot)<=&
          dg_dc_gs_electron_count_tolerance.and.&
          1d0-projector_change>=continuation_controls%minimum_projector_overlap
        real_space_residual=huge(1d0);interface_action_residuals=huge(1d0)
        symmetry_residual=huge(1d0);projector_symmetry_residual=huge(1d0)
        low_energy_occupied_symmetry_defect=huge(1d0);occupied_symmetry_defect=huge(1d0)
        target_symmetry_defect=huge(1d0);target_energy_symmetry_defect=huge(1d0)
        input_density_symmetry_defect=huge(1d0);output_density_symmetry_defect=huge(1d0)
        density_symmetry_defect=huge(1d0);physical_symmetry_defect=huge(1d0)
        hamiltonian_hermiticity=huge(1d0);hamiltonian_scale=1d0;hamiltonian_finite=.false.
        call schedule_dg_hybrid_candidate_checks(stage_schedule,cheap_candidate,run_expensive)
        if(run_expensive)then
          full_action_values=interior_kinetic_action+interior_nonlocal_action
          do p=1,size(ow_core_ids)
            full_action_values(:,p)=full_action_values(:,p)+local_potential(p)*interior_values(:,p)
          enddo
          call evaluate_dg_hybrid_real_space_residual(dc%icomm_tot,int(product(int(dc%lg_tot%num,8))),&
            ow_core_ids,interior_weights,size(effective_ids),row_ids,interior_values,full_action_values,&
            coefficients,eigenvalues,real_space_residual,local_ok,continuation_message)
          if(.not.local_ok)then;write(0,'(a)')trim(continuation_message);error stop 'DG strong residual evaluation failed';endif
          call reconstruct_dg_hybrid_production_interface_actions(dc%icomm_tot,size(effective_ids),row_ids,&
            coefficients,production_faces,dg_dc_gs_sipg_penalty_factor,interface_component_actions,&
            local_ok,continuation_message)
          if(.not.local_ok)then;write(0,'(a)')trim(continuation_message);error stop 'DG SIPG action reconstruction failed';endif
          call evaluate_dg_hybrid_face_action_residuals(dc%icomm_tot,size(effective_ids),row_ids,&
            interface_component_rows,iterate%hamiltonian_rows,coefficients,sc_epsilon,trial_lambda,&
            interface_component_actions,interface_action_residuals,&
            local_ok,continuation_message)
          if(.not.local_ok)then;write(0,'(a)')trim(continuation_message);error stop 'DG SIPG action residual evaluation failed';endif
          if(size(basis_representation,3)==0)then
            symmetry_residual=0d0;projector_symmetry_residual=0d0
            low_energy_occupied_symmetry_defect=0d0;target_symmetry_defect=0d0
            target_energy_symmetry_defect=0d0;occupied_symmetry_defect=0d0
            input_density_symmetry_defect=0d0;output_density_symmetry_defect=0d0
            density_symmetry_defect=0d0;symmetry_evaluation_ok=.true.
          else
            call measure_dg_hybrid_operator_covariance(dc%icomm_tot,row_ids,iterate%hamiltonian_rows,&
              basis_representation,symmetry_residual,local_ok)
            if(.not.local_ok)error stop 'DG continuation Hamiltonian covariance measurement failed'
            call measure_dg_hybrid_projector_covariance(dc%icomm_tot,row_ids,projector_rows,&
              basis_representation,projector_symmetry_residual,local_ok)
            if(.not.local_ok)error stop 'DG continuation occupied-projector covariance measurement failed'
            call evaluate_dg_hybrid_distributed_low_energy_symmetry(dc%icomm_tot,row_ids,&
              fixed_payload%metric_rows,basis_representation,solver_coefficients,solver_eigenvalues,&
              occupied_symmetry_rank,extended_target_rank,dg_ow_symmetry_tolerance,&
              low_energy_occupied_symmetry_defect,target_symmetry_defect,target_energy_symmetry_defect,&
              symmetry_evaluation_ok,continuation_message)
            if(.not.symmetry_evaluation_ok)then
              write(0,'(a)')trim(continuation_message)
              error stop 'DG continuation low-energy symmetry evaluation failed'
            endif
            occupied_symmetry_defect=max(low_energy_occupied_symmetry_defect,projector_symmetry_residual)
            call measure_dg_hybrid_core_density_covariance(dc%icomm_tot,rho_out,core_symmetry_maps_arg,&
              output_density_symmetry_defect,local_ok,continuation_message)
            if(.not.local_ok)then
              write(0,'(a)')trim(continuation_message);error stop 'DG continuation output-density symmetry failed'
            endif
            call measure_dg_hybrid_core_density_covariance(dc%icomm_tot,rho_in,core_symmetry_maps_arg,&
              input_density_symmetry_defect,local_ok,continuation_message)
            if(.not.local_ok)then
              write(0,'(a)')trim(continuation_message);error stop 'DG continuation input-density symmetry failed'
            endif
            density_symmetry_defect=max(input_density_symmetry_defect,output_density_symmetry_defect)
          endif
          physical_symmetry_defect=max(occupied_symmetry_defect,density_symmetry_defect)
          if(dg_hybrid_symmetry_energy_window==-1d0)physical_symmetry_defect=&
            max(physical_symmetry_defect,target_symmetry_defect,target_energy_symmetry_defect)
          call ow_distributed_hermiticity(dc%icomm_tot,row_ids,iterate%hamiltonian_rows,&
            hamiltonian_hermiticity,hamiltonian_scale,hamiltonian_finite)
        endif
        stage_converged=cheap_candidate.and.real_space_residual<=stage_report%tolerances(1).and.&
          all(interface_action_residuals<=stage_report%tolerances(1)).and.&
          physical_symmetry_defect<=dg_ow_symmetry_tolerance.and.&
          occupation_kernel_ok.and.hamiltonian_finite.and.&
          hamiltonian_hermiticity<=dg_dc_gs_hermiticity_tolerance*max(1d0,hamiltonian_scale)
        call complete_dg_hybrid_stage_solve(stage_schedule,stage_converged,trial_lambda,&
          refresh_scheduled,final_refresh_performed)
        if(refresh_scheduled)then
          rho_in=rho_out;previous_interface_state=interface_state;stage_converged=.false.
          cycle stage_pass
        endif
        if(stage_converged)exit stage_pass
        rho_in=rho_in+continuation_controller%controls%density_damping*(rho_out-rho_in)
        previous_interface_state=interface_state
      enddo stage_pass
      if(.not.stage_converged)reject_trial=.true.
      if(reject_trial)then
        call reject_dg_hybrid_trial(dc%icomm_tot,continuation_controller,trial_state,&
          'inner fixed point did not converge',local_ok,continuation_message)
        if(.not.local_ok)then;write(0,'(a)')trim(continuation_message);error stop 'DG continuation rollback failed';endif
        rho_in=trial_state%density
        previous_interface_state=trial_state%trace
        call propose_dg_hybrid_trial(dc%icomm_tot,continuation_controller,trial_state,local_ok,continuation_message)
        if(.not.local_ok)error stop 'DG continuation retry proposal failed'
        trial_lambda=continuation_controller%trial_lambda
        cycle
      endif
      previous_interface_state=interface_state
      call set_dg_hybrid_trial_state(trial_state,rho_in,local_potential,&
        occupation_result%occupations(1:occupation_result%noccupied),eigenvalues,&
        projector_rows,interface_state,iteration,fixed_payload%fingerprint,solver_fingerprint)
      if(.not.continuation_controller%valid)then
        call initialize_dg_hybrid_controller(dc%icomm_tot,continuation_controls,0d0,trial_state,&
          size(production_faces),continuation_controller,local_ok,continuation_message)
        if(.not.local_ok)error stop 'DG continuation controller initialization failed'
        accepted_lambda=0d0;lambda_zero_accepted=.true.
      else
        stage_report%residuals=[residuals%r_h,residuals%r_rho,residuals%r_t,residuals%r_s]
        stage_report%projector_overlap=max(0d0,1d0-projector_change)
        stage_report%iteration=iteration
        stage_report%electron_ok=abs(electron_count-dc%elec_num_tot)<=dg_dc_gs_electron_count_tolerance
        stage_report%occupation_ok=occupation_kernel_ok
        stage_report%hermitian_ok=hamiltonian_finite.and.&
          hamiltonian_hermiticity<=dg_dc_gs_hermiticity_tolerance*max(1d0,hamiltonian_scale)
        stage_report%symmetry_ok=physical_symmetry_defect<=dg_ow_symmetry_tolerance
        stage_report%real_space_ok=real_space_residual<=stage_report%tolerances(1).and.&
          all(interface_action_residuals<=stage_report%tolerances(1)).and.&
          residuals%r_t<=stage_report%tolerances(3)
        stage_report%finite_ok=hamiltonian_finite.and.all(ieee_is_finite(solver_eigenvalues)).and.&
          all(ieee_is_finite([occupied_symmetry_defect,target_symmetry_defect,target_energy_symmetry_defect,&
          density_symmetry_defect])).and.ieee_is_finite(electron_count).and.&
          (.not.meaningful_gap.or.ieee_is_finite(occupied_unoccupied_gap))
        stage_report%gap_shrinking=meaningful_gap.and.accepted_gap<huge(1d0).and.&
          occupied_unoccupied_gap<accepted_gap
        call decide_dg_hybrid_stage(dc%icomm_tot,continuation_controller,trial_state,stage_report,&
          accept_stage,local_ok,continuation_message)
        if(.not.local_ok.or..not.accept_stage)error stop 'DG continuation converged stage was not accepted'
        accepted_lambda=continuation_controller%accepted_lambda
        if(accepted_lambda==0d0)lambda_zero_accepted=.true.
      endif
      if(meaningful_gap)accepted_gap=occupied_unoccupied_gap
      if(accepted_lambda==1d0)exit
      call propose_dg_hybrid_trial(dc%icomm_tot,continuation_controller,trial_state,local_ok,continuation_message)
      if(.not.local_ok)error stop 'DG continuation trial proposal failed'
      rho_in=trial_state%density;previous_interface_state=trial_state%trace
      trial_lambda=continuation_controller%trial_lambda
    enddo
    if(.not.final_refresh_performed)error stop 'DG continuation lambda-one state was not fully refreshed'
    call ow_fingerprint_distributed_matrix(dc%icomm_tot,row_ids,iterate%hamiltonian_rows,&
      final_operator_fingerprint,local_ok)
    if(.not.local_ok)error stop 'DG continuation final Hamiltonian fingerprint failed'

    ! These physical gates are unconditional for the final candidate.  The
    ! density below is the authoritative output reconstructed with Task-9
    ! occupations, never a symmetrized replacement.
    rho_in=rho_out
    if(size(basis_representation,3)==0)then
      projector_symmetry_residual=0d0;density_symmetry_defect=0d0
    else
      call measure_dg_hybrid_projector_covariance(dc%icomm_tot,row_ids,projector_rows,&
        basis_representation,projector_symmetry_residual,local_ok)
      if(.not.local_ok)error stop 'DG final occupied-projector covariance measurement failed'
      call measure_dg_hybrid_core_density_covariance(dc%icomm_tot,rho_in,core_symmetry_maps_arg,&
        density_symmetry_defect,local_ok,continuation_message)
      if(.not.local_ok)error stop 'DG final density covariance measurement failed'
    endif
    occupied_gate_ok=projector_symmetry_residual-dg_ow_symmetry_tolerance<=0d0
    density_gate_ok=density_symmetry_defect<=dg_ow_symmetry_tolerance
    call record_dg_hybrid_unconditional_gates(dc%icomm_tot,candidate_acceptance,&
      occupied_gate_ok,density_gate_ok,projector_symmetry_residual,density_symmetry_defect,&
      dg_ow_symmetry_tolerance,local_ok,continuation_message)
    if(.not.local_ok)then
      write(0,'(a)')trim(continuation_message)
      error stop 'DG final occupied-projector or density gate failed'
    endif

    call collect_dg_hybrid_full_rows(dc%icomm_tot,size(effective_ids),row_ids,&
      fixed_payload%metric_rows,full_metric,local_ok)
    if(.not.local_ok)error stop 'DG final metric collection failed'
    call collect_dg_hybrid_full_rows(dc%icomm_tot,size(effective_ids),row_ids,&
      complete_eigensystem%coefficients,full_coefficients,local_ok)
    if(.not.local_ok)error stop 'DG complete eigenvector collection failed'
    ! The analysis routines use a generator set, but the checkpoint operation
    ! catalog is a complete advertised set and therefore carries identity
    ! explicitly as operation one.
    allocate(checkpoint_basis_representation(size(effective_ids),size(effective_ids),&
      size(basis_representation,3)+1))
    checkpoint_basis_representation=(0d0,0d0)
    do p=1,size(effective_ids)
      checkpoint_basis_representation(p,p,1)=(1d0,0d0)
    enddo
    if(size(basis_representation,3)>0)&
      checkpoint_basis_representation(:,:,2:)=basis_representation
    allocate(checkpoint_cartesian_rotations(3,3,size(cartesian_rotations_arg,3)+1))
    checkpoint_cartesian_rotations=0d0
    do p=1,3
      checkpoint_cartesian_rotations(p,p,1)=1d0
    enddo
    ! The production point action satisfies the covector convention measured
    ! by measure_dg_spatial_gradient_covariance: R^T grad(g r)=D^T grad(r).
    ! Task 11 stores the component matrix in D^H X_a D=sum_b Q_ab X_b,
    ! hence canonical momentum uses Q=R^T.
    do p=1,size(cartesian_rotations_arg,3)
      checkpoint_cartesian_rotations(:,:,p+1)=transpose(cartesian_rotations_arg(:,:,p))
    enddo
    dynamic_requested_rank=occupation_result%noccupied
    if(dg_hybrid_symmetry_energy_window==-1d0)dynamic_requested_rank=&
      max(occupation_result%noccupied,min(symmetry_target_rank_arg,size(effective_ids)))
    call certify_dg_hybrid_energy_window(dc%icomm_tot,full_metric,checkpoint_basis_representation,&
      full_coefficients,complete_eigensystem%eigenvalues,occupation_result%noccupied,&
      occupation_result%e_homo,dg_hybrid_symmetry_energy_window,&
      dynamic_requested_rank,&
      dg_dc_gs_final_orbital_tolerance,dg_ow_symmetry_tolerance,projector_symmetry_residual,&
      density_symmetry_defect,spectral_certification,local_ok,continuation_message)
    if(.not.local_ok)then
      write(0,'(a)')trim(continuation_message)
      error stop 'DG final LCFO energy-window certification failed'
    endif
    extended_target_rank=spectral_certification%certified_rank
    occupied_symmetry_defect=max(spectral_certification%occupied_subspace_defect,&
      projector_symmetry_residual)
    target_symmetry_defect=spectral_certification%target_subspace_defect
    target_energy_symmetry_defect=spectral_certification%target_energy_defect
    physical_symmetry_defect=max(occupied_symmetry_defect,target_symmetry_defect,&
      target_energy_symmetry_defect,density_symmetry_defect)
    call record_dg_hybrid_spectral_certification(dc%icomm_tot,candidate_acceptance,&
      spectral_certification%requested_rank,spectral_certification%certified_rank,&
      spectral_certification%boundary_cluster_rank,spectral_certification%proof_state_present,&
      spectral_certification%compatibility_dynamic_rank,dg_hybrid_symmetry_energy_window==-1d0,&
      spectral_certification%fingerprint,local_ok,continuation_message)
    if(.not.local_ok)error stop 'DG spectral certification receipt failed'

    certified_rank=spectral_certification%certified_rank
    allocate(checkpoint_coordinates(3,size(ow_core_ids)))
    do p=1,size(ow_core_ids)
      checkpoint_coordinates(1,p)=real(modulo(ow_core_ids(p)-1_8,int(dc%lg_tot%num(1),8)),8)*dc%system_tot%hgs(1)
      checkpoint_coordinates(2,p)=real(modulo((ow_core_ids(p)-1_8)/int(dc%lg_tot%num(1),8),&
        int(dc%lg_tot%num(2),8)),8)*dc%system_tot%hgs(2)
      checkpoint_coordinates(3,p)=real((ow_core_ids(p)-1_8)/int(dc%lg_tot%num(1)*dc%lg_tot%num(2),8),8)*&
        dc%system_tot%hgs(3)
    enddo
    call assemble_dg_cell_wrapped_position(dc%icomm_tot,ow_core_ids,interior_weights,checkpoint_coordinates,&
      [0d0,0d0,0d0],real(dc%lg_tot%num,8)*dc%system_tot%hgs,interior_values,checkpoint_position,&
      checkpoint_payload%position_convention_fingerprint,local_ok,continuation_message)
    if(.not.local_ok)error stop 'DG continuation periodic position assembly failed'
    allocate(c_cert_rows,source=complete_eigensystem%coefficients(:,1:certified_rank))
    call build_dg_hybrid_certified_representation(full_metric,checkpoint_basis_representation,&
      full_coefficients(:,1:certified_rank),certified_representation)
    allocate(certified_scalar_operators(certified_rank,certified_rank,5),&
      certified_vector_operators(certified_rank,certified_rank,3,&
        rt_dg_hybrid_vector_canonical_momentum),&
      certified_tensor_operators(certified_rank,certified_rank,3,3,0))
    call collect_dg_hybrid_full_rows(dc%icomm_tot,size(effective_ids),row_ids,&
      fixed_payload%kinetic_rows,full_component,local_ok)
    if(.not.local_ok)error stop 'DG certified kinetic source collection failed'
    call project_dg_hybrid_operator_to_rt(full_coefficients(:,1:certified_rank),full_component,projected_position)
    certified_scalar_operators(:,:,1)=projected_position
    call collect_dg_hybrid_full_rows(dc%icomm_tot,size(effective_ids),row_ids,&
      fixed_payload%nonlocal_rows,full_component,local_ok)
    if(.not.local_ok)error stop 'DG certified nonlocal source collection failed'
    call project_dg_hybrid_operator_to_rt(full_coefficients(:,1:certified_rank),full_component,projected_position)
    certified_scalar_operators(:,:,2)=projected_position
    call collect_dg_hybrid_full_rows(dc%icomm_tot,size(effective_ids),row_ids,&
      iterate%local_rows,full_component,local_ok)
    if(.not.local_ok)error stop 'DG certified local source collection failed'
    call project_dg_hybrid_operator_to_rt(full_coefficients(:,1:certified_rank),full_component,projected_position)
    certified_scalar_operators(:,:,3)=projected_position
    call collect_dg_hybrid_full_rows(dc%icomm_tot,size(effective_ids),row_ids,&
      fixed_payload%interface_rows,full_component,local_ok)
    if(.not.local_ok)error stop 'DG certified SIPG source collection failed'
    call project_dg_hybrid_operator_to_rt(full_coefficients(:,1:certified_rank),full_component,projected_position)
    certified_scalar_operators(:,:,4)=projected_position
    call collect_dg_hybrid_full_rows(dc%icomm_tot,size(effective_ids),row_ids,&
      iterate%hamiltonian_rows,full_component,local_ok)
    if(.not.local_ok)error stop 'DG certified Hamiltonian source collection failed'
    call project_dg_hybrid_operator_to_rt(full_coefficients(:,1:certified_rank),full_component,projected_position)
    certified_scalar_operators(:,:,5)=projected_position
    ! Vector item one is the canonical-momentum (velocity-gauge) coupling.
    ! Unlike cell-wrapped position it obeys a homogeneous crystallographic
    ! vector law even for operations with an affine translation.
    call assemble_dg_hybrid_canonical_momentum(dc%icomm_tot,interior_weights,interior_values,&
      interior_gradients,checkpoint_momentum,local_ok,continuation_message)
    if(.not.local_ok)error stop 'DG certified canonical-momentum assembly failed'
    do p=1,3
      call project_dg_hybrid_operator_to_rt(full_coefficients(:,1:certified_rank),&
        checkpoint_momentum(p,:,:),projected_position)
      certified_vector_operators(:,:,p,rt_dg_hybrid_vector_canonical_momentum)=projected_position
    enddo
    if(allocated(dg_certified_localizer_values))deallocate(dg_certified_localizer_values)
    if(allocated(dg_certified_localizer_weights))deallocate(dg_certified_localizer_weights)
    if(allocated(dg_certified_localizer_grid_ids))deallocate(dg_certified_localizer_grid_ids)
    allocate(dg_certified_localizer_values,source=interior_values)
    allocate(dg_certified_localizer_weights,source=interior_weights)
    allocate(dg_certified_localizer_grid_ids,source=ow_core_ids)
    call build_dg_hybrid_certified_rt_basis(dc%icomm_tot,size(effective_ids),row_ids,&
      fixed_payload%metric_rows,c_cert_rows,&
      complete_eigensystem%eigenvalues(1:certified_rank),occupation_result%noccupied,&
      certified_representation,checkpoint_cartesian_rotations,certified_scalar_operators,&
      certified_vector_operators,certified_tensor_operators,dg_ow_symmetry_tolerance,&
      dg_hybrid_certified_localizer_candidate,certified_rt_basis,local_ok,continuation_message)
    deallocate(dg_certified_localizer_values,dg_certified_localizer_weights,dg_certified_localizer_grid_ids)
    if(.not.local_ok)then
      write(0,'(a)')trim(continuation_message)
      error stop 'DG certified RT basis construction failed'
    endif
    call record_dg_hybrid_certified_rt_basis(dc%icomm_tot,candidate_acceptance,&
      certified_rt_basis%certified_rank,certified_rt_basis%b_rt_fingerprint,&
      certified_rt_basis%operator_fingerprint,local_ok,continuation_message)
    if(.not.local_ok)error stop 'DG certified RT basis receipt failed'

    call validate_dg_hybrid_ground_state(dc%icomm_tot,size(effective_ids),occupation_result%noccupied,&
      row_ids,coefficients,occupation_result%occupations(1:occupation_result%noccupied),eigenvalues,&
      dc%elec_num_tot,fixed_payload%basis_fingerprint,&
      fixed_payload%metric_fingerprint,final_operator_fingerprint,&
      checkpoint_payload%position_convention_fingerprint,&
      dg_dc_gs_final_orbital_tolerance,final_ground_state,final_state_workspace,final_state_fingerprint,&
      local_ok,continuation_message)
    if(.not.local_ok)then;write(0,'(a)')trim(continuation_message);error stop 'DG continuation final state publication failed';endif
    final_ground_state%converged=.true.;final_ground_state%final_eigensolve_count=1
    final_ground_state%spectral_certification=spectral_certification
    allocate(final_density,source=rho_in);allocate(final_trace,source=interface_state)
    allocate(final_hamiltonian_rows,source=iterate%hamiltonian_rows)
    call fingerprint_rt_dg_hybrid_component(dc%icomm_tot,row_ids,fixed_payload%kinetic_rows,kinetic_fingerprint,local_ok)
    if(.not.local_ok)error stop 'DG continuation kinetic checkpoint fingerprint failed'
    call fingerprint_rt_dg_hybrid_component(dc%icomm_tot,row_ids,fixed_payload%nonlocal_rows,nonlocal_fingerprint,local_ok)
    if(.not.local_ok)error stop 'DG continuation nonlocal checkpoint fingerprint failed'
    call fingerprint_rt_dg_hybrid_component(dc%icomm_tot,row_ids,iterate%local_rows,local_fingerprint,local_ok)
    if(.not.local_ok)error stop 'DG continuation local checkpoint fingerprint failed'
    call fingerprint_rt_dg_hybrid_component(dc%icomm_tot,row_ids,fixed_payload%interface_rows,sipg_fingerprint,local_ok)
    if(.not.local_ok)error stop 'DG continuation SIPG checkpoint fingerprint failed'
    call checkpoint_grid_real_fingerprint(dc%icomm_tot,ow_core_ids,dc_seed_density,&
      seed_fingerprint,local_ok)
    if(.not.local_ok)error stop 'DG continuation seed checkpoint fingerprint failed'
    checkpoint_payload%valid=.true.;checkpoint_payload%final_refresh_complete=final_refresh_performed
    checkpoint_payload%analysis_complete=.true.;checkpoint_payload%identity_only=size(basis_representation,3)==0
    checkpoint_payload%operation_count=size(checkpoint_basis_representation,3)
    checkpoint_payload%nonidentity_operation_count=size(basis_representation,3)
    checkpoint_payload%global_count=size(effective_ids);checkpoint_payload%noccupied=occupation_result%noccupied
    checkpoint_payload%state_fingerprint=final_state_fingerprint
    checkpoint_payload%metric_fingerprint=fixed_payload%metric_fingerprint
    checkpoint_payload%operator_structure_fingerprint=fixed_payload%fingerprint
    checkpoint_payload%operator_value_fingerprint=final_operator_fingerprint
    checkpoint_payload%kinetic_fingerprint=kinetic_fingerprint
    checkpoint_payload%nonlocal_fingerprint=nonlocal_fingerprint
    checkpoint_payload%local_fingerprint=local_fingerprint
    checkpoint_payload%sipg_fingerprint=sipg_fingerprint
    checkpoint_payload%basis_fingerprint=fixed_payload%basis_fingerprint
    checkpoint_payload%face_fingerprint=fixed_payload%interface_fingerprint
    checkpoint_payload%dc_seed_fingerprint=seed_fingerprint
    checkpoint_payload%continuation_fingerprint=ieor(final_state_fingerprint,continuation_controller%accepted_state%operator_value_fingerprint)
    if(checkpoint_payload%continuation_fingerprint==0_8)checkpoint_payload%continuation_fingerprint=1_8
    checkpoint_payload%scope_fingerprint=scope_fingerprint_arg
    checkpoint_payload%selection_fingerprint=selection_fingerprint_arg
    checkpoint_payload%analysis_fingerprint=checkpoint_analysis_fingerprint(checkpoint_basis_representation)
    checkpoint_payload%pseudopotential_fingerprint=pseudopotential_fingerprint_arg
    allocate(checkpoint_payload%row_ids,source=row_ids)
    allocate(checkpoint_payload%metric_row_offsets,source=metric_offsets_arg)
    allocate(checkpoint_payload%metric_column_ids,source=metric_columns_arg)
    allocate(checkpoint_payload%operator_row_offsets,source=operator_offsets_arg)
    allocate(checkpoint_payload%operator_column_ids,source=operator_columns_arg)
    allocate(checkpoint_payload%metric_rows,source=fixed_payload%metric_rows)
    allocate(checkpoint_payload%kinetic_rows,source=fixed_payload%kinetic_rows)
    allocate(checkpoint_payload%nonlocal_rows,source=fixed_payload%nonlocal_rows)
    allocate(checkpoint_payload%local_rows,source=iterate%local_rows)
    allocate(checkpoint_payload%sipg_rows,source=fixed_payload%interface_rows)
    allocate(checkpoint_payload%hamiltonian_rows,source=iterate%hamiltonian_rows)
    allocate(checkpoint_payload%coefficients,source=final_ground_state%coefficients)
    allocate(checkpoint_payload%symmetry_representation,source=checkpoint_basis_representation)
    allocate(checkpoint_payload%position_rows(3,size(row_ids),size(checkpoint_position,3)))
    do i=1,size(row_ids)
      checkpoint_payload%position_rows(:,i,:)=checkpoint_position(:,int(row_ids(i)),:)
    enddo
    checkpoint_payload%global_grid_count=product(dc%lg_tot%num)
    allocate(checkpoint_payload%occupations,source=final_ground_state%occupations)
    allocate(checkpoint_payload%eigenvalues,source=final_ground_state%eigenvalues)
    allocate(checkpoint_payload%grid_ids,source=ow_core_ids)
    allocate(checkpoint_payload%grid_weights,source=interior_weights)
    allocate(checkpoint_payload%partition_ids,source=interior_fragment)
    allocate(checkpoint_payload%basis_values,source=interior_values)
    allocate(checkpoint_payload%density,source=rho_in)
    allocate(checkpoint_payload%requested_ids,source=requested_ids_arg)
    allocate(checkpoint_payload%effective_ids,source=selection_effective_ids_arg)
    allocate(checkpoint_payload%added_ids,source=added_ids_arg)
    allocate(checkpoint_payload%closure_parent,source=closure_parent_arg)
    allocate(checkpoint_payload%closure_reason,source=closure_reason_arg)
    allocate(checkpoint_payload%closure_action,source=closure_action_arg)
    allocate(checkpoint_payload%scope_selectors(size(scope_selectors_arg)),checkpoint_payload%xc_types(size(xc_func%xctype)))
    checkpoint_payload%scope_selectors=scope_selectors_arg
    checkpoint_payload%xc_types=xc_func%xctype
    allocate(checkpoint_payload%continuation_receipt(16))
    checkpoint_payload%continuation_receipt=[accepted_lambda,residuals%r_h,residuals%r_rho,residuals%r_t,residuals%r_s,&
      electron_count,real(symmetry_target_rank_arg,8),real(extended_target_rank,8),symmetry_residual,&
      occupied_symmetry_defect,target_symmetry_defect,target_energy_symmetry_defect,density_symmetry_defect,&
      projector_symmetry_residual,real_space_residual,maxval(interface_action_residuals)]
    allocate(checkpoint_payload%pseudopotential_receipt(6),checkpoint_payload%energy_receipt(7))
    checkpoint_payload%pseudopotential_receipt=[real(dc%system_tot%nion,8),pp%zion,real(pp%lmax,8),&
      real(pp%nrmax,8),real(ppg%Nlma,8),real(size(effective_ids),8)*real(size(effective_ids),8)]
    allocate(energy_local_coefficients(size(effective_ids),occupation_result%noccupied),&
      energy_global_coefficients(size(effective_ids),occupation_result%noccupied))
    energy_local_coefficients=(0d0,0d0)
    do energy_row=1,size(row_ids)
      energy_local_coefficients(int(row_ids(energy_row)),:)=final_ground_state%coefficients(energy_row,:)
    enddo
    call MPI_Allreduce(energy_local_coefficients,energy_global_coefficients,size(energy_global_coefficients),&
      MPI_DOUBLE_COMPLEX,MPI_SUM,dc%icomm_tot,ierr_local)
    if(ierr_local/=MPI_SUCCESS)error stop 'DG continuation energy coefficient redistribution failed'
    local_energy_parts=0d0
    do energy_row=1,size(row_ids);do energy_state=1,occupation_result%noccupied
      local_energy_parts(1)=local_energy_parts(1)+final_ground_state%occupations(energy_state)*real(&
        conjg(final_ground_state%coefficients(energy_row,energy_state))*&
        sum((fixed_payload%kinetic_rows(energy_row,:)+fixed_payload%interface_rows(energy_row,:))*&
        energy_global_coefficients(:,energy_state)))
      local_energy_parts(2)=local_energy_parts(2)+final_ground_state%occupations(energy_state)*real(&
        conjg(final_ground_state%coefficients(energy_row,energy_state))*&
        sum(fixed_payload%nonlocal_rows(energy_row,:)*energy_global_coefficients(:,energy_state)))
    enddo;enddo
    do p=1,size(ow_core_ids)
      energy_gx=int(modulo(ow_core_ids(p)-1_8,int(dc%lg_tot%num(1),8)))+1
      energy_gy=int(modulo((ow_core_ids(p)-1_8)/int(dc%lg_tot%num(1),8),int(dc%lg_tot%num(2),8)))+1
      energy_gz=int((ow_core_ids(p)-1_8)/int(dc%lg_tot%num(1)*dc%lg_tot%num(2),8))+1
      energy_ix=findloc(dc%jxyz_tot(:,1),energy_gx,dim=1)
      energy_iy=findloc(dc%jxyz_tot(:,2),energy_gy,dim=1)
      energy_iz=findloc(dc%jxyz_tot(:,3),energy_gz,dim=1)
      if(energy_ix<1.or.energy_iy<1.or.energy_iz<1)error stop 'DG continuation energy grid mapping failed'
      local_energy_parts(3)=local_energy_parts(3)+eexc_tmp(energy_ix,energy_iy,energy_iz)*interior_weights(p)
    enddo
    call MPI_Allreduce(local_energy_parts,global_energy_parts,3,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.any(.not.ieee_is_finite(global_energy_parts)))&
      error stop 'DG continuation energy decomposition failed'
    checkpoint_energy%E_kin=global_energy_parts(1)
    checkpoint_energy%E_ion_nloc=global_energy_parts(2)
    checkpoint_energy%E_xc=global_energy_parts(3)
    call calc_Total_Energy_periodic(dc%mg_tot,ewald,dc%system_tot,dc%info_tot,pp,dc%ppg_tot,&
      dc%fg_tot,dc%poisson_tot,.true.,checkpoint_energy)
    final_energy_receipt=[checkpoint_energy%E_tot,checkpoint_energy%E_kin,checkpoint_energy%E_h,&
      checkpoint_energy%E_xc,checkpoint_energy%E_ion_ion,checkpoint_energy%E_ion_loc,checkpoint_energy%E_ion_nloc]
    if(any(.not.ieee_is_finite(final_energy_receipt)).or.&
        .not.(abs(final_energy_receipt(1)-sum(final_energy_receipt(2:7)))<=&
        100d0*epsilon(1d0)*max(1d0,abs(final_energy_receipt(1)))))&
      error stop 'DG continuation final energy receipt is inconsistent'
    checkpoint_payload%energy_receipt=final_energy_receipt
    checkpoint_payload%energy_fingerprint=checkpoint_real_fingerprint(checkpoint_payload%energy_receipt)
    allocate(checkpoint_payload%nonlocal_ids,source=ow_core_ids)
    allocate(checkpoint_payload%nonlocal_owner(size(ow_core_ids)),source=rank_local)
    allocate(checkpoint_payload%nonlocal_values,source=interior_nonlocal_action)
    allocate(checkpoint_face_owner(size(production_faces)))
    do p=1,size(production_faces)
      checkpoint_face_owner(p)=merge(rank_local,huge(0),production_faces(p)%frozen)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,checkpoint_face_owner,size(checkpoint_face_owner),MPI_INTEGER,MPI_MIN,&
      dc%icomm_tot,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.any(checkpoint_face_owner==huge(0)))&
      error stop 'DG continuation checkpoint face ownership failed'
    owned_face_count=count(checkpoint_face_owner==rank_local)
    face_point_count=0;face_value_count=0;face_basis_count=0;face_weight_count=0;face_observable_count=0
    do p=1,size(production_faces)
      if(checkpoint_face_owner(p)/=rank_local)cycle
      face_point_count=face_point_count+size(production_faces(p)%point_ids_minus)+size(production_faces(p)%point_ids_plus)
      face_value_count=face_value_count+size(production_faces(p)%value_minus)+size(production_faces(p)%value_plus)+&
        size(production_faces(p)%derivative_minus)+size(production_faces(p)%derivative_plus)
      face_basis_count=face_basis_count+size(production_faces(p)%basis_ids_minus)+size(production_faces(p)%basis_ids_plus)
      face_weight_count=face_weight_count+size(production_faces(p)%weights)
      face_observable_count=face_observable_count+size(production_faces(p)%weights)**2
    enddo
    allocate(checkpoint_payload%face_ids(owned_face_count),checkpoint_payload%face_metadata(8,owned_face_count),&
      checkpoint_payload%face_normals(3,owned_face_count),checkpoint_payload%face_offsets(owned_face_count+1),&
      checkpoint_payload%face_weight_offsets(owned_face_count+1),&
      checkpoint_payload%face_basis_offsets(owned_face_count+1),&
      checkpoint_payload%face_value_offsets(owned_face_count+1),&
      checkpoint_payload%face_observable_offsets(owned_face_count+1),&
      checkpoint_payload%face_basis_ids(face_basis_count),&
      checkpoint_payload%face_point_ids(face_point_count),checkpoint_payload%face_weights(face_weight_count),&
      checkpoint_payload%face_values(1,face_value_count),&
      checkpoint_payload%interface_observables(3,face_observable_count))
    face_point_position=0;face_value_position=0;face_basis_position=0;face_weight_position=0
    face_observable_position=0;interface_cursor=0;face_slot=0
    checkpoint_payload%face_offsets(1)=1;checkpoint_payload%face_weight_offsets(1)=1
    checkpoint_payload%face_basis_offsets(1)=1;checkpoint_payload%face_value_offsets(1)=1
    checkpoint_payload%face_observable_offsets(1)=1
    do p=1,size(production_faces)
      if(.not.production_faces(p)%frozen)cycle
      if(checkpoint_face_owner(p)/=rank_local)then
        interface_cursor=interface_cursor+size(production_faces(p)%weights)**2
        cycle
      endif
      face_slot=face_slot+1;checkpoint_payload%face_ids(face_slot)=production_faces(p)%global_face_id
      checkpoint_payload%face_metadata(:,face_slot)=[production_faces(p)%minus_fragment,production_faces(p)%plus_fragment,&
        production_faces(p)%periodic_shift,size(production_faces(p)%weights),size(production_faces(p)%basis_ids_minus),&
        size(production_faces(p)%basis_ids_plus)]
      checkpoint_payload%face_normals(:,face_slot)=production_faces(p)%canonical_normal
      checkpoint_payload%face_point_ids(face_point_position+1:face_point_position+size(production_faces(p)%point_ids_minus))=&
        production_faces(p)%point_ids_minus;face_point_position=face_point_position+size(production_faces(p)%point_ids_minus)
      checkpoint_payload%face_point_ids(face_point_position+1:face_point_position+size(production_faces(p)%point_ids_plus))=&
        production_faces(p)%point_ids_plus;face_point_position=face_point_position+size(production_faces(p)%point_ids_plus)
      checkpoint_payload%face_offsets(face_slot+1)=face_point_position+1
      checkpoint_payload%face_basis_ids(face_basis_position+1:face_basis_position+size(production_faces(p)%basis_ids_minus))=&
        production_faces(p)%basis_ids_minus;face_basis_position=face_basis_position+size(production_faces(p)%basis_ids_minus)
      checkpoint_payload%face_basis_ids(face_basis_position+1:face_basis_position+size(production_faces(p)%basis_ids_plus))=&
        production_faces(p)%basis_ids_plus;face_basis_position=face_basis_position+size(production_faces(p)%basis_ids_plus)
      checkpoint_payload%face_basis_offsets(face_slot+1)=face_basis_position+1
      checkpoint_payload%face_weights(face_weight_position+1:face_weight_position+size(production_faces(p)%weights))=&
        production_faces(p)%weights;face_weight_position=face_weight_position+size(production_faces(p)%weights)
      checkpoint_payload%face_weight_offsets(face_slot+1)=face_weight_position+1
      call pack_checkpoint_face_values(production_faces(p),checkpoint_payload%face_values,face_value_position)
      checkpoint_payload%face_value_offsets(face_slot+1)=face_value_position+1
      face_observable_count=size(production_faces(p)%weights)**2
      checkpoint_payload%interface_observables(:,face_observable_position+1:&
        face_observable_position+face_observable_count)=transpose(interface_state(&
        interface_cursor+1:interface_cursor+face_observable_count,:))
      interface_cursor=interface_cursor+face_observable_count
      face_observable_position=face_observable_position+face_observable_count
      checkpoint_payload%face_observable_offsets(face_slot+1)=face_observable_position+1
    enddo

    ! Named checkpoint-v3 construction provenance.
    allocate(construction_ids(size(effective_ids)),construction_generations(size(effective_ids)),&
      construction_ordering(size(effective_ids)),construction_ownership(size(effective_ids)))
    construction_ids=int(effective_ids,8);construction_generations=1
    do p=1,size(effective_ids)
      construction_ordering(p)=p
    enddo
    construction_ownership=0
    do p=1,size(row_ids);construction_ownership(int(row_ids(p)))=rank_local+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,construction_ownership,size(construction_ownership),MPI_INTEGER,&
      MPI_SUM,dc%icomm_tot,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.any(construction_ownership<=0))&
      error stop 'DG construction ownership receipt failed'
    catalog_ids_fingerprint=checkpoint_integer64_fingerprint(construction_ids)
    catalog_generation_fingerprint=checkpoint_integer_fingerprint(construction_generations)
    catalog_ordering_fingerprint=checkpoint_integer_fingerprint(construction_ordering)
    catalog_ownership_fingerprint=checkpoint_integer_fingerprint(construction_ownership)
    catalog_provenance_fingerprint=fixed_payload%basis_fingerprint
    catalog_fingerprint=int(z'6A09E667F3BCC909',8)
    catalog_fingerprint=ieor(ishftc(catalog_fingerprint,7),catalog_ids_fingerprint)
    catalog_fingerprint=ieor(ishftc(catalog_fingerprint,7),catalog_generation_fingerprint)
    catalog_fingerprint=ieor(ishftc(catalog_fingerprint,7),catalog_ordering_fingerprint)
    catalog_fingerprint=ieor(ishftc(catalog_fingerprint,7),catalog_ownership_fingerprint)
    catalog_fingerprint=ieor(ishftc(catalog_fingerprint,7),catalog_provenance_fingerprint)
    if(catalog_fingerprint==0_8)catalog_fingerprint=1_8
    checkpoint_payload%construction_catalog%valid=.true.
    checkpoint_payload%construction_catalog%global_count=size(effective_ids)
    allocate(checkpoint_payload%construction_catalog%ids,source=construction_ids)
    allocate(checkpoint_payload%construction_catalog%generations,source=construction_generations)
    allocate(checkpoint_payload%construction_catalog%ordering,source=construction_ordering)
    allocate(checkpoint_payload%construction_catalog%ownership,source=construction_ownership)
    checkpoint_payload%construction_catalog%ids_fingerprint=catalog_ids_fingerprint
    checkpoint_payload%construction_catalog%generation_fingerprint=catalog_generation_fingerprint
    checkpoint_payload%construction_catalog%ordering_fingerprint=catalog_ordering_fingerprint
    checkpoint_payload%construction_catalog%ownership_fingerprint=catalog_ownership_fingerprint
    checkpoint_payload%construction_catalog%provenance_fingerprint=catalog_provenance_fingerprint
    checkpoint_payload%construction_catalog%catalog_fingerprint=catalog_fingerprint
    checkpoint_payload%catalog_fingerprint=catalog_fingerprint

    ! The fixed-rank certified embedding remains row-owned in construction
    ! space.  U and every RT matrix are independently cyclic-sharded.
    call collect_dg_hybrid_full_rows(dc%icomm_tot,size(effective_ids),row_ids,&
      certified_rt_basis%b_rt,full_b_rt,local_ok)
    if(.not.local_ok)error stop 'DG certified RT embedding collection failed'
    rt_metric=certified_rt_basis%metric_rt
    rt_kinetic=certified_rt_basis%scalar_operators_rt(:,:,1)
    rt_nonlocal=certified_rt_basis%scalar_operators_rt(:,:,2)
    rt_local=certified_rt_basis%scalar_operators_rt(:,:,3)
    rt_sipg=certified_rt_basis%scalar_operators_rt(:,:,4)
    rt_hamiltonian=rt_kinetic+rt_nonlocal+rt_local+rt_sipg
    allocate(rt_basis_values(certified_rank,size(interior_values,2)))
    rt_basis_values=matmul(transpose(full_b_rt),interior_values)
    rt_row_count=0
    do p=1,certified_rank
      if(mod(p-1,nproc_local)==rank_local)rt_row_count=rt_row_count+1
    enddo
    allocate(rt_row_ids(rt_row_count),rt_row_owner_keys(certified_rank))
    rt_row_position=0
    do p=1,certified_rank
      rt_row_owner_keys(p)=mod(p-1,nproc_local)+1
      if(mod(p-1,nproc_local)/=rank_local)cycle
      rt_row_position=rt_row_position+1;rt_row_ids(rt_row_position)=int(p,8)
    enddo

    checkpoint_payload%certified_basis%valid=.true.
    checkpoint_payload%certified_basis%localization_converged=certified_rt_basis%localization_converged
    checkpoint_payload%certified_basis%localization_symmetry_constrained=&
      certified_rt_basis%localization_symmetry_constrained
    checkpoint_payload%certified_basis%construction_count=size(effective_ids)
    checkpoint_payload%certified_basis%certified_count=spectral_certification%certified_rank
    checkpoint_payload%certified_basis%occupied_count=occupation_result%noccupied
    checkpoint_payload%certified_basis%localization_iterations=certified_rt_basis%localization_iterations
    allocate(checkpoint_payload%certified_basis%construction_row_ids,source=row_ids)
    allocate(checkpoint_payload%certified_basis%transformation_row_ids,source=rt_row_ids)
    allocate(checkpoint_payload%certified_basis%c_cert,source=certified_rt_basis%c_cert)
    allocate(checkpoint_payload%certified_basis%b_rt,source=certified_rt_basis%b_rt)
    allocate(checkpoint_payload%certified_basis%u_rt(rt_row_count,certified_rank))
    do p=1,rt_row_count
      checkpoint_payload%certified_basis%u_rt(p,:)=certified_rt_basis%u_rt(int(rt_row_ids(p)),:)
    enddo
    allocate(checkpoint_payload%certified_basis%initial_occupied_amplitudes,&
      source=certified_rt_basis%initial_occupied_amplitudes)
    allocate(checkpoint_payload%certified_basis%certified_eigenvalues,&
      source=certified_rt_basis%certified_eigenvalues)
    allocate(checkpoint_payload%certified_basis%occupations,&
      source=occupation_result%occupations(1:occupation_result%noccupied))
    allocate(checkpoint_payload%certified_basis%centers,source=certified_rt_basis%centers)
    allocate(checkpoint_payload%certified_basis%spreads_before,source=certified_rt_basis%spreads_before)
    allocate(checkpoint_payload%certified_basis%spreads_after,source=certified_rt_basis%spreads_after)
    checkpoint_payload%certified_basis%spread_before_total=certified_rt_basis%spread_before_total
    checkpoint_payload%certified_basis%spread_after_total=certified_rt_basis%spread_after_total
    checkpoint_payload%certified_basis%spread_improvement=certified_rt_basis%spread_improvement
    checkpoint_payload%certified_basis%transform_unitarity_defect=certified_rt_basis%transform_unitarity_defect
    checkpoint_payload%certified_basis%certified_metric_defect=certified_rt_basis%certified_metric_defect
    checkpoint_payload%certified_basis%rt_metric_defect=certified_rt_basis%rt_metric_defect
    checkpoint_payload%certified_basis%embedding_defect=certified_rt_basis%embedding_defect
    checkpoint_payload%certified_basis%projector_invariance_defect=certified_rt_basis%projector_invariance_defect
    checkpoint_payload%certified_basis%target_symmetry_defect_before=&
      certified_rt_basis%target_symmetry_defect_before
    checkpoint_payload%certified_basis%target_symmetry_defect_after=&
      certified_rt_basis%target_symmetry_defect_after
    checkpoint_payload%certified_basis%energy_symmetry_defect_before=&
      certified_rt_basis%energy_symmetry_defect_before
    checkpoint_payload%certified_basis%energy_symmetry_defect_after=&
      certified_rt_basis%energy_symmetry_defect_after
    checkpoint_payload%certified_basis%symmetry_defect_invariance=certified_rt_basis%symmetry_defect_invariance
    checkpoint_payload%certified_basis%scalar_covariance_defect=certified_rt_basis%scalar_covariance_defect
    checkpoint_payload%certified_basis%vector_covariance_defect=certified_rt_basis%vector_covariance_defect
    checkpoint_payload%certified_basis%tensor_covariance_defect=certified_rt_basis%tensor_covariance_defect
    initial_state_fingerprint=checkpoint_complex_matrix_fingerprint(&
      certified_rt_basis%initial_occupied_amplitudes)
    checkpoint_payload%certified_basis%c_cert_fingerprint=certified_rt_basis%c_cert_fingerprint
    checkpoint_payload%certified_basis%u_rt_fingerprint=certified_rt_basis%localization_fingerprint
    checkpoint_payload%certified_basis%b_rt_fingerprint=certified_rt_basis%b_rt_fingerprint
    checkpoint_payload%certified_basis%initial_state_fingerprint=initial_state_fingerprint
    checkpoint_payload%certified_basis%transformation_fingerprint=certified_rt_basis%localization_fingerprint
    checkpoint_payload%certified_basis%operator_fingerprint=certified_rt_basis%operator_fingerprint
    checkpoint_payload%certified_basis%fingerprint=certified_rt_basis%fingerprint

    checkpoint_payload%electron_count%valid=.true.
    checkpoint_payload%electron_count%expected_count=dc%elec_num_tot
    checkpoint_payload%electron_count%actual_count=occupation_result%electron_count
    checkpoint_payload%electron_count%tolerance=dg_dc_gs_electron_count_tolerance
    checkpoint_payload%electron_count%defect=abs(dc%elec_num_tot-occupation_result%electron_count)
    checkpoint_payload%electron_count%omitted_tail=occupation_result%omitted_occupation_tail
    checkpoint_payload%electron_count%chemical_potential=occupation_result%chemical_potential
    checkpoint_payload%electron_count%fingerprint=occupation_result%fingerprint

    checkpoint_payload%rt_space%valid=.true.
    checkpoint_payload%rt_space%rank=spectral_certification%certified_rank
    checkpoint_payload%rt_space%operation_count=size(certified_rt_basis%representation_rt,3)
    checkpoint_payload%rt_space%scalar_count=size(certified_rt_basis%scalar_operators_rt,3)
    checkpoint_payload%rt_space%vector_count=size(certified_rt_basis%vector_operators_rt,4)
    checkpoint_payload%rt_space%tensor_count=size(certified_rt_basis%tensor_operators_rt,5)
    allocate(checkpoint_payload%rt_space%row_ids,source=rt_row_ids)
    allocate(checkpoint_payload%rt_space%row_owner_keys,source=rt_row_owner_keys)
    allocate(checkpoint_payload%rt_space%grid_owner_keys(size(ow_core_ids)),source=rank_local+1)
    allocate(checkpoint_payload%rt_space%metric_rows(rt_row_count,certified_rank),&
      checkpoint_payload%rt_space%kinetic_rows(rt_row_count,certified_rank),&
      checkpoint_payload%rt_space%nonlocal_rows(rt_row_count,certified_rank),&
      checkpoint_payload%rt_space%local_rows(rt_row_count,certified_rank),&
      checkpoint_payload%rt_space%sipg_rows(rt_row_count,certified_rank),&
      checkpoint_payload%rt_space%hamiltonian_rows(rt_row_count,certified_rank),&
      checkpoint_payload%rt_space%scalar_operator_rows(rt_row_count,certified_rank,&
        size(certified_rt_basis%scalar_operators_rt,3)),&
      checkpoint_payload%rt_space%vector_operator_rows(rt_row_count,certified_rank,3,&
        size(certified_rt_basis%vector_operators_rt,4)),&
      checkpoint_payload%rt_space%tensor_operator_rows(rt_row_count,certified_rank,3,3,&
        size(certified_rt_basis%tensor_operators_rt,5)))
    do p=1,rt_row_count
      construction_index=int(rt_row_ids(p))
      checkpoint_payload%rt_space%metric_rows(p,:)=rt_metric(construction_index,:)
      checkpoint_payload%rt_space%kinetic_rows(p,:)=rt_kinetic(construction_index,:)
      checkpoint_payload%rt_space%nonlocal_rows(p,:)=rt_nonlocal(construction_index,:)
      checkpoint_payload%rt_space%local_rows(p,:)=rt_local(construction_index,:)
      checkpoint_payload%rt_space%sipg_rows(p,:)=rt_sipg(construction_index,:)
      checkpoint_payload%rt_space%hamiltonian_rows(p,:)=rt_hamiltonian(construction_index,:)
      checkpoint_payload%rt_space%scalar_operator_rows(p,:,:)=&
        certified_rt_basis%scalar_operators_rt(construction_index,:,:)
      checkpoint_payload%rt_space%vector_operator_rows(p,:,:,:)=&
        certified_rt_basis%vector_operators_rt(construction_index,:,:,:)
      checkpoint_payload%rt_space%tensor_operator_rows(p,:,:,:,:)=&
        certified_rt_basis%tensor_operators_rt(construction_index,:,:,:,:)
    enddo
    allocate(checkpoint_payload%rt_space%representation,source=certified_rt_basis%representation_rt)
    allocate(checkpoint_payload%rt_space%cartesian_rotations,source=certified_rt_basis%cartesian_rotations)
    allocate(checkpoint_payload%rt_space%basis_values,source=rt_basis_values)
    allocate(checkpoint_payload%rt_space%density,source=rho_in)
    call fingerprint_rt_dg_hybrid_component(dc%icomm_tot,rt_row_ids,&
      checkpoint_payload%rt_space%metric_rows,rt_metric_fingerprint,local_ok)
    if(local_ok)call fingerprint_rt_dg_hybrid_component(dc%icomm_tot,rt_row_ids,&
      checkpoint_payload%rt_space%kinetic_rows,rt_kinetic_fingerprint,local_ok)
    if(local_ok)call fingerprint_rt_dg_hybrid_component(dc%icomm_tot,rt_row_ids,&
      checkpoint_payload%rt_space%nonlocal_rows,rt_nonlocal_fingerprint,local_ok)
    if(local_ok)call fingerprint_rt_dg_hybrid_component(dc%icomm_tot,rt_row_ids,&
      checkpoint_payload%rt_space%local_rows,rt_local_fingerprint,local_ok)
    if(local_ok)call fingerprint_rt_dg_hybrid_component(dc%icomm_tot,rt_row_ids,&
      checkpoint_payload%rt_space%sipg_rows,rt_sipg_fingerprint,local_ok)
    if(local_ok)call fingerprint_rt_dg_hybrid_component(dc%icomm_tot,rt_row_ids,&
      checkpoint_payload%rt_space%hamiltonian_rows,rt_hamiltonian_fingerprint,local_ok)
    if(.not.local_ok)error stop 'DG RT component fingerprint failed'
    call checkpoint_grid_complex_fingerprint(dc%icomm_tot,ow_core_ids,rt_basis_values,&
      rt_basis_fingerprint,local_ok)
    if(local_ok)call checkpoint_grid_real_fingerprint(dc%icomm_tot,ow_core_ids,rho_in,&
      rt_density_fingerprint,local_ok)
    if(.not.local_ok)error stop 'DG RT grid fingerprint failed'
    call checkpoint_grid_integer_fingerprint(dc%icomm_tot,ow_core_ids,&
      checkpoint_payload%rt_space%grid_owner_keys,rt_grid_ownership_fingerprint,local_ok)
    if(.not.local_ok)error stop 'DG RT grid-ownership fingerprint failed'
    rt_ownership_fingerprint=checkpoint_integer_fingerprint(rt_row_owner_keys)
    rt_ownership_fingerprint=ieor(ishftc(rt_ownership_fingerprint,7),rt_grid_ownership_fingerprint)
    if(rt_ownership_fingerprint==0_8)rt_ownership_fingerprint=1_8
    rt_scalar_fingerprint=checkpoint_analysis_fingerprint(certified_rt_basis%scalar_operators_rt)
    rt_vector_fingerprint=checkpoint_complex_rank4_fingerprint(certified_rt_basis%vector_operators_rt)
    rt_tensor_fingerprint=checkpoint_complex_rank5_fingerprint(certified_rt_basis%tensor_operators_rt)
    rt_representation_fingerprint=checkpoint_analysis_fingerprint(certified_rt_basis%representation_rt)
    rt_rotation_fingerprint=checkpoint_real_rank3_fingerprint(certified_rt_basis%cartesian_rotations)
    rt_payload_fingerprint=int(z'5BE0CD19137E2179',8)
    rt_payload_fingerprint=ieor(ishftc(rt_payload_fingerprint,7),rt_metric_fingerprint)
    rt_payload_fingerprint=ieor(ishftc(rt_payload_fingerprint,7),rt_kinetic_fingerprint)
    rt_payload_fingerprint=ieor(ishftc(rt_payload_fingerprint,7),rt_nonlocal_fingerprint)
    rt_payload_fingerprint=ieor(ishftc(rt_payload_fingerprint,7),rt_local_fingerprint)
    rt_payload_fingerprint=ieor(ishftc(rt_payload_fingerprint,7),rt_sipg_fingerprint)
    rt_payload_fingerprint=ieor(ishftc(rt_payload_fingerprint,7),rt_hamiltonian_fingerprint)
    rt_payload_fingerprint=ieor(ishftc(rt_payload_fingerprint,7),rt_basis_fingerprint)
    rt_payload_fingerprint=ieor(ishftc(rt_payload_fingerprint,7),rt_density_fingerprint)
    rt_payload_fingerprint=ieor(ishftc(rt_payload_fingerprint,7),rt_ownership_fingerprint)
    rt_payload_fingerprint=ieor(ishftc(rt_payload_fingerprint,7),rt_scalar_fingerprint)
    rt_payload_fingerprint=ieor(ishftc(rt_payload_fingerprint,7),rt_vector_fingerprint)
    rt_payload_fingerprint=ieor(ishftc(rt_payload_fingerprint,7),rt_tensor_fingerprint)
    rt_payload_fingerprint=ieor(ishftc(rt_payload_fingerprint,7),rt_representation_fingerprint)
    rt_payload_fingerprint=ieor(ishftc(rt_payload_fingerprint,7),rt_rotation_fingerprint)
    if(rt_payload_fingerprint==0_8)rt_payload_fingerprint=1_8
    checkpoint_payload%rt_space%metric_fingerprint=rt_metric_fingerprint
    checkpoint_payload%rt_space%kinetic_fingerprint=rt_kinetic_fingerprint
    checkpoint_payload%rt_space%nonlocal_fingerprint=rt_nonlocal_fingerprint
    checkpoint_payload%rt_space%local_fingerprint=rt_local_fingerprint
    checkpoint_payload%rt_space%sipg_fingerprint=rt_sipg_fingerprint
    checkpoint_payload%rt_space%hamiltonian_fingerprint=rt_hamiltonian_fingerprint
    checkpoint_payload%rt_space%basis_fingerprint=rt_basis_fingerprint
    checkpoint_payload%rt_space%density_fingerprint=rt_density_fingerprint
    checkpoint_payload%rt_space%ownership_fingerprint=rt_ownership_fingerprint
    checkpoint_payload%rt_space%scalar_fingerprint=rt_scalar_fingerprint
    checkpoint_payload%rt_space%vector_fingerprint=rt_vector_fingerprint
    checkpoint_payload%rt_space%tensor_fingerprint=rt_tensor_fingerprint
    checkpoint_payload%rt_space%representation_fingerprint=rt_representation_fingerprint
    checkpoint_payload%rt_space%fingerprint=rt_payload_fingerprint

    checkpoint_payload%energy_window%valid=.true.
    checkpoint_payload%energy_window%compatibility_dynamic_rank=&
      spectral_certification%compatibility_dynamic_rank
    checkpoint_payload%energy_window%proof_state_present=spectral_certification%proof_state_present
    checkpoint_payload%energy_window%mode=merge(rt_dg_hybrid_energy_window_legacy_dynamic,&
      rt_dg_hybrid_energy_window_explicit,spectral_certification%compatibility_dynamic_rank)
    checkpoint_payload%energy_window%construction_rank=size(effective_ids)
    checkpoint_payload%energy_window%solved_rank=size(effective_ids)
    checkpoint_payload%energy_window%occupied_rank=occupation_result%noccupied
    checkpoint_payload%energy_window%requested_rank=spectral_certification%requested_rank
    checkpoint_payload%energy_window%certified_rank=spectral_certification%certified_rank
    checkpoint_payload%energy_window%extension_states=spectral_certification%extension_states
    checkpoint_payload%energy_window%boundary_cluster_rank=spectral_certification%boundary_cluster_rank
    checkpoint_payload%energy_window%proof_status=merge(1,0,spectral_certification%proof_state_present)
    checkpoint_payload%energy_window%window_size=spectral_certification%energy_window
    checkpoint_payload%energy_window%e_homo=spectral_certification%e_homo
    checkpoint_payload%energy_window%requested_cutoff=spectral_certification%requested_cutoff
    checkpoint_payload%energy_window%certified_cutoff=spectral_certification%certified_cutoff
    checkpoint_payload%energy_window%extension_energy=spectral_certification%extension_energy
    checkpoint_payload%energy_window%proof_energy=spectral_certification%proof_energy
    checkpoint_payload%energy_window%fingerprint=spectral_certification%fingerprint

    checkpoint_payload%symmetry_receipt%valid=.true.
    checkpoint_payload%symmetry_receipt%worst_operation=max(1,spectral_certification%worst_operation)
    checkpoint_payload%symmetry_receipt%occupied_subspace_defect=spectral_certification%occupied_subspace_defect
    checkpoint_payload%symmetry_receipt%occupied_projector_defect=projector_symmetry_residual
    checkpoint_payload%symmetry_receipt%target_subspace_defect=spectral_certification%target_subspace_defect
    checkpoint_payload%symmetry_receipt%target_energy_defect=spectral_certification%target_energy_defect
    checkpoint_payload%symmetry_receipt%density_defect=density_symmetry_defect
    checkpoint_payload%symmetry_receipt%scalar_covariance_defect=certified_rt_basis%scalar_covariance_defect
    checkpoint_payload%symmetry_receipt%vector_covariance_defect=certified_rt_basis%vector_covariance_defect
    checkpoint_payload%symmetry_receipt%tensor_covariance_defect=certified_rt_basis%tensor_covariance_defect
    checkpoint_payload%symmetry_receipt%final_basis_defect=max(certified_rt_basis%rt_metric_defect,&
      certified_rt_basis%embedding_defect,certified_rt_basis%projector_invariance_defect)
    checkpoint_payload%symmetry_receipt%worst_operation_defect=spectral_certification%worst_operation_defect
    checkpoint_payload%symmetry_receipt%maximum_physical_defect=max(&
      spectral_certification%maximum_physical_defect,certified_rt_basis%scalar_covariance_defect,&
      certified_rt_basis%vector_covariance_defect,certified_rt_basis%tensor_covariance_defect,&
      checkpoint_payload%symmetry_receipt%final_basis_defect)
    symmetry_receipt_fingerprint=checkpoint_real_fingerprint([&
      checkpoint_payload%symmetry_receipt%occupied_subspace_defect,&
      checkpoint_payload%symmetry_receipt%occupied_projector_defect,&
      checkpoint_payload%symmetry_receipt%target_subspace_defect,&
      checkpoint_payload%symmetry_receipt%target_energy_defect,&
      checkpoint_payload%symmetry_receipt%density_defect,&
      checkpoint_payload%symmetry_receipt%scalar_covariance_defect,&
      checkpoint_payload%symmetry_receipt%vector_covariance_defect,&
      checkpoint_payload%symmetry_receipt%tensor_covariance_defect,&
      checkpoint_payload%symmetry_receipt%final_basis_defect,&
      checkpoint_payload%symmetry_receipt%worst_operation_defect,&
      checkpoint_payload%symmetry_receipt%maximum_physical_defect])
    symmetry_receipt_fingerprint=ieor(ishftc(symmetry_receipt_fingerprint,7),&
      spectral_certification%fingerprint)
    symmetry_receipt_fingerprint=ieor(ishftc(symmetry_receipt_fingerprint,7),&
      certified_rt_basis%fingerprint)
    symmetry_receipt_fingerprint=ieor(ishftc(symmetry_receipt_fingerprint,7),&
      int(checkpoint_payload%symmetry_receipt%worst_operation,8))
    if(symmetry_receipt_fingerprint==0_8)symmetry_receipt_fingerprint=1_8
    checkpoint_payload%symmetry_receipt%fingerprint=symmetry_receipt_fingerprint

    checkpoint_payload%handoff_receipts%valid=.true.
    checkpoint_payload%handoff_receipts%position_fingerprint=&
      checkpoint_payload%position_convention_fingerprint
    checkpoint_payload%handoff_receipts%nonlocal_fingerprint=checkpoint_payload%nonlocal_fingerprint
    checkpoint_payload%handoff_receipts%face_fingerprint=checkpoint_payload%face_fingerprint
    checkpoint_payload%handoff_receipts%pseudopotential_fingerprint=&
      checkpoint_payload%pseudopotential_fingerprint
    checkpoint_payload%handoff_receipts%transformation_fingerprint=&
      checkpoint_payload%certified_basis%transformation_fingerprint
    handoff_fingerprint=ieor(checkpoint_payload%position_convention_fingerprint,&
      ishftc(checkpoint_payload%nonlocal_fingerprint,7))
    handoff_fingerprint=ieor(handoff_fingerprint,ishftc(checkpoint_payload%face_fingerprint,13))
    handoff_fingerprint=ieor(handoff_fingerprint,&
      ishftc(checkpoint_payload%pseudopotential_fingerprint,17))
    handoff_fingerprint=ieor(handoff_fingerprint,&
      ishftc(checkpoint_payload%certified_basis%transformation_fingerprint,19))
    if(handoff_fingerprint==0_8)handoff_fingerprint=1_8
    checkpoint_payload%handoff_receipts%fingerprint=handoff_fingerprint

    call fingerprint_rt_dg_hybrid_ground_state_payload(dc%icomm_tot,checkpoint_payload,&
      checkpoint_fingerprint,payload_ready,continuation_message)
    call authorize_dg_hybrid_v3_publication(dc%icomm_tot,candidate_acceptance,&
      rt_dg_hybrid_ground_state_checkpoint_version,checkpoint_payload%rt_space%rank,&
      payload_ready,local_ok,continuation_message)
    if(.not.local_ok.or..not.payload_ready)then
      write(0,'(a)')trim(continuation_message)
      error stop 'DG checkpoint-v3 publication authorization failed'
    endif
    call write_rt_dg_hybrid_ground_state_checkpoint(dc%icomm_tot,'./hybrid_dg_ground_state.chk',checkpoint_payload,&
      checkpoint_fingerprint,local_ok,continuation_message)
    if(.not.local_ok)then;write(0,'(a)')trim(continuation_message);error stop 'DG continuation complete checkpoint failed';endif
    window_mode='explicit'
    if(spectral_certification%compatibility_dynamic_rank)window_mode='legacy_dynamic'
    if(rank_local==0)write(*,'(a,2(a,a),2(a,i0),6(a,es24.16),3(a,i0))')&
      '[HYBRID-LCFO-WINDOW] mode=',trim(window_mode),&
      ' compatibility=',merge('1','0',spectral_certification%compatibility_dynamic_rank),&
      ' requested_rank=',spectral_certification%requested_rank,&
      ' certified_rank=',spectral_certification%certified_rank,&
      ' delta_e=',spectral_certification%energy_window,&
      ' homo=',spectral_certification%e_homo,&
      ' requested_cutoff=',spectral_certification%requested_cutoff,&
      ' certified_cutoff=',spectral_certification%certified_cutoff,&
      ' extension_energy=',spectral_certification%extension_energy,&
      ' proof_energy=',spectral_certification%proof_energy,&
      ' boundary_rank=',spectral_certification%boundary_cluster_rank,&
      ' extension_states=',spectral_certification%extension_states,&
      ' proof_state=',merge(1,0,spectral_certification%proof_state_present)
    if(rank_local==0)write(*,'(a,5(a,es16.8),a,i0)')&
      '[HYBRID-LCFO-SYMMETRY] occupied=',occupied_symmetry_defect,&
      ' target=',target_symmetry_defect,' energy=',target_energy_symmetry_defect,&
      ' density=',density_symmetry_defect,&
      ' maximum=',checkpoint_payload%symmetry_receipt%maximum_physical_defect,&
      ' worst_operation=',checkpoint_payload%symmetry_receipt%worst_operation
    if(rank_local==0)write(*,'(a,3(a,i0),a,es16.8,a,i0,a,es16.8,a,i0)')&
      '[HYBRID-RT-BASIS] construction_rank=',size(effective_ids),&
      ' certified_rank=',spectral_certification%certified_rank,&
      ' rt_rank=',checkpoint_payload%rt_space%rank,&
      ' localization_spread=',certified_rt_basis%spread_after_total,&
      ' embedding_fingerprint=',certified_rt_basis%b_rt_fingerprint,&
      ' operator_covariance=',max(certified_rt_basis%scalar_covariance_defect,&
        certified_rt_basis%vector_covariance_defect,certified_rt_basis%tensor_covariance_defect),&
      ' checkpoint_version=',rt_dg_hybrid_ground_state_checkpoint_version
    if(rank_local==0)write(*,'(a,4(a,es16.8))')'[OW-GS] fully refreshed DG continuation converged',&
      ' lambda=',accepted_lambda,' h_residual=',residuals%r_h,' density_residual=',residuals%r_rho,&
      ' interface_residual=',residuals%r_t
    if(rank_local==0)write(*,'(a,6(a,i0),12(a,es16.8),a,i0)')&
      '[HYBRID-GS-ACCEPTANCE] seed_identity=',merge(1,0,seed_identity_accepted),&
      ' lambda_zero=',merge(1,0,lambda_zero_accepted),' lambda_one=',merge(1,0,accepted_lambda==1d0),&
      ' final_refresh=',merge(1,0,final_refresh_performed),&
      ' requested_target_rank=',spectral_certification%requested_rank,&
      ' extended_target_rank=',spectral_certification%certified_rank,&
      ' r_h=',residuals%r_h,&
      ' r_rho=',residuals%r_rho,' r_t=',residuals%r_t,' r_s=',residuals%r_s,&
      ' electron=',abs(electron_count-occupation_result%electron_count),' symmetry=',physical_symmetry_defect,&
      ' real_space=',max(real_space_residual,maxval(interface_action_residuals)),&
      ' occupied_defect=',occupied_symmetry_defect,' target_defect=',target_symmetry_defect,&
      ' target_energy_defect=',target_energy_symmetry_defect,' density_defect=',density_symmetry_defect,&
      ' full_operator_defect=',symmetry_residual,' payload_fingerprint=',checkpoint_fingerprint
  end subroutine run_dg_hybrid_concrete_continuation

  subroutine dg_hybrid_certified_localizer_candidate(comm_arg,global_count_arg,certified_rank,&
      row_ids_arg,c_cert_arg,require_unconstrained,transform,centers,spreads_before,spreads_after,&
      iterations,symmetry_constrained,converged,callback_ok,callback_message)
    integer,intent(in)::comm_arg,global_count_arg,certified_rank
    integer(8),intent(in)::row_ids_arg(:)
    complex(8),intent(in)::c_cert_arg(:,:)
    logical,intent(in)::require_unconstrained
    complex(8),allocatable,intent(out)::transform(:,:)
    real(8),allocatable,intent(out)::centers(:,:),spreads_before(:),spreads_after(:)
    integer,intent(out)::iterations
    logical,intent(out)::symmetry_constrained,converged,callback_ok
    character(*),intent(out)::callback_message
    complex(8),allocatable::full_c(:,:),certified_values(:,:),anchors(:,:),m_matrix(:,:,:),a_matrix(:,:)
    real(8),allocatable::fractional(:,:),eigenvalues_arg(:),atom_positions(:,:)
    real(8)::lattice_inverse(3,3),reciprocal_lattice(3,3),determinant,spread_receipt(3)
    integer,allocatable::nncell(:,:)
    character(8),allocatable::atom_symbols(:)
    integer::nntot,point,atom,ierr_local,ierr_status,local_bad,global_bad
    integer(8)::coordinator_bytes,workspace_peak,byte_limit,nxy_local
    logical::local_ok

    callback_ok=.false.;callback_message='';iterations=0
    symmetry_constrained=.false.;converged=.false.
    local_bad=merge(0,1,require_unconstrained)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr_status)
    if(ierr_status/=MPI_SUCCESS.or.global_bad/=0)then
      callback_message='certified RT localizer requires the unconstrained mode';return
    endif
    local_bad=merge(0,1,allocated(dg_certified_localizer_values).and.&
      allocated(dg_certified_localizer_weights).and.allocated(dg_certified_localizer_grid_ids))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr_status)
    if(ierr_status/=MPI_SUCCESS.or.global_bad/=0)then
      callback_message='certified RT localizer context is unavailable';return
    endif
    local_bad=merge(0,1,size(dg_certified_localizer_values,1)==global_count_arg.and.&
      size(dg_certified_localizer_values,2)==size(dg_certified_localizer_grid_ids).and.&
      size(dg_certified_localizer_weights)==size(dg_certified_localizer_grid_ids).and.&
      size(c_cert_arg,2)==certified_rank)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr_status)
    if(ierr_status/=MPI_SUCCESS.or.global_bad/=0)then
      callback_message='certified RT localizer context shape mismatch';return
    endif
    call collect_dg_hybrid_full_rows(comm_arg,global_count_arg,row_ids_arg,c_cert_arg,full_c,local_ok)
    if(.not.local_ok)then;callback_message='certified coefficient collection failed';return;endif
    allocate(certified_values(certified_rank,size(dg_certified_localizer_grid_ids)))
    certified_values=matmul(transpose(full_c),dg_certified_localizer_values)
    allocate(anchors,source=certified_values)
    allocate(fractional(3,size(dg_certified_localizer_grid_ids)))
    nxy_local=int(dc%lg_tot%num(1),8)*int(dc%lg_tot%num(2),8)
    do point=1,size(dg_certified_localizer_grid_ids)
      fractional(1,point)=real(modulo(dg_certified_localizer_grid_ids(point)-1_8,&
        int(dc%lg_tot%num(1),8)),8)/real(dc%lg_tot%num(1),8)
      fractional(2,point)=real(modulo((dg_certified_localizer_grid_ids(point)-1_8)/&
        int(dc%lg_tot%num(1),8),int(dc%lg_tot%num(2),8)),8)/real(dc%lg_tot%num(2),8)
      fractional(3,point)=real((dg_certified_localizer_grid_ids(point)-1_8)/nxy_local,8)/&
        real(dc%lg_tot%num(3),8)
    enddo
    call measure_dg_hybrid_periodic_spreads(comm_arg,certified_values,&
      dg_certified_localizer_weights,fractional,dc%system_tot%primitive_a,spreads_before,local_ok)
    if(.not.local_ok)then;callback_message='certified initial spread measurement failed';return;endif
    call invert_ow_lattice(dc%system_tot%primitive_a,lattice_inverse,determinant,local_ok)
    if(.not.local_ok)then;callback_message='certified RT Wannier lattice is singular';return;endif
    reciprocal_lattice=2d0*pi*transpose(lattice_inverse)
    allocate(atom_symbols(dc%system_tot%nion),atom_positions(3,dc%system_tot%nion))
    atom_positions=dc%system_tot%Rion
    do atom=1,dc%system_tot%nion
      if(dc%system_tot%kion(atom)<1.or.dc%system_tot%kion(atom)>size(pp%atom_symbol))then
        callback_message='certified RT atom species is outside the pseudopotential table';return
      endif
      atom_symbols(atom)=pp%atom_symbol(dc%system_tot%kion(atom))
    enddo
    call setup_dg_w90_gamma_library(comm_arg,'hybrid_certified_rt_mlwf',&
      dc%system_tot%primitive_a,reciprocal_lattice,atom_symbols,atom_positions,&
      certified_rank,certified_rank,wannier_num_iter,dg_ow_w90_initial_projection,&
      DG_W90_UNCONSTRAINED,nntot,nncell,local_ok,callback_message)
    if(.not.local_ok)return
    byte_limit=8_8*1024_8*1024_8*1024_8
    call assemble_dg_w90_gamma_matrices(comm_arg,certified_values,anchors,&
      dg_certified_localizer_weights,fractional,nncell,byte_limit,m_matrix,a_matrix,&
      coordinator_bytes,workspace_peak,local_ok,callback_message)
    if(.not.local_ok)return
    allocate(eigenvalues_arg(certified_rank));eigenvalues_arg=0d0
    call run_dg_w90_gamma_library(comm_arg,'hybrid_certified_rt_mlwf',&
      dc%system_tot%primitive_a,reciprocal_lattice,atom_symbols,atom_positions,&
      m_matrix,a_matrix,eigenvalues_arg,sum(spreads_before),dg_ow_symmetry_tolerance,&
      wannier_num_iter,transform,centers,spreads_after,spread_receipt,local_ok,&
      callback_message,convergence_iterations_out=iterations,require_nonincreasing_spread=.false.)
    if(.not.local_ok)return
    centers=matmul(lattice_inverse,centers)
    call apply_dg_w90_gamma_transform(comm_arg,dg_certified_localizer_grid_ids,&
      certified_values,transform=transform,centers=centers,tolerance=dg_ow_symmetry_tolerance,&
      ok=local_ok,message=callback_message,spreads=spreads_after)
    if(.not.local_ok)return
    call MPI_Allreduce(merge(0,1,all(ieee_is_finite(spreads_before)).and.&
      all(ieee_is_finite(spreads_after))),ierr_local,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr_status)
    if(ierr_status/=MPI_SUCCESS.or.ierr_local/=0)then
      callback_message='certified RT localization returned a nonfinite spread receipt';return
    endif
    symmetry_constrained=.false.;converged=.true.;callback_ok=.true.;callback_message=''
  end subroutine dg_hybrid_certified_localizer_candidate

  subroutine measure_dg_hybrid_periodic_spreads(comm_arg,values,weights,fractional,lattice,spreads,callback_ok)
    integer,intent(in)::comm_arg
    complex(8),intent(in)::values(:,:)
    real(8),intent(in)::weights(:),fractional(:,:),lattice(3,3)
    real(8),allocatable,intent(out)::spreads(:)
    logical,intent(out)::callback_ok
    complex(8),allocatable::local_phase(:,:),global_phase(:,:)
    real(8),allocatable::local_norm(:),global_norm(:),local_spread(:),global_spread(:)
    real(8)::angle,delta_fractional(3),delta_cartesian(3),probability
    integer::state,point,axis,ierr_local,local_bad,global_bad,local_states,minimum_states,maximum_states
    callback_ok=.false.
    local_bad=merge(0,1,size(values,2)==size(weights).and.&
      all(shape(fractional)==[3,size(weights)]).and.all(weights>=0d0).and.&
      all(ieee_is_finite(weights)).and.all(ieee_is_finite(fractional)).and.&
      all(ieee_is_finite(lattice)).and.all(ieee_is_finite(real(values))).and.&
      all(ieee_is_finite(aimag(values))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.global_bad/=0)return
    local_states=size(values,1)
    call MPI_Allreduce(local_states,minimum_states,1,MPI_INTEGER,MPI_MIN,comm_arg,ierr_local)
    if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(local_states,maximum_states,1,MPI_INTEGER,MPI_MAX,&
      comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.minimum_states/=maximum_states.or.minimum_states<1)return
    allocate(local_phase(3,size(values,1)),global_phase(3,size(values,1)),&
      local_norm(size(values,1)),global_norm(size(values,1)),local_spread(size(values,1)),&
      global_spread(size(values,1)),spreads(size(values,1)))
    local_phase=(0d0,0d0);local_norm=0d0
    do point=1,size(weights);do state=1,size(values,1)
      probability=weights(point)*abs(values(state,point))**2
      local_norm(state)=local_norm(state)+probability
      do axis=1,3
        angle=2d0*pi*fractional(axis,point)
        local_phase(axis,state)=local_phase(axis,state)+probability*cmplx(cos(angle),sin(angle),8)
      enddo
    enddo;enddo
    call MPI_Allreduce(local_phase,global_phase,size(local_phase),MPI_DOUBLE_COMPLEX,MPI_SUM,comm_arg,ierr_local)
    if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(local_norm,global_norm,size(local_norm),&
      MPI_DOUBLE_PRECISION,MPI_SUM,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.any(global_norm<=tiny(1d0)))then;callback_ok=.false.;return;endif
    local_spread=0d0
    do point=1,size(weights);do state=1,size(values,1)
      do axis=1,3
        angle=atan2(aimag(global_phase(axis,state)),real(global_phase(axis,state)))/(2d0*pi)
        delta_fractional(axis)=modulo(fractional(axis,point)-angle+0.5d0,1d0)-0.5d0
      enddo
      delta_cartesian=matmul(lattice,delta_fractional)
      probability=weights(point)*abs(values(state,point))**2
      local_spread(state)=local_spread(state)+probability*sum(delta_cartesian**2)
    enddo;enddo
    call MPI_Allreduce(local_spread,global_spread,size(local_spread),MPI_DOUBLE_PRECISION,MPI_SUM,&
      comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS)then;callback_ok=.false.;return;endif
    spreads=max(0d0,global_spread/global_norm)
    callback_ok=all(ieee_is_finite(spreads))
  end subroutine measure_dg_hybrid_periodic_spreads

  subroutine collect_dg_hybrid_full_rows(comm_arg,global_count_arg,row_ids_arg,row_values,full_values,callback_ok)
    integer,intent(in)::comm_arg,global_count_arg
    integer(8),intent(in)::row_ids_arg(:)
    complex(8),intent(in)::row_values(:,:)
    complex(8),allocatable,intent(out)::full_values(:,:)
    logical,intent(out)::callback_ok
    complex(8),allocatable::local_values(:,:)
    integer::row,ierr_local,local_bad,global_bad,local_shape(2),minimum_shape(2),maximum_shape(2)
    callback_ok=.false.
    local_bad=merge(0,1,global_count_arg>0.and.size(row_values,1)==size(row_ids_arg).and.&
      all(row_ids_arg>=1_8).and.all(row_ids_arg<=int(global_count_arg,8)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.global_bad/=0)return
    local_shape=[global_count_arg,size(row_values,2)]
    call MPI_Allreduce(local_shape,minimum_shape,2,MPI_INTEGER,MPI_MIN,comm_arg,ierr_local)
    if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(local_shape,maximum_shape,2,MPI_INTEGER,MPI_MAX,&
      comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.any(minimum_shape/=maximum_shape))return
    allocate(local_values(global_count_arg,size(row_values,2)),&
      full_values(global_count_arg,size(row_values,2)))
    local_values=(0d0,0d0)
    do row=1,size(row_ids_arg);local_values(int(row_ids_arg(row)),:)=row_values(row,:);enddo
    call MPI_Allreduce(local_values,full_values,size(full_values),MPI_DOUBLE_COMPLEX,MPI_SUM,&
      comm_arg,ierr_local)
    callback_ok=ierr_local==MPI_SUCCESS.and.all(ieee_is_finite(real(full_values))).and.&
      all(ieee_is_finite(aimag(full_values)))
  end subroutine collect_dg_hybrid_full_rows

  subroutine build_dg_hybrid_certified_representation(metric,representation,c_cert,certified_representation)
    complex(8),intent(in)::metric(:,:),representation(:,:,:),c_cert(:,:)
    complex(8),allocatable,intent(out)::certified_representation(:,:,:)
    integer::operation
    allocate(certified_representation(size(c_cert,2),size(c_cert,2),size(representation,3)))
    do operation=1,size(representation,3)
      certified_representation(:,:,operation)=matmul(conjg(transpose(c_cert)),&
        matmul(metric,matmul(representation(:,:,operation),c_cert)))
    enddo
  end subroutine build_dg_hybrid_certified_representation

  subroutine assemble_dg_hybrid_canonical_momentum(comm_arg,weights,values,gradients,momentum,&
      callback_ok,callback_message)
    integer,intent(in)::comm_arg
    real(8),intent(in)::weights(:)
    complex(8),intent(in)::values(:,:),gradients(:,:,:)
    complex(8),allocatable,intent(out)::momentum(:,:,:)
    logical,intent(out)::callback_ok
    character(*),intent(out)::callback_message
    complex(8),allocatable::local_derivative(:,:,:),derivative(:,:,:)
    integer::n,point,row,column,axis,ierr_local,local_bad,global_bad,minimum_n,maximum_n
    callback_ok=.false.;callback_message='';n=size(values,1);local_bad=0
    if(n<1.or.size(values,2)/=size(weights).or.&
        any(shape(gradients)/=[3,n,size(weights)]).or.any(weights<=0d0).or.&
        any(.not.ieee_is_finite(weights)).or.any(.not.ieee_is_finite(real(values))).or.&
        any(.not.ieee_is_finite(aimag(values))).or.any(.not.ieee_is_finite(real(gradients))).or.&
        any(.not.ieee_is_finite(aimag(gradients))))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.global_bad/=0)then
      callback_message='invalid canonical-momentum construction payload';return
    endif
    call MPI_Allreduce(n,minimum_n,1,MPI_INTEGER,MPI_MIN,comm_arg,ierr_local)
    if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(n,maximum_n,1,MPI_INTEGER,MPI_MAX,&
      comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.minimum_n/=maximum_n)then
      callback_message='inconsistent canonical-momentum basis rank';return
    endif
    allocate(local_derivative(3,n,n),derivative(3,n,n),momentum(3,n,n))
    local_derivative=(0d0,0d0)
    do point=1,size(weights);do column=1,n;do row=1,n;do axis=1,3
      local_derivative(axis,row,column)=local_derivative(axis,row,column)+weights(point)*&
        conjg(values(row,point))*gradients(axis,column,point)
    enddo;enddo;enddo;enddo
    call MPI_Allreduce(local_derivative,derivative,3*n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,&
      comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS)then
      callback_message='canonical-momentum reduction failed';return
    endif
    do axis=1,3
      momentum(axis,:,:)=-cmplx(0d0,1d0,8)*0.5d0*&
        (derivative(axis,:,:)-conjg(transpose(derivative(axis,:,:))))
    enddo
    if(any(.not.ieee_is_finite(real(momentum))).or.any(.not.ieee_is_finite(aimag(momentum))))then
      callback_message='canonical-momentum projection is nonfinite';return
    endif
    callback_ok=.true.
  end subroutine assemble_dg_hybrid_canonical_momentum

  subroutine project_dg_hybrid_operator_to_rt(basis,operator,projected)
    complex(8),intent(in)::basis(:,:),operator(:,:)
    complex(8),allocatable,intent(out)::projected(:,:)
    allocate(projected(size(basis,2),size(basis,2)))
    projected=matmul(conjg(transpose(basis)),matmul(operator,basis))
  end subroutine project_dg_hybrid_operator_to_rt

  integer(8) function checkpoint_integer64_fingerprint(values) result(hash)
    integer(8),intent(in)::values(:)
    integer::p
    hash=int(z'6A09E667F3BCC909',8)
    hash=ieor(ishftc(hash,7),int(size(values),8))
    do p=1,size(values);hash=ieor(ishftc(hash,11),ieor(values(p),int(p,8)));enddo
    if(hash==0_8)hash=1_8
  end function checkpoint_integer64_fingerprint

  integer(8) function checkpoint_integer_fingerprint(values) result(hash)
    integer,intent(in)::values(:)
    integer::p
    hash=int(z'BB67AE8584CAA73B',8)
    hash=ieor(ishftc(hash,7),int(size(values),8))
    do p=1,size(values);hash=ieor(ishftc(hash,11),ieor(int(values(p),8),int(p,8)));enddo
    if(hash==0_8)hash=1_8
  end function checkpoint_integer_fingerprint

  integer(8) function checkpoint_complex_matrix_fingerprint(values) result(hash)
    complex(8),intent(in)::values(:,:)
    integer::p,q
    integer(8)::bits
    hash=int(z'3C6EF372FE94F82B',8)
    hash=ieor(ishftc(hash,7),int(size(values,1),8));hash=ieor(ishftc(hash,7),int(size(values,2),8))
    do q=1,size(values,2);do p=1,size(values,1)
      bits=transfer(real(values(p,q)),bits);hash=ieor(ishftc(hash,11),bits)
      bits=transfer(aimag(values(p,q)),bits);hash=ieor(ishftc(hash,13),bits)
    enddo;enddo
    if(hash==0_8)hash=1_8
  end function checkpoint_complex_matrix_fingerprint

  integer(8) function checkpoint_complex_rank4_fingerprint(values) result(hash)
    complex(8),intent(in)::values(:,:,:,:)
    integer::p,q,r,s
    integer(8)::bits
    hash=int(z'A54FF53A5F1D36F1',8)
    hash=ieor(hash,int(size(values),8));hash=ieor(hash,ishftc(int(size(values,4),8),17))
    do s=1,size(values,4);do r=1,size(values,3);do q=1,size(values,2);do p=1,size(values,1)
      bits=transfer(real(values(p,q,r,s)),bits);hash=ieor(ishftc(hash,7),bits)
      bits=transfer(aimag(values(p,q,r,s)),bits);hash=ieor(ishftc(hash,11),bits)
    enddo;enddo;enddo;enddo
    if(hash==0_8)hash=1_8
  end function checkpoint_complex_rank4_fingerprint

  integer(8) function checkpoint_complex_rank5_fingerprint(values) result(hash)
    complex(8),intent(in)::values(:,:,:,:,:)
    integer::p,q,r,s,t
    integer(8)::bits
    hash=int(z'510E527FADE682D1',8)
    hash=ieor(hash,int(size(values),8));hash=ieor(hash,ishftc(int(size(values,5),8),17))
    do t=1,size(values,5);do s=1,size(values,4);do r=1,size(values,3)
      do q=1,size(values,2);do p=1,size(values,1)
        bits=transfer(real(values(p,q,r,s,t)),bits);hash=ieor(ishftc(hash,7),bits)
        bits=transfer(aimag(values(p,q,r,s,t)),bits);hash=ieor(ishftc(hash,11),bits)
      enddo;enddo
    enddo;enddo;enddo
    if(hash==0_8)hash=1_8
  end function checkpoint_complex_rank5_fingerprint

  integer(8) function checkpoint_real_rank3_fingerprint(values) result(hash)
    real(8),intent(in)::values(:,:,:)
    integer::p,q,r
    integer(8)::bits
    hash=int(z'1F83D9ABFB41BD6B',8)
    hash=ieor(ishftc(hash,7),int(size(values,1),8))
    hash=ieor(ishftc(hash,7),int(size(values,2),8))
    hash=ieor(ishftc(hash,7),int(size(values,3),8))
    do r=1,size(values,3);do q=1,size(values,2);do p=1,size(values,1)
      bits=transfer(values(p,q,r),bits);hash=ieor(ishftc(hash,11),bits)
    enddo;enddo;enddo
    if(hash==0_8)hash=1_8
  end function checkpoint_real_rank3_fingerprint

  subroutine checkpoint_grid_complex_fingerprint(comm_arg,grid_ids_arg,values,fingerprint,callback_ok)
    integer,intent(in)::comm_arg
    integer(8),intent(in)::grid_ids_arg(:)
    complex(8),intent(in)::values(:,:)
    integer(8),intent(out)::fingerprint
    logical,intent(out)::callback_ok
    integer::p,q,ierr_local,local_bad,global_bad,local_rows,minimum_rows,maximum_rows
    integer(8)::local_hash,bits,local_count,global_count
    callback_ok=.false.;fingerprint=0_8
    local_bad=merge(0,1,size(values,2)==size(grid_ids_arg))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.global_bad/=0)return
    local_rows=size(values,1)
    call MPI_Allreduce(local_rows,minimum_rows,1,MPI_INTEGER,MPI_MIN,comm_arg,ierr_local)
    if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(local_rows,maximum_rows,1,MPI_INTEGER,MPI_MAX,&
      comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.minimum_rows/=maximum_rows)return
    local_hash=0_8
    do p=1,size(grid_ids_arg);do q=1,size(values,1)
      bits=transfer(real(values(q,p)),bits)
      local_hash=ieor(local_hash,ishftc(ieor(bits,grid_ids_arg(p)),mod(7*q,63)))
      bits=transfer(aimag(values(q,p)),bits)
      local_hash=ieor(local_hash,ishftc(ieor(bits,ishftc(grid_ids_arg(p),17)),mod(11*q,63)))
    enddo;enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm_arg,ierr_local)
    local_count=int(size(grid_ids_arg),8);global_count=0_8
    if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER8,MPI_SUM,&
      comm_arg,ierr_local)
    fingerprint=ieor(fingerprint,ishftc(int(size(values,1),8),31))
    fingerprint=ieor(fingerprint,ishftc(global_count,23))
    if(fingerprint==0_8)fingerprint=1_8
    callback_ok=ierr_local==MPI_SUCCESS
  end subroutine checkpoint_grid_complex_fingerprint

  subroutine checkpoint_grid_real_fingerprint(comm_arg,grid_ids_arg,values,fingerprint,callback_ok)
    integer,intent(in)::comm_arg
    integer(8),intent(in)::grid_ids_arg(:)
    real(8),intent(in)::values(:)
    integer(8),intent(out)::fingerprint
    logical,intent(out)::callback_ok
    integer::p,ierr_local,local_bad,global_bad
    integer(8)::local_hash,bits,local_count,global_count
    callback_ok=.false.;fingerprint=0_8
    local_bad=merge(0,1,size(values)==size(grid_ids_arg))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.global_bad/=0)return
    local_hash=0_8
    do p=1,size(values)
      bits=transfer(values(p),bits)
      local_hash=ieor(local_hash,ishftc(ieor(bits,ishftc(grid_ids_arg(p),17)),&
        int(modulo(grid_ids_arg(p),63_8))))
    enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm_arg,ierr_local)
    local_count=int(size(values),8);global_count=0_8
    if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER8,MPI_SUM,&
      comm_arg,ierr_local)
    fingerprint=ieor(fingerprint,ishftc(global_count,31))
    if(fingerprint==0_8)fingerprint=1_8
    callback_ok=ierr_local==MPI_SUCCESS
  end subroutine checkpoint_grid_real_fingerprint

  subroutine checkpoint_grid_integer_fingerprint(comm_arg,grid_ids_arg,values,fingerprint,callback_ok)
    integer,intent(in)::comm_arg,values(:)
    integer(8),intent(in)::grid_ids_arg(:)
    integer(8),intent(out)::fingerprint
    logical,intent(out)::callback_ok
    integer::p,ierr_local,local_bad,global_bad
    integer(8)::local_hash,entry,local_count,global_count
    callback_ok=.false.;fingerprint=0_8
    local_bad=merge(0,1,size(values)==size(grid_ids_arg).and.all(values>0))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.global_bad/=0)return
    local_hash=0_8
    do p=1,size(values)
      entry=ieor(grid_ids_arg(p),ishftc(int(values(p),8),17))
      local_hash=ieor(local_hash,ishftc(entry,int(modulo(grid_ids_arg(p),63_8))))
    enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm_arg,ierr_local)
    local_count=int(size(values),8);global_count=0_8
    if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER8,MPI_SUM,&
      comm_arg,ierr_local)
    fingerprint=ieor(fingerprint,ishftc(global_count,31))
    if(fingerprint==0_8)fingerprint=1_8
    callback_ok=ierr_local==MPI_SUCCESS
  end subroutine checkpoint_grid_integer_fingerprint

  subroutine pack_checkpoint_face_values(face,values,position)
    type(s_dg_hybrid_production_face_trace),intent(in)::face
    complex(8),intent(inout)::values(:,:)
    integer,intent(inout)::position
    integer::count
    count=size(face%value_minus)
    if(count>0)values(1,position+1:position+count)=reshape(face%value_minus,[count]);position=position+count
    count=size(face%value_plus)
    if(count>0)values(1,position+1:position+count)=reshape(face%value_plus,[count]);position=position+count
    count=size(face%derivative_minus)
    if(count>0)values(1,position+1:position+count)=reshape(face%derivative_minus,[count]);position=position+count
    count=size(face%derivative_plus)
    if(count>0)values(1,position+1:position+count)=reshape(face%derivative_plus,[count]);position=position+count
  end subroutine pack_checkpoint_face_values

  subroutine build_checkpoint_topology_graphs(global_ids,row_ids,basis_fragment,faces,&
      metric_offsets,metric_columns,operator_offsets,operator_columns)
    integer,intent(in)::global_ids(:),basis_fragment(:)
    integer(8),intent(in)::row_ids(:)
    type(s_dg_hybrid_production_face_trace),intent(in)::faces(:)
    integer,allocatable,intent(out)::metric_offsets(:),metric_columns(:),operator_offsets(:),operator_columns(:)
    integer::row,column,metric_count,operator_count,row_position
    ! The broken-volume metric has support only within a fragment.  The
    ! nonlocal operator contract permits every retained basis pair, so the
    ! shared Hamiltonian envelope is deliberately dense and retains explicit
    ! zero edges from kinetic, local, and SIPG components.
    metric_count=0
    do row=1,size(row_ids)
      row_position=findloc(global_ids,int(row_ids(row)),dim=1)
      metric_count=metric_count+count(basis_fragment==basis_fragment(row_position))
    enddo
    operator_count=size(row_ids)*size(global_ids)
    allocate(metric_offsets(size(row_ids)+1),metric_columns(metric_count),&
      operator_offsets(size(row_ids)+1),operator_columns(operator_count))
    metric_count=0;operator_count=0;metric_offsets(1)=1;operator_offsets(1)=1
    do row=1,size(row_ids)
      row_position=findloc(global_ids,int(row_ids(row)),dim=1)
      do column=1,size(global_ids)
        if(basis_fragment(column)==basis_fragment(row_position))then
          metric_count=metric_count+1;metric_columns(metric_count)=global_ids(column)
        endif
        operator_count=operator_count+1;operator_columns(operator_count)=global_ids(column)
      enddo
      metric_offsets(row+1)=metric_count+1;operator_offsets(row+1)=operator_count+1
    enddo
  end subroutine build_checkpoint_topology_graphs

  subroutine selection_added_members(requested,effective,added)
    integer,intent(in)::requested(:),effective(:)
    integer,allocatable,intent(out)::added(:)
    integer::i,position
    allocate(added(count([(count(requested==effective(i))==0,i=1,size(effective))])))
    position=0
    do i=1,size(effective)
      if(any(requested==effective(i)))cycle
      position=position+1;added(position)=effective(i)
    enddo
  end subroutine selection_added_members

  integer(8) function checkpoint_analysis_fingerprint(representation) result(fingerprint)
    complex(8),intent(in)::representation(:,:,:)
    integer::i,j,k;integer(8)::bits
    fingerprint=int(size(representation,3),8)
    do k=1,size(representation,3);do j=1,size(representation,2);do i=1,size(representation,1)
      bits=transfer(real(representation(i,j,k),8),bits);fingerprint=ieor(ishftc(fingerprint,9),bits)
      bits=transfer(aimag(representation(i,j,k)),bits);fingerprint=ieor(ishftc(fingerprint,9),bits)
    enddo;enddo;enddo
    if(fingerprint==0_8)fingerprint=1_8
  end function checkpoint_analysis_fingerprint

  integer(8) function checkpoint_real_fingerprint(values) result(fingerprint)
    real(8),intent(in)::values(:)
    integer(8)::bits
    integer::i
    fingerprint=int(z'9E3779B97F4A7C15',8)
    do i=1,size(values)
      bits=transfer(values(i),bits)
      fingerprint=ieor(ishftc(fingerprint,11),ieor(bits,int(i,8)))
    enddo
    if(fingerprint==0_8)fingerprint=1_8
  end function checkpoint_real_fingerprint

  subroutine set_dg_hybrid_trial_state(state,density,potential,occupations_arg,eigenvalues_arg,&
      projector,trace,epoch,operator_fingerprint,solver_fingerprint)
    type(s_dg_hybrid_trial_state),intent(out)::state
    real(8),intent(in)::density(:),potential(:),occupations_arg(:),eigenvalues_arg(:)
    complex(8),intent(in)::projector(:,:),trace(:,:)
    integer,intent(in)::epoch
    integer(8),intent(in)::operator_fingerprint,solver_fingerprint
    allocate(state%density,source=density);allocate(state%potential,source=potential)
    allocate(state%occupations,source=occupations_arg);allocate(state%eigenvalues,source=eigenvalues_arg)
    allocate(state%mixing_history,source=density);allocate(state%projector,source=projector)
    allocate(state%trace,source=trace)
    state%density_epoch=epoch;state%operator_epoch=epoch;state%projector_epoch=epoch
    state%trace_epoch=epoch;state%derived_epoch=epoch
    state%operator_structure_fingerprint=operator_fingerprint
    state%operator_value_fingerprint=ieor(operator_fingerprint,solver_fingerprint)
    if(state%operator_value_fingerprint==0_8)state%operator_value_fingerprint=1_8
    state%trace_cache_valid=.true.
  end subroutine set_dg_hybrid_trial_state

  subroutine extract_dg_hybrid_core_local_potential(core_ids,core_potential,callback_ok)
    integer(8),intent(in)::core_ids(:)
    real(8),intent(out)::core_potential(:)
    logical,intent(out)::callback_ok
    integer::p,gx,gy,gz,ix,iy,iz
    callback_ok=.false.
    if(size(core_potential)/=size(core_ids))return
    do p=1,size(core_ids)
      gx=int(modulo(core_ids(p)-1_8,int(dc%lg_tot%num(1),8)))+1
      gy=int(modulo((core_ids(p)-1_8)/int(dc%lg_tot%num(1),8),int(dc%lg_tot%num(2),8)))+1
      gz=int((core_ids(p)-1_8)/int(dc%lg_tot%num(1)*dc%lg_tot%num(2),8))+1
      ix=findloc(dc%jxyz_tot(:,1),gx,dim=1);iy=findloc(dc%jxyz_tot(:,2),gy,dim=1)
      iz=findloc(dc%jxyz_tot(:,3),gz,dim=1)
      if(ix<1.or.iy<1.or.iz<1)return
      core_potential(p)=v_local(1)%f(ix,iy,iz)
    enddo
    callback_ok=all(ieee_is_finite(core_potential))
  end subroutine extract_dg_hybrid_core_local_potential

  subroutine form_dg_hybrid_coefficient_actions(comm_arg,row_ids_arg,hrows_arg,srows_arg,&
      coefficients_arg,eigenvalues_arg,hc_arg,sc_epsilon_arg,callback_ok)
    integer,intent(in)::comm_arg
    integer(8),intent(in)::row_ids_arg(:)
    complex(8),intent(in)::hrows_arg(:,:),srows_arg(:,:),coefficients_arg(:,:)
    real(8),intent(in)::eigenvalues_arg(:)
    complex(8),allocatable,intent(out)::hc_arg(:,:),sc_epsilon_arg(:,:)
    logical,intent(out)::callback_ok
    complex(8),allocatable::local_coefficients(:,:),global_coefficients(:,:)
    integer::i,ierr_local,global_count
    global_count=size(hrows_arg,2);callback_ok=.false.
    allocate(local_coefficients(global_count,size(coefficients_arg,2)),&
      global_coefficients(global_count,size(coefficients_arg,2)))
    local_coefficients=(0d0,0d0)
    do i=1,size(row_ids_arg);local_coefficients(int(row_ids_arg(i)),:)=coefficients_arg(i,:);enddo
    call MPI_Allreduce(local_coefficients,global_coefficients,size(global_coefficients),&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS)return
    allocate(hc_arg(size(coefficients_arg,1),size(coefficients_arg,2)),&
      sc_epsilon_arg(size(coefficients_arg,1),size(coefficients_arg,2)))
    hc_arg=matmul(hrows_arg,global_coefficients);sc_epsilon_arg=matmul(srows_arg,global_coefficients)
    do i=1,size(eigenvalues_arg);sc_epsilon_arg(:,i)=eigenvalues_arg(i)*sc_epsilon_arg(:,i);enddo
    callback_ok=all(ieee_is_finite(real(hc_arg))).and.all(ieee_is_finite(aimag(hc_arg))).and.&
      all(ieee_is_finite(real(sc_epsilon_arg))).and.all(ieee_is_finite(aimag(sc_epsilon_arg)))
  end subroutine form_dg_hybrid_coefficient_actions

  subroutine measure_dg_hybrid_operator_covariance(comm_arg,row_ids_arg,hrows_arg,representation_arg,&
      residual,callback_ok)
    integer,intent(in)::comm_arg
    integer(8),intent(in)::row_ids_arg(:)
    complex(8),intent(in)::hrows_arg(:,:),representation_arg(:,:,:)
    real(8),intent(out)::residual
    logical,intent(out)::callback_ok
    integer::rank_local,nproc_local,ierr_local,nowned,global_count,noperation,r,operation,i,j,k,offset,nrows
    integer,allocatable::counts(:),displacements(:)
    integer(8),allocatable::all_rows(:)
    complex(8),allocatable::h_times_d(:,:),partial(:,:),reduced(:,:)
    real(8)::local_defect,global_defect,local_scale,global_scale
    callback_ok=.false.;residual=huge(1d0);nowned=size(row_ids_arg);global_count=size(hrows_arg,2)
    noperation=size(representation_arg,3)
    if(any(shape(hrows_arg)/=[nowned,global_count]).or.&
        any(shape(representation_arg)/=[global_count,global_count,noperation]).or.noperation<1)return
    call MPI_Comm_rank(comm_arg,rank_local,ierr_local);if(ierr_local/=MPI_SUCCESS)return
    call MPI_Comm_size(comm_arg,nproc_local,ierr_local);if(ierr_local/=MPI_SUCCESS)return
    allocate(counts(nproc_local),displacements(nproc_local))
    call MPI_Allgather(nowned,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm_arg,ierr_local)
    displacements(1)=0
    do r=2,nproc_local;displacements(r)=displacements(r-1)+counts(r-1);enddo
    allocate(all_rows(sum(counts)))
    call MPI_Allgatherv(row_ids_arg,nowned,MPI_INTEGER8,all_rows,counts,displacements,MPI_INTEGER8,&
      comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.size(all_rows)/=global_count)return
    local_defect=0d0;local_scale=max(1d0,maxval(abs(hrows_arg)))
    do operation=1,noperation
      allocate(h_times_d(nowned,global_count))
      h_times_d=matmul(hrows_arg,representation_arg(:,:,operation))
      do r=0,nproc_local-1
        nrows=counts(r+1);offset=displacements(r+1)
        allocate(partial(nrows,global_count),reduced(nrows,global_count));partial=(0d0,0d0)
        do i=1,nrows
          do k=1,nowned
            partial(i,:)=partial(i,:)+conjg(representation_arg(int(row_ids_arg(k)),&
              int(all_rows(offset+i)),operation))*h_times_d(k,:)
          enddo
        enddo
        call MPI_Reduce(partial,reduced,nrows*global_count,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm_arg,ierr_local)
        if(ierr_local/=MPI_SUCCESS)return
        if(rank_local==r.and.nrows>0)then
          local_defect=max(local_defect,maxval(abs(reduced-hrows_arg)))
          local_scale=max(local_scale,maxval(abs(reduced)))
        endif
        deallocate(partial,reduced)
      enddo
      deallocate(h_times_d)
    enddo
    call MPI_Allreduce(local_defect,global_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm_arg,ierr_local)
    if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
      comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS)return
    residual=global_defect/max(1d0,global_scale);callback_ok=ieee_is_finite(residual)
  end subroutine measure_dg_hybrid_operator_covariance

  subroutine measure_dg_hybrid_projector_covariance(comm_arg,row_ids_arg,projector_rows_arg,&
      representation_arg,residual,callback_ok)
    integer,intent(in)::comm_arg
    integer(8),intent(in)::row_ids_arg(:)
    complex(8),intent(in)::projector_rows_arg(:,:),representation_arg(:,:,:)
    real(8),intent(out)::residual
    logical,intent(out)::callback_ok
    integer::rank_local,nproc_local,ierr_local,nowned,global_count,noperation,r,operation,i,k,offset,nrows
    integer,allocatable::counts(:),displacements(:)
    integer(8),allocatable::all_rows(:)
    complex(8),allocatable::q_times_d(:,:),partial(:,:),reduced(:,:)
    real(8)::local_defect,global_defect,local_scale,global_scale
    callback_ok=.false.;residual=huge(1d0);nowned=size(row_ids_arg);global_count=size(projector_rows_arg,2)
    noperation=size(representation_arg,3)
    if(any(shape(projector_rows_arg)/=[nowned,global_count]).or.&
        any(shape(representation_arg)/=[global_count,global_count,noperation]).or.noperation<1)return
    call MPI_Comm_rank(comm_arg,rank_local,ierr_local);if(ierr_local/=MPI_SUCCESS)return
    call MPI_Comm_size(comm_arg,nproc_local,ierr_local);if(ierr_local/=MPI_SUCCESS)return
    allocate(counts(nproc_local),displacements(nproc_local))
    call MPI_Allgather(nowned,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm_arg,ierr_local)
    displacements(1)=0
    do r=2,nproc_local;displacements(r)=displacements(r-1)+counts(r-1);enddo
    allocate(all_rows(sum(counts)))
    call MPI_Allgatherv(row_ids_arg,nowned,MPI_INTEGER8,all_rows,counts,displacements,MPI_INTEGER8,&
      comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.size(all_rows)/=global_count)return
    local_defect=0d0;local_scale=max(1d0,maxval(abs(projector_rows_arg)))
    do operation=1,noperation
      allocate(q_times_d(nowned,global_count))
      q_times_d=matmul(projector_rows_arg,representation_arg(:,:,operation))
      do r=0,nproc_local-1
        nrows=counts(r+1);offset=displacements(r+1)
        allocate(partial(nrows,global_count),reduced(nrows,global_count));partial=(0d0,0d0)
        do i=1,nrows
          do k=1,nowned
            partial(i,:)=partial(i,:)+representation_arg(int(all_rows(offset+i)),&
              int(row_ids_arg(k)),operation)*projector_rows_arg(k,:)
          enddo
        enddo
        call MPI_Reduce(partial,reduced,nrows*global_count,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm_arg,ierr_local)
        if(ierr_local/=MPI_SUCCESS)return
        if(rank_local==r.and.nrows>0)then
          local_defect=max(local_defect,maxval(abs(reduced-q_times_d)))
          local_scale=max(local_scale,maxval(abs(reduced)),maxval(abs(q_times_d)))
        endif
        deallocate(partial,reduced)
      enddo
      deallocate(q_times_d)
    enddo
    call MPI_Allreduce(local_defect,global_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm_arg,ierr_local)
    if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
      comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS)return
    residual=global_defect/max(1d0,global_scale);callback_ok=ieee_is_finite(residual)
  end subroutine measure_dg_hybrid_projector_covariance

  subroutine evaluate_dg_hybrid_distributed_low_energy_symmetry(comm_arg,row_ids_arg,metric_rows_arg,&
      representation_arg,coefficients_arg,eigenvalues_arg,occupied_rank_arg,target_rank_arg,tolerance_arg,&
      occupied_defect_arg,target_defect_arg,energy_defect_arg,callback_ok,message_arg)
    integer,intent(in)::comm_arg,occupied_rank_arg,target_rank_arg
    integer(8),intent(in)::row_ids_arg(:)
    complex(8),intent(in)::metric_rows_arg(:,:),representation_arg(:,:,:),coefficients_arg(:,:)
    real(8),intent(in)::eigenvalues_arg(:),tolerance_arg
    real(8),intent(out)::occupied_defect_arg,target_defect_arg,energy_defect_arg
    logical,intent(out)::callback_ok
    character(*),intent(out)::message_arg
    integer::global_count,state_count,operation_count,i,row,ierr_local,local_bad,global_bad,allocation_status
    integer::local_contract(5),minimum_contract(5),maximum_contract(5)
    integer(8)::tolerance_bits,minimum_tolerance_bits,maximum_tolerance_bits
    integer,allocatable::ownership_count(:)
    complex(8),allocatable::global_metric(:,:),global_coefficients(:,:)
    logical::symmetry_ok

    callback_ok=.false.;message_arg='';occupied_defect_arg=huge(1d0)
    target_defect_arg=huge(1d0);energy_defect_arg=huge(1d0)
    global_count=size(metric_rows_arg,2);state_count=size(coefficients_arg,2)
    operation_count=size(representation_arg,3);local_bad=0
    local_contract=[global_count,state_count,operation_count,occupied_rank_arg,target_rank_arg]
    call MPI_Allreduce(local_contract,minimum_contract,5,MPI_INTEGER,MPI_MIN,comm_arg,ierr_local)
    if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(local_contract,maximum_contract,5,MPI_INTEGER,MPI_MAX,&
      comm_arg,ierr_local)
    tolerance_bits=transfer(tolerance_arg,tolerance_bits)
    if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(tolerance_bits,minimum_tolerance_bits,1,MPI_INTEGER8,MPI_MIN,&
      comm_arg,ierr_local)
    if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(tolerance_bits,maximum_tolerance_bits,1,MPI_INTEGER8,MPI_MAX,&
      comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS)then
      message_arg='distributed low-energy symmetry agreement failed';return
    endif
    if(any(minimum_contract/=maximum_contract).or.minimum_tolerance_bits/=maximum_tolerance_bits)then
      message_arg='rank-disagreeing distributed low-energy symmetry contract';return
    endif
    if(global_count<1.or.state_count<1.or.operation_count<1.or.&
        size(metric_rows_arg,1)/=size(row_ids_arg).or.size(coefficients_arg,1)/=size(row_ids_arg).or.&
        size(eigenvalues_arg)/=state_count.or.&
        any(shape(representation_arg)/=[global_count,global_count,operation_count]).or.&
        occupied_rank_arg<1.or.occupied_rank_arg>target_rank_arg.or.target_rank_arg>state_count.or.&
        any(row_ids_arg<1_8).or.any(row_ids_arg>int(global_count,8)))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.global_bad/=0)then
      message_arg='invalid distributed low-energy symmetry contract';return
    endif
    allocate(ownership_count(global_count),global_metric(global_count,global_count),&
      global_coefficients(global_count,state_count),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.global_bad/=0)then
      message_arg='distributed low-energy symmetry allocation failed';return
    endif
    ownership_count=0;global_metric=(0d0,0d0);global_coefficients=(0d0,0d0)
    do i=1,size(row_ids_arg)
      row=int(row_ids_arg(i));ownership_count(row)=ownership_count(row)+1
      global_metric(row,:)=metric_rows_arg(i,:);global_coefficients(row,:)=coefficients_arg(i,:)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_count,MPI_INTEGER,MPI_SUM,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.any(ownership_count/=1))then
      message_arg='distributed low-energy symmetry rows are not owned exactly once';return
    endif
    call MPI_Allreduce(MPI_IN_PLACE,global_metric,size(global_metric),MPI_DOUBLE_COMPLEX,MPI_SUM,&
      comm_arg,ierr_local)
    if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(MPI_IN_PLACE,global_coefficients,size(global_coefficients),&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS)then
      message_arg='distributed low-energy symmetry row collection failed';return
    endif
    call evaluate_dg_hybrid_low_energy_symmetry(comm_arg,global_metric,representation_arg,global_coefficients,&
      eigenvalues_arg,occupied_rank_arg,target_rank_arg,tolerance_arg,occupied_defect_arg,target_defect_arg,&
      energy_defect_arg,symmetry_ok,message_arg)
    if(symmetry_ok)then
      callback_ok=.true.;return
    endif
    if(trim(message_arg)=='low-energy LCFO eigenspace is not symmetry closed'.and.&
        all(ieee_is_finite([occupied_defect_arg,target_defect_arg,energy_defect_arg])))then
      callback_ok=.true.;return
    endif
  end subroutine evaluate_dg_hybrid_distributed_low_energy_symmetry

  subroutine measure_dg_hybrid_core_density_covariance(comm_arg,density_arg,maps_arg,&
      residual_arg,callback_ok,message_arg)
    integer,intent(in)::comm_arg
    real(8),intent(in)::density_arg(:)
    integer(8),intent(in)::maps_arg(:,:)
    real(8),intent(out)::residual_arg
    logical,intent(out)::callback_ok
    character(*),intent(out)::message_arg
    complex(8),allocatable::density_probe(:,:),mapped_density(:,:)
    real(8)::local_values(2),global_values(2)
    integer::operation_count,minimum_operation_count,maximum_operation_count,operation,&
      ierr_local,local_bad,global_bad,allocation_status

    callback_ok=.false.;message_arg='';residual_arg=huge(1d0);operation_count=size(maps_arg,2)
    call MPI_Allreduce(operation_count,minimum_operation_count,1,MPI_INTEGER,MPI_MIN,comm_arg,ierr_local)
    if(ierr_local==MPI_SUCCESS)call MPI_Allreduce(operation_count,maximum_operation_count,1,MPI_INTEGER,MPI_MAX,&
      comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS)then
      message_arg='mapped-core density operation agreement failed';return
    endif
    if(minimum_operation_count/=maximum_operation_count)then
      message_arg='rank-disagreeing mapped-core density operation count';return
    endif
    local_bad=merge(0,1,size(density_arg)>0.and.size(maps_arg,1)==size(density_arg).and.&
      operation_count>0.and.all(ieee_is_finite(density_arg)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.global_bad/=0)then
      message_arg='invalid mapped-core density covariance contract';return
    endif
    allocate(density_probe(1,size(density_arg)),mapped_density(1,size(density_arg)),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr_local)
    if(ierr_local/=MPI_SUCCESS.or.global_bad/=0)then
      message_arg='mapped-core density covariance allocation failed';return
    endif
    density_probe(1,:)=cmplx(density_arg,0d0,8);residual_arg=0d0
    do operation=1,operation_count
      call exchange_dg_point_permuted_orbital_rows(comm_arg,density_probe,maps_arg(:,operation),&
        mapped_density,callback_ok,message_arg)
      if(.not.callback_ok)return
      local_values=[maxval(abs(mapped_density-density_probe)),&
        max(maxval(abs(mapped_density)),maxval(abs(density_probe)))]
      call MPI_Allreduce(local_values,global_values,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm_arg,ierr_local)
      if(ierr_local/=MPI_SUCCESS)then
        callback_ok=.false.;message_arg='mapped-core density covariance reduction failed';return
      endif
      residual_arg=max(residual_arg,global_values(1)/max(tiny(1d0),global_values(2)))
    enddo
    callback_ok=ieee_is_finite(residual_arg)
    if(.not.callback_ok)message_arg='mapped-core density covariance defect is not finite'
  end subroutine measure_dg_hybrid_core_density_covariance

  subroutine assemble_dg_hybrid_divided_core_density(core_density,electron_count,callback_ok)
    real(8),intent(out)::core_density(:),electron_count
    logical,intent(out)::callback_ok
    integer::p,position,ierr_local,local_bad,global_bad
    real(8)::local_electron_count

    callback_ok=.false.;core_density=0d0;electron_count=0d0
    local_bad=0
    if(.not.allocated(divided_fragment_density).or.size(core_density)/=size(ow_core_ids))then
      local_bad=1
    else
      do p=1,size(ow_core_ids)
        position=findloc(divided_fragment_basis%buffer_point_ids,ow_core_ids(p),dim=1)
        if(position<1)then
          local_bad=1
        else
          core_density(p)=divided_fragment_density(position)
        endif
      enddo
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr_local)
    if(ierr_local/=MPI_SUCCESS)return
    if(global_bad/=0)return
    local_electron_count=sum(ow_core_weights*core_density)
    call MPI_Allreduce(local_electron_count,electron_count,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
      dc%icomm_tot,ierr_local)
    callback_ok=ierr_local==MPI_SUCCESS.and.all(ieee_is_finite(core_density)).and.&
      ieee_is_finite(electron_count)
  end subroutine assemble_dg_hybrid_divided_core_density

  subroutine mix_dg_hybrid_divided_density(iteration,input_density,new_density,mixed_density,callback_ok)
    integer,intent(in)::iteration
    real(8),intent(in)::input_density(:),new_density(:)
    real(8),intent(out)::mixed_density(:)
    logical,intent(out)::callback_ok
    real(8),allocatable::input_total(:,:,:),new_total(:,:,:)
    integer::p,ix_local,iy_local,iz_local,mixing_iteration,mixing_basis_generation
    integer(8)::mixing_inventory_fingerprint
    character(16)::selected_mixing_method
    character(64)::reset_reason
    character(512)::mixing_message
    real(8)::effective_mixrate
    logical::reset_required,mixing_ok

    callback_ok=.false.;mixed_density=0d0
    call gather_dg_hybrid_divided_core_density(input_density,input_total,callback_ok)
    if(.not.callback_ok)return
    call gather_dg_hybrid_divided_core_density(new_density,new_total,callback_ok)
    if(.not.callback_ok)return
    mixing_basis_generation=divided_mixing_basis_generation
    mixing_inventory_fingerprint=divided_mixing_inventory_fingerprint
    call prepare_dg_hybrid_divided_mixing(dc%icomm_tot,dg_hybrid_divided_mixing,method_mixing,&
      mixing%mixrate,mixing_basis_generation,mixing_inventory_fingerprint,&
      divided_mixing_rollback_pending,divided_mixing_state,selected_mixing_method,&
      effective_mixrate,mixing_iteration,reset_required,reset_reason,mixing_ok,mixing_message)
    if(.not.mixing_ok)then
      if(dc%id_tot==0)write(0,'(a)')trim(mixing_message)
      return
    endif
    if(dc%id_tot==0)write(*,'(a,a,a,es16.8,a,i0,a,l1,a,a,a)')&
      '[DG-HYBRID-MIXING] method=',trim(selected_mixing_method),' mixrate=',effective_mixrate,&
      ' history_length=',divided_mixing_state%history_length,' reset=',reset_required,&
      ' reset_reason=',trim(reset_reason),' seed_fingerprint_unchanged=T'
    dc%rho_tot_s(1)%f=input_total
    call copy_density(mixing_iteration,dc%system_tot%nspin,dc%mg_tot,dc%rho_tot_s,mixing)
    dc%rho_tot_s(1)%f=new_total
    select case(selected_mixing_method)
    case('simple')
      call simple_mixing(dc%mg_tot,dc%system_tot,1d0-effective_mixrate,effective_mixrate,dc%rho_tot_s,mixing)
    case('broyden')
      call wrapper_broyden(dc%info_tot%icomm_r,dc%mg_tot,dc%system_tot,dc%rho_tot_s,mixing_iteration,mixing)
    case('pulay')
      call pulay(dc%mg_tot,dc%info_tot,dc%system_tot,dc%rho_tot_s,mixing_iteration,mixing)
    case default
      return
    end select
    do p=1,size(ow_core_ids)
      ix_local=int(modulo(ow_core_ids(p)-1_8,int(dc%lg_tot%num(1),8)))+1
      iy_local=int(modulo((ow_core_ids(p)-1_8)/int(dc%lg_tot%num(1),8),int(dc%lg_tot%num(2),8)))+1
      iz_local=int((ow_core_ids(p)-1_8)/(int(dc%lg_tot%num(1),8)*int(dc%lg_tot%num(2),8)))+1
      mixed_density(p)=dc%rho_tot_s(1)%f(ix_local,iy_local,iz_local)
    enddo
    callback_ok=all(ieee_is_finite(mixed_density))
    call accept_dg_hybrid_divided_mixing(dc%icomm_tot,mixing_iteration,divided_mixing_state,&
      mixing_ok,mixing_message,local_accept_ok=callback_ok)
    if(.not.mixing_ok.and.dc%id_tot==0)write(0,'(a)')trim(mixing_message)
    callback_ok=mixing_ok
    divided_mixing_rollback_pending=.not.mixing_ok
  end subroutine mix_dg_hybrid_divided_density

  subroutine solve_final_dg_hybrid_divided_lcfo(comm_arg,global_count_arg,nstate_arg,row_ids_arg,&
      hrows_arg,srows_arg,tolerance_arg,coefficients_arg,eigenvalues_arg,maximum_residual_arg,&
      orthogonality_defect_arg,projector_defect_arg,workspace_peak_bytes_arg,fingerprint_arg,&
      callback_ok,callback_message)
    integer,intent(in)::comm_arg,global_count_arg,nstate_arg
    integer(8),intent(in)::row_ids_arg(:)
    complex(8),intent(in)::hrows_arg(:,:),srows_arg(:,:)
    real(8),intent(in)::tolerance_arg
    complex(8),allocatable,intent(out)::coefficients_arg(:,:)
    real(8),intent(out)::eigenvalues_arg(:),maximum_residual_arg,orthogonality_defect_arg,projector_defect_arg
    integer(8),intent(out)::workspace_peak_bytes_arg,fingerprint_arg
    logical,intent(out)::callback_ok
    character(*),intent(out)::callback_message
    call solve_dg_hybrid_generalized_scalapack(comm_arg,global_count_arg,nstate_arg,row_ids_arg,hrows_arg,&
      srows_arg,tolerance_arg,coefficients_arg,eigenvalues_arg,maximum_residual_arg,orthogonality_defect_arg,&
      projector_defect_arg,workspace_peak_bytes_arg,fingerprint_arg,callback_ok,callback_message)
  end subroutine solve_final_dg_hybrid_divided_lcfo

  ! Apply SALMON's established total-system Hamiltonian to one bounded tile.
  ! The callback contract supplies values in the current row-owned physical-ID
  ! order.  We explicitly map those IDs into dc%mg_tot instead of assuming that
  ! the two local array orders happen to coincide.
  subroutine apply_ow_full_cell_hpsi_tile(tile_in,tile_out,callback_ok)
    complex(8),intent(in)::tile_in(:,:)
    complex(8),intent(out)::tile_out(:,:)
    logical,intent(out)::callback_ok
    type(s_parallel_info)::tile_info
    type(s_orbital)::tile_psi,tile_hpsi
    type(s_sendrecv_grid)::tile_srg
    type(s_scalar),allocatable::zero_vlocal(:)
    complex(8),allocatable::mg_tile_in(:,:),mg_tile_out(:,:),projector_local(:,:),projector_world(:,:),&
      atom_uVpsi(:,:,:,:,:),atom_uVpsi_reduced(:,:,:,:,:)
    integer::width,p,ix,iy,iz,io,local_index,nx,ny,ierr,local_bad,global_bad,allocation_status,rank,&
      full_nonfinite_count,local_only_nonfinite_count,shape_bad,parallel_bad,owned_grid_count,ilma,ia,j
    integer(8)::physical_id,nxy,redistribution_workspace
    logical::full_output_finite,local_only_finite,redistribution_ok
    character(256)::redistribution_message
    real(8)::tile_input_peak,vlocal_peak,local_projector_difference,global_projector_difference,&
      local_atom_scale,global_atom_scale,local_world_scale,global_world_scale
    callback_ok=.false.;tile_out=(0d0,0d0);local_bad=0;full_output_finite=.true.
    call MPI_Comm_rank(dc%icomm_tot,rank,ierr)
    width=size(tile_in,1)
    shape_bad=merge(1,0,size(tile_in,2)/=size(ow_core_ids).or.any(shape(tile_out)/=shape(tile_in)))
    parallel_bad=merge(1,0,dc%info_tot%isize_o/=1.or.dc%info_tot%isize_k/=1.or.&
      dc%info_tot%numk/=1.or.dc%info_tot%numm/=1.or.dc%system_tot%nspin/=1)
    local_bad=max(shape_bad,parallel_bad)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(local_bad/=0)write(0,'(a,i0,2(a,i0),a,2(i0,1x),a,2(i0,1x),6(a,i0))')&
        '[OW-HPSI-CONTRACT-FAILURE] rank=',rank,' shape_bad=',shape_bad,' parallel_bad=',parallel_bad,&
        ' tile_shape=',shape(tile_in),' output_shape=',shape(tile_out),&
        ' core_count=',size(ow_core_ids),' isize_o=',dc%info_tot%isize_o,&
        ' isize_k=',dc%info_tot%isize_k,' numk=',dc%info_tot%numk,&
        ' numm=',dc%info_tot%numm,' nspin=',dc%system_tot%nspin
      return
    endif
    nx=dc%lg_tot%num(1);ny=dc%lg_tot%num(2);nxy=int(nx,8)*int(ny,8)
    owned_grid_count=product(dc%mg_tot%ie-dc%mg_tot%is+1)
    if(.not.allocated(ow_hpsi_grid_ids))then
      allocate(ow_hpsi_grid_ids(owned_grid_count),stat=allocation_status)
      local_bad=merge(0,1,allocation_status==0)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
      local_index=0
      do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)
      do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
      do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
        local_index=local_index+1
        ow_hpsi_grid_ids(local_index)=int(ix,8)+int(nx,8)*(int(iy-1,8)+int(ny,8)*int(iz-1,8))
      enddo
      enddo
      enddo
    endif
    if(.not.ow_hpsi_redistribution%initialized)then
      call initialize_dg_full_cell_redistribution(dc%icomm_tot,int(ow_global_grid_count),ow_core_ids,&
        ow_hpsi_grid_ids,ow_hpsi_redistribution,redistribution_workspace,redistribution_ok,&
        redistribution_message)
      if(.not.redistribution_ok)then
        if(rank==0)write(0,'(2a)')'[OW-HPSI-REDISTRIBUTION-FAILURE] ',trim(redistribution_message)
        return
      endif
      ow_hpsi_redistribution_workspace=redistribution_workspace
      if(rank==0)write(*,'(a,i0)')'[OW-GS-DIAGNOSTIC] hpsi_redistribution_workspace_bytes=',&
        ow_hpsi_redistribution_workspace
    endif
    allocate(mg_tile_in(width,owned_grid_count),mg_tile_out(width,owned_grid_count),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call apply_dg_full_cell_redistribution_forward(ow_hpsi_redistribution,tile_in,mg_tile_in,&
      redistribution_ok,redistribution_message)
    if(.not.redistribution_ok)then
      if(rank==0)write(0,'(2a)')'[OW-HPSI-REDISTRIBUTION-FAILURE] ',trim(redistribution_message)
      return
    endif
    tile_info%im_s=1;tile_info%im_e=1;tile_info%numm=1
    tile_info%ik_s=1;tile_info%ik_e=1;tile_info%numk=1
    tile_info%io_s=1;tile_info%io_e=width;tile_info%numo=width
    tile_info%if_divide_rspace=dc%info_tot%if_divide_rspace
    tile_info%if_divide_orbit=.false.
    tile_info%icomm_r=dc%info_tot%icomm_r
    tile_info%icomm_rko=dc%info_tot%icomm_rko
    allocate(tile_psi%zwf(dc%mg_tot%is_array(1):dc%mg_tot%ie_array(1),&
      dc%mg_tot%is_array(2):dc%mg_tot%ie_array(2),dc%mg_tot%is_array(3):dc%mg_tot%ie_array(3),&
      1,1:width,1:1,1:1),tile_hpsi%zwf(dc%mg_tot%is_array(1):dc%mg_tot%ie_array(1),&
      dc%mg_tot%is_array(2):dc%mg_tot%ie_array(2),dc%mg_tot%is_array(3):dc%mg_tot%ie_array(3),&
      1,1:width,1:1,1:1),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(tile_psi%zwf))deallocate(tile_psi%zwf)
      if(allocated(tile_hpsi%zwf))deallocate(tile_hpsi%zwf)
      return
    endif
    tile_psi%zwf=(0d0,0d0);tile_hpsi%zwf=(0d0,0d0)
    call init_sendrecv_grid(tile_srg,dc%mg_tot,width,dc%info_tot%icomm_rko,dc%srg_tot%neig)
    local_index=0
    do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)
    do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
    do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
      local_index=local_index+1
      do io=1,width;tile_psi%zwf(ix,iy,iz,1,io,1,1)=mg_tile_in(io,local_index);enddo
    enddo
    enddo
    enddo
    if(.not.ow_projector_stage_diagnosed)then
      allocate(projector_local(width,dc%ppg_tot%Nlma),projector_world(width,dc%ppg_tot%Nlma))
      projector_local=0d0
      do ilma=1,dc%ppg_tot%Nlma
        ia=dc%ppg_tot%ia_tbl(ilma)
        do j=1,dc%ppg_tot%mps(ia)
          ix=dc%ppg_tot%jxyz(1,j,ia);iy=dc%ppg_tot%jxyz(2,j,ia);iz=dc%ppg_tot%jxyz(3,j,ia)
          do io=1,width
            projector_local(io,ilma)=projector_local(io,ilma)+&
              conjg(dc%ppg_tot%zekr_uV(j,ilma,1))*tile_psi%zwf(ix,iy,iz,1,io,1,1)
          enddo
        enddo
        projector_local(:,ilma)=projector_local(:,ilma)*dc%ppg_tot%rinv_uvu(ilma)
      enddo
      call MPI_Allreduce(projector_local,projector_world,size(projector_local),MPI_DOUBLE_COMPLEX,&
        MPI_SUM,dc%icomm_tot,ierr)
      call calc_uVpsi_rdivided(1,tile_info,dc%ppg_tot,tile_psi,atom_uVpsi,atom_uVpsi_reduced)
      local_projector_difference=0d0;local_atom_scale=0d0;local_world_scale=0d0
      do ilma=1,dc%ppg_tot%Nlma
        ia=dc%ppg_tot%ia_tbl(ilma)
        if(.not.dc%ppg_tot%ireferred_atom(ia))cycle
        local_projector_difference=max(local_projector_difference,&
          maxval(abs(atom_uVpsi_reduced(1,1:width,1,1,ilma)-projector_world(:,ilma))))
        local_atom_scale=max(local_atom_scale,maxval(abs(atom_uVpsi_reduced(1,1:width,1,1,ilma))))
        local_world_scale=max(local_world_scale,maxval(abs(projector_world(:,ilma))) )
      enddo
      call MPI_Allreduce(local_projector_difference,global_projector_difference,1,&
        MPI_DOUBLE_PRECISION,MPI_MAX,dc%icomm_tot,ierr)
      call MPI_Allreduce(local_atom_scale,global_atom_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,dc%icomm_tot,ierr)
      call MPI_Allreduce(local_world_scale,global_world_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,dc%icomm_tot,ierr)
      if(rank==0)write(*,'(a,3(a,es16.8))')&
        '[OW-GS-DIAGNOSTIC] projector_atom_comm/world_overlap_difference',&
        ' difference=',global_projector_difference,' atom_comm_scale=',global_atom_scale,&
        ' world_scale=',global_world_scale
      deallocate(projector_local,projector_world,atom_uVpsi,atom_uVpsi_reduced)
      ow_projector_stage_diagnosed=.true.
    endif
    if(local_index==owned_grid_count)then
      select case(ow_full_cell_component_mode)
      case(1)
        allocate(zero_vlocal(1));call allocate_scalar(dc%mg_tot,zero_vlocal(1))
        zero_vlocal(1)%f=0d0
        call hpsi(tile_psi,tile_hpsi,tile_info,dc%mg_tot,zero_vlocal,dc%system_tot,stencil,&
          tile_srg,dc%ppg_tot,include_nonlocal=.false.)
      case(2)
        do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)
        do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
        do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
          do io=1,width
            tile_hpsi%zwf(ix,iy,iz,1,io,1,1)=&
              dc%vloc_tot(1)%f(ix,iy,iz)*tile_psi%zwf(ix,iy,iz,1,io,1,1)
          enddo
        enddo
        enddo
        enddo
      case default
        call hpsi(tile_psi,tile_hpsi,tile_info,dc%mg_tot,dc%vloc_tot,dc%system_tot,stencil,&
          tile_srg,dc%ppg_tot,include_nonlocal=.false.)
        call apply_ow_world_reduced_nonlocal(tile_psi,tile_hpsi,width,redistribution_ok,&
          redistribution_message)
        if(.not.redistribution_ok)then
          if(rank==0)write(0,'(2a)')'[OW-HPSI-NONLOCAL-FAILURE] ',trim(redistribution_message)
          return
        endif
        full_output_finite=all(ieee_is_finite(real(tile_hpsi%zwf))).and.&
          all(ieee_is_finite(aimag(tile_hpsi%zwf)))
        if(.not.full_output_finite)then
          full_nonfinite_count=count(.not.ieee_is_finite(real(tile_hpsi%zwf)))+&
            count(.not.ieee_is_finite(aimag(tile_hpsi%zwf)))
          tile_input_peak=maxval(abs(tile_psi%zwf))
          vlocal_peak=maxval(abs(dc%vloc_tot(1)%f))
          tile_hpsi%zwf=(0d0,0d0)
          call hpsi(tile_psi,tile_hpsi,tile_info,dc%mg_tot,dc%vloc_tot,dc%system_tot,stencil,&
            tile_srg,dc%ppg_tot,include_nonlocal=.false.)
          local_only_finite=all(ieee_is_finite(real(tile_hpsi%zwf))).and.&
            all(ieee_is_finite(aimag(tile_hpsi%zwf)))
          local_only_nonfinite_count=count(.not.ieee_is_finite(real(tile_hpsi%zwf)))+&
            count(.not.ieee_is_finite(aimag(tile_hpsi%zwf)))
          write(0,'(a,i0,2(a,l1),2(a,i0),2(a,es16.8))')'[OW-HPSI-FAILURE] rank=',rank,&
            ' full_finite=',full_output_finite,' local_only_finite=',local_only_finite,&
            ' full_nonfinite=',full_nonfinite_count,' local_only_nonfinite=',local_only_nonfinite_count,&
            ' input_peak=',tile_input_peak,' vlocal_peak=',vlocal_peak
        endif
      end select
      local_index=0
      do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)
      do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
      do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
        local_index=local_index+1
        do io=1,width;mg_tile_out(io,local_index)=tile_hpsi%zwf(ix,iy,iz,1,io,1,1);enddo
      enddo
      enddo
      enddo
      call apply_dg_full_cell_redistribution_reverse(ow_hpsi_redistribution,mg_tile_out,tile_out,&
        redistribution_ok,redistribution_message)
      callback_ok=redistribution_ok.and.full_output_finite.and.&
        all(ieee_is_finite(real(tile_out))).and.all(ieee_is_finite(aimag(tile_out)))
    endif
    call dealloc_cache(tile_srg)
    if(allocated(zero_vlocal))then
      if(allocated(zero_vlocal(1)%f))deallocate(zero_vlocal(1)%f)
      deallocate(zero_vlocal)
    endif
    if(allocated(tile_psi%zwf))deallocate(tile_psi%zwf)
    if(allocated(tile_hpsi%zwf))deallocate(tile_hpsi%zwf)
    if(allocated(mg_tile_in))deallocate(mg_tile_in)
    if(allocated(mg_tile_out))deallocate(mg_tile_out)
  end subroutine apply_ow_full_cell_hpsi_tile

  subroutine apply_ow_world_reduced_nonlocal(tile_psi,tile_hpsi,width,ok,message)
    type(s_orbital),intent(in)::tile_psi
    type(s_orbital),intent(inout)::tile_hpsi
    integer,intent(in)::width
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::local_overlap(:,:),global_overlap(:,:)
    complex(8)::coefficient
    integer::ilma,ia,j,ix,iy,iz,io,ierr,allocation_status,local_bad,global_bad
    ok=.false.;message='';local_bad=0
    allocate(local_overlap(width,dc%ppg_tot%Nlma),global_overlap(width,dc%ppg_tot%Nlma),&
      stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='cannot allocate world-reduced nonlocal overlap';return
    endif
    local_overlap=0d0
    do ilma=1,dc%ppg_tot%Nlma
      ia=dc%ppg_tot%ia_tbl(ilma)
      do j=1,dc%ppg_tot%mps(ia)
        ix=dc%ppg_tot%jxyz(1,j,ia);iy=dc%ppg_tot%jxyz(2,j,ia);iz=dc%ppg_tot%jxyz(3,j,ia)
        do io=1,width
          local_overlap(io,ilma)=local_overlap(io,ilma)+&
            conjg(dc%ppg_tot%zekr_uV(j,ilma,1))*tile_psi%zwf(ix,iy,iz,1,io,1,1)
        enddo
      enddo
    enddo
    call MPI_Allreduce(local_overlap,global_overlap,size(local_overlap),MPI_DOUBLE_COMPLEX,&
      MPI_SUM,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)then;message='world-reduced nonlocal overlap failed';return;endif
    do ilma=1,dc%ppg_tot%Nlma
      ia=dc%ppg_tot%ia_tbl(ilma)
      do j=1,dc%ppg_tot%mps(ia)
        ix=dc%ppg_tot%jxyz(1,j,ia);iy=dc%ppg_tot%jxyz(2,j,ia);iz=dc%ppg_tot%jxyz(3,j,ia)
        do io=1,width
          coefficient=dc%ppg_tot%rinv_uvu(ilma)*global_overlap(io,ilma)
          tile_hpsi%zwf(ix,iy,iz,1,io,1,1)=tile_hpsi%zwf(ix,iy,iz,1,io,1,1)+&
            coefficient*dc%ppg_tot%zekr_uV(j,ilma,1)
        enddo
      enddo
    enddo
    ok=.true.
  end subroutine apply_ow_world_reduced_nonlocal

  subroutine ow_distributed_hermiticity(comm,row_ids,rows,defect,scale,finite)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::rows(:,:)
    real(8),intent(out)::defect,scale
    logical,intent(out)::finite
    integer,parameter::row_batch_size=32
    integer::rank,nproc,ierr,r,i,j,nrows,local_finite,global_finite,local_bad,global_bad,&
      batch_first,batch_count
    integer,allocatable::counts(:),displacements(:)
    integer(8),allocatable::all_ids(:)
    complex(8),allocatable::block(:,:)
    real(8)::local_defect,local_scale
    local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    allocate(counts(nproc),displacements(nproc))
    call MPI_Allgather(size(row_ids),1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      defect=huge(1d0);scale=huge(1d0);finite=.false.;return
    endif
    displacements(1)=0
    do r=2,nproc;displacements(r)=displacements(r-1)+counts(r-1);enddo
    allocate(all_ids(sum(counts)))
    call MPI_Allgatherv(row_ids,size(row_ids),MPI_INTEGER8,all_ids,counts,displacements,&
      MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0.or.ierr/=MPI_SUCCESS)then
      defect=huge(1d0);scale=huge(1d0);finite=.false.;return
    endif
    local_defect=0d0;local_scale=0d0
    if(size(rows)>0)local_scale=maxval(abs(rows))
    local_finite=merge(1,0,all(ieee_is_finite(real(rows))).and.all(ieee_is_finite(aimag(rows))))
    do r=0,nproc-1
      nrows=counts(r+1)
      do batch_first=1,nrows,row_batch_size
        batch_count=min(row_batch_size,nrows-batch_first+1)
        allocate(block(batch_count,size(rows,2)));block=(0d0,0d0)
        if(rank==r)block=rows(batch_first:batch_first+batch_count-1,:)
        call MPI_Bcast(block,batch_count*size(rows,2),MPI_DOUBLE_COMPLEX,r,comm,ierr)
        if(ierr/=MPI_SUCCESS)local_bad=1
        do i=1,size(row_ids);do j=1,batch_count
          local_defect=max(local_defect,abs(rows(i,int(all_ids(&
            displacements(r+1)+batch_first+j-1)))-conjg(block(j,int(row_ids(i))))))
        enddo;enddo
        deallocate(block)
      enddo
    enddo
    call MPI_Allreduce(local_defect,defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(local_scale,scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(local_finite,global_finite,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    finite=global_finite==1.and.local_bad==0.and.ierr==MPI_SUCCESS
  end subroutine

  subroutine ow_fingerprint_distributed_matrix(comm,row_ids,matrix,fingerprint,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::matrix(:,:)
    integer(8),intent(out)::fingerprint
    logical,intent(out)::ok
    integer::i,j,ierr,local_bad,global_bad
    integer(8)::local_hash,real_bits,imaginary_bits,entry_hash
    local_bad=0
    if(size(matrix,1)/=size(row_ids))local_bad=1
    if(.not.all(ieee_is_finite(real(matrix))).or..not.all(ieee_is_finite(aimag(matrix))))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;fingerprint=0_8;ok=.false.;return;endif
    local_hash=0_8
    do j=1,size(matrix,2);do i=1,size(row_ids)
      real_bits=transfer(real(matrix(i,j),8),real_bits)
      imaginary_bits=transfer(aimag(matrix(i,j)),imaginary_bits)
      entry_hash=ieor(ishftc(real_bits,modulo(int(row_ids(i)),63)),ishftc(imaginary_bits,modulo(j+11,63)))
      entry_hash=ieor(entry_hash,ishftc(ieor(row_ids(i),ishft(int(j,8),21)),17))
      local_hash=ieor(local_hash,entry_hash)
    enddo;enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;fingerprint=0_8;ok=.false.;return;endif
    fingerprint=ieor(fingerprint,ishftc(int(size(matrix,2),8),29))
    fingerprint=ieor(fingerprint,int(z'6A09E667F3BCC909',8))
    if(fingerprint==0_8)fingerprint=1_8
    ok=.true.
  end subroutine ow_fingerprint_distributed_matrix

  integer(8) function ow_collective_operator_fingerprint(comm)
    integer,intent(in)::comm
    integer(8)::local_fingerprint,ranked_fingerprint
    integer::rank,ierr
    call MPI_Comm_rank(comm,rank,ierr)
    local_fingerprint=dg_dc_operator_fingerprint(.false.)
    ranked_fingerprint=ishftc(local_fingerprint,modulo(rank,63))
    call MPI_Allreduce(ranked_fingerprint,ow_collective_operator_fingerprint,1,MPI_INTEGER8,&
      MPI_BXOR,comm,ierr)
    ow_collective_operator_fingerprint=ieor(ow_collective_operator_fingerprint,&
      int(z'BB67AE8584CAA73B',8))
    if(ow_collective_operator_fingerprint==0_8)ow_collective_operator_fingerprint=1_8
  end function

  subroutine invert_ow_metric(metric,inverse,ok,message)
    complex(8),intent(in)::metric(:,:)
    complex(8),allocatable,intent(out)::inverse(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::n,info,i,j
    external::zpotrf,zpotri
    n=size(metric,1);ok=.false.;message=''
    if(n<1.or.size(metric,2)/=n)then;message='invalid metric extent';return;endif
    allocate(inverse,source=metric)
    call zpotrf('U',n,inverse,n,info)
    if(info/=0)then;message='metric is not positive definite';return;endif
    call zpotri('U',n,inverse,n,info)
    if(info/=0)then;message='metric inversion failed';return;endif
    do j=1,n;do i=j+1,n;inverse(i,j)=conjg(inverse(j,i));enddo;enddo
    ok=.true.
  end subroutine

  subroutine compute_ow_periodic_spread(comm,spread_max)
    integer,intent(in)::comm
    real(8),intent(out)::spread_max
    complex(8),allocatable::local_moment(:,:),global_moment(:,:)
    real(8),allocatable::local_norm(:),global_norm(:)
    complex(8)::phase
    integer::axis,wannier,point,ierr,coordinate_index
    allocate(local_moment(3,size(ow_core_values,1)),global_moment(3,size(ow_core_values,1)))
    allocate(local_norm(size(ow_core_values,1)),global_norm(size(ow_core_values,1)))
    local_moment=(0d0,0d0);local_norm=0d0
    do point=1,size(ow_core_ids)
      do wannier=1,size(ow_core_values,1)
        local_norm(wannier)=local_norm(wannier)+ow_core_weights(point)*&
          abs(ow_core_values(wannier,point))**2
      enddo
      do axis=1,3
        select case(axis)
        case(1)
          coordinate_index=int(modulo(ow_core_ids(point)-1_8,int(dc%lg_tot%num(1),8)))
        case(2)
          coordinate_index=int(modulo((ow_core_ids(point)-1_8)/int(dc%lg_tot%num(1),8),&
            int(dc%lg_tot%num(2),8)))
        case default
          coordinate_index=int((ow_core_ids(point)-1_8)/&
            (int(dc%lg_tot%num(1),8)*int(dc%lg_tot%num(2),8)))
        end select
        phase=exp(cmplx(0d0,2d0*pi*real(coordinate_index,8)/real(dc%lg_tot%num(axis),8),8))
        do wannier=1,size(ow_core_values,1)
          local_moment(axis,wannier)=local_moment(axis,wannier)+ow_core_weights(point)*&
            abs(ow_core_values(wannier,point))**2*phase
        enddo
      enddo
    enddo
    call MPI_Allreduce(local_moment,global_moment,size(local_moment),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_norm,global_norm,size(local_norm),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    spread_max=0d0
    do wannier=1,size(global_norm)
      if(global_norm(wannier)>0d0)spread_max=max(spread_max,&
        sum(max(0d0,1d0-abs(global_moment(:,wannier)/global_norm(wannier))**2)))
    enddo
  end subroutine

  subroutine write_ow_ground_state_evidence(metric_spectrum,occupied_required,nproc,rank,spread_max,&
      complete_shell_channels,core_atoms)
    real(8),intent(in)::metric_spectrum(:),spread_max
    integer,intent(in)::occupied_required,nproc,rank,complete_shell_channels,core_atoms
    real(8)::s_defect,s_scale,total_energy
    logical::finite_s
    integer::local_target,local_owned_rows,max_owned_rows,ierr
    integer(8)::local_overlap_bytes,max_overlap_bytes,max_hamiltonian_bytes
    call ow_distributed_hermiticity(dc%icomm_tot,ow_row_ids,ow_srows,s_defect,s_scale,finite_s)
    if(.not.finite_s.or.s_defect>dg_dc_gs_hermiticity_tolerance*max(1d0,s_scale))&
      error stop 'row-owned overlapping-Wannier metric evidence gate failed'
    local_owned_rows=size(ow_row_ids)
    local_overlap_bytes=int(size(ow_srows),8)*16_8
    call MPI_Allreduce(local_owned_rows,max_owned_rows,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)error stop 'row-owned matrix row-count evidence reduction failed'
    call MPI_Allreduce(local_overlap_bytes,max_overlap_bytes,1,MPI_INTEGER8,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)error stop 'row-owned overlap byte evidence reduction failed'
    call MPI_Allreduce(ow_diag_h_local_bytes,max_hamiltonian_bytes,1,MPI_INTEGER8,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)error stop 'row-owned Hamiltonian byte evidence reduction failed'
    if(rank/=0)return
    total_energy=sum(ow_checkpoint%occupations*ow_state%eigenvalues)
    local_target=size(ow_basis%center_owner_rank)/nproc
    write(*,'(a,i0)')'[OW-GS-EVIDENCE] mpi_ranks=',nproc
    write(*,'(a,i0)')'[OW-GS-EVIDENCE] omp_threads=',omp_get_max_threads()
    write(*,'(a,i0)')'[OW-GS-EVIDENCE] candidate_per_fragment=',&
      merge(dg_ow_candidate_states_per_fragment,system%no,dg_ow_candidate_states_per_fragment>0)
    write(*,'(a,i0)')'[OW-GS-EVIDENCE] target_per_fragment=',local_target
    write(*,'(a,i0)')'[OW-GS-EVIDENCE] core_atoms_per_fragment=',core_atoms
    write(*,'(a,i0)')'[OW-GS-EVIDENCE] global_target=',size(ow_basis%center_owner_rank)
    write(*,'(a,i0)')'[OW-GS-EVIDENCE] checkpoint_format_version=',3
    write(*,'(a,i0)')'[OW-GS-EVIDENCE] matrix_owned_rows_max=',max_owned_rows
    write(*,'(a,i0)')'[OW-GS-EVIDENCE] overlap_local_bytes_max=',max_overlap_bytes
    write(*,'(a,i0)')'[OW-GS-EVIDENCE] hamiltonian_local_bytes_max=',max_hamiltonian_bytes
    write(*,'(a,i0)')'[OW-GS-EVIDENCE] complete_shell_channels=',complete_shell_channels
    write(*,'(a,i0)')'[OW-GS-EVIDENCE] complete_shell_residual_rank=',&
      local_target-occupied_required
    write(*,'(a,i0)')'[OW-GS-EVIDENCE] bond_center_orbit_closed=',1
    write(*,'(a,i0)')'[OW-GS-EVIDENCE] occupied_included=',occupied_required
    write(*,'(a,i0)')'[OW-GS-EVIDENCE] occupied_required=',occupied_required
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] metric_cholesky_pivot_min=',minval(metric_spectrum)
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] metric_cholesky_pivot_max=',maxval(metric_spectrum)
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] metric_pivot_condition=',ow_checkpoint%metric_condition
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] occupied_inclusion_residual=',&
      ow_basis%occupied_inclusion_residual
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] complete_shell_inclusion_residual=',&
      ow_basis%projection_inclusion_residual
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] center_defect=',ow_basis%symmetry_closure_residual
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] spread_max=',spread_max
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] tail_value_norm=',ow_basis%boundary_value_max
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] tail_gradient_norm=',ow_basis%boundary_gradient_max
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] s_hermiticity=',s_defect
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] h_hermiticity=',ow_diag_h_hermiticity
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] t_hermiticity=',ow_diag_t_hermiticity
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] vnl_hermiticity=',ow_diag_vnl_hermiticity
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] vlocal_hermiticity=',ow_diag_vlocal_hermiticity
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] density_residual=',ow_result%density_residual
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] unmixed_density_residual=',ow_result%unmixed_density_residual
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] coefficient_residual=',ow_result%coefficient_residual
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] s_orthogonality=',ow_result%orthogonality_defect
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] trace_charge=',ow_result%trace_charge
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] integrated_charge=',ow_result%integrated_charge
    write(*,'(a,es24.16)')'[OW-GS-EVIDENCE] total_energy=',total_energy
  end subroutine

  subroutine assemble_ow_nonlocal_rows(comm,matrix_rows,ownership_count,ok,message)
    integer,intent(in)::comm
    complex(8),allocatable,intent(out)::matrix_rows(:,:)
    integer,intent(out)::ownership_count
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::local_overlap(:,:),owned_overlap(:,:)
    integer,allocatable::local_atom_ids(:),local_ordinals(:)
    integer(8),allocatable::projector_ids(:)
    real(8),allocatable::local_matrix_strength(:),local_action_strength(:),matrix_strength(:)
    logical,allocatable::complete(:,:)
    logical::global_ok
    integer::ilma,ia,j,ix,iy,iz,p,nwann,total_projectors,&
      canonical_index(3),ordinal

    nwann=size(ow_core_values,1)
    ok=.true.
    allocate(local_overlap(nwann,ppg%Nlma),local_atom_ids(ppg%Nlma),local_ordinals(ppg%Nlma),&
      local_matrix_strength(ppg%Nlma),local_action_strength(ppg%Nlma))
    local_overlap=(0d0,0d0);local_atom_ids=0;local_ordinals=0
    do ilma=1,ppg%Nlma
      ia=ppg%ia_tbl(ilma)
      call map_dc_atom_to_physical_atom(ia,local_atom_ids(ilma),ok)
      ordinal=count(ppg%ia_tbl(1:ilma)==ia);local_ordinals(ilma)=ordinal
      local_matrix_strength(ilma)=system%hvol*ppg%rinv_uvu(ilma)
      local_action_strength(ilma)=ppg%rinv_uvu(ilma)
      do j=1,ppg%mps(ia)
        ix=ppg%jxyz(1,j,ia);iy=ppg%jxyz(2,j,ia);iz=ppg%jxyz(3,j,ia)
        canonical_index=[dc_to_canonical_index(ix,ow_core_size(1),ow_buffer(1)),&
          dc_to_canonical_index(iy,ow_core_size(2),ow_buffer(2)),&
          dc_to_canonical_index(iz,ow_core_size(3),ow_buffer(3))]
        if(any(canonical_index<1).or.any(canonical_index>ow_box_size))then
          ok=.false.;cycle
        endif
        p=canonical_index(1)+ow_box_size(1)*((canonical_index(2)-1)+&
          ow_box_size(2)*(canonical_index(3)-1))
        local_overlap(:,ilma)=local_overlap(:,ilma)+ppg%uV(j,ilma)*&
          sqrt(ow_partition_weight(p))*ow_box_values(:,p)
      enddo
    enddo
    call comm_logical_and(ok,global_ok,comm);ok=global_ok
    if(.not.ok)then;message='cannot canonicalize complete physical atom/projector support';return;endif
    call collect_dg_overlapping_wannier_projector_overlaps(comm,nwann,local_atom_ids,local_ordinals,&
      local_matrix_strength,local_action_strength,local_overlap,projector_ids,matrix_strength,&
      owned_overlap,total_projectors,ok,message)
    if(.not.ok)return
    allocate(complete(nwann,size(projector_ids)))
    complete=.true.
    call assemble_dg_overlapping_wannier_nonlocal_rows(comm,nwann,ow_row_ids,projector_ids,matrix_strength,&
      owned_overlap,complete,int(total_projectors,8),matrix_rows,ownership_count,ok,message)
  end subroutine

  subroutine assemble_dg_hybrid_divided_nonlocal_rows(fragment_basis,row_ids,global_count,&
      matrix_rows,nonlocal_action,ownership_count,ok,message)
    type(s_dg_hybrid_fragment_basis),intent(in)::fragment_basis
    integer(8),intent(in)::row_ids(:)
    integer,intent(in)::global_count
    complex(8),intent(out)::matrix_rows(:,:)
    complex(8),allocatable,intent(out)::nonlocal_action(:,:)
    integer,intent(out)::ownership_count
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::local_overlap(:,:),owned_overlap(:,:),complete_overlap(:,:),assembled_matrix_rows(:,:)
    integer,allocatable::local_atom_ids(:),local_ordinals(:),complete_atom_ids(:),complete_ordinals(:)
    integer,allocatable::support_core_positions(:),support_projector_positions(:)
    integer(8),allocatable::projector_ids(:)
    real(8),allocatable::local_matrix_strength(:),local_action_strength(:),owned_matrix_strength(:),&
      complete_action_strength(:)
    logical,allocatable::complete(:,:)
    complex(8),allocatable::support_projector_values(:)
    integer::ilma,ia,j,ix,iy,iz,ix_tot,iy_tot,iz_tot,basis,position,ordinal,total_projectors,q,core_position,&
      local_bad,global_bad,ierr,support_count
    integer(8)::point_id

    ok=.false.;message='';ownership_count=0;local_bad=0
    if(global_count<1.or..not.allocated(fragment_basis%global_ids).or.&
        .not.allocated(fragment_basis%buffer_point_ids).or..not.allocated(fragment_basis%buffer_values))then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(global_bad/=0)then;message='invalid divided nonlocal fragment basis';return;endif
    allocate(local_overlap(global_count,ppg%Nlma),local_atom_ids(ppg%Nlma),&
      local_ordinals(ppg%Nlma),local_matrix_strength(ppg%Nlma),local_action_strength(ppg%Nlma))
    local_overlap=(0d0,0d0);local_atom_ids=0;local_ordinals=0
    local_matrix_strength=0d0;local_action_strength=0d0
    do ilma=1,ppg%Nlma
      ia=ppg%ia_tbl(ilma)
      call map_dc_atom_to_physical_atom(ia,local_atom_ids(ilma),ok)
      if(.not.ok)then;local_bad=1;cycle;endif
      ordinal=count(ppg%ia_tbl(1:ilma)==ia);local_ordinals(ilma)=ordinal
      local_matrix_strength(ilma)=system%hvol*ppg%rinv_uvu(ilma)
      local_action_strength(ilma)=ppg%rinv_uvu(ilma)
      do j=1,ppg%mps(ia)
        ix=ppg%jxyz(1,j,ia);iy=ppg%jxyz(2,j,ia);iz=ppg%jxyz(3,j,ia)
        ix_tot=dc%jxyz_tot(ix,1);iy_tot=dc%jxyz_tot(iy,2);iz_tot=dc%jxyz_tot(iz,3)
        point_id=1_8+int(ix_tot-1,8)+int(dc%lg_tot%num(1),8)*(&
          int(iy_tot-1,8)+int(dc%lg_tot%num(2),8)*int(iz_tot-1,8))
        position=findloc(fragment_basis%buffer_point_ids,point_id,dim=1)
        if(position<=0)then;local_bad=1;cycle;endif
        do basis=1,size(fragment_basis%global_ids)
          if(fragment_basis%global_ids(basis)<1_8.or.fragment_basis%global_ids(basis)>int(global_count,8))then
            local_bad=1;cycle
          endif
          local_overlap(int(fragment_basis%global_ids(basis)),ilma)=&
            local_overlap(int(fragment_basis%global_ids(basis)),ilma)+&
            ppg%uV(j,ilma)*fragment_basis%buffer_values(position,basis)
        enddo
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(global_bad/=0)then;message='divided basis omits or misidentifies nonlocal projector support';return;endif
    call collect_dg_overlapping_wannier_projector_overlaps(dc%icomm_tot,global_count,local_atom_ids,&
      local_ordinals,local_matrix_strength,local_action_strength,local_overlap,projector_ids,&
      owned_matrix_strength,owned_overlap,total_projectors,ok,message,complete_atom_ids,complete_ordinals,&
      complete_action_strength=complete_action_strength,complete_overlap=complete_overlap)
    if(.not.ok)return
    allocate(support_core_positions(sum(ppg%mps(ppg%ia_tbl))),&
      support_projector_positions(sum(ppg%mps(ppg%ia_tbl))),&
      support_projector_values(sum(ppg%mps(ppg%ia_tbl))))
    support_count=0
    do ilma=1,ppg%Nlma
      q=0
      do position=1,size(complete_atom_ids)
        if(complete_atom_ids(position)==local_atom_ids(ilma).and.&
            complete_ordinals(position)==local_ordinals(ilma))then;q=position;exit;endif
      enddo
      if(q==0)then;local_bad=1;cycle;endif
      ia=ppg%ia_tbl(ilma)
      do j=1,ppg%mps(ia)
        ix=ppg%jxyz(1,j,ia);iy=ppg%jxyz(2,j,ia);iz=ppg%jxyz(3,j,ia)
        ix_tot=dc%jxyz_tot(ix,1);iy_tot=dc%jxyz_tot(iy,2);iz_tot=dc%jxyz_tot(iz,3)
        point_id=1_8+int(ix_tot-1,8)+int(dc%lg_tot%num(1),8)*(&
          int(iy_tot-1,8)+int(dc%lg_tot%num(2),8)*int(iz_tot-1,8))
        core_position=findloc(ow_core_ids,point_id,dim=1)
        if(core_position<=0)cycle
        support_count=support_count+1
        support_core_positions(support_count)=core_position
        support_projector_positions(support_count)=q
        support_projector_values(support_count)=ppg%uV(j,ilma)
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(global_bad/=0)then;message='complete nonlocal projector payload lacks local identity';return;endif
    call apply_dg_overlapping_wannier_nonlocal_action(dc%icomm_tot,global_count,size(ow_core_ids),&
      support_core_positions(1:support_count),support_projector_positions(1:support_count),&
      support_projector_values(1:support_count),complete_action_strength,complete_overlap,&
      nonlocal_action,ok,message)
    if(.not.ok)return
    allocate(complete(global_count,size(projector_ids)));complete=.true.
    call assemble_dg_overlapping_wannier_nonlocal_rows(dc%icomm_tot,global_count,row_ids,projector_ids,&
      owned_matrix_strength,owned_overlap,complete,int(total_projectors,8),assembled_matrix_rows,&
      ownership_count,ok,message)
    if(.not.ok)return
    local_bad=0
    if(.not.allocated(assembled_matrix_rows))then
      local_bad=1
    else if(any(shape(assembled_matrix_rows)/=[size(row_ids),global_count]))then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      ok=.false.;message='invalid assembled divided nonlocal row extent';return
    endif
    matrix_rows=assembled_matrix_rows
  end subroutine assemble_dg_hybrid_divided_nonlocal_rows

  subroutine map_dc_atom_to_physical_atom(local_atom,physical_atom,ok)
    integer,intent(in)::local_atom
    integer,intent(out)::physical_atom
    logical,intent(out)::ok
    real(8)::position(3),delta(3),best,tolerance
    integer::atom,sx,sy,sz,fx,fy,fz
    best=huge(best);physical_atom=0
    do atom=1,dc%system_tot%nion
      if(system%kion(local_atom)/=dc%system_tot%kion(atom))cycle
      do fz=-1,1;do fy=-1,1;do fx=-1,1
      position=system%Rion(:,local_atom)+dc%rxyz_frag(:,dc%i_frag)+&
        fx*system%primitive_a(:,1)+fy*system%primitive_a(:,2)+fz*system%primitive_a(:,3)
      do sz=-1,1;do sy=-1,1;do sx=-1,1
        delta=position-dc%system_tot%Rion(:,atom)-sx*dc%system_tot%primitive_a(:,1)-&
          sy*dc%system_tot%primitive_a(:,2)-sz*dc%system_tot%primitive_a(:,3)
        if(sum(delta*delta)<best)then;best=sum(delta*delta);physical_atom=atom;endif
      enddo;enddo;enddo
      enddo;enddo;enddo
    enddo
    tolerance=1024d0*epsilon(1d0)*max(1d0,maxval(abs(dc%system_tot%primitive_a)))**2
    ok=physical_atom>0.and.best<=tolerance
  end subroutine

  integer function dc_to_canonical_index(index,core_count,buffer_count)
    integer,intent(in)::index,core_count,buffer_count
    integer::wrapped_index
    wrapped_index=modulo(index-1,core_count+2*buffer_count)+1
    if(wrapped_index<=core_count+buffer_count)then
      dc_to_canonical_index=buffer_count+wrapped_index
    else
      dc_to_canonical_index=wrapped_index-(core_count+buffer_count)
    endif
  end function

  subroutine ow_hybrid_update_potential(input_density,callback_ok)
    real(8),intent(in)::input_density(:);logical,intent(out)::callback_ok
    callback_ok=allocated(ow_hybrid_density)
    if(callback_ok)callback_ok=size(input_density)==size(ow_hybrid_density)
    if(callback_ok)ow_hybrid_density=input_density
  end subroutine ow_hybrid_update_potential

  subroutine ow_hybrid_assemble_hamiltonian(iteration,callback_ok)
    integer,intent(in)::iteration;logical,intent(out)::callback_ok
    character(256)::callback_message
    call ow_build_hamiltonian(dc%icomm_tot,ow_hybrid_density,ow_hybrid_hrows,&
      ow_hybrid_potential,ow_hybrid_operator_fingerprint,callback_ok,callback_message,&
      update_auxiliary_pencil=.false.)
    if(.not.callback_ok.and.nproc_id_global==0)write(0,'(a)')trim(callback_message)
  end subroutine ow_hybrid_assemble_hamiltonian

  subroutine ow_hybrid_solve_occupied(iteration,total_energy,residual,electron_defect,symmetry_defect,callback_ok)
    integer,intent(in)::iteration
    real(8),intent(out)::total_energy,residual,electron_defect,symmetry_defect
    logical,intent(out)::callback_ok
    real(8)::projector_defect
    integer(8)::workspace_peak,fingerprint
    character(256)::callback_message
    call solve_dg_hybrid_generalized_scalapack(dc%icomm_tot,size(ow_hybrid_hrows,2),&
      size(ow_hybrid_occupations),ow_row_ids,ow_hybrid_hrows,ow_srows,dg_dc_gs_final_orbital_tolerance,&
      ow_hybrid_coefficients,ow_hybrid_eigenvalues,residual,ow_hybrid_orthogonality,projector_defect,&
      workspace_peak,fingerprint,callback_ok,callback_message)
    if(callback_ok)then
      total_energy=sum(ow_hybrid_occupations*ow_hybrid_eigenvalues) ! occupied band-energy indicator, not total DFT energy
      electron_defect=abs(sum(ow_hybrid_occupations)-dc%elec_num_tot)
      symmetry_defect=ow_hybrid_symmetry_defect;ow_hybrid_eigensystem_residual=residual
    else
      total_energy=huge(1d0);electron_defect=huge(1d0);symmetry_defect=huge(1d0)
    endif
    if(.not.callback_ok.and.nproc_id_global==0)write(0,'(a)')trim(callback_message)
  end subroutine ow_hybrid_solve_occupied

  subroutine ow_hybrid_reconstruct_density(output_density,callback_ok)
    real(8),intent(out)::output_density(:);logical,intent(out)::callback_ok
    real(8),allocatable::density(:)
    real(8)::electron_count
    integer(8)::workspace_peak,fingerprint
    character(256)::callback_message
    call reconstruct_dg_hybrid_density(dc%icomm_tot,int(ow_global_grid_count),ow_core_ids,ow_core_weights,&
      size(ow_core_values,1),ow_row_ids,ow_hybrid_coefficients,ow_hybrid_occupations,&
      ow_hybrid_basis_provider,min(16,size(ow_core_values,1)),min(16,size(ow_hybrid_occupations)),&
      ow_symmetry_fingerprint,dg_dc_gs_final_density_tolerance,density,electron_count,workspace_peak,&
      fingerprint,callback_ok,callback_message)
    if(callback_ok)then
      callback_ok=allocated(density)
      if(callback_ok)callback_ok=size(output_density)==size(density)
      if(callback_ok)output_density=density
      if(callback_ok)then
        callback_ok=abs(electron_count-dc%elec_num_tot)<=dg_dc_gs_electron_count_tolerance
        if(.not.callback_ok)callback_message='reconstructed hybrid density violates the electron-count gate'
      endif
    endif
    if(.not.callback_ok.and.nproc_id_global==0)write(0,'(a)')trim(callback_message)
  end subroutine ow_hybrid_reconstruct_density

  subroutine ow_hybrid_basis_provider(first_column,column_count,tile_values,callback_ok)
    integer,intent(in)::first_column,column_count
    complex(8),intent(out)::tile_values(:,:);logical,intent(out)::callback_ok
    callback_ok=first_column>=1.and.column_count>=1
    if(callback_ok)callback_ok=first_column<=size(ow_core_values,1)
    if(callback_ok)callback_ok=column_count<=size(ow_core_values,1)-first_column+1
    if(callback_ok)callback_ok=all(shape(tile_values)==[column_count,size(ow_core_values,2)])
    if(callback_ok)tile_values=ow_core_values(first_column:first_column+column_count-1,:)
  end subroutine ow_hybrid_basis_provider

  subroutine ow_hybrid_density_mix(iteration,input_density,output_density,reset_history,reduce_rate,&
      mixed_density,callback_ok)
    integer,intent(in)::iteration
    real(8),intent(in)::input_density(:),output_density(:)
    logical,intent(in)::reset_history,reduce_rate
    real(8),intent(out)::mixed_density(:);logical,intent(out)::callback_ok
    character(256)::callback_message
    integer::new_history_count
    if(reduce_rate)ow_hybrid_mixing_rate=max(1d-3,0.5d0*ow_hybrid_mixing_rate)
    if(reset_history)then
      ow_hybrid_density_history(:,1)=input_density
      ow_hybrid_density_history(:,2)=input_density
      ow_hybrid_history_count=0
    endif
    call mix_dg_overlapping_wannier_density_history(dc%icomm_tot,ow_hybrid_mixing_rate,input_density,&
      output_density,ow_hybrid_density_history,ow_hybrid_history_count,mixed_density,&
      ow_hybrid_new_history,new_history_count,callback_ok,callback_message)
    if(callback_ok)then
      ow_hybrid_density_history=ow_hybrid_new_history;ow_hybrid_history_count=new_history_count
    endif
    if(.not.callback_ok.and.nproc_id_global==0)write(0,'(a)')trim(callback_message)
  end subroutine ow_hybrid_density_mix

  subroutine ow_hybrid_density_to_dc(values,callback_ok)
    real(8),intent(in)::values(:);logical,intent(out)::callback_ok
    integer::p,ix,iy,iz
    logical::global_ok
    callback_ok=size(values)==size(ow_core_ids)
    if(.not.callback_ok)return
    do p=1,size(ow_core_ids)
      ix=int(modulo(ow_core_ids(p)-1_8,int(dc%lg_tot%num(1),8)))+1
      iy=int(modulo((ow_core_ids(p)-1_8)/int(dc%lg_tot%num(1),8),int(dc%lg_tot%num(2),8)))+1
      iz=int((ow_core_ids(p)-1_8)/(int(dc%lg_tot%num(1),8)*int(dc%lg_tot%num(2),8)))+1
      if(ix<dc%mg_tot%is(1).or.ix>dc%mg_tot%ie(1).or.iy<dc%mg_tot%is(2).or.iy>dc%mg_tot%ie(2).or.&
        iz<dc%mg_tot%is(3).or.iz>dc%mg_tot%ie(3))then;callback_ok=.false.;cycle;endif
      dc%rho_tot_s(1)%f(ix,iy,iz)=values(p)
    enddo
    call comm_logical_and(callback_ok,global_ok,dc%icomm_tot);callback_ok=global_ok
  end subroutine ow_hybrid_density_to_dc

  subroutine ow_hybrid_density_from_dc(values,callback_ok)
    real(8),intent(out)::values(:);logical,intent(out)::callback_ok
    integer::p,ix,iy,iz
    callback_ok=size(values)==size(ow_core_ids)
    if(.not.callback_ok)return
    do p=1,size(ow_core_ids)
      ix=int(modulo(ow_core_ids(p)-1_8,int(dc%lg_tot%num(1),8)))+1
      iy=int(modulo((ow_core_ids(p)-1_8)/int(dc%lg_tot%num(1),8),int(dc%lg_tot%num(2),8)))+1
      iz=int((ow_core_ids(p)-1_8)/(int(dc%lg_tot%num(1),8)*int(dc%lg_tot%num(2),8)))+1
      values(p)=dc%rho_tot_s(1)%f(ix,iy,iz)
    enddo
  end subroutine ow_hybrid_density_from_dc

  subroutine ow_mix_density(comm,mixing_rate,current_density,raw_density,history,history_count,&
      mixed_density,new_history,new_history_count,ok,message)
    integer,intent(in)::comm,history_count
    real(8),intent(in)::mixing_rate,current_density(:),raw_density(:),history(:,:)
    real(8),intent(out)::mixed_density(:),new_history(:,:)
    integer,intent(out)::new_history_count
    logical,intent(out)::ok
    character(*),intent(out)::message
    call mix_dg_overlapping_wannier_density_history(comm,mixing_rate,current_density,raw_density,&
      history,history_count,mixed_density,new_history,new_history_count,ok,message)
  end subroutine

  subroutine ow_transaction(action,ok,message)
    integer,intent(in)::action
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::is,ix,iy,iz,jx,jy,jz
    ok=.true.;message=''
    select case(action)
    case(0)
      if(allocated(ow_density_snapshot))deallocate(ow_density_snapshot)
      allocate(ow_density_snapshot(size(dc%rho_tot_s(1)%f,1),size(dc%rho_tot_s(1)%f,2),&
        size(dc%rho_tot_s(1)%f,3),system%nspin))
      ow_density_snapshot=0d0
      do is=1,system%nspin
      do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)
      do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
      do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
        jx=ix-lbound(dc%rho_tot_s(is)%f,1)+1
        jy=iy-lbound(dc%rho_tot_s(is)%f,2)+1
        jz=iz-lbound(dc%rho_tot_s(is)%f,3)+1
        ow_density_snapshot(jx,jy,jz,is)=dc%rho_tot_s(is)%f(ix,iy,iz)
      enddo
      enddo
      enddo
      enddo
      ow_potential_epoch_snapshot=dg_gs_potential_epoch
      ow_transaction_active=.true.
    case(-1)
      if(.not.ow_transaction_active)then;ok=.false.;message='missing overlapping-Wannier rollback snapshot';return;endif
      call dg_dc_update_potential_from_density(ow_density_snapshot,ok,message)
      if(ok)dg_gs_potential_epoch=ow_potential_epoch_snapshot
      ow_transaction_active=.false.
    case(1)
      ok=ow_transaction_active
      if(.not.ok)message='missing overlapping-Wannier prepared transaction'
    end select
  end subroutine

  subroutine ow_commit_transaction()
    ow_transaction_active=.false.
    if(allocated(ow_density_snapshot))deallocate(ow_density_snapshot)
  end subroutine

  subroutine promote_and_project_ow_matrices(metric,hamiltonian,hamiltonian_components,position,velocity,&
      local_integer_rotations,local_cartesian_rotations,local_product,local_representation,&
      pre_projection_defect,post_projection_defect,promoted_group_order,ok,message)
    complex(8),intent(inout)::metric(:,:),hamiltonian(:,:),position(:,:,:),velocity(:,:,:)
    complex(8),intent(in)::hamiltonian_components(:,:,:)
    integer,intent(in)::local_integer_rotations(:,:,:),local_product(:,:)
    real(8),intent(in)::local_cartesian_rotations(:,:,:)
    complex(8),intent(in)::local_representation(:,:,:)
    real(8),intent(out)::pre_projection_defect,post_projection_defect
    integer,intent(out)::promoted_group_order
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::root_integer_rotations(:,:,:),common_integer_rotations(:,:,:),&
      local_match(:),common_product(:,:),promoted(:),promoted_product(:,:)
    real(8),allocatable::root_cartesian_rotations(:,:,:),common_cartesian_rotations(:,:,:),&
      promoted_rotations(:,:,:),cross_block_scalar_residual(:,:),cross_block_vector_residual(:,:),&
      cross_block_scalar_field_residual(:,:),cross_block_vector_field_residual(:,:),&
      cross_block_hamiltonian_component_residual(:,:),&
      scalar_residual(:),vector_residual(:)
    complex(8),allocatable::raw_representation(:,:,:),promoted_raw(:,:,:),representation(:,:,:),&
      common_local_representation(:,:,:),scalars(:,:,:),vectors(:,:,:,:),&
      projected_scalars(:,:,:),projected_vectors(:,:,:,:),transformed_block(:,:),target_block(:,:)
    real(8),allocatable::fragment_centers(:,:)
    integer,allocatable::fragment_ids(:),fragment_permutation(:,:)
    logical,allocatable::fragment_exact(:,:)
    integer::rank,nproc,ierr,root_count,common_count,local_index,present_local,present_global,&
      iop,jop,kop,i,j,k,nlocal,nwann,source,target,mapped_source,mapped_target,block,axis,field
    real(8)::raw_unitarity_defect,unitarity_defect,closure_defect,scalar_scale,vector_scale,&
      scalar_defect,vector_defect
    logical::residual_ok

    ok=.false.;message='';promoted_group_order=0;pre_projection_defect=huge(1d0)
    post_projection_defect=huge(1d0)
    call MPI_Comm_rank(dc%icomm_tot,rank,ierr);call MPI_Comm_size(dc%icomm_tot,nproc,ierr)
    nwann=size(metric,1);nlocal=size(local_representation,1)
    if(size(metric,2)/=nwann.or.any(shape(hamiltonian)/=[nwann,nwann]).or. &
        any(shape(position)/=[3,nwann,nwann]).or.any(shape(velocity)/=[3,nwann,nwann]).or. &
        any(shape(hamiltonian_components)/=[nwann,nwann,3]).or. &
        size(local_representation,2)/=nlocal.or.size(local_representation,3)/= &
        size(local_integer_rotations,3).or.nwann/=nlocal*nproc.or. &
        size(local_integer_rotations,1)/=3.or.size(local_integer_rotations,2)/=3.or. &
        any(shape(local_cartesian_rotations)/=shape(local_integer_rotations)).or. &
        any(shape(local_product)/=[size(local_representation,3),size(local_representation,3)]))then
      message='invalid local/global crystallographic promotion dimensions';return
    end if
    root_count=0;if(rank==0)root_count=size(local_integer_rotations,3)
    call MPI_Bcast(root_count,1,MPI_INTEGER,0,dc%icomm_tot,ierr)
    allocate(root_integer_rotations(3,3,root_count),root_cartesian_rotations(3,3,root_count))
    if(rank==0)then
      root_integer_rotations=local_integer_rotations;root_cartesian_rotations=local_cartesian_rotations
    end if
    call MPI_Bcast(root_integer_rotations,9*root_count,MPI_INTEGER,0,dc%icomm_tot,ierr)
    call MPI_Bcast(root_cartesian_rotations,9*root_count,MPI_DOUBLE_PRECISION,0,dc%icomm_tot,ierr)
    allocate(local_match(root_count));local_match=0;common_count=0
    do iop=1,root_count
      do jop=1,size(local_integer_rotations,3)
        if(all(local_integer_rotations(:,:,jop)==root_integer_rotations(:,:,iop)))then
          local_match(iop)=jop;exit
        end if
      end do
      present_local=merge(1,0,local_match(iop)>0)
      call MPI_Allreduce(present_local,present_global,1,MPI_INTEGER,MPI_MIN,dc%icomm_tot,ierr)
      if(present_global==1)common_count=common_count+1
    end do
    if(common_count<1)then;message='fragment point groups have no common identity';return;end if
    allocate(common_integer_rotations(3,3,common_count),common_cartesian_rotations(3,3,common_count),&
      common_local_representation(nlocal,nlocal,common_count),common_product(common_count,common_count))
    common_count=0
    do iop=1,root_count
      present_local=merge(1,0,local_match(iop)>0)
      call MPI_Allreduce(present_local,present_global,1,MPI_INTEGER,MPI_MIN,dc%icomm_tot,ierr)
      if(present_global/=1)cycle
      common_count=common_count+1;common_integer_rotations(:,:,common_count)=root_integer_rotations(:,:,iop)
      common_cartesian_rotations(:,:,common_count)=root_cartesian_rotations(:,:,iop)
      common_local_representation(:,:,common_count)=local_representation(:,:,local_match(iop))
    end do
    allocate(fragment_ids(nproc),fragment_centers(3,nproc))
    call MPI_Allgather(dc%i_frag,1,MPI_INTEGER,fragment_ids,1,MPI_INTEGER,dc%icomm_tot,ierr)
    do source=1,nproc
      fragment_centers(:,source)=modulo((real(dc%ixyz_frag(:,fragment_ids(source))-1,8)+&
        0.5d0*real(ow_core_size,8))/real(dc%lg_tot%num,8),1d0)
    end do
    call build_dg_fragment_permuted_representation(common_local_representation,&
      real(common_integer_rotations,8),fragment_centers,dg_ow_symmetry_tolerance,raw_representation,&
      fragment_permutation,ok,message)
    if(.not.ok)return
    do iop=1,common_count;do jop=1,common_count
      kop=0
      do k=1,common_count
        if(all(matmul(common_integer_rotations(:,:,iop),common_integer_rotations(:,:,jop))== &
            common_integer_rotations(:,:,k)))then;kop=k;exit;end if
      end do
      if(kop==0)then;message='common fragment point rotations are not closed';return;end if
      common_product(iop,jop)=kop
    end do;end do
    allocate(scalars(nwann,nwann,2),vectors(nwann,nwann,3,2))
    scalars(:,:,1)=metric;scalars(:,:,2)=hamiltonian
    do i=1,3
      vectors(:,:,i,1)=position(i,:,:);vectors(:,:,i,2)=velocity(i,:,:)
    end do
    call evaluate_dg_covariance_residuals_by_operation(raw_representation,common_cartesian_rotations,&
      scalars,vectors,scalar_residual,vector_residual,residual_ok,message)
    if(.not.residual_ok)return
    allocate(fragment_exact(nproc,common_count),&
      cross_block_scalar_residual(nproc*nproc,common_count),&
      cross_block_vector_residual(nproc*nproc,common_count),&
      cross_block_scalar_field_residual(2,common_count),&
      cross_block_vector_field_residual(2,common_count),&
      cross_block_hamiltonian_component_residual(3,common_count),&
      transformed_block(nlocal,nlocal),target_block(nlocal,nlocal));fragment_exact=.true.
    cross_block_scalar_residual=0d0;cross_block_vector_residual=0d0
    cross_block_scalar_field_residual=0d0;cross_block_vector_field_residual=0d0
    cross_block_hamiltonian_component_residual=0d0
    do iop=1,common_count
      block=0
      do source=1,nproc;do target=1,nproc
        block=block+1;mapped_source=fragment_permutation(source,iop)
        mapped_target=fragment_permutation(target,iop)
        do field=1,size(scalars,3)
          scalar_scale=max(1d0,maxval(abs(scalars(:,:,field))))
          transformed_block=matmul(conjg(transpose(raw_representation(&
            (mapped_source-1)*nlocal+1:mapped_source*nlocal,(source-1)*nlocal+1:source*nlocal,iop))),&
            matmul(scalars((mapped_source-1)*nlocal+1:mapped_source*nlocal,&
              (mapped_target-1)*nlocal+1:mapped_target*nlocal,field),raw_representation(&
              (mapped_target-1)*nlocal+1:mapped_target*nlocal,&
              (target-1)*nlocal+1:target*nlocal,iop)))
          scalar_defect=maxval(abs(transformed_block-scalars((source-1)*nlocal+1:source*nlocal,&
            (target-1)*nlocal+1:target*nlocal,field)))/scalar_scale
          cross_block_scalar_residual(block,iop)=max(cross_block_scalar_residual(block,iop),scalar_defect)
          cross_block_scalar_field_residual(field,iop)=max(&
            cross_block_scalar_field_residual(field,iop),scalar_defect)
        end do
        do field=1,size(hamiltonian_components,3)
          scalar_scale=max(1d0,maxval(abs(hamiltonian_components(:,:,field))))
          transformed_block=matmul(conjg(transpose(raw_representation(&
            (mapped_source-1)*nlocal+1:mapped_source*nlocal,(source-1)*nlocal+1:source*nlocal,iop))),&
            matmul(hamiltonian_components((mapped_source-1)*nlocal+1:mapped_source*nlocal,&
              (mapped_target-1)*nlocal+1:mapped_target*nlocal,field),raw_representation(&
              (mapped_target-1)*nlocal+1:mapped_target*nlocal,&
              (target-1)*nlocal+1:target*nlocal,iop)))
          scalar_defect=maxval(abs(transformed_block-hamiltonian_components(&
            (source-1)*nlocal+1:source*nlocal,(target-1)*nlocal+1:target*nlocal,field)))/scalar_scale
          cross_block_hamiltonian_component_residual(field,iop)=max(&
            cross_block_hamiltonian_component_residual(field,iop),scalar_defect)
        end do
        do field=1,size(vectors,4)
          vector_scale=max(1d0,maxval(abs(vectors(:,:,:,field))))
          do axis=1,3
          transformed_block=matmul(conjg(transpose(raw_representation(&
            (mapped_source-1)*nlocal+1:mapped_source*nlocal,(source-1)*nlocal+1:source*nlocal,iop))),&
            matmul(vectors((mapped_source-1)*nlocal+1:mapped_source*nlocal,&
              (mapped_target-1)*nlocal+1:mapped_target*nlocal,axis,field),raw_representation(&
              (mapped_target-1)*nlocal+1:mapped_target*nlocal,&
              (target-1)*nlocal+1:target*nlocal,iop)))
          target_block=(0d0,0d0)
          do j=1,3
            target_block=target_block+common_cartesian_rotations(axis,j,iop)*vectors(&
              (source-1)*nlocal+1:source*nlocal,(target-1)*nlocal+1:target*nlocal,j,field)
          end do
          vector_defect=maxval(abs(transformed_block-target_block))/vector_scale
          cross_block_vector_residual(block,iop)=max(cross_block_vector_residual(block,iop),vector_defect)
          cross_block_vector_field_residual(field,iop)=max(&
            cross_block_vector_field_residual(field,iop),vector_defect)
          end do
        end do
      end do;end do
    end do
    call promote_dg_exact_global_subgroup(common_product,fragment_exact,cross_block_scalar_residual,&
      cross_block_vector_residual,dg_ow_symmetry_tolerance,promoted,ok,message)
    if(.not.ok)return
    if(rank==0)then
      do iop=1,common_count
        write(*,'(a,i0,6(a,es12.4),a,l1)')&
          '[OW-GS-DIAGNOSTIC] common_point_operation=',iop,&
          ' cross_block_scalar_residual=',maxval(cross_block_scalar_residual(:,iop)),&
          ' cross_block_vector_residual=',maxval(cross_block_vector_residual(:,iop)),&
          ' S_residual=',cross_block_scalar_field_residual(1,iop),&
          ' H_residual=',cross_block_scalar_field_residual(2,iop),&
          ' X_residual=',cross_block_vector_field_residual(1,iop),&
          ' V_residual=',cross_block_vector_field_residual(2,iop),&
          ' promoted=',any(promoted==iop)
        write(*,'(a,i0,3(a,es12.4))')&
          '[OW-GS-DIAGNOSTIC] common_point_operation_components=',iop,&
          ' T_residual=',cross_block_hamiltonian_component_residual(1,iop),&
          ' Vlocal_residual=',cross_block_hamiltonian_component_residual(2,iop),&
          ' Vnonlocal_residual=',cross_block_hamiltonian_component_residual(3,iop)
      end do
    end if
    promoted_group_order=size(promoted)
    allocate(promoted_raw(nwann,nwann,promoted_group_order),&
      promoted_rotations(3,3,promoted_group_order),promoted_product(promoted_group_order,promoted_group_order))
    promoted_raw=raw_representation(:,:,promoted);promoted_rotations=common_cartesian_rotations(:,:,promoted)
    do i=1,promoted_group_order;do j=1,promoted_group_order
      kop=common_product(promoted(i),promoted(j));k=findloc(promoted,kop,dim=1)
      if(k==0)then;message='promoted crystallographic subgroup is not closed';ok=.false.;return;end if
      promoted_product(i,j)=k
    end do;end do
    call build_dg_fragment_group_representation(metric,promoted_raw,promoted_product,&
      dg_ow_symmetry_tolerance,representation,raw_unitarity_defect,unitarity_defect,&
      closure_defect,ok,message)
    if(.not.ok)return
    call project_dg_fragment_covariant_operators(representation,promoted_rotations,scalars,vectors,&
      dg_ow_symmetry_tolerance,projected_scalars,projected_vectors,pre_projection_defect,&
      post_projection_defect,ok,message)
    if(.not.ok)return
    metric=projected_scalars(:,:,1);hamiltonian=projected_scalars(:,:,2)
    do i=1,3
      position(i,:,:)=projected_vectors(:,:,i,1);velocity(i,:,:)=projected_vectors(:,:,i,2)
    end do
  end subroutine promote_and_project_ow_matrices

  subroutine prepare_ow_global_point_action(local_ids,target_ids,integer_rotations,rotations,&
      fractional_translations,product_table,translation_subgroup,point_representatives,&
      point_product,translation_cocycle,inversion_present,ok,message)
    integer(8),intent(in)::local_ids(:)
    integer(8),allocatable,intent(out)::target_ids(:,:)
    integer,allocatable,intent(out)::integer_rotations(:,:,:)
    real(8),allocatable,intent(out)::rotations(:,:,:)
    real(8),allocatable,intent(out)::fractional_translations(:,:)
    integer,allocatable,intent(out)::product_table(:,:)
    integer,allocatable,intent(out)::translation_subgroup(:),point_representatives(:),&
      point_product(:,:),translation_cocycle(:,:)
    logical,intent(out)::inversion_present,ok
    character(*),intent(out)::message
    type(t_sawf_crystallographic_catalog)::catalog
    type(t_sawf_operation_index)::operation_index
    type(t_sawf_symop),allocatable::selected_catalog_operations(:)
    real(8),allocatable::fractional_positions(:,:)
    integer,allocatable::species(:),selected(:),mapped_owner(:),mapped_local(:),mapped_wrap(:,:)
    integer(8),allocatable::all_ids(:,:),mapped_ids(:)
    real(8)::lattice_inverse(3,3),determinant,common_center(3),common_center_residual
    integer::rank,nproc,ierr,nlocal,atom,operation,axis,translation_grid(3),nselected,&
      g,h
    logical::inverse_ok,map_ok,matched,duplicate_operation,have_common_center,center_solver_ok
    character(256)::detail

    ok=.false.;inversion_present=.false.;message=''
    call MPI_Comm_rank(dc%icomm_tot,rank,ierr);call MPI_Comm_size(dc%icomm_tot,nproc,ierr)
    nlocal=size(local_ids)
    if(nlocal<1)then;message='global point action has no local core IDs';return;end if
    call invert_ow_lattice(dc%system_tot%primitive_a,lattice_inverse,determinant,inverse_ok)
    if(.not.inverse_ok)then;message='global point-action lattice is singular';return;end if
    allocate(fractional_positions(3,dc%system_tot%nion),species(dc%system_tot%nion))
    do atom=1,dc%system_tot%nion
      fractional_positions(:,atom)=modulo(matmul(lattice_inverse,dc%system_tot%Rion(:,atom)),1d0)
      species(atom)=dc%system_tot%kion(atom)
    end do
    call load_sawf_crystallographic_catalog_auto(dc%system_tot%primitive_a,fractional_positions,&
      species,dg_ow_symmetry_tolerance,catalog,map_ok,detail)
    if(.not.map_ok)then;message='global point-action catalog: '//trim(detail);return;end if
    allocate(selected(size(catalog%operations)));nselected=0
    do operation=1,size(catalog%operations)
      do axis=1,3
        translation_grid(axis)=nint(catalog%fractional_translation(axis,operation)*&
          real(dc%lg_tot%num(axis),8))
      end do
      if(maxval(abs(real(translation_grid,8)/real(dc%lg_tot%num,8)-&
          catalog%fractional_translation(:,operation)-anint(real(translation_grid,8)/&
          real(dc%lg_tot%num,8)-catalog%fractional_translation(:,operation))))>&
          dg_ow_symmetry_tolerance)cycle
      duplicate_operation=.false.
      do g=1,nselected
        if(all(catalog%integer_rotation(:,:,selected(g))==&
            catalog%integer_rotation(:,:,operation)).and.maxval(abs(&
            catalog%fractional_translation(:,selected(g))-catalog%fractional_translation(:,operation)-&
            anint(catalog%fractional_translation(:,selected(g))-&
            catalog%fractional_translation(:,operation))))<=dg_ow_symmetry_tolerance)then
          duplicate_operation=.true.;exit
        end if
      end do
      if(duplicate_operation)cycle
      nselected=nselected+1;selected(nselected)=operation
      if(all(catalog%integer_rotation(:,:,operation)==&
          reshape([-1,0,0,0,-1,0,0,0,-1],[3,3])))inversion_present=.true.
    end do
    if(nselected<1)then;message='global point-action catalog has no grid-commensurate operation';return;end if
    call solve_dg_affine_common_fixed_point(catalog%integer_rotation(:,:,selected(1:nselected)),&
      catalog%fractional_translation(:,selected(1:nselected)),dg_ow_symmetry_tolerance,&
      have_common_center,common_center,common_center_residual,center_solver_ok,detail)
    if(.not.center_solver_ok)then;message='global affine center solve: '//trim(detail);return;end if
    if(rank==0)write(*,'(a,i0,2(a,l1),a,3(es12.4,1x),a,es12.4)')&
      '[OW-GS-DIAGNOSTIC] global_affine_group_order=',nselected,&
      ' inversion=',inversion_present,' common_center=',have_common_center,&
      ' center_fractional=',common_center,' center_residual=',common_center_residual
    allocate(all_ids(nlocal,nproc),target_ids(nlocal,nselected),integer_rotations(3,3,nselected),&
      rotations(3,3,nselected),fractional_translations(3,nselected))
    call MPI_Allgather(local_ids,nlocal,MPI_INTEGER8,all_ids,nlocal,MPI_INTEGER8,&
      dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)then;message='global point-action owner gather failed';return;end if
    do g=1,nselected
      operation=selected(g);integer_rotations(:,:,g)=catalog%integer_rotation(:,:,operation)
      rotations(:,:,g)=catalog%operations(operation)%R
      fractional_translations(:,g)=catalog%fractional_translation(:,operation)
      call build_dg_pointwise_affine_owner_map(dc%lg_tot%num,local_ids,all_ids,&
        catalog%integer_rotation(:,:,operation),catalog%fractional_translation(:,operation),&
        dg_ow_symmetry_tolerance,mapped_ids,mapped_owner,mapped_local,mapped_wrap,map_ok,detail)
      if(.not.map_ok)then;message='global point-action map: '//trim(detail);return;end if
      target_ids(:,g)=int(mapped_owner,8)*int(nlocal,8)+int(mapped_local,8)
      deallocate(mapped_ids,mapped_owner,mapped_local,mapped_wrap)
    end do
    allocate(selected_catalog_operations(nselected),source=catalog%operations(selected(1:nselected)))
    call build_sawf_operation_index(selected_catalog_operations,dg_ow_symmetry_tolerance,&
      operation_index,map_ok,detail)
    if(.not.map_ok)then;message='global affine operation index: '//trim(detail);return;endif
    allocate(product_table(nselected,nselected));product_table=0
    do g=1,nselected;do h=1,nselected
      call lookup_sawf_operation_product(operation_index,selected_catalog_operations,g,h,&
        product_table(g,h),map_ok,detail)
      if(.not.map_ok)then;message='global affine product lookup: '//trim(detail);return;endif
    enddo;enddo
    call factor_dg_affine_translation_cocycle(integer_rotations,fractional_translations,&
      product_table,dg_ow_symmetry_tolerance,translation_subgroup,point_representatives,&
      point_product,translation_cocycle,map_ok,detail)
    if(.not.map_ok)then;message='global affine factorization: '//trim(detail);return;end if
    if(rank==0)write(*,'(2(a,i0))')'[OW-GS-DIAGNOSTIC] translation_subgroup_order=',&
      size(translation_subgroup),' point_cogroup_order=',size(point_representatives)
    ok=.true.
  end subroutine prepare_ow_global_point_action

  subroutine prepare_ow_fixed_center_group(local_ids,operations,target_ids,product_table,&
      fixed_center,inversion_present,group_fingerprint,ok,message)
    integer(8),intent(in)::local_ids(:)
    type(t_sawf_symop),allocatable,intent(out)::operations(:)
    integer(8),allocatable,intent(out)::target_ids(:,:)
    integer,allocatable,intent(out)::product_table(:,:)
    real(8),intent(out)::fixed_center(3)
    logical,intent(out)::inversion_present,ok
    integer(8),intent(out)::group_fingerprint
    character(*),intent(out)::message
    type(t_sawf_crystallographic_catalog)::catalog
    type(t_sawf_operation_index)::operation_index
    real(8),allocatable::fractional_positions(:,:)
    integer,allocatable::species(:),selected(:),candidate_selected(:),mapped_owner(:),mapped_local(:),mapped_wrap(:,:)
    integer(8),allocatable::all_ids(:,:),mapped_ids(:)
    real(8)::lattice_inverse(3,3),determinant,fixed_residual(3),center_min(3),center_max(3),&
      candidate_center(3),best_norm,candidate_norm,best_lex,candidate_lex
    integer::rank,nproc,ierr,nlocal,atom,operation,inversion_operation,identity_operation,&
      selected_count,selected_count_min,selected_count_max,candidate_count,best_count,g,h,axis,&
      candidate_operation,sx,sy,sz,translation_grid(3)
    integer(8)::fingerprint_min,fingerprint_max
    logical::inverse_ok,map_ok,duplicate_rotation
    character(256)::detail

    ok=.false.;message='';fixed_center=0d0;inversion_present=.false.;group_fingerprint=0_8
    call MPI_Comm_rank(dc%icomm_tot,rank,ierr);call MPI_Comm_size(dc%icomm_tot,nproc,ierr)
    nlocal=size(local_ids)
    if(nlocal<1)then;message='fixed-center group has no local core IDs';return;endif
    call invert_ow_lattice(dc%system_tot%primitive_a,lattice_inverse,determinant,inverse_ok)
    if(.not.inverse_ok)then;message='fixed-center lattice is singular';return;endif
    allocate(fractional_positions(3,dc%system_tot%nion),species(dc%system_tot%nion))
    do atom=1,dc%system_tot%nion
      fractional_positions(:,atom)=modulo(matmul(lattice_inverse,dc%system_tot%Rion(:,atom)),1d0)
      species(atom)=dc%system_tot%kion(atom)
    enddo
    call load_sawf_crystallographic_catalog_auto(dc%system_tot%primitive_a,fractional_positions,&
      species,dg_ow_symmetry_tolerance,catalog,map_ok,detail)
    if(.not.map_ok)then;message='fixed-center catalog: '//trim(detail);return;endif
    inversion_operation=0;identity_operation=0
    do operation=1,size(catalog%operations)
      do axis=1,3
        translation_grid(axis)=nint(catalog%fractional_translation(axis,operation)*&
          real(dc%lg_tot%num(axis),8))
      enddo
      if(maxval(abs(real(translation_grid,8)/real(dc%lg_tot%num,8)-&
          catalog%fractional_translation(:,operation)-anint(real(translation_grid,8)/&
          real(dc%lg_tot%num,8)-catalog%fractional_translation(:,operation))))>&
          dg_ow_symmetry_tolerance)cycle
      if(all(catalog%integer_rotation(:,:,operation)==reshape([1,0,0,0,1,0,0,0,1],[3,3])).and.&
          maxval(abs(catalog%fractional_translation(:,operation)-&
          anint(catalog%fractional_translation(:,operation))))<=dg_ow_symmetry_tolerance)&
        identity_operation=operation
    enddo
    if(identity_operation==0)then
      message='fixed-center group requires grid-commensurate identity and inversion';return
    endif
    allocate(candidate_selected(size(catalog%operations)))
    best_count=0;best_norm=huge(1d0);best_lex=huge(1d0)
    do candidate_operation=1,size(catalog%operations)
      if(any(catalog%integer_rotation(:,:,candidate_operation)/=&
          reshape([-1,0,0,0,-1,0,0,0,-1],[3,3])))cycle
      do axis=1,3
        translation_grid(axis)=nint(catalog%fractional_translation(axis,candidate_operation)*&
          real(dc%lg_tot%num(axis),8))
      enddo
      if(maxval(abs(real(translation_grid,8)/real(dc%lg_tot%num,8)-&
          catalog%fractional_translation(:,candidate_operation)-anint(real(translation_grid,8)/&
          real(dc%lg_tot%num,8)-catalog%fractional_translation(:,candidate_operation))))>&
          dg_ow_symmetry_tolerance)cycle
      do sx=0,1;do sy=0,1;do sz=0,1
      candidate_center=modulo(0.5d0*catalog%fractional_translation(:,candidate_operation)+&
        0.5d0*[real(sx,8),real(sy,8),real(sz,8)],1d0)
      candidate_count=0
      do operation=1,size(catalog%operations)
        do axis=1,3
          translation_grid(axis)=nint(catalog%fractional_translation(axis,operation)*&
            real(dc%lg_tot%num(axis),8))
        enddo
        if(maxval(abs(real(translation_grid,8)/real(dc%lg_tot%num,8)-&
            catalog%fractional_translation(:,operation)-anint(real(translation_grid,8)/&
            real(dc%lg_tot%num,8)-catalog%fractional_translation(:,operation))))>&
            dg_ow_symmetry_tolerance)cycle
        fixed_residual=catalog%fractional_translation(:,operation)-candidate_center+&
          matmul(real(catalog%integer_rotation(:,:,operation),8),candidate_center)
        fixed_residual=fixed_residual-anint(fixed_residual)
        if(maxval(abs(fixed_residual))>dg_ow_symmetry_tolerance)cycle
        duplicate_rotation=.false.
        do g=1,candidate_count
          if(all(catalog%integer_rotation(:,:,candidate_selected(g))==&
              catalog%integer_rotation(:,:,operation)))then
            duplicate_rotation=.true.;exit
          endif
        enddo
        if(duplicate_rotation)cycle
        candidate_count=candidate_count+1;candidate_selected(candidate_count)=operation
      enddo
      candidate_norm=sum(min(candidate_center,1d0-candidate_center)**2)
      candidate_lex=candidate_center(1)+1d-3*candidate_center(2)+1d-6*candidate_center(3)
      if(candidate_count>best_count.or.(candidate_count==best_count.and.&
          (candidate_norm<best_norm-dg_ow_symmetry_tolerance.or.&
          (abs(candidate_norm-best_norm)<=dg_ow_symmetry_tolerance.and.candidate_lex<best_lex))))then
        best_count=candidate_count;best_norm=candidate_norm;best_lex=candidate_lex
        inversion_operation=candidate_operation;fixed_center=candidate_center
      endif
      enddo;enddo;enddo
    enddo
    deallocate(candidate_selected)
    if(inversion_operation==0)then
      message='fixed-center group requires grid-commensurate identity and inversion';return
    endif
    allocate(selected(size(catalog%operations)));selected_count=1;selected(1)=identity_operation
    do operation=1,size(catalog%operations)
      if(operation==identity_operation)cycle
      do axis=1,3
        translation_grid(axis)=nint(catalog%fractional_translation(axis,operation)*&
          real(dc%lg_tot%num(axis),8))
      enddo
      if(maxval(abs(real(translation_grid,8)/real(dc%lg_tot%num,8)-&
          catalog%fractional_translation(:,operation)-anint(real(translation_grid,8)/&
          real(dc%lg_tot%num,8)-catalog%fractional_translation(:,operation))))>&
          dg_ow_symmetry_tolerance)cycle
      fixed_residual=catalog%fractional_translation(:,operation)-fixed_center+&
        matmul(real(catalog%integer_rotation(:,:,operation),8),fixed_center)
      fixed_residual=fixed_residual-anint(fixed_residual)
      if(maxval(abs(fixed_residual))>dg_ow_symmetry_tolerance)cycle
      duplicate_rotation=.false.
      do g=1,selected_count
        if(all(catalog%integer_rotation(:,:,selected(g))==catalog%integer_rotation(:,:,operation)))then
          duplicate_rotation=.true.;exit
        endif
      enddo
      if(duplicate_rotation)cycle
      selected_count=selected_count+1;selected(selected_count)=operation
    enddo
    if(selected_count<2.or.selected_count>48)then
      message='fixed-center crystallographic subgroup order is outside [2,48]';return
    endif
    allocate(operations(selected_count),source=catalog%operations(selected(1:selected_count)))
    inversion_present=.false.
    do g=1,selected_count
      if(all(operations(g)%W==reshape([-1,0,0,0,-1,0,0,0,-1],[3,3])))inversion_present=.true.
    enddo
    if(.not.inversion_present)then;message='fixed-center subgroup lost inversion';return;endif
    call build_sawf_operation_index(operations,dg_ow_symmetry_tolerance,operation_index,map_ok,detail)
    if(.not.map_ok)then;message='fixed-center operation index: '//trim(detail);return;endif
    allocate(product_table(selected_count,selected_count));product_table=0
    do g=1,selected_count;do h=1,selected_count
      call lookup_sawf_operation_product(operation_index,operations,g,h,product_table(g,h),map_ok,detail)
      if(.not.map_ok)then;message='fixed-center product: '//trim(detail);return;endif
    enddo;enddo
    group_fingerprint=fingerprint_dg_exact_fragment_symmetry(&
      catalog%integer_rotation(:,:,selected(1:selected_count)),product_table,dg_ow_symmetry_tolerance,&
      catalog%fractional_translation(:,selected(1:selected_count)))
    if(group_fingerprint==0_8)then;message='fixed-center group fingerprint is zero';return;endif
    call MPI_Allreduce(selected_count,selected_count_min,1,MPI_INTEGER,MPI_MIN,dc%icomm_tot,ierr)
    call MPI_Allreduce(selected_count,selected_count_max,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    call MPI_Allreduce(group_fingerprint,fingerprint_min,1,MPI_INTEGER8,MPI_MIN,dc%icomm_tot,ierr)
    call MPI_Allreduce(group_fingerprint,fingerprint_max,1,MPI_INTEGER8,MPI_MAX,dc%icomm_tot,ierr)
    call MPI_Allreduce(fixed_center,center_min,3,MPI_DOUBLE_PRECISION,MPI_MIN,dc%icomm_tot,ierr)
    call MPI_Allreduce(fixed_center,center_max,3,MPI_DOUBLE_PRECISION,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.selected_count_min/=selected_count_max.or.&
        fingerprint_min/=fingerprint_max.or.maxval(abs(center_max-center_min))>dg_ow_symmetry_tolerance)then
      message='fixed-center subgroup receipts differ across ranks';return
    endif
    allocate(all_ids(nlocal,nproc),target_ids(nlocal,selected_count))
    call MPI_Allgather(local_ids,nlocal,MPI_INTEGER8,all_ids,nlocal,MPI_INTEGER8,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fixed-center owner gather failed';return;endif
    do g=1,selected_count
      operation=selected(g)
      call build_dg_pointwise_affine_owner_map(dc%lg_tot%num,local_ids,all_ids,&
        catalog%integer_rotation(:,:,operation),catalog%fractional_translation(:,operation),&
        dg_ow_symmetry_tolerance,mapped_ids,mapped_owner,mapped_local,mapped_wrap,map_ok,detail)
      if(.not.map_ok)then;message='fixed-center point map: '//trim(detail);return;endif
      target_ids(:,g)=int(mapped_owner,8)*int(nlocal,8)+int(mapped_local,8)
      deallocate(mapped_ids,mapped_owner,mapped_local,mapped_wrap)
    enddo
    if(rank==0)write(*,'(a,i0,a,l1,a,3(es12.4,1x),a,i0)')&
      '[OW-GS-DIAGNOSTIC] fixed_center_group_order=',selected_count,&
      ' inversion=',inversion_present,' center_fractional=',fixed_center,&
      ' fingerprint=',group_fingerprint
    ok=.true.;message=''
  end subroutine prepare_ow_fixed_center_group

  subroutine project_ow_exact_global_group(metric,hamiltonian,position,velocity,promoted,group_order,&
      fixed_center,fixed_group_fingerprint,workspace_peak_bytes,&
      pre_projection_defect,post_projection_defect,ok,message)
    complex(8),intent(inout)::metric(:,:),hamiltonian(:,:),position(:,:,:),velocity(:,:,:)
    logical,intent(out)::promoted,ok
    integer,intent(out)::group_order
    real(8),intent(out)::fixed_center(3)
    integer(8),intent(out)::fixed_group_fingerprint,workspace_peak_bytes
    real(8),intent(out)::pre_projection_defect,post_projection_defect
    character(*),intent(out)::message
    type(t_sawf_crystallographic_catalog)::catalog
    real(8),allocatable::fractional_positions(:,:),rotations(:,:,:)
    integer,allocatable::species(:),product_table(:,:),selected_operations(:),point_map(:),&
      orbit_operations(:),rank_fragment(:),source_to_target(:),fragment_maps(:,:),&
      fragment_orbit(:),orbit_representative(:)
    integer(8),allocatable::core_symmetry_map(:,:),all_core_ids(:,:),mapped_core_ids(:)
    integer,allocatable::mapped_core_owner(:),mapped_core_local(:),mapped_core_wrap(:,:)
    complex(8),allocatable::raw(:,:,:),symmetry_overlap(:,:,:),&
      representation(:,:,:),metric_inverse(:,:),&
      scalars(:,:,:),vectors(:,:,:,:),&
      projected_scalars(:,:,:),projected_vectors(:,:,:,:)
    real(8)::lattice_inverse(3,3),determinant,raw_unitarity_defect,unitarity_defect,closure_defect,&
      cell_length(3),inversion_center(3),center_fractional(3),fixed_residual(3)
    integer::rank,nproc,ierr,nwann,ncore,atom,operation,inversion_operation,&
      axis,translation_grid(3),&
      g,h,k,selected_count,product_rotation(3,3),source_rank,target_rank,&
      source_fragment,target_fragment,base_operation,max_targets,&
      representative_fragment,representative_rank,valid_operation_count,&
      point,source_grid(3),mapped_grid(3),location(2),mapping_bad,global_mapping_bad
    real(8)::relative_cartesian_rotation(3,3)
    logical::inverse_ok,representation_ok,grid_ok,fragment_ok,center_available
    real(8)::grid_residual,center_grid(3),mapped_coordinate
    character(256)::detail

    promoted=.false.;group_order=0;fixed_center=0d0;fixed_group_fingerprint=0_8
    workspace_peak_bytes=0_8;ok=.false.;message='';pre_projection_defect=0d0;post_projection_defect=0d0
    call MPI_Comm_rank(dc%icomm_tot,rank,ierr);call MPI_Comm_size(dc%icomm_tot,nproc,ierr)
    nwann=size(metric,1)
    if(size(metric,2)/=nwann.or.any(shape(hamiltonian)/=[nwann,nwann]).or.&
        any(shape(position)/=[3,nwann,nwann]).or.any(shape(velocity)/=[3,nwann,nwann]))then
      message='invalid exact global group projection dimensions';return
    end if
    call invert_ow_lattice(dc%system_tot%primitive_a,lattice_inverse,determinant,inverse_ok)
    if(.not.inverse_ok)then;message='global inversion lattice is singular';return;end if
    allocate(fractional_positions(3,dc%system_tot%nion),species(dc%system_tot%nion))
    do atom=1,dc%system_tot%nion
      fractional_positions(:,atom)=modulo(matmul(lattice_inverse,dc%system_tot%Rion(:,atom)),1d0)
      species(atom)=dc%system_tot%kion(atom)
    end do
    call load_sawf_crystallographic_catalog_auto(dc%system_tot%primitive_a,fractional_positions,&
      species,dg_ow_symmetry_tolerance,catalog,representation_ok,detail)
    if(.not.representation_ok)then;message='global inversion catalog: '//trim(detail);return;end if
    inversion_operation=0
    do operation=1,size(catalog%operations)
      if(any(catalog%integer_rotation(:,:,operation)/=&
          reshape([-1,0,0,0,-1,0,0,0,-1],[3,3])))cycle
      do axis=1,3
        translation_grid(axis)=nint(catalog%fractional_translation(axis,operation)*&
          real(dc%lg_tot%num(axis),8))
      end do
      if(maxval(abs(real(translation_grid,8)/real(dc%lg_tot%num,8)-&
          catalog%fractional_translation(:,operation)-anint(real(translation_grid,8)/&
          real(dc%lg_tot%num,8)-catalog%fractional_translation(:,operation))))>&
          dg_ow_symmetry_tolerance)cycle
      inversion_operation=operation
      exit
    end do
    if(inversion_operation==0)then;ok=.true.;return;end if
    center_fractional=0.5d0*catalog%fractional_translation(:,inversion_operation)
    allocate(selected_operations(size(catalog%operations)));selected_count=0
    do operation=1,size(catalog%operations)
      do axis=1,3
        translation_grid(axis)=nint(catalog%fractional_translation(axis,operation)*&
          real(dc%lg_tot%num(axis),8))
      end do
      if(maxval(abs(real(translation_grid,8)/real(dc%lg_tot%num,8)-&
          catalog%fractional_translation(:,operation)-anint(real(translation_grid,8)/&
          real(dc%lg_tot%num,8)-catalog%fractional_translation(:,operation))))>&
          dg_ow_symmetry_tolerance)cycle
      fixed_residual=catalog%fractional_translation(:,operation)-center_fractional+&
        matmul(real(catalog%integer_rotation(:,:,operation),8),center_fractional)
      fixed_residual=fixed_residual-anint(fixed_residual)
      if(maxval(abs(fixed_residual))>dg_ow_symmetry_tolerance)cycle
      do g=1,selected_count
        if(all(catalog%integer_rotation(:,:,selected_operations(g))==&
            catalog%integer_rotation(:,:,operation)))exit
      end do
      if(g<=selected_count)cycle
      selected_count=selected_count+1;selected_operations(selected_count)=operation
    end do
    if(selected_count<2)then;message='global exact point group contains only identity';return;end if
    allocate(raw(nwann,nwann,selected_count),&
      product_table(selected_count,selected_count),rotations(3,3,selected_count));raw=(0d0,0d0)
    do g=1,selected_count
      operation=selected_operations(g);rotations(:,:,g)=catalog%operations(operation)%R
    end do
    do g=1,selected_count;do h=1,selected_count
      ! raw(:,:,g) is the pullback f(r)->f(R_g r+t_g), so matrix
      ! multiplication composes the underlying spatial maps in reverse order.
      product_rotation=matmul(catalog%integer_rotation(:,:,selected_operations(h)),&
        catalog%integer_rotation(:,:,selected_operations(g)));product_table(g,h)=0
      do k=1,selected_count
        if(all(product_rotation==catalog%integer_rotation(:,:,selected_operations(k))))then
          product_table(g,h)=k;exit
        end if
      end do
      if(product_table(g,h)==0)then;message='global exact point rotations are not closed';return;end if
    end do;end do
    fixed_center=modulo(center_fractional,1d0)
    fixed_group_fingerprint=fingerprint_dg_exact_fragment_symmetry(&
      catalog%integer_rotation(:,:,selected_operations(1:selected_count)),product_table,&
      dg_ow_symmetry_tolerance,catalog%fractional_translation(:,selected_operations(1:selected_count)))
    if(fixed_group_fingerprint==0_8)then;message='fixed-center point-group fingerprint is zero';return;endif
    workspace_peak_bytes=16_8*int(nwann,8)*int(nwann,8)*int(selected_count,8)
    ncore=size(ow_core_weights)
    if(ncore<1.or.size(ow_core_values,1)/=nwann.or.size(ow_core_values,2)/=ncore)then
      message='global exact group distributed Wannier core is unavailable';return
    end if
    allocate(core_symmetry_map(ncore,selected_count),all_core_ids(ncore,nproc))
    call MPI_Allgather(ow_core_ids,ncore,MPI_INTEGER8,all_core_ids,ncore,MPI_INTEGER8,&
      dc%icomm_tot,ierr)
    do g=1,selected_count
      operation=selected_operations(g)
      call build_dg_pointwise_affine_owner_map(dc%lg_tot%num,ow_core_ids,all_core_ids,&
        catalog%integer_rotation(:,:,operation),catalog%fractional_translation(:,operation),&
        dg_ow_symmetry_tolerance,mapped_core_ids,mapped_core_owner,mapped_core_local,&
        mapped_core_wrap,representation_ok,detail)
      if(.not.representation_ok)then;message='full-system pointwise action: '//trim(detail);return;end if
      core_symmetry_map(:,g)=int(mapped_core_owner,8)*int(ncore,8)+int(mapped_core_local,8)
      deallocate(mapped_core_ids,mapped_core_owner,mapped_core_local,mapped_core_wrap)
    end do
    call assemble_dg_distributed_basis_symmetry_overlap(dc%icomm_tot,ow_core_values,&
      ow_core_weights,core_symmetry_map,symmetry_overlap,representation_ok,detail)
    if(.not.representation_ok)then
      message='global exact distributed Wannier action: '//trim(detail);return
    end if
    call invert_ow_metric(metric,metric_inverse,representation_ok,detail)
    if(.not.representation_ok)then;message='global exact action metric inverse: '//trim(detail);return;end if
    do g=1,selected_count
      raw(:,:,g)=matmul(metric_inverse,symmetry_overlap(:,:,g))
    end do
    call build_dg_fragment_group_representation(metric,raw,product_table,dg_ow_symmetry_tolerance,&
      representation,raw_unitarity_defect,unitarity_defect,closure_defect,representation_ok,detail,2d0)
    if(rank==0)write(*,'(a,3(a,es12.4))')&
      '[OW-GS-DIAGNOSTIC] global_exact_group_representation',&
      ' raw_unitarity_defect=',raw_unitarity_defect,&
      ' unitarity_defect=',unitarity_defect,' closure_defect=',closure_defect
    if(.not.representation_ok)then;message='global inversion representation: '//trim(detail);return;end if
    allocate(scalars(nwann,nwann,2),vectors(nwann,nwann,3,2))
    workspace_peak_bytes=max(workspace_peak_bytes,16_8*int(nwann,8)*int(nwann,8)*&
      int(2*selected_count+16,8))
    scalars(:,:,1)=metric;scalars(:,:,2)=hamiltonian
    cell_length=real(dc%lg_tot%num,8)*dc%system_tot%hgs
    inversion_center=0.5d0*catalog%fractional_translation(:,inversion_operation)*cell_length
    do axis=1,3
      vectors(:,:,axis,1)=position(axis,:,:)-inversion_center(axis)*metric
      vectors(:,:,axis,2)=velocity(axis,:,:)
    end do
    call project_dg_fragment_covariant_operators(representation,rotations,scalars,vectors,&
      dg_ow_symmetry_tolerance,projected_scalars,projected_vectors,pre_projection_defect,&
      post_projection_defect,representation_ok,detail,2d0)
    if(.not.representation_ok)then;message='global inversion projection: '//trim(detail);return;end if
    metric=projected_scalars(:,:,1);hamiltonian=projected_scalars(:,:,2)
    do axis=1,3
      position(axis,:,:)=projected_vectors(:,:,axis,1)+inversion_center(axis)*metric
      velocity(axis,:,:)=projected_vectors(:,:,axis,2)
    end do
    promoted=.true.;group_order=selected_count;ok=.true.
  end subroutine project_ow_exact_global_group

  subroutine populate_ow_checkpoint(occupations,condition_number,closure_residual,operator_fingerprint,&
      localization_initial_spread,localization_final_spread,localization_maximum_gradient,&
      localization_iterations,localization_converged,mlwf_input_fingerprint,&
      mlwf_transform_fingerprint,mlwf_spreads,mlwf_coordinator_bytes,mlwf_workspace_peak_bytes,&
      mlwf_coordinator_byte_limit,mlwf_symmetry_receipts,mlwf_canonical,affine_workspace_peak_bytes,&
      affine_integer_rotations,affine_fractional_translations,occupied_subspace_distance,&
      occupied_electron_count_drift,occupied_density_interior_difference,&
      occupied_density_boundary_difference,occupied_closure_before,occupied_closure_after,&
      occupied_selected_edge,occupied_rejected_edge,occupied_cluster_gap,&
      occupied_selected_block_dimension,occupied_adaptation_workspace_peak_bytes)
    real(8),intent(in)::occupations(:),condition_number,closure_residual
    integer(8),intent(in)::operator_fingerprint
    real(8),intent(in)::localization_initial_spread,localization_final_spread,&
      localization_maximum_gradient
    integer,intent(in)::localization_iterations
    logical,intent(in)::localization_converged
    integer(8),intent(in)::mlwf_input_fingerprint,mlwf_transform_fingerprint,&
      mlwf_coordinator_bytes,mlwf_workspace_peak_bytes,mlwf_coordinator_byte_limit
    real(8),intent(in)::mlwf_spreads(3),mlwf_symmetry_receipts(3)
    logical,intent(in)::mlwf_canonical
    integer(8),intent(in)::affine_workspace_peak_bytes
    integer,intent(in)::affine_integer_rotations(:,:,:)
    real(8),intent(in)::affine_fractional_translations(:,:)
    real(8),intent(in)::occupied_subspace_distance,occupied_electron_count_drift,&
      occupied_density_interior_difference,occupied_density_boundary_difference,&
      occupied_closure_before,occupied_closure_after,occupied_selected_edge,&
      occupied_rejected_edge,occupied_cluster_gap
    integer,intent(in)::occupied_selected_block_dimension
    integer(8),intent(in)::occupied_adaptation_workspace_peak_bytes
    integer::rank,i,j,nowned,nbox,nproc,ierr,axis,point,ownership_count,operation
    integer(8)::tail_count8,occupation_hash,redistribution_local_hash,redistribution_global_hash,word
    integer(8)::fixed_center_group_fingerprint,point_projection_workspace_peak_bytes
    integer(8),allocatable::all_tail_ids(:)
    complex(8),allocatable::hrows(:,:),metric(:,:),hamiltonian(:,:),metric_inverse(:,:),zero_hamiltonian(:,:),&
      position(:,:,:),derivative(:,:,:),canonical_momentum(:,:,:),&
      velocity(:,:,:),nonlocal_velocity(:,:,:),&
      residual_rows(:,:)
    real(8),allocatable::coordinates(:,:)
    real(8)::origin(3),cell_length(3),local_residual_norm,global_residual_norm,&
      local_h_norm,global_h_norm,local_s_norm,global_s_norm,published_coefficient_residual
    real(8)::pre_projection_defect,post_projection_defect,inversion_pre_defect,inversion_post_defect
    real(8)::fixed_center_fractional(3)
    integer::promoted_group_order,global_exact_group_order
    logical::ok,global_inversion_promoted,global_exact_group_promoted
    character(256)::message
    call MPI_Comm_rank(dc%icomm_tot,rank,i)
    call MPI_Comm_size(dc%icomm_tot,nproc,ierr)
    nowned=count(ow_basis%center_owner_rank==rank)
    if(size(ow_basis%physical_grid_ids)>huge(nbox)/nproc)&
      error stop 'overlapping-Wannier checkpoint box extent overflow'
    nbox=size(ow_basis%physical_grid_ids)*nproc
    tail_count8=int(nowned,8)*int(nbox,8)
    if(tail_count8>int(huge(nbox),8))error stop 'overlapping-Wannier checkpoint tail extent overflow'
    allocate(all_tail_ids(nbox))
    call MPI_Allgather(ow_basis%physical_grid_ids,size(ow_basis%physical_grid_ids),MPI_INTEGER8,&
      all_tail_ids,size(ow_basis%physical_grid_ids),MPI_INTEGER8,dc%icomm_tot,ierr)
    if(.not.allocated(ow_published_hrows))error stop 'one-shot published Hamiltonian is unavailable'
    allocate(hrows,source=ow_published_hrows)
    allocate(residual_rows(size(hrows,1),size(ow_state%coefficients,2)))
    residual_rows=matmul(hrows,ow_state%coefficients)
    do j=1,size(residual_rows,2)
      residual_rows(:,j)=residual_rows(:,j)-&
        ow_state%eigenvalues(j)*matmul(ow_srows,ow_state%coefficients(:,j))
    enddo
    published_coefficient_residual=0d0
    do j=1,size(residual_rows,2)
      local_residual_norm=sum(abs(residual_rows(:,j))**2)
      local_h_norm=sum(abs(matmul(hrows,ow_state%coefficients(:,j)))**2)
      local_s_norm=sum(abs(matmul(ow_srows,ow_state%coefficients(:,j)))**2)
      call MPI_Allreduce(local_residual_norm,global_residual_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        dc%icomm_tot,ierr)
      call MPI_Allreduce(local_h_norm,global_h_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
      call MPI_Allreduce(local_s_norm,global_s_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
      published_coefficient_residual=max(published_coefficient_residual,&
        sqrt(max(0d0,global_residual_norm))/max(tiny(1d0),sqrt(max(0d0,global_h_norm))+&
        abs(ow_state%eigenvalues(j))*sqrt(max(0d0,global_s_norm))))
    enddo
    if(published_coefficient_residual>dg_dc_gs_final_orbital_tolerance)&
      error stop 'published overlapping-Wannier H0 is inconsistent with accepted GS coefficients'
    allocate(metric(size(ow_srows,2),size(ow_srows,2)),hamiltonian(size(ow_srows,2),size(ow_srows,2)))
    metric=(0d0,0d0);hamiltonian=(0d0,0d0)
    do i=1,size(ow_row_ids)
      metric(int(ow_row_ids(i)),:)=ow_srows(i,:)
      hamiltonian(int(ow_row_ids(i)),:)=hrows(i,:)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,metric,size(metric),MPI_DOUBLE_COMPLEX,MPI_SUM,dc%icomm_tot,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,hamiltonian,size(hamiltonian),MPI_DOUBLE_COMPLEX,MPI_SUM,dc%icomm_tot,ierr)
    call invert_ow_metric(metric,metric_inverse,ok,message)
    if(.not.ok)then;write(0,'(a)')trim(message);error stop 'overlapping-Wannier metric inverse publication gate failed';endif
    allocate(zero_hamiltonian(size(metric,1),size(metric,2)));zero_hamiltonian=(0d0,0d0)
    allocate(coordinates(3,size(ow_core_ids)))
    do point=1,size(ow_core_ids)
      coordinates(1,point)=real(modulo(ow_core_ids(point)-1_8,int(dc%lg_tot%num(1),8)),8)*dc%system_tot%hgs(1)
      coordinates(2,point)=real(modulo((ow_core_ids(point)-1_8)/int(dc%lg_tot%num(1),8),&
        int(dc%lg_tot%num(2),8)),8)*dc%system_tot%hgs(2)
      coordinates(3,point)=real((ow_core_ids(point)-1_8)/&
        (int(dc%lg_tot%num(1),8)*int(dc%lg_tot%num(2),8)),8)*dc%system_tot%hgs(3)
    enddo
    origin=0d0;cell_length=real(dc%lg_tot%num,8)*dc%system_tot%hgs
    call assemble_dg_overlapping_wannier_observables(dc%icomm_tot,size(metric,1),ow_core_ids,&
      ow_core_weights,coordinates,origin,cell_length,'cell_wrapped',ow_core_values,ow_core_gradients,&
      metric,metric_inverse,hamiltonian,zero_hamiltonian,ow_global_grid_count,&
      dg_dc_gs_hermiticity_tolerance,position,derivative,canonical_momentum,velocity,&
      nonlocal_velocity,ownership_count,ok,message)
    if(.not.ok.or.int(ownership_count,8)/=ow_global_grid_count)then
      write(0,'(a)')trim(message)
      error stop 'overlapping-Wannier observable publication gate failed'
    endif
    call project_ow_exact_global_group(metric,hamiltonian,position,velocity,global_exact_group_promoted,&
      global_exact_group_order,fixed_center_fractional,fixed_center_group_fingerprint,&
      point_projection_workspace_peak_bytes,inversion_pre_defect,inversion_post_defect,ok,message)
    global_inversion_promoted=global_exact_group_promoted
    if(.not.ok)then
      write(0,'(a)')trim(message)
      error stop 'overlapping-Wannier exact global inversion publication gate failed'
    end if
    promoted_group_order=global_exact_group_order
    pre_projection_defect=inversion_pre_defect;post_projection_defect=inversion_post_defect
    if(rank==0)write(*,'(a,l1,2(a,es12.4))')&
      '[OW-GS-DIAGNOSTIC] global_inversion_promoted=',global_inversion_promoted,&
      ' pre_projection_defect=',inversion_pre_defect,&
      ' post_projection_defect=',inversion_post_defect
    if(rank==0)write(*,'(a,l1,a,i0)')&
      '[OW-GS-DIAGNOSTIC] global_exact_group_promoted=',global_exact_group_promoted,&
      ' global_exact_group_order=',global_exact_group_order
    ! GP4 already symmetrized the published row-owned pencil.  The dense exact-group
    ! projection above is diagnostic only; do not replace or re-solve the accepted H/S/C.
    write(*,'(a,i0,a,es12.4,a,es12.4)')'[OW-GS-DIAGNOSTIC] promoted_point_group_order=',&
      promoted_group_order,' pre_projection_defect=',pre_projection_defect,&
      ' post_projection_defect=',post_projection_defect
    ow_checkpoint=s_dg_overlapping_wannier_checkpoint()
    ow_checkpoint%basis_generation=ow_basis%generation;ow_checkpoint%geometry_generation=1
    ow_checkpoint%basis_fingerprint=ow_state%basis_fingerprint
    ow_checkpoint%operator_fingerprint=operator_fingerprint
    ow_checkpoint%global_lcfo_fingerprint=ieor(ow_state%basis_fingerprint,&
      ishftc(int(ow_basis%retained_rank,8),17))
    if(ow_checkpoint%global_lcfo_fingerprint==0_8)ow_checkpoint%global_lcfo_fingerprint=1_8
    occupation_hash=int(z'6A09E667F3BCC909',8)
    do i=1,size(occupations)
      word=transfer(occupations(i),word)
      occupation_hash=ieor(ishftc(occupation_hash,7),ieor(word,int(i,8)))
    end do
    ow_checkpoint%occupation_block_fingerprint=occupation_hash
    if(ow_checkpoint%occupation_block_fingerprint==0_8)ow_checkpoint%occupation_block_fingerprint=1_8
    ow_checkpoint%affine_cocycle_fingerprint=ow_symmetry_fingerprint
    if(ow_checkpoint%affine_cocycle_fingerprint==0_8)&
      error stop 'missing affine-cocycle checkpoint provenance'
    redistribution_local_hash=ieor(int(z'BB67AE8584CAA73B',8),int(rank+1,8))
    if(rank==0)then
      do i=1,size(ow_basis%center_owner_rank)
        redistribution_local_hash=ieor(ishftc(redistribution_local_hash,11),&
          int(ow_basis%center_owner_rank(i)+1,8))
      end do
    end if
    do i=1,size(ow_core_ids)
      redistribution_local_hash=ieor(ishftc(redistribution_local_hash,13),ow_core_ids(i))
    end do
    do i=1,size(all_tail_ids)
      redistribution_local_hash=ieor(ishftc(redistribution_local_hash,17),all_tail_ids(i))
    end do
    call MPI_Allreduce(redistribution_local_hash,redistribution_global_hash,1,MPI_INTEGER8,&
      MPI_BXOR,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)error stop 'redistribution checkpoint provenance reduction failed'
    ow_checkpoint%redistribution_fingerprint=redistribution_global_hash
    if(ow_checkpoint%redistribution_fingerprint==0_8)ow_checkpoint%redistribution_fingerprint=1_8
    call compute_dg_overlapping_wannier_matrix_fingerprints(dc%icomm_tot,ow_row_ids,hrows,&
      position(:,ow_row_ids,:),velocity(:,ow_row_ids,:),ow_checkpoint%hamiltonian_fingerprint,&
      ow_checkpoint%observable_fingerprint,ok)
    if(.not.ok)error stop 'overlapping-Wannier matrix fingerprint publication gate failed'
    ow_checkpoint%field_coupling_convention='cell_wrapped_length_velocity'
    ow_checkpoint%mlwf_backend='wannier90';ow_checkpoint%mlwf_version='3.1.0'
    ow_checkpoint%mlwf_input_fingerprint=mlwf_input_fingerprint
    ow_checkpoint%mlwf_transform_fingerprint=mlwf_transform_fingerprint
    ow_checkpoint%mlwf_spreads=mlwf_spreads
    ow_checkpoint%mlwf_coordinator_bytes=mlwf_coordinator_bytes
    ow_checkpoint%mlwf_workspace_peak_bytes=mlwf_workspace_peak_bytes
    ow_checkpoint%mlwf_coordinator_byte_limit=mlwf_coordinator_byte_limit
    ow_checkpoint%mlwf_symmetry_receipts=mlwf_symmetry_receipts
    ow_checkpoint%mlwf_canonical=mlwf_canonical
    ow_checkpoint%affine_group_order=size(affine_integer_rotations,3)
    ow_checkpoint%translation_subgroup_order=0
    do operation=1,ow_checkpoint%affine_group_order
      if(all(affine_integer_rotations(:,:,operation)==reshape([1,0,0,0,1,0,0,0,1],[3,3])))&
        ow_checkpoint%translation_subgroup_order=ow_checkpoint%translation_subgroup_order+1
    end do
    if(ow_checkpoint%translation_subgroup_order<=0.or.&
       mod(ow_checkpoint%affine_group_order,ow_checkpoint%translation_subgroup_order)/=0)&
      error stop 'invalid affine factor orders at checkpoint publication'
    ow_checkpoint%point_cogroup_order=ow_checkpoint%affine_group_order/&
      ow_checkpoint%translation_subgroup_order
    ow_checkpoint%fixed_center_group_order=global_exact_group_order
    ow_checkpoint%fixed_center_inversion_present=global_inversion_promoted
    ow_checkpoint%fixed_center_fractional=fixed_center_fractional
    ow_checkpoint%fixed_center_group_fingerprint=fixed_center_group_fingerprint
    ow_checkpoint%affine_proof_workspace_peak_bytes=affine_workspace_peak_bytes
    ow_checkpoint%point_projection_workspace_peak_bytes=point_projection_workspace_peak_bytes
    ow_checkpoint%occupied_subspace_distance=occupied_subspace_distance
    ow_checkpoint%occupied_electron_count_drift=occupied_electron_count_drift
    ow_checkpoint%occupied_density_interior_difference=occupied_density_interior_difference
    ow_checkpoint%occupied_density_boundary_difference=occupied_density_boundary_difference
    ow_checkpoint%occupied_density_interior_tolerance=dg_dc_gs_final_density_tolerance
    ow_checkpoint%occupied_density_boundary_tolerance=10d0*dg_dc_gs_final_density_tolerance
    ow_checkpoint%occupied_closure_before=occupied_closure_before
    ow_checkpoint%occupied_closure_after=occupied_closure_after
    ow_checkpoint%occupied_selected_edge=occupied_selected_edge
    ow_checkpoint%occupied_rejected_edge=occupied_rejected_edge
    ow_checkpoint%occupied_cluster_gap=occupied_cluster_gap
    ow_checkpoint%occupied_selected_block_dimension=occupied_selected_block_dimension
    ow_checkpoint%occupied_adaptation_workspace_peak_bytes=occupied_adaptation_workspace_peak_bytes
    allocate(ow_checkpoint%center_owner,source=ow_basis%center_owner_rank)
    allocate(ow_checkpoint%overlap,source=ow_srows)
    allocate(ow_checkpoint%hamiltonian0,source=hrows)
    allocate(ow_checkpoint%position,source=position(:,ow_row_ids,:))
    allocate(ow_checkpoint%velocity,source=velocity(:,ow_row_ids,:))
    allocate(ow_checkpoint%overlap_row_ids,source=ow_row_ids)
    allocate(ow_checkpoint%coefficients,source=ow_state%coefficients)
    allocate(ow_checkpoint%occupations,source=occupations)
    allocate(ow_checkpoint%core_physical_ids,source=ow_core_ids)
    allocate(ow_checkpoint%density,source=ow_state%density)
    allocate(ow_checkpoint%tail_center(nowned),ow_checkpoint%tail_generation(nowned),&
      ow_checkpoint%tail_offsets(nowned+1),ow_checkpoint%tail_physical_ids(int(tail_count8)))
    i=0;ow_checkpoint%tail_offsets(1)=1
    do j=1,size(ow_basis%center_owner_rank)
      if(ow_basis%center_owner_rank(j)/=rank)cycle
      i=i+1;ow_checkpoint%tail_center(i)=j;ow_checkpoint%tail_generation(i)=ow_basis%generation
      ow_checkpoint%tail_offsets(i+1)=i*nbox+1
      ow_checkpoint%tail_physical_ids((i-1)*nbox+1:i*nbox)=all_tail_ids
    enddo
    ow_checkpoint%density_residual=ow_result%density_residual
    ow_checkpoint%unmixed_density_residual=ow_result%unmixed_density_residual
    ow_checkpoint%coefficient_residual=ow_result%coefficient_residual
    ow_checkpoint%orthogonality_defect=ow_result%orthogonality_defect
    ow_checkpoint%metric_condition=condition_number
    ow_checkpoint%charge_error=ow_result%integrated_charge-ow_result%trace_charge
    ow_checkpoint%density_tolerance=dg_dc_gs_final_density_tolerance
    ow_checkpoint%coefficient_tolerance=dg_dc_gs_final_orbital_tolerance
    ow_checkpoint%orthogonality_tolerance=10d0*dg_dc_gs_final_orbital_tolerance
    ow_checkpoint%charge_tolerance=dg_dc_gs_electron_count_tolerance
    ow_checkpoint%condition_limit=1d0/dg_dc_metric_rank_tolerance
    ow_checkpoint%symmetry_closure_residual=closure_residual
    ow_checkpoint%symmetry_tolerance=dg_ow_symmetry_tolerance
    ow_checkpoint%localization_initial_spread=localization_initial_spread
    ow_checkpoint%localization_final_spread=localization_final_spread
    ow_checkpoint%localization_maximum_gradient=localization_maximum_gradient
    ow_checkpoint%localization_iterations=localization_iterations
    ow_checkpoint%localization_converged=localization_converged
    ow_checkpoint%gs_acceptance_receipts=[0d0,0d0,&
      abs(sum(occupations)-dc%elec_num_tot)/dg_dc_gs_electron_count_tolerance,&
      ow_result%coefficient_residual/dg_dc_gs_final_orbital_tolerance,&
      ow_result%density_residual/dg_dc_gs_final_density_tolerance,&
      ow_result%coefficient_residual/dg_dc_gs_final_orbital_tolerance,&
      merge(inversion_post_defect/dg_ow_symmetry_tolerance,0d0,global_inversion_promoted),&
      0d0,0d0,max(post_projection_defect,&
        merge(inversion_post_defect,0d0,global_inversion_promoted))/dg_ow_symmetry_tolerance,&
      published_coefficient_residual/dg_dc_gs_final_orbital_tolerance]
    ow_checkpoint%gs_acceptance_tolerance=1d0
    ow_checkpoint%accepted=ow_result%converged.and.localization_converged
  end subroutine

  subroutine dg_dc_update_potential_from_density(density_arg,ok,message)
    real(8),intent(in)::density_arg(:,:,:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::is,ix,iy,iz,jx,jy,jz
    logical::global_layout
    global_layout=all(shape(density_arg(:,:,:,1))==dc%lg_tot%num)
    do is=1,system%nspin
      do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3);do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2);do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
        if(global_layout)then;jx=ix;jy=iy;jz=iz
        else
          jx=ix-lbound(dc%rho_tot_s(is)%f,1)+1;jy=iy-lbound(dc%rho_tot_s(is)%f,2)+1
          jz=iz-lbound(dc%rho_tot_s(is)%f,3)+1
        endif
        dc%rho_tot_s(is)%f(ix,iy,iz)=density_arg(jx,jy,jz,is)
      enddo;enddo;enddo
    enddo
    do is=1,system%nspin
      do iz=mg%is(3),mg%ie(3);do iy=mg%is(2),mg%ie(2);do ix=mg%is(1),mg%ie(1)
        if(global_layout)then;jx=dc%jxyz_tot(ix,1);jy=dc%jxyz_tot(iy,2);jz=dc%jxyz_tot(iz,3)
        else
          jx=ix-lbound(rho_s(is)%f,1)+1;jy=iy-lbound(rho_s(is)%f,2)+1;jz=iz-lbound(rho_s(is)%f,3)+1
        endif
        rho_s(is)%f(ix,iy,iz)=density_arg(jx,jy,jz,is)
      enddo;enddo;enddo
    enddo
    dc%rho_tot%f=0d0
    do is=1,system%nspin;dc%rho_tot%f=dc%rho_tot%f+dc%rho_tot_s(is)%f;enddo
    call finish_dg_dc_potential_update(ok,message)
  end subroutine dg_dc_update_potential_from_density

  subroutine dg_dc_update_potential_from_distributed_density(point_ids_arg,density_arg,ok,message)
    integer(8),intent(in)::point_ids_arg(:)
    real(8), intent(in) :: density_arg(:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: is,ix,iy,iz,jx,jy,jz,p
    call load_dg_hybrid_distributed_dc_density(dc,point_ids_arg,density_arg,ok,message)
    if(.not.ok)return
    do is=1,system%nspin;rho_s(is)%f=0d0;enddo
    do p=1,size(point_ids_arg)
      jx=int(modulo(point_ids_arg(p)-1_8,int(dc%lg_tot%num(1),8)))+1
      jy=int(modulo((point_ids_arg(p)-1_8)/int(dc%lg_tot%num(1),8),int(dc%lg_tot%num(2),8)))+1
      jz=int((point_ids_arg(p)-1_8)/int(dc%lg_tot%num(1)*dc%lg_tot%num(2),8))+1
      ix=findloc(dc%jxyz_tot(:,1),jx,dim=1);iy=findloc(dc%jxyz_tot(:,2),jy,dim=1)
      iz=findloc(dc%jxyz_tot(:,3),jz,dim=1)
      if(ix>0.and.iy>0.and.iz>0)rho_s(1)%f(ix,iy,iz)=density_arg(p)
    enddo
    call finish_dg_dc_potential_update(ok,message)
  end subroutine dg_dc_update_potential_from_distributed_density

  subroutine finish_dg_dc_potential_update(ok,message)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::is,ix,iy,iz
    logical::rho_finite,hartree_finite,vxc_finite,vlocal_finite,nlcc_finite
    real(8)::density_minimum,density_maximum,nlcc_minimum,nlcc_maximum
    type(s_scalar),allocatable::fragment_hartree(:)
    call hartree(dc%lg_tot,dc%mg_tot,dc%info_tot,dc%system_tot,dc%fg_tot,dc%poisson_tot, &
      dc%srg_scalar_tot,stencil,dc%rho_tot,dc%Vh_tot)
    allocate(fragment_hartree(system%nspin))
    do is=1,system%nspin
      call allocate_scalar(mg,fragment_hartree(is))
      dc%vloc_tot(is)%f=dc%Vh_tot%f
    enddo
    call calc_vlocal_fragment_dcdft(system%nspin,mg,fragment_hartree,dc)
    call exchange_correlation(system,xc_func,mg,srg_scalar,srg,rho_s,&
      pp,ppn,info,spsi,stencil,Vxc,energy%E_xc)
    call update_vlocal(mg,system%nspin,fragment_hartree(1),Vpsl,Vxc,v_local)
    dg_gs_potential_epoch=dg_gs_potential_epoch+1_8
    rho_finite=.true.;hartree_finite=.true.;vxc_finite=.true.;vlocal_finite=.true.
    density_minimum=huge(1d0);density_maximum=-huge(1d0)
    nlcc_finite=all(ieee_is_finite(ppn%rho_nlcc))
    nlcc_minimum=minval(ppn%rho_nlcc)
    nlcc_maximum=maxval(ppn%rho_nlcc)
    do is=1,system%nspin
    do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)
    do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
    do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
      rho_finite=rho_finite.and.ieee_is_finite(dc%rho_tot%f(ix,iy,iz))
      hartree_finite=hartree_finite.and.ieee_is_finite(dc%Vh_tot%f(ix,iy,iz))
      if(ieee_is_finite(dc%rho_tot%f(ix,iy,iz)))then
        density_minimum=min(density_minimum,dc%rho_tot%f(ix,iy,iz))
        density_maximum=max(density_maximum,dc%rho_tot%f(ix,iy,iz))
      endif
    enddo
    enddo
    enddo
    enddo
    do is=1,system%nspin
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      vxc_finite=vxc_finite.and.ieee_is_finite(Vxc(is)%f(ix,iy,iz))
      vlocal_finite=vlocal_finite.and.ieee_is_finite(v_local(is)%f(ix,iy,iz))
    enddo
    enddo
    enddo
    enddo
    ok=rho_finite.and.hartree_finite.and.vxc_finite.and.vlocal_finite
    if(ok) then
      message=''
    else
      write(message,'(a,4(l1,1x),a,2(es12.4,1x),a,l1,a,2(es12.4,1x))')&
        'DG DC GS: finite rho/hartree/vxc/vlocal=',rho_finite,hartree_finite,&
        vxc_finite,vlocal_finite,' density_min/max=',density_minimum,density_maximum,&
        ' nlcc_finite=',nlcc_finite,' nlcc_min/max=',nlcc_minimum,nlcc_maximum
    end if
    do is=1,size(fragment_hartree)
      if(allocated(fragment_hartree(is)%f))deallocate(fragment_hartree(is)%f)
    enddo
    deallocate(fragment_hartree)
  end subroutine finish_dg_dc_potential_update

  integer function canonical_to_dc_index(index,core_count,buffer_count)
    integer, intent(in) :: index,core_count,buffer_count
    if(index<=buffer_count) then
      canonical_to_dc_index=core_count+buffer_count+index
    else if(index<=buffer_count+core_count) then
      canonical_to_dc_index=index-buffer_count
    else
      canonical_to_dc_index=core_count+index-(buffer_count+core_count)
    end if
  end function canonical_to_dc_index

  integer(8) function dg_dc_geometry_fingerprint()
    integer(8) :: local_hash
    integer :: ii,jj
    local_hash=not(0_8)
    call hash_integer(local_hash,dc%lg_tot%num(1))
    call hash_integer(local_hash,dc%lg_tot%num(2))
    call hash_integer(local_hash,dc%lg_tot%num(3))
    call hash_real(local_hash,dc%system_tot%Hvol)
    do jj=1,3
    do ii=1,3
      call hash_real(local_hash,dc%system_tot%primitive_a(ii,jj))
      call hash_real(local_hash,dc%system_tot%primitive_b(ii,jj))
    end do
    end do
    do jj=1,dc%n_frag
    do ii=1,3
      call hash_integer(local_hash,dc%nxyz_domain_frag(ii,jj))
      call hash_integer(local_hash,dc%ixyz_frag(ii,jj))
    end do
    end do
    do jj=1,dc%system_tot%nion
      call hash_integer(local_hash,dc%system_tot%kion(jj))
      do ii=1,3
        call hash_real(local_hash,dc%system_tot%Rion(ii,jj))
      end do
    end do
    if(local_hash==0_8) local_hash=1_8
    dg_dc_geometry_fingerprint=local_hash
  end function dg_dc_geometry_fingerprint

  integer(8) function dg_dc_operator_fingerprint(include_potential)
    logical,intent(in),optional::include_potential
    integer(8) :: local_hash,point_hash,global_potential_hash,point_owner
    integer :: ix,iy,iz,is,ii,jj,kk,mpi_error
    logical::hash_potential
    hash_potential=.true.;if(present(include_potential))hash_potential=include_potential
    local_hash=0_8;global_potential_hash=0_8
    if(hash_potential)then
    do is=1,dc%system_tot%nspin
    do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)
    do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
    do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
      point_owner=int(ix-dc%lg_tot%is(1)+1,8)+int(dc%lg_tot%num(1),8)* &
        (int(iy-dc%lg_tot%is(2),8)+int(dc%lg_tot%num(2),8)* &
        (int(iz-dc%lg_tot%is(3),8)+int(dc%lg_tot%num(3),8)*int(is-1,8)))
      point_hash=not(0_8)
      call hash_integer8(point_hash,point_owner)
      call hash_real(point_hash,dc%vloc_tot(is)%f(ix,iy,iz))
      local_hash=ieor(local_hash,point_hash)
    end do
    end do
    end do
    end do
#ifdef USE_MPI
    call MPI_Allreduce(local_hash,global_potential_hash,1,MPI_INTEGER8,MPI_BXOR,dc%icomm_tot,mpi_error)
    if(mpi_error/=MPI_SUCCESS) stop 'DG DC operator fingerprint reduction failed'
#else
    global_potential_hash=local_hash
#endif
    endif
    local_hash=not(0_8)
    call hash_integer8(local_hash,global_potential_hash)
    call hash_integer(local_hash,dc%system_tot%nspin)
    call hash_real(local_hash,dg_dc_gs_sipg_penalty_factor)
    call hash_real(local_hash,dg_dc_gs_target_lambda)
    call hash_real(local_hash,stencil%coef_lap0)
    do jj=1,3
    do ii=1,4
      call hash_real(local_hash,stencil%coef_lap(ii,jj))
      call hash_real(local_hash,stencil%coef_nab(ii,jj))
    end do
    end do
    call hash_integer(local_hash,pp%lmax)
    call hash_integer(local_hash,pp%nrmax)
    if(allocated(pp%zps)) then
      do ii=1,size(pp%zps)
        call hash_integer(local_hash,pp%zps(ii))
        call hash_real(local_hash,pp%rloc(ii))
        call hash_real(local_hash,pp%rps(ii))
      end do
    end if
    if(allocated(pp%rad)) then
      do jj=1,size(pp%rad,2)
      do ii=1,size(pp%rad,1)
        call hash_real(local_hash,pp%rad(ii,jj))
      end do
      end do
    end if
    if(allocated(pp%radnl)) then
      do jj=1,size(pp%radnl,2)
      do ii=1,size(pp%radnl,1)
        call hash_real(local_hash,pp%radnl(ii,jj))
      end do
      end do
    end if
    if(allocated(pp%nrps))then
      do ii=1,size(pp%nrps);call hash_integer(local_hash,pp%nrps(ii));enddo
    endif
    if(allocated(pp%nrps_ao))then
      do ii=1,size(pp%nrps_ao);call hash_integer(local_hash,pp%nrps_ao(ii));enddo
    endif
    if(allocated(pp%mlps))then
      do ii=1,size(pp%mlps);call hash_integer(local_hash,pp%mlps(ii));enddo
    endif
    if(allocated(pp%nproj))then
      do jj=1,size(pp%nproj,2)
      do ii=lbound(pp%nproj,1),ubound(pp%nproj,1)
        call hash_integer(local_hash,pp%nproj(ii,jj))
      enddo
      enddo
    endif
    if(allocated(pp%inorm))then
      do jj=1,size(pp%inorm,2)
      do ii=lbound(pp%inorm,1),ubound(pp%inorm,1)
        call hash_integer(local_hash,pp%inorm(ii,jj))
      enddo
      enddo
    endif
    if(allocated(pp%vloctbl)) then
      do jj=1,size(pp%vloctbl,2)
      do ii=1,size(pp%vloctbl,1)
        call hash_real(local_hash,pp%vloctbl(ii,jj))
      end do
      end do
    end if
    if(allocated(pp%udvtbl).and.allocated(pp%nrps).and.allocated(pp%nproj)) then
      do kk=1,size(pp%udvtbl,3)
      do ii=1,pp%nrps(kk)
      do jj=0,sum(pp%nproj(:,kk))-1
        call hash_real(local_hash,pp%udvtbl(ii,jj,kk))
      end do
      end do
      end do
    end if
    call hash_character(local_hash,trim(xc))
    if(local_hash==0_8) local_hash=1_8
    dg_dc_operator_fingerprint=local_hash
  end function dg_dc_operator_fingerprint

  subroutine fingerprint_ow_spatial_frame(comm,row_ids,values,point_weights,tolerance,&
      fingerprint,gram_defect,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::values(:,:)
    real(8),intent(in)::point_weights(:),tolerance
    integer(8),intent(out)::fingerprint
    real(8),intent(out)::gram_defect
    logical,intent(out)::ok
    complex(8),allocatable::gram(:,:)
    integer::i,j,p,nstate,ierr,local_bad,global_bad
    integer(8)::local_hash,quantized
    real(8)::quantum,local_max,global_max

    nstate=size(values,1);fingerprint=0_8;gram_defect=huge(1d0);local_bad=0
    if(nstate<1.or.size(values,2)/=size(row_ids).or.size(point_weights)/=size(row_ids).or.&
        .not.ieee_is_finite(tolerance).or.tolerance<=0d0.or.&
        .not.all(ieee_is_finite(real(values))).or..not.all(ieee_is_finite(aimag(values))).or.&
        .not.all(ieee_is_finite(point_weights)).or.any(point_weights<0d0))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ok=.false.;return;endif
    allocate(gram(nstate,nstate),stat=local_bad)
    call MPI_Allreduce(merge(0,1,local_bad==0),global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(gram))deallocate(gram);ok=.false.;return
    endif
    gram=(0d0,0d0)
    do p=1,size(row_ids)
      do j=1,nstate;do i=1,nstate
        gram(i,j)=gram(i,j)+point_weights(p)*conjg(values(i,p))*values(j,p)
      enddo;enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,gram,nstate*nstate,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;deallocate(gram);ok=.false.;return;endif
    gram_defect=0d0
    do j=1,nstate;do i=1,nstate
      if(i==j)then
        gram_defect=max(gram_defect,abs(gram(i,j)-1d0))
      else
        gram_defect=max(gram_defect,abs(gram(i,j)))
      endif
    enddo;enddo
    quantum=100d0*tolerance;local_max=0d0
    if(size(row_ids)>0)local_max=maxval(abs(values))
    call MPI_Allreduce(local_max,global_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS.and.gram_defect<=10d0*tolerance.and.&
      global_max<=0.25d0*real(huge(0_8),8)*quantum)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;deallocate(gram);ok=.false.;return;endif
    local_hash=0_8
    do p=1,size(row_ids);do i=1,nstate
      quantized=nint(real(values(i,p),8)/quantum,8)
      local_hash=ieor(local_hash,ieor(ishftc(row_ids(p),7),ieor(ishftc(int(i,8),17),quantized)))
      quantized=nint(aimag(values(i,p))/quantum,8)
      local_hash=ieor(local_hash,ieor(ishftc(row_ids(p),11),ieor(ishftc(int(i,8),23),quantized)))
    enddo;enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    fingerprint=ieor(fingerprint,int(z'243F6A8885A308D3',8))
    if(fingerprint==0_8)fingerprint=1_8
    ok=ierr==MPI_SUCCESS;deallocate(gram)
  end subroutine fingerprint_ow_spatial_frame

  subroutine fingerprint_ow_w90_matrices(comm,m_matrix,a_matrix,fingerprint,ok)
    integer,intent(in)::comm
    complex(8),intent(in)::m_matrix(:,:,:),a_matrix(:,:)
    integer(8),intent(out)::fingerprint
    logical,intent(out)::ok
    integer::rank,i,j,k,ierr
    call MPI_Comm_rank(comm,rank,ierr)
    ok=ierr==MPI_SUCCESS;fingerprint=int(z'510E527FADE682D1',8)
    if(rank==0.and.ok)then
      call hash_integer(fingerprint,size(m_matrix,1));call hash_integer(fingerprint,size(m_matrix,2))
      call hash_integer(fingerprint,size(m_matrix,3));call hash_integer(fingerprint,size(a_matrix,1))
      call hash_integer(fingerprint,size(a_matrix,2))
      do k=1,size(m_matrix,3);do j=1,size(m_matrix,2);do i=1,size(m_matrix,1)
        call hash_real(fingerprint,real(m_matrix(i,j,k),8))
        call hash_real(fingerprint,aimag(m_matrix(i,j,k)))
      enddo;enddo;enddo
      do j=1,size(a_matrix,2);do i=1,size(a_matrix,1)
        call hash_real(fingerprint,real(a_matrix(i,j),8))
        call hash_real(fingerprint,aimag(a_matrix(i,j)))
      enddo;enddo
      if(fingerprint==0_8)fingerprint=1_8
    endif
    call MPI_Bcast(fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    ok=ok.and.ierr==MPI_SUCCESS.and.fingerprint/=0_8
  end subroutine fingerprint_ow_w90_matrices

  subroutine fingerprint_ow_w90_transform(comm,transform,fingerprint,ok)
    integer,intent(in)::comm
    complex(8),intent(in)::transform(:,:)
    integer(8),intent(out)::fingerprint
    logical,intent(out)::ok
    integer(8)::minimum_hash,maximum_hash
    integer::i,j,ierr
    fingerprint=int(z'9B05688C2B3E6C1F',8)
    call hash_integer(fingerprint,size(transform,1));call hash_integer(fingerprint,size(transform,2))
    do j=1,size(transform,2);do i=1,size(transform,1)
      call hash_real(fingerprint,real(transform(i,j),8));call hash_real(fingerprint,aimag(transform(i,j)))
    enddo;enddo
    if(fingerprint==0_8)fingerprint=1_8
    call MPI_Allreduce(fingerprint,minimum_hash,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    ok=ierr==MPI_SUCCESS
    call MPI_Allreduce(fingerprint,maximum_hash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    ok=ok.and.ierr==MPI_SUCCESS.and.minimum_hash==maximum_hash
    fingerprint=minimum_hash
  end subroutine fingerprint_ow_w90_transform

  subroutine hash_integer(hash,value)
    integer(8), intent(inout) :: hash
    integer, intent(in) :: value
    integer(8) :: bits
    integer :: ibyte
    bits=int(value,8)
    do ibyte=0,7
      call hash_byte(hash,int(ibits(bits,8*ibyte,8)))
    end do
    if(hash==0_8) hash=1_8
  end subroutine hash_integer

  subroutine hash_integer8(hash,value)
    integer(8), intent(inout) :: hash
    integer(8), intent(in) :: value
    integer :: ibyte
    do ibyte=0,7
      call hash_byte(hash,int(ibits(value,8*ibyte,8)))
    end do
    if(hash==0_8) hash=1_8
  end subroutine hash_integer8

  subroutine hash_real(hash,value)
    integer(8), intent(inout) :: hash
    real(8), intent(in) :: value
    integer(8) :: bits
    integer :: ibyte
    bits=transfer(value,bits)
    do ibyte=0,7
      call hash_byte(hash,int(ibits(bits,8*ibyte,8)))
    end do
    if(hash==0_8) hash=1_8
  end subroutine hash_real

  subroutine hash_character(hash,value)
    integer(8), intent(inout) :: hash
    character(*), intent(in) :: value
    integer :: ii
    do ii=1,len(value)
      call hash_byte(hash,iachar(value(ii:ii)))
    end do
  end subroutine hash_character

  subroutine hash_byte(hash,value)
    integer(8), intent(inout) :: hash
    integer, intent(in) :: value
    integer(8), parameter :: polynomial=int(z'C96C5795D7870F42',8)
    integer :: ibit
    hash=ieor(hash,int(iand(value,255),8))
    do ibit=1,8
      if(btest(hash,0)) then
        hash=ieor(shiftr(hash,1),polynomial)
      else
        hash=shiftr(hash,1)
      end if
    end do
  end subroutine hash_byte

end subroutine main_dft

subroutine prepare_dg_hybrid_divided_production_basis(comm_arg,fragment_count_arg,fragment_id_arg,&
    grid_num_arg,fragment_origin_arg,fragment_size_arg,hgs_arg,ncore_arg,nbox_arg,nxy_arg,&
    global_count_arg,physical_ids_arg,core_ids_arg,core_weights_arg,core_values_arg,box_values_arg,&
    wannier_owner_arg,fragment_basis_arg,nwannier_arg,noperation_arg,pencil_maps_arg,raw_partition_arg,&
    reciprocal_lattice_arg,reciprocal_rotations_arg,pw_cutoff_arg,pw_max_arg,tolerance_arg,basis_fingerprint_arg,&
    pw_fingerprint_arg,buffer_fingerprint_arg,fragment_fingerprint_arg,selection_arg,callback_ok,callback_message)
  use dg_hybrid_windowed_pw_types,only:s_dg_hybrid_basis_catalog,s_dg_hybrid_production_selection
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  use dg_hybrid_production_pw_basis,only:analyze_dg_hybrid_lcfo_selection,freeze_dg_hybrid_production_selection
  use dg_hybrid_continuation_state,only:close_dg_hybrid_selection
  use dg_hybrid_window_distribution,only:redistribute_dg_hybrid_fragment_windows
  use dg_hybrid_projected_fragment_pipeline,only:build_dg_hybrid_projected_fragment_basis
  implicit none
  integer,intent(in)::comm_arg,fragment_count_arg,fragment_id_arg,grid_num_arg(3),&
    fragment_origin_arg(3,fragment_count_arg),fragment_size_arg(3,fragment_count_arg),ncore_arg,nbox_arg,&
    nwannier_arg,noperation_arg,pw_max_arg
  integer(8),intent(in)::nxy_arg,global_count_arg,physical_ids_arg(nbox_arg),core_ids_arg(ncore_arg),&
    pencil_maps_arg(ncore_arg,noperation_arg),basis_fingerprint_arg
  real(8),intent(in)::hgs_arg(3),core_weights_arg(ncore_arg),raw_partition_arg(nbox_arg),&
    reciprocal_lattice_arg(3,3),reciprocal_rotations_arg(3,3,noperation_arg),pw_cutoff_arg,tolerance_arg
  complex(8),intent(in)::core_values_arg(nwannier_arg,ncore_arg),box_values_arg(nwannier_arg,nbox_arg)
  integer,intent(in)::wannier_owner_arg(nwannier_arg)
  type(s_dg_hybrid_fragment_basis),intent(out)::fragment_basis_arg
  integer(8),intent(out)::pw_fingerprint_arg,buffer_fingerprint_arg,fragment_fingerprint_arg
  type(s_dg_hybrid_production_selection),intent(out)::selection_arg
  logical,intent(out)::callback_ok
  character(*),intent(out)::callback_message
  integer,allocatable::fragment_ids(:),core_fragment_ids(:),row_action(:,:)
  real(8),allocatable::box_windows(:,:),core_windows(:,:),buffer_windows(:,:),&
    core_coordinates(:,:),buffer_coordinates(:,:),g_vectors(:,:)
  type(s_dg_hybrid_basis_catalog)::pw_catalog
  integer,allocatable::effective_packet_ids(:),closure_parent(:),closure_action(:)
  integer(8)::selection_fingerprint
  integer::point,fragment,grid_x,grid_y,grid_z
  integer(8)::pw_workspace,buffer_window_workspace,fragment_workspace

  callback_ok=.false.;callback_message=''
  allocate(fragment_ids(1),core_fragment_ids(ncore_arg),&
    row_action(ncore_arg,size(pencil_maps_arg,2)),box_windows(1,nbox_arg))
  allocate(core_coordinates(3,ncore_arg),buffer_coordinates(3,nbox_arg))
  fragment_ids(1)=fragment_id_arg
  do point=1,ncore_arg
    core_fragment_ids(point)=0
    grid_x=int(modulo(core_ids_arg(point)-1_8,int(grid_num_arg(1),8)))
    grid_y=int(modulo((core_ids_arg(point)-1_8)/int(grid_num_arg(1),8),int(grid_num_arg(2),8)))
    grid_z=int((core_ids_arg(point)-1_8)/nxy_arg)
    do fragment=1,fragment_count_arg
      if(grid_x<fragment_origin_arg(1,fragment).or.&
        grid_x>=fragment_origin_arg(1,fragment)+fragment_size_arg(1,fragment))cycle
      if(grid_y<fragment_origin_arg(2,fragment).or.&
        grid_y>=fragment_origin_arg(2,fragment)+fragment_size_arg(2,fragment))cycle
      if(grid_z<fragment_origin_arg(3,fragment).or.&
        grid_z>=fragment_origin_arg(3,fragment)+fragment_size_arg(3,fragment))cycle
      core_fragment_ids(point)=fragment
      exit
    enddo
  enddo
  if(any(core_fragment_ids==0))then
    callback_message='divided Hybrid core fragment ownership failed';return
  endif
  row_action=int(pencil_maps_arg)
  if(any(row_action<1).or.any(row_action>int(global_count_arg)))then
    callback_message='divided Hybrid physical row-action reconstruction failed';return
  endif
  box_windows(1,:)=raw_partition_arg
  do point=1,ncore_arg
    core_coordinates(1,point)=real(modulo(core_ids_arg(point)-1_8,int(grid_num_arg(1),8)),8)*hgs_arg(1)
    core_coordinates(2,point)=real(modulo((core_ids_arg(point)-1_8)/int(grid_num_arg(1),8),&
      int(grid_num_arg(2),8)),8)*hgs_arg(2)
    core_coordinates(3,point)=real((core_ids_arg(point)-1_8)/nxy_arg,8)*hgs_arg(3)
  enddo
  do point=1,nbox_arg
    buffer_coordinates(1,point)=real(modulo(physical_ids_arg(point)-1_8,int(grid_num_arg(1),8)),8)*hgs_arg(1)
    buffer_coordinates(2,point)=real(modulo((physical_ids_arg(point)-1_8)/int(grid_num_arg(1),8),&
      int(grid_num_arg(2),8)),8)*hgs_arg(2)
    buffer_coordinates(3,point)=real((physical_ids_arg(point)-1_8)/nxy_arg,8)*hgs_arg(3)
  enddo
  call analyze_dg_hybrid_lcfo_selection(comm_arg,int(global_count_arg),fragment_count_arg,&
    fragment_ids,physical_ids_arg,box_windows,core_ids_arg,core_fragment_ids,core_coordinates,row_action,&
    reciprocal_lattice_arg,reciprocal_rotations_arg,basis_fingerprint_arg,pw_cutoff_arg,pw_max_arg,16,tolerance_arg,&
    core_windows,g_vectors,&
    selection_arg,pw_workspace,pw_fingerprint_arg,callback_ok,callback_message)
  if(.not.callback_ok)return
  call close_dg_hybrid_selection(comm_arg,selection_arg%requested_packet_ids,selection_arg%packet_ids,&
    selection_arg%packet_action,effective_packet_ids,closure_parent,closure_action,selection_fingerprint,&
    callback_ok,callback_message)
  if(.not.callback_ok)return
  call freeze_dg_hybrid_production_selection(comm_arg,selection_arg,effective_packet_ids,pw_catalog,&
    pw_fingerprint_arg,callback_ok,callback_message)
  if(.not.callback_ok)return
  call redistribute_dg_hybrid_fragment_windows(comm_arg,int(global_count_arg),fragment_count_arg,&
    fragment_ids,physical_ids_arg,box_windows,physical_ids_arg,buffer_windows,buffer_window_workspace,&
    buffer_fingerprint_arg,callback_ok,callback_message)
  if(.not.callback_ok)return
  call build_dg_hybrid_projected_fragment_basis(comm_arg,int(global_count_arg),fragment_count_arg,&
    fragment_id_arg,core_ids_arg,core_weights_arg,core_values_arg,core_coordinates,core_windows,&
    physical_ids_arg,box_values_arg,buffer_coordinates,buffer_windows,pw_catalog,g_vectors,wannier_owner_arg,&
    16,tolerance_arg,basis_fingerprint_arg,fragment_basis_arg,fragment_workspace,&
    fragment_fingerprint_arg,callback_ok,callback_message)
end subroutine prepare_dg_hybrid_divided_production_basis
