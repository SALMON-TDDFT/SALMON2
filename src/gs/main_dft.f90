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
#if defined(USE_MPI) && defined(USE_SCALAPACK)
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
  dg_dc_seed_mode,dg_dc_seed_directory,dg_fragment_wf_checkpoint_mode,&
  dg_fragment_wf_checkpoint_directory,dg_fragment_w90_initial_projection,&
  dg_hybrid_symmetry_energy_window,temperature,&
  dg_hybrid_fragment_cg_steps
use dg_dc_seed_checkpoint,only:s_dg_dc_seed_contract,s_dg_dc_seed_payload,&
  DG_DC_SEED_ABSENT,DG_DC_SEED_VALID,DG_DC_SEED_INVALID,&
  build_dg_dc_seed_contract,probe_dg_dc_seed,&
  read_dg_dc_seed,write_dg_dc_seed,restore_dg_dc_seed_payload,resolve_dg_dc_seed_mode
use dg_canonical_pp_fingerprint,only:canonical_pp_fingerprint,canonical_pp_digest,canonical_pp_valence_sum
use dg_overlapping_wannier_construction, only: s_dg_overlapping_wannier_construction, &
  construct_dg_overlapping_wannier_basis,verify_dg_overlapping_wannier_periodic_closure,&
  replicate_dg_fragment_wannier_representative,verify_dg_fragment_wannier_streaming_closure,&
  verify_dg_fragment_center_orbit,verify_dg_uniform_fragment_target_rank
use dg_overlapping_wannier_construction, only: orthonormalize_dg_distributed_seed_space
use dg_overlapping_wannier_construction, only: compose_dg_buffered_orbital_tile_to_physical_grid
#ifdef USE_EIGENEXA
use dg_overlapping_wannier_construction, only: measure_dg_rank_fixed_symmetry_residuals_eigenexa
use dg_overlapping_wannier_construction, only: build_dg_group_averaged_occupied_candidates_eigenexa
use dg_overlapping_wannier_construction, only: build_dg_cocycle_averaged_occupied_candidates_eigenexa
use dg_overlapping_wannier_construction, only: split_dg_translation_character_sector_eigenexa
#endif
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
  validate_dg_factored_point_cogroup_gauge,&
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
use dg_overlapping_wannier_nonlocal,only:&
  assemble_dg_overlapping_wannier_nonlocal_rows,collect_dg_overlapping_wannier_projector_overlaps
use dg_overlapping_wannier_scf, only: s_dg_overlapping_wannier_scf_state, &
  s_dg_overlapping_wannier_scf_result, &
  compute_dg_overlapping_wannier_scf_fingerprint,mix_dg_overlapping_wannier_density_history
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
  s_dg_hybrid_fragment_candidate_catalog,fragment_seed,fragment_pw
use dg_hybrid_divided_operator,only:freeze_dg_hybrid_single_owner_payload
use dg_hybrid_production_face_traces,only:s_dg_hybrid_production_face_trace,&
  freeze_dg_hybrid_basis_directory,materialize_dg_hybrid_production_face_collection,&
  assemble_dg_hybrid_production_interface_component_rows,materialize_dg_hybrid_production_interior,&
  reconstruct_dg_hybrid_production_interface_state,reconstruct_dg_hybrid_production_interface_actions
use dg_hybrid_broken_volume,only:assemble_dg_hybrid_broken_volume_rows,assemble_dg_hybrid_local_potential_rows
use dg_hybrid_variational_payload,only:s_dg_hybrid_fixed_payload,&
  freeze_dg_hybrid_variational_payload,compose_dg_hybrid_variational_hamiltonian,&
  write_dg_hybrid_variational_payload_bundle
use dg_hybrid_publication_policy,only:s_dg_hybrid_candidate_acceptance,&
  validate_dg_hybrid_v5_publication_rank_policy
use dg_hybrid_terminal_refinement,only:s_dg_hybrid_terminal_refinement_controls,&
  s_dg_hybrid_terminal_refinement_state,s_dg_hybrid_terminal_refinement_receipt,&
  s_dg_hybrid_terminal_operator_guard,initialize_dg_hybrid_terminal_refinement,&
  observe_dg_hybrid_terminal_refinement,initialize_dg_hybrid_terminal_operator_guard,&
  validate_dg_hybrid_terminal_operator_guard
use dg_hybrid_localization_first,only:s_dg_hybrid_localization_receipt,&
  prepare_dg_hybrid_localization_first_seed,build_dg_hybrid_localization_receipt
use plusU_global,only:PLUS_U_ON
use dg_hybrid_projected_fragment_pipeline,only:build_dg_hybrid_projected_local_fragment_basis
use dg_hybrid_interface_continuation,only:s_dg_hybrid_interface_continuation,&
  initialize_dg_hybrid_interface_continuation,accept_dg_hybrid_interface_point
use dg_hybrid_schwarz_state,only:s_dg_hybrid_schwarz_state,initialize_dg_hybrid_schwarz_state,&
  validate_dg_hybrid_schwarz_dynamic_receipt
use dg_hybrid_schwarz_operator,only:s_dg_hybrid_schwarz_schedule,&
  build_dg_hybrid_schwarz_schedule,apply_dg_hybrid_schwarz_hamiltonian,apply_dg_hybrid_schwarz_rows
use dg_hybrid_schwarz_solver,only:advance_dg_hybrid_schwarz_epoch,&
  assign_dg_hybrid_schwarz_occupations
use dg_nonlocal_projector_range,only:s_dg_nonlocal_range_receipt,analyze_dg_nonlocal_projector_range
use dg_hybrid_generalized_eigensystem,only:&
  solve_dg_hybrid_generalized_scalapack,solve_dg_hybrid_generalized_once_and_publish,&
  solve_dg_hybrid_generalized_complete_once
use dg_hybrid_ground_state_types,only:s_dg_hybrid_ground_state,&
  validate_dg_hybrid_ground_state
use rt_dg_hybrid_occupied_checkpoint,only:write_rt_dg_hybrid_occupied_checkpoint
use rt_dg_hybrid_checkpoint_v5,only:s_rt_dg_hybrid_v5_shard,&
  collective_rt_dg_hybrid_publication_precondition,&
  collective_rt_dg_hybrid_publication_mapping_precondition,publish_rt_dg_hybrid_checkpoint_v5,&
  s_rt_dg_hybrid_v5_publication_authorization
use rt_dg_hybrid_refinement_receipt,only:s_rt_dg_hybrid_refinement_receipt,&
  write_rt_dg_hybrid_refinement_receipt
use rt_dg_hybrid_initialization,only:fingerprint_rt_dg_hybrid_scope,&
  fingerprint_rt_dg_hybrid_sparse_structure
use rt_dg_hybrid_system_identity,only:fingerprint_rt_dg_hybrid_system
use rt_dg_hybrid_structural_graph,only:build_rt_dg_hybrid_structural_graph
use rt_dg_hybrid_sparse_projection,only:project_rt_dg_hybrid_point_csr_edges
use rt_dg_hybrid_sparse_exchange,only:s_rt_dg_sparse_exchange,build_rt_dg_sparse_exchange,&
  apply_rt_dg_sparse_rows_tiled
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
use dg_overlapping_wannier_symmetry,only:&
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
  export_dg_w90_replay_bundle,DG_W90_CONSTRAINED
use lcfo_wannier_sawf, only: t_sawf_crystallographic_catalog,t_sawf_symop,&
  load_sawf_crystallographic_catalog_auto
use lcfo_wannier_sawf_dmn,only:t_sawf_dmn_writer,t_sawf_operation_index,&
  begin_sawf_dmn,append_sawf_dmn_operation,finish_sawf_dmn,abort_sawf_dmn,&
  build_sawf_operation_index,lookup_sawf_operation_product,&
  convert_sawf_pullback_to_active_representation
use lcfo_wannier_sawf_band, only: validate_sawf_fragment_symmetry_map,&
  build_sawf_fragment_buffer_point_map
#endif
#ifdef USE_EIGENEXA
#if defined(USE_MPI) && defined(USE_SCALAPACK)
use eigenexa_module, only: init_eigenexa_mod=>init_eigenexa,finalize_eigenexa
#endif
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
#if defined(USE_MPI) && defined(USE_SCALAPACK)
use lcfo_flux
#endif
use lcfo_soi
implicit none
integer :: ix,iy,iz
integer :: Miter,iatom,jj,nspin
integer(8) :: dg_gs_potential_epoch
real(8) :: sum1
character(100) :: comment_line
character(1024) :: variational_payload_capture_prefix

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
type(s_unfold) :: unfold
logical :: rion_update
logical :: flag_opt_conv
logical :: local_basis_route_active
integer :: Miopt, iopt,nopt_max,i
integer :: iter_band_kpt, iter_band_kpt_end, iter_band_kpt_stride
logical :: is_checkpoint_iter, is_shutdown_time
integer :: ilevel_print

#if defined(USE_MPI) && defined(USE_SCALAPACK)
type(s_dg_dc_seed_contract) :: dg_dc_seed_contract
type(s_dg_dc_seed_payload) :: dg_dc_seed_payload
logical :: dg_dc_seed_run_scf,dg_dc_seed_load,dg_dc_seed_publish,dg_dc_seed_fatal
logical :: dg_dc_seed_scf_skipped,dg_dc_seed_ok,dg_dc_seed_collective_ok
integer :: dg_dc_seed_status
integer :: dg_dc_seed_rwf_bounds(14),dg_dc_seed_rho_bounds(6),dg_dc_seed_vloc_bounds(6)
integer(int64) :: dg_dc_seed_publication_id
integer(int64) :: dg_dc_seed_immutable_inputs(6),dg_dc_seed_ownership_map(8)
real(8) :: dg_dc_seed_electron_tolerance
character(512) :: dg_dc_seed_message,dg_dc_seed_probe_message
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
type(s_dg_hybrid_schwarz_state) :: bounded_schwarz_state
type(s_dg_hybrid_schwarz_schedule) :: bounded_schwarz_schedule
real(8),allocatable :: ow_hybrid_divided_total_density(:,:,:)
character(16) :: ow_hybrid_divided_convergence
real(8) :: ow_hybrid_divided_threshold
integer(8) :: divided_solver_fingerprint=0_8
real(8) :: divided_fragment_residual=huge(1d0),&
  divided_fragment_orthogonality=huge(1d0)
type(s_dg_hybrid_fixed_payload) :: bounded_fixed_payload
complex(8),allocatable :: bounded_interior_values(:,:),bounded_local_potential_rows(:,:),&
  bounded_fragment_h(:,:),bounded_fragment_s(:,:)
complex(8),allocatable :: bounded_schwarz_candidate_vectors(:,:)
real(8),allocatable :: bounded_schwarz_candidate_energies(:)
integer(8),allocatable :: bounded_schwarz_candidate_ids(:)
real(8),allocatable :: bounded_core_weights(:)
integer,allocatable :: bounded_basis_fragment(:),bounded_basis_local_slot(:),&
  bounded_basis_generation(:),bounded_interior_fragment(:)
integer(8),allocatable :: bounded_core_ids(:)
integer(8) :: bounded_directory_fingerprint=0_8
integer(8) :: bounded_face_fingerprint=0_8,bounded_mapping_fingerprint=0_8,&
  bounded_candidate_fingerprint=0_8
real(8) :: bounded_interface_scale=1d0
integer :: bounded_last_peer_exchange_count=0,bounded_last_accepted_cg_steps=0
real(8) :: ow_hybrid_symmetry_defect=huge(1d0)
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
logical :: ow_direct_nonlocal_compared=.false.
logical :: ow_projector_stage_diagnosed=.false.

interface
end interface
#endif

#if defined(USE_MPI) && defined(USE_SCALAPACK)
dg_dc_seed_ok=.true.;dg_dc_seed_collective_ok=.true.
if(trim(dg_dc_seed_mode)/='off')then
  dg_dc_seed_ok=yn_dc=='y'.and.yn_dg_dc_overlapping_wannier=='y'.and.&
    trim(theory)=='dft'.and.iperiodic==3.and.yn_spinorbit=='n'.and.yn_opt/='y'.and.&
    .not.PLUS_U_ON.and.yn_hse/='y'.and.yn_fix_func/='y'.and.yn_jm/='y'
  call comm_logical_and(dg_dc_seed_ok,dg_dc_seed_collective_ok,nproc_group_global)
  if(.not.dg_dc_seed_collective_ok)&
    error stop 'DG DC seed mode is enabled outside its supported conventional-DC scope'
endif
#endif

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
call init_dft(nproc_group_global,info,lg,mg,system,stencil,fg,poisson,srg,srg_scalar,ofl,unfold)
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

#if defined(USE_MPI) && defined(USE_SCALAPACK)
dg_dc_seed_run_scf=.true.;dg_dc_seed_load=.false.;dg_dc_seed_publish=.false.
dg_dc_seed_fatal=.false.;dg_dc_seed_scf_skipped=.false.;dg_dc_seed_ok=.true.
dg_dc_seed_collective_ok=.true.
dg_dc_seed_status=DG_DC_SEED_ABSENT;dg_dc_seed_publication_id=0_int64
dg_dc_seed_message='';dg_dc_seed_contract=s_dg_dc_seed_contract()
dg_dc_seed_probe_message=''
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
  dg_dc_seed_probe_message=dg_dc_seed_message
  call resolve_dg_dc_seed_mode(dg_dc_seed_mode,dg_dc_seed_status,&
    dg_dc_seed_run_scf,dg_dc_seed_load,dg_dc_seed_publish,dg_dc_seed_fatal,&
    dg_dc_seed_scf_skipped,dg_dc_seed_ok,dg_dc_seed_message)
  if(dg_dc_seed_status==DG_DC_SEED_INVALID.and.len_trim(dg_dc_seed_probe_message)>0)&
    dg_dc_seed_message=dg_dc_seed_probe_message
case('auto')
  call atomic_create_directory(trim(dg_dc_seed_directory),dc%icomm_tot,dc%id_tot)
  call probe_dg_dc_seed(dc%icomm_tot,trim(dg_dc_seed_directory),dg_dc_seed_contract,&
    dc%system_tot%hvol,dc%elec_num_tot,dg_dc_seed_electron_tolerance,threshold,&
    dg_dc_seed_status,dg_dc_seed_publication_id,dg_dc_seed_message)
  dg_dc_seed_probe_message=dg_dc_seed_message
  call resolve_dg_dc_seed_mode(dg_dc_seed_mode,dg_dc_seed_status,&
    dg_dc_seed_run_scf,dg_dc_seed_load,dg_dc_seed_publish,dg_dc_seed_fatal,&
    dg_dc_seed_scf_skipped,dg_dc_seed_ok,dg_dc_seed_message)
  if(dg_dc_seed_status==DG_DC_SEED_INVALID.and.len_trim(dg_dc_seed_probe_message)>0)&
    dg_dc_seed_message=dg_dc_seed_probe_message
case default
  dg_dc_seed_fatal=.true.;dg_dc_seed_message='unknown DG DC seed mode'
end select
if(dg_dc_seed_fatal)then
  if(dc%id_tot==0)then
    write(error_unit,'(a,a)')'[DG-DC-SEED-ERROR] ',trim(dg_dc_seed_message)
    flush(error_unit)
  endif
  error stop 'DG DC seed is absent, invalid, or incompatible'
endif
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
#endif

Miopt = 0
nopt_max = 1
if(yn_opt=='y') call initialization_opt(Miopt,opt,system,flag_opt_conv,nopt_max,ofl)

call timer_end(LOG_INIT_GS)

#if defined(USE_MPI) && defined(USE_SCALAPACK)
if(yn_dc == 'y' .and. yn_dc_lcfo_wannier == 'y' .and. dc_lcfo_wannier_import_only_requested()) then
  if(comm_is_root(nproc_id_global)) &
    write(*,'(1x,a)') '[DC-LCFO-W90-IMPORT] import-only mode: skip SCF and reuse external Wannier90 outputs'
  call dc_lcfo_wannier_import_only(dc)
  call timer_end(LOG_TOTAL)
  return
end if
#endif

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
#if defined(USE_MPI) && defined(USE_SCALAPACK)
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
                        band, ilevel_print,dg_dc_seed_publish,&
                        dg_dc_seed_electron_tolerance,dc)
endif
#else
call scf_iteration_dft( Miter,rion_update,sum1,  &
                        system,energy,ewald,  &
                        lg,mg,info,poisson,fg,cg,mixing,stencil,srg,srg_scalar, &
                        spsi,shpsi,sttpsi,rho,rho_jm,rho_s,V_local,Vh,Vxc,Vpsl,xc_func, &
                        pp,ppg,ppn,band,ilevel_print,.false.,0d0,dc)
#endif


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
if(yn_out_tm  == 'y'.or. yn_out_tm_bin == 'y'.or.yn_out_gs_sgm_eps=='y') then
   select case(iperiodic)
   case(3)
      call write_k_data(system,stencil)  !need? (probably remove later)
      call write_tm_data(spsi,system,info,mg,stencil,srg,ppg,energy)
   case(0)
     write(*,*) "error: yn_out_tm='y','yn_out_tm_bin='y',yn_out_gs_sgm_eps='y' & iperiodic=0"
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

#if defined(USE_MPI) && defined(USE_SCALAPACK)
      if(yn_dg_dc_overlapping_wannier/='y' .and. &
         (is_checkpoint_iter .or. is_shutdown_time)) then
#else
      if(is_checkpoint_iter .or. is_shutdown_time) then
#endif
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
#if defined(USE_MPI) && defined(USE_SCALAPACK)
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
      if(.not.dg_dc_seed_ok)then
        if(dc%id_tot==0)then
          write(error_unit,'(a,a)')'[DG-DC-SEED-ERROR] ',trim(dg_dc_seed_message)
          flush(error_unit)
        endif
        error stop 'failed to publish converged DG DC seed state'
      endif
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
      error stop 'bare overlapping-Wannier GS is retired; select the divided Hybrid route'
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
#else
  if(yn_dc_lcfo == 'y') then
    if(yn_spinorbit == 'y') then
      call dc_lcfo_soi(lg,mg,system,info,stencil,ppg,energy,v_local,spsi,shpsi,sttpsi,srg,dc)
    else
      call dc_lcfo(lg,mg,system,info,stencil,ppg,energy,v_local,spsi,shpsi,sttpsi,srg,dc)
    end if
  end if
#endif
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

#if defined(USE_MPI) && defined(USE_SCALAPACK)
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
    integer(int64)::pp_fingerprint
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
    pp_fingerprint=canonical_pp_fingerprint(pp)
    if(pp_fingerprint==0_int64)then
      hash=0_int64;return
    endif
    call hash_integer8(hash,pp_fingerprint)
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




  subroutine run_dg_hybrid_continuation_ground_state_for_main
    if(dc%id_tot==0)write(*,'(a)')&
      '[DG-HYBRID-ROUTE] continuation selector uses divided local-plus-terminal LCFO route'
    call run_dg_hybrid_divided_ground_state_for_main
  end subroutine run_dg_hybrid_continuation_ground_state_for_main

  subroutine run_dg_hybrid_divided_ground_state_for_main
    character(1024)::diagnostic_prefix,diagnostic_root_prefix
    integer::diagnostic_status,diagnostic_length
    logical::density_diagnostic
    complex(8),allocatable::diagnostic_box(:,:)
    real(8),allocatable::diagnostic_core(:,:,:),diagnostic_occupations(:),diagnostic_conventional(:),&
      diagnostic_frozen_density(:),diagnostic_frozen_potential(:)
    integer(int64),parameter::dg_hybrid_max_interface_points=1000000_int64
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
    type(s_dg_hybrid_interface_continuation)::interface_continuation
    type(s_dg_hybrid_schwarz_state)::accepted_schwarz_state
    type(s_dg_hybrid_terminal_refinement_controls)::terminal_refinement_controls
    type(s_dg_hybrid_terminal_refinement_state)::terminal_refinement_state
    type(s_dg_hybrid_terminal_refinement_receipt)::terminal_refinement_receipt
    type(s_dg_hybrid_terminal_operator_guard)::terminal_operator_guard
    integer::nproc,rank,ierr,status,p,q,axis,index3(3),raw_grid(3),core_grid(3),global_point_count,&
      local_basis_count,total_basis_count,face_count,initial_count,guard_count,candidate_count,&
      pw_candidate_count,global_column
    integer,allocatable::fragment_ids(:),core_fragment_ids(:),row_action(:,:),basis_counts(:),&
      basis_displacements(:),effective_basis_ids(:),basis_owner(:),basis_fragment(:),&
      fragment_origins(:,:),fragment_sizes(:,:),projector_offsets(:),selected_seeds(:)
    integer(int64)::raw_count,byte_limit,pw_workspace,window_workspace,basis_workspace,&
      pw_fingerprint,window_fingerprint,basis_fingerprint,frame_fingerprint,face_fingerprint,&
      global_basis_fingerprint,global_frame_fingerprint,metric_fingerprint,interface_fingerprint,directory_fingerprint,&
      final_operator_fingerprint,&
      final_state_workspace,final_state_fingerprint,final_solver_workspace,&
      final_solver_fingerprint,final_checkpoint_fingerprint,final_provenance(6),&
      attempted_continuation_fingerprint,terminal_fingerprints(2),terminal_fingerprints_min(2),&
      terminal_fingerprints_max(2),terminal_dynamic_receipt
    integer(int64)::support_fingerprints(3),fragment_wf_publication_id
    integer(int64),allocatable::candidate_grid_ids(:),core_ids(:),gathered_basis_ids(:),projector_grid_ids(:)
    complex(8),allocatable::buffer_candidates(:,:),projector_candidates(:,:),reference_frame(:,:),&
      interior_values(:,:),interior_gradients(:,:,:),interior_kinetic_action(:,:),kinetic_rows(:,:),&
      metric_rows(:,:),local_potential_rows(:,:),nonlocal_rows(:,:),interface_components(:,:,:),&
      interface_rows(:,:),schwarz_coupling_rows(:,:),seed_coefficients(:,:),&
      final_local_potential_rows(:,:),final_hrows(:,:),final_srows(:,:)
    complex(8),allocatable::final_solved_coefficients(:,:)
    real(8),allocatable::core_lower(:,:),core_extent(:,:),atom_positions(:,:),raw_weight(:),&
      raw_gradient(:,:),partition_weight(:),partition_gradient(:,:),box_windows(:,:),&
      core_coordinates(:,:),buffer_coordinates(:,:),core_windows(:,:),buffer_windows(:,:),&
      g_vectors(:,:),core_weights(:),projector_weights(:),unit_potential(:),local_potential(:),&
      initial_density(:),final_occupations(:)
    real(8),allocatable::terminal_density_input(:),terminal_density_output(:),terminal_density_mixed(:),&
      terminal_density_history(:,:),terminal_density_new_history(:,:),terminal_solve_local_potential(:)
    real(8),allocatable::final_solved_eigenvalues(:)
    complex(8),allocatable::projector_support_values(:)
    real(8)::axis_weight(3),axis_gradient(3),coordinate,sum_defect,gradient_defect,&
      denominator,convergence_value,electron_defect,terminal_electron_defect,accepted_interface_scale,&
      terminal_density_change,terminal_total_energy,terminal_previous_total_energy,terminal_energy_change
    real(8)::volume_diagnostics(4),local_potential_diagnostics(2),final_scf_receipts(5),&
      final_residual,final_orthogonality,final_projector_defect
    real(8)::reciprocal_rotation(3,3,1),fragment_lattice(3,3),fragment_reciprocal_lattice(3,3)
    integer,allocatable::payload_owner(:),payload_fragment(:),payload_local_slot(:),payload_generation(:),&
      interior_fragment(:)
    character(8),allocatable::atom_symbols(:)
    integer::scf_iterations,final_state_count,terminal_history_count,terminal_new_history_count
    integer(int64)::continuation_point,continuation_point_limit
    logical::ok,collective_ok,point_ok,diagnostic_ok,terminal_state_ok,fragment_wf_checkpoint_hit,&
      terminal_request_another,terminal_have_previous_energy
    character(512)::message,fragment_wf_checkpoint_reason

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
    if(trim(dg_fragment_wf_checkpoint_mode)/='off')&
      call atomic_create_directory(trim(dg_fragment_wf_checkpoint_directory),dc%icomm_tot,dc%id_tot)
    call build_dg_hybrid_fragment_wannier_from_dc_seed(dc%icomm_tot,dc%icomm_frag,info%icomm_o,&
      dc%i_frag,1,'dgfw',raw_grid,[1,1,1],raw_grid,&
      spsi%rwf,energy%esp,system%rocc,system%hvol,candidate_grid_ids,buffer_candidates,&
      projector_candidates,dg_dc_metric_rank_tolerance,fragment_lattice,fragment_reciprocal_lattice,&
      atom_symbols,atom_positions,wannier_num_iter,dg_ow_localization_gradient_tolerance,&
      byte_limit,fragment_cache,ok,message,initial_projection=dg_fragment_w90_initial_projection,&
      checkpoint_mode=dg_fragment_wf_checkpoint_mode,&
      checkpoint_directory=dg_fragment_wf_checkpoint_directory,&
      dc_seed_publication_id=dg_dc_seed_publication_id,&
      mapping_fingerprint=dg_dc_seed_contract%ownership_fingerprint,&
      immutable_fingerprint=dg_dc_seed_contract%immutable_fingerprint,&
      checkpoint_hit=fragment_wf_checkpoint_hit,&
      checkpoint_publication_id=fragment_wf_publication_id,&
      checkpoint_reason=fragment_wf_checkpoint_reason)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local DC-to-Wannier construction failed'
    endif
    if(rank==0)write(*,'(a,a,a,l1,a,i0,a,a)')'[DG-FRAGMENT-WF] mode=',&
      trim(dg_fragment_wf_checkpoint_mode),' checkpoint_hit=',fragment_wf_checkpoint_hit,&
      ' publication_id=',fragment_wf_publication_id,' reason=',trim(fragment_wf_checkpoint_reason)
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
    call MPI_Allreduce(frame_fingerprint,global_frame_fingerprint,1,MPI_INTEGER8,MPI_BXOR,&
      dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)error stop 'fragment-local reference-frame fingerprint reduction failed'
    global_frame_fingerprint=ieor(global_frame_fingerprint,int(z'3C6EF372FE94F82B',int64))
    if(global_frame_fingerprint==0_int64)global_frame_fingerprint=1_int64
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
    allocate(schwarz_coupling_rows(size(projected_basis%global_ids),total_basis_count))
    schwarz_coupling_rows=(0d0,0d0)
    where(abs(kinetic_rows)>0d0.or.abs(nonlocal_rows)>0d0.or.&
      abs(interface_rows)>0d0.or.abs(metric_rows)>0d0)
      schwarz_coupling_rows=(1d0,0d0)
    end where
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
    if(rank==0)write(*,'(a,i0)')'[DG-HYBRID-DIVIDED] projected_basis_fingerprint=',&
      global_basis_fingerprint
    call freeze_dg_hybrid_single_owner_payload(dc%icomm_tot,dc%n_frag,projected_basis,metric_rows,&
      kinetic_rows,nonlocal_rows,interface_rows,global_basis_fingerprint,metric_fingerprint,&
      interface_fingerprint,fixed_payload,payload_owner,payload_fragment,payload_local_slot,&
      payload_generation,directory_fingerprint,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'fragment-local variational payload freeze failed'
    endif
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
      schwarz_coupling_rows,directory_fingerprint,face_fingerprint,bounded_mapping_fingerprint,&
      bounded_schwarz_schedule,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'immutable Schwarz neighbor schedule construction failed'
    endif

    ! Publish the immutable frame/operator context used by sibling callbacks.
    ! Each callback is MPI_COMM_SELF inside its fragment; only the common
    ! occupation and fixed-density interface continuation communicate over dc%icomm_tot.
    divided_fragment_basis=projected_basis
    bounded_fixed_payload=fixed_payload
    bounded_interior_values=interior_values
    bounded_local_potential_rows=local_potential_rows
    bounded_core_weights=core_weights
    bounded_basis_fragment=payload_fragment
    bounded_basis_local_slot=payload_local_slot
    bounded_basis_generation=payload_generation
    bounded_interior_fragment=interior_fragment
    bounded_core_ids=core_ids
    bounded_directory_fingerprint=directory_fingerprint
    bounded_face_fingerprint=face_fingerprint
    if(allocated(ow_core_ids))deallocate(ow_core_ids)
    if(allocated(ow_core_weights))deallocate(ow_core_weights)
    allocate(ow_core_ids,source=core_ids);allocate(ow_core_weights,source=core_weights)
    ow_global_grid_count=int(global_point_count,int64)
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
    call update_dg_hybrid_divided_potential(initial_density,ok)
    if(.not.ok)error stop 'fixed ordinary-DC density potential update failed'
    diagnostic_prefix=''
    call get_environment_variable('SALMON_DG_DENSITY_DIAGNOSTIC_PREFIX',diagnostic_prefix,&
      length=diagnostic_length,status=diagnostic_status)
    call comm_logical_and(diagnostic_status==0.or.diagnostic_status==1,ok,dc%icomm_tot)
    if(.not.ok)error stop 'density diagnostic environment is invalid or truncated'
    diagnostic_root_prefix=diagnostic_prefix
    call comm_bcast(diagnostic_root_prefix,dc%icomm_tot,0)
    call comm_logical_and(diagnostic_prefix==diagnostic_root_prefix,ok,dc%icomm_tot)
    if(.not.ok)error stop 'density diagnostic prefix differs across ranks'
    density_diagnostic=diagnostic_status==0.and.diagnostic_length>0
    call comm_logical_and(density_diagnostic,collective_ok,dc%icomm_tot)
    call comm_logical_and(.not.density_diagnostic,ok,dc%icomm_tot)
    if(.not.collective_ok.and..not.ok)error stop 'density diagnostic must be enabled on every rank'
    if(density_diagnostic)then
      ! Same V_local as the first DG solve; no density feedback in DC-LCFO.
      call dc_lcfo(lg,mg,system,info,stencil,ppg,energy,v_local,spsi,shpsi,sttpsi,srg,dc,&
        retained_count=dc%nstate_tot,retained_box_contribution=diagnostic_box,&
        retained_occupations=diagnostic_occupations,write_files=.false.,&
        retained_core_density=diagnostic_core)
      allocate(diagnostic_conventional(size(core_ids)),diagnostic_frozen_potential(size(core_ids)))
      do p=1,size(core_ids)
        q=core_selection%core_row_slots(p)-1
        index3(1)=modulo(q,raw_grid(1))+1;q=q/raw_grid(1)
        index3(2)=modulo(q,raw_grid(2))+1;index3(3)=q/raw_grid(2)+1
        diagnostic_conventional(p)=diagnostic_core(index3(1),index3(2),index3(3))
      enddo
      deallocate(diagnostic_core,diagnostic_box,diagnostic_occupations)
      call extract_dg_hybrid_core_local_potential(core_ids,diagnostic_frozen_potential,ok)
      if(.not.ok)error stop 'density diagnostic fixed potential extraction failed'
    endif
    call initialize_dg_hybrid_interface_continuation(dc%icomm_tot,projected_basis%generation,&
      bounded_mapping_fingerprint,mixing%mixrate,interface_continuation,ok,message,full_from_start=.true.)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'DG interface continuation initialization failed'
    endif
    accepted_schwarz_state=bounded_schwarz_state
    accepted_interface_scale=0d0
    continuation_point=0_int64
    point_ok=interface_continuation%rate>=1d0/&
      real(dg_hybrid_max_interface_points-2_int64,8)
    call comm_logical_and(point_ok,collective_ok,dc%icomm_tot)
    if(.not.collective_ok)then
      if(rank==0)write(error_unit,'(a,es16.8,a,i0)')&
        '[DG-HYBRID-DIVIDED] interface increment is too small: rate=',&
        interface_continuation%rate,' maximum_points=',dg_hybrid_max_interface_points
      error stop 'DG interface continuation exceeds the finite production point budget'
    endif
    continuation_point_limit=ceiling(1d0/interface_continuation%rate,kind=int64)+2_int64
    do while(.not.interface_continuation%finished)
      continuation_point=continuation_point+1_int64
      if(continuation_point>continuation_point_limit)then
        bounded_schwarz_state=accepted_schwarz_state
        if(rank==0)write(error_unit,'(a,es16.8)')&
          '[DG-HYBRID-DIVIDED] continuation point limit at accepted lambda=',&
          accepted_interface_scale
        error stop 'DG interface continuation exceeded its defensive point limit'
      endif
      bounded_interface_scale=interface_continuation%lambda
      attempted_continuation_fingerprint=interface_continuation%fingerprint
      call solve_dg_hybrid_schwarz_fragments(int(continuation_point),point_ok)
      point_ok=point_ok.and.ieee_is_finite(divided_fragment_residual).and.&
        ieee_is_finite(divided_fragment_orthogonality).and.&
        ieee_is_finite(bounded_schwarz_state%electron_defect)
      call comm_logical_and(point_ok,collective_ok,dc%icomm_tot)
      if(.not.collective_ok)then
        bounded_schwarz_state=accepted_schwarz_state
        call record_dg_hybrid_interface_continuation_diagnostic(accepted_interface_scale,&
          attempted_continuation_fingerprint,interface_continuation,.false.,diagnostic_ok)
        if(rank==0)write(error_unit,'(a,es16.8,2a)')&
          '[DG-HYBRID-DIVIDED] rejected point; last accepted lambda=',&
          accepted_interface_scale,' ',trim(message)
        error stop 'DG interface continuation point failed collectively'
      endif
      call record_dg_hybrid_interface_continuation_diagnostic(bounded_interface_scale,&
        attempted_continuation_fingerprint,interface_continuation,.true.,diagnostic_ok)
      if(.not.diagnostic_ok)then
        bounded_schwarz_state=accepted_schwarz_state
        if(rank==0)write(error_unit,'(a,es16.8)')&
          '[DG-HYBRID-DIVIDED] diagnostic/acceptance rollback; last accepted lambda=',&
          accepted_interface_scale
        error stop 'DG interface continuation diagnostic certification failed'
      endif
      accepted_schwarz_state=bounded_schwarz_state
      accepted_interface_scale=bounded_interface_scale
    enddo
    scf_iterations=interface_continuation%accepted_steps
    convergence_value=divided_fragment_residual
    electron_defect=bounded_schwarz_state%electron_defect
    if(rank==0)write(*,'(a,2(a,es12.4),a,i0)')'[DG-HYBRID-DIVIDED] selected basis prepared',&
      ' partition_sum_defect=',sum_defect,' partition_gradient_defect=',gradient_defect,&
      ' global_basis_count=',size(projected_basis%global_ids)
    if(rank==0)write(*,'(a,i0,2(a,es12.4))')'[OW-GS] fixed-density DG continuation points=',&
      scf_iterations,' residual=',convergence_value,' electron_defect=',electron_defect

    ! Project the fixed ordinary-DC potential, compose the complete
    ! row-distributed DG operator, and diagonalize it exactly once.
    allocate(final_local_potential_rows(size(projected_basis%global_ids),total_basis_count),&
      final_hrows(size(projected_basis%global_ids),total_basis_count),&
      final_srows(size(projected_basis%global_ids),total_basis_count))
    call extract_dg_hybrid_core_local_potential(core_ids,local_potential,ok)
    if(.not.ok)error stop 'terminal divided Hybrid potential extraction failed'
    call assemble_dg_hybrid_local_potential_rows(dc%icomm_tot,total_basis_count,&
      projected_basis%global_ids,payload_fragment,core_ids,interior_fragment,core_weights,&
      interior_values,local_potential,final_local_potential_rows,local_potential_diagnostics,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'terminal divided Hybrid potential projection failed'
    endif
    call validate_dg_hybrid_schwarz_dynamic_receipt(dc%icomm_tot,bounded_schwarz_state,&
      terminal_dynamic_receipt,terminal_state_ok,message)
    terminal_fingerprints=[interface_continuation%fingerprint,bounded_schwarz_state%fingerprint]
    call MPI_Allreduce(terminal_fingerprints,terminal_fingerprints_min,2,MPI_INTEGER8,MPI_MIN,&
      dc%icomm_tot,ierr)
    terminal_state_ok=terminal_state_ok.and.ierr==MPI_SUCCESS
    call MPI_Allreduce(terminal_fingerprints,terminal_fingerprints_max,2,MPI_INTEGER8,MPI_MAX,&
      dc%icomm_tot,ierr)
    terminal_state_ok=terminal_state_ok.and.ierr==MPI_SUCCESS.and.&
      all(terminal_fingerprints==terminal_fingerprints_min).and.&
      all(terminal_fingerprints==terminal_fingerprints_max).and.&
      interface_continuation%valid.and.interface_continuation%finished.and.&
      interface_continuation%lambda==1d0.and.bounded_interface_scale==1d0.and.&
      accepted_interface_scale==1d0.and.interface_continuation%fingerprint/=0_int64.and.&
      bounded_schwarz_state%valid.and.bounded_schwarz_state%fingerprint/=0_int64.and.&
      interface_continuation%basis_generation==bounded_schwarz_state%basis_generation.and.&
      interface_continuation%mapping_fingerprint==bounded_schwarz_state%mapping_fingerprint
    call comm_logical_and(terminal_state_ok,collective_ok,dc%icomm_tot)
    if(.not.collective_ok)error stop 'terminal LCFO requires a consistent finished lambda-one continuation state'
    final_hrows=bounded_fixed_payload%kinetic_rows+bounded_fixed_payload%nonlocal_rows+&
      bounded_fixed_payload%interface_rows+final_local_potential_rows
    final_srows=bounded_fixed_payload%metric_rows
    call ow_fingerprint_distributed_matrix(dc%icomm_tot,projected_basis%global_ids,final_hrows,&
      final_operator_fingerprint,ok)
    if(.not.ok)error stop 'terminal divided Hybrid operator fingerprint failed'
    if(.not.allocated(bounded_schwarz_state%occupations))&
      error stop 'terminal divided Hybrid occupations are unavailable'
    ! PZHEEVD computes the complete construction-basis spectrum.  The normal
    ! path executes this loop once.  Only a density/energy mismatch requests
    ! one of at most three additional local-potential refinements.
    final_state_count=total_basis_count
    if(final_state_count<1)error stop 'terminal divided Hybrid occupied inventory is empty'
    allocate(final_occupations(final_state_count))
    final_occupations=0d0
    final_occupations(:min(final_state_count,size(bounded_schwarz_state%occupations)))=&
      bounded_schwarz_state%wspin*bounded_schwarz_state%occupations(&
        :min(final_state_count,size(bounded_schwarz_state%occupations)))
    allocate(terminal_density_input,source=initial_density)
    allocate(terminal_density_output(size(initial_density)),terminal_density_mixed(size(initial_density)),&
      terminal_density_history(size(initial_density),2),terminal_density_new_history(size(initial_density),2),&
      terminal_solve_local_potential(size(initial_density)))
    terminal_density_history(:,1)=initial_density;terminal_density_history(:,2)=initial_density
    terminal_density_new_history=terminal_density_history;terminal_history_count=0
    terminal_have_previous_energy=.false.;terminal_previous_total_energy=0d0
    terminal_refinement_controls%maximum_additional_solves=3
    terminal_refinement_controls%density_tolerance=dg_dc_gs_final_density_tolerance
    terminal_refinement_controls%energy_tolerance=ow_hybrid_divided_threshold
    call initialize_dg_hybrid_terminal_refinement(dc%icomm_tot,terminal_refinement_controls,&
      terminal_refinement_state,ok,message)
    if(.not.ok)error stop 'terminal divided Hybrid refinement initialization failed'
    call initialize_dg_hybrid_terminal_operator_guard(dc%icomm_tot,bounded_fixed_payload%metric_rows,&
      bounded_fixed_payload%kinetic_rows,bounded_fixed_payload%nonlocal_rows,&
      bounded_fixed_payload%interface_rows,payload_generation,payload_owner,&
      bounded_fixed_payload%fingerprint,initial_density,terminal_operator_guard,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'terminal divided Hybrid immutable operator guard initialization failed'
    endif
terminal_lcfo_refinement: do
      call validate_dg_hybrid_terminal_operator_guard(dc%icomm_tot,bounded_fixed_payload%metric_rows,&
        bounded_fixed_payload%kinetic_rows,bounded_fixed_payload%nonlocal_rows,&
        bounded_fixed_payload%interface_rows,payload_generation,payload_owner,&
        bounded_fixed_payload%fingerprint,initial_density,terminal_operator_guard,ok,message)
      if(.not.ok)then
        if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
        error stop 'terminal divided Hybrid immutable operator changed before solve'
      endif
      terminal_solve_local_potential=local_potential
      call solve_dg_hybrid_generalized_once_and_publish(dc%icomm_tot,total_basis_count,final_state_count,&
        projected_basis%global_ids,final_hrows,final_srows,dg_dc_gs_final_orbital_tolerance,&
        final_occupations,dc%elec_num_tot,global_basis_fingerprint,metric_fingerprint,&
        final_operator_fingerprint,global_frame_fingerprint,solve_final_dg_hybrid_divided_lcfo,&
        ow_hybrid_ground_state,final_state_workspace,final_state_fingerprint,final_residual,&
        final_orthogonality,final_projector_defect,final_solver_workspace,final_solver_fingerprint,&
        ok,message,electronic_temperature=max(0d0,temperature),&
        occupation_electron_tolerance=dg_dc_gs_electron_count_tolerance,&
        solved_coefficients=final_solved_coefficients,solved_eigenvalues=final_solved_eigenvalues)
      if(.not.ok)then
        if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
        error stop 'terminal divided Hybrid LCFO solve failed'
      endif
      call validate_dg_hybrid_terminal_operator_guard(dc%icomm_tot,bounded_fixed_payload%metric_rows,&
        bounded_fixed_payload%kinetic_rows,bounded_fixed_payload%nonlocal_rows,&
        bounded_fixed_payload%interface_rows,payload_generation,payload_owner,&
        bounded_fixed_payload%fingerprint,initial_density,terminal_operator_guard,ok,message)
      if(.not.ok)then
        if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
        error stop 'terminal divided Hybrid immutable operator changed during solve'
      endif
      call reconstruct_dg_hybrid_terminal_density(projected_basis%global_ids,interior_values,&
        ow_hybrid_ground_state,terminal_density_output,ok,message)
      if(.not.ok)error stop 'terminal divided Hybrid density reconstruction failed'
      if(density_diagnostic.and..not.allocated(diagnostic_frozen_density))then
        if(any(local_potential/=diagnostic_frozen_potential))&
          error stop 'density diagnostic potential changed before first DG solve'
        allocate(diagnostic_frozen_density,source=terminal_density_output)
      endif
      call measure_dg_hybrid_terminal_density_change(dc%icomm_tot,core_weights,terminal_density_input,&
        terminal_density_output,terminal_density_change,ok,message)
      if(.not.ok)error stop 'terminal divided Hybrid density-change measurement failed'
      call evaluate_dg_hybrid_terminal_total_energy(core_ids,core_weights,terminal_density_output,&
        terminal_solve_local_potential,ow_hybrid_ground_state,terminal_total_energy,ok,message)
      if(.not.ok)error stop 'terminal divided Hybrid total-energy evaluation failed'
      if(terminal_have_previous_energy)then
        terminal_energy_change=abs(terminal_total_energy-terminal_previous_total_energy)
      else
        ! A good DC seed may legitimately terminate after its first LCFO solve;
        ! there is no preceding terminal energy against which to form a change.
        terminal_energy_change=0d0
      endif
      call observe_dg_hybrid_terminal_refinement(dc%icomm_tot,terminal_refinement_state,&
        terminal_density_change,terminal_energy_change,.true.,terminal_request_another,&
        terminal_refinement_receipt,ok,message)
      if(.not.ok)error stop 'terminal divided Hybrid refinement observation failed'
      if(rank==0)write(*,'(a,i0,3(a,es16.8),a,l1)')'[OW-GS] terminal LCFO solve=',&
        terminal_refinement_receipt%total_solve_count,' density_change=',terminal_density_change,&
        ' total_energy_change=',terminal_energy_change,' total_energy=',terminal_total_energy,&
        ' converged=',terminal_refinement_receipt%converged
      if(.not.terminal_request_another)exit terminal_lcfo_refinement
      if(dg_dc_gs_density_mix_rate==1d0)then
        terminal_density_mixed=terminal_density_output
        terminal_density_new_history=terminal_density_history
        terminal_density_new_history(:,1)=terminal_density_output;terminal_new_history_count=1;ok=.true.
      else
        call mix_dg_overlapping_wannier_density_history(dc%icomm_tot,dg_dc_gs_density_mix_rate,&
          terminal_density_input,terminal_density_output,terminal_density_history,terminal_history_count,&
          terminal_density_mixed,terminal_density_new_history,terminal_new_history_count,ok,message)
      endif
      if(.not.ok)error stop 'terminal divided Hybrid density-history mixing failed'
      terminal_density_history=terminal_density_new_history
      terminal_history_count=terminal_new_history_count
      terminal_density_input=terminal_density_mixed
      terminal_previous_total_energy=terminal_total_energy;terminal_have_previous_energy=.true.
      call dg_dc_update_potential_from_distributed_density(core_ids,terminal_density_input,ok,message)
      if(.not.ok)error stop 'terminal divided Hybrid refined potential update failed'
      call extract_dg_hybrid_core_local_potential(core_ids,local_potential,ok)
      if(.not.ok)error stop 'terminal divided Hybrid refined potential extraction failed'
      call assemble_dg_hybrid_local_potential_rows(dc%icomm_tot,total_basis_count,&
        projected_basis%global_ids,payload_fragment,core_ids,interior_fragment,core_weights,&
        interior_values,local_potential,final_local_potential_rows,local_potential_diagnostics,ok,message)
      if(.not.ok)error stop 'terminal divided Hybrid refined potential projection failed'
      final_hrows=bounded_fixed_payload%kinetic_rows+bounded_fixed_payload%nonlocal_rows+&
        bounded_fixed_payload%interface_rows+final_local_potential_rows
      call ow_fingerprint_distributed_matrix(dc%icomm_tot,projected_basis%global_ids,final_hrows,&
        final_operator_fingerprint,ok)
      if(.not.ok)error stop 'terminal divided Hybrid refined operator fingerprint failed'
    enddo terminal_lcfo_refinement
    if(density_diagnostic)call write_dg_density_diagnostic(diagnostic_prefix,core_ids,core_weights,&
      initial_density,diagnostic_conventional,diagnostic_frozen_density,terminal_density_output,&
      diagnostic_frozen_potential)
    if(.not.terminal_refinement_receipt%converged.and.&
       .not.terminal_refinement_receipt%publish_last_valid)&
      error stop 'terminal divided Hybrid refinement ended without a publishable state'
    ow_hybrid_ground_state%final_eigensolve_count=terminal_refinement_receipt%total_solve_count
    ow_hybrid_ground_state%additional_refinement_count=&
      terminal_refinement_receipt%additional_refinement_count
    ow_hybrid_ground_state%refinement_converged=terminal_refinement_receipt%converged
    ow_hybrid_ground_state%refinement_exhausted=terminal_refinement_receipt%exhausted
    ow_hybrid_ground_state%terminal_density_change=terminal_refinement_receipt%density_change
    ow_hybrid_ground_state%terminal_energy_change=terminal_refinement_receipt%energy_change
    terminal_electron_defect=abs(sum(ow_hybrid_ground_state%occupations)-dc%elec_num_tot)
    if(rank==0.and.ow_hybrid_ground_state%final_eigensolve_count==1)write(*,'(a,4(a,es16.8))')&
      '[OW-GS] fixed-density/non-self-consistent divided WF+PW LCFO solved once',&
      ' residual=',final_residual,' orthogonality=',final_orthogonality,&
      ' projector=',final_projector_defect,' electron_defect=',terminal_electron_defect
    if(rank==0)write(*,'(a,i0,4(a,es16.8))')&
      '[OW-GS] divided WF+PW terminal LCFO total_solve_count=',&
      ow_hybrid_ground_state%final_eigensolve_count,&
      ' residual=',final_residual,' orthogonality=',final_orthogonality,&
      ' projector=',final_projector_defect,' electron_defect=',terminal_electron_defect
    final_provenance=[pw_fingerprint,window_fingerprint,global_basis_fingerprint,&
      bounded_fixed_payload%fingerprint,final_operator_fingerprint,final_solver_fingerprint]
    final_scf_receipts=[convergence_value,final_residual,final_orthogonality,&
      final_projector_defect,electron_defect]
    call write_rt_dg_hybrid_occupied_checkpoint(dc%icomm_tot,'./overlapping_wannier_occupied.chk',&
      ow_hybrid_ground_state%global_count,ow_hybrid_ground_state%owned_row_ids,&
      ow_hybrid_ground_state%coefficients,ow_hybrid_ground_state%occupations,&
      ow_hybrid_ground_state%eigenvalues,pw_fingerprint,global_basis_fingerprint,&
      final_provenance,final_operator_fingerprint,final_state_fingerprint,final_scf_receipts,&
      max(dg_dc_gs_final_orbital_tolerance,dg_dc_gs_electron_count_tolerance,&
      maxval(final_scf_receipts)),final_checkpoint_fingerprint,ok,message)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED] ',trim(message)
      error stop 'terminal divided Hybrid occupied checkpoint failed'
    endif
    call publish_dg_hybrid_divided_v5(projected_basis%global_ids,payload_owner,payload_generation,&
      core_ids,core_weights,interior_fragment,interior_values,interior_gradients,&
      bounded_fixed_payload%metric_rows,bounded_fixed_payload%kinetic_rows,&
      bounded_fixed_payload%nonlocal_rows,final_local_potential_rows,bounded_fixed_payload%interface_rows,&
      final_hrows,final_solved_coefficients,final_solved_eigenvalues,ow_hybrid_ground_state,&
      global_basis_fingerprint,global_frame_fingerprint,metric_fingerprint,&
      bounded_fixed_payload%interface_fingerprint,interface_continuation%fingerprint,&
      final_operator_fingerprint,final_residual,final_orthogonality,final_projector_defect,&
      terminal_electron_defect,0,.true.,ok,message,terminal_refinement=terminal_refinement_receipt)
    if(.not.ok)then
      if(rank==0)write(error_unit,'(a,a)')'[DG-HYBRID-DIVIDED-V5] ',trim(message)
      error stop 'terminal divided Hybrid v5 publication failed'
    endif
  end subroutine run_dg_hybrid_divided_ground_state_for_main

  subroutine write_dg_density_diagnostic(prefix,ids,weights,seed,conventional,frozen,relaxed,potential)
    character(*),intent(in)::prefix
    integer(8),intent(in)::ids(:)
    real(8),intent(in)::weights(:),seed(:),conventional(:),frozen(:),relaxed(:),potential(:)
    integer::rank,nproc,ierr,unit,status,close_status,p
    logical::ok,all_ok
    character(1200)::filename
    call MPI_Comm_rank(dc%icomm_tot,rank,ierr)
    call MPI_Comm_size(dc%icomm_tot,nproc,ierr)
    write(filename,'(a,a,i8.8)')trim(prefix),'.rank-',rank
    ! Diagnostic-only export, not a restart cache. Never overwrite evidence.
    open(newunit=unit,file=trim(filename),status='new',action='write',iostat=status)
    ok=status==0
    if(ok)then
      write(unit,'(a)',iostat=status)'SALMON_DG_DENSITY_DIAGNOSTIC_V1'
      if(status==0)write(unit,*,iostat=status)nproc,rank,dc%i_frag,size(ids),product(dc%lg_tot%num)
      if(status==0)write(unit,*,iostat=status)dg_dc_seed_publication_id,&
        dg_dc_seed_contract%ownership_fingerprint,dg_dc_seed_contract%immutable_fingerprint
      if(status==0)write(unit,*,iostat=status)dc%elec_num_tot,temperature
      do p=1,size(ids)
        if(status/=0)exit
        write(unit,'(i0,6(1x,es26.17e3))',iostat=status)ids(p),weights(p),seed(p),conventional(p),&
          frozen(p),relaxed(p),potential(p)
      enddo
      if(status==0)write(unit,'(a)',iostat=status)'END_DENSITY_DIAGNOSTIC'
      close(unit,iostat=close_status)
      ok=status==0.and.close_status==0
    endif
    call comm_logical_and(ok,all_ok,dc%icomm_tot)
    if(.not.all_ok)error stop 'density diagnostic export failed; partial files are not reusable'
    if(rank==0)write(*,'(a)')'[DG-DENSITY-DIAGNOSTIC] exported; absolute DC error unavailable'
  end subroutine write_dg_density_diagnostic

  subroutine reconstruct_dg_hybrid_terminal_density(row_ids,basis_values,state,density,ok,message)
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::basis_values(:,:)
    type(s_dg_hybrid_ground_state),intent(in)::state
    real(8),intent(out)::density(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::orbital_values(:)
    integer::i,p
    logical::local_ok,collective_ok

    ok=.false.;message='';density=0d0
    local_ok=state%valid.and.allocated(state%coefficients).and.allocated(state%occupations)
    if(local_ok)local_ok=size(density)==size(basis_values,2).and.&
      size(row_ids)==size(state%coefficients,1).and.&
      size(state%coefficients,2)==size(state%occupations)
    if(local_ok)local_ok=all(row_ids>=1_8).and.all(row_ids<=int(size(basis_values,1),8))
    call comm_logical_and(local_ok,collective_ok,dc%icomm_tot)
    if(.not.collective_ok)then;message='invalid terminal density reconstruction contract';return;endif
    allocate(orbital_values(size(state%occupations)))
    do p=1,size(density)
      orbital_values=(0d0,0d0)
      do i=1,size(row_ids)
        orbital_values=orbital_values+basis_values(int(row_ids(i)),p)*state%coefficients(i,:)
      enddo
      density(p)=sum(state%occupations*abs(orbital_values)**2)
    enddo
    local_ok=all(ieee_is_finite(density)).and.all(density>=0d0)
    call comm_logical_and(local_ok,collective_ok,dc%icomm_tot)
    if(.not.collective_ok)then;message='nonfinite terminal reconstructed density';return;endif
    ok=.true.;message=''
  end subroutine reconstruct_dg_hybrid_terminal_density

  subroutine measure_dg_hybrid_terminal_density_change(comm,weights,input_density,output_density,&
      density_change,ok,message)
    integer,intent(in)::comm
    real(8),intent(in)::weights(:),input_density(:),output_density(:)
    real(8),intent(out)::density_change
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8)::local_values(2),global_values(2)
    integer::ierr,local_bad,global_bad

    density_change=huge(1d0);ok=.false.;message=''
    local_bad=merge(0,1,size(weights)==size(input_density).and.&
      size(output_density)==size(input_density).and.all(ieee_is_finite(weights)).and.&
      all(weights>0d0).and.all(ieee_is_finite(input_density)).and.&
      all(ieee_is_finite(output_density)).and.all(input_density>=0d0).and.all(output_density>=0d0))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid terminal density-change contract';return;endif
    local_values(1)=sum(weights*(output_density-input_density)**2)
    local_values(2)=max(sum(weights*input_density**2),sum(weights*output_density**2))
    call MPI_Allreduce(local_values,global_values,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='terminal density-change reduction failed';return;endif
    density_change=sqrt(max(0d0,global_values(1))/max(tiny(1d0),global_values(2)))
    ok=ieee_is_finite(density_change)
    if(ok)then;message='';else;message='nonfinite terminal density change';endif
  end subroutine measure_dg_hybrid_terminal_density_change

  subroutine evaluate_dg_hybrid_terminal_total_energy(grid_ids,grid_weights,density,solve_local_potential,&
      state,total_energy,ok,message)
    integer(8),intent(in)::grid_ids(:)
    real(8),intent(in)::grid_weights(:),density(:),solve_local_potential(:)
    type(s_dg_hybrid_ground_state),intent(in)::state
    real(8),intent(out)::total_energy
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_dft_energy)::trial_energy
    real(8)::band_energy,local_parts(2),global_parts(2)
    integer::p,gx,gy,gz,ix,iy,iz,ierr

    total_energy=huge(1d0);ok=.false.;message=''
    if(size(grid_ids)/=size(density).or.size(grid_weights)/=size(density).or.&
       size(solve_local_potential)/=size(density).or..not.state%valid.or.&
       .not.allocated(state%occupations).or..not.allocated(state%eigenvalues))then
      message='invalid terminal total-energy contract';return
    endif
    if(size(state%occupations)/=size(state%eigenvalues))then
      message='terminal total-energy occupation/eigenvalue shape mismatch';return
    endif
    if(any(.not.ieee_is_finite(density)).or.any(.not.ieee_is_finite(solve_local_potential)).or.&
       any(.not.ieee_is_finite(grid_weights)))then;message='nonfinite terminal total-energy input';return;endif
    band_energy=sum(state%occupations*state%eigenvalues)
    local_parts=0d0
    local_parts(1)=sum(grid_weights*density*solve_local_potential)
    call dg_dc_update_potential_from_distributed_density(grid_ids,density,ok,message)
    if(.not.ok)return
    do p=1,size(grid_ids)
      gx=int(modulo(grid_ids(p)-1_8,int(dc%lg_tot%num(1),8)))+1
      gy=int(modulo((grid_ids(p)-1_8)/int(dc%lg_tot%num(1),8),int(dc%lg_tot%num(2),8)))+1
      gz=int((grid_ids(p)-1_8)/int(dc%lg_tot%num(1)*dc%lg_tot%num(2),8))+1
      ix=findloc(dc%jxyz_tot(:,1),gx,dim=1);iy=findloc(dc%jxyz_tot(:,2),gy,dim=1)
      iz=findloc(dc%jxyz_tot(:,3),gz,dim=1)
      if(ix<1.or.iy<1.or.iz<1)then;message='terminal total-energy grid mapping failed';ok=.false.;return;endif
      local_parts(2)=local_parts(2)+eexc_tmp(ix,iy,iz)*grid_weights(p)
    enddo
    call MPI_Allreduce(local_parts,global_parts,2,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.any(.not.ieee_is_finite(global_parts)).or..not.ieee_is_finite(band_energy))then
      message='terminal total-energy reduction failed';ok=.false.;return
    endif
    trial_energy%E_tot=0d0;trial_energy%E_kin=band_energy-global_parts(1)
    trial_energy%E_h=0d0;trial_energy%E_xc=global_parts(2);trial_energy%E_ion_ion=0d0
    trial_energy%E_ion_loc=0d0;trial_energy%E_ion_nloc=0d0;trial_energy%E_U=0d0
    trial_energy%E_tot0=0d0;trial_energy%elec_num=sum(state%occupations)
    trial_energy%elec_num_raw=trial_energy%elec_num;trial_energy%pw_weight_raw=0d0
    call calc_Total_Energy_periodic(dc%mg_tot,ewald,dc%system_tot,dc%info_tot,pp,dc%ppg_tot,&
      dc%fg_tot,dc%poisson_tot,.true.,trial_energy)
    total_energy=trial_energy%E_tot;ok=ieee_is_finite(total_energy)
    if(ok)then;message='';else;message='nonfinite terminal total energy';endif
  end subroutine evaluate_dg_hybrid_terminal_total_energy

  subroutine publish_dg_hybrid_divided_v5(row_ids,row_owner,row_generation,grid_ids,grid_weights,&
      grid_fragment,basis_values,basis_gradients,metric_rows,kinetic_rows,nonlocal_rows,local_rows,&
      sipg_rows,hamiltonian_rows,solved_coefficients,solved_eigenvalues,occupied_state,&
      basis_fingerprint,dc_seed_fingerprint,metric_fingerprint,face_fingerprint,&
      continuation_fingerprint,operator_fingerprint,stationarity_defect,metric_defect,&
      projector_defect,electron_defect,certified_rank_receipt,publication_authorized,ok,message,&
      publication_receipt,terminal_refinement)
    integer(8),intent(in)::row_ids(:),grid_ids(:)
    integer,intent(in)::row_owner(:),row_generation(:),grid_fragment(:)
    real(8),intent(in)::grid_weights(:),solved_eigenvalues(:),stationarity_defect,metric_defect,&
      projector_defect,electron_defect
    complex(8),intent(in)::basis_values(:,:),basis_gradients(:,:,:),metric_rows(:,:),kinetic_rows(:,:),&
      nonlocal_rows(:,:),local_rows(:,:),sipg_rows(:,:),hamiltonian_rows(:,:),solved_coefficients(:,:)
    type(s_dg_hybrid_ground_state),intent(in)::occupied_state
    integer(8),intent(in)::basis_fingerprint,dc_seed_fingerprint,metric_fingerprint,face_fingerprint,&
      continuation_fingerprint,operator_fingerprint
    integer,intent(in)::certified_rank_receipt
    logical,intent(in)::publication_authorized
    type(s_dg_hybrid_candidate_acceptance),intent(in),optional::publication_receipt
    type(s_dg_hybrid_terminal_refinement_receipt),intent(in),optional::terminal_refinement
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_rt_dg_hybrid_v5_shard)::payload
    type(s_rt_dg_hybrid_refinement_receipt)::companion_receipt
    integer,allocatable::metric_offsets(:),metric_columns(:),operator_offsets(:),operator_columns(:)
    complex(8),allocatable::empty_position(:,:,:),orbital_values(:)
    complex(8),allocatable::energy_action(:,:),energy_component_values(:)
    real(8),allocatable::density(:),coordinates(:,:),coordinate_component(:)
    integer::n,nocc,nrow,npoint,rank,nproc,ierr,i,j,p,a,edge,nnz_basis,slot,requested_rank,certified_rank,&
      global_projector_count,energy_state,energy_gx,energy_gy,energy_gz,energy_ix,energy_iy,energy_iz,&
      energy_payload_count
    integer,allocatable::energy_column_slots(:)
    integer(8)::structure_fingerprint
    integer(8)::energy_workspace_peak
    real(8)::electron_count,window,cutoff,cluster_tolerance,local_energy_parts(3),global_energy_parts(3)
    logical::precondition_ok,local_ok,energy_exchange_ok
    character(512)::local_message,energy_exchange_message
    type(s_rt_dg_sparse_exchange)::energy_exchange
    type(s_rt_dg_hybrid_v5_publication_authorization)::authorization
    type(s_dft_energy)::checkpoint_energy
    type(s_dft_system)::identity_system
    type(s_dg_hybrid_candidate_acceptance)::rank_policy_receipt

    ok=.false.;message='';n=size(solved_eigenvalues);nocc=occupied_state%noccupied
    nrow=size(row_ids);npoint=size(grid_ids)
    call MPI_Comm_rank(dc%icomm_tot,rank,ierr);call MPI_Comm_size(dc%icomm_tot,nproc,ierr)
    precondition_ok=ierr==MPI_SUCCESS.and.n>=2.and.nocc>=1.and.nocc<n.and.nproc==dc%n_frag.and.dc%i_frag==rank+1
    precondition_ok=precondition_ok.and.size(row_owner)==n.and.size(row_generation)==n.and.&
      size(grid_weights)==npoint.and.size(grid_fragment)==npoint.and.all(grid_fragment==dc%i_frag)
    precondition_ok=precondition_ok.and.all(shape(basis_values)==[n,npoint]).and.&
      all(shape(metric_rows)==[nrow,n]).and.all(shape(kinetic_rows)==[nrow,n]).and.&
      all(shape(nonlocal_rows)==[nrow,n]).and.all(shape(local_rows)==[nrow,n]).and.&
      all(shape(sipg_rows)==[nrow,n]).and.all(shape(hamiltonian_rows)==[nrow,n])
    precondition_ok=precondition_ok.and.occupied_state%valid.and.occupied_state%converged.and.&
      occupied_state%global_count==n.and.occupied_state%final_eigensolve_count>=1.and.&
      occupied_state%final_eigensolve_count<=4.and.&
      allocated(occupied_state%owned_row_ids).and.allocated(occupied_state%coefficients).and.&
      allocated(occupied_state%occupations).and.allocated(occupied_state%eigenvalues)
    if(precondition_ok)precondition_ok=size(occupied_state%owned_row_ids)==nrow.and.&
      all(shape(occupied_state%coefficients)==[nrow,nocc]).and.size(occupied_state%occupations)==nocc.and.&
      size(occupied_state%eigenvalues)==nocc
    call collective_rt_dg_hybrid_publication_precondition(dc%icomm_tot,precondition_ok,n,nocc,local_ok,local_message)
    if(.not.local_ok)then;message='terminal divided v5 publication precondition failed: '//trim(local_message);return;endif
    call collective_rt_dg_hybrid_publication_mapping_precondition(dc%icomm_tot,n,row_ids,row_owner,&
      occupied_state%owned_row_ids,precondition_ok,local_ok,local_message)
    if(.not.local_ok)then;message='terminal divided v5 row mapping failed: '//trim(local_message);return;endif

    window=max(0d0,dg_hybrid_symmetry_energy_window);cutoff=solved_eigenvalues(nocc)+window;requested_rank=nocc
    do while(requested_rank<n.and.solved_eigenvalues(requested_rank+1)<=cutoff);requested_rank=requested_rank+1;enddo
    certified_rank=requested_rank
    do while(certified_rank<n)
      cluster_tolerance=max(dg_ow_symmetry_tolerance,64d0*epsilon(1d0))*&
        max(1d0,abs(solved_eigenvalues(certified_rank)),abs(solved_eigenvalues(certified_rank+1)))
      if(solved_eigenvalues(certified_rank+1)-solved_eigenvalues(certified_rank)>cluster_tolerance)exit
      certified_rank=certified_rank+1
    enddo
    if(certified_rank_receipt>0)then
      if(certified_rank_receipt<requested_rank.or.certified_rank_receipt>n)then
        message='terminal divided v5 certified-rank receipt is inconsistent';return
      endif
      certified_rank=certified_rank_receipt
    endif
    rank_policy_receipt=s_dg_hybrid_candidate_acceptance()
    if(present(publication_receipt))rank_policy_receipt=publication_receipt
    call validate_dg_hybrid_v5_publication_rank_policy(dc%icomm_tot,rank_policy_receipt,&
      dg_hybrid_symmetry_energy_window,requested_rank,certified_rank,n,local_ok,local_message)
    if(.not.local_ok)then;message=trim(local_message);return;endif

    ! Pointwise support supplies every potentially nonzero local-potential and
    ! position edge.  Fixed matrices add their exact structural support.  No
    ! spectral rotation and no global R-by-R work array is formed here.
    allocate(empty_position(3,0,0))
    call build_rt_dg_hybrid_structural_graph(dc%icomm_tot,n,row_ids,basis_values,metric_rows,kinetic_rows,&
      nonlocal_rows,local_rows,sipg_rows,hamiltonian_rows,empty_position,metric_offsets,metric_columns,&
      operator_offsets,operator_columns,local_ok,local_message,basis_owners=row_owner,local_owner=rank)
    deallocate(empty_position)
    if(.not.local_ok)then;message='terminal divided v5 structural graph failed: '//trim(local_message);return;endif
    call fingerprint_rt_dg_hybrid_sparse_structure(dc%icomm_tot,n,row_ids,operator_offsets,operator_columns,&
      basis_fingerprint,int(z'43454C4C57524150',8),structure_fingerprint,local_ok,local_message)
    if(.not.local_ok)then;message='terminal divided v5 structure fingerprint failed';return;endif

    allocate(payload%row_ids,source=row_ids);allocate(payload%metric_offsets,source=metric_offsets)
    allocate(payload%metric_columns,source=metric_columns);allocate(payload%metric_values(size(metric_columns)))
    allocate(payload%operator_offsets,source=operator_offsets);allocate(payload%operator_columns,source=operator_columns)
    allocate(payload%operator_values(size(operator_columns)),payload%kinetic_values(size(operator_columns)),&
      payload%nonlocal_values(size(operator_columns)),payload%local_values(size(operator_columns)),&
      payload%sipg_values(size(operator_columns)),payload%position_values(3,size(operator_columns)))
    do i=1,nrow
      do edge=metric_offsets(i),metric_offsets(i+1)-1
        j=metric_columns(edge);payload%metric_values(edge)=metric_rows(i,j)
      enddo
      do edge=operator_offsets(i),operator_offsets(i+1)-1
        j=operator_columns(edge);payload%operator_values(edge)=hamiltonian_rows(i,j)
        payload%kinetic_values(edge)=kinetic_rows(i,j);payload%nonlocal_values(edge)=nonlocal_rows(i,j)
        payload%local_values(edge)=local_rows(i,j);payload%sipg_values(edge)=sipg_rows(i,j)
      enddo
    enddo
    nnz_basis=count((basis_values/=(0d0,0d0)).and.spread(row_owner==rank,2,npoint))
    allocate(payload%grid_ids,source=grid_ids);allocate(payload%grid_weights,source=grid_weights)
    allocate(payload%basis_point_offsets(npoint+1),payload%basis_support_ids(nnz_basis),&
      payload%basis_support_values(nnz_basis))
    slot=0;payload%basis_point_offsets(1)=1
    do p=1,npoint
      do j=1,n
        if(row_owner(j)==rank.and.basis_values(j,p)/=(0d0,0d0))then
          slot=slot+1;payload%basis_support_ids(slot)=j;payload%basis_support_values(slot)=basis_values(j,p)
        endif
      enddo
      payload%basis_point_offsets(p+1)=slot+1
    enddo
    allocate(coordinates(3,npoint),coordinate_component(npoint))
    do p=1,npoint
      coordinates(1,p)=real(modulo(grid_ids(p)-1_8,int(dc%lg_tot%num(1),8)),8)*dc%system_tot%hgs(1)
      coordinates(2,p)=real(modulo((grid_ids(p)-1_8)/int(dc%lg_tot%num(1),8),int(dc%lg_tot%num(2),8)),8)*&
        dc%system_tot%hgs(2)
      coordinates(3,p)=real((grid_ids(p)-1_8)/int(dc%lg_tot%num(1)*dc%lg_tot%num(2),8),8)*dc%system_tot%hgs(3)
    enddo
    do a=1,3
      coordinate_component=coordinates(a,:)
      call project_rt_dg_hybrid_point_csr_edges(dc%icomm_tot,n,row_ids,operator_offsets,operator_columns,&
        grid_ids,grid_weights,payload%basis_point_offsets,payload%basis_support_ids,&
        payload%basis_support_values,coordinate_component,payload%position_values(a,:),local_ok,local_message)
      if(.not.local_ok)then;message='terminal divided v5 sparse position projection failed: '//trim(local_message);return;endif
    enddo
    allocate(density(npoint),orbital_values(nocc));density=0d0
    do p=1,npoint
      orbital_values=(0d0,0d0)
      do i=1,nrow
        orbital_values=orbital_values+basis_values(int(row_ids(i)),p)*occupied_state%coefficients(i,:)
      enddo
      density(p)=sum(occupied_state%occupations*abs(orbital_values)**2)
    enddo
    call MPI_Allreduce(sum(density*grid_weights),electron_count,1,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)then;message='terminal divided v5 density electron reduction failed';return;endif
    ! Synchronize only the Hartree/XC potential used for the physical energy
    ! receipt with the immutable LCFO density being published.  This does not
    ! update the density, diagonalize again, or alter the localized basis and
    ! coefficients; RT performs the identical t=0 refresh.
    call dg_dc_update_potential_from_distributed_density(grid_ids,density,local_ok,local_message)
    if(.not.local_ok)then;message='terminal divided v5 energy potential refresh failed: '//trim(local_message);return;endif
    allocate(payload%density,source=density)
    allocate(payload%initial_occupied_amplitudes,source=occupied_state%coefficients)
    allocate(payload%occupations,source=occupied_state%occupations)
    allocate(payload%eigenvalues,source=occupied_state%eigenvalues)
    allocate(payload%scope_selectors(8),payload%xc_types(size(xc_func%xctype)))
    payload%scope_selectors=[1,1,1,0,0,0,0,0];payload%xc_types=xc_func%xctype
    allocate(payload%acceptance_receipts(8));payload%acceptance_receipts=[stationarity_defect,metric_defect,&
      projector_defect,electron_defect,electron_count,dc%elec_num_tot,real(requested_rank,8),real(certified_rank,8)]
    global_projector_count=dc%ppg_tot%Nlma
    allocate(payload%pseudopotential_receipt(6));payload%pseudopotential_receipt=[real(dc%system_tot%nion,8),&
      canonical_pp_valence_sum(pp),real(pp%lmax,8),real(pp%nrmax,8),real(global_projector_count,8),real(n,8)]
    ! Evaluate the same physical energy decomposition consumed by RT without
    ! gathering the distributed occupied coefficient rows.  The sparse halo
    ! is bounded by the frozen operator graph and is released with this scope.
    call build_rt_dg_sparse_exchange(dc%icomm_tot,n,structure_fingerprint,row_ids,&
      operator_columns,energy_exchange,local_ok,local_message)
    if(.not.local_ok)then;message='terminal divided v5 energy halo construction failed: '//trim(local_message);return;endif
    allocate(energy_action(nrow,nocc),energy_component_values(size(operator_columns)),&
      energy_column_slots(size(operator_columns)))
    energy_column_slots=[(i,i=1,size(energy_column_slots))]
    energy_component_values=payload%kinetic_values+payload%sipg_values
    call apply_rt_dg_sparse_rows_tiled(dc%icomm_tot,energy_exchange,operator_offsets,&
      energy_component_values,energy_column_slots,occupied_state%coefficients,energy_action,16,&
      energy_workspace_peak,energy_payload_count,energy_exchange_ok,energy_exchange_message)
    if(.not.energy_exchange_ok)then
      message='terminal divided v5 energy coefficient exchange failed: '//trim(energy_exchange_message);return
    endif
    local_energy_parts=0d0
    do energy_state=1,nocc
      local_energy_parts(1)=local_energy_parts(1)+occupied_state%occupations(energy_state)*real(sum(&
        conjg(occupied_state%coefficients(:,energy_state))*energy_action(:,energy_state)),8)
    enddo
    call apply_rt_dg_sparse_rows_tiled(dc%icomm_tot,energy_exchange,operator_offsets,&
      payload%nonlocal_values,energy_column_slots,occupied_state%coefficients,energy_action,16,&
      energy_workspace_peak,energy_payload_count,energy_exchange_ok,energy_exchange_message)
    if(.not.energy_exchange_ok)then
      message='terminal divided v5 nonlocal energy action failed: '//trim(energy_exchange_message);return
    endif
    do energy_state=1,nocc
      local_energy_parts(2)=local_energy_parts(2)+occupied_state%occupations(energy_state)*real(sum(&
        conjg(occupied_state%coefficients(:,energy_state))*energy_action(:,energy_state)),8)
    enddo
    do p=1,npoint
      energy_gx=int(modulo(grid_ids(p)-1_8,int(dc%lg_tot%num(1),8)))+1
      energy_gy=int(modulo((grid_ids(p)-1_8)/int(dc%lg_tot%num(1),8),int(dc%lg_tot%num(2),8)))+1
      energy_gz=int((grid_ids(p)-1_8)/int(dc%lg_tot%num(1)*dc%lg_tot%num(2),8))+1
      energy_ix=findloc(dc%jxyz_tot(:,1),energy_gx,dim=1)
      energy_iy=findloc(dc%jxyz_tot(:,2),energy_gy,dim=1)
      energy_iz=findloc(dc%jxyz_tot(:,3),energy_gz,dim=1)
      if(energy_ix<1.or.energy_iy<1.or.energy_iz<1)then
        message='terminal divided v5 energy grid mapping failed';return
      endif
      local_energy_parts(3)=local_energy_parts(3)+eexc_tmp(energy_ix,energy_iy,energy_iz)*grid_weights(p)
    enddo
    call MPI_Allreduce(local_energy_parts,global_energy_parts,3,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.any(.not.ieee_is_finite(global_energy_parts)))then
      message='terminal divided v5 energy decomposition reduction failed';return
    endif
    checkpoint_energy%E_kin=global_energy_parts(1);checkpoint_energy%E_ion_nloc=global_energy_parts(2)
    checkpoint_energy%E_xc=global_energy_parts(3)
    call calc_Total_Energy_periodic(dc%mg_tot,ewald,dc%system_tot,dc%info_tot,pp,dc%ppg_tot,&
      dc%fg_tot,dc%poisson_tot,.true.,checkpoint_energy)
    allocate(payload%energy_receipt(7))
    payload%energy_receipt=[checkpoint_energy%E_tot,checkpoint_energy%E_kin,checkpoint_energy%E_h,&
      checkpoint_energy%E_xc,checkpoint_energy%E_ion_ion,checkpoint_energy%E_ion_loc,checkpoint_energy%E_ion_nloc]
    if(any(.not.ieee_is_finite(payload%energy_receipt)).or.&
        abs(payload%energy_receipt(1)-sum(payload%energy_receipt(2:7)))>&
        100d0*epsilon(1d0)*max(1d0,abs(payload%energy_receipt(1))))then
      message='terminal divided v5 final energy receipt is inconsistent';return
    endif
    if(rank==0)write(*,'(a,7(a,es16.8))')'[HYBRID-GS-ENERGY-RECEIPT]',&
      ' total=',payload%energy_receipt(1),' kinetic=',payload%energy_receipt(2),&
      ' hartree=',payload%energy_receipt(3),' xc=',payload%energy_receipt(4),&
      ' ion_ion=',payload%energy_receipt(5),' ion_local=',payload%energy_receipt(6),&
      ' ion_nonlocal=',payload%energy_receipt(7)
    payload%global_count=n;payload%global_grid_count=product(dc%lg_tot%num);payload%nocc=nocc
    payload%certified_rank=certified_rank;payload%fragment_id=rank+1
    payload%basis_fingerprint=basis_fingerprint;payload%operator_fingerprint=operator_fingerprint
    payload%operator_structure_fingerprint=structure_fingerprint
    payload%scope_fingerprint=fingerprint_rt_dg_hybrid_scope(payload%scope_selectors,payload%xc_types)
    ! Physical system identity is independent of solver state count and thermal
    ! redistribution; LCFO occupations remain authenticated in the payload.
    identity_system=dc%system_tot
    if(identity_system%nspin/=1.or.identity_system%nk/=1)then
      message='terminal divided v5 electronic identity shape is unsupported';return
    endif
    payload%pseudopotential_digest=canonical_pp_digest(pp)
    payload%system_fingerprint=fingerprint_rt_dg_hybrid_system(identity_system,dc%lg_tot%num,.true.,&
      dc%ppg_tot%Nlma,payload%pseudopotential_digest,xc_func%xctype,.false.,.false.,.false.,.false.,.false.,&
      nint(sum(occupied_state%occupations)),[0,0])
    payload%pseudopotential_fingerprint=canonical_pp_fingerprint(pp)
    payload%payload_fingerprint=ieor(ieor(basis_fingerprint,operator_fingerprint),&
      ieor(dc_seed_fingerprint,ieor(face_fingerprint,continuation_fingerprint)))
    if(payload%payload_fingerprint==0_8)payload%payload_fingerprint=1_8
    authorization%valid=publication_authorized;authorization%checkpoint_version=5
    authorization%published_rank=n;authorization%basis_fingerprint=basis_fingerprint
    authorization%operator_fingerprint=operator_fingerprint
    call publish_rt_dg_hybrid_checkpoint_v5(dc%icomm_tot,'./hybrid_dg_ground_state.chk',n,nocc,&
      row_ids,row_owner,occupied_state%owned_row_ids,payload,authorization,precondition_ok,local_ok,local_message)
    if(.not.local_ok)then;message='terminal divided v5 write failed: '//trim(local_message);return;endif
    if(present(terminal_refinement))then
      companion_receipt%version=1;companion_receipt%fragment_id=rank+1
      companion_receipt%v5_publication_fingerprint=payload%payload_fingerprint
      companion_receipt%total_solve_count=terminal_refinement%total_solve_count
      companion_receipt%additional_refinement_count=terminal_refinement%additional_refinement_count
      companion_receipt%density_change=terminal_refinement%density_change
      companion_receipt%energy_change=terminal_refinement%energy_change
      companion_receipt%converged=terminal_refinement%converged
      companion_receipt%exhausted=terminal_refinement%exhausted
      if(terminal_refinement%exhausted)then
        companion_receipt%exit_reason='maximum-additional-solves-exhausted'
      else if(terminal_refinement%total_solve_count==1)then
        companion_receipt%exit_reason='initial-lcfo-converged'
      else
        companion_receipt%exit_reason='refined-lcfo-converged'
      endif
      call write_rt_dg_hybrid_refinement_receipt(dc%icomm_tot,&
        './hybrid_dg_ground_state.chk.refinement',companion_receipt,local_ok,local_message)
      if(.not.local_ok)then;message='terminal refinement receipt write failed: '//trim(local_message);return;endif
      if(rank==0.and.terminal_refinement%exhausted)write(error_unit,'(a)')&
        '[DG-HYBRID-REFINEMENT-WARNING] maximum additional LCFO solves exhausted; publishing last finite valid state'
    endif
    if(rank==0)write(*,'(a,6(a,i0),4(a,es16.8),a,i0)')'[HYBRID-GS-HANDOFF] route=divided-terminal-lcfo-v5',&
      ' construction_rank=',n,' solved_rank=',n,' certified_rank=',certified_rank,' rt_rank=',n,&
      ' occupied_rank=',nocc,' projector_count=',global_projector_count,' stationarity=',stationarity_defect,&
      ' metric=',metric_defect,' projector=',projector_defect,' electron=',electron_defect,' writer_count=',1
    ok=.true.;message=''
  end subroutine publish_dg_hybrid_divided_v5


  subroutine solve_dg_hybrid_schwarz_fragments(iteration,callback_ok)
    integer,intent(in)::iteration
    logical,intent(out)::callback_ok
    real(8),allocatable::core_potential(:)
    real(8),allocatable::occupations(:),energies(:)
    real(8)::diagnostics(2)
    integer(8)::potential_fingerprint
    integer::local_iterations,extensions,rank_local,ierr_local,status_local
    logical::converged,rolled_back,extended,collective_callback_ok,all_extended,all_not_extended
    character(512)::solver_message

    bounded_last_accepted_cg_steps=0
    call MPI_Comm_rank(dc%icomm_tot,rank_local,ierr_local)
    callback_ok=ierr_local==MPI_SUCCESS.and.bounded_schwarz_state%valid.and.&
      bounded_schwarz_schedule%valid.and.allocated(bounded_core_ids)
    call comm_logical_and(callback_ok,collective_callback_ok,dc%icomm_tot)
    if(.not.collective_callback_ok)then
      callback_ok=.false.
      if(ierr_local==MPI_SUCCESS.and.rank_local==0)&
        write(error_unit,'(a)')'Schwarz callback preflight failed collectively'
      return
    endif
    allocate(core_potential(size(bounded_core_ids)),stat=status_local)
    callback_ok=status_local==0
    call comm_logical_and(callback_ok,collective_callback_ok,dc%icomm_tot)
    if(.not.collective_callback_ok)then
      callback_ok=.false.
      if(rank_local==0)write(error_unit,'(a)')'Schwarz core-potential allocation failed collectively'
      return
    endif
    call extract_dg_hybrid_core_local_potential(bounded_core_ids,core_potential,callback_ok)
    call comm_logical_and(callback_ok,collective_callback_ok,dc%icomm_tot)
    if(.not.collective_callback_ok)then
      callback_ok=.false.
      if(rank_local==0)write(error_unit,'(a)')'Schwarz core-potential extraction failed collectively'
      return
    endif
    call assemble_dg_hybrid_local_potential_rows(dc%icomm_tot,bounded_fixed_payload%global_basis_count,&
      divided_fragment_basis%global_ids,bounded_basis_fragment,bounded_core_ids,&
      bounded_interior_fragment,bounded_core_weights,bounded_interior_values,core_potential,&
      bounded_local_potential_rows,diagnostics,callback_ok,solver_message)
    call comm_logical_and(callback_ok,collective_callback_ok,dc%icomm_tot)
    if(.not.collective_callback_ok)then
      callback_ok=.false.
      if(rank_local==0)write(error_unit,'(a,a)')&
        'bounded fragment potential projection failed collectively: ',trim(solver_message)
      return
    endif
    call ow_fingerprint_distributed_matrix(dc%icomm_tot,divided_fragment_basis%global_ids,&
      bounded_local_potential_rows,potential_fingerprint,callback_ok)
    call comm_logical_and(callback_ok,collective_callback_ok,dc%icomm_tot)
    if(.not.collective_callback_ok)then
      callback_ok=.false.
      if(rank_local==0)write(error_unit,'(a)')'Schwarz potential fingerprint failed collectively'
      return
    endif
    call assemble_dg_hybrid_schwarz_local_preconditioner_blocks(callback_ok)
    call comm_logical_and(callback_ok,collective_callback_ok,dc%icomm_tot)
    if(.not.collective_callback_ok)then
      callback_ok=.false.
      if(rank_local==0)write(error_unit,'(a)')&
        'Schwarz local preconditioner block assembly failed collectively'
      return
    endif
    bounded_last_peer_exchange_count=0
    call advance_dg_hybrid_schwarz_epoch(dc%icomm_tot,bounded_schwarz_state%basis_generation,&
      dg_hybrid_fragment_cg_steps,dg_dc_gs_intermediate_orbital_tolerance,&
      dg_dc_gs_orthogonality_tolerance,dg_dc_gs_allowed_residual_growth,&
      apply_dg_hybrid_schwarz_h,apply_dg_hybrid_schwarz_s,&
      apply_dg_hybrid_schwarz_preconditioner,bounded_schwarz_state,local_iterations,&
      divided_fragment_residual,divided_fragment_orthogonality,converged,rolled_back,&
      callback_ok,solver_message)
    bounded_last_accepted_cg_steps=local_iterations
    callback_ok=callback_ok.and..not.rolled_back
    call comm_logical_and(callback_ok,collective_callback_ok,dc%icomm_tot)
    if(.not.collective_callback_ok)then
      callback_ok=.false.
      if(rank_local==0)write(error_unit,'(a,a)')&
        'bounded Schwarz update failed or rolled back collectively: ',trim(solver_message)
      return
    endif
    callback_ok=local_iterations<=dg_hybrid_fragment_cg_steps
    call comm_logical_and(callback_ok,collective_callback_ok,dc%icomm_tot)
    if(.not.collective_callback_ok)then
      callback_ok=.false.
      if(rank_local==0)write(error_unit,'(a)')'Schwarz CG step cap was exceeded collectively'
      return
    endif
    extensions=0
    do
      call assign_dg_hybrid_schwarz_occupations(dc%icomm_tot,bounded_schwarz_state%basis_generation,&
        300d0,2d0,dc%elec_num_tot,1d-10,dg_dc_gs_subspace_tolerance,&
        bounded_schwarz_candidate_ids,bounded_schwarz_candidate_energies,&
        bounded_schwarz_candidate_vectors,apply_dg_hybrid_schwarz_h,apply_dg_hybrid_schwarz_s,&
        bounded_schwarz_state,occupations,energies,extended,callback_ok,solver_message)
      call comm_logical_and(callback_ok,collective_callback_ok,dc%icomm_tot)
      if(.not.collective_callback_ok)then
        callback_ok=.false.
        if(rank_local==0)write(error_unit,'(a,a)')&
          'global Schwarz occupations failed collectively: ',trim(solver_message)
        return
      endif
      call comm_logical_and(extended,all_extended,dc%icomm_tot)
      call comm_logical_and(.not.extended,all_not_extended,dc%icomm_tot)
      if(.not.all_extended.and..not.all_not_extended)then
        callback_ok=.false.
        if(rank_local==0)write(error_unit,'(a)')&
          'global Schwarz occupation extension decision disagrees across ranks'
        return
      endif
      if(all_not_extended)exit
      extensions=extensions+1
    enddo
    system%mu=bounded_schwarz_state%chemical_potential
    if(rank_local==0)write(*,'(a,i0,3(a,i0),3(a,es12.4))')'[DG-HYBRID-SCHWARZ] epoch=',iteration,&
      ' neighbor_exchanges=',bounded_last_peer_exchange_count,' accepted_cg_steps=',local_iterations,&
      ' common_extensions=',extensions,' residual=',divided_fragment_residual,&
      ' electron_defect=',bounded_schwarz_state%electron_defect,&
      ' temperature=',bounded_schwarz_state%temperature
    callback_ok=local_iterations<=dg_hybrid_fragment_cg_steps
  end subroutine solve_dg_hybrid_schwarz_fragments

  subroutine record_dg_hybrid_interface_continuation_diagnostic(&
      diagnostic_state_lambda,continuation_fingerprint,continuation_state,local_accept,diagnostic_ok)
    real(8),intent(in)::diagnostic_state_lambda
    integer(int64),intent(in)::continuation_fingerprint
    type(s_dg_hybrid_interface_continuation),intent(inout)::continuation_state
    logical,intent(in)::local_accept
    logical,intent(out)::diagnostic_ok
    complex(8),allocatable::hcoeff(:,:),scoeff(:,:),scaled_interface_action(:,:)
    real(8),allocatable::local_h_diagonal(:),global_h_diagonal(:),&
      local_s_diagonal(:),global_s_diagonal(:)
    real(8)::local_interface_norm_squared,global_interface_norm_squared,&
      rayleigh_energy_trace,scaled_interface_action_norm,residual_record,&
      orthogonality_record,electron_defect_record,diagnostic_state_lambda_record
    integer(int64)::continuation_fingerprint_min,continuation_fingerprint_max
    integer(int64)::dynamic_receipt
    integer::nstate,j,status_local,ierr_local,rank_local,peer_exchanges,minimum_steps,maximum_steps,&
      minimum_accept_request,maximum_accept_request
    real(8)::minimum_state_lambda,maximum_state_lambda
    logical::local_ok,collective_diagnostic_ok,measurement_available,accept_ok,record_ok,dynamic_ok
    character(16)::measurement_status,record_status
    character(512)::diagnostic_message

    call MPI_Comm_rank(dc%icomm_tot,rank_local,ierr_local)
    local_ok=ierr_local==MPI_SUCCESS.and.bounded_schwarz_state%valid.and.&
      allocated(bounded_schwarz_state%coefficients).and.allocated(bounded_schwarz_state%occupations)
    call comm_logical_and(local_ok,collective_diagnostic_ok,dc%icomm_tot)
    measurement_available=collective_diagnostic_ok
    nstate=bounded_schwarz_state%trial_count
    if(measurement_available)then
      allocate(hcoeff(bounded_schwarz_state%local_basis_count,nstate),&
        scoeff(bounded_schwarz_state%local_basis_count,nstate),&
        scaled_interface_action(bounded_schwarz_state%local_basis_count,nstate),&
        local_h_diagonal(nstate),global_h_diagonal(nstate),local_s_diagonal(nstate),&
        global_s_diagonal(nstate),stat=status_local)
      call comm_logical_and(status_local==0,collective_diagnostic_ok,dc%icomm_tot)
      measurement_available=collective_diagnostic_ok
    endif
    if(measurement_available)then
      call apply_dg_hybrid_schwarz_hamiltonian(dc%icomm_tot,bounded_schwarz_schedule,&
        bounded_schwarz_state%basis_generation,bounded_directory_fingerprint,bounded_face_fingerprint,&
        bounded_mapping_fingerprint,divided_fragment_basis%global_ids,bounded_basis_fragment,&
        bounded_basis_local_slot,bounded_fixed_payload%kinetic_rows,bounded_fixed_payload%nonlocal_rows,&
        bounded_fixed_payload%interface_rows,bounded_local_potential_rows,bounded_interface_scale,&
        bounded_schwarz_state%coefficients,hcoeff,peer_exchanges,local_ok,diagnostic_message)
      call comm_logical_and(local_ok,collective_diagnostic_ok,dc%icomm_tot)
      measurement_available=collective_diagnostic_ok
    endif
    if(measurement_available)then
      call apply_dg_hybrid_schwarz_s(bounded_schwarz_state%coefficients,scoeff,local_ok)
      call comm_logical_and(local_ok,collective_diagnostic_ok,dc%icomm_tot)
      measurement_available=collective_diagnostic_ok
    endif
    if(measurement_available)then
      call apply_dg_hybrid_schwarz_rows(dc%icomm_tot,bounded_schwarz_schedule,&
        bounded_schwarz_state%basis_generation,bounded_directory_fingerprint,bounded_face_fingerprint,&
        bounded_mapping_fingerprint,divided_fragment_basis%global_ids,bounded_basis_fragment,&
        bounded_basis_local_slot,bounded_interface_scale*bounded_fixed_payload%interface_rows,&
        bounded_schwarz_state%coefficients,scaled_interface_action,peer_exchanges,local_ok,diagnostic_message)
      call comm_logical_and(local_ok,collective_diagnostic_ok,dc%icomm_tot)
      measurement_available=collective_diagnostic_ok
    endif
    if(measurement_available)then
      do j=1,nstate
        local_h_diagonal(j)=real(sum(conjg(bounded_schwarz_state%coefficients(:,j))*hcoeff(:,j)),8)
        local_s_diagonal(j)=real(sum(conjg(bounded_schwarz_state%coefficients(:,j))*scoeff(:,j)),8)
      enddo
      call MPI_Allreduce(local_h_diagonal,global_h_diagonal,nstate,MPI_DOUBLE_PRECISION,MPI_SUM,&
        dc%icomm_tot,ierr_local)
      measurement_available=ierr_local==MPI_SUCCESS
      call MPI_Allreduce(local_s_diagonal,global_s_diagonal,nstate,MPI_DOUBLE_PRECISION,MPI_SUM,&
        dc%icomm_tot,ierr_local)
      measurement_available=measurement_available.and.ierr_local==MPI_SUCCESS
      local_interface_norm_squared=sum(abs(scaled_interface_action)**2)
      call MPI_Allreduce(local_interface_norm_squared,global_interface_norm_squared,1,&
        MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr_local)
      measurement_available=measurement_available.and.ierr_local==MPI_SUCCESS.and.&
        all(global_s_diagonal>dg_dc_metric_rank_tolerance)
    endif
    if(measurement_available)then
      ! The trace is occupation weighted and normalizes each diagonal Rayleigh
      ! quotient by its live distributed S norm.  The separate orthogonality
      ! defect reports the off-diagonal S overlap that can still bias this trace.
      rayleigh_energy_trace=bounded_schwarz_state%wspin*sum(&
        bounded_schwarz_state%occupations*(global_h_diagonal/global_s_diagonal))
      ! Each rank owns every coefficient-space output row exactly once.  Summing
      ! those row-local squares gives the row-owned Frobenius action norm
      ! ||lambda H_interface C||_F without claiming basis invariance.
      scaled_interface_action_norm=sqrt(max(0d0,global_interface_norm_squared))
      measurement_available=ieee_is_finite(rayleigh_energy_trace).and.&
        ieee_is_finite(scaled_interface_action_norm)
    endif
    if(measurement_available)then
      call validate_dg_hybrid_schwarz_dynamic_receipt(dc%icomm_tot,bounded_schwarz_state,&
        dynamic_receipt,dynamic_ok,diagnostic_message)
      measurement_available=dynamic_ok
    endif
    ! Certify every field used by the rank-zero record before allowing the
    ! continuation controller to advance.  A metadata failure is therefore a
    ! measured rollback, not an acceptance whose controller state must be undone.
    record_ok=ieee_is_finite(diagnostic_state_lambda).and.diagnostic_state_lambda>=0d0.and.&
      diagnostic_state_lambda<=1d0
    call MPI_Allreduce(bounded_last_accepted_cg_steps,minimum_steps,1,MPI_INTEGER,MPI_MIN,&
      dc%icomm_tot,ierr_local)
    record_ok=record_ok.and.ierr_local==MPI_SUCCESS
    call MPI_Allreduce(bounded_last_accepted_cg_steps,maximum_steps,1,MPI_INTEGER,MPI_MAX,&
      dc%icomm_tot,ierr_local)
    record_ok=record_ok.and.ierr_local==MPI_SUCCESS.and.minimum_steps==maximum_steps.and.&
      minimum_steps>=0.and.maximum_steps<=dg_hybrid_fragment_cg_steps
    call MPI_Allreduce(continuation_fingerprint,continuation_fingerprint_min,1,MPI_INTEGER8,MPI_MIN,&
      dc%icomm_tot,ierr_local)
    record_ok=record_ok.and.ierr_local==MPI_SUCCESS
    call MPI_Allreduce(continuation_fingerprint,continuation_fingerprint_max,1,MPI_INTEGER8,MPI_MAX,&
      dc%icomm_tot,ierr_local)
    record_ok=record_ok.and.ierr_local==MPI_SUCCESS.and.&
      continuation_fingerprint_min==continuation_fingerprint_max.and.continuation_fingerprint_min/=0_int64.and.&
      continuation_fingerprint==continuation_state%fingerprint
    call MPI_Allreduce(merge(1,0,local_accept),minimum_accept_request,1,MPI_INTEGER,MPI_MIN,&
      dc%icomm_tot,ierr_local)
    record_ok=record_ok.and.ierr_local==MPI_SUCCESS
    call MPI_Allreduce(merge(1,0,local_accept),maximum_accept_request,1,MPI_INTEGER,MPI_MAX,&
      dc%icomm_tot,ierr_local)
    record_ok=record_ok.and.ierr_local==MPI_SUCCESS.and.minimum_accept_request==maximum_accept_request
    call MPI_Allreduce(diagnostic_state_lambda,minimum_state_lambda,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
      dc%icomm_tot,ierr_local)
    record_ok=record_ok.and.ierr_local==MPI_SUCCESS
    call MPI_Allreduce(diagnostic_state_lambda,maximum_state_lambda,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
      dc%icomm_tot,ierr_local)
    record_ok=record_ok.and.ierr_local==MPI_SUCCESS.and.ieee_is_finite(minimum_state_lambda).and.&
      ieee_is_finite(maximum_state_lambda).and.minimum_state_lambda==maximum_state_lambda.and.&
      minimum_state_lambda>=0d0.and.maximum_state_lambda<=1d0
    call comm_logical_and(record_ok,collective_diagnostic_ok,dc%icomm_tot)
    measurement_available=measurement_available.and.collective_diagnostic_ok
    if(measurement_available)then
      measurement_status='valid'
    else
      measurement_status='unavailable'
      rayleigh_energy_trace=huge(1d0)
      scaled_interface_action_norm=huge(1d0)
    endif
    if(local_accept.and.measurement_available)then
      call accept_dg_hybrid_interface_point(dc%icomm_tot,bounded_schwarz_state%basis_generation,&
        bounded_mapping_fingerprint,.true.,continuation_state,accept_ok,diagnostic_message)
    else
      call accept_dg_hybrid_interface_point(dc%icomm_tot,bounded_schwarz_state%basis_generation,&
        bounded_mapping_fingerprint,.false.,continuation_state,accept_ok,diagnostic_message)
      accept_ok=.false.
    endif
    diagnostic_ok=local_accept.and.measurement_available.and.accept_ok
    record_status=merge('accepted','rollback',diagnostic_ok)
    ! diagnostic_state_lambda identifies the coefficient/occupation snapshot
    ! probed by the attempted operator lambda; on rollback it is not promoted.
    ! Preserve a parseable finite record even when a rejected solver scalar is NaN/Inf.
    residual_record=finite_diagnostic_value(divided_fragment_residual)
    orthogonality_record=finite_diagnostic_value(divided_fragment_orthogonality)
    electron_defect_record=finite_diagnostic_value(bounded_schwarz_state%electron_defect)
    diagnostic_state_lambda_record=finite_diagnostic_value(diagnostic_state_lambda)
    if(rank_local==0)write(*,'(2(a,es24.16),a,i0,5(a,es24.16),5a,i0)')&
      '[DG-HYBRID-CONTINUATION] lambda=',bounded_interface_scale,&
      ' diagnostic_state_lambda=',diagnostic_state_lambda_record,&
      ' accepted_cg_steps=',minimum_steps,' residual=',residual_record,&
      ' orthogonality_defect=',orthogonality_record,&
      ' electron_defect=',electron_defect_record,&
      ' rayleigh_energy_trace=',rayleigh_energy_trace,&
      ' scaled_interface_action_norm=',scaled_interface_action_norm,&
      ' measurement_status=',trim(measurement_status),' status=',trim(record_status),&
      ' continuation_fingerprint=',continuation_fingerprint_min
  end subroutine record_dg_hybrid_interface_continuation_diagnostic

  pure real(8) function finite_diagnostic_value(value)result(record_value)
    real(8),intent(in)::value
    if(ieee_is_finite(value))then
      record_value=value
    else
      record_value=huge(1d0)
    endif
  end function finite_diagnostic_value

  subroutine assemble_dg_hybrid_schwarz_local_preconditioner_blocks(callback_ok)
    logical,intent(out)::callback_ok
    integer::global_column,local_slot,local_count,status_local
    local_count=size(divided_fragment_basis%global_ids);callback_ok=.false.
    if(allocated(bounded_fragment_h))deallocate(bounded_fragment_h)
    if(allocated(bounded_fragment_s))deallocate(bounded_fragment_s)
    allocate(bounded_fragment_h(local_count,local_count),bounded_fragment_s(local_count,local_count),&
      stat=status_local)
    if(status_local/=0)return
    bounded_fragment_h=0d0;bounded_fragment_s=0d0
    do global_column=1,bounded_fixed_payload%global_basis_count
      if(bounded_basis_fragment(global_column)/=dc%i_frag)cycle
      local_slot=bounded_basis_local_slot(global_column)
      bounded_fragment_h(:,local_slot)=bounded_fixed_payload%kinetic_rows(:,global_column)+&
        bounded_fixed_payload%nonlocal_rows(:,global_column)+&
        bounded_interface_scale*bounded_fixed_payload%interface_rows(:,global_column)+&
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
      bounded_fixed_payload%interface_rows,bounded_local_potential_rows,bounded_interface_scale,input,candidate,&
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

















  subroutine update_dg_hybrid_divided_potential(core_density,callback_ok)
    real(8),intent(in)::core_density(:)
    logical,intent(out)::callback_ok
    character(256)::potential_message

    call dg_dc_update_potential_from_distributed_density(ow_core_ids,core_density,callback_ok,potential_message)
    if(.not.callback_ok)write(0,'(2a)')'divided potential update: ',trim(potential_message)
  end subroutine update_dg_hybrid_divided_potential



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
#endif

end subroutine main_dft

#if defined(USE_MPI) && defined(USE_SCALAPACK)
#endif
