#!/usr/bin/env python3
from pathlib import Path
import re
import sys

root = Path(__file__).resolve().parents[2]
source = (root / "src/gs/main_dft.f90").read_text()

flag = "yn_dg_hybrid_continuation_scf == 'y'"
assert flag in source, "missing explicit DG continuation production branch"
branch = source[source.index(flag):]
branch = branch[:branch.index("endif")]
assert "run_dg_hybrid_continuation_ground_state_for_main" in branch

ground_state_name = "subroutine run_dg_overlapping_wannier_ground_state_for_main"
ground_state_match = re.search(
    rf"{ground_state_name}(?P<body>.*?)end\s+subroutine",
    source,
    re.I | re.S,
)
assert ground_state_match, "missing overlapping-Wannier production body"
ground_state = ground_state_match.group("body")
ground_state_lower = ground_state.lower()
branch_begin = "! hybrid_localization_first_branch_begin"
arm_begin = "! hybrid_localization_first_arm_begin"
arm_end = "! hybrid_localization_first_arm_end"
legacy_begin = "! hybrid_constrained_legacy_arm_begin"
legacy_end = "! hybrid_constrained_legacy_arm_end"
branch_end = "! hybrid_localization_first_branch_end"
for marker in (branch_begin, arm_begin, arm_end, legacy_begin, legacy_end, branch_end):
    assert ground_state_lower.count(marker) == 1, f"missing unique route marker {marker}"

raw_seed_position = ground_state_lower.index(
    "global_seed_values(nstate+projector_tile_first:nstate+projector_tile_last,:)="
)
crystal_map_position = ground_state_lower.index("call prepare_ow_global_point_action")
branch_position = ground_state_lower.index(branch_begin)
assert raw_seed_position < crystal_map_position < branch_position, (
    "Hybrid localization must branch once after raw occupied+s+p assembly and crystal-map construction"
)

localization_arm = ground_state_lower[
    ground_state_lower.index(arm_begin) + len(arm_begin):ground_state_lower.index(arm_end)
]
legacy_arm = ground_state_lower[
    ground_state_lower.index(legacy_begin) + len(legacy_begin):ground_state_lower.index(legacy_end)
]
common_tail = ground_state_lower[ground_state_lower.index(branch_end) + len(branch_end):]

for token in (
    "call find_dg_group_identity",
    "call select_dg_group_generators",
    "call measure_dg_grid_map_stencil_defect",
    "call prepare_dg_hybrid_localization_first_seed",
    "call materialize_ow_distributed_core_to_buffer",
    "call reindex_dg_point_maps_between_row_layouts",
    "ow_core_ids=initial_core_ids",
    "call setup_dg_w90_gamma_library",
    "dg_w90_unconstrained",
    "call assemble_dg_w90_gamma_matrices",
    "call run_dg_w90_gamma_library",
    "call apply_dg_w90_gamma_transform",
    "call build_dg_hybrid_localization_receipt",
    "[hybrid-wf-localization] symmetry_constraint=off",
):
    assert token in localization_arm, f"localization-first arm omits {token}"

localization_order = [
    localization_arm.index("call prepare_dg_hybrid_localization_first_seed"),
    localization_arm.index("call materialize_ow_distributed_core_to_buffer"),
    localization_arm.index("call reindex_dg_point_maps_between_row_layouts"),
    localization_arm.index("ow_core_ids=initial_core_ids"),
    localization_arm.index("call setup_dg_w90_gamma_library"),
    localization_arm.index("call assemble_dg_w90_gamma_matrices"),
    localization_arm.index("call run_dg_w90_gamma_library"),
    localization_arm.index("call apply_dg_w90_gamma_transform"),
    localization_arm.index("call build_dg_hybrid_localization_receipt"),
]
assert localization_order == sorted(localization_order), (
    "localization-first route must prepare, reindex, localize, apply, and certify one fixed-rank gauge in order"
)
assert "global_seed_values(1:nstate,:)=adapted" not in localization_arm.replace(" ", ""), (
    "localization-first route must not replace the raw occupied seed block by a symmetry-adapted block"
)
assert "if(.not.center_diagnostic_ok)" not in localization_arm.replace(" ", ""), (
    "individual-WF and center symmetry diagnostics must not become acceptance gates"
)
for forbidden_gate in (
    "global_retained_group_closure_defect>dg_ow_symmetry_tolerance",
    "maxval(lcfo_total_symmetry_residual)>dg_ow_symmetry_tolerance",
):
    assert forbidden_gate not in localization_arm.replace(" ", ""), (
        "complete construction-basis symmetry is diagnostic-only before LCFO"
    )

for forbidden in (
    "prepare_ow_fixed_center_group",
    "build_dg_group_averaged_occupied_candidates_eigenexa",
    "build_dg_cocycle_averaged_occupied_candidates_eigenexa",
    "build_dg_periodic_spectral_basins",
    "prepare_dg_spectral_basin_operators",
    "propagate_dg_spectral_basin_orbit_channels",
    "build_dg_finite_abelian_character_table",
    "split_dg_translation_character_sector_eigenexa",
    "prepare_dg_translation_character_action",
    "build_dg_translation_character_intertwining_phase_prepared",
    "begin_sawf_dmn",
    "append_sawf_dmn_operation",
    "finish_sawf_dmn",
    "validate_dg_w90_generator_covariance",
    "inherit_dg_w90_affine_receipts",
    "validate_dg_factored_point_cogroup_gauge",
    "verify_dg_wannier_center_affine_orbits",
):
    assert forbidden not in localization_arm, f"localization-first arm still invokes {forbidden}"

assert "dg_w90_constrained" in legacy_arm and "dg_w90_unconstrained" not in legacy_arm, (
    "legacy overlapping-Wannier arm must remain symmetry constrained"
)
assert "ow_pencil_generator_maps,source=global_symmetry_map" in common_tail.replace(" ", ""), (
    "localization-first route must retain the complete crystal maps for LCFO/operator diagnostics"
)

for token in [
    "freeze_dg_hybrid_basis_directory",
    "materialize_dg_hybrid_production_interior",
    "assemble_dg_hybrid_broken_volume_rows",
    "assemble_dg_hybrid_local_potential_rows",
    "assemble_dg_hybrid_divided_nonlocal_rows",
    "freeze_dg_hybrid_variational_payload",
    "materialize_dg_hybrid_production_face_collection",
    "assemble_dg_hybrid_production_interface_component_rows",
    "reconstruct_dg_hybrid_production_interface_actions",
    "evaluate_dg_hybrid_face_action_residuals",
]:
    assert token in source, f"production continuation does not connect {token}"

capture_variable = "SALMON_DG_VARIATIONAL_PAYLOAD_CAPTURE"
assert capture_variable in source, "continuation omits opt-in variational payload capture"
capture_call = "call write_dg_hybrid_variational_payload_bundle"
assert capture_call in source, "continuation omits variational payload bundle writer"
capture_position = source.index(capture_call)
freeze_position = source.index("call freeze_dg_hybrid_variational_payload")
assert capture_position < freeze_position, "variational payload must be captured before freeze"
capture_guard = source[:capture_position].rsplit("if(", 1)[-1]
assert "payload_capture" in capture_guard and ">0" in capture_guard, (
    "variational payload capture must be disabled for an empty environment value"
)
between_capture_and_freeze = source[capture_position:freeze_position]
assert "error stop 'divided Hybrid variational payload capture failed'" in between_capture_and_freeze, (
    "capture failure must stop before payload freeze"
)
assert "[HYBRID-VARIATIONAL-PAYLOAD-CAPTURE]" in source, (
    "continuation omits variational payload capture receipt"
)
for token in (
    "divided_full_basis_closure_defect",
    "divided_full_basis_closure_ok",
    "[HYBRID-RETAINED-BASIS-SYMMETRY]",
):
    assert token.lower() in source.lower(), f"continuation omits full-basis diagnostic {token}"
assert "if(.not.divided_full_basis_closure_ok)error stop" not in source.replace(" ", "").lower(), (
    "full retained-basis nonclosure must not stop LCFO"
)

nonlocal_name = "subroutine assemble_dg_hybrid_divided_nonlocal_rows"
nonlocal_body = source[source.index(nonlocal_name):]
nonlocal_body = nonlocal_body[:nonlocal_body.index("end subroutine assemble_dg_hybrid_divided_nonlocal_rows")]
for token in (
    "ppg%uv",
    "fragment_basis%buffer_values",
    "dc%jxyz_tot",
    "map_dc_atom_to_physical_atom",
    "collect_dg_overlapping_wannier_projector_overlaps",
    "assemble_dg_overlapping_wannier_nonlocal_rows",
):
    assert token in nonlocal_body.lower(), "production nonlocal path omits " + token
assert "complex(8),intent(out)::matrix_rows(:,:)" in nonlocal_body, (
    "production nonlocal output must use a caller-owned fixed-extent matrix"
)
assert "assembled_matrix_rows" in nonlocal_body, (
    "production nonlocal assembly must validate a temporary result before publication"
)
assert "shape(assembled_matrix_rows)/=[size(row_ids),global_count]" in nonlocal_body.replace(" ", ""), (
    "production nonlocal assembly omits exact extent validation"
)
assert "matrix_rows=assembled_matrix_rows" in nonlocal_body.replace(" ", ""), (
    "production nonlocal assembly does not publish the validated temporary"
)
nonlocal_call_position = source.index("call assemble_dg_hybrid_divided_nonlocal_rows")
nonlocal_allocation = "allocate(dg_hybrid_nonlocal_rows(size(divided_lcfo_row_ids),size(divided_effective_ids)))"
assert nonlocal_allocation in source[:nonlocal_call_position].replace("&\n", ""), (
    "continuation must allocate the exact nonlocal row extent before assembly"
)

if "--interface-only" in sys.argv:
    print("PASS production continuation SIPG interface connection")
    raise SystemExit(0)

setup_name = "if(yn_dg_hybrid_divided_scf=='y'.or.yn_dg_hybrid_continuation_scf=='y')then"
setup_start = source.index(setup_name)
setup = source[setup_start:source.index("return", setup_start)]
continuation_name = "if(yn_dg_hybrid_continuation_scf=='y')then"
continuation = setup[setup.index(continuation_name):]
continuation_only = continuation.split("else", 1)[0]
assert "run_dg_hybrid_divided_scf" not in continuation_only, (
    "continuation must start directly from the exact converged DC density"
)
assert "assemble_dg_hybrid_lcfo_rows" not in continuation_only, (
    "continuation metric must not be obtained through ordinary hpsi projection"
)
assert "apply_dg_hybrid_divided_fragment_hpsi" not in continuation_only, (
    "ordinary real-space hpsi is forbidden in the DG variational continuation"
)
assert "dg_hybrid_unit_local_potential" in continuation_only, (
    "continuation must construct S from the frozen basis volume integral"
)
assert "size(divided_lcfo_hrows,2)" not in continuation_only, (
    "continuation must not query an unallocated legacy LCFO Hamiltonian"
)

reuse = source[source.index("if(ok.and.reusable"):].split("endif", 1)[0]
assert "yn_dg_hybrid_continuation_scf/='y'" in reuse, (
    "an old overlapping-Wannier checkpoint must not bypass continuation"
)

preparation_name = "subroutine prepare_dg_hybrid_divided_production_basis"
preparation_start = source.index(preparation_name)
preparation = source[preparation_start:source.index(
    "end subroutine prepare_dg_hybrid_divided_production_basis", preparation_start
)]
assert "analyze_dg_hybrid_lcfo_selection" in preparation, (
    "Hybrid production preparation does not defer symmetry recovery to LCFO"
)
assert "basis_fingerprint_arg" in preparation[
    preparation.index("call analyze_dg_hybrid_lcfo_selection"):
], "LCFO-deferred preparation omits authoritative Wannier provenance"
assert "pw_max_arg" in preparation[
    preparation.index("call analyze_dg_hybrid_lcfo_selection"):
], "LCFO-deferred preparation omits the post-completion PW capacity guard"
assert "call analyze_dg_hybrid_production_selection" not in preparation, (
    "Hybrid production preparation still invokes strict fragment covariance analysis"
)
for token in (
    "fragment_origin_arg",
    "fragment_size_arg",
    "core_ids_arg",
    "core_fragment_ids(point)=fragment",
):
    assert token in preparation, f"Hybrid production preparation omits DC-derived DG input {token}"
handoff = "[HYBRID-LCFO-SYMMETRY-HANDOFF]"
assert handoff in source, "Hybrid production preparation omits LCFO symmetry-handoff receipt"
receipt = source[source.index(handoff):source.index(handoff) + 300]
assert "wannier_fingerprint=" in receipt and "basis_fingerprint" in receipt, (
    "LCFO symmetry-handoff receipt omits the nonzero Wannier fingerprint"
)
assert "production_fingerprint=" in receipt and "divided_pw_fingerprint" in receipt, (
    "LCFO symmetry-handoff receipt omits the nonzero production fingerprint"
)
pw_receipt_marker = "[HYBRID-PW-CUTOFF]"
assert pw_receipt_marker in source, "Hybrid continuation omits the completed PW cutoff receipt"
pw_receipt = source[source.index(pw_receipt_marker):source.index(pw_receipt_marker) + 400]
for token in ("requested=", "effective=", "shell_added=", "orbit_added="):
    assert token in pw_receipt, f"Hybrid PW cutoff receipt omits {token}"
assert "wannier_pw_max" in source[:source.index("contains")], (
    "Hybrid continuation does not import the user PW capacity guard"
)
production_call_start = source.index("call prepare_dg_hybrid_divided_production_basis")
production_call = source[production_call_start:production_call_start + 1800]
assert "ow_reciprocal_operation_maps" in production_call, (
    "Hybrid PW preparation still uses only the LCFO diagnostic generator maps"
)
assert "global_point_rotations(:,:,global_point_representatives)" in production_call.replace(" ", ""), (
    "Hybrid PW preparation omits the complete authoritative point-cogroup rotations"
)

driver_name = "subroutine run_dg_hybrid_concrete_continuation"
assert driver_name in source, "missing contained concrete continuation driver"
start = source.index(driver_name)
implementation = source[start:source.index("end subroutine run_dg_hybrid_concrete_continuation", start)]
assert "384" not in implementation, "production continuation hard-codes the Si64 target state count"

required = [
    "dc_seed_density",
    "fixed_payload",
    "production_faces",
    "dg_dc_update_potential_from_distributed_density",
    "assemble_dg_hybrid_local_potential_rows",
    "compose_dg_hybrid_variational_hamiltonian",
    "reconstruct_dg_hybrid_occupied_state",
    "reconstruct_dg_hybrid_production_interface_state",
    "reject_dg_hybrid_trial",
    "final_refresh_performed",
    "validate_dg_hybrid_ground_state",
    "final_density",
    "final_trace",
    "final_hamiltonian_rows",
    "measure_dg_hybrid_projector_covariance",
    "measure_dg_hybrid_core_density_covariance",
    "dg_hybrid_continuation_state_count",
    "select_dg_hybrid_symmetry_target",
    "evaluate_dg_hybrid_distributed_low_energy_symmetry",
    "symmetry_target_rank_arg",
    "core_symmetry_maps_arg",
]
for token in required:
    assert token in implementation, f"continuation implementation does not use {token}"

for forbidden in [
    "dg_hybrid_production_continuation_adapter",
    "build_dg_hybrid_production_catalog",
    "s_dg_hybrid_continuation_callbacks",
    "solve_dg_hybrid_generalized_once_and_publish",
    "write_overlapping_wannier_occupied_checkpoint",
]:
    assert forbidden not in implementation, f"continuation implementation uses forbidden {forbidden}"

assert implementation.index("accepted_lambda==1d0") < implementation.index(
    "validate_dg_hybrid_ground_state"
), "lambda-one state is published before the final refresh gate"
assert "real(8),parameter::density_damping" not in implementation
assert "continuation_controller%controls%density_damping" in implementation
assert "schedule_dg_hybrid_candidate_checks" in implementation
assert "begin_dg_hybrid_stage_solve" in implementation
assert "complete_dg_hybrid_stage_solve" in implementation
assert "stage_report%gap_shrinking=.false." not in implementation, (
    "adaptive lambda must use the measured occupied-unoccupied gap"
)
assert "stage_report%hermitian_ok=.true." not in implementation
assert "ow_distributed_hermiticity" in implementation
assert "continuation_state_count=nstate+1" not in implementation
state_count_call = implementation[implementation.index("call dg_hybrid_continuation_state_count"):]
state_count_call = state_count_call[:state_count_call.index("local_ok)")]
assert "symmetry_target_rank_arg" in state_count_call, (
    "continuation state-count policy does not receive the material-dependent target rank"
)
assert "solver_coefficients,solver_eigenvalues" in implementation.replace("&\n", ""), (
    "continuation discards the target-window eigensystem before symmetry analysis"
)
for token in ("target_window_solve: do", "unresolved degenerate cluster"):
    assert token in implementation, f"continuation omits adaptive degeneracy-boundary handling {token}"
assert "continuation_state_count=continuation_state_count+1" in implementation.replace(" ", ""), (
    "continuation does not enlarge an unresolved target window"
)
assert "coefficients,occupied_occupations,rho_out" in implementation.replace("&\n", ""), (
    "density reconstruction no longer uses only occupied coefficients and physical occupations"
)
assert "occupied_unoccupied_gap>dg_dc_gs_final_orbital_tolerance" not in implementation
assert "stage_report%occupation_ok=occupation_kernel_ok" in implementation
assert "hamiltonian_hermiticity<=dg_dc_gs_hermiticity_tolerance" in implementation
assert "max(occupied_symmetry_defect,target_symmetry_defect,target_energy_symmetry_defect," in implementation.replace(
    "&\n", ""
), "post-LCFO acceptance does not gate all physical symmetry defects"
assert "symmetry_residual<=dg_ow_symmetry_tolerance" not in implementation.replace(" ", ""), (
    "full retained-Hamiltonian covariance still gates post-LCFO acceptance"
)
continuation_call = source[source.index("call run_dg_hybrid_concrete_continuation"):start]
for token in ("ntarget", "ow_pencil_generator_maps"):
    assert token in continuation_call, f"production route does not pass {token} into concrete continuation"
assert "checkpoint_payload%energy_receipt=[energy%E_tot" not in implementation, (
    "checkpoint energy receipt still serializes the pre-continuation energy"
)
for receipt_token in (
    "[HYBRID-GS-ACCEPTANCE]",
    "seed_identity=",
    "lambda_zero=",
    "requested_target_rank=",
    "extended_target_rank=",
    "occupied_defect=",
    "target_defect=",
    "target_energy_defect=",
    "density_defect=",
    "payload_fingerprint=",
):
    assert receipt_token in implementation, f"production continuation omits acceptance receipt {receipt_token}"
for token in (
    "energy_global_coefficients",
    "fixed_payload%kinetic_rows",
    "fixed_payload%nonlocal_rows",
    "eexc_tmp(energy_ix,energy_iy,energy_iz)",
    "calc_Total_Energy_periodic(dc%mg_tot,ewald,dc%system_tot",
    "checkpoint_energy%E_ion_ion",
    "fixed_payload%kinetic_rows(energy_row,:)+fixed_payload%interface_rows(energy_row,:)",
):
    assert token in implementation, f"final DG energy decomposition omits {token}"

distributed_helper_name = "subroutine evaluate_dg_hybrid_distributed_low_energy_symmetry"
distributed_helper = source[source.index(distributed_helper_name):]
distributed_helper = distributed_helper[:distributed_helper.index(
    "end subroutine evaluate_dg_hybrid_distributed_low_energy_symmetry"
)]
metric_collective = distributed_helper.index("MPI_IN_PLACE,global_metric")
for token in ("minimum_contract", "maximum_contract", "minimum_tolerance_bits", "maximum_tolerance_bits"):
    assert token in distributed_helper[:metric_collective], (
        f"distributed low-energy symmetry omits pre-collective rank agreement {token}"
    )

density_helper_name = "subroutine measure_dg_hybrid_core_density_covariance"
density_helper = source[source.index(density_helper_name):]
density_helper = density_helper[:density_helper.index(
    "end subroutine measure_dg_hybrid_core_density_covariance"
)]
assert "exchange_dg_point_permuted_orbital_rows" in density_helper, (
    "mapped-core density symmetry does not apply the authoritative distributed point map directly"
)
density_measurement = density_helper.index("call exchange_dg_point_permuted_orbital_rows")
for token in ("minimum_operation_count", "maximum_operation_count"):
    assert token in density_helper[:density_measurement], (
        f"mapped-core density symmetry omits operation-count rank agreement {token}"
    )
assert "MPI_MAX" in density_helper, "mapped-core density symmetry does not retain the worst point defect"
assert "maxval(abs(mapped_density-density_probe))" in density_helper.replace(" ", ""), (
    "mapped-core density symmetry does not measure the direct worst-point mismatch"
)
assert "call measure_dg_spatial_basis_covariance" not in density_helper, (
    "mapped-core density symmetry still uses a grid-diluted L2 residual"
)

gs_input = (root / "tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_dg_continuation.in").read_text()
rt_input = (root / "tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_dg_zero_field_rt.in").read_text()
acceptance_runner = (root / "tests/dg/run_dg_hybrid_si64_continuation_rt.py").read_text()
for token in (
    "requested_target_rank", "extended_target_rank", "occupied_defect", "target_defect",
    "target_energy_defect", "density_defect", "full_operator_defect",
):
    assert token in acceptance_runner, f"Si64 acceptance parser omits {token}"
assert "physical_symmetry_defects" in acceptance_runner, (
    "Si64 acceptance parser does not gate the four post-LCFO physical defects"
)
assert "yn_dg_hybrid_continuation_scf='y'" in gs_input
assert "yn_dg_hybrid_divided_scf='n'" in gs_input
assert "yn_rt_dg_hybrid_continuation='y'" in rt_input
assert "ae_shape1='none'" in rt_input
print("PASS concrete DG continuation production route contract")
