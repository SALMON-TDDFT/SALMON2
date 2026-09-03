#!/usr/bin/env python3
from pathlib import Path
import re
import sys

root = Path(__file__).resolve().parents[2]
source = (root / "src/gs/main_dft.f90").read_text()
controller_source = (root / "src/gs/dc/dg_hybrid_continuation_controller.f90").read_text()
symmetry_source = (root / "src/gs/dc/dg_hybrid_low_energy_symmetry.f90").read_text()
certified_basis_source = (root / "src/gs/dc/dg_hybrid_certified_rt_basis.f90").read_text()


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text.lower())

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
assert localization_arm.count("if(size(global_affine_generators)>0)then") >= 2, (
    "identity-only localization still enters zero-operation symmetry diagnostics"
)
assert "lcfo_symmetry_worst_operation=global_identity_operation" in localization_arm, (
    "identity-only localization does not publish a vacuous identity diagnostic"
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
assert "if(size(ow_pencil_generator_maps,2)>0)then" in source.replace(" ", ""), (
    "common generator-only diagnostics are not guarded for identity-only systems"
)
assert "allocate(ow_pencil_generator_representation(ntarget,ntarget,0))" in source.replace(" ", ""), (
    "identity-only systems do not retain an explicit empty nonidentity-generator representation"
)
assert "if(size(ow_pencil_generator_maps,2)>0.and.ok)then" in source.replace(" ", ""), (
    "common map/gradient diagnostics still receive zero generator operations"
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

shared_setup_flag = "if(yn_dg_hybrid_divided_scf=='y'.or.yn_dg_hybrid_continuation_scf=='y')then"
shared_setup_start = ground_state_lower.index(shared_setup_flag)
shared_route = ground_state_lower[shared_setup_start:]
shared_branch_position = shared_route.index("if(yn_dg_hybrid_continuation_scf=='y')then")
shared_prefix = shared_route[:shared_branch_position]
shared_payload_steps = [
    "call freeze_dg_hybrid_basis_directory",
    "call materialize_dg_hybrid_production_interior",
    "call assemble_dg_hybrid_broken_volume_rows",
    "call assemble_dg_hybrid_divided_nonlocal_rows",
    "call materialize_dg_hybrid_production_face_collection",
    "call assemble_dg_hybrid_production_interface_component_rows",
    "call freeze_dg_hybrid_variational_payload",
]
for token in shared_payload_steps:
    assert token in shared_prefix, (
        f"fixed variational payload step remains branch-local: {token}"
    )
shared_positions = [shared_prefix.index(token) for token in shared_payload_steps]
assert shared_positions == sorted(shared_positions), (
    "shared basis/interior/projector/face/payload construction is out of order"
)
assert shared_route.count("call freeze_dg_hybrid_variational_payload") == 1, (
    "divided/reference routes must freeze exactly one shared variational payload"
)
assert "[hybrid-shared-variational-payload]" in shared_prefix, (
    "shared immutable variational payload fingerprint receipt is missing"
)
assert "dg_hybrid_fixed_payload%fingerprint" in shared_prefix, (
    "shared variational payload receipt omits its immutable fingerprint"
)
continuation_arm = shared_route[shared_branch_position:shared_route.index("else", shared_branch_position)]
divided_arm = shared_route[shared_route.index("else", shared_branch_position):]
for arm_name, arm in (("continuation", continuation_arm), ("divided", divided_arm)):
    assert "divided_fixed_payload_fingerprint" in arm, (
        f"{arm_name} route does not consume the shared immutable payload fingerprint"
    )
    assert "call freeze_dg_hybrid_variational_payload" not in arm, (
        f"{arm_name} route refreezes its variational payload"
    )
shared_reuse_guard = source[source.index("if(ok.and.reusable"):].split("endif", 1)[0]
assert "yn_dg_hybrid_continuation_scf/='y'" in shared_reuse_guard, (
    "an old overlapping-Wannier checkpoint must not bypass continuation"
)
assert "yn_dg_hybrid_divided_scf/='y'" in shared_reuse_guard, (
    "an old overlapping-Wannier checkpoint must not bypass divided payload construction"
)

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
assert "dg_hybrid_unit_local_potential" in shared_prefix, (
    "shared divided/reference setup must construct S from the frozen basis volume integral"
)
assert "size(divided_lcfo_hrows,2)" not in continuation_only, (
    "continuation must not query an unallocated legacy LCFO Hamiltonian"
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
implementation_lower = implementation.lower()
implementation_compact = compact(implementation)

final_acceptance_calls = (
    "call solve_dg_hybrid_generalized_complete_once",
    "call record_dg_hybrid_complete_lcfo_solve",
    "call derive_dg_hybrid_occupation_policy",
    "call record_dg_hybrid_occupation_policy",
    "call record_dg_hybrid_unconditional_gates",
    "call certify_dg_hybrid_energy_window",
    "call record_dg_hybrid_spectral_certification",
    "call build_dg_hybrid_certified_rt_basis",
    "call record_dg_hybrid_certified_rt_basis",
    "call authorize_dg_hybrid_v3_publication",
    "call write_rt_dg_hybrid_ground_state_checkpoint",
)
for token in final_acceptance_calls:
    assert token in implementation_lower, f"final LCFO acceptance omits {token}"
assert implementation.count("call solve_dg_hybrid_generalized_complete_once") == 1, (
    "final candidate must publish one complete LCFO eigensystem"
)
final_acceptance_order = [implementation_lower.index(token) for token in final_acceptance_calls]
assert final_acceptance_order == sorted(final_acceptance_order), (
    "final acceptance must record the full solve, Task 9 occupations, unconditional gates, "
    "Task 10 certification, Task 11 basis, v3 authorization, and write in order"
)
assert "symmetrize_dg_distributed_pencil_rows" not in implementation, (
    "final LCFO acceptance must not symmetrize its Hamiltonian, coefficients, occupations, or density"
)

solve_position = implementation_lower.index("call solve_dg_hybrid_generalized_complete_once")
occupation_position = implementation_lower.index("call derive_dg_hybrid_occupation_policy")
certification_position = implementation_lower.index("call certify_dg_hybrid_energy_window")
basis_position = implementation_lower.index("call build_dg_hybrid_certified_rt_basis")
authorize_position = implementation_lower.index("call authorize_dg_hybrid_v3_publication")
write_position = implementation_lower.index("call write_rt_dg_hybrid_ground_state_checkpoint")
final_acceptance = implementation_lower[solve_position:write_position]
occupation_receipt_position = implementation_lower.index("call record_dg_hybrid_occupation_policy")
gate_position = implementation_lower.index("call record_dg_hybrid_unconditional_gates")
final_projector_measurement = implementation_lower.find(
    "call measure_dg_hybrid_projector_covariance", occupation_receipt_position
)
final_density_measurement = implementation_lower.find(
    "call measure_dg_hybrid_core_density_covariance", occupation_receipt_position
)
assert final_projector_measurement >= 0, "final candidate omits occupied-projector covariance measurement"
assert final_density_measurement >= 0, "final candidate omits density covariance measurement"
assert occupation_receipt_position < final_projector_measurement < gate_position, (
    "occupied-projector covariance is not an unconditional final-candidate gate"
)
assert occupation_receipt_position < final_density_measurement < gate_position, (
    "density covariance is not an unconditional final-candidate gate"
)
assert implementation_compact.count("if(size(basis_representation,3)==0)then") >= 2, (
    "identity-only systems are not handled as a vacuous physical symmetry gate"
)
assert "complete_eigensystem%eigenvalues" in final_acceptance[:certification_position - solve_position], (
    "Task 9 and Task 10 must consume the complete LCFO spectrum"
)
for token in (
    "occupation_result%noccupied",
    "occupation_result%occupations",
    "occupation_result%e_homo",
    "complete_eigensystem%coefficients",
):
    assert token in final_acceptance[occupation_position - solve_position:], (
        f"final LCFO state is not reconstructed from the authoritative Task 9 result: {token}"
    )
assert "occupied_occupations" not in final_acceptance[occupation_position - solve_position:], (
    "final acceptance must not reuse the pre-LCFO occupied-state count or occupations"
)

named_v3_fields = (
    "checkpoint_payload%construction_catalog%valid=.true.",
    "checkpoint_payload%certified_basis%valid=.true.",
    "checkpoint_payload%electron_count%valid=.true.",
    "checkpoint_payload%energy_window%valid=.true.",
    "checkpoint_payload%symmetry_receipt%valid=.true.",
    "checkpoint_payload%rt_space%valid=.true.",
    "checkpoint_payload%handoff_receipts%valid=.true.",
)
for token in named_v3_fields:
    assert token in implementation_compact, f"v3 publication omits {token}"
assert max(implementation_compact.index(token) for token in named_v3_fields) < implementation_compact.index(
    "callauthorize_dg_hybrid_v3_publication"
), "v3 publication is authorized before every named payload section is valid"
assert "occupied_occupations" not in implementation_lower, (
    "pre-LCFO occupations still influence the final continuation candidate"
)
assert "allocate(checkpoint_payload%effective_ids,source=selection_effective_ids_arg)" in implementation_compact, (
    "checkpoint selection receipt no longer preserves the PW packet-ID domain"
)
assert "construction_ids=int(effective_ids,8)" in implementation_compact, (
    "v3 construction catalog uses PW packet IDs instead of construction-column IDs"
)
assert "construction_ids=int(selection_effective_ids_arg,8)" not in implementation_compact
for fingerprint in (
    "catalog_ids_fingerprint", "catalog_generation_fingerprint", "catalog_ordering_fingerprint",
    "catalog_ownership_fingerprint", "catalog_provenance_fingerprint",
):
    assert re.search(
        rf"catalog_fingerprint=ieor\(ishftc\(catalog_fingerprint,7\),{fingerprint}\)",
        implementation_compact,
    ), f"aggregate construction-catalog fingerprint omits {fingerprint}"
assert "checkpoint_payload%construction_catalog%catalog_fingerprint=catalog_fingerprint" in implementation_compact
assert "checkpoint_payload%catalog_fingerprint=catalog_fingerprint" in implementation_compact
assert "checkpoint_payload%catalog_fingerprint=fixed_payload%basis_fingerprint" not in implementation_compact
assert "real(size(checkpoint_payload%nonlocal_rows),8)" not in implementation_compact, (
    "replicated pseudopotential receipt contains a rank-local matrix extent"
)
assert "callcheckpoint_grid_real_fingerprint(dc%icomm_tot,ow_core_ids,dc_seed_density" in implementation_compact, (
    "DC-seed fingerprint is not keyed by global grid identity"
)
grid_real_fingerprint = compact(
    source.split("subroutine checkpoint_grid_real_fingerprint", 1)[1].split(
        "end subroutine checkpoint_grid_real_fingerprint", 1
    )[0]
)
assert "mod(13*p,63)" not in grid_real_fingerprint, (
    "grid-real fingerprint still depends on rank-local point ordering"
)
assert "grid_ids_arg(p),63_8" in grid_real_fingerprint and "global_count" in grid_real_fingerprint, (
    "grid-real fingerprint omits global-ID keys or the collective grid count"
)

for token in (
    "checkpoint_payload%rt_space%rank=spectral_certification%certified_rank",
    "checkpoint_payload%certified_basis%certified_count=spectral_certification%certified_rank",
    "checkpoint_payload%energy_window%construction_rank=size(effective_ids)",
    "checkpoint_payload%energy_window%solved_rank=size(effective_ids)",
    "checkpoint_payload%energy_window%occupied_rank=occupation_result%noccupied",
    "checkpoint_payload%energy_window%requested_rank=spectral_certification%requested_rank",
    "checkpoint_payload%energy_window%certified_rank=spectral_certification%certified_rank",
    "checkpoint_payload%energy_window%boundary_cluster_rank=spectral_certification%boundary_cluster_rank",
    "rt_dg_hybrid_ground_state_checkpoint_version",
):
    assert token in implementation_compact, f"v3 certified-rank publication omits {token}"
for forbidden_rank in (
    "checkpoint_payload%rt_space%rank=size(effective_ids)",
    "checkpoint_payload%certified_basis%certified_count=size(effective_ids)",
    "checkpoint_payload%rt_space%rank=continuation_state_count",
):
    assert forbidden_rank not in implementation_compact, (
        "construction-only directions leaked into the published RT dimensions"
    )

assert authorize_position < write_position, "checkpoint write precedes v3 authorization"
for receipt in ("[HYBRID-LCFO-WINDOW]", "[HYBRID-LCFO-SYMMETRY]", "[HYBRID-RT-BASIS]"):
    assert receipt in implementation, f"final LCFO acceptance omits receipt {receipt}"
window_receipt_position = implementation.index("[HYBRID-LCFO-WINDOW]")
symmetry_receipt_position = implementation.index("[HYBRID-LCFO-SYMMETRY]")
rt_receipt_position = implementation.index("[HYBRID-RT-BASIS]")
window_receipt_block = implementation[window_receipt_position:symmetry_receipt_position]
symmetry_receipt_block = implementation[symmetry_receipt_position:rt_receipt_position]
rt_receipt_block = implementation[rt_receipt_position:implementation.index("[OW-GS]", rt_receipt_position)]
for field in (
    "mode=", "compatibility=", "delta_e=", "homo=", "requested_cutoff=", "requested_rank=",
    "certified_cutoff=", "certified_rank=", "boundary_rank=", "extension_states=", "extension_energy=",
    "proof_state=", "proof_energy=",
):
    assert field in window_receipt_block, f"production LCFO-window receipt omits {field}"
for field in ("occupied=", "target=", "energy=", "density=", "worst_operation="):
    assert field in symmetry_receipt_block, f"production LCFO-symmetry receipt omits {field}"
for field in (
    "construction_rank=", "certified_rank=", "rt_rank=", "localization_spread=",
    "embedding_fingerprint=", "operator_covariance=", "checkpoint_version=",
):
    assert field in rt_receipt_block, f"production RT-basis receipt omits {field}"
assert "if(energy_window==-1d0)then" in compact(symmetry_source), (
    "exact -1 legacy dynamic-rank branch is missing"
)
assert "energy_window<0d0.and.energy_window/=-1d0" in compact(symmetry_source), (
    "negative non-sentinel Hybrid windows are not rejected exactly"
)
assert (
    "[HYBRID-SYMMETRY-COMPATIBILITY-WARNING] energy_window=-1 uses dynamic requested-rank selection"
    in symmetry_source
), (
    "legacy dynamic-rank branch omits its explicit warning"
)
for token in (
    "dg_hybrid_symmetry_energy_window==-1d0",
    "rt_dg_hybrid_energy_window_legacy_dynamic",
    "rt_dg_hybrid_energy_window_explicit",
):
    assert token in implementation_compact, f"final publication omits exact -1 compatibility receipt {token}"
for token in (
    "s_dg_hybrid_candidate_acceptance",
    "initialize_dg_hybrid_candidate_acceptance",
    "record_dg_hybrid_complete_lcfo_solve",
    "record_dg_hybrid_occupation_policy",
    "record_dg_hybrid_unconditional_gates",
    "record_dg_hybrid_spectral_certification",
    "record_dg_hybrid_certified_rt_basis",
    "authorize_dg_hybrid_v3_publication",
    "legacy_warning_required",
    "legacy_warning_observed",
    "published_rt_rank",
):
    assert token in controller_source, f"controller omits ordered final-acceptance contract {token}"

projection_helper_match = re.search(
    r"subroutine\s+(?P<name>[a-z0-9_]*certified[a-z0-9_]*representation[a-z0-9_]*)\b"
    r"(?P<body>.*?)end\s+subroutine",
    source,
    re.I | re.S,
)
assert projection_helper_match, "missing C^dagger S O_R C certified-representation projector"
projection_helper_name = projection_helper_match.group("name").lower()
projection_helper = compact(projection_helper_match.group("body"))
assert f"call{projection_helper_name}" in implementation_compact, (
    "final acceptance does not call the certified-representation projector"
)
assert "conjg(transpose(" in projection_helper, (
    "certified representation omits the C^dagger projection"
)
assert re.search(r"matmul\([a-z0-9_]*metric[a-z0-9_]*,", projection_helper), (
    "certified representation omits the construction metric S"
)
assert re.search(
    r"matmul\([a-z0-9_]*representation[a-z0-9_]*\(:,:,operation\),",
    projection_helper,
), "certified representation omits the construction action O_R"
operator_route_compact = implementation_compact.replace("&", "")
position_assembly_position = operator_route_compact.index("callassemble_dg_cell_wrapped_position")
assert "checkpoint_position(p,:,:),projected_position" not in operator_route_compact, (
    "affine cell-wrapped position was incorrectly sent through the homogeneous-vector gate"
)
for token in (
    "rt_dg_hybrid_vector_canonical_momentum",
    "callassemble_dg_hybrid_canonical_momentum",
    "interior_gradients,checkpoint_momentum",
    "checkpoint_momentum(p,:,:),projected_position",
    "certified_vector_operators(:,:,p,rt_dg_hybrid_vector_canonical_momentum)=projected_position",
    "checkpoint_payload%rt_space%vector_count=size(certified_rt_basis%vector_operators_rt,4)",
    "certified_rt_basis%vector_operators_rt(construction_index,:,:,:)",
    "checkpoint_cartesian_rotations(:,:,p+1)=transpose(cartesian_rotations_arg(:,:,p))",
):
    assert token in operator_route_compact, f"canonical-momentum vector handoff omits {token}"
momentum_gate_position = operator_route_compact.index("callassemble_dg_hybrid_canonical_momentum")
assert position_assembly_position < momentum_gate_position < operator_route_compact.index(
    "callbuild_dg_hybrid_certified_rt_basis"
), "canonical momentum is not gated before certified RT publication"
assert "fixed_payload%interface_fingerprint,dg_dc_gs_final_orbital_tolerance" not in operator_route_compact
assert (
    "final_operator_fingerprint,checkpoint_payload%position_convention_fingerprint,"
    "dg_dc_gs_final_orbital_tolerance"
) in operator_route_compact, "final ground state is not bound to the actual position convention"

for token in (
    "certified_scalar_operators(certified_rank,certified_rank,5)",
    "fixed_payload%kinetic_rows,full_component",
    "fixed_payload%nonlocal_rows,full_component",
    "iterate%local_rows,full_component",
    "fixed_payload%interface_rows,full_component",
    "iterate%hamiltonian_rows,full_component",
    "rt_kinetic=certified_rt_basis%scalar_operators_rt(:,:,1)",
    "rt_nonlocal=certified_rt_basis%scalar_operators_rt(:,:,2)",
    "rt_local=certified_rt_basis%scalar_operators_rt(:,:,3)",
    "rt_sipg=certified_rt_basis%scalar_operators_rt(:,:,4)",
):
    assert token in operator_route_compact, f"component-wise fixed-operator gate omits {token}"

for token in (
    "size(basis_representation,3)+1",
    "checkpoint_basis_representation(p,p,1)=(1d0,0d0)",
    "checkpoint_basis_representation(:,:,2:)=basis_representation",
    "checkpoint_cartesian_rotations(p,p,1)=1d0",
    "checkpoint_cartesian_rotations(:,:,p+1)=transpose(cartesian_rotations_arg(:,:,p))",
    "full_metric,checkpoint_basis_representation",
    "certified_representation,checkpoint_cartesian_rotations",
    "checkpoint_payload%nonidentity_operation_count=size(basis_representation,3)",
):
    assert token in operator_route_compact, f"checkpoint symmetry catalog omits explicit identity: {token}"

for token in (
    "callcheckpoint_grid_integer_fingerprint(dc%icomm_tot,ow_core_ids,",
    "checkpoint_payload%rt_space%grid_owner_keys,rt_grid_ownership_fingerprint",
    "rt_ownership_fingerprint=ieor(ishftc(rt_ownership_fingerprint,7),rt_grid_ownership_fingerprint)",
):
    assert token in operator_route_compact, f"RT ownership receipt omits {token}"
for fingerprint in (
    "rt_metric_fingerprint", "rt_kinetic_fingerprint", "rt_nonlocal_fingerprint",
    "rt_local_fingerprint", "rt_sipg_fingerprint", "rt_hamiltonian_fingerprint",
    "rt_basis_fingerprint", "rt_density_fingerprint", "rt_ownership_fingerprint",
    "rt_scalar_fingerprint", "rt_vector_fingerprint", "rt_tensor_fingerprint",
    "rt_representation_fingerprint", "rt_rotation_fingerprint",
):
    assert re.search(
        rf"rt_payload_fingerprint=ieor\(ishftc\(rt_payload_fingerprint,7\),{fingerprint}\)",
        operator_route_compact,
    ), f"aggregate RT fingerprint omits {fingerprint}"
for field in (
    "occupied_subspace_defect", "occupied_projector_defect", "target_subspace_defect",
    "target_energy_defect", "density_defect", "scalar_covariance_defect",
    "vector_covariance_defect", "tensor_covariance_defect", "final_basis_defect",
    "worst_operation_defect", "maximum_physical_defect",
):
    assert (
        f"checkpoint_payload%symmetry_receipt%{field}" in operator_route_compact[
            operator_route_compact.index("symmetry_receipt_fingerprint=checkpoint_real_fingerprint"):
            operator_route_compact.index("checkpoint_payload%symmetry_receipt%fingerprint")
        ]
    ), f"aggregate symmetry receipt omits {field}"
handoff_fingerprint_block = operator_route_compact[
    operator_route_compact.index("handoff_fingerprint=ieor("):
    operator_route_compact.index("checkpoint_payload%handoff_receipts%fingerprint")
]
assert "checkpoint_payload%pseudopotential_fingerprint" in handoff_fingerprint_block, (
    "aggregate handoff fingerprint omits pseudopotential provenance"
)
assert "dynamic_requested_rank" in operator_route_compact
assert "min(extended_target_rank,size(effective_ids))" not in operator_route_compact, (
    "legacy Task 10 receipt uses the stage-expanded rank instead of the requested rank"
)

basis_call_end = implementation_lower.index(
    "call record_dg_hybrid_certified_rt_basis", basis_position
)
basis_call = compact(implementation[basis_position:basis_call_end])
for token in (
    "size(effective_ids)",
    "fixed_payload%metric_rows",
    "certified_rank",
    "certified_representation",
):
    assert token in basis_call, f"Task 11 production call omits certified input {token}"
assert "certified_rank=spectral_certification%certified_rank" in implementation_compact, (
    "Task 11 rank is not the Task 10 certified rank"
)
assert "complete_eigensystem%coefficients" in compact(
    implementation[certification_position:basis_position]
), "Task 11 C_cert is not taken from the complete LCFO eigensystem"
assert "result%b_rt=right_transform_rows(c_cert_rows,transform)" in compact(certified_basis_source), (
    "Task 11 no longer enforces B_rt = C_cert U_rt"
)
assert re.search(
    r"matmul\(transpose\([^)]*b_rt[^)]*\),interior_values\)",
    implementation_compact,
), "final real-space RT basis is not reconstructed from B_rt"

localizer_match = re.search(
    r"subroutine\s+(?P<name>[a-z0-9_]*certified[a-z0-9_]*localiz[a-z0-9_]*)\b"
    r"(?P<header>.*?require_unconstrained.*?)"
    r"(?P<body>.*?)end\s+subroutine",
    source,
    re.I | re.S,
)
assert localizer_match, "missing production fixed-rank certified-space localizer"
localizer_name = localizer_match.group("name").lower()
localizer = compact(localizer_match.group("header") + localizer_match.group("body"))
assert localizer_name in basis_call, "Task 11 is not connected to the production certified localizer"
for token in (
    "require_unconstrained",
    "dg_w90_unconstrained",
    "certified_rank,certified_rank",
    "spreads_before",
    "spreads_after",
    "symmetry_constrained=.false.",
    "require_nonincreasing_spread=.false.",
):
    assert token in localizer, f"certified localizer is not unconstrained and fixed rank: {token}"

for token in (
    "rt_row_ids",
    "checkpoint_payload%certified_basis%transformation_row_ids",
    "checkpoint_payload%rt_space%row_ids",
    "checkpoint_payload%rt_space%row_owner_keys",
    "checkpoint_payload%rt_space%grid_owner_keys",
):
    assert token in implementation_lower, f"v3 producer omits unique RT row ownership: {token}"
assert re.search(r"mod\([a-z0-9_]+-1,[a-z0-9_]+\)==rank_local", implementation_compact), (
    "v3 producer does not assign each replicated RT row to one cyclic owner"
)
assert re.search(
    r"transformation_row_ids(?:\([^)]*\))?(?:,source=|=)rt_row_ids",
    implementation_compact,
), "U_rt transformation rows are not uniquely row-owned"
assert re.search(
    r"rt_space%row_ids(?:\([^)]*\))?(?:,source=|=)rt_row_ids",
    implementation_compact,
), "RT operator rows are not uniquely row-owned"

for offsets in (
    "face_offsets",
    "face_weight_offsets",
    "face_basis_offsets",
    "face_value_offsets",
    "face_observable_offsets",
):
    assert f"checkpoint_payload%{offsets}" in implementation_lower, (
        f"v3 face handoff omits {offsets}"
    )
for token in (
    "checkpoint_payload%face_offsets(face_slot+1)=face_point_position+1",
    "checkpoint_payload%face_weight_offsets(face_slot+1)=face_weight_position+1",
    "checkpoint_payload%face_basis_offsets(face_slot+1)=face_basis_position+1",
    "checkpoint_payload%face_value_offsets(face_slot+1)=face_value_position+1",
    "checkpoint_payload%face_observable_offsets(face_slot+1)=face_observable_position+1",
    "face_observable_count",
    "interface_cursor",
    "transpose(interface_state(",
):
    assert token in implementation_compact, f"v3 canonical face packing omits {token}"
assert "allocate(checkpoint_payload%interface_observables,source=interface_state)" not in implementation_compact, (
    "v3 face observables still publish every local duplicate instead of canonical-owner entries"
)

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
assert "accepted_lambda=0d0;lambda_zero_accepted=.true." in implementation_compact, (
    "controller initialization does not record the accepted lambda-zero state"
)
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
for obsolete_prefix_route in (
    "target_window_solve",
    "continuation_state_count=continuation_state_count+1",
    "callselect_dg_hybrid_symmetry_target",
    "calldg_hybrid_continuation_state_count",
):
    assert obsolete_prefix_route not in implementation_compact, (
        f"continuation still uses the obsolete prefix eigensolve route: {obsolete_prefix_route}"
    )
stage_target_setup = implementation_compact[
    implementation_compact.index("occupied_symmetry_rank=occupation_result%noccupied"):
    implementation_compact.index("if(allocated(coefficients))")
]
assert "extended_target_rank=occupation_result%noccupied" in stage_target_setup, (
    "explicit energy-window continuation still starts from a legacy target rank"
)
legacy_target_guard = "if(dg_hybrid_symmetry_energy_window==-1d0)then"
assert legacy_target_guard in stage_target_setup, (
    "dynamic requested-rank setup is not restricted to the exact -1 compatibility mode"
)
legacy_guard_position = stage_target_setup.index(legacy_target_guard)
legacy_guard_end = stage_target_setup.index("endif", legacy_guard_position)
assert "symmetry_target_rank_arg" not in stage_target_setup[:legacy_guard_position]
assert "symmetry_target_rank_arg" in stage_target_setup[legacy_guard_position:legacy_guard_end], (
    "exact -1 compatibility mode no longer receives its dynamic requested rank"
)
stage_symmetry_gate = implementation_compact[
    implementation_compact.index("density_symmetry_defect=max(input_density_symmetry_defect"):
    implementation_compact.index("callow_distributed_hermiticity")
].replace("&", "")
assert "physical_symmetry_defect=max(occupied_symmetry_defect,density_symmetry_defect)" in stage_symmetry_gate, (
    "pre-certification SCF still requires an arbitrary target prefix to be symmetry closed"
)
assert (
    "if(dg_hybrid_symmetry_energy_window==-1d0)physical_symmetry_defect="
    "max(physical_symmetry_defect,target_symmetry_defect,target_energy_symmetry_defect)"
) in stage_symmetry_gate, "legacy dynamic-rank target gate was not isolated to exact -1 mode"
assert "stage_report%symmetry_ok=physical_symmetry_defect<=dg_ow_symmetry_tolerance" in implementation_compact, (
    "stage acceptance bypasses the mode-aware physical symmetry gate"
)
assert "spreads_before(1)=spreads_before(1)+" not in implementation_compact, (
    "certified localizer rewrites the measured pre-localization spread"
)
assert "sum(spreads_before)+dg_ow_symmetry_tolerance>=sum(spreads_after)" not in implementation_compact, (
    "certified localizer rejects a signed spread receipt using a material-independent threshold"
)
complete_call_tail = implementation_lower[solve_position:occupation_position]
assert "solve_dg_hybrid_generalized_scalapack" in complete_call_tail, (
    "complete-once wrapper is not connected to the production generalized eigensolver"
)
assert "size(effective_ids)" in compact(complete_call_tail), (
    "complete-once solve is not dimensioned by the full construction basis"
)
assert "occupied_unoccupied_gap>dg_dc_gs_final_orbital_tolerance" not in implementation
assert "stage_report%occupation_ok=occupation_kernel_ok" in implementation
assert "hamiltonian_hermiticity<=dg_dc_gs_hermiticity_tolerance" in implementation
post_certification = implementation_compact[
    implementation_compact.index("callcertify_dg_hybrid_energy_window"):
].replace("&", "")
assert (
    "physical_symmetry_defect=max(occupied_symmetry_defect,target_symmetry_defect,"
    "target_energy_symmetry_defect,density_symmetry_defect)"
) in post_certification, "post-LCFO acceptance does not gate all physical symmetry defects"
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
for token in (
    "LCFO_WINDOW_RECEIPT",
    "LCFO_SYMMETRY_RECEIPT",
    "RT_BASIS_RECEIPT",
    "LEGACY_WARNING",
    "certified_rank",
    "certified_rt_rank",
    "embedding_fingerprint",
    "operator_covariance",
    "checkpoint_version",
):
    assert token in acceptance_runner, f"Si64 parser omits certified-v3 evidence {token}"
assert "384" not in acceptance_runner, "Si64 parser still embeds a material-specific LCFO rank fixture"
assert "yn_dg_hybrid_continuation_scf='y'" in gs_input
assert "yn_dg_hybrid_divided_scf='n'" in gs_input
assert "yn_rt_dg_hybrid_continuation='y'" in rt_input
assert "ae_shape1='none'" in rt_input
print("PASS concrete DG continuation production route contract")
