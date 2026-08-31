#!/usr/bin/env python3
from pathlib import Path
import sys

root = Path(__file__).resolve().parents[2]
source = (root / "src/gs/main_dft.f90").read_text()

flag = "yn_dg_hybrid_continuation_scf == 'y'"
assert flag in source, "missing explicit DG continuation production branch"
branch = source[source.index(flag):]
branch = branch[:branch.index("endif")]
assert "run_dg_hybrid_continuation_ground_state_for_main" in branch

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
assert "yn_dg_hybrid_continuation_scf='y'" in gs_input
assert "yn_dg_hybrid_divided_scf='n'" in gs_input
assert "yn_rt_dg_hybrid_continuation='y'" in rt_input
assert "ae_shape1='none'" in rt_input
print("PASS concrete DG continuation production route contract")
