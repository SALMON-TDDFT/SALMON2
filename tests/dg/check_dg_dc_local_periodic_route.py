#!/usr/bin/env python3
import re
from pathlib import Path

root = Path(__file__).resolve().parents[2]
global_source = (root / "src/io/salmon_global.f90").read_text()
input_source = (root / "src/io/inputoutput.f90").read_text()
scf_source = (root / "src/gs/scf_iteration_dft.f90").read_text()
main_source = (root / "src/gs/main_dft.f90").read_text()
main_driver_source = main_source.split("contains", 1)[0]
scf_kernel_source = (root / "src/gs/scf_iteration.f90").read_text()
cg_source = (root / "src/gs/conjugate_gradient.f90").read_text()
dc_source = (root / "src/gs/dc/dcdft.f90").read_text()
handoff_source = (root / "src/gs/dc/dg_dc_handoff.f90").read_text()
buffer_projector_source = (root / "src/common/dg_buffer_window_projector.f90").read_text()
ground_state_source = (root / "src/gs/dc/dg_dc_ground_state.f90").read_text()
adapter_source = (root / "src/gs/dc/dg_dc_ground_state_adapter.f90").read_text()
cmake_source = (root / "src/gs/dc/CMakeLists.txt").read_text()

controls = {
    "yn_dg_dc_local_periodic": r"yn_dg_dc_local_periodic\s*=\s*'n'",
    "dg_dc_handoff_min_iter": r"dg_dc_handoff_min_iter\s*=\s*[1-9]\d*",
    "dg_dc_handoff_tolerance": r"dg_dc_handoff_tolerance\s*=\s*[0-9.d+\-]+",
    "dg_dc_candidate_orbitals_per_atom": r"dg_dc_candidate_orbitals_per_atom\s*=\s*40\b",
    "dg_dc_metric_rank_tolerance": r"dg_dc_metric_rank_tolerance\s*=\s*[0-9.d+\-]+",
}
for name, default_pattern in controls.items():
    assert re.search(rf"\b{name}\b", global_source, re.I), f"missing global control {name}"
    assert re.search(rf"namelist/dc/.*\b{name}\b", input_source, re.I | re.S), f"{name} absent from &dc"
    assert re.search(default_pattern, input_source, re.I), f"missing default for {name}"
    assert re.search(rf"call\s+comm_bcast\s*\(\s*{name}\b", input_source, re.I), f"missing broadcast for {name}"

assert re.search(r"call\s+yn_argument_check\s*\(\s*yn_dg_dc_local_periodic\s*\)", input_source, re.I)
for name in controls.keys() - {"yn_dg_dc_local_periodic"}:
    assert re.search(rf"{name}\s*(<=|<)\s*", input_source, re.I), f"{name} must be validated even disabled"

assert "dg_dc_handoff.f90" in cmake_source
assert "dg_buffer_window_projector.f90" in (root / "src/common/CMakeLists.txt").read_text()
assert "occupation" not in buffer_projector_source.lower(), \
    "the full configured window projector must not apply occupation weighting"
for projector_contract in [
    "S_dg_buffer_projector_diagnostics",
    "build_dg_buffer_window_projector",
    "retained_rank",
    "minimum_retained_singular_value",
    "projection_residual",
    "escape_norm",
]:
    assert projector_contract.lower() in buffer_projector_source.lower(), \
        f"buffer-window projector is missing {projector_contract}"
assert re.search(r"call\s+evaluate_dg_dc_handoff", scf_source, re.I)
handoff_accept = re.search(
    r"if\s*\(\s*dg_handoff_accept\s*\)\s*then(.*?)exit\s+DFT_Iteration",
    scf_source, re.I | re.S,
)
assert handoff_accept and re.search(
    r"call\s+timer_end\s*\(\s*LOG_WRITE_GS_RESULTS\s*\)", handoff_accept.group(1), re.I
), "accepted DC handoff must close the active results timer"
assert not re.search(r"call\s+materialize_dg_dc_candidates", main_driver_source, re.I), \
    "coefficient DG route must not materialize a global nodal state"
assert not re.search(r"call\s+solve_nodal_ground_state_cg_mpi", main_driver_source, re.I), \
    "production DG route must use distributed local-basis coefficients"
assert re.search(
    r"nstate\s*=\s*dg_dc_candidate_orbitals_per_atom\s*\*\s*natom", dc_source, re.I
), "enabled DC solve must create the complete configured untruncated candidate pool"
assert re.search(r"call\s+preserve_dg_dc_density_potential", scf_source, re.I)
assert re.search(r"call\s+discard_dc_mixing_history", scf_source, re.I)
assert "candidate_metric_rank" in handoff_source and "maxloc" in handoff_source
assert re.search(r"call\s+gram_schmidt\s*\(", scf_kernel_source, re.I), \
    "existing Gram-Schmidt must remain the orthogonalization owner"
for token in [
    "DG_DC_PRECONVERGENCE",
    "DG_DC_CONTINUATION",
    "DG_DC_FULL_SCF",
    "DG_DC_UNMIXED_FIXED_POINT",
    "DG_DC_ACCEPTED",
    "DG_DC_FAILED",
    "lambda_step=0.5d0*lambda_step",
    "allocate(accepted_state,source=state",
    "accepted_density=density",
    "restore_nodal_state(state,accepted_state)",
    "density=accepted_density",
]:
    assert token.replace(" ", "").lower() in ground_state_source.replace(" ", "").lower(), token

gs_controls = [
    "dg_dc_gs_intermediate_orbital_tolerance",
    "dg_dc_gs_intermediate_density_tolerance",
    "dg_dc_gs_final_orbital_tolerance",
    "dg_dc_gs_final_density_tolerance",
    "dg_dc_gs_subspace_tolerance",
    "dg_dc_gs_initial_lambda_step",
    "dg_dc_gs_minimum_lambda_step",
    "dg_dc_gs_maximum_lambda_step",
    "dg_dc_gs_allowed_residual_growth",
    "dg_dc_gs_density_mix_rate",
    "dg_dc_gs_sipg_penalty_factor",
    "dg_dc_gs_target_lambda",
    "dg_dc_gs_hermiticity_tolerance",
    "dg_dc_gs_orthogonality_tolerance",
    "dg_dc_gs_face_balance_tolerance",
    "dg_dc_gs_electron_count_tolerance",
    "dg_dc_gs_minimum_projector_overlap",
    "dg_dc_gs_maximum_scf_iterations",
    "dg_dc_gs_maximum_eigensolver_iterations",
    "dg_dc_gs_maximum_rollbacks",
]
for name in gs_controls:
    assert re.search(rf"\b{name}\b", global_source, re.I), name
    assert re.search(rf"namelist/dc/.*\b{name}\b", input_source, re.I | re.S), name
    assert re.search(rf"call\s+comm_bcast\s*\(\s*{name}\b", input_source, re.I), name

assert "dg_dc_ground_state.f90" in cmake_source
assert "dg_dc_ground_state_adapter.f90" in cmake_source
assert re.search(r"call\s+run_dg_dc_local_basis_ground_state_for_main\s*\(", main_driver_source, re.I), \
    "production DG route must enter the fragment-local-basis coefficient solver"
local_basis_driver = re.search(
    r"subroutine\s+run_dg_dc_local_basis_ground_state_for_main\b(.*?)"
    r"end\s+subroutine\s+run_dg_dc_local_basis_ground_state_for_main",
    main_source, re.I | re.S,
)
assert local_basis_driver, "production local-basis adapter must be defined"
for required_call in [
    "orthonormalize_dg_dc_fragment_core_basis",
    "transform_dg_dc_fragment_buffer_basis",
    "project_dg_dc_local_basis_volume",
    "assemble_dg_dc_local_basis_interface_rows",
    "compose_dg_dc_distributed_hamiltonian_rows",
    "solve_dg_dc_local_basis_bands_cg",
    "reconstruct_dg_dc_local_basis_density",
    "validate_dg_dc_local_basis_density",
    "calc_rho_total_dcdft",
    "dg_dc_update_potential_from_density",
]:
    assert re.search(rf"call\s+{required_call}\s*\(", local_basis_driver.group(1), re.I), \
        f"production local-basis adapter must call {required_call}"
assert not re.search(r"nodal|lcfo|eigenexa", local_basis_driver.group(1), re.I), \
    "production local-basis adapter must not fall back to nodal, LCFO, or EigenExa paths"
assert re.search(r"call\s+comm_logical_and\s*\(\s*ok\b", local_basis_driver.group(1), re.I), \
    "rank-local stages must gate failure collectively"
assert "81d0" not in local_basis_driver.group(1), \
    "production SIPG penalty must come from a validated control"
assert re.search(r"dg_dc_local_basis_state", local_basis_driver.group(1), re.I), \
    "converged coefficients/eigenvalues must persist in an explicit DG-only state"
assert re.search(r"spsi_shape\s*\(\s*1\s*:\s*3\s*\)\s*==\s*sttpsi_shape",
                 local_basis_driver.group(1), re.I), \
    "production adapter must reject mismatched full-buffer orbital shapes"
for persisted_field in ["full_fragment_basis", "basis_transform", "basis_offsets", "fragment_ids"]:
    assert re.search(rf"candidate_state%{persisted_field}", local_basis_driver.group(1), re.I), \
        f"DG-only result must transactionally retain {persisted_field}"
assert re.search(r"iteration_operator_fingerprint\s*=\s*dg_dc_operator_fingerprint\s*\(\s*\)",
                 local_basis_driver.group(1), re.I), \
    "persisted state must identify the Hamiltonian used by the final coefficient solve"
for continuation_control in [
    "dg_dc_gs_initial_lambda_step",
    "dg_dc_gs_minimum_lambda_step",
    "dg_dc_gs_maximum_lambda_step",
    "dg_dc_gs_allowed_residual_growth",
    "dg_dc_gs_maximum_rollbacks",
]:
    assert continuation_control in local_basis_driver.group(1), \
        f"production local-basis continuation must use {continuation_control}"
assert re.search(r"trial_lambda", local_basis_driver.group(1), re.I), \
    "production local-basis Hamiltonian must use a variable continuation scale"
assert re.search(r"rollback", local_basis_driver.group(1), re.I), \
    "production local-basis continuation must rollback rejected trials"
assert re.search(r"unmixed", local_basis_driver.group(1), re.I), \
    "production local-basis route must perform an unmixed fixed-point gate"
assert re.search(r"yn_dg_dc_local_periodic\s*/=\s*'y'.*?checkpoint_gs",
                 main_source, re.I | re.S), \
    "pre-DG structure checkpoints must be disabled for the opt-in route"
assert re.search(r"yn_dg_dc_local_periodic\s*==\s*'y'.*?time_shutdown\s*>\s*0",
                 input_source, re.I | re.S), \
    "the no-checkpoint route must reject shutdown-checkpoint mode up front"
assert re.search(r"local_basis_route_active\s*=\s*\.true\.", main_source, re.I), \
    "the DG-only result must mark the standard publication path inactive"
assert re.search(r"if\s*\(\s*\.not\.\s*local_basis_route_active\s*\)\s*then.*?"
                 r"checkpoint_gs", main_source, re.I | re.S), \
    "standard checkpoint publication must be gated off for the DG-only state"
assert "solve_dg_dc_local_basis_bands_reference" not in main_source, \
    "production DG route must not call the replicated root reference eigensolver"
assert not re.search(r"call\s+publish_dg_dc_ground_state_for_main\s*\(", main_driver_source, re.I)
assert re.search(r"call\s+calc_rho_total_dcdft\s*\(", main_source, re.I)
assert re.search(r"call\s+calc_vlocal_fragment_dcdft\s*\(", main_source, re.I)
assert re.search(r"hpsi\s*=\s*hpsi_volume_nonlocal\s*\+\s*lambda\s*\*\s*hpsi_sipg", adapter_source, re.I)
assert re.search(r"yn_dg_dc_local_periodic\s*==\s*'y'", main_source, re.I)
assert re.search(r"yn_dg_dc_local_periodic\s*==\s*'y'.*yn_dc_lcfo", main_source, re.I | re.S)
assert "dc_lcfo(" in main_source, "normal DC-LCFO route must remain present"
assert "dc_lcfo_flux(" in main_source, "normal flux route must remain present"

print("PASS DG DC local-periodic route contract")
