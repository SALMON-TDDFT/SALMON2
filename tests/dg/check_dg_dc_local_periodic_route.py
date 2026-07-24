#!/usr/bin/env python3
import re
from pathlib import Path

root = Path(__file__).resolve().parents[2]
global_source = (root / "src/io/salmon_global.f90").read_text()
input_source = (root / "src/io/inputoutput.f90").read_text()
scf_source = (root / "src/gs/scf_iteration_dft.f90").read_text()
main_source = (root / "src/gs/main_dft.f90").read_text()
dc_source = (root / "src/gs/dc/dcdft.f90").read_text()
handoff_source = (root / "src/gs/dc/dg_dc_handoff.f90").read_text()
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
assert re.search(r"call\s+evaluate_dg_dc_handoff", scf_source, re.I)
assert re.search(r"exit\s+DFT_Iteration", scf_source, re.I)
assert re.search(r"call\s+materialize_dg_dc_candidates", main_source, re.I)
assert "MPI_COMM_SELF" not in main_source, "production materialization must be globally collective"
assert re.search(r"materialize_dg_dc_candidates.*dc%icomm_tot", main_source, re.I | re.S)
assert re.search(
    r"nstate\s*=\s*dg_dc_candidate_orbitals_per_atom\s*\*\s*natom", dc_source, re.I
), "enabled DC solve must create the complete configured untruncated candidate pool"
assert re.search(r"call\s+preserve_dg_dc_density_potential", scf_source, re.I)
assert re.search(r"call\s+discard_dc_mixing_history", scf_source, re.I)
assert "candidate_metric_rank" in handoff_source and "maxloc" in handoff_source
for token in [
    "DG_DC_PRECONVERGENCE",
    "DG_DC_CONTINUATION",
    "DG_DC_FULL_SCF",
    "DG_DC_UNMIXED_FIXED_POINT",
    "DG_DC_ACCEPTED",
    "DG_DC_FAILED",
    "lambda_step=0.5d0*lambda_step",
    "accepted_state=state",
    "accepted_density=density",
    "state=accepted_state",
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
assert re.search(r"call\s+run_dg_dc_ground_state\s*\(", main_source, re.I)
assert re.search(r"call\s+calc_rho_total_dcdft\s*\(", main_source, re.I)
assert re.search(r"call\s+calc_vlocal_fragment_dcdft\s*\(", main_source, re.I)
assert re.search(r"hpsi\s*=\s*hpsi_volume_nonlocal\s*\+\s*lambda\s*\*\s*hpsi_sipg", adapter_source, re.I)
assert re.search(r"yn_dg_dc_local_periodic\s*==\s*'y'", main_source, re.I)
assert re.search(r"yn_dg_dc_local_periodic\s*==\s*'y'.*yn_dc_lcfo", main_source, re.I | re.S)
assert "dc_lcfo(" in main_source, "normal DC-LCFO route must remain present"
assert "dc_lcfo_flux(" in main_source, "normal flux route must remain present"

print("PASS DG DC local-periodic route contract")
