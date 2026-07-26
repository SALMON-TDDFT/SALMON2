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
direct_sipg_source = (root / "src/common/dg_dc_direct_sipg.f90").read_text()
scf_driver_source = (root / "src/gs/scf_iteration.f90").read_text()
dc_source = (root / "src/gs/dc/dcdft.f90").read_text()
handoff_source = (root / "src/gs/dc/dg_dc_handoff.f90").read_text()
buffer_projector_source = (root / "src/common/dg_buffer_window_projector.f90").read_text()
buffer_core_faces_source = (root / "src/gs/dc/dg_dc_buffer_core_faces.f90").read_text()
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
assert "dg_dc_buffer_core_faces.f90" in cmake_source
for face_contract in [
    "physical_grid_ids",
    "local_buffer_indices",
    "neighbor_core_indices",
    "neighbor_core_values",
    "neighbor_core_normals",
    "configured_states",
    "MPI_Allreduce",
    "project_dg_dc_buffer_core_face",
]:
    assert face_contract.lower() in buffer_core_faces_source.lower(), \
        f"buffer/core face mapping is missing {face_contract}"
assert re.search(r"allocate\s*\(\s*faces\s*\(\s*6\s*\)",buffer_core_faces_source,re.I), \
    "DG mapping must construct exactly six signed physical faces"
assert "27" not in buffer_core_faces_source, \
    "edge/corner DC neighbor metadata must not become DG faces"
assert re.search(
    r"candidate_normals\s*=\s*matmul\s*\(\s*face%neighbor_core_normals\s*,\s*"
    r"candidate_coefficients\s*\)",
    buffer_core_faces_source,re.I,
), "value and normal traces must use the same frozen projector coefficients"
for preparation_call in [
    "build_dg_dc_buffer_core_faces",
    "validate_dg_dc_buffer_core_faces",
    "exchange_dg_dc_buffer_core_state_window",
    "project_dg_dc_buffer_core_face",
]:
    assert re.search(rf"call\s+{preparation_call}\s*\(", adapter_source, re.I), \
        f"ground-state adapter must actually call {preparation_call}"
assert re.search(
    r"subroutine\s+prepare_dg_dc_buffer_core_projectors\b",
    adapter_source,
    re.I,
), "outer-SCF projector preparation entry point is missing"
assert re.search(r"call\s+evaluate_dg_dc_handoff", scf_source, re.I)
assert re.search(
    r"call\s+solve_orbitals\s*\(.*?\bdc\s*,\s*dg_dc_handoff_runtime%interface_scale\s*\)",
    scf_source,re.I|re.S,
), "the existing DC orbital solve must receive the current DG scale"
assert re.search(
    r"call\s+gscg_rwf\s*\(.*?\bdc\s*,\s*dg_scale\s*\)",
    scf_driver_source,re.I|re.S,
), "the existing orthogonalized CG must own the direct DG action"
assert re.search(
    r"direct_dg_active\s*=.*present\s*\(\s*dg_scale\s*\).*dg_scale\s*>\s*0",
    scf_driver_source,re.I|re.S,
), "the orbital solver must identify a nonzero direct-DG Hamiltonian"
assert re.search(
    r"if\s*\(\s*yn_subspace_diagonalization\s*==\s*'y'.*?"
    r"\.not\.\s*direct_dg_active\s*\)\s*then",
    scf_driver_source,re.I|re.S,
), "volume-only subspace diagonalization must be disabled for a nonzero DG Hamiltonian"
assert len(re.findall(r"call\s+apply_dc_direct_dg_hpsi_rwf\s*\(",cg_source,re.I))>=3, \
    "the DG action must be recomputed for xk, pk, and refreshed xk"
direct_cg_operator = re.search(
    r"subroutine\s+apply_dc_direct_dg_hpsi_rwf\b(.*?)"
    r"end\s+subroutine\s+apply_dc_direct_dg_hpsi_rwf",
    cg_source,re.I|re.S,
)
assert direct_cg_operator
direct_projector_prepare = re.search(
    r"subroutine\s+prepare_dc_direct_dg_projectors_rwf\b(.*?)"
    r"end\s+subroutine\s+prepare_dc_direct_dg_projectors_rwf",
    cg_source,re.I|re.S,
)
assert direct_projector_prepare
assert re.search(
    r"call\s+prepare_dg_dc_buffer_core_projectors",direct_projector_prepare.group(1),re.I
), "outer-SCF freeze must use the validated Task 2 buffer/core projector API"
assert re.search(
    r"call\s+build_dg_dc_buffer_core_faces",direct_projector_prepare.group(1),re.I
), "production projector preparation must construct canonical six-face metadata"
assert "frozen_projector_generation" in direct_cg_operator.group(1)
assert re.search(
    r"call\s+apply_dg_dc_frozen_projected_face_plane",direct_cg_operator.group(1),re.I
), "production Hamiltonian must consume reconstructed remote value and normal traces"
assert re.search(
    r"call\s+comm_summation\s*\(\s*owned_buffer_face_values\s*,\s*global_buffer_values",
    direct_cg_operator.group(1),re.I,
), "each frozen action must assemble the complete distributed state window"
assert re.search(
    r"call\s+comm_summation\s*\(\s*core_trace_owned\s*,\s*core_trace_values",
    direct_cg_operator.group(1),re.I,
), "normal traces must be assembled from raw core owners"
assert re.search(r"do\s+iface\s*=\s*1\s*,\s*6",direct_cg_operator.group(1),re.I)
assert "same-state" not in direct_cg_operator.group(1).lower()
assert not re.search(
    r"neighbor_value\s*=\s*cmplx\s*\(\s*buf_recv\(.*?\bilo\b",
    direct_cg_operator.group(1),re.I|re.S,
), "production Hamiltonian must not use the same local state index"
gscg_rwf_source = re.search(
    r"subroutine\s+gscg_rwf\b(.*?)end\s+subroutine\s+gscg_rwf",
    cg_source,re.I|re.S,
)
assert gscg_rwf_source and not re.search(
    r"build_dg_buffer_window_projector|prepare_dg_dc_buffer_core_projectors",
    gscg_rwf_source.group(1),re.I,
), "the frozen projector must not be rebuilt inside the orbital CG sweep"
assert not re.search(r"call\s+apply_dc_lcfo_flux_hpsi_rwf",direct_cg_operator.group(1),re.I)
assert "total_value" in direct_sipg_source and "total_normal" in direct_sipg_source
assert re.search(
    r"use\s+dg_dc_direct_sipg\s*,\s*only\s*:\s*apply_dg_dc_frozen_projected_face_plane",
    direct_cg_operator.group(1),re.I,
)
assert re.search(
    r"call\s+calc_eigen_energy.*?call\s+apply_dc_direct_dg_hpsi_rwf"
    r".*?call\s+recompute_direct_dg_eigenvalues",
    scf_source,re.I|re.S,
), "energy evaluation must restore the same direct DG action"
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
assert not re.search(r"call\s+run_dg_dc_local_basis_ground_state_for_main\s*\(", main_driver_source, re.I), \
    "production DG route must stay in the existing DC orbital CG"
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
