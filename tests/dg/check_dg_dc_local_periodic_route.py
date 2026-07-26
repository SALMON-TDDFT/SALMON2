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
checkpoint_source = (root / "src/common/dg_ground_state_checkpoint.f90").read_text()
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
assert re.search(
    r"yn_dg_dc_local_periodic\s*==\s*'y'\s*\.and\.\s*yn_spinorbit\s*==\s*'y'.*?"
    r"sawf_input_fatal",
    input_source,re.I|re.S,
), "real-orbital-only DG direct route must reject spin-orbit input before initialization"
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
assert re.search(
    r"active_core\s*=\s*\.not\.\s*any\s*\(\s*local_core_size\s*<=\s*0\s*\).*?"
    r"nowned\s*=\s*merge\s*\(.*?0\s*,\s*active_core\s*\)",
    direct_projector_prepare.group(1),re.I|re.S,
), "buffer-only r-space ranks must remain in collectives with an empty state directory"
assert "frozen_projector_generation" in direct_cg_operator.group(1), \
    "production direct-DG action must accept the frozen outer-SCF projector generation"
assert re.search(
    r"call\s+apply_dg_dc_frozen_projected_face_plane",direct_cg_operator.group(1),re.I
), "production Hamiltonian must consume reconstructed remote value and normal traces"
assert re.search(
    r"call\s+comm_summation\s*\(\s*owned_buffer_face_values\s*,\s*global_buffer_values",
    direct_cg_operator.group(1),re.I,
), "each frozen action must assemble the complete distributed nstate_frag buffer window"
assert re.search(
    r"call\s+comm_summation\s*\(\s*core_trace_owned\s*,\s*core_trace_values",
    direct_cg_operator.group(1),re.I,
), "normal traces must be assembled from raw core owners"
assert len(re.findall(
    r"if\s*\(\s*all\s*\(\s*il(?:ift)?\s*>=\s*mg%is\s*\).*?"
    r"hpsi%rwf",
    direct_cg_operator.group(1),re.I|re.S,
)) >= 2, "face values and normal lifts must be written only by their raw-grid owner"
assert re.search(r"do\s+iface\s*=\s*1\s*,\s*6",direct_cg_operator.group(1),re.I), \
    "the production Hamiltonian must apply all six signed physical faces"
assert "same-state" not in direct_cg_operator.group(1).lower(), \
    "same-state neighbor diagnostics must not remain in the production Hamiltonian"
assert not re.search(
    r"neighbor_value\s*=\s*cmplx\s*\(\s*buf_recv\(.*?\bilo\b",
    direct_cg_operator.group(1),re.I|re.S,
), "production Hamiltonian must not use the neighbor state with the same local io index"
gscg_rwf_source = re.search(
    r"subroutine\s+gscg_rwf\b(.*?)end\s+subroutine\s+gscg_rwf",
    cg_source,re.I|re.S,
)
assert gscg_rwf_source and not re.search(
    r"build_dg_buffer_window_projector|prepare_dg_dc_buffer_core_projectors",
    gscg_rwf_source.group(1),re.I,
), "the frozen projector must not be rebuilt inside the orbital CG sweep"
assert not re.search(r"call\s+apply_dc_lcfo_flux_hpsi_rwf",direct_cg_operator.group(1),re.I), \
    "the direct SIPG action must not delegate to the LCFO flux approximation"
for component in ["total_value", "total_normal"]:
    assert component in direct_sipg_source, \
        f"direct CG SIPG action is missing {component}"
assert "dg_dc_gs_sipg_penalty_factor" in direct_cg_operator.group(1)
assert "evaluate_dg_nodal_sipg_face" in direct_sipg_source
assert re.search(
    r"use\s+dg_dc_direct_sipg\s*,\s*only\s*:\s*apply_dg_dc_frozen_projected_face_plane",
    direct_cg_operator.group(1),re.I,
), "production CG must use the numerically tested route-neutral SIPG evaluator"
assert re.search(
    r"call\s+apply_dg_dc_frozen_projected_face_plane",direct_cg_operator.group(1),re.I
), "production CG must call the tested SIPG evaluator"
assert re.search(
    r"call\s+calc_eigen_energy.*?call\s+apply_dc_direct_dg_hpsi_rwf"
    r".*?call\s+recompute_direct_dg_eigenvalues",
    scf_source,re.I|re.S,
), "energy evaluation must restore the same direct DG action after conventional hpsi"
assert re.search(
    r"yn_dg_dc_local_periodic\s*==\s*'y'.*?"
    r"comm_logical_and\s*\(\s*system%if_real_orbital\s*\.and\.\s*allocated\s*\(\s*spsi%rwf\s*\)",
    main_source,re.I|re.S,
), "DG direct route must collectively preflight real allocated orbitals before SCF"
direct_diagnostics = re.search(
    r"subroutine\s+diagnose_direct_dg_orbitals\b(.*?)"
    r"end\s+subroutine\s+diagnose_direct_dg_orbitals",
    scf_source,re.I|re.S,
)
assert direct_diagnostics
assert re.search(
    r"do\s+jo_local\s*=\s*1\s*,\s*system%no.*?"
    r"call\s+comm_bcast\s*\(.*?info%icomm_o\s*,\s*info%irank_io\s*\(\s*jo_local\s*\)\s*\)",
    direct_diagnostics.group(1),re.I|re.S,
), "full-scale diagnostics must exchange every state-axis shard before forming global matrices"
assert re.search(
    r"do\s+io_local\s*=\s*info%io_s\s*,\s*info%io_e",
    direct_diagnostics.group(1),re.I,
), "full-scale diagnostics may directly index only the rank-local state interval"
for guard in ["info%ik_s /= 1", "info%im_s /= 1", "psi orbital range mismatch",
              "hpsi orbital range mismatch", "spin range mismatch"]:
    assert guard.lower() in direct_cg_operator.group(1).lower(), \
        f"direct CG decomposition preflight is missing {guard}"
assert re.search(
    r"allocate\s*\(\s*owned_buffer_face_values.*?stat\s*=",
    direct_cg_operator.group(1),re.I|re.S,
), "direct CG projected-face allocation must be checked"
assert re.search(r"call\s+comm_summation\s*\(\s*failure_local\s*,\s*failure_global",
                 direct_cg_operator.group(1),re.I), \
    "rank-local SIPG failures must be agreed before the next halo exchange"
assert re.search(r"use_direct_dg\s*=.*yn_dg_dc_local_periodic\s*==\s*'y'",cg_source,re.I), \
    "the direct DG hook must be independent of LCFO controls"
assert re.search(
    r"if\s*\(\s*use_direct_dg\s*\)\s*use_fixed_flux\s*=\s*\.false\.",
    cg_source,re.I,
), "direct DG must suppress the legacy LCFO-flux diagnostics and operator path"
assert re.search(r"dg_fragment_orbitals_are_finite\s*\(\s*\)", scf_source, re.I)
assert re.search(r"handoff rejected:.*residual\(max\).*charge_error\(max\)",
                 scf_source, re.I | re.S), \
    "below-threshold handoff rejection must expose collective diagnostics"
assert re.search(
    r"if\s*\(\s*dg_handoff_accept\s*\.and\.\s*\.not\.\s*dg_was_accepted\s*\)\s*then"
    r".*?call\s+initialize_dg_dc_direct_continuation"
    r".*?call\s+timer_end\s*\(\s*LOG_WRITE_GS_RESULTS\s*\)",
    scf_source, re.I | re.S,
), "accepted DC handoff must continue through the normal results-timer close"
assert not re.search(r"call\s+materialize_dg_dc_candidates", main_driver_source, re.I), \
    "coefficient DG route must not materialize a global nodal state"
assert not re.search(r"call\s+solve_nodal_ground_state_cg_mpi", main_driver_source, re.I), \
    "production DG route must use distributed local-basis coefficients"
assert re.search(
    r"direct_ground_state_complete.*?call\s+publish_dg_dc_direct_ground_state_for_main",
    main_driver_source,re.I|re.S,
), "direct checkpoint publication must be guarded by full-scale fixed-point completion"
assert re.search(
    r"subroutine\s+publish_dg_dc_direct_ground_state_for_main.*?"
    r"call\s+publish_dg_dc_direct_checkpoint",
    main_source,re.I|re.S,
), "production must publish the direct fragment-orbital checkpoint"
assert re.search(r"checkpoint%nstate\s*=\s*system%no",main_source,re.I), \
    "direct checkpoint must persist the global nstate_frag count"
assert re.search(r"checkpoint%state_start\s*=\s*info%io_s",main_source,re.I), \
    "direct checkpoint must persist the rank-local state-axis origin"
assert re.search(
    r"orbital_core_start\s*=\s*max\s*\(\s*mg%is\s*,\s*\[\s*1\s*,\s*1\s*,\s*1\s*\]\s*\).*?"
    r"checkpoint%orbital_core_origin\s*=\s*checkpoint%fragment_origin",
    main_source,re.I|re.S,
), "direct checkpoint orbital ownership must come from fragment mg, not total-density decomposition"
assert not re.search(
    r"checkpoint%orbital_core_(origin|size)\s*=\s*checkpoint%density_(origin|size)",
    main_source,re.I,
), "fragment orbital ownership must remain independent of total-density ownership"
assert re.search(
    r"checkpoint%occupations\s*=\s*system%rocc\s*\(\s*info%io_s\s*:\s*info%io_e",
    main_source,re.I,
), "direct checkpoint occupations must match the rank-local orbital shard"
assert not re.search(
    r"nstate\s*=\s*dg_dc_candidate_orbitals_per_atom\s*\*", dc_source, re.I
), "DC fragment state count must come from namelist nstate_frag, not a fixed per-atom multiplier"
assert re.search(
    r"dc%nstate_frag\s*=\s*nstate_frag", dc_source, re.I
), "DC fragment state count must preserve the namelist value"
assert not re.search(
    r"(system%no|configured_count)\s*==\s*natom\s*\*\s*"
    r"(dg_dc_candidate_orbitals_per_atom|state%candidate_orbitals_per_atom)",
    main_source + handoff_source, re.I,
), "DG handoff must accept the namelist nstate_frag instead of imposing states-per-atom"
assert re.search(r"call\s+preserve_dg_dc_density_potential", scf_source, re.I)
for rollback_field in ["rho_tot_s", "Vh_tot", "Vxc_tot", "vloc_tot"]:
    assert re.search(
        rf"dg_continuation_rollback.*?dc%{rollback_field}",
        scf_source,re.I|re.S,
    ), f"continuation rollback does not restore {rollback_field}"
mixing_reset = re.search(
    r"subroutine\s+discard_dc_mixing_history\b(.*?)"
    r"end\s+subroutine\s+discard_dc_mixing_history",
    scf_source,re.I|re.S,
)
assert mixing_reset, "DG continuation must define an explicit mixing-history reset"
assert re.search(
    r"subroutine\s+discard_dc_mixing_history\s*\(\s*history\s*,\s*rho_s\s*,\s*Vh\s*,\s*Vxc\s*\)",
    scf_source,re.I,
), "mixing reset must be seeded from the accepted density and potential"
assert not re.search(
    r"history%(rho(_s)?|Vh|Vxc)_(in|out).*?=\s*0d0",
    mixing_reset.group(1),re.I|re.S,
), "mixing reset must not replace the accepted state with a zero-density history"
assert len(re.findall(
    r"call\s+discard_dc_mixing_history\s*\(\s*mixing\s*,\s*dc%rho_tot_s\s*,"
    r"\s*dc%Vh_tot\s*,\s*dc%Vxc_tot\s*\)",
    scf_source,re.I,
)) >= 3, "handoff, stage advance, and rollback must all reseed mixing from the restored state"
assert "direct_dg_solve_count" in handoff_source, \
    "handoff state must record completed direct-DG orbital solves"
assert re.search(
    r"call\s+solve_orbitals.*?interface_scale\s*\).*?"
    r"direct_dg_solve_count\s*=\s*dg_dc_handoff_runtime%direct_dg_solve_count\s*\+\s*1",
    scf_source,re.I|re.S,
), "a nonzero-DG CG pass must be recorded immediately after the orbital solve"
assert re.search(
    r"sum1\s*<\s*threshold.*?direct_ground_state_complete",
    scf_source,re.I|re.S,
), "the DG route must not declare convergence before the full-scale fixed point"
assert re.search(
    r"end\s*do\s+DFT_Iteration.*?yn_dg_dc_local_periodic.*?\.not\.\s*"
    r"dg_dc_handoff_runtime%direct_ground_state_complete.*?stop",
    scf_source,re.I|re.S,
), "nscf exhaustion must not return any DG state before the full-scale fixed point"
assert re.search(
    r"if\s*\(\s*yn_dg_dc_local_periodic\s*==\s*'y'\s*\)\s*then.*?"
    r"direct_ground_state_complete.*?publish_dg_dc_direct_ground_state_for_main.*?"
    r"else\s+if\s*\(\s*yn_dc_lcfo_flux",
    main_source,re.I|re.S,
), "the DG publication branch must not fall through to LCFO or normal DC output"
assert re.search(
    r"dg_continuation_rollback.*?dg_dc_projector_generation\s*="
    r"dg_dc_projector_generation\s*\+\s*1.*?"
    r"dg_dc_frozen_operator_fingerprint\s*=\s*0_8",
    scf_source,re.I|re.S,
), "rollback must explicitly invalidate the rejected frozen projector"
assert re.search(
    r"call\s+collective_and\s*\(\s*unmixed\s*,\s*communicator\s*,\s*collective_unmixed\s*\).*?"
    r"state%interface_scale\s*<\s*1d0.*?else\s+if\s*\(\s*collective_unmixed\s*\)\s*then"
    r".*?state%direct_ground_state_complete\s*=\s*\.true\.",
    handoff_source,re.I|re.S,
), "direct completion must require collective agreement on a full-scale unmixed fixed point"
for diagnostic_gate in [
    "direct_orbital_residual", "direct_orthogonality_defect",
    "direct_hermiticity_defect", "direct_face_balance_defect", "charge_error",
]:
    assert re.search(
        rf"{diagnostic_gate}\s*<=\s*dg_dc_gs_",scf_source,re.I
    ), f"direct completion gate is missing measured {diagnostic_gate}"
assert re.search(
    r"allocate\s*\(\s*density_candidate.*?stat\s*=\s*allocation_status",
    handoff_source,re.I|re.S,
), "handoff snapshot allocation must be checked"
assert re.search(
    r"call\s+collective_and\s*\(\s*local_ok\s*,\s*communicator\s*,\s*ok\s*\)",
    handoff_source,re.I,
), "handoff snapshot allocation failure must be agreed collectively"
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
assert re.search(r"call\s+publish_dg_dc_ground_state_for_main\s*\(",local_basis_driver.group(1),re.I), \
    "accepted local-basis fixed point must publish the DG-only checkpoint"
assert re.search(r"spsi_shape\s*\(\s*1\s*:\s*3\s*\)\s*==\s*sttpsi_shape",
                 local_basis_driver.group(1), re.I), \
    "production adapter must reject mismatched full-buffer orbital shapes"
for persisted_field in ["full_fragment_basis", "basis_transform", "basis_offsets", "fragment_ids"]:
    assert re.search(rf"candidate_state%{persisted_field}", local_basis_driver.group(1), re.I), \
        f"DG-only result must transactionally retain {persisted_field}"
assert re.search(r"iteration_operator_fingerprint\s*=\s*dg_dc_operator_fingerprint\s*\(\s*\)",
                 local_basis_driver.group(1), re.I), \
    "persisted state must identify the Hamiltonian used by the final coefficient solve"
assert re.search(r"hamiltonian_operator_fingerprint\s*=\s*iteration_operator_fingerprint",
                 local_basis_driver.group(1),re.I)
assert re.search(r"candidate_state%operator_fingerprint\s*=\s*dg_dc_operator_fingerprint\s*\(\s*\)",
                 local_basis_driver.group(1),re.I), \
    "checkpoint identity must fingerprint the updated potential actually published"
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
assert re.search(r"yn_dg_dc_local_periodic\s*/=\s*'y'.*?checkpoint_gs",
                 scf_source, re.I | re.S), \
    "in-loop standard checkpoints must be disabled for the opt-in route"
assert re.search(
    r"dg_handoff_count\s*=\s*product\s*\(\s*dc%mg_tot%ie\s*-\s*"
    r"dc%mg_tot%is\s*\+\s*1\s*\)",
    scf_source, re.I,
), "handoff snapshots must allocate only the locally owned grid points"
assert re.search(r"yn_dg_dc_local_periodic\s*==\s*'y'.*?time_shutdown\s*>\s*0",
                 input_source, re.I | re.S), \
    "the no-checkpoint route must reject shutdown-checkpoint mode up front"
assert re.search(
    r"yn_dg_dc_local_periodic\s*==\s*'y'.*?trim\s*\(\s*theory\s*\)\s*/=\s*'dft'",
    input_source,re.I|re.S,
), "the DG ground-state flag must reject conventional RT theories"
assert re.search(
    r"yn_dg_dc_local_periodic\s*==\s*'y'.*?yn_dc_lcfo.*?yn_eigenexa",
    input_source,re.I|re.S,
), "the DG route must reject LCFO, EigenExa, and WPW fallback controls"
assert re.search(r"local_basis_route_active\s*=\s*\.true\.", main_source, re.I), \
    "the DG-only result must mark the standard publication path inactive"
assert re.search(r"if\s*\(\s*\.not\.\s*local_basis_route_active\s*\)\s*then.*?"
                 r"checkpoint_gs", main_source, re.I | re.S), \
    "standard checkpoint publication must be gated off for the DG-only state"
assert "solve_dg_dc_local_basis_bands_reference" not in main_source, \
    "production DG route must not call the replicated root reference eigensolver"
assert not re.search(r"call\s+publish_dg_dc_ground_state_for_main\s*\(", main_driver_source, re.I)
assert re.search(r"type\s*,\s*public\s*::\s*s_dg_dc_direct_checkpoint_state",
                 checkpoint_source,re.I), \
    "DG checkpoint must define a direct fragment-orbital payload"
direct_checkpoint_type = re.search(
    r"type\s*,\s*public\s*::\s*s_dg_dc_direct_checkpoint_state(.*?)end\s+type",
    checkpoint_source,re.I|re.S,
)
assert direct_checkpoint_type
for field in ["fragment_orbitals", "occupations", "nstate", "nspin", "core_size",
              "full_spatial_shape", "geometry_fingerprint", "operator_fingerprint"]:
    assert field in direct_checkpoint_type.group(1), f"direct checkpoint is missing {field}"
assert "coefficient_rows" not in direct_checkpoint_type.group(1), \
    "direct checkpoint must not persist projected coefficient rows"
for field in ["projector_generation", "projector_fingerprint", "projector_retained_rank",
              "projector_projection_residual", "projector_escape_norm",
              "face_action_norm", "face_pair_balance", "accepted", "unmixed_verified"]:
    assert field in direct_checkpoint_type.group(1), f"direct checkpoint is missing {field}"
assert re.search(
    r"checkpoint%geometry_fingerprint\s*=\s*dg_dc_geometry_fingerprint\s*\(\s*\).*?"
    r"checkpoint%operator_fingerprint\s*=\s*dg_dc_operator_fingerprint\s*\(\s*\).*?"
    r"hamiltonian_operator_fingerprint\s*=\s*ieor",
    main_source,re.I|re.S,
), "direct checkpoint must identify the final geometry, potential, and frozen projector"
assert re.search(r"call\s+calc_rho_total_dcdft\s*\(", main_source, re.I)
assert re.search(r"call\s+calc_vlocal_fragment_dcdft\s*\(", main_source, re.I)
assert re.search(r"hpsi\s*=\s*hpsi_volume_nonlocal\s*\+\s*lambda\s*\*\s*hpsi_sipg", adapter_source, re.I)
assert re.search(r"yn_dg_dc_local_periodic\s*==\s*'y'", main_source, re.I)
assert re.search(r"yn_dg_dc_local_periodic\s*==\s*'y'.*yn_dc_lcfo", main_source, re.I | re.S)
assert "dc_lcfo(" in main_source, "normal DC-LCFO route must remain present"
assert "dc_lcfo_flux(" in main_source, "normal flux route must remain present"

print("PASS DG DC local-periodic route contract")
