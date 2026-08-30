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
    "assemble_dg_hybrid_production_interface_rows",
]:
    assert token in source, f"production continuation does not connect {token}"

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

driver_name = "subroutine run_dg_hybrid_concrete_continuation"
assert driver_name in source, "missing contained concrete continuation driver"
start = source.index(driver_name)
implementation = source[start:source.index("end subroutine run_dg_hybrid_concrete_continuation", start)]

required = [
    "dc_seed_density",
    "fixed_payload",
    "production_faces",
    "dg_dc_update_potential_from_density",
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
    "nstate+1",
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
assert "occupied_unoccupied_gap>dg_dc_gs_final_orbital_tolerance" in implementation
assert "hamiltonian_hermiticity<=dg_dc_gs_hermiticity_tolerance" in implementation

print("PASS concrete DG continuation production route contract")
