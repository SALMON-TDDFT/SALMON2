#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
common_state = root / "src/common/dg_nodal_state.f90"
common_interfaces = root / "src/common/dg_nodal_interfaces.f90"
sipg = root / "src/common/dg_nodal_sipg.f90"
rt_types = root / "src/rt/dg/rt_dg_nodal_types.f90"
solver = root / "src/rt/dg/rt_dg_nodal_ground_state_solver.f90"
cg = root / "src/rt/dg/rt_dg_nodal_cg.f90"
density = root / "src/rt/dg/rt_dg_nodal_density.f90"

assert common_state.exists(), "missing GS/RT-neutral nodal state"
assert common_interfaces.exists(), "missing GS/RT-neutral nodal callback interfaces"
assert sipg.exists(), "missing complete Hermitian nodal SIPG face operator"

state_body = common_state.read_text().lower()
for token in (
    "type, public :: s_dg_nodal_common_state",
    "psi_core",
    "occupations",
    "core_owner",
    "geometry_fingerprint",
    "operator_fingerprint",
    "initialize_dg_nodal_common_state",
    "release_dg_nodal_common_state",
    "validate_dg_nodal_common_state_mpi",
):
    assert token in state_body, f"missing common-state contract: {token}"

interfaces_body = common_interfaces.read_text().lower()
for token in (
    "nodal_complete_h_action",
    "nodal_density_builder",
    "nodal_subspace_rotation",
    "nodal_collective_validator",
):
    assert token in interfaces_body, f"missing neutral callback: {token}"

rt_body = rt_types.read_text().lower()
assert "use dg_nodal_state" in rt_body
assert "type :: s_dg_nodal_state" not in rt_body, "RT wrapper must not own a second nodal state"
assert "dg_ground_state_ready" not in state_body, "readiness must have one source of truth"

for path in (solver, cg, density):
    body = path.read_text().lower()
    assert "use dg_nodal_interfaces" in body, f"{path.name} must consume neutral callback declarations"

common_cmake = (root / "src/common/CMakeLists.txt").read_text()
assert "dg_nodal_state.f90" in common_cmake
assert "dg_nodal_interfaces.f90" in common_cmake
assert "dg_nodal_sipg.f90" in common_cmake

sipg_body = sipg.read_text().lower()
for token in (
    "consistency_value",
    "symmetry_normal",
    "penalty_value",
    "hermiticity_defect",
    "internal_cancellation_defect",
    "jump_norm",
    "penalty_energy",
    "trace_epoch",
    "canonical_owner",
    "physical_face",
    "auxiliary_periodic_wrap",
):
    assert token in sipg_body, f"missing SIPG contract: {token}"

print("PASS nodal ground-state route depends on the GS/RT-neutral common core")
