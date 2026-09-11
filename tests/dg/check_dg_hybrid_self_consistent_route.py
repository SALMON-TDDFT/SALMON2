#!/usr/bin/env python3
"""Production route contract for the self-consistent hybrid ground state."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MAIN = (ROOT / "src/gs/main_dft.f90").read_text(errors="replace").lower()
MAIN_COMPACT = "".join(MAIN.replace("&", "").split())
SCF_PATH = ROOT / "tests/dg/legacy_support/dg_hybrid_scf.f90"
STATE_PATH = ROOT / "src/gs/dc/dg_hybrid_ground_state_types.f90"
CHECKPOINT = (ROOT / "src/rt/dg/rt_dg_hybrid_occupied_checkpoint.f90").read_text(
    errors="replace"
).lower()

assert SCF_PATH.is_file(), "missing self-consistent hybrid GS controller"
assert STATE_PATH.is_file(), "missing distributed hybrid occupied-state contract"

scf = SCF_PATH.read_text(errors="replace").lower()
state = STATE_PATH.read_text(errors="replace").lower()

for token in (
    "hybrid_basis_fingerprint",
    "metric_fingerprint",
    "update_potential",
    "assemble_hamiltonian",
    "solve_occupied_states",
    "reconstruct_density",
    "mix_density",
    "density_residual",
    "energy_residual",
    "eigensystem_residual",
    "electron_count_defect",
    "converged",
):
    assert token in scf, f"missing hybrid SCF route contract: {token}"

for token in ("coefficients", "occupations", "eigenvalues", "noccupied"):
    assert token in state, f"missing multi-orbital state field: {token}"

forbidden_inside_loop = (
    "select_dg_hybrid_wannier_blocks",
    "build_dg_hybrid_windowed_pw_catalog",
    "project_dg_hybrid_wannier_complement",
)
for token in forbidden_inside_loop:
    assert token not in scf, f"hybrid SCF must not reselect its fixed basis: {token}"

# The controller is retained as test support, not a selectable production route.
assert "call run_dg_hybrid_self_consistent_ground_state" not in MAIN
inputs = (ROOT / "src/io/inputoutput.f90").read_text().lower()
assert "yn_dg_hybrid_scf is retired" in inputs
assert "call run_dg_hybrid_divided_ground_state_for_main" in MAIN

assert "write_rt_dg_hybrid_occupied_checkpoint" in CHECKPOINT and "coefficients(:,:)" in CHECKPOINT.replace(" ", ""), (
    "hybrid RT checkpoint must store the occupied coefficient matrix"
)
assert "occupations" in CHECKPOINT, "hybrid RT checkpoint must store occupations"

print("self-consistent hybrid GS route contract: PASS")
