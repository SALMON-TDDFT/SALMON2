#!/usr/bin/env python3
"""Production route contract for the self-consistent hybrid ground state."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MAIN = (ROOT / "src/gs/main_dft.f90").read_text(errors="replace").lower()
SCF_PATH = ROOT / "src/gs/dc/dg_hybrid_scf.f90"
STATE_PATH = ROOT / "src/gs/dc/dg_hybrid_ground_state_types.f90"
CHECKPOINT = (ROOT / "src/rt/dg/rt_dg_hybrid_checkpoint.f90").read_text(
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

main_tokens = (
    "yn_dg_hybrid_scf=='y'",
    "call run_dg_hybrid_self_consistent_ground_state",
    "hybrid_ground_state%converged=.true.",
    "call write_rt_dg_hybrid_occupied_checkpoint",
)
positions = [MAIN.find(token) for token in main_tokens]
assert all(position >= 0 for position in positions), (
    "production must converge the hybrid GS before checkpoint publication"
)
assert positions == sorted(positions), "hybrid GS/checkpoint production order is invalid"

for token in (
    "ow_hybrid_density=ow_initial_occupied_density",
    "if(ok.and.reusable.and.yn_dg_hybrid_scf/='y')then",
    "min(dg_dc_gs_electron_count_tolerance,dg_ow_symmetry_tolerance)",
    "call mix_dg_overlapping_wannier_density_history",
    "ow_hybrid_new_history=ow_hybrid_density_history;ow_hybrid_history_count=0",
    "update_auxiliary_pencil=.false.",
    "hybrid_metric_fingerprint",
    "spectral_catalog_fingerprint",
    "band_energy",
):
    assert token in MAIN, f"production hybrid SCF is missing required safety contract: {token}"

assert "if(callback_ok.and.size(output_density)==size(density))" not in MAIN, (
    "density callback must not rely on Fortran short-circuit evaluation"
)
assert "call pulay(dc%mg_tot" not in MAIN, (
    "hybrid SCF must use the previously validated bounded two-point Anderson mixer, not full Pulay"
)

assert "write_rt_dg_hybrid_occupied_checkpoint" in CHECKPOINT and "coefficients(:,:)" in CHECKPOINT.replace(" ", ""), (
    "hybrid RT checkpoint must store the occupied coefficient matrix"
)
assert "occupations" in CHECKPOINT, "hybrid RT checkpoint must store occupations"

print("self-consistent hybrid GS route contract: PASS")
