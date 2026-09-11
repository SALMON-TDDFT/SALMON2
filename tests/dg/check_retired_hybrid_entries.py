"""Retired Hybrid selectors must never enter legacy production drivers."""
from pathlib import Path

root = Path(__file__).resolve().parents[2]
gs = (root / "src/gs/main_dft.f90").read_text().lower()
rt = (root / "src/rt/main_tddft.f90").read_text().lower()
inputs = (root / "src/io/inputoutput.f90").read_text().lower()
assert "call run_dg_hybrid_self_consistent_ground_state" not in gs
assert "subroutine run_dg_overlapping_wannier_ground_state_for_main" not in gs
assert "call run_dg_overlapping_wannier_ground_state_for_main" not in gs
assert "bare overlapping-wannier gs is retired" in inputs
for name in ("populate_ow_checkpoint", "prepare_ow_fixed_center_group", "ow_build_hamiltonian",
             "write_ow_ground_state_evidence", "restore_ow_checkpoint_density"):
    assert f"subroutine {name}" not in gs, f"unreferenced bare-OW helper remains: {name}"
assert "run_dg_overlapping_wannier_coefficient_rt" not in rt
assert "read_dg_overlapping_wannier_checkpoint" not in rt
for name in ("ow_hybrid_update_potential", "ow_hybrid_assemble_hamiltonian",
             "ow_hybrid_solve_occupied", "ow_hybrid_reconstruct_density",
             "ow_hybrid_basis_provider", "ow_hybrid_density_mix"):
    assert name not in gs, f"unused old SCF callback remains: {name}"
assert "dg_hybrid_scf.f90" not in (root / "src/gs/dc/CMakeLists.txt").read_text()
assert "rt_dg_overlapping_wannier.f90" not in (root / "src/rt/CMakeLists.txt").read_text()
assert "use dg_hybrid_continuation_controller" not in gs
assert "use dg_hybrid_publication_policy" in gs
assert "dg_hybrid_continuation_controller.f90" not in (root / "src/gs/dc/CMakeLists.txt").read_text()
assert "build_dg_hybrid_scope_receipt" not in gs
for name in ("divided_lcfo_hrows", "divided_production_selection", "hybrid_scf_fingerprint",
             "dg_certified_localizer_values", "divided_full_basis_closure_ok",
             "local_occupied_values", "hybrid_converged_density", "spectral_complement_row_ids",
             "hybrid_localization_seed_fingerprint"):
    assert name not in gs, f"unused retired DG workspace remains: {name}"
for retired_module in ("dg_hybrid_divided_scf.f90", "dg_hybrid_real_space_residual.f90",
                       "dg_hybrid_continuation_state.f90"):
    assert retired_module not in (root / "src/gs/dc/CMakeLists.txt").read_text()
for flag in ("yn_dg_hybrid_scf", "yn_dg_overlapping_wannier_rt", "yn_dg_overlapping_wannier_rt_restart"):
    assert f"{flag} is retired" in inputs, f"missing retirement error for {flag}"
assert "call run_dg_hybrid_divided_ground_state_for_main" in gs
assert "call run_dg_hybrid_continuation_rt" in rt
for name in ("dg_hybrid_scf", "dg_hybrid_divided_scf", "dg_hybrid_real_space_residual",
             "dg_hybrid_continuation_state", "dg_hybrid_continuation_controller"):
    assert not (root / f"src/gs/dc/{name}.f90").exists()
    assert (root / f"tests/dg/legacy_support/{name}.f90").is_file()
assert not (root / "src/rt/dg/rt_dg_overlapping_wannier.f90").exists()
assert (root / "tests/dg/legacy_support/rt_dg_overlapping_wannier.f90").is_file()
print("retired Hybrid entry contract: PASS")
