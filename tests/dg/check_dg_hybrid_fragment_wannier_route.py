#!/usr/bin/env python3
"""Task 8 production gate; helpers alone do not establish main-route wiring.

This gate intentionally remains RED until the divided production route is
separated from legacy complete-system construction and wired to the new kernels.
It is a source contract, not a replacement for MPI/physical regression tests.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
source = (ROOT / "src/gs/main_dft.f90").read_text().lower()
source = re.sub(r"!.*", "", source).replace("&", "")
entry = "run_dg_hybrid_divided_ground_state_for_main"
assert re.search(r"if\s*\(\s*yn_dg_hybrid_divided_scf\s*==\s*'y'\s*\)\s*then\s*call\s+" + entry, source), (
    "Task 8 incomplete: divided main must dispatch to its own DC-seed construction entry"
)
assert "divided hybrid scf must dispatch to the schwarz production entry" in source, (
    "the legacy overlapping-Wannier entry must reject accidental divided-SCF routing"
)
match = re.search(r"\bsubroutine\s+" + entry + r"\b(.*?)\bend subroutine\s+" + entry, source, re.S)
assert match, "missing separated divided production routine"
route = match.group(1)
assert "nproc==dc%n_frag" in re.sub(r"\s+", "", route), (
    "divided production entry must require MPI size equal to fragment count"
)
required = (
    "build_dg_hybrid_fragment_wannier_from_dc_seed",
    "select_dg_hybrid_core_wannier",
    "prepare_dg_hybrid_selected_catalog",
    "build_dg_hybrid_projected_local_fragment_basis",
    "prepare_dg_hybrid_selected_trial",
    "export_dg_hybrid_selected_basis_frame",
    "prepare_dg_hybrid_schwarz_candidate_inventory",
    "initialize_dg_hybrid_schwarz_state",
    "build_dg_hybrid_schwarz_schedule",
    "initialize_dg_hybrid_interface_continuation",
    "solve_dg_hybrid_generalized_once_and_publish",
)
for name in required:
    assert re.search(r"\bcall\s+" + name + r"\b", route), f"missing production call: {name}"
assert route.index("call build_dg_hybrid_fragment_wannier_from_dc_seed") < route.index(
    "call initialize_dg_hybrid_interface_continuation"
)
for name in ("dc_lcfo", "run_dg_overlapping_wannier_ground_state_for_main",
             "setup_dg_w90_gamma_library", "run_dg_w90_gamma_library",
             "apply_dg_hybrid_divided_fragment_hpsi", "solve_dg_hybrid_fragment_spectrum"):
    assert not re.search(r"\bcall\s+" + name + r"\b", route), f"legacy production fallback: {name}"
assert not re.search(r"\b384\b", route), "material-specific state count in production"
solver_name = "solve_dg_hybrid_schwarz_fragments"
solver_match = re.search(r"\bsubroutine\s+" + solver_name + r"\b(.*?)\bend subroutine\s+" + solver_name,
                         source, re.S)
assert solver_match, "missing Schwarz fragment solve callback"
solver = solver_match.group(1)
assert solver_name in route, "divided SCF does not receive the Schwarz fragment callback"
assert "advance_dg_hybrid_schwarz_epoch" in solver
assert "assign_dg_hybrid_schwarz_occupations" in solver
assert "dg_hybrid_fragment_cg_steps" in solver, "Schwarz production route ignores the user CG cap"
assert "300d0" in solver, "common-mu occupation epoch omits the 300 K electron temperature"
for obsolete in (
    "solve_dg_hybrid_bounded_fragments",
    "solve_dg_hybrid_divided_fragments",
    "apply_dg_hybrid_divided_fragment_hpsi",
    "apply_dg_hybrid_divided_fragment_metric",
    "assemble_dg_hybrid_divided_core_density",
):
    assert not re.search(r"\bsubroutine\s+" + obsolete + r"\b", source), (
        f"obsolete divided production callback remains: {obsolete}"
    )
print("fragment-local DC-to-Wannier production route: PASS")
