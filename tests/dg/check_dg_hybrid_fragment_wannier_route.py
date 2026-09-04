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
match = re.search(r"\bsubroutine\s+" + entry + r"\b(.*?)\bend subroutine\s+" + entry, source, re.S)
assert match, "missing separated divided production routine"
route = match.group(1)
assert "nproc==dc%n_frag" in re.sub(r"\s+", "", route), (
    "divided production entry must require MPI size equal to fragment count"
)
required = (
    "build_dg_hybrid_fragment_wannier_from_dc_seed",
    "extract_dg_hybrid_fragment_self_block",
    "initialize_dg_hybrid_fragment_subspace",
    "advance_dg_hybrid_fragment_epoch",
    "run_dc_fragment_occupation_epoch",
    "measure_dg_hybrid_fragment_core_norms",
    "run_dg_hybrid_divided_scf",
)
for name in required:
    assert re.search(r"\bcall\s+" + name + r"\b", route), f"missing production call: {name}"
assert route.index("call build_dg_hybrid_fragment_wannier_from_dc_seed") < route.index("call run_dg_hybrid_divided_scf")
for name in ("dc_lcfo", "run_dg_overlapping_wannier_ground_state_for_main",
             "setup_dg_w90_gamma_library", "run_dg_w90_gamma_library",
             "apply_dg_hybrid_divided_fragment_hpsi", "solve_dg_hybrid_fragment_spectrum"):
    assert not re.search(r"\bcall\s+" + name + r"\b", route), f"legacy production fallback: {name}"
assert not re.search(r"\b384\b", route), "material-specific state count in production"
print("fragment-local DC-to-Wannier production route: PASS")
