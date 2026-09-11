#!/usr/bin/env python3
"""Current DG dispatch; former bare-OW driver assertions are retired.

Standalone numerical module tests remain separate from production wiring.
"""
from pathlib import Path
import runpy

root = Path(__file__).resolve().parents[2]
main = (root / "src/gs/main_dft.f90").read_text().lower()
inputs = (root / "src/io/inputoutput.f90").read_text().lower()
assert "run_dg_overlapping_wannier_ground_state_for_main" not in main
assert "bare overlapping-wannier gs is retired" in inputs
for call in ("call dc_lcfo(", "call dc_lcfo_flux", "call dc_lcfo_soi"):
    assert call in main, f"ordinary DC/LCFO route removed: {call}"
assert "lcfo.f90" in (root / "src/gs/dc/CMakeLists.txt").read_text()
for gate in ("check_retired_hybrid_entries.py", "check_dg_hybrid_terminal_refinement_route.py",
             "check_dg_hybrid_fragment_wannier_route.py"):
    runpy.run_path(str(root / "tests/dg" / gate), run_name="__main__")
print("current DG/ordinary LCFO dispatch contract: PASS")
