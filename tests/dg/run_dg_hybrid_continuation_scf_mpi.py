#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
obsolete = root / "src/gs/dc/dg_hybrid_continuation_scf.f90"
assert not obsolete.exists(), "obsolete callback continuation module remains"

source = (root / "src/gs/main_dft.f90").read_text().lower()
assert "use dg_hybrid_continuation_scf" not in source, (
    "production still imports the obsolete callback continuation module"
)
assert "run_dg_hybrid_continuation_scf_fixture" not in source, (
    "production still calls the obsolete fixture-only continuation entry point"
)
print("PASS obsolete callback continuation layer removed")
