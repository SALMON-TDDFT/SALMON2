#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
global_text = (root / "src/io/salmon_global.f90").read_text()
input_text = (root / "src/io/inputoutput.f90").read_text()
scf_text = (root / "src/gs/scf_iteration.f90").read_text()

assert "nstate_freeze_gs" in global_text
assert "nstate_freeze_gs" in input_text
assert "call comm_bcast(nstate_freeze_gs" in input_text
assert "frozen_rwf" in scf_text and "frozen_zwf" in scf_text
assert "nstate_freeze_gs > 0" in scf_text
assert "yn_subspace_diagonalization == 'y' .and. nstate_freeze_gs == 0" in scf_text
print("PASS GS frozen occupied-orbital input and solve wiring")
