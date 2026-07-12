#!/usr/bin/env python3
from pathlib import Path


source = (Path(__file__).resolve().parents[2] / "src/rt/main_tddft.f90").read_text()
needle = "'[DG-DCDFT-SEED-HRES-BPW-SCF]', full_h_seed_hres"
assert needle in source
assert "Projected Wannier seed is not an eigenstate of the full DG Hamiltonian" in source
assert "stop 'DG-Fragment RT: projected Wannier seed is not stationary'" in source
print("PASS projected Wannier RT seed has a full-DG stationarity gate")
