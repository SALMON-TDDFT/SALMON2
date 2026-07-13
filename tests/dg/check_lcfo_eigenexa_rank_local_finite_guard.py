#!/usr/bin/env python3
from pathlib import Path


source = (Path(__file__).resolve().parents[2] / "src/gs/dc/lcfo_flux.f90").read_text()
diag = source.split('if(dc%id_tot==0) write(*,*) "h_div: done"', 1)[1].split(
    "call eigen_sx(n, n, h_div", 1
)[0]

required = [
    "ieee_is_finite(h_div)",
    "DC-LCFO-EIGENEXA-LOCAL-H",
    "rank=",
    "local_i=",
    "local_j=",
]
missing = [token for token in required if token not in diag]
if missing:
    raise SystemExit(
        "rank-local EigenExa Hamiltonian finite guard missing: " + ", ".join(missing)
    )

print("PASS LCFO EigenExa rejects non-finite rank-local Hamiltonian blocks")
