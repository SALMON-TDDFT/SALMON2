#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
htext = (root / "src/rt/dg/rt_dg_fragment_hamiltonian.f90").read_text()
ftext = (root / "src/rt/dg/rt_dg_fragment.f90").read_text()

assert "DC-exported overlap matrix calculated" in htext
assert "allocate(s_local(nrow,nrow)" in ftext
assert "call eigen_zheev(s_dense, seval, svec)" in ftext
assert "horth = matmul(conjg(transpose(sinvhalf)), matmul(h_dense, sinvhalf))" in ftext
assert "evec = matmul(sinvhalf, zorth)" in ftext
print("PASS Full-H seed uses the explicit DG overlap generalized eigenproblem")
