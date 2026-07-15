#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
src = (ROOT / "src/rt/dg/rt_dg_wpw_owner_exchange.f90").read_text().lower()

for token in ("iso_fortran_env", "int64", "validate_exchange_sizes", "huge(0)",
              "collective_validation_handshake"):
    assert token in src, f"missing collective MPI count overflow guard: {token}"
assert src.count("call validate_exchange_sizes") >= 2
for token in ("mpi_abort", "abort_mpi_collective", "mpi collective failed"):
    assert token in src, f"MPI collective failure is not communicator-fatal: {token}"
assert src.count("call abort_mpi_collective") >= 10, "some owner-exchange collective failures still return locally"
print("PASS owner exchange preflights MPI integer counts collectively")
