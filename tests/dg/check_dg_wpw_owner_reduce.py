#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
src = ROOT / "src/rt/dg/rt_dg_wpw_owner_exchange.f90"

assert src.exists(), "missing owner-targeted WP partial-sum exchange"
text = src.read_text().lower()
for token in (
    "reduce_w_partial_to_owners",
    "mpi_alltoallv",
    "mpi_allreduce",
    "row_owner_function",
    "owned_row_ids",
    "local_info",
    "ieee_is_finite",
):
    assert token in text, f"missing owner exchange token: {token}"
for forbidden in ("comm_summation", "n_global_basis", "yw_global"):
    assert forbidden not in text, f"owner exchange uses global reduction/storage: {forbidden}"
assert "row_owner(:)" not in text, "owner reduction requires O(N) replicated owner metadata"
print("PASS owner-targeted WP partial-sum exchange contract")
