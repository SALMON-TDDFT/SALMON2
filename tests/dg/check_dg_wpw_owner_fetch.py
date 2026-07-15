#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
src = (ROOT / "src/rt/dg/rt_dg_wpw_owner_exchange.f90").read_text().lower()

assert "fetch_rows_from_owners" in src
body = src.split("subroutine fetch_rows_from_owners", 1)[1]
body = body.split("end subroutine fetch_rows_from_owners", 1)[0]
for token in ("requested_row_ids", "owned_row_ids", "mpi_alltoallv", "local_info", "ieee_is_finite"):
    assert token in body, f"missing owner-fetch token: {token}"
assert "procedure(row_owner_function)" in body, "owner lookup must be computed, not globally replicated"
for forbidden in ("comm_bcast", "comm_summation", "n_global_basis", "x_global"):
    assert forbidden not in body, f"owner fetch replicates global data: {forbidden}"
assert "row_owner(:)" not in body, "owner fetch requires O(N) replicated owner metadata"
print("PASS owner-targeted arbitrary-vector row fetch contract")
