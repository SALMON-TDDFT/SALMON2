#!/usr/bin/env python3
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
source = (ROOT / "src/rt/dg/rt_dg_wpw_candidate_halo.f90").read_text()

assert "MPI_Wtime" in source, "candidate route must measure wall-clock time"
assert "[DG-WPW-ROUTE]" in source, "candidate route must publish a production timing record"
for field in ("staged_local=", "routed_local=", "wp_owned=", "pp_owned=", "seconds="):
    assert field in source, f"candidate route timing record is missing {field}"

print("PASS WPW candidate route exposes bounded performance diagnostics")
