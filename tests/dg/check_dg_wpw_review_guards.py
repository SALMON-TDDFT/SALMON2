#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
exchange = (ROOT / "src/rt/dg/rt_dg_wpw_owner_exchange.f90").read_text().lower()
builder = (ROOT / "src/rt/dg/rt_dg_wpw_sparse_builder.f90").read_text().lower()
assert "merge(global_info" not in exchange, "MPI failure reads possibly undefined global_info"
assert "has_duplicate_pairs" not in builder, "sparse candidate validation is quadratic"
assert "pairs_strictly_increasing" in builder
print("PASS review guards avoid undefined MPI status and quadratic sparse validation")
