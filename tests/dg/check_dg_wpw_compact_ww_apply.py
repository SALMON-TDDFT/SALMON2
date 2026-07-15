#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
src = (ROOT / "src/rt/dg/rt_dg_fragment_ops.f90").read_text().lower()

assert "public :: apply_matrix_blocks_batch_compact" in src
assert "subroutine apply_matrix_blocks_batch_compact" in src
body = src.split("subroutine apply_matrix_blocks_batch_compact", 1)[1]
body = body.split("end subroutine apply_matrix_blocks_batch_compact", 1)[0]
for token in ("compact_row_ids", "block_ids", "find_compact_row", "info"):
    assert token in body, f"missing compact WW token: {token}"
for forbidden in ("n_mat_max, n_mat_max", "allocate(h_global", "comm_summation"):
    assert forbidden not in body, f"compact WW apply uses global storage/reduction: {forbidden}"
print("PASS compact WW/face block apply contract")
