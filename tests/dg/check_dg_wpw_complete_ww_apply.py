#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ops = (ROOT / "src/rt/dg/rt_dg_fragment_ops.f90").read_text().lower()
adapter = (ROOT / "src/rt/dg/rt_dg_wpw_matrix_free_adapter.f90").read_text().lower()

assert "subroutine apply_complex_matrix_blocks_batch_compact" in ops
assert "public :: apply_complex_matrix_blocks_batch_compact" in ops
hbody = adapter.split("subroutine apply_h_wpw_distributed", 1)[1]
hbody = hbody.split("end subroutine apply_h_wpw_distributed", 1)[0]
for token in ("complex_matrix_block_info", "ww_nl_blocks", "ww_nl_block_ids",
              "ww_real_local_sipg", "ww_complex_nonlocal"):
    assert token in adapter, f"Hamiltonian adapter omits complex/nonlocal WW data: {token}"
assert "apply_complex_matrix_blocks_batch_compact" in adapter
assert "public :: apply_h_wpw_distributed" not in adapter
assert "public :: apply_s_wpw_distributed" not in adapter
print("PASS distributed Hamiltonian includes compact complex/nonlocal WW action")
