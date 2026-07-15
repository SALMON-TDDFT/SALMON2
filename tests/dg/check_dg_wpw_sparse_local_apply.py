#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
src = ROOT / "src/rt/dg/rt_dg_wpw_sparse_blocks.f90"

assert src.exists(), "missing sparse WPW local block algebra"
text = src.read_text().lower()
for token in (
    "type, public :: s_dg_wpw_sparse_blocks",
    "wp_w_row_ids",
    "wp_pw_col_ids",
    "pp_pw_row_ids",
    "pp_pw_col_ids",
    "apply_wp_owned_columns",
    "apply_pp_owned_rows",
    "yw_partial",
    "yp_owned",
):
    assert token in text, f"missing sparse block token: {token}"
for forbidden in ("h_global", "s_global", "matmul("):
    assert forbidden not in text, f"sparse local apply uses dense operation: {forbidden}"
print("PASS sparse WPW local block algebra contract")
