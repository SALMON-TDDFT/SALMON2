#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
src = ROOT / "src/rt/dg/rt_dg_wpw_column_layout.f90"

assert src.exists(), "missing production (K,G) column layout"
text = src.read_text().lower()
for token in (
    "type, public :: s_dg_wpw_column_layout",
    "basis_kind",
    "windowed_kg",
    "n_global_columns",
    "n_g_modes",
    "pw_fragment_ids",
    "pw_g_ids",
    "pw_owner",
    "owned_column_ids",
    "initialize_wpw_column_layout",
    "wpw_column_id",
    "wpw_column_pair",
):
    assert token in text, f"missing column-layout token: {token}"
for forbidden in ("n_plane_waves", "k_pw", "s_dg_fragment_rt"):
    assert forbidden not in text, f"production layout is coupled to legacy G-only state: {forbidden}"
print("PASS production (K,G) column-layout contract")
