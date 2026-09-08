#!/usr/bin/env python3
from pathlib import Path

root = Path(".")
lcfo = (root / "src" / "gs" / "dc" / "lcfo_flux.f90").read_text()
inp = (root / "src" / "io" / "inputoutput.f90").read_text()

for token in [
    "direct_amn_bond_block",
    "direct_amn_bond_global",
    "call transform_wannier_position_gauge",
    "basis source=direct_amn_bond_projectors_block",
    "basis source=direct_amn_bond_projectors_global",
    "v_matrix(1:num_bands_chk,1:num_wann_chk) = v_direct",
    "v_matrix(1:num_bands_chk,1:num_wann_chk) = v_direct_global",
]:
    if token not in lcfo:
        raise SystemExit(f"missing lcfo token: {token}")

for token in ["direct_amn_bond_block", "direct_amn_bond_global"]:
    if token not in inp:
        raise SystemExit(f"input validation does not allow {token}")
