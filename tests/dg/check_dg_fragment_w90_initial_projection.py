#!/usr/bin/env python3
"""Regression gate for fragment-local Wannier90 initial-projection routing."""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
main = (ROOT / "src/gs/main_dft.f90").read_text().lower()
builder = (ROOT / "src/gs/dc/dg_hybrid_fragment_wannier.f90").read_text().lower()

entry = "run_dg_hybrid_divided_ground_state_for_main"
route_match = re.search(
    rf"\bsubroutine\s+{entry}\b(.*?)\bend subroutine\s+{entry}", main, re.S
)
assert route_match, "missing divided Hybrid production entry"
route = route_match.group(1)
fragment_call = re.search(
    r"call\s+build_dg_hybrid_fragment_wannier_from_dc_seed\s*\((.*?)\)",
    route,
    re.S,
)
assert fragment_call and "dg_ow_w90_initial_projection" in fragment_call.group(1), (
    "fragment Wannier90 construction ignores the user-selected initial projection"
)

construct = builder.split("subroutine construct_fragment_wannier", 1)[1].split(
    "end subroutine construct_fragment_wannier", 1
)[0]
compact_construct = re.sub(r"\s+|&", "", construct)
assert "num_iter,trim(initial_projection),dg_w90_unconstrained" in compact_construct, (
    "fragment Wannier90 setup hard-codes its initial projection"
)

contract = builder.split("subroutine validate_fragment_contract", 1)[1].split(
    "end subroutine validate_fragment_contract", 1
)[0]
compact_contract = re.sub(r"\s+|&", "", contract)
assert "hash_character(replica_hash,trim(initial_projection))" in compact_contract
assert "hash_character(basis_fingerprint,trim(initial_projection))" in compact_contract

print("fragment Wannier90 initial-projection route: PASS")
