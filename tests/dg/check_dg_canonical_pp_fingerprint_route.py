#!/usr/bin/env python3
"""Source contract for canonical PP provenance in saved DC and Hybrid handoff."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
MAIN = (ROOT / "src/gs/main_dft.f90").read_text()
HELPER = (ROOT / "src/gs/dc/dg_canonical_pp_fingerprint.f90").read_text()
CMAKE = (ROOT / "src/gs/dc/CMakeLists.txt").read_text()


def extent(source: str, kind: str, name: str) -> str:
    match = re.search(
        rf"\b{kind}\s+{name}\b(?P<body>.*?)\bend\s+{kind}\s+{name}\b",
        source,
        re.IGNORECASE | re.DOTALL,
    )
    assert match, f"missing {kind} {name}"
    return match.group("body")


assert "dg_canonical_pp_fingerprint.f90" in CMAKE
assert "SALMON-DG-CANONICAL-PP" in HELPER
assert re.search(r"canonical_pp_schema\s*=\s*1_int64", HELPER, re.IGNORECASE)

seed_body = extent(MAIN, "function", "dg_dc_seed_operator_input_fingerprint")
assert "canonical_pp_fingerprint(pp)" in seed_body
assert re.search(
    r"if\s*\(\s*pp_fingerprint\s*==\s*0_int64\s*\)\s*then\s*"
    r"hash\s*=\s*0_int64\s*;?\s*return\s*;?\s*endif",
    seed_body,
    re.IGNORECASE | re.DOTALL,
), "invalid canonical PP input is not handed to the collective seed-contract gate"
assert "pp%zion" not in seed_body.lower()
assert "hash_pp_info_for_dg_dc_seed" not in MAIN.lower()

hybrid_start = MAIN.lower().index("subroutine run_dg_overlapping_wannier_ground_state_for_main")
hybrid_end = MAIN.lower().index("subroutine build_ow_complete_sp_projectors", hybrid_start)
hybrid_body = MAIN[hybrid_start:hybrid_end]
assert "pseudopotential_fingerprint=canonical_pp_fingerprint(pp)" in hybrid_body.replace(" ", "").lower()
for reduction in ("MPI_MIN", "MPI_MAX"):
    assert re.search(
        rf"MPI_Allreduce\s*\(\s*pseudopotential_fingerprint\s*,.*?{reduction}",
        hybrid_body,
        re.IGNORECASE | re.DOTALL,
    ), f"Hybrid PP provenance lacks collective {reduction} agreement"
publisher_body = extent(MAIN, "subroutine", "publish_dg_hybrid_divided_v4")
assert "canonical_pp_valence_sum(pp)" in publisher_body
assert "canonical_pp_digest(pp)" in publisher_body
assert "pp%zion" not in MAIN.lower()

print("PASS canonical PP provenance production route contract")
