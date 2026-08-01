#!/usr/bin/env python3
"""Production wiring contract for exact buffered-fragment crystallographic symmetry."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
MAIN = (ROOT / "src/gs/main_dft.f90").read_text().lower()
SYMMETRY = (ROOT / "src/gs/dc/dg_overlapping_wannier_symmetry.f90").read_text().lower()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


for routine in (
    "load_sawf_crystallographic_catalog_auto",
    "select_dg_exact_fragment_subgroup",
    "build_dg_fragment_group_representation",
    "promote_dg_exact_global_subgroup",
    "project_dg_fragment_covariant_operators",
):
    require(re.search(rf"\bcall\s+{routine}\b", MAIN) is not None,
            f"production overlapping-Wannier route does not call {routine}")

require("dc%system_tot%rion" in MAIN and "dc%system_tot%kion" in MAIN,
        "fragment symmetry discovery must use instantaneous coordinates and species")
require(re.search(r"fragment.*atom.*(mask|index)", MAIN) is not None,
        "production must construct a separate buffered-fragment atom set")
require(re.search(r"local_symmetry_map\s*\([^)]*,[^)]*\)", MAIN) is not None,
        "production must pass nontrivial local point maps into Wannier construction")
require("point_group_symbol" in MAIN and "space_group_number" in MAIN,
        "production diagnostics must identify each exact fragment group")
require("pre_projection_defect" in MAIN and "post_projection_defect" in MAIN,
        "production evidence must report covariance before and after projection")
require("cross_block_scalar_residual" in MAIN and "cross_block_vector_residual" in MAIN,
        "global promotion must audit scalar and vector cross-fragment blocks")
require(MAIN.index("call promote_dg_exact_global_subgroup") <
        MAIN.index("call project_dg_fragment_covariant_operators"),
        "local fragment symmetry must not project global matrices before promotion")
require(MAIN.index("call project_dg_fragment_covariant_operators") <
        MAIN.index("call compute_dg_overlapping_wannier_matrix_fingerprints"),
        "S/H/X/V must be projected before V3 matrix fingerprints are published")

for forbidden in ("thermal_symmetry_tolerance", "parent_point_group", "force_parent_symmetry"):
    require(forbidden not in MAIN and forbidden not in SYMMETRY,
            f"production must not restore thermally broken parent symmetry: {forbidden}")

print("PASS exact buffered-fragment symmetry is wired into production publication")
