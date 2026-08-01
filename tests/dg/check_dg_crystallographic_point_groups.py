#!/usr/bin/env python3
"""Contract for the complete nonmagnetic crystallographic point-group catalog."""

from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "tests/dg/fixtures/crystallographic_point_groups.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


groups = json.loads(CATALOG.read_text())
require(len(groups) == 32, "catalog must contain all 32 crystallographic point groups")
require([group["number"] for group in groups] == list(range(1, 33)),
        "point-group numbers must be contiguous and deterministic")
require(len({group["symbol"] for group in groups}) == 32,
        "point-group symbols must be unique")
require(all(isinstance(group["order"], int) and 1 <= group["order"] <= 48
            for group in groups), "invalid crystallographic point-group order")
require({group["symbol"] for group in groups} >=
        {"1", "-1", "2/m", "mmm", "4/mmm", "-3m", "6/mmm", "m-3m"},
        "catalog omits a crystallographic crystal system or inversion class")

c_source = (ROOT / "src/gs/dc/lcfo_wannier_spglib.c").read_text()
for token in ("spg_get_dataset", "pointgroup_symbol", "spacegroup_number"):
    require(token in c_source, f"spglib wrapper does not publish {token}")

fortran_source = (ROOT / "src/gs/dc/lcfo_wannier_sawf.f90").read_text().lower()
require(re.search(r"type\s*(?:::\s*)?t_sawf_crystallographic_catalog", fortran_source) is not None,
        "missing normalized crystallographic catalog type")
for token in ("point_group_symbol", "space_group_number", "hall_number",
              "integer_rotation", "fractional_translation"):
    require(token in fortran_source, f"normalized catalog omits {token}")
require("load_sawf_crystallographic_catalog_auto" in fortran_source,
        "missing automatic crystallographic catalog loader")

print("PASS complete 32-point-group crystallographic catalog contract")
