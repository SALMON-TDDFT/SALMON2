#!/usr/bin/env python3
"""Focused source and geometry checks for rank-local SAWF closure diagnostics."""

from collections import Counter
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
FLUX = ROOT / "src/gs/dc/lcfo_flux.f90"
BAND = ROOT / "src/gs/dc/lcfo_wannier_sawf_band.f90"


def owner_2x2x2(point):
    return tuple(0 if coordinate < 16 else 1 for coordinate in point)


def histogram(transform):
    counts = Counter()
    for z in range(16):
        for y in range(16):
            for x in range(16):
                counts[owner_2x2x2(transform((x, y, z)))] += 1
    return sorted(counts.values())


identity_histogram = histogram(lambda point: point)
assert identity_histogram == [4096], identity_histogram

# W=-I, tau=1/8 on a 32^3 mesh gives i' = 4-i (mod 32), in zero-based indices.
inversion_histogram = histogram(
    lambda point: tuple((4 - coordinate) % 32 for coordinate in point)
)
assert inversion_histogram == sorted([125] + [275] * 3 + [605] * 3 + [1331]), (
    inversion_histogram
)

# Mathematical oracle for the logged quantity: ||(I-BB^T)G||_F^2/||G||_F^2.
# An identity-mapped complete target basis closes exactly, while dropping one
# target direction gives finite leakage.
def projected_leakage(transformed, target_basis):
    # Columns are basis functions and these fixtures are orthonormal.
    coefficient = [
        [sum(target_basis[r][i] * transformed[r][j] for r in range(len(transformed)))
         for j in range(len(transformed[0]))]
        for i in range(len(target_basis[0]))
    ]
    residual2 = 0.0
    norm2 = 0.0
    for r, row in enumerate(transformed):
        for j, value in enumerate(row):
            projection = sum(target_basis[r][i] * coefficient[i][j]
                             for i in range(len(target_basis[0])))
            residual2 += (value - projection) ** 2
            norm2 += value ** 2
    return residual2 / norm2


eye2 = [[1.0, 0.0], [0.0, 1.0]]
assert projected_leakage(eye2, eye2) < 1.0e-15
assert projected_leakage(eye2, [[1.0], [0.0]]) > 0.49

source = FLUX.read_text()
entry = source.split("type :: t_sawf_fragment_state_entry", 1)[1].split(
    "end type t_sawf_fragment_state_entry", 1
)[0]
assert re.search(r"real\(8\),\s*allocatable\s*::\s*basis\s*\(:,:\)", entry, re.I), (
    "rank-local closure diagnostic must retain the raw LCFO basis in the bounded cache"
)

loader = source.split("subroutine load_sawf_fragment_state_entry", 1)[1].split(
    "end subroutine load_sawf_fragment_state_entry", 1
)[0]
assert re.search(r"entry%basis\s*=", loader, re.I), "raw basis is not retained"

clearer = source.split("subroutine clear_sawf_fragment_state_entry", 1)[1].split(
    "end subroutine clear_sawf_fragment_state_entry", 1
)[0]
assert "deallocate(entry%basis)" in clearer.lower(), "cached raw basis is not released"

builder = source.split("subroutine build_sawf_dmn_operation_fragment_local", 1)[1].split(
    "end subroutine build_sawf_dmn_operation_fragment_local", 1
)[0]
assert "diagnose_sawf_fragment_basis_closure" in builder
assert builder.count("[DC-LCFO-SAWF-CLOSURE-LOCAL]") == 1, (
    "expected one summary write site per operation/source fragment"
)
assert re.search(r"histogram=", builder, re.I)
assert re.search(r"aggregate_leakage=", builder, re.I)
assert re.search(r"max_block_leakage=", builder, re.I)
assert re.search(r"histogram_entry_cap\s*=\s*\d+", builder, re.I)
assert re.search(r"truncated:", builder, re.I)
assert "target_histogram(target_frag)>0" in builder.replace(" ", "").lower()

band_source = BAND.read_text()
diagnostic = band_source.split("subroutine diagnose_sawf_fragment_basis_closure", 1)[1].split(
    "end subroutine diagnose_sawf_fragment_basis_closure", 1
)[0]
assert "matmul(transpose(target_basis)" in diagnostic.lower()
assert "transformed_basis" in diagnostic.lower()
assert "residual" in diagnostic.lower()
assert re.search(r"residual_norm2\s*=", diagnostic, re.I)
assert re.search(r"transformed_norm2\s*=", diagnostic, re.I)

# Logging must happen after the grid map loops, never once per grid point.
point_loops_end = builder.find("nmapped=0")
log_site = builder.find("[DC-LCFO-SAWF-CLOSURE-LOCAL]")
assert point_loops_end >= 0 and log_site > point_loops_end

print("PASS SAWF rank-local fragment closure diagnostic source and geometry")
