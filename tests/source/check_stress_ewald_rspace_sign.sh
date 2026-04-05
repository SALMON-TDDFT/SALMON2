#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
stress_f90="$repo_root/src/common/stress.f90"

python3 - <<'PY' "$stress_f90"
from pathlib import Path
import re
import sys

text = Path(sys.argv[1]).read_text()
match = re.search(
    r"subroutine calc_stress_ewa\b(.*?)end subroutine calc_stress_ewa",
    text,
    re.S,
)
if not match:
    raise SystemExit("calc_stress_ewa block not found")

block = match.group(1)

required = [
    "E_ewa_R_loc = E_ewa_R_loc + 0.5d0 * pp%zps(kion(ia)) * pp%zps(kion(ib))",
    "strs_R(a,b) = strs_R(a,b) + fact * rab(a) * rab(b)",
]
for needle in required:
    if needle not in block:
        raise SystemExit(f"missing expected statement: {needle}")

expected = "fact = -0.5d0 * pp%zps(kion(ia)) * pp%zps(kion(ib)) / (V * r_abs**3)"
if expected not in block:
    raise SystemExit(f"expected Ewald R-space stress sign not found: {expected}")
PY
