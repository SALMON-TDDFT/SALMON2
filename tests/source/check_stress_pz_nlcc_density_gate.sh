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
    r"subroutine calc_stress_xc_builtin_pz\b(.*?)end subroutine calc_stress_xc_builtin_pz",
    text,
    re.S,
)
if not match:
    raise SystemExit("calc_stress_xc_builtin_pz block not found")

block = match.group(1)

required = [
    "rho_xc = rho_s(ispin)%f(ix,iy,iz)",
    "trho_xc = rho_xc",
    "trho_xc = trho_xc + ppn%rho_nlcc(ix,iy,iz)",
    "E_vx_loc = E_vx_loc + rho_xc * vxc_x",
    "E_vc_loc = E_vc_loc + rho_xc * vxc_c",
]
for needle in required:
    if needle not in block:
        raise SystemExit(f"missing expected statement: {needle}")

forbidden = [
    "rho_xc = rho_xc + 0.5d0 * ppn%rho_nlcc(ix,iy,iz)",
]
for needle in forbidden:
    if needle in block:
        raise SystemExit(f"forbidden statement present: {needle}")
PY
