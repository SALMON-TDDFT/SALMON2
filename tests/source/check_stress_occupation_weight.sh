#!/bin/sh
set -eu

stress_f90="$(CDPATH= cd -- "$(dirname "$0")/../../src/common" && pwd)/stress.f90"

python3 - <<'PY' "$stress_f90"
import pathlib
import re
import sys

text = pathlib.Path(sys.argv[1]).read_text()

for name in ("calc_stress_kin", "calc_stress_nl"):
    match = re.search(rf"subroutine {name}\(.*?end subroutine {name}", text, re.S)
    if match is None:
        raise SystemExit(1)
    body = match.group(0)
    if "rtmp = 2d0 * system%rocc" in body:
        raise SystemExit(1)
    if "rtmp = system%rocc(io,ik,ispin) * system%wtk(ik) * system%Hvol" not in body:
        raise SystemExit(1)
PY
