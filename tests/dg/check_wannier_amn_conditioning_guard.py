#!/usr/bin/env python3
from pathlib import Path


root = Path(__file__).resolve().parents[2]
src = (root / "src" / "gs" / "dc" / "lcfo_flux.f90").read_text()
global_src = (root / "src" / "io" / "salmon_global.f90").read_text()
io_src = (root / "src" / "io" / "inputoutput.f90").read_text()

for name in ["wannier_amn_svd_tol", "wannier_amn_reject_tol"]:
    if name not in global_src:
        raise SystemExit(f"missing global namelist variable: {name}")
    if name not in io_src:
        raise SystemExit(f"missing input/broadcast/log handling for: {name}")

if "subroutine diagnose_wannier_amn_conditioning" not in src:
    raise SystemExit("missing AMN conditioning diagnostic helper")

if "call diagnose_wannier_amn_conditioning(nband_wann, wannier_num_wann)" not in src:
    raise SystemExit("Wannier AMN conditioning must be checked after writing .amn")

for marker in [
    "[DC-LCFO-WANNIER-AMN]",
    "min_sv_rel=",
    "near_null=",
    "wannier_amn_reject_tol",
    "DC-LCFO Wannier export: AMN projection matrix is rank deficient",
]:
    if marker not in src:
        raise SystemExit(f"missing AMN conditioning marker: {marker}")

print("Wannier AMN conditioning diagnostics are wired into DC-LCFO export")
