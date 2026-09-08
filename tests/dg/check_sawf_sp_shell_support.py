#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
sawf = (root / "src/gs/dc/lcfo_wannier_sawf.f90").read_text()
projection = (root / "src/gs/dc/dg_overlapping_wannier_projection.f90").read_text()
flux = (root / "src/gs/dc/lcfo_flux.f90").read_text()

required_sawf = [
    "sawf_projection_shell_lmax",
    "dg_complete_shell_target_count",
    "build_dg_complete_shell_manifest",
    "validate_dg_complete_shell_manifest",
    "max_l",
    "':s;p'",
    "':s;p;d'",
]
missing = [token for token in required_sawf if token not in sawf]
if missing:
    raise SystemExit("SAWF s+p shell support missing: " + ", ".join(missing))
for token in ("case(1);channels_per_atom=4", "case(2);channels_per_atom=9"):
    if token not in projection.replace(" ", ""):
        raise SystemExit("authoritative complete-shell count missing: " + token)

required_flux = [
    "call sawf_projection_shell_lmax",
    "wannier_num_wann",
    "projection_lmax",
]
missing = [token for token in required_flux if token not in flux]
if missing:
    raise SystemExit("lcfo_flux shell selection missing: " + ", ".join(missing))

print("PASS SAWF complete s+p and s+p+d shell selection wiring")
