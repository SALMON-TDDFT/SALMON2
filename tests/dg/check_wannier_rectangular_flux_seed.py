#!/usr/bin/env python3
from pathlib import Path


root = Path(__file__).resolve().parents[2]
src = (root / "src" / "gs" / "dc" / "lcfo_flux.f90").read_text()

required = [
    "subroutine build_wannier_flux_eigen_seed_from_transform",
    "h_wann(iw,jw) = h_wann(iw,jw) + conjg(v_matrix(ib,iw)) * esp_in(ib,ispin) * v_matrix(ib,jw)",
    "call eigen_zheev(h_wann, eval, eigvec)",
    "[DC-LCFO-W90-SEED] rectangular transform projected and diagonalized:",
]

for marker in required:
    if marker not in src:
        raise SystemExit(f"missing rectangular Wannier flux seed marker: {marker}")

for forbidden in [
    "rectangular/disentangled transform is not supported yet",
    "DC-LCFO Wannier flux seed: non-square transform",
    "skip flux eigen seed for rectangular transform",
]:
    if forbidden in src:
        raise SystemExit(f"rectangular Wannier flux seed path still has obsolete guard: {forbidden}")

print("Rectangular Wannier90 transforms build Flux eigen seed by diagonalizing U^dag E U")
