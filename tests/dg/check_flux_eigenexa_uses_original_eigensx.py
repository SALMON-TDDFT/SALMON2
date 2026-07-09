#!/usr/bin/env python3
from pathlib import Path


root = Path(__file__).resolve().parents[2]
src = (root / "src" / "gs" / "dc" / "lcfo_flux.f90").read_text()
ref = (root / "src" / "gs" / "dc" / "lcfo.f90").read_text()

if "call eigen_sx(n, n, h_div, nx" not in ref:
    raise SystemExit("reference LCFO EigenExa path no longer uses eigen_sx")

if "call eigen_sx(n, n, h_div, nx, esp_tot(1:n,ispin), v_div, nx)" not in src:
    raise SystemExit("Flux-LCFO EigenExa path must use eigen_sx like the original LCFO path")

if "call eigen_s(n, n, h_div, nx, esp_tot(1:n,ispin), v_div, nx)" in src:
    raise SystemExit("Flux-LCFO EigenExa path still calls eigen_s instead of eigen_sx")

print("Flux-LCFO EigenExa diagonalization uses original LCFO eigen_sx call")
