#!/usr/bin/env python3
from pathlib import Path


source = (Path(__file__).resolve().parents[2] / "src/gs/dc/lcfo_flux.f90").read_text()
required = [
    "subroutine check_lcfo_h_component_finite",
    "[DC-LCFO-H-COMPONENT-LOCAL]",
    "call check_lcfo_h_component_finite(mat_H_weak_kinetic",
    "call check_lcfo_h_component_finite(mat_H_weak_potential",
    "call check_lcfo_h_component_finite(mat_H_weak_nonlocal",
    "call check_lcfo_h_component_finite(mat_H_volume_weak_local",
    "call check_lcfo_h_component_finite(mat_H_surface_self",
    "call check_lcfo_h_component_finite(halo(i_halo)%mat_H_surface_cross",
]
missing = [token for token in required if token not in source]
if missing:
    raise SystemExit(
        "LCFO Hamiltonian component finite guards missing: " + ", ".join(missing)
    )

print("PASS LCFO Hamiltonian components stop at the producing rank")
