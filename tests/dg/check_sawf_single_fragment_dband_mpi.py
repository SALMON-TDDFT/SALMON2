#!/usr/bin/env python3
from pathlib import Path

source = (Path(__file__).resolve().parents[2] / "src/gs/dc/lcfo_flux.f90").read_text()

required = [
    "dc%id_frag/=0 .and. dc%n_frag/=1",
    "modulo(p-1,dc%isize_frag)/=dc%id_frag",
    "expected_local_points",
]
missing = [token for token in required if token not in source]
if missing:
    raise SystemExit(
        "SAWF single-fragment D_band MPI partition missing: " + ", ".join(missing)
    )

print("PASS SAWF single-fragment D_band grid points are MPI partitioned")
