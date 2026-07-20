#!/usr/bin/env python3
from pathlib import Path
import os
import subprocess

ROOT = Path(__file__).resolve().parents[2]
FC = os.environ.get("MPIFC", "mpifort")
EXE = Path("/tmp/test_dg_wpw_core_wannier_seed_mpi")
MOD = Path("/tmp/test_dg_wpw_core_wannier_seed_mod")
MOD.mkdir(exist_ok=True)

sources = [
    ROOT / "src/gs/dc/lcfo_wannier_sawf_seed.f90",
    ROOT / "src/gs/dc/dg_wpw_wannier_tail_halo.f90",
    ROOT / "tests/dg/test_dg_wpw_core_wannier_seed_mpi.f90",
]
subprocess.run(
    [FC, "-J", str(MOD), "-I", str(MOD), *(str(path) for path in sources),
     "-llapack", "-lblas", "-o", str(EXE)],
    cwd="/tmp",
    check=True,
)
subprocess.run(["mpirun", "-np", "2", str(EXE)], cwd=ROOT, check=True)
