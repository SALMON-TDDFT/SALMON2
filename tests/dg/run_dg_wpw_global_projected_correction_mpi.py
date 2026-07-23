#!/usr/bin/env python3
from pathlib import Path
import os
import subprocess

ROOT = Path(__file__).resolve().parents[2]
MOD = Path("/tmp/test_dg_wpw_global_projected_correction_mod")
EXE = Path("/tmp/test_dg_wpw_global_projected_correction_mpi")
MOD.mkdir(exist_ok=True)

subprocess.run([
    os.environ.get("MPIFC", "mpifort"),
    "-J", str(MOD), "-I", str(MOD),
    str(ROOT / "src/common/dg_wpw_global_projected_correction.f90"),
    str(ROOT / "tests/dg/test_dg_wpw_global_projected_correction_mpi.f90"),
    "-llapack", "-lblas", "-o", str(EXE),
], check=True, cwd="/tmp")
subprocess.run(["mpirun", "-np", "2", str(EXE)], check=True, cwd=ROOT)
