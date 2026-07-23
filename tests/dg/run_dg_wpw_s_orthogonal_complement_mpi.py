#!/usr/bin/env python3
from pathlib import Path
import os
import subprocess

ROOT = Path(__file__).resolve().parents[2]
BUILD = Path(os.environ.get("SALMON_BUILD_DIR", ROOT / "build-mpi-eigenexa"))
FC = os.environ.get("MPIFC", "mpifort")
EXE = Path("/tmp/test_dg_wpw_s_orthogonal_complement_mpi")
MOD = Path("/tmp/test_dg_wpw_s_orthogonal_complement_mod")
MOD.mkdir(exist_ok=True)

sources = [
    ROOT / "src/common/dg_wpw_owner_exchange.f90",
    ROOT / "src/common/dg_wpw_bounded_operator.f90",
    ROOT / "src/common/dg_wpw_s_orthogonal_complement.f90",
    ROOT / "tests/dg/test_dg_wpw_s_orthogonal_complement_mpi.f90",
]
cmd = [FC, "-J", str(MOD), "-I", str(MOD), "-I", str(BUILD),
       *(str(path) for path in sources), "-llapack", "-lblas", "-o", str(EXE)]
subprocess.run(cmd, cwd="/tmp", check=True)
subprocess.run(["mpirun", "-np", "2", str(EXE)], cwd=ROOT, check=True)
