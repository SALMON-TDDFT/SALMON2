#!/usr/bin/env python3
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="dg_dc_rt_handoff_") as tmp:
    tmp = Path(tmp)
    exe = tmp / "test_dg_dc_rt_handoff_mpi"
    subprocess.run([
        "mpifort", "-std=f2008", "-J", str(tmp),
        str(ROOT / "src/common/dg_wpw_owner_exchange.f90"),
        str(ROOT / "src/common/dg_wpw_bounded_operator.f90"),
        str(ROOT / "src/common/dg_wpw_checkpoint.f90"),
        str(ROOT / "src/rt/dg/rt_dg_wpw_checkpoint_handoff.f90"),
        str(ROOT / "tests/dg/test_dg_dc_rt_handoff_mpi.f90"),
        "-o", str(exe),
    ], check=True, cwd=tmp)
    subprocess.run(["mpirun", "-np", "2", str(exe)], check=True, cwd=tmp)
