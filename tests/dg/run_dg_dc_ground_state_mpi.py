#!/usr/bin/env python3
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[2]
sources = [
    root / "src/common/dg_nodal_state.f90",
    root / "src/gs/dc/dg_dc_ground_state.f90",
    root / "tests/dg/test_dg_dc_ground_state_mpi.f90",
]
for source in sources:
    assert source.exists(), f"missing DG-DC ground-state source: {source}"

mpifort = shutil.which("mpifort")
mpiexec = shutil.which("mpiexec")
assert mpifort and mpiexec

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    (tmp / "config.h").write_text("")
    exe = tmp / "test_dg_dc_ground_state_mpi"
    subprocess.run(
        [mpifort, "-cpp", "-DUSE_MPI", "-I", str(tmp), "-fcheck=all", "-fbacktrace",
         *map(str, sources), "-o", str(exe)],
        cwd=tmp, check=True,
    )
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    result = subprocess.run(
        [mpiexec, "-n", "2", str(exe)], text=True, capture_output=True, env=env
    )
    if result.returncode:
        raise AssertionError(
            f"MPI fixture failed ({result.returncode})\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    assert "PASS adaptive self-consistent DG continuation" in result.stdout

print("PASS DG-DC ground-state continuation MPI fixture")
