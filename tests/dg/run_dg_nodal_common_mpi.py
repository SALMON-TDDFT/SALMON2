#!/usr/bin/env python3
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[2]
sources = [
    root / "src/common/dg_nodal_state.f90",
    root / "src/common/dg_nodal_interfaces.f90",
    root / "tests/dg/test_dg_nodal_common_mpi.f90",
]
for source in sources:
    assert source.exists(), f"missing nodal common source: {source}"

mpifort = shutil.which("mpifort")
mpiexec = shutil.which("mpiexec")
assert mpifort and mpiexec, "MPI compiler and launcher are required"

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    (tmp / "config.h").write_text("")
    exe = tmp / "test_dg_nodal_common_mpi"
    subprocess.run(
        [
            mpifort,
            "-cpp",
            "-DUSE_MPI",
            "-I",
            str(tmp),
            "-fcheck=all",
            "-fbacktrace",
            *map(str, sources),
            "-o",
            str(exe),
        ],
        check=True,
        cwd=tmp,
    )
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    completed = subprocess.run(
        [mpiexec, "-n", "2", str(exe)],
        text=True,
        capture_output=True,
        env=env,
    )
    if completed.returncode != 0:
        raise AssertionError(
            f"MPI fixture failed ({completed.returncode})\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    assert "PASS GS/RT-neutral nodal common state" in completed.stdout

print("PASS nodal common state MPI fixture")
