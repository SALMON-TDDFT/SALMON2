#!/usr/bin/env python3
from pathlib import Path
import os
import shutil
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
mpifort = shutil.which("mpifort")
mpiexec = shutil.which("mpiexec")
assert mpifort and mpiexec, "MPI compiler and launcher are required"

with tempfile.TemporaryDirectory(prefix="ow-projection-") as temporary:
    build = Path(temporary)
    (build / "config.h").write_text("")
    executable = build / "projection"
    subprocess.run(
        [
            mpifort,
            "-cpp",
            "-DUSE_MPI",
            "-I",
            str(build),
            "-J",
            str(build),
            "-fcheck=all",
            "-ffpe-trap=invalid,zero,overflow",
            "-fbacktrace",
            str(root / "src/gs/dc/dg_overlapping_wannier_projection.f90"),
            str(root / "tests/dg/test_dg_overlapping_wannier_projection_mpi.f90"),
            "-o",
            str(executable),
        ],
        check=True,
    )
    environment = {**os.environ, "OMP_NUM_THREADS": "1"}
    for ranks in (1, 2, 4, 8):
        completed = subprocess.run(
            [mpiexec, "-n", str(ranks), str(executable)],
            text=True,
            capture_output=True,
            env=environment,
            timeout=30,
        )
        if completed.returncode:
            raise SystemExit(completed.stdout + completed.stderr)

print("PASS complete-s+p projection manifest fixture on 1, 2, 4, and 8 ranks")
