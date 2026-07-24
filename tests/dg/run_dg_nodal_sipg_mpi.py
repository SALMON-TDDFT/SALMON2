#!/usr/bin/env python3
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[2]
sources = [
    root / "src/common/dg_nodal_sipg.f90",
    root / "tests/dg/test_dg_nodal_sipg_mpi.f90",
]
for source in sources:
    assert source.exists(), f"missing SIPG source: {source}"
fc = shutil.which("mpifort")
run = shutil.which("mpiexec")
assert fc and run
with tempfile.TemporaryDirectory() as directory:
    directory = Path(directory)
    (directory / "config.h").write_text("")
    exe = directory / "test_dg_nodal_sipg_mpi"
    subprocess.run(
        [fc, "-cpp", "-DUSE_MPI", "-fcheck=all", "-fbacktrace", "-I", str(directory),
         *map(str, sources), "-o", str(exe)],
        check=True,
        cwd=directory,
    )
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    completed = subprocess.run([run, "-n", "2", str(exe)], text=True, capture_output=True, env=env)
    if completed.returncode:
        raise AssertionError(f"SIPG fixture failed\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}")
    assert "PASS analytic complete Hermitian nodal SIPG face operator" in completed.stdout
print("PASS nodal SIPG MPI fixture")
