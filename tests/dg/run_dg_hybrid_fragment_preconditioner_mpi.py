#!/usr/bin/env python3
"""Distributed fixed-frame action and bounded-update covariance fixtures."""
from pathlib import Path
import os
import shlex
import shutil
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
if shutil.which("pkg-config"):
    lapack = shlex.split(subprocess.check_output(["pkg-config", "--libs", "openblas"], text=True))
else:
    prefix = subprocess.check_output(["brew", "--prefix", "openblas"], text=True).strip()
    lapack = [f"-L{prefix}/lib", "-lopenblas"]
with tempfile.TemporaryDirectory(prefix="hybrid-fragment-preconditioner-") as name:
    build = Path(name)
    exe = build / "preconditioner"
    subprocess.run([
        shutil.which("mpifort"), "-std=f2008", "-J", str(build), "-I", str(build),
        "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-fbacktrace",
        str(root / "src/gs/occupation_kernel.f90"),
        str(root / "src/gs/dc/dc_fragment_occupation.f90"),
        str(root / "src/gs/dc/dg_hybrid_fragment_preconditioner.f90"),
        str(root / "src/gs/dc/dg_hybrid_fragment_subspace.f90"),
        str(root / "tests/dg/test_dg_hybrid_fragment_preconditioner_mpi.f90"),
        *lapack, "-o", str(exe),
    ], check=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = env["OPENBLAS_NUM_THREADS"] = "1"
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for ranks in (1, 2, 4, 8):
        run = subprocess.run([shutil.which("mpiexec"), "-n", str(ranks), str(exe)],
                             env=env, capture_output=True, text=True, timeout=60)
        assert run.returncode == 0, (ranks, run.stdout, run.stderr)
        assert f"PASS fragment preconditioner on {ranks} ranks" in run.stdout, run.stdout
        print(run.stdout.strip())
