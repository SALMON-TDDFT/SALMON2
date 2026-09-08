#!/usr/bin/env python3
from pathlib import Path
import os
import shutil
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="hybrid-real-space-residual-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    exe = build / "test_real_space_residual"
    subprocess.run([
        shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-std=f2008",
        "-ffree-line-length-none", "-I", str(build), "-J", str(build),
        "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-fbacktrace",
        str(root / "src/gs/dc/dg_hybrid_real_space_residual.f90"),
        str(root / "tests/dg/test_dg_hybrid_real_space_residual_mpi.f90"),
        "-o", str(exe),
    ], check=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for ranks in (1, 2, 4):
        run = subprocess.run([shutil.which("mpiexec"), "-n", str(ranks), str(exe)],
                             capture_output=True, text=True, env=env)
        assert run.returncode == 0, (ranks, run.stdout, run.stderr)
        assert f"PASS hybrid real-space residual on {ranks} ranks" in run.stdout
print("PASS hybrid real-space residual on 1, 2, and 4 ranks")
