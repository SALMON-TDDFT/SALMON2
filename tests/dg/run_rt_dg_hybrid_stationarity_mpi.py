#!/usr/bin/env python3
from pathlib import Path
import os
import shutil
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="hybrid-stationarity-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    exe = build / "hybrid_stationarity"
    subprocess.run([
        shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-I", str(build),
        "-J", str(build), "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
        "-fbacktrace", str(root / "src/common/dg_hybrid_continuation_residuals.f90"),
        str(root / "src/common/dg_hybrid_total_energy.f90"),
        str(root / "src/rt/dg/rt_dg_hybrid_stationarity.f90"),
        str(root / "tests/dg/test_rt_dg_hybrid_stationarity_mpi.f90"),
        "-o", str(exe),
    ], check=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for ranks in (1, 2, 4):
        run = subprocess.run([shutil.which("mpiexec"), "-n", str(ranks), str(exe)],
                             capture_output=True, text=True, env=env, timeout=90)
        assert run.returncode == 0, (ranks, run.stdout, run.stderr)
        assert f"PASS hybrid RT stationarity on {ranks} ranks" in run.stdout
print("PASS hybrid RT stationarity on 1, 2, and 4 ranks")
