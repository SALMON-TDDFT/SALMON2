#!/usr/bin/env python3
from pathlib import Path
import os
import shlex
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]


def lapack_libraries():
    configured = os.environ.get("SALMON_LAPACK_LIBS")
    if configured:
        return shlex.split(configured)
    if shutil.which("pkg-config"):
        for package in ("openblas", "lapack"):
            probe = subprocess.run(["pkg-config", "--libs", package], capture_output=True, text=True)
            if probe.returncode == 0:
                return shlex.split(probe.stdout)
    if shutil.which("brew"):
        probe = subprocess.run(["brew", "--prefix", "openblas"], capture_output=True, text=True)
        if probe.returncode == 0:
            return [f"-L{probe.stdout.strip()}/lib", "-lopenblas"]
    return ["-llapack", "-lblas"]


with tempfile.TemporaryDirectory(prefix="hybrid-schwarz-solver-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    executable = build / "hybrid_schwarz_solver"
    subprocess.run(
        [
            shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-std=f2008",
            "-ffree-line-length-none", "-I", str(build), "-J", str(build),
            "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-fbacktrace",
            str(ROOT / "src/gs/dc/dg_hybrid_schwarz_state.f90"),
            str(ROOT / "src/gs/dc/dg_hybrid_schwarz_solver.f90"),
            str(ROOT / "tests/dg/test_dg_hybrid_schwarz_solver_mpi.f90"),
            *lapack_libraries(), "-o", str(executable),
        ], check=True,
    )
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for rank_count in (2, 4, 8):
        run = subprocess.run(
            [shutil.which("mpiexec"), "-n", str(rank_count), str(executable)],
            capture_output=True, text=True, env=environment, timeout=30,
        )
        assert run.returncode == 0, (rank_count, run.stdout, run.stderr)
        assert f"PASS hybrid Schwarz solver on {rank_count} ranks" in run.stdout

print("PASS hybrid Schwarz solver on 2, 4, and 8 ranks")
