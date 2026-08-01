#!/usr/bin/env python3
"""Run scalar and polar-vector fragment covariance projection on MPI ranks."""

from pathlib import Path
import os
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-symmetry-projection-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    executable = build / "symmetry-projection"
    result = subprocess.run([
        shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-I", str(build), "-J", str(build),
        "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
        str(ROOT / "src/gs/dc/dg_overlapping_wannier_symmetry.f90"),
        str(ROOT / "tests/dg/test_dg_overlapping_wannier_symmetry_projection_mpi.f90"),
        "-llapack", "-lblas", "-o", str(executable),
    ], capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError("symmetry-projection compile failed:\n" + result.stdout + result.stderr)
    environment = os.environ.copy()
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for ranks in (1, 2, 4, 8):
        result = subprocess.run([shutil.which("mpiexec"), "-n", str(ranks), str(executable)],
                                capture_output=True, text=True, env=environment)
        if result.returncode:
            raise RuntimeError(f"{ranks}-rank symmetry projection failed:\n{result.stdout}{result.stderr}")
        if f"PASS overlapping-Wannier symmetry projection on {ranks} ranks" not in result.stdout:
            raise RuntimeError(f"missing {ranks}-rank PASS marker:\n{result.stdout}")
print("PASS overlapping-Wannier symmetry projection fixture on 1, 2, 4, and 8 ranks")
