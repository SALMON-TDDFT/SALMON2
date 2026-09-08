#!/usr/bin/env python3
"""Run retained-basis crystallographic group representations on MPI ranks."""

from pathlib import Path
import os
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-point-group-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    executable = build / "point-group"
    result = subprocess.run([
        shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-I", str(build), "-J", str(build),
        "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
        str(ROOT / "src/gs/dc/dg_overlapping_wannier_symmetry.f90"),
        str(ROOT / "tests/dg/test_dg_overlapping_wannier_point_group_mpi.f90"),
        "-llapack", "-lblas", "-o", str(executable),
    ], capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError("point-group compile failed:\n" + result.stdout + result.stderr)
    environment = os.environ.copy()
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for ranks in (1, 2, 4, 8):
        result = subprocess.run([shutil.which("mpiexec"), "-n", str(ranks), str(executable)],
                                capture_output=True, text=True, env=environment)
        if result.returncode:
            raise RuntimeError(f"{ranks}-rank point group failed:\n{result.stdout}{result.stderr}")
        if f"PASS overlapping-Wannier point group on {ranks} ranks" not in result.stdout:
            raise RuntimeError(f"missing {ranks}-rank PASS marker:\n{result.stdout}")
print("PASS overlapping-Wannier point-group fixture on 1, 2, 4, and 8 ranks")
