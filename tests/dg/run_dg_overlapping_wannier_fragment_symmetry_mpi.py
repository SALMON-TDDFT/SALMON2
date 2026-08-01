#!/usr/bin/env python3
"""Compile and run exact buffered-fragment subgroup selection on MPI ranks."""

from pathlib import Path
import os
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-fragment-symmetry-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    executable = build / "fragment-symmetry"
    compile_result = subprocess.run([
        shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-I", str(build), "-J", str(build),
        "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
        str(ROOT / "src/gs/dc/dg_overlapping_wannier_symmetry.f90"),
        str(ROOT / "tests/dg/test_dg_overlapping_wannier_fragment_symmetry_mpi.f90"),
        "-o", str(executable),
    ], capture_output=True, text=True)
    if compile_result.returncode:
        raise RuntimeError("fragment-symmetry compile failed:\n" + compile_result.stdout + compile_result.stderr)
    environment = os.environ.copy()
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for ranks in (1, 2, 4, 8):
        result = subprocess.run([shutil.which("mpiexec"), "-n", str(ranks), str(executable)],
                                capture_output=True, text=True, env=environment)
        if result.returncode:
            raise RuntimeError(f"{ranks}-rank fragment symmetry failed:\n{result.stdout}{result.stderr}")
        expected = f"PASS exact buffered-fragment symmetry on {ranks} ranks"
        if expected not in result.stdout:
            raise RuntimeError(f"missing {ranks}-rank PASS marker:\n{result.stdout}")
print("PASS exact buffered-fragment symmetry fixture on 1, 2, 4, and 8 ranks")
