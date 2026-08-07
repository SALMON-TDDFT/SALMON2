#!/usr/bin/env python3
"""Compile and run buffer-local periodic Wannier localization checks."""

from pathlib import Path
import os
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-localization-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    executable = build / "localization"
    compile_result = subprocess.run([
        shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-I", str(build), "-J", str(build),
        "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
        str(ROOT / "src/gs/dc/dg_overlapping_wannier_localization.f90"),
        str(ROOT / "tests/dg/test_dg_overlapping_wannier_localization_mpi.f90"),
        "-o", str(executable),
    ], capture_output=True, text=True)
    if compile_result.returncode:
        raise RuntimeError("localization compile failed:\n" + compile_result.stdout + compile_result.stderr)
    environment = os.environ.copy()
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for ranks in (1, 2, 4, 8):
        result = subprocess.run([shutil.which("mpiexec"), "-n", str(ranks), str(executable)],
                                capture_output=True, text=True, env=environment)
        if result.returncode:
            raise RuntimeError(f"{ranks}-rank localization failed:\n{result.stdout}{result.stderr}")
        expected = f"PASS buffer-local periodic Wannier spread on {ranks} ranks"
        if expected not in result.stdout:
            raise RuntimeError(f"missing {ranks}-rank PASS marker:\n{result.stdout}")
print("PASS buffer-local periodic Wannier spread fixture on 1, 2, 4, and 8 ranks")
