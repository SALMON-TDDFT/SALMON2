#!/usr/bin/env python3
"""Build and run overlapping-Wannier metadata fixtures."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MPIFORT = shutil.which("mpifort")
MPIEXEC = shutil.which("mpiexec")
assert MPIFORT and MPIEXEC, "MPI compiler and launcher are required"

with tempfile.TemporaryDirectory(prefix="dg-overlapping-wannier-metadata-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    executable = build / "test_dg_overlapping_wannier_metadata_mpi"
    subprocess.run(
        [
            MPIFORT,
            "-cpp",
            "-DUSE_MPI",
            "-I",
            str(build),
            "-J",
            str(build),
            "-fcheck=all",
            "-ffpe-trap=invalid,zero,overflow",
            "-fbacktrace",
            "-g",
            str(ROOT / "src/gs/dc/dg_overlapping_wannier_types.f90"),
            str(ROOT / "tests/dg/test_dg_overlapping_wannier_metadata_mpi.f90"),
            "-o",
            str(executable),
        ],
        check=True,
        cwd=build,
    )
    environment = os.environ.copy()
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    fingerprints: list[int] = []
    for process_count in (1, 2, 4, 8):
        completed = subprocess.run(
            [MPIEXEC, "-n", str(process_count), str(executable)],
            text=True,
            capture_output=True,
            env=environment,
            timeout=90,
        )
        if completed.returncode != 0:
            raise AssertionError(
                f"{process_count}-rank fixture failed ({completed.returncode})\n"
                f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
            )
        marker = f"PASS overlapping-Wannier metadata on {process_count} ranks"
        assert marker in completed.stdout, completed.stdout
        match = re.search(r"fingerprint=(-?\d+)", completed.stdout)
        assert match, completed.stdout
        fingerprints.append(int(match.group(1)))
    assert len(set(fingerprints)) == 1, fingerprints

print("PASS overlapping-Wannier metadata fixture on 1, 2, 4, and 8 ranks")
