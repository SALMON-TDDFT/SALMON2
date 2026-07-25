#!/usr/bin/env python3
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCES = [
    ROOT / "src/common/dg_buffer_window_projector.f90",
    ROOT / "src/gs/dc/dg_dc_buffer_core_faces.f90",
    ROOT / "tests/dg/test_dg_dc_buffer_core_faces_mpi.f90",
]
mpifort = shutil.which("mpifort")
mpiexec = shutil.which("mpiexec")
assert mpifort and mpiexec, "MPI compiler and launcher are required"

with tempfile.TemporaryDirectory(prefix="dg-buffer-core-faces-") as temporary_directory:
    build = Path(temporary_directory)
    (build / "config.h").write_text("")
    executable = build / "test_dg_dc_buffer_core_faces_mpi"
    subprocess.run(
        [
            mpifort,
            "-cpp",
            "-DUSE_MPI",
            "-I",
            str(build),
            "-J",
            str(build),
            "-fcheck=all",
            "-fbacktrace",
            "-g",
            *map(str, SOURCES),
            "-llapack",
            "-lblas",
            "-o",
            str(executable),
        ],
        check=True,
        cwd=build,
    )
    environment = os.environ.copy()
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for process_count in (1, 2, 8):
        completed = subprocess.run(
            [mpiexec, "-n", str(process_count), str(executable)],
            text=True,
            capture_output=True,
            env=environment,
            timeout=60,
        )
        if completed.returncode != 0:
            raise AssertionError(
                f"{process_count}-rank fixture failed ({completed.returncode})\n"
                f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
            )
        assert "PASS buffer-to-neighbor-core physical face mapping" in completed.stdout

print("PASS DG DC buffer/core face MPI fixture on 1, 2, and 8 ranks")
