#!/usr/bin/env python3
from pathlib import Path
import os
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-w90-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    exe = build / "w90_adapter"
    subprocess.run(
        [
            shutil.which("mpifort"),
            "-cpp",
            "-I",
            str(build),
            "-J",
            str(build),
            "-fcheck=all",
            "-ffpe-trap=invalid,zero,overflow",
            str(ROOT / "src/gs/dc/dg_overlapping_wannier_w90.f90"),
            str(ROOT / "tests/dg/test_dg_overlapping_wannier_w90_mpi.f90"),
            "-o",
            str(exe),
        ],
        check=True,
    )
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for ranks in (1, 2, 4, 8):
        result = subprocess.run(
            [shutil.which("mpiexec"), "-n", str(ranks), str(exe)],
            capture_output=True,
            text=True,
            env=env,
        )
        assert result.returncode == 0, (ranks, result.stdout, result.stderr)
        assert "PASS Wannier90 MLWF adapter validation" in result.stdout
    library = os.environ.get("SALMON_WANNIER90_LIB")
    if library:
        actual = build / "w90_library"
        subprocess.run(
            [
                shutil.which("mpifort"),
                "-cpp",
                "-DUSE_MPI",
                "-DUSE_WANNIER90",
                "-I",
                str(build),
                "-J",
                str(build),
                str(ROOT / "src/gs/dc/dg_overlapping_wannier_w90.f90"),
                str(ROOT / "tests/dg/test_dg_overlapping_wannier_w90_mpi.f90"),
                library,
                "-framework",
                "Accelerate",
                "-o",
                str(actual),
            ],
            check=True,
        )
        for ranks in (1, 2, 4, 8):
            result = subprocess.run(
                [shutil.which("mpiexec"), "-n", str(ranks), str(actual)],
                cwd=build,
                capture_output=True,
                text=True,
                env=env,
            )
            assert result.returncode == 0, (ranks, result.stdout, result.stderr)
            assert "PASS Wannier90 MLWF adapter validation" in result.stdout
print("PASS Wannier90 MLWF adapter validation on 1, 2, 4, and 8 ranks")
