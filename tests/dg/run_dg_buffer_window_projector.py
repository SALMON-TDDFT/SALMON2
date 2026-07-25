#!/usr/bin/env python3
from pathlib import Path
import os
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
fc = os.environ.get("MPIFC", "mpifort")

with tempfile.TemporaryDirectory(prefix="dg-buffer-projector-") as temporary_directory:
    build = Path(temporary_directory)
    exe = build / "test_dg_buffer_window_projector"
    subprocess.run(
        [
            fc,
            "-J",
            str(build),
            "-I",
            str(build),
            str(ROOT / "src/common/dg_buffer_window_projector.f90"),
            str(ROOT / "tests/dg/test_dg_buffer_window_projector.f90"),
            "-fcheck=all",
            "-llapack",
            "-lblas",
            "-o",
            str(exe),
        ],
        check=True,
        cwd=build,
    )
    subprocess.run([str(exe)], check=True, cwd=ROOT, timeout=10)
