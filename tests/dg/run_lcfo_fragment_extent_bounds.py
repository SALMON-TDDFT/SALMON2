#!/usr/bin/env python3
"""Compile and run the optimized-fragment extent bounds fixture."""

from pathlib import Path
import shutil
import subprocess
import tempfile


root = Path(__file__).resolve().parents[2]
compiler = shutil.which("gfortran")
if compiler is None:
    raise SystemExit("gfortran is required")
with tempfile.TemporaryDirectory(prefix="lcfo-fragment-extent-") as name:
    executable = Path(name) / "extent_bounds"
    subprocess.run(
        [compiler, "-O0", "-fcheck=all", "-fbacktrace", str(root / "tests/dg/test_lcfo_fragment_extent_bounds.f90"), "-o", str(executable)],
        check=True,
    )
    subprocess.run([str(executable)], check=True)
