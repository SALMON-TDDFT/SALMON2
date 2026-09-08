#!/usr/bin/env python3
"""Compile/run complex spinor optimized-extent bounds coverage."""

from pathlib import Path
import shutil
import subprocess
import tempfile


root = Path(__file__).resolve().parents[2]
compiler = shutil.which("gfortran")
if compiler is None:
    raise SystemExit("gfortran is required")
with tempfile.TemporaryDirectory(prefix="lcfo-soi-fragment-extent-") as holder:
    executable = Path(holder) / "extent_bounds"
    command = [compiler, "-O0", "-fcheck=all", "-fbacktrace",
               str(root / "tests/dg/test_lcfo_soi_fragment_extent_bounds.f90"),
               "-o", str(executable)]
    print("compile:", " ".join(command), flush=True)
    subprocess.run(command, check=True)
    subprocess.run([str(executable)], check=True)
