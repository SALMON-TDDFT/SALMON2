#!/usr/bin/env python3
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]

with tempfile.TemporaryDirectory(prefix="dg_wpw_checkpoint_") as tmp:
    tmp = Path(tmp)
    exe = tmp / "test_dg_wpw_checkpoint_roundtrip"
    subprocess.run([
        "gfortran", "-std=f2008", "-J", str(tmp),
        str(ROOT / "src/common/dg_wpw_checkpoint.f90"),
        str(ROOT / "tests/dg/test_dg_wpw_checkpoint_roundtrip.f90"),
        "-o", str(exe),
    ], check=True, cwd=tmp)
    subprocess.run([str(exe)], check=True, cwd=tmp)
