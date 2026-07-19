#!/usr/bin/env python3
from pathlib import Path
import subprocess, tempfile

root = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory() as td:
    exe = Path(td) / "test"
    subprocess.run(["gfortran", "-std=f2008", "-Wall", "-Wextra", "-fcheck=all", "-J", td,
                    str(root/"src/rt/dg/rt_dg_wpw_length_gauge.f90"),
                    str(root/"tests/dg/test_wpw_length_gauge_observable.f90"),
                    "-o", str(exe)], check=True, cwd=td)
    subprocess.run([str(exe)], check=True)
print("WPW length-gauge observable: PASS")
