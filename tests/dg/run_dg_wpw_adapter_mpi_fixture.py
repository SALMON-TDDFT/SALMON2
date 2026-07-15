#!/usr/bin/env python3
"""Build the WPW adapter fixture against the exact SALMON object set and run it."""
from pathlib import Path
import shlex
import subprocess

ROOT = Path(__file__).resolve().parents[2]
BUILD = ROOT / "build-mpi-eigenexa"
OBJ = Path("/tmp/test_dg_wpw_adapter_mpi.o")
EXE = Path("/tmp/test_dg_wpw_adapter_mpi")

subprocess.run([
    "/opt/homebrew/bin/mpifort", "-I", str(BUILD), "-c",
    str(ROOT / "tests/dg/test_dg_wpw_adapter_mpi.f90"), "-o", str(OBJ),
], check=True)

link_dir = BUILD / "src"
args = shlex.split((link_dir / "CMakeFiles/salmon.dir/link.txt").read_text())
args.insert(1, str(OBJ))
args.remove("CMakeFiles/salmon.dir/main.f90.o")
out_index = args.index("-o") + 1
args[out_index] = str(EXE)
subprocess.run(args, cwd=link_dir, check=True)
subprocess.run(["mpirun", "-np", "2", str(EXE)], check=True)
