#!/usr/bin/env python3
import argparse
from pathlib import Path
import shlex
import subprocess

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument("--build-dir", default="build-mpi-eigenexa")
args = parser.parse_args()
build = Path(args.build_dir)
if not build.is_absolute():
    build = ROOT / build

obj = Path("/tmp/test_dg_wpw_face_trace_scanner.o")
exe = Path("/tmp/test_dg_wpw_face_trace_scanner")
subprocess.run(["/opt/homebrew/bin/mpifort", "-I", str(build), "-c",
                str(ROOT / "tests/dg/test_dg_wpw_face_trace_scanner.f90"),
                "-o", str(obj)], check=True)
link_dir = build / "src"
link_args = shlex.split((link_dir / "CMakeFiles/salmon.dir/link.txt").read_text())
link_args.insert(1, str(obj))
link_args.remove("CMakeFiles/salmon.dir/main.f90.o")
link_args[link_args.index("-o") + 1] = str(exe)
subprocess.run(link_args, cwd=link_dir, check=True)
subprocess.run([str(exe)], check=True)
