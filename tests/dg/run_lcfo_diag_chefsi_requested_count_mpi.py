#!/usr/bin/env python3
"""Build and run actual CheFSI against a one-rank ScaLAPACK fixture."""

from argparse import ArgumentParser
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
PRE_FIX_BLOB = "130a9c12c2065f495e8c57e6826e614e82c76bb8"
SCALAPACK = Path("/opt/homebrew/opt/scalapack/lib/libscalapack.dylib")
OPENBLAS = Path("/opt/homebrew/opt/openblas/lib/libopenblas.dylib")

parser = ArgumentParser()
parser.add_argument("--source", choices=("current", "pre-fix"), default="current")
args = parser.parse_args()

compiler = shutil.which("mpifort")
launcher = shutil.which("mpirun")
if compiler is None or launcher is None:
    raise SystemExit("mpifort and mpirun are required")
for library in (SCALAPACK, OPENBLAS):
    if not library.exists():
        raise SystemExit(f"required library is missing: {library}")

with tempfile.TemporaryDirectory(prefix="lcfo-chefsi-runtime-") as name:
    build = Path(name)
    solver = ROOT / "src/gs/dc/lcfo_diag_chefsi.f90"
    old_interface = False
    if args.source == "pre-fix":
        solver = build / "lcfo_diag_chefsi.pre-fix.f90"
        solver.write_bytes(
            subprocess.check_output(["git", "cat-file", "blob", PRE_FIX_BLOB], cwd=ROOT)
        )
        old_interface = True

    flags = [
        "-cpp",
        "-O0",
        "-g",
        "-fcheck=all",
        "-fbacktrace",
        "-ffree-line-length-none",
        "-J",
        str(build),
        "-I",
        str(build),
    ]
    stubs = ROOT / "tests/dg/lcfo_diag_chefsi_runtime_stubs.f90"
    driver = ROOT / "tests/dg/test_lcfo_diag_chefsi_requested_count_mpi.f90"
    objects = [build / "stubs.o", build / "solver.o", build / "driver.o"]
    commands = [
        [compiler, *flags, "-c", str(stubs), "-o", str(objects[0])],
        [compiler, *flags, "-c", str(solver), "-o", str(objects[1])],
        [
            compiler,
            *flags,
            *(["-DOLD_CHEFSI_INTERFACE"] if old_interface else []),
            "-c",
            str(driver),
            "-o",
            str(objects[2]),
        ],
        [
            compiler,
            "-fcheck=all",
            "-fbacktrace",
            *(str(item) for item in objects),
            str(SCALAPACK),
            str(OPENBLAS),
            "-o",
            str(build / "test_chefsi"),
        ],
    ]

    print(f"source={args.source}")
    print(f"source_identity={PRE_FIX_BLOB if old_interface else solver}")
    print(subprocess.check_output([compiler, "--version"], text=True).splitlines()[0])
    print(subprocess.check_output([launcher, "--version"], text=True).splitlines()[0])
    print(f"scalapack={SCALAPACK.resolve()}")
    print(f"openblas={OPENBLAS.resolve()}")
    print("bounds_check=-fcheck=all")
    for command in commands:
        print("compile:", " ".join(command))
        subprocess.run(command, cwd=build, check=True)
    completed = subprocess.run(
        [launcher, "-np", "1", str(build / "test_chefsi")],
        cwd=build,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    print(completed.stdout, end="")
    illegal_diagnostic = re.compile(
        r"parameter\s+number|illegal\s+value|\bPDORMTR\b|\bXERBLA\b",
        re.IGNORECASE,
    )
    if illegal_diagnostic.search(completed.stdout):
        raise SystemExit("ScaLAPACK/XERBLA illegal-parameter diagnostic detected")
