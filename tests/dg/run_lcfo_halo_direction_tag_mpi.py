#!/usr/bin/env python3
"""Compile/run the LCFO same-dvec basis-halo tag MPI regression."""

from argparse import ArgumentParser
from pathlib import Path
import os
import shutil
import subprocess
import tempfile


parser = ArgumentParser()
parser.add_argument("--legacy", action="store_true")
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
fc = shutil.which("mpifort") or shutil.which("mpif90")
launcher = shutil.which("mpirun") or shutil.which("mpiexec")
if not fc or not launcher:
    raise SystemExit("MPI Fortran compiler/launcher is required")

with tempfile.TemporaryDirectory(prefix="lcfo-halo-tag-") as holder:
    build = Path(holder)
    sources = [
        root / "tests/dg/dc_fragment_geometry_test_structures.f90",
        root / "src/gs/dc/dc_fragment_geometry.f90",
        root / "tests/dg/test_lcfo_halo_direction_tag_mpi.F90",
    ]
    objects = []
    for index, source in enumerate(sources):
        obj = build / f"source-{index}.o"
        command = [fc, "-cpp", "-O0", "-g", "-fcheck=all", "-fbacktrace",
                   "-J", str(build), "-I", str(build)]
        if args.legacy:
            command.append("-DLEGACY_TAG")
        command.extend(["-c", str(source), "-o", str(obj)])
        print("compile:", " ".join(command), flush=True)
        subprocess.run(command, check=True)
        objects.append(obj)
    executable = build / "test_lcfo_halo_direction_tag"
    command = [fc, "-fcheck=all", "-fbacktrace", *map(str, objects), "-o", str(executable)]
    print("link:", " ".join(command), flush=True)
    subprocess.run(command, check=True)
    environment = os.environ.copy()
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for ranks in (2, 3):
        command = [launcher, "-np", str(ranks), str(executable)]
        print("run:", " ".join(command), flush=True)
        subprocess.run(command, check=True, env=environment, timeout=30)
