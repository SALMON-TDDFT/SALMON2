#!/usr/bin/env python3
"""Compile the production face-topology helper and exercise it on two ranks."""

from pathlib import Path
import os
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
FC = shutil.which("mpifort") or shutil.which("mpif90")
MPIEXEC = shutil.which("mpirun") or shutil.which("mpiexec")
if not FC or not MPIEXEC:
    raise SystemExit("MPI Fortran compiler/launcher is required")

with tempfile.TemporaryDirectory(prefix="dc-fragment-face-") as holder:
    build = Path(holder)
    sources = [
        ROOT / "tests/dg/dc_fragment_geometry_test_structures.f90",
        ROOT / "src/gs/dc/dc_fragment_geometry.f90",
        ROOT / "tests/dg/test_dc_fragment_face_neighbor_mpi.f90",
    ]
    objects = []
    for index, source in enumerate(sources):
        obj = build / f"source-{index}.o"
        command = [FC, "-O0", "-g", "-fcheck=all", "-fbacktrace", "-J", str(build),
                   "-I", str(build), "-c", str(source), "-o", str(obj)]
        print("compile:", " ".join(command), flush=True)
        subprocess.run(command, check=True)
        objects.append(obj)
    executable = build / "test_dc_fragment_face_neighbor"
    command = [FC, "-fcheck=all", "-fbacktrace", *map(str, objects), "-o", str(executable)]
    print("link:", " ".join(command), flush=True)
    subprocess.run(command, check=True)
    environment = os.environ.copy()
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    command = [MPIEXEC, "-np", "2", str(executable)]
    print("run:", " ".join(command), flush=True)
    subprocess.run(command, check=True, env=environment)
