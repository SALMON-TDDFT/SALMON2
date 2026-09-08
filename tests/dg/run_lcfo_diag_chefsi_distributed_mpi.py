#!/usr/bin/env python3
"""Build and run an actual two-rank distributed CheFSI fixture."""

from pathlib import Path
import shutil
import subprocess
import tempfile


root = Path(__file__).resolve().parents[2]
compiler = shutil.which("mpifort")
launcher = shutil.which("mpirun")
scalapack = Path("/opt/homebrew/opt/scalapack/lib/libscalapack.dylib")
openblas = Path("/opt/homebrew/opt/openblas/lib/libopenblas.dylib")
if compiler is None or launcher is None or not scalapack.exists() or not openblas.exists():
    raise SystemExit("MPI, ScaLAPACK, and OpenBLAS are required")

with tempfile.TemporaryDirectory(prefix="lcfo-chefsi-distributed-") as name:
    build = Path(name)
    flags = ["-cpp", "-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffree-line-length-none", "-J", str(build), "-I", str(build)]
    sources = (
        root / "tests/dg/lcfo_diag_chefsi_runtime_stubs.f90",
        root / "src/gs/dc/lcfo_diag_chefsi.f90",
        root / "tests/dg/test_lcfo_diag_chefsi_distributed_mpi.f90",
    )
    objects = []
    print(subprocess.check_output([compiler, "--version"], text=True).splitlines()[0])
    print(subprocess.check_output([launcher, "--version"], text=True).splitlines()[0])
    print(f"scalapack={scalapack.resolve()}")
    print(f"openblas={openblas.resolve()}")
    print("bounds_check=-fcheck=all")
    for index, source in enumerate(sources):
        obj = build / f"source-{index}.o"
        command = [compiler, *flags, "-c", str(source), "-o", str(obj)]
        print("compile:", " ".join(command))
        subprocess.run(command, check=True)
        objects.append(obj)
    executable = build / "test_chefsi_distributed"
    link = [compiler, "-fcheck=all", "-fbacktrace", *(str(item) for item in objects), str(scalapack), str(openblas), "-o", str(executable)]
    print("link:", " ".join(link))
    subprocess.run(link, check=True)
    result = subprocess.run([launcher, "-np", "2", str(executable)], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(result.stdout, end="")
    if "illegal value" in result.stdout.lower():
        raise SystemExit("ScaLAPACK emitted an illegal-parameter diagnostic")
    if "CheFSI cycle:" not in result.stdout:
        raise SystemExit("distributed fixture did not exercise CheFSI filtering")
    if result.returncode:
        raise SystemExit(result.returncode)
