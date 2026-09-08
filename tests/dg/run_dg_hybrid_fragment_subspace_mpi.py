#!/usr/bin/env python3
"""Small distributed bounded-update tests; never starts a material calculation."""
from pathlib import Path
import os
import shlex
import shutil
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
if shutil.which("pkg-config"):
    lapack = shlex.split(subprocess.check_output(["pkg-config", "--libs", "openblas"], text=True))
else:
    prefix = subprocess.check_output(["brew", "--prefix", "openblas"], text=True).strip()
    lapack = [f"-L{prefix}/lib", "-lopenblas"]
with tempfile.TemporaryDirectory(prefix="hybrid-fragment-subspace-") as name:
    build = Path(name)
    executable = build / "fragment_subspace"
    subprocess.run([
        shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-std=f2008",
        "-J", str(build), "-I", str(build), "-fcheck=all",
        "-ffpe-trap=invalid,zero,overflow", "-fbacktrace",
        str(root / "src/gs/occupation_kernel.f90"),
        str(root / "src/gs/dc/dc_fragment_occupation.f90"),
        str(root / "src/gs/dc/dg_hybrid_fragment_subspace.f90"),
        str(root / "tests/dg/test_dg_hybrid_fragment_subspace_mpi.f90"),
        *lapack, "-o", str(executable),
    ], check=True)
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment["OPENBLAS_NUM_THREADS"] = "1"
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for ranks in (1, 2, 4, 8):
        run = subprocess.run([shutil.which("mpiexec"), "-n", str(ranks), str(executable)],
                             env=environment, capture_output=True, text=True, timeout=60)
        assert run.returncode == 0, (ranks, run.stdout, run.stderr)
        assert f"PASS hybrid fragment subspace on {ranks} ranks" in run.stdout, run.stdout
        print(run.stdout.strip())
