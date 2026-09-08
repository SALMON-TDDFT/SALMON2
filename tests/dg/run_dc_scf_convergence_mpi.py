#!/usr/bin/env python3
from pathlib import Path
import os
import shutil
import subprocess
import tempfile


root = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="dc-scf-convergence-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    executable = build / "dc_scf_convergence"
    subprocess.run(
        [
            shutil.which("mpifort"),
            "-cpp",
            "-DUSE_MPI",
            "-I",
            str(build),
            "-J",
            str(build),
            "-fcheck=all",
            "-ffpe-trap=invalid,zero,overflow",
            "-fbacktrace",
            str(root / "src/gs/dc/dc_scf_convergence.f90"),
            str(root / "tests/dg/test_dc_scf_convergence_mpi.f90"),
            "-o",
            str(executable),
        ],
        check=True,
    )
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for rank_count in (1, 2, 4):
        run = subprocess.run(
            [shutil.which("mpiexec"), "-n", str(rank_count), str(executable)],
            capture_output=True,
            text=True,
            env=environment,
        )
        assert run.returncode == 0, (rank_count, run.stdout, run.stderr)
        assert f"PASS DC SCF convergence on {rank_count} ranks" in run.stdout
print("PASS DC SCF convergence on 1, 2, and 4 ranks")
