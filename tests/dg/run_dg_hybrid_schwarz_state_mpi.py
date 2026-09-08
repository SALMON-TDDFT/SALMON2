#!/usr/bin/env python3
from pathlib import Path
import os
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]

with tempfile.TemporaryDirectory(prefix="hybrid-schwarz-state-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    executable = build / "hybrid_schwarz_state"
    subprocess.run(
        [
            shutil.which("mpifort"),
            "-cpp",
            "-DUSE_MPI",
            "-std=f2008",
            "-ffree-line-length-none",
            "-I",
            str(build),
            "-J",
            str(build),
            "-fcheck=all",
            "-ffpe-trap=invalid,zero,overflow",
            "-fbacktrace",
            str(ROOT / "src/gs/dc/dg_hybrid_schwarz_state.f90"),
            str(ROOT / "tests/dg/test_dg_hybrid_schwarz_state_mpi.f90"),
            "-o",
            str(executable),
        ],
        check=True,
    )
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for rank_count in (2, 4, 8):
        run = subprocess.run(
            [shutil.which("mpiexec"), "-n", str(rank_count), str(executable)],
            capture_output=True,
            text=True,
            env=environment,
            timeout=30,
        )
        assert run.returncode == 0, (rank_count, run.stdout, run.stderr)
        assert f"PASS hybrid Schwarz state on {rank_count} ranks" in run.stdout

print("PASS hybrid Schwarz state on 2, 4, and 8 ranks")
