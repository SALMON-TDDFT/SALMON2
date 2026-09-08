#!/usr/bin/env python3
from pathlib import Path
import os
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (ROOT / "src/gs/dc/dg_hybrid_schwarz_operator.f90").read_text().lower()
if "subroutine apply_dg_hybrid_schwarz_rows" in SOURCE:
    APPLY = SOURCE.split("subroutine apply_dg_hybrid_schwarz_rows", 1)[1].split(
        "end subroutine apply_dg_hybrid_schwarz_rows", 1
    )[0]
    assert "mpi_allgather" not in APPLY and "mpi_allgatherv" not in APPLY
    assert "mpi_irecv" in APPLY and "mpi_isend" in APPLY and "mpi_waitall" in APPLY

with tempfile.TemporaryDirectory(prefix="hybrid-schwarz-operator-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    executable = build / "hybrid_schwarz_operator"
    subprocess.run(
        [
            shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-std=f2008",
            "-ffree-line-length-none", "-I", str(build), "-J", str(build),
            "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-fbacktrace",
            str(ROOT / "src/gs/dc/dg_hybrid_schwarz_operator.f90"),
            str(ROOT / "tests/dg/test_dg_hybrid_schwarz_operator_mpi.f90"),
            "-o", str(executable),
        ],
        check=True,
    )
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for rank_count in (2, 4, 8):
        run = subprocess.run(
            [shutil.which("mpiexec"), "-n", str(rank_count), str(executable)],
            capture_output=True, text=True, env=environment, timeout=30,
        )
        assert run.returncode == 0, (rank_count, run.stdout, run.stderr)
        assert f"PASS hybrid Schwarz schedule on {rank_count} ranks" in run.stdout

print("PASS hybrid Schwarz schedule on 2, 4, and 8 ranks")
