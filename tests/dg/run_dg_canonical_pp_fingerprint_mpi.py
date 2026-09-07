#!/usr/bin/env python3
"""Build and exercise the canonical pseudopotential fingerprint."""

from pathlib import Path
import os
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    compiler = shutil.which("mpifort")
    launcher = shutil.which("mpiexec")
    if compiler is None or launcher is None:
        raise SystemExit("mpifort and mpiexec are required")
    with tempfile.TemporaryDirectory(prefix="dg-canonical-pp-fingerprint-") as name:
        build = Path(name)
        (build / "config.h").write_text("")
        executable = build / "test_dg_canonical_pp_fingerprint"
        subprocess.run(
            [
                compiler,
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
                str(ROOT / "src/common/structures.f90"),
                str(ROOT / "src/gs/dc/dg_canonical_pp_fingerprint.f90"),
                str(ROOT / "tests/dg/test_dg_canonical_pp_fingerprint_mpi.f90"),
                "-o",
                str(executable),
            ],
            check=True,
        )
        environment = os.environ.copy()
        environment["OMP_NUM_THREADS"] = "1"
        environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
        for ranks in (1, 2, 4):
            completed = subprocess.run(
                [launcher, "-n", str(ranks), str(executable)],
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=30,
            )
            assert completed.returncode == 0, (ranks, completed.stdout)
            assert f"PASS canonical PP fingerprint on {ranks} ranks" in completed.stdout
    print("PASS canonical PP fingerprint on 1, 2, and 4 ranks")


if __name__ == "__main__":
    main()
