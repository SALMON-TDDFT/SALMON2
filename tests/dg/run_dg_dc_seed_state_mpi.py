#!/usr/bin/env python3
"""Build and exercise the conventional-DC seed state helpers."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    compiler = shutil.which("mpifort")
    launcher = shutil.which("mpiexec")
    if compiler is None or launcher is None:
        raise SystemExit("mpifort and mpiexec are required")

    with tempfile.TemporaryDirectory(prefix="dg-dc-seed-state-") as name:
        build = Path(name)
        (build / "config.h").write_text("")
        executable = build / "test_dg_dc_seed_state"
        subprocess.run(
            [
                compiler,
                "-cpp",
                "-DUSE_MPI",
                "-I",
                str(build),
                "-J",
                str(build),
                "-fcheck=all",
                "-ffpe-trap=invalid,zero,overflow",
                "-fbacktrace",
                str(ROOT / "src/gs/dc/dg_dc_seed_checkpoint.f90"),
                str(ROOT / "tests/dg/test_dg_dc_seed_state_mpi.f90"),
                "-o",
                str(executable),
            ],
            check=True,
        )

        environment = os.environ.copy()
        environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
        for ranks in (1, 2, 4, 8):
            case = build / f"ranks-{ranks}"
            seed = case / "seed"
            seed.mkdir(parents=True)
            environment["DG_DC_SEED_DIRECTORY"] = str(seed)
            completed = subprocess.run(
                [launcher, "-n", str(ranks), str(executable)],
                cwd=case,
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=30,
            )
            assert completed.returncode == 0, (ranks, completed.stdout)
            assert f"PASS DG DC seed state ranks={ranks}" in completed.stdout, (
                completed.stdout
            )

    print("PASS conventional DC seed state on 1, 2, 4, and 8 ranks")


if __name__ == "__main__":
    main()
