#!/usr/bin/env python3
"""Build and exercise the localization-first seed and receipt contract."""

from __future__ import annotations

import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "src/gs/dc/dg_hybrid_localization_first.f90"


def lapack_libraries() -> list[str]:
    configured = os.environ.get("SALMON_LAPACK_LIBS")
    if configured:
        return shlex.split(configured)
    if shutil.which("pkg-config"):
        openblas = subprocess.run(
            ["pkg-config", "--libs", "openblas"], capture_output=True, text=True
        )
        if openblas.returncode == 0:
            return shlex.split(openblas.stdout)
        lapack = subprocess.run(
            ["pkg-config", "--libs", "lapack"], capture_output=True, text=True
        )
        if lapack.returncode == 0:
            return shlex.split(lapack.stdout)
    if shutil.which("brew"):
        prefix = subprocess.run(
            ["brew", "--prefix", "openblas"], capture_output=True, text=True
        )
        if prefix.returncode == 0:
            return [f"-L{prefix.stdout.strip()}/lib", "-lopenblas"]
    return ["-llapack", "-lblas"]


def require_narrow_metric_route() -> None:
    source = MODULE.read_text().lower()
    start = source.index("subroutine prepare_dg_hybrid_localization_first_seed")
    end = source.index("end subroutine prepare_dg_hybrid_localization_first_seed", start)
    body = source[start:end]
    assert "call orthonormalize_dg_distributed_seed_space" in body
    forbidden = (
        "group_averaged",
        "character_sector",
        "symmetry_closed",
        "orbit_completion",
        "per_wf",
        "prun",
    )
    for token in forbidden:
        assert token not in body, f"localization-first seed preparation calls {token}"


def main() -> None:
    require_narrow_metric_route()
    compiler = shutil.which("mpifort")
    launcher = shutil.which("mpiexec")
    if compiler is None or launcher is None:
        raise SystemExit("mpifort and mpiexec are required")
    with tempfile.TemporaryDirectory(prefix="dg-localization-first-") as name:
        build = Path(name)
        (build / "config.h").write_text("")
        executable = build / "test_dg_hybrid_localization_first"
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
                str(ROOT / "src/gs/dc/dg_overlapping_wannier_construction.f90"),
                str(MODULE),
                str(ROOT / "tests/dg/test_dg_hybrid_localization_first_mpi.f90"),
                *lapack_libraries(),
                "-o",
                str(executable),
            ],
            check=True,
            timeout=120,
        )
        environment = os.environ.copy()
        environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
        receipts: list[tuple[int, int, float]] = []
        for ranks in (1, 2, 4, 8):
            completed = subprocess.run(
                [launcher, "-n", str(ranks), str(executable)],
                cwd=build,
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=60,
            )
            assert completed.returncode == 0, (ranks, completed.stdout)
            assert f"PASS localization-first contract on {ranks} ranks" in completed.stdout
            match = re.search(
                r"LOCALIZATION_FIRST ranks=\d+ seed=(-?\d+) "
                r"transform=(-?\d+) spread_total=\s*([+\-0-9.Ee]+)",
                completed.stdout,
            )
            assert match, completed.stdout
            receipts.append((int(match.group(1)), int(match.group(2)), float(match.group(3))))
        assert len(set(value[0] for value in receipts)) == 1, receipts
        assert len(set(value[1] for value in receipts)) == 1, receipts
        assert len(set(value[2] for value in receipts)) == 1, receipts
    print("PASS localization-first contract on 1, 2, 4, and 8 ranks")


if __name__ == "__main__":
    main()
