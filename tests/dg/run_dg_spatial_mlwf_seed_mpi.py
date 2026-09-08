#!/usr/bin/env python3
"""Compile and run the spatial-MLWF sparse-seed contract."""

from pathlib import Path
import argparse
import os
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument("--ranks", default="1,2,4")
args = parser.parse_args()
ranks = [int(value) for value in args.ranks.split(",")]

with tempfile.TemporaryDirectory(prefix="spatial-mlwf-seed-") as name:
    build = Path(name)
    (build / "config.h").write_text("#define USE_MPI\n")
    executable = build / "seed"
    compiler = shutil.which("mpifort")
    if compiler is None:
        raise RuntimeError("mpifort is required")
    compile_result = subprocess.run([
        compiler, "-cpp", "-I", str(build), "-J", str(build),
        "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
        str(ROOT / "src/gs/dc/dg_spatial_mlwf_seed.f90"),
        str(ROOT / "tests/dg/test_dg_spatial_mlwf_seed_mpi.f90"),
        "-llapack", "-lblas", "-o", str(executable),
    ], capture_output=True, text=True)
    if compile_result.returncode:
        raise RuntimeError("seed compile failed:\n" + compile_result.stdout + compile_result.stderr)
    environment = os.environ.copy()
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for rank_count in ranks:
        result = subprocess.run(
            [shutil.which("mpiexec"), "-n", str(rank_count), str(executable)],
            capture_output=True, text=True, env=environment,
        )
        if result.returncode:
            raise RuntimeError(f"{rank_count}-rank seed failed:\n{result.stdout}{result.stderr}")
        marker = f"PASS spatial MLWF sparse seed on {rank_count} ranks"
        if marker not in result.stdout:
            raise RuntimeError(f"missing {rank_count}-rank PASS marker:\n{result.stdout}")
print("PASS spatial MLWF sparse seed fixture on " + ", ".join(map(str, ranks)) + " ranks")
