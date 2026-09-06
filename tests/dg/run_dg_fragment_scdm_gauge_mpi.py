#!/usr/bin/env python3
"""Compile and run the deterministic fragment-SCDM gauge contract."""

from pathlib import Path
import os
import re
import shlex
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
if os.environ.get("SALMON_LAPACK_LIBS"):
    lapack = shlex.split(os.environ["SALMON_LAPACK_LIBS"])
elif shutil.which("pkg-config") and subprocess.run(
    ["pkg-config", "--exists", "openblas"], check=False
).returncode == 0:
    lapack = shlex.split(subprocess.check_output(["pkg-config", "--libs", "openblas"], text=True))
elif shutil.which("brew"):
    prefix = subprocess.check_output(["brew", "--prefix", "openblas"], text=True).strip()
    lapack = [f"-L{prefix}/lib", "-lopenblas"]
else:
    lapack = ["-llapack", "-lblas"]

with tempfile.TemporaryDirectory(prefix="fragment-scdm-gauge-") as name:
    build = Path(name)
    (build / "config.h").write_text("#define USE_MPI\n")
    executable = build / "fragment_scdm_gauge"
    compile_result = subprocess.run(
        [
            shutil.which("mpifort"),
            "-cpp",
            "-I",
            str(build),
            "-J",
            str(build),
            "-O0",
            "-g",
            "-fcheck=all",
            "-ffpe-trap=invalid,zero,overflow",
            "-fbacktrace",
            str(ROOT / "src/gs/dc/dg_fragment_scdm_gauge.f90"),
            str(ROOT / "tests/dg/test_dg_fragment_scdm_gauge_mpi.f90"),
            *lapack,
            "-o",
            str(executable),
        ],
        capture_output=True,
        text=True,
    )
    if compile_result.returncode:
        raise RuntimeError("fragment SCDM compile failed:\n" + compile_result.stdout + compile_result.stderr)
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    fingerprints = []
    for rank_count in (1, 2, 4, 8):
        result = subprocess.run(
            [shutil.which("mpiexec"), "-n", str(rank_count), str(executable)],
            capture_output=True,
            text=True,
            env=environment,
        )
        if result.returncode:
            raise RuntimeError(
                f"{rank_count}-rank fragment SCDM test failed:\n{result.stdout}{result.stderr}"
            )
        marker = f"PASS fragment SCDM gauge on {rank_count} ranks"
        if marker not in result.stdout:
            raise RuntimeError(f"missing {rank_count}-rank PASS marker:\n{result.stdout}")
        match = re.search(r"SCDM_GAUGE ranks=\d+ fingerprint=(-?\d+)", result.stdout)
        if not match:
            raise RuntimeError(f"missing {rank_count}-rank fingerprint:\n{result.stdout}")
        fingerprints.append(int(match.group(1)))
    if len(set(fingerprints)) != 1:
        raise RuntimeError(f"layout-dependent SCDM fingerprints: {fingerprints}")
print("PASS fragment SCDM gauge on 1, 2, 4, and 8 ranks")
