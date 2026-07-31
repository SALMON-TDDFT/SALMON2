#!/usr/bin/env python3
from pathlib import Path
import os
import re
import shutil
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-coefficient-rt-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    exe = build / "coefficient_rt"
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
            str(root / "src/rt/dg/rt_dg_overlapping_wannier.f90"),
            str(root / "tests/dg/test_rt_dg_overlapping_wannier_mpi.f90"),
            "-L/opt/homebrew/opt/openblas/lib",
            "-lopenblas",
            "-o",
            str(exe),
        ],
        check=True,
    )
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    fingerprints = []
    for ranks in (1, 2, 4, 8):
        env["OW_RT_RESTART_PREFIX"] = str(build / f"restart-{ranks}.bin")
        result = subprocess.run(
            [shutil.which("mpiexec"), "-n", str(ranks), str(exe)],
            capture_output=True,
            text=True,
            env=env,
        )
        assert result.returncode == 0, (ranks, result.stdout, result.stderr)
        match = re.search(
            rf"PASS overlapping-Wannier coefficient RT on {ranks} ranks fingerprint=([0-9A-F]+)",
            result.stdout,
        )
        assert match, result.stdout
        fingerprints.append(match.group(1))
    assert len(set(fingerprints)) == 1, fingerprints
print("PASS overlapping-Wannier coefficient RT fixture on 1, 2, 4, and 8 ranks")
