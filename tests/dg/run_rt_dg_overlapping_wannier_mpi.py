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
        observable_prefix = build / f"observables-{ranks}"
        env["OW_RT_OBSERVABLE_PREFIX"] = str(observable_prefix)
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
        one_shot = Path(f"{observable_prefix}-one.dat")
        split = Path(f"{observable_prefix}-split.dat")
        assert one_shot.read_bytes() == split.read_bytes(), (
            "one-shot and restart-split observable evidence differ",
            one_shot.read_text(),
            split.read_text(),
        )
        lines = one_shot.read_text().splitlines()
        assert sum("step time Ex Ey Ez Px Py Pz Jx Jy Jz" in line for line in lines) == 1
        data = [line.split() for line in lines if line and not line.startswith("#")]
        assert len(data) == 3, data
        assert all(len(row) == 11 for row in data), data
        steps = [int(row[0]) for row in data]
        times = [float(row[1]) for row in data]
        assert steps == [0, 1, 2], steps
        assert all(right > left for left, right in zip(times, times[1:])), times
        assert not list(build.glob(f"observables-{ranks}-*.rank*.dat")), (
            "observable publication must be rank-zero-only"
        )
    assert len(set(fingerprints)) == 1, fingerprints
print("PASS overlapping-Wannier coefficient RT fixture on 1, 2, 4, and 8 ranks")
