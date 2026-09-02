#!/usr/bin/env python3
from pathlib import Path
import os
import re
import shutil
import subprocess
import tempfile


root = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="dc-fragment-occupation-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    executable = build / "dc_fragment_occupation"
    subprocess.run(
        [
            shutil.which("mpifort"),
            "-cpp",
            "-DUSE_MPI",
            "-std=f2008",
            "-I",
            str(build),
            "-J",
            str(build),
            "-fcheck=all",
            "-ffpe-trap=invalid,zero,overflow",
            "-fbacktrace",
            str(root / "src/gs/occupation_kernel.f90"),
            str(root / "src/gs/dc/dc_fragment_occupation.f90"),
            str(root / "tests/dg/test_dc_fragment_occupation_mpi.f90"),
            "-o",
            str(executable),
        ],
        check=True,
    )
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    fingerprints = []
    for rank_count in (1, 2, 4):
        run = subprocess.run(
            [shutil.which("mpiexec"), "-n", str(rank_count), str(executable)],
            capture_output=True,
            text=True,
            env=environment,
            timeout=30,
        )
        assert run.returncode == 0, (rank_count, run.stdout, run.stderr)
        assert f"PASS DC fragment occupation on {rank_count} ranks" in run.stdout
        match = re.search(
            r"DC_FRAGMENT_OCCUPATION ranks=\d+ fingerprint=(-?\d+)", run.stdout
        )
        assert match, run.stdout
        fingerprints.append(int(match.group(1)))
    assert len(set(fingerprints)) == 1, fingerprints
print("PASS DC fragment occupation on 1, 2, and 4 ranks")
