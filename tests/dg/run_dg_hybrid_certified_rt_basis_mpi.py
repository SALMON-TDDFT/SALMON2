#!/usr/bin/env python3
"""Build and exercise the certified localized Hybrid RT basis contract."""

from pathlib import Path
import os
import re
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "src/gs/dc/dg_hybrid_certified_rt_basis.f90"
if MODULE.exists():
    source = MODULE.read_text().lower()
    builder = source[source.index("subroutine build_dg_hybrid_certified_rt_basis"):]
    builder = builder[:builder.index("end subroutine build_dg_hybrid_certified_rt_basis")]
    for forbidden in ("where(abs(", "pack(abs(", "drop_small", "sparsify", "pruning_threshold"):
        assert forbidden not in builder, f"certified RT builder independently prunes operator elements: {forbidden}"
    assert "call localizer(" in builder, "certified RT builder does not invoke the second localizer"

with tempfile.TemporaryDirectory(prefix="hybrid-certified-rt-basis-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    executable = build / "test_dg_hybrid_certified_rt_basis_mpi"
    subprocess.run([
        shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-std=f2008", "-ffree-line-length-none",
        "-I", str(build), "-J", str(build), "-fcheck=all",
        "-ffpe-trap=invalid,zero,overflow", "-fbacktrace",
        str(ROOT / "src/gs/dc/dg_hybrid_ground_state_types.f90"),
        str(ROOT / "src/gs/dc/dg_hybrid_low_energy_symmetry.f90"),
        str(MODULE),
        str(ROOT / "tests/dg/test_dg_hybrid_certified_rt_basis_mpi.f90"),
        "-o", str(executable),
    ], check=True)
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    fingerprints = []
    for ranks in (1, 2, 4, 8):
        completed = subprocess.run(
            [shutil.which("mpiexec"), "-n", str(ranks), str(executable)],
            capture_output=True, text=True, env=environment, timeout=120,
        )
        assert completed.returncode == 0, (ranks, completed.stdout, completed.stderr)
        assert f"PASS certified Hybrid RT basis on {ranks} ranks" in completed.stdout
        match = re.search(
            r"HYBRID_CERTIFIED_RT_BASIS ranks=\d+ c=(-?\d+) u=(-?\d+) b=(-?\d+) "
            r"operators=(-?\d+) fingerprint=(-?\d+)", completed.stdout,
        )
        assert match, completed.stdout
        fingerprints.append(tuple(int(value) for value in match.groups()))
    assert len(set(fingerprints)) == 1, fingerprints

print("PASS certified Hybrid RT basis on 1, 2, 4, and 8 ranks")
