#!/usr/bin/env python3
from pathlib import Path
import os
import re
import shutil
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
source = (root / "src/gs/dc/dg_hybrid_low_energy_symmetry.f90").read_text().lower()
adaptive = source[source.index("subroutine certify_dg_hybrid_energy_window"):]
adaptive = adaptive[:adaptive.index("end subroutine certify_dg_hybrid_energy_window")]
assert not re.search(r"\b(?:128|384)\b", adaptive), "material-specific rank literal entered energy-window policy"
assert "energy_window==-1d0" in adaptive, "exact -1 dynamic-rank compatibility branch is missing"
assert "compatibility" in adaptive and "warning" in adaptive, "dynamic-rank compatibility warning is missing"
assert "precompute_dg_hybrid_symmetry_prefix_defects" in adaptive, "prefix defects are not precomputed once"
assert "call evaluate_dg_hybrid_low_energy_symmetry" not in adaptive, (
    "adaptive cluster search repeats the full dense symmetry evaluator")
with tempfile.TemporaryDirectory(prefix="hybrid-low-energy-symmetry-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    exe = build / "test_dg_hybrid_low_energy_symmetry_mpi"
    subprocess.run([
        shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-std=f2008", "-ffree-line-length-none",
        "-I", str(build), "-J", str(build), "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-fbacktrace",
        str(root / "src/gs/dc/dg_hybrid_ground_state_types.f90"),
        str(root / "src/gs/dc/dg_hybrid_low_energy_symmetry.f90"),
        str(root / "tests/dg/test_dg_hybrid_low_energy_symmetry_mpi.f90"), "-o", str(exe),
    ], check=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    fingerprints = []
    for nrank in (1, 2, 4, 8):
        run = subprocess.run([shutil.which("mpiexec"), "-n", str(nrank), str(exe)],
                             capture_output=True, text=True, env=env)
        assert run.returncode == 0, (nrank, run.stdout, run.stderr)
        assert f"PASS low-energy Hybrid symmetry on {nrank} ranks" in run.stdout
        assert "compatibility" in (run.stdout + run.stderr).lower(), (run.stdout, run.stderr)
        match = re.search(r"HYBRID_SPECTRAL_CERTIFICATION ranks=\d+ fingerprint=(-?\d+)", run.stdout)
        assert match, run.stdout
        fingerprints.append(int(match.group(1)))
    assert len(set(fingerprints)) == 1, fingerprints
print("PASS low-energy Hybrid symmetry on 1, 2, 4, and 8 ranks")
