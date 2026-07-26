#!/usr/bin/env python3
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[2]
sources = [
    root / "src/common/dg_nodal_state.f90",
    root / "src/common/dg_nodal_sipg.f90",
    root / "src/gs/dc/dg_dc_ground_state.f90",
    root / "src/gs/dc/dg_dc_local_basis_ground_state.f90",
    root / "src/common/dg_ground_state_checkpoint.f90",
    root / "tests/dg/test_dg_ground_state_checkpoint_mpi.f90",
]
for source in sources:
    assert source.exists(), f"missing DG ground-state checkpoint source: {source}"

mpifort = shutil.which("mpifort")
mpiexec = shutil.which("mpiexec")
assert mpifort and mpiexec

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    (tmp / "config.h").write_text("")
    exe = tmp / "test_dg_ground_state_checkpoint_mpi"
    subprocess.run(
        [mpifort, "-cpp", "-DUSE_MPI", "-I", str(tmp), "-fcheck=all",
         "-fbacktrace", *map(str, sources), "-llapack", "-lblas", "-o", str(exe)],
        cwd=tmp, check=True,
    )
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    result = subprocess.run(
        [mpiexec, "-n", "2", str(exe), str(tmp / "accepted_dg_gs")],
        text=True, capture_output=True, env=env,
    )
    if result.returncode:
        raise AssertionError(
            f"MPI fixture failed ({result.returncode})\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    assert "PASS verified DG ground-state checkpoint" in result.stdout

print("PASS DG ground-state checkpoint MPI fixture")
