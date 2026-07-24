#!/usr/bin/env python3
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[2]
adapter = root / "src/gs/dc/dg_dc_ground_state_adapter.f90"
assert adapter.exists(), f"missing DG-DC production adapter: {adapter}"

sources = [
    root / "src/common/dg_nodal_state.f90",
    root / "src/common/dg_nodal_sipg.f90",
    root / "src/gs/dc/dg_dc_ground_state.f90",
    adapter,
    root / "tests/dg/test_dg_dc_ground_state_adapter_mpi.f90",
]
mpifort = shutil.which("mpifort")
mpiexec = shutil.which("mpiexec")
assert mpifort and mpiexec

with tempfile.TemporaryDirectory() as tmp_name:
    tmp = Path(tmp_name)
    (tmp / "config.h").write_text("")
    exe = tmp / "test_dg_dc_ground_state_adapter_mpi"
    subprocess.run(
        [mpifort, "-cpp", "-DUSE_MPI", "-I", str(tmp), "-fcheck=all",
         "-fbacktrace", *map(str, sources), "-o", str(exe)],
        cwd=tmp, check=True,
    )
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    result = subprocess.run(
        [mpiexec, "-n", "2", str(exe)], text=True, capture_output=True, env=env
    )
    if result.returncode:
        raise AssertionError(
            f"MPI fixture failed ({result.returncode})\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    assert "PASS DG-DC production adapter" in result.stdout

print("PASS DG-DC production adapter MPI fixture")
