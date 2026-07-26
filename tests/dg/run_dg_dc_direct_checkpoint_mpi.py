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
    root / "tests/dg/test_dg_dc_direct_checkpoint_mpi.f90",
]
mpifort = shutil.which("mpifort")
mpiexec = shutil.which("mpiexec")
assert mpifort and mpiexec
with tempfile.TemporaryDirectory() as tmp_name:
    tmp = Path(tmp_name)
    (tmp / "config.h").write_text("")
    exe = tmp / "test_dg_dc_direct_checkpoint_mpi"
    subprocess.run(
        [mpifort, "-cpp", "-DUSE_MPI", "-I", str(tmp), "-fcheck=all", "-fbacktrace",
         *map(str, sources), "-framework", "Accelerate", "-o", str(exe)], check=True, cwd=tmp,
    )
    case = tmp / "checkpoint"
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    completed = subprocess.run([mpiexec, "-n", "2", str(exe), str(case)],
                               text=True, capture_output=True, env=env)
    if completed.returncode:
        raise AssertionError(
            f"MPI fixture failed ({completed.returncode})\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    assert "PASS direct DC checkpoint MPI fixture" in completed.stdout
print("PASS direct DC checkpoint MPI fixture")
