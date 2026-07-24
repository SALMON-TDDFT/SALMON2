#!/usr/bin/env python3
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[2]
module = root / "src/gs/dc/dg_dc_local_basis_ground_state.f90"
fixture = root / "tests/dg/test_dg_dc_local_basis_layout_mpi.f90"
if not module.exists():
    raise AssertionError("missing DG-DC local-basis ground-state module")

fc = shutil.which("mpifort")
launcher = shutil.which("mpiexec")
assert fc and launcher
with tempfile.TemporaryDirectory() as tmp_name:
    tmp = Path(tmp_name)
    (tmp / "config.h").write_text("")
    executable = tmp / "test_dg_dc_local_basis_layout_mpi"
    subprocess.run(
        [
            fc,
            "-cpp",
            "-DUSE_MPI",
            "-I",
            str(tmp),
            str(module),
            str(fixture),
            "-o",
            str(executable),
        ],
        check=True,
        cwd=tmp,
    )
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    result = subprocess.run(
        [launcher, "-n", "2", str(executable)],
        text=True,
        capture_output=True,
        env=env,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"MPI fixture failed ({result.returncode})\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    assert "PASS DG-DC local basis layout keeps global bands independent" in result.stdout

print("PASS DG-DC local-basis layout MPI fixture")
