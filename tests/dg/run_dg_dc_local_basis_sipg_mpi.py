#!/usr/bin/env python3
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[2]
sources = [
    root / "src/common/dg_nodal_sipg.f90",
    root / "src/gs/dc/dg_dc_local_basis_ground_state.f90",
    root / "tests/dg/test_dg_dc_local_basis_sipg_mpi.f90",
]
fc = shutil.which("mpifort")
launcher = shutil.which("mpiexec")
assert fc and launcher and all(path.exists() for path in sources)
with tempfile.TemporaryDirectory() as tmp_name:
    tmp = Path(tmp_name)
    (tmp / "config.h").write_text("")
    executable = tmp / "test_dg_dc_local_basis_sipg_mpi"
    subprocess.run(
        [fc, "-cpp", "-DUSE_MPI", "-I", str(tmp), *map(str, sources), "-o", str(executable)],
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
    assert "PASS DG-DC local-basis SIPG pair assembly" in result.stdout

print("PASS DG-DC local-basis SIPG MPI fixture")
