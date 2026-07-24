#!/usr/bin/env python3
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[2]
sources = [
    root / "src/common/dg_nodal_sipg.f90",
    root / "src/gs/dc/dg_dc_local_basis_ground_state.f90",
    root / "tests/dg/test_dg_dc_local_basis_sipg_mpi.f90",
]
axis_exchange_source = root / "tests/dg/test_dg_dc_local_basis_axis_exchange_mpi.f90"
fc = shutil.which("mpifort")
launcher = shutil.which("mpiexec")
assert fc and launcher and all(path.exists() for path in sources)
linear_algebra = ["-framework", "Accelerate"] if sys.platform == "darwin" else ["-llapack", "-lblas"]
with tempfile.TemporaryDirectory() as tmp_name:
    tmp = Path(tmp_name)
    (tmp / "config.h").write_text("")
    executable = tmp / "test_dg_dc_local_basis_sipg_mpi"
    subprocess.run(
        [fc, "-cpp", "-DUSE_MPI", "-I", str(tmp), *map(str, sources), *linear_algebra, "-o", str(executable)],
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
    axis_executable = tmp / "test_dg_dc_local_basis_axis_exchange_mpi"
    subprocess.run(
        [
            fc,
            "-cpp",
            "-DUSE_MPI",
            "-I",
            str(tmp),
            str(sources[0]),
            str(sources[1]),
            str(axis_exchange_source),
            *linear_algebra,
            "-o",
            str(axis_executable),
        ],
        check=True,
        cwd=tmp,
    )
    axis_result = subprocess.run(
        [launcher, "-n", "4", str(axis_executable)],
        text=True,
        capture_output=True,
        env=env,
    )
    if axis_result.returncode != 0:
        raise AssertionError(
            f"four-rank axis fixture failed ({axis_result.returncode})\n"
            f"stdout:\n{axis_result.stdout}\nstderr:\n{axis_result.stderr}"
        )
    assert "PASS DG-DC local-basis four-rank axis exchange" in axis_result.stdout

print("PASS DG-DC local-basis SIPG MPI fixture")
