#!/usr/bin/env python3
from pathlib import Path
import os
import shlex
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
if os.environ.get("SALMON_LAPACK_LIBS"):
    lapack_libs = shlex.split(os.environ["SALMON_LAPACK_LIBS"])
elif shutil.which("pkg-config") and subprocess.run(
    ["pkg-config", "--exists", "openblas"], check=False
).returncode == 0:
    lapack_libs = shlex.split(
        subprocess.check_output(["pkg-config", "--libs", "openblas"], text=True)
    )
else:
    generic_lapack = None
    if shutil.which("pkg-config"):
        probe = subprocess.run(["pkg-config", "--libs", "lapack"], capture_output=True, text=True)
        if probe.returncode == 0:
            generic_lapack = shlex.split(probe.stdout)
    if generic_lapack:
        lapack_libs = generic_lapack
    elif shutil.which("brew"):
        probe = subprocess.run(["brew", "--prefix", "openblas"], capture_output=True, text=True)
        if probe.returncode == 0:
            lapack_libs = [f"-L{probe.stdout.strip()}/lib", "-lopenblas"]
        else:
            lapack_libs = ["-llapack", "-lblas"]
    else:
        lapack_libs = ["-llapack", "-lblas"]
with tempfile.TemporaryDirectory(prefix="ow-w90-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    exe = build / "w90_adapter"
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
            str(ROOT / "src/gs/dc/lcfo_wannier_sawf_seed.f90"),
            str(ROOT / "src/gs/dc/dg_overlapping_wannier_w90.f90"),
            str(ROOT / "tests/dg/test_dg_overlapping_wannier_w90_mpi.f90"),
            *lapack_libs,
            "-o",
            str(exe),
        ],
        check=True,
    )
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    fingerprints = []
    sector_fingerprints = []
    position_tuple_fingerprints = []
    joint_center_fingerprints = []
    for ranks in (1, 2, 4, 8):
        result = subprocess.run(
            [shutil.which("mpiexec"), "-n", str(ranks), str(exe)],
            cwd=build,
            capture_output=True,
            text=True,
            env=env,
        )
        assert result.returncode == 0, (ranks, result.stdout, result.stderr)
        assert "PASS Wannier90 MLWF adapter validation" in result.stdout
        line = next(line for line in result.stdout.splitlines() if line.startswith("W90_MATRIX_FINGERPRINT"))
        fingerprints.append(float(line.split()[1]))
        sector_line = next(line for line in result.stdout.splitlines() if line.startswith("W90_SECTOR_FINGERPRINT"))
        sector_fingerprints.append(int(sector_line.split()[1]))
        position_line = next(
            line for line in result.stdout.splitlines() if line.startswith("W90_POSITION_TUPLE_FINGERPRINT")
        )
        position_tuple_fingerprints.append(int(position_line.split()[1]))
        joint_line = next(
            line for line in result.stdout.splitlines() if line.startswith("W90_JOINT_CENTER_FINGERPRINT")
        )
        joint_center_fingerprints.append(int(joint_line.split()[1]))
    assert max(fingerprints) - min(fingerprints) < 1.0e-12, fingerprints
    assert len(set(sector_fingerprints)) == 1, sector_fingerprints
    assert len(set(position_tuple_fingerprints)) == 1, position_tuple_fingerprints
    assert len(set(joint_center_fingerprints)) == 1, joint_center_fingerprints
    library = os.environ.get("SALMON_WANNIER90_LIB")
    if library:
        actual = build / "w90_library"
        subprocess.run(
            [
                shutil.which("mpifort"),
                "-cpp",
                "-DUSE_MPI",
                "-DUSE_WANNIER90",
                "-I",
                str(build),
                "-J",
                str(build),
                str(ROOT / "src/gs/dc/lcfo_wannier_sawf_seed.f90"),
                str(ROOT / "src/gs/dc/dg_overlapping_wannier_w90.f90"),
                str(ROOT / "tests/dg/test_dg_overlapping_wannier_w90_mpi.f90"),
                library,
                *lapack_libs,
                "-o",
                str(actual),
            ],
            check=True,
        )
        library_fingerprints = []
        for ranks in (1, 2, 4, 8):
            result = subprocess.run(
                [shutil.which("mpiexec"), "-n", str(ranks), str(actual)],
                cwd=build,
                capture_output=True,
                text=True,
                env=env,
            )
            assert result.returncode == 0, (ranks, result.stdout, result.stderr)
            assert "PASS Wannier90 MLWF adapter validation" in result.stdout
            line = next(line for line in result.stdout.splitlines() if line.startswith("W90_MATRIX_FINGERPRINT"))
            library_fingerprints.append(float(line.split()[1]))
        assert max(library_fingerprints) - min(library_fingerprints) < 1.0e-12, library_fingerprints
print("PASS Wannier90 MLWF adapter validation on 1, 2, 4, and 8 ranks")
