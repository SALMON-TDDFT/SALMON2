#!/usr/bin/env python3
"""Selection geometry plus real raw-cache validation, with only W90 stubbed."""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
compiler = shutil.which("mpifort")
launcher = shutil.which("mpiexec")
assert compiler and launcher, "MPI compiler and launcher required"
if os.environ.get("SALMON_LAPACK_LIBS"):
    libraries = shlex.split(os.environ["SALMON_LAPACK_LIBS"])
elif shutil.which("pkg-config") and subprocess.run(
        ["pkg-config", "--exists", "openblas"]).returncode == 0:
    libraries = shlex.split(subprocess.check_output(
        ["pkg-config", "--libs", "openblas"], text=True))
elif shutil.which("brew"):
    prefix = subprocess.check_output(["brew", "--prefix", "openblas"], text=True).strip()
    libraries = [f"-L{prefix}/lib", "-lopenblas"]
else:
    libraries = ["-llapack", "-lblas"]

with tempfile.TemporaryDirectory(prefix="hybrid-selection-") as temporary:
    build = Path(temporary)
    (build / "config.h").write_text("\n".join(
        f"#define {name}" for name in ("SYSTEM_HAS_POSIX", "SYSTEM_HAS_POSIX_STAT",
        "SYSTEM_HAS_POSIX_ACCESS", "SYSTEM_HAS_POSIX_MKDIR", "SYSTEM_HAS_POSIX_NFTW")))
    sources = ["src/io/posix.c", "src/gs/dc/lcfo_wannier_sawf_seed.f90",
               "src/gs/dc/dg_overlapping_wannier_w90.f90",
               "src/gs/dc/dg_hybrid_fragment_wannier.f90"]
    selection = "src/gs/dc/dg_hybrid_fragment_selection.f90"
    if (ROOT / selection).exists():
        sources.append(selection)
    sources += ["src/math/phys_constants.f90",
                "src/gs/occupation_kernel.f90",
                "src/gs/dc/dc_fragment_occupation.f90",
                "src/common/dg_hybrid_windowed_pw_types.f90",
                "src/common/dg_hybrid_reciprocal_catalog.f90",
                "src/common/dg_hybrid_windowed_pw_basis.f90",
                "src/common/dg_hybrid_wannier_complement.f90",
                "src/gs/dc/dg_hybrid_fragment_basis.f90",
                "src/gs/dc/dg_hybrid_fragment_basis_stream.f90",
                "src/gs/dc/dg_hybrid_broken_volume.f90",
                "src/gs/dc/dg_hybrid_sipg_operator.f90",
                "src/gs/dc/dg_overlapping_wannier_nonlocal.f90",
                "src/gs/dc/dg_hybrid_variational_payload.f90",
                "src/gs/dc/dg_hybrid_divided_operator.f90",
                "src/gs/dc/dg_hybrid_fragment_subspace.f90",
                "src/gs/dc/dg_hybrid_fragment_solver.f90",
                "src/gs/dc/dg_hybrid_fragment_preconditioner.f90",
                "src/gs/dc/dg_hybrid_projected_fragment_pipeline.f90",
                "src/gs/dc/dg_hybrid_fragment_admission.f90",
                "src/gs/dc/dg_hybrid_fragment_thermal.f90"]
    sources += ["tests/dg/test_dg_hybrid_fragment_wannier_mpi.f90",
                "tests/dg/test_dg_hybrid_fragment_selection_mpi.f90"]
    executable = build / "selection"
    result = subprocess.run([compiler, "-cpp", "-DUSE_MPI", "-DUSE_WANNIER90",
        "-DW90_TEST_STUBS", "-DDG_W90_STUBS_ONLY", "-ffree-line-length-none",
        "-g", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-fbacktrace",
        "-I", str(build), "-J", str(build), *(str(ROOT / p) for p in sources),
        *libraries, "-o", str(executable)], capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for count in (1, 2, 4, 8):
        run_dir = build / f"ranks-{count}"
        run_dir.mkdir()
        result = subprocess.run([launcher, "-n", str(count), str(executable)],
            cwd=run_dir, env=environment, capture_output=True, text=True, timeout=60)
        assert result.returncode == 0, (count, result.stdout, result.stderr)
        assert f"PASS core-center selection on {count} ranks" in result.stdout, result.stdout
        if count > 1:
            assert f"PASS thermal DG handoff on {count} ranks" in result.stdout, result.stdout
            print(next(line for line in result.stdout.splitlines() if line.startswith("PASS thermal DG")))
            failure_marker = f"PASS thermal DG capacity/tail failures on {count} ranks"
            assert failure_marker in result.stdout, result.stdout
            print(failure_marker)
            transaction_marker = f"PASS thermal DG transaction on {count} ranks"
            assert transaction_marker in result.stdout, result.stdout
            print(transaction_marker)
            assert f"PASS explicit cutoff projection on {count} ranks" in result.stdout, result.stdout
            print(next(line for line in result.stdout.splitlines() if line.startswith("PASS explicit cutoff")))
print("PASS core-center selection on 1, 2, 4, and 8 ranks")
