#!/usr/bin/env python3
from pathlib import Path
import argparse
import os
import shlex
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
FRAGMENT_SOURCE = ROOT / "src/gs/dc/dg_hybrid_fragment_wannier.f90"
TEST_SOURCE = ROOT / "tests/dg/test_dg_hybrid_fragment_wannier_mpi.f90"

parser = argparse.ArgumentParser()
parser.add_argument("--core-metric-audit", action="store_true",
                    help="Exercise the pending DC-to-CG handoff with the production core-only metric")
arguments = parser.parse_args()


def lapack_libraries():
    configured = os.environ.get("SALMON_LAPACK_LIBS")
    if configured:
        return shlex.split(configured)
    if shutil.which("pkg-config"):
        openblas = subprocess.run(
            ["pkg-config", "--libs", "openblas"], capture_output=True, text=True
        )
        if openblas.returncode == 0:
            return shlex.split(openblas.stdout)
        lapack = subprocess.run(
            ["pkg-config", "--libs", "lapack"], capture_output=True, text=True
        )
        if lapack.returncode == 0:
            return shlex.split(lapack.stdout)
    if shutil.which("brew"):
        openblas = subprocess.run(
            ["brew", "--prefix", "openblas"], capture_output=True, text=True
        )
        if openblas.returncode == 0:
            return [f"-L{openblas.stdout.strip()}/lib", "-lopenblas"]
    return ["-llapack", "-lblas"]


with tempfile.TemporaryDirectory(prefix="hybrid-fragment-wannier-") as name:
    assert TEST_SOURCE.is_file(), f"missing fragment-Wannier fixture: {TEST_SOURCE}"
    build = Path(name)
    (build / "config.h").write_text(
        "#define SYSTEM_HAS_POSIX\n"
        "#define SYSTEM_HAS_POSIX_STAT\n"
        "#define SYSTEM_HAS_POSIX_ACCESS\n"
        "#define SYSTEM_HAS_POSIX_MKDIR\n"
        "#define SYSTEM_HAS_POSIX_NFTW\n"
    )
    executable = build / "hybrid_fragment_wannier"
    sources = [
        ROOT / "src/io/posix.c",
        ROOT / "src/gs/dc/lcfo_wannier_sawf_seed.f90",
        ROOT / "src/gs/dc/dg_overlapping_wannier_w90.f90",
        ROOT / "src/gs/dc/dg_fragment_scdm_gauge.f90",
        ROOT / "src/gs/dc/dg_fragment_wf_checkpoint.f90",
        ROOT / "src/gs/dc/dg_hybrid_fragment_subspace.f90",
        ROOT / "src/common/dg_hybrid_windowed_pw_types.f90",
        ROOT / "src/common/dg_hybrid_windowed_pw_basis.f90",
        ROOT / "src/common/dg_hybrid_wannier_complement.f90",
        ROOT / "src/gs/dc/dg_hybrid_fragment_basis.f90",
        ROOT / "src/gs/dc/dg_hybrid_fragment_basis_stream.f90",
        ROOT / "src/gs/dc/dg_hybrid_projected_fragment_pipeline.f90",
        ROOT / "src/gs/dc/dg_hybrid_variational_payload.f90",
        ROOT / "src/gs/dc/dg_hybrid_broken_volume.f90",
        ROOT / "src/gs/dc/dg_hybrid_divided_operator.f90",
        ROOT / "src/gs/dc/dg_hybrid_fragment_preconditioner.f90",
        ROOT / "src/gs/dc/dg_hybrid_fragment_solver.f90",
    ]
    # During RED the wished-for module is absent, so compiling the fixture itself
    # gives the useful "cannot open ...mod" diagnostic.  Once implemented, the
    # exact same permanent runner compiles the production source before the test.
    if FRAGMENT_SOURCE.exists():
        sources.append(FRAGMENT_SOURCE)
    sources.append(TEST_SOURCE)
    compiler = shutil.which("mpifort")
    launcher = shutil.which("mpiexec")
    assert compiler, "mpifort is required for the fragment-Wannier MPI fixture"
    assert launcher, "mpiexec is required for the fragment-Wannier MPI fixture"
    compile_result = subprocess.run(
        [
            compiler,
            "-cpp",
            "-DUSE_MPI",
            "-DUSE_WANNIER90",
            "-DW90_TEST_STUBS",
            "-std=f2008",
            "-fall-intrinsics",
            "-ffree-line-length-none",
            "-I",
            str(build),
            "-J",
            str(build),
            "-fcheck=all",
            "-ffpe-trap=invalid,zero,overflow",
            "-fbacktrace",
            *(str(source) for source in sources),
            *lapack_libraries(),
            "-o",
            str(executable),
        ],
        capture_output=True,
        text=True,
        timeout=120,
    )
    if compile_result.returncode != 0:
        diagnostic = compile_result.stdout + compile_result.stderr
        if not FRAGMENT_SOURCE.exists():
            assert "dg_hybrid_fragment_wannier.mod" in diagnostic.lower(), diagnostic
        raise RuntimeError("fragment-Wannier compile failed:\n" + diagnostic)

    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    cases = [(2, True)] if arguments.core_metric_audit else [(2, False), (4, False), (8, False), (2, True)]
    for rank_count, audit in cases:
        run_directory = build / (f"audit-ranks-{rank_count}" if audit else f"ranks-{rank_count}")
        run_directory.mkdir()
        run = subprocess.run(
            [launcher, "-n", str(rank_count), str(executable),
             *(["--core-metric-audit"] if audit else [])],
            cwd=run_directory,
            capture_output=True,
            text=True,
            env=environment,
            timeout=60,
        )
        assert run.returncode == 0, (rank_count, run.stdout, run.stderr)
        if audit:
            for fragment in (1, 2):
                assert run.stdout.count(f"PASS expected core-null rejection fragment={fragment}") == 1, run.stdout
        assert (
            f"PASS hybrid fragment Wannier on {rank_count} ranks" in run.stdout
        ), run.stdout
        artifact_root = run_directory / "fragment-wannier-artifacts"
        expected_generations = (7, 8, 9, 11, 12, 13, 27, 28, 29)
        expected_directories = {
            artifact_root / f"f{fragment_id:06d}" / f"g{generation:08d}"
            for fragment_id in (1, 2)
            for generation in expected_generations
        }
        actual_directories = set(artifact_root.glob("f*/g*"))
        assert actual_directories == expected_directories, actual_directories
        assert not list(artifact_root.rglob("*.dmn")), "unconstrained mode emitted .dmn"

print("PASS core-metric audit on 2 ranks" if arguments.core_metric_audit else
      "PASS hybrid fragment Wannier on 2, 4, and 8 ranks")
if not arguments.core_metric_audit:
    print("PASS expected core-null rejection on 2 ranks")
