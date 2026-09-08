#!/usr/bin/env python3
from pathlib import Path
import os
import shlex
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests/dg/test_dg_hybrid_projected_fragment_pipeline_mpi.f90"
DEPENDENCIES = [
    ROOT / "src/common/dg_hybrid_windowed_pw_types.f90",
    ROOT / "src/common/dg_hybrid_windowed_pw_basis.f90",
    ROOT / "src/common/dg_hybrid_wannier_complement.f90",
    ROOT / "src/gs/dc/dg_hybrid_fragment_basis.f90",
    ROOT / "src/gs/dc/dg_hybrid_fragment_basis_stream.f90",
]
PIPELINE = ROOT / "src/gs/dc/dg_hybrid_projected_fragment_pipeline.f90"
EXPECTED_CONTRACT = "s_dg_hybrid_projection_factorization_receipt"


def lapack_libraries():
    configured = os.environ.get("SALMON_LAPACK_LIBS")
    if configured:
        return shlex.split(configured)
    if shutil.which("pkg-config"):
        for package in ("openblas", "lapack"):
            probe = subprocess.run(
                ["pkg-config", "--libs", package], capture_output=True, text=True
            )
            if probe.returncode == 0:
                return shlex.split(probe.stdout)
    if shutil.which("brew"):
        probe = subprocess.run(
            ["brew", "--prefix", "openblas"], capture_output=True, text=True
        )
        if probe.returncode == 0:
            return [f"-L{probe.stdout.strip()}/lib", "-lopenblas"]
    return ["-llapack", "-lblas"]

compiler = shutil.which("mpifort")
launcher = shutil.which("mpiexec")
assert compiler, "mpifort is required for the projected-fragment pipeline fixture"
assert launcher, "mpiexec is required for the projected-fragment pipeline fixture"

flags = [
    "-cpp",
    "-DUSE_MPI",
    "-std=f2008",
    "-ffree-line-length-none",
    "-fcheck=all",
    "-ffpe-trap=invalid,zero,overflow",
    "-fbacktrace",
]

with tempfile.TemporaryDirectory(prefix="hybrid-projected-fragment-") as name:
    build = Path(name)
    syntax_build = build / "syntax"
    production_build = build / "production"
    syntax_build.mkdir()
    production_build.mkdir()
    (syntax_build / "config.h").write_text("")
    (production_build / "config.h").write_text("")

    for source in DEPENDENCIES:
        subprocess.run(
            [
                compiler,
                *flags,
                "-I",
                str(syntax_build),
                "-J",
                str(syntax_build),
                "-c",
                str(source),
                "-o",
                str(syntax_build / f"{source.stem}.o"),
            ],
            check=True,
        )
    syntax = subprocess.run(
        [
            compiler,
            *flags,
            "-DDG_GENERALIZED_PIPELINE_CONTRACT_SYNTAX",
            "-I",
            str(syntax_build),
            "-J",
            str(syntax_build),
            "-c",
            str(FIXTURE),
            "-o",
            str(syntax_build / "fixture.o"),
        ],
        capture_output=True,
        text=True,
    )
    if syntax.returncode != 0:
        raise RuntimeError(
            "projected-fragment factorization receipt syntax preflight failed:\n"
            + syntax.stdout
            + syntax.stderr
        )

    executable = production_build / "hybrid_projected_fragment"
    compile_result = subprocess.run(
        [
            compiler,
            *flags,
            "-I",
            str(production_build),
            "-J",
            str(production_build),
            *(str(source) for source in DEPENDENCIES),
            str(PIPELINE),
            str(FIXTURE),
            *lapack_libraries(),
            "-o",
            str(executable),
        ],
        capture_output=True,
        text=True,
    )
    if compile_result.returncode != 0:
        diagnostic = compile_result.stdout + compile_result.stderr
        assert EXPECTED_CONTRACT in diagnostic.lower(), (
            "pipeline compile failed without naming the wished-for generation receipt:\n"
            + diagnostic
        )
        raise RuntimeError(
            "CONTRACT/SYNTAX PREFLIGHT PASSED; EXPECTED RED: production is missing "
            + EXPECTED_CONTRACT
            + "\n"
            + diagnostic
        )

    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    # Four ranks exercise both active fragment publishers and idle ranks first,
    # so RED diagnostics cover publication semantics before shorter layouts.
    for rank_count in (4, 2, 8):
        run = subprocess.run(
            [launcher, "-n", str(rank_count), str(executable)],
            capture_output=True,
            text=True,
            env=environment,
            timeout=60,
        )
        assert run.returncode == 0, (rank_count, run.stdout, run.stderr)
        assert f"PASS hybrid projected fragment pipeline on {rank_count} ranks" in run.stdout

print("PASS hybrid projected fragment pipeline on 2, 4, and 8 ranks")
