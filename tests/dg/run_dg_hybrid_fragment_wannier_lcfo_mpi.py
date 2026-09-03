#!/usr/bin/env python3
from pathlib import Path
import os
import re
import shlex
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests/dg/test_dg_hybrid_fragment_wannier_lcfo_mpi.f90"
SOURCES = [
    ROOT / "src/common/dg_hybrid_windowed_pw_types.f90",
    ROOT / "src/common/dg_hybrid_windowed_pw_basis.f90",
    ROOT / "src/common/dg_hybrid_wannier_complement.f90",
    ROOT / "src/gs/dc/dg_hybrid_fragment_basis.f90",
    ROOT / "src/gs/dc/dg_hybrid_fragment_basis_stream.f90",
    ROOT / "src/gs/dc/dg_hybrid_projected_fragment_pipeline.f90",
]
EXPECTED_GENERALIZED_SYMBOLS = (
    "compute_dg_hybrid_generalized_wannier_projection_tile",
    "build_dg_hybrid_complete_union_map",
    "s_dg_hybrid_dual_basis_catalog",
    "finalize_dg_hybrid_dual_basis_catalog",
    "expected_fragment_wannier_ranks",
)


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
assert compiler, "mpifort is required for the generalized fragment-union fixture"
assert launcher, "mpiexec is required for the generalized fragment-union fixture"
assert FIXTURE.is_file(), f"missing generalized fragment-union fixture: {FIXTURE}"

common_flags = [
    "-cpp",
    "-DUSE_MPI",
    "-std=f2008",
    "-ffree-line-length-none",
    "-fcheck=all",
    "-ffpe-trap=invalid,zero,overflow",
    "-fbacktrace",
]

with tempfile.TemporaryDirectory(prefix="hybrid-fragment-wannier-lcfo-") as name:
    build = Path(name)
    syntax_build = build / "syntax"
    production_build = build / "production"
    syntax_build.mkdir()
    production_build.mkdir()
    (syntax_build / "config.h").write_text("")
    (production_build / "config.h").write_text("")

    # Compile the fixture against its declaration-only wished-for contract first.
    # This separates a fixture syntax/interface defect from the intentional RED
    # caused by production not exporting the generalized APIs yet.
    fragment_object = syntax_build / "fragment_basis.o"
    subprocess.run(
        [
            compiler,
            *common_flags,
            "-I",
            str(syntax_build),
            "-J",
            str(syntax_build),
            "-c",
            str(ROOT / "src/gs/dc/dg_hybrid_fragment_basis.f90"),
            "-o",
            str(fragment_object),
        ],
        check=True,
    )
    syntax_result = subprocess.run(
        [
            compiler,
            *common_flags,
            "-DDG_GENERALIZED_CONTRACT_SYNTAX",
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
    if syntax_result.returncode != 0:
        raise RuntimeError(
            "generalized fragment-union fixture contract/syntax preflight failed:\n"
            + syntax_result.stdout
            + syntax_result.stderr
        )

    executable = production_build / "hybrid_fragment_wannier_lcfo"
    compile_result = subprocess.run(
        [
            compiler,
            *common_flags,
            "-I",
            str(production_build),
            "-J",
            str(production_build),
            *(str(source) for source in SOURCES),
            str(FIXTURE),
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
        lowered = diagnostic.lower()
        matched = [name for name in EXPECTED_GENERALIZED_SYMBOLS if name in lowered]
        assert matched, (
            "compile failed without naming a wished-for generalized API/type:\n" + diagnostic
        )
        forbidden = ("syntax error", "unexpected end", "unterminated", "invalid character")
        assert not any(token in lowered for token in forbidden), (
            "fixture syntax/build defect obscured the intended RED:\n" + diagnostic
        )
        raise RuntimeError(
            "CONTRACT/SYNTAX PREFLIGHT PASSED; EXPECTED RED: production is missing "
            "generalized fragment-union API/type "
            + ", ".join(matched)
            + "\n"
            + diagnostic
        )

    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    generalized_fingerprints = set()
    map_fingerprints = set()
    catalog_fingerprints = {}
    receipt_pattern = re.compile(
        r"HYBRID_FRAGMENT_WANNIER_LCFO ranks=(\d+) metric_rank=(\d+) "
        r"generalized_fingerprint=(-?\d+) map_fingerprint=(-?\d+) "
        r"catalog_fingerprint=(-?\d+)"
    )
    for rank_count in (1, 2, 4, 8):
        run = subprocess.run(
            [launcher, "-n", str(rank_count), str(executable)],
            capture_output=True,
            text=True,
            env=environment,
            timeout=60,
        )
        assert run.returncode == 0, (rank_count, run.stdout, run.stderr)
        assert (
            f"PASS hybrid fragment Wannier LCFO on {rank_count} ranks" in run.stdout
        ), run.stdout
        match = receipt_pattern.search(run.stdout)
        assert match, (rank_count, run.stdout)
        assert int(match.group(1)) == rank_count
        assert int(match.group(2)) == 5
        assert all(int(match.group(i)) != 0 for i in (3, 4, 5))
        generalized_fingerprints.add(int(match.group(3)))
        map_fingerprints.add(int(match.group(4)))
        catalog_fingerprints[rank_count] = int(match.group(5))
    assert len(generalized_fingerprints) == 1, generalized_fingerprints
    assert len(map_fingerprints) == 1, map_fingerprints
    # Fragment publishers remain ranks 0 and 1 for 2/4/8, but the exact
    # communicator topology is part of reuse provenance.
    assert len({catalog_fingerprints[ranks] for ranks in (2, 4, 8)}) == 3, catalog_fingerprints
    assert catalog_fingerprints[1] != catalog_fingerprints[2], catalog_fingerprints

print("PASS hybrid fragment Wannier LCFO on 1, 2, 4, and 8 ranks")
