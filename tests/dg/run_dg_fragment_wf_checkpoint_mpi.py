#!/usr/bin/env python3
"""Build and exercise strict rank-local fragment-WF checkpoints."""

from pathlib import Path
import os
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
SCENARIOS = (
    "roundtrip", "rank_count_mismatch", "rank_fragment_permutation",
    "dc_seed_mismatch", "dc_seed_fingerprint_mismatch", "grid_mismatch", "cell_mismatch",
    "fragment_geometry_mismatch", "pseudopotential_mismatch", "boundary_mismatch",
    "inventory_mismatch", "ordering_mismatch", "selection_mismatch", "local_layout_mismatch",
    "generation_mismatch", "gauge_version_mismatch", "gauge_mode_mismatch",
    "gauge_fingerprint_mismatch", "missing_manifest", "missing_peer",
    "truncated_payload", "corrupt_payload", "unknown_version", "corrupt_manifest",
    "interrupted_publication", "interrupted_update_preserves_committed",
    "mode_disagreement", "auto_invalid", "status_disagreement",
)

with tempfile.TemporaryDirectory(prefix="dg-fragment-wf-checkpoint-") as name:
    build = Path(name)
    (build / "config.h").write_text("#define USE_MPI\n")
    executable = build / "test_dg_fragment_wf_checkpoint"
    compile_result = subprocess.run([
        shutil.which("mpifort"), "-cpp", "-I", str(build), "-J", str(build),
        "-O0", "-g", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-fbacktrace",
        str(ROOT / "src/gs/dc/dg_fragment_wf_checkpoint.f90"),
        str(ROOT / "tests/dg/test_dg_fragment_wf_checkpoint_mpi.f90"),
        "-o", str(executable),
    ], capture_output=True, text=True)
    if compile_result.returncode:
        raise RuntimeError("fragment-WF checkpoint compile failed:\n" + compile_result.stdout + compile_result.stderr)
    environment = os.environ.copy()
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for ranks in (1, 2, 4, 8):
        for scenario in SCENARIOS:
            case = build / f"ranks-{ranks}" / scenario
            checkpoint = case / "checkpoint"
            checkpoint.mkdir(parents=True)
            environment["DG_FRAGMENT_WF_DIRECTORY"] = str(checkpoint)
            environment["DG_FRAGMENT_WF_SCENARIO"] = scenario
            result = subprocess.run(
                [shutil.which("mpiexec"), "-n", str(ranks), str(executable)],
                cwd=case, env=environment, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, timeout=30,
            )
            if result.returncode:
                raise RuntimeError(f"{ranks}-rank {scenario} failed:\n{result.stdout}")
            marker = f"PASS fragment-WF checkpoint scenario={scenario} ranks={ranks}"
            if marker not in result.stdout:
                raise RuntimeError(f"missing PASS marker for {ranks} {scenario}:\n{result.stdout}")
    persisted = build / "persisted-rank-count-mismatch"
    checkpoint = persisted / "checkpoint"
    checkpoint.mkdir(parents=True)
    environment["DG_FRAGMENT_WF_DIRECTORY"] = str(checkpoint)
    environment["DG_FRAGMENT_WF_SCENARIO"] = "write_only"
    written = subprocess.run([shutil.which("mpiexec"), "-n", "2", str(executable)],
        cwd=persisted, env=environment, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, timeout=30)
    if written.returncode:
        raise RuntimeError("persisted 2-rank write failed:\n" + written.stdout)
    environment["DG_FRAGMENT_WF_SCENARIO"] = "read_existing_rank_mismatch"
    rejected = subprocess.run([shutil.which("mpiexec"), "-n", "4", str(executable)],
        cwd=persisted, env=environment, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, timeout=30)
    if rejected.returncode:
        raise RuntimeError("persisted rank-count mismatch was not rejected:\n" + rejected.stdout)
print("PASS strict fragment-WF checkpoints on 1, 2, 4, and 8 ranks")
