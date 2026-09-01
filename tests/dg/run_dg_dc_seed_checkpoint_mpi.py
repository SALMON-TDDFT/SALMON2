#!/usr/bin/env python3
"""Build and exercise strict rank-sharded DC seed checkpoints."""

from pathlib import Path
import os
import re
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
SCENARIOS = (
    "absent",
    "roundtrip",
    "missing_manifest",
    "missing_shard",
    "truncated_shard",
    "corrupt_manifest",
    "corrupt_shard",
    "interrupted_publication",
    "interrupted_update_preserves_committed",
    "mixed_publication_read",
    "rank_directory_write_mismatch",
    "rank_count_mismatch",
    "rank_fragment_mismatch",
    "local_bound_mismatch",
    "ownership_map_mismatch",
    "immutable_mismatch",
    "electron_count_mismatch",
    "residual_threshold_mismatch",
)


def assert_bounded_memory_validation_routes() -> None:
    source = (ROOT / "src/gs/dc/dg_dc_seed_checkpoint.f90").read_text()
    write_body = re.search(
        r"subroutine write_dg_dc_seed\b.*?end subroutine write_dg_dc_seed",
        source,
        re.DOTALL,
    )
    probe_body = re.search(
        r"subroutine probe_dg_dc_seed\b.*?end subroutine probe_dg_dc_seed",
        source,
        re.DOTALL,
    )
    assert write_body and probe_body
    assert "verified_payload" not in write_body.group(0), (
        "writer validation must not allocate a second full seed payload"
    )
    assert "read_shard_file" not in write_body.group(0), (
        "writer validation must use bounded-memory shard streaming"
    )
    assert "type(s_dg_dc_seed_payload)::payload" not in probe_body.group(0), (
        "probe must not materialize the full rank-local seed payload"
    )


def main() -> None:
    assert_bounded_memory_validation_routes()
    compiler = shutil.which("mpifort")
    launcher = shutil.which("mpiexec")
    if compiler is None or launcher is None:
        raise SystemExit("mpifort and mpiexec are required")

    with tempfile.TemporaryDirectory(prefix="dg-dc-seed-checkpoint-") as name:
        build = Path(name)
        (build / "config.h").write_text("")
        executable = build / "test_dg_dc_seed_checkpoint"
        subprocess.run(
            [
                compiler,
                "-cpp",
                "-DUSE_MPI",
                "-I",
                str(build),
                "-J",
                str(build),
                "-fcheck=all",
                "-ffpe-trap=invalid,zero,overflow",
                "-fbacktrace",
                str(ROOT / "src/gs/dc/dg_dc_seed_checkpoint.f90"),
                str(ROOT / "tests/dg/test_dg_dc_seed_checkpoint_mpi.f90"),
                "-o",
                str(executable),
            ],
            check=True,
        )

        environment = os.environ.copy()
        environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
        for ranks in (1, 2, 4, 8):
            for scenario in SCENARIOS:
                case = build / f"ranks-{ranks}" / scenario
                seed = case / "seed"
                seed.mkdir(parents=True)
                if scenario in (
                    "mixed_publication_read",
                    "rank_directory_write_mismatch",
                ):
                    (seed / "a").mkdir()
                    (seed / "b").mkdir()
                environment["DG_DC_SEED_DIRECTORY"] = str(seed)
                environment["DG_DC_SEED_SCENARIO"] = scenario
                completed = subprocess.run(
                    [launcher, "-n", str(ranks), str(executable)],
                    cwd=case,
                    env=environment,
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    timeout=30,
                )
                assert completed.returncode == 0, (
                    ranks,
                    scenario,
                    completed.stdout,
                )
                marker = f"PASS DG DC seed scenario={scenario} ranks={ranks}"
                assert marker in completed.stdout, completed.stdout
                match = re.search(r"publication_id=(-?\d+)", completed.stdout)
                assert match, completed.stdout
                failed_write_may_not_reserve_an_id = scenario == (
                    "interrupted_update_preserves_committed"
                ) or (scenario == "rank_directory_write_mismatch" and ranks > 1)
                if scenario != "absent" and not failed_write_may_not_reserve_an_id:
                    assert int(match.group(1)) != 0, completed.stdout

        case = build / "persisted-rank-count-mismatch"
        seed = case / "seed"
        seed.mkdir(parents=True)
        environment["DG_DC_SEED_DIRECTORY"] = str(seed)
        environment["DG_DC_SEED_SCENARIO"] = "write_only"
        written = subprocess.run(
            [launcher, "-n", "2", str(executable)],
            cwd=case,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=30,
        )
        assert written.returncode == 0, written.stdout
        assert "PASS DG DC seed scenario=write_only ranks=2" in written.stdout

        environment["DG_DC_SEED_SCENARIO"] = "read_existing_rank_mismatch"
        rejected = subprocess.run(
            [launcher, "-n", "4", str(executable)],
            cwd=case,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=30,
        )
        assert rejected.returncode == 0, rejected.stdout
        assert (
            "PASS DG DC seed scenario=read_existing_rank_mismatch ranks=4"
            in rejected.stdout
        )

    print("PASS strict DG DC seed checkpoints on 1, 2, 4, and 8 ranks")


if __name__ == "__main__":
    main()
