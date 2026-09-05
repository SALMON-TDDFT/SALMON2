#!/usr/bin/env python3
"""Run Si64 divided WF+PW SCF followed by exactly one distributed LCFO solve."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import struct
import subprocess
import time
from pathlib import Path


EXPECTED_PUBLICATION_ID = 7047888166118007469
EXPECTED_MPI_SIZE = 8
EXPECTED_MAPPING_FINGERPRINT = 254086644876463474
DEFAULT_SEED_DIRECTORY = Path("/tmp/si64-task8-bounded-smoke-20260905/dc-seed")


def validate_seed_directory(seed_directory: Path) -> dict[str, int]:
    """Reject the comparison before MPI startup unless the exact seed is present."""
    manifest = seed_directory / "dg_dc_seed.manifest"
    if not manifest.is_file():
        raise RuntimeError(f"missing reusable seed manifest: {manifest}")
    manifest_header = manifest.read_bytes()[:48]
    if len(manifest_header) != 48:
        raise RuntimeError("truncated reusable seed manifest")
    magic, version, mpi_size, publication_id = struct.unpack("<32siiq", manifest_header)
    if magic.rstrip() != b"SALMON_DG_DC_SEED_MANIFEST_V1":
        raise RuntimeError("unexpected reusable seed manifest magic")
    if (version, mpi_size, publication_id) != (
        1,
        EXPECTED_MPI_SIZE,
        EXPECTED_PUBLICATION_ID,
    ):
        raise RuntimeError("reusable seed manifest identity changed")
    publication_hex = f"{publication_id:016X}"
    mapping_fingerprint = None
    for rank in range(EXPECTED_MPI_SIZE):
        shard = seed_directory / f"dg_dc_seed.{publication_hex}.rank{rank:08d}.shard"
        header = shard.read_bytes()[:80] if shard.is_file() else b""
        if len(header) != 80:
            raise RuntimeError(f"missing or truncated reusable seed shard for rank {rank}")
        values = struct.unpack("<32siiiiqqqq", header)
        shard_magic, shard_version, shard_mpi, shard_rank, fragment_id = values[:5]
        shard_publication, _, ownership_fingerprint, _ = values[5:]
        if (
            shard_magic.rstrip() != b"SALMON_DG_DC_SEED_SHARD_V1"
            or (shard_version, shard_mpi, shard_rank, fragment_id)
            != (1, EXPECTED_MPI_SIZE, rank, rank + 1)
            or shard_publication != publication_id
        ):
            raise RuntimeError(f"reusable seed rank-fragment contract changed at rank {rank}")
        if rank == 0:
            mapping_fingerprint = ownership_fingerprint
    if mapping_fingerprint != EXPECTED_MAPPING_FINGERPRINT:
        raise RuntimeError("reusable seed mapping fingerprint changed")
    return {
        "publication_id": publication_id,
        "mpi_size": mpi_size,
        "mapping_fingerprint": mapping_fingerprint,
    }


def render_input(template: str, mixing_method: str, seed_directory: Path) -> str:
    selector = re.compile(
        r"^(?P<prefix>\s*dg_hybrid_divided_mixing\s*=\s*)['\"](?:simple|pulay)['\"]\s*$",
        re.IGNORECASE | re.MULTILINE,
    )
    if len(selector.findall(template)) != 1:
        raise RuntimeError("expected exactly one Hybrid divided mixing selector")
    rendered = selector.sub(lambda match: match.group("prefix") + f"'{mixing_method}'", template)
    if re.search(r"\bdg_dc_seed_(?:mode|directory)\s*=", rendered, re.IGNORECASE):
        raise RuntimeError("Si64 divided fixture unexpectedly contains a seed override")
    dc_start = rendered.lower().index("&dc")
    dc_end = rendered.index("/", dc_start)
    seed_controls = (
        " dg_dc_seed_mode='read'\n"
        f" dg_dc_seed_directory='{seed_directory}'\n"
    )
    return rendered[:dc_end] + seed_controls + rendered[dc_end:]


def seed_receipt(text: str) -> dict[str, int]:
    matches = re.findall(
        r"\[DG-DC-SEED\]\s+mode=(\w+)\s+publication_id=(-?\d+)\s+"
        r"scf_skipped=([TF])\s+mpi_size=(\d+)\s+mapping_fingerprint=(-?\d+)",
        text,
    )
    if len(matches) != 1:
        raise RuntimeError("expected exactly one reusable DC seed receipt")
    mode, publication, skipped, mpi_size, mapping = matches[0]
    receipt = {
        "publication_id": int(publication),
        "mpi_size": int(mpi_size),
        "mapping_fingerprint": int(mapping),
    }
    if mode != "read" or skipped != "T" or receipt != {
        "publication_id": EXPECTED_PUBLICATION_ID,
        "mpi_size": EXPECTED_MPI_SIZE,
        "mapping_fingerprint": EXPECTED_MAPPING_FINGERPRINT,
    }:
        raise RuntimeError("runtime DC seed receipt changed or ordinary DC was not skipped")
    if "DC #SCF =" in text:
        raise RuntimeError("ordinary DC SCF executed during a seed-reuse comparison")
    return receipt


def process_rows() -> list[tuple[int, int, int, str]]:
    result = subprocess.run(
        ["ps", "-axo", "pid=,ppid=,rss=,command="], check=True, capture_output=True, text=True
    )
    rows = []
    for line in result.stdout.splitlines():
        fields = line.strip().split(None, 3)
        if len(fields) == 4:
            rows.append((int(fields[0]), int(fields[1]), int(fields[2]), fields[3]))
    return rows


def descendants(root: int, rows: list[tuple[int, int, int, str]]) -> set[int]:
    selected = {root}
    changed = True
    while changed:
        changed = False
        for pid, parent, _, _ in rows:
            if parent in selected and pid not in selected:
                selected.add(pid)
                changed = True
    return selected


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--mpi-ranks", type=int, default=8)
    parser.add_argument("--binary", type=Path, default=root / "build-hybrid-release/salmon")
    parser.add_argument("--mixing", choices=("simple", "pulay"), required=True)
    parser.add_argument("--seed-directory", type=Path, default=DEFAULT_SEED_DIRECTORY)
    parser.add_argument("--result-dir", type=Path)
    args = parser.parse_args()
    if args.mpi_ranks != 8:
        raise RuntimeError("Si64 divided LCFO validation requires exactly 8 MPI ranks")
    binary = args.binary.resolve(strict=True)
    seed_directory = args.seed_directory.resolve(strict=True)
    expected_seed = validate_seed_directory(seed_directory)
    result_dir = (
        args.result_dir
        if args.result_dir is not None
        else Path("/tmp") / f"si64-task8-schwarz-{args.mixing}-{time.strftime('%Y%m%d-%H%M%S')}"
    ).resolve()
    if result_dir.exists():
        raise RuntimeError(f"result directory must be fresh: {result_dir}")
    fixture = root / "tests/dg/data/si64_overlapping_wannier_rt"
    result_dir.mkdir(parents=True)
    template = (fixture / "input_hybrid_divided_lcfo.in").read_text()
    rendered_input = render_input(template, args.mixing, seed_directory)
    (result_dir / "inputfile").write_text(rendered_input)
    shutil.copy2(fixture / "atom.dat", result_dir / "atom.dat")
    shutil.copy2(root / "samples/exercise_04_bulkSi_gs/Si_rps.dat", result_dir / "Si_rps.dat")

    log_path = result_dir / "run.log"
    peak_rss: dict[int, int] = {}
    started = time.monotonic()
    with (result_dir / "inputfile").open("rb") as source, log_path.open("wb") as log:
        process = subprocess.Popen(
            [shutil.which("mpirun") or "mpirun", "-np", "8", str(binary)],
            cwd=result_dir,
            stdin=source,
            stdout=log,
            stderr=subprocess.STDOUT,
            env={**os.environ, "OMP_NUM_THREADS": "1"},
        )
        while process.poll() is None:
            rows = process_rows()
            family = descendants(process.pid, rows)
            for pid, _, rss, command in rows:
                if pid in family and binary.name in command:
                    peak_rss[pid] = max(peak_rss.get(pid, 0), rss)
            time.sleep(2.0)
        return_code = process.returncode

    text = log_path.read_text(errors="replace")
    divided = re.findall(
        r"\[OW-GS\] divided WF\+PW SCF converged iterations=(\d+)\s+density=\s*([0-9Ee+\-.]+)", text
    )
    final = re.findall(
        r"\[OW-GS\] divided WF\+PW LCFO solved once\s+residual=\s*([0-9Ee+\-.]+)\s+"
        r"orthogonality=\s*([0-9Ee+\-.]+)\s+projector=\s*([0-9Ee+\-.]+)", text
    )
    schwarz = re.findall(
        r"\[DG-HYBRID-SCHWARZ\]\s+epoch=(\d+)\s+neighbor_exchanges=(\d+)\s+"
        r"accepted_cg_steps=(\d+)\s+common_extensions=(\d+)\s+residual=\s*([0-9Ee+\-.]+)\s+"
        r"electron_defect=\s*([0-9Ee+\-.]+)\s+temperature=\s*([0-9Ee+\-.]+)",
        text,
    )
    convergence_history = re.findall(
        r"\[DG-HYBRID-DIVIDED-SCF\]\s+iteration=(\d+)\s+convergence=\s*([0-9Ee+\-.]+)\s+"
        r"electron_defect=\s*([0-9Ee+\-.]+)",
        text,
    )
    runtime_seed = seed_receipt(text)
    cg_steps = [int(row[2]) for row in schwarz]
    neighbor_exchanges = [int(row[1]) for row in schwarz]
    has_invalid_number = re.search(r"(?<![A-Za-z])(?:nan|[+\-]?inf)(?![A-Za-z])", text, re.IGNORECASE) is not None
    has_rollback_or_fallback = (
        re.search(r"\b(?:rolled_back|rollback|fallback)\b", text, re.IGNORECASE) is not None
    )
    scf_log_position = text.rfind("[OW-GS] divided WF+PW SCF converged")
    lcfo_log_position = text.find("[OW-GS] divided WF+PW LCFO solved once")
    post_lcfo_text = text[lcfo_log_position:] if lcfo_log_position >= 0 else ""
    evidence = {
        "return_code": return_code,
        "mpi_ranks": 8,
        "omp_threads": 1,
        "mixing": args.mixing,
        "elapsed_seconds": time.monotonic() - started,
        "peak_rss_kib_by_process": sorted(peak_rss.values(), reverse=True),
        "expected_seed": expected_seed,
        "runtime_seed": runtime_seed,
        "ordinary_dc_scf_executed": "DC #SCF =" in text,
        "wannier90_accepted": "overlapping_wannier_mlwf.chk" in " ".join(p.name for p in result_dir.iterdir()),
        "schwarz_epochs": schwarz,
        "density_convergence_history": convergence_history,
        "full_neighbor_dg_from_epoch_1": bool(
            schwarz and int(schwarz[0][0]) == 1 and int(schwarz[0][1]) > 0
        ),
        "cg_steps_within_1_to_3": bool(cg_steps) and all(1 <= step <= 3 for step in cg_steps),
        "neighbor_exchange_minimum": min(neighbor_exchanges) if neighbor_exchanges else None,
        "invalid_number_seen": has_invalid_number,
        "rollback_or_fallback_seen": has_rollback_or_fallback,
        "repeated_full_diagonalization": len(final) > 1,
        "divided_scf": divided[-1] if divided else None,
        "final_lcfo": final[-1] if final else None,
        "final_lcfo_count": text.count("[OW-GS] divided WF+PW LCFO solved once"),
        "lcfo_after_divided_scf": 0 <= scf_log_position < lcfo_log_position,
        "post_lcfo_density_updates": post_lcfo_text.count("[DG-HYBRID-DIVIDED-SCF] iteration="),
        "occupied_checkpoint": (result_dir / "overlapping_wannier_occupied.chk").is_file(),
    }
    (result_dir / "divided_lcfo_evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    if evidence["expected_seed"] != evidence["runtime_seed"] or evidence["ordinary_dc_scf_executed"]:
        raise RuntimeError("Si64 run did not reuse the exact authoritative DC seed")
    if (
        not evidence["full_neighbor_dg_from_epoch_1"]
        or not evidence["cg_steps_within_1_to_3"]
        or evidence["invalid_number_seen"]
        or evidence["rollback_or_fallback_seen"]
        or evidence["repeated_full_diagonalization"]
        or evidence["post_lcfo_density_updates"] != 0
    ):
        raise RuntimeError("Si64 Schwarz invariants failed; see divided_lcfo_evidence.json")
    if divided and (len(final) != 1 or not evidence["lcfo_after_divided_scf"]):
        raise RuntimeError("converged divided SCF did not perform exactly one terminal LCFO solve")
    if return_code != 0:
        raise RuntimeError(f"Si64 divided LCFO failed with exit {return_code}; see {log_path}")
    if (
        not divided
        or len(final) != 1
        or evidence["final_lcfo_count"] != 1
        or not evidence["lcfo_after_divided_scf"]
        or evidence["post_lcfo_density_updates"] != 0
    ):
        raise RuntimeError("Si64 run did not complete one divided SCF followed by exactly one LCFO solve")
    if not evidence["occupied_checkpoint"]:
        raise RuntimeError("Si64 divided LCFO checkpoint was not published")
    print(result_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
