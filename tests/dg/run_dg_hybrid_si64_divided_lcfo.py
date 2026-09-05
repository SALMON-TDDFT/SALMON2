#!/usr/bin/env python3
"""Run Si64 divided WF+PW SCF followed by exactly one distributed LCFO solve."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path


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
    parser.add_argument("--binary", type=Path, default=root / "build-hybrid-commit/salmon")
    parser.add_argument("--result-dir", type=Path, default=root / "verification-si64-divided-lcfo")
    args = parser.parse_args()
    if args.mpi_ranks != 8:
        raise RuntimeError("Si64 divided LCFO validation requires exactly 8 MPI ranks")
    binary = args.binary.resolve(strict=True)
    result_dir = args.result_dir.resolve()
    if result_dir.exists():
        raise RuntimeError(f"result directory must be fresh: {result_dir}")
    fixture = root / "tests/dg/data/si64_overlapping_wannier_rt"
    result_dir.mkdir(parents=True)
    shutil.copy2(fixture / "input_hybrid_divided_lcfo.in", result_dir / "inputfile")
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
    scf_log_position = text.rfind("[OW-GS] divided WF+PW SCF converged")
    lcfo_log_position = text.find("[OW-GS] divided WF+PW LCFO solved once")
    post_lcfo_text = text[lcfo_log_position:] if lcfo_log_position >= 0 else ""
    evidence = {
        "return_code": return_code,
        "mpi_ranks": 8,
        "omp_threads": 1,
        "elapsed_seconds": time.monotonic() - started,
        "peak_rss_kib_by_process": sorted(peak_rss.values(), reverse=True),
        "conventional_dc_converged": "#GS converged at" in text,
        "wannier90_accepted": "overlapping_wannier_mlwf.chk" in " ".join(p.name for p in result_dir.iterdir()),
        "divided_scf": divided[-1] if divided else None,
        "final_lcfo": final[-1] if final else None,
        "final_lcfo_count": text.count("[OW-GS] divided WF+PW LCFO solved once"),
        "lcfo_after_divided_scf": 0 <= scf_log_position < lcfo_log_position,
        "post_lcfo_density_updates": post_lcfo_text.count("[DG-HYBRID-DIVIDED-SCF] iteration="),
        "occupied_checkpoint": (result_dir / "overlapping_wannier_occupied.chk").is_file(),
    }
    (result_dir / "divided_lcfo_evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    if return_code != 0:
        raise RuntimeError(f"Si64 divided LCFO failed with exit {return_code}; see {log_path}")
    if not evidence["conventional_dc_converged"] or not evidence["wannier90_accepted"]:
        raise RuntimeError("Si64 run ended without accepted conventional DC/Wannier90 seed")
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
