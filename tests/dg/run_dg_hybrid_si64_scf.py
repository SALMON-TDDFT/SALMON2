#!/usr/bin/env python3
"""Run the default-off Si64 hybrid ScaLAPACK SCF route with RSS monitoring."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path


def process_table() -> list[tuple[int, int, int, str]]:
    result = subprocess.run(
        ["ps", "-axo", "pid=,ppid=,rss=,command="],
        check=True,
        capture_output=True,
        text=True,
    )
    rows: list[tuple[int, int, int, str]] = []
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
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("result_dir", type=Path)
    parser.add_argument("--ranks", type=int, default=8)
    parser.add_argument("--timeout-minutes", type=float, default=0.0)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    if args.ranks < 1:
        raise RuntimeError("MPI rank count must be positive")
    binary = args.binary.resolve(strict=True)
    result_dir = args.result_dir.resolve()
    if result_dir.exists():
        raise RuntimeError(f"result directory must be fresh: {result_dir}")
    fixtures = args.repo / "tests/dg/data/si64_overlapping_wannier_rt"
    result_dir.mkdir(parents=True)
    for name in ("atom.dat", "input_hybrid_scf.in"):
        shutil.copy2(fixtures / name, result_dir / ("inputfile" if name.startswith("input_") else name))
    shutil.copy2(args.repo / "samples/exercise_04_bulkSi_gs/Si_rps.dat", result_dir / "Si_rps.dat")

    log_path = result_dir / "run.log"
    peak_rss_kib: dict[int, int] = {}
    started = time.monotonic()
    with (result_dir / "inputfile").open("rb") as source, log_path.open("wb") as log:
        process = subprocess.Popen(
            ["mpirun", "-np", str(args.ranks), str(binary)],
            cwd=result_dir,
            stdin=source,
            stdout=log,
            stderr=subprocess.STDOUT,
            env={**os.environ, "OMP_NUM_THREADS": "1"},
        )
        while process.poll() is None:
            rows = process_table()
            family = descendants(process.pid, rows)
            for pid, _, rss, command in rows:
                if pid in family and (str(binary) in command or Path(binary).name in command):
                    peak_rss_kib[pid] = max(peak_rss_kib.get(pid, 0), rss)
            if args.timeout_minutes > 0 and time.monotonic() - started > 60.0 * args.timeout_minutes:
                process.terminate()
                try:
                    process.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    process.kill()
                raise RuntimeError("Si64 hybrid SCF exceeded the requested timeout")
            time.sleep(2.0)
        return_code = process.returncode

    log_text = log_path.read_text(errors="replace")
    convergence = re.findall(
        r"\[OW-GS\] hybrid SCF converged iterations=(\d+).*?density=\s*([0-9Ee+\-.]+).*?"
        r"band_energy_change=\s*([0-9Ee+\-.]+).*?eigensystem=\s*([0-9Ee+\-.]+).*?"
        r"electrons=\s*([0-9Ee+\-.]+).*?symmetry=\s*([0-9Ee+\-.]+)",
        log_text,
    )
    iteration_receipts = re.findall(
        r"\[HYBRID-SCF\] iteration=(\d+).*?density=\s*([0-9Ee+\-.]+).*?"
        r"band_energy_change=\s*([0-9Ee+\-.]+).*?eigensystem=\s*([0-9Ee+\-.]+).*?"
        r"electrons=\s*([0-9Ee+\-.]+).*?symmetry=\s*([0-9Ee+\-.]+).*?"
        r"band_energy=\s*([0-9Ee+\-.]+)",
        log_text,
    )
    evidence = {
        "return_code": return_code,
        "mpi_ranks": args.ranks,
        "omp_threads": 1,
        "elapsed_seconds": time.monotonic() - started,
        "peak_rss_kib_by_process": sorted(peak_rss_kib.values(), reverse=True),
        "peak_rss_kib_max": max(peak_rss_kib.values(), default=0),
        "convergence": convergence[-1] if convergence else None,
        "iteration_receipts": iteration_receipts,
        "occupied_checkpoint": (result_dir / "overlapping_wannier_occupied.chk").is_file(),
    }
    (result_dir / "hybrid_scf_evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    if return_code != 0:
        raise RuntimeError(f"Si64 hybrid SCF failed with exit {return_code}; see {log_path}")
    if not convergence or not evidence["occupied_checkpoint"]:
        raise RuntimeError("Si64 run ended without converged hybrid SCF checkpoint publication")
    if "[OW-GS] reused accepted route checkpoint" in log_text:
        raise RuntimeError("hybrid route incorrectly reused the legacy checkpoint")
    print(result_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
