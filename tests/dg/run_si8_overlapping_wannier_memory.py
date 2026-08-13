#!/usr/bin/env python3
"""Run Si8 on eight MPI ranks while recording nonintrusive memory evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import time


PAGE_SIZE = re.compile(r"page size of\s+(\d+)\s+bytes", re.IGNORECASE)
PAGE_LINE = re.compile(r"^Pages\s+(free|inactive|speculative|purgeable):\s+(\d+)\.", re.MULTILINE)


def parse_available_bytes(text: str) -> int:
    match = PAGE_SIZE.search(text)
    if not match:
        raise ValueError("vm_stat page size is missing")
    pages = {name: int(value) for name, value in PAGE_LINE.findall(text)}
    required = {"free", "inactive", "speculative", "purgeable"}
    if not required.issubset(pages):
        raise ValueError("vm_stat available-page fields are missing")
    return int(match.group(1)) * sum(pages[name] for name in required)


def parse_rank_rss_kib(text: str, root_pid: int) -> dict[int, int]:
    rows: dict[int, tuple[int, int, str]] = {}
    for line in text.splitlines():
        fields = line.strip().split(None, 3)
        if len(fields) != 4:
            continue
        try:
            pid, ppid, rss = map(int, fields[:3])
        except ValueError:
            continue
        rows[pid] = (ppid, rss, fields[3])
    descendants = {root_pid}
    changed = True
    while changed:
        changed = False
        for pid, (ppid, _, _) in rows.items():
            if ppid in descendants and pid not in descendants:
                descendants.add(pid)
                changed = True
    return {
        pid: rss
        for pid, (_, rss, command) in rows.items()
        if pid in descendants and Path(command.split()[0]).name == "salmon"
    }


def latest_phase(text: str) -> str:
    phases = [line.strip() for line in text.splitlines() if "[OW-GS-DIAGNOSTIC]" in line]
    return phases[-1] if phases else "startup"


def safety_reason(
    available_bytes: int,
    rank_rss_kib: dict[int, int],
    available_floor_bytes: int,
    rank_ceiling_kib: int,
) -> str | None:
    if available_bytes < available_floor_bytes:
        return "available_memory"
    if rank_rss_kib and max(rank_rss_kib.values()) > rank_ceiling_kib:
        return "rank_rss"
    return None


def build_launch_command(binary: Path, which=shutil.which) -> list[str]:
    command = ["/usr/bin/nice", "-n", "15", "mpirun", "-np", "8", str(binary)]
    taskpolicy = which("taskpolicy")
    return [taskpolicy, "-b", *command] if taskpolicy else command


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tail(path: Path, maximum: int = 1024 * 1024) -> str:
    if not path.exists():
        return ""
    with path.open("rb") as stream:
        stream.seek(max(0, path.stat().st_size - maximum))
        return stream.read().decode(errors="replace")


def terminate_group(process: subprocess.Popen[bytes]) -> None:
    try:
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=30)
    except ProcessLookupError:
        return
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait(timeout=10)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("result_root", type=Path)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--interval", type=float, default=30.0)
    parser.add_argument("--available-floor-gib", type=float, default=8.0)
    parser.add_argument("--rank-ceiling-gib", type=float, default=3.0)
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    result = args.result_root.resolve()
    if result.exists():
        raise RuntimeError(f"result directory must be fresh: {result}")
    result.mkdir(parents=True)
    fixture = args.repo / "tests/dg/data/si8_overlapping_wannier"
    pseudo = args.repo / "samples/benchmark/bulk_Si/Si_rps.dat"
    shutil.copy2(fixture / "inputfile.in", result / "inputfile")
    shutil.copy2(fixture / "atom.dat", result / "atom.dat")
    shutil.copy2(pseudo, result / "Si_rps.dat")
    provenance = {
        path.name: sha256(path)
        for path in (binary, result / "inputfile", result / "atom.dat", result / "Si_rps.dat")
    }
    (result / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    log_path = result / "run.log"
    monitor_path = result / "memory.csv"
    env = {
        **os.environ,
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "GFORTRAN_UNBUFFERED_ALL": "1",
    }
    command = build_launch_command(binary)
    available_floor = int(args.available_floor_gib * (1 << 30))
    rank_ceiling = int(args.rank_ceiling_gib * (1 << 20))
    stopped_reason = ""
    with (result / "inputfile").open("rb") as input_stream, log_path.open("wb") as log_stream, \
            monitor_path.open("w", newline="") as monitor_stream:
        writer = csv.writer(monitor_stream)
        writer.writerow(["elapsed_s", "available_bytes", "rank_count", "total_rss_kib", "max_rss_kib", "phase"])
        monitor_stream.flush()
        started = time.monotonic()
        process = subprocess.Popen(
            command, cwd=result, stdin=input_stream, stdout=log_stream,
            stderr=subprocess.STDOUT, env=env, start_new_session=True,
        )
        try:
            while process.poll() is None:
                vm_text = subprocess.run(["vm_stat"], text=True, capture_output=True, check=True).stdout
                ps_text = subprocess.run(
                    ["ps", "-axo", "pid=,ppid=,rss=,command="],
                    text=True, capture_output=True, check=True,
                ).stdout
                available = parse_available_bytes(vm_text)
                rss = parse_rank_rss_kib(ps_text, process.pid)
                phase = latest_phase(read_tail(log_path))
                writer.writerow([
                    f"{time.monotonic() - started:.3f}", available, len(rss),
                    sum(rss.values()), max(rss.values(), default=0), phase,
                ])
                monitor_stream.flush()
                stopped_reason = safety_reason(available, rss, available_floor, rank_ceiling) or ""
                if stopped_reason:
                    (result / "SAFETY_STOP").write_text(stopped_reason + "\n")
                    terminate_group(process)
                    break
                time.sleep(max(1.0, args.interval))
        except BaseException:
            terminate_group(process)
            raise
        return_code = process.wait()
    if stopped_reason:
        raise RuntimeError(f"Si8 memory safety gate stopped the run: {stopped_reason}")
    if return_code != 0:
        raise RuntimeError(f"Si8 MPI8 run failed with exit {return_code}")
    print(f"Si8 MPI8 memory diagnostic completed: {result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
