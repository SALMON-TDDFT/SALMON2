#!/usr/bin/env python3
"""Run the fixed Si64 reference and overlapping-Wannier Cartesian matrix."""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import time
from pathlib import Path


DECOMPOSITIONS = {"2x2x2": "2, 2, 2"}
assert tuple(DECOMPOSITIONS) == ("2x2x2",), "Task 9 runs only the fixed 2x2x2 physical decomposition"
BOX_PROFILES = {
    "buffer5": (5, 5, 5),
    "buffer6": (6, 6, 6),
}
WINDOWS = {"c192-sp48": 192, "c256-sp48": 256}
TIME_RSS = re.compile(r"^\s*(\d+)\s+maximum resident set size\s*$", re.MULTILINE)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def materialize(template: Path, destination: Path, replacements: dict[str, str]) -> None:
    text = template.read_text()
    for key, value in replacements.items():
        text = text.replace(f"@{key}@", value)
    unresolved = re.findall(r"@[A-Z_]+@", text)
    if unresolved:
        raise RuntimeError(f"unresolved fixture placeholders: {unresolved}")
    destination.write_text(text)


def run_once(binary: Path, workdir: Path, log_name: str) -> tuple[float, float]:
    time_log = workdir / f"{log_name}.time"
    log = workdir / log_name
    command = [
        "/usr/bin/time",
        "-l",
        "-o",
        str(time_log),
        "mpirun",
        "-np",
        "8",
        str(binary),
    ]
    started = time.monotonic()
    with (workdir / "inputfile").open("rb") as input_stream, log.open("wb") as output:
        result = subprocess.run(
            command,
            cwd=workdir,
            stdin=input_stream,
            stdout=output,
            stderr=subprocess.STDOUT,
            env={**os.environ, "OMP_NUM_THREADS": "1"},
            check=False,
        )
    runtime = time.monotonic() - started
    if result.returncode != 0:
        raise RuntimeError(f"{workdir.name}/{log_name} failed with exit {result.returncode}")
    timing = time_log.read_text(errors="replace")
    match = TIME_RSS.search(timing)
    if not match:
        raise RuntimeError(f"maximum RSS missing from {time_log}")
    return runtime, int(match.group(1)) / (1024.0 * 1024.0)


def prepare_common(row: Path, atom: Path, pseudo: Path) -> None:
    row.mkdir(parents=True)
    shutil.copy2(atom, row / "atom.dat")
    shutil.copy2(pseudo, row / "C_rps.dat")


def append_external_evidence(log: Path, runtime: float, memory: float, manifest: Path) -> None:
    with log.open("a") as stream:
        stream.write(f"[OW-GS-EVIDENCE] runtime_s={runtime:.16e}\n")
        stream.write(f"[OW-GS-EVIDENCE] peak_memory_mib={memory:.16e}\n")
        stream.write(f"[OW-GS-EVIDENCE] checkpoint_manifest_sha256={sha256(manifest)}\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("result_root", type=Path)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--skip-reference", action="store_true")
    parser.add_argument("--minimum-row", action="store_true")
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    result_root = args.result_root.resolve()
    if result_root.exists():
        raise RuntimeError(f"result root must be fresh: {result_root}")
    result_root.mkdir(parents=True)
    fixture = args.repo / "tests/dg/data/si64_overlapping_wannier"
    atom = fixture / "atom.dat"
    pseudo = args.repo / "samples/exercise_01_C2H2_gs/C_rps.dat"
    if not args.skip_reference:
        reference = result_root / "normal-reference"
        prepare_common(reference, atom, pseudo)
        shutil.copy2(fixture / "inputfile_reference.in", reference / "inputfile")
        run_once(binary, reference, "run.log")
    for decomposition, decomposition_value in DECOMPOSITIONS.items():
        for box_profile, buffer_width in BOX_PROFILES.items():
            for window, candidate in WINDOWS.items():
                if args.minimum_row and (decomposition, box_profile, window) != (
                    "2x2x2",
                    "buffer6",
                    "c192-sp48",
                ):
                    continue
                row = result_root / (
                    f"decomp-{decomposition}_box-{box_profile}_window-{window}"
                )
                prepare_common(row, atom, pseudo)
                materialize(
                    fixture / "inputfile_ow.in",
                    row / "inputfile",
                    {
                        "DECOMPOSITION": decomposition_value,
                        "BUFFER": ", ".join(str(value) for value in buffer_width),
                        "CANDIDATE": str(candidate),
                        "TARGET": "0",
                    },
                )
                runtime, memory = run_once(binary, row, "run.log")
                manifest = row / "overlapping_wannier_gs.manifest"
                if not manifest.is_file():
                    raise RuntimeError(f"route checkpoint missing from {row}")
                manifest_digest_before = sha256(manifest)
                append_external_evidence(row / "run.log", runtime, memory, manifest)
                run_once(binary, row, "restart.log")
                if sha256(manifest) != manifest_digest_before:
                    raise RuntimeError(f"route checkpoint changed during restart reuse: {row}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
