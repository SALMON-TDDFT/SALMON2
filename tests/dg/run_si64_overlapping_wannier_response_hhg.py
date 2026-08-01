#!/usr/bin/env python3
"""Run the immutable genuine-Si64 polarization response/HHG matrix."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


AXES = {"x": "1d0,0d0,0d0", "y": "0d0,1d0,0d0", "z": "0d0,0d0,1d0"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def render(template: Path, destination: Path, replacements: dict[str, str]) -> None:
    text = template.read_text()
    for key, value in replacements.items():
        text = text.replace(f"@{key}@", value)
    unresolved = re.findall(r"@[A-Z_]+@", text)
    if unresolved:
        raise RuntimeError(f"unresolved input placeholders: {unresolved}")
    destination.write_text(text)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("checkpoint_dir", type=Path)
    parser.add_argument("result_root", type=Path)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    checkpoint = args.checkpoint_dir.resolve(strict=True)
    root = args.result_root.resolve()
    if root.exists():
        raise RuntimeError(f"result root must be fresh: {root}")
    atoms = (checkpoint / "atom.dat").read_text()
    atom_lines = [line for line in atoms.splitlines() if line.strip()]
    if len(atom_lines) != 64 or any("'Si'" not in line for line in atom_lines):
        raise RuntimeError("checkpoint provenance is not 64 silicon atoms")
    if atoms != (args.repo / "tests/dg/data/si64_overlapping_wannier_rt/atom.dat").read_text():
        raise RuntimeError("checkpoint Si64 coordinates differ from the tracked production fixture")
    if not (checkpoint / "Si_rps.dat").is_file():
        raise RuntimeError("genuine Si pseudopotential is missing")
    if sha256(checkpoint / "Si_rps.dat") != sha256(args.repo / "samples/exercise_04_bulkSi_gs/Si_rps.dat"):
        raise RuntimeError("checkpoint silicon pseudopotential differs from the tracked SALMON sample")
    manifest_path = checkpoint / "overlapping_wannier_gs.manifest"
    if b"SALMON_OW_GS_CHECKPOINT_V3" not in manifest_path.read_bytes():
        raise RuntimeError("accepted V3 checkpoint manifest is missing")
    root.mkdir(parents=True)
    fixtures = args.repo / "tests/dg/data/si64_overlapping_wannier_rt"
    cases: list[tuple[str, str, str, str, str]] = [
        ("fieldoff", "input_fieldoff.in", "x", "0.5d0", "800"),
        ("fieldoff-half-dt", "input_fieldoff.in", "x", "0.25d0", "1600"),
        ("impulse-x", "input_impulse.in", "x", "0.5d0", "800"),
        ("impulse-x-half", "input_impulse.in", "x", "0.5d0", "800"),
        ("impulse-y", "input_impulse.in", "y", "0.5d0", "800"),
        ("impulse-z", "input_impulse.in", "z", "0.5d0", "800"),
        ("laser-weak-x", "input_laser_weak.in", "x", "0.5d0", "800"),
        ("laser-hhg-x", "input_laser_hhg.in", "x", "0.5d0", "800"),
        ("laser-hhg-y", "input_laser_hhg.in", "y", "0.5d0", "800"),
        ("laser-hhg-z", "input_laser_hhg.in", "z", "0.5d0", "800"),
        ("laser-hhg-x-half-dt", "input_laser_hhg.in", "x", "0.25d0", "1600"),
    ]
    checkpoint_files = sorted(checkpoint.glob("overlapping_wannier_gs*"))
    checkpoint_hashes = {item.name: sha256(item) for item in checkpoint_files}
    checkpoint_digest = sha256(manifest_path)
    for name, template, axis, dt, nt in cases:
        case = root / name; case.mkdir()
        shutil.copy2(checkpoint / "atom.dat", case / "atom.dat")
        shutil.copy2(checkpoint / "Si_rps.dat", case / "Si_rps.dat")
        for item in checkpoint_files:
            shutil.copy2(item, case / item.name)
        impulse = "5d-3" if name == "impulse-x-half" else "1d-2"
        render(fixtures / template, case / "inputfile", {
            "AXIS": AXES[axis], "DT": dt, "NT": nt, "IMPULSE": impulse,
        })
        with (case / "inputfile").open("rb") as source, (case / "run.log").open("wb") as log:
            result = subprocess.run(["mpirun", "-np", "8", str(binary)], cwd=case,
                stdin=source, stdout=log, stderr=subprocess.STDOUT,
                env={**os.environ, "OMP_NUM_THREADS": "1"})
        if result.returncode:
            raise RuntimeError(f"{name} failed with exit {result.returncode}")
        observable = case / "overlapping_wannier_rt_observables.dat"
        restart = case / "overlapping_wannier_rt.restart"
        if not restart.is_file() or any(sha256(checkpoint / name) != expected
                                        for name, expected in checkpoint_hashes.items()):
            raise RuntimeError(f"{name}: restart missing or source V3 checkpoint changed")
        evidence = {
            "material": "Si", "atomic_number": 14, "atom_count": 64,
            "checkpoint_magic": "SALMON_OW_GS_CHECKPOINT_V3",
            "checkpoint_manifest_sha256": checkpoint_digest,
            "observable_sha256": sha256(observable), "axis": axis,
            "rt_restart_sha256": sha256(restart),
            "input_sha256": sha256(case / "inputfile"),
            "binary_sha256": sha256(binary),
            "dt_au": float(dt.removesuffix("d0")), "nt": int(nt),
        }
        (case / "manifest.json").write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
        if name.startswith("laser-hhg") or name == "laser-weak-x":
            background = root / ("fieldoff-half-dt" if "half-dt" in name else "fieldoff")
            subprocess.run([sys.executable, str(args.repo / "tools/analyze_overlapping_wannier_spectra.py"), "hhg",
                "--input", observable, "--background", background / "overlapping_wannier_rt_observables.dat",
                "--axis", axis, "--window", "exponential", "--damping-time", "1000",
                "--carrier-ev", "1.55", "--output", case / "hhg-spectrum.tsv",
                "--summary", case / "hhg-summary.json"], check=True)
            summary_path = case / "hhg-summary.json"
            summary = json.loads(summary_path.read_text())
            summary["fundamental_bin_error"] = min(
                abs(index * summary["frequency_resolution_au"] * 27.211386245988 / 1.55 - 1.0)
                for index in range(summary["sample_count"] // 2 + 1)
            )
            summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    if any(sha256(checkpoint / name) != expected
           for name, expected in checkpoint_hashes.items()):
        raise RuntimeError("source V3 checkpoint files changed during response/HHG matrix")
    print(root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
