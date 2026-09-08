#!/usr/bin/env python3
"""Replay exported Gamma-point Wannier90 inputs without rerunning SALMON."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time


COMMON_SUFFIXES = (".win", ".mmn", ".amn", ".eig")
CONSTRAINED_SUFFIXES = (".win", ".dmn", ".mmn", ".amn", ".eig")


def replace_keyword(text: str, keyword: str, value: str) -> str:
    lines = text.splitlines()
    replacement = f"{keyword} = {value}"
    for index, line in enumerate(lines):
        if line.strip().lower().startswith(keyword.lower()):
            lines[index] = replacement
            break
    else:
        lines.append(replacement)
    return "\n".join(lines) + "\n"


def remove_keyword(text: str, keyword: str) -> str:
    lines = [line for line in text.splitlines()
             if not line.strip().lower().startswith(keyword.lower())]
    return "\n".join(lines) + "\n"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("--seed", default="overlapping_wannier_mlwf")
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument(
        "--symmetry-mode", choices=("constrained", "unconstrained"), required=True
    )
    parser.add_argument("--symmetrize-eps")
    parser.add_argument("--site-symmetry", choices=("true", "false"))
    step = parser.add_mutually_exclusive_group()
    step.add_argument("--trial-step")
    step.add_argument("--fixed-step")
    parser.add_argument("--num-cg-steps", type=int)
    parser.add_argument("--num-iter", type=int)
    args = parser.parse_args()

    suffixes = (
        CONSTRAINED_SUFFIXES
        if args.symmetry_mode == "constrained"
        else COMMON_SUFFIXES
    )
    inputs = [args.bundle / f"{args.seed}{suffix}" for suffix in suffixes]
    missing = [path for path in inputs if not path.is_file()]
    if missing:
        print("missing replay input: " + ", ".join(map(str, missing)), file=sys.stderr)
        return 2
    if not args.executable.is_file():
        print(f"Wannier90 executable is missing: {args.executable}", file=sys.stderr)
        return 2

    args.output_directory.mkdir(parents=True, exist_ok=True)
    input_hashes = {path.name: sha256(path) for path in inputs}
    with tempfile.TemporaryDirectory(prefix="salmon-w90-replay-") as name:
        run_directory = Path(name)
        for path in inputs:
            shutil.copy2(path, run_directory / path.name)
        win = run_directory / f"{args.seed}.win"
        win_text = win.read_text()
        effective_symmetrize_eps = args.symmetrize_eps
        effective_site_symmetry = args.site_symmetry
        if args.symmetry_mode == "unconstrained":
            win_text = remove_keyword(win_text, "symmetrize_eps")
            win_text = remove_keyword(win_text, "site_symmetry")
            win_text = replace_keyword(win_text, "site_symmetry", "false")
            effective_symmetrize_eps = None
            effective_site_symmetry = "false"
        else:
            if args.symmetrize_eps:
                win_text = replace_keyword(win_text, "symmetrize_eps", args.symmetrize_eps)
            if args.site_symmetry:
                win_text = replace_keyword(win_text, "site_symmetry", args.site_symmetry)
        if args.trial_step:
            win_text = remove_keyword(win_text, "fixed_step")
            win_text = replace_keyword(win_text, "trial_step", args.trial_step)
        if args.fixed_step:
            win_text = remove_keyword(win_text, "trial_step")
            win_text = replace_keyword(win_text, "fixed_step", args.fixed_step)
        if args.num_cg_steps is not None:
            win_text = replace_keyword(win_text, "num_cg_steps", str(args.num_cg_steps))
        if args.num_iter is not None:
            win_text = replace_keyword(win_text, "num_iter", str(args.num_iter))
        win.write_text(win_text)
        started = time.monotonic()
        result = subprocess.run(
            [str(args.executable), args.seed], cwd=run_directory,
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        elapsed = time.monotonic() - started
        wout = run_directory / f"{args.seed}.wout"
        if wout.is_file():
            shutil.copy2(wout, args.output_directory / wout.name)
        (args.output_directory / "wannier90-replay.log").write_text(result.stdout)

    receipt = {
        "command": [str(args.executable), args.seed],
        "elapsed_seconds": elapsed,
        "input_sha256": input_hashes,
        "returncode": result.returncode,
        "seed": args.seed,
        "symmetry_mode": args.symmetry_mode,
        "symmetrize_eps": effective_symmetrize_eps,
        "site_symmetry": effective_site_symmetry,
        "trial_step": args.trial_step,
        "fixed_step": args.fixed_step,
        "num_cg_steps": args.num_cg_steps,
        "num_iter": args.num_iter,
    }
    (args.output_directory / "replay-receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
