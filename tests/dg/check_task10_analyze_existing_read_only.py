#!/usr/bin/env python3
"""Prove both Task10 analyze-existing CLIs leave their raw result trees unchanged."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import stat
import subprocess
import sys
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def tree_receipt(root: Path) -> list[dict[str, int | str]]:
    receipt = []
    for path in sorted(root.rglob("*"), key=lambda item: item.relative_to(root).as_posix()):
        info = path.lstat()
        item: dict[str, int | str] = {
            "path": path.relative_to(root).as_posix(),
            "mode": f"{stat.S_IFMT(info.st_mode):o}:{stat.S_IMODE(info.st_mode):04o}",
            "size": info.st_size,
            "mtime_ns": info.st_mtime_ns,
        }
        if stat.S_ISREG(info.st_mode):
            item["sha256"] = sha256(path)
        elif stat.S_ISLNK(info.st_mode):
            item["target"] = os.readlink(path)
        receipt.append(item)
    return receipt


def run_read_only(runner: Path, result_dir: Path) -> dict:
    before = tree_receipt(result_dir)
    completed = subprocess.run(
        [sys.executable, str(runner), "--result-dir", str(result_dir), "--analyze-existing"],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    after = tree_receipt(result_dir)
    if after != before:
        raise RuntimeError(f"analyze-existing modified its raw result tree: {result_dir}")
    payload = json.loads(completed.stdout)
    if payload.get("status") != "PASS" or payload.get("result_dir") != str(result_dir):
        raise RuntimeError(f"analyze-existing did not emit a complete JSON validation: {runner}")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--si64-result-dir", type=Path, required=True)
    parser.add_argument("--smoke-result-dir", type=Path, required=True)
    args = parser.parse_args()
    repository = args.repository.resolve(strict=True)
    si64_result = args.si64_result_dir.resolve(strict=True)
    smoke_result = args.smoke_result_dir.resolve(strict=True)
    results = {
        "si64": run_read_only(
            repository / "tests/dg/run_dg_hybrid_si64_scdm_reuse.py", si64_result),
        "smoke": run_read_only(
            repository / "tests/dg/run_dg_fragment_wf_production_smoke.py", smoke_result),
    }
    print(json.dumps({"status": "PASS", "results": results}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
