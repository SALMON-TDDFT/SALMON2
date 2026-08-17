#!/usr/bin/env python3
from pathlib import Path
import json
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
DRIVER = ROOT / "tests/dg/replay_dg_wannier90_bundle.py"

with tempfile.TemporaryDirectory(prefix="w90-replay-test-") as name:
    root = Path(name)
    bundle = root / "bundle"
    output = root / "output"
    bundle.mkdir()
    seed = "sample"
    for suffix, payload in {
        ".win": "num_wann = 1\nsymmetrize_eps = 1.0d-10\n",
        ".dmn": "dmn\n",
        ".mmn": "mmn\n",
        ".amn": "amn\n",
        ".eig": "1 1 0.0\n",
    }.items():
        (bundle / f"{seed}{suffix}").write_text(payload)
    fake = root / "wannier90.x"
    fake.write_text("#!/bin/sh\ngrep symmetrize_eps \"$1.win\" > \"$1.wout\"\n")
    fake.chmod(0o755)
    result = subprocess.run(
        ["python3", str(DRIVER), str(bundle), "--seed", seed,
         "--executable", str(fake), "--output-directory", str(output),
         "--symmetrize-eps", "2.5d-8"],
        text=True, capture_output=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "2.5d-8" in (output / f"{seed}.wout").read_text()
    assert "1.0d-10" in (bundle / f"{seed}.win").read_text()
    receipt = json.loads((output / "replay-receipt.json").read_text())
    assert receipt["returncode"] == 0
    assert set(receipt["input_sha256"]) == {f"{seed}{s}" for s in (".win", ".dmn", ".mmn", ".amn", ".eig")}

    (bundle / f"{seed}.amn").unlink()
    missing = subprocess.run(
        ["python3", str(DRIVER), str(bundle), "--seed", seed,
         "--executable", str(fake), "--output-directory", str(output)],
        text=True, capture_output=True,
    )
    assert missing.returncode != 0
    assert "missing replay input" in missing.stderr.lower()

print("PASS standalone Wannier90 replay driver")
