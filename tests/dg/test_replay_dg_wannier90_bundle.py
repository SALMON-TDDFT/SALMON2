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
        ".win": (
            "num_wann = 1\n"
            "site_symmetry = true\n"
            "symmetrize_eps = 1.0d-10\n"
        ),
        ".dmn": "dmn\n",
        ".mmn": "mmn\n",
        ".amn": "amn\n",
        ".eig": "1 1 0.0\n",
    }.items():
        (bundle / f"{seed}{suffix}").write_text(payload)
    fake = root / "wannier90.x"
    fake.write_text(
        "#!/bin/sh\n"
        "if grep -qi 'site_symmetry *= *true' \"$1.win\" && [ ! -e \"$1.dmn\" ]; then exit 92; fi\n"
        "if grep -qi 'site_symmetry *= *false' \"$1.win\" && [ -e \"$1.dmn\" ]; then exit 91; fi\n"
        "cat \"$1.win\" > \"$1.wout\"\n"
    )
    fake.chmod(0o755)
    result = subprocess.run(
        ["python3", str(DRIVER), str(bundle), "--seed", seed,
         "--executable", str(fake), "--output-directory", str(output),
         "--symmetry-mode", "constrained", "--symmetrize-eps", "2.5d-8"],
        text=True, capture_output=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    replay_text = (output / f"{seed}.wout").read_text().lower()
    assert "2.5d-8" in replay_text
    assert "site_symmetry = true" in replay_text
    assert "1.0d-10" in (bundle / f"{seed}.win").read_text()
    receipt = json.loads((output / "replay-receipt.json").read_text())
    assert receipt["returncode"] == 0
    assert receipt["symmetry_mode"] == "constrained"
    assert set(receipt["input_sha256"]) == {
        f"{seed}{suffix}"
        for suffix in (".win", ".dmn", ".mmn", ".amn", ".eig")
    }

    with (bundle / f"{seed}.win").open("a") as stream:
        stream.write("site_symmetry = true\nsymmetrize_eps = 2.0d-10\n")

    unconstrained_output = root / "unconstrained-output"
    unconstrained = subprocess.run(
        ["python3", str(DRIVER), str(bundle), "--seed", seed,
         "--executable", str(fake), "--output-directory", str(unconstrained_output),
         "--symmetry-mode", "unconstrained", "--symmetrize-eps", "9.9d-9"],
        text=True, capture_output=True,
    )
    assert unconstrained.returncode == 0, unconstrained.stdout + unconstrained.stderr
    unconstrained_text = (unconstrained_output / f"{seed}.wout").read_text().lower()
    assert "site_symmetry = false" in unconstrained_text
    assert unconstrained_text.count("site_symmetry") == 1
    assert "symmetrize_eps" not in unconstrained_text
    unconstrained_receipt = json.loads(
        (unconstrained_output / "replay-receipt.json").read_text()
    )
    assert unconstrained_receipt["symmetry_mode"] == "unconstrained"
    assert unconstrained_receipt["symmetrize_eps"] is None
    assert set(unconstrained_receipt["input_sha256"]) == {
        f"{seed}{suffix}" for suffix in (".win", ".mmn", ".amn", ".eig")
    }
    assert "site_symmetry = true" in (bundle / f"{seed}.win").read_text().lower()
    assert "symmetrize_eps = 1.0d-10" in (bundle / f"{seed}.win").read_text().lower()

    (bundle / f"{seed}.dmn").unlink()

    missing_mode = subprocess.run(
        ["python3", str(DRIVER), str(bundle), "--seed", seed,
         "--executable", str(fake), "--output-directory", str(root / "missing-mode")],
        text=True, capture_output=True,
    )
    assert missing_mode.returncode != 0
    assert "--symmetry-mode" in missing_mode.stderr

    missing_dmn = subprocess.run(
        ["python3", str(DRIVER), str(bundle), "--seed", seed,
         "--executable", str(fake), "--output-directory", str(root / "missing-dmn"),
         "--symmetry-mode", "constrained"],
        text=True, capture_output=True,
    )
    assert missing_dmn.returncode != 0
    assert "missing replay input" in missing_dmn.stderr.lower()

    (bundle / f"{seed}.amn").unlink()
    missing = subprocess.run(
        ["python3", str(DRIVER), str(bundle), "--seed", seed,
         "--executable", str(fake), "--output-directory", str(output),
         "--symmetry-mode", "unconstrained"],
        text=True, capture_output=True,
    )
    assert missing.returncode != 0
    assert "missing replay input" in missing.stderr.lower()

print("PASS standalone Wannier90 replay driver")
