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
            "fixed_step = 0.01\n"
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
         "--symmetry-mode", "constrained",
         "--symmetrize-eps", "2.5d-8", "--site-symmetry", "true",
         "--trial-step", "0.25", "--num-cg-steps", "2", "--num-iter", "75"],
        text=True, capture_output=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "2.5d-8" in (output / f"{seed}.wout").read_text()
    replay_text = (output / f"{seed}.wout").read_text().lower()
    assert "site_symmetry = true" in replay_text
    assert "trial_step = 0.25" in replay_text
    assert "fixed_step" not in replay_text
    assert "num_cg_steps = 2" in replay_text
    assert "num_iter = 75" in replay_text
    assert "1.0d-10" in (bundle / f"{seed}.win").read_text()
    receipt = json.loads((output / "replay-receipt.json").read_text())
    assert receipt["returncode"] == 0
    assert receipt["symmetry_mode"] == "constrained"
    assert receipt["site_symmetry"] == "true"
    assert receipt["trial_step"] == "0.25"
    assert receipt["num_cg_steps"] == 2
    assert receipt["num_iter"] == 75
    assert set(receipt["input_sha256"]) == {f"{seed}{s}" for s in (".win", ".dmn", ".mmn", ".amn", ".eig")}

    fixed_output = root / "fixed-output"
    fixed = subprocess.run(
        ["python3", str(DRIVER), str(bundle), "--seed", seed,
         "--executable", str(fake), "--output-directory", str(fixed_output),
         "--symmetry-mode", "constrained",
         "--fixed-step", "0.05", "--num-iter", "20"],
        text=True, capture_output=True,
    )
    assert fixed.returncode == 0, fixed.stdout + fixed.stderr
    fixed_text = (fixed_output / f"{seed}.wout").read_text().lower()
    assert "fixed_step = 0.05" in fixed_text
    assert "trial_step" not in fixed_text
    assert json.loads((fixed_output / "replay-receipt.json").read_text())["fixed_step"] == "0.05"

    with (bundle / f"{seed}.win").open("a") as stream:
        stream.write("site_symmetry = true\nsymmetrize_eps = 2.0d-10\n")

    unconstrained_output = root / "unconstrained-output"
    unconstrained = subprocess.run(
        ["python3", str(DRIVER), str(bundle), "--seed", seed,
         "--executable", str(fake), "--output-directory", str(unconstrained_output),
         "--symmetry-mode", "unconstrained",
         "--symmetrize-eps", "9.9d-9", "--site-symmetry", "true",
         "--trial-step", "0.125", "--num-cg-steps", "3", "--num-iter", "50"],
        text=True, capture_output=True,
    )
    assert unconstrained.returncode == 0, unconstrained.stdout + unconstrained.stderr
    unconstrained_text = (unconstrained_output / f"{seed}.wout").read_text().lower()
    assert "site_symmetry = false" in unconstrained_text
    assert unconstrained_text.count("site_symmetry") == 1
    assert "symmetrize_eps" not in unconstrained_text
    assert "trial_step = 0.125" in unconstrained_text
    assert "fixed_step" not in unconstrained_text
    unconstrained_receipt = json.loads(
        (unconstrained_output / "replay-receipt.json").read_text()
    )
    assert unconstrained_receipt["symmetry_mode"] == "unconstrained"
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
