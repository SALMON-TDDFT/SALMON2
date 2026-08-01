#!/usr/bin/env python3
"""Focused RED/GREEN test for polarization-derived HHG figure generation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np


repo = Path(__file__).resolve().parents[2]
tool = repo / "tools/plot_overlapping_wannier_hhg.py"

with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    ideal = root / "ideal"
    displaced = root / "displaced"
    order = np.linspace(0.0, 12.0, 241)
    for matrix, scale in ((ideal, 1.0), (displaced, 0.7)):
        for axis_index, axis in enumerate("xyz"):
            case = matrix / f"laser-hhg-{axis}"
            case.mkdir(parents=True)
            power = scale * (1.0e-12 + np.exp(-((order - 3.0) / 0.18) ** 2))
            power *= 1.0 + 0.02 * axis_index
            spectrum = np.column_stack((order * 0.05696, order * 1.55, order, power))
            np.savetxt(case / "hhg-spectrum.tsv", spectrum)
            (case / "hhg-summary.json").write_text(json.dumps({
                "spectrum_source": "polarization",
                "carrier_energy_ev": 1.55,
            }))
    prefix = root / "polarization-hhg"
    subprocess.run([
        sys.executable, str(tool), "--ideal-root", str(ideal),
        "--displaced-root", str(displaced), "--output-prefix", str(prefix),
        "--maximum-order", "12",
    ], check=True)
    for suffix in (".png", ".pdf", ".json"):
        assert prefix.with_suffix(suffix).is_file()
        assert prefix.with_suffix(suffix).stat().st_size > 100
    evidence = json.loads(prefix.with_suffix(".json").read_text())
    assert evidence["spectrum_source"] == "polarization"
    assert evidence["panels"] == ["ideal_xyz", "ideal_displaced_x"]
    assert evidence["maximum_order"] == 12.0
    assert evidence["input_sha256"]["ideal_x"] == hashlib.sha256(
        (ideal / "laser-hhg-x/hhg-spectrum.tsv").read_bytes()).hexdigest()

print("PASS polarization-derived HHG figure generation")

