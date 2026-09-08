#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
ANALYZER = ROOT / "tools/analyze_overlapping_wannier_spectra.py"
HEADER = """# SALMON_OW_COEFFICIENT_RT_OBSERVABLES_V1
# sign electronic_charge=-1
# units time=au electric_field=au polarization=electron/bohr^2 current=electron/(bohr^2*au)
# volume_au   1.00000000000000000E+000
# provenance 2 1 1 11 12 13 14 cell_wrapped_length_velocity
# step time Ex Ey Ez Px Py Pz Jx Jy Jz
"""


def write_signal(path: Path, time: np.ndarray, polarization: np.ndarray) -> None:
    omega = 2.0 * np.pi * 17 / (len(time) * (time[1] - time[0]))
    current = omega * np.cos(omega * time)
    rows = np.column_stack(
        [np.arange(len(time)), time, np.zeros((len(time), 3)), polarization,
         current, np.zeros(len(time)), np.zeros(len(time))]
    )
    path.write_text(HEADER + "\n".join(" ".join(f"{value:.17e}" for value in row) for row in rows) + "\n")


def run(*args: object, expect_ok: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        [sys.executable, str(ANALYZER), *(str(arg) for arg in args)],
        text=True, capture_output=True,
    )
    assert (result.returncode == 0) == expect_ok, (result.stdout, result.stderr)
    return result


with tempfile.TemporaryDirectory(prefix="ow-spectra-") as temporary:
    work = Path(temporary)
    count, dt, bin_index = 1024, 0.2, 17
    time = np.arange(count) * dt
    omega = 2.0 * np.pi * bin_index / (count * dt)
    background = np.column_stack([0.02 + 1e-5 * time, np.zeros(count), np.zeros(count)])
    driven = background.copy()
    driven[:, 0] += 3e-4 * np.sin(omega * time)
    write_signal(work / "background.dat", time, background)
    write_signal(work / "driven.dat", time, driven)

    linear_tsv, linear_json = work / "linear.tsv", work / "linear.json"
    run("linear", "--input", work / "driven.dat", "--background", work / "background.dat",
        "--axis", "x", "--window", "hann", "--impulse-area", "2e-4",
        "--output", linear_tsv, "--summary", linear_json)
    linear = np.loadtxt(linear_tsv, comments="#")
    assert abs(linear[np.argmax(linear[:, 4]), 0] - omega) <= 2.0 * np.pi / (count * dt)
    summary = json.loads(linear_json.read_text())
    assert summary["spectrum_source"] == "polarization"
    assert summary["field_off_subtracted"] is True
    assert summary["impulse_area_au"] == 2e-4
    assert summary["frequency_resolution_au"] == 2.0 * np.pi / (count * dt)
    linear_tsv_bytes, linear_json_bytes = linear_tsv.read_bytes(), linear_json.read_bytes()
    run("linear", "--input", work / "driven.dat", "--background", work / "background.dat",
        "--axis", "x", "--window", "hann", "--impulse-area", "2e-4",
        "--output", linear_tsv, "--summary", linear_json)
    assert linear_tsv.read_bytes() == linear_tsv_bytes
    assert linear_json.read_bytes() == linear_json_bytes

    hhg_tsv, hhg_json = work / "hhg.tsv", work / "hhg.json"
    carrier_ev = omega * 27.211386245988
    run("hhg", "--input", work / "driven.dat", "--background", work / "background.dat",
        "--axis", "x", "--window", "exponential", "--damping-time", "1000",
        "--carrier-ev", carrier_ev, "--output", hhg_tsv, "--summary", hhg_json)
    hhg = np.loadtxt(hhg_tsv, comments="#")
    peak = hhg[np.argmax(hhg[1:, 3]) + 1]
    assert abs(peak[2] - 1.0) <= 1.0 / bin_index
    assert np.allclose(hhg[:, 3], hhg[:, 0] ** 4 * hhg[:, 4])
    assert json.loads(hhg_json.read_text())["formula"] == "omega^4*abs(P_omega)^2"

    bad = work / "nonuniform.dat"
    bad_time = time.copy(); bad_time[20] += 0.01
    write_signal(bad, bad_time, driven)
    run("hhg", "--input", bad, "--background", work / "background.dat", "--axis", "x",
        "--window", "hann", "--carrier-ev", carrier_ev, "--output", work / "bad.tsv",
        "--summary", work / "bad.json", expect_ok=False)
    truncated = work / "truncated.dat"; truncated.write_text(HEADER + "0 0 0\n")
    run("linear", "--input", truncated, "--background", work / "background.dat", "--axis", "x",
        "--window", "hann", "--impulse-area", "1", "--output", work / "truncated.tsv",
        "--summary", work / "truncated.json", expect_ok=False)
    nonfinite = work / "nonfinite.dat"
    nonfinite_lines = (work / "driven.dat").read_text().splitlines()
    first_data = next(index for index, line in enumerate(nonfinite_lines) if not line.startswith("#"))
    tokens = nonfinite_lines[first_data].split(); tokens[5] = "nan"
    nonfinite_lines[first_data] = " ".join(tokens)
    nonfinite.write_text("\n".join(nonfinite_lines) + "\n")
    run("linear", "--input", nonfinite, "--background", work / "background.dat", "--axis", "x",
        "--window", "hann", "--impulse-area", "1", "--output", work / "nan.tsv",
        "--summary", work / "nan.json", expect_ok=False)

print("PASS overlapping-Wannier polarization spectrum analyzer")
