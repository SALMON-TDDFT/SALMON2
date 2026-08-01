#!/usr/bin/env python3
"""Validate the immutable genuine-Si64 response/HHG production matrix."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


HARTREE_EV = 27.211386245988
FIELD_OFF_P_DRIFT_MAX_E_PER_BOHR2 = 1.0e-9
CURRENT_DERIVATIVE_RELATIVE_RMS_MAX = 5.0e-2
IMPULSE_LINEARITY_RELATIVE_MAX = 2.0e-2
CUBIC_EQUIVALENCE_RELATIVE_MAX = 5.0e-2
TRANSVERSE_POLARIZATION_RATIO_MAX = 5.0e-2
TIMESTEP_RELATIVE_MAX = 5.0e-2
HARMONIC_BIN_TOLERANCE = 1.0
FREQUENCY_RESOLUTION_MAX_EV = 0.45
NYQUIST_MIN_EV = 50.0
ODD_TO_EVEN_HARMONIC_POWER_MIN = 3.0
NONLINEAR_THIRD_HARMONIC_POWER_RATIO_MIN = 2.0e4
FORBIDDEN_MARKERS = ("dg_wpw", "time_evolution_dg_fragment", "initialization_rt")
REQUIRED_CASES = (
    "fieldoff", "fieldoff-half-dt", "impulse-x", "impulse-x-half", "impulse-y", "impulse-z",
    "laser-weak-x", "laser-hhg-x", "laser-hhg-y", "laser-hhg-z", "laser-hhg-x-half-dt",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_case(root: Path, name: str) -> tuple[np.ndarray, dict[str, object]]:
    case = root / name
    manifest = json.loads((case / "manifest.json").read_text())
    require(manifest["material"] == "Si", f"{name}: material is not Si")
    require(manifest["atomic_number"] == 14, f"{name}: atomic number is not 14")
    require(manifest["atom_count"] == 64, f"{name}: atom count is not 64")
    require(manifest["checkpoint_magic"] == "SALMON_OW_GS_CHECKPOINT_V3",
            f"{name}: checkpoint is not V3")
    observable = case / "overlapping_wannier_rt_observables.dat"
    require(manifest["observable_sha256"] == digest(observable),
            f"{name}: observable hash mismatch")
    require(manifest["rt_restart_sha256"] == digest(case / "overlapping_wannier_rt.restart"),
            f"{name}: RT restart hash mismatch")
    require(manifest["input_sha256"] == digest(case / "inputfile"),
            f"{name}: input hash mismatch")
    data = np.loadtxt(observable)
    require(data.ndim == 2 and data.shape[1] == 11 and np.isfinite(data).all(),
            f"{name}: invalid observable table")
    require(np.array_equal(data[:, 0], np.arange(len(data))),
            f"{name}: non-contiguous observable steps")
    log = (case / "run.log").read_text(errors="replace").lower()
    require(not any(marker in log for marker in FORBIDDEN_MARKERS),
            f"{name}: forbidden route marker in run log")
    return data, manifest


def relative_rms(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.sqrt(np.mean((left - right) ** 2)) / max(np.max(np.abs(right)), 1e-30))


def harmonic_power(spectrum: np.ndarray, order: int) -> float:
    index = int(np.argmin(np.abs(spectrum[:, 2] - order)))
    return float(spectrum[index, 3])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_root", type=Path)
    args = parser.parse_args()
    root = args.result_root.resolve(strict=True)
    cases = {name: load_case(root, name) for name in REQUIRED_CASES}
    require(len({item[1]["checkpoint_manifest_sha256"] for item in cases.values()}) == 1,
            "cases do not share one V3 checkpoint")
    require(len({item[1]["binary_sha256"] for item in cases.values()}) == 1,
            "cases do not share one executable")
    off = cases["fieldoff"][0]
    require(np.max(np.abs(off[:, 5:8] - off[0, 5:8])) <= FIELD_OFF_P_DRIFT_MAX_E_PER_BOHR2,
            "field-off polarization drift exceeds tolerance")

    impulse = cases["impulse-x"][0]
    delta = impulse[:, 5:8] - off[:, 5:8]
    derivative = (delta[2:, 0] - delta[:-2, 0]) / (impulse[2:, 1] - impulse[:-2, 1])
    current = impulse[1:-1, 8] - off[1:-1, 8]
    require(relative_rms(current, derivative) <= CURRENT_DERIVATIVE_RELATIVE_RMS_MAX,
            "current disagrees with centered polarization derivative")
    half = cases["impulse-x-half"][0][:, 5] - off[:, 5]
    require(relative_rms(delta[:, 0], 2.0 * half) <= IMPULSE_LINEARITY_RELATIVE_MAX,
            "impulse response is not amplitude-linear")
    amplitudes = []
    for axis, name in enumerate(("impulse-x", "impulse-y", "impulse-z")):
        response = cases[name][0][:, 5 + axis] - off[:, 5 + axis]
        amplitudes.append(np.max(np.abs(response)))
    require((max(amplitudes) - min(amplitudes)) / max(amplitudes) <= CUBIC_EQUIVALENCE_RELATIVE_MAX,
            "ideal Si64 Cartesian impulse amplitudes violate cubic equivalence")
    require(np.max(np.abs(delta[:, 1:])) / np.max(np.abs(delta[:, 0])) <=
            TRANSVERSE_POLARIZATION_RATIO_MAX, "impulse transverse response exceeds tolerance")

    weak = cases["laser-weak-x"][0][:, 5] - off[:, 5]
    strong = cases["laser-hhg-x"][0][:, 5] - off[:, 5]
    require(np.max(np.abs(strong)) > np.max(np.abs(weak)),
            "strong-field polarization does not exceed weak-field polarization")
    strong_full = cases["laser-hhg-x"][0][:, 5:8] - off[:, 5:8]
    require(np.max(np.abs(strong_full[:, 1:])) / np.max(np.abs(strong_full[:, 0])) <=
            TRANSVERSE_POLARIZATION_RATIO_MAX, "laser transverse response exceeds tolerance")
    summaries = [json.loads((root / f"laser-hhg-{axis}" / "hhg-summary.json").read_text())
                 for axis in "xyz"]
    require(all(item["spectrum_source"] == "polarization" for item in summaries),
            "HHG spectrum source is not polarization")
    require(all(item["fundamental_bin_error"] <= HARMONIC_BIN_TOLERANCE for item in summaries),
            "laser fundamental is outside the accepted spectral bin")
    require(all(item["frequency_resolution_au"] * HARTREE_EV <= FREQUENCY_RESOLUTION_MAX_EV
                for item in summaries), "frequency resolution is insufficient")
    require(all(item["nyquist_au"] * HARTREE_EV >= NYQUIST_MIN_EV for item in summaries),
            "Nyquist energy is insufficient")
    spectra = [np.loadtxt(root / f"laser-hhg-{axis}" / "hhg-spectrum.tsv") for axis in "xyz"]
    odd_powers = [harmonic_power(spectrum, 3) + harmonic_power(spectrum, 5) for spectrum in spectra]
    even_powers = [harmonic_power(spectrum, 2) + harmonic_power(spectrum, 4) for spectrum in spectra]
    require(all(odd / max(even, 1e-300) >= ODD_TO_EVEN_HARMONIC_POWER_MIN
                for odd, even in zip(odd_powers, even_powers)),
            "ideal centrosymmetric Si64 lacks odd-harmonic dominance")
    require((max(odd_powers) - min(odd_powers)) / max(odd_powers) <=
            CUBIC_EQUIVALENCE_RELATIVE_MAX, "HHG axes violate cubic equivalence")
    weak_spectrum = np.loadtxt(root / "laser-weak-x" / "hhg-spectrum.tsv")
    require(harmonic_power(spectra[0], 3) /
            max(harmonic_power(weak_spectrum, 3), 1e-300) >=
            NONLINEAR_THIRD_HARMONIC_POWER_RATIO_MIN,
            "third-harmonic scaling does not establish nonlinear response")
    base = np.loadtxt(root / "laser-hhg-x" / "hhg-spectrum.tsv")
    half_dt = np.loadtxt(root / "laser-hhg-x-half-dt" / "hhg-spectrum.tsv")
    interpolated = np.interp(base[:, 0], half_dt[:, 0], half_dt[:, 3])
    band = base[:, 2] <= 15.0
    require(relative_rms(base[band, 3], interpolated[band]) <= TIMESTEP_RELATIVE_MAX,
            "HHG spectrum fails time-step convergence")
    print("PASS genuine-Si64 overlapping-Wannier response/HHG matrix")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
