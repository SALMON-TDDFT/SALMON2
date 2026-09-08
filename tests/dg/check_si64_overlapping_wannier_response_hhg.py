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
IMPULSE_LINEARITY_RELATIVE_MAX = 2.0e-2
CUBIC_EQUIVALENCE_RELATIVE_MAX = 5.0e-2
TRANSVERSE_POLARIZATION_RATIO_MAX = 5.0e-2
TIMESTEP_RELATIVE_MAX = 5.0e-2
EVEN_TO_THIRD_POWER_MAX = 1.0e-3
HARMONIC_BIN_TOLERANCE = 1.0
FREQUENCY_RESOLUTION_MAX_EV = 0.45
NYQUIST_MIN_EV = 50.0
EXPECTED_LASER_CYCLES = 10.0
EXPECTED_POST_PULSE_CYCLES = 2.0
EXPECTED_LASER_DT_AU = 2.0
EXPECTED_LASER_REFERENCE_DT_AU = 1.0
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
    require(manifest.get("structure") == "ideal", f"{name}: structure is not ideal")
    require(int(manifest.get("global_point_group_order", 0)) > 1,
            f"{name}: missing full-system point-group evidence")
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


def harmonic_local_morphology(spectrum: np.ndarray, order: int) -> dict[str, float | str]:
    """Classify the nearest harmonic bin against its immediate spectral neighbors."""
    index = int(np.argmin(np.abs(spectrum[:, 2] - order)))
    require(0 < index < len(spectrum) - 1, f"harmonic {order} lacks neighboring bins")
    center = float(spectrum[index, 3])
    left = float(spectrum[index - 1, 3])
    right = float(spectrum[index + 1, 3])
    kind = "peak" if center > max(left, right) else "dip" if center < min(left, right) else "slope"
    return {
        "kind": kind,
        "neighbor_ratio": center / max(left, right, 1e-300),
        "harmonic_order": float(spectrum[index, 2]),
        "power": center,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_root", type=Path)
    args = parser.parse_args()
    root = args.result_root.resolve(strict=True)
    cases = {name: load_case(root, name) for name in REQUIRED_CASES}
    require(all(item[1]["structure"] == "ideal" for item in cases.values()),
            "ideal matrix is not marked ideal")
    require(len({item[1]["checkpoint_manifest_sha256"] for item in cases.values()}) == 1,
            "cases do not share one V3 checkpoint")
    require(len({item[1]["binary_sha256"] for item in cases.values()}) == 1,
            "cases do not share one executable")
    require(all(item[1]["global_point_group_order"] ==
                cases["fieldoff"][1]["global_point_group_order"] for item in cases.values()),
            "ideal Si64 cases do not share one full-system point group")
    require(all(item[1].get("global_inversion_promoted") is True for item in cases.values()),
            "ideal Si64 V3 route did not promote exact fragment-spanning inversion")
    laser_manifests = [cases[f"laser-hhg-{axis}"][1] for axis in "xyz"]
    require(all(item.get("laser_cycles") == EXPECTED_LASER_CYCLES and
                item.get("post_pulse_cycles") == EXPECTED_POST_PULSE_CYCLES and
                item["dt_au"] == EXPECTED_LASER_DT_AU for item in laser_manifests),
            "ideal HHG matrix does not use the accepted long production pulse")
    reference_manifest = cases["laser-hhg-x-half-dt"][1]
    require(reference_manifest.get("laser_cycles") == EXPECTED_LASER_CYCLES and
            reference_manifest.get("post_pulse_cycles") == EXPECTED_POST_PULSE_CYCLES and
            reference_manifest["dt_au"] == EXPECTED_LASER_REFERENCE_DT_AU and
            reference_manifest["total_duration_au"] == laser_manifests[0]["total_duration_au"],
            "ideal HHG reference does not share the accepted physical duration")
    off = cases["fieldoff"][0]
    require(np.max(np.abs(off[:, 5:8] - off[0, 5:8])) <= FIELD_OFF_P_DRIFT_MAX_E_PER_BOHR2,
            "field-off polarization drift exceeds tolerance")

    impulse = cases["impulse-x"][0]
    delta = impulse[:, 5:8] - off[:, 5:8]
    derivative = (delta[2:, 0] - delta[:-2, 0]) / (impulse[2:, 1] - impulse[:-2, 1])
    current = impulse[1:-1, 8] - off[1:-1, 8]
    current_derivative_mismatch = relative_rms(current, derivative)
    half = cases["impulse-x-half"][0][:, 5] - off[:, 5]
    require(relative_rms(delta[:, 0], 2.0 * half) <= IMPULSE_LINEARITY_RELATIVE_MAX,
            "impulse response is not amplitude-linear")
    amplitudes = []
    for axis, name in enumerate(("impulse-x", "impulse-y", "impulse-z")):
        response = cases[name][0][:, 5 + axis] - off[:, 5 + axis]
        amplitudes.append(np.max(np.abs(response)))
    require((max(amplitudes) - min(amplitudes)) / max(amplitudes) <= CUBIC_EQUIVALENCE_RELATIVE_MAX,
            "ideal Si64 Cartesian impulse amplitudes violate cubic equivalence")
    impulse_transverse_ratio = float(
        np.max(np.abs(delta[:, 1:])) / np.max(np.abs(delta[:, 0])))
    require(impulse_transverse_ratio <= TRANSVERSE_POLARIZATION_RATIO_MAX,
            "impulse transverse response exceeds tolerance")

    weak = cases["laser-weak-x"][0][:, 5] - off[:, 5]
    strong = cases["laser-hhg-x"][0][:, 5] - off[:, 5]
    require(np.max(np.abs(strong)) > np.max(np.abs(weak)),
            "strong-field polarization does not exceed weak-field polarization")
    strong_full = cases["laser-hhg-x"][0][:, 5:8] - off[:, 5:8]
    laser_transverse_ratio = float(
        np.max(np.abs(strong_full[:, 1:])) / np.max(np.abs(strong_full[:, 0])))
    require(laser_transverse_ratio <= TRANSVERSE_POLARIZATION_RATIO_MAX,
            "laser transverse response exceeds tolerance")
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
    morphologies = [{order: harmonic_local_morphology(spectrum, order) for order in (2, 3, 4)}
                    for spectrum in spectra]
    require(all(item[3]["kind"] == "peak" for item in morphologies),
            "ideal Si64 third harmonic is not a local spectral peak")
    require(all(harmonic_power(spectrum, 2) / max(harmonic_power(spectrum, 3), 1e-300) <=
                EVEN_TO_THIRD_POWER_MAX and
                harmonic_power(spectrum, 4) / max(harmonic_power(spectrum, 3), 1e-300) <=
                EVEN_TO_THIRD_POWER_MAX for spectrum in spectra),
            "ideal centrosymmetric Si64 does not suppress H2/H4 relative to H3")
    odd_powers = [harmonic_power(spectrum, 3) + harmonic_power(spectrum, 5) for spectrum in spectra]
    require((max(odd_powers) - min(odd_powers)) / max(odd_powers) <=
            CUBIC_EQUIVALENCE_RELATIVE_MAX, "HHG axes violate cubic equivalence")
    weak_spectrum = np.loadtxt(root / "laser-weak-x" / "hhg-spectrum.tsv")
    nonlinear_third_harmonic_power_ratio = harmonic_power(spectra[0], 3) / max(
        harmonic_power(weak_spectrum, 3), 1e-300)
    base = np.loadtxt(root / "laser-hhg-x" / "hhg-spectrum.tsv")
    half_dt = np.loadtxt(root / "laser-hhg-x-half-dt" / "hhg-spectrum.tsv")
    interpolated = np.interp(base[:, 0], half_dt[:, 0], half_dt[:, 3])
    band = base[:, 2] <= 15.0
    require(relative_rms(base[band, 3], interpolated[band]) <= TIMESTEP_RELATIVE_MAX,
            "HHG spectrum fails time-step convergence")

    diagnostics = {
        "ideal_cartesian_amplitudes": [float(value) for value in amplitudes],
        "ideal_harmonic_morphologies": morphologies,
        "ideal_current_derivative_relative_rms": current_derivative_mismatch,
        "impulse_transverse_ratio": impulse_transverse_ratio,
        "laser_transverse_ratio": laser_transverse_ratio,
        "nonlinear_third_harmonic_power_ratio": nonlinear_third_harmonic_power_ratio,
    }
    print(json.dumps(diagnostics, sort_keys=True))
    print("PASS ideal genuine-Si64 overlapping-Wannier response/HHG matrix")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
