#!/usr/bin/env python3
"""Analyze overlapping-Wannier RT spectra from polarization, never current."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


HARTREE_EV = 27.211386245988
MAGIC = "# SALMON_OW_COEFFICIENT_RT_OBSERVABLES_V1"
COLUMNS = "# step time Ex Ey Ez Px Py Pz Jx Jy Jz"
AXIS_COLUMN = {"x": 5, "y": 6, "z": 7}


def load_series(path: Path) -> tuple[np.ndarray, tuple[str, ...]]:
    lines = path.read_text().splitlines()
    if not lines or lines[0] != MAGIC or COLUMNS not in lines:
        raise ValueError(f"{path}: invalid observable header")
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] != 11 or data.shape[0] < 3 or not np.isfinite(data).all():
        raise ValueError(f"{path}: expected at least three finite 11-column samples")
    steps = data[:, 0]
    if not np.array_equal(steps, np.arange(len(data)) + steps[0]) or np.any(steps != np.floor(steps)):
        raise ValueError(f"{path}: steps must be consecutive integers")
    differences = np.diff(data[:, 1])
    if differences[0] <= 0 or not np.allclose(differences, differences[0], rtol=1e-12, atol=1e-14):
        raise ValueError(f"{path}: time samples must be strictly increasing and uniform")
    column_index = lines.index(COLUMNS)
    return data, tuple(lines[: column_index + 1])


def window_values(kind: str, count: int, dt: float, damping_time: float | None) -> np.ndarray:
    if kind == "hann":
        return np.hanning(count)
    if damping_time is None or not np.isfinite(damping_time) or damping_time <= 0:
        raise ValueError("exponential window requires a finite positive --damping-time")
    return np.exp(-np.arange(count) * dt / damping_time)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("linear", "hhg"))
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--background", type=Path, required=True)
    parser.add_argument("--axis", choices=tuple(AXIS_COLUMN), required=True)
    parser.add_argument("--window", choices=("hann", "exponential"), required=True)
    parser.add_argument("--damping-time", type=float)
    parser.add_argument("--impulse-area", type=float)
    parser.add_argument("--carrier-ev", type=float)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    driven, driven_header = load_series(args.input)
    background, background_header = load_series(args.background)
    if driven_header != background_header:
        raise ValueError("driven and field-off observable headers/provenance must match")
    if driven.shape != background.shape or not np.array_equal(driven[:, :2], background[:, :2]):
        raise ValueError("driven and field-off series must have identical step/time grids")
    dt = float(driven[1, 1] - driven[0, 1])
    count = len(driven)
    polarization = driven[:, AXIS_COLUMN[args.axis]] - background[:, AXIS_COLUMN[args.axis]]
    window = window_values(args.window, count, dt, args.damping_time)
    transform = dt * np.fft.rfft(polarization * window)
    omega = 2.0 * np.pi * np.fft.rfftfreq(count, d=dt)
    energy_ev = omega * HARTREE_EV
    abs_squared = np.abs(transform) ** 2
    resolution = 2.0 * np.pi / (count * dt)
    summary: dict[str, object] = {
        "axis": args.axis,
        "dt_au": dt,
        "field_off_subtracted": True,
        "frequency_resolution_au": resolution,
        "nyquist_au": np.pi / dt,
        "sample_count": count,
        "spectrum_source": "polarization",
        "window": args.window,
        "window_coherent_gain": float(np.mean(window)),
    }
    if args.mode == "linear":
        if args.impulse_area is None or not np.isfinite(args.impulse_area) or args.impulse_area == 0:
            raise ValueError("linear mode requires finite nonzero --impulse-area")
        response = transform / args.impulse_area
        table = np.column_stack([omega, energy_ev, response.real, response.imag, np.abs(response) ** 2])
        header = "omega_au energy_ev response_re response_im abs_response_squared"
        summary.update(mode="linear", impulse_area_au=args.impulse_area,
                       normalization="P_omega/integrated_impulse_field")
    else:
        if args.carrier_ev is None or not np.isfinite(args.carrier_ev) or args.carrier_ev <= 0:
            raise ValueError("hhg mode requires finite positive --carrier-ev")
        harmonic_order = energy_ev / args.carrier_ev
        intensity = omega ** 4 * abs_squared
        table = np.column_stack([omega, energy_ev, harmonic_order, intensity, abs_squared])
        header = "omega_au energy_ev harmonic_order omega4_abs_p_squared abs_p_squared"
        summary.update(mode="hhg", carrier_ev=args.carrier_ev,
                       formula="omega^4*abs(P_omega)^2")
    np.savetxt(args.output, table, fmt="%.17e", header=header)
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
