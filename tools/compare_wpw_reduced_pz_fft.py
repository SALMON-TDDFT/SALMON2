#!/usr/bin/env python3
"""Compare production/reduced WPW Pz diagnostic time series.

Input files are written by the DG-WPW reduced diagnostic:
  1: step, 2: time[a.u.], 3: keep_n, 4: Pz_prod,
  5: Pz_reduced, 6: Pz_diff, 7: rel_Pz_diff
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


AU_TO_EV = 27.211386245988


def load_pz(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 7:
        raise ValueError(f"{path}: expected at least 7 columns, got {data.shape[1]}")
    return data


def fft_amp(signal: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    centered = signal - signal[0]
    window = np.hanning(centered.size) if centered.size > 2 else np.ones_like(centered)
    spec = np.fft.rfft(centered * window)
    freq_au = 2.0 * np.pi * np.fft.rfftfreq(centered.size, d=dt)
    return freq_au * AU_TO_EV, np.abs(spec)


def summarize(path: Path) -> None:
    data = load_pz(path)
    time = data[:, 1]
    keep_n = int(data[0, 2])
    pz_prod = data[:, 3]
    pz_reduced = data[:, 4]
    pz_diff = data[:, 5]
    if time.size > 1:
        dt = float(np.median(np.diff(time)))
    else:
        dt = 1.0
    energy_ev, amp_prod = fft_amp(pz_prod, dt)
    _, amp_reduced = fft_amp(pz_reduced, dt)
    _, amp_diff = fft_amp(pz_diff, dt)
    if energy_ev.size > 1:
        peak_idx = 1 + int(np.argmax(amp_prod[1:]))
    else:
        peak_idx = 0
    amp_delta = np.abs(amp_prod - amp_reduced)
    print(f"file={path}")
    print(f"  keep_n={keep_n} n={time.size} dt_au={dt:.16e}")
    print(f"  max_abs_Pz_diff={np.max(np.abs(pz_diff)):.16e}")
    print(f"  rms_abs_Pz_diff={np.sqrt(np.mean(pz_diff**2)):.16e}")
    print(f"  fft_prod_peak_eV={energy_ev[peak_idx]:.16e}")
    print(f"  fft_prod_peak_amp={amp_prod[peak_idx]:.16e}")
    print(f"  fft_reduced_peak_amp_at_prod_peak={amp_reduced[peak_idx]:.16e}")
    print(f"  fft_amp_diff_at_prod_peak={amp_delta[peak_idx]:.16e}")
    print(f"  fft_max_amp_diff={np.max(amp_delta):.16e}")
    print(f"  fft_max_diff_signal_amp={np.max(amp_diff):.16e}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path)
    args = parser.parse_args()
    for path in args.files:
        summarize(path)


if __name__ == "__main__":
    main()
