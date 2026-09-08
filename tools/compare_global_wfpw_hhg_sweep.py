#!/usr/bin/env python3
"""Compare a Global WF+PW basis sweep against a Full TDDFT reference.

The manifest is a TSV/CSV table with at least a ``path`` column. Optional
columns such as ``basis_id``, ``WF_count``, ``PW_count``, ``PW_cutoff_or_shell``,
``projector_count``, ``dt``, ``propagator_kind``, ``observable_source``,
``gauge``, and ``volume_normalization`` are copied to the output table.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable

import numpy as np

from compare_dg_polarization_response_fft import (
    load_polarization,
    parse_harmonic_orders,
    response_spectrum,
)


def read_manifest(path: Path) -> list[dict[str, str]]:
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open(newline="") as fp:
        rows = list(csv.DictReader(fp, delimiter=delimiter))
    if not rows:
        raise ValueError(f"{path}: empty manifest")
    required = ("path", "basis_id", "WF_count", "PW_count", "PW_cutoff_or_shell")
    missing = [key for key in required if key not in rows[0]]
    if missing:
        raise ValueError(f"{path}: manifest missing required columns: {','.join(missing)}")
    return rows


def select_signal_window(
    time: np.ndarray,
    signal: np.ndarray,
    axis: int,
    skip_time_au: float | None,
    tmax_au: float | None,
    detrend: str,
) -> tuple[np.ndarray, np.ndarray]:
    start = 0
    if skip_time_au is not None:
        start = int(np.searchsorted(time, time[0] + skip_time_au, side="left"))
    end = time.size
    if tmax_au is not None:
        end = int(np.searchsorted(time, time[start] + tmax_au, side="right"))
    if end - start < 2:
        raise ValueError("need at least two samples in time-domain comparison window")
    t = time[start:end] - time[start]
    y = signal[start:end, axis].copy()
    y -= y[0]
    if detrend == "linear":
        slope, intercept = np.polyfit(t, y, 1)
        y -= slope * t + intercept
    elif detrend == "mean":
        y -= np.mean(y)
    elif detrend != "none":
        raise ValueError(f"invalid detrend mode: {detrend}")
    return t, y


def rel_time_rms(
    ref_time: np.ndarray,
    ref_signal: np.ndarray,
    cmp_time: np.ndarray,
    cmp_signal: np.ndarray,
    axis: int,
    skip_time_au: float | None,
    tmax_au: float | None,
    detrend: str,
) -> float:
    tref, yref = select_signal_window(ref_time, ref_signal, axis, skip_time_au, tmax_au, detrend)
    tcmp, ycmp = select_signal_window(cmp_time, cmp_signal, axis, skip_time_au, tmax_au, detrend)
    if tref[0] < tcmp[0] - 1.0e-12 or tref[-1] > tcmp[-1] + 1.0e-12:
        mask = (tref >= tcmp[0] - 1.0e-12) & (tref <= tcmp[-1] + 1.0e-12)
        tref = tref[mask]
        yref = yref[mask]
    yinterp = np.interp(tref, tcmp, ycmp)
    ref_norm = max(float(np.sqrt(np.mean(yref**2))), 1.0e-300)
    return float(np.sqrt(np.mean((yinterp - yref) ** 2))) / ref_norm


def cutoff_energy(energy: np.ndarray, amp: np.ndarray, mask: np.ndarray, level: float) -> float:
    idx = np.where(mask & (amp >= level))[0]
    return float(energy[idx[-1]]) if idx.size else 0.0


def harmonic_peak_positions(
    energy: np.ndarray,
    amp: np.ndarray,
    fundamental_ev: float,
    orders: Iterable[int],
    half_width_ev: float,
) -> str:
    if fundamental_ev <= 0.0:
        return ""
    values = []
    for order in orders:
        center = order * fundamental_ev
        mask = np.abs(energy - center) <= half_width_ev
        if not np.any(mask):
            idx = int(np.argmin(np.abs(energy - center)))
            values.append(f"{order}:{energy[idx]:.8e}")
            continue
        local = np.where(mask)[0]
        peak = local[int(np.argmax(amp[local]))]
        values.append(f"{order}:{energy[peak]:.8e}")
    return ";".join(values)


def harmonic_intensity_errors(
    energy: np.ndarray,
    ref_amp: np.ndarray,
    amp: np.ndarray,
    fundamental_ev: float,
    orders: Iterable[int],
    half_width_ev: float,
) -> str:
    if fundamental_ev <= 0.0:
        return ""
    values = []
    for order in orders:
        center = order * fundamental_ev
        mask = np.abs(energy - center) <= half_width_ev
        if not np.any(mask):
            idx = int(np.argmin(np.abs(energy - center)))
            mask = np.zeros_like(energy, dtype=bool)
            mask[idx] = True
        ref_peak = max(float(np.max(ref_amp[mask])), 1.0e-300)
        cmp_peak = float(np.max(amp[mask]))
        values.append(f"{order}:{(cmp_peak - ref_peak) / ref_peak:.8e}")
    return ";".join(values)


def harmonic_intensity_error_map(
    energy: np.ndarray,
    ref_amp: np.ndarray,
    amp: np.ndarray,
    fundamental_ev: float,
    orders: Iterable[int],
    half_width_ev: float,
) -> dict[int, float]:
    if fundamental_ev <= 0.0:
        return {}
    values = {}
    for order in orders:
        center = order * fundamental_ev
        mask = np.abs(energy - center) <= half_width_ev
        if not np.any(mask):
            idx = int(np.argmin(np.abs(energy - center)))
            mask = np.zeros_like(energy, dtype=bool)
            mask[idx] = True
        ref_peak = max(float(np.max(ref_amp[mask])), 1.0e-300)
        cmp_peak = float(np.max(amp[mask]))
        values[order] = (cmp_peak - ref_peak) / ref_peak
    return values


def parse_optional_float(value: str) -> float | None:
    if not value:
        return None
    return float(value)


def fmt(value: str | float) -> str:
    if isinstance(value, float):
        return f"{value:.16e}"
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--reference-source", choices=("dg_polarization", "rt_current_integral", "rt_current_integral_minus"), default="rt_current_integral")
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--axis", choices=("x", "y", "z"), default="z")
    parser.add_argument("--damping", type=float, default=0.0)
    parser.add_argument("--pad-factor", type=int, default=4)
    parser.add_argument("--skip-time-au", type=float, default=None)
    parser.add_argument("--tmax-au", type=float, default=None)
    parser.add_argument("--detrend", choices=("none", "mean", "linear"), default="linear")
    parser.add_argument("--emin", type=float, default=0.0)
    parser.add_argument("--emax", type=float, default=0.0)
    parser.add_argument("--hhg-emin", type=float, default=-1.0)
    parser.add_argument("--hhg-emax", type=float, default=0.0)
    parser.add_argument("--cutoff-threshold-rel", type=float, default=1.0e-4)
    parser.add_argument("--fundamental-ev", type=float, default=0.0)
    parser.add_argument("--harmonic-orders", default="")
    parser.add_argument("--harmonic-half-width-ev", type=float, default=0.25)
    parser.add_argument("--peak-shift-tol-ev", type=float, default=-1.0)
    parser.add_argument("--cutoff-shift-tol-ev", type=float, default=-1.0)
    parser.add_argument("--hhg-window-rms-tol", type=float, default=-1.0)
    parser.add_argument("--norm-charge-drift-tol", type=float, default=-1.0)
    args = parser.parse_args()

    axis = {"x": 0, "y": 1, "z": 2}[args.axis]
    orders = parse_harmonic_orders(args.harmonic_orders)

    ref_time, ref_signal = load_polarization(args.reference, args.reference_source)
    energy, ref_amp = response_spectrum(
        ref_time,
        ref_signal,
        axis,
        args.damping,
        args.pad_factor,
        0,
        args.skip_time_au,
        args.tmax_au,
        args.detrend,
    )
    mask = np.ones_like(energy, dtype=bool)
    if args.emin >= 0.0:
        mask &= energy >= args.emin
    if args.emax > 0.0:
        mask &= energy <= args.emax
    if np.count_nonzero(mask) == 0:
        raise ValueError("empty energy comparison window")
    hhg_mask = mask.copy()
    hhg_emin = args.hhg_emin if args.hhg_emin >= 0.0 else args.emin
    hhg_emax = args.hhg_emax if args.hhg_emax > 0.0 else args.emax
    if hhg_emin >= 0.0:
        hhg_mask &= energy >= hhg_emin
    if hhg_emax > 0.0:
        hhg_mask &= energy <= hhg_emax
    if np.count_nonzero(hhg_mask) == 0:
        raise ValueError("empty HHG energy comparison window")

    masked_idx = np.where(mask)[0]
    ref_peak = masked_idx[int(np.argmax(ref_amp[masked_idx]))]
    ref_norm = max(float(np.sqrt(np.mean(ref_amp[mask] ** 2))), 1.0e-300)
    hhg_ref_norm = max(float(np.sqrt(np.mean(ref_amp[hhg_mask] ** 2))), 1.0e-300)
    cutoff_ref = cutoff_energy(
        energy,
        ref_amp,
        hhg_mask,
        args.cutoff_threshold_rel * max(float(ref_amp[ref_peak]), 1.0e-300),
    )

    harmonic_columns = tuple(f"H{order:02d}_rel_intensity_error" for order in orders)
    columns = (
        "case",
        "method",
        "basis",
        "basis_id",
        "WF_count",
        "PW_count",
        "PW_cutoff_or_shell",
        "projector_count",
        "dt",
        "propagator_kind",
        "observable_source",
        "gauge",
        "volume_normalization",
        "field_abs",
        "axis",
        "block",
        "peak_shift_eV",
        "rel_peak_height_error",
        "rel_rms_error",
        "hhg_window_rms_error",
        "HHG_peak_positions",
        "HHG_plateau_metric",
        "HHG_cutoff",
        "cutoff_shift_eV",
        "harmonic_intensities",
        "polarization_RMS",
        "current_RMS",
        "energy_drift",
        "norm_charge_drift",
        *harmonic_columns,
        "bad",
    )
    print("\t".join(columns))

    for row in read_manifest(args.manifest):
        out = {key: row.get(key, "") for key in columns}
        out["axis"] = out["axis"] or args.axis
        try:
            path = Path(row["path"])
            source = row.get("source") or row.get("observable_source") or "dg_polarization"
            time, signal = load_polarization(path, source)
            ecmp, acmp = response_spectrum(
                time,
                signal,
                axis,
                args.damping,
                args.pad_factor,
                0,
                args.skip_time_au,
                args.tmax_au,
                args.detrend,
            )
            window_energy = energy[mask]
            if window_energy[0] < ecmp[0] - 1.0e-12 or window_energy[-1] > ecmp[-1] + 1.0e-12:
                raise ValueError("candidate energy grid does not cover comparison window")
            amp = np.interp(energy, ecmp, acmp)
            diff = amp - ref_amp
            cmp_peak = masked_idx[int(np.argmax(amp[masked_idx]))]
            peak_shift = float(energy[cmp_peak] - energy[ref_peak])
            rel_peak_height = float(amp[ref_peak] - ref_amp[ref_peak]) / max(float(abs(ref_amp[ref_peak])), 1.0e-300)
            hhg_rms = float(np.sqrt(np.mean(diff[hhg_mask] ** 2))) / hhg_ref_norm
            rms = float(np.sqrt(np.mean(diff[mask] ** 2))) / ref_norm
            cutoff_cmp = cutoff_energy(
                energy,
                amp,
                hhg_mask,
                args.cutoff_threshold_rel * max(float(amp[cmp_peak]), 1.0e-300),
            )
            out.update(
                {
                    "peak_shift_eV": peak_shift,
                    "rel_peak_height_error": rel_peak_height,
                    "rel_rms_error": rms,
                    "hhg_window_rms_error": hhg_rms,
                    "HHG_peak_positions": harmonic_peak_positions(
                        energy,
                        amp,
                        args.fundamental_ev,
                        orders,
                        args.harmonic_half_width_ev,
                    ),
                    "HHG_plateau_metric": float(np.mean(amp[hhg_mask]) / max(float(np.mean(ref_amp[hhg_mask])), 1.0e-300)),
                    "HHG_cutoff": cutoff_cmp,
                    "cutoff_shift_eV": cutoff_cmp - cutoff_ref,
                    "harmonic_intensities": harmonic_intensity_errors(
                        energy,
                        ref_amp,
                        amp,
                        args.fundamental_ev,
                        orders,
                        args.harmonic_half_width_ev,
                    ),
                    "polarization_RMS": rel_time_rms(
                        ref_time,
                        ref_signal,
                        time,
                        signal,
                        axis,
                        args.skip_time_au,
                        args.tmax_au,
                        args.detrend,
                    ),
                    "bad": "F",
                }
            )
            for order, value in harmonic_intensity_error_map(
                energy,
                ref_amp,
                amp,
                args.fundamental_ev,
                orders,
                args.harmonic_half_width_ev,
            ).items():
                out[f"H{order:02d}_rel_intensity_error"] = value
            bad_reasons = []
            if args.peak_shift_tol_ev >= 0.0 and abs(peak_shift) > args.peak_shift_tol_ev:
                bad_reasons.append("peak_shift")
            if args.cutoff_shift_tol_ev >= 0.0 and abs(cutoff_cmp - cutoff_ref) > args.cutoff_shift_tol_ev:
                bad_reasons.append("cutoff_shift")
            if args.hhg_window_rms_tol >= 0.0 and hhg_rms > args.hhg_window_rms_tol:
                bad_reasons.append("hhg_window_rms")
            norm_charge_drift = parse_optional_float(row.get("norm_charge_drift", ""))
            if (
                args.norm_charge_drift_tol >= 0.0
                and norm_charge_drift is not None
                and abs(norm_charge_drift) > args.norm_charge_drift_tol
            ):
                bad_reasons.append("norm_charge_drift")
            if bad_reasons:
                out["bad"] = "T:" + ",".join(bad_reasons)
        except Exception as exc:  # noqa: BLE001 - sweep output should identify bad cases.
            out["bad"] = f"T:{exc}"
        print("\t".join(fmt(out.get(col, "")) for col in columns))


if __name__ == "__main__":
    main()
