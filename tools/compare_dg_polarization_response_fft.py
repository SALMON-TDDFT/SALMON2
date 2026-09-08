#!/usr/bin/env python3
"""Compare polarization-derived response spectra from P(t) time-series files.

The comparison reads DG length-gauge polarization files directly, or integrates
SALMON matter current from *_rt.data as a fallback for conventional TDDFT. It
uses the same preprocessing for all inputs and reports peak/RMS differences.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


AU_TO_EV = 27.211386245988


def cumulative_trapezoid(time: np.ndarray, values: np.ndarray) -> np.ndarray:
    out = np.zeros_like(values)
    if time.size <= 1:
        return out
    dt = np.diff(time)[:, None]
    out[1:, :] = np.cumsum(0.5 * (values[1:, :] + values[:-1, :]) * dt, axis=0)
    return out


def load_polarization(path: Path, source: str) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if source == "dg_polarization":
        if data.shape[1] < 4:
            raise ValueError(f"{path}: expected at least 4 columns in polarization data, got {data.shape[1]}")
        return data[:, 0], data[:, 1:4]
    if source == "rt_current_integral":
        if data.shape[1] < 16:
            raise ValueError(f"{path}: expected at least 16 columns in *_rt.data, got {data.shape[1]}")
        time = data[:, 0]
        # SALMON current is dP/dt up to convention. Use the same sign for all
        # current-integrated references; sign checks belong to Stage 1g.
        return time, cumulative_trapezoid(time, data[:, 13:16])
    if source == "rt_current_integral_minus":
        time, pol = load_polarization(path, "rt_current_integral")
        return time, -pol
    raise ValueError(f"unknown source: {source}")


def response_spectrum(
    time: np.ndarray,
    signal: np.ndarray,
    axis: int,
    damping: float,
    pad_factor: int,
    skip: int,
    skip_time_au: float | None,
    tmax_au: float | None,
    detrend: str,
) -> tuple[np.ndarray, np.ndarray]:
    if skip_time_au is not None:
        if skip_time_au < 0.0:
            raise ValueError(f"invalid skip_time_au={skip_time_au}")
        skip = int(np.searchsorted(time, time[0] + skip_time_au, side="left"))
    if skip < 0 or skip >= time.size:
        raise ValueError(f"invalid skip={skip} for n={time.size}")
    end = time.size
    if tmax_au is not None:
        if tmax_au <= 0.0:
            raise ValueError(f"invalid tmax_au={tmax_au}")
        end = int(np.searchsorted(time, time[skip] + tmax_au, side="right"))
    if end - skip < 2:
        raise ValueError("need at least two samples after skip/tmax selection")
    t = time[skip:end] - time[skip]
    y = signal[skip:end, axis].copy()
    if y.size < 2:
        raise ValueError("need at least two samples after skip")
    y -= y[0]
    if detrend == "linear":
        slope, intercept = np.polyfit(t, y, 1)
        y -= slope * t + intercept
    elif detrend == "mean":
        y -= np.mean(y)
    elif detrend != "none":
        raise ValueError(f"invalid detrend mode: {detrend}")
    if damping > 0.0:
        y *= np.exp(-damping * t)
    if y.size > 2:
        y *= np.hanning(y.size)
    nfft = int(2 ** np.ceil(np.log2(max(1, y.size * max(1, pad_factor)))))
    dt = float(np.median(np.diff(t)))
    spec = np.fft.rfft(y, n=nfft) * dt
    omega_ev = 2.0 * np.pi * np.fft.rfftfreq(nfft, d=dt) * AU_TO_EV
    return omega_ev, np.abs(spec)


def parse_meta(text: str | None) -> dict[str, str]:
    meta = {
        "case": "",
        "observable": "Pz",
        "field_type": "laser_pulse",
        "gauge": "",
        "method": "",
        "basis": "",
        "fragment": "",
        "propagator": "",
        "dt_au": "",
        "Tmax_au": "",
        "field_abs": "",
        "axis": "",
        "block": "",
    }
    if not text:
        return meta
    for item in text.split(","):
        if not item.strip():
            continue
        if "=" not in item:
            raise ValueError(f"metadata item must be key=value: {item}")
        key, value = item.split("=", 1)
        key = key.strip()
        if key not in meta:
            raise ValueError(f"unknown metadata key: {key}")
        meta[key] = value.strip()
    return meta


def parse_harmonic_orders(text: str) -> list[int]:
    if not text.strip():
        return []
    orders = []
    for item in text.split(","):
        if not item.strip():
            continue
        order = int(item)
        if order <= 0:
            raise ValueError(f"harmonic order must be positive: {item}")
        orders.append(order)
    return orders


def harmonic_error_text(
    energy: np.ndarray,
    ref_amp: np.ndarray,
    amp: np.ndarray,
    fundamental_ev: float,
    orders: list[int],
    half_width_ev: float,
) -> str:
    if fundamental_ev <= 0.0 or not orders:
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


def format_table_row(row: dict[str, str | float | int]) -> str:
    columns = (
        "case",
        "method",
        "basis",
        "fragment",
        "propagator",
        "dt_au",
        "Tmax_au",
        "field_abs",
        "axis",
        "block",
        "observable",
        "peak_shift_eV",
        "peak_height_error",
        "rel_peak_height_error",
        "rms_error",
        "rel_rms_error",
        "low_energy_rel_rms_error",
        "hhg_window_rms_error",
        "cutoff_shift_eV",
        "harmonic_rel_intensity_errors",
        "bad",
    )
    out = []
    for col in columns:
        value = row.get(col, "")
        if isinstance(value, float):
            out.append(f"{value:.16e}")
        else:
            out.append(str(value))
    return "\t".join(out)


def summarize(
    ref_path: Path,
    cmp_paths: list[Path],
    ref_source: str,
    cmp_sources: list[str],
    ref_meta: dict[str, str],
    cmp_meta: list[dict[str, str]],
    axis: int,
    damping: float,
    pad_factor: int,
    skip: int,
    skip_time_au: float | None,
    tmax_au: float | None,
    detrend: str,
    emin: float,
    emax: float,
    hhg_emin: float,
    hhg_emax: float,
    cutoff_threshold_rel: float,
    fundamental_ev: float,
    harmonic_orders: list[int],
    harmonic_half_width_ev: float,
    table: bool,
) -> None:
    ref_time, ref_signal = load_polarization(ref_path, ref_source)
    energy, ref_amp = response_spectrum(ref_time, ref_signal, axis, damping, pad_factor, skip, skip_time_au, tmax_au, detrend)
    mask = np.ones_like(energy, dtype=bool)
    if emin >= 0.0:
        mask &= energy >= emin
    if emax > 0.0:
        mask &= energy <= emax
    if np.count_nonzero(mask) == 0:
        raise ValueError("empty energy comparison window")
    masked_idx = np.where(mask)[0]
    if masked_idx.size > 1:
        peak_idx = masked_idx[int(np.argmax(ref_amp[masked_idx]))]
    else:
        peak_idx = masked_idx[0]
    ref_norm = max(float(np.sqrt(np.mean(ref_amp[mask] ** 2))), 1.0e-300)
    hhg_mask = mask.copy()
    if hhg_emin >= 0.0:
        hhg_mask &= energy >= hhg_emin
    if hhg_emax > 0.0:
        hhg_mask &= energy <= hhg_emax
    if np.count_nonzero(hhg_mask) == 0:
        raise ValueError("empty HHG energy comparison window")
    hhg_ref_norm = max(float(np.sqrt(np.mean(ref_amp[hhg_mask] ** 2))), 1.0e-300)
    cutoff_ref_level = cutoff_threshold_rel * max(float(ref_amp[peak_idx]), 1.0e-300)
    cutoff_ref_candidates = np.where(hhg_mask & (ref_amp >= cutoff_ref_level))[0]
    cutoff_ref = float(energy[cutoff_ref_candidates[-1]]) if cutoff_ref_candidates.size else 0.0
    if table:
        print(
            "case\tmethod\tbasis\tfragment\tpropagator\tdt_au\tTmax_au\tfield_abs\taxis\tblock\t"
            "observable\t"
            "peak_shift_eV\tpeak_height_error\t"
            "rel_peak_height_error\trms_error\trel_rms_error\tlow_energy_rel_rms_error\t"
            "hhg_window_rms_error\tcutoff_shift_eV\tharmonic_rel_intensity_errors\tbad"
        )
    else:
        print(f"reference={ref_path}")
        if any(ref_meta.values()):
            print("reference_meta=" + ",".join(f"{key}={value}" for key, value in ref_meta.items()))
        print(f"axis={axis} n={ref_time.size} dt_au={np.median(np.diff(ref_time)):.16e}")
        print(f"reference_source={ref_source}")
        print(f"damping_au={damping:.16e} pad_factor={pad_factor} skip={skip}")
        print(f"detrend={detrend}")
        if skip_time_au is not None:
            print(f"skip_time_au={skip_time_au:.16e}")
        if tmax_au is not None:
            print(f"tmax_au={tmax_au:.16e}")
        print(f"ref_peak_eV={energy[peak_idx]:.16e}")
        print(f"ref_peak_amp={ref_amp[peak_idx]:.16e}")
    for path, source, meta in zip(cmp_paths, cmp_sources, cmp_meta, strict=True):
        time, signal = load_polarization(path, source)
        energy_cmp, amp_cmp = response_spectrum(time, signal, axis, damping, pad_factor, skip, skip_time_au, tmax_au, detrend)
        energy_window = energy[mask]
        if energy_window[0] < energy_cmp[0] - 1.0e-12 or energy_window[-1] > energy_cmp[-1] + 1.0e-12:
            raise ValueError(f"{path}: candidate energy grid does not cover comparison energy window")
        amp = np.interp(energy, energy_cmp, amp_cmp)
        diff = amp - ref_amp
        abs_diff = np.abs(diff)
        rms = float(np.sqrt(np.mean(diff[mask] ** 2)))
        rel_rms = rms / ref_norm
        hhg_rms = float(np.sqrt(np.mean(diff[hhg_mask] ** 2))) / hhg_ref_norm
        peak_shift = 0.0
        cmp_peak_idx = masked_idx[int(np.argmax(amp[masked_idx]))]
        peak_shift = float(energy[cmp_peak_idx] - energy[peak_idx])
        peak_height_error = float(amp[peak_idx] - ref_amp[peak_idx])
        rel_peak_height_error = peak_height_error / max(float(abs(ref_amp[peak_idx])), 1.0e-300)
        low = masked_idx[: max(1, min(8, masked_idx.size))]
        low_ref = max(float(np.sqrt(np.mean(ref_amp[low] ** 2))), 1.0e-300)
        low_err = float(np.sqrt(np.mean(diff[low] ** 2))) / low_ref
        cutoff_cmp_level = cutoff_threshold_rel * max(float(amp[cmp_peak_idx]), 1.0e-300)
        cutoff_cmp_candidates = np.where(hhg_mask & (amp >= cutoff_cmp_level))[0]
        cutoff_cmp = float(energy[cutoff_cmp_candidates[-1]]) if cutoff_cmp_candidates.size else 0.0
        cutoff_shift = cutoff_cmp - cutoff_ref
        harmonic_errors = harmonic_error_text(
            energy,
            ref_amp,
            amp,
            fundamental_ev,
            harmonic_orders,
            harmonic_half_width_ev,
        )
        if table:
            if not meta["Tmax_au"] and tmax_au is not None:
                meta["Tmax_au"] = f"{tmax_au:.16e}"
            if not meta["axis"]:
                meta["axis"] = "xyz"[axis]
            row = {
                **meta,
                "case": meta["case"] or path.stem,
                "peak_shift_eV": peak_shift,
                "peak_height_error": peak_height_error,
                "rel_peak_height_error": rel_peak_height_error,
                "rms_error": rms,
                "rel_rms_error": rel_rms,
                "low_energy_rel_rms_error": low_err,
                "hhg_window_rms_error": hhg_rms,
                "cutoff_shift_eV": cutoff_shift,
                "harmonic_rel_intensity_errors": harmonic_errors,
                "bad": "F",
            }
            print(format_table_row(row))
        else:
            print(f"candidate={path}")
            if any(meta.values()):
                print("  candidate_meta=" + ",".join(f"{key}={value}" for key, value in meta.items()))
            print(f"  candidate_dt_au={np.median(np.diff(time)):.16e}")
            print(f"  peak_eV={energy[cmp_peak_idx]:.16e}")
            print(f"  peak_position_shift_eV={peak_shift:.16e}")
            print(f"  peak_height_error={peak_height_error:.16e}")
            print(f"  rel_peak_height_error={rel_peak_height_error:.16e}")
            print(f"  rms_spectrum_error={rms:.16e}")
            print(f"  rel_rms_spectrum_error={rel_rms:.16e}")
            print(f"  max_abs_spectrum_error={float(np.max(abs_diff[mask])):.16e}")
            print(f"  low_energy_rel_rms_error={low_err:.16e}")
            print(f"  hhg_window_rms_error={hhg_rms:.16e}")
            print(f"  cutoff_shift_eV={cutoff_shift:.16e}")
            if harmonic_errors:
                print(f"  harmonic_rel_intensity_errors={harmonic_errors}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("candidates", nargs="+", type=Path)
    parser.add_argument("--reference-source", choices=("dg_polarization", "rt_current_integral", "rt_current_integral_minus"), default="dg_polarization")
    parser.add_argument("--candidate-source", action="append", default=[], choices=("dg_polarization", "rt_current_integral", "rt_current_integral_minus"), help="source for one candidate; repeat in candidate order")
    parser.add_argument("--reference-meta", default="", help="comma-separated key=value metadata")
    parser.add_argument(
        "--candidate-meta",
        action="append",
        default=[],
        help="comma-separated key=value metadata for one candidate; repeat in candidate order",
    )
    parser.add_argument("--axis", choices=("x", "y", "z"), default="z")
    parser.add_argument("--damping", type=float, default=0.0, help="exponential damping in a.u.^-1")
    parser.add_argument("--pad-factor", type=int, default=4)
    parser.add_argument("--skip", type=int, default=0)
    parser.add_argument("--skip-time-au", type=float, default=None, help="physical start time in a.u.; overrides --skip")
    parser.add_argument("--tmax-au", type=float, default=None, help="physical window length in a.u. after skip")
    parser.add_argument("--detrend", choices=("none", "mean", "linear"), default="none")
    parser.add_argument("--emin", type=float, default=0.0, help="minimum energy in eV")
    parser.add_argument("--emax", type=float, default=0.0, help="maximum energy in eV; <=0 means no cap")
    parser.add_argument("--hhg-emin", type=float, default=-1.0, help="minimum HHG-window energy in eV; <0 reuses --emin")
    parser.add_argument("--hhg-emax", type=float, default=0.0, help="maximum HHG-window energy in eV; <=0 reuses --emax")
    parser.add_argument("--cutoff-threshold-rel", type=float, default=1.0e-4, help="relative spectral threshold for cutoff estimate")
    parser.add_argument("--fundamental-ev", type=float, default=0.0, help="laser fundamental energy for harmonic-order diagnostics")
    parser.add_argument("--harmonic-orders", default="", help="comma-separated harmonic orders, for example 3,5,7")
    parser.add_argument("--harmonic-half-width-ev", type=float, default=0.25, help="half-width around each harmonic peak")
    parser.add_argument("--table", action="store_true", help="print compact TSV table")
    args = parser.parse_args()
    axis = {"x": 0, "y": 1, "z": 2}[args.axis]
    if args.candidate_meta and len(args.candidate_meta) != len(args.candidates):
        raise SystemExit("--candidate-meta must be repeated once per candidate")
    if args.candidate_source and len(args.candidate_source) != len(args.candidates):
        raise SystemExit("--candidate-source must be repeated once per candidate")
    candidate_meta = [parse_meta(text) for text in args.candidate_meta]
    while len(candidate_meta) < len(args.candidates):
        candidate_meta.append(parse_meta(None))
    candidate_sources = list(args.candidate_source)
    while len(candidate_sources) < len(args.candidates):
        candidate_sources.append("dg_polarization")
    summarize(
        args.reference,
        args.candidates,
        args.reference_source,
        candidate_sources,
        parse_meta(args.reference_meta),
        candidate_meta,
        axis,
        args.damping,
        args.pad_factor,
        args.skip,
        args.skip_time_au,
        args.tmax_au,
        args.detrend,
        args.emin,
        args.emax,
        args.hhg_emin if args.hhg_emin >= 0.0 else args.emin,
        args.hhg_emax if args.hhg_emax > 0.0 else args.emax,
        args.cutoff_threshold_rel,
        args.fundamental_ev,
        parse_harmonic_orders(args.harmonic_orders),
        args.harmonic_half_width_ev,
        args.table,
    )


if __name__ == "__main__":
    main()
