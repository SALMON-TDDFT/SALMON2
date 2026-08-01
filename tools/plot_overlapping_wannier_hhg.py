#!/usr/bin/env python3
"""Plot overlapping-Wannier HHG spectra derived strictly from polarization."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_spectrum(root: Path, axis: str) -> tuple[np.ndarray, Path]:
    case = root / f"laser-hhg-{axis}"
    summary = json.loads((case / "hhg-summary.json").read_text())
    if summary.get("spectrum_source") != "polarization":
        raise RuntimeError(f"{case}: HHG spectrum is not polarization-derived")
    path = case / "hhg-spectrum.tsv"
    spectrum = np.loadtxt(path)
    if spectrum.ndim != 2 or spectrum.shape[1] < 4 or not np.isfinite(spectrum).all():
        raise RuntimeError(f"{path}: invalid HHG spectrum")
    return spectrum, path


def plot_panel(axis: plt.Axes, curves: list[tuple[np.ndarray, str]], maximum_order: float) -> None:
    positive = np.concatenate([curve[:, 3][curve[:, 3] > 0.0] for curve, _ in curves])
    if positive.size == 0:
        raise RuntimeError("HHG spectra contain no positive power")
    floor = max(float(np.min(positive)) * 0.1, float(np.max(positive)) * 1.0e-14)
    for spectrum, label in curves:
        selected = spectrum[:, 2] <= maximum_order
        axis.semilogy(spectrum[selected, 2], np.maximum(spectrum[selected, 3], floor), label=label)
    for order in range(1, int(maximum_order) + 1):
        axis.axvline(order, color="0.86", linewidth=0.6, zorder=0)
    axis.set_xlim(0.0, maximum_order)
    axis.set_xlabel("Harmonic order")
    axis.set_ylabel(r"Polarization HHG power  $\omega^4 |P(\omega)|^2$")
    axis.grid(True, which="both", axis="y", alpha=0.2)
    axis.legend(frameon=False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ideal-root", type=Path, required=True)
    parser.add_argument("--displaced-root", type=Path)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--maximum-order", type=float, default=15.0)
    args = parser.parse_args()
    if args.maximum_order <= 1.0:
        raise RuntimeError("maximum harmonic order must exceed one")

    ideal: dict[str, np.ndarray] = {}
    hashes: dict[str, str] = {}
    for axis_name in "xyz":
        ideal[axis_name], path = load_spectrum(args.ideal_root.resolve(strict=True), axis_name)
        hashes[f"ideal_{axis_name}"] = digest(path)

    panels = ["ideal_xyz"]
    panel_count = 1 if args.displaced_root is None else 2
    figure, axes = plt.subplots(panel_count, 1, figsize=(8.2, 4.8 * panel_count), squeeze=False)
    plot_panel(axes[0, 0], [(ideal[name], f"ideal {name}") for name in "xyz"], args.maximum_order)
    axes[0, 0].set_title("Ideal structure: Cartesian polarization spectra")

    if args.displaced_root is not None:
        displaced, path = load_spectrum(args.displaced_root.resolve(strict=True), "x")
        hashes["displaced_x"] = digest(path)
        plot_panel(axes[1, 0], [(ideal["x"], "ideal x"), (displaced, "fixed-displaced x")],
                   args.maximum_order)
        axes[1, 0].set_title("Exact instantaneous symmetry comparison")
        panels.append("ideal_displaced_x")

    figure.tight_layout()
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    png = args.output_prefix.with_suffix(".png")
    pdf = args.output_prefix.with_suffix(".pdf")
    figure.savefig(png, dpi=180)
    figure.savefig(pdf)
    plt.close(figure)
    evidence = {
        "spectrum_source": "polarization",
        "panels": panels,
        "maximum_order": args.maximum_order,
        "input_sha256": hashes,
        "png_sha256": digest(png),
        "pdf_sha256": digest(pdf),
    }
    args.output_prefix.with_suffix(".json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

