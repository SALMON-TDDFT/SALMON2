#!/usr/bin/env python3
"""Plot DG-Wannier+BPW Ecut/Sperp sweep CSV files as dependency-free SVGs."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


DEFAULT_METRICS = (
    "np",
    "fsum_avg",
    "c_comm_sum",
    "min_sperp",
    "max_sperp",
    "herm",
    "selected_nraw",
    "selected_g2_max",
)

COLORS = (
    "#1f77b4",
    "#d62728",
    "#2ca02c",
    "#9467bd",
    "#ff7f0e",
    "#17becf",
    "#8c564b",
    "#7f7f7f",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot BPW Ecut sweep CSV diagnostics.")
    parser.add_argument("csv", nargs="+", help="CSV files from dg_bpw_ecut_sweep.py")
    parser.add_argument("--out-dir", default="bpw_sweep_plots", help="Output directory")
    parser.add_argument("--metrics", nargs="+", default=list(DEFAULT_METRICS))
    parser.add_argument(
        "--log-y",
        nargs="*",
        default=("fsum_avg", "c_comm_sum", "herm"),
        help="Metrics to plot on log10 y scale.",
    )
    parser.add_argument("--width", type=int, default=900)
    parser.add_argument("--height", type=int, default=560)
    return parser.parse_args()


def read_rows(paths: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        with open(path, newline="") as csv_file:
            for row in csv.DictReader(csv_file):
                row = dict(row)
                row["_source"] = path
                rows.append(row)
    return rows


def as_float(row: dict[str, str], key: str) -> float | None:
    try:
        value = row.get(key, "")
        if value == "":
            return None
        parsed = float(value)
        if math.isnan(parsed) or math.isinf(parsed):
            return None
        return parsed
    except ValueError:
        return None


def group_rows(rows: list[dict[str, str]], metric: str):
    grouped = defaultdict(list)
    for row in rows:
        x = as_float(row, "ecut")
        y = as_float(row, metric)
        if x is None or y is None:
            continue
        case_label = row.get("case_label", "") or "case"
        sperp_tol = row.get("sperp_tol", "") or "default"
        grouped[(case_label, sperp_tol)].append((x, y))
    for key in grouped:
        grouped[key] = sorted(grouped[key])
    return grouped


def nice_range(values: list[float], log_y: bool) -> tuple[float, float]:
    if not values:
        return 0.0, 1.0
    if log_y:
        positive = [value for value in values if value > 0.0]
        if not positive:
            return -1.0, 1.0
        lo = math.floor(math.log10(min(positive)))
        hi = math.ceil(math.log10(max(positive)))
        if lo == hi:
            hi = lo + 1.0
        return lo, hi
    lo = min(values)
    hi = max(values)
    if lo == hi:
        pad = 1.0 if lo == 0.0 else abs(lo) * 0.1
        return lo - pad, hi + pad
    pad = 0.06 * (hi - lo)
    return lo - pad, hi + pad


def svg_text(x: float, y: float, text: str, size: int = 13, anchor: str = "middle") -> str:
    escaped = (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}">{escaped}</text>'


def plot_metric(rows: list[dict[str, str]], metric: str, out_path: Path, width: int, height: int, log_y: bool) -> None:
    grouped = group_rows(rows, metric)
    margin_l, margin_r, margin_t, margin_b = 90, 30, 52, 76
    plot_w = width - margin_l - margin_r
    plot_h = height - margin_t - margin_b

    xs = [x for points in grouped.values() for x, _ in points]
    ys = [y for points in grouped.values() for _, y in points]
    if not xs or not ys:
        return
    x_min, x_max = min(xs), max(xs)
    if x_min == x_max:
        x_min -= 0.5
        x_max += 0.5
    y_min, y_max = nice_range(ys, log_y)

    def sx(x: float) -> float:
        return margin_l + (x - x_min) / (x_max - x_min) * plot_w

    def sy(y: float) -> float:
        if log_y:
            if y <= 0.0:
                return margin_t + plot_h
            y = math.log10(y)
        return margin_t + (y_max - y) / (y_max - y_min) * plot_h

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:-apple-system,BlinkMacSystemFont,"Helvetica Neue",Arial,sans-serif;fill:#222}'
        '.axis{stroke:#333;stroke-width:1.2}.grid{stroke:#ddd;stroke-width:1}.line{fill:none;stroke-width:2.2}'
        '.pt{stroke:white;stroke-width:1}</style>',
        svg_text(width / 2, 28, f"{metric} vs Ecut", 18),
        f'<line class="axis" x1="{margin_l}" y1="{margin_t + plot_h}" x2="{margin_l + plot_w}" y2="{margin_t + plot_h}"/>',
        f'<line class="axis" x1="{margin_l}" y1="{margin_t}" x2="{margin_l}" y2="{margin_t + plot_h}"/>',
    ]

    for i in range(6):
        t = i / 5.0
        x_val = x_min + t * (x_max - x_min)
        x_pos = sx(x_val)
        parts.append(f'<line class="grid" x1="{x_pos:.1f}" y1="{margin_t}" x2="{x_pos:.1f}" y2="{margin_t + plot_h}"/>')
        parts.append(svg_text(x_pos, margin_t + plot_h + 24, f"{x_val:g}", 12))

        y_val = y_min + t * (y_max - y_min)
        y_pos = margin_t + (1.0 - t) * plot_h
        label = f"1e{int(y_val)}" if log_y else f"{y_val:.3g}"
        parts.append(f'<line class="grid" x1="{margin_l}" y1="{y_pos:.1f}" x2="{margin_l + plot_w}" y2="{y_pos:.1f}"/>')
        parts.append(svg_text(margin_l - 10, y_pos + 4, label, 12, "end"))

    parts.append(svg_text(margin_l + plot_w / 2, height - 22, "Ecut [a.u.]", 14))
    y_label = f"log10({metric})" if log_y else metric
    parts.append(f'<text x="22" y="{margin_t + plot_h / 2:.1f}" font-size="14" text-anchor="middle" '
                 f'transform="rotate(-90 22 {margin_t + plot_h / 2:.1f})">{y_label}</text>')

    legend_x = margin_l + 10
    legend_y = margin_t + 18
    for idx, (key, points) in enumerate(sorted(grouped.items())):
        color = COLORS[idx % len(COLORS)]
        coords = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in points)
        parts.append(f'<polyline class="line" points="{coords}" stroke="{color}"/>')
        for x, y in points:
            parts.append(f'<circle class="pt" cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="4" fill="{color}"/>')
        label = f"{key[0]} / Sperp={key[1]}"
        ly = legend_y + idx * 20
        parts.append(f'<line x1="{legend_x}" y1="{ly - 4}" x2="{legend_x + 22}" y2="{ly - 4}" stroke="{color}" stroke-width="2.2"/>')
        parts.append(svg_text(legend_x + 30, ly, label, 12, "start"))

    parts.append("</svg>")
    out_path.write_text("\n".join(parts))


def main() -> int:
    args = parse_args()
    rows = read_rows(args.csv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_metrics = set(args.log_y)
    for metric in args.metrics:
        plot_metric(rows, metric, out_dir / f"{metric}.svg", args.width, args.height, metric in log_metrics)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
