#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plotter for SALMON SBE real-time output files.

Drop this script into a calculation directory (next to SYSNAME_sbe_rt.data,
SYSNAME_sbe_rt_energy.data, SYSNAME_sbe_nex.data, SYSNAME_sbe_nex_k.data, ...)
and run it. It scans the directory for these files, plots:
  - SYSNAME_sbe_rt_energy.data : total energy vs time
  - SYSNAME_sbe_nex.data       : number of excited electrons/holes vs time
  - SYSNAME_sbe_nex_k.data     : per-k Houston-basis population, plotted as
                                 three 2D heatmap slices (kx-ky, kx-kz, ky-kz)
                                 -- one PNG per saved time, all times encoded
                                 in the file names.

No interactive windows are opened (Agg backend); everything is saved as PNG
into an output directory.
"""

import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

plt.switch_backend('Agg')  # Non-interactive backend


def parse_header(header_line):
    """Extract column names, ignoring units in brackets."""
    pattern = r'\d+:([^\[\s]+)(?:\[[^\]]*\])?'
    return re.findall(pattern, header_line)


def find_header(filepath):
    """Find header line with column numbers. Returns (header_text, line_index)."""
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('#') and re.search(r'#\s*1\s*:\s*\S', line):
                return line.strip(), i
    raise ValueError("Header line with column numbers not found.")


def load_columns(filepath):
    """Load a simple whitespace-separated SBE .data file into a 2D array."""
    header_line, _ = find_header(filepath)
    column_names = parse_header(header_line)

    rows = []
    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            try:
                rows.append([float(x) for x in stripped.split()])
            except ValueError:
                continue

    data = np.array(rows) if rows else np.empty((0, len(column_names)))
    return column_names, data


def plot_xy(time, values, time_name, col_name, output_path, dpi=150):
    if len(time) == 0:
        print(f"  (skip) no data for {col_name}")
        return

    plt.figure(figsize=(10, 6))
    plt.plot(time, values, linewidth=1.0)
    plt.xlabel(time_name)
    plt.ylabel(col_name)
    plt.title(f'{col_name} vs {time_name}')
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.tight_layout()

    safe_name = re.sub(r'[^\w\-]', '_', col_name)
    safe_time = re.sub(r'[^\w\-]', '_', time_name)
    out_file = output_path / f'{safe_name}_vs_{safe_time}.png'
    plt.savefig(out_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"  saved {out_file.name}")


def plot_energy_and_nex(filepath, output_dir, dpi=150):
    print(f"Processing {filepath.name} ...")
    cols, data = load_columns(filepath)
    if data.size == 0:
        print("  (skip) no data lines found")
        return
    time_name = cols[0]
    time = data[:, 0]
    for j in range(1, len(cols)):
        plot_xy(time, data[:, j], time_name, cols[j], output_dir, dpi=dpi)


def read_nex_k_blocks(filepath):
    """
    Parse SYSNAME_sbe_nex_k.data: a sequence of blocks, each starting with a
    "# t = <value> <unit>" comment line, followed by `nk` data lines
    "ik kx ky kz population".
    Returns a list of (t_value, t_unit, kpoints[nk,3], population[nk]).
    """
    blocks = []
    t_value = None
    t_unit = ''
    kx, ky, kz, pop = [], [], [], []

    time_re = re.compile(r'#\s*t\s*=\s*([-\d.eEdD+]+)\s*(\S*)')

    def flush():
        nonlocal t_value, t_unit, kx, ky, kz, pop
        if t_value is not None and kx:
            blocks.append((t_value, t_unit,
                           np.array([kx, ky, kz]).T,
                           np.array(pop)))
        t_value, t_unit = None, ''
        kx, ky, kz, pop = [], [], [], []

    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            m = time_re.match(stripped)
            if m:
                flush()
                t_value = float(m.group(1))
                t_unit = m.group(2)
                continue
            if stripped.startswith('#'):
                continue
            parts = stripped.split()
            if len(parts) < 5:
                continue
            try:
                kx.append(float(parts[1]))
                ky.append(float(parts[2]))
                kz.append(float(parts[3]))
                pop.append(float(parts[4]))
            except ValueError:
                continue
    flush()
    return blocks


def heatmap_slice(ax, k_a, k_b, values, label_a, label_b, title):
    """Bin scattered (k_a, k_b, value) points onto a regular grid and draw a heatmap."""
    if len(k_a) == 0:
        ax.set_title(title + " (no data)")
        return

    ua = np.unique(np.round(k_a, 10))
    ub = np.unique(np.round(k_b, 10))
    if len(ua) > 1 and len(ub) > 1:
        grid = np.full((len(ub), len(ua)), np.nan)
        ia = np.searchsorted(ua, np.round(k_a, 10))
        ib = np.searchsorted(ub, np.round(k_b, 10))
        for i in range(len(values)):
            grid[ib[i], ia[i]] = values[i]
        im = ax.imshow(grid, origin='lower', aspect='auto',
                       extent=[ua.min(), ua.max(), ub.min(), ub.max()],
                       cmap='viridis')
    else:
        im = ax.scatter(k_a, k_b, c=values, cmap='viridis', s=20)

    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_xlabel(label_a)
    ax.set_ylabel(label_b)
    ax.set_title(title)


def plot_nex_k(filepath, output_dir, dpi=150):
    print(f"Processing {filepath.name} ...")
    blocks = read_nex_k_blocks(filepath)
    if not blocks:
        print("  (skip) no data blocks found")
        return

    for t_value, t_unit, kpoints, pop in blocks:
        fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

        heatmap_slice(axes[0], kpoints[:, 0], kpoints[:, 1], pop,
                      'kx [a.u.]', 'ky [a.u.]', 'population_lcb: kx-ky (kz~0 slice)')
        heatmap_slice(axes[1], kpoints[:, 0], kpoints[:, 2], pop,
                      'kx [a.u.]', 'kz [a.u.]', 'population_lcb: kx-kz (ky~0 slice)')
        heatmap_slice(axes[2], kpoints[:, 1], kpoints[:, 2], pop,
                      'ky [a.u.]', 'kz [a.u.]', 'population_lcb: ky-kz (kx~0 slice)')

        fig.suptitle(f'Houston-basis lowest-conduction-band population, t = {t_value:.6f} {t_unit}')
        fig.tight_layout(rect=[0, 0, 1, 0.94])

        out_file = output_dir / f'nex_k_t_{t_value:.6f}{t_unit}.png'
        fig.savefig(out_file, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        print(f"  saved {out_file.name}")


def main():
    parser = argparse.ArgumentParser(
        description='Plot SALMON SBE real-time output files (energy, nex, per-k population).')
    parser.add_argument('-i', '--input-dir', default='.', help='Directory with SYSNAME_sbe_*.data files')
    parser.add_argument('-o', '--output', default='sbe_plots', help='Output directory for PNGs')
    parser.add_argument('--dpi', type=int, default=150, help='Image resolution')
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    energy_files = sorted(input_dir.glob('*_sbe_rt_energy.data'))
    nex_files = sorted(input_dir.glob('*_sbe_nex.data'))
    nex_k_files = sorted(input_dir.glob('*_sbe_nex_k.data'))

    for f in energy_files:
        plot_energy_and_nex(f, output_dir, dpi=args.dpi)
    for f in nex_files:
        plot_energy_and_nex(f, output_dir, dpi=args.dpi)
    for f in nex_k_files:
        plot_nex_k(f, output_dir, dpi=args.dpi)

    if not (energy_files or nex_files or nex_k_files):
        print(f"No SYSNAME_sbe_rt_energy.data / SYSNAME_sbe_nex.data / SYSNAME_sbe_nex_k.data "
              f"files found in {input_dir.resolve()}")
        return

    print(f"\nDone. Output directory: {output_dir.resolve()}")


if __name__ == '__main__':
    main()
