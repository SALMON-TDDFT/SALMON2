#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

BOHR_PER_ANG = 1.0 / 0.529177249
GPA_PER_HA_BOHR3 = 29421.02648


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare primitive nonorthogonal kinetic stress against finite differences."
    )
    parser.add_argument("--base-dir", required=True, help="Directory containing m*/p* subdirectories")
    parser.add_argument("--mode", choices=("scale", "shear"), required=True)
    parser.add_argument("--tolerance", type=float, default=0.1, help="Allowed absolute mismatch in GPa")
    return parser.parse_args()


def read_first_data_row(path: Path) -> list[float]:
    for line in path.read_text().splitlines():
        if line and not line.startswith("#"):
            return [float(x) for x in line.split()]
    raise RuntimeError(f"no data rows found in {path}")


def read_conventional_a_ang(case_dir: Path) -> float:
    input_text = (case_dir / "inputfile").read_text()
    match = re.search(r"al_vec1\s*=\s*([0-9.E+\-]+)d0,\s*([0-9.E+\-]+)d0,\s*([0-9.E+\-]+)d0", input_text)
    if not match:
        raise RuntimeError(f"al_vec1 not found in {case_dir / 'inputfile'}")
    return 2.0 * float(match.group(1))


def read_sector_component(info_path: Path, sector: str, component: str) -> float:
    table_pattern = re.compile(
        rf"^\s*{re.escape(sector)}\s+([0-9.E+\-]+)\s+([0-9.E+\-]+)\s+([0-9.E+\-]+)\s+([0-9.E+\-]+)\s+([0-9.E+\-]+)\s+([0-9.E+\-]+)\s+([0-9.E+\-]+)\s*$",
        re.M,
    )
    match = table_pattern.search(info_path.read_text())
    if not match:
        raise RuntimeError(f"sector row '{sector}' not found in {info_path}")
    columns = {
        "xx": 1,
        "yy": 2,
        "zz": 3,
        "xy": 4,
        "yz": 5,
        "xz": 6,
        "P": 7,
    }
    return float(match.group(columns[component]))


def scale_fd(rows: list[tuple[float, Path]]) -> tuple[float, float]:
    if len(rows) != 3:
        raise RuntimeError("scale mode expects exactly 3 points")
    entries = []
    for _, case_dir in rows:
        vals = read_first_data_row(next(case_dir.glob("*_stress_energy.data")))
        e_kin = vals[2]  # col 3: E_kin [Ha]
        p_kin = vals[12]  # col 13: P_kin [GPa]
        a_ang = read_conventional_a_ang(case_dir)
        entries.append((a_ang, e_kin, p_kin))
    entries.sort(key=lambda item: item[0])
    a_bohr = [a * BOHR_PER_ANG for a, _, _ in entries]
    e_kin = [e for _, e, _ in entries]
    analytic = entries[1][2]
    dE_da = (e_kin[2] - e_kin[0]) / (a_bohr[2] - a_bohr[0])
    dV_da = 3.0 * a_bohr[1] ** 2 / 4.0
    fd = -(dE_da / dV_da) * GPA_PER_HA_BOHR3
    return analytic, fd


def shear_fd(rows: list[tuple[float, Path]]) -> tuple[float, float]:
    if len(rows) != 3:
        raise RuntimeError("shear mode expects exactly 3 points")
    entries = []
    for delta, case_dir in rows:
        vals = read_first_data_row(next(case_dir.glob("*_stress_energy.data")))
        e_kin = vals[2]  # col 3: E_kin [Ha]
        info_path = next(case_dir.glob("*_info.data"))
        xy_kin = read_sector_component(info_path, "Kinetic", "xy")
        a_ang = read_conventional_a_ang(case_dir)
        entries.append((delta, e_kin, xy_kin, info_path, a_ang))
    entries.sort(key=lambda item: item[0])
    a_center_bohr = entries[1][4] * BOHR_PER_ANG
    v_prim = a_center_bohr ** 3 / 4.0
    dE_ddelta = (entries[2][1] - entries[0][1]) / (entries[2][0] - entries[0][0])
    fd = dE_ddelta / v_prim * GPA_PER_HA_BOHR3
    analytic = entries[1][2]
    return analytic, fd


def collect_cases(base_dir: Path, mode: str) -> list[tuple[float, Path]]:
    names = [("m001", -0.005), ("p000", 0.0), ("p001", 0.005)] if mode == "scale" else [("m005", -0.005), ("p000", 0.0), ("p005", 0.005)]
    rows = []
    for name, value in names:
        case_dir = base_dir / name
        if not case_dir.is_dir():
            raise RuntimeError(f"missing case directory: {case_dir}")
        rows.append((value, case_dir))
    return rows


def main() -> int:
    args = parse_args()
    base_dir = Path(args.base_dir)
    rows = collect_cases(base_dir, args.mode)
    if args.mode == "scale":
        analytic, fd = scale_fd(rows)
        label = "P_kin"
    else:
        analytic, fd = shear_fd(rows)
        label = "sigma_kin_xy"
    diff = analytic - fd
    print(f"{label}: analytic={analytic:.6f} GPa fd={fd:.6f} GPa diff={diff:.6f} GPa")
    return 0 if abs(diff) <= args.tolerance else 1


if __name__ == "__main__":
    sys.exit(main())
