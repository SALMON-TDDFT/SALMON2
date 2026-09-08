#!/usr/bin/env python3
"""Convert dg_bpw_auto_report.dat files into a compact CSV table."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


SUMMARY_RE = re.compile(r"(\w+)=\s*([^\s]+)")

FIELDS = [
    "case",
    "path",
    "candidate_BPW_number",
    "selected_raw_BPW_number",
    "selected_perp_BPW_number",
    "selected_G2_max",
    "selected_Emax",
    "Sperp_cutoff_used",
    "singular_value_min_selected",
    "singular_value_max_selected",
    "singular_value_condition",
    "projectability_min_selected",
    "projectability_max_rejected",
    "projectability_gap",
    "hermiticity_max",
    "C_comm_sum",
    "fsum_proxy_wp_avg",
    "fsum_proxy_total_avg",
    "warnings",
    "recommendation",
]


def parse_report(path: Path) -> dict[str, str]:
    data: dict[str, str] = {"path": str(path)}
    warning_lines = 0
    in_warnings = False

    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line == "warnings:":
            in_warnings = True
            continue
        if line == "end_warnings":
            in_warnings = False
            continue
        if in_warnings:
            warning_lines += 1
            continue
        if line.startswith("SUMMARY "):
            for key, value in SUMMARY_RE.findall(line):
                data[f"summary_{key}"] = value
            continue
        if "=" in line:
            key, value = line.split("=", 1)
            data[key.strip()] = value.strip()

    data.setdefault("warnings", data.get("summary_warnings", str(warning_lines)))
    return data


def case_name(path: Path) -> str:
    if path.name == "dg_bpw_auto_report.dat":
        return path.parent.name
    return path.stem.replace("_dg_bpw_auto_report", "")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs="+", type=Path)
    parser.add_argument("--output", "-o", type=Path, default=Path("-"))
    args = parser.parse_args()

    rows = []
    for report in args.reports:
        row = parse_report(report)
        row["case"] = case_name(report)
        rows.append(row)

    out = None
    try:
        if str(args.output) == "-":
            out = None
            writer = csv.DictWriter(__import__("sys").stdout, fieldnames=FIELDS, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        else:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
                writer.writeheader()
                writer.writerows(rows)
    finally:
        if out is not None:
            out.close()


if __name__ == "__main__":
    main()
