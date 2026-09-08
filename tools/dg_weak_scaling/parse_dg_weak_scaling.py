#!/usr/bin/env python3
"""Summarize RT-DG weak-scaling logs.

The script is intentionally log-format tolerant: missing fields are left blank
so that partially failed Fugaku runs can still be compared.
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics
from pathlib import Path


FIELDS = [
    "case",
    "status",
    "mpi_ranks",
    "ranks_per_fragment_subgroup",
    "n_frag",
    "states_per_fragment",
    "n_mat",
    "nstate_tot",
    "dist_eig_n",
    "dist_eig_nkeep",
    "dist_eig_solved",
    "total_time_s",
    "momentum_self_s",
    "momentum_grad_s",
    "rk_iter_count",
    "rk_mean_s",
    "stage_mean_s",
    "deriv_mean_s",
    "postres_per_occ_last",
    "full_rel_initial",
    "full_rel_final",
]


PATTERNS = {
    "mpi_ranks": re.compile(r"MPI parallelization:\s*(\d+)\s+processes"),
    "subgroup": re.compile(r"MPI ranks per fragment subgroup:\s*(\d+)"),
    "n_frag": re.compile(r"Number of fragments:\s*(\d+)"),
    "states_per_fragment": re.compile(r"States per fragment:\s*(\d+)"),
    "core_eig": re.compile(r"cap=\s*(\d+)\s+n_mat=\s*(\d+)\s+nstate_tot=\s*(\d+)"),
    "dist_setup": re.compile(r"solve(?: setup)?:\s*n=\s*(\d+)\s+nkeep=\s*(\d+)"),
    "dist_final_solved": re.compile(r"\[DG-DIST-EIG-FINAL\]\s+distributed_solved=\s*([TF])"),
    "dist_initial": re.compile(r"\[DG-DIST-EIG-INITIAL\]\s+full_rel=\s*([+\-\d.Ee]+)"),
    "dist_final": re.compile(r"\[DG-DIST-EIG-FINAL\]\s+full_rel=\s*([+\-\d.Ee]+)"),
    "total_time": re.compile(r"total calculation time,\s*([+\-\d.Ee]+)"),
    "momentum": re.compile(
        r"momentum timing:\s+halo_exchange=\s*([+\-\d.Ee]+)\s+grad=\s*([+\-\d.Ee]+)\s+self=\s*([+\-\d.Ee]+)"
    ),
    "rk_total": re.compile(r"rk timing:\s*itt=\s*\d+.*?total=\s*([+\-\d.Ee]+)"),
    "rk_stage_update": re.compile(r"\[DG-RK\]\s+stage\s+\d+\s+density/H update done time=\s*([+\-\d.Ee]+)"),
    "rk_deriv": re.compile(r"\[DG-RK\]\s+stage\s+\d+\s+derivative done time=\s*([+\-\d.Ee]+)"),
    "postres": re.compile(r"postResPerOcc=.*?\n?\s*([+\-\d.Ee]+)\s*$"),
}


def parse_float(text: str) -> float | None:
    try:
        return float(text.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def mean(values: list[float]) -> str:
    if not values:
        return ""
    return f"{statistics.fmean(values):.6g}"


def parse_log(path: Path) -> dict[str, str]:
    text = path.read_text(errors="replace")
    row: dict[str, str] = {field: "" for field in FIELDS}
    row["case"] = path.parent.name if path.name == "run.log" else path.stem
    row["status"] = "fatal" if re.search(r"\[FATAL\]|SIGSEGV|Program received signal|error\(s\) in input", text) else "ok"

    for key, pattern in PATTERNS.items():
        if key in {"momentum", "rk_total", "rk_stage_update", "rk_deriv", "postres"}:
            continue
        matches = list(pattern.finditer(text))
        if not matches:
            continue
        m = matches[-1]
        if key == "core_eig":
            row["states_per_fragment"] = row["states_per_fragment"] or m.group(1)
            row["n_mat"] = m.group(2)
            row["nstate_tot"] = m.group(3)
        elif key == "dist_setup":
            row["dist_eig_n"] = m.group(1)
            row["dist_eig_nkeep"] = m.group(2)
        elif key == "dist_final_solved":
            row["dist_eig_solved"] = m.group(1)
        elif key == "dist_initial":
            row["full_rel_initial"] = m.group(1)
        elif key == "dist_final":
            row["full_rel_final"] = m.group(1)
        else:
            field = {
                "subgroup": "ranks_per_fragment_subgroup",
                "total_time": "total_time_s",
            }.get(key, key)
            row[field] = m.group(1)

    momentum_matches = list(PATTERNS["momentum"].finditer(text))
    if momentum_matches:
        m = momentum_matches[-1]
        row["momentum_grad_s"] = m.group(2)
        row["momentum_self_s"] = m.group(3)

    rk_total = [v for v in (parse_float(m.group(1)) for m in PATTERNS["rk_total"].finditer(text)) if v is not None]
    stage = [v for v in (parse_float(m.group(1)) for m in PATTERNS["rk_stage_update"].finditer(text)) if v is not None]
    deriv = [v for v in (parse_float(m.group(1)) for m in PATTERNS["rk_deriv"].finditer(text)) if v is not None]
    if rk_total:
        row["rk_iter_count"] = str(len(rk_total))
        row["rk_mean_s"] = mean(rk_total)
    elif stage or deriv:
        row["rk_iter_count"] = str(max(len(stage), len(deriv)))
        row["rk_mean_s"] = mean([s + d for s, d in zip(stage, deriv)])
    row["stage_mean_s"] = mean(stage)
    row["deriv_mean_s"] = mean(deriv)

    postres = list(PATTERNS["postres"].finditer(text))
    if postres:
        row["postres_per_occ_last"] = postres[-1].group(1)

    return row


def print_markdown(rows: list[dict[str, str]]) -> None:
    cols = [
        "case",
        "status",
        "mpi_ranks",
        "ranks_per_fragment_subgroup",
        "n_frag",
        "n_mat",
        "nstate_tot",
        "dist_eig_solved",
        "total_time_s",
        "rk_mean_s",
        "stage_mean_s",
        "deriv_mean_s",
        "full_rel_final",
    ]
    print("| " + " | ".join(cols) + " |")
    print("| " + " | ".join(["---"] * len(cols)) + " |")
    for row in rows:
        print("| " + " | ".join(row.get(col, "") for col in cols) + " |")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", nargs="+", type=Path, help="run.log files or directories containing run.log")
    parser.add_argument("--format", choices=("csv", "markdown"), default="csv")
    args = parser.parse_args()

    paths: list[Path] = []
    for item in args.logs:
        if item.is_dir():
            paths.extend(sorted(item.glob("**/run.log")))
        else:
            paths.append(item)
    rows = [parse_log(path) for path in paths]

    if args.format == "markdown":
      print_markdown(rows)
    else:
      writer = csv.DictWriter(__import__("sys").stdout, fieldnames=FIELDS)
      writer.writeheader()
      writer.writerows(rows)


if __name__ == "__main__":
    main()
