#!/usr/bin/env python3
"""Run DG-Wannier+BPW shell_ecut diagnostics and collect a CSV summary.

The script intentionally keeps BPW selection as a plain cutoff sweep.  It
does not optimize commutators or f-sums; those quantities are parsed only as
diagnostics from SALMON stdout.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
from pathlib import Path


DEFAULT_ECUTS = ("0.5", "1.0", "1.5", "2.0", "3.0", "4.0")


SELECT_RE = re.compile(
    r"\[DG-BPW-SELECT\]\s+mode=(?P<mode>\S+).*?"
    r"selected_shells=\s*(?P<selected_shells>\d+).*?"
    r"selected_nraw=\s*(?P<selected_nraw>\d+).*?"
    r"selected_G2_max=\s*(?P<selected_g2_max>\d+).*?"
    r"selected_Emax=\s*(?P<selected_emax>[-+0-9.Ee]+)"
)
SHELL_RE = re.compile(
    r"shell_id=\s*(?P<shell_id>\d+)\s+G2=\s*(?P<g2>\d+)\s+"
    r"E_shell=\s*(?P<e_shell>[-+0-9.Ee]+)\s+nraw=\s*(?P<nraw>\d+)\s+"
    r"np=\s*(?P<np>\d+)\s+fWP_avg_occ=\s*(?P<fwp_avg_occ>[-+0-9.Ee]+)"
)
FWP_COMM_RE = re.compile(
    r"fWP_xyz_occ=\s*(?P<fwp_x>[-+0-9.Ee]+)\s+"
    r"(?P<fwp_y>[-+0-9.Ee]+)\s+(?P<fwp_z>[-+0-9.Ee]+)\s+"
    r"C_comm_xyz=\s*(?P<c_x>[-+0-9.Ee]+)\s+"
    r"(?P<c_y>[-+0-9.Ee]+)\s+(?P<c_z>[-+0-9.Ee]+)"
)
FSUM_RE = re.compile(
    r"fsum_xyz=\s*(?P<fsum_x>[-+0-9.Ee]+)\s+"
    r"(?P<fsum_y>[-+0-9.Ee]+)\s+(?P<fsum_z>[-+0-9.Ee]+)\s+"
    r"min_Sperp=\s*(?P<min_sperp>[-+0-9.Ee]+)\s+"
    r"max_Sperp=\s*(?P<max_sperp>[-+0-9.Ee]+)"
)
HERMIT_RE = re.compile(r"hermiticity max\|Z-Z\^H\|=\s*(?P<herm>[-+0-9.Ee]+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sweep SALMON_DG_BPW_ECUT for DG-Wannier+BPW diagnostics."
    )
    parser.add_argument("--salmon", required=True, help="Path to the salmon executable.")
    parser.add_argument("--input", required=True, help="SALMON input file.")
    parser.add_argument("--workdir", default=".", help="Run directory.")
    parser.add_argument("--output", default="bpw_ecut_sweep.csv", help="CSV output path.")
    parser.add_argument("--log-dir", default="", help="Optional directory for per-Ecut stdout logs.")
    parser.add_argument("--ecut", nargs="+", default=list(DEFAULT_ECUTS), help="Ecut values in a.u.")
    parser.add_argument(
        "--sperp-tol",
        nargs="+",
        default=[""],
        help="Optional SALMON_DG_MIXED_SPERP_TOL values. Empty default uses SALMON default.",
    )
    parser.add_argument("--np", type=int, default=8, help="MPI ranks.")
    parser.add_argument("--omp", type=int, default=2, help="OMP_NUM_THREADS.")
    parser.add_argument(
        "--case-label",
        default="",
        help="Optional label for the Wannier set/input case, written to the CSV.",
    )
    parser.add_argument("--extra-env", action="append", default=[], metavar="KEY=VALUE")
    parser.add_argument(
        "--allow-diagnostic-stop",
        action="store_true",
        help=(
            "Accept a non-zero SALMON exit after BPW shell diagnostics were printed. "
            "This is useful for nt=0 diagnostic sweeps that only need mixed-Z setup output."
        ),
    )
    return parser.parse_args()


def parse_key_values(items: list[str]) -> dict[str, str]:
    env: dict[str, str] = {}
    for item in items:
        if "=" not in item:
            raise SystemExit(f"--extra-env must be KEY=VALUE, got {item!r}")
        key, value = item.split("=", 1)
        env[key] = value
    return env


def parse_summary(stdout: str, ecut: str, sperp_tol: str, case_label: str) -> dict[str, str]:
    selected: dict[str, str] = {
        "case_label": case_label,
        "ecut": ecut,
        "sperp_tol": sperp_tol or "default",
    }
    select_match = SELECT_RE.search(stdout)
    if select_match:
        selected.update(select_match.groupdict())

    shell_rows: list[dict[str, str]] = []
    pending: dict[str, str] | None = None
    for line in stdout.splitlines():
        shell_match = SHELL_RE.search(line)
        if shell_match:
            if pending is not None:
                shell_rows.append(pending)
            pending = shell_match.groupdict()
            continue
        if pending is None:
            continue
        fwp_match = FWP_COMM_RE.search(line)
        if fwp_match:
            pending.update(fwp_match.groupdict())
            continue
        fsum_match = FSUM_RE.search(line)
        if fsum_match:
            pending.update(fsum_match.groupdict())
            continue
    if pending is not None:
        shell_rows.append(pending)

    selected_shell = selected.get("selected_shells")
    row = None
    if selected_shell:
        if selected_shell != "0":
            row = next((item for item in shell_rows if item.get("shell_id") == selected_shell), None)
        else:
            row = {
                "shell_id": "0",
                "g2": selected.get("selected_g2_max", "0"),
                "e_shell": selected.get("selected_emax", "0.0"),
                "nraw": selected.get("selected_nraw", "0"),
                "np": "0",
                "fwp_avg_occ": "0.0",
                "fwp_x": "0.0",
                "fwp_y": "0.0",
                "fwp_z": "0.0",
                "c_x": "0.0",
                "c_y": "0.0",
                "c_z": "0.0",
                "fsum_x": "0.0",
                "fsum_y": "0.0",
                "fsum_z": "0.0",
            }
    if row is None and shell_rows:
        row = shell_rows[-1]
    if row:
        selected.update(row)

    herm_match = HERMIT_RE.search(stdout)
    if herm_match:
        selected.update(herm_match.groupdict())

    add_derived_columns(selected)
    return selected


def add_derived_columns(row: dict[str, str]) -> None:
    def get_float(key: str) -> float | None:
        try:
            return float(row[key])
        except (KeyError, ValueError):
            return None

    c_vals = [get_float("c_x"), get_float("c_y"), get_float("c_z")]
    if all(value is not None for value in c_vals):
        row["c_comm_sum"] = f"{sum(value for value in c_vals if value is not None):.16e}"

    fsum_vals = [get_float("fsum_x"), get_float("fsum_y"), get_float("fsum_z")]
    if all(value is not None for value in fsum_vals):
        row["fsum_avg"] = f"{sum(value for value in fsum_vals if value is not None) / 3.0:.16e}"


def main() -> int:
    args = parse_args()
    workdir = Path(args.workdir)
    salmon = Path(args.salmon)
    input_path = Path(args.input)
    output = Path(args.output)
    log_dir = Path(args.log_dir) if args.log_dir else None
    if log_dir:
        log_dir.mkdir(parents=True, exist_ok=True)

    extra_env = parse_key_values(args.extra_env)
    rows: list[dict[str, str]] = []
    for sperp_tol in args.sperp_tol:
        for ecut in args.ecut:
            env = os.environ.copy()
            env.update(extra_env)
            env["OMP_NUM_THREADS"] = str(args.omp)
            env["SALMON_DG_MIXED_Z"] = "1"
            env["SALMON_DG_BPW_SELECT_MODE"] = "shell_ecut"
            env["SALMON_DG_BPW_ECUT"] = str(ecut)
            env["SALMON_DG_BPW_COMM_DIAG"] = "1"
            env["SALMON_DG_BPW_FSUM_DIAG"] = "1"
            if sperp_tol:
                env["SALMON_DG_MIXED_SPERP_TOL"] = str(sperp_tol)

            with input_path.open("rb") as input_file:
                result = subprocess.run(
                    ["mpirun", "-np", str(args.np), str(salmon)],
                    cwd=workdir,
                    env=env,
                    stdin=input_file,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            stdout = result.stdout.decode(errors="replace")
            if log_dir:
                safe_ecut = str(ecut).replace(".", "p").replace("-", "m")
                safe_tol = str(sperp_tol or "default").replace(".", "p").replace("-", "m")
                (log_dir / f"bpw_ecut_{safe_ecut}_sperp_{safe_tol}.log").write_text(stdout)
            has_diag = SELECT_RE.search(stdout) and "BPW shell diag:" in stdout
            if result.returncode != 0 and not (args.allow_diagnostic_stop and has_diag):
                sys.stderr.write(stdout)
                raise SystemExit(
                    f"SALMON failed for Ecut={ecut}, Sperp_tol={sperp_tol or 'default'} "
                    f"with exit code {result.returncode}"
                )
            if result.returncode != 0:
                sys.stderr.write(
                    f"warning: accepted diagnostic output for Ecut={ecut}, "
                    f"Sperp_tol={sperp_tol or 'default'} after SALMON exit code "
                    f"{result.returncode}\n"
                )
            rows.append(parse_summary(stdout, str(ecut), str(sperp_tol), args.case_label))

    fields = [
        "case_label",
        "ecut",
        "sperp_tol",
        "mode",
        "selected_shells",
        "selected_nraw",
        "selected_g2_max",
        "selected_emax",
        "shell_id",
        "g2",
        "e_shell",
        "nraw",
        "np",
        "fwp_avg_occ",
        "fwp_x",
        "fwp_y",
        "fwp_z",
        "c_x",
        "c_y",
        "c_z",
        "c_comm_sum",
        "fsum_x",
        "fsum_y",
        "fsum_z",
        "fsum_avg",
        "min_sperp",
        "max_sperp",
        "herm",
    ]
    with output.open("w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
