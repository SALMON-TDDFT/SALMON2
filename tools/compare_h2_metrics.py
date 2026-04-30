#!/usr/bin/env python3
import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple


@dataclass
class CaseResult:
    label: str
    ne0: Optional[float]
    ne_last: Optional[float]
    ne_drift: Optional[float]
    n_energy_rows: int
    n_rt_rows: int
    max_abs_rt_col2_4: Optional[float]
    stop_itt_nan: Optional[int]
    max_raw_coef: Optional[float]
    max_dcoef_dt_h0: Optional[float]


FLOAT_RE = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"


def parse_numeric_rows(path: Path) -> List[List[float]]:
    rows: List[List[float]] = []
    if not path.exists():
        return rows
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            try:
                row = [float(x) for x in s.split()]
            except ValueError:
                continue
            if row:
                rows.append(row)
    return rows


def parse_stop_itt_from_log(path: Path) -> Optional[int]:
    if not path.exists():
        return None
    pat = re.compile(r"\[NaN\]\s+coef:.*?itt=\s*(\d+)")
    best: Optional[int] = None
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            m = pat.search(line)
            if not m:
                continue
            v = int(m.group(1))
            best = v if best is None else min(best, v)
    return best


def parse_rhs_metrics(path: Path) -> Tuple[Optional[float], Optional[float]]:
    if not path.exists():
        return None, None
    pat_raw = re.compile(rf"stage=raw_coef.*?max_abs=\s*({FLOAT_RE})")
    pat_h0 = re.compile(rf"stage=dcoef_dt_h0.*?max_abs=\s*({FLOAT_RE})")
    raw_max: Optional[float] = None
    h0_max: Optional[float] = None
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            m1 = pat_raw.search(line)
            if m1:
                v = float(m1.group(1))
                raw_max = v if raw_max is None else max(raw_max, v)
            m2 = pat_h0.search(line)
            if m2:
                v = float(m2.group(1))
                h0_max = v if h0_max is None else max(h0_max, v)
    return raw_max, h0_max


def format_float(x: Optional[float]) -> str:
    if x is None:
        return "-"
    if x == 0.0:
        return "0"
    ax = abs(x)
    if ax >= 1e4 or ax < 1e-3:
        return f"{x:.3e}"
    return f"{x:.8f}"


def analyze_case(label: str, energy: Path, rt: Path, log: Path) -> CaseResult:
    erows = parse_numeric_rows(energy)
    rrows = parse_numeric_rows(rt)

    ne0 = erows[0][3] if erows and len(erows[0]) >= 4 else None
    ne_last = erows[-1][3] if erows and len(erows[-1]) >= 4 else None
    ne_drift = (ne_last - ne0) if (ne0 is not None and ne_last is not None) else None

    max_abs_rt_col2_4: Optional[float] = None
    for row in rrows:
        if len(row) >= 4:
            for v in row[1:4]:
                av = abs(v)
                max_abs_rt_col2_4 = av if max_abs_rt_col2_4 is None else max(max_abs_rt_col2_4, av)

    stop_itt = parse_stop_itt_from_log(log)
    raw_max, h0_max = parse_rhs_metrics(log)

    return CaseResult(
        label=label,
        ne0=ne0,
        ne_last=ne_last,
        ne_drift=ne_drift,
        n_energy_rows=len(erows),
        n_rt_rows=len(rrows),
        max_abs_rt_col2_4=max_abs_rt_col2_4,
        stop_itt_nan=stop_itt,
        max_raw_coef=raw_max,
        max_dcoef_dt_h0=h0_max,
    )


def main() -> int:
    p = argparse.ArgumentParser(
        description=(
            "Compare H2 run metrics across cases. "
            "Use repeated --case LABEL ENERGY RT LOG."
        )
    )
    p.add_argument(
        "--case",
        nargs=4,
        metavar=("LABEL", "ENERGY", "RT", "LOG"),
        action="append",
        required=True,
    )
    args = p.parse_args()

    results: List[CaseResult] = []
    for label, energy_s, rt_s, log_s in args.case:
        results.append(analyze_case(label, Path(energy_s), Path(rt_s), Path(log_s)))

    header = [
        "label",
        "Ne0",
        "Ne_last",
        "dNe",
        "energy_rows",
        "rt_rows",
        "max|rt[2:4]|",
        "NaN_itt",
        "max_raw_coef",
        "max_dcoef_dt_h0",
    ]
    print("\t".join(header))
    for r in results:
        print(
            "\t".join(
                [
                    r.label,
                    format_float(r.ne0),
                    format_float(r.ne_last),
                    format_float(r.ne_drift),
                    str(r.n_energy_rows),
                    str(r.n_rt_rows),
                    format_float(r.max_abs_rt_col2_4),
                    "-" if r.stop_itt_nan is None else str(r.stop_itt_nan),
                    format_float(r.max_raw_coef),
                    format_float(r.max_dcoef_dt_h0),
                ]
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())