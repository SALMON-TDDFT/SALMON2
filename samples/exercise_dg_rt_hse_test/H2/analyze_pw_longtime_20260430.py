from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import math
import re

BASE = Path(__file__).resolve().parent

@dataclass
class Case:
    name: str
    nt: int
    kick: bool
    sysname: str
    log: str

CASES = [
    Case("nt200_nokick", 200, False, "H2_contract_nt200_t0", "run_pw_nt200_nokick.log"),
    Case("nt200_kick",   200, True,  "H2_contract_nt200_kick", "run_pw_nt200_kick.log"),
    Case("nt400_nokick", 400, False, "H2_contract_nt400_t0", "run_pw_nt400_nokick.log"),
    Case("nt400_kick",   400, True,  "H2_contract_nt400_kick", "run_pw_nt400_kick.log"),
]

F2 = lambda x: float(x.replace("D", "E").replace("d", "e"))


def read_table(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    if not path.exists():
        return rows
    for line in path.read_text(errors="ignore").splitlines():
        s = line.strip()
        if (not s) or s.startswith("#"):
            continue
        try:
            rows.append([F2(tok) for tok in s.split()])
        except Exception:
            pass
    return rows


def tail_mean(rows: list[list[float]], idx: int, n: int = 20) -> float:
    vals = [r[idx] for r in rows if len(r) > idx]
    if not vals:
        return math.nan
    t = vals[-min(n, len(vals)):]
    return sum(t) / len(t)


def max_abs(rows: list[list[float]], idx: int) -> float:
    vals = [abs(r[idx]) for r in rows if len(r) > idx]
    return max(vals) if vals else math.nan


def final_val(rows: list[list[float]], idx: int) -> float:
    vals = [r[idx] for r in rows if len(r) > idx]
    return vals[-1] if vals else math.nan


def log_info(path: Path) -> dict:
    text = path.read_text(errors="ignore") if path.exists() else ""
    def has(pat: str) -> bool:
        return re.search(pat, text) is not None
    floor_lines = re.findall(r"\[RHO-FLOOR\].*", text)
    has_nan_error = (
        has(r"STOP\s+NaN")
        or has(r"NaN in ")
        or has(r"\[FP-DOMAIN\].*NaN")
        or has(r"\brt abort\b")
    )
    return {
        "end_salmon": has(r"(^|\n)\s*end SALMON\s*(\n|$)"),
        "has_nan": has_nan_error,
        "has_npw0_fatal": has(r"DG-FATAL unsupported n_pw==0"),
        "pw_on_seen": has(r"use_plane_wave_basis\s*=\s*T"),
        "npw_seen": has(r"n_plane_waves\s*=\s*[1-9]"),
        "mixed_ready_seen": has(r"mixed_basis_ready\s*=\s*T"),
        "rho_floor_count": len(floor_lines),
        "rho_floor_last": floor_lines[-1] if floor_lines else "",
    }


def summarize_case(case: Case) -> dict:
    cur = read_table(BASE / f"{case.sysname}_dg_current_decomp.data")
    eng = read_table(BASE / f"{case.sysname}_rt_energy.data")
    info = log_info(BASE / case.log)

    return {
        "case": case.name,
        "nt": case.nt,
        "kick": case.kick,
        "end_SALMON": info["end_salmon"],
        "NaN_or_STOP": info["has_nan"],
        "n_pw0_fatal": info["has_npw0_fatal"],
        "pw_on_seen": info["pw_on_seen"],
        "n_plane_waves_seen": info["npw_seen"],
        "mixed_basis_ready_seen": info["mixed_ready_seen"],
        "rho_floor_count": info["rho_floor_count"],
        "rho_floor_last": info["rho_floor_last"],
        "last_it": final_val(cur, 0),
        "max_abs_J_para_x": max_abs(cur, 2),
        "max_abs_J_dia_x": max_abs(cur, 5),
        "max_abs_J_total_x": max_abs(cur, 8),
        "tail_mean_J_para_x": tail_mean(cur, 2),
        "tail_mean_J_total_x": tail_mean(cur, 8),
        "final_J_para_x": final_val(cur, 2),
        "final_J_total_x": final_val(cur, 8),
        "final_rho_drift": final_val(cur, 11),
        "final_coef_occ_norm_drift": final_val(cur, 12),
        "final_cdagsc_i": final_val(cur, 14),
        "final_occvirt_leakage": final_val(cur, 16),
        "max_abs_rho_drift": max_abs(cur, 11),
        "max_abs_coef_occ_norm_drift": max_abs(cur, 12),
        "max_abs_cdagsc_i": max_abs(cur, 14),
        "max_abs_occvirt_leakage": max_abs(cur, 16),
        "final_Ne_raw": final_val(eng, 2),
        "final_Ne_raw_minus_Ne0": final_val(eng, 3),
        "max_abs_Ne_raw_minus_Ne0": max_abs(eng, 3),
    }


def fmt(v: float) -> str:
    if isinstance(v, bool):
        return "T" if v else "F"
    if isinstance(v, str):
        return v
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "nan"
    return f"{v:.6e}"


def main() -> None:
    rows = [summarize_case(c) for c in CASES]

    csv_path = BASE / "pw_longtime_20260430_summary.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    by = {r["case"]: r for r in rows}
    delta_rows = []
    for nt in (200, 400):
        a = by[f"nt{nt}_nokick"]
        b = by[f"nt{nt}_kick"]
        delta_rows.append({
            "nt": nt,
            "d_final_J_total_x": b["final_J_total_x"] - a["final_J_total_x"],
            "d_final_J_para_x": b["final_J_para_x"] - a["final_J_para_x"],
            "d_tail_mean_J_total_x": b["tail_mean_J_total_x"] - a["tail_mean_J_total_x"],
            "d_tail_mean_J_para_x": b["tail_mean_J_para_x"] - a["tail_mean_J_para_x"],
            "d_final_Ne_raw": b["final_Ne_raw"] - a["final_Ne_raw"],
            "d_final_rho_drift": b["final_rho_drift"] - a["final_rho_drift"],
            "d_final_coef_occ_norm_drift": b["final_coef_occ_norm_drift"] - a["final_coef_occ_norm_drift"],
            "d_final_cdagsc_i": b["final_cdagsc_i"] - a["final_cdagsc_i"],
            "d_final_occvirt_leakage": b["final_occvirt_leakage"] - a["final_occvirt_leakage"],
        })

    delta_csv = BASE / "pw_longtime_20260430_delta.csv"
    with delta_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(delta_rows[0].keys()))
        w.writeheader()
        w.writerows(delta_rows)

    md_path = BASE / "pw_longtime_20260430_summary.md"
    with md_path.open("w") as f:
        f.write("|Case|nt|kick|end SALMON|NaN|Ne_raw(final)|max|J_total_x||max|J_para_x||drift metrics (final)|\n")
        f.write("|---|---:|:---:|:---:|:---:|---:|---:|---:|---|\n")
        for r in rows:
            drift = (
                f"rho={fmt(r['final_rho_drift'])}, "
                f"coef={fmt(r['final_coef_occ_norm_drift'])}, "
                f"CdagSC-I={fmt(r['final_cdagsc_i'])}, "
                f"occvirt={fmt(r['final_occvirt_leakage'])}"
            )
            f.write(
                f"|{r['case']}|{r['nt']}|{'kick' if r['kick'] else 'no-kick'}|"
                f"{'T' if r['end_SALMON'] else 'F'}|{'T' if r['NaN_or_STOP'] else 'F'}|"
                f"{fmt(r['final_Ne_raw'])}|{fmt(r['max_abs_J_total_x'])}|{fmt(r['max_abs_J_para_x'])}|{drift}|\n"
            )

        f.write("\n")
        f.write("|nt|ΔJ_total(final)|ΔJ_para(final)|tail ΔJ_total|tail ΔJ_para|ΔNe_raw|Δrho_drift|Δcoef_occ_norm_drift|ΔCdagSC-I|Δoccvirt_leakage|\n")
        f.write("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for d in delta_rows:
            f.write(
                f"|{d['nt']}|{fmt(d['d_final_J_total_x'])}|{fmt(d['d_final_J_para_x'])}|"
                f"{fmt(d['d_tail_mean_J_total_x'])}|{fmt(d['d_tail_mean_J_para_x'])}|{fmt(d['d_final_Ne_raw'])}|"
                f"{fmt(d['d_final_rho_drift'])}|{fmt(d['d_final_coef_occ_norm_drift'])}|{fmt(d['d_final_cdagsc_i'])}|"
                f"{fmt(d['d_final_occvirt_leakage'])}|\n"
            )

    print(csv_path)
    print(delta_csv)
    print(md_path)


if __name__ == "__main__":
    main()
