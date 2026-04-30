#!/usr/bin/env bash
set -u

ROOT='/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2'
SAMP="$ROOT/samples/exercise_dg_rt_hse_test/H2"
EXE="$ROOT/build/salmon"
OUT="$SAMP/pw_ab_sweep_20260429_nt40"
mkdir -p "$OUT/inputs" "$OUT/logs" "$OUT/data"

mk_input() {
  local src="$1" dst="$2" sys="$3" npw="$4" kick="$5"
  cp "$src" "$dst"
    perl -0777 -i -pe "s/sysname\s*=\s*'[^']*'/sysname = '$sys'/g; s/n_plane_waves_dg\s*=\s*[0-9]+/n_plane_waves_dg = $npw/g; s/I_wcm2_1\s*=\s*[-0-9.dD+]+/I_wcm2_1 = $kick/g; s/nproc_rgrid\(1:3\)\s*=\s*[^\n]*/nproc_rgrid(1:3) = 1,1,1/g" "$dst"
}

for npw in 10 20 32; do
    mk_input "$SAMP/input_contract_h2_nt40_mixed_npw1_nokick" "$OUT/inputs/input_nokick_pw${npw}.in" "H2_nt40_nokick_pw${npw}" "$npw" "0.0d0"
  mk_input "$SAMP/input_contract_h2_nt40_mixed_npw1_kick" "$OUT/inputs/input_kick_pw${npw}.in" "H2_nt40_kick_pw${npw}" "$npw" "1.0d12"
done

run_case() {
  local tag="$1" in_rel="$2" sfp="$3" mfp="$4" sys="$5"
  local log="$OUT/logs/${tag}.log"
    echo "[RUN] $tag"

  pushd "$SAMP" >/dev/null
  rm -f "${sys}"_* "$log"
  env SALMON_DG_SFP_MODE="$sfp" SALMON_DG_MFP_MODE="$mfp" SALMON_DG_HPP_MODE=d SALMON_DG_REORTH_MIXED_OCC=0 OMP_NUM_THREADS=1 \
      mpiexec -n 2 "$EXE" < "$in_rel" > "$log" 2>&1
  local rc=$?
  if [[ $rc -eq 0 ]]; then
    cp "${sys}_dg_current_decomp.data" "$OUT/data/${tag}_dg_current_decomp.data"
    cp "${sys}_rt_energy.data" "$OUT/data/${tag}_rt_energy.data"
    cp "${sys}_rt.data" "$OUT/data/${tag}_rt.data"
  fi
  popd >/dev/null
  echo "$tag,$rc" >> "$OUT/run_status.csv"
    echo "[DONE] $tag rc=$rc"
  return 0
}

: > "$OUT/run_status.csv"
for npw in 10 20 32; do
    run_case "A_nokick_pw${npw}" "pw_ab_sweep_20260429_nt40/inputs/input_nokick_pw${npw}.in" l 0 "H2_nt40_nokick_pw${npw}"
    run_case "B_nokick_pw${npw}" "pw_ab_sweep_20260429_nt40/inputs/input_nokick_pw${npw}.in" f f "H2_nt40_nokick_pw${npw}"
    run_case "A_kick_pw${npw}"   "pw_ab_sweep_20260429_nt40/inputs/input_kick_pw${npw}.in"   l 0 "H2_nt40_kick_pw${npw}"
    run_case "B_kick_pw${npw}"   "pw_ab_sweep_20260429_nt40/inputs/input_kick_pw${npw}.in"   f f "H2_nt40_kick_pw${npw}"
done

python3 - <<'PY'
from pathlib import Path
import re, csv, math
out = Path('/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/samples/exercise_dg_rt_hse_test/H2/pw_ab_sweep_20260429_nt40')
logd = out/'logs'
datad = out/'data'

mixed_re = re.compile(r"Mixed-S composition: pw_eff=\s*([0-9.EeDd+\-]+).*pw_dom_modes=([0-9]+)\s*/\s*([0-9]+)")


def f2(x):
    return float(x.replace('D','E').replace('d','e'))

def read_table(p):
    rows=[]
    if not p.exists():
        return rows
    for ln in p.read_text().splitlines():
        s=ln.strip()
        if not s or s.startswith('#'):
            continue
        try:
            rows.append([f2(t) for t in s.split()])
        except Exception:
            pass
    return rows

def max_abs_col(rows, idx):
    vals=[abs(r[idx]) for r in rows if len(r)>idx]
    return max(vals) if vals else math.nan

def max_abs_diff_col(a,b,idx):
    n=min(len(a),len(b))
    vals=[abs(a[i][idx]-b[i][idx]) for i in range(n) if len(a[i])>idx and len(b[i])>idx]
    return max(vals) if vals else math.nan

def scan_mixed(logp):
    vals=[]
    dom=[]
    if not logp.exists():
        return math.nan, '', 0
    for ln in logp.read_text(errors='ignore').splitlines():
        m=mixed_re.search(ln)
        if m:
            vals.append(f2(m.group(1)))
            dom.append(f"{int(m.group(2))}/{int(m.group(3))}")
    if not vals:
        return math.nan, '', 0
    return max(vals), dom[-1], len(vals)

rows=[]
for npw in (10,20,32):
    for kind in ('nokick','kick'):
        A = f'A_{kind}_pw{npw}'
        B = f'B_{kind}_pw{npw}'

        A_cur = read_table(datad/f'{A}_dg_current_decomp.data')
        B_cur = read_table(datad/f'{B}_dg_current_decomp.data')
        A_eng = read_table(datad/f'{A}_rt_energy.data')
        B_eng = read_table(datad/f'{B}_rt_energy.data')

        A_pweff, A_dom, A_hits = scan_mixed(logd/f'{A}.log')
        B_pweff, B_dom, B_hits = scan_mixed(logd/f'{B}.log')

        rec = {
            'npw': npw,
            'kind': kind,
            'A_pw_eff_max': A_pweff,
            'B_pw_eff_max': B_pweff,
            'd_pw_eff_max': abs(A_pweff-B_pweff) if math.isfinite(A_pweff) and math.isfinite(B_pweff) else math.nan,
            'A_pw_dom_last': A_dom,
            'B_pw_dom_last': B_dom,
            'A_pw_lines': A_hits,
            'B_pw_lines': B_hits,
            'A_Jpara_x_maxabs': max_abs_col(A_cur,2),
            'B_Jpara_x_maxabs': max_abs_col(B_cur,2),
            'd_Jpara_x_maxabs': max_abs_diff_col(A_cur,B_cur,2),
            'A_Jdia_x_maxabs': max_abs_col(A_cur,5),
            'B_Jdia_x_maxabs': max_abs_col(B_cur,5),
            'd_Jdia_x_maxabs': max_abs_diff_col(A_cur,B_cur,5),
            'A_Jtot_x_maxabs': max_abs_col(A_cur,8),
            'B_Jtot_x_maxabs': max_abs_col(B_cur,8),
            'd_Jtot_x_maxabs': max_abs_diff_col(A_cur,B_cur,8),
            'A_Edrift_maxabs': max_abs_col(A_eng,2),
            'B_Edrift_maxabs': max_abs_col(B_eng,2),
            'd_Edrift_maxabs': max_abs_diff_col(A_eng,B_eng,2),
            'A_Nedrift_maxabs': max_abs_col(A_eng,4),
            'B_Nedrift_maxabs': max_abs_col(B_eng,4),
            'd_Nedrift_maxabs': max_abs_diff_col(A_eng,B_eng,4),
            'A_PWweight_maxabs': max_abs_col(A_eng,5),
            'B_PWweight_maxabs': max_abs_col(B_eng,5),
            'd_PWweight_maxabs': max_abs_diff_col(A_eng,B_eng,5),
        }
        rows.append(rec)

fields=list(rows[0].keys()) if rows else []
with (out/'summary.csv').open('w', newline='') as f:
    w=csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

with (out/'summary.md').open('w') as f:
    f.write('|npw|kind|d_pw_eff_max|A_dom|B_dom|d_Jtot_x|max|d_Edrift|max|d_Nedrift|max|\n')
    f.write('|---:|:---:|---:|:---|:---|---:|---:|---:|\n')
    def fmt(v):
        return f"{v:.6e}" if math.isfinite(v) else 'nan'
    for r in rows:
        f.write(f"|{r['npw']}|{r['kind']}|{fmt(r['d_pw_eff_max'])}|{r['A_pw_dom_last']}|{r['B_pw_dom_last']}|{fmt(r['d_Jtot_x_maxabs'])}|{fmt(r['d_Edrift_maxabs'])}|{fmt(r['d_Nedrift_maxabs'])}|\n")
PY
