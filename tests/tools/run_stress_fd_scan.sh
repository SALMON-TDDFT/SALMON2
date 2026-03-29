#!/bin/bash
# run_stress_fd_scan.sh -- sector-resolved FD validation for SALMON stress
# MUST specify either --setup-submit or --fd-only
#
# Phase 1: run_stress_fd_scan.sh --setup-submit [options]
# Phase 2: run_stress_fd_scan.sh --fd-only --base-dir DIR --scan-dir NAME

set -euo pipefail

MODE=""
BASE_DIR=""
SCAN_DIR="stress_scan_fd_C"
INP="Si_gs.inp"
PSEUDO="Si_rps.dat"
A_CENTER="5.43"
N_POINTS="7"
DA="0.01"
BINARY=""

usage() {
  echo "Usage:"
  echo "  $0 --setup-submit --base-dir DIR --scan-dir NAME --inp FILE --pseudo FILE"
  echo "              --a-center FLOAT --n-points INT --da FLOAT --binary PATH"
  echo "  $0 --fd-only --base-dir DIR --scan-dir NAME"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --setup-submit) MODE="setup-submit"; shift ;;
    --fd-only) MODE="fd-only"; shift ;;
    --base-dir) BASE_DIR="$2"; shift 2 ;;
    --scan-dir) SCAN_DIR="$2"; shift 2 ;;
    --inp) INP="$2"; shift 2 ;;
    --pseudo) PSEUDO="$2"; shift 2 ;;
    --a-center) A_CENTER="$2"; shift 2 ;;
    --n-points) N_POINTS="$2"; shift 2 ;;
    --da) DA="$2"; shift 2 ;;
    --binary) BINARY="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

[[ -z "$MODE" ]] && { echo "ERROR: must specify --setup-submit or --fd-only"; usage; }
[[ -z "$BASE_DIR" ]] && { echo "ERROR: --base-dir required"; usage; }

SCAN_PATH="$BASE_DIR/$SCAN_DIR"

if [[ "$MODE" == "setup-submit" ]]; then
  [[ -z "$BINARY" ]] && { echo "ERROR: --binary required for --setup-submit"; usage; }
  [[ ! -f "$BASE_DIR/$INP" ]] && { echo "ERROR: $BASE_DIR/$INP not found"; exit 1; }
  [[ ! -f "$BASE_DIR/$PSEUDO" ]] && { echo "ERROR: $BASE_DIR/$PSEUDO not found"; exit 1; }
  [[ ! -f "$BASE_DIR/calc.sh" ]] && { echo "ERROR: $BASE_DIR/calc.sh not found"; exit 1; }

  mkdir -p "$SCAN_PATH"
  : > "$SCAN_PATH/scan_jobs.txt"
  echo "Scan directory: $SCAN_PATH"

  HALF=$(( N_POINTS / 2 ))
  declare -a A_LIST
  for (( i=-HALF; i<=HALF; i++ )); do
    A=$(python3 -c "a_center=float('$A_CENTER'); i_val=int('$i'); da=float('$DA'); print(f'{a_center + i_val * da:.2f}')")
    A_LIST+=("$A")
  done

  for A in "${A_LIST[@]}"; do
    TAG=$(echo "$A" | tr -d '.')
    DIR="$SCAN_PATH/a${TAG}"
    mkdir -p "$DIR"
    cp "$BASE_DIR/$INP" "$BASE_DIR/$PSEUDO" "$BASE_DIR/calc.sh" "$DIR/"

    SYSNAME="Si_a${TAG}"
    perl -0pi -e "
      s/al\(1:3\)\s*=\s*[0-9.]+d0,\s*[0-9.]+d0,\s*[0-9.]+d0/al(1:3) = ${A}d0, ${A}d0, ${A}d0/;
      s/sysname\s*=\s*'[^']*'/sysname = '${SYSNAME}'/;
    " "$DIR/$INP"
    perl -0pi -e "s|EXEC=.*|EXEC=$BINARY|" "$DIR/calc.sh"

    (
      cd "$DIR"
      JOB_ID=$(pjsub calc.sh | grep -oP '(?<=Job )\d+')
      echo "$A $TAG $JOB_ID" >> "$SCAN_PATH/scan_jobs.txt"
      echo "  Submitted a=$A (tag=$TAG) job=$JOB_ID"
    )
  done

  echo
  echo "All jobs submitted. Job IDs saved to: $SCAN_PATH/scan_jobs.txt"
  echo "Monitor with: pjstat | grep -f <(awk '{print \$3}' $SCAN_PATH/scan_jobs.txt)"
  echo
  echo "When all jobs complete, run:"
  echo "  $0 --fd-only --base-dir $BASE_DIR --scan-dir $SCAN_DIR"
elif [[ "$MODE" == "fd-only" ]]; then
  [[ ! -f "$SCAN_PATH/scan_jobs.txt" ]] && { echo "ERROR: $SCAN_PATH/scan_jobs.txt not found. Run --setup-submit first."; exit 1; }

  mapfile -t ROWS < <(sort -k1 -n "$SCAN_PATH/scan_jobs.txt")
  N="${#ROWS[@]}"

  echo "Collecting *_stress_energy.data from $N scan points..."

  TMP=$(mktemp)
  for ROW in "${ROWS[@]}"; do
    A=$(echo "$ROW" | awk '{print $1}')
    TAG=$(echo "$ROW" | awk '{print $2}')
    DATA_FILE="$SCAN_PATH/a${TAG}/Si_a${TAG}_stress_energy.data"
    if [[ ! -f "$DATA_FILE" ]]; then
      echo "  WARNING: $DATA_FILE not found -- job may not have completed"
      continue
    fi
    DATA_ROW=$(grep -v '^#' "$DATA_FILE" | head -1)
    echo "$A $DATA_ROW" >> "$TMP"
  done

  echo "Running FD calculation..."

  python3 - "$TMP" "$SCAN_PATH/stress_fd_results.txt" <<'PY'
import sys

tmp_file = sys.argv[1]
out_file = sys.argv[2]

rows = []
for line in open(tmp_file):
    parts = line.split()
    a_ang = float(parts[0])
    cols = [float(x) for x in parts[1:]]
    rows.append((a_ang, cols))

rows.sort(key=lambda r: r[0])
n = len(rows)
if n < 3:
    print(f"ERROR: need at least 3 scan points, got {n}")
    sys.exit(1)

mid = n // 2
a_center_ang = rows[mid][0]

bohr_per_ang = 1.0 / 0.529177249
a_bohr = [r[0] * bohr_per_ang for r in rows]
h = a_bohr[1] - a_bohr[0]

coeff7 = [1/60, -9/60, 45/60, 0, -45/60, 9/60, -1/60]

def fd7(vals, step):
    if len(vals) == 7:
        return (coeff7[0]*vals[0] + coeff7[1]*vals[1] + coeff7[2]*vals[2]
              + coeff7[4]*vals[4] + coeff7[5]*vals[5] + coeff7[6]*vals[6]) / step
    return (vals[2] - vals[0]) / (2 * step)

def p_from_deda(dE_da_ha, a_bohr_center):
    dE_dV = dE_da_ha / (3 * a_bohr_center**2)
    return -dE_dV * 29421.02648

ha_per_ev = 1.0 / 27.211396132

sector_defs = [
    ('kin',    3,  13, 'Kinetic'),
    ('har',    4,  14, 'Hartree'),
    ('xc',     5,  15, 'XC'),
    ('loc_sr', 9,  None, 'Local SR'),
    ('loc_lr', 10, None, 'Local LR'),
    ('nl',     7,  24, 'Nonlocal'),
    ('ewald',  8,  None, 'Ewald (ion_ion)'),
    ('total',  2,  25, 'Total'),
]

lines_out = []
lines_out.append("=== Stress FD Validation ===")
lines_out.append(f"a_center = {a_center_ang:.2f} A")
lines_out.append(f"n_points = {n}")
lines_out.append("")
lines_out.append(f"{'sector':<15} {'P_analytic[GPa]':>17} {'P_FD[GPa]':>12} {'diff[GPa]':>11}  status")
lines_out.append("-" * 65)

for sec, ecol, pcol, label in sector_defs:
    if ecol == 2:
        e_vals = [r[1][ecol - 1] * ha_per_ev for r in rows]
    else:
        e_vals = [r[1][ecol - 1] for r in rows]

    dE_da = fd7(e_vals, h) if n == 7 else (e_vals[mid + 1] - e_vals[mid - 1]) / (2 * h)
    p_fd = p_from_deda(dE_da, a_bohr[mid])

    if sec == 'loc_sr':
        p_an = rows[mid][1][15] + rows[mid][1][16]
    elif sec == 'loc_lr':
        p_an = rows[mid][1][17] + rows[mid][1][18]
    elif sec == 'ewald':
        p_an = rows[mid][1][19] + rows[mid][1][20] + rows[mid][1][21] + rows[mid][1][22]
    elif pcol is not None:
        p_an = rows[mid][1][pcol - 1]
    else:
        p_an = None

    if p_an is not None:
        diff = p_an - p_fd
        ok = abs(diff) < 0.5
        status = "PASS" if ok else "FAIL"
        lines_out.append(f"{label:<15} {p_an:>17.3f} {p_fd:>12.3f} {diff:>11.3f}  {status}")
    else:
        lines_out.append(f"{label:<15} {'(no analytic)':>17} {p_fd:>12.3f} {'':>11}  (FD only)")

lines_out.append("")

output = "\n".join(lines_out)
print(output)
open(out_file, 'w').write(output + "\n")
print(f"Results saved to: {out_file}")
PY

  rm -f "$TMP"
fi
