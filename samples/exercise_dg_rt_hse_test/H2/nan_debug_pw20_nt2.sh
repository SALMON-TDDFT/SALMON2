#!/bin/bash
# NaN debug for B_kick_pw20 vs A_kick_pw20 at nt=2
# Purpose: Isolate where NaN first appears in coefficient space

set -e

SALMON_DIR='/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2'
EXE="${SALMON_DIR}/build/salmon"
TEMPLATE_DIR="$(cd "$(dirname "$BASH_SOURCE")"; pwd)"
TEMPLATE_BASE="${TEMPLATE_DIR}/input_contract_h2_nt40_mixed_npw1_kick"
WORK_DIR="$(pwd)/nan_debug_pw20_nt2"

mkdir -p "$WORK_DIR"/{logs,data}
cd "$WORK_DIR"

# Copy required pseudopotential
cp "${TEMPLATE_DIR}/H_rps.dat" . 2>/dev/null || echo "Warning: H_rps.dat not found, assuming it exists locally"

# Test matrix: 2 cases (A=legacy/0, B=full/f) × 2 modes (nokick, kick)
# Here we focus on kick cases with pw=20
declare -a CASES=(
  "A_kick_pw20:l:0"     # SFP_MODE=l MFP_MODE=0
  "B_kick_pw20:f:f"     # SFP_MODE=f MFP_MODE=f
)

mk_input_nt2() {
  local work_dir="$1" case_name="$2" sfp_mode="$3" mfp_mode="$4"
  local out_file="${work_dir}/${case_name}.in"
  
  # Copy base and modify for nt=2, npw=20, kick
  perl -0777 -pe "
    s/Nt\s*=\s*\d+/Nt = 2/;
    s/n_plane_waves_dg\s*=\s*\d+/n_plane_waves_dg = 20/;
    s/y_amplitude\s*=\s*[0-9.dD+-]+/y_amplitude = 1.0d12/;
    s/nproc_rgrid\(1:3\)\s*=\s*[^\/]+/nproc_rgrid(1:3) = 1,1,1/;
  " "$TEMPLATE_BASE" > "$out_file"
  
  echo "Generated: $out_file"
}

run_case_pair() {
  local work_dir="$1"
  
  for spec in "${CASES[@]}"; do
    IFS=':' read -r case_name sfp_mode mfp_mode <<< "$spec"
    
    in_file="${work_dir}/${case_name}.in"
    log="${work_dir}/logs/${case_name}.log"
    rc=0
    
    echo "[SETUP] $case_name (SFP_MODE=$sfp_mode MFP_MODE=$mfp_mode)"
    mk_input_nt2 "$work_dir" "$case_name" "$sfp_mode" "$mfp_mode"
    
    echo "[RUN] $case_name ..."
    env \
      SALMON_DG_SFP_MODE="$sfp_mode" \
      SALMON_DG_MFP_MODE="$mfp_mode" \
      SALMON_DG_HPP_MODE="d" \
      SALMON_DG_REORTH_MIXED_OCC="0" \
      OMP_NUM_THREADS=1 \
      mpiexec -n 2 "$EXE" < "$in_file" > "$log" 2>&1 || rc=$?
    
    echo "[DONE] $case_name (rc=$rc)"
    echo ""
  done
}

# Main
echo "=== NaN Debug Session: B_kick_pw20 vs A_kick_pw20 (nt=2) ==="
echo "Workspace: $WORK_DIR"
echo ""

run_case_pair "$WORK_DIR"

# Summary
echo ""
echo "=== Log Locations ==="
ls -lh "${WORK_DIR}/logs/"*.log

echo ""
echo "=== Error Grep ==="
for log in "${WORK_DIR}"/logs/*.log; do
  echo "File: $(basename "$log")"
  grep -E 'STOP|NaN|Error' "$log" | head -5 || echo "  (no error markers)"
done
