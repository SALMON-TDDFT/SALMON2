#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

build_dir="${SALMON_BUILD_DIR:-$repo_root/build-mpi-eigenexa-wannier-lib}"
case_dir="${SALMON_CASE_DIR:-$build_dir/samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac}"
salmon_bin="${SALMON_BIN:-$build_dir/salmon}"
input_field_off="${SALMON_INPUT_FIELD_OFF:-/tmp/inputfile_rt_w90_mixed_nofield_nt100}"
input_field_on="${SALMON_INPUT_FIELD_ON:-/tmp/inputfile_rt_w90_mixed_z_nt100}"
np="${SALMON_NP:-8}"
omp_threads="${OMP_NUM_THREADS:-2}"
out_dir="${SALMON_OBS_REGRESSION_OUT_DIR:-/tmp}"
tol="${SALMON_OBS_REGRESSION_TOL:-1.0e-10}"

off_log="$out_dir/salmon_mixed_z_obs_regression_field_off.log"
on_log="$out_dir/salmon_mixed_z_obs_regression_field_on.log"

require_file() {
  local path="$1"
  if [[ ! -e "$path" ]]; then
    echo "missing required path: $path" >&2
    exit 2
  fi
}

require_file "$salmon_bin"
require_file "$case_dir"
require_file "$input_field_off"
require_file "$input_field_on"

run_case() {
  local input="$1"
  local log="$2"
  (
    cd "$case_dir"
    env OMP_NUM_THREADS="$omp_threads" \
      SALMON_DG_MIXED_Z=1 \
      SALMON_DG_EXPDIAG_GLOBAL_FLUX=1 \
      SALMON_DG_MIXED_Z_USE_LOCAL_SPLIT_PROPAGATION=1 \
      SALMON_DG_MIXED_Z_LOCAL_RHO_WRITEBACK_PROD_TARGET=1 \
      SALMON_DG_MIXED_Z_LOCAL_PZ_WRITEBACK_TOTAL=1 \
      SALMON_DG_MIXED_Z_LOCAL_CURRENT_WRITEBACK_TOTAL=1 \
      mpirun -np "$np" "$salmon_bin" < "$input" > "$log" 2>&1
  )
}

run_case "$input_field_off" "$off_log"
run_case "$input_field_on" "$on_log"

python3 - "$tol" "$off_log" "$on_log" <<'PY'
import math
import re
import sys

tol = float(sys.argv[1])
cases = [("field-off", sys.argv[2]), ("field-on", sys.argv[3])]
checks = [
    (
        "rho",
        "DG-MIXEDZ-LOCAL-RHO-WRITEBACK-PROD-TARGET-CMP",
        [
            "int_diff_candidate_minus_before",
            "max_abs_candidate_minus_before",
            "rms_candidate_minus_before",
            "after_minus_candidate_int",
        ],
    ),
    (
        "pz",
        "DG-MIXEDZ-LOCAL-PZ-WRITEBACK-TOTAL-CMP",
        [
            "diff_candidate_minus_before",
            "after_minus_candidate",
            "after_minus_before",
        ],
    ),
    (
        "current",
        "DG-MIXEDZ-LOCAL-CURRENT-WRITEBACK-TOTAL-CMP",
        [
            "diff_candidate_minus_before_norm",
            "after_minus_candidate_norm",
            "after_minus_before_norm",
        ],
    ),
]
legacy_tags = [
    "DG-MIXEDZ-LOCAL-SPLIT",
    "DG-MIXEDZ-LOCAL-FRAG-",
    "DG-MIXEDZ-LOCAL-PZ-FIELDON",
    "DG-MIXEDZ-LOCAL-RHO-BLOCK-CMP",
    "DG-MIXEDZ-LOCAL-RHO-WW-GRID-CMP",
    "DG-MIXEDZ-LOCAL-CURRENT-PATH-CMP",
]

failed = False
for label, path in cases:
    data = open(path, "rb").read()
    text = data.decode(errors="ignore")
    print(f"=== {label} {path}")
    basic_bad = {
        "missing_end": "end SALMON" not in text,
        "FATAL": text.count("FATAL"),
        "NaN": text.count("NaN"),
        "SIGBUS": text.count("SIGBUS"),
        "NUL": data.count(b"\x00"),
    }
    print("basic", basic_bad)
    if any(v for v in basic_bad.values()):
        failed = True

    legacy_count = 0
    for tag in legacy_tags:
        count = text.count("[" + tag)
        legacy_count += count
        print(f"legacy {tag} count={count}")
    if legacy_count != 0:
        failed = True

    for name, tag, keys in checks:
        rows = []
        for line in text.splitlines():
            if f"[{tag}]" not in line:
                continue
            vals = {k: float(v) for k, v in re.findall(r"([A-Za-z0-9_]+)=\s*([-+0-9.Ee]+)", line)}
            for k in ["replacement_ready", "replacement_applied", "production_value_modified", "bad"]:
                m = re.search(k + r"=([FT])", line)
                vals[k] = m.group(1) if m else "NA"
            m = re.search(r"replacement_block_reason=\s*(\S+)", line)
            vals["reason"] = m.group(1) if m else "NA"
            rows.append(vals)

        bad = sum(row.get("bad") == "T" for row in rows)
        applied = sum(row.get("replacement_applied") == "T" for row in rows)
        modified = sum(row.get("production_value_modified") == "T" for row in rows)
        reasons = sorted({row.get("reason") for row in rows})
        print(f"{name} nline={len(rows)} bad={bad} applied={applied} modified={modified} reasons={reasons}")
        if not rows or bad or applied != len(rows) or modified != len(rows) or reasons != ["none"]:
            failed = True

        for key in keys:
            vals = [abs(float(row.get(key, math.nan))) for row in rows]
            max_abs = max(vals) if vals else math.nan
            print(f"  {key} max_abs={max_abs:.6e}")
            if not math.isfinite(max_abs) or max_abs > tol:
                failed = True

if failed:
    sys.exit(1)
PY
