#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
script="$repo_root/tests/tools/run_primitive_si_r2scan_stress_smoke.sh"

test -x "$script"
grep -Fq "al_vec1 =" "$script"
grep -Fq "al_vec2 =" "$script"
grep -Fq "al_vec3 =" "$script"
grep -Fq '"Si_primitive_r2scan_gs" "r2scan" "no"' "$script"
grep -Fq '"Si_primitive_r2scan_stress" "r2scan" "yes"' "$script"
grep -Fq -- '--stress-output-level LEVEL' "$script"
grep -Fq 'STRESS_OUTPUT_LEVEL="high"' "$script"
grep -Fq 'case "$STRESS_OUTPUT_LEVEL" in' "$script"
grep -Fq 'low|middle|high) ;;' "$script"
grep -Fq "stress_output_level = '\${STRESS_OUTPUT_LEVEL}'" "$script"
grep -Fq "nproc_rgrid(1) = 1" "$script"
grep -Fq "nproc_rgrid(2) = 1" "$script"
grep -Fq "nproc_rgrid(3) = 1" "$script"
