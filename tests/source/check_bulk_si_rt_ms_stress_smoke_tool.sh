#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
script="$repo_root/tests/tools/run_bulk_si_rt_ms_stress_smoke.sh"
cmake_lists="$repo_root/tests/source/CMakeLists.txt"

grep -Fq 'check_bulk_si_rt_ms_stress_smoke_tool.sh' "$cmake_lists"
test -x "$script"
grep -Fq 'bulk Si RT/MS stress smoke helper' "$script"
grep -Fq 'data_for_restart' "$script"
grep -Fq 'restart/' "$script"
grep -Fq "theory = 'tddft_response'" "$script"
grep -Fq "theory = 'multi_scale_maxwell_tddft'" "$script"
grep -Fq "xc = 'PZ'" "$script"
grep -Fq "xc = 'r2scan'" "$script"
grep -Fq 'primitive' "$script"
grep -Fq 'conventional' "$script"
grep -Fq "stress_output_level = '" "$script"
grep -Fq "yn_out_stress = 'y'" "$script"
grep -Fq 'ln -s' "$script"
grep -Fq 'RT and MS stress smoke completed.' "$script"
