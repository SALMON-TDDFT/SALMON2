#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
stress_f90="$repo_root/src/common/stress.f90"
input_pp_f90="$repo_root/src/atom/pp/input_pp.f90"

check_block_contains() {
  file_path=$1
  block_start=$2
  block_end=$3
  pattern=$4
  if ! sed -n "/$block_start/,/$block_end/p" "$file_path" | rg -q "$pattern"; then
    echo "missing pattern in $file_path: $pattern" >&2
    return 1
  fi
}

check_block_absent() {
  file_path=$1
  block_start=$2
  block_end=$3
  pattern=$4
  if sed -n "/$block_start/,/$block_end/p" "$file_path" | rg -q "$pattern"; then
    echo "unexpected legacy pattern in $file_path: $pattern" >&2
    return 1
  fi
}

rg -q "use prep_pp_sub, only: eval_local_sr_shared_u" "$stress_f90"
rg -q "call build_local_sr_shared_u\\(pp, ik\\)" "$input_pp_f90"

check_block_contains "$stress_f90" "subroutine calc_stress_loc_sr_rs" "end subroutine calc_stress_loc_sr_rs" "call eval_local_sr_shared_u\\(pp,ik,r_abs,u_r,du_r,intr\\)"
check_block_absent "$stress_f90" "subroutine calc_stress_loc_sr_rs" "end subroutine calc_stress_loc_sr_rs" "pp%dvloctbl"
check_block_absent "$stress_f90" "subroutine calc_stress_loc_sr_rs" "end subroutine calc_stress_loc_sr_rs" "pp%vloctbl"

check_block_contains "$stress_f90" "function eval_local_sr_vg_from_table" "end function eval_local_sr_vg_from_table" "call eval_local_sr_shared_u\\(pp,ik,r1,u_r,du_r,i-1\\)"
check_block_absent "$stress_f90" "function eval_local_sr_vg_from_table" "end function eval_local_sr_vg_from_table" "pp%vloctbl"

echo "shared-u localSR source checks passed"
