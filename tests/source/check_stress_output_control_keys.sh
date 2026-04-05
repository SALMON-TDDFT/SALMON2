#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
global_f90="$repo_root/src/io/salmon_global.f90"
input_f90="$repo_root/src/io/inputoutput.f90"
write_f90="$repo_root/src/io/write.f90"
stress_f90="$repo_root/src/common/stress.f90"
structures_f90="$repo_root/src/common/structures.f90"

assert_absent() {
  pattern="$1"
  file="$2"
  if rg -F -q "$pattern" "$file"; then
    echo "unexpected pattern still present in $file: $pattern" >&2
    exit 1
  fi
}

grep -Fq "stress_output_level" "$global_f90"
grep -Fq "yn_out_stress_terms" "$global_f90"
grep -Fq "yn_out_stress_details" "$global_f90"
grep -Fq "yn_out_stress_numerics" "$global_f90"
grep -Fq "stress_l_decomp" "$global_f90"

grep -Fq "stress_output_level" "$input_f90"
grep -Fq "yn_out_stress_terms" "$input_f90"
grep -Fq "yn_out_stress_details" "$input_f90"
grep -Fq "yn_out_stress_numerics" "$input_f90"
grep -Fq "stress_l_decomp" "$input_f90"
grep -Fq "call normalize_stress_output_controls" "$input_f90"
grep -Fq "stress_output_level must be 'low', 'middle', or 'high'" "$input_f90"
grep -Fq "stress_l_decomp must be 'no', 'species', or 'atom'" "$input_f90"
grep -Fq "stress_l_decomp='atom' is not implemented yet" "$input_f90"
grep -Fq "stress_l_decomp requires yn_out_stress_details='y'" "$input_f90"
grep -Fq "yn_out_stress_decomp is no longer supported" "$input_f90"
grep -Fq "stress_fd_detail is no longer supported" "$input_f90"

grep -Fq "stress_output_level = 'high'" "$input_f90"
grep -Fq "yn_out_stress_terms = ' '" "$input_f90"
grep -Fq "yn_out_stress_details = ' '" "$input_f90"
grep -Fq "yn_out_stress_numerics = ' '" "$input_f90"
grep -Fq "stress_l_decomp = '__unset__'" "$input_f90"
grep -Fq "stress_l_decomp = 'species'" "$input_f90"

grep -Fq "yn_out_stress_numerics" "$stress_f90"
grep -Fq "stress_l_decomp" "$stress_f90"
grep -Fq "stress_nl_species_l" "$stress_f90"
grep -Fq "stress_nl_species_l" "$structures_f90"

assert_absent "stress_fd_detail" "$write_f90"
assert_absent "yn_out_stress_decomp" "$write_f90"
assert_absent "stress_l_decomp must be 'no', 'total', 'species', or 'atom'" "$input_f90"
assert_absent "stress_l_decomp = 'total'" "$input_f90"
assert_absent "case('no','total','species','atom')" "$input_f90"
assert_absent "stress_nl_l" "$write_f90"
assert_absent "stress_nl_l" "$stress_f90"
assert_absent "stress_nl_l" "$structures_f90"
