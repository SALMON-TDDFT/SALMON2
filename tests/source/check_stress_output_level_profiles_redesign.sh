#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
write_f90="$repo_root/src/io/write.f90"

assert_absent() {
  pattern="$1"
  if rg -F -q "$pattern" "$write_f90"; then
    echo "unexpected pattern still present: $pattern" >&2
    exit 1
  fi
}

grep -Fq "terms_on = (yn_out_stress_terms == 'y')" "$write_f90"
grep -Fq "details_on = (yn_out_stress_details == 'y')" "$write_f90"
grep -Fq "numerics_on = (yn_out_stress_numerics == 'y')" "$write_f90"

grep -Fq 'Stress tensor [GPa]' "$write_f90"
grep -Fq 'Stress residual diagnostics [Hartree]' "$write_f90"
grep -Fq 'Stress decomposition tensor [GPa]' "$write_f90"
grep -Fq 'Stress decomposition detail [GPa]' "$write_f90"

grep -Fq "call write_stress_tensor_row_gpa(fh, 'Kinetic'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'Total'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'XC-local'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'XC-grad'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'XC-tau'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'Local-SR'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'Local-LR'" "$write_f90"
grep -Fq 'call write_nl_l_channel_tensor_rows_gpa(' "$write_f90"
grep -Fq 'call write_nl_species_l_channel_tensor_rows_gpa(' "$write_f90"

assert_absent 'stress_detail_level' "$write_f90"
assert_absent 'Stress decomposition pressure [GPa]'
