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

grep -Fq "call begin_compact_column_header(fh)" "$write_f90"
grep -Fq "call write_compact_column_header(fh, line_cols, col, 'iter')" "$write_f90"
grep -Fq "call write_compact_column_header(fh, line_cols, col, 'E_total [eV]')" "$write_f90"
grep -Fq "call write_compact_column_header(fh, line_cols, col, 'P_total [GPa]')" "$write_f90"
grep -Fq "call write_compact_column_header(fh, line_cols, col, 'P_xc_grad_payload [GPa]')" "$write_f90"
grep -Fq "call write_compact_column_header(fh, line_cols, col, 'P_loc_sr_rs [GPa]')" "$write_f90"
grep -Fq "call end_compact_column_header(fh)" "$write_f90"
grep -Fq "call write_nl_l_pressure_column_headers" "$write_f90"
grep -Fq "call write_nl_species_l_pressure_column_headers" "$write_f90"
grep -Fq "call write_nl_l_pressure_values" "$write_f90"
grep -Fq "call write_nl_species_l_pressure_values" "$write_f90"

assert_absent "# stress_output_level = "
assert_absent "# yn_out_stress_terms = "
assert_absent "# yn_out_stress_details = "
assert_absent "# yn_out_stress_numerics = "
assert_absent "# stress_l_decomp = "
assert_absent "# col "
assert_absent "(middle+, else 0)"
assert_absent "(high, else 0)"
assert_absent "real(8) :: p(13)"
assert_absent "(i6,11e20.12,13f20.8)"
assert_absent "stress_detail_level"
