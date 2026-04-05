#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
write_f90="$repo_root/src/io/write.f90"
rt_init_f90="$repo_root/src/rt/initialization_rt.f90"
rt_step_f90="$repo_root/src/rt/time_evolution_step.f90"
ms_main_f90="$repo_root/src/ms/main_ms.f90"

assert_absent() {
  pattern="$1"
  file="$2"
  if rg -F -q "$pattern" "$file"; then
    echo "unexpected pattern still present in $file: $pattern" >&2
    exit 1
  fi
}

grep -Fq "subroutine write_stress_rt(iter, ofl, dt, system, energy, pp)" "$write_f90"
grep -Fq "yn_out_stress_terms" "$write_f90"
grep -Fq "yn_out_stress_details" "$write_f90"
grep -Fq "yn_out_stress_numerics" "$write_f90"
grep -Fq "stress_l_decomp" "$write_f90"
grep -Fq 'write(ofl%fh_stress,10) "Stress", "Stress tensor time series"' "$write_f90"
grep -Fq 'write(ofl%fh_stress,10) "Pressure", "-Tr(stress)/3"' "$write_f90"
grep -Fq "call begin_compact_column_header(ofl%fh_stress)" "$write_f90"
grep -Fq "call write_compact_column_header(ofl%fh_stress, line_cols, col, 'time [" "$write_f90"
grep -Fq "call end_compact_column_header(ofl%fh_stress)" "$write_f90"
grep -Fq "write_nl_species_l_pressure_column_headers" "$write_f90"
grep -Fq "write_nl_species_l_pressure_values" "$write_f90"
grep -Fq "call write_stress_rt(-1, ofl, dt, system, energy, pp)" "$rt_init_f90"
grep -Fq "call write_stress_rt(stress_label_iter, ofl, dt, system, energy, pp)" "$rt_init_f90"
grep -Fq "call write_stress_rt(itt, ofl, dt, system, energy, pp)" "$rt_step_f90"
grep -Fq "call write_stress_rt(stress_label_iter, ofl, dt, system, energy, pp)" "$ms_main_f90"

assert_absent "yn_out_stress_decomp" "$write_f90"
assert_absent "# stress_output_level = " "$write_f90"
assert_absent "# yn_out_stress_terms = " "$write_f90"
assert_absent "# yn_out_stress_details = " "$write_f90"
assert_absent "# yn_out_stress_numerics = " "$write_f90"
assert_absent "# stress_l_decomp = " "$write_f90"
assert_absent "# col " "$write_f90"
