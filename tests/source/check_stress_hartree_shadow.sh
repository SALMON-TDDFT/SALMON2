#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
structures_f90="$repo_root/src/common/structures.f90"
stress_f90="$repo_root/src/common/stress.f90"
main_dft_f90="$repo_root/src/gs/main_dft.f90"
write_f90="$repo_root/src/io/write.f90"

grep -Fq "stress_har_shadow(3,3)" "$structures_f90"
grep -Fq "subroutine calc_stress_har_shadow(" "$stress_f90"
grep -Fq "call calc_stress_har_shadow(" "$main_dft_f90"
grep -Fq "system%stress_har_shadow = -strs_sum" "$stress_f90"
grep -Fq "Tr(har_shadow)*V + E_h" "$write_f90"
grep -Fq "P_har_shadow [GPa]" "$write_f90"
grep -Fq "P_har_shadow+loc_fullobj+ewald [GPa]" "$write_f90"
