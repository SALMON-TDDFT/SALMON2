#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
stress_f90="$repo_root/src/common/stress.f90"
write_f90="$repo_root/src/io/write.f90"

grep -Fq "subroutine calc_stress_loc_sr_rs(" "$stress_f90"
grep -Fq "system%stress_loc_sr_rs = -strs_sum / V" "$stress_f90"
grep -Fq "P_loc_sr [GPa]" "$write_f90"
grep -Fq "P_loc_sr_rs [GPa]" "$write_f90"
grep -Fq "P_loc_sr_rs - P_loc_sr [GPa]" "$write_f90"

! rg -F -q "P_loc_sr_grad_G [GPa]" "$write_f90"
! rg -F -q "P_loc_sr_rs-G [GPa]" "$write_f90"
! rg -F -q "P_har+loc_sr_rs+lr_G+diag+ewa [GPa]" "$write_f90"
