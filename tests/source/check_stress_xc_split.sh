#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
write_f90="$repo_root/src/io/write.f90"
stress_f90="$repo_root/src/common/stress.f90"
structures_f90="$repo_root/src/common/structures.f90"

grep -Fq 'real(8) :: stress_x(3,3)' "$structures_f90"
grep -Fq 'real(8) :: stress_c(3,3)' "$structures_f90"
grep -Fq 'real(8) :: stress_x_e_vx' "$structures_f90"
grep -Fq 'real(8) :: stress_c_e_vc' "$structures_f90"

grep -Fq 'system%stress_x = 0d0' "$stress_f90"
grep -Fq 'system%stress_c = 0d0' "$stress_f90"
grep -Fq 'system%stress_x_e_vx = 0d0' "$stress_f90"
grep -Fq 'system%stress_c_e_vc = 0d0' "$stress_f90"
grep -Fq 'system%stress_xc = system%stress_x + system%stress_c' "$stress_f90"
grep -Fq 'system%stress_xc_e_vxc = system%stress_x_e_vx + system%stress_c_e_vc' "$stress_f90"

grep -Fq "call write_stress_tensor_row_gpa(fh, '  XC'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, '    X'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, '    C'" "$write_f90"
