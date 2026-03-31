#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
write_f90="$repo_root/src/io/write.f90"

grep -Fq 'Stress decomposition tensor [GPa]' "$write_f90"
grep -Fq "write(fh,'(1x,a8,1x,7a16)') 'sector', 'xx', 'yy', 'zz', 'xy', 'yz', 'xz', 'P'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'Kinetic'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'Hartree'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'XC'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'Local'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'Nonlocal'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'Ewald'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'Total'" "$write_f90"
