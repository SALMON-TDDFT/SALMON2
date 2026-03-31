#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
write_f90="$repo_root/src/io/write.f90"

grep -Fq 'Stress decomposition tensor [GPa]' "$write_f90"
grep -Fq "write(fh,'(1x,a12,1x,7a16)') 'sector', 'xx', 'yy', 'zz', 'xy', 'yz', 'xz', 'P'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'Kinetic'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'Hartree'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'XC'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'Local'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'Nonlocal'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'Ewald'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'Total'" "$write_f90"
grep -Fq 'Stress decomposition pressure [GPa]' "$write_f90"
grep -Fq 'Stress decomposition detail [GPa]' "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, '  Local'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, '    Local-SR'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, '    Local-LR'" "$write_f90"
grep -Fq "call write_stress_tensor_row_gpa(fh, '  Nonlocal'" "$write_f90"
grep -Fq 'call write_nl_l_channel_tensor_rows_gpa(' "$write_f90"
grep -Fq 'real(8),allocatable :: stress_nl_l(:,:,:)' "$repo_root/src/common/structures.f90"
grep -Fq 'if(pp%inorm(l,ik) == 0) cycle' "$repo_root/src/common/stress.f90"
