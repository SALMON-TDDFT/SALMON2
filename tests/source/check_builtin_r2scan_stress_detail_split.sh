#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
stress_src="$repo_root/src/common/stress.f90"
structures_src="$repo_root/src/common/structures.f90"
write_src="$repo_root/src/io/write.f90"

grep -Fq "stress_xc_local(3,3)" "$structures_src"
grep -Fq "stress_xc_grad(3,3)" "$structures_src"
grep -Fq "stress_xc_tau(3,3)" "$structures_src"
grep -Fq "system%stress_xc_local(a,a) = -(E_vxc - energy%E_xc) / V" "$stress_src"
grep -Fq "system%stress_xc_grad = strs_grad_sum / V" "$stress_src"
grep -Fq "system%stress_xc_tau = strs_tau_sum / V" "$stress_src"
grep -Fq "system%stress_xc = system%stress_xc_local + system%stress_xc_grad + system%stress_xc_tau" "$stress_src"
grep -Fq "call write_stress_tensor_row_gpa(fh, '    XC-local'" "$write_src"
grep -Fq "call write_stress_tensor_row_gpa(fh, '    XC-grad'" "$write_src"
grep -Fq "call write_stress_tensor_row_gpa(fh, '    XC-tau'" "$write_src"
