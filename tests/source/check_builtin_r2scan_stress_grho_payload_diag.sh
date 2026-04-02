#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
stress_src="$repo_root/src/common/stress.f90"
structures_src="$repo_root/src/common/structures.f90"
write_src="$repo_root/src/io/write.f90"

grep -Fq "type(s_vector) :: grho" "$structures_src"
grep -Fq "real(8) :: stress_xc_grad_payload(3,3)" "$structures_src"
grep -Fq "system%stress_xc_grad_payload = 0d0" "$stress_src"
grep -Fq "if (allocated(system%xc_payload%grho%v)) then" "$stress_src"
grep -Fq "system%stress_xc_grad_payload(a,b) = system%stress_xc_grad_payload(a,b) +" "$stress_src"
grep -Fq "call write_stress_tensor_row_gpa(fh, 'XC-grad-pyld'" "$write_src"
