#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)
structures_f90="$repo_root/src/common/structures.f90"
stress_f90="$repo_root/src/common/stress.f90"
write_f90="$repo_root/src/io/write.f90"

grep -F "real(8) :: stress_loc_fullobj_grad(3,3)" "$structures_f90" >/dev/null
grep -F "real(8) :: stress_loc_fullobj_diag(3,3)" "$structures_f90" >/dev/null

grep -F "system%stress_loc_fullobj_grad = strs_grad_sum" "$stress_f90" >/dev/null
grep -F "system%stress_loc_fullobj_diag(a,a) = energy%E_ion_loc / V" "$stress_f90" >/dev/null

grep -F 'P_loc_fullobj_grad [GPa]' "$write_f90" >/dev/null
grep -F 'P_loc_fullobj_diag [GPa]' "$write_f90" >/dev/null
grep -F 'P_loc_fullobj [GPa]' "$write_f90" >/dev/null
grep -F 'P_har+loc_fullobj+ewald [GPa]' "$write_f90" >/dev/null
