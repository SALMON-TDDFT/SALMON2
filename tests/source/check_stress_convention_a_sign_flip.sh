#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
stress_f90="$repo_root/src/common/stress.f90"

grep -Fq 'Nielsen-Martin formula: sigma_{ab} = +(1/V) dE/d(eps_{ab})' "$stress_f90"
grep -Fq 'system%stress_har = -strs_sum' "$stress_f90"
grep -Fq 'system%stress_xc(a,a) = -(E_vxc - energy%E_xc) / V' "$stress_f90"
grep -Fq 'system%stress_loc_grad = -strs_grad_sum' "$stress_f90"
grep -Fq 'system%stress_loc_fullobj_grad = -strs_grad_sum' "$stress_f90"
grep -Fq 'system%stress_loc_sr_grad = -strs_sr_sum' "$stress_f90"
grep -Fq 'system%stress_loc_lr_grad = -strs_lr_sum' "$stress_f90"
grep -Fq 'system%stress_loc_sr_scr_grad = -(strs_grad_sum - strs_lr_scr_sum)' "$stress_f90"
grep -Fq 'system%stress_loc_lr_scr_grad = -strs_lr_scr_sum' "$stress_f90"
grep -Fq 'system%stress_loc_sr_diag(a,a) = -E_sr / V' "$stress_f90"
grep -Fq 'system%stress_loc_lr_diag(a,a) = -E_lr / V' "$stress_f90"
grep -Fq 'system%stress_loc_sr_scr_diag(a,a) = -E_sr_scr / V' "$stress_f90"
grep -Fq 'system%stress_loc_lr_scr_diag(a,a) = -E_lr_scr / V' "$stress_f90"
grep -Fq 'system%stress_loc_diag(a,a) = -strs_diag_sum(a,a)' "$stress_f90"
grep -Fq 'system%stress_loc_fullobj_diag(a,a) = -energy%E_ion_loc / V' "$stress_f90"
grep -Fq 'system%stress_loc = -strs_sum' "$stress_f90"
grep -Fq 'system%stress_nl = -strs_sum' "$stress_f90"
