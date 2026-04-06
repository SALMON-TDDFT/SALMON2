#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
xc_f90="$repo_root/src/xc/salmon_xc.f90"
scf_f90="$repo_root/src/gs/scf_iteration.f90"
mixing_f90="$repo_root/src/gs/mixing.f90"

grep -Fq "tau_override" "$xc_f90"
grep -Fq "tau_override" "$scf_f90"
grep -Fq "mix_xc_operator_payload" "$mixing_f90"
grep -Fq "vtau" "$mixing_f90"
