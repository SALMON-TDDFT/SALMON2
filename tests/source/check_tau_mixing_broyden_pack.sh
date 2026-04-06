#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
structures_f90="$repo_root/src/common/structures.f90"
mixing_f90="$repo_root/src/gs/mixing.f90"

grep -Fq "tau_in" "$structures_f90"
grep -Fq "tau_out" "$structures_f90"
grep -Fq "tau_metric_weight" "$mixing_f90"
grep -Fq "wrapper_broyden" "$mixing_f90"
grep -Fq "pulay" "$mixing_f90"
grep -Fq "pack" "$mixing_f90"
grep -Fq "unpack" "$mixing_f90"
grep -Fq "sqrt(tau_metric_weight)" "$mixing_f90"
