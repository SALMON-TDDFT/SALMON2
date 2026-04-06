#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
input_f90="$repo_root/src/io/inputoutput.f90"
global_f90="$repo_root/src/io/salmon_global.f90"

grep -Fq "yn_tau_mixing" "$global_f90"
grep -Fq "tau_mixrate" "$global_f90"
grep -Fq "tau_metric_weight" "$global_f90"

grep -Fq "yn_tau_mixing" "$input_f90"
grep -Fq "tau_mixrate" "$input_f90"
grep -Fq "tau_metric_weight" "$input_f90"
grep -Fq "variables.log" "$input_f90"
