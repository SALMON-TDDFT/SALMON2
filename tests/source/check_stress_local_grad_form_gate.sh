#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
stress_f90="$repo_root/src/common/stress.f90"

grep -F "coeff_lr = -2d0 * dble(conjg(rho_e) * V_lr_sum) / G2" "$stress_f90" >/dev/null
grep -F "strs_sr(a,b) = strs_sr(a,b) + 2d0 * dble(conjg(rho_e) * dVsr_dG2_sum) * g(a) * g(b)" "$stress_f90" >/dev/null

if grep -Fq "coeff_lr = 2d0 * dble(conjg(rho_e) * V_lr_sum) / G2" "$stress_f90"; then
  echo "local LR gradient sign still uses the pre-QE mapping" >&2
  exit 1
fi

if grep -Fq "strs_sr(a,b) = strs_sr(a,b) - 2d0 * dble(conjg(rho_e) * dVsr_dG2_sum) * g(a) * g(b)" "$stress_f90"; then
  echo "local SR gradient sign still uses the pre-QE mapping" >&2
  exit 1
fi
