#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
stress_f90="$repo_root/src/common/stress.f90"

if grep -Eq 'E_sr[[:space:]]*/[[:space:]]*3d0' "$stress_f90"; then
  echo "empirical local stress diagonal weight detected"
  exit 1
fi

grep -F "strs_sum(a,a) = strs_sum(a,a) + (E_sr + E_lr) / V" "$stress_f90" >/dev/null
grep -F "system%stress_loc_diag(a,a) = (E_sr + E_lr) / V" "$stress_f90" >/dev/null
grep -F "strs_sum = strs_sum / V" "$stress_f90" >/dev/null

if grep -Fq "strs_sum = strs_sum / V**2" "$stress_f90"; then
  echo "local stress gradient normalization still uses /V**2" >&2
  exit 1
fi
