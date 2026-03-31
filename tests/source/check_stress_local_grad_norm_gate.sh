#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
stress_f90="$repo_root/src/common/stress.f90"

if grep -Eq 'strs_sum[[:space:]]*=[[:space:]]*strs_sum[[:space:]]*/[[:space:]]*V\*\*2' "$stress_f90"; then
  echo "local stress gradient still divides by V**2" >&2
  exit 1
fi

grep -Eq 'strs_sum[[:space:]]*=[[:space:]]*strs_sum[[:space:]]*/[[:space:]]*V([[:space:]]|$)' "$stress_f90"
