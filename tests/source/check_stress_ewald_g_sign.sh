#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
stress_f90="$repo_root/src/common/stress.f90"

grep -F "strs_G(a,b) = strs_G(a,b) + fact * 2d0 * g(a) * g(b) / G2 * (1d0 + G2/(4d0*aEwald))" "$stress_f90" >/dev/null
grep -F "strs_G(a,a) = strs_G(a,a) - fact" "$stress_f90" >/dev/null
grep -F "strs_G_sum(a,a) = strs_G_sum(a,a) + pi * Qtot**2 / (2d0*aEwald*V**2)" "$stress_f90" >/dev/null
