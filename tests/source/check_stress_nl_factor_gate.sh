#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
stress_f90="$repo_root/src/common/stress.f90"

grep -F "contrib = 2d0 * rtmp * dble(conjg(uVpsi_ilma) * r_uVpsi_b(b))" "$stress_f90" >/dev/null

if grep -Fq "contrib = rtmp * dble(conjg(uVpsi_ilma) * r_uVpsi_b(b))" "$stress_f90"; then
  echo "nonlocal stress gradient still lacks the Hermitian factor 2" >&2
  exit 1
fi
