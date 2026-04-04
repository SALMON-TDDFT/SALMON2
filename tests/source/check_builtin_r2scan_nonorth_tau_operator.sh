#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)

if grep -Fq 'stop "error: tau operator support requires orthogonal stencil"' \
  "$repo_root/src/common/hamiltonian.f90"; then
  echo "unexpected orthogonal-only tau operator guard still present" >&2
  exit 1
fi

grep -Fq 'call add_xc_tau_operator_nonorth_complex' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'stencil%coef_F(4)' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'check_builtin_r2scan_nonorth_tau_operator.sh' "$repo_root/tests/source/CMakeLists.txt"
