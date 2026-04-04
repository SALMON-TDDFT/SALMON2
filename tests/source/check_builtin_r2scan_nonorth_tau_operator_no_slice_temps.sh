#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
ham_f90="$repo_root/src/common/hamiltonian.f90"
cmake_lists="$repo_root/tests/source/CMakeLists.txt"

grep -Fq 'check_builtin_r2scan_nonorth_tau_operator_no_slice_temps.sh' "$cmake_lists"
grep -Fq 'subroutine add_xc_tau_operator_nonorth_complex' "$ham_f90"

if grep -Fq 'dpsi(1,:,:,:)' "$ham_f90" || \
   grep -Fq 'dpsi(2,:,:,:)' "$ham_f90" || \
   grep -Fq 'dpsi(3,:,:,:)' "$ham_f90"; then
  echo "nonorth tau operator still passes noncontiguous dpsi slices to helper calls" >&2
  exit 1
fi
