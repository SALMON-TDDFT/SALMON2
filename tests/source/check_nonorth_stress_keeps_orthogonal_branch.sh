#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
ham_f90="$repo_root/src/common/hamiltonian.f90"
cmake_lists="$repo_root/tests/source/CMakeLists.txt"

grep -Fq 'check_nonorth_stress_keeps_orthogonal_branch.sh' "$cmake_lists"
grep -Fq 'if (.not. stencil%if_orthogonal) then' "$ham_f90"
grep -Fq 'stop "error: nonorthogonal tau operator with real orbitals is unsupported"' "$ham_f90"
grep -Fq 'call add_xc_tau_operator_nonorth_complex' "$ham_f90"
grep -Fq '!$omp parallel do collapse(4) default(none)' "$ham_f90"
