#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)

grep -Fq 'type s_xc_operator_payload' "$repo_root/src/common/structures.f90"
grep -Fq 'logical :: use_tau_operator' "$repo_root/src/common/structures.f90"
grep -Fq 'logical :: use_laplacian_operator' "$repo_root/src/common/structures.f90"
grep -Fq 'payload%use_tau_operator' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'type(s_xc_operator_payload), intent(out), optional :: xc_payload' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'type(s_xc_operator_payload), intent(out), optional :: payload' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'xc_payload%use_tau_operator = .false.' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'check_builtin_r2scan_operator_gate.sh' "$repo_root/tests/source/CMakeLists.txt"
