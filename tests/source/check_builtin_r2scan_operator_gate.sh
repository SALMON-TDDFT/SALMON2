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
grep -Fq 'type(s_xc_operator_payload) :: xc_payload' "$repo_root/src/gs/scf_iteration.f90"
grep -Fq 'call exchange_correlation(' "$repo_root/src/gs/scf_iteration.f90"
grep -Fq 'xc_payload' "$repo_root/src/gs/initialization_dft.f90"
grep -Fq 'xc_payload' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'subroutine update_vlocal(' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'subroutine add_xc_tau_operator' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'check_builtin_r2scan_operator_gate.sh' "$repo_root/tests/source/CMakeLists.txt"
