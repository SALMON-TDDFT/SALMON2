#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)

grep -Fq 'type s_xc_operator_payload' "$repo_root/src/common/structures.f90"
grep -Fq 'logical :: use_tau_operator' "$repo_root/src/common/structures.f90"
grep -Fq 'logical :: use_laplacian_operator' "$repo_root/src/common/structures.f90"
grep -Fq 'logical :: vtau_has_shadow_values = .false.' "$repo_root/src/common/structures.f90"
grep -Fq 'payload%use_tau_operator' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'type(s_xc_operator_payload), intent(out), optional :: xc_payload' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'type(s_xc_operator_payload), intent(out), optional :: payload' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'xc_payload%use_tau_operator = .false.' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'xc_payload%vtau_has_shadow_values = .false.' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'payload%vtau_has_shadow_values = .false.' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'type(s_xc_operator_payload) :: xc_payload' "$repo_root/src/gs/scf_iteration.f90"
grep -Fq 'call exchange_correlation(' "$repo_root/src/gs/scf_iteration.f90"
grep -Fq 'xc_payload' "$repo_root/src/gs/initialization_dft.f90"
grep -Fq 'xc_payload' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'xc_payload=system%xc_payload' "$repo_root/src/gs/conjugate_gradient.f90"
grep -Fq 'xc_payload=system%xc_payload' "$repo_root/src/so/subspace_diagonalization_so.f90"
grep -Fq 'subroutine update_vlocal(' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'subroutine add_xc_tau_operator' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'xc_payload%vtau%f' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'if (.not. xc_payload%vtau_has_shadow_values) then' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'error: tau operator support is unavailable for OpenACC builds' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'if (.not. xc_payload%use_tau_operator) return' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'if (present(xc_payload)) then' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'else if (system%xc_payload%use_tau_operator) then' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'call add_xc_tau_operator(htpsi,tpsi,info,mg,system,stencil,srg,ppg,xc_payload)' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'call add_xc_tau_operator(htpsi,tpsi,info,mg,system,stencil,srg,ppg,system%xc_payload)' "$repo_root/src/common/hamiltonian.f90"
grep -Fq 'Vlocal(is)%f(ix,iy,iz) = Vpsl%f(ix,iy,iz) + Vh%f(ix,iy,iz) + Vxc(is)%f(ix,iy,iz)' "$repo_root/src/common/hamiltonian.f90"
update_vlocal_block=$(sed -n '/subroutine update_vlocal(/,/end subroutine update_vlocal/p' "$repo_root/src/common/hamiltonian.f90")
printf '%s\n' "$update_vlocal_block" | grep -Fq 'Vlocal(is)%f(ix,iy,iz) = Vpsl%f(ix,iy,iz) + Vh%f(ix,iy,iz) + Vxc(is)%f(ix,iy,iz)'
if printf '%s\n' "$update_vlocal_block" | grep -Eq 'xc_payload|vtau'; then
  exit 1
fi
grep -Fq 'check_builtin_r2scan_operator_gate.sh' "$repo_root/tests/source/CMakeLists.txt"
