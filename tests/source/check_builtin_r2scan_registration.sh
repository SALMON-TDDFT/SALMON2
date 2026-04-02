#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)

# Keep this aligned with the current baseline skeleton only.
grep -Fq 'builtin_r2scan.f90' "$repo_root/src/xc/CMakeLists.txt"
grep -Fq "case ('r2scan')" "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'salmon_xctype_r2scan' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'use builtin_r2scan' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'call exec_builtin_r2scan()' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'subroutine exec_builtin_r2scan()' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'stop "r2SCAN supports only nspin=1"' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'rho_1d = reshape(rho, (/nl/))' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'grho_s_1d = reshape(grho(:, :, :, :), (/nl, 3/)) * 0.5' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'tau_s_1d = reshape(tau(:, :, :), (/nl/)) * 0.5' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'call exc_cor_r2scan' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'payload%use_tau_operator = .true.' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'payload%vtau%f' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'rdedd(ix,iy,iz,idir)' "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'module builtin_r2scan' "$repo_root/src/xc/builtin_r2scan.f90"
grep -Fq 'subroutine exc_cor_r2scan' "$repo_root/src/xc/builtin_r2scan.f90"
grep -Fq 'real(8), intent(out) :: vtau(nl)' "$repo_root/src/xc/builtin_r2scan.f90"
grep -Fq 'real(8), intent(out) :: vgrad(nl)' "$repo_root/src/xc/builtin_r2scan.f90"
grep -Fq 'check_builtin_r2scan_registration.sh' "$repo_root/tests/source/CMakeLists.txt"
