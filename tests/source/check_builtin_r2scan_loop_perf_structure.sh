#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
builtin_src="$repo_root/src/xc/builtin_r2scan.f90"
xc_src="$repo_root/src/xc/salmon_xc.f90"
cmake_lists="$repo_root/tests/source/CMakeLists.txt"

assert_absent() {
  pattern="$1"
  file="$2"
  if rg -F -q "$pattern" "$file"; then
    echo "unexpected pattern still present in $file: $pattern" >&2
    exit 1
  fi
}

grep -Fq 'check_builtin_r2scan_loop_perf_structure.sh' "$cmake_lists"
grep -Fq '!$omp parallel do default(none)' "$builtin_src"
grep -Fq 'real(8), intent(in) :: grho_norm(nl)' "$builtin_src"
grep -Fq 'call eval_r2scan_point(rho(i), rho_s(i), grho_norm(i), tau_s(i),' "$builtin_src"
grep -Fq 'real(8), intent(in) :: grho_norm' "$builtin_src"

grep -Fq 'grho_norm_1d = reshape(sqrt(grho(:,:,:,1)**2 + grho(:,:,:,2)**2 + grho(:,:,:,3)**2), (/nl/))' "$xc_src"
grep -Fq 'call exc_cor_r2scan(nl, rho_1d, rho_s_1d, grho_norm_1d, tau_s_1d,' "$xc_src"
grep -Fq 'if (grho_norm_1d(i) <= grad_floor) cycle' "$xc_src"
grep -Fq 'where (grho_norm_1d > grad_floor)' "$xc_src"

assert_absent 'grho_s(i,:)' "$builtin_src"
assert_absent 'real(8) :: grho_norm_3d(nx, ny, nz)' "$xc_src"
assert_absent 'real(8) :: vgrad_3d(nx, ny, nz)' "$xc_src"
