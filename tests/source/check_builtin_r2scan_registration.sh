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
grep -Fq 'module builtin_r2scan' "$repo_root/src/xc/builtin_r2scan.f90"
grep -Fq 'check_builtin_r2scan_registration.sh' "$repo_root/tests/source/CMakeLists.txt"
