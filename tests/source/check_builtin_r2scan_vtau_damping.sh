#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
xc_src="$repo_root/src/xc/salmon_xc.f90"
cmake_lists="$repo_root/tests/source/CMakeLists.txt"

grep -Fq 'check_builtin_r2scan_vtau_damping.sh' "$cmake_lists"
grep -Fq 'r2scan_vtau_damping_alpha = 0.5d0' "$xc_src"
grep -Fq "use salmon_global, only: yn_spinorbit, calc_mode" "$xc_src"
grep -Fq "trim(calc_mode) == 'GS'" "$xc_src"
grep -Fq 'allocate(vtau_prev_local(mg%num(1), mg%num(2), mg%num(3)))' "$xc_src"
grep -Fq 'vtau_prev_local = payload%vtau%f(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))' "$xc_src"
grep -Fq 'payload%vtau%f = r2scan_vtau_damping_alpha * payload%vtau%f' "$xc_src"
grep -Fq '+ (1d0 - r2scan_vtau_damping_alpha) * vtau_prev_local' "$xc_src"
