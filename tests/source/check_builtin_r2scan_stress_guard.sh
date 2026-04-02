#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)

if grep -Fq "stress tensor does not support xc='r2scan'" "$repo_root/src/common/stress.f90"; then
  echo "unexpected r2scan stress guard still present" >&2
  exit 1
fi

grep -Fq "call calc_stress_xc(system, pp, info, mg, stencil, srg, ppn, rho_s, Vxc, energy, xc_func, tpsi, field_state, srg_scalar)" "$repo_root/src/common/stress.f90"
grep -Fq "call calc_stress_xc_builtin_r2scan" "$repo_root/src/common/stress.f90"
grep -Fq "system%xc_payload%rdedd%v" "$repo_root/src/common/stress.f90"
grep -Fq "system%xc_payload%vtau%f" "$repo_root/src/common/stress.f90"
grep -Fq "r2SCAN supports only GS calculations" "$repo_root/src/xc/salmon_xc.f90"
grep -Fq 'check_builtin_r2scan_stress_guard.sh' "$repo_root/tests/source/CMakeLists.txt"
