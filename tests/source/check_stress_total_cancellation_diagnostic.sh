#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)
structures_f90="$repo_root/src/common/structures.f90"
stress_f90="$repo_root/src/common/stress.f90"
write_f90="$repo_root/src/io/write.f90"

if ! rg -q "stress_xc_e_vxc" "$structures_f90"; then
  echo "missing stress_xc_e_vxc field in structures.f90" >&2
  exit 1
fi

if ! rg -q "system%stress_xc_e_vxc *= *E_vxc" "$stress_f90"; then
  echo "missing E_vxc capture in calc_stress_xc" >&2
  exit 1
fi

if ! rg -F -q 'Total cancellation diagnostics [Hartree]' "$write_f90"; then
  echo "missing total cancellation diagnostics heading" >&2
  exit 1
fi

if ! rg -F -q 'Tr(xc)*V - 3(E_vxc-E_xc)' "$write_f90"; then
  echo "missing XC virial residual line" >&2
  exit 1
fi

if ! rg -F -q 'RHS_known(kin+har+xc+nl)' "$write_f90"; then
  echo "missing known-block theorem RHS line" >&2
  exit 1
fi

if ! rg -q 'Remainder_after_known' "$write_f90"; then
  echo "missing remainder-after-known line" >&2
  exit 1
fi
