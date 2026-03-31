#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
write_f90="$repo_root/src/io/write.f90"

grep -Fq 'Stress residual diagnostics [Hartree]' "$write_f90"
grep -Fq 'Tr(kin)*V + 2E_kin' "$write_f90"
grep -Fq 'virial_har = stress_tensor_trace(system%stress_har) * system%det_a + energy%E_h' "$write_f90"
grep -Fq 'Tr(har)*V + E_h' "$write_f90"
grep -Fq 'Tr(xc)*V + 3(E_vxc-E_xc)' "$write_f90"
grep -Fq 'Tr(loc_lr)*V + E_lr' "$write_f90"
grep -Fq 'Tr(ewa)*V + E_ion_ion' "$write_f90"
grep -Fq 'Residual_kin_har_xc' "$write_f90"

! rg -F -q 'Energy decomposition [Hartree]' "$write_f90"
! rg -F -q 'Total cancellation diagnostics [Hartree]' "$write_f90"
! rg -F -q 'Local/Ewald diagnostics' "$write_f90"
! rg -F -q 'NL virial diagnostics [Hartree]' "$write_f90"
! rg -F -q 'Kinetic stress diagnostics [Hartree]' "$write_f90"
