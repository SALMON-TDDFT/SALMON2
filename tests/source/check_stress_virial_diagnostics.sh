#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
write_f90="$repo_root/src/io/write.f90"

grep -Fq 'Stress virial diagnostics [Hartree]' "$write_f90"
grep -Fq 'Tr(kin)*V + 2E_kin' "$write_f90"
grep -Fq 'Tr(har)*V + E_h' "$write_f90"
grep -Fq 'Kinetic stress diagnostics [Hartree]' "$write_f90"
grep -Fq 'E_kin(from stress trace)' "$write_f90"
grep -Fq 'E_kin(grad2/2)' "$write_f90"
grep -Fq 'E_kin(cross)' "$write_f90"
grep -Fq 'E_kin(k2/2)' "$write_f90"
grep -Fq 'E_kin(reconstructed)' "$write_f90"
grep -Fq 'E_kin(reference)' "$write_f90"
grep -Fq 'Local/Ewald diagnostics' "$write_f90"
grep -Fq 'E_loc_sr' "$write_f90"
grep -Fq 'E_loc_lr' "$write_f90"
grep -Fq 'P_ewald_G [GPa]' "$write_f90"
grep -Fq 'P_ewald_R [GPa]' "$write_f90"
