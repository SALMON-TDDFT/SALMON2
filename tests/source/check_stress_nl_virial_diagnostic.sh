#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
write_f90="$repo_root/src/io/write.f90"

grep -Fq 'NL virial diagnostics [Hartree]' "$write_f90"
grep -Fq 'Tr(nl)*V              =' "$write_f90"
grep -Fq '3E_nl                 =' "$write_f90"
grep -Fq 'Tr(nl_grad)*V         =' "$write_f90"
grep -Fq 'P_nl_diag [GPa]       =' "$write_f90"
grep -Fq 'P_nl_grad [GPa]       =' "$write_f90"
