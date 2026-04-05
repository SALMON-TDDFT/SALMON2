#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
stress_src="$repo_root/src/common/stress.f90"
structures_src="$repo_root/src/common/structures.f90"
write_src="$repo_root/src/io/write.f90"

grep -Fq "system%xc_payload%rdedd%v" "$stress_src"
grep -Fq "system%xc_payload%vtau%f" "$stress_src"
grep -Fq "stress_xc_dbg_grho_local_payload_maxdiff" "$structures_src"
grep -Fq "stress_xc_dbg_grho_direct_payload_maxdiff" "$structures_src"
grep -Fq "stress_xc_dbg_grho_direct_local_maxdiff" "$structures_src"
grep -Fq "max|grho_payload-grho_local|" "$write_src"
grep -Fq "max|grho_payload-grho_direct|" "$write_src"
grep -Fq "max|grho_direct-grho_local|" "$write_src"
grep -Fq "int rdedd.grho_local dV" "$write_src"
grep -Fq "int rdedd.grho_payload dV" "$write_src"
grep -Fq "int rho div(rdedd) dV" "$write_src"
