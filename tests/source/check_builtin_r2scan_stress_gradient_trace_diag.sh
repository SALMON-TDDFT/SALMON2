#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
stress_src="$repo_root/src/common/stress.f90"
structures_src="$repo_root/src/common/structures.f90"
write_src="$repo_root/src/io/write.f90"

grep -Fq "stress_xc_dbg_grho_local_payload_maxdiff" "$structures_src"
grep -Fq "stress_xc_dbg_grho_direct_payload_maxdiff" "$structures_src"
grep -Fq "stress_xc_dbg_grho_direct_local_maxdiff" "$structures_src"
grep -Fq "stress_xc_dbg_rdedd_dot_grho_local" "$structures_src"
grep -Fq "stress_xc_dbg_rdedd_dot_grho_payload" "$structures_src"
grep -Fq "stress_xc_dbg_rho_div_rdedd" "$structures_src"
grep -Fq "system%stress_xc_dbg_grho_direct_payload_maxdiff = maxdiff_out(2)" "$stress_src"
grep -Fq "system%stress_xc_dbg_grho_direct_local_maxdiff = maxdiff_out(3)" "$stress_src"
grep -Fq "system%stress_xc_dbg_rdedd_dot_grho_local = diag_out(1)" "$stress_src"
grep -Fq "system%stress_xc_dbg_rdedd_dot_grho_payload = diag_out(2)" "$stress_src"
grep -Fq "system%stress_xc_dbg_rho_div_rdedd = diag_out(3)" "$stress_src"
grep -Fq "max|grho_payload-grho_local|" "$write_src"
grep -Fq "max|grho_payload-grho_direct|" "$write_src"
grep -Fq "max|grho_direct-grho_local|" "$write_src"
grep -Fq "int rdedd.grho_payload + rho div(rdedd) dV" "$write_src"
