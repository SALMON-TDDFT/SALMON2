#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
stress_src="$repo_root/src/common/stress.f90"
structures_src="$repo_root/src/common/structures.f90"
write_src="$repo_root/src/io/write.f90"

grep -Fq "system%xc_payload%rdedd%v" "$stress_src"
grep -Fq "system%xc_payload%vtau%f" "$stress_src"
grep -Fq "stress_xc_dbg_rdedd_refresh_maxdiff" "$structures_src"
grep -Fq "stress_xc_dbg_rdedd_refresh_maxdiff" "$stress_src"
grep -Fq "call exchange_correlation(" "$stress_src"
grep -Fq "max|rdedd_payload-rdedd_xc_refresh|" "$write_src"
! rg -F -q "stress_xc_dbg_rdedd_payload_maxdiff" "$structures_src"
! rg -F -q "calc_r2scan_rdedd_payload_maxdiff" "$stress_src"
! rg -F -q "calc_r2scan_tau_density" "$stress_src"
! rg -F -q "max|rdedd_payload-rdedd_rebuilt|" "$write_src"
