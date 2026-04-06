#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
xc_src="$repo_root/src/xc/salmon_xc.f90"
struct_src="$repo_root/src/common/structures.f90"
mixing_src="$repo_root/src/gs/mixing.f90"
scf_src="$repo_root/src/gs/scf_iteration.f90"
cmake_lists="$repo_root/tests/source/CMakeLists.txt"

grep -Fq 'check_builtin_r2scan_vtau_damping.sh' "$cmake_lists"
grep -Fq 'type(s_scalar) :: Vtau_in, Vtau_out' "$struct_src"
grep -Fq 'subroutine mix_xc_operator_payload(mg,system,c1,c2,mixing)' "$mixing_src"
grep -Fq 'if (.not. system%xc_payload%use_tau_operator) return' "$mixing_src"
grep -Fq 'if (.not. allocated(mixing%Vtau_in%f)) then' "$mixing_src"
grep -Fq 'mixing%Vtau_out%f = system%xc_payload%vtau%f' "$mixing_src"
grep -Fq 'system%xc_payload%vtau%f = c1*mixing%Vtau_in%f + c2*mixing%Vtau_out%f' "$mixing_src"
grep -Fq 'call mix_xc_operator_payload(mg,system,1.d0-mixing%mixrate,mixing%mixrate,mixing)' "$scf_src"
! grep -Fq 'r2scan_vtau_damping_alpha' "$xc_src"
! grep -Fq 'vtau_prev_local' "$xc_src"
