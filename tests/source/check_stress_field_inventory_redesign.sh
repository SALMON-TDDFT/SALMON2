#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
structures_src="$repo_root/src/common/structures.f90"

assert_absent() {
  pattern="$1"
  if rg -F -q "$pattern" "$structures_src"; then
    echo "unexpected pattern still present: $pattern" >&2
    exit 1
  fi
}

grep -Fq "real(8) :: stress_tensor(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_kin(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_har(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_xc(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_loc(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_nl(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_ewa(3,3)" "$structures_src"

grep -Fq "real(8) :: stress_xc_e_vxc" "$structures_src"
grep -Fq "real(8) :: stress_loc_sr_grad(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_loc_lr_grad(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_loc_sr_diag(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_loc_lr_diag(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_loc_sr_energy" "$structures_src"
grep -Fq "real(8) :: stress_loc_lr_energy" "$structures_src"
grep -Fq "real(8) :: stress_ewa_g(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_ewa_r(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_ewa_energy_G" "$structures_src"
grep -Fq "real(8) :: stress_ewa_energy_R" "$structures_src"

grep -Fq "real(8) :: stress_har_shadow(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_xc_local(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_xc_grad(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_xc_grad_payload(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_xc_grad_vsigma(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_xc_tau(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_x(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_c(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_loc_sr_scr_grad(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_loc_lr_scr_grad(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_loc_sr_scr_diag(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_loc_lr_scr_diag(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_loc_sr_rs(3,3)" "$structures_src"
grep -Fq "real(8),allocatable :: stress_nl_species_l(:,:,:,:)" "$structures_src"
grep -Fq "real(8) :: stress_loc_sr_scr_energy" "$structures_src"
grep -Fq "real(8) :: stress_loc_lr_scr_energy" "$structures_src"
grep -Fq "real(8) :: stress_ewa_g_grad(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_ewa_g_diag(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_ewa_g_self(3,3)" "$structures_src"
grep -Fq "real(8) :: stress_xc_dbg_grho_local_payload_maxdiff" "$structures_src"
grep -Fq "real(8) :: stress_xc_dbg_grho_direct_payload_maxdiff" "$structures_src"
grep -Fq "real(8) :: stress_xc_dbg_grho_direct_local_maxdiff" "$structures_src"
grep -Fq "real(8) :: stress_xc_dbg_rdedd_dot_grho_local" "$structures_src"
grep -Fq "real(8) :: stress_xc_dbg_rdedd_dot_grho_payload" "$structures_src"
grep -Fq "real(8) :: stress_xc_dbg_rho_div_rdedd" "$structures_src"

assert_absent "stress_loc_fd"
assert_absent "stress_loc_grad"
assert_absent "stress_loc_diag"
assert_absent "stress_loc_fullobj_grad"
assert_absent "stress_loc_fullobj_diag"
assert_absent "stress_x_e_vx"
assert_absent "stress_c_e_vc"
assert_absent "stress_kin_dbg_grad2"
assert_absent "stress_kin_dbg_cross"
assert_absent "stress_kin_dbg_k2"
assert_absent "stress_xc_dbg_rdedd_refresh_maxdiff"
assert_absent "stress_xc_dbg_grho_refresh_maxdiff"
assert_absent "stress_xc_dbg_rho_box_direct_maxdiff"
assert_absent "stress_xc_dbg_rho_box_direct_active_maxdiff"
assert_absent "stress_xc_dbg_grho_direct_local_early_maxdiff"
assert_absent "stress_xc_dbg_grho_direct_local_bulk_maxdiff"
assert_absent "real(8),allocatable :: stress_nl_l(:,:,:)"
