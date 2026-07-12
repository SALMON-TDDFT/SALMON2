#!/usr/bin/env python3
from pathlib import Path

lcfo = Path("src/gs/dc/lcfo_flux.f90").read_text()
rtio = Path("src/rt/dg/rt_dg_fragment_io.f90").read_text()
inp = Path("src/io/inputoutput.f90").read_text()

for token in [
    "direct_amn_global",
    "basis source=direct_amn_projectors_global",
    "call read_wannier90_amn_direct_transform_import(dc, num_bands_chk, num_wann_chk, owner_frag, &",
    "call diagnose_wannier_transform_subspace_overlap(num_bands_chk, num_wann_chk, &",
    "call build_direct_wannier_buffer_position_from_lcfo_import(dc, num_bands_chk, num_wann_chk, &",
    "position source=direct_amn_buffer_local_r",
    "W90_vs_direct_amn_global",
]:
    if token not in lcfo:
        raise SystemExit(f"missing generic direct AMN global token: {token}")

direct_block = lcfo.split(
    "if(trim(dg_wannier_symmetry_gauge) == 'direct_amn_global') then\n"
    "          call read_wannier90_amn_direct_transform_import(dc, num_bands_chk, num_wann_chk, owner_frag, &",
    1,
)[1].split(
    "if(direct_gauge_requested .and. .not. direct_gauge_applied)", 1
)[0]
for bad_route in ["position source=direct_amn_lcfo_grid", "position source=direct_amn_phase_log"]:
    if bad_route in direct_block:
        raise SystemExit("direct_amn_global must use buffer-local continuous r, not sawtooth or phase-log position")

for bad_token in [
    "call set_wannier_centers_from_position_diagonal(num_wann_chk, center_bohr, aa_global)",
    "owner_frag(iw) = find_owner_fragment_from_center_import(dc, center_bohr(1:3,iw))",
    "call rebalance_wannier_owner_fragments_import(dc, center_bohr, owner_frag, num_wann_chk)",
]:
    if bad_token in direct_block:
        raise SystemExit("direct_amn_global must not reinterpret buffer-local position as global centers/owners")

if "call build_direct_wannier_buffer_position_from_lcfo_import(dc, num_bands_chk, num_wann_chk, &" in direct_block:
    if "owner_frag, center_bohr, v_direct_global, aa_direct, ok_position)" not in direct_block:
        raise SystemExit("buffer-local position builder must keep physical global centers/owners separate")

if "direct_gauge_applied = .true." not in direct_block:
    raise SystemExit("direct_amn_global must mark the direct AMN route as applied")

sym_block = lcfo.split("if(direct_gauge_requested .and. .not. direct_gauge_applied)", 1)[1].split(
    "if(allocated(symops)) deallocate(symops)", 1
)[0]
for token in [
    "requested direct AMN gauge was not applied",
    "trim(dg_wannier_symmetry_gauge) == 'local_inversion_position'",
    "call symmetrize_fragment_wannier_position_import(dc, center_bohr, owner_frag, num_wann_chk, &",
    "call diagnose_global_wannier_pbc_operator_symmetry(dc, center_bohr, num_wann_chk, &",
]:
    if token not in sym_block:
        raise SystemExit(f"missing direct AMN diagnostic control-flow token: {token}")

if "call symmetrize_global_wannier_position_import" in sym_block:
    raise SystemExit("direct AMN gauges must not trigger an automatic inversion projection")

for token in [
    "requested direct AMN gauge requires Wannier90 position data",
    "direct_amn_global operator diagnostic deferred:",
    "unitary symmetry representation required",
]:
    if token not in lcfo:
        raise SystemExit(f"missing direct-gauge diagnostic safety token: {token}")

post_direct_block = sym_block.split(
    'write(*,\'(1x,a,a)\') "[DC-LCFO-W90-SYM] position sym mode="', 1
)[1]
if "if(trim(dg_wannier_symmetry_gauge) == 'direct_amn_global') then" not in post_direct_block:
    raise SystemExit("direct_amn_global must not run a center-permutation diagnostic in a mismatched gauge")

diagnostic_block = lcfo.split(
    "subroutine diagnose_global_wannier_pbc_operator_symmetry", 1
)[1].split(
    "end subroutine diagnose_global_wannier_pbc_operator_symmetry", 1
)[0]
if '" z_odd_available=F"' not in diagnostic_block:
    raise SystemExit("unavailable position diagnostics must not be reported as z_odd_res=0")

for token in [
    "direct AMN bond gauge requires bond_centers projection",
    "requested direct_amn_global but AMN read failed",
]:
    if token not in lcfo:
        raise SystemExit(f"direct AMN route must fail loudly instead of silently no-oping: {token}")

if "direct_amn_global" not in inp:
    raise SystemExit("input validation does not allow generic direct_amn_global")

for token in [
    'binfile_wf_wannier_seed = "wavefunctions_wannier_seed.bin"',
    'binfile_w90g_persistent = "wannier90_global_basis_persistent.bin"',
    'binfile_w90seed_persistent = "wannier_flux_eigen_seed_persistent.bin"',
    "filename = wannier_wavefunction_seed_filename_import(ifrag)",
    "filename = wannier_wavefunction_seed_filename_import(1)",
    "filename = trim(base_directory)//binfile_wf_wannier_seed",
    "binfile_w90g_persistent",
    "binfile_w90seed_persistent",
]:
    if token not in lcfo:
        raise SystemExit(f"missing persistent Wannier seed protection token: {token}")

for token in [
    'binfile_wf_wannier_seed = "wavefunctions_wannier_seed.bin"',
    'binfile_w90g_persistent = "wannier90_global_basis_persistent.bin"',
    'binfile_w90seed_persistent = "wannier_flux_eigen_seed_persistent.bin"',
    "character(256) function dg_wavefunction_seed_filename",
    "character(256) function dg_total_persistent_filename",
    "filename = dg_wavefunction_seed_filename(bdir_frag, ifrag, binfile_wf, binfile_wf_wannier_seed)",
    "filename = dg_wavefunction_seed_filename(bdir_frag, 1, binfile_wf, binfile_wf_wannier_seed)",
    "filename = dg_total_persistent_filename(binfile_w90g, binfile_w90g_persistent)",
    "filename = dg_total_persistent_filename(binfile_w90seed, binfile_w90seed_persistent)",
    "use salmon_global, only: wannier_num_bands, wannier_num_wann",
    "DG-Fragment RT: incompatible Wannier90 global basis dimensions",
    "DG-Fragment RT: Wannier90 global basis exceeds wavefunction seed states",
]:
    if token not in rtio:
        raise SystemExit(f"missing RT persistent Wannier seed protection token: {token}")

if "Wannier90 global basis has more bands than RT wavefunction seed states" not in rtio:
    raise SystemExit("RT must fail loudly when a 256-state wavefunction seed is paired with a larger Wannier basis")

if "target_wann exceeds unique bond-center candidates" not in lcfo:
    raise SystemExit("bond-center projection map must reject duplicated bond-center projectors")

if "direct AMN projection matrix is rank deficient" not in lcfo:
    raise SystemExit("direct AMN import must reject rank-deficient projector spaces")

if "eval(k) > lambda_cut" not in lcfo:
    raise SystemExit("direct AMN Lowdin cutoff must use lambda_cut")
