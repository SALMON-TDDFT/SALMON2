#!/usr/bin/env python3
"""Contract for reusing authoritative DC controls in divided Hybrid SCF."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GLOBAL = (ROOT / "src/io/salmon_global.f90").read_text(errors="replace").lower()
INPUT = (ROOT / "src/io/inputoutput.f90").read_text(errors="replace").lower()
DCDTF = (ROOT / "src/gs/dc/dcdft.f90").read_text(errors="replace").lower()
MAIN = (ROOT / "src/gs/main_dft.f90").read_text(errors="replace").lower()

assert "character(1)   :: yn_dg_hybrid_divided_scf" in GLOBAL
assert "yn_dg_hybrid_divided_scf = 'n'" in INPUT, "divided route must default off"
assert "call comm_bcast(yn_dg_hybrid_divided_scf" in INPUT
assert "call yn_argument_check(yn_dg_hybrid_divided_scf)" in INPUT

for token in (
    "subroutine prepare_dg_hybrid_divided_dc_controls",
    "convergence_mode=convergence",
    "density_threshold=threshold",
    "allocate(initial_total_density(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3)",
    "do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)",
    "initial_total_density(ix,iy,iz)=dc%rho_tot%f(ix,iy,iz)",
    "call mpi_allreduce(mpi_in_place,initial_total_density,size(initial_total_density)",
    "subroutine load_dg_hybrid_distributed_dc_density",
    "call mpi_alltoallv",
):
    assert token in DCDTF, f"DC adapter is missing authoritative control: {token}"

assert "yn_dg_hybrid_divided_scf" in MAIN
branch_start = "if(yn_dg_hybrid_divided_scf=='y'.or.yn_dg_hybrid_continuation_scf=='y')then"
assert branch_start in MAIN, "missing default-off divided DC preparation branch"
branch = MAIN[MAIN.index(branch_start) :].split("endif", 1)[0]
assert "prepare_dg_hybrid_divided_dc_controls" in branch
assert "call prepare_dg_hybrid_divided_production_basis" in branch
production_call = MAIN.index("call prepare_dg_hybrid_divided_production_basis")
assert "subroutine prepare_dg_hybrid_divided_production_basis" in MAIN
main_end = MAIN.index("end subroutine main_dft")
prep_start = MAIN.index("subroutine prepare_dg_hybrid_divided_production_basis")
assert prep_start > main_end, (
    "production-basis preparation must not be an internal procedure of the giant "
    "main_dft frame"
)
prep = MAIN[MAIN.index("subroutine prepare_dg_hybrid_divided_production_basis") :].split(
    "end subroutine", 1
)[0]
for token in (
    "analyze_dg_hybrid_lcfo_selection",
    "freeze_dg_hybrid_production_selection",
    "redistribute_dg_hybrid_fragment_windows",
    "build_dg_hybrid_projected_fragment_basis",
):
    assert token in prep, f"divided route is missing production WF+PW construction: {token}"
assert prep.index("analyze_dg_hybrid_lcfo_selection") < prep.index(
    "freeze_dg_hybrid_production_selection"
)
lcfo_call = prep[prep.index("call analyze_dg_hybrid_lcfo_selection") :]
assert "basis_fingerprint_arg" in lcfo_call.split("if(.not.callback_ok)return", 1)[0], (
    "LCFO-deferred production preparation omits authoritative Wannier provenance"
)
assert prep.index("freeze_dg_hybrid_production_selection") < prep.index(
    "build_dg_hybrid_projected_fragment_basis"
)
assert "build_dg_hybrid_production_pw_basis" not in prep
assert "core_fragment_ids=dc%i_frag" not in prep, (
    "pencil-owned core points must not inherit one rank-local fragment id"
)
for token in (
    "fragment_origin_arg",
    "fragment_size_arg",
    "do fragment=1,fragment_count_arg",
):
    assert token in prep, f"missing physical-point fragment ownership reconstruction: {token}"
assert "row_action(size(ow_pencil_generator_maps,1)" not in prep, (
    "production row action must cover global physical points, not local pencil rows"
)
for token in (
    "row_action(ncore_arg,size(pencil_maps_arg,2))",
    "row_action=int(pencil_maps_arg)",
):
    assert token in prep, f"missing physical-ID row-action preservation: {token}"
assert "all_core_ids_arg" not in prep, (
    "global physical point IDs must not be reinterpreted as gathered-row ordinals"
)
assert "row_action=0" not in prep, (
    "the fully populated local row action must not use the crashing redundant bulk memset"
)
assert "core_fragment_ids=0" not in prep
assert "do point=1,ncore_arg" in prep
assert "core_fragment_ids(point)=0" in prep
assert "dc%rho_tot" in DCDTF
density_loader = DCDTF[DCDTF.index("subroutine load_dg_hybrid_distributed_dc_density") :]
density_loader = density_loader.split("end subroutine load_dg_hybrid_distributed_dc_density", 1)[0]
for token in (
    "point_multiplicity",
    "any(point_multiplicity/=1)",
    "distributed dc density catalog is not exactly once",
    "distributed dc density count exchange failed",
):
    assert token in density_loader, f"distributed density loader lacks catalog/error contract: {token}"
count_exchange = density_loader.index("call mpi_alltoall(send_counts")
assert count_exchange < density_loader.index("if(ierr/=mpi_success", count_exchange) < density_loader.index(
    "recv_displs(1)=0", count_exchange
)
callback_start = "subroutine apply_dg_hybrid_divided_fragment_hpsi"
assert callback_start in MAIN, "missing divided fragment Hamiltonian callback"
callback = MAIN[MAIN.index(callback_start) :].split("end subroutine", 1)[0]
for token in ("call hpsi", "mg", "v_local", "system", "ppg"):
    assert token in callback, f"divided Hamiltonian callback is missing fragment object: {token}"
for forbidden in ("dc%mg_tot", "dc%vloc_tot", "dc%system_tot", "dc%ppg_tot"):
    assert forbidden not in callback, f"divided Hamiltonian callback used total-system operator: {forbidden}"
solve_start = "subroutine solve_dg_hybrid_divided_fragments"
assert solve_start in MAIN, "missing divided fragment eigensolver adapter"
solve_callback = MAIN[MAIN.index(solve_start) :].split("end subroutine", 1)[0]
assert "solve_dg_hybrid_fragment_basis(dc%icomm_frag" in solve_callback
assert "solve_dg_hybrid_fragment_basis(dc%icomm_tot" not in solve_callback
mix_start = "subroutine mix_dg_hybrid_divided_density"
assert mix_start in MAIN, "missing divided DC density mixer adapter"
mix_callback = MAIN[MAIN.index(mix_start) :].split("end subroutine", 1)[0]
assert "call copy_density" in mix_callback
assert "select case(method_mixing)" in mix_callback
assert "mixing%mixrate" in mix_callback
for token in ("simple_mixing", "wrapper_broyden", "pulay"):
    assert token in mix_callback, f"divided adapter does not reuse DC mixer: {token}"
for forbidden in (
    "dg_hybrid_divided_density_tolerance",
    "dg_dc_gs_final_density_tolerance",
    "post_lcfo_density",
):
    assert forbidden not in branch, f"divided route introduced a new density gate: {forbidden}"

potential_name = "subroutine finish_dg_dc_potential_update"
potential = MAIN[MAIN.index(potential_name) :].split(
    "end subroutine finish_dg_dc_potential_update", 1
)[0]
assert "call hartree(dc%lg_tot" in potential, "hybrid Hartree must use the established total FFT grid"
assert "call exchange_correlation(system,xc_func,mg" in potential, (
    "hybrid XC must be evaluated on the fragment density and halo"
)
assert "call exchange_correlation(dc%system_tot" not in potential, (
    "hybrid XC must not create a second full-system evaluation"
)
assert "call calc_vlocal_fragment_dcdft" in potential, "total Hartree is not returned to fragments"

continuation = MAIN.split("subroutine run_dg_hybrid_concrete_continuation", 1)[1].split(
    "end subroutine run_dg_hybrid_concrete_continuation", 1
)[0]
assert "dg_dc_update_potential_from_distributed_density" in continuation
assert "allocate(density4(dc%lg_tot%num" not in continuation
assert "gather_dg_hybrid_divided_core_density(rho_in" not in continuation

print("divided Hybrid DC controls contract: PASS")
