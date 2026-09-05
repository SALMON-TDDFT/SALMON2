#!/usr/bin/env python3
"""Contract for reusing authoritative DC controls in divided Hybrid SCF."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
GLOBAL = (ROOT / "src/io/salmon_global.f90").read_text(errors="replace").lower()
INPUT = (ROOT / "src/io/inputoutput.f90").read_text(errors="replace").lower()
DCDTF = (ROOT / "src/gs/dc/dcdft.f90").read_text(errors="replace").lower()
MAIN = (ROOT / "src/gs/main_dft.f90").read_text(errors="replace").lower()
BROKEN = (ROOT / "src/gs/dc/dg_hybrid_broken_volume.f90").read_text(errors="replace").lower()
MIXING = (ROOT / "src/gs/dc/dg_hybrid_divided_mixing.f90").read_text(errors="replace").lower()

assert "character(1)   :: yn_dg_hybrid_divided_scf" in GLOBAL
assert "yn_dg_hybrid_divided_scf = 'n'" in INPUT, "divided route must default off"
assert "call comm_bcast(yn_dg_hybrid_divided_scf" in INPUT
assert "call yn_argument_check(yn_dg_hybrid_divided_scf)" in INPUT
assert "character(16)  :: dg_hybrid_divided_mixing" in GLOBAL
dc_namelist = INPUT[INPUT.index("namelist/dc/") : INPUT.index("!! == default for &unit")]
assert "dg_hybrid_divided_mixing" in dc_namelist
assert "dg_hybrid_divided_mixing = 'pulay'" in INPUT
assert "call string_lowercase(dg_hybrid_divided_mixing)" in INPUT
assert "call comm_bcast(dg_hybrid_divided_mixing" in INPUT
assert "'dg_hybrid_divided_mixing',trim(dg_hybrid_divided_mixing)" in INPUT
assert "select case(trim(dg_hybrid_divided_mixing))" in INPUT
assert "case('inherit','simple','pulay','broyden')" in INPUT

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

dc_occupation = DCDTF[DCDTF.index("subroutine ne2mu_dcdft") :].split(
    "end subroutine ne2mu_dcdft", 1
)[0]
assert "use occupation_kernel, only: solve_weighted_state_occupations" in dc_occupation
assert "call solve_weighted_state_occupations" in dc_occupation
assert "state_weights=reshape(ne_frag_orb,[state_count])" in dc_occupation
assert "state_weights=0.5d0*state_weights" in dc_occupation
assert "rocc=reshape(solved_occupations,[system%no,system%nspin,dc%n_frag])" in dc_occupation
assert "system%rocc(1:system%no,1,1:system%nspin)=rocc(:,:,dc%i_frag)" in dc_occupation
assert dc_occupation.index("call calc_ne_each") < dc_occupation.index(
    "call comm_summation(wrk1,esp"
) < dc_occupation.index("state_weights=reshape(ne_frag_orb,[state_count])")
assert dc_occupation.index("state_weights=reshape(ne_frag_orb,[state_count])") < dc_occupation.index(
    "call solve_weighted_state_occupations"
) < dc_occupation.index("rocc=reshape(solved_occupations") < dc_occupation.index(
    "system%rocc(1:system%no,1,1:system%nspin)=rocc(:,:,dc%i_frag)"
)
assert "dc%elec_num_tot,max(0d0,temperature),wspin" in dc_occupation
assert "subroutine mu2ne" not in dc_occupation
assert "subroutine ne2mu_core" not in dc_occupation

assert "yn_dg_hybrid_divided_scf" in MAIN
assert "ok=nproc==dc%n_frag.and..not.dc%optimized_fragment_geometry" in MAIN, (
    "production DC handoff must retain exactly one MPI rank per fragment"
)
rank_guard = MAIN.index("ok=nproc==dc%n_frag.and..not.dc%optimized_fragment_geometry")
assert "call comm_logical_and(ok,reusable,dc%icomm_tot)" in MAIN[rank_guard : rank_guard + 300]
assert "requires one rank per valid dc fragment" in MAIN[rank_guard : rank_guard + 400]
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
entry = "subroutine run_dg_hybrid_divided_ground_state_for_main"
assert entry in MAIN, "missing separated bounded divided production entry"
bounded = MAIN[MAIN.index(entry) :].split("end subroutine run_dg_hybrid_divided_ground_state_for_main", 1)[0]
for token in (
    "initialize_dg_hybrid_schwarz_state",
    "build_dg_hybrid_schwarz_schedule",
    "solve_dg_hybrid_schwarz_fragments",
):
    assert token in bounded, f"Schwarz divided production entry is missing {token}"
assert "extract_dg_hybrid_fragment_self_block" not in bounded, (
    "production divided SCF still extracts a self block instead of applying full DG rows"
)
for forbidden in (
    "apply_dg_hybrid_divided_fragment_hpsi",
    "solve_dg_hybrid_fragment_spectrum",
    "determine_dc_fragment_occupations",
):
    assert not re.search(r"\bcall\s+" + forbidden + r"\b", bounded), (
        f"bounded divided production entry still calls legacy {forbidden}"
    )
bounded_solver_name = "subroutine solve_dg_hybrid_schwarz_fragments"
assert bounded_solver_name in MAIN, "missing bounded Schwarz fragment callback"
bounded_solver = MAIN[MAIN.index(bounded_solver_name) :].split(
    "end subroutine solve_dg_hybrid_schwarz_fragments", 1
)[0]
assert "advance_dg_hybrid_schwarz_epoch" in bounded_solver
assert "assign_dg_hybrid_schwarz_occupations" in bounded_solver
assert "dg_hybrid_fragment_cg_steps" in bounded_solver and "temperature" in bounded_solver
assert "apply_dg_hybrid_schwarz_h" in bounded_solver
assert "apply_dg_hybrid_schwarz_s" in bounded_solver
assert "local_iterations<=dg_hybrid_fragment_cg_steps" in bounded_solver
for token in (
    "apply_dg_hybrid_schwarz_hamiltonian",
    "apply_dg_hybrid_schwarz_rows",
    "neighbor_exchanges=",
    "accepted_cg_steps=",
    "common_extensions=",
):
    assert token in MAIN, f"production Schwarz diagnostics/route is missing {token}"
local_projection = BROKEN[BROKEN.index("subroutine assemble_dg_hybrid_local_potential_rows") :].split(
    "end subroutine assemble_dg_hybrid_local_potential_rows", 1
)[0]
assert "partial(i,columns)=matmul" in re.sub(r"\s+", "", local_projection), (
    "density-epoch local-potential projection must use the fragment block matrix product"
)
for token in (
    "single_fragment_owner",
    "local_rows(:,columns)=matmul",
):
    assert token in re.sub(r"\s+", "", local_projection), (
        "one-rank/one-fragment density epochs must avoid redundant global reductions: " + token
    )
thermal = (ROOT / "src/gs/dc/dg_hybrid_fragment_thermal.f90").read_text(errors="replace").lower()
assert "run_dc_fragment_occupation_epoch" in thermal
assert "reconstruct_dg_hybrid_fragment_density" in thermal
mix_start = "subroutine mix_dg_hybrid_divided_density"
assert mix_start in MAIN, "missing divided DC density mixer adapter"
mix_callback = MAIN[MAIN.index(mix_start) :].split("end subroutine", 1)[0]
for token in (
    "prepare_dg_hybrid_divided_mixing",
    "accept_dg_hybrid_divided_mixing",
    "mixing_iteration",
    "reset_reason",
    "history_length",
):
    assert token in mix_callback, f"divided mixer lacks persistent lifecycle: {token}"
assert "call copy_density" in mix_callback
assert "dg_hybrid_divided_mixing,method_mixing" in mix_callback
assert "select case(selected_mixing_method)" in mix_callback
assert "mixing%mixrate" in mix_callback
for token in ("simple_mixing", "wrapper_broyden", "pulay"):
    assert token in mix_callback, f"divided adapter does not reuse DC mixer: {token}"
assert "copy_density(mixing_iteration" in mix_callback
assert "wrapper_broyden" in mix_callback and "mixing_iteration,mixing" in mix_callback
assert "pulay" in mix_callback and "mixing_iteration,mixing" in mix_callback
for token in (
    "basis_generation",
    "common_inventory",
    "collective_rollback",
    "unsupported divided hybrid mixing method",
):
    assert token in MIXING, f"divided mixing lifecycle lacks contract: {token}"
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

single_owner = re.search(
    r"if\s*\(\s*yn_dg_hybrid_divided_scf\s*==\s*'y'\s*\)\s*then\s*"
    r"call\s+freeze_dg_hybrid_single_owner_payload\b(.*?)endif", MAIN, re.S
)
assert single_owner, "divided main must freeze its single-owner column directory with the fixed payload"
for token in ("divided_fragment_basis", "divided_basis_local_slot", "divided_basis_generation",
              "divided_basis_directory_fingerprint", "dg_hybrid_fixed_payload"):
    assert token in single_owner.group(1), f"single-owner payload handoff lacks {token}"
assert "call freeze_dg_hybrid_variational_payload" in single_owner.group(1).split("else", 1)[1], (
    "continuation must retain its existing payload publication path"
)

print("divided Hybrid DC controls contract: PASS")
