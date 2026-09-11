#!/usr/bin/env python3
"""Run a real divided-Hybrid GS producer and consume its distributed v5 checkpoint in Exp RT."""

from pathlib import Path
import hashlib
import math
import os
import re
import shutil
import subprocess
import tempfile


root = Path(__file__).resolve().parents[2]
main_source = (root / "src/rt/main_tddft.f90").read_text()
main_dft_source = (root / "src/gs/main_dft.f90").read_text()

formal_publisher = main_dft_source.split(
    "subroutine publish_dg_hybrid_divided_v5", 1
)[1].split("end subroutine publish_dg_hybrid_divided_v5", 1)[0]
def require_localized_publisher(body: str) -> None:
    compact = re.sub(r"\s+", "", body.lower())
    assert "allocate(payload%initial_occupied_amplitudes,source=occupied_state%coefficients)" in compact, \
        "terminal LCFO eigenvectors are not stored by owned localized-basis rows"
    assert "full_coefficients" not in compact and "collect_dg_hybrid_full_rows" not in compact, \
        "formal v5 publisher reconstructs a global coefficient matrix"
    assert "allocate(payload%basis_point_offsets(npoint+1)" in compact and "payload%basis_support_values(slot)=basis_values(j,p)" in compact, \
        "formal v5 publisher does not retain the localized construction basis as point CSR"
    assert "global_projector_count=dc%ppg_tot%nlma" in compact and \
        "real(global_projector_count,8)" in compact, \
        "formal v3 pseudopotential receipt does not use the authoritative full-system projector count"
    assert "mpi_allreduce(ppg%nlma,global_projector_count" not in compact, \
        "formal v3 pseudopotential receipt double-counts overlapping fragment-buffer projectors"
    assert "payload%energy_receipt=[checkpoint_energy%e_tot" in compact, \
        "formal v5 checkpoint does not contain a physical GS energy decomposition"
    assert "callcalc_total_energy_periodic" in compact, \
        "formal v5 energy reference is not produced by SALMON total-energy evaluation"
    assert "payload%energy_receipt=0d0" not in compact, \
        "all-zero energy receipt silently disables GS-to-RT energy identity"
    assert "payload%system_fingerprint=fingerprint_rt_dg_hybrid_system" in compact, \
        "formal v5 checkpoint is not bound to the GS physical system"
    assert "identity_system=dc%system_tot" in compact and \
        "fingerprint_rt_dg_hybrid_system(identity_system," in compact and \
        "nint(sum(occupied_state%occupations)),[0,0]" in compact and "identity_system%rocc" not in compact, \
        "formal checkpoint does not use the common physical electronic specification"
    assert "payload%pseudopotential_fingerprint=canonical_pp_fingerprint(pp)" in compact, \
        "formal v5 checkpoint is not bound to the canonical GS pseudopotential"
    assert "payload%pseudopotential_digest=canonical_pp_digest(pp)" in compact, \
        "formal checkpoint does not independently authenticate all PP tables"


require_localized_publisher(formal_publisher)
for old, replacement in (
    ("allocate(payload%initial_occupied_amplitudes,source=occupied_state%coefficients)",
     "allocate(payload%initial_occupied_amplitudes,source=solved_coefficients)"),
    ("allocate(payload%basis_point_offsets(npoint+1)", "allocate(payload%removed_point_offsets(npoint+1)"),
    ("global_projector_count=dc%ppg_tot%Nlma", "global_projector_count=ppg%Nlma"),
    ("payload%energy_receipt=[checkpoint_energy%E_tot", "payload%energy_receipt=0d0"),
    ("payload%pseudopotential_fingerprint=canonical_pp_fingerprint(pp)",
     "payload%pseudopotential_fingerprint=1_8"),
    ("payload%pseudopotential_digest=canonical_pp_digest(pp)",
     "payload%pseudopotential_digest=1_8"),
):
    mutated = formal_publisher.replace(old, replacement, 1)
    assert mutated != formal_publisher, old
    try:
        require_localized_publisher(mutated)
    except AssertionError:
        pass
    else:
        raise AssertionError(f"localized-basis mutation survived: {old}")


def require_hybrid_route_contract(source: str) -> None:
    continuation = source.split("subroutine run_dg_hybrid_continuation_rt()", 1)[1].split(
        "end subroutine run_dg_hybrid_continuation_rt", 1
    )[0]
    projection = source.split("subroutine project_salmon_local_rows", 1)[1].split(
        "end subroutine project_salmon_local_rows", 1
    )[0]
    assert "[HYBRID-RT-ROUTE] propagator=EXP potential=PP+HARTREE+XC" in continuation
    assert "[HYBRID-RT-STEP] step=" in continuation and "path=PP+HARTREE+XC" in continuation
    assert "[HYBRID-RT-STATE]" in continuation
    assert "[HYBRID-RT-PROJECTION]" in continuation
    assert "[HYBRID-RT-CHECKPOINT-STATIONARITY]" in continuation
    assert "[HYBRID-RT-REFRESH-STATIONARITY]" in continuation
    assert "establish_fixed_density_reference=.true." in continuation.lower()
    assert "all(hybrid_state%energy_receipt==0d0)" in continuation.lower()
    assert "hybrid dg rt physical energy receipt is invalid" in continuation.lower()
    assert "energy%e_ion_ion=hybrid_state%energy_receipt(5)" in continuation.lower()
    density_update_source = (root / "src/rt/dg/rt_dg_hybrid_density_update.f90").read_text().lower()
    assert "new_h=new_h+state%hamiltonian_reference_correction" in density_update_source
    assert "call propagate_rt_dg_hybrid_length_gauge" in continuation
    assert "canonical_pp_fingerprint(pp)" in continuation.lower()
    assert "fingerprint_rt_dg_hybrid_system" in continuation.lower()
    for token in ("call hartree", "call exchange_correlation_density", "call update_vlocal"):
        assert token in projection.lower(), token


def h4_gs_input() -> str:
    return """&calculation
 theory='dft'
 yn_dc='y'
/
&control
 sysname='h4_hybrid_smoke'
 yn_reset_step_restart='y'
 write_gs_restart_data='no'
/
&units
 unit_system='a.u.'
/
&dc
 num_fragment=1,1,4
 num_rgrid_buffer=0,0,4
 nproc_rgrid_tot=1,1,4
 nstate_frag=20
 energy_cut=10d0
 yn_dc_lcfo='n'
 yn_dc_lcfo_flux='n'
 yn_dc_lcfo_diag='n'
 yn_dc_fragment_optimization='n'
 yn_dg_dc_overlapping_wannier='y'
 yn_dg_hybrid_scf='n'
 yn_dg_hybrid_divided_scf='y'
 dg_hybrid_divided_mixing='pulay'
 dg_dc_seed_mode='off'
 dg_fragment_wf_checkpoint_mode='off'
 dg_fragment_w90_initial_projection='scdm'
 dg_ow_candidate_states_per_fragment=4
 dg_ow_target_wanniers_per_fragment=0
 dg_dc_gs_density_mix_rate=0.2d0
 dg_dc_gs_maximum_scf_iterations=12
 dg_dc_gs_final_density_tolerance=1d-4
 dg_dc_gs_final_orbital_tolerance=1d-4
 dg_dc_gs_electron_count_tolerance=1d-6
 dg_ow_boundary_value_tolerance=1d-3
 dg_ow_boundary_gradient_tolerance=1d-3
 dg_ow_symmetry_tolerance=1d-8
 dg_ow_localization_gradient_tolerance=2d-2
 wannier_num_iter=500
 wannier_pw_cutoff=0.5d0
/
&parallel
 nproc_k=1
 nproc_ob=1
 nproc_rgrid=1,1,1
 yn_eigenexa='n'
 yn_scalapack='y'
/
&system
 yn_periodic='y'
 al=12d0,12d0,18d0
 nelem=1
 nstate=2
 nelec=4
 natom=4
 temperature_k=300d0
/
&pseudo
 izatom(1)=1
 file_pseudo(1)='H_rps.dat'
 lloc_ps(1)=0
/
&functional
 xc='PZ'
/
&rgrid
 num_rgrid=16,16,24
/
&kgrid
 num_kgrid=1,1,1
/
&scf
 nscf_init_redistribution=2
 nscf_init_no_diagonal=0
 nscf=60
 ncg=3
 method_mixing='simple'
 mixrate=0.2d0
 threshold=1d-4
 yn_preconditioning='y'
/
&atomic_coor
 'H' 0d0 0d0 -4.7d0 1
 'H' 0d0 0d0 -3.3d0 1
 'H' 0d0 0d0  3.3d0 1
 'H' 0d0 0d0  4.7d0 1
/
"""


def h4_rt_input() -> str:
    return """&calculation
 theory='tddft_response'
/
&control
 sysname='h4_hybrid_rt_smoke'
/
&units
 unit_system='a.u.'
/
&parallel
 nproc_k=1
 nproc_ob=1
 nproc_rgrid=1,1,4
 yn_eigenexa='n'
/
&system
 yn_periodic='y'
 al=12d0,12d0,18d0
 nelem=1
 nstate=2
 nelec=4
 natom=4
/
&pseudo
 izatom(1)=1
 file_pseudo(1)='H_rps.dat'
 lloc_ps(1)=0
/
&functional
 xc='PZ'
/
&rgrid
 num_rgrid=16,16,24
/
&kgrid
 num_kgrid=1,1,1
/
&tgrid
 dt=0.02d0
 nt=2
/
&propagation
 yn_rt_dg_hybrid_continuation='y'
 yn_dg_length_gauge='y'
/
&emfield
 ae_shape1='impulse'
 e_impulse=1d-4
 epdir_re1=0d0,0d0,1d0
/
&atomic_coor
 'H' 0d0 0d0 -4.7d0 1
 'H' 0d0 0d0 -3.3d0 1
 'H' 0d0 0d0  3.3d0 1
 'H' 0d0 0d0  4.7d0 1
/
"""


def h4_rt_zero_input() -> str:
    return h4_rt_input().replace("ae_shape1='impulse'", "ae_shape1='none'").replace("e_impulse=1d-4", "e_impulse=0d0")


require_hybrid_route_contract(main_source)
divided_route = main_dft_source.split("subroutine run_dg_hybrid_divided_ground_state_for_main", 1)[1].split(
    "end subroutine run_dg_hybrid_divided_ground_state_for_main", 1
)[0]
divided_publisher = main_dft_source.split("subroutine publish_dg_hybrid_divided_v5", 1)[1].split(
    "end subroutine publish_dg_hybrid_divided_v5", 1
)[0]
assert "yn_dg_hybrid_divided_scf=='y'.or.yn_dg_hybrid_continuation_scf=='y'" in main_dft_source.lower()
nonlocal_builder = main_dft_source.lower().split(
    "subroutine assemble_dg_hybrid_divided_nonlocal_rows", 1
)[1].split("end subroutine assemble_dg_hybrid_divided_nonlocal_rows", 1)[0]
assert "complex(8),intent(out)::nonlocal_action(:,:)" in nonlocal_builder
assert "any(shape(nonlocal_action)/=[global_count,size(ow_core_ids)])" in nonlocal_builder
assert "allocate(dg_hybrid_interior_nonlocal_action(size(divided_effective_ids),size(ow_core_ids)))" in main_dft_source.lower()
assert "run_dg_hybrid_concrete_continuation" not in main_dft_source.lower()
assert "enddo terminal_lcfo_refinement" in divided_route.lower()
assert "complete_action_strength(q)*support_projector_values(position)*complete_overlap(basis,q)" in nonlocal_builder
assert divided_route.lower().count("call publish_dg_hybrid_divided_v5") == 1
assert divided_publisher.lower().count("call publish_rt_dg_hybrid_checkpoint_v5") == 1
assert "call write_rt_dg_hybrid_checkpoint_v5" not in divided_publisher.lower()
assert "solved_coefficients=final_solved_coefficients" in divided_route.lower()
assert "solved_eigenvalues=final_solved_eigenvalues" in divided_route.lower()
for old, replacement in (
    ("[HYBRID-RT-ROUTE]", "[REMOVED-HYBRID-RT-ROUTE]"),
    ("[HYBRID-RT-STEP]", "[REMOVED-HYBRID-RT-STEP]"),
    ("potential=PP+HARTREE+XC", "potential=CONVENTIONAL"),
    ("[HYBRID-RT-CHECKPOINT-STATIONARITY]", "[REMOVED-CHECKPOINT-STATIONARITY]"),
    ("[HYBRID-RT-REFRESH-STATIONARITY]", "[REMOVED-REFRESH-STATIONARITY]"),
    ("establish_fixed_density_reference=.true.", "establish_fixed_density_reference=.false."),
):
    mutated = main_source.replace(old, replacement)
    try:
        require_hybrid_route_contract(mutated)
    except AssertionError:
        pass
    else:
        raise AssertionError(f"route contract mutation survived: {old}")

env = os.environ.copy()
env["OMP_NUM_THREADS"] = "1"
env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
mpiexec = shutil.which("mpiexec")
assert mpiexec, "MPI launcher is required"
configured_build = os.environ.get("SALMON_DG_PRODUCTION_BUILD")
assert configured_build, "production smoke requires the current validated SALMON build"
salmon_build = Path(configured_build).resolve()
cache = (salmon_build / "CMakeCache.txt").read_text()
for option in ("USE_MPI:BOOL=ON", "USE_SCALAPACK:BOOL=ON", "USE_SPGLIB:BOOL=ON"):
    assert option in cache, option
assert "USE_EIGENEXA:BOOL=OFF" in cache
salmon = salmon_build / "salmon"
assert salmon.exists(), "current-source SALMON build did not produce the executable"
binary_hash = hashlib.sha256(salmon.read_bytes()).hexdigest()

with tempfile.TemporaryDirectory(prefix="hybrid-production-smoke-") as name:
    work = Path(name)
    shutil.copy2(root / "samples/exercise_01_C2H2_gs/H_rps.dat", work / "H_rps.dat")
    gs = subprocess.run(
        [mpiexec, "-n", "4", str(salmon)], input=h4_gs_input(), cwd=work, env=env,
        capture_output=True, text=True, timeout=180,
    )
    assert gs.returncode == 0, (gs.stdout, gs.stderr)
    assert gs.stdout.count("end SALMON") == 4, gs.stdout
    for receipt in (
        "[DG-DC-SEED]", "[DG-FRAGMENT-WF]", "[DG-HYBRID-DIVIDED-SEED]",
        "[OW-GS] fixed-density/non-self-consistent divided WF+PW LCFO solved once",
    ):
        assert receipt in gs.stdout, receipt
    handoff = re.findall(
        r"\[HYBRID-GS-HANDOFF\] route=divided-terminal-lcfo-v5 construction_rank=(\d+) solved_rank=(\d+) "
        r"certified_rank=(\d+) rt_rank=(\d+) occupied_rank=(\d+) projector_count=(\d+).*writer_count=(\d+)", gs.stdout,
    )
    assert len(handoff) == 1, gs.stdout
    construction, solved, certified, rt_rank, occupied, projector_count, writer_count = map(int, handoff[0])
    assert construction == solved == rt_rank and construction > certified >= occupied > 0 and writer_count == 1
    assert projector_count == 12, (
        "two-fragment overlapping buffers did not retain the 12-projector full H4 system identity",
        projector_count,
    )
    checkpoint = work / "hybrid_dg_ground_state.chk.manifest"
    shards = sorted(work.glob("hybrid_dg_ground_state.chk.v5.*.rank*.shard"))
    assert checkpoint.is_file() and checkpoint.stat().st_size > 64 and len(shards) == 4, (
        "formal divided terminal LCFO did not publish its v5 manifest and rank shards",
        gs.stdout[-6000:], gs.stderr,
    )
    assert checkpoint.read_bytes()[:32].rstrip(b" \0") == b"SALMON_HYBRID_DG_MANIFEST_V5"

    rt = subprocess.run(
        [mpiexec, "-n", "4", str(salmon)], input=h4_rt_input(), cwd=work, env=env,
        capture_output=True, text=True, timeout=180,
    )
    assert rt.returncode == 0, (gs.stdout, rt.stdout, rt.stderr)
    assert rt.stdout.count("[HYBRID-RT-ROUTE] propagator=EXP potential=PP+HARTREE+XC") == 1
    assert rt.stdout.count("[HYBRID-RT-HANDOFF]") == 1
    assert rt.stdout.count("[HYBRID-RT-POTENTIAL] update=0 path=PP+HARTREE+XC") == 1
    energy_identity = re.search(
        r"\[HYBRID-RT-ENERGY-IDENTITY\]\s+checkpoint=\s*(\S+)\s+refreshed=\s*(\S+)\s+defect=\s*(\S+)",
        rt.stdout,
    )
    assert energy_identity, "GS-to-RT physical energy identity was skipped"
    energy_values = [float(value) for value in energy_identity.groups()]
    assert all(math.isfinite(value) for value in energy_values)
    assert energy_values[2] <= 1e-10 * max(1.0, abs(energy_values[0]))
    state_match = re.search(
        r"\[HYBRID-RT-STATE\]\s+basis_norm=\s*([^ ]+)\s+density_norm=\s*([^ ]+)\s+"
        r"occupation_sum=\s*([^ ]+)\s+hamiltonian_norm=\s*([^ ]+)\s+kinetic_norm=\s*([^ ]+)\s+"
        r"nonlocal_norm=\s*([^ ]+)\s+local_norm=\s*([^ ]+)\s+sipg_norm=\s*([^ ]+)\s+global_nnz=(\d+)",
        rt.stdout,
    )
    assert state_match, rt.stdout
    state_values = [float(value) for value in state_match.groups()[:8]]
    assert all(math.isfinite(value) and value > 0.0 for value in state_values), state_values
    global_nnz = int(state_match.group(9))
    assert 0 < global_nnz < construction * construction, (
        "four-fragment DG operator unexpectedly became globally dense", global_nnz, construction,
    )
    assert global_nnz <= 3 * construction * math.ceil(construction / 4), (
        "DG CSR exceeds the local-plus-two-neighbor block bound", global_nnz, construction,
    )
    projection_match = re.search(
        r"\[HYBRID-RT-PROJECTION\]\s+local_potential_norm=\s*([^ ]+)\s+"
        r"hamiltonian_norm=\s*([^ ]+)\s+initial_delta=\s*([^ ]+)\s+path=PP\+HARTREE\+XC", rt.stdout,
    )
    assert projection_match, rt.stdout
    assert all(math.isfinite(float(value)) and float(value) > 0.0 for value in projection_match.groups())
    steps = re.findall(
        r"\[HYBRID-RT-STEP\] step=(\d+) metric_norm=\s*([^ ]+) orbital_energy=\s*([^ ]+) "
        r"polarization_norm=\s*([^ ]+) density_norm=\s*([^ ]+) update_count=(\d+) path=PP\+HARTREE\+XC",
        rt.stdout,
    )
    assert steps and [int(step[0]) for step in steps] == [1, 2], rt.stdout
    assert all(math.isfinite(float(value)) for step in steps for value in step[1:5])
    assert all(float(step[1]) > 0.0 and float(step[4]) > 0.0 for step in steps)
    assert "[DG-OW-RT]" not in rt.stdout
    assert "production DG requires a build with MPI and ScaLAPACK support" not in rt.stdout + rt.stderr

    zero = subprocess.run(
        [mpiexec, "-n", "4", str(salmon)], input=h4_rt_zero_input(), cwd=work, env=env,
        capture_output=True, text=True, timeout=180,
    )
    assert zero.returncode == 0, (zero.stdout, zero.stderr)
    checkpoint_stationarity = re.search(
        r"\[HYBRID-RT-CHECKPOINT-STATIONARITY\]\s+orbital_residual=\s*(\S+)\s+metric_defect=\s*(\S+)",
        zero.stdout,
    )
    assert checkpoint_stationarity, zero.stdout
    assert all(float(value) <= 1e-4 for value in checkpoint_stationarity.groups()), checkpoint_stationarity.groups()
    refresh_stationarity = re.search(
        r"\[HYBRID-RT-REFRESH-STATIONARITY\]\s+h_residual=\s*(\S+)\s+hamiltonian_delta=\s*(\S+)",
        zero.stdout,
    )
    assert refresh_stationarity, zero.stdout
    assert float(refresh_stationarity.group(1)) <= 1e-4, refresh_stationarity.groups()
    zero_steps = re.findall(
        r"\[HYBRID-RT-STATIONARITY\] step=(\d+) density=\s*(\S+) energy=\s*(\S+) "
        r"projector=\s*(\S+) electron=\s*(\S+) h_residual=\s*(\S+) current_bound=\s*(\S+)", zero.stdout,
    )
    assert [int(step[0]) for step in zero_steps] == [1, 2], zero.stdout
    assert all(math.isfinite(float(value)) and float(value) <= 1e-4 for step in zero_steps for value in step[1:])
    checkpoint_hash = hashlib.sha256(checkpoint.read_bytes() + b"".join(s.read_bytes() for s in shards)).hexdigest()
    print(f"producer_binary_sha256={binary_hash}")
    print(f"checkpoint_sha256={checkpoint_hash}")
    for line in rt.stdout.splitlines():
        if line.startswith("[HYBRID-RT-"):
            print(f"ranks=4 {line}")
    for line in zero.stdout.splitlines():
        if line.startswith("[HYBRID-RT-STATIONARITY]") or line.startswith("[HYBRID-RT-REFRESH-STATIONARITY]"):
            print(f"ranks=4 zero-field {line}")

print("PASS actual production divided H4 GS-to-v5-to-Exp-RT smoke")
