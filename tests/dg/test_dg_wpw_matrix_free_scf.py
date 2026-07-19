#!/usr/bin/env python3
from pathlib import Path
import os
import subprocess

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "src/gs/dc/dg_wpw_matrix_free_scf.f90"
MAIN = ROOT / "src/gs/main_dft.f90"
LCFO = ROOT / "src/gs/dc/lcfo_flux.f90"
GLOBAL = ROOT / "src/io/salmon_global.f90"
INPUT = ROOT / "src/io/inputoutput.f90"

assert SOURCE.exists(), "missing bounded matrix-free DG-DC consumer"
source = SOURCE.read_text().lower()
main = MAIN.read_text().lower()
lcfo = LCFO.read_text().lower()
controls = (GLOBAL.read_text() + INPUT.read_text()).lower()

for token in (
    "run_dg_wpw_matrix_free_scf",
    "apply_h_batch",
    "apply_s_batch",
    "global_gram_batch",
    "n_occ + extra_states",
    "density_residual",
    "potential_residual",
    "energy_residual",
    "projector_residual",
    "generalized_residual",
    "metric_orthonormality",
    "charge_error",
    "gap",
    "fixed_point",
):
    assert token in source, f"missing matrix-free SCF contract: {token}"

for forbidden in (
    "allocate(h_global",
    "allocate(s_global",
    "allocate(density_matrix",
    "run_dg_wpw_fixed_scf",
):
    assert forbidden not in source + main, f"dense production path detected: {forbidden}"

for stage in (
    "initial_apply_begin", "initial_apply_end", "initial_eigh_end",
    "expanded_apply_end", "reduced_eigh_end", "inner_end",
):
    assert f"call trace_solve_window('{stage}'" in source, \
        f"missing solve-window boundary diagnostic: {stage}"
assert "dg_wpw_preconditioner" in source and "present(precondition)" in source, \
    "matrix-free solver must accept an optional bounded preconditioner callback"

for token in (
    "yn_dg_wpw_production",
    "dg_wpw_extra_states",
    "dg_wpw_gap_threshold",
    "dg_wpw_metric_cutoff",
    "dg_wpw_scf_mix",
    "dg_wpw_scf_max_iter",
    "dg_wpw_scf_residual_tolerance",
):
    assert token in controls, f"missing DG-DC namelist control: {token}"

mod = Path("/tmp/test_dg_wpw_matrix_free_scf_mod")
mod.mkdir(exist_ok=True)
exe = Path("/tmp/test_dg_wpw_matrix_free_scf_mpi")
fc = os.environ.get("MPIFC", "mpifort")
subprocess.run([fc,"-J",str(mod),"-I",str(mod),
    str(ROOT/"src/common/dg_generalized_algebra.f90"),
    str(ROOT/"src/common/dg_wpw_matrix_free_operator.f90"),str(SOURCE),
    str(ROOT/"tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90"),"-llapack","-lblas","-o",str(exe)],check=True,cwd="/tmp")
subprocess.run(["mpirun","-np","2",str(exe)],check=True,cwd=ROOT)
assert "call dc_lcfo_flux" in main, "main_dft does not reach the LCFO production owner"
assert "call run_wpw_production_scf" in lcfo, "lcfo_flux does not reach the stepwise production consumer"
print("PASS bounded matrix-free DG-DC production route contract")
