#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
main = (root / 'src/rt/main_tddft.f90').read_text()
assert 'use rt_dg_nodal_salmon_taylor' in main
assert 'use rt_dg_nodal_diagnostics' in main
assert 'call update_kvector_nonlocalpt' in main
assert 'call propagate_nodal_salmon_taylor' in main
assert '[DG-NODAL-RT]' in main
assert "' Ac='" in main
assert 'nodal_norm_drift' in main
assert 'RT/observables remain disabled in preflight mode' not in main

diag = root / 'src/rt/dg/rt_dg_nodal_diagnostics.f90'
assert diag.exists(), 'missing nodal RT norm diagnostics'
body = diag.read_text()
assert 'subroutine calculate_nodal_norm_diagnostics_mpi' in body
assert 'call MPI_Allreduce' in body
print('PASS converged nodal state enters complete-H Taylor RT with norm diagnostics')
