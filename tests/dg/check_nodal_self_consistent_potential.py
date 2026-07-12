#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
source = root / 'src/rt/dg/rt_dg_nodal_density.f90'
assert source.exists(), 'missing nodal all-system density reconstruction'
body = source.read_text()
assert 'subroutine reconstruct_nodal_density_mpi' in body
assert 'system%rocc' in body
assert '/system%hvol' in body.replace(' ', '')
assert 'call MPI_Allreduce' in body
assert 'call sym_rho' in body

main = (root / 'src/rt/main_tddft.f90').read_text()
assert 'use rt_dg_nodal_density' in main
assert 'call reconstruct_nodal_density_mpi' in main
assert 'call hartree(' in main
assert 'call exchange_correlation(' in main
assert 'call build_nodal_local_potential' in main
propagate = main.index('call propagate_nodal_salmon_taylor')
density = main.index('call reconstruct_nodal_density_mpi', propagate)
potential = main.index('call build_nodal_local_potential', density)
assert propagate < density < potential
print('PASS nodal Taylor RT rebuilds rho, Vh, and Vxc every step')
