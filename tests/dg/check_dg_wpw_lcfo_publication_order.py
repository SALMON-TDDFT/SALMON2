#!/usr/bin/env python3
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
src=(ROOT/'src/gs/dc/lcfo_flux.f90').read_text().lower()
calc=src.index('call calc_hamiltonian_matrix')
imp=src.index('call import_wpw_lcfo_ww_components',calc)
build=src.index('call build_dg_wpw_production_operator',calc)
bind=src.index('call bind_dg_wpw_hs_callbacks',calc)
assert calc < imp < build < bind
assert 'call build_dg_wpw_production_operator' not in src[:calc]
assert 'call bind_dg_wpw_hs_callbacks' not in src[:calc]
release=src.index('call release_dg_wpw_bounded_operator')
comm_free=src.index('call mpi_comm_free(wpw_production_comm')
assert bind < release < comm_free
print('PASS LCFO WW import precedes atomic production operator publication')
