#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
layout = (ROOT / 'src/rt/dg/rt_dg_wpw_column_layout.f90').read_text().lower()
context = (ROOT / 'src/common/dg_wpw_production_context.f90').read_text().lower()
lcfo = (ROOT / 'src/gs/dc/lcfo_flux.f90').read_text().lower()
builder = (ROOT / 'src/rt/dg/rt_dg_wpw_production_builder.f90').read_text().lower()
sparse = (ROOT / 'src/rt/dg/rt_dg_wpw_sparse_builder.f90').read_text().lower()

for token in (
    'fragment_root',
    'initialize_wpw_fragment_root_layout',
    'wpw_fragment_root_owner',
    'initialize_dg_wpw_fragment_root_context',
    'owned_fragment_id',
):
    if token not in layout + context:
        raise SystemExit(f'missing fragment-root ownership contract: {token}')

for forbidden in ('global_owner', 'owner_by_fragment', 'fragment_owner_ids'):
    if forbidden in layout + context:
        raise SystemExit(f'forbidden O(N) owner metadata: {forbidden}')

print('PASS fragment-root ownership source contract')

for token in (
    'mpi_comm_split',
    'mpi_undefined',
    'initialize_dg_wpw_fragment_root_context',
):
    if token not in lcfo:
        raise SystemExit(f'lcfo_flux does not create the fragment-root production communicator: {token}')

if 'initialize_dg_wpw_production_context(wpw_context' in lcfo:
    raise SystemExit('lcfo_flux still initializes the arithmetic all-rank production context')
if 'mpi_allreduce' not in lcfo:
    raise SystemExit('fragment-root initialization failure is not reduced over all ranks')
if "if(wpw_mode_info/=0)stop 'dg wpw fragment-root context initialization failed'" in lcfo:
    raise SystemExit('fragment root can stop before all-rank initialization failure reduction')

print('PASS lcfo_flux fragment-root communicator source contract')

if 'initialize_wpw_fragment_root_layout' not in builder:
    raise SystemExit('production builder recreates no fragment-root column layout')
if 'wpw_fragment_root_owner' not in sparse:
    raise SystemExit('sparse builder does not validate fragment-root candidate ownership')
print('PASS fragment-root sparse builder source contract')
