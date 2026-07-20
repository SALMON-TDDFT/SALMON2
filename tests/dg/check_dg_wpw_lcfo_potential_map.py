#!/usr/bin/env python3
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
map_src=ROOT/'src/gs/dc/dg_wpw_lcfo_potential_map.f90'
volume_src=ROOT/'src/gs/dc/dg_wpw_lcfo_volume_operator.f90'
assert map_src.exists(), 'missing transactional LCFO potential map'
assert volume_src.exists(), 'missing shared LCFO volume rebuild kernel'
text=map_src.read_text().lower();volume=volume_src.read_text().lower()
lcfo=(ROOT/'src/gs/dc/lcfo_flux.f90').read_text().lower()
for token in (
  'type,public::s_dg_wpw_map_candidate',
  'rho_in','rho_raw','rho_candidate','potential_candidate','energy_candidate',
  'operator_candidate','publication_epoch','mixer_history','residual_history',
  'subroutine run_dg_wpw_lcfo_potential_map',
  'subroutine rollback_dg_wpw_lcfo_potential_map',
  'subroutine publish_dg_wpw_lcfo_potential_map',
  'density_mix_only',
): assert token in text, f'missing transactional map contract: {token}'
for token in ('rebuild_dg_wpw_ww_volume','rebuild_dg_wpw_wp_pp_volume','fixed_kinetic','fixed_nonlocal','fixed_face'):
  assert token in volume, f'missing volume rebuild contract: {token}'
for forbidden in ('global_h(', 'global_s(', 'all_fragment_owner', 'do ifrag=1,n_frag',
                  'mpi_allgather', 'mpi_allgatherv', 'rho_global(:)'):
  assert forbidden not in text+volume, f'forbidden production construct: {forbidden}'
assert 'grid_id=wpw_core_global_grid_id(grid)' in lcfo, \
  'production quadrature does not retain its established zero-origin global grid ID convention'
assert 'global_grid_id=wpw_core_global_grid_id(global_point-1)' in lcfo, \
  'Wannier tail points are not converted from jxyz_tot to the core accumulator grid ID convention'
assert 'point_capacity=wpw_local_core_point_count' in lcfo, \
  'production quadrature does not preallocate bounded local core-point storage'
print('PASS transactional LCFO density/Hartree/LDA potential-map source contract')
