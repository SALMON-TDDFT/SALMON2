#!/usr/bin/env python3
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
src=(ROOT/'src/rt/dg/rt_dg_plane_wave.f90').read_text().lower()
start=src.index('subroutine init_plane_wave_basis')
end=src.index('end subroutine init_plane_wave_basis',start)
body=src[start:end]
assert 'use dg_wpw_g_modes' in body, 'RT does not use the GS-neutral G-mode selector'
assert 'select_dg_wpw_g_modes' in body, 'RT does not call the GS-neutral G-mode selector'
for duplicate in ('shell_keys','n_pw_candidate','cubic_k_order_less(k_indices'):
  assert duplicate not in body, f'RT retains duplicate G-mode selection: {duplicate}'
print('PASS GS and RT share deterministic WPW G-mode selection')
