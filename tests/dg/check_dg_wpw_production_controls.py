#!/usr/bin/env python3
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
global_src=(ROOT/'src/io/salmon_global.f90').read_text().lower()
io_src=(ROOT/'src/io/inputoutput.f90').read_text().lower()
lcfo=(ROOT/'src/gs/dc/lcfo_flux.f90').read_text().lower()
rt_window=(ROOT/'src/rt/dg/rt_dg_wpw_window.f90').read_text().lower()

for name in ('dg_wpw_window_buffer','dg_wpw_window_width'):
    if name not in global_src:
        raise SystemExit(f'missing production WPW namelist state: {name}')
    if io_src.count(name)<5:
        raise SystemExit(f'production WPW control is not defaulted/read/broadcast/logged/validated: {name}')
    if name not in lcfo:
        raise SystemExit(f'lcfo_flux does not consume production WPW control: {name}')

print('PASS production WPW window controls are explicit namelist state')

if 'get_environment_variable' in rt_window:
    raise SystemExit('RT WPW window still has result-changing environment-only controls')
for name in ('dg_wpw_window_buffer','dg_wpw_window_width'):
    if name not in rt_window:
        raise SystemExit(f'RT WPW window does not consume shared namelist state: {name}')
print('PASS GS and RT reject environment-only WPW window controls')

for token in ('build_wpw_support_geometry','build_dg_wpw_fragment_support','wpw_support_fragment_ids',
              'dg_wpw_support_overlap_box','wpw_support_records'):
    if token not in lcfo:
        raise SystemExit(f'lcfo_flux does not prepare bounded WPW support geometry: {token}')
print('PASS lcfo_flux prepares bounded structured WPW support geometry')
