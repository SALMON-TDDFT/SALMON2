#!/usr/bin/env python3
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
main=(ROOT/'src/rt/main_tddft.f90').read_text().lower()
inp=(ROOT/'src/io/inputoutput.f90').read_text().lower()
glob=(ROOT/'src/io/salmon_global.f90').read_text().lower()
for token in ('yn_dg_wpw_checkpoint_rt','dg_wpw_checkpoint_manifest','dg_wpw_checkpoint_rank_prefix',
              'dg_wpw_checkpoint_identity_tolerance'):
    assert token in inp and token in glob, f'missing explicit RT checkpoint control {token}'
assert 'load_rt_dg_wpw_checkpoint_handoff' in main, 'RT entry point does not load WPW checkpoint'
assert "yn_restart == 'y' .and. yn_dg_wpw_checkpoint_rt /= 'y'" in main, \
    'checkpoint-backed RT can still project conventional restart orbitals'
assert 'wpw_checkpoint_handoff%valid' in main, 'RT entry point does not require validated handoff'
print('PASS RT entry point uses checkpoint identity route without conventional projection')
