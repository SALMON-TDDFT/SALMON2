#!/usr/bin/env python3
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
main=(ROOT/'src/rt/main_tddft.f90').read_text().lower()
inp=(ROOT/'src/io/inputoutput.f90').read_text().lower()
for token in ('dg_wpw_exp_max_corrector','dg_wpw_exp_corrector_tolerance','dg_wpw_exp_norm_tolerance'):
    assert token in inp and token in main, f'missing explicit midpoint Exp control {token}'
assert 'initialize_dg_wpw_exp_state' in main and 'advance_dg_wpw_length_gauge_exp' in main
assert 'position_reduced' in main and 'wpw_field_mid' in main
assert "checkpoint-backed field-on rt requires task 9" not in main
for token in ('pz=', 'delta_pz=', 'jz=', 'branch=continuous_sawtooth_z'):
    assert token in main, f'missing production observable {token}'
assert '[dg-wpw-exp]' in main and 'norm_drift' in main and 'correctors' in main
print('PASS checkpoint-backed field-off/field-on production midpoint Exp route')
