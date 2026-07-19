#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
source = (ROOT / 'src/rt/dg/rt_dg_wpw_window.f90').read_text()

required = [
    'use dg_wpw_windows, only: evaluate_dg_wpw_normalized_windows',
    'call evaluate_dg_wpw_normalized_windows',
]
for token in required:
    if token not in source:
        raise SystemExit(f'RT WPW window does not share the production evaluator: {token}')

print('PASS RT and production share the normalized WPW window evaluator')
