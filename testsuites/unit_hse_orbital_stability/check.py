#!/usr/bin/env python3
"""Check a native frozen-H solve diagnostic; never generate or collect GS data."""
import argparse
import math
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument('log', type=Path)
a = p.parse_args()
s = a.log.read_text()
rows = [list(map(float, line.split(':', 1)[1].split())) for line in s.splitlines()
        if line.startswith('HSE_STATIONARITY gram_error/')]
if len(rows) != 3 or 'end SALMON' not in s:
    raise SystemExit('FAIL incomplete fixed-H diagnostic')
if any(len(r) != 5 or not all(map(math.isfinite, r)) for r in rows):
    raise SystemExit('FAIL invalid diagnostic metrics')
if any(r[0] > 1e-10 or r[1] > 1e-10 for r in rows):
    raise SystemExit('FAIL orthogonality or projected Hermiticity')
if rows[0][2] <= 1e-6 or rows[1][2] >= 1e-6 or rows[2][2] >= 1e-6:
    raise SystemExit(f'FAIL stationary occupied density not recovered: {[r[2] for r in rows]}')
print('PASS frozen-H orbital stability:', [r[2] for r in rows])
