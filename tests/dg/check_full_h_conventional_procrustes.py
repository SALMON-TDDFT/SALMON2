#!/usr/bin/env python3
from pathlib import Path

root=Path(__file__).resolve().parents[2]
source=(root/'src/rt/dg/rt_dg_fragment.f90').read_text()
assert 'subroutine build_full_h_conventional_alignment' in source
assert 'qblock=matmul(oblock,metric_invhalf)' in source
assert "'[DG-FULL-H-PROCRUSTES] block='" in source
assert 'full_h_conventional_alignment(1:nb,1:nb,1)' in source
print('PASS Full-H/conventional block Procrustes alignment is wired')
