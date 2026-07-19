#!/usr/bin/env python3
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
text=(ROOT/'src/gs/dc/dg_wpw_rank_local_quadrature.f90').read_text().lower()
start=text.index('subroutine add_dg_wpw_core_point')
end=text.index('end subroutine',start)
body=text[start:end]
assert 'allocate(' not in body, 'streaming core-point accumulation allocates per point'
assert 'wpw_volume_weak_pair' in body, 'streaming accumulator does not evaluate weak pairs directly'
print('PASS core-point accumulator reuses persistent storage')
