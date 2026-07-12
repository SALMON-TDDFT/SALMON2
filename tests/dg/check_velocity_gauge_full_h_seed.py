#!/usr/bin/env python3
from pathlib import Path

source=(Path(__file__).resolve().parents[2]/'src/rt/main_tddft.f90').read_text()
assert "yn_dg_full_h_eigen_seed == 'y' .and. yn_dg_length_gauge == 'n'" in source
assert "'[DG-FULL-H-SEED-HRES-VG]'" in source
assert "velocity-gauge Full DG seed is not stationary" in source
print('PASS velocity-gauge propagation receives a stationary Full-DG seed')
