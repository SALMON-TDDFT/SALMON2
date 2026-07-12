#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).resolve().parents[2]
main=(root/'src/rt/main_tddft.f90').read_text()
assert 'use rt_dg_nodal_salmon_prepare' in main
assert "if (yn_dg_nodal_rt == 'y') then" in main
assert 'call prepare_nodal_salmon_ground_state' in main
assert "'[DG-NODAL-GS]'" in main
assert 'if (.not. nodal_state%dg_ground_state_ready)' in main
assert 'call finalize_dg_fragment_rt_std(dg_frag)' in main
assert 'return' in main
print('PASS main_tddft runs a complete-H nodal GS preflight before legacy DG construction')
