#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).resolve().parents[2]
glob=(root/'src/io/salmon_global.f90').read_text()
inp=(root/'src/io/inputoutput.f90').read_text()
for name in ['yn_dg_nodal_rt','dg_nodal_gs_relax_step','dg_nodal_gs_max_iter','dg_nodal_gs_tol','dg_nodal_taylor_order']:
    assert name in glob, f'missing global {name}'
    assert name in inp, f'missing input wiring {name}'
assert "if(yn_dg_nodal_rt=='y' .and. yn_dg_fragment_rt/='y')" in inp
assert "time_integrator_dg_fragment/='taylor4pc'" in inp
assert 'dt > 0.02d0' in inp
assert 'call yn_argument_check(yn_dg_nodal_rt)' in inp
print('PASS nodal real-space controls are namelist-driven and enforce dt<=0.02 Taylor')
