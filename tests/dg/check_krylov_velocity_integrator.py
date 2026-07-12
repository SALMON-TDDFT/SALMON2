#!/usr/bin/env python3
from pathlib import Path

root=Path(__file__).resolve().parents[2]
fragment=(root/'src/rt/dg/rt_dg_fragment.f90').read_text()
iteration=(root/'src/rt/dg/rt_dg_iteration.f90').read_text()
inputs=(root/'src/io/inputoutput.f90').read_text()
krylov=(root/'src/rt/dg/rt_dg_integrator_krylov.f90')

assert krylov.exists(), 'missing Krylov integrator implementation'
body=krylov.read_text()
assert 'subroutine time_evolution_krylov' in body
assert 'call calculate_time_derivative' in body
assert 'deriv(nlocal,nstate_prop,1)' in body
assert 'call calculate_time_derivative(dg_frag,system,mg,ppg,Ac_mid,deriv,1,nstate_prop)' in body
assert 'subroutine global_dots' in body
assert 'logical,allocatable :: active(:)' in body
assert 'call eigen_zheev' in body
assert "case('krylov')" in fragment
assert 'case(6)' in iteration and 'call time_evolution_krylov' in iteration
assert "time_integrator_dg_fragment/='krylov'" in inputs
print('PASS block-sparse velocity-gauge Krylov integrator is wired')
