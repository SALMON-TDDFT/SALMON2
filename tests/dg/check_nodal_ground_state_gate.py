#!/usr/bin/env python3
from pathlib import Path

root=Path(__file__).resolve().parents[2]
types=(root/'src/rt/dg/rt_dg_nodal_types.f90').read_text()
taylor=(root/'src/rt/dg/rt_dg_nodal_taylor.f90').read_text()
assert 'logical :: dg_ground_state_ready = .false.' in types
assert 'real(8) :: dg_ground_state_residual' in types
assert 'subroutine accept_nodal_dg_ground_state' in types
assert "if (.not. state%dg_ground_state_ready) stop" in taylor
assert 'state%dg_ground_state_residual' in taylor
print('PASS nodal RT requires a verified eigenstate of its own DG Hamiltonian')
