#!/usr/bin/env python3
from pathlib import Path
import re
ROOT=Path(__file__).resolve().parents[2]
text=(ROOT/'src/gs/dc/lcfo_flux.f90').read_text().lower()
start=text.index("if(yn_dg_wpw_production=='y') then")
end=text.index('    endif\n    call check_lcfo_basis_potential_inputs_finite',start)
block=text[start:end]
assert "stop 'dg wpw" not in block and 'stop "dg wpw' not in block, \
  'WPW production block contains rank-local STOP'
assert 'subroutine wpw_collective_require' in text, 'missing collective WPW failure gate'
gate=text[text.index('subroutine wpw_collective_require'):]
assert re.search(r'call\s+mpi_allreduce',gate), 'collective gate does not reduce local failures'
assert re.search(r'call\s+mpi_abort',gate), 'collective gate is not communicator-fatal'
print('PASS WPW production failures are communicator-fatal')
