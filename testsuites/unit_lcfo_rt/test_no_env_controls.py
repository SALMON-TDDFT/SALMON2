"""Guard the user contract: custom HSE/LCFO configuration comes from input files."""
from pathlib import Path
repo=Path(__file__).resolve().parents[2]
files=['src/rt/lcfo_rt_basis.f90','src/xc/exx_lcfo_rt.f90','src/xc/lcfo_rt_wannier.f90',
       'src/xc/lcfo_seed.f90','src/xc/exx_k_exchange.f90','src/xc/exx_native.f90',
       'src/common/hse_reference_export.f90']
for name in files:
    assert 'get_environment_variable' not in (repo/name).read_text().lower(), name
print('HSE/LCFO environment control reads removed')
