"""The v5 endpoint and auxiliary occupied stream have distinct owners."""
from pathlib import Path

root = Path(__file__).resolve().parents[2]
old = root / 'src/rt/dg/rt_dg_hybrid_checkpoint.f90'
assert not old.exists(), 'ambiguous checkpoint module remains'
v5 = (root / 'src/rt/dg/rt_dg_hybrid_checkpoint_v5.f90').read_text().lower()
occupied = (root / 'src/rt/dg/rt_dg_hybrid_occupied_checkpoint.f90').read_text().lower()
for name in ('publish_rt_dg_hybrid_checkpoint_v5',
             'collective_rt_dg_hybrid_publication_precondition',
             'collective_rt_dg_hybrid_publication_mapping_precondition'):
    assert f'subroutine {name}(' in v5
    assert name not in occupied
assert 's_rt_dg_hybrid_v5_publication_authorization' in v5
assert 'use rt_dg_hybrid_checkpoint_v5' not in occupied
assert 'write_rt_dg_hybrid_occupied_checkpoint' in occupied
assert 'read_rt_dg_hybrid_occupied_checkpoint' in occupied
assert 'occupied_version=2' in occupied
assert 'schema_version=5' in v5
cmake = (root / 'src/rt/CMakeLists.txt').read_text()
assert 'dg/rt_dg_hybrid_checkpoint.f90' not in cmake
assert 'dg/rt_dg_hybrid_occupied_checkpoint.f90' in cmake
print('checkpoint responsibility separation: PASS')
