"""Do not retain unreachable DG response construction inside generic output."""
from pathlib import Path

root = Path(__file__).resolve().parents[2]
writer = (root/'src/io/write.f90').read_text().lower()
rt = (root/'src/rt/main_tddft.f90').read_text().lower()
for name in ('dg_polarization_history', 'write_dg_polarization_data',
             'write_dg_polarization_response_3d'):
    assert name not in writer+rt, f'unreachable DG response processing remains: {name}'
# Preserve the used compute-then-output observable route and conventional output.
assert 'call evaluate_dg_overlapping_wannier_observables' in rt
assert 'call write_dg_overlapping_wannier_rt_observable_sample' in rt
for name in ('write_response_0d', 'write_response_3d', 'write_pulse_0d', 'write_pulse_3d'):
    assert f'subroutine {name}(' in writer
print('DG write output scope: PASS')
