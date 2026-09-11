"""Do not retain unreachable DG response construction inside generic output."""
from pathlib import Path

root = Path(__file__).resolve().parents[2]
writer = (root/'src/io/write.f90').read_text().lower()
rt = (root/'src/rt/main_tddft.f90').read_text().lower()
for name in ('dg_polarization_history', 'write_dg_polarization_data',
             'write_dg_polarization_response_3d'):
    assert name not in writer+rt, f'unreachable DG response processing remains: {name}'
# Preserve the current v5 compute-then-diagnostic route and conventional output.
hybrid_rt = rt.split('subroutine run_dg_hybrid_continuation_rt', 1)[1].split(
    'end subroutine run_dg_hybrid_continuation_rt', 1)[0]
tokens = ('call propagate_rt_dg_hybrid_length_gauge',
          'call reconstruct_rt_dg_hybrid_density',
          'call update_rt_dg_hybrid_density', "'[hybrid-rt-step] step='")
# The initial density update precedes propagation; inspect the time-step loop.
step_loop = hybrid_rt.split('do step=1,nt', 1)[1]
def ordered_output(body):
    positions = [body.find(token) for token in tokens]
    return all(position >= 0 for position in positions) and positions == sorted(positions)
assert ordered_output(step_loop)
for token in tokens:
    assert not ordered_output(step_loop.replace(token, 'omitted', 1)), token
for name in ('write_response_0d', 'write_response_3d', 'write_pulse_0d', 'write_pulse_3d'):
    assert f'subroutine {name}(' in writer
print('DG write output scope: PASS')
