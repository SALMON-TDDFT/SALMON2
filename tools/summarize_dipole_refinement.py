"""Summarize fixed study cases without executing or mutating a solver.

Usage: python tools/summarize_dipole_refinement.py NEW_JSON COARSE TIME SPACE FINE FINE_TIME
"""
import json
from pathlib import Path
import sys
import numpy as np

out = Path(sys.argv[1])
labels = ('coarse','time','space','fine','fine-time')
if len(sys.argv) != 7:
    raise ValueError('new output JSON and five study directories required')
cases, arrays = {}, {}
for label, directory in zip(labels,sys.argv[2:]):
    root = Path(directory)
    data = json.loads((root/'report.json').read_text())
    with np.load(root/'comparison.npz',allow_pickle=False) as a:
        arrays[label] = {key:a[key].copy() for key in ('time','expected_Ez','actual_Ez','E')}
        data['self_E_max'] = float(np.max(np.abs(a['E'][:,0,2])))
    sample = arrays[label]
    if not all(np.isfinite(v).all() for v in sample.values()):
        raise ValueError('nonfinite reference or measured arrays')
    norm = np.linalg.norm(sample['expected_Ez'])
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError('positive finite reference norm required')
    recomputed = np.linalg.norm(sample['actual_Ez']-sample['expected_Ez'])/norm
    np.testing.assert_allclose(recomputed,data['relative_L2_E'],rtol=1.e-13,atol=1.e-15)
    cases[label] = data
for label in labels:
    for key in ('source_position','observation_position','pulse_start','pulse_duration','observation_end'):
        if cases[label][key] != cases['coarse'][key]:
            raise ValueError('physical configuration changed: '+key)
spatial = ('time','space','fine')
if len({cases[label]['dt'] for label in spatial}) != 1:
    raise ValueError('spatial comparison requires identical dt')
if not cases['time']['dx'] > cases['space']['dx'] > cases['fine']['dx']:
    raise ValueError('spatial comparison requires strictly decreasing dx')
for label in spatial[1:]:
    np.testing.assert_array_equal(arrays[label]['time'],arrays['time']['time'])
    np.testing.assert_array_equal(arrays[label]['expected_Ez'],arrays['time']['expected_Ez'])
for key in ('relative_L2_E','source_work_relative_error'):
    if not cases['fine'][key] < cases['space'][key] < cases['time'][key]:
        raise AssertionError('spatial refinement did not improve '+key)
changes = {}
for fine,coarse in (('time','coarse'),('fine-time','space')):
    a,b = arrays[fine],arrays[coarse]
    if cases[fine]['dx'] != cases[coarse]['dx']:
        raise ValueError('temporal comparison must retain dx')
    np.testing.assert_allclose(a['time'][::2],b['time'],rtol=0,atol=1.e-14)
    changes[fine] = float(np.linalg.norm(a['actual_Ez'][::2]-b['actual_Ez'])/np.linalg.norm(b['expected_Ez']))
summary = dict(cases=cases,temporal_half_step_field_change=changes,
               scope='finite-window point-dipole refinement; source work is not measured surface flux',
               self_field_feedback_validated=False)
with out.open('x') as stream:
    json.dump(summary,stream,indent=2)
    stream.write('\n')
print('spatial refinement and temporal comparisons recorded:',out)
