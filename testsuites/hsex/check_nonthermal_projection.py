#!/usr/bin/env python3
"""Offline assertions for the small-Si adaptive TD-projection campaign.

No SALMON computation or production W orchestration is performed in Python.
Keep wavefunctions, GS data and checkpoints on scratch. See the developer
manual for the reference/captured/initial5/cap/grown/restart input variants.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import struct
import numpy as np
from test_internal_si import compare_w, compare_current, read_w, completed_log


def captures(directory):
    text = (directory/'run.log').read_text()
    rows = re.findall(r'HSEX_PROJECTION_CAPTURE states=(\d+) rt_occupied=(\d+) missing_norm=\s*(\S+)', text)
    if not rows:
        raise AssertionError(f'Capture diagnostics missing: {directory}')
    result = [(int(n), int(rt), float(loss)) for n, rt, loss in rows]
    if any(rt != 4 or not np.isfinite(loss) for _, rt, loss in result):
        raise AssertionError('Propagated count changed or nonfinite capture')
    return result


def checkpoint_count(path):
    data = path.read_bytes()
    # Two scalar sequential-unformatted records: nk, no (4-byte GNU markers).
    if struct.unpack_from('=i', data)[0] != 4 or struct.unpack_from('=i', data, 12)[0] != 4:
        raise AssertionError('Unexpected checkpoint record markers')
    return struct.unpack_from('=i', data, 16)[0]


def check(root):
    result = {}
    for name in ('reference', 'captured', 'initial5', 'grown', 'restart'):
        completed_log(root/name/'run.log', require_team=name != 'restart')
    result['unchanged_W'] = compare_w(root/'captured/hsex-updated-000000001.bin',
                                     root/'reference/hsex-updated-000000001.bin')
    result['unchanged_current'] = compare_current(root/'captured/Si_rt.data',
                                                 root/'reference/Si_rt.data', [1, 2, 3])
    result['restart_current'] = compare_current(root/'restart/Si_rt.data',
                                               root/'grown/Si_rt.data', [2, 3])
    initial, grown = captures(root/'initial5'), captures(root/'grown')
    if grown[0][0] != 5 or grown[0][2] <= 1e-10 or grown[-1][0] <= 5 or grown[-1][2] > 1e-10:
        raise AssertionError('Forced expansion did not satisfy capture')
    if 'HSEX_PROJECTION_GROW old=5 new=10 rt_occupied=4' not in (root/'grown/run.log').read_text():
        raise AssertionError('Expected auxiliary-only growth missing')
    for name in ('captured', 'initial5', 'grown'):
        text = (root/name/'run.log').read_text()
        endpoint = np.array([float(v) for v in re.search(r'HSEX_PROJECTION_A endpoint:\s*([^\n]+)', text)[1].split()])
        trajectory = np.loadtxt(root/name/'Si_rt.data')
        np.testing.assert_allclose(endpoint, trajectory[0, 7:10], rtol=0, atol=1e-14)
        if np.linalg.norm(endpoint) == 0:
            raise AssertionError('Test does not exercise nonzero endpoint A')
        if checkpoint_count(root/name/'checkpoint_rt_000001/info.bin') != 4:
            raise AssertionError('RT checkpoint contains auxiliary propagated states')
    header, values, tail = read_w(root/'grown/hsex-updated-000000001.bin')
    step, nb = struct.unpack_from('=2i', tail)
    nk = int(np.prod(np.frombuffer(header, '=i4', 9, 4)[3:6]))
    if step != 1 or nb != grown[-1][0] or len(tail) != 8 + 16*nb*nk:
        raise AssertionError('Expanded W snapshot metadata invalid')
    pop = np.frombuffer(tail, '=f8', count=nb*nk, offset=8+8*nb*nk).reshape(nk, nb)
    if not np.isfinite(values).all() or not np.isfinite(pop).all():
        raise AssertionError('Nonfinite W or occupations')
    # Compare raw coefficients in W to the projection output, rather than
    # merely testing that a normalized trace equals the electron count.
    lines = (root/'grown/Si_ovlp.data').read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == '1')
    raw = np.array([[float(x) for x in line.split()[1:]] for line in lines[start+1:start+1+nk]])
    np.testing.assert_allclose(pop, raw, rtol=0, atol=2e-14)
    if not np.any((pop > 0) & (pop < 1e-7)):
        raise AssertionError('Projected nonthermal small populations missing')
    cap = root/'cap'
    if 'HSEX_PROJECTION_REJECT cap=5' not in (cap/'run.log').read_text():
        raise AssertionError('Collective cap rejection missing')
    if int((cap/'exit-status').read_text()) == 0 or (cap/'hsex-updated-000000001.bin').exists():
        raise AssertionError('Insufficient projection published W')
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    if sha(cap/'hsex-w.bin') != sha(root/'initial5/hsex-w.bin'):
        raise AssertionError('Cap rejection modified previous W')
    if (root/'salmon-release').exists():
        completed_log(root/'restart-release/run.log')
        result['release_restart_current'] = compare_current(root/'restart-release/Si_rt.data',
                                                           root/'grown/Si_rt.data', [2, 3])
        release_cap = root/'cap-release'
        if 'HSEX_PROJECTION_REJECT cap=5' not in (release_cap/'run.log').read_text():
            raise AssertionError('Final build did not reject insufficient capture')
        if int((release_cap/'exit-status').read_text()) == 0:
            raise AssertionError('Final build accepted insufficient capture')
        if (release_cap/'hsex-updated-000000001.bin').exists():
            raise AssertionError('Final build published W despite insufficient capture')
        if sha(release_cap/'hsex-w.bin') != sha(root/'initial5/hsex-w.bin'):
            raise AssertionError('Final cap rejection modified previous W')
        if checkpoint_count(root/'restart-release/checkpoint_rt_000002/info.bin') != 4:
            raise AssertionError('Final build changed propagated orbital count')
        result['release_sha256'] = sha(root/'salmon-release')
    result.update(initial_capture=initial, grown_capture=grown,
                  raw_population_max_error=float(np.max(abs(pop-raw))),
                  projected_electrons=float(pop.sum()/nk), auxiliary_states=nb,
                  propagated_states=4, cap_rejection_preserves_W=True,
                  same_nonzero_endpoint_A=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check(args.work)
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text)
