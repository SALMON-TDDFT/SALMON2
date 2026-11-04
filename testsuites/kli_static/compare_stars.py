"""Compare a complete reduced-mesh Si test to the matching full-mesh trajectory."""
import argparse
from pathlib import Path
import re
import numpy as np


def gaps(run):
    log = (run / 'warm-fast' / 'run.log').read_text()
    assert '#GS converged' in log and 'end SALMON' in log
    return np.array([float(re.findall(pattern, log)[-1]) for pattern in
                     (r'Fundamental gap\[eV\]\s+(\S+)', r'BG between same k-point\[eV\]\s+(\S+)')])


def current(run, case, steps):
    rows = []
    with (run / case / 'Si_rt.data').open() as stream:
        for line in stream:
            if line.strip() and not line.startswith('#'):
                rows.append(line)
            if len(rows) == steps:
                break
    assert len(rows) == steps
    return np.loadtxt(rows, ndmin=2)[:, [0, 13, 14, 15]]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('full', type=Path)
    parser.add_argument('reduced', type=Path)
    parser.add_argument('--steps', type=int, default=20)
    args = parser.parse_args()
    counts = re.search(r'KLI representative/full k points:\s+(\d+)\s+(\d+)',
                       (args.reduced / 'warm-fast' / 'run.log').read_text())
    assert counts and int(counts[1]) < int(counts[2])
    np.testing.assert_allclose(gaps(args.full), gaps(args.reduced), rtol=0, atol=1e-5)
    for mode in ('frozen8', 'dynamic8'):
        assert 'end SALMON' in (args.reduced / mode / 'run.log').read_text()
        full = current(args.full, mode, args.steps)
        reduced = current(args.reduced, mode, args.steps)
        np.testing.assert_array_equal(full[:, 0], reduced[:, 0])
        np.testing.assert_allclose(full[:, 1:], reduced[:, 1:], atol=1e-9, rtol=0)
        print(mode, 'max |J_reduced - J_full| [a.u.]', abs(full[:, 1:] - reduced[:, 1:]).max())
