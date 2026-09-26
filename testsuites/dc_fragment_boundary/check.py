"""Check the actual MPI fragment initializer using the Si64 axial fixture.

Usage: python3 check.py /absolute/path/to/salmon /new/run/directory
Runs one SCF iteration; does not establish SCF convergence.
"""
from pathlib import Path
import os
import re
import shutil
import subprocess
import sys

exe = Path(sys.argv[1]).resolve()
run = Path(sys.argv[2]).resolve()
fixture = Path(__file__).resolve().parents[2] / 'samples/dc_hse/si64-chain'
run.mkdir(parents=True, exist_ok=False)
for name in ('atom.dat', 'Si_rps.dat'):
    shutil.copy2(fixture / name, run / name)
input_text = (fixture / 'inputfile.reference').read_text()
input_text = re.sub(r'\bnscf\s*=\s*1000\b', 'nscf=1', input_text)
input_text = re.sub(r'\bcheckpoint_interval\s*=\s*10\b', 'checkpoint_interval=1', input_text)
(run / 'inputfile').write_text(input_text)
env = dict(os.environ, OMP_NUM_THREADS='2', OMP_MAX_ACTIVE_LEVELS='1', OPENBLAS_NUM_THREADS='1')
with (run / 'inputfile').open('rb') as inp, (run / 'outputfile').open('wb') as out:
    subprocess.run(['mpiexec', '-n', '8', str(exe)], cwd=run, stdin=inp,
                   stdout=out, stderr=subprocess.STDOUT, env=env, check=True, timeout=90)
text = (run / 'outputfile').read_text()
rows = sorted(tuple(map(int, m)) for m in re.findall(
    r'fragment, natom, nelec, ixyz_frag:\s+(\d+)\s+(\d+)\s+(\d+)', text))
expected = [(i, 16, 64) for i in range(1, 9)]
assert rows == expected, f'fragment/atoms/initial electrons: {rows}; expected {expected}'
sites = []
for i in range(1, 9):
    path = run / f'data_dcdft/fragments/{i:06d}/checkpoint_gs_000001/atomic_coor.txt'
    coordinates = [list(map(float, line.split()[1:4])) for line in path.read_text().splitlines()]
    periodic_sites = set()
    for xyz in coordinates:
        fractional = [x / 2.565 for x in xyz]
        assert all(abs(x - round(x)) < 1e-11 for x in fractional), xyz
        periodic_sites.add(tuple(round(x) % n for x, n in zip(fractional, (8, 4, 4))))
    assert len(periodic_sites) == 16, (i, periodic_sites)
    sites.append(periodic_sites)
assert all(s == sites[0] for s in sites), 'periodic fragment geometries are not equivalent'
print('PASS: eight equivalent periodic fragments, each with 16 distinct atoms and initially 64 electrons')
