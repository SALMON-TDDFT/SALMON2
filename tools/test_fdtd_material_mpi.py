"""Compare magnetic interfaces across native MPI decompositions.
Usage: python tools/test_fdtd_material_mpi.py MPI_EXECUTABLE SERIAL_FIXTURE NEW_OUTPUT
SERIAL_FIXTURE is magnetic_dielectric_0 from test_fdtd_material_interfaces.py.
The binary must already pass the source/object audit. This tests the native
launcher; it does not add a public multi-rank Python adapter.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('executable'); p.add_argument('fixture'); p.add_argument('output')
a = p.parse_args()
root = Path(a.output).resolve(); root.mkdir(parents=True, exist_ok=False)
text = (Path(a.fixture)/'inputfile').read_text()
# Put the lower y/z surface on the rank corner. At some H sites only
# the diagonal adjacent cell is material, with permeability != vacuum.
assert text.count('ori_s(1,:) = -2, 0, 0') == 1
assert text.count('ori_s(2,:) = 2, 0, 0') == 1
text = text.replace('ori_s(1,:) = -2, 0, 0', 'ori_s(1,:) = -2, 8, 8')
text = text.replace('ori_s(2,:) = 2, 0, 0', 'ori_s(2,:) = 2, 8, 8')
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
reference = None
report = {}
try:
    for dims in [(1,1,1),(1,2,2),(2,1,2),(2,2,1)]:
        label = 'x'.join(map(str,dims)); run = root/label; run.mkdir()
        inp = text + '\n&parallel\n nproc_k=1\n nproc_ob=1\n nproc_rgrid='+','.join(map(str,dims))+'\n/\n'
        (run/'inputfile').write_text(inp)
        with (run/'inputfile').open() as stdin, (run/'stdout.log').open('w') as log:
            subprocess.run(['mpiexec','-n',str(np.prod(dims)),str(Path(a.executable).resolve())],
                           stdin=stdin,stdout=log,stderr=subprocess.STDOUT,cwd=run,env=env,check=True,timeout=300)
        assert 'SALMONPY_FDTD_INTERFACE_V2' in (run/'stdout.log').read_text()
        values = [np.loadtxt(f) for f in sorted(run.glob('obs*_at_point_rt.data'))]
        assert len(values)==2 and all(len(v)==160 and np.isfinite(v).all() for v in values)
        if reference is None: reference = values
        diff = max(float(np.max(abs(x-y))) for x,y in zip(reference,values))
        report[label] = {'max_abs':diff,'bitwise_equal': all(np.array_equal(x,y) for x,y in zip(reference,values))}
        for x,y in zip(reference,values): np.testing.assert_allclose(x,y,rtol=1e-12,atol=1e-14)
        print('PASS',label,diff,flush=True)
finally:
    (root/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
