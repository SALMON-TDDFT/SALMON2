"""Bounded local regressions for the experimental embedded MPI coupling.
python tools/run_coupled_mpi_tests.py EMBEDDED_BUILD_JSON SERIAL_RUNTIME GS_BUNDLE NEW_OUTPUT
"""
import argparse,json,os,sys
from pathlib import Path
import numpy as np
from run_parallel_mpi_tests import run_bounded

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('build');p.add_argument('runtime');p.add_argument('bundle');p.add_argument('output')
a=p.parse_args();out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False)
root=Path(__file__).resolve().parents[1]
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
cases=[]
for name,n,extra in [('orbital',5,[]),('unequal',4,[]),('realspace',5,['--rgrid'])]:
    cases.append((name,n,[str(root/'tools/test_embedded_tddft_mpi.py'),a.build,a.runtime,a.bundle,str(out/name),*extra]))
for mode,n,axis in [('zero',2,0),('zero',3,0),('ld',3,0),('current',3,0),('current',3,1),('current',3,2),('tddft',5,0),('tddft',6,0),('mixed',6,0)]:
    name=f'{mode}-{n}-axis{axis}'
    cases.append((name,n,[str(root/'tools/test_coupled_mpi.py'),mode,a.build,a.runtime,a.bundle,str(out/name),str(axis)]))
cases.append(('sample',6,['-m','salmonpy.examples.mpi_coupled_tddft',a.build,a.bundle,'rt.inp',str(out/'sample')]))
for name,n,args in cases:
    with (out/(name+'.log')).open('w') as log:
        run_bounded(['mpiexec','-n',str(n),sys.executable,*args],env=env,log=log,timeout=60)
    print(name,'PASS',flush=True)
full=json.loads((out/'tddft-6-axis0/history.json').read_text())
zero=json.loads((out/'zero-3-axis0/history.json').read_text())
single=json.loads((out/'tddft-5-axis0/history.json').read_text())
for key in ('A_over_c','J_charge','E_at_points'):
    np.testing.assert_array_equal([v[key] for v in full],[v[key] for v in single])
delta=float(np.max(np.abs(np.array([v['A_over_c'] for v in full])-np.array([v['A_over_c'] for v in zero]))))
assert delta>1e-10,'feedback must differ from zero-current propagation'
(out/'summary.json').write_text(json.dumps({'cases_passed':len(cases),'field_rank_parity':'exact',
    'feedback_max_A_difference':delta,'scope':'connection regressions, not material convergence'},indent=2)+'\n')
print('feedback and rank parity PASS',delta)
