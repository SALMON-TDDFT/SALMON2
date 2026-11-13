"""Bounded native/Python MPI regression runner, requiring a host-MPI bridge.
python tools/run_parallel_mpi_tests.py BRIDGE_BUILD_JSON NEW_OUTPUT
"""
import argparse,json,os,subprocess,sys
from pathlib import Path
def run_bounded(command, *, env, log, timeout=60):
    """Terminate the launcher/process group on timeout, including stuck ranks."""
    import signal
    with subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                          start_new_session=True) as process:
        try:
            code=process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGTERM)
            try: process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL);process.wait()
            raise
        if code: raise subprocess.CalledProcessError(code,command)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('build');p.add_argument('output');p.add_argument('--tddft-runtime');p.add_argument('--tddft-bundle');a=p.parse_args()
    if bool(a.tddft_runtime)!=bool(a.tddft_bundle): p.error('supply both TDDFT runtime and bundle')
    if os.name!='posix': p.error('this MPI regression launcher currently targets macOS/Linux')
    j=json.loads(Path(a.build).read_text());out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False)
    root=Path(__file__).resolve().parents[1]
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    cases=[('custom',6,[str(root/'tools/test_distributed_mpi.py'),j['library'],j['sha256']]),
           ('ld',4,[str(root/'tools/test_ld_distributed_mpi.py')]),
           ('sample',5,['-m','salmonpy.examples.mpi_response',str(out/'sample.json'),'--bridge',j['library'],'--sha256',j['sha256']])]
    if a.tddft_runtime:
        cases += [('tddft',3,[str(root/'tools/test_tddft_distributed_mpi.py'),a.tddft_runtime,a.tddft_bundle]),
                  ('tddft-sample',3,['-m','salmonpy.examples.mpi_tddft',a.tddft_runtime,a.tddft_bundle,'rt.inp',str(out/'tddft-sample.json')])]
    for name,count,args in cases:
        with (out/(name+'.log')).open('w') as log:
            run_bounded(['mpiexec','-n',str(count),sys.executable,*args],env=env,log=log)
        if name == 'sample':
            import numpy as np
            sample = json.loads((out/'sample.json').read_text())
            assert sample['point_ids'] == [10, 20], 'one point per material group'
            assert len(sample['history']) == 4
            for step, frame in enumerate(sample['history'], 1):
                A = np.array([[point*.01*step, 0., 0.] for point in (10, 20)])
                expected = np.array([2., 3.])[:, None]*A + np.array([10., 20.])[:, None]*(.1*step)
                np.testing.assert_allclose(frame['A_over_c'], A, rtol=0, atol=1e-15)
                np.testing.assert_allclose(frame['J_charge'], expected, rtol=0, atol=1e-14)
        print(name,'PASS',flush=True)


if __name__=='__main__': main()
