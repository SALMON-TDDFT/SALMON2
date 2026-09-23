"""Extend the existing strong-pulse inputs to 12000 steps; change nt only."""
from pathlib import Path
import json, os, subprocess
import numpy as np
root=Path(__file__).resolve().parents[3]
source=root/'calculations/si_tdcdft_k4/instant_smoke'
work=root/'calculations/si_tdcdft_k4/instant_long'
work.mkdir(exist_ok=True)
exe=os.environ.get('SALMON_EXE','/private/tmp/salmon-tdcdft-mpi/salmon')
env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1')
results={}
for name in ['strong_fixed','strong_screened']:
    d=work/name;d.mkdir(exist_ok=True)
    for target,src in [('restart',root/'calculations/si_tdcdft_k4/gs/data_for_restart'),
                       ('Si_rps.dat',root/'testsuites/pseudo/Si_rps.dat')]:
        if not (d/target).exists(): (d/target).symlink_to(src)
    old=(source/name/'inputfile').read_text()
    assert old.count('nt = 1200')==1
    inp=old.replace('nt = 1200','nt = 12000')
    (d/'inputfile').write_text(inp)
    print('Running',name,flush=True)
    with (d/'outputfile').open('w') as out:
        result=subprocess.run(['/opt/homebrew/bin/mpiexec','-n','4',exe],input=inp,text=True,
                              cwd=d,env=env,stdout=out,stderr=subprocess.STDOUT)
    results[name]={'exit_code':result.returncode}
    if (d/'Si_rt_xc.data').exists():
        x=np.loadtxt(d/'Si_rt_xc.data')
        results[name].update(rows=len(x),finite=bool(np.isfinite(x).all()))
        previous=np.loadtxt(source/name/'Si_rt_xc.data')
        if len(x)>=len(previous):
            delta=float(np.max(abs(x[:len(previous)]-previous)))
            results[name]['short_prefix_max_difference']=delta
            np.testing.assert_allclose(x[:len(previous)],previous,rtol=1e-10,atol=1e-12)
    print(name,results[name],flush=True)
    (root/'docs/results/si-instant-screening/long_status.json').write_text(json.dumps(results,indent=2)+'\n')
