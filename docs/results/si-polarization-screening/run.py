"""Same laser test for the E_xc=alpha P closure; retain all old comparison data."""
from pathlib import Path
import json,os,subprocess,sys
import numpy as np
root=Path(__file__).resolve().parents[3]
source=root/'calculations/si_tdcdft_k4/instant_long/strong_screened/inputfile'
work=root/'calculations/si_tdcdft_k4/polarization'
work.mkdir(exist_ok=True)
case=sys.argv[1] if len(sys.argv)>1 else 'long'
assert case in ('long','fine')
nt,dt=(12000,.08) if case=='long' else (6000,.04)
d=work/case;d.mkdir(exist_ok=True)
for name,target in [('restart',root/'calculations/si_tdcdft_k4/gs/data_for_restart'),
                    ('Si_rps.dat',root/'testsuites/pseudo/Si_rps.dat')]:
    if not (d/name).exists(): (d/name).symlink_to(target)
text=source.read_text().replace("tdcdft_screening='instant'","tdcdft_screening='polarization'")
if case=='fine': text=text.replace('nt = 12000','nt = 6000').replace('dt = 0.08','dt = 0.04')
(d/'inputfile').write_text(text)
exe=os.environ.get('SALMON_EXE','/private/tmp/salmon-tdcdft-mpi/salmon')
env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1')
with (d/'outputfile').open('w') as out:
    result=subprocess.run(['/opt/homebrew/bin/mpiexec','-n','4',exe],input=text,text=True,cwd=d,env=env,
                          stdout=out,stderr=subprocess.STDOUT)
x=np.loadtxt(d/'Si_rt_xc.data');j=np.loadtxt(d/'Si_rt.data')
status=dict(exit_code=result.returncode,nt=nt,dt_au=dt,rows=len(x),finite=bool(np.isfinite(x).all() and np.isfinite(j).all()))
(root/f'docs/results/si-polarization-screening/{case}_status.json').write_text(json.dumps(status,indent=2)+'\n')
print(status,flush=True)
assert result.returncode==0 and len(x)==nt and status['finite']
