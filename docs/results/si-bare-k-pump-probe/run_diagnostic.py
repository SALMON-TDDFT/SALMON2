"""Check whether a delayed impulse reactivates full instantaneous screening."""
from pathlib import Path
import os,subprocess,json,sys
root=Path(__file__).resolve().parents[3]
case=sys.argv[1]
assert case in ('diagnostic_full','diagnostic_half')
d=root/'calculations/si_tdcdft_k4/bare_k_probe'/case
d.mkdir(parents=True,exist_ok=True)
for name,target in [('restart',root/'calculations/si_tdcdft_k4/gs/data_for_restart'),('Si_rps.dat',root/'testsuites/pseudo/Si_rps.dat')]:
    if not (d/name).exists(): (d/name).symlink_to(target)
s=(root/'calculations/si_tdcdft_k4/polarization/long/inputfile').read_text()
s=s[s.index('&calculation'):].replace('tdcdft_alpha=0.2','tdcdft_alpha=1.0').replace('tdcdft_screen_reference=0.19591042214256474','tdcdft_screen_reference=0.0')
nt=1200
s=s.replace('nt = 12000',f'nt = {nt}')
s=s.replace("ae_shape1 = 'Acos2'","ae_shape1 = 'Acos2'\n  ae_shape2 = 'impulse'\n  epdir_re2 = 0., 0., 1.").replace('T1_T2 = 620.1205995','T1_T2 = 50')
if case=='diagnostic_half':s=s.replace('e_impulse = 0.001','e_impulse = 0.0005')
(d/'inputfile').write_text(s)
env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1')
with (d/'outputfile').open('w') as f:
    r=subprocess.run(['/opt/homebrew/bin/mpiexec','-n','4','/private/tmp/salmon-tdcdft-mpi/salmon'],input=s,text=True,cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
status=dict(case=case,exit_code=r.returncode,requested_steps=nt,completed='end SALMON' in (d/'outputfile').read_text())
(root/'docs/results/si-bare-k-pump-probe'/f'{case}_status.json').write_text(json.dumps(status,indent=2)+'\n')
print(status,flush=True)
