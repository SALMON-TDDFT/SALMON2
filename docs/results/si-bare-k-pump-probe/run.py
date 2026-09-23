"""Frozen pump-screening response at t_probe=80 au, same Si pump and 4^3 k grid."""
from pathlib import Path
import os,subprocess,json,sys,fcntl
root=Path(__file__).resolve().parents[3]
case=sys.argv[1]
assert case in ('pumped','pumped_half','pumped_small','ground_one','ground_screened')
d=root/'calculations/si_tdcdft_k4/bare_k_probe'/case;d.mkdir(parents=True,exist_ok=True)
lock=(d/'run.lock').open('w')
fcntl.flock(lock,fcntl.LOCK_EX)
statusfile=root/'docs/results/si-bare-k-pump-probe'/f'{case}_status.json'
if statusfile.exists() and json.loads(statusfile.read_text()).get('completed'):
 print(f'{case}: existing completed trajectory retained',flush=True);sys.exit(0)
for name,target in [('restart',root/'calculations/si_tdcdft_k4/gs/data_for_restart'),('Si_rps.dat',root/'testsuites/pseudo/Si_rps.dat')]:
 if not (d/name).exists():(d/name).symlink_to(target)
s=(root/'calculations/si_tdcdft_k4/bare_k/strong_long/inputfile').read_text()
s=s.replace("tdcdft_screening='polarization'","tdcdft_screening='polarization'\n  tdcdft_screen_stop=60")
s=s.replace("ae_shape1 = 'Acos2'","ae_shape1 = 'Acos2'\n  ae_shape2 = 'impulse'\n  epdir_re2 = 0., 0., 1.").replace('T1_T2 = 620.1205995','T1_T2 = 50')
if case=='pumped_small':s=s.replace('e_impulse = 0.001','e_impulse = 0.0001')
if case=='pumped_half':s=s.replace('e_impulse = 0.001','e_impulse = 0.0005')
if case.startswith('ground'):s=s.replace('I_wcm2_1 = 10000000000000.0','I_wcm2_1 = 0')
if case=='ground_screened':s=s.replace('tdcdft_alpha=1.0','tdcdft_alpha=0.0002062912827269997')
(d/'inputfile').write_text(s)
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
with (d/'outputfile').open('w') as f:
 r=subprocess.run(['/opt/homebrew/bin/mpiexec','-n','4','/private/tmp/salmon-tdcdft-mpi/salmon'],input=s,text=True,cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
status=dict(case=case,exit_code=r.returncode,requested_steps=12000,completed='end SALMON' in (d/'outputfile').read_text())
(root/'docs/results/si-bare-k-pump-probe'/f'{case}_status.json').write_text(json.dumps(status,indent=2)+'\n')
print(status,flush=True)
assert status['completed'] and r.returncode==0
