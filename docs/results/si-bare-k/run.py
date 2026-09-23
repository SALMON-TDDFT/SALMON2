"""Alpha0=1, K0=0 trial using the existing instantaneous polarization closure."""
from pathlib import Path
import os,subprocess,json,sys
root=Path(__file__).resolve().parents[3]
case=sys.argv[1]
assert case in ('weak','strong','strong_long','strong_matched_gate')
d=root/'calculations/si_tdcdft_k4/bare_k'/case
d.mkdir(parents=True,exist_ok=True)
for name,target in [('restart',root/'calculations/si_tdcdft_k4/gs/data_for_restart'),('Si_rps.dat',root/'testsuites/pseudo/Si_rps.dat')]:
    if not (d/name).exists(): (d/name).symlink_to(target)
s=(root/'calculations/si_tdcdft_k4/polarization/long/inputfile').read_text()
s=s[s.index('&calculation'):].replace('tdcdft_alpha=0.2','tdcdft_alpha=1.0').replace('tdcdft_screen_reference=0.19591042214256474','tdcdft_screen_reference=0.0')
nt=12000 if case=='strong_long' else 1200
s=s.replace('nt = 12000',f'nt = {nt}')
if case=='weak': s=s.replace('I_wcm2_1 = 10000000000000.0','I_wcm2_1 = 100000000.0')
if case=='strong_matched_gate': s=s.replace('tdcdft_screen_floor=2e-5','tdcdft_screen_floor=0.006324555320336759')
(d/'inputfile').write_text(s)
env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1')
with (d/'outputfile').open('w') as f:
    r=subprocess.run(['/opt/homebrew/bin/mpiexec','-n','4','/private/tmp/salmon-tdcdft-mpi/salmon'],input=s,text=True,cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
status=dict(case=case,exit_code=r.returncode,requested_steps=nt,completed='end SALMON' in (d/'outputfile').read_text())
(root/'docs/results/si-bare-k'/f'{case}_status.json').write_text(json.dumps(status,indent=2)+'\n')
print(status,flush=True)
