"""Produce Bloch snapshots; dynamics/alpha unchanged from the bare-K trial."""
from pathlib import Path
import subprocess,os,json,sys
root=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent
case=sys.argv[1];assert case in ('none','weak','strong')
work=Path('/private/tmp/salmon-si-time-wannier');d=work/case;d.mkdir(parents=True,exist_ok=True)
for name,target in [('restart',root/'calculations/si_tdcdft_k4/gs/data_for_restart'),('Si_rps.dat',root/'testsuites/pseudo/Si_rps.dat')]:
 if not (d/name).exists():(d/name).symlink_to(target)
s=(root/'calculations/si_tdcdft_k4/bare_k/strong_long/inputfile').read_text()
s=s.replace('nt = 12000','nt = 1600').replace("sysname = 'Si'","sysname = 'Si'\n  checkpoint_interval=200")
s=s.replace('&analysis',"&analysis\n  projection_option='gs'\n  out_projection_step=200")
intensity={'none':'0','weak':'100000000.0','strong':'10000000000000.0'}[case]
s=s.replace('I_wcm2_1 = 10000000000000.0','I_wcm2_1 = '+intensity)
(d/'inputfile').write_text(s);(out/f'{case}.inp').write_text(s)
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
with (d/'outputfile').open('w') as f:
 r=subprocess.run(['/opt/homebrew/bin/mpiexec','-n','4','/private/tmp/salmon-tdcdft-mpi/salmon'],input=s,text=True,cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
status=dict(case=case,exit_code=r.returncode,complete='end SALMON' in (d/'outputfile').read_text(),raw_directory=str(d))
(out/f'{case}_status.json').write_text(json.dumps(status,indent=2)+'\n');print(status,flush=True)
assert status['complete'] and r.returncode==0
