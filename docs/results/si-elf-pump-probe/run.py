"""Native ELF pump-probe: select cases with ELF_PROBE_CASES; no alpha freezing."""
from pathlib import Path
import os,sys,subprocess,json,time,concurrent.futures,fcntl
out=Path(__file__).resolve().parent;root=out.parents[2]
work=Path('/private/tmp/salmon-si-elf-pump-probe');work.mkdir(exist_ok=True)
exe=Path(sys.argv[1]).resolve();launcher=sys.argv[2:]
base=(out.parent/'si-elf-feedback/strong.inp').read_text().replace('nt = 1600','nt = 12000').replace('checkpoint_interval=200','checkpoint_interval=2000').replace("projection_option='gs'","projection_option='no'")
base=base.replace("ae_shape1 = 'Acos2'","ae_shape1 = 'Acos2'\n  ae_shape2 = 'impulse'\n  epdir_re2 = 0., 0., 1.\n  T1_T2 = 50")
amplitudes={'pump':0.,'plus':1e-4,'minus':-1e-4,'half_plus':5e-5,'half_minus':-5e-5,'ground':1e-4,'ground_half':5e-5,'ground_zero':0.,'diagnostic_halfdt':0.,'diagnostic_stride1':0.}
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
def run(name):
 inp=base.replace('e_impulse = 0.001',f'e_impulse = {amplitudes[name]:.17g}')
 if name.startswith('ground'):inp=inp.replace('I_wcm2_1 = 10000000000000.0','I_wcm2_1 = 0')
 if name=='diagnostic_halfdt':inp=inp.replace('nt = 12000','nt = 6240').replace('dt = 0.08','dt = 0.04').replace('tdcdft_elf_stride=10','tdcdft_elf_stride=20')
 if name=='diagnostic_stride1':inp=inp.replace('nt = 12000','nt = 3120').replace('tdcdft_elf_stride=10','tdcdft_elf_stride=1')
 p=work/name;p.mkdir(exist_ok=True)
 with (p/'run.lock').open('w') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  if (p/'status.json').exists():
   if (p/'inputfile').read_text()!=inp:raise RuntimeError(f'{name}: input changed; choose a fresh work directory')
   return json.loads((p/'status.json').read_text())
  for link,target in [('restart',root/'calculations/si_tdcdft_k4/gs/data_for_restart'),('Si_rps.dat',root/'testsuites/pseudo/Si_rps.dat')]:
   if not (p/link).exists():(p/link).symlink_to(target)
  (p/'inputfile').write_text(inp);(out/f'{name}.inp').write_text(inp);start=time.monotonic()
  with (p/'outputfile').open('w') as log:
   result=subprocess.run(launcher+[str(exe)],input=inp,text=True,cwd=p,stdout=log,stderr=subprocess.STDOUT,env=env)
  status=dict(case=name,exit_code=result.returncode,completed=result.returncode==0 and 'end SALMON' in (p/'outputfile').read_text(),elapsed_seconds=time.monotonic()-start)
  (p/'status.json').write_text(json.dumps(status,indent=2)+'\n');(out/f'{name}_status.json').write_text(json.dumps(status,indent=2)+'\n')
  return status
names=os.environ.get('ELF_PROBE_CASES','pump,plus,ground').split(',')
assert all(x in amplitudes for x in names)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 for job in concurrent.futures.as_completed([pool.submit(run,x) for x in names]):print(json.dumps(job.result()),flush=True)
