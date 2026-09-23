"""ELF + Proca pump-probe. Cases: g1/g4/g10/g40 with screen/pump/plus/minus/halfplus/halfminus/ground."""
from pathlib import Path
import os,sys,json,subprocess,time,fcntl,concurrent.futures
out=Path(__file__).resolve().parent;root=out.parents[3];work=Path('/private/tmp/salmon-si-elf-proca-probe');work.mkdir(exist_ok=True)
exe=Path(sys.argv[1]).resolve();launcher=sys.argv[2:]
base=(out.parent/'pump.inp').read_text().replace("tdcdft='lrc'","tdcdft='proca'")
gammas={'g1':.0001,'g4':.0004,'g10':.001,'g40':.004}
amplitudes={'screen':0.,'pump':0.,'plus':1e-4,'minus':-1e-4,'halfplus':5e-5,'halfminus':-5e-5,'ground':1e-4,'groundhalf':5e-5}
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
def run(name):
 coefficient,kind=name.split('_');gamma=gammas[coefficient];amplitude=amplitudes[kind]
 inp=base.replace('tdcdft_restoring=0',f'tdcdft_restoring={gamma:.17g}').replace('e_impulse = 0','e_impulse = '+format(amplitude,'.17g'))
 steps=5000 if kind=='screen' else 12000
 inp=inp.replace('nt = 12000',f'nt = {steps}')
 if kind.startswith('ground'):inp=inp.replace('I_wcm2_1 = 10000000000000.0','I_wcm2_1 = 0')
 p=work/name;p.mkdir(exist_ok=True)
 with (p/'run.lock').open('w') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  if (p/'status.json').exists():
   if (p/'inputfile').read_text()!=inp:raise RuntimeError(f'{name}: saved input differs')
   return json.loads((p/'status.json').read_text())
  for link,target in [('restart',root/'calculations/si_tdcdft_k4/gs/data_for_restart'),('Si_rps.dat',root/'testsuites/pseudo/Si_rps.dat')]:
   if not (p/link).exists():(p/link).symlink_to(target)
  (p/'inputfile').write_text(inp);(out/f'{name}.inp').write_text(inp);start=time.monotonic()
  with (p/'outputfile').open('w') as f:result=subprocess.run(launcher+[str(exe)],input=inp,text=True,cwd=p,stdout=f,stderr=subprocess.STDOUT,env=env)
  status=dict(case=name,gamma=gamma,beta=0.,probe_amplitude=amplitude,requested_steps=steps,exit_code=result.returncode,
    completed=result.returncode==0 and 'end SALMON' in (p/'outputfile').read_text(),elapsed_seconds=time.monotonic()-start)
  (p/'status.json').write_text(json.dumps(status,indent=2)+'\n');(out/f'{name}_status.json').write_text(json.dumps(status,indent=2)+'\n');return status
names=os.environ.get('ELF_PROCA_CASES','g1_screen,g4_screen,g10_screen').split(',')
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 for job in concurrent.futures.as_completed([pool.submit(run,n) for n in names]):print(json.dumps(job.result()),flush=True)
