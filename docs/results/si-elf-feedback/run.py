"""Run the native ELF feedback comparison; EXE [MPI launcher prefix]."""
from pathlib import Path
import os,sys,subprocess,concurrent.futures,time
out=Path(__file__).resolve().parent;root=out.parents[2]
exe=Path(sys.argv[1]).resolve();launcher=sys.argv[2:]
work=Path('/private/tmp/salmon-si-elf-feedback');work.mkdir(exist_ok=True)
base=(out.parent/'si-time-wannier/strong_dense.inp').read_text()
base=base.replace('checkpoint_interval=10','checkpoint_interval=200').replace('tdcdft_alpha=1.0','tdcdft_alpha=0.2').replace("tdcdft_screening='polarization'","tdcdft_screening='elf'\n  tdcdft_elf_stride=10")
for line in ('  tdcdft_screen_omega=0.2\n','  tdcdft_screen_reference=0.0\n','  tdcdft_screen_strength=1\n','  tdcdft_screen_floor=2e-5\n','  T1_T2 = 620.1205995\n'):
 base=base.replace(line,'')
cases={name:base.replace('I_wcm2_1 = 10000000000000.0',f'I_wcm2_1 = {intensity}') for name,intensity in [('none',0),('weak',1e8),('strong',1e13)]}
cases['fixed_strong']=cases['strong'].replace("tdcdft_screening='elf'","tdcdft_screening='none'")
imp=base.replace("theory='tddft_pulse'","theory='tddft_response'").replace("ae_shape1 = 'Acos2'","ae_shape1 = 'impulse'").replace('e_impulse = 0.001','e_impulse = 0.0001')
cases['impulse']=imp;cases['fixed_impulse']=imp.replace("tdcdft_screening='elf'","tdcdft_screening='none'")
cases['strong_stride1']=cases['strong'].replace('tdcdft_elf_stride=10','tdcdft_elf_stride=1')
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
def run(item):
 name,inp=item;p=work/name;p.mkdir(exist_ok=True);(out/f'{name}.inp').write_text(inp)
 if (p/'outputfile').exists() and 'end SALMON' in (p/'outputfile').read_text():
  if (p/'inputfile').read_text()!=inp:raise RuntimeError(f'{name}: saved input differs; choose a fresh work directory')
  return name,'already complete'
 for link,target in [('restart',root/'calculations/si_tdcdft_k4/gs/data_for_restart'),('Si_rps.dat',root/'testsuites/pseudo/Si_rps.dat')]:
  if not (p/link).exists():(p/link).symlink_to(target)
 (p/'inputfile').write_text(inp);t=time.monotonic()
 with (p/'outputfile').open('w') as f:
  r=subprocess.run(launcher+[str(exe)],input=inp,text=True,cwd=p,stdout=f,stderr=subprocess.STDOUT,env=env)
 if r.returncode or 'end SALMON' not in (p/'outputfile').read_text():raise RuntimeError(f'{name}: {p}/outputfile')
 return name,round(time.monotonic()-t,1)
# Cases may be selected through ELF_CASES without changing launcher arguments.
items=[(k,v) for k,v in cases.items() if k in os.environ.get('ELF_CASES',','.join(cases)).split(',')]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 for job in concurrent.futures.as_completed([pool.submit(run,x) for x in items]):print(job.result(),flush=True)
