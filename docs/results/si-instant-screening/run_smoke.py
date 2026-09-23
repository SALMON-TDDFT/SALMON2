"""Short 4^3-k laser comparison; numerical trial, not a converged transient spectrum."""
from pathlib import Path
import json,os,subprocess,sys
import numpy as np
root=Path(__file__).resolve().parents[3]
work=root/'calculations/si_tdcdft_k4/instant_smoke'
work.mkdir(exist_ok=True)
base=(root/'samples/exercise_si_tdcdft/Si_pump.inp').read_text()
base=base.replace('nt = 30000','nt = 1200').replace('tw1 = 413.413733','tw1 = 60')
base=base.replace('omega1 = 0.05879892','omega1 = 0.2')
base=base.replace("tdcdft='proca'","tdcdft='lrc'").replace('tdcdft_restoring=0.0001','tdcdft_restoring=0')
env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1')
exe=os.environ.get('SALMON_EXE','/private/tmp/salmon-tdcdft-mpi/salmon')
def run(name,intensity,strength,reference=0,floor=2e-5):
    d=work/name;d.mkdir(exist_ok=True)
    for target,src in [('restart',root/'calculations/si_tdcdft_k4/gs/data_for_restart'),
                       ('Si_rps.dat',root/'testsuites/pseudo/Si_rps.dat')]:
        if not (d/target).exists(): (d/target).symlink_to(src)
    inp=base.replace('I_wcm2_1 = 1e10',f'I_wcm2_1 = {intensity}')
    inp=inp.replace('tdcdft_alpha=0.2',f'''tdcdft_alpha=0.2
  tdcdft_screening='instant'
  tdcdft_screen_omega=0.2
  tdcdft_screen_reference={reference:.17g}
  tdcdft_screen_strength={strength}
  tdcdft_screen_floor={floor:.17g}''')
    (d/'inputfile').write_text(inp)
    print('Running',name,flush=True)
    with (d/'outputfile').open('w') as out:
        subprocess.run(['/opt/homebrew/bin/mpiexec','-n','4',exe],input=inp,text=True,cwd=d,env=env,
                       stdout=out,stderr=subprocess.STDOUT,check=True)
    x=np.loadtxt(d/'Si_rt_xc.data');j=np.loadtxt(d/'Si_rt.data')
    assert x.shape==(1200,12) and np.isfinite(x).all() and np.isfinite(j).all()
    return x,j
if '--floor-check' in sys.argv:
    metrics=json.loads((root/'docs/results/si-instant-screening/metrics.json').read_text())
    active=np.loadtxt(work/'strong_screened/Si_rt_xc.data')
    aj=np.loadtxt(work/'strong_screened/Si_rt.data')
    results=[]
    for floor in [1e-5,4e-5]:
        x,j=run('strong_floor_'+str(floor),1e13,1,metrics['reference_K'],floor)
        results.append(dict(floor=floor,alpha_min=float(x[:,7].min()),alpha_end=float(x[-1,7]),
            relative_current_change_from_default=float(np.linalg.norm(j[:,13:16]-aj[:,13:16])/np.linalg.norm(aj[:,13:16]))))
    (root/'docs/results/si-instant-screening/floor_sensitivity.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2),flush=True)
else:
    weak,wj=run('weak_reference',1e8,0)
    k0=max(0,float(weak[:,8].max()))+1e-12
    fixed,fj=run('strong_fixed',1e13,0,k0)
    active,aj=run('strong_screened',1e13,1,k0)
    check,cj=run('weak_screened',1e8,1,k0)
    metrics=dict(k_grid=[4,4,4],dt_au=.08,nt=1200,pulse_duration_au=60,omega_au=.2,
                 reference_K=k0,field_floor_au=2e-5,weak_intensity=1e8,strong_intensity=1e13,
                 weak_max_current_difference=float(np.max(abs(wj[:,13:16]-cj[:,13:16]))),
                 strong_relative_current_change=float(np.linalg.norm(aj[:,13:16]-fj[:,13:16])/np.linalg.norm(fj[:,13:16])),
                 alpha_min=float(active[:,7].min()),alpha_max=float(active[:,7].max()),
                 fraction_screened=float(np.mean(active[:,7]<.2)),
                 max_Axc_fixed=float(np.max(abs(fixed[:,1:4]))),max_Axc_screened=float(np.max(abs(active[:,1:4]))))
    (root/'docs/results/si-instant-screening/metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    print(json.dumps(metrics,indent=2),flush=True)
