"""Postprocess completed native impulse traces with the established response convention."""
from pathlib import Path
import sys,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'samples/exercise_si_tdcdft'))
from analyze import read_current,response,peak_metrics
OUT=Path(__file__).resolve().parent
BASE=ROOT/'calculations/si_hse_native'
energy=np.arange(.01,8.001,.01)
traces={}
for name,directory in [('HSE06','impulse_hse06_fixedalpha'),('TDCDFT alpha=0.2','impulse_tdcdft_fixedalpha')]:
 d=BASE/directory
 if (d/'status.txt').read_text().strip()!='completed':raise RuntimeError(f'{name} incomplete')
 a=read_current(d/'Si_rt.data')
 if not np.isclose(a[-1,0],880):raise RuntimeError(f'{name} wrong end time')
 if not np.allclose(a[:,3],1e-4,atol=1e-14,rtol=0):raise RuntimeError('Wrong impulse')
 traces[name]=a
metrics={};fig,axs=plt.subplots(2,1,figsize=(8,7),sharex=True)
for T in (440,660,880):
 columns=[energy];metrics[str(T)]={}
 for name,a in traces.items():
  data=a[a[:,0]<=T+1e-9]
  eps=response(data[:,0],data[:,15],1e-4,energy)
  columns.extend([eps.real,eps.imag])
  metrics[str(T)][name]={f'{lo}-{hi}eV':peak_metrics(energy,eps.imag,lo,hi) for lo,hi in [(2,4),(4,6)]}
  if T==880:
   axs[0].plot(energy,eps.real,label=name);axs[1].plot(energy,eps.imag,label=name)
 np.savetxt(OUT/f'spectra_{T}au.csv',np.column_stack(columns),delimiter=',',header='energy_eV,HSE_Re,HSE_Im,TDCDFT_Re,TDCDFT_Im')
for ax,label in zip(axs,['Re epsilon','Im epsilon']):ax.set_ylabel(label);ax.grid(alpha=.25);ax.legend()
axs[1].set_xlabel('Energy (eV)');fig.suptitle('Si impulse: HSE06 vs fixed-alpha TDCDFT (21.29 fs)');fig.tight_layout()
fig.savefig(OUT/'comparison.png',dpi=180);fig.savefig(OUT/'comparison.pdf')
(OUT/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
print('Analysis completed; inspect plots and finite-window dependence before interpreting peaks.')
