"""Render the recorded exchange timing/accuracy measurements; no simulations."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(__file__).resolve().parent
fig,ax=plt.subplots(1,2,figsize=(10.5,4),layout='constrained')
for name,label in [('hse','HSE ground state'),('initial','PZ ground state'),('strong','TDCDFT pumped snapshot')]:
 d=json.loads((p/f'{name}_optimized_benchmark.json').read_text());r=[x for x in d['rows'] if x['action_relative_error']>0]
 ax[0].plot([x['median_seconds'] for x in r],[max(x['action_relative_error'],1e-16) for x in r],'o-',label=label)
 ax[1].plot([x['median_seconds'] for x in r],[abs(x['hse_energy_error_meV_per_atom']) for x in r],'o-',label=label)
ax[0].set(yscale='log',xlabel='Exchange action wall time (s)',ylabel='Relative action error',ylim=(1e-5,.2),title='Pair-screening accuracy (nonzero cutoffs)')
ax[1].set(yscale='log',xlabel='Exchange action wall time (s)',ylabel='HSE exchange error (meV/atom)',title='Same snapshot, full-pair reference')
for a in ax:a.grid(alpha=.25);a.legend(fontsize=8)
fig.suptitle('Si8 / 4³ k points — single-thread FFTW, median of 3 runs',fontsize=12)
fig.savefig(p/'exchange_efficiency.png',dpi=180)
