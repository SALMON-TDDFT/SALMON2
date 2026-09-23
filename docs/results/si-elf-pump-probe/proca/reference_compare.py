"""Compare completed unpumped references at common delay/window."""
from analyze_results import *

def main():
 baseline=np.loadtxt(out.parent/'ground_spectrum.csv',delimiter=',',skiprows=1)
 spectra={'gamma=0':baseline[:,2]};metrics={}
 for name in ('g4_ground','g10_ground'):
  if not (work/name/'status.json').exists():continue
  r,x,status=common.load(name)
  if not status['valid_trajectory']:raise RuntimeError(f'{name}: invalid reference')
  eps=response(r[:,0],r[:,15],eta,energy,probe);spectra[f"gamma={status['gamma']:g}"]=eps.imag
  m=peak_metrics(energy,eps.imag,2,4);m['apparent_FWHM_eV']=apparent_width(energy,eps.imag)
  m['max_norm_error']=status['max_norm_error'];m['gamma']=status['gamma'];m['windows']={}
  for span in (600.,720.,880.):
   keep=r[:,0]<=probe+span+1e-8;ee=response(r[keep,0],r[keep,15],eta,energy,probe);m['windows'][str(span)]=peak_metrics(energy,ee.imag,2,4)
  metrics[name]=m
  np.savetxt(out/f'{name}_spectrum.csv',np.column_stack((energy,eps.real,eps.imag)),delimiter=',',header='energy_eV,Re_epsilon,Im_epsilon',comments='')
 metrics['zero_gamma']=peak_metrics(baseline[:,0],baseline[:,2],2,4)
 (out/'reference_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
 fig,ax=plt.subplots(2,1,figsize=(9,7),layout='constrained')
 for label,y in spectra.items():
  for a,limits in zip(ax,[(2,4.5),(.1,1.5)]):
   show=(energy>=limits[0])&(energy<=limits[1]);a.plot(energy[show],y[show],label=label)
 ax[0].set(xlim=(2,4.5),xlabel='Energy (eV)',ylabel='Im epsilon',title='Optical peak; unpumped references')
 ax[1].set(xlim=(.1,1.5),xlabel='Energy (eV)',ylabel='Im epsilon',title='Auxiliary low-frequency response; not an exciton peak')
 for a in ax:a.legend();a.axhline(0,color='0.5',lw=.6);a.grid(alpha=.2)
 fig.suptitle('Si 4x4x4 k: effect of Proca stabilization on the reference')
 fig.savefig(out/'reference_comparison.png',dpi=160);fig.savefig(out/'reference_comparison.pdf')
 print(json.dumps(metrics,indent=2))
if __name__=='__main__':main()
