"""Short stability screen, explicitly separate from full spectral validation."""
from analyze_results import *

def main():
 summary={};traces={}
 for name in ('g1_screen','g4_screen','g10_screen'):
  r,x,status=common.load(name,5000);summary[name]=status;traces[name]=(r,x)
 (out/'gamma_scan.json').write_text(json.dumps(summary,indent=2)+'\n')
 fig,ax=plt.subplots(3,1,figsize=(9,8),sharex=True,layout='constrained')
 for name,(r,x) in traces.items():
  stop=summary[name]['first_logged_norm_error_over_1e_4_fs'];keep=r[:,0]*fs<(stop if stop is not None else np.inf)
  t=r[keep,0]*fs;xx=x[keep];label=f"gamma={summary[name]['gamma']:g}"
  ax[0].plot(t,xx[:,7],label=label);ax[1].plot(t,xx[:,3],label=label);ax[2].plot(t,r[keep,15],label=label)
 ax[0].set(ylabel='alpha');ax[0].legend();ax[1].set(ylabel='Axc,z / c (a.u.)');ax[2].set(ylabel='Jz (a.u.)',xlabel='Time (fs)')
 for a in ax:a.axvline(60*fs,color='0.5',ls=':');a.grid(alpha=.2)
 fig.suptitle('Si strong pump: ELF-Proca stability screen, beta=0\nGamma=1e-4 trace ends before norm error; other runs reach9.68 fs')
 fig.savefig(out/'gamma_scan.png',dpi=160);fig.savefig(out/'gamma_scan.pdf')
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
