"""Record why a valid pumped spectrum cannot be extracted; compare short controls."""
from analyze_probe import *

def main():
 data={};summary={}
 for name,steps in [('pump',12000),('plus',12000),('ground',12000),('diagnostic_halfdt',6240),('diagnostic_stride1',3120)]:
  r,x,s=load(name,steps);data[name]=(r,x);summary[name]=s
 r,x=data['pump'];rp,xp=data['plus'];pre=r[:,0]<=probe
 summary['checks']=dict(pre_probe_current_difference=float(np.max(abs(rp[:len(r)][pre,15]-r[pre,15]))),
   max_external_field_after_pulse=float(abs(r[r[:,0]>60.08,1:7]).max()),
   available_post_probe_before_abort_fs=(summary['pump']['failure_time_fs']-probe*fs),
   clean_post_probe_before_norm_threshold_fs=summary['pump']['first_logged_norm_error_over_1e_4_fs']-probe*fs)
 clean_time=summary['checks']['clean_post_probe_before_norm_threshold_fs']
 summary['checks']['approximate_clean_window_resolution_eV']=4.135667696/clean_time
 # Compare only through 6.0375fs, before the coarse run's recorded norm failure.
 n=3120;ref=r[:n];refx=x[:n]
 for name,stride in [('diagnostic_halfdt',2),('diagnostic_stride1',1)]:
  a,b=data[name];a=a[stride-1::stride];b=b[stride-1::stride]
  np.testing.assert_allclose(a[:,0],ref[:,0],atol=1e-8,rtol=0)
  summary[name]['comparison_to_dt008_stride10']=dict(
   end_time_fs=float(a[-1,0]*fs),current_relative_l2=float(np.linalg.norm(a[:,15]-ref[:,15])/np.linalg.norm(ref[:,15])),
   Axc_relative_l2=float(np.linalg.norm(b[:,3]-refx[:,3])/np.linalg.norm(refx[:,3])),max_alpha_difference=float(abs(b[:,7]-refx[:,7]).max()),
   final_Axc=float(b[-1,3]),reference_final_Axc=float(refx[-1,3]),final_current=float(a[-1,15]),reference_final_current=float(ref[-1,15]))
 # A valid equilibrium reference is useful; no failed pumped transform is made.
 gr,gx=data['ground'];assert summary['ground']['valid_trajectory']
 assert abs(gr[gr[:,0]<=probe,15]).max()<1e-9
 eps=response(gr[:,0],gr[:,15],eta,energy,probe)
 m=peak_metrics(energy,eps.imag,2,4);m['apparent_FWHM_eV']=apparent_width(energy,eps.imag);m['windows']={}
 for span in (600.,720.,880.):
  keep=gr[:,0]<=probe+span+1e-8;ee=response(gr[keep,0],gr[keep,15],eta,energy,probe)
  m['windows'][str(span)]=peak_metrics(energy,ee.imag,2,4)
 summary['ground_spectrum']=m
 np.savetxt(out/'ground_spectrum.csv',np.column_stack((energy,eps.real,eps.imag)),delimiter=',',header='energy_eV,Re_epsilon,Im_epsilon',comments='')
 summary['pumped_spectrum']=dict(valid=False,reason='Both pump-only and pump+probe fail norm before the required observation window; no pumped peak extracted.')
 for time in (3.,4.,5.,6.):
  i=np.argmin(abs(r[:,0]*fs-time));summary[f'pump_at_{time}fs']=dict(time_fs=float(r[i,0]*fs),alpha=float(x[i,7]),Axc=float(x[i,3]),current=float(r[i,15]),Exc=float(x[i,6]),P=float(x[i,11]))
 (out/'diagnostic_metrics.json').write_text(json.dumps(summary,indent=2)+'\n')
 fig,ax=plt.subplots(4,1,figsize=(9,10),sharex=True,layout='constrained')
 styles={'pump':('tab:red','-','Pump, dt=.08, ELF stride10'), 'plus':('tab:orange','--','Pump + probe'),
  'diagnostic_halfdt':('tab:blue',':','Pump, dt=.04, ELF stride20'), 'diagnostic_stride1':('tab:green','-.','Pump, dt=.08, ELF stride1')}
 for name,(color,style,label) in styles.items():
  a,b=data[name];t=a[:,0]*fs
  for axis,col in [(ax[0],7),(ax[1],3)]:axis.plot(t,b[:,col],color=color,ls=style,label=label)
  ax[2].plot(t,a[:,15],color=color,ls=style)
  z=np.load(out/f'{name}_trace.npz');ax[3].semilogy(z['norm_time_fs'],np.maximum(abs(z['electron_count']-32),1e-10),color=color,ls=style)
 ax[0].set(ylabel='alpha');ax[0].legend(fontsize=8);ax[1].set(ylabel='Axc,z / c (a.u.)');ax[2].set(ylabel='Jz (a.u.)');ax[3].set(ylabel='|N - 32|',xlabel='Time (fs)')
 bad=summary['pump']['first_logged_norm_error_over_1e_4_fs']
 for a in ax:
  a.axvline(60*fs,color='0.5',ls=':',lw=1);a.axvline(probe*fs,color='0.5',ls='--',lw=1)
  a.axvspan(bad,6.55,color='0.5',alpha=.15);a.set_xlim(0,6.55);a.grid(alpha=.2)
 ax[2].set_ylim(-.13,.004);ax[3].axhline(1e-4,color='0.5',ls=':')
 fig.suptitle('Si 4x4x4 k: strong pump cannot reach the spectral observation window\nShading: electron-count error > 1e-4; late alpha is not physically reliable')
 fig.savefig(out/'instability.png',dpi=160);fig.savefig(out/'instability.pdf')
 fig,ax=plt.subplots(figsize=(8,4),layout='constrained');ax.plot(energy,eps.imag)
 ax.set(xlim=(2,4.5),xlabel='Energy (eV)',ylabel='Im epsilon',title='Unpumped Si: dynamic ELF alpha0=0.2\nValid reference only; no pumped peak extracted');ax.grid(alpha=.2)
 fig.savefig(out/'ground_spectrum.png',dpi=160)
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
