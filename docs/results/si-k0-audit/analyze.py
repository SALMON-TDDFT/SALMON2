"""Audit K0 and field gating using existing fixed-alpha data; no new TDDFT trajectory."""
from pathlib import Path
import json,os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent
w=root/'calculations/si_tdcdft_k4/instant_smoke'
fs=.02418884326505;omega=.2;alpha0=.2;floor=2e-5;k0=.19591042214256474

def load(name):
    r=np.loadtxt(w/name/'Si_rt.data');x=np.loadtxt(w/name/'Si_rt_xc.data')
    np.testing.assert_allclose(r[:,0],x[:,0],atol=1e-8,rtol=0)
    a=r[:,1:4];dt=r[1,0]-r[0,0]
    # Known compact external pulse has a(0)=0 and remains zero after this file's last time.
    ext=np.vstack([np.zeros((1,3)),a,np.zeros((1,3))])
    e=-(ext[2:]-ext[:-2])/(2*dt)
    return r[:,0],a,e,r[:,13:16],x[:,9:12],x[:,8]

def estimate(a,e,j,p,threshold,reference):
    denom=np.sum(a*a+(e/omega)**2,axis=1);num=np.sum(a*j-e*p,axis=1)
    k=np.divide(num,denom,out=np.zeros_like(num),where=denom>0)
    valid=denom>threshold**2
    alphas=np.empty_like(k);held=alpha0
    for n in range(len(k)):
        if valid[n]:held=alpha0/(1+4*np.pi*max(k[n]-reference,0)/omega**2)
        alphas[n]=held
    return k,denom,num,valid,alphas

t,a,e,j,p,saved=load('weak_reference')
k,d,num,valid,weak_alpha=estimate(a,e,j,p,floor,k0)
np.testing.assert_allclose(k[valid],saved[valid],rtol=1e-10,atol=1e-12)
scale=np.sqrt(1e5)
ks,ds,ns,vs,scaled_alpha=estimate(a*scale,e*scale,j*scale,p*scale,floor,k0)
np.testing.assert_allclose(ks[d>0],k[d>0],rtol=1e-12,atol=1e-12)
assert weak_alpha.min()==alpha0 and scaled_alpha.min()<1e-4
# A common relative threshold is invariant under pure amplitude rescaling.
rel=.1
kr,dr,nr,vr,ar=estimate(a,e,j,p,rel*np.sqrt(d.max()),k0)
krs,drs,nrs,vrs,ars=estimate(a*scale,e*scale,j*scale,p*scale,rel*np.sqrt(ds.max()),k0)
np.testing.assert_array_equal(vr,vrs)
np.testing.assert_allclose(ar,ars,rtol=1e-12,atol=1e-12)
_,sa,se,sj,sp,_=load('strong_fixed')
np.testing.assert_allclose(sa,scale*a,rtol=1e-12,atol=1e-14)
sk,sd,sn,sv,_=estimate(sa,se,sj,sp,rel*np.sqrt(np.sum(sa*sa+(se/omega)**2,axis=1).max()),k0)
np.testing.assert_array_equal(vr,sv)
# Finite-window equilibrium optical reference, not a static or exact susceptibility.
spectrum=np.loadtxt(root/'docs/results/si-k4/lrc.csv',delimiter=',',skiprows=1)
energy=omega*27.21138505
re=np.interp(energy,spectrum[:,0],spectrum[:,1]);im=np.interp(energy,spectrum[:,0],spectrum[:,2])
optical_k=omega**2*(1-re)/(4*np.pi)
weighted_k=float(num[vr].sum()/d[vr].sum())
candidates={'old_maximum_absolute_gate':k0,'weak_field_weighted_fit_relative_gate':weighted_k,
            'equilibrium_optical_real_response':float(optical_k),
            'weak_maximum_relative_gate':float(k[vr].max())}
checks={}
for name,value in candidates.items():
    checks[name]=dict(K0=value,
        weak_min_alpha=float(estimate(a,e,j,p,rel*np.sqrt(d.max()),value)[-1].min()),
        strong_fixed_offline_min_alpha=float(estimate(sa,se,sj,sp,rel*np.sqrt(sd.max()),value)[-1].min()))
changed=np.flatnonzero(scaled_alpha<.1999)
metrics=dict(omega_au=omega,pulse_end_fs=60*fs,absolute_floor=floor,relative_floor=rel,
    weak_peak_field_norm=float(np.sqrt(d.max())),strong_equivalent_peak_field_norm=float(np.sqrt(ds.max())),
    weak_last_valid_fs=float(t[np.flatnonzero(valid)[-1]]*fs),
    scaled_linear_last_valid_fs=float(t[np.flatnonzero(vs)[-1]]*fs),
    weak_alpha_final=float(weak_alpha[-1]),scaled_linear_alpha_final=float(scaled_alpha[-1]),
    scaled_linear_first_correction_fs=float(t[changed[0]]*fs),
    relative_gate_alpha_max_rescaling_difference=float(np.max(abs(ar-ars))),
    weak_max_K_relative_gate=float(k[vr].max()),strong_max_K_relative_gate=float(sk[sv].max()),
    weak_weighted_K=weighted_k,strong_weighted_K=float(sn[sv].sum()/sd[sv].sum()),
    equilibrium_energy_eV=energy,equilibrium_Re_epsilon=float(re),equilibrium_Im_epsilon=float(im),
    candidate_reference_checks=checks)
(out/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
fig,ax=plt.subplots(2,1,figsize=(9,7),sharex=True,layout='constrained')
ax[0].plot(t[vr]*fs,k[vr],label='Weak fixed-alpha response')
ax[0].plot(t[sv]*fs,sk[sv],label='Strong fixed-alpha response')
ax[0].axhline(k0,color='0.5',ls='--',label='Original K0')
ax[0].set_ylabel('K (common relative gate 0.1)');ax[0].legend(fontsize=9)
ax[1].plot(t*fs,weak_alpha,label='Original weak data, absolute gate')
ax[1].plot(t*fs,scaled_alpha,label='Same linear data x316, absolute gate')
ax[1].plot(t*fs,ars,ls='--',label='Same linear data x316, relative gate 0.1')
ax[1].set_ylabel('Offline alpha estimate');ax[1].set_xlabel('Time (fs)');ax[1].legend(fontsize=9)
for z in ax:
    z.axvline(60*fs,color='0.6',ls=':');z.grid(alpha=.2);z.set_xlim(0,1.6)
fig.suptitle('K0 audit: amplitude scaling alone triggers the old correction\nNo nonlinear electronic evolution in the scaled-data test')
fig.savefig(out/'audit.png',dpi=160);fig.savefig(out/'audit.pdf')
print(json.dumps(metrics,indent=2))
