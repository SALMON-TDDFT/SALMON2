"""Offline reduced-W/RT comparison. Screening and propagation stay in Fortran."""
from pathlib import Path
import argparse,hashlib,json,re,os
import numpy as np
EV=27.211386245988
FS=.024188843265857

def read_w(path):
    data=path.read_bytes();dims=np.frombuffer(data,'=i4',count=9,offset=4)
    assert data[:4]==b'HXT1' and len(dims)==9
    ng=int(np.prod(dims[:3]));nk=int(np.prod(dims[3:6]));na=int(dims[6])**3
    end=40+8*(9+3*nk+3*int(dims[8])+ng+1+3*nk)
    w=np.frombuffer(data,'=f8',offset=end)
    assert len(w)==nk*na*na and np.isfinite(w).all()
    return data[:end],w,hashlib.sha256(data).hexdigest()

def trajectory(path,steps=None):
    data=np.loadtxt(path)
    assert data.ndim==2 and data.shape[1]==16 and np.isfinite(data).all()
    if steps is not None:assert len(data)==steps,(path,len(data),steps)
    assert np.max(abs(np.diff(data[:,0])-.05))<1e-10
    return data

def epsilon(time,current,energies,kick=1e-4):
    window=1-3*(time/time[-1])**2+2*(time/time[-1])**3
    frequencies=energies/EV
    sigma=np.array([np.sum(current*window*np.exp(1j*w*time))*.05/kick for w in frequencies])
    return 1+4j*np.pi*sigma/frequencies

def peaks(energies,eps):
    maxima=np.flatnonzero((eps.imag[1:-1]>eps.imag[:-2])&(eps.imag[1:-1]>=eps.imag[2:]))+1
    return [dict(energy_eV=float(energies[i]),height=float(eps.imag[i])) for i in maxima
            if 1<=energies[i]<=6 and eps.imag[i]>0]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['screen','replay','spectra'])
    parser.add_argument('root',type=Path)
    parser.add_argument('--reference-kick',type=Path,default=Path('/private/tmp/si-lda-w-k4-20261121/tbw/Si_rt.data'))
    parser.add_argument('--reference-zero',type=Path,default=Path('/private/tmp/si-lda-w-k4-20261121/tbw0/Si_rt.data'))
    args=parser.parse_args();root=args.root;out=root/'analysis';out.mkdir(exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR',str(out/'mpl-cache'))
    if args.mode=='replay':
        report={};passed=True
        for stage,old in [('reference-replay',args.reference_kick),('reference-zero-replay',args.reference_zero)]:
            assert 'end SALMON' in (root/stage/'run.log').read_text()
            new=trajectory(root/stage/'Si_rt.data',1024);previous=trajectory(old)[:len(new)]
            np.testing.assert_array_equal(new[:,0],previous[:,0])
            error=float(np.max(abs(new[:,13:16]-previous[:,13:16])))
            report[stage]=dict(steps=len(new),all_J_max_absolute_difference=error,source=str(old))
            passed &= error<1e-11
        _,new_w,_=read_w(root/'gs-w/hsex-reference.bin')
        _,old_w,_=read_w(args.reference_kick.parent/'hsex-w.bin')
        report['reference_W_relative_difference_from_old']=float(np.linalg.norm(new_w-old_w)/np.linalg.norm(old_w))
        passed &= report['reference_W_relative_difference_from_old']<1e-8
        report['passed']=bool(passed)
        (out/'reference-replay.json').write_text(json.dumps(report,indent=2)+'\n')
        # Independent check of the transform against SALMON's own short-run file.
        emitted=np.loadtxt(root/'reference-replay/Si_response.data')
        raw=trajectory(root/'reference-replay/Si_rt.data',1024)
        wanted=epsilon(raw[:,0],raw[:,15],emitted[:,0]*EV)
        actual=emitted[:,9]+1j*emitted[:,12]
        np.testing.assert_allclose(wanted,actual,atol=1e-8,rtol=1e-10)
        report['Fourier_max_absolute_difference_from_SALMON']=float(np.max(abs(wanted-actual)))
        (out/'reference-replay.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report,indent=2));return
    text=(root/'gs-w/run.log').read_text()
    assert 'end SALMON' in text and 'HSEX_REDUCED_REFERENCE_W_PASS' in text
    metadata,a,sha_a=read_w(root/'gs-w/hsex-reference.bin');metadata_b,b,sha_b=read_w(root/'gs-w/hsex-w.bin')
    assert metadata==metadata_b,'Reference and approximate W must use identical GS/H metadata'
    tensors=np.array([list(map(float,l.split()[-3:])) for l in text.splitlines()
                      if l.startswith('HSEX_W_EPS_LONGWAVE diagonal')])
    assert tensors.shape==(2,3)
    np.testing.assert_allclose(tensors[0],tensors[1],atol=1e-9,rtol=1e-10)
    timings=[list(map(float,l.split(':')[1].split())) for l in text.splitlines()
             if l.startswith('HSEX_W_SECONDS response/assembly/transform:')]
    rows=[list(map(float,l.split()[-9:])) for l in text.splitlines() if l.startswith('HSEX_REDUCED_W_K ')]
    assert len(rows)==64 and sorted(int(r[1]) for r in rows)==list(range(1,65))
    screen=dict(reference_W_sha256=sha_a,reduced_W_sha256=sha_b,identical_GS_metadata=True,
                W_relative_Frobenius_difference=float(np.linalg.norm(b-a)/np.linalg.norm(a)),
                W_max_absolute_difference=float(np.max(abs(b-a))),
                exact_optical_epsilon_diagonal=tensors.tolist(),W_seconds=timings,
                unfinished_batches=int(sum(r[6] for r in rows)),total_batches=int(sum(r[5] for r in rows)),
                max_orbital_response_residual=max(r[8] for r in rows),
                direct_gap_eV=float(re.search(r'BG between same k-point\[eV\]\s+(\S+)',text)[1]))
    (out/'screening.json').write_text(json.dumps(screen,indent=2)+'\n')
    if args.mode=='screen':print(json.dumps(screen,indent=2));return
    for stage in ['reduced','reduced0']:
        assert 'end SALMON' in (root/stage/'run.log').read_text(),stage+' unfinished'
    new=trajectory(root/'reduced/Si_rt.data',16384);zero=trajectory(root/'reduced0/Si_rt.data',16384)
    old=trajectory(args.reference_kick,16384);oldzero=trajectory(args.reference_zero,16384)
    for data in [zero,old,oldzero]:np.testing.assert_array_equal(new[:,0],data[:,0])
    reused=not args.reference_kick.is_relative_to(root)
    if reused:assert json.loads((out/'reference-replay.json').read_text())['passed']
    time=new[:,0];currents=[old[:,15]-oldzero[:,15],new[:,15]-zero[:,15]]
    energies=np.arange(.005,8.0001,.005);visible=(energies>=1)&(energies<=6)
    spectra=[epsilon(time,j,energies) for j in currents]
    report=dict(screening=screen,steps=len(time),duration_fs=float(time[-1]*FS),
                nominal_fourier_resolution_eV=float(2*np.pi/time[-1]*EV),
                window='SALMON cubic: 1-3(t/T)^2+2(t/T)^3; sigma0 correction off; zero-field subtracted',
                reference_reused=reused,reference_paths=[str(args.reference_kick),str(args.reference_zero)],
                current_relative_L2_difference=float(np.linalg.norm(currents[1]-currents[0])/np.linalg.norm(currents[0])),
                current_max_absolute_difference=float(np.max(abs(currents[1]-currents[0]))),
                zero_field_max_absolute_Jz=[float(np.max(abs(oldzero[:,15]))),float(np.max(abs(zero[:,15])))],
                Im_relative_L2_difference_1_6_eV=float(np.linalg.norm((spectra[1].imag-spectra[0].imag)[visible])/np.linalg.norm(spectra[0].imag[visible])),
                Re_relative_L2_difference_1_6_eV=float(np.linalg.norm((spectra[1].real-spectra[0].real)[visible])/np.linalg.norm(spectra[0].real[visible])),cases={})
    for label,eps,current in zip(['reference','reduced'],spectra,currents):
        p=peaks(energies,eps);n=len(time)*3//4;p75=peaks(energies,epsilon(time[:n],current[:n],energies))
        stable=[v for v in p if v['height']>=.05*max(eps.imag[visible]) and any(abs(v['energy_eV']-w['energy_eV'])<=.1 for w in p75)]
        report['cases'][label]=dict(first_stable_peak=stable[0] if stable else None,
                                   peaks_1_6_eV=p,peaks_75percent_window=p75)
    ref=report['cases']['reference']['first_stable_peak'];red=report['cases']['reduced']['first_stable_peak']
    if ref is not None and red is not None:
        report['first_peak_energy_difference_eV']=red['energy_eV']-ref['energy_eV']
        report['first_peak_height_relative_difference']=red['height']/ref['height']-1
        report['comparison_within_declared_tolerance']=bool(abs(report['first_peak_energy_difference_eV'])<=.01 and abs(report['first_peak_height_relative_difference'])<=.01 and report['Im_relative_L2_difference_1_6_eV']<=.01)
    else:report['comparison_within_declared_tolerance']=False
    report['comparison_tolerance']='Same finite window: first peak <=0.01eV and height <=1%, Im L2 <=1%; not a physical convergence bound'
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.2))
    for eps,label,color in zip(spectra,['Reference W','Reduced W: rank384, symmetry'],['#333333','#d55e00']):
        style='--' if label.startswith('Reduced') else '-'
        axes[0].plot(energies,eps.real,label=label,color=color,ls=style,lw=1.5)
        axes[1].plot(energies,eps.imag,label=label,color=color,ls=style,lw=1.5)
    for ax,y in zip(axes,['Re epsilon zz','Im epsilon zz']):
        ax.set(xlim=(1,6),xlabel='Energy (eV)',ylabel=y);ax.grid(alpha=.25)
        ax.axvline(screen['direct_gap_eV'],color='.6',ls=':',lw=.8)
    allreal=np.r_[spectra[0].real[visible],spectra[1].real[visible]]
    pad=.05*np.ptp(allreal);axes[0].set_ylim(allreal.min()-pad,allreal.max()+pad)
    axes[1].set_ylim(0,1.06*max(max(z.imag[visible]) for z in spectra))
    axes[0].legend(fontsize=9);fig.suptitle('Si Nr12^3, Nk=Nq4^3: fixed W, occupied-only TB-mBJ RT\nSame 19.82 fs cubic window; zero-field subtracted',fontsize=11)
    fig.tight_layout();fig.savefig(out/'reference-vs-reduced-spectrum.png',dpi=180);fig.savefig(out/'reference-vs-reduced-spectrum.pdf')
    fig,axes=plt.subplots(2,1,figsize=(9,6),sharex=True)
    for j,label,color in zip(currents,['Reference W','Reduced W'],['#333333','#d55e00']):
        axes[0].plot(time*FS,j,label=label,color=color,ls='--' if label.startswith('Reduced') else '-')
    axes[1].plot(time*FS,currents[1]-currents[0],color='#d55e00');axes[0].legend()
    axes[0].set_ylabel('Jz (a.u.)');axes[1].set_ylabel('Delta Jz (a.u.)');axes[1].set_xlabel('Time (fs)')
    for ax in axes:ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(out/'reference-vs-reduced-current.png',dpi=180)
    np.savetxt(out/'spectra.tsv',np.column_stack([energies,spectra[0].real,spectra[0].imag,spectra[1].real,spectra[1].imag]),header='eV Re_reference Im_reference Re_reduced Im_reduced')
    np.savetxt(out/'currents.tsv',np.column_stack([time*FS,*currents,currents[1]-currents[0],oldzero[:,15],zero[:,15]]),header='fs J_reference J_reduced delta_J zero_reference zero_reduced')
    (out/'metrics.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
