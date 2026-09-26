"""Independent eigen-residual and periodic axial-tail diagnostics (atomic units)."""
from pathlib import Path
import numpy as np


def eigen_diagnostics(psi,hpsi,dv,occupation):
    psi=np.asarray(psi);hpsi=np.asarray(hpsi)
    norm=np.sum(abs(psi)**2,axis=0)*dv
    eigen=np.sum(psi.conj()*hpsi,axis=0)*dv/norm
    residual=np.sqrt(np.sum(abs(hpsi-psi*eigen)**2,axis=0)*dv/norm)
    overlap=psi.conj().T@psi*dv
    return dict(max_residual_Ha=float(residual.max()),
                occupied_rms_residual_Ha=float(np.sqrt(np.sum(occupation*residual**2)/np.sum(occupation))),
                max_orthogonality_error=float(np.max(abs(overlap-np.eye(len(norm))))),
                max_rayleigh_imaginary_Ha=float(np.max(abs(eigen.imag))),
                eigenvalues_Ha=eigen.real.tolist(),residuals_Ha=residual.tolist())


def read_eigen_pair(path):
    with Path(path).open('rb') as f:
        magic=f.read(4);endian='<' if magic==b'\x04\x03\x02\x01' else '>'
        f.seek(0);header=np.fromfile(f,endian+'i4',5)
        if len(header)!=5 or header[0]!=16909060 or header[1]!=1 or header[4]!=1:
            raise ValueError('invalid Gamma eigen diagnostic header')
        ng,no=map(int,header[2:4]);dv=np.fromfile(f,endian+'f8',1)[0]
        occ=np.fromfile(f,endian+'f8',no)
        psi=np.fromfile(f,endian+'c16',ng*no).reshape((ng,no),order='F')
        hp=np.fromfile(f,endian+'c16',ng*no).reshape((ng,no),order='F')
        if f.read(1):raise ValueError('trailing data')
    return psi,hp,dv,occ


def axial_tail(q,spacing,radii):
    """Integrate over y,z; radius is a minimum-image x half-width, not sphere."""
    q=np.asarray(q);spacing=np.asarray(spacing);radii=np.asarray(radii)
    weight=np.sum(abs(q)**2,axis=(2,3));norm=weight.sum(axis=1)
    x=np.arange(q.shape[1])*spacing[0];length=q.shape[1]*spacing[0]
    moment=weight@np.exp(2j*np.pi*x/length)
    center=(np.angle(moment)%(2*np.pi))*length/(2*np.pi)
    distance=abs((x[None,:]-center[:,None]+length/2)%length-length/2)
    tail=np.array([np.sum(weight*(distance>r),axis=1)/norm for r in radii])
    return dict(center_bohr=center,center_reliability=abs(moment)/norm,
                tail_fraction=tail,radii_bohr=radii,axial_rms_bohr=np.sqrt(np.sum(weight*distance**2,axis=1)/norm))


def fock_action(source,target,spacing,omega,batch=8):
    """Full periodic Gamma kernel, no alpha. Each source defines one common operator."""
    axes=[2*np.pi*np.fft.fftfreq(n,d=h) for n,h in zip(source.shape[1:],spacing)]
    g2=sum(g*g for g in np.meshgrid(*axes,indexing='ij'))
    kernel=np.full(source.shape[1:],np.pi/omega**2)
    np.divide(4*np.pi*(-np.expm1(-g2/(4*omega**2))),g2,out=kernel,where=g2>0)
    result=np.zeros_like(target,dtype=complex)
    for q in source:
        for first in range(0,len(target),batch):
            sl=slice(first,first+batch)
            potential=np.fft.ifftn(np.fft.fftn(q.conj()*target[sl],axes=(-3,-2,-1))*kernel,axes=(-3,-2,-1))
            result[sl]-=q*potential
    return result


def support_sweep(q,spacing,omega,radii,natom,min_center_reliability=0.):
    """Frozen density diagnostic. No renormalization and no SCF certification."""
    q=np.asarray(q,dtype=complex);spacing=np.asarray(spacing,dtype=float)
    radii=np.asarray(radii,dtype=float)
    if q.ndim!=4 or not np.isfinite(q).all() or np.any(np.sum(abs(q)**2,axis=(1,2,3))==0):
        raise ValueError('finite nonzero factors required')
    if spacing.shape!=(3,) or np.any(spacing<=0) or not np.isfinite(spacing).all():
        raise ValueError('positive finite spacing required')
    if omega<=0 or not np.isfinite(omega) or natom<1 or not np.isfinite(radii).all() or np.any(radii<0):
        raise ValueError('invalid sweep controls')
    if not np.isfinite(min_center_reliability) or not 0<=min_center_reliability<=1:
        raise ValueError('center reliability threshold must lie in [0,1]')
    tail=axial_tail(q,spacing,radii)
    protected=tail['center_reliability']<min_center_reliability
    x=np.arange(q.shape[1])*spacing[0];length=q.shape[1]*spacing[0]
    distance=abs((x[None,:]-tail['center_bohr'][:,None]+length/2)%length-length/2)
    dv=float(np.prod(spacing));alpha=.25
    reference=fock_action(q,q,spacing,omega)
    energy=alpha*dv*np.vdot(q,reference).real
    reports=[]
    for i,radius in enumerate(radii):
        mask=(distance<=radius)|protected[:,None]
        candidate=q if radius>=length/2 else q*mask[:,:,None,None]
        if np.array_equal(candidate,q):
            action=reference;self_action=reference
        else:
            action=fock_action(candidate,q,spacing,omega)
            self_action=fock_action(candidate,candidate,spacing,omega)
        expectation=alpha*dv*np.vdot(q,action).real
        self_energy=alpha*dv*np.vdot(candidate,self_action).real
        metric=dv*q.reshape(len(q),-1).conj()@action.reshape(len(q),-1).T
        reports.append(dict(radius_bohr=float(radius),radius_angstrom=float(radius*.529177210903),
            center_reliability_threshold=float(min_center_reliability),
            full_support_factor_count=int(protected.sum()),
            minimum_center_reliability=float(tail['center_reliability'].min()),
            max_truncated_factor_tail_fraction=float(tail['tail_fraction'][i,~protected].max()) if np.any(~protected) else 0.,
            max_tail_fraction=float(tail['tail_fraction'][i].max()),
            discarded_norm_fraction=float(1-np.vdot(candidate,candidate).real/np.vdot(q,q).real),
            full_exchange_Ha=float(energy),truncated_density_exchange_Ha=float(self_energy),
            exchange_error_Ha=float(self_energy-energy),exchange_error_meV_per_atom=float((self_energy-energy)*27211.386245988/natom),
            original_density_exchange_expectation_Ha=float(expectation),
            relative_action_error=float(np.linalg.norm(action-reference)/np.linalg.norm(reference)),
            target_metric_antihermitian_relative=float(np.linalg.norm(metric-metric.conj().T)/max(np.linalg.norm(metric),1e-300)),
            certified_for_scf=False))
    return reports


def main():
    import argparse,json
    from read_snapshot import read_snapshot
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('snapshot',type=Path);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--radii',nargs='+',type=float,required=True)
    p.add_argument('--natom',type=int,required=True);p.add_argument('--allow-unconverged',action='store_true')
    p.add_argument('--min-center-reliability',type=float,default=0.)
    a=p.parse_args();s=read_snapshot(a.snapshot)
    if np.prod(s['mesh'])!=1:raise ValueError('Gamma required')
    if (not s['converged'] or s['localization_status']!=0) and not a.allow_unconverged:
        raise ValueError('Both SCF and localization convergence required for reference sweep')
    q=s['q'].reshape(tuple(s['n'])+(s['q'].shape[1],),order='F').transpose(3,0,1,2)
    result=dict(scope='Frozen whole-periodic-cell density; axial x support, no renormalization. Not a total-SCF-energy or force certificate.',
        snapshot=str(a.snapshot.resolve()),scf_converged=s['converged'],localization_status=s['localization_status'],
        reports=support_sweep(q,s['spacing'],s['omega'],a.radii,a.natom,a.min_center_reliability))
    a.output.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
