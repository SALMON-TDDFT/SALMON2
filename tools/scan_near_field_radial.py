"""Saved-field diagnostic: local volumes, held-out points, radial-complete solid harmonics vs Cartesian polynomials.

Polynomial fits are an error diagnostic, not a new TDDFT input representation.
Numerical fits and residuals are evaluated by Fortran. Never starts a solver.
"""
import argparse
import ctypes as c
import json
import hashlib
from pathlib import Path
from unittest.mock import patch
import warnings
import numpy as np
from salmonpy import _material_handles
from salmonpy.field_sampling import SampledField
from salmonpy.near_field import NearFieldSeries


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('samples',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--kernels',type=Path,required=True)
    p.add_argument('--polynomial-library',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    lib=c.CDLL(str(a.polynomial_library))
    poly=lib.sp_polynomial_probe
    poly.argtypes=[c.c_int]*3+[c.c_void_p]*6
    poly.restype=None
    def polynomial(points,values,probes,degree):
        points=np.ascontiguousarray(points);values=np.ascontiguousarray(values)
        probes=np.ascontiguousarray(probes);out=np.empty(len(probes))
        rank=c.c_int();status=c.c_int()
        poly(len(points),len(probes),degree,points.ctypes.data,values.ctypes.data,probes.ctypes.data,
             out.ctypes.data,c.byref(rank),c.byref(status))
        if status.value:raise ValueError(f'polynomial rank {rank.value}, status {status.value}')
        return out
    rows=[]
    with patch.object(_material_handles,'_library_path',return_value=a.kernels):
        api=_material_handles.default_sbe_api()
        compare=api.library.sp_convergence_waveforms
        compare.argtypes=[c.c_int,c.c_void_p,c.c_void_p,c.c_double,c.c_double,c.c_void_p,c.c_void_p]
        compare.restype=c.c_int
        def metric(pred,ref):
            pred=np.ascontiguousarray(pred);ref=np.ascontiguousarray(ref);stats=np.empty(4);ok=c.c_int()
            if compare(ref.size,pred.ctypes.data,ref.ctypes.data,0.,.05,stats.ctypes.data,c.byref(ok)):
                raise RuntimeError('Fortran comparison failed')
            return dict(zip(('absolute_rms','relative_rms','reference_rms','maximum_absolute'),map(float,stats)))
        # Independent oracle: missing r^2 term, exactly recovered by complete degree2.
        rng=np.random.default_rng(991)
        x=rng.uniform(-1,1,(180,3));q=rng.uniform(-1,1,(90,3))
        v=np.einsum('ij,ij->i',x,x);truth=np.einsum('ij,ij->i',q,q)
        probe=metric(polynomial(x,v,q,2),truth)
        if probe['maximum_absolute']>1.e-12:raise RuntimeError('radial polynomial oracle failed')
        samples=SampledField('E',[0.,.02],x,np.tile(v[None,:,None],(2,1,3)),'E:au')
        zero=SampledField('H',[0.,.02],x,np.zeros((2,len(x),3)),'H:au')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',UserWarning)
            harmonic=NearFieldSeries.from_samples(samples,zero,center=(0,0,0),radius=1.,lmax=2)
        radial_missing=metric(harmonic.fields(q,step_index=0)[0][:,0],truth)
        if radial_missing['relative_rms']>1.e-12:raise RuntimeError('radial-complete solid-harmonic oracle failed')
        with np.load(a.samples,allow_pickle=False) as data:
            positions=data['positions'];center=np.array([5.17,4.77,5.31])
            for halfwidth in (2.,4.,6.,8.):
                selected=np.all(np.abs(positions-[5.,5.,5.])<=halfwidth+1.e-10,axis=1)
                local=positions[selected];held=np.random.default_rng(729).permutation(len(local))
                validation=held[:max(1,len(local)//4)];training=held[len(validation):]
                # fixed normalization and split across all orders/bases for this volume
                scale=halfwidth*2.;train=(local[training]-center)/scale;test=(local[validation]-center)/scale
                for degree in (1,2,3,4,6,8):
                    for step in (40,71,82):
                        E=data['E'][step,selected];H=data['H'][step,selected]
                        row=dict(halfwidth_bohr=halfwidth,side_angstrom=2*halfwidth*.529177210903,
                                 lmax=degree,step=step,training_points=len(training),held_out_points=len(validation))
                        for basis in ('radial_solid_harmonic','complete_polynomial'):
                            try:
                                if basis=='radial_solid_harmonic':
                                    with warnings.catch_warnings():
                                        warnings.simplefilter('ignore',UserWarning)
                                        fitted=NearFieldSeries.from_samples(
                                            SampledField('E',[0.,.02],local[training],[E[training],E[training]],'E:au'),
                                            SampledField('H',[0.,.02],local[training],[H[training],H[training]],'H:au'),
                                            center=center,radius=scale,lmax=degree)
                                    predicted=fitted.fields(local[validation],step_index=0)
                                else:
                                    predicted=[np.stack([polynomial(train,f[training,k],test,degree) for k in range(3)],axis=1) for f in (E,H)]
                                row[basis]={name:metric(pred,ref[validation]) for name,pred,ref in zip(('E','H'),predicted,(E,H))}
                            except (ValueError,RuntimeError) as exc:row[basis]={'unavailable':str(exc)}
                        rows.append(row)
    a.output.write_text(json.dumps(dict(schema='salmonpy.radial-basis-probe/2',source=str(a.samples),
        center_bohr=center.tolist(),steps=[40,71,82],radial_oracle=probe,solid_radial_oracle=radial_missing,
        source_sha256=hashlib.sha256(a.samples.read_bytes()).hexdigest(),
        kernels_sha256=hashlib.sha256(a.kernels.read_bytes()).hexdigest(),
        polynomial_library_sha256=hashlib.sha256(a.polynomial_library.read_bytes()).hexdigest(),
        polynomial_source_sha256=hashlib.sha256(Path(__file__).with_name('near_field_polynomial_probe.f90').read_bytes()).hexdigest(),
        scope='LD sphere transient; not Tip/Bow Tie convergence; diagnostic polynomial is not a TDDFT coupling',rows=rows),indent=2)+'\n')
    for width in (2.,4.,6.,8.):
        print('halfwidth',width)
        for row in rows:
            if row['halfwidth_bohr']==width and row['step']==71:
                print(row['lmax'],*[f"{basis}: "+str(row[basis].get('E',row[basis])) for basis in ('radial_solid_harmonic','complete_polynomial')])


if __name__=='__main__':main()
