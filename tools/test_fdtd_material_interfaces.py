"""Native material-label invariance and singleton compatibility.
Usage: PYTHONPATH=. python tools/test_fdtd_material_interfaces.py BEFORE AFTER NEW_OUTPUT
Requires audited binaries. One thread, finite synthetic structures, no R/T claim.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from salmonpy.fdtd import FDTD, Dielectric, PlaneWave, Impulse, SalmonExecutable
from salmonpy.geometry import Grid, Geometry, Shape
from salmonpy.fields import AcosPulse
from salmonpy.ld import LorentzDrude

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('before');p.add_argument('after');p.add_argument('output')
p.add_argument('--quick',action='store_true',help='only two-layer ID permutation (red test)')
a=p.parse_args();out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False)
engines={k:SalmonExecutable(v,hashlib.sha256(Path(v).read_bytes()).hexdigest()) for k,v in [('before',a.before),('after',a.after)]}
ld1=LorentzDrude(plasma_frequency=20,strength=[1],resonance=[0],damping=[2])
ld2=LorentzDrude(plasma_frequency=15,strength=[.7,.3],resonance=[25,35],damping=[2,3])
grid=Grid((36,36,36),(2,2,2))
report={}
def make(materials,axis=0,swap=False,single=False,impulse=False):
    size=[16.,16.,16.];size[axis]=4.
    left=[0.,0.,0.];right=list(left);left[axis]=-2.;right[axis]=2.
    if single:
        shape=Shape('box',(8.,8.,8.));shapes=(shape,)
    else:
        shapes=(Shape('box',size,left,2 if swap else 1),Shape('box',size,right,1 if swap else 2))
    pol=[0.,0.,0.];pol[(axis+2)%3]=1.
    source=Impulse(polarization=pol) if impulse else PlaneWave(AcosPulse(.001,20.,.15,polarization=pol),axis,-14.)
    observation=[0.,0.,0.];observation[axis]=12.
    return FDTD(Geometry(grid,shapes),materials,.003,160,source=source,
                observations=(tuple(observation),(0.,0.,0.)),name='interface')
def run(label,calc,engine='after'):
    start=time.monotonic();r=engines[engine].run(calc,out/label,threads=1,timeout=300)
    report[label]={'seconds':time.monotonic()-start}
    return r

def equal(left,right,label):
    for i,(x,y) in enumerate(zip(left.observations,right.observations),1):
        diff=float(np.max(abs(x.values-y.values)))
        report[label+f'_obs{i}']={'max_abs':diff,'bitwise_equal':bool(np.array_equal(x.values,y.values))}
        np.testing.assert_array_equal(x.values,y.values,err_msg=label)
    print('PASS',label,flush=True)
try:
    axes=(0,) if a.quick else (0,1,2)
    profiles={'ld_ld':(ld1,ld2)} if a.quick else {
        'ld_ld':(ld1,ld2), 'ld_dielectric':(ld1,Dielectric(2.3,1.2,.02)),
        'magnetic_dielectric':(Dielectric(2,1.1,.01),Dielectric(3,1.8,.02))}
    for name,materials in profiles.items():
        for axis in axes:
            label=f'{name}_{axis}'
            one=run(label,make(materials,axis))
            two=run(label+'_relabel',make(materials[::-1],axis,swap=True))
            equal(one,two,label)
    if not a.quick:
        from dataclasses import replace
        calc=make((ld1,ld2))
        equal(run('shape_forward',calc),
              run('shape_reverse',replace(calc,geometry=Geometry(grid,calc.geometry.shapes[::-1]))),
              'shape_order')
        duplicate=make((ld2,ld2))
        one_id=replace(duplicate,materials=(ld2,),geometry=Geometry(grid,tuple(
            Shape('box',s.size,s.center,1) for s in duplicate.geometry.shapes)))
        equal(run('duplicate_ids',duplicate),run('one_id',one_id),'duplicate_material')
        for name,materials in [('ld1',(ld1,)),('ld2',(ld2,)),('dielectric',(Dielectric(2.3,1.2,.02),))]:
            for impulse in ([False,True] if name.startswith('ld') else [False]):
                label=name+('_impulse' if impulse else '_pulse')
                calc=make(materials,single=True,impulse=impulse)
                equal(run(label+'_before',calc,'before'),run(label+'_after',calc),label)
        for axis in axes:
            one=run(f'impulse_{axis}',make((ld1,ld2),axis,impulse=True))
            two=run(f'impulse_{axis}_relabel',make((ld2,ld1),axis,swap=True,impulse=True))
            equal(one,two,f'impulse_{axis}')
finally:
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
