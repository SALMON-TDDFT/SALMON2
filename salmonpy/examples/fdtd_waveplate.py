"""Synthetic waveplate: preview by default; --run explicitly starts FDTD.

python -m salmonpy.examples.fdtd_waveplate BUILD_JSON NEW_DIRECTORY --plate quarter --run
Use --spacing 4 for a spatial refinement.
This is a lossless principal-axis validation slab, not a fitted real material.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from salmonpy.ld import TensorLorentzDrude
from salmonpy.fdtd import FDTD,PlaneWave,Boundaries,SalmonExecutable,FDTDMPI
from salmonpy.fields import AcosPulse
from salmonpy.geometry import Grid,Shape,Geometry
from salmonpy.subpixel import SubpixelSmoothing
from salmonpy.analysis import waveplate_diagnostic,read_table


def configuration(plate='quarter',spacing=8.,angle=0.,subpixel=False):
    if plate not in ('quarter','half'):raise ValueError('quarter or half required')
    if not np.isfinite(spacing) or spacing<=0:raise ValueError('positive spacing required')
    if not np.isfinite(angle):raise ValueError('finite rotation angle in radians required')
    cs,sn=np.cos(angle),np.sin(angle)
    rotation=np.array([[cs,-sn,0.],[sn,cs,0.],[0.,0.,1.]])
    thickness=432. if plate=='quarter' else 864.
    grid=Grid((4,4,2*int(np.ceil(2048./spacing))),(spacing,)*3)
    shape=Shape('box',(4*spacing,4*spacing,thickness))
    geometry=Geometry(grid,(shape,))
    material=TensorLorentzDrude(plasma_frequency=0.,strength=[np.zeros((3,3))],resonance=[.4],damping=[0.],epsilon_infinity=rotation@np.diag([2.25,4.,2.25])@rotation.T)
    dt=.0025*spacing
    calculation=FDTD(geometry,materials=(material,),dt=dt,steps=int(np.ceil(180./dt)),
        boundaries=Boundaries(('periodic',)*4+('pml','pml')),
        source=PlaneWave(AcosPulse(amplitude=.001,omega=1.,duration=60.,polarization=rotation@np.array([2**-.5,2**-.5,0.])),2,-1400.),
        observations=((0.,0.,1400.),),observation_stride=2,
        smoothing=SubpixelSmoothing() if subpixel else None)
    actual=float(np.count_nonzero(geometry.material_ids()[2,2,:])*spacing)
    return calculation,thickness if subpixel else actual


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build_record',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--plate',choices=['quarter','half'],default='quarter')
    p.add_argument('--angle-degrees',type=float,default=0.,help='rotate principal axes about lab z; incident principal components remain equal')
    p.add_argument('--spacing',type=float,default=8.)
    p.add_argument('--ranks',type=int,default=1)
    p.add_argument('--subpixel',action='store_true',help='LD-only interface preparation at the input pulse carrier')
    p.add_argument('--run',action='store_true',help='explicitly launch FDTD; omit for preview')
    a=p.parse_args();calculation,thickness=configuration(a.plate,a.spacing,np.deg2rad(a.angle_degrees),a.subpixel)
    decomposition=FDTDMPI((1,1,a.ranks))
    if not a.run:
        print(json.dumps(dict(mode='preview_only',grid=calculation.geometry.grid.counts,dt=calculation.dt,steps=calculation.steps,
            thickness_bohr=thickness,smoothing_reference_omega=calculation.smoothing_reference,assets=list(calculation.input_assets()),message='No solver launched. Use --run explicitly inside an allocated compute job.'),indent=2))
        return
    runner=SalmonExecutable.from_mpi_build(a.build_record)
    result=runner.run(calculation,a.output,mpi=decomposition,threads=1,timeout=300.)
    table=read_table(result.directory/'obs1_at_point_rt.data')
    electric=np.column_stack([table.column('E_'+axis) for axis in 'xyz'])
    report=waveplate_diagnostic(table.column('Time'),electric,omega=1.,indices=(1.5,2.),thickness=thickness,units='atomic_gaussian',angle=np.deg2rad(a.angle_degrees))
    report['phasor']=[[float(z.real),float(z.imag)] for z in report['phasor']]
    for key in ('stokes_normalized','analytic_stokes_normalized'):report[key]=report[key].tolist()
    report.update(plate=a.plate,thickness_bohr=thickness,synthetic_material=True,angle_degrees=a.angle_degrees,smoothing_reference_omega=calculation.smoothing_reference)
    (a.output/'polarization.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
