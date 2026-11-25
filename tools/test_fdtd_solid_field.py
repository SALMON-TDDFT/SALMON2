"""Collect collocated 3-D LD/FDTD E samples for explicit one-way Ylm fitting.

Run with two host MPI ranks: one field rank, one labels-only rank. No electrons,
GS, implicit solver launch, field feedback or per-step world collective.
"""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import ParallelPlan,ParallelContext,SolverGroup
from salmonpy.fdtd import FDTD,PlaneWave
from salmonpy.geometry import Grid,Geometry,Shape
from salmonpy.fields import AcosPulse
from salmonpy.ld import LorentzDrude
from salmonpy.coupled import NativeEHField


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--eh',action='store_true',help='also collect and Fortran-center H')
    p.add_argument('--kernels',type=Path,help='explicit checksum-recorded calculation library for H centering')
    p.add_argument('--vacuum',action='store_true',help='same source and sampling without the LD sphere')
    p.add_argument('--warmup',type=int,default=0,help='explicit field steps before the selected electron time window')
    p.add_argument('--sample-indices',type=int,nargs='+',default=[18,20,22],help='grid indices used on each sampling axis')
    p.add_argument('--grid-size',type=int,default=40)
    p.add_argument('--sample-x',type=int,nargs='+')
    p.add_argument('--sample-y',type=int,nargs='+')
    p.add_argument('--sample-z',type=int,nargs='+')
    a=p.parse_args();world=MPI.COMM_WORLD
    if a.warmup<0:raise ValueError('nonnegative warmup required')
    if world.size!=2:raise ValueError('this short probe uses exactly two host ranks')
    axes=[values if values is not None else a.sample_indices for values in (a.sample_x,a.sample_y,a.sample_z)]
    for values in axes:
        if len(set(values))!=len(values) or any(i<1 or i>a.grid_size-2 for i in values):
            raise ValueError('distinct interior sample indices required')
    points=dict(enumerate(itertools.product(*axes)))
    plan=ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('labels',(1,),point_ids=tuple(points))))
    error=None
    if world.rank==0:
        try:a.output.mkdir(parents=True,exist_ok=False)
        except Exception as exc:error=str(exc)
    error=world.bcast(error,root=0)
    if error:raise ValueError(error)
    with ParallelContext.create(world,plan) as context:
        if world.rank==0:
            try:
                material=LorentzDrude(plasma_frequency=2.,strength=[1.],resonance=[0.],damping=[.2])
                geometry=Geometry(Grid((a.grid_size,)*3,(2.,)*3),(Shape('ellipsoid',(6.,6.,6.),(-10.,0.,0.)),))
                if a.vacuum:geometry=Geometry(geometry.grid,())
                calculation=FDTD(geometry,() if a.vacuum else (material,),.003,82+a.warmup,
                    source=PlaneWave(AcosPulse(.01,20.,.15),0,-16.),observations=((1.,1.,1.),))
                frames=[];times=[];left_frames=[];right_frames=[];Hframes=[]
                if a.eh:
                    import ctypes as c
                    from salmonpy._material_handles import SBEAPI,default_sbe_api
                    api=default_sbe_api() if a.kernels is None else SBEAPI(a.kernels)
                    fn=api.library.sp_molecular_interval_field
                    fn.argtypes=[c.c_void_p]*3;fn.restype=c.c_int
                    def average(left,right):
                        values=np.empty_like(left)
                        for q in range(len(left)):
                            api.call('sp_molecular_interval_field',left[q].ctypes.data,right[q].ctypes.data,values[q].ctypes.data)
                        return values
                with NativeEHField(context,calculation,a.build,a.output/'field',points=points,grid=(1,1,1)) as field:
                    for _ in range(a.warmup):field.advance(np.zeros((len(points),3)))
                    if a.eh:
                        previous=field.sample_magnetic_field().values[0] if a.warmup else np.zeros((len(points),3))
                    for _ in range(82):
                        field.advance(np.zeros((len(points),3)))
                        sample=field.sample_fields()['E']
                        frames.append(sample.values[0]);times.append(sample.times[0])
                        if a.eh:
                            current=field.sample_magnetic_field().values[0]
                            left_frames.append(previous);right_frames.append(current);Hframes.append(average(previous,current))
                            previous=current
                    positions=sample.positions
                    # NativeEHField admits only zero-initialized plane-wave runs.
                    # Record that declared step-0 condition explicitly.
                    native_times=np.array(times)
                    if a.warmup:
                        time_origin=float(native_times[0]);times=np.arange(len(times))*.003
                    else:
                        time_origin=0.;frames=[np.zeros_like(frames[0]),*frames];times=[0.,*times]
                        native_times=np.array(times)
                    extra={}
                    if a.eh:
                        if not a.warmup:
                            initial=np.zeros_like(Hframes[0])
                            Hframes=[initial,*Hframes];left_frames=[initial,*left_frames];right_frames=[initial,*right_frames]
                        extra=dict(E=frames,H=Hframes,H_left=left_frames,H_right=right_frames,
                                   H_left_times=np.asarray(times)-.0015,H_right_times=np.asarray(times)+.0015)
                    np.savez(a.output/'samples.npz',times=times,positions=positions,values=frames,
                             native_times=native_times,time_origin=time_origin,**extra)
                    report=dict(schema='salmonpy.fdtd-solid-samples/1',steps=82,dt=.003,
                        field_ranks=1,point_count=len(points),H_included=a.eh,H_centering='Fortran half-step average' if a.eh else None,declared_initial_E_zero=a.warmup==0,time_origin=time_origin,warmup_steps=a.warmup,quantity='E',unit='E:au',position_unit='bohr',
                        max_abs_E=float(np.max(np.abs(frames))),field_library_sha256=field.manifest['library_sha256'],
                        source='plane wave',environment='vacuum' if a.vacuum else 'LD sphere',mode='one-way saved samples; no molecular radiation feedback')
                    regular=all(len(v)>=2 and v==list(range(v[0],v[-1]+1,v[1]-v[0])) and v[1]>v[0] for v in axes)
                    report.update(grid_counts=[len(v) for v in axes] if regular else None,
                        grid_origin=positions[0].tolist() if regular else None,
                        grid_spacing=[2.*(v[1]-v[0]) for v in axes] if regular else None,frames=len(times))
                    if report['max_abs_E']<=1.e-18:raise AssertionError('field did not reach the sampling region')
                    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
            except Exception as exc:error=f'{type(exc).__name__}: {exc}'
        error=world.bcast(error,root=0)
        if error:raise RuntimeError(error)
    if world.rank==0:print(f'PASS: 3D LD/FDTD, {len(points)} collocated sites, {len(times)} E/H frames')


if __name__=='__main__':main()
