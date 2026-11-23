"""Explicit field-only constitutive energy and restart regression; no GS/SCF.

mpiexec -n 3 python tools/test_diagonal_field_energy_mpi.py BUILD NEW_OUTPUT \
    --field-ranks 2 --profile diagonal-ld --observe
Use --checkpoint NEW and --stop-after32, then --restore NEW in a new MPI job.
This development fixture stores direct numerical states; it is not a new API.
"""
import argparse
import ctypes as c
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time
import numpy as np
from salmonpy.distributed import SolverGroup,ParallelPlan,ParallelContext
from salmonpy.coupled import NativeEHField
from salmonpy.fdtd import FDTD,PlaneWave,Dielectric
from salmonpy.geometry import Grid,Geometry,Shape
from salmonpy.fields import AcosPulse
from salmonpy.ld import LorentzDrude,TensorLorentzDrude
from salmonpy.native_steps import _NativeBackup
from salmonpy._native_checkpoint import NativeRestoreCandidate
from salmonpy.api import _bind_field_operations


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--field-ranks',type=int,default=2)
    p.add_argument('--split-axis',type=int,choices=(0,1,2),default=0)
    p.add_argument('--profile',choices=('legacy','dielectric','scalar-ld','diagonal-ld','diagonal-unit-x','tensor-pole','tensor-background','subpixel','subpixel-aligned','ohmic','ohmic-unit','ohmic-ld'),default='diagonal-ld')
    p.add_argument('--owned-current',action='store_true')
    p.add_argument('--driven-current',action='store_true',help='manufactured constant j after step1; test only')
    p.add_argument('--execution-mode',choices=('checked','performance'),default='checked')
    p.add_argument('--grid-count',type=int,default=40)
    p.add_argument('--observation-stride',type=int,default=1)
    p.add_argument('--fail-resource-probe',action='store_true',help='explicit measurement-failure fixture; computation still runs')
    p.add_argument('--energy-lower',type=int,nargs=3)
    p.add_argument('--observe',action='store_true');p.add_argument('--checkpoint',type=Path)
    p.add_argument('--restore',type=Path);p.add_argument('--stop-after',type=int,default=64)
    a=p.parse_args()
    if a.observation_stride<1:p.error('observation-stride must be positive')
    if a.grid_count<40 or a.grid_count%2:p.error('grid-count must be even and at least40')
    if not 1<=a.stop_after<=64:p.error('stop-after must be in1..64')
    if a.field_ranks<1 or a.grid_count%a.field_ranks:p.error('field-ranks must divide grid-count')
    from mpi4py import MPI
    world=MPI.COMM_WORLD
    center=a.grid_count//2
    lower=tuple(a.energy_lower) if a.energy_lower is not None else (center-7,)*3
    upper=(center+7,)*3
    phases={name:0. for name in ('setup','advance','backup_capture','observation','backup_release','loop','checkpoint','close')}
    energy_queries=0;role=None
    setup_start=time.perf_counter()
    if world.size!=a.field_ranks+1:raise ValueError('field ranks plus one bookkeeping rank required')
    if world.rank==0:a.output.mkdir(parents=True,exist_ok=False)
    world.Barrier()
    plan=ParallelPlan((SolverGroup('field',tuple(range(a.field_ranks)),role='field'),
                      SolverGroup('bookkeeping',(a.field_ranks,),point_ids=(10,))))
    history=[];steps=0
    if a.restore:
        data=json.loads((a.restore/'complete.json').read_text())
        if data['profile']!=a.profile:raise ValueError('profile mismatch')
        steps=data['step'];history=data['history']
    with ParallelContext.create(world,plan) as context:
        field=None;role=context.group.role
        try:
            if context.group.role=='field':
                if a.profile in ('dielectric','ohmic','ohmic-ld'):material=Dielectric(2.5,mu=1.4,conductivity=.2 if a.profile.startswith('ohmic') else 0.)
                elif a.profile=='ohmic-unit':material=Dielectric(conductivity=.2)
                elif a.profile in ('diagonal-ld','diagonal-unit-x','tensor-pole','tensor-background'):material=TensorLorentzDrude(plasma_frequency=20.,
                    strength=[[[1.,.25,.1],[.25,2.,.2],[.1,.2,3.]]] if a.profile=='tensor-pole' else [[[1.,0.,0.],[0.,2.,0.],[0.,0.,3.]]],resonance=[0.],damping=[2.],
                    epsilon_infinity=[[2.,.3,.1],[.3,3.,.2],[.1,.2,4.]] if a.profile=='tensor-background' else [[1. if a.profile=='diagonal-unit-x' else 2.,0.,0.],[0.,3.,0.],[0.,0.,4.]])
                else:material=LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.],
                    epsilon_infinity=4. if a.profile in ('scalar-ld','subpixel','subpixel-aligned') else 1.)
                from salmonpy.subpixel import SubpixelSmoothing
                shape=Shape('ellipsoid',(9.,11.,13.),(6.,.4,.7)) if a.profile=='subpixel' else Shape('box',(6.,12.,12.),(6.,0.,0.))
                materials=(material,);shapes=(shape,)
                if a.profile=='ohmic-ld':
                    materials+=(LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.],epsilon_infinity=4.),)
                    shapes+=(Shape('box',(4.,10.,10.),(-6.,0.,0.),material=2),)
                calc=FDTD(Geometry(Grid((a.grid_count,)*3,(2.,)*3),shapes),
                    materials,.003,64,source=PlaneWave(AcosPulse(.1,20.,.15),0,-16.),observations=((-1.,1.,1.),),
                    smoothing=SubpixelSmoothing() if a.profile.startswith('subpixel') else None,
                    observation_stride=a.observation_stride)
                field=NativeEHField(context,calc,a.build,a.output/'field',points={10:(center-1,center,center)},grid=tuple(a.field_ranks if d==a.split_axis else 1 for d in range(3)),execution_mode=a.execution_mode)
                _bind_field_operations(field.lib,('sp_eh_checkpoint','sp_eh_step_indices'))
                if a.restore:
                    candidate=NativeRestoreCandidate(field,'sp_eh_checkpoint',a.restore/f'field-{world.rank}.bin')
                    candidate.stage();indices=np.empty(2,dtype=np.int32)
                    assert field._call(field.lib.sp_eh_step_indices,indices)==0 and indices[1]==steps
                    if context.is_leader:
                        for source in (a.restore/'streams').glob('*.data'):shutil.copyfile(source,field.directory/source.name)
                    candidate.publish();field.A=np.array(data['A']);field.E=np.array(data['E']);field.step_index=steps
                observer=field.energy_observer(lower,upper) if a.observe else None
                backup=_NativeBackup(field,context.solver_comm,'sp_eh_step_transaction',field.directory,
                    ('A','E','step_index','failed','closed')) if a.observe else None
                phases['setup']=time.perf_counter()-setup_start
                loop_start=time.perf_counter()
                for n in range(steps,a.stop_after):
                    sampled=(n+1)%a.observation_stride==0 or n+1==a.stop_after
                    query=observer is not None and sampled
                    start=time.perf_counter()
                    token=backup.capture() if query else None
                    phases['backup_capture']+=time.perf_counter()-start
                    current=np.full((1,3),.001 if a.driven_current and n>0 else 0.)
                    start=time.perf_counter()
                    if a.owned_current:field.advance_owned(current)
                    else:field.advance(current)
                    phases['advance']+=time.perf_counter()-start
                    start=time.perf_counter()
                    budget=observer.sample_step() if query else None
                    phases['observation']+=time.perf_counter()-start
                    if query:energy_queries+=1
                    if context.is_leader and query and not budget['available']:
                        print('energy observation unavailable for supported profile',flush=True)
                        world.Abort(2)
                    if context.is_leader and sampled:
                        history.append(dict(step=n+1,A=field.A.tolist(),E=field.E.tolist(),budget=budget))
                    start=time.perf_counter()
                    if query:backup.release(token)
                    phases['backup_release']+=time.perf_counter()-start
                phases['loop']=time.perf_counter()-loop_start
            if a.checkpoint:
                if world.rank==0:a.checkpoint.mkdir(parents=True,exist_ok=False)
                world.Barrier()
                if field is not None:
                    assert field._call(field.lib.sp_eh_checkpoint,str(a.checkpoint/f'field-{world.rank}.bin').encode(),0)==0
                    if context.is_leader:
                        (a.checkpoint/'streams').mkdir()
                        for source in field.directory.glob('*.data'):shutil.copyfile(source,a.checkpoint/'streams'/source.name)
                world.Barrier()
                if world.rank==0:
                    (a.checkpoint/'complete.json').write_text(json.dumps(dict(profile=a.profile,step=field.step_index,
                        A=field.A.tolist(),E=field.E.tolist(),history=history),indent=2)+'\n')
            if field is not None:
                checkpoint_start=time.perf_counter()
                assert field._call(field.lib.sp_eh_checkpoint,str(a.output/f'final-{world.rank}.bin').encode(),0)==0
                phases['checkpoint']=time.perf_counter()-checkpoint_start
        finally:
            if field is not None:
                start=time.perf_counter();field.close();phases['close']=time.perf_counter()-start
    from salmonpy.memory import process_memory
    rss,basis=(None,'injected resource probe failure') if a.fail_resource_probe else process_memory()
    measurements=world.gather(dict(rank=world.rank,role=role,peak_rss_bytes=rss,memory_basis=basis,
        timings_seconds=phases,energy_queries=energy_queries),root=0)
    if world.rank==0:
        inventory=None;inventory_error=None
        try:
            files=[(str(f.relative_to(a.output)),f.stat().st_size) for f in a.output.rglob('*') if f.is_file()]
            inventory=dict(file_count=len(files),bytes=sum(size for _,size in files),
                native_point_stream_bytes=sum(size for name,size in files if name.endswith('_at_point_rt.data')),
                native_point_stream_files=sum(name.endswith('_at_point_rt.data') for name,_ in files))
        except OSError as exc:inventory_error=type(exc).__name__+': '+str(exc)
        (a.output/'report.json').write_text(json.dumps(dict(profile=a.profile,history=history,
            observed=a.observe,field_ranks=a.field_ranks,grid_count=a.grid_count,
            observation_stride=a.observation_stride,execution_mode=a.execution_mode,calculation_status='completed',
            resource_measurement_status='complete' if inventory_error is None and all(x['peak_rss_bytes'] is not None for x in measurements) else 'partial',
            rank_measurements=measurements,output_inventory=inventory,output_inventory_error=inventory_error,
            timing_scope='per-rank inclusive wall seconds; advance and observation include MPI wait; communication time not isolated',
            rss_scope='per-process high-water RSS; sum is not a simultaneous job peak; children excluded',
            output_scope='files present before report.json; includes final checkpoints/static inputs; external launcher log excluded',
            energy_scope='sampled individual steps only; sparse samples are not a cumulative energy balance'),indent=2)+'\n')


if __name__=='__main__':
    try:main()
    except SystemExit:raise
    except BaseException:
        import traceback
        traceback.print_exc()
        mpi=sys.modules.get('mpi4py.MPI')
        if mpi is not None and mpi.Is_initialized() and not mpi.Is_finalized():mpi.COMM_WORLD.Abort(1)
        raise
