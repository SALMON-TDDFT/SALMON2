"""Ordinary layer/material-ID assembly with LD, SBE and different HSE inputs.

Explicit existing GS catalogs only; no build/download/SCF. 3D interface fixture,
not physical convergence. TDDFT production normally prefers one-dimensional setups.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from salmonpy.fast_hse import FastHSEMaterialProvider
from salmonpy.fast_hse_mpi import NativeFastHSEBackend
from salmonpy.material_response import MaterialBinding,ResponseSession
from salmonpy.materials import bind_material
from salmonpy.distributed import ParallelContext
from salmonpy.layers import assemble_native_layers
from salmonpy.coupled import NativeEHField,CoupledStepper
from salmonpy.geometry import Grid
from salmonpy.fdtd import PlaneWave
from salmonpy.fields import AcosPulse
from salmonpy.ld import LorentzDrude
from salmonpy.examples.mpi_mixed_multiscale import models
from salmonpy.examples._tddft_materials import load_materials
from salmonpy.examples._direct_checkpoint import checkpoint,digest


def main():
    from mpi4py import MPI
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('field_build','hse_build','materials','output'):p.add_argument(name,type=Path)
    p.add_argument('--layers',type=Path,default=Path(__file__).with_name('mixed_hse_models_layers.json'))
    p.add_argument('--separate',action='store_true',help='independent material MPI groups for comparison')
    p.add_argument('--steps',type=int,default=32)
    p.add_argument('--save-at',type=int)
    p.add_argument('--resume',type=Path)
    p.add_argument('--compare',type=Path)
    p.add_argument('--checkpoint-interval-steps',type=int,help='direct checkpoint every N accepted steps; continue')
    p.add_argument('--checkpoint-interval-seconds',type=float,help='periodic direct checkpoint in elapsed seconds; continue')
    p.add_argument('--stop-after-seconds',type=float,help='save a direct checkpoint and exit after this elapsed wall time')
    args=p.parse_args();world=MPI.COMM_WORLD
    from salmonpy.checkpointing import CheckpointSchedule
    schedule=CheckpointSchedule(interval_steps=args.checkpoint_interval_steps,
        interval_seconds=args.checkpoint_interval_seconds,stop_after_seconds=args.stop_after_seconds,comm=world)
    if args.steps<1:p.error('positive steps required')
    if args.save_at is not None and not 0<args.save_at<args.steps:p.error('save-at must lie inside the run')
    output=args.output.resolve()
    error=None
    if world.rank==0:
        try:output.mkdir(parents=True,exist_ok=False)
        except OSError as exc:error=str(exc)
    error=world.bcast(error,root=0)
    if error:raise ValueError(error)
    inputs=load_materials(args.materials)
    grid=Grid((40,40,40),(2.,2.,2.));dt=.001
    _,sbe=models()
    available={'sbe-point':bind_material(sbe,key='sbe-point',identity='synthetic-two-band-v1'),
        'drude':LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.])}
    rank_settings={'sbe-point':1}
    for mid,spec in inputs.items():
        if spec.decomposition[1:]!=(1,1,1,1):p.error('HSE requires k-only decomposition')
        available[mid]=MaterialBinding(mid,
            FastHSEMaterialProvider(spec.bundle,spec.input_file,max_steps=args.steps),
            NativeFastHSEBackend(args.hse_build,output,allow_symmetry=True,registry_model=not args.separate))
        rank_settings[mid]=spec.ranks
    layout=json.loads(args.layers.read_text())
    selected={entry['material_id'] for entry in layout['layers']}
    chosen=[mid for mid in inputs if mid in selected]
    group=None if args.separate else {'hse-models':chosen}
    assembly=assemble_native_layers(layout,available,grid=grid,dt=dt,steps=args.steps,
        source=PlaneWave(AcosPulse(.001,20.,.15),0,-16.),field_grid=(1,1,1),
        ranks_per_point=rank_settings,rank_budget=world.size,model_groups=group)
    if assembly.record()['required_ranks']!=world.size:p.error('launch exactly '+str(assembly.record()['required_ranks'])+' ranks')
    signature={'steps':args.steps,'assembly':assembly.record(),
        'field_library':json.loads(args.field_build.read_text())['sha256'],
        'hse_library':json.loads(args.hse_build.read_text())['sha256'],
        'assets':{mid:{str(path.relative_to(spec.bundle)):digest(path) for path in spec.bundle.rglob('*') if path.is_file()} for mid,spec in inputs.items()}}
    history=[]
    with ParallelContext.create(world,assembly.parallel_plan) as context:
        field=response=None
        try:
            field=NativeEHField(context,assembly.calculation,args.field_build,output/'field',
                               points=assembly.points,grid=(1,1,1)) if context.group.role=='field' else None
            response=ResponseSession(context,assembly.bindings,dt=dt)
            coupled=CoupledStepper(response,field)
            start=checkpoint(world,coupled,args.resume,signature,restore=True) if args.resume else 0
            stop=args.save_at or args.steps
            from salmonpy.native_steps import native_step_targets
            targets=native_step_targets(coupled)
            def timed_save(step):
                checkpoint(world,coupled,output/f'checkpoint-{step:08d}',signature,frames=history)
            for index in range(start,stop):
                frame=coupled.step_protected(targets)
                if frame is not None:history.append(frame)
                if schedule.after_step(index+1,timed_save):
                    stop=index+1
                    break
            if args.save_at and not schedule.stopped:checkpoint(world,coupled,output/'checkpoint',signature,frames=history)
        finally:
            if response is not None:response.close()
            if field is not None:field.close()
    if world.rank==0:
        arrays={name:np.asarray([frame[name] for frame in history]) for name in ('A_over_c','J_charge','E_at_points')}
        np.savez(output/'history.npz',**arrays)
        if args.compare:
            with np.load(args.compare/'history.npz') as reference:
                for key,values in arrays.items():np.testing.assert_array_equal(values,reference[key][start:stop])
        (output/'result.json').write_text(json.dumps({'status':'pass','steps':stop,'restored_step':start,
            'stopped_for_checkpoint':schedule.stopped,'native_steps_this_job':stop-start,'replayed_steps':0,'ranks':world.size,
            'separate':args.separate,'point_materials':assembly.record()['point_materials'] if 'point_materials' in assembly.record() else {str(k):b.key for k,b in assembly.bindings.items()}},indent=2)+'\n')
        print('LD/SBE/different-HSE layer assembly and protected continuation PASS')


if __name__=='__main__':main()
