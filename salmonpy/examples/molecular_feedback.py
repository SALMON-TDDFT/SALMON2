"""Single molecular dipole coupled to caller-owned FDTD groups.

mpiexec -n 6 python -m salmonpy.examples.molecular_feedback FIELD_BUILD \
    E_WORKER VALIDATION_JSON EXISTING_GS NEW_OUTPUT --field-grid 2 1 1

This explicit calculation reuses a prepared GS; it never builds or runs SCF.
All values are atomic units. The geometry is the small C2H2 wiring example,
not a converged molecular spectrum. The own-vacuum policy removes vacuum
radiation reaction. Electronic rollback uses a direct live-worker savepoint, without replay.
"""
import argparse
from contextlib import ExitStack
import json
from math import prod
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('field_build','worker','validation','bundle','output'):
        parser.add_argument(name,type=Path)
    parser.add_argument('--input',default='input.inp')
    parser.add_argument('--field-grid',nargs=3,type=int,default=(1,1,1))
    parser.add_argument('--steps',type=int,default=64)
    parser.add_argument('--dt',type=float,default=.003)
    parser.add_argument('--execution-mode',choices=('checked','performance'),default='checked',help='field, molecule and outer feedback mode; trial/restore and explicit validation remain')
    parser.add_argument('--ld',action='store_true',help='add the verified synthetic Drude box')
    parser.add_argument('--checkpoint',type=Path,help='new direct bundle directory')
    parser.add_argument('--resume',type=Path,help='direct bundle to load before updating')
    parser.add_argument('--stop-after',type=int,help='accepted FDTD step at which to save and exit')
    parser.add_argument('--checkpoint-interval-steps',type=int,help='direct checkpoint every N accepted steps; continue')
    parser.add_argument('--checkpoint-interval-seconds',type=float,help='periodic direct checkpoint in elapsed seconds; continue')
    parser.add_argument('--stop-after-seconds',type=float,help='save a direct checkpoint and exit after this elapsed wall time')
    args=parser.parse_args()
    from mpi4py import MPI
    from salmonpy.distributed import SolverGroup,ParallelPlan,ParallelContext
    from salmonpy.coupled import NativeEHField
    from salmonpy.geometry import Grid,Geometry,Shape
    from salmonpy.fdtd import FDTD,PlaneWave
    from salmonpy.fields import AcosPulse
    from salmonpy.ld import LorentzDrude
    from salmonpy.isolated_tddft import load_isolated_runtime,IsolatedTDDFTInitialState
    from salmonpy.dipole_source import DipoleMap,MolecularFeedback
    world=MPI.COMM_WORLD
    from salmonpy.checkpointing import CheckpointSchedule
    schedule=CheckpointSchedule(interval_steps=args.checkpoint_interval_steps,
        interval_seconds=args.checkpoint_interval_seconds,stop_after_seconds=args.stop_after_seconds,comm=world)
    nfield=prod(args.field_grid);group_size=nfield+1
    if any(n<1 for n in args.field_grid) or world.size!=2*group_size:
        raise ValueError('launch 2*(NX*NY*NZ+1) ranks with positive field-grid dimensions')
    pair=world.rank//group_size
    local_world=world.Split(pair,world.rank)
    plan=ParallelPlan((SolverGroup('field',tuple(range(nfield)),role='field'),
                      SolverGroup('bookkeeping',(nfield,),point_ids=(10,))))
    try:
        with ExitStack() as states,ParallelContext.create(local_world,plan) as context:
            state=None;error=None
            if world.rank==0:
                try:
                    args.output.mkdir(parents=True,exist_ok=False)
                    runtime=load_isolated_runtime(args.worker,args.validation)
                    model=states.enter_context(IsolatedTDDFTInitialState.from_salmon(
                        args.bundle,args.input,runtime=runtime,max_steps=args.steps,threads=1))
                    state=states.enter_context(model.new_state(dt=args.dt,position=(-1.,1.,1.),execution_mode=args.execution_mode))
                except Exception as exc:error=str(exc)
            error=world.bcast(error,root=0)
            if error:raise RuntimeError(error)
            field=None
            try:
                error=None
                try:
                    if context.group.role=='field':
                        shapes=();materials=()
                        if args.ld and pair==0:
                            shapes=(Shape('box',(6.,12.,12.),(6.,0.,0.)),)
                            materials=(LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.]),)
                        calculation=FDTD(Geometry(Grid((40,40,40),(2.,)*3),shapes),materials,args.dt,args.steps,
                            source=PlaneWave(AcosPulse(.1 if pair==0 else 0.,20.,.15),0,-16.),
                            observations=((-1.,1.,1.),))
                        field=NativeEHField(context,calculation,args.field_build,args.output/('field-'+str(pair)),
                            points={10:(19,20,20)},grid=tuple(args.field_grid),execution_mode=args.execution_mode)
                except Exception as exc:error=str(exc)
                errors=world.allgather(error)
                if any(errors):raise RuntimeError(str(errors))
                observer=None if field is None else field.energy_observer((13,13,13),(27,27,27))
                coupling=MolecularFeedback(world,field,state=state,
                    mapping=DipoleMap(('molecule',),(10,),[[1.]],[8.]),
                    total_root=0,reference_root=group_size,policy='environment_scattered_only',observer=observer,execution_mode=args.execution_mode)
                restored_steps=None
                if args.resume:
                    saved=coupling.restore_checkpoint(args.resume)
                    records=saved['history'] if world.rank==0 else None
                    error='sample resume requires saved frame history' if world.rank==0 and records is None else None
                    error=world.bcast(error,root=0)
                    if error:raise ValueError(error)
                    restored_steps=saved['verified_steps'] if world.rank==0 else None
                else:
                    first=coupling.start()
                    records=[first.record()] if world.rank==0 else None
                stop=args.steps if args.stop_after is None else args.stop_after
                if not coupling.step_index<=stop<=args.steps:raise ValueError('stop-after must be between current FDTD step and capacity')
                restore_statistics=state.worker_statistics() if world.rank==0 and args.resume else None
                def timed_save(step):
                    coupling.save_checkpoint(args.output/f'checkpoint-{step:08d}',history=records)
                for _ in range(stop-coupling.step_index):
                    frame=coupling.step()
                    if world.rank==0:records.append(frame.record())
                    if schedule.after_step(coupling.step_index,timed_save):break
                if args.checkpoint:coupling.save_checkpoint(args.checkpoint,history=records)
                statistics=state.worker_statistics() if world.rank==0 else None
            finally:
                if field is not None:field.close()
        # Finish collective cleanup before root-only file I/O.
        if world.rank==0:
            with (args.output/'molecular-feedback.json').open('x') as out:
                json.dump(dict(stopped_for_checkpoint=schedule.stopped,policy='environment_scattered_only',field_grid=args.field_grid,
                    units='atomic',execution_mode=args.execution_mode,
                    molecular_execution_mode=state.execution_mode,initial_hold=args.dt,frames=records,worker_statistics=statistics,restore_statistics=restore_statistics,restored_steps=restored_steps),out,indent=2)
            print(args.output/'molecular-feedback.json')
    finally:local_world.Free()


if __name__=='__main__':main()
