"""Development spatial electric/orbital-magnetic molecular/FDTD feedback and direct restart.

Run one field group for external controls, two for feedback. Each group has
FIELD_RANKS solver ranks and one bookkeeping rank. This is a finite-scope
validation driver, not a general production failure coordinator.
"""
import argparse
from contextlib import ExitStack,contextmanager
import json
import os
from pathlib import Path
import shutil
import numpy as np


@contextmanager
def trial_guard(world):
    try:yield
    except Exception:
        import traceback
        traceback.print_exc();world.Abort(1);raise


def main(argv=None,*,mode=None):
    parser=argparse.ArgumentParser(description=("One launch: spatial LD/FDTD E/H drives a molecule; no molecular current return. Fixed C2H2 demonstration geometry." if mode=='one-way' else __doc__))
    if mode is None:parser.add_argument('mode',choices=('external','feedback','one-way'))
    else:parser.set_defaults(mode=mode)
    for name in ('field_build','worker','validation','bundle','output'):parser.add_argument(name,type=Path)
    parser.add_argument('--steps',type=int,default=8)
    parser.add_argument('--field-ranks',type=int,default=1)
    parser.add_argument('--ld',action='store_true')
    parser.add_argument('--orbital-magnetic',action='store_true',help=argparse.SUPPRESS if mode=='one-way' else 'Yee E/H input and orbital magnetic-current return')
    parser.add_argument('--baseline',type=Path,help=argparse.SUPPRESS if mode=='one-way' else 'reference output')
    parser.add_argument('--inject-reject-step',type=int,help=argparse.SUPPRESS if mode=='one-way' else 'development rejection step')
    parser.add_argument('--save-step',type=int)
    parser.add_argument('--resume',type=Path)
    parser.add_argument('--kernels',type=Path,help=argparse.SUPPRESS)
    parser.add_argument('--stop-after',type=int,help='save all states and exit at this molecular step')
    args=parser.parse_args(argv)
    if args.mode=='one-way':args.orbital_magnetic=True
    if args.stop_after is not None and not 1<=args.stop_after<=args.steps:
        parser.error('stop-after must be between one and steps')
    if args.mode!='feedback' and args.inject_reject_step is not None:
        parser.error('rejection injection belongs to the feedback validation route')
    if args.kernels:
        import salmonpy._material_handles as handles
        handles._library_path=lambda:args.kernels.resolve()
    from mpi4py import MPI
    from salmonpy.distributed import ParallelPlan,ParallelContext,SolverGroup
    from salmonpy.coupled import NativeEHField
    from salmonpy.geometry import Grid,Geometry,Shape
    from salmonpy.fdtd import FDTD,PlaneWave
    from salmonpy.fields import AcosPulse
    from salmonpy.ld import LorentzDrude
    from salmonpy.near_field import GridFieldSeries
    from salmonpy.isolated_tddft import load_isolated_runtime,IsolatedTDDFTInitialState
    from salmonpy.native_steps import _NativeBackup
    from salmonpy.api import _bind_field_operations
    from salmonpy.convergence import _comparison
    from salmonpy._native_checkpoint import NativeRestoreCandidate
    from salmonpy.dipole_source import _MolecularStreams
    world=MPI.COMM_WORLD;feedback=args.mode=='feedback';quantum=args.mode!='external';size=args.field_ranks+1
    if args.field_ranks not in (1,2) or world.size!=size*(2 if feedback else 1):
        raise ValueError('launch (field-ranks+1) ranks per group')
    if args.steps<1:raise ValueError('positive step count required')
    if args.save_step is not None and not 1<=args.save_step<=args.steps:
        raise ValueError('save-step must be between one and steps')
    if args.save_step is not None and args.stop_after is not None and args.save_step!=args.stop_after:
        raise ValueError('choose save-step or stop-after; this sample writes one checkpoint directory')
    pair=world.rank//size;local_world=world.Split(pair,world.rank)
    dt=.003;start=np.array([8,8,8] if args.orbital_magnetic else [10,12,9])
    counts=(24,24,24) if args.orbital_magnetic else (17,17,21);position=(-2.54,1.62,1.94)
    indices=np.array([[x,y,z] for z in range(start[2],start[2]+counts[2])
                      for y in range(start[1],start[1]+counts[1])
                      for x in range(start[0],start[0]+counts[0])],dtype=np.int32)
    shape=(counts[2],counts[1],counts[0],3)
    if world.rank==0:args.output.mkdir(parents=True,exist_ok=False)
    world.Barrier()
    plan=ParallelPlan((SolverGroup('field',tuple(range(args.field_ranks)),role='field'),
                      SolverGroup('labels',(args.field_ranks,),point_ids=(10,))))
    with ExitStack() as stack,ParallelContext.create(local_world,plan) as context,trial_guard(world):
        field=None;backup=None;bounds=None
        if context.group.role=='field':
            shapes=();materials=()
            if args.ld and pair==0:
                shapes=(Shape('box',(4.,12.,12.),(17.,0.,0.)),)
                materials=(LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.]),)
            try:
                calculation=FDTD(Geometry(Grid((40,40,40),(2.,)*3),shapes),materials,dt,args.steps+1,
                    source=PlaneWave(AcosPulse(.1 if pair==0 else 0.,20.,.15),0,-20.),observations=((-1.,1.,1.),))
                field=NativeEHField(context,calculation,args.field_build,args.output/f'field-{pair}',
                                    points={10:(19,20,20)},grid=(args.field_ranks,1,1))
            except BaseException:
                import traceback
                traceback.print_exc();world.Abort(1);raise
            if feedback:
                backup=_NativeBackup(field,context.solver_comm,'sp_eh_step_transaction',field.directory,('A','E','step_index','failed','closed'))
        geometry=(field.component_field_geometry() if args.orbital_magnetic else field.electric_edge_geometry()) if field is not None else None
        geometry=world.bcast(geometry if world.rank==0 else None,root=0)
        origin=geometry['origin']+start*geometry['spacing']
        def window(left,right,n):
            if args.orbital_magnetic:
                return GridFieldSeries(origin,geometry['spacing'],[0.,dt],np.array([left[0],right[0]]),
                    np.array([left[1],right[1]]),molecule_position=position,time_origin=n*dt)
            return GridFieldSeries(np.stack((origin,origin)),geometry['spacing'],[0.,dt],
                np.array([left,right]).reshape((2,*shape)),np.zeros((2,*shape)),
                molecule_position=position,time_origin=n*dt)
        def selected(total,own,n):return window(total[0],total[1],n).difference(window(own[0],own[1],n))
        state=None;model=None
        if world.rank==0 and quantum:
            source=args.output/'molecule-input';source.mkdir()
            for name in ('C_rps.dat','H_rps.dat','input.inp'):shutil.copy2(args.bundle/name,source/name)
            shutil.copytree(args.bundle/'data_for_restart',source/'data_for_restart')
            empty=np.zeros((args.steps+2,*shape))
            GridFieldSeries(origin if args.orbital_magnetic else np.stack((origin,origin)),geometry['spacing'],np.arange(args.steps+2)*dt,empty,empty,
                            molecule_position=position).write(source/'salmonpy-near-eh.dat')
            runtime=load_isolated_runtime(args.worker,args.validation)
            model=stack.enter_context(IsolatedTDDFTInitialState.from_salmon(source,'input.inp',runtime=runtime,
                                                                          max_steps=args.steps,near_field=True))
            state=stack.enter_context(model.new_state(dt=dt,position=position))
        restored_steps=None
        traces=observations=None
        if args.resume:
            meta=json.loads((args.resume/'metadata.json').read_text())
            expected=dict(mode=args.mode,steps=args.steps,field_ranks=args.field_ranks,ld=args.ld,orbital_magnetic=args.orbital_magnetic,dt=dt,
                          field_sha256=json.loads(args.field_build.read_text())['sha256'])
            if any(meta[k]!=v for k,v in expected.items()):raise ValueError('checkpoint configuration differs')
            if field is not None:
                candidate=NativeRestoreCandidate(field,'sp_eh_checkpoint',args.resume/f'field-{world.rank}.bin')
                candidate.stage()
                streams=None
                if context.is_leader:
                    streams=_MolecularStreams(field,args.resume/f'streams-{world.rank}');streams.publish()
                candidate.publish()
                if streams:streams.close()
                with np.load(args.resume/f'python-{world.rank}.npz',allow_pickle=False) as saved:
                    field.A=saved['A'].copy();field.E=saved['E'].copy()
                field.step_index=meta['field_step']
                _bind_field_operations(field.lib,('sp_eh_step_indices',))
                actual_steps=np.empty(2,dtype=np.int32)
                assert field._call(field.lib.sp_eh_step_indices,actual_steps)==0
                assert tuple(actual_steps)==(meta['field_step'],-1)
            if world.rank==0:
                if state is not None:
                    state.restore_direct_checkpoint(args.resume/'molecule')
                    assert state.worker_statistics()['updates']==0
                    assert state.worker_statistics()['step']==meta['molecular_step']
                with np.load(args.resume/'history.npz',allow_pickle=False) as saved:
                    traces={k:list(saved[k]) for k in ('total','own','selected')};observations=list(saved['observations'])
                restored_steps=(meta['field_step'],meta['molecular_step'])
                if args.mode=='one-way':shutil.copy2(args.resume/'field-step-0.npz',args.output/'field-step-0.npz')
            nstart=meta['molecular_step']+1
        else:
            if args.mode=='one-way':
                zero=None
                if field is not None:
                    value=field.component_fields(indices).reshape((2,*shape))
                    if context.is_leader:zero=value
                zeros=world.gather(zero,root=0)
                if world.rank==0:np.savez(args.output/'field-step-0.npz',EH=zeros[0])
            if field is not None:field.advance(np.zeros((1,3)))
            nstart=1
        if field is not None:
            _bind_field_operations(field.lib,('sp_eh_plane_read','sp_eh_checkpoint','sp_eh_step_indices'))
            if feedback:
                info=np.empty(21,dtype=np.int32);dummy=np.empty(0)
                assert field._call(field.lib.sp_eh_plane_read,1,20,0,info.ctypes.data,dummy.ctypes.data)==0
                bounds=context.solver_comm.allgather((info[11:14]-info[5:8],info[14:17]-info[5:8]))
        def collect():
            data=None
            if field is not None:
                values=field.component_fields(indices).reshape((2,*shape)) if args.orbital_magnetic else field.electric_edges(indices).reshape(shape)
                if context.is_leader:data=values
            return world.gather(data,root=0)
        initial=collect()
        if world.rank==0:
            if traces is None:
                own=initial[size] if feedback else np.zeros_like(initial[0])
                if feedback:
                    first_drive=selected((initial[0],initial[0]),(own,own),1)
                    first=np.stack((first_drive.E[0],first_drive.H[0])) if args.orbital_magnetic else first_drive.E[0]
                else:first=initial[0]
                traces=dict(total=[initial[0]],own=[own],selected=[first])
                observations=[] if state is None else [state._engine.current.copy()]
            previous_total=traces['total'][-1];previous_own=traces['own'][-1]
            density=state.charge_density() if feedback else None
        rejected=False;drive_errors=[];magnetic_peak=0.
        for n in range(nstart,args.steps+1):
            while True:
                token=None if backup is None else backup.capture()
                if field is not None:field.advance(np.zeros((1,3)))
                trial=collect()
                if not feedback:
                    if backup:backup.release(token)
                    if world.rank==0:
                        if quantum:
                            state.step_grid_field(window(previous_total,trial[0],n),step_index=n)
                            observations.append(state._engine.current.copy())
                        previous_total=trial[0]
                        traces['total'].append(trial[0]);traces['own'].append(np.zeros_like(trial[0]));traces['selected'].append(trial[0])
                    break
                payloads=None;quantum_token=None
                if world.rank==0:
                    drive=selected((previous_total,trial[0]),(previous_own,trial[size]),n)
                    previous_magnetic=state.orbital_magnetic_response() if args.orbital_magnetic else None
                    quantum_token=state.savepoint();state.step_grid_field(drive,step_index=n)
                    next_density=state.charge_density()
                    current=drive.charge_current(next_density.positions,density.charge_density,next_density.charge_density,
                                                volume_elements=next_density.volume_elements,dt=dt)
                    if args.orbital_magnetic:
                        magnetic=state.orbital_magnetic_response()
                        magnetic_peak=max(magnetic_peak,float(np.max(np.abs(magnetic.current_density(previous=previous_magnetic)))))
                        current=magnetic.total_current(current,previous=previous_magnetic)
                    payloads=[current if rank in (0,size) else None for rank in range(world.size)]
                group_current=world.scatter(payloads,root=0)
                if field is not None:
                    backup.rollback(token);backup.release(token);token=backup.capture()
                    packets=None
                    if context.is_leader:
                        values=group_current.reshape(-1,3);packets=[]
                        for lower,upper in bounds:
                            keep=np.all((indices>=lower)&(indices<=upper),axis=1)
                            packets.append((np.ascontiguousarray(indices[keep]),np.ascontiguousarray(values[keep])))
                    owned_indices,owned_current=context.solver_comm.scatter(packets,root=0)
                    field.advance_edges(owned_indices,owned_current)
                accepted=collect();decision=None
                if world.rank==0:
                    actual=selected((previous_total,accepted[0]),(previous_own,accepted[size]),n)
                    error=_comparison(np.stack((drive.E,drive.H)) if args.orbital_magnetic else drive.E,
                                      np.stack((actual.E,actual.H)) if args.orbital_magnetic else actual.E,1e-13,1e-10)
                    if not error['within_tolerance']:raise AssertionError('selected drive changed after radiation: '+str(error))
                    decision=args.inject_reject_step==n and not rejected
                    drive_errors.append(error['absolute_error'])
                reject=world.bcast(decision,root=0)
                if reject:
                    if backup:backup.rollback(token);backup.release(token)
                    if world.rank==0:
                        state.restore_savepoint(quantum_token);state.release_savepoint(quantum_token)
                        np.testing.assert_array_equal(state.charge_density().charge_density,density.charge_density)
                    rejected=True
                    continue
                if backup:backup.release(token)
                if world.rank==0:
                    state.release_savepoint(quantum_token);density=next_density
                    previous_total=accepted[0];previous_own=accepted[size]
                    traces['total'].append(previous_total);traces['own'].append(previous_own);traces['selected'].append(np.stack((actual.E[1],actual.H[1])) if args.orbital_magnetic else actual.E[1])
                    observations.append(state._engine.current.copy())
                    assert state.worker_statistics()['step']==n
                break
            if n==args.save_step or n==args.stop_after:
                checkpoint=args.output/'checkpoint'
                if world.rank==0:checkpoint.mkdir()
                world.Barrier()
                if field is not None:
                    assert field._call(field.lib.sp_eh_checkpoint,os.fsencode(checkpoint/f'field-{world.rank}.bin'),0)==0
                    np.savez(checkpoint/f'python-{world.rank}.npz',A=field.A,E=field.E)
                    if context.is_leader:
                        folder=checkpoint/f'streams-{world.rank}';folder.mkdir()
                        for file in field.directory.glob('*.data'):shutil.copy2(file,folder/file.name)
                if world.rank==0:
                    if state is not None:state.save_direct_checkpoint(checkpoint/'molecule')
                    np.savez(checkpoint/'history.npz',**traces,observations=observations)
                    if args.mode=='one-way':shutil.copy2(args.output/'field-step-0.npz',checkpoint/'field-step-0.npz')
                    meta=dict(mode=args.mode,steps=args.steps,field_ranks=args.field_ranks,ld=args.ld,orbital_magnetic=args.orbital_magnetic,dt=dt,field_step=n+1,molecular_step=n,
                              field_sha256=json.loads(args.field_build.read_text())['sha256'])
                    (checkpoint/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
                world.Barrier()
            if n==args.stop_after:break
        if world.rank==0:
            environment_difference=None;trajectory_control=False
            if args.baseline:
                with np.load(args.baseline/'history.npz',allow_pickle=False) as saved:
                    if quantum and 'observations' in saved and saved['observations'].size:
                        trajectory_control=True
                        for name in traces:np.testing.assert_array_equal(traces[name],saved[name])
                        np.testing.assert_array_equal(observations,saved['observations'])
                    elif feedback and not args.ld:
                        np.testing.assert_allclose(traces['selected'],saved['selected'],rtol=1e-10,atol=1e-13)
                    elif feedback:
                        environment_difference=_comparison(traces['selected'],saved['selected'],0.,0.)
                        assert environment_difference['absolute_error']>0.
            np.savez(args.output/'history.npz',**traces,observations=observations,
                     field_step=np.arange(1,len(traces['total'])+1),molecular_step=np.arange(len(observations)),
                     electronic_field_half_step=2*np.arange(len(observations))+1,
                     electronic_field_time=(np.arange(len(observations))+.5)*dt)
            if state is not None:
                folder=args.output/'molecule-streams';folder.mkdir()
                for file in (Path(state._engine._work.name)/'run').glob('*.data'):shutil.copy2(file,folder/file.name)
            summary=dict(electronic_clock_offset='1/2',used_grid_frame='n; n+1 is prefetch only for middlepoint',mode=args.mode,ld=args.ld,molecular_steps=state.worker_statistics()['step'] if state is not None else None,field_steps=field.step_index,
                worker_statistics=state.worker_statistics() if state is not None else None,
                feedback_enabled=feedback,self_vacuum_reference=feedback,
                field_ranks=args.field_ranks,position_bohr=position,patch_counts=counts,
                policy='environment_scattered_only' if feedback else 'prescribed_EH',electric_only_feedback=feedback and not args.orbital_magnetic,magnetic_return=feedback and args.orbital_magnetic,
                magnetic_current_peak=magnetic_peak,
                maximum_drive_change=max(drive_errors,default=0.),rejected_and_retried=rejected,
                restored_steps=restored_steps,restore_extra_evolutions=0 if args.resume else None,
                baseline_matched=bool(args.baseline) and (trajectory_control or not args.ld),
                environment_difference=environment_difference,physical_convergence=False)
            (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
        if field is not None:field.close()
    local_world.Free()


def entry(argv=None,*,mode=None):
    """Terminate all ranks on an unhandled example error; never hang peers."""
    try:main(argv,mode=mode)
    except Exception:
        import traceback
        from mpi4py import MPI
        traceback.print_exc();MPI.COMM_WORLD.Abort(1);raise


if __name__=='__main__':entry()
