"""One MPI launch: caller-owned FDTD group -> MPI-parallel molecular TDDFT group.

Fixed C2H2/LD connection example. No nested launcher, child process, current return
or own-vacuum field. Direct restart requires the same MPI layout and build.
General deep-error retry is not provided by this dedicated-rank example.
"""
import argparse
from pathlib import Path
import shutil
import json
import os
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('build_record','bundle','output'):parser.add_argument(name,type=Path)
    parser.add_argument('--field-ranks',type=int,default=2)
    parser.add_argument('--molecule-ranks',type=int,default=2)
    parser.add_argument('--steps',type=int,default=4)
    parser.add_argument('--ld',action='store_true')
    parser.add_argument('--stop-after',type=int)
    parser.add_argument('--resume',type=Path)
    args=parser.parse_args()
    if args.stop_after is not None and not 1<=args.stop_after<=args.steps:parser.error('invalid stop-after')
    from mpi4py import MPI
    from salmonpy.distributed import ParallelPlan,ParallelContext,SolverGroup
    from salmonpy.coupled import NativeEHField
    from salmonpy.geometry import Grid,Geometry,Shape
    from salmonpy.fdtd import FDTD,PlaneWave
    from salmonpy.fields import AcosPulse
    from salmonpy.ld import LorentzDrude
    from salmonpy.near_field import GridFieldSeries
    from salmonpy.isolated_tddft import run_molecular_mpi,MolecularMPIController
    from salmonpy.outputs import save_arrays
    from salmonpy.tddft import _digest
    world=MPI.COMM_WORLD
    if args.field_ranks not in (1,2) or args.molecule_ranks<1 or args.steps<1:
        parser.error('one/two field ranks, positive molecular ranks and steps required')
    if world.size!=args.field_ranks+args.molecule_ranks:
        parser.error('MPI ranks must equal field-ranks + molecule-ranks')
    transport=world.Dup();work=args.output.resolve()/'molecule-work'
    plan=ParallelPlan((SolverGroup('field',tuple(range(args.field_ranks)),role='field'),
        SolverGroup('molecule',tuple(range(args.field_ranks,world.size)),point_ids=(10,))))
    field=None;controller=None
    try:
        if world.rank==0:args.output.mkdir(parents=True,exist_ok=False)
        world.Barrier()
        with ParallelContext.create(world,plan) as context:
            geometry=None
            if context.group.role=='field':
                shapes=(Shape('box',(4.,12.,12.),(17.,0.,0.)),) if args.ld else ()
                materials=(LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.]),) if args.ld else ()
                calculation=FDTD(Geometry(Grid((40,40,40),(2.,)*3),shapes),materials,.003,args.steps+1,
                    source=PlaneWave(AcosPulse(.1,20.,.15),0,-20.),observations=((-1.,1.,1.),))
                field=NativeEHField(context,calculation,args.build_record,args.output/'field',
                    points={10:(19,20,20)},grid=(args.field_ranks,1,1))
                geometry=field.component_field_geometry()
            geometry=world.bcast(geometry if world.rank==0 else None,root=0)
            start=np.array([8,8,8]);shape=(24,24,24,3);dt=.003
            origin=geometry['origin']+start*geometry['spacing'];position=(-2.54,1.62,1.94)
            indices=np.array([[x,y,z] for z in range(8,32) for y in range(8,32) for x in range(8,32)],dtype=np.int32)
            if world.rank==args.field_ranks:
                work.mkdir()
                for name in ('input.inp','C_rps.dat','H_rps.dat'):shutil.copy2(args.bundle/name,work/name)
                shutil.copytree(args.bundle/'data_for_restart',work/'data_for_restart')
                empty=np.zeros((args.steps+2,*shape))
                GridFieldSeries(origin,geometry['spacing'],np.arange(args.steps+2)*dt,empty,empty,
                    molecule_position=position).write(work/'salmonpy-near-eh.dat')
            world.Barrier()
            if context.group.role!='field':
                run_molecular_mpi(args.build_record,context.solver_comm,transport,controller_rank=0,
                    work_directory=work,decomposition=(1,args.molecule_ranks,1,1,1))
            else:
                if world.rank==0:
                    controller=MolecularMPIController(transport,solver_root=args.field_ranks,
                        work_directory=work,dt=dt,max_steps=args.steps,decomposition=(1,args.molecule_ranks,1,1,1))
                    frames=[controller.observables];field_frames=[]
                # Retain actual field-step0, then the established one-step hold.
                zero=field.component_fields(indices).reshape((2,*shape))
                if not args.resume:field.advance(np.zeros((1,3)))
                previous=field.component_fields(indices).reshape((2,*shape))
                if world.rank==0:field_frames=[zero,previous]
                start_step=0
                if args.resume:
                    from salmonpy._native_checkpoint import NativeRestoreCandidate
                    from salmonpy.dipole_source import _MolecularStreams
                    from salmonpy.api import _bind_field_operations
                    meta=json.loads((args.resume/'complete.json').read_text())
                    record=json.loads(args.build_record.read_text())
                    expected=dict(steps=args.steps,field_ranks=args.field_ranks,molecule_ranks=args.molecule_ranks,
                                  ld=args.ld,library_sha256=record['sha256'])
                    if (meta.get('schema')!='salmonpy.one-way-molecular-mpi/1' or
                            any(meta.get(k)!=v for k,v in expected.items())):
                        raise ValueError('joint restart configuration differs')
                    error=None
                    try:
                        files=meta['files']
                        for i,(name,sha) in enumerate(sorted(files.items())):
                            path=Path(name)
                            if path.is_absolute() or '..' in path.parts or '\\' in name:
                                raise ValueError('invalid joint checkpoint filename')
                            item=args.resume/path
                            if i%context.solver_comm.size==context.solver_comm.rank:
                                if item.is_symlink() or not item.is_file() or _digest(item)!=sha:
                                    raise ValueError('joint checkpoint checksum differs: '+name)
                    except Exception as exc:error=str(exc)
                    errors=context.solver_comm.allgather(error)
                    if any(errors):raise ValueError(next(error for error in errors if error))
                    candidate=NativeRestoreCandidate(field,'sp_eh_checkpoint',args.resume/f'field-{world.rank}.bin')
                    candidate.stage()
                    streams=None
                    if context.is_leader:
                        streams=_MolecularStreams(field,args.resume/'field-streams');streams.publish()
                    candidate.publish()
                    if streams:streams.close()
                    with np.load(args.resume/f'field-python-{world.rank}.npz') as data:
                        field.A=data['A'].copy();field.E=data['E'].copy()
                    start_step=meta['molecular_step'];field.step_index=start_step+1
                    _bind_field_operations(field.lib,('sp_eh_step_indices',))
                    actual=np.empty(2,dtype=np.int32)
                    if field._call(field.lib.sp_eh_step_indices,actual) or tuple(actual)!=(start_step+1,-1):
                        raise ValueError('restored field step differs')
                    previous=field.component_fields(indices).reshape((2,*shape))
                    if world.rank==0:
                        controller.restore_direct_checkpoint(args.resume/'molecule')
                        restore_statistics=controller.worker_statistics()
                        if restore_statistics!={'step':start_step,'updates':0}:raise ValueError('restore evolved molecular state')
                        with np.load(args.resume/'history.npz') as saved:
                            field_frames=list(saved['fields']);saved_values=saved['observations']
                        from salmonpy.isolated_tddft import DipoleObservables
                        frames=[DipoleObservables._from_worker(i*dt,value,checked=True) for i,value in enumerate(saved_values)]
                        for name in ('dipole','electronic_change','electron_number'):
                            np.testing.assert_array_equal(getattr(controller.observables,name),getattr(frames[-1],name))
                stop=args.steps if args.stop_after is None else args.stop_after
                if stop<start_step:raise ValueError('stop-after precedes restored step')
                for n in range(start_step+1,stop+1):
                    field.advance(np.zeros((1,3)))
                    current=field.component_fields(indices).reshape((2,*shape))
                    if world.rank==0:
                        window=GridFieldSeries(origin,geometry['spacing'],[0.,dt],
                            np.array([previous[0],current[0]]),np.array([previous[1],current[1]]),
                            molecule_position=position,time_origin=n*dt)
                        frames.append(controller.step_grid_field(window,step_index=n))
                        field_frames.append(current)
                    previous=current
                if world.rank==0:
                    stats=controller.worker_statistics()
                if args.stop_after is not None:
                    from salmonpy.api import _bind_field_operations
                    checkpoint=args.output/'checkpoint'
                    if context.is_leader:checkpoint.mkdir()
                    context.solver_comm.Barrier()
                    _bind_field_operations(field.lib,('sp_eh_checkpoint',))
                    if field._call(field.lib.sp_eh_checkpoint,os.fsencode(checkpoint/f'field-{world.rank}.bin'),0):
                        raise RuntimeError('field direct save failed')
                    np.savez(checkpoint/f'field-python-{world.rank}.npz',A=field.A,E=field.E)
                    context.solver_comm.Barrier()
                    if world.rank==0:
                        folder=checkpoint/'field-streams';folder.mkdir()
                        for path in field.directory.glob('*.data'):shutil.copy2(path,folder/path.name)
                        controller.save_direct_checkpoint(checkpoint/'molecule')
                        values=np.array([[*f.electronic_change,*f.dipole,f.electron_number,
                            0. if f.internal_energy is None else f.internal_energy,
                            float(f.internal_energy is not None)] for f in frames])
                        np.savez(checkpoint/'history.npz',fields=field_frames,observations=values)
                        files={}
                        for item in sorted(checkpoint.rglob('*')):
                            if not item.is_file():continue
                            if item.is_relative_to(checkpoint/'molecule') and item.name!='manifest.json':continue
                            files[item.relative_to(checkpoint).as_posix()]=_digest(item)
                        meta=dict(schema='salmonpy.one-way-molecular-mpi/1',files=files,steps=args.steps,field_ranks=args.field_ranks,molecule_ranks=args.molecule_ranks,
                            ld=args.ld,library_sha256=json.loads(args.build_record.read_text())['sha256'],molecular_step=stop)
                        (checkpoint/'complete.json').write_text(json.dumps(meta,indent=2)+'\n')
                if world.rank==0:controller.close()
                field.close();field=None
        # Both groups have returned and detached; the caller's MPI remains live.
        alive=world.allreduce(1)
        if alive!=world.size:raise RuntimeError('caller MPI lost a rank')
        if world.rank==0:
            save_arrays(args.output/'response',{
                'molecular_step':np.arange(stop+1),
                'time':np.array([f.time for f in frames]),
                'dipole':np.array([f.dipole for f in frames]),
                'electronic_change':np.array([f.electronic_change for f in frames]),
                'electron_number':np.array([f.electron_number for f in frames]),
                'E':np.array([f[0] for f in field_frames]),'H':np.array([f[1] for f in field_frames]),
                'field_step':np.arange(len(field_frames))},metadata={
                'coupling':'one-way-host-MPI','material_id':'molecule','point_id':10,'field_ranks':args.field_ranks,
                'molecule_ranks':args.molecule_ranks,'decomposition':[1,args.molecule_ranks,1,1,1],
                'feedback':False,'self_vacuum_reference':False,'caller_MPI_alive':True,
                'worker_statistics':stats,'restore_statistics':restore_statistics if args.resume else None,'dt':dt,'physical_convergence':False,
                'units':{'time':'atomic_time','dipole':'e*bohr','electronic_change':'e*bohr','E':'E:au','H':'H:au'}})
            print('One-way FDTD / molecular MPI response saved:',args.output/'response')
    except Exception:
        import traceback
        traceback.print_exc();world.Abort(1);raise
    finally:
        transport.Free()


if __name__=='__main__':main()
