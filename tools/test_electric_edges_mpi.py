"""Explicit raw E-component current tests, with optional real TDDFT density history.

Launch FIELD_RANKS+1 ranks; the last rank only labels the macro point.
The density-history case is prescribed molecular radiation with an explicit first
field hold, not a self-consistent spatial-field molecular feedback driver.
"""
import argparse
import json
import os
from pathlib import Path
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('macro','edges','density','response'))
    parser.add_argument('build',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--response-history',type=Path)
    parser.add_argument('--response-component',choices=('combined','electric','magnetic'),default='combined')
    parser.add_argument('--baseline',type=Path);parser.add_argument('--density-history',type=Path)
    parser.add_argument('--save-step',type=int)
    parser.add_argument('--resume',type=Path)
    parser.add_argument('--kernels',type=Path,help='explicit host calculation library for density response')
    args=parser.parse_args()
    if args.kernels is not None:
        import salmonpy._material_handles as handles
        handles._library_path=lambda:args.kernels.resolve()
    from mpi4py import MPI
    from salmonpy.distributed import ParallelPlan,ParallelContext,SolverGroup
    from salmonpy.coupled import NativeEHField
    from salmonpy.fdtd import FDTD,PlaneWave
    from salmonpy.geometry import Grid,Geometry
    from salmonpy.fields import AcosPulse
    from salmonpy.api import _bind_field_operations
    from salmonpy.native_steps import _NativeBackup
    from salmonpy.near_field import GridFieldSeries
    world=MPI.COMM_WORLD;nfield=world.size-1
    if nfield not in (1,2):raise ValueError('test supports one/two field ranks plus one bookkeeping rank')
    if args.mode=='density' and args.density_history is None:raise ValueError('density history required')
    if args.mode=='response' and args.response_history is None:raise ValueError('response history required')
    plan=ParallelPlan((SolverGroup('field',tuple(range(nfield)),role='field'),SolverGroup('labels',(nfield,),point_ids=(10,))))
    if world.rank==0:args.output.mkdir(parents=True,exist_ok=False)
    world.Barrier()
    summary=None
    with ParallelContext.create(world,plan) as context:
        if context.group.role=='field':
            comm=context.solver_comm
            density=None
            if args.mode=='density':
                with np.load(args.density_history,allow_pickle=False) as saved:
                    density={k:saved[k].copy() for k in saved.files}
                steps=len(density['charge_density']);dt=float(density['dt'])
            elif args.mode=='response':
                with np.load(args.response_history,allow_pickle=False) as saved:
                    response={k:saved[k].copy() for k in saved.files}
                steps=len(response[args.response_component])+1;dt=float(response['dt'])
            else:steps=16;dt=.003
            calculation=FDTD(Geometry(Grid((40,40,40),(2.,)*3)),(),dt,steps,
                source=PlaneWave(AcosPulse(0.,20.,.15),0,-16.),observations=((-1.,1.,1.),))
            with NativeEHField(context,calculation,args.build,args.output/'field',points={10:(19,20,20)},grid=(nfield,1,1)) as field:
                geometry=field.electric_edge_geometry()
                expected=np.full((3,3),-39.)+np.eye(3)
                np.testing.assert_array_equal(geometry['origin'],expected)
                query=np.array([[x,y,z] for z in (19,20,21) for y in (19,20,21) for x in (18,19,20)],dtype=np.int32)
                np.testing.assert_array_equal(field.electric_edges(query),0.)
                restored_fields=restored_raw=None
                _bind_field_operations(field.lib,('sp_eh_checkpoint','sp_eh_step_indices'))
                start=2;restored_step=None
                if args.resume:
                    saved=json.loads((args.resume/'metadata.json').read_text())
                    expected=dict(mode=args.mode,response_component=args.response_component,field_ranks=nfield,steps=steps,dt=dt,
                                  library_sha256=json.loads(args.build.read_text())['sha256'])
                    if any(saved[k]!=v for k,v in expected.items()):raise ValueError('checkpoint configuration differs')
                    code=field._call(field.lib.sp_eh_checkpoint,os.fsencode(args.resume/f'field-{comm.rank}.bin'),1)
                    assert code==0
                    with np.load(args.resume/f'python-{comm.rank}.npz',allow_pickle=False) as data:
                        field.A=data['A'].copy();field.E=data['E'].copy()
                        restored_fields=list(data['fields']);restored_raw=list(data['raw'])
                    field.step_index=saved['step'];start=field.step_index+1;restored_step=field.step_index
                    actual=np.empty(2,dtype=np.int32)
                    assert field._call(field.lib.sp_eh_step_indices,actual)==0
                    assert tuple(actual)==(restored_step,-1)

                # Fresh runs hold once; restore goes directly to the saved step.
                if not args.resume:field.advance(np.zeros((1,3)))
                _bind_field_operations(field.lib,('sp_eh_plane_read','sp_eh_step_edges'))
                info=np.empty(21,dtype=np.int32);dummy=np.empty(0)
                status=field._call(field.lib.sp_eh_plane_read,1,20,0,info.ctypes.data,dummy.ctypes.data)
                assert status==0
                lower=info[11:14]-info[5:8];upper=info[14:17]-info[5:8]
                bounds=comm.allgather((lower.tolist(),upper.tolist()))
                mapping=None
                if args.mode=='density' and comm.rank==0:
                    shape=(2,40,40,40,3);zeros=np.zeros(shape)
                    mapping=GridFieldSeries(np.stack((geometry['origin'],geometry['origin'])),geometry['spacing'],[0.,dt],zeros,zeros,
                                            molecule_position=(-2.54,1.62,1.94))
                fields=[field.E.copy()] if restored_fields is None else restored_fields
                raw=[field.electric_edges(query)] if restored_raw is None else restored_raw
                counts=[]
                _bind_field_operations(field.lib,('sp_eh_checkpoint','sp_eh_step_indices'))
                backup=_NativeBackup(field,comm,'sp_eh_step_transaction',field.directory,('A','E','step_index','failed','closed'))
                invalid_preserved=False;rollback_exact=False
                for n in range(start,steps+1):
                    vector=np.array([.0001,-.0002,.0003]) if n<=12 else np.zeros(3)
                    packets=None
                    if comm.rank==0:
                        if args.mode=='density':
                            current=mapping.charge_current(density['positions'],density['charge_density'][n-2],density['charge_density'][n-1],
                                volume_elements=density['volume_elements'],dt=dt).reshape(-1,3)
                            full_indices=np.array([[x,y,z] for z in range(40) for y in range(40) for x in range(40)],dtype=np.int32)
                            active=np.any(current!=0.,axis=1);full_indices=full_indices[active];current=current[active]
                        elif args.mode=='response':
                            grid_current=response[args.response_component][n-2]
                            nz,ny,nx,_=grid_current.shape;begin=response['start']
                            full_indices=np.array([[x+begin[0],y+begin[1],z+begin[2]] for z in range(nz)
                                                   for y in range(ny) for x in range(nx)],dtype=np.int32)
                            current=grid_current.reshape(-1,3);active=np.any(current!=0.,axis=1)
                            full_indices=full_indices[active];current=current[active]
                        else:
                            full_indices=[];current=[]
                            for d in range(3):
                                for side in (0,1):
                                    index=np.array([19,20,20]);index[d]-=side
                                    item=np.zeros(3);item[d]=.5*vector[d]
                                    full_indices.append(index);current.append(item)
                            full_indices=np.array(full_indices,dtype=np.int32);current=np.array(current)
                        packets=[]
                        for lo,hi in bounds:
                            selected=np.all((full_indices>=lo)&(full_indices<=hi),axis=1)
                            packets.append((np.ascontiguousarray(full_indices[selected]),np.ascontiguousarray(current[selected])))
                    indices,current=comm.scatter(packets,root=0)
                    counts.append(comm.allgather(len(indices)))
                    if n==2 and args.mode not in ('macro','response'):
                        before=field.electric_edges(query);bad=current.copy()
                        if comm.rank==0:bad[0,0]=np.nan
                        output=np.empty_like(field.E)
                        code=field._call(field.lib.sp_eh_step_edges,len(field.point_ids),len(indices),indices,bad,field.A,output)
                        assert code==10
                        np.testing.assert_array_equal(field.electric_edges(query),before)
                        assert field.step_index==1 and not field.failed
                        invalid_preserved=True
                        # Finite duplicated inputs overflow inside the solver,
                        # exercising the existing deep checked recovery gate.
                        from salmonpy.errors import NumericalError
                        overflow_indices=np.array([[19,20,20],[19,20,20]],dtype=np.int32) if comm.rank==0 else np.empty((0,3),dtype=np.int32)
                        overflow_current=np.full((len(overflow_indices),3),np.finfo(float).max)
                        try:field.advance_edges(overflow_indices,overflow_current)
                        except NumericalError:pass
                        else:raise AssertionError('deep overflow was not rejected')
                        assert field.step_index==1 and not field.failed
                        np.testing.assert_array_equal(field.electric_edges(query),before)
                    token=backup.capture() if n==3 else None
                    def advance():
                        if args.mode=='macro':field.advance(vector[None,:])
                        else:field.advance_edges(indices,current)
                    advance()
                    if token is not None:
                        accepted_E=field.E.copy();accepted_A=field.A.copy();accepted_raw=field.electric_edges(query)
                        backup.rollback(token);advance()
                        np.testing.assert_array_equal(field.E,accepted_E);np.testing.assert_array_equal(field.A,accepted_A)
                        np.testing.assert_array_equal(field.electric_edges(query),accepted_raw)
                        backup.release(token);rollback_exact=True
                    fields.append(field.E.copy());raw.append(field.electric_edges(query))
                    if n==args.save_step:
                        checkpoint=args.output/'checkpoint'
                        if comm.rank==0:checkpoint.mkdir()
                        comm.Barrier()
                        assert field._call(field.lib.sp_eh_checkpoint,os.fsencode(checkpoint/f'field-{comm.rank}.bin'),0)==0
                        np.savez(checkpoint/f'python-{comm.rank}.npz',A=field.A,E=field.E,fields=fields,raw=raw)
                        if comm.rank==0:
                            metadata=dict(mode=args.mode,response_component=args.response_component,field_ranks=nfield,steps=steps,dt=dt,step=n,
                                          library_sha256=json.loads(args.build.read_text())['sha256'])
                            (checkpoint/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')

                assert np.any(np.asarray(raw)!=0.)
                if comm.rank==0:
                    fields=np.array(fields);raw=np.array(raw)
                    if args.baseline:
                        with np.load(args.baseline/'history.npz',allow_pickle=False) as expected:
                            np.testing.assert_array_equal(fields,expected['E']);np.testing.assert_array_equal(raw,expected['raw_E'])
                    np.savez(args.output/'history.npz',E=fields,raw_E=raw)
                    summary=dict(mode=args.mode,field_ranks=nfield,steps=steps,initial_hold=dt,
                        baseline_exact=bool(args.baseline),invalid_current_preserves_state=invalid_preserved,deep_overflow_preserves_state=invalid_preserved,
                        direct_live_rollback_retry_exact=rollback_exact,owned_entry_counts=counts,
                        radiation_nonzero=True,component_geometry_exact=True,physical_convergence=False,restored_step=restored_step,
                        full_spatial_feedback=False,restore_setup_evolution_calls=0 if args.resume else None)
                    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    summary=world.bcast(summary,root=0)
    if world.rank==0:print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
