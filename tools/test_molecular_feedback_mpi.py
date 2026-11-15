"""Development oscillator feedback regression; not an isolated TDDFT solver.

external: 2 ranks, illumination-only baseline. feedback: 4 ranks, two native
field contexts. Each context owns one field rank plus one bookkeeping rank.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def oscillator_step(p, v, E, dt, omega):
    """Implicit midpoint with constant prescribed drive on one interval."""
    new_p = ((1-(omega*dt/2)**2)*p+dt*v+.5*dt*dt*E)/(1+(omega*dt/2)**2)
    new_v = 2*(new_p-p)/dt-v
    energy_old = .5*np.sum(v*v+omega*omega*p*p)
    energy_new = .5*np.sum(new_v*new_v+omega*omega*new_p*new_p)
    work = float(np.dot(E,new_p-p))
    return new_p,new_v,energy_new-energy_old-work


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('external','feedback'))
    parser.add_argument('build',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--baseline',type=Path)
    parser.add_argument('--dt',type=float,default=.003)
    parser.add_argument('--steps',type=int,default=64)
    args=parser.parse_args()
    if not np.isfinite(args.dt) or args.dt<=0 or args.steps<2:parser.error('positive dt and steps>=2 required')
    from mpi4py import MPI
    from salmonpy.coupled import NativeEHField
    from salmonpy.distributed import ParallelPlan,ParallelContext,SolverGroup
    from salmonpy.fdtd import FDTD,PlaneWave
    from salmonpy.geometry import Grid,Geometry
    from salmonpy.fields import AcosPulse
    from salmonpy.field_input import FieldValue
    from salmonpy.dipole_source import select_molecular_field
    world=MPI.COMM_WORLD
    feedback=args.mode=='feedback'
    if world.size!=(4 if feedback else 2):raise ValueError('wrong rank count')
    if feedback and args.baseline is None:raise ValueError('--baseline required')
    pair=world.rank//2
    sub=world.Split(pair,world.rank)
    dt=args.dt;steps=args.steps;volume=8.;omega=2.
    plan=ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('labels',(1,),point_ids=(10,))))
    try:
        if world.rank==0:args.output.mkdir(parents=True,exist_ok=False)
        world.Barrier()
        field=None
        with ParallelContext.create(sub,plan) as context:
            try:
                if context.group.role=='field':
                    grid=Grid((40,40,40),(2.,)*3)
                    calculation=FDTD(Geometry(grid),(),dt,steps,
                        source=PlaneWave(AcosPulse(.001 if pair==0 else 0.,20.,.15),0,-16.),
                        observations=((-1.,1.,1.),))
                    field=NativeEHField(context,calculation,args.build,args.output/('field-'+str(pair)),
                        points={10:(19,20,20)},grid=(1,1,1))
                p=np.zeros(3);v=np.zeros(3);current=np.zeros((1,3))
                previous=None;records=[];balance=[];work=[];ps=[p.copy()];js=[]
                for n in range(1,steps+1):
                    if field is not None:field.advance(current)
                    gathered=world.gather(None if field is None else field.E[0].copy(),root=0)
                    if world.rank==0:
                        total=gathered[0];own=gathered[2] if feedback else np.zeros(3)
                        if feedback:
                            provenance=dict(footprint_id='native-cell-10',sample_map_id='native-E-average',
                                            own_source_history_id='accepted-p')
                            drive=select_molecular_field(FieldValue(total,n*dt,[-1,1,1],'E:au'),
                                molecule_id='oscillator',policy='environment_scattered_only',provenance=provenance,
                                reference=FieldValue(own,n*dt,[-1,1,1],'E:au'),
                                reference_provenance=dict(provenance,source_ids=['oscillator'],medium='vacuum',incident=False)).field.value
                        else:drive=total
                        records.append([total.copy(),own.copy(),drive.copy()])
                        if previous is not None:
                            old_total,old_own,old_drive,dp=previous
                            work.append([float(np.dot(.5*(old_total+total),dp)),
                                         float(np.dot(.5*(old_own+own),dp)),float(np.dot(old_drive,dp))])
                        if feedback and n<steps:
                            new_p,new_v,residual=oscillator_step(p,v,drive,dt,omega)
                            dp=new_p-p;balance.append(residual);previous=(total,own,drive,dp)
                            current=(dp/(volume*dt))[None,:];p,v=new_p,new_v
                            ps.append(p.copy());js.append(current[0].copy())
                    current=world.bcast(current if world.rank==0 else None,root=0)
                if world.rank==0:
                    records=np.asarray(records)
                    report=dict(mode=args.mode,dt=dt,steps=steps,initial_hold=dt,
                                full_field_energy_balance=False,physical_validation=False,
                                library_sha256=field.manifest['library_sha256'])
                    report.update(omega=omega,cell_volume=volume,policy='environment_scattered_only' if feedback else 'illumination_only',
                        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                        native_input_sha256={str(k):hashlib.sha256((args.output/('field-'+str(k))/'inputfile').read_bytes()).hexdigest() for k in range(2 if feedback else 1)})
                    if feedback:
                        report['baseline_sha256']={name:hashlib.sha256((args.baseline/name).read_bytes()).hexdigest() for name in ('report.json','history.npz')}
                        baseline_report=json.loads((args.baseline/'report.json').read_text())
                        assert baseline_report['mode']=='external' and baseline_report['dt']==dt and baseline_report['steps']==steps
                        assert baseline_report['library_sha256']==field.manifest['library_sha256']
                        with np.load(args.baseline/'history.npz',allow_pickle=False) as baseline:
                            expected=baseline['fields'][:,0,:]
                        np.testing.assert_allclose(records[:,2,:],expected,rtol=1e-10,atol=1e-14)
                        assert np.max(np.abs(records[:,1,:]))>1e-20,'trivial zero self field'
                        assert max(abs(x) for x in balance)<1e-12,'oscillator work balance'
                        work=np.asarray(work)
                        report.update(drive_error_max=float(np.max(np.abs(records[:,2,:]-expected))),
                            self_field_max=float(np.max(np.abs(records[:,1,:]))),
                            oscillator_balance_max=float(np.max(np.abs(balance))),
                            total_source_work=float(work[:,0].sum()),reference_source_work=float(work[:,1].sum()),
                            drive_work=float(work[:,2].sum()),
                            temporal_splitting_work=float((work[:,0]-work[:,1]-work[:,2]).sum()))
                    np.savez(args.output/'history.npz',fields=records,E_time=np.arange(1,steps+1)*dt,
                             p=ps,j=js,p_time=np.arange(1,len(ps)+1)*dt,
                             j_time=(np.arange(len(js))+1.5)*dt,work=work,balance=balance)
                    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
                    print('PASS',args.mode,flush=True)
            except BaseException:
                # Abort before collective field/context cleanup on a subset of ranks.
                import traceback
                traceback.print_exc();world.Abort(1);raise
            finally:
                if field is not None:field.close()
    except BaseException:
        world.Abort(1);raise
    finally:sub.Free()


if __name__=='__main__':main()
