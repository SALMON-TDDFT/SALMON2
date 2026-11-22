"""Development oscillator feedback regression; not an isolated TDDFT solver.

external: 2 ranks, illumination-only baseline. feedback: 4 ranks, two native
field contexts. Each context owns one field rank plus one bookkeeping rank.
replay: 2 ranks, accepted currents through ordinary FDTD with no trial/rollback.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def oscillator_step(p, v, E, dt, omega):
    """Prescribed-drive control using the common Fortran test kernel."""
    from salmonpy.dipole_source import _midpoint_oscillator
    q,w,_,stats=_midpoint_oscillator(p,v,E,E,dt=dt,omega=omega)
    return q,w,float(stats[5])


def centered_loop(world,field,steps,dt,omega,volume):
    """Trial/rollback/correct current; only accepted field output survives."""
    from salmonpy.native_steps import _NativeBackup
    from salmonpy.dipole_source import (_midpoint_oscillator,_molecular_work,
                                       _selected_molecular_field,_molecular_summary)
    backup=None if field is None else _NativeBackup(field,field.context.solver_comm,
        'sp_eh_step_transaction',field.directory,('A','E','step_index','failed','closed'))
    zero=np.zeros((1,3));p=np.zeros(3);v=np.zeros(3)
    if field is not None:field.advance(zero)
    gathered=world.gather(None if field is None else field.E[0].copy(),root=0)
    records=[];ps=[];js=[];statistics=[]
    if world.rank==0:
        left_total=gathered[0];left_own=gathered[2]
        records.append([left_total,left_own,_selected_molecular_field(left_total,left_own)])
        ps.append(p.copy())
    for n in range(1,steps):
        token=None if backup is None else backup.capture()
        if field is not None:field.advance(zero)
        trial=world.gather(None if field is None else field.E[0].copy(),root=0)
        if world.rank==0:
            q,w,current,_=_midpoint_oscillator(p,v,left_total,trial[0],dt=dt,omega=omega,volume=volume,
                                             own_left=left_own,own_right=trial[2])
        else:current=None
        current=world.bcast(current,root=0)
        if backup is not None:backup.rollback(token)
        if field is not None:field.advance(current[None,:])
        accepted=world.gather(None if field is None else field.E[0].copy(),root=0)
        if world.rank==0:
            checked_p,checked_v,_,_=_midpoint_oscillator(p,v,left_total,accepted[0],dt=dt,omega=omega,
                volume=volume,own_left=left_own,own_right=accepted[2])
            np.testing.assert_allclose(q,checked_p,rtol=1e-10,atol=1e-15)
            np.testing.assert_allclose(w,checked_v,rtol=1e-10,atol=1e-15)
            stat=_molecular_work(p,v,q,w,left_total,accepted[0],omega=omega,
                own_left=left_own,own_right=accepted[2])
            statistics.append(stat)
            p,v=q,w;left_total=accepted[0];left_own=accepted[2]
            records.append([left_total,left_own,_selected_molecular_field(left_total,left_own)])
            ps.append(p.copy());js.append(current.copy())
        if backup is not None:backup.release(token)
    if world.rank!=0:return None
    return np.asarray(records),np.asarray(ps),np.asarray(js),np.asarray(statistics),_molecular_summary(statistics)



def isolated_loop(world, field, args):
    """Private E-worker coupling; no quantum energy-balance claim.

    Field time = electronic time + dt after the explicit first hold. Initial
    selected E must be zero, matching the existing worker's initial E. The worker
    receives endpoint E and retains its own integration order without averaging
    it again in Python. Restoring electrons currently replays accepted history.
    """
    from salmonpy.native_steps import _NativeBackup
    from salmonpy.dipole_source import DipoleMap, _IsolatedDipoleTrial, _selected_molecular_field
    from salmonpy.isolated_tddft import load_isolated_runtime, IsolatedTDDFTInitialState
    from contextlib import ExitStack
    backup=None if field is None else _NativeBackup(field,field.context.solver_comm,
        'sp_eh_step_transaction',field.directory,('A','E','step_index','failed','closed'))
    zero=np.zeros((1,3))
    if field is not None:field.advance(zero)
    initial=world.gather(None if field is None else field.E[0].copy(),root=0)
    with ExitStack() as stack:
        bridge=None
        if world.rank==0:
            selected=_selected_molecular_field(initial[0],initial[2])
            if np.any(selected):
                raise ValueError('initial molecular E must be zero; increase source separation or delay')
            runtime=load_isolated_runtime(args.isolated_executable,args.isolated_validation)
            model=stack.enter_context(IsolatedTDDFTInitialState.from_salmon(
                args.isolated_bundle,args.isolated_input,runtime=runtime,
                max_steps=args.steps,threads=1))
            state=stack.enter_context(model.new_state(dt=args.dt,position=(-1.,1.,1.)))
            mapping=DipoleMap(('molecule',),(10,),[[1.]],[8.])
            bridge=_IsolatedDipoleTrial(state,mapping)
            ps=[state.observables.dipole.copy()];js=[];drives=[];norms=[state.electron_number]
            fields=[[initial[0],initial[2],selected]]
        n=1
        rejected=False
        while n<args.steps:
            token=None if backup is None else backup.capture()
            if field is not None:field.advance(zero)
            predicted=world.gather(None if field is None else field.E[0].copy(),root=0)
            current=bridge.propose(predicted[0],predicted[2]) if world.rank==0 else None
            current=world.bcast(current,root=0)
            if backup is not None:
                backup.rollback(token)
                # Native backups are single-use: current is move-restored.
                # Retire the predictor token and protect the accepted attempt.
                backup.release(token)
                token=backup.capture()
            if field is not None:field.advance(current)
            accepted=world.gather(None if field is None else field.E[0].copy(),root=0)
            injected=args.inject_reject_step==n and not rejected
            error=None
            if world.rank==0:
                try:
                    if injected:raise ValueError('injected endpoint rejection')
                    bridge.validate(accepted[0],accepted[2])
                except Exception as exc:error=str(exc)
            error=world.bcast(error,root=0)
            if error:
                if backup is not None:backup.rollback(token)
                if bridge is not None:bridge.cancel()
                if backup is not None:backup.release(token)
                if injected:
                    rejected=True
                    continue
                raise ValueError(error)
            if world.rank==0:
                provided=bridge._drive.copy()
                bridge.commit()
                selected=_selected_molecular_field(accepted[0],accepted[2])
                ps.append(state.observables.dipole.copy());js.append(current[0].copy());norms.append(state.electron_number)
                drives.append(provided);fields.append([accepted[0],accepted[2],selected])
            if backup is not None:backup.release(token)
            n+=1
        if world.rank==0:
            # Independent response to the accepted E trace, with no field trials.
            control=stack.enter_context(model.new_state(dt=args.dt,position=(-1.,1.,1.)))
            for E,p,Ne in zip(drives,ps[1:],norms[1:]):
                np.testing.assert_array_equal(control.step(E).dipole,p)
                np.testing.assert_equal(control.electron_number,Ne)
            np.savez(args.output/'history.npz',fields=fields,p=ps,j=js,drive=drives,electron_number=norms,
                E_time=np.arange(1,args.steps+1)*args.dt,
                p_time=np.arange(1,args.steps+1)*args.dt,
                electronic_time=np.arange(args.steps)*args.dt,
                j_time=(np.arange(len(js))+1.5)*args.dt)
            from salmonpy.convergence import _norm_report
            norm_report=_norm_report(norms,norms[0],0.,1.e-5)
            report=dict(mode='feedback',schedule='centered-trial-rollback',
                norm_diagnostic=norm_report,norm_loss_target=1.e-5,
                dt=args.dt,steps=args.steps,material_environment=args.environment,
                library_sha256=field.manifest['library_sha256'],
                isolated_executable_sha256=runtime.executable.sha256,
                TDDFT_attached=not json.loads(args.isolated_validation.read_text()).get('synthetic',False),
                response_kind='isolated-E-worker',
                injected_rejection_retried=rejected,
                synthetic_worker=json.loads(args.isolated_validation.read_text()).get('synthetic',False),
                independent_drive_replay_exact=True,
                initial_hold=args.dt,electronic_time_offset=args.dt,
                rollback='isolated accepted-history replay',full_field_energy_balance=False,
                physical_validation=False,policy='environment_scattered_only',
                isolated_input_sha256=hashlib.sha256((args.isolated_bundle/args.isolated_input).read_bytes()).hexdigest(),
                isolated_validation_sha256=hashlib.sha256(args.isolated_validation.read_bytes()).hexdigest(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
            (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
            print('PASS isolated-feedback interface; physical convergence not tested',flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('external','feedback','replay'))
    parser.add_argument('build',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--baseline',type=Path)
    parser.add_argument('--isolated-executable',type=Path)
    parser.add_argument('--isolated-validation',type=Path)
    parser.add_argument('--isolated-bundle',type=Path)
    parser.add_argument('--isolated-input',default='rt.inp')
    parser.add_argument('--inject-reject-step',type=int,help='isolated development test: reject once and retry this interval')
    parser.add_argument('--trajectory',type=Path,help='accepted centered history for a single-pass current replay')
    parser.add_argument('--environment',choices=('vacuum','ld'),default='vacuum')
    parser.add_argument('--schedule',choices=('held','centered'),default='centered',help='held is an explicit old-algorithm control, not a fallback')
    parser.add_argument('--dt',type=float,default=.003)
    parser.add_argument('--steps',type=int,default=64)
    args=parser.parse_args()
    if not np.isfinite(args.dt) or args.dt<=0 or args.steps<2:parser.error('positive dt and steps>=2 required')
    from mpi4py import MPI
    from salmonpy.coupled import NativeEHField
    from salmonpy.distributed import ParallelPlan,ParallelContext,SolverGroup
    from salmonpy.fdtd import FDTD,PlaneWave
    from salmonpy.geometry import Grid,Geometry,Shape
    from salmonpy.ld import LorentzDrude
    from salmonpy.fields import AcosPulse
    from salmonpy._sbe_handles import default_api
    from salmonpy.field_input import FieldValue
    from salmonpy.dipole_source import select_molecular_field
    world=MPI.COMM_WORLD
    feedback=args.mode=='feedback'
    isolated=args.isolated_executable is not None
    if args.inject_reject_step is not None and (not isolated or not 1<=args.inject_reject_step<args.steps):
        raise ValueError('injection needs an isolated interval in 1..steps-1')
    if isolated and (not feedback or args.schedule!='centered' or args.isolated_validation is None or args.isolated_bundle is None):
        raise ValueError('isolated coupling needs feedback, centered, validation and bundle')
    if world.size!=(4 if feedback else 2):raise ValueError('wrong rank count')
    if feedback and not isolated and args.baseline is None:raise ValueError('--baseline required')
    error=None
    if feedback and not isolated and world.rank==0:
        try:
            control=json.loads((args.baseline/'report.json').read_text())
            runtime=json.loads(args.build.read_text())
            if (control['mode']!='external' or control['dt']!=args.dt or control['steps']!=args.steps
                    or control['library_sha256']!=runtime['sha256'] or control['environment']!=args.environment):
                raise ValueError('baseline mode/dt/steps/library/environment mismatch')
        except Exception as exc:error=str(exc)
    replay_current=replay_field=None
    if args.mode=='replay' and world.rank==0:
        try:
            if args.trajectory is None:raise ValueError('--trajectory required for replay')
            original=json.loads((args.trajectory/'report.json').read_text())
            if (original['schedule']!='centered-trial-rollback' or original['dt']!=args.dt
                    or original['steps']!=args.steps or original['material_environment']!=args.environment
                    or original['library_sha256']!=json.loads(args.build.read_text())['sha256']):
                raise ValueError('replay declaration mismatch')
            with np.load(args.trajectory/'history.npz',allow_pickle=False) as saved:
                replay_current=saved['j'].copy();replay_field=saved['fields'][:,0,:].copy()
            if replay_current.shape!=(args.steps-1,3):raise ValueError('replay current extent differs')
        except Exception as exc:error=str(exc)
    error=world.bcast(error,root=0)
    if error:raise ValueError(error)
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
                    shapes=();materials=()
                    if args.environment=='ld' and pair==0:
                        shapes=(Shape('box',(6.,12.,12.),(6.,0.,0.)),)
                        materials=(LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.]),)
                    calculation=FDTD(Geometry(grid,shapes),materials,dt,steps,
                        source=PlaneWave(AcosPulse(.001 if pair==0 else 0.,20.,.15),0,-16.),
                        observations=((-1.,1.,1.),))
                    field=NativeEHField(context,calculation,args.build,args.output/('field-'+str(pair)),
                        points={10:(19,20,20)},grid=(1,1,1))
                if isolated:
                    isolated_loop(world,field,args)
                    return
                if feedback and args.schedule=='centered':
                    result=centered_loop(world,field,steps,dt,omega,volume)
                    if world.rank==0:
                        records,ps,js,stats,totals=result
                        with np.load(args.baseline/'history.npz',allow_pickle=False) as baseline:
                            expected=baseline['fields'][:,0,:]
                        if args.environment=='vacuum':
                            np.testing.assert_allclose(records[:,2,:],expected,rtol=1e-10,atol=1e-14)
                        assert np.max(abs(records[:,1,:]))>1e-20
                        feedback_max=float(np.max(abs(records[:,2,:]-expected)))
                        if args.environment=='ld':assert feedback_max>1e-20,'material feedback must be nonzero'
                        assert totals[3]<1e-12 and abs(totals[4])<1e-12
                        assert max(totals[8:11])<1e-10
                        report=dict(mode='feedback',schedule='centered-trial-rollback',dt=dt,steps=steps,
                            initial_hold=dt,policy='environment_scattered_only',omega=omega,cell_volume=volume,
                            library_sha256=field.manifest['library_sha256'],
                            oscillator_balance_max=float(totals[3]),temporal_splitting_work=float(totals[4]),
                            total_source_work=float(totals[0]),reference_source_work=float(totals[1]),
                            drive_work=float(totals[2]),held_work_negative_control=float(totals[5]),
                            accumulated_oscillator_balance=float(totals[7]),
                            relative_step_balance_max=float(totals[8]),
                            relative_accumulated_balance=float(totals[9]),
                            relative_temporal_work_error=float(totals[10]),
                            full_field_energy_balance=False,physical_validation=False,
                            material_environment=args.environment,TDDFT_attached=False,
                            scattered_feedback_max=feedback_max,
                            molecular_kernel_sha256=default_api().sha256,
                            baseline_sha256={name:hashlib.sha256((args.baseline/name).read_bytes()).hexdigest() for name in ('report.json','history.npz')},
                            native_input_sha256={str(k):hashlib.sha256((args.output/('field-'+str(k))/'inputfile').read_bytes()).hexdigest() for k in range(2)},
                            script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
                        np.savez(args.output/'history.npz',fields=records,p=ps,j=js,statistics=stats,
                            E_time=np.arange(1,steps+1)*dt,p_time=np.arange(1,steps+1)*dt,
                            j_time=(np.arange(len(js))+1.5)*dt)
                        (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
                        print('PASS centered-feedback',flush=True)
                    return
                p=np.zeros(3);v=np.zeros(3);current=np.zeros((1,3))
                previous=None;records=[];balance=[];work=[];ps=[p.copy()];js=[]
                for n in range(1,steps+1):
                    if args.mode=='replay' and world.rank==0:
                        if n>=2:current=replay_current[n-2][None,:]
                        js.append(current[0].copy())
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
                    report=dict(mode=args.mode,environment=args.environment,dt=dt,steps=steps,initial_hold=dt,
                                full_field_energy_balance=False,physical_validation=False,
                                library_sha256=field.manifest['library_sha256'])
                    if args.mode=='replay':
                        np.testing.assert_array_equal(records[:,0,:],replay_field)
                        report['accepted_field_replay_exact']=True
                        report['trajectory_sha256']={name:hashlib.sha256((args.trajectory/name).read_bytes()).hexdigest() for name in ('report.json','history.npz')}
                    report.update(omega=omega,cell_volume=volume,policy='environment_scattered_only' if feedback else 'accepted_current_replay' if args.mode=='replay' else 'illumination_only',
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
                             j_time=(np.arange(len(js))+(0.5 if args.mode=='replay' else 1.5))*dt,work=work,balance=balance)
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
