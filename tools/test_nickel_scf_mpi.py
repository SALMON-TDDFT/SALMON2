"""Ni SOC SCF snapshots/interleaving versus independent controls, all arrays in RAM.

mpiexec -n 3*SIZE python tools/test_nickel_scf_mpi.py BUILD_JSON NEW_OUTPUT --potential UPF --size SIZE
No GS orbital/density/restart files are collected. Only input/log/scalar evidence is written.
"""
import argparse
from dataclasses import replace,fields
import json
from pathlib import Path
from types import SimpleNamespace
import sys
import numpy as np
from mpi4py import MPI

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from salmonpy.dft import DFT,Pseudopotential,SCFSettings,_restart_output_text
from salmonpy.native_gs import NativeGSBackend
from salmonpy.pseudopotentials import catalog_reference
from test_nickel_soc import nickel_model


class InMemoryDFT(DFT):
    def input_text(self,scf=SCFSettings()):
        return _restart_output_text(super().input_text(scf),False)


def frame(state):
    return {name:state.array(name).values for name in
            ('density','occupations','eigenvalues','orbitals','scf_history')}


def compare(left,right):
    errors={}
    for name in left:
        if left[name].shape!=right[name].shape or not np.isfinite(right[name]).all():
            raise AssertionError('invalid '+name)
        errors[name]=float(np.max(np.abs(left[name]-right[name])))
        if not np.array_equal(left[name],right[name]):
            raise AssertionError('SCF differs after restore/interleaving: '+name+' '+str(errors[name]))
    return errors


def alignment_study(session,model,args,comm,work):
    """Long SCF comparison; electronic arrays stay in RAM, no RT or checkpoint."""
    import csv,time,warnings,os
    from dataclasses import asdict
    from salmonpy.errors import SCFReviewRequired
    if args.stage_limit<1 or args.batch<1:raise ValueError('positive stage limit and batch required')
    from salmonpy.ld import _positive
    _positive(args.preliminary_residual,'preliminary residual')
    _positive(args.preliminary_moment_change,'preliminary moment change')
    settings=SCFSettings(max_iterations=1+4*args.stage_limit,tolerance=1e-9,
        mixing_method='simple_dm',mixing_rate=.3,adaptive_mixing=True,
        initial_wavefunction='random',random_seed=0)
    schedule=([0.,0.,.002],[0.,0.,.0005],[0.,0.,.0001],[0.,0.,0.])
    configuration={'schema':'salmonpy.nickel-spin-alignment/1','stage_limit':args.stage_limit,
        'batch':args.batch,'charge_tolerance':1e-9,'moment_tolerance_muB':1e-5,
        'preliminary_residual':args.preliminary_residual,'preliminary_moment_change_muB':args.preliminary_moment_change,
        'SCF':asdict(settings),'field_schedule_au':schedule,'MPI_ranks':comm.size,
        'grid':[16,16,16],'kmesh':[2,2,2],'electrons':18,'states':36,'temperature_k':300,
        'runtime_sha256':session.backend.sha256,'OMP_NUM_THREADS':os.environ.get('OMP_NUM_THREADS'),
        'wave_density_restart_saved':False,
        'scope':'fixed coarse input SCF stability; no grid/k/temperature or physical-material convergence certification'}
    control_stages=None;control_rows=[]
    if args.control_record is not None:
        import hashlib
        ref=args.control_record.resolve();previous_config=json.loads((ref.parent/'configuration.json').read_text())
        for name in ('SCF','runtime_sha256','MPI_ranks','grid','kmesh','electrons','states','temperature_k','OMP_NUM_THREADS','stage_limit'):
            if previous_config[name]!=configuration[name]:raise ValueError('control configuration differs: '+name)
        control_stages=json.loads(ref.read_text())
        if len(control_stages)!=1 or any(control_stages[0]['field_au']):raise ValueError('zero-field control required')
        configuration['reused_zero_control']={'path':str(ref),'sha256':hashlib.sha256(ref.read_bytes()).hexdigest()}
        if comm.rank==0:
            with (ref.parent/'history.csv').open() as f:
                for row in csv.DictReader(f):
                    if row['case']!='zero-control':continue
                    control_rows.append({k:v if k=='case' else None if v=='' else float(v) for k,v in row.items()})
    if comm.rank==0:(work/'configuration.json').write_text(json.dumps(configuration,indent=2)+'\n')
    started=time.perf_counter();rows=[];results={}
    columns=['case','stage','iteration','Bx_au','By_au','Bz_au','residual','energy_hartree',
        'zeeman_energy_hartree','energy_without_auxiliary_hartree','mx_muB','my_muB','mz_muB',
        'moment_norm_muB','moment_change_muB','sample_interval','alignment_cosine_z','mixing_rate','mixing_changes','elapsed_seconds']
    stream=(work/'history.csv').open('w',newline='') if comm.rank==0 else None
    writer=csv.DictWriter(stream,fieldnames=columns) if stream else None
    if writer:writer.writeheader();stream.flush()
    def sample(state,label,stage,previous):
        info=state.info;spin=state.scf_spin_field;policy=state.scf_mixing
        moment=np.asarray(spin['spin_magnetic_moment_muB']);length=float(np.linalg.norm(moment))
        difference=None if previous is None else float(np.max(np.abs(moment-np.asarray(previous['moment']))))
        interval=0 if previous is None else info.iterations-previous['iteration']
        row=dict(zip(columns,[label,stage,info.iterations,*spin['field_au'],
            info.residual if info.residual<1e100 else None,info.energy_hartree,
            spin['zeeman_energy_hartree'],spin['energy_without_auxiliary_hartree'],*moment.tolist(),
            length,difference,interval,None if length<1e-12 else float(moment[2]/length),
            policy['rate'],policy['changes'],time.perf_counter()-started]))
        if comm.rank==0:
            writer.writerow(row);stream.flush();rows.append(row)
            (work/'progress.json').write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
            print('ALIGNMENT_SCF',json.dumps(row,allow_nan=False),flush=True)
        return info,{'moment':moment.tolist(),'iteration':info.iterations},difference
    def run_case(state,label,fields):
        stages=[]
        for stage,field in enumerate(fields):
            state.set_scf_spin_field(field)
            start=state.iterations;stop=start+args.stage_limit;previous=None
            preliminary=label=='field-ramp' and any(field)
            residual_target=args.preliminary_residual if preliminary else settings.tolerance
            moment_target=args.preliminary_moment_change if preliminary else 1e-5
            info,previous,difference=sample(state,label,stage,previous)
            def ready():
                charge_ok=info.residual<=residual_target if preliminary else info.converged
                return charge_ok and difference is not None and difference<=moment_target
            while not ready() and not info.converged and info.iterations<stop:
                # Resolve changes of moment near charge convergence, not just at25-step intervals.
                count=1 if info.residual<1e-6 else args.batch
                try:state.continue_scf(steps=min(count,stop-info.iterations))
                except SCFReviewRequired:pass
                info,previous,difference=sample(state,label,stage,previous)
            stable=difference is not None and difference<=moment_target
            status={'stage':stage,'field_au':field,'start_iteration':start,'end_iteration':info.iterations,
                'charge_converged':info.converged,'moment_stable_over_last_sample':stable,
                'stage_target_met':ready(),'preliminary_stage':preliminary,'residual_target':residual_target,'moment_change_target_muB':moment_target,
                'moment_change_muB':difference,'sample_interval':rows[-1]['sample_interval'] if comm.rank==0 else None,
                'info':asdict(info),'spin':state.scf_spin_field,'mixing':state.scf_mixing}
            stages.append(status)
            if comm.rank==0:(work/(label+'-stages.json')).write_text(json.dumps(stages,indent=2,allow_nan=False)+'\n')
            if not ready():
                if comm.rank==0:print('SCF_REVIEW_REQUIRED',label,stage,'No automatic field reduction after an unresolved stage.',flush=True)
                break
        state.write_scf_report(work/(label+'-residuals'))
        return {'stages':stages,'completed_schedule':len(stages)==len(fields) and all(x['stage_target_met'] for x in stages)}
    try:
        with session.run(model,scf=settings,steps=1,material_id='ni-zero-control') as zero:
            with session.run(model,scf=settings,steps=1,material_id='ni-field-ramp') as field:
                if args.control_record is not None:
                    import hashlib
                    matched=False
                    for path in args.control_record.resolve().parent.glob('gs-*/gs.inp'):
                        pp=path.parent/'pseudo1.upf'
                        if path.read_text()==zero._input_text and pp.is_file() and hashlib.sha256(pp.read_bytes()).hexdigest()==model.pseudopotentials[0].sha256:
                            matched=True;break
                    if not matched:raise ValueError('saved control input/pseudopotential differs')
                left=frame(zero);right=frame(field);error=None;initial_difference={};initial_exact=True
                try:
                    for name in left:
                        if not np.isfinite(left[name]).all() or not np.isfinite(right[name]).all():raise AssertionError('nonfinite initial '+name)
                        np.testing.assert_allclose(left[name],right[name],rtol=2e-12,atol=2e-14,equal_nan=False)
                        initial_difference[name]=float(np.max(np.abs(left[name]-right[name])))
                        initial_exact=initial_exact and np.array_equal(left[name],right[name])
                except Exception as exc:error=str(exc)
                errors=comm.allgather(error)
                if any(errors):raise AssertionError('initial state comparison: '+str(errors))
                initial_records=comm.gather({'bitwise_equal':initial_exact,'differences':initial_difference},root=0)
                del left,right

                results['zero_control']=run_case(zero,'zero-control',[[0.,0.,0.]]) if control_stages is None else {'stages':control_stages,'completed_schedule':bool(control_stages[0]['charge_converged'] and control_stages[0]['moment_stable_over_last_sample']),'reused_from':str(args.control_record.resolve())}
                results['field_ramp']=run_case(field,'field-ramp',schedule)
    finally:
        if stream:stream.close()
    if comm.rank==0:
        unexpected=[str(p) for p in work.rglob('*') if p.is_file() and p.suffix in ('.bin','.npz')]
        if unexpected:raise AssertionError('unexpected electronic data export: '+str(unexpected))
        result={'status':'completed' if all(x['completed_schedule'] for x in results.values()) else 'SCF review required',
            'cases':results,'initial_state_comparison':{'rtol':2e-12,'atol':2e-14,'ranks':initial_records},'GS_datasets_collected':False,
            'TDDFT_started':False,'physical_convergence_verified':False,'elapsed_seconds':time.perf_counter()-started}
        (work/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        try:
            import matplotlib
            matplotlib.use('Agg')
            from matplotlib import pyplot as plt
            fig,axes=plt.subplots(3,1,figsize=(8,9),sharex=True)
            for label in ('zero-control','field-ramp'):
                all_selected=[r for r in rows+control_rows if r['case']==label]
                selected=[r for r in all_selected if r['residual'] is not None]
                axes[0].semilogy([r['iteration'] for r in selected],[max(r['residual'],1e-30) for r in selected],label=label)
                for component in ('mx_muB','my_muB','mz_muB'):
                    axes[1].plot([r['iteration'] for r in selected],[r[component] for r in selected],label=label+' '+component)
                axes[2].step([r['iteration'] for r in all_selected],[r['Bz_au'] for r in all_selected],where='post',label=label)
            axes[0].axhline(1e-9,color='gray',linestyle=':');axes[0].set_ylabel('SCF residual')
            axes[1].set_ylabel('Spin moment (mu_B)');axes[2].set_ylabel('Auxiliary Bz (au)');axes[2].set_xlabel('Accepted SCF iteration (field applies to next step)')
            for ax in axes:ax.legend(fontsize=7);ax.grid(alpha=.2)
            fig.tight_layout();fig.savefig(work/'scf-alignment.png',dpi=150);plt.close(fig)
            result['plot_status']='complete'
        except Exception as exc:
            result['plot_status']='unavailable';result['plot_error']=str(exc)
            warnings.warn('SCF data retained; plot unavailable: '+str(exc))
        (work/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print('Ni SOC long SCF study:',result['status'],work/'result.json',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--potential',type=Path,required=True)
    parser.add_argument('--size',type=int,default=1)
    parser.add_argument('--alignment-study',action='store_true',help='long memory-only zero-field/auxiliary-field comparison; launch size ranks')
    parser.add_argument('--stage-limit',type=int,default=1000)
    parser.add_argument('--preliminary-residual',type=float,default=1e-5)
    parser.add_argument('--preliminary-moment-change',type=float,default=1e-3)
    parser.add_argument('--control-record',type=Path,help='explicit completed zero-control-stages.json from the same configuration')
    parser.add_argument('--batch',type=int,default=25)
    parser.add_argument('--adaptive-mixing',action='store_true')
    parser.add_argument('--extend-limit',type=int,help='validate live and terminal-cap SCF budget extension')
    args=parser.parse_args();world=MPI.COMM_WORLD
    expected=args.size if args.alignment_study else 3*args.size
    if args.size<1 or world.size!=expected:raise ValueError('launch exactly '+str(expected)+' ranks')
    group=0 if args.alignment_study else world.rank//args.size
    comm=world.Dup() if args.alignment_study else world.Split(group,world.rank)
    work=args.output if args.alignment_study else args.output/str(group)
    if comm.rank==0:work.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    reference=catalog_reference('pseudodojo-pbe-0.4-fr-standard-ni')
    potential=Pseudopotential.from_file(args.potential,atomic_number=28,local_channel=2,reference=reference)
    base=nickel_model(potential)
    a=InMemoryDFT(**{f.name:getattr(base,f.name) for f in fields(DFT)})
    b=replace(a,cell=a.cell*1.02)
    sa=SCFSettings(max_iterations=12,tolerance=1e-12,mixing_method='simple_dm',mixing_rate=.2,adaptive_mixing=args.adaptive_mixing,
                   initial_wavefunction='random',random_seed=0)
    sb=replace(sa,mixing_rate=.35,adaptive_mixing=False,random_seed=7)
    context=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    backend=NativeGSBackend(args.build,work,decomposition=(args.size,1,1,1,1))
    own={};report={}
    if args.alignment_study:
        with backend.open(context) as session:
            alignment_study(session,a,args,comm,work)
        comm.Free()
        return
    with backend.open(context) as session:
        if group==0:
            with session.run(a,scf=sa,steps=2,material_id='ni-a') as first:
                start=frame(first)
                initial_policy=first.scf_mixing if args.adaptive_mixing else None
                if args.extend_limit is not None:
                    before=first.scf_status
                    first.extend_scf_limit(args.extend_limit)
                    report['extension_boundary_errors']=compare(start,frame(first))
                    if first.scf_status['executed_iterations']!=before['executed_iterations']:
                        raise AssertionError('budget extension executed SCF')
                snapshot=first.snapshot()
                try:
                    with session.run(b,scf=sb,steps=2,material_id='ni-b') as second:
                        if second.scf_status['max_iterations']!=sb.max_iterations:raise AssertionError('sibling cap changed')
                        first.continue_scf(steps=2)
                        own['a']=frame(first)
                        before=first.scf_status['executed_iterations']
                        first.restore(snapshot)
                        if args.adaptive_mixing:
                            if first.scf_mixing!=initial_policy:raise AssertionError('adaptive state not restored')
                            if second.scf_mixing['enabled']:raise AssertionError('sibling adaptive mode changed')
                        report['restored_boundary_errors']=compare(start,frame(first))
                        report['restore_scf_iterations']=first.scf_status['executed_iterations']-before
                        if report['restore_scf_iterations']!=0:raise AssertionError('restore recalculated SCF')
                        first.continue_scf(steps=2)
                        report['continuation_errors']=compare(own['a'],frame(first))
                        second.continue_scf(steps=2)
                        own['b']=frame(second)
                        if first.scf_status['initializations']!=1 or second.scf_status['initializations']!=1:
                            raise AssertionError('state was reinitialized')
                finally:snapshot.close()
                if args.extend_limit is not None:
                    from salmonpy.errors import SCFReviewRequired,InputError
                    short=replace(sa,max_iterations=2)
                    with session.run(a,scf=short,steps=1,material_id='ni-small-budget') as capped:
                        try:capped.continue_scf(steps=1)
                        except SCFReviewRequired:pass
                        if capped.iterations!=2:raise AssertionError('old cap was not reached')
                        boundary=frame(capped);counts=capped.scf_status
                        capped.extend_scf_limit(args.extend_limit)
                        compare(boundary,frame(capped))
                        if capped.scf_status['executed_iterations']!=counts['executed_iterations']:
                            raise AssertionError('terminal extension executed SCF')
                        try:capped.extend_scf_limit(args.extend_limit-1)
                        except InputError:pass
                        else:raise AssertionError('decreasing cap accepted')
                        capped.continue_scf(steps=2)
                        report['terminal_extension_continuation_errors']=compare(own['a'],frame(capped))
                        report['extended_limit']=capped.scf_status['max_iterations']

        else:
            label='b' if group==1 else 'a'
            with session.run(b if group==1 else a,scf=sb if group==1 else sa,
                             steps=4,material_id='ni-'+label) as state:
                own[label]=frame(state)
    if group==0:
        report['independent_b_errors']=compare(own['b'],world.recv(source=args.size+comm.rank,tag=31))
        report['independent_a_errors']=compare(own['a'],world.recv(source=2*args.size+comm.rank,tag=32))
        report['status']='passed'
    else:
        world.send(own['b' if group==1 else 'a'],dest=comm.rank,tag=31 if group==1 else 32)
    reports=world.gather(report if group==0 else None,root=0)
    if world.rank==0:
        if any(p.is_file() for p in args.output.rglob('data_for_restart/*')):
            raise AssertionError('GS data unexpectedly exported')
        result={'schema':'salmonpy.nickel-soc-scf-memory/1','status':'passed','solver_ranks':args.size,
                'kmesh':[2,2,2],'states':36,'electronic_temperature_k':300,'semicore_shells':['3s','3p'],'valence_electrons':18,
                'reports':[r for r in reports if r], 'gs_exported':False,
                'scope':'paused memory-only SCF snapshot and independent states, not physical convergence or disk restart'}
        (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
        print('Ni SOC full-SCF memory restore PASS',args.output/'result.json')
    comm.Free()


if __name__=='__main__':
    main()
