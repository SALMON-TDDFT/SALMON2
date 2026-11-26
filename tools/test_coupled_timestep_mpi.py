"""Explicit TDDFT+E/H timestep study; never run as required physics CI.

build-analysis NEW_DIRECTORY: compile/audit the small Fortran analysis module.
run BUILD EXISTING_PZ_FIXTURE NEW_OUTPUT --dt .003 --steps 5334:
launch one existing MPI job with --field-ranks2 plus two electronic k ranks.
--dt/--steps refer to ELECTRONIC updates. --field-steps M gives field dt=dt/M.
Compare field-only and electronic-only refinement at identical duration. Existing GS is
referenced through symlinks, not copied or calculated. No electronic checkpoint.
compare COARSE FINE ANALYSIS_BUILD NEW_JSON: read completed runs, never start MPI.
"""
import argparse
import ctypes as c
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import numpy as np

NAMES=('incoming_energy','outgoing_energy','field_modified_change','field_physical_change',
       'field_current_work','midpoint_electronic_work','electronic_energy_at_field_end',
       'field_modified_residual','field_physical_residual','total_modified_residual','total_physical_residual',
       'coupling_work_defect','electronic_work_defect','electron_loss_fraction','electron_gain_fraction',
       'current_peak','electric_peak','initial_cell_energy','electronic_energy_at_endpoint','late_current_peak',
       'field_modified_relative_incident','field_physical_relative_incident',
       'total_modified_relative_incident','total_physical_relative_incident',
       'coupling_defect_relative_work','electronic_defect_relative_work','max_step_field_residual','net_outward_energy','staggered_point_work','field_point_work_mismatch','staggered_work_endpoint_term','staggered_work_identity_residual','held_current_work_term')


def package_root():
    root=Path(__file__).resolve().parents[1]
    return root


def build_analysis(output):
    output=output.resolve();output.mkdir(exist_ok=False)
    root=package_root();source=root/'tests/native/test_coupled_timestep.f90'
    checker=root/'tools/check_forbidden_eigensolvers.py'
    subprocess.run([sys.executable,str(checker),'--source-root',str(source.parent)],check=True)
    obj=output/'analysis.o'
    subprocess.run(['gfortran','-std=f2008','-O2','-fPIC','-J'+str(output),'-c',str(source),'-o',str(obj)],check=True)
    subprocess.run([sys.executable,str(checker),'--object-tree',str(output)],check=True)
    library=output/('analysis.dylib' if sys.platform=='darwin' else 'analysis.so')
    subprocess.run(['gfortran','-dynamiclib' if sys.platform=='darwin' else '-shared',str(obj),'-o',str(library)],check=True)
    lib=c.CDLL(str(library));lib.study_self_test.restype=c.c_int
    assert lib.study_self_test()==0,'Fortran timestep/energy self-test failed'
    (output/'build.json').write_text(json.dumps({'schema':'salmonpy.coupled-timestep-analysis/2','library':str(library),
        'sha256':hashlib.sha256(library.read_bytes()).hexdigest(),
        'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'source_object_audits':'passed','finite_NaN_Inf_and_clock_tests':'passed'},indent=2)+'\n')


class ExistingGSProvider:
    def __init__(self,fixture,output,steps,dt,n_hamil):
        self.fixture=fixture.resolve();self.output=output.resolve();self.steps=steps;self.dt=dt;self.n_hamil=n_hamil
        self.identity='existing-PZ-Si/'+hashlib.sha256((self.fixture/'rt.inp').read_bytes()).hexdigest()
    def prepare(self,context):
        from salmonpy.bands import _input_groups
        from salmonpy._hse_geometry import _numeric
        groups=_input_groups((self.fixture/'rt.inp').read_text())
        # This is deliberately a fixed Si regression, not a general material runner.
        for axis in range(3):
            expected=[0.,0.,0.];expected[axis]=10.26
            if not np.array_equal(_numeric(groups['system'],'al_vec'+str(axis+1),3),expected):
                raise ValueError('this fixed Si audit requires explicit10.26bohr cubic cell vectors')
        work=self.output/'electrons';comm=context.solver_comm
        if comm.rank==0:
            work.mkdir(exist_ok=False)
            groups=_input_groups((self.fixture/'rt.inp').read_text())
            assert "xc='pz'" in groups['functional'].replace(' ','')
            groups['tgrid']=f'\n dt={self.dt:.17g}, nt={self.steps}\n'
            groups['analysis']='\n out_rt_energy_step=1\n'
            groups['propagation']=groups.get('propagation','')+f'\n n_hamil={self.n_hamil}\n'
            (work/'rt.inp').write_text('\n'.join('&'+k+v+'/' for k,v in groups.items())+'\n')
            for name in ('data_for_restart','pseudo1_rps.dat'):
                (work/name).symlink_to(self.fixture/name,target_is_directory=name=='data_for_restart')
        comm.Barrier();return work


class ExistingGSBackend:
    """Manual fixture using real shared state handles without GS work copies."""
    execution='collective';reduction='leader'
    def __init__(self,build,steps):
        self.record=json.loads(build.read_text());self.steps=steps
        self.identity='shared-PZ/'+self.record['sha256']
    def open_group(self,work,dt,context):
        from types import SimpleNamespace
        from salmonpy.distributed import NativeParallelContext
        from salmonpy._shared_tddft import NativeTDDFTModel,_configure
        owner=SimpleNamespace(solver_comm=context.solver_comm,closed=False,_dependents=0)
        native=NativeParallelContext(owner,self.record['library'],self.record['sha256'])
        lib=c.CDLL(self.record['library']);_configure(lib)
        return NativeTDDFTModel(native,lib,work,
            SimpleNamespace(max_steps=self.steps,manifest=(),input_file='rt.inp'),dt,
            (context.solver_comm.size,1,1,1,1),self.identity,Path('.'))


def run(a):
    from mpi4py import MPI
    from salmonpy.distributed import ParallelPlan,SolverGroup,ParallelContext
    from salmonpy.material_response import ResponseSession,MaterialBinding
    from salmonpy.coupled import NativeEHField,FieldSubcycledStepper
    from salmonpy.fdtd import FDTD,PlaneWave,Boundaries
    from salmonpy.fields import AcosPulse
    from salmonpy.geometry import Grid,Geometry
    from salmonpy.native_steps import _NativeBackup
    from salmonpy.examples.molecular_fdtd_mpi_response import abort_on_error
    world=MPI.COMM_WORLD;out=a.output.resolve();start=time.perf_counter()
    if world.size!=a.field_ranks+2:raise ValueError('field ranks plus two TDDFT k ranks required')
    if a.steps<2 or a.dt<=0:raise ValueError('positive dt and at least two updates required')
    if world.rank==0:out.mkdir(exist_ok=False)
    world.Barrier()
    plan=ParallelPlan((SolverGroup('field',tuple(range(a.field_ranks)),role='field'),
                       SolverGroup('si',(a.field_ranks,a.field_ranks+1),point_ids=(10,))))
    provider=ExistingGSProvider(a.fixture,out,a.steps,a.dt,a.n_hamil)
    binding=MaterialBinding('si',provider,ExistingGSBackend(a.build,a.steps))
    history=[];norms=[];initial_current=None;quantum_current=[];token=None
    field_dt=a.dt/a.field_steps;field_count=a.steps*a.field_steps
    with ParallelContext.create(world,plan) as context,abort_on_error(world):
        field=None;response=None
        try:
            if context.group.role=='field':
                calculation=FDTD(Geometry(Grid((40,40,40),(a.spacing,)*3)),(),field_dt,field_count,
                    boundaries=Boundaries(('pml','pml','periodic','periodic','periodic','periodic')),
                    source=PlaneWave(AcosPulse(a.amplitude,a.omega,a.pulse_duration),0,-8*a.spacing))
                field=NativeEHField(context,calculation,a.build,out/'field',points={10:(20,20,20)},grid=(a.field_ranks,1,1))
                observer=field.energy_observer((13,13,13),(27,27,27))
                backup=_NativeBackup(field,context.solver_comm,'sp_eh_step_transaction',field.directory,
                                     ('A','E','step_index','failed','closed'))
            response=ResponseSession(context,{10:binding},dt=a.dt)
            stepper=FieldSubcycledStepper(response,field,field_steps=a.field_steps)
            if context.group.role!='field':
                state=response.states[10];norms.append(state.electron_number)
                initial_current=state.initial_current.tolist()
            def before_field(field):
                nonlocal token
                token=backup.capture()
            def after_field(field):
                budget=observer.sample_step();error=None
                if context.is_leader:
                    try:
                        if not budget['available']:raise RuntimeError(budget)
                        history.append({'E':field.E[0].tolist(),'J_used':stepper.J[0].tolist(),
                            'budget':[budget[name] for name in observer._names]})
                    except Exception as exc:error=str(exc)
                errors=context.solver_comm.allgather(error)
                if any(errors):raise RuntimeError(str(errors))
                backup.release(token)
            for n in range(a.steps):
                frame=stepper.step(before_field_step=before_field,after_field_step=after_field)
                if field is not None and context.is_leader:quantum_current.append(frame['J_charge'][0].tolist())
                if field is None:norms.append(state.electron_number)
                if world.rank==0 and (n+1)%512==0:print(f'accepted electronic {n+1}/{a.steps}; field {(n+1)*a.field_steps}/{field_count}',flush=True)
        except BaseException:
            # Abort before collective cleanup if a rank-local observation failed.
            import traceback
            traceback.print_exc();world.Abort(1)
            raise
        finally:
            if response:response.close()
            if field:field.close()
    records=world.gather({'norms':norms,'initial_current':initial_current} if world.rank==a.field_ranks else None,root=0)
    if world.rank==0:
        q=records[a.field_ranks]
        table=np.loadtxt(next((out/'electrons').glob('*_rt_energy.data')))
        assert len(table)==a.steps+1
        np.savez(out/'history.npz',E=np.array([v['E'] for v in history]),
            J=np.array([q['initial_current']]+quantum_current),
            J_used=np.array([v['J_used'] for v in history]),
            budget=np.array([v['budget'] for v in history]),
            quantum_energy=table[:,1],electron_number=np.array(q['norms']))
        (out/'run.json').write_text(json.dumps({'dt':field_dt,'steps':field_count,'field_dt':field_dt,'electronic_dt':a.dt,
            'electronic_steps':a.steps,'field_steps_per_electronic_step':a.field_steps,
            'current_policy':'hold-last-accepted','n_hamil':a.n_hamil,'electronic_time_origin':.5*field_dt,'field_end_time':a.dt*a.steps,
            'pulse_duration':a.pulse_duration,'omega':a.omega,'amplitude':a.amplitude,
            'field_counts':[40]*3,'field_spacing':[a.spacing]*3,'point':[20]*3,'box':[[13]*3,[27]*3],
            'macro_volume':a.spacing**3,'cell_volume':10.26**3,'field_ranks':a.field_ranks,'electronic_ranks':2,
            'wall_seconds':time.perf_counter()-start,'GS_copied':False,'SCF_runs':0,
            'library_sha256':json.loads(a.build.read_text())['sha256'],
            'MPI_alive':not MPI.Is_finalized(),'fixture':str(a.fixture.resolve())},indent=2)+'\n')
        print('RUN COMPLETE',out,flush=True)


def compare(a):
    rec=json.loads((a.analysis/'build.json').read_text());path=Path(rec['library'])
    assert hashlib.sha256(path.read_bytes()).hexdigest()==rec['sha256']
    if rec.get('schema')!='salmonpy.coupled-timestep-analysis/2':raise ValueError('rebuild matching analysis protocol2 explicitly')
    lib=c.CDLL(str(path));ptr=c.c_void_p
    lib.study_reduce.argtypes=[c.c_int,c.c_int,c.c_double,c.c_double,c.c_double]+[ptr]*7
    lib.study_reduce.restype=c.c_int
    lib.study_compare.argtypes=[c.c_int]*4+[c.c_double]*4+[ptr]*5;lib.study_compare.restype=c.c_int
    runs=[]
    for directory in (a.coarse,a.fine):
        meta=json.loads((directory/'run.json').read_text());data=dict(np.load(directory/'history.npz'))
        from salmonpy.bands import _input_groups
        from salmonpy._hse_geometry import _numeric
        actual=_input_groups((directory/'electrons/rt.inp').read_text())
        meta['actual_n_hamil']=int(_numeric(actual['propagation'],'n_hamil',1,True,default=4)[0])
        arrays=[np.ascontiguousarray(data[k]) for k in ('budget','E','J_used','J','quantum_energy','electron_number')]
        summary=np.empty(33)
        assert lib.study_reduce(meta['electronic_steps'],meta['field_steps_per_electronic_step'],meta['electronic_dt'],meta['macro_volume']/meta['cell_volume'],meta['macro_volume'],
            *(v.ctypes.data for v in arrays),summary.ctypes.data)==0
        runs.append((meta,data,dict(zip(NAMES,map(float,summary)))))
    coarse,fine=runs
    assert coarse[0]['field_end_time']==fine[0]['field_end_time']
    for k in ('pulse_duration','omega','amplitude','field_counts','field_spacing','point','box','macro_volume','cell_volume','library_sha256','fixture','actual_n_hamil'):
        assert coarse[0][k]==fine[0][k],k
    arrays=[np.ascontiguousarray(d[k]) for d,k in ((coarse[1],'J'),(fine[1],'J'),(coarse[1],'E'),(fine[1],'E'))]
    stats=np.empty(6)
    assert lib.study_compare(coarse[0]['electronic_steps'],fine[0]['electronic_steps'],coarse[0]['steps'],fine[0]['steps'],
        coarse[0]['electronic_dt'],fine[0]['electronic_dt'],coarse[0]['field_dt'],fine[0]['field_dt'],
        *(v.ctypes.data for v in arrays),stats.ctypes.data)==0
    lib.study_difference.argtypes=[c.c_double,c.c_double,ptr];lib.study_difference.restype=c.c_int
    energy_comparison={}
    for key in ('field_current_work','electronic_energy_at_field_end'):
        pair=np.empty(2);status=lib.study_difference(coarse[2][key],fine[2][key],pair.ctypes.data)
        if status==1:energy_comparison[key]={'available':False,'reason':'zero fine reference'}
        else:
            assert status==0
            energy_comparison[key]={'absolute_difference':float(pair[0]),'relative_difference':float(pair[1])}
    result={'coarse':coarse[2],'fine':fine[2],'energy_comparison':energy_comparison,
        'comparison':dict(zip(('current_rms_abs','current_rms_relative','current_max_abs','E_rms_abs','E_rms_relative','E_max_abs'),map(float,stats))),
        'time_alignment':'E at common field times; electronic J at n*dt_e+dt_field/2 explicitly interpolated in Fortran; no extrapolation',
        'timesteps':{'coarse':{k:coarse[0][k] for k in ('field_dt','electronic_dt','field_steps_per_electronic_step','actual_n_hamil')},
                     'fine':{k:fine[0][k] for k in ('field_dt','electronic_dt','field_steps_per_electronic_step','actual_n_hamil')}},
        'field_end_time':coarse[0]['field_end_time'],
        'advisory':{'current_relative_RMS_target':.05,'electron_loss_fraction_target':1e-5,
            'current_below_target':float(stats[1])<=.05,
            'both_norms_below_target':all(run[2]['electron_loss_fraction']<=1e-5 for run in runs),
            'current_has_signal':fine[2]['current_peak']>0.,
            'note':'targets report diagnostics, never stop valid propagation'},
        'order_claim':'two levels measure a difference, not an observed convergence order',
        'scope':'fixed Si seed, temporal coupling check; spatial/k/physical convergence not certified'}
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='mode',required=True)
    q=sub.add_parser('build-analysis');q.add_argument('output',type=Path)
    q=sub.add_parser('run')
    for k in ('build','fixture','output'):q.add_argument(k,type=Path)
    q.add_argument('--dt',type=float,default=.003,help='electronic time step');q.add_argument('--steps',type=int,default=5334,help='electronic updates')
    q.add_argument('--field-steps',type=int,default=1,help='field updates per electronic update; field dt=dt/field-steps')
    q.add_argument('--n-hamil',type=int,default=4,help='Taylor expansion order in this fixed PZ midpoint fixture')
    q.add_argument('--amplitude',type=float,default=.001,help='pulse amplitude;0 is an unforced electronic stability control')
    q.add_argument('--spacing',type=float,default=2.,help='Maxwell grid spacing in bohr; does not change electronic grid')
    q.add_argument('--field-ranks',type=int,default=2);q.add_argument('--omega',type=float,default=.5)
    q.add_argument('--pulse-duration',type=float,default=8.)
    q=sub.add_parser('compare')
    for k in ('coarse','fine','analysis','output'):q.add_argument(k,type=Path)
    a=p.parse_args()
    if a.mode=='run':run(a)
    elif a.mode=='compare':compare(a)
    else:build_analysis(a.output)


if __name__=='__main__':main()
