"""Native E/H group numerical retry: BUILD NEW_OUTPUT [--baseline|--faults].

Use 1/2 ranks, --split for two independent two-rank groups, OMP/BLAS1 and
external timeout. Fault hooks must be compiled separately, never in production.
"""
import ctypes as c
import json
from pathlib import Path
from types import SimpleNamespace
import sys
import numpy as np
from mpi4py import MPI
from salmonpy.coupled import NativeEHField
from salmonpy.fdtd import FDTD,PlaneWave
from salmonpy.fields import AcosPulse
from salmonpy.geometry import Grid,Geometry,Shape
from salmonpy.ld import LorentzDrude
from salmonpy.errors import NumericalError


def main():
    build,destination=sys.argv[1:3]
    split='--split' in sys.argv;faults='--faults' in sys.argv
    epsilon=float(sys.argv[sys.argv.index('--epsilon-infinity')+1]) if '--epsilon-infinity' in sys.argv else 1.0
    resume_midpoint='--resume-midpoint' in sys.argv
    if resume_midpoint and not '--baseline' in sys.argv:raise ValueError('resume-midpoint uses baseline mode')
    no_chdir='--no-chdir' in sys.argv
    pending_close='--pending-on-close' in sys.argv;saved_close='--saved-on-close' in sys.argv
    inspect_mpi='--inspect-mpi' in sys.argv
    inspect_close='--inspect-close' in sys.argv or pending_close or saved_close or inspect_mpi
    stop=int(sys.argv[sys.argv.index('--stop-at')+1]) if '--stop-at' in sys.argv else 32
    if not 1<=stop<=32 or (pending_close and saved_close):raise ValueError('invalid close probe configuration')
    stride=int(sys.argv[sys.argv.index('--observation-stride')+1]) if '--observation-stride' in sys.argv else 1
    world=MPI.COMM_WORLD
    if split:assert world.size==4
    else:assert world.size in (1,2)
    group=world.rank//2 if split else 0
    comm=world.Split(group,world.rank) if split else world
    baseline='--baseline' in sys.argv or group==1
    out=Path(destination).resolve()
    if split:out=out/f'group-{group}'
    if comm.rank==0:out.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    ctx=SimpleNamespace(solver_comm=comm,is_leader=comm.rank==0,closed=False,_dependents=0,
        group=SimpleNamespace(role='field'),plan=SimpleNamespace(point_ids=(1,)))
    anchor=Path.cwd()
    from unittest.mock import patch
    guard=patch('os.chdir',side_effect=AssertionError('FDTD host chdir forbidden')) if no_chdir else None
    if guard:guard.start()
    mpi_before=None;mpi_initialized=None;mpi_after_native=None;mpi_after_detach=None;host_mpi=None
    if inspect_mpi:
        import hashlib
        receipt=json.loads(Path(build).read_text())
        assert hashlib.sha256(Path(receipt['library']).read_bytes()).hexdigest()==receipt['sha256']
        probe_lib=c.CDLL(receipt['library'],mode=c.RTLD_GLOBAL)
        mpi_query=getattr(probe_lib,'sp_eh_mpi_resources',None)
        old_query=getattr(probe_lib,'sp_tddft_resources',None)
        def mpi_counts():
            values=np.zeros(3 if mpi_query else 4,dtype=np.int32)
            fn=mpi_query or old_query;fn.argtypes=[c.c_void_p];fn.restype=c.c_int
            assert fn(values.ctypes.data)==0
            return values[:3].tolist() if mpi_query else [int(values[0]),int(values[1]),None]
        mpi_before=mpi_counts();assert mpi_before[:2]==[0,0]
        if '--host-mpi-library' in sys.argv:
            directory=Path(sys.argv[sys.argv.index('--host-mpi-library')+1])
            hr=json.loads((directory/'result.json').read_text())
            assert hr['runtime_sha256']==receipt['sha256']
            assert hashlib.sha256(Path(hr['library']).read_bytes()).hexdigest()==hr['sha256']
            host_mpi=c.CDLL(hr['library']).test_host_mpi;host_mpi.argtypes=[c.c_int];host_mpi.restype=c.c_int
    grid=Grid((40,40,40),(2.,2.,2.));dt=.003;steps=32
    ld=(LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.],epsilon_infinity=epsilon),
        LorentzDrude(plasma_frequency=15.,strength=[.7,.3],resonance=[25.,35.],damping=[2.,3.],epsilon_infinity=1.0 if epsilon==1.0 else epsilon+.75))
    geometry=Geometry(grid,(Shape('box',(2.,4.,4.),(7.,1.,1.),material=1),Shape('box',(2.,4.,4.),(9.,1.,1.),material=2)))
    calc=FDTD(geometry,ld,dt,steps,source=PlaneWave(AcosPulse(.001,20.,.15),0,-16.),observations=((1.,1.,1.),),observation_stride=stride)
    admission_cases=[];geometry_mpi=[]
    original_call=NativeEHField._call
    def admission_call(self,function,*args):
        if getattr(function,'__name__','')=='sp_eh_init':
            path,count,indices=args
            trials=[('rank-local-zero-points',path,0 if comm.rank==0 else count),
                    ('rank-local-missing-input',b'nonexistent-admission-input' if comm.rank==0 else path,count),
                    ('rank-local-long-path',b'x'*4096 if comm.rank==0 else path,count)]
            for name,trial_path,trial_count in trials:
                code=original_call(self,function,trial_path,trial_count,indices)
                assert comm.allgather(int(code))==[1]*comm.size
                if inspect_mpi:assert mpi_counts()[0]==0
                admission_cases.append(name)
            if '--geometry-admission' in sys.argv:
                import re
                text=(self.directory/'inputfile').read_text()
                specifications=[
                    ('grid-overdetermined',text.replace('&maxwell',"&maxwell\n num_rgrid_em=40,40,40",1),1101),
                    ('boundary-conflict',text.replace("boundary_em(1,1) = 'abc'","boundary_em(1,1) = 'periodic'",1),1102),
                    ('time-overdetermined',text.replace('&maxwell',"&maxwell\n at_em=1",1),1103),
                    ('cfl',re.sub(r'\bdt_em\s*=\s*[-+0-9.eEdD]+','dt_em = 1.0',text),1104),
                ]
                for name,trial,expected in specifications:
                    assert trial!=text,name
                    for attempt in range(2):
                        if comm.rank==0:(self.directory/'geometry-rejected-input').write_text(trial)
                        comm.Barrier()
                        code=original_call(self,function,b'geometry-rejected-input',count,indices)
                        assert comm.allgather(int(code))==[expected]*comm.size,(name,code)
                        flags=np.ones(3,dtype=np.int32)
                        fn=self.lib.sp_eh_resources;fn.argtypes=[c.c_void_p];fn.restype=c.c_int
                        assert fn(flags.ctypes.data)==0 and flags.tolist()==[0,0,0]
                        if inspect_mpi:
                            counts=mpi_counts()
                            assert counts[0]==0 and (not mpi_query or counts[2]==0),counts
                            if attempt:assert counts[1]==geometry_mpi[-1]['counts'][1],counts
                            geometry_mpi.append({'case':name,'attempt':attempt,'counts':counts})
                        admission_cases.append('geometry-'+name+'-'+str(attempt))
            if '--parsed-admission' in sys.argv:
                rejected=self.directory/'restart-rejected-input'
                if comm.rank==0:
                    text=(self.directory/'inputfile').read_text()
                    rejected.write_text(text.replace('&control',"&control\n yn_restart = 'y'",1))
                comm.Barrier()
                for attempt in range(2):
                    code=original_call(self,function,b'restart-rejected-input',count,indices)
                    assert comm.allgather(int(code))==[1]*comm.size
                    if inspect_mpi:assert mpi_counts()[0]==0
                    admission_cases.append('parsed-restart-rejected-'+str(attempt))
            if '--malformed-admission' in sys.argv:
                rejected=self.directory/'malformed-rejected-input'
                if comm.rank==0:
                    text=(self.directory/'inputfile').read_text()
                    rejected.write_text(text.replace('&maxwell',"&maxwell\n nt_em = 'not-an-integer'",1))
                comm.Barrier()
                for attempt in range(2):
                    code=original_call(self,function,b'malformed-rejected-input',count,indices)
                    assert all(value!=0 for value in comm.allgather(int(code)))
                    if inspect_mpi:assert mpi_counts()[0]==0
                    admission_cases.append('malformed-integer-rejected-'+str(attempt))
            if '--unit-admission' in sys.argv:
                rejected=self.directory/'unit-rejected-input'
                if comm.rank==0:
                    text=(self.directory/'inputfile').read_text()
                    rejected.write_text(text.replace("unit_system = 'au'", "unit_system = 'unsupported-units'", 1))
                comm.Barrier()
                for attempt in range(2):
                    code=original_call(self,function,b'unit-rejected-input',count,indices)
                    assert comm.allgather(int(code))==[1001]*comm.size
                    if inspect_mpi:assert mpi_counts()[0]==0
                    admission_cases.append('invalid-units-rejected-'+str(attempt))
            if '--coordinate-admission' in sys.argv:
                specifications=[
                    ('missing-source', "natom=1", 1003),
                    ('duplicate-source', "natom=1, file_atom_coor='missing.xyz', file_atom_red_coor='missing.red'", 1003),
                    ('reduced-isolated', "natom=1, file_atom_red_coor='missing.red'", 1004),
                    ('missing-file', "natom=1, file_atom_coor='missing.xyz'", 1005),
                    ('bad-row', "natom=1, file_atom_coor='bad-coordinates.txt'", 1006),
                    ('truncated', "natom=2, file_atom_coor='short-coordinates.txt'", 1006),
                    ('bad-reduced', "natom=1, file_atom_red_coor='bad-coordinates.txt'", 1006),
                    ('bad-optimization-row', "natom=1, file_atom_coor='short-coordinates.txt'", 1006)]
                for name,settings,expected in specifications:
                    if comm.rank==0:
                        text=(self.directory/'inputfile').read_text().replace('&system','&system\n '+settings,1)
                        if name=='bad-reduced':text=text.replace("yn_periodic = 'n'", "yn_periodic = 'y'",1)
                        if name=='bad-optimization-row':text=text.replace('&calculation',"&calculation\n yn_opt='y'",1)
                        (self.directory/'coordinate-rejected-input').write_text(text)
                        (self.directory/'bad-coordinates.txt').write_text('Si invalid 0 0 1\n')
                        (self.directory/'short-coordinates.txt').write_text('Si 0 0 0 1\n')
                    comm.Barrier()
                    code=original_call(self,function,b'coordinate-rejected-input',count,indices)
                    assert comm.allgather(int(code))==[expected]*comm.size,(name,code)
                    if inspect_mpi:assert mpi_counts()[0]==0
                    admission_cases.append('coordinates-'+name)
            if '--semantic-admission' in sys.argv:
                specifications=[
                    ('yn-maxwell', "maxwell:yn_wf_em='z'",1013),
                    ('yn-shape', "maxwell:yn_make_shape='z'",1013)]
                for name,settings,expected in specifications:
                    if comm.rank==0:
                        text=(self.directory/'inputfile').read_text()
                        if ':' in settings:
                            group,entry=settings.split(':',1)
                            if name=='yn-shape':text=text.replace("yn_make_shape = 'y'","yn_make_shape = 'z'",1)
                            else:text=text.replace('&'+group,'&'+group+'\n '+entry,1)
                        else:text+='\n'+settings+'\n'
                        (self.directory/'semantic-rejected-input').write_text(text)
                    comm.Barrier()
                    code=original_call(self,function,b'semantic-rejected-input',count,indices)
                    assert comm.allgather(int(code))==[expected]*comm.size,(name,code)
                    if inspect_mpi:assert mpi_counts()[0]==0
                    admission_cases.append('semantic-'+name)
                if comm.rank==0:
                    text=(self.directory/'inputfile').read_text()
                    if '--inactive-electronic-settings' in sys.argv:
                        text=text.replace('&system','&system\n nelem=1',1)
                        text+=("\n&pseudo\n file_pseudo(1)='missing-electronic-only_rps.dat'\n/\n"
                               "&functional\n hse_mlwf_interval=11, exx_mlwf_interval=12,\n"
                               " hse_mlwf_maxiter=13, exx_mlwf_maxiter=14,\n"
                               " hse_mlwf_tolerance=1d-7, exx_mlwf_tolerance=2d-7,\n"
                               " exx_k_algorithm='user_route', exx_factor_history_frames=0,\n"
                               " exx_factor_n_teach=0, exx_factor_exact_interval=0\n/\n")
                    (self.directory/'inputfile').write_text(
                        text.replace('&maxwell',"&maxwell\n yn_wf_em='N'",1)
                        +"\n&propagation\n propagator='user_extension'\n/\n"
                        +"&scf\n method_init_wf='user_initializer'\n/\n")
                comm.Barrier()


        return original_call(self,function,*args)
    with patch.object(NativeEHField,'_call',admission_call) if '--geometry-admission' in sys.argv or '--admission' in sys.argv or '--parsed-admission' in sys.argv or '--malformed-admission' in sys.argv or '--unit-admission' in sys.argv or '--coordinate-admission' in sys.argv or '--semantic-admission' in sys.argv else patch.object(NativeEHField,'_call',original_call):
        f=NativeEHField(ctx,calc,build,out/'field',points={1:(20,20,20)},grid=(1,comm.size,1))
    if inspect_mpi:
        mpi_initialized=mpi_counts();assert mpi_initialized[0]>0
        if mpi_query:assert mpi_initialized[2]==mpi_initialized[0]
        if host_mpi:assert host_mpi(0)==0
    if not baseline:
        f.lib.sp_eh_recovery_protocol.restype=c.c_int
        assert f.lib.sp_eh_recovery_protocol()==1
        assert hasattr(f.lib,'sp_eh_test_fault')==faults
    cp=f.lib.sp_eh_checkpoint;cp.argtypes=[c.c_char_p,c.c_int];cp.restype=c.c_int
    def save(tag):
        p=out/f'{tag}-{comm.rank}.bin';assert f._call(cp,str(p).encode(),0)==0;return p.read_bytes()
    def files():return {p.name:p.read_bytes() for p in f.directory.glob('*.data')}
    cases=[];rows=[]
    if not baseline and comm.size>1:
        before=save('mixed-entry-before')
        E=np.empty_like(f.A);J=np.zeros_like(f.A)
        floats=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
        f.lib.sp_eh_step_checked.argtypes=[c.c_int,floats,floats,floats]
        f.lib.sp_eh_step_checked.restype=c.c_int
        if comm.rank==0:code=f._call(f.lib.sp_eh_step_checked,1,J,f.A,E)
        else:code=f._call(f.lib.sp_eh_step,1,J,E)
        assert code==11
        assert save('mixed-entry-after')==before
        cases.append(('admission','mixed-entry'))
    for step in range(1,stop+1):
        J=np.array([[0.,.002*np.sin(.2*step),0.]]) if step>1 else np.zeros((1,3))
        if not baseline and step==4:
            before=save('scalar-admission-before');outputs=files()
            for label,bad in [('nan',float('nan')),('positive-inf',float('inf')),('negative-inf',-float('inf'))]:
                for which in ('J','A'):
                    trial_J=J.copy();trial_A=f.A.copy()
                    if comm.rank==comm.size-1:
                        (trial_J if which=='J' else trial_A)[0,2]=bad
                    E=np.empty_like(f.A)
                    assert f._call(f.lib.sp_eh_step_checked,1,trial_J,trial_A,E)==10
                    assert save('scalar-'+which+'-'+label+'-after')==before
                    assert files()==outputs and f.step_index==step-1
                    cases.append(('scalar-admission',which,label))
        if not baseline and step in (1,12,32):
            attempts=list(range(20,24)) if faults else []
            for stage in attempts:
                A=f.A.copy();E=f.E.copy();before=save(f'{step}-{stage}-before');outputs=files()
                fn=f.lib.sp_eh_test_fault;fn.argtypes=[c.c_int,c.c_int];fn.restype=c.c_int
                assert fn(stage,comm.size-1)==0
                try:f.advance(J)
                except NumericalError:pass
                else:raise AssertionError('numerical error not returned')
                assert not f.failed and f.step_index==step-1
                np.testing.assert_array_equal(f.A,A);np.testing.assert_array_equal(f.E,E)
                assert save(f'{step}-{stage}-after')==before
                assert files()==outputs
                cases.append((step,stage))
            if step==12:
                # Real finite-input overflow of A-dt*sample(E), before ordinary output.
                A=f.A.copy();E=f.E.copy();f.A[:]=np.finfo(float).max
                before=save('integrated-A-before');outputs=files()
                try:f.advance(np.full((1,3),np.finfo(float).max))
                except NumericalError:pass
                else:raise AssertionError('integrated A overflow not returned')
                assert not f.failed and f.step_index==step-1
                assert save('integrated-A-after')==before and files()==outputs
                np.testing.assert_array_equal(f.E,E);f.A=A
                cases.append((step,'integrated-A-overflow'))
        rows.append(np.r_[f.advance(J).ravel(),f.E.ravel(),f.sample_magnetic_field().values.ravel()])
        if resume_midpoint and step==16:
            save('midpoint')
            midpoint=dict(A=f.A.copy(),E=f.E.copy(),step_index=f.step_index)
        if resume_midpoint and step==18:
            from salmonpy._native_checkpoint import NativeRestoreCandidate
            candidate=NativeRestoreCandidate(f,'sp_eh_checkpoint',out/f'midpoint-{comm.rank}.bin')
            with patch.object(f,'advance',side_effect=AssertionError('restoration must not advance')):
                candidate.stage();candidate.publish()
            for name,value in midpoint.items():setattr(f,name,value.copy() if isinstance(value,np.ndarray) else value)
            for repeated in (17,18):
                incoming=np.array([[0.,.002*np.sin(.2*repeated),0.]])
                actual=np.r_[f.advance(incoming).ravel(),f.E.ravel(),f.sample_magnetic_field().values.ravel()]
                np.testing.assert_array_equal(actual,rows[repeated-1])

    if not baseline and f.step_index==steps:
        before=save('capacity-before');E=np.empty_like(f.A)
        code=f._call(f.lib.sp_eh_step,1,np.zeros_like(f.A),E)
        assert code==13 and save('capacity-after')==before
        cases.append(('admission','capacity'))
    np.save(out/f'fields-{comm.rank}.npy',rows)
    allocation_before=None;allocation_after=None
    if inspect_close:
        resources=getattr(f.lib,'sp_eh_resources',None)
        if resources is None:raise RuntimeError('matching native field-resource diagnostics required')
        resources.argtypes=[c.c_void_p];resources.restype=c.c_int
        def flags():
            values=np.zeros(3,dtype=np.int32);assert resources(values.ctypes.data)==0
            return values.tolist()
        assert flags()==[1,0,0]
        if pending_close:
            save('close-candidate')
            outputs=files()
            assert f._call(cp,str(out/f'close-candidate-{comm.rank}.bin').encode(),2)==0
            assert files()==outputs
        if saved_close:
            transaction=f.lib.sp_eh_step_transaction;transaction.argtypes=[c.c_int];transaction.restype=c.c_int
            assert f._call(transaction,0)==0
        allocation_before=flags()
        assert allocation_before==[1,int(saved_close),int(pending_close)]
        outputs_before_close=files()
    if inspect_mpi:
        original_call=f._call
        def watched(fn,*args):
            nonlocal mpi_after_native
            status=original_call(fn,*args)
            if getattr(fn,'__name__','')=='sp_eh_close' and status==0:
                mpi_after_native=mpi_counts()
                (out/f'mpi-native-close-{comm.rank}.json').write_text(json.dumps(mpi_after_native))
                assert mpi_after_native[0]==int(host_mpi is not None),'field children survived native close'
                if mpi_query:assert mpi_after_native[2]==0
                comm.Barrier()  # Caller parent must remain valid before detach.
                if host_mpi:
                    assert host_mpi(1)==0 and mpi_counts()[0]==0
            return status
        f._call=watched
    f.close()
    if inspect_mpi:
        mpi_after_detach=mpi_counts();assert mpi_after_detach[:2]==[0,0]
        comm.Barrier()
    if inspect_close:
        allocation_after=flags();assert allocation_after==[0,0,0]
        assert f.lib.sp_eh_close()==1 and flags()==[0,0,0]
        if stop<steps:assert files()==outputs_before_close
    assert Path.cwd()==anchor
    if guard:guard.stop()
    if comm.rank==0:
        (out/'validation.json').write_text(json.dumps(dict(status='passed',epsilon_infinity=[epsilon,1.0 if epsilon==1.0 else epsilon+.75],checkpoint_restore_without_reinitialization=resume_midpoint,validation_continuation_repeated_steps=2 if resume_midpoint else 0,new_runtime_no_chdir=no_chdir,admission_cases=admission_cases,geometry_mpi=geometry_mpi,cases=cases,ranks=comm.size,steps=stop,capacity=steps,replayed_steps=0,allocation_before_close=allocation_before,allocation_after_close=allocation_after,mpi_before=mpi_before,mpi_initialized=mpi_initialized,mpi_after_native_close=mpi_after_native,mpi_after_detach=mpi_after_detach,other_owner_preserved=bool(host_mpi)),indent=2)+'\n')
        print('E/H deep recovery PASS',cases)


if __name__=='__main__':main()
