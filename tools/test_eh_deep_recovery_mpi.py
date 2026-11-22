"""Native E/H group numerical retry: BUILD NEW_OUTPUT [--baseline|--faults].

Use 1/2/7 ranks; --auto-layout-admission uses7 for unsupported automatic
factorization and explicit-layout retry. --split for two independent two-rank groups, OMP/BLAS1 and
external timeout. --tensor exercises reciprocal off-diagonal LD using the same
recovery/close suite; use --baseline --resume-midpoint on a production build.
Fault hooks must be compiled separately, never in production.
--host-unit-library DIR --host-unit-lifetime verifies caller-owned file units
through initialization, stepping, direct restore and finalization.
--observation-planes also exercises legacy plane snapshots/integrals.
--material-input-mode diagonal|tensor|subpixel and --material-input-admission
exercise reader rejection, public lifetime release and valid-input retry.
--coefficient-admission uses those modes for late cold-init1111/1112/1114 returns.
--allocation-admission injects partial core-3D allocation failures (test build only).
--ld-allocation-admission tests partial LD pole/current-history allocation failures.
--aux-allocation-admission tests material/boundary coefficient vectors.
--material-allocation-failures tests partial LD tensor-input array allocation failures.
--ld-coefficient-allocation-failures tests LD material-list and coefficient scratch failures.
--point-output-failure MODE --point-output-reference FILE checks optional point I/O and repair.
--integral-output-failure MODE uses the same faults for xy integrals; --integrals-only enables all integral controls.
--info-output-failure open|write tests optional final metadata failure and cleanup.
--snapshot-output-failure open|write tests plane snapshot failure,repair and direct continuation.
--checkpoint-save-failures checks rank-local required-save failure,state retention and repair.
--compact-checkpoint uses a smaller valid geometry for checkpoint-I/O controls,not API defaults.
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
from salmonpy.ld import LorentzDrude,TensorLorentzDrude
from salmonpy.errors import NumericalError


def main():
    build,destination=sys.argv[1:3]
    split='--split' in sys.argv;faults='--faults' in sys.argv
    material_mode=sys.argv[sys.argv.index('--material-input-mode')+1] if '--material-input-mode' in sys.argv else None
    if material_mode not in (None,'diagonal','tensor','subpixel'):raise ValueError('unknown material-input fixture')
    tensor='--tensor' in sys.argv or material_mode in ('diagonal','tensor')
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
    snapshot_failure_mode=sys.argv[sys.argv.index('--snapshot-output-failure')+1] if '--snapshot-output-failure' in sys.argv else None
    if snapshot_failure_mode not in (None,'open','write'):raise ValueError('unknown snapshot-output failure fixture')
    if snapshot_failure_mode and (not resume_midpoint or 18%stride):
        raise ValueError('snapshot recovery fixture requires midpoint resume and a stride dividing18')
    info_failure_mode=sys.argv[sys.argv.index('--info-output-failure')+1] if '--info-output-failure' in sys.argv else None
    if info_failure_mode not in (None,'open','write'):raise ValueError('unknown info-output failure fixture')
    integral_failure='--integral-output-failure' in sys.argv
    failure_flag='--integral-output-failure' if integral_failure else '--point-output-failure'
    point_failure_mode=sys.argv[sys.argv.index(failure_flag)+1] if failure_flag in sys.argv else None
    if point_failure_mode not in (None,'header-open','header-write','append-open','append-write'):
        raise ValueError('unknown observation-output failure fixture')
    point_reference=Path(sys.argv[sys.argv.index('--point-output-reference')+1]) if point_failure_mode else None
    point_output_name='obs1_xy_integral_rt.data' if integral_failure else 'obs1_at_point_rt.data'
    output_header_lines=5 if integral_failure else 4
    world=MPI.COMM_WORLD
    if split:assert world.size==4
    else:assert world.size in (1,2,7)
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
    compact='--compact-checkpoint' in sys.argv
    grid=Grid((24,24,24) if compact else (40,40,40),(2.,2.,2.));dt=.003;steps=32
    sample_points={1:(12,10,10) if compact else (20,20,20)}
    ld=(LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.],epsilon_infinity=epsilon),
        LorentzDrude(plasma_frequency=15.,strength=[.7,.3],resonance=[25.,35.],damping=[2.,3.],epsilon_infinity=1.0 if epsilon==1.0 else epsilon+.75))
    if tensor:
        cs,sn=np.cos(.4),np.sin(.4)
        rotation=np.eye(3) if material_mode=='diagonal' else np.array([[cs,-sn,0.],[sn,cs,0.],[0.,0.,1.]])
        ld=tuple(TensorLorentzDrude(plasma_frequency=m.config.plasma_frequency,
            strength=[rotation@np.diag([f,1.1*f,1.2*f])@rotation.T for f in m.config.strength],
            resonance=m.config.resonance,damping=m.config.damping,
            epsilon_infinity=rotation@np.diag([m.config.epsilon_infinity,1.2*m.config.epsilon_infinity,1.4*m.config.epsilon_infinity])@rotation.T)
            for m in ld)
    box_size=(2.,2.,2.) if compact else (2.,4.,4.)
    box_yz=(3.,3.) if compact else (1.,1.)
    geometry=Geometry(grid,(Shape('box',box_size,(-1. if compact else 7.,*box_yz),material=1),
        Shape('box',box_size,(1. if compact else 9.,*box_yz),material=2)))
    smoothing=None
    if material_mode=='subpixel':
        from salmonpy.subpixel import SubpixelSmoothing
        smoothing=SubpixelSmoothing(reference_omega=50.)
    calc=FDTD(geometry,ld,dt,steps,source=PlaneWave(AcosPulse(.001,20.,.15),0,-4. if compact else -16.),observations=((1.,1.,1.),),observation_stride=stride,smoothing=smoothing)
    admission_cases=[];geometry_mpi=[]
    host_unit=None
    if '--host-unit-library' in sys.argv:
        import hashlib
        receipt=json.loads((Path(sys.argv[sys.argv.index('--host-unit-library')+1])/'result.json').read_text())
        assert hashlib.sha256(Path(receipt['library']).read_bytes()).hexdigest()==receipt['sha256']
        host_lib=c.CDLL(receipt['library']);host_unit=host_lib.test_host_unit
        host_unit.argtypes=[c.c_char_p,c.c_int];host_unit.restype=c.c_int
        assert host_unit(str(out/f'host-{comm.rank}.data').encode(),0)==0
    original_call=NativeEHField._call
    def automatic_input(text):
        import re
        for key,value in [('nproc_k','0'),('nproc_ob','0'),('nproc_rgrid','0,0,0')]:
            text=re.sub(r'\b'+key+r'\s*=.*',key+' = '+value,text)
        return text
    def shape_text(text,name):
        import re
        text=text.replace("yn_make_shape = 'y'","yn_make_shape = 'n'",1)
        if re.search(r'\bshape_file\s*=',text):
            return re.sub(r'\bshape_file\s*=.*',"shape_file='"+name+"'",text)
        return text.replace('&maxwell',"&maxwell\n shape_file='"+name+"'",1)
    if '--shape-input' in sys.argv:
        shape_source=Path(sys.argv[sys.argv.index('--shape-input')+1]).read_bytes()
        inherited_shape_call=original_call
        def external_shape_call(self,function,*call_args):
            if getattr(function,'__name__','')=='sp_eh_init':
                if comm.rank==0 and call_args[0]==b'inputfile':
                    path=self.directory/'inputfile'
                    if 'shape_file' not in path.read_text():path.write_text(shape_text(path.read_text(),'external-shape.cube'))
                    (self.directory/'external-shape.cube').write_bytes(shape_source)
                comm.Barrier()
            return inherited_shape_call(self,function,*call_args)
        original_call=external_shape_call
    if '--observation-planes' in sys.argv or '--integrals-only' in sys.argv or integral_failure or snapshot_failure_mode:
        inherited_plane_call=original_call
        def observation_plane_call(self,function,*call_args):
            if getattr(function,'__name__','')=='sp_eh_init':
                if comm.rank==0 and call_args[0]==b'inputfile':
                    path=self.directory/'inputfile'
                    path.write_text(path.read_text().replace('&maxwell',
                        ("&maxwell\n yn_obs_plane_em(1)='y', yn_obs_plane_integral_em(1)='y'" if '--observation-planes' in sys.argv or snapshot_failure_mode else "&maxwell\n yn_obs_plane_integral_em(1)='y'"),1))
                comm.Barrier()
            return inherited_plane_call(self,function,*call_args)
        original_call=observation_plane_call
    if '--orbital-mapping' in sys.argv:
        inherited_call=original_call
        def implemented_mapping_call(self,function,*call_args):
            if getattr(function,'__name__','')=='sp_eh_init':
                if comm.rank==0 and call_args[0]==b'inputfile':
                    path=self.directory/'inputfile'
                    path.write_text(path.read_text().replace('&parallel',"&parallel\n process_allocation='orbital_sequential'",1))
                comm.Barrier()
            return inherited_call(self,function,*call_args)
        original_call=implemented_mapping_call
    if point_failure_mode and point_failure_mode.startswith('header'):
        inherited_point_call=original_call
        def point_io_init_call(self,function,*call_args):
            if getattr(function,'__name__','')=='sp_eh_init':
                if comm.rank==0:
                    path=self.directory/point_output_name
                    if point_failure_mode=='header-open':path.mkdir()
                    else:
                        path.write_bytes(b'');path.chmod(0o400)
                comm.Barrier()
            return inherited_point_call(self,function,*call_args)
        original_call=point_io_init_call
    def admission_call(self,function,*args):
        nonlocal host_unit
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
            if '--auto-layout-admission' in sys.argv:
                assert comm.size==7,'this unsupported-prime admission fixture uses seven ranks'
                if comm.rank==0:
                    (self.directory/'automatic-rejected-input').write_text(automatic_input((self.directory/'inputfile').read_text()))
                comm.Barrier()
                for attempt in range(2):
                    code=original_call(self,function,b'automatic-rejected-input',count,indices)
                    assert comm.allgather(int(code))==[1117]*comm.size,code
                    if inspect_mpi:assert mpi_counts()[0]==0
                    admission_cases.append('automatic-layout-retry-'+str(attempt))
            if '--geometry-admission' in sys.argv:
                import re
                text=(self.directory/'inputfile').read_text()
                specifications=[
                    ('grid-overdetermined',text.replace('&maxwell',"&maxwell\n num_rgrid_em=40,40,40",1),1101),
                    ('boundary-conflict',text.replace("boundary_em(1,1) = 'abc'","boundary_em(1,1) = 'periodic'",1),1102),
                    ('time-overdetermined',text.replace('&maxwell',"&maxwell\n at_em=1",1),1103),
                    ('cfl',re.sub(r'\bdt_em\s*=\s*[-+0-9.eEdD]+','dt_em = 1.0',text),1104),
                    ('small-local-grid',re.sub(r'\bal_em\s*=.*',f'al_em = 80, {2*(4*(comm.size-1)+3)}, 80',text),1106),
                    ('maxwell-k-layout',re.sub(r'\bnproc_k\s*=\s*1\b','nproc_k = 2',text),1105),
                    ('maxwell-orbital-layout',re.sub(r'\bnproc_ob\s*=\s*1\b','nproc_ob = 2',text),1105),
                    ('mpi-count-mismatch',re.sub(r'\bnproc_rgrid\s*=.*',f'nproc_rgrid = 1, {comm.size+1}, 1',text),1105),
                    ('unimplemented-mapping',text.replace('&parallel',"&parallel\n process_allocation='user_mapping_extension'",1),1107),
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
            if '--shape-input-admission' in sys.argv:
                text=(self.directory/'inputfile').read_text()
                header="material IDs\nnot electronic density\n1 0 0 0\n40 2 0 0\n40 0 2 0\n40 0 0 2\n1 1 0 0 0\n"
                samples=[('missing',None),('header',b'truncated\n'),
                    ('dimensions',header.replace('40 2 0 0','39 2 0 0').encode()),
                    ('truncated',(header+'0 0 0 0 0 0\n').encode()),
                    ('invalid-number',(header+'bad 0 0 0 0 0\n').encode()),
                    ('nan',(header+'NaN 0 0 0 0 0\n').encode()),
                    ('infinity',(header+'Inf 0 0 0 0 0\n').encode()),
                    ('nan-grid-count',header.replace('40 2 0 0','NaN 2 0 0').encode()),
                    ('nan-tail',(header+('0 0 0 0 0 0\n'*(40**3//6))+'0 0 0 NaN\n').encode()),
                    ('infinity-tail',(header+('0 0 0 0 0 0\n'*(40**3//6))+'0 0 0 Inf\n').encode())]
                for name,payload in samples:
                    for attempt in range(2):
                        filename='missing-shape.cube' if payload is None else 'rejected-shape.cube'
                        if comm.rank==0:
                            (self.directory/'shape-rejected-input').write_text(shape_text(text,filename))
                            if payload is not None:(self.directory/filename).write_bytes(payload)
                        comm.Barrier()
                        code=original_call(self,function,b'shape-rejected-input',count,indices)
                        assert comm.allgather(int(code))==[1115]*comm.size,(name,code)
                        if inspect_mpi:assert mpi_counts()[0]==0
                        admission_cases.append('shape-'+name+'-'+str(attempt))
                if comm.size>1:
                    # Different file visibility/data on one owner must not strand peers.
                    from salmonpy._native_group import _directory_call
                    rank_dir=self.directory/f'shape-rank-{comm.rank}';rank_dir.mkdir()
                    (rank_dir/'inputfile').write_text(shape_text(text,'local-shape.cube'))
                    payload=(header+'0 0 0 0 0 0\n') if comm.rank==comm.size-1 else (header+('0 0 0 0 0 0\n'*(40**3//6))+'0 0 0 0\n')
                    (rank_dir/'local-shape.cube').write_text(payload)
                    for attempt in range(2):
                        code=_directory_call(self.lib,comm,rank_dir,function,b'inputfile',count,indices,minimum_protocol=3)
                        assert comm.allgather(int(code))==[1115]*comm.size
                        if inspect_mpi:assert mpi_counts()[0]==0
                        admission_cases.append('shape-rank-local-failure-'+str(attempt))
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
                               "&functional\n exx_mlwf_interval=12,\n"
                               " exx_mlwf_maxiter=14,\n"
                               " exx_mlwf_tolerance=2d-7,\n"
                               " exx_k_algorithm='user_route', exx_factor_history_frames=0,\n"
                               " exx_factor_n_teach=0, exx_factor_exact_interval=0\n/\n")
                    (self.directory/'inputfile').write_text(
                        text.replace('&maxwell',"&maxwell\n yn_wf_em='N'",1)
                        +"\n&propagation\n propagator='user_extension'\n/\n"
                        +"&scf\n method_init_wf='user_initializer'\n/\n")
                comm.Barrier()


        if getattr(function,'__name__','')=='sp_eh_init' and host_unit is not None and '--host-unit-lifetime' not in sys.argv:
            assert host_unit(b'',1)==0
            assert host_unit(b'',2)==0
            assert (out/f'host-{comm.rank}.data').read_text()=='host-before\nhost-after\n'
            host_unit=None
            admission_cases.append('host-unit-preserved-through-rejected-inputs')
        return original_call(self,function,*args)
    if '--material-input-admission' in sys.argv:
        if material_mode is None:raise ValueError('choose material-input-mode for material admission')
        from salmonpy.errors import InputError
        from salmonpy import embedded_tddft
        asset={'diagonal':'salmonpy-ld-diagonal.dat','tensor':'salmonpy-ld-tensor.dat',
               'subpixel':'salmonpy-subpixel.dat'}[material_mode]
        expected=1113 if material_mode=='subpixel' else 1110
        def corrupt(text,variant):
            if variant=='missing':return None
            if variant=='header':return 'INVALID_MATERIAL_ASSET\n'
            if variant=='truncated':return text.splitlines()[0]+'\n'
            lines=text.splitlines();row={'diagonal':3,'tensor':4,'subpixel':1}[material_mode]
            values=lines[row].split();values[0]='NaN';lines[row]=' '.join(values)
            return '\n'.join(lines)+'\n'
        for variant in ('missing','header','truncated','nonfinite'):
            def rejected_material_call(self,function,*call_args):
                if getattr(function,'__name__','')=='sp_eh_init':
                    if comm.rank==0:
                        path=self.directory/asset;value=corrupt(path.read_text(),variant)
                        if value is None:path.unlink()
                        else:path.write_text(value)
                    comm.Barrier()  # Fixture mutates assets after public staging.
                return original_call(self,function,*call_args)
            for attempt in range(2):
                try:
                    with patch.object(NativeEHField,'_call',rejected_material_call):
                        NativeEHField(ctx,calc,build,out/f'rejected-material-{variant}-{attempt}',
                                      points=sample_points,grid=(1,comm.size,1))
                except InputError as exc:
                    assert str(expected) in str(exc),exc
                else:raise AssertionError('invalid material input was accepted')
                assert ctx._dependents==0
                if inspect_mpi:assert mpi_counts()[:2]==[0,0]
                assert not embedded_tddft._used,'clean native material rejection still consumes Python lifetime'
                comm.Barrier()  # Caller survives; correct input will be staged next.
                admission_cases.append('material-'+material_mode+'-'+variant+'-'+str(attempt))
    if '--coefficient-admission' in sys.argv:
        if material_mode is None:raise ValueError('choose material-input-mode for coefficient admission')
        import re
        from salmonpy.errors import InputError
        from salmonpy import embedded_tddft
        expected={'diagonal':1111,'tensor':1112,'subpixel':1114}[material_mode]
        def rejected_coefficient_call(self,function,*call_args):
            if getattr(function,'__name__','')=='sp_eh_init':
                if comm.rank==0:
                    path=self.directory/'inputfile';text=path.read_text()
                    if material_mode=='diagonal':
                        text=re.sub(r'omega_p_ld\(1\)\s*=.*','omega_p_ld(1) = 1e200',text)
                    elif material_mode=='tensor':
                        text=text.replace("yn_periodic = 'n'","yn_periodic = 'y'",1)
                        text=text.replace('&maxwell',"&maxwell\n bloch_k_em=0,.01,0",1)
                        text=re.sub(r'source_loc1\s*=.*','source_loc1 = 16,0,0',text)
                        text=re.sub(r'ori_s\(1,:\)\s*=.*','ori_s(1,:) = 32,40,40',text)
                        text=re.sub(r'ori_s\(2,:\)\s*=.*','ori_s(2,:) = 36,40,40',text)
                        for axis in (2,3):
                            for face in (1,2):
                                text=text.replace(f"boundary_em({axis},{face}) = 'abc'",f"boundary_em({axis},{face}) = 'periodic'")
                    else:
                        asset=self.directory/'salmonpy-subpixel.dat';lines=asset.read_text().splitlines()
                        # First material oscillator rate: finite input, overflow in native recurrence.
                        lines[7]='1e200, 2';asset.write_text('\n'.join(lines)+'\n')
                    path.write_text(text)
                comm.Barrier()  # Validation-only input mutation after public staging.
            return original_call(self,function,*call_args)
        for attempt in range(2):
            try:
                with patch.object(NativeEHField,'_call',rejected_coefficient_call):
                    NativeEHField(ctx,calc,build,out/f'rejected-coefficient-{attempt}',
                                  points=sample_points,grid=(1,comm.size,1))
            except InputError as exc:
                assert str(expected) in str(exc),exc
            else:raise AssertionError('invalid coefficient/mode initialization was accepted')
            assert ctx._dependents==0
            if inspect_mpi:assert mpi_counts()[:2]==[0,0]
            assert not embedded_tddft._used,'clean coefficient rejection still consumes Python lifetime'
            comm.Barrier()
            admission_cases.append('coefficient-'+material_mode+'-'+str(attempt))
    if '--allocation-admission' in sys.argv or '--ld-allocation-admission' in sys.argv or '--aux-allocation-admission' in sys.argv or '--material-allocation-failures' in sys.argv or '--ld-coefficient-allocation-failures' in sys.argv:
        # Deterministic STAT failure, confined to a separately built test runtime.
        from salmonpy.errors import InputError
        from salmonpy import embedded_tddft
        allocation_stages=[]
        if '--allocation-admission' in sys.argv:allocation_stages.extend((24,25))
        if '--ld-allocation-admission' in sys.argv:allocation_stages.extend((26,27))
        if '--aux-allocation-admission' in sys.argv:allocation_stages.extend((28,29))
        if '--material-allocation-failures' in sys.argv:
            if not tensor and material_mode!='diagonal':raise ValueError('material allocation tests require tensor or diagonal LD')
            allocation_stages.extend((30,) if material_mode=='diagonal' else (30,31,32))
        if '--ld-coefficient-allocation-failures' in sys.argv:
            allocation_stages.extend((33,34,35,36) if tensor else (33,34,35))
        for stage in allocation_stages:
            def allocation_init_call(self,function,*call_args):
                if getattr(function,'__name__','')=='sp_eh_init':
                    hook=self.lib.sp_eh_test_fault
                    hook.argtypes=[c.c_int,c.c_int];hook.restype=c.c_int
                    assert hook(stage,comm.size-1)==0,'allocation fault stage unavailable'
                return original_call(self,function,*call_args)
            for attempt in range(2):
                try:
                    with patch.object(NativeEHField,'_call',allocation_init_call):
                        NativeEHField(ctx,calc,build,out/f'allocation-{stage}-{attempt}',
                                      points=sample_points,grid=(1,comm.size,1))
                except InputError as exc:
                    assert str(1120 if stage in (28,29,30,31,32,33,34,35,36) else 1119 if stage in (26,27) else 1118) in str(exc),exc
                else:raise AssertionError('rank-local allocation failure was not returned')
                assert ctx._dependents==0
                if inspect_mpi:assert mpi_counts()[:2]==[0,0]
                assert not embedded_tddft._used,'allocation rejection consumes Python lifetime'
                comm.Barrier()
                admission_cases.append(f'allocation-{stage}-retry-{attempt}')
    if '--grid-init-retry' in sys.argv:
        # Check public lifetime release, beyond the raw C entry alone.
        import re
        from salmonpy.errors import InputError
        retry_inputs=[
            ('grid-specification',1101,lambda text:text.replace('&maxwell',"&maxwell\n num_rgrid_em=40,40,40",1)),
            ('boundary-pairing',1102,lambda text:text.replace("boundary_em(1,1) = 'abc'","boundary_em(1,1) = 'periodic'",1)),
            ('time-specification',1103,lambda text:text.replace('&maxwell',"&maxwell\n at_em=1",1)),
            ('cfl',1104,lambda text:re.sub(r'\bdt_em\s*=\s*[-+0-9.eEdD]+','dt_em = 1.0',text)),
            ('grid',1106,lambda text:re.sub(r'\bal_em\s*=.*',f'al_em = 80, {2*(4*(comm.size-1)+3)}, 80',text)),
            ('partition',1105,lambda text:re.sub(r'\bnproc_rgrid\s*=.*',f'nproc_rgrid = 1, {comm.size+1}, 1',text)),
            ('mapping',1107,lambda text:text.replace('&parallel',"&parallel\n process_allocation='user_mapping_extension'",1)),
        ]
        if '--auto-layout-admission' in sys.argv:
            retry_inputs.append(('automatic-layout',1117,automatic_input))
        if '--shape-input-admission' in sys.argv:
            retry_inputs.append(('shape',1115,lambda text:shape_text(text,'missing-shape.cube')))
        for name,expected,transform in retry_inputs:
            def rejected_init_call(self,function,*call_args):
                if getattr(function,'__name__','')=='sp_eh_init':
                    if comm.rank==0:
                        path=self.directory/'inputfile';path.write_text(transform(path.read_text()))
                    comm.Barrier()
                return original_call(self,function,*call_args)
            for attempt in range(2):
                try:
                    with patch.object(NativeEHField,'_call',rejected_init_call):
                        NativeEHField(ctx,calc,build,out/f'rejected-{name}-{attempt}',points=sample_points,grid=(1,comm.size,1))
                except InputError as exc:
                    assert str(expected) in str(exc),exc
                else:raise AssertionError('invalid MPI setup was not rejected: '+name)
                assert ctx._dependents==0
                if inspect_mpi:assert mpi_counts()[:2]==[0,0]
                from salmonpy import embedded_tddft
                assert not embedded_tddft._used,'clean early initialization rejection still consumes Python lifetime'
                comm.Barrier()  # Caller MPI survives rejected constructor/detach.
                admission_cases.append('public-'+name+'-retry-'+str(attempt))
    with patch.object(NativeEHField,'_call',admission_call) if '--auto-layout-admission' in sys.argv or '--shape-input-admission' in sys.argv or '--geometry-admission' in sys.argv or '--admission' in sys.argv or '--parsed-admission' in sys.argv or '--malformed-admission' in sys.argv or '--unit-admission' in sys.argv or '--coordinate-admission' in sys.argv or '--semantic-admission' in sys.argv else patch.object(NativeEHField,'_call',original_call):
        f=NativeEHField(ctx,calc,build,out/'field',points=sample_points,grid=(1,comm.size,1))
    if host_unit is not None and '--host-unit-lifetime' not in sys.argv:
        assert host_unit(b'',1)==0
        assert host_unit(b'',2)==0
        assert (out/f'host-{comm.rank}.data').read_text()=='host-before\nhost-after\n'
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
    def files():return {p.name:p.read_bytes() for p in f.directory.glob('*.data') if p.is_file()}
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
        if '--checkpoint-save-failures' in sys.argv and step==8:
            from salmonpy._native_checkpoint import NativeRestoreCandidate
            reference=save('checkpoint-save-before')
            saved_outputs=files()
            original_A=f.A.copy();original_E=f.E.copy();original_step=f.step_index
            verification=[]
            def verify_save_state(label):
                import tempfile,hashlib
                with tempfile.TemporaryDirectory(prefix='checkpoint-verify-',dir=out) as scratch:
                    path=Path(scratch)/f'rank-{comm.rank}.bin'
                    assert f._call(cp,str(path).encode(),0)==0
                    data=path.read_bytes()
                    assert data==reference
                    verification.append({'check':label,'sha256':hashlib.sha256(data).hexdigest()})
            for kind in ('existing','directory','missing-parent'):
                target=out/f'save-{kind}-{comm.rank}.bin'
                if comm.rank==comm.size-1:
                    if kind=='existing':target.write_bytes(b'previous checkpoint must remain unchanged')
                    elif kind=='directory':target.mkdir()
                    else:target=out/f'absent-parent-{comm.rank}'/target.name
                code=f._call(cp,str(target).encode(),0)
                assert code!=0,'rank-local save failure was not returned collectively'
                if comm.rank==comm.size-1 and kind=='existing':
                    assert target.read_bytes()==b'previous checkpoint must remain unchanged'
                verify_save_state(kind+'-after')
                assert files()==saved_outputs and not f.failed and f.step_index==original_step
                np.testing.assert_array_equal(f.A,original_A);np.testing.assert_array_equal(f.E,original_E)
                # Preserve partial failed-save files as evidence; they are not an archive.
                if target.exists():target.rename(out/f'failed-save-{kind}-{comm.rank}')
                if not target.parent.exists():target.parent.mkdir()
                assert f._call(cp,str(target).encode(),0)==0
                assert target.read_bytes()==reference
                candidate=NativeRestoreCandidate(f,'sp_eh_checkpoint',target)
                with patch.object(f,'advance',side_effect=AssertionError('save/retry/restore must not advance')):
                    candidate.stage();candidate.publish()
                verify_save_state(kind+'-restored')
                assert files()==saved_outputs and f.step_index==original_step
                comm.Barrier()
                admission_cases.append('checkpoint-save-failure-repair-'+kind)
            (out/f'checkpoint-save-verification-{comm.rank}.json').write_text(json.dumps(verification,indent=2)+'\n')
        if snapshot_failure_mode and step==18:
            if comm.rank==0:
                path=f.directory/'obs1_ex_xy_18.data'
                if snapshot_failure_mode=='open':path.mkdir()
                else:
                    path.write_bytes(b'');path.chmod(0o400)
            comm.Barrier()
        if point_failure_mode and point_failure_mode.startswith('append') and step==4:
            if comm.rank==0:
                path=f.directory/point_output_name
                if point_failure_mode=='append-open':
                    path.rename(out/'point-before-fault.data');path.mkdir()
                else:path.chmod(0o400)
            comm.Barrier()
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
        if point_failure_mode and step==6:
            if comm.rank==0:
                path=f.directory/point_output_name
                if path.is_dir():path.rmdir()
                elif path.exists():path.chmod(0o600)
                if point_failure_mode.startswith('header'):
                    path.write_bytes(b''.join(point_reference.read_bytes().splitlines(keepends=True)[:output_header_lines]))
                elif point_failure_mode=='append-open':(out/'point-before-fault.data').rename(path)
            comm.Barrier()
        if resume_midpoint and step==16:
            save('midpoint')
            midpoint=dict(A=f.A.copy(),E=f.E.copy(),step_index=f.step_index)
        if resume_midpoint and step==18:
            if snapshot_failure_mode:
                if comm.rank==0:
                    path=f.directory/'obs1_ex_xy_18.data'
                    if path.is_dir():path.rmdir()
                    else:path.chmod(0o600)
                comm.Barrier()
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
    if info_failure_mode:
        if comm.rank==0:
            path=f.directory/'obs0_info.data'
            if info_failure_mode=='open':path.mkdir()
            else:
                path.write_bytes(b'');path.chmod(0o400)
        comm.Barrier()
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
    if host_unit is not None and '--host-unit-lifetime' in sys.argv:
        assert host_unit(b'',1)==0, 'host file unit was closed by native output'
        assert host_unit(b'',2)==0
        assert (out/f'host-{comm.rank}.data').read_text()=='host-before\nhost-after\n'
        admission_cases.append('host-unit-preserved-through-init-step-restore-finalize')
    assert Path.cwd()==anchor
    if guard:guard.stop()
    if comm.rank==0:
        (out/'validation.json').write_text(json.dumps(dict(status='passed',compact_checkpoint_fixture=compact,snapshot_output_failure=snapshot_failure_mode,info_output_failure=info_failure_mode,point_output_failure=point_failure_mode,observation_output_kind='integral' if integral_failure else 'point',process_allocation='orbital_sequential' if '--orbital-mapping' in sys.argv else 'grid_sequential',reciprocal_tensor=tensor and material_mode!='diagonal',material_input_mode=material_mode,epsilon_infinity=[m.config.epsilon_infinity for m in ld],checkpoint_restore_without_reinitialization=resume_midpoint,validation_continuation_repeated_steps=2 if resume_midpoint else 0,new_runtime_no_chdir=no_chdir,admission_cases=admission_cases,geometry_mpi=geometry_mpi,cases=cases,ranks=comm.size,steps=stop,capacity=steps,replayed_steps=0,allocation_before_close=allocation_before,allocation_after_close=allocation_after,mpi_before=mpi_before,mpi_initialized=mpi_initialized,mpi_after_native_close=mpi_after_native,mpi_after_detach=mpi_after_detach,other_owner_preserved=bool(host_mpi)),indent=2)+'\n')
        print('E/H deep recovery PASS',cases)


if __name__=='__main__':main()
