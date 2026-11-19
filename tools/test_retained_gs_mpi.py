"""Three groups: retained A+B vs independent A/B, no GS restart collection.

mpiexec -n 6 python tools/test_retained_gs_mpi.py BUILD NEW_OUTPUT --size 2 --split spatial
Only digest/metadata evidence is written; requested numerical arrays stay in memory.
"""
import argparse
import ctypes as c
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import warnings
import numpy as np
from mpi4py import MPI
from salmonpy.dft import DFT,Pseudopotential,SCFSettings
from salmonpy.native_gs import NativeGSBackend
from salmonpy.errors import InputError,RuntimeFailure


def import_controls(args):
    """Three-group import comparison; no SCF or saved numerical-array dataset."""
    w=MPI.COMM_WORLD
    if w.size!=3*args.size:raise ValueError('launch 3*size ranks')
    if w.rank==0:args.output.mkdir(parents=True,exist_ok=False)
    w.Barrier();group=w.rank//args.size;comm=w.Split(group,w.rank)
    work=args.output/str(group)
    if comm.rank==0:work.mkdir()
    comm.Barrier();ctx=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    backend=NativeGSBackend(args.build,work,decomposition=(args.size,1,1,1,1))
    def snapshot(state):
        info=state.info
        assert info.energy_hartree is None and info.residual is None and info.converged is None
        result={}
        for q in ('density','occupations','orbitals'):
            a=state.array(q);assert not a.values.flags.writeable
            result[q]=(a.lower_indices,a.values.copy())
        return result
    def equal(a,b):
        for q in a:
            assert a[q][0]==b[q][0];np.testing.assert_array_equal(a[q][1],b[q][1])
    bundles={1:args.import_bundle,2:args.second_bundle or args.import_bundle}
    frames={}
    with backend.open(ctx) as session:
        if group:
            state=session.load(bundles[group],material_id='material-'+str(group))
            frames[group]=snapshot(state);state.close()
        else:
            a=session.load(bundles[1],material_id='material-1');before=snapshot(a)
            b=session.load(bundles[2],material_id='material-2');equal(before,snapshot(a))
            frames={1:before,2:snapshot(b)}
            for operation in (lambda:session.run(None),
                              lambda:a.export_archive(work/'must-not-exist')):
                try:operation()
                except InputError:pass
                else:raise AssertionError('computed-only operation was accepted')
            equal(before,snapshot(a));equal(frames[2],snapshot(b))
            for q in ('forces','eigenvalues'):
                try:a.array(q)
                except InputError:pass
                else:raise AssertionError('uncomputed quantity advertised')
            old=a._token;a.close();equal(frames[2],snapshot(b))
            fresh=session.load(bundles[1],material_id='material-1');equal(before,snapshot(fresh))
            assert fresh._token!=old
            assert session.lib.sp_gs_close(old)!=0
            if comm.size>1:
                bad=fresh._token if comm.rank==0 else b._token
                stats=np.zeros(4,np.int32);values=np.zeros(3)
                assert session.lib.sp_gs_info(bad,stats.ctypes.data,values.ctypes.data)!=0
                equal(before,snapshot(fresh));equal(frames[2],snapshot(b))
            fresh.close();b.close()
        assert session.resources()['states']==0
        assert session.resources()['live_child_communicators']==0
    local=comm.gather(frames,root=0)
    all_frames=w.gather(local if comm.rank==0 else None,root=0)
    if w.rank==0:
        mixed=all_frames[0];first=all_frames[args.size];second=all_frames[2*args.size]
        for rank in range(args.size):
            equal(mixed[rank][1],first[rank][1]);equal(mixed[rank][2],second[rank][2])
        result={'status':'pass','ranks':w.size,'independent_arrays_exact':True,
                'fresh_import_exact':True,'close_isolation':True,'stale_tokens_rejected':True,
                'collective_mixed_tokens_rejected':args.size>1,'SCF_calls':0}
        (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--size',type=int,choices=(1,2),default=2)
    parser.add_argument('--current-steps',type=int,default=0,help='also exercise computed GS to independent RT states and direct resume')
    parser.add_argument('--baseline',type=Path,help='explicit independently built control runtime')
    parser.add_argument('--check-export-assembly',action='store_true',help='intercept shared writer and inspect assembled arrays in RAM; do not write GS archives')
    parser.add_argument('--split',choices=('orbital','spatial'),default='orbital')
    parser.add_argument('--xc-a',choices=('pz','libxc_pz','libxc_pbe'),default='pz')
    parser.add_argument('--xc-b',choices=('pz','libxc_pz','libxc_pbe'),default='libxc_pbe')
    parser.add_argument('--gamma',action='store_true')
    parser.add_argument('--tolerance',type=float,default=1e-12)
    parser.add_argument('--b-iterations',type=int,default=4)
    parser.add_argument('--b-tolerance',type=float)
    parser.add_argument('--import-bundle',type=Path,help='existing synthetic GS; import-only mode, never SCF')
    parser.add_argument('--second-bundle',type=Path)
    args=parser.parse_args()
    if args.current_steps and args.current_steps<4:parser.error('--current-steps must be at least four')
    if args.import_bundle:
        import_controls(args);return
    world=MPI.COMM_WORLD
    if world.size!=3*args.size:raise ValueError('launch exactly 3*size ranks')
    if world.rank==0:args.output.mkdir(parents=True,exist_ok=False)
    world.Barrier();g=world.rank//args.size;comm=world.Split(g,world.rank)
    work=args.output/str(g)
    if comm.rank==0:work.mkdir()
    comm.Barrier();context=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    root=Path(__file__).resolve().parents[1]
    pseudo=root/'salmonpy/examples/si_intro_data/Si_rps.dat'
    if not pseudo.exists():pseudo=root.parent/'examples/si_intro_data/Si_rps.dat'
    pp=Pseudopotential.from_file(pseudo,atomic_number=14,local_channel=2)
    a=DFT((14,),np.eye(3)*8,[[.3,.4,.45]],(pp,),4,2,(8,8,8),
          (1,1,1) if args.gamma else (2,1,1),xc=args.xc_a)
    b=DFT((14,14),np.eye(3)*9,[[.3,.4,.45],[.7,.6,.55]],(pp,),8,4,(12,12,12),(1,2,1),xc=args.xc_b)
    dec=(1,args.size,1,1,1) if args.split=='orbital' else (1,1,args.size,1,1)
    backend=NativeGSBackend(args.baseline if g and args.baseline else args.build,work,decomposition=dec)
    scf=SCFSettings(max_iterations=4,tolerance=args.tolerance)
    scf_b=SCFSettings(max_iterations=args.b_iterations,tolerance=args.tolerance if args.b_tolerance is None else args.b_tolerance)
    frames={};evidence={};responses={};warnings.simplefilter('always')

    def snapshot(state):
        result={'info':state.info}
        for q in ('density','orbitals','occupations','eigenvalues','forces'):
            arr=state.array(q)
            assert not arr.values.flags.writeable
            try:arr.values.setflags(write=True)
            except ValueError:pass
            else:raise AssertionError('array became writable')
            result[q]=(arr.lower_indices,arr.values)
        if args.check_export_assembly:
            import salmonpy._gs_archive as archive
            original=archive._write;seen=[]
            def inspect(gs,dest,block,waves,occ,eigen,info,model):
                assert waves.shape==(int(np.prod(model.kmesh)),model.states,*reversed(model.grid))
                assert np.isfinite(waves).all() and np.isfinite(occ).all() and np.isfinite(eigen).all()
                seen.append({'shape':list(waves.shape),'wave_sha256':hashlib.sha256(waves.tobytes()).hexdigest()})
                return None
            archive._write=inspect
            try:state.export_archive(work/('unwritten-'+str(state._token)))
            finally:archive._write=original
            assert not (work/('unwritten-'+str(state._token))).exists()
            if comm.rank==0:assert len(seen)==1
        return result

    def equal(first,second):
        assert first['info']==second['info']
        for q in ('density','orbitals','occupations','eigenvalues','forces'):
            assert first[q][0]==second[q][0]
            np.testing.assert_array_equal(first[q][1],second[q][1])

    def begin_response(gs,key):
        if not args.current_steps:return None
        model=gs.tddft_model(dt=.001,max_steps=args.current_steps+1)
        assert model.material_id==gs.material_id
        first=model.new_state('first');sibling=model.new_state('sibling')
        def drive(i):return 1e-4*np.array([np.sin(.31*i),np.cos(.23*i)-1,np.sin(.17*i)])
        trace=[];half=args.current_steps//2
        for i in range(1,half+1):
            trace.append(np.r_[first.step(drive(i)),sibling.step(-drive(i))])
        path=work/f'rt-{key}-{comm.rank}.bin'
        first.checkpoint(path,False)
        restored=model.new_state('restored');restored.checkpoint(path,True)
        path.unlink()  # transient RT restart control; never retain a GS dataset
        assert restored.steps==first.steps==half
        return model,first,sibling,restored,trace,drive,half

    def finish_response(rt,key):
        if rt is None:return
        model,first,sibling,restored,trace,drive,half=rt
        try:
            # Seed remains owned after the original computed GS has closed.
            fresh=model.new_state('fresh');fresh.close()
            for i in range(half+1,args.current_steps+1):
                trace.append(np.r_[first.step(drive(i)),sibling.step(-drive(i))])
                np.testing.assert_array_equal(restored.step(drive(i)),first.current)
                assert restored.steps==first.steps==i
            responses[key]=np.asarray(trace)
            assert not list(model.work.rglob('wfn.bin'))
        finally:
            restored.close();first.close();sibling.close();model.close()

    with backend.open(context) as session:
        if g==0:
            first=session.run(a,scf=scf,material_id='silicon-a')
            before=snapshot(first)
            rt_a=begin_response(first,'a')
            second=session.run(b,scf=scf_b,material_id='silicon-b')
            equal(before,snapshot(first));frames['a']=before;frames['b']=snapshot(second)
            rt_b=begin_response(second,'b')
            assert first.material_id=='silicon-a' and second.material_id=='silicon-b'
            counts=session.resources();assert counts['states']==2 and counts['live_child_communicators']>0
            # A failed collective input agreement must not change either owner.
            if args.size>1:
                called=[]
                def must_not_enter():called.append(True)
                missing=work/'rank-local-missing' if comm.rank==args.size-1 else work
                try:session._call(missing,must_not_enter)
                except InputError:pass
                else:raise AssertionError('rank-local missing work accepted')
                assert not called
                equal(before,snapshot(first));equal(frames['b'],snapshot(second))
                try:session.run(a,scf=scf,material_id='same' if comm.rank==0 else 'different')
                except InputError:pass
                else:raise AssertionError('different material identities accepted')
                equal(before,snapshot(first));equal(frames['b'],snapshot(second))
                stats=np.zeros(4,dtype=np.int32);vals=np.zeros(3)
                bad=first._token if comm.rank==0 else second._token
                code=session.lib.sp_gs_info(bad,stats.ctypes.data,vals.ctypes.data)
                assert code!=0
            old=first._token;first.close();first.close()
            assert session.resources()['states']==1
            equal(frames['b'],snapshot(second))
            remaining_children=session.resources()['live_child_communicators']
            # Freed slots must not authorize an old token; live B persists while A runs again.
            for repeat in range(3):
                fresh=session.run(a,scf=scf,material_id='silicon-a')
                equal(frames['a'],snapshot(fresh));equal(frames['b'],snapshot(second))
                stats=np.zeros(4,dtype=np.int32);vals=np.zeros(3)
                assert session.lib.sp_gs_info(old,stats.ctypes.data,vals.ctypes.data)!=0
                fresh.close()
                now=session.resources();assert now['states']==1
                assert now['live_child_communicators']==remaining_children
                assert now['communicator_slots']<=counts['communicator_slots']
            second.close()
            finish_response(rt_a,'a');finish_response(rt_b,'b')
            assert session.resources()['live_child_communicators']==0
        else:
            key='a' if g==1 else 'b'
            state=session.run(a if g==1 else b,scf=scf if g==1 else scf_b,material_id='silicon-'+key)
            frames[key]=snapshot(state)
            response=begin_response(state,key);state.close();finish_response(response,key)
    assert context._dependents==0
    gathered=world.allgather(frames)
    rt_frames=world.allgather(responses)
    for offset in range(args.size):
        for key,index in (('a',args.size+offset),('b',2*args.size+offset)):
            first,second=gathered[offset][key],gathered[index][key]
            if args.current_steps:np.testing.assert_array_equal(rt_frames[offset][key],rt_frames[index][key])
            if args.baseline:
                for name in ('iterations','converged','atoms','bands'):assert getattr(first['info'],name)==getattr(second['info'],name)
                np.testing.assert_allclose([first['info'].energy_hartree,first['info'].residual,first['info'].tolerance],
                                          [second['info'].energy_hartree,second['info'].residual,second['info'].tolerance],rtol=3e-12,atol=3e-13)
                for q in ('density','orbitals','occupations','eigenvalues','forces'):
                    assert first[q][0]==second[q][0]
                    np.testing.assert_allclose(first[q][1],second[q][1],rtol=3e-12,atol=3e-13)
            else:equal(first,second)
    for key,frame in frames.items():
        evidence[key]={'info':vars(frame['info']),'arrays':{q:{'lower_indices':frame[q][0],
            'shape':frame[q][1].shape,'sha256':hashlib.sha256(frame[q][1].tobytes()).hexdigest()}
            for q in ('density','orbitals','occupations','eigenvalues','forces')}}
    (work/f'rank-{comm.rank}.json').write_text(json.dumps(evidence,indent=2)+'\n')
    world.Barrier()
    if world.rank==0:
        (args.output/'result.json').write_text(json.dumps({'status':'pass','size':args.size,'split':args.split,
            'xc_a':args.xc_a,'xc_b':args.xc_b,'gamma_a':args.gamma,'b_iterations':args.b_iterations,'b_tolerance':args.b_tolerance,
            'independent_arrays_exact':True,
            'retained_siblings_unchanged':True,'close_reuse_stale_tokens':True,'MPI_resources_released':True,
            'SCF_nonconvergence_warns_and_retains':args.tolerance==1e-12,
            'native_convergence_flags':{key:frame['info'].converged for key,frame in frames.items()},'GS_restart_collected':False,'computed_GS_to_RT_steps':args.current_steps,
            'RT_independent_and_direct_resume_exact':bool(args.current_steps),
            'RT_replayed_steps':0,'RT_after_source_GS_closed':bool(args.current_steps)},indent=2)+'\n')
    comm.Free()


if __name__=='__main__':
    try:main()
    except BaseException:
        import traceback;traceback.print_exc();MPI.COMM_WORLD.Abort(1)
