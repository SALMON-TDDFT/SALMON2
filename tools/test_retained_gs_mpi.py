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


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--size',type=int,choices=(1,2),default=2)
    parser.add_argument('--split',choices=('orbital','spatial'),default='orbital')
    parser.add_argument('--xc-a',choices=('pz','libxc_pz','libxc_pbe'),default='pz')
    parser.add_argument('--xc-b',choices=('pz','libxc_pz','libxc_pbe'),default='libxc_pbe')
    parser.add_argument('--gamma',action='store_true')
    parser.add_argument('--tolerance',type=float,default=1e-12)
    parser.add_argument('--b-iterations',type=int,default=4)
    parser.add_argument('--b-tolerance',type=float)
    args=parser.parse_args();world=MPI.COMM_WORLD
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
    backend=NativeGSBackend(args.build,work,decomposition=dec)
    scf=SCFSettings(max_iterations=4,tolerance=args.tolerance)
    scf_b=SCFSettings(max_iterations=args.b_iterations,tolerance=args.tolerance if args.b_tolerance is None else args.b_tolerance)
    frames={};evidence={};warnings.simplefilter('always')

    def snapshot(state):
        result={'info':state.info}
        for q in ('density','orbitals','occupations','eigenvalues','forces'):
            arr=state.array(q)
            assert not arr.values.flags.writeable
            try:arr.values.setflags(write=True)
            except ValueError:pass
            else:raise AssertionError('array became writable')
            result[q]=(arr.lower_indices,arr.values)
        return result

    def equal(first,second):
        assert first['info']==second['info']
        for q in ('density','orbitals','occupations','eigenvalues','forces'):
            assert first[q][0]==second[q][0]
            np.testing.assert_array_equal(first[q][1],second[q][1])

    with backend.open(context) as session:
        if g==0:
            first=session.run(a,scf=scf,material_id='silicon-a')
            before=snapshot(first)
            second=session.run(b,scf=scf_b,material_id='silicon-b')
            equal(before,snapshot(first));frames['a']=before;frames['b']=snapshot(second)
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
            second.close();assert session.resources()['live_child_communicators']==0
        else:
            key='a' if g==1 else 'b'
            state=session.run(a if g==1 else b,scf=scf if g==1 else scf_b,material_id='silicon-'+key)
            frames[key]=snapshot(state)
    assert context._dependents==0
    gathered=world.allgather(frames)
    for offset in range(args.size):
        equal(gathered[offset]['a'],gathered[args.size+offset]['a'])
        equal(gathered[offset]['b'],gathered[2*args.size+offset]['b'])
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
            'native_convergence_flags':{key:frame['info'].converged for key,frame in frames.items()},'GS_restart_collected':False},indent=2)+'\n')
    comm.Free()


if __name__=='__main__':
    try:main()
    except BaseException:
        import traceback;traceback.print_exc();MPI.COMM_WORLD.Abort(1)
