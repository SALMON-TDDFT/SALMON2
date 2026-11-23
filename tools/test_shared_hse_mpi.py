"""Real Python shared-model regression: BUILD FULL_K_FIXTURE NEW_OUTPUT.

Run with 2*--group-size MPI ranks. --mode single-a/single-b supplies independent
controls for mixed. Explicit runtime/inputs; no build/download. Synthetic test.
"""
import argparse
import ctypes as c
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.fast_hse import FastHSEMaterialProvider
from salmonpy.fast_hse_mpi import NativeFastHSEBackend
from salmonpy.native_steps import native_step_targets
from salmonpy.fast_hse_mpi import _agree
from salmonpy.errors import InputError,RuntimeFailure


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('build','fixture','output'):parser.add_argument(key,type=Path)
    parser.add_argument('--group-size',type=int,default=2)
    parser.add_argument('--mode',choices=('mixed','single-a','single-b'),default='mixed')
    parser.add_argument('--legacy',action='store_true')
    args=parser.parse_args();world=MPI.COMM_WORLD
    if args.legacy and args.mode=='mixed':parser.error('legacy is isolated reference only')
    assert world.size==2*args.group_size
    group=world.rank//args.group_size;comm=world.Split(group,world.rank)
    out=args.output.resolve()/f'group-{group}'
    if comm.rank==0:out.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    points=(args.mode[-1],) if args.legacy else ('a','b')
    context=SimpleNamespace(solver_comm=comm,point_ids=points,closed=False,_dependents=0)
    initial=FastHSEMaterialProvider(args.fixture,max_steps=34).prepare(context)
    backend=NativeFastHSEBackend(args.build,out)
    model=None if args.legacy else backend.open_model(initial,.001,context)
    if model is not None:assert context._dependents==1
    states={};rows={'a':[],'b':[]}
    def make(owner):
        if args.legacy:
            return backend.create_state(initial,owner,.001 if owner=='a' else .002,context)
        return model.new_state(owner,dt=.001 if owner=='a' else .002,capacity=32 if owner=='a' else 34)
    def value(owner,step):
        a=np.zeros(3)
        if owner=='a':a[2]=.0001*np.sin(.3*step)
        else:a[:2]=.0001*np.array([np.cos(.2*step)-1,np.sin(.2*step)])
        return a
    def advance(state,owner,step):
        j=state.step(value(owner,step));assert state.steps==step
        return np.r_[j,state.energy]
    def backup(state):
        response=SimpleNamespace(context=context,states={state.point_id:state},agree=lambda error:_agree(comm,error))
        return native_step_targets(SimpleNamespace(response=response,field=None))[0]
    def path(owner):return out/f'{owner}-saved-{comm.rank}.bin'
    if args.mode!='single-b':states['a']=make('a')
    if args.mode!='single-a':states['b']=make('b')
    if args.mode=='mixed':
        assert states['a'].work==states['b'].work==model.work
        assert states['a'].handle!=states['b'].handle
        assert set(states['a'].output_files).isdisjoint(states['b'].output_files)
    for step in range(1,33):
        if args.mode=='mixed' and step in (6,10):
            guards={owner:backup(state) for owner,state in states.items()}
            tokens={owner:guard.capture() for owner,guard in guards.items()}
        for owner,state in states.items():rows[owner].append(advance(state,owner,step))
        if args.mode=='mixed' and step==6:
            # Release B to flush its accepted step6; A rollback must not touch B.
            guards['b'].release(tokens['b'])
            before={n:(model.work/n).read_bytes() for n in states['b'].output_files}
            b_current=states['b'].current.copy()
            guards['a'].rollback(tokens['a'])
            np.testing.assert_array_equal(states['b'].current,b_current)
            assert {n:(model.work/n).read_bytes() for n in states['b'].output_files}==before
            np.testing.assert_array_equal(advance(states['a'],'a',step),rows['a'][-1])
            guards['a'].release(tokens['a'])
        if args.mode=='mixed' and step==10:
            for owner in states:
                guards[owner].rollback(tokens[owner])
                np.testing.assert_array_equal(advance(states[owner],owner,step),rows[owner][-1])
                guards[owner].release(tokens[owner])
        if step==8:
            for owner,state in states.items():state.checkpoint(path(owner),False)
    if 'b' in states:
        for step in (33,34):rows['b'].append(advance(states['b'],'b',step))
    if args.mode=='mixed':
        original={owner:(state.A.copy(),state.current.copy(),state.steps) for owner,state in states.items()}
        bad=make('a')
        bad_path=out/f'broken-{comm.rank}.bin'
        bad_path.write_bytes(path('a').read_bytes()[:256])
        with np.testing.assert_raises(RuntimeFailure):bad.checkpoint(bad_path,True)
        bad.close()
        resumed={}
        for owner in states:
            resumed[owner]=make(owner)
            resumed[owner].checkpoint(path(owner),True)
            assert resumed[owner].steps==8
        for owner,state in resumed.items():
            for step in range(9,33 if owner=='a' else 34):
                np.testing.assert_array_equal(advance(state,owner,step),rows[owner][step-1])
        for owner,state in states.items():
            A,j,index=original[owner]
            np.testing.assert_array_equal(state.A,A);np.testing.assert_array_equal(state.current,j)
            assert state.steps==index
            state.close()
        resumed['a'].close()
        model.close()
        assert context._dependents==1
        with np.testing.assert_raises(InputError):model.new_state('closed-model')
        np.testing.assert_array_equal(advance(resumed['b'],'b',34),rows['b'][33])
        resumed['b'].close()
        assert context._dependents==0
    else:
        owner=next(iter(states));states[owner].close()
        # All numerical states closed, but shared model can allocate again.
        for _ in range(0 if args.legacy else 2):
            fresh=make(owner)
            for step in range(1,4):
                np.testing.assert_array_equal(advance(fresh,owner,step),rows[owner][step-1])
            fresh.close()
        if model is not None:model.close()
        assert context._dependents==0
    for owner,trace in rows.items():
        if trace:np.save(out/f'{owner}-trace-{comm.rank}.npy',trace)
    if comm.rank==0:
        (out/'validation.json').write_text(json.dumps({'status':'pass','mode':args.mode,'group_size':comm.size,
            'work_copy_count':len(list(out.glob('native-fast-hse-*/bundle'))),'replayed_steps':0,
            'source':str(args.fixture),'shared_model_identity':model.identity if model else backend.identity,'legacy':args.legacy},indent=2)+'\n')
        assert len(list(out.glob('native-fast-hse-*/bundle')))==1
        print('Python shared fast-PR HSE PASS',args.mode,group,flush=True)
    comm.Free()


if __name__=='__main__':main()
