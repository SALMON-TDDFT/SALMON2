"""Shared reduced/full-k HSE numerical probe: BUILD FIXTURE NEW_OUTPUT.

Launch 3*--group-size ranks. Group0 shares A/B; groups1/2 use independent legacy
states as controls. Every context sees a symmetry-preserving Cartesian z drive.
Explicit local build/inputs only; synthetic ownership test, not convergence.
"""
import argparse
import ctypes as c,re,hashlib
from salmonpy.distributed import NativeParallelContext
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.fast_hse import FastHSEMaterialProvider
from salmonpy.fast_hse_mpi import NativeFastHSEBackend, _history_info


def timer_probe(a,world):
    """Read existing GS in place; never copy/export GS or electronic checkpoints."""
    from salmonpy._native_group import _directory_call
    from salmonpy.fast_hse_mpi import _State
    group=world.rank//a.group_size;comm=world.Split(group,world.rank)
    record=json.loads((a.build if group==0 else a.baseline).read_text())
    work=a.output.resolve()/str(group)
    heterogeneous=bool(a.second_fixture)
    fixtures=[a.fixture,Path(a.second_fixture)] if heterogeneous else [a.fixture,a.fixture]
    input_works=[work/'a',work/'b'] if heterogeneous and group==0 else [work,work]
    if comm.rank==0:
        work.mkdir(parents=True,exist_ok=False)
        for destination,fixture in (zip(input_works,fixtures) if group==0 and heterogeneous else [(work,fixtures[group-1] if group else a.fixture)]):
            if destination!=work:destination.mkdir()
            text=(fixture/'rt.inp').read_text()
            text=re.sub(r'(?i)\bnproc_k\s*=\s*\d+',f'nproc_k={a.group_size}',text)
            (destination/'data_for_restart').symlink_to(fixture.resolve()/'data_for_restart',target_is_directory=True)
            for asset in fixture.iterdir():
                if asset.is_file() and asset.name!='rt.inp':
                    (destination/asset.name).symlink_to(asset.resolve())
            (destination/'rt.inp').write_text(text)
    comm.Barrier()
    owner=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    native=NativeParallelContext(owner,record['library'],record['sha256']);lib=c.CDLL(record['library'])
    state_work={};model_work={};initial_work=work
    def call(fn,*args):
        destination=initial_work
        if args and isinstance(args[0],int):
            destination=state_work.get(args[0],model_work.get(args[0],destination))
        return _directory_call(lib,comm,destination,fn,*args,minimum_protocol=2)
    vec=np.ctypeslib.ndpointer(dtype=np.float64,shape=(3,),flags='C_CONTIGUOUS')
    token=c.c_int64;zero=np.zeros(3);j=np.zeros(3);energy=c.c_double()
    def configure(name,args):
        fn=getattr(lib,'sp_fast_hse_'+name);fn.argtypes=args;fn.restype=c.c_int;return fn
    def drive(i,k):return np.array([0.,0.,1e-4*np.sin(.3*i)*(1 if k==1 else -.7)])
    traces={};timer_checks=0;handles=[]
    host_timer=None;host_timer_checks=0
    if a.host_timer_library and group==0:
        hr=json.loads((a.host_timer_library/'result.json').read_text())
        assert hr['runtime_sha256']==record['sha256']
        assert hashlib.sha256(Path(hr['library']).read_bytes()).hexdigest()==hr['sha256']
        c.CDLL(record['library'],mode=c.RTLD_GLOBAL)
        host_timer=c.CDLL(hr['library']).test_host_timers
        host_timer.argtypes=[c.c_int];host_timer.restype=c.c_int
        assert host_timer(0)==0
    def check_host_timer():
        nonlocal host_timer_checks
        if host_timer:
            assert host_timer(1)==0,'HSE work changed host timer counters'
            host_timer_checks+=1
    try:
        if group or a.legacy_timers:
            init=configure('init_symmetry' if a.symmetry else 'init',[c.c_char_p,c.c_double,c.c_int,vec])
            step=configure('step',[vec,vec,vec,c.POINTER(c.c_double)]);close=configure('close',[])
            history=configure('history_info',[c.c_void_p])
            assert call(init,b'rt.inp',.001,32,j)==0
            check_host_timer()
            old=zero.copy();rows=[];point=group or 1
            for i in range(1,33):
                new=drive(i,point)
                if group==0 and i==4:
                    stale=old.copy();stale[0]+=1.
                    assert call(step,stale,new,j,c.byref(energy))!=0
                    check_host_timer()
                assert call(step,old,new,j,c.byref(energy))==0
                check_host_timer()
                hist=np.zeros(8,dtype=np.int32);assert call(history,hist.ctypes.data)==0
                rows.append(np.r_[j,energy.value,hist]);old=new
            traces[point]=np.array(rows)
            if group==0 and a.legacy_timers:
                q=object.__new__(_State);q.closed=False;q.failed=False;q.lib=lib
                q._native=SimpleNamespace(context=owner);q._call=call
                assert q.timing_seconds['density']>0
            assert call(close)==0;check_host_timer()
        else:
            init=configure('handle_init',[c.c_char_p,c.c_double,c.c_int,c.c_int,c.POINTER(token),vec])
            new=configure('model_new',[c.c_double,c.c_int,c.POINTER(token),vec]);step=configure('handle_step',[token,vec,vec,vec,c.POINTER(c.c_double)])
            history=configure('handle_history',[token,c.c_void_p])
            model_token=token();model_tokens=[]
            register=configure('model_register',[c.c_char_p,c.c_double,c.c_int,c.c_int,c.POINTER(token),c.POINTER(token),vec])
            spawn=configure('model_spawn',[token,c.c_double,c.c_int,c.POINTER(token),vec])
            close=configure('handle_close',[token]);transaction=configure('handle_transaction',[token,c.c_int])
            def create(first=False):
                nonlocal initial_work
                h=token()
                if heterogeneous and len(model_tokens)<2:
                    initial_work=input_works[len(model_tokens)]
                    code=call(register,b'rt.inp',.001,32,1,c.byref(model_token),c.byref(h),j)
                    model_tokens.append(model_token.value);model_work[model_token.value]=initial_work
                elif first:
                    code=call(register,b'rt.inp',.001,32,int(a.symmetry),c.byref(model_token),c.byref(h),j) if a.timer_registry else call(init,b'rt.inp',.001,32,int(a.symmetry),c.byref(h),j)
                else:
                    code=call(spawn,model_token.value,.001,32,c.byref(h),j) if a.timer_registry else call(new,.001,32,c.byref(h),j)
                assert code==0,code
                if not heterogeneous and a.timer_registry and first:model_tokens.append(model_token.value)
                state_work[h.value]=initial_work if len(handles)<2 else state_work.get(handles[0],initial_work)
                handles.append(h.value);return h.value
            first=create(True);second=create()
            def query(h):
                state=object.__new__(_State);state.closed=False;state.failed=False;state.handle=h
                state._native=SimpleNamespace(context=owner);state.lib=lib;state._call=call
                return state.timing_seconds
            assert query(first)==query(first)
            last={1:zero.copy(),2:zero.copy()};rows={1:[],2:[]}
            for i in range(1,33):
                for k,h,sibling in [(1,first,second),(2,second,first)]:
                    before=query(sibling)
                    previous=query(h)
                    A=drive(i,k);assert call(step,h,last[k],A,j,c.byref(energy))==0
                    hist=np.zeros(8,dtype=np.int32);assert call(history,h,hist.ctypes.data)==0
                    rows[k].append(np.r_[j,energy.value,hist]);last[k]=A
                    assert query(sibling)==before
                    after=query(h);assert all(after[n]>=previous[n] for n in after)
                    timer_checks+=1
                if i==8:
                    assert call(transaction,first,0)==0
                    before=query(second)
                    outputs={p:p.read_bytes() for p in work.rglob('*.data')}
                    assert call(step,first,last[1],drive(9,1),j,c.byref(energy))==0
                    assert call(transaction,first,1)==0
                    for path,data in outputs.items():path.write_bytes(data)
                    assert call(transaction,first,3)==0
                    assert call(transaction,first,2)==0
                    assert query(second)==before
            assert query(first)['density']>0 and query(second)['density']>0
            before=query(second);assert call(close,first)==0;handles.remove(first);assert query(second)==before
            switch_history=np.zeros(8,dtype=np.int32);assert call(history,second,switch_history.ctypes.data)==0
            if heterogeneous:
                model_token.value=model_tokens[0];initial_work=model_work[model_token.value]
            fresh=create();assert query(fresh)['density']<query(second)['density']
            # New generation must not inherit old elapsed counters; rejected creation also preserves siblings.
            before=query(second);collision=(model_work[model_tokens[0]] if heterogeneous else work)/'native_hse_h3_g1_rt.data'
            if comm.rank==0:collision.write_text('keep')
            comm.Barrier();h=token()
            code=call(spawn,model_token.value,.001,32,c.byref(h),j) if a.timer_registry else call(new,.001,32,c.byref(h),j)
            assert code==comm.size,('collision admission status',code,h.value);assert query(second)==before
            other_history=np.zeros(8,dtype=np.int32);assert call(history,second,other_history.ctypes.data)==0
            assert call(step,fresh,zero,drive(1,1),j,c.byref(energy))==0
            fresh_history=np.zeros(8,dtype=np.int32);assert call(history,fresh,fresh_history.ctypes.data)==0
            np.testing.assert_array_equal(np.r_[j,energy.value,fresh_history],rows[1][0])
            if heterogeneous:
                assert call(close,fresh)==0;handles.remove(fresh)
                release=configure('model_close',[token]);old_model=model_tokens.pop(0)
                assert call(release,old_model)==0
                before_reuse=query(second);reused=create()
                assert call(step,reused,zero,drive(1,2),j,c.byref(energy))==0
                reused_history=np.zeros(8,dtype=np.int32);assert call(history,reused,reused_history.ctypes.data)==0
                np.testing.assert_array_equal(np.r_[j,energy.value,reused_history],rows[2][0])
                assert query(second)==before_reuse
                fresh=reused
            traces={k:np.array(v) for k,v in rows.items()}
            diagnostics={'first_and_sibling_isolated_checks':timer_checks,'survivor':query(second),'fresh':query(fresh)}
            if comm.rank==0:(work/'timers.json').write_text(json.dumps(diagnostics))
    finally:
        for h in reversed(handles):assert call(close,h)==0
        if group==0 and a.timer_registry:
            release=configure('model_close',[token])
            for model in model_tokens:assert call(release,model)==0
        native.close()
        check_host_timer()
    gathered=world.allgather(traces)
    for k in ((1,) if a.legacy_timers else (1,2)):np.testing.assert_array_equal(gathered[0][k],gathered[k*a.group_size][k])
    if world.rank==0:(a.output/'result.json').write_text(json.dumps({'status':'passed','group_size':a.group_size,'symmetry':a.symmetry,'registry':a.timer_registry,'heterogeneous':heterogeneous,'fresh_after_switch_exact':True,'model_generation_reuse_exact':heterogeneous,'independent_old_current_energy_history_exact':True,'sibling_timer_checks':timer_checks,'GS_copies':0,'SCF_runs':0,'disk_electronic_checkpoints':0,'memory_rollback_exact':None if a.legacy_timers else True,'host_timer_checks':host_timer_checks,'legacy_timers':a.legacy_timers})+'\n')
    comm.Free()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build','fixture','output'):p.add_argument(name,type=Path)
    p.add_argument('--group-size',type=int,default=2)
    p.add_argument('--timers-only',action='store_true')
    p.add_argument('--legacy-timers',action='store_true')
    p.add_argument('--host-timer-library',type=Path)
    p.add_argument('--timer-registry',action='store_true')
    p.add_argument('--second-fixture',type=Path)
    p.add_argument('--symmetry',action='store_true')
    p.add_argument('--steps',type=int,default=32)
    p.add_argument('--baseline',type=Path)
    p.add_argument('--no-chdir',action='store_true')
    args=p.parse_args();world=MPI.COMM_WORLD
    assert world.size==3*args.group_size and args.steps>=16
    if args.timers_only:return timer_probe(args,world)
    group=world.rank//args.group_size;comm=world.Split(group,world.rank)
    out=args.output.resolve();work=out/str(group)
    if comm.rank==0:work.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    context=SimpleNamespace(solver_comm=comm,point_ids=(1,2) if group==0 else (group,),closed=False,_dependents=0)
    anchor=Path.cwd()
    from unittest.mock import patch
    guard=patch('os.chdir',side_effect=AssertionError('shared HSE host chdir forbidden')) if args.no_chdir and group==0 else None
    if guard:guard.start()
    initial=FastHSEMaterialProvider(args.fixture,max_steps=args.steps).prepare(context)
    backend=NativeFastHSEBackend(args.baseline if group and args.baseline else args.build,work,allow_symmetry=args.symmetry)
    states=[];model=None;traces={};histories={};resumes={}
    def drive(index,which):
        return np.array([0.,0.,1e-4*np.sin(.3*index)*(1 if which==1 else -.7)])
    try:
        if group==0:
            model=backend.open_model(initial,.001,context)
            states=[model.new_state(1),model.new_state(2)]
        else:states=[backend.create_state(initial,group,.001,context)]
        if group==0 and args.symmetry:
            from salmonpy.errors import InputError
            before=[(s.current.copy(),s.A.copy(),s.steps) for s in states]
            try:states[1].step([1e-4,0.,0.])
            except InputError as exc:assert 'breaks' in str(exc)
            else:raise AssertionError('symmetry-breaking A was accepted')
            for state,(j,A,step) in zip(states,before):
                np.testing.assert_array_equal(state.current,j);np.testing.assert_array_equal(state.A,A)
                assert state.steps==step and not state.failed
        for state in states:traces[state.point_id if group==0 else group]=[]
        # Independent control create_state has no point_id attribute.
        for index in range(1,args.steps+1):
            for offset,state in enumerate(states):
                which=offset+1 if group==0 else group
                traces[which].append(np.r_[state.step(drive(index,which)),state.energy])
                if group==0 and index==7+which:
                    state.checkpoint(work/f'{which}-{comm.rank}.bin',False)
        if group==0:
            for state in states:histories[state.point_id]=state._call(_history_info,state.lib).tolist()
            saved_b=(states[1].A.copy(),states[1].current.copy(),states[1].energy,states[1].steps)
            for which in (1,2):
                fresh=model.new_state(which);states.append(fresh)
                fresh.checkpoint(work/f'{which}-{comm.rank}.bin',True)
                trace=[]
                for index in range(8+which,args.steps+1):trace.append(np.r_[fresh.step(drive(index,which)),fresh.energy])
                np.testing.assert_array_equal(trace,traces[which][7+which:])
                resumes[which]=trace
            np.testing.assert_array_equal(states[1].A,saved_b[0]);np.testing.assert_array_equal(states[1].current,saved_b[1])
            assert (states[1].energy,states[1].steps)==saved_b[2:]
        if group==0:
            for state in reversed(states):state.close()
            states=[]
            fresh=model.new_state(3);states.append(fresh)
            np.testing.assert_array_equal(np.r_[fresh.step(drive(1,1)),fresh.energy],traces[1][0])
            model.close()  # A surviving lease still works after creation is disabled.
            np.testing.assert_array_equal(np.r_[fresh.step(drive(2,1)),fresh.energy],traces[1][1])
        np.savez(work/f'trace-{comm.rank}.npz',**{f'point{k}':v for k,v in traces.items()})
    finally:
        for state in reversed(states):state.close()
        if model is not None:model.close()
    assert context._dependents==0
    assert Path.cwd()==anchor
    if guard:guard.stop()
    gathered=world.allgather(traces)
    for which in (1,2):
        ref=gathered[which*args.group_size][which]
        for r in range(args.group_size):np.testing.assert_array_equal(gathered[r][which],ref)
    if world.rank==0:
        assert not np.array_equal(traces[1],traces[2])
        (out/'result.json').write_text(json.dumps({'status':'pass','symmetry':args.symmetry,
            'group_size':args.group_size,'steps':args.steps,'replayed_steps':0,
            'independent_controls_exact':True,'new_runtime_no_chdir':args.no_chdir,'previous_runtime_compared':bool(args.baseline),'resume_exact':True,'histories':histories},indent=2)+'\n')
        print('Shared symmetry HSE: controls and direct resume PASS',flush=True)
    comm.Free()


if __name__=='__main__':
    try:main()
    except BaseException:
        import traceback
        traceback.print_exc()
        MPI.COMM_WORLD.Abort(1)
