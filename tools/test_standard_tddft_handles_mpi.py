"""Same-model standard TDDFT C handles versus independent singleton controls.

BUILD OLD_BUILD GS NEW_OUTPUT; launch3*--group-size MPI ranks. Synthetic native
ownership/transaction/direct-restart regression, not material convergence.
"""
import argparse,ctypes as c,json,os,shutil
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import NativeParallelContext


def timer_probe(a,world):
    """Read existing GS in place; never copy/export GS or electronic checkpoints."""
    from salmonpy._native_group import _directory_call
    from salmonpy._shared_tddft import _SharedState
    group=world.rank//a.group_size;comm=world.Split(group,world.rank)
    record=json.loads((a.build if group==0 else a.old_build).read_text())
    work=a.output.resolve()/str(group)
    if comm.rank==0:
        work.mkdir(parents=True,exist_ok=False)
        text=(a.fixture/'rt.inp').read_text()
        (work/'data_for_restart').symlink_to(a.fixture.resolve()/'data_for_restart',target_is_directory=True)
        (work/'pseudo1_rps.dat').symlink_to(a.fixture.resolve()/'pseudo1_rps.dat')
        (work/'rt.inp').write_text(text)
    comm.Barrier()
    owner=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    native=NativeParallelContext(owner,record['library'],record['sha256']);lib=c.CDLL(record['library'])
    def call(fn,*args):return _directory_call(lib,comm,work,fn,*args,minimum_protocol=2)
    vec=np.ctypeslib.ndpointer(dtype=np.float64,shape=(3,),flags='C_CONTIGUOUS')
    dec_t=np.ctypeslib.ndpointer(dtype=np.int32,shape=(5,),flags='C_CONTIGUOUS')
    token=c.c_int64;zero=np.zeros(3);j=np.zeros(3)
    dec=np.array([1,a.group_size,1,1,1] if a.split=='orbital' else [1,1,a.group_size,1,1],dtype=np.int32)
    def configure(name,args):
        fn=getattr(lib,'sp_tddft_'+name);fn.argtypes=args;fn.restype=c.c_int;return fn
    def drive(i,k):return 1e-4*np.array([np.sin(.31*i),np.cos(.23*i)-1,np.sin(.17*i)])*(1 if k==1 else -.7)
    traces={};timer_checks=0;handles=[]
    host=None
    if a.host_unit_library and (group==0 or a.host_all_groups):
        import hashlib
        record_host=json.loads((a.host_unit_library/'result.json').read_text())
        path=Path(record_host['library'])
        assert hashlib.sha256(path.read_bytes()).hexdigest()==record_host['sha256']
        host=c.CDLL(str(path)).test_host_unit;host.argtypes=[c.c_char_p,c.c_int];host.restype=c.c_int
        assert host(os.fsencode(work/f'host-{comm.rank}.txt'),0)==0
    host_timer=None;host_timer_checks=0
    if a.host_timer_library and group==0:
        import hashlib
        hr=json.loads((a.host_timer_library/'result.json').read_text())
        assert hr['runtime_sha256']==record['sha256']
        assert hashlib.sha256(Path(hr['library']).read_bytes()).hexdigest()==hr['sha256']
        c.CDLL(record['library'],mode=c.RTLD_GLOBAL)  # Test-only Fortran module access on macOS.
        host_timer=c.CDLL(hr['library']).test_host_timers;host_timer.argtypes=[c.c_int];host_timer.restype=c.c_int
        assert host_timer(0)==0
    def check_host_timer():
        nonlocal host_timer_checks
        if host_timer:
            assert host_timer(1)==0,'native work changed host timer counters'
            host_timer_checks+=1
    try:
        if group or a.legacy_timers:
            init=configure('init',[c.c_char_p,c.c_double,c.c_int,vec,dec_t,vec])
            step=configure('step',[vec,vec,vec]);close=configure('close',[])
            assert call(init,b'rt.inp',.003,32,zero,dec,j)==0
            check_host_timer()
            old=zero.copy();rows=[];point=group or 1
            for i in range(1,33):
                new=drive(i,point)
                if group==0 and i==4:
                    bad=old.copy();bad[0]+=1.
                    assert call(step,bad,new,j)==10;check_host_timer()
                assert call(step,old,new,j)==0;rows.append(-j.copy());old=new
                check_host_timer()
            traces[point]=np.array(rows)
            if group==0 and a.legacy_timers:
                from salmonpy.embedded_tddft import _NativeState
                q=object.__new__(_NativeState);q.closed=False;q.failed=False;q.lib=lib;q.native=SimpleNamespace(context=owner);q._call=call
                assert q.timing_seconds['density']>0
            assert call(close)==0;check_host_timer()
        else:
            init=configure('handle_init',[c.c_char_p,c.c_double,c.c_int,vec,dec_t,c.POINTER(token),vec])
            new=configure('model_new',[c.POINTER(token),vec]);step=configure('handle_step',[token,vec,vec,vec])
            close=configure('handle_close',[token]);transaction=configure('handle_transaction',[token,c.c_int])
            def create(first=False):
                h=token();code=call(init,b'rt.inp',.003,32,zero,dec,c.byref(h),j) if first else call(new,c.byref(h),j)
                assert code==0,code
                handles.append(h.value);return h.value
            first=create(True);second=create()
            def query(h):
                state=object.__new__(_SharedState);state.closed=False;state.failed=False;state.handle=h
                state.native=SimpleNamespace(model=SimpleNamespace(lib=lib),context=owner);state._call=call
                return state.timing_seconds
            assert query(first)==query(first)
            last={1:zero.copy(),2:zero.copy()};rows={1:[],2:[]}
            for i in range(1,33):
                for k,h,sibling in [(1,first,second),(2,second,first)]:
                    before=query(sibling)
                    previous=query(h)
                    A=drive(i,k);assert call(step,h,last[k],A,j)==0
                    rows[k].append(j.copy());last[k]=A
                    assert query(sibling)==before
                    after=query(h);assert all(after[n]>=previous[n] for n in after)
                    timer_checks+=1
                if i==8:
                    assert call(transaction,first,0)==0
                    before=query(second)
                    outputs={p:p.read_bytes() for p in work.glob('*.data')}
                    assert call(step,first,last[1],drive(9,1),j)==0
                    assert call(transaction,first,1)==0
                    for path,data in outputs.items():path.write_bytes(data)
                    assert call(transaction,first,3)==0
                    assert call(transaction,first,2)==0
                    assert query(second)==before
            assert query(first)['density']>0 and query(second)['density']>0
            before=query(second);assert call(close,first)==0;handles.remove(first);assert query(second)==before
            fresh=create();assert query(fresh)['density']<query(second)['density']
            # New generation must not inherit old elapsed counters; rejected creation also preserves siblings.
            before=query(second);collision=work/'native_tddft_h3_g1_rt.data'
            if comm.rank==0:collision.write_text('keep')
            comm.Barrier();h=token();assert call(new,c.byref(h),j)==1;assert query(second)==before
            traces={k:np.array(v) for k,v in rows.items()}
            diagnostics={'first_and_sibling_isolated_checks':timer_checks,'survivor':query(second),'fresh':query(fresh)}
            if comm.rank==0:(work/'timers.json').write_text(json.dumps(diagnostics))
    finally:
        for h in reversed(handles):assert call(close,h)==0
        native.close()
        check_host_timer()
        if host:
            assert host(b'',1)==0,'native close disconnected host unit30'
            assert host(b'',2)==0
            assert (work/f'host-{comm.rank}.txt').read_text()=='host-before\nhost-after\n'
    gathered=world.allgather(traces)
    for k in ((1,) if a.legacy_timers else (1,2)):np.testing.assert_array_equal(gathered[0][k],gathered[k*a.group_size][k])
    if world.rank==0:(a.output/'result.json').write_text(json.dumps({'status':'passed','group_size':a.group_size,'split':a.split,'independent_controls_exact':True,'previous_runtime_compared':a.build.resolve()!=a.old_build.resolve(),'sibling_timer_checks':timer_checks,'GS_copies':0,'SCF_runs':0,'disk_electronic_checkpoints':0,'memory_rollback_exact':not a.legacy_timers,'host_timer_checks':host_timer_checks,'legacy_timers':a.legacy_timers,'host_unit30_preserved':bool(a.host_unit_library)})+'\n')
    comm.Free()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('build','old_build','fixture','output'):p.add_argument(key,type=Path)
    p.add_argument('--group-size',type=int,choices=(1,2),default=2)
    p.add_argument('--timers-only',action='store_true')
    p.add_argument('--host-unit-library',type=Path)
    p.add_argument('--host-all-groups',action='store_true')
    p.add_argument('--host-timer-library',type=Path)
    p.add_argument('--legacy-timers',action='store_true')
    p.add_argument('--split',choices=('orbital','rgrid'),default='orbital')
    a=p.parse_args();world=MPI.COMM_WORLD;assert world.size==3*a.group_size
    if a.timers_only:return timer_probe(a,world)
    group=world.rank//a.group_size;comm=world.Split(group,world.rank)
    record=json.loads((a.build if group==0 else a.old_build).read_text());lib=c.CDLL(record['library'])
    if group==0:assert hasattr(lib,'sp_tddft_handle_protocol'),'explicit native handles rebuild required'
    work=a.output.resolve()/str(group)
    if comm.rank==0:shutil.copytree(a.fixture,work)
    comm.Barrier();os.chdir(work)
    owner=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    native=NativeParallelContext(owner,record['library'],record['sha256'])
    vec=np.ctypeslib.ndpointer(dtype=np.float64,shape=(3,),flags='C_CONTIGUOUS')
    dec_t=np.ctypeslib.ndpointer(dtype=np.int32,shape=(5,),flags='C_CONTIGUOUS')
    token=c.c_int64;v=np.zeros(3);dec=np.array([1,a.group_size,1,1,1] if a.split=='orbital' else [1,1,a.group_size,1,1],dtype=np.int32)
    specs={'init':[c.c_char_p,c.c_double,c.c_int,vec,dec_t,vec],'step':[vec,vec,vec],
           'checkpoint':[c.c_char_p,c.c_int],'close':[]}
    for name,args in specs.items():fn=getattr(lib,'sp_tddft_'+name);fn.argtypes=args;fn.restype=c.c_int
    def drive(i,point):
        return 1e-4*np.array([np.sin(.31*i),np.cos(.23*i)-1,np.sin(.17*i)])*(1 if point==1 else -.7)
    traces={};initials={};states=[];counts=[]
    try:
        if group:
            j=np.zeros(3);assert lib.sp_tddft_init(b'rt.inp',.003,32,v,dec,j)==0
            initials[group]=-j.copy();values=[];old=v.copy()
            for i in range(1,33):
                new=drive(i,group);assert lib.sp_tddft_step(old,new,j)==0;values.append(-j.copy());old=new
            traces[group]=np.array(values);assert lib.sp_tddft_close()==0
        else:
            hspec={'protocol':[],'init':[c.c_char_p,c.c_double,c.c_int,vec,dec_t,c.POINTER(token),vec],
                   'step':[token,vec,vec,vec],'info':[token,c.POINTER(c.c_int),vec,vec],
                   'checkpoint':[token,c.c_char_p,c.c_int],'transaction':[token,c.c_int],
                   'label':[token,c.c_void_p,c.c_int],'close':[token]}
            h={}
            for name,args in hspec.items():
                fn=getattr(lib,'sp_tddft_handle_'+name);fn.argtypes=args;fn.restype=c.c_int;h[name]=fn
            assert h['protocol']()==1
            new=lib.sp_tddft_model_new;new.argtypes=[c.POINTER(token),vec];new.restype=c.c_int
            resource=lib.sp_tddft_resources;resource.argtypes=[np.ctypeslib.ndpointer(dtype=np.int32,shape=(4,),flags='C_CONTIGUOUS')];resource.restype=c.c_int
            def stats():
                values=np.zeros(4,dtype=np.int32);assert resource(values)==0;counts.append(values.tolist());return values
            def create(first=False):
                handle=token();j=np.zeros(3)
                code=h['init'](b'rt.inp',.003,32,v,dec,c.byref(handle),j) if first else new(c.byref(handle),j)
                assert code==0,('new',code)
                states.append(handle.value);return handle.value,j
            def info(handle):
                step=c.c_int();A=np.zeros(3);j=np.zeros(3);code=h['info'](handle,c.byref(step),A,j)
                assert code==0;return step.value,A,j
            def checkpoint(handle,name,mode):
                return h['checkpoint'](handle,os.fsencode(work/f'{name}-{comm.rank}.bin'),mode)
            def step(handle,i,point):
                _,old,_=info(handle);j=np.zeros(3);assert h['step'](handle,old,drive(i,point),j)==0;return j.copy()
            first,j=create(True);initials[1]=j.copy();second,j=create();initials[2]=j.copy()
            baseline=stats()[0];assert baseline>0
            extras=[]
            for unused in range(4):
                handle,initial=create();extras.append(handle)
                np.testing.assert_array_equal(initial,initials[1])
            assert stats()[0]>baseline
            for handle in reversed(extras):
                assert h['close'](handle)==0;states.remove(handle)
            assert stats()[0]==baseline
            label=c.create_string_buffer(257)
            assert h['label'](first,label,1 if comm.rank==0 else 257)==13
            if comm.size>1:
                st=c.c_int();A=np.zeros(3);j=np.zeros(3)
                assert h['info'](first if comm.rank==0 else second,c.byref(st),A,j)==11
            values={1:[],2:[]}
            for i in range(1,33):
                values[1].append(step(first,i,1));values[2].append(step(second,i,2))
                if i==4:
                    assert checkpoint(first,'a4',0)==0 and checkpoint(second,'b4',0)==0
                    assert h['transaction'](second,0)==0
                    label=c.create_string_buffer(257);assert h['label'](first,label,len(label))==0
                    names=[work/(label.value.decode()+tail) for tail in ('_rt.data','_rt_energy.data')]
                    assert h['transaction'](first,0)==0
                    saved={path:path.read_bytes() for path in names}
                    assert h['transaction'](second,2)==0
                    assert checkpoint(second,'b4',2)==0
                    step(first,5,1)
                    assert h['transaction'](first,1)==0
                    for path,data in saved.items():path.write_bytes(data)
                    assert h['transaction'](first,3)==0 and h['transaction'](first,2)==0
                    assert info(first)[0]==4 and info(second)[0]==4
                    assert checkpoint(first,'a4',2)==0  # Both owners have a pending candidate.
                    assert checkpoint(first,'unused',4)==0
                    old=info(second)[1];j=np.zeros(3)
                    assert h['step'](second,old,drive(5,2),j)==14
                    assert checkpoint(second,'unused',4)==0
                if i in (8,9):assert checkpoint(first if i==8 else second,'a8' if i==8 else 'b9',0)==0
            traces={k:np.array(value) for k,value in values.items()}
            assert h['transaction'](first,0)==0  # Survives sibling creation/load/close.
            stale=[]
            for point,split in ((1,8),(2,9)):
                sibling=info(second)
                fresh,_=create()
                for old_handle in stale:
                    st=c.c_int();A=np.zeros(3);j=np.zeros(3)
                    assert h['info'](old_handle,c.byref(st),A,j)==12
                assert checkpoint(fresh,'a8' if point==1 else 'b9',1)==0
                continued=[step(fresh,i,point) for i in range(split+1,33)]
                np.testing.assert_array_equal(continued,traces[point][split:])
                assert h['close'](fresh)==0;states.remove(fresh)
                stale.append(fresh)
                st=c.c_int();A=np.zeros(3);j=np.zeros(3);assert h['info'](fresh,c.byref(st),A,j)==12
                now=info(second);assert now[0]==sibling[0]
                np.testing.assert_array_equal(now[1],sibling[1]);np.testing.assert_array_equal(now[2],sibling[2])
            assert h['transaction'](first,2)==0
            assert checkpoint(second,'b9',2)==0
            assert h['close'](first)==0;states.remove(first)
            assert stats()[2]==1  # Closing A does not cancel B's pending candidate.
            assert h['close'](second)==0;states.remove(second)
            assert stats()[0]==0
            idle,_=create();np.testing.assert_array_equal(step(idle,1,1),traces[1][0]);assert h['close'](idle)==0;states.remove(idle)
            assert stats()[0]==0
            assert lib.sp_tddft_checkpoint(b'unused',4)==2
    finally:
        for handle in states:lib.sp_tddft_handle_close(handle)
        native.close()
    assert owner._dependents==0 and comm.allreduce(1)==comm.size
    values=world.allgather((traces,initials))
    for point in (1,2):
        ref,init=values[point*a.group_size]
        for r in range(a.group_size):
            np.testing.assert_array_equal(values[r][0][point],ref[point]);np.testing.assert_array_equal(values[r][1][point],init[point])
    for point,trace in traces.items():
        assert np.isfinite(trace).all() and np.max(abs(trace))>1e-12
        np.save(work/f'point{point}-{comm.rank}.npy',trace)
    if world.rank==0:
        (a.output/'result.json').write_text(json.dumps({'status':'pass','group_size':a.group_size,'split':a.split,
            'independent_old_exact':True,'resume_exact':True,'transactions_isolated':True,'stale_token_rejected':True,
            'replayed_steps':0,'counts':counts},indent=2)+'\n');print('Standard shared native handles and independent old controls PASS',flush=True)
    comm.Free()


if __name__=='__main__':
    try:main()
    except BaseException:
        import traceback
        traceback.print_exc();MPI.COMM_WORLD.Abort(1)
