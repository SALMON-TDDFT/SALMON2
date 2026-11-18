"""Two registered models vs independent controls. Three MPI groups, --size ranks each."""
import argparse
import ctypes as c
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.embedded_tddft import NativeTDDFTBackend,NativeTDDFTMaterialProvider
from salmonpy.errors import InputError


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('build','baseline','fixture_a','fixture_b','output'):p.add_argument(key,type=Path)
    p.add_argument('--size',type=int,default=2)
    p.add_argument('--b-split',choices=('orbital','rgrid'),default='rgrid')
    p.add_argument('--reverse',action='store_true')
    p.add_argument('--capacity-test',action='store_true')
    a=p.parse_args();world=MPI.COMM_WORLD
    assert world.size==3*a.size
    group=world.rank//a.size;comm=world.Split(group,world.rank)
    root=a.output.resolve()/str(group)
    if comm.rank==0:root.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    context=SimpleNamespace(solver_comm=comm,point_ids=(1,2) if group==0 else (group,),closed=False,_dependents=0)
    decs={1:(1,a.size,1,1,1),2:(1,1,1,1,a.size) if a.b_split=='rgrid' else (1,a.size,1,1,1)}
    dt={1:.003,2:.004};caps={1:48,2:64}
    fixtures={1:a.fixture_a,2:a.fixture_b}
    def drive(i,key):
        return 1e-4*np.array([np.sin(.31*i),np.cos(.23*i)-1,np.sin(.17*i)])*(1 if key==1 else -.7)
    steps_b=56 if a.capacity_test else 40
    traces={};models={};states={};registry=None
    if group:
        key=group
        provider=NativeTDDFTMaterialProvider(fixtures[key],'rt.inp',max_steps=caps[key])
        backend=NativeTDDFTBackend(a.baseline,root,decomposition=decs[key])
        state=backend.create_state(provider.prepare(context),key,dt[key],context)
        try:
            values=[]
            for i in range(1,(steps_b if key==2 else 40)+1):values.append(state.advance(state.A.copy(),drive(i,key)))
            traces[key]=np.asarray(values)
        finally:state.close()
    else:
        backend=NativeTDDFTBackend(a.build,root,decomposition=decs[1])
        registry=backend.open_registry(context)
        initials={k:NativeTDDFTMaterialProvider(v,'rt.inp',max_steps=caps[k]).prepare(context) for k,v in fixtures.items()}
        try:
            for k in ((2,1) if a.reverse else (1,2)):
                models[k]=registry.open_model(initials[k],dt[k],decomposition=decs[k])
            if comm.size>1:
                chosen=models[1] if comm.rank==0 else models[2]
                try:chosen.new_state('mixed-model')
                except InputError:pass
                else:raise AssertionError('mixed Python model selection accepted')
                assert models[1].lib.pending is not None and models[2].lib.pending is not None
            for k in (1,2):states[k]=models[k].new_state('material-'+str(k));traces[k]=[]
            for i in range(1,33):
                for k in (1,2):
                    state=states[k];traces[k].append(state.step(drive(i,k)))
                    assert state.time==i*dt[k]
                    if (k,i) in ((1,8),(2,9)):state.checkpoint(root/f'{k}-{comm.rank}.bin')
            for k,split in ((1,8),(2,9)):
                fresh=models[k].new_state('resumed-'+str(k))
                try:
                    fresh.checkpoint(root/f'{k}-{comm.rank}.bin',restore=True)
                    for i in range(split+1,33):np.testing.assert_array_equal(fresh.step(drive(i,k)),traces[k][i-1])
                finally:fresh.close()
            old_model_token=models[1].model_token
            models[1].close()
            with_error=False
            try:models[1].new_state()
            except InputError:with_error=True
            assert with_error
            traces[1].append(states[1].step(drive(33,1)))
            states[1].close();states.pop(1)
            for i in range(33,steps_b+1):traces[2].append(states[2].step(drive(i,2)))
            states[2].close();states.pop(2)
            idle=models[2].new_state('idle')
            try:np.testing.assert_array_equal(idle.step(drive(1,2)),traces[2][0])
            finally:idle.close()
            models[2].close()
            # Reuse registry slots and distinguish model/state/stale/mixed tokens.
            first=registry.open_model(initials[1],dt[1],decomposition=decs[1]);models[3]=first
            second=registry.open_model(initials[2],dt[2],decomposition=decs[2]);models[4]=second
            child=first.new_state('survives-registry-close');states[3]=child
            np.testing.assert_array_equal(child.step(drive(1,1)),traces[1][0])
            handle=c.c_int64();current=np.zeros(3)
            assert registry.lib.sp_tddft_model_spawn(old_model_token,c.byref(handle),current)==12
            assert registry.lib.sp_tddft_model_spawn(child.handle,c.byref(handle),current)==12
            if comm.size>1:
                token=first.model_token if comm.rank==0 else second.model_token
                assert registry.lib.sp_tddft_model_spawn(token,c.byref(handle),current)==11
            # Last default model differs from active owner; failed register must retain both.
            mt,st=c.c_int64(),c.c_int64()
            code=registry._call(first.work,registry.lib.sp_tddft_model_register,b'rt.inp',float('nan'),16,
                                np.zeros(3),np.asarray(decs[1],dtype=np.int32),c.byref(mt),c.byref(st),current)
            assert code==10 and mt.value==0 and st.value==0
            sibling=second.new_state('after-rejected-register')
            try:np.testing.assert_array_equal(sibling.step(drive(1,2)),traces[2][0])
            finally:sibling.close()
            queued=registry.open_model(initials[2],dt[2],decomposition=decs[2]);models[5]=queued
            class RankOrder(list):
                def discard(self,value):
                    if value in self:self.remove(value)
            ordered=list(registry._models)
            ordered.sort(key=lambda model:model.model_token,reverse=bool(comm.rank%2))
            registry._models=RankOrder(ordered)
            registry.close()
            np.testing.assert_array_equal(child.step(drive(2,1)),traces[1][1])
            child.close();states.pop(3)
            assert len(list(root.glob('model-*/bundle/data_for_restart')))==5
            traces={k:np.asarray(v) for k,v in traces.items()}
        finally:
            for state in reversed(tuple(states.values())):state.close()
            for model in tuple(models.values()):model.close()
            registry.close()
    assert context._dependents==0
    all_traces=world.allgather(traces)
    for k in (1,2):np.testing.assert_array_equal(all_traces[0][k],all_traces[k*a.size][k][:len(all_traces[0][k])])
    np.savez(root/f'trace-{comm.rank}.npz',**{f'p{k}':v for k,v in traces.items()})
    if world.rank==0:
        (a.output/'result.json').write_text(json.dumps({'status':'pass','size':a.size,'b_split':a.b_split,
            'reverse':a.reverse,'accepted_B_steps':steps_b,'independent_exact':True,'resume_exact':True,'replayed_steps':0,
            'model_close_preserves_siblings':True,'stale_and_mixed_model_tokens_rejected':True,
            'registry_close_order_independent':True,'failed_register_preserves_models':True})+'\n')
    comm.Free()


if __name__=='__main__':
    try:main()
    except BaseException:
        import traceback
        traceback.print_exc();MPI.COMM_WORLD.Abort(1)
