"""Test repeated shared HSE state close/recreate without retaining MPI children.

BUILD FIXTURE NEW_OUTPUT; launch 2*--group-size MPI ranks. Group0 repeatedly
creates/closes a sibling while group1 supplies an independent propagation control.
Explicit current runtime only. No build/download or material-convergence claim.
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


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('build','fixture','output'):p.add_argument(key,type=Path)
    p.add_argument('--group-size',type=int,choices=(1,2),default=2)
    p.add_argument('--cycles',type=int,default=16)
    p.add_argument('--symmetry',action='store_true')
    a=p.parse_args();world=MPI.COMM_WORLD
    assert world.size==2*a.group_size and a.cycles>=2
    group=world.rank//a.group_size;comm=world.Split(group,world.rank)
    work=a.output.resolve()/str(group)
    if comm.rank==0:work.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    ctx=SimpleNamespace(solver_comm=comm,point_ids=(1,2) if group==0 else (1,),closed=False,_dependents=0)
    backend=NativeFastHSEBackend(a.build,work,allow_symmetry=a.symmetry)
    initial=FastHSEMaterialProvider(a.fixture,max_steps=a.cycles+2).prepare(ctx)
    states=[];model=None;trace=[];counts=[]
    try:
        model=backend.open_model(initial,.001,ctx)
        query=getattr(model.lib,'sp_fast_hse_context_mpi_protocol',None)
        assert query is not None,'explicit rebuild required: context-owned MPI protocol missing'
        query.argtypes=[];query.restype=c.c_int;assert query()==1
        stats=model.lib.sp_fast_hse_mpi_resources
        vector=np.ctypeslib.ndpointer(dtype=np.int32,shape=(2,),flags='C_CONTIGUOUS')
        stats.argtypes=[vector];stats.restype=c.c_int
        def count():
            value=np.zeros(2,dtype=np.int32)
            assert stats(value)==0
            counts.append(value.tolist());return value.copy()
        survivor=model.new_state(1);states.append(survivor)
        baseline=count();assert baseline[0]>0
        for step in range(1,a.cycles+1):
            if group==0:
                sibling=model.new_state(2);states.append(sibling)
                peak=count();assert peak[0]>baseline[0]
                sibling.step([0.,0.,-7e-5*np.sin(.2*step)])
                sibling.close();states.remove(sibling)
                after=count();assert after[0]==baseline[0],(baseline,after)
                assert after[1]<=2*baseline[1],(baseline,after)
            trace.append(np.r_[survivor.step([0.,0.,1e-4*np.sin(.3*step)]),survivor.energy])
        survivor.close();states.remove(survivor)
        assert count()[0]==0,'closed context children retained'
        fresh=model.new_state(3);states.append(fresh)
        np.testing.assert_array_equal(np.r_[fresh.step([0.,0.,1e-4*np.sin(.3)]),fresh.energy],trace[0])
        fresh.close();states.remove(fresh);assert count()[0]==0
    finally:
        for state in reversed(states):state.close()
        if model is not None:model.close()
    assert ctx._dependents==0
    assert comm.allreduce(1)==comm.size  # Host communicator still alive.
    values=world.allgather(trace)
    for value in values:np.testing.assert_array_equal(value,values[0])
    if world.rank==0:
        (a.output/'result.json').write_text(json.dumps({'status':'pass','cycles':a.cycles,
            'group_size':a.group_size,'symmetry':a.symmetry,'counts':counts,
            'independent_control_exact':True,'host_communicator_alive':True},indent=2)+'\n')
        print('HSE context MPI retirement, bounded registry and independent control PASS',flush=True)
    comm.Free()


if __name__=='__main__':
    try:main()
    except BaseException:
        import traceback
        traceback.print_exc()
        MPI.COMM_WORLD.Abort(1)
