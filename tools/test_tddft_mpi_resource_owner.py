"""Raw standard TDDFT close returns initialization MPI resources to host.

BUILD FIXTURE NEW_OUTPUT; launch2 MPI ranks. --baseline records old-runtime
numerics without new diagnostics. No build/download. One native lifetime remains.
"""
import argparse
import ctypes as c
import json
import os
from pathlib import Path
import shutil
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import NativeParallelContext


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build','fixture','output'):p.add_argument(name,type=Path)
    p.add_argument('--baseline',action='store_true')
    p.add_argument('--split',choices=('orbital','rgrid'),default='orbital')
    a=p.parse_args();comm=MPI.COMM_WORLD;assert comm.size==2
    record=json.loads(a.build.read_text());lib=c.CDLL(record['library'])
    if not a.baseline:
        q=getattr(lib,'sp_tddft_context_mpi_protocol',None)
        assert q is not None,'explicit rebuild required: standard TDDFT MPI ownership protocol missing'
        q.argtypes=[];q.restype=c.c_int;assert q()==1
    root=a.output.resolve()
    if comm.rank==0:shutil.copytree(a.fixture,root/'bundle')
    comm.Barrier()
    owner=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    native=NativeParallelContext(owner,record['library'],record['sha256'])
    vector=np.ctypeslib.ndpointer(dtype=np.float64,shape=(3,),flags='C_CONTIGUOUS')
    integers=np.ctypeslib.ndpointer(dtype=np.int32,shape=(5,),flags='C_CONTIGUOUS')
    lib.sp_tddft_init.argtypes=[c.c_char_p,c.c_double,c.c_int,vector,integers,vector]
    lib.sp_tddft_init.restype=c.c_int
    lib.sp_tddft_step.argtypes=[vector,vector,vector];lib.sp_tddft_step.restype=c.c_int
    lib.sp_tddft_checkpoint.argtypes=[c.c_char_p,c.c_int];lib.sp_tddft_checkpoint.restype=c.c_int
    lib.sp_tddft_close.restype=c.c_int
    states=[]
    def counts():
        if a.baseline:return None
        values=np.zeros(4,dtype=np.int32)
        fn=lib.sp_tddft_resources
        fn.argtypes=[np.ctypeslib.ndpointer(dtype=np.int32,shape=(4,),flags='C_CONTIGUOUS')];fn.restype=c.c_int
        assert fn(values)==0
        states.append(values.tolist());return values
    before=counts()
    if before is not None:np.testing.assert_array_equal(before,[0,0,0,0])
    previous=Path.cwd();initialized=False
    try:
        os.chdir(root/'bundle')
        old=np.zeros(3);j=np.zeros(3)
        decomposition=np.array([1,2,1,1,1] if a.split=='orbital' else [1,1,2,1,1],dtype=np.int32)
        assert lib.sp_tddft_init(b'rt.inp',.003,32,old,decomposition,j)==0
        initialized=True;live=counts()
        if live is not None:assert live[0]>0 and live[2]==0 and live[3]==1
        trace=[-j.copy()]
        for index in range(1,33):
            new=np.array([0.,0.,1e-4*np.sin(.3*index)])
            assert lib.sp_tddft_step(old,new,j)==0
            old=new;trace.append(-j.copy())
        path=os.fsencode(root/f'state-{comm.rank}.bin')
        assert lib.sp_tddft_checkpoint(path,0)==0
        assert lib.sp_tddft_checkpoint(path,2)==0  # Close with a staged, uncommitted payload.
        staged=counts()
        if staged is not None:assert staged[2]==1 and staged[0]==live[0]
        assert lib.sp_tddft_close()==0;initialized=False
        closed=counts()
        if closed is not None:assert closed[0]==0 and closed[2]==0 and closed[3]==0
        # Resource retirement does not silently enable unsafe process reinitialization.
        assert lib.sp_tddft_init(b'rt.inp',.003,32,old,decomposition,j)==12
        assert comm.allreduce(1)==comm.size
        np.save(root/f'current-{comm.rank}.npy',trace)
    finally:
        if initialized:lib.sp_tddft_close()
        os.chdir(previous);native.close()
    assert owner._dependents==0
    if comm.rank==0:
        (root/'result.json').write_text(json.dumps({'status':'pass','split':a.split,'baseline':a.baseline,
            'resource_counts':states,'host_communicator_alive':True,'one_lifetime_preserved':True,
            'staged_payload_released':not a.baseline},indent=2)+'\n')
        print('Standard TDDFT close, staged payload and host MPI PASS',flush=True)


if __name__=='__main__':
    try:main()
    except BaseException:
        import traceback
        traceback.print_exc()
        MPI.COMM_WORLD.Abort(1)
