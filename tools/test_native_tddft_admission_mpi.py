"""Two-rank raw TDDFT admission/retry checks. BUILD GS NEW_OUTPUT [--baseline].

Use an MPI timeout and OMP/BLAS=1. No Python input validator shields the raw ABI.
The baseline route uses only valid calls and records the identical drive.
"""
import ctypes as c
import json
import os
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import NativeParallelContext


def main():
    build, fixture, destination = sys.argv[1:4]
    baseline = '--baseline' in sys.argv[4:]
    world = MPI.COMM_WORLD
    split = '--split' in sys.argv[4:]
    assert world.size == (4 if split else 2)
    group = world.rank//2 if split else 0
    comm = world.Split(group,world.rank) if split else world
    baseline = baseline or (split and group == 1)
    record = json.loads(Path(build).read_text())
    lib = c.CDLL(record['library'])
    if not baseline:
        query = lib.sp_tddft_error_protocol
        query.restype = c.c_int
        assert query() == 1
    root = Path(destination).resolve()
    if split:root = root/('group-'+str(group))
    if comm.rank == 0:
        root.mkdir(parents=True, exist_ok=False)
        shutil.copytree(fixture, root/'bundle')
    comm.Barrier()
    if '--python-init-retry' in sys.argv[4:]:
        from salmonpy.embedded_tddft import NativeTDDFTMaterialProvider, NativeTDDFTBackend
        from salmonpy.errors import InputError
        context=SimpleNamespace(solver_comm=comm,point_ids=(1,))
        initial=NativeTDDFTMaterialProvider(root/'bundle','rt.inp',max_steps=8).prepare(context)
        backend=NativeTDDFTBackend(build,root,decomposition=(1,2,1,1,1))
        try:
            backend.create_state(initial,1,.003*(comm.rank+1),context)
        except InputError:
            pass
        else:
            raise AssertionError('different initial dt accepted')
        state=backend.create_state(initial,1,.003,context)
        currents=[]
        for step in range(1,9):
            value=np.array([.001*np.sin(.2*step),0.,0.])
            currents.append(-state.advance(state.A.copy(),value))
        state.close()
        np.save(root/f'currents-{comm.rank}.npy',currents)
        if comm.rank==0:
            (root/'validation.json').write_text(json.dumps({'status':'passed','ranks':2,'python_init_retry':True,'steps':8},indent=2)+'\n')
        return
    owner = SimpleNamespace(solver_comm=comm, closed=False, _dependents=0)
    native = NativeParallelContext(owner, record['library'], record['sha256'])
    vector = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, shape=(3,), flags='C_CONTIGUOUS')
    layout = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, shape=(5,), flags='C_CONTIGUOUS')
    lib.sp_tddft_init.argtypes = [c.c_char_p,c.c_double,c.c_int,vector,layout,vector]
    lib.sp_tddft_init.restype = c.c_int
    lib.sp_tddft_step.argtypes = [vector,vector,vector]
    lib.sp_tddft_step.restype = c.c_int
    lib.sp_tddft_checkpoint.argtypes = [c.c_char_p,c.c_int]
    lib.sp_tddft_checkpoint.restype = c.c_int
    lib.sp_tddft_close.restype = c.c_int
    os.chdir(root/'bundle')
    A = np.zeros(3); J = np.zeros(3); dec = np.array([1,2,1,1,1],dtype=np.int32)
    rejected = []
    def expect(code, wanted, name):
        assert comm.allgather(code) == [wanted]*2, (name,code)
        rejected.append(name)
    if not baseline:
        expect(lib.sp_tddft_step(A,A,J),12,'step-before-init')
        expect(lib.sp_tddft_init(b'rt.inp',float('nan') if comm.rank else .003,8,A,dec,J),10,'init-nonfinite')
        expect(lib.sp_tddft_init(b'rt.inp',.003,0,A,dec,J),10,'init-capacity')
        bad = dec.copy();bad[0] = 0
        expect(lib.sp_tddft_init(b'rt.inp',.003,8,A,bad,J),10,'init-layout')
        expect(lib.sp_tddft_init(b'missing.inp',.003,8,A,dec,J),10,'init-missing-file')
        expect(lib.sp_tddft_init(b'rt.inp',.003*(comm.rank+1),8,A,dec,J),11,'init-disagreement')
    assert lib.sp_tddft_init(b'rt.inp',.003,8,A,dec,J) == 0
    def bytes_now(label):
        path = root/f'{label}-{comm.rank}.bin'
        assert lib.sp_tddft_checkpoint(os.fsencode(path),0) == 0
        return path.read_bytes()
    def reject_step(old,new,status,name):
        before=bytes_now(name+'-before')
        expect(lib.sp_tddft_step(old,new,J),status,name)
        assert bytes_now(name+'-after') == before, name
    if not baseline:
        before=bytes_now('init-reuse-before')
        expect(lib.sp_tddft_init(b'rt.inp',.003,8,A,dec,J),12,'init-reuse')
        assert bytes_now('init-reuse-after')==before
        seed=root/f'pending-seed-{comm.rank}.bin'
        bytes_now('pending-seed')
        assert lib.sp_tddft_checkpoint(os.fsencode(seed),2)==0
        reject_step(A,A,14,'step-pending-candidate')
        assert lib.sp_tddft_checkpoint(b'',4)==0
        reject_step(A,np.array([np.nan if comm.rank else 0.,0.,0.]),10,'step-nonfinite')
        reject_step(np.ones(3),A,10,'step-stale')
        reject_step(A,np.array([.001*comm.rank,0.,0.]),11,'step-disagreement')
    currents=[]
    for step in range(1,9):
        new=np.array([.001*np.sin(.2*step),0.,0.])
        assert lib.sp_tddft_step(A,new,J) == 0
        currents.append(J.copy());A=new
    if not baseline:
        reject_step(A,A,13,'step-exhausted')
    np.save(root/f'currents-{comm.rank}.npy',currents)
    assert lib.sp_tddft_close() == 0
    if not baseline:
        expect(lib.sp_tddft_step(A,A,J),12,'step-after-close')
    native.close()
    if comm.rank == 0:
        (root/'validation.json').write_text(json.dumps({'status':'passed','ranks':2,'cases':rejected,'steps':8,'baseline':baseline},indent=2)+'\n')
        print('native TDDFT admission PASS',rejected)


if __name__ == '__main__':
    main()
