"""Real deep-error/retry test. BUILD GS NEW_OUTPUT [--baseline] [--faults].

Two ranks; OMP/BLAS1 and launcher timeout required. --faults requires a separately
compiled SALMONPY_RECOVERY_TESTING build; ordinary builds must not export its hook.
"""
import ctypes as c
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.embedded_tddft import NativeTDDFTBackend, NativeTDDFTMaterialProvider
from salmonpy.errors import NumericalError


def main():
    build,gs,destination=sys.argv[1:4]
    baseline='--baseline' in sys.argv[4:];faults='--faults' in sys.argv[4:]
    fault_step=int(sys.argv[sys.argv.index('--fault-step')+1]) if '--fault-step' in sys.argv else 4
    world=MPI.COMM_WORLD
    split='--split' in sys.argv[4:]
    assert world.size==(4 if split else 2)
    group=world.rank//2 if split else 0
    comm=world.Split(group,world.rank) if split else world
    partition=sys.argv[sys.argv.index('--partition')+1] if '--partition' in sys.argv else 'orbital'
    decomposition={'orbital':(1,2,1,1,1),'k':(2,1,1,1,1),'rgrid':(1,1,2,1,1)}[partition]
    baseline=baseline or (split and group==1)
    lib=c.CDLL(json.loads(Path(build).read_text())['library'])
    if not baseline:
        lib.sp_tddft_recovery_protocol.restype=c.c_int
        assert lib.sp_tddft_recovery_protocol()==1
        assert hasattr(lib,'sp_tddft_test_fault')==faults
    out=Path(destination).resolve()
    if split:out=out/('group-'+str(group))
    if comm.rank==0:out.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    ctx=SimpleNamespace(solver_comm=comm,point_ids=(1,))
    initial=NativeTDDFTMaterialProvider(gs,'rt.inp',max_steps=16).prepare(ctx)
    state=NativeTDDFTBackend(build,out,decomposition=decomposition).create_state(initial,1,.003,ctx)
    checkpoint=state.lib.sp_tddft_checkpoint;checkpoint.argtypes=[c.c_char_p,c.c_int];checkpoint.restype=c.c_int
    def save(tag):
        p=out/f'{tag}-{comm.rank}.bin'
        assert state._call(checkpoint,os.fsencode(p),0)==0
        return p.read_bytes()
    def trace_files():
        return {str(p.relative_to(state.work)):p.read_bytes() for p in state.work.glob('*.data')}
    rows=[];cases=[]
    for step in range(1,17):
        value=np.array([.001*np.sin(.2*step),0.,0.])
        if not baseline and step==fault_step:
            attempts=[('overflow',None)]
            if faults:attempts += [('stage-'+str(stage),stage) for stage in (20,21,22,23,24)]
            for tag,stage in attempts:
                before=save(tag+'-before');files=trace_files();A=state.A.copy();index=state.steps
                if stage is not None:
                    fn=state.lib.sp_tddft_test_fault;fn.argtypes=[c.c_int,c.c_int];fn.restype=c.c_int
                    assert fn(stage,1)==0
                drive=np.array([1e100,0.,0.]) if stage is None else value
                try:state.advance(A,drive)
                except NumericalError:pass
                else:raise AssertionError('deep error was not returned')
                assert state.steps==index
                np.testing.assert_array_equal(state.A,A)
                assert save(tag+'-after')==before,tag
                assert trace_files()==files,tag+' wrote failed output'
                cases.append(tag)
        rows.append(state.advance(state.A.copy(),value))
    np.save(out/f'current-{comm.rank}.npy',rows)
    state.close()
    if comm.rank==0:
        traces=out/'accepted-traces';traces.mkdir()
        import shutil
        for path in state.work.glob('*.data'):shutil.copyfile(path,traces/path.name)
    if comm.rank==0:
        (out/'validation.json').write_text(json.dumps({'status':'passed','cases':cases,'steps':16,'ranks':2,'replayed_steps':0,'partition':partition},indent=2)+'\n')
        print('deep recovery PASS',cases)


if __name__=='__main__':main()
