"""Fast-PR HSE deep-error/retry test. BUILD GS NEW_OUTPUT [--baseline] [--faults].

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
from salmonpy.fast_hse_mpi import NativeFastHSEBackend
from salmonpy.fast_hse import FastHSEMaterialProvider
from salmonpy.errors import NumericalError


def main():
    build,gs,destination=sys.argv[1:4]
    baseline='--baseline' in sys.argv[4:];faults='--faults' in sys.argv[4:]
    fault_step=int(sys.argv[sys.argv.index('--fault-step')+1]) if '--fault-step' in sys.argv else 4
    overflow_A=float(sys.argv[sys.argv.index('--overflow-A')+1]) if '--overflow-A' in sys.argv else 1e200
    world=MPI.COMM_WORLD
    symmetry='--symmetry' in sys.argv[4:]
    axis=2 if symmetry else 0
    split='--split' in sys.argv[4:]
    assert world.size==(4 if split else 2)
    group=world.rank//2 if split else 0
    steps=int(sys.argv[sys.argv.index('--steps')+1]) if '--steps' in sys.argv else 16
    comm=world.Split(group,world.rank) if split else world
    partition='k-only'
    baseline=baseline or (split and group==1)
    lib=c.CDLL(json.loads(Path(build).read_text())['library'])
    if not baseline:
        lib.sp_fast_hse_recovery_protocol.restype=c.c_int
        assert lib.sp_fast_hse_recovery_protocol()==1
        assert hasattr(lib,'sp_fast_hse_test_fault')==faults
    out=Path(destination).resolve()
    if split:out=out/('group-'+str(group))
    if comm.rank==0:out.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    ctx=SimpleNamespace(solver_comm=comm,point_ids=(1,))
    initial=FastHSEMaterialProvider(gs,'rt.inp',max_steps=steps).prepare(ctx)
    state=NativeFastHSEBackend(build,out,allow_symmetry=symmetry).create_state(initial,1,.001,ctx)
    checkpoint=state.lib.sp_fast_hse_checkpoint;checkpoint.argtypes=[c.c_char_p,c.c_int];checkpoint.restype=c.c_int
    def save(tag):
        p=out/f'{tag}-{comm.rank}.bin'
        assert state._call(checkpoint,os.fsencode(p),0)==0
        return p.read_bytes()
    def trace_files():
        return {str(p.relative_to(state.work)):p.read_bytes() for p in state.work.glob('*.data')}
    rows=[];cases=[]
    ptcn='hse_ptcn' in (Path(gs)/'rt.inp').read_text().lower()
    for step in range(1,steps+1):
        value=np.zeros(3);value[axis]=.001*np.sin(.2*step)
        if not baseline and step==fault_step:
            attempts=[('overflow',None)]
            if faults and ptcn:attempts.append(('ptcn-solver',32))
            if faults:attempts += [('stage-'+str(stage),stage) for stage in (20,21,22,23,24,30,33)]
            for tag,stage in attempts:
                before=save(tag+'-before');files=trace_files();A=state.A.copy();index=state.steps;old_current=state.current.copy();old_energy=state.energy
                if stage is not None:
                    fn=state.lib.sp_fast_hse_test_fault;fn.argtypes=[c.c_int,c.c_int];fn.restype=c.c_int
                    assert fn(stage,1)==0
                drive=value.copy()
                if stage is None:drive[axis]=overflow_A
                try:state.advance(A,drive)
                except NumericalError:pass
                else:raise AssertionError('deep error was not returned')
                assert state.steps==index
                np.testing.assert_array_equal(state.A,A)
                assert not state.failed
                np.testing.assert_array_equal(state.current,old_current)
                assert state.energy==old_energy
                assert save(tag+'-after')==before,tag
                assert trace_files()==files,tag+' wrote failed output'
                cases.append(tag)
        rows.append(np.r_[state.advance(state.A.copy(),value),state.energy])
    np.save(out/f'current-{comm.rank}.npy',rows)
    if '--capacity-probe' in sys.argv:
        # Bypass Python's capacity guard to verify the raw ABI never reports
        # success without propagation. This deliberately ends the native lifetime.
        current=np.empty(3);energy=c.c_double()
        code=state._call(state.lib.sp_fast_hse_step,state.A.copy(),np.array([.123,0.,0.]),current,c.byref(energy))
        assert code!=0, 'exhausted raw ABI falsely reported success'
        cases.append('raw-capacity-rejection')
    state.close()
    if comm.rank==0:
        traces=out/'accepted-traces';traces.mkdir()
        import shutil
        for path in state.work.glob('*.data'):shutil.copyfile(path,traces/path.name)
    if comm.rank==0:
        (out/'validation.json').write_text(json.dumps({'status':'passed','cases':cases,'steps':steps,'ranks':2,'replayed_steps':0,'partition':partition},indent=2)+'\n')
        print('deep recovery PASS',cases)


if __name__=='__main__':main()
