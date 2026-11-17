"""Native E/H group numerical retry: BUILD NEW_OUTPUT [--baseline|--faults].

Use 1/2 ranks, --split for two independent two-rank groups, OMP/BLAS1 and
external timeout. Fault hooks must be compiled separately, never in production.
"""
import ctypes as c
import json
from pathlib import Path
from types import SimpleNamespace
import sys
import numpy as np
from mpi4py import MPI
from salmonpy.coupled import NativeEHField
from salmonpy.fdtd import FDTD,PlaneWave
from salmonpy.fields import AcosPulse
from salmonpy.geometry import Grid,Geometry,Shape
from salmonpy.ld import LorentzDrude
from salmonpy.errors import NumericalError


def main():
    build,destination=sys.argv[1:3]
    split='--split' in sys.argv;faults='--faults' in sys.argv
    stride=int(sys.argv[sys.argv.index('--observation-stride')+1]) if '--observation-stride' in sys.argv else 1
    world=MPI.COMM_WORLD
    if split:assert world.size==4
    else:assert world.size in (1,2)
    group=world.rank//2 if split else 0
    comm=world.Split(group,world.rank) if split else world
    baseline='--baseline' in sys.argv or group==1
    out=Path(destination).resolve()
    if split:out=out/f'group-{group}'
    if comm.rank==0:out.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    ctx=SimpleNamespace(solver_comm=comm,is_leader=comm.rank==0,closed=False,_dependents=0,
        group=SimpleNamespace(role='field'),plan=SimpleNamespace(point_ids=(1,)))
    grid=Grid((40,40,40),(2.,2.,2.));dt=.003;steps=32
    ld=(LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.]),
        LorentzDrude(plasma_frequency=15.,strength=[.7,.3],resonance=[25.,35.],damping=[2.,3.]))
    geometry=Geometry(grid,(Shape('box',(2.,4.,4.),(7.,1.,1.),material=1),Shape('box',(2.,4.,4.),(9.,1.,1.),material=2)))
    calc=FDTD(geometry,ld,dt,steps,source=PlaneWave(AcosPulse(.001,20.,.15),0,-16.),observations=((1.,1.,1.),),observation_stride=stride)
    f=NativeEHField(ctx,calc,build,out/'field',points={1:(20,20,20)},grid=(1,comm.size,1))
    if not baseline:
        f.lib.sp_eh_recovery_protocol.restype=c.c_int
        assert f.lib.sp_eh_recovery_protocol()==1
        assert hasattr(f.lib,'sp_eh_test_fault')==faults
    cp=f.lib.sp_eh_checkpoint;cp.argtypes=[c.c_char_p,c.c_int];cp.restype=c.c_int
    def save(tag):
        p=out/f'{tag}-{comm.rank}.bin';assert f._call(cp,str(p).encode(),0)==0;return p.read_bytes()
    def files():return {p.name:p.read_bytes() for p in f.directory.glob('*.data')}
    cases=[];rows=[]
    if not baseline and comm.size>1:
        before=save('mixed-entry-before')
        E=np.empty_like(f.A);J=np.zeros_like(f.A)
        floats=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
        f.lib.sp_eh_step_checked.argtypes=[c.c_int,floats,floats,floats]
        f.lib.sp_eh_step_checked.restype=c.c_int
        if comm.rank==0:code=f._call(f.lib.sp_eh_step_checked,1,J,f.A,E)
        else:code=f._call(f.lib.sp_eh_step,1,J,E)
        assert code==11
        assert save('mixed-entry-after')==before
        cases.append(('admission','mixed-entry'))
    for step in range(1,steps+1):
        J=np.array([[0.,.002*np.sin(.2*step),0.]]) if step>1 else np.zeros((1,3))
        if not baseline and step in (1,12,32):
            attempts=list(range(20,24)) if faults else []
            for stage in attempts:
                A=f.A.copy();E=f.E.copy();before=save(f'{step}-{stage}-before');outputs=files()
                fn=f.lib.sp_eh_test_fault;fn.argtypes=[c.c_int,c.c_int];fn.restype=c.c_int
                assert fn(stage,comm.size-1)==0
                try:f.advance(J)
                except NumericalError:pass
                else:raise AssertionError('numerical error not returned')
                assert not f.failed and f.step_index==step-1
                np.testing.assert_array_equal(f.A,A);np.testing.assert_array_equal(f.E,E)
                assert save(f'{step}-{stage}-after')==before
                assert files()==outputs
                cases.append((step,stage))
            if step==12:
                # Real finite-input overflow of A-dt*sample(E), before ordinary output.
                A=f.A.copy();E=f.E.copy();f.A[:]=np.finfo(float).max
                before=save('integrated-A-before');outputs=files()
                try:f.advance(np.full((1,3),np.finfo(float).max))
                except NumericalError:pass
                else:raise AssertionError('integrated A overflow not returned')
                assert not f.failed and f.step_index==step-1
                assert save('integrated-A-after')==before and files()==outputs
                np.testing.assert_array_equal(f.E,E);f.A=A
                cases.append((step,'integrated-A-overflow'))
        rows.append(np.r_[f.advance(J).ravel(),f.E.ravel(),f.sample_magnetic_field().values.ravel()])
    if not baseline:
        before=save('capacity-before');E=np.empty_like(f.A)
        code=f._call(f.lib.sp_eh_step,1,np.zeros_like(f.A),E)
        assert code==13 and save('capacity-after')==before
        cases.append(('admission','capacity'))
    np.save(out/f'fields-{comm.rank}.npy',rows)
    f.close()
    if comm.rank==0:
        (out/'validation.json').write_text(json.dumps(dict(status='passed',cases=cases,ranks=comm.size,steps=steps,replayed_steps=0),indent=2)+'\n')
        print('E/H deep recovery PASS',cases)


if __name__=='__main__':main()
