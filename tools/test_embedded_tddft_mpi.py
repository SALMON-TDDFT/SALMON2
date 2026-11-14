"""Five ranks: field + two native two-rank TDDFT groups, one point each.
Args: EMBEDDED_BUILD_JSON SERIAL_RUNTIME GS_BUNDLE SHARED_NEW_OUTPUT
"""
import json,sys
from pathlib import Path
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import ParallelPlan,SolverGroup,ParallelContext
from salmonpy.material_response import MaterialBinding,ResponseSession
from salmonpy.embedded_tddft import NativeTDDFTMaterialProvider,NativeTDDFTBackend
from salmonpy.runtime import load_runtime
from salmonpy.tddft import TDDFTInitialState

world=MPI.COMM_WORLD
if world.size not in (4,5):raise ValueError('four or five ranks')
record,runtime,bundle,out=sys.argv[1:5];out=Path(out)
rgrid='--rgrid' in sys.argv[5:]
kgrid='--kgrid' in sys.argv[5:]
try:
    if world.rank==0:out.mkdir(parents=True,exist_ok=False)
    world.Barrier()
    drives=np.array([[[step*1e-4,0.,0.],[0.,step*2e-4,0.]] for step in range(1,5)])
    expected=None
    if world.rank==0:
        exe=load_runtime(runtime,required=['tddft','single-rank'])
        expected=np.empty_like(drives)
        with TDDFTInitialState.from_salmon(bundle,'rt.inp',executable=exe,max_steps=4,threads=1) as initial:
            for p in range(2):
                with initial.new_state(dt=.02) as state:
                    for s in range(4):expected[s,p]=state.step(drives[s,p])  # TDDFTState returns charge current j
    plan=ParallelPlan((SolverGroup('field',(0,),role='field'),
                       SolverGroup('left',tuple(range(1,world.size-2)),point_ids=(1,)),SolverGroup('right',(world.size-2,world.size-1),point_ids=(2,))))
    provider=NativeTDDFTMaterialProvider(bundle,'rt.inp',max_steps=4)
    bindings={}
    for point,size in ((1,world.size-3),(2,2)):
        decomp=(size,1,1,1,1) if kgrid else ((1,1,size,1,1) if rgrid else (1,size,1,1,1))
        backend=NativeTDDFTBackend(record,out,decomposition=decomp)
        bindings[point]=MaterialBinding('point-'+str(point),provider,backend)
    with ParallelContext.create(world,plan) as context:
        session=ResponseSession(context,bindings,dt=.02)
        actual=[]
        for state in session.states.values():
            try:state.advance(np.ones(3),np.zeros(3))
            except ValueError as exc:assert 'stale A' in str(exc)
            else:raise AssertionError('stale A must be rejected before native advance')
        try:
            for A in drives:
                J=session.advance(A if world.rank==0 else None)
                if world.rank==0:actual.append(J.tolist())
        finally:session.close()
        try:
            ResponseSession(context,bindings,dt=.02)
        except ValueError as exc:
            assert 'one native electronic lifetime' in str(exc)
        else:raise AssertionError('reopening native module state must be rejected')
    # This would fail/hang if native code finalized MPI or used world collectives.
    assert not MPI.Is_finalized()
    assert world.allreduce(1)==world.size
    if world.rank==0:
        actual=np.array(actual)
        np.testing.assert_allclose(actual,expected,rtol=2e-10,atol=2e-12)
        (out/'comparison.json').write_text(json.dumps({'max_abs_error':float(np.max(np.abs(actual-expected))),
             'actual':actual.tolist(),'reference':expected.tolist(),'mpi_alive_after_close':True},indent=2))
        print('embedded two-group TDDFT PASS',float(np.max(np.abs(actual-expected))))
except BaseException:
    import traceback
    traceback.print_exc();world.Abort(1)
