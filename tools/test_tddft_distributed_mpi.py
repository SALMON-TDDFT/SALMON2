"""mpiexec -n 3 python tools/test_tddft_distributed_mpi.py SERIAL_RUNTIME BUNDLE
Four independent points in two one-rank material groups, native TDDFT workers.
Synthetic GS fixture from test_tddft_entrypoint.py; not a converged calculation.
"""
from dataclasses import replace
from salmonpy.fdtd import SalmonExecutable
import sys
import tempfile,shutil
from pathlib import Path
import numpy as np
from mpi4py import MPI
from salmonpy import material_response as m
assert hasattr(m,'TDDFTMaterialProvider'), 'TDDFT material provider is missing'
from salmonpy.distributed import ParallelPlan,SolverGroup,ParallelContext
from salmonpy.runtime import load_runtime
from salmonpy.tddft import TDDFTInitialState
assert MPI.COMM_WORLD.size==3
# A changed bundle is rejected before a worker is started.
with tempfile.TemporaryDirectory() as tmp:
    bundle=Path(tmp)/'bundle';shutil.copytree(sys.argv[2],bundle)
    changed=m.TDDFTMaterialProvider(bundle,'rt.inp',runtime=sys.argv[1],max_steps=8,threads=1)
    inp=bundle/'rt.inp';inp.write_text(inp.read_text()+'\n! changed after identity\n')
    try: changed.prepare(None)
    except ValueError as exc: assert 'changed' in str(exc)
    else: raise AssertionError('changed initial state accepted')
    changed.close()
provider=m.TDDFTMaterialProvider(sys.argv[2],'rt.inp',runtime=sys.argv[1],max_steps=8,threads=1)
plan=ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('a',(1,),point_ids=(1,2)),SolverGroup('b',(2,),point_ids=(3,4))))
binding=m.MaterialBinding('crystal',provider,m.TDDFTResponseBackend())
bad_plan=ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('bad',(1,2),point_ids=(1,2))))
with ParallelContext.create(MPI.COMM_WORLD,bad_plan) as context:
    try: m.ResponseSession(context,{p:binding for p in bad_plan.point_ids},dt=.02)
    except ValueError as exc: assert 'one rank' in str(exc)
    else: raise AssertionError('serial TDDFT accepted a multi-rank material group')
with ParallelContext.create(MPI.COMM_WORLD,plan) as context:
    with TDDFTInitialState.from_salmon(sys.argv[2],'rt.inp',executable=load_runtime(sys.argv[1],required=['tddft','single-rank']),max_steps=8,threads=1) as baseline:
        try: m.TDDFTResponseBackend().create_state(baseline,1,.02,None)
        except TypeError: pass
        else: raise AssertionError('unverified initial state accepted by serial backend')
        # The backend must launch the verified bundle, not a same-hash path
        # supplied by a custom loader (which could resolve different libraries).
        with tempfile.TemporaryDirectory() as tmp:
            relocated=replace(baseline,executable=SalmonExecutable(Path(tmp)/'not-present',baseline.executable.sha256))
            guarded=m.TDDFTResponseBackend().create_state(m.TDDFTPreparedState(relocated,sys.argv[1]),1,.02,None)
            guarded.close()
        references=[baseline.new_state(dt=.02) for _ in plan.point_ids] if context.group.role=='field' else []
        try:
            histories=[]
            for repeat in range(2):
                session=m.ResponseSession(context,{p:binding for p in plan.point_ids},dt=.02)
                if repeat==1: provider.close()  # Live point states own their copies.
                for step in range(1,4):
                    A=np.array([[p*step*1e-4,0.,0.] for p in plan.point_ids])
                    actual=session.advance(A if context.group.role=='field' else None)
                    if actual is not None:
                        if repeat==0:
                            expected=np.array([state.step(a) for state,a in zip(references,A)])  # charge current j
                            histories.append(expected.copy())
                        else: expected=histories[step-1]
                        np.testing.assert_array_equal(actual,expected)
                session.close()
        finally:
            for state in references:state.close()
provider.close()
if MPI.COMM_WORLD.rank==0: print('four-point TDDFT charge J, independent repeats and serial equality PASS')
