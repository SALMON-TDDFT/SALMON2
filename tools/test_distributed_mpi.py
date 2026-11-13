"""Run with mpiexec -n 6; tests groups, user response and native subgroup bridge.
Usage: python tools/test_distributed_mpi.py BRIDGE_LIBRARY SHA256
No full FDTD/TDDFT coupling is claimed by this manufactured response test.
"""
import sys
import numpy as np
from mpi4py import MPI
from salmonpy import distributed as d
assert hasattr(d,'NativeParallelContext'), 'Fortran communicator attachment missing'
from salmonpy.material_response import MaterialBinding, ResponseSession

class Provider:
    identity='test-initial/1'
    def prepare(self,context): return 2.
class Backend:
    identity='test-memory/1'; execution='collective'; reduction='sum'
    def create_state(self,initial,point_id,dt,context):
        class State:
            def __init__(self): self.j=np.zeros(3)
            def advance(self,old,new):
                self.j=self.j+initial*(new-old)+dt*point_id
                return self.j/context.solver_comm.size
            def close(self): pass
        return State()

world=MPI.COMM_WORLD
assert world.size==6
plan=d.ParallelPlan((d.SolverGroup('field',(0,),role='field'),
                    d.SolverGroup('one',(1,2),point_ids=(10,11)),
                    d.SolverGroup('two',(3,4,5),point_ids=(20,21))))
# An inconsistent plan must fail collectively before Split.
try:
    d.ParallelContext.create(world,plan if world.rank else d.ParallelPlan(plan.groups[::-1]))
except ValueError: pass
else: raise AssertionError('inconsistent plans accepted')
with d.ParallelContext.create(world,plan) as ctx:
    native=d.NativeParallelContext(ctx,sys.argv[1],sys.argv[2])
    try: ctx.close()
    except RuntimeError: pass
    else: raise AssertionError('parent communicator freed while native context is attached')
    try: d.NativeParallelContext(ctx,sys.argv[1],sys.argv[2])
    except RuntimeError: pass
    else: raise AssertionError('duplicate native attachment accepted')
    assert native.info()==(ctx.solver_comm.rank,ctx.solver_comm.size)
    assert native.sum(float(ctx.solver_comm.rank))==ctx.solver_comm.size*(ctx.solver_comm.size-1)/2
    for _ in range(ctx.solver_comm.size): native.barrier()
    native.close()
    assert not MPI.Is_finalized()
    # A field group may have no native attachment while material groups do.
    material_native=d.NativeParallelContext(ctx,sys.argv[1],sys.argv[2]) if ctx.group.role=='material' else None
    try: ctx.close()
    except RuntimeError: pass
    else: raise AssertionError('one group freed shared leaders while another was busy')
    if material_native is not None: material_native.close()
    # Different states from the same prepared material, independent per point.
    binding=MaterialBinding('custom',Provider(),Backend())
    session=ResponseSession(ctx,{p:binding for p in plan.point_ids},dt=.1)
    for step in range(1,5):
        values=np.array([[step*.01*p,0.,0.] for p in plan.point_ids])
        result=session.advance(values if ctx.group.role=='field' and ctx.is_leader else None)
        if ctx.group.role=='field':
            expected=np.array([[2*step*.01*p+step*.1*p,step*.1*p,step*.1*p] for p in plan.point_ids])
            np.testing.assert_allclose(result,expected,rtol=1e-14,atol=1e-14)
    # Invalid drive is rejected by every group; no one hangs in a collective.
    try: session.advance(np.full((4,3),np.nan) if ctx.group.role=='field' else None)
    except ValueError: pass
    else: raise AssertionError('invalid drive accepted')
    session.close()
    collision={p:binding for p in plan.point_ids}
    collision[11]=MaterialBinding('custom',Provider(),Backend())
    try: ResponseSession(ctx,collision,dt=.1)
    except ValueError: pass
    else: raise AssertionError('ambiguous material key accepted')
    class Strided(Backend):
        identity='test-strided/1'
        def create_state(self,initial,point_id,dt,context):
            class State:
                def advance(self,old,new): return np.arange(6.,dtype=float)[::2]/context.solver_comm.size
                def close(self): pass
            return State()
    binding=MaterialBinding('strided',Provider(),Strided())
    session=ResponseSession(ctx,{p:binding for p in plan.point_ids},dt=.1)
    result=session.advance(np.zeros((4,3)) if ctx.group.role=='field' else None)
    if ctx.group.role=='field': np.testing.assert_allclose(result,np.tile([0.,2.,4.],(4,1)))
    session.close()
    closed=[]
    class Failing(Backend):
        identity='test-factory-error/1'
        def create_state(self,initial,point_id,dt,context):
            if point_id==11: raise ValueError('factory failure')
            class State:
                def close(self): closed.append(point_id)
            return State()
    binding=MaterialBinding('failing',Provider(),Failing())
    try: ResponseSession(ctx,{p:binding for p in plan.point_ids},dt=.1)
    except ValueError: pass
    else: raise AssertionError('factory error ignored')
    assert closed==[p for p in ctx.group.point_ids if p!=11], closed
    class Overflow(Backend):
        identity='test-overflow/1'
        def create_state(self,initial,point_id,dt,context):
            class State:
                def advance(self,old,new): return np.full(3,np.finfo(float).max)
                def close(self): pass
            return State()
    binding=MaterialBinding('overflow',Provider(),Overflow())
    session=ResponseSession(ctx,{p:binding for p in plan.point_ids},dt=.1)
    for attempt in range(2):
        try: session.advance(np.zeros((4,3)) if ctx.group.role=='field' else None)
        except ValueError: pass
        else: raise AssertionError('overflow/failed session accepted')
    assert session.failed and session.step_index==0
    session.close()


world.Barrier()
if world.rank==0: print('Python/Fortran groups, independent custom states and A/J exchange PASS')
