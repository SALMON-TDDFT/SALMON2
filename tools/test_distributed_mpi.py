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
    # A mismatch on one subgroup member must reject everyone before attachment.
    from unittest.mock import patch
    from types import SimpleNamespace
    class Function:
        def __init__(self,fn):self.fn=fn
        def __call__(self,*args):return self.fn(*args)
    for variant in ('missing','version','width','status',
                    'sp_attach','sp_detach','sp_info','sp_barrier','sp_sum','sp_library'):
        attached=[]
        def attach(handle):
            attached.append(handle)
            raise AssertionError('mismatched host must never attach')
        def descriptor(capacity,words):
            values=[1,32,64,64,32,64,1,1,1,3]
            if ctx.solver_comm.rank==ctx.solver_comm.size-1:
                if variant=='status':return 13
                if variant=='version':values[0]=2
                if variant=='width':values[2]=32
            for i,v in enumerate(values):words[i]=v
            return 0
        def library(buffer,capacity,active):
            buffer.value=MPI.Get_library_version().encode()
            active._obj.value=1
        fake=SimpleNamespace(sp_host_abi=Function(descriptor),sp_attach=Function(attach),
            sp_detach=Function(lambda:0),sp_info=Function(lambda *a:None),
            sp_barrier=Function(lambda:None),sp_sum=Function(lambda x:x),sp_library=Function(library))
        if variant=='missing' and ctx.solver_comm.rank==ctx.solver_comm.size-1:del fake.sp_host_abi
        if variant.startswith('sp_') and ctx.solver_comm.rank==ctx.solver_comm.size-1:delattr(fake,variant)
        with patch('ctypes.CDLL',return_value=fake):
            try:d.NativeParallelContext(ctx,sys.argv[1],sys.argv[2])
            except ValueError:pass
            else:raise AssertionError('invalid host ABI accepted')
        assert not attached and ctx._dependents==0
        assert MPI.Is_initialized() and not MPI.Is_finalized()
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
