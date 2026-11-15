"""Three-rank sample regression plus collective callback-failure handling."""
import sys
from mpi4py import MPI
from salmonpy.examples import mpi_field_inputs as demo
from salmonpy.errors import InputError

mode=sys.argv[1]
sys.argv=[sys.argv[0]]
if mode=='callback-failure':
    original=demo.Backend.create_state
    def create(self,*args,**kwargs):
        state=original(self,*args,**kwargs)
        if MPI.COMM_WORLD.rank==2:
            def fail(fields): raise RuntimeError('injected E/H failure')
            state.step_fields=fail
        return state
    demo.Backend.create_state=create
elif mode=='output-mismatch':
    original=demo.Backend.create_state
    def create(self,initial,point_id,dt,context):
        state=original(self,initial,point_id,dt,context)
        ports=(demo.InputPort('A'),) if point_id==10 else (demo.InputPort('E'),demo.InputPort('H'))
        state.timing=demo.TimeSpec(dt,ports,(demo.OutputPort('J_matter'),))
        return state
    demo.Backend.create_state=create
elif mode!='success':raise ValueError('unknown scenario')
failed=False
try:demo.main()
except InputError as exc:
    expected={'callback-failure':'injected E/H failure','output-mismatch':'backend output timing'}
    if mode not in expected or expected[mode] not in str(exc): raise
    failed=True
assert MPI.COMM_WORLD.allreduce(int(failed))==(0 if mode=='success' else 3)
if MPI.COMM_WORLD.rank==0:print('PASS',mode)
