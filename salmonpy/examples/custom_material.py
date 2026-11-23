"""One reusable synthetic material, serial or point-distributed MPI.

python -m salmonpy.examples.custom_material
mpiexec -n 3 python -m salmonpy.examples.custom_material --mpi NEW_JSON
This memory law demonstrates the interface; it is not a physical material fit.
"""
from dataclasses import dataclass,asdict
import json
import numpy as np
from salmonpy.fields import FieldSamples
from salmonpy.errors import InputError,NumericalError
from salmonpy.timing import TimeSpec,InputPort,OutputPort


@dataclass(frozen=True)
class ExampleMaterial:
    gain: float = 2.
    decay: float = .8
    backend: str = 'example-memory-material/1'

    def __post_init__(self):
        if not np.isfinite(self.gain) or not np.isfinite(self.decay):raise InputError('finite coefficients required')

    @property
    def identity(self):
        return json.dumps(asdict(self),sort_keys=True,allow_nan=False)

    def new_state(self, *, dt, time0=0., initial_A=(0.,0.,0.)):
        return _State(self,dt,time0,initial_A)


class _State:
    def __init__(self,model,dt,time0,A):
        self.model=model;self.dt=dt;self.time0=time0;self.time=time0;self.n=0
        self.current=np.zeros(3);self.A=np.array(A,dtype=float,copy=True);self.closed=False
        self.timing=TimeSpec(dt,(InputPort('A'),),(OutputPort('j'),),t0=time0,
                             history='previous A and j',integrator='synthetic memory recurrence')

    def step(self,A):
        if self.closed:raise InputError('state is closed')
        new=np.array(A,dtype=float,copy=True)
        if new.shape!=(3,) or not np.isfinite(new).all():raise InputError('finite three-vector A required')
        with np.errstate(over='raise',invalid='raise'):
            current=self.model.decay*self.current+self.model.gain*(new-self.A)
        if not np.isfinite(current).all():raise NumericalError('nonfinite material current')
        # Every check precedes accepted-state mutation.
        self.current=current;self.A=new;self.n+=1;self.time=self.time0+self.n*self.dt
        return self.current.copy()

    def close(self):self.closed=True


def waveform(point=10):
    n=np.arange(9);values=np.column_stack((n,np.sin(.3*n),np.cos(.2*n)-1))*point*1e-4
    return FieldSamples(n*.02,values,'A')


def bindings():
    from salmonpy.materials import bind_material
    a=ExampleMaterial(gain=1.);b=ExampleMaterial(gain=2.)
    one=bind_material(a,key='one',identity=a.identity)
    two=bind_material(b,key='two',identity=b.identity)
    return {10:one,20:two,30:one,40:two}


def plan_for(size):
    from salmonpy.distributed import ParallelPlan,SolverGroup
    if size not in (2,3):raise InputError('example supports2 or3 MPI ranks')
    points=(10,20,30,40)
    return ParallelPlan((SolverGroup('field',(0,),role='field'),)+tuple(
        SolverGroup('material'+str(rank),(rank,),point_ids=points[rank-1::size-1])
        for rank in range(1,size)))


def run_mpi(world):
    from salmonpy.distributed import ParallelContext
    from salmonpy.material_response import ResponseSession
    plan=plan_for(world.size);history=[]
    with ParallelContext.create(world,plan) as context:
        session=ResponseSession(context,bindings(),dt=.02)
        try:
            timing=session.timing_records()
            initial=session.initial_current
            for n in range(1,9):
                A=np.array([waveform(p).values[n] for p in plan.point_ids]) if world.rank==0 else None
                current=session.advance(A)
                if world.rank==0:history.append(current.tolist())
            communication=session.communication_report()
        finally:session.close()
    return {'point_ids':plan.point_ids,'history':history,'initial_current':initial,
            'timing':timing,'communication':communication,
            'scope':'synthetic serial model distributed across points; not internal MPI parallelism'} if world.rank==0 else None


def main():
    import argparse
    from salmonpy.materials import check_material
    from salmonpy.crystal import respond
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mpi',action='store_true');parser.add_argument('output',nargs='?')
    args=parser.parse_args()
    if args.mpi:
        from mpi4py import MPI
        result=run_mpi(MPI.COMM_WORLD)
        if MPI.COMM_WORLD.rank:return
    else:
        model=ExampleMaterial();drive=waveform()
        result={'check':asdict(check_material(model,drive)),
                'current':respond(model,drive).current.tolist(),'scope':'synthetic demonstration only'}
    if args.output:
        from pathlib import Path
        with Path(args.output).open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    else:print(json.dumps(result,indent=2,allow_nan=False))

if __name__=='__main__':main()
