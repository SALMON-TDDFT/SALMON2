"""Three ranks: distinct A and E/H contracts with explicit charge-current models.

mpiexec -n 3 python -m salmonpy.examples.mpi_field_inputs
Synthetic interface demonstration, not molecular TDDFT or physical radiation.
"""
import argparse
import numpy as np
from salmonpy.field_input import FieldValue, FieldPacket, FieldContract
from salmonpy.timing import TimeSpec, InputPort, OutputPort
from salmonpy.material_response import MaterialBinding, ResponseSession
from salmonpy.distributed import SolverGroup, ParallelPlan, ParallelContext
from salmonpy.errors import InputError


class Provider:
    identity='field-demo-provider/1'
    def prepare(self,context): return None


class Backend:
    identity='field-demo-current/1'
    execution='collective'
    reduction='leader'
    def create_state(self,initial,point_id,dt,context):
        class State:
            def __init__(self): self.calls=0
            def step_fields(self,packet):
                self.calls+=1
                if point_id==10:
                    assert packet.E is None and packet.H is None
                    return -.2*packet.A.value
                assert packet.A is None
                return .3*packet.E.value+.1*packet.H.value
            def close(self): pass
        return State()


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    from mpi4py import MPI
    world=MPI.COMM_WORLD
    if world.size!=3: raise ValueError('exactly three ranks required')
    plan=ParallelPlan((SolverGroup('field',(0,),role='field'),
                       SolverGroup('A-material',(1,),point_ids=(10,)),
                       SolverGroup('EH-material',(2,),point_ids=(20,))))
    dt=.02
    contracts={10:FieldContract(TimeSpec(dt,(InputPort('A'),),(OutputPort('j'),)),
                               {'A':'A/c:au'},(0,0,0),'velocity'),
               20:FieldContract(TimeSpec(dt,(InputPort('E'),InputPort('H')),(OutputPort('j'),)),
                               {'E':'E:au','H':'H:au'},(1,0,0),'custom')}
    binding=MaterialBinding('demo',Provider(),Backend())
    with ParallelContext.create(world,plan) as context:
        session=ResponseSession(context,{10:binding,20:binding},dt=dt,field_contracts=contracts)
        try:
            timings=session.timing_records()
            assert set(timings[10]['inputs'])=={'A'}
            assert set(timings[20]['inputs'])=={'E','H'}
            def packets(t,missing=False):
                return {10:FieldPacket(A=FieldValue([1,0,0],time=t,position=(0,0,0),unit='A/c:au')),
                        20:FieldPacket(E=FieldValue([0,2,0],time=t,position=(1,0,0),unit='E:au'),
                         H=None if missing else FieldValue([0,0,3],time=t,position=(1,0,0),unit='H:au'))}
            for bad in (packets(dt,True),packets(2*dt)):
                try:session.step_fields(bad if world.rank==0 else None)
                except InputError:pass
                else:raise AssertionError('invalid packet accepted')
                assert session.step_index==0
                assert all(s.calls==0 for s in session.states.values())
            for n in range(1,4):
                result=session.step_fields(packets(n*dt) if world.rank==0 else None)
                if world.rank==0:np.testing.assert_allclose(result,[[-.2,0,0],[0,.6,.3]])
            assert all(s.calls==3 for s in session.states.values())
        finally:session.close()
    if world.rank==0:print('PASS MPI A/E/H contracts, preflight rejection and explicit current response')


if __name__=='__main__':main()
