"""Manual six-rank owner-route contract; no solver, native build or GS required.
mpiexec -n 6 python tools/test_owner_routing_mpi.py NEW_REPORT.json
The reporting gather after context close is test-only, outside step traffic.
"""
import json,sys
from pathlib import Path
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import ParallelContext,ParallelPlan,SolverGroup
from salmonpy.material_response import MaterialBinding,ResponseSession
from salmonpy.field_input import FieldContract,FieldValue,FieldPacket
from salmonpy.timing import TimeSpec,InputPort,OutputPort

class Provider:
    identity='routing-initial/1'
    def prepare(self,context):return None
class Backend:
    identity='routing-custom-current/1';execution='collective';reduction='sum'
    def create_state(self,initial,point,dt,context):
        class State:
            current=np.zeros(3)
            def __init__(self):self.steps=0
            def step(self,A):
                self.steps+=1
                self.current=((point+.5)*A+self.steps*.01)/context.solver_comm.size
                return self.current
            def step_fields(self,packet):return self.step(packet.A.value)
            def close(self):pass
        return State()
class ObservedComm:
    def __init__(self,comm,scope,rows):self.comm=comm;self.scope=scope;self.rows=rows
    def __getattr__(self,name):return getattr(self.comm,name)
    def Bcast(self,value,*args,**kw):
        result=self.comm.Bcast(value,*args,**kw)
        self.rows.append((self.scope,'Bcast',list(value.shape)));return result
    def bcast(self,value,*args,**kw):
        result=self.comm.bcast(value,*args,**kw)
        if isinstance(result,dict):self.rows.append((self.scope,'bcast-dict',sorted(result)))
        return result
    def Scatterv(self,send,recv,*args,**kw):
        result=self.comm.Scatterv(send,recv,*args,**kw)
        self.rows.append((self.scope,'Scatterv',list(recv.shape)));return result
    def scatter(self,send,*args,**kw):
        result=self.comm.scatter(send,*args,**kw)
        self.rows.append((self.scope,'scatter-dict',sorted(result)));return result
world=MPI.COMM_WORLD
if world.size!=6:raise ValueError('launch six ranks')
plan=ParallelPlan((SolverGroup('left',(0,1),point_ids=(51,7)),SolverGroup('field',(2,3),role='field'),
                   SolverGroup('right',(4,5),point_ids=(99,11,30))))
violations=[];records=[];dt=.01
for packet in (False,True):
 with ParallelContext.create(world,plan) as context:
    traffic=[]
    context.solver_comm=ObservedComm(context.solver_comm,'subgroup',traffic)
    if context.is_leader:context.leaders=ObservedComm(context.leaders,'leaders',traffic)
    contracts={p:FieldContract(TimeSpec(dt,(InputPort('A'),),(OutputPort('j'),)),{'A':'A/c:au'},(p,0,0),'velocity') for p in plan.point_ids}
    binding=MaterialBinding('custom',Provider(),Backend())
    session=ResponseSession(context,{p:binding for p in plan.point_ids},dt=dt,field_contracts=contracts if packet else None)
    traffic.clear()
    for step in range(1,5):
        values=np.array([[p*step*.001,-p*.002,step*.003] for p in plan.point_ids])
        if packet:
            inputs={p:FieldPacket(A=FieldValue(v,time=step*dt,position=(p,0,0),unit='A/c:au')) for p,v in zip(plan.point_ids,values)}
            current=session.step_fields(inputs if context.group.role=='field' and context.is_leader else None)
        else:current=session.step(values if context.group.role=='field' and context.is_leader else None)
        if current is not None:
            expected=np.array([(p+.5)*v+step*.01 for p,v in zip(plan.point_ids,values)])
            if not np.array_equal(current,expected):violations.append('current/control mismatch')
    for scope,method,payload in traffic:
        if method in ('Bcast','Scatterv'):
            if payload!=[len(context.group.point_ids),3]:violations.append((scope,method,payload))
            if scope=='leaders' and method=='Bcast':violations.append('full leader Bcast remains')
        else:
            if payload!=list(context.group.point_ids):violations.append((scope,method,payload))
    records.append(dict(rank=world.rank,group=context.group.name,packet=packet,owned=list(context.group.point_ids),traffic=traffic))
    session.close()
all_results=world.gather(dict(violations=violations,records=records),root=0)
if world.rank==0:
    report=dict(scope='synthetic owner routing; no material physics or solver',ranks=6,point_ids=list(plan.point_ids),results=all_results)
    with Path(sys.argv[1]).open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    if any(row['violations'] for row in all_results):raise AssertionError([row['violations'] for row in all_results])
    print('PASS owner-only A and packet payloads, independent currents and host MPI')
