"""Count outer agreements while running an explicitly selected MPI driver.

Default driver is the ordinary mpi_mixed_native sample. --driver SCRIPT allows
other standalone MPI contract drivers. Counts exclude setup/teardown and native
subgroup communication. All assertions run after the driver closes its contexts.
"""
import argparse
import json
from pathlib import Path
import runpy
import sys
from unittest.mock import patch
from mpi4py import MPI

def main():
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--expected-votes',type=int,required=True)
    parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--expected-leader-ops',type=int)
    parser.add_argument('--driver',type=Path)
    parser.add_argument('--failures-only',action='store_true')
    args,remaining=parser.parse_known_args()
    from salmonpy.material_response import ResponseSession
    from salmonpy.coupled import CoupledStepper,PacketCoupledStepper
    original=ResponseSession.agree
    counts=[0];steps=[];failed_steps=[]
    class Observed:
        def __init__(self,comm,events,parent=False):self.comm=comm;self.events=events;self.parent=parent
        def __getattr__(self,name):
            value=getattr(self.comm,name)
            if name not in ('allgather','bcast','scatter','alltoall','Scatterv','Gatherv','Alltoallv'):return value
            def call(*a,**kw):
                self.events.append({'operation':name,**({'root':kw.get('root',0)} if self.parent else {})})
                return value(*a,**kw)
            return call
    def agree(self,error):
        counts[0]+=1
        return original(self,error)
    def wrap(method,kind):
        def step(self):
            before=counts[0];context=self.response.context;leaders=context.leaders
            parent=getattr(context,'_error_comm',None);events=[];broadcasts=[];success=False
            if context.is_leader:context.leaders=Observed(leaders,events)
            if parent is not None:context._error_comm=Observed(parent,broadcasts,parent=True)
            try:
                result=method(self);success=True
                return result
            finally:
                context.leaders=leaders
                if parent is not None:context._error_comm=parent
                record={'kind':kind,'step':self.step_index,'agreements':counts[0]-before,
                        'is_leader':context.is_leader,'leader_operations':events,'failure_broadcasts':broadcasts}
                (steps if success else failed_steps).append(record)
        return step
    sys.argv=[str(args.driver) if args.driver else sys.argv[0],*remaining]
    with patch.object(ResponseSession,'agree',agree),patch.object(ResponseSession,'_agree',agree), \
         patch.object(CoupledStepper,'step',wrap(CoupledStepper.step,'A')), \
         patch.object(PacketCoupledStepper,'step',wrap(PacketCoupledStepper.step,'packets')):
        if args.driver:runpy.run_path(str(args.driver),run_name='__main__')
        else:
            from salmonpy.examples.mpi_mixed_native import main as sample
            sample()
    rows=MPI.COMM_WORLD.gather({'rank':MPI.COMM_WORLD.rank,'steps':steps,'failed_steps':failed_steps},root=0)
    failure=None
    if MPI.COMM_WORLD.rank==0:
        valid=all((row['failed_steps'] and not row['steps']) if args.failures_only else
                  (row['steps'] and all(s['agreements']==args.expected_votes for s in row['steps'])) for row in rows)
        if args.expected_leader_ops is not None:
            valid=valid and all(len(s['leader_operations'])==args.expected_leader_ops
                for row in rows for s in row['steps'] if s['is_leader'])
        valid=valid and all(not s['failure_broadcasts'] for row in rows for s in row['steps'])
        args.report.write_text(json.dumps({'passed':valid,'expected':args.expected_votes,'ranks':rows},indent=2)+'\n')
        if not valid:failure='outer agreement count differs from the requested contract'
    failure=MPI.COMM_WORLD.bcast(failure,root=0)
    if failure:raise AssertionError(failure)
    if MPI.COMM_WORLD.rank==0:print('PASS outer agreement count',args.expected_votes)

if __name__=='__main__':main()
