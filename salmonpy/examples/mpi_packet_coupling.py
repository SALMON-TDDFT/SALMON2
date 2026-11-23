"""Explicit packet feedback with synthetic currents (not a real material).

3 ranks for Multiscale1D; 3/4 for native (1/2 field ranks + 2 material ranks).
"""
import argparse,json
from pathlib import Path
import numpy as np
from salmonpy.coupled import MultiscaleField,NativeEHField,CoupledStepper,PacketCoupledStepper
from salmonpy.distributed import ParallelPlan,ParallelContext,SolverGroup
from salmonpy.material_response import MaterialBinding,ResponseSession,packetize_A_binding,TDDFTMaterialProvider,TDDFTResponseBackend
from salmonpy.field_input import FieldContract
from salmonpy.timing import TimeSpec,InputPort,OutputPort
from salmonpy.maxwell import Multiscale1D
from salmonpy.fields import AcosPulse
from salmonpy.geometry import Grid,Geometry
from salmonpy.fdtd import FDTD,PlaneWave
from salmonpy.errors import InputError


class Provider:
    identity='packet-demo/1'
    def prepare(self,context):return None


class Backend:
    identity='packet-demo-current/1';execution='leader_only';reduction='leader'
    def create_state(self,initial,point_id,dt,context):
        class State:
            def step(self,A):return -.01*A
            def step_fields(self,packet):
                if packet.A is not None:return self.step(packet.A.value)
                return .01*packet.E.value+.001*packet.H.value
            def close(self):pass
        return State()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--field',choices=('multiscale','native'),default='multiscale')
    parser.add_argument('--build',type=Path)
    parser.add_argument('--ports',choices=('a','eh'),default='a')
    parser.add_argument('--adapt-a',action='store_true',help='explicitly wrap the legacy A backend')
    parser.add_argument('--runtime',type=Path,help='managed TDDFT runtime (requires --bundle)')
    parser.add_argument('--bundle',type=Path,help='TDDFT bundle containing rt.inp')
    parser.add_argument('--legacy',action='store_true',help='A-only reference path')
    parser.add_argument('--reject-phase',action='store_true',help='test clock refusal before update')
    args=parser.parse_args()
    if args.ports=='eh' and (args.field!='native' or args.legacy):parser.error('E/H requires native packet path')
    if (args.field=='native')!=(args.build is not None):parser.error('--build required only for native')
    if args.adapt_a and (args.ports!='a' or args.legacy):parser.error('--adapt-a requires A packet mode')
    if (args.runtime is None)!=(args.bundle is None):parser.error('--runtime and --bundle required together')
    if args.runtime and not (args.adapt_a or args.legacy):parser.error('TDDFT requires --adapt-a or --legacy')
    from mpi4py import MPI
    world=MPI.COMM_WORLD;nf=world.size-2
    if nf not in (1,2) or (args.field=='multiscale' and nf!=1):raise InputError('invalid rank allocation')
    dt=.003;steps=64
    plan=ParallelPlan((SolverGroup('field',tuple(range(nf)),role='field'),
                       SolverGroup('left',(nf,),point_ids=(10,)),SolverGroup('right',(nf+1,),point_ids=(20,))))
    with ParallelContext.create(world,plan) as context:
        try:
            field=None;native=args.field=='native'
            positions=np.array([[-1.,1.,1.],[5.,1.,1.]]) if native else np.array([[2.,0.,0.],[4.,0.,0.]])
            phase=.5 if native else 0.
            ports=(InputPort('A',offset=phase),) if args.ports=='a' else (InputPort('E'),InputPort('H',offset=.5))
            units={'A':'A/c:au'} if args.ports=='a' else {'E':'E:au','H':'H:au'}
            contracts={p:FieldContract(TimeSpec(dt,ports,(OutputPort('j',offset=(0. if native else .5) if args.reject_phase else phase),)),
                                          units,positions[i],'custom') for i,p in enumerate(plan.point_ids)}
            if context.group.role=='field':
                if native:
                    calc=FDTD(Geometry(Grid((40,40,40),(2.,)*3)),(),dt,steps,
                              source=PlaneWave(AcosPulse(.001,20.,.15),0,-16.),observations=tuple(map(tuple,positions)))
                    field=NativeEHField(context,calc,args.build,args.output/'field',points={10:(19,20,20),20:(22,20,20)},grid=(nf,1,1))
                else:
                    model=Multiscale1D(2,2.,vacuum=(10,10),incident=AcosPulse(.001,20.,.15))
                    field=MultiscaleField(context,model.new_state(dt=dt,steps=steps))
            provider=TDDFTMaterialProvider(args.bundle,'rt.inp',runtime=args.runtime,max_steps=steps,threads=1) if args.runtime else Provider()
            backend=TDDFTResponseBackend() if args.runtime else Backend()
            if args.adapt_a and not args.runtime:
                original_create=backend.create_state
                def create(*values):
                    state=original_create(*values)
                    state.timing=TimeSpec(values[2],(InputPort('A'),),(OutputPort('j'),))
                    return state
                backend.create_state=create
            binding=MaterialBinding('response',provider,backend)
            if args.adapt_a:binding=packetize_A_binding(binding,field_phase='yee_half_step' if native else 'endpoint')
            response=ResponseSession(context,{p:binding for p in plan.point_ids},dt=dt,
                                     field_contracts=None if args.legacy else contracts)
            try:
                if args.reject_phase:
                    try:PacketCoupledStepper(response,field)
                    except InputError:pass
                    else:raise AssertionError('invalid current phase accepted')
                    assert response.step_index==0
                    if field is not None:assert field.step_index==0
                    if world.rank==0:print('PASS phase rejected before update')
                    return
                driver=CoupledStepper(response,field) if args.legacy else PacketCoupledStepper(response,field)
                declarations=response.timing_records()
                currents=[];electric=[];potential=[]
                for _ in range(steps):
                    result=driver.step()
                    if world.rank==0:
                        currents.append(result['J_charge']);electric.append(field.E.copy());potential.append(field.A.copy())
                if world.rank==0:
                    args.output.mkdir(parents=True,exist_ok=native)
                    np.savez(args.output/'history.npz',j=currents,E=electric,A=potential)
                    (args.output/'timing.json').write_text(json.dumps({'coupling':driver.timing,'materials':declarations},indent=2))
                    print('PASS packet coupling' if not args.legacy else 'PASS legacy reference',args.output)
            finally:
                response.close()
                if args.runtime:provider.close()
                if field is not None:field.close()
        except BaseException:
            # Abort before context-manager cleanup can strand other solver groups.
            import traceback
            traceback.print_exc()
            world.Abort(1)
            raise


if __name__=='__main__':
    # Parsing/help remains usable without importing MPI.
    main()
