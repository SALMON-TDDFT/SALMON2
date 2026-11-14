"""Native E/H stepping/coupling regressions, launched by a bounded MPI controller.
Args: MODE BUILD_JSON SERIAL_RUNTIME GS_BUNDLE NEW_OUTPUT
MODE: zero (field+dummy), current (field+dummy), tddft (field+2x2 material ranks).
Use 2/3 ranks for zero/current and 5/6 for TDDFT (one/two field ranks).
"""
import json,sys
from pathlib import Path
import numpy as np
from mpi4py import MPI
from salmonpy.fdtd import FDTD,PlaneWave
from salmonpy.fields import AcosPulse
from salmonpy.geometry import Grid,Geometry,Shape
from salmonpy.ld import LorentzDrude
from salmonpy.distributed import ParallelPlan,SolverGroup,ParallelContext
from salmonpy.embedded_tddft import NativeTDDFTMaterialProvider,NativeTDDFTBackend
from salmonpy.material_response import MaterialBinding,ResponseSession
from salmonpy.coupled import NativeEHField,CoupledStepper
from salmonpy.runtime import load_runtime
from salmonpy.tddft import TDDFTInitialState

world=MPI.COMM_WORLD
mode,record,runtime,bundle,out=sys.argv[1:6];out=Path(out)
field_axis=int(sys.argv[6]) if len(sys.argv)>6 else 0
field_size=world.size-(4 if mode in ('tddft','mixed') else 1)
if field_size not in (1,2):raise ValueError('one or two field ranks required')
try:
    if world.rank==0:out.mkdir(parents=True,exist_ok=False)
    world.Barrier()
    dt=.003;steps=2 if mode=='current' else 64
    grid=Grid((40,40,40),(2.,2.,2.));points={1:(19,20,20),2:(20,20,20)}
    observations=tuple(tuple(grid.coordinates[d][p[d]] for d in range(3)) for p in points.values())
    geometry=Geometry(grid,(Shape('box',(4.,4.,4.),(9.,1.,1.)),)) if mode=='ld' else Geometry(grid)
    materials=(LorentzDrude(plasma_frequency=20.,strength=[.7,.3],resonance=[0.,25.],damping=[2.,3.]),) if mode=='ld' else ()
    if mode=='mixed':
        geometry=Geometry(grid,(Shape('box',(2.,4.,4.),(7.,1.,1.),material=1),Shape('box',(2.,4.,4.),(9.,1.,1.),material=2)))
        materials=(LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.]),
                   LorentzDrude(plasma_frequency=15.,strength=[.7,.3],resonance=[25.,35.],damping=[2.,3.]))
    calculation=FDTD(geometry,materials,dt,steps,
        source=PlaneWave(AcosPulse(0. if mode=='current' else .001,20.,.15),0,-16.),observations=observations)
    if mode in ('zero','ld') and world.rank==0:
        reference=load_runtime(runtime).run(calculation,out/'batch',threads=1)
    groups=[SolverGroup('field',tuple(range(field_size)),role='field')]
    if mode in ('tddft','mixed'):
        groups += [SolverGroup('left',(field_size,field_size+1),point_ids=(1,)),
                   SolverGroup('right',(field_size+2,field_size+3),point_ids=(2,))]
    else:groups += [SolverGroup('dummy',(field_size,),point_ids=(1,2))]
    plan=ParallelPlan(tuple(groups));history=[]
    with ParallelContext.create(world,plan) as context:
        field=None;session=None
        try:
            if context.group.role=='field':
                decomposition=[1,1,1];decomposition[field_axis]=field_size
                field=NativeEHField(context,calculation,record,out/'field',points=points,grid=decomposition)
            if mode in ('tddft','mixed'):
                provider=NativeTDDFTMaterialProvider(bundle,'rt.inp',max_steps=steps)
                backend=NativeTDDFTBackend(record,out,decomposition=(1,2,1,1,1))
                binding=MaterialBinding('same-GS',provider,backend)
                session=ResponseSession(context,{1:binding,2:binding},dt=dt)
                coupled=CoupledStepper(session,field)
                for n in range(steps):
                    frame=coupled.step()
                    if world.rank==0:history.append({k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in frame.items()})
            elif field is not None:
                J=np.zeros((2,3))
                for n in range(steps):
                    if mode=='current' and n==1:J[:]=[[1e-4,3e-4,-2e-4],[2e-4,-1e-4,4e-4]]
                    A=field.advance(J)
                    if world.rank==0:history.append({'A_over_c':A.tolist(),'E_at_points':field.E.tolist()})
                if mode=='current':
                    expected=-4*np.pi*dt*np.array([[1e-4,1.5e-4,-1e-4],[1.25e-4,-.5e-4,2e-4]])
                    np.testing.assert_allclose(field.E,expected,rtol=2e-14,atol=1e-18)
        finally:
            if session is not None:session.close()
            if field is not None:field.close()
    assert world.allreduce(1)==world.size
    if world.rank==0:
        if mode in ('zero','ld'):
            for i in (1,2):
                x=np.loadtxt(out/'batch'/f'obs{i}_at_point_rt.data')
                y=np.loadtxt(out/'field'/f'obs{i}_at_point_rt.data')
                np.testing.assert_array_equal(x,y)
        if mode in ('tddft','mixed'):
            replay=np.empty((steps,2,3))
            with TDDFTInitialState.from_salmon(bundle,'rt.inp',executable=load_runtime(runtime),max_steps=steps,threads=1) as initial:
                for p in range(2):
                    with initial.new_state(dt=dt) as state:
                        for n,frame in enumerate(history):replay[n,p]=state.step(frame['A_over_c'][p])  # TDDFTState returns charge current j
            actual=np.array([f['J_charge'] for f in history])
            np.testing.assert_allclose(actual,replay,rtol=2e-10,atol=2e-12)
            previous=np.zeros((2,3));old_A=np.zeros((2,3))
            for n,frame in enumerate(history,1):
                np.testing.assert_array_equal(frame['J_previous'],previous)
                np.testing.assert_allclose(frame['A_over_c'],old_A-dt*np.array(frame['E_at_points']),rtol=0,atol=1e-18)
                assert frame['electronic_field_clock']==(n+.5)*dt
                previous=np.array(frame['J_charge']);old_A=np.array(frame['A_over_c'])
            print('replay max error',float(np.max(np.abs(actual-replay))))
        (out/'history.json').write_text(json.dumps(history,indent=2))
        print(mode,'PASS',field_size,'field ranks')
except BaseException:
    import traceback
    traceback.print_exc();world.Abort(1)
