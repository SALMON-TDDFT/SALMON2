"""Small retained-GS ownership demo; four SCF iterations, not converged material data.

mpiexec -n 2 python -m salmonpy.examples.retained_gs BUILD NEW_OUTPUT
For existing fast GS: add --import-bundle GS [--second-bundle OTHER_GS].
Optional --current-steps 16 uses retained GS→RT without GS files (declared RT profile).
Computed mode performs four SCF iterations first; imported mode performs no SCF.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from salmonpy.dft import DFT,Pseudopotential,SCFSettings
from salmonpy.native_gs import NativeGSBackend


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--xc',choices=('pz','libxc_pz','libxc_pbe','hse06'),help='computed-mode XC override; hse06 requires a compute-capable fast-HSE build')
    parser.add_argument('--import-bundle',type=Path,help='read existing compatible GS without SCF')
    parser.add_argument('--second-bundle',type=Path)
    parser.add_argument('--current-steps',type=int,default=0,help='optional current steps; requires a declared retained-GS RT profile')
    parser.add_argument('--export-archives',action='store_true',
                        help='explicitly save regular-mesh GS/2 arrays on subgroup root; not a full SCF checkpoint')
    args=parser.parse_args()
    if args.import_bundle and args.xc is not None:
        parser.error('--xc selects computed GS; imported GS retains its recorded XC')
    if args.current_steps<0:
        parser.error('--current-steps requires a nonnegative count')
    from mpi4py import MPI
    world=MPI.COMM_WORLD;error=None
    if world.size not in (1,2):raise ValueError('this small demo uses one or two ranks')
    if world.rank==0:
        try:args.output.mkdir(parents=True,exist_ok=False)
        except Exception as exc:error=str(exc)
    error=world.bcast(error,root=0)
    if error:raise ValueError(error)
    def response_from(gs):
        if not args.current_steps:return None
        model=gs.tddft_model(dt=.001,max_steps=args.current_steps)
        gs.close()  # RT seed remains owned independently of the source GS.
        first=model.new_state('pulse-a');second=model.new_state('pulse-b')
        try:
            currents=[]
            for i in range(1,args.current_steps+1):
                A=1e-4*np.array([np.sin(.23*i),np.cos(.17*i)-1,np.sin(.13*i)])
                currents.append(first.step(A).tolist());second.step(-A)
            return {'material_id':model.material_id,'dt_au':.001,
                    'charge_current_au':currents,'replayed_steps':0}
        finally:
            first.close();second.close();model.close()

    if args.import_bundle:
        context=SimpleNamespace(solver_comm=world,closed=False,_dependents=0)
        backend=NativeGSBackend(args.build,args.output,decomposition=(world.size,1,1,1,1))
        with backend.open(context) as session:
            states=[session.load(bundle,material_id=name) for bundle,name in zip(
                (args.import_bundle,args.second_bundle or args.import_bundle),('material-a','material-b'))]
            before=states[1].array('density').values
            reports=[{'material_id':state.material_id,**vars(state.info),
                      'density_shape':list(state.array('density').values.shape)} for state in states]
            if args.export_archives:
                for state in states:state.export_archive(args.output/(state.material_id+'-gs2'))
            response=response_from(states[0])
            states[0].close()
            np.testing.assert_array_equal(before,states[1].array('density').values)
            if world.rank==0:
                (args.output/'result.json').write_text(json.dumps({'status':'pass','materials':reports,
                    'SCF_calls':0,'response':response,'scope':'existing GS import; convergence/energy not established'},indent=2)+'\n')
        assert context._dependents==0
        return
    pp=Pseudopotential.from_file(Path(__file__).with_name('si_intro_data')/'Si_rps.dat',atomic_number=14,local_channel=2)
    models=[DFT((14,),np.eye(3)*8,[[.3,.4,.45]],(pp,),4,2,(8,8,8),(2,1,1),xc=args.xc or 'pz'),
            DFT((14,14),np.eye(3)*9,[[.3,.4,.45],[.7,.6,.55]],(pp,),8,4,(12,12,12),(1,2,1),xc=args.xc or 'libxc_pbe')]
    context=SimpleNamespace(solver_comm=world,closed=False,_dependents=0)
    decomposition=(world.size,1,1,1,1) if args.xc=='hse06' else (1,world.size,1,1,1)
    backend=NativeGSBackend(args.build,args.output,decomposition=decomposition)
    with backend.open(context) as session:
        states=[session.run(model,scf=SCFSettings(max_iterations=4),material_id=name)
                for model,name in zip(models,('silicon-a','silicon-b'))]
        reports=[{'material_id':s.material_id,**vars(s.info),
                  'eigenvalues_hartree':s.array('eigenvalues').values.tolist()} for s in states]
        if args.export_archives:
            for state in states:
                state.export_archive(args.output/(state.material_id+'-gs2'))
        response=response_from(states[0])
        states[0].close()
        assert states[1].info.energy_hartree==reports[1]['energy_hartree']
        if world.rank==0:
            (args.output/'result.json').write_text(json.dumps({'status':'pass','materials':reports,
                'response':response,'scope':'four-iteration retained ownership demonstration; not converged physics'},indent=2)+'\n')
    assert context._dependents==0


if __name__=='__main__':main()
