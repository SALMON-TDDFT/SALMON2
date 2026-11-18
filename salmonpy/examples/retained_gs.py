"""Small retained-GS ownership demo; four SCF iterations, not converged material data.

mpiexec -n 2 python -m salmonpy.examples.retained_gs BUILD NEW_OUTPUT
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from salmonpy.dft import DFT,Pseudopotential,SCFSettings
from salmonpy.native_gs import NativeGSBackend


def main():
    from mpi4py import MPI
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--export-archives',action='store_true',
                        help='explicitly save Gamma GS/2 arrays; serial only, not a full SCF checkpoint')
    args=parser.parse_args();world=MPI.COMM_WORLD;error=None
    if args.export_archives and world.size!=1:
        raise ValueError('--export-archives requires one rank')
    if world.size not in (1,2):raise ValueError('this small demo uses one or two ranks')
    if world.rank==0:
        try:args.output.mkdir(parents=True,exist_ok=False)
        except Exception as exc:error=str(exc)
    error=world.bcast(error,root=0)
    if error:raise ValueError(error)
    pp=Pseudopotential.from_file(Path(__file__).with_name('si_intro_data')/'Si_rps.dat',atomic_number=14,local_channel=2)
    models=[DFT((14,),np.eye(3)*8,[[.3,.4,.45]],(pp,),4,2,(8,8,8),(1,1,1) if args.export_archives else (2,1,1),xc='pz'),
            DFT((14,14),np.eye(3)*9,[[.3,.4,.45],[.7,.6,.55]],(pp,),8,4,(12,12,12),(1,1,1) if args.export_archives else (1,2,1),xc='libxc_pbe')]
    context=SimpleNamespace(solver_comm=world,closed=False,_dependents=0)
    backend=NativeGSBackend(args.build,args.output,decomposition=(1,world.size,1,1,1))
    with backend.open(context) as session:
        states=[session.run(model,scf=SCFSettings(max_iterations=4),material_id=name)
                for model,name in zip(models,('silicon-a','silicon-b'))]
        reports=[{'material_id':s.material_id,**vars(s.info),
                  'eigenvalues_hartree':s.array('eigenvalues').values.tolist()} for s in states]
        if args.export_archives:
            for state in states:
                state.export_archive(args.output/(state.material_id+'-gs2'))
        states[0].close()
        assert states[1].info.energy_hartree==reports[1]['energy_hartree']
        if world.rank==0:
            (args.output/'result.json').write_text(json.dumps({'status':'pass','materials':reports,
                'scope':'retained ownership demonstration; not converged physics'},indent=2)+'\n')
    assert context._dependents==0


if __name__=='__main__':main()
