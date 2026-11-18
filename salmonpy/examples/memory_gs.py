"""Read existing GS, make independent transverse TDDFT states, no SCF.

mpiexec -n 2 python -m salmonpy.examples.memory_gs BUILD EXISTING_GS NEW_WORK
The source GS must match its rt.inp and remain unchanged through load.
"""
import argparse
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from salmonpy.native_gs import NativeGSBackend


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build',type=Path)
    parser.add_argument('initial',type=Path)
    parser.add_argument('work',type=Path)
    parser.add_argument('--describe-input',action='store_true',help='show captured declared input metadata; requires metadata extra')
    args=parser.parse_args()
    from mpi4py import MPI
    comm=MPI.COMM_WORLD
    if comm.rank==0:args.work.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    context=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    backend=NativeGSBackend(args.build,args.work,decomposition=(1,1,comm.size,1,1))
    with backend.open(context) as session:
        gs=session.load(args.initial,material_id='silicon')
        model=gs.tddft_model(dt=.003,max_steps=24)
        description=None;error=None
        if args.describe_input:
            if comm.rank==0:
                try:description=gs.input_metadata()
                except Exception as exc:error=str(exc)
            error=comm.bcast(error,root=0)
            if error:raise RuntimeError('input metadata description failed: '+error)
        gs.close()
        if description is not None and comm.rank==0:
            import json
            print(json.dumps(description.as_dict(),ensure_ascii=False,indent=2))
        first=model.new_state('first');second=model.new_state('second')
        try:
            for step in range(1,17):
                A=1e-4*np.array([np.sin(.31*step),np.cos(.23*step)-1,np.sin(.17*step)])
                reference=first.step(A)
                np.testing.assert_array_equal(second.step(A),reference)
            if comm.rank==0:print('Two independent states agree after 16 steps; no SCF run.')
        finally:
            first.close();second.close();model.close()


if __name__=='__main__':main()
