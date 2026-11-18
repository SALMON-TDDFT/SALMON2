"""Two different TDDFT inputs in one process per MPI rank, one shared attachment.

mpiexec -n 2 python -m salmonpy.examples.tddft_models BUILD INITIAL_A INITIAL_B NEW_OUTPUT
Only verified existing bundles/runtime; no DFT, download or build is started.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from salmonpy.embedded_tddft import NativeTDDFTBackend,NativeTDDFTMaterialProvider


def main():
    from mpi4py import MPI
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build','initial_a','initial_b','output'):p.add_argument(name,type=Path)
    p.add_argument('--steps',type=int,default=16)
    p.add_argument('--dt-a',type=float,default=.003)
    p.add_argument('--dt-b',type=float,default=.004)
    a=p.parse_args();world=MPI.COMM_WORLD;error=None
    if world.rank==0:
        try:a.output.mkdir(parents=True,exist_ok=False)
        except Exception as exc:error=str(exc)
    error=world.bcast(error,root=0)
    if error:raise ValueError(error)
    context=SimpleNamespace(solver_comm=world,closed=False,_dependents=0)
    backend=NativeTDDFTBackend(a.build,a.output,decomposition=(1,world.size,1,1,1))
    providers=[NativeTDDFTMaterialProvider(path,'rt.inp',max_steps=a.steps+2) for path in (a.initial_a,a.initial_b)]
    registry=backend.open_registry(context);models=[];states=[];currents=[[],[]]
    try:
        for name,provider,dt in zip(('material-a','material-b'),providers,(a.dt_a,a.dt_b)):
            model=registry.open_model(provider.prepare(context),dt)
            models.append(model);states.append(model.new_state(name))
        for i in range(1,a.steps+1):
            for k,state in enumerate(states):currents[k].append(state.step([0.,1e-4*np.sin(.2*i),0.]))
        if world.rank==0:
            np.savez(a.output/'response.npz',time_a=np.arange(1,a.steps+1)*a.dt_a,
                     time_b=np.arange(1,a.steps+1)*a.dt_b,j_a=currents[0],j_b=currents[1])
            (a.output/'result.json').write_text(json.dumps({'status':'pass','models':[
                {'material_id':name,'identity':model.identity,'dt':model.dt,'capacity':model.capacity,
                 'decomposition':model.decomposition} for name,model in zip(('material-a','material-b'),models)],
                'units':{'time':'atomic time','j':'charge-current-density:au'},
                'scope':'ownership example; not a physical convergence certificate'},indent=2)+'\n')
    finally:
        for state in reversed(states):state.close()
        for model in reversed(models):model.close()
        registry.close()
    assert context._dependents==0


if __name__=='__main__':main()
