"""Explicit MPI shared-model fast-PR HSE example; requires a prepared matching GS.

mpiexec -n 2 python shared_fast_hse.py BUILD_JSON GS NEW_OUTPUT
No build/download. This checks composition/restart, not physical convergence.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from salmonpy.fast_hse import FastHSEMaterialProvider
from salmonpy.fast_hse_mpi import NativeFastHSEBackend


def main():
    from mpi4py import MPI
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('build','initial','output'):parser.add_argument(name,type=Path)
    parser.add_argument('--steps',type=int,default=16)
    parser.add_argument('--split',type=int,default=8)
    parser.add_argument('--dt',type=float,default=.001)
    parser.add_argument('--symmetry',action='store_true',help='explicit same-symmetry reduced-k model; both example drives along z')
    args=parser.parse_args();comm=MPI.COMM_WORLD
    if not 0<args.split<args.steps:parser.error('require 0 < split < steps')
    output=args.output.resolve();error=None
    if comm.rank==0:
        try:output.mkdir(parents=True,exist_ok=False)
        except OSError as exc:error=str(exc)
    error=comm.bcast(error,root=0)
    if error:raise ValueError(error)
    context=SimpleNamespace(solver_comm=comm,point_ids=('point-a','point-b'),closed=False,_dependents=0)
    initial=FastHSEMaterialProvider(args.initial,max_steps=args.steps).prepare(context)
    model=NativeFastHSEBackend(args.build,output,allow_symmetry=args.symmetry).open_model(initial,args.dt,context)
    states=[];continuous=[];second=[];resumed=[]
    saved=output/f'point-a-{comm.rank}.bin'
    def drive(index):
        t=index*args.dt
        return np.array([0.,0.,1e-4*np.sin(.3*t)])
    try:
        with model:
            a=model.new_state('point-a');states.append(a)
            b=model.new_state('point-b');states.append(b)
            for index in range(1,args.steps+1):
                continuous.append(np.r_[a.step(drive(index)),a.energy])
                bfield=np.array([0.,0.,-7e-5*np.sin(.2*index*args.dt)]) if args.symmetry else np.array([1e-4*np.sin(.2*index*args.dt),0.,0.])
                second.append(np.r_[b.step(bfield),b.energy])
                if index==args.split:a.checkpoint(saved,False)
            fresh=model.new_state('point-a');states.append(fresh)
            fresh.checkpoint(saved,True)
            for index in range(args.split+1,args.steps+1):
                resumed.append(np.r_[fresh.step(drive(index)),fresh.energy])
            np.testing.assert_allclose(resumed,np.asarray(continuous)[args.split:],rtol=1e-10,atol=1e-12)
            np.savez(output/f'trace-{comm.rank}.npz',point_a=continuous,point_b=second,resumed=resumed)
    finally:
        # Model exit prevents creation. Each state owns a lease and closes explicitly.
        for state in reversed(states):state.close()
        model.close()
    assert context._dependents==0
    if comm.rank==0:
        (output/'result.json').write_text(json.dumps({'status':'pass','steps':args.steps,'split':args.split,
            'replayed_steps':0,'model_identity':model.identity,'point_ids':['point-a','point-b'],
            'A_unit':'A/c atomic units','current_unit':'charge current density atomic units',
            'energy_unit':'Hartree','caveat':'Composition/restart check, not material convergence.'},indent=2)+'\n')
        print('Shared model: independent states and direct resume agree')


if __name__=='__main__':main()
