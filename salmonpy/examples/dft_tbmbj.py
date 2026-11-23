"""Prepare periodic Si states with built-in TB-mBJ (no Libxc requirement).

python -m salmonpy.examples.dft_tbmbj RUNTIME SI_RPS NEW_GS [--cval 1.2]
An explicit GS calculation; use a compute node. Values are starter settings,
not converged optical properties. Energy/forces are diagnostic, not ASE fitting.
"""
import argparse
import numpy as np
from salmonpy.dft import DFT,Pseudopotential,SCFSettings
from salmonpy.runtime import load_runtime

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('runtime');p.add_argument('pseudo');p.add_argument('directory')
p.add_argument('--cval',type=float,default=None,help='fixed c; omitted means native automatic c')
p.add_argument('--threads',type=int,default=1)
a=p.parse_args()
positions=[[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0],
           [.25,.25,.25],[.25,.75,.75],[.75,.25,.75],[.75,.75,.25]]
model=DFT([14]*8,np.eye(3)*10.26,positions,
    [Pseudopotential.from_file(a.pseudo,atomic_number=14,local_channel=2)],
    electrons=32,states=32,grid=(16,16,16),kmesh=(2,2,2),xc='tbmbj',cval=a.cval)
engine=load_runtime(a.runtime,required=['dft'])
result=model.run(engine,a.directory,
    scf=SCFSettings(max_iterations=300,tolerance=1e-9),threads=a.threads)
print(result.directory)
print('SCF converged:',result.converged,'iterations:',result.iterations)
# result.tddft_initial_state(executable=engine) retains the same cval in RT input.
