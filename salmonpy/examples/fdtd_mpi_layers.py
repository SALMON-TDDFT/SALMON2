"""Three LD materials inside existing MPI FDTD; no dedicated LD ranks.
Run from serial Python, NOT under mpiexec:
python -m salmonpy.examples.fdtd_mpi_layers MPI_BUILD_JSON NEW_DIRECTORY --grid 2 2 1
Synthetic finite structure, not converged multilayer-film optics. Atomic units.
"""
import argparse
import numpy as np
from salmonpy.fdtd import FDTD,FDTDMPI,SalmonExecutable,PlaneWave
from salmonpy.geometry import Grid,Geometry,Shape
from salmonpy.ld import LorentzDrude
from salmonpy.fields import AcosPulse

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('build_record');p.add_argument('directory')
p.add_argument('--grid',type=int,nargs=3,default=(2,2,1),metavar=('PX','PY','PZ'))
p.add_argument('--threads',type=int,default=1)
a=p.parse_args()
materials=(LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.]),
           LorentzDrude(plasma_frequency=15.,strength=[.7,.3],resonance=[25.,35.],damping=[2.,3.]),
           LorentzDrude(plasma_frequency=12.,strength=[.5,.3,.2],resonance=[0.,20.,30.],damping=[1.,2.,3.]))
geometry=Geometry(Grid((48,40,40),(2.,2.,2.)),tuple(
    Shape('box',(8.,16.,16.),(x,0.,0.),material=i) for i,x in enumerate((-8.,0.,8.),1)))
calculation=FDTD(geometry,materials,.003,160,
    source=PlaneWave(AcosPulse(.001,20.,.15),0,-24.),observations=((-16.,0.,0.),(0.,0.,0.),(16.,0.,0.)))
mpi=FDTDMPI(tuple(a.grid))
ids,counts=np.unique(geometry.material_ids(),return_counts=True)
print('MPI ranks:',mpi.ranks,'OpenMP threads/rank:',a.threads)
print('Material-cell counts (0=vacuum):',dict(zip(ids.tolist(),counts.tolist())))
print('Memory has not been measured; these grid/pulse values are a connection example.')
result=SalmonExecutable.from_mpi_build(a.build_record).run(calculation,a.directory,mpi=mpi,threads=a.threads)
print(result.directory)
