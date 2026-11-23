"""Write FDTD inputs only; never launch SALMON or submit a job.

python -m salmonpy.examples.prepare_fdtd NEW_DIRECTORY --mpi-grid 2 1 1
This small LD fixture demonstrates staging, not converged optical physics.
"""
import argparse
import json
from pathlib import Path
from salmonpy.fdtd import FDTD,FDTDMPI
from salmonpy.geometry import Geometry,Grid,Shape
from salmonpy.ld import LorentzDrude
from salmonpy.preparation import prepare_fdtd


def calculation():
    material=LorentzDrude(plasma_frequency=.3,strength=[1.],resonance=[0.],damping=[.01])
    return FDTD(Geometry(Grid((40,40,40),(2.,2.,2.)),(Shape('ellipsoid',(8.,8.,8.)),)),
                (material,),dt=.003,steps=32,observations=((0.,0.,0.),))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--mpi-grid',type=int,nargs=3,default=(2,1,1))
    args=parser.parse_args()
    mpi=FDTDMPI(tuple(args.mpi_grid))
    directory=prepare_fdtd(calculation(),args.directory,mpi=mpi)
    print(json.dumps({'status':'prepared-not-run','directory':str(directory),
                      'ranks':mpi.ranks,'scope':'input files only; no solver, launcher or job submission'},indent=2))


if __name__=='__main__':main()
