"""Two custom dispersive layers (atomic units, synthetic materials).
python -m salmonpy.examples.fdtd_layers MANAGED_RUNTIME NEW_DIRECTORY
Requires runtime interface v2. One MPI rank, one thread; no compilation.
This finite structure illustrates the API, not converged infinite-film optics.
"""
import argparse
from salmonpy.runtime import load_runtime
from salmonpy.ld import LorentzDrude
from salmonpy.fdtd import FDTD, PlaneWave
from salmonpy.geometry import Grid, Geometry, Shape
from salmonpy.fields import AcosPulse

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('runtime'); p.add_argument('directory')
a = p.parse_args()
# Replace these synthetic coefficients with validated material data.
# Frequencies and damping are in atomic units. Each material can have a
# different number of poles. Dielectric(epsilon=..., mu=..., conductivity=...) is
# also accepted in this tuple for non-dispersive isotropic materials.
metal = LorentzDrude(plasma_frequency=20., strength=[1.], resonance=[0.], damping=[2.])
resonant = LorentzDrude(plasma_frequency=15., strength=[.7, .3],
                       resonance=[25., 35.], damping=[2., 3.])
geometry = Geometry(Grid((36,36,36),(2.,2.,2.)), (
    Shape('box',(4.,16.,16.),(-2.,0.,0.),material=1),
    Shape('box',(4.,16.,16.),(2.,0.,0.),material=2)))
calculation = FDTD(geometry, (metal, resonant), dt=.003, steps=160,
    source=PlaneWave(AcosPulse(.001,20.,.15,polarization=(0.,0.,1.)),axis=0,position=-14.),
    observations=((12.,0.,0.),(0.,0.,0.)))
result = load_runtime(a.runtime, required=['fdtd','material-id-independent']).run(
    calculation, a.directory, threads=1)
for table in result.observations:
    print(table.columns)
    print(table.values[-1])
