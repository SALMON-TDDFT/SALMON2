"""Transverse-periodic slab, x PML; synthetic optical input, atomic units.

python -m salmonpy.examples.fdtd_periodic_slab RUNTIME NEW_OUTPUT [--ld]
Requires the runtime's per-axis FDTD boundary implementation. No build/download.
"""
import argparse
from salmonpy.runtime import load_runtime
from salmonpy.fdtd import FDTD, Boundaries, Dielectric, PlaneWave
from salmonpy.geometry import Grid, Geometry, Shape
from salmonpy.fields import AcosPulse
from salmonpy.ld import LorentzDrude

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('runtime');p.add_argument('directory');p.add_argument('--ld',action='store_true')
p.add_argument('--epsilon-infinity',type=float,default=1.0,help='instantaneous LD background, used with --ld')
a=p.parse_args()
material=(LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[25.],damping=[2.],epsilon_infinity=a.epsilon_infinity)
          if a.ld else Dielectric(epsilon=2.25))
calculation=FDTD(Geometry(Grid((256,8,8),(.5,.5,.5)),(Shape('box',(20.,4.,4.)),)),
    (material,),.0015,4000,
    boundaries=Boundaries(('pml','pml','periodic','periodic','periodic','periodic')),
    source=PlaneWave(AcosPulse(.001,10.,2.),0,-40.),observations=((30.,-.25,-.25),(30.,.25,.25)))
result=load_runtime(a.runtime,required=['fdtd']).run(calculation,a.directory,threads=1)
print(result.directory)
print('For vacuum-normalized complex Fresnel comparison, use tools/test_fdtd_periodic.py.')
