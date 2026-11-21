"""Run in the installed Python environment: python -m salmonpy.examples.ld_response"""
import argparse
import numpy as np
from salmonpy import LorentzDrude

parser=argparse.ArgumentParser(description="Synthetic LD response; optional reciprocal anisotropic tensor")
parser.add_argument("--anisotropic",action="store_true")
args=parser.parse_args()

material = LorentzDrude(
    plasma_frequency=0.5,
    strength=[0.8, 0.2],
    resonance=[0.0, 0.3],
    damping=[0.01, 0.02],
    epsilon_infinity=2.25,  # instantaneous background; P/J below are dispersive poles
)
if args.anisotropic:
    from salmonpy.ld import TensorLorentzDrude
    material = TensorLorentzDrude(plasma_frequency=.5,
        strength=[[[.8,.1,0.],[.1,.4,.05],[0.,.05,.3]],
                  [[.2,0.,0.],[0.,.3,0.],[0.,0.,.1]]],
        resonance=[0.,.3],damping=[.01,.02],
        epsilon_infinity=[[2.25,.1,0.],[.1,2.,0.],[0.,0.,2.5]])
dt = 0.1  # atomic time units
field = np.zeros((1001, 3))
time = np.arange(len(field))*dt
field[:, 0] = 1e-3*np.exp(-((time-50)/10)**2)*np.cos(0.2*(time-50))
response = material.respond(field, dt=dt)
print('Final time (a.u.):', response.time[-1])
print('Final current (a.u.):', response.current[-1])
print('epsilon tensor at omega=0.2 a.u.:', material.permittivity_tensor(0.2))

print('Final displacement D (a.u.):', material.displacement(field, response.polarization)[-1])
