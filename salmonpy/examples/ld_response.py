"""Run from checkout root: python3 -m salmonpy.examples.ld_response"""
import numpy as np
from salmonpy import LorentzDrude

material = LorentzDrude(
    plasma_frequency=0.5,
    strength=[0.8, 0.2],
    resonance=[0.0, 0.3],
    damping=[0.01, 0.02],
)
dt = 0.1  # atomic time units
field = np.zeros((1001, 3))
time = np.arange(len(field))*dt
field[:, 0] = 1e-3*np.exp(-((time-50)/10)**2)*np.cos(0.2*(time-50))
response = material.respond(field, dt=dt)
print('Final time (a.u.):', response.time[-1])
print('Final current (a.u.):', response.current[-1])
print('epsilon at omega=0.2 a.u.:', material.permittivity(0.2))
