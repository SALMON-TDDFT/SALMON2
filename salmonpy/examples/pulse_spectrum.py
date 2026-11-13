"""python -m salmonpy.examples.pulse_spectrum (no native solver required)."""
import numpy as np
from salmonpy.fields import AcosPulse, sample_pulses
from salmonpy.analysis import spectrum

pulses=[AcosPulse(.01,.2,100.,start=10.),AcosPulse(.01,.2,100.,start=130.)]
a=sample_pulses(pulses,np.arange(0,250,.1))
e=a.electric_field()
s=spectrum(e.time,e.values[:,2],window='hann')
print(dict(s.config))
print('Peak angular frequency:',s.omega[np.argmax(abs(s.values))])
