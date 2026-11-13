"""Validate the documented four-segment Si native band run (after completion)."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from salmonpy.bands import BandStructure

root=Path(sys.argv[1]);b=BandStructure.from_salmon(root/'band.dat')
vertices=np.array([[0,0,0],[.5,0,0],[.5,.5,0],[.5,.5,.5],[0,0,0]])
expected=np.concatenate([a+(z-a)*np.arange(4)[:,None]/4 for a,z in zip(vertices[:-1],vertices[1:])] + [vertices[-1:]])*2*np.pi/10.26
np.testing.assert_allclose(b.k_cartesian[:17],expected,atol=1e-15)
np.testing.assert_allclose(b.k_cartesian[17:],0,atol=1e-15)
# Each of the three blocks must reach convergence; native emits cycle afterwards.
blocks=(root/'stdout.log').read_text().split('iter_band_kpt=')[1:]
assert len(blocks)==3 and all('cycle!!!' in block for block in blocks)
# The repeated Gamma endpoint must reproduce the first Gamma spectrum.
np.testing.assert_allclose(b.energies[0,:24], b.energies[16,:24],atol=1e-4,rtol=0)
print('PASS: Cartesian path, endpoint/padding, three converged blocks, Gamma closure')
print('Gamma closure max (Ha):',np.max(abs(b.energies[0,:24]-b.energies[16,:24])))
