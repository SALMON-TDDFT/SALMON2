"""Run after the documented 8-atom Si band.inp calculation (17 physical k points)."""
import argparse
import numpy as np
from salmonpy.bands import BandStructure, HARTREE_EV

parser = argparse.ArgumentParser()
parser.add_argument('band_dat')
parser.add_argument('figure')
args = parser.parse_args()
raw = BandStructure.from_salmon(args.band_dat)
# Four segments, each with four intervals; exclude native batch padding.
bands = BandStructure(raw.k_cartesian[:17], raw.energies[:17, :24])
# 32 electrons, unpolarized: 16 occupied bands. Path maximum, not global VBM.
reference = np.max(bands.energies[:, :16])*HARTREE_EV
ax = bands.plot(reference_ev=reference,
                ticks=[(0,'Γ'),(4,'X'),(8,'M'),(12,'R'),(16,'Γ')],
                title='Si: 8-atom cubic cell (folded bands)\nPZ, coarse 12³ grid; path valence maximum = 0 eV')
ax.set_ylim(-14,8)
ax.figure.savefig(args.figure, dpi=160)
