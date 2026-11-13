import tempfile
import unittest
import importlib.util
from pathlib import Path
import numpy as np
from salmonpy.bands import BandStructure

class BandsTest(unittest.TestCase):
    def test_blocks_and_spins(self):
        text = ('Number_of_Bands: 2\nNumber_of_kpt_in_each_block: 1\nNumber_of_blocks: 2\n'
                '1 9 9 9 0 0 0\n1 1 -1 -2\n1 2 3 4\n'
                '2 9 9 9 .3 .4 0\n1 1 -3 -4\n1 2 5 6\n')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'band.dat';p.write_text(text)
            b=BandStructure.from_salmon(p)
            np.testing.assert_allclose(b.distance,[0,.5])
            self.assertEqual(b.energies.shape,(2,2,2))
            self.assertEqual(b.energies[1,0,1],-4)
            with self.assertRaises(ValueError): b.energies.setflags(write=True)
            p.write_text(text.rsplit('\n',2)[0])
            with self.assertRaises(ValueError): BandStructure.from_salmon(p)
    def test_invalid_data(self):
        with self.assertRaises(ValueError): BandStructure([[0,0,0],[1,0,0]], [[[float('nan')]],[[1]]])
    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "optional plot dependency")
    def test_plot(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        b=BandStructure([[0,0,0],[1,0,0]], [[[1]],[[2]]])
        ax=b.plot(reference_ev=1,ticks=[(0,'Γ'),(1,'X')])
        self.assertAlmostEqual(ax.lines[0].get_ydata()[0],26.211386245988)
        plt.close(ax.figure)
        with self.assertRaises(ValueError): b.plot(reference_ev=0,ticks=[(2,'X')])
