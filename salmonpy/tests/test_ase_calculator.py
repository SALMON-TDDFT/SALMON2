import importlib.util
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np

@unittest.skipUnless(importlib.util.find_spec('ase'),'optional ASE dependency')
class ASECalculatorTest(unittest.TestCase):
    def test_cache_invalidation_units_and_failure(self):
        from ase import Atoms
        from ase.calculators.calculator import PropertyNotImplementedError
        from salmonpy.ase_calculator import SalmonCalculator
        from salmonpy.dft import DFT,Pseudopotential
        from salmonpy.fdtd import SalmonExecutable
        from salmonpy.bands import HARTREE_EV
        from salmonpy.ase import BOHR_ANGSTROM
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as tmp:
            calc=SalmonCalculator(executable=SalmonExecutable(Path(tmp)/'fake','0'*64),directory=Path(tmp)/'runs',
                 pseudopotentials=[Pseudopotential(14,2,'_rps.dat',b'fixture')],electrons=4,states=8,grid=[12]*3,kmesh=[2]*3)
            atoms=Atoms('Si',positions=[[0,0,0]],cell=np.eye(3)*5,pbc=True);atoms.calc=calc
            fake=SimpleNamespace(energy_hartree=-1.,forces_hartree_per_bohr=np.array([[.1,.2,.3]]))
            with patch.object(DFT,'run',return_value=fake) as run:
                self.assertAlmostEqual(atoms.get_potential_energy(),-HARTREE_EV)
                f=atoms.get_forces();np.testing.assert_allclose(f,[[.1,.2,.3]]*np.array(HARTREE_EV/BOHR_ANGSTROM))
                self.assertEqual(run.call_count,1)
                f[:]=0
                self.assertNotEqual(atoms.get_forces()[0,0],0)
                atoms.positions[0,0]+=.1;atoms.get_forces();self.assertEqual(run.call_count,2)
                atoms.cell[0,0]+=.1;atoms.get_potential_energy();self.assertEqual(run.call_count,3)
                dirs=[str(c.args[1]) for c in run.call_args_list];self.assertEqual(len(set(dirs)),3)
                with self.assertRaises(PropertyNotImplementedError):atoms.get_stress()
                with self.assertRaises(ValueError):calc.set(grid=(16,16,16))
            atoms.positions[0,0]+=.1
            with patch.object(DFT,'run',side_effect=RuntimeError('unconverged')):
                with self.assertRaisesRegex(RuntimeError,'unconverged'):atoms.get_forces()
            self.assertEqual(calc.results,{})
            self.assertIsNone(calc.last_ground_state)
