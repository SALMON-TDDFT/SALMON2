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
                 pseudopotentials=[Pseudopotential(14,2,'_rps.dat',b'2000 .01 2 4.0\n')],electrons=4,states=8,grid=[12]*3,kmesh=[2]*3)
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

    def test_neutral_inference_and_mismatch_before_native_run(self):
        from ase import Atoms
        from salmonpy.ase_calculator import SalmonCalculator
        from salmonpy.dft import DFT, Pseudopotential
        from salmonpy.fdtd import SalmonExecutable
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as tmp:
            args = dict(executable=SalmonExecutable(Path(tmp)/'fake','0'*64), directory=Path(tmp)/'runs',
                        pseudopotentials=[Pseudopotential(14,2,'_rps.dat',b'2000 .01 2 4.0\n')],
                        states=8, grid=(12,)*3, kmesh=(2,)*3)
            atoms = Atoms('Si', cell=np.eye(3)*5, pbc=True)
            atoms.calc = SalmonCalculator(**args)
            fake = SimpleNamespace(energy_hartree=-1., forces_hartree_per_bohr=np.zeros((1,3)))
            with patch.object(DFT, 'run', autospec=True, return_value=fake) as run:
                atoms.get_potential_energy()
                self.assertEqual(run.call_args.args[0].electrons, 4)
            atoms.calc = SalmonCalculator(electrons=6, **args)
            with patch.object(DFT, 'run') as run, self.assertRaisesRegex(ValueError, 'neutral'):
                atoms.get_potential_energy()
            run.assert_not_called()

    def test_opt_in_warm_start_uses_only_compatible_success(self):
        from ase import Atoms
        from types import SimpleNamespace
        from salmonpy.dft import DFT, Pseudopotential
        from salmonpy.ase_calculator import SalmonCalculator
        from salmonpy.fdtd import SalmonExecutable
        with tempfile.TemporaryDirectory() as tmp:
            calc = SalmonCalculator(executable=SalmonExecutable(Path(tmp)/'fake','0'*64), directory=tmp,
                    pseudopotentials=[Pseudopotential(14,2,'_rps.dat',b'20 .01 2 4\n')],
                    states=8, grid=(12,)*3, kmesh=(2,)*3, warm_start=True)
            atoms = Atoms('Si', cell=np.eye(3)*5, pbc=True); atoms.calc=calc
            fake = SimpleNamespace(energy_hartree=-1.,forces_hartree_per_bohr=np.zeros((1,3)))
            with patch.object(DFT,'run',return_value=fake) as run:
                atoms.get_potential_energy()
                self.assertNotIn('initial_guess',run.call_args.kwargs)
                atoms.positions[0,0]+=.01; atoms.get_potential_energy()
                self.assertIs(run.call_args.kwargs['initial_guess'],fake)
                atoms.cell[0,0]+=.01; atoms.get_potential_energy()
                self.assertNotIn('initial_guess',run.call_args.kwargs)
            atoms.positions[0,0]+=.01
            with patch.object(DFT,'run',side_effect=RuntimeError('unconverged')):
                with self.assertRaises(RuntimeError): atoms.get_potential_energy()
            atoms.positions[0,0]+=.01
            with patch.object(DFT,'run',return_value=fake) as run:
                atoms.get_potential_energy()
                self.assertNotIn('initial_guess',run.call_args.kwargs)

    def test_libxc_and_nc_upf_scope(self):
        from salmonpy.ase_calculator import SalmonCalculator
        from salmonpy.dft import Pseudopotential
        from salmonpy.fdtd import SalmonExecutable
        with tempfile.TemporaryDirectory() as tmp:
            options = dict(executable=SalmonExecutable(Path(tmp)/'fake','0'*64), directory=tmp,
                           states=8, grid=(12,)*3, kmesh=(2,)*3)
            for xc, suffix, data in [('libxc_pz','_rps.dat',b'20 .01 2 4\n'),
                                    ('libxc_pbe','.upf',b'<UPF><PP_HEADER pseudo_type="NC" z_valence="4"/></UPF>')]:
                calc = SalmonCalculator(xc=xc, pseudopotentials=[Pseudopotential(14,2,suffix,data)], **options)
                self.assertEqual(calc.settings['xc'],xc)
            with self.assertRaises(ValueError):
                SalmonCalculator(xc='tbmbj',pseudopotentials=[Pseudopotential(14,2,'_rps.dat',b'20 .01 2 4')],**options)

    def test_nlcc_blocks_forces_but_preserves_energy_and_cache(self):
        from ase import Atoms
        from ase.calculators.calculator import PropertyNotImplementedError
        from salmonpy.ase_calculator import SalmonCalculator
        from salmonpy.dft import DFT, Pseudopotential
        from salmonpy.fdtd import SalmonExecutable
        from types import SimpleNamespace
        headers = ('core_correction="T"', 'core_correction="true"',
                   'core_correction=".TRUE."', 'core_correction="1"', 'core_correction="unknown"')
        data = [f'<UPF><PP_HEADER pseudo_type="NC" z_valence="4" {h}/></UPF>'.encode()
                for h in headers]
        # SALMON reads PP_NLCC itself; a false/missing header must not bypass protection.
        data += [b'<UPF><PP_HEADER pseudo_type="NC" z_valence="4" core_correction="F"/>'
                 b'<PP_NLCC size="1">1.0</PP_NLCC></UPF>',
                 b'<UPF><PP_HEADER pseudo_type="NC" z_valence="4"/><PP_NLCC/></UPF>']
        for raw in data:
            for xc in ('pz','libxc_pz','libxc_pbe'):
                with self.subTest(raw=raw,xc=xc),tempfile.TemporaryDirectory() as tmp:
                    calc=SalmonCalculator(executable=SalmonExecutable(Path(tmp)/'fake','0'*64),
                        directory=tmp,pseudopotentials=[Pseudopotential(14,2,'.upf',raw)],
                        states=8,grid=(12,)*3,kmesh=(2,)*3,xc=xc)
                    atoms=Atoms('Si',cell=np.eye(3)*5,pbc=True);atoms.calc=calc
                    # No force attribute: an energy request must not read forces.
                    fake=SimpleNamespace(energy_hartree=-1.)
                    with patch.object(DFT,'run',return_value=fake) as run:
                        with self.assertRaisesRegex(PropertyNotImplementedError,'NLCC'):
                            atoms.get_forces()
                        run.assert_not_called()
                        energy=atoms.get_potential_energy()
                        self.assertTrue(np.isfinite(energy));self.assertEqual(run.call_count,1)
                        self.assertNotIn('forces',calc.results)
                        self.assertNotIn('forces',calc.implemented_properties)
                        with self.assertRaisesRegex(PropertyNotImplementedError,'NLCC'):
                            atoms.get_forces()
                        with self.assertRaisesRegex(PropertyNotImplementedError,'NLCC'):
                            calc.calculate(atoms,properties=['energy','forces'])
                        self.assertEqual(run.call_count,1)
                        self.assertEqual(atoms.get_potential_energy(),energy)

    def test_non_nlcc_upf_still_returns_forces(self):
        from ase import Atoms
        from salmonpy.ase_calculator import SalmonCalculator
        from salmonpy.dft import DFT, Pseudopotential
        from salmonpy.fdtd import SalmonExecutable
        from types import SimpleNamespace
        for flag in ('', 'core_correction="F"', 'core_correction="false"',
                     'core_correction=".FALSE."', 'core_correction="0"'):
            with self.subTest(flag=flag),tempfile.TemporaryDirectory() as tmp:
                data=f'<UPF><PP_HEADER pseudo_type="NC" z_valence="4" {flag}/></UPF>'.encode()
                calc=SalmonCalculator(executable=SalmonExecutable(Path(tmp)/'fake','0'*64),
                    directory=tmp,pseudopotentials=[Pseudopotential(14,2,'.upf',data)],
                    states=8,grid=(12,)*3,kmesh=(2,)*3,xc='libxc_pbe')
                atoms=Atoms('Si',cell=np.eye(3)*5,pbc=True);atoms.calc=calc
                fake=SimpleNamespace(energy_hartree=-1.,forces_hartree_per_bohr=np.zeros((1,3)))
                with patch.object(DFT,'run',return_value=fake):
                    np.testing.assert_array_equal(atoms.get_forces(),np.zeros((1,3)))
