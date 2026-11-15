import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
from salmonpy.errors import InputError


@unittest.skipUnless(importlib.util.find_spec('ase'), 'optional ASE')
class StructureTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('salmonpy.structure'), 'structure preparation missing')
        from salmonpy import structure
        self.api=structure
        from ase import Atoms
        self.atoms=Atoms('Si2',positions=[[0,0,0],[1,2,3]],cell=[5,6,7],pbc=True)

    def test_atoms_are_detached_without_wrapping(self):
        self.atoms.positions[0]=[-1,0,0]
        value=self.api.StructureInput.from_atoms(self.atoms)
        self.atoms.positions[:]=0
        np.testing.assert_array_equal(value.atoms.positions[0],[-1,0,0])
        changed=value.atoms; changed.positions[:]=9
        np.testing.assert_array_equal(value.atoms.positions[0],[-1,0,0])
        self.assertEqual(value.record()['atom_count'],2)

    def test_file_frame_and_hash(self):
        from ase.io import write
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'super@cell.xyz'
            second=self.atoms.repeat((2,1,1))
            write(path,[self.atoms,second],format='extxyz')
            value=self.api.StructureInput.from_file(path,format='extxyz',index=1)
            self.assertEqual(len(value.atoms),4)
            np.testing.assert_allclose(value.atoms.cell.array,second.cell.array)
            self.assertEqual(value.record()['sha256'],hashlib.sha256(path.read_bytes()).hexdigest())
            self.assertEqual(value.record()['index'],1)
            self.assertNotIn(tmp,str(value.record()))
            with self.assertRaises(InputError): self.api.StructureInput.from_file(path,format='auto')
            with self.assertRaises(InputError): self.api.StructureInput.from_file(path,format='extxyz',index=-1)

    def test_grid_ceiling_alignment_and_units(self):
        selection=self.api.grid_from_spacing(self.atoms.cell.array,spacing=0.4,unit='angstrom')
        self.assertEqual(selection.counts,(13,15,18))
        np.testing.assert_allclose(selection.actual_spacing_angstrom,np.array([5,6,7])/selection.counts)
        self.assertTrue(all(h<=.4 for h in selection.actual_spacing_angstrom))
        aligned=self.api.grid_from_spacing(self.atoms.cell.array,spacing=.4,unit='angstrom',multiples=(4,2,4))
        self.assertEqual(aligned.counts,(16,16,20))
        from salmonpy.ase import BOHR_ANGSTROM
        same=self.api.grid_from_spacing(self.atoms.cell.array,spacing=.4/BOHR_ANGSTROM,unit='bohr')
        self.assertEqual(same.counts,selection.counts)
        with self.assertRaises(AttributeError): selection.counts=(1,1,1)

    def test_grid_rejects_invalid_geometry_and_settings(self):
        cell=self.atoms.cell.array
        for width in (0,-1,float('nan'),True,[.1,.2],1e-300):
            with self.assertRaises(InputError): self.api.grid_from_spacing(cell,spacing=width,unit='angstrom')
        for bad in (np.zeros((3,3)), [[5,1,0],[0,6,0],[0,0,7]],np.diag([-5,6,7])):
            with self.assertRaises(InputError): self.api.grid_from_spacing(bad,spacing=.4,unit='angstrom')
        with self.assertRaises(InputError): self.api.grid_from_spacing(cell,spacing=.4,unit='nm')
        with self.assertRaises(InputError): self.api.grid_from_spacing(cell,spacing=.4,unit='angstrom',multiples=(0,1,1))

    def test_dft_width_matches_explicit_counts(self):
        from salmonpy.ase import dft_from_atoms
        from salmonpy.dft import Pseudopotential
        with tempfile.TemporaryDirectory() as tmp:
            # Use the same minimal KY fixture pattern as existing DFT tests.
            path=Path(tmp)/'Si_rps.dat'
            path.write_text('1 0.1 2 4\n')
            pseudo=Pseudopotential.from_file(path,atomic_number=14,local_channel=2)
            args=dict(pseudopotentials=[pseudo],states=4,kmesh=(1,1,1))
            old=dft_from_atoms(self.atoms,grid=(13,15,18),**args)
            new=dft_from_atoms(self.atoms,grid_spacing=.4,grid_spacing_unit='angstrom',**args)
            self.assertEqual(old.input_text(),new.input_text())
            for extra in ({'grid':(13,15,18),'grid_spacing':.4,'grid_spacing_unit':'angstrom'},
                          {'grid_spacing':.4},{}):
                with self.assertRaises(InputError): dft_from_atoms(self.atoms,**args,**extra)

    def test_cif_preserves_expanded_cell(self):
        from ase.io import write
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'cell.cif'
            write(path,self.atoms,format='cif')
            value=self.api.StructureInput.from_file(path,format='cif')
            np.testing.assert_allclose(value.atoms.cell.array,self.atoms.cell.array,atol=1e-12)
            np.testing.assert_allclose(value.atoms.positions,self.atoms.positions,atol=1e-12)

    def test_grid_requires_scalar_or_three_values(self):
        with self.assertRaises(InputError):
            self.api.grid_from_spacing(self.atoms.cell.array,spacing=[.4],unit='angstrom')

    def test_nested_metadata_is_independent(self):
        self.atoms.info['nested']={'values':[1]}
        value=self.api.StructureInput.from_atoms(self.atoms)
        self.atoms.info['nested']['values'][0]=2
        self.assertEqual(value.atoms.info['nested']['values'],[1])
        copied=value.atoms; copied.info['nested']['values'][0]=3
        self.assertEqual(value.atoms.info['nested']['values'],[1])

    def test_partial_occupancy_is_rejected_by_loader_and_dft(self):
        from salmonpy.ase import dft_from_atoms
        self.atoms.info['occupancy']={'0':{'Si':.6,'Ge':.4}}
        with self.assertRaisesRegex(InputError,'occupancy'):
            self.api.StructureInput.from_atoms(self.atoms)
        with self.assertRaisesRegex(InputError,'occupancy'):
            dft_from_atoms(self.atoms,pseudopotentials=[],states=4,kmesh=(1,1,1),grid=(8,8,8))
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'mixed.cif'
            path.write_text('''data_test
_cell_length_a 5
_cell_length_b 5
_cell_length_c 5
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Si1 Si 0 0 0 0.6
Ge1 Ge 0 0 0 0.4
''')
            with self.assertRaisesRegex(InputError,'occupancy'):
                self.api.StructureInput.from_file(path,format='cif')
