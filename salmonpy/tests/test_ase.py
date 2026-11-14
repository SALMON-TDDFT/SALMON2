import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np

HAS_ASE=importlib.util.find_spec('ase') and importlib.util.find_spec('f90nml')

@unittest.skipUnless(HAS_ASE,'optional ASE dependencies')
class ASEBridgeTest(unittest.TestCase):
    def test_structure_units_and_mismatch(self):
        from salmonpy.ase import atoms_from_initial, validate_structure, BOHR_ANGSTROM
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as tmp:
            text="""&calculation theory='dft' /
&system yn_periodic='y', al_vec1=4,0,0, al_vec2=1,5,0, al_vec3=0,1,6, natom=2 /
&pseudo izatom(1)=14 /
&atomic_red_coor
'Si' 0 0 0 1
'Si' .25 .25 .25 1
/
"""
            initial=SimpleNamespace(_input=text,_snapshot=Path(tmp))
            atoms=atoms_from_initial(initial)
            np.testing.assert_allclose(atoms.cell.array, np.array([[4,0,0],[1,5,0],[0,1,6]])*BOHR_ANGSTROM)
            np.testing.assert_allclose(atoms.get_scaled_positions(),[[0,0,0],[.25,.25,.25]])
            validate_structure(initial,atoms)
            changed=atoms.copy();changed.positions[1,0]+=.01
            with self.assertRaisesRegex(ValueError,'positions'):validate_structure(initial,changed)
            changed=atoms.copy();changed.numbers[0]=6
            with self.assertRaisesRegex(ValueError,'species'):validate_structure(initial,changed)
            changed=atoms.copy();changed.cell[0,0]+=.01
            with self.assertRaisesRegex(ValueError,'cell'):validate_structure(initial,changed)

    def test_disconnected_path_and_energy_axes(self):
        from ase.cell import Cell
        from salmonpy.ase import explicit_path, to_ase_band_structure, BOHR_ANGSTROM
        from salmonpy.bands import BandStructure,HARTREE_EV
        cell=Cell([[4.,0,0],[1.,5.,0],[0,1.,6.]])
        specials={'G':[0,0,0],'X':[.5,0,0],'M':[.5,.5,0],'R':[.5,.5,.5]}
        path=cell.bandpath('GX,MR',npoints=12,special_points=specials)
        copied=path.kpts.copy();p=explicit_path(path)
        np.testing.assert_array_equal(p.fractional,copied)
        reciprocal=2*np.pi*np.linalg.inv(cell.array/BOHR_ANGSTROM)
        energies=np.arange(len(copied)*3*2).reshape(len(copied),3,2)*.01
        bands=BandStructure(p.cartesian(reciprocal),energies)
        result=to_ase_band_structure(bands,path,reference_ev=2.5)
        np.testing.assert_array_equal(result.energies,energies.transpose(2,0,1)*HARTREE_EV)
        self.assertEqual(result.reference,2.5)
        self.assertEqual(result.path.path,'GX,MR')
        path.kpts[0,0]=.25
        np.testing.assert_array_equal(result.path.kpts,copied)
        np.testing.assert_array_equal(p.fractional,copied)
        with self.assertRaisesRegex(ValueError,'k'):to_ase_band_structure(bands,path,reference_ev=0)

    def test_cartesian_angstrom_file_and_periodic_equivalence(self):
        from salmonpy.ase import atoms_from_initial, validate_structure
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'atoms.txt').write_text("'Si' 1 2 3 1\n")
            initial=SimpleNamespace(_snapshot=root,_input="""&calculation theory='dft' /
&units unit_system='A_eV_fs' /
&system yn_periodic='y', al=4,5,6, natom=1, file_atom_coor='atoms.txt' /
&pseudo izatom(1)=14 /
""")
            atoms=atoms_from_initial(initial)
            np.testing.assert_allclose(atoms.positions,[[1,2,3]])
            np.testing.assert_allclose(atoms.cell.array,np.diag([4,5,6]))
            atoms.positions[0]+=atoms.cell[1]
            validate_structure(initial,atoms)
            atoms.pbc[2]=False
            with self.assertRaisesRegex(ValueError,'periodic'):validate_structure(initial,atoms)

    def test_ase_bohr_roundoff_accepted_consistently(self):
        from ase.cell import Cell
        from ase.units import Bohr
        from salmonpy.ase import to_ase_band_structure
        from salmonpy.bands import BandStructure
        path=Cell(np.eye(3)*10.26*Bohr).bandpath('GX',npoints=5)
        native=BandStructure(path.kpts*(2*np.pi/10.26),np.zeros((len(path.kpts),2,1)))
        result=to_ase_band_structure(native,path,reference_ev=0)
        np.testing.assert_array_equal(result.path.kpts,path.kpts)

    def test_dft_model_from_atoms_is_independent(self):
        from ase import Atoms
        from salmonpy.ase import dft_from_atoms,BOHR_ANGSTROM
        from salmonpy.dft import Pseudopotential
        atoms=Atoms('Si',positions=[[0,0,0]],cell=np.eye(3)*5,pbc=True)
        pseudo=Pseudopotential(14,2,'_rps.dat',b'2000 .01 2 4.0\n')
        args=dict(pseudopotentials=[pseudo],electrons=4,states=8,grid=[12]*3,kmesh=[2]*3)
        model=dft_from_atoms(atoms,**args)
        atoms.positions[0,0]=1
        np.testing.assert_array_equal(model.fractional_positions,[[0,0,0]])
        np.testing.assert_allclose(model.cell,np.eye(3)*5/BOHR_ANGSTROM)
        atoms.set_initial_magnetic_moments([1])
        with self.assertRaises(ValueError):dft_from_atoms(atoms,**args)
        atoms.set_initial_magnetic_moments([0]);atoms.pbc[2]=False
        with self.assertRaises(ValueError):dft_from_atoms(atoms,**args)

class OptionalDependencyTest(unittest.TestCase):
    def test_import_does_not_load_ase(self):
        import subprocess,sys
        subprocess.run([sys.executable,'-c',"import sys; import salmonpy.ase; assert 'ase' not in sys.modules; assert 'f90nml' not in sys.modules"],check=True)

@unittest.skipUnless(HAS_ASE, 'optional ASE dependencies')
class ASEOpticalTest(unittest.TestCase):
    def inputs(self):
        from ase import Atoms
        from salmonpy.dft import Pseudopotential
        atoms = Atoms('Si2', scaled_positions=[[0,0,0],[.25,.25,.25]], cell=np.eye(3)*5, pbc=True)
        options = dict(pseudopotentials=[Pseudopotential(14,2,'_rps.dat',b'2000 .01 2 4.0\n')],
                       states=8, grid=(12,12,12), kmesh=(2,2,2))
        return atoms, options

    def test_neutral_count_and_explicit_mismatch(self):
        from salmonpy.ase import dft_from_atoms
        atoms, options = self.inputs()
        self.assertEqual(dft_from_atoms(atoms, **options).electrons, 8)
        self.assertEqual(dft_from_atoms(atoms, electrons=8, **options).electrons, 8)
        for count in (4, 10, True, 8.5):
            with self.subTest(count=count), self.assertRaises(ValueError):
                dft_from_atoms(atoms, electrons=count, **options)

    def test_optical_entry_composes_existing_runner(self):
        from unittest.mock import patch
        from salmonpy.ase import tddft_from_atoms
        from salmonpy.dft import DFT, GroundStateResult, SCFSettings
        from salmonpy.fdtd import SalmonExecutable
        atoms, options = self.inputs()
        with tempfile.TemporaryDirectory() as tmp:
            engine = SalmonExecutable(Path(tmp)/'fake', '0'*64)
            gs = GroundStateResult(Path(tmp)/'gs', -1., 10)
            marker = object()
            scf = SCFSettings(tolerance=1e-10)
            with patch.object(DFT, 'run', return_value=gs) as run, patch.object(GroundStateResult, 'tddft_initial_state', return_value=marker) as load:
                initial = tddft_from_atoms(atoms, gs.directory, executable=engine, max_steps=20,
                                          threads=2, timeout=60., scf=scf, **options)
                self.assertIs(initial, marker)
                run.assert_called_once_with(engine, gs.directory, scf=scf, threads=2, timeout=60.)
                load.assert_called_once_with(executable=engine, max_steps=20, threads=2, timeout=60.)
            for extra in ({'max_steps':0}, {'threads':0}, {'electrons':4}):
                with patch.object(DFT, 'run') as run, self.assertRaises(ValueError):
                    tddft_from_atoms(atoms, gs.directory, executable=engine, **options, **extra)
                run.assert_not_called()
