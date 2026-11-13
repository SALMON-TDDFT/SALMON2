from pathlib import Path
import tempfile
import unittest
import numpy as np
from salmonpy.dft import DFT, Pseudopotential, SCFSettings, _read_status, GroundStateResult

class DFTTest(unittest.TestCase):
    def test_model_copies_and_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'Si_rps.dat';path.write_bytes(b'original')
            pseudo=Pseudopotential.from_file(path,atomic_number=14,local_channel=2)
            cell=np.eye(3)*10.26;positions=np.array([[0.,0.,0.]])
            model=DFT([14],cell,positions,[pseudo],electrons=4,states=8,grid=[12]*3,kmesh=[2]*3)
            cell[0,0]=99;positions[0,0]=.5;path.write_bytes(b'changed')
            self.assertEqual(model.cell[0,0],10.26)
            self.assertEqual(model.fractional_positions[0,0],0)
            self.assertEqual(pseudo.data,b'original')
            text=model.input_text(SCFSettings(max_iterations=100,tolerance=1e-8))
            self.assertIn("theory='dft'",text)
            self.assertIn('nscf=100',text)
            self.assertIn('_rps.dat',text)
            with self.assertRaises(ValueError):model.cell.setflags(write=True)
            with self.assertRaises(ValueError):DFT([14],np.eye(3),[[0,0,0]],[pseudo],electrons=3,states=8,grid=[12]*3,kmesh=[2]*3)

    def test_native_status_and_incomplete_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'gs_status.dat'
            p.write_text('SALMONPY_GS_V1\n1 37\n1e-10 1e-9 -31.25\nrho_dne\n')
            self.assertEqual(_read_status(p)['iterations'],37)
            p.write_text(p.read_text().replace('1 37','0 37'))
            with self.assertRaisesRegex(RuntimeError,'converg'):_read_status(p)
            with self.assertRaises((ValueError,FileNotFoundError)):GroundStateResult.open(tmp)

    def test_extra_native_files_rejected_on_reopen(self):
        import json,hashlib
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);native=root/'native';(native/'data_for_restart').mkdir(parents=True)
            for name in ['info.bin','occupation.bin','wfn.bin']:(native/'data_for_restart'/name).write_bytes(b'x')
            (native/'gs.inp').write_text('input')
            (native/'gs_status.dat').write_text('SALMONPY_GS_V1\n1 37\n1e-10 1e-9 -31.25\nrho_dne\n')
            files={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in native.rglob('*') if p.is_file()}
            (root/'ground-state.json').write_text(json.dumps({'schema':'salmonpy.gs/1','status':'complete','energy_hartree':-31.25,'iterations':37,'files':files}))
            GroundStateResult.open(root)
            (native/'hubbard_param.dat').write_text('unexpected physical input')
            with self.assertRaisesRegex(ValueError,'files'):GroundStateResult.open(root)

    def test_forces_shape_units_and_readonly(self):
        from salmonpy.dft import _read_forces
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'forces'
            p.write_text('SALMONPY_FORCES_V1\n2\n.1 -.2 0\n-.1 .2 0\n')
            f=_read_forces(p,2)
            np.testing.assert_array_equal(f,[[.1,-.2,0],[-.1,.2,0]])
            with self.assertRaises(ValueError):f.setflags(write=True)
            with self.assertRaises(ValueError):_read_forces(p,3)
            p.write_text('SALMONPY_FORCES_V1\n1\nnan 0 0\n')
            with self.assertRaises(ValueError):_read_forces(p,1)
