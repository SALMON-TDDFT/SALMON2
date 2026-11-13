import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from salmonpy.bands import BandPath, BandInitialState, _check_report, _validate_gs_input
from salmonpy.fdtd import SalmonExecutable

class BandRunnerTest(unittest.TestCase):
    def test_skew_path_and_ticks(self):
        p=BandPath([[0,0,0],[.5,.25,0],[0,0,0]],[2,1],['Γ','Q','Γ'])
        b=np.array([[2,1,0],[0,3,.5],[0,0,4.]])
        np.testing.assert_allclose(p.cartesian(b),[[0,0,0],[.625,.375,0],[1.25,.75,0],[0,0,0]])
        self.assertEqual(p.ticks,((0,'Γ'),(2,'Q'),(3,'Γ')))
        with self.assertRaises(ValueError):p.vertices.setflags(write=True)
        for intervals in ([0,1],[1.5,1],[True,1]):
            with self.assertRaises(ValueError):BandPath(p.vertices,intervals,p.labels)

    def test_report_nonconvergence_and_exact_division(self):
        # Columns of reciprocal lattice, then start/count/converged block records.
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'report'
            p.write_text('SALMONPY_BANDS_V1\n1 0 0\n0 1 0\n0 0 1\n2 3\n1 2 1\n3 2 1\n')
            np.testing.assert_equal(_check_report(p,4,3),np.eye(3))
            p.write_text(p.read_text().replace('3 2 1','3 2 0'))
            with self.assertRaisesRegex(RuntimeError,'converg'):_check_report(p,4,3)
            p.write_text('SALMONPY_BANDS_V1\n1 0 0\n0 1 0\n0 0 1\n2 3\n1 2 1\n')
            with self.assertRaises(RuntimeError):_check_report(p,4,3)

    def test_snapshot_independence_and_input_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);src=root/'gs';src.mkdir();(src/'data_for_restart').mkdir()
            for name in ['info.bin','occupation.bin','wfn.bin']:(src/'data_for_restart'/name).write_bytes(b'gs')
            (src/'Si.dat').write_text('pseudo')
            inp="&calculation theory='dft' /\n&system yn_periodic='y', nstate=4 /\n&functional xc='PZ' /\n&pseudo file_pseudo(1)='Si.dat' /\n&atomic_red_coor\n'Si' 0 0 0 1\n/\n"
            (src/'gs.inp').write_text(inp)
            exe=root/'salmon';exe.write_text('fake')
            engine=SalmonExecutable(exe,hashlib.sha256(exe.read_bytes()).hexdigest())
            with BandInitialState.from_salmon(src,'gs.inp',executable=engine) as gs:
                (src/'Si.dat').write_text('changed after loading')
                self.assertEqual((gs._snapshot/'Si.dat').read_text(),'pseudo')
                path=BandPath([[0,0,0],[.5,0,0]],[3],['Γ','X'])
                def fail(*args,**kwargs):
                    work=Path(kwargs['cwd']);(work/'Si.dat').write_text('mutated run')
                    raise RuntimeError('deliberate failure')
                with patch('subprocess.run',side_effect=fail):
                    for name in ['run1','run2']:
                        with self.assertRaisesRegex(RuntimeError,'deliberate'):gs.run(path,root/name,bands=3)
                self.assertEqual((gs._snapshot/'Si.dat').read_text(),'pseudo')
                self.assertEqual((src/'Si.dat').read_text(),'changed after loading')
            (src/'gs.inp').write_text(inp.replace('Si.dat','../outside.dat'))
            with self.assertRaises(ValueError):BandInitialState.from_salmon(src,'gs.inp',executable=engine)

    def test_commented_assets_and_inactive_coordinate_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'Si.dat').write_text('pseudo')
            text="&calculation theory='dft' /\n&system yn_periodic='y', nstate=4, file_atom_coor='none', file_atom_red_coor='none' /\n&pseudo\n! file_pseudo(1)='old.dat'\nfile_pseudo(1)='Si.dat' /\n"
            _validate_gs_input(text,root)
