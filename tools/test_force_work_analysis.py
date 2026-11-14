"""Numerical sign/quadrature and scalar-only retention regressions."""
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import measure_nlcc_force_work as probe

class ForceWorkAnalysis(unittest.TestCase):
    def test_cubic_quadrature_and_conservative_force(self):
        x=np.linspace(0,2,17)
        self.assertAlmostEqual(probe.simpson(x**3,2),4)
        rows=[dict(energy_hartree=v*v,path_force_hartree_per_bohr=-2*v) for v in x]
        self.assertAlmostEqual(probe.summarize(rows,2)['consistency_residual_hartree'],0)
    def test_known_force_defect_sign(self):
        x=np.linspace(0,2,17)
        rows=[dict(energy_hartree=v*v,path_force_hartree_per_bohr=-2*v+.01) for v in x]
        result=probe.summarize(rows,2)
        self.assertAlmostEqual(result['consistency_residual_hartree'],.02)
        self.assertAlmostEqual(result['mean_force_consistency_defect_hartree_per_bohr'],-.01)
    def test_nonfinite_or_odd_intervals_rejected(self):
        for y in ([0,1],[0,np.nan,1],[0,np.inf,1]):
            with self.assertRaises(ValueError):probe.simpson(y,1)
    def test_resumed_spot_samples_sort_before_integration(self):
        x=np.linspace(0,1,17)
        rows=[dict(index=i,energy_hartree=v*v,path_force_hartree_per_bohr=-2*v) for i,v in enumerate(x)]
        unordered=[rows[0],rows[16]]+rows[1:16]
        self.assertAlmostEqual(probe.summarize(probe.ordered_samples(unordered,16),1)['consistency_residual_hartree'],0)
        with self.assertRaises(ValueError):probe.ordered_samples(rows+[rows[0]],16)
    def test_json_replacement_keeps_previous_samples_if_interrupted(self):
        with tempfile.TemporaryDirectory() as root:
            path=Path(root)/'samples.json';probe.write_json(path,{'old':1})
            with patch.object(Path,'replace',side_effect=OSError('interrupted')):
                with self.assertRaises(OSError):probe.write_json(path,{'new':2})
            import json
            self.assertEqual(json.loads(path.read_text()),{'old':1})

    def test_only_requested_reports_survive_success_or_timeout(self):
        for fail in (False,True):
            with self.subTest(fail=fail),tempfile.TemporaryDirectory() as root:
                dest=Path(root)/'out';scratch=[]
                def native(*args,**kw):
                    w=kw['cwd'];scratch.append(w)
                    for name in ('gs_eigen.data','gs_info.data','gs_k.data','wfn.bin',
                                 'gs_status.dat','gs_forces.dat'):
                        (w/name).write_text('fixture')
                    if fail:raise subprocess.TimeoutExpired('fixture',1)
                    return subprocess.CompletedProcess(args,0)
                with patch.object(probe.subprocess,'run',side_effect=native):
                    if fail:
                        with self.assertRaises(subprocess.TimeoutExpired):
                            probe.run_scalar_only(Path('/fake'),'input',b'potential',dest,{})
                    else:probe.run_scalar_only(Path('/fake'),'input',b'potential',dest,{})
                self.assertEqual({p.name for p in dest.iterdir()},set(probe.SCALAR_FILES))
                self.assertFalse(scratch[0].exists())

    def test_mpi_scalar_export_filter_also_runs_on_failure(self):
        with tempfile.TemporaryDirectory() as root:
            dest=Path(root)/'out';scratch=[]
            def mpi(command,text,**kw):
                work=kw['cwd'];scratch.append(work)
                for name in ('gs_eigen.data','gs_info.data','gs_k.data','gs_status.dat'):
                    (work/name).write_text('fixture')
                raise subprocess.TimeoutExpired('mpi-fixture',1)
            with patch('salmonpy.fdtd._run_fdtd_mpi',side_effect=mpi):
                with self.assertRaises(subprocess.TimeoutExpired):
                    probe.run_scalar_only(Path('/fake'),'input',b'potential',dest,{},launcher=('mpiexec','-n','2'))
            self.assertEqual({p.name for p in dest.iterdir()},{'gs.inp','stdout.log','gs_status.dat'})
            self.assertFalse(scratch[0].exists())

if __name__=='__main__':unittest.main()
