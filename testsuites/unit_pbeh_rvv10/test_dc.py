"""DC total-density integration; SALMON_TEST_MPIEXEC enables 2/4-rank checks."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE'), 'SALMON_TEST_EXE is required')
class DCTest(unittest.TestCase):
    def run_case(self, ranks=1, conventional=False, xc='pbeh40_rvv10', radius=0, initial_only=False, method='gauss'):
        inp = Path(__file__).with_name('dc_hydrogen.inp').read_text()
        inp = inp.replace("xc='pbeh40_rvv10'", f"xc='{xc}'")
        if radius>0:
            inp=inp.replace('hse_mlwf_maxiter=20',f'exx_mlwf_maxiter=100\n exx_mlwf_interval=5\n exx_mlwf_tolerance=1d-7\n exx_mlwf_radius={radius}')
        if ranks < 4:
            inp = inp.replace('nproc_k=2', 'nproc_k=1')
            inp = inp.replace('nproc_rgrid_tot=4,1,1', f'nproc_rgrid_tot={ranks},1,1')
        if ranks == 1:
            inp = inp.replace('num_fragment=2,1,1', 'num_fragment=1,1,1')
            inp = inp.replace('num_rgrid_buffer=4,0,0', 'num_rgrid_buffer=0,0,0')
        if conventional:
            inp = inp.replace("yn_dc='y'", "yn_dc='n'")
        if initial_only:
            inp=inp.replace('nscf=200',f"nscf=1\n method_init_wf='{method}'")
        command = [os.environ['SALMON_TEST_EXE']]
        if ranks > 1:
            command = [os.environ['SALMON_TEST_MPIEXEC'], '-n', str(ranks)] + command
        with tempfile.TemporaryDirectory() as tmp:
            shutil.copy(ROOT/'testsuites/pseudo/H_rps.dat', tmp)
            run = subprocess.run(command, input=inp, cwd=tmp, text=True,
                capture_output=True, timeout=180,
                env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))
            if initial_only:
                self.assertNotEqual(run.returncode,0)
                self.assertIn('localization not converged; SCF result rejected',run.stdout+run.stderr)
                return float(re.findall(r'DC #SCF.*Total Energy =\s*([\d.E+-]+)',run.stdout)[0])
            self.assertEqual(run.returncode, 0, run.stdout[-1000:]+run.stderr)
            self.assertIn('end SALMON', run.stdout)
            if conventional:
                self.assertIn('#GS converged', run.stdout)
            else:
                differences = re.findall(r'DC #SCF.*diff =\s*([\d.E+-]+)', run.stdout)
                self.assertTrue(differences)
                self.assertLess(float(differences[-1]), 1e-8)
            info = next(Path(tmp).rglob('*_info.data')).read_text()
            energy = float(re.search(r'Total energy \(eV\) =\s*([\d.E+-]+)', info)[1])
            return energy

    def test_one_fragment_matches_conventional(self):
        # Exercise both dispersion energy accounting and its SCF potential.
        dc = self.run_case()
        conventional = self.run_case(conventional=True)
        self.assertLess(abs(dc-conventional), 2e-6)
        without_dispersion = self.run_case(xc='pbeh40')
        self.assertGreater(abs(dc-without_dispersion), 1e-4)

    @unittest.skipUnless(os.environ.get('SALMON_TEST_MPIEXEC'), 'MPI launcher is required')
    def test_finite_radius_initial_rank_invariance_and_localization_guard(self):
        for method in ['gauss','random']:
            two=self.run_case(ranks=2,radius=4,initial_only=True,method=method)
            four=self.run_case(ranks=4,radius=4,initial_only=True,method=method)
            self.assertLess(abs(two-four),2e-6)

    @unittest.skipUnless(os.environ.get('SALMON_TEST_MPIEXEC'), 'MPI launcher is required')
    def test_two_fragments_rank_invariance(self):
        two = self.run_case(ranks=2)
        four = self.run_case(ranks=4)
        self.assertLess(abs(two-four), 2e-6)
        # These buffers cover the full cell, so decomposition should reproduce
        # the one-fragment result as well (no fragment truncation here).
        self.assertLess(abs(two-self.run_case()), 2e-6)

if __name__ == '__main__':
    unittest.main()
