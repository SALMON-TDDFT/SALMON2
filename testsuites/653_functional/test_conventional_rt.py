"""Ordinary hybrid GS -> native RT, with real MPI layouts and saved-input guards.

Set SALMON_TEST_EXE and SALMON_TEST_MPIEXEC to run these small integration tests.
Each producer is generated afresh, without DC or fractional occupations.
"""
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FUNCTIONALS = ('pbe0', 'pbeh40', 'pbeh40_rvv10')


@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE') and os.environ.get('SALMON_TEST_MPIEXEC'),
                     'SALMON_TEST_EXE and SALMON_TEST_MPIEXEC required')
class ConventionalHybridRT(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='conventional-hybrid-rt-')
        cls.addClassCleanup(cls.temp.cleanup)
        cls.root = Path(cls.temp.name)
        cls.producers = {}

    def execute(self, name, inp, ranks=1, threads=1):
        folder = self.root / name
        folder.mkdir()
        shutil.copy(ROOT / 'testsuites/pseudo/H_rps.dat', folder)
        (folder / 'inputfile').write_text(inp)
        command = [os.environ['SALMON_TEST_EXE']]
        if ranks > 1:
            command = [os.environ['SALMON_TEST_MPIEXEC'], '-n', str(ranks)] + command
        result = subprocess.run(command, input=inp, cwd=folder, text=True, capture_output=True,
                                timeout=180, env=dict(os.environ, OMP_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS='1'))
        (folder / 'outputfile').write_text(result.stdout + result.stderr)
        if os.environ.get('SALMON_TEST_SAVE_DIR'):
            shutil.copytree(folder, Path(os.environ['SALMON_TEST_SAVE_DIR']) / self.root.name / name)
        return folder, result

    def ground_state(self, functional, multik=False):
        key = (functional, multik)
        if key not in self.producers:
            template = ROOT / 'testsuites/431_H_pbe0_conventional_gs/inputfile'
            inp = template.read_text().replace("xc='pbe0'", f"xc='{functional}'")
            inp = inp.replace('nproc_k=2', 'nproc_k=1')
            if not multik:
                inp = inp.replace('num_kgrid=1,2,1', 'num_kgrid=1,1,1')
            folder, result = self.execute(f'gs_{functional}_{int(multik)}', inp)
            self.assertEqual(result.returncode, 0, result.stdout[-2500:] + result.stderr)
            self.assertIn('#GS converged at', result.stdout)
            self.assertIn('end SALMON', result.stdout)
            self.assertTrue((folder / 'data_for_restart/hybrid_gs.bin').is_file())
            if multik:
                self.assertIn('EXX_DISTRIBUTED_K:', result.stdout)
                self.assertNotIn('EXX_WANNIER refresh/', result.stdout)
            self.producers[key] = folder
        return self.producers[key]

    def rt_input(self, functional, producer, multik=False, ranks=1, field='zero'):
        inp = (ROOT / 'testsuites/432_H_pbe0_conventional_rt/inputfile').read_text()
        inp = inp.replace("xc='pbe0'", f"xc='{functional}'")
        inp = re.sub(r"directory_read_data='[^']+'",
                     f"directory_read_data='{producer / 'data_for_restart'}/'", inp)
        inp = inp.replace('nt=16', 'nt=8')
        if multik:
            inp = inp.replace('nproc_k=2', f'nproc_k={ranks}')
        else:
            inp = inp.replace('num_kgrid=1,2,1', 'num_kgrid=1,1,1').replace('nproc_k=2', 'nproc_k=1')
            inp = inp.replace('nproc_rgrid=1,1,1', f'nproc_rgrid=1,{ranks},1')
        if field == 'zero':
            inp = inp.replace('e_impulse=1d-4', 'e_impulse=0d0')
        elif field == 'pulse':
            inp = inp.replace("theory='tddft_response'", "theory='tddft_pulse'")
            inp = inp.replace("ae_shape1='impulse'", "ae_shape1='Acos2'")
            inp = inp.replace('e_impulse=1d-4',
                              'E_amplitude1=.001d0\n omega1=10d0\n tw1=.16d0\n t1_start=0d0\n phi_CEP1=0d0')
        return inp

    def assert_rt(self, folder, result):
        self.assertEqual(result.returncode, 0, result.stdout[-2500:] + result.stderr)
        self.assertIn('end SALMON', result.stdout)
        self.assertNotIn('DC-LCFO wavefunction reconstruction', result.stdout)
        self.assertNotIn('Native LCFO RT active', result.stdout)
        current = np.loadtxt(next(folder.glob('*_rt.data')))
        energy = np.loadtxt(next(folder.glob('*_rt_energy.data')))
        self.assertEqual(current.shape, (8, 16))
        self.assertEqual(len(energy), 9)
        self.assertTrue(np.isfinite(current).all())
        self.assertTrue(np.isfinite(energy).all())
        self.assertAlmostEqual(current[-1, 0], .16, places=10)
        charges = re.findall(r'^\s*\d+\s+[\d.]+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s*$', result.stdout, re.M)
        self.assertTrue(charges, 'Missing electron norm')
        self.assertTrue(all(np.isfinite(float(value)) for value in charges))
        self.assertLess(max(abs(float(value) - 4.) for value in charges), 1e-7)
        return current, energy

    def check_layouts(self, field, multik):
        for functional in FUNCTIONALS:
            with self.subTest(functional=functional, multik=multik, field=field):
                producer = self.ground_state(functional, multik)
                reference = None
                for ranks in (1, 2):
                    inp = self.rt_input(functional, producer, multik, ranks, field)
                    folder, result = self.execute(f'{field}_{functional}_{int(multik)}_{ranks}', inp, ranks)
                    current, energy = self.assert_rt(folder, result)
                    if multik:
                        self.assertIn('EXX_DISTRIBUTED_K:', result.stdout)
                        self.assertNotIn('EXX_WANNIER refresh/', result.stdout)
                    if field == 'zero':
                        self.assertLess(np.max(np.abs(energy[:, 1] - energy[0, 1])), 1e-7)
                        self.assertLess(np.max(np.abs(current[:, 13:16])), 1e-6)
                    else:
                        self.assertGreater(np.max(np.abs(current[:, 13:16])), 1e-12)
                    if reference is not None:
                        np.testing.assert_allclose(current, reference[0], atol=1e-8, rtol=1e-7)
                        np.testing.assert_allclose(energy, reference[1], atol=1e-8, rtol=1e-7)
                    reference = (current, energy)

    def test_multik_fractional_source_keeps_occupation_aware_route(self):
        inp = (ROOT / 'testsuites/431_H_pbe0_conventional_gs/inputfile').read_text()
        inp = inp.replace('nproc_k=2', 'nproc_k=1')
        inp = inp.replace('nstate=2', 'nstate=3\n temperature_k=300d0')
        folder, result = self.execute('gs_thermal_extra_state', inp)
        self.assertEqual(result.returncode, 0, result.stdout[-2000:] + result.stderr)
        self.assertIn('#GS converged at', result.stdout)
        self.assertIn('EXX_WANNIER refresh/', result.stdout)
        self.assertNotIn('EXX_DISTRIBUTED_K:', result.stdout)

    def test_hse_multik_snapshot_remains_available(self):
        inp = (ROOT / 'testsuites/431_H_pbe0_conventional_gs/inputfile').read_text()
        inp = inp.replace('nproc_k=2', 'nproc_k=1')
        inp = inp.replace("xc='pbe0'", "xc='hse06'\n yn_hse_wannier='y'\n yn_hse_wannier_snapshot='y'")
        folder, result = self.execute('gs_hse_snapshot', inp)
        self.assertEqual(result.returncode, 0, result.stdout[-2000:] + result.stderr)
        self.assertIn('#GS converged at', result.stdout)
        self.assertTrue(list(folder.rglob('hse_wannier_snapshot.bin')))

    def test_gamma_zero_field_and_spatial_mpi(self):
        self.check_layouts('zero', False)

    def test_full_k_mesh_impulse_and_k_mpi(self):
        self.check_layouts('impulse', True)

    def test_full_k_mesh_acos2_and_k_mpi(self):
        self.check_layouts('pulse', True)

    def test_gamma_acos2_and_spatial_mpi(self):
        self.check_layouts('pulse', False)

    def test_gamma_adaptive_source_ace_and_spatial_mpi(self):
        for functional in FUNCTIONALS:
            with self.subTest(functional=functional):
                producer = self.ground_state(functional)
                reference = None
                for ranks in (1, 2):
                    inp = self.rt_input(functional, producer, ranks=ranks, field='impulse')
                    inp = inp.replace("exx_mlwf_interval=5",
                                      "exx_mlwf_norm_fraction=.999d0\n exx_ace_support='source'\n exx_mlwf_interval=5")
                    folder, result = self.execute(f'adaptive_{functional}_{ranks}', inp, ranks)
                    current, energy = self.assert_rt(folder, result)
                    self.assertIn('EXX_SUPPORT_ACE accepted: T', result.stdout)
                    self.assertNotIn('EXX_SUPPORT_ACE accepted: F', result.stdout)
                    self.assertGreater(np.max(np.abs(current[:, 13:16])), 1e-12)
                    if reference is not None:
                        np.testing.assert_allclose(current, reference[0], atol=1e-8, rtol=1e-7)
                        np.testing.assert_allclose(energy, reference[1], atol=1e-8, rtol=1e-7)
                    reference = (current, energy)

    def test_openmp_native_rt_parity(self):
        for functional in FUNCTIONALS:
            for multik in (False, True):
                with self.subTest(functional=functional, multik=multik):
                    producer = self.ground_state(functional, multik)
                    inp = self.rt_input(functional, producer, multik, 2, 'pulse')
                    if not multik:
                        inp = inp.replace("exx_mlwf_radius=0d0",
                                          "exx_mlwf_radius=0d0\n exx_mlwf_norm_fraction=.999d0\n exx_ace_support='source'")
                    reference = None
                    for threads in (1, 2):
                        folder, result = self.execute(f'omp_{functional}_{multik}_{threads}', inp, 2, threads)
                        current, energy = self.assert_rt(folder, result)
                        if reference is not None:
                            np.testing.assert_allclose(current, reference[0], atol=1e-8, rtol=1e-7)
                            np.testing.assert_allclose(energy, reference[1], atol=1e-8, rtol=1e-7)
                        reference = (current, energy)

    def test_openmp_dc_core_exchange_parity(self):
        template = (ROOT / 'testsuites/425_H_pbeh40_dc_gs/inputfile').read_text()
        for functional in ('hse06',) + FUNCTIONALS:
            with self.subTest(functional=functional):
                reference = None
                for threads in (1, 2):
                    inp = template.replace("xc='pbeh40'", f"xc='{functional}'")
                    _, result = self.execute(f'omp_dc_{functional}_{threads}', inp, 2, threads)
                    self.assertEqual(result.returncode, 0, result.stdout[-2500:] + result.stderr)
                    self.assertIn('end SALMON', result.stdout)
                    residual = re.findall(r'DC #SCF.*diff =\s*(\S+)', result.stdout)
                    self.assertLess(float(residual[-1]), 1e-10)
                    exchange = float(re.findall(r'DC_HSE_CORE exchange Ha =\s*(\S+)', result.stdout)[-1])
                    charge = float(re.findall(r'integral\(rho_tot\)=\s*(\S+)', result.stdout)[-1])
                    self.assertTrue(np.isfinite(exchange))
                    self.assertLess(abs(charge - 4.), 1e-10)
                    if reference is not None:
                        self.assertLess(abs(exchange - reference), 1e-9)
                    reference = exchange

    def test_metadata_rejections_before_wavefunction_read(self):
        original = self.ground_state('pbeh40')
        for mode in ('missing', 'truncated', 'functional', 'cutoff', 'ions', 'k_mesh', 'cell'):
            with self.subTest(mode=mode):
                producer = self.root / ('bad_metadata_' + mode)
                shutil.copytree(original, producer)
                metadata = producer / 'data_for_restart/hybrid_gs.bin'
                if mode == 'missing':
                    metadata.unlink()
                elif mode == 'truncated':
                    metadata.write_bytes(metadata.read_bytes()[:11])
                inp = self.rt_input('pbeh40', producer)
                if mode == 'functional':
                    inp = inp.replace("xc='pbeh40'", "xc='pbe0'")
                elif mode == 'cutoff':
                    inp = inp.replace('pbeh_coulomb_radius=4d0', 'pbeh_coulomb_radius=3d0')
                elif mode == 'ions':
                    inp = inp.replace('3.3d0', '3.31d0')
                elif mode == 'k_mesh':
                    inp = inp.replace('num_kgrid=1,1,1', 'num_kgrid=1,2,1')
                elif mode == 'cell':
                    inp = inp.replace('al=16d0,8d0,8d0', 'al=17d0,8d0,8d0')
                # A missing wfn makes the validation ordering observable.
                (producer / 'data_for_restart/wfn.bin').unlink()
                _, result = self.execute('reject_' + mode, inp)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('Conventional hybrid GS metadata missing, malformed, or mismatched',
                              result.stdout + result.stderr)
                self.assertNotIn('end SALMON', result.stdout)

    def test_rvv10_parameter_mismatch(self):
        producer = self.ground_state('pbeh40_rvv10')
        inp = self.rt_input('pbeh40_rvv10', producer)
        inp = inp.replace("xc='pbeh40_rvv10'", "xc='pbeh40_rvv10'\n rvv10_b=6d0")
        _, result = self.execute('reject_rvv10', inp)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Conventional hybrid GS metadata missing, malformed, or mismatched',
                      result.stdout + result.stderr)

    def test_saved_occupation_payload_is_checked(self):
        original = self.ground_state('pbeh40')
        producer = self.root / 'bad_occupation'
        shutil.copytree(original, producer)
        occupation = producer / 'data_for_restart/occupation.bin'
        data = bytearray(occupation.read_bytes())
        # Compiler-native sequential records: locate the two-real payload via
        # its matching record markers instead of assuming one marker width.
        for marker_size, marker_format in ((4, '=i'), (8, '=q')):
            if len(data) == 16 + 2 * marker_size and struct.unpack_from(marker_format, data)[0] == 16:
                self.assertEqual(struct.unpack_from(marker_format, data, len(data) - marker_size)[0], 16)
                struct.pack_into('=d', data, marker_size, 1.5)
                break
        else:
            self.fail('Unexpected occupation record layout')
        occupation.write_bytes(data)
        _, result = self.execute('reject_occupation', self.rt_input('pbeh40', producer))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Conventional hybrid GS occupation payload mismatch', result.stdout + result.stderr)
        self.assertNotIn('end SALMON', result.stdout)


if __name__ == '__main__':
    unittest.main()
