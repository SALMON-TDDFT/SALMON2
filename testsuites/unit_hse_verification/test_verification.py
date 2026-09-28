#!/usr/bin/env python3
"""Check the GS/RT verifiers with synthetic data, including rejected failures."""
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class VerificationTest(unittest.TestCase):
    def run_case(self, case, files, ok):
        with tempfile.TemporaryDirectory() as tmp:
            for name, content in files.items():
                path = Path(tmp) / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(content)
            script = ROOT / 'testsuites' / case / 'verification'
            # Also exercise absence of the Python 3-only math helpers.
            runner = ('import math,runpy,sys; sys.modules["pathlib"]=None; '
                      'del math.isfinite; del math.isclose; '
                      'runpy.run_path(' + repr(str(script)) + ',run_name="__main__")')
            result = subprocess.run([sys.executable, '-c', runner], cwd=tmp,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(result.returncode == 0, ok, result.stderr.decode())

    def test_gs(self):
        files = {'outputfile': '#GS converged at 87 : 0.43519093E-08\nend SALMON\n',
                 'variables.log': 'hse_omega(bohr^-1)=0.11\n',
                 'data_for_restart/wfn.bin': 'fixture',
                 'Si_eigen.data': '1 -1.0\n',
                 'Si_k.data': ''.join('%d 0 0 0 0.015625\n' % i for i in range(1, 65))}
        self.run_case('420_bulk_Si_hse_gs', files, True)
        for name, bad in [('outputfile', 'end SALMON\n'),
                          ('outputfile', '#GS converged at 87 : 1e-3\nend SALMON\n'),
                          ('Si_eigen.data', '1 nan\n')]:
            broken = dict(files); broken[name] = bad
            self.run_case('420_bulk_Si_hse_gs', broken, False)

    def test_rt(self):
        rt = []
        spectrum = []
        for step in range(1, 65):
            row = [0.] * 16
            row[0], row[3], row[15] = step * .16, 1e-4, 1e-6
            rt.append(' '.join(map(str, row)))
        for i in range(1, 1001):
            row = [0.] * 13
            row[0], row[3], row[6] = i * .001, 1e-3, 2e-3
            row[9] = 1 - 4 * math.pi / row[0] * row[6]
            row[12] = 4 * math.pi / row[0] * row[3]
            spectrum.append(' '.join(map(str, row)))
        files = {'outputfile': '10 1.600 0 0 0 32.000 -1\nend SALMON\n',
                 'variables.log': 'propagator=hse_taylor4\n',
                 'Si_rt.data': '\n'.join(rt),
                 'Si_rt_energy.data': '\n'.join('%s -1' % (i * .16) for i in range(65)),
                 'Si_response.data': '\n'.join(spectrum)}
        self.run_case('421_bulk_Si_hse_rt', files, True)
        for name, bad in [('outputfile', files['outputfile'].replace('32.000', '31.000')),
                          ('Si_rt_energy.data', files['Si_rt_energy.data'].replace('10.24 -1', '10.24 -2')),
                          ('Si_response.data', files['Si_response.data'].replace('0.001', 'nan', 1)),
                          ('Si_response.data', files['Si_response.data'].replace(spectrum[0], ' '.join(['0.001'] + ['0'] * 12), 1))]:
            broken = dict(files); broken[name] = bad
            self.run_case('421_bulk_Si_hse_rt', broken, False)


if __name__ == '__main__':
    unittest.main()
