#!/usr/bin/env python3
"""Production V2 reader tests using synthetic paired regular Pauli primitives."""
from pathlib import Path
import importlib.util
import math
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('legacy_reader_tests', ROOT / 'testsuites/450_ky_soc_nlcc/test_reader.py')
legacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)


def synthetic_v2():
    lines = (ROOT / 'testsuites/pseudo/Se_rps.dat').read_text().splitlines()
    n, h = int(lines[0].split()[0]), float(lines[0].split()[1])
    nc = int(lines[3].split()[2])
    length = 0.529177210903
    lines[1] = 'KY_SOC_NLCC_V2'
    aux = []
    for i in range(n + 1):
        r = i * h / length
        a = math.exp(-r*r)
        da = -2*r*a
        b = 0.05*r*a
        db = 0.05*(1-2*r*r)*a
        rho = 2*(a*a+b*b)/(4*math.pi)
        drho = 4*(a*da+b*db)/(4*math.pi)
        dda = (4*r*r-2)*a
        ddb = 0.05*(4*r*r*r-6*r)*a
        ddrho = 4*(da*da+a*dda+db*db+b*ddb)/(4*math.pi)
        angular = b/r if r else db
        tau = 2*(da*da+db*db+2*angular*angular)/(4*math.pi)
        current = -4*b*b/(8*math.pi*r) if r else 0.0
        row = lines[4+nc+i].split()
        row[2:5] = [f'{rho/length**3:.16e}',f'{drho/length**4:.16e}',f'{ddrho/length**5:.16e}']
        lines[4+nc+i] = ' '.join(row)
        aux.append(' '.join(f'{x:.16e}' for x in [r,rho,tau,current,a,da,b,db]))
    return '\n'.join(lines[:-1] + ['CORE_AUX_PAULI_V1',
        'PAULI_CLOSED_SHELL_V1 ATOMIC_UNITS TAU_NO_HALF',f'{n+1} 2 2',
        '1 0 0 2','1 1 -2 2',*aux,'END_KY_SOC_NLCC_V2'])+'\n'


class Core23Reader(legacy.Reader):
    def trial(self, text):
        p = self.d / 'synthetic_v2_rps.dat'
        p.write_text(text)
        return self.run_reader(p)

    def test_v2_valid_regular_joint_primitives(self):
        q = self.trial(synthetic_v2())
        self.assertEqual(q.returncode,0,q.stderr+q.stdout)

    def test_v2_requires_model_units_complete_block(self):
        text = synthetic_v2()
        for bad in [text.replace('PAULI_CLOSED_SHELL_V1','UNKNOWN_CORE'),
                    text.replace('TAU_NO_HALF','TAU_HALF'),
                    text.replace('CORE_AUX_PAULI_V1','MISSING_BLOCK'),
                    text.replace('1 1 -2 2','1 1 1 2'),
                    text.replace('1 1 -2 2','2 1 -2 2')]:
            q = self.trial(bad)
            self.assertNotEqual(q.returncode,0,q.stdout)

    def test_v2_scalar_nonfinite_and_closure_rejected(self):
        lines = synthetic_v2().splitlines()
        row = lines.index('CORE_AUX_PAULI_V1')+6
        for column,value in [(1,'NaN'),(2,'Infinity'),(3,'NaN'),(4,'NaN'),(5,'NaN'),(2,'1.23456789')]:
            bad = list(lines)
            fields = bad[row].split();fields[column]=value;bad[row]=' '.join(fields)
            q = self.trial('\n'.join(bad)+'\n')
            self.assertNotEqual(q.returncode,0,(column,value,q.stdout))


if __name__ == '__main__':
    unittest.main()
