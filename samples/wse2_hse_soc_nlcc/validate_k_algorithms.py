#!/usr/bin/env python3
"""Compare native k algorithms on one completed SOC/NLCC symmetry fixture."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--salmon', type=Path, required=True)
    ap.add_argument('--reference', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    binary, source, root = args.salmon.resolve(), args.reference.resolve(), args.output.resolve()
    assert json.loads((source/'results.json').read_text())['result'] == 'PASS'
    digest = hashlib.sha256(binary.read_bytes()).hexdigest()
    assert (source/'execution-binary.sha256').read_text().strip() == digest
    root.mkdir()
    report = {'binary_sha256': digest, 'source_fixture': str(source)}

    def run(name, text, sym, selected=None, rejection=None):
        d = root/name
        d.mkdir()
        (d/'inputfile').write_text(text)
        (d/'sym.dat').write_text(sym.read_text())
        with (d/'inputfile').open() as inp, (d/'run.log').open('w') as log:
            result = subprocess.run([str(binary)], cwd=d, stdin=inp, stdout=log,
                stderr=subprocess.STDOUT, env=dict(os.environ, OMP_NUM_THREADS='2',
                OPENBLAS_NUM_THREADS='1'), timeout=1800)
        log = (d/'run.log').read_text()
        if rejection:
            assert rejection in log and 'end SALMON' not in log, log[-2000:]
            report[name] = {'rejected': rejection}
        else:
            assert result.returncode == 0 and 'end SALMON' in log, log[-2000:]
            assert 'selected='+selected in log, log[-2000:]
        return d, log

    # Same initial density and input as auto GS, now with explicit pair selection.
    text = (source/'gs-off/inputfile').read_text().replace("xc='hse06'", "xc='hse06',exx_k_algorithm='pair'")
    d, log = run('gs-pair', text, source/'gs-off/sym.dat', selected='spinor_pair')
    residual = float(re.findall(r'iter and int_x.*?=\s*\d+\s+([0-9.Ee+-]+)', log)[-1])
    assert residual < 1e-8 and '#GS converged at' in log
    energy = float(re.findall(r'Total Energy=\s*([0-9.Ee+-]+)', log)[-1])
    delta = abs(energy-json.loads((source/'results.json').read_text())['gs-off']['printed_energy_eV'])
    assert delta < 1e-4
    report['gs-pair'] = {'residual': residual, 'auto_energy_difference_eV': delta}

    def rows(path):
        data = [[float(t) for t in line.split()] for line in path.read_text().splitlines()
                if line.strip() and not line.lstrip().startswith('#')]
        assert data and all(math.isfinite(x) for row in data for x in row), path
        assert all(len(row) == 16 for row in data), path
        return data

    for axis in ('x', 'y'):
        original = source/('rt-'+axis+'-on')
        text = (original/'inputfile').read_text()
        # Each comparison consumes exactly the same verified GS checkpoint.
        text = re.sub(r"directory_read_data='([^']+)'",
            lambda m: "directory_read_data='%s'" % (original/m.group(1)).resolve(), text)
        reference = rows(original/'Si_rt.data')
        for mode in ('auto', 'reference', 'pair', 'legacy-pair'):
            setting = "exx_k_algorithm='%s'" % mode
            if mode == 'legacy-pair':
                setting = "yn_exx_spinor_pair_fft='y'"
            selected = 'spinor_pair' if mode in ('pair', 'legacy-pair') else 'reference'
            d, _ = run('rt-'+axis+'-'+mode, text.replace("xc='hse06'", "xc='hse06',"+setting),
                       original/'sym.dat', selected=selected)
            actual = rows(d/'Si_rt.data')
            assert len(actual) == len(reference) == 20
            delta = max(abs(a-b) for r, s in zip(actual, reference) for a, b in zip(r, s))
            assert delta < 1e-10
            variables = (d/'variables.log').read_text()
            assert 'exx_k_algorithm_reason=' in variables
            report['rt-'+axis+'-'+mode] = {'steps': len(actual), 'max_auto_difference': delta}

    text = (source/'gs-off/inputfile').read_text()
    for name, setting, rejection in (
        ('conflict', "exx_k_algorithm='reference',yn_exx_spinor_pair_fft='y'",
         'Conflicting exx_k_algorithm and legacy exchange switches'),
        ('invalid', "exx_k_algorithm='invalid'", 'Invalid exx_k_algorithm or legacy yes/no value'),
        ('soc-orbits', "exx_k_algorithm='orbits'", 'FFT orbits require native CPU dense-k hybrid symmetry')):
        run(name, text.replace("xc='hse06'", "xc='hse06',"+setting), source/'gs-off/sym.dat', rejection=rejection)
    run('manual-theory', text.replace("theory='dft'", "theory='dft_md'").replace(
        "xc='hse06'", "xc='hse06',exx_k_algorithm='pair'"), source/'gs-off/sym.dat',
        rejection='Manual EXX k algorithm requires dft or fixed-ion RT')
    run('scalar-pair-symmetry', text.replace("spin='noncollinear'", "spin='unpolarized'")
        .replace("yn_spinorbit='y'", "yn_spinorbit='n'").replace('nstate=8', 'nstate=4')
        .replace("yn_symmetry='n'", "yn_symmetry='yyy'")
        .replace("xc='hse06'", "xc='hse06',exx_k_algorithm='pair'"), source/'gs-off/sym.dat',
        rejection='Manual scalar pair algorithm requires crystal symmetry off')
    report['result'] = 'PASS'
    (root/'results.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
