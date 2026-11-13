"""Check excluded inputs against a built executable without loading GS data."""
from pathlib import Path
import os
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[2]
executable = Path(sys.argv[1]).resolve()
template = (root / 'testsuites/421_bulk_Si_hse_rt/inputfile').read_text()
template = template.replace("'Si_rps.dat'", "'deliberately_missing_pseudo.dat'")
template = template.replace('num_kgrid = 4, 4, 4', 'num_kgrid = 1, 1, 1')
template = template.replace('nproc_k=4', 'nproc_k=1')
cases = {
    'source': ("exx_ace_support='source'", 'source-support ACE is deferred'),
    'finite_range': ('hse_sr_tolerance=1d-3', 'Finite-range Gamma exchange is deferred'),
    'nan_range': ('hse_sr_tolerance=NaN', 'Finite-range Gamma exchange is deferred'),
    'infinite_range': ('hse_sr_tolerance=Infinity', 'Finite-range Gamma exchange is deferred'),
}
for name, (setting, message) in cases.items():
    with tempfile.TemporaryDirectory(prefix='hse-scope-') as directory:
        result = subprocess.run(
            [str(executable)], input=template.replace("xc = 'hse06'", "xc = 'hse06'\n  " + setting),
            cwd=directory, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=30, env={**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
        assert result.returncode != 0 and message in result.stdout, (name, result.stdout[-5000:])
        assert not list(Path(directory).rglob('*wfn*')), 'Input rejection wrote wavefunction data'
        print(name + ': rejected before loading GS data')

for name, setting in {
    'group': "hse_sr_group_fft='direct'",
    'wave': "yn_hse_sr_wave_first='y'",
    'return_sum': "yn_hse_sr_return_sum='y'",
    'sr_radius': 'hse_sr_radius=1d0',
    'ordinary_gamma': '',
}.items():
    with tempfile.TemporaryDirectory(prefix='hse-scope-') as directory:
        result = subprocess.run(
            [str(executable)], input=template.replace("xc = 'hse06'", "xc = 'hse06'\n  " + setting),
            cwd=directory, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=30, env={**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
        if name == 'ordinary_gamma':
            assert 'Unprepared ps_format' in result.stdout, result.stdout[-5000:]
            assert 'error(s) in input' not in result.stdout, result.stdout[-5000:]
            print(name + ': passed input validation, stopped at deliberately invalid pseudo fixture')
        else:
            # Upstream namelist failures use STOP without a nonzero exit code.
            assert 'error(s) in input' in result.stdout, result.stdout[-5000:]
            print(name + ': removed namelist input rejected')
        assert not list(Path(directory).rglob('*wfn*')), 'Input fixture wrote wavefunction data'
