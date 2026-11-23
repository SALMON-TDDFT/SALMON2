"""Run the short LD and synthetic SBE examples from an installed wheel.

No compiler, MPI, SALMON executable, network access or additional input is used.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from .errors import RuntimeFailure
from ._version import __version__


def run(directory=None):
    """Run both examples with one thread; preserve logs and a completion report."""
    if directory is None:
        output = Path(tempfile.mkdtemp(prefix='salmonpy-first-'))
    else:
        output = Path(directory).expanduser().resolve()
        output.mkdir(parents=True, exist_ok=False)
    report = {'version': __version__, 'status': 'running', 'output_directory': str(output),
              'purpose': 'installation check with synthetic models; not material-property validation',
              'threads': 1, 'calculations': []}
    environment = dict(os.environ)
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        environment[name] = '1'
    for example, extra in (('ld_response', []), ('sbe_response', [str(output/'sbe')])):
        started = time.perf_counter()
        log = output/(example + '.log')
        entry = {'example': example, 'log': str(log)}
        try:
            with log.open('w', encoding='utf-8') as stream:
                result = subprocess.run([sys.executable, '-I', '-X', 'utf8', '-m', 'salmonpy.examples.' + example, *extra],
                    stdout=stream, stderr=subprocess.STDOUT, env=environment, timeout=60)
            if result.returncode:
                raise RuntimeFailure(f'{example} failed; read {log}. Install a matching SALMONpy wheel; do not build automatically.')
        except (OSError, subprocess.TimeoutExpired, RuntimeFailure) as exc:
            entry.update(status='failed', error=str(exc), elapsed_seconds=time.perf_counter()-started)
            report['calculations'].append(entry)
            report['status'] = 'failed'
            (output/'summary.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
            raise RuntimeFailure(f'First calculation failed; see {output}/summary.json') from exc
        entry.update(status='passed', elapsed_seconds=time.perf_counter()-started)
        report['calculations'].append(entry)
    report['status'] = 'passed'
    (output/'summary.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    return report


def main():
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, 'reconfigure'):
            stream.reconfigure(encoding='utf-8')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='new output directory; otherwise create a temporary directory')
    args = parser.parse_args()
    try:
        report = run(args.output)
    except (OSError, RuntimeFailure) as exc:
        parser.exit(1, str(exc)+'\n')
    print('SALMONpy first calculations: PASS')
    print('LD: 1000 updates; synthetic SBE: 500 updates; one thread per calculation.')
    print('Output:', report['output_directory'])
    print('This confirms execution, not converged material properties.')


if __name__ == '__main__':
    main()
