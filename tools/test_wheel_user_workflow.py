"""Install a supplied binary wheel and run five user examples in a fresh venv.

This source-checkout tool uses only the standard library. It never builds SALMONpy.
Run with the Python version matching the wheel; output must be outside the checkout.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import venv


def run_case(command, directory, label, env, *, clock=time.perf_counter):
    """Keep calculation outcome independent of optional timing instrumentation."""
    timing = {'status': 'unavailable', 'wall_seconds': None}
    try:
        start = clock()
    except Exception as exc:
        timing['error'] = str(exc)
        start = None
    record = {'label': label, 'command': list(map(str, command)),
              'status': 'not_started', 'exit_code': None,
              'timing': timing,
              'peak_memory': {'status': 'not_measured', 'bytes': None,
                              'scope': 'job including child processes'}}
    try:
        with (directory / (label + '.stdout.txt')).open('w') as stdout, \
             (directory / (label + '.stderr.txt')).open('w') as stderr:
            result = subprocess.run(command, cwd=directory, env=env,
                                    stdout=stdout, stderr=stderr, check=False)
        record.update(exit_code=result.returncode,
                      status='passed' if result.returncode == 0 else 'failed')
    except OSError as exc:
        record.update(status='launch_failed', error=str(exc))
    if start is not None:
        try:
            timing.update(status='measured', wall_seconds=clock() - start)
        except Exception as exc:
            timing['error'] = str(exc)
    (directory / (label + '.json')).write_text(json.dumps(record, indent=2) + '\n')
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('wheel', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--wheelhouse', type=Path,
                        help='offline directory containing compatible dependency wheels')
    parser.add_argument('--memory-budget', required=True, help='user allowance, or unknown')
    args = parser.parse_args()
    wheel = args.wheel.expanduser().resolve()
    output = args.output.expanduser().resolve()
    checkout = Path(__file__).resolve().parents[1]
    if not wheel.is_file() or wheel.suffix != '.whl':
        parser.error('provide an existing .whl file; source builds are not supported')
    if output.is_relative_to(checkout):
        parser.error('use a new output directory outside the checkout')
    output.mkdir(parents=True, exist_ok=False)
    report = {'schema': 'salmonpy.wheel-user-test/1', 'status': 'started',
              'wheel': str(wheel), 'wheel_sha256': hashlib.sha256(wheel.read_bytes()).hexdigest(),
              'python': sys.version, 'platform': platform.platform(),
              'memory_budget': args.memory_budget, 'memory_enforced': False,
              'cpu_policy': 'sequential, one numerical-library thread',
              'scope': 'wheel installation and five example smoke tests; not scientific convergence',
              'independent_numerical_validation': {'status': 'not_performed',
                  'note': 'Only checks present in the official examples run.'},
              'cases': []}
    def save():
        (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    env = dict(os.environ)
    for name in ('PYTHONPATH', 'PYTHONHOME'):
        env.pop(name, None)
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        env[name] = '1'
    save()
    try:
        venv.EnvBuilder(with_pip=True).create(output / 'venv')
        py = output / 'venv' / ('Scripts/python.exe' if os.name == 'nt' else 'bin/python')
        install = [str(py), '-I', '-m', 'pip', '--isolated', 'install',
                   '--only-binary=:all:', '--disable-pip-version-check']
        if args.wheelhouse:
            install += ['--no-index', '--find-links', str(args.wheelhouse.expanduser().resolve())]
        install.append(str(wheel))
        record = run_case(install, output, 'install', env)
        report['cases'].append(record); save()
        if record['status'] != 'passed':
            report['status'] = 'installation_failed'; save(); return 1
        probe = """import json,sys,platform
from pathlib import Path
from importlib.metadata import version
import salmonpy,numpy
p=Path(salmonpy.__file__).resolve()
assert p.is_relative_to(Path(sys.prefix).resolve()), str(p)
print(json.dumps(dict(salmonpy=version('salmonpy'),numpy=numpy.__version__,
    module_path=str(p),python=sys.version,platform=platform.platform()),indent=2))
"""
        record = run_case([str(py), '-I', '-c', probe], output, 'installed_identity', env)
        report['cases'].append(record); save()
        if record['status'] != 'passed':
            report['status'] = 'identity_failed'; save(); return 1
        for name, tail in [('ld_response', []), ('sbe_response', [str(output / 'sbe')]),
                           ('pulse_spectrum', []), ('geometry_preview', [str(output / 'grid.npz')]),
                           ('outputs_roundtrip', [str(output / 'archive')])]:
            record = run_case([str(py), '-I', '-m', 'salmonpy.examples.' + name, *tail],
                              output, name, env)
            report['cases'].append(record); save()
            print(name, record['status'], flush=True)
        report['status'] = 'passed' if all(r['status'] == 'passed' for r in report['cases']) else 'failed'
    except BaseException as exc:
        report.update(status='interrupted' if isinstance(exc, KeyboardInterrupt) else 'setup_failed',
                      error=str(exc))
        raise
    finally:
        save()
    print(output / 'report.json')
    return 0 if report['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
