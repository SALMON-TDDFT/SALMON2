"""Install a binary wheel in a fresh venv and run the complete Si starter.

No source builds, implicit runtime downloads, or checkout imports in the venv.
The runtime argument is the exact outer Actions artifact ZIP from the guide.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import venv
from test_wheel_user_workflow import run_case


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('wheel', type=Path)
    parser.add_argument('runtime_artifact', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--wheelhouse', type=Path)
    parser.add_argument('--memory-budget-gib', type=float, required=True)
    parser.add_argument('--threads', type=int, choices=range(1, 5), default=1)
    args = parser.parse_args()
    import math
    if not math.isfinite(args.memory_budget_gib) or args.memory_budget_gib <= 0:
        parser.error('finite positive memory budget required')
    wheel = args.wheel.expanduser().resolve()
    artifact = args.runtime_artifact.expanduser().resolve()
    output = args.output.expanduser().resolve()
    checkout = Path(__file__).resolve().parents[1]
    if not wheel.is_file() or wheel.suffix != '.whl' or not artifact.is_file():
        parser.error('existing binary wheel and runtime artifact required')
    if output.is_relative_to(checkout):
        parser.error('use a new output directory outside the checkout')
    output.mkdir(parents=True, exist_ok=False)
    report = {'schema': 'salmonpy.wheel-si-intro/1', 'status': 'started',
              'wheel_sha256': hashlib.sha256(wheel.read_bytes()).hexdigest(),
              'runtime_artifact_sha256': hashlib.sha256(artifact.read_bytes()).hexdigest(),
              'python': sys.version, 'platform': platform.platform(),
              'memory_budget_gib': args.memory_budget_gib, 'memory_enforced': False,
              'threads': args.threads, 'cases': [],
              'scope': 'fixed Si introductory regression, not optical convergence'}
    def save():
        (output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    env = dict(os.environ)
    for name in ('PYTHONPATH', 'PYTHONHOME', 'SALMONPY_TDDFT_WORKER'):
        env.pop(name, None)
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        env[name] = '1'
    def case(label, command):
        result = run_case(command, output, label, env)
        report['cases'].append(result)
        save()
        if result['status'] != 'passed':
            raise RuntimeError(f'{label} failed; inspect its logs and phase record')
        print(label, 'PASS', flush=True)
    save()
    try:
        venv.EnvBuilder(with_pip=True).create(output/'venv')
        python = output/'venv'/('Scripts/python.exe' if os.name == 'nt' else 'bin/python')
        install = [str(python), '-I', '-m', 'pip', '--isolated', 'install',
                   '--only-binary=:all:', '--disable-pip-version-check']
        if args.wheelhouse:
            install += ['--no-index', '--find-links', str(args.wheelhouse.resolve())]
        case('install', install+[str(wheel)])
        probe = """import json,sys
from pathlib import Path
import salmonpy,numpy
p=Path(salmonpy.__file__).resolve()
assert p.is_relative_to(Path(sys.prefix).resolve()), str(p)
print(json.dumps({'module':str(p),'numpy':numpy.__version__,'python':sys.version}))
"""
        case('identity', [str(python), '-I', '-c', probe])
        case('doctor', [str(python), '-I', '-m', 'salmonpy.doctor', '--json'])
        command = [str(python), '-I', '-m', 'salmonpy.examples.si_intro']
        case('prepare', command+['prepare', str(output/'si'), '--runtime-artifact', str(artifact)])
        case('calculate', command+['run', str(output/'si'), '--threads', str(args.threads),
                                  '--memory-budget-gib', str(args.memory_budget_gib)])
        report['starter_report'] = json.loads((output/'si/result/report.json').read_text())
        report['status'] = 'passed'
    except BaseException as exc:
        report.update(status='failed', error=str(exc))
        raise
    finally:
        save()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
