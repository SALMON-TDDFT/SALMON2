"""Audit/link a local native pseudopotential input configuration probe. No network or CI.

BUILD_RECORD NEW_OUTPUT [--threads N]; all outputs must be external and new.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    record = json.loads(args.build.read_text())
    engine = Path(record['executable'])
    build = engine.parent
    if (record.get('schema') != 'salmonpy.hse-rt-state-build/1'
            or record.get('status') != 'built-and-audited'
            or hashlib.sha256(engine.read_bytes()).hexdigest() != record['sha256']):
        raise ValueError('unchanged audited build required')
    config = (build / 'config.h').read_text()
    if re.search(r'^#define (USE_MPI|USE_OPENACC)\b', config, re.M):
        raise ValueError('serial CPU build required')
    if args.threads < 1:
        raise ValueError('positive thread count required')
    output = args.output.resolve()
    if any(output.is_relative_to(path.resolve())
           for path in [root, build, Path(record['source_root'])]):
        raise ValueError('new external output required')
    output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(root / 'salmonpy/examples/si_intro_data/Si_rps.dat', output / 'Si_rps.dat')
    probe = output / 'probe.f90'
    shutil.copy2(root / 'tests/native/test_pp_grid_config.f90', probe)
    checker = root / 'tools/check_forbidden_eigensolvers.py'
    flags = (build / 'src/CMakeFiles/salmon.dir/flags.make').read_text()
    options = [v for v in shlex.split(re.search(r'^Fortran_FLAGS\s*=\s*(.*)$', flags, re.M)[1])
               if not v.startswith('-J')]
    includes = shlex.split(re.search(r'^Fortran_INCLUDES\s*=\s*(.*)$', flags, re.M)[1])
    fc = re.search(r'^CMAKE_Fortran_COMPILER:[^=]*=(.*)$',
                   (build / 'CMakeCache.txt').read_text(), re.M)[1].strip()
    commands = []
    with (output / 'build.log').open('w') as log:
        def run(command, cwd=output):
            commands.append([str(v) for v in command])
            subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=120,
                           env={**os.environ, 'OMP_NUM_THREADS': str(args.threads),
                                'OPENBLAS_NUM_THREADS': '1'})
        run([sys.executable, checker, '--source-root', record['source_root'] + '/src'])
        run([sys.executable, checker, '--source-root', output])
        run([fc, *options, '-g', '-fcheck=all', *includes, '-J' + str(output), '-c', probe, '-o', output / 'probe.o'])
        run([sys.executable, checker, '--object-tree', build])
        run([sys.executable, checker, '--object-tree', output])
        link = shlex.split((build / 'src/CMakeFiles/salmon.dir/link.txt').read_text())
        link[link.index('-o') + 1] = str(output / 'probe')
        link = [str(output / 'probe.o') if v == 'CMakeFiles/salmon.dir/main.f90.o' else v for v in link]
        run(link, build / 'src')
        run([output / 'probe'])
        for mode, expected in [('bad_grid', 'PP grid geometry is fixed'),
                               ('bad_layout', 'PP grid radial layout is fixed'),
                               ('bad_mapping', 'PP grid species mapping is fixed'),
                               ('bad_source', 'PP grid radial input is fixed')]:
            case = output / mode
            case.mkdir()
            shutil.copy2(output / 'Si_rps.dat', case / 'Si_rps.dat')
            command = [str(output / 'probe'), mode]
            commands.append(command)
            result = subprocess.run(command, cwd=case, capture_output=True, text=True,
                                    timeout=120, env={**os.environ, 'OMP_NUM_THREADS': str(args.threads),
                                                     'OPENBLAS_NUM_THREADS': '1'})
            (case / 'run.log').write_text(result.stdout + result.stderr)
            if result.returncode == 0 or expected not in result.stdout + result.stderr:
                raise AssertionError(f'{mode} did not reject the incompatible update')
    (output / 'result.json').write_text(json.dumps({
        'status': 'passed', 'scope': 'scalar pseudopotential grid ownership; not heterogeneous TDDFT',
        'engine_sha256': record['sha256'], 'threads': args.threads,
        'source_sha256': hashlib.sha256(probe.read_bytes()).hexdigest(), 'commands': commands,
    }, indent=2) + '\n')
    print('Pseudopotential grid configuration probe PASS:', output)


if __name__ == '__main__':
    main()
