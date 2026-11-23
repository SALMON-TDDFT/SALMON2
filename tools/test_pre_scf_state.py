"""Compile/run the owned pre-SCF probe against an existing audited serial HSE build.

BUILD_RECORD NEW_OUTPUT. Explicit local compilation only; no downloads.
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
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    record = json.loads(args.build.read_text())
    executable = Path(record['executable'])
    if (record.get('schema') != 'salmonpy.hse-rt-state-build/1'
            or record.get('status') != 'built-and-audited'
            or hashlib.sha256(executable.read_bytes()).hexdigest() != record['sha256']):
        raise ValueError('unchanged audited build required')
    build = executable.parent
    output = args.output.resolve()
    if any(output.is_relative_to(p.resolve()) for p in (root, build, Path(record['source_root']))):
        raise ValueError('new external output directory required')
    config = (build / 'config.h').read_text()
    if not re.search(r'^#define USE_HSE\b', config, re.M):
        raise ValueError('HSE build required')
    if re.search(r'^#define (USE_MPI|USE_OPENACC)\b', config, re.M):
        raise ValueError('serial CPU build required')
    output.mkdir(parents=True, exist_ok=False)
    source = output / 'src'
    source.mkdir()
    shutil.copy2(root / 'tests/native/test_pre_scf_state.f90', source / 'probe.f90')
    flags = (build / 'src/CMakeFiles/salmon.dir/flags.make').read_text()
    options = [v for v in shlex.split(re.search(r'^Fortran_FLAGS\s*=\s*(.*)$', flags, re.M)[1])
               if not v.startswith('-J')]
    includes = shlex.split(re.search(r'^Fortran_INCLUDES\s*=\s*(.*)$', flags, re.M)[1])
    compiler = re.search(r'^CMAKE_Fortran_COMPILER:[^=]*=(.*)$',
                         (build / 'CMakeCache.txt').read_text(), re.M)[1].strip()
    checker = root / 'tools/check_forbidden_eigensolvers.py'
    commands = []
    with (output / 'build.log').open('w') as log:
        def run(command, **kwargs):
            commands.append([str(v) for v in command])
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT, **kwargs)
        run([sys.executable, checker, '--source-root', Path(record['source_root']) / 'src'])
        run([sys.executable, checker, '--source-root', source])
        run([compiler, *options, *includes, '-J' + str(output), '-c', source / 'probe.f90',
             '-o', output / 'probe.o'], cwd=output)
        run([sys.executable, checker, '--object-tree', build])
        run([sys.executable, checker, '--object-tree', output])
        link = shlex.split((build / 'src/CMakeFiles/salmon.dir/link.txt').read_text())
        link[link.index('-o') + 1] = str(output / 'probe')
        link = [str(output / 'probe.o') if v == 'CMakeFiles/salmon.dir/main.f90.o' else v for v in link]
        run(link, cwd=build / 'src')
        run([output / 'probe'], cwd=output,
            env={**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
    (output / 'result.json').write_text(json.dumps({'status': 'passed', 'commands': commands,
        'probe_sha256': hashlib.sha256((source / 'probe.f90').read_bytes()).hexdigest(),
        'build_record': str(args.build.resolve())}, indent=2) + '\n')


if __name__ == '__main__':
    main()
