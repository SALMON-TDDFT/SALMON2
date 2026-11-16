"""Explicitly link a serial fast-PR HSE command worker from audited build objects.

BUILD_JSON NEW_OUTPUT. No source preparation, downloads, or implicit compilation.
Only the independent fast-PR explicit-routing profile is accepted. The output
record identifies this worker, its adapter sources and the original audited build.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys


FAST_COMMIT = 'c9e8ce431f43eefd4a88998134adb4796f8ec837'


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            value.update(block)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    record = json.loads(args.build.read_text())
    executable = Path(record['executable']).resolve()
    build = executable.parent
    source = Path(record['source_root']).resolve()
    marker_path = source/'salmonpy-hse-rt-state.json'
    marker = json.loads(marker_path.read_text())
    if (record.get('schema') != 'salmonpy.hse-rt-state-build/1'
            or record.get('status') != 'built-and-audited'
            or digest(executable) != record.get('sha256')
            or record.get('base_commit') != FAST_COMMIT
            or marker.get('upstream_commit') != FAST_COMMIT
            or not marker.get('explicit_rt_routing')):
        raise ValueError('audited explicit-routing fast-PR HSE build required')
    config = (build/'config.h').read_text()
    if (not re.search(r'^#define USE_HSE\b', config, re.M)
            or re.search(r'^#define (USE_MPI|USE_OPENACC|USE_FFTW)\b', config, re.M)):
        raise ValueError('serial CPU fast-PR HSE/FFTE build required')
    out = args.output.resolve()
    if any(out.is_relative_to(parent) for parent in (root, build, source)):
        raise ValueError('new external output directory required')
    adapters = root/'salmonpy/native_adapters'
    if not adapters.is_dir():
        adapters = root.parent/'native_adapters'
    names = ('hse_checkpoint', 'tddft_context', 'fast_hse_worker')
    for name in names:
        if not (adapters/(name+'.f90')).is_file():
            raise ValueError('missing adapter source: '+name)
    flags = (build/'src/CMakeFiles/salmon.dir/flags.make').read_text()
    options = [v for v in shlex.split(re.search(r'^Fortran_FLAGS\s*=\s*(.*)$', flags, re.M)[1])
               if not v.startswith('-J')]
    includes = shlex.split(re.search(r'^Fortran_INCLUDES\s*=\s*(.*)$', flags, re.M)[1])
    compiler = re.search(r'^CMAKE_Fortran_COMPILER:[^=]*=(.*)$',
                        (build/'CMakeCache.txt').read_text(), re.M)[1].strip()
    link = shlex.split((build/'src/CMakeFiles/salmon.dir/link.txt').read_text())
    main_object = 'CMakeFiles/salmon.dir/main.f90.o'
    if link.count(main_object) != 1 or '-o' not in link:
        raise ValueError('unsupported original link command')
    original_objects = [(build/'src'/part).resolve() for part in link if part.endswith('.o') and part != main_object]
    if not original_objects or any(not path.is_file() for path in original_objects):
        raise ValueError('audited build object files are required')
    out.mkdir(parents=True, exist_ok=False)
    staged = out/'src'
    staged.mkdir()
    for name in names:
        shutil.copy2(adapters/(name+'.f90'), staged/(name+'.f90'))
    worker = out/'fast_hse_worker'
    commands = []
    checker = root/'tools/check_forbidden_eigensolvers.py'
    with (out/'build.log').open('w') as log:
        def run(command, **kwargs):
            command = [str(item) for item in command]
            commands.append(command)
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT,
                           timeout=300, **kwargs)
        run([sys.executable, checker, '--source-root', source/'src'])
        run([sys.executable, checker, '--source-root', staged])
        for name in names:
            run([compiler, *options, '-I'+str(out), *includes, '-J'+str(out),
                 '-c', staged/(name+'.f90'), '-o', out/(name+'.o')], cwd=out)
        run([sys.executable, checker, '--object-tree', build])
        run([sys.executable, checker, '--object-tree', out])
        link[link.index('-o')+1] = str(worker)
        link = [str(out/'fast_hse_worker.o') if item == main_object else item for item in link]
        link[1:1] = [str(out/'tddft_context.o'), str(out/'hse_checkpoint.o')]
        run(link, cwd=build/'src')
    receipt = {
        'schema': 'salmonpy.fast-hse-worker/1', 'status': 'built-and-audited',
        'executable': str(worker), 'sha256': digest(worker),
        'origin_sha256': record['sha256'],
        'build_origin': {'record_path': str(args.build.resolve()), 'record_sha256': digest(args.build),
                         'record': record, 'source_marker_sha256': digest(marker_path)},
        'source_hashes': {path.name: digest(path) for path in sorted(staged.iterdir())},
        'linked_object_hashes': {str(path): digest(path) for path in original_objects},
        'commands': commands,
        'protocol': 'salmonpy.fast-hse-worker/1',
        'scope': 'serial fast-PR HSE only; one process per state; Unix inputfile plus command stdin',
    }
    (out/'build.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(out/'build.json')


if __name__ == '__main__':
    main()
